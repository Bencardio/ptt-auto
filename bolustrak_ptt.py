# -*- coding: utf-8 -*-
"""
bolustrak_ptt.py
================

Transit AP -> VG/aorte mesuré sur la série 2D_BOLUSTRAK (fluoroscopie de
déclenchement de l'ARM aortique, Philips), pour les patients sans perfusion.

Cette série s'arrête au déclenchement de l'ARM, souvent AVANT le pic dans le
VG/aorte. On mesure donc le transit sur des repères qui ne demandent pas le pic :
  - « pied » : intersection de la tangente au point de pente maximale avec la
    ligne de base (temps d'arrivée du bolus) ;
  - « pente max » : instant où la courbe monte le plus vite ;
  - « pic » : seulement si le pic est atteint dans les deux régions.
PTT = temps(VG/aorte) - temps(AP/cœur droit). Pour chaque patient, l'outil
dit si la pente maximale et le pic du VG ont bien été atteints avant la fin.

ATTENTION : mesure différente du PTT de perfusion (pic à pic VD -> VG) ; non
validée comme équivalente. Outil de recherche, non destiné au diagnostic.

Utilisation (dans le dossier qui contient aussi ptt_pipeline.py) :
    python bolustrak_ptt.py D:\\PTT_projet\\dicom
    python bolustrak_ptt.py D:\\PTT_projet\\dicom --patient P001
    python bolustrak_ptt.py D:\\PTT_projet\\dicom --patient P001 --serie 1201
Le dossier donné contient un sous-dossier par patient. Résultats :
    resultats_bolustrak/resultats_bolustrak.csv  (mis à jour après chaque patient)
    resultats_bolustrak/controles/<patient>.png  (image de contrôle)

Version 0.2.1 (série trop courte = non mesurable ; ligne de base après la mise à l'équilibre ; pic AP non atteint = non mesurable)
"""

import argparse
import sys
import traceback
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
try:
    import ptt_pipeline as pp          # lecture DICOM et temps : mêmes fonctions
except ImportError:
    print("\nERREUR : placez bolustrak_ptt.py dans le même dossier que ptt_pipeline.py\n")
    sys.exit(1)

import pandas as pd
from scipy import ndimage
from scipy.signal import savgol_filter
import matplotlib.pyplot as plt

# =============================================================================
# PARAMÈTRES MODIFIABLES
# =============================================================================

# Mots qui, dans le nom de la série, désignent la fluoroscopie de déclenchement.
MOTS_BOLUSTRAK = ["bolustrak", "bolus trak", "bolus", "care bolus", "smartprep", "fluoro", "timing"]

# Nombre minimal d'images dans le temps.
MIN_IMAGES = 20

# Les premières images sont plus claires (mise à l'équilibre du signal, ~2 s).
# Le début de la ligne de base est détecté automatiquement quand le signal
# moyen s'est stabilisé (à TOLERANCE_EQUILIBRE près), au plus tard à
# EQUILIBRE_MAX_S. La ligne de base dure ensuite DUREE_BASE_S.
TOLERANCE_EQUILIBRE = 0.03
EQUILIBRE_MAX_S = 6.0
DUREE_BASE_S = 3.0

# Nombre minimal d'images après la ligne de base pour chercher le passage du produit.
MIN_IMAGES_APRES_BASE = 10

# Écart minimal (s) entre les temps d'arrivée du groupe précoce (cœur droit/AP)
# et du groupe tardif (VG/aorte) : en dessous, il n'y a pas deux groupes
# distincts (par exemple une AP dont l'arrivée s'étale).
ECART_MIN_S = 3.0

# Un pixel n'est retenu que si son rehaussement maximal dépasse cette fraction
# du rehaussement des pixels les plus brillants (écarte le parenchyme pulmonaire).
SEUIL_REHAUSSEMENT = 0.25

# Bornes de plausibilité du transit (s).
PTT_MIN_S, PTT_MAX_S = 2.0, 20.0

# =============================================================================


def choisir_serie_bolustrak(series, numero=None):
    if numero is not None:
        for s in series.values():
            if s["serie_num"] == numero:
                return s, "série %d imposée" % numero
        return None, "série %d introuvable" % numero
    cand = [s for s in series.values()
            if pp.contient(s["description"], MOTS_BOLUSTRAK)
            and s["n_coupes"] == 1 and s["images_par_coupe"] >= MIN_IMAGES]
    if not cand:
        trop_courtes = [s for s in series.values() if pp.contient(s["description"], MOTS_BOLUSTRAK)]
        if trop_courtes:
            return None, "série BolusTrak trop courte (%d images)" % max(s["images_par_coupe"] for s in trop_courtes)
        return None, "aucune série BolusTrak"
    cand.sort(key=lambda s: s["images_par_coupe"], reverse=True)
    return cand[0], "série BolusTrak" + (" (la plus longue de %d)" % len(cand) if len(cand) > 1 else "")


def debut_equilibre(pile, temps):
    """Indice de la 1re image après la mise à l'équilibre du signal."""
    moy = pile.mean(axis=0)
    corps = moy > 0.2 * np.percentile(moy, 99)
    m = pile[:, corps].mean(axis=1)
    ref_m = (temps >= 3.0) & (temps <= EQUILIBRE_MAX_S)
    ref = float(np.median(m[ref_m])) if ref_m.sum() >= 2 else float(np.median(m[2:8]))
    stable = np.abs(m - ref) <= TOLERANCE_EQUILIBRE * abs(ref)
    i_max = int(np.searchsorted(temps, EQUILIBRE_MAX_S))
    for i in range(1, min(i_max, len(m) - 3) + 1):
        if stable[i:i + 3].all():
            return i
    return max(1, min(i_max, len(m) - 4))


def fenetre_base(temps, i_eq):
    i = np.arange(len(temps))
    m = (i >= i_eq) & (temps <= temps[i_eq] + DUREE_BASE_S)
    if m.sum() < 3:
        m = (i >= i_eq) & (i < i_eq + 3)
    return m


def lisser(y, dt):
    """Lissage léger (Savitzky-Golay sur ~2 s), sans décaler le temps."""
    w = max(5, int(round(2.0 / max(dt, 0.05))) | 1)
    if w > len(y) - 2:
        w = (len(y) - 2) - (1 - (len(y) - 2) % 2)
    return savgol_filter(y, w, 2) if w >= 5 else y.copy()


def kmeans_1d(v, n_iter=30):
    """Deux groupes de valeurs (k-moyennes en 1D) : renvoie les 2 centres triés."""
    c1, c2 = np.percentile(v, 15), np.percentile(v, 85)
    for _ in range(n_iter):
        g2 = np.abs(v - c2) < np.abs(v - c1)
        if g2.all() or (~g2).all():
            break
        n1, n2 = v[~g2].mean(), v[g2].mean()
        if n1 == c1 and n2 == c2:
            break
        c1, c2 = n1, n2
    return min(c1, c2), max(c1, c2)


def plus_grande_composante(masque):
    lab, n = ndimage.label(masque)
    if n == 0:
        return masque
    tailles = ndimage.sum(masque, lab, range(1, n + 1))
    return lab == (1 + int(np.argmax(tailles)))


def detecter_regions(pile, temps, i_eq):
    """AP/cœur droit (arrivée précoce) et VG/aorte (arrivée tardive).

    1. temps d'arrivée de chaque pixel rehaussé ;
    2. deux groupes de temps d'arrivée : précoce (cœur droit) et tardif (VG/aorte) ;
    3. dans le groupe précoce, s'il existe un sous-groupe plus précoce d'au moins
       1 s (veine cave supérieure, tronc brachio-céphalique), il est écarté."""
    dt = float(np.median(np.diff(temps)))
    lisse = ndimage.gaussian_filter(pile, sigma=(0, 1.2, 1.2))
    base_m = fenetre_base(temps, i_eq)
    i_fin_base = int(np.nonzero(base_m)[0].max()) + 1
    if len(temps) - i_fin_base < MIN_IMAGES_APRES_BASE:
        return None, "série trop courte : pas d'images après la ligne de base"
    base = lisse[base_m].mean(axis=0)
    bruit = lisse[base_m].std(axis=0)
    e = lisse - base[None]
    # rehaussement cherché APRÈS la ligne de base : les pixels qui ne sont
    # brillants qu'au début (mise à l'équilibre) ne sont pas retenus
    emax = e[i_fin_base:].max(axis=0)
    # seuil assez haut pour écarter le parenchyme pulmonaire (rehaussement faible)
    seuil = max(SEUIL_REHAUSSEMENT * np.percentile(emax, 99.5), 1e-9)
    cand = (emax > seuil) & (emax > 6 * np.median(bruit))
    cand = ndimage.binary_opening(cand, iterations=1)
    if cand.sum() < 30:
        return None, "trop peu de pixels rehaussés"

    # temps d'arrivée : 1re image >= 30 % du maximum du pixel
    rel = e / np.where(emax > 0, emax, 1)[None]
    rel[:i_fin_base] = 0
    depasse = rel >= 0.30
    ta = temps[np.argmax(depasse, axis=0)]
    c1, c2 = kmeans_1d(ta[cand])
    coupure = (c1 + c2) / 2
    ecart = float(c2 - c1)
    if ecart < ECART_MIN_S:
        return None, ("pas d'arrivée distincte dans le VG/aorte (écart %.1f s < %g s) : "
                      "série trop courte ?" % (ecart, ECART_MIN_S))

    precoce = cand & (ta < coupure)
    tardif = cand & (ta > coupure)
    if tardif.sum() < 10:
        return None, "pas d'arrivée visible dans le VG/aorte (série trop courte ?)"
    if precoce.sum() < 10:
        return None, "pas d'arrivée visible dans le cœur droit"
    # groupe précoce : veines (VCS, tronc brachio-céphalique) puis cœur droit/AP.
    # S'il y a deux sous-groupes séparés d'au moins 1 s, on garde le plus tardif.
    tp = ta[precoce]
    v1, v2 = kmeans_1d(tp)
    g2 = np.abs(tp - v2) < np.abs(tp - v1)
    if v2 - v1 >= 1.0 and g2.sum() >= 10:
        coeur = np.zeros_like(precoce)
        coeur[precoce] = g2
    else:
        coeur = precoce
    droit = plus_grande_composante(coeur)
    gauche = plus_grande_composante(tardif)
    veineux = precoce & ~droit
    # garder le cœur de chaque région (pixels les plus rehaussés)
    for m in (droit, gauche):
        if m.sum() > 20:
            m &= emax >= np.percentile(emax[m], 40)
    return {"droit": droit, "gauche": gauche, "veineux": veineux,
            "emax": emax, "coupure": coupure, "ecart_arrivee": ecart, "dt": dt}, ""


def reperes_courbe(t, y, base_m, i_eq):
    """Pied, pente maximale et pic d'une courbe (peut-être tronquée)."""
    dt = float(np.median(np.diff(t)))
    b0, sd = float(y[base_m].mean()), float(y[base_m].std() or 1e-9)
    e = y - b0
    e[:i_eq] = 0          # images avant l'équilibre : remplacées par la ligne de base
    ys = lisser(e, dt)
    n = len(ys)
    d = np.gradient(ys, t)
    d[:int(np.nonzero(base_m)[0].max()) + 1] = 0
    i_max = int(np.argmax(ys))
    # pente maximale cherchée avant le maximum observé (ignore la recirculation)
    i_s = int(np.argmax(d[:i_max + 1])) if i_max > 0 else int(np.argmax(d))
    pente = float(d[i_s])
    r = {"amplitude": float(ys[i_max]), "bruit": sd, "i_pente": i_s, "i_pic": i_max,
         "t_pente": float(t[i_s]), "pente": pente, "lisse": ys}
    if pente <= 0:
        r["t_pied"] = np.nan
    else:
        # tangente au point de pente maximale ∩ ligne de base (0)
        r["t_pied"] = float(t[i_s] - ys[i_s] / pente)
    # la pente maximale est fiable si la courbe ralentit ensuite nettement
    r["pente_atteinte"] = bool(i_s <= n - 3 and d[-1] < 0.6 * pente)
    # pic atteint : maximum avant les 2 dernières images et redescente nette,
    # ou plateau (saturation) d'au moins 3 images
    plateau = int(np.sum(ys[i_max:] >= 0.95 * ys[i_max]))
    r["pic_atteint"] = bool(i_max <= n - 3 and (ys[-1] < 0.85 * ys[i_max] or plateau >= 3))
    if r["pic_atteint"]:
        res = pp.analyser_courbe(t, ys)
        r["t_pic"] = float(res.get("t_pic_ajuste", res.get("t_pic_brut", t[i_max]))) if res.get("ok") else float(t[i_max])
    else:
        r["t_pic"] = np.nan
    return r


def analyser_patient(pid, dossier, numero_serie, dossier_png):
    ligne = {"ID": pid, "statut": "échec"}
    series = pp.resumer_series(pp.lire_entetes_patient(dossier))
    s, msg = choisir_serie_bolustrak(series, numero_serie)
    ligne["choix_serie"] = msg
    if s is None:
        ligne["statut"] = "sans BolusTrak"
        ligne["message"] = msg
        return ligne
    ligne.update({"serie_num": s["serie_num"], "serie_description": s["description"]})
    imgs = next(iter(s["coupes"].values()))
    pile, temps, imgs = pp.charger_pixels(imgs)
    ligne["n_images"] = len(imgs)
    ligne["source_temps"] = imgs[0].get("source_temps")
    if temps is None:
        ligne["statut"] = "non mesurable"
        ligne["message"] = "pas d'heure par image utilisable : temps inconnu"
        return ligne
    ligne["duree_serie_s"] = round(float(temps[-1]), 1)
    ligne["intervalle_images_s"] = round(float(np.median(np.diff(temps))), 2)

    i_eq = debut_equilibre(pile, temps)
    ligne["debut_ligne_base_s"] = round(float(temps[i_eq]), 1)
    apres = len(temps) - (int(np.nonzero(fenetre_base(temps, i_eq))[0].max()) + 1)
    if apres < MIN_IMAGES_APRES_BASE:
        ligne["statut"] = "non mesurable"
        ligne["message"] = ("série trop courte (%.1f s, %d images) : pas de passage du produit "
                            "après la ligne de base ; s'agit-il bien d'un BolusTrak ?"
                            % (float(temps[-1]), len(temps)))
        return ligne
    reg, err = detecter_regions(pile, temps, i_eq)
    if reg is None:
        ligne["message"] = err
        if "VG" in err:
            ligne["statut"] = "non mesurable"
        return ligne
    base_m = fenetre_base(temps, i_eq)
    courbes, rep = {}, {}
    for cle in ("droit", "gauche"):
        courbes[cle] = pile[:, reg[cle]].mean(axis=1)
        rep[cle] = reperes_courbe(temps, courbes[cle], base_m, i_eq)
    D, G = rep["droit"], rep["gauche"]

    ptt_pied = G["t_pied"] - D["t_pied"]
    ptt_pente = G["t_pente"] - D["t_pente"]
    ptt_pic = G["t_pic"] - D["t_pic"] if (D["pic_atteint"] and G["pic_atteint"]) else np.nan
    alertes = []
    if not G["pente_atteinte"]:
        alertes.append("pente max VG non atteinte (série trop courte) : transit sous-estimé")
    if not G["pic_atteint"]:
        alertes.append("pic VG non atteint")
    if not D["pic_atteint"]:
        alertes.append("pic AP/cœur droit non atteint (arrivée trop tardive pour la durée de la série)")
    elif np.isfinite(D["t_pic"]) and G["t_pied"] < D["t_pic"] - 1.0:
        alertes.append("arrivée VG avant le pic du cœur droit : régions à vérifier")
    if not (PTT_MIN_S <= ptt_pied <= PTT_MAX_S):
        alertes.append("transit (pied) hors bornes %g-%g s" % (PTT_MIN_S, PTT_MAX_S))
    if G["amplitude"] < 5 * G["bruit"] or D["amplitude"] < 5 * D["bruit"]:
        alertes.append("rehaussement faible")
    ligne.update({
        "PTT_pied_s": round(ptt_pied, 2), "PTT_pente_s": round(ptt_pente, 2),
        "PTT_pic_s": round(ptt_pic, 2) if np.isfinite(ptt_pic) else None,
        "t_pied_AP_s": round(D["t_pied"], 2), "t_pied_VG_s": round(G["t_pied"], 2),
        "pente_VG_atteinte": int(G["pente_atteinte"]), "pic_VG_atteint": int(G["pic_atteint"]),
        "pic_AP_atteint": int(D["pic_atteint"]),
        "pixels_AP": int(reg["droit"].sum()), "pixels_VG": int(reg["gauche"].sum()),
        "alertes": " ; ".join(alertes),
    })
    benignes = ("pic VG non atteint",)
    if not G["pente_atteinte"] or not D["pic_atteint"]:
        ligne["statut"] = "non mesurable"
    elif all(x in benignes for x in alertes):
        ligne["statut"] = "ok"
    else:
        ligne["statut"] = "à vérifier"
    dessiner(pid, pile, temps, reg, courbes, rep, ligne, dossier_png, base_m)
    return ligne


def dessiner(pid, pile, temps, reg, courbes, rep, ligne, dossier_png, base_m):
    dossier_png.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(1, 2, figsize=(13, 5.2))
    ax[0].imshow(pile[-1], cmap="gray")
    for cle, coul in (("veineux", "gray"), ("droit", "deepskyblue"), ("gauche", "red")):
        if reg[cle].any():
            ax[0].contour(reg[cle], levels=[0.5], colors=coul, linewidths=1.2)
    ax[0].set_title("%s - %s (dernière image)\nbleu = AP/cœur droit, rouge = VG/aorte, gris = veineux écarté"
                    % (pid, ligne.get("serie_description", "")), fontsize=9)
    ax[0].axis("off")
    for cle, coul, nom in (("droit", "deepskyblue", "AP/cœur droit"), ("gauche", "red", "VG/aorte")):
        r = rep[cle]
        y = courbes[cle] - courbes[cle][base_m].mean()
        ax[1].plot(temps, y, ".", color=coul, label=nom + " (mesuré)")
        ax[1].plot(temps, r["lisse"], "-", color=coul, lw=1)
        if np.isfinite(r["t_pied"]):
            ax[1].axvline(r["t_pied"], color=coul, ls="--", lw=1)
        ax[1].axvline(r["t_pente"], color=coul, ls=":", lw=1)
        if np.isfinite(r["t_pic"]):
            ax[1].axvline(r["t_pic"], color=coul, ls="-.", lw=1)
    ax[1].axhline(0, color="k", lw=0.5)
    ax[1].axvspan(temps[base_m][0], temps[base_m][-1], color="0.92", zorder=0)
    ax[1].set_xlabel("temps (s)")
    ax[1].set_ylabel("signal - ligne de base")
    ax[1].legend(fontsize=8, loc="upper left")
    pic = ligne.get("PTT_pic_s")
    ax[1].set_title("Transit pied = %.2f s | pente max = %.2f s | pic = %s\nstatut : %s%s"
                    % (ligne["PTT_pied_s"], ligne["PTT_pente_s"], "%.2f s" % pic if pic is not None else "non atteint",
                       ligne["statut"], ("\n" + ligne["alertes"]) if ligne["alertes"] else ""), fontsize=9)
    ax[1].text(0.01, 0.01, "tirets = pied ; pointillés = pente max ; tirets-points = pic",
               transform=ax[1].transAxes, fontsize=7, color="gray")
    fig.tight_layout()
    fig.savefig(dossier_png / (pid + ".png"), dpi=90)
    plt.close(fig)


COLONNES = ["ID", "statut", "PTT_pied_s", "PTT_pente_s", "PTT_pic_s", "pente_VG_atteinte",
            "pic_VG_atteint", "pic_AP_atteint", "t_pied_AP_s", "t_pied_VG_s", "debut_ligne_base_s", "duree_serie_s",
            "intervalle_images_s", "n_images", "serie_num", "serie_description", "source_temps",
            "pixels_AP", "pixels_VG", "alertes", "choix_serie", "message"]


def main():
    ap = argparse.ArgumentParser(description="Transit AP -> VG/aorte sur la série BolusTrak")
    ap.add_argument("dossier", help="dossier contenant un sous-dossier par patient")
    ap.add_argument("--patient", default=None, help="un seul patient (nom du sous-dossier)")
    ap.add_argument("--serie", type=int, default=None, help="numéro de série à imposer (avec --patient)")
    ap.add_argument("--sortie", default=None, help="dossier des résultats (défaut : resultats_bolustrak à côté du script)")
    a = ap.parse_args()
    base = Path(a.dossier)
    sortie = Path(a.sortie) if a.sortie else Path(__file__).resolve().parent / "resultats_bolustrak"
    sortie.mkdir(parents=True, exist_ok=True)
    patients = sorted([d for d in base.iterdir() if d.is_dir()])
    if a.patient:
        patients = [d for d in patients if d.name == a.patient]
    lignes = []
    for i, d in enumerate(patients, 1):
        print("[%d/%d] %s ..." % (i, len(patients), d.name), flush=True)
        try:
            ligne = analyser_patient(d.name, d, a.serie, sortie / "controles")
        except Exception as e:
            traceback.print_exc()
            ligne = {"ID": d.name, "statut": "échec", "message": "erreur : %s" % e}
        print("    %s  %s" % (ligne.get("statut"), ligne.get("PTT_pied_s", ligne.get("message", ""))))
        lignes.append(ligne)
        df = pd.DataFrame(lignes)
        df = df[[c for c in COLONNES if c in df.columns] + [c for c in df.columns if c not in COLONNES]]
        df.to_csv(sortie / "resultats_bolustrak.csv", index=False, sep=";", encoding="utf-8-sig")
    print("\nRésultats :", sortie / "resultats_bolustrak.csv")


if __name__ == "__main__":
    main()
