# -*- coding: utf-8 -*-
"""
site_ptt.py
===========

Petit « site » LOCAL pour analyser le PTT / PBVi par glisser-déposer.

    python site_ptt.py

Le navigateur s'ouvre sur http://127.0.0.1:8765 . Le site ne fonctionne que
sur cet ordinateur : il n'est pas accessible depuis le réseau et aucune image
ne quitte la machine. Tous les calculs sont ceux de ptt_pipeline.py (qui doit
se trouver dans le même dossier) ; ce fichier n'ajoute que l'interface.

Aucune bibliothèque supplémentaire n'est nécessaire : celles du pipeline
(pydicom, numpy, scipy, pandas, matplotlib, openpyxl) suffisent.

Dossiers créés à côté de ce fichier :
    site_travail/     copies des DICOM déposés (effaçables depuis le site)
    resultats_site/   un sous-dossier par analyse : resultats_PTT.csv (perfusion)
                      ou resultats_bolustrak.csv (sous-dossier ..._bolustrak),
                      .xlsx, controles/*.png, inventaire_series.csv

Outil de recherche, non destiné au diagnostic clinique.
"""

import base64
import hashlib
import pickle
import io
import json
import re
import shutil
import sys
import threading
import time
import traceback
import webbrowser
import zipfile
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, unquote, urlparse

BASE = Path(__file__).resolve().parent
sys.path.insert(0, str(BASE))
try:
    import ptt_pipeline as pp
except SystemExit:  # bibliothèque manquante : le pipeline a déjà affiché quoi faire
    input("Appuyez sur Entrée pour fermer.")
    raise
except ImportError:
    print("\nERREUR : ptt_pipeline.py doit être dans le même dossier que site_ptt.py :", BASE)
    input("Appuyez sur Entrée pour fermer.")
    sys.exit(1)

try:
    import bolustrak_ptt as bt     # analyse 2D BolusTrak (facultative)
except Exception as _e:            # fichier absent : le site marche sans ce mode
    bt = None
    print("Mode BolusTrak indisponible (%s)." % _e)
# Une relance manuelle tourne en parallèle de l'analyse de la cohorte :
# matplotlib (pyplot) n'aime pas tracer depuis deux fils à la fois.
_trace_verrou = threading.Lock()


def _un_trace_a_la_fois(f):
    def g(*a, **k):
        with _trace_verrou:
            return f(*a, **k)
    return g


for _m, _noms in ((pp, ("dessiner_controle", "dessiner_controle_segment")), (bt, ("dessiner",))):
    for _n in _noms:
        if _m is not None and hasattr(_m, _n):
            setattr(_m, _n, _un_trace_a_la_fois(getattr(_m, _n)))

if bt is not None and hasattr(pp, "MOTS_A_LIRE_EN_PLUS"):
    # séries BolusTrak toujours lues en entier (lecture rapide du pipeline)
    pp.MOTS_A_LIRE_EN_PLUS = sorted(set(pp.MOTS_A_LIRE_EN_PLUS) | set(bt.MOTS_BOLUSTRAK))

import numpy as np
import pandas as pd

PORT = 8765
TRAVAIL = BASE / "site_travail"
RESULTATS = BASE / "resultats_site"
PAGE = BASE / "page_site_ptt.html"

COLONNES_TABLEAU = ["ID", "presence_sequence_perfusion", "statut", "PTT_principal_s", "PTT_pic_s", "PTT_centroide_s",
                    "PTTc_pic", "PTTc_centroide", "PBVi_VD_pic_mL_m2", "PBVi_VG_pic_mL_m2",
                    "RR_ms", "serie_num", "coupe_retenue", "alertes", "message"]


# -----------------------------------------------------------------------------
# Lecture des en-têtes une seule fois par patient (inventaire + analyse)
# -----------------------------------------------------------------------------
_lire_entetes_orig = pp.lire_entetes_patient
_cache_entetes = {}
_cache_verrou = threading.Lock()


def lire_entetes_cache(dossier, complet=False):
    cle = str(Path(dossier).resolve())
    with _cache_verrou:
        if cle in _cache_entetes and (_cache_entetes[cle][1] or not complet):
            return _cache_entetes[cle][0]
    images = _lire_entetes_orig(dossier, complet=complet)
    with _cache_verrou:
        _cache_entetes[cle] = (images, complet)
    return images


pp.lire_entetes_patient = lire_entetes_cache
pp.GARDER_SERIE_RETENUE = True


def oublier_entetes(dossier):
    with _cache_verrou:
        _cache_entetes.pop(str(Path(dossier).resolve()), None)


# -----------------------------------------------------------------------------
# Éditeur de ROI : pixels d'une coupe gardés en mémoire (une seule à la fois)
# -----------------------------------------------------------------------------
_coupe_cache = {"cle": None, "val": None}
_coupe_verrou = threading.Lock()


def charger_coupe(lot, pid, serie, coupe):
    cle = (lot.id, pid, serie, coupe)
    with _coupe_verrou:
        if _coupe_cache["cle"] == cle:
            return _coupe_cache["val"]
    s = lot.series_retenues.get(pid)
    if s is None or s["serie_num"] != serie:
        # autre série que celle retenue : il faut relire les en-têtes du patient
        series = pp.resumer_series(pp.lire_entetes_patient(lot.patients[pid]["dossier"], complet=True))
        s = next((x for x in series.values() if x["serie_num"] == serie), None)
    if s is None:
        raise ValueError("série %s introuvable pour ce patient" % serie)
    positions = sorted(s["coupes"])
    if not 1 <= coupe <= len(positions):
        raise ValueError("coupe %s demandée, %d disponibles" % (coupe, len(positions)))
    pile, temps, _ = pp.charger_pixels(s["coupes"][positions[coupe - 1]])
    unite = "s"
    if temps is None:      # pas d'heure par image : axe en numéros d'image
        temps, unite = np.arange(pile.shape[0], dtype=float), "image"
    fin, _ = pp.fenetre_premier_passage(pile, temps)
    try:
        mm_px = float(s["coupes"][positions[coupe - 1]][0]["espacement"][0])
    except Exception:
        mm_px = None
    val = {"pile": pile, "temps": temps, "unite": unite, "fin": fin, "n_coupes": len(positions), "mm_px": mm_px}
    with _coupe_verrou:
        _coupe_cache.update(cle=cle, val=val)
    return val


def images_coupe(v):
    """Toutes les images de la coupe en niveaux de gris 8 bits (même fenêtrage),
    en un seul bloc base64 : la page les dessine elle-même, sans PNG."""
    pile = v["pile"]
    bas, haut = np.percentile(pile, [1, 99.7])
    gris = np.clip((pile - bas) * (255.0 / max(haut - bas, 1e-6)), 0, 255).astype(np.uint8)
    return base64.b64encode(gris.tobytes()).decode()


def courbes_roi(v, vd, vg, fenetre=None):
    pile, t, fin = v["pile"], v["temps"], v["fin"]
    deb = 0
    if fenetre:
        deb, fin = pp.indices_fenetre(t, fenetre, fin)
    res = {"t": [round(float(x), 2) for x in t], "fin": fin, "debut": deb, "fin_auto": v["fin"],
           "unite": v["unite"]}
    for nom, c in (("vd", vd), ("vg", vg)):
        if not c:
            continue
        m = pp.masque_disque(pile.shape[1:], *c)
        if m.sum() < 1:
            continue
        y = pile[:, m].mean(axis=1) - pile[:pp.N_IMAGES_BASE, m].mean()
        res[nom] = [round(float(x), 1) for x in y]
        a = pp.analyser_courbe(t[deb:fin], y[deb:fin])
        if a.get("ok"):
            res["pic_" + nom] = round(a["t_pic_brut"], 2)
    if "pic_vd" in res and "pic_vg" in res:
        res["ptt_pic"] = round(res["pic_vg"] - res["pic_vd"], 2)
    return res


def lire_fenetre(f):
    """[début, fin] de la fenêtre d'analyse en secondes (None = automatique)."""
    if not f:
        return None
    if not isinstance(f, (list, tuple)) or len(f) != 2:
        raise ValueError("fenêtre : 2 nombres attendus (début, fin en s)")
    f = [None if x in (None, "") else float(str(x).replace(",", ".")) for x in f]
    if f[0] is not None and f[1] is not None and f[1] <= f[0]:
        raise ValueError("fenêtre : la fin doit être après le début")
    return f if any(x is not None for x in f) else None


def lire_cercle(c):
    """[ligne, colonne, rayon] en pixels, ou None."""
    if c is None:
        return None
    if not isinstance(c, (list, tuple)) or len(c) != 3:
        raise ValueError("ROI : 3 nombres attendus (ligne, colonne, rayon)")
    c = [float(x) for x in c]
    if c[2] <= 0:
        raise ValueError("ROI : rayon nul")
    return [round(x, 1) for x in c]


# -----------------------------------------------------------------------------
# Lots d'analyse
# -----------------------------------------------------------------------------
def nom_sur(texte, defaut="patient"):
    """Nom de patient / de dossier sans caractère dangereux."""
    t = re.sub(r"[^A-Za-z0-9_.\- ]+", "_", str(texte)).strip(" .")
    return t[:80] or defaut


# Deux analyses distinctes, aux résultats séparés (dossiers et fichiers différents) :
#   perfusion : PTT VD -> VG (ptt_pipeline.py)
#   bolustrak : transit AP -> VG/aorte sur la série 2D_BOLUSTRAK (bolustrak_ptt.py)
MODES = {
    "perfusion": {"fichier": "resultats_PTT", "suffixe": "",
                  "statuts_faits": ("ok", "à vérifier", "sans perfusion")},
    "bolustrak": {"fichier": "resultats_bolustrak", "suffixe": "_bolustrak",
                  "statuts_faits": ("ok", "à vérifier", "sans BolusTrak", "non mesurable")},
    # courbes VD/VG exportées de Segment (« ROI Analysis », un .txt par patient)
    "segment": {"fichier": "resultats_PTT_segment", "suffixe": "_segment",
                "statuts_faits": ("ok", "à vérifier")},
}


def fichier_txt_segment(chemin):
    """Export Segment d'un patient : le fichier lui-même, ou le premier .txt du dossier."""
    p = Path(chemin)
    if p.is_file():
        return p
    txt = sorted(x for x in p.rglob("*") if x.is_file() and x.suffix.lower() == ".txt")
    if not txt:
        raise ValueError("aucun fichier .txt (export Segment) pour ce patient")
    return txt[0]


def ajouter_export_segment(lot, f, copie):
    """Un .txt = un patient (nom du fichier) ; un classeur Excel = un patient
    par onglet (nom de l'onglet). Autres fichiers ignorés."""
    f = Path(f)
    ext = f.suffix.lower()
    if ext == ".txt":
        lot.patients[nom_sur(f.stem)] = {"dossier": f, "fichiers": 1, "copie": copie}
    elif ext in pp.EXTENSIONS_CLASSEUR:
        feuilles = pp.feuilles_classeur(f)
        for feuille in feuilles:
            # onglet unique au nom par défaut (Sheet1, Feuil1...) : patient = nom du fichier
            defaut = len(feuilles) == 1 and str(feuille).strip().lower() in pp.NOMS_ONGLET_PAR_DEFAUT
            pid = nom_sur(f.stem if defaut else str(feuille))
            lot.patients[pid] = {"dossier": f, "feuille": feuille, "fichiers": 1, "copie": copie}
            lot.lignes.pop(pid, None)


def rr_des_analyses_perfusion():
    """RR (ms) lu dans les en-têtes DICOM et intervalle réel entre images (s),
    par ID patient, d'après les analyses de perfusion déjà faites sur ce PC :
    servent aux courbes Segment (temps exporté parfois inutilisable, et
    intervalle entre images différent du RR chez Philips)."""
    rr, pas = {}, {}
    for f in sorted(RESULTATS.glob("*/" + ETAT), key=lambda x: x.stat().st_mtime):
        try:
            e = json.loads(f.read_text(encoding="utf-8"))
        except Exception:
            continue
        if e.get("mode", "perfusion") != "perfusion":
            continue
        for pid, l in e.get("lignes", {}).items():
            if l.get("RR_ms") and str(l.get("source_RR", "")).startswith("en-tête"):
                rr[pid] = float(l["RR_ms"])
            if l.get("intervalle_images_s") and "reconstruit" not in str(l.get("source_temps", "")):
                pas[pid] = float(l["intervalle_images_s"])
    return rr, pas


class Lot:
    _compteur = 0

    def __init__(self, mode="perfusion", lid=None):
        Lot._compteur += 1
        self.mode = mode
        self.prefixe = prefixe_par_defaut()   # préfixe des codes pseudonymes de cette analyse
        self.id = lid or (time.strftime("%Y%m%d_%H%M%S") + "_%d" % Lot._compteur + MODES[mode]["suffixe"])
        self.dossier_travail = TRAVAIL / self.id
        self.dossier_resultats = RESULTATS / self.id
        self.patients = {}           # pid -> {"dossier": Path, "fichiers": int, "copie": bool}
        self.clinique = None         # Path du fichier clinique déposé
        self.lignes = {}             # pid -> ligne de résultats (dict)
        self.inventaire = {}         # pid -> liste des séries
        self.etat = "attente"        # attente | en_cours | termine
        self.en_cours = None
        self.fait = 0
        self.a_faire = 0
        self.debut = None
        self.journal = []
        self.etape = ""
        self.intervalle_segment = None   # s entre images, saisi sur le site (courbes Segment)
        self.series_retenues = {}    # pid -> série retenue (fichiers par coupe) : éditeur de ROI rapide
        self.file = []               # file d'attente : [(pid, choix manuel ou None)]
        self.file_relance = []       # relances manuelles pendant l'analyse (traitées en parallèle)
        self.relance = None          # patient en cours de relance manuelle
        self.principal_actif = False
        self.version = {}            # pid -> n° de relance : un ancien calcul ne l'écrase pas
        self.envoi = {"debut": None, "fin": None, "octets": 0}
        self.verrou = threading.Lock()

    def resume(self):
        lignes = [self.lignes[p] for p in sorted(self.patients) if p in self.lignes]
        return {
            "lot": self.id,
            "mode": self.mode,
            "etat": self.etat,
            "en_cours": self.en_cours,
            "relance": self.relance,
            "prefixe": self.prefixe,
            "relances_en_attente": [p for p, _ in self.file_relance],
            "fait": self.fait,
            "a_faire": self.a_faire,
            "duree_s": round(time.time() - self.debut) if self.debut else 0,
            "patients": [{"ID": p, "fichiers": d["fichiers"], "copie": d["copie"]}
                         for p, d in self.patients.items()],
            "clinique": self.clinique.name if self.clinique else (
                "cliniques.xlsx du dossier du site" if clinique_par_defaut() else None),
            "resultats": [propre({k: v for k, v in l.items() if k != "_sequences"}) for l in lignes],
            "inventaire": self.inventaire,
            "dossier_resultats": str(self.dossier_resultats),
            "journal": self.journal[-20:],
            "etape": self.etape,
            "envoi_s": round(self.envoi["fin"] - self.envoi["debut"], 1) if self.envoi["fin"] else None,
            "envoi_mo": round(self.envoi["octets"] / 1048576, 1),
        }


LOTS = {}
LOTS_VERROU = threading.Lock()
ANALYSE_VERROU = threading.Lock()


def lot_depuis(qs_ou_json):
    lid = qs_ou_json.get("lot")
    if isinstance(lid, list):
        lid = lid[0]
    with LOTS_VERROU:
        return LOTS.get(lid)


def propre(ligne):
    """Valeurs JSON-compatibles (NaN -> None, types numpy -> Python)."""
    out = {}
    for k, v in ligne.items():
        if hasattr(v, "item"):
            v = v.item()
        if isinstance(v, float) and v != v:
            v = None
        out[k] = v
    return out


def clinique_par_defaut():
    for nom in ("cliniques.xlsx", "cliniques.csv"):
        if (BASE / nom).exists():
            return BASE
    return None


def choisir_dossier_natif():
    """Fenêtre Windows/Mac « Choisir un dossier » ouverte par le serveur : le
    chemin est connu, les DICOM sont lus sur place, sans copie ni envoi.
    Lancée dans un processus séparé (tkinter n'aime pas les fils d'exécution)."""
    import subprocess
    code = ("import tkinter as tk; from tkinter import filedialog; r = tk.Tk(); r.withdraw(); "
            "r.attributes('-topmost', True); "
            "print(filedialog.askdirectory(title='Dossier DICOM à analyser') or '')")
    try:
        sortie = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True,
                                timeout=900, encoding="utf-8")
        if sortie.returncode != 0:
            return {"chemin": "", "erreur": "fenêtre de choix indisponible sur ce PC"}
        return {"chemin": sortie.stdout.strip()}
    except Exception as e:
        return {"chemin": "", "erreur": "fenêtre de choix indisponible (%s)" % e}


def inventaire_patient(dossier, mode="perfusion"):
    series = pp.resumer_series(pp.lire_entetes_patient(dossier))
    if mode == "bolustrak":
        choisi, _ = bt.choisir_serie_bolustrak(series)
    else:
        choisi, _ = pp.choisir_serie(series)
    sortie = []
    for s in sorted(series.values(), key=lambda s: s["serie_num"]):
        sortie.append({
            "SeriesNumber": s["serie_num"], "SeriesDescription": s["description"],
            "n_coupes": s["n_coupes"], "images_par_coupe": s["images_par_coupe"],
            "heure_presente": bool(s["temps_dispo"]),
            "auto": bool(choisi is not None and s is choisi),
        })
        s.pop("_candidats", None)
    return sortie


def contexte_analyse(lot):
    if lot.clinique:
        clin = pp.lire_clinique(lot.clinique.parent)
    elif clinique_par_defaut():
        clin = pp.lire_clinique(BASE)
    else:
        clin = {}
    rr_segment, pas_segment = rr_des_analyses_perfusion() if lot.mode == "segment" else ({}, {})
    return {"clin": clin, "rr_segment": rr_segment, "pas_segment": pas_segment,
            "png": lot.dossier_resultats / "controles"}


def lancer_analyse(lot):
    """Analyse séquentielle de la file d'attente (dans un fil d'exécution séparé).
    Les relances manuelles faites pendant ce temps (autre série, autre coupe,
    ROI redessinées) passent par lancer_relances, en parallèle."""
    with ANALYSE_VERROU:
        lot.debut = time.time()
        lot.fait = 0
        png = lot.dossier_resultats / "controles"
        png.mkdir(parents=True, exist_ok=True)
        ctx = contexte_analyse(lot)
        while True:
            with lot.verrou:
                if not lot.file:
                    lot.en_cours = None
                    lot.etape = ""
                    lot.principal_actif = False
                    if not lot.relance and not lot.file_relance:
                        lot.etat = "termine"
                    break
                pid, c = lot.file.pop(0)
                lot.en_cours = pid
                version = lot.version.get(pid, 0)
            traiter_patient(lot, pid, c, ctx, version, principal=True)
            lot.fait += 1


def lancer_relances(lot):
    """Relances manuelles demandées pendant l'analyse de la cohorte : faites
    tout de suite, sans attendre la fin du patient en cours."""
    ctx = contexte_analyse(lot)
    ctx["png"].mkdir(parents=True, exist_ok=True)
    while True:
        with lot.verrou:
            if not lot.file_relance:
                lot.relance = None
                if not lot.principal_actif and not lot.file:
                    lot.etat = "termine"
                break
            pid, c = lot.file_relance.pop(0)
            lot.relance = pid
            version = lot.version.get(pid, 0)
        traiter_patient(lot, pid, c, ctx, version, principal=False)


def traiter_patient(lot, pid, c, ctx, version, principal):
    """Analyse d'un patient et mise à jour des résultats. `version` : si une
    relance manuelle de ce patient a été demandée entre-temps, ce calcul est
    jeté (il ne doit pas écraser le choix de l'utilisateur)."""
    clin, png = ctx["clin"], ctx["png"]
    rr_segment, pas_segment = ctx["rr_segment"], ctx["pas_segment"]
    choix = {pid: c} if c else {}
    def etape(t):
        if principal:
            lot.etape = t
    dossier = lot.patients[pid]["dossier"]
    t0 = time.time()
    if lot.mode == "segment":
        try:
            etape("calcul sur les courbes Segment")
            feuille = lot.patients[pid].get("feuille")
            ligne = pp.analyser_export_segment(pid, dossier if feuille else fichier_txt_segment(dossier),
                                               clin, rr_segment, png, feuille,
                                               pas_segment.get(pid) or lot.intervalle_segment,
                                               (c or {}).get("fenetre"))
        except Exception as e:
            ligne = {"ID": pid, "statut": "échec", "message": "export Segment illisible : %s" % e,
                     "source_mesure": "Segment"}
        if clin.get(pid, {}).get("Groupe") is not None:
            ligne["Groupe"] = clin[pid]["Groupe"]
        if lot.lignes.get(pid, {}).get("commentaire"):
            ligne["commentaire"] = lot.lignes[pid]["commentaire"]
        if lot.version.get(pid, 0) != version:
            return
        lot.lignes[pid] = finaliser(ligne, lot)
        lot.journal.append("%s : %s %s" % (pid, ligne.get("statut"), ligne.get("PTT_principal_s", "") or ""))
        enregistrer(lot)
        return
    try:
        etape("lecture des en-têtes DICOM")
        deja = lot.series_retenues.get(pid)
        rapide = bool(c and c.get("serie") is not None and deja is not None
                      and deja.get("serie_num") == c["serie"])
        if rapide:
            # même série que l'analyse précédente (ROI redessinées, autre coupe) :
            # ses en-têtes sont déjà connus, inutile de relire tout l'examen
            images = [im for liste in deja["coupes"].values() for im in liste]
            with _cache_verrou:
                _cache_entetes[str(Path(dossier).resolve())] = (images, True)
        else:
            # série imposée : lecture complète (elle peut être dans un dossier sauté)
            images = pp.lire_entetes_patient(dossier, complet=bool(c and c.get("serie") is not None))
        t_entetes = time.time() - t0
        etape("calcul (%d fichiers DICOM lus en %.0f s)" % (len(images), t_entetes))
        if pid not in lot.inventaire:
            lot.inventaire[pid] = inventaire_patient(dossier, lot.mode)
        if lot.mode == "bolustrak":
            ligne = bt.analyser_patient(pid, dossier, (c or {}).get("serie"), png)
            if not rapide:
                ligne.update(pp.informations_examen(images))
                ligne["_sequences"] = pp.inventaire_sequences(pp.resumer_series(images))
        else:
            ligne = pp.analyser_patient(pid, dossier, clin, choix, png)
        if rapide:
            # seule la série retenue a été relue : informations d'examen et
            # séquences réalisées reprises de l'analyse précédente
            ancienne = lot.lignes.get(pid, {})
            ligne.pop("_sequences", None)
            for k in pp.COLONNES_EXAMEN + getattr(pp, "CLES_IDENTITE", []) + ["_sequences"]:
                if ancienne.get(k) is not None:
                    ligne[k] = ancienne[k]
        ligne["n_fichiers_lus"] = len(images)
        ligne["duree_entetes_s"] = round(t_entetes, 1)
        ligne["duree_calcul_s"] = round(time.time() - t0 - t_entetes, 1)
        oublier_entetes(dossier)   # mémoire : ne pas garder des centaines de patients
    except Exception as e:
        oublier_entetes(dossier)
        ligne = {"ID": pid, "statut": "échec", "message": "erreur : %s" % e}
        lot.dossier_resultats.mkdir(parents=True, exist_ok=True)
        with open(lot.dossier_resultats / "journal_erreurs.txt", "a", encoding="utf-8") as fj:
            fj.write("=== %s (%s) ===\n%s\n" % (pid, time.strftime("%Y-%m-%d %H:%M"),
                                              traceback.format_exc()))
    if clin.get(pid, {}).get("Groupe") is not None:
        ligne["Groupe"] = clin[pid]["Groupe"]
    if pid in choix:
        for k in ("choix_serie", "message"):
            if isinstance(ligne.get(k), str):
                ligne[k] = ligne[k].replace("(choix_series.csv)", "(choix manuel sur le site)")
        ligne["choix_manuel"] = "série %s" % choix[pid].get("serie") + (
            "" if lot.mode == "bolustrak" else ", coupe %s" % (choix[pid].get("coupe") or "auto")) + (
            ", ROI dessinées à la main" if choix[pid].get("roi") else "")
    if lot.version.get(pid, 0) != version:
        return     # relance manuelle demandée pendant ce calcul : on garde la sienne
    serie_ret = ligne.pop("_serie", None)
    if serie_ret:
        lot.series_retenues[pid] = serie_ret
        sauver_serie(lot, pid, serie_ret)
    if lot.lignes.get(pid, {}).get("commentaire"):
        ligne["commentaire"] = lot.lignes[pid]["commentaire"]   # le commentaire survit à une relance
    lot.lignes[pid] = finaliser(ligne, lot)   # nouveau calcul : validation manuelle remise à zéro
    lot.journal.append("%s : %s %s" % (pid, ligne.get("statut"),
                                       ligne.get("PTT_principal_s", ligne.get("PTT_pied_s", "")) or ""))
    enregistrer(lot)




# -----------------------------------------------------------------------------
# Codes pseudonymes : un code stable par patient (PTT-001, PTT-002...), gardé
# dans pseudonymes.csv à côté du site (même code d'une analyse à l'autre).
# Ce fichier contient les noms : il ne doit pas quitter le PC de l'hôpital.
# -----------------------------------------------------------------------------
PREFIXE_PSEUDO = "PTT-"                      # préfixe par défaut (modifiable sur le site)
FICHIER_PSEUDO = BASE / "pseudonymes.csv"
FICHIER_PREFIXE = BASE / "prefixe_pseudo.txt"   # dernier préfixe choisi sur le site
_pseudo_verrou = threading.Lock()
_pseudos = None     # {(préfixe, ID): code}


def prefixe_par_defaut():
    try:
        v = FICHIER_PREFIXE.read_text(encoding="utf-8").strip()
        return v if prefixe_valide(v) else PREFIXE_PSEUDO
    except OSError:
        return PREFIXE_PSEUDO


def prefixe_valide(v):
    return bool(re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_\-]{0,19}", str(v or "")))


def _prefixe_du_code(code):
    return re.sub(r"\d+$", "", str(code or ""))


def code_pseudo(pid, prefixe=None):
    """Code stable par patient et par préfixe (PTT-001, PTT-IMM-001...) :
    un même patient garde son code d'une analyse à l'autre ; chaque préfixe
    a sa propre numérotation."""
    global _pseudos
    prefixe = prefixe or PREFIXE_PSEUDO
    with _pseudo_verrou:
        if _pseudos is None:
            _pseudos = {}
            if FICHIER_PSEUDO.exists():
                try:
                    df = pd.read_csv(FICHIER_PSEUDO, sep=";", dtype=str, encoding="utf-8-sig")
                    df = df.astype(object).where(pd.notna(df), None)
                    for r in df.to_dict("records"):
                        if r.get("ID") and r.get("code_pseudo"):
                            _pseudos[(r.get("prefixe") or _prefixe_du_code(r["code_pseudo"]), r["ID"])] = r["code_pseudo"]
                except Exception:
                    _pseudos = {}
        if (prefixe, pid) not in _pseudos:
            nums = [int(c[len(prefixe):]) for (p, _), c in _pseudos.items()
                    if p == prefixe and str(c)[len(prefixe):].isdigit()]
            _pseudos[(prefixe, pid)] = "%s%03d" % (prefixe, max(nums, default=0) + 1)
            tmp = FICHIER_PSEUDO.with_suffix(".tmp")
            pd.DataFrame([{"ID": i, "code_pseudo": c, "prefixe": p} for (p, i), c in _pseudos.items()]).to_csv(
                tmp, sep=";", index=False, encoding="utf-8-sig")
            tmp.replace(FICHIER_PSEUDO)
        return _pseudos[(prefixe, pid)]


def finaliser(ligne, lot=None):
    """La validation manuelle (ok / échec) remplace le statut ; le statut donné
    par le programme reste dans statut_auto (décocher le rétablit)."""
    if ligne.get("ID"):
        ligne["code_pseudo"] = code_pseudo(str(ligne["ID"]), getattr(lot, "prefixe", None))
    ligne.setdefault("statut_auto", ligne.get("statut"))
    ligne["statut"] = ligne.get("validation_manuelle") or ligne["statut_auto"]
    return ligne


def fichier_csv(lot):
    return lot.dossier_resultats / (MODES[lot.mode]["fichier"] + ".csv")


_ecriture_verrou = threading.RLock()


def completer_identites(lot):
    """Patients analysés avec une version antérieure (sans date d'examen ni
    nom DICOM) : un seul en-tête DICOM est lu par patient pour compléter
    l'identité et la date d'examen. Rien n'est relu pour les autres."""
    if lot.mode == "segment":
        return 0
    import os
    n = 0
    for pid, l in list(lot.lignes.items()):
        if l.get("_nom_dicom") or l.get("_naissance") or l.get("date_examen"):
            continue
        d = lot.patients.get(pid)
        if not d or not Path(d["dossier"]).is_dir():
            continue
        essais = 0
        for racine, _, noms in os.walk(str(d["dossier"])):
            for nom in sorted(noms):
                if nom.upper() == "DICOMDIR" or nom.startswith("."):
                    continue
                essais += 1
                images = pp._entetes_fichier(Path(racine) / nom)
                if images:
                    info = pp.informations_examen(images[:1])
                    for k in ("date_examen", "heure_examen", "age_ans", "sexe", "identifiant_dicom",
                              "_naissance", "_nom_dicom") + (("date_naissance",) if pp.EXTRAIRE_DATE_NAISSANCE else ()):
                        if info.get(k) is not None and l.get(k) is None:
                            l[k] = info[k]
                    n += 1
                    break
                if essais >= 20:
                    break
            if essais >= 20 or l.get("date_examen") or l.get("_nom_dicom"):
                break
    if n:
        enregistrer(lot)
    return n


def enregistrer(lot):
    """Résultats (resultats_PTT.csv, même format que le pipeline, ou
    resultats_bolustrak.csv) + inventaire + identités, après chaque patient."""
    with _ecriture_verrou:
        _enregistrer(lot)


def _enregistrer(lot):
    lot.dossier_resultats.mkdir(parents=True, exist_ok=True)
    lignes = [lot.lignes[p] for p in sorted(lot.patients) if p in lot.lignes]
    if lignes and lot.mode == "bolustrak":
        tableau(lot).to_csv(fichier_csv(lot), index=False, sep=";", decimal=",", encoding="utf-8-sig")
    elif lignes:
        pp.ecrire_resultats(lignes, fichier_csv(lot))
    if lignes and lot.mode != "segment":
        try:
            pp.ecrire_sequences(lignes, lot.dossier_resultats)
        except Exception:
            pass
    if lignes:
        try:
            pp.ecrire_identites(lignes, lot.dossier_resultats)
        except Exception:
            traceback.print_exc()
    inv = [dict(ID=p, **s) for p, l in lot.inventaire.items() for s in l]
    if inv:
        pd.DataFrame(inv).to_csv(lot.dossier_resultats / "inventaire_series.csv",
                                 index=False, sep=";", encoding="utf-8-sig")
    sauver_etat(lot)


# -----------------------------------------------------------------------------
# Sauvegarde du lot (pour le rouvrir après fermeture du site)
# -----------------------------------------------------------------------------
ETAT = "etat_lot.json"


def fichier_serie(lot, pid):
    return lot.dossier_resultats / ".series" / (hashlib.md5(pid.encode("utf-8")).hexdigest() + ".pkl")


def sauver_etat(lot):
    """etat_lot.json : patients, dossiers DICOM, résultats (validations et ROI
    manuelles comprises). Écrit à côté du CSV après chaque changement."""
    etat = {"version": 1, "id": lot.id, "mode": lot.mode, "sauve_le": time.strftime("%Y-%m-%d %H:%M"),
            "clinique": str(lot.clinique) if lot.clinique else None,
            "patients": {p: {"dossier": str(d["dossier"]), "fichiers": d["fichiers"], "copie": d["copie"],
                             **({"feuille": d["feuille"]} if d.get("feuille") else {})}
                         for p, d in lot.patients.items()},
            "lignes": {p: propre(l) for p, l in lot.lignes.items()},
            "inventaire": lot.inventaire, "intervalle_segment": lot.intervalle_segment,
            "prefixe": lot.prefixe}
    f = lot.dossier_resultats / ETAT
    tmp = f.with_suffix(".tmp")
    tmp.write_text(json.dumps(etat, ensure_ascii=False, default=str), encoding="utf-8")
    tmp.replace(f)


def sauver_serie(lot, pid, serie):
    """Série retenue d'un patient (liste des fichiers par coupe) : éditeur de ROI
    rapide même après réouverture. Un fichier par patient."""
    try:
        f = fichier_serie(lot, pid)
        f.parent.mkdir(parents=True, exist_ok=True)
        with open(f, "wb") as fh:
            pickle.dump(serie, fh, protocol=pickle.HIGHEST_PROTOCOL)
    except Exception:
        pass


def lister_lots_sauves():
    sortie = []
    for f in sorted(RESULTATS.glob("*/" + ETAT), key=lambda x: x.stat().st_mtime, reverse=True)[:30]:
        try:
            e = json.loads(f.read_text(encoding="utf-8"))
        except Exception:
            continue
        lignes = e.get("lignes", {}).values()
        sortie.append({"id": e["id"], "mode": e.get("mode", "perfusion"), "sauve_le": e.get("sauve_le"),
                       "n_patients": len(e.get("patients", {})), "n_resultats": len(lignes),
                       "a_verifier": sum(1 for l in lignes if l.get("statut") == "à vérifier"),
                       "ouvert": e["id"] in LOTS})
    return sortie


def rouvrir_lot(lid):
    with LOTS_VERROU:
        if lid in LOTS:
            return LOTS[lid]
    f = RESULTATS / lid / ETAT
    if not f.exists() or f.parent.parent != RESULTATS:
        raise ValueError("analyse introuvable")
    e = json.loads(f.read_text(encoding="utf-8"))
    lot = Lot(e.get("mode", "perfusion"), lid=e["id"])
    lot.clinique = Path(e["clinique"]) if e.get("clinique") and Path(e["clinique"]).exists() else None
    for p, d in e.get("patients", {}).items():
        dossier = Path(d["dossier"])
        lot.patients[p] = {"dossier": dossier, "copie": d.get("copie", False),
                           "fichiers": d.get("fichiers", 0) if dossier.exists() else 0}
        if d.get("feuille"):
            lot.patients[p]["feuille"] = d["feuille"]
    lot.lignes = e.get("lignes", {})
    lot.inventaire = e.get("inventaire", {})
    lot.intervalle_segment = e.get("intervalle_segment")
    # analyses d'avant le choix du préfixe : celui de leurs codes
    lot.prefixe = e.get("prefixe") or next((_prefixe_du_code(l.get("code_pseudo")) for l in lot.lignes.values()
                                            if l.get("code_pseudo")), PREFIXE_PSEUDO)
    for p in lot.lignes:
        fs = fichier_serie(lot, p)
        if fs.exists():
            try:
                with open(fs, "rb") as fh:
                    lot.series_retenues[p] = pickle.load(fh)
            except Exception:
                pass
    lot.etat = "termine"
    lot.fait = lot.a_faire = len(lot.lignes)
    with LOTS_VERROU:
        LOTS[lot.id] = lot
    return lot


def tableau(lot):
    lignes = [lot.lignes[p] for p in sorted(lot.patients) if p in lot.lignes]
    df = pd.DataFrame(lignes)
    premieres = (["ID", "code_pseudo"] + [c for c in bt.COLONNES if c != "ID"] + ["validation_manuelle", "statut_auto", "commentaire"]
                 if lot.mode == "bolustrak" else pp.PREMIERES_COLONNES)
    cols = [c for c in premieres if c in df] + [c for c in df.columns
                                                if c not in premieres and not c.startswith("_")]
    df = pp.ajouter_nom_prenom(df[cols]) if len(df) else df
    for c in df.columns:   # 0/1, numéros de série... restent des entiers malgré les cases vides
        v = [x for x in df[c] if x is not None and x == x]
        if v and all(isinstance(x, int) and not isinstance(x, bool) for x in v):
            df[c] = df[c].astype("Int64")
    return df


# -----------------------------------------------------------------------------
# Serveur HTTP (127.0.0.1 uniquement)
# -----------------------------------------------------------------------------
class Gestionnaire(BaseHTTPRequestHandler):
    server_version = "SitePTT/1.0"

    def log_message(self, fmt, *args):  # console silencieuse
        pass

    # --- réponses ---------------------------------------------------------
    def envoyer(self, code, corps, type_="application/json; charset=utf-8", nom=None):
        if isinstance(corps, (dict, list)):
            corps = json.dumps(corps, ensure_ascii=False).encode("utf-8")
        elif isinstance(corps, str):
            corps = corps.encode("utf-8")
        self.send_response(code)
        self.send_header("Content-Type", type_)
        self.send_header("Content-Length", str(len(corps)))
        self.send_header("Cache-Control", "no-store")
        if nom:
            self.send_header("Content-Disposition", 'attachment; filename="%s"' % nom)
        self.end_headers()
        self.wfile.write(corps)

    def erreur(self, message, code=400):
        self.envoyer(code, {"erreur": message})

    def lire_json(self):
        n = int(self.headers.get("Content-Length") or 0)
        return json.loads(self.rfile.read(n).decode("utf-8") or "{}") if n else {}

    # --- GET -----------------------------------------------------------------
    def do_GET(self):
        u = urlparse(self.path)
        qs = parse_qs(u.query)
        if u.path in ("/", "/index.html"):
            return self.envoyer(200, PAGE.read_bytes(), "text/html; charset=utf-8")
        if u.path == "/web/pont.js":   # pont de la version en ligne : inutile ici
            return self.envoyer(200, b"", "text/javascript; charset=utf-8")
        lot = lot_depuis(qs)
        if u.path == "/api/etat":
            if not lot:
                return self.erreur("lot inconnu", 404)
            return self.envoyer(200, lot.resume())
        if u.path == "/api/png":
            if not lot:
                return self.erreur("lot inconnu", 404)
            pid = qs.get("patient", [""])[0]
            if pid not in lot.patients:
                return self.erreur("patient inconnu", 404)
            f = lot.dossier_resultats / "controles" / ("%s.png" % pid)
            if not f.exists():
                return self.erreur("pas d'image de contrôle", 404)
            return self.envoyer(200, f.read_bytes(), "image/png")
        if u.path == "/api/export":
            if not lot or not lot.lignes:
                return self.erreur("aucun résultat à exporter", 404)
            fmt = qs.get("format", ["csv"])[0]
            if fmt == "csv":
                enregistrer(lot)
                f = fichier_csv(lot)
                return self.envoyer(200, f.read_bytes(), "text/csv; charset=utf-8",
                                    "%s_%s.csv" % (MODES[lot.mode]["fichier"], lot.id))
            if fmt == "identites":
                completer_identites(lot)
                enregistrer(lot)
                f = lot.dossier_resultats / "identites_patients.csv"
                return self.envoyer(200, f.read_bytes(), "text/csv; charset=utf-8",
                                    "identites_patients_%s.csv" % lot.id)
            if fmt == "xlsx":
                f = fichier_csv(lot).with_suffix(".xlsx")
                completer_identites(lot)
                with pd.ExcelWriter(f) as w:
                    tableau(lot).to_excel(w, sheet_name="Resultats", index=False)
                    pp.identites([lot.lignes[p] for p in sorted(lot.patients) if p in lot.lignes]).to_excel(
                        w, sheet_name="Identites", index=False)
                    inv = [dict(ID=p, identifiant_dicom=lot.lignes.get(p, {}).get("identifiant_dicom"), **s)
                           for p, l in lot.inventaire.items() for s in l]
                    if inv:
                        pp.ajouter_nom_prenom(pd.DataFrame(inv)).to_excel(w, sheet_name="Inventaire_series", index=False)
                    for nom, feuille in (("sequences_realisees.csv", "Sequences_realisees"),
                                         ("sequences_detail.csv", "Sequences_detail")):
                        fs = lot.dossier_resultats / nom
                        if lot.mode != "segment" and fs.exists():
                            pd.read_csv(fs, sep=";", encoding="utf-8-sig").to_excel(
                                w, sheet_name=feuille, index=False)
                return self.envoyer(200, f.read_bytes(),
                                    "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                                    "%s_%s.xlsx" % (MODES[lot.mode]["fichier"], lot.id))
            if fmt == "zip":
                tampon = io.BytesIO()
                with zipfile.ZipFile(tampon, "w", zipfile.ZIP_DEFLATED) as z:
                    enregistrer(lot)
                    for f in sorted(lot.dossier_resultats.rglob("*")):
                        if f.is_file() and not f.name.endswith(".tmp"):
                            z.write(f, f.relative_to(lot.dossier_resultats))
                return self.envoyer(200, tampon.getvalue(), "application/zip",
                                    "%s_%s.zip" % (MODES[lot.mode]["fichier"], lot.id))
            return self.erreur("format inconnu")
        return self.erreur("page introuvable", 404)

    # --- POST ----------------------------------------------------------------
    def do_POST(self):
        # En-tête obligatoire : un autre site ouvert dans le navigateur ne peut
        # pas l'envoyer sans autorisation, donc ne peut pas piloter ce serveur.
        if self.headers.get("X-Site-PTT") != "1":
            return self.erreur("requête refusée", 403)
        u = urlparse(self.path)
        qs = parse_qs(u.query)
        try:
            if u.path == "/api/nouveau":
                mode = self.lire_json().get("mode") or "perfusion"
                if mode not in MODES:
                    return self.erreur("type d'analyse inconnu")
                if mode == "bolustrak" and bt is None:
                    return self.erreur("bolustrak_ptt.py est absent du dossier du site")
                lot = Lot(mode)
                with LOTS_VERROU:
                    LOTS[lot.id] = lot
                return self.envoyer(200, {"lot": lot.id, "prefixe": lot.prefixe})

            if u.path == "/api/lots":
                return self.envoyer(200, {"lots": lister_lots_sauves()})
            if u.path == "/api/rouvrir":
                try:
                    lot = rouvrir_lot(str(self.lire_json().get("id") or ""))
                except ValueError as e:
                    return self.erreur(str(e), 404)
                return self.envoyer(200, {"lot": lot.id, "mode": lot.mode, "prefixe": lot.prefixe})

            lot = lot_depuis(qs)
            if u.path == "/api/fichiers":
                return self.recevoir_paquet(lot)
            if u.path == "/api/choisir_dossier":
                return self.envoyer(200, choisir_dossier_natif())
            if u.path == "/api/clinique":
                return self.recevoir_clinique(lot, qs)

            data = self.lire_json()
            lot = lot or lot_depuis(data)
            if not lot:
                return self.erreur("lot inconnu (rechargez la page)", 404)
            if u.path == "/api/dossier_local":
                return self.dossier_local(lot, data)
            if u.path == "/api/analyser":
                return self.analyser(lot, data)
            if u.path == "/api/prefixe":
                v = str(data.get("prefixe") or "").strip()
                if not prefixe_valide(v):
                    return self.erreur("préfixe : lettres, chiffres, - ou _ (20 caractères au plus), ex. PTT-IMM-")
                lot.prefixe = v
                FICHIER_PREFIXE.write_text(v, encoding="utf-8")
                with lot.verrou:
                    for l in lot.lignes.values():
                        if l.get("ID"):
                            l["code_pseudo"] = code_pseudo(str(l["ID"]), v)
                if lot.lignes:
                    enregistrer(lot)
                return self.envoyer(200, {"ok": True, "prefixe": v})
            if u.path == "/api/commenter":
                pid, texte = data.get("patient"), str(data.get("texte") or "").strip()[:500]
                if pid not in lot.lignes:
                    return self.erreur("patient inconnu", 404)
                with lot.verrou:
                    if texte:
                        lot.lignes[pid]["commentaire"] = texte
                    else:
                        lot.lignes[pid].pop("commentaire", None)
                enregistrer(lot)
                return self.envoyer(200, {"ok": True})
            if u.path == "/api/valider":
                pid, v = data.get("patient"), data.get("valeur") or None
                if pid not in lot.lignes:
                    return self.erreur("patient inconnu", 404)
                if v not in (None, "ok", "échec"):
                    return self.erreur("valeur attendue : ok, échec ou vide")
                with lot.verrou:
                    if v:
                        lot.lignes[pid]["validation_manuelle"] = v
                    else:
                        lot.lignes[pid].pop("validation_manuelle", None)
                    finaliser(lot.lignes[pid], lot)
                enregistrer(lot)
                return self.envoyer(200, {"ok": True})
            if u.path in ("/api/roi_images", "/api/roi_courbes"):
                return self.roi(lot, data, u.path.endswith("images"))
            if u.path == "/api/effacer_copies":
                if lot.etat == "en_cours":
                    return self.erreur("analyse en cours")
                for pid, d in lot.patients.items():
                    if d["copie"]:
                        oublier_entetes(d["dossier"])
                shutil.rmtree(lot.dossier_travail, ignore_errors=True)
                for d in lot.patients.values():
                    if d["copie"]:
                        d["fichiers"] = 0
                return self.envoyer(200, {"ok": True})
            return self.erreur("action inconnue", 404)
        except Exception as e:
            traceback.print_exc()
            return self.erreur("erreur interne : %s" % e, 500)

    def recevoir_paquet(self, lot):
        """Plusieurs fichiers dans une seule requête (beaucoup plus rapide que
        une requête par fichier) : l'en-tête X-Fichiers donne, pour chaque
        fichier, [patient, chemin relatif, taille] ; le corps les met bout à bout."""
        if not lot:
            return self.erreur("lot inconnu (rechargez la page)", 404)
        if lot.etat == "en_cours":
            return self.erreur("analyse en cours")
        liste = json.loads(unquote(self.headers.get("X-Fichiers", "[]")))
        if lot.envoi["debut"] is None:
            lot.envoi["debut"] = time.time()
        for patient, chemin, taille in liste:
            pid = nom_sur(patient)
            morceaux = [nom_sur(m, "x") for m in str(chemin).replace("\\", "/").split("/")
                        if m not in ("", ".", "..")] or ["fichier"]
            dest = (lot.dossier_travail / pid).joinpath(*morceaux)
            dest.parent.mkdir(parents=True, exist_ok=True)
            reste = int(taille)
            with open(dest, "wb") as f:
                while reste > 0:
                    bloc = self.rfile.read(min(reste, 1 << 20))
                    if not bloc:
                        return self.erreur("envoi interrompu")
                    f.write(bloc)
                    reste -= len(bloc)
            with lot.verrou:
                d = lot.patients.setdefault(pid, {"dossier": lot.dossier_travail / pid,
                                                  "fichiers": 0, "copie": True})
                d["fichiers"] += 1
                lot.envoi["octets"] += int(taille)
            lot.lignes.pop(pid, None)
            lot.inventaire.pop(pid, None)
        for pid in {nom_sur(x[0]) for x in liste}:
            oublier_entetes(lot.dossier_travail / pid)
        if lot.mode == "segment":
            # classeur Excel déposé : un patient par onglet
            for pid in sorted({nom_sur(x[0]) for x in liste}):
                d = lot.patients.get(pid)
                if not d or d.get("feuille"):
                    continue
                classeurs = [x for x in Path(d["dossier"]).rglob("*")
                             if x.is_file() and x.suffix.lower() in pp.EXTENSIONS_CLASSEUR]
                if classeurs:
                    del lot.patients[pid]
                    try:
                        for c in classeurs:
                            ajouter_export_segment(lot, c, copie=True)
                    except Exception as e:
                        return self.erreur("classeur Excel illisible : %s" % e)
        lot.envoi["fin"] = time.time()
        return self.envoyer(200, {"ok": True, "n": len(liste)})

    def recevoir_clinique(self, lot, qs):
        if not lot:
            return self.erreur("lot inconnu (rechargez la page)", 404)
        nom = qs.get("nom", ["cliniques.xlsx"])[0].lower()
        ext = ".csv" if nom.endswith(".csv") else ".xlsx"
        d = lot.dossier_travail / "_clinique"
        shutil.rmtree(d, ignore_errors=True)
        d.mkdir(parents=True)
        f = d / ("cliniques" + ext)
        n = int(self.headers.get("Content-Length") or 0)
        f.write_bytes(self.rfile.read(n))
        try:
            clin = pp.lire_clinique(d)
        except Exception as e:
            shutil.rmtree(d, ignore_errors=True)
            return self.erreur("fichier clinique illisible : %s (colonne ID obligatoire)" % e)
        lot.clinique = f
        return self.envoyer(200, {"ok": True, "n_patients": len(clin)})

    def dossier_local(self, lot, data):
        """Dossier déjà présent sur ce PC (cohorte entière) : rien n'est copié."""
        if lot.etat == "en_cours":
            return self.erreur("analyse en cours")
        brut = str(data.get("chemin", "")).strip().strip('"')
        if not brut:
            return self.erreur("indiquez un dossier")
        p = Path(brut).expanduser()
        if not p.is_dir():
            return self.erreur("dossier introuvable : %s" % brut)
        if lot.mode == "segment":
            # un fichier .txt (export Segment) = un patient, nommé comme le fichier
            # ... et un classeur Excel = un patient par onglet
            avant = set(lot.patients)
            for x in sorted(p.rglob("*")):
                if not x.is_file() or x.name.startswith((".", "~$")):
                    continue
                try:
                    ajouter_export_segment(lot, x, copie=False)
                except Exception as e:
                    return self.erreur("classeur illisible (%s) : %s" % (x.name, e))
            nouveaux = sorted(set(lot.patients) - avant)
            if not nouveaux:
                return self.erreur("aucun export Segment (.txt) ni classeur Excel (.xlsx) dans %s" % p)
            return self.envoyer(200, {"ok": True, "n_patients": len(nouveaux), "patients": nouveaux})
        if data.get("un_seul_patient"):
            sous = [p]
        else:
            if (p / "dicom").is_dir():
                p = p / "dicom"
            sous = sorted(x for x in p.iterdir() if x.is_dir() and not x.name.startswith("."))
            if not sous:
                return self.erreur("aucun sous-dossier patient dans %s" % p)
        for x in sous:
            lot.patients[nom_sur(x.name)] = {"dossier": x, "fichiers": None, "copie": False}
        return self.envoyer(200, {"ok": True, "n_patients": len(sous),
                                  "patients": [nom_sur(x.name) for x in sous]})

    def roi(self, lot, data, avec_images):
        pid = data.get("patient")
        if lot.mode != "perfusion":
            return self.erreur("éditeur de ROI réservé aux perfusions")
        if pid not in lot.patients:
            return self.erreur("patient inconnu", 404)
        d = lot.patients[pid]
        if d["copie"] and not d["fichiers"]:
            return self.erreur("les copies DICOM de ce patient ont été effacées : redéposez-le")
        try:
            serie, coupe = int(data.get("serie")), int(data.get("coupe") or 1)
            v = charger_coupe(lot, pid, serie, coupe)
            res = courbes_roi(v, lire_cercle(data.get("vd")), lire_cercle(data.get("vg")),
                              lire_fenetre(data.get("fenetre")))
        except (TypeError, ValueError) as e:
            return self.erreur(str(e))
        if avec_images:
            res.update(gris=images_coupe(v), n_images=int(v["pile"].shape[0]), n_coupes=v["n_coupes"], coupe=coupe, serie=serie,
                       hauteur=int(v["pile"].shape[1]), largeur=int(v["pile"].shape[2]), mm_px=v.get("mm_px"))
        return self.envoyer(200, res)

    def analyser(self, lot, data):
        if not lot.patients:
            return self.erreur("aucun patient : déposez d'abord des fichiers DICOM"
                               + (" ou des exports Segment (.txt)" if lot.mode == "segment" else ""))
        lot.dossier_resultats.mkdir(parents=True, exist_ok=True)
        sauver_etat(lot)   # la liste des patients est retrouvée même si le site est fermé pendant l'analyse
        pids = [p for p in (data.get("patients") or list(lot.patients)) if p in lot.patients]
        if lot.mode == "segment" and "intervalle_s" in data:
            v = str(data.get("intervalle_s") or "").strip().replace(",", ".")
            try:
                lot.intervalle_segment = float(v) if v else None
            except ValueError:
                return self.erreur("intervalle entre images : nombre en secondes attendu (ex. 1,2)")
            if lot.intervalle_segment is not None and not 0.05 <= lot.intervalle_segment <= 10:
                return self.erreur("intervalle entre images : entre 0,05 et 10 s")
        choix = {}
        for pid, c in (data.get("choix") or {}).items():
            try:
                serie = int(c.get("serie")) if c.get("serie") not in (None, "") else None
                coupe = int(c.get("coupe")) if c.get("coupe") not in (None, "") else None
            except (TypeError, ValueError):
                return self.erreur("série et coupe doivent être des nombres entiers")
            if serie is not None or coupe is not None:
                choix[pid] = {"serie": serie, "coupe": coupe}
            try:
                fen = lire_fenetre(c.get("fenetre"))
            except (TypeError, ValueError) as e:
                return self.erreur(str(e))
            if fen:
                choix.setdefault(pid, {"serie": serie, "coupe": coupe})["fenetre"] = fen
            if c.get("roi") and lot.mode == "perfusion":
                try:
                    roi = {"vd": lire_cercle(c["roi"].get("vd")), "vg": lire_cercle(c["roi"].get("vg"))}
                except (TypeError, ValueError, AttributeError) as e:
                    return self.erreur(str(e))
                if not (roi["vd"] and roi["vg"] and serie is not None and coupe is not None):
                    return self.erreur("dessinez les deux ROI (VD et VG) avant de recalculer")
                choix[pid]["roi"] = roi
        if data.get("reprendre"):
            pids = [p for p in pids if lot.lignes.get(p, {}).get("statut") not in MODES[lot.mode]["statuts_faits"]
                    and not lot.lignes.get(p, {}).get("validation_manuelle")]   # échec validé à la main : on n'y touche plus
        relance = False
        with lot.verrou:
            en_cours = lot.principal_actif
            deja = {p for p, _ in lot.file}
            if en_cours and len(pids) <= 5:
                # relance d'un patient pendant l'analyse de la cohorte : faite tout
                # de suite, en parallèle (sans attendre la fin du patient en cours)
                for p in pids:
                    lot.version[p] = lot.version.get(p, 0) + 1   # tout calcul plus ancien sera ignoré
                lot.file_relance = [x for x in lot.file_relance if x[0] not in pids] + [(p, choix.get(p)) for p in pids]
                n_avant = len(lot.file)
                lot.file = [x for x in lot.file if x[0] not in pids]
                lot.a_faire -= n_avant - len(lot.file)
                relance = lot.relance is None
            else:
                lot.file += [(p, choix.get(p)) for p in pids if p not in deja]
                if not en_cours:
                    lot.fait = 0
                    lot.principal_actif = True
                lot.a_faire = lot.fait + len(lot.file) + (1 if en_cours else 0)
            lot.etat = "en_cours"
            if relance:
                lot.relance = pids[0] if pids else None   # réservé : un seul fil de relance à la fois
        if relance:
            threading.Thread(target=lancer_relances, args=(lot,), daemon=True).start()
        elif not en_cours:
            threading.Thread(target=lancer_analyse, args=(lot,), daemon=True).start()
        return self.envoyer(200, {"ok": True, "n": len(pids), "en_file": en_cours})


def main():
    port = PORT
    serveur = None
    for essai in range(PORT, PORT + 20):
        try:
            serveur = ThreadingHTTPServer(("127.0.0.1", essai), Gestionnaire)
            port = essai
            break
        except OSError:
            continue
    if serveur is None:
        print("ERREUR : aucun port libre entre %d et %d." % (PORT, PORT + 19))
        input("Appuyez sur Entrée pour fermer.")
        return
    url = "http://127.0.0.1:%d/" % port
    print("=" * 66)
    print(" Site PTT / PBVi (local) : %s" % url)
    print(" Les images restent sur cet ordinateur.")
    print(" Résultats enregistrés dans : %s" % RESULTATS)
    print(" Pour arrêter : fermez cette fenêtre (ou Ctrl+C).")
    print("=" * 66)
    if "--sans-navigateur" not in sys.argv:
        threading.Timer(1.0, lambda: webbrowser.open(url)).start()
    try:
        serveur.serve_forever()
    except KeyboardInterrupt:
        print("\nSite arrêté.")


if __name__ == "__main__":
    main()
