# -*- coding: utf-8 -*-
"""
ptt_pipeline.py
================

Calcul automatisé du temps de transit pulmonaire (PTT) et du volume sanguin
pulmonaire indexé (PBVi) à partir des séries DICOM de perfusion de premier
passage (IRM cardiaque), en rétrospectif.

Principe (indicateur-dilution, d'après Ricci 2019, Nelsson 2021, Seraphim 2021) :
  1. lecture des images DICOM de chaque patient ;
  2. choix de la série de perfusion de repos et de la coupe où le VD et le VG
     sont visibles ;
  3. détection automatique des cavités VD (rehaussement précoce) et VG
     (rehaussement tardif) ;
  4. courbes signal-temps, ajustement gamma-variate (élimine la recirculation) ;
  5. PTT = temps(VG) - temps(VD), par centre de gravité ET par pic à pic ;
  6. PTTc = PTT / RR ; nombre de battements = PTT / RR ;
     PBVi = battements x VES indexé du VD (Ricci) et PTT x débit cardiaque /
     surface corporelle (Seraphim) ;
  7. un tableau de résultats (CSV) et une image de contrôle (PNG) par patient.

Utilisation (dans l'invite de commandes, depuis le dossier du projet) :
    python ptt_pipeline.py inventaire          -> liste les séries trouvées
    python ptt_pipeline.py analyse             -> calcule PTT/PBVi pour tous
    python ptt_pipeline.py analyse --patient P001   -> un seul patient
    python ptt_pipeline.py analyse --reprendre      -> continue un calcul
                                                       interrompu
    python ptt_pipeline.py importer                 -> recopie les résultats
                                 dans une copie de la base Excel de recueil
    python ptt_pipeline.py segment                  -> lit les exports
                                 « ROI Analysis » de Segment (dossier
                                 exports_segment/) et calcule PTT/PBVi
    python ptt_pipeline.py diagnostic               -> heures DICOM de chaque
                                 série de perfusion (pour comprendre un échec)
    python ptt_pipeline.py importer --source fusion -> importe l'analyse
                                 automatique, remplacée par la mesure Segment
                                 pour les patients mesurés dans Segment
    Le fichier de résultats est enregistré après CHAQUE patient : une
    coupure (veille, redémarrage, Ctrl+C) ne fait perdre que le patient en
    cours.

Organisation attendue du dossier du projet :
    PTT_projet/
        ptt_pipeline.py          (ce fichier)
        cliniques.xlsx           (FACULTATIF : sans lui, le PTT et le
                                  PTTc sont calculés ; le PBVi est alors
                                  calculé dans la base Excel de recueil)
        Base_IRM_multiparametrique.xlsx (facultatif : pour le mode importer)
        choix_series.csv         (facultatif : choix manuel de la série)
        exports_segment/         (facultatif : exports Segment, un .txt par
                                  patient nommé avec son code, ex. CMH-012.txt)
        dicom/
            P001/   (tous les fichiers DICOM du patient, sous-dossiers permis)
            P002/
            ...
    Les résultats sont écrits dans PTT_projet/resultats/

Version 1.4.6 - outil de recherche, non destiné au diagnostic clinique.
"""

import argparse
import math
import sys
import time
import traceback
import warnings
from pathlib import Path

import numpy as np

warnings.filterwarnings("ignore")

try:
    import pydicom
    import pandas as pd
    from scipy import ndimage
    from scipy.optimize import curve_fit
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError as err:  # message clair pour un utilisateur novice
    print("\nERREUR : une bibliothèque Python manque :", err)
    print("Tapez dans l'invite de commandes :")
    print("  python -m pip install pydicom numpy scipy pandas matplotlib openpyxl\n")
    sys.exit(1)

VERSION = "1.4.8"

# =============================================================================
# 1. PARAMÈTRES MODIFIABLES
#    (ne modifiez que ce qui est entre les guillemets ou les nombres)
# =============================================================================

# Mots qui, dans le nom de la série (SeriesDescription), désignent une
# perfusion. La comparaison ignore majuscules/minuscules.
MOTS_PERFUSION = ["perf", "first pass", "firstpass", "fpp", "dyn", "srtfe",
                  "fgre tc", "time course"]

# Si AUCUNE série ne porte un nom de perfusion (MOTS_PERFUSION), faut-il
# chercher une série « dynamique » de remplacement ? False = non : le patient
# n'est pas analysé, statut « sans perfusion » et presence_sequence_perfusion = 0
# (évite d'analyser une mauvaise séquence et fait gagner du temps).
# True = ancien comportement (série dynamique de remplacement, à vérifier).
AUTORISER_SERIE_DE_REMPLACEMENT = False

# Utilisé par le site : garder la série retenue (liste des fichiers par coupe)
# dans le résultat, pour rouvrir les images sans relire tous les en-têtes.
GARDER_SERIE_RETENUE = False

# Lecture des en-têtes : nombre de fichiers lus en même temps (1 = un par un).
N_FILS_LECTURE = 8

# Lecture rapide : ne pas lire en entier les sous-dossiers qui ne contiennent
# qu'une série sans intérêt (ciné, LGE, T1...). False = tout lire (plus lent).
LECTURE_RAPIDE = True

# Séries à toujours lire en entier, en plus des perfusions (ex. BolusTrak).
MOTS_A_LIRE_EN_PLUS = ["bolustrak", "bolus trak", "bolus", "care bolus", "smartprep", "fluoro", "timing"]

# Séries dérivées à ne jamais utiliser (versions basse résolution, cartes...).
MOTS_DERIVES_EXCLUS = ["low res", "r_squared", "map"]

# Quand une perfusion existe en version recalée (MoCo) ET en version d'origine
# (préfixe « ORIG »), laquelle prendre ? True = la version recalée (les ROI
# restent dans les cavités malgré la respiration) ; False = l'origine.
PREFERER_MOCO = True

# Nombre de battements entre deux images de perfusion quand il n'est pas lisible
# dans les noms de séries (« Perf 2RR ») : sert seulement si le RR doit être
# déduit de l'intervalle entre images. None = détection par le nom, sinon 1.
RR_PAR_IMAGE = None

# Si plusieurs perfusions sont possibles, celles dont le rehaussement est
# inférieur à cette fraction du meilleur sont écartées (injection ratée).
SEUIL_REHAUSSEMENT_RELATIF = 0.3

# Mots qui désignent des séries à NE JAMAIS prendre pour une perfusion, même si
# elles comportent beaucoup d'images (TI scout / Look-Locker, cartographies,
# ciné, flux, rehaussement tardif...). Ignorés si le nom contient un mot de
# MOTS_PERFUSION.
MOTS_EXCLUS = ["scout", "look", "locker", "ti_", "ti ", "molli", "shmolli", "sasha", "map",
               "t1", "t2", "cine", "ciné", "trufi", "flow", "flux", "phase", "venc",
               "lge", "psir", "de ", "viab", "tagging", "localizer", "loc"]

# Nombre minimal d'images par coupe pour une série dont le NOM indique une
# perfusion (plus bas que MIN_IMAGES_DYNAMIQUES : perfusions acquises tous les
# 2 battements, ou courtes).
MIN_IMAGES_PERF_NOMMEE = 15

# Mots qui désignent une perfusion de STRESS (série à éviter).
MOTS_STRESS = ["stress", "adeno", "regadeno", "dobu", "effort"]

# Mots qui désignent explicitement une perfusion de REPOS (série à préférer).
MOTS_REPOS = ["rest", "repos", "rust"]

# Plusieurs perfusions dont l'une n'a que 3 ou 4 coupes (souvent « 3 sl » dans
# le nom) et les autres davantage (en général 8) : on analyse celle à peu de
# coupes, sur sa coupe du milieu (2e sur 3). False = ancien comportement.
PREFERER_PERFUSION_PEU_DE_COUPES = True
N_COUPES_MAX_PREFEREE = 4

# S'il reste plusieurs séries dynamiques sans indication repos/stress :
#   "derniere" = prendre la dernière acquise (le repos suit habituellement le
#                stress) ; "premiere" = prendre la première.
ORDRE_SI_AMBIGU = "derniere"

# Nombre minimal d'images dans le temps pour une coupe dynamique.
MIN_IMAGES_DYNAMIQUES = 30

# Nombre d'images de début utilisées comme ligne de base (avant le contraste).
N_IMAGES_BASE = 3

# Méthode principale pour le PTT : "pic" (pic à pic, recommandé pour les
# perfusions simple bolus où le sommet de la courbe VG est souvent écrêté) ou
# "centroide" (centre de gravité, utilisé par Seraphim 2021 avec une AIF dédiée).
METHODE_PRINCIPALE = "pic"

# Si l'heure de chaque image manque dans le DICOM : reconstruire le temps à
# partir du RR (1 image par battement supposée) ? True = oui, False = non.
AUTORISER_TEMPS_RECONSTRUIT = True

# Bornes de plausibilité du PTT (secondes) : en dehors, le résultat est signalé.
PTT_MIN_S, PTT_MAX_S = 2.0, 25.0

# Fenêtre du PREMIER PASSAGE : seules les images entre le début et l'arrivée du
# contraste + FENETRE_PREMIER_PASSAGE_S secondes servent à trouver VD/VG et à
# mesurer les courbes. Évite d'être trompé par une 2e injection pendant la série,
# par la recirculation ou par les organes abdominaux qui se rehaussent tard.
FENETRE_PREMIER_PASSAGE_S = 35.0

# Distance (mm) au-delà de laquelle un couple VD/VG devient peu plausible
# (les deux cavités sont voisines).
DISTANCE_CAVITES_OK_MM = 70.0

# Taille des ROI automatiques (des disques) : rayon = (1 - FRACTION_NOYAU_ROI) x
# rayon du plus grand cercle qui tient dans la cavité. Plus la valeur est
# grande, plus la ROI est petite. Les pixels
# moins rehaussés que FRACTION_REHAUSSEMENT_ROI x le rehaussement typique de la
# cavité (piliers, trabéculations, myocarde) sont aussi retirés.
FRACTION_NOYAU_ROI = 0.65
# Les ROI VD et VG automatiques ont la MÊME taille : le rayon commun est celui de
# la plus petite des deux cavités, multiplié par FACTEUR_TAILLE_ROI (1.0 = taille
# habituelle ; 1.3 = ROI 30 % plus larges, 0.7 = 30 % plus petites). Les pixels
# hors de la cavité ne sont jamais pris.
FACTEUR_TAILLE_ROI = 1.0
FRACTION_REHAUSSEMENT_ROI = 0.6

# Écart centroïde-pic (s) au-delà duquel une alerte est levée. Seuil empirique
# (contrôle qualité), aucun seuil n'est publié.
ECART_PIC_CENTROIDE_MAX_S = 1.5

# Dispersion : rapport des largeurs à mi-hauteur (LMH) VG / VD au-delà duquel la
# courbe VG est dite « étalée ». Seuil empirique, à ajuster sur vos témoins.
RAPPORT_LMH_ETALE = 1.8
# Courbes « propres » : R2 gamma >= R2_COURBE_PROPRE, pas de plateau (saturation).
R2_COURBE_PROPRE = 0.90

# Sauvegarde des courbes signal-temps VD/VG de chaque patient (dossier
# resultats/courbes/, un fichier .txt par patient). Ces fichiers se relisent
# sans les DICOM par le mode « segment » (site : « Courbes Segment ») : toutes
# les mesures (pics, centroïdes, LMH...) sont alors recalculées sur les mêmes
# courbes, avec la même fenêtre de premier passage et le même RR.
SAUVER_COURBES = True

# Informations d'examen lues dans les en-têtes DICOM (date d'examen, âge, taille,
# poids...). La date de naissance est une donnée identifiante : False = seul
# l'âge à l'examen est écrit (à choisir avant de partager les résultats).
EXTRAIRE_DATE_NAISSANCE = True

# Tableau des séquences réalisées (sequences_realisees.csv, 0/1 par type).
# Le type est reconnu d'après le nom de la série ; pour les cinés sans incidence
# dans le nom, d'après l'orientation des coupes (signalé « orientation » dans
# sequences_detail.csv, à vérifier). Mots cherchés (majuscules ignorées) :
MOTS_SEQUENCES = {
    "localisateur": ["survey", "scout", "localizer", "localiser", "loc", "3plane", "tri-pilot",
                     "senseref", "coilsurvey", "refscan", "calibration", "shim", "b0map"],
    "ti_scout": ["look", "locker", "looklocker", "ti scout", "tiscout", "ti_scout", "ti-scout", "ll"],
    "ecv": ["ecv"],
    "t2_etoile": ["t2*", "t2star", "t2 star", "t2_star", "t2etoile"],
    "t2_mapping": ["t2 map", "t2map", "t2_map", "t2-map", "t2 mapping", "t2prep map", "gradse", "t2 quant"],
    "t1_mapping": ["molli", "shmolli", "sasha", "savior", "t1 map", "t1map", "t1_map", "t1-map",
                   "t1 mapping", "mapit t1"],
    "t2_oedeme": ["stir", "t2w", "t2 w", "t2-w", "t2_w", "tse t2", "t2 tse", "t2_tse", "t2-tse",
                  "tirm", "black blood", "acut2e", "t2 bb", "t2_bb"],
    "angio": ["bolustrak", "bolus trak", "bolus", "care bolus", "smartprep", "fluoro", "timing",
              "ce-aorta", "ce aorta", "ce_aorta", "aorta", "cemra", "ce-mra", "ce_mra", "mra",
              "angio", "mip"],
    "flux_4d": ["4d flow", "4dflow", "4d_flow", "4d-flow", "flow 4d"],
    "flux": ["qflow", "q flow", "q-flow", "q_flow", "flow", "flux", "venc", "phase contrast", "pc"],
    "tagging": ["tag", "tagging", "3dtag", "tfepi", "fepi", "spamm", "grid", "dense"],
    "ege": ["ege", "early", "precoce", "précoce"],
    "lge": ["lge", "psir", "viab", "viability", "late", "delayed", "tardif", "scar",
            "3dscar", "ir_tfe", "ir-tfe", "irtfe", "ir_btfe", "ir-btfe", "irbtfe", "mag ir", "de_"],
    "cine": ["cine", "ciné", "sbtfe", "btfe", "b-tfe", "b tfe", "bffe", "ssfp", "trufi", "truefisp", "true fisp",
             "fiesta", "real time", "realtime"],
    "morphologie": ["haste", "ssfse", "sshtse", "ssh tse", "tse", "bb", "t1w", "t1 tse", "dixon",
                    "axial", "axiale"],
    # incidences des cinés (dans le nom)
    "2cav": ["2ch", "2 ch", "2-ch", "2_ch", "2cv", "2 cv", "2cav", "2 cav", "2c", "vla",
             "2chamber", "2 chamber", "deux cavites", "deux cavités"],
    "3cav": ["3ch", "3 ch", "3-ch", "3_ch", "3cv", "3 cv", "3cav", "3 cav", "3c", "lvot", "ccvg",
             "3chamber", "3 chamber", "trois cavites", "trois cavités"],
    "4cav": ["4ch", "4 ch", "4-ch", "4_ch", "4cv", "4 cv", "4cav", "4 cav", "4c", "hla",
             "4chamber", "4 chamber", "quatre cavites", "quatre cavités"],
    "petit_axe": ["sa", "sax", "short", "short axis", "petit axe", "pa", "sa stack", "stack"],
    "aortique": ["ao", "aorte", "aortic", "aortique", "valve ao", "ao valve", "av", "asc ao"],
    "vd": ["rv", "vd", "rvot", "ccvd", "tricuspid", "tricuspide", "rv inflow"],
    "autre_incidence": ["ax", "axial", "trans", "cor", "coronal", "sag", "sagittal", "mitral",
                        "mv", "pulmonary", "pulmonaire", "ap"],
    # injection de gadolinium (cartographies T1)
    "post_gado": ["post", "gd", "gado", "gad", "+c", "+gd", "apres", "après", "after", "postc",
                  "post-contrast", "postcontrast", "15min", "15 min", "10min", "10 min"],
    "natif": ["pre", "pré", "native", "natif", "nat", "precontrast", "pre-contrast"],
}

# Ordre et intitulé des colonnes de sequences_realisees.csv
COLONNES_SEQUENCES = [
    ("cine_2cav", "ciné 2 cavités"), ("cine_3cav", "ciné 3 cavités"),
    ("cine_4cav", "ciné 4 cavités"), ("cine_petit_axe", "ciné petit axe"),
    ("cine_aortique", "ciné valve aortique / aorte"), ("cine_vd", "ciné VD / CCVD"),
    ("cine_autre", "ciné autre incidence ou non déterminée"),
    ("tagging", "tagging"),
    ("t1_natif", "cartographie T1 native"), ("t1_post", "cartographie T1 post-gadolinium"),
    ("t1_indetermine", "cartographie T1 (avant/après gado non déterminé)"),
    ("ecv", "carte ECV"), ("t2_mapping", "cartographie T2"), ("t2_etoile", "T2*"),
    ("t2_oedeme", "T2 pondéré / STIR (œdème)"),
    ("perfusion", "perfusion (toute)"), ("perfusion_stress", "perfusion stress"),
    ("perfusion_repos", "perfusion repos"),
    ("ti_scout", "TI scout / Look-Locker"), ("ege", "rehaussement précoce (EGE)"),
    ("lge", "rehaussement tardif (LGE, tous types)"),
    ("flux", "flux (contraste de phase 2D)"), ("flux_4d", "flux 4D"),
    ("angio", "angiographie / bolus test"), ("morphologie", "morphologie (single-shot, sang noir...)"),
]

# =============================================================================
# 2. FONCTIONS UTILITAIRES
# =============================================================================


def temps_dicom_en_secondes(valeur):
    """Convertit une heure DICOM 'HHMMSS.ffffff' (ou date+heure) en secondes."""
    if valeur is None:
        return None
    s = str(valeur).strip()
    if not s:
        return None
    if len(s) >= 14 and s[:8].isdigit():  # format DT : AAAAMMJJHHMMSS.ffffff
        s = s[8:]
    try:
        hh = int(s[0:2])
        mm = int(s[2:4]) if len(s) >= 4 else 0
        ss = float(s[4:]) if len(s) > 4 else 0.0
        return hh * 3600 + mm * 60 + ss
    except ValueError:
        return None


def _get(ds, nom, defaut=None):
    try:
        v = getattr(ds, nom)
        return v if v != "" else defaut
    except Exception:
        return defaut


def position_coupe(ipp, iop, slice_location):
    """Position de la coupe le long de sa normale (mm), arrondie à 0,5 mm."""
    try:
        if ipp is not None and iop is not None:
            r = np.array(iop[:3], dtype=float)
            c = np.array(iop[3:], dtype=float)
            n = np.cross(r, c)
            return round(float(np.dot(np.array(ipp, dtype=float), n)) * 2) / 2
    except Exception:
        pass
    if slice_location is not None:
        try:
            return round(float(slice_location) * 2) / 2
        except Exception:
            pass
    return 0.0


def lister_fichiers(dossier):
    """Tous les fichiers d'un dossier et de ses sous-dossiers (os.walk : sous
    Windows, pas d'accès disque supplémentaire par fichier)."""
    import os
    sortie = []
    for racine, _, noms in os.walk(str(dossier)):
        r = Path(racine)
        sortie += [r / n for n in noms if n.upper() != "DICOMDIR" and not n.startswith(".")]
    return sortie


def _entetes_fichier(f):
    """En-têtes d'un fichier DICOM : liste d'images (1 pour un fichier classique,
    plusieurs pour un DICOM « Enhanced »), vide si le fichier n'est pas une image."""
    images = []
    try:
        ds = pydicom.dcmread(str(f), stop_before_pixels=True, force=True)
    except Exception:
        return []
    if _get(ds, "Rows") is None:
        return []
    serie = str(_get(ds, "SeriesInstanceUID", f.parent.name))
    info_commune = {
        "fichier": str(f),
        "serie_uid": serie,
        "serie_num": int(_get(ds, "SeriesNumber", 0) or 0),
        "description": str(_get(ds, "SeriesDescription", "") or ""),
        "rr_ms_tag": None,
        "fc_tag": None,
        "espacement": _get(ds, "PixelSpacing", None),
        "protocole": str(_get(ds, "ProtocolName", "") or ""),
        "iop": _orientation(_get(ds, "ImageOrientationPatient")),
        "examen": _infos_examen_ds(ds),
    }
    ni = _get(ds, "NominalInterval", None)
    if ni:
        info_commune["rr_ms_tag"] = float(ni)
    hr = _get(ds, "HeartRate", None)
    if hr:
        info_commune["fc_tag"] = float(hr)

    nframes = int(_get(ds, "NumberOfFrames", 1) or 1)
    if nframes <= 1:
        t_acq = temps_dicom_en_secondes(_get(ds, "AcquisitionTime", None))
        t_cont = temps_dicom_en_secondes(_get(ds, "ContentTime", None))
        t = t_acq if t_acq is not None else t_cont
        img = dict(info_commune)
        img.update({
            "trame": None,
            "position": position_coupe(_get(ds, "ImagePositionPatient"),
                                       _get(ds, "ImageOrientationPatient"),
                                       _get(ds, "SliceLocation")),
            "temps": t,
            "temps_acq": t_acq,
            "temps_contenu": t_cont,
            "trigger_ms": _get(ds, "TriggerTime", None),
            "instance": int(_get(ds, "InstanceNumber", 0) or 0),
        })
        images.append(img)
    else:  # DICOM « Enhanced » multi-images (fréquent chez Philips)
        iop_commun = None
        try:
            sh = ds.SharedFunctionalGroupsSequence[0]
            iop_commun = sh.PlaneOrientationSequence[0].ImageOrientationPatient
        except Exception:
            pass
        try:
            pf = ds.PerFrameFunctionalGroupsSequence
        except Exception:
            pf = [None] * nframes
        for k in range(nframes):
            g = pf[k] if k < len(pf) else None
            ipp = iop = t = rr = None
            try:
                ipp = g.PlanePositionSequence[0].ImagePositionPatient
            except Exception:
                pass
            try:
                iop = g.PlaneOrientationSequence[0].ImageOrientationPatient
            except Exception:
                iop = iop_commun
            if iop is None:
                iop = _get(ds, "ImageOrientationPatient")
            try:
                fc = g.FrameContentSequence[0]
                t = temps_dicom_en_secondes(_get(fc, "FrameAcquisitionDateTime")
                                            or _get(fc, "FrameReferenceDateTime"))
            except Exception:
                pass
            try:
                rr = float(g.CardiacSynchronizationSequence[0].RRIntervalTimeNominal)
            except Exception:
                pass
            img = dict(info_commune)
            img.update({
                "trame": k,
                "iop": _orientation(iop),
                "position": position_coupe(ipp, iop, None),
                "temps": t,
                "trigger_ms": None,
                "instance": k + 1,
            })
            if rr:
                img["rr_ms_tag"] = rr
            images.append(img)
    return images


def _orientation(iop):
    try:
        v = [round(float(x), 4) for x in iop]
        return tuple(v) if len(v) == 6 else None
    except Exception:
        return None


def _date_dicom(v):
    """'AAAAMMJJ' -> 'JJ/MM/AAAA' (None si absente ou illisible)."""
    s = str(v or "").strip()[:8]
    if len(s) == 8 and s.isdigit() and s != "00000000":
        return "%s/%s/%s" % (s[6:8], s[4:6], s[0:4])
    return None


def _infos_examen_ds(ds):
    """Informations patient / examen d'un en-tête (valeurs brutes, une fois par fichier)."""
    def num(nom):
        try:
            v = float(_get(ds, nom))
            return v if v > 0 else None
        except Exception:
            return None
    heure = str(_get(ds, "StudyTime", "") or "")
    return {
        "date": _date_dicom(_get(ds, "StudyDate") or _get(ds, "SeriesDate") or _get(ds, "AcquisitionDate")),
        "heure": "%s:%s" % (heure[:2], heure[2:4]) if len(heure) >= 4 and heure[:4].isdigit() else None,
        "naissance": _date_dicom(_get(ds, "PatientBirthDate")),
        "age": str(_get(ds, "PatientAge", "") or "").strip() or None,
        "sexe": str(_get(ds, "PatientSex", "") or "").strip().upper() or None,
        "id_patient": str(_get(ds, "PatientID", "") or "").strip() or None,
        "nom_patient": str(_get(ds, "PatientName", "") or "").strip() or None,
        "taille": num("PatientSize"),
        "poids": num("PatientWeight"),
        "fabricant": str(_get(ds, "Manufacturer", "") or "").strip() or None,
        "modele": str(_get(ds, "ManufacturerModelName", "") or "").strip() or None,
        "champ": num("MagneticFieldStrength"),
        "agent": str(_get(ds, "ContrastBolusAgent", "") or "").strip() or None,
        "volume": num("ContrastBolusVolume"),
        "debit": _premier_nombre(_get(ds, "ContrastFlowRate")),
        "duree_inj": _premier_nombre(_get(ds, "ContrastFlowDuration")),
        "dose": num("ContrastBolusTotalDose"),
        "concentration": num("ContrastBolusIngredientConcentration"),
        "voie": str(_get(ds, "ContrastBolusRoute", "") or "").strip() or None,
        "debut_inj": _heure_hhmmss(_get(ds, "ContrastBolusStartTime")),
    }


def _premier_nombre(v):
    """Débit / durée d'injection : un nombre ou une liste (un par phase)."""
    try:
        if isinstance(v, (list, tuple)) or hasattr(v, "__iter__") and not isinstance(v, str):
            v = list(v)[0]
        v = float(v)
        return v if v > 0 else None
    except Exception:
        return None


def _heure_hhmmss(v):
    t = temps_dicom_en_secondes(v)
    if t is None:
        return None
    return "%02d:%02d:%02d" % (t // 3600, (t % 3600) // 60, t % 60)


def _plus_frequent(valeurs):
    v = [x for x in valeurs if x not in (None, "")]
    if not v:
        return None
    from collections import Counter
    return Counter(v).most_common(1)[0][0]


# Colonnes d'informations d'examen (fichier de résultats)
COLONNES_EXAMEN = ["identifiant_dicom", "date_examen", "heure_examen", "date_naissance", "age_ans", "sexe",
                   "taille_cm_dicom", "poids_kg_dicom", "IMC_dicom", "SC_dicom_m2",
                   "fabricant_IRM", "modele_IRM", "champ_T", "produit_contraste",
                   "volume_contraste_mL", "dose_totale_contraste_mL", "concentration_contraste_mg_mL",
                   "debit_injection_mL_s", "duree_injection_s", "heure_debut_injection", "voie_injection"]


def informations_examen(images):
    """Date d'examen et de naissance (JJ/MM/AAAA), âge, sexe, taille, poids,
    IMC, surface corporelle (Mosteller), IRM et produit de contraste, quand ils
    sont présents dans les en-têtes DICOM."""
    ex = [im.get("examen") or {} for im in images]
    if not any(ex):
        return {}
    perf = [im.get("examen") or {} for im in images if contient(im.get("description", ""), MOTS_PERFUSION)]
    val = lambda k, l=ex: _plus_frequent(e.get(k) for e in l)
    out = {"date_examen": val("date"), "heure_examen": val("heure"), "sexe": val("sexe"),
           "fabricant_IRM": val("fabricant"), "modele_IRM": val("modele"), "champ_T": val("champ"),
           "produit_contraste": val("agent", perf) or val("agent"),
           "volume_contraste_mL": val("volume", perf) or val("volume"),
           "debit_injection_mL_s": val("debit", perf) or val("debit"),
           "dose_totale_contraste_mL": val("dose", perf) or val("dose"),
           "concentration_contraste_mg_mL": val("concentration", perf) or val("concentration"),
           "duree_injection_s": val("duree_inj", perf) or val("duree_inj"),
           "heure_debut_injection": val("debut_inj", perf) or val("debut_inj"),
           "voie_injection": val("voie", perf) or val("voie")}
    naissance = val("naissance")
    age = None
    if naissance and out["date_examen"]:
        jn, mn, an = (int(x) for x in naissance.split("/"))
        je, me, ae = (int(x) for x in out["date_examen"].split("/"))
        age = ae - an - ((me, je) < (mn, jn))
    elif val("age"):
        a = val("age")
        try:
            n_, u = int(a[:-1]), a[-1].upper()
            age = n_ if u == "Y" else (round(n_ / 12, 1) if u == "M" else (round(n_ / 52, 1) if u == "W" else 0))
        except Exception:
            age = None
    out["age_ans"] = age
    ident = val("id_patient")
    if ident and _est_identifiant(ident):
        out["identifiant_dicom"] = ident
    if EXTRAIRE_DATE_NAISSANCE:
        out["date_naissance"] = naissance
    # identité complète (fichier identites_patients.csv seulement, jamais dans les résultats)
    out["_naissance"] = naissance
    out["_nom_dicom"] = val("nom_patient")
    taille = val("taille")
    if taille:
        taille = taille * 100 if taille < 3 else taille       # DICOM : mètres
        out["taille_cm_dicom"] = round(taille, 1)
    poids = val("poids")
    if poids:
        out["poids_kg_dicom"] = round(poids, 1)
    if taille and poids and 100 <= taille <= 230 and 20 <= poids <= 300:
        out["IMC_dicom"] = round(poids / (taille / 100) ** 2, 1)
        out["SC_dicom_m2"] = round(math.sqrt(taille * poids / 3600), 3)   # Mosteller
    return {k: v for k, v in out.items() if v is not None}


def _lire_en_parallele(fichiers):
    """Lecture des en-têtes de plusieurs fichiers à la fois (beaucoup plus rapide
    sur un disque réseau, un disque USB ou avec un antivirus), ordre conservé."""
    if N_FILS_LECTURE <= 1 or len(fichiers) < 20:
        return [_entetes_fichier(f) for f in fichiers]
    from concurrent.futures import ThreadPoolExecutor
    with ThreadPoolExecutor(N_FILS_LECTURE) as ex:
        return list(ex.map(_entetes_fichier, fichiers))


def _series_utiles(description):
    return contient(description, MOTS_PERFUSION + MOTS_A_LIRE_EN_PLUS)


def lire_entetes_patient(dossier_patient, complet=False):
    """Lit les en-têtes DICOM (sans les pixels) et renvoie une liste d'images.

    Chaque image est un dict : série, description, position, temps, RR...
    Gère les fichiers classiques (1 image/fichier) et multi-images (Enhanced).

    Lecture rapide (LECTURE_RAPIDE = True, sauf complet=True), prévue pour les
    exports rangés une séquence par sous-dossier : dans chaque sous-dossier d'au
    moins 3 fichiers, seuls le premier et le dernier fichier sont lus (tous les
    dossiers en même temps). S'ils appartiennent à une même série qui n'est pas
    une perfusion ni un BolusTrak (ciné, LGE, T1, flux...), le reste du dossier
    n'est pas lu ; la série reste visible dans l'inventaire avec la mention
    « [lecture partielle] ». Parmi les perfusions, les séries dérivées (cartes,
    « low res ») et, s'il existe une perfusion de repos ou sans précision, les
    perfusions de stress ne sont pas lues en entier non plus : elles ne servent
    pas au calcul. Une série imposée à la main (choix manuel) relit tout.
    """
    fichiers = lister_fichiers(dossier_patient)
    if complet or not LECTURE_RAPIDE or AUTORISER_SERIE_DE_REMPLACEMENT:
        return [im for lot in _lire_en_parallele(fichiers) for im in lot]
    groupes = {}
    for f in fichiers:
        groupes.setdefault(f.parent, []).append(f)
    a_lire, images, echantillons = [], [], []
    for dossier, fs in groupes.items():
        fs = sorted(fs, key=lambda x: x.name)
        if len(fs) < 3:
            a_lire += fs
        else:
            echantillons.append((dossier, fs, [fs[0], fs[-1]]))
    # 1er et dernier fichier de chaque dossier, tous lus en même temps
    plats = [x for _, _, ech in echantillons for x in ech]
    lus_tous = _lire_en_parallele_force(plats)
    utiles, k = [], 0
    for dossier, fs, ech in echantillons:
        lus = lus_tous[k:k + len(ech)]
        k += len(ech)
        uids = {im["serie_uid"] for lot in lus for im in lot}
        descs = {im["description"] for lot in lus for im in lot}
        une_serie = all(lus) and len(uids) == 1 and all(len(lot) == 1 for lot in lus)
        if une_serie and not any(_series_utiles(d) for d in descs) and not _series_utiles(dossier.name):
            _marquer_partiel(lus, len(fs), images)      # dossier d'une seule série inutile : sauté
        else:
            utiles.append((dossier, fs, ech, lus, une_serie, " ".join(descs)))
    # perfusions qui ne servent pas au calcul : dérivées, et stress s'il existe un repos
    def perf_simple(d):
        return (contient(d, MOTS_PERFUSION) and not contient(d, MOTS_DERIVES_EXCLUS)
                and not contient(d, MOTS_A_LIRE_EN_PLUS))
    repos_dispo = any(u[4] and perf_simple(u[5]) and not contient(u[5], MOTS_STRESS) for u in utiles)
    for dossier, fs, ech, lus, une_serie, desc in utiles:
        inutile = une_serie and perf_simple(desc) and not contient(desc, MOTS_A_LIRE_EN_PLUS) and (
            False if not contient(desc, MOTS_STRESS) else repos_dispo)
        if une_serie and contient(desc, MOTS_PERFUSION) and contient(desc, MOTS_DERIVES_EXCLUS):
            inutile = True
        if inutile:
            _marquer_partiel(lus, len(fs), images)
        else:
            images += [im for lot in lus for im in lot]
            a_lire += [f for f in fs if f not in ech]
    images += [im for lot in _lire_en_parallele(a_lire) for im in lot]
    return images


def _marquer_partiel(lus, n_fichiers, images):
    for lot in lus:
        for im in lot:
            im["description"] += " [lecture partielle : %d fichiers]" % n_fichiers
            images.append(im)


def _lire_en_parallele_force(fichiers):
    """Comme _lire_en_parallele, mais en parallèle dès 2 fichiers."""
    if N_FILS_LECTURE <= 1 or len(fichiers) < 2:
        return [_entetes_fichier(f) for f in fichiers]
    from concurrent.futures import ThreadPoolExecutor
    with ThreadPoolExecutor(N_FILS_LECTURE) as ex:
        return list(ex.map(_entetes_fichier, fichiers))


def resumer_series(images):
    """Regroupe les images par série et par position de coupe."""
    series = {}
    for im in images:
        s = series.setdefault(im["serie_uid"], {
            "serie_uid": im["serie_uid"], "serie_num": im["serie_num"],
            "description": im["description"], "coupes": {}, "temps_min": None,
            "protocole": im.get("protocole", ""), "iop": im.get("iop")})
        s["coupes"].setdefault(im["position"], []).append(im)
        if im["temps"] is not None:
            s["temps_min"] = im["temps"] if s["temps_min"] is None else min(s["temps_min"], im["temps"])
    for s in series.values():
        n = [len(v) for v in s["coupes"].values()]
        s["n_images"] = sum(n)
        s["n_coupes"] = len(n)
        s["images_par_coupe"] = int(np.median(n)) if n else 0
        s["temps_dispo"] = all(im["temps"] is not None
                               for v in s["coupes"].values() for im in v)
    return fusionner_series_eclatees(series)


def fusionner_series_eclatees(series):
    """Certains exports Philips mettent chaque dynamique de la perfusion dans sa
    propre série (n° 120100001, 120100002, ... : 1 image par coupe chacune). Ces
    séries de même nom et de numéros qui se suivent sont regroupées en une seule."""
    groupes = {}
    for uid, s in series.items():
        if s["images_par_coupe"] <= 2 and contient(s["description"], MOTS_PERFUSION):
            groupes.setdefault(s["description"].strip().lower(), []).append(uid)
    for uids in groupes.values():
        if len(uids) < MIN_IMAGES_PERF_NOMMEE:
            continue
        uids.sort(key=lambda u: series[u]["serie_num"])
        nums = [series[u]["serie_num"] for u in uids]
        if nums[-1] - nums[0] > 2 * len(nums):        # pas une suite de numéros
            continue
        premiere = series[uids[0]]
        f = {"serie_uid": "fusion:" + premiere["serie_uid"], "serie_num": nums[0],
             "description": premiere["description"], "coupes": {}, "temps_min": None,
             "protocole": premiere.get("protocole", ""), "iop": premiere.get("iop"),
             "fusion": "%d séries %d à %d fusionnées (1 série par dynamique)" % (len(uids), nums[0], nums[-1])}
        for u in uids:
            s = series.pop(u)
            for pos, ims in s["coupes"].items():
                f["coupes"].setdefault(pos, []).extend(ims)
            if s["temps_min"] is not None:
                f["temps_min"] = s["temps_min"] if f["temps_min"] is None else min(f["temps_min"], s["temps_min"])
        n = [len(v) for v in f["coupes"].values()]
        f["n_images"], f["n_coupes"] = sum(n), len(n)
        f["images_par_coupe"] = int(np.median(n))
        f["temps_dispo"] = all(im["temps"] is not None for v in f["coupes"].values() for im in v)
        series[f["serie_uid"]] = f
    return series


def contient(texte, mots):
    t = texte.lower()
    return any(m.lower() in t for m in mots)


def presence_perfusion(series):
    """1 si au moins une série (non dérivée) porte un nom de perfusion, sinon 0."""
    return int(any(contient(s["description"], MOTS_PERFUSION)
                   and not contient(s["description"], MOTS_DERIVES_EXCLUS)
                   for s in series.values()))


def _mots_presents(texte, cle):
    """Un mot de MOTS_SEQUENCES[cle] est-il dans le texte ? Les mots courts
    (<= 3 caractères, ex. « sa », « de », « pc ») doivent être isolés."""
    import re
    for m in MOTS_SEQUENCES.get(cle, []):
        m = m.lower()
        if len(m) <= 3 or m in ("loc", "tag", "flow", "stack", "post", "late"):
            if re.search(r"(?<![a-z0-9])" + re.escape(m) + r"(?![a-z0-9])", texte):
                return True
        elif m in texte:
            return True
    return False


def _vecteur_normal(iop):
    if not iop:
        return None
    r, c = np.array(iop[:3]), np.array(iop[3:])
    n = np.cross(r, c)
    nn = np.linalg.norm(n)
    return n / nn if nn > 0.5 else None


def _incidence_grand_axe(normales, axe):
    """Incidences 2/3/4 cavités d'après l'orientation des coupes (repère patient
    DICOM : x vers la gauche, y vers l'arrière, z vers la tête).

    axe = grand axe du VG (normale du petit axe). Chaque coupe grand axe coupe le
    petit axe selon une droite ; son angle est mesuré depuis la paroi antérieure
    (projection de la tête) en tournant vers le septum (projection de la droite
    du patient) : 2 cavités ~ 0°, 3 cavités ~ 60°, 4 cavités ~ 120° (segments
    AHA espacés de 60°). Les incidences sont attribuées ensemble pour coller au
    mieux à ces angles. Renvoie une liste d'étiquettes (ou None)."""
    import itertools
    n = axe / np.linalg.norm(axe)
    if n @ np.array([1.0, -1.0, -1.0]) < 0:          # vers la pointe (gauche, avant, bas)
        n = -n
    a = np.array([0.0, 0.0, 1.0]) - n[2] * n
    if np.linalg.norm(a) < 0.2:
        return [None] * len(normales)
    a /= np.linalg.norm(a)
    s = np.cross(n, a)
    if s[0] > 0:                                       # vers la droite du patient
        s = -s
    angles = []
    for m in normales:
        d = np.cross(n, m)
        if np.linalg.norm(d) < 1e-6:
            angles.append(None)
            continue
        angles.append(math.degrees(math.atan2(d @ s, d @ a)) % 180)
    ref = {"cine_2cav": 0.0, "cine_3cav": 60.0, "cine_4cav": 120.0}
    dist = lambda x, y: min(abs(x - y) % 180, 180 - abs(x - y) % 180)
    # orientations regroupées (même incidence répétée : écart < 15°)
    groupes = []
    for i, ang in enumerate(angles):
        if ang is None:
            continue
        for g in groupes:
            if dist(g["angle"], ang) < 15:
                g["membres"].append(i)
                break
        else:
            groupes.append({"angle": ang, "membres": [i]})
    etiq = [None] * len(normales)
    if groupes and len(groupes) <= 3:
        meilleur = min(itertools.permutations(ref, len(groupes)),
                       key=lambda p: sum(dist(g["angle"], ref[c]) for g, c in zip(groupes, p)))
        for g, c in zip(groupes, meilleur):
            if dist(g["angle"], ref[c]) <= 45:
                for i in g["membres"]:
                    etiq[i] = c
    else:
        for g in groupes:
            c = min(ref, key=lambda k: dist(g["angle"], ref[k]))
            if dist(g["angle"], ref[c]) <= 30:
                for i in g["membres"]:
                    etiq[i] = c
    return etiq


def inventaire_sequences(series):
    """Types de séquences réalisées pour un patient.

    Renvoie {"colonnes": {clé: 0/1}, "non_classees": texte,
             "detail": [une ligne par série : numéro, nom, type, méthode]}.
    Plusieurs séries du même type (LGE répétés, plusieurs cinés dans la même
    incidence...) comptent une seule fois."""
    import re
    lignes = []
    for s in sorted(series.values(), key=lambda s: (s["serie_num"], s.get("temps_min") or 0)):
        nom = re.sub(r"\s*\[lecture partielle.*\]\s*$", "", s["description"])
        texte = (nom + " " + (s.get("protocole") or "")).lower().replace("_", " ")
        texte_brut = (nom + " " + (s.get("protocole") or "")).lower()
        mp = lambda c: _mots_presents(texte, c) or _mots_presents(texte_brut, c)
        cat, methode = None, "nom"
        ipc = s.get("images_par_coupe", 0)
        partielle = re.search(r"\[lecture partielle : (\d+) fichiers\]", s["description"])
        if partielle:          # seuls quelques fichiers lus : nombre d'images estimé
            ipc = max(ipc, int(partielle.group(1)) // max(s.get("n_coupes", 1), 1))
        if mp("ti_scout") and not contient(texte, MOTS_PERFUSION):
            cat = "ti_scout"
        elif mp("localisateur"):
            cat = "localisateur"
        elif mp("ecv"):
            cat = "ecv"
        elif mp("t2_etoile"):
            cat = "t2_etoile"
        elif mp("t2_mapping"):
            cat = "t2_mapping"
        elif mp("t1_mapping"):
            cat = "t1"
        elif contient(texte, MOTS_PERFUSION) and not mp("flux") and not mp("cine"):
            cat = ("perfusion_stress" if contient(texte, MOTS_STRESS) else
                   "perfusion_repos" if contient(texte, MOTS_REPOS) else "perfusion")
        elif mp("t2_oedeme"):
            cat = "t2_oedeme"
        elif mp("flux_4d"):
            cat = "flux_4d"
        elif mp("flux"):
            cat = "flux"
        elif mp("angio"):
            cat = "angio"
        elif mp("tagging"):
            cat = "tagging"
        elif mp("ege"):
            cat = "ege"
        elif mp("lge"):
            cat = "lge"
        elif mp("cine") and ipc >= 5:
            cat = "cine"
        elif mp("morphologie") or mp("cine"):
            cat = "morphologie"
        elif ipc >= 10 and s.get("n_coupes", 1) <= 20:
            cat, methode = "cine", "nombre d'images (nom non reconnu)"
        if cat == "cine":
            for c in ("2cav", "3cav", "4cav", "petit_axe", "aortique", "vd", "autre_incidence"):
                if mp(c):
                    cat = {"autre_incidence": "cine_autre"}.get(c, "cine_" + c)
                    break
        # pile de coupes : >= 3 positions, ou >= 2 quand seuls quelques fichiers ont été lus
        multi = s.get("n_coupes", 1) >= (2 if partielle else 3)
        lignes.append({"s": s, "nom": nom, "texte": texte, "categorie": cat, "methode": methode,
                       "ipc": ipc, "multi": multi})

    # --- cinés sans incidence dans le nom : orientation des coupes -----------------
    normale = lambda l: _vecteur_normal(l["s"].get("iop"))
    axe, source_axe = None, None
    for choix_axe in (lambda l: l["categorie"] == "cine_petit_axe",
                      lambda l: l["categorie"] == "cine" and l["multi"],
                      lambda l: l["categorie"].startswith("perfusion") if l["categorie"] else False):
        cands = [normale(l) for l in lignes if choix_axe(l) and normale(l) is not None]
        if cands:
            axe = cands[0]
            source_axe = "petit axe"
            break
    a_classer = [l for l in lignes if l["categorie"] == "cine" and normale(l) is not None]
    if axe is None:
        uniques = [normale(l) for l in a_classer if not l["multi"]]
        meilleur = None
        for i in range(len(uniques)):
            for j in range(i + 1, len(uniques)):
                c = np.cross(uniques[i], uniques[j])
                if meilleur is None or np.linalg.norm(c) > np.linalg.norm(meilleur):
                    meilleur = c
        if meilleur is not None and np.linalg.norm(meilleur) > 0.5:
            axe, source_axe = meilleur / np.linalg.norm(meilleur), "intersection des grands axes"
    grands_axes = []
    for l in lignes:
        if l["categorie"] != "cine":
            continue
        m = normale(l)
        multi = l["multi"]
        if m is None or axe is None:
            l["categorie"] = "cine_petit_axe" if (multi and m is None or multi and axe is None) else "cine_autre"
            l["methode"] = "pile de coupes" if l["categorie"] == "cine_petit_axe" else "incidence non déterminée"
            continue
        cos = abs(float(m @ axe))
        if cos >= 0.8:
            l["categorie"], l["methode"] = "cine_petit_axe", "orientation (à vérifier)"
        elif cos <= 0.5 and not multi:
            grands_axes.append(l)
        else:
            l["categorie"], l["methode"] = "cine_autre", "orientation (à vérifier)"
    if grands_axes:
        etiq = _incidence_grand_axe([normale(l) for l in grands_axes], axe)
        for l, e in zip(grands_axes, etiq):
            l["categorie"] = e or "cine_autre"
            l["methode"] = "orientation (%s, à vérifier)" % source_axe

    # --- cartographies T1 : avant ou après gadolinium -------------------------------
    avec_contraste = ("perfusion", "perfusion_stress", "perfusion_repos", "angio", "lge", "ege",
                      "ti_scout", "ecv")
    t_contraste = [l["s"]["temps_min"] for l in lignes
                   if l["categorie"] in avec_contraste and l["s"].get("temps_min") is not None]
    t0 = min(t_contraste) if t_contraste else None
    gado = any(l["categorie"] in avec_contraste for l in lignes)
    for l in lignes:
        if l["categorie"] != "t1":
            continue
        if mp_post(l["texte"]):
            l["categorie"] = "t1_post"
        elif mp_natif(l["texte"]):
            l["categorie"] = "t1_natif"
        elif not gado:
            l["categorie"], l["methode"] = "t1_natif", "aucune série avec gadolinium"
        elif t0 is not None and l["s"].get("temps_min") is not None:
            l["categorie"] = "t1_post" if l["s"]["temps_min"] > t0 else "t1_natif"
            l["methode"] = "heure par rapport à la 1re série avec gadolinium"
        else:
            l["categorie"], l["methode"] = "t1_indetermine", "heure d'acquisition absente"
    # --- perfusions sans repos/stress dans le nom : 2 séries -> stress puis repos --
    perfs = [l for l in lignes if l["categorie"] == "perfusion"
             and not contient(l["texte"], MOTS_DERIVES_EXCLUS)]
    nommees = [l for l in lignes if l["categorie"] in ("perfusion_stress", "perfusion_repos")]
    if len(perfs) == 2 and not nommees and all(l["s"].get("temps_min") is not None for l in perfs):
        perfs.sort(key=lambda l: l["s"]["temps_min"])
        for l, c in zip(perfs, ("perfusion_stress", "perfusion_repos")):
            l["categorie"], l["methode"] = c, "2 perfusions : la 1re = stress (ordre habituel, à vérifier)"

    presentes = {l["categorie"] for l in lignes}
    colonnes = {c: int(c in presentes) for c, _ in COLONNES_SEQUENCES}
    colonnes["perfusion"] = int(any(c and c.startswith("perfusion") for c in presentes))
    non_classees = sorted({l["nom"] for l in lignes if l["categorie"] is None})
    detail = []
    for l in lignes:
        s = l["s"]
        t = s.get("temps_min")
        detail.append({
            "SeriesNumber": s["serie_num"], "SeriesDescription": l["nom"],
            "n_coupes": s.get("n_coupes"), "images_par_coupe": l["ipc"],
            "heure": "%02d:%02d" % (t // 3600, (t % 3600) // 60) if t is not None else None,
            "type": dict(COLONNES_SEQUENCES).get(l["categorie"], l["categorie"] or "non classée"),
            "reconnu_par": l["methode"] if l["categorie"] else ""})
    return {"colonnes": colonnes, "non_classees": " | ".join(non_classees), "detail": detail}


def mp_post(texte):
    return _mots_presents(texte, "post_gado")


def mp_natif(texte):
    return _mots_presents(texte, "natif")


def ecrire_sequences(lignes, dossier):
    """sequences_realisees.csv (une ligne par patient, 0/1 par type) et
    sequences_detail.csv (une ligne par série). Les patients sans inventaire
    dans `lignes` (reprise d'un calcul) gardent leur ligne précédente."""
    dossier = Path(dossier)
    f1, f2 = dossier / "sequences_realisees.csv", dossier / "sequences_detail.csv"
    anciens1, anciens2 = {}, {}
    for f, d, groupe in ((f1, anciens1, False), (f2, anciens2, True)):
        if f.exists():
            try:
                df = pd.read_csv(f, sep=";", encoding="utf-8-sig", dtype=str)
                df = df.astype(object).where(pd.notna(df), None)
                for r in df.to_dict("records"):
                    if groupe:
                        d.setdefault(str(r["ID"]), []).append(r)
                    else:
                        d[str(r["ID"])] = r
            except Exception:
                pass
    r1, r2 = [], []
    for l in lignes:
        pid = str(l.get("ID"))
        inv = l.get("_sequences")
        if isinstance(inv, dict) and "colonnes" in inv:
            base = {"ID": pid}
            for k in ("code_pseudo", "Groupe", "date_examen", "identifiant_dicom"):
                if l.get(k) is not None:
                    base[k] = l[k]
            r1.append(dict(base, **{lib: inv["colonnes"].get(c, 0) for c, lib in COLONNES_SEQUENCES},
                           **{"series_non_classees": inv.get("non_classees") or ""}))
            r2 += [dict({"ID": pid, "identifiant_dicom": l.get("identifiant_dicom")}, **d)
                   for d in inv.get("detail", [])]
        else:
            if pid in anciens1:
                r1.append(anciens1[pid])
            r2 += anciens2.get(pid, [])
    if not r1:
        return None
    for f, r in ((f1, r1), (f2, r2)):
        tmp = f.with_name(f.stem + "_en_cours.tmp")
        ajouter_nom_prenom(pd.DataFrame(r)).to_csv(tmp, index=False, sep=";", encoding="utf-8-sig")
        try:
            tmp.replace(f)
        except PermissionError:
            tmp.replace(f.with_name(f.stem + "_secours.csv"))
    return f1


def choisir_serie(series, choix_manuel=None):
    """Choisit la série de perfusion de repos. Renvoie (serie, message)."""
    if choix_manuel is not None:
        for s in series.values():
            if s["serie_num"] == choix_manuel:
                return s, f"série {choix_manuel} imposée (choix_series.csv)"
        return None, f"série {choix_manuel} demandée mais introuvable"

    utilisables = [s for s in series.values() if not contient(s["description"], MOTS_DERIVES_EXCLUS)]
    nommees = [s for s in utilisables if contient(s["description"], MOTS_PERFUSION)
               and s["images_par_coupe"] >= MIN_IMAGES_PERF_NOMMEE]
    if nommees:
        perf = nommees
    elif not AUTORISER_SERIE_DE_REMPLACEMENT:
        if presence_perfusion(series):
            return None, ("série(s) de perfusion trouvée(s) mais trop courte(s) "
                          "(< %d images par coupe)" % MIN_IMAGES_PERF_NOMMEE)
        return None, ("aucune série au nom de perfusion (%s) : patient non analysé"
                      % "/".join(MOTS_PERFUSION))
    else:
        # aucune série nommée « perfusion » : séries dynamiques, sauf celles dont
        # le nom désigne autre chose (TI scout, cartographie, ciné, flux...)
        perf = [s for s in utilisables if s["images_par_coupe"] >= MIN_IMAGES_DYNAMIQUES
                and not contient(s["description"], MOTS_EXCLUS)]
        if not perf:
            return None, ("aucune série de perfusion reconnue (nom contenant %s) : "
                          "indiquer la série dans choix_series.csv" % "/".join(MOTS_PERFUSION))
    # doublons « ORIG X » / « X » (version d'origine et version recalée MoCo)
    def sans_orig(d):
        d = d.strip()
        return d[5:].strip().lower() if d.lower().startswith("orig ") else d.lower()
    origines = [s for s in perf if s["description"].strip().lower().startswith("orig ")]
    derivees = [s for s in perf if s not in origines]
    noms_derivees = {sans_orig(s["description"]) for s in derivees}
    if PREFERER_MOCO:
        perf = derivees + [s for s in origines if sans_orig(s["description"]) not in noms_derivees]
    else:
        noms_orig = {sans_orig(s["description"]) for s in origines}
        perf = origines + [s for s in derivees if sans_orig(s["description"]) not in noms_orig]
    note = ""
    if PREFERER_PERFUSION_PEU_DE_COUPES:
        # perfusion à 3-4 coupes (« 3 sl ») préférée aux perfusions à 8 coupes
        hors_stress = [s for s in perf if not contient(s["description"], MOTS_STRESS)]
        peu = [s for s in hors_stress if peu_de_coupes(s)]
        if peu and len(peu) < len(hors_stress):
            perf = peu + [s for s in perf if contient(s["description"], MOTS_STRESS)]
            note = " ; perfusion à %d coupes préférée aux perfusions à plus de coupes" % len(peu[0]["coupes"])
            if len(peu) == 1:
                n_autres = len(hors_stress) - 1
                return peu[0], "série à %d coupes préférée à %d perfusion%s à plus de coupes" % (
                    len(peu[0]["coupes"]), n_autres, "s" if n_autres > 1 else "")
    repos = [s for s in perf if contient(s["description"], MOTS_REPOS)
             and not contient(s["description"], MOTS_STRESS)]
    if len(repos) == 1:
        return repos[0], "série de repos identifiée par son nom" + note
    non_stress = [s for s in perf if not contient(s["description"], MOTS_STRESS)]
    candidats = repos or non_stress or perf
    if len(candidats) == 1:
        msg = "série unique" if len(perf) == 1 else "seule série non-stress"
        return candidats[0], msg + note
    cle = lambda s: (s["temps_min"] if s["temps_min"] is not None else s["serie_num"])
    candidats = sorted(candidats, key=cle)
    choisi = candidats[-1] if ORDRE_SI_AMBIGU == "derniere" else candidats[0]
    choisi["_candidats"] = candidats
    return choisi, "AMBIGU : %d séries possibles, %s retenue (à vérifier)" % (
        len(candidats), "dernière" if ORDRE_SI_AMBIGU == "derniere" else "première") + note


def peu_de_coupes(serie):
    """Perfusion à 3-4 coupes : nombre de coupes, ou « 3 sl » / « 4sl » dans le nom."""
    import re
    n = len(serie.get("coupes") or {}) or serie.get("n_coupes") or 0
    if 3 <= n <= N_COUPES_MAX_PREFEREE:
        return True
    return n <= N_COUPES_MAX_PREFEREE and bool(
        re.search(r"(?<!\d)[34]\s*(sl|slices?|coupes?)\b", str(serie.get("description", "")).lower()))


def rehaussement_serie(serie):
    """Intensité du rehaussement d'une série (coupe du milieu) : 99e percentile
    de (maximum dans le temps - ligne de base) / ligne de base. Une série sans
    injection (bolus raté, série de test) donne une valeur faible."""
    positions = sorted(serie["coupes"].keys())
    imgs = serie["coupes"][positions[len(positions) // 2]]
    pile, _, _ = charger_pixels(imgs)
    base = pile[:N_IMAGES_BASE].mean(axis=0)
    hausse = pile.max(axis=0) - base
    utiles = base > np.percentile(base, 50)          # tissus/sang, pas le fond
    if not utiles.any():
        return 0.0
    return float(np.percentile(hausse[utiles] / (base[utiles] + 1e-6), 99))


def _temps_valable(valeurs):
    """Une série de temps est utilisable si chaque image a sa valeur, si les
    valeurs sont (presque toutes) différentes et couvrent plus d'une seconde."""
    if any(v is None for v in valeurs) or len(valeurs) < 3:
        return False
    v = np.array(valeurs, dtype=float)
    return (np.ptp(v) > 1.0) and (len(np.unique(np.round(v, 3))) >= 0.8 * len(v))


def charger_pixels(liste_images):
    """Charge les pixels d'une coupe, triés dans le temps.
    Renvoie (pile, temps en s depuis la 1re image ou None, images triées).
    La source du temps retenue est notée dans images[0]["source_temps"].
    Sources essayées : AcquisitionTime, ContentTime, TriggerTime. Certaines
    séries dérivées (recalage MoCo) recopient la même heure sur toutes les
    images : elles sont alors rejetées et le temps sera reconstruit."""
    def lire_pixels(f):
        ds = pydicom.dcmread(f, force=True)
        arr = ds.pixel_array.astype(np.float32)
        pente = float(_get(ds, "RescaleSlope", 1) or 1)
        ord0 = float(_get(ds, "RescaleIntercept", 0) or 0)
        return arr * pente + ord0

    uniques = list(dict.fromkeys(im["fichier"] for im in liste_images))
    if N_FILS_LECTURE > 1 and len(uniques) >= 20:   # plusieurs fichiers lus à la fois
        from concurrent.futures import ThreadPoolExecutor
        with ThreadPoolExecutor(N_FILS_LECTURE) as ex:
            cache = dict(zip(uniques, ex.map(lire_pixels, uniques)))
    else:
        cache = {f: lire_pixels(f) for f in uniques}
    paires = []
    for im in liste_images:
        arr = cache[im["fichier"]]
        a = arr[im["trame"]] if im["trame"] is not None else arr
        paires.append(((im["serie_num"], im["instance"]), a, im))
    paires.sort(key=lambda p: p[0])

    def heures(cle, facteur=1.0):
        vals = []
        for _, _, im in paires:
            v = im.get(cle)
            try:
                vals.append(float(v) * facteur if v is not None else None)
            except (TypeError, ValueError):
                vals.append(None)
        return vals

    sources = [("AcquisitionTime", heures("temps_acq")),
               ("ContentTime", heures("temps_contenu")),
               ("TriggerTime", heures("trigger_ms", 0.001))]
    if all(p[2].get("temps_acq") is None and p[2].get("temps_contenu") is None for p in paires):
        sources.insert(0, ("temps de l'image", heures("temps")))   # DICOM Enhanced
    temps, source = None, "aucune heure par image utilisable"
    for nom, vals in sources:
        if _temps_valable(vals):
            v = np.array(vals, dtype=float)
            v = np.where(v < v[0] - 3600, v + 86400, v)            # passage de minuit
            ordre = np.argsort(v, kind="stable")
            paires = [paires[i] for i in ordre]
            temps = v[ordre] - v[ordre][0]
            source = nom
            break
    pile = np.stack([p[1] for p in paires], axis=0)
    imgs = [p[2] for p in paires]
    imgs[0]["source_temps"] = source
    return pile, temps, imgs


# =============================================================================
# 3. DÉTECTION DES CAVITÉS ET COURBES
# =============================================================================


def _carte_rehaussement(pile):
    """Pile légèrement lissée (médiane sur 3 images dans le temps : supprime les
    images isolées bougées ; gaussienne dans l'espace) et rehaussement par pixel.
    Renvoie (rehausse, rmax, coherent) ; coherent = pixels dont la courbe varie
    de façon régulière (bolus) et non en dents de scie (mouvement, paroi)."""
    lisse = ndimage.median_filter(pile, size=(3, 1, 1)) if pile.shape[0] >= 5 else pile
    lisse = ndimage.gaussian_filter(lisse, sigma=(0, 1, 1))
    base = lisse[:N_IMAGES_BASE].mean(axis=0)
    rehausse = lisse - base[None]
    rmax = rehausse.max(axis=0)
    brut = ndimage.gaussian_filter(pile, sigma=(0, 1, 1))
    saut = np.percentile(np.abs(np.diff(brut, axis=0)), 75, axis=0)
    coherent = rmax > 3.0 * (saut + 1e-6)
    return rehausse, rmax, coherent


def fenetre_premier_passage(pile, temps):
    """Indice de fin (exclu) de la fenêtre du premier passage et heure d'arrivée
    du contraste. Arrivée = première image où au moins 20 pixels nettement
    rehaussés atteignent 30 % de leur maximum."""
    n = pile.shape[0]
    if n < 8:
        return n, None
    rehausse, rmax, coherent = _carte_rehaussement(pile)
    forts = coherent & (rmax > 0.2 * np.percentile(rmax[coherent], 99.5)) if coherent.any() else None
    if forts is None or forts.sum() < 20:
        return n, None
    arrivee = np.argmax(rehausse[:, forts] >= 0.3 * rmax[forts][None], axis=0)
    compte = np.bincount(arrivee, minlength=n).cumsum()
    i_a = int(np.argmax(compte >= 20))
    t_a = float(temps[i_a])
    fin = int(np.searchsorted(temps, t_a + FENETRE_PREMIER_PASSAGE_S, side="right"))
    fin = min(n, max(fin, i_a + 8, N_IMAGES_BASE + 8))
    return fin, t_a


def detecter_cavites(pile, temps, espacement_mm):
    """Détecte VD (rehaussement précoce) et VG (tardif). Renvoie dict ou None.
    La pile doit déjà être limitée à la fenêtre du premier passage."""
    rehausse, rmax, coherent = _carte_rehaussement(pile)
    i_ttp = np.argmax(rehausse, axis=0)
    ttp = temps[i_ttp]

    seuil = 0.35 * np.percentile(rmax[coherent] if coherent.sum() >= 20 else rmax, 99.5)
    # pixel retenu : nettement rehaussé, courbe régulière, et pic AVANT la fin de
    # la fenêtre (un pixel encore en montée à la fin = organe tardif, pas une cavité)
    cand = (rmax > seuil) & coherent & (i_ttp < pile.shape[0] - 2)
    cand = ndimage.binary_opening(cand, iterations=1)
    if cand.sum() < 20:
        return None, "trop peu de pixels rehaussés"

    # séparation en 2 groupes selon le temps au pic (k-moyennes 1D simple)
    t_c = ttp[cand]
    c1, c2 = np.percentile(t_c, 20), np.percentile(t_c, 80)
    for _ in range(50):
        g1 = np.abs(t_c - c1) <= np.abs(t_c - c2)
        n1, n2 = t_c[g1], t_c[~g1]
        if len(n1) == 0 or len(n2) == 0:
            break
        nc1, nc2 = n1.mean(), n2.mean()
        if abs(nc1 - c1) < 1e-3 and abs(nc2 - c2) < 1e-3:
            break
        c1, c2 = nc1, nc2
    if abs(c2 - c1) < 1.5:
        return None, "pas de séparation nette précoce/tardif (écart %.1f s)" % abs(c2 - c1)
    precoce_t, tardif_t = min(c1, c2), max(c1, c2)
    limite = (precoce_t + tardif_t) / 2

    try:
        dy, dx = (float(espacement_mm[0]), float(espacement_mm[1]))
    except Exception:
        dy = dx = 1.5

    def composantes(masque, k=4):
        lab, n = ndimage.label(masque)
        if n == 0:
            return []
        poids = ndimage.sum(rmax, lab, index=range(1, n + 1))
        taille = ndimage.sum(np.ones_like(rmax), lab, index=range(1, n + 1))
        res = []
        for j in np.argsort(poids)[::-1][:k]:
            if taille[j] < 10:
                continue
            comp = ndimage.binary_fill_holes(lab == j + 1)
            res.append((comp, float(poids[j]), ndimage.center_of_mass(comp)))
        return res

    def noyau(comp):
        # ROI ronde : disque centré au point le plus profond de la cavité (le plus
        # loin des parois), calculé sur le sang seulement (piliers, trabéculations
        # et myocarde peu rehaussés exclus), donc loin de tout volume partiel.
        sang = comp & (rmax >= FRACTION_REHAUSSEMENT_ROI * np.percentile(rmax[comp], 90))
        zone = sang if sang.sum() >= 12 else comp
        dist = ndimage.distance_transform_edt(zone)
        dist_l = ndimage.gaussian_filter(dist, 1.0) * zone     # centre stable si plateau
        cy, cx = np.unravel_index(int(np.argmax(dist_l)), dist.shape)
        r = max(1.5, (1.0 - FRACTION_NOYAU_ROI) * float(dist.max()))
        return comp, zone, cy, cx, r

    def disque_commun(c, r):
        comp, zone, cy, cx, _ = c
        disque = masque_disque(comp.shape, cy, cx, r) & zone
        return disque if disque.sum() >= 6 else comp

    cvd = composantes(cand & (ttp <= limite))
    cvg = composantes(cand & (ttp > limite))
    if not cvd or not cvg:
        return None, "cavité VD ou VG non trouvée"
    # meilleur couple : deux cavités bien rehaussées ET voisines
    meilleur = None
    for a in cvd:
        for b in cvg:
            d = math.hypot((a[2][0] - b[2][0]) * dy, (a[2][1] - b[2][1]) * dx)
            proximite = 1.0 if d <= DISTANCE_CAVITES_OK_MM else math.exp(-(d - DISTANCE_CAVITES_OK_MM) / 25.0)
            score = math.sqrt(a[1] * b[1]) * proximite
            if meilleur is None or score > meilleur[0]:
                meilleur = (score, a[0], b[0], d)
    _, vd, vg, dist = meilleur
    nvd, nvg = noyau(vd), noyau(vg)
    r = max(1.5, min(nvd[4], nvg[4]) * FACTEUR_TAILLE_ROI)     # même taille VD et VG
    vd, vg = disque_commun(nvd, r), disque_commun(nvg, r)
    return {"vd": vd, "vg": vg, "distance_mm": dist, "rayon_px": r, "centres": ((nvd[2], nvd[3]), (nvg[2], nvg[3])),
            "aire_vd": int(vd.sum()), "aire_vg": int(vg.sum()),
            "ecart_groupes_s": tardif_t - precoce_t}, "ok"


def masque_disque(forme, cy, cx, r):
    yy, xx = np.ogrid[:forme[0], :forme[1]]
    return (yy - float(cy)) ** 2 + (xx - float(cx)) ** 2 <= float(r) ** 2


def cercle_equivalent(masque):
    """(ligne du centre, colonne du centre, rayon) du disque de même aire."""
    cy, cx = ndimage.center_of_mass(masque)
    return [round(float(cy), 1), round(float(cx), 1), round(math.sqrt(masque.sum() / math.pi), 1)]


def _cercle_auto(cav, i, cle):
    """Cercle (ligne, colonne, rayon) d'une ROI automatique pour l'éditeur."""
    if cav.get("centres") and cav.get("rayon_px"):
        cy, cx = cav["centres"][i]
        return [round(float(cy), 1), round(float(cx), 1), round(float(cav["rayon_px"]), 1)]
    return cercle_equivalent(cav[cle])


def cavites_manuelles(roi, forme, espacement_mm):
    """ROI dessinées à la main : roi = {"vd": [ligne, colonne, rayon], "vg": [...]}
    en pixels de l'image."""
    vd, vg = masque_disque(forme, *roi["vd"]), masque_disque(forme, *roi["vg"])
    if vd.sum() < 3 or vg.sum() < 3:
        return None, "ROI manuelle trop petite (moins de 3 pixels)"
    try:
        dy, dx = (float(espacement_mm[0]), float(espacement_mm[1]))
    except Exception:
        dy = dx = 1.5
    dist = math.hypot((roi["vd"][0] - roi["vg"][0]) * dy, (roi["vd"][1] - roi["vg"][1]) * dx)
    return {"vd": vd, "vg": vg, "distance_mm": dist, "aire_vd": int(vd.sum()),
            "aire_vg": int(vg.sum()), "ecart_groupes_s": None}, "ok"


def gamma(t, ym, t0, tm, a):
    """Gamma-variate normalisée : pic ym au temps t0 + tm."""
    x = np.clip((t - t0) / tm, 1e-9, None)
    y = ym * x ** a * np.exp(a * (1 - x))
    return np.where(t > t0, y, 0.0)


def analyser_courbe(t, y):
    """Centre de gravité et pic d'une courbe de premier passage."""
    res = {"ok": False}
    if len(y) < 8 or np.nanmax(y) <= 0:
        return res
    i_pic = int(np.argmax(y))
    ymax = float(y[i_pic])
    # pic brut affiné par interpolation parabolique
    if 0 < i_pic < len(y) - 1:
        y0, y1, y2 = y[i_pic - 1], y[i_pic], y[i_pic + 1]
        den = (y0 - 2 * y1 + y2)
        dec = 0.5 * (y0 - y2) / den if den != 0 else 0.0
        dec = float(np.clip(dec, -0.5, 0.5))
        dt = (t[i_pic + 1] - t[i_pic - 1]) / 2
        t_pic = float(t[i_pic] + dec * dt)
    else:
        t_pic = float(t[i_pic])
    # sommet écrêté (saturation) : on prend le milieu du plateau >= 95 % du max
    haut = y >= 0.95 * ymax
    if haut.sum() >= 2:
        g = d = i_pic
        while g - 1 >= 0 and haut[g - 1]:
            g -= 1
        while d + 1 < len(y) and haut[d + 1]:
            d += 1
        if d > g:
            t_pic = float((t[g] + t[d]) / 2)
    res["t_pic_brut"] = t_pic

    # début : dernier point < 10 % du maximum avant le pic
    debut = 0
    for i in range(i_pic, -1, -1):
        if y[i] < 0.10 * ymax:
            debut = i
            break
    # fin : premier point < 50 % après le pic, ou remontée (recirculation)
    fin = len(y) - 1
    for i in range(i_pic + 1, len(y)):
        if y[i] < 0.50 * ymax or (i > i_pic + 1 and y[i] > y[i - 1] * 1.05):
            fin = i
            break
    # et pas plus loin que 1,5 fois la durée de montée après le pic : sans lavage
    # net (plateau de recirculation), la queue fausserait le centre de gravité
    t_lim = t_pic + max(1.5 * (t_pic - float(t[debut])), 3.0)
    fin = min(fin, max(i for i in range(len(t)) if t[i] <= t_lim) if t[0] <= t_lim else fin)
    fin = max(fin, min(i_pic + 2, len(y) - 1))
    fen = slice(debut, fin + 1)
    tf, yf = t[fen], y[fen]

    # plateau (saturation) : nombre d'images >= 95 % du maximum
    res["n_plateau"] = int(np.sum(y >= 0.95 * ymax))

    # En cas de saturation (sommet écrêté), on ajuste la gamma-variate sur les
    # flancs montant et descendant seulement : la courbe reconstruit le sommet.
    if res["n_plateau"] >= 2:
        garde = yf < 0.95 * ymax
        if garde.sum() >= 5:
            tf, yf = tf[garde], yf[garde]
            res["sommet_exclu"] = True

    try:
        p0 = [ymax, float(t[debut]) - 0.5, max(t_pic - float(t[debut]), 0.5), 2.0]
        bornes = ([0.3 * ymax, float(t[0]) - 5, 0.2, 0.3],
                  [3.0 * ymax, t_pic - 0.05, 15.0, 15.0])
        p0[1] = float(np.clip(p0[1], bornes[0][1] + 1e-3, bornes[1][1] - 1e-3))
        p0[2] = float(np.clip(p0[2], bornes[0][2] + 1e-3, bornes[1][2] - 1e-3))
        popt, _ = curve_fit(gamma, tf, yf, p0=p0, bounds=bornes, maxfev=20000)
        ym, t0, tm, a = popt
        pred = gamma(tf, *popt)
        ss_res = float(np.sum((yf - pred) ** 2))
        ss_tot = float(np.sum((yf - yf.mean()) ** 2)) or 1e-9
        res.update({"ok": True, "ajustement": "gamma",
                    "t_centroide": float(t0 + tm * (a + 1) / a),
                    "t_pic_ajuste": float(t0 + tm),
                    "r2": 1 - ss_res / ss_tot, "params": popt})
    except Exception:
        w = np.clip(yf, 0, None)
        if w.sum() > 0:
            res.update({"ok": True, "ajustement": "brut",
                        "t_centroide": float(np.sum(tf * w) / np.sum(w)),
                        "t_pic_ajuste": t_pic, "r2": float("nan"), "params": None})
    res["fenetre"] = (debut, fin)
    res["lmh"] = largeur_mi_hauteur(t, y, res)
    # sommet écrêté (saturation) : plateau long par rapport à la largeur de la
    # courbe, ou gamma ajustée nettement au-dessus du maximum mesuré. Une courbe
    # simplement large a aussi plusieurs images >= 95 % sans être saturée : sur
    # une gamma, la durée >= 95 % vaut environ 0,25 x la largeur à mi-hauteur.
    p, dt = res.get("params"), float(np.median(np.diff(t)))
    lmh = res["lmh"]
    res["ecrete"] = bool(res["n_plateau"] >= 3 and (
        (lmh and (res["n_plateau"] - 1) * dt / lmh >= 0.35) or (p is not None and p[0] > 1.15 * ymax)))
    return res


def largeur_mi_hauteur(t, y, res=None):
    """Largeur à mi-hauteur (s) mesurée sur la courbe elle-même (interpolation
    linéaire), dans la fenêtre analysée. None si la courbe ne redescend pas sous
    50 % du maximum (plateau, recirculation, série trop courte) : la gamma
    ajustée n'est pas utilisée car, au-delà des points mesurés, sa queue est
    extrapolée."""
    i = int(np.argmax(y)); demi = y[i] / 2
    g = i
    while g > 0 and y[g - 1] >= demi:
        g -= 1
    d = i
    while d < len(y) - 1 and y[d + 1] >= demi:
        d += 1
    if g == 0 or d == len(y) - 1:
        return None
    tg = t[g - 1] + (demi - y[g - 1]) / (y[g] - y[g - 1]) * (t[g] - t[g - 1])
    td = t[d] + (y[d] - demi) / (y[d] - y[d + 1]) * (t[d + 1] - t[d])
    return float(td - tg)


def bilan_dispersion(a_vd, a_vg, ptt_c, ptt_p):
    """Largeurs à mi-hauteur, rapport VG/VD et alerte d'écart centroïde-pic
    qualifiée : courbes propres + VG étalé = dispersion physiologique probable ;
    courbes de mauvaise qualité = erreur de mesure probable."""
    l_vd, l_vg = a_vd.get("lmh"), a_vg.get("lmh")
    rapport = l_vg / l_vd if l_vd and l_vg else None
    defauts = []
    for nom, a in (("VD", a_vd), ("VG", a_vg)):
        if a.get("ajustement") != "gamma" or not a.get("r2", 0) >= R2_COURBE_PROPRE:
            defauts.append("ajustement %s médiocre" % nom)
        if a.get("ecrete"):
            defauts.append("sommet %s écrêté" % nom)
        if not a.get("lmh"):
            defauts.append("%s ne redescend pas sous 50 %%" % nom)
    propres = not defauts
    etale = bool(propres and rapport >= RAPPORT_LMH_ETALE)
    cols = {"LMH_VD_s": round(l_vd, 2) if l_vd else None,
            "LMH_VG_s": round(l_vg, 2) if l_vg else None,
            "rapport_LMH_VG_VD": round(rapport, 2) if rapport else None,
            "ecart_centroide_pic_s": round(ptt_c - ptt_p, 2),
            "courbes_propres": int(propres),
            "dispersion_VG": "VG étalé (courbes propres)" if etale else ""}
    alertes = []
    ecart = ptt_c - ptt_p
    if abs(ecart) > ECART_PIC_CENTROIDE_MAX_S:
        txt = "écart centroïde-pic %+.1f s" % ecart
        if etale and ecart > 0:
            txt += " : courbes propres et VG étalé (rapport LMH %.1f), dispersion physiologique probable" % rapport
        elif not propres:
            txt += " : %s, erreur de mesure probable" % ", ".join(defauts)
        else:
            txt += " : courbes propres sans VG étalé (rapport LMH %.1f), vérifier les ROI" % rapport
        alertes.append(txt)
    return cols, alertes


# =============================================================================
# 4. ANALYSE D'UN PATIENT
# =============================================================================


def rr_de_la_serie(images_coupe, temps):
    """Intervalle RR (ms) et sa source."""
    vals = [im["rr_ms_tag"] for im in images_coupe if im.get("rr_ms_tag")]
    if vals:
        return float(np.median(vals)), "en-tête DICOM (NominalInterval)"
    fcs = [im["fc_tag"] for im in images_coupe if im.get("fc_tag")]
    if fcs:
        return 60000.0 / float(np.median(fcs)), "en-tête DICOM (HeartRate)"
    if temps is not None and len(temps) > 3:
        d = np.median(np.diff(temps)) * 1000
        if 300 <= d <= 2000:
            return float(d), "intervalle entre images (1 image/battement supposée)"
    return None, "indisponible"


def analyser_patient(id_patient, dossier, clin, choix, dossier_png):
    ligne = {"ID": id_patient, "statut": "échec", "message": ""}
    # série imposée : lecture complète (elle peut être dans un dossier sauté)
    images = lire_entetes_patient(dossier, complet=choix.get(id_patient, {}).get("serie") is not None)
    if not images:
        ligne["message"] = "aucune image DICOM lisible"
        return ligne
    series = resumer_series(images)
    ligne["presence_sequence_perfusion"] = presence_perfusion(series)
    ligne.update(informations_examen(images))
    try:
        ligne["_sequences"] = inventaire_sequences(series)
    except Exception as e:
        ligne["_sequences"] = {"colonnes": {}, "non_classees": "inventaire impossible : %s" % e, "detail": []}
    choix_serie = choix.get(id_patient, {}).get("serie")
    choix_coupe = choix.get(id_patient, {}).get("coupe")
    roi_man = choix.get(id_patient, {}).get("roi")      # ROI dessinées à la main (site)
    fen_man = choix.get(id_patient, {}).get("fenetre")  # [début, fin] (s) choisis à la main (site)
    if (not ligne["presence_sequence_perfusion"] and choix_serie is None
            and not AUTORISER_SERIE_DE_REMPLACEMENT):
        # pas de perfusion nommée : on passe au patient suivant, sans pixel lu
        ligne["statut"] = "sans perfusion"
        ligne["message"] = ligne["choix_serie"] = (
            "aucune série au nom de perfusion (%s) : patient non analysé" % "/".join(MOTS_PERFUSION))
        return ligne
    serie, msg = choisir_serie(series, choix_serie)
    ligne["choix_serie"] = msg
    if serie is None:
        ligne["message"] = msg
        return ligne
    # Plusieurs perfusions possibles : on écarte celles sans rehaussement
    # (injection ratée ou série de test), puis on applique la règle habituelle.
    cands = serie.pop("_candidats", None)
    if cands and len(cands) > 1:
        try:
            rh = {id(c): rehaussement_serie(c) for c in cands}
            meilleur_rh = max(rh.values())
            gardes = [c for c in cands if rh[id(c)] >= SEUIL_REHAUSSEMENT_RELATIF * meilleur_rh]
            ecartees = [str(c["serie_num"]) for c in cands if c not in gardes]
            if ecartees:
                serie = gardes[-1] if ORDRE_SI_AMBIGU == "derniere" else gardes[0]
                msg = "série %d retenue ; série(s) %s écartée(s) : peu ou pas de rehaussement" % (
                    serie["serie_num"], ", ".join(ecartees))
                if len(gardes) > 1:
                    msg += " ; AMBIGU entre %d séries rehaussées (à vérifier)" % len(gardes)
                ligne["choix_serie"] = msg
        except Exception as e:
            ligne["choix_serie"] = msg + " (contrôle du rehaussement impossible : %s)" % e
    ligne["serie_num"] = serie["serie_num"]
    ligne["serie_description"] = serie["description"]
    if serie.get("fusion"):
        ligne["choix_serie"] = (ligne.get("choix_serie") or "") + " ; " + serie["fusion"]

    positions = sorted(serie["coupes"].keys())
    ordre = list(range(len(positions)))
    coupe_milieu = None
    if PREFERER_PERFUSION_PEU_DE_COUPES and not choix_coupe and 3 <= len(positions) <= N_COUPES_MAX_PREFEREE:
        # perfusion à 3-4 coupes : la coupe du milieu (2 sur 3) est analysée en
        # premier et retenue dès qu'elle est exploitable
        coupe_milieu = (len(positions) - 1) // 2
        ordre = [coupe_milieu] + [k for k in ordre if k != coupe_milieu]
    if choix_coupe:
        if 1 <= choix_coupe <= len(positions):
            ordre = [choix_coupe - 1]
        else:
            ligne["message"] = "coupe %s demandée, %d disponibles" % (choix_coupe, len(positions))
            return ligne

    import re as _re
    n_rr = RR_PAR_IMAGE
    if n_rr is None:
        n_rr = 1
        for s_ in series.values():
            d_ = s_["description"].lower()
            m_rr = _re.search(r"(\d)\s*rr", d_)
            if m_rr and contient(d_, MOTS_PERFUSION):
                n_rr = int(m_rr.group(1))
                break
    # série d'origine jumelle (« ORIG <même nom> ») d'une série recalée
    jumelle = None
    nom_serie = serie["description"].strip().lower()
    for s_ in series.values():
        d_ = s_["description"].strip().lower()
        if d_.startswith("orig ") and d_[5:].strip() == nom_serie:
            if jumelle is None or abs(s_["serie_num"] - 40000 - serie["serie_num"]) == 0:
                jumelle = s_

    meilleur = None
    raisons = []
    for k in ordre:
        imgs = serie["coupes"][positions[k]]
        if len(imgs) < MIN_IMAGES_PERF_NOMMEE:
            raisons.append("coupe %d : %d images seulement" % (k + 1, len(imgs)))
            continue
        pile, temps, imgs_triees = charger_pixels(imgs)
        temps_reconstruit = False
        if temps is None and jumelle is not None:
            # série recalée (MoCo) sans heure par image : on reprend les heures
            # de la série d'origine (« ORIG ... »), coupe par coupe, dans l'ordre
            pos_j = sorted(jumelle["coupes"].keys())
            if len(pos_j) == len(positions):
                imgs_j = sorted(jumelle["coupes"][pos_j[k]], key=lambda im: im["instance"])
                if len(imgs_j) == pile.shape[0]:
                    for nom, cle, f in (("AcquisitionTime", "temps_acq", 1.0),
                                        ("ContentTime", "temps_contenu", 1.0)):
                        vals = [im.get(cle) for im in imgs_j]
                        if _temps_valable(vals):
                            v = np.array(vals, dtype=float)
                            v = np.where(v < v[0] - 3600, v + 86400, v)
                            if np.all(np.diff(v) > 0):
                                temps = v - v[0]
                                imgs_triees[0]["source_temps"] = ("%s de la série d'origine %d"
                                                                  % (nom, jumelle["serie_num"]))
                                break
        if temps is None:
            # Pas d'heure par image : si le RR est connu, on suppose 1 image par
            # battement (cas habituel des perfusions) et on reconstruit le temps.
            rr_tmp, _ = rr_de_la_serie(imgs_triees, None)
            if rr_tmp is None and AUTORISER_TEMPS_RECONSTRUIT:
                c_ = clin.get(id_patient, {})
                rr_tmp = c_.get("RR") or (60000.0 / c_["FC"] if c_.get("FC") else None)
            if rr_tmp is None or not AUTORISER_TEMPS_RECONSTRUIT:
                raisons.append("coupe %d : pas d'heure utilisable par image et RR inconnu "
                               "(renseigner FC ou RR_ms dans cliniques.xlsx)" % (k + 1))
                continue
            temps = np.arange(pile.shape[0]) * rr_tmp * n_rr / 1000.0
            temps_reconstruit = True
            imgs_triees[0]["source_temps"] = ("reconstruit : 1 image tous les %d RR, RR = %.0f ms"
                                              % (n_rr, rr_tmp))
        # premier passage seulement (2e injection, recirculation, organes tardifs exclus)
        fin, t_arr = fenetre_premier_passage(pile, temps)
        fin_auto, deb = fin, 0
        if fen_man:     # fenêtre d'analyse des courbes choisie à la main
            deb, fin = indices_fenetre(temps, fen_man, fin)
        if roi_man:
            cav, m = cavites_manuelles(roi_man, pile.shape[1:], imgs_triees[0]["espacement"])
        else:           # ROI automatiques : toujours placées sur le premier passage entier
            cav, m = detecter_cavites(pile[:fin_auto], temps[:fin_auto], imgs_triees[0]["espacement"])
        if cav is None:
            raisons.append("coupe %d : %s" % (k + 1, m))
            continue
        y_vd = (pile[:, cav["vd"]].mean(axis=1) - pile[:N_IMAGES_BASE, cav["vd"]].mean())
        y_vg = (pile[:, cav["vg"]].mean(axis=1) - pile[:N_IMAGES_BASE, cav["vg"]].mean())
        a_vd, a_vg = analyser_courbe(temps[deb:fin], y_vd[deb:fin]), analyser_courbe(temps[deb:fin], y_vg[deb:fin])
        if not (a_vd["ok"] and a_vg["ok"]):
            raisons.append("coupe %d : courbe non analysable" % (k + 1))
            continue
        ptt = a_vg["t_pic_brut"] - a_vd["t_pic_brut"]
        plaus = PTT_MIN_S <= ptt <= PTT_MAX_S
        r2_min = min(a_vd.get("r2", 0), a_vg.get("r2", 0))
        r2_min = 0.0 if r2_min != r2_min else r2_min           # nan -> 0
        score = (1 if plaus else 0, 1 if cav["distance_mm"] <= 120 else 0,
                 1 if r2_min >= 0.8 else 0, min(cav["aire_vd"], cav["aire_vg"]))
        cand = dict(k=k, pile=pile, temps=temps, imgs=imgs_triees, cav=cav,
                    temps_reconstruit=temps_reconstruit, fin=fin, debut=deb, t_arrivee=t_arr,
                    y_vd=y_vd, y_vg=y_vg, a_vd=a_vd, a_vg=a_vg, score=score)
        if meilleur is None or cand["score"] > meilleur["score"]:
            meilleur = cand
        if k == coupe_milieu:
            break       # coupe du milieu exploitable : c'est elle qu'on garde
    if coupe_milieu is not None and meilleur is not None:
        if meilleur["k"] == coupe_milieu:
            ligne["choix_serie"] = (ligne.get("choix_serie") or "") + " ; coupe du milieu (%d/%d)" % (
                coupe_milieu + 1, len(positions))
        else:
            ligne["choix_serie"] = (ligne.get("choix_serie") or "") + (
                " ; coupe du milieu (%d/%d) inexploitable, coupe %d retenue" % (
                    coupe_milieu + 1, len(positions), meilleur["k"] + 1))
    if meilleur is None:
        ligne["message"] = " ; ".join(raisons) or "aucune coupe exploitable"
        return ligne

    b = meilleur
    rr_ms, src_rr = rr_de_la_serie(b["imgs"], b["temps"])
    if rr_ms and src_rr.startswith("intervalle") and n_rr > 1:
        # série acquise tous les n battements (« Perf 2RR ») : 1 image = n RR
        rr_ms, src_rr = rr_ms / n_rr, "intervalle entre images / %d (série « %dRR »)" % (n_rr, n_rr)
    fc_clin = clin.get(id_patient, {}).get("FC")
    if rr_ms is None and fc_clin:
        rr_ms, src_rr = 60000.0 / fc_clin, "FC du fichier clinique"
    a_vd, a_vg = b["a_vd"], b["a_vg"]
    ptt_c = a_vg["t_centroide"] - a_vd["t_centroide"]
    ptt_p = a_vg["t_pic_brut"] - a_vd["t_pic_brut"]

    alertes = []
    ptt_princ = ptt_p if METHODE_PRINCIPALE == "pic" else ptt_c
    if not (PTT_MIN_S <= ptt_princ <= PTT_MAX_S):
        alertes.append("PTT hors bornes")
    ptt_autre, nom_autre = (ptt_c, "centroïde") if METHODE_PRINCIPALE == "pic" else (ptt_p, "pic")
    if not (PTT_MIN_S <= ptt_autre <= PTT_MAX_S):
        alertes.append("PTT %s hors bornes (%.1f s)" % (nom_autre, ptt_autre))
    if a_vg.get("ecrete"):
        alertes.append("plateau VG (saturation probable)")
    if a_vd.get("ecrete"):
        alertes.append("plateau VD (saturation probable)")
    disp, al = bilan_dispersion(a_vd, a_vg, ptt_c, ptt_p)
    alertes += al
    if b["cav"]["distance_mm"] > 120:
        alertes.append("cavités éloignées (%.0f mm)" % b["cav"]["distance_mm"])
    for nom, a in (("VD", a_vd), ("VG", a_vg)):
        if a.get("ajustement") == "brut":
            alertes.append("ajustement gamma échoué %s" % nom)
        elif a.get("r2", 1) < 0.90:
            alertes.append("ajustement %s médiocre (R2=%.2f)" % (nom, a["r2"]))
    if rr_ms is None:
        alertes.append("RR inconnu")
    elif src_rr.startswith("intervalle entre images (") and rr_ms > INTERVALLE_SUSPECT_MS:
        alertes.append("intervalle entre images %.0f ms : 1 image tous les 2 battements ? "
                       "renseigner RR_ms ou FC dans cliniques.xlsx" % rr_ms)
    if "AMBIGU" in msg:
        alertes.append("choix de série ambigu")
    if b.get("temps_reconstruit"):
        alertes.append("temps reconstruit à partir du RR (heures DICOM inutilisables)")

    ligne.update({
        "statut": "ok" if not alertes else "à vérifier",
        "coupe_retenue": b["k"] + 1,
        "n_coupes_serie": len(positions),
        "n_images": len(b["temps"]),
        "ROI": "manuelles" if roi_man else "automatiques",
        "fenetre_premier_passage_s": round(float(b["temps"][b["fin"] - 1]), 1),
        "fenetre_analyse_s": "%.1f-%.1f" % (b["temps"][b.get("debut", 0)], b["temps"][b["fin"] - 1]),
        "fenetre_manuelle": "oui" if fen_man else "",
        "intervalle_images_s": round(float(np.median(np.diff(b["temps"]))), 3),
        "source_temps": b["imgs"][0].get("source_temps"),
        "aire_roi_VD_px": b["cav"]["aire_vd"],
        "aire_roi_VG_px": b["cav"]["aire_vg"],
        "RR_ms": round(rr_ms, 1) if rr_ms else None,
        "source_RR": src_rr,
        "PTT_centroide_s": round(ptt_c, 2),
        "PTT_pic_s": round(ptt_p, 2),
        "R2_gamma_VD": round(a_vd.get("r2", float("nan")), 3),
        "R2_gamma_VG": round(a_vg.get("r2", float("nan")), 3),
        "plateau_VG_images": a_vg["n_plateau"],
        "alertes": " ; ".join(alertes),
    })
    ligne.update(disp)
    if rr_ms:
        rr_s = rr_ms / 1000
        ligne["FC_bpm"] = round(60.0 / rr_s, 1)
        for suf, ptt in (("centroide", ptt_c), ("pic", ptt_p)):
            # PTTc = PTT / RR = PTT x FC / 60 = nombre de cycles cardiaques de transit
            ligne["PTTc_%s" % suf] = round(ptt / rr_s, 2)
            # PTT corrigé de la FC (formule de Bazett, comme le QTc) : PTT / racine(RR en s)
            ligne["PTT_bazett_%s" % suf] = round(ptt / math.sqrt(rr_s), 2)

    # PBVi (surface corporelle du fichier clinique, sinon celle des en-têtes DICOM)
    c = clin.get(id_patient, {})
    sc = c.get("SC") or ligne.get("SC_dicom_m2")
    if sc:
        ligne["SC_m2"] = round(sc, 3)
        ligne["source_SC"] = "fichier clinique" if c.get("SC") else "en-têtes DICOM"
        fc_ref = c.get("FC") or (60000.0 / rr_ms if rr_ms else None)
        for suf, ptt in (("centroide", ptt_c), ("pic", ptt_p)):
            if c.get("VES_VD") and rr_ms:
                ligne["PBVi_VD_%s_mL_m2" % suf] = round((ptt / (rr_ms / 1000)) * c["VES_VD"] / sc, 1)
            if c.get("VES_VG") and fc_ref:
                ligne["PBVi_VG_%s_mL_m2" % suf] = round(ptt * (c["VES_VG"] * fc_ref / 60.0) / sc, 1)

    ligne["PTT_principal_s"] = ligne["PTT_centroide_s"] if METHODE_PRINCIPALE == "centroide" else ligne["PTT_pic_s"]
    # géométrie des ROI (pour l'éditeur du site ; non écrite dans le CSV)
    ligne["_roi"] = {"serie": serie["serie_num"], "coupe": b["k"] + 1,
                     "image": b.get("debut", 0) + int(np.argmax(b["y_vg"][b.get("debut", 0):b["fin"]])),
                     "fenetre": [round(float(b["temps"][b.get("debut", 0)]), 2), round(float(b["temps"][b["fin"] - 1]), 2)]
                     if fen_man else None,
                     "vd": list(roi_man["vd"]) if roi_man else _cercle_auto(b["cav"], 0, "vd"),
                     "vg": list(roi_man["vg"]) if roi_man else _cercle_auto(b["cav"], 1, "vg")}
    b["roi_manuelles"] = bool(roi_man)
    if GARDER_SERIE_RETENUE:
        ligne["_serie"] = serie
    ligne["message"] = msg
    if SAUVER_COURBES:
        try:
            ligne["fichier_courbes"] = ecrire_courbes(id_patient, b, rr_ms, src_rr, ligne, roi_man,
                                                      Path(dossier_png).parent / "courbes")
        except Exception as e:
            ligne["message"] += " ; courbes non sauvegardées (%s)" % e
    try:
        dessiner_controle(id_patient, b, ligne, dossier_png)
    except Exception as e:
        ligne["message"] += " ; image de contrôle non créée (%s)" % e
    return ligne


def nom_fichier_sur(texte):
    import re
    return re.sub(r'[<>:"/\\|?*\x00-\x1f]', "_", str(texte)).strip(" .") or "patient"


def ecrire_courbes(id_patient, b, rr_ms, src_rr, ligne, roi_man, dossier):
    """Courbes VD/VG (signal moins ligne de base, toutes les images de la coupe)
    dans un .txt tabulé relisible par lire_export_segment. Les lignes « # »
    gardent ce qu'il faut pour refaire le calcul à l'identique (RR, fenêtre)."""
    dossier = Path(dossier)
    dossier.mkdir(parents=True, exist_ok=True)
    f = dossier / (nom_fichier_sur(id_patient) + ".txt")
    roi = ligne.get("_roi") or {}
    meta = [("version", "ptt_pipeline %s" % VERSION),
            ("RR_ms", "%.1f" % rr_ms if rr_ms else ""),
            ("source_RR", src_rr or ""),
            ("fin_premier_passage", int(b["fin"])),
            ("debut_analyse", int(b.get("debut", 0))),
            ("temps_reconstruit", int(bool(b.get("temps_reconstruit")))),
            ("serie", ligne.get("serie_num", "")),
            ("coupe", b["k"] + 1),
            ("ROI", "manuelles" if roi_man else "automatiques"),
            ("ROI_VD_ligne_colonne_rayon_px", ",".join(str(v) for v in roi.get("vd", []))),
            ("ROI_VG_ligne_colonne_rayon_px", ",".join(str(v) for v in roi.get("vg", [])))]
    lignes = ["# Courbes signal-temps (moins la ligne de base) enregistrées par ptt_pipeline.py.",
              "# A relire avec le mode « segment » (site : Courbes Segment) pour refaire les calculs sans les DICOM."]
    lignes += ["# %s\t%s" % (k, v) for k, v in meta]
    lignes += ["", "Temps (s)\tVD\tVG"]
    lignes += ["%.4f\t%.4f\t%.4f" % (t, a, c) for t, a, c in zip(b["temps"], b["y_vd"], b["y_vg"])]
    tmp = f.with_suffix(".tmp")
    tmp.write_text("\n".join(lignes) + "\n", encoding="utf-8")
    tmp.replace(f)
    return f.name


def indices_fenetre(temps, fenetre, fin_auto):
    """Fenêtre d'analyse choisie à la main : [début, fin] en secondes (None = auto).
    Renvoie (indice de début, indice de fin exclu) ; au moins 8 images gardées."""
    t_deb, t_fin = (list(fenetre) + [None, None])[:2]
    n = len(temps)
    deb = int(np.searchsorted(temps, float(t_deb) - 1e-6)) if t_deb is not None else 0
    fin = int(np.searchsorted(temps, float(t_fin) + 1e-6, side="right")) if t_fin is not None else fin_auto
    deb, fin = max(0, min(deb, n - 1)), max(1, min(fin, n))
    if fin - deb < 8:
        return 0, fin_auto
    return deb, fin


def dessiner_controle(id_patient, b, ligne, dossier_png):
    """Image PNG de contrôle : ROI sur l'image au pic VG + courbes."""
    t = b["temps"]
    deb = b.get("debut", 0)
    i_vg = deb + int(np.argmax(b["y_vg"][deb:b.get("fin", len(t))]))
    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    ax[0].imshow(b["pile"][i_vg], cmap="gray")
    ax[0].contour(b["cav"]["vd"], levels=[0.5], colors="deepskyblue", linewidths=1.5)
    ax[0].contour(b["cav"]["vg"], levels=[0.5], colors="red", linewidths=1.5)
    ax[0].set_title("%s - coupe %d - image au pic VG\nbleu = VD, rouge = VG%s" % (
        id_patient, b["k"] + 1, " (ROI manuelles)" if b.get("roi_manuelles") else ""))
    ax[0].axis("off")

    ax[1].plot(t, b["y_vd"], "o", ms=3, color="deepskyblue", label="VD (mesuré)")
    ax[1].plot(t, b["y_vg"], "o", ms=3, color="red", label="VG (mesuré)")
    fin = b.get("fin", len(t))
    if fin < len(t):   # images hors du premier passage : affichées, non utilisées
        ax[1].axvspan(t[fin - 1], t.max(), color="gray", alpha=0.15, lw=0)
        ax[1].text(0.99, 0.5, "hors premier passage\n(non utilisé)", transform=ax[1].transAxes,
                   fontsize=7, color="gray", ha="right")
    if deb > 0:        # début exclu à la main
        ax[1].axvspan(t.min(), t[deb], color="gray", alpha=0.15, lw=0)
    tt = np.linspace(t[deb], t[fin - 1], 500)
    for a, col in ((b["a_vd"], "deepskyblue"), (b["a_vg"], "red")):
        if a.get("params") is not None:
            ax[1].plot(tt, gamma(tt, *a["params"]), "-", color=col, lw=1.5)
        ax[1].axvline(a["t_centroide"], color=col, ls="--", lw=1)
        ax[1].axvline(a["t_pic_brut"], color=col, ls=":", lw=1)
    ax[1].set_xlabel("temps (s)")
    ax[1].set_ylabel("signal - ligne de base")
    txt = "PTT centroïde = %.2f s | PTT pic = %.2f s" % (ligne["PTT_centroide_s"], ligne["PTT_pic_s"])
    if ligne.get("PTTc_centroide") is not None:
        txt += "\nPTTc = %.2f | RR = %s ms" % (ligne["PTTc_centroide"], ligne["RR_ms"])
    if ligne.get("alertes"):
        txt += "\nALERTES : " + ligne["alertes"]
    ax[1].set_title(txt, fontsize=9)
    ax[1].legend(fontsize=8, loc="upper right")
    ax[1].text(0.01, 0.01, "tirets = centre de gravité ; pointillés = pic",
               transform=ax[1].transAxes, fontsize=7, color="gray")
    fig.tight_layout()
    fig.savefig(Path(dossier_png) / ("%s.png" % id_patient), dpi=110)
    plt.close(fig)


# =============================================================================
# 5. FICHIERS D'ENTRÉE (clinique, choix manuels)
# =============================================================================


def lire_table(chemin):
    if chemin.suffix.lower() in (".xlsx", ".xls"):
        return pd.read_excel(chemin, dtype={"ID": str})
    try:
        return pd.read_csv(chemin, sep=None, engine="python", dtype={"ID": str})
    except Exception:
        return pd.read_csv(chemin, sep=";", dtype={"ID": str})


def lire_clinique(base):
    for nom in ("cliniques.xlsx", "cliniques.csv"):
        p = base / nom
        if p.exists():
            df = lire_table(p)
            clin = {}
            for _, r in df.iterrows():
                d = {}
                def num(col):
                    try:
                        v = float(str(r[col]).replace(",", "."))
                        return v if not math.isnan(v) else None
                    except Exception:
                        return None
                d["VES_VG"] = num("VES_VG_mL") if "VES_VG_mL" in df else None
                d["VES_VD"] = num("VES_VD_mL") if "VES_VD_mL" in df else None
                d["FC"] = num("FC") if "FC" in df else None
                d["RR"] = num("RR_ms") if "RR_ms" in df else None
                h = num("Taille_cm") if "Taille_cm" in df else None
                w = num("Poids_kg") if "Poids_kg" in df else None
                d["SC"] = math.sqrt(h * w / 3600) if h and w else None  # Mosteller
                d["Groupe"] = r["Groupe"] if "Groupe" in df else None
                clin[str(r["ID"]).strip()] = d
            print("Fichier clinique lu : %s (%d patients)" % (nom, len(clin)))
            return clin
    print("Pas de fichier clinique (cliniques.xlsx) : le PBVi ne sera pas calculé.")
    return {}


def lire_choix(base):
    p = base / "choix_series.csv"
    choix = {}
    if p.exists():
        df = lire_table(p)
        for _, r in df.iterrows():
            d = {}
            try:
                d["serie"] = int(r["SeriesNumber"])
            except Exception:
                d["serie"] = None
            try:
                d["coupe"] = int(r["Coupe"]) if "Coupe" in df and str(r["Coupe"]) not in ("", "nan") else None
            except Exception:
                d["coupe"] = None
            choix[str(r["ID"]).strip()] = d
        print("Choix manuels lus : %d patients (choix_series.csv)" % len(choix))
    return choix


# =============================================================================
# 6. PROGRAMME PRINCIPAL
# =============================================================================


def patients(base, un_seul=None):
    d = base / "dicom"
    if not d.exists():
        print("\nERREUR : le dossier 'dicom' est introuvable dans", base)
        sys.exit(1)
    ids = sorted(p.name for p in d.iterdir() if p.is_dir())
    if un_seul:
        ids = [i for i in ids if i == un_seul]
        if not ids:
            print("ERREUR : patient %s introuvable dans dicom/" % un_seul)
            sys.exit(1)
    return [(i, d / i) for i in ids]


def mode_inventaire(base, un_seul):
    sortie = base / "resultats"
    sortie.mkdir(exist_ok=True)
    lignes, lignes_seq = [], []
    for i, (pid, dossier) in enumerate(patients(base, un_seul), 1):
        print("[%d] inventaire %s ..." % (i, pid))
        images = lire_entetes_patient(dossier)
        series = resumer_series(images)
        lignes_seq.append(dict({"ID": pid, "_sequences": inventaire_sequences(series)},
                               **informations_examen(images)))
        choisi, msg = choisir_serie(series)
        for s in sorted(series.values(), key=lambda s: s["serie_num"]):
            lignes.append({
                "ID": pid, "presence_sequence_perfusion": presence_perfusion(series),
                "SeriesNumber": s["serie_num"],
                "SeriesDescription": s["description"],
                "n_images": s["n_images"], "n_coupes": s["n_coupes"],
                "images_par_coupe": s["images_par_coupe"],
                "heure_acquisition_presente": "oui" if s["temps_dispo"] else "NON",
                "serie_choisie_automatiquement": "OUI" if (choisi and s is choisi) else "",
                "raison_du_choix": msg if (choisi and s is choisi) else "",
            })
    df = pd.DataFrame(lignes)
    f = sortie / "inventaire_series.csv"
    df.to_csv(f, index=False, sep=";", encoding="utf-8-sig")
    print("\nInventaire écrit dans :", f)
    if ecrire_sequences(lignes_seq, sortie):
        print("Séquences réalisées (0/1 par type) :", sortie / "sequences_realisees.csv")


PREMIERES_COLONNES = ["ID", "code_pseudo", "Groupe", "date_examen", "date_naissance", "age_ans", "sexe",
                      "taille_cm_dicom", "poids_kg_dicom", "SC_dicom_m2", "presence_sequence_perfusion", "statut", "validation_manuelle", "statut_auto", "commentaire", "PTT_principal_s", "PTT_centroide_s", "PTT_pic_s",
                      "PTTc_centroide", "PTTc_pic", "PTT_bazett_centroide", "PTT_bazett_pic", "FC_bpm", "PBVi_VD_centroide_mL_m2",
                      "PBVi_VG_centroide_mL_m2", "PBVi_VD_pic_mL_m2", "PBVi_VG_pic_mL_m2",
                      "ecart_centroide_pic_s", "LMH_VD_s", "LMH_VG_s", "rapport_LMH_VG_VD", "dispersion_VG",
                      "RR_ms", "source_RR", "alertes", "message"]


def lire_resultats_existants(f):
    """Relit un resultats_PTT.csv déjà produit (liste de dictionnaires)."""
    if not f.exists():
        return []
    try:
        df = pd.read_csv(f, sep=";", decimal=",", encoding="utf-8-sig", dtype={"ID": str})
    except Exception as e:
        print("ATTENTION : impossible de relire %s (%s) ; on repart de zéro." % (f.name, e))
        return []
    df = df.astype(object).where(pd.notna(df), None)
    return df.to_dict("records")


def _est_identifiant(mot):
    """Identifiant hospitalier : commence et finit par une lettre, avec au moins
    un chiffre entre les deux (ex. A1234567B)."""
    import re
    return bool(re.fullmatch(r"[A-Za-z][A-Za-z0-9]*[0-9][A-Za-z0-9]*[A-Za-z]", str(mot or "")))


def nom_prenom(identifiant):
    """Identifiant du patient (nom du dossier) -> (nom, prénom, identifiant) :
    le nom est le premier mot, le prénom tout ce qui suit (prénoms composés
    compris), sauf l'identifiant hospitalier (lettre + chiffres + lettre), mis à
    part. Espaces, « _ » et « ^ » (format DICOM) servent de séparateurs."""
    import re
    mots = [m for m in re.split(r"[\s_^]+", str(identifiant or "").strip()) if m]
    ident = next((m for m in mots if _est_identifiant(m)), None)
    mots = [m for m in mots if m != ident]
    if not mots:
        return None, None, ident
    return mots[0], (" ".join(mots[1:]) or None), ident


CLES_IDENTITE = ["_naissance", "_nom_dicom"]
COLONNES_IDENTITE = ["ID", "code_pseudo", "nom", "prenom", "identifiant", "date_naissance", "age_ans", "sexe",
                     "date_examen", "heure_examen", "source_nom"]


def _nom_dicom_utilisable(v):
    """Nom DICOM « NOM^PRENOM » (None si vide ou anonymisé)."""
    import re
    v = str(v or "").strip()
    if not v or re.search(r"anonym|^x+$|^patient\b|^\d+$", v, re.I):
        return None
    morceaux = [m.strip() for m in v.split("^")]
    if not morceaux[0] or not re.search(r"[A-Za-zÀ-ÿ]", morceaux[0]):
        return None
    return morceaux[0], (" ".join(m for m in morceaux[1:3] if m) or None)


def identites(lignes):
    """Tableau d'identité : un patient par ligne, nom / prénom (en-tête DICOM,
    sinon nom du dossier), identifiant, date de naissance et date d'examen.
    Données nominatives : ce fichier ne doit pas quitter l'hôpital."""
    out = []
    for l in lignes:
        nom, prenom, ident = nom_prenom(l.get("ID"))
        dicom = _nom_dicom_utilisable(l.get("_nom_dicom"))
        source = "dossier"
        if dicom:
            nom, prenom, source = dicom[0], dicom[1], "DICOM"
        idd = l.get("identifiant_dicom")
        out.append({"ID": l.get("ID"), "code_pseudo": l.get("code_pseudo"), "nom": nom, "prenom": prenom,
                    "identifiant": ident or (idd if isinstance(idd, str) and idd else None),
                    "date_naissance": l.get("_naissance") or l.get("date_naissance"),
                    "age_ans": l.get("age_ans"), "sexe": l.get("sexe"),
                    "date_examen": l.get("date_examen"), "heure_examen": l.get("heure_examen"),
                    "source_nom": source if nom else None})
    return _age_entier(pd.DataFrame(out, columns=COLONNES_IDENTITE))


def _age_entier(df):
    """Âge en années entières sans « ,0 » (cases vides permises)."""
    if len(df):
        a = pd.to_numeric(df["age_ans"], errors="coerce")
        if a.notna().any() and (a.dropna() % 1 == 0).all():
            df["age_ans"] = a.astype("Int64")
    return df


def ecrire_identites(lignes, dossier):
    """identites_patients.csv à côté des résultats (remplacé à chaque patient)."""
    f = Path(dossier) / "identites_patients.csv"
    df = identites(lignes).astype(object)
    if f.exists() and len(df):
        # reprise d'un calcul (lignes relues du CSV, sans le nom DICOM) : on
        # complète avec ce que le fichier précédent contenait déjà
        try:
            ancien = pd.read_csv(f, sep=";", encoding="utf-8-sig", dtype=str)
            ancien = ancien.astype(object).where(pd.notna(ancien), None)
            ancien = {str(r["ID"]): r for r in ancien.to_dict("records")}
            for i, pid in enumerate(df["ID"]):
                r = ancien.get(str(pid))
                if not r:
                    continue
                if r.get("source_nom") == "DICOM" and df.at[i, "source_nom"] != "DICOM":
                    for c in ("nom", "prenom", "source_nom"):
                        df.at[i, c] = r.get(c)
                for c in COLONNES_IDENTITE:
                    if c in r and r[c] is not None and pd.isna(df.at[i, c]):
                        df.at[i, c] = r[c]
        except Exception:
            pass
    tmp = f.with_name(f.stem + "_en_cours.tmp")
    df = _age_entier(df)
    df.to_csv(tmp, index=False, sep=";", encoding="utf-8-sig")
    try:
        tmp.replace(f)
    except PermissionError:      # ouvert dans Excel
        try:
            tmp.replace(f.with_name(f.stem + "_secours.csv"))
        except PermissionError:
            pass
    return f


def ajouter_nom_prenom(df):
    """Colonnes « nom », « prenom » et « identifiant » juste après la colonne ID
    (identifiant : pris dans le nom du dossier, sinon dans l'en-tête DICOM)."""
    if "ID" not in df.columns or not len(df):
        return df
    dicom = list(df["identifiant_dicom"]) if "identifiant_dicom" in df.columns else [None] * len(df)
    df = df.drop(columns=[c for c in ("nom", "prenom", "identifiant", "identifiant_dicom") if c in df.columns])
    np_ = [nom_prenom(x) for x in df["ID"]]
    i = list(df.columns).index("ID") + 1
    df.insert(i, "nom", [a for a, _, _ in np_])
    df.insert(i + 1, "prenom", [b for _, b, _ in np_])
    df.insert(i + 2, "identifiant", [c or (d if isinstance(d, str) and d else None)
                                     for (_, _, c), d in zip(np_, dicom)])
    return df


def ecrire_resultats(lignes, f):
    """Écrit le CSV complet (appelé après CHAQUE patient).
    Écriture dans un fichier temporaire puis remplacement : une coupure
    pendant l'écriture ne peut pas abîmer le fichier précédent."""
    df = pd.DataFrame(lignes)
    cols = [c for c in PREMIERES_COLONNES if c in df] + [c for c in df.columns
                                                         if c not in PREMIERES_COLONNES and not c.startswith("_")]
    df = ajouter_nom_prenom(df[cols])
    tmp = f.with_name(f.stem + "_en_cours.tmp")
    df.to_csv(tmp, index=False, sep=";", decimal=",", encoding="utf-8-sig")
    try:
        tmp.replace(f)
        return f
    except PermissionError:
        # Cas typique sous Windows : resultats_PTT.csv est ouvert dans Excel.
        secours = f.with_name(f.stem + "_secours.csv")
        try:
            tmp.replace(secours)
        except PermissionError:
            return tmp
        return secours


def mode_analyse(base, un_seul, reprendre=False):
    sortie = base / "resultats"
    png = sortie / "controles"
    png.mkdir(parents=True, exist_ok=True)
    f = sortie / "resultats_PTT.csv"
    clin = lire_clinique(base)
    choix = lire_choix(base)
    liste = patients(base, un_seul)

    # --- point de départ : résultats déjà présents ? -----------------------
    anciennes = lire_resultats_existants(f)
    if reprendre or un_seul:
        # on garde les résultats existants ; le(s) patient(s) analysé(s)
        # remplacent leur ancienne ligne
        lignes = anciennes
        if reprendre and not un_seul:
            faits = {str(l.get("ID")) for l in anciennes
                     if l.get("statut") in ("ok", "à vérifier", "sans perfusion")}
            avant = len(liste)
            liste = [(pid, d) for pid, d in liste if pid not in faits]
            print("Reprise : %d patient(s) déjà analysé(s) sautés, %d à faire "
                  "(les échecs éventuels sont retentés)." % (avant - len(liste), len(liste)))
    else:
        lignes = []
        if f.exists():
            horodatage = time.strftime("%Y%m%d_%H%M")
            sauvegarde = f.with_name("resultats_PTT_sauvegarde_%s.csv" % horodatage)
            try:
                f.replace(sauvegarde)
                print("L'ancien fichier de résultats est conservé sous :", sauvegarde.name)
            except PermissionError:
                print("\nERREUR : resultats_PTT.csv est ouvert (probablement dans Excel).")
                print("Fermez-le puis relancez la commande.\n")
                sys.exit(1)

    index = {str(l.get("ID")): k for k, l in enumerate(lignes)}
    debut = time.time()
    ecrit = f
    try:
        for i, (pid, dossier) in enumerate(liste, 1):
            print("[%d/%d] %s ..." % (i, len(liste), pid), end=" ", flush=True)
            try:
                ligne = analyser_patient(pid, dossier, clin, choix, png)
            except Exception as e:
                ligne = {"ID": pid, "statut": "échec", "message": "erreur : %s" % e}
                with open(sortie / "journal_erreurs.txt", "a", encoding="utf-8") as fj:
                    fj.write("=== %s (%s) ===\n%s\n" % (pid, time.strftime("%Y-%m-%d %H:%M"), traceback.format_exc()))
            if clin.get(pid, {}).get("Groupe") is not None:
                ligne["Groupe"] = clin[pid]["Groupe"]
            if pid in index:
                lignes[index[pid]] = ligne
            else:
                index[pid] = len(lignes)
                lignes.append(ligne)
            ecrit = ecrire_resultats(lignes, f)          # sauvegarde après chaque patient
            ecrire_sequences(lignes, sortie)
            ecrire_identites(lignes, sortie)
            reste = (time.time() - debut) / i * (len(liste) - i)
            print(ligne["statut"], "-", ligne.get("PTT_principal_s", ""),
                  ligne.get("alertes", "") or ligne.get("message", ""),
                  "| reste env. %d min" % round(reste / 60))
    except KeyboardInterrupt:
        print("\n\nInterrompu par l'utilisateur. Les patients déjà traités sont enregistrés.")
        print("Pour continuer plus tard :  python ptt_pipeline.py analyse --reprendre\n")
        return

    if not lignes:
        print("\nAucun patient à analyser.")
        return
    df = pd.DataFrame(lignes)
    n_ok = int((df["statut"] == "ok").sum())
    n_v = int((df["statut"] == "à vérifier").sum())
    n_e = int((df["statut"] == "échec").sum())
    n_s = int((df["statut"] == "sans perfusion").sum())
    print("\nTerminé en %d min : %d ok, %d à vérifier, %d échecs, %d sans perfusion (total du fichier)."
          % (round((time.time() - debut) / 60), n_ok, n_v, n_e, n_s))
    print("Résultats :", ecrit)
    if ecrit != f:
        print("ATTENTION : resultats_PTT.csv était ouvert (Excel ?). Les résultats sont dans",
              ecrit.name, "; fermez Excel et renommez ce fichier en resultats_PTT.csv.")
    print("Images de contrôle :", png)


# =============================================================================
# MESURES FAITES DANS SEGMENT (mode « segment »)
# =============================================================================
# Chaque export « ROI Analysis » de Segment est enregistré en .txt dans
# PTT_projet/exports_segment/, sous le code d'étude du patient (CMH-012.txt).
# Le programme lit tous les fichiers, applique la même méthode que l'analyse
# automatique (ligne de base, pic, fenêtre, gamma-variate) et écrit
# resultats/resultats_PTT_segment.csv (+ une image de contrôle par patient).
# La 2e ligne de l'export (nom du patient) n'est jamais lue ni recopiée.

MOTS_ROI_VD = ["right", "vd", "rv", "droit"]
MOTS_ROI_VG = ["left", "vg", "lv", "gauche"]

# Au-delà de cet intervalle entre images (ms), on soupçonne 1 image tous les
# 2 battements : le RR déduit des images est alors signalé.
INTERVALLE_SUSPECT_MS = 1200

# Courbes Segment : intervalle réel entre images (s) quand l'axe du temps de
# l'export est inutilisable (Segment exporte parfois un pas de quelques ms par
# image). None = repris de l'analyse automatique du même patient, sinon
# 1 image par battement (RR) supposée.
INTERVALLE_IMAGES_SEGMENT_S = None
# Fichiers de courbes enregistrés par ce programme : False = garder la fenêtre de
# premier passage trouvée sur les images (calcul identique) ; True = la
# recalculer sur les courbes (utile si vous modifiez la règle de la fenêtre).
RECALCULER_FENETRE_COURBES = False
NOMS_ONGLET_PAR_DEFAUT = ("sheet1", "feuil1", "feuille1", "blad1", "tabelle1", "hoja1", "foglio1")


def _nombre(texte):
    """'312.45' ou '312,45' -> 312.45 ; sinon None."""
    s = str(texte).strip().replace(",", ".")
    if not s:
        return None
    try:
        v = float(s)
        return v if math.isfinite(v) else None
    except ValueError:
        return None


def _contient_mot(texte, mots):
    """Recherche de mots entiers (évite que 'lv' soit trouvé dans 'value')."""
    import re
    t = re.sub(r"[^a-z]+", " ", str(texte).lower())
    return any((" %s " % m) in (" %s " % t) for m in mots)


EXTENSIONS_CLASSEUR = (".xlsx", ".xlsm")


def feuilles_classeur(chemin):
    """Onglets non vides d'un classeur Excel (un onglet = un patient)."""
    import openpyxl
    wb = openpyxl.load_workbook(chemin, read_only=True, data_only=True)
    try:
        return [ws.title for ws in wb.worksheets if (ws.max_row or 0) >= 9]
    finally:
        wb.close()


def _lignes_feuille(chemin, feuille):
    import openpyxl
    wb = openpyxl.load_workbook(chemin, read_only=True, data_only=True)
    try:
        if feuille not in wb.sheetnames:
            raise ValueError("onglet « %s » introuvable dans %s" % (feuille, Path(chemin).name))
        return [["" if v is None else (repr(v) if isinstance(v, float) else str(v)) for v in r]
                for r in wb[feuille].iter_rows(values_only=True)]
    finally:
        wb.close()


def _est_entete_temps(cel):
    c = str(cel).strip().lower()
    return c.startswith("time") or c.startswith("temps") or c in ("t", "t(s)", "t (s)", "t[s]", "t [s]", "t (ms)", "t [ms]")


def lire_export_segment(chemin, feuille=None):
    """Lit un export « ROI Analysis » de Segment (texte tabulé), ou un onglet
    d'un classeur Excel (feuille = nom de l'onglet) : une colonne de temps puis
    les colonnes VD et VG (nommées RV/VD/droit et LV/VG/gauche).
    Renvoie (t_s, y_vd_brut, y_vg_brut, infos) ; lève ValueError si illisible."""
    if Path(chemin).suffix.lower() in EXTENSIONS_CLASSEUR:
        lignes = _lignes_feuille(chemin, feuille)
    else:
        contenu = None
        for enc in ("utf-8-sig", "cp1252", "latin-1"):
            try:
                contenu = Path(chemin).read_text(encoding=enc)
                break
            except UnicodeDecodeError:
                continue
        if contenu is None:
            raise ValueError("encodage du fichier illisible")
        lignes = [l.rstrip("\r\n").split("\t") for l in contenu.split("\n")]

    # ligne d'en-tête : première cellule commençant par "time" (ou "temps")
    i_tete = None
    for i, cel in enumerate(lignes):
        if cel and _est_entete_temps(cel[0]):
            i_tete = i
            break
    if i_tete is None:
        raise ValueError("ligne d'en-tête 'Time' (ou 'Temps') introuvable en 1re colonne")
    tete = [str(c).strip() for c in lignes[i_tete]]
    etiq = [str(c).strip() for c in lignes[i_tete - 1]] if i_tete >= 1 else []
    etiq += [""] * (len(tete) - len(etiq))

    t0l = tete[0].lower().replace(" ", "")
    # unité du temps : indiquée dans l'en-tête, sinon devinée (pas médian > 50 = ms)
    unite = 1.0 if ("[s]" in t0l or "(s)" in t0l) else 0.001 if "ms" in t0l else None
    cols_mean = [j for j, h in enumerate(tete) if "mean" in h.lower()]

    def colonne_roi(mots):
        for j, lab in enumerate(etiq):
            if lab and _contient_mot(lab, mots):
                # colonne « Mean Intensity » du bloc de cette ROI
                for k in cols_mean:
                    if k >= j:
                        return k
                return j
        for j, h in enumerate(tete):            # étiquette dans l'en-tête lui-même
            if j and _contient_mot(h, mots):
                return j
        return None

    infos = {"alertes": [], "meta": {}}
    for cel in lignes[:i_tete]:          # fichier de courbes du programme : lignes « # clé<TAB>valeur »
        if cel and str(cel[0]).startswith("# ") and len(cel) >= 2:
            infos["meta"][str(cel[0])[2:].strip()] = str(cel[1]).strip()
    j_vd, j_vg = colonne_roi(MOTS_ROI_VD), colonne_roi(MOTS_ROI_VG)
    if j_vd is None or j_vg is None or j_vd == j_vg:
        if len(cols_mean) < 2:
            raise ValueError("ROI VD et VG non identifiées et moins de 2 colonnes 'Mean Intensity'")
        j_vd, j_vg = cols_mean[0], cols_mean[1]
        infos["alertes"].append("ROI non identifiées par leur nom (1re = VD, 2e = VG supposées)")
    infos["col_VD"], infos["col_VG"] = j_vd + 1, j_vg + 1

    t, vd, vg = [], [], []
    for cel in lignes[i_tete + 1:]:
        if not cel or not cel[0].strip():
            if t:
                break
            continue
        vals = [_nombre(cel[j]) if j < len(cel) else None for j in (0, j_vd, j_vg)]
        if None in vals:
            if t:
                break
            continue
        t.append(vals[0])
        vd.append(vals[1])
        vg.append(vals[2])
    if len(t) < 8:
        raise ValueError("moins de 8 images lisibles dans l'export")
    t = np.array(t)
    if unite is None:
        unite = 0.001 if np.median(np.diff(t)) > 50 else 1.0
    t = t * unite
    # un premier passage dure plusieurs dizaines de secondes : une série entière
    # de moins de 5 s veut dire que l'axe du temps exporté n'est pas le vrai temps
    infos["temps_inutilisable"] = bool(len(t) >= 15 and t[-1] - t[0] < 5.0)
    if np.any(np.diff(t) <= 0):
        raise ValueError("axe du temps non croissant")
    return t, np.array(vd), np.array(vg), infos


def fenetre_courbes(t, y_vd, y_vg):
    """Fin (exclue) du premier passage sur des courbes VD/VG : arrivée du
    contraste dans le VD + FENETRE_PREMIER_PASSAGE_S, et coupure avant une
    nouvelle montée du VG (2e injection, mouvement respiratoire...)."""
    n = len(t)
    if n < 8 or np.max(y_vd) <= 0:
        return n
    i_a = int(np.argmax(y_vd >= 0.3 * np.max(y_vd)))
    fin = int(np.searchsorted(t, t[i_a] + FENETRE_PREMIER_PASSAGE_S, side="right"))
    fin = min(n, max(fin, i_a + 8, N_IMAGES_BASE + 8))
    # premier pic du VG après l'arrivée, puis coupure si le VG remonte nettement
    yg = y_vg[:fin]
    if yg.max() > 0:
        i = i_a + int(np.argmax(yg[i_a:] >= 0.5 * yg[i_a:].max()))
        while i + 1 < fin and yg[i + 1] >= yg[i]:
            i += 1
        pic, i_min = yg[i], i
        for j in range(i + 1, fin):
            if yg[j] < yg[i_min]:
                i_min = j
            elif yg[j] - yg[i_min] > 0.3 * (pic - yg[i_min]) and yg[i_min] < 0.7 * pic:
                fin = i_min + 1
                break
    return max(fin, min(n, i_a + 8))


def analyser_export_segment(pid, chemin, clin, rr_dicom, dossier_png, feuille=None, intervalle_s=None,
                            fenetre=None):
    """intervalle_s : intervalle réel entre images (s), utilisé si l'axe du temps
    de l'export est inutilisable (sinon INTERVALLE_IMAGES_SEGMENT_S, sinon RR)."""
    ligne = {"ID": pid, "statut": "échec", "message": "", "source_mesure": "Segment"}
    t, vd, vg, infos = lire_export_segment(chemin, feuille)
    meta = infos.get("meta") or {}
    c = clin.get(pid, {})
    source_temps = "export Segment"
    if meta.get("version", "").startswith("ptt_pipeline"):
        ligne["source_mesure"] = "courbes sauvegardées (%s, série %s coupe %s, ROI %s)" % (
            meta["version"], meta.get("serie", "?"), meta.get("coupe", "?"), meta.get("ROI", "?"))
        source_temps = "fichier de courbes du programme"
    if infos.get("temps_inutilisable"):
        duree = t[-1] - t[0]
        pas = intervalle_s or INTERVALLE_IMAGES_SEGMENT_S
        if pas:
            source_temps = "reconstruit : %d images x %.3f s (le temps Segment couvrait %.2f s)" % (len(t), pas, duree)
        else:
            rr_tmp = c.get("RR") or rr_dicom.get(pid) or (60000.0 / c["FC"] if c.get("FC") else None)
            if not rr_tmp:
                ligne["message"] = ("temps de l'export Segment inutilisable (%.2f s pour %d images) : indiquez "
                                    "l'intervalle réel entre images, ou le RR/la FC du patient" % (duree, len(t)))
                return ligne
            pas = rr_tmp / 1000.0
            source_temps = "reconstruit : 1 image par battement supposée (RR %.0f ms ; le temps Segment couvrait %.2f s)" % (rr_tmp, duree)
            infos["alertes"].append("temps reconstruit (1 image/battement supposée) : vérifier l'intervalle entre images")
        t = np.arange(len(t)) * float(pas)
    y_vd = vd - vd[:N_IMAGES_BASE].mean()
    y_vg = vg - vg[:N_IMAGES_BASE].mean()
    try:            # fenêtre du premier passage déjà déterminée sur les images
        fin = min(len(t), int(meta["fin_premier_passage"])) if RECALCULER_FENETRE_COURBES is False else None
    except (KeyError, ValueError):
        fin = None
    if not fin or fin < 8:
        fin = fenetre_courbes(t, y_vd, y_vg)
    try:
        deb = int(meta.get("debut_analyse") or 0)
    except ValueError:
        deb = 0
    if fenetre:          # fenêtre choisie à la main (secondes)
        deb, fin = indices_fenetre(t, fenetre, fin)
        ligne["fenetre_manuelle"] = "oui"
    if fin - deb < 8:
        deb = 0
    a_vd, a_vg = analyser_courbe(t[deb:fin], y_vd[deb:fin]), analyser_courbe(t[deb:fin], y_vg[deb:fin])
    if not (a_vd["ok"] and a_vg["ok"]):
        ligne["message"] = "courbe VD ou VG non analysable"
        return ligne
    try:
        rr_fichier = float(meta.get("RR_ms") or 0) or None
    except ValueError:
        rr_fichier = None

    # RR : base clinique (RR_ms) > en-tête DICOM (analyse automatique) >
    #      intervalle entre images > FC clinique
    intervalle = float(np.median(np.diff(t)) * 1000)
    if c.get("RR"):
        rr_ms, src_rr = c["RR"], "RR_ms du fichier clinique"
    elif rr_fichier:
        rr_ms, src_rr = rr_fichier, "fichier de courbes (%s)" % (meta.get("source_RR") or "analyse d'origine")
    elif rr_dicom.get(pid):
        rr_ms, src_rr = rr_dicom[pid], "en-tête DICOM (repris de l'analyse automatique)"
    elif 300 <= intervalle <= 2000:
        rr_ms, src_rr = intervalle, "intervalle entre images (1 image/battement supposée)"
    elif c.get("FC"):
        rr_ms, src_rr = 60000.0 / c["FC"], "FC du fichier clinique"
    else:
        rr_ms, src_rr = None, "indisponible"

    ptt_c = a_vg["t_centroide"] - a_vd["t_centroide"]
    ptt_p = a_vg["t_pic_brut"] - a_vd["t_pic_brut"]
    alertes = list(infos["alertes"])
    if meta.get("temps_reconstruit") == "1":
        alertes.append("temps reconstruit à partir du RR (heures DICOM inutilisables)")
    ptt_princ = ptt_p if METHODE_PRINCIPALE == "pic" else ptt_c
    if not (PTT_MIN_S <= ptt_princ <= PTT_MAX_S):
        alertes.append("PTT hors bornes")
    if a_vg.get("ecrete"):
        alertes.append("plateau VG (saturation probable)")
    if a_vd.get("ecrete"):
        alertes.append("plateau VD (saturation probable)")
    disp, al = bilan_dispersion(a_vd, a_vg, ptt_c, ptt_p)
    alertes += al
    if a_vg["t_pic_brut"] <= a_vd["t_pic_brut"]:
        alertes.append("pic VG avant le pic VD (ROI inversées ?)")
    for nom, a in (("VD", a_vd), ("VG", a_vg)):
        if a.get("ajustement") == "brut":
            alertes.append("ajustement gamma échoué %s" % nom)
        elif a.get("r2", 1) < 0.90:
            alertes.append("ajustement %s médiocre (R2=%.2f)" % (nom, a["r2"]))
    if rr_ms is None:
        alertes.append("RR inconnu")
    elif src_rr.startswith("intervalle") and intervalle > INTERVALLE_SUSPECT_MS:
        alertes.append("intervalle entre images %.0f ms : 1 image tous les 2 battements ? "
                       "renseigner RR_ms dans cliniques.xlsx" % intervalle)

    ligne.update({
        "statut": "ok" if not alertes else "à vérifier",
        "fichier_export": Path(chemin).name + (" / onglet %s" % feuille if feuille else ""),
        "n_images": len(t),
        "intervalle_images_ms": round(intervalle, 1),
        "source_temps": source_temps,
        "fenetre_premier_passage_s": round(float(t[fin - 1]), 1),
        "fenetre_analyse_s": "%.1f-%.1f" % (t[deb], t[fin - 1]),
        "colonne_VD": infos["col_VD"], "colonne_VG": infos["col_VG"],
        "RR_ms": round(rr_ms, 1) if rr_ms else None,
        "source_RR": src_rr,
        "PTT_centroide_s": round(ptt_c, 2),
        "PTT_pic_s": round(ptt_p, 2),
        "R2_gamma_VD": round(a_vd.get("r2", float("nan")), 3),
        "R2_gamma_VG": round(a_vg.get("r2", float("nan")), 3),
        "plateau_VG_images": a_vg["n_plateau"],
        "alertes": " ; ".join(alertes),
    })
    ligne.update(disp)
    if rr_ms:
        ligne["FC_bpm"] = round(60000.0 / rr_ms, 1)
        for suf, ptt in (("centroide", ptt_c), ("pic", ptt_p)):
            ligne["PTTc_%s" % suf] = round(ptt / (rr_ms / 1000), 2)
            ligne["PTT_bazett_%s" % suf] = round(ptt / math.sqrt(rr_ms / 1000), 2)
    sc = c.get("SC")
    if sc:
        ligne["SC_m2"] = round(sc, 3)
        fc_ref = c.get("FC") or (60000.0 / rr_ms if rr_ms else None)
        for suf, ptt in (("centroide", ptt_c), ("pic", ptt_p)):
            if c.get("VES_VD") and rr_ms:
                ligne["PBVi_VD_%s_mL_m2" % suf] = round((ptt / (rr_ms / 1000)) * c["VES_VD"] / sc, 1)
            if c.get("VES_VG") and fc_ref:
                ligne["PBVi_VG_%s_mL_m2" % suf] = round(ptt * (c["VES_VG"] * fc_ref / 60.0) / sc, 1)
    ligne["PTT_principal_s"] = ligne["PTT_centroide_s"] if METHODE_PRINCIPALE == "centroide" else ligne["PTT_pic_s"]

    try:
        dessiner_controle_segment(pid, t, y_vd, y_vg, a_vd, a_vg, ligne, dossier_png, fin, deb)
    except Exception as e:
        ligne["message"] += "image de contrôle non créée (%s)" % e
    return ligne


def dessiner_controle_segment(pid, t, y_vd, y_vg, a_vd, a_vg, ligne, dossier_png, fin=None, deb=0):
    fig, ax = plt.subplots(figsize=(8, 5))
    if deb:
        ax.axvspan(t[0], t[deb], color="grey", alpha=0.15)
    if fin is not None and fin < len(t):
        ax.axvspan(t[fin - 1], t[-1], color="grey", alpha=0.15, label="hors premier passage (non utilisé)")
    ax.plot(t, y_vd, "o", ms=3, color="deepskyblue", label="VD (Segment)")
    ax.plot(t, y_vg, "o", ms=3, color="red", label="VG (Segment)")
    tt = np.linspace(t.min(), t.max(), 500)
    for a, col in ((a_vd, "deepskyblue"), (a_vg, "red")):
        if a.get("params") is not None:
            ax.plot(tt, gamma(tt, *a["params"]), "-", color=col, lw=1.5)
        d, f = a["fenetre"]
        ax.axvspan(t[d], t[f], color=col, alpha=0.06)
        ax.axvline(a["t_centroide"], color=col, ls="--", lw=1)
        ax.axvline(a["t_pic_brut"], color=col, ls=":", lw=1)
    txt = "%s (Segment) | PTT pic = %.2f s | PTT centroïde = %.2f s" % (
        pid, ligne["PTT_pic_s"], ligne["PTT_centroide_s"])
    if ligne.get("PTTc_pic") is not None:
        txt += "\nPTTc pic = %.2f | RR = %s ms (%s)" % (ligne["PTTc_pic"], ligne["RR_ms"], ligne["source_RR"])
    if ligne.get("alertes"):
        txt += "\nALERTES : " + ligne["alertes"]
    ax.set_title(txt, fontsize=8)
    ax.set_xlabel("temps (s)")
    ax.set_ylabel("signal - ligne de base")
    ax.legend(fontsize=8, loc="upper right")
    ax.text(0.01, 0.01, "zones colorées = fenêtres de premier passage ; tirets = centroïde ; pointillés = pic",
            transform=ax.transAxes, fontsize=7, color="gray")
    fig.tight_layout()
    fig.savefig(Path(dossier_png) / ("%s.png" % pid), dpi=110)
    plt.close(fig)


def comparer_auto_segment(sortie, lignes_seg):
    """Compare, patient par patient, l'analyse automatique et la mesure Segment."""
    f_auto = sortie / "resultats_PTT.csv"
    auto = {str(l.get("ID")): l for l in lire_resultats_existants(f_auto)}
    comp = []
    for s in lignes_seg:
        a = auto.get(str(s.get("ID")))
        if not a or s.get("statut") == "échec" or a.get("statut") == "échec":
            continue
        d = {"ID": s["ID"], "statut_auto": a.get("statut"), "statut_segment": s.get("statut")}
        for m in ("PTT_pic_s", "PTT_centroide_s"):
            try:
                va, vs = float(a.get(m)), float(s.get(m))
            except (TypeError, ValueError):
                continue
            d[m + "_auto"], d[m + "_segment"] = va, vs
            d[m + "_diff_auto_moins_segment"] = round(va - vs, 2)
        comp.append(d)
    if not comp:
        return
    df = pd.DataFrame(comp)
    f = sortie / "comparaison_auto_segment.csv"
    df.to_csv(f, index=False, sep=";", decimal=",", encoding="utf-8-sig")
    print("\nComparaison automatique / Segment (%d patients) : %s" % (len(df), f.name))
    for m in ("PTT_pic_s", "PTT_centroide_s"):
        col = m + "_diff_auto_moins_segment"
        if col in df and df[col].notna().sum() >= 3:
            x = df[col].dropna().astype(float)
            biais, et = x.mean(), x.std(ddof=1)
            print("  %-16s biais = %+.2f s ; limites d'agrément = %.2f à %.2f s (n = %d)"
                  % (m, biais, biais - 1.96 * et, biais + 1.96 * et, len(x)))


def mode_segment(base, un_seul, reprendre=False):
    dossier = base / "exports_segment"
    if not dossier.exists():
        print("\nERREUR : le dossier 'exports_segment' est introuvable dans", base)
        print("Créez-le et enregistrez-y les exports de Segment (un fichier .txt par patient,")
        print("nommé avec le code d'étude, par exemple CMH-012.txt), ou un classeur Excel")
        print("avec un onglet par patient (nom de l'onglet = code d'étude).\n")
        sys.exit(1)
    # (identifiant, fichier, onglet) : un .txt = un patient ; un classeur Excel =
    # un patient par onglet (nom de l'onglet = identifiant)
    fichiers = []
    for p in sorted(dossier.iterdir()):
        if p.suffix.lower() in (".txt", ".csv"):
            fichiers.append((p.stem.strip(), p, None))
        elif p.suffix.lower() in EXTENSIONS_CLASSEUR and not p.name.startswith("~$"):
            noms = feuilles_classeur(p)
            for nom in noms:
                ident = p.stem.strip() if (len(noms) == 1 and str(nom).strip().lower() in NOMS_ONGLET_PAR_DEFAUT) else str(nom).strip()
                fichiers.append((ident, p, nom))
    if un_seul:
        fichiers = [x for x in fichiers if x[0] == un_seul]
    if not fichiers:
        print("\nAucun export trouvé dans exports_segment/%s.\n" % ((" pour " + un_seul) if un_seul else ""))
        sys.exit(1)
    sortie = base / "resultats"
    png = sortie / "controles_segment"
    png.mkdir(parents=True, exist_ok=True)
    f = sortie / "resultats_PTT_segment.csv"
    clin = lire_clinique(base)
    rr_dicom, pas_dicom = {}, {}
    for l in lire_resultats_existants(sortie / "resultats_PTT.csv"):
        try:
            if l.get("intervalle_images_s"):
                pas_dicom[str(l["ID"])] = float(str(l["intervalle_images_s"]).replace(",", "."))
        except (TypeError, ValueError):
            pass
        if str(l.get("source_RR") or "").startswith("en-tête DICOM") and l.get("RR_ms"):
            try:
                rr_dicom[str(l["ID"])] = float(l["RR_ms"])
            except (TypeError, ValueError):
                pass
    if rr_dicom:
        print("RR de l'en-tête DICOM disponible (analyse automatique) pour %d patient(s)." % len(rr_dicom))

    lignes = lire_resultats_existants(f) if (reprendre or un_seul) else []
    if reprendre and not un_seul:
        faits = {str(l.get("ID")) for l in lignes if l.get("statut") in ("ok", "à vérifier")}
        fichiers = [x for x in fichiers if x[0] not in faits]
        print("Reprise : %d export(s) à traiter." % len(fichiers))
    elif not (reprendre or un_seul) and f.exists():
        sauvegarde = f.with_name("resultats_PTT_segment_sauvegarde_%s.csv" % time.strftime("%Y%m%d_%H%M"))
        try:
            f.replace(sauvegarde)
            print("L'ancien fichier est conservé sous :", sauvegarde.name)
        except PermissionError:
            print("\nERREUR : resultats_PTT_segment.csv est ouvert (Excel ?). Fermez-le et relancez.\n")
            sys.exit(1)
    index = {str(l.get("ID")): k for k, l in enumerate(lignes)}
    ecrit = f
    for i, (pid, p, feuille) in enumerate(fichiers, 1):
        print("[%d/%d] %s ..." % (i, len(fichiers), pid), end=" ", flush=True)
        try:
            ligne = analyser_export_segment(pid, p, clin, rr_dicom, png, feuille, pas_dicom.get(pid))
        except ValueError as e:
            ligne = {"ID": pid, "statut": "échec", "source_mesure": "Segment", "message": str(e)}
        except Exception as e:
            ligne = {"ID": pid, "statut": "échec", "source_mesure": "Segment", "message": "erreur : %s" % e}
            with open(sortie / "journal_erreurs.txt", "a", encoding="utf-8") as fj:
                fj.write("=== segment %s (%s) ===\n%s\n" % (pid, time.strftime("%Y-%m-%d %H:%M"),
                                                        traceback.format_exc()))
        if clin.get(pid, {}).get("Groupe") is not None:
            ligne["Groupe"] = clin[pid]["Groupe"]
        if pid in index:
            lignes[index[pid]] = ligne
        else:
            index[pid] = len(lignes)
            lignes.append(ligne)
        ecrit = ecrire_resultats(lignes, f)
        print(ligne["statut"], "-", ligne.get("PTT_principal_s", ""),
              ligne.get("alertes", "") or ligne.get("message", ""))
    df = pd.DataFrame(lignes)
    print("\nTerminé : %d ok, %d à vérifier, %d échecs." % (
        int((df["statut"] == "ok").sum()), int((df["statut"] == "à vérifier").sum()),
        int((df["statut"] == "échec").sum())))
    print("Résultats :", ecrit)
    print("Images de contrôle :", png)
    comparer_auto_segment(sortie, lignes)



# =============================================================================
# REPORT DES RÉSULTATS DANS LA BASE DE RECUEIL EXCEL (mode « importer »)
# =============================================================================
# Colonne du CSV -> colonne (ligne 2) de la feuille « Base »
CORRESPONDANCE_BASE = {
    "serie_num": "Perf_serie_num",
    "coupe_retenue": "Perf_coupe",
    "RR_ms": "RR_perf_ms",
    "source_RR": "Source_RR",
    "PTT_pic_s": "PTT_pic_s",
    "PTT_centroide_s": "PTT_centroide_s",
    "alertes": "Alertes_script",
}


def mode_importer(base, chemin_base, source="auto"):
    """Recopie les résultats de resultats_PTT.csv dans une COPIE de la base
    Excel, en faisant correspondre les patients par leur code (ID_etude).
    Le fichier d'origine n'est jamais modifié."""
    try:
        from openpyxl import load_workbook
    except ImportError:
        print("\nERREUR : la bibliothèque openpyxl manque. Tapez :")
        print("  python -m pip install openpyxl\n")
        sys.exit(1)
    f_auto = base / "resultats" / "resultats_PTT.csv"
    f_seg = base / "resultats" / "resultats_PTT_segment.csv"
    besoin = {"auto": [f_auto], "segment": [f_seg], "fusion": [f_auto, f_seg]}[source]
    for f_csv in besoin:
        if not f_csv.exists():
            print("\nERREUR : %s introuvable. Lancez d'abord l'analyse correspondante.\n" % f_csv)
            sys.exit(1)
    if chemin_base:
        f_base = Path(chemin_base)
        if not f_base.is_absolute():
            f_base = base / f_base
    else:
        candidats = sorted(p for p in base.glob("Base_IRM*.xlsx") if "_PTT_" not in p.name)
        if not candidats:
            print("\nERREUR : aucune base Excel (Base_IRM*.xlsx) dans", base)
            print("Copiez la base dans le dossier du projet, ou indiquez-la avec --base.\n")
            sys.exit(1)
        f_base = candidats[0]
    if not f_base.exists():
        print("\nERREUR : base Excel introuvable :", f_base, "\n")
        sys.exit(1)

    def lire_ok(f_csv, etiquette):
        d = pd.read_csv(f_csv, sep=";", decimal=",", encoding="utf-8-sig", dtype={"ID": str})
        d = d[d["statut"].isin(["ok", "à vérifier"])].copy()
        d["source_mesure"] = etiquette
        return d
    if source == "auto":
        res = lire_ok(f_auto, "automatique")
    elif source == "segment":
        res = lire_ok(f_seg, "Segment")
    else:
        ra, rs_ = lire_ok(f_auto, "automatique"), lire_ok(f_seg, "Segment")
        res = pd.concat([ra[~ra["ID"].isin(rs_["ID"])], rs_], ignore_index=True)
        print("Fusion : %d patient(s) mesuré(s) dans Segment remplacent l'analyse automatique."
              % int(ra["ID"].isin(rs_["ID"]).sum()))
    print("Base lue :", f_base.name, "| résultats exploitables :", len(res), "(source : %s)" % source)

    wb = load_workbook(f_base)
    if "Base" not in wb.sheetnames:
        print("\nERREUR : la feuille « Base » est absente de", f_base.name, "\n")
        sys.exit(1)
    ws = wb["Base"]
    col = {}
    for c in range(1, ws.max_column + 1):
        v = ws.cell(row=2, column=c).value
        if v:
            col[str(v).strip()] = c
    manquantes = [b for b in list(CORRESPONDANCE_BASE.values()) + ["ID_etude"] if b not in col]
    if manquantes:
        print("\nERREUR : colonnes absentes de la base :", ", ".join(manquantes), "\n")
        sys.exit(1)
    lignes_id = {}
    for r in range(4, ws.max_row + 1):
        v = ws.cell(row=r, column=col["ID_etude"]).value
        if v not in (None, ""):
            lignes_id.setdefault(str(v).strip(), []).append(r)

    n_ok, absents, doublons = 0, [], []
    for _, l in res.iterrows():
        pid = str(l["ID"]).strip()
        rows = lignes_id.get(pid, [])
        if not rows:
            absents.append(pid)
            continue
        if len(rows) > 1:
            doublons.append(pid)
            continue
        r = rows[0]
        for c_csv, c_base in CORRESPONDANCE_BASE.items():
            v = l.get(c_csv)
            if v is None or (isinstance(v, float) and math.isnan(v)):
                v = None
            elif isinstance(v, (np.integer,)):
                v = int(v)
            elif isinstance(v, (np.floating,)):
                v = float(v)
            ws.cell(row=r, column=col[c_base]).value = v
        if "Source_PTT" in col:
            ws.cell(row=r, column=col["Source_PTT"]).value = l.get("source_mesure")
        elif l.get("source_mesure") == "Segment":
            al0 = ws.cell(row=r, column=col["Alertes_script"]).value
            ws.cell(row=r, column=col["Alertes_script"]).value = ("[mesure Segment] " + (al0 or "")).strip()
        if "Plateau_VG" in col:
            al = str(l.get("alertes") or "")
            ws.cell(row=r, column=col["Plateau_VG"]).value = "Oui" if "plateau VG" in al else "Non"
        n_ok += 1

    try:
        wb.calculation.fullCalcOnLoad = True    # Excel recalcule tout à l'ouverture
    except Exception:
        pass
    sortie = f_base.with_name(f_base.stem + "_PTT_" + time.strftime("%Y%m%d_%H%M") + ".xlsx")
    wb.save(sortie)
    print("\n%d patient(s) reporté(s) dans : %s" % (n_ok, sortie.name))
    if absents:
        print("ATTENTION : %d code(s) du CSV absent(s) de la colonne ID_etude de la base :" % len(absents))
        print("  " + ", ".join(absents))
    if doublons:
        print("ATTENTION : code(s) présent(s) plusieurs fois dans la base (non reportés) :")
        print("  " + ", ".join(doublons))
    print("Le fichier d'origine %s n'a pas été modifié.\n" % f_base.name)


def mode_diagnostic(base, un_seul):
    """Pour chaque série de perfusion (y compris ORIG / Low Res) : les heures
    DICOM disponibles par image. Sert à comprendre un échec. Aucun pixel lu."""
    sortie = base / "resultats"
    sortie.mkdir(exist_ok=True)
    lignes = []
    for pid, dossier in patients(base, un_seul):
        print("diagnostic %s ..." % pid)
        series = resumer_series(lire_entetes_patient(dossier))
        for s_ in sorted(series.values(), key=lambda x: x["serie_num"]):
            if not contient(s_["description"], MOTS_PERFUSION):
                continue
            for k, pos in enumerate(sorted(s_["coupes"].keys()), 1):
                imgs = sorted(s_["coupes"][pos], key=lambda im: im["instance"])
                d = {"ID": pid, "SeriesNumber": s_["serie_num"], "SeriesDescription": s_["description"],
                     "coupe": k, "n_images": len(imgs)}
                for nom, cle, f in (("AcquisitionTime", "temps_acq", 1), ("ContentTime", "temps_contenu", 1),
                                    ("TriggerTime", "trigger_ms", 0.001)):
                    vals = [im.get(cle) for im in imgs]
                    ok = [float(v) * f for v in vals if v is not None]
                    d[nom + "_presentes"] = len(ok)
                    d[nom + "_distinctes"] = len(set(round(v, 3) for v in ok))
                    d[nom + "_etendue_s"] = round(max(ok) - min(ok), 2) if ok else None
                    d[nom + "_utilisable"] = "oui" if _temps_valable(vals if f == 1 else
                                                                      [float(v) * f if v is not None else None for v in vals]) else "non"
                d["HeartRate"] = imgs[0].get("fc_tag")
                d["NominalInterval_ms"] = imgs[0].get("rr_ms_tag")
                lignes.append(d)
    f = sortie / "diagnostic_temps.csv"
    pd.DataFrame(lignes).to_csv(f, index=False, sep=";", decimal=",", encoding="utf-8-sig")
    print("\nDiagnostic écrit dans :", f)



def main():
    ap = argparse.ArgumentParser(description="Calcul du PTT/PBVi à partir des perfusions DICOM")
    ap.add_argument("mode", nargs="?", default="analyse", choices=["inventaire", "analyse", "importer", "segment", "diagnostic"])
    ap.add_argument("--patient", default=None, help="analyser un seul patient (nom du dossier)")
    ap.add_argument("--dossier", default=None, help="dossier du projet (par défaut : celui du script)")
    ap.add_argument("--base", default=None,
                    help="mode importer : base Excel de recueil (par défaut : Base_IRM*.xlsx du dossier du projet)")
    ap.add_argument("--source", default="auto", choices=["auto", "segment", "fusion"],
                    help="mode importer : auto (défaut), segment, ou fusion (Segment prioritaire)")
    ap.add_argument("--reprendre", action="store_true",
                    help="saute les patients déjà présents (ok / à vérifier) dans resultats_PTT.csv")
    args = ap.parse_args()
    base = Path(args.dossier).resolve() if args.dossier else Path(__file__).resolve().parent
    print("Dossier du projet :", base)
    if args.mode == "inventaire":
        mode_inventaire(base, args.patient)
    elif args.mode == "importer":
        mode_importer(base, args.base, args.source)
    elif args.mode == "diagnostic":
        mode_diagnostic(base, args.patient)
    elif args.mode == "segment":
        mode_segment(base, args.patient, args.reprendre)
    else:
        mode_analyse(base, args.patient, args.reprendre)


if __name__ == "__main__":
    main()
