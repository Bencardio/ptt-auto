# -*- coding: utf-8 -*-
"""
web_ptt.py
==========

Version en ligne du site PTT : fait tourner site_ptt.py (et donc
ptt_pipeline.py, sans rien changer aux calculs) dans le navigateur, avec
Pyodide. Ce fichier remplace seulement ce qui dépend d'un vrai serveur :

- les requêtes /api/... de la page arrivent ici au lieu de passer par HTTP ;
- les fils d'exécution (analyse de la cohorte, relances) sont remplacés par
  pas() : un patient à la fois, appelé en boucle par web/worker.js, pour que
  la page reste réactive entre deux patients ;
- les DICOM sont lus sur place dans le dossier choisi (montage /dicom/...),
  sans copie ni envoi ;
- les résultats, pseudonymes.csv, prefixe_pseudo.txt et cliniques.xlsx sont
  dans le « dossier de travail » choisi sur le PC (monté sur /travail).

Rien ne quitte l'ordinateur : tout est calculé dans l'onglet du navigateur.
"""

import io
import os
import sys
import time
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")

import ptt_pipeline as pp

pp.N_FILS_LECTURE = 1          # pas de fils d'exécution dans le navigateur

import site_ptt as s

PRINCIPAUX = {}    # lot.id -> contexte d'analyse (analyse de la cohorte en cours)
RELANCES = {}      # lot.id -> contexte d'analyse (relances manuelles en attente)


class _FilDiffere:
    """Remplace threading.Thread dans site_ptt : au lieu de démarrer un fil,
    inscrit le travail ; pas() le fait ensuite patient par patient."""

    def __init__(self, target=None, args=(), daemon=None, **_):
        self.target, self.args = target, args

    def start(self):
        lot = self.args[0]
        if self.target is s.lancer_analyse:
            lot.debut = time.time()
            lot.fait = 0
            ctx = s.contexte_analyse(lot)
            ctx["png"].mkdir(parents=True, exist_ok=True)
            PRINCIPAUX[lot.id] = ctx
        elif self.target is s.lancer_relances:
            ctx = s.contexte_analyse(lot)
            ctx["png"].mkdir(parents=True, exist_ok=True)
            RELANCES[lot.id] = ctx
        else:
            self.target(*self.args)


class _Threading:
    """Module threading de site_ptt : Thread différé, le reste inchangé."""

    def __init__(self, vrai):
        self._vrai = vrai

    def __getattr__(self, nom):
        return _FilDiffere if nom == "Thread" else getattr(self._vrai, nom)


s.threading = _Threading(s.threading)


def demarrer(dossier_travail):
    """Le dossier de travail remplace le dossier du site (BASE)."""
    base = Path(dossier_travail)
    base.mkdir(parents=True, exist_ok=True)
    s.BASE = base
    s.TRAVAIL = base / "site_travail"
    s.RESULTATS = base / "resultats_site"
    s.FICHIER_PSEUDO = base / "pseudonymes.csv"
    s.FICHIER_PREFIXE = base / "prefixe_pseudo.txt"
    s._pseudos = None
    s.LOTS.clear()
    PRINCIPAUX.clear()
    RELANCES.clear()
    return {"version": pp.VERSION, "dossier": str(base)}


def travail_en_attente():
    return bool(PRINCIPAUX or RELANCES)


def pas():
    """Un patient (relance manuelle d'abord, comme sur le site local, où elle
    passe devant la cohorte). Renvoie True s'il reste du travail."""
    for lid, ctx in list(RELANCES.items()):
        lot = s.LOTS.get(lid)
        if lot is None:
            RELANCES.pop(lid)
            continue
        with lot.verrou:
            if not lot.file_relance:
                lot.relance = None
                if not lot.principal_actif and not lot.file:
                    lot.etat = "termine"
                RELANCES.pop(lid)
                continue
            pid, c = lot.file_relance.pop(0)
            lot.relance = pid
            version = lot.version.get(pid, 0)
        s.traiter_patient(lot, pid, c, ctx, version, principal=False)
        return True
    for lid, ctx in list(PRINCIPAUX.items()):
        lot = s.LOTS.get(lid)
        if lot is None:
            PRINCIPAUX.pop(lid)
            continue
        with lot.verrou:
            if not lot.file:
                lot.en_cours = None
                lot.etape = ""
                lot.principal_actif = False
                if not lot.relance and not lot.file_relance:
                    lot.etat = "termine"
                PRINCIPAUX.pop(lid)
                continue
            pid, c = lot.file.pop(0)
            lot.en_cours = pid
            version = lot.version.get(pid, 0)
        s.traiter_patient(lot, pid, c, ctx, version, principal=True)
        lot.fait += 1
        with lot.verrou:
            # la page n'est mise à jour qu'entre deux patients : on annonce le suivant
            if lot.file:
                lot.en_cours = lot.file[0][0]
                lot.etape = "lecture des en-têtes DICOM"
        return True
    return False


# -----------------------------------------------------------------------------
# Requêtes de la page (au lieu du serveur HTTP)
# -----------------------------------------------------------------------------
class _Entetes(dict):
    def get(self, cle, defaut=None):
        for k, v in self.items():
            if k.lower() == cle.lower():
                return v
        return defaut


class _Gestionnaire(s.Gestionnaire):
    def __init__(self, chemin, entetes, corps):
        self.path = chemin
        self.headers = _Entetes(entetes or {})
        self.headers["Content-Length"] = str(len(corps or b""))
        self.rfile = io.BytesIO(corps or b"")
        self.wfile = io.BytesIO()
        self.code = 500
        self.entetes_reponse = {}

    def send_response(self, code, message=None):
        self.code = code

    def send_header(self, cle, valeur):
        self.entetes_reponse[cle] = valeur

    def end_headers(self):
        pass


def _sortie(code, entetes, corps):
    """Réponse pour web/worker.js : le corps passe en Uint8Array JavaScript."""
    try:
        from pyodide.ffi import to_js
        corps = to_js(memoryview(bytes(corps)))
    except ImportError:     # hors navigateur (tests)
        pass
    return code, entetes, corps


def _json(code, obj):
    import json
    return _sortie(code, {"Content-Type": "application/json; charset=utf-8"},
                   json.dumps(obj, ensure_ascii=False).encode("utf-8"))


def requete(methode, chemin, entetes, corps):
    """(code, en-têtes, octets) : même réponse que le serveur local."""
    if hasattr(corps, "to_bytes"):          # Uint8Array JavaScript
        corps = corps.to_bytes()
    if isinstance(corps, str):
        corps = corps.encode("utf-8")
    elif corps is not None and not isinstance(corps, (bytes, bytearray)):
        corps = bytes(corps)
    route = chemin.split("?", 1)[0]
    if route == "/api/dossier_local":
        import json
        data = json.loads((corps or b"{}").decode("utf-8") or "{}")
        if not str(data.get("chemin", "")).startswith("/dicom/"):
            return _json(400, {"erreur": "dans la version en ligne, le navigateur ne peut pas ouvrir un dossier à partir de "
                                         "son chemin : utilisez le bouton « Choisir un dossier sur ce PC et analyser »"})
    g = _Gestionnaire(chemin, entetes, corps)
    if methode == "GET":
        g.do_GET()
    else:
        g.do_POST()
    return _sortie(g.code, g.entetes_reponse, g.wfile.getvalue())


def ajouter_depot(lid, base):
    """Fichiers glissés-déposés (montés sur place sous base/<patient>/...) :
    même effet que l'envoi au serveur local (recevoir_paquet), sans copie."""
    lot = s.LOTS.get(lid)
    if lot is None:
        return _json(404, {"erreur": "lot inconnu (rechargez la page)"})
    if lot.etat == "en_cours":
        return _json(400, {"erreur": "analyse en cours"})
    base = Path(base)
    pids = []
    for d in sorted(x for x in base.iterdir() if x.is_dir()):
        pid = s.nom_sur(d.name)
        n = sum(1 for f in d.rglob("*") if f.is_file())
        lot.patients[pid] = {"dossier": d, "fichiers": n, "copie": False}
        lot.lignes.pop(pid, None)
        lot.inventaire.pop(pid, None)
        s.oublier_entetes(d)
        pids.append(pid)
    if lot.mode == "segment":
        # classeur Excel déposé : un patient par onglet (comme recevoir_paquet)
        for pid in pids:
            d = lot.patients.get(pid)
            if not d or d.get("feuille"):
                continue
            classeurs = [x for x in Path(d["dossier"]).rglob("*")
                         if x.is_file() and x.suffix.lower() in pp.EXTENSIONS_CLASSEUR]
            if classeurs:
                del lot.patients[pid]
                try:
                    for c in classeurs:
                        s.ajouter_export_segment(lot, c, copie=False)
                except Exception as e:
                    return _json(400, {"erreur": "classeur Excel illisible : %s" % e})
    return _json(200, {"ok": True, "n": len(pids)})
