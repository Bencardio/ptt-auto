# -*- coding: utf-8 -*-
"""
mise_a_jour.py
==============

Met l'outil à jour depuis GitHub au lancement (appelé par Lancer_site_PTT.bat
et Lancer_site_PTT_Mac.command, ou à la main : python mise_a_jour.py).

- Télécharge la dernière version du dépôt DEPOT (branche BRANCHE), sans git.
- Ne remplace que les fichiers du programme (FICHIERS_PROGRAMME) ; les dossiers
  de travail (resultats_site, site_travail, pseudonymes.csv, cliniques.xlsx...)
  ne sont jamais touchés.
- Les fichiers remplacés sont gardés dans anciennes_versions/<date_heure>/.
- Sans internet (ou GitHub bloqué par le réseau), le site démarre avec la
  version déjà installée.
"""

import io
import json
import shutil
import sys
import time
import urllib.request
import zipfile
from pathlib import Path

DEPOT = "Bencardio/ptt-auto"
BRANCHE = "main"
FICHIERS_PROGRAMME = ["ptt_pipeline.py", "site_ptt.py", "bolustrak_ptt.py", "page_site_ptt.html",
                      "mise_a_jour.py", "LISEZMOI.md", "Methode_calcul_PTT.docx", "exemple_classeur_segment.xlsx"]
# Les lanceurs (.bat, .command) ne sont pas remplacés pendant qu'ils tournent.

BASE = Path(__file__).resolve().parent
FICHIER_VERSION = BASE / ".version_installee"
DELAI_S = 10


def lire(url):
    req = urllib.request.Request(url, headers={"User-Agent": "ptt-auto-mise-a-jour"})
    with urllib.request.urlopen(req, timeout=DELAI_S) as r:
        return r.read()


def main():
    try:
        info = json.loads(lire("https://api.github.com/repos/%s/commits/%s" % (DEPOT, BRANCHE)))
        sha = info["sha"]
    except Exception as e:
        print("Mise à jour : GitHub injoignable (%s), version installée conservée." % e.__class__.__name__)
        return
    installee = FICHIER_VERSION.read_text().strip() if FICHIER_VERSION.exists() else ""
    if sha == installee:
        print("Mise à jour : déjà la dernière version (%s)." % sha[:7])
        return
    try:
        archive = zipfile.ZipFile(io.BytesIO(lire("https://codeload.github.com/%s/zip/%s" % (DEPOT, sha))))
    except Exception as e:
        print("Mise à jour : téléchargement impossible (%s), version installée conservée." % e.__class__.__name__)
        return
    sauvegarde = BASE / "anciennes_versions" / time.strftime("%Y%m%d_%H%M%S")
    changes = []
    for nom in archive.namelist():
        fichier = nom.split("/", 1)[-1]
        if fichier not in FICHIERS_PROGRAMME:
            continue
        nouveau = archive.read(nom)
        cible = BASE / fichier
        if cible.exists() and cible.read_bytes() == nouveau:
            continue
        if cible.exists():
            sauvegarde.mkdir(parents=True, exist_ok=True)
            shutil.copy2(cible, sauvegarde / fichier)
        cible.write_bytes(nouveau)
        changes.append(fichier)
    FICHIER_VERSION.write_text(sha)
    msg = info.get("commit", {}).get("message", "").splitlines()[0] if info.get("commit") else ""
    if changes:
        print("Mise à jour installée (%s) : %s" % (sha[:7], ", ".join(changes)))
        if msg:
            print("  Dernière modification : %s" % msg)
    else:
        print("Mise à jour : fichiers déjà à jour (%s)." % sha[:7])


if __name__ == "__main__":
    try:
        main()
    except Exception as e:  # la mise à jour ne doit jamais empêcher le site de démarrer
        print("Mise à jour ignorée : %s" % e)
    sys.exit(0)
