#!/bin/bash
# Double-cliquez sur ce fichier (Mac) pour ouvrir le site PTT dans le navigateur.
# Au lancement, l'outil se met à jour depuis GitHub (mise_a_jour.py).
cd "$(dirname "$0")"
python3 mise_a_jour.py
python3 site_ptt.py
