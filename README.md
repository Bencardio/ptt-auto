# Outil PTT / PBVi (perfusion de premier passage, IRM cardiaque)

Outil de recherche (non certifié, non destiné au diagnostic) : temps de transit
pulmonaire (PTT pic / centroïde, PTTc, Bazett), largeurs à mi-hauteur, PBVi, sur
les perfusions de premier passage, les séries 2D BolusTrak ou des courbes
exportées de Segment.

**Ouvrir l'outil : https://bencardio.github.io/ptt-auto/**

- Le programme s'exécute entièrement dans le navigateur (Python via Pyodide) :
  les DICOM sont lus sur l'ordinateur de l'utilisateur et ne sont envoyés à aucun serveur.
- Version locale (Python installé sur le PC) : `Lancer_site_PTT.bat` / `Lancer_site_PTT_Mac.command`,
  avec mise à jour automatique au lancement.
- Mode d'emploi complet : `LISEZMOI.md` ; méthode de calcul : `Methode_calcul_PTT.docx`.
- Ce dépôt ne contient que le programme, jamais de données patient (voir `.gitignore`).
