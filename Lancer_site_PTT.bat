@echo off
REM Double-cliquez sur ce fichier pour ouvrir le site PTT dans le navigateur.
REM Au lancement, l'outil se met a jour depuis GitHub (mise_a_jour.py).
cd /d "%~dp0"
where py >nul 2>nul
if %errorlevel%==0 (
  py mise_a_jour.py
  py site_ptt.py
) else (
  python mise_a_jour.py
  python site_ptt.py
)
pause
