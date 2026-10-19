/* Version en ligne du site PTT : Python (Pyodide) dans un Web Worker.
   Les calculs sont ceux de ptt_pipeline.py / site_ptt.py, chargés tels quels
   depuis le dépôt. Les DICOM sont lus sur place (WORKERFS), les résultats
   écrits dans le dossier de travail choisi (mountNativeFS). Rien n'est envoyé. */
"use strict";
let py = null, natif = null, travail = false, boucle = false;
const file = [];     // messages reçus pendant le démarrage

function dire(msg) { postMessage({type: "etat", msg}); }

async function demarrer({pyodide, base, roues}) {
  importScripts(pyodide + "pyodide.js");
  dire("Chargement de Python dans le navigateur…");
  py = await loadPyodide({indexURL: pyodide});
  dire("Chargement des bibliothèques (numpy, scipy, pandas, matplotlib)…");
  await py.loadPackage(["numpy", "scipy", "pandas", "matplotlib", "micropip"]);
  const micropip = py.pyimport("micropip");
  await micropip.install(roues.map(r => base + "web/wheels/" + r));
  dire("Chargement de l'outil PTT…");
  py.FS.mkdirTree("/app/web");
  for (const f of ["ptt_pipeline.py", "site_ptt.py", "bolustrak_ptt.py", "web/web_ptt.py"]) {
    const r = await fetch(base + f, {cache: "no-cache"});
    if (!r.ok) throw new Error(`${f} introuvable (${r.status})`);
    py.FS.writeFile("/app/" + f, new Uint8Array(await r.arrayBuffer()));
  }
  py.FS.mkdirTree("/dicom");
  await py.runPythonAsync(`
import sys, os
os.chdir("/app")
sys.path[:0] = ["/app", "/app/web"]
sys.argv = ["site_ptt.py", "--sans-navigateur"]
import web_ptt
`);
  postMessage({type: "pret", version: py.runPython("web_ptt.pp.VERSION")});
}

async function choisirTravail(poignee) {
  if (natif) { try { await natif.syncfs(); } catch (e) {} }
  if (poignee) {
    try { py.FS.unmount("/travail"); } catch (e) {}
    try { py.FS.rmdir("/travail"); } catch (e) {}
    natif = await py.mountNativeFS("/travail", poignee);
  } else {
    natif = null;
    py.FS.mkdirTree("/travail");
  }
  travail = true;
  return py.runPython("import json; json.dumps(web_ptt.demarrer('/travail'))");
}

async function sauver() { if (natif) { try { await natif.syncfs(); } catch (e) { dire("Écriture des résultats impossible : " + e.message); } } }

function monter(base, fichiers) {
  // base = /dicom/<nom> ; fichiers = [{nom: "patient/serie/IM_0001", blob}]
  try { py.FS.unmount(base); } catch (e) {}
  try { py.FS.rmdir(base); } catch (e) {}
  py.FS.mkdirTree(base);
  py.FS.mount(py.FS.filesystems.WORKERFS, {blobs: fichiers.map(f => ({name: f.nom, data: f.blob}))}, base);
}

function appelPython(nom, ...args) {
  const f = py.globals.get("web_ptt")[nom];
  const r = f(...args);
  const code = r.get(0), e = r.get(1), corps = r.get(2);
  const entetes = e.toJs({dict_converter: Object.fromEntries});
  e.destroy(); r.destroy(); f.destroy();
  return {code, entetes, corps};
}

async function lancerBoucle() {
  if (boucle) return;
  boucle = true;
  const tour = async () => {
    let encore = false;
    try { encore = py.runPython("web_ptt.pas()"); }
    catch (e) { dire("Erreur pendant l'analyse : " + e.message); }
    await sauver();
    if (encore) setTimeout(tour, 0); else boucle = false;
  };
  setTimeout(tour, 0);
}

async function traiter(m) {
  if (m.type === "travail") return {retour: JSON.parse(await choisirTravail(m.poignee))};
  if (m.type === "monter") { monter(m.base, m.fichiers); return {retour: true}; }
  if (m.type === "depot") {
    monter(m.base, m.fichiers);
    return appelPython("ajouter_depot", m.lot, m.base);
  }
  if (m.type === "api") {
    const corps = m.corps == null ? null : (typeof m.corps === "string" ? m.corps : new Uint8Array(m.corps));
    const r = appelPython("requete", m.methode, m.chemin, py.toPy(m.entetes || {}), corps);
    if (m.methode !== "GET") await sauver();
    if (py.runPython("web_ptt.travail_en_attente()")) lancerBoucle();
    return r;
  }
  throw new Error("message inconnu : " + m.type);
}

onmessage = async e => {
  const m = e.data;
  if (m.type === "demarrer") {
    try { await demarrer(m); } catch (er) { postMessage({type: "echec", msg: er.message}); return; }
    while (file.length) await repondre(file.shift());
    return;
  }
  if (!py) { file.push(m); return; }
  await repondre(m);
};

async function repondre(m) {
  try {
    const r = await traiter(m);
    if (r && r.corps instanceof Uint8Array) {
      postMessage({type: "reponse", id: m.id, ...r, corps: r.corps.buffer}, [r.corps.buffer]);
    } else postMessage({type: "reponse", id: m.id, ...r});
  } catch (er) {
    postMessage({type: "reponse", id: m.id, code: 500, entetes: {"Content-Type": "application/json"},
                 corps: new TextEncoder().encode(JSON.stringify({erreur: "erreur interne : " + er.message})).buffer});
  }
}
