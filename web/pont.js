/* Pont de la version en ligne du site PTT (GitHub Pages).
   La page page_site_ptt.html est la même que sur le site local : ses requêtes
   /api/... sont envoyées ici à Python (Pyodide, web/worker.js) au lieu d'un
   serveur. Les DICOM sont lus sur place dans le dossier choisi, les résultats
   écrits dans le « dossier de travail » choisi : rien ne quitte le PC. */
"use strict";
(() => {
  const PYODIDE = new URLSearchParams(location.search).get("pyodide") || "https://cdn.jsdelivr.net/pyodide/v0.27.7/full/";
  const BASE = new URL(".", location.href).href;
  const ROUES = ["pydicom-3.0.1-py3-none-any.whl", "openpyxl-3.1.5-py2.py3-none-any.whl", "et_xmlfile-2.0.0-py3-none-any.whl"];
  const IGNORES = /^(\.|thumbs\.db$|desktop\.ini$)/i;
  const vraiFetch = window.fetch.bind(window);

  const travailleur = new Worker(BASE + "web/worker.js");
  let compteur = 0, pret = false, travailPret = false;
  const attentes = {}, attendre = [];
  let depot = [];          // fichiers glissés-déposés, montés au lancement de l'analyse
  let nDepot = 0;

  travailleur.onmessage = e => {
    const m = e.data;
    if (m.type === "etat") return bandeau(m.msg);
    if (m.type === "echec") return bandeau("Impossible de charger l'outil : " + m.msg + ". Rechargez la page.", true);
    if (m.type === "pret") { pret = true; versionPipeline = m.version; return demanderTravail(); }
    if (m.type === "reponse" && attentes[m.id]) { attentes[m.id](m); delete attentes[m.id]; }
  };
  let versionPipeline = "";
  function envoyer(msg, transfert) {
    return new Promise(ok => { msg.id = ++compteur; attentes[msg.id] = ok; travailleur.postMessage(msg, transfert || []); });
  }
  travailleur.postMessage({type: "demarrer", pyodide: PYODIDE, base: BASE, roues: ROUES});

  /* ---------- bandeau du haut : chargement et dossier de travail ---------- */
  const css = document.createElement("style");
  css.textContent = `#pttWeb{position:sticky;top:0;z-index:50;background:#1f2a44;color:#fff;padding:10px 16px;font:14px system-ui,sans-serif;
    display:flex;gap:10px;flex-wrap:wrap;align-items:center}#pttWeb button{font:inherit;padding:5px 10px;border-radius:6px;border:0;cursor:pointer}
    #pttWeb .p{background:#4c8dff;color:#fff}#pttWeb .s{background:#e8ecf4;color:#1d2330}#pttWeb.erreur{background:#8b1d1d}
    #pttWeb small{opacity:.85}.ptt-web-cache{display:none!important}`;
  document.head.appendChild(css);
  const B = document.createElement("div");
  B.id = "pttWeb";
  const placer = () => { if (!B.isConnected) document.body.prepend(B); };
  if (document.body) placer(); else document.addEventListener("DOMContentLoaded", placer);
  function bandeau(html, erreur) { placer(); B.classList.toggle("erreur", !!erreur); B.innerHTML = html; }
  bandeau("Version en ligne : chargement de Python dans le navigateur (environ 1 minute la première fois, ensuite plus rapide)…");

  /* ---------- poignée du dossier de travail gardée d'une visite à l'autre ---------- */
  function idb(action, valeur) {
    return new Promise(ok => {
      try {
        const r = indexedDB.open("ptt-web", 1);
        r.onupgradeneeded = () => r.result.createObjectStore("poignees");
        r.onerror = () => ok(null);
        r.onsuccess = () => {
          try {
            const t = r.result.transaction("poignees", action === "lire" ? "readonly" : "readwrite").objectStore("poignees");
            const q = action === "lire" ? t.get("travail") : t.put(valeur, "travail");
            q.onsuccess = () => ok(q.result || null); q.onerror = () => ok(null);
          } catch (e) { ok(null); }
        };
      } catch (e) { ok(null); }
    });
  }

  async function taille(poignee, limite) {
    // un dossier de résultats reste petit : on refuse un dossier qui contient les DICOM
    let n = 0, octets = 0;
    const tour = async d => {
      for await (const h of d.values()) {
        if (n > limite.n || octets > limite.o) return;
        if (h.kind === "directory") await tour(h);
        else { n++; if (n % 200 === 0 || n < 50) octets += (await h.getFile()).size; }
      }
    };
    await tour(poignee);
    return n > limite.n || octets > limite.o;
  }

  async function utiliserTravail(poignee) {
    if (travailPret) {          // changement de dossier : la page repart de zéro sur le nouveau dossier
      if (poignee) await idb("ecrire", poignee);
      sessionStorage.setItem("ptt-sans-dossier", poignee ? "" : "1");
      return location.reload();
    }
    if (poignee) {
      bandeau(`Ouverture du dossier de travail « ${esc(poignee.name)} »…`);
      if (await taille(poignee, {n: 20000, o: 1.5e9}))
        return demanderTravail(`Le dossier « ${esc(poignee.name)} » est trop gros pour un dossier de résultats (il contient sans doute les DICOM). Créez un dossier vide, par exemple « PTT_resultats », et choisissez-le.`);
      await idb("ecrire", poignee);
    }
    const r = await envoyer({type: "travail", poignee: poignee || null});
    travailPret = true;
    bandeau(`<span>Version en ligne · pipeline ${esc(versionPipeline)} · les DICOM restent sur ce PC.</span>
      <small>${poignee ? `Dossier de travail : <b>${esc(poignee.name)}</b> (résultats, pseudonymes.csv, cliniques.xlsx)`
                       : "<b>Sans dossier de travail</b> : téléchargez les résultats avant de fermer l'onglet."}</small>
      <button class="s" id="pttChanger">Changer</button>`);
    document.getElementById("pttChanger").onclick = () => demanderTravail();
    while (attendre.length) attendre.shift()();
  }

  async function demanderTravail(message) {
    const ancienne = await idb("lire");
    const peut = "showDirectoryPicker" in window;
    if (!message && !travailPret) {
      // permission déjà accordée (Chrome : « Autoriser à chaque visite ») : rien à cliquer
      if (ancienne && sessionStorage.getItem("ptt-sans-dossier") !== "1") {
        try { if ((await ancienne.queryPermission({mode: "readwrite"})) === "granted") return utiliserTravail(ancienne); } catch (e) {}
      }
      if (sessionStorage.getItem("ptt-sans-dossier") === "1") return utiliserTravail(null);
    }
    bandeau(`<span>${message || "Choisissez le <b>dossier de travail</b> : les résultats, pseudonymes.csv et cliniques.xlsx y sont enregistrés (un dossier vide la première fois, par ex. « PTT_resultats », pas celui des DICOM)."}</span>
      ${ancienne && peut ? `<button class="p" id="pttReprendre">Reprendre « ${esc(ancienne.name)} »</button>` : ""}
      ${peut ? `<button class="${ancienne ? "s" : "p"}" id="pttChoisir">Choisir un dossier…</button>` : `<small>Ce navigateur ne peut pas écrire dans un dossier (utilisez Chrome ou Edge).</small>`}
      <button class="s" id="pttSans">Continuer sans dossier</button>`, !!message);
    const b = id => document.getElementById(id);
    if (b("pttReprendre")) b("pttReprendre").onclick = async () => {
      try {
        if ((await ancienne.requestPermission({mode: "readwrite"})) !== "granted") return;
        await utiliserTravail(ancienne);
      } catch (e) { demanderTravail("Dossier introuvable : choisissez-le à nouveau."); }
    };
    if (b("pttChoisir")) b("pttChoisir").onclick = async () => {
      let p; try { p = await showDirectoryPicker({id: "ptt-travail", mode: "readwrite"}); } catch (e) { return; }
      utiliserTravail(p);
    };
    b("pttSans").onclick = () => utiliserTravail(null);
  }

  /* ---------- choix d'un dossier DICOM : lu sur place ---------- */
  async function lister(poignee) {
    const sortie = [];
    const tour = async (d, chemin) => {
      for await (const h of d.values()) {
        if (IGNORES.test(h.name)) continue;
        const c = chemin ? chemin + "/" + h.name : h.name;
        if (h.kind === "directory") await tour(h, c);
        else {
          sortie.push({nom: c, blob: await h.getFile()});
          if (sortie.length % 500 === 0) bandeauTemp(`Inventaire du dossier : ${sortie.length} fichiers…`);
        }
      }
    };
    await tour(poignee, "");
    return sortie;
  }
  function listerParInput() {
    return new Promise(ok => {
      const i = document.createElement("input");
      i.type = "file"; i.webkitdirectory = true; i.multiple = true;
      i.onchange = () => {
        const f = [...i.files];
        if (!f.length) return ok(null);
        const racine = (f[0].webkitRelativePath || "").split("/")[0] || "dossier";
        ok({nom: racine, fichiers: f.filter(x => !x.webkitRelativePath.split("/").some(p => IGNORES.test(p)))
                                  .map(x => ({nom: x.webkitRelativePath.split("/").slice(1).join("/"), blob: x}))});
      };
      i.oncancel = () => ok(null);
      i.click();
    });
  }
  let texteBandeau = "";
  function bandeauTemp(t) { if (!texteBandeau) texteBandeau = B.innerHTML; B.innerHTML = esc(t); }
  function finTemp() { if (texteBandeau) { B.innerHTML = texteBandeau; texteBandeau = ""; const c = document.getElementById("pttChanger"); if (c) c.onclick = () => demanderTravail(); } }
  const sur = n => n.replace(/[^\w\-. ]/g, "_").slice(0, 60) || "dossier";

  async function choisirDossierDicom() {
    let nom, fichiers;
    if ("showDirectoryPicker" in window) {
      let p; try { p = await showDirectoryPicker({id: "ptt-dicom", mode: "read"}); } catch (e) { return {chemin: ""}; }
      nom = p.name; fichiers = await lister(p);
    } else {
      const r = await listerParInput();
      if (!r) return {chemin: ""};
      nom = r.nom; fichiers = r.fichiers;
    }
    finTemp();
    if (!fichiers.length) return {chemin: "", erreur: "dossier vide"};
    const base = "/dicom/" + sur(nom);
    await envoyer({type: "monter", base, fichiers});
    return {chemin: base};
  }

  /* ---------- requêtes /api/... : à Python au lieu d'un serveur ---------- */
  function reponse(code, corps, entetes) {
    return new Response(corps, {status: code, headers: entetes || {"Content-Type": "application/json; charset=utf-8"}});
  }
  const json = (code, obj) => reponse(code, JSON.stringify(obj));

  async function api(url, opt) {
    if (!pret || !travailPret) await new Promise(ok => attendre.push(ok));
    const u = new URL(url, location.href), route = u.pathname.replace(/^.*\/api\//, "/api/");
    const chemin = route + u.search;
    const methode = (opt && opt.method) || "GET";
    const entetes = Object.assign({}, (opt && opt.headers) || {});
    let corps = opt && opt.body;
    if (route === "/api/choisir_dossier") return json(200, await choisirDossierDicom());
    if (route === "/api/fichiers") {
      // glisser-déposer : pas de copie, les fichiers seront lus sur place
      const liste = JSON.parse(decodeURIComponent(entetes["X-Fichiers"] || "[]"));
      let o = 0;
      for (const [pid, rel, t] of liste) { depot.push({nom: sur(pid) + "/" + rel, blob: corps.slice(o, o + t)}); o += t; }
      return json(200, {ok: true, n: liste.length});
    }
    if (route === "/api/analyser" && depot.length) {
      const lot = u.searchParams.get("lot") || JSON.parse(corps).lot;
      const r = await envoyer({type: "depot", base: "/dicom/depot_" + (++nDepot), fichiers: depot, lot});
      depot = [];
      if (r.code !== 200) return reponse(r.code, r.corps, r.entetes);
    }
    if (corps instanceof Blob) corps = await corps.arrayBuffer();
    const r = await envoyer({type: "api", methode, chemin, entetes, corps: corps ?? null},
                            corps instanceof ArrayBuffer ? [corps] : []);
    return reponse(r.code, r.corps, r.entetes);
  }

  window.fetch = (url, opt) => {
    const s = typeof url === "string" ? url : url.url;
    if (/^\/api\//.test(s) || s.startsWith(BASE + "api/")) return api(s, opt);
    return vraiFetch(url, opt);
  };

  function esc(t) { return String(t ?? "").replace(/[&<>"]/g, c => ({"&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;"}[c])); }

  window.PTT_WEB = {
    _test: {utiliserTravail},
    async telecharger(url) {
      const r = await fetch(url);
      if (!r.ok) { const j = await r.json().catch(() => ({})); return alert(j.erreur || "téléchargement impossible"); }
      const cd = r.headers.get("Content-Disposition") || "", m = /filename="([^"]+)"/.exec(cd);
      const a = document.createElement("a");
      a.href = URL.createObjectURL(await r.blob());
      a.download = m ? m[1] : "resultats";
      document.body.appendChild(a); a.click(); a.remove();
      setTimeout(() => URL.revokeObjectURL(a.href), 60000);
    },
  };

  // images de contrôle : <img data-png="/api/png?..."> fabriquées par Python
  new MutationObserver(() => {
    document.querySelectorAll("img[data-png]").forEach(async img => {
      const u = img.dataset.png; img.removeAttribute("data-png");
      const r = await fetch(u);
      if (!r.ok) return img.remove();
      img.src = URL.createObjectURL(await r.blob());
    });
  }).observe(document.documentElement, {childList: true, subtree: true});

  // options sans objet en ligne : chemin à coller (le navigateur ne lit pas un dossier par son chemin), copies DICOM
  document.addEventListener("DOMContentLoaded", () => {
    const c = document.getElementById("chemin");
    if (c) c.closest("div[style]").classList.add("ptt-web-cache");
    const e = document.getElementById("btnEffacer");
    if (e) e.classList.add("ptt-web-cache");
    const n = document.getElementById("btnNatif");
    if (n && n.nextElementSibling) n.nextElementSibling.textContent = "Les fichiers sont lus sur place par le navigateur : rien n'est copié ni envoyé.";
  });
})();
