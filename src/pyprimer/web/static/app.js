/* PyPrimer web interface. */
(() => {
  "use strict";

  const $ = (sel, root = document) => root.querySelector(sel);
  const $$ = (sel, root = document) => Array.from(root.querySelectorAll(sel));
  const PREFS_KEY = "pyprimer.prefs.v1";
  const MAP_WIDTH = 80;

  const state = {
    defaults: null,
    prefs: null,
    results: {}, // mode -> { data, rows, request, selected, sort, enzymes, cloning }
  };

  // ------------------------------------------------------------------ utils

  const escapeHtml = (s) =>
    String(s ?? "").replace(/[&<>"']/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" })[c]);

  const fmt = (v, digits = 2) => (typeof v === "number" ? v.toFixed(digits) : v ?? "");

  async function api(path, body, { method = body ? "POST" : "GET", raw = false } = {}) {
    const response = await fetch(`/api${path}`, {
      method,
      headers: body ? { "Content-Type": "application/json" } : {},
      body: body ? JSON.stringify(body) : undefined,
    });
    if (!response.ok) {
      let message = `Request failed (${response.status})`;
      try {
        message = (await response.json()).error || message;
      } catch { /* not JSON */ }
      throw new Error(message);
    }
    return raw ? response : response.json();
  }

  function storageGet(key) {
    try {
      return JSON.parse(localStorage.getItem(key));
    } catch {
      return null;
    }
  }

  function storageSet(key, value) {
    try {
      localStorage.setItem(key, JSON.stringify(value));
    } catch { /* storage unavailable */ }
  }

  const clone = (obj) => JSON.parse(JSON.stringify(obj));

  function getPath(obj, path) {
    return path.split(".").reduce((o, k) => (o == null ? undefined : o[k]), obj);
  }

  function setPath(obj, path, value) {
    const keys = path.split(".");
    let o = obj;
    keys.slice(0, -1).forEach((k) => {
      if (o[k] == null || typeof o[k] !== "object") o[k] = {};
      o = o[k];
    });
    o[keys[keys.length - 1]] = value;
  }

  function readInput(el) {
    if (el.type === "checkbox") return el.checked;
    if (el.dataset.type === "intlist") {
      const parts = el.value.split(/[\s,;]+/).filter(Boolean);
      if (!parts.length) return null;
      const nums = parts.map(Number);
      if (nums.some((n) => !Number.isInteger(n))) throw new Error("Exon boundaries must be whole numbers");
      return nums;
    }
    if (el.type === "number") return el.value === "" ? null : Number(el.value);
    return el.value;
  }

  function writeInput(el, value) {
    if (value === undefined) return;
    if (el.type === "checkbox") el.checked = Boolean(value);
    else if (Array.isArray(value)) el.value = value.join(", ");
    else el.value = value ?? "";
    if (el.tagName === "TEXTAREA") updateLength(el);
  }

  function readForm(form) {
    const out = {};
    $$("[data-path]", form).forEach((el) => setPath(out, el.dataset.path, readInput(el)));
    return out;
  }

  function fillForm(form, obj) {
    $$("[data-path]", form).forEach((el) => writeInput(el, getPath(obj, el.dataset.path)));
  }

  function field(panel, path) {
    return $(`[data-path="${path}"]`, panel);
  }

  function setStatus(root, message, kind = "") {
    const el = $('[data-role="status"]', root);
    if (!el) return;
    el.innerHTML = message ? `<p class="msg ${kind}">${escapeHtml(message)}</p>` : "";
  }

  function cleanSequence(text) {
    let body = text.trim();
    if (body.startsWith(">")) {
      const lines = body.split(/\r?\n/);
      const end = lines.findIndex((l, i) => i > 0 && l.startsWith(">"));
      body = lines.slice(1, end === -1 ? undefined : end).join("");
    }
    return body.replace(/[\s\d-]/g, "");
  }

  function updateLength(textarea) {
    const label = $(`[data-len-for="${textarea.dataset.path}"]`, textarea.closest(".panel"));
    if (label) label.textContent = `${cleanSequence(textarea.value).length.toLocaleString()} bp`;
  }

  function reverseComplement(seq) {
    const map = { A: "T", C: "G", G: "C", T: "A", a: "t", c: "g", g: "c", t: "a" };
    return seq.split("").reverse().map((c) => map[c] || c).join("");
  }

  function download(filename, text, type = "text/plain") {
    const url = URL.createObjectURL(new Blob([text], { type }));
    const a = Object.assign(document.createElement("a"), { href: url, download: filename });
    document.body.append(a);
    a.click();
    a.remove();
    setTimeout(() => URL.revokeObjectURL(url), 1000);
  }

  async function busy(button, fn) {
    if (button) button.disabled = true;
    try {
      return await fn();
    } finally {
      if (button) button.disabled = false;
    }
  }

  // ------------------------------------------------------------ preferences

  function defaultPrefs() {
    const d = state.defaults;
    return {
      conditions: clone(d.conditions),
      gc: clone(d.pcr.gc),
      repeats: clone(d.pcr.repeats),
      bis_repeats: clone(d.bisulphite.repeats),
      cpg: clone(d.cpg),
    };
  }

  function loadPrefs() {
    const saved = storageGet(PREFS_KEY);
    const prefs = defaultPrefs();
    if (saved && typeof saved === "object") {
      for (const key of Object.keys(prefs)) Object.assign(prefs[key], saved[key] || {});
    }
    return prefs;
  }

  function initPrefs() {
    const panel = $("#tab-prefs");
    const form = $("form", panel);
    fillForm(form, state.prefs);
    form.addEventListener("change", () => {
      state.prefs = readForm(form);
      storageSet(PREFS_KEY, state.prefs);
      setStatus(panel, "Saved in this browser");
    });
    $('[data-action="reset-prefs"]', panel).addEventListener("click", () => {
      state.prefs = defaultPrefs();
      fillForm(form, state.prefs);
      storageSet(PREFS_KEY, state.prefs);
      setStatus(panel, "Defaults restored");
    });
  }

  // ------------------------------------------------------------------- tabs

  function initTabs() {
    const tabs = $$(".tabs [role=tab]");
    const show = (name) => {
      tabs.forEach((t) => t.setAttribute("aria-selected", String(t.dataset.tab === name)));
      $$(".panel").forEach((p) => (p.hidden = p.id !== `tab-${name}`));
      storageSet("pyprimer.tab", name);
    };
    tabs.forEach((t) => t.addEventListener("click", () => show(t.dataset.tab)));
    const saved = storageGet("pyprimer.tab");
    if (saved && $(`#tab-${saved}`)) show(saved);
  }

  // ---------------------------------------------------------- design panels

  const PAIR_COLUMNS = [
    ["Forward", "forward", "seq"],
    ["Pos", "forward_pos", "num"],
    ["Len", "forward_len", "num"],
    ["Tm", "forward_tm", "tm"],
    ["Reverse", "reverse", "seq"],
    ["Pos", "reverse_end", "num"],
    ["Len", "reverse_len", "num"],
    ["Tm", "reverse_tm", "tm"],
    ["Amp", "amplicon", "num"],
    ["Ext. dimer ΔG", "dimer_dg", "tm"],
    ["Full dimer ΔG", "full_dimer_dg", "tm"],
  ];
  const SEQ_COLUMNS = [
    ["Primer", "seq", "seq"],
    ["Pos", "pos", "num"],
    ["Len", "length", "num"],
    ["Tm", "tm", "tm"],
    ["Self-dimer ΔG", "dimer_dg", "tm"],
  ];

  function buildRequest(mode, form) {
    const values = readForm(form);
    const prefs = state.prefs;
    const settings = { ...values.settings };
    if (mode === "bisulphite") settings.repeats = prefs.bis_repeats;
    else {
      settings.gc = prefs.gc;
      settings.repeats = prefs.repeats;
    }
    const req = { settings, conditions: prefs.conditions };
    if (mode === "qpcr") {
      Object.assign(req, { mrna: values.mrna, genomic: values.genomic, exon_boundaries: values.exon_boundaries });
    } else if (mode === "sequencing") {
      Object.assign(req, { sequence: values.sequence, start: values.start, end: values.end });
    } else {
      const t = values.target;
      const anySet = Object.values(t).some((v) => v !== null);
      Object.assign(req, { sequence: values.sequence, target: anySet ? t : null, search: values.search });
    }
    return req;
  }

  function initDesignPanel(panel) {
    const mode = panel.dataset.mode;
    const form = $("form", panel);
    const d = state.defaults;
    fillForm(form, { settings: d[mode], search: "outwards" });
    state.results[mode] = { sort: null };

    form.addEventListener("submit", (event) => {
      event.preventDefault();
      runDesign(mode, panel, form, event.submitter);
    });

    $$("[data-action=features]", panel).forEach((btn) =>
      btn.addEventListener("click", () => busy(btn, () => findFeatures(panel, btn.dataset.kind))),
    );
    $('[data-action="clear-target"]', panel)?.addEventListener("click", () => {
      ["outer_start", "inner_start", "inner_end", "outer_end"].forEach((k) => (field(panel, `target.${k}`).value = ""));
    });
    $('[data-action="map-exons"]', panel)?.addEventListener("click", (e) => busy(e.currentTarget, () => mapExons(panel)));
    $('[data-action="export"]', panel).addEventListener("click", () => exportTable(mode));
  }

  async function runDesign(mode, panel, form, button) {
    let request;
    try {
      request = buildRequest(mode, form);
    } catch (err) {
      setStatus(panel, err.message, "error");
      return;
    }
    const result = state.results[mode];
    setStatus(panel, "Searching for primers…", "busy");
    $('[data-role="details"]', panel).hidden = true;
    await busy(button, async () => {
      try {
        const data = await api(`/design/${mode}`, request);
        Object.assign(result, { data, request, rows: data.pairs || data.primers, selected: null, sort: null, enzymes: null });
        if (data.target && mode !== "qpcr" && mode !== "sequencing") fillForm(form, { target: data.target });
        if (mode === "qpcr" && data.exon_boundaries) {
          field(panel, "exon_boundaries").value = data.exon_boundaries.join(", ");
          if (data.exons) renderExons(panel, data.exons);
        }
        setStatus(panel, data.messages.join(" · "), result.rows.length ? "" : "error");
        renderTable(mode, panel);
        renderOverview(mode, panel);
      } catch (err) {
        setStatus(panel, err.message, "error");
      }
    });
  }

  async function findFeatures(panel, kind) {
    const textarea = field(panel, panel.dataset.mode === "qpcr" ? "mrna" : "sequence");
    try {
      const data = await api("/features", { sequence: textarea.value, kind, cpg: state.prefs.cpg });
      if (!data.regions.length) {
        setStatus(panel, kind === "orf" ? "No open reading frame found" : "No CpG island found", "error");
        return;
      }
      const header = textarea.value.trim().startsWith(">") ? textarea.value.trim().split(/\r?\n/)[0] + "\n" : "";
      textarea.value = header + data.highlighted;
      updateLength(textarea);
      const t = data.target;
      if (panel.dataset.mode === "sequencing") {
        fillForm($("form", panel), { start: t.inner_start, end: t.inner_end });
      } else {
        fillForm($("form", panel), { target: t });
      }
      const label = kind === "orf" ? "ORF" : `${data.regions.length} CpG island(s)`;
      const desc = data.regions.map((r) => `${r.start}–${r.end}`).join(", ");
      setStatus(panel, `${label}: ${desc} (shown in upper case)`);
    } catch (err) {
      setStatus(panel, err.message, "error");
    }
  }

  async function mapExons(panel) {
    setStatus(panel, "Aligning mRNA to genomic sequence…", "busy");
    try {
      const data = await api("/exons", { mrna: field(panel, "mrna").value, genomic: field(panel, "genomic").value });
      field(panel, "exon_boundaries").value = data.boundaries.join(", ");
      renderExons(panel, data.exons);
      const tail = data.unaligned_tail ? ` (${data.unaligned_tail} unaligned bases at the 3′ end, e.g. poly-A)` : "";
      setStatus(panel, `Found ${data.exons.length} exons${tail}`);
    } catch (err) {
      setStatus(panel, err.message, "error");
    }
  }

  function renderExons(panel, exons) {
    $('[data-role="exons"]', panel).innerHTML = `<table><thead><tr><th>Exon</th><th>mRNA</th><th>Genomic</th><th>Length</th></tr></thead><tbody>${exons
      .map((e, i) => `<tr><td>${i + 1}</td><td>${e.mrna_start + 1}–${e.mrna_end}</td><td>${e.genomic_start + 1}–${e.genomic_end}</td><td>${e.mrna_end - e.mrna_start}</td></tr>`)
      .join("")}</tbody></table>`;
  }

  // ---------------------------------------------------------- results table

  function columnsFor(mode) {
    return mode === "sequencing" ? SEQ_COLUMNS : PAIR_COLUMNS;
  }

  function renderTable(mode, panel) {
    const result = state.results[mode];
    const rows = result.rows;
    const wrap = $('[data-role="table-wrap"]', panel);
    const actions = $('[data-role="table-actions"]', panel);
    wrap.hidden = actions.hidden = !rows.length;
    if (!rows.length) return;
    const total = result.data.total;
    $('[data-role="count"]', panel).textContent =
      `${total.toLocaleString()} ${mode === "sequencing" ? "primers" : "primer pairs"}` +
      (result.data.truncated ? ` (showing ${rows.length.toLocaleString()})` : "");

    const cols = columnsFor(mode);
    const table = $('[data-role="table"]', panel);
    const sort = result.sort;
    table.innerHTML =
      `<thead><tr>${cols
        .map(([label], i) => {
          const s = sort && sort.col === i ? (sort.dir > 0 ? "ascending" : "descending") : "none";
          return `<th scope="col" data-col="${i}" aria-sort="${s}">${escapeHtml(label)}</th>`;
        })
        .join("")}</tr></thead>` +
      `<tbody>${rows
        .map(
          (row, r) =>
            `<tr data-row="${r}" aria-selected="${result.selected === row}">${cols
              .map(([, key, kind]) => `<td class="${kind === "seq" ? "seq" : "num"}">${escapeHtml(kind === "tm" ? fmt(row[key]) : row[key])}</td>`)
              .join("")}</tr>`,
        )
        .join("")}</tbody>`;

    $$("th", table).forEach((th) =>
      th.addEventListener("click", () => {
        const col = Number(th.dataset.col);
        const key = cols[col][1];
        const dir = result.sort && result.sort.col === col ? -result.sort.dir : 1;
        result.sort = { col, dir };
        result.rows = [...result.rows].sort((a, b) => (a[key] < b[key] ? -dir : a[key] > b[key] ? dir : 0));
        renderTable(mode, panel);
      }),
    );
    $$("tbody tr", table).forEach((tr) =>
      tr.addEventListener("click", () => selectRow(mode, panel, result.rows[Number(tr.dataset.row)])),
    );
  }

  function exportTable(mode) {
    const result = state.results[mode];
    const cols = columnsFor(mode);
    const lines = [cols.map(([label]) => label).join("\t")];
    result.rows.forEach((row) => lines.push(cols.map(([, key]) => row[key]).join("\t")));
    download(`primers_${mode}.tsv`, lines.join("\n") + "\n", "text/tab-separated-values");
  }

  // ----------------------------------------------------------------- details

  function detailsRequest(mode, result) {
    const r = result.request;
    const req = {
      mode,
      sequence: mode === "qpcr" ? r.mrna : r.sequence,
      conditions: state.prefs.conditions,
      cpg: state.prefs.cpg,
      title: document.title,
    };
    if (mode === "sequencing") req.primer = result.selected;
    else req.pair = result.selected;
    if (mode === "pcr" && result.cloning) req.cloning = result.cloning;
    return req;
  }

  async function selectRow(mode, panel, row) {
    const result = state.results[mode];
    result.selected = row;
    $$('[data-role="table"] tbody tr', panel).forEach((tr) =>
      tr.setAttribute("aria-selected", String(result.rows[Number(tr.dataset.row)] === row)),
    );
    renderOverview(mode, panel);
    await refreshDetails(mode, panel);
  }

  async function refreshDetails(mode, panel) {
    const result = state.results[mode];
    const box = $('[data-role="details"]', panel);
    try {
      const data = await api("/details", detailsRequest(mode, result));
      if (mode === "pcr" && !result.enzymes) {
        const t = result.data.target;
        const region = cleanSequence(result.request.sequence).slice(t.outer_start, t.outer_end + 1);
        result.enzymes = await api("/enzymes", { region, exclude_found: true }).catch(() => []);
      }
      box.innerHTML = renderDetails(mode, result, data);
      box.hidden = false;
      bindDetails(mode, panel, box, data);
    } catch (err) {
      box.hidden = false;
      box.innerHTML = `<p class="msg error">${escapeHtml(err.message)}</p>`;
    }
  }

  function primerCard(name, cls, analysis, alt, altLabel) {
    if (!analysis) return "";
    const warnings = [];
    if (analysis.run) warnings.push("Contains a run of identical bases");
    if (analysis.repeat) warnings.push(`Contains a repeat (${escapeHtml(analysis.repeat)})`);
    return `<div class="primer-card">
      <div class="name">${name}</div>
      <div class="seq ${cls}">5′-${escapeHtml(analysis.seq)}-3′</div>
      ${alt ? `<div class="alt" title="${escapeHtml(altLabel)}">${escapeHtml(altLabel)}: 5′-${escapeHtml(alt)}-3′</div>` : ""}
      <dl class="stats">
        <div><dt>Tm</dt><dd>${fmt(analysis.tm)} °C</dd></div>
        <div><dt>Length</dt><dd>${analysis.length} nt</dd></div>
        <div><dt>GC</dt><dd>${fmt(analysis.gc, 0)} %</dd></div>
        <div><dt>ΔH°</dt><dd>${fmt(analysis.dh)} kcal/mol</dd></div>
        <div><dt>ΔS°</dt><dd>${fmt(analysis.ds)} eu</dd></div>
        <div><dt>ΔG°</dt><dd>${fmt(analysis.dg)} kcal/mol</dd></div>
      </dl>
      ${warnings.map((w) => `<p class="warning">${w}</p>`).join("")}
    </div>`;
  }

  function dimerList(title, dimers, emptyText) {
    const items = dimers.length
      ? dimers
          .map((d) => `<div class="dimer ${d.dg < -5 ? "strong" : ""}"><div>${escapeHtml(d.label)}: <span class="dg">${fmt(d.dg)} kcal/mol</span></div><pre>${escapeHtml(d.diagram)}</pre></div>`)
          .join("")
      : `<p class="muted">${emptyText}</p>`;
    return `<div><h3>${title}</h3><div class="dimers">${items}</div></div>`;
  }

  function analysisHtml(analysis, alts = {}, altLabel = "") {
    const temp = state.prefs.conditions.dimer_temp_c;
    return `
      ${analysis.errors.map((e) => `<p class="msg error">${escapeHtml(e)}</p>`).join("")}
      <div class="primer-cards">
        ${primerCard("Forward primer", "fwd", analysis.forward, alts.forward, altLabel)}
        ${primerCard("Reverse primer", "rev", analysis.reverse, alts.reverse, altLabel)}
      </div>
      ${dimerList(`Extensible 3′ primer-dimers at ${temp} °C`, analysis.extensible, "None found")}
      ${dimerList("More stable non-extensible dimers", analysis.non_extensible, "None found")}`;
  }

  function blastUrl(primers) {
    const query = primers.filter(([, s]) => s).map(([n, s]) => `>${n}\n${s}`).join("\n");
    const params = new URLSearchParams({ PAGE_TYPE: "BlastSearch", PROGRAM: "blastn", QUERY: query, WORD_SIZE: "7", EXPECT: "1000" });
    return `https://blast.ncbi.nlm.nih.gov/Blast.cgi?${params}`;
  }

  function renderDetails(mode, result, data) {
    const row = result.selected;
    const title =
      mode === "sequencing"
        ? `Primer at ${row.pos}`
        : `Pair: ${row.forward_pos}–${row.reverse_end} · ${row.amplicon} bp amplicon`;
    const alts = data.cloned || data.degenerate || {};
    const altLabel = data.cloned ? "With cloning tail" : data.degenerate ? "Degenerate (CpG)" : "";
    const blast = blastUrl([["forward", row.forward || row.seq], ["reverse", row.reverse]]);
    return `
      <div class="detail-head">
        <h2>${escapeHtml(title)}</h2><span class="spacer"></span>
        <button type="button" class="btn ghost small" data-action="copy">Copy primers</button>
        <button type="button" class="btn ghost small" data-action="report">Download report</button>
        <a class="btn ghost small" href="${escapeHtml(blast)}" target="_blank" rel="noopener noreferrer">BLAST at NCBI ↗</a>
      </div>
      ${analysisHtml(data.analysis, alts, altLabel)}
      ${mode === "pcr" ? cloningHtml(result) : ""}
      <div>
        <h3>Sequence map</h3>
        <div class="legend">
          <span style="--c: var(--fwd)">Forward primer</span>
          ${mode === "sequencing" ? "" : '<span style="--c: var(--rev)">Reverse primer</span>'}
          ${mode === "pcr" || mode === "bisulphite" ? '<span style="--c: var(--inner)">Inner range</span><span style="--c: var(--outer)">Outer range</span>' : ""}
          <span style="--c: var(--feature-alt)">${mode === "bisulphite" ? "CpG island" : "ORF"}</span>
          ${mode === "bisulphite" ? '<span style="--c: var(--cpg)">CpG</span>' : ""}
          ${mode === "qpcr" ? '<span style="--c: var(--cpg)">Exon junction</span>' : ""}
        </div>
        <div class="seqmap">${renderMap(mode, result, data)}</div>
      </div>`;
  }

  function cloningHtml(result) {
    const c = result.cloning || { anchor: state.defaults.anchor, forward_enzyme: "", reverse_enzyme: "", forward_frame: null, reverse_frame: null };
    const options = (selected) =>
      `<option value="">none</option>` +
      (result.enzymes || [])
        .map((e) => `<option value="${escapeHtml(e.name)}" ${e.name === selected ? "selected" : ""}>${escapeHtml(e.name)} – ${escapeHtml(e.pattern)}</option>`)
        .join("");
    const frames = (selected) =>
      `<option value="">ignore frame</option>` + [0, 1, 2].map((f) => `<option value="${f}" ${selected === f ? "selected" : ""}>${f === 0 ? "in frame" : `+${f}`}</option>`).join("");
    return `<form class="cloning" data-role="cloning">
      <h3>Cloning sites</h3>
      <p class="help">Add restriction sites (six-cutters that do not cut the amplified region) and a 5′ anchor. Choose the frame of the site in the vector to keep a fusion in frame.</p>
      <div class="grid">
        <label class="field"><span>Forward site</span><select name="forward_enzyme">${options(c.forward_enzyme)}</select></label>
        <label class="field"><span>Frame</span><select name="forward_frame">${frames(c.forward_frame)}</select></label>
        <label class="field"><span>Reverse site</span><select name="reverse_enzyme">${options(c.reverse_enzyme)}</select></label>
        <label class="field"><span>Frame</span><select name="reverse_frame">${frames(c.reverse_frame)}</select></label>
        <label class="field"><span>5′ anchor</span><input name="anchor" class="mono" value="${escapeHtml(c.anchor)}"></label>
      </div>
    </form>`;
  }

  function bindDetails(mode, panel, box, data) {
    const result = state.results[mode];
    const row = result.selected;
    $('[data-action="copy"]', box).addEventListener("click", async (e) => {
      const alts = data.cloned || data.degenerate || {};
      const fwd = alts.forward || row.forward || row.seq;
      const rev = alts.reverse || row.reverse;
      const text = rev ? `${fwd}\t${rev}` : fwd;
      try {
        await navigator.clipboard.writeText(text);
        e.target.textContent = "Copied";
      } catch {
        window.prompt("Copy the primers:", text);
      }
    });
    $('[data-action="report"]', box).addEventListener("click", (e) =>
      busy(e.target, async () => {
        try {
          const response = await api("/report", detailsRequest(mode, result), { raw: true });
          const name = (response.headers.get("Content-Disposition") || "").match(/filename="([^"]+)"/);
          download(name ? name[1] : "primer_report.txt", await response.text());
        } catch (err) {
          setStatus(panel, err.message, "error");
        }
      }),
    );
    const cloning = $('[data-role="cloning"]', box);
    if (cloning) {
      cloning.addEventListener("change", () => {
        const f = new FormData(cloning);
        const frame = (v) => (v === "" ? null : Number(v));
        result.cloning = {
          forward_enzyme: f.get("forward_enzyme"),
          reverse_enzyme: f.get("reverse_enzyme"),
          anchor: f.get("anchor").toUpperCase().replace(/[^ACGT]/g, ""),
          forward_frame: frame(f.get("forward_frame")),
          reverse_frame: frame(f.get("reverse_frame")),
        };
        refreshDetails(mode, panel);
      });
      cloning.addEventListener("submit", (e) => e.preventDefault());
    }
  }

  // ------------------------------------------------------------ sequence map

  function renderMap(mode, result, data) {
    const m = data.map;
    const seq = m.sequence;
    const n = seq.length;
    const cls = Array.from({ length: n }, () => []);
    const mark = (start, end, c) => {
      for (let i = Math.max(0, start); i < Math.min(n, end); i++) cls[i].push(c);
    };
    data.regions.forEach((r) => mark(r.start, r.end, "s-feature"));
    const t = result.data.target;
    if (t && t.outer_start != null && (mode === "pcr" || mode === "bisulphite")) {
      mark(t.outer_start, t.outer_end + 1, "s-outer");
      mark(t.inner_start, t.inner_end + 1, "s-inner");
    }
    if (mode === "bisulphite") {
      for (let i = 0; i < n - 1; i++) if ((seq[i] === "C" || seq[i] === "c") && (seq[i + 1] === "G" || seq[i + 1] === "g")) mark(i, i + 2, "s-mark");
    }
    if (mode === "qpcr") (result.data.exon_boundaries || []).forEach((b) => mark(b - 1, b + 1, "s-mark"));
    const row = result.selected;
    if (mode === "sequencing") mark(row.pos, row.pos + row.length, "s-fwd");
    else {
      mark(row.forward_pos, row.forward_pos + row.forward_len, "s-fwd");
      mark(row.reverse_end - row.reverse_len + 1, row.reverse_end + 1, "s-rev");
    }

    const priority = ["s-fwd", "s-rev", "s-mark", "s-inner", "s-outer"];
    const classOf = (i) => {
      const c = cls[i];
      const main = priority.find((p) => c.includes(p)) || "";
      return [main, c.includes("s-feature") && main !== "s-fwd" && main !== "s-rev" && main !== "s-mark" ? "s-feature" : ""].filter(Boolean).join(" ");
    };

    const blocks = [];
    for (let start = 0; start < n; start += MAP_WIDTH) {
      const end = Math.min(start + MAP_WIDTH, n);
      let seqHtml = "";
      let runStart = start;
      for (let i = start + 1; i <= end; i++) {
        if (i === end || classOf(i) !== classOf(runStart)) {
          const c = classOf(runStart);
          const text = escapeHtml(seq.slice(runStart, i));
          seqHtml += c ? `<span class="${c}">${text}</span>` : text;
          runStart = i;
        }
      }
      const slice = (line) => escapeHtml((line || "").slice(start, end).replace(/\s+$/, ""));
      blocks.push(
        `<div class="block"><div class="ruler">${slice(m.ruler)}</div>` +
          `<div class="primers">${slice(m.primers) || " "}</div>` +
          `<div>${seqHtml}</div>` +
          `<div class="feat">${slice(m.features) || " "}</div></div>`,
      );
    }
    return blocks.join("");
  }

  // --------------------------------------------------------------- overview

  function renderOverview(mode, panel) {
    const svg = $('[data-role="overview"]', panel);
    const result = state.results[mode];
    const n = result.data && result.data.length;
    if (!n) {
      svg.toggleAttribute("hidden", true);
      return;
    }
    svg.toggleAttribute("hidden", false);
    const W = svg.clientWidth || 800;
    svg.setAttribute("viewBox", `0 0 ${W} 62`);
    const pad = 14;
    const x = (pos) => pad + (Math.max(0, Math.min(pos, n)) / n) * (W - 2 * pad);
    const parts = [`<rect class="dna" x="${pad}" y="26" width="${W - 2 * pad}" height="10" rx="3"/>`];
    const t = result.data.target;
    if (t && t.outer_start != null && (mode === "pcr" || mode === "bisulphite")) {
      parts.push(`<rect class="outer" x="${x(t.outer_start)}" y="20" width="${Math.max(1, x(t.outer_end) - x(t.outer_start))}" height="22" rx="3"/>`);
      parts.push(`<rect class="inner" x="${x(t.inner_start)}" y="23" width="${Math.max(1, x(t.inner_end) - x(t.inner_start))}" height="16" rx="2"/>`);
    }
    if (mode === "qpcr") (result.data.exon_boundaries || []).forEach((b) => parts.push(`<line class="boundary" x1="${x(b)}" x2="${x(b)}" y1="18" y2="44"/>`));
    if (mode === "sequencing") {
      result.rows.forEach((p) => parts.push(`<path class="fwd" d="M${x(p.pos)} 14 l7 4 l-7 4 z" opacity="${p === result.selected ? 1 : 0.45}"/>`));
    } else if (result.selected) {
      const p = result.selected;
      parts.push(`<path class="fwd" d="M${x(p.forward_pos)} 12 h${Math.max(4, x(p.forward_pos + p.forward_len) - x(p.forward_pos))} l5 5 h-${Math.max(4, x(p.forward_pos + p.forward_len) - x(p.forward_pos)) + 5} z"/>`);
      const rs = x(p.reverse_end - p.reverse_len + 1);
      const rw = Math.max(4, x(p.reverse_end) - rs);
      parts.push(`<path class="rev" d="M${rs - 5} 45 h${rw + 5} v5 h-${rw} z"/>`);
    }
    parts.push(`<text x="${pad}" y="58">0</text><text x="${W - pad}" y="58" text-anchor="end">${n.toLocaleString()} bp</text>`);
    svg.innerHTML = parts.join("");
  }

  // ---------------------------------------------------------- analysis tab

  function initAnalysis() {
    const panel = $("#tab-analyse");
    const form = $("form", panel);
    form.addEventListener("submit", async (event) => {
      event.preventDefault();
      const values = readForm(form);
      const box = $('[data-role="details"]', panel);
      setStatus(panel, "");
      await busy(event.submitter, async () => {
        try {
          const data = await api("/analyse", { ...values, conditions: state.prefs.conditions, repeats: state.prefs.repeats });
          const blast = blastUrl([["forward", data.forward?.seq], ["reverse", data.reverse?.seq]]);
          box.innerHTML = `<div class="detail-head"><h2>Primer analysis</h2><span class="spacer"></span>
            <a class="btn ghost small" href="${escapeHtml(blast)}" target="_blank" rel="noopener noreferrer">BLAST at NCBI ↗</a></div>
            ${analysisHtml(data)}`;
          box.hidden = false;
        } catch (err) {
          box.hidden = true;
          setStatus(panel, err.message, "error");
        }
      });
    });
  }

  // ---------------------------------------------------------- AutoDimer tab

  const AUTODIMER_KEY = "pyprimer.autodimer.v1";
  const AUTODIMER_COLUMNS = [
    // label, key, numeric
    ["Oligos", "labels", false],
    ["Score", "score", true],
    ["Matches", "matches", true],
    ["Tm (°C)", "tm", true],
    ["ΔG (kcal/mol)", "dg", true],
    ["Alignment", "lineup", false],
  ];
  const AUTODIMER_KINDS = {
    dimers: { name: "primer-dimers", scorePath: "settings.min_score", lower: (s) => s - 1 },
    // Hairpin stems are short, so drop straight to a typical hairpin score.
    hairpins: { name: "hairpins", scorePath: "settings.hairpin_min_score", lower: (s) => Math.min(s - 1, 3) },
  };

  function countOligos(text) {
    const lines = text.split(/\r?\n/).map((l) => l.trim()).filter(Boolean);
    const labels = lines.filter((l) => l.startsWith(">")).length;
    return labels || lines.length;
  }

  function initAutoDimer() {
    const panel = $("#tab-autodimer");
    const form = $("form", panel);
    const textarea = field(panel, "sequences");
    const counter = $('[data-role="oligo-count"]', panel);
    let result = null; // { data, settings, sections: { dimers: {rows, sort}, hairpins: {rows, sort} } }

    const updateCount = () => {
      const n = countOligos(textarea.value);
      counter.textContent = `${n} sequence${n === 1 ? "" : "s"}`;
    };
    const save = () => storageSet(AUTODIMER_KEY, readForm(form));
    fillForm(form, { sequences: "", settings: state.defaults.autodimer });
    const saved = storageGet(AUTODIMER_KEY);
    if (saved && typeof saved === "object") fillForm(form, saved);
    updateCount();
    textarea.addEventListener("input", () => (updateCount(), save()));
    form.addEventListener("change", save);
    $('[data-action="autodimer-defaults"]', panel).addEventListener("click", () => {
      fillForm(form, { settings: state.defaults.autodimer });
      save();
    });

    const value = (row, key) => (key === "labels" ? row.label ?? `${row.label_1} vs ${row.label_2}` : row[key]);
    const shown = (v, missing) => (v === null ? missing : v.toFixed(1));
    const labelsHtml = (row) => {
      if (row.label !== undefined) return escapeHtml(row.label);
      if (row.self_dimer) return `${escapeHtml(row.label_1)}<br><span class="vs">self-dimer</span>`;
      return `${escapeHtml(row.label_1)}<br><span class="vs">vs</span> ${escapeHtml(row.label_2)}`;
    };

    function renderSection(kind) {
      const section = $(`[data-section="${kind}"]`, panel);
      const { rows, sort } = result.sections[kind];
      const { comparisons } = result.data[kind];
      const score = getPath(result.settings, AUTODIMER_KINDS[kind].scorePath);
      section.hidden = false;
      $('[data-role="section-count"]', section).textContent =
        `${rows.length.toLocaleString()} hit${rows.length === 1 ? "" : "s"} · score ≥ ${score} · ` +
        `${comparisons.toLocaleString()} ${kind === "dimers" ? "pairs" : "oligos"} screened`;
      const empty = $('[data-role="section-empty"]', section);
      const wrap = $('[data-role="table-wrap"]', section);
      empty.hidden = rows.length > 0;
      wrap.hidden = !rows.length;
      if (!rows.length) {
        empty.innerHTML = `No ${AUTODIMER_KINDS[kind].name} found at score ${escapeHtml(score)}. <a href="#" data-action="lower-score">Try a lower score</a>`;
        $('[data-action="lower-score"]', empty).addEventListener("click", (event) => {
          event.preventDefault();
          const input = field(panel, AUTODIMER_KINDS[kind].scorePath);
          input.value = AUTODIMER_KINDS[kind].lower(score);
          save();
          form.requestSubmit();
        });
        return;
      }
      const table = $('[data-role="table"]', section);
      table.innerHTML =
        `<thead><tr>${AUTODIMER_COLUMNS.map(([label, , numeric], i) => {
          const s = sort && sort.col === i ? (sort.dir > 0 ? "ascending" : "descending") : "none";
          const text = i === 4 ? `ΔG at ${result.settings.settings.dg_temp_c} °C (kcal/mol)` : label;
          return `<th scope="col" data-col="${i}" aria-sort="${s}" ${numeric ? "" : 'class="text"'}>${escapeHtml(text)}</th>`;
        }).join("")}</tr></thead>` +
        `<tbody>${rows
          .map(
            (row) => `<tr>
              <td class="labels">${labelsHtml(row)}</td>
              <td class="num">${row.score}</td>
              <td class="num">${row.matches}</td>
              <td class="num">${shown(row.tm, "< 0")}</td>
              <td class="num">${shown(row.dg, "> 0")}</td>
              <td class="seq"><pre>${escapeHtml(row.lineup)}</pre></td>
            </tr>`,
          )
          .join("")}</tbody>`;
      $$("th", table).forEach((th) =>
        th.addEventListener("click", () => {
          const col = Number(th.dataset.col);
          const key = AUTODIMER_COLUMNS[col][1];
          const current = result.sections[kind];
          const dir = current.sort && current.sort.col === col ? -current.sort.dir : 1;
          // "< 0" and "> 0" (null) sort as the least and most stable values.
          const missing = key === "tm" ? -Infinity : key === "dg" ? Infinity : "";
          const v = (row) => value(row, key) ?? missing;
          current.sort = { col, dir };
          current.rows = [...current.rows].sort((a, b) => (v(a) < v(b) ? -dir : v(a) > v(b) ? dir : 0));
          renderSection(kind);
        }),
      );
    }

    function render() {
      $('[data-role="table-actions"]', panel).hidden = false;
      $('[data-role="count"]', panel).textContent =
        `${result.data.sequences} sequence${result.data.sequences === 1 ? "" : "s"} screened`;
      Object.keys(AUTODIMER_KINDS).forEach(renderSection);
    }

    form.addEventListener("submit", async (event) => {
      event.preventDefault();
      setStatus(panel, "");
      $('[data-role="table-actions"]', panel).hidden = true;
      $$("[data-section]", panel).forEach((s) => (s.hidden = true));
      const button = $("button[type=submit]", form);
      await busy(button, async () => {
        try {
          const values = readForm(form);
          setStatus(panel, "Screening…", "busy");
          const data = await api("/autodimer", values);
          setStatus(panel, "");
          result = {
            data,
            settings: values,
            sections: {
              dimers: { rows: data.dimers.hits, sort: null },
              hairpins: { rows: data.hairpins.hits, sort: null },
            },
          };
          render();
        } catch (err) {
          setStatus(panel, err.message, "error");
        }
        revealResults();
      });
    });

    // The button sits below a tall form, so bring the results (at the top of
    // their column, possibly under the sticky header) into view.
    function revealResults() {
      const box = $(".results", panel);
      const header = $(".topbar").offsetHeight;
      const top = box.getBoundingClientRect().top;
      if (top < header || top > window.innerHeight - 60) {
        window.scrollTo(0, Math.max(0, window.scrollY + top - header - 12));
      }
    }

    $('[data-action="save-text"]', panel).addEventListener("click", () => {
      if (result) download("autodimer.txt", result.data.text);
    });
  }

  // ---------------------------------------------------------- sequence boxes

  function initSequenceBoxes() {
    $$("textarea[data-path]").forEach((ta) => ta.addEventListener("input", () => updateLength(ta)));
    $$("input[type=file][data-load-into]").forEach((input) =>
      input.addEventListener("change", async () => {
        const file = input.files[0];
        if (!file) return;
        const ta = field(input.closest(".panel"), input.dataset.loadInto);
        ta.value = await file.text();
        ta.dispatchEvent(new Event("input"));
        input.value = "";
      }),
    );
    $$("[data-action=revcomp]").forEach((btn) =>
      btn.addEventListener("click", () => {
        const ta = field(btn.closest(".panel"), btn.dataset.target);
        const seq = cleanSequence(ta.value);
        const header = ta.value.trim().startsWith(">") ? ta.value.trim().split(/\r?\n/)[0] + " (reverse complement)\n" : "";
        ta.value = header + reverseComplement(seq);
        updateLength(ta);
      }),
    );
  }

  // ----------------------------------------------------------------- Ensembl

  function initEnsembl() {
    const dialog = $("#ensembl-dialog");
    const status = (msg, kind = "") => {
      $('[data-role="ensembl-status"]', dialog).innerHTML = msg ? `<span class="msg ${kind}">${escapeHtml(msg)}</span>` : "";
    };
    let target = null;

    $$("[data-action=ensembl]").forEach((btn) =>
      btn.addEventListener("click", () => {
        target = { panel: btn.closest(".panel"), path: btn.dataset.target };
        $('[data-role="ensembl-gene"]', dialog).innerHTML = "";
        status("");
        dialog.showModal();
        $("input[name=q]", dialog).focus();
      }),
    );

    const search = async () => {
      const species = $("input[name=species]", dialog).value;
      const q = $("input[name=q]", dialog).value;
      status("Searching Ensembl…", "busy");
      try {
        const gene = await api(`/ensembl/gene?${new URLSearchParams({ species, q })}`);
        status("");
        renderGene(gene);
      } catch (err) {
        status(err.message, "error");
      }
    };
    $('[data-action="ensembl-search"]', dialog).addEventListener("click", search);
    $("input[name=q]", dialog).addEventListener("keydown", (e) => {
      if (e.key === "Enter") {
        e.preventDefault();
        search();
      }
    });

    function renderGene(gene) {
      const qpcr = target.panel.dataset.mode === "qpcr";
      const box = $('[data-role="ensembl-gene"]', dialog);
      box.innerHTML = `
        <p><strong>${escapeHtml(gene.name)}</strong> · ${escapeHtml(gene.id)} · ${escapeHtml(gene.location)}<br>
        <span class="muted">${escapeHtml(gene.description)}</span></p>
        <div class="transcripts">${gene.transcripts
          .map(
            (t, i) => `<label><input type="radio" name="transcript" value="${i}" ${i === 0 ? "checked" : ""}>
              <span>${escapeHtml(t.name)} <span class="muted">${escapeHtml(t.id)} · ${escapeHtml(t.biotype)} · ${t.exon_lengths.length} exons</span></span>
              ${t.canonical ? '<span class="badge">canonical</span>' : ""}</label>`,
          )
          .join("")}</div>
        <div class="row">
          ${qpcr ? "" : `<label class="field"><span>Sequence</span><select name="seqtype">
            <option value="cdna">cDNA (transcript)</option><option value="cds">CDS (transcript)</option>
            <option value="genomic">Genomic (gene)</option></select></label>`}
          <span class="spacer"></span>
          <button type="button" class="btn primary" data-action="ensembl-load">${qpcr ? "Load cDNA and exon junctions" : "Load sequence"}</button>
        </div>`;
      $('[data-action="ensembl-load"]', box).addEventListener("click", (e) =>
        busy(e.target, async () => {
          const t = gene.transcripts[Number($("input[name=transcript]:checked", box).value)];
          const type = qpcr ? "cdna" : $("select[name=seqtype]", box).value;
          const id = type === "genomic" ? gene.id : t.id;
          status("Downloading sequence…", "busy");
          try {
            const data = await api(`/ensembl/sequence?${new URLSearchParams({ id, type })}`);
            const ta = field(target.panel, target.path);
            ta.value = `>${type === "genomic" ? gene.name : t.name} ${id} ${type}\n${data.sequence}`;
            updateLength(ta);
            if (qpcr) field(target.panel, "exon_boundaries").value = t.exon_boundaries.join(", ");
            dialog.close();
          } catch (err) {
            status(err.message, "error");
          }
        }),
      );
    }
  }

  // ------------------------------------------------------------------- init

  async function init() {
    initTabs();
    try {
      state.defaults = await api("/defaults");
    } catch (err) {
      document.querySelector("main").insertAdjacentHTML("afterbegin", `<p class="msg error">Could not load settings: ${escapeHtml(err.message)}</p>`);
      return;
    }
    state.prefs = loadPrefs();
    initPrefs();
    $$(".panel[data-mode]").forEach((panel) => {
      if (["pcr", "bisulphite", "qpcr", "sequencing"].includes(panel.dataset.mode)) initDesignPanel(panel);
    });
    initAnalysis();
    initSequenceBoxes();
    initAutoDimer();
    initEnsembl();
    window.addEventListener("resize", () => {
      for (const mode of Object.keys(state.results)) {
        const panel = $(`#tab-${mode}`);
        if (!panel.hidden) renderOverview(mode, panel);
      }
    });
  }

  init();
})();
