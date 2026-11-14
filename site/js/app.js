// tcga_survivome: cohort builder + survival analysis UI.
import { initDb, query, ident, lit } from './db.js';
import { kmFit, kmMedian, survAt, medianFollowUp, logRank, coxFit, quantileSorted } from './stats.js';
import { renderKM } from './kmplot.js';

const $ = (id) => document.getElementById(id);
function h(tag, props = {}, ...kids) {
  const e = document.createElement(tag);
  for (const [k, v] of Object.entries(props)) {
    if (k === 'class') e.className = v;
    else if (k.startsWith('on')) e.addEventListener(k.slice(2), v);
    else if (v !== false && v != null) e.setAttribute(k, v === true ? '' : v);
  }
  for (const kid of kids.flat()) if (kid != null) e.append(kid);
  return e;
}

const GROUP_LABELS = {
  clinical: 'Clinical', therapy: 'Therapy', mutation: 'Mutation burden and MSI', dna_repair: 'DNA repair',
  chromosomal_instability: 'Chromosomal instability', tumor: 'Tumor', subtype: 'Subtypes', immune: 'Immune',
  viral: 'Viral', qc: 'Quality review', arm_level_cna: 'Arm-level copy number',
};
const ENDPOINTS = {
  OS: 'Overall survival', DSS: 'Disease-specific survival',
  PFI: 'Progression-free interval', DFI: 'Disease-free interval',
};
const COLORS = [1, 2, 3, 4, 5, 6, 7, 8].map((i) => `var(--s${i})`);
const MAX_GROUPS = COLORS.length;
const MIN_LEVEL_N = 10; // smallest category kept in a Cox model
const collator = new Intl.Collator('en', { numeric: true });

let dict = {};       // column name -> dictionary entry
let variables = [];  // columns that can be filtered / compared
const levelCache = new Map();

const defaults = () => ({
  filters: [{ col: 'cancer_type', values: [] }],
  excludeQC: true, endpoint: 'OS', groupBy: 'stage', split: 'median', cutoff: '', withinType: false,
  unit: 'months', censorAt: '', stratify: false, showCI: true, showCensor: true,
  coxGroups: true, covs: [],
});
let state = defaults();

// ---------- formatting ----------

const isNum = (col) => dict[col].type === 'numeric';
const sqlExpr = (col) => (dict[col].type === 'boolean' ? `CAST(${ident(col)} AS VARCHAR)` : ident(col));
function fmtLevel(col, v) {
  if (v == null) return '(missing)';
  if (dict[col].group === 'arm_level_cna') return { '-1': 'Loss', 0: 'Neutral', 1: 'Gain' }[v] ?? v;
  if (dict[col].type === 'boolean') return v === 'true' ? 'Yes' : 'No';
  return String(v);
}
function fmtNum(v) {
  if (v == null || Number.isNaN(v)) return '–';
  if (v === 0) return '0';
  const a = Math.abs(v);
  if (a >= 1000) return Math.round(v).toLocaleString('en');
  return String(+v.toPrecision(3));
}
function fmtP(p) {
  if (p == null || Number.isNaN(p)) return '–';
  if (p < 1e-4) return p.toExponential(1).replace('e-', ' × 10⁻').replace(/\d+$/, (d) => [...d].map((c) => '⁰¹²³⁴⁵⁶⁷⁸⁹'[c]).join(''));
  if (p < 0.001) return p.toFixed(4);
  return p.toFixed(3);
}
const pct = (v) => (Number.isNaN(v) ? '–' : `${Math.round(v * 100)}%`);
const fmtHR = (v) => (!Number.isFinite(v) ? '–' : v >= 100 ? '>100' : v < 0.01 ? '<0.01' : v.toFixed(2));

// ---------- URL state ----------

function saveState() {
  history.replaceState(null, '', '#' + encodeURIComponent(JSON.stringify(state, (k, v) => (k === 'ui' ? undefined : v))));
}
function loadState() {
  if (location.hash.length < 2) return;
  try {
    const s = JSON.parse(decodeURIComponent(location.hash.slice(1)));
    const ok = (c) => c === '' || variables.some((v) => v.column === c);
    state = { ...defaults(), ...s };
    state.filters = (s.filters || []).filter((f) => ok(f.col) && f.col);
    state.covs = (s.covs || []).filter((c) => ok(c.col) && c.col);
    if (!ok(state.groupBy)) state.groupBy = '';
    if (!ENDPOINTS[state.endpoint]) state.endpoint = 'OS';
  } catch { state = defaults(); }
}

// ---------- cohort SQL ----------

function filterSql(f) {
  if (isNum(f.col)) {
    const parts = [];
    if (f.min !== '' && f.min != null && Number.isFinite(+f.min)) parts.push(`${ident(f.col)} >= ${+f.min}`);
    if (f.max !== '' && f.max != null && Number.isFinite(+f.max)) parts.push(`${ident(f.col)} <= ${+f.max}`);
    return parts.join(' AND ');
  }
  if (!f.values?.length) return '';
  const e = sqlExpr(f.col), named = f.values.filter((v) => v != null), parts = [];
  if (named.length) parts.push(`${e} IN (${named.map(lit).join(', ')})`);
  if (named.length < f.values.length) parts.push(`${e} IS NULL`);
  return `(${parts.join(' OR ')})`;
}
function whereSql(skip = -1) {
  const parts = [];
  if (state.excludeQC) parts.push('NOT qc_do_not_use AND NOT qc_pathology_excluded AND NOT redacted');
  state.filters.forEach((f, i) => { if (i !== skip) { const s = filterSql(f); if (s) parts.push(s); } });
  return parts.length ? parts.join(' AND ') : 'TRUE';
}

// ---------- controls ----------

function fillVarSelect(sel, placeholder) {
  sel.replaceChildren(h('option', { value: '' }, placeholder));
  const groups = new Map();
  for (const v of variables) {
    if (!groups.has(v.group)) groups.set(v.group, h('optgroup', { label: GROUP_LABELS[v.group] || v.group }));
    groups.get(v.group).append(h('option', { value: v.column }, v.label));
  }
  sel.append(...groups.values());
}

async function levelsOf(col) {
  if (!levelCache.has(col)) {
    const rows = await query(`SELECT DISTINCT ${sqlExpr(col)} AS v FROM s`);
    const vals = rows.map((r) => r.v).filter((v) => v != null).sort((a, b) => collator.compare(fmtLevel(col, a), fmtLevel(col, b)));
    if (rows.some((r) => r.v == null)) vals.push(null);
    levelCache.set(col, vals);
  }
  return levelCache.get(col);
}

async function renderFilters() {
  const list = $('filter-list');
  list.replaceChildren();
  for (const [i, f] of state.filters.entries()) {
    const d = dict[f.col];
    const card = h('div', { class: 'filter' },
      h('div', { class: 'filter-head' }, h('b', {}, d.label),
        h('button', { type: 'button', title: 'Remove filter', 'aria-label': `Remove ${d.label} filter`,
          onclick: () => { state.filters.splice(i, 1); renderFilters().then(update); } }, '×')));
    f.ui = {};
    if (isNum(f.col)) {
      f.ui.sub = h('p', { class: 'filter-sub' });
      const mk = (key, ph) => h('input', { type: 'number', step: 'any', placeholder: ph, value: f[key] ?? '', 'aria-label': `${d.label} ${ph}`,
        onchange: (ev) => { f[key] = ev.target.value; update(); } });
      card.append(f.ui.sub, h('div', { class: 'range' }, mk('min', 'min'), '–', mk('max', 'max')));
    } else {
      const levels = await levelsOf(f.col);
      f.values ??= [];
      f.ui.counts = new Map();
      const box = h('div', { class: 'levels' });
      for (const v of levels) {
        const ct = h('span', { class: 'ct' });
        f.ui.counts.set(v, ct);
        box.append(h('label', { 'data-name': fmtLevel(f.col, v).toLowerCase() },
          h('input', { type: 'checkbox', checked: f.values.includes(v),
            onchange: (ev) => {
              f.values = f.values.filter((x) => x !== v);
              if (ev.target.checked) f.values.push(v);
              update();
            } }),
          h('span', { class: 'lv', title: fmtLevel(f.col, v) }, fmtLevel(f.col, v)), ct));
      }
      f.ui.sub = h('p', { class: 'filter-sub' });
      card.append(f.ui.sub);
      if (levels.length > 12) {
        card.append(h('input', { type: 'search', placeholder: 'Search…', 'aria-label': `Search ${d.label}`,
          oninput: (ev) => {
            const q = ev.target.value.trim().toLowerCase();
            for (const row of box.children) row.hidden = !!q && !row.dataset.name.includes(q);
          } }));
      }
      card.append(box);
    }
    list.append(card);
  }
}

async function refreshFilterCounts() {
  await Promise.all(state.filters.map(async (f, i) => {
    const w = whereSql(i);
    if (isNum(f.col)) {
      const [r] = await query(`SELECT min(${ident(f.col)}) lo, max(${ident(f.col)}) hi, count(${ident(f.col)}) n FROM s WHERE ${w}`);
      const unit = dict[f.col].unit ? ` ${dict[f.col].unit}` : '';
      f.ui.sub.textContent = r.n ? `${fmtNum(r.lo)} to ${fmtNum(r.hi)}${unit} · ${r.n.toLocaleString('en')} with a value` : 'No values in this cohort';
    } else {
      const rows = await query(`SELECT ${sqlExpr(f.col)} AS v, count(*) n FROM s WHERE ${w} GROUP BY 1`);
      const m = new Map(rows.map((r) => [r.v, r.n]));
      for (const [v, span] of f.ui.counts) {
        const n = m.get(v) || 0;
        span.textContent = n.toLocaleString('en');
        span.parentElement.style.opacity = n ? 1 : 0.45;
      }
      f.ui.sub.textContent = f.values.length ? `${f.values.length} selected` : 'Any (tick to restrict)';
    }
  }));
}

function renderCovs() {
  $('cov-list').replaceChildren(...state.covs.map((c, i) => {
    const chip = h('span', { class: 'chip' }, dict[c.col].label);
    if (isNum(c.col)) {
      const sel = h('select', { 'aria-label': `Scale for ${dict[c.col].label}`, onchange: (ev) => { c.scale = ev.target.value; update(); } },
        h('option', { value: 'unit' }, `per ${dict[c.col].unit ? dict[c.col].unit.replace(/s$/, '') : 'unit'}`),
        h('option', { value: 'sd' }, 'per SD'),
        h('option', { value: 'log2' }, 'per doubling'));
      sel.value = c.scale || 'unit';
      chip.append(sel);
    }
    chip.append(h('button', { type: 'button', 'aria-label': `Remove ${dict[c.col].label}`,
      onclick: () => { state.covs.splice(i, 1); renderCovs(); update(); } }, '×'));
    return chip;
  }));
}

function syncControls() {
  $('exclude-qc').checked = state.excludeQC;
  $('endpoint').value = state.endpoint;
  $('group-by').value = state.groupBy;
  $('split').value = state.split;
  $('cutoff').value = state.cutoff;
  $('within-type').checked = state.withinType;
  $('unit').value = state.unit;
  $('censor-at').value = state.censorAt;
  $('stratify').checked = state.stratify;
  $('show-ci').checked = state.showCI;
  $('show-censor').checked = state.showCensor;
  $('cox-groups').checked = state.coxGroups;
  const numeric = !!state.groupBy && isNum(state.groupBy);
  $('split-wrap').hidden = !numeric;
  $('cutoff-wrap').hidden = !(numeric && state.split === 'custom');
  $('within-wrap').hidden = !(numeric && state.split !== 'custom');
}

function bind(id, key, { checkbox = false, event = 'change' } = {}) {
  $(id).addEventListener(event, (ev) => {
    state[key] = checkbox ? ev.target.checked : ev.target.value;
    syncControls();
    update();
  });
}

// ---------- grouping ----------

/** Assigns each row a group code (-1 = not shown). Returns { labels, codes, notes }. */
function makeGroups(rows) {
  const col = state.groupBy, notes = [];
  if (!col) return { labels: ['All patients'], codes: rows.map(() => 0), notes };
  const name = dict[col].label;
  const missing = rows.filter((r) => r.g == null).length;
  if (missing) notes.push(`${missing.toLocaleString('en')} patients have no value for ${name} and are not shown.`);

  if (!isNum(col)) {
    const counts = new Map();
    for (const r of rows) if (r.g != null) counts.set(r.g, (counts.get(r.g) || 0) + 1);
    let levels = [...counts.keys()];
    if (levels.length > MAX_GROUPS) {
      levels.sort((a, b) => counts.get(b) - counts.get(a));
      const dropped = levels.slice(MAX_GROUPS).reduce((s, l) => s + counts.get(l), 0);
      notes.push(`${name} has ${levels.length} categories here. The ${MAX_GROUPS} largest are shown; ${dropped.toLocaleString('en')} patients in smaller ones are not. Filter the cohort to compare others.`);
      levels = levels.slice(0, MAX_GROUPS);
    }
    levels.sort((a, b) => collator.compare(fmtLevel(col, a), fmtLevel(col, b)));
    const code = new Map(levels.map((l, i) => [l, i]));
    return { labels: levels.map((l) => fmtLevel(col, l)), codes: rows.map((r) => code.get(r.g) ?? -1), notes };
  }

  // numeric: cut into bins
  const split = state.split;
  const qs = { median: [0.5], tertile: [1 / 3, 2 / 3], quartile: [0.25, 0.5, 0.75], extremes: [0.25, 0.75] }[split];
  const cutsFor = (vals) => {
    const sorted = vals.slice().sort((a, b) => a - b);
    return qs.map((q) => quantileSorted(sorted, q));
  };
  let cutOf, shown = null; // shown: cut values to print in labels
  if (split === 'custom') {
    let c = parseFloat(state.cutoff);
    if (!Number.isFinite(c)) {
      c = quantileSorted(rows.filter((r) => r.g != null).map((r) => r.g).sort((a, b) => a - b), 0.5);
      notes.push(`Enter a cut point for ${name}. Using the median (${fmtNum(c)}) for now.`);
    }
    shown = [c]; cutOf = () => shown;
  } else if (state.withinType) {
    const byType = new Map();
    for (const r of rows) if (r.g != null) { if (!byType.has(r.ct)) byType.set(r.ct, []); byType.get(r.ct).push(r.g); }
    const cuts = new Map([...byType].map(([k, v]) => [k, cutsFor(v)]));
    cutOf = (r) => cuts.get(r.ct);
  } else {
    shown = cutsFor(rows.filter((r) => r.g != null).map((r) => r.g));
    cutOf = () => shown;
  }
  const nb = (split === 'custom' ? 1 : qs.length) + 1;
  const c = shown ? shown.map(fmtNum) : null;
  const rel = ' (within type)';
  let labels;
  if (nb === 2) labels = c ? [`Low (≤ ${c[0]})`, `High (> ${c[0]})`] : ['Low (≤ type median)', 'High (> type median)'];
  else if (split === 'tertile') labels = c ? [`Low (≤ ${c[0]})`, `Middle (${c[0]} to ${c[1]})`, `High (> ${c[1]})`] : ['Low third', 'Middle third', 'High third'].map((s) => s + rel);
  else if (split === 'quartile') labels = c ? [`Q1 (≤ ${c[0]})`, `Q2 (${c[0]} to ${c[1]})`, `Q3 (${c[1]} to ${c[2]})`, `Q4 (> ${c[2]})`] : ['Q1 (lowest)', 'Q2', 'Q3', 'Q4 (highest)'].map((s) => s + rel);
  else labels = c ? [`Bottom quarter (≤ ${c[0]})`, null, `Top quarter (> ${c[1]})`] : [`Bottom quarter${rel}`, null, `Top quarter${rel}`];

  let codes = rows.map((r) => {
    if (r.g == null) return -1;
    const cuts = cutOf(r);
    let b = 0;
    for (const cut of cuts) if (r.g > cut) b++;
    return labels[b] == null ? -1 : b;
  });
  // drop empty bins (heavy ties, e.g. a median of 0) and renumber
  const used = [...new Set(codes.filter((x) => x >= 0))].sort((a, b) => a - b);
  if (used.length < labels.filter((l) => l != null).length) {
    notes.push(`Many patients share the same ${name} value, so some groups are empty.`);
  }
  const remap = new Map(used.map((u, i) => [u, i]));
  codes = codes.map((x) => (x < 0 ? -1 : remap.get(x)));
  return { labels: used.map((u) => labels[u]), codes, notes };
}

// ---------- analysis ----------

let runId = 0;
let lastSvg = null;

async function update() {
  const id = ++runId;
  saveState();
  const where = whereSql();
  const ep = state.endpoint;
  const covCols = state.covs.map((c, i) => `, ${sqlExpr(c.col)} AS c${i}`).join('');
  const gSel = state.groupBy ? `, ${sqlExpr(state.groupBy)} AS g` : '';
  const [[{ n: cohortN }], rows] = await Promise.all([
    query(`SELECT count(*) n FROM s WHERE ${where}`),
    query(`SELECT cancer_type AS ct, ${ident(ep + '_time')} AS t, ${ident(ep)} AS e${gSel}${covCols}
           FROM s WHERE ${where} AND ${ident(ep)} IS NOT NULL AND ${ident(ep + '_time')} >= 0`),
    refreshFilterCounts(),
  ]);
  if (id !== runId) return;

  // time scale and administrative censoring
  const div = state.unit === 'years' ? 365.25 : 30.4375;
  const cap = parseFloat(state.censorAt);
  for (const r of rows) {
    r.t /= div;
    if (cap > 0 && r.t > cap) { r.t = cap; r.e = 0; }
  }
  const nTypes = new Set(rows.map((r) => r.ct)).size;
  $('cohort-n').textContent = cohortN.toLocaleString('en');
  $('cohort-n-sub').textContent = `patients · ${rows.length.toLocaleString('en')} with ${ENDPOINTS[ep].toLowerCase()} data · ${nTypes} cancer type${nTypes === 1 ? '' : 's'}`;

  const { labels, codes, notes } = makeGroups(rows);
  const keep = rows.map((_, i) => i).filter((i) => codes[i] >= 0);
  const A = keep.map((i) => rows[i]);
  const g = keep.map((i) => codes[i]);
  const time = A.map((r) => r.t), event = A.map((r) => r.e);
  const k = labels.length;
  const typeCode = new Map();
  const strata = state.stratify ? A.map((r) => { if (!typeCode.has(r.ct)) typeCode.set(r.ct, typeCode.size); return typeCode.get(r.ct); }) : null;

  const series = labels.map((label, j) => {
    const idx = g.map((v, i) => (v === j ? i : -1)).filter((i) => i >= 0);
    const fit = kmFit(idx.map((i) => time[i]), idx.map((i) => event[i]));
    return { label, color: COLORS[j], fit };
  });
  const totalEvents = event.reduce((s, e) => s + e, 0);

  // notes
  const byVar = state.groupBy ? dict[state.groupBy] : null;
  if (byVar?.group === 'therapy' || state.filters.some((f) => dict[f.col].group === 'therapy' && filterSql(f))) {
    notes.push('Therapy records are incomplete: “no record” does not mean untreated. Patients must also survive long enough to receive a treatment, which biases these comparisons (immortal-time bias).');
  }
  if (k > 1 && nTypes > 1 && !state.stratify) {
    notes.push(`This cohort spans ${nTypes} cancer types with very different baseline survival. Consider stratifying tests by cancer type, or filtering to one type.`);
  }
  if (A.length && totalEvents < 10) notes.push(`Only ${totalEvents} events in this cohort. Estimates will be unreliable.`);
  $('notes').replaceChildren(...notes.map((n) => h('div', { class: 'note' }, n)));

  // title, test, legend
  $('km-title').textContent = ENDPOINTS[ep] + (byVar ? ` by ${byVar.label.replace(/^[A-Z](?=[a-z])/, (c) => c.toLowerCase())}` : '');
  const test = $('km-test');
  test.replaceChildren();
  let lr = null;
  if (k > 1 && totalEvents > 0) {
    lr = logRank(time, event, g, k, strata);
    test.append(`${state.stratify ? 'Stratified log-rank' : 'Log-rank'} test: `, h('strong', {}, `p = ${fmtP(lr.p)}`),
      ` (χ² = ${fmtNum(lr.chisq)}, ${lr.df} df)`);
  }
  if (A.length) {
    test.append(`${lr ? ' · ' : ''}Median follow-up ${fmtNum(medianFollowUp(time, event))} ${state.unit}`);
  }
  $('km-legend').replaceChildren(...(k > 1 ? series.map((s) =>
    h('span', {}, h('i', { class: 'key', style: `border-top-color:${s.color}` }), s.label)) : []));

  // plot
  if (!A.length) {
    $('km-plot').replaceChildren(h('div', { class: 'empty' }, 'No patients match these filters with data for this endpoint.'));
    lastSvg = null;
  } else {
    lastSvg = renderKM($('km-plot'), $('km-tip'), series, {
      unit: state.unit, xmax: Math.max(...time), showCI: state.showCI, showCensor: state.showCensor,
      yLabel: ep === 'OS' || ep === 'DSS' ? 'Surviving' : 'Event-free',
    });
  }
  $('btn-svg').hidden = !lastSvg;

  // summary table
  let hr = null;
  if (k > 1 && totalEvents > 0) {
    const dummies = labels.slice(1).map((_, j) => g.map((v) => (v === j + 1 ? 1 : 0)));
    hr = coxFit(time, event, dummies, strata);
  }
  const yr = state.unit === 'years' ? 1 : 12;
  const head = ['Group', 'Patients', 'Events', `Median, ${state.unit} (95% CI)`, '1-year', '3-year', '5-year'];
  if (hr) head.push('Hazard ratio (95% CI)');
  const body = series.map((s, j) => {
    const m = kmMedian(s.fit);
    const med = Number.isNaN(m.median) ? 'Not reached' : `${fmtNum(m.median)} (${fmtNum(m.lower)} to ${Number.isNaN(m.upper) ? 'not reached' : fmtNum(m.upper)})`;
    const cells = [
      h('td', {}, k > 1 ? h('i', { class: 'key', style: `border-top-color:${s.color}` }) : null, s.label),
      h('td', {}, s.fit.n.toLocaleString('en')), h('td', {}, s.fit.events.toLocaleString('en')), h('td', {}, med),
      ...[1, 3, 5].map((y) => h('td', {}, pct(survAt(s.fit, y * yr)))),
    ];
    if (hr) {
      const ok = hr.converged && j > 0 && hr.se[j - 1] < 10;
      cells.push(h('td', {}, j === 0 ? 'Reference' : ok ? `${fmtHR(hr.hr[j - 1])} (${fmtHR(hr.lower[j - 1])} to ${fmtHR(hr.upper[j - 1])})` : 'Not estimable'));
    }
    return h('tr', {}, cells);
  });
  $('km-table').replaceChildren(h('thead', {}, h('tr', {}, head.map((t) => h('th', {}, t)))), h('tbody', {}, body));

  renderByType(A, g, labels, nTypes);
  renderCox(A, g, labels, strata);
}

// ---------- results for each cancer type ----------

let byType = { cols: [], rows: [], sort: null };

/** Benjamini–Hochberg adjusted p-values (NaN entries are ignored). */
function adjustBH(ps) {
  const idx = ps.map((p, i) => i).filter((i) => !Number.isNaN(ps[i])).sort((a, b) => ps[b] - ps[a]);
  const out = ps.map(() => NaN), m = idx.length;
  let min = 1;
  idx.forEach((i, rank) => { min = Math.min(min, (ps[i] * m) / (m - rank)); out[i] = min; });
  return out;
}

/** Repeats the current analysis inside each cancer type of the cohort. */
function renderByType(A, g, labels, nTypes) {
  const box = $('bytype');
  box.hidden = nTypes < 2;
  if (nTypes < 2) return;
  const k = labels.length, yr = state.unit === 'years' ? 1 : 12;
  const members = new Map();
  A.forEach((r, i) => { if (!members.has(r.ct)) members.set(r.ct, []); members.get(r.ct).push(i); });

  const rows = [...members].map(([ct, idx]) => {
    const time = idx.map((i) => A[i].t), event = idx.map((i) => A[i].e), grp = idx.map((i) => g[i]);
    const row = { ct, n: idx.length, events: event.reduce((s, e) => s + e, 0), p: NaN, hr: NaN, lo: NaN, hi: NaN };
    const fits = labels.map((_, j) => {
      const sel = grp.map((v, i) => (v === j ? i : -1)).filter((i) => i >= 0);
      return kmFit(sel.map((i) => time[i]), sel.map((i) => event[i]));
    });
    row.sizes = fits.map((f) => f.n);
    row.medians = fits.map((f) => (f.n ? kmMedian(f) : { median: NaN, lower: NaN, upper: NaN }));
    row.landmarks = [1, 3, 5].map((y) => survAt(fits[0], y * yr));
    // test only the groups that have patients in this cancer type
    const present = labels.map((_, j) => j).filter((j) => fits[j].n > 0);
    if (present.length > 1 && row.events > 0) {
      const code = new Map(present.map((j, i) => [j, i]));
      row.p = logRank(time, event, grp.map((v) => code.get(v)), present.length).p;
      if (k === 2) {
        const fit = coxFit(time, event, [grp.map((v) => (v === 1 ? 1 : 0))]);
        if (fit.converged && fit.se[0] < 10) { row.hr = fit.hr[0]; row.lo = fit.lower[0]; row.hi = fit.upper[0]; }
      }
    }
    return row;
  });
  const q = adjustBH(rows.map((r) => r.p));
  rows.forEach((r, i) => { r.q = q[i]; });

  const medianText = (m) => (Number.isNaN(m.median) ? 'not reached' : fmtNum(m.median));
  const cols = [
    { title: 'Cancer type', val: (r) => r.ct, text: (r) => r.ct },
    { title: 'Patients', val: (r) => r.n, text: (r) => r.n.toLocaleString('en') },
    { title: 'Events', val: (r) => r.events, text: (r) => r.events.toLocaleString('en') },
  ];
  if (k === 1) {
    cols.push({ title: `Median, ${state.unit} (95% CI)`, val: (r) => r.medians[0].median,
      text: (r) => { const m = r.medians[0]; return Number.isNaN(m.median) ? 'not reached' : `${fmtNum(m.median)} (${fmtNum(m.lower)} to ${Number.isNaN(m.upper) ? 'not reached' : fmtNum(m.upper)})`; } });
    ['1-year', '3-year', '5-year'].forEach((t, j) => cols.push({ title: t, val: (r) => r.landmarks[j], text: (r) => pct(r.landmarks[j]) }));
  } else {
    if (k <= 4) labels.forEach((l, j) => cols.push({ title: `${l}: patients; median ${state.unit}`, val: (r) => r.medians[j].median,
      text: (r) => (r.sizes[j] ? `${r.sizes[j].toLocaleString('en')}; ${medianText(r.medians[j])}` : '–') }));
    if (k === 2) cols.push({ title: 'Hazard ratio (95% CI)', val: (r) => r.hr,
      text: (r) => (Number.isNaN(r.hr) ? '–' : `${fmtHR(r.hr)} (${fmtHR(r.lo)} to ${fmtHR(r.hi)})`) });
    cols.push({ title: 'Log-rank p', val: (r) => r.p, text: (r) => fmtP(r.p), raw: true });
    cols.push({ title: 'Adjusted p', val: (r) => r.q, text: (r) => fmtP(r.q), raw: true });
  }
  const pCol = cols.findIndex((c) => c.title === 'Log-rank p');
  byType = { cols, rows, sort: pCol >= 0 ? { col: pCol, dir: 1 } : { col: 0, dir: 1 } };

  $('bytype-summary').textContent = `Results for each of the ${nTypes} cancer types`;
  $('bytype-hint').textContent = k === 1
    ? 'The same endpoint within each cancer type. Click a column heading to sort, or a cancer type to restrict the cohort to it.'
    : 'The same comparison repeated within each cancer type, using the groups defined above' +
      (state.groupBy && isNum(state.groupBy) && !state.withinType && state.split !== 'custom' ? ' (cut points come from the whole cohort; tick “Compute cut points within each cancer type” to split each type at its own values)' : '') +
      '. Adjusted p corrects for testing several cancer types (Benjamini–Hochberg); shaded rows have adjusted p below 0.05. ' +
      (k === 2 ? `Hazard ratios are ${labels[1]} relative to ${labels[0]}. ` : '') +
      'Click a column heading to sort, or a cancer type to restrict the cohort to it.';
  drawByType();
}

function drawByType() {
  const { cols, rows, sort } = byType;
  const c = cols[sort.col];
  const sorted = rows.slice().sort((a, b) => {
    const x = c.val(a), y = c.val(b);
    if (typeof x === 'string') return sort.dir * collator.compare(x, y);
    if (Number.isNaN(x) || Number.isNaN(y)) return Number.isNaN(x) - Number.isNaN(y); // missing last
    return sort.dir * (x - y);
  });
  $('bytype-table').replaceChildren(
    h('thead', {}, h('tr', {}, cols.map((col, i) => h('th', {
      'aria-sort': i === sort.col ? (sort.dir > 0 ? 'ascending' : 'descending') : null,
      onclick: () => { byType.sort = { col: i, dir: i === sort.col ? -sort.dir : 1 }; drawByType(); },
    }, col.title + (i === sort.col ? (sort.dir > 0 ? ' ▲' : ' ▼') : ''))))),
    h('tbody', {}, sorted.map((r) => h('tr', { class: r.q < 0.05 ? 'sig' : null }, cols.map((col, i) => h('td', {},
      i === 0 ? h('button', { type: 'button', class: 'linklike', title: `Restrict the cohort to ${r.ct}`, onclick: () => focusType(r.ct) }, r.ct) : col.text(r)))))));
}

async function focusType(ct) {
  let f = state.filters.find((x) => x.col === 'cancer_type');
  if (!f) { f = { col: 'cancer_type', values: [] }; state.filters.unshift(f); }
  f.values = [ct];
  await renderFilters();
  await update();
  window.scrollTo({ top: 0, behavior: 'smooth' });
}

function downloadByType() {
  const { cols, rows } = byType;
  const cell = (v) => (/[",\n]/.test(v) ? `"${v.replace(/"/g, '""')}"` : v);
  download('tcga_survivome_by_cancer_type.csv', 'text/csv',
    [cols.map((c) => cell(c.title)).join(','), ...rows.map((r) => cols.map((c) => (c.raw ? (Number.isNaN(c.val(r)) ? '' : String(c.val(r))) : cell(c.text(r)))).join(','))].join('\n'));
}

// ---------- Cox model ----------

function renderCox(A, g, labels, strata) {
  const out = $('cox-out');
  const msgs = [];
  const useGroups = state.coxGroups && labels.length > 1;
  const covs = state.covs.map((c, i) => ({ ...c, key: `c${i}` })).filter((c) => {
    if (useGroups && c.col === state.groupBy) { msgs.push(`${dict[c.col].label} is already in the model as the comparison groups.`); return false; }
    if (strata && c.col === 'cancer_type') { msgs.push('Cancer type is handled by stratification, so it is not also a covariate.'); return false; }
    return true;
  });
  if (!useGroups && !covs.length) {
    out.replaceChildren(h('p', { class: 'hint' }, 'Add a variable, or choose a comparison above, to fit a model.'));
    return;
  }

  // complete cases; drop categories too small to estimate
  let idx = A.map((_, i) => i).filter((i) => covs.every((c) => {
    const v = A[i][c.key];
    return v != null && !(isNum(c.col) && c.scale === 'log2' && v < 0);
  }));
  for (const c of covs.filter((c) => !isNum(c.col))) {
    const counts = new Map();
    for (const i of idx) counts.set(A[i][c.key], (counts.get(A[i][c.key]) || 0) + 1);
    idx = idx.filter((i) => counts.get(A[i][c.key]) >= MIN_LEVEL_N);
  }
  const dropped = A.length - idx.length;

  // design matrix
  const terms = []; // { label, ref } for reference rows, { label, col } for estimated ones
  const cols = [];
  const add = (label, values) => { terms.push({ label, at: cols.length }); cols.push(values); };
  if (useGroups) {
    const present = [...new Set(idx.map((i) => g[i]))].sort((a, b) => a - b);
    const name = state.groupBy ? dict[state.groupBy].label : 'Group';
    terms.push({ label: `${name}: ${labels[present[0]]}`, ref: true });
    for (const j of present.slice(1)) add(`${name}: ${labels[j]}`, idx.map((i) => (g[i] === j ? 1 : 0)));
  }
  for (const c of covs) {
    const d = dict[c.col];
    const vals = idx.map((i) => A[i][c.key]);
    if (isNum(c.col)) {
      if (c.scale === 'log2') add(`${d.label} (per doubling of value + 1)`, vals.map((v) => Math.log2(v + 1)));
      else if (c.scale === 'sd') {
        const mean = vals.reduce((s, v) => s + v, 0) / vals.length;
        const sd = Math.sqrt(vals.reduce((s, v) => s + (v - mean) ** 2, 0) / Math.max(1, vals.length - 1)) || 1;
        add(`${d.label} (per SD = ${fmtNum(sd)})`, vals.map((v) => (v - mean) / sd));
      } else add(`${d.label} (per ${d.unit ? d.unit.replace(/s$/, '') : 'unit'})`, vals);
    } else {
      const counts = new Map();
      for (const v of vals) counts.set(v, (counts.get(v) || 0) + 1);
      const levels = [...counts.keys()].sort((a, b) => counts.get(b) - counts.get(a));
      const ref = levels[0];
      terms.push({ label: `${d.label}: ${fmtLevel(c.col, ref)}`, ref: true });
      for (const l of levels.slice(1).sort((a, b) => collator.compare(fmtLevel(c.col, a), fmtLevel(c.col, b)))) {
        add(`${d.label}: ${fmtLevel(c.col, l)}`, vals.map((v) => (v === l ? 1 : 0)));
      }
    }
  }
  const time = idx.map((i) => A[i].t), event = idx.map((i) => A[i].e);
  const nEv = event.reduce((s, e) => s + e, 0);
  if (!cols.length || nEv < cols.length + 1) {
    out.replaceChildren(h('p', { class: 'hint' }, `Not enough data to fit this model (${idx.length} patients with complete data, ${nEv} events).`));
    return;
  }
  if (cols.length > 40) {
    out.replaceChildren(h('p', { class: 'hint' }, `This model has ${cols.length} terms. Remove a variable with many categories.`));
    return;
  }
  const fit = coxFit(time, event, cols, strata ? idx.map((i) => strata[i]) : null);

  // forest scale (log), shared by all rows
  const est = terms.filter((t) => !t.ref).map((t) => t.at).filter((a) => fit.se[a] < 10);
  const lo = Math.max(0.05, Math.min(0.5, ...est.map((a) => fit.lower[a])));
  const hi = Math.min(20, Math.max(2, ...est.map((a) => fit.upper[a])));
  const FW = 220, pad = 8;
  const fx = (v) => pad + ((Math.log(Math.min(hi, Math.max(lo, v))) - Math.log(lo)) / (Math.log(hi) - Math.log(lo))) * (FW - 2 * pad);
  const svgRow = (inner) => {
    const s = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
    s.setAttribute('viewBox', `0 0 ${FW} 26`); s.setAttribute('width', FW); s.setAttribute('height', 26);
    s.innerHTML = `<line x1="${fx(1)}" x2="${fx(1)}" y1="0" y2="26" style="stroke:var(--axis)"/>` + inner;
    return s;
  };
  const rows = terms.map((t) => {
    if (t.ref) return h('tr', { class: 'ref' }, h('td', {}, t.label), h('td', {}, '1 (reference)'), h('td', {}), h('td', {}), h('td', { class: 'forest' }, svgRow('')));
    const a = t.at, ok = fit.se[a] < 10 && Number.isFinite(fit.se[a]);
    if (!ok) return h('tr', {}, h('td', {}, t.label), h('td', { colspan: 3 }, 'Not estimable'), h('td', { class: 'forest' }, svgRow('')));
    const mark = `<line x1="${fx(fit.lower[a])}" x2="${fx(fit.upper[a])}" y1="13" y2="13" style="stroke:var(--s1)" stroke-width="2" stroke-linecap="round"/>` +
      `<circle cx="${fx(fit.hr[a])}" cy="13" r="4.5" style="fill:var(--s1);stroke:var(--surface)" stroke-width="2"/>`;
    const cell = h('td', { class: 'forest' }, svgRow(mark));
    cell.title = `Hazard ratio ${fmtHR(fit.hr[a])} (95% CI ${fmtHR(fit.lower[a])} to ${fmtHR(fit.upper[a])})`;
    return h('tr', {}, h('td', {}, t.label), h('td', {}, fmtHR(fit.hr[a])),
      h('td', {}, `${fmtHR(fit.lower[a])} to ${fmtHR(fit.upper[a])}`), h('td', {}, fmtP(fit.p[a])), cell);
  });
  const tickVals = [0.05, 0.1, 0.25, 0.5, 1, 2, 4, 10, 20].filter((v) => v >= lo && v <= hi);
  const axis = svgRow(tickVals.map((v) => `<text x="${fx(v)}" y="18" text-anchor="middle" style="fill:var(--muted);font-size:11px">${v}</text>`).join(''));
  rows.push(h('tr', {}, h('td', { colspan: 4, style: 'color:var(--muted);font-size:12px' }, 'Hazard ratio, log scale. Right of 1 means higher risk.'), h('td', { class: 'forest' }, axis)));

  const foot = [`${fit.n.toLocaleString('en')} patients, ${fit.nEvent.toLocaleString('en')} events.`];
  if (dropped) foot.push(`${dropped.toLocaleString('en')} left out (missing values, or categories with fewer than ${MIN_LEVEL_N} patients).`);
  foot.push(`Likelihood-ratio test: χ² = ${fmtNum(fit.lrChisq)}, ${fit.lrDf} df, p = ${fmtP(fit.lrP)}.`);
  if (strata) foot.push('Stratified by cancer type.');
  if (!fit.converged) msgs.push('The model did not converge, usually because a category has no events. Remove a variable or widen the cohort.');
  out.replaceChildren(
    h('div', { class: 'table-scroll' }, h('table', { class: 'data' },
      h('thead', {}, h('tr', {}, ['Term', 'Hazard ratio', '95% CI', 'p', ''].map((t) => h('th', {}, t)))),
      h('tbody', {}, rows))),
    h('p', { class: 'hint' }, foot.join(' ')),
    ...msgs.map((m) => h('div', { class: 'note' }, m)));
}

// ---------- downloads ----------

function download(name, mime, content) {
  const a = h('a', { href: URL.createObjectURL(new Blob([content], { type: mime })), download: name });
  document.body.append(a); a.click(); a.remove();
  setTimeout(() => URL.revokeObjectURL(a.href), 1000);
}
async function downloadCsv() {
  const rows = await query(`SELECT * FROM s WHERE ${whereSql()} ORDER BY cancer_type, patient_id`);
  if (!rows.length) return;
  const keys = Object.keys(rows[0]);
  const cell = (v) => (v == null ? '' : /[",\n]/.test(String(v)) ? `"${String(v).replace(/"/g, '""')}"` : String(v));
  download('tcga_survivome_cohort.csv', 'text/csv', [keys.join(','), ...rows.map((r) => keys.map((k) => cell(r[k])).join(','))].join('\n'));
}
function downloadSvg() {
  if (!lastSvg) return;
  const css = getComputedStyle(document.documentElement);
  const clone = lastSvg.cloneNode(true);
  clone.setAttribute('xmlns', 'http://www.w3.org/2000/svg');
  const style = document.createElementNS('http://www.w3.org/2000/svg', 'style');
  style.textContent = 'text{font:12px system-ui,-apple-system,"Segoe UI",sans-serif;fill:var(--ink-2)}';
  clone.prepend(style);
  // resolve theme variables so the file stands alone
  const src = new XMLSerializer().serializeToString(clone)
    .replace(/var\((--[\w-]+)\)/g, (_, v) => css.getPropertyValue(v).trim() || '#000');
  download('tcga_survivome_km.svg', 'image/svg+xml', src);
}

// ---------- start ----------

async function main() {
  const [d] = await Promise.all([
    fetch('data/data_dictionary.json').then((r) => { if (!r.ok) throw new Error('Could not load the data dictionary'); return r.json(); }),
    initDb('data/survivome.parquet'),
  ]);
  for (const c of d.columns) dict[c.column] = c;
  variables = d.columns.filter((c) => ['categorical', 'numeric', 'boolean'].includes(c.type));
  loadState();

  fillVarSelect($('add-filter'), '+ Add a filter…');
  fillVarSelect($('group-by'), 'No comparison (one curve)');
  fillVarSelect($('add-cov'), '+ Add a variable…');
  syncControls();
  renderCovs();
  await renderFilters();

  bind('exclude-qc', 'excludeQC', { checkbox: true });
  bind('endpoint', 'endpoint'); bind('group-by', 'groupBy'); bind('split', 'split'); bind('cutoff', 'cutoff');
  bind('within-type', 'withinType', { checkbox: true });
  bind('unit', 'unit'); bind('censor-at', 'censorAt');
  bind('stratify', 'stratify', { checkbox: true });
  bind('show-ci', 'showCI', { checkbox: true }); bind('show-censor', 'showCensor', { checkbox: true });
  bind('cox-groups', 'coxGroups', { checkbox: true });
  $('unit').addEventListener('change', () => { state.censorAt = ''; syncControls(); }, { capture: true });

  $('add-filter').addEventListener('change', async (ev) => {
    const col = ev.target.value; ev.target.value = '';
    if (!col) return;
    state.filters.push(isNum(col) ? { col, min: '', max: '' } : { col, values: [] });
    await renderFilters(); update();
  });
  $('add-cov').addEventListener('change', (ev) => {
    const col = ev.target.value; ev.target.value = '';
    if (!col || state.covs.some((c) => c.col === col)) return;
    state.covs.push({ col, scale: 'unit' });
    renderCovs(); update();
  });
  $('clear-filters').addEventListener('click', async () => {
    state.filters = defaults().filters;
    await renderFilters(); update();
  });
  $('btn-csv').addEventListener('click', downloadCsv);
  $('btn-svg').addEventListener('click', downloadSvg);
  $('btn-bytype-csv').addEventListener('click', downloadByType);
  $('btn-link').addEventListener('click', async (ev) => {
    try { await navigator.clipboard.writeText(location.href); ev.target.textContent = 'Link copied'; }
    catch { ev.target.textContent = 'Copy the address bar'; }
    setTimeout(() => { ev.target.textContent = 'Copy link to this view'; }, 1800);
  });

  $('loading').hidden = true;
  $('app').hidden = false;
  await update();
}

main().catch((err) => {
  console.error(err);
  const l = $('loading');
  l.hidden = false; l.classList.add('error');
  l.textContent = `Something went wrong while loading: ${err.message}`;
});
