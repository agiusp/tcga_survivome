// Kaplan–Meier plot as inline SVG: step curves, confidence bands, censoring marks,
// a numbers-at-risk table aligned to the time axis, and a crosshair tooltip.
import { survAt, nRiskAt } from './stats.js';

const NS = 'http://www.w3.org/2000/svg';
function el(tag, attrs = {}, parent) {
  const e = document.createElementNS(NS, tag);
  for (const k in attrs) e.setAttribute(k, attrs[k]);
  if (parent) parent.appendChild(e);
  return e;
}
function text(str, attrs, parent) {
  const t = el('text', attrs, parent);
  t.textContent = str;
  return t;
}
const short = (s, n) => (s.length > n ? s.slice(0, n - 1) + '…' : s);

function tickStep(xmax, unit) {
  const steps = unit === 'years' ? [0.25, 0.5, 1, 2, 5, 10, 20] : [1, 2, 3, 6, 12, 24, 36, 48, 60, 120, 240];
  return steps.find((s) => xmax / s <= 10) ?? steps[steps.length - 1];
}

/**
 * @param container element that receives the <svg>
 * @param tip tooltip element (absolutely positioned inside the same wrapper)
 * @param series [{ label, color, fit }]
 * @param opts { unit, xmax, showCI, showCensor, yLabel }
 */
export function renderKM(container, tip, series, opts) {
  const W = 860, mL = 158, mR = 18, mT = 10, plotH = 330;
  const riskTop = mT + plotH + 50, rowH = 17;
  const H = riskTop + 8 + series.length * rowH + 8;
  const pw = W - mL - mR;
  const xmax = opts.xmax > 0 ? opts.xmax : 1;
  const x = (t) => mL + (Math.min(t, xmax) / xmax) * pw;
  const y = (s) => mT + (1 - s) * plotH;

  container.replaceChildren();
  const svg = el('svg', { viewBox: `0 0 ${W} ${H}`, role: 'img', 'aria-label': `Kaplan–Meier curves: ${opts.yLabel}` }, container);

  // grid and axes
  for (const p of [0, 0.25, 0.5, 0.75, 1]) {
    el('line', { x1: mL, x2: W - mR, y1: y(p), y2: y(p), style: `stroke:var(${p === 0 ? '--axis' : '--grid'})`, 'stroke-width': 1 }, svg);
    text(`${p * 100}%`, { x: mL - 8, y: y(p) + 4, 'text-anchor': 'end' }, svg);
  }
  const step = tickStep(xmax, opts.unit);
  const ticks = [];
  for (let t = 0; t <= xmax + 1e-9; t += step) ticks.push(t);
  for (const t of ticks) {
    el('line', { x1: x(t), x2: x(t), y1: y(0), y2: y(0) + 5, style: 'stroke:var(--axis)' }, svg);
    text(String(+t.toFixed(2)), { x: x(t), y: y(0) + 19, 'text-anchor': 'middle' }, svg);
  }
  text(`Time since diagnosis (${opts.unit})`, { x: mL + pw / 2, y: y(0) + 38, 'text-anchor': 'middle', class: 'axis-title' }, svg);
  const yl = text(opts.yLabel, { x: 0, y: 0, 'text-anchor': 'middle', class: 'axis-title', transform: `translate(${mL - 46} ${mT + plotH / 2}) rotate(-90)` }, svg);
  yl.setAttribute('aria-hidden', 'true');

  // curves
  for (const s of series) {
    const f = s.fit;
    if (!f.time.length) continue;
    const last = f.time[f.time.length - 1];
    if (opts.showCI) {
      let up = `M${x(0)} ${y(1)}`, lo = '', prevLo = 1, end = 0;
      for (let i = 0; i < f.time.length; i++) {
        if (!f.nEvent[i]) continue;
        if (Number.isNaN(f.upper[i])) break;
        up += `H${x(f.time[i])}V${y(f.upper[i])}`;
        lo = `V${y(prevLo)}H${x(f.time[i])}` + lo; // traced backwards
        prevLo = f.lower[i]; end = i;
      }
      const stop = Number.isNaN(f.upper[f.time.length - 1]) ? f.time[end] : last;
      el('path', { d: `${up}H${x(stop)}V${y(prevLo)}${lo}V${y(1)}Z`, style: `fill:${s.color}`, 'fill-opacity': 0.13, stroke: 'none' }, svg);
    }
    let d = `M${x(0)} ${y(1)}`;
    for (let i = 0; i < f.time.length; i++) if (f.nEvent[i]) d += `H${x(f.time[i])}V${y(f.surv[i])}`;
    d += `H${x(last)}`;
    el('path', { d, fill: 'none', style: `stroke:${s.color}`, 'stroke-width': 2, 'stroke-linejoin': 'round' }, svg);
    if (opts.showCensor) {
      let c = '';
      for (let i = 0; i < f.time.length; i++) {
        if (f.nCensor[i]) c += `M${x(f.time[i])} ${y(f.surv[i]) - 3.5}v7`;
      }
      if (c) el('path', { d: c, fill: 'none', style: `stroke:${s.color}`, 'stroke-width': 1, 'stroke-opacity': 0.75 }, svg);
    }
  }

  // numbers at risk
  text('Number at risk', { x: 4, y: riskTop, class: 'axis-title' }, svg);
  series.forEach((s, k) => {
    const yy = riskTop + 8 + (k + 1) * rowH - 4;
    el('line', { x1: 4, x2: 18, y1: yy - 4, y2: yy - 4, style: `stroke:${s.color}`, 'stroke-width': 3, 'stroke-linecap': 'round' }, svg);
    const lab = text(short(s.label, 17), { x: 24, y: yy, class: 'risk-label' }, svg);
    el('title', {}, lab).textContent = s.label;
    for (const t of ticks) text(String(nRiskAt(s.fit, t)), { x: x(t), y: yy, 'text-anchor': 'middle' }, svg);
  });

  // hover layer
  const cross = el('line', { y1: mT, y2: y(0), style: 'stroke:var(--muted)', 'stroke-width': 1, visibility: 'hidden' }, svg);
  const hit = el('rect', { x: mL, y: mT, width: pw, height: plotH, fill: 'transparent' }, svg);
  const hide = () => { cross.setAttribute('visibility', 'hidden'); tip.hidden = true; };
  hit.addEventListener('pointerleave', hide);
  hit.addEventListener('pointermove', (ev) => {
    const box = svg.getBoundingClientRect(), scale = box.width / W;
    const px = (ev.clientX - box.left) / scale;
    const t = Math.max(0, Math.min(xmax, ((px - mL) / pw) * xmax));
    cross.setAttribute('x1', x(t)); cross.setAttribute('x2', x(t));
    cross.setAttribute('visibility', 'visible');
    tip.replaceChildren();
    const head = document.createElement('div');
    head.className = 'tip-t';
    head.textContent = `At ${t.toFixed(opts.unit === 'years' ? 2 : 1)} ${opts.unit}`;
    tip.appendChild(head);
    for (const s of series) {
      const row = document.createElement('div'); row.className = 'row';
      const key = document.createElement('span'); key.className = 'key'; key.style.borderTopColor = s.color;
      const v = document.createElement('b');
      const sv = survAt(s.fit, t);
      v.textContent = Number.isNaN(sv) ? '–' : `${(sv * 100).toFixed(1)}%`;
      const nm = document.createElement('span'); nm.className = 'nm';
      nm.textContent = `${short(s.label, 28)} · ${nRiskAt(s.fit, t)} at risk`;
      row.append(key, v, nm); tip.appendChild(row);
    }
    tip.hidden = false;
    const wrap = container.parentElement.getBoundingClientRect();
    const left = ev.clientX - wrap.left;
    const flip = left + tip.offsetWidth + 24 > wrap.width;
    tip.style.left = `${flip ? Math.max(0, left - tip.offsetWidth - 14) : left + 14}px`;
    tip.style.top = `${Math.max(0, ev.clientY - wrap.top - tip.offsetHeight / 2)}px`;
  });
  return svg;
}
