// Survival statistics in plain JS: Kaplan–Meier, log-rank, Cox (Efron ties).
// No dependencies; runs in the browser and in Node (see tests/).
// Results are checked against R's `survival` package in tests/stats.test.mjs.

const Z95 = 1.959963984540054;

// ---------- distributions ----------

function gammaLn(x) {
  const c = [76.18009172947146, -86.50532032941677, 24.01409824083091,
    -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5];
  let y = x, tmp = x + 5.5;
  tmp -= (x + 0.5) * Math.log(tmp);
  let ser = 1.000000000190015;
  for (let j = 0; j < 6; j++) ser += c[j] / ++y;
  return -tmp + Math.log(2.5066282746310005 * ser / x);
}

/** Regularized upper incomplete gamma Q(a, x). */
export function gammaQ(a, x) {
  if (!(x > 0)) return 1;
  if (x < a + 1) {
    let ap = a, sum = 1 / a, del = sum;
    for (let n = 0; n < 1000; n++) {
      ap += 1; del *= x / ap; sum += del;
      if (Math.abs(del) < Math.abs(sum) * 1e-15) break;
    }
    return 1 - sum * Math.exp(-x + a * Math.log(x) - gammaLn(a));
  }
  const tiny = 1e-300;
  let b = x + 1 - a, c = 1 / tiny, d = 1 / b, h = d;
  for (let i = 1; i < 1000; i++) {
    const an = -i * (i - a);
    b += 2;
    d = an * d + b; if (Math.abs(d) < tiny) d = tiny;
    c = b + an / c; if (Math.abs(c) < tiny) c = tiny;
    d = 1 / d;
    const del = d * c; h *= del;
    if (Math.abs(del - 1) < 1e-15) break;
  }
  return Math.exp(-x + a * Math.log(x) - gammaLn(a)) * h;
}

/** Upper-tail chi-square p-value. */
export const pChisq = (x, df) => (x >= 0 && df > 0 ? gammaQ(df / 2, x / 2) : NaN);
/** Two-sided normal p-value. */
export const pNorm2 = (z) => gammaQ(0.5, (z * z) / 2);

// ---------- linear algebra ----------

/** Inverse of a symmetric matrix (array of rows); null if singular. */
export function invert(A) {
  const n = A.length;
  const M = A.map((r, i) => [...r, ...Array.from({ length: n }, (_, j) => (i === j ? 1 : 0))]);
  const scale = Math.max(...A.map((r, i) => Math.abs(r[i])), 1e-300);
  for (let c = 0; c < n; c++) {
    let piv = c;
    for (let r = c + 1; r < n; r++) if (Math.abs(M[r][c]) > Math.abs(M[piv][c])) piv = r;
    if (Math.abs(M[piv][c]) < 1e-11 * scale) return null;
    [M[c], M[piv]] = [M[piv], M[c]];
    const d = M[c][c];
    for (let j = 0; j < 2 * n; j++) M[c][j] /= d;
    for (let r = 0; r < n; r++) {
      if (r === c) continue;
      const f = M[r][c];
      if (f !== 0) for (let j = 0; j < 2 * n; j++) M[r][j] -= f * M[c][j];
    }
  }
  return M.map((r) => r.slice(n));
}

// ---------- Kaplan–Meier ----------

/**
 * Kaplan–Meier estimate with Greenwood variance and log-transformed 95% CI
 * (R's survfit default, conf.type = "log").
 * @param {ArrayLike<number>} time
 * @param {ArrayLike<number>} event 1 = event, 0 = censored
 */
export function kmFit(time, event) {
  const n = time.length;
  const idx = Array.from({ length: n }, (_, i) => i).sort((a, b) => time[a] - time[b]);
  const fit = {
    n, events: 0, time: [], nRisk: [], nEvent: [], nCensor: [], surv: [], lower: [], upper: [],
    sortedTimes: idx.map((i) => time[i]),
  };
  let s = 1, v = 0, i = 0;
  while (i < n) {
    const t = time[idx[i]];
    let d = 0, c = 0, j = i;
    while (j < n && time[idx[j]] === t) { if (event[idx[j]]) d++; else c++; j++; }
    const r = n - i;
    if (d > 0) {
      s *= (r - d) / r;
      v += r > d ? d / (r * (r - d)) : Infinity;
      fit.events += d;
    }
    const ok = s > 0 && Number.isFinite(v);
    fit.time.push(t); fit.nRisk.push(r); fit.nEvent.push(d); fit.nCensor.push(c);
    fit.surv.push(s);
    fit.lower.push(ok ? Math.exp(Math.log(s) - Z95 * Math.sqrt(v)) : NaN);
    fit.upper.push(ok ? Math.min(1, Math.exp(Math.log(s) + Z95 * Math.sqrt(v))) : NaN);
    i = j;
  }
  return fit;
}

function lastIndexAtOrBefore(arr, t) {
  let lo = 0, hi = arr.length; // first index with arr[i] > t
  while (lo < hi) { const m = (lo + hi) >> 1; if (arr[m] <= t) lo = m + 1; else hi = m; }
  return lo - 1;
}

/** Step-function value of a KM curve ('surv' | 'lower' | 'upper') at time t. */
export function survAt(fit, t, key = 'surv') {
  const i = lastIndexAtOrBefore(fit.time, t);
  if (i < 0) return 1;
  // Past the last observation the curve is undefined unless it has reached 0.
  if (t > fit.time[fit.time.length - 1] && fit.surv[fit.surv.length - 1] > 0) return NaN;
  return fit[key][i];
}

/** Number still at risk just before time t (time >= t). */
export function nRiskAt(fit, t) {
  const a = fit.sortedTimes;
  let lo = 0, hi = a.length; // first index with a[i] >= t
  while (lo < hi) { const m = (lo + hi) >> 1; if (a[m] < t) lo = m + 1; else hi = m; }
  return a.length - lo;
}

function quantileOf(fit, key, q) {
  const target = 1 - q, tol = 1e-9;
  const y = fit[key];
  for (let i = 0; i < y.length; i++) {
    const v = Number.isNaN(y[i]) ? 0 : y[i];
    if (v <= target + tol) {
      if (Math.abs(v - target) <= tol) {
        // Curve sits exactly on the target: R averages with the next drop.
        for (let j = i + 1; j < y.length; j++) {
          const w = Number.isNaN(y[j]) ? 0 : y[j];
          if (w < target - tol) return (fit.time[i] + fit.time[j]) / 2;
        }
      }
      return fit.time[i];
    }
  }
  return NaN;
}

/** Median survival with 95% CI (NaN where not reached). */
export function kmMedian(fit) {
  return {
    median: quantileOf(fit, 'surv', 0.5),
    lower: quantileOf(fit, 'lower', 0.5),
    upper: quantileOf(fit, 'upper', 0.5),
  };
}

/** Median follow-up by the reverse Kaplan–Meier method. */
export function medianFollowUp(time, event) {
  const rev = Array.from(event, (e) => (e ? 0 : 1));
  return kmMedian(kmFit(time, rev)).median;
}

// ---------- log-rank ----------

function splitByStratum(n, strata) {
  if (!strata) return [Array.from({ length: n }, (_, i) => i)];
  const m = new Map();
  for (let i = 0; i < n; i++) {
    if (!m.has(strata[i])) m.set(strata[i], []);
    m.get(strata[i]).push(i);
  }
  return [...m.values()];
}

/**
 * k-sample log-rank test (optionally stratified), as in R's survdiff(rho = 0).
 * @param group integer codes 0..k-1
 */
export function logRank(time, event, group, k, strata = null) {
  const O = new Array(k).fill(0), E = new Array(k).fill(0);
  const V = Array.from({ length: k }, () => new Array(k).fill(0));
  for (const members of splitByStratum(time.length, strata)) {
    members.sort((a, b) => time[a] - time[b]);
    const atRisk = new Array(k).fill(0);
    for (const i of members) atRisk[group[i]]++;
    let nTot = members.length, i = 0;
    while (i < members.length) {
      const t = time[members[i]];
      const d = new Array(k).fill(0), gone = new Array(k).fill(0);
      let dTot = 0, j = i;
      while (j < members.length && time[members[j]] === t) {
        const g = group[members[j]];
        gone[g]++;
        if (event[members[j]]) { d[g]++; dTot++; }
        j++;
      }
      if (dTot > 0) {
        const f = nTot > 1 ? (dTot * (nTot - dTot)) / (nTot - 1) : 0;
        for (let a = 0; a < k; a++) {
          O[a] += d[a];
          E[a] += (dTot * atRisk[a]) / nTot;
          const pa = atRisk[a] / nTot;
          for (let b = 0; b < k; b++) V[a][b] += f * pa * ((a === b ? 1 : 0) - atRisk[b] / nTot);
        }
      }
      for (let a = 0; a < k; a++) atRisk[a] -= gone[a];
      nTot -= j - i;
      i = j;
    }
  }
  const df = k - 1;
  let chisq = NaN;
  if (df >= 1) {
    const Vi = invert(V.slice(0, df).map((r) => r.slice(0, df)));
    if (Vi) {
      chisq = 0;
      for (let a = 0; a < df; a++) for (let b = 0; b < df; b++) chisq += (O[a] - E[a]) * Vi[a][b] * (O[b] - E[b]);
    }
  }
  return { observed: O, expected: E, chisq, df, p: pChisq(chisq, df) };
}

// ---------- Cox proportional hazards ----------

/**
 * Cox model by Newton–Raphson with Efron's handling of ties (R's coxph default).
 * @param cols array of covariate columns, each of length n
 * @param strata optional stratum code per subject
 */
export function coxFit(time, event, cols, strata = null, { maxIter = 40, eps = 1e-9 } = {}) {
  const n = time.length, p = cols.length;
  const means = cols.map((c) => { let s = 0; for (let i = 0; i < n; i++) s += c[i]; return s / n; });
  const X = cols.map((c, a) => Float64Array.from(c, (v) => v - means[a]));
  const order = Array.from({ length: n }, (_, i) => i);
  const sKey = strata || new Array(n).fill(0);
  order.sort((a, b) => (sKey[a] < sKey[b] ? -1 : sKey[a] > sKey[b] ? 1 : time[b] - time[a]));
  let nEvent = 0;
  for (let i = 0; i < n; i++) if (event[i]) nEvent++;

  function evaluate(beta) {
    let l = 0;
    const g = new Float64Array(p), H = new Float64Array(p * p);
    let S0 = 0;
    const S1 = new Float64Array(p), S2 = new Float64Array(p * p);
    const D1 = new Float64Array(p), D2 = new Float64Array(p * p), x = new Float64Array(p), m = new Float64Array(p);
    let i = 0, cur;
    while (i < n) {
      const st = sKey[order[i]];
      if (st !== cur) { cur = st; S0 = 0; S1.fill(0); S2.fill(0); }
      const t = time[order[i]];
      let D0 = 0, d = 0;
      D1.fill(0); D2.fill(0);
      while (i < n && sKey[order[i]] === st && time[order[i]] === t) {
        const s = order[i];
        let eta = 0;
        for (let a = 0; a < p; a++) { x[a] = X[a][s]; eta += beta[a] * x[a]; }
        const w = Math.exp(eta), ev = event[s] ? 1 : 0;
        S0 += w;
        if (ev) { D0 += w; d++; l += eta; }
        for (let a = 0; a < p; a++) {
          const wx = w * x[a];
          S1[a] += wx;
          if (ev) { D1[a] += wx; g[a] += x[a]; }
          for (let b = a; b < p; b++) {
            const wxx = wx * x[b];
            S2[a * p + b] += wxx;
            if (ev) D2[a * p + b] += wxx;
          }
        }
        i++;
      }
      for (let e = 0; e < d; e++) {
        const f = e / d, s0 = S0 - f * D0;
        l -= Math.log(s0);
        for (let a = 0; a < p; a++) { m[a] = (S1[a] - f * D1[a]) / s0; g[a] -= m[a]; }
        for (let a = 0; a < p; a++) for (let b = a; b < p; b++) {
          H[a * p + b] += (S2[a * p + b] - f * D2[a * p + b]) / s0 - m[a] * m[b];
        }
      }
    }
    const Hm = Array.from({ length: p }, (_, a) => Array.from({ length: p }, (_, b) => (a <= b ? H[a * p + b] : H[b * p + a])));
    return { l, g, H: Hm };
  }

  let beta = new Float64Array(p);
  let cur = evaluate(beta);
  const loglik0 = cur.l;
  let converged = p === 0, iterations = 0, singular = false;
  while (!converged && iterations < maxIter) {
    iterations++;
    const Hi = invert(cur.H);
    if (!Hi) { singular = true; break; }
    const step = Hi.map((row) => row.reduce((s, v, b) => s + v * cur.g[b], 0));
    let next, nb, halvings = 0;
    for (;;) {
      nb = beta.map((v, a) => v + step[a]);
      next = evaluate(nb);
      if (Number.isFinite(next.l) && next.l >= cur.l - 1e-12) break;
      if (++halvings > 25) break;
      for (let a = 0; a < p; a++) step[a] /= 2;
    }
    const change = Math.abs(next.l - cur.l);
    beta = nb; cur = next;
    if (change <= eps * (Math.abs(cur.l) + eps)) converged = true;
  }
  const Vm = singular ? null : invert(cur.H);
  const se = Array.from({ length: p }, (_, a) => (Vm ? Math.sqrt(Vm[a][a]) : NaN));
  const coef = Array.from(beta);
  const lr = 2 * (cur.l - loglik0);
  return {
    n, nEvent, coef, se,
    z: coef.map((b, a) => b / se[a]),
    p: coef.map((b, a) => pNorm2(b / se[a])),
    hr: coef.map(Math.exp),
    lower: coef.map((b, a) => Math.exp(b - Z95 * se[a])),
    upper: coef.map((b, a) => Math.exp(b + Z95 * se[a])),
    loglik0, loglik: cur.l, lrChisq: lr, lrDf: p, lrP: pChisq(lr, p),
    converged: converged && !singular, singular, iterations,
  };
}

// ---------- helpers for grouping ----------

/** Quantile, type 7 (R default), of a sorted numeric array. */
export function quantileSorted(sorted, q) {
  if (!sorted.length) return NaN;
  const h = (sorted.length - 1) * q, lo = Math.floor(h), hi = Math.ceil(h);
  return sorted[lo] + (sorted[hi] - sorted[lo]) * (h - lo);
}
