// Checks site/js/stats.js against R's survival package (fixtures from make_fixtures.R).
// Run: node --test tests/stats.test.mjs
import test from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { kmFit, kmMedian, logRank, coxFit, pChisq } from '../site/js/stats.js';

const fx = JSON.parse(readFileSync(new URL('./fixtures.json', import.meta.url)));
const { time, event, g, age, male, stratum } = fx.data;
const close = (a, b, tol = 1e-8, msg = '') => {
  if (a == null || Number.isNaN(a)) { assert.ok(b == null || Number.isNaN(b), `${msg}: expected ${a}, got ${b}`); return; }
  assert.ok(Math.abs(a - b) <= tol * Math.max(1, Math.abs(a)), `${msg}: expected ${a}, got ${b}`);
};
const dummies = [1, 2, 3].map((k) => g.map((v) => (v === k ? 1 : 0)));

test('chi-square tail', () => {
  fx.pchisq.x.forEach((x, i) => close(fx.pchisq.p[i], pChisq(x, fx.pchisq.df[i]), 1e-9, `pchisq ${x}`));
});

test('Kaplan–Meier matches survfit', () => {
  fx.km.forEach((ref, k) => {
    const idx = g.map((v, i) => (v === k ? i : -1)).filter((i) => i >= 0);
    const fit = kmFit(idx.map((i) => time[i]), idx.map((i) => event[i]));
    assert.equal(fit.time.length, ref.time.length);
    ref.time.forEach((t, i) => {
      close(t, fit.time[i], 1e-12, 'time');
      close(ref.surv[i], fit.surv[i], 1e-10, `surv g${k}`);
      close(ref.lower[i], fit.lower[i], 1e-9, `lower g${k}`);
      close(ref.upper[i], fit.upper[i], 1e-9, `upper g${k}`);
      assert.equal(fit.nRisk[i], ref.n_risk[i]);
    });
    const m = kmMedian(fit);
    close(ref.median, m.median, 1e-10, `median g${k}`);
    close(ref.median_lower, m.lower, 1e-10, `median lower g${k}`);
    close(ref.median_upper, m.upper, 1e-10, `median upper g${k}`);
  });
});

test('log-rank matches survdiff', () => {
  const r = logRank(time, event, g, 4);
  close(fx.logrank.chisq, r.chisq, 1e-8, 'chisq');
  fx.logrank.obs.forEach((o, i) => close(o, r.observed[i], 1e-10, 'obs'));
  fx.logrank.exp.forEach((e, i) => close(e, r.expected[i], 1e-8, 'exp'));
  close(fx.logrank_strat.chisq, logRank(time, event, g, 4, stratum).chisq, 1e-8, 'stratified chisq');
});

test('Cox matches coxph (Efron)', () => {
  for (const [ref, strata] of [[fx.cox, null], [fx.cox_strat, stratum]]) {
    const r = coxFit(time, event, [age, male, ...dummies], strata);
    assert.ok(r.converged);
    ref.coef.forEach((b, i) => close(b, r.coef[i], 1e-6, `coef ${i}`));
    ref.se.forEach((s, i) => close(s, r.se[i], 1e-6, `se ${i}`));
    close(ref.loglik[0], r.loglik0, 1e-9, 'loglik0');
    close(ref.loglik[1], r.loglik, 1e-9, 'loglik');
    close(ref.lr, r.lrChisq, 1e-7, 'LR');
  }
});
