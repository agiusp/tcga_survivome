# How the survival statistics are computed

Everything on the site is computed in your browser by `site/js/stats.js`, a small
dependency-free module. The methods follow the defaults of R's `survival` package
(`survfit`, `survdiff`, `coxph`), and the results are tested against it
(see [Validation](#validation)).

## 1. Which patients enter an analysis

1. **Cohort.** Patients matching every filter in the Cohort panel. Within one categorical
   filter the ticked values are alternatives (OR); different filters combine (AND). With
   "Exclude cases flagged by TCGA quality review" on, patients with the AWG do-not-use flag,
   the pathology-review exclusion or a redaction are removed.
2. **Endpoint data.** Patients need both an event indicator and a non-negative time for the
   chosen endpoint. The cohort count shows how many have it.
3. **Comparison groups.** Patients with no value for the "Compare by" variable are left out,
   and the page says how many.

### Endpoints

Endpoints come from the TCGA Clinical Data Resource (Liu et al., *Cell* 2018), which
defines them in full.

| Endpoint | Event |
|---|---|
| Overall survival (OS) | Death from any cause |
| Disease-specific survival (DSS) | Death from the cancer |
| Progression-free interval (PFI) | Progression, recurrence, new primary tumor, or death with tumor |
| Disease-free interval (DFI) | Recurrence or new tumor in patients who were disease-free after initial treatment |

Patients without the event are censored at last follow-up. Competing events (for example a
death from another cause under DSS) are treated as censoring, as coded in TCGA-CDR. No
competing-risks method is applied.

### Time

Times are stored in days since diagnosis and shown as months (days / 30.4375) or years
(days / 365.25). If you set "Censor follow-up at", any patient followed beyond that time is
censored at it. This changes the curves, the tests and the Cox model, not just the plot range.

## 2. Comparison groups

**Categorical variables** use their categories as groups. The plot holds at most 8 groups.
If there are more, the 8 with the most patients are kept and the rest are left out, with a
note.

**Numeric variables** are cut into groups:

| Split | Groups |
|---|---|
| Median | at or below the median, above it |
| Tertiles, quartiles | 3 or 4 groups at the 1/3 and 2/3, or 25/50/75%, quantiles |
| Top vs. bottom quartile | at or below the 25% quantile vs. above the 75% quantile; the middle half is left out |
| Custom cut point | at or below your value, above it |

- A value equal to a cut point goes into the lower group.
- Quantiles use linear interpolation (R's default, type 7).
- Cut points are computed from the patients in the current cohort who have endpoint data.
- With "Compute cut points within each cancer type", each cancer type gets its own cut
  points, so "high" means high relative to that patient's cancer type. This matters for
  measures such as TMB whose typical values differ a lot between cancer types.
- If many patients share one value (for example a median of 0), a group can be empty. Empty
  groups are dropped, with a note.

## 3. Kaplan–Meier curves

For each group, at each distinct event time $t_i$ with $d_i$ events among $n_i$ patients
still at risk:

$$\hat S(t) = \prod_{t_i \le t} \left(1 - \frac{d_i}{n_i}\right)$$

Patients censored at $t_i$ count as at risk at $t_i$.

- **Confidence bands** are pointwise 95% intervals from Greenwood's variance on the log
  scale (R's `conf.type = "log"`):
  $\hat S(t)\,\exp\!\left(\pm 1.96 \sqrt{\sum_{t_i \le t} \frac{d_i}{n_i (n_i - d_i)}}\right)$,
  with the upper limit capped at 100%. They are not defined once the curve reaches 0.
- **Censoring marks** are the small vertical ticks, drawn at each time a patient was censored.
- **Number at risk** at time $t$ counts patients whose follow-up time is at least $t$.
- **Median survival** is the first time the curve is at or below 50%. If the curve sits
  exactly on 50%, it is the midpoint between that time and the next drop. Its 95% interval
  is read the same way off the lower and upper confidence curves. "Not reached" means the
  curve, or that confidence curve, never gets down to 50%.
- **1-, 3- and 5-year survival** is the curve's value at 12, 36 and 60 months. It shows "–"
  past the end of a group's follow-up.
- **Median follow-up** uses the reverse Kaplan–Meier method (censoring is treated as the
  event and events as censoring), over all patients shown.

## 4. Log-rank test

The test compares the whole curves across the $k$ groups (unweighted, R's
`survdiff(rho = 0)`). At each event time, with $d$ events among $n$ at risk in total and
$n_j$ at risk in group $j$:

- expected events in group $j$: $E_j = d\,n_j / n$
- covariance: $V_{jl} = \dfrac{d\,(n-d)}{n-1}\,\dfrac{n_j}{n}\left(\delta_{jl} - \dfrac{n_l}{n}\right)$

Observed minus expected counts and the covariances are summed over event times. The
statistic $(O-E)^\top V^{-1} (O-E)$, using any $k-1$ groups, is compared with a chi-square
distribution on $k-1$ degrees of freedom.

**Stratified by cancer type.** With that option on, $O-E$ and $V$ are computed inside each
cancer type and then added. Patients are only ever compared with patients who have the same
cancer type, so differences in baseline survival between cancer types don't drive the result.
The plotted curves are not stratified, so a stratified p-value can disagree with what the
pooled curves appear to show.

## 5. Cox proportional hazards models

Two places use a Cox model:

- **Hazard ratios in the group table**: a model with the comparison groups only, each
  against the first group. Unadjusted, but stratified if the option is on.
- **The Cox panel**: the comparison groups (optional) plus any variables you add.

### Fitting

- Partial likelihood with **Efron's approximation** for tied event times (R's default).
- Maximized by Newton–Raphson from all coefficients at 0, with step-halving if a step lowers
  the likelihood. It stops when the log-likelihood changes by less than 1 part in 10⁹, or
  after 40 iterations.
- **Stratification by cancer type** gives each cancer type its own baseline hazard; the
  hazard ratios are shared. Cancer type then can't also be a covariate.
- Standard errors come from the inverse of the observed information matrix.

### What is reported

| Output | Computation |
|---|---|
| Hazard ratio | $e^{\hat\beta}$ |
| 95% CI | $e^{\hat\beta \pm 1.96\,\mathrm{SE}}$ |
| p (per term) | Wald test: $z = \hat\beta / \mathrm{SE}$, two-sided normal |
| Likelihood-ratio test | twice the gain in log partial likelihood over a model with no terms; chi-square on the number of terms |

### How variables are coded

- **Numeric**, your choice of scale:
  - per unit: the raw value;
  - per SD: divided by the standard deviation among the patients in the model;
  - per doubling: $\log_2(\text{value} + 1)$. Negative values are left out.
- **Categorical**: one indicator per category against a reference. The reference is the most
  common category (for the comparison groups, the first group).
- **Missing values**: complete-case analysis. A patient missing any variable in the model is
  left out, and the count is reported under the table.
- **Small categories**: categories with fewer than 10 patients in the model are left out,
  since their hazard ratios can't be estimated reliably.

### When a result is withheld

- "Not estimable" appears when a coefficient's standard error is extreme (10 or more on the
  log scale). This usually means a category has no events.
- A warning appears if the model did not converge.
- The model is not fitted with fewer events than terms + 1, or with more than 40 terms.

## 6. Results for each cancer type

When the cohort contains more than one cancer type, a table under the plot repeats the
analysis inside each cancer type. It is meant for spotting cancer types worth a closer
look; clicking one restricts the cohort to it.

- **Groups** are the ones already defined for the plot. Numeric cut points therefore come
  from the whole cohort unless "Compute cut points within each cancer type" is on, in which
  case each type is split at its own values.
- **Per-type statistics**: patients, events, and for each group its size and Kaplan–Meier
  median (shown for up to 4 groups). With no comparison selected, the table instead gives
  the median with its 95% interval and 1-, 3- and 5-year survival.
- **Log-rank p** is the unstratified test within that cancer type, over the groups that have
  patients there. It is blank when only one group is present or there are no events.
- **Hazard ratio** (two groups only) is from an unadjusted Cox model, second group relative
  to the first. It is blank when it can't be estimated.
- **Adjusted p** is the Benjamini–Hochberg false-discovery-rate adjustment across the cancer
  types that were tested. Rows with adjusted p below 0.05 are shaded. It corrects only for
  the number of cancer types in this one table, not for other variables, cut points or
  endpoints you have tried.
- The "Stratify tests by cancer type" and "Censor follow-up at" options: stratification has
  no effect inside a single cancer type; censoring applies as in the main analysis.

Small cancer types have few events (see the README), so a large p-value there often means
too little data, not no effect.

## 7. What is not done

Keep these in mind when reading results:

- **No check of the proportional-hazards assumption** (for example Schoenfeld residuals).
- **No correction for multiple testing.** Trying many variables or cut points and reporting
  the best one overstates significance.
- **No competing-risks analysis**, no left truncation, no time-varying covariates.
- **Treatment comparisons are biased.** Therapy is recorded during follow-up but treated as
  known at diagnosis, so patients had to survive long enough to be treated (immortal-time
  bias). "No record" also doesn't mean untreated.
- **Data limitations** are listed in the [README](../README.md#known-limitations): short
  follow-up and few events for some cancer types, and incomplete staging and therapy data.

## Validation

`tests/stats.test.mjs` compares `stats.js` with R's `survival` package on 2,548 patients
(BRCA, LUAD, COAD and KIRC with stage and age recorded):

| Check | R reference | Agreement required |
|---|---|---|
| Kaplan–Meier curve, 95% limits, numbers at risk, median and its interval, per stage | `survfit`, `quantile` | 1e-9 or tighter |
| Log-rank across 4 stages, plain and stratified by cancer type | `survdiff` | 1e-8 |
| Cox coefficients, standard errors, log-likelihoods, likelihood-ratio statistic, plain and stratified | `coxph` (Efron) | 1e-6 or tighter |
| Chi-square tail probabilities | `pchisq` | 1e-9 |

```sh
Rscript tests/make_fixtures.R      # regenerate the R reference values
node --test tests/stats.test.mjs
```
