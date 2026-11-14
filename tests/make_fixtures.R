# Reference results from R's survival package for tests/stats.test.mjs.
# Run from the repo root: Rscript tests/make_fixtures.R
suppressMessages({library(arrow); library(data.table); library(survival); library(jsonlite)})
d <- as.data.table(read_parquet("data/processed/survivome.parquet"))
d <- d[cancer_type %in% c("BRCA", "LUAD", "COAD", "KIRC") & !is.na(OS) & !is.na(OS_time) &
         !is.na(stage) & !is.na(age), .(cancer_type, time = OS_time, event = OS, stage, age, sex)]
d[, g := as.integer(factor(stage)) - 1L]
d[, male := as.integer(sex == "male")]

km <- lapply(sort(unique(d$g)), function(k) {
  f <- survfit(Surv(time, event) ~ 1, d[g == k])
  q <- quantile(f, 0.5)
  list(time = f$time, surv = f$surv, lower = f$lower, upper = f$upper, n_risk = f$n.risk,
       median = unname(q$quantile), median_lower = unname(q$lower), median_upper = unname(q$upper))
})
lr <- survdiff(Surv(time, event) ~ g, d)
lrs <- survdiff(Surv(time, event) ~ g + strata(cancer_type), d)
cx <- coxph(Surv(time, event) ~ age + male + factor(g), d)
cxs <- coxph(Surv(time, event) ~ age + male + factor(g) + strata(cancer_type), d)
cox_out <- function(m) list(coef = unname(coef(m)), se = unname(sqrt(diag(vcov(m)))),
                            loglik = m$loglik, lr = unname(summary(m)$logtest["test"]))
out <- list(
  data = list(time = d$time, event = d$event, g = d$g, age = d$age, male = d$male,
              stratum = as.integer(factor(d$cancer_type))),
  km = km,
  logrank = list(chisq = lr$chisq, obs = as.numeric(lr$obs), exp = as.numeric(lr$exp)),
  logrank_strat = list(chisq = lrs$chisq),
  cox = cox_out(cx), cox_strat = cox_out(cxs),
  pchisq = list(x = c(0.5, 3.84, 20, 80), df = c(1, 1, 3, 5),
                p = pchisq(c(0.5, 3.84, 20, 80), c(1, 1, 3, 5), lower.tail = FALSE))
)
write_json(out, "tests/fixtures.json", digits = NA, auto_unbox = TRUE, na = "null")
cat("rows:", nrow(d), "\n")
