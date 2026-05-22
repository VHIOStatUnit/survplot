# ── 0. Packages & source ──────────────────────────────────────────────────────
pkgs <- c("survival", "survminer", "tidyverse", "ggplot2",
          "gridExtra", "patchwork", "grid", "gtable",
          "ggpubr", "ggtext", "emmeans")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}

# ── 1. Data preparation ───────────────────────────────────────────────────────
# lung: NCCTG lung cancer trial, shipped with {survival}
#   time   : survival time (days)
#   status : 1 = censored, 2 = dead
#   sex    : 1 = Male, 2 = Female
#   ph.ecog: ECOG performance score 0-4 (0 = good)
#   age    : age in years

df <- lung %>%
  mutate(
    sex        = factor(sex,    levels = c(1, 2), labels = c("Male", "Female")),
    ecog       = factor(ph.ecog, levels = c(0, 1, 2, 3),
                        labels = c("ECOG 0", "ECOG 1", "ECOG 2", "ECOG 3")),
    age_grp    = factor(ifelse(age < 63, "Age < 63", "Age ≥ 63")),
    status_bin = as.numeric(status == 2)
  ) %>%
  drop_na(time, status_bin, sex, ecog, age_grp, age)


# ── 2. Building Cox Model ─────────────────────────────────────────────────────
cox.model <- coxph(Surv(time, status_bin) ~ sex * age_grp, data = df)


# ── 3. Creating data frame for the hazard ratio ───────────────────────────────
df_forest_cox<- data.frame(emmeans(cox.model, ~ sex * age_grp))


df_forest_cox <- df_forest_cox |>
  arrange(sex, age_grp) |>
  mutate(exp.coef = exp(emmean),
         exp.low.ci = exp(asymp.LCL),
         exp.up.ci = exp(asymp.UCL),
         "HR (95% CI)" = sprintf("%.2f (%.2f, %.2f)", exp.coef, exp.low.ci, exp.up.ci)) |>
  mutate(Subgroups = paste(sex, "+", age_grp))

df_forest_cox <- df_forest_cox |> slice(-1) |>
  arrange(sex, age_grp) |>
  select(Subgroups, exp.coef, SE, exp.low.ci, exp.up.ci, `HR (95% CI)`)

df.plot <- df_forest_cox |> 
  select(Subgroups, exp.coef, exp.low.ci, exp.up.ci) |>
  rename(HR = exp.coef,
         conf.low = exp.low.ci,
         conf.high = exp.up.ci)

df.plot
