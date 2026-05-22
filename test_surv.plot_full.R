# =============================================================================
# Comprehensive test script for surv.plot() and ggsurvhrplot()
# Dataset: built-in `lung` (survival package) — no external data needed
#
# HOW TO RUN
#   1. Place this file in the same folder as surv_plot_functions.R
#   2. Open an R session and run:  source("test_surv_plot_full.R")
#   3. Each test prints a heading, draws the plot, and pauses 1 s.
#      If running interactively, plots appear in your graphics device.
#      To save all plots to PDF instead, set  SAVE_PDF <- TRUE  below.
# =============================================================================

SAVE_PDF <- FALSE          # set TRUE to write all plots to surv_plot_tests.pdf
if (SAVE_PDF) pdf("~/Desktop/VHIO/R code/survplot_test.pdf", width = 10, height = 8)

# ── 0. Packages & source ──────────────────────────────────────────────────────
pkgs <- c("survival", "survminer", "tidyverse", "ggplot2",
          "gridExtra", "patchwork", "grid", "gtable",
          "ggpubr", "ggtext")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}

source("~/Desktop/VHIO/R code/survplot_median.R")

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

# Convenience wrappers ─────────────────────────────────────────────────────────
run_test <- function(label, expr) {
  cat(sprintf("\n%s\n%s\n", label, strrep("=", nchar(label))))
  p <- tryCatch(eval(expr), error = function(e) {
    cat("  ERROR:", conditionMessage(e), "\n"); NULL
  })
  if (!is.null(p)) print(p)
  Sys.sleep(1)
  invisible(p)
}

# Common ggsurvplot.args reused across many tests
base_gsp <- list(
  risk.table       = "nrisk_cumevents",
  risk.table.title = "No. at Risk (No. Events)",
  xlab             = "Time (days)",
  ylab             = "Survival probability (%)",
  break.time.by    = 90,
  xlim             = c(0, 900),
  surv.scale       = "percent",
  conf.int         = TRUE
)

pal2 <- c("dodgerblue3", "firebrick")   # 2-group palette
pal4 <- c("dodgerblue3", "firebrick",   # 4-group palette
          "forestgreen",  "darkorange")


# =============================================================================
# BLOCK A — Model paths (Cox / AFT / data.table)
# =============================================================================

# ── A1. Default Cox, no extras ────────────────────────────────────────────────
run_test("A1: Cox model | defaults (no risk table, no median)", quote(
  surv.plot(
    surv.formula    = Surv(time, status_bin) ~ sex,
    dat             = df,
    labels          = c("Male", "Female"),
    palette         = pal2,
    ggsurvplot.args = list(
      xlab      = "Time (days)",
      ylab      = "Survival probability (%)",
      surv.scale = "percent"
    )
  )
))

# ── A2. Cox with risk table ───────────────────────────────────────────────────
run_test("A2: Cox model | with risk table", quote(
  surv.plot(
    surv.formula    = Surv(time, status_bin) ~ sex,
    dat             = df,
    labels          = c("Male", "Female"),
    palette         = pal2,
    ggsurvplot.args = base_gsp
  )
))

# ── A3. Cox with risk table + median.show ─────────────────────────────────────
run_test("A3: Cox model | risk table + median.show = TRUE", quote(
  surv.plot(
    surv.formula    = Surv(time, status_bin) ~ sex,
    dat             = df,
    labels          = c("Male", "Female"),
    palette         = pal2,
    ggsurvplot.args = base_gsp,
    tableplot.args  = list(
      median.show  = TRUE,
      table_width  = 0.52,
      table_height = 0.18
    )
  )
))

# ── A4. Cox without p-value (p.value suppressed via data.table path) ──────────
# NOTE: Cox path always shows p-value (forced inside ggsurvhrplot).
# To suppress it, use the data.table path with p.value_df = FALSE (see B-tests).

# ── A5. AFT (Weibull) model ───────────────────────────────────────────────────
run_test("A5: AFT (Weibull) model | with risk table", quote(
  surv.plot(
    surv.formula    = Surv(time, status_bin) ~ sex,
    dat             = df,
    labels          = c("Male", "Female"),
    palette         = pal2,
    aft.model       = TRUE,
    ggsurvplot.args = base_gsp
  )
))

# ── A6. AFT + median.show ─────────────────────────────────────────────────────
run_test("A6: AFT model | median.show = TRUE", quote(
  surv.plot(
    surv.formula    = Surv(time, status_bin) ~ sex,
    dat             = df,
    labels          = c("Male", "Female"),
    palette         = pal2,
    aft.model       = TRUE,
    ggsurvplot.args = base_gsp,
    tableplot.args  = list(
      median.show  = TRUE,
      table_width  = 0.52,
      table_height = 0.18
    )
  )
))

# ── A7. data.table path — p-value shown ───────────────────────────────────────
run_test("A7: data.table path | p.value_df = TRUE", quote(
  surv.plot(
    surv.formula    = Surv(time, status_bin) ~ sex,
    dat             = df,
    table.df        = data.frame(HR = 1.45, conf.low = 1.10,
                                 conf.high = 1.91, p.val = 0.008),
    p.value_df      = TRUE,
    labels          = c("Male", "Female"),
    palette         = pal2,
    ggsurvplot.args = base_gsp
  )
))

# ── A8. data.table path — p-value hidden ─────────────────────────────────────
run_test("A8: data.table path | p.value_df = FALSE (HR only)", quote(
  surv.plot(
    surv.formula    = Surv(time, status_bin) ~ sex,
    dat             = df,
    table.df        = data.frame(HR = 1.45, conf.low = 1.10, conf.high = 1.91),
    p.value_df      = FALSE,
    labels          = c("Male", "Female"),
    palette         = pal2,
    ggsurvplot.args = base_gsp
  )
))

# ── A9. data.table path + median.show ────────────────────────────────────────
run_test("A9: data.table path | p.value_df = TRUE + median.show = TRUE", quote(
  surv.plot(
    surv.formula    = Surv(time, status_bin) ~ sex,
    dat             = df,
    table.df        = data.frame(HR = 1.45, conf.low = 1.10,
                                 conf.high = 1.91, p.val = 0.008),
    p.value_df      = TRUE,
    labels          = c("Male", "Female"),
    palette         = pal2,
    ggsurvplot.args = base_gsp,
    tableplot.args  = list(
      median.show  = TRUE,
      table_width  = 0.60,
      table_height = 0.18
    )
  )
))

# ── A10. data.table path — median only, no p-value ───────────────────────────
run_test("A10: data.table path | p.value_df = FALSE + median.show = TRUE", quote(
  surv.plot(
    surv.formula    = Surv(time, status_bin) ~ sex,
    dat             = df,
    table.df        = data.frame(HR = 1.45, conf.low = 1.10, conf.high = 1.91),
    p.value_df      = FALSE,
    labels          = c("Male", "Female"),
    palette         = pal2,
    ggsurvplot.args = base_gsp,
    tableplot.args  = list(
      median.show  = TRUE,
      table_width  = 0.52,
      table_height = 0.18
    )
  )
))


# =============================================================================
# BLOCK B — labels.table.unadj (unadjusted risk table)
# =============================================================================

# ── B1. Cox + unadjusted risk table ──────────────────────────────────────────
run_test("B1: Cox model | labels.table.unadj = TRUE", quote(
  surv.plot(
    surv.formula       = Surv(time, status_bin) ~ sex,
    dat                = df,
    labels             = c("Male", "Female"),
    palette            = pal2,
    labels.table.unadj = TRUE,
    ggsurvplot.args    = modifyList(base_gsp, list(
      risk.table.title = "Adjusted No. at Risk (No. Events)"
    ))
  )
))

# ── B2. Cox + unadjusted risk table + median.show ────────────────────────────
run_test("B2: Cox model | labels.table.unadj = TRUE + median.show = TRUE", quote(
  surv.plot(
    surv.formula       = Surv(time, status_bin) ~ sex,
    dat                = df,
    labels             = c("Male", "Female"),
    palette            = pal2,
    labels.table.unadj = TRUE,
    ggsurvplot.args    = modifyList(base_gsp, list(
      risk.table.title = "Adjusted No. at Risk (No. Events)"
    )),
    tableplot.args = list(
      median.show  = TRUE,
      table_width  = 0.52,
      table_height = 0.18
    )
  )
))


# =============================================================================
# BLOCK C — Position corners
# =============================================================================

corners <- c("topright", "topleft", "bottomright", "bottomleft")

for (pos in corners) {
  run_test(sprintf("C: position = '%s'", pos), bquote(
    surv.plot(
      surv.formula    = Surv(time, status_bin) ~ sex,
      dat             = df,
      labels          = c("Male", "Female"),
      palette         = pal2,
      ggsurvplot.args = base_gsp,
      tableplot.args  = list(position = .(pos))
    )
  ))
}


# =============================================================================
# BLOCK D — x_margin / y_margin sign behaviour
# =============================================================================

# ── D1. Negative x_margin → shift LEFT from topright corner ──────────────────
run_test("D1: x_margin = -0.10 (shift LEFT)", quote(
  surv.plot(
    surv.formula    = Surv(time, status_bin) ~ sex,
    dat             = df,
    labels          = c("Male", "Female"),
    palette         = pal2,
    ggsurvplot.args = base_gsp,
    tableplot.args  = list(position = "topright", x_margin = -0.10)
  )
))

# ── D2. Positive x_margin → shift RIGHT from topleft corner ──────────────────
run_test("D2: x_margin = +0.10 (shift RIGHT)", quote(
  surv.plot(
    surv.formula    = Surv(time, status_bin) ~ sex,
    dat             = df,
    labels          = c("Male", "Female"),
    palette         = pal2,
    ggsurvplot.args = base_gsp,
    tableplot.args  = list(position = "topleft", x_margin = 0.10)
  )
))

# ── D3. Negative y_margin → shift DOWN from topright corner ──────────────────
run_test("D3: y_margin = -0.15 (shift DOWN)", quote(
  surv.plot(
    surv.formula    = Surv(time, status_bin) ~ sex,
    dat             = df,
    labels          = c("Male", "Female"),
    palette         = pal2,
    ggsurvplot.args = base_gsp,
    tableplot.args  = list(position = "topright", y_margin = -0.15)
  )
))

# ── D4. Positive y_margin → shift UP from bottomright corner ─────────────────
run_test("D4: y_margin = +0.15 (shift UP)", quote(
  surv.plot(
    surv.formula    = Surv(time, status_bin) ~ sex,
    dat             = df,
    labels          = c("Male", "Female"),
    palette         = pal2,
    ggsurvplot.args = base_gsp,
    tableplot.args  = list(position = "bottomright", y_margin = 0.15)
  )
))

# ── D5. Combined: negative x + negative y ────────────────────────────────────
run_test("D5: x_margin = -0.08, y_margin = -0.10 (left + down from topright)", quote(
  surv.plot(
    surv.formula    = Surv(time, status_bin) ~ sex,
    dat             = df,
    labels          = c("Male", "Female"),
    palette         = pal2,
    ggsurvplot.args = base_gsp,
    tableplot.args  = list(position = "topright", x_margin = -0.08, y_margin = -0.10)
  )
))

# ── D6. Combined: positive x + positive y ────────────────────────────────────
run_test("D6: x_margin = +0.08, y_margin = +0.10 (right + up from bottomleft)", quote(
  surv.plot(
    surv.formula    = Surv(time, status_bin) ~ sex,
    dat             = df,
    labels          = c("Male", "Female"),
    palette         = pal2,
    ggsurvplot.args = base_gsp,
    tableplot.args  = list(position = "bottomleft", x_margin = 0.08, y_margin = 0.10)
  )
))


# =============================================================================
# BLOCK E — table_width / table_height
# =============================================================================

# ── E1. Narrow table ─────────────────────────────────────────────────────────
run_test("E1: table_width = 0.25, table_height = 0.12 (compact)", quote(
  surv.plot(
    surv.formula    = Surv(time, status_bin) ~ sex,
    dat             = df,
    labels          = c("Male", "Female"),
    palette         = pal2,
    ggsurvplot.args = base_gsp,
    tableplot.args  = list(table_width = 0.25, table_height = 0.12)
  )
))

# ── E2. Wide table ───────────────────────────────────────────────────────────
run_test("E2: table_width = 0.60, table_height = 0.22 (spacious)", quote(
  surv.plot(
    surv.formula    = Surv(time, status_bin) ~ sex,
    dat             = df,
    labels          = c("Male", "Female"),
    palette         = pal2,
    ggsurvplot.args = base_gsp,
    tableplot.args  = list(table_width = 0.60, table_height = 0.22)
  )
))


# =============================================================================
# BLOCK F — Weights
# =============================================================================

# Simulate IPTW-style weights (random here, just to exercise the code path)
set.seed(42)
df$weight_sim <- runif(nrow(df), 0.5, 2.0)

run_test("F1: Cox model with survey weights", quote(
  surv.plot(
    surv.formula    = Surv(time, status_bin) ~ sex,
    dat             = df,
    weights         = df$weight_sim,
    labels          = c("Male", "Female"),
    palette         = pal2,
    ggsurvplot.args = base_gsp,
    tableplot.args  = list(median.show = TRUE, table_width = 0.52)
  )
))

run_test("F2: Cox + weights + unadjusted risk table", quote(
  surv.plot(
    surv.formula       = Surv(time, status_bin) ~ sex,
    dat                = df,
    weights            = df$weight_sim,
    labels             = c("Male", "Female"),
    palette            = pal2,
    labels.table.unadj = TRUE,
    ggsurvplot.args    = modifyList(base_gsp, list(
      risk.table.title = "Adjusted No. at Risk (No. Events)"
    )),
    tableplot.args = list(median.show = TRUE, table_width = 0.52)
  )
))


# =============================================================================
# BLOCK G — Three-group model (ECOG 0 / 1 / 2)
# =============================================================================

df3 <- df %>% filter(ecog %in% c("ECOG 0", "ECOG 1", "ECOG 2")) %>%
  droplevels()
pal3 <- c("dodgerblue3", "firebrick", "forestgreen")

run_test("G1: Three groups (ECOG) | Cox + median.show", quote(
  surv.plot(
    surv.formula    = Surv(time, status_bin) ~ ecog,
    dat             = df3,
    labels          = c("ECOG 0", "ECOG 1", "ECOG 2"),
    palette         = pal3,
    ggsurvplot.args = modifyList(base_gsp, list(conf.int = FALSE)),
    tableplot.args  = list(
      median.show  = TRUE,
      position     = "topright",
      table_width  = 0.52,
      table_height = 0.26
    )
  )
))

run_test("G2: Three groups (ECOG) | data.table path + median.show", quote(
  surv.plot(
    surv.formula    = Surv(time, status_bin) ~ ecog,
    dat             = df3,
    table.df        = data.frame(
      HR        = c(1.50, 2.10),
      conf.low  = c(0.95, 1.35),
      conf.high = c(2.35, 3.28),
      p.val     = c(0.08, 0.001)
    ),
    p.value_df      = TRUE,
    labels          = c("ECOG 0", "ECOG 1", "ECOG 2"),
    palette         = pal3,
    ggsurvplot.args = modifyList(base_gsp, list(conf.int = FALSE)),
    tableplot.args  = list(
      median.show  = TRUE,
      position     = "topright",
      table_width  = 0.60,
      table_height = 0.26
    )
  )
))


# =============================================================================
# BLOCK H — survfit.args and model.args pass-through
# =============================================================================

# ── H1. survfit.args: start.time ─────────────────────────────────────────────
run_test("H1: survfit.args — start.time = 30 (landmark analysis at day 30)", quote(
  surv.plot(
    surv.formula    = Surv(time, status_bin) ~ sex,
    dat             = df,
    labels          = c("Male", "Female"),
    palette         = pal2,
    survfit.args    = list(start.time = 30),
    ggsurvplot.args = base_gsp
  )
))

# ── H2. model.args: robust standard errors for Cox ───────────────────────────
run_test("H2: model.args — robust SEs in Cox (ties = 'efron')", quote(
  surv.plot(
    surv.formula    = Surv(time, status_bin) ~ sex,
    dat             = df,
    labels          = c("Male", "Female"),
    palette         = pal2,
    model.args      = list(ties = "efron"),
    ggsurvplot.args = base_gsp
  )
))


# =============================================================================
# BLOCK I — ggsurvplot.args pass-through
# =============================================================================

# ── I1. No confidence interval ────────────────────────────────────────────────
run_test("I1: conf.int = FALSE", quote(
  surv.plot(
    surv.formula    = Surv(time, status_bin) ~ sex,
    dat             = df,
    labels          = c("Male", "Female"),
    palette         = pal2,
    ggsurvplot.args = modifyList(base_gsp, list(conf.int = FALSE))
  )
))

# ── I2. Cumulative events instead of survival ─────────────────────────────────
run_test("I2: fun = 'event' (cumulative incidence)", quote(
  surv.plot(
    surv.formula    = Surv(time, status_bin) ~ sex,
    dat             = df,
    labels          = c("Male", "Female"),
    palette         = pal2,
    ggsurvplot.args = modifyList(base_gsp, list(
      fun       = "event",
      ylab      = "Cumulative incidence (%)",
      conf.int  = FALSE
    )),
    tableplot.args = list(position = "bottomright")
  )
))

# ── I3. Log-log scale ─────────────────────────────────────────────────────────
run_test("I3: fun = 'cloglog' (log-log scale)", quote(
  surv.plot(
    surv.formula    = Surv(time, status_bin) ~ sex,
    dat             = df,
    labels          = c("Male", "Female"),
    palette         = pal2,
    ggsurvplot.args = list(
      fun       = "cloglog",
      xlab      = "log(Time)",
      ylab      = "log(-log(S(t)))",
      risk.table = FALSE,
      conf.int  = FALSE
    ),
    tableplot.args = list(position = "bottomright")
  )
))

# ── I4. Custom break.time.by ──────────────────────────────────────────────────
run_test("I4: break.time.by = 180 (6-month ticks)", quote(
  surv.plot(
    surv.formula    = Surv(time, status_bin) ~ sex,
    dat             = df,
    labels          = c("Male", "Female"),
    palette         = pal2,
    ggsurvplot.args = modifyList(base_gsp, list(break.time.by = 180))
  )
))


# =============================================================================
# BLOCK J — Edge / error cases
# =============================================================================

# ── J1. median.show = TRUE without survfit → should error gracefully ──────────
cat("\nJ1: median.show = TRUE called directly on ggsurvhrplot without survfit\n",
    strrep("=", 65), "\n")
tryCatch(
  ggsurvhrplot(
    cox_model   = coxph(Surv(time, status_bin) ~ sex, data = df),
    ggsurvplot  = ggsurvplot(survfit(Surv(time, status_bin) ~ sex, data = df),
                             data = df),
    labels      = c("Male", "Female"),
    col         = pal2,
    median.show = TRUE      # survfit not supplied — should error
  ),
  error = function(e) cat("  Expected error caught:", conditionMessage(e), "\n")
)

# ── J2. data.table path with p.value = TRUE but missing p.val column ──────────
cat("\nJ2: data.table missing p.val column when p.value_df = TRUE\n",
    strrep("=", 60), "\n")
tryCatch(
  surv.plot(
    surv.formula    = Surv(time, status_bin) ~ sex,
    dat             = df,
    table.df        = data.frame(HR = 1.45, conf.low = 1.10, conf.high = 1.91),
    p.value_df      = TRUE,     # p.val column missing → error expected
    labels          = c("Male", "Female"),
    palette         = pal2,
    ggsurvplot.args = base_gsp
  ),
  error = function(e) cat("  Expected error caught:", conditionMessage(e), "\n")
)

# ── J3. Sanity-check: median values match raw survfit output ──────────────────
cat("\nJ3: Sanity check — median from survfit directly\n",
    strrep("=", 50), "\n")
svf_check <- survfit(Surv(time, status_bin) ~ sex, data = df)
print(summary(svf_check)$table[, c("median", "0.95LCL", "0.95UCL")])
cat("  These values should match the 'Median (95%CI)' column in tests A3, A6, A9.\n")


# =============================================================================
# Done
# =============================================================================
if (SAVE_PDF) {
  dev.off()
  cat("\nAll plots saved to surv_plot_tests.pdf\n")
} else {
  cat("\nAll tests complete.\n")
}
