# ggsurvplot-extended

> Publication-ready Kaplan–Meier survival plots with embedded statistic tables — built on `{survival}`, `{survminer}`, and `{patchwork}`.

---

## Overview

This repository provides three R functions that extend the standard `ggsurvplot()` workflow into a single, reproducible call:

| Function | Role |
|---|---|
| `surv.plot()` | **Main entry point.** Fits the survival model, builds the KM curve, and assembles the final composite figure. |
| `ggsurvhrplot()` | **Lower-level compositor.** Combines a `ggsurvplot` object, a risk table, and an inset statistic table (HR / median) into one `patchwork` figure. Called internally by `surv.plot()` but can also be used directly. |
| `customize_labels()` | **Internal helper.** Applies font styling to `ggsurvplot` risk-table panels using `{ggtext}` markdown elements. |

The typical output is a composite figure with:

- A Kaplan–Meier curve with optional confidence bands and median survival reference line
- An inset table showing HR (95 % CI), p-value, and/or median survival per group
- A risk table (number at risk / cumulative events) below the curve

---

## Installation

No package installation is required — source the file directly:

```r
source("surv.plot_functions.R")
```

### Dependencies

Install any missing packages with `install.packages()`:

```r
install.packages(c(
  "survival", "survminer", "tidyverse", "ggplot2",
  "gridExtra", "patchwork", "grid", "gtable",
  "ggpubr", "ggtext"
))
```

| Package | Purpose |
|---|---|
| `survival` | `survfit()`, `coxph()`, `survreg()` |
| `survminer` | `ggsurvplot()`, `customize_labels()` internals |
| `ggplot2` | Plot construction and `annotation_custom()` |
| `patchwork` | Stacking KM curve and risk table vertically |
| `gridExtra` | `tableGrob()` for the inset statistic table |
| `gtable` | Adding the colour-segment column to the inset table |
| `ggpubr` | Font parsing inside `customize_labels()` |
| `ggtext` | Markdown-aware axis labels |

---

## Quick start

```r
library(survival)
source("surv_plot_functions.R")

df <- lung |>
  transform(
    sex        = factor(sex, levels = c(1, 2), labels = c("Male", "Female")),
    status_bin = as.numeric(status == 2)
  ) |>
  na.omit()

surv.plot(
  surv.formula = Surv(time, status_bin) ~ sex,
  dat          = df,
  labels       = c("Male", "Female"),
  palette      = c("dodgerblue3", "firebrick"),

  ggsurvplot.args = list(
    risk.table    = "nrisk_cumevents",
    xlab          = "Time (days)",
    ylab          = "Survival probability (%)",
    surv.scale    = "percent",
    conf.int      = TRUE,
    break.time.by = 90,
    xlim          = c(0, 900)
  ),

  tableplot.args = list(
    median.show = TRUE,
    position    = "topright"
  )
)
```

---

## Function reference

---

### `surv.plot()`

The main user-facing function. Fits the survival curve and statistical model, applies plot customisation, and returns a composite `patchwork` figure.

```r
surv.plot(
  surv.formula       = NULL,
  weights            = NULL,
  dat                = NULL,
  table.df           = NULL,
  labels             = NULL,
  p.value_df         = FALSE,
  labels.table.unadj = FALSE,
  aft.model          = FALSE,
  palette            = NULL,
  survfit.args       = list(),
  model.args         = list(),
  tableplot.args     = list(),
  ggsurvplot.args    = list()
)
```

#### Arguments

| Argument | Type | Default | Description |
|---|---|---|---|
| `surv.formula` | `formula` | `NULL` | A `Surv()` formula, e.g. `Surv(time, status) ~ group`. Used for both `survfit()` and the statistical model. **Required.** |
| `weights` | `numeric vector` | `NULL` | Optional observation weights (e.g. IPTW). Passed to both `survfit()` and `coxph()` / `survreg()`. |
| `dat` | `data.frame` | `NULL` | The analysis dataset. **Required.** |
| `table.df` | `data.frame` | `NULL` | Pre-computed statistics to display in the inset table instead of fitting a model. Must contain columns `HR`, `conf.low`, `conf.high`, and optionally `p.val` (required when `p.value_df = TRUE`). One row per non-reference group. |
| `labels` | `character vector` | `NULL` | Display labels for each group, in the same order as the factor levels. Used in the legend, the risk table row names, and the inset table row names. |
| `p.value_df` | `logical` | `FALSE` | When `table.df` is supplied, controls whether the `p.val` column is shown in the inset table. Ignored when a model is fitted (p-value is always shown in that case). |
| `labels.table.unadj` | `logical` | `FALSE` | When `TRUE`, replaces the risk table with an **unweighted** risk table (always based on a plain `survfit()` regardless of weights). Useful for weighted analyses where the weighted risk table would be misleading. |
| `aft.model` | `logical` | `FALSE` | When `FALSE` (default), fits a Cox proportional hazards model via `coxph()`. When `TRUE`, fits a Weibull accelerated failure time model via `survreg(dist = "weibull")` and back-transforms coefficients to the hazard ratio scale. |
| `palette` | `character vector` | `NULL` | Colour palette passed to `ggsurvplot()` and used for the colour segments in the inset table. One colour per group. |
| `survfit.args` | `list` | `list()` | Additional named arguments forwarded to `survfit()`, e.g. `list(start.time = 30)` for a landmark analysis. |
| `model.args` | `list` | `list()` | Additional named arguments forwarded to `coxph()` or `survreg()`, e.g. `list(ties = "efron")`. |
| `tableplot.args` | `list` | `list()` | Named arguments forwarded to `ggsurvhrplot()` to control the inset table (position, size, margins, median display). See `ggsurvhrplot()` reference below for the full list. |
| `ggsurvplot.args` | `list` | `list()` | Named arguments forwarded directly to `ggsurvplot()`, e.g. `list(conf.int = TRUE, risk.table = "nrisk_cumevents", xlim = c(0, 900))`. See `?ggsurvplot` for all options. |

#### Return value

A `patchwork` composite object containing the KM plot (with embedded inset table) stacked above the risk table. Print it with `print()` or pass it to `ggsave()`.

#### Notes

- The KM plot always includes a horizontal dashed reference line at 50 % survival.
- When `table.df = NULL`, the model is fitted automatically. When `table.df` is supplied, no model is fitted and `table.df` is used directly.
- `survfit.args` affects only the KM curve (and median extraction). The statistical model always uses the full dataset.

---

### `ggsurvhrplot()`

Lower-level compositor. Takes an already-built `ggsurvplot` object and adds an inset statistic table, then stacks the risk table below using `patchwork`. Called internally by `surv.plot()` via `tableplot.args`, but can be used directly when you need full manual control.

```r
ggsurvhrplot(
  cox_model    = NULL,
  aft_model    = NULL,
  data.table   = NULL,
  p.value      = TRUE,
  survfit      = NULL,
  median.show  = FALSE,
  ggsurvplot   = NULL,
  ggsurvtable  = NULL,
  labels       = c("Trt A", "Trt B"),
  col          = c("dodgerblue3", "firebrick"),
  position     = "topright",
  x_margin     = 0.02,
  y_margin     = 0.02,
  table_width  = 0.38,
  table_height = 0.16
)
```

#### Arguments

**Statistics source** — exactly one of the first three must be supplied, automatically done by calling `surv.plot` according to the specification of `aft.model` or `table.df`:

| Argument | Type | Default | Description |
|---|---|---|---|
| `cox_model` | `coxph` object | `NULL` | A fitted Cox model. HR, 95 % CI, and p-value are extracted automatically. p-value is always displayed. |
| `aft_model` | `survreg` object | `NULL` | A fitted Weibull AFT model. Coefficients are back-transformed to the HR scale using the estimated scale parameter. p-value is always displayed. |
| `data.table` | `data.frame` | `NULL` | Pre-computed statistics. Required columns: `HR`, `conf.low`, `conf.high`. Optional column: `p.val` (required when `p.value = TRUE`). One row per non-reference group; the reference row (`Ref.`) is prepended automatically. |

**Table content:**

| Argument | Type | Default | Description |
|---|---|---|---|
| `p.value` | `logical` | `TRUE` | Whether to display the p-value column. Only consulted when `data.table` is provided; for `cox_model` and `aft_model` the p-value is always shown. |
| `survfit` | `survfit` object | `NULL` | The `survfit` object used to extract median survival times. **Required when `median.show = TRUE`**; ignored otherwise. |
| `median.show` | `logical` | `FALSE` | When `TRUE`, adds a `Median (95% CI)` column to the inset table. Medians are extracted from `survfit`. Groups whose median is never reached are displayed as `NR`. |

**Plot inputs:**

| Argument | Type | Default | Description |
|---|---|---|---|
| `ggsurvplot` | `ggsurvplot` object | `NULL` | The `ggsurvplot` object whose `$plot` panel will receive the inset table. **Required.** |
| `ggsurvtable` | `ggplot` object | `NULL` | The risk table panel to stack below the KM plot. If `NULL`, falls back to `ggsurvplot$table`. Pass `NULL` explicitly to suppress the risk table entirely. |

**Appearance:**

| Argument | Type | Default | Description |
|---|---|---|---|
| `labels` | `character vector` | `c("Trt A", "Trt B")` | Row labels for the inset table, one per group including the reference. Must match the order of groups in the survival model / `data.table`. |
| `col` | `character vector` | `c("dodgerblue3", "firebrick")` | Colours used for the line segments in the leftmost column of the inset table. Should match the KM curve palette. |

**Inset table placement:**

| Argument | Type | Default | Description |
|---|---|---|---|
| `position` | `character` | `"topright"` | Corner of the plot panel to anchor the inset table. One of `"topright"`, `"topleft"`, `"bottomright"`, `"bottomleft"`. |
| `x_margin` | `numeric` | `0.02` | Horizontal shift as a fraction of the x-axis range. **Positive → right, negative → left**, regardless of which corner is chosen. A value of `0` places the table flush against the 2 % built-in padding. |
| `y_margin` | `numeric` | `0.02` | Vertical shift as a fraction of the y-axis range. **Positive → up, negative → down**, regardless of which corner is chosen. |
| `table_width` | `numeric` | `0.38` | Width of the inset table as a fraction of the x-axis range. Increase when adding extra columns (e.g. `0.52` with `median.show = TRUE`, `0.60` with three columns). |
| `table_height` | `numeric` | `0.16` | Height of the inset table as a fraction of the y-axis range. Increase for more groups (e.g. `0.26` for three groups). |

#### Return value

A `patchwork` object: KM plot with the inset table embedded (via `annotation_custom`) stacked above the risk table with a 5:1 height ratio.

#### Notes on `x_margin` / `y_margin` sign convention

The sign is **position-independent** — it always means the same direction regardless of which corner is chosen:

```
x_margin > 0  →  shift RIGHT      x_margin < 0  →  shift LEFT
y_margin > 0  →  shift UP         y_margin < 0  →  shift DOWN
```

Example — nudge the table slightly left and down from the top-right corner:

```r
tableplot.args = list(position = "topright", x_margin = -0.06, y_margin = -0.10)
```

---

### `customize_labels()`

Internal helper that applies `{ggtext}` markdown font styling to one or more `ggplot` panels (typically `ggsurvplot$table`). You will not normally call this directly.

```r
customize_labels(
  p,
  font.title     = NULL,
  font.subtitle  = NULL,
  font.caption   = NULL,
  font.x         = NULL,
  font.y         = NULL,
  font.xtickslab = NULL,
  font.ytickslab = NULL
)
```

Each font argument accepts a vector of the form `c(size, face)` or `c(size, face, colour)`, parsed by `ggpubr:::.parse_font()`. For example, `c(10, "bold")` sets 10 pt bold text.

---

### `format_number()`

Internal helper that formats a numeric value for display in the inset table. Numbers with five or more digits before the decimal point are shown in scientific notation (`1e+05`); all others are rounded to two decimal places.

---

## Common recipes

### Weighted analysis with unadjusted risk table

```r
surv.plot(
  surv.formula       = Surv(time, status_bin) ~ group,
  dat                = df,
  weights            = df$iptw,
  labels             = c("Control", "Treatment"),
  palette            = c("dodgerblue3", "firebrick"),
  labels.table.unadj = TRUE,          # show unweighted counts below
  ggsurvplot.args    = list(
    risk.table.title = "Adjusted No. at Risk (No. Events)",
    surv.scale       = "percent",
    conf.int         = TRUE
  ),
  tableplot.args = list(median.show = TRUE, table_width = 0.52)
)
```

### Pre-computed HR from an adjusted model

```r
# Fit your own adjusted model externally, then pass the results in
adj_hr <- data.frame(HR = 1.45, conf.low = 1.10, conf.high = 1.91, p.val = 0.008)

surv.plot(
  surv.formula = Surv(time, status_bin) ~ group,
  dat          = df,
  table.df     = adj_hr,
  p.value_df   = TRUE,
  labels       = c("Control", "Treatment"),
  palette      = c("dodgerblue3", "firebrick"),
  ggsurvplot.args = list(surv.scale = "percent", conf.int = TRUE),
  tableplot.args  = list(median.show = TRUE, table_width = 0.52)
)
```

### AFT (Weibull) model

```r
surv.plot(
  surv.formula = Surv(time, status_bin) ~ group,
  dat          = df,
  labels       = c("Control", "Treatment"),
  palette      = c("dodgerblue3", "firebrick"),
  aft.model    = TRUE,
  ggsurvplot.args = list(surv.scale = "percent"),
  tableplot.args  = list(median.show = TRUE)
)
```

### Landmark analysis starting at day 30

```r
surv.plot(
  surv.formula = Surv(time, status_bin) ~ group,
  dat          = df,
  labels       = c("Control", "Treatment"),
  palette      = c("dodgerblue3", "firebrick"),
  survfit.args = list(start.time = 30),
  ggsurvplot.args = list(surv.scale = "percent", xlim = c(30, 900))
)
```

### Save to file

```r
p <- surv.plot(...)
ggsave("km_plot.pdf", plot = p, width = 10, height = 8)
ggsave("km_plot.png", plot = p, width = 10, height = 8, dpi = 300)
```

---

## File structure

```
.
├── surv.plot_functions.R   # All functions (surv.plot, ggsurvhrplot, helpers)
├── create_df.table.R       # An example how to create the HR table for interactions
├── test_surv.plot_full.R   # 27 tests covering every argument and code path
└── README.md
```

---

## Dependencies and versions

Developed and tested with:

| Package | Version |
|---|---|
| R | ≥ 4.2 |
| survival | ≥ 3.5 |
| survminer | ≥ 0.4.9 |
| ggplot2 | ≥ 3.4 |
| patchwork | ≥ 1.2 |
| gridExtra | ≥ 2.3 |
| gtable | ≥ 0.3 |
| ggpubr | ≥ 0.6 |
| ggtext | ≥ 0.1.2 |

---

## License

MIT
