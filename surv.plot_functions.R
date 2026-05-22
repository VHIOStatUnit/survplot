# =============================================================================
# Customised ggsurvplot
# =============================================================================
surv.plot <- function(
    surv.formula        = NULL,
    weights             = NULL,
    dat                 = NULL,
    table.df            = NULL,
    labels              = NULL,
    p.value_df          = FALSE,
    labels.table.unadj  = FALSE,
    aft.model           = FALSE,
    palette             = NULL,
    # --- explicit pass-through args for survfit() ---
    survfit.args        = list(),
    # --- explicit pass-through args for coxph() / survreg() ---
    model.args          = list(),
    # --- explicit pass-through args for ggsurvhrplot() ---
    tableplot.args         = list(),
    # --- all remaining args go to ggsurvplot() ---
    ggsurvplot.args      = list()) {
  
  # ── 0. Load packages ──────────────────────────────────────────────────────
  for (pkg in c("survminer", "tidyverse", "ggplot2", "survival", "gridExtra", "patchwork", "grid")) {
    if (!base::paste0("package:", pkg) %in% search()) {
      library(pkg, character.only = TRUE)
    }
  }
  
  # ── 1. Build survfit / model ───────────────────────────────────────────────
  base_survfit_args <- list(formula = surv.formula, data = dat)
  if (!is.null(weights)) base_survfit_args$weights <- weights
  svf <- do.call(survival::survfit, c(base_survfit_args, survfit.args))
  
  base_model_args <- list(formula = surv.formula, data = dat)
  if (!is.null(weights)) base_model_args$weights <- weights
  
  if (!aft.model) {
    cox.model <- do.call(survival::coxph, c(base_model_args, model.args))
    aft.mod   <- NULL
  } else {
    aft.mod   <- do.call(survival::survreg, c(base_model_args, list(dist = "weibull"), model.args))
    cox.model <- NULL
  }
  
  # ── 2. Main ggsurvplot ────────────────────────────────────────────────────
  ggsurvplot_args <- c(
    list(fit         = svf,
         data        = dat,
         legend.labs = labels,
         palette     = palette),
    ggsurvplot.args
  )
  
  psurv <- do.call(survminer::ggsurvplot, ggsurvplot_args)
  
  # ── 3. Customise the KM plot panel ────────────────────────────────────────
  psurv$plot <- psurv$plot +
    ggplot2::annotate("segment",
             x = 0, xend = 45, y = 50, yend = 50,
             linetype = "dashed", color = "black", linewidth = 0.2) +
    ggplot2::theme(
      legend.title   = ggplot2::element_text(face = "bold"),
      axis.line      = ggplot2::element_line(linewidth = 0.2),
      axis.title     = ggplot2::element_text(face = "bold"),
      axis.title.x   = ggplot2::element_text(face = "bold", margin = ggplot2::margin(t = 10)),
      axis.title.y   = ggplot2::element_text(face = "bold", margin = ggplot2::margin(r = 10)),
      plot.margin    = ggplot2::margin(6, 6, 6, 27)
    ) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(linewidth = 0.8)))
  
  # ── 4. Customise the risk table panel ─────────────────────────────────────
  if (!is.null(psurv$table)) {
    psurv$table <- customize_labels(
      psurv$table,
      font.title     = c(10, "bold"),
      font.ytickslab = c(10)
    )
    psurv$table$theme$plot.title$margin <- ggplot2::margin(0, 0, 6, -40)
    psurv$table$theme$text$hjust        <- 0
  }
  
  # ── 5. Optional: unadjusted risk table ────────────────────────────────────
  if (labels.table.unadj) {
    
    svf_unadj <- do.call(survival::survfit, list(formula = surv.formula, data = dat))
    
    if ("risk.table" %in% names(ggsurvplot.args)) {
      ggsurvplot.args$risk.table <- NULL
    }
    
    if ("risk.table.title" %in% names(ggsurvplot.args)) {
      ggsurvplot.args$risk.table.title <- NULL
    }
    
    ggsurvplot_unadj_args <- c(
      list(fit                  = svf_unadj,
           data                 = dat,
           risk.table           = "nrisk_cumevents",
           risk.table.title     = "Unadjusted No. at Risk (No. Events)",
           legend.labs          = labels,
           palette              = palette),
      ggsurvplot.args
    )
    psurv_unadj <- do.call(survminer::ggsurvplot, ggsurvplot_unadj_args)
    
    psurv_unadj$table <- customize_labels(
      psurv_unadj$table,
      font.title     = c(10, "bold"),
      font.ytickslab = c(10)
    )
    psurv_unadj$table$theme$plot.title$margin <- ggplot2::margin(0, 0, 6, -40)
    psurv_unadj$table$theme$text$hjust        <- 0
    
    active_table <- psurv_unadj$table
    
  } else {
    active_table <- psurv$table
  }
  
  # ── 6. Assemble the final composite plot ──────────────────────────────────
  # Build the common args that always go to ggsurvhrplot.
  # Old positional tab_* args are gone — user controls placement via
  # position / x_margin / y_margin / table_width / table_height in tableplot.args
  base_hrplot_args <- list(
    ggsurvplot  = psurv,
    ggsurvtable = active_table,
    labels      = labels,
    col         = palette
  )
  
  if (is.null(table.df)) {
    # Cox or AFT model path
    do.call(
      ggsurvhrplot,
      c(base_hrplot_args,
        list(cox_model = cox.model,
             aft_model = aft.mod,
             survfit   = svf),        # <-- pass survfit for median extraction
        tableplot.args)
    )
    
  } else {
    # Pre-supplied data.frame path
    do.call(
      ggsurvhrplot,
      c(base_hrplot_args,
        list(data.table = table.df,
             p.value    = p.value_df,
             survfit    = svf),        # <-- pass survfit for median extraction
        tableplot.args)
    )
  }
}


### =============================================================================
### Function to add the ggsurvplot KM, table and Cox model together
### =============================================================================

# -----------------------------------------------------------
# Helping function to extract the statistic of the Cox model
# -----------------------------------------------------------
# If the number of digits is greater than 4, the number will be illustrated
# in academic format
format_number <- function(x) {
  digits_before <- nchar(as.character(floor(abs(x))))
  
  if (digits_before > 4) {
    format(x, scientific = TRUE, digits = 1)
  } else {
    round(x, 2)
  }
}


ggsurvhrplot <- function(cox_model    = NULL,
                         aft_model    = NULL,
                         data.table   = NULL,      # alternative to cox/aft model
                         p.value      = TRUE,      # only used when data.table is provided
                         survfit      = NULL,      # survfit object for median extraction
                         median.show  = FALSE,     # add a Median (95% CI) column to the table
                         ggsurvplot   = NULL,
                         ggsurvtable  = NULL,
                         labels       = c("Trt A", "Trt B"),
                         col          = c("dodgerblue3", "firebrick"),
                         position     = "topright",
                         x_margin     = 0.02,
                         y_margin     = 0.02,
                         table_width  = 0.38,
                         table_height = 0.16) {
  
  position <- match.arg(position, c("topright", "topleft", "bottomright", "bottomleft"))
  
  # =====================================================
  # Input validation
  # =====================================================
  n_sources <- sum(!is.null(cox_model), !is.null(aft_model), !is.null(data.table))
  if (n_sources == 0) stop("Provide one of: cox_model, aft_model, or data.table.")
  if (n_sources  > 1) stop("Provide only one of: cox_model, aft_model, or data.table.")
  
  if (median.show && is.null(survfit)) {
    stop("A survfit object must be supplied via the 'survfit' argument when median.show = TRUE.")
  }
  
  # =====================================================
  # Extract statistics — Cox model
  # =====================================================
  if (!is.null(cox_model)) {
    cox_summary <- summary(cox_model)
    hr      <- sapply(cox_summary$coefficients[, "exp(coef)"], format_number)
    ci_low  <- sapply(cox_summary$conf.int[, "lower .95"],     format_number)
    ci_high <- sapply(cox_summary$conf.int[, "upper .95"],     format_number)
    pval     <- cox_summary$coefficients[, "Pr(>|z|)"]
    pval_txt <- ifelse(pval < 0.01, "<0.01", round(pval, 2))
    p.value  <- TRUE   # always show p-value from model
  }
  
  # =====================================================
  # Extract statistics — AFT model
  # =====================================================
  if (!is.null(aft_model)) {
    aft_summary <- summary(aft_model)
    names_var   <- names(coef(aft_model))[-1]
    scale       <- aft_model$scale
    coef_val <- vapply(names_var, function(v) coef(aft_model)[v],          numeric(1))
    se       <- vapply(names_var, function(v) sqrt(vcov(aft_model)[v, v]), numeric(1))
    hr      <- sapply(exp(-coef_val / scale),                format_number)
    ci_low  <- sapply(exp(-(coef_val + 1.96 * se) / scale), format_number)
    ci_high <- sapply(exp(-(coef_val - 1.96 * se) / scale), format_number)
    pval     <- vapply(names_var, function(v) aft_summary$table[v, "p"], numeric(1))
    pval_txt <- ifelse(pval < 0.01, "<0.01", round(pval, 2))
    p.value  <- TRUE
  }
  
  # =====================================================
  # Extract statistics — data.frame input
  # Expected columns: HR, conf.low, conf.high, (p.val optional)
  # First row is treated as reference and prepended automatically
  # =====================================================
  if (!is.null(data.table)) {
    hr      <- sapply(data.table$HR,        format_number)
    ci_low  <- sapply(data.table$conf.low,  format_number)
    ci_high <- sapply(data.table$conf.high, format_number)
    
    if (p.value) {
      if (!"p.val" %in% names(data.table))
        stop("data.table must contain a 'p.val' column when p.value = TRUE.")
      pval     <- data.table$p.val
      pval_txt <- ifelse(pval < 0.01, "<0.01", round(pval, 2))
    }
  }
  
  # =====================================================
  # Extract median survival times (if requested)
  # =====================================================
  if (median.show) {
    sv_sum <- summary(survfit)$table
    
    # summary()$table returns a matrix (multi-group) or named vector (single group)
    if (is.matrix(sv_sum)) {
      med_vals  <- sv_sum[, "median"]
      med_lo    <- sv_sum[, "0.95LCL"]
      med_hi    <- sv_sum[, "0.95UCL"]
    } else {
      # single-group survfit — wrap in a 1-row structure
      med_vals  <- sv_sum["median"]
      med_lo    <- sv_sum["0.95LCL"]
      med_hi    <- sv_sum["0.95UCL"]
    }
    
    # Helper: format a single median value / CI bound (NA → "NR")
    fmt_med <- function(x) ifelse(is.na(x), "NR", format_number(x))
    
    med_txt <- mapply(
      function(m, lo, hi)
        base::paste0(fmt_med(m), " (", fmt_med(lo), " - ", fmt_med(hi), ")"),
      med_vals, med_lo, med_hi
    )
  }
  
  # =====================================================
  # Build table data frame
  # =====================================================
  hr_col <- c("Ref.", base::paste0(hr, " (", ci_low, " - ", ci_high, ")"))
  
  if (p.value && median.show) {
    table_df <- data.frame(
      `Median (95%CI)` = med_txt,
      `HR (95%CI)` = hr_col,
      `p-value`    = c("-", pval_txt),
      check.names  = FALSE
    )
  } else if (p.value) {
    table_df <- data.frame(
      `HR (95%CI)` = hr_col,
      `p-value`    = c("-", pval_txt),
      check.names  = FALSE
    )
  } else if (median.show) {
    table_df <- data.frame(
      `Median (95%CI)` = med_txt,
      `HR (95%CI)`     = hr_col,
      check.names      = FALSE
    )
  } else {
    table_df <- data.frame(
      `HR (95%CI)` = hr_col,
      check.names  = FALSE
    )
  }
  
  tbl <- gridExtra::tableGrob(
    table_df,
    rows  = labels,
    theme = gridExtra::ttheme_minimal(
      core    = list(fg_params = list(cex = 0.7)),
      colhead = list(fg_params = list(fontface = "bold",  cex = 0.8)),
      rowhead = list(fg_params = list(fontface = "plain", cex = 0.75))
    )
  )
  
  # =====================================================
  # Inject color segments as extra leftmost gtable column
  # =====================================================
  n_groups <- length(labels)
  
  line_grobs <- c(
    list(grid::nullGrob()),
    lapply(seq_len(n_groups), function(i) {
      grid::segmentsGrob(
        x0 = unit(0.1, "npc"), x1 = unit(0.9, "npc"),
        y0 = unit(0.5, "npc"), y1 = unit(0.5, "npc"),
        gp = grid::gpar(col = col[i], lwd = 2.5)
      )
    })
  )
  
  tbl <- gtable::gtable_add_cols(tbl, widths = unit(1.5, "char"), pos = 0)
  for (row_i in seq_along(line_grobs)) {
    tbl <- gtable::gtable_add_grob(
      tbl,
      grobs = line_grobs[[row_i]],
      t = row_i, b = row_i,
      l = 1,     r = 1
    )
  }
  
  # =====================================================
  # Compute annotation_custom coordinates from axis ranges
  #
  # x_margin / y_margin sign convention (position-independent):
  #   positive x_margin → shift table RIGHT
  #   negative x_margin → shift table LEFT
  #   positive y_margin → shift table UP
  #   negative y_margin → shift table DOWN
  # =====================================================
  built   <- ggplot2::ggplot_build(ggsurvplot$plot)
  x_range <- built$layout$panel_params[[1]]$x.range
  y_range <- built$layout$panel_params[[1]]$y.range
  
  x_span <- diff(x_range)
  y_span <- diff(y_range)
  
  xm <- x_span * x_margin   # signed offset in data units
  ym <- y_span * y_margin   # signed offset in data units
  tw <- x_span * table_width
  th <- y_span * table_height
  
  if (position == "topright") {
    # anchor: top-right corner; positive xm pushes further right, positive ym pushes further up
    xmax <- x_range[2] - (x_span * 0.02) + xm;  xmin <- xmax - tw
    ymax <- y_range[2] - (y_span * 0.02) + ym;  ymin <- ymax - th
    
  } else if (position == "topleft") {
    # anchor: top-left corner; positive xm pushes further right, positive ym pushes further up
    xmin <- x_range[1] + (x_span * 0.02) + xm;  xmax <- xmin + tw
    ymax <- y_range[2] - (y_span * 0.02) + ym;  ymin <- ymax - th
    
  } else if (position == "bottomright") {
    # anchor: bottom-right corner; positive xm pushes further right, positive ym pushes further up
    xmax <- x_range[2] - (x_span * 0.02) + xm;  xmin <- xmax - tw
    ymin <- y_range[1] + (y_span * 0.02) + ym;  ymax <- ymin + th
    
  } else if (position == "bottomleft") {
    # anchor: bottom-left corner; positive xm pushes further right, positive ym pushes further up
    xmin <- x_range[1] + (x_span * 0.02) + xm;  xmax <- xmin + tw
    ymin <- y_range[1] + (y_span * 0.02) + ym;  ymax <- ymin + th
  }
  
  # =====================================================
  # Inject table and combine with risk table
  # =====================================================
  final_plot <- ggsurvplot$plot +
    ggplot2::annotation_custom(
      grob = tbl,
      xmin = xmin, xmax = xmax,
      ymin = ymin, ymax = ymax
    )
  
  bottom <- if (!is.null(ggsurvtable)) ggsurvtable else ggsurvplot$table
  if (!is.null(bottom)) final_plot / bottom + patchwork::plot_layout(heights = c(5, 1)) else final_plot + patchwork::plot_layout(heights = c(6))
}


# =============================================================================
# A function to customize the labels of the ggsurvplot table
# =============================================================================
customize_labels <- function (p, font.title = NULL,
                              font.subtitle = NULL, font.caption = NULL,
                              font.x = NULL, font.y = NULL, font.xtickslab = NULL,
                              font.ytickslab = NULL){
  original.p <- p
  if(ggplot2::is_ggplot(original.p)) list.plots <- list(original.p)
  else if(is.list(original.p)) list.plots <- original.p
  else stop("Can't handle an object of class ", class (original.p))
  .set_font <- function(font){
    font <- ggpubr:::.parse_font(font)
    ggtext::element_markdown (size = font$size, face = font$face, colour = font$color)
  }
  for(i in 1:length(list.plots)){
    p <- list.plots[[i]]
    if(ggplot2::is_ggplot(p)){
      if (!is.null(font.title)) p <- p + ggplot2::theme(plot.title = .set_font(font.title))
      if (!is.null(font.subtitle)) p <- p + ggplot2::theme(plot.subtitle = .set_font(font.subtitle))
      if (!is.null(font.caption)) p <- p + ggplot2::theme(plot.caption = .set_font(font.caption))
      if (!is.null(font.x)) p <- p + ggplot2::theme(axis.title.x = .set_font(font.x))
      if (!is.null(font.y)) p <- p + ggplot2::theme(axis.title.y = .set_font(font.y))
      if (!is.null(font.xtickslab)) p <- p + ggplot2::theme(axis.text.x = .set_font(font.xtickslab))
      if (!is.null(font.ytickslab)) p <- p + ggplot2::theme(axis.text.y = .set_font(font.ytickslab))
      list.plots[[i]] <- p
    }
  }
  if(ggplot2::is_ggplot(original.p)) list.plots[[1]]
  else list.plots
}