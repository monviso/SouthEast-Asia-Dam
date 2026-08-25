# Step3D_plot_significant_GAM_effects_ggplot_v5.R
#
# Re-plot all statistically significant GAM smooth terms from selected Step3A
# models using ggplot2 + patchwork. This version avoids plot.gam() and fixes
# the common issue where smooth objects do not expose first.para/last.para by
# recovering coefficient indices from coefficient names.
#
# Outputs one multi-panel figure for each selected model:
#   - Tlake
#   - Toutlet
#   - Tdiff
#
# 2D distributed-lag smooths are plotted as raster + contours.
# 1D smooths are plotted as partial-effect lines with 95% confidence bands.

required_packages <- c(
  "mgcv", "dplyr", "readr", "tibble", "stringr", "purrr",
  "ggplot2", "patchwork", "scales"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Install missing packages first: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

library(mgcv)
library(dplyr)
library(readr)
library(tibble)
library(stringr)
library(purrr)
library(ggplot2)
library(patchwork)
library(scales)

# -------------------------------------------------------------------------
# User settings
# -------------------------------------------------------------------------

step3A_dir <- "E:/OneDrive/Documents/satData_AppEEars/model_outputs/step3A_main_models"
fit_search_dir <- step3A_dir

output_dir <- file.path(step3A_dir, "step3D_significant_effect_figures_ggplot_v5")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

alpha <- 0.05

# Current selected models from your LODO results.
# Set use_manual_selection <- FALSE to choose lowest RMSE automatically from
# model_comparison_all.csv / model_comparison_lag.csv.
use_manual_selection <- FALSE

manual_selected_models <- tibble::tribble(
  ~model,    ~response, ~lag_window,
  "Tlake",   "Tm_l",    "lag_24h",
  "Toutlet", "Tm_o",    "lag_4h",
  "Tdiff",   "T_diff",  "lag_24h"
)

# Keep these modest to avoid long rendering times.
surface_n_x <- 55
surface_n_y <- 40
smooth_n <- 130
range_trim <- 0.02

figure_width <- 32
panel_height <- 10
max_cols <- 3
png_dpi <- 300

# To avoid the elapsed-time problem, only PNG is saved by default.
# Switch to TRUE when the PNGs look right and you want vector output.
save_pdf <- FALSE
save_individual_panels <- FALSE

# Use the same colour range for all 2D lag-effect surfaces.
# This makes colour intensity comparable across meteorological-lag panels and
# lets patchwork collect a single shared legend.
use_global_2d_colour_scale <- TRUE

# Scope for the shared colour range:
#   "within_model" = one colour scale per model figure, calculated only from
#                    the significant 2D lag surfaces of that specific model.
#   "all_models"   = one colour scale per model figure.
# Default is "within_model" because this keeps the colour bar comparable within
# each response model, without being dominated by another response.
global_2d_colour_scope <- "within_model"

# -------------------------------------------------------------------------
# Labels and helpers
# -------------------------------------------------------------------------

clean_file_label <- function(x) {
  x |>
    as.character() |>
    stringr::str_replace_all("[^A-Za-z0-9_]+", "_") |>
    stringr::str_replace_all("^_|_$", "")
}

normalise_term <- function(x) {
  x |>
    as.character() |>
    stringr::str_replace_all("\\s+", "") |>
    stringr::str_replace_all("`", "")
}

pretty_model_label <- function(model, lag_window) {
  response_label <- dplyr::case_when(
    model == "Tlake" ~ "Lake WST",
    model == "Toutlet" ~ "Outlet WST",
    model == "Tdiff" ~ " WST difference",
    TRUE ~ model
  )
  
  lag_label <- ifelse(
    lag_window == "0h",
    "0 h",
    stringr::str_replace(lag_window, "lag_", "")
  )
  
  paste0(response_label, " (", lag_label, " model)")
}

pretty_term_label <- function(term) {
  term0 <- normalise_term(term)
  
  dplyr::case_when(
    str_detect(term0, "^te\\(Air_T_l,lag") ~ "Air temp.",
    str_detect(term0, "^te\\(prec_l,lag") ~ "Precipitation",
    str_detect(term0, "^te\\(wind_l,lag") ~ "Wind speed",
    str_detect(term0, "^te\\(shortW_l,lag") ~ "Shortwave rad.",
    str_detect(term0, "^te\\(pres_l,lag") ~ "Atm. pres.",
    
    str_detect(term0, "^s\\(Air_T") ~ "Air temp.",
    str_detect(term0, "^s\\(prec") ~ "Precipitation",
    str_detect(term0, "^s\\(wind") ~ "Wind speed",
    str_detect(term0, "^s\\(shortW") ~ "Shortwave rad.",
    str_detect(term0, "^s\\(pres") ~ "Atm. pres.",
    
    str_detect(term0, "^s\\(doy") ~ "Day of year",
    str_detect(term0, "^s\\(m_LOC") ~ "Month",
    str_detect(term0, "^s\\(H_LOC") ~ "Local solar hour",
    
    str_detect(term0, "^s\\(Dam_hgt") ~ "Dam height",
    str_detect(term0, "^s\\(elevation") ~ "Elevation",
    str_detect(term0, "^s\\(mean_depth_m") ~ "Approximate mean depth",
    str_detect(term0, "^s\\(log_res_area") ~ "Reservoir area",
    str_detect(term0, "^s\\(log_res_capacity") ~ "Reservoir capacity",
    str_detect(term0, "^s\\(Res_area_kmsq") ~ "Reservoir area",
    str_detect(term0, "^s\\(Res_capacity_km3") ~ "Reservoir capacity",
    str_detect(term0, "^s\\(Res_capaci_kmcube") ~ "Reservoir capacity",
    
    TRUE ~ term
  )
}

axis_label <- function(var) {
  dplyr::case_when(
    var %in% c("Air_T", "Air_T_l") ~ "Air temperature (°C)",
    var %in% c("prec", "prec_l") ~ "Precipitation (mm)",
    var %in% c("wind", "wind_l") ~ "Wind speed (m s-1)",
    var %in% c("shortW", "shortW_l") ~ "Shortwave radiation (J m-2)",
    var %in% c("pres", "pres_l") ~ "Atmospheric pressure (Pa)",
    var == "lag" ~ "Lag (hours)",
    var == "doy" ~ "Day of year",
    var == "m_LOC" ~ "Month",
    var == "H_LOC" ~ "Local solar hour",
    var == "Dam_hgt" ~ "Dam height (m)",
    var == "elevation" ~ "Elevation (m a.s.l.)",
    var == "mean_depth_m" ~ "Approximate mean depth (m)",
    var == "log_res_area" ~ "Reservoir area (km2)",
    var == "log_res_capacity" ~ "Reservoir capacity (km3)",
    var == "Res_area_kmsq" ~ "Reservoir area (km2)",
    var == "Res_capacity_km3" ~ "Reservoir capacity (km3)",
    var == "Res_capaci_kmcube" ~ "Reservoir capacity (km3)",
    TRUE ~ var
  )
}

term_order_rank <- function(term) {
  term0 <- normalise_term(term)
  dplyr::case_when(
    str_detect(term0, "Air_T") ~ 1,
    str_detect(term0, "prec") ~ 2,
    str_detect(term0, "wind") ~ 3,
    str_detect(term0, "shortW") ~ 4,
    str_detect(term0, "pres") ~ 5,
    str_detect(term0, "doy|m_LOC|H_LOC") ~ 6,
    str_detect(term0, "Dam_hgt") ~ 7,
    str_detect(term0, "elevation") ~ 8,
    str_detect(term0, "mean_depth") ~ 9,
    str_detect(term0, "area|log_res_area") ~ 10,
    str_detect(term0, "capaci|capacity") ~ 11,
    TRUE ~ 99
  )
}

format_p <- function(p) {
  format.pval(p, digits = 2, eps = 1e-4)
}

find_p_col <- function(tab) {
  p_cols <- grep("p", names(tab), ignore.case = TRUE, value = TRUE)
  if (length(p_cols) == 0) return(NA_character_)
  p_cols[length(p_cols)]
}

match_smooth_index <- function(term, fit) {
  smooth_labels <- vapply(fit$smooth, function(x) x$label, character(1))
  
  idx <- match(term, smooth_labels)
  if (!is.na(idx)) return(idx)
  
  idx <- match(normalise_term(term), normalise_term(smooth_labels))
  if (!is.na(idx)) return(idx)
  
  NA_integer_
}

get_smooth_vars <- function(smooth_obj) {
  if (!is.null(smooth_obj$term) && length(smooth_obj$term) > 0) {
    return(as.character(smooth_obj$term))
  }
  
  if (!is.null(smooth_obj$term.names) && length(smooth_obj$term.names) > 0) {
    return(as.character(smooth_obj$term.names))
  }
  
  if (!is.null(smooth_obj$margin) && length(smooth_obj$margin) > 0) {
    vars <- unlist(lapply(smooth_obj$margin, get_smooth_vars), use.names = FALSE)
    vars <- unique(as.character(vars))
    if (length(vars) > 0) return(vars)
  }
  
  # Fallback from label, e.g. "te(Air_T_l,lag)" or "s(doy)".
  lab <- normalise_term(smooth_obj$label)
  inside <- stringr::str_match(lab, "^[a-z]+\\((.*)\\)$")[, 2]
  if (!is.na(inside)) {
    vars <- strsplit(inside, ",", fixed = TRUE)[[1]]
    vars <- vars[!str_detect(vars, "=")]
    vars <- vars[nzchar(vars)]
    if (length(vars) > 0) return(vars)
  }
  
  character(0)
}

coef_indices_for_smooth <- function(fit, smooth_obj) {
  first <- smooth_obj$first.para
  last <- smooth_obj$last.para
  
  if (!is.null(first) && !is.null(last) && length(first) == 1 && length(last) == 1) {
    return(seq.int(first, last))
  }
  
  cn <- names(stats::coef(fit))
  lab <- smooth_obj$label
  
  # Most mgcv smooth coefficients are named like "s(doy).1" or
  # "te(Air_T_l,lag).23".
  idx <- which(startsWith(cn, paste0(lab, ".")))
  
  # Fallback for minor whitespace/name differences.
  if (length(idx) == 0) {
    cn_norm <- normalise_term(cn)
    lab_norm <- normalise_term(lab)
    idx <- which(startsWith(cn_norm, paste0(lab_norm, ".")))
  }
  
  if (length(idx) == 0) {
    stop(
      "Could not identify coefficient columns for smooth label: ", lab,
      call. = FALSE
    )
  }
  
  idx
}

smooth_vcov <- function(fit, smooth_obj) {
  idx <- coef_indices_for_smooth(fit, smooth_obj)
  
  list(
    idx = idx,
    beta = stats::coef(fit)[idx],
    Vp = fit$Vp[idx, idx, drop = FALSE]
  )
}

get_model_var <- function(fit, var) {
  if (!is.null(fit$model) && var %in% names(fit$model)) {
    return(fit$model[[var]])
  }
  
  if (!is.null(fit$var.summary) && var %in% names(fit$var.summary)) {
    return(fit$var.summary[[var]])
  }
  
  NULL
}

finite_values <- function(x) {
  if (is.null(x)) return(numeric(0))
  
  if (is.matrix(x) || is.data.frame(x)) {
    x <- as.vector(as.matrix(x))
  }
  
  x <- suppressWarnings(as.numeric(x))
  x[is.finite(x)]
}

var_range <- function(fit, var, trim = range_trim) {
  vals <- finite_values(get_model_var(fit, var))
  
  if (var == "lag") {
    vals <- vals[is.finite(vals)]
    if (length(vals) > 0) return(range(vals, na.rm = TRUE))
    return(c(0, 48))
  }
  
  if (var == "doy") return(c(1, 366))
  if (var == "m_LOC") return(c(1, 12))
  if (var == "H_LOC") return(c(0, 24))
  
  if (length(vals) == 0) {
    stop("Could not find finite model values for variable: ", var, call. = FALSE)
  }
  
  if (trim > 0 && length(unique(vals)) > 10) {
    rr <- as.numeric(stats::quantile(vals, probs = c(trim, 1 - trim), na.rm = TRUE))
  } else {
    rr <- range(vals, na.rm = TRUE)
  }
  
  if (!all(is.finite(rr)) || diff(rr) == 0) {
    rr <- range(vals, na.rm = TRUE)
    if (diff(rr) == 0) rr <- rr + c(-0.5, 0.5)
  }
  
  rr
}

seq_for_var <- function(fit, var, n = smooth_n) {
  if (var == "lag") {
    vals <- sort(unique(finite_values(get_model_var(fit, var))))
    vals <- vals[is.finite(vals)]
    if (length(vals) > 0 && length(vals) <= n) return(vals)
  }
  
  rr <- var_range(fit, var)
  seq(rr[1], rr[2], length.out = n)
}

eval_smooth_direct <- function(fit, smooth_obj, newdata) {
  X <- mgcv::PredictMat(smooth_obj, newdata)
  vc <- smooth_vcov(fit, smooth_obj)
  
  if (ncol(X) != length(vc$beta)) {
    stop(
      "Design matrix has ", ncol(X), " columns but smooth has ",
      length(vc$beta), " coefficients for ", smooth_obj$label,
      call. = FALSE
    )
  }
  
  fit_val <- as.numeric(X %*% vc$beta)
  se_val <- sqrt(rowSums((X %*% vc$Vp) * X))
  
  tibble(
    fit = fit_val,
    se = se_val,
    lower = fit_val - 1.96 * se_val,
    upper = fit_val + 1.96 * se_val
  )
}

get_significant_smooths <- function(fit, alpha = 0.05) {
  sm <- summary(fit)
  
  if (is.null(sm$s.table)) return(tibble())
  
  smooth_tab <- as.data.frame(sm$s.table) |>
    rownames_to_column("term") |>
    as_tibble()
  
  p_col <- find_p_col(smooth_tab)
  if (is.na(p_col)) return(tibble())
  
  smooth_tab |>
    mutate(
      p_value = as.numeric(.data[[p_col]]),
      select_idx = vapply(term, match_smooth_index, integer(1), fit = fit),
      pretty_term = vapply(term, pretty_term_label, character(1)),
      term_type = if_else(str_detect(normalise_term(term), "^te\\("), "lag_surface", "smooth"),
      order_rank = vapply(term, term_order_rank, numeric(1))
    ) |>
    filter(
      is.finite(p_value),
      p_value < alpha,
      !str_detect(term, "dam_id"),
      !is.na(select_idx)
    ) |>
    arrange(term_type != "lag_surface", order_rank, p_value)
}

find_fit_file <- function(model, lag_window, search_dir = fit_search_dir) {
  all_rds <- list.files(search_dir, pattern = "\\.rds$", recursive = TRUE, full.names = TRUE)
  
  if (length(all_rds) == 0) {
    stop("No .rds model files found under: ", search_dir, call. = FALSE)
  }
  
  if (lag_window == "0h") {
    exact_names <- c(
      paste0("fit_zerolag_", model, ".rds"),
      paste0("fit_0h_", model, ".rds")
    )
  } else {
    exact_names <- c(
      paste0("fit_", lag_window, "_", model, ".rds"),
      paste0("fit_", model, "_", lag_window, ".rds")
    )
  }
  
  exact_hits <- all_rds[basename(all_rds) %in% exact_names]
  if (length(exact_hits) > 0) return(exact_hits[1])
  
  model_re <- clean_file_label(model)
  lag_re <- clean_file_label(lag_window)
  
  if (lag_window == "0h") {
    hits <- all_rds[
      str_detect(basename(all_rds), fixed(model_re)) &
        str_detect(basename(all_rds), regex("zerolag|0h", ignore_case = TRUE))
    ]
  } else {
    hits <- all_rds[
      str_detect(basename(all_rds), fixed(model_re)) &
        str_detect(basename(all_rds), fixed(lag_re))
    ]
  }
  
  if (length(hits) == 0) {
    stop(
      "Could not find fitted model RDS for model = ", model,
      ", lag_window = ", lag_window,
      ". Search folder: ", search_dir,
      call. = FALSE
    )
  }
  
  hits[1]
}

select_models_from_comparison <- function(step3A_dir) {
  candidates <- c(
    file.path(step3A_dir, "best_model_by_response.csv"),
    file.path(step3A_dir, "model_comparison_all.csv"),
    file.path(step3A_dir, "model_comparison_lag.csv")
  )
  
  candidates <- candidates[file.exists(candidates)]
  
  if (length(candidates) == 0) {
    warning("No comparison file found. Using manual_selected_models.")
    return(manual_selected_models)
  }
  
  comp <- read_csv(candidates[1], show_col_types = FALSE)
  
  if (!all(c("model", "lag_window") %in% names(comp))) {
    warning("Comparison file does not contain model and lag_window. Using manual_selected_models.")
    return(manual_selected_models)
  }
  
  metric_col <- if ("RMSE" %in% names(comp)) "RMSE" else if ("AIC" %in% names(comp)) "AIC" else NA_character_
  
  if (is.na(metric_col)) {
    warning("No RMSE or AIC column found. Using manual_selected_models.")
    return(manual_selected_models)
  }
  
  comp |>
    filter(model %in% c("Tlake", "Toutlet", "Tdiff"), is.finite(.data[[metric_col]])) |>
    group_by(model) |>
    arrange(.data[[metric_col]], .by_group = TRUE) |>
    slice(1) |>
    ungroup() |>
    mutate(
      response = case_when(
        model == "Tlake" ~ "Tm_l",
        model == "Toutlet" ~ "Tm_o",
        model == "Tdiff" ~ "T_diff",
        TRUE ~ NA_character_
      )
    ) |>
    select(model, response, lag_window)
}

make_panel_data_1d <- function(fit, smooth_obj, x_var) {
  x_grid <- seq_for_var(fit, x_var, n = smooth_n)
  
  plot_x <- x_grid
  if (x_var %in% c("log_res_area", "log_res_capacity")) {
    plot_x <- expm1(x_grid)
  }
  
  nd <- data.frame(x_grid)
  names(nd) <- x_var
  
  eff <- eval_smooth_direct(fit, smooth_obj, nd)
  
  bind_cols(
    tibble(x = plot_x, x_predict = x_grid),
    eff
  )
}

make_panel_data_2d <- function(fit, smooth_obj, x_var, y_var) {
  x_grid <- seq_for_var(fit, x_var, n = surface_n_x)
  y_grid <- seq_for_var(fit, y_var, n = surface_n_y)
  
  grd <- expand.grid(
    x = x_grid,
    y = y_grid,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  ) |>
    as_tibble()
  
  nd <- as.data.frame(grd)
  names(nd) <- c(x_var, y_var)
  
  eff <- eval_smooth_direct(fit, smooth_obj, nd)
  
  bind_cols(grd, eff)
}

plot_1d_effect <- function(panel_df, term, p_value, x_label) {
  ggplot(panel_df, aes(x = x, y = fit)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.22) +
    geom_line(linewidth = 0.75) +
    geom_hline(yintercept = 0, linewidth = 0.25, linetype = "dashed") +
    labs(
      title = paste0(pretty_term_label(term), "  (p = ", format_p(p_value), ")"),
      x = x_label,
      y = "Partial effect (°C)"
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 15),
      panel.grid.minor = element_blank(),
      axis.title =  element_text(size = 15),axis.text = element_text(size=12),
      legend.text = element_text(size=12), legend.title = element_text(size = 15)
    )
}

plot_2d_effect <- function(panel_df, term, p_value, x_label, y_label,
                           colour_limits = NULL) {
  if (is.null(colour_limits)) {
    lim <- suppressWarnings(max(abs(panel_df$fit), na.rm = TRUE))
    if (!is.finite(lim) || lim == 0) lim <- 1
    colour_limits <- c(-lim, lim)
  }
  
  ggplot(panel_df, aes(x = x, y = y, fill = fit)) +
    geom_raster(interpolate = TRUE) +
    geom_contour(aes(z = fit), colour = "grey25", linewidth = 0.16, alpha = 0.55, bins = 8) +
    scale_fill_gradient2(
      low = "#313695",
      mid = "#f7f7f7",
      high = "#a50026",
      midpoint = 0,
      limits = colour_limits,
      oob = scales::squish,
      name = "Partial\neffect (°C)"
    ) +
    coord_cartesian(expand = FALSE) +
    labs(
      title = paste0(pretty_term_label(term), "  (p = ", format_p(p_value), ")"),
      x = x_label,
      y = y_label
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 15),
      panel.grid = element_blank(),
      legend.position = "right",
      axis.title =  element_text(size = 15),axis.text = element_text(size=12),
      legend.text = element_text(size=12), legend.title = element_text(size = 15)
      
    )
}

make_effect_plot <- function(fit, sig_row, colour_limits_2d = NULL) {
  smooth_obj <- fit$smooth[[sig_row$select_idx]]
  vars <- get_smooth_vars(smooth_obj)
  
  if (length(vars) == 0) {
    stop("Could not identify variables for smooth: ", smooth_obj$label, call. = FALSE)
  }
  
  if (length(vars) == 1) {
    x_var <- vars[1]
    panel_df <- make_panel_data_1d(fit, smooth_obj, x_var)
    
    p <- plot_1d_effect(
      panel_df = panel_df,
      term = sig_row$term,
      p_value = sig_row$p_value,
      x_label = axis_label(x_var)
    )
    
    return(list(
      plot = p,
      data = panel_df,
      status = "ok",
      error = NA_character_,
      plot_type = "1d",
      x_var = x_var,
      y_var = NA_character_
    ))
  }
  
  # Use first two variables for 2D surfaces. For current lagged GAMs this is
  # meteorological variable and lag.
  x_var <- vars[1]
  y_var <- vars[2]
  panel_df <- make_panel_data_2d(fit, smooth_obj, x_var, y_var)
  
  p <- plot_2d_effect(
    panel_df = panel_df,
    term = sig_row$term,
    p_value = sig_row$p_value,
    x_label = axis_label(x_var),
    y_label = axis_label(y_var),
    colour_limits = colour_limits_2d
  )
  
  return(list(
    plot = p,
    data = panel_df,
    status = "ok",
    error = NA_character_,
    plot_type = "2d",
    x_var = x_var,
    y_var = y_var
  ))
}

surface_absmax_for_sig_term <- function(fit, sig_row) {
  smooth_obj <- fit$smooth[[sig_row$select_idx]]
  vars <- get_smooth_vars(smooth_obj)
  
  if (length(vars) < 2) {
    return(NA_real_)
  }
  
  panel_df <- make_panel_data_2d(fit, smooth_obj, vars[1], vars[2])
  suppressWarnings(max(abs(panel_df$fit), na.rm = TRUE))
}

compute_2d_colour_limits <- function(selected_models,
                                     scope = c("all_models", "within_model")) {
  scope <- match.arg(scope)
  
  limits_by_key <- list()
  all_abs <- numeric(0)
  
  for (ii in seq_len(nrow(selected_models))) {
    model <- selected_models$model[ii]
    response <- selected_models$response[ii]
    lag_window <- selected_models$lag_window[ii]
    
    message("Pre-scanning 2D colour range for ", model, " / ", lag_window)
    
    fit_file <- find_fit_file(model, lag_window, fit_search_dir)
    fit <- readRDS(fit_file)
    
    sig_terms <- get_significant_smooths(fit, alpha = alpha) |>
      mutate(
        model = model,
        response = response,
        lag_window = lag_window,
        fit_file = fit_file,
        .before = 1
      )
    
    sig_2d <- sig_terms |> filter(term_type == "lag_surface")
    
    vals <- numeric(0)
    
    if (nrow(sig_2d) > 0) {
      vals <- purrr::map_dbl(
        seq_len(nrow(sig_2d)),
        function(jj) {
          tryCatch(
            surface_absmax_for_sig_term(fit, sig_2d[jj, ]),
            error = function(e) NA_real_
          )
        }
      )
      vals <- vals[is.finite(vals)]
    }
    
    key <- paste(model, lag_window, sep = "_")
    
    if (length(vals) > 0) {
      lim <- max(vals, na.rm = TRUE)
      if (!is.finite(lim) || lim == 0) lim <- 1
      limits_by_key[[key]] <- c(-lim, lim)
      all_abs <- c(all_abs, vals)
    } else {
      limits_by_key[[key]] <- NULL
    }
  }
  
  if (scope == "all_models") {
    global_lim <- suppressWarnings(max(all_abs, na.rm = TRUE))
    if (!is.finite(global_lim) || global_lim == 0) global_lim <- 1
    
    for (ii in seq_len(nrow(selected_models))) {
      key <- paste(selected_models$model[ii], selected_models$lag_window[ii], sep = "_")
      limits_by_key[[key]] <- c(-global_lim, global_lim)
    }
  }
  
  limits_by_key
}

empty_plot <- function(title, message) {
  ggplot() +
    annotate("text", x = 0, y = 0, label = message, size = 4.2) +
    xlim(-1, 1) +
    ylim(-1, 1) +
    labs(title = title) +
    theme_void() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
}

# -------------------------------------------------------------------------
# Main execution
# -------------------------------------------------------------------------

selected_models <- if (use_manual_selection) {
  manual_selected_models
} else {
  select_models_from_comparison(step3A_dir)
}

selected_models <- selected_models |>
  mutate(
    model = as.character(model),
    response = as.character(response),
    lag_window = as.character(lag_window)
  )

colour_limits_2d_by_model <- list()

if (use_global_2d_colour_scale) {
  colour_limits_2d_by_model <- compute_2d_colour_limits(
    selected_models,
    scope = global_2d_colour_scope
  )
  
  colour_keys <- paste(selected_models$model, selected_models$lag_window, sep = "_")
  
  lim_table <- tibble(
    model = selected_models$model,
    lag_window = selected_models$lag_window,
    colour_min = purrr::map_dbl(
      colour_keys,
      ~ if (is.null(colour_limits_2d_by_model[[.x]])) NA_real_ else colour_limits_2d_by_model[[.x]][1]
    ),
    colour_max = purrr::map_dbl(
      colour_keys,
      ~ if (is.null(colour_limits_2d_by_model[[.x]])) NA_real_ else colour_limits_2d_by_model[[.x]][2]
    )
  )
  
  write_csv(lim_table, file.path(output_dir, "shared_2d_colour_limits_ggplot_v5.csv"))
}

all_sig_terms <- list()
panel_status <- list()
figure_log <- list()

for (i in seq_len(nrow(selected_models))) {
  model <- selected_models$model[i]
  response <- selected_models$response[i]
  lag_window <- selected_models$lag_window[i]
  
  message("Processing ", model, " / ", lag_window)
  
  fit_file <- find_fit_file(model, lag_window, fit_search_dir)
  fit <- readRDS(fit_file)
  
  sig_terms <- get_significant_smooths(fit, alpha = alpha) |>
    mutate(
      model = model,
      response = response,
      lag_window = lag_window,
      fit_file = fit_file,
      .before = 1
    )
  
  all_sig_terms[[paste(model, lag_window, sep = "_")]] <- sig_terms
  
  safe_model <- clean_file_label(paste(model, lag_window, sep = "_"))
  
  if (nrow(sig_terms) == 0) {
    panels <- list(empty_plot(
      title = pretty_model_label(model, lag_window),
      message = paste0("No significant smooth terms at alpha = ", alpha)
    ))
  } else {
    panels <- vector("list", nrow(sig_terms))
    
    for (j in seq_len(nrow(sig_terms))) {
      sig_row <- sig_terms[j, ]
      
      colour_limits_this_model <- colour_limits_2d_by_model[[safe_model]]
      
      out <- tryCatch(
        make_effect_plot(
          fit,
          sig_row,
          colour_limits_2d = colour_limits_this_model
        ),
        error = function(e) {
          list(
            plot = empty_plot(
              title = pretty_term_label(sig_row$term),
              message = paste0("Plot failed:\n", conditionMessage(e))
            ),
            data = tibble(),
            status = "failed",
            error = conditionMessage(e),
            plot_type = NA_character_,
            x_var = NA_character_,
            y_var = NA_character_
          )
        }
      )
      
      panels[[j]] <- out$plot
      
      panel_status[[length(panel_status) + 1]] <- tibble(
        model = model,
        response = response,
        lag_window = lag_window,
        term = sig_row$term,
        pretty_term = sig_row$pretty_term,
        p_value = sig_row$p_value,
        select_idx = sig_row$select_idx,
        plot_type = out$plot_type,
        x_var = out$x_var,
        y_var = out$y_var,
        status = out$status,
        error = out$error
      )
      
      if (save_individual_panels) {
        indiv_file <- file.path(
          output_dir,
          paste0(
            "panel_",
            safe_model,
            "__",
            sprintf("%02d_", j),
            clean_file_label(sig_row$pretty_term),
            ".png"
          )
        )
        
        ggsave(indiv_file, out$plot, width = 16, height = 5, units = "cm", dpi = png_dpi)
      }
    }
  }
  
  n_col <- min(max_cols, length(panels))
  n_row <- ceiling(length(panels) / n_col)
  
  # Patchwork compatibility note:
  # Some older patchwork/ggplot2 combinations do not support the `& theme(...)`
  # operator for applying a theme to all panels, producing:
  #   Can't find method for generic `&(e1, e2)`.
  # Therefore we avoid `& theme(...)` here and keep layout/theme control inside
  # each panel. The global annotation is added without a theme argument.
  p_combined <- patchwork::wrap_plots(
    panels,
    ncol = n_col,
    guides = "collect"
  ) +
    patchwork::plot_annotation(
      title = pretty_model_label(model, lag_window),
      tag_levels = "A"
    )
  
  fig_height <- max(4.5, panel_height * n_row + 0.7)
  
  png_file <- file.path(output_dir, paste0("figure_effects_", safe_model, "_ggplot_v5.png"))
  pdf_file <- file.path(output_dir, paste0("figure_effects_", safe_model, "_ggplot_v5.pdf"))
  
  ggsave(
    png_file,
    p_combined,
    width = figure_width,
    height = fig_height,
    dpi = png_dpi,
    limitsize = FALSE,
    units="cm"
  )
  
  if (save_pdf) {
    ggsave(
      pdf_file,
      p_combined,
      width = figure_width,
      height = fig_height,
      limitsize = FALSE,
      units="cm"
    )
  } else {
    pdf_file <- NA_character_
  }
  
  figure_log[[length(figure_log) + 1]] <- tibble(
    model = model,
    response = response,
    lag_window = lag_window,
    fit_file = fit_file,
    n_significant_terms = nrow(sig_terms),
    ncol = n_col,
    nrow = n_row,
    png_file = png_file,
    pdf_file = pdf_file
  )
}

sig_out <- bind_rows(all_sig_terms)
panel_status_out <- bind_rows(panel_status)
figure_log_out <- bind_rows(figure_log)

write_csv(sig_out, file.path(output_dir, "significant_terms_selected_models_ggplot_v5.csv"))
write_csv(panel_status_out, file.path(output_dir, "panel_status_ggplot_v5.csv"))
write_csv(figure_log_out, file.path(output_dir, "figure_log_ggplot_v5.csv"))

message("Done.")
message("Outputs written to: ", output_dir)
message("Significant terms table: ", file.path(output_dir, "significant_terms_selected_models_ggplot_v5.csv"))
message("Panel status table: ", file.path(output_dir, "panel_status_ggplot_v5.csv"))
message("Figure log: ", file.path(output_dir, "figure_log_ggplot_v5.csv"))