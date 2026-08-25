# Step 3A - Main ECOSTRESS GAM comparison: no local-hour primary models
#
# Purpose:
#   Compare no-lag and distributed-lag GAMs for:
#     1) reservoir/lake WST   (Tlake  = Tm_l)
#     2) downstream outlet WST (Toutlet = Tm_o)
#     3) outlet-lake WST difference (Tdiff = Tm_o - Tm_l)
#
# Design choice:
#   H_LOC is intentionally excluded from the primary models. This keeps the
#   main inference focused on lagged meteorological forcing + dam properties.
#   H_LOC / sampling-time sensitivity is handled in Step3B.
#
# Required Step 2 outputs in input_dir:
#   dt0_nasa_power_zerolag.csv
#   gam_lag_data_nasa_power.rds
# Optional fallback if the RDS is missing:
#   lag_long_nasa_power_0_48h.csv
#
# Main outputs in output_dir:
#   model_sample_sizes.csv
#   model_comparison_zerolag.csv
#   model_comparison_lag.csv
#   model_comparison_all.csv
#   best_model_by_response.csv
#   lodo_predictions/*.csv
#   lodo_by_dam/*.csv
#   model_summaries/*.csv
#   fits/*.rds

required_packages <- c("mgcv", "dplyr", "purrr", "readr", "tibble", "tidyr", "lubridate")

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
library(purrr)
library(readr)
library(tibble)
library(tidyr)
library(lubridate)

# -------------------------------------------------------------------------
# Local paths
# -------------------------------------------------------------------------

input_dir <- "~/model_inputs"
output_root <- "~/model_outputs"
output_dir <- file.path(output_root, "step3A_main_models")

# Set FALSE if you do not want partial-effect plots. The model summary CSVs
# are always written.
make_partial_effect_plots <- TRUE

# Main lag windows from Step 2.
lag_sets <- list(
  lag_4h = 0:4,
  lag_12h = 0:12,
  lag_24h = 0:24,
  lag_48h = 0:48
)

# Dam-level covariates used in the main models. These are skipped if absent.
# log_res_capacity is derived below from Res_capaci_kmcube.
main_dam_predictors <- c(
  "Dam_hgt",
  "elevation",
  "log_res_area",
  "log_res_capacity"
)

# Minimum observations required for a LODO fold to be fitted.
min_train_rows <- 25L
min_test_rows <- 1L

# -------------------------------------------------------------------------
# Output folders
# -------------------------------------------------------------------------

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
summary_dir <- file.path(output_dir, "model_summaries")
effect_dir <- file.path(output_dir, "significant_partial_effects")
fit_dir <- file.path(output_dir, "fits")
pred_dir <- file.path(output_dir, "lodo_predictions")
by_dam_dir <- file.path(output_dir, "lodo_by_dam")
log_dir <- file.path(output_dir, "logs")

dir.create(summary_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(effect_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fit_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(pred_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(by_dam_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)

plot_log_path <- file.path(log_dir, "partial_effect_plot_log.csv")
lodo_failure_log_path <- file.path(log_dir, "lodo_failure_log.csv")

# -------------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------------

clean_file_label <- function(x) {
  gsub("[^A-Za-z0-9_]+", "_", x)
}

safe_numeric <- function(x) {
  suppressWarnings(as.numeric(x))
}

safe_cor <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if (length(x) < 3 || stats::sd(x) == 0 || stats::sd(y) == 0) {
    return(NA_real_)
  }
  stats::cor(x, y)
}

prediction_metrics <- function(predictions) {
  predictions |>
    filter(is.finite(observed), is.finite(predicted)) |>
    mutate(
      error = predicted - observed,
      abs_error = abs(error),
      sq_error = error^2
    ) |>
    summarise(
      n_dams = n_distinct(dam_id),
      n_observations = n(),
      MAE = mean(abs_error),
      RMSE = sqrt(mean(sq_error)),
      bias = mean(error),
      r = safe_cor(observed, predicted),
      R2_CV = if_else(
        sum((observed - mean(observed))^2) > 0,
        1 - sum(sq_error) / sum((observed - mean(observed))^2),
        NA_real_
      ),
      .groups = "drop"
    )
}

by_dam_metrics <- function(predictions) {
  predictions |>
    filter(is.finite(observed), is.finite(predicted)) |>
    mutate(
      error = predicted - observed,
      abs_error = abs(error),
      sq_error = error^2
    ) |>
    group_by(dam_id) |>
    summarise(
      n = n(),
      MAE = mean(abs_error),
      RMSE = sqrt(mean(sq_error)),
      bias = mean(error),
      r = safe_cor(observed, predicted),
      .groups = "drop"
    )
}

model_metrics <- function(fit) {
  sm <- summary(fit)
  tibble(
    n = length(fit$y),
    AIC = AIC(fit),
    dev_expl = sm$dev.expl,
    adj_r2 = sm$r.sq,
    scale = fit$sig2
  )
}

add_derived_dam_terms <- function(dat) {
  for (cc in c("Dam_hgt", "elevation", "Res_capaci_kmcube", "Res_area_kmsq")) {
    if (cc %in% names(dat)) dat[[cc]] <- safe_numeric(dat[[cc]])
  }
  
  if (!"Res_capaci_kmcube" %in% names(dat)) dat$Res_capaci_kmcube <- NA_real_
  if (!"Res_area_kmsq" %in% names(dat)) dat$Res_area_kmsq <- NA_real_
  
  dat |>
    mutate(
      # Some legacy tables used large raw values for capacity. Keep the earlier
      # practical conversion rule for consistency with previous scripts.
      Res_capacity_km3 = case_when(
        !is.na(Res_capaci_kmcube) & Res_capaci_kmcube > 20 ~ Res_capaci_kmcube / 1000,
        TRUE ~ Res_capaci_kmcube
      ),
      mean_depth_m = if_else(
        !is.na(Res_capacity_km3) & !is.na(Res_area_kmsq) & Res_area_kmsq > 0,
        1000 * Res_capacity_km3 / Res_area_kmsq,
        NA_real_
      ),
      log_res_area = if_else(!is.na(Res_area_kmsq) & Res_area_kmsq >= 0,
                             log1p(Res_area_kmsq), NA_real_),
      log_res_capacity = if_else(!is.na(Res_capacity_km3) & Res_capacity_km3 >= 0,
                                 log1p(Res_capacity_km3), NA_real_)
    )
}

add_basic_fields <- function(dt0) {
  if (!"has_lake" %in% names(dt0)) dt0$has_lake <- !is.na(dt0$Tm_l)
  if (!"has_outlet" %in% names(dt0)) dt0$has_outlet <- !is.na(dt0$Tm_o)
  if (!"paired" %in% names(dt0)) dt0$paired <- dt0$has_lake & dt0$has_outlet
  
  if (!"datetime_utc" %in% names(dt0) && "date.time" %in% names(dt0)) {
    dt0$datetime_utc <- ymd_hms(dt0$date.time, tz = "UTC")
  } else if ("datetime_utc" %in% names(dt0)) {
    dt0$datetime_utc <- ymd_hms(dt0$datetime_utc, tz = "UTC")
  }
  
  if (!"doy" %in% names(dt0) && "datetime_utc" %in% names(dt0)) {
    dt0$doy <- yday(dt0$datetime_utc)
  }
  
  if (!"H_LOC" %in% names(dt0)) {
    if ("Longitude_ref" %in% names(dt0) && "datetime_utc" %in% names(dt0)) {
      H_UTC <- hour(dt0$datetime_utc) + minute(dt0$datetime_utc) / 60 + second(dt0$datetime_utc) / 3600
      dt0$H_LOC <- (H_UTC + safe_numeric(dt0$Longitude_ref) / 15) %% 24
    } else {
      dt0$H_LOC <- NA_real_
    }
  }
  
  dt0 |>
    mutate(
      obs_id = if ("obs_id" %in% names(dt0)) obs_id else row_number(),
      dam = as.character(dam),
      dam_id = factor(dam),
      T_diff = if_else(paired & !is.na(Tm_l) & !is.na(Tm_o), Tm_o - Tm_l, NA_real_)
    ) |>
    add_derived_dam_terms()
}

# -------------------------------------------------------------------------
# Build lag list fallback, if Step 2 RDS is absent
# -------------------------------------------------------------------------

matrix_from_lag_long <- function(lag_long, obs_ids, value_col, lags) {
  lag_cols <- paste0("lag_", lags)
  
  wide <- lag_long |>
    filter(obs_id %in% obs_ids, lag %in% lags) |>
    select(obs_id, lag, value = all_of(value_col)) |>
    mutate(lag = paste0("lag_", lag)) |>
    pivot_wider(names_from = lag, values_from = value)
  
  wide <- tibble(obs_id = obs_ids) |>
    left_join(wide, by = "obs_id") |>
    arrange(match(obs_id, obs_ids))
  
  for (cc in lag_cols) {
    if (!cc %in% names(wide)) wide[[cc]] <- NA_real_
  }
  
  as.matrix(wide[, lag_cols, drop = FALSE])
}

make_mgcv_lag_data <- function(dt0, lag_long, response, lags, require_paired = FALSE) {
  keep <- !is.na(dt0[[response]])
  if (require_paired) keep <- keep & !is.na(dt0$Tm_l) & !is.na(dt0$Tm_o)
  
  base <- dt0[keep, , drop = FALSE] |>
    arrange(obs_id)
  
  obs_ids <- base$obs_id
  out <- as.list(base)
  
  out$lag <- matrix(
    rep(lags, times = nrow(base)),
    nrow = nrow(base),
    byrow = TRUE
  )
  
  out$Air_T_l <- matrix_from_lag_long(lag_long, obs_ids, "Air_T_l", lags)
  out$prec_l <- matrix_from_lag_long(lag_long, obs_ids, "prec_l", lags)
  out$wind_l <- matrix_from_lag_long(lag_long, obs_ids, "wind_l", lags)
  out$pres_l <- matrix_from_lag_long(lag_long, obs_ids, "pres_l", lags)
  out$shortW_l <- matrix_from_lag_long(lag_long, obs_ids, "shortW_l", lags)
  
  out
}

build_gam_lag_data_from_csv <- function(dt0, lag_long_path) {
  message("gam_lag_data_nasa_power.rds not found; rebuilding from lag_long CSV.")
  lag_long <- read_csv(lag_long_path, show_col_types = FALSE)
  
  list(
    lake = map(lag_sets, ~ make_mgcv_lag_data(dt0, lag_long, "Tm_l", lags = .x)),
    outlet = map(lag_sets, ~ make_mgcv_lag_data(dt0, lag_long, "Tm_o", lags = .x)),
    diff = map(
      lag_sets,
      ~ make_mgcv_lag_data(dt0, lag_long, "T_diff", lags = .x, require_paired = TRUE)
    )
  )
}

# -------------------------------------------------------------------------
# Formula builders
# -------------------------------------------------------------------------

smooth_or_linear_terms <- function(dat, vars, k_max = 6) {
  vars <- vars[vars %in% names(dat)]
  out <- character()
  
  for (v in vars) {
    x <- dat[[v]]
    if (is.matrix(x)) x <- as.vector(x)
    x <- safe_numeric(x)
    n_unique <- length(unique(x[is.finite(x)]))
    
    if (n_unique >= 5) {
      k_use <- max(4L, min(k_max, n_unique - 1L))
      out <- c(out, paste0("s(", v, ", bs = 'ts', k = ", k_use, ")"))
    } else if (n_unique >= 2) {
      # Too few unique values for a smooth; retain as a simple parametric term.
      out <- c(out, v)
    }
  }
  
  out
}

cyclic_terms_no_hour <- function(dat) {
  out <- character()
  if ("doy" %in% names(dat)) {
    out <- c(out, "s(doy, bs = 'cc', k = 12)")
  } else if ("m_LOC" %in% names(dat)) {
    out <- c(out, "s(m_LOC, bs = 'cc', k = 12)")
  }
  out
}

meteo_zerolag_terms <- function(dat) {
  smooth_or_linear_terms(dat, c("Air_T", "prec", "wind", "shortW", "pres"), k_max = 6)
}

meteo_lag_terms <- function(lags) {
  k_lag <- max(3L, min(6L, length(lags)))
  c(
    paste0("te(Air_T_l, lag, k = c(6, ", k_lag, "), bs = c('ts', 'ts'))"),
    paste0("te(prec_l, lag, k = c(6, ", k_lag, "), bs = c('ts', 'ts'))"),
    paste0("te(wind_l, lag, k = c(6, ", k_lag, "), bs = c('ts', 'ts'))"),
    paste0("te(shortW_l, lag, k = c(6, ", k_lag, "), bs = c('ts', 'ts'))"),
    paste0("te(pres_l, lag, k = c(6, ", k_lag, "), bs = c('ts', 'ts'))")
  )
}

dam_terms <- function(dat) {
  smooth_or_linear_terms(dat, main_dam_predictors, k_max = 6)
}

make_formula <- function(response, terms, add_re = TRUE) {
  terms <- terms[nzchar(terms)]
  if (add_re) terms <- c(terms, "s(dam_id, bs = 're')")
  
  if (length(terms) == 0) {
    as.formula(paste(response, "~ 1"))
  } else {
    as.formula(paste(response, "~", paste(terms, collapse = " + ")))
  }
}

make_zerolag_formula <- function(response, dat, add_re = TRUE) {
  terms <- c(
    meteo_zerolag_terms(dat),
    cyclic_terms_no_hour(dat),
    dam_terms(dat)
  )
  make_formula(response, terms, add_re = add_re)
}

make_lag_formula <- function(response, lags, dat, add_re = TRUE) {
  terms <- c(
    meteo_lag_terms(lags),
    cyclic_terms_no_hour(dat),
    dam_terms(dat)
  )
  make_formula(response, terms, add_re = add_re)
}

cyclic_knots <- function(dat) {
  out <- list()
  if ("doy" %in% names(dat)) out$doy <- c(0.5, 366.5)
  if (!"doy" %in% names(dat) && "m_LOC" %in% names(dat)) out$m_LOC <- c(0.5, 12.5)
  out
}

complete_rows_df <- function(data, formula) {
  vars <- all.vars(formula)
  vars <- vars[vars %in% names(data)]
  scalar_vars <- vars[!vapply(data[vars], is.matrix, logical(1))]
  
  if (length(scalar_vars) == 0) return(rep(TRUE, nrow(data)))
  stats::complete.cases(data[, scalar_vars, drop = FALSE])
}

complete_rows_list <- function(data, formula) {
  vars <- all.vars(formula)
  vars <- vars[vars %in% names(data)]
  
  response <- vars[1]
  first <- data[[response]]
  n <- if (!is.null(dim(first))) nrow(first) else length(first)
  keep <- rep(TRUE, n)
  
  for (v in vars) {
    x <- data[[v]]
    if (!is.null(dim(x))) {
      if (nrow(x) != n) stop("Variable ", v, " has inconsistent rows.", call. = FALSE)
      if (is.numeric(x) || is.integer(x)) {
        keep <- keep & rowSums(!is.finite(x)) == 0
      } else {
        keep <- keep & rowSums(is.na(x)) == 0
      }
    } else if (length(x) == n) {
      keep <- keep & !is.na(x)
      if (is.numeric(x) || is.integer(x)) keep <- keep & is.finite(x)
    }
  }
  
  keep
}

subset_gam_list <- function(dat, rows) {
  lapply(dat, function(x) {
    if (is.matrix(x)) {
      x[rows, , drop = FALSE]
    } else if (length(x) == length(rows)) {
      x[rows]
    } else {
      x
    }
  })
}

add_dam_terms_to_lag_data <- function(x, dt0) {
  for (nm in names(x)) {
    if (is.null(x[[nm]]$obs_id)) next
    
    obs_ids <- x[[nm]]$obs_id
    add_cols <- dt0 |>
      select(any_of(c(
        "obs_id", "Dam_hgt", "elevation", "Res_capaci_kmcube",
        "Res_capacity_km3", "Res_area_kmsq", "mean_depth_m",
        "log_res_area", "log_res_capacity"
      ))) |>
      filter(obs_id %in% obs_ids) |>
      arrange(match(obs_id, obs_ids))
    
    if (nrow(add_cols) != length(obs_ids)) {
      stop("Could not align dam covariates for ", nm, call. = FALSE)
    }
    
    for (cc in setdiff(names(add_cols), "obs_id")) {
      x[[nm]][[cc]] <- add_cols[[cc]]
    }
  }
  x
}

save_full_model_summary <- function(fit, model_id) {
  sm <- summary(fit)
  safe_id <- clean_file_label(model_id)
  
  global <- tibble(
    model_id = model_id,
    n = length(fit$y),
    AIC = AIC(fit),
    dev_expl = sm$dev.expl,
    adj_r2 = sm$r.sq,
    scale = fit$sig2,
    family = fit$family$family,
    method = fit$method
  )
  
  write_csv(global, file.path(summary_dir, paste0("global_", safe_id, ".csv")))
  
  if (!is.null(sm$p.table)) {
    write_csv(
      as.data.frame(sm$p.table) |>
        rownames_to_column("term") |>
        as_tibble() |>
        mutate(model_id = model_id, .before = 1),
      file.path(summary_dir, paste0("parametric_terms_", safe_id, ".csv"))
    )
  }
  
  if (!is.null(sm$s.table)) {
    smooth <- as.data.frame(sm$s.table) |>
      rownames_to_column("term") |>
      as_tibble() |>
      mutate(model_id = model_id, .before = 1)
    
    write_csv(smooth, file.path(summary_dir, paste0("smooth_terms_", safe_id, ".csv")))
    
    if (isTRUE(make_partial_effect_plots)) {
      save_significant_partial_effects(fit, smooth, model_id)
    }
  }
}

save_significant_partial_effects <- function(fit, smooth_table, model_id, alpha = 0.05) {
  p_col <- grep("p", names(smooth_table), ignore.case = TRUE, value = TRUE)
  p_col <- p_col[length(p_col)]
  if (length(p_col) == 0 || is.na(p_col)) return(invisible(NULL))
  
  sig_terms <- smooth_table |>
    filter(is.finite(.data[[p_col]]), .data[[p_col]] < alpha, !grepl("dam_id", term))
  
  if (nrow(sig_terms) == 0) return(invisible(NULL))
  
  smooth_labels <- vapply(fit$smooth, function(x) x$label, character(1))
  safe_model_id <- clean_file_label(model_id)
  plot_log <- list()
  
  for (j in seq_len(nrow(sig_terms))) {
    term <- sig_terms$term[j]
    select_idx <- match(term, smooth_labels)
    if (is.na(select_idx)) next
    
    out_file <- file.path(
      effect_dir,
      paste0(safe_model_id, "__", clean_file_label(term), ".png")
    )
    
    status <- "saved"
    error_message <- NA_character_
    
    tryCatch(
      {
        png(out_file, width = 1800, height = 1400, res = 180)
        plot(fit, select = select_idx, shade = TRUE, pages = 1, scheme = 2,
             main = paste(model_id, term))
        dev.off()
      },
      error = function(e) {
        status <<- "failed"
        error_message <<- conditionMessage(e)
        if (names(dev.cur()) != "null device") dev.off()
      }
    )
    
    if (names(dev.cur()) != "null device") dev.off()
    
    plot_log[[length(plot_log) + 1]] <- tibble(
      model_id = model_id,
      term = term,
      file = out_file,
      status = status,
      error = error_message
    )
  }
  
  if (length(plot_log) > 0) {
    log_df <- bind_rows(plot_log)
    if (file.exists(plot_log_path)) {
      write_csv(bind_rows(read_csv(plot_log_path, show_col_types = FALSE), log_df), plot_log_path)
    } else {
      write_csv(log_df, plot_log_path)
    }
  }
  
  invisible(NULL)
}

append_lodo_failure <- function(model_id, held_out, message_text) {
  row <- tibble(model_id = model_id, held_out = held_out, message = message_text)
  if (file.exists(lodo_failure_log_path)) {
    write_csv(bind_rows(read_csv(lodo_failure_log_path, show_col_types = FALSE), row), lodo_failure_log_path)
  } else {
    write_csv(row, lodo_failure_log_path)
  }
}

lodo_gam_df <- function(data, formula, response, model_id, dam_col = "dam_id", knots = cyclic_knots(data)) {
  data <- data |>
    filter(!is.na(.data[[response]]))
  
  dams <- unique(as.character(data[[dam_col]]))
  predictions <- vector("list", length(dams))
  
  for (i in seq_along(dams)) {
    held_out <- dams[i]
    train <- data[as.character(data[[dam_col]]) != held_out, , drop = FALSE]
    test <- data[as.character(data[[dam_col]]) == held_out, , drop = FALSE]
    
    if (nrow(train) < min_train_rows || nrow(test) < min_test_rows) {
      append_lodo_failure(model_id, held_out, "Too few train/test rows")
      next
    }
    
    pred <- tryCatch(
      {
        fit <- gam(formula, data = train, family = gaussian(), method = "REML", knots = knots)
        as.numeric(predict(fit, newdata = test, type = "response"))
      },
      error = function(e) {
        append_lodo_failure(model_id, held_out, conditionMessage(e))
        rep(NA_real_, nrow(test))
      }
    )
    
    predictions[[i]] <- tibble(
      obs_id = if ("obs_id" %in% names(test)) test$obs_id else NA_integer_,
      dam_id = held_out,
      observed = test[[response]],
      predicted = pred
    )
  }
  
  bind_rows(predictions)
}

lodo_gam_list <- function(data, formula, response, model_id, dam_col = "dam_id", knots = cyclic_knots(data)) {
  keep <- !is.na(data[[response]])
  data <- subset_gam_list(data, keep)
  
  dams <- unique(as.character(data[[dam_col]]))
  predictions <- vector("list", length(dams))
  
  for (i in seq_along(dams)) {
    held_out <- dams[i]
    train_rows <- as.character(data[[dam_col]]) != held_out
    test_rows <- as.character(data[[dam_col]]) == held_out
    
    train <- subset_gam_list(data, train_rows)
    test <- subset_gam_list(data, test_rows)
    
    n_train <- sum(train_rows)
    n_test <- sum(test_rows)
    
    if (n_train < min_train_rows || n_test < min_test_rows) {
      append_lodo_failure(model_id, held_out, "Too few train/test rows")
      next
    }
    
    pred <- tryCatch(
      {
        fit <- gam(formula, data = train, family = gaussian(), method = "REML", knots = knots)
        as.numeric(predict(fit, newdata = test, type = "response"))
      },
      error = function(e) {
        append_lodo_failure(model_id, held_out, conditionMessage(e))
        rep(NA_real_, n_test)
      }
    )
    
    predictions[[i]] <- tibble(
      obs_id = if (!is.null(test$obs_id)) test$obs_id else NA_integer_,
      dam_id = held_out,
      observed = test[[response]],
      predicted = pred
    )
  }
  
  bind_rows(predictions)
}

# -------------------------------------------------------------------------
# Read inputs
# -------------------------------------------------------------------------

dt0_path <- file.path(input_dir, "dt0_nasa_power_zerolag.csv")
gam_lag_path <- file.path(input_dir, "gam_lag_data_nasa_power.rds")
lag_long_path <- file.path(input_dir, "lag_long_nasa_power_0_48h.csv")

if (!file.exists(dt0_path)) stop("Missing Step 2 output: ", dt0_path, call. = FALSE)

dt0 <- read_csv(dt0_path, show_col_types = FALSE) |>
  select(-matches("^\\.\\.\\.|^Unnamed")) |>
  add_basic_fields()

if (file.exists(gam_lag_path)) {
  gam_lag_data <- readRDS(gam_lag_path)
} else if (file.exists(lag_long_path)) {
  gam_lag_data <- build_gam_lag_data_from_csv(dt0, lag_long_path)
  saveRDS(gam_lag_data, file.path(output_dir, "gam_lag_data_rebuilt_from_lag_long.rds"))
} else {
  stop("Missing both gam_lag_data_nasa_power.rds and lag_long_nasa_power_0_48h.csv.", call. = FALSE)
}

gam_lag_data$lake <- add_dam_terms_to_lag_data(gam_lag_data$lake, dt0)
gam_lag_data$outlet <- add_dam_terms_to_lag_data(gam_lag_data$outlet, dt0)
gam_lag_data$diff <- add_dam_terms_to_lag_data(gam_lag_data$diff, dt0)

available_dam_predictors <- main_dam_predictors[main_dam_predictors %in% names(dt0)]
message("Dam predictors used if complete: ", paste(available_dam_predictors, collapse = ", "))

model_sample_sizes <- tibble(
  model = c("Tlake", "Toutlet", "Tdiff"),
  response = c("Tm_l", "Tm_o", "T_diff"),
  n = c(sum(!is.na(dt0$Tm_l)), sum(!is.na(dt0$Tm_o)), sum(!is.na(dt0$T_diff))),
  n_dams = c(
    n_distinct(dt0$dam_id[!is.na(dt0$Tm_l)]),
    n_distinct(dt0$dam_id[!is.na(dt0$Tm_o)]),
    n_distinct(dt0$dam_id[!is.na(dt0$T_diff)])
  )
)
write_csv(model_sample_sizes, file.path(output_dir, "model_sample_sizes.csv"))
print(model_sample_sizes)

# -------------------------------------------------------------------------
# Fit zero-lag models
# -------------------------------------------------------------------------

zero_specs <- tribble(
  ~model, ~response,
  "Tlake", "Tm_l",
  "Toutlet", "Tm_o",
  "Tdiff", "T_diff"
)

zero_results <- vector("list", nrow(zero_specs))

for (i in seq_len(nrow(zero_specs))) {
  model_name <- zero_specs$model[i]
  response <- zero_specs$response[i]
  lag_window <- "0h"
  model_id <- paste(model_name, lag_window, sep = "_")
  
  message("Fitting zero-lag model: ", model_id)
  
  dat <- dt0 |>
    filter(!is.na(.data[[response]]))
  
  form <- make_zerolag_formula(response, dat, add_re = TRUE)
  form_cv <- make_zerolag_formula(response, dat, add_re = FALSE)
  dat <- dat[complete_rows_df(dat, form), , drop = FALSE]
  
  fit <- gam(form, data = dat, family = gaussian(), method = "REML", knots = cyclic_knots(dat))
  save_full_model_summary(fit, model_id)
  saveRDS(fit, file.path(fit_dir, paste0("fit_", model_id, ".rds")))
  
  preds <- lodo_gam_df(data = dat, formula = form_cv, response = response, model_id = model_id, knots = cyclic_knots(dat))
  pred_file <- file.path(pred_dir, paste0("lodo_predictions_", model_id, ".csv"))
  by_dam_file <- file.path(by_dam_dir, paste0("lodo_by_dam_", model_id, ".csv"))
  write_csv(preds, pred_file)
  write_csv(by_dam_metrics(preds), by_dam_file)
  
  zero_results[[i]] <- model_metrics(fit) |>
    bind_cols(prediction_metrics(preds)) |>
    mutate(
      model = model_name,
      response = response,
      lag_window = lag_window,
      model_type = "zerolag",
      formula = paste(deparse(form), collapse = " "),
      cv_formula = paste(deparse(form_cv), collapse = " "),
      prediction_file = pred_file,
      by_dam_file = by_dam_file,
      fit_file = file.path(fit_dir, paste0("fit_", model_id, ".rds"))
    ) |>
    select(model, response, lag_window, model_type, everything())
}

zero_comparison <- bind_rows(zero_results)
write_csv(zero_comparison, file.path(output_dir, "model_comparison_zerolag.csv"))

# -------------------------------------------------------------------------
# Fit distributed-lag models
# -------------------------------------------------------------------------

lag_specs <- tribble(
  ~model, ~response, ~data_name,
  "Tlake", "Tm_l", "lake",
  "Toutlet", "Tm_o", "outlet",
  "Tdiff", "T_diff", "diff"
)

lag_results <- list()

for (i in seq_len(nrow(lag_specs))) {
  model_name <- lag_specs$model[i]
  response <- lag_specs$response[i]
  data_name <- lag_specs$data_name[i]
  
  for (lag_window in names(gam_lag_data[[data_name]])) {
    model_id <- paste(model_name, lag_window, sep = "_")
    message("Fitting lag model: ", model_id)
    
    dat <- gam_lag_data[[data_name]][[lag_window]]
    lags <- sort(unique(as.vector(dat$lag)))
    form <- make_lag_formula(response, lags, dat, add_re = TRUE)
    form_cv <- make_lag_formula(response, lags, dat, add_re = FALSE)
    dat <- subset_gam_list(dat, complete_rows_list(dat, form))
    
    fit <- gam(form, data = dat, family = gaussian(), method = "REML", knots = cyclic_knots(dat))
    save_full_model_summary(fit, model_id)
    saveRDS(fit, file.path(fit_dir, paste0("fit_", model_id, ".rds")))
    
    preds <- lodo_gam_list(data = dat, formula = form_cv, response = response, model_id = model_id, knots = cyclic_knots(dat))
    pred_file <- file.path(pred_dir, paste0("lodo_predictions_", model_id, ".csv"))
    by_dam_file <- file.path(by_dam_dir, paste0("lodo_by_dam_", model_id, ".csv"))
    write_csv(preds, pred_file)
    write_csv(by_dam_metrics(preds), by_dam_file)
    
    lag_results[[model_id]] <- model_metrics(fit) |>
      bind_cols(prediction_metrics(preds)) |>
      mutate(
        model = model_name,
        response = response,
        lag_window = lag_window,
        model_type = "distributed_lag",
        formula = paste(deparse(form), collapse = " "),
        cv_formula = paste(deparse(form_cv), collapse = " "),
        prediction_file = pred_file,
        by_dam_file = by_dam_file,
        fit_file = file.path(fit_dir, paste0("fit_", model_id, ".rds"))
      ) |>
      select(model, response, lag_window, model_type, everything())
  }
}

lag_comparison <- bind_rows(lag_results)
write_csv(lag_comparison, file.path(output_dir, "model_comparison_lag.csv"))

model_comparison_all <- bind_rows(zero_comparison, lag_comparison) |>
  arrange(model, RMSE, AIC)

write_csv(model_comparison_all, file.path(output_dir, "model_comparison_all.csv"))

best_model_by_response <- model_comparison_all |>
  filter(is.finite(RMSE)) |>
  group_by(model, response) |>
  arrange(RMSE, AIC, .by_group = TRUE) |>
  slice(1) |>
  ungroup()

write_csv(best_model_by_response, file.path(output_dir, "best_model_by_response.csv"))

message("Done Step3A main model comparison.")
message("Outputs written to: ", output_dir)
message("Best models by LODO RMSE:")
print(best_model_by_response |> select(model, response, lag_window, n_observations, MAE, RMSE, bias, r, R2_CV, AIC, dev_expl))
