# Step 3B - ECOSTRESS sampling coverage + local-hour / sampling sensitivity
#
# Purpose:
#   Supporting-material diagnostics for the reviewer concern that ECOSTRESS
#   observations are sparse and unevenly distributed in dam, month/day and hour.
#
# What this script does:
#   1) Summarises observation availability by dam, month, day-of-year and local hour.
#   2) Calculates evenness / concentration metrics for sampling distribution.
#   3) Compares selected Step3A models with and without H_LOC/cyclic terms.
#   4) Runs optional sampling sensitivity for selected models:
#        - full unweighted
#        - full weighted by dam-month-hour sampling density
#        - balanced one-observation-per-dam-month-hour-bin subset
#
# Required Step 2 outputs:
#   dt0_nasa_power_zerolag.csv
#   gam_lag_data_nasa_power.rds
#
# Required Step 3A output, if available:
#   model_outputs/step3A_main_models/best_model_by_response.csv
# If absent, default lag windows are used.

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

has_ggplot2 <- requireNamespace("ggplot2", quietly = TRUE)

# -------------------------------------------------------------------------
# Local paths and switches
# -------------------------------------------------------------------------

input_dir <- "~/model_inputs"
output_root <- "~/model_outputs"
step3A_dir <- file.path(output_root, "step3A_main_models")
output_dir <- file.path(output_root, "step3B_sampling_coverage_hour_sensitivity")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_dir, "coverage_tables"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_dir, "plots"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_dir, "hour_structure_sensitivity"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_dir, "sampling_sensitivity"), showWarnings = FALSE, recursive = TRUE)

run_hour_structure_sensitivity <- TRUE
run_sampling_sensitivity <- TRUE

# Structures for the local-hour/cyclic sensitivity.
# no_hour is the main-inference structure used in Step3A.
structure_set <- c("no_hour", "with_hour", "no_cyclic", "cyclic_only")

main_dam_predictors <- c(
  "Dam_hgt",
  "elevation",
  "log_res_area",
  "log_res_capacity"
)

min_train_rows <- 25L

# -------------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------------

safe_numeric <- function(x) suppressWarnings(as.numeric(x))

safe_cor <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if (length(x) < 3 || stats::sd(x) == 0 || stats::sd(y) == 0) return(NA_real_)
  stats::cor(x, y)
}

clean_file_label <- function(x) gsub("[^A-Za-z0-9_]+", "_", x)

calc_gini <- function(x) {
  x <- x[is.finite(x) & x >= 0]
  if (length(x) == 0 || sum(x) == 0) return(NA_real_)
  x <- sort(x)
  n <- length(x)
  (2 * sum(seq_len(n) * x) / (n * sum(x))) - (n + 1) / n
}

calc_shannon_evenness <- function(x) {
  x <- x[is.finite(x) & x > 0]
  if (length(x) <= 1 || sum(x) == 0) return(NA_real_)
  p <- x / sum(x)
  -sum(p * log(p)) / log(length(x))
}

prediction_metrics <- function(predictions) {
  predictions |>
    filter(is.finite(observed), is.finite(predicted)) |>
    mutate(error = predicted - observed, abs_error = abs(error), sq_error = error^2) |>
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

add_derived_dam_terms <- function(dat) {
  for (cc in c("Dam_hgt", "elevation", "Res_capaci_kmcube", "Res_area_kmsq")) {
    if (cc %in% names(dat)) dat[[cc]] <- safe_numeric(dat[[cc]])
  }
  
  if (!"Res_capaci_kmcube" %in% names(dat)) dat$Res_capaci_kmcube <- NA_real_
  if (!"Res_area_kmsq" %in% names(dat)) dat$Res_area_kmsq <- NA_real_
  
  dat |>
    mutate(
      Res_capacity_km3 = case_when(
        !is.na(Res_capaci_kmcube) & Res_capaci_kmcube > 20 ~ Res_capaci_kmcube / 1000,
        TRUE ~ Res_capaci_kmcube
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
  
  if (!"doy" %in% names(dt0) && "datetime_utc" %in% names(dt0)) dt0$doy <- yday(dt0$datetime_utc)
  if (!"m_LOC" %in% names(dt0) && "datetime_utc" %in% names(dt0)) dt0$m_LOC <- month(dt0$datetime_utc)
  
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
      T_diff = if_else(paired & !is.na(Tm_l) & !is.na(Tm_o), Tm_o - Tm_l, NA_real_),
      obs_domain_clean = case_when(
        !is.na(obs_domain) ~ as.character(obs_domain),
        paired ~ "paired",
        has_lake & !has_outlet ~ "lake_only",
        has_outlet & !has_lake ~ "outlet_only",
        TRUE ~ "missing_domain"
      ),
      hour_bin = if_else(is.finite(H_LOC), floor(H_LOC) %% 24, NA_real_),
      month_bin = if_else(!is.na(m_LOC), as.integer(m_LOC), as.integer(month(datetime_utc))),
      doy_bin = if_else(!is.na(doy), as.integer(doy), as.integer(yday(datetime_utc)))
    ) |>
    add_derived_dam_terms()
}

# -------------------------------------------------------------------------
# Formula helpers for sensitivity models
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
      out <- c(out, v)
    }
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

cyclic_terms_for_structure <- function(dat, structure) {
  out <- character()
  if (structure %in% c("with_hour", "cyclic_only") && "H_LOC" %in% names(dat)) {
    out <- c(out, "s(H_LOC, bs = 'cc', k = 8)")
  }
  if (structure %in% c("no_hour", "with_hour", "cyclic_only")) {
    if ("doy" %in% names(dat)) {
      out <- c(out, "s(doy, bs = 'cc', k = 12)")
    } else if ("m_LOC" %in% names(dat)) {
      out <- c(out, "s(m_LOC, bs = 'cc', k = 12)")
    }
  }
  out
}

dam_terms <- function(dat) {
  smooth_or_linear_terms(dat, main_dam_predictors, k_max = 6)
}

make_sensitivity_formula <- function(response, dat, lag_window, structure) {
  terms <- character()
  
  if (structure != "cyclic_only") {
    if (lag_window == "0h") {
      terms <- c(terms, meteo_zerolag_terms(dat))
    } else {
      lags <- sort(unique(as.vector(dat$lag)))
      terms <- c(terms, meteo_lag_terms(lags))
    }
  }
  
  terms <- c(
    terms,
    cyclic_terms_for_structure(dat, structure),
    dam_terms(dat)
  )
  
  terms <- terms[nzchar(terms)]
  if (length(terms) == 0) as.formula(paste(response, "~ 1")) else as.formula(paste(response, "~", paste(terms, collapse = " + ")))
}

cyclic_knots <- function(dat, structure = "no_hour") {
  out <- list()
  if (structure %in% c("with_hour", "cyclic_only") && "H_LOC" %in% names(dat)) out$H_LOC <- c(0, 24)
  if (structure %in% c("no_hour", "with_hour", "cyclic_only")) {
    if ("doy" %in% names(dat)) out$doy <- c(0.5, 366.5)
    if (!"doy" %in% names(dat) && "m_LOC" %in% names(dat)) out$m_LOC <- c(0.5, 12.5)
  }
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
  n <- if (!is.null(dim(data[[response]]))) nrow(data[[response]]) else length(data[[response]])
  keep <- rep(TRUE, n)
  
  for (v in vars) {
    x <- data[[v]]
    if (!is.null(dim(x))) {
      keep <- keep & rowSums(!is.finite(x)) == 0
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

add_cols_to_lag_data <- function(x, dt0, cols) {
  for (nm in names(x)) {
    if (is.null(x[[nm]]$obs_id)) next
    obs_ids <- x[[nm]]$obs_id
    add_cols <- dt0 |>
      select(any_of(c("obs_id", cols))) |>
      filter(obs_id %in% obs_ids) |>
      arrange(match(obs_id, obs_ids))
    if (nrow(add_cols) != length(obs_ids)) stop("Could not align columns for ", nm, call. = FALSE)
    for (cc in setdiff(names(add_cols), "obs_id")) x[[nm]][[cc]] <- add_cols[[cc]]
  }
  x
}

lodo_df <- function(data, formula, response, knots, weights_col = NULL) {
  data <- data |>
    filter(!is.na(.data[[response]]))
  dams <- unique(as.character(data$dam_id))
  out <- vector("list", length(dams))
  
  for (i in seq_along(dams)) {
    held <- dams[i]
    train <- data[as.character(data$dam_id) != held, , drop = FALSE]
    test <- data[as.character(data$dam_id) == held, , drop = FALSE]
    if (nrow(train) < min_train_rows || nrow(test) == 0) next
    
    pred <- tryCatch({
      if (!is.null(weights_col) && weights_col %in% names(train)) {
        fit <- gam(formula, data = train, family = gaussian(), method = "REML", knots = knots,
                   weights = train[[weights_col]])
      } else {
        fit <- gam(formula, data = train, family = gaussian(), method = "REML", knots = knots)
      }
      as.numeric(predict(fit, newdata = test, type = "response"))
    }, error = function(e) rep(NA_real_, nrow(test)))
    
    out[[i]] <- tibble(obs_id = test$obs_id, dam_id = held, observed = test[[response]], predicted = pred)
  }
  
  bind_rows(out)
}

lodo_list <- function(data, formula, response, knots, weights_col = NULL) {
  keep <- !is.na(data[[response]])
  data <- subset_gam_list(data, keep)
  dams <- unique(as.character(data$dam_id))
  out <- vector("list", length(dams))
  
  for (i in seq_along(dams)) {
    held <- dams[i]
    train_rows <- as.character(data$dam_id) != held
    test_rows <- as.character(data$dam_id) == held
    train <- subset_gam_list(data, train_rows)
    test <- subset_gam_list(data, test_rows)
    n_test <- sum(test_rows)
    if (sum(train_rows) < min_train_rows || n_test == 0) next
    
    pred <- tryCatch({
      if (!is.null(weights_col) && !is.null(train[[weights_col]])) {
        fit <- gam(formula, data = train, family = gaussian(), method = "REML", knots = knots,
                   weights = train[[weights_col]])
      } else {
        fit <- gam(formula, data = train, family = gaussian(), method = "REML", knots = knots)
      }
      as.numeric(predict(fit, newdata = test, type = "response"))
    }, error = function(e) rep(NA_real_, n_test))
    
    out[[i]] <- tibble(obs_id = test$obs_id, dam_id = held, observed = test[[response]], predicted = pred)
  }
  
  bind_rows(out)
}

get_response_data <- function(model_name, response, data_name, lag_window, dt0, gam_lag_data) {
  if (lag_window == "0h") {
    dt0 |>
      filter(!is.na(.data[[response]]))
  } else {
    gam_lag_data[[data_name]][[lag_window]]
  }
}

run_lodo_model <- function(dat, response, lag_window, structure, weights_col = NULL) {
  form <- make_sensitivity_formula(response, dat, lag_window, structure)
  knots <- cyclic_knots(dat, structure)
  
  if (lag_window == "0h") {
    dat <- dat[complete_rows_df(dat, form), , drop = FALSE]
    lodo_df(dat, form, response, knots, weights_col = weights_col)
  } else {
    dat <- subset_gam_list(dat, complete_rows_list(dat, form))
    lodo_list(dat, form, response, knots, weights_col = weights_col)
  }
}

# -------------------------------------------------------------------------
# Read inputs
# -------------------------------------------------------------------------

dt0_path <- file.path(input_dir, "dt0_nasa_power_zerolag.csv")
gam_lag_path <- file.path(input_dir, "gam_lag_data_nasa_power.rds")

if (!file.exists(dt0_path)) stop("Missing: ", dt0_path, call. = FALSE)
if (!file.exists(gam_lag_path)) stop("Missing: ", gam_lag_path, call. = FALSE)

dt0 <- read_csv(dt0_path, show_col_types = FALSE) |>
  select(-matches("^\\.\\.\\.|^Unnamed")) |>
  add_basic_fields()

gam_lag_data <- readRDS(gam_lag_path)
cols_to_add <- c(
  "H_LOC", "doy", "m_LOC", "hour_bin", "month_bin", "doy_bin",
  "Dam_hgt", "elevation", "Res_capaci_kmcube", "Res_capacity_km3",
  "Res_area_kmsq", "log_res_area", "log_res_capacity",
  "observation_source", "obs_domain_clean", "filled_from_additional"
)

gam_lag_data$lake <- add_cols_to_lag_data(gam_lag_data$lake, dt0, cols_to_add)
gam_lag_data$outlet <- add_cols_to_lag_data(gam_lag_data$outlet, dt0, cols_to_add)
gam_lag_data$diff <- add_cols_to_lag_data(gam_lag_data$diff, dt0, cols_to_add)

# -------------------------------------------------------------------------
# Coverage tables
# -------------------------------------------------------------------------

coverage_dir <- file.path(output_dir, "coverage_tables")

coverage_overall <- dt0 |>
  summarise(
    n_rows = n(),
    n_dams = n_distinct(dam),
    n_lake = sum(!is.na(Tm_l)),
    n_outlet = sum(!is.na(Tm_o)),
    n_paired = sum(!is.na(T_diff)),
    min_datetime = min(datetime_utc, na.rm = TRUE),
    max_datetime = max(datetime_utc, na.rm = TRUE),
    .groups = "drop"
  )
write_csv(coverage_overall, file.path(coverage_dir, "coverage_overall.csv"))

coverage_by_domain <- dt0 |>
  count(obs_domain_clean, observation_source, name = "n") |>
  arrange(desc(n))
write_csv(coverage_by_domain, file.path(coverage_dir, "coverage_by_domain_and_source.csv"))

coverage_by_dam <- dt0 |>
  group_by(dam, obs_domain_clean) |>
  summarise(
    n = n(),
    n_lake = sum(!is.na(Tm_l)),
    n_outlet = sum(!is.na(Tm_o)),
    n_paired = sum(!is.na(T_diff)),
    .groups = "drop"
  )
write_csv(coverage_by_dam, file.path(coverage_dir, "coverage_by_dam.csv"))

coverage_by_month <- dt0 |>
  group_by(month_bin, obs_domain_clean) |>
  summarise(n = n(), n_dams = n_distinct(dam), .groups = "drop")
write_csv(coverage_by_month, file.path(coverage_dir, "coverage_by_month.csv"))

coverage_by_hour <- dt0 |>
  group_by(hour_bin, obs_domain_clean) |>
  summarise(n = n(), n_dams = n_distinct(dam), .groups = "drop")
write_csv(coverage_by_hour, file.path(coverage_dir, "coverage_by_local_hour.csv"))

coverage_by_doy <- dt0 |>
  group_by(doy_bin, obs_domain_clean) |>
  summarise(n = n(), n_dams = n_distinct(dam), .groups = "drop")
write_csv(coverage_by_doy, file.path(coverage_dir, "coverage_by_day_of_year.csv"))

coverage_month_hour <- dt0 |>
  group_by(month_bin, hour_bin, obs_domain_clean) |>
  summarise(n = n(), n_dams = n_distinct(dam), .groups = "drop")
write_csv(coverage_month_hour, file.path(coverage_dir, "coverage_by_month_local_hour.csv"))

coverage_dam_month_hour <- dt0 |>
  group_by(dam, month_bin, hour_bin, obs_domain_clean) |>
  summarise(n = n(), .groups = "drop")
write_csv(coverage_dam_month_hour, file.path(coverage_dir, "coverage_by_dam_month_local_hour.csv"))

make_evenness_row <- function(dat, label, group_vars, n_possible = NA_real_) {
  counts <- dat |>
    count(across(all_of(group_vars)), name = "n")
  
  tibble(
    distribution = label,
    n_observations = sum(counts$n),
    n_occupied_bins = nrow(counts),
    n_possible_bins = n_possible,
    coverage_fraction = if_else(is.finite(n_possible) && n_possible > 0, nrow(counts) / n_possible, NA_real_),
    shannon_evenness = calc_shannon_evenness(counts$n),
    gini = calc_gini(counts$n),
    max_bin_share = max(counts$n) / sum(counts$n)
  )
}

n_dams <- n_distinct(dt0$dam)
evenness_all <- bind_rows(
  make_evenness_row(dt0, "dam", c("dam"), n_dams),
  make_evenness_row(dt0, "month", c("month_bin"), 12),
  make_evenness_row(dt0, "local_hour", c("hour_bin"), 24),
  make_evenness_row(dt0, "day_of_year", c("doy_bin"), 366),
  make_evenness_row(dt0, "month_local_hour", c("month_bin", "hour_bin"), 12 * 24),
  make_evenness_row(dt0, "dam_month_local_hour", c("dam", "month_bin", "hour_bin"), n_dams * 12 * 24)
)
write_csv(evenness_all, file.path(coverage_dir, "coverage_evenness_all_observations.csv"))

evenness_by_response <- bind_rows(
  make_evenness_row(filter(dt0, !is.na(Tm_l)), "Tlake_dam_month_hour", c("dam", "month_bin", "hour_bin"), n_dams * 12 * 24) |> mutate(model = "Tlake"),
  make_evenness_row(filter(dt0, !is.na(Tm_o)), "Toutlet_dam_month_hour", c("dam", "month_bin", "hour_bin"), n_dams * 12 * 24) |> mutate(model = "Toutlet"),
  make_evenness_row(filter(dt0, !is.na(T_diff)), "Tdiff_dam_month_hour", c("dam", "month_bin", "hour_bin"), n_dams * 12 * 24) |> mutate(model = "Tdiff")
) |>
  select(model, everything())
write_csv(evenness_by_response, file.path(coverage_dir, "coverage_evenness_by_model_dataset.csv"))

# Optional quick plots.
if (has_ggplot2) {
  p1 <- ggplot2::ggplot(coverage_by_hour, ggplot2::aes(x = hour_bin, y = n)) +
    ggplot2::geom_col() +
    ggplot2::facet_wrap(~obs_domain_clean, scales = "free_y") +
    ggplot2::labs(x = "Local solar hour", y = "Number of observations")
  ggplot2::ggsave(file.path(output_dir, "plots", "coverage_by_local_hour.png"), p1, width = 9, height = 5, dpi = 200)
  
  p2 <- ggplot2::ggplot(coverage_month_hour, ggplot2::aes(x = hour_bin, y = month_bin, fill = n)) +
    ggplot2::geom_tile() +
    ggplot2::facet_wrap(~obs_domain_clean) +
    ggplot2::labs(x = "Local solar hour", y = "Month", fill = "n")
  ggplot2::ggsave(file.path(output_dir, "plots", "coverage_month_by_hour_heatmap.png"), p2, width = 10, height = 5, dpi = 200)
}
p1
p2
# -------------------------------------------------------------------------
# Determine selected models from Step3A
# -------------------------------------------------------------------------

best_path <- file.path(step3A_dir, "best_model_by_response.csv")

if (file.exists(best_path)) {
  selected_models <- read_csv(best_path, show_col_types = FALSE) |>
    select(model, response, lag_window)
} else {
  warning("Step3A best_model_by_response.csv not found; using default lag windows.")
  selected_models <- tribble(
    ~model, ~response, ~lag_window,
    "Tlake", "Tm_l", "lag_24h",
    "Toutlet", "Tm_o", "lag_4h",
    "Tdiff", "T_diff", "lag_24h"
  )
}

selected_models <- selected_models |>
  mutate(
    data_name = case_when(
      model == "Tlake" ~ "lake",
      model == "Toutlet" ~ "outlet",
      model == "Tdiff" ~ "diff",
      TRUE ~ NA_character_
    )
  )

write_csv(selected_models, file.path(output_dir, "selected_models_for_sensitivity.csv"))

# -------------------------------------------------------------------------
# Hour/cyclic structure sensitivity
# -------------------------------------------------------------------------

if (run_hour_structure_sensitivity) {
  message("Running local-hour/cyclic structure sensitivity...")
  hour_results <- list()
  sens_dir <- file.path(output_dir, "hour_structure_sensitivity")
  
  for (i in seq_len(nrow(selected_models))) {
    model_name <- selected_models$model[i]
    response <- selected_models$response[i]
    data_name <- selected_models$data_name[i]
    lag_window <- selected_models$lag_window[i]
    
    for (structure in structure_set) {
      dat <- get_response_data(model_name, response, data_name, lag_window, dt0, gam_lag_data)
      preds <- run_lodo_model(dat, response, lag_window, structure)
      
      pred_file <- file.path(
        sens_dir,
        paste0("lodo_predictions_", model_name, "_", lag_window, "_", structure, ".csv")
      )
      write_csv(preds, pred_file)
      
      hour_results[[paste(model_name, lag_window, structure, sep = "_")]] <- prediction_metrics(preds) |>
        mutate(
          model = model_name,
          response = response,
          selected_lag_window = lag_window,
          structure = structure,
          prediction_file = pred_file
        ) |>
        select(model, response, selected_lag_window, structure, everything())
    }
  }
  
  hour_sensitivity <- bind_rows(hour_results) |>
    arrange(model, RMSE)
  
  write_csv(hour_sensitivity, file.path(sens_dir, "model_comparison_hour_cyclic_structure_sensitivity.csv"))
}

# -------------------------------------------------------------------------
# Sampling density sensitivity: full, weighted, balanced
# -------------------------------------------------------------------------

make_response_weights <- function(dt0, response) {
  dat <- dt0 |>
    filter(!is.na(.data[[response]])) |>
    mutate(
      sampling_bin = paste(dam_id, month_bin, hour_bin, sep = "__")
    ) |>
    group_by(sampling_bin) |>
    mutate(bin_n = n()) |>
    ungroup() |>
    mutate(
      sampling_weight = 1 / bin_n,
      sampling_weight = sampling_weight / mean(sampling_weight, na.rm = TRUE)
    ) |>
    select(obs_id, sampling_weight, bin_n)
  
  dt0 |>
    select(-any_of(c("sampling_weight", "bin_n"))) |>
    left_join(dat, by = "obs_id") |>
    mutate(sampling_weight = if_else(is.na(sampling_weight), 1, sampling_weight))
}

select_balanced_ids <- function(dt0, model_name, response) {
  support_col <- case_when(
    model_name == "Tlake" ~ "Npixel_l",
    model_name == "Toutlet" ~ "Npixel_o",
    model_name == "Tdiff" ~ "Npixel_pair",
    TRUE ~ "Npixel_l"
  )
  
  tmp <- dt0 |>
    mutate(
      Npixel_pair = pmin(
        if_else(is.na(Npixel_l), 0, Npixel_l),
        if_else(is.na(Npixel_o), 0, Npixel_o)
      ),
      support_value = if (support_col %in% names(dt0)) safe_numeric(.data[[support_col]]) else 1
    ) |>
    filter(!is.na(.data[[response]]), !is.na(dam_id), !is.na(month_bin), !is.na(hour_bin)) |>
    arrange(dam_id, month_bin, hour_bin, desc(support_value), obs_id) |>
    group_by(dam_id, month_bin, hour_bin) |>
    slice(1) |>
    ungroup()
  
  tmp$obs_id
}

subset_response_data_by_ids <- function(dat, obs_ids, lag_window) {
  if (lag_window == "0h") {
    dat |>
      filter(obs_id %in% obs_ids)
  } else {
    subset_gam_list(dat, dat$obs_id %in% obs_ids)
  }
}

add_weights_to_response_data <- function(dat, dt0_weighted, lag_window) {
  wtab <- dt0_weighted |>
    select(obs_id, sampling_weight, bin_n)
  
  if (lag_window == "0h") {
    dat |>
      select(-any_of(c("sampling_weight", "bin_n"))) |>
      left_join(wtab, by = "obs_id") |>
      mutate(sampling_weight = if_else(is.na(sampling_weight), 1, sampling_weight))
  } else {
    idx <- match(dat$obs_id, wtab$obs_id)
    dat$sampling_weight <- wtab$sampling_weight[idx]
    dat$bin_n <- wtab$bin_n[idx]
    dat$sampling_weight[is.na(dat$sampling_weight)] <- 1
    dat
  }
}

if (run_sampling_sensitivity) {
  message("Running sampling density sensitivity...")
  sampling_results <- list()
  sampling_dir <- file.path(output_dir, "sampling_sensitivity")
  
  for (i in seq_len(nrow(selected_models))) {
    model_name <- selected_models$model[i]
    response <- selected_models$response[i]
    data_name <- selected_models$data_name[i]
    lag_window <- selected_models$lag_window[i]
    structure <- "no_hour"
    
    base_dat <- get_response_data(model_name, response, data_name, lag_window, dt0, gam_lag_data)
    
    # 1. Full unweighted.
    preds_full <- run_lodo_model(base_dat, response, lag_window, structure)
    write_csv(preds_full, file.path(sampling_dir, paste0("lodo_", model_name, "_", lag_window, "_full_unweighted.csv")))
    sampling_results[[paste(model_name, "full_unweighted", sep = "_")]] <- prediction_metrics(preds_full) |>
      mutate(model = model_name, response = response, lag_window = lag_window, sampling_mode = "full_unweighted")
    
    # 2. Full weighted by dam-month-hour density.
    dt0_weighted <- make_response_weights(dt0, response)
    weighted_dat <- add_weights_to_response_data(base_dat, dt0_weighted, lag_window)
    preds_weighted <- run_lodo_model(weighted_dat, response, lag_window, structure, weights_col = "sampling_weight")
    write_csv(preds_weighted, file.path(sampling_dir, paste0("lodo_", model_name, "_", lag_window, "_full_weighted.csv")))
    sampling_results[[paste(model_name, "full_weighted", sep = "_")]] <- prediction_metrics(preds_weighted) |>
      mutate(model = model_name, response = response, lag_window = lag_window, sampling_mode = "full_weighted")
    
    # 3. Balanced subset: one observation per dam-month-hour bin.
    balanced_ids <- select_balanced_ids(dt0, model_name, response)
    balanced_dat <- subset_response_data_by_ids(base_dat, balanced_ids, lag_window)
    preds_balanced <- run_lodo_model(balanced_dat, response, lag_window, structure)
    write_csv(preds_balanced, file.path(sampling_dir, paste0("lodo_", model_name, "_", lag_window, "_balanced_dam_month_hour.csv")))
    sampling_results[[paste(model_name, "balanced", sep = "_")]] <- prediction_metrics(preds_balanced) |>
      mutate(
        model = model_name,
        response = response,
        lag_window = lag_window,
        sampling_mode = "balanced_one_per_dam_month_hour",
        n_balanced_obs_ids = length(balanced_ids)
      )
    
    write_csv(
      tibble(model = model_name, response = response, lag_window = lag_window, obs_id = balanced_ids),
      file.path(sampling_dir, paste0("balanced_obs_ids_", model_name, "_", lag_window, ".csv"))
    )
  }
  
  sampling_sensitivity <- bind_rows(sampling_results) |>
    select(model, response, lag_window, sampling_mode, everything()) |>
    arrange(model, sampling_mode)
  
  write_csv(sampling_sensitivity, file.path(sampling_dir, "sampling_sensitivity_lodo_metrics.csv"))
}

message("Done Step3B sampling coverage and sensitivity.")
message("Outputs written to: ", output_dir)
