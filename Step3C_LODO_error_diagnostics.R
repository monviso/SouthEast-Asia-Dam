# Fit 0-lag and distributed-lag GAMs for reservoir lake temperature,
# outlet temperature, and outlet-lake temperature difference.
#
# Requires outputs from get_climate_predictors.R:
#   /workspace/model_inputs/dt0_nasa_power_zerolag.csv
#   /workspace/model_inputs/gam_lag_data_nasa_power.rds
#
# Main outputs:
#   /workspace/model_outputs/model_comparison_zerolag.csv
#   /workspace/model_outputs/model_comparison_lag.csv
#   /workspace/model_outputs/lodo_predictions_*.csv
#   /workspace/model_outputs/variance_partition_*.csv
#   /workspace/model_outputs/venn_*.pdf

required_packages <- c("mgcv", "dplyr", "purrr", "readr", "tibble")

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

input_dir <- "E:/OneDrive/Documents/satData_AppEEars/model_inputs"
output_dir <- "E:/OneDrive/Documents/satData_AppEEars/model_outputs"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
summary_dir <- file.path(output_dir, "model_summaries")
effect_dir <- file.path(output_dir, "significant_partial_effects")
plot_log_path <- file.path(output_dir, "partial_effect_plot_log.csv")
dir.create(summary_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(effect_dir, showWarnings = FALSE, recursive = TRUE)

dt0_path <- file.path(input_dir, "dt0_nasa_power_zerolag.csv")
lag_path <- file.path(input_dir, "gam_lag_data_nasa_power.rds")

dt0 <- read_csv(dt0_path, show_col_types = FALSE)

if (!"has_lake" %in% names(dt0)) {
  dt0$has_lake <- !is.na(dt0$Tm_l)
}
if (!"has_outlet" %in% names(dt0)) {
  dt0$has_outlet <- !is.na(dt0$Tm_o)
}
if (!"paired" %in% names(dt0)) {
  dt0$paired <- dt0$has_lake & dt0$has_outlet
}

dt0 <- dt0 |>
  mutate(
    dam_id = factor(dam),
    T_diff = if_else(paired, Tm_o - Tm_l, NA_real_)
  )

gam_lag_data <- readRDS(lag_path)

# -------------------------------------------------------------------------
# Optional dam-characteristic handling
# -------------------------------------------------------------------------
# If your zero-lag file already contains dam characteristics, this script uses
# them automatically. If not, join them before fitting, for example:
#
# dam_covariates <- read_csv("/path/to/dam_characteristics.csv")
# dt0 <- dt0 |> left_join(dam_covariates, by = "dam")
#
# Then rebuild gam_lag_data from get_climate_predictors.R or join the same
# columns into each list element before fitting.

add_derived_dam_terms <- function(dat) {
  if (all(c("Res_capaci_kmcube", "Res_area_kmsq") %in% names(dat))) {
    dat <- dat |>
      mutate(
        Res_capacity_km3 = case_when(
          !is.na(Res_capaci_kmcube) & Res_capaci_kmcube > 20 ~
            Res_capaci_kmcube / 1000,
          TRUE ~ Res_capaci_kmcube
        ),
        mean_depth_m = if_else(
          Res_area_kmsq > 0,
          1000 * Res_capacity_km3 / Res_area_kmsq,
          NA_real_
        ),
        log_res_area = log1p(Res_area_kmsq)
      )
  }
  dat
}

dt0 <- add_derived_dam_terms(dt0)

add_dam_terms_to_lag_data <- function(x, dt0) {
  for (nm in names(x)) {
    if (is.null(x[[nm]]$obs_id)) next

    obs_ids <- x[[nm]]$obs_id

    add_cols <- dt0 |>
      select(any_of(c(
        "obs_id", "Dam_hgt", "elevation", "Res_capaci_kmcube",
        "Res_capacity_km3", "Res_area_kmsq", "mean_depth_m", "log_res_area"
      ))) |>
      filter(obs_id %in% obs_ids) |>
      arrange(match(obs_id, obs_ids))

    if (ncol(add_cols) <= 1) next

    if (nrow(add_cols) != length(obs_ids)) {
      stop(
        "Could not align dam covariates for ", nm,
        ": expected ", length(obs_ids), " rows, got ", nrow(add_cols),
        ". Check obs_id consistency between dt0 and gam_lag_data.",
        call. = FALSE
      )
    }

    for (cc in setdiff(names(add_cols), "obs_id")) {
      x[[nm]][[cc]] <- add_cols[[cc]]
    }
  }
  x
}

gam_lag_data$lake <- add_dam_terms_to_lag_data(gam_lag_data$lake, dt0)
gam_lag_data$outlet <- add_dam_terms_to_lag_data(gam_lag_data$outlet, dt0)
gam_lag_data$diff <- add_dam_terms_to_lag_data(gam_lag_data$diff, dt0)

dam_predictors <- c(
  "Dam_hgt",
  "elevation",
  "mean_depth_m",
  "log_res_area"
)

dam_predictors <- dam_predictors[dam_predictors %in% names(dt0)]

if (length(dam_predictors) == 0) {
  message(
    "No dam-characteristic columns found. Predictive models will use ",
    "cyclic + meteorological terms + dam random effect only."
  )
}

# -------------------------------------------------------------------------
# Formula builders
# -------------------------------------------------------------------------

smooth_terms <- function(vars, k = 6) {
  paste0("s(", vars, ", bs = 'ts', k = ", k, ")")
}

cyclic_terms <- function(dat) {
  out <- character()
  if ("H_LOC" %in% names(dat)) {
    out <- c(out, "s(H_LOC, bs = 'cc', k = 8)")
  }
  if ("doy" %in% names(dat)) {
    out <- c(out, "s(doy, bs = 'cc', k = 12)")
  } else if ("m_LOC" %in% names(dat)) {
    out <- c(out, "s(m_LOC, bs = 'cc', k = 12)")
  }
  out
}

dam_terms <- function(dat) {
  vars <- dam_predictors[dam_predictors %in% names(dat)]
  smooth_terms(vars, k = 6)
}

meteo_zerolag_terms <- function(dat) {
  vars <- c("Air_T", "prec", "wind", "shortW", "pres")
  vars <- vars[vars %in% names(dat)]
  smooth_terms(vars, k = 6)
}

meteo_lag_terms <- function(lags) {
  k_lag <- min(6L, length(lags))
  c(
    paste0("te(Air_T_l, lag, k = c(6, ", k_lag, "), bs = c('ts', 'ts'))"),
    paste0("te(prec_l, lag, k = c(6, ", k_lag, "), bs = c('ts', 'ts'))"),
    paste0("te(wind_l, lag, k = c(6, ", k_lag, "), bs = c('ts', 'ts'))"),
    paste0("te(shortW_l, lag, k = c(6, ", k_lag, "), bs = c('ts', 'ts'))"),
    paste0("te(pres_l, lag, k = c(6, ", k_lag, "), bs = c('ts', 'ts'))")
  )
}

make_formula <- function(response, terms, add_re = TRUE) {
  if (add_re) {
    terms <- c(terms, "s(dam_id, bs = 're')")
  }

  as.formula(
    paste(response, "~", paste(terms, collapse = " + "))
  )
}

make_zerolag_formula <- function(response, dat = dt0, add_re = TRUE) {
  terms <- c(
    meteo_zerolag_terms(dat),
    cyclic_terms(dat),
    dam_terms(dat)
  )
  make_formula(response, terms, add_re = add_re)
}

make_lag_formula <- function(response, lags, dat, add_re = TRUE) {
  terms <- c(
    meteo_lag_terms(lags),
    cyclic_terms(dat),
    dam_terms(dat)
  )
  make_formula(response, terms, add_re = add_re)
}

cyclic_knots <- function(dat) {
  out <- list()
  if ("H_LOC" %in% names(dat)) out$H_LOC <- c(0, 24)
  if ("doy" %in% names(dat)) out$doy <- c(0.5, 366.5)
  if (!"doy" %in% names(dat) && "m_LOC" %in% names(dat)) {
    out$m_LOC <- c(0.5, 12.5)
  }
  out
}

# -------------------------------------------------------------------------
# Metrics and leave-one-dam-out CV
# -------------------------------------------------------------------------

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

clean_file_label <- function(x) {
  gsub("[^A-Za-z0-9_]+", "_", x)
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

  write_csv(
    global,
    file.path(summary_dir, paste0("global_", safe_id, ".csv"))
  )

  if (!is.null(sm$p.table)) {
    parametric <- as.data.frame(sm$p.table) |>
      rownames_to_column("term") |>
      as_tibble() |>
      mutate(model_id = model_id, .before = 1)

    write_csv(
      parametric,
      file.path(summary_dir, paste0("parametric_terms_", safe_id, ".csv"))
    )
  }

  if (!is.null(sm$s.table)) {
    smooth <- as.data.frame(sm$s.table) |>
      rownames_to_column("term") |>
      as_tibble() |>
      mutate(model_id = model_id, .before = 1)

    write_csv(
      smooth,
      file.path(summary_dir, paste0("smooth_terms_", safe_id, ".csv"))
    )

    save_significant_partial_effects(fit, smooth, model_id)
  }
}

save_significant_partial_effects <- function(fit, smooth_table, model_id,
                                             alpha = 0.05) {
  p_col <- grep("p", names(smooth_table), ignore.case = TRUE, value = TRUE)
  p_col <- p_col[length(p_col)]

  if (length(p_col) == 0 || is.na(p_col)) {
    return(invisible(NULL))
  }

  sig_terms <- smooth_table |>
    filter(
      is.finite(.data[[p_col]]),
      .data[[p_col]] < alpha,
      !grepl("dam_id", term)
    )

  if (nrow(sig_terms) == 0) {
    return(invisible(NULL))
  }

  smooth_labels <- vapply(fit$smooth, function(x) x$label, character(1))
  safe_model_id <- clean_file_label(model_id)
  plot_log <- list()

  for (j in seq_len(nrow(sig_terms))) {
    term <- sig_terms$term[j]
    select_idx <- match(term, smooth_labels)

    if (is.na(select_idx)) next

    file_name <- paste0(
      safe_model_id,
      "__",
      clean_file_label(term),
      ".png"
    )

    out_file <- file.path(effect_dir, file_name)

    status <- "saved"
    error_message <- NA_character_

    tryCatch(
      {
        png(out_file, width = 1800, height = 1400, res = 180)
        plot(
          fit,
          select = select_idx,
          shade = TRUE,
          pages = 1,
          scheme = 2,
          main = paste(model_id, term)
        )
        dev.off()
      },
      error = function(e) {
        status <<- "failed"
        error_message <<- conditionMessage(e)
        if (names(dev.cur()) != "null device") {
          dev.off()
        }
      }
    )

    if (names(dev.cur()) != "null device") {
      dev.off()
    }

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
      write_csv(
        bind_rows(read_csv(plot_log_path, show_col_types = FALSE), log_df),
        plot_log_path
      )
    } else {
      write_csv(log_df, plot_log_path)
    }
  }
}

prediction_metrics <- function(predictions) {
  safe_cor <- function(x, y) {
    if (length(x) < 3 || sd(x) == 0 || sd(y) == 0) {
      return(NA_real_)
    }
    cor(x, y)
  }

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
      R2_CV = 1 - sum(sq_error) / sum((observed - mean(observed))^2),
      .groups = "drop"
    )
}

by_dam_metrics <- function(predictions) {
  safe_cor <- function(x, y) {
    if (length(x) < 3 || sd(x) == 0 || sd(y) == 0) {
      return(NA_real_)
    }
    cor(x, y)
  }

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

complete_rows_df <- function(data, formula) {
  vars <- all.vars(formula)
  vars <- vars[vars %in% names(data)]
  scalar_vars <- vars[!vapply(data[vars], is.matrix, logical(1))]

  if (length(scalar_vars) == 0) {
    return(rep(TRUE, nrow(data)))
  }

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
      if (nrow(x) != n) {
        stop(
          "Variable ", v, " has ", nrow(x),
          " rows but expected ", n, ".",
          call. = FALSE
        )
      }

      if (is.numeric(x) || is.integer(x)) {
        keep <- keep & rowSums(!is.finite(x)) == 0
      } else {
        keep <- keep & rowSums(is.na(x)) == 0
      }
    } else if (length(x) == n) {
      keep <- keep & !is.na(x)
    } else if (length(x) == n * length(unique(as.vector(data$lag)))) {
      xmat <- matrix(x, nrow = n, byrow = FALSE)
      keep <- keep & rowSums(is.na(xmat)) == 0
    } else if (length(x) != 1) {
      stop(
        "Variable ", v, " has length ", length(x),
        " but expected ", n, ".",
        call. = FALSE
      )
    }
  }

  keep
}

lodo_gam_df <- function(data, formula, response,
                        dam_col = "dam_id",
                        knots = cyclic_knots(data),
                        re_term = NULL) {
  data <- data |>
    filter(!is.na(.data[[response]]))

  dams <- unique(as.character(data[[dam_col]]))
  predictions <- vector("list", length(dams))

  for (i in seq_along(dams)) {
    held_out <- dams[i]

    train <- data[as.character(data[[dam_col]]) != held_out, , drop = FALSE]
    test <- data[as.character(data[[dam_col]]) == held_out, , drop = FALSE]

    train[[dam_col]] <- droplevels(factor(train[[dam_col]]))

    fit <- gam(
      formula,
      data = train,
      family = gaussian(),
      method = "REML",
      knots = knots
    )

    newdata <- test

    if (is.null(re_term)) {
      pred <- predict(fit, newdata = newdata, type = "response")
    } else {
      pred <- predict(
        fit,
        newdata = newdata,
        type = "response",
        exclude = re_term
      )
    }

    predictions[[i]] <- tibble(
      obs_id = if ("obs_id" %in% names(test)) test$obs_id else NA_integer_,
      dam_id = held_out,
      observed = test[[response]],
      predicted = as.numeric(pred)
    )
  }

  bind_rows(predictions)
}

subset_gam_list <- function(dat, rows) {
  lapply(dat, function(x) {
    if (is.matrix(x)) {
      x[rows, , drop = FALSE]
    } else {
      x[rows]
    }
  })
}

lodo_gam_list <- function(data, formula, response,
                          dam_col = "dam_id",
                          knots = cyclic_knots(data),
                          re_term = NULL) {
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

    train[[dam_col]] <- droplevels(factor(train[[dam_col]]))

    fit <- gam(
      formula,
      data = train,
      family = gaussian(),
      method = "REML",
      knots = knots
    )

    newdata <- test

    if (is.null(re_term)) {
      pred <- predict(fit, newdata = newdata, type = "response")
    } else {
      pred <- predict(
        fit,
        newdata = newdata,
        type = "response",
        exclude = re_term
      )
    }

    predictions[[i]] <- tibble(
      obs_id = if ("obs_id" %in% names(test)) test$obs_id else NA_integer_,
      dam_id = held_out,
      observed = test[[response]],
      predicted = as.numeric(pred)
    )
  }

  bind_rows(predictions)
}

# -------------------------------------------------------------------------
# Fit 0-lag models
# -------------------------------------------------------------------------

zero_specs <- tribble(
  ~model, ~response,
  "Tlake", "Tm_l",
  "Toutlet", "Tm_o",
  "Tdiff", "T_diff"
)

zero_results <- vector("list", nrow(zero_specs))

model_sample_sizes <- tibble(
  response = c("Tm_l", "Tm_o", "T_diff"),
  n = c(
    sum(!is.na(dt0$Tm_l)),
    sum(!is.na(dt0$Tm_o)),
    sum(dt0$paired & !is.na(dt0$T_diff))
  ),
  n_dams = c(
    n_distinct(dt0$dam_id[!is.na(dt0$Tm_l)]),
    n_distinct(dt0$dam_id[!is.na(dt0$Tm_o)]),
    n_distinct(dt0$dam_id[dt0$paired & !is.na(dt0$T_diff)])
  )
)

write_csv(model_sample_sizes, file.path(output_dir, "model_sample_sizes.csv"))
print(model_sample_sizes)

for (i in seq_len(nrow(zero_specs))) {
  model_name <- zero_specs$model[i]
  response <- zero_specs$response[i]

  message("Fitting zero-lag model: ", model_name)

  dat <- dt0 |>
    filter(!is.na(.data[[response]]))

  form <- make_zerolag_formula(response, dat, add_re = TRUE)
  form_cv <- make_zerolag_formula(response, dat, add_re = FALSE)

  dat <- dat[complete_rows_df(dat, form), , drop = FALSE]

  fit <- gam(
    form,
    data = dat,
    family = gaussian(),
    method = "REML",
    knots = cyclic_knots(dat)
  )

  save_full_model_summary(
    fit,
    model_id = paste0("zerolag_", model_name)
  )

  preds <- lodo_gam_df(
    data = dat,
    formula = form_cv,
    response = response,
    knots = cyclic_knots(dat)
  )

  write_csv(
    preds,
    file.path(output_dir, paste0("lodo_predictions_zerolag_", model_name, ".csv"))
  )

  write_csv(
    by_dam_metrics(preds),
    file.path(output_dir, paste0("lodo_by_dam_zerolag_", model_name, ".csv"))
  )

  zero_results[[i]] <- model_metrics(fit) |>
    bind_cols(prediction_metrics(preds)) |>
    mutate(
      model = model_name,
      response = response,
      lag_window = "0h",
      formula = paste(deparse(form), collapse = " ")
    ) |>
    select(model, response, lag_window, everything())

  saveRDS(fit, file.path(output_dir, paste0("fit_zerolag_", model_name, ".rds")))
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
    message("Fitting lag model: ", model_name, " / ", lag_window)

    dat <- gam_lag_data[[data_name]][[lag_window]]
    lags <- sort(unique(as.vector(dat$lag)))
    form <- make_lag_formula(response, lags, dat, add_re = TRUE)
    form_cv <- make_lag_formula(response, lags, dat, add_re = FALSE)
    dat <- subset_gam_list(dat, complete_rows_list(dat, form))

    fit <- gam(
      form,
      data = dat,
      family = gaussian(),
      method = "REML",
      knots = cyclic_knots(dat)
    )

    save_full_model_summary(
      fit,
      model_id = paste0(lag_window, "_", model_name)
    )

    preds <- lodo_gam_list(
      data = dat,
      formula = form_cv,
      response = response,
      knots = cyclic_knots(dat)
    )

    pred_file <- file.path(
      output_dir,
      paste0("lodo_predictions_", lag_window, "_", model_name, ".csv")
    )

    write_csv(preds, pred_file)

    write_csv(
      by_dam_metrics(preds),
      file.path(output_dir, paste0("lodo_by_dam_", lag_window, "_", model_name, ".csv"))
    )

    lag_results[[paste(model_name, lag_window, sep = "_")]] <- model_metrics(fit) |>
      bind_cols(prediction_metrics(preds)) |>
      mutate(
        model = model_name,
        response = response,
        lag_window = lag_window,
        formula = paste(deparse(form), collapse = " ")
      ) |>
      select(model, response, lag_window, everything())

    saveRDS(
      fit,
      file.path(output_dir, paste0("fit_", lag_window, "_", model_name, ".rds"))
    )
  }
}

lag_comparison <- bind_rows(lag_results)
write_csv(lag_comparison, file.path(output_dir, "model_comparison_lag.csv"))

# -------------------------------------------------------------------------
# Diagnostic model: where is Tdiff prediction error largest?
# -------------------------------------------------------------------------
# This is a second-stage diagnostic model. It uses the absolute LODO error of
# the selected Tdiff model as the response. This is better than using in-sample
# residuals because it targets conditions where the model generalises poorly to
# held-out dams.

fit_tdiff_error_diagnostic <- function(lag_comparison, gam_lag_data) {
  if (!"Tdiff" %in% lag_comparison$model) {
    warning("No Tdiff lag models available; skipping Tdiff error diagnostic.")
    return(invisible(NULL))
  }

  best <- lag_comparison |>
    filter(model == "Tdiff", is.finite(RMSE)) |>
    arrange(RMSE) |>
    slice(1)

  if (nrow(best) == 0) {
    warning("No finite Tdiff RMSE found; skipping Tdiff error diagnostic.")
    return(invisible(NULL))
  }

  lag_window <- best$lag_window[1]
  pred_file <- file.path(
    output_dir,
    paste0("lodo_predictions_", lag_window, "_Tdiff.csv")
  )

  if (!file.exists(pred_file)) {
    warning("Prediction file not found for Tdiff error diagnostic: ", pred_file)
    return(invisible(NULL))
  }

  preds <- read_csv(pred_file, show_col_types = FALSE) |>
    mutate(
      Tdiff_error = predicted - observed,
      Tdiff_abs_error = abs(Tdiff_error),
      Tdiff_abs_error_eps = Tdiff_abs_error + 0.01
    ) |>
    filter(!is.na(obs_id), is.finite(Tdiff_abs_error_eps))

  dat <- gam_lag_data$diff[[lag_window]]
  match_idx <- match(dat$obs_id, preds$obs_id)

  dat$Tdiff_error <- preds$Tdiff_error[match_idx]
  dat$Tdiff_abs_error <- preds$Tdiff_abs_error[match_idx]
  dat$Tdiff_abs_error_eps <- preds$Tdiff_abs_error_eps[match_idx]

  lags <- sort(unique(as.vector(dat$lag)))
  form <- make_lag_formula(
    "Tdiff_abs_error_eps",
    lags = lags,
    dat = dat,
    add_re = TRUE
  )

  dat <- subset_gam_list(dat, complete_rows_list(dat, form))

  message(
    "Fitting Tdiff absolute-error diagnostic model for best lag window: ",
    lag_window
  )

  fit <- gam(
    form,
    data = dat,
    family = Gamma(link = "log"),
    method = "REML",
    knots = cyclic_knots(dat)
  )

  saveRDS(
    fit,
    file.path(output_dir, paste0("fit_Tdiff_abs_error_", lag_window, ".rds"))
  )

  save_full_model_summary(
    fit,
    model_id = paste0("Tdiff_abs_error_", lag_window)
  )

  scalar_vars <- names(dat)[!vapply(dat, function(x) !is.null(dim(x)), logical(1))]

  write_csv(
    as_tibble(dat[scalar_vars]) |>
      select(
        any_of(c(
          "obs_id", "dam", "dam_id", "T_diff", "Tdiff_error",
          "Tdiff_abs_error", "Tdiff_abs_error_eps", "H_LOC", "doy",
          "Dam_hgt", "elevation", "mean_depth_m", "log_res_area"
        ))
      ),
    file.path(output_dir, paste0("Tdiff_abs_error_data_", lag_window, ".csv"))
  )

  diagnostic_summary <- tibble(
    selected_lag_window = lag_window,
    source_model_RMSE = best$RMSE[1],
    source_model_MAE = best$MAE[1],
    n = length(fit$y),
    error_model_AIC = AIC(fit),
    error_model_dev_expl = summary(fit)$dev.expl
  )

  write_csv(
    diagnostic_summary,
    file.path(output_dir, "Tdiff_abs_error_diagnostic_summary.csv")
  )

  invisible(fit)
}

fit_tdiff_error_diagnostic(lag_comparison, gam_lag_data)

# -------------------------------------------------------------------------
# Variance partition / Venn-style commonality analysis
# -------------------------------------------------------------------------
# This fits all seven non-empty combinations of predictor groups and uses
# deviance explained to compute approximate unique/shared components.
#
# Important interpretation:
#   - This is not causal attribution.
#   - Negative shared components can occur with correlated predictors,
#     suppressor effects, or smoothing penalties.
#   - The default excludes s(dam_id, bs='re') because the requested group is
#     dam characteristics, not dam identity. Set include_site_re = TRUE only
#     if you want dam identity included as part of the dam group.

fit_group_subset <- function(data, response, groups, selected_groups,
                             include_site_re = FALSE) {
  terms <- unlist(groups[selected_groups], use.names = FALSE)

  if (length(terms) == 0) {
    stop("No terms available for selected groups: ",
         paste(selected_groups, collapse = ", "))
  }

  if (include_site_re && "dam" %in% selected_groups) {
    terms <- c(terms, "s(dam_id, bs = 're')")
  }

  form <- make_formula(response, terms, add_re = FALSE)
  data <- data[complete_rows_df(data, form), , drop = FALSE]

  gam(
    form,
    data = data,
    family = gaussian(),
    method = "REML",
    knots = cyclic_knots(data)
  )
}

variance_partition <- function(data, response,
                               include_site_re = FALSE,
                               prefix = response) {
  data <- data |>
    filter(!is.na(.data[[response]]))

  groups <- list(
    cyclic = cyclic_terms(data),
    meteo = meteo_zerolag_terms(data),
    dam = dam_terms(data)
  )

  groups <- groups[lengths(groups) > 0]

  if (length(groups) < 2) {
    warning("Need at least two non-empty predictor groups for partitioning.")
    return(NULL)
  }

  group_names <- names(groups)

  combos <- unlist(
    map(seq_along(group_names), ~ combn(group_names, .x, simplify = FALSE)),
    recursive = FALSE
  )

  subset_results <- map_dfr(combos, function(g) {
    fit <- fit_group_subset(
      data = data,
      response = response,
      groups = groups,
      selected_groups = g,
      include_site_re = include_site_re
    )

    tibble(
      groups = paste(g, collapse = "+"),
      n_groups = length(g),
      dev_expl = summary(fit)$dev.expl,
      AIC = AIC(fit)
    )
  })

  write_csv(
    subset_results,
    file.path(output_dir, paste0("variance_partition_subsets_", prefix, ".csv"))
  )

  # Commonality coefficients for exactly three groups.
  if (length(group_names) == 3) {
    R <- function(g) {
      key <- paste(g, collapse = "+")
      subset_results$dev_expl[match(key, subset_results$groups)]
    }

    A <- group_names[1]
    B <- group_names[2]
    C <- group_names[3]

    R_A <- R(A)
    R_B <- R(B)
    R_C <- R(C)
    R_AB <- R(c(A, B))
    R_AC <- R(c(A, C))
    R_BC <- R(c(B, C))
    R_ABC <- R(c(A, B, C))

    regions <- tibble(
      region = c(
        A,
        B,
        C,
        paste(A, B, sep = "&"),
        paste(A, C, sep = "&"),
        paste(B, C, sep = "&"),
        paste(A, B, C, sep = "&"),
        "unexplained"
      ),
      value = c(
        R_ABC - R_BC,
        R_ABC - R_AC,
        R_ABC - R_AB,
        R_AC + R_BC - R_C - R_ABC,
        R_AB + R_BC - R_B - R_ABC,
        R_AB + R_AC - R_A - R_ABC,
        R_A + R_B + R_C - R_AB - R_AC - R_BC + R_ABC,
        1 - R_ABC
      )
    )

    write_csv(
      regions,
      file.path(output_dir, paste0("variance_partition_regions_", prefix, ".csv"))
    )

    draw_venn_base(
      regions = regions,
      group_names = group_names,
      output_pdf = file.path(output_dir, paste0("venn_", prefix, ".pdf")),
      title = paste("Approximate variance partition:", prefix)
    )
  }

  subset_results
}

draw_venn_base <- function(regions, group_names, output_pdf, title) {
  fmt <- function(region) {
    val <- regions$value[match(region, regions$region)]
    ifelse(
      is.na(val),
      "",
      paste0(region, "\n", sprintf("%.1f%%", 100 * val))
    )
  }

  A <- group_names[1]
  B <- group_names[2]
  C <- group_names[3]

  pdf(output_pdf, width = 7, height = 6)
  op <- par(mar = c(1, 1, 3, 1))
  plot.new()
  plot.window(xlim = c(0, 10), ylim = c(0, 9), asp = 1)
  title(main = title, cex.main = 1.1)

  symbols(4, 5.2, circles = 2.2, inches = FALSE, add = TRUE,
          bg = adjustcolor("#4E79A7", alpha.f = 0.25),
          fg = "#4E79A7", lwd = 2)
  symbols(6, 5.2, circles = 2.2, inches = FALSE, add = TRUE,
          bg = adjustcolor("#F28E2B", alpha.f = 0.25),
          fg = "#F28E2B", lwd = 2)
  symbols(5, 3.5, circles = 2.2, inches = FALSE, add = TRUE,
          bg = adjustcolor("#59A14F", alpha.f = 0.25),
          fg = "#59A14F", lwd = 2)

  text(2.5, 7.3, A, col = "#4E79A7", font = 2)
  text(7.5, 7.3, B, col = "#F28E2B", font = 2)
  text(5, 1.0, C, col = "#59A14F", font = 2)

  text(3.1, 5.6, fmt(A), cex = 0.8)
  text(6.9, 5.6, fmt(B), cex = 0.8)
  text(5, 2.7, fmt(C), cex = 0.8)
  text(5, 6.1, fmt(paste(A, B, sep = "&")), cex = 0.8)
  text(4.1, 4.1, fmt(paste(A, C, sep = "&")), cex = 0.8)
  text(5.9, 4.1, fmt(paste(B, C, sep = "&")), cex = 0.8)
  text(5, 4.9, fmt(paste(A, B, C, sep = "&")), cex = 0.8)

  unexplained <- regions$value[match("unexplained", regions$region)]
  text(
    5, 0.3,
    paste0("Unexplained: ", sprintf("%.1f%%", 100 * unexplained)),
    cex = 0.85
  )

  par(op)
  dev.off()
}

variance_partition(dt0, "Tm_l", prefix = "Tlake_zerolag")
variance_partition(dt0, "Tm_o", prefix = "Toutlet_zerolag")
variance_partition(dt0, "T_diff", prefix = "Tdiff_zerolag")

message("Done. Outputs written to: ", output_dir)
