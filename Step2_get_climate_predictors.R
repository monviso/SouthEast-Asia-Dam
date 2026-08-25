# Prepare climatic predictors for ECOSTRESS dam temperature models.
#
# Input expected:
#   sites_complete_HQ_obs.csv with columns:
#   date.time, dam, Latitude_l, Longitude_l, Latitude_o, Longitude_o, Tm_l, Tm_o
#
# Climate source:
#   NASA POWER hourly point API, requested in UTC.
#
# Important:
#   ECOSTRESS timestamps are assumed to be UTC. Local hour/day predictors are
#   computed separately from the UTC timestamp. By default this script uses
#   local solar time from longitude, which is usually preferable for thermal
#   diel-cycle modelling across multiple Southeast Asian time zones.

required_packages <- c(
  "dplyr", "tidyr", "purrr", "readr", "lubridate", "jsonlite", "stringr"
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

library(dplyr)
library(tidyr)
library(purrr)
library(readr)
library(lubridate)
library(jsonlite)
library(stringr)

input_candidates <- c( "E:/OneDrive/Documents/satData_AppEEars/model_inputs/sites_complete_HQ_obs_ECOSTRESS.csv"
)

input_csv <- input_candidates[file.exists(input_candidates)][1]

if (is.na(input_csv)) {
  stop("Could not find an input CSV in /workspace/.cache.", call. = FALSE)
}

power_cache_dir <- "E:/OneDrive/Documents/satData_AppEEars/climate_cache_nasa_power"
output_dir <- "E:/OneDrive/Documents/satData_AppEEars/model_inputs/"

dir.create(power_cache_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# NASA POWER hourly parameters.
# T2M                 = 2 m air temperature, deg C
# PRECTOTCORR         = corrected precipitation, mm/hour
# WS2M                = 2 m wind speed, m/s
# PS                  = surface pressure, kPa
# ALLSKY_SFC_SW_DWN   = all-sky shortwave radiation at surface
power_parameters <- c(
  "T2M",
  "PRECTOTCORR",
  "WS2M",
  "PS",
  "ALLSKY_SFC_SW_DWN"
)

max_lag_hours <- 48L
lag_sets <- list(
  lag_4h = 0:4,
  lag_12h = 0:12,
  lag_24h = 0:24,
  lag_48h = 0:48
)

safe_name <- function(x) {
  x |>
    str_replace_all("[^A-Za-z0-9]+", "_") |>
    str_replace_all("^_|_$", "")
}

# More aggressive dam key used only to reconcile names that differ by spaces
# or punctuation, e.g. "Huai Kum" vs "HuaiKum".
canonical_dam_key <- function(x) {
  x |>
    as.character() |>
    str_to_lower() |>
    str_replace_all("[^a-z0-9]+", "")
}

median_or_na <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (all(is.na(x))) {
    NA_real_
  } else {
    median(x, na.rm = TRUE)
  }
}

first_non_missing_chr <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) {
    NA_character_
  } else {
    x[1]
  }
}

parse_power_hourly <- function(power_json) {
  pars <- power_json$properties$parameter
  time_keys <- names(pars[[1]])
  
  out <- tibble(
    datetime_utc = ymd_h(time_keys, tz = "UTC")
  )
  
  for (p in names(pars)) {
    vals <- unlist(pars[[p]][time_keys], use.names = FALSE)
    out[[p]] <- as.numeric(vals)
  }
  
  out |>
    mutate(across(-datetime_utc, ~ na_if(.x, -999))) |>
    transmute(
      datetime_utc,
      Air_T = T2M,
      prec = PRECTOTCORR,
      wind = WS2M,
      pres = PS,
      shortW = ALLSKY_SFC_SW_DWN
    )
}

download_power_hourly <- function(dam, lat, lon, start_date, end_date,
                                  cache_dir = power_cache_dir,
                                  parameters = power_parameters) {
  if (!is.finite(lat) || !is.finite(lon)) {
    stop(
      "Missing/invalid NASA POWER coordinates for dam: ", dam,
      " (lat=", lat, ", lon=", lon, ")",
      call. = FALSE
    )
  }
  
  cache_file <- file.path(
    cache_dir,
    paste0(safe_name(dam), "_", start_date, "_", end_date, ".csv")
  )
  
  if (file.exists(cache_file)) {
    return(read_csv(cache_file, show_col_types = FALSE))
  }
  
  base_url <- "https://power.larc.nasa.gov/api/temporal/hourly/point"
  query <- list(
    parameters = paste(parameters, collapse = ","),
    community = "AG",
    longitude = sprintf("%.5f", lon),
    latitude = sprintf("%.5f", lat),
    start = format(start_date, "%Y%m%d"),
    end = format(end_date, "%Y%m%d"),
    format = "JSON",
    `time-standard` = "UTC"
  )
  
  url <- paste0(
    base_url, "?",
    paste(
      paste0(names(query), "=", utils::URLencode(unlist(query), reserved = TRUE)),
      collapse = "&"
    )
  )
  
  message("Downloading NASA POWER for ", dam)
  power_json <- jsonlite::fromJSON(url, simplifyVector = FALSE)
  clim <- parse_power_hourly(power_json)
  
  write_csv(clim, cache_file)
  clim
}

add_local_predictors <- function(dat, lon_col = "Longitude_l",
                                 time_col = "datetime_utc",
                                 mode = c("solar", "timezone")) {
  mode <- match.arg(mode)
  
  if (mode == "solar") {
    # Solar local time: offset = longitude / 15 hours.
    # This avoids country timezone boundaries and is physically meaningful for
    # diurnal heating. Example: 105E is UTC+7 solar hours.
    dat |>
      mutate(
        local_offset_hours = .data[[lon_col]] / 15,
        datetime_local = .data[[time_col]] +
          seconds(local_offset_hours * 3600),
        H_LOC = hour(datetime_local) +
          minute(datetime_local) / 60 +
          second(datetime_local) / 3600,
        doy = yday(datetime_local),
        m_LOC = month(datetime_local),
        local_date = as.Date(datetime_local)
      )
  } else {
    if (!requireNamespace("lutz", quietly = TRUE)) {
      stop(
        "Install package 'lutz' for timezone mode, or use mode = 'solar'.",
        call. = FALSE
      )
    }
    
    dat$tz_name <- lutz::tz_lookup_coords(
      lat = dat$Latitude_l,
      lon = dat[[lon_col]],
      method = "accurate"
    )
    
    local_list <- map2(
      dat[[time_col]],
      dat$tz_name,
      ~ with_tz(.x, tzone = .y)
    )
    
    dat$H_LOC <- map_dbl(
      local_list,
      ~ hour(.x) + minute(.x) / 60 + second(.x) / 3600
    )
    dat$doy <- map_int(local_list, yday)
    dat$m_LOC <- map_int(local_list, month)
    dat$local_date <- as.Date(map_chr(local_list, as.character))
    dat$datetime_local_chr <- map_chr(local_list, as.character)
    dat
  }
}

extract_lagged_climate <- function(site_obs, site_clim, lags = 0:max_lag_hours) {
  target <- tidyr::crossing(
    obs_id = site_obs$obs_id,
    lag = lags
  ) |>
    left_join(
      site_obs |> select(obs_id, datetime_utc),
      by = "obs_id"
    ) |>
    mutate(target_time_utc = datetime_utc - hours(lag))
  
  climate_vars <- c("Air_T", "prec", "wind", "pres", "shortW")
  
  for (v in climate_vars) {
    target[[paste0(v, "_l")]] <- approx(
      x = as.numeric(site_clim$datetime_utc),
      y = site_clim[[v]],
      xout = as.numeric(target$target_time_utc),
      rule = 2,
      ties = "ordered"
    )$y
  }
  
  target |>
    select(obs_id, lag, target_time_utc, ends_with("_l"))
}

make_zero_lag_data <- function(obs, lag_long) {
  met0 <- lag_long |>
    filter(lag == 0) |>
    transmute(
      obs_id,
      Air_T = Air_T_l,
      prec = prec_l,
      wind = wind_l,
      pres = pres_l,
      shortW = shortW_l
    )
  
  obs |>
    left_join(met0, by = "obs_id") |>
    mutate(
      T_diff = Tm_o - Tm_l,
      Res_capacity_km3 = case_when(
        !is.na(Res_capaci_kmcube) & Res_capaci_kmcube > 20 ~
          Res_capaci_kmcube / 1000,
        TRUE ~ Res_capaci_kmcube
      ),
      mean_depth_m = if_else(
        !is.na(Res_capacity_km3) &
          !is.na(Res_area_kmsq) &
          Res_area_kmsq > 0,
        1000 * Res_capacity_km3 / Res_area_kmsq,
        NA_real_
      ),
      log_res_area = log1p(Res_area_kmsq),
      dam_id = factor(dam)
    )
}

matrix_from_lag_long <- function(lag_long, obs_ids, value_col, lags) {
  wide <- lag_long |>
    filter(obs_id %in% obs_ids, lag %in% lags) |>
    select(obs_id, lag, value = all_of(value_col)) |>
    mutate(lag = paste0("lag_", lag)) |>
    pivot_wider(names_from = lag, values_from = value) |>
    arrange(match(obs_id, obs_ids))
  
  lag_cols <- paste0("lag_", lags)
  as.matrix(wide[, lag_cols, drop = FALSE])
}

make_mgcv_lag_data <- function(dt0, lag_long, response,
                               lags = 0:24,
                               require_paired = FALSE) {
  keep <- !is.na(dt0[[response]])
  
  if (require_paired) {
    keep <- keep & !is.na(dt0$Tm_l) & !is.na(dt0$Tm_o)
  }
  
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

obs_raw <- read_csv(input_csv, show_col_types = FALSE) |>
  select(-matches("^\\.\\.\\.|^Unnamed"))

for (cc in c("Latitude_ref", "Longitude_ref", "latitude", "longitude",
             "Latitude_l", "Longitude_l", "Latitude_o", "Longitude_o",
             "Res_capaci_kmcube", "Res_area_kmsq")) {
  if (!cc %in% names(obs_raw)) {
    obs_raw[[cc]] <- NA_real_
  }
}

# Some rows from the additional ECOSTRESS file have dam names without spaces
# and no coordinates, e.g. "HuaiKum", while the main S2-mask rows use
# "Huai Kum" and do have coordinates. Build a canonical key that ignores
# spaces/punctuation, then fill missing coordinates from any matching row.
obs_pre <- obs_raw |>
  mutate(
    obs_id = row_number(),
    dam_original = as.character(dam),
    dam_key = canonical_dam_key(dam_original),
    datetime_utc = ymd_hms(date.time, tz = "UTC"),
    Latitude_ref = coalesce(Latitude_ref, latitude, Latitude_l, Latitude_o),
    Longitude_ref = coalesce(Longitude_ref, longitude, Longitude_l, Longitude_o)
  )

site_lookup <- obs_pre |>
  group_by(dam_key) |>
  summarise(
    dam_preferred = first_non_missing_chr(
      dam_original[!is.na(Latitude_ref) & !is.na(Longitude_ref)]
    ),
    dam_fallback = first_non_missing_chr(dam_original),
    Latitude_site = median_or_na(Latitude_ref),
    Longitude_site = median_or_na(Longitude_ref),
    n_rows_site_key = n(),
    names_combined = paste(sort(unique(dam_original)), collapse = " | "),
    .groups = "drop"
  ) |>
  mutate(
    dam_preferred = coalesce(dam_preferred, dam_fallback)
  ) |>
  select(-dam_fallback)

write_csv(
  site_lookup,
  file.path(output_dir, "site_coordinate_lookup_after_name_reconciliation.csv")
)

obs <- obs_pre |>
  left_join(site_lookup, by = "dam_key") |>
  mutate(
    dam = dam_preferred,
    Latitude_ref = coalesce(Latitude_ref, Latitude_site),
    Longitude_ref = coalesce(Longitude_ref, Longitude_site)
  )

sites <- obs |>
  group_by(dam_key, dam) |>
  summarise(
    lat = median_or_na(Latitude_ref),
    lon = median_or_na(Longitude_ref),
    start_date = as.Date(min(datetime_utc, na.rm = TRUE) - hours(max_lag_hours + 2)),
    end_date = as.Date(max(datetime_utc, na.rm = TRUE) + hours(2)),
    n_obs = n(),
    names_combined = paste(sort(unique(dam_original)), collapse = " | "),
    .groups = "drop"
  )

missing_sites <- sites |>
  filter(!is.finite(lat) | !is.finite(lon))

if (nrow(missing_sites) > 0) {
  write_csv(
    missing_sites,
    file.path(output_dir, "sites_missing_coordinates_for_nasa_power.csv")
  )
  
  warning(
    "Dropping ", nrow(missing_sites),
    " site(s) with no valid coordinates before NASA POWER download: ",
    paste(missing_sites$names_combined, collapse = "; "),
    call. = FALSE
  )
  
  obs <- obs |>
    filter(!dam_key %in% missing_sites$dam_key)
  
  sites <- sites |>
    filter(dam_key %in% unique(obs$dam_key))
}

obs <- obs |>
  add_local_predictors(lon_col = "Longitude_ref", mode = "solar")

climate_by_dam <- pmap(
  sites,
  function(dam_key, dam, lat, lon, start_date, end_date, n_obs, names_combined) {
    download_power_hourly(
      dam = dam,
      lat = lat,
      lon = lon,
      start_date = start_date,
      end_date = end_date
    ) |>
      mutate(
        dam_key = dam_key,
        dam = dam
      )
  }
)

climate_hourly <- bind_rows(climate_by_dam)

lag_long <- obs |>
  group_split(dam_key) |>
  map_dfr(function(site_obs) {
    this_dam_key <- unique(site_obs$dam_key)
    site_clim <- climate_hourly |>
      filter(dam_key == this_dam_key) |>
      arrange(datetime_utc)
    extract_lagged_climate(site_obs, site_clim, lags = 0:max_lag_hours)
  })

dt0 <- make_zero_lag_data(obs, lag_long)

# mgcv list-data objects for distributed lag models.
# Use require_paired = TRUE for T_diff because it needs both Tm_l and Tm_o.
gam_lag_data <- list(
  lake = map(lag_sets, ~ make_mgcv_lag_data(dt0, lag_long, "Tm_l", lags = .x)),
  outlet = map(lag_sets, ~ make_mgcv_lag_data(dt0, lag_long, "Tm_o", lags = .x)),
  diff = map(
    lag_sets,
    ~ make_mgcv_lag_data(
      dt0,
      lag_long,
      "T_diff",
      lags = .x,
      require_paired = TRUE
    )
  )
)

write_csv(dt0, file.path(output_dir, "dt0_nasa_power_zerolag.csv"))
write_csv(lag_long, file.path(output_dir, "lag_long_nasa_power_0_48h.csv"))
saveRDS(gam_lag_data, file.path(output_dir, "gam_lag_data_nasa_power.rds"))
saveRDS(climate_hourly, file.path(output_dir, "climate_hourly_nasa_power.rds"))

message("Wrote:")
message("  ", file.path(output_dir, "dt0_nasa_power_zerolag.csv"))
message("  ", file.path(output_dir, "lag_long_nasa_power_0_48h.csv"))
message("  ", file.path(output_dir, "gam_lag_data_nasa_power.rds"))
message("  ", file.path(output_dir, "climate_hourly_nasa_power.rds"))
