# -----------------------------------------------------------------------------
# Update ECOSTRESS dam figures from the revised observation table
# -----------------------------------------------------------------------------
# Outputs:
#   1) Map with retained HQ observations coloured by number of observations,
#      candidate dams without HQ observations in orange, and Sungai Selangor in blue.
#   2) Temperature distribution plot for Lake WST, Outlet WST, and Tdiff.
#
# Main input is the revised ECOSTRESS observation table, e.g.:
#   sites_complete_HQ_obs_ECOSTRESS.csv
#
# Optional input for orange "no HQ observation" points:
#   a candidate dam shapefile, e.g. site_full_ext.shp, OR a CSV with dam/name + lat/lon.
#   If no candidate file is found, the map is still made for HQ sites, but orange points
#   cannot be drawn because the observation table only contains retained observations.
# -----------------------------------------------------------------------------

required_packages <- c(
  "dplyr", "readr", "ggplot2", "stringr", "tidyr", "lubridate",
  "sf", "maps", "viridis", "scales"
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
library(readr)
library(ggplot2)
library(stringr)
library(tidyr)
library(lubridate)
library(sf)
library(maps)
library(viridis)
library(scales)

# -----------------------------------------------------------------------------
# Local paths: edit only this block if needed
# -----------------------------------------------------------------------------

input_csv <- "E:/OneDrive/Documents/satData_AppEEars/model_inputs/sites_complete_HQ_obs_ECOSTRESS.csv"
inp<-read.csv(input_csv)
mean(inp$Tm_l, na.rm=T); sd(inp$Tm_l, na.rm=T)
mean(inp$Tm_o, na.rm=T); sd(inp$Tm_o, na.rm=T)
mean(inp$T_diff, na.rm=T); sd(inp$T_diff, na.rm=T)
median(inp$T_diff, na.rm=T)
# Candidate dam file is only needed to plot orange "no HQ observation" dams.
# Option A: shapefile with lake/river polygons and a dam-name column.
candidate_shp <- "E:/OneDrive/Documents/satData_AppEEars/csv_dam/shp/sites_full.shp"
cand<-read_sf(candidate_shp)
# Option B: CSV with dam/name plus latitude/longitude columns.
#candidate_csv <- "E:/OneDrive/Documents/satData_AppEEars/model_inputs/site_coordinate_lookup_after_name_reconciliation.csv"

output_dir <- "E:/OneDrive/Documents/satData_AppEEars/figures"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# If TRUE, point size on map also varies with number of observations.
scale_point_size_by_n <- FALSE

# Temperature thresholds are optional; leave as NULL or set, e.g. c(28, 30, 32)
temperature_thresholds <- NULL

# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------

normalise_dam_key <- function(x) {
  x |>
    as.character() |>
    str_to_lower() |>
    str_replace_all("[^a-z0-9]+", "") |>
    str_trim()
}

first_existing <- function(x) {
  x[file.exists(x)][1]
}

coalesce_numeric <- function(...) {
  dplyr::coalesce(!!!rlang::list2(...))
}

find_first_col <- function(dat, candidates) {
  hit <- candidates[candidates %in% names(dat)][1]
  if (length(hit) == 0 || is.na(hit)) NA_character_ else hit
}

# -----------------------------------------------------------------------------
# Read revised observation table
# -----------------------------------------------------------------------------

if (!file.exists(input_csv)) {
  stop("Input CSV not found: ", input_csv, call. = FALSE)
}

obs_raw <- read_csv(input_csv, show_col_types = FALSE) |>
  select(-matches("^\\.\\.\\.|^Unnamed"))

# Make sure optional columns exist before mutate()
for (cc in c(
  "dam", "datetime_utc", "date.time",
  "year", "month", "H_LOC",
  "has_lake", "has_outlet", "paired", "T_diff",
  "Latitude_ref", "Longitude_ref",
  "latitude", "longitude",
  "Latitude_l", "Longitude_l", "Latitude_o", "Longitude_o"
)) {
  if (!cc %in% names(obs_raw)) {
    obs_raw[[cc]] <- NA
  }
}

# Create/repair key columns outside mutate()
if (all(is.na(obs_raw$dam_key))) {
  obs_raw$dam_key <- normalise_dam_key(obs_raw$dam)
}

if (all(is.na(obs_raw$datetime_utc))) {
  obs_raw$datetime_utc <- suppressWarnings(ymd_hms(obs_raw$date.time, tz = "UTC"))
} else {
  obs_raw$datetime_utc <- suppressWarnings(ymd_hms(obs_raw$datetime_utc, tz = "UTC"))
}

lon_tmp <- coalesce_numeric(
  obs_raw$Longitude_ref,
  obs_raw$longitude,
  obs_raw$Longitude_l,
  obs_raw$Longitude_o
)

lat_tmp <- coalesce_numeric(
  obs_raw$Latitude_ref,
  obs_raw$latitude,
  obs_raw$Latitude_l,
  obs_raw$Latitude_o
)

if (all(is.na(obs_raw$H_LOC))) {
  obs_raw$H_LOC <- (
    lubridate::hour(obs_raw$datetime_utc) +
      lubridate::minute(obs_raw$datetime_utc) / 60 +
      lon_tmp / 15
  ) %% 24
}

obs <- obs_raw |>
  mutate(
    dam = as.character(dam),
    year = if_else(is.na(year), lubridate::year(datetime_utc), as.integer(year)),
    month = if_else(is.na(month), lubridate::month(datetime_utc), as.integer(month)),
    has_lake = if_else(is.na(has_lake), !is.na(Tm_l), as.logical(has_lake)),
    has_outlet = if_else(is.na(has_outlet), !is.na(Tm_o), as.logical(has_outlet)),
    paired = if_else(is.na(paired), !is.na(Tm_l) & !is.na(Tm_o), as.logical(paired)),
    T_diff = if_else(is.na(T_diff), Tm_o - Tm_l, as.numeric(T_diff)),
    Latitude_plot = lat_tmp,
    Longitude_plot = lon_tmp
  )

# One row per retained dam in the observation table.
obs_sites <- obs |>
  group_by(dam_key) |>
  summarise(
    dam = dplyr::first(na.omit(dam)),
    Latitude = median(Latitude_plot, na.rm = TRUE),
    Longitude = median(Longitude_plot, na.rm = TRUE),
    n_obs = n_distinct(datetime_utc),
    n_lake = sum(!is.na(Tm_l)),
    n_outlet = sum(!is.na(Tm_o)),
    n_paired = sum(!is.na(Tm_l) & !is.na(Tm_o)),
    .groups = "drop"
  ) |>
  filter(is.finite(Latitude), is.finite(Longitude)) |>
  mutate(
    site_class = if_else(
      normalise_dam_key(dam) == normalise_dam_key("Sungai Selangor"),
      "Sungai Selangor validation site",
      "Retained ECOSTRESS WST observations"
    )
  )

# -----------------------------------------------------------------------------
# Candidate dams for orange no-HQ points
# -----------------------------------------------------------------------------

candidate_sites <- NULL

if (file.exists(candidate_shp)) {
  message("Reading candidate dam shapefile: ", candidate_shp)
  
  cand_gdf <- sf::st_read(candidate_shp, quiet = TRUE)
  
  name_col <- find_first_col(cand_gdf, c("name", "dam", "Dam", "DAM", "dam_name", "Name"))
  if (is.na(name_col)) {
    warning("Candidate shapefile found, but no dam-name column was recognised; orange no-HQ points will not be drawn.")
  } else {
    candidate_sites <- cand_gdf |>
      mutate(
        dam = as.character(.data[[name_col]]),
        dam_key = normalise_dam_key(dam)
      ) |>
      group_by(dam_key, dam) |>
      summarise(geometry = sf::st_union(geometry), .groups = "drop") |>
      sf::st_point_on_surface() |>
      sf::st_transform(4326)
    
    coords <- sf::st_coordinates(candidate_sites)
    candidate_sites <- candidate_sites |>
      sf::st_drop_geometry() |>
      mutate(Longitude = coords[, 1], Latitude = coords[, 2])
  }
  
} else if (file.exists(candidate_csv)) {
  message("Reading candidate dam CSV: ", candidate_csv)
  
  cand <- read_csv(candidate_csv, show_col_types = FALSE)
  name_col <- find_first_col(cand, c("dam_preferred", "dam", "name", "Dam", "DAM", "dam_name", "Name"))
  lat_col <- find_first_col(cand, c("Latitude_site", "latitude", "Latitude", "lat", "Lat"))
  lon_col <- find_first_col(cand, c("Longitude_site", "longitude", "Longitude", "lon", "Lon"))
  
  if (any(is.na(c(name_col, lat_col, lon_col)))) {
    warning("Candidate CSV found, but dam/lat/lon columns were not recognised; orange no-HQ points will not be drawn.")
  } else {
    candidate_sites <- cand |>
      transmute(
        dam = as.character(.data[[name_col]]),
        dam_key = normalise_dam_key(dam),
        Latitude = as.numeric(.data[[lat_col]]),
        Longitude = as.numeric(.data[[lon_col]])
      ) |>
      group_by(dam_key) |>
      summarise(
        dam = dplyr::first(na.omit(dam)),
        Latitude = median(Latitude, na.rm = TRUE),
        Longitude = median(Longitude, na.rm = TRUE),
        .groups = "drop"
      ) |>
      filter(is.finite(Latitude), is.finite(Longitude))
  }
}

if (is.null(candidate_sites)) {
  warning(
    "No candidate dam file available. The map will show only dams with retained observations; ",
    "orange no-HQ sites require a candidate shapefile/CSV."
  )
  candidate_sites <- obs_sites |>
    select(dam_key, dam, Latitude, Longitude)
}

map_sites <- candidate_sites |>
  left_join(
    obs_sites |>
      select(dam_key, n_obs, n_lake, n_outlet, n_paired),
    by = "dam_key"
  ) |>
  mutate(
    has_hq = !is.na(n_obs) & n_obs > 0,
    is_selangor = dam_key == normalise_dam_key("Sungai Selangor"),
    map_class = case_when(
      is_selangor ~ "Sungai Selangor validation site",
      has_hq ~ "Retained ECOSTRESS WST observations",
      TRUE ~ "No retained HQ observation"
    )
  )

write_csv(map_sites, file.path(output_dir, "figure1_map_site_observation_counts.csv"))

# -----------------------------------------------------------------------------
# Figure 1: map
# -----------------------------------------------------------------------------

world <- map_data("world") |>
  filter(long >= 88, long <= 130, lat >= -12, lat <= 30)

no_hq_sites <- map_sites |> filter(map_class == "No retained HQ observation")
hq_sites <- map_sites |> filter(map_class == "Retained ECOSTRESS WST observations")
selangor_site <- map_sites |> filter(map_class == "Sungai Selangor validation site")

p_map <- ggplot() +
  geom_polygon(
    data = world,
    aes(x = long, y = lat, group = group),
    fill = "grey92", colour = "grey75", linewidth = 0.2
  ) +
  geom_point(
    data = no_hq_sites,
    aes(x = Longitude, y = Latitude, colour = "No retained HQ observation"),
    shape = 21, fill = "red3", size = 2.0, stroke = 0.8
  ) +
  {
    if (scale_point_size_by_n) {
      geom_point(
        data = hq_sites,
        aes(x = Longitude, y = Latitude, fill = n_obs, size = n_obs),
        shape = 21, colour = "black", stroke = 0.35, alpha = 0.95
      )
    } else {
      geom_point(
        data = hq_sites,
        aes(x = Longitude, y = Latitude, fill = n_obs),
        shape = 21, colour = "black", size = 3.5, stroke = 0.35, alpha = 0.95
      )
    }
  } +
  geom_point(
    data = selangor_site,
    aes(x = Longitude, y = Latitude, colour = "Sungai Selangor validation site"),
    shape = 21,  fill = "green4", size = 3.0, stroke = 0.35
  ) +
  coord_sf(xlim = c(88, 130), ylim = c(-12, 30), expand = FALSE) +
  scale_fill_viridis_c(
    option = "C",
    name = "Retained WST observations",
    breaks = pretty_breaks(n = 5)
  ) +
  scale_colour_manual(
    name = NULL,
    values = c(
      "No retained HQ observation" = "red3",
      "Sungai Selangor validation site" = "green4"
    )
  ) +
  #scale_size_continuous(
    #name = "Retained WST\nobservations",
    #range = c(2.5, 7),
    #breaks = pretty_breaks(n = 4)
  #) +
  labs(
    #title = "ECOSTRESS WST observations across Southeast Asian dams",
    #subtitle = "Retained dams coloured by number of high-quality WST observation times",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    legend.position = "bottom",
    legend.box = "vertical",
    plot.title = element_text(face = "bold"),
    legend.text = element_text(size=15), legend.title = element_text(size=15), 
    legend.key.width = unit(3, "cm"), legend.key.height = unit(1, "cm")
  )
p_map
if (scale_point_size_by_n) {
  p_map <- p_map + guides(
    fill = guide_colourbar(order = 1),
    size = guide_legend(order = 2),
    colour = guide_legend(order = 3, override.aes = list(fill = c("red3", "green4"), size = 6))
  )
} else {
  p_map <- p_map + guides(
    fill = guide_colourbar(order = 1),
    colour = guide_legend(order = 2, override.aes = list(fill = c("red3", "green4"), size = 6))
  )
}
p_map
ggsave(
  file.path(output_dir, "figure1_dam_map_observation_count.png"),
  p_map,
  width = 9.5,
  height = 7.0,
  dpi = 300
)

ggsave(
  file.path(output_dir, "figure1_dam_map_observation_count.pdf"),
  p_map,
  width = 9.5,
  height = 7.0
)

# -----------------------------------------------------------------------------
# Figure 2: Lake / Outlet / Tdiff distributions
# -----------------------------------------------------------------------------

# Long table for lake and outlet WST.
wst_long <- obs |>
  select(dam, dam_key, datetime_utc, Tm_l, Tm_o, T_diff, paired) |>
  pivot_longer(
    cols = c(Tm_l, Tm_o),
    names_to = "variable",
    values_to = "value"
  ) |>
  filter(is.finite(value)) |>
  mutate(
    metric = recode(
      variable,
      Tm_l = "Reservoir WST",
      Tm_o = "Outlet WST"
    ),
    axis_group = "WST"
  ) |>
  select(dam, dam_key, datetime_utc, metric, value, axis_group)

# Tdiff only for paired observations.
tdiff_long <- obs |>
  filter(paired, is.finite(T_diff)) |>
  transmute(
    dam,
    dam_key,
    datetime_utc,
    metric = "Outlet - reservoir WST",
    value = T_diff,
    axis_group = "Tdiff"
  )

temp_long <- bind_rows(wst_long, tdiff_long) |>
  mutate(
    metric = factor(
      metric,
      levels = c("Reservoir WST", "Outlet WST", "Outlet - reservoir WST")
    )
  )

write_csv(temp_long, file.path(output_dir, "figure2_temperature_distribution_data_long.csv"))

summary_temp <- temp_long |>
  group_by(metric) |>
  summarise(
    n = n(),
    n_dams = n_distinct(dam_key),
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    median = median(value, na.rm = TRUE),
    q25 = quantile(value, 0.25, na.rm = TRUE),
    q75 = quantile(value, 0.75, na.rm = TRUE),
    min = min(value, na.rm = TRUE),
    max = max(value, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(summary_temp, file.path(output_dir, "figure2_temperature_distribution_summary.csv"))

# Recommended version: facets, because WST and Tdiff have different physical meaning.
p_temp_facet <- ggplot(temp_long, aes(x = metric, y = value, fill = metric)) +
  geom_hline(
    data = tibble(axis_group = "Tdiff", yint = 0),
    aes(yintercept = yint),
    inherit.aes = FALSE,
    linewidth = 0.4,
    linetype = "dashed",
    colour = "grey35"
  ) +
  geom_boxplot(outlier.alpha = 0.45, width = 0.55) +
  geom_jitter(width = 0.08, alpha = 0.05, size = 0.9) +
  facet_wrap(
    ~ axis_group,
    scales = "free_y",
    #labeller = as_labeller(c(
      #WST = "Water surface temperature (°C)",
      #Tdiff = "Outlet - reservoir WST (°C)"
    #))
  ) +
  scale_fill_manual(
    values = c(
      "Reservoir WST" = "royalblue3",
      "Outlet WST" = "forestgreen",
      "Outlet - reservoir WST" = "grey55"
    )
  ) +
  labs(
    title = "Observed ECOSTRESS water-surface temperature distributions",
    x = NULL,
    y = NULL,
    fill = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )
p_temp_facet
ggsave(
  file.path(output_dir, "figure2_wst_and_tdiff_boxplots_faceted.png"),
  p_temp_facet,
  width = 9.5,
  height = 5.8,
  dpi = 300
)

ggsave(
  file.path(output_dir, "figure2_wst_and_tdiff_boxplots_faceted.pdf"),
  p_temp_facet,
  width = 9.5,
  height = 5.8
)

# Optional version: all three in one panel. This avoids secondary-axis scaling,
# but Tdiff will be visually compressed compared with WST.
p_temp_single <- ggplot(temp_long, aes(x = metric, y = value, fill = metric)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.35, colour = "grey45") +
  geom_boxplot(outlier.alpha = 0.45, width = 0.6) +
  geom_jitter(width = 0.08, alpha = 0.12, size = 0.8) +
  scale_fill_manual(
    values = c(
      "Reservoir WST" = "royalblue3",
      "Outlet WST" = "forestgreen",
      "Outlet - reservoir WST" = "grey55"
    )
  ) +
  labs(
    title = "Observed ECOSTRESS WST and outlet-reservoir temperature difference",
    x = NULL,
    y = "Temperature / temperature difference (°C)",
    fill = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold")
  )
p_temp_single
ggsave(
  file.path(output_dir, "figure2_wst_and_tdiff_boxplots_single_axis.png"),
  p_temp_single,
  width = 8.5,
  height = 5.5,
  dpi = 300
)

# Optional secondary-axis style figure.
# This is provided only for visual comparison. The faceted version is recommended
# for the paper because WST and Tdiff are different quantities.
range_wst <- range(c(obs$Tm_l, obs$Tm_o), na.rm = TRUE)
range_td <- range(obs$T_diff[obs$paired], na.rm = TRUE)
scale_factor <- diff(range_wst) / diff(range_td)
offset <- mean(range_wst) - mean(range_td) * scale_factor

tdiff_scaled <- tdiff_long |>
  mutate(value_scaled = value * scale_factor + offset)
wst_long$metric<-factor(wst_long$metric, levels = c("Outlet WST", "Reservoir WST", "Outlet - reservoir WST"), labels = c("Outlet", "Lake", "Outlet-Lake diff"))
tdiff_scaled$metric<-factor(tdiff_scaled$metric, levels=c("Outlet WST", "Reservoir WST", "Outlet - reservoir WST"), labels = c("Outlet", "Lake","Outlet-Lake diff"))
p_temp_secondary <- ggplot() +
  geom_boxplot(
    data = wst_long,
    aes(x = metric, y = value, fill = metric),
    outlier.alpha = 0.45,
    width = 0.55
  ) +
  geom_boxplot(
    data = tdiff_scaled,
    aes(x = metric, y = value_scaled, fill = metric),
    outlier.alpha = 0.45,
    width = 0.55
  ) +
  scale_y_continuous(
    name = "WST (°C)",
    sec.axis = sec_axis(
      trans = ~ (. - offset) / scale_factor,
      name = " WST diff (°C)"
    )
  ) +
  scale_fill_manual(
    values = c(
      "Lake" = "royalblue3",
      "Outlet" = "forestgreen",
      "Outlet-Lake diff" = "grey55"
    )
  ) +
  
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(), axis.text = element_text(size=12),
    axis.title.y = element_text(size=15)
  ) +
  geom_vline(xintercept = 2.5)
p_temp_secondary
ggsave(
  file.path(output_dir, "figure2_wst_and_tdiff_boxplots_secondary_axis.png"),
  p_temp_secondary,
  width = 9.5,
  height = 5.8,
  dpi = 300
)

message("Wrote outputs to: ", output_dir)
message("Map data: ", file.path(output_dir, "figure1_map_site_observation_counts.csv"))
message("Temperature summary: ", file.path(output_dir, "figure2_temperature_distribution_summary.csv"))
message("Recommended figures:")
message("  ", file.path(output_dir, "figure1_dam_map_observation_count.png"))
message("  ", file.path(output_dir, "figure2_wst_and_tdiff_boxplots_faceted.png"))
