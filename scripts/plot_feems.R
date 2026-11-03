#!/usr/bin/env Rscript
# plot_feems.R — figures for feems_assign_nonbreeding.py output (LOSH)
#
# Usage:
#   Rscript plot_feems.R <run_dir> [calib_dir]
#
#   run_dir    tempered run (e.g. feems_LOSH_fst_T251): nonbreeding_assignment_long.csv,
#              nonbreeding_assignment_summary.csv, cv_error.csv, run_info.json
#   calib_dir  run that used --loo (e.g. feems_LOSH_fst): loo_validation.csv,
#              loo_temperature_calibration.csv. Defaults to run_dir.
#
# Figures go to <run_dir>/figures/:
#   fig1_cv_surface.png             CV error over lambda x lambda_q
#   fig2_temperature_calibration.png  LOO 95%-set coverage vs temperature T
#   fig3_loo_error.png              LOO error histogram + true -> MAP map
#   fig4_origin_density.png         mean posterior over all nonbreeding birds
#   fig5_assignment_certainty.png   MAP prob and 95%-set size, sampled vs unsampled MAP
#   fig6_individual_posteriors.pdf  one page per nonbreeding bird
#   nonbreeding_95cs_R.csv          per-bird 95% set recomputed here (for checking)
#
# Packages: data.table, ggplot2, sf, patchwork; optional rnaturalearth (+ rnaturalearthdata)
# or maps for the basemap, jsonlite for run_info.json.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(sf)
  library(patchwork)
})
sf_use_s2(FALSE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript plot_feems.R <run_dir> [calib_dir]")
run_dir   <- args[1]
calib_dir <- if (length(args) >= 2) args[2] else run_dir
fig_dir   <- file.path(run_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# North America Albers equal-area; edit lon_0 / lat_0 if your range is centred elsewhere
CRS_PLOT <- "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-100 +datum=WGS84 +units=m +no_defs"

find_file <- function(name) {
  for (d in unique(c(run_dir, calib_dir))) {
    f <- file.path(d, name)
    if (file.exists(f)) return(f)
  }
  NA_character_
}
pick_col <- function(dt, choices) {
  hit <- intersect(choices, names(dt))
  if (length(hit)) hit[1] else NA_character_
}
save_fig <- function(p, name, w = 8, h = 6) {
  ggsave(file.path(fig_dir, name), p, width = w, height = h, dpi = 300)
  message("[fig] ", name)
}
theme_set(theme_bw(base_size = 11) +
            theme(panel.grid.minor = element_blank(),
                  strip.background = element_rect(fill = "grey95", colour = NA)))

# ---- read assignment output --------------------------------------------------
f_long <- find_file("nonbreeding_assignment_long.csv")
f_summ <- find_file("nonbreeding_assignment_summary.csv")
if (is.na(f_long)) stop("nonbreeding_assignment_long.csv not found in ", run_dir)
long <- fread(f_long)
summ <- if (!is.na(f_summ)) fread(f_summ) else NULL
long[, deme := as.character(deme)]

T_used <- NA_real_
f_info <- find_file("run_info.json")
if (!is.na(f_info) && requireNamespace("jsonlite", quietly = TRUE)) {
  info <- jsonlite::fromJSON(f_info)
  if (!is.null(info$temp)) T_used <- info$temp
}

# ---- deme table + sampled flag ----------------------------------------------
demes <- unique(long[, .(deme, deme_lon, deme_lat)])
samp_col <- pick_col(long, c("sampled", "is_sampled", "n_breeding", "n_samples", "n_breed"))
f_loo <- find_file("loo_validation.csv")
loo   <- if (!is.na(f_loo)) fread(f_loo) else NULL
if (!is.na(samp_col)) {
  s <- unique(long[, c("deme", samp_col), with = FALSE])
  demes[s, sampled := as.numeric(get(samp_col)) > 0, on = "deme"]
} else if (!is.null(loo)) {
  # every deme held out in leave-one-deme-out had breeding samples
  demes[, sampled := deme %in% as.character(loo$true_deme)]
} else {
  demes[, sampled := NA]
}
message(sprintf("[demes] %d demes, %s sampled", nrow(demes),
                ifelse(all(is.na(demes$sampled)), "unknown", sum(demes$sampled, na.rm = TRUE))))

# ---- spatial objects: deme cells (Voronoi, clipped) + basemap ---------------
deme_pts <- st_transform(st_as_sf(demes, coords = c("deme_lon", "deme_lat"), crs = 4326), CRS_PLOT)
xy <- st_coordinates(deme_pts)
nn <- apply(as.matrix(dist(xy)), 1, function(r) min(r[r > 0]))
spacing_m <- median(nn)
message(sprintf("[grid] median deme spacing %.0f km", spacing_m / 1000))

vor <- st_collection_extract(st_voronoi(st_union(deme_pts)), "POLYGON")
clip <- st_union(st_buffer(deme_pts, spacing_m * 0.6))
cells <- st_sf(geometry = st_intersection(st_sfc(vor, crs = CRS_PLOT), clip))
cells <- st_join(cells, deme_pts, join = st_contains)
cells <- cells[!is.na(cells$deme), ]

bbox <- st_bbox(st_buffer(st_union(deme_pts), spacing_m))
get_basemap <- function() {
  if (requireNamespace("rnaturalearth", quietly = TRUE) &&
      requireNamespace("rnaturalearthdata", quietly = TRUE)) {
    w <- rnaturalearth::ne_countries(scale = 50, returnclass = "sf")
  } else if (requireNamespace("maps", quietly = TRUE)) {
    w <- st_as_sf(maps::map("world", plot = FALSE, fill = TRUE))
  } else return(NULL)
  w <- st_sf(geometry = st_make_valid(st_transform(st_geometry(w), CRS_PLOT)))
  st_crop(w, bbox)
}
get_states <- function() {
  if (!requireNamespace("maps", quietly = TRUE)) return(NULL)
  s <- st_sf(geometry = st_make_valid(st_transform(st_geometry(
         st_as_sf(maps::map("state", plot = FALSE, fill = TRUE))), CRS_PLOT)))
  st_crop(s, bbox)
}
world  <- suppressWarnings(get_basemap())
states <- suppressWarnings(get_states())

base_map <- function() {
  g <- ggplot()
  if (!is.null(world))  g <- g + geom_sf(data = world, fill = "grey92", colour = "grey60", linewidth = 0.25)
  if (!is.null(states)) g <- g + geom_sf(data = states, fill = NA, colour = "grey75", linewidth = 0.15)
  g
}
finish_map <- function(g) {
  g + coord_sf(crs = CRS_PLOT, xlim = bbox[c("xmin", "xmax")], ylim = bbox[c("ymin", "ymax")],
               expand = FALSE) +
    theme(axis.title = element_blank())
}
sampled_pts <- deme_pts[deme_pts$deme %in% demes[sampled == TRUE, deme], ]

# ---- per-bird 95% credible set (recomputed; same rule as the LOO) -----------
setorder(long, sample_id, -post_prob)
long[, better_mass := cumsum(post_prob) - post_prob, by = sample_id]
long[, in95 := better_mass < 0.95]
cs <- long[, .(map_deme = deme[1], map_prob = post_prob[1],
               map_lon = deme_lon[1], map_lat = deme_lat[1],
               n_demes_95cs = sum(in95)), by = sample_id]
cs[demes, map_sampled := i.sampled, on = c(map_deme = "deme")]

nb_lon <- if (!is.null(summ)) pick_col(summ, c("lon", "sample_lon", "nb_lon", "longitude")) else NA
nb_lat <- if (!is.null(summ)) pick_col(summ, c("lat", "sample_lat", "nb_lat", "latitude")) else NA
if (!is.na(nb_lon) && !is.na(nb_lat)) {
  cs[summ, `:=`(nb_lon = get(paste0("i.", nb_lon)), nb_lat = get(paste0("i.", nb_lat))),
     on = "sample_id"]
} else {
  cs[, `:=`(nb_lon = NA_real_, nb_lat = NA_real_)]
}
fwrite(cs, file.path(fig_dir, "nonbreeding_95cs_R.csv"))
if (!is.null(summ) && "n_demes_95cs" %in% names(summ)) {
  chk <- merge(cs[, .(sample_id, r = n_demes_95cs)],
               summ[, .(sample_id, py = n_demes_95cs)], by = "sample_id")
  message(sprintf("[check] 95%% set size matches python summary for %.0f%% of birds",
                  100 * mean(chk$r == chk$py)))
}

# ---- fig1: CV surface --------------------------------------------------------
f_cv <- find_file("cv_error.csv")
if (!is.na(f_cv)) {
  cv <- fread(f_cv, header = TRUE)  # numeric header, so force it
  setnames(cv, 1, "lamb_q")
  cvl <- melt(cv, id.vars = "lamb_q", variable.name = "lamb", value.name = "cv_err")
  cvl[, lamb := as.numeric(as.character(lamb))]
  cvl[, pct_above_min := 100 * (cv_err / min(cv_err, na.rm = TRUE) - 1)]
  best <- cvl[which.min(cv_err)]
  p1 <- ggplot(cvl, aes(log10(lamb), log10(lamb_q), fill = pct_above_min)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.1f", pct_above_min)), size = 2.6) +
    geom_point(data = best, shape = 4, size = 5, stroke = 1.2, colour = "black") +
    scale_fill_viridis_c(option = "magma", direction = -1, name = "% above\nmin CV error") +
    labs(x = expression(log[10]~lambda), y = expression(log[10]~lambda[q]),
         title = "FEEMS cross-validation error",
         subtitle = sprintf("min at lambda = %.3g, lambda_q = %.3g (x); values = %% above minimum",
                            best$lamb, best$lamb_q)) +
    coord_fixed()
  save_fig(p1, "fig1_cv_surface.png", 8, 6)
}

# ---- fig2: temperature calibration ------------------------------------------
f_cal <- find_file("loo_temperature_calibration.csv")
if (!is.na(f_cal)) {
  cal <- fread(f_cal)
  T_cal <- cal[coverage >= 0.95, min(T)]
  p2 <- ggplot(cal, aes(T, coverage)) +
    geom_hline(yintercept = 0.95, linetype = 2, colour = "grey40") +
    geom_vline(xintercept = T_cal, linetype = 3, colour = "firebrick") +
    geom_line() + geom_point(size = 1.8) +
    annotate("text", x = T_cal, y = 0.05, hjust = -0.1, colour = "firebrick",
             label = sprintf("T = %.0f", T_cal)) +
    scale_x_log10() + scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    labs(x = "Temperature T (log scale)", y = "True deme inside 95% set",
         title = "Leave-one-deme-out calibration of assignment posteriors")
  save_fig(p2, "fig2_temperature_calibration.png", 7, 4.5)
}

# ---- fig3: LOO error ---------------------------------------------------------
if (!is.null(loo)) {
  med <- median(loo$err_km)
  p3a <- ggplot(loo, aes(err_km)) +
    geom_histogram(binwidth = spacing_m / 1000 / 2, boundary = 0,
                   fill = "grey70", colour = "white") +
    geom_vline(xintercept = med, colour = "firebrick") +
    annotate("text", x = med, y = Inf, vjust = 1.5, hjust = -0.1, colour = "firebrick",
             label = sprintf("median %.0f km\n(%.1f deme spacings)", med, med / (spacing_m / 1000))) +
    labs(x = "Distance, sample location to MAP deme (km)", y = "Breeding birds",
         title = "Leave-one-deme-out error")
  seg <- copy(loo)
  t_xy <- st_coordinates(st_transform(st_as_sf(seg, coords = c("true_lon", "true_lat"), crs = 4326), CRS_PLOT))
  m_xy <- st_coordinates(st_transform(st_as_sf(seg, coords = c("map_lon", "map_lat"), crs = 4326), CRS_PLOT))
  seg[, c("x0", "y0", "x1", "y1") := .(t_xy[, 1], t_xy[, 2], m_xy[, 1], m_xy[, 2])]
  p3b <- finish_map(
    base_map() +
      geom_segment(data = seg, aes(x0, y0, xend = x1, yend = y1, colour = err_km),
                   linewidth = 0.4, alpha = 0.7,
                   arrow = arrow(length = unit(1.2, "mm"), type = "closed")) +
      geom_point(data = seg, aes(x0, y0), size = 0.8) +
      scale_colour_viridis_c(option = "plasma", name = "Error (km)") +
      labs(title = "Sample location -> MAP deme"))
  save_fig(p3a + p3b + plot_layout(widths = c(1, 1.3)), "fig3_loo_error.png", 13, 5.5)
}

# ---- fig4: aggregate origin density -----------------------------------------
dens <- long[, .SD[1], by = sample_id][, .(n_map = .N), by = .(deme)][  # long is sorted by -post_prob
  long[, .(mean_post = mean(post_prob)), by = deme], on = "deme"]
dens[is.na(n_map), n_map := 0]
cells_d <- merge(cells, dens, by = "deme")
ttl <- if (!is.na(T_used)) sprintf(" (T = %.0f)", T_used) else ""
p4 <- base_map() +
  geom_sf(data = cells_d, aes(fill = mean_post), colour = "white", linewidth = 0.1, alpha = 0.9) +
  scale_fill_viridis_c(option = "viridis", name = "Mean posterior\nacross birds") +
  geom_sf(data = sampled_pts, shape = 21, fill = NA, colour = "black", size = 1.4, stroke = 0.4)
if (!all(is.na(cs$nb_lon))) {
  nbp <- st_transform(st_as_sf(cs[!is.na(nb_lon)], coords = c("nb_lon", "nb_lat"), crs = 4326), CRS_PLOT)
  p4 <- p4 + geom_sf(data = nbp, shape = 17, colour = "firebrick", size = 1.3, alpha = 0.7)
}
p4 <- finish_map(p4) +
  labs(title = paste0("Estimated breeding origins of nonbreeding LOSH", ttl),
       subtitle = sprintf("%d birds; open circles = sampled breeding demes; red triangles = nonbreeding sample sites",
                          uniqueN(long$sample_id)))
save_fig(p4, "fig4_origin_density.png", 9, 7)

# ---- fig5: certainty ---------------------------------------------------------
cs[, map_type := fifelse(is.na(map_sampled), "unknown",
                          fifelse(map_sampled, "MAP deme sampled", "MAP deme unsampled"))]
p5a <- ggplot(cs, aes(map_prob, fill = map_type)) +
  geom_histogram(bins = 25, boundary = 0, colour = "white", position = "stack") +
  scale_fill_manual(values = c("MAP deme sampled" = "#2c7fb8", "MAP deme unsampled" = "#fdae61",
                               unknown = "grey60"), name = NULL) +
  labs(x = "Posterior probability of MAP deme", y = "Nonbreeding birds")
p5b <- ggplot(cs, aes(n_demes_95cs, fill = map_type)) +
  geom_histogram(binwidth = 1, boundary = 0.5, colour = "white", position = "stack") +
  scale_fill_manual(values = c("MAP deme sampled" = "#2c7fb8", "MAP deme unsampled" = "#fdae61",
                               unknown = "grey60"), name = NULL) +
  labs(x = "Demes in 95% credible set", y = "Nonbreeding birds")
p5 <- (p5a + p5b) + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
p5 <- p5 + plot_annotation(title = paste0("Assignment certainty", ttl),
                           subtitle = sprintf("median MAP prob %.2f; median 95%% set %d demes (~%.0f km grid)",
                                              median(cs$map_prob), as.integer(median(cs$n_demes_95cs)),
                                              spacing_m / 1000))
save_fig(p5, "fig5_assignment_certainty.png", 10, 4.5)

# ---- fig6: one page per bird -------------------------------------------------
pdf(file.path(fig_dir, "fig6_individual_posteriors.pdf"), width = 7, height = 6)
for (id in cs[order(-map_prob), sample_id]) {
  b <- long[sample_id == id]
  cb <- merge(cells, b[, .(deme, post_prob, in95)], by = "deme")
  r <- cs[sample_id == id]
  g <- base_map() +
    geom_sf(data = cb, aes(fill = post_prob), colour = NA, alpha = 0.9) +
    geom_sf(data = st_union(cb[cb$in95, ]), fill = NA, colour = "black", linewidth = 0.6) +
    geom_sf(data = sampled_pts, shape = 21, fill = NA, colour = "grey30", size = 1, stroke = 0.3) +
    scale_fill_viridis_c(option = "rocket", direction = -1, name = "Posterior",
                         limits = c(0, max(b$post_prob)))
  if (!is.na(r$nb_lon)) {
    nb1 <- st_transform(st_as_sf(r, coords = c("nb_lon", "nb_lat"), crs = 4326), CRS_PLOT)
    g <- g + geom_sf(data = nb1, shape = 24, fill = "deepskyblue", colour = "black", size = 2.8)
  }
  g <- finish_map(g) +
    labs(title = id,
         subtitle = sprintf("MAP prob %.2f; 95%% set = %d demes (outlined); MAP deme %s%s",
                            r$map_prob, r$n_demes_95cs,
                            ifelse(isTRUE(r$map_sampled), "sampled", "unsampled"), ttl))
  print(g)
}
invisible(dev.off())
message("[fig] fig6_individual_posteriors.pdf (", nrow(cs), " pages)")
message("[done] figures in ", fig_dir)
