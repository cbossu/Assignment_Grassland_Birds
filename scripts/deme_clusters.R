#!/usr/bin/env Rscript
# deme_clusters.R — label each FEEMS deme with its genoscape cluster
#
# Usage:
#   Rscript deme_clusters.R <genoscape.tif> <nonbreeding_assignment_long.csv> <out.csv> \
#                           [buffer_km=100] [color_names.csv]
#
# Each deme gets the most common genoscape colour within buffer_km of its centre
# (deme centres often fall in range gaps, so a single pixel is not enough).
# Colours are matched to the main genoscape colours; anti-aliased edge pixels are
# snapped to the nearest main colour. Demes with no genoscape pixels get "none".
#
# color_names.csv (optional): two columns, hex,cluster  e.g.
#   hex,cluster
#   #E69F00,West
#   #009E73,Central
#   #56B4E9,East
# Run once without it: the script prints the colours it found so you can name them.

suppressPackageStartupMessages({ library(terra); library(data.table) })

a <- commandArgs(trailingOnly = TRUE)
if (length(a) < 3) stop("Usage: Rscript deme_clusters.R <genoscape.tif> <long.csv> <out.csv> [buffer_km] [color_names.csv]")
f_gs <- a[1]; f_long <- a[2]; f_out <- a[3]
buf_km  <- if (length(a) >= 4) as.numeric(a[4]) else 100
f_names <- if (length(a) >= 5) a[5] else NA

# ---- genoscape pixels -> hex colours ------------------------------------------
g <- rast(f_gs)
if (nlyr(g) >= 4) g <- mask(g[[1:3]], g[[4]] > 0, maskvalues = FALSE)
g <- g[[1:3]]
mx <- max(minmax(g)[2, ], na.rm = TRUE)
scale255 <- if (mx <= 1) 255 else 1

v <- values(g, mat = TRUE)
v <- v[complete.cases(v), , drop = FALSE]
hex_all <- rgb(v[, 1] * scale255, v[, 2] * scale255, v[, 3] * scale255, maxColorValue = 255)
tab <- sort(table(hex_all), decreasing = TRUE)
main <- names(tab)[tab / sum(tab) >= 0.01]          # main colours: >= 1% of range pixels
# drop white/near-white backgrounds if the raster stores them instead of NA
main <- main[colSums(col2rgb(main)) < 3 * 245]
cat("[colours] main genoscape colours (share of occupied pixels):\n")
print(round(100 * tab[main] / sum(tab), 1))

main_rgb <- t(col2rgb(main))
snap <- function(hex) {                             # nearest main colour, per unique hex
  u <- unique(hex)
  if (length(u) == 0) return(character(0))
  x <- t(col2rgb(u))
  nearest <- apply(x, 1, function(p) main[which.min(colSums((t(main_rgb) - p)^2))])
  setNames(nearest, u)[hex]
}

# ---- demes ---------------------------------------------------------------------
long  <- fread(f_long)
demes <- unique(long[, .(deme, deme_lon, deme_lat)])
pts   <- vect(demes, geom = c("deme_lon", "deme_lat"), crs = "EPSG:4326")
if (!same.crs(pts, g)) pts <- project(pts, crs(g))
bufs  <- buffer(pts, width = buf_km * 1000)

ex <- extract(g, bufs)                               # all cells in each buffer, with ID
setDT(ex)
ex <- ex[complete.cases(ex)]
names(ex)[2:4] <- c("r", "g", "b")
ex[, hex := snap(rgb(r * scale255, g * scale255, b * scale255, maxColorValue = 255))]
lab <- ex[, .N, by = .(ID, hex)][order(ID, -N)][, .(hex = hex[1], share = N[1] / sum(N)), by = ID]

demes[, ID := .I]
demes[lab, `:=`(hex = i.hex, share = i.share), on = "ID"]
demes[is.na(hex), `:=`(hex = "none", share = NA_real_)]

if (!is.na(f_names)) {
  nm <- fread(f_names, strip.white = TRUE)
  nm[, `:=`(hex = toupper(trimws(gsub("\r", "", hex))), cluster = trimws(gsub("\r", "", cluster)))]
  nm[!startsWith(hex, "#"), hex := paste0("#", hex)]
  demes[, hex := toupper(hex)]
  miss <- setdiff(toupper(main), nm$hex)
  if (length(miss)) warning("colours with no name in ", f_names, ": ", paste(miss, collapse = ", "),
                            "  (names file has: ", paste(nm$hex, collapse = ", "), ")")
  demes[nm, cluster := i.cluster, on = "hex"]
  demes[is.na(cluster), cluster := ifelse(hex == "NONE", "none", hex)]
} else {
  demes[, cluster := hex]
}
demes[, ID := NULL]
fwrite(demes, f_out)
cat(sprintf("[demes] %d demes labelled (%d with no genoscape pixels within %g km)\n",
            nrow(demes), sum(demes$cluster == "none"), buf_km))
print(demes[, .N, by = cluster][order(-N)])
