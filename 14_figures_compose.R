# Compose Figures 1 and 2 (2 x 2 panels).
# Arranges the four panel PNGs already produced by 04_model_patient_level.R and
# 05_model_aggregate_level.R into single two-by-two figures. No analysis and no model
# is re-run, so no number can change; output is written at 300 dpi.

# Mapping (Figure 1 = Analysis 1, Figure 2 = Analysis 2):
#   aggregate_level_* = Analysis 1, hospital-wide aggregate, species-proxy outcome -> Figure 1
#   patient_level_* = Analysis 2, individual-level prescription only, AST-confirmed MDR -> Figure 2
# The aggregate_level_* panels come from 05_model_aggregate_level.R, the patient_level_* panels from 04_model_patient_level.R.
# This script only re-arranges the bitmaps and writes 300 dpi metadata.

suppressPackageStartupMessages({
  library(png)
  library(grid)
})
dir.create("output/ml", showWarnings = FALSE, recursive = TRUE)

# Draw one panel: the PNG is scaled to preserve its aspect ratio and centred in the viewport.
# When tag_band > 0 the image is shrunk by (1 - tag_band) and placed in the lower part
# of the viewport, leaving the top band free for the A/B tag so the tag cannot sit on
# the panel's own title or axis labels.
draw_panel <- function(file, tag, tag_band = 0, fontsize = 16, tag_y = 2) {
  img <- readPNG(file)
  vp <- current.viewport()
  vp_w <- convertWidth(unit(1, "npc"), "npc", valueOnly = TRUE)
  vp_h <- convertHeight(unit(1, "npc"), "npc", valueOnly = TRUE)
  ar_img <- dim(img)[2] / dim(img)[1]        # width / height
  ar_vp <- vp_w / vp_h
  if (ar_img > ar_vp) {
    w <- 1; h <- ar_vp / ar_img
  } else {
    h <- 1; w <- ar_img / ar_vp
  }
  if (tag_band > 0) {
    s <- 1 - tag_band
    w <- w * s
    h <- h * s
    y_c <- s / 2                          # image in lower [0, s]; top band [s, 1] reserved
  } else {
    y_c <- 0.5
  }
  grid.raster(img, x = 0.5, y = y_c, width = unit(w, "npc"), height = unit(h, "npc"))
  grid.text(tag, x = unit(2, "mm"), y = unit(1, "npc") - unit(tag_y, "mm"),
            just = c("left", "top"),
            gp = gpar(fontsize = fontsize, fontface = "bold"))
}

compose <- function(files, tags, out_png, out_pdf, width_in = 7, height_in = 7, res = 300,
                    tag_band = 0, fontsize = 16, tag_y = 2) {
  stopifnot(length(files) == 4, all(file.exists(files)))
  pad <- 0.015
  gw <- (1 - 3 * pad) / 2
  gh <- (1 - 3 * pad) / 2

# viewport origin is bottom-left: A top-left, B top-right, C bottom-left, D bottom-right
  slots <- list(
    c(pad,               1 - pad - gh),
    c(2 * pad + gw,      1 - pad - gh),
    c(pad,               pad),
    c(2 * pad + gw,      pad)
  )

  png(out_png, width = width_in * res, height = height_in * res, res = res, bg = "white")
  grid.newpage()
  for (i in seq_along(files)) {
    pushViewport(viewport(x = slots[[i]][1], y = slots[[i]][2],
                          width = gw, height = gh, just = c("left", "bottom")))
    draw_panel(files[i], tags[i], tag_band, fontsize, tag_y)
    popViewport()
  }
  dev.off()

  pdf(out_pdf, width = width_in, height = height_in)
  grid.newpage()
  for (i in seq_along(files)) {
    pushViewport(viewport(x = slots[[i]][1], y = slots[[i]][2],
                          width = gw, height = gh, just = c("left", "bottom")))
    draw_panel(files[i], tags[i], tag_band, fontsize, tag_y)
    popViewport()
  }
  dev.off()

  cat(sprintf("%s  %dx%d px @ %d dpi\n", basename(out_png),
              width_in * res, height_in * res, res))
}

base_ml <- "output/ml"

# Figure 1 - Analysis 1 (hospital-wide, species-proxy outcome) -> aggregate_level_*
# tag_band = 0.09 reserves 9% at the top for the A/B/C/D tags so they clear the panel titles
# font size 14, slightly smaller than Figure 2's 16, to fit the tag band
compose(
  files = file.path(base_ml, c("aggregate_level_roc.png", "aggregate_level_shap.png",
                               "aggregate_level_calibration.png", "aggregate_level_lasso_cv.png")),
  tags = c("A", "B", "C", "D"),
  out_png = file.path(base_ml, "figure1_analysis1.png"),
  out_pdf = file.path(base_ml, "figure1_analysis1.pdf"),
  tag_band = 0.09,
  fontsize = 14
)

# Figure 2 - Analysis 2 (individual-patient, AST-confirmed MDR) -> patient_level_*
# tag_band = 0.09 reserves the top band so the tag clears the panel C title
# font size 14, matching Figure 1
compose(
  files = file.path(base_ml, c("patient_level_roc.png", "patient_level_shap.png",
                               "patient_level_calibration.png", "patient_level_lasso_cv.png")),
  tags = c("A", "B", "C", "D"),
  out_png = file.path(base_ml, "figure2_analysis2.png"),
  out_pdf = file.path(base_ml, "figure2_analysis2.pdf"),
  tag_band = 0.09,
  fontsize = 14
)
