library(patchwork)
library(ggplot2)
library(cowplot)
library(rsvg)
library(grid)
library(magick)
library(ComplexHeatmap)
load("../240328_RNAseq2nd/analysis/figure_RNAseq.rda")
load("../240415_baQFA/integrated_analysis/figure_integrated_analysis.rda")
heatmap_grob <- grid.grabExpr(draw(heatmap, newpage=FALSE), )
heatmap_gg <- ggdraw() + draw_grob(heatmap_grob)
fig1a_svg <- image_read_svg("./fig1a.svg")
fig1a_grob <- rasterGrob(image = as.raster(fig1a_svg))
fig1a <- ggdraw() +
  draw_grob(fig1a_grob)
(fig1a + plot_integrated) / (mds_plot + heatmap_gg)


