# Shared plot palette and display labels for the benchmark figures.
#
# Sourced by plot_roc_partial.R and plot_fp_universe_shift.R.
# Defined ONCE here because these three scripts drifted apart three times:
# DEXSeq was blue in one and orange in two, GrASE_BH had a colour in one and was
# missing from every TOOLS list, and the disp() relabel reached only one script.
# Any new benchmark figure should source this rather than copy the vectors.
#
# Conventions:
#   GrASE       reds, darkening with the effect-size gate (dpi0 -> dpi0.2)
#   GrASE_BH    yellow-orange -- a GrASE variant, but an FDR-architecture change
#               (plain BH) rather than an effect-size gate, so off the red ramp
#   DEXSeq      blue
#   MAJIQ       greens, darkening with confidence threshold
#   rMATS       purples, light -> dark with the dPSI gate
COL <- c("GrASE_BH"            = "#fec44f",
         "GrASE"               = "#fb6a4a",
         "GrASE_merged_all"    = "#fb6a4a",
         "GrASE_dpi0.1"        = "#de2d26",
         "GrASE_merged_dpi0.1" = "#de2d26",
         "GrASE_dpi0.2"        = "#a50f15",
         "GrASE_merged_dpi0.2" = "#a50f15",
         "GrASE_internal"        = "#67000d",
         "GrASE_merged_internal" = "#67000d",
         "DEXSeq"              = "#3182bd",
         "MAJIQ_C0.10"         = "#a1d99b",
         "MAJIQ_C0.20"         = "#238b45",
         "rMATS"               = "#bcbddc",
         "rMATS_dpsi0.1"       = "#756bb1",
         "rMATS_dpsi0.2"       = "#3f007d")

# Merged is the only GrASE family plotted, so the "_merged" infix carries no
# information in a legend. The ungated run is GrASE_dpi0, never plain "GrASE":
# exontest.R ships --min_dpi 0.1, so the unsuffixed name would misrepresent the
# default. Presentation only -- tool names in the written tables are unchanged.
disp <- function(x)
  sub("^GrASE_merged_", "GrASE_", sub("^GrASE_merged_all$", "GrASE_dpi0", x))

# Legend label: GrASE_BH spelled out, everything else via disp().
plab <- function(x) ifelse(x == "GrASE_BH", "GrASE (plain BH)", disp(x))
