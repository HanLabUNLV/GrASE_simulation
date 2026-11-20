# Shared plot palette and display labels for the benchmark figures.
#
# Sourced by plot_roc_partial.R and plot_fp_universe_shift.R.
# Defined ONCE here because these three scripts drifted apart three times:
# DEXSeq was blue in one and orange in two, GrASE_BH had a colour in one and was
# missing from every TOOLS list, and the disp() relabel reached only one script.
# Any new benchmark figure should source this rather than copy the vectors.
#
# Conventions:
#   GrASE       reds, darkening with the delta_pi gate (dpi0 -> dpi0.2)
#   GrASE_eff   pinks, darkening with the |effect_size| gate (eff0.3 -> eff1.0)
#   GrASE_merged_default  black -- exontest.R's default call rule (both gates)
#   GrASE_nestedBH  yellow-orange -- a GrASE variant, but an FDR-architecture change
#               (nested BH; the other GrASE entries use plain BH) rather than an
#               effect-size gate, so off the red ramp
#   DEXSeq      blue
#   MAJIQ       greens, darkening with confidence threshold
#   rMATS       purples, light -> dark with the dPSI gate
COL <- c("GrASE_merged_default" = "#000000",
         "GrASE_nestedBH"      = "#fec44f",
         "GrASE"               = "#fb6a4a",
         "GrASE_merged_all"    = "#fb6a4a",
         "GrASE_dpi0.1"        = "#de2d26",
         "GrASE_merged_dpi0.1" = "#de2d26",
         "GrASE_dpi0.2"        = "#a50f15",
         "GrASE_merged_dpi0.2" = "#a50f15",
         "GrASE_merged_eff0.3" = "#fcc5c0",
         "GrASE_merged_eff0.4" = "#fa9fb5",
         "GrASE_merged_eff0.5" = "#f768a1",
         "GrASE_merged_eff1.0" = "#ae017e",
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

# Legend label: GrASE_nestedBH spelled out, everything else via disp().
plab <- function(x) ifelse(x == "GrASE_nestedBH", "GrASE (nested BH)",
                     ifelse(x == "GrASE_merged_default", "GrASE (default)", disp(x)))

# Main-text figures show GrASE as one |effect_size| sweep on plain BH; the nested
# BH and delta_pi curves go to the supplementary copy (*.supp.pdf). The sweep is
# eff0 (no effect-size gate) -> eff0.3 -> eff0.4 -> eff0.5 (exontest.R's default) -> eff1.0, on a
# red ramp that darkens with the gate. GrASE_merged_eff0.5 is the same call set
# as GrASE_merged_default and is not drawn twice.
MAIN_TOOLS <- c("GrASE_merged_all", "GrASE_merged_eff0.3", "GrASE_merged_eff0.4",
                "GrASE_merged_default", "GrASE_merged_eff1.0",
                "DEXSeq", "MAJIQ_C0.10", "MAJIQ_C0.20",
                "rMATS", "rMATS_dpsi0.1", "rMATS_dpsi0.2")
COL_MAIN <- COL
COL_MAIN[c("GrASE_merged_all", "GrASE_merged_eff0.3", "GrASE_merged_eff0.4",
           "GrASE_merged_default", "GrASE_merged_eff1.0")] <-
  c("#fcbba1", "#fc9272", "#fb6a4a", "#de2d26", "#67000d")
plab_main <- function(x) ifelse(x == "GrASE_merged_all", "GrASE_eff0",
                          ifelse(x == "GrASE_merged_default", "GrASE_eff0.5 (default)", plab(x)))
