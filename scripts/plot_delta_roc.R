#!/usr/bin/env Rscript
#
# plot_delta_roc.R
#
# Plots precision-recall tradeoff across delta values for the lfc_diff_net filter.
# Reads the full-universe roc_data.full.txt AND the matching restricted-universe
# roc_data.restricted.txt (same directory), and draws them side by side.
#
# Usage:
#   Rscript plot_delta_roc.R <roc_data.full.txt> <out_pdf>
# The restricted file is derived by replacing ".full." with ".restricted." in the
# full path; if it is missing, only the full universe is plotted.

library(ggplot2)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  roc_file <- "posthoc_lfc_filter.betabinom_EBmap.delta0/roc_data.full.txt"
  out_pdf  <- "scripts/plots/delta_roc.pdf"
} else {
  roc_file <- args[1]
  out_pdf  <- args[2]
}

# read one universe's roc table, tag it, keep DTE/DTU
read_universe <- function(path, universe) {
  if (!file.exists(path)) {
    message("Skipping (not found): ", path)
    return(NULL)
  }
  df <- read.table(path, header = TRUE, sep = "\t")
  df <- df[df$sim_type %in% c("DTE", "DTU"), ]
  df$universe <- universe
  df
}

restricted_file <- sub("\\.full\\.", ".restricted.", roc_file)
if (identical(restricted_file, roc_file)) {
  # path did not contain ".full." -- try inserting the universe token
  restricted_file <- sub("roc_data", "roc_data.restricted", roc_file)
}

roc <- bind_rows(
  read_universe(roc_file, "Full universe"),
  read_universe(restricted_file, "Restricted universe")
)
if (is.null(roc) || nrow(roc) == 0) stop("No roc data could be read.")

# order universe facets Full then Restricted
roc$universe <- factor(roc$universe,
                       levels = c("Full universe", "Restricted universe"))

# Replace -Inf delta with a label for plotting
roc$delta_num <- roc$delta
roc$delta_lab <- ifelse(is.infinite(roc$delta), "-Inf (no filter)",
                        as.character(roc$delta))

# Identify key points
roc$point_type <- "other"
roc$point_type[is.infinite(roc$delta)] <- "baseline"
roc$point_type[roc$delta == 0] <- "delta0"

# Order for line drawing (within each universe x sim_type)
roc <- roc[order(roc$universe, roc$sim_type, roc$delta_num), ]

label_rows <- roc[roc$delta %in% c(-Inf, 0, 0.5, 1, 2) &
                  roc$sim_type %in% c("DTE", "DTU"), ]
label_rows$label <- ifelse(is.infinite(label_rows$delta), "baseline",
                    paste0("delta=", label_rows$delta))

p <- ggplot(roc, aes(x = micro_recall, y = micro_precision,
                     color = sim_type, group = sim_type)) +
  geom_line(linewidth = 0.8) +
  geom_point(data = roc[roc$point_type == "baseline", ],
             shape = 1, size = 4, stroke = 1.2) +
  geom_point(data = roc[roc$point_type == "delta0", ],
             shape = 18, size = 4) +
  geom_point(data = roc[roc$point_type == "other", ],
             shape = 16, size = 2) +
  facet_wrap(~ universe) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1)) +
  scale_color_manual(values = c("DTE" = "#d73027", "DTU" = "#4575b4"),
                     name = NULL) +
  labs(
    title = "lfc_diff_net filter: precision-recall tradeoff",
    subtitle = "Sweep of delta threshold on lfc_diff_net > delta  |  padj < 0.01",
    x = "Recall (micro)",
    y = "Precision (micro)",
    caption = "Diamond = delta=0  Circle = baseline (no filter)"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "top")

# Try ggrepel; fall back to geom_text if not installed
have_repel <- requireNamespace("ggrepel", quietly = TRUE)
if (have_repel) {
  p <- p + ggrepel::geom_text_repel(data = label_rows,
                                    aes(label = label),
                                    size = 3, show.legend = FALSE,
                                    min.segment.length = 0.2,
                                    box.padding = 0.3)
} else {
  p <- p + geom_text(data = label_rows, aes(label = label),
                     size = 3, hjust = -0.1, show.legend = FALSE)
}

n_uni <- length(unique(roc$universe))
ggsave(out_pdf, p, width = if (n_uni > 1) 10 else 6, height = 5)
message("Written: ", out_pdf, "  (", n_uni, " universe(s))")
