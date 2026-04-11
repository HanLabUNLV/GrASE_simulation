#!/usr/bin/env Rscript
#
# plot_delta_roc.R
#
# Plots precision-recall tradeoff across delta values for the lfc_diff_net filter.
# Reads roc_data.restricted.txt from the posthoc_lfc_filter output directory.
#
# Usage:
#   Rscript plot_delta_roc.R <roc_data_restricted.txt> <out_pdf>

library(ggplot2)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  roc_file <- "posthoc_lfc_filter.delta0/roc_data.restricted.txt"
  out_pdf  <- "scripts/plots/delta_roc.pdf"
} else {
  roc_file <- args[1]
  out_pdf  <- args[2]
}

roc <- read.table(roc_file, header = TRUE, sep = "\t")

# Keep only DTE and DTU
roc <- roc[roc$sim_type %in% c("DTE", "DTU"), ]

# Replace -Inf delta with a label for plotting
roc$delta_num <- roc$delta
roc$delta_lab <- ifelse(is.infinite(roc$delta), "-Inf (no filter)", as.character(roc$delta))

# Identify key points
roc$point_type <- "other"
roc$point_type[is.infinite(roc$delta)] <- "baseline"
roc$point_type[roc$delta == 0] <- "delta0"

# Order for line drawing
roc <- roc[order(roc$sim_type, roc$delta_num), ]

label_rows <- roc[roc$delta %in% c(-Inf, 0, 0.5, 1, 2) & roc$sim_type %in% c("DTE", "DTU"), ]
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
  ggrepel::geom_text_repel(data = label_rows,
                            aes(label = label),
                            size = 3, show.legend = FALSE,
                            min.segment.length = 0.2,
                            box.padding = 0.3) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1)) +
  scale_color_manual(values = c("DTE" = "#d73027", "DTU" = "#4575b4"),
                     name = NULL) +
  labs(
    title = "lfc_diff_net filter: precision-recall tradeoff",
    subtitle = "Sweep of delta threshold on lfc_diff_net > delta  |  padj < 0.05",
    x = "Recall (restricted, micro)",
    y = "Precision (restricted, micro)",
    caption = "Diamond = delta=0  Circle = baseline (no filter)"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "top")

# Try ggrepel; fall back to geom_text if not installed
tryCatch({
  library(ggrepel)
}, error = function(e) {
  p <<- p +
    geom_text(data = label_rows, aes(label = label),
              size = 3, hjust = -0.1, show.legend = FALSE)
})

ggsave(out_pdf, p, width = 6, height = 5)
message("Written: ", out_pdf)
