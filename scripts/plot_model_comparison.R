library(ggplot2)
library(tidyr)
library(dplyr)
library(gridExtra)

# Output directory for plots (absolute path, independent of working directory)
OUTDIR <- "/mnt/data1/home/mirahan/GrASE_simulation/scripts/plots/model_comparison"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# Define paths
f_EBmap <- "~/GrASE_simulation/bipartition.test/test_bipartition.internal_betabinom_EBmap.txt"
f_EBapprox <- "~/GrASE_simulation/bipartition.test/test_bipartition.internal_betabinom_EBapprox.txt"
f_MLE  <- "~/GrASE_simulation/bipartition.test/test_bipartition.internal_betabinom_MLE.txt"
f_wilcoxon <- "~/GrASE_simulation/bipartition.test/test_bipartition.internal_wilcoxon.txt"
# Raw per-event phi table carries the LRT identifiability flag (BB significantly
# overdispersed vs binomial). Used to color the dispersion scatter.
f_phi <- "~/GrASE_simulation/bipartition.test/phi.glmmtmb.internal.txt"

cat("Reading data...\n")
t_EBapprox <- read.table(f_EBapprox, header=TRUE, stringsAsFactors=FALSE)
t_EBmap <- read.table(f_EBmap, header=TRUE, stringsAsFactors=FALSE)
t_MLE  <- read.table(f_MLE,  header=TRUE, stringsAsFactors=FALSE)
t_wilcoxon <- read.table(f_wilcoxon, header=TRUE, stringsAsFactors=FALSE)

# Identifiability flag per gene/event/comparison (may be absent in older runs).
ident_df <- tryCatch({
  tp <- read.table(f_phi, header=TRUE, stringsAsFactors=FALSE)
  if ("identifiable" %in% names(tp)) {
    tp %>% select(gene, event, comparison, identifiable) %>%
      mutate(identifiable = as.logical(identifiable))
  } else NULL
}, error = function(e) NULL)

cat("Processing...\n")
# Select and label
df_EBapprox <- t_EBapprox %>% 
  select(gene, event,comparison, p.value, phi, effect_size) %>% 
  mutate(method = "EBapprox")

df_EBmap <- t_EBmap %>% 
  select(gene, event,comparison, p.value, phi, effect_size) %>% 
  mutate(method = "EBmap")

df_MLE <- t_MLE %>% 
  select(gene, event,comparison, p.value, phi, effect_size) %>% 
  mutate(method = "MLE")

df_wilcoxon <- t_wilcoxon %>% 
  select(gene, event,comparison, p.value, phi, effect_size) %>% 
  mutate(method = "Wilcoxon", phi = as.numeric(phi))

# Combine
all_data <- bind_rows(df_EBapprox, df_EBmap, df_MLE, df_wilcoxon)
all_data[is.na(all_data$effect_size),'effect_size'] <- 0

# Clean and transform
# We limit phi to visualization range since MLE goes to 10^69
MAX_PHI_PLOT <- 20000 

all_data <- all_data %>%
  mutate(
    p.value = suppressWarnings(as.numeric(p.value)),
    # Handle extremely small p-values (avoid log(0))
    p_clamped = pmax(p.value, 1e-300),
    log_p = -log10(p_clamped),
    
    phi = suppressWarnings(as.numeric(phi)),
    # Clamp phi for plotting visualization
    phi_plot = ifelse(is.na(phi) | phi > MAX_PHI_PLOT, MAX_PHI_PLOT, phi),
    log_phi = log10(phi_plot + 0.1),

    effect_size = suppressWarnings(as.numeric(effect_size))
  )

# Widen for direct comparison
wide_df <- all_data %>%
  select(gene, event,comparison, method, log_phi, log_p, effect_size) %>%
  pivot_wider(
    names_from = method,
    values_from = c(log_phi, log_p, effect_size),
    names_sep = "_"
  )

# Attach the identifiability flag (per gene/event/comparison) for coloring.
if (!is.null(ident_df)) {
  wide_df <- wide_df %>% left_join(ident_df, by = c("gene", "event", "comparison"))
  wide_df$ident_lab <- ifelse(is.na(wide_df$identifiable), "no phi estimate",
                       ifelse(wide_df$identifiable, "identifiable (BB > binomial)",
                                                    "non-identifiable (~binomial)"))
  wide_df$ident_lab <- factor(wide_df$ident_lab,
    levels = c("identifiable (BB > binomial)", "non-identifiable (~binomial)", "no phi estimate"))
}

cat("Plotting...\n")

pdf(file.path(OUTDIR, "internalp1.pdf"), width=12, height=12)
# P1: Phi Comparison (EBmap vs EBapprox), colored by identifiability.
# The off-diagonal "square" is the non-identifiable events: with no dispersion
# signal, EBmap (MAP) and EBapprox (variance-weighted shrinkage) fill in
# different moderated values, so they scatter instead of lying on the 1:1 line.
if (!is.null(ident_df)) {
  p1 <- ggplot(wide_df, aes(x=log_phi_EBmap, y=log_phi_EBapprox, color=ident_lab)) +
    geom_point(alpha=0.25, size=0.5) +
    geom_abline(color="black", linetype="dashed") +
    scale_color_manual(values=c("identifiable (BB > binomial)"   = "#2166ac",
                                "non-identifiable (~binomial)"    = "#d6604d",
                                "no phi estimate"                 = "grey60"),
                       name="dispersion", drop=FALSE) +
    guides(color = guide_legend(override.aes = list(alpha=1, size=2))) +
    theme_bw() +
    theme(legend.position="top") +
    labs(title="Dispersion (Log10 Phi): EBmap vs EBapprox",
         subtitle="Dashed line = 1:1 identity. Off-diagonal square = non-identifiable events (no dispersion signal).")
} else {
  p1 <- ggplot(wide_df, aes(x=log_phi_EBmap, y=log_phi_EBapprox)) +
    geom_point(alpha=0.2, size=0.5) +
    geom_abline(color="red", linetype="dashed") +
    theme_bw() +
    labs(title="Dispersion (Log10 Phi): EBmap vs EBapprox",
         subtitle="Red line = 1:1 identity")
}
print(p1)
dev.off()

pdf(file.path(OUTDIR, "internalp2.pdf"), width=12, height=12)
# P2: Phi Comparison (EBmap vs MLE), colored by identifiability.
if (!is.null(ident_df)) {
  p2 <- ggplot(wide_df, aes(x=log_phi_EBmap, y=log_phi_MLE, color=ident_lab)) +
    geom_point(alpha=0.25, size=0.5) +
    geom_abline(color="black", linetype="dashed") +
    scale_color_manual(values=c("identifiable (BB > binomial)"   = "#2166ac",
                                "non-identifiable (~binomial)"    = "#d6604d",
                                "no phi estimate"                 = "grey60"),
                       name="dispersion", drop=FALSE) +
    guides(color = guide_legend(override.aes = list(alpha=1, size=2))) +
    theme_bw() +
    theme(legend.position="top") +
    labs(title="Dispersion (Log10 Phi): EBmap vs MLE",
         subtitle="Dashed line = 1:1 identity. MLE leaves non-identifiable phi at the binomial-limit boundary.")
} else {
  p2 <- ggplot(wide_df, aes(x=log_phi_EBmap, y=log_phi_MLE)) +
    geom_point(alpha=0.2, size=0.5) +
    geom_abline(color="red", linetype="dashed") +
    theme_bw() +
    labs(title="Dispersion (Log10 Phi): EBmap vs MLE",
         subtitle="MLE often estimates infinite (Binomial) phi")
}
print(p2)
dev.off()

pdf(file.path(OUTDIR, "internalp3.pdf"), width=12, height=12)
# P3: P-value Comparison (EBmap vs MLE)
p3 <- ggplot(wide_df, aes(x=log_p_EBmap, y=log_p_MLE)) +
  geom_point(alpha=0.2, size=0.5) +
  geom_abline(color="red", linetype="dashed") +
  theme_bw() + 
  labs(title="Significance (-Log10 P): EBmap vs MLE")
print(p3)
dev.off()

pdf(file.path(OUTDIR, "internalp4.pdf"), width=12, height=12)
# P4: P-value Comparison (EBmap vs EBapprox)
p4 <- ggplot(wide_df, aes(x=log_p_EBmap, y=log_p_EBapprox)) +
  geom_point(alpha=0.2, size=0.5) +
  geom_abline(color="red", linetype="dashed") +
  theme_bw() + 
  labs(title="Significance (-Log10 P): EBmap vs EBapprox")
print(p4)
dev.off()

pdf(file.path(OUTDIR, "internalp5.pdf"), width=12, height=12)
# P5: P-value Comparison (EBmap vs Wilcoxon)
p5 <- ggplot(wide_df, aes(x=log_p_EBmap, y=log_p_Wilcoxon)) +
  geom_point(alpha=0.2, size=0.5) +
  geom_abline(color="red", linetype="dashed") +
  theme_bw() + 
  labs(title="Significance (-Log10 P): EBmap vs Wilcoxon")
print(p5)
dev.off()

pdf(file.path(OUTDIR, "internalp6.pdf"), width=12, height=12)
# P6: P-value Comparison (EBapprox vs Wilcoxon)
p6 <- ggplot(wide_df, aes(x=log_p_EBapprox, y=log_p_Wilcoxon)) +
  geom_point(alpha=0.2, size=0.5) +
  geom_abline(color="red", linetype="dashed") +
  theme_bw() + 
  labs(title="Significance (-Log10 P): EBapprox vs Wilcoxon")
print(p6)
dev.off()

pdf(file.path(OUTDIR, "internalp7.pdf"), width=12, height=12)
# P7: Effect Size Comparison (EBmap vs EBapprox)
p7 <- ggplot(wide_df, aes(x=effect_size_EBmap, y=effect_size_EBapprox)) +
  geom_point(alpha=0.2, size=0.5) +
  geom_abline(color="red", linetype="dashed") +
  theme_bw() + 
  labs(title="Effect Size: EBmap vs EBapprox")
print(p7)
dev.off()

pdf(file.path(OUTDIR, "internalp8.pdf"), width=12, height=12)
# P8: Effect Size Comparison (EBmap vs MLE)
# Filter outliers beyond +/- 100 for better visualization
p8_data <- wide_df %>% 
  filter(abs(effect_size_MLE) <= 100 & abs(effect_size_EBmap) <= 100)

p8 <- ggplot(p8_data, aes(x=effect_size_EBmap, y=effect_size_MLE)) +
  geom_point(alpha=0.2, size=0.5) +
  geom_abline(color="red", linetype="dashed") +
  theme_bw() + 
  labs(title="Effect Size: EBmap vs MLE (Outliers > +/- 100 removed)")
print(p8)
dev.off()

pdf(file.path(OUTDIR, "internalp9.pdf"), width=12, height=12)
# P9: Effect Size Comparison (EBmap vs Wilcoxon)
p9 <- ggplot(wide_df, aes(x=effect_size_EBmap, y=effect_size_Wilcoxon)) +
  geom_point(alpha=0.2, size=0.5) +
  geom_abline(color="red", linetype="dashed") +
  theme_bw() + 
  labs(title="Effect Size: EBmap vs Wilcoxon")
print(p9)
dev.off()

pdf(file.path(OUTDIR, "internalp10.pdf"), width=12, height=12)
# P10: Effect Size Comparison (EBapprox vs Wilcoxon)
p10 <- ggplot(wide_df, aes(x=effect_size_EBapprox, y=effect_size_Wilcoxon)) +
  geom_point(alpha=0.2, size=0.5) +
  geom_abline(color="red", linetype="dashed") +
  theme_bw() + 
  labs(title="Effect Size: EBapprox vs Wilcoxon")
print(p10)
dev.off()

pdf(file.path(OUTDIR, "model_comparisons.pdf"), width=12, height=25)
grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, nrow=5)
dev.off()

cat("Success. Plots saved to model_comparisons.pdf\n")
