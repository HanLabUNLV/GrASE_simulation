library(ggplot2)
library(tidyr)
library(dplyr)
library(gridExtra)

# Define paths
f_prior <- "~/GrASE_simulation/bipartition.test/test_bipartition.internal_glmmTMB_prior.txt"
f_fixed <- "~/GrASE_simulation/bipartition.test/test_bipartition.internal_glmmTMB_fixedEB.txt"
f_noprior  <- "~/GrASE_simulation/bipartition.test/test_bipartition.internal_glmmTMB_no_prior.txt"
f_wilcoxon <- "~/GrASE_simulation/bipartition.test/test_bipartition.internal_wilcoxon.txt"

cat("Reading data...\n")
t_fixed <- read.table(f_fixed, header=TRUE, stringsAsFactors=FALSE)
t_prior <- read.table(f_prior, header=TRUE, stringsAsFactors=FALSE)
t_noprior  <- read.table(f_noprior,  header=TRUE, stringsAsFactors=FALSE)
t_wilcoxon <- read.table(f_wilcoxon, header=TRUE, stringsAsFactors=FALSE)

cat("Processing...\n")
# Select and label
df_fixed <- t_fixed %>% 
  select(gene, event,comparison, p.value, phi, effect_size) %>% 
  mutate(method = "FixedEB")

df_prior <- t_prior %>% 
  select(gene, event,comparison, p.value, phi, effect_size) %>% 
  mutate(method = "Prior")

df_noprior <- t_noprior %>% 
  select(gene, event,comparison, p.value, phi, effect_size) %>% 
  mutate(method = "noprior")

df_wilcoxon <- t_wilcoxon %>% 
  select(gene, event,comparison, p.value, phi, effect_size) %>% 
  mutate(method = "Wilcoxon", phi = as.numeric(phi))

# Combine
all_data <- bind_rows(df_fixed, df_prior, df_noprior, df_wilcoxon)
all_data[is.na(all_data$effect_size),'effect_size'] <- 0

# Clean and transform
# We limit phi to visualization range since noprior goes to 10^69
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

cat("Plotting...\n")

pdf("plots/internalp1.pdf", width=12, height=25)
# P1: Phi Comparison (Prior vs FixedEB)
p1 <- ggplot(wide_df, aes(x=log_phi_Prior, y=log_phi_FixedEB)) +
  geom_point(alpha=0.2, size=0.5) +
  geom_abline(color="red", linetype="dashed") +
  theme_bw() + 
  labs(title="Dispersion (Log10 Phi): Prior vs Fixed EB", 
       subtitle="Red line = 1:1 identity")
print(p1)
dev.off()

pdf("plots/internalp2.pdf", width=12, height=25)
# P2: Phi Comparison (Prior vs noprior)
p2 <- ggplot(wide_df, aes(x=log_phi_Prior, y=log_phi_noprior)) +
  geom_point(alpha=0.2, size=0.5) +
  geom_abline(color="red", linetype="dashed") +
  theme_bw() + 
  labs(title="Dispersion (Log10 Phi): Prior vs noprior",
       subtitle="noprior often estimates infinite (Binomial) phi")
print(p2)
dev.off()

pdf("plots/internalp3.pdf", width=12, height=25)
# P3: P-value Comparison (Prior vs noprior)
p3 <- ggplot(wide_df, aes(x=log_p_Prior, y=log_p_noprior)) +
  geom_point(alpha=0.2, size=0.5) +
  geom_abline(color="red", linetype="dashed") +
  theme_bw() + 
  labs(title="Significance (-Log10 P): Prior vs noprior")
print(p3)
dev.off()

pdf("plots/internalp4.pdf", width=12, height=25)
# P4: P-value Comparison (Prior vs Fixed EB)
p4 <- ggplot(wide_df, aes(x=log_p_Prior, y=log_p_FixedEB)) +
  geom_point(alpha=0.2, size=0.5) +
  geom_abline(color="red", linetype="dashed") +
  theme_bw() + 
  labs(title="Significance (-Log10 P): Prior vs Fixed EB")
print(p4)
dev.off()

pdf("plots/internalp5.pdf", width=12, height=25)
# P5: P-value Comparison (Prior vs Wilcoxon)
p5 <- ggplot(wide_df, aes(x=log_p_Prior, y=log_p_Wilcoxon)) +
  geom_point(alpha=0.2, size=0.5) +
  geom_abline(color="red", linetype="dashed") +
  theme_bw() + 
  labs(title="Significance (-Log10 P): Prior vs Wilcoxon")
print(p5)
dev.off()

pdf("plots/internalp6.pdf", width=12, height=25)
# P6: P-value Comparison (FixedEB vs Wilcoxon)
p6 <- ggplot(wide_df, aes(x=log_p_FixedEB, y=log_p_Wilcoxon)) +
  geom_point(alpha=0.2, size=0.5) +
  geom_abline(color="red", linetype="dashed") +
  theme_bw() + 
  labs(title="Significance (-Log10 P): FixedEB vs Wilcoxon")
print(p6)
dev.off()

pdf("plots/internalp7.pdf", width=12, height=25)
# P7: Effect Size Comparison (Prior vs Fixed EB)
p7 <- ggplot(wide_df, aes(x=effect_size_Prior, y=effect_size_FixedEB)) +
  geom_point(alpha=0.2, size=0.5) +
  geom_abline(color="red", linetype="dashed") +
  theme_bw() + 
  labs(title="Effect Size: Prior vs Fixed EB")
print(p7)
dev.off()

pdf("plots/internalp8.pdf", width=12, height=25)
# P8: Effect Size Comparison (Prior vs noprior)
# Filter outliers beyond +/- 100 for better visualization
p8_data <- wide_df %>% 
  filter(abs(effect_size_noprior) <= 100 & abs(effect_size_Prior) <= 100)

p8 <- ggplot(p8_data, aes(x=effect_size_Prior, y=effect_size_noprior)) +
  geom_point(alpha=0.2, size=0.5) +
  geom_abline(color="red", linetype="dashed") +
  theme_bw() + 
  labs(title="Effect Size: Prior vs noprior (Outliers > +/- 100 removed)")
print(p8)
dev.off()

pdf("plots/internalp9.pdf", width=12, height=25)
# P9: Effect Size Comparison (Prior vs Wilcoxon)
p9 <- ggplot(wide_df, aes(x=effect_size_Prior, y=effect_size_Wilcoxon)) +
  geom_point(alpha=0.2, size=0.5) +
  geom_abline(color="red", linetype="dashed") +
  theme_bw() + 
  labs(title="Effect Size: Prior vs Wilcoxon")
print(p9)
dev.off()

pdf("plots/internalp10.pdf", width=12, height=25)
# P10: Effect Size Comparison (FixedEB vs Wilcoxon)
p10 <- ggplot(wide_df, aes(x=effect_size_FixedEB, y=effect_size_Wilcoxon)) +
  geom_point(alpha=0.2, size=0.5) +
  geom_abline(color="red", linetype="dashed") +
  theme_bw() + 
  labs(title="Effect Size: FixedEB vs Wilcoxon")
print(p10)
dev.off()

pdf("plots/model_comparisons.pdf", width=12, height=25)
grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, nrow=5)
dev.off()

cat("Success. Plots saved to model_comparisons.pdf\n")
