library(ggplot2)
library(tidyr)
library(dplyr)
library(gridExtra)

# Define paths
f_EBmap <- "~/GrASE_simulation/bipartition.test/test_bipartition.internal_betabinom_EBmap.txt"
f_EBapprox <- "~/GrASE_simulation/bipartition.test/test_bipartition.internal_betabinom_EBapprox.txt"
f_MLE  <- "~/GrASE_simulation/bipartition.test/test_bipartition.internal_betabinom_MLE.txt"
f_wilcoxon <- "~/GrASE_simulation/bipartition.test/test_bipartition.internal_wilcoxon.txt"

cat("Reading data...\n")
t_EBapprox <- read.table(f_EBapprox, header=TRUE, stringsAsFactors=FALSE)
t_EBmap <- read.table(f_EBmap, header=TRUE, stringsAsFactors=FALSE)
t_MLE  <- read.table(f_MLE,  header=TRUE, stringsAsFactors=FALSE)
t_wilcoxon <- read.table(f_wilcoxon, header=TRUE, stringsAsFactors=FALSE)

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

cat("Plotting...\n")

pdf("plots/internalp1.pdf", width=12, height=25)
# P1: Phi Comparison (EBmap vs EBapprox)
p1 <- ggplot(wide_df, aes(x=log_phi_EBmap, y=log_phi_EBapprox)) +
  geom_point(alpha=0.2, size=0.5) +
  geom_abline(color="red", linetype="dashed") +
  theme_bw() + 
  labs(title="Dispersion (Log10 Phi): EBmap vs EBapprox", 
       subtitle="Red line = 1:1 identity")
print(p1)
dev.off()

pdf("plots/internalp2.pdf", width=12, height=25)
# P2: Phi Comparison (EBmap vs MLE)
p2 <- ggplot(wide_df, aes(x=log_phi_EBmap, y=log_phi_MLE)) +
  geom_point(alpha=0.2, size=0.5) +
  geom_abline(color="red", linetype="dashed") +
  theme_bw() + 
  labs(title="Dispersion (Log10 Phi): EBmap vs MLE",
       subtitle="MLE often estimates infinite (Binomial) phi")
print(p2)
dev.off()

pdf("plots/internalp3.pdf", width=12, height=25)
# P3: P-value Comparison (EBmap vs MLE)
p3 <- ggplot(wide_df, aes(x=log_p_EBmap, y=log_p_MLE)) +
  geom_point(alpha=0.2, size=0.5) +
  geom_abline(color="red", linetype="dashed") +
  theme_bw() + 
  labs(title="Significance (-Log10 P): EBmap vs MLE")
print(p3)
dev.off()

pdf("plots/internalp4.pdf", width=12, height=25)
# P4: P-value Comparison (EBmap vs EBapprox)
p4 <- ggplot(wide_df, aes(x=log_p_EBmap, y=log_p_EBapprox)) +
  geom_point(alpha=0.2, size=0.5) +
  geom_abline(color="red", linetype="dashed") +
  theme_bw() + 
  labs(title="Significance (-Log10 P): EBmap vs EBapprox")
print(p4)
dev.off()

pdf("plots/internalp5.pdf", width=12, height=25)
# P5: P-value Comparison (EBmap vs Wilcoxon)
p5 <- ggplot(wide_df, aes(x=log_p_EBmap, y=log_p_Wilcoxon)) +
  geom_point(alpha=0.2, size=0.5) +
  geom_abline(color="red", linetype="dashed") +
  theme_bw() + 
  labs(title="Significance (-Log10 P): EBmap vs Wilcoxon")
print(p5)
dev.off()

pdf("plots/internalp6.pdf", width=12, height=25)
# P6: P-value Comparison (EBapprox vs Wilcoxon)
p6 <- ggplot(wide_df, aes(x=log_p_EBapprox, y=log_p_Wilcoxon)) +
  geom_point(alpha=0.2, size=0.5) +
  geom_abline(color="red", linetype="dashed") +
  theme_bw() + 
  labs(title="Significance (-Log10 P): EBapprox vs Wilcoxon")
print(p6)
dev.off()

pdf("plots/internalp7.pdf", width=12, height=25)
# P7: Effect Size Comparison (EBmap vs EBapprox)
p7 <- ggplot(wide_df, aes(x=effect_size_EBmap, y=effect_size_EBapprox)) +
  geom_point(alpha=0.2, size=0.5) +
  geom_abline(color="red", linetype="dashed") +
  theme_bw() + 
  labs(title="Effect Size: EBmap vs EBapprox")
print(p7)
dev.off()

pdf("plots/internalp8.pdf", width=12, height=25)
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

pdf("plots/internalp9.pdf", width=12, height=25)
# P9: Effect Size Comparison (EBmap vs Wilcoxon)
p9 <- ggplot(wide_df, aes(x=effect_size_EBmap, y=effect_size_Wilcoxon)) +
  geom_point(alpha=0.2, size=0.5) +
  geom_abline(color="red", linetype="dashed") +
  theme_bw() + 
  labs(title="Effect Size: EBmap vs Wilcoxon")
print(p9)
dev.off()

pdf("plots/internalp10.pdf", width=12, height=25)
# P10: Effect Size Comparison (EBapprox vs Wilcoxon)
p10 <- ggplot(wide_df, aes(x=effect_size_EBapprox, y=effect_size_Wilcoxon)) +
  geom_point(alpha=0.2, size=0.5) +
  geom_abline(color="red", linetype="dashed") +
  theme_bw() + 
  labs(title="Effect Size: EBapprox vs Wilcoxon")
print(p10)
dev.off()

pdf("plots/model_comparisons.pdf", width=12, height=25)
grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, nrow=5)
dev.off()

cat("Success. Plots saved to model_comparisons.pdf\n")
