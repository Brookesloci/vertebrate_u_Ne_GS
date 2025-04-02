library(ape)
library(brms)
library(dplyr)
library(tidyr)
library(geiger)
library(evobiR)
library(ggplot2)
library(ggpubr)



time_tree <- read.tree(file="tree_calibrated_UCE_morecal.nwk")
my_tree = drop.tip(time_tree, c("Macaca_mulatta","Elephas_maximus", 
                                "Cercocebus_lunulatus", "Muntiacus_reevesi", 
                                "Oryzias_latipes", "Callithrix_jacchus"))
d <- read.csv("my_Ne_u_data_Nc.csv")
rownames(d)<-d$Species.name


tip <- c("Arctocephalus_gazella","Betta_splendens","Chrysemys_picta","Fukomys_damarensis","Larus_argentatus","Larus_marinus","Pogona_vitticeps","Rhea_pennata","Rousettus_aegyptiacus","Sphaerodactylus_inigoi","Syngnathus_scovelli","Thamnophis_sirtalis","Tupaia_belangeri" )
time_tree_u <-drop.tip(my_tree,tip)
d <- d %>% drop_na(Modeled.rate.per.generation..m_generation_modeled.)

name.check(time_tree_u,d)
d <-ReorderData(time_tree_u,d)
rownames(d) <- NULL

## construct covariance matrix

A <- ape::vcv.phylo(time_tree_u)

m1_d <- data.frame(d$Modeled.rate.per.generation..m_generation_modeled.,d$Ne..harmonic.mean..psmc., d$Species.name)
colnames(m1_d) <- c("u", "Ne","phylo")

m1_d$Ne <- scale(log10(m1_d$Ne))
m1_d$u <- scale(log10(m1_d$u))
m1 <- brm( u ~ Ne + (1|gr(phylo, cov=A)), data = m1_d, family = gaussian(), data2 = list(A = A), prior = c( prior(normal(0, 1), "b"), prior(normal(0, 1), "Intercept"), prior(student_t(3, 0, 2.5), "sd"), prior(student_t(3, 0, 2.5), "sigma") ), control = list(adapt_delta = 0.999), iter = 10000 ,save_pars = save_pars(all = TRUE))

summary(m1)
plot(m1)

m2_d <- data.frame(d$Modeled.rate.per.generation..m_generation_modeled.,d$Ne..harmonic.mean..psmc., d$Species.name, d$Generation.time..years.)
colnames(m2_d) <- c("u", "Ne","phylo","GT")

m2_d$Ne <- scale(log10(m2_d$Ne))
m2_d$u <- scale(log10(m2_d$u))
m2_d$GT <- scale(log10(m2_d$GT))

m2 <- brm( u ~ Ne + GT + (1|gr(phylo, cov=A)), data = m2_d, family = gaussian(), data2 = list(A = A), prior = c( prior(normal(0, 1), "b"), prior(normal(0, 1), "Intercept"), prior(student_t(3, 0, 2.5), "sd"), prior(student_t(3, 0, 2.5), "sigma") ), control = list(adapt_delta = 0.999), iter = 10000 ,save_pars = save_pars(all = TRUE))
summary(m2)
plot(m2)

options(future.globals.maxSize = 2000 * 1024^2) 

loo_m1 <- loo(m1, moment_match = TRUE, reloo = TRUE)
loo_m2 <- loo(m2, moment_match = TRUE, reloo = TRUE)
loo_compare(loo_m2, loo_m1)

waic_m1 <- waic(m1)
waic_m2 <- waic(m2)
waic_1 <- waic_m1$estimates["waic", "Estimate"]
waic_1
waic_2 <- waic_m2$estimates["waic", "Estimate"]
waic_2
delta_waic <- waic_1 - waic_2
delta_waic

loo_m1_k <- loo(m1, moment_match = TRUE, reloo = TRUE)$diagnostics$pareto_k
loo_m2_k <- loo(m2, moment_match = TRUE, reloo = TRUE)$diagnostics$pareto_k


high_k_m1 <- data.frame(Species = m1_d$phylo, Pareto_k = loo_m1_k) %>%
  filter(Pareto_k > 0.7)


high_k_m2 <- data.frame(Species = m2_d$phylo, Pareto_k = loo_m2_k) %>%
  filter(Pareto_k > 0.7)


m2_d$Taxonomic.group <- d$Taxonomic.group
m1_d$Taxonomic.group <- d$Taxonomic.group
my_colors <- c(
  "Bird" = "#EBCC2A",
  "Reptile" = "#3B9AB2",
  "Mammal" = "#35274A",
  "Fish" = "#F2300F"
)




conditional_effects_ne_m1_data <- conditional_effects(m1, effects = "Ne", plot = FALSE)[[1]]

# effects plot for Ne from m1
marginal_ne_m1_plot <- ggplot(conditional_effects_ne_m1_data, aes(x = Ne, y = estimate__)) +
  geom_line(color = "red") +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.2, fill = "lightcoral") +
  geom_point(data = m1_d, aes(x = Ne, y = u, color = Taxonomic.group), size = 2, alpha = 0.6) +
  scale_color_manual(values = my_colors, name = "Taxonomic Group") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "none"
  ) +
  xlab("Standardized log10(Ne)") +
  ylab("Standardized log10(u)") +
  ggtitle("Effect of Ne on u")

marginal_ne_m1_plot


# Extract marginal effects data for Ne
conditional_effects_ne_data <- conditional_effects(m2, effects = "Ne", plot = FALSE)[[1]]

# Marginal effects plot for Ne
marginal_ne_plot <- ggplot(conditional_effects_ne_data, aes(x = Ne, y = estimate__)) +
  geom_line(color = "blue") +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.2, fill = "lightblue") +
  geom_point(data = m2_d, aes(x = Ne, y = u, color = Taxonomic.group), size = 2, alpha = 0.6) +
  scale_color_manual(values = my_colors, name = "Taxonomic Group") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "none"
  ) +
  xlab("Standardized log10(Ne)") +
  ylab("Standardized log10(u)") +
  ggtitle("Marginal Effect of Ne on u")
marginal_ne_plot

# Marginal effects plot for GT
conditional_effects_gt_data <- conditional_effects(m2, effects = "GT", plot = FALSE)[[1]]

marginal_gt_plot <- ggplot(conditional_effects_gt_data, aes(x = GT, y = estimate__)) +
  geom_line(color = "green") +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.2, fill = "lightgreen") +
  geom_point(data = m2_d, aes(x = GT, y = u, color = Taxonomic.group), size = 2, alpha = 0.6) +
  scale_color_manual(values = my_colors, name = "Taxonomic Group") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "none"
  ) +
  xlab("Standardized log10(GT)") +
  ylab("Standardized log10(u)") +
  ggtitle("Marginal Effect of GT on u")


marginal_gt_plot



combined_figure <- ggarrange(marginal_ne_m1_plot, marginal_ne_plot, marginal_gt_plot,
                             ncol = 2, nrow = 2, labels = c("A", "B", "C"),
                             common.legend = FALSE)
print(combined_figure)
