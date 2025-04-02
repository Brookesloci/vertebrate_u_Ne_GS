library(phylopath)
library(ape)
library(dplyr)
library(evobiR)
library(geiger)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(car)
library(dagitty)
library(nlme)
library(ggpubr)
library(nlme)




time_tree <- read.tree(file="tree_calibrated_UCE_morecal.nwk")
my_tree = drop.tip(time_tree, c("Macaca_mulatta","Elephas_maximus", 
                           "Cercocebus_lunulatus", "Muntiacus_reevesi", 
                           "Oryzias_latipes", "Callithrix_jacchus"))
d <- read.csv("primary_dataset.csv")
rownames(d)<-d$Species.name


tip <- c("Arctocephalus_gazella","Betta_splendens","Chrysemys_picta","Fukomys_damarensis","Larus_argentatus","Larus_marinus","Pogona_vitticeps","Rhea_pennata","Rousettus_aegyptiacus","Sphaerodactylus_inigoi","Syngnathus_scovelli","Thamnophis_sirtalis","Tupaia_belangeri" )
time_tree_u <-drop.tip(my_tree,tip)
d <- d %>% drop_na(Modeled.rate.per.generation..m_generation_modeled.)

name.check(time_tree_u,d)
d <-ReorderData(time_tree_u,d)



#d$MA_over_LS <- d$Average.maturation.time/d$Lifespan.in.the.wild..years.

d$Modeled.rate.per.generation..m_generation_modeled. <- log10(d$Modeled.rate.per.generation..m_generation_modeled.)
d$Ne..harmonic.mean..psmc. <- log10(d$Ne..harmonic.mean..psmc.)
d$Generation.time..years. <- log10(d$Generation.time..years.)
d$Average.maturation.time <- log10(d$Average.maturation.time)
d$Lifespan.in.the.wild..years. <- log10(d$Lifespan.in.the.wild..years.)
d$Mass..grams. <- log10(d$Mass..grams.)



###check raw variables
scatter_plot <- function(data, x, y, x_label, y_label) { 
  ggplot(data, aes_string(x = x, y = y)) + 
    geom_point() + 
    theme_minimal() + 
    labs(x = x_label, y = y_label) +
    ggtitle(paste(y_label, "vs", x_label)) + 
    theme(plot.title = element_text(size = 12, face = "bold")) 
}

histogram_plot <- function(data, var, var_label) { 
  ggplot(data, aes_string(x = var)) + 
    geom_histogram(bins = 30, fill = "blue", color = "black", alpha = 0.7) + 
    theme_minimal() + 
    labs(x = var_label) +
    ggtitle(paste("Histogram of", var_label)) + 
    theme(plot.title = element_text(size = 12, face = "bold")) 
}

variables <- c("Modeled.rate.per.generation..m_generation_modeled.", "Generation.time..years.", "Average.maturation.time", "Lifespan.in.the.wild..years.", "Mass..grams.", "Ne..harmonic.mean..psmc.")
labels <- c("u", "GT", "MA", "LS", "BM", "Ne")

plots <- list() 

for (i in 1:length(variables)) { 
  var <- variables[i] 
  label <- labels[i] 
  scatter_plot_ne <- scatter_plot(d, "Ne..harmonic.mean..psmc.", var, "Ne", label) 
  scatter_plot_u <- scatter_plot(d, "Modeled.rate.per.generation..m_generation_modeled.", var, "u", label) 
  histogram_var <- histogram_plot(d, var, label) 
  plots <- c(plots, list(scatter_plot_ne, scatter_plot_u, histogram_var)) 
}

grid.arrange(grobs = plots, ncol = 3)

options(upper.bound = 1, lower.bound = 0)


## DAG business

DAG1 <- dagitty("dag{
                OPG -> Ne
                OPG -> u
                OPG -> GT
                GT -> Ne
                GT -> u
                MS -> GT
                MS -> LS
                BM -> LS
                BM -> Ne
                Ne -> u
                MS -> Ne
                u -> LS
                BM -> u
                MA -> LS
                MA -> GT
                BM -> MA
                BM -> MS
                OPG -> MA}")

plot(DAG1)

impliedConditionalIndependencies(DAG1)

adjustmentSets(DAG1,
               exposure = "Ne",
               outcome = "u",
               effect = "direct")

# Define models
paths(DAG1, "Ne", "u")

path1_d = data.frame(d$Ne..harmonic.mean..psmc.,d$Modeled.rate.per.generation..m_generation_modeled., d$Generation.time..years.,
                        as.factor(d$Mating.strategy..binary.), 
                        d$Mass..grams.,
                        d$Lifespan.in.the.wild..years., 
                        as.factor(d$Number.of.offspring.per.generation),
                        d$Average.maturation.time)
rownames(path1_d) <- d$Species.name
colnames(path1_d) <- c("Ne","u", "GT", "MS", "BM", "LS", "OPG", "MA")


m1 <- define_model_set( 
  null = c(), 
  direct = c(u ~ Ne), 
  indirect = c(u ~GT + BM + OPG, Ne ~ GT + BM + OPG + MS, GT ~ MS + MA +OPG, BM ~ MS + MA, OPG ~ MA),
  both = c(u~Ne,u ~GT + BM + OPG, Ne ~ GT + BM + OPG + MS, GT ~ MS + MA, BM ~ MS + MA, OPG ~ MA),
  .common = c(LS ~ u + BM + MA + MS)
)

p1 <- phylo_path(m1,path1_d,time_tree_u)
plot_model_set(m1)
p1
s1 <- summary(p1)
s1

plot(s1)
best_model1 <- best(p1, boot = 1000)
coef_plot(best_model1)

coef <- best_model1$coef
lower <- best_model1$lower
upper <- best_model1$upper

# Get row and column names (excluding the diagonal)
rows <- rownames(coef)
cols <- colnames(coef)

# Initialize a data frame to store results
results <- data.frame(From = character(), To = character(), Estimate = numeric(), Lower = numeric(), Upper = numeric(), Significant = logical())

# Iterate through the coefficient matrix (excluding the diagonal)
for (i in 1:nrow(coef)) {
  for (j in 1:ncol(coef)) {
    if (i != j) { # Avoid diagonal elements (self-effects)
      from <- rows[i]
      to <- cols[j]
      estimate <- coef[i, j]
      lower_bound <- lower[i, j]
      upper_bound <- upper[i, j]
      significant <- (lower_bound > 0 || upper_bound < 0) # Check if CI overlaps zero
      
      results <- rbind(results, data.frame(From = from, To = to, Estimate = estimate, Lower = lower_bound, Upper = upper_bound, Significant = significant))
    }
  }
}

print(results)
write.csv(results, file = "best_model_coefficients.csv", row.names = FALSE)



both_model <- choice(p1,"both", boot=1000)

coef_plot(both_model)

coef1 <- both_model$coef
lower1 <- both_model$lower
upper1 <- both_model$upper

# Get row and column names (excluding the diagonal)
rows1 <- rownames(coef1)
cols1 <- colnames(coef1)

# Initialize a data frame to store results
results1 <- data.frame(From = character(), To = character(), Estimate = numeric(), Lower = numeric(), Upper = numeric(), Significant = logical())

# Iterate through the coefficient matrix (excluding the diagonal)
for (i in 1:nrow(coef1)) {
  for (j in 1:ncol(coef1)) {
    if (i != j) { # Avoid diagonal elements (self-effects)
      from <- rows1[i]
      to <- cols1[j]
      estimate <- coef1[i, j]
      lower_bound <- lower1[i, j]
      upper_bound <- upper1[i, j]
      significant <- (lower_bound > 0 || upper_bound < 0) # Check if CI overlaps zero
      
      results1 <- rbind(results, data.frame(From = from, To = to, Estimate = estimate, Lower = lower_bound, Upper = upper_bound, Significant = significant))
    }
  }
}

print(results1)
write.csv(results, file = "S2_best_model_coefficients.csv", row.names = FALSE)

modify_plot <- function(p, node_text_size = 8, edge_text_size = 6) { # Added edge_text_size argument
  for (i in seq_along(p$layers)) {
    if (inherits(p$layers[[i]]$geom, "GeomText")) {
      p$layers[[i]]$aes_params$size <- node_text_size # Node text size
    } else if (inherits(p$layers[[i]]$geom, "GeomLabel")) { # Check for edge labels
      p$layers[[i]]$aes_params$size <- edge_text_size # Edge text size
    }
  }
  return(p)
}



plot_best <- modify_plot(plot_best, node_text_size = 6, edge_text_size = 5)  # Reduced edge label size
plot_both <- modify_plot(plot_both, node_text_size = 6, edge_text_size = 5)  # Reduced edge label size


plot_best <- plot(best_model1) + 
  ggtitle("(a) Best-Supported Indirect Effects Model") +
  theme(plot.title = element_text(size = 12, hjust = 0.5))

plot_best <- modify_plot(plot_best, node_text_size = 6) # Call the function to change node text size

plot_both <- plot(both_model) + 
  ggtitle("(b) Combined Model") +
  theme(plot.title = element_text(size = 12, hjust = 0.5))

plot_both <- modify_plot(plot_both, node_text_size = 6) # Call the function to change node text size


combined_plot <- ggarrange(
  plot_best, 
  plot_both, 
  ncol = 1,
  # Reduce the space between plots
  heights = c(5, 5)
)

combined_plot

ggsave("combined_path_diagram.pdf", combined_plot, width = 10, height = 10, units = "in", dpi = 300) 







coef_plot(best_model1, error_bar = "se", order_by = "strength", to="u") + ggplot2::coord_flip()






coef_plot(both_model, error_bar = "se", order_by = "strength", to ="u") + ggplot2::coord_flip()


coef_table_indirect <- best_model1$coef
coef_table_both <- both_model$coef
indirect_SE <- best_model1$se
both_SE <- both_model$se


pgls1_d = data.frame(d$Ne..harmonic.mean..psmc.,d$Modeled.rate.per.generation..m_generation_modeled., d$Generation.time..years.,
                     d$Mass..grams.,
                     as.factor(d$Number.of.offspring.per.generation))
rownames(pgls1_d) <- d$Species.name
colnames(pgls1_d) <- c("Ne","u", "GT", "BM", "OPG")


summary(best_model1)


na_counts_base <- colSums(is.na(d))
print(na_counts_base)








