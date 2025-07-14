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
library(tibble)
library(caper)



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

d$Nc_u <- d$Nc * d$Modeled.rate.per.generation..m_generation_modeled.


d$expected_pi <- 4* d$Nc_u

d$UCE_mu <- (d$UCE.substitution.rate/1000000)* d$Generation.time..years.

d$Ne_pi <- d$diversity/(4*d$UCE_mu)

d$Ne_x_u <- d$Ne..harmonic.mean..psmc. * d$Modeled.rate.per.generation..m_generation_modeled.
d$piNe_x_u <- d$Ne_pi * d$Modeled.rate.per.generation..m_generation_modeled.

d$piNe_Nc <- d$Ne_pi/d$Nc
d$Ne_Nc <- d$Ne..harmonic.mean..psmc./d$Nc

d$piNe_Nc <- log10(d$piNe_Nc)
d$Ne_Nc <- log10(d$Ne_Nc)


d$piNe_x_u <- log10(d$piNe_x_u)
d$Ne_x_u <- log(d$Ne_x_u)
d$genome_size <- log10(d$genome_size)
d$Modeled.rate.per.generation..m_generation_modeled. <- log10(d$Modeled.rate.per.generation..m_generation_modeled.)
d$Ne..harmonic.mean..psmc. <- log10(d$Ne..harmonic.mean..psmc.)
d$Generation.time..years. <- log10(d$Generation.time..years.)
d$Average.maturation.time <- log10(d$Average.maturation.time)
d$Lifespan.in.the.wild..years. <- log10(d$Lifespan.in.the.wild..years.)
d$Mass..grams. <- log10(d$Mass..grams.)
d$Ne_pi <- log10(d$Ne_pi)
d$Mating.strategy..binary. <- as.factor(d$Mating.strategy..binary.)
d$Number.of.offspring.per.generation <- as.factor(d$Number.of.offspring.per.generation)
d$Nc_u <- log10(d$Nc_u)
d$Nc <- log10(d$Nc)
d$expected_pi <- log10(d$expected_pi)
d$diversity <- log10(d$diversity)


my_colors <- c(
  "Bird" = "#EBCC2A",
  "Reptile" = "#3B9AB2",
  "Mammal" = "#35274A",
  "Fish" = "#F2300F"
)

columns_to_plot <- c("Mass..grams.", "Average.maturation.time", "Lifespan.in.the.wild..years.", 
                     "Modeled.rate.per.generation..m_generation_modeled.", 
                     "Ne..harmonic.mean..psmc.","diversity") 

plot_list <- list()

# Loop through each column and create a scatter plot against Generation.time..years.
for (col in columns_to_plot) {
  p <- ggplot(d, aes_string(x = "Generation.time..years.", y = col, color = "Taxonomic.group")) +
    geom_point() +
    labs(x = "Generation time (years)", y = col) +
    scale_color_manual(values = my_colors) +  # Apply the color palette here
    theme_minimal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                            panel.background = element_blank(), axis.line = element_line(colour = "black"))
  plot_list <- c(plot_list, list(p))
}

# Combine all plots into one figure with a common legend
combined_plot <- ggarrange(plotlist = plot_list, ncol = 2, nrow = ceiling(length(plot_list)/2),
                           common.legend = TRUE, legend = "right", labels = c("A", "B", "C", "D", "E", "F"))

combined_plot



plot(d$Average.maturation.time ~d$Mating.strategy..binary.)
plot(d$Modeled.rate.per.generation..m_generation_modeled.~d$Mating.strategy..binary.)
plot(d$Mass..grams.~d$Number.of.offspring.per.generation)



path1_d = data.frame(d$Ne..harmonic.mean..psmc.,d$Modeled.rate.per.generation..m_generation_modeled., d$Generation.time..years.,
                     d$Mating.strategy..binary., 
                     d$Mass..grams.,
                     d$Lifespan.in.the.wild..years., 
                     d$Average.maturation.time,
                     d$Number.of.offspring.per.generation)
rownames(path1_d) <- d$Species.name
colnames(path1_d) <- c("Ne","u", "GT", "MS", "BM", "LS", "MA","OPG")

m1<- define_model_set(
  A = c(u ~ Ne + MS + BM + LS + OPG + GT +MA)
)

plot_model_set(m1)

p1 <- phylo_path(m1,path1_d,time_tree_u, lower.bound = 0.000, upper.bound = 1.00)

s1 <- summary(p1)
s1

dsep_p1 <- as.data.frame (p1$d_sep)
dsep_p1
dsep_p1$A.model <- NULL
dsep_p1$significant <- dsep_p1$A.p < 0.05
dsep_p1[] <- lapply(dsep_p1, function(x) if (is.list(x)) as.character(x) else x)
#write.csv(dsep_p1, "modelA_dsep_results.csv", row.names = FALSE)

best1 <- best(p1)
plot(best1)


m2 <- define_model_set(
  B = c(u ~ Ne + MS + BM + LS + OPG + GT +MA, Ne ~ GT, GT ~ MA, GT ~ LS)
)

plot_model_set(m2)

p2 <- phylo_path(m2,path1_d,time_tree_u, lower.bound = 0.000, upper.bound = 1.00)

s2 <- summary(p2)
s2

dsep_p2 <- as.data.frame (p2$d_sep)
dsep_p2$B.model <- NULL
dsep_p2$significant <- dsep_p2$B.p < 0.05
dsep_p2[] <- lapply(dsep_p2, function(x) if (is.list(x)) as.character(x) else x)
dsep_p2

#write.csv(dsep_p2, "modelB_dsep_results.csv", row.names = FALSE)


best2 <- best(p2)
plot(best2)


m3 <- define_model_set(
  C = c(
    u ~ Ne + MS + BM + LS + GT + MA,
    Ne ~ GT,
    GT ~ MA + LS,
    LS ~ BM,
    MA ~ BM,
    OPG ~ BM + MA + LS + GT
  )
)

plot_model_set(m3)

p3<- phylo_path(m3,path1_d,time_tree_u, lower.bound = 0.000, upper.bound = 1.00)

s3 <- summary(p3)
s3

dsep_p3 <- as.data.frame (p3$d_sep)
dsep_p3$C.model <- NULL
dsep_p3$significant <- dsep_p3$C.p < 0.05
dsep_p3[] <- lapply(dsep_p3, function(x) if (is.list(x)) as.character(x) else x)
dsep_p3
#write.csv(dsep_p3, "modelC_dsep_results.csv", row.names = FALSE)

best3 <- best(p3)
plot(best3)


m4 <- define_model_set(
  D = c(
    u ~ Ne + MS + BM + LS + GT + MA,
    Ne ~ GT,
    GT ~ MA + LS,
    LS ~ BM,
    MA ~ BM,
    OPG ~ BM + MA + LS + GT,
    LS ~ MA,
    MS ~ BM,
    MS ~ MA
  )
)


plot_model_set(m4)

p4<- phylo_path(m4,path1_d,time_tree_u, lower.bound = 0.000, upper.bound = 1.00)

s4 <- summary(p4)
s4

best4 <- best(p4)
plot(best4)

dsep_p4 <- as.data.frame (p4$d_sep)
dsep_p4$D.model <- NULL
dsep_p4$significant <- dsep_p4$D.p < 0.05
dsep_p4[] <- lapply(dsep_p4, function(x) if (is.list(x)) as.character(x) else x)
dsep_p4
#write.csv(dsep_p4, "modelD_dsep_results.csv", row.names = FALSE)


boot_modelD <- best(p4, boot = 1000)
coef_plot(boot_modelD)

coefD<- boot_modelD$coef
lowerD <- boot_modelD$lower
upperD <- boot_modelD$upper

rowsD <- rownames(coefD)
colsD <- colnames(coefD)

# Initialize a data frame to store results
resultsD <- data.frame(From = character(), To = character(), Estimate = numeric(), Lower = numeric(), Upper = numeric(), Significant = logical())

# Iterate through the coefficient matrix (excluding the diagonal)
for (i in 1:nrow(coefD)) {
  for (j in 1:ncol(coefD)) {
    if (i != j) { # Avoid diagonal elements (self-effects)
      from <- rowsD[i]
      to <- colsD[j]
      estimate <- coefD[i, j]
      lower_bound <- lowerD[i, j]
      upper_bound <- upperD[i, j]
      significant <- (lower_bound > 0 || upper_bound < 0) # Check if CI overlaps zero
      
      resultsD <- rbind(resultsD, data.frame(From = from, To = to, Estimate = estimate, Lower = lower_bound, Upper = upper_bound, Significant = significant))
    }
  }
}

print(resultsD)
#write.csv(resultsD, file = "modelD_coefficients.csv", row.names = FALSE)


m5 <- define_model_set(
  E = c(
    u ~ MS + BM + LS + GT + MA,
    Ne ~ GT,
    GT ~ MA + LS,
    LS ~ BM,
    MA ~ BM,
    OPG ~ BM + MA + LS + GT,
    LS ~ MA,
    MS ~ BM,
    MS ~ MA
  )
)


plot_model_set(m5)

p5<- phylo_path(m5,path1_d,time_tree_u, lower.bound = 0.000, upper.bound = 1.00)

s5 <- summary(p5)
s5

best5 <- best(p5)
plot(best5)

dsep_p5 <- as.data.frame (p5$d_sep)
dsep_p5$E.model <- NULL
dsep_p5$significant <- dsep_p5$E.p < 0.05
dsep_p5[] <- lapply(dsep_p5, function(x) if (is.list(x)) as.character(x) else x)
dsep_p5
#write.csv(dsep_p5, "modelE_dsep_results.csv", row.names = FALSE)


boot_modelE <- best(p5, boot = 1000)
coef_plot(boot_modelE)

coefE<- boot_modelE$coef
lowerE <- boot_modelE$lower
upperE <- boot_modelE$upper

rowsE <- rownames(coefE)
colsE <- colnames(coefE)

# Initialize a data frame to store results
resultsE <- data.frame(From = character(), To = character(), Estimate = numeric(), Lower = numeric(), Upper = numeric(), Significant = logical())

# Iterate through the coefficient matrix (excluding the diagonal)
for (i in 1:nrow(coefE)) {
  for (j in 1:ncol(coefE)) {
    if (i != j) { # Avoid diagonal elements (self-effects)
      from <- rowsE[i]
      to <- colsE[j]
      estimate <- coefE[i, j]
      lower_bound <- lowerE[i, j]
      upper_bound <- upperE[i, j]
      significant <- (lower_bound > 0 || upper_bound < 0) # Check if CI overlaps zero
      
      resultsE <- rbind(resultsE, data.frame(From = from, To = to, Estimate = estimate, Lower = lower_bound, Upper = upper_bound, Significant = significant))
    }
  }
}

print(resultsE)
#write.csv(resultsE, file = "modelE_coefficients.csv", row.names = FALSE)




full_models <- define_model_set(
  A = c(u ~ Ne + MS + BM + LS + OPG + GT +MA),
  B = c(u ~ Ne + MS + BM + LS + OPG + GT +MA, Ne ~ GT, GT ~ MA, GT ~ LS),
  C = c(
    u ~ Ne + MS + BM + LS + GT + MA,
    Ne ~ GT,
    GT ~ MA + LS,
    LS ~ BM,
    MA ~ BM,
    OPG ~ BM + MA + LS + GT),
  D = c(
    u ~ Ne + MS + BM + LS + GT + MA,
    Ne ~ GT,
    GT ~ MA + LS,
    LS ~ BM,
    MA ~ BM,
    OPG ~ BM + MA + LS + GT,
    LS ~ MA,
    MS ~ BM,
    MS ~ MA
  ),
  E = c(
    u ~ MS + BM + LS + GT + MA,
    Ne ~ GT,
    GT ~ MA + LS,
    LS ~ BM,
    MA ~ BM,
    OPG ~ BM + MA + LS + GT,
    LS ~ MA,
    MS ~ BM,
    MS ~ MA)
  )

fig2 <- plot_model_set(full_models)

fig2 + theme(strip.text = element_text(size = 16)) 


fit_all <- phylo_path(full_models, data = path1_d, tree = time_tree_u,lower.bound = 0.000, upper.bound = 1.00)

summary_all <- summary(fit_all)
print(summary_all)


plot_bestD <- plot(boot_modelD) + 
  ggtitle("(A) Model D Standardized Path Coefficents") +
  theme(plot.title = element_text(size = 12, hjust = 0.5))

plot_bestD

coef_df <- as.data.frame(boot_modelD$coef)
lower_df <- as.data.frame(boot_modelD$lower)
upper_df <- as.data.frame(boot_modelD$upper)

plot_data <- coef_df %>%
  rownames_to_column("From") %>%
  pivot_longer(-From, names_to = "To", values_to = "Estimate") %>%
  left_join(
    lower_df %>% rownames_to_column("From") %>%
      pivot_longer(-From, names_to = "To", values_to = "Lower"),
    by = c("From", "To")
  ) %>%
  left_join(
    upper_df %>% rownames_to_column("From") %>%
      pivot_longer(-From, names_to = "To", values_to = "Upper"),
    by = c("From", "To")
  ) %>%
  filter(From != To & !(Estimate == 0 & Lower == 0 & Upper == 0)) %>%
  mutate(
    Pair = paste(From, "→", To),
    Response = To,
    Significant = Lower > 0 | Upper < 0,
    Direction = case_when(
      Significant & Estimate > 0 ~ "Positive",
      Significant & Estimate < 0 ~ "Negative",
      TRUE ~ "Non-significant"
    )
  )

# Desired response order
response_order <- c("u", "Ne", "GT", "MS", "OPG", "LS", "MA", "BM")

# Order by response variable grouping
plot_data <- plot_data %>%
  mutate(Response = factor(Response, levels = response_order)) %>%
  arrange(Response, To, From) %>%
  mutate(Pair = factor(Pair, levels = unique(Pair)))

# Create the plot
plotDcoef <- ggplot(plot_data, aes(x = Pair, y = Estimate, color = Direction)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2) +
  scale_color_manual(
    values = c("Positive" = "blue", "Negative" = "red", "Non-significant" = "gray50")
  ) +
  theme_minimal() +
  labs(
    title = "(B) Model D Path Coefficients Ordered by Response Variable",
    x = "Life History Trait Pairing",
    y = "Standardized Path Coefficient"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    plot.title = element_text(hjust = 0.5, size = 12),
    legend.position = "none"
  )
plotDcoef

plot_bestD <- plot_bestD + 
  theme(
    plot.margin = margin(5, 10, 5, 10),
    plot.title = element_text(size = 11),
    axis.text.x = element_blank(),        # hides tick labels
    axis.ticks.x = element_blank(),       # hides tick marks
    axis.title.x = element_blank(),        # hides axis title
    axis.text.y= element_blank(),
  )
plot_bestD

fig3 <- ggarrange(
  plot_bestD, 
  plotDcoef, 
  nrow = 2,
  heights = c(1.5, 1) 
)

fig3

#### redo with Ne estimated from pi


path2_d = data.frame(d$Ne_pi,d$Modeled.rate.per.generation..m_generation_modeled., d$Generation.time..years.,
                     d$Mating.strategy..binary., 
                     d$Mass..grams.,
                     d$Lifespan.in.the.wild..years., 
                     d$Average.maturation.time,
                     d$Number.of.offspring.per.generation)
rownames(path2_d) <- d$Species.name
colnames(path2_d) <- c("Ne","u", "GT", "MS", "BM", "LS", "MA","OPG")

path2_d<- path2_d %>% drop_na(Ne)

tip2 <- c( "Orcinus_orca","Panthera_tigris","Paralichthys_olivaceus", "Rangifer_tarandus")   
tree_pi_Ne = drop.tip(time_tree_u, tip2)

name.check(tree_pi_Ne,path2_d)
path2_d<-ReorderData(tree_pi_Ne,path2_d)


mF<- define_model_set(
  F = c(u ~ Ne + MS + BM + LS + OPG + GT +MA)
)

plot_model_set(mF)

pF <- phylo_path(mF,path2_d,tree_pi_Ne, lower.bound = 0.000, upper.bound = 1.00)

sF <- summary(pF)
sF

dsep_pF <- as.data.frame (pF$d_sep)
dsep_pF
dsep_pF$F.model <- NULL
dsep_pF$significant <- dsep_pF$F.p < 0.05
dsep_pF[] <- lapply(dsep_pF, function(x) if (is.list(x)) as.character(x) else x)
#write.csv(dsep_pF, "modelF_dsep_results.csv", row.names = FALSE)

bestF <- best(pF)
plot(bestF)


mG <- define_model_set(
  G = c(u ~ Ne + MS + BM + LS + OPG + GT +MA, Ne ~ GT, GT ~ MA, GT ~ LS)
)

plot_model_set(mG)

pG <- phylo_path(mG,path2_d,tree_pi_Ne, lower.bound = 0.000, upper.bound = 1.00)

sG <- summary(pG)
sG

dsep_pG <- as.data.frame (pG$d_sep)
dsep_pG$G.model <- NULL
dsep_pG$significant <- dsep_pG$G.p < 0.05
dsep_pG[] <- lapply(dsep_pG, function(x) if (is.list(x)) as.character(x) else x)
dsep_pG

#write.csv(dsep_pG, "modelG_dsep_results.csv", row.names = FALSE)


bestG <- best(pG)
plot(bestG)


mH <- define_model_set(
  H = c(
    u ~ Ne + MS + BM + LS + GT + MA,
    Ne ~ GT,
    GT ~ MA + LS,
    LS ~ BM,
    MA ~ BM,
    OPG ~ BM + MA + LS + GT
  )
)

plot_model_set(mH)

pH<- phylo_path(mH,path2_d,tree_pi_Ne, lower.bound = 0.000, upper.bound = 1.00)

sH <- summary(pH)
sH

dsep_pH <- as.data.frame (pH$d_sep)
dsep_pH$H.model <- NULL
dsep_pH$significant <- dsep_pH$H.p < 0.05
dsep_pH[] <- lapply(dsep_pH, function(x) if (is.list(x)) as.character(x) else x)
dsep_pH
#write.csv(dsep_pH, "modelH_dsep_results.csv", row.names = FALSE)

bestH <- best(pH)
plot(bestH)


mI <- define_model_set(
  I = c(
    u ~ Ne + MS + BM + LS + GT + MA,
    Ne ~ GT,
    GT ~ MA + LS,
    LS ~ BM,
    MA ~ BM,
    OPG ~ BM + MA + LS + GT,
    LS ~ MA,
    MS ~ BM,
    MS ~ MA
  )
)


plot_model_set(mI)

pI<- phylo_path(mI,path2_d,tree_pi_Ne, lower.bound = 0.000, upper.bound = 1.00)

sI <- summary(pI)
sI

bestI <- best(pI)
plot(bestI)

dsep_pI <- as.data.frame (pI$d_sep)
dsep_pI$I.model <- NULL
dsep_pI$significant <- dsep_pI$I.p < 0.05
dsep_pI[] <- lapply(dsep_pI, function(x) if (is.list(x)) as.character(x) else x)
dsep_pI
#write.csv(dsep_pI, "modelI_dsep_results.csv", row.names = FALSE)


boot_modelI <- best(pI, boot = 1000)
coef_plot(boot_modelI)

coefI<- boot_modelI$coef
lowerI <- boot_modelI$lower
upperI <- boot_modelI$upper

rowsI <- rownames(coefI)
colsI <- colnames(coefI)

# Initialize a data frame to store results
resultsI <- data.frame(From = character(), To = character(), Estimate = numeric(), Lower = numeric(), Upper = numeric(), Significant = logical())

# Iterate through the coefficient matrix (excluding the diagonal)
for (i in 1:nrow(coefI)) {
  for (j in 1:ncol(coefI)) {
    if (i != j) { # Avoid diagonal elements (self-effects)
      from <- rowsI[i]
      to <- colsI[j]
      estimate <- coefI[i, j]
      lower_bound <- lowerI[i, j]
      upper_bound <- upperI[i, j]
      significant <- (lower_bound > 0 || upper_bound < 0) # Check if CI overlaps zero
      
      resultsI <- rbind(resultsI, data.frame(From = from, To = to, Estimate = estimate, Lower = lower_bound, Upper = upper_bound, Significant = significant))
    }
  }
}

print(resultsI)
#write.csv(resultsI, file = "modelI_coefficients.csv", row.names = FALSE)


mJ <- define_model_set(
  J = c(
    u ~ MS + BM + LS + GT + MA,
    Ne ~ GT,
    GT ~ MA + LS,
    LS ~ BM,
    MA ~ BM,
    OPG ~ BM + MA + LS + GT,
    LS ~ MA,
    MS ~ BM,
    MS ~ MA
  )
)


plot_model_set(mJ)

pJ<- phylo_path(mJ,path2_d,tree_pi_Ne, lower.bound = 0.000, upper.bound = 1.00)

sJ <- summary(pJ)
sJ

bestJ <- best(pJ)
plot(bestJ)

dsep_pJ <- as.data.frame (pJ$d_sep)
dsep_pJ$J.model <- NULL
dsep_pJ$significant <- dsep_pJ$J.p < 0.05
dsep_pJ[] <- lapply(dsep_pJ, function(x) if (is.list(x)) as.character(x) else x)
dsep_pJ
#write.csv(dsep_pJ, "modelJ_dsep_results.csv", row.names = FALSE)


boot_modelJ <- best(pJ, boot = 1000)
coef_plot(boot_modelJ)

coefJ<- boot_modelJ$coef
lowerJ <- boot_modelJ$lower
upperJ <- boot_modelJ$upper

rowsJ <- rownames(coefJ)
colsJ <- colnames(coefJ)

# Initialize a data frame to store results
resultsJ <- data.frame(From = character(), To = character(), Estimate = numeric(), Lower = numeric(), Upper = numeric(), Significant = logical())

# Iterate through the coefficient matrix (excluding the diagonal)
for (i in 1:nrow(coefJ)) {
  for (j in 1:ncol(coefJ)) {
    if (i != j) { # Avoid diagonal elements (self-effects)
      from <- rowsJ[i]
      to <- colsJ[j]
      estimate <- coefJ[i, j]
      lower_bound <- lowerJ[i, j]
      upper_bound <- upperJ[i, j]
      significant <- (lower_bound > 0 || upper_bound < 0) # Check if CI overlaps zero
      
      resultsJ <- rbind(resultsJ, data.frame(From = from, To = to, Estimate = estimate, Lower = lower_bound, Upper = upper_bound, Significant = significant))
    }
  }
}

print(resultsJ)
#write.csv(resultsJ, file = "modelJ_coefficients.csv", row.names = FALSE)


full_models2 <- define_model_set(
  F = c(u ~ Ne + MS + BM + LS + OPG + GT +MA),
  G = c(u ~ Ne + MS + BM + LS + OPG + GT +MA, Ne ~ GT, GT ~ MA, GT ~ LS),
  H = c(
    u ~ Ne + MS + BM + LS + GT + MA,
    Ne ~ GT,
    GT ~ MA + LS,
    LS ~ BM,
    MA ~ BM,
    OPG ~ BM + MA + LS + GT),
  I = c(
    u ~ Ne + MS + BM + LS + GT + MA,
    Ne ~ GT,
    GT ~ MA + LS,
    LS ~ BM,
    MA ~ BM,
    OPG ~ BM + MA + LS + GT,
    LS ~ MA,
    MS ~ BM,
    MS ~ MA
  ),
  J = c(
    u ~ MS + BM + LS + GT + MA,
    Ne ~ GT,
    GT ~ MA + LS,
    LS ~ BM,
    MA ~ BM,
    OPG ~ BM + MA + LS + GT,
    LS ~ MA,
    MS ~ BM,
    MS ~ MA)
)

plot_model_set(full_models2)


fit_all2 <- phylo_path(full_models2, data = path2_d, tree = tree_pi_Ne,lower.bound = 0.000, upper.bound = 1.00)

summary_all2 <- summary(fit_all2)
print(summary_all2)


##### explore GS ~ Ne + u

path3_d = data.frame(d$Ne..harmonic.mean..psmc.,d$Modeled.rate.per.generation..m_generation_modeled., d$Generation.time..years.,
                     d$Mating.strategy..binary., 
                     d$Mass..grams.,
                     d$Lifespan.in.the.wild..years., 
                     d$Average.maturation.time,
                     d$Number.of.offspring.per.generation, 
                     d$genome_size)
rownames(path3_d) <- d$Species.name
colnames(path3_d) <- c("Ne","u", "GT", "MS", "BM", "LS", "MA","OPG","GS")

mK<- define_model_set(
  K = c(GS ~ Ne + MS + BM + LS + OPG + GT +MA + u)
)

plot_model_set(mK)


pK <- phylo_path(mK,path3_d,time_tree_u, lower.bound = 0.000, upper.bound = 1.00)

sK <- summary(pK)
sK

dsep_pK <- as.data.frame (pK$d_sep)
dsep_pK
dsep_pK$K.model <- NULL
dsep_pK$significant <- dsep_pK$K.p < 0.05
dsep_pK[] <- lapply(dsep_pK, function(x) if (is.list(x)) as.character(x) else x)
#write.csv(dsep_pK, "modelK_dsep_results.csv", row.names = FALSE)

bestK <- best(pK)
plot(bestK)


### add 3 bio grounded links MA -> GT, GT -> Ne, LS -> GT plus 3 top effects on u GT -> u, MA -> u, MS -> u
mL<- define_model_set(
  L = c(GS ~ Ne + MS + BM + LS + OPG + GT +MA + u, GT ~ MA + LS, Ne ~ GT, u ~ GT + MA + MS))

plot_model_set(mL)

pL <- phylo_path(mL,path3_d,time_tree_u, lower.bound = 0.000, upper.bound = 1.00)

sL <- summary(pL)
sL

dsep_pL <- as.data.frame (pL$d_sep)
dsep_pL
dsep_pL$L.model <- NULL
dsep_pL$significant <- dsep_pL$L.p < 0.05
dsep_pL[] <- lapply(dsep_pL, function(x) if (is.list(x)) as.character(x) else x)
#write.csv(dsep_pL, "modelL_dsep_results.csv", row.names = FALSE)


bestL <- best(pL)
plot(bestL)



#plus other sig effects from model D: BM -> LS, MA -> LS, BM -> MA, MA -> MS, BM -> MS, & OPG ~ GT + BM + MA, remove OPG -> GS path
mM<- define_model_set(
  M = c(GS ~ Ne + MS + BM + LS + GT +MA + u, GT ~ MA + LS, Ne ~ GT, u ~ GT + MA + MS, LS ~ BM + MA, MA ~ BM, MS ~ MA + BM, OPG ~ GT + BM + MA))

plot_model_set(mM)

pM <- phylo_path(mM,path3_d,time_tree_u, lower.bound = 0.000, upper.bound = 1.00)

sM <- summary(pM)
sM

dsep_pM <- as.data.frame (pM$d_sep)
dsep_pM
dsep_pM$M.model <- NULL
dsep_pM$significant <- dsep_pM$M.p < 0.05
dsep_pM[] <- lapply(dsep_pM, function(x) if (is.list(x)) as.character(x) else x)
#write.csv(dsep_pM, "modelM_dsep_results.csv", row.names = FALSE)

boot_modelM <- best(pM, boot = 1000)


coef_plot(boot_modelM)

coefM<- boot_modelM$coef
lowerM <- boot_modelM$lower
upperM <- boot_modelM$upper

rowsM <- rownames(coefM)
colsM <- colnames(coefM)

# Initialize a data frame to store results
resultsM <- data.frame(From = character(), To = character(), Estimate = numeric(), Lower = numeric(), Upper = numeric(), Significant = logical())

# Iterate through the coefficient matrix (excluding the diagonal)
for (i in 1:nrow(coefM)) {
  for (j in 1:ncol(coefM)) {
    if (i != j) { # Avoid diagonal elements (self-effects)
      from <- rowsM[i]
      to <- colsM[j]
      estimate <- coefM[i, j]
      lower_bound <- lowerM[i, j]
      upper_bound <- upperM[i, j]
      significant <- (lower_bound > 0 || upper_bound < 0) # Check if CI overlaps zero
      
      resultsM <- rbind(resultsM, data.frame(From = from, To = to, Estimate = estimate, Lower = lower_bound, Upper = upper_bound, Significant = significant))
    }
  }
}

print(resultsM)
#write.csv(resultsM, file = "modelM_coefficients.csv", row.names = FALSE)


plot_bestM <- plot(boot_modelM) + 
  ggtitle("(A) Model M Standardized Path Coefficents") +
  theme(plot.title = element_text(size = 12, hjust = 0.5))

plot_bestM


coefM_df <- as.data.frame(boot_modelM$coef)
lowerM_df <- as.data.frame(boot_modelM$lower)
upperM_df <- as.data.frame(boot_modelM$upper)

plot_dataM <- coefM_df %>%
  rownames_to_column("From") %>%
  pivot_longer(-From, names_to = "To", values_to = "Estimate") %>%
  left_join(
    lowerM_df %>% rownames_to_column("From") %>%
      pivot_longer(-From, names_to = "To", values_to = "Lower"),
    by = c("From", "To")
  ) %>%
  left_join(
    upperM_df %>% rownames_to_column("From") %>%
      pivot_longer(-From, names_to = "To", values_to = "Upper"),
    by = c("From", "To")
  ) %>%
  filter(From != To & !(Estimate == 0 & Lower == 0 & Upper == 0)) %>%
  mutate(
    Pair = paste(From, "→", To),
    Response = To,
    Significant = Lower > 0 | Upper < 0,
    Direction = case_when(
      Significant & Estimate > 0 ~ "Positive",
      Significant & Estimate < 0 ~ "Negative",
      TRUE ~ "Non-significant"
    )
  )

# Desired response order
response_orderM <- c("GS","u", "Ne", "GT", "MS", "OPG", "LS", "MA", "BM")

# Order by response variable grouping
plot_dataM <- plot_dataM %>%
  mutate(Response = factor(Response, levels = response_order)) %>%
  arrange(Response, To, From) %>%
  mutate(Pair = factor(Pair, levels = unique(Pair)))

# Create the plot
plotMcoef <- ggplot(plot_dataM, aes(x = Pair, y = Estimate, color = Direction)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2) +
  scale_color_manual(
    values = c("Positive" = "blue", "Negative" = "red", "Non-significant" = "gray50")
  ) +
  theme_minimal() +
  labs(
    title = "(B) Model M Path Coefficients Ordered by Response Variable",
    x = "Life History Trait Pairing",
    y = "Standardized Path Coefficient"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    plot.title = element_text(hjust = 0.5, size = 12),
    legend.position = "none"
  )
plotMcoef

plot_bestM <- plot_bestM + 
  theme(
    plot.margin = margin(5, 10, 5, 10),
    plot.title = element_text(size = 11),
    axis.text.x = element_blank(),        # hides tick labels
    axis.ticks.x = element_blank(),       # hides tick marks
    axis.title.x = element_blank(),        # hides axis title
    axis.text.y= element_blank(),
  )
plot_bestM

fig4 <- ggarrange(
  plot_bestM, 
  plotMcoef, 
  nrow = 2,
  heights = c(1.5, 1) 
)

fig4



mN<- define_model_set(
  N = c(GS ~  MS + BM + LS + GT +MA , GT ~ MA + LS, Ne ~ GT, u ~ GT + MA + MS, LS ~ BM + MA, MA ~ BM, MS ~ MA + BM, OPG ~ GT + BM + MA))

plot_model_set(mN)

pN<- phylo_path(mN,path3_d,time_tree_u, lower.bound = 0.000, upper.bound = 1.00)

sN <- summary(pN)
sN

dsep_pN <- as.data.frame (pN$d_sep)
dsep_pN
dsep_pN$N.model <- NULL
dsep_pN$significant <- dsep_pN$N.p < 0.05
dsep_pN[] <- lapply(dsep_pN, function(x) if (is.list(x)) as.character(x) else x)
#write.csv(dsep_pN, "modelN_dsep_results.csv", row.names = FALSE)

boot_modelN <- best(pN, boot = 1000)
coef_plot(boot_modelN)

coefN<- boot_modelN$coef
lowerN <- boot_modelN$lower
upperN <- boot_modelN$upper

rowsN <- rownames(coefN)
colsN <- colnames(coefN)

# Initialize a data frame to store results
resultsN <- data.frame(From = character(), To = character(), Estimate = numeric(), Lower = numeric(), Upper = numeric(), Significant = logical())

# Iterate through the coefficient matrix (excluding the diagonal)
for (i in 1:nrow(coefN)) {
  for (j in 1:ncol(coefN)) {
    if (i != j) { # Avoid diagonal elements (self-effects)
      from <- rowsN[i]
      to <- colsN[j]
      estimate <- coefN[i, j]
      lower_bound <- lowerN[i, j]
      upper_bound <- upperN[i, j]
      significant <- (lower_bound > 0 || upper_bound < 0) # Check if CI overlaps zero
      
      resultsN <- rbind(resultsN, data.frame(From = from, To = to, Estimate = estimate, Lower = lower_bound, Upper = upper_bound, Significant = significant))
    }
  }
}

print(resultsN)
#write.csv(resultsN, file = "modelN_coefficients.csv", row.names = FALSE)



full_models3 <- define_model_set(
  K = c(GS ~ Ne + MS + BM + LS + OPG + GT +MA + u),
  L = c(GS ~ Ne + MS + BM + LS + OPG + GT +MA + u, GT ~ MA + LS, Ne ~ GT, u ~ GT + MA + MS),
  M = c(GS ~ Ne + MS + BM + LS + GT +MA + u, GT ~ MA + LS, Ne ~ GT, u ~ GT + MA + MS, LS ~ BM + MA, MA ~ BM, MS ~ MA + BM, OPG ~ GT + BM + MA),
  N = c(GS ~  MS + BM + LS + GT +MA , GT ~ MA + LS, Ne ~ GT, u ~ GT + MA + MS, LS ~ BM + MA, MA ~ BM, MS ~ MA + BM, OPG ~ GT + BM + MA))

plot_model_set(full_models3)

#fixfig2 <- plot_model_set(full_models3, rotation = 120)
#fixfig2

#plot_data2 <- ggplot_build(fixfig)$data[[2]]
#node_positions2 <- unique(plot_data[c("x", "y", "label")])

#custom_layout2 <- data.frame(
  #name = node_positions$label,
  #x = node_positions$x,
  #y = node_positions$y
#)

# Print the custom layout to see the default coordinates
#print(custom_layout2)

#custom_layout2 <- data.frame(
  #name = c("u", "MA", "GT", "OPG", "LS", "BM", "MS", "Ne", "GS"),
 # x = c(-0.9378208, 0.2393997, 0.5485132, 1.4732415, 1.5817125, 0.6203072, -1.2395277, -0.6375485, -1.1427627),
 # y = c(1.4703548, 0.4387891, -1.2430324, -0.3938475, 0.4235073, 1.2552865, 0.5182317, -1.8405801, -0.2748826)
#)


#plot_model_set(full_models3,  manual_layout = custom_layout)

fit_all3 <- phylo_path(full_models3, data = path3_d, tree = time_tree_u,lower.bound = 0.000, upper.bound = 1.00)

summary_all3 <- summary(fit_all3)
print(summary_all3)





##### check GS ~ Neu with pi Ne


path4_d = data.frame(d$Ne_pi   ,d$Modeled.rate.per.generation..m_generation_modeled., d$Generation.time..years.,
                     d$Mating.strategy..binary., 
                     d$Mass..grams.,
                     d$Lifespan.in.the.wild..years., 
                     d$Average.maturation.time,
                     d$Number.of.offspring.per.generation, 
                     d$genome_size)
rownames(path4_d) <- d$Species.name
colnames(path4_d) <- c("Ne","u", "GT", "MS", "BM", "LS", "MA","OPG","GS")

path4_d<- path4_d %>% drop_na(Ne)

name.check(tree_pi_Ne,path4_d)
path4_d<-ReorderData(tree_pi_Ne,path4_d)



mO<- define_model_set(
  O = c(GS ~ Ne + MS + BM + LS + OPG + GT +MA + u)
)

plot_model_set(mO)


pO <- phylo_path(mO,path4_d,tree_pi_Ne, lower.bound = 0.000, upper.bound = 1.00)

sO <- summary(pO)
sO

dsep_pO <- as.data.frame (pO$d_sep)
dsep_pO
dsep_pO$O.model <- NULL
dsep_pO$significant <- dsep_pO$O.p < 0.05
dsep_pO[] <- lapply(dsep_pO, function(x) if (is.list(x)) as.character(x) else x)
#write.csv(dsep_pO, "modelO_dsep_results.csv", row.names = FALSE)

bestO <- best(pO)
plot(bestO)

mP<- define_model_set(
  P = c(GS ~ Ne + MS + BM + LS + OPG + GT +MA + u, GT ~ MA + LS, Ne ~ GT, u ~ GT + MA + MS))

plot_model_set(mP)

pP <- phylo_path(mP,path4_d,tree_pi_Ne, lower.bound = 0.000, upper.bound = 1.00)

sP <- summary(pP)
sP

dsep_pP <- as.data.frame (pP$d_sep)
dsep_pP
dsep_pP$P.model <- NULL
dsep_pP$significant <- dsep_pP$P.p < 0.05
dsep_pP[] <- lapply(dsep_pP, function(x) if (is.list(x)) as.character(x) else x)
#write.csv(dsep_pP, "modelP_dsep_results.csv", row.names = FALSE)


bestP <- best(pP)
plot(bestP)

mQ<- define_model_set(
  Q = c(GS ~ Ne + MS + BM + LS + GT +MA + u, GT ~ MA + LS, Ne ~ GT, u ~ GT + MA + MS, LS ~ BM + MA, MA ~ BM, MS ~ MA + BM, OPG ~ GT + BM + MA))

plot_model_set(mQ)

pQ <- phylo_path(mQ,path4_d,tree_pi_Ne, lower.bound = 0.000, upper.bound = 1.00)

sQ <- summary(pQ)
sQ

dsep_pQ <- as.data.frame (pQ$d_sep)
dsep_pQ
dsep_pQ$Q.model <- NULL
dsep_pQ$significant <- dsep_pQ$Q.p < 0.05
dsep_pQ[] <- lapply(dsep_pQ, function(x) if (is.list(x)) as.character(x) else x)
#write.csv(dsep_pQ, "modelQ_dsep_results.csv", row.names = FALSE)

boot_modelQ <- best(pQ, boot = 1000)
coef_plot(boot_modelQ)

coefQ<- boot_modelQ$coef
lowerQ <- boot_modelQ$lower
upperQ <- boot_modelQ$upper

rowsQ <- rownames(coefQ)
colsQ <- colnames(coefQ)

# Initialize a data frame to store results
resultsQ <- data.frame(From = character(), To = character(), Estimate = numeric(), Lower = numeric(), Upper = numeric(), Significant = logical())

# Iterate through the coefficient matrix (excluding the diagonal)
for (i in 1:nrow(coefQ)) {
  for (j in 1:ncol(coefQ)) {
    if (i != j) { # Avoid diagonal elements (self-effects)
      from <- rowsQ[i]
      to <- colsQ[j]
      estimate <- coefQ[i, j]
      lower_bound <- lowerQ[i, j]
      upper_bound <- upperQ[i, j]
      significant <- (lower_bound > 0 || upper_bound < 0) # Check if CI overlaps zero
      
      resultsQ <- rbind(resultsQ, data.frame(From = from, To = to, Estimate = estimate, Lower = lower_bound, Upper = upper_bound, Significant = significant))
    }
  }
}

print(resultsQ)
#write.csv(resultsQ, file = "modelQ_coefficients.csv", row.names = FALSE)


mR<- define_model_set(
  R = c(GS ~  MS + BM + LS + GT +MA , GT ~ MA + LS, Ne ~ GT, u ~ GT + MA + MS, LS ~ BM + MA, MA ~ BM, MS ~ MA + BM, OPG ~ GT + BM + MA))

plot_model_set(mR)

pR<- phylo_path(mR,path4_d,tree_pi_Ne, lower.bound = 0.000, upper.bound = 1.00)

sR <- summary(pR)
sR

dsep_pR <- as.data.frame (pR$d_sep)
dsep_pR
dsep_pR$R.model <- NULL
dsep_pR$significant <- dsep_pR$R.p < 0.05
dsep_pR[] <- lapply(dsep_pN, function(x) if (is.list(x)) as.character(x) else x)
#write.csv(dsep_pR, "modelR_dsep_results.csv", row.names = FALSE)

boot_modelR <- best(pR, boot = 1000)
coef_plot(boot_modelR)

coefR<- boot_modelR$coef
lowerR <- boot_modelR$lower
upperR <- boot_modelR$upper

rowsR <- rownames(coefR)
colsR <- colnames(coefR)

# Initialize a data frame to store results
resultsR <- data.frame(From = character(), To = character(), Estimate = numeric(), Lower = numeric(), Upper = numeric(), Significant = logical())

# Iterate through the coefficient matrix (excluding the diagonal)
for (i in 1:nrow(coefR)) {
  for (j in 1:ncol(coefR)) {
    if (i != j) { # Avoid diagonal elements (self-effects)
      from <- rowsR[i]
      to <- colsR[j]
      estimate <- coefR[i, j]
      lower_bound <- lowerR[i, j]
      upper_bound <- upperR[i, j]
      significant <- (lower_bound > 0 || upper_bound < 0) # Check if CI overlaps zero
      
      resultsR <- rbind(resultsR, data.frame(From = from, To = to, Estimate = estimate, Lower = lower_bound, Upper = upper_bound, Significant = significant))
    }
  }
}

print(resultsR)
#write.csv(resultsR, file = "modelR_coefficients.csv", row.names = FALSE)



full_models4 <- define_model_set(
  O = c(GS ~ Ne + MS + BM + LS + OPG + GT +MA + u),
  P = c(GS ~ Ne + MS + BM + LS + OPG + GT +MA + u, GT ~ MA + LS, Ne ~ GT, u ~ GT + MA + MS),
  Q = c(GS ~ Ne + MS + BM + LS + GT +MA + u, GT ~ MA + LS, Ne ~ GT, u ~ GT + MA + MS, LS ~ BM + MA, MA ~ BM, MS ~ MA + BM, OPG ~ GT + BM + MA),
  R = c(GS ~  MS + BM + LS + GT +MA , GT ~ MA + LS, Ne ~ GT, u ~ GT + MA + MS, LS ~ BM + MA, MA ~ BM, MS ~ MA + BM, OPG ~ GT + BM + MA))



plot_model_set(full_models4)


fit_all4 <- phylo_path(full_models4, data = path4_d, tree= tree_pi_Ne,lower.bound = 0.000, upper.bound = 1.00)

summary_all4 <- summary(fit_all4)
print(summary_all4)


####### Neu ~ GS PGLS

d1_vars <- c("Species.name","Taxonomic.group","Modeled.rate.per.generation..m_generation_modeled." ,"genome_size","Ne..harmonic.mean..psmc.")
d1 <- d[d1_vars]

name.check(time_tree_u,d1)

d1 <-ReorderData(time_tree_u,d1)


vert1 <- comparative.data(phy=time_tree_u, data=d1, names.col = Species.name, vcv = TRUE, na.omit = TRUE, warn.dropped = TRUE)

GS_Neu_pgls <- pgls( genome_size~ Ne..harmonic.mean..psmc. * Modeled.rate.per.generation..m_generation_modeled., data=vert1, lambda = "ML")
summary(GS_Neu_pgls)
d1$Ne_mu_log <- d1$Ne..harmonic.mean..psmc. + d1$Modeled.rate.per.generation..m_generation_modeled.  # Correct: log10(Ne × μ)

coef <- coefficients(GS_Neu_pgls)
int <- coef[1]        # Intercept
b_Ne <- coef[2]       # Effect of Ne
b_u <- coef[3]        # Effect of μ
b_inter <- coef[4]    # Interaction term

# Hold u constant at its mean for Ne effect
u_mean <- mean(d1$Modeled.rate.per.generation..m_generation_modeled., na.rm = TRUE)
d1$fit_Ne <- int + b_Ne * d1$Ne..harmonic.mean..psmc. + b_u * u_mean + b_inter * d1$Ne..harmonic.mean..psmc. * u_mean

# Hold Ne constant at its mean for μ effect
Ne_mean <- mean(d1$Ne..harmonic.mean..psmc., na.rm = TRUE)
d1$fit_u <- int + b_Ne * Ne_mean + b_u * d1$Modeled.rate.per.generation..m_generation_modeled. + b_inter * Ne_mean * d1$Modeled.rate.per.generation..m_generation_modeled.

# Use full prediction for Neμ
d1$fit_Ne_mu <- int + b_Ne * d1$Ne..harmonic.mean..psmc. + b_u * d1$Modeled.rate.per.generation..m_generation_modeled. + b_inter * (d1$Ne..harmonic.mean..psmc. * d1$Modeled.rate.per.generation..m_generation_modeled.)



r1 <-ggplot(d1) +
  geom_point(aes(x = Ne_mu_log, y = genome_size, color = Taxonomic.group)) +
  geom_line(aes(x = Ne_mu_log, y = fit_Ne_mu), color = "black", linetype = "solid") +
  geom_line(aes(x = Ne_mu_log, y = fit_Ne), color = "blue", linetype = "dashed") +
  geom_line(aes(x = Ne_mu_log, y = fit_u), color = "red", linetype = "dotted") +
  xlab("log₁₀(Ne × μ)") +
  ylab("Genome size (log₁₀ bp)") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), 
        legend.position = "none")


r1

d2_vars <- c("Species.name","Taxonomic.group","Modeled.rate.per.generation..m_generation_modeled." ,"genome_size","Ne_pi")
d2 <- d[d2_vars]

d2<- d2 %>% drop_na(Ne_pi)

d2 <-ReorderData(tree_pi_Ne,d2)




vert2 <- comparative.data(phy=tree_pi_Ne, data=d2, names.col = Species.name, vcv = TRUE, na.omit = TRUE, warn.dropped = TRUE)

GS_piNeu_pgls <- pgls( genome_size~ Modeled.rate.per.generation..m_generation_modeled.*Ne_pi, data=vert2, lambda = "ML")
summary(GS_piNeu_pgls)



##show LP


d3_vars <- c("Species.name","Taxonomic.group","diversity" ,"Nc", "expected_pi")
d3 <- d[d3_vars]

d3 <- d3 %>% drop_na(Nc)
d3 <- d3 %>% drop_na(diversity)


tip3 <- c("Ailurus_fulgens","Amphiprion_ocellaris","Aptenodytes_forsteri","Betta_splendens","Canis_lupus_familiaris","Chauna_torquata","Chrysemys_picta","Coleonyx_brevis","Eublepharis_macularius","Felis_catus","Fukomys_damarensis","Larus_marinus","Moschus_berezovskii","Neovison_vison","Orcinus_orca","Phoenicopterus_roseus","Platalea_ajaja","Pogona_vitticeps","Rhea_pennata","Rousettus_aegyptiacus","Saxicola_maurus","Sphaerodactylus_inigoi","Syngnathus_scovelli","Thamnophis_sirtalis","Tursiops_truncatus","Vicugna_pacos","Arctocephalus_gazella","Larus_argentatus", "Tupaia_belangeri" ,"Panthera_tigris" , "Paralichthys_olivaceus","Rangifer_tarandus")
time_tree_Nc_pi <-drop.tip(time_tree_u,tip3)
name.check(time_tree_Nc_pi,d3)


d3 <-ReorderData(time_tree_Nc_pi,d3)

vert3 <- comparative.data(phy=time_tree_Nc_pi, data=d3, names.col = Species.name, vcv = TRUE, na.omit = TRUE, warn.dropped = TRUE)

pi_Nc_pgls <- pgls(diversity ~Nc, data=vert3, lambda = "ML")
summary(pi_Nc_pgls)

coef3 <- coefficients(pi_Nc_pgls)
int3 <-coef3[1] 
B3 <- coef3[2]

r3 <- ggplot(data = d3, aes(x = Nc, y = diversity)) +
  xlab("Census population size (log10)") + 
  ylab("nucleotide silent-site diversity (log10)") +
  geom_point() +
  scale_color_manual(values = my_colors) +
  geom_abline(intercept = int3, slope = B3) +
  geom_line(data = d3, aes(x = Nc, y = expected_pi), color = "red") +
  geom_point(aes(color = Taxonomic.group)) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), 
        legend.position = "none")

r3


#show Ne/Nc/Nc

d4_vars <- c("Species.name","Taxonomic.group","Nc", "piNe_Nc","Ne_Nc")
d4 <- d[d4_vars]

d4 <- d4 %>% drop_na(Nc)
d4 <- d4 %>% drop_na(piNe_Nc)


name.check(time_tree_Nc_pi,d4)


d4 <-ReorderData(time_tree_Nc_pi,d4)

vert4 <- comparative.data(phy=time_tree_Nc_pi, data=d4, names.col = Species.name, vcv = TRUE, na.omit = TRUE, warn.dropped = TRUE)

NeNc_Nc_pgls <- pgls(Ne_Nc ~Nc, data=vert4, lambda = "ML")
summary(NeNc_Nc_pgls)

coef4 <- coefficients(NeNc_Nc_pgls)
int4 <-coef4[1] 
B4 <- coef4[2]

r4 <- ggplot(data = d4, aes(x = Nc, y = Ne_Nc)) +
  xlab("Census population size (log10)") + 
  ylab("Ne/Nc (log10)") +
  geom_point() +
  scale_color_manual(values = my_colors) +
  geom_abline(intercept = int4, slope = B4) +
  geom_point(aes(color = Taxonomic.group)) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), 
        legend.position = "none")

r4

piNeNc_Nc_pgls <- pgls (piNe_Nc ~ Nc, data = vert4, lambda = "ML")
summary(piNeNc_Nc_pgls)
coef5 <- coefficients(piNeNc_Nc_pgls)
int5 <-coef5[1] 
B5 <- coef5[2]


r5 <- ggplot(data = d4, aes(x = Nc, y = piNe_Nc)) +
  xlab("Census population size (log10)") + 
  ylab("pi-based Ne/Nc (log10)") +
  geom_point() +
  scale_color_manual(values = my_colors) +
  geom_abline(intercept = int5, slope = B5) +
  geom_point(aes(color = Taxonomic.group)) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), 
        legend.position = "none")

r5


#### Ne/Nc ~ u
d5_vars <- c("Species.name","Taxonomic.group","Modeled.rate.per.generation..m_generation_modeled.", "piNe_Nc","Ne_Nc")
d5 <- d[d5_vars]

d5 <- d5 %>% drop_na(piNe_Nc)



name.check(time_tree_Nc_pi,d5)


d5 <-ReorderData(time_tree_Nc_pi,d5)

vert5 <- comparative.data(phy=time_tree_Nc_pi, data=d5, names.col = Species.name, vcv = TRUE, na.omit = TRUE, warn.dropped = TRUE)

NeNc_u_pgls <- pgls(Ne_Nc ~Modeled.rate.per.generation..m_generation_modeled., data=vert5, lambda = "ML")
summary(NeNc_u_pgls)

coef6 <- coefficients(NeNc_u_pgls)
int6 <-coef6[1] 
B6 <- coef6[2]

r6 <- ggplot(data = d5, aes(x = Modeled.rate.per.generation..m_generation_modeled., y = Ne_Nc)) +
  xlab("Mutation rate per site per generation (log10)") + 
  ylab("Ne/Nc (log10)") +
  geom_point() +
  scale_color_manual(values = my_colors) +
  geom_abline(intercept = int6, slope = B6) +
  geom_point(aes(color = Taxonomic.group)) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), 
        legend.position = "none")

r6

piNe_Nc_u_pgls <- pgls (piNe_Nc ~ Modeled.rate.per.generation..m_generation_modeled., data = vert5, lambda = "ML")
summary(piNe_Nc_u_pgls)
coef7 <- coefficients(piNe_Nc_u_pgls)
int7 <-coef7[1] 
B7 <- coef7[2]


r7 <- ggplot(data = d5, aes(x = Modeled.rate.per.generation..m_generation_modeled., y = piNe_Nc)) +
  xlab("Mutation rate per site per generation (log10)") + 
  ylab("pi-based Ne/Nc (log10)") +
  geom_point() +
  scale_color_manual(values = my_colors) +
  geom_abline(intercept = int7, slope = B7) +
  geom_point(aes(color = Taxonomic.group)) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), 
        legend.position = "none")

r7







figure6 <- ggarrange(r1, r4,r3, r6,
                     labels = c("A", "B", "C","D"),
                     common.legend = TRUE, legend = "bottom",
                     ncol = 2, nrow = 2)
figure6



###Ne ~ Nc


d6_vars <- c("Species.name","Taxonomic.group","Ne..harmonic.mean..psmc.", "Ne_pi", "Nc")
d6 <- d[d6_vars]

d6 <- d6 %>% drop_na(Nc)
d6 <- d6 %>% drop_na(Ne_pi)


name.check(time_tree_Nc_pi,d6)


d6 <-ReorderData(time_tree_Nc_pi,d6)

vert6 <- comparative.data(phy=time_tree_Nc_pi, data=d6, names.col = Species.name, vcv = TRUE, na.omit = TRUE, warn.dropped = TRUE)

Ne_Nc_pgls <- pgls(Ne..harmonic.mean..psmc.~ Nc, data=vert6, lambda = "ML")
summary(Ne_Nc_pgls)

coef8 <- coefficients(Ne_Nc_pgls)
int8 <-coef8[1] 
B8 <- coef8[2]

r8 <- ggplot(data = d6, aes(x = Nc ,y = Ne..harmonic.mean..psmc.)) +
  xlab("Census population size (log10)") + 
  ylab("PSMC effective population size (log10)") +
  geom_point() +
  scale_color_manual(values = my_colors) +
  geom_abline(intercept = int8, slope = B8) +
  geom_point(aes(color = Taxonomic.group)) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), 
        legend.position = "none")

r8


piNe_Nc_pgls <- pgls (Ne_pi~ Nc , data = vert6, lambda = "ML")
summary(piNe_Nc_pgls)

coef9 <- coefficients(piNe_Nc_pgls)
int9 <-coef9[1] 
B9 <- coef9[2]


r9 <- ggplot(data = d6, aes(x = Nc, y = Ne_pi)) +
  xlab("Census population size (log10)") + 
  ylab("pi-based Ne (log10)") +
  geom_point() +
  scale_color_manual(values = my_colors) +
  geom_abline(intercept = int9, slope = B9) +
  geom_point(aes(color = Taxonomic.group)) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), 
        legend.position = "none")

r9







figure6 <- ggarrange(r1, r8, r3,r4, r6,
                     labels = c("A", "B", "C","D","E"),
                     common.legend = TRUE, legend = "bottom",
                     ncol = 3, nrow = 2)
figure6





























