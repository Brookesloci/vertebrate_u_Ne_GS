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
library(car)
library(phylolm)
library(visreg)
library(purrr)
library(MASS)
library(stats)
library(MuMIn)

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
    scale_color_manual(values = my_colors) +  
    theme_minimal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                            panel.background = element_blank(), axis.line = element_line(colour = "black"))
  plot_list <- c(plot_list, list(p))
}


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

plot_A <-plot_model_set(m1)

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
print(best1)


m2 <- define_model_set(
  B = c(u ~ Ne + MS + BM + LS + OPG + GT +MA, Ne ~ GT, GT ~ MA, GT ~ LS)
)

plot_B <-plot_model_set(m2)

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

coef_df <- as.data.frame(best2$coef)
se_df <- as.data.frame(best2$se)


coef_df$From <- rownames(coef_df)
se_df$From <- rownames(se_df)

# Reshape to long format
coef_long <- reshape(coef_df, varying = list(names(coef_df)[-ncol(coef_df)]),
                     v.names = "Coefficient", timevar = "To", 
                     times = names(coef_df)[-ncol(coef_df)], direction = "long")

se_long <- reshape(se_df, varying = list(names(se_df)[-ncol(se_df)]),
                   v.names = "StdError", timevar = "To", 
                   times = names(se_df)[-ncol(se_df)], direction = "long")

# Merge the two long data frames
combined_df <- merge(coef_long, se_long, by = c("From", "To"), sort = FALSE)


combined_df <- combined_df[, c("From", "To", "Coefficient", "StdError")]


#write.csv(combined_df, "ModelB_combined.csv", row.names = FALSE)





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

plot_C <- plot_model_set(m3)

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


plot_D <-plot_model_set(m4)

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



mD2 <- define_model_set(
  D2 = c(
    u ~ Ne + MS + GT,
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


plot_E <-plot_model_set(mD2)

pD2<- phylo_path(mD2,path1_d,time_tree_u, lower.bound = 0.000, upper.bound = 1.00)

sD2 <- summary(pD2)
sD2

bestD2 <- best(pD2)
plot(bestD2)

dsep_pD2 <- as.data.frame (pD2$d_sep)
dsep_pD2$D2.model <- NULL
dsep_pD2$significant <- dsep_pD2$D2.p < 0.05
dsep_pD2[] <- lapply(dsep_pD2, function(x) if (is.list(x)) as.character(x) else x)
dsep_pD2
#write.csv(dsep_pD2, "modelD2_dsep_results.csv", row.names = FALSE)


boot_modelD2 <- best(pD2, boot = 1000)
coef_plot(boot_modelD2)

coefD2<- boot_modelD2$coef
lowerD2 <- boot_modelD2$lower
upperD2 <- boot_modelD2$upper

rowsD2 <- rownames(coefD2)
colsD2 <- colnames(coefD2)

# Initialize a data frame to store results
resultsD2 <- data.frame(From = character(), To = character(), Estimate = numeric(), Lower = numeric(), Upper = numeric(), Significant = logical())

# Iterate through the coefficient matrix (excluding the diagonal)
for (i in 1:nrow(coefD2)) {
  for (j in 1:ncol(coefD2)) {
    if (i != j) { # Avoid diagonal elements (self-effects)
      from <- rowsD2[i]
      to <- colsD2[j]
      estimate <- coefD2[i, j]
      lower_bound <- lowerD2[i, j]
      upper_bound <- upperD2[i, j]
      significant <- (lower_bound > 0 || upper_bound < 0) # Check if CI overlaps zero
      
      resultsD2 <- rbind(resultsD2, data.frame(From = from, To = to, Estimate = estimate, Lower = lower_bound, Upper = upper_bound, Significant = significant))
    }
  }
}

print(resultsD2)
#write.csv(resultsD2, file = "modelD2_coefficients.csv", row.names = FALSE)



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
    u ~ Ne + MS + GT,
    Ne ~ GT,
    GT ~ MA + LS,
    LS ~ BM,
    MA ~ BM,
    OPG ~ BM + MA + LS + GT,
    LS ~ MA,
    MS ~ BM,
    MS ~ MA)
  )



fit_all <- phylo_path(full_models, data = path1_d, tree = time_tree_u,lower.bound = 0.000, upper.bound = 1.00)

summary_all <- summary(fit_all)
print(summary_all)

m5 <- define_model_set(
  E = c(
    u ~ MS + GT,
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


plot_F <- plot_model_set(m5)
plot_F
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



full_models2 <- define_model_set(
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
    u ~ Ne + MS + GT,
    Ne ~ GT,
    GT ~ MA + LS,
    LS ~ BM,
    MA ~ BM,
    OPG ~ BM + MA + LS + GT,
    LS ~ MA,
    MS ~ BM,
    MS ~ MA
    ),
   F = c(
    u ~ MS + GT,
    Ne ~ GT,
    GT ~ MA + LS,
    LS ~ BM,
    MA ~ BM,
    OPG ~ BM + MA + LS + GT,
    LS ~ MA,
    MS ~ BM,
    MS ~ MA)
)



fit_all2 <- phylo_path(full_models2, data = path1_d, tree = time_tree_u,lower.bound = 0.000, upper.bound = 1.00)

summary_all2 <- summary(fit_all2)
print(summary_all2)


plot_list <- list(plot_A, plot_B, plot_C, plot_D, plot_E, plot_F)


plot_list <- lapply(plot_list, function(p) {
  p + theme(plot.margin = margin(2, 2, 2, 2))  # Adjust values as needed
})


ggarrange(plotlist = plot_list,
          ncol = 2, nrow = 3,
          labels = c("A", "B", "C", "D", "E", "F"),
          align = "hv")





plot_bestD2 <- plot(boot_modelD2) + 
  ggtitle("(A) Model E Standardized Path Coefficents") +
  theme(plot.title = element_text(size = 12, hjust = 0.5))

plot_bestD2

coef_df <- as.data.frame(boot_modelD2$coef)
lower_df <- as.data.frame(boot_modelD2$lower)
upper_df <- as.data.frame(boot_modelD2$upper)

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
plotD2coef <- ggplot(plot_data, aes(x = Pair, y = Estimate, color = Direction)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2) +
  scale_color_manual(
    values = c("Positive" = "blue", "Negative" = "red", "Non-significant" = "gray50")
  ) +
  theme_minimal() +
  labs(
    title = "(B) Model E Path Coefficients Ordered by Response Variable",
    x = "Life History Trait Pairing",
    y = "Standardized Path Coefficient"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    plot.title = element_text(hjust = 0.5, size = 12),
    legend.position = "none"
  )
plotD2coef

plot_bestE <- plot_bestD2 + 
  theme(
    plot.margin = margin(5, 10, 5, 10),
    plot.title = element_text(size = 11),
    axis.text.x = element_blank(),        # hides tick labels
    axis.ticks.x = element_blank(),       # hides tick marks
    axis.title.x = element_blank(),        # hides axis title
    axis.text.y= element_blank(),
  )
plot_bestD2

fig3 <- ggarrange(
  plot_bestD2, 
  plotD2coef, 
  nrow = 2,
  heights = c(1.5, 1) 
)

fig3


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

mG<- define_model_set(
  G = c(GS ~ Ne + MS + BM + LS + OPG + GT +MA + u)
)

plot_model_set(mG)


pG <- phylo_path(mG,path3_d,time_tree_u, lower.bound = 0.000, upper.bound = 1.00)

sG <- summary(pG)
sG

dsep_pG <- as.data.frame (pG$d_sep)
dsep_pG
dsep_pG$G.model <- NULL
dsep_pG$significant <- dsep_pG$G.p < 0.05
dsep_pG[] <- lapply(dsep_pG, function(x) if (is.list(x)) as.character(x) else x)
#write.csv(dsep_pG, "modelG_dsep_results.csv", row.names = FALSE)

bestG <- best(pG)
plot(bestG)

print(bestG)



#### add life history structure from model E



mH <- define_model_set(
  H = c(
    u ~ MS + GT,
    Ne ~ GT,
    GT ~ MA + LS,
    LS ~ BM,
    MA ~ BM,
    OPG ~ BM + MA + LS + GT,
    LS ~ MA,
    MS ~ BM,
    MS ~ MA,
    GS ~ Ne + MS + BM + LS + OPG + GT +MA + u)
)

plot_model_set(mH)

pH <- phylo_path(mH,path3_d,time_tree_u, lower.bound = 0.000, upper.bound = 1.00)

sH <- summary(pH)
sH

dsep_pH <- as.data.frame (pH$d_sep)
dsep_pH
dsep_pH$H.model <- NULL
dsep_pH$significant <- dsep_pH$H.p < 0.05
dsep_pH[] <- lapply(dsep_pH, function(x) if (is.list(x)) as.character(x) else x)
#write.csv(dsep_pH, "modelH_dsep_results.csv", row.names = FALSE)


bestH <- best(pH)
plot(bestH)



boot_modelH <- best(pH, boot = 1000)


coef_plot(boot_modelH)

coefH<- boot_modelH$coef
lowerH <- boot_modelH$lower
upperH <- boot_modelH$upper

rowsH <- rownames(coefH)
colsH <- colnames(coefH)

# Initialize a data frame to store results
resultsH <- data.frame(From = character(), To = character(), Estimate = numeric(), Lower = numeric(), Upper = numeric(), Significant = logical())

# Iterate through the coefficient matrix (excluding the diagonal)
for (i in 1:nrow(coefH)) {
  for (j in 1:ncol(coefH)) {
    if (i != j) { # Avoid diagonal elements (self-effects)
      from <- rowsH[i]
      to <- colsH[j]
      estimate <- coefH[i, j]
      lower_bound <- lowerH[i, j]
      upper_bound <- upperH[i, j]
      significant <- (lower_bound > 0 || upper_bound < 0) # Check if CI overlaps zero
      
      resultsH <- rbind(resultsH, data.frame(From = from, To = to, Estimate = estimate, Lower = lower_bound, Upper = upper_bound, Significant = significant))
    }
  }
}

print(resultsH)
#write.csv(resultsH, file = "modelH_coefficients.csv", row.names = FALSE)


plot_bestH <- plot(boot_modelH) + 
  ggtitle("(A) Model H Standardized Path Coefficents") +
  theme(plot.title = element_text(size = 12, hjust = 0.5))

plot_bestH


coefH_df <- as.data.frame(boot_modelH$coef)
lowerH_df <- as.data.frame(boot_modelH$lower)
upperH_df <- as.data.frame(boot_modelH$upper)

plot_dataH <- coefH_df %>%
  rownames_to_column("From") %>%
  pivot_longer(-From, names_to = "To", values_to = "Estimate") %>%
  left_join(
    lowerH_df %>% rownames_to_column("From") %>%
      pivot_longer(-From, names_to = "To", values_to = "Lower"),
    by = c("From", "To")
  ) %>%
  left_join(
    upperH_df %>% rownames_to_column("From") %>%
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
response_order <- c("GS","u", "Ne", "GT", "MS", "OPG", "LS", "MA", "BM")

# Order by response variable grouping
plot_dataH <- plot_dataH %>%
  mutate(Response = factor(Response, levels = response_order)) %>%
  arrange(Response, To, From) %>%
  mutate(Pair = factor(Pair, levels = unique(Pair)))

# Create the plot
plotHcoef <- ggplot(plot_dataH, aes(x = Pair, y = Estimate, color = Direction)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2) +
  scale_color_manual(
    values = c("Positive" = "blue", "Negative" = "red", "Non-significant" = "gray50")
  ) +
  theme_minimal() +
  labs(
    title = "(B) Model H Path Coefficients Ordered by Response Variable",
    x = "Life History Trait Pairing",
    y = "Standardized Path Coefficient"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    plot.title = element_text(hjust = 0.5, size = 12),
    legend.position = "none"
  )
plotHcoef

plot_bestH <- plot_bestH + 
  theme(
    plot.title = element_text(size = 11),
    axis.text.x = element_blank(),        # hides tick labels
    axis.ticks.x = element_blank(),       # hides tick marks
    axis.title.x = element_blank(),        # hides axis title
    axis.text.y= element_blank(),
  )

plot_bestH$data[, c("name", "x", "y")]

plot_bestH$data <- plot_bestH$data |>
  dplyr::mutate(
    x = dplyr::case_when(
      name == "OPG" ~ x + 2,
      name == "MS"  ~ x + 1.5,  # move right
      TRUE ~ x
    ),
    y = dplyr::case_when(
      name == "OPG" ~ y + 2.5,
      name == "MS"  ~ y - 0.5,  # move down
      TRUE ~ y
    )
  )

plot_bestH

fig4 <- ggarrange(
  plot_bestH, 
  plotHcoef, 
  nrow = 2,
  heights = c(1.5, 1) 
)

fig4



mI <- define_model_set(
  I = c(
    u ~ MS + GT,
    Ne ~ GT,
    GT ~ MA + LS,
    LS ~ BM,
    MA ~ BM,
    OPG ~ BM + MA + LS + GT,
    LS ~ MA,
    MS ~ BM,
    MS ~ MA,
    GS ~ MS + BM + LS + OPG + GT +MA)
)
plot_model_set(mI)

pI<- phylo_path(mI,path3_d,time_tree_u, lower.bound = 0.000, upper.bound = 1.00)

sI <- summary(pI)
sI

dsep_pI <- as.data.frame (pI$d_sep)
dsep_pI
dsep_pI$I.model <- NULL
dsep_pI$significant <- dsep_pI$I.p < 0.05
dsep_pI[] <- lapply(dsep_pI, function(x) if (is.list(x)) as.character(x) else x)
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



full_models3 <- define_model_set(
  G = c(GS ~ Ne + MS + BM + LS + OPG + GT +MA + u),
  H = c(u ~ MS + GT,
          Ne ~ GT,
          GT ~ MA + LS,
          LS ~ BM,
          MA ~ BM,
          OPG ~ BM + MA + LS + GT,
          LS ~ MA,
          MS ~ BM,
          MS ~ MA,
          GS ~ Ne + MS + BM + LS + OPG + GT +MA + u),   
  I = c(u ~ MS + GT,
    Ne ~ GT,
    GT ~ MA + LS,
    LS ~ BM,
    MA ~ BM,
    OPG ~ BM + MA + LS + GT,
    LS ~ MA,
    MS ~ BM,
    MS ~ MA,
    GS ~ MS + BM + LS + OPG + GT +MA))
  

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




#### redo with Ne estimated from pi



path4_d = data.frame(d$Ne_pi,d$Modeled.rate.per.generation..m_generation_modeled., d$Generation.time..years.,
                     d$Mating.strategy..binary., 
                     d$Mass..grams.,
                     d$Lifespan.in.the.wild..years., 
                     d$Average.maturation.time,
                     d$Number.of.offspring.per.generation, 
                     d$genome_size)
rownames(path4_d) <- d$Species.name
colnames(path4_d) <- c("Ne","u", "GT", "MS", "BM", "LS", "MA","OPG","GS")


species_in_data <- rownames(path4_d)


species_in_data <- rownames(path4_d)
tree_pi_Ne <- drop.tip(time_tree, setdiff(time_tree$tip.label, species_in_data))



name.check(tree_pi_Ne,path4_d)
path4_d<-ReorderData(tree_pi_Ne,path4_d)


mH_pi <- define_model_set(
  H_pi = c(
    u ~ MS + GT,
    Ne ~ GT,
    GT ~ MA + LS,
    LS ~ BM,
    MA ~ BM,
    OPG ~ BM + MA + LS + GT,
    LS ~ MA,
    MS ~ BM,
    MS ~ MA,
    GS ~ Ne + MS + BM + LS + OPG + GT +MA + u)
)

plot_model_set(mH_pi)

pH_pi <- phylo_path(mH_pi,path4_d,tree_pi_Ne, lower.bound = 0.000, upper.bound = 1.00)

sH_pi <- summary(pH_pi)
sH_pi

dsep_pH_pi <- as.data.frame (pH_pi$d_sep)
dsep_pH_pi
dsep_pH_pi$H_pi.model <- NULL
dsep_pH_pi$significant <- dsep_pH_pi$H_pi.p < 0.05
dsep_pH_pi[] <- lapply(dsep_pH_pi, function(x) if (is.list(x)) as.character(x) else x)
#write.csv(dsep_pH_pi, "modelH_pi_dsep_results.csv", row.names = FALSE)


bestH_pi <- best(pH_pi)
plot(bestH_pi)



boot_modelH_pi <- best(pH_pi, boot = 1000)


coef_plot(boot_modelH_pi)

coefH_pi<- boot_modelH_pi$coef
lowerH_pi <- boot_modelH_pi$lower
upperH_pi <- boot_modelH_pi$upper

rowsH_pi <- rownames(coefH_pi)
colsH_pi <- colnames(coefH_pi)

# Initialize a data frame to store results
resultsH_pi <- data.frame(From = character(), To = character(), Estimate = numeric(), Lower = numeric(), Upper = numeric(), Significant = logical())

# Iterate through the coefficient matrix (excluding the diagonal)
for (i in 1:nrow(coefH_pi)) {
  for (j in 1:ncol(coefH_pi)) {
    if (i != j) { # Avoid diagonal elements (self-effects)
      from <- rowsH_pi[i]
      to <- colsH_pi[j]
      estimate <- coefH_pi[i, j]
      lower_bound <- lowerH_pi[i, j]
      upper_bound <- upperH_pi[i, j]
      significant <- (lower_bound > 0 || upper_bound < 0) # Check if CI overlaps zero
      
      resultsH_pi <- rbind(resultsH_pi, data.frame(From = from, To = to, Estimate = estimate, Lower = lower_bound, Upper = upper_bound, Significant = significant))
    }
  }
}

print(resultsH_pi)
#write.csv(resultsH_pi, file = "modelH_pi_coefficients.csv", row.names = FALSE)


plot_bestH_pi <- plot(boot_modelH_pi) + 
  ggtitle("(A) Model H_pi Standardized Path Coefficents") +
  theme(plot.title = element_text(size = 12, hjust = 0.5))

plot_bestH_pi


coefH_pi_df <- as.data.frame(boot_modelH_pi$coef)
lowerH_pi_df <- as.data.frame(boot_modelH_pi$lower)
upperH_pi_df <- as.data.frame(boot_modelH_pi$upper)

plot_dataH_pi <- coefH_pi_df %>%
  rownames_to_column("From") %>%
  pivot_longer(-From, names_to = "To", values_to = "Estimate") %>%
  left_join(
    lowerH_pi_df %>% rownames_to_column("From") %>%
      pivot_longer(-From, names_to = "To", values_to = "Lower"),
    by = c("From", "To")
  ) %>%
  left_join(
    upperH_pi_df %>% rownames_to_column("From") %>%
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


# Order by response variable grouping
plot_dataH_pi <- plot_dataH_pi %>%
  mutate(Response = factor(Response, levels = response_order)) %>%
  arrange(Response, To, From) %>%
  mutate(Pair = factor(Pair, levels = unique(Pair)))

# Create the plot
plotH_picoef <- ggplot(plot_dataH_pi, aes(x = Pair, y = Estimate, color = Direction)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2) +
  scale_color_manual(
    values = c("Positive" = "blue", "Negative" = "red", "Non-significant" = "gray50")
  ) +
  theme_minimal() +
  labs(
    title = "(B) Model H_pi Path Coefficients Ordered by Response Variable",
    x = "Life History Trait Pairing",
    y = "Standardized Path Coefficient"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    plot.title = element_text(hjust = 0.5, size = 12),
    legend.position = "none"
  )
plotH_picoef

plot_bestH_pi <- plot_bestH_pi + 
  theme(
    plot.margin = margin(5, 10, 5, 10),
    plot.title = element_text(size = 11),
    axis.text.x = element_blank(),        # hides tick labels
    axis.ticks.x = element_blank(),       # hides tick marks
    axis.title.x = element_blank(),        # hides axis title
    axis.text.y= element_blank(),
  )




plot_bestH_pi$data <- plot_bestH_pi$data |>
  dplyr::mutate(
    x = dplyr::case_when(
      name == "OPG" ~ x + 2,
      name == "MS"  ~ x + 1.5,  # move right
      TRUE ~ x
    ),
    y = dplyr::case_when(
      name == "OPG" ~ y + 2.5,
      name == "MS"  ~ y - 0.5,  # move down
      TRUE ~ y
    )
  )



plot_bestH_pi










figS3 <- ggarrange(
  plot_bestH_pi, 
  plotH_picoef, 
  nrow = 2,
  heights = c(1.5, 1) 
)

figS3




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


r1 <- ggplot(d1, aes(x = Ne_mu_log, y = genome_size, color = Taxonomic.group)) +
  geom_point(show.legend = TRUE) +
  scale_color_manual(values = my_colors, drop = FALSE) +
  geom_line(aes(y = fit_Ne_mu), color = "black", linetype = "solid") +
  geom_line(aes(y = fit_Ne), color = "blue", linetype = "dashed") +
  geom_line(aes(y = fit_u), color = "red", linetype = "dotted") +
  xlab("log₁₀(Ne × μ)") +
  ylab("Genome size (log₁₀ bp)") +
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "bottom")

r1

d2_vars <- c("Species.name","Taxonomic.group","Modeled.rate.per.generation..m_generation_modeled." ,"genome_size","Ne_pi")
d2 <- d[d2_vars]

d2<- d2 %>% drop_na(Ne_pi)

d2 <-ReorderData(tree_pi_Ne,d2)




vert2 <- comparative.data(phy=tree_pi_Ne, data=d2, names.col = Species.name, vcv = TRUE, na.omit = TRUE, warn.dropped = TRUE)

GS_piNeu_pgls <- pgls( genome_size~ Modeled.rate.per.generation..m_generation_modeled.*Ne_pi, data=vert2, lambda = "ML")
summary(GS_piNeu_pgls)



####### condition on GT



mu_col <- "Modeled.rate.per.generation..m_generation_modeled."
ne_col <- "Ne..harmonic.mean..psmc."
gt_col <- "Generation.time..years."
grp_col <- "Taxonomic.group"

# ---- Data matched to tree ----
vars <- c("Species.name", grp_col, mu_col, ne_col, gt_col)
d_u  <- d[, vars, drop = FALSE]
name.check(time_tree_u, d_u)
d_u <- ReorderData(time_tree_u, d_u)

vert_u <- comparative.data(
  phy         = time_tree_u,
  data        = d_u,
  names.col   = Species.name,
  vcv         = TRUE,
  na.omit     = TRUE,
  warn.dropped = TRUE
)


fml <- as.formula(paste0("`", mu_col, "` ~ `", ne_col, "` + `", gt_col, "`"))
fit_full <- pgls(fml, data = vert_u, lambda = "ML")
s_full   <- summary(fit_full)
print(s_full)

lambda_hat <- as.numeric(s_full$param["lambda"])
cat("Estimated Pagel's λ =", round(lambda_hat, 3), "\n")


if (lambda_hat <= 1e-8) {
  message("λ ≈ 0: refitting GLS without phylogenetic correlation.")
  gls_obj <- nlme::gls(fml, data = vert_u$data, method = "ML")
} else {
  gls_obj <- nlme::gls(
    fml,
    data = vert_u$data,
    method = "ML",
    correlation = ape::corPagel(value = lambda_hat, phy = time_tree_u, fixed = TRUE)
  )
}


resolve_term <- function(pattern, coef_names) {
  hit <- grep(pattern, coef_names, ignore.case = TRUE)
  if (!length(hit)) stop("Term not found for pattern: ", pattern,
                         "\nChoices: ", paste(coef_names, collapse = ", "))
  coef_names[hit[1]]
}
partial_r2_from_t <- function(tval, df_res) (tval^2) / (tval^2 + df_res)

boot_partial_R2 <- function(gls_fit, pgls_fit, term_pattern, B = 10000, seed = 1L) {
  set.seed(seed)
  
  # Coefs / VCOV from GLS (works for class 'gls')
  coefs <- coef(gls_fit)                   # named vector
  V     <- as.matrix(stats::vcov(gls_fit)) # vcov.gls
  
  # Residual df = n - p
  n_used <- nobs(gls_fit)                  # safer than nrow(gls_fit$data)
  p_par  <- length(coefs)
  df_res <- n_used - p_par
  
  # Map term name across PGLS summary rows and GLS coef names
  s_fit    <- summary(pgls_fit)
  row_name <- resolve_term(term_pattern, rownames(s_fit$coef))
  col_name <- resolve_term(term_pattern, names(coefs))
  
  # Point estimate via reported t from PGLS
  t_obs  <- s_fit$coef[row_name, "t value"]
  r2_hat <- partial_r2_from_t(t_obs, df_res)
  
  # Parametric coefficient draws ~ MVN(beta_hat, V)
  draws   <- MASS::mvrnorm(B, mu = coefs, Sigma = V)
  se_term <- sqrt(V[col_name, col_name])
  t_draw  <- draws[, col_name] / se_term
  r2_draw <- partial_r2_from_t(t_draw, df_res)
  ci      <- stats::quantile(r2_draw, c(0.025, 0.975), na.rm = TRUE)
  
  list(term = col_name, r2 = as.numeric(r2_hat), ci = as.numeric(ci), draws = r2_draw)
}


ne_pat <- "Ne\\.\\.harmonic\\.mean\\.\\.psmc\\.|^`?Ne"
gt_pat <- "Generation\\.time\\.\\.years\\.|Generation\\.time"

boot_ne <- boot_partial_R2(gls_obj, fit_full, ne_pat, B = 10000, seed = 42)
boot_gt <- boot_partial_R2(gls_obj, fit_full, gt_pat, B = 10000, seed = 43)

r2_tab <- data.frame(
  term       = c("Ne", "GT"),
  R2_partial = c(boot_ne$r2, boot_gt$r2),
  CI_low     = c(boot_ne$ci[1], boot_gt$ci[1]),
  CI_high    = c(boot_ne$ci[2], boot_gt$ci[2])
)
print(r2_tab, row.names = FALSE, digits = 4)


cf_names <- names(coef(fit_full))
ne_name  <- resolve_term(ne_pat, cf_names)
gt_name  <- resolve_term(gt_pat, cf_names)

beta_Ne <- as.numeric(coef(fit_full)[ne_name])
beta_GT <- as.numeric(coef(fit_full)[gt_name])

dat_used <- vert_u$data
x_Ne  <- dat_used[[ne_name]]
x_GT  <- dat_used[[gt_name]]
grp   <- dat_used[[grp_col]]

e_full  <- residuals(fit_full, type = "response")
pres_Ne <- e_full + beta_Ne * x_Ne
pres_GT <- e_full + beta_GT * x_GT

df_ne_all <- data.frame(Ne = x_Ne, pres = pres_Ne, grp = grp)
df_gt_all <- data.frame(GT = x_GT, pres = pres_GT, grp = grp)

df_ne_lines <- subset(df_ne_all, grp %in% c("Bird","Mammal","Fish"))
df_gt_lines <- subset(df_gt_all, grp %in% c("Bird","Mammal","Fish"))

alltaxa_ne <- data.frame(slope = beta_Ne, intercept = 0, grp = "All taxa")
alltaxa_gt <- data.frame(slope = beta_GT, intercept = 0, grp = "All taxa")

my_colors1 <- c(
  "All taxa" = "#6C3BAA",
  "Bird"     = "#EBCC2A",
  "Reptile"  = "#3B9AB2",
  "Mammal"   = "#35274A",
  "Fish"     = "#F2300F"
)

p_ne <- ggplot(df_ne_all, aes(Ne, pres, color = grp)) +
  geom_point(size = 2, alpha = 0.9) +
  geom_smooth(data = df_ne_lines, method = "lm", se = FALSE,
              linewidth = 0.9, fullrange = TRUE) +
  geom_abline(data = alltaxa_ne,
              aes(slope = slope, intercept = intercept, color = grp),
              linewidth = 1.0, key_glyph = draw_key_path) +
  scale_color_manual(values = my_colors1, drop = FALSE, name = "Group") +
  labs(x = expression(log[10](N[e])),
       y = expression("Partial residuals of " * mu * " (| GT)"),
       title = expression(mu %~% N[e] + GT ~ ":" ~ "effect of " ~ N[e])) +
  theme_classic(base_size = 12) +
  guides(color = guide_legend(override.aes = list(linetype = "solid")))

p_gt <- ggplot(df_gt_all, aes(GT, pres, color = grp)) +
  geom_point(size = 2, alpha = 0.9) +
  geom_smooth(data = df_gt_lines, method = "lm", se = FALSE,
              linewidth = 0.9, fullrange = TRUE) +
  geom_abline(data = alltaxa_gt,
              aes(slope = slope, intercept = intercept, color = grp),
              linewidth = 1.0, key_glyph = draw_key_path) +
  scale_color_manual(values = my_colors1, drop = FALSE, name = "Group") +
  labs(x = expression(log[10](Generation ~ time ~ (years))),
       y = expression("Partial residuals of " * mu * " (| " * N[e] * ")"),
       title = expression(mu %~% N[e] + GT ~ ":" ~ "effect of GT")) +
  theme_classic(base_size = 12) +
  guides(color = guide_legend(override.aes = list(linetype = "solid")))

# ---- Assemble figure with bootstrap R² caption ----
fmt_pct <- function(x) sprintf("%.1f%%", 100 * x)
cap_text <- sprintf(
  "Partial R^2 (bootstrap 95%% CI): Ne = %s [%s–%s]; GT = %s [%s–%s].",
  fmt_pct(r2_tab$R2_partial[r2_tab$term == "Ne"]),
  fmt_pct(r2_tab$CI_low[r2_tab$term == "Ne"]),
  fmt_pct(r2_tab$CI_high[r2_tab$term == "Ne"]),
  fmt_pct(r2_tab$R2_partial[r2_tab$term == "GT"]),
  fmt_pct(r2_tab$CI_low[r2_tab$term == "GT"]),
  fmt_pct(r2_tab$CI_high[r2_tab$term == "GT"])
)

fig_min_adj <- ggarrange(p_ne, p_gt, ncol = 2, common.legend = TRUE, legend = "bottom")
fig_min_adj <- annotate_figure(fig_min_adj, bottom = text_grob(cap_text, size = 10))
print(fig_min_adj)





########## mammals only
d_mammals <- d %>% filter(Taxonomic.group == "Mammal")
mammal_species <- d_mammals$Species.name

# Drop non-mammalian species from the tree
tree_mammals <- drop.tip(time_tree_u, setdiff(time_tree_u$tip.label, mammal_species))

name.check(tree_mammals, d_mammals)
d_mammals <- ReorderData(tree_mammals, d_mammals)

pathM_d = data.frame(d_mammals$Ne..harmonic.mean..psmc.,d_mammals$Modeled.rate.per.generation..m_generation_modeled., d_mammals$Generation.time..years.,
                     d_mammals$Mating.strategy..binary., 
                     d_mammals$Mass..grams.,
                     d_mammals$Lifespan.in.the.wild..years., 
                     d_mammals$Average.maturation.time,
                     d_mammals$Number.of.offspring.per.generation)
rownames(pathM_d) <- d_mammals$Species.name
colnames(pathM_d) <- c("Ne","u", "GT", "MS", "BM", "LS", "MA","OPG")



pM2 <- phylo_path(m2,pathM_d,tree_mammals, lower.bound = 0.000, upper.bound = 1.00)

s2M <- summary(pM2)
s2M

dsep_pM2 <- as.data.frame (pM2$d_sep)
dsep_pM2$B.model <- NULL
dsep_pM2$significant <- dsep_pM2$B.p < 0.05
dsep_pM2[] <- lapply(dsep_pM2, function(x) if (is.list(x)) as.character(x) else x)
dsep_pM2

#write.csv(dsep_pM2, "mammals_modelB_dsep_results.csv", row.names = FALSE)



mE_M <- define_model_set(
  E_M = c(
    u ~ Ne + MS + GT,
    Ne ~ GT,
    GT ~ MA + LS,
    LS ~ BM,
    MA ~ BM,
    OPG ~ BM + MA + LS + GT,
    LS ~ MA,
    MS ~ BM,
    MS ~ MA)
)

plot_model_set(mE_M)

pE_M <- phylo_path(mE_M,pathM_d,tree_mammals, lower.bound = 0.000, upper.bound = 1.00)

sE_M <- summary(pE_M)
sE_M

dsep_pE_M <- as.data.frame (pE_M$d_sep)
dsep_pE_M
dsep_pE_M$E_M.model <- NULL
dsep_pE_M$significant <- dsep_pE_M$E_M.p < 0.05
dsep_pE_M[] <- lapply(dsep_pE_M, function(x) if (is.list(x)) as.character(x) else x)
#write.csv(dsep_pE_M, "modelE_M_dsep_results.csv", row.names = FALSE)


bestE_M <- best(pE_M)
plot(bestE_M)
print(bestE_M)


boot_modelE_M <- best(pE_M, boot = 1000)


coef_plot(boot_modelE_M)

coefE_M<- boot_modelE_M$coef
lowerE_M <- boot_modelE_M$lower
upperE_M <- boot_modelE_M$upper

rowsE_M <- rownames(coefE_M)
colsE_M <- colnames(coefE_M)

# Initialize a data frame to store results
resultsE_M <- data.frame(From = character(), To = character(), Estimate = numeric(), Lower = numeric(), Upper = numeric(), Significant = logical())

# Iterate through the coefficient matrix (excluding the diagonal)
for (i in 1:nrow(coefE_M)) {
  for (j in 1:ncol(coefE_M)) {
    if (i != j) { # Avoid diagonal elements (self-effects)
      from <- rowsE_M[i]
      to <- colsE_M[j]
      estimate <- coefE_M[i, j]
      lower_bound <- lowerE_M[i, j]
      upper_bound <- upperE_M[i, j]
      significant <- (lower_bound > 0 || upper_bound < 0) # Check if CI overlaps zero
      
      resultsE_M <- rbind(resultsE_M, data.frame(From = from, To = to, Estimate = estimate, Lower = lower_bound, Upper = upper_bound, Significant = significant))
    }
  }
}

print(resultsE_M)
#write.csv(resultsE_M, file = "modelE_M_coefficients.csv", row.names = FALSE)


plot_bestE_M <- plot(boot_modelE_M) + 
  ggtitle("(A) Model E Mammals Only Standardized Path Coefficents") +
  theme(plot.title = element_text(size = 12, hjust = 0.5))

plot_bestE_M


coefE_M_df <- as.data.frame(boot_modelE_M$coef)
lowerE_M_df <- as.data.frame(boot_modelE_M$lower)
upperE_M_df <- as.data.frame(boot_modelE_M$upper)

plot_dataE_M <- coefE_M_df %>%
  rownames_to_column("From") %>%
  pivot_longer(-From, names_to = "To", values_to = "Estimate") %>%
  left_join(
    lowerE_M_df %>% rownames_to_column("From") %>%
      pivot_longer(-From, names_to = "To", values_to = "Lower"),
    by = c("From", "To")
  ) %>%
  left_join(
    upperE_M_df %>% rownames_to_column("From") %>%
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


# Order by response variable grouping
plot_dataE_M <- plot_dataE_M %>%
  mutate(Response = factor(Response, levels = response_order)) %>%
  arrange(Response, To, From) %>%
  mutate(Pair = factor(Pair, levels = unique(Pair)))

# Create the plot
plotE_Mcoef <- ggplot(plot_dataE_M, aes(x = Pair, y = Estimate, color = Direction)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2) +
  scale_color_manual(
    values = c("Positive" = "blue", "Negative" = "red", "Non-significant" = "gray50")
  ) +
  theme_minimal() +
  labs(
    title = "(B) Model E Mammals Only Path Coefficients Ordered by Response Variable",
    x = "Life History Trait Pairing",
    y = "Standardized Path Coefficient"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    plot.title = element_text(hjust = 0.5, size = 12),
    legend.position = "none"
  )
plotE_Mcoef

plot_bestE_M <- plot_bestE_M + 
  theme(
    plot.margin = margin(5, 10, 5, 10),
    plot.title = element_text(size = 11),
    axis.text.x = element_blank(),        # hides tick labels
    axis.ticks.x = element_blank(),       # hides tick marks
    axis.title.x = element_blank(),        # hides axis title
    axis.text.y= element_blank(),
  )
plot_bestE_M

figS3 <- ggarrange(
  plot_bestE_M, 
  plotE_Mcoef, 
  nrow = 2,
  heights = c(1.5, 1) 
)

figS3


#### switch GT -> LS for reviewer




mE_rev <- define_model_set(
  E_rev = c(
    u  ~ MS + GT,        
    Ne ~ GT,
    GT ~ MA,
    LS ~ BM + MA + GT,     
    OPG ~ BM + MA + LS + GT,
    MS ~ BM,
    MS ~ MA
  )
)

pE_rev  <- phylo_path(mE_rev, path1_d, time_tree_u, lower.bound = 0.000, upper.bound = 1.00)
sE_rev  <- summary(pE_rev)
bestE_rev <- best(pE_rev)


comp_E <- data.frame(
  model = c("MD2 (LS → GT)", "E_rev (GT → LS)"),
  C     = c(sD2$C,          sE_rev$C),
  p     = c(sD2$p,          sE_rev$p),
  CICc  = c(sD2$CICc,       sE_rev$CICc)
)
comp_E$Delta_CICc <- comp_E$CICc - min(comp_E$CICc)
print(comp_E, row.names = FALSE, digits = 4)












































