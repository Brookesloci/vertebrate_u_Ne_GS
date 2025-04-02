library(ggplot2)
library(dplyr)
library(ape)
library(geiger)
library(evobiR)
library(phytools)
library(ggpmisc)
library(caper)
library(ggpubr)
library(tidyverse)
library(nlme)
library(olsrr)



d <- read.csv("primary_dataset.csv")
time_tree <- read.tree(file="tree_calibrated_UCE_morecal.nwk")
rownames(d)<-d$Species.name
tip <- c("Callithrix_jacchus", "Cercocebus_lunulatus", "Elephas_maximus", "Macaca_mulatta", "Muntiacus_reevesi","Oryzias_latipes")

time_tree <-drop.tip(time_tree,tip)

name.check(time_tree,d)
d <-ReorderData(time_tree,d)

d$Ne_u <- d$Ne..harmonic.mean..psmc. * d$Modeled.rate.per.generation..m_generation_modeled.
d$Nc_u <- d$Nc * d$Modeled.rate.per.generation..m_generation_modeled.
d$expected_pi <- 4* d$Nc_u
d$NeNc <- d$Ne..harmonic.mean..psmc./d$Nc

d$genome_size <- log10(d$genome_size)
d$Ne..harmonic.mean..psmc. <- log10(d$Ne..harmonic.mean..psmc.)
d$Nc <- log10(d$Nc)
d$Modeled.rate.per.generation..m_generation_modeled. <- log10(d$Modeled.rate.per.generation..m_generation_modeled.)
d$Ne_u <- log10(d$Ne_u)
d$UCE.substitution.rate <- log10(d$UCE.substitution.rate)
d$NeNc <- log10(d$NeNc)
d$expected_pi <- log10(d$expected_pi)
d$Mass..grams. <- log10(d$Mass..grams.)
d$Generation.time..years. <- log10(d$Generation.time..years.)
d$Average.maturation.time <- log10(d$Average.maturation.time)
d$Lifespan.in.the.wild..years. <- log10(d$Lifespan.in.the.wild..years.)

my_colors <- c(
  "Bird" = "#EBCC2A",
  "Reptile" = "#3B9AB2",
  "Mammal" = "#35274A",
  "Fish" = "#F2300F"
)

columns_to_plot <- c("Mass..grams.", "Average.maturation.time", "Lifespan.in.the.wild..years.", 
                     "Modeled.rate.per.generation..m_generation_modeled.", 
                     "Ne..harmonic.mean..psmc.","UCE.substitution.rate") 

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

columns_to_plot <- c("Mass..grams.", "Average.maturation.time", "Lifespan.in.the.wild..years.",
                     "Modeled.rate.per.generation..m_generation_modeled.",
                     "Ne..harmonic.mean..psmc.","UCE.substitution.rate", "Generation.time..years.")


# Fit a linear model with Modeled.rate as the response variable
model <- lm(Modeled.rate.per.generation..m_generation_modeled. ~., data = d[, columns_to_plot])


##lambda for individual traits
phylosig(time_tree,d$genome_size,method = "lambda", test= TRUE)
phylosig(time_tree,d$Ne..harmonic.mean..psmc., method = "lambda", test=TRUE)

plot(d$genome_size ~ d$Ne..harmonic.mean..psmc.)

pagel_GS_Ne <- gls(genome_size~Ne..harmonic.mean..psmc., data=d,correlation=corPagel(value=0.8,phy=time_tree,fixed=FALSE))

Pagel0_GS_Ne <- gls(genome_size ~ Ne..harmonic.mean..psmc., correlation = corPagel(value=0, phy=time_tree, fixed = TRUE), data=d )

anova(pagel_GS_Ne,Pagel0_GS_Ne)

summary(pagel_GS_Ne)

GS_Ne.lm <- lm(genome_size~Ne..harmonic.mean..psmc., data = d)

summary(GS_Ne.lm)

d_vars <- c("Species.name","Ne..harmonic.mean..psmc.","Taxonomic.group","genome_size")
d1 <- d[d_vars]

vert1 <- comparative.data(phy=time_tree, data=d1, names.col = Species.name, vcv = TRUE, na.omit = TRUE, warn.dropped = TRUE)

GS_Ne_pgls <- pgls(genome_size ~ Ne..harmonic.mean..psmc., data=vert1, lambda = "ML")
summary(GS_Ne_pgls)

coef <- coefficients(GS_Ne_pgls)
##Intercept
int1 <-coef[1] 
## predictor 1
B1 <- coef[2]


r1 <- ggplot(data=d, aes(x=Ne..harmonic.mean..psmc., y=genome_size)) + 
  xlab("Harmonic mean Ne (log10)") + ylab("genome size (log10 bp)") +
  geom_point()  + geom_abline(intercept=int1, slope=B1)+ geom_point(aes(color=Taxonomic.group)) + 
  scale_color_manual(values = my_colors)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                                        panel.background = element_blank(), axis.line = element_line(colour = "black") + theme(legend.position = "bottom") )
r1


##census pop ~ Ne

d2_vars <- c("Species.name","Ne..harmonic.mean..psmc.","Taxonomic.group","Nc")
d2 <- d[d2_vars]

d2 <- d2 %>% drop_na(Nc)


tip1 <- c("Ailurus_fulgens","Amphiprion_ocellaris","Aptenodytes_forsteri","Betta_splendens","Canis_lupus_familiaris","Chauna_torquata","Chrysemys_picta","Coleonyx_brevis","Eublepharis_macularius","Felis_catus","Fukomys_damarensis","Larus_marinus","Moschus_berezovskii","Neovison_vison","Orcinus_orca","Phoenicopterus_roseus","Platalea_ajaja","Pogona_vitticeps","Rhea_pennata","Rousettus_aegyptiacus","Saxicola_maurus","Sphaerodactylus_inigoi","Syngnathus_scovelli","Thamnophis_sirtalis","Tursiops_truncatus","Vicugna_pacos") 
time_tree_Nc <-drop.tip(time_tree,tip1)

name.check(time_tree_Nc,d2)
d2 <-ReorderData(time_tree_Nc,d2)


plot(d2$Ne..harmonic.mean..psmc.~d2$Nc)
phylosig(time_tree_Nc,d2$Nc, method = "lambda", test=TRUE)


pagel_Nc_Ne <- gls(Ne..harmonic.mean..psmc.~Nc, data=d2,correlation=corPagel(value=0.8,phy=time_tree_Nc,fixed=FALSE))

intervals(pagel_Nc_Ne,which="var-cov")

Pagel0_Nc_Ne <- gls(Ne..harmonic.mean..psmc. ~ Nc, correlation = corPagel(value=0, phy=time_tree_Nc, fixed = TRUE), data=d2 )

anova(pagel_Nc_Ne,Pagel0_Nc_Ne)

summary(pagel_Nc_Ne)

Nc_Ne.lm <- lm(Ne..harmonic.mean..psmc.~Nc, data = d2)

summary(Nc_Ne.lm)

vert2 <- comparative.data(phy=time_tree_Nc, data=d2, names.col = Species.name, vcv = TRUE, na.omit = TRUE, warn.dropped = TRUE)

Nc_Ne_pgls <- pgls(Ne..harmonic.mean..psmc.~Nc, data=vert2, lambda = "ML")
summary(Nc_Ne_pgls)

coef1 <- coefficients(Nc_Ne_pgls)
##Intercept
int2 <-coef1[1] 
## predictor 1
B2 <- coef1[2]


r2 <- ggplot(data=d2, aes(x=Nc, y=Ne..harmonic.mean..psmc.)) + 
  ylab("Harmonic mean Ne (log10)") + xlab("Nc (log10") +   geom_smooth(method="lm", se=FALSE) +
  scale_color_manual(values = my_colors) +
  geom_point() + geom_abline(intercept=int2, slope=B2)+ geom_point(aes(color=Taxonomic.group)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                                       panel.background = element_blank(), axis.line = element_line(colour = "black") + theme(legend.position = "bottom") )
r2

###Nc vs GS

d3_vars <- c("Species.name","genome_size","Taxonomic.group","Nc")
d3 <- d[d3_vars]

d3 <- d3 %>% drop_na(Nc)
name.check(time_tree_Nc,d3)
d3 <-ReorderData(time_tree_Nc,d3)


plot(d3$genome_size~d3$Nc)

Nc_GS.lm <- lm(genome_size~Nc, data = d3)

summary(Nc_GS.lm)

vert3 <- comparative.data(phy=time_tree_Nc, data=d3, names.col = Species.name, vcv = TRUE, na.omit = TRUE, warn.dropped = TRUE)

Nc_GS_pgls <- pgls(genome_size~Nc, data=vert3, lambda = "ML")
summary(Nc_GS_pgls)

pagel_Nc_GS <- gls(genome_size~Nc, data=d3,correlation=corPagel(value=0.8,phy=time_tree_Nc,fixed=FALSE))

intervals(pagel_Nc_GS,which="var-cov")

Pagel0_Nc_GS <- gls(genome_size ~ Nc, correlation = corPagel(value=0, phy=time_tree_Nc, fixed = TRUE), data=d3 )

anova(pagel_Nc_GS,Pagel0_Nc_GS)

summary(pagel_Nc_GS)

coef2 <- coefficients(Nc_GS_pgls)
##Intercept
int3 <-coef2[1] 
## predictor 1
B3 <- coef2[2]


r3 <- ggplot(data=d3, aes(x=Nc, y=genome_size)) + 
  ylab("genome size (log10 bp") + xlab("Nc (log10)")  +
  scale_color_manual(values = my_colors)+
  geom_point() + geom_abline(intercept=int3, slope=B3)+ geom_point(aes(color=Taxonomic.group)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                                       panel.background = element_blank(), axis.line = element_line(colour = "black") + theme(legend.position = "bottom") )
r3


### Ne u business



d4_vars <- c("Species.name","Ne..harmonic.mean..psmc.","Taxonomic.group","Modeled.rate.per.generation..m_generation_modeled.")
d4 <- d[d4_vars]

d4 <- d4 %>% drop_na(Modeled.rate.per.generation..m_generation_modeled.)

tip4 <- c("Arctocephalus_gazella","Betta_splendens","Chrysemys_picta","Fukomys_damarensis","Larus_argentatus","Larus_marinus","Pogona_vitticeps","Rhea_pennata","Rousettus_aegyptiacus","Sphaerodactylus_inigoi","Syngnathus_scovelli","Thamnophis_sirtalis","Tupaia_belangeri" )
time_tree_u <-drop.tip(time_tree,tip4)

name.check(time_tree_u,d4)

d4 <-ReorderData(time_tree_u,d4)

phylosig(time_tree_u,d4$Modeled.rate.per.generation..m_generation_modeled., method = "lambda", test=TRUE)


plot(d4$Modeled.rate.per.generation..m_generation_modeled. ~ d4$Ne..harmonic.mean..psmc.)

pagel_u_Ne <- gls(Modeled.rate.per.generation..m_generation_modeled. ~Ne..harmonic.mean..psmc., data=d4,correlation=corPagel(value=0.8,phy=time_tree_u,fixed=FALSE))

intervals(pagel_u_Ne,which="var-cov")

Pagel0_u_Ne <- gls(Modeled.rate.per.generation..m_generation_modeled. ~ Ne..harmonic.mean..psmc., correlation = corPagel(value=0, phy=time_tree_u, fixed = TRUE), data=d4 )

anova(pagel_u_Ne,Pagel0_u_Ne)

summary(pagel_u_Ne)

u_Ne.lm <- lm(Modeled.rate.per.generation..m_generation_modeled.~Ne..harmonic.mean..psmc., data = d4)

summary(u_Ne.lm)

vert4 <- comparative.data(phy=time_tree_u, data=d4, names.col = Species.name, vcv = TRUE, na.omit = TRUE, warn.dropped = TRUE)

u_Ne_pgls <- pgls( Modeled.rate.per.generation..m_generation_modeled.~ Ne..harmonic.mean..psmc., data=vert4)
summary(u_Ne_pgls)

coef3 <- coefficients(u_Ne_pgls)
##Intercept
int4 <-coef3[1] 
## predictor 1
B4 <- coef3[2]


r4 <- ggplot(data=d4, aes(x=Ne..harmonic.mean..psmc., y=Modeled.rate.per.generation..m_generation_modeled.)) + 
  xlab("Harmonic mean Ne (log10)") + ylab("Mutation rate per site per generation (10-8) (log10)") +
  scale_color_manual(values = my_colors)+
  geom_point()  + geom_abline(intercept=int4, slope=B4)+ geom_point(aes(color=Taxonomic.group)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                                                                              panel.background = element_blank(), axis.line = element_line(colour = "black") + theme(legend.position = "bottom") )
r4


##u vs GS



d6_vars <- c("Species.name","genome_size","Taxonomic.group","Modeled.rate.per.generation..m_generation_modeled.")
d6 <- d[d6_vars]

d6 <- d6 %>% drop_na(Modeled.rate.per.generation..m_generation_modeled.)

name.check(time_tree_u,d6)
d6 <-ReorderData(time_tree_u,d6)

plot(d6$genome_size ~d6$Modeled.rate.per.generation..m_generation_modeled. )

pagel_u_GS <- gls(genome_size ~Modeled.rate.per.generation..m_generation_modeled., data=d6,correlation=corPagel(value=0.8,phy=time_tree_u,fixed=FALSE))

intervals(pagel_u_GS,which="var-cov")

Pagel0_u_GS <- gls(genome_size ~Modeled.rate.per.generation..m_generation_modeled. , correlation = corPagel(value=0, phy=time_tree_u, fixed = TRUE), data=d6 )

anova(pagel_u_GS,Pagel0_u_GS)

summary(pagel_u_GS)

u_GS.lm <- lm(genome_size ~ Modeled.rate.per.generation..m_generation_modeled., data = d6)

summary(u_GS.lm)

vert6 <- comparative.data(phy=time_tree_u, data=d6, names.col = Species.name, vcv = TRUE, na.omit = TRUE, warn.dropped = TRUE)

u_GS_pgls <- pgls( genome_size ~ Modeled.rate.per.generation..m_generation_modeled., data=vert6)
summary(u_GS_pgls)

coef5 <- coefficients(u_GS_pgls)
##Intercept
int6 <-coef5[1] 
## predictor 1
B6 <- coef5[2]


r4a <- ggplot(data=d6, aes(x=Modeled.rate.per.generation..m_generation_modeled., y=genome_size)) + 
  xlab("Mutation rate per site per generation (10-8) (log10)") + ylab("genome size (log10 bp)") +
  scale_color_manual(values = my_colors)+
  geom_point()  + geom_abline(intercept=int6, slope=B6)+ geom_point(aes(color=Taxonomic.group)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                                        panel.background = element_blank(), axis.line = element_line(colour = "black") + theme(legend.position = "bottom") )
r4a


d5_vars <- c("Species.name","Taxonomic.group","Ne_u" ,"genome_size")
d5 <- d[d5_vars]

d5 <- d5 %>% drop_na(Ne_u)

name.check(time_tree_u,d5)

d5 <-ReorderData(time_tree_u,d5)
plot(d5$genome_size~d5$Ne_u)


phylosig(time_tree_u,d5$Ne_u, method = "lambda", test=TRUE)

pagel_GS_u_Ne <- gls(genome_size ~Ne_u, data=d5,correlation=corPagel(value=0.8,phy=time_tree_u,fixed=FALSE))

Pagel0_GS_u_Ne <- gls(genome_size~ Ne_u, correlation = corPagel(value=0, phy=time_tree_u, fixed = TRUE), data=d5 )

anova(pagel_GS_u_Ne,Pagel0_GS_u_Ne)

summary(pagel_GS_u_Ne)

GS_u_Ne.lm <- lm(genome_size~Ne_u, data = d5)

summary(GS_u_Ne.lm)

vert5 <- comparative.data(phy=time_tree_u, data=d5, names.col = Species.name, vcv = TRUE, na.omit = TRUE, warn.dropped = TRUE)

GS_u_Ne_pgls <- pgls( genome_size~ Ne_u, data=vert5)
summary(GS_u_Ne_pgls)

coef4 <- coefficients(GS_u_Ne_pgls)
##Intercept
int5 <-coef4[1] 
## predictor 1
B5 <- coef4[2]


r5 <- ggplot(data=d5, aes(x=Ne_u, y=genome_size)) + 
  xlab("Ne * u (log10)") + ylab("genome size (log10 bp)")+
  scale_color_manual(values = my_colors) +
  geom_point()  + geom_abline(intercept=int5, slope=B5)+ geom_point(aes(color=Taxonomic.group)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                                        panel.background = element_blank(), axis.line = element_line(colour = "black") + theme(legend.position = "bottom") )
r5

figure2 <- ggarrange(r1, r2, r3,r4a, r5,
                     labels = c("A", "B", "C","D","E"),
                     common.legend = TRUE, legend = "bottom",
                     ncol = 2, nrow = 3)
figure2




##show LP
d7_vars <- c("Species.name","Taxonomic.group","UCE.substitution.rate" ,"Nc", "expected_pi")
d7 <- d[d7_vars]

d7 <- d7 %>% drop_na(Nc)
d7 <- d7 %>% drop_na(expected_pi)
tip5 <- c("Arctocephalus_gazella","Larus_argentatus", "Tupaia_belangeri"  )

time_tree_Nc_u <-drop.tip(time_tree_Nc,tip5)


name.check(time_tree_Nc_u,d7)

d7 <-ReorderData(time_tree_Nc_u,d7)

pagel_pi_Nc <- gls(UCE.substitution.rate ~Nc, data=d7,correlation=corPagel(value=0.8,phy=time_tree_Nc_u,fixed=FALSE))

Pagel0_pi_Nc <- gls(UCE.substitution.rate ~Nc, correlation = corPagel(value=0, phy=time_tree_Nc_u, fixed = TRUE), data=d7)

anova(pagel_pi_Nc,Pagel0_pi_Nc)

summary(pagel_pi_Nc)

pi_Nc.lm <- lm(UCE.substitution.rate ~Nc, data = d7)

summary(pi_Nc.lm)

vert7 <- comparative.data(phy=time_tree_Nc_u, data=d7, names.col = Species.name, vcv = TRUE, na.omit = TRUE, warn.dropped = TRUE)

pi_Nc_pgls <- pgls(UCE.substitution.rate ~Nc, data=vert7)
summary(pi_Nc_pgls)

coef6 <- coefficients(pi_Nc_pgls)
##Intercept
int7 <-coef6[1] 
## predictor 1
B7 <- coef6[2]

r6 <- ggplot(data = d7, aes(x = Nc, y = UCE.substitution.rate)) +
  xlab("Census population size (log10)") + 
  ylab("UCE substitution rate (log10)") +
  geom_point() +
  scale_color_manual(values = my_colors) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_abline(intercept = int7, slope = B7) +
  geom_line(data = d7, aes(x = Nc, y = expected_pi), color = "red") +
  geom_point(aes(color = Taxonomic.group)) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), 
        legend.position = "none")

r6



###Ne/Nc ~ u


d8_vars <- c("Species.name","Taxonomic.group","NeNc","Modeled.rate.per.generation..m_generation_modeled.")

d8 <- d[d8_vars]

d8 <- d8 %>% drop_na(NeNc)
d8 <- d8 %>% drop_na(Modeled.rate.per.generation..m_generation_modeled.)





d8 <-ReorderData(time_tree_Nc_u,d8)


phylosig(time_tree_Nc_u,d8$NeNc, method = "lambda", test=TRUE)

pagel_NeNc_u <- gls(NeNc ~Modeled.rate.per.generation..m_generation_modeled., data=d8,correlation=corPagel(value=0.8,phy=time_tree_Nc_u,fixed=FALSE))

Pagel0_NeNc_u <- gls(NeNc~ Modeled.rate.per.generation..m_generation_modeled., correlation = corPagel(value=0, phy=time_tree_Nc_u, fixed = TRUE), data=d8)

anova(pagel_NeNc_u,Pagel0_NeNc_u)

summary(pagel_NeNc_u)

NeNc_u.lm <- lm(NeNc~Modeled.rate.per.generation..m_generation_modeled., data = d8)

summary(NeNc_u.lm)

vert8 <- comparative.data(phy=time_tree_Nc_u, data=d8, names.col = Species.name, vcv = TRUE, na.omit = TRUE, warn.dropped = TRUE)

NeNc_u_pgls <- pgls( NeNc~ Modeled.rate.per.generation..m_generation_modeled., data=vert8)
summary(NeNc_u_pgls)

coef7 <- coefficients(NeNc_u_pgls)
##Intercept
int8 <-coef7[1] 
## predictor 1
B8 <- coef7[2]


r7 <- ggplot(data=d8, aes(x=Modeled.rate.per.generation..m_generation_modeled., y=NeNc)) + 
  xlab(" Mutation rate (log10)") + ylab("(Ne/Nc) (log10)") + geom_abline(intercept=int8, slope=B8) +
  scale_color_manual(values = my_colors)+
  geom_point()  + geom_smooth(method="lm", se=FALSE) + geom_point(aes(color=Taxonomic.group)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                                      panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position = "none") 
r7



###Ne/Nc ~ Nc


d9_vars <- c("Species.name","Taxonomic.group","NeNc","Nc")

d9 <- d[d9_vars]

d9 <- d9 %>% drop_na(NeNc)
name.check(time_tree_Nc,d9)


d9 <-ReorderData(time_tree_Nc,d9)

plot(d9$NeNc~d9$Nc)
pagel_NeNc_Nc <- gls(NeNc ~Nc, data=d9,correlation=corPagel(value=0.8,phy=time_tree_Nc,fixed=FALSE))

Pagel0_NeNc_Nc <- gls(NeNc~ Nc, correlation = corPagel(value=0, phy=time_tree_Nc, fixed = TRUE), data=d9)

anova(pagel_NeNc_Nc,Pagel0_NeNc_Nc)

summary(pagel_NeNc_Nc)

NeNc_Nc.lm <- lm(NeNc~Nc, data = d9)

summary(NeNc_Nc.lm)

vert9 <- comparative.data(phy=time_tree_Nc, data=d9, names.col = Species.name, vcv = TRUE, na.omit = TRUE, warn.dropped = TRUE)

NeNc_Nc_pgls <- pgls( NeNc~ Nc, data=vert9)
summary(NeNc_Nc_pgls)

coef8 <- coefficients(NeNc_Nc_pgls)
##Intercept
int9 <-coef8[1] 
## predictor 1
B9 <- coef8[2]


r8 <- ggplot(data=d9, aes(x=Nc, y=NeNc)) +xlab(" Nc (log10)") + ylab("(Ne/Nc) (log10)") + geom_abline(intercept=int9, slope=B9)+
  geom_point()  +
  scale_color_manual(values = my_colors)+ geom_smooth(method="lm", se=FALSE) + geom_point(aes(color=Taxonomic.group)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                                      panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position = "none") 
r8

figure3 <- ggarrange(r6, r8, r7,
                     labels = c("A", "B", "C"),
                     common.legend = FALSE,
                     ncol = 2, nrow = 2)
figure3



























