#### 0 - Clear the memory and import libraries ####
rm(list=ls())
library(data.table)
library(tidyverse)
library(zoo)
library(plotrix)
library(sf)
library(rnaturalearth)
library(rworldmap)
library(ggrepel)
library(cowplot)
library(ggpubr)
library(Hmisc)
library(corrplot)
library(FSA)

#### 1 - Import data ####
  #### Import Site information ####

master <- fread("data_table.csv")
cases_all <- fread("cases.csv")

  #### Record and use LoD and LoQ for later use ####
# N1,  E gene and CrAssphage LoQs/ LoDs are expressed in copies per 100ml of sewage
# LoDs/ LoQs assume standard curves used 1x 10^x, rather than 5x 10^x as was the case
# LoDs/ LoQs for N1 were 2x higher for the first week of the study. Only one point in
# week 1 is between LoQ and LoD and would be anyway if the lower LoQ and LoD thresholds
# are used that assume 100 ml sample volume

N1_LoQ <- 59*20 
N1_LoD <- 9*20
E_LoQ <- 126*20
E_LoD <- 22 * 20
C_LoQ <- 20 * 20
C_LoD <- 2 * 20

#### Handling of samples between LoD and LoQ for correlation analysis ####
# - Create a data table for correlation analysis
corr_dt <- master
# - Set N1 levels below LoD to zero
corr_dt$gc_100ml[corr_dt$gc_100ml < N1_LoD] <- 0
corr_dt$gc_100ml[corr_dt$gc_100ml >= N1_LoD & corr_dt$gc_100ml < N1_LoQ] <- 0.1
corr_dt$E_gene_100ml[corr_dt$E_gene_100ml < E_LoD] <- 0
corr_dt$E_gene_100ml[corr_dt$E_gene_100mll >= E_LoD & corr_dt$E_gene_100ml < E_LoQ] <- 0.1

#### Chemistry data ####

# Flow
shapiro.test((corr_dt$Flow/ corr_dt$Population_Equivalent))
# p <2.2e-16
kruskal.test((corr_dt$Flow/ corr_dt$Population_Equivalent), corr_dt$Site_Name)
# p <2.2e-16
dunnTest((corr_dt$Flow/ corr_dt$Population_Equivalent) ~  corr_dt$Site_Name)
# Gwynedd-Manchester not sig, all others highly sig
cor.test((corr_dt$Flow/ corr_dt$Population_Equivalent), corr_dt$gc_100ml, method = "spearman", use = "complete.obs")
# p = 0.1108

# NH4
shapiro.test(corr_dt$NH4)
# p = 5.986e-07
kruskal.test(corr_dt$NH4, corr_dt$Site_Name)
# p = 4.402e-05
cor.test(corr_dt$NH4, corr_dt$gc_100ml, method = "spearman", use = "complete.obs")
# p = 0.8238

# PO4
shapiro.test(corr_dt$PO4)
# p = 5.114e-05
kruskal.test(corr_dt$PO4, corr_dt$Site_Name)
# p = 0.0006494
dunnTest(corr_dt$PO4 ~ corr_dt$Site_Name)
# Cardiff - Liverpool (0.006), Cardiff - Wirral (0.02), Liverpool-Wrexham (0.017), Wirral-Wrexham (0.0490)
cor.test(corr_dt$PO4, corr_dt$gc_100ml, method = "spearman", use = "complete.obs")
# p = 0.1462

# pH
shapiro.test(corr_dt$pH)
# p = 0.009681
kruskal.test(corr_dt$pH, corr_dt$Site_Name)
# p = 0.004938
dunnTest(corr_dt$pH ~ corr_dt$Site_Name)
# Gwynedd-Wirral (0.0300)
cor.test(corr_dt$pH, corr_dt$gc_100ml, method = "spearman", use = "complete.obs")
# p = 0.8141

# EC
shapiro.test(corr_dt$EC)
# p = 1.549e-07
kruskal.test(corr_dt$EC, corr_dt$Site_Name)
# p-value = 7.178e-08
cor.test(corr_dt$EC, corr_dt$gc_100ml, method = "spearman", use = "complete.obs")
# p = 0.5206

# NO3
shapiro.test(corr_dt$NO3)
# p = 2.678e-15
kruskal.test(corr_dt$NO3, corr_dt$Site_Name)
# p = 0.003202
cor.test(corr_dt$NO3, corr_dt$gc_100ml, method = "spearman", use = "complete.obs")
# p = 0.06433

#### Extraction Efficiency Analysis ####
summary(corr_dt$Crass_Recovery)
# mean = 25.9%
std.error(corr_dt$Crass_Recovery)
# SEM = 4.27
colSums(!is.na(corr_dt))
# n = 20

shapiro.test(corr_dt$Crass_Recovery)
# p = 1.69e-5
kruskal.test(corr_dt$Crass_Recovery, corr_dt$Site_Name)
# p = 0.07278
kruskal.test(corr_dt$Crass_Recovery, corr_dt$Week)
# p = 0.6529

cor.test(corr_dt$Crass_Recovery, corr_dt$gc_100ml, method = "spearman", use = "complete.obs")
p = 0.5077

#### Positive tests descriptive stats ####
# Classify N1 quantities
master$cov_class[is.na(master$gc_100ml) == F] <- "Quantified"
master$cov_class[master$gc_100ml < N1_LoQ] <- "Detected"
master$cov_class[master$gc_100ml < N1_LoD] <- "Low amplification"
master$cov_class[master$gc_100ml == 0] <- "No amplification"

# Add label if CoV was detected
master$cov_detect <- T
master$cov_detect[master$gc_100ml < N1_LoD] <- F

# Create summary table of tests
pos_tests <- master[is.na(master$cov_class) == F]%>%
  group_by(Site_Name) %>% dplyr::count(cov_class)

pos_tests <- pos_tests[pos_tests$cov_class != "",]

# Calculate percentages
pos_tests$percent <- (pos_tests$n/15) * 100

# Order values
pos_tests$cov_class <- ordered(pos_tests$cov_class, levels = c("Quantified", "Detected", "Low amplification", "No amplification"))

# Plot figure
ggplot(pos_tests) + theme_bw() + geom_bar(aes(x = Site_Name, y = percent, fill = cov_class),
                        position = "stack", stat = "identity") +
  theme(legend.position = "bottom") +
  labs(x = "WWTP", y = "Proportion of test results", fill = element_blank()) +
  theme(text = element_text(size = 15)) + 
  scale_fill_manual(values = c("red3", "orange1", "grey50","grey"))
sum(pos_tests$n)
# Calculate the mean and sem of positive tests
Detected <- pos_tests[pos_tests$cov_class %in% c("Quantified", "Detected"),]
DetectedAgregate <- aggregate(Detected$percent, by = list(Site = Detected$Site_Name), sum)
summary(DetectedAgregate)
# mean = 64%
std.error(DetectedAgregate$x)
# 6.8 sem

# Calculate the mean and sem of quantified tests
Quantified <- pos_tests[pos_tests$cov_class %in% c("Quantified"),]
QuantifiedAgregate <- aggregate(Quantified$percent, by = list(Site = Quantified$Site_Name), sum)
summary(QuantifiedAgregate)
# mean = 28.9%
std.error(QuantifiedAgregate$x)
# 2.2 sem

# Calculate the mean and sem of the minimum cases for a positive test
Detected_all <- master[master$cov_class %in% c("Quantified", "Detected")]
Detected_all <- aggregate(Detected_all$Daily_Cases_100k, by = list(Site = Detected_all$Site_Name), min)
summary(Detected_all)
# mean minimum cases = 1.2
std.error(Detected_all$x)
# SEM = 0.26

#### E gene data analysis ####

shapiro.test((corr_dt$E_gene_100ml))
# p 1.811e-14
kruskal.test((corr_dt$E_gene_100ml), corr_dt$Site_Name)
# p = 0.1414
cor.test(corr_dt$E_gene_100ml, corr_dt$gc_100ml, method = "spearman", use = "complete.obs")
# rho = 0.538, p = 3.969e-5


E_LoQ/ N1_LoQ
# 2.14
E_LoD/ N1_LoD
#2.44

# Classify E gene values
master$E_class[is.na(master$E_gene_100ml) == F] <- "Quantified"
master$E_class[master$E_gene_100ml < E_LoQ] <- "Detected"
master$E_class[master$E_gene_100ml < E_LoD] <- "Low amplification"
master$E_class[master$E_gene_100ml == 0] <- "No amplification"

# Create a separate table of tests that contained both E and N1
Detected_markers <- master[master$cov_class %in% c("Quantified", "Detected", "Low amplification", "No amplification")]
Detected_markers <- Detected_markers[Detected_markers$E_class %in% c("Quantified", "Detected", "Low amplification", "No amplification")]

# Count N1 classes
Detected_N1 <- Detected_markers %>%
  dplyr::count(cov_class)

# Label as N1
Detected_N1$assay <- "N1"

# Count E classes
Detected_E <- Detected_markers %>%
  dplyr::count(E_class)

# Label as E and rename to allow merging
Detected_E$assay <- "E"
Detected_E <- Detected_E %>% dplyr::rename(cov_class = E_class)

# Merge anc calculate percentages
Detected_marker_summary <- full_join(Detected_N1, Detected_E)
Detected_marker_summary$percentage <- NA
sumN <- sum(Detected_marker_summary$n[Detected_marker_summary$assay == "N1"])
sumE <- sum(Detected_marker_summary$n[Detected_marker_summary$assay == "E"])
Detected_marker_summary$percentage[Detected_marker_summary$assay == "N1"] <- (Detected_marker_summary$n[Detected_marker_summary$assay == "N1"]/sumN)*100
Detected_marker_summary$percentage[Detected_marker_summary$assay == "E"] <- (Detected_marker_summary$n[Detected_marker_summary$assay == "E"]/sumE)*100

# Set order
Detected_marker_summary$cov_class <- ordered(Detected_marker_summary$cov_class, levels = c("Quantified", "Detected", "Low amplification", "No amplification"))

# Plot figure
ggplot(Detected_marker_summary) + theme_bw() + geom_bar(aes(x = assay, y = percentage, fill = cov_class),
                                     position = "stack", stat = "identity", width = 0.5) +
  theme(legend.position = "bottom") +
  labs(x = "WWTP", y = "Proportion of test results", fill = element_blank()) +
  theme(text = element_text(size = 15)) + 
  scale_x_discrete(labels = c("E Sarbeco", "N1 CDC")) +
  scale_fill_manual(values = c("red3", "orange1", "grey50","grey"))

#### Data plotting ####
  #### Figure 2 - Plot Chem data ####

# Create data table of just chemistry and CoV data
corr_dt <- master[is.na(master$gc_100ml) == F]
corr_dt <- select(corr_dt, -c(Crass_Recovery, E_gene_100ml))
corr_dt <- unique(corr_dt)

textsize <- 14

# Plot flow sub-plot
Flow_plt <- ggplot(corr_dt) +
  geom_boxplot(aes(x = Site_Name, y = Flow/ Population_Equivalent, fill = Site_Name),
               outlier.shape = NA) +
  geom_jitter(aes(x = Site_Name, y = Flow/ Population_Equivalent), width = 0.1) +
  theme_bw() +
  scale_fill_brewer(palette = "Dark2") +
  labs(x = element_blank(), y = expression("Daily flow/PE"~(m^{3}~PE^{-1})), fill = "WWTP") +
  scale_x_discrete(guide = guide_axis(n.dodge=3)) +
  ylim(min = 0, max = 0.525) +
  theme(text = element_text(size = textsize)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

# Plot NH4 sub-plot
NH4_plt <- ggplot(corr_dt) +
  geom_boxplot(aes(x = Site_Name, y = NH4, fill = Site_Name),
               outlier.shape = NA) +
  geom_jitter(aes(x = Site_Name, y = NH4), width = 0.1) +
  theme_bw() +
  scale_fill_brewer(palette = "Dark2") +
  labs(x = element_blank(), y = expression(NH[4]^{"+"}~"(mg"~L^{"-1"}~")"), fill = "WWTP") +
  scale_x_discrete(guide = guide_axis(n.dodge=3)) +
  ylim(min = 0, max = 90) +
  theme(text = element_text(size = textsize)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

# Plot NO3 sub-plot (not used)
NO3_plt <- ggplot(corr_dt) +
  geom_boxplot(aes(x = Site_Name, y = NO3, fill = Site_Name),
               outlier.shape = NA) +
  geom_jitter(aes(x = Site_Name, y = NO3), width = 0.1) +
  theme_bw() +
  scale_fill_brewer(palette = "Dark2") +
  labs(x = element_blank(), y = expression(NO[3]^{"-"}~concentration~"(mg"~L^{"-1"}~")"), fill = "WWTP") +
  scale_x_discrete(guide = guide_axis(n.dodge=3)) +
  ylim(min = 0, max = 90) +
  theme(text = element_text(size = textsize)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

# Plot PO4 sub-plot
PO4_plt <- ggplot(corr_dt) +
  geom_boxplot(aes(x = Site_Name, y = PO4, fill = Site_Name),
               outlier.shape = NA) +
  geom_jitter(aes(x = Site_Name, y = PO4), width = 0.1) +
  theme_bw() +
  scale_fill_brewer(palette = "Dark2") +
  labs(x = element_blank(), y = expression(MRP~"(mg"~L^{"-1"}~")"), fill = "WWTP") +
  scale_x_discrete(guide = guide_axis(n.dodge=3)) +
  ylim(min = 0, max = 8.5) +
  theme(text = element_text(size = textsize)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

# Plot EC sub-plot
EC_plt <- ggplot(corr_dt) +
  geom_boxplot(aes(x = Site_Name, y = EC, fill = Site_Name),
               outlier.shape = NA) +
  geom_jitter(aes(x = Site_Name, y = EC), width = 0.1) +
  theme_bw() +
  scale_fill_brewer(palette = "Dark2") +
  labs(x = element_blank(), y = expression(paste("EC (",mu,"S cm"^{-1}~")")),
       fill = "WWTP") +
  scale_x_discrete(guide = guide_axis(n.dodge=3)) +
  ylim(min = 0, max = 3500) +
  theme(text = element_text(size = textsize)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

# Plot pH sub-plot
pH_plt <- ggplot(corr_dt) +
  geom_boxplot(aes(x = Site_Name, y = pH, fill = Site_Name),
               outlier.shape = NA) +
  geom_jitter(aes(x = Site_Name, y = pH), width = 0.1) +
  theme_bw() +
  scale_fill_brewer(palette = "Dark2") +
  labs(x = element_blank(), y = "pH", fill = "WWTP") +
  scale_x_discrete(guide = guide_axis(n.dodge=3)) +
  theme(text = element_text(size = textsize)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

# Plot CrAssphage sub-plot
Crass_plt <- ggplot(corr_dt) + geom_boxplot(aes(x = Site_Name, y = Crass_gc_100, fill = Site_Name)) +
  geom_jitter(aes(x = Site_Name, y = Crass_gc_100), width = 0.1) +
  ylim(min = 0, max = 7500) +
  theme_bw() +
  theme(text = element_text(size = textsize)) +
  ylab(expression("CrAssphage (gc 100 ml"^{-1}~")")) +
  xlab(element_blank()) +
  theme(legend.position = "none") +
  scale_x_discrete(guide = guide_axis(n.dodge=3)) +
  scale_fill_brewer(palette = "Dark2") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

# Generate figure - note, sometimes this needs rerunning
ggarrange(Flow_plt, NH4_plt, PO4_plt, pH_plt, EC_plt,  Crass_plt,
          nrow = 2, ncol = 3, common.legend = T,
          legend = "bottom", align = "hv")

  #### Figure 3 - Plot CoV levels ####

# Set scaling coefficient for multiple axes
coeff <- 0.0025

# Plot figure
ggplot(master) + 
  geom_line(data = cases_all, aes(x = Date, y = Daily_Cases_100k/ coeff), colour = "grey") +
  geom_point(aes(x = Date, y = Deaths_100k / coeff), shape = 17, colour = "#808080") +
  geom_point(aes(x = Date, y = gc_100ml, colour = Site_Name, shape = cov_detect)) +
  scale_shape_manual(values = c(2,17)) +
  scale_y_continuous( name = "SARS-CoV-2 genome copies / 100 mL",
                      sec.axis = sec_axis(~.*coeff, name = "Deaths per 100,000 per week/\nCases per 100,000 per day")) +
  facet_wrap(~ Site_Name) +
  xlim(as.Date("01/03/2020", format = "%d/%m/%Y"), as.Date("01/08/2020", format = "%d/%m/%Y")) +
  theme_bw() +
  geom_hline(yintercept = N1_LoQ, color = "grey", linetype = "dashed") +
  geom_hline(yintercept = N1_LoD, color = "grey", linetype = "dotted") +
  geom_vline(xintercept = as.Date("2020-03-16"), color = "grey", linetype = "longdash") +
  theme(legend.position="bottom") +
  scale_colour_brewer(palette = "Dark2") +
  theme(text = element_text(size = 15)) + 
  labs(colour = "WWTP", shape = "Amplification")

# - Plot gc/ 100 mL against deaths for supplemental information
ggplot(corr_dt) + geom_point(aes(x = Deaths_100k, y = gc_100ml, colour = Site_Name),
                             shape = 17, size = 3) +
  theme_bw() +
  ylab(element_text("SARS-CoV-2\nGenome copies/ 100 mL")) +
  xlab("Deaths per 100,000") +
  theme(text = element_text(size = 15)) +
  labs(colour = "WWTP") + 
  scale_colour_brewer(palette = "Dark2") +
  theme_bw() +
  theme(legend.position="bottom") + facet_wrap(~Site_Name)

# - Plot gc/ 100 mL against cases for supplemental information
ggplot(corr_dt) + geom_point(aes(x = Daily_Cases_100k, y = gc_100ml, colour = Site_Name),
                             shape = 17, size = 3) +
  theme_bw() +
  xlab("Daily positive tests per 100,000") +
  ylab(element_text("SARS-CoV-2\nGenome copies/ 100 mL")) +
  theme(text = element_text(size = 15)) +
  scale_colour_brewer(palette = "Dark2") +
  theme(legend.position="bottom") +
  labs(colour = "WWTP") + facet_wrap(~Site_Name)

# - Plot deaths against cases for supplemental information
ggplot(corr_dt) + geom_point(aes(x = Daily_Cases_100k,
                                 y = Deaths_100k,
                                 colour = Site_Name),
                             shape = 17, size = 3) +
  theme_bw() +
  ylab(element_text("Deaths per 100,000")) +
  xlab("Daily positive tests per 100,000") +
  theme(text = element_text(size = 15)) +
  scale_colour_brewer(palette = "Dark2") +
  theme(legend.position="bottom") +
  labs(colour = "WWTP") + facet_wrap(~Site_Name)

# Calculate weekly cases for correlations
corr_dt <- corr_dt %>% group_by(Site_Name, Week) %>%
  mutate(Weekly_Cases = sum(Daily_Cases_100k))

corr_dt <- as.data.table(corr_dt)

corr_dt <- corr_dt[is.na(corr_dt$gc_100ml) == F,]

#### Plot Correllogram ####
# - Create a function for generating correllograms per site
CorrMe <- function(site){
  dt<- corr_dt
  dt <- dt[dt$Site_Name == site]
  dt <- select(dt, c(Deaths_100k, Rolling_Sum_100k, gc_100ml, Week))
  dt <- dplyr::rename(dt, Deaths = Deaths_100k, Cases = Rolling_Sum_100k, CoV = gc_100ml)
  dt <- sapply(dt, as.numeric)
  mx <- rcorr(as.matrix(dt), type = "spearman")
  corrplot(mx$r, type = "upper", order = "alphabet",
          p.mat = mx$P, sig.level = 0.05, insig = "blank",
          col=c("black","white"),
          bg = "#DCDCDC",
          method = "pie",
#          col = brewer.pal(n = 8, name = "RdBu"),
          tl.col="black", tl.srt=45, tl.cex = 1.5, title = site,
          mar=c(5,0,2,0),
          cex.main = 2,
          cl.cex = 1.5)
}

# - Create grid for plots
par(mfrow=c(2,3))

sapply(c("Cardiff", "Gwynedd", "Liverpool","Manchester","The Wirral","Wrexham"),function(x) CorrMe(x))


study_wk = c(1,2,3,4,5,6,7)
wk_num <- study_wk + 13
myDate=paste("2020/",as.character(wk_num), "/1", sep = "")
as.Date("2020/13/1",format="%Y/%V/%u")

as.Date("2020-1-1", format ="%Y-%W-%u")

strftime("2020-03-30", format = "%V")


#### Test correction factors on correlations ####
# Create a separate data frame
corr_test_dt <- corr_dt
# Adjust Daily cases by ratio of population equivalent to local authority population
corr_test_dt$Daily_Cases_100k_adj <- corr_test_dt$Daily_Cases_100k * (corr_test_dt$Population_Equivalent/corr_test_dt$Local_Authority_Population)

# Adjust Daily cases by ratio of population equivalent to local authority population
corr_test_dt$Deaths_100k_adj <- corr_test_dt$Deaths_100k * (corr_test_dt$Population_Equivalent/corr_test_dt$Local_Authority_Population)

# Normalise CoV wastewater quantities by flow
corr_test_dt$gc_100ml_adj <- corr_test_dt$gc_100ml / corr_test_dt$Flow

# - Create a function for generating correllograms per site
CorrAdj <- function(site){
  dt<- corr_test_dt
  dt <- dt[dt$Site_Name == site]
  dt <- select(dt, c(Deaths_100k, Deaths_100k_adj,
                     Daily_Cases_100k, Daily_Cases_100k_adj,
                     gc_100ml,
                     gc_100ml_adj))
  dt <- dplyr::rename(dt,
               Deaths = Deaths_100k,
               Tests = Daily_Cases_100k,
               CoV = gc_100ml,
               `Deaths*` = Deaths_100k_adj,
               `Tests*` = Daily_Cases_100k_adj,
               `CoV*` = gc_100ml_adj)
  dt <- sapply(dt, as.numeric)
  mx <- rcorr(as.matrix(dt), type = "spearman")
  corrplot(mx$r, type = "upper", order = "alphabet",
           p.mat = mx$P, sig.level = 0.05, insig = "blank",
           col=c("black","white"),
           bg = "#DCDCDC",
           method = "pie",
           #          col = brewer.pal(n = 8, name = "RdBu"),
           tl.col="black", tl.srt=45, tl.cex = 1.5, title = site,
           mar=c(5,0,2,0),
           cex.main = 2,
           cl.cex = 1.5)
}

# - Create grid for plots
par(mfrow=c(2,3))

# - Generate correllograms
sapply(c("Cardiff", "Gwynedd", "Liverpool","Manchester","The Wirral","Wrexham"),function(x) CorrAdj(x))

#### Lag Analysis ####
# Note - this code is ugly, but it works!
library(plyr)
# Create a vector of sites, offset from wastewater sampling day and days of cases to sum to consider
site_vector <- c("Liverpool", "Manchester", "The Wirral")
window_size_vector <- c(1:14)
window_size <- 5

# Create a copy of the data to perform the lag analysis on
lag_dt <- master
lag_dt$gc_100ml[lag_dt$gc_100ml < N1_LoD] <- 0
lag_dt$gc_100ml[lag_dt$gc_100ml >= N1_LoD & lag_dt$gc_100ml < N1_LoQ] <- 0.1

# Create a function to calculate the Spearman p-value for correlation between wastewater signal and cases
win_size_p <- function(window_size, site){
  dt <- lag_dt[lag_dt$Site_Name == site]
  dt <- dt[order(Date)]
  dt <- select(dt, c(Date, gc_100ml,Daily_Cases_100k))
  dt$rollingsum <- rollsum(shift(dt$Daily_Cases_100k, n = offset, type = direction), window_size, fill = NA, align = "right")
  correlation <- cor.test(dt$gc_100ml, dt$rollingsum, method = "spearman", use = "complete.obs")
  p <- correlation[3]
  rho <- correlation[4]
  return(p)
}

# Create a function to calculate the Spearman rho for correlation between wastewater signal and cases
win_size_rho <- function(window_size, site){
  dt <- lag_dt[lag_dt$Site_Name == site]
  dt <- dt[order(Date)]
  dt <- select(dt, c(Date, gc_100ml,Daily_Cases_100k))
  dt$rollingsum <- rollsum(shift(dt$Daily_Cases_100k, n = offset, type = direction), window_size, fill = NA, align = "right")
  correlation <- cor.test(dt$gc_100ml, dt$rollingsum, method = "spearman", use = "complete.obs")
  p <- correlation[3]
  rho <- correlation[4]
  return(rho)
}

# Create a function to calculate the range of possible p-values and rho over a range of window sizes
win_size_func <- function(site){
  dt_p = ldply(mapply(win_size_p, window_size_vector, site, SIMPLIFY = T), data.table)
  dt_rho <- ldply(mapply(win_size_rho, window_size_vector, site, SIMPLIFY = T), data.table)
  dt_p$win_size <- window_size_vector
  dt_rho$win_size <- window_size_vector
  dt_p$p_value <- dt_p$V1
  dt_rho$rho <- dt_rho$V1
  dt_p <- select(dt_p, win_size, p_value)
  dt_rho <- select(dt_rho, win_size, rho)
  dt_full <- full_join(dt_p, dt_rho)
  dt_full$site <- site
  dt_full$win_offset <- offset
  return(dt_full)
}

# Create a data table to store the results
dt_all <- data.table()

# Set direction of lag/ lead to test
direction <- "lead"

# Calculate for zero offset
offset <- 0
for(site in site_vector){
  dt_site <- win_size_func(site)
  dt_all <- rbind(dt_all, dt_site)
}

# Calculate for leading offsets
offset_vector <- c(1:14)
for(offset in offset_vector){
  for(site in site_vector){
    dt_site <- win_size_func(site)
    dt_all <- rbind(dt_all, dt_site)
  }
}

# Change direction to lag
direction <- "lag"
offset_vector <- c(1:5)
# Calculate for lagging offsets
for(offset in offset_vector){
  for(site in site_vector){
    dt_site <- win_size_func(site)
    dt_site$win_offset <- dt_site$win_offset*-1
    dt_all <- rbind(dt_all, dt_site)
  }
}

# Create columns of rho values
dt_all$rho_adj <- dt_all$rho
dt_all$rho_sig <- dt_all$rho

# Reject rho where p >= 0.05
dt_all$rho_sig[dt_all$p_value >= 0.05] <- NA

# Calculate FDR corrected p values
dt_all$p_adj <- as.numeric(p.adjust(dt_all$p_value, method = "fdr", n = length(dt_all$p_value)))

# Reject rho where FDR corrected p >= 0.05
dt_all$rho_adj[dt_all$p_adj >= 0.05] <- NA

# Plot figure

#dt_all$rho_adj <- round(dt_all$rho_adj/0.05)*0.05

ggplot(dt_all) +
  geom_tile(aes(x = win_offset, y = win_size, fill = rho_adj)) +
  theme_bw() + 
  scale_fill_continuous(type = "viridis") +
  coord_equal() +
  facet_wrap(~site) +
#  xlim(min = -5, max = 14) +
#  ylim(min = 0, max = 14) +
  labs( x = "Number of days offset from wastewater sampling", y = "Days to sum cases over", fill = "Spearman's rho") +
  theme(legend.position = "bottom") +
  theme(text = element_text(size = 15))
  




