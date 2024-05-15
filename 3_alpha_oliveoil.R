# Title: Cross-sectional associations of olive oil consumption with human gut microbiota diversity and composition
# Author: Jiaqi Ni [ORCID](https://orcid.org/0000-0001-5625-0373)
# Description: This is the analysis script for alpha diversity analyses

########################################################### STARTUP ################################################################

# Set up working environment  ------------------------------------------------------------------
rm(list=ls())
getwd()
setwd("/home/jiaqi/Jiaqi_Ni/OliveMicrobiomeCognition_Microbiome")
options(max.print=1000000)
options(scipen = 999)

if (sys.nframe() == 0L) {
  this.dir <- "/home/jiaqi/Jiaqi_Ni/OliveMicrobiomeCognition_Microbiome"
} else {
  freeze <- ls()
  this.dir <- source.dir
}

# Load dependencies -------------------------------------------------------
library(mia)
library(miaTime)
library(miaViz)
library(scater) # Load package to plot reducedDim
library(vegan)
library(patchwork)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(gridExtra)
library(cowplot)
library(Maaslin2)
library(ggalluvial)
library(tidyverse)
library(modelsummary)
library(RColorBrewer)
library(cowplot)
library(grid)
# Data preparation ---------------------------------------------------
# * load data ---------------------------------------------------
ps <- readRDS("./DATA/jiaqi_Phyloseq_110324.rds") #baseline microbiota N=656
## convert phyloseq to TSE
tse <- makeTreeSummarizedExperimentFromPhyloseq(ps) 
tse
# * calculate alpha diversity measures ---------------------------------------------------
#Richness
tse <- mia::estimateRichness(tse, 
                             index = c("observed", "chao1", "ace", "hill"), 
                             name=c("Observed", "Chao1", "ACE", "Hill"))

#Evenness
tse <- estimateEvenness(tse, 
                        index = c("simpson", "pielou", "camargo", "simpson_evenness", "evar", "bulla"),
                        name = c("Simpson", "Pielou", "Camargo", "Simpson_evenness", "Evar", "Bulla"))

#Diversity
tse <- estimateDiversity(tse, 
                         index = c("shannon", "gini_simpson", "inverse_simpson", "coverage", "fisher", "faith", "log_modulo_skewness"),
                         name = c("Shannon", "GiniSimpson",  "InverseSimpson",  "Coverage", "Fisher", "Faith", "LogModSkewness"))

#Dominance
tse <- estimateDominance(tse, 
                         index = c("dbp", "dmn", "absolute", "relative", "simpson_lambda", "core_abundance", "gini"),
                         name  = c("BergerParker", "McNaughton", "Absolute", "Relative", "SimpsonLambda", "CoreAbundance", "Gini"))

#Divergence
tse <- mia::estimateDivergence(tse,
                               assay.type = "counts",
                               reference = "median",
                               FUN = vegan::vegdist)

# * extract data for dataframe ---------------------------------------------------
df <- as.data.frame(colData(tse))
df<-subset(df, select = -c(Chao1_se,ACE, ACE_se))
which(colnames(df)=="Observed")
alphaDiversity <-colnames(df)[169:ncol(df)]
# Define and fix continuous and categorical variables
catVars <- c("sexo_s1", "nodo", "area", "grupo_int_v00", "edu_level_v00", "civil_status_v00",
             "ter_ajust_olivatot_v00", "ter_ajust_ac_olivavir_v00", "ter_ajust_re_oliva_v00", 
              "diab_prev_s1", "hta_v00", "colest_v00", "depression_v00",  "antihypertensive_v00", "antidepressant_v00", 
             "antidiabetic_v00", "hypolipidemic_v00", "ter_exercise_day_v00", 
             "smoke_status_v00","ter_MED12score_v00","ter_alcoholg_v00")
df[,catVars] <- lapply(df[,catVars], as.factor)

conVars <- c("edad_s1", "olivatot_v00",  "ac_olivavir_v00", "re_oliva_v00",  
             "ajust_olivatot_v00", "ajust_ac_olivavir_v00", "ajust_re_oliva_v00", 
             "ptrend_ter_ajust_olivatot_v00", "ptrend_ter_ajust_ac_olivavir_v00", "ptrend_ter_ajust_re_oliva_v00",
             "imc_v00", "waist_cir_v00", "exercise_day_v00", "alcoholg_v00", "qalcoholg_v00", "MED12score_v00")
df[,conVars] <- lapply(df[,conVars], as.numeric)

########################################################### ANALYSIS ################################################################

# MEAN COMPARISONS ---------------------------------------------------
# * normality test ---------------------------------------------------
shapiro_results <- list()
for(i in alphaDiversity) {
  variables <- df[[i]]
  if(length(unique(variables))>1){
    shapiro_result <- shapiro.test(variables)
    shapiro_results[[i]] <- shapiro_result
  } else {
    cat("Skipped Shapiro-wilk test for", i, "due to no variability. \n\n")
  }
  
}
for(i in alphaDiversity) {
  if (!is.null(shapiro_results[[i]])){
    cat("Shapiro-Wilk test for", i, ":\n")
    print(shapiro_results[[i]])
    cat("\n")
  }
  
}
#BOX+VIOLIN Plots ------------------------------------------------------------- 
# * total olive oil ---------------------------------------------------
comb <- split(t(combn(levels(df$ter_ajust_olivatot_v00), 2)), 
              seq(nrow(t(combn(levels(df$ter_ajust_olivatot_v00), 2)))))

alphaplots<- list()
for(i in alphaDiversity) {
  p <-ggplot(df, aes(x=ter_ajust_olivatot_v00, y = !!sym(i), fill= ter_ajust_olivatot_v00)) + 
    geom_violin(trim = F, width=.5, alpha=.3)+ 
    geom_boxplot(width=.25, alpha=.3) +
    geom_jitter(width = .05,  alpha=.5) +
    stat_compare_means(comparisons = comb, method= "wilcox.test", map_signif_level = FALSE, size = 5, family = "Times New Roman") +
    labs(x = "",
         y = i) +
    theme_test() +
    guides(fill=guide_legend(title = "TOO tertiles")) + 
    scale_fill_hue(labels=c("T1","T2","T3")) +
    theme(text = element_text(size = 30, family = "Times New Roman"),
          legend.text = element_text(size = 30)) +
    scale_fill_brewer(palette="YlOrBr")
  alphaplots[[i]] <- p
}

for(i in seq_along(alphaplots)) {
  print(alphaplots[[i]])
}

#plot 
figTOO <- wrap_plots(alphaplots$Chao1, alphaplots$Simpson, alphaplots$Shannon, alphaplots$InverseSimpson) +
  plot_layout(guides = "collect") #generate multi-panel plot
figTOO
# * extra virgin olive oil ---------------------------------------------------
comb <- split(t(combn(levels(df$ter_ajust_ac_olivavir_v00), 2)), 
              seq(nrow(t(combn(levels(df$ter_ajust_ac_olivavir_v00), 2)))))
alphaplots<- list()
for(i in alphaDiversity) {
  p <-ggplot(df, aes(x=ter_ajust_ac_olivavir_v00, y = !!sym(i), fill= ter_ajust_ac_olivavir_v00)) + 
    geom_violin(trim = F, width=.5, alpha=.3)+ 
    geom_boxplot(width=.25, alpha=.3) +
    geom_jitter(width = .05,  alpha=.5) +
    stat_compare_means(comparisons = comb, method= "wilcox.test", map_signif_level = FALSE, size = 5, family = "Times New Roman") +
    labs(x = "",
         y = i) +
    theme_test() +    
    guides(fill=guide_legend(title = "VOO tertiles")) + 
    scale_fill_hue(labels=c("T1","T2","T3")) +
    theme(text = element_text(size = 30, family = "Times New Roman"),
          legend.text = element_text(size = 30)) +
    scale_fill_brewer(palette="Greens")
  
  alphaplots[[i]] <- p
}

for(i in seq_along(alphaplots)) {
  print(alphaplots[[i]])
}

#plot 
figVOO <-wrap_plots(alphaplots$Chao1, alphaplots$Simpson, alphaplots$Shannon, alphaplots$InverseSimpson) +
  plot_layout(guides = "collect") #generate multi-panel plot

# * common olive oil ---------------------------------------------------
comb <- split(t(combn(levels(df$ter_ajust_re_oliva_v00), 2)), 
              seq(nrow(t(combn(levels(df$ter_ajust_re_oliva_v00), 2)))))
alphaplots<- list()
for(i in alphaDiversity) {
  p <-ggplot(df, aes(x=ter_ajust_re_oliva_v00, y = !!sym(i), fill= ter_ajust_re_oliva_v00)) + 
    geom_violin(trim = F, width=.5, alpha=.3)+ 
    geom_boxplot(width=.25, alpha=.3) +
    geom_jitter(width = .05,  alpha=.5) +
    stat_compare_means(comparisons = comb, method= "wilcox.test", map_signif_level = FALSE, size = 5, family = "Times New Roman") +
    labs(x = "",
         y = i) +
    theme_test() +
    guides(fill=guide_legend(title = "COO tertiles")) + 
    scale_fill_hue(labels=c("T1","T2","T3")) +
    theme(text = element_text(size = 30, family = "Times New Roman"),
          legend.text = element_text(size = 30)) +
    scale_fill_brewer(palette="Reds")
  
  alphaplots[[i]] <- p
}

for(i in seq_along(alphaplots)) {
  print(alphaplots[[i]])
}
#plot 
figCOO <- wrap_plots(alphaplots$Chao1, alphaplots$Simpson, alphaplots$Shannon, alphaplots$InverseSimpson) +
  plot_layout(guides = "collect") #generate multi-panel plot

## COMBINE 3 FIGURES
combined_alpha <- plot_grid(figTOO, figVOO, figCOO, nrow = 3)
combined_alpha

## Save figure
png(file="/home/jiaqi/Jiaqi_Ni/OliveMicrobiomeCognition_Microbiome/Output/alpha.png",width=2100,height=2500, pointsize=50)
combined_alpha
dev.off()


# LINEAR REGRESSION ---------------------------------------------------
# Basic model ------------------------------------------------------------- 
# baseline age, sex
## Set variables for Loop
alphaindex <-c("Chao1", "Simpson", "Shannon", "InverseSimpson")

### * total olive oil  (cat.) -------------------------------------------------------------
baselinebasicresults <- list()
for (i in alphaindex){
  formula <- as.formula(paste(i, "~ ter_ajust_olivatot_v00 +  edad_s1 + sexo_s1"))
  baselinebasicmodels <- lm ( formula, data = df)
  baselinebasicresults[[i]] <- summary(baselinebasicmodels)
  
}

for (i in alphaindex){
  print(paste("Regression results for", i))
  print(baselinebasicresults[[i]])
  
}

modelsummary(baselinebasicresults, 
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "p.value")
### total olive oil  (p-trend)
baselinebasicresults <- list()
for (i in alphaindex){
  formula <- as.formula(paste(i, "~ ptrend_ter_ajust_olivatot_v00 +  edad_s1 + sexo_s1"))
  baselinebasicmodels <- lm ( formula, data = df)
  baselinebasicresults[[i]] <- summary(baselinebasicmodels)
  
}

for (i in alphaindex){
  print(paste("Regression results for", i))
  print(baselinebasicresults[[i]])
  
}

modelsummary(baselinebasicresults, 
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "p.value")

### * total olive oil  (cont.) -------------------------------------------------------------
baselinebasicresults <- list()
for (i in alphaindex){
  formula <- as.formula(paste(i, "~ ajust_tspolivatot_v00 +  edad_s1 + sexo_s1"))
  baselinebasicmodels <- lm ( formula, data = df)
  baselinebasicresults[[i]] <- summary(baselinebasicmodels)
  
}

for (i in alphaindex){
  print(paste("Regression results for", i))
  print(baselinebasicresults[[i]])
  
}

modelsummary(baselinebasicresults, 
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "p.value")
### * extra virgin olive oil  (cat.) -------------------------------------------------------------
baselinebasicresults <- list()
for (i in alphaindex){
  formula <- as.formula(paste(i, "~ ter_ajust_ac_olivavir_v00 +  edad_s1 + sexo_s1"))
  baselinebasicmodels <- lm ( formula, data = df)
  baselinebasicresults[[i]] <- summary(baselinebasicmodels)
  
}

for (i in alphaindex){
  print(paste("Regression results for", i))
  print(baselinebasicresults[[i]])
  
}

modelsummary(baselinebasicresults, 
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "p.value")
### extra virgin olive oil  (p-trend)
baselinebasicresults <- list()
for (i in alphaindex){
  formula <- as.formula(paste(i, "~ ptrend_ter_ajust_ac_olivavir_v00 +  edad_s1 + sexo_s1"))
  baselinebasicmodels <- lm ( formula, data = df)
  baselinebasicresults[[i]] <- summary(baselinebasicmodels)
  
}

for (i in alphaindex){
  print(paste("Regression results for", i))
  print(baselinebasicresults[[i]])
  
}

modelsummary(baselinebasicresults, 
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "p.value")
### * extra virgin olive oil  (cont.) -------------------------------------------------------------
baselinebasicresults <- list()
for (i in alphaindex){
  formula <- as.formula(paste(i, "~ ajust_tspac_olivavir_v00 +  edad_s1 + sexo_s1"))
  baselinebasicmodels <- lm ( formula, data = df)
  baselinebasicresults[[i]] <- summary(baselinebasicmodels)
  
}

for (i in alphaindex){
  print(paste("Regression results for", i))
  print(baselinebasicresults[[i]])
  
}

modelsummary(baselinebasicresults, 
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "p.value")
### * common olive oil  (cat.) -------------------------------------------------------------
baselinebasicresults <- list()
for (i in alphaindex){
  formula <- as.formula(paste(i, "~ ter_ajust_re_oliva_v00 +  edad_s1 + sexo_s1"))
  baselinebasicmodels <- lm ( formula, data = df)
  baselinebasicresults[[i]] <- summary(baselinebasicmodels)
  
}

for (i in alphaindex){
  print(paste("Regression results for", i))
  print(baselinebasicresults[[i]])
  
}

modelsummary(baselinebasicresults, 
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "p.value")
### refined olive oil  (p-trend)
baselinebasicresults <- list()
for (i in alphaindex){
  formula <- as.formula(paste(i, "~ ptrend_ter_ajust_re_oliva_v00 +  edad_s1 + sexo_s1"))
  baselinebasicmodels <- lm ( formula, data = df)
  baselinebasicresults[[i]] <- summary(baselinebasicmodels)
  
}

for (i in alphaindex){
  print(paste("Regression results for", i))
  print(baselinebasicresults[[i]])
  
}

modelsummary(baselinebasicresults, 
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "p.value")
### * common olive oil  (cont.) -------------------------------------------------------------
baselinebasicresults <- list()
for (i in alphaindex){
  formula <- as.formula(paste(i, "~ ajust_tspre_oliva_v00 +  edad_s1 + sexo_s1"))
  baselinebasicmodels <- lm ( formula, data = df)
  baselinebasicresults[[i]] <- summary(baselinebasicmodels)
  
}

for (i in alphaindex){
  print(paste("Regression results for", i))
  print(baselinebasicresults[[i]])
  
}

modelsummary(baselinebasicresults, 
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "p.value")
# Socio-demographically adjusted model ------------------------------------------------------------- 
# baseline age, sex, geographical area, educational level, civil status
### * total olive oil  (cat.) -------------------------------------------------------------
baselinebasicresults <- list()
for (i in alphaindex){
  formula <- as.formula(paste(i, "~ ter_ajust_olivatot_v00 +  edad_s1 + sexo_s1 + area + edu_level_v00 + civil_status_v00 "))
  baselinebasicmodels <- lm ( formula, data = df)
  baselinebasicresults[[i]] <- summary(baselinebasicmodels)
  
}

for (i in alphaindex){
  print(paste("Regression results for", i))
  print(baselinebasicresults[[i]])
  
}

modelsummary(baselinebasicresults, 
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "p.value")
### total olive oil  (p-trend)
baselinebasicresults <- list()
for (i in alphaindex){
  formula <- as.formula(paste(i, "~ ptrend_ter_ajust_olivatot_v00 +  edad_s1 + sexo_s1 + area + edu_level_v00 + civil_status_v00"))
  baselinebasicmodels <- lm ( formula, data = df)
  baselinebasicresults[[i]] <- summary(baselinebasicmodels)
  
}

for (i in alphaindex){
  print(paste("Regression results for", i))
  print(baselinebasicresults[[i]])
  
}

modelsummary(baselinebasicresults, 
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "p.value")

### * total olive oil  (cont.) -------------------------------------------------------------
baselinebasicresults <- list()
for (i in alphaindex){
  formula <- as.formula(paste(i, "~ ajust_tspolivatot_v00 +  edad_s1 + sexo_s1 + area + edu_level_v00 + civil_status_v00"))
  baselinebasicmodels <- lm ( formula, data = df)
  baselinebasicresults[[i]] <- summary(baselinebasicmodels)
  
}

for (i in alphaindex){
  print(paste("Regression results for", i))
  print(baselinebasicresults[[i]])
  
}

modelsummary(baselinebasicresults, 
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "p.value")
### * extra virgin olive oil  (cat.) -------------------------------------------------------------
baselinebasicresults <- list()
for (i in alphaindex){
  formula <- as.formula(paste(i, "~ ter_ajust_ac_olivavir_v00 +  edad_s1 + sexo_s1 + area + edu_level_v00 + civil_status_v00"))
  baselinebasicmodels <- lm ( formula, data = df)
  baselinebasicresults[[i]] <- summary(baselinebasicmodels)
  
}

for (i in alphaindex){
  print(paste("Regression results for", i))
  print(baselinebasicresults[[i]])
  
}

modelsummary(baselinebasicresults, 
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "p.value")
### extra virgin olive oil  (p-trend)
baselinebasicresults <- list()
for (i in alphaindex){
  formula <- as.formula(paste(i, "~ ptrend_ter_ajust_ac_olivavir_v00 +  edad_s1 + sexo_s1 + area + edu_level_v00 + civil_status_v00"))
  baselinebasicmodels <- lm ( formula, data = df)
  baselinebasicresults[[i]] <- summary(baselinebasicmodels)
  
}

for (i in alphaindex){
  print(paste("Regression results for", i))
  print(baselinebasicresults[[i]])
  
}

modelsummary(baselinebasicresults, 
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "p.value")
### * extra virgin olive oil  (cont.) -------------------------------------------------------------
baselinebasicresults <- list()
for (i in alphaindex){
  formula <- as.formula(paste(i, "~ ajust_tspac_olivavir_v00 +  edad_s1 + sexo_s1 + area + edu_level_v00 + civil_status_v00"))
  baselinebasicmodels <- lm ( formula, data = df)
  baselinebasicresults[[i]] <- summary(baselinebasicmodels)
  
}

for (i in alphaindex){
  print(paste("Regression results for", i))
  print(baselinebasicresults[[i]])
  
}

modelsummary(baselinebasicresults, 
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "p.value")
### * common olive oil  (cat.) -------------------------------------------------------------
baselinebasicresults <- list()
for (i in alphaindex){
  formula <- as.formula(paste(i, "~ ter_ajust_re_oliva_v00 +  edad_s1 + sexo_s1 + area + edu_level_v00 + civil_status_v00"))
  baselinebasicmodels <- lm ( formula, data = df)
  baselinebasicresults[[i]] <- summary(baselinebasicmodels)
  
}

for (i in alphaindex){
  print(paste("Regression results for", i))
  print(baselinebasicresults[[i]])
  
}

modelsummary(baselinebasicresults, 
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "p.value")
### refined olive oil  (p-trend)
baselinebasicresults <- list()
for (i in alphaindex){
  formula <- as.formula(paste(i, "~ ptrend_ter_ajust_re_oliva_v00 +  edad_s1 + sexo_s1 + area + edu_level_v00 + civil_status_v00"))
  baselinebasicmodels <- lm ( formula, data = df)
  baselinebasicresults[[i]] <- summary(baselinebasicmodels)
  
}

for (i in alphaindex){
  print(paste("Regression results for", i))
  print(baselinebasicresults[[i]])
  
}

modelsummary(baselinebasicresults, 
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "p.value")
### * common olive oil  (cont.) -------------------------------------------------------------
baselinebasicresults <- list()
for (i in alphaindex){
  formula <- as.formula(paste(i, "~ ajust_tspre_oliva_v00 +  edad_s1 + sexo_s1 + area + edu_level_v00 + civil_status_v00"))
  baselinebasicmodels <- lm ( formula, data = df)
  baselinebasicresults[[i]] <- summary(baselinebasicmodels)
  
}

for (i in alphaindex){
  print(paste("Regression results for", i))
  print(baselinebasicresults[[i]])
  
}

modelsummary(baselinebasicresults, 
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "p.value")
# Fully adjusted model ------------------------------------------------------------- 
#baseline age, sex, geographical area, educational level, civil status
#BMI, physical activity, smoking status, alcohol consumption, depressive symptomatology, diabetes prevalence, hypertension prevalence, hypercholesterolemia prevalence, 
#adherence to Mediterranean diet (12 scores)
### * total olive oil  (cat.) -------------------------------------------------------------
baselinebasicresults <- list()
for (i in alphaindex){
  formula <- as.formula(paste(i, "~ ter_ajust_olivatot_v00 +  edad_s1 + sexo_s1 + area + edu_level_v00 + civil_status_v00 + 
                              imc_v00 + exercise_day_v00 + smoke_status_v00 + alcoholg_v00 + qalcoholg_v00 + 
                              depression_v00 + diab_prev_s1 + hta_v00 + colest_v00 + MED12score_v00"))
  baselinebasicmodels <- lm (formula, data = df)
  baselinebasicresults[[i]] <- summary(baselinebasicmodels)
  
}

for (i in alphaindex){
  print(paste("Regression results for", i))
  print(baselinebasicresults[[i]])
  
}

modelsummary(baselinebasicresults, 
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "p.value")

### total olive oil  (p-trend)
baselinebasicresults <- list()
for (i in alphaindex){
  formula <- as.formula(paste(i, "~ ptrend_ter_ajust_olivatot_v00 +  edad_s1 + sexo_s1 + area + edu_level_v00 + civil_status_v00 + 
                              imc_v00 + exercise_day_v00 + smoke_status_v00 + alcoholg_v00 + qalcoholg_v00 + 
                              depression_v00 + diab_prev_s1 + hta_v00 + colest_v00 +  MED12score_v00"))
  baselinebasicmodels <- lm ( formula, data = df)
  baselinebasicresults[[i]] <- summary(baselinebasicmodels)
  
}

for (i in alphaindex){
  print(paste("Regression results for", i))
  print(baselinebasicresults[[i]])
  
}

modelsummary(baselinebasicresults, 
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "p.value")

### * total olive oil  (cont.) -------------------------------------------------------------
baselinebasicresults <- list()
for (i in alphaindex){
  formula <- as.formula(paste(i, "~ ajust_tspolivatot_v00 +  edad_s1 + sexo_s1 + area + edu_level_v00 + civil_status_v00 + 
                              imc_v00 + exercise_day_v00 + smoke_status_v00 + alcoholg_v00 + qalcoholg_v00 + 
                              depression_v00 + diab_prev_s1 + hta_v00 + colest_v00 +  MED12score_v00"))
  baselinebasicmodels <- lm ( formula, data = df)
  baselinebasicresults[[i]] <- summary(baselinebasicmodels)
  
}

for (i in alphaindex){
  print(paste("Regression results for", i))
  print(baselinebasicresults[[i]])
  
}

modelsummary(baselinebasicresults, 
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "p.value")
### * extra virgin olive oil  (cat.) -------------------------------------------------------------
baselinebasicresults <- list()
for (i in alphaindex){
  formula <- as.formula(paste(i, "~ ter_ajust_ac_olivavir_v00 +  edad_s1 + sexo_s1 + area + edu_level_v00 + civil_status_v00 + 
                              imc_v00 + exercise_day_v00 + smoke_status_v00 + alcoholg_v00 + qalcoholg_v00 + 
                              depression_v00 + diab_prev_s1 + hta_v00 + colest_v00 +  MED12score_v00"))
  baselinebasicmodels <- lm ( formula, data = df)
  baselinebasicresults[[i]] <- summary(baselinebasicmodels)
  
}

for (i in alphaindex){
  print(paste("Regression results for", i))
  print(baselinebasicresults[[i]])
  
}

modelsummary(baselinebasicresults, 
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "p.value")
### extra virgin olive oil  (p-trend)
baselinebasicresults <- list()
for (i in alphaindex){
  formula <- as.formula(paste(i, "~ ptrend_ter_ajust_ac_olivavir_v00 + edad_s1 + sexo_s1 + area + edu_level_v00 + civil_status_v00 + 
                              imc_v00 + exercise_day_v00 + smoke_status_v00 + alcoholg_v00 + qalcoholg_v00 + 
                              depression_v00 + diab_prev_s1 + hta_v00 + colest_v00 +  MED12score_v00"))
  baselinebasicmodels <- lm ( formula, data = df)
  baselinebasicresults[[i]] <- summary(baselinebasicmodels)
  
}

for (i in alphaindex){
  print(paste("Regression results for", i))
  print(baselinebasicresults[[i]])
  
}

modelsummary(baselinebasicresults, 
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "p.value")
### * extra virgin olive oil  (cont.) -------------------------------------------------------------
baselinebasicresults <- list()
for (i in alphaindex){
  formula <- as.formula(paste(i, "~ ajust_tspac_olivavir_v00 +  edad_s1 + sexo_s1 + area + edu_level_v00 + civil_status_v00 + 
                              imc_v00 + exercise_day_v00 + smoke_status_v00 + alcoholg_v00 + qalcoholg_v00 + 
                              depression_v00 + diab_prev_s1 + hta_v00 + colest_v00 +  MED12score_v00"))
  baselinebasicmodels <- lm ( formula, data = df)
  baselinebasicresults[[i]] <- summary(baselinebasicmodels)
  
}

for (i in alphaindex){
  print(paste("Regression results for", i))
  print(baselinebasicresults[[i]])
  
}

modelsummary(baselinebasicresults, 
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "p.value")

### * common olive oil  (cat.) -------------------------------------------------------------
baselinebasicresults <- list()
for (i in alphaindex){
  formula <- as.formula(paste(i, "~ ter_ajust_re_oliva_v00 +  edad_s1 + sexo_s1 + area + edu_level_v00 + civil_status_v00 + 
                              imc_v00 + exercise_day_v00 + smoke_status_v00 + alcoholg_v00 + qalcoholg_v00 + 
                              depression_v00 + diab_prev_s1 + hta_v00 + colest_v00 + MED12score_v00"))
  baselinebasicmodels <- lm ( formula, data = df)
  baselinebasicresults[[i]] <- summary(baselinebasicmodels)
  
}

for (i in alphaindex){
  print(paste("Regression results for", i))
  print(baselinebasicresults[[i]])
  
}

modelsummary(baselinebasicresults, 
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "p.value")
### refined olive oil  (p-trend)
baselinebasicresults <- list()
for (i in alphaindex){
  formula <- as.formula(paste(i, "~ ptrend_ter_ajust_re_oliva_v00 +  edad_s1 + sexo_s1 + area + edu_level_v00 + civil_status_v00 + 
                              imc_v00 + exercise_day_v00 + smoke_status_v00 + alcoholg_v00 + qalcoholg_v00 + 
                              depression_v00 + diab_prev_s1 + hta_v00 + colest_v00 + MED12score_v00"))
  baselinebasicmodels <- lm ( formula, data = df)
  baselinebasicresults[[i]] <- summary(baselinebasicmodels)
  
}

for (i in alphaindex){
  print(paste("Regression results for", i))
  print(baselinebasicresults[[i]])
  
}

modelsummary(baselinebasicresults, 
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "p.value")
### * common olive oil  (cont.) -------------------------------------------------------------
baselinebasicresults <- list()
for (i in alphaindex){
  formula <- as.formula(paste(i, "~ ajust_tspre_oliva_v00 +  edad_s1 + sexo_s1 + area + edu_level_v00 + civil_status_v00 + 
                              imc_v00 + exercise_day_v00 + smoke_status_v00 + alcoholg_v00 + qalcoholg_v00 + 
                              depression_v00 + diab_prev_s1 + hta_v00 + colest_v00 +MED12score_v00"))
  baselinebasicmodels <- lm ( formula, data = df)
  baselinebasicresults[[i]] <- summary(baselinebasicmodels)
  
}

for (i in alphaindex){
  print(paste("Regression results for", i))
  print(baselinebasicresults[[i]])
  
}

modelsummary(baselinebasicresults, 
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "p.value")

