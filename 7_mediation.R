# Title: Mediation analysis
# Author: Jiaqi Ni [ORCID](https://orcid.org/0000-0001-5625-0373)
# Description: This is the analysis script for assessing the mediation effect of gut microbiota in the associations between olive oil consumption and cognitive change

########################################################### STARTUP ################################################################
# Set up working environment  ------------------------------------------------------------------
rm(list=ls())
getwd()
setwd("/home/jiaqi/Jiaqi_Ni/OliveMicrobiomeCognition_Microbiome/")
options(max.print=1000000)
options(scipen = 999)
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
library(mediation)
library(lmtest)
library(HIMA)
library(glue)
library(diagram)


########################################################### OVERALL COMPOSITION ################################################################

# Data preparation ---------------------------------------------------

# * load data ---------------------------------------------------
ps <- readRDS("/home/jiaqi/Jiaqi_Ni/OliveMicrobiomeCognition_Microbiome/DATA/jiaqi_Phyloseq_110324.rds") #baseline microbiota N=656
tse <- makeTreeSummarizedExperimentFromPhyloseq(ps) # convert phyloseq to TSE
tse
# * calculate alpha diversity  ---------------------------------------------------
#Richness
tse <- mia::estimateRichness(tse, 
                             index = c( "chao1"), 
                             name=c( "Chao1"))

#Evenness
tse <- estimateEvenness(tse, 
                        index = c("simpson"),
                        name = c("Simpson"))

#Diversity
tse <- estimateDiversity(tse, 
                         index = c("shannon", "inverse_simpson"),
                         name = c("Shannon", "InverseSimpson"))


# * calculate beta diversity  (Aitchison distance) ---------------------------------------------------
##CLR transform
tse <- transformAssay(x = tse, assay.type = "counts", method = "clr", 
                      pseudocount = 1, name = "clr")

# Calculate Aitchison distance (euclidean distance of CLR-transformed species abundances)
beta <- vegan::vegdist(t(assays(tse)$clr), method = "euclidean")
# Take a PCA and choose the first 2 principal components
pca <- cmdscale(beta, k = 2, eig = TRUE)
##with alpha indexes
df <- as.data.frame(colData(tse))
df<-subset(df, select = -c(Chao1_se))

##with pc1 & pc2
pcadf <- as.data.frame(pca$points)
names(pcadf) <- c("PC1", "PC2")
pcadf$SampleIDcounts <- rownames(pcadf)

##final dataframe
df <- merge(df, pcadf, by = "SampleIDcounts")
which(colnames(df)=="Chao1")
df[, 169:ncol(df)] <- scale(df[, 169:ncol(df)]) #scale the continuous indices

# Define and fix continuous and categorical variables
catVars <- c("sexo_s1",  "nodo", "area", "grupo_int_v00", "edu_level_v00", "civil_status_v00",
             "ter_ajust_olivatot_v00", "ter_ajust_ac_olivavir_v00", "ter_ajust_re_oliva_v00", 
              "diab_prev_s1", "hta_v00", "colest_v00", "depression_v00",  "antihypertensive_v00", "antidepressant_v00", 
             "antidiabetic_v00", "hypolipidemic_v00", "ter_exercise_day_v00", 
             "smoke_status_v00","ter_MED12score_v00","ter_alcoholg_v00")
df[,catVars] <- lapply(df[,catVars], as.factor)

conVars <- c("edad_s1", "olivatot_v00",  "ac_olivavir_v00", "re_oliva_v00",  "ajust_tspolivatot_v00", "ajust_tspac_olivavir_v00", "ajust_tspre_oliva_v00",
             "ajust_olivatot_v00", "ajust_ac_olivavir_v00", "ajust_re_oliva_v00", 
             "ptrend_ter_ajust_olivatot_v00", "ptrend_ter_ajust_ac_olivavir_v00", "ptrend_ter_ajust_re_oliva_v00",
             "imc_v00", "waist_cir_v00", "exercise_day_v00", "alcoholg_v00", "qalcoholg_v00", "MED12score_v00",
             "Chao1","Simpson" ,"Shannon","InverseSimpson","PC1","PC2")
df[,conVars] <- lapply(df[,conVars], as.numeric)
# vectors of exposure, mediator, and outcome names 
exposures <- c("ajust_tspolivatot_v00", "ajust_tspac_olivavir_v00", "ajust_tspre_oliva_v00") 
mediators <- c("Chao1", "Simpson", "Shannon", "InverseSimpson", "PC1", "PC2") 


# Fully adjusted model ---------------------------------------------------
# baseline age, sex, baseline cognitive test score, intervention group, geographical area, educational level, civil status
#BMI, physical activity, smoking status, alcohol consumption, depressive symptomatology, diabetes prevalence, hypertension prevalence, hypercholesterolemia prevalence, 
#adherence to Mediterranean diet (12 scores)

##GCF---------------------------------------------------
# Create an empty data frame to store results 
mediation_results_gcf <- data.frame(   
  exposure = character(),   
  mediator = character(),   
  ACME_estimate = numeric(),   ACME_lower95 = numeric(),   ACME_upper95 = numeric(),   ACME_p = numeric(),   
  ADE_estimate = numeric(),   ADE_lower95 = numeric(),   ADE_upper95 = numeric(),   ADE_p = numeric(),   
  total_estimate = numeric(),   total_lower95 = numeric(),   total_upper95 = numeric(),   total_p = numeric(),   
  prop_estimate = numeric(),   prop_lower95 = numeric(),   prop_upper95 = numeric(),   prop_p = numeric(),   stringsAsFactors = FALSE ) 


set.seed(2023) 
# Nested loops for each combination of exposure, mediator, and outcome 
for (exposure in exposures) {   
  for (mediator in mediators) {     
    med_formula <- as.formula(paste0(mediator, "~ ", exposure,  " + edad_s1 + sexo_s1 + edu_level_v00 + civil_status_v00 + area + 
              imc_v00 + exercise_day_v00 + smoke_status_v00 + alcoholg_v00 + qalcoholg_v00 + 
              depression_v00 + diab_prev_s1 + hta_v00 + colest_v00 + MED12score_v00"))       
    out_formula <- as.formula(paste0("zGCF_2yc ~ ", mediator, " + ", exposure, " + ", " zGCF_v00 + edad_s1 + sexo_s1 + grupo_int_v00 + edu_level_v00 + civil_status_v00 + area + 
              imc_v00 + exercise_day_v00 + smoke_status_v00 + alcoholg_v00 + qalcoholg_v00 + 
              depression_v00 + diab_prev_s1 + hta_v00 + colest_v00 + MED12score_v00"))             
    # Fit the mediation models      
    med_fit <- lm(med_formula, data = df)       
    out_fit <- lm(out_formula, data = df)             
    # Perform mediation analysis   
    med_out <- mediate(med_fit, out_fit, treat = exposure, mediator = mediator, sims = 1000)             
    # Store the results in the data frame       
    mediation_results_gcf <- rbind(mediation_results_gcf, data.frame(         
      exposure = exposure,         
      mediator = mediator,         
      ACME_estimate = summary(med_out)$d.avg,         
      ACME_lower95 = summary(med_out)$d.avg.ci[1],         
      ACME_upper95 = summary(med_out)$d.avg.ci[2],         
      ACME_p = summary(med_out)$d.avg.p,         
      ADE_estimate = summary(med_out)$z.avg,         
      ADE_lower95 = summary(med_out)$z.avg.ci[1],         
      ADE_upper95 = summary(med_out)$z.avg.ci[2],         
      ADE_p = summary(med_out)$z.avg.p,         
      total_estimate = summary(med_out)$tau.coef,         
      total_lower95 = summary(med_out)$tau.ci[1],         
      total_upper95 = summary(med_out)$tau.ci[2],         
      total_p = summary(med_out)$tau.p,         
      prop_estimate = summary(med_out)$n.avg,         
      prop_lower95 = summary(med_out)$n.avg.ci[1],         
      prop_upper95 = summary(med_out)$n.avg.ci[2],         
      prop_p = summary(med_out)$n.avg.p       
    ))     
  }   
} 

mediation_results_gcf$Outcome <- factor("Global cognitive function")


##General cognitive function---------------------------------------------------
# Create an empty data frame to store results 
mediation_results_gen <- data.frame(   
  exposure = character(),   
  mediator = character(),   
  ACME_estimate = numeric(),   ACME_lower95 = numeric(),   ACME_upper95 = numeric(),   ACME_p = numeric(),   
  ADE_estimate = numeric(),   ADE_lower95 = numeric(),   ADE_upper95 = numeric(),   ADE_p = numeric(),   
  total_estimate = numeric(),   total_lower95 = numeric(),   total_upper95 = numeric(),   total_p = numeric(),   
  prop_estimate = numeric(),   prop_lower95 = numeric(),   prop_upper95 = numeric(),   prop_p = numeric(),   stringsAsFactors = FALSE ) 


set.seed(2023) 
# Nested loops for each combination of exposure, mediator, and outcome 
for (exposure in exposures) {   
  for (mediator in mediators) {     
    med_formula <- as.formula(paste0(mediator, "~ ", exposure,  " + edad_s1 + sexo_s1 + edu_level_v00 + civil_status_v00 + area + 
              imc_v00 + exercise_day_v00 + smoke_status_v00 + alcoholg_v00 + qalcoholg_v00 + 
              depression_v00 + diab_prev_s1 + hta_v00 + colest_v00 + MED12score_v00"))       
    out_formula <- as.formula(paste0("zgenCF_2yc ~ ", mediator, " + ", exposure, " + ", " zgenCF_v00 + edad_s1 + sexo_s1 + grupo_int_v00 + edu_level_v00 + civil_status_v00 + area + 
              imc_v00 + exercise_day_v00 + smoke_status_v00 + alcoholg_v00 + qalcoholg_v00 + 
              depression_v00 + diab_prev_s1 + hta_v00 + colest_v00 + MED12score_v00"))             
    # Fit the mediation models    
    med_fit <- lm(med_formula, data = df)       
    out_fit <- lm(out_formula, data = df)             
    # Perform mediation analysis  
    med_out <- mediate(med_fit, out_fit, treat = exposure, mediator = mediator, sims = 1000)             
    # Store the results in the data frame       
    mediation_results_gen <- rbind(mediation_results_gen, data.frame(         
      exposure = exposure,         
      mediator = mediator,         
      ACME_estimate = summary(med_out)$d.avg,         
      ACME_lower95 = summary(med_out)$d.avg.ci[1],         
      ACME_upper95 = summary(med_out)$d.avg.ci[2],         
      ACME_p = summary(med_out)$d.avg.p,         
      ADE_estimate = summary(med_out)$z.avg,         
      ADE_lower95 = summary(med_out)$z.avg.ci[1],         
      ADE_upper95 = summary(med_out)$z.avg.ci[2],         
      ADE_p = summary(med_out)$z.avg.p,         
      total_estimate = summary(med_out)$tau.coef,         
      total_lower95 = summary(med_out)$tau.ci[1],         
      total_upper95 = summary(med_out)$tau.ci[2],         
      total_p = summary(med_out)$tau.p,         
      prop_estimate = summary(med_out)$n.avg,         
      prop_lower95 = summary(med_out)$n.avg.ci[1],         
      prop_upper95 = summary(med_out)$n.avg.ci[2],         
      prop_p = summary(med_out)$n.avg.p       
    ))     
  }   
} 

mediation_results_gen$Outcome <- factor("General cognitive function")

##Executive function---------------------------------------------------
# Create an empty data frame to store results 
mediation_results_ex <- data.frame(   
  exposure = character(),   
  mediator = character(),   
  ACME_estimate = numeric(),   ACME_lower95 = numeric(),   ACME_upper95 = numeric(),   ACME_p = numeric(),   
  ADE_estimate = numeric(),   ADE_lower95 = numeric(),   ADE_upper95 = numeric(),   ADE_p = numeric(),   
  total_estimate = numeric(),   total_lower95 = numeric(),   total_upper95 = numeric(),   total_p = numeric(),   
  prop_estimate = numeric(),   prop_lower95 = numeric(),   prop_upper95 = numeric(),   prop_p = numeric(),   stringsAsFactors = FALSE ) 


set.seed(2023) 
# Nested loops for each combination of exposure, mediator, and outcome 
for (exposure in exposures) {   
  for (mediator in mediators) {     
    med_formula <- as.formula(paste0(mediator, "~ ", exposure,  " + edad_s1 + sexo_s1 + edu_level_v00 + civil_status_v00 + area + 
              imc_v00 + exercise_day_v00 + smoke_status_v00 + alcoholg_v00 + qalcoholg_v00 + 
              depression_v00 + diab_prev_s1 + hta_v00 + colest_v00 + MED12score_v00"))       
    out_formula <- as.formula(paste0("zExF_2yc ~ ", mediator, " + ", exposure, " + ", "zExF_v00 + edad_s1 + sexo_s1 + grupo_int_v00 + edu_level_v00 + civil_status_v00 + area + 
              imc_v00 + exercise_day_v00 + smoke_status_v00 + alcoholg_v00 + qalcoholg_v00 + 
              depression_v00 + diab_prev_s1 + hta_v00 + colest_v00 + MED12score_v00"))             
    # Fit the mediation models 
    med_fit <- lm(med_formula, data = df)       
    out_fit <- lm(out_formula, data = df)             
    # Perform mediation analysis    
    med_out <- mediate(med_fit, out_fit, treat = exposure, mediator = mediator, sims = 1000)             
    # Store the results in the data frame       
    mediation_results_ex <- rbind(mediation_results_ex, data.frame(         
      exposure = exposure,         
      mediator = mediator,         
      ACME_estimate = summary(med_out)$d.avg,         
      ACME_lower95 = summary(med_out)$d.avg.ci[1],         
      ACME_upper95 = summary(med_out)$d.avg.ci[2],         
      ACME_p = summary(med_out)$d.avg.p,         
      ADE_estimate = summary(med_out)$z.avg,         
      ADE_lower95 = summary(med_out)$z.avg.ci[1],         
      ADE_upper95 = summary(med_out)$z.avg.ci[2],         
      ADE_p = summary(med_out)$z.avg.p,         
      total_estimate = summary(med_out)$tau.coef,         
      total_lower95 = summary(med_out)$tau.ci[1],         
      total_upper95 = summary(med_out)$tau.ci[2],         
      total_p = summary(med_out)$tau.p,         
      prop_estimate = summary(med_out)$n.avg,         
      prop_lower95 = summary(med_out)$n.avg.ci[1],         
      prop_upper95 = summary(med_out)$n.avg.ci[2],         
      prop_p = summary(med_out)$n.avg.p       
    ))     
  }   
} 

mediation_results_ex$Outcome <- factor("Executive function")

##Attention---------------------------------------------------
# Create an empty data frame to store results 
mediation_results_att <- data.frame(   
  exposure = character(),   
  mediator = character(),   
  ACME_estimate = numeric(),   ACME_lower95 = numeric(),   ACME_upper95 = numeric(),   ACME_p = numeric(),   
  ADE_estimate = numeric(),   ADE_lower95 = numeric(),   ADE_upper95 = numeric(),   ADE_p = numeric(),   
  total_estimate = numeric(),   total_lower95 = numeric(),   total_upper95 = numeric(),   total_p = numeric(),   
  prop_estimate = numeric(),   prop_lower95 = numeric(),   prop_upper95 = numeric(),   prop_p = numeric(),   stringsAsFactors = FALSE ) 


set.seed(2023) 
# Nested loops for each combination of exposure, mediator, and outcome 
for (exposure in exposures) {   
  for (mediator in mediators) {     
    med_formula <- as.formula(paste0(mediator, "~ ", exposure,  " + edad_s1 + sexo_s1 + edu_level_v00 + civil_status_v00 + area + 
              imc_v00 + exercise_day_v00 + smoke_status_v00 + alcoholg_v00 + qalcoholg_v00 + 
              depression_v00 + diab_prev_s1 + hta_v00 + colest_v00 + MED12score_v00"))       
    out_formula <- as.formula(paste0("zattention_2yc ~ ", mediator, " + ", exposure, " + ", "zattention_v00 + edad_s1 + sexo_s1 + grupo_int_v00 + edu_level_v00 + civil_status_v00 + area + 
              imc_v00 + exercise_day_v00 + smoke_status_v00 + alcoholg_v00 + qalcoholg_v00 + 
              depression_v00 + diab_prev_s1 + hta_v00 + colest_v00 + MED12score_v00"))             
    # Fit the mediation models      
    med_fit <- lm(med_formula, data = df)       
    out_fit <- lm(out_formula, data = df)             
    # Perform mediation analysis      
    med_out <- mediate(med_fit, out_fit, treat = exposure, mediator = mediator, sims = 1000)             
    # Store the results in the data frame       
    mediation_results_att <- rbind(mediation_results_att, data.frame(         
      exposure = exposure,         
      mediator = mediator,         
      ACME_estimate = summary(med_out)$d.avg,         
      ACME_lower95 = summary(med_out)$d.avg.ci[1],         
      ACME_upper95 = summary(med_out)$d.avg.ci[2],         
      ACME_p = summary(med_out)$d.avg.p,         
      ADE_estimate = summary(med_out)$z.avg,         
      ADE_lower95 = summary(med_out)$z.avg.ci[1],         
      ADE_upper95 = summary(med_out)$z.avg.ci[2],         
      ADE_p = summary(med_out)$z.avg.p,         
      total_estimate = summary(med_out)$tau.coef,         
      total_lower95 = summary(med_out)$tau.ci[1],         
      total_upper95 = summary(med_out)$tau.ci[2],         
      total_p = summary(med_out)$tau.p,         
      prop_estimate = summary(med_out)$n.avg,         
      prop_lower95 = summary(med_out)$n.avg.ci[1],         
      prop_upper95 = summary(med_out)$n.avg.ci[2],         
      prop_p = summary(med_out)$n.avg.p       
    ))     
  }   
}
mediation_results_att$Outcome <- factor("Attention")
##Language---------------------------------------------------
# Create an empty data frame to store results 
mediation_results_lang <- data.frame(   
  exposure = character(),   
  mediator = character(),   
  ACME_estimate = numeric(),   ACME_lower95 = numeric(),   ACME_upper95 = numeric(),   ACME_p = numeric(),   
  ADE_estimate = numeric(),   ADE_lower95 = numeric(),   ADE_upper95 = numeric(),   ADE_p = numeric(),   
  total_estimate = numeric(),   total_lower95 = numeric(),   total_upper95 = numeric(),   total_p = numeric(),   
  prop_estimate = numeric(),   prop_lower95 = numeric(),   prop_upper95 = numeric(),   prop_p = numeric(),   stringsAsFactors = FALSE ) 


set.seed(2023) 
# Nested loops for each combination of exposure, mediator, and outcome 
for (exposure in exposures) {   
  for (mediator in mediators) {     
    med_formula <- as.formula(paste0(mediator, "~ ", exposure,  " + edad_s1 + sexo_s1 + edu_level_v00 + civil_status_v00 + area + 
              imc_v00 + exercise_day_v00 + smoke_status_v00 + alcoholg_v00 + qalcoholg_v00 + 
              depression_v00 + diab_prev_s1 + hta_v00 + colest_v00 + MED12score_v00"))       
    out_formula <- as.formula(paste0("zlanguage_2yc ~ ", mediator, " + ", exposure, " + ", "zlanguage_v00 + edad_s1 + sexo_s1 + grupo_int_v00 + edu_level_v00 + civil_status_v00 + area + 
              imc_v00 + exercise_day_v00 + smoke_status_v00 + alcoholg_v00 + qalcoholg_v00 + 
              depression_v00 + diab_prev_s1 + hta_v00 + colest_v00 + MED12score_v00"))             
    # Fit the mediation models 
    med_fit <- lm(med_formula, data = df)       
    out_fit <- lm(out_formula, data = df)             
    # Perform mediation analysis     
    med_out <- mediate(med_fit, out_fit, treat = exposure, mediator = mediator, sims = 1000)             
    # Store the results in the data frame       
    mediation_results_lang <- rbind(mediation_results_lang, data.frame(         
      exposure = exposure,         
      mediator = mediator,         
      ACME_estimate = summary(med_out)$d.avg,         
      ACME_lower95 = summary(med_out)$d.avg.ci[1],         
      ACME_upper95 = summary(med_out)$d.avg.ci[2],         
      ACME_p = summary(med_out)$d.avg.p,         
      ADE_estimate = summary(med_out)$z.avg,         
      ADE_lower95 = summary(med_out)$z.avg.ci[1],         
      ADE_upper95 = summary(med_out)$z.avg.ci[2],         
      ADE_p = summary(med_out)$z.avg.p,         
      total_estimate = summary(med_out)$tau.coef,         
      total_lower95 = summary(med_out)$tau.ci[1],         
      total_upper95 = summary(med_out)$tau.ci[2],         
      total_p = summary(med_out)$tau.p,         
      prop_estimate = summary(med_out)$n.avg,         
      prop_lower95 = summary(med_out)$n.avg.ci[1],         
      prop_upper95 = summary(med_out)$n.avg.ci[2],         
      prop_p = summary(med_out)$n.avg.p       
    ))     
  }   
} 

mediation_results_lang$Outcome <- factor("Language")

mediation_results <- rbind(mediation_results_gcf, mediation_results_gen, mediation_results_ex, mediation_results_att, mediation_results_lang)
mediation_results_sig <- mediation_results[which((mediation_results$ACME_lower95 >= 0 & mediation_results$ACME_upper95 >= 0) |  mediation_results$ACME_lower95 <= 0 & mediation_results$ACME_upper95 <= 0),]


#################################################  OLIVE OIL-RELATED TAXA ###################################################

# Data preparation (bacterial genus)---------------------------------------------------
# * load data ---------------------------------------------------
ps <- readRDS("/home/jiaqi/Jiaqi_Ni/OliveMicrobiomeCognition_Microbiome/DATA/jiaqi_Phyloseq_110324.rds") #baseline microbiota N=656
tse <- makeTreeSummarizedExperimentFromPhyloseq(ps) # convert phyloseq to TSE
tse
# * filter taxa ---------------------------------------------------
# Genus level
tse_genus <- mergeFeaturesByRank(tse, rank = "Genus",
                                 agglomerateTree = TRUE) #agglomerated OTU to genus 
tse_genus <- transformAssay(tse_genus, assay.type = "counts", method = "relabundance") #add relative abundance assay
tse_genus <- subsetByPrevalentFeatures(tse_genus,
                                       prevalence = 10 / 100,
                                       detection = 1/1000, as_relative = TRUE, include.lowest = TRUE) #include only features with total relative abundance >= 0.001 in at least 10% of the samples
tse_genus <- transformAssay(x = tse_genus, assay.type = "relabundance", method = "clr", 
                            pseudocount = 1, name = "clr") #add CLR transformed relative abundances assay 
taxa_table_g <- as.data.frame(t(assays(tse_genus)$clr))
colnames(taxa_table_g) <- sub("^Genus:", "", colnames(taxa_table_g))

taxa_sig <- c("Acidaminococcus", "Adlercreutzia", "Akkermansia" ,"Bacteroides", "Blautia", "CAG-56", "Clostridium_sensu_stricto_1",   "Collinsella"                  
              , "Dorea"   ,                      "Eubacterium_hallii_group"     
              , "Mogibacterium"  ,               "Phascolarctobacterium"        
              , "Romboutsia" ,                   "Streptococcus"                
              , "___5"     ,                     "uncultured_8"                 
              , "Senegalimassilia"  ,            "Christensenellaceae_R-7_group"
              , "Faecalibacterium")                

taxa_table_g_sig<- taxa_table_g[, taxa_sig]
columns<- taxa_table_g_sig[, sapply(taxa_table_g_sig, is.numeric)] #extract the columns to scale 
taxa_table_g_sig_s<-as.data.frame(scale(columns)) #scale taxa 
taxa_table_g_sig_s$SampleIDcounts <- rownames(taxa_table_g_sig_s)
# Meta data
meta_table <- as.data.frame(colData(tse))
df <-  merge(meta_table, taxa_table_g_sig_s,  by = "SampleIDcounts")
# Define and fix continuous and categorical variables
catVars <- c("sexo_s1",  "nodo", "area", "grupo_int_v00", "edu_level_v00", "civil_status_v00",
             "ter_ajust_olivatot_v00", "ter_ajust_ac_olivavir_v00", "ter_ajust_re_oliva_v00", 
             "diab_prev_s1", "hta_v00", "colest_v00", "depression_v00",  "antihypertensive_v00", "antidepressant_v00", 
             "antidiabetic_v00", "hypolipidemic_v00", "ter_exercise_day_v00", 
             "smoke_status_v00","ter_MED12score_v00","ter_alcoholg_v00")
df[,catVars] <- lapply(df[,catVars], as.factor)

conVars <- c("edad_s1", "olivatot_v00",  "ac_olivavir_v00", "re_oliva_v00", "ajust_tspolivatot_v00", "ajust_tspac_olivavir_v00", "ajust_tspre_oliva_v00", 
             "ajust_olivatot_v00", "ajust_ac_olivavir_v00", "ajust_re_oliva_v00", 
             "ptrend_ter_ajust_olivatot_v00", "ptrend_ter_ajust_ac_olivavir_v00", "ptrend_ter_ajust_re_oliva_v00",
             "imc_v00", "waist_cir_v00", "exercise_day_v00", "alcoholg_v00", "qalcoholg_v00", "MED12score_v00",
             "Acidaminococcus", "Adlercreutzia", "Akkermansia" ,"Bacteroides", "Blautia", "CAG-56", "Clostridium_sensu_stricto_1",   "Collinsella"                  
             , "Dorea"   ,                      "Eubacterium_hallii_group"     
             , "Mogibacterium"  ,               "Phascolarctobacterium"        
             , "Romboutsia" ,                   "Streptococcus"                
             , "___5"     ,                     "uncultured_8"                 
             , "Senegalimassilia"  ,            "Christensenellaceae_R-7_group"
             , "Faecalibacterium")
df[,conVars] <- lapply(df[,conVars], as.numeric)

colnames(df)[which(colnames(df)=="CAG-56")] <- "CAG_56"
colnames(df)[which(colnames(df)=="Christensenellaceae_R-7_group")] <- "Christensenellaceae_R_7_group"
colnames(df)[which(colnames(df)=="___5")] <- "u___5"
# Example vectors of exposure, mediator, and outcome names 
exposures <- c("ajust_tspolivatot_v00", "ajust_tspac_olivavir_v00", "ajust_tspre_oliva_v00") 
mediators <- c("Acidaminococcus", "Adlercreutzia", "Akkermansia" ,"Bacteroides", "Blautia", "CAG_56", "Clostridium_sensu_stricto_1",   "Collinsella"                  
               , "Dorea"   ,                      "Eubacterium_hallii_group"     
               , "Mogibacterium"  ,               "Phascolarctobacterium"        
               , "Romboutsia" ,                   "Streptococcus"                
               , "u___5"     ,                     "uncultured_8"                 
               , "Senegalimassilia"  ,            "Christensenellaceae_R_7_group"
               , "Faecalibacterium") 

# Fully adjusted model ---------------------------------------------------
# baseline age, sex, baseline cognitive test score, intervention group, geographical area, educational level, civil status
#BMI, physical activity, smoking status, alcohol consumption, depressive symptomatology, diabetes prevalence, hypertension prevalence, hypercholesterolemia prevalence, 
#adherence to Mediterranean diet (12 scores)

##GCF---------------------------------------------------
# Create an empty data frame to store results 
mediation_results_gcf <- data.frame(   
  exposure = character(),   
  mediator = character(),   
  ACME_estimate = numeric(),   ACME_lower95 = numeric(),   ACME_upper95 = numeric(),   ACME_p = numeric(),   
  ADE_estimate = numeric(),   ADE_lower95 = numeric(),   ADE_upper95 = numeric(),   ADE_p = numeric(),   
  total_estimate = numeric(),   total_lower95 = numeric(),   total_upper95 = numeric(),   total_p = numeric(),   
  prop_estimate = numeric(),   prop_lower95 = numeric(),   prop_upper95 = numeric(),   prop_p = numeric(),   stringsAsFactors = FALSE ) 


set.seed(2023) 
# Nested loops for each combination of exposure, mediator, and outcome 
for (exposure in exposures) {   
  for (mediator in mediators) {     
    med_formula <- as.formula(paste0(mediator, "~ ", exposure,  " + edad_s1 + sexo_s1 + edu_level_v00 + civil_status_v00 + area + 
              imc_v00 + exercise_day_v00 + smoke_status_v00 + alcoholg_v00 + qalcoholg_v00 + 
              depression_v00 + diab_prev_s1 + hta_v00 + colest_v00 + MED12score_v00"))       
    out_formula <- as.formula(paste0("zGCF_2yc ~ ", mediator, " + ", exposure, " + ", " zGCF_v00 + edad_s1 + sexo_s1 + grupo_int_v00 + edu_level_v00 + civil_status_v00 + area + 
              imc_v00 + exercise_day_v00 + smoke_status_v00 + alcoholg_v00 + qalcoholg_v00 + 
              depression_v00 + diab_prev_s1 + hta_v00 + colest_v00 + MED12score_v00"))             
    # Fit the mediation models   
    med_fit <- lm(med_formula, data = df)       
    out_fit <- lm(out_formula, data = df)             
    # Perform mediation analysis 
    med_out <- mediate(med_fit, out_fit, treat = exposure, mediator = mediator, sims = 1000, boot = T)             
    # Store the results in the data frame       
    mediation_results_gcf <- rbind(mediation_results_gcf, data.frame(         
      exposure = exposure,         
      mediator = mediator,         
      ACME_estimate = summary(med_out)$d.avg,         
      ACME_lower95 = summary(med_out)$d.avg.ci[1],         
      ACME_upper95 = summary(med_out)$d.avg.ci[2],         
      ACME_p = summary(med_out)$d.avg.p,         
      ADE_estimate = summary(med_out)$z.avg,         
      ADE_lower95 = summary(med_out)$z.avg.ci[1],         
      ADE_upper95 = summary(med_out)$z.avg.ci[2],         
      ADE_p = summary(med_out)$z.avg.p,         
      total_estimate = summary(med_out)$tau.coef,         
      total_lower95 = summary(med_out)$tau.ci[1],         
      total_upper95 = summary(med_out)$tau.ci[2],         
      total_p = summary(med_out)$tau.p,         
      prop_estimate = summary(med_out)$n.avg,         
      prop_lower95 = summary(med_out)$n.avg.ci[1],         
      prop_upper95 = summary(med_out)$n.avg.ci[2],         
      prop_p = summary(med_out)$n.avg.p       
    ))     
  }   
} 

mediation_results_gcf$Outcome <- factor("Global cognitive function")


##General cognitive function---------------------------------------------------
# Create an empty data frame to store results 
mediation_results_gen <- data.frame(   
  exposure = character(),   
  mediator = character(),   
  ACME_estimate = numeric(),   ACME_lower95 = numeric(),   ACME_upper95 = numeric(),   ACME_p = numeric(),   
  ADE_estimate = numeric(),   ADE_lower95 = numeric(),   ADE_upper95 = numeric(),   ADE_p = numeric(),   
  total_estimate = numeric(),   total_lower95 = numeric(),   total_upper95 = numeric(),   total_p = numeric(),   
  prop_estimate = numeric(),   prop_lower95 = numeric(),   prop_upper95 = numeric(),   prop_p = numeric(),   stringsAsFactors = FALSE ) 


set.seed(2023) 
# Nested loops for each combination of exposure, mediator, and outcome 
for (exposure in exposures) {   
  for (mediator in mediators) {     
    med_formula <- as.formula(paste0(mediator, "~ ", exposure,  " + edad_s1 + sexo_s1 + edu_level_v00 + civil_status_v00 + area + 
              imc_v00 + exercise_day_v00 + smoke_status_v00 + alcoholg_v00 + qalcoholg_v00 + 
              depression_v00 + diab_prev_s1 + hta_v00 + colest_v00 + MED12score_v00"))       
    out_formula <- as.formula(paste0("zgenCF_2yc ~ ", mediator, " + ", exposure, " + ", " zgenCF_v00 + edad_s1 + sexo_s1 + grupo_int_v00 + edu_level_v00 + civil_status_v00 + area + 
              imc_v00 + exercise_day_v00 + smoke_status_v00 + alcoholg_v00 + qalcoholg_v00 + 
              depression_v00 + diab_prev_s1 + hta_v00 + colest_v00 + MED12score_v00"))             
    # Fit the mediation models    
    med_fit <- lm(med_formula, data = df)       
    out_fit <- lm(out_formula, data = df)             
    # Perform mediation analysis    
    med_out <- mediate(med_fit, out_fit, treat = exposure, mediator = mediator, sims = 1000, boot = T)             
    # Store the results in the data frame       
    mediation_results_gen <- rbind(mediation_results_gen, data.frame(         
      exposure = exposure,         
      mediator = mediator,         
      ACME_estimate = summary(med_out)$d.avg,         
      ACME_lower95 = summary(med_out)$d.avg.ci[1],         
      ACME_upper95 = summary(med_out)$d.avg.ci[2],         
      ACME_p = summary(med_out)$d.avg.p,         
      ADE_estimate = summary(med_out)$z.avg,         
      ADE_lower95 = summary(med_out)$z.avg.ci[1],         
      ADE_upper95 = summary(med_out)$z.avg.ci[2],         
      ADE_p = summary(med_out)$z.avg.p,         
      total_estimate = summary(med_out)$tau.coef,         
      total_lower95 = summary(med_out)$tau.ci[1],         
      total_upper95 = summary(med_out)$tau.ci[2],         
      total_p = summary(med_out)$tau.p,         
      prop_estimate = summary(med_out)$n.avg,         
      prop_lower95 = summary(med_out)$n.avg.ci[1],         
      prop_upper95 = summary(med_out)$n.avg.ci[2],         
      prop_p = summary(med_out)$n.avg.p       
    ))     
  }   
} 

mediation_results_gen$Outcome <- factor("General cognitive function")

##Executive function---------------------------------------------------
# Create an empty data frame to store results 
mediation_results_ex <- data.frame(   
  exposure = character(),   
  mediator = character(),   
  ACME_estimate = numeric(),   ACME_lower95 = numeric(),   ACME_upper95 = numeric(),   ACME_p = numeric(),   
  ADE_estimate = numeric(),   ADE_lower95 = numeric(),   ADE_upper95 = numeric(),   ADE_p = numeric(),   
  total_estimate = numeric(),   total_lower95 = numeric(),   total_upper95 = numeric(),   total_p = numeric(),   
  prop_estimate = numeric(),   prop_lower95 = numeric(),   prop_upper95 = numeric(),   prop_p = numeric(),   stringsAsFactors = FALSE ) 


set.seed(2023) 
# Nested loops for each combination of exposure, mediator, and outcome 
for (exposure in exposures) {   
  for (mediator in mediators) {     
    med_formula <- as.formula(paste0(mediator, "~ ", exposure,  " + edad_s1 + sexo_s1 + edu_level_v00 + civil_status_v00 + area + 
              imc_v00 + exercise_day_v00 + smoke_status_v00 + alcoholg_v00 + qalcoholg_v00 + 
              depression_v00 + diab_prev_s1 + hta_v00 + colest_v00 + MED12score_v00"))       
    out_formula <- as.formula(paste0("zExF_2yc ~ ", mediator, " + ", exposure, " + ", "zExF_v00 + edad_s1 + sexo_s1 + grupo_int_v00 + edu_level_v00 + civil_status_v00 + area + 
              imc_v00 + exercise_day_v00 + smoke_status_v00 + alcoholg_v00 + qalcoholg_v00 + 
              depression_v00 + diab_prev_s1 + hta_v00 + colest_v00 + MED12score_v00"))             
    # Fit the mediation models    
    med_fit <- lm(med_formula, data = df)       
    out_fit <- lm(out_formula, data = df)             
    # Perform mediation analysis        
    med_out <- mediate(med_fit, out_fit, treat = exposure, mediator = mediator, sims = 1000, boot = T)             
    # Store the results in the data frame       
    mediation_results_ex <- rbind(mediation_results_ex, data.frame(         
      exposure = exposure,         
      mediator = mediator,         
      ACME_estimate = summary(med_out)$d.avg,         
      ACME_lower95 = summary(med_out)$d.avg.ci[1],         
      ACME_upper95 = summary(med_out)$d.avg.ci[2],         
      ACME_p = summary(med_out)$d.avg.p,         
      ADE_estimate = summary(med_out)$z.avg,         
      ADE_lower95 = summary(med_out)$z.avg.ci[1],         
      ADE_upper95 = summary(med_out)$z.avg.ci[2],         
      ADE_p = summary(med_out)$z.avg.p,         
      total_estimate = summary(med_out)$tau.coef,         
      total_lower95 = summary(med_out)$tau.ci[1],         
      total_upper95 = summary(med_out)$tau.ci[2],         
      total_p = summary(med_out)$tau.p,         
      prop_estimate = summary(med_out)$n.avg,         
      prop_lower95 = summary(med_out)$n.avg.ci[1],         
      prop_upper95 = summary(med_out)$n.avg.ci[2],         
      prop_p = summary(med_out)$n.avg.p       
    ))     
  }   
} 

mediation_results_ex$Outcome <- factor("Executive function")

##Attention---------------------------------------------------
# Create an empty data frame to store results 
mediation_results_att <- data.frame(   
  exposure = character(),   
  mediator = character(),   
  ACME_estimate = numeric(),   ACME_lower95 = numeric(),   ACME_upper95 = numeric(),   ACME_p = numeric(),   
  ADE_estimate = numeric(),   ADE_lower95 = numeric(),   ADE_upper95 = numeric(),   ADE_p = numeric(),   
  total_estimate = numeric(),   total_lower95 = numeric(),   total_upper95 = numeric(),   total_p = numeric(),   
  prop_estimate = numeric(),   prop_lower95 = numeric(),   prop_upper95 = numeric(),   prop_p = numeric(),   stringsAsFactors = FALSE ) 


set.seed(2023) 
# Nested loops for each combination of exposure, mediator, and outcome 
for (exposure in exposures) {   
  for (mediator in mediators) {     
    med_formula <- as.formula(paste0(mediator, "~ ", exposure,  " + edad_s1 + sexo_s1 + edu_level_v00 + civil_status_v00 + area + 
              imc_v00 + exercise_day_v00 + smoke_status_v00 + alcoholg_v00 + qalcoholg_v00 + 
              depression_v00 + diab_prev_s1 + hta_v00 + colest_v00 + MED12score_v00"))       
    out_formula <- as.formula(paste0("zattention_2yc ~ ", mediator, " + ", exposure, " + ", "zattention_v00 + edad_s1 + sexo_s1 + grupo_int_v00 + edu_level_v00 + civil_status_v00 + area + 
              imc_v00 + exercise_day_v00 + smoke_status_v00 + alcoholg_v00 + qalcoholg_v00 + 
              depression_v00 + diab_prev_s1 + hta_v00 + colest_v00 + MED12score_v00"))             
    # Fit the mediation models      
    med_fit <- lm(med_formula, data = df)       
    out_fit <- lm(out_formula, data = df)             
    # Perform mediation analysis       
    med_out <- mediate(med_fit, out_fit, treat = exposure, mediator = mediator, sims = 1000, boot = T)             
    # Store the results in the data frame       
    mediation_results_att <- rbind(mediation_results_att, data.frame(         
      exposure = exposure,         
      mediator = mediator,         
      ACME_estimate = summary(med_out)$d.avg,         
      ACME_lower95 = summary(med_out)$d.avg.ci[1],         
      ACME_upper95 = summary(med_out)$d.avg.ci[2],         
      ACME_p = summary(med_out)$d.avg.p,         
      ADE_estimate = summary(med_out)$z.avg,         
      ADE_lower95 = summary(med_out)$z.avg.ci[1],         
      ADE_upper95 = summary(med_out)$z.avg.ci[2],         
      ADE_p = summary(med_out)$z.avg.p,         
      total_estimate = summary(med_out)$tau.coef,         
      total_lower95 = summary(med_out)$tau.ci[1],         
      total_upper95 = summary(med_out)$tau.ci[2],         
      total_p = summary(med_out)$tau.p,         
      prop_estimate = summary(med_out)$n.avg,         
      prop_lower95 = summary(med_out)$n.avg.ci[1],         
      prop_upper95 = summary(med_out)$n.avg.ci[2],         
      prop_p = summary(med_out)$n.avg.p       
    ))     
  }   
}
mediation_results_att$Outcome <- factor("Attention")

##Language---------------------------------------------------
# Create an empty data frame to store results 
mediation_results_lang <- data.frame(   
  exposure = character(),   
  mediator = character(),   
  ACME_estimate = numeric(),   ACME_lower95 = numeric(),   ACME_upper95 = numeric(),   ACME_p = numeric(),   
  ADE_estimate = numeric(),   ADE_lower95 = numeric(),   ADE_upper95 = numeric(),   ADE_p = numeric(),   
  total_estimate = numeric(),   total_lower95 = numeric(),   total_upper95 = numeric(),   total_p = numeric(),   
  prop_estimate = numeric(),   prop_lower95 = numeric(),   prop_upper95 = numeric(),   prop_p = numeric(),   stringsAsFactors = FALSE ) 


set.seed(2023) 
# Nested loops for each combination of exposure, mediator, and outcome 
for (exposure in exposures) {   
  for (mediator in mediators) {     
    med_formula <- as.formula(paste0(mediator, "~ ", exposure,  " + edad_s1 + sexo_s1 + edu_level_v00 + civil_status_v00 + area + 
              imc_v00 + exercise_day_v00 + smoke_status_v00 + alcoholg_v00 + qalcoholg_v00 + 
              depression_v00 + diab_prev_s1 + hta_v00 + colest_v00 + MED12score_v00"))       
    out_formula <- as.formula(paste0("zlanguage_2yc ~ ", mediator, " + ", exposure, " + ", "zlanguage_v00 + edad_s1 + sexo_s1 + grupo_int_v00 + edu_level_v00 + civil_status_v00 + area + 
              imc_v00 + exercise_day_v00 + smoke_status_v00 + alcoholg_v00 + qalcoholg_v00 + 
              depression_v00 + diab_prev_s1 + hta_v00 + colest_v00 + MED12score_v00"))             
    # Fit the mediation models        
    med_fit <- lm(med_formula, data = df)       
    out_fit <- lm(out_formula, data = df)             
    # Perform mediation analysis   
    med_out <- mediate(med_fit, out_fit, treat = exposure, mediator = mediator, sims = 1000, boot = T)             
    # Store the results in the data frame       
    mediation_results_lang <- rbind(mediation_results_lang, data.frame(         
      exposure = exposure,         
      mediator = mediator,         
      ACME_estimate = summary(med_out)$d.avg,         
      ACME_lower95 = summary(med_out)$d.avg.ci[1],         
      ACME_upper95 = summary(med_out)$d.avg.ci[2],         
      ACME_p = summary(med_out)$d.avg.p,         
      ADE_estimate = summary(med_out)$z.avg,         
      ADE_lower95 = summary(med_out)$z.avg.ci[1],         
      ADE_upper95 = summary(med_out)$z.avg.ci[2],         
      ADE_p = summary(med_out)$z.avg.p,         
      total_estimate = summary(med_out)$tau.coef,         
      total_lower95 = summary(med_out)$tau.ci[1],         
      total_upper95 = summary(med_out)$tau.ci[2],         
      total_p = summary(med_out)$tau.p,         
      prop_estimate = summary(med_out)$n.avg,         
      prop_lower95 = summary(med_out)$n.avg.ci[1],         
      prop_upper95 = summary(med_out)$n.avg.ci[2],         
      prop_p = summary(med_out)$n.avg.p       
    ))     
  }   
} 

mediation_results_lang$Outcome <- factor("Language")

mediation_results_taxa <- rbind(mediation_results_gcf, mediation_results_gen, mediation_results_ex, mediation_results_att, mediation_results_lang)
mediation_results_taxa_sig <- mediation_results_taxa[which((mediation_results_taxa$ACME_lower95 >= 0 & mediation_results_taxa$ACME_upper95 >= 0) |  
                                                             mediation_results_taxa$ACME_lower95 <= 0 & mediation_results_taxa$ACME_upper95 <= 0),]

mediation_results_full_sig<-rbind(mediation_results_sig, mediation_results_taxa_sig)
mediation_results_full_sig_sig <- mediation_results_full_sig[which((mediation_results_full_sig$total_lower95 >= 0 & mediation_results_full_sig$total_upper95 >= 0) |  
                                                                     mediation_results_full_sig$total_lower95 <= 0 & mediation_results_full_sig$total_upper95 <= 0),]


#################################################  EXPOSURE-MEDIATOR INTERACTIONS ###################################################

# * TOTAL-Adlercreutzia-General cognitive function ---------------------------------------------------
yx.fit <- lm(data = df, formula= zgenCF_2yc ~ ajust_olivatot_v00 + zgenCF_v00 +  edad_s1 + sexo_s1 
             + area + edu_level_v00 + civil_status_v00 + grupo_int_v00 + imc_v00 + exercise_day_v00 + smoke_status_v00 
             + alcoholg_v00 + qalcoholg_v00 + depression_v00 + diab_prev_s1 + hta_v00 + colest_v00 + MED12score_v00 )
summary(yx.fit)
med.fit <-  lm(data = df, formula = Adlercreutzia ~ ajust_olivatot_v00 +  edad_s1 + sexo_s1 
               + area + edu_level_v00 + civil_status_v00 + imc_v00 + exercise_day_v00 + smoke_status_v00 
               + alcoholg_v00 + qalcoholg_v00 + depression_v00 + diab_prev_s1 + hta_v00 + colest_v00 +
                 MED12score_v00)
summary(med.fit)
out.fit <- lm(data = df, formula= zgenCF_2yc ~ Adlercreutzia * ajust_olivatot_v00 + zgenCF_v00 +  edad_s1 + sexo_s1 
              + area + edu_level_v00 + civil_status_v00 + grupo_int_v00 + imc_v00 + exercise_day_v00 + smoke_status_v00 
              + alcoholg_v00 + qalcoholg_v00 + depression_v00 + diab_prev_s1 + hta_v00 + colest_v00 + MED12score_v00 )

summary(out.fit)

set.seed(2023)
med.out <- mediate(med.fit, out.fit, treat = "ajust_olivatot_v00", mediator = "Adlercreutzia", sims=1000, boot = T)

test.TMint(med.out, conf.level = .95)


# * Extra virgin-Adlercreutzia-General cognitive function ---------------------------------------------------

yx.fit <- lm(data = df, formula= zgenCF_2yc ~ ajust_ac_olivavir_v00 + zgenCF_v00 +  edad_s1 + sexo_s1 
             + area + edu_level_v00 + civil_status_v00 + grupo_int_v00 + imc_v00 + exercise_day_v00 + smoke_status_v00 
             + alcoholg_v00 + qalcoholg_v00 + depression_v00 + diab_prev_s1 + hta_v00 + colest_v00 + MED12score_v00 )
summary(yx.fit)
med.fit <-  lm(data = df, formula = Adlercreutzia ~ ajust_ac_olivavir_v00 +  edad_s1 + sexo_s1 
               + area + edu_level_v00 + civil_status_v00 + imc_v00 + exercise_day_v00 + smoke_status_v00 
               + alcoholg_v00 + qalcoholg_v00 + depression_v00 + diab_prev_s1 + hta_v00 + colest_v00 +
                 MED12score_v00)
summary(med.fit)
out.fit <- lm(data = df, formula= zgenCF_2yc ~ Adlercreutzia * ajust_ac_olivavir_v00 + zgenCF_v00 +  edad_s1 + sexo_s1 
              + area + edu_level_v00 + civil_status_v00 + grupo_int_v00 + imc_v00 + exercise_day_v00 + smoke_status_v00 
              + alcoholg_v00 + qalcoholg_v00 + depression_v00 + diab_prev_s1 + hta_v00 + colest_v00 + MED12score_v00 )

summary(out.fit)
set.seed(2023)
med.out <- mediate(med.fit, out.fit, treat = "ajust_ac_olivavir_v00", mediator = "Adlercreutzia", sims=1000, boot = T)

test.TMint(med.out, conf.level = .95)



# * load data 
ps <- readRDS("/home/jiaqi/Jiaqi_Ni/OliveMicrobiomeCognition_Microbiome/DATA/jiaqi_Phyloseq_110324.rds") #baseline microbiota N=656
tse <- makeTreeSummarizedExperimentFromPhyloseq(ps) # convert phyloseq to TSE
tse
# * calculate alpha diversity  
#Richness
tse <- mia::estimateRichness(tse, 
                             index = c( "chao1"), 
                             name=c( "Chao1"))

#Evenness
tse <- estimateEvenness(tse, 
                        index = c("simpson"),
                        name = c("Simpson"))

#Diversity
tse <- estimateDiversity(tse, 
                         index = c("shannon", "inverse_simpson"),
                         name = c("Shannon", "InverseSimpson"))


# * calculate beta diversity  (Aitchison distance) 
##CLR transform
tse <- transformAssay(x = tse, assay.type = "counts", method = "clr", 
                      pseudocount = 1, name = "clr")

# Calculate Aitchison distance (euclidean distance of CLR-transformed species abundances)
beta <- vegan::vegdist(t(assays(tse)$clr), method = "euclidean")
# Take a PCA and choose the first 2 principal components
pca <- cmdscale(beta, k = 2, eig = TRUE)
##with alpha indexes
df <- as.data.frame(colData(tse))
df<-subset(df, select = -c(Chao1_se))

##with pc1 & pc2
pcadf <- as.data.frame(pca$points)
names(pcadf) <- c("PC1", "PC2")
pcadf$SampleIDcounts <- rownames(pcadf)

##final dataframe
df <- merge(df, pcadf, by = "SampleIDcounts")
which(colnames(df)=="Chao1")
df[, 169:ncol(df)] <- scale(df[, 169:ncol(df)]) #scale the continuous indices

# Define and fix continuous and categorical variables
catVars <- c("sexo_s1",  "nodo", "area", "grupo_int_v00", "edu_level_v00", "civil_status_v00",
             "ter_ajust_olivatot_v00", "ter_ajust_ac_olivavir_v00", "ter_ajust_re_oliva_v00", 
             "diab_prev_s1", "hta_v00", "colest_v00", "depression_v00",  "antihypertensive_v00", "antidepressant_v00", 
             "antidiabetic_v00", "hypolipidemic_v00", "ter_exercise_day_v00", 
             "smoke_status_v00","ter_MED12score_v00","ter_alcoholg_v00")
df[,catVars] <- lapply(df[,catVars], as.factor)

conVars <- c("edad_s1", "olivatot_v00",  "ac_olivavir_v00", "re_oliva_v00",  "ajust_tspolivatot_v00", "ajust_tspac_olivavir_v00", "ajust_tspre_oliva_v00",
             "ajust_olivatot_v00", "ajust_ac_olivavir_v00", "ajust_re_oliva_v00", 
             "ptrend_ter_ajust_olivatot_v00", "ptrend_ter_ajust_ac_olivavir_v00", "ptrend_ter_ajust_re_oliva_v00",
             "imc_v00", "waist_cir_v00", "exercise_day_v00", "alcoholg_v00", "qalcoholg_v00", "MED12score_v00",
             "Chao1","Simpson" ,"Shannon","InverseSimpson","PC1","PC2")
df[,conVars] <- lapply(df[,conVars], as.numeric)
# Example vectors of exposure, mediator, and outcome names 
exposures <- c("ajust_tspolivatot_v00", "ajust_tspac_olivavir_v00", "ajust_tspre_oliva_v00") 
mediators <- c("Chao1", "Simpson", "Shannon", "InverseSimpson", "PC1", "PC2") 

# * TOTAL-PC1-EXECUTIVE ---------------------------------------------------

yx.fit <- lm(data = df, formula= zExF_2yc ~ ajust_olivatot_v00 + zExF_v00 +  edad_s1 + sexo_s1 
             + area + edu_level_v00 + civil_status_v00 + grupo_int_v00 + imc_v00 + exercise_day_v00 + smoke_status_v00 
             + alcoholg_v00 + qalcoholg_v00 + depression_v00 + diab_prev_s1 + hta_v00 + colest_v00 + MED12score_v00 )
summary(yx.fit)
med.fit <-  lm(data = df, formula = PC1 ~ ajust_olivatot_v00 +  edad_s1 + sexo_s1 
               + area + edu_level_v00 + civil_status_v00 + imc_v00 + exercise_day_v00 + smoke_status_v00 
               + alcoholg_v00 + qalcoholg_v00 + depression_v00 + diab_prev_s1 + hta_v00 + colest_v00 +
                 MED12score_v00)
summary(med.fit)
out.fit <- lm(data = df, formula= zExF_2yc ~ PC1 * ajust_olivatot_v00 + zExF_v00 +  edad_s1 + sexo_s1 
              + area + edu_level_v00 + civil_status_v00 + grupo_int_v00 + imc_v00 + exercise_day_v00 + smoke_status_v00 
              + alcoholg_v00 + qalcoholg_v00 + depression_v00 + diab_prev_s1 + hta_v00 + colest_v00 + MED12score_v00 )

summary(out.fit)

set.seed(2023)
med.out <- mediate(med.fit, out.fit, treat = "ajust_olivatot_v00", mediator = "PC1", sims=1000)

test.TMint(med.out, conf.level = .95)

# * VIRGIN-PC1-EXECUTIVE ---------------------------------------------------

yx.fit <- lm(data = df, formula= zExF_2yc ~ ajust_ac_olivavir_v00 + zExF_v00 +  edad_s1 + sexo_s1 
             + area + edu_level_v00 + civil_status_v00 + grupo_int_v00 + imc_v00 + exercise_day_v00 + smoke_status_v00 
             + alcoholg_v00 + qalcoholg_v00 + depression_v00 + diab_prev_s1 + hta_v00 + colest_v00 + MED12score_v00 )
summary(yx.fit)
med.fit <-  lm(data = df, formula = PC1 ~ ajust_ac_olivavir_v00 +  edad_s1 + sexo_s1 
               + area + edu_level_v00 + civil_status_v00 + imc_v00 + exercise_day_v00 + smoke_status_v00 
               + alcoholg_v00 + qalcoholg_v00 + depression_v00 + diab_prev_s1 + hta_v00 + colest_v00 +
                 MED12score_v00)
summary(med.fit)
out.fit <- lm(data = df, formula= zExF_2yc ~ PC1 * ajust_ac_olivavir_v00 + zExF_v00 +  edad_s1 + sexo_s1 
              + area + edu_level_v00 + civil_status_v00 + grupo_int_v00 + imc_v00 + exercise_day_v00 + smoke_status_v00 
              + alcoholg_v00 + qalcoholg_v00 + depression_v00 + diab_prev_s1 + hta_v00 + colest_v00 + MED12score_v00 )

summary(out.fit)

set.seed(2023)
med.out <- mediate(med.fit, out.fit, treat = "ajust_ac_olivavir_v00", mediator = "PC1", sims=1000)

test.TMint(med.out, conf.level = .95)

# * VIRGIN-PC1-LANGUAGE ---------------------------------------------------

yx.fit <- lm(data = df, formula= zlanguage_2yc ~ ajust_ac_olivavir_v00 + zlanguage_v00 +  edad_s1 + sexo_s1 
             + area + edu_level_v00 + civil_status_v00 + grupo_int_v00 + imc_v00 + exercise_day_v00 + smoke_status_v00 
             + alcoholg_v00 + qalcoholg_v00 + depression_v00 + diab_prev_s1 + hta_v00 + colest_v00 + MED12score_v00 )
summary(yx.fit)
med.fit <-  lm(data = df, formula = PC1 ~ ajust_ac_olivavir_v00 +  edad_s1 + sexo_s1 
               + area + edu_level_v00 + civil_status_v00 + imc_v00 + exercise_day_v00 + smoke_status_v00 
               + alcoholg_v00 + qalcoholg_v00 + depression_v00 + diab_prev_s1 + hta_v00 + colest_v00 +
                 MED12score_v00)
summary(med.fit)
out.fit <- lm(data = df, formula= zlanguage_2yc ~ PC1 * ajust_ac_olivavir_v00 + zlanguage_v00 +  edad_s1 + sexo_s1 
              + area + edu_level_v00 + civil_status_v00 + grupo_int_v00 + imc_v00 + exercise_day_v00 + smoke_status_v00 
              + alcoholg_v00 + qalcoholg_v00 + depression_v00 + diab_prev_s1 + hta_v00 + colest_v00 + MED12score_v00 )

summary(out.fit)

set.seed(2023)
med.out <- mediate(med.fit, out.fit, treat = "ajust_ac_olivavir_v00", mediator = "PC1", sims=1000)

test.TMint(med.out, conf.level = .95)
