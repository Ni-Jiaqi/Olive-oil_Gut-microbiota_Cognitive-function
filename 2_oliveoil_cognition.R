# Title: olive oil consumption and changes in cognitive function 
# Author: Jiaqi Ni [ORCID](https://orcid.org/0000-0001-5625-0373)
# Description: This is the analysis script for assess prospective associations of olive oil consumption with changes in cognitive function in 2 years. 


# Set up working environment  ------------------------------------------------------------------
rm(list=ls())
getwd()
setwd("/home/jiaqi/Jiaqi_Ni/OliveMicrobiomeCognition_Microbiome")
options(max.print=1000000)
options(scipen = 999)

# Load dependencies -------------------------------------------------------
library(modelsummary)

# Load and prepare data ---------------------------------------------------
ps_df<- readRDS("DATA/metadata_FINAL_110324.rds") #load metadata
# Define and fix continuous and categorical variables
catVars <- c("sexo_s1",  "nodo", "area", "grupo_int_v00", "edu_level_v00", "civil_status_v00",
             "ter_ajust_olivatot_v00", "ter_ajust_ac_olivavir_v00", "ter_ajust_re_oliva_v00", 
             "diab_prev_s1", "hta_v00", "colest_v00", "depression_v00",  "antihypertensive_v00", "antidepressant_v00", 
             "antidiabetic_v00", "hypolipidemic_v00", "ter_exercise_day_v00", 
             "smoke_status_v00","ter_MED12score_v00","ter_alcoholg_v00")
ps_df[,catVars] <- lapply(ps_df[,catVars], as.factor)

conVars <- c("edad_s1", "olivatot_v00",  "ac_olivavir_v00", "re_oliva_v00",  
             "ajust_olivatot_v00", "ajust_ac_olivavir_v00", "ajust_re_oliva_v00", "ajust_tspolivatot_v00", "ajust_tspac_olivavir_v00", "ajust_tspre_oliva_v00",
             "ptrend_ter_ajust_olivatot_v00", "ptrend_ter_ajust_ac_olivavir_v00", "ptrend_ter_ajust_re_oliva_v00",
             "MMSE_v00", "MMSE_v02", "CDT_v00", "CDT_v02", "VFTa_v00", "VFTa_v02", "VFTp_v00", "VFTp_v02", 
             "TMTa_v00", "TMTa_v02", "TMTb_v00", "TMTb_v02", "DSTf_v00", "DSTf_v02", "DSTb_v00", "DSTb_v02", 
             "zMMSE_v00", "zMMSE_v02", "zCDT_v00", "zCDT_v02", "zVFTa_v00", "zVFTa_v02", "zVFTp_v00", "zVFTp_v02", 
             "zTMTa_v00", "zTMTa_v02", "zTMTb_v00", "zTMTb_v02", "zDSTf_v00", "zDSTf_v02", "zDSTb_v00", "zDSTb_v02", 
             "GCF_v00", "GCF_v02", "genCF_v00", "genCF_v02", "ExF_v00", "ExF_v02", "attention_v00", "attention_v02", 
             "language_v00", "language_v02", "zGCF_v00", "zGCF_v02", "zgenCF_v00", "zgenCF_v02", "zExF_v00", "zExF_v02", 
             "zattention_v00", "zattention_v02", "zlanguage_v00", "zlanguage_v02", 
             "zMMSE_2yc", "zCDT_2yc", "zVFTa_2yc", "zVFTp_2yc", "zTMTa_2yc", "zTMTb_2yc", "zDSTf_2yc", "zDSTb_2yc", 
             "zGCF_2yc", "zgenCF_2yc", "zExF_2yc", "zattention_2yc", "zlanguage_2yc", "MED12score_v00", 
             "imc_v00", "waist_cir_v00", "exercise_day_v00", "alcoholg_v00", "qalcoholg_v00",
             "hc_v00", "prot_v00", "gratot_v00", "mo_v00", "po_v00", "sa_v00",
             "porc_hc_v00", "porc_pr_v00", "porc_gr_v00", "porc_mo_v00", "porc_po_v00", "porc_sa_v00", 
             "fibra_v00", "col_v00", "fit_v00", "trans_v00", "linoleico_v00", "linolenico_v00", "omega3_v00", "n3marinos_v00", 
             "verdutot_v00", "frutatot_v00", "legumbre_v00", "cereal_v00", "lacteos_v00", "carnicos_v00", "pescados_v00", "fsecos_v00", "gallet_v00")
ps_df[,conVars] <- lapply(ps_df[,conVars], as.numeric)

################################################## MAIN ANALYSIS ######################################################
# Basic model ------------------------------------------------------------- 
# baseline age, sex, baseline cognitive test score
## Set variables for Loop

cognitionscore_baseline <-c("zMMSE_v00", "zCDT_v00", "zVFTa_v00", "zVFTp_v00", 
                            "zTMTa_v00", "zTMTb_v00", "zDSTf_v00", "zDSTb_v00", 
                            "zGCF_v00", "zgenCF_v00", "zExF_v00", "zattention_v00", "zlanguage_v00")

cognitionscore_change <-c("zMMSE_2yc", "zCDT_2yc", "zVFTa_2yc", "zVFTp_2yc", 
                          "zTMTa_2yc", "zTMTb_2yc", "zDSTf_2yc", "zDSTb_2yc", 
                          "zGCF_2yc", "zgenCF_2yc", "zExF_2yc", "zattention_2yc", "zlanguage_2yc")

### * total olive oil  (cat.) -------------------------------------------------------------
comon <- c("ter_ajust_olivatot_v00" ,   "edad_s1",  "sexo_s1")

changebasicresults <- list()

for (i in 1:length(cognitionscore_baseline)){
  baseline_var<-cognitionscore_baseline[i]
  change_var<-cognitionscore_change[i]
  formula <- as.formula(paste(change_var, "~", baseline_var, "+", paste(comon, collapse  = "+") ))
  changebasicmodels <- lm ( formula, data = ps_df)
  results_key<-paste(change_var, baseline_var, sep = "_")
  changebasicresults[[results_key]] <- summary(changebasicmodels)
  
}   

for (j in names(changebasicresults)){
  print(paste("Regression results for", j))
  print(changebasicresults[[j]])
  
}

modelsummary(changebasicresults, 
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "p.value")

### total olive oil  (p-trend)
comon <- c("ptrend_ter_ajust_olivatot_v00" ,   "edad_s1",  "sexo_s1")


changebasicresults <- list()

for (i in 1:length(cognitionscore_baseline)){
  baseline_var<-cognitionscore_baseline[i]
  change_var<-cognitionscore_change[i]
  formula <- as.formula(paste(change_var, "~", baseline_var, "+", paste(comon, collapse  = "+") ))
  changebasicmodels <- lm ( formula, data = ps_df)
  results_key<-paste(baseline_var, change_var, sep = "_")
  changebasicresults[[results_key]] <- summary(changebasicmodels)
  
}   

for (j in names(changebasicresults)){
  print(paste("Regression results for", j))
  print(changebasicresults[[j]])
  
}

modelsummary(changebasicresults, 
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "p.value")


### * total olive oil  (cont.) -------------------------------------------------------------
comon <- c("ajust_tspolivatot_v00" ,   "edad_s1",  "sexo_s1")

changebasicresults <- list()

for (i in 1:length(cognitionscore_baseline)){
  baseline_var<-cognitionscore_baseline[i]
  change_var<-cognitionscore_change[i]
  formula <- as.formula(paste(change_var, "~", baseline_var, "+", paste(comon, collapse  = "+") ))
  changebasicmodels <- lm ( formula, data = ps_df)
  results_key<-paste(baseline_var, change_var, sep = "_")
  changebasicresults[[results_key]] <- summary(changebasicmodels)
  
}   

for (j in names(changebasicresults)){
  print(paste("Regression results for", j))
  print(changebasicresults[[j]])
  
}
modelsummary(changebasicresults, 
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "p.value")

### * extra virgin olive oil  (cat.) -------------------------------------------------------------
comon <- c("ter_ajust_ac_olivavir_v00" ,   "edad_s1",  "sexo_s1")

changebasicresults <- list()

for (i in 1:length(cognitionscore_baseline)){
  baseline_var<-cognitionscore_baseline[i]
  change_var<-cognitionscore_change[i]
  formula <- as.formula(paste(change_var, "~", baseline_var, "+", paste(comon, collapse  = "+") ))
  changebasicmodels <- lm ( formula, data = ps_df)
  results_key<-paste(change_var, baseline_var, sep = "_")
  changebasicresults[[results_key]] <- summary(changebasicmodels)
  
}   

for (j in names(changebasicresults)){
  print(paste("Regression results for", j))
  print(changebasicresults[[j]])
  
}

modelsummary(changebasicresults, 
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "p.value")

### extra virgin olive oil  (p-trend)
comon <- c("ptrend_ter_ajust_ac_olivavir_v00" ,   "edad_s1",  "sexo_s1")


changebasicresults <- list()

for (i in 1:length(cognitionscore_baseline)){
  baseline_var<-cognitionscore_baseline[i]
  change_var<-cognitionscore_change[i]
  formula <- as.formula(paste(change_var, "~", baseline_var, "+", paste(comon, collapse  = "+") ))
  changebasicmodels <- lm ( formula, data = ps_df)
  results_key<-paste(baseline_var, change_var, sep = "_")
  changebasicresults[[results_key]] <- summary(changebasicmodels)
  
}   

for (j in names(changebasicresults)){
  print(paste("Regression results for", j))
  print(changebasicresults[[j]])
  
}

modelsummary(changebasicresults, 
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "p.value")


### * extra virgin olive oil  (cont.) -------------------------------------------------------------
comon <- c("ajust_tspac_olivavir_v00" ,   "edad_s1",  "sexo_s1")

changebasicresults <- list()

for (i in 1:length(cognitionscore_baseline)){
  baseline_var<-cognitionscore_baseline[i]
  change_var<-cognitionscore_change[i]
  formula <- as.formula(paste(change_var, "~", baseline_var, "+", paste(comon, collapse  = "+") ))
  changebasicmodels <- lm ( formula, data = ps_df)
  results_key<-paste(baseline_var, change_var, sep = "_")
  changebasicresults[[results_key]] <- summary(changebasicmodels)
  
}   

for (j in names(changebasicresults)){
  print(paste("Regression results for", j))
  print(changebasicresults[[j]])
  
}
modelsummary(changebasicresults, 
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "p.value")

### * refined olive oil  (cat.) -------------------------------------------------------------
comon <- c("ter_ajust_re_oliva_v00" ,   "edad_s1",  "sexo_s1")

changebasicresults <- list()

for (i in 1:length(cognitionscore_baseline)){
  baseline_var<-cognitionscore_baseline[i]
  change_var<-cognitionscore_change[i]
  formula <- as.formula(paste(change_var, "~", baseline_var, "+", paste(comon, collapse  = "+") ))
  changebasicmodels <- lm ( formula, data = ps_df)
  results_key<-paste(change_var, baseline_var, sep = "_")
  changebasicresults[[results_key]] <- summary(changebasicmodels)
  
}   

for (j in names(changebasicresults)){
  print(paste("Regression results for", j))
  print(changebasicresults[[j]])
  
}

modelsummary(changebasicresults, 
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "p.value")

### refined olive oil  (p-trend)
comon <- c("ptrend_ter_ajust_re_oliva_v00" ,   "edad_s1",  "sexo_s1")


changebasicresults <- list()

for (i in 1:length(cognitionscore_baseline)){
  baseline_var<-cognitionscore_baseline[i]
  change_var<-cognitionscore_change[i]
  formula <- as.formula(paste(change_var, "~", baseline_var, "+", paste(comon, collapse  = "+") ))
  changebasicmodels <- lm ( formula, data = ps_df)
  results_key<-paste(baseline_var, change_var, sep = "_")
  changebasicresults[[results_key]] <- summary(changebasicmodels)
  
}   

for (j in names(changebasicresults)){
  print(paste("Regression results for", j))
  print(changebasicresults[[j]])
  
}

modelsummary(changebasicresults, 
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "p.value")


### * refined olive oil  (cont.) -------------------------------------------------------------
comon <- c("ajust_tspre_oliva_v00" ,   "edad_s1",  "sexo_s1")

changebasicresults <- list()

for (i in 1:length(cognitionscore_baseline)){
  baseline_var<-cognitionscore_baseline[i]
  change_var<-cognitionscore_change[i]
  formula <- as.formula(paste(change_var, "~", baseline_var, "+", paste(comon, collapse  = "+") ))
  changebasicmodels <- lm ( formula, data = ps_df)
  results_key<-paste(baseline_var, change_var, sep = "_")
  changebasicresults[[results_key]] <- summary(changebasicmodels)
  
}   

for (j in names(changebasicresults)){
  print(paste("Regression results for", j))
  print(changebasicresults[[j]])
  
}
modelsummary(changebasicresults, 
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "p.value")


# Socio-demographically adjusted model ------------------------------------------------------------- 
# baseline age, sex, baseline cognitive test score, PREDIMED-Plus randomized groups, geographical area, educational level, civil status

### * total olive oil  (cat.) -------------------------------------------------------------
comon <- c("ter_ajust_olivatot_v00", "edad_s1",  "sexo_s1", "grupo_int_v00" , "area", "edu_level_v00" , "civil_status_v00")

changebasicresults <- list()

for (i in 1:length(cognitionscore_baseline)){
  baseline_var<-cognitionscore_baseline[i]
  change_var<-cognitionscore_change[i]
  formula <- as.formula(paste(change_var, "~", baseline_var, "+", paste(comon, collapse  = "+") ))
  changebasicmodels <- lm ( formula, data = ps_df)
  results_key<-paste(change_var, baseline_var, sep = "_")
  changebasicresults[[results_key]] <- summary(changebasicmodels)
  
}   

for (j in names(changebasicresults)){
  print(paste("Regression results for", j))
  print(changebasicresults[[j]])
  
}

modelsummary(changebasicresults, 
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "p.value")

### total olive oil  (p-trend)
comon <- c("ptrend_ter_ajust_olivatot_v00" ,   "edad_s1",  "sexo_s1", "grupo_int_v00" ,  "area", "edu_level_v00" , "civil_status_v00")

changebasicresults <- list()
for (i in 1:length(cognitionscore_baseline)){
  baseline_var<-cognitionscore_baseline[i]
  change_var<-cognitionscore_change[i]
  formula <- as.formula(paste(change_var, "~", baseline_var, "+", paste(comon, collapse  = "+") ))
  changebasicmodels <- lm ( formula, data = ps_df)
  results_key<-paste(baseline_var, change_var, sep = "_")
  changebasicresults[[results_key]] <- summary(changebasicmodels)
  
}   

for (j in names(changebasicresults)){
  print(paste("Regression results for", j))
  print(changebasicresults[[j]])
  
}

modelsummary(changebasicresults, 
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "p.value")


### * total olive oil  (cont.) -------------------------------------------------------------
comon <- c("ajust_tspolivatot_v00" , "edad_s1",  "sexo_s1", "grupo_int_v00" ,  "area", "edu_level_v00" , "civil_status_v00")

changebasicresults <- list()
for (i in 1:length(cognitionscore_baseline)){
  baseline_var<-cognitionscore_baseline[i]
  change_var<-cognitionscore_change[i]
  formula <- as.formula(paste(change_var, "~", baseline_var, "+", paste(comon, collapse  = "+") ))
  changebasicmodels <- lm ( formula, data = ps_df)
  results_key<-paste(baseline_var, change_var, sep = "_")
  changebasicresults[[results_key]] <- summary(changebasicmodels)
  
}   

for (j in names(changebasicresults)){
  print(paste("Regression results for", j))
  print(changebasicresults[[j]])
  
}
modelsummary(changebasicresults, 
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "p.value")

### * extra virgin olive oil  (cat.) -------------------------------------------------------------
comon <- c("ter_ajust_ac_olivavir_v00" ,  "edad_s1",  "sexo_s1", "grupo_int_v00" ,  "area", "edu_level_v00" , "civil_status_v00")

changebasicresults <- list()
for (i in 1:length(cognitionscore_baseline)){
  baseline_var<-cognitionscore_baseline[i]
  change_var<-cognitionscore_change[i]
  formula <- as.formula(paste(change_var, "~", baseline_var, "+", paste(comon, collapse  = "+") ))
  changebasicmodels <- lm ( formula, data = ps_df)
  results_key<-paste(change_var, baseline_var, sep = "_")
  changebasicresults[[results_key]] <- summary(changebasicmodels)
  
}   

for (j in names(changebasicresults)){
  print(paste("Regression results for", j))
  print(changebasicresults[[j]])
  
}

modelsummary(changebasicresults, 
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "p.value")

### extra virgin olive oil  (p-trend)
comon <- c("ptrend_ter_ajust_ac_olivavir_v00" ,   "edad_s1",  "sexo_s1", "grupo_int_v00" ,  "area", "edu_level_v00" , "civil_status_v00")

changebasicresults <- list()
for (i in 1:length(cognitionscore_baseline)){
  baseline_var<-cognitionscore_baseline[i]
  change_var<-cognitionscore_change[i]
  formula <- as.formula(paste(change_var, "~", baseline_var, "+", paste(comon, collapse  = "+") ))
  changebasicmodels <- lm ( formula, data = ps_df)
  results_key<-paste(baseline_var, change_var, sep = "_")
  changebasicresults[[results_key]] <- summary(changebasicmodels)
  
}   

for (j in names(changebasicresults)){
  print(paste("Regression results for", j))
  print(changebasicresults[[j]])
  
}

modelsummary(changebasicresults, 
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "p.value")


### * extra virgin olive oil  (cont.) -------------------------------------------------------------
comon <- c("ajust_tspac_olivavir_v00" ,   "edad_s1",  "sexo_s1", "grupo_int_v00" ,  "area", "edu_level_v00" , "civil_status_v00")

changebasicresults <- list()
for (i in 1:length(cognitionscore_baseline)){
  baseline_var<-cognitionscore_baseline[i]
  change_var<-cognitionscore_change[i]
  formula <- as.formula(paste(change_var, "~", baseline_var, "+", paste(comon, collapse  = "+") ))
  changebasicmodels <- lm ( formula, data = ps_df)
  results_key<-paste(baseline_var, change_var, sep = "_")
  changebasicresults[[results_key]] <- summary(changebasicmodels)
  
}   

for (j in names(changebasicresults)){
  print(paste("Regression results for", j))
  print(changebasicresults[[j]])
  
}
modelsummary(changebasicresults, 
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "p.value")

### * refined olive oil  (cat.) -------------------------------------------------------------
comon <- c("ter_ajust_re_oliva_v00" , "edad_s1",  "sexo_s1", "grupo_int_v00" ,  "area", "edu_level_v00" , "civil_status_v00")

changebasicresults <- list()
for (i in 1:length(cognitionscore_baseline)){
  baseline_var<-cognitionscore_baseline[i]
  change_var<-cognitionscore_change[i]
  formula <- as.formula(paste(change_var, "~", baseline_var, "+", paste(comon, collapse  = "+") ))
  changebasicmodels <- lm ( formula, data = ps_df)
  results_key<-paste(change_var, baseline_var, sep = "_")
  changebasicresults[[results_key]] <- summary(changebasicmodels)
  
}   

for (j in names(changebasicresults)){
  print(paste("Regression results for", j))
  print(changebasicresults[[j]])
  
}

modelsummary(changebasicresults, 
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "p.value")

### refined olive oil  (p-trend)
comon <- c("ptrend_ter_ajust_re_oliva_v00" ,  "edad_s1",  "sexo_s1", "grupo_int_v00" ,  "area", "edu_level_v00" , "civil_status_v00")


changebasicresults <- list()
for (i in 1:length(cognitionscore_baseline)){
  baseline_var<-cognitionscore_baseline[i]
  change_var<-cognitionscore_change[i]
  formula <- as.formula(paste(change_var, "~", baseline_var, "+", paste(comon, collapse  = "+") ))
  changebasicmodels <- lm ( formula, data = ps_df)
  results_key<-paste(baseline_var, change_var, sep = "_")
  changebasicresults[[results_key]] <- summary(changebasicmodels)
  
}   

for (j in names(changebasicresults)){
  print(paste("Regression results for", j))
  print(changebasicresults[[j]])
  
}

modelsummary(changebasicresults, 
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "p.value")


### * refined olive oil  (cont.) -------------------------------------------------------------
comon <- c("ajust_tspre_oliva_v00" ,  "edad_s1",  "sexo_s1", "grupo_int_v00" ,  "area", "edu_level_v00" , "civil_status_v00")

changebasicresults <- list()
for (i in 1:length(cognitionscore_baseline)){
  baseline_var<-cognitionscore_baseline[i]
  change_var<-cognitionscore_change[i]
  formula <- as.formula(paste(change_var, "~", baseline_var, "+", paste(comon, collapse  = "+") ))
  changebasicmodels <- lm ( formula, data = ps_df)
  results_key<-paste(baseline_var, change_var, sep = "_")
  changebasicresults[[results_key]] <- summary(changebasicmodels)
  
}   

for (j in names(changebasicresults)){
  print(paste("Regression results for", j))
  print(changebasicresults[[j]])
  
}
modelsummary(changebasicresults, 
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "p.value")


# Fully adjusted model ------------------------------------------------------------- 
#baseline age, sex, baseline cognitive test score, PREDIMED-Plus randomized groups, geographical area, educational level, civil status
#BMI, physical activity, smoking status, alcohol consumption, depressive symptomatology, diabetes prevalence, hypertension prevalence, hypercholesterolemia prevalence, 
#adherence to Mediterranean diet (12 scores)

### * total olive oil  (cat.) -------------------------------------------------------------
comon <- c("ter_ajust_olivatot_v00", "edad_s1",  "sexo_s1", "grupo_int_v00" ,"area", "edu_level_v00" , "civil_status_v00",
           "imc_v00" , "exercise_day_v00",  "smoke_status_v00" , "alcoholg_v00" , "qalcoholg_v00" , 
           "depression_v00" , "diab_prev_s1" ,"hta_v00", "colest_v00" ,  "MED12score_v00")

changebasicresults <- list()

for (i in 1:length(cognitionscore_baseline)){
  baseline_var<-cognitionscore_baseline[i]
  change_var<-cognitionscore_change[i]
  formula <- as.formula(paste(change_var, "~", baseline_var, "+", paste(comon, collapse  = "+") ))
  changebasicmodels <- lm ( formula, data = ps_df)
  results_key<-paste(change_var, baseline_var, sep = "_")
  changebasicresults[[results_key]] <- summary(changebasicmodels)
  
}   

for (j in names(changebasicresults)){
  print(paste("Regression results for", j))
  print(changebasicresults[[j]])
  
}

modelsummary(changebasicresults, 
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "p.value")

### total olive oil  (p-trend)
comon <- c("ptrend_ter_ajust_olivatot_v00" ,   "edad_s1",  "sexo_s1", "grupo_int_v00" , "area", "edu_level_v00" , "civil_status_v00",
           "imc_v00" , "exercise_day_v00",  "smoke_status_v00" , "alcoholg_v00" , "qalcoholg_v00" , 
           "depression_v00" , "diab_prev_s1" ,"hta_v00", "colest_v00" ,  "MED12score_v00")

changebasicresults <- list()
for (i in 1:length(cognitionscore_baseline)){
  baseline_var<-cognitionscore_baseline[i]
  change_var<-cognitionscore_change[i]
  formula <- as.formula(paste(change_var, "~", baseline_var, "+", paste(comon, collapse  = "+") ))
  changebasicmodels <- lm ( formula, data = ps_df)
  results_key<-paste(baseline_var, change_var, sep = "_")
  changebasicresults[[results_key]] <- summary(changebasicmodels)
  
}   

for (j in names(changebasicresults)){
  print(paste("Regression results for", j))
  print(changebasicresults[[j]])
  
}

modelsummary(changebasicresults, 
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "p.value")


### * total olive oil  (cont.) -------------------------------------------------------------
comon <- c("ajust_tspolivatot_v00" , "edad_s1",  "sexo_s1", "grupo_int_v00" , "area", "edu_level_v00" , "civil_status_v00",
           "imc_v00" , "exercise_day_v00",  "smoke_status_v00" , "alcoholg_v00" , "qalcoholg_v00" , 
           "depression_v00" , "diab_prev_s1" ,"hta_v00", "colest_v00" ,  "MED12score_v00")

changebasicresults <- list()
for (i in 1:length(cognitionscore_baseline)){
  baseline_var<-cognitionscore_baseline[i]
  change_var<-cognitionscore_change[i]
  formula <- as.formula(paste(change_var, "~", baseline_var, "+", paste(comon, collapse  = "+") ))
  changebasicmodels <- lm ( formula, data = ps_df)
  results_key<-paste(baseline_var, change_var, sep = "_")
  changebasicresults[[results_key]] <- summary(changebasicmodels)
  
}   

for (j in names(changebasicresults)){
  print(paste("Regression results for", j))
  print(changebasicresults[[j]])
  
}
modelsummary(changebasicresults, 
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "p.value")


### * extra virgin olive oil  (cat.) -------------------------------------------------------------
comon <- c("ter_ajust_ac_olivavir_v00" ,  "edad_s1",  "sexo_s1", "grupo_int_v00" , "area", "edu_level_v00" , "civil_status_v00",
           "imc_v00" , "exercise_day_v00",  "smoke_status_v00" , "alcoholg_v00" , "qalcoholg_v00" , 
           "depression_v00" , "diab_prev_s1" ,"hta_v00", "colest_v00" ,  "MED12score_v00")

changebasicresults <- list()
for (i in 1:length(cognitionscore_baseline)){
  baseline_var<-cognitionscore_baseline[i]
  change_var<-cognitionscore_change[i]
  formula <- as.formula(paste(change_var, "~", baseline_var, "+", paste(comon, collapse  = "+") ))
  changebasicmodels <- lm ( formula, data = ps_df)
  results_key<-paste(change_var, baseline_var, sep = "_")
  changebasicresults[[results_key]] <- summary(changebasicmodels)
  
}   

for (j in names(changebasicresults)){
  print(paste("Regression results for", j))
  print(changebasicresults[[j]])
  
}

modelsummary(changebasicresults, 
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "p.value")

### extra virgin olive oil  (p-trend)
comon <- c("ptrend_ter_ajust_ac_olivavir_v00" ,   "edad_s1",  "sexo_s1", "grupo_int_v00" , "area", "edu_level_v00" , "civil_status_v00",
           "imc_v00" , "exercise_day_v00",  "smoke_status_v00" , "alcoholg_v00" , "qalcoholg_v00" , 
           "depression_v00" , "diab_prev_s1" ,"hta_v00", "colest_v00" ,  "MED12score_v00")

changebasicresults <- list()
for (i in 1:length(cognitionscore_baseline)){
  baseline_var<-cognitionscore_baseline[i]
  change_var<-cognitionscore_change[i]
  formula <- as.formula(paste(change_var, "~", baseline_var, "+", paste(comon, collapse  = "+") ))
  changebasicmodels <- lm ( formula, data = ps_df)
  results_key<-paste(baseline_var, change_var, sep = "_")
  changebasicresults[[results_key]] <- summary(changebasicmodels)
  
}   

for (j in names(changebasicresults)){
  print(paste("Regression results for", j))
  print(changebasicresults[[j]])
  
}

modelsummary(changebasicresults, 
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "p.value")

### * extra virgin olive oil  (cont.) -------------------------------------------------------------
comon <- c("ajust_tspac_olivavir_v00" ,   "edad_s1",  "sexo_s1", "grupo_int_v00" ,"area", "edu_level_v00" , "civil_status_v00",
           "imc_v00" , "exercise_day_v00",  "smoke_status_v00" , "alcoholg_v00" , "qalcoholg_v00" , 
           "depression_v00" , "diab_prev_s1" ,"hta_v00", "colest_v00" ,  "MED12score_v00")

changebasicresults <- list()
for (i in 1:length(cognitionscore_baseline)){
  baseline_var<-cognitionscore_baseline[i]
  change_var<-cognitionscore_change[i]
  formula <- as.formula(paste(change_var, "~", baseline_var, "+", paste(comon, collapse  = "+") ))
  changebasicmodels <- lm ( formula, data = ps_df)
  results_key<-paste(baseline_var, change_var, sep = "_")
  changebasicresults[[results_key]] <- summary(changebasicmodels)
  
}   

for (j in names(changebasicresults)){
  print(paste("Regression results for", j))
  print(changebasicresults[[j]])
  
}
modelsummary(changebasicresults, 
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "p.value")


### * refined olive oil  (cat.) -------------------------------------------------------------
comon <- c("ter_ajust_re_oliva_v00" , "edad_s1",  "sexo_s1", "grupo_int_v00" ,"area", "edu_level_v00" , "civil_status_v00",
           "imc_v00" , "exercise_day_v00",  "smoke_status_v00" , "alcoholg_v00" , "qalcoholg_v00" , 
           "depression_v00" , "diab_prev_s1" ,"hta_v00", "colest_v00" ,  "MED12score_v00")

changebasicresults <- list()
for (i in 1:length(cognitionscore_baseline)){
  baseline_var<-cognitionscore_baseline[i]
  change_var<-cognitionscore_change[i]
  formula <- as.formula(paste(change_var, "~", baseline_var, "+", paste(comon, collapse  = "+") ))
  changebasicmodels <- lm ( formula, data = ps_df)
  results_key<-paste(change_var, baseline_var, sep = "_")
  changebasicresults[[results_key]] <- summary(changebasicmodels)
  
}   

for (j in names(changebasicresults)){
  print(paste("Regression results for", j))
  print(changebasicresults[[j]])
  
}

modelsummary(changebasicresults, 
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "p.value")

### refined olive oil  (p-trend)
comon <- c("ptrend_ter_ajust_re_oliva_v00" ,  "edad_s1",  "sexo_s1", "grupo_int_v00" ,"area", "edu_level_v00" , "civil_status_v00",
           "imc_v00" , "exercise_day_v00",  "smoke_status_v00" , "alcoholg_v00" , "qalcoholg_v00" , 
           "depression_v00" , "diab_prev_s1" ,"hta_v00", "colest_v00" ,  "MED12score_v00")

changebasicresults <- list()
for (i in 1:length(cognitionscore_baseline)){
  baseline_var<-cognitionscore_baseline[i]
  change_var<-cognitionscore_change[i]
  formula <- as.formula(paste(change_var, "~", baseline_var, "+", paste(comon, collapse  = "+") ))
  changebasicmodels <- lm ( formula, data = ps_df)
  results_key<-paste(baseline_var, change_var, sep = "_")
  changebasicresults[[results_key]] <- summary(changebasicmodels)
  
}   

for (j in names(changebasicresults)){
  print(paste("Regression results for", j))
  print(changebasicresults[[j]])
  
}

modelsummary(changebasicresults, 
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "p.value")


### * refined olive oil  (cont.) -------------------------------------------------------------
comon <- c("ajust_tspre_oliva_v00" ,  "edad_s1",  "sexo_s1", "grupo_int_v00" ,"area", "edu_level_v00" , "civil_status_v00",
           "imc_v00" , "exercise_day_v00",  "smoke_status_v00" , "alcoholg_v00" , "qalcoholg_v00" , 
           "depression_v00" , "diab_prev_s1" ,"hta_v00", "colest_v00" ,  "MED12score_v00")

changebasicresults <- list()
for (i in 1:length(cognitionscore_baseline)){
  baseline_var<-cognitionscore_baseline[i]
  change_var<-cognitionscore_change[i]
  formula <- as.formula(paste(change_var, "~", baseline_var, "+", paste(comon, collapse  = "+") ))
  changebasicmodels <- lm ( formula, data = ps_df)
  results_key<-paste(baseline_var, change_var, sep = "_")
  changebasicresults[[results_key]] <- summary(changebasicmodels)
  
}   

for (j in names(changebasicresults)){
  print(paste("Regression results for", j))
  print(changebasicresults[[j]])
  
}
modelsummary(changebasicresults, 
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "p.value")

