# Title: Producing descriptive statistics
# Author: Jiaqi Ni [ORCID](https://orcid.org/0000-0001-5625-0373)
# Description: This is the analysis script for descriptive statistics

# Set up working environment  ------------------------------------------------------------------
rm(list=ls())
getwd()
setwd("/home/jiaqi/Jiaqi_Ni/OliveMicrobiomeCognition_Microbiome")
options(max.print=1000000)
options(scipen = 999)

# Load dependencies -------------------------------------------------------
library(dplyr)
library(tableone)
library(openxlsx)

# Data preparation -----------------------------------------------------------
ps_df<- readRDS("DATA/metadata_FINAL_110324.rds") #load metadata
# Define and fix continuous and categorical variables
## Vector of variables to summarize
myVars <- c("sexo_s1", "edad_s1",  "nodo",  "area","grupo_int_v00", "edu_level_v00", "civil_status_v00", 
            "olivatot_v00",  "ac_olivavir_v00", "re_oliva_v00",  "otheroil_v00", "ajust_otheroil_v00" , "ajust_tspotheroil_v00",
            "ajust_olivatot_v00", "ajust_ac_olivavir_v00", "ajust_re_oliva_v00", 
            "ter_ajust_olivatot_v00", "ter_ajust_ac_olivavir_v00", "ter_ajust_re_oliva_v00", 
            "ptrend_ter_ajust_olivatot_v00", "ptrend_ter_ajust_ac_olivavir_v00", "ptrend_ter_ajust_re_oliva_v00", 
            "MMSE_v00", "MMSE_v02", "CDT_v00", "CDT_v02", "VFTa_v00", "VFTa_v02", "VFTp_v00", "VFTp_v02", 
            "TMTa_v00", "TMTa_v02", "TMTb_v00", "TMTb_v02", "DSTf_v00", "DSTf_v02", "DSTb_v00", "DSTb_v02", 
            "zMMSE_v00", "zMMSE_v02", "zCDT_v00", "zCDT_v02", "zVFTa_v00", "zVFTa_v02", "zVFTp_v00", "zVFTp_v02", 
            "zTMTa_v00", "zTMTa_v02", "zTMTb_v00", "zTMTb_v02", "zDSTf_v00", "zDSTf_v02", "zDSTb_v00", "zDSTb_v02", 
            "GCF_v00", "GCF_v02", "genCF_v00", "genCF_v02", "ExF_v00", "ExF_v02", "attention_v00", "attention_v02", 
            "language_v00", "language_v02", "zGCF_v00", "zGCF_v02", "zgenCF_v00", "zgenCF_v02", "zExF_v00", "zExF_v02", 
            "zattention_v00", "zattention_v02", "zlanguage_v00", "zlanguage_v02", 
            "zMMSE_2yc", "zCDT_2yc", "zVFTa_2yc", "zVFTp_2yc", "zTMTa_2yc", "zTMTb_2yc", "zDSTf_2yc", "zDSTb_2yc", 
            "zGCF_2yc", "zgenCF_2yc", "zExF_2yc", "zattention_2yc", "zlanguage_2yc", 
            "imc_v00",  "waist_cir_v00", "diab_prev_s1", "hta_v00", "colest_v00", "depression_v00",  
            "antihypertensive_v00", "antidiabetic_v00", "hypolipidemic_v00", "antidepressant_v00", 
            "exercise_day_v00", "ter_exercise_day_v00", "smoke_status_v00", "MED14score_v00", "MED12score_v00", "ter_MED12score_v00", 
            "energiat_v00", "alcoholg_v00", "ter_alcoholg_v00", "qalcoholg_v00", "hc_v00", "prot_v00", "gratot_v00", "mo_v00", "po_v00", "sa_v00",
            "porc_hc_v00", "porc_pr_v00", "porc_gr_v00", "porc_mo_v00", "porc_po_v00", "porc_sa_v00", 
            "fibra_v00", "col_v00", "fit_v00", "trans_v00", "linoleico_v00", "linolenico_v00", "omega3_v00", "n3marinos_v00", 
            "verdutot_v00", "frutatot_v00", "legumbre_v00", "cereal_v00", "lacteos_v00", "carnicos_v00", "pescados_v00", "fsecos_v00", "gallet_v00",
            "ajust_hc_v00", "ajust_prot_v00", "ajust_gratot_v00", "ajust_mo_v00", "ajust_po_v00", "ajust_sa_v00", "ajust_porc_hc_v00", "ajust_porc_pr_v00", "ajust_porc_gr_v00", "ajust_porc_mo_v00", "ajust_porc_po_v00", "ajust_porc_sa_v00", 
            "ajust_fibra_v00", "ajust_verdutot_v00", "ajust_frutatot_v00", "ajust_legumbre_v00", 
            "ajust_cereal_v00", "ajust_lacteos_v00", "ajust_carnicos_v00", "ajust_pescados_v00", "ajust_fsecos_v00", "ajust_gallet_v00")

catVars <- c("sexo_s1",  "nodo", "area", "grupo_int_v00", "edu_level_v00", "civil_status_v00", 
             "ter_ajust_olivatot_v00", "ter_ajust_ac_olivavir_v00", "ter_ajust_re_oliva_v00", 
              "diab_prev_s1", "hta_v00", "colest_v00", "depression_v00",  "antihypertensive_v00", "antidepressant_v00", 
             "antidiabetic_v00", "hypolipidemic_v00","ter_exercise_day_v00",
             "smoke_status_v00","ter_MED12score_v00","ter_alcoholg_v00")
ps_df[,catVars] <- lapply(ps_df[,catVars], as.factor)

conVars <- c("edad_s1", "olivatot_v00",  "ac_olivavir_v00", "re_oliva_v00",  
             "ajust_olivatot_v00", "ajust_ac_olivavir_v00", "ajust_re_oliva_v00", "otheroil_v00", "ajust_otheroil_v00" , "ajust_tspotheroil_v00",
             "ptrend_ter_ajust_olivatot_v00", "ptrend_ter_ajust_ac_olivavir_v00", "ptrend_ter_ajust_re_oliva_v00",
             "MMSE_v00", "MMSE_v02", "CDT_v00", "CDT_v02", "VFTa_v00", "VFTa_v02", "VFTp_v00", "VFTp_v02", 
             "TMTa_v00", "TMTa_v02", "TMTb_v00", "TMTb_v02", "DSTf_v00", "DSTf_v02", "DSTb_v00", "DSTb_v02", 
             "zMMSE_v00", "zMMSE_v02", "zCDT_v00", "zCDT_v02", "zVFTa_v00", "zVFTa_v02", "zVFTp_v00", "zVFTp_v02", 
             "zTMTa_v00", "zTMTa_v02", "zTMTb_v00", "zTMTb_v02", "zDSTf_v00", "zDSTf_v02", "zDSTb_v00", "zDSTb_v02", 
             "GCF_v00", "GCF_v02", "genCF_v00", "genCF_v02", "ExF_v00", "ExF_v02", "attention_v00", "attention_v02", 
             "language_v00", "language_v02", "zGCF_v00", "zGCF_v02", "zgenCF_v00", "zgenCF_v02", "zExF_v00", "zExF_v02", 
             "zattention_v00", "zattention_v02", "zlanguage_v00", "zlanguage_v02", 
             "zMMSE_2yc", "zCDT_2yc", "zVFTa_2yc", "zVFTp_2yc", "zTMTa_2yc", "zTMTb_2yc", "zDSTf_2yc", "zDSTb_2yc", 
             "zGCF_2yc", "zgenCF_2yc", "zExF_2yc", "zattention_2yc", "zlanguage_2yc", 
             "imc_v00", "waist_cir_v00", "exercise_day_v00", "alcoholg_v00", "qalcoholg_v00", "MED12score_v00", 
             "hc_v00", "prot_v00", "gratot_v00", "mo_v00", "po_v00", "sa_v00",
             "porc_hc_v00", "porc_pr_v00", "porc_gr_v00", "porc_mo_v00", "porc_po_v00", "porc_sa_v00", 
             "fibra_v00", "col_v00", "fit_v00", "trans_v00", "linoleico_v00", "linolenico_v00", "omega3_v00", "n3marinos_v00", 
             "verdutot_v00", "frutatot_v00", "legumbre_v00", "cereal_v00", "lacteos_v00", "carnicos_v00", "pescados_v00", "fsecos_v00", "gallet_v00",
             "ajust_hc_v00", "ajust_prot_v00", "ajust_gratot_v00", "ajust_mo_v00", "ajust_po_v00", "ajust_sa_v00", "ajust_porc_hc_v00", "ajust_porc_pr_v00", "ajust_porc_gr_v00", "ajust_porc_mo_v00", "ajust_porc_po_v00", "ajust_porc_sa_v00", 
             "ajust_fibra_v00", "ajust_verdutot_v00", "ajust_frutatot_v00", "ajust_legumbre_v00", 
             "ajust_cereal_v00", "ajust_lacteos_v00", "ajust_carnicos_v00", "ajust_pescados_v00", "ajust_fsecos_v00", "ajust_gallet_v00")
ps_df[,conVars] <- lapply(ps_df[,conVars], as.numeric)
# Descriptive statistics -----------------------------------------------------------
## Create Table 1 grouped by ter_ajust_olivatot_v00
table1 <- CreateTableOne(vars = myVars, strata = c("ter_ajust_olivatot_v00"), data = ps_df, factorVars = catVars, addOverall = TRUE)
table1$CatTable #for categorical variables
table1$ContTable #for continuous variables
# Save results -----------------------------------------------------------
table1df <- as.data.frame(print(table1, showAllLevels = TRUE, formatOptions = list(big.mark = ","), smd = TRUE))
write.xlsx(table1df, file = "Output/descriptives.xlsx", rowNames = TRUE)

