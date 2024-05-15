# Title: Data preparation 
# Author: Jiaqi Ni [ORCID](https://orcid.org/0000-0001-5625-0373)
# Description: This is the script for subsetting and preparing data for downstream analyses


# Set up working environment  ------------------------------------------------------------------
rm(list=ls())
getwd()
setwd("/home/jiaqi/Jiaqi_Ni/OliveMicrobiomeCognition_Microbiome")
options(max.print=1000000)
options(scipen = 999)

# Load dependencies -------------------------------------------------------
library(haven) #for reading STATA data
library(dplyr)
library(summarytools)
library(labelled)
library(purrr)
library(phyloseq)

# Load raw data -----------------------------------------------------------

PMplus.all <- read_dta("./DATA/202012220958_PREDIMEDplus_2anys_2020-12-22.dta") #load PREDIMED-Plus original metadata 
ps <- readRDS("./DATA/predimed_phyloseq_all.rds") #load original phyloseq object 

# Prepare the final dataset -----------------------------------------------
# * Flowchart -------------------------------------------------
# The flowchart starts
# * * Keep participants with baseline gut microbiota data  -------------------------------------------------
ps <- subset_samples(ps, Timepoint_hr == "Baseline")
# From the Sample Data of this phyloseq object we get "paciente" from the column 24 (OriginalID) 
paciente <- gsub("P20|\\.A\\.Amplicon|\\.B\\.Amplicon|\\.A\\.r\\.Amplicon|\\.A\\.d\\.Amplicon|\\.B\\.r\\.Amplicon|P2", "", sample_data(ps)$OriginalID)
sample_data(ps)$paciente <- paciente
#Keep the columns of this metadata that are necessary
sample_data(ps)$SampleIDcounts <- rownames(sample_data(ps)) #create a new variable with rowname 
keep.cols <- c("paciente", "SampleIdentifier", "SampleID", "OriginalID", "SampleIDs","SampleIDcounts", "Timepoint_hr")
sample_data(ps) <- sample_data(ps) [, keep.cols]

# * * Clean phyloseq samples -------------------------------------------------
# Remove duplicated or replicated sample with lower DNA concentration
to.move <- c("105.B.d.Amplicon", "133.A.Amplicon", "134.A.d.Amplicon", "173.A.d.Amplicon", "19.A.Amplicon", "20.B.Amplicon", "225.B.d.Amplicon", "233.A.Amplicon", "293.A.d.Amplicon", "383.A.Amplicon", "395.A.r.Amplicon", "477.A.r.Amplicon", "512.B.Amplicon", "519.A.r.Amplicon", "568.B.r.Amplicon", "600.A.d.Amplicon", "667.B.Amplicon", "682.B.Amplicon", "74.B.Amplicon", "743.B.r.Amplicon", "762.A.d.Amplicon", "812.B.r.Amplicon", "819.B.Amplicon", "820.B.Amplicon", "886.A.r.Amplicon") 
ps <- prune_samples(!(sample_names(ps) %in% to.move), ps)
# Remove samples from antibiotic users
to.move2 <- c("80.A.Amplicon", "80.B.Amplicon", "89.A.Amplicon", "89.B.Amplicon", "100.A.Amplicon", "100.B.Amplicon", "172.A.Amplicon", "172.B.Amplicon","377.A.Amplicon", "377.B.Amplicon", "405.A.Amplicon", "405.B.Amplicon","472.A.Amplicon", "472.B.Amplicon","600.A.Amplicon", "600.B.Amplicon","859.A.Amplicon", "859.B.Amplicon") 
ps <- prune_samples(!(sample_names(ps) %in% to.move2), ps)
# retrieve feature counts table
counts = as.data.frame(otu_table(ps))
# retrieve the original metadata
metadata = as.data.frame(as.matrix(sample_data(ps)))
# Merge the 2 data frames based on the "paciente" columnn, the original PM-plus metadata and the prepared metadata
merged <- merge(metadata, PMplus.all, by = "paciente", all.x= TRUE) #all.x= TRUE ensures that all rows from the original sample_data are retained even if there are missing matches in the "add.samdf"

PMplus.all <- merged 
PMplus.ok<-PMplus.all[!is.na(PMplus.all$grupo_int_v00),] #Filter out Individuals that were not randomized (actual beginning of the flowchart)
PMplus.ok<-PMplus.ok[!is.na(PMplus.ok$energiat_v00),] #Filter out Individuals w/o FFQ at baseline
PMplus.ok$energia_v00ok <- with(PMplus.ok, (sexo_s1==0 & (energiat_v00>=800 & energiat_v00<=4000)| sexo_s1==1 & (energiat_v00>=500 & energiat_v00<=3500))) #Create Extreme Kcal Intakes
PMplus.ok<-PMplus.ok[PMplus.ok$energia_v00ok==TRUE,] #Filter out Extreme Kcal Intakes 
PMplus.ok<-PMplus.ok[PMplus.ok$olivatot_v00<=100,] #Filter out Individuals reporting olive oil intake greater than 100g/d 
PMplus.ok <- PMplus.ok %>%
  rename(
    MMSE_v00 = puntuacion_mmse_v00, 
    MMSE_v02 = puntuacion_mmse_v02, 
    CDT_v00 = puntua_reloj_v00, 
    CDT_v02 = puntua_reloj_v02,
    VFTa_v00 = total_animales_v00, 
    VFTa_v02 = total_animales_v02, 
    VFTp_v00 = total_letra_p_v00, 
    VFTp_v02 = total_letra_p_v02, 
    TMTa_v00 = tiempo_totala_v00, 
    TMTa_v02 = tiempo_totala_v02, 
    TMTb_v00 = tiempo_totalb_v00, 
    TMTb_v02 = tiempo_totalb_v02, 
    DSTf_v00 = orden_directo_v00, 
    DSTf_v02 = orden_directo_v02, 
    DSTb_v00 = orden_inverso_v00, 
    DSTb_v02 = orden_inverso_v02 
  ) # rename variables
PMplus.ok <- subset(PMplus.ok, !(is.na(MMSE_v00) | is.na(MMSE_v02) | is.na(CDT_v00) | is.na(CDT_v02) | 
                                   is.na(VFTa_v00) | is.na(VFTa_v02) | is.na(VFTp_v00) | is.na(VFTp_v02) |
                                   is.na(TMTa_v00) | is.na(TMTa_v02) | is.na(TMTb_v00) | is.na(TMTb_v02) | 
                                   is.na(DSTf_v00) | is.na(DSTf_v02) | is.na(DSTb_v00) | is.na(DSTb_v02))) #Filter out Individuals without cognitive test scores at baseline & year 2 
# * Prepare metadata -------------------------------------------------
# * * Exposure -------------------------------------------------
summary(PMplus.ok[c("olivatot_v00", "ac_oliva_v00", "ac_olivavir_v00", "ac_orujo_v00")])
PMplus.ok$re_oliva_v00 <- PMplus.ok$ac_oliva_v00 + PMplus.ok$ac_orujo_v00 #generate common olive oil variable
PMplus.ok$otheroil_v00 <- PMplus.ok$ac_maiz_v00 + PMplus.ok$ac_girasol_v00 + PMplus.ok$ac_soja_v00 + PMplus.ok$ac_mezcla_v00 # generate variable for other types of oil

summary(PMplus.ok[c("olivatot_v00", "ac_olivavir_v00", "re_oliva_v00", "otheroil_v00")])
PMplus.ok$tspolivatot_v00 <- PMplus.ok$olivatot_v00/10 #per 10g 
PMplus.ok$tspac_olivavir_v00 <- PMplus.ok$ac_olivavir_v00/10 #per 10g 
PMplus.ok$tspre_oliva_v00 <- PMplus.ok$re_oliva_v00/10 #per 10g 
PMplus.ok$tspotheroil_v00 <- PMplus.ok$otheroil_v00/10 #per 10g 
summary(PMplus.ok[c("tspolivatot_v00", "tspac_olivavir_v00", "tspre_oliva_v00", "tspotheroil_v00")])


# Energy adjusted using the residual methods
var.list <- c("olivatot_v00", "ac_olivavir_v00", "re_oliva_v00", "otheroil_v00", "tspolivatot_v00", "tspac_olivavir_v00", "tspre_oliva_v00", "tspotheroil_v00") # List of variables that need to be adjusted by energy intake 
for (var in var.list) { # perform the loop over the variables 
  # run the regression 
  model<-lm(paste(var, "energiat_v00", sep = "~"), data = PMplus.ok )
  # get the residuals from the regression 
  residuals <- residuals(model)
  # calculate the mean energy intake 
  energy_mean <- mean(PMplus.ok$energiat_v00, na.rm = TRUE)
  # expected food intake for a person with the mean energy intake
  expected_var <- coef(model)[1] + coef(model)[2] * energy_mean 
  # create the new variable "ajust_var" and add the residuals and expected food intake
  PMplus.ok[paste0("ajust_", var)] <-expected_var + residuals 
}
summary(PMplus.ok[c("ajust_olivatot_v00", "ajust_ac_olivavir_v00", "ajust_re_oliva_v00","ajust_otheroil_v00" , "ajust_tspolivatot_v00", "ajust_tspac_olivavir_v00", "ajust_tspre_oliva_v00","ajust_tspotheroil_v00")])
table(PMplus.ok$ajust_olivatot_v00)
table(PMplus.ok$ajust_ac_olivavir_v00)
table(PMplus.ok$ajust_re_oliva_v00)
# Tertiles of energy-adjusted olive oil consumption
var.adjusted.list <- c("ajust_olivatot_v00", "ajust_ac_olivavir_v00", "ajust_re_oliva_v00")
for (var in var.adjusted.list) {
  PMplus.ok <- PMplus.ok %>%
    mutate(!!paste0("ter_", var) := ntile(!!sym(var), 3))
}

summary(PMplus.ok[c("ter_ajust_olivatot_v00", "ter_ajust_ac_olivavir_v00", "ter_ajust_re_oliva_v00")])
PMplus.ok$ter_ajust_olivatot_v00 <- as.factor(PMplus.ok$ter_ajust_olivatot_v00) #convert the variable as factor
PMplus.ok$ter_ajust_ac_olivavir_v00 <- as.factor(PMplus.ok$ter_ajust_ac_olivavir_v00)
PMplus.ok$ter_ajust_re_oliva_v00 <- as.factor(PMplus.ok$ter_ajust_re_oliva_v00)


summary_stats_oil <- PMplus.ok %>%
  group_by(ter_ajust_olivatot_v00) %>%
  summarise(count = n(),
            mean = mean(olivatot_v00),
            sd = sd(olivatot_v00),
            min = min(olivatot_v00),
            q25 = quantile(olivatot_v00, 0.25),
            median = median(olivatot_v00),
            q75 = quantile(olivatot_v00, 0.75),
            max = max(olivatot_v00))
print(summary_stats_oil) #description of total olive oil consumption per group

summary_stats_extravirgin <- PMplus.ok %>%
  group_by(ter_ajust_ac_olivavir_v00) %>%
  summarise(count = n(),
            mean = mean(ac_olivavir_v00),
            sd = sd(ac_olivavir_v00),
            min = min(ac_olivavir_v00),
            q25 = quantile(ac_olivavir_v00, 0.25),
            median = median(ac_olivavir_v00),
            q75 = quantile(ac_olivavir_v00, 0.75),
            max = max(ac_olivavir_v00))
print(summary_stats_extravirgin) #description of extra virgin olive oil consumption per group

summary_stats_refined <- PMplus.ok %>%
  group_by(ter_ajust_re_oliva_v00) %>%
  summarise(count = n(),
            mean = mean(re_oliva_v00),
            sd = sd(re_oliva_v00),
            min = min(re_oliva_v00),
            q25 = quantile(re_oliva_v00, 0.25),
            median = median(re_oliva_v00),
            q75 = quantile(re_oliva_v00, 0.75),
            max = max(re_oliva_v00))

print(summary_stats_refined) #description of extra virgin olive oil consumption per group

# Create p-trend of categories of energy-adjusted olive oil consumption
## Total olive oil
PMplus.ok <- PMplus.ok %>% mutate(ptrend_ter_ajust_olivatot_v00 = NA_real_)
#Iterate through tertiles and calculate medians
for (q in 1:3){
  PMplus.ok_subset <- PMplus.ok %>%
    filter(ter_ajust_olivatot_v00 == q) %>%
    summarise(median_value = median(ajust_olivatot_v00))
  PMplus.ok$ptrend_ter_ajust_olivatot_v00[PMplus.ok$ter_ajust_olivatot_v00 == q] <- PMplus.ok_subset$median_value
}

## Extra virgin olive oil
PMplus.ok <- PMplus.ok %>% mutate(ptrend_ter_ajust_ac_olivavir_v00 = NA_real_)
#Iterate through tertiles and calculate medians
for (q in 1:3){
  PMplus.ok_subset <- PMplus.ok %>%
    filter(ter_ajust_ac_olivavir_v00 == q) %>%
    summarise(median_value = median(ajust_ac_olivavir_v00))
  PMplus.ok$ptrend_ter_ajust_ac_olivavir_v00[PMplus.ok$ter_ajust_ac_olivavir_v00 == q] <- PMplus.ok_subset$median_value
}

## Refined olive oil
PMplus.ok <- PMplus.ok %>% mutate(ptrend_ter_ajust_re_oliva_v00 = NA_real_)
#Iterate through tertiles and calculate medians
for (q in 1:3){
  PMplus.ok_subset <- PMplus.ok %>%
    filter(ter_ajust_re_oliva_v00 == q) %>%
    summarise(median_value = median(ajust_re_oliva_v00))
  PMplus.ok$ptrend_ter_ajust_re_oliva_v00[PMplus.ok$ter_ajust_re_oliva_v00 == q] <- PMplus.ok_subset$median_value
}
# * * Outcome -------------------------------------------------
# Generate Z-SCORES (to mean and SD of Baseline) for each test at each visit
#MMSE
##baseline
MMSEmean_v00 <- mean(PMplus.ok$MMSE_v00, na.rm = TRUE)
MMSEsd_v00 <- sd(PMplus.ok$MMSE_v00, na.rm = TRUE)
PMplus.ok <- PMplus.ok %>%
  mutate(zMMSE_v00 = (MMSE_v00 - MMSEmean_v00) / MMSEsd_v00)
##2-year
PMplus.ok <- PMplus.ok %>%
  mutate(zMMSE_v02 = (MMSE_v02 - MMSEmean_v00) / MMSEsd_v00)

#CDT
##baseline
CDTmean_v00 <- mean(PMplus.ok$CDT_v00, na.rm = TRUE)
CDTsd_v00 <- sd(PMplus.ok$CDT_v00, na.rm = TRUE)
PMplus.ok <- PMplus.ok %>%
  mutate(zCDT_v00 = (CDT_v00 - CDTmean_v00) / CDTsd_v00)
##2-year
PMplus.ok <- PMplus.ok %>%
  mutate(zCDT_v02 = (CDT_v02 - CDTmean_v00) / CDTsd_v00)

#VFT-a
##baseline
VFTamean_v00 <- mean(PMplus.ok$VFTa_v00, na.rm = TRUE)
VFTasd_v00 <- sd(PMplus.ok$VFTa_v00, na.rm = TRUE)
PMplus.ok <- PMplus.ok %>%
  mutate(zVFTa_v00 = (VFTa_v00 - VFTamean_v00) / VFTasd_v00)
##2-year
PMplus.ok <- PMplus.ok %>%
  mutate(zVFTa_v02 = (VFTa_v02 - VFTamean_v00) / VFTasd_v00)


#VFT-p
##baseline
VFTpmean_v00 <- mean(PMplus.ok$VFTp_v00, na.rm = TRUE)
VFTpsd_v00 <- sd(PMplus.ok$VFTp_v00, na.rm = TRUE)
PMplus.ok <- PMplus.ok %>%
  mutate(zVFTp_v00 = (VFTp_v00 - VFTpmean_v00) / VFTpsd_v00)
##2-year
PMplus.ok <- PMplus.ok %>%
  mutate(zVFTp_v02 = (VFTp_v02 - VFTpmean_v00) / VFTpsd_v00)


#TMT-a
##baseline
TMTamean_v00 <- mean(PMplus.ok$TMTa_v00, na.rm = TRUE)
TMTasd_v00 <- sd(PMplus.ok$TMTa_v00, na.rm = TRUE)
PMplus.ok <- PMplus.ok %>%
  mutate(zTMTa_v00 = (TMTa_v00 - TMTamean_v00) / TMTasd_v00)
##2-year
PMplus.ok <- PMplus.ok %>%
  mutate(zTMTa_v02 = (TMTa_v02 - TMTamean_v00) / TMTasd_v00)


#TMT-b
##baseline
TMTbmean_v00 <- mean(PMplus.ok$TMTb_v00, na.rm = TRUE)
TMTbsd_v00 <- sd(PMplus.ok$TMTb_v00, na.rm = TRUE)
PMplus.ok <- PMplus.ok %>%
  mutate(zTMTb_v00 = (TMTb_v00 - TMTbmean_v00) / TMTbsd_v00)
##2-year
PMplus.ok <- PMplus.ok %>%
  mutate(zTMTb_v02 = (TMTb_v02 - TMTbmean_v00) / TMTbsd_v00)


#DST-f
##baseline
DSTfmean_v00 <- mean(PMplus.ok$DSTf_v00, na.rm = TRUE)
DSTfsd_v00 <- sd(PMplus.ok$DSTf_v00, na.rm = TRUE)
PMplus.ok <- PMplus.ok %>%
  mutate(zDSTf_v00 = (DSTf_v00 - DSTfmean_v00) / DSTfsd_v00)
##2-year
PMplus.ok <- PMplus.ok %>%
  mutate(zDSTf_v02 = (DSTf_v02 - DSTfmean_v00) / DSTfsd_v00)


#DST-b
##baseline
DSTbmean_v00 <- mean(PMplus.ok$DSTb_v00, na.rm = TRUE)
DSTbsd_v00 <- sd(PMplus.ok$DSTb_v00, na.rm = TRUE)
PMplus.ok <- PMplus.ok %>%
  mutate(zDSTb_v00 = (DSTb_v00 - DSTbmean_v00) / DSTbsd_v00)
##2-year
PMplus.ok <- PMplus.ok %>%
  mutate(zDSTb_v02 = (DSTb_v02 - DSTbmean_v00) / DSTbsd_v00)

# Generate Cognitive Domains
#Global Cognitive Function 
##baseline
PMplus.ok <- PMplus.ok %>%
  mutate(GCF_v00 = (zMMSE_v00 + zCDT_v00 + zVFTa_v00 + zVFTp_v00 - zTMTa_v00 - zTMTb_v00 + zDSTf_v00 + zDSTb_v00) / 8) #TMTa and TMTb, lower score better cognitive performance 
##2-year
PMplus.ok <- PMplus.ok %>%
  mutate(GCF_v02 = (zMMSE_v02 + zCDT_v02 + zVFTa_v02 + zVFTp_v02 - zTMTa_v02 - zTMTb_v02 + zDSTf_v02 + zDSTb_v02) / 8)


#General Cognitive Function 
##baseline
PMplus.ok <- PMplus.ok %>%
  mutate(genCF_v00 = (zMMSE_v00 + zCDT_v00) / 2)
##2-year
PMplus.ok <- PMplus.ok %>%
  mutate(genCF_v02 = (zMMSE_v02 + zCDT_v02) / 2)


#Executive Function 
##baseline
PMplus.ok <- PMplus.ok %>%
  mutate(ExF_v00 = (zVFTa_v00 + zVFTp_v00 - zTMTb_v00 + zDSTb_v00) / 4)
##2-year
PMplus.ok <- PMplus.ok %>%
  mutate(ExF_v02 = (zVFTa_v02 + zVFTp_v02 - zTMTb_v02 + zDSTb_v02) / 4)


#Attention
##baseline
PMplus.ok <- PMplus.ok %>%
  mutate(attention_v00 = (- zTMTa_v00 + zDSTf_v00) / 2)
##2-year
PMplus.ok <- PMplus.ok %>%
  mutate(attention_v02 = (- zTMTa_v02 + zDSTf_v02) / 2)


#Language
##baseline
PMplus.ok <- PMplus.ok %>%
  mutate(language_v00 = (zVFTa_v00 + zVFTp_v00) / 2)
##2-year
PMplus.ok <- PMplus.ok %>%
  mutate(language_v02 = (zVFTa_v02 + zVFTp_v02) / 2)

# Generate Z-SCORES (to mean and SD of Baseline) for each Cognitive Domain at each visit
#Global Cognitive Function
##baseline
GCFmean_v00 <- mean(PMplus.ok$GCF_v00, na.rm = TRUE)
GCFsd_v00 <- sd(PMplus.ok$GCF_v00, na.rm = TRUE)
PMplus.ok <- PMplus.ok %>%
  mutate(zGCF_v00 = (GCF_v00 - GCFmean_v00) / GCFsd_v00)

##2-year
PMplus.ok <- PMplus.ok %>%
  mutate(zGCF_v02 = (GCF_v02 - GCFmean_v00) / GCFsd_v00)


#General Cognitive Function
##baseline
genCFmean_v00 <- mean(PMplus.ok$genCF_v00, na.rm = TRUE)
genCFsd_v00 <- sd(PMplus.ok$genCF_v00, na.rm = TRUE)
PMplus.ok <- PMplus.ok %>%
  mutate(zgenCF_v00 = (genCF_v00 - genCFmean_v00) / genCFsd_v00)

##2-year
PMplus.ok <- PMplus.ok %>%
  mutate(zgenCF_v02 = (genCF_v02 - genCFmean_v00) / genCFsd_v00)


#Executive Function
##baseline
ExFmean_v00 <- mean(PMplus.ok$ExF_v00, na.rm = TRUE)
ExFsd_v00 <- sd(PMplus.ok$ExF_v00, na.rm = TRUE)
PMplus.ok <- PMplus.ok %>%
  mutate(zExF_v00 = (ExF_v00 - ExFmean_v00) / ExFsd_v00)

##2-year
PMplus.ok <- PMplus.ok %>%
  mutate(zExF_v02 = (ExF_v02 - ExFmean_v00) / ExFsd_v00)


#Attention
##baseline
attentionmean_v00 <- mean(PMplus.ok$attention_v00, na.rm = TRUE)
attentionsd_v00 <- sd(PMplus.ok$attention_v00, na.rm = TRUE)
PMplus.ok <- PMplus.ok %>%
  mutate(zattention_v00 = (attention_v00 - attentionmean_v00) / attentionsd_v00)

##2-year
PMplus.ok <- PMplus.ok %>%
  mutate(zattention_v02 = (attention_v02 - attentionmean_v00) / attentionsd_v00)


#Language
##baseline
languagemean_v00 <- mean(PMplus.ok$language_v00, na.rm = TRUE)
languagesd_v00 <- sd(PMplus.ok$language_v00, na.rm = TRUE)
PMplus.ok <- PMplus.ok %>%
  mutate(zlanguage_v00 = (language_v00 - languagemean_v00) / languagesd_v00)
##2-year
PMplus.ok <- PMplus.ok %>%
  mutate(zlanguage_v02 = (language_v02 - languagemean_v00) / languagesd_v00)

# Generate 2-YEAR z-score changes
PMplus.ok <- PMplus.ok %>%
  mutate(zMMSE_2yc = zMMSE_v02 - zMMSE_v00)
PMplus.ok <- PMplus.ok %>%
  mutate(zCDT_2yc = zCDT_v02 - zCDT_v00)
PMplus.ok <- PMplus.ok %>%
  mutate(zVFTa_2yc = zVFTa_v02 - zVFTa_v00)
PMplus.ok <- PMplus.ok %>%
  mutate(zVFTp_2yc = zVFTp_v02 - zVFTp_v00)
PMplus.ok <- PMplus.ok %>%
  mutate(zTMTa_2yc = zTMTa_v02 - zTMTa_v00)
PMplus.ok <- PMplus.ok %>%
  mutate(zTMTb_2yc = zTMTb_v02 - zTMTb_v00)
PMplus.ok <- PMplus.ok %>%
  mutate(zDSTf_2yc = zDSTf_v02 - zDSTf_v00)
PMplus.ok <- PMplus.ok %>%
  mutate(zDSTb_2yc = zDSTb_v02 - zDSTb_v00)
PMplus.ok <- PMplus.ok %>%
  mutate(zGCF_2yc = zGCF_v02 - zGCF_v00)
PMplus.ok <- PMplus.ok %>%
  mutate(zgenCF_2yc = zgenCF_v02 - zgenCF_v00)
PMplus.ok <- PMplus.ok %>%
  mutate(zExF_2yc = zExF_v02 - zExF_v00)
PMplus.ok <- PMplus.ok %>%
  mutate(zattention_2yc = zattention_v02 - zattention_v00)
PMplus.ok <- PMplus.ok %>%
  mutate(zlanguage_2yc = zlanguage_v02 - zlanguage_v00)
# * * Covariate -------------------------------------------------
#Socio-demographic COVARIATES
#Age
class(PMplus.ok$edad_s1)
summary(PMplus.ok$edad_s1) 

#Sex
class(PMplus.ok$sexo_s1)
PMplus.ok$sexo_s1 <- as.factor(PMplus.ok$sexo_s1)
table(PMplus.ok$sexo_s1)

#Intervention group //***TAKE INTO ACCOUNT IN LONGITUDINAL STUDIES***
class(PMplus.ok$grupo_int_v00)
PMplus.ok$grupo_int_v00 <- as.factor(PMplus.ok$grupo_int_v00)
table(PMplus.ok$grupo_int_v00)

#create area 
PMplus.ok <- PMplus.ok %>%
  mutate(area = nodo)
PMplus.ok$area[PMplus.ok$nodo %in% c("6","7")] <-0
PMplus.ok$area[PMplus.ok$nodo %in% c("12","21")] <-1
class(PMplus.ok$area)
PMplus.ok$area <- as.factor(PMplus.ok$area)
PMplus.ok$area <- factor(PMplus.ok$area,
                              levels = c(0,1),
                              labels = c("Catalunya" , "Valencia"))
table(PMplus.ok$area, useNA = "ifany")

##nodo for gut microbiota analysis 
class(PMplus.ok$nodo)
PMplus.ok$nodo <- as.factor(PMplus.ok$nodo)
table(PMplus.ok$nodo)

#Education
class(PMplus.ok$escola_v00)
table(PMplus.ok$escola_v00, useNA = "ifany")
PMplus.ok <- PMplus.ok %>%
  mutate(edu_level_v00 = escola_v00)
PMplus.ok$edu_level_v00[PMplus.ok$escola_v00 %in% c(4,5)] <-0
PMplus.ok$edu_level_v00[PMplus.ok$escola_v00 %in% c(3)] <-1
PMplus.ok$edu_level_v00[PMplus.ok$escola_v00 %in% c(1,2)] <-2
class(PMplus.ok$edu_level_v00)
PMplus.ok$edu_level_v00 <- as.factor(PMplus.ok$edu_level_v00)
PMplus.ok$edu_level_v00 <- factor(PMplus.ok$edu_level_v00,
                                  levels = c(0,1,2),
                                  labels = c("Primary school or less" , "High school" , "College"))
table(PMplus.ok$edu_level_v00, useNA = "ifany")


#Civil status
class(PMplus.ok$civil_v00)
table(PMplus.ok$civil_v00, useNA = "ifany")
PMplus.ok <- PMplus.ok %>%
  mutate(civil_status_v00 = civil_v00)
PMplus.ok$civil_status_v00[PMplus.ok$civil_v00 %in% c(1,4,5,6)] <-0
PMplus.ok$civil_status_v00[PMplus.ok$civil_v00 %in% c(2, NA)] <-1
PMplus.ok$civil_status_v00[PMplus.ok$civil_v00 %in% c(3)] <-2
class(PMplus.ok$civil_status_v00)
PMplus.ok$civil_status_v00 <- as.factor(PMplus.ok$civil_status_v00)
PMplus.ok$civil_status_v00 <- factor(PMplus.ok$civil_status_v00,
                                     levels = c(0,1,2),
                                     labels = c("Single, divorced or separated" , "Married" , "Widower"))
table(PMplus.ok$civil_status_v00, useNA = "ifany")

#Anthropometric variables & Personal history of illness (including Depressive symptomatology) & Lifestyle COVARIATES
#BMI
class(PMplus.ok$imc_v00)
summary(PMplus.ok$imc_v00)

#Waist circumference
PMplus.ok$waist_cir_v00 <- rowMeans(PMplus.ok[, c("cintura1_v00", "cintura2_v00")], na.rm = TRUE)
summary(PMplus.ok$waist_cir_v00)
class(PMplus.ok$waist_cir_v00)

#Diabetes status
class(PMplus.ok$diab_prev_s1)
PMplus.ok$diab_prev_s1 <- as.factor(PMplus.ok$diab_prev_s1)
table(PMplus.ok$diab_prev_s1, useNA = "ifany")

#Hypertension
class(PMplus.ok$hta_s1)
table(PMplus.ok$hta_s1, useNA = "ifany")
PMplus.ok <- PMplus.ok %>%
  mutate(hta_v00 = hta_s1)
PMplus.ok$hta_v00[PMplus.ok$hta_s1 %in% c(0)] <-0
PMplus.ok$hta_v00[PMplus.ok$hta_s1 %in% c(1,9)] <-1
class(PMplus.ok$hta_v00)
PMplus.ok$hta_v00 <- as.factor(PMplus.ok$hta_v00)
PMplus.ok$hta_v00 <- factor(PMplus.ok$hta_v00,
                            levels = c(0,1),
                            labels = c("No" , "Yes"))
table(PMplus.ok$hta_v00, useNA = "ifany")

#Hypercholesterolemia
class(PMplus.ok$colest_s1)
table(PMplus.ok$colest_s1, useNA = "ifany")
PMplus.ok <- PMplus.ok %>%
  mutate(colest_v00 = colest_s1)
PMplus.ok$colest_v00[PMplus.ok$colest_s1 %in% c(0)] <-0
PMplus.ok$colest_v00[PMplus.ok$colest_s1 %in% c(1,9)] <-1
class(PMplus.ok$colest_v00)
PMplus.ok$colest_v00 <- as.factor(PMplus.ok$colest_v00)
PMplus.ok$colest_v00 <- factor(PMplus.ok$colest_v00,
                               levels = c(0,1),
                               labels = c("No" , "Yes"))
table(PMplus.ok$colest_v00, useNA = "ifany")

#Depressive symptomatology
class(PMplus.ok$bdi_total_v00)
table(PMplus.ok$bdi_total_v00, useNA = "ifany")
PMplus.ok$depression_v00 <- NA
PMplus.ok$depression_v00[PMplus.ok$bdi_total_v00 >=14] <- 1
PMplus.ok$depression_v00[PMplus.ok$bdi_total_v00 <14 | is.na(PMplus.ok$bdi_total_v00)] <- 0
PMplus.ok$depression_v00 <- factor(PMplus.ok$depression_v00,
                                   levels = c(0,1),
                                   labels = c("No" , "Yes"))
table(PMplus.ok$depression_v00, useNA = "ifany")

#Medication use 
##Antihypertensive
PMplus.ok <- PMplus.ok %>%
  mutate(med.hypertensive_v00 = idr_v00+ araii_v00 +ieca_v00+ tazd_v00 +antagonistas_calcio_v00 +diur_asa_v00+ diur_potasio_v00+ betabloq_v00+ otros_hipotens_v00)
PMplus.ok$antihypertensive_v00 <- NA
PMplus.ok$antihypertensive_v00[PMplus.ok$med.hypertensive_v00 >=1] <- 1
PMplus.ok$antihypertensive_v00[PMplus.ok$med.hypertensive_v00 <1 ] <- 0
PMplus.ok$antihypertensive_v00 <- factor(PMplus.ok$antihypertensive_v00,
                                         levels = c(0,1),
                                         labels = c("No" , "Yes"))
table(PMplus.ok$antihypertensive_v00, useNA = "ifany")
table(PMplus.ok$tto_ta_v00, useNA = "ifany")
##Antidiabetic
PMplus.ok <- PMplus.ok %>%
  mutate(med.antidiabetic_v00 = otroshipo_noinsu_v00 +insulinas_v00+ metformina_v00+ otras_biguanidas_v00 +sulfonilureas_v00+ inhi_alfa_glucosidasa_v00+ tiazolinidionas_v00+ secretagogos_v00+ analogos_glp1_v00+ idpp4_v00 +islgt2_v00)
PMplus.ok$antidiabetic_v00 <- NA
PMplus.ok$antidiabetic_v00[PMplus.ok$med.antidiabetic_v00 >=1] <- 1
PMplus.ok$antidiabetic_v00[PMplus.ok$med.antidiabetic_v00 <1 ] <- 0
PMplus.ok$antidiabetic_v00 <- factor(PMplus.ok$antidiabetic_v00,
                                     levels = c(0,1),
                                     labels = c("No" , "Yes"))
table(PMplus.ok$antidiabetic_v00, useNA = "ifany")
table(PMplus.ok$tto_dm_v00, useNA = "ifany")
##Hypolipidemic
PMplus.ok <- PMplus.ok %>%
  mutate(med.hypolipidemic_v00 = estatina_v00 + otrhipo_v00)
PMplus.ok$hypolipidemic_v00 <- NA
PMplus.ok$hypolipidemic_v00[PMplus.ok$med.hypolipidemic_v00 >=1] <- 1
PMplus.ok$hypolipidemic_v00[PMplus.ok$med.hypolipidemic_v00 <1 ] <- 0
PMplus.ok$hypolipidemic_v00 <- factor(PMplus.ok$hypolipidemic_v00,
                                      levels = c(0,1),
                                      labels = c("No" , "Yes"))
table(PMplus.ok$hypolipidemic_v00, useNA = "ifany")
table(PMplus.ok$tto_col_v00, useNA = "ifany")

##Antidepressant
class(PMplus.ok$tto_tranq_v00)
table(PMplus.ok$tto_tranq_v00, useNA = "ifany")
PMplus.ok$antidepressant_v00 <- factor(PMplus.ok$tto_tranq_v00,
                                     levels = c(0,1),
                                     labels = c("No" , "Yes"))

table(PMplus.ok$antidepressant_v00, useNA = "ifany")

#Clinical data
summary(PMplus.ok[c("glucosa_v00", "coltot_v00", "trigli_v00", "hba1c_v00")])

#Physical activity
summary(PMplus.ok$geaf_tot_v00)
PMplus.ok <- PMplus.ok %>%
  mutate(exercise_day_v00 = geaf_tot_v00/7) #change it from per week to per day, use the per day value
summary(PMplus.ok$exercise_day_v00)

#Smoking status
class(PMplus.ok$fuma_s1)
table(PMplus.ok$fuma_s1, useNA = "ifany")
PMplus.ok <- PMplus.ok %>%
  mutate(smoke_status_v00 = fuma_s1)
PMplus.ok$smoke_status_v00[PMplus.ok$fuma_s1 %in% c(5,9)] <-0
PMplus.ok$smoke_status_v00[PMplus.ok$fuma_s1 %in% c(2,3,4)] <-1
PMplus.ok$smoke_status_v00[PMplus.ok$fuma_s1 %in% c(1)] <-2
class(PMplus.ok$smoke_status_v00)
PMplus.ok$smoke_status_v00 <- as.factor(PMplus.ok$smoke_status_v00)
PMplus.ok$smoke_status_v00 <- factor(PMplus.ok$smoke_status_v00,
                                     levels = c(0,1,2),
                                     labels = c("Never smoker" , "Former smoker" , "Smoker"))
table(PMplus.ok$smoke_status_v00, useNA = "ifany")

#Energy (caloric) intake, kcal/day
summary(PMplus.ok$energiat_v00)

#Alcohol consumption, g/day
summary(PMplus.ok$alcoholg_v00)
PMplus.ok$qalcoholg_v00 <-PMplus.ok$alcoholg_v00*PMplus.ok$alcoholg_v00

#Adherence to Med. diet (14-point) baseline:
#/*Mediterranean diet score:
#Adherence based on the 14 point Questionnaire for Med Diet Adherence used in the PREDIMED-Plus trial.
#Foods considered / Alimentos a considerar: 
#1)  1 point = Olive oil is the primary oil used in cooking.
##1. Same as 17 Point Q1 (p17_1_v00). Generation of score for the use of olive oil as primary oil.
class(PMplus.ok$p17_1_v00)
table(PMplus.ok$p17_1_v00)
#2)  1 point = >= 4 Tbsp olive oil / day.
##2. Generation of score for the amount of olive oil consumed.
PMplus.ok <- PMplus.ok %>%
  mutate(oliveoil_v00 = ac_oliva_v00 + ac_olivavir_v00)
PMplus.ok$p14_oliveoil_v00 <- NA
PMplus.ok$p14_oliveoil_v00[PMplus.ok$oliveoil_v00 >=40] <- 1
PMplus.ok$p14_oliveoil_v00[PMplus.ok$oliveoil_v00 <40] <- 0
class(PMplus.ok$p14_oliveoil_v00)
table(PMplus.ok$p14_oliveoil_v00, useNA = "ifany")
#PMplus.ok$p14_oliveoil_v00 <- factor(PMplus.ok$p14_oliveoil_v00,
#                                  levels = c(0,1),
#                                  labels = c("less than 4 Tbsp per day" , "at least 4 Tbsp per day"))
#table(PMplus.ok$p14_oliveoil_v00, useNA = "ifany")
#PMplus.ok$p14_oliveoil_v00<- as.numeric(PMplus.ok$p14_oliveoil_v00)

#3)  1 point = >= 2 servings of vegetables / day (1 servings = 200 g).
##3. Same as 17 Point Q2 (p17_2_v00).
class(PMplus.ok$p17_2_v00)
table(PMplus.ok$p17_2_v00)
#4)  1 point = >= 3 servings of fruits / day.
##4. Same as 17 Point Q3 (p17_3_v00).
class(PMplus.ok$p17_3_v00)
table(PMplus.ok$p17_3_v00)
#5)  1 point = < 1 serving of red meat / day (1 servings = 100 to 150 g).
PMplus.ok <- PMplus.ok %>%
  mutate(redmeat_v00 = c_ternera_v00 + c_cerdo_v00 + c_cordero_v00 + higad_v00 + visceras_v00 + j_serrano_v00 + j_cocido_v00 + embutidos_v00 + pates_v00 + hamburguesa_v00 + bacon_v00)
PMplus.ok$p14_redmeat_v00 <- NA
PMplus.ok$p14_redmeat_v00[PMplus.ok$redmeat_v00 >=150] <- 0
PMplus.ok$p14_redmeat_v00[PMplus.ok$redmeat_v00 <150] <- 1
#PMplus.ok$p14_redmeat_v00 <- factor(PMplus.ok$p14_redmeat_v00,
#                                  levels = c(0,1),
#                                  labels = c("at least 150 g per day" , "less than 150 g per day"))
#table(PMplus.ok$p14_redmeat_v00, useNA = "ifany")
#PMplus.ok$p14_redmeat_v00<- as.numeric(PMplus.ok$p14_redmeat_v00)

#6)  1 point = < 1 serving of butter/margarine / day (1 serving = 12 g).
PMplus.ok <- PMplus.ok %>%
  mutate(butter_v00 = mantequillas_v00 + margarinas_v00 + mantecacer_v00 + nata_crema_v00)
PMplus.ok$p14_butter_v00 <- NA
PMplus.ok$p14_butter_v00[PMplus.ok$butter_v00 >=12] <- 0
PMplus.ok$p14_butter_v00[PMplus.ok$butter_v00 <12] <- 1
#PMplus.ok$p14_butter_v00 <- factor(PMplus.ok$p14_butter_v00,
#                                  levels = c(0,1),
#                                  labels = c("at least 12 g per day" , "less than 12 g per day"))
#table(PMplus.ok$p14_butter_v00, useNA = "ifany")
#PMplus.ok$p14_butter_v00<-as.numeric(PMplus.ok$p14_butter_v00)

#7)  1 point = < 1 serving of SSBs / day.
PMplus.ok <- PMplus.ok %>%
  mutate(ssb_v00 = refrescos_v00 + z_botella_v00)
PMplus.ok$p14_ssb_v00 <- NA
PMplus.ok$p14_ssb_v00[PMplus.ok$ssb_v00 >=200] <- 0
PMplus.ok$p14_ssb_v00[PMplus.ok$ssb_v00 <200] <- 1
#PMplus.ok$p14_ssb_v00 <- factor(PMplus.ok$p14_ssb_v00,
#                                  levels = c(0,1),
#                                  labels = c("at least 200 g per day" , "less than 200 g per day"))
#table(PMplus.ok$p14_ssb_v00, useNA = "ifany")
#PMplus.ok$p14_ssb_v00 <- as.numeric(PMplus.ok$p14_ssb_v00)

#8)  1 point = >= 7 glasses of wine / week.
PMplus.ok <- PMplus.ok %>%
  mutate(wine_v00 = (v_tintojov_v00/100 + v_tintoanej_v00/100)*7)
PMplus.ok$p14_wine_v00 <- NA
PMplus.ok$p14_wine_v00[PMplus.ok$wine_v00 >=7] <- 1
PMplus.ok$p14_wine_v00[PMplus.ok$wine_v00 <7] <- 0
#PMplus.ok$p14_wine_v00 <- factor(PMplus.ok$p14_wine_v00,
#                                  levels = c(0,1),
#                                  labels = c("less than 7 glasses per week" , "at least 7 glasses per week"))
#table(PMplus.ok$p14_wine_v00, useNA = "ifany")
#PMplus.ok$p14_wine_v00 <- as.numeric(PMplus.ok$p14_wine_v00)

#9)  1 point = >= 3 servings of legumes / week (1 serving = 150 g).
##9. Same as 17 Point Q7 (p17_7_v00).
class(PMplus.ok$p17_7_v00)
table(PMplus.ok$p17_7_v00)
#10) 1 point = >= 3 servings of fish/shellfish / week (1 serving = 100 to 150 g of fish or 4 to 5 pieces of shellfish).
##10. Same as 17 Point Q8 (p17_8_v00).
class(PMplus.ok$p17_8_v00)
table(PMplus.ok$p17_8_v00)
#11) 1 point = < 2 servings of commercial baked goods / week.
##11.Generation of score for commercial baked goods (sweets)
PMplus.ok <- PMplus.ok %>%
  mutate(sweets_v00 = (bizcocho_v00/50 + croissant_v00/50 + donut_v00/50 + magdalena_v00/40 + pastel_v00/50 + churro_v00/100 + mazapan_v00/90)*7)
PMplus.ok$p14_sweets_v00 <- NA
PMplus.ok$p14_sweets_v00[PMplus.ok$sweets_v00 >=2] <- 0
PMplus.ok$p14_sweets_v00[PMplus.ok$sweets_v00 <2] <- 1
#PMplus.ok$p14_sweets_v00 <- factor(PMplus.ok$p14_sweets_v00,
#                                  levels = c(0,1),
#                                  labels = c("at least 2 times per week" , "less than 2 times per week"))
#table(PMplus.ok$p14_sweets_v00, useNA = "ifany")
#PMplus.ok$p14_sweets_v00 <- as.numeric(PMplus.ok$p14_sweets_v00)

#12) 1 point = >= 3 servings of nuts / week (1 serving = 30 g).
##12. Generation of score for nuts
PMplus.ok <- PMplus.ok %>%
  mutate(nuts_v00 = (almendras_v00 + pistachos_v00 + nuez_v00 + f_secos_v00)/30*7)
PMplus.ok$p14_nuts_v00 <- NA
PMplus.ok$p14_nuts_v00[PMplus.ok$nuts_v00 >=3] <- 1
PMplus.ok$p14_nuts_v00[PMplus.ok$nuts_v00 <3] <- 0
#PMplus.ok$p14_nuts_v00 <- factor(PMplus.ok$p14_nuts_v00,
#                                  levels = c(0,1),
#                                  labels = c("less than 3 servings per week" , "at least 3 servings per week"))
#table(PMplus.ok$p14_nuts_v00, useNA = "ifany")
#PMplus.ok$p14_nuts_v00 <- as.numeric(PMplus.ok$p14_nuts_v00)

#13) 1 point = White meat is preferred over red meat.
##13. Same as 17 Point Q11 (p17_11_v00).
class(PMplus.ok$p17_11_v00)
table(PMplus.ok$p17_11_v00)

#14) 1 point = >= 2 servings a week of cooked vegetables, pasta, rice or other dishes seasoned with tomato sauce, garlic, onion or leek cooked over low heat with olive oil (sofrito).
##14. Same as 17 Point Q12 (p17_12_v00).
class(PMplus.ok$p17_12_v00)
table(PMplus.ok$p17_12_v00)
#*/

#*Creation of MED Diet 14-Score
PMplus.ok <- PMplus.ok %>%
  mutate(MED14score_v00 = p17_1_v00 + p14_oliveoil_v00 + p17_2_v00 + p17_3_v00 + p14_redmeat_v00 + p14_butter_v00 + p14_ssb_v00 + p14_wine_v00 + p17_7_v00 + p17_8_v00 + p14_sweets_v00 + p14_nuts_v00 + p17_11_v00 + p17_12_v00)

class(PMplus.ok$MED14score_v00)
summary(PMplus.ok$MED14score_v00)


class(PMplus.ok$p14_total_v00)
summary(PMplus.ok$p14_total_v00)

#*Creation of Adherence to Med. diet (12-point) without olive oil
PMplus.ok <- PMplus.ok %>%
  mutate(MED12score_v00 = p17_2_v00 + p17_3_v00 + p14_redmeat_v00 + p14_butter_v00 + p14_ssb_v00 + p14_wine_v00 + p17_7_v00 + p17_8_v00 + p14_sweets_v00 + p14_nuts_v00 + p17_11_v00 + p17_12_v00)

class(PMplus.ok$MED12score_v00)
summary(PMplus.ok$MED12score_v00)

cov.list <- c("exercise_day_v00", "alcoholg_v00", "MED12score_v00")
for (var in cov.list) {
  PMplus.ok <- PMplus.ok %>%
    mutate(!!paste0("ter_", var) := ntile(!!sym(var), 3))
}

PMplus.ok$ter_exercise_day_v00 <- as.factor(PMplus.ok$ter_exercise_day_v00)
PMplus.ok$ter_alcoholg_v00 <- as.factor(PMplus.ok$ter_alcoholg_v00)
PMplus.ok$ter_MED12score_v00 <- as.factor(PMplus.ok$ter_MED12score_v00)
# Dietary variables 
# Energy adjusted using the residual methods
var.list <- c("hc_v00", "prot_v00", "gratot_v00", "mo_v00", "po_v00", "sa_v00", "porc_hc_v00", "porc_pr_v00", "porc_gr_v00", "porc_mo_v00", "porc_po_v00", "porc_sa_v00", 
              "fibra_v00", "verdutot_v00", "frutatot_v00", "legumbre_v00", 
              "cereal_v00", "lacteos_v00", "carnicos_v00", "pescados_v00", "fsecos_v00", "gallet_v00") # List of variables that need to be adjusted by energy intake 
for (var in var.list) { # perform the loop over the variables 
  # run the regression 
  model<-lm(paste(var, "energiat_v00", sep = "~"), data = PMplus.ok )
  # get the residuals from the regression 
  residuals <- residuals(model)
  # calculate the mean energy intake 
  energy_mean <- mean(PMplus.ok$energiat_v00, na.rm = TRUE)
  # expected food intake for a person with the mean energy intake
  expected_var <- coef(model)[1] + coef(model)[2] * energy_mean 
  # create the new variable "ajust_var" and add the residuals and expected food intake
  PMplus.ok[paste0("ajust_", var)] <-expected_var + residuals 
}
# * * Keep variables of interests -------------------------------------------------
metadata_final <- dplyr::select (PMplus.ok, paciente, SampleIdentifier, SampleID, OriginalID, SampleIDs, SampleIDcounts, Timepoint_hr, 
                                 olivatot_v00, ac_olivavir_v00,  re_oliva_v00, otheroil_v00, ajust_olivatot_v00, ajust_ac_olivavir_v00, ajust_re_oliva_v00, ajust_otheroil_v00, ajust_tspolivatot_v00, ajust_tspac_olivavir_v00, ajust_tspre_oliva_v00, ajust_tspotheroil_v00,
                                 ter_ajust_olivatot_v00, ter_ajust_ac_olivavir_v00, ter_ajust_re_oliva_v00,
                                 ptrend_ter_ajust_olivatot_v00, ptrend_ter_ajust_ac_olivavir_v00, ptrend_ter_ajust_re_oliva_v00, 
                                 MMSE_v00, MMSE_v02, CDT_v00, CDT_v02, VFTa_v00, VFTa_v02, VFTp_v00, VFTp_v02, TMTa_v00, TMTa_v02, 
                                 TMTb_v00, TMTb_v02, DSTf_v00, DSTf_v02, DSTb_v00, DSTb_v02, 
                                 zMMSE_v00, zMMSE_v02, zCDT_v00, zCDT_v02, zVFTa_v00, zVFTa_v02, zVFTp_v00, zVFTp_v02, 
                                 zTMTa_v00, zTMTa_v02, zTMTb_v00, zTMTb_v02, zDSTf_v00, zDSTf_v02, zDSTb_v00, zDSTb_v02, 
                                 GCF_v00, GCF_v02, genCF_v00, genCF_v02, ExF_v00, ExF_v02, attention_v00, attention_v02, language_v00, language_v02, 
                                 zGCF_v00, zGCF_v02, zgenCF_v00, zgenCF_v02, zExF_v00, zExF_v02, zattention_v00, zattention_v02, zlanguage_v00, zlanguage_v02, 
                                 zMMSE_2yc, zCDT_2yc, zVFTa_2yc, zVFTp_2yc, zTMTa_2yc, zTMTb_2yc, zDSTf_2yc, zDSTb_2yc, zGCF_2yc, zgenCF_2yc, zExF_2yc, zattention_2yc, zlanguage_2yc, 
                                 sexo_s1, edad_s1,  nodo, area, grupo_int_v00, edu_level_v00, civil_status_v00, 
                                 imc_v00,  waist_cir_v00, diab_prev_s1, hta_v00, colest_v00, depression_v00,  
                                 antihypertensive_v00, antidiabetic_v00, hypolipidemic_v00, antidepressant_v00, 
                                 exercise_day_v00, ter_exercise_day_v00, smoke_status_v00, MED14score_v00, MED12score_v00, ter_MED12score_v00, 
                                 energiat_v00, alcoholg_v00, ter_alcoholg_v00, qalcoholg_v00, 
                                 hc_v00, prot_v00, gratot_v00, mo_v00, po_v00, sa_v00, porc_hc_v00, porc_pr_v00, porc_gr_v00, porc_mo_v00, porc_po_v00, porc_sa_v00, 
                                 fibra_v00, col_v00, fit_v00, trans_v00, linoleico_v00, linolenico_v00, omega3_v00, n3marinos_v00, verdutot_v00, frutatot_v00, legumbre_v00, 
                                 cereal_v00, lacteos_v00, carnicos_v00, pescados_v00, fsecos_v00, gallet_v00, 
                                 ajust_hc_v00, ajust_prot_v00, ajust_gratot_v00, ajust_mo_v00, ajust_po_v00, ajust_sa_v00, ajust_porc_hc_v00, ajust_porc_pr_v00, ajust_porc_gr_v00, ajust_porc_mo_v00, ajust_porc_po_v00, ajust_porc_sa_v00, 
                                 ajust_fibra_v00, ajust_verdutot_v00, ajust_frutatot_v00, ajust_legumbre_v00, 
                                 ajust_cereal_v00, ajust_lacteos_v00, ajust_carnicos_v00, ajust_pescados_v00, ajust_fsecos_v00, ajust_gallet_v00) 
# SAVE FINAL METADATA  -----------------------------------------------
#SAVE THE SAMPLE DATA IN .rds FILE 
rownames(metadata_final)<- metadata_final$SampleIDcounts
saveRDS(metadata_final, file = "DATA/metadata_FINAL_110324.rds")
write.csv(metadata_final, file = "DATA/metadata.csv")

# SAVE FINAL Phyloseq object  -----------------------------------------------
sample_data(ps) <- metadata_final #update the sample_data component of phyloseq object with the merged meta data 
saveRDS(ps, file = "DATA/jiaqi_Phyloseq_110324.rds")

# SAVE OTU table  -----------------------------------------------
otudata_final <- as.data.frame(otu_table(ps))
write.csv(otudata_final, file = "DATA/ASV_table.csv")
# SAVE Taxonomic table  -----------------------------------------------
taxadata_final <- as.data.frame(tax_table(ps))
write.csv(taxadata_final, file = "DATA/taxonomic_table.csv")


## Clear the global environment
if (sys.nframe() != 0L) {
  rm(list = setdiff(ls(), c(freeze))) # remove objects created by the script when sourced from main.R to preserve memory
  invisible(gc()) # Do a silent garbage collect
}

