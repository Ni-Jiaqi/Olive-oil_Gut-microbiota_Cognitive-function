# Title: Cross-sectional associations of olive oil consumption with human gut microbiota diversity and composition
# Author: Jiaqi Ni [ORCID](https://orcid.org/0000-0001-5625-0373)
# Description: This is the analysis script for differential abundance analysis 

########################################################### STARTUP ################################################################
# Set up working environment  ------------------------------------------------------------------
rm(list=ls())
options(max.print=1000000)
options(scipen = 999)

getwd()
setwd("/home/jiaqi/Jiaqi_Ni/OliveMicrobiomeCognition_Microbiome/Output/maaslin/")
if (sys.nframe() == 0L) {
  this.dir <- "/home/jiaqi/Jiaqi_Ni/OliveMicrobiomeCognition_Microbiome/"
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
library(grid)
library(gridExtra)
library(ggExtra)
library(ggside)
library(cowplot)
library(sjPlot)
library(Maaslin2)
library(ggalluvial)
library(tidyverse)
library(modelsummary)
require(phyloseq)
# Data preparation ---------------------------------------------------

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
# * extract data for dataframe ---------------------------------------------------
taxa_table_g <- as.data.frame(t(assays(tse_genus)$clr))
# Meta data
meta_table <- as.data.frame(colData(tse))
# Define and fix continuous and categorical variables
catVars <- c("sexo_s1",  "nodo", "area", "grupo_int_v00", "edu_level_v00", "civil_status_v00",
             "ter_ajust_olivatot_v00", "ter_ajust_ac_olivavir_v00", "ter_ajust_re_oliva_v00", 
              "diab_prev_s1", "hta_v00", "colest_v00", "depression_v00",  "antihypertensive_v00", "antidepressant_v00", 
             "antidiabetic_v00", "hypolipidemic_v00", "ter_exercise_day_v00", 
             "smoke_status_v00","ter_MED12score_v00","ter_alcoholg_v00")
meta_table[,catVars] <- lapply(meta_table[,catVars], as.factor)

conVars <- c("edad_s1", "olivatot_v00",  "ac_olivavir_v00", "re_oliva_v00",  
             "ajust_olivatot_v00", "ajust_ac_olivavir_v00", "ajust_re_oliva_v00", 
             "ptrend_ter_ajust_olivatot_v00", "ptrend_ter_ajust_ac_olivavir_v00", "ptrend_ter_ajust_re_oliva_v00",
             "imc_v00", "waist_cir_v00", "exercise_day_v00", "alcoholg_v00", "qalcoholg_v00", "MED12score_v00")
meta_table[,conVars] <- lapply(meta_table[,conVars], as.numeric)
########################################################### ANALYSIS ################################################################
#  MAASLIN 2 ---------------------------------------------------
# Define function for easier running
# Maaslin will be run with the following parameters:
# - Fixed effects:  olive oil tertiles, baseline age, sex, geographical area (Catalunya, Valencia), 
#                  educational level (primary, secondary, or college), civil status (single, divorced or separated, married, widower), 
#                  BMI (kg/m2), physical activity (METs/min/day), smoking status (current, former, or never), alcohol consumption in g/day (and adding the quadratic term), 
#                  depressive symptomatology (yes/no), diabetes prevalence (yes/no), hypertension prevalence (yes/no), and hypercholesterolemia prevalence (yes/no), 
#                  adherence to Mediterranean diet 12 scores
# - Filtering: minimum 10 % prevalence and 0.1 % abundance (NIH)
# - Transformation: CLR transformed relative abundances
# - Standardization of continuous variables: NO
# - FDR correction: Benjamini-Hochberg (q = 0.05 significance threshold)

# - significant taxa (q<0.25) ------
run_maaslin_full <- function(taxa.table, meta.table, fixed.effects, folder.name, reference.group) {
  fit.data <- Maaslin2(
    input_data = taxa.table, # Define the taxa data frame
    input_metadata = meta.table, # Define the phenotype data frame from where to pull covariates and nutritional data
    fixed_effects = fixed.effects, # Define the variables to use on the right-hand-side of the model formula
    output = folder.name, # Define the name for the output folder
    min_abundance = -Inf, # Don't set any filtering thresholds for abundance or prevalence (we have already done this manually above)
    min_prevalence = -Inf,
    normalization = "NONE", # Don't do any normalization (we have done this already manually above with CLR)
    transform = "NONE", # Don't make any transformations
    standardize = F, # T: Convert continuous variables into Z-scores (same as "scale()")
    analysis_method = "LM", # Use linear regression models
    correction = "BH", # Use Benjamini-Hochberg false-discovery-rate correction
    reference = reference.group,
    max_significance = 0.25, # Set the significance threshold for the FDR-corrected p-values
    plot_heatmap = FALSE, # Don't plot heatmaps
    plot_scatter = FALSE) # Don't plot scatter plots
}

# Run the analyses
run_maaslin_full(taxa_table_g, meta_table, fixed.effects = c("ajust_tspolivatot_v00","edad_s1","sexo_s1"), folder.name = "Basic_total", reference.group = "NULL")
run_maaslin_full(taxa_table_g, meta_table, fixed.effects = c("ajust_tspolivatot_v00","edad_s1","sexo_s1", "area", "edu_level_v00" , "civil_status_v00"), folder.name = "Sociodemo_total", reference.group = "NULL")
run_maaslin_full(taxa_table_g, meta_table, fixed.effects = c("ajust_tspolivatot_v00","edad_s1","sexo_s1", "area", "edu_level_v00" , "civil_status_v00",
                                                             "imc_v00" , "exercise_day_v00", "smoke_status_v00" , "alcoholg_v00" , "qalcoholg_v00" , 
                                                             "depression_v00" , "diab_prev_s1" , "hta_v00" , "colest_v00" ,
                                                             "MED12score_v00"), folder.name = "Full_total", reference.group = "NULL")

run_maaslin_full(taxa_table_g, meta_table, fixed.effects = c("ajust_tspac_olivavir_v00","edad_s1","sexo_s1"), folder.name = "Basic_ex", reference.group = "NULL")
run_maaslin_full(taxa_table_g, meta_table, fixed.effects = c("ajust_tspac_olivavir_v00","edad_s1","sexo_s1", "area", "edu_level_v00" , "civil_status_v00"), folder.name = "Sociodemo_ex", reference.group = "NULL")
run_maaslin_full(taxa_table_g, meta_table, fixed.effects = c("ajust_tspac_olivavir_v00","edad_s1","sexo_s1", "area", "edu_level_v00" , "civil_status_v00",
                                                             "imc_v00" , "exercise_day_v00", "smoke_status_v00" , "alcoholg_v00" , "qalcoholg_v00" , 
                                                             "depression_v00" , "diab_prev_s1" , "hta_v00" , "colest_v00" ,
                                                             "MED12score_v00"), folder.name = "Full_ex", reference.group = "NULL")

run_maaslin_full(taxa_table_g, meta_table, fixed.effects = c("ajust_tspre_oliva_v00","edad_s1","sexo_s1"), folder.name = "Basic_re", reference.group = "NULL")
run_maaslin_full(taxa_table_g, meta_table, fixed.effects = c("ajust_tspre_oliva_v00","edad_s1","sexo_s1", "area", "edu_level_v00" , "civil_status_v00"), folder.name = "Sociodemo_re", reference.group = "NULL")
run_maaslin_full(taxa_table_g, meta_table, fixed.effects = c("ajust_tspre_oliva_v00","edad_s1","sexo_s1", "area", "edu_level_v00" , "civil_status_v00",
                                                             "imc_v00" , "exercise_day_v00", "smoke_status_v00" , "alcoholg_v00" , "qalcoholg_v00" , 
                                                             "depression_v00" , "diab_prev_s1" , "hta_v00" , "colest_v00" ,
                                                             "MED12score_v00"), folder.name = "Full_re", reference.group = "NULL")

# Load results into memory
g_results_1 <- read.table(paste0(this.dir,"Output/maaslin/Basic_total/significant_results.tsv"), sep = "\t", header = TRUE)
g_results_2 <- read.table(paste0(this.dir,"Output/maaslin/Sociodemo_total/significant_results.tsv"), sep = "\t", header = TRUE)
g_results_3 <- read.table(paste0(this.dir,"Output/maaslin/Full_total/significant_results.tsv"), sep = "\t", header = TRUE)
g_results_4 <- read.table(paste0(this.dir,"Output/maaslin/Basic_ex/significant_results.tsv"), sep = "\t", header = TRUE)
g_results_5 <- read.table(paste0(this.dir,"Output/maaslin/Sociodemo_ex/significant_results.tsv"), sep = "\t", header = TRUE)
g_results_6 <- read.table(paste0(this.dir,"Output/maaslin/Full_ex/significant_results.tsv"), sep = "\t", header = TRUE)
g_results_7 <- read.table(paste0(this.dir,"Output/maaslin/Basic_re/significant_results.tsv"), sep = "\t", header = TRUE)
g_results_8 <- read.table(paste0(this.dir,"Output/maaslin/Sociodemo_re/significant_results.tsv"), sep = "\t", header = TRUE)
g_results_9 <- read.table(paste0(this.dir,"Output/maaslin/Full_re/significant_results.tsv"), sep = "\t", header = TRUE)


# combine results from maasline regressions across different dietary variables
sig_results<- function(dir, exposure){
  result<-read.table(file = dir, sep = '\t', header = TRUE, check.names=FALSE)
  sig_result <- subset(result, metadata==exposure, select =c(feature, metadata))
  return(sig_result)
}
sigtoo <- sig_results(dir= paste0(this.dir,"Output/maaslin/Full_total/significant_results.tsv"),   exposure = "ajust_tspolivatot_v00")
sigvoo <- sig_results(dir= paste0(this.dir,"Output/maaslin/Full_ex/significant_results.tsv"),     exposure = "ajust_tspac_olivavir_v00")
sigcoo <- sig_results(dir= paste0(this.dir,"Output/maaslin/Full_re/significant_results.tsv"),   exposure = "ajust_tspre_oliva_v00")



# - all taxa  ------
run_maaslin_full <- function(taxa.table, meta.table, fixed.effects, folder.name, reference.group) {
  fit.data <- Maaslin2(
    input_data = taxa.table, # Define the taxa data frame
    input_metadata = meta.table, # Define the phenotype data frame from where to pull covariates and nutritional data
    fixed_effects = fixed.effects, # Define the variables to use on the right-hand-side of the model formula
    output = folder.name, # Define the name for the output folder
    min_abundance = -Inf, # Don't set any filtering thresholds for abundance or prevalence (we have already done this manually above)
    min_prevalence = -Inf,
    normalization = "NONE", # Don't do any normalization (we have done this already manually above with CLR)
    transform = "NONE", # Don't make any transformations
    standardize = F, # T: Convert continuous variables into Z-scores (same as "scale()")
    analysis_method = "LM", # Use linear regression models
    correction = "BH", # Use Benjamini-Hochberg false-discovery-rate correction
    reference = reference.group,
    max_significance = 1, # Set the significance threshold for the FDR-corrected p-values
    plot_heatmap = FALSE, # Don't plot heatmaps
    plot_scatter = FALSE) # Don't plot scatter plots
}

# Run the analyses
run_maaslin_full(taxa_table_g, meta_table, fixed.effects = c("ajust_tspolivatot_v00","edad_s1","sexo_s1"), folder.name = "Basic_total", reference.group = "NULL")
run_maaslin_full(taxa_table_g, meta_table, fixed.effects = c("ajust_tspolivatot_v00","edad_s1","sexo_s1", "area", "edu_level_v00" , "civil_status_v00"), folder.name = "Sociodemo_total", reference.group = "NULL")
run_maaslin_full(taxa_table_g, meta_table, fixed.effects = c("ajust_tspolivatot_v00","edad_s1","sexo_s1", "area", "edu_level_v00" , "civil_status_v00",
                                                             "imc_v00" , "exercise_day_v00", "smoke_status_v00" , "alcoholg_v00" , "qalcoholg_v00" , 
                                                             "depression_v00" , "diab_prev_s1" , "hta_v00" , "colest_v00" ,
                                                             "MED12score_v00"), folder.name = "Full_total", reference.group = "NULL")

run_maaslin_full(taxa_table_g, meta_table, fixed.effects = c("ajust_tspac_olivavir_v00","edad_s1","sexo_s1"), folder.name = "Basic_ex", reference.group = "NULL")
run_maaslin_full(taxa_table_g, meta_table, fixed.effects = c("ajust_tspac_olivavir_v00","edad_s1","sexo_s1", "area", "edu_level_v00" , "civil_status_v00"), folder.name = "Sociodemo_ex", reference.group = "NULL")
run_maaslin_full(taxa_table_g, meta_table, fixed.effects = c("ajust_tspac_olivavir_v00","edad_s1","sexo_s1", "area", "edu_level_v00" , "civil_status_v00",
                                                             "imc_v00" , "exercise_day_v00", "smoke_status_v00" , "alcoholg_v00" , "qalcoholg_v00" , 
                                                             "depression_v00" , "diab_prev_s1" , "hta_v00" , "colest_v00" ,
                                                             "MED12score_v00"), folder.name = "Full_ex", reference.group = "NULL")

run_maaslin_full(taxa_table_g, meta_table, fixed.effects = c("ajust_tspre_oliva_v00","edad_s1","sexo_s1"), folder.name = "Basic_re", reference.group = "NULL")
run_maaslin_full(taxa_table_g, meta_table, fixed.effects = c("ajust_tspre_oliva_v00","edad_s1","sexo_s1", "area", "edu_level_v00" , "civil_status_v00"), folder.name = "Sociodemo_re", reference.group = "NULL")
run_maaslin_full(taxa_table_g, meta_table, fixed.effects = c("ajust_tspre_oliva_v00","edad_s1","sexo_s1", "area", "edu_level_v00" , "civil_status_v00",
                                                             "imc_v00" , "exercise_day_v00", "smoke_status_v00" , "alcoholg_v00" , "qalcoholg_v00" , 
                                                             "depression_v00" , "diab_prev_s1" , "hta_v00" , "colest_v00" ,
                                                             "MED12score_v00"), folder.name = "Full_re", reference.group = "NULL")

# Load results into memory
g_results_1 <- read.table(paste0(this.dir,"Output/maaslin/Basic_total/all_results.tsv"), sep = "\t", header = TRUE)
g_results_2 <- read.table(paste0(this.dir,"Output/maaslin/Sociodemo_total/all_results.tsv"), sep = "\t", header = TRUE)
g_results_3 <- read.table(paste0(this.dir,"Output/maaslin/Full_total/all_results.tsv"), sep = "\t", header = TRUE)
g_results_4 <- read.table(paste0(this.dir,"Output/maaslin/Basic_ex/all_results.tsv"), sep = "\t", header = TRUE)
g_results_5 <- read.table(paste0(this.dir,"Output/maaslin/Sociodemo_ex/all_results.tsv"), sep = "\t", header = TRUE)
g_results_6 <- read.table(paste0(this.dir,"Output/maaslin/Full_ex/all_results.tsv"), sep = "\t", header = TRUE)
g_results_7 <- read.table(paste0(this.dir,"Output/maaslin/Basic_re/all_results.tsv"), sep = "\t", header = TRUE)
g_results_8 <- read.table(paste0(this.dir,"Output/maaslin/Sociodemo_re/all_results.tsv"), sep = "\t", header = TRUE)
g_results_9 <- read.table(paste0(this.dir,"Output/maaslin/Full_re/all_results.tsv"), sep = "\t", header = TRUE)


all_results <- function(dir, exposure, label){
  result<-read.table(file =dir, sep = '\t', header = TRUE, check.names=FALSE)
  all_result <- subset(result, metadata==exposure, select =-c(value, N))
  all_result$meta <- label
  return(all_result)
}



sigtooall <- all_results(dir=paste0(this.dir,"Output/maaslin/Full_total/all_results.tsv"),   exposure = "ajust_tspolivatot_v00", label="TOO")
sigvooall <- all_results(dir=paste0(this.dir,"Output/maaslin/Full_ex/all_results.tsv"),     exposure = "ajust_tspac_olivavir_v00", label="VOO")
sigcooall <- all_results(dir=paste0(this.dir,"Output/maaslin/Full_re/all_results.tsv"),   exposure = "ajust_tspre_oliva_v00", label="COO")

join1<-full_join(sigtoo, sigvoo, by="feature")
join2<-full_join(join1, sigcoo, by="feature")

sigfeature<-subset(join2, select=feature)

join1<-left_join(sigfeature, sigtooall, by="feature")
join2<-left_join(sigfeature, sigvooall, by="feature")
join3<-left_join(sigfeature, sigcooall, by="feature")


bind4<-rbind(join1, join2, join3)
bind4$FeatureID <- sub("^Genus\\.", "Genus\\:", bind4$feature)
bind4$FeatureID <- sub("Genus:CAG.56", "Genus:CAG-56", bind4$FeatureID)
bind4$FeatureID <- sub("Genus:Christensenellaceae_R.7_group", "Genus:Christensenellaceae_R-7_group", bind4$FeatureID)

# read in table of species matched to phyla

phylum <- mia::meltAssay(tse_genus,
                         add_row_data = T,
                         add_col_data = T,
                         assay.type = "clr")


phylum <- phylum[, c("FeatureID", "Phylum")]
phylum_S<- phylum[!duplicated(phylum[, c("FeatureID", "Phylum")]), ]

# format specie & phylum names
bind4<-left_join(bind4, phylum_S, by="FeatureID")
bind4$feature <- substring(bind4$feature, 7)
bind4<-bind4[rev(order(bind4$meta, bind4$Phylum, bind4$feature)),]
bind4$name<-gsub("[^[:alnum:][:blank:]+?&/\\-]", " ", bind4$feature)
bug_names<-bind4[1:(nrow(bind4)/3), (ncol(bind4)-1):ncol(bind4)]

# create heatmap for diet-taxonomy associations
level_x_order <- factor(bind4$meta, level = c('TOO', 'VOO', 'COO'))
level_y_order <- factor(bind4$name, levels = bug_names$name)
bind4$stars <- cut(bind4$qval, breaks=c(-Inf, 0.01, 0.05, 0.1, 0.25, Inf), label=c("****", "***", "**", "*", ""))


p1<-ggplot(bind4, aes(level_x_order, level_y_order)) +
  geom_tile(aes(fill = coef), color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", space = "Lab" ) +
  geom_text(aes(label=stars), color="black", size=35, show.legend = TRUE) +
  xlab(NULL)  +
  theme(legend.title = element_text(size = 50, family = "Times New Roman"),
        legend.text = element_text(size = 50, family = "Times New Roman"),
        legend.position = "left",
        plot.title = element_text(size=50, family = "Times New Roman"),
        axis.title=element_text(size=50,face="bold", family = "Times New Roman"),
        axis.title.y= element_blank(),
        axis.text.y=element_text(size=50, face="bold.italic", family = "Times New Roman"),
        axis.text.x = element_text(angle = 0, hjust = 1, size=50, face="bold", family = "Times New Roman")) +
  labs(fill = "β coefficient") +
  guides(fill = guide_colorbar(barwidth = 1.5,
                               barheight = 10))
# create left bar for categorizing species to phyla
bind4$category<-as.factor(bind4$Phylum)
p2<-ggplot(bind4, aes(level_x_order, level_y_order))+
  scale_y_discrete(position = "right")+
  geom_tile(aes(fill = category)) +
  xlab(NULL)  +
  theme(legend.title = element_text(size = 50, family = "Times New Roman"),
        legend.text = element_text(size = 50, family = "Times New Roman"),
        axis.title.y= element_blank(),
        axis.title=element_text(size=50,face="bold", family = "Times New Roman"),
        axis.text.x = element_text(angle = 0, hjust = 1, size=50, family = "Times New Roman")) +
  labs(fill = "Phyla")

png(file="/home/jiaqi/Jiaqi_Ni/OliveMicrobiomeCognition_Microbiome/Output/maaslin_taxa.png",width=4500,height=3000, pointsize=50)
grid.newpage()
grid.draw(cbind(ggplotGrob(p1), ggplotGrob(p2), size = "last"))
dev.off()
