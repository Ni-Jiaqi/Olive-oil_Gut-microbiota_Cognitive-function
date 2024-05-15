# Title: Prospective associations of human gut microbiota diversity and composition with changes in cognitive function in 2 years. 
# Author: Jiaqi Ni [ORCID](https://orcid.org/0000-0001-5625-0373)
# Description: This is the analysis script for associations of alpha, beta diversity and per-taxa composition with changes in cognitive function in 2 years 

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
library(sjPlot)
library(Maaslin2)
library(ggalluvial)
library(tidyverse)
library(modelsummary)
library(grid)

###################################### OVERALL COMPOSITION ######################################
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
which(colnames(df)=="Chao1")
alphaDiversity <-colnames(df)[169:ncol(df)]
##with pc1 & pc2
pcadf <- as.data.frame(pca$points)
names(pcadf) <- c("PC1", "PC2")
pcadf$SampleIDcounts <- rownames(pcadf)
##final dataframe
df <- merge(df, pcadf, by = "SampleIDcounts")
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
             "imc_v00", "waist_cir_v00", "exercise_day_v00", "alcoholg_v00", "qalcoholg_v00", "MED12score_v00",
             "Chao1","Simpson" ,"Shannon","InverseSimpson","PC1","PC2")
df[,conVars] <- lapply(df[,conVars], as.numeric)
# Fully adjusted model ---------------------------------------------------
# baseline age, sex, baseline cognitive test score, intervention group, geographical area, educational level, civil status
#BMI, physical activity, smoking status, alcohol consumption, depressive symptomatology, diabetes prevalence, hypertension prevalence, hypercholesterolemia prevalence, 
#adherence to Mediterranean diet (12 scores)
which(colnames(df)=="Chao1")

##GCF
lm_taxa_gcf_f <- data.frame(row.names = colnames(df[,169:ncol(df)]), 
                            taxon = colnames(df[,169:ncol(df)]), 
                            betaestimate = NA, lower95 = NA, upper95 = NA, p = NA, q = NA)
for (i in 169:ncol(df)){
  taxon <- colnames(df)[i]
  fit <- lm(zGCF_2yc ~ scale(df[,i]) +  zGCF_v00 + edad_s1 + sexo_s1 + grupo_int_v00 + edu_level_v00 + civil_status_v00 + area + 
              imc_v00 + exercise_day_v00 + smoke_status_v00 + alcoholg_v00 + qalcoholg_v00 + 
              depression_v00 + diab_prev_s1 + hta_v00 + colest_v00 + MED12score_v00, data = df)
  
  lm_taxa_gcf_f[taxon,"betaestimate"] <- summary(fit)$coefficients[2,1]
  lm_taxa_gcf_f[taxon,"lower95"] <- confint(fit)[2,1]
  lm_taxa_gcf_f[taxon,"upper95"] <- confint(fit)[2,2]
  lm_taxa_gcf_f[taxon,"p"] <- summary(fit)$coefficients[2,4]
}
lm_taxa_gcf_f$Outcome <- factor("Global cognitive function")
##General cognitive function
lm_taxa_gencf_f <- data.frame(row.names = colnames(df[,169:ncol(df)]), 
                              taxon = colnames(df[,169:ncol(df)]), 
                              betaestimate = NA, lower95 = NA, upper95 = NA, p = NA, q = NA)
for (i in 169:ncol(df)){
  taxon <- colnames(df)[i]
  fit <- lm(zgenCF_2yc ~ scale(df[,i]) +  zgenCF_v00 + edad_s1 + sexo_s1 + grupo_int_v00 + edu_level_v00 + civil_status_v00 + area + 
              imc_v00 + exercise_day_v00 + smoke_status_v00 + alcoholg_v00 + qalcoholg_v00 + 
              depression_v00 + diab_prev_s1 + hta_v00 + colest_v00 + MED12score_v00, data = df)
  
  lm_taxa_gencf_f[taxon,"betaestimate"] <- summary(fit)$coefficients[2,1]
  lm_taxa_gencf_f[taxon,"lower95"] <- confint(fit)[2,1]
  lm_taxa_gencf_f[taxon,"upper95"] <- confint(fit)[2,2]
  lm_taxa_gencf_f[taxon,"p"] <- summary(fit)$coefficients[2,4]
}
lm_taxa_gencf_f$Outcome <- factor("General cognitive function")
##Executive function
lm_taxa_ex_f <- data.frame(row.names = colnames(df[,169:ncol(df)]), 
                           taxon = colnames(df[,169:ncol(df)]), 
                           betaestimate = NA, lower95 = NA, upper95 = NA, p = NA, q = NA)
for (i in 169:ncol(df)){
  taxon <- colnames(df)[i]
  fit <- lm(zExF_2yc ~ scale(df[,i]) +  zExF_v00 + edad_s1 + sexo_s1 + grupo_int_v00 + edu_level_v00 + civil_status_v00 + area + 
              imc_v00 + exercise_day_v00 + smoke_status_v00 + alcoholg_v00 + qalcoholg_v00 + 
              depression_v00 + diab_prev_s1 + hta_v00 + colest_v00 + MED12score_v00, data = df)
  
  lm_taxa_ex_f[taxon,"betaestimate"] <- summary(fit)$coefficients[2,1]
  lm_taxa_ex_f[taxon,"lower95"] <- confint(fit)[2,1]
  lm_taxa_ex_f[taxon,"upper95"] <- confint(fit)[2,2]
  lm_taxa_ex_f[taxon,"p"] <- summary(fit)$coefficients[2,4]
}
lm_taxa_ex_f$Outcome <- factor("Executive function")
##Attention
lm_taxa_att_f <- data.frame(row.names = colnames(df[,169:ncol(df)]), 
                            taxon = colnames(df[,169:ncol(df)]), 
                            betaestimate = NA, lower95 = NA, upper95 = NA, p = NA, q = NA)
for (i in 169:ncol(df)){
  taxon <- colnames(df)[i]
  fit <- lm(zattention_2yc ~ scale(df[,i]) +  zattention_v00 + edad_s1 + sexo_s1 + grupo_int_v00 + edu_level_v00 + civil_status_v00 + area + 
              imc_v00 + exercise_day_v00 + smoke_status_v00 + alcoholg_v00 + qalcoholg_v00 + 
              depression_v00 + diab_prev_s1 + hta_v00 + colest_v00 + MED12score_v00, data = df)
  
  lm_taxa_att_f[taxon,"betaestimate"] <- summary(fit)$coefficients[2,1]
  lm_taxa_att_f[taxon,"lower95"] <- confint(fit)[2,1]
  lm_taxa_att_f[taxon,"upper95"] <- confint(fit)[2,2]
  lm_taxa_att_f[taxon,"p"] <- summary(fit)$coefficients[2,4]
}
lm_taxa_att_f$Outcome <- factor("Attention")
##Language
lm_taxa_lang_f <- data.frame(row.names = colnames(df[,169:ncol(df)]), 
                             taxon = colnames(df[,169:ncol(df)]), 
                             betaestimate = NA, lower95 = NA, upper95 = NA, p = NA, q = NA)
for (i in 169:ncol(df)){
  taxon <- colnames(df)[i]
  fit <- lm(zlanguage_2yc ~ scale(df[,i]) +  zlanguage_v00 + edad_s1 + sexo_s1 + grupo_int_v00 + edu_level_v00 + civil_status_v00 + area + 
              imc_v00 + exercise_day_v00 + smoke_status_v00 + alcoholg_v00 + qalcoholg_v00 + 
              depression_v00 + diab_prev_s1 + hta_v00 + colest_v00 + MED12score_v00, data = df)
  
  lm_taxa_lang_f[taxon,"betaestimate"] <- summary(fit)$coefficients[2,1]
  lm_taxa_lang_f[taxon,"lower95"] <- confint(fit)[2,1]
  lm_taxa_lang_f[taxon,"upper95"] <- confint(fit)[2,2]
  lm_taxa_lang_f[taxon,"p"] <- summary(fit)$coefficients[2,4]
}
lm_taxa_lang_f$Outcome <- factor("Language")


overallcog_total <- rbind(lm_taxa_gcf_f, lm_taxa_gencf_f, lm_taxa_ex_f, lm_taxa_att_f, lm_taxa_lang_f)
overallcog_total$q <- p.adjust(overallcog_total$p, method = "BH")
overallcog_total$stars <- cut(overallcog_total$p, breaks=c(-Inf, 0.005, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))

overallcog_total <- overallcog_total %>% 
  mutate(taxon = recode(taxon, "Chao1"="Alpha diversity_Chao1", "InverseSimpson"="Alpha diversity_Inverse Simpson", "Shannon"="Alpha diversity_Shannon", "Simpson"="Alpha diversity_Simpson","PC1"="Beta diversity_PC1","PC2"="Beta diversity_PC2"))

################################################# PER TAXA ###################################################
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
colnames(taxa_table_g) <- sub("^Genus:", "Genus_", colnames(taxa_table_g))
colnames(taxa_table_g)<-gsub("[^[:alnum:][:blank:]+?&/\\-]", " ", colnames(taxa_table_g))
colnames(taxa_table_g) <- sub("^Genus ", "Genus_", colnames(taxa_table_g))
taxa_table_g$SampleIDcounts <- rownames(taxa_table_g)
taxa_sig <- c("Genus_Acidaminococcus", "Genus_Adlercreutzia", "Genus_Akkermansia" ,"Genus_Bacteroides", "Genus_Blautia", "Genus_CAG-56", "Genus_Clostridium sensu stricto 1",   "Genus_Collinsella"                  
              , "Genus_Dorea"   ,                      "Genus_Eubacterium hallii group"     
              , "Genus_Mogibacterium"  ,               "Genus_Phascolarctobacterium"        
              , "Genus_Romboutsia" ,                   "Genus_Streptococcus"                
              , "Genus_   5"     ,                     "Genus_uncultured 8"                 
              , "Genus_Senegalimassilia"  ,            "Genus_Christensenellaceae R-7 group"
              , "Genus_Faecalibacterium", "SampleIDcounts")                


taxa_table_g_sig<- taxa_table_g[, taxa_sig]
df <-  merge(df, taxa_table_g_sig,  by = "SampleIDcounts")
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

which(colnames(df)=="Genus_Acidaminococcus")
# Fully adjusted model ---------------------------------------------------
# baseline age, sex, baseline cognitive test score, intervention group, geographical area, educational level, civil status
#BMI, physical activity, smoking status, alcohol consumption, depressive symptomatology, diabetes prevalence, hypertension prevalence, hypercholesterolemia prevalence, 
#adherence to Mediterranean diet (12 scores)

##GCF
lm_taxa_gcf_f <- data.frame(row.names = colnames(df[,175:ncol(df)]), 
                            taxon = colnames(df[,175:ncol(df)]), 
                            betaestimate = NA, lower95 = NA, upper95 = NA, p = NA, q = NA)
for (i in 175:ncol(df)){
  taxon <- colnames(df)[i]
  fit <- lm(zGCF_2yc ~ scale(df[,i]) +  zGCF_v00 + edad_s1 + sexo_s1 + grupo_int_v00 + edu_level_v00 + civil_status_v00 + area + 
              imc_v00 + exercise_day_v00 + smoke_status_v00 + alcoholg_v00 + qalcoholg_v00 + 
              depression_v00 + diab_prev_s1 + hta_v00 + colest_v00 + MED12score_v00, data = df)
  
  lm_taxa_gcf_f[taxon,"betaestimate"] <- summary(fit)$coefficients[2,1]
  lm_taxa_gcf_f[taxon,"lower95"] <- confint(fit)[2,1]
  lm_taxa_gcf_f[taxon,"upper95"] <- confint(fit)[2,2]
  lm_taxa_gcf_f[taxon,"p"] <- summary(fit)$coefficients[2,4]
}
lm_taxa_gcf_f$Outcome <- factor("Global cognitive function")
##General cognitive function
lm_taxa_gencf_f <- data.frame(row.names = colnames(df[,175:ncol(df)]), 
                              taxon = colnames(df[,175:ncol(df)]), 
                              betaestimate = NA, lower95 = NA, upper95 = NA, p = NA, q = NA)
for (i in 175:ncol(df)){
  taxon <- colnames(df)[i]
  fit <- lm(zgenCF_2yc ~ scale(df[,i]) +  zgenCF_v00 + edad_s1 + sexo_s1 + grupo_int_v00 + edu_level_v00 + civil_status_v00 + area + 
              imc_v00 + exercise_day_v00 + smoke_status_v00 + alcoholg_v00 + qalcoholg_v00 + 
              depression_v00 + diab_prev_s1 + hta_v00 + colest_v00 + MED12score_v00, data = df)
  
  lm_taxa_gencf_f[taxon,"betaestimate"] <- summary(fit)$coefficients[2,1]
  lm_taxa_gencf_f[taxon,"lower95"] <- confint(fit)[2,1]
  lm_taxa_gencf_f[taxon,"upper95"] <- confint(fit)[2,2]
  lm_taxa_gencf_f[taxon,"p"] <- summary(fit)$coefficients[2,4]
}
lm_taxa_gencf_f$Outcome <- factor("General cognitive function")
##Executive function
lm_taxa_ex_f <- data.frame(row.names = colnames(df[,175:ncol(df)]), 
                           taxon = colnames(df[,175:ncol(df)]), 
                           betaestimate = NA, lower95 = NA, upper95 = NA, p = NA, q = NA)
for (i in 175:ncol(df)){
  taxon <- colnames(df)[i]
  fit <- lm(zExF_2yc ~ scale(df[,i]) +  zExF_v00 + edad_s1 + sexo_s1 + grupo_int_v00 + edu_level_v00 + civil_status_v00 + area + 
              imc_v00 + exercise_day_v00 + smoke_status_v00 + alcoholg_v00 + qalcoholg_v00 + 
              depression_v00 + diab_prev_s1 + hta_v00 + colest_v00 + MED12score_v00, data = df)
  
  lm_taxa_ex_f[taxon,"betaestimate"] <- summary(fit)$coefficients[2,1]
  lm_taxa_ex_f[taxon,"lower95"] <- confint(fit)[2,1]
  lm_taxa_ex_f[taxon,"upper95"] <- confint(fit)[2,2]
  lm_taxa_ex_f[taxon,"p"] <- summary(fit)$coefficients[2,4]
}
lm_taxa_ex_f$Outcome <- factor("Executive function")
##Attention
lm_taxa_att_f <- data.frame(row.names = colnames(df[,175:ncol(df)]), 
                            taxon = colnames(df[,175:ncol(df)]), 
                            betaestimate = NA, lower95 = NA, upper95 = NA, p = NA, q = NA)
for (i in 175:ncol(df)){
  taxon <- colnames(df)[i]
  fit <- lm(zattention_2yc ~ scale(df[,i]) +  zattention_v00 + edad_s1 + sexo_s1 + grupo_int_v00 + edu_level_v00 + civil_status_v00 + area + 
              imc_v00 + exercise_day_v00 + smoke_status_v00 + alcoholg_v00 + qalcoholg_v00 + 
              depression_v00 + diab_prev_s1 + hta_v00 + colest_v00 + MED12score_v00, data = df)
  
  lm_taxa_att_f[taxon,"betaestimate"] <- summary(fit)$coefficients[2,1]
  lm_taxa_att_f[taxon,"lower95"] <- confint(fit)[2,1]
  lm_taxa_att_f[taxon,"upper95"] <- confint(fit)[2,2]
  lm_taxa_att_f[taxon,"p"] <- summary(fit)$coefficients[2,4]
}
lm_taxa_att_f$Outcome <- factor("Attention")
##Language
lm_taxa_lang_f <- data.frame(row.names = colnames(df[,175:ncol(df)]), 
                             taxon = colnames(df[,175:ncol(df)]), 
                             betaestimate = NA, lower95 = NA, upper95 = NA, p = NA, q = NA)
for (i in 175:ncol(df)){
  taxon <- colnames(df)[i]
  fit <- lm(zlanguage_2yc ~ scale(df[,i]) +  zlanguage_v00 + edad_s1 + sexo_s1 + grupo_int_v00 + edu_level_v00 + civil_status_v00 + area + 
              imc_v00 + exercise_day_v00 + smoke_status_v00 + alcoholg_v00 + qalcoholg_v00 + 
              depression_v00 + diab_prev_s1 + hta_v00 + colest_v00 + MED12score_v00, data = df)
  
  lm_taxa_lang_f[taxon,"betaestimate"] <- summary(fit)$coefficients[2,1]
  lm_taxa_lang_f[taxon,"lower95"] <- confint(fit)[2,1]
  lm_taxa_lang_f[taxon,"upper95"] <- confint(fit)[2,2]
  lm_taxa_lang_f[taxon,"p"] <- summary(fit)$coefficients[2,4]
}
lm_taxa_lang_f$Outcome <- factor("Language")

taxacog_total <- rbind(lm_taxa_gcf_f, lm_taxa_gencf_f, lm_taxa_ex_f, lm_taxa_att_f, lm_taxa_lang_f)
taxacog_total$q <- p.adjust(taxacog_total$p, method = "BH")
taxacog_total$stars <- cut(taxacog_total$p, breaks=c(-Inf, 0.005, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))




gutcognition <- rbind(overallcog_total, taxacog_total)
gutcognition$padjust <- p.adjust(gutcognition$p, method = "BH")


# read in table of species matched to phyla

phylum <- mia::meltAssay(tse_genus,
                         add_row_data = T,
                         add_col_data = T,
                         assay.type = "clr")


phylum <- phylum[, c("FeatureID", "Phylum")]
phylum_S<- phylum[!duplicated(phylum[, c("FeatureID", "Phylum")]), ]
phylum_S$taxon<-gsub("[^[:alnum:][:blank:]+?&/\\-]", " ", phylum_S$FeatureID)
phylum_S$taxon <- sub("^Genus ", "Genus_", phylum_S$taxon)

gutcognition<-left_join(gutcognition, phylum_S, by="taxon")

gutcognition<-gutcognition[order(gutcognition$Outcome, gutcognition$Phylum, gutcognition$taxon),]


########################################################## PLOT ##################################################
#gutcognition <- gutcognition %>% arrange(gutcognition$Outcome, gutcognition$taxon)
feature_names <- unique(gutcognition$taxon)
# create heatmap for olive oil-taxonomy associations
level_x_order <- factor(gutcognition$Outcome, level = c('Global cognitive function', 'General cognitive function', 'Executive function', 'Attention', 'Language'))
level_y_order <- factor(gutcognition$taxon, levels = feature_names)


p1<-ggplot(gutcognition, aes(level_x_order, level_y_order)) +
  geom_tile(aes(fill = betaestimate), color = "white",
            lwd = 1.5,
            linetype = 1) +
  coord_fixed(ratio = 0.5)+
  scale_fill_gradient2(low = "blue", high = "red", space = "Lab" ) +
  geom_text(aes(label=stars), color="black", size=35, show.legend = TRUE,vjust = 0.7) +
  xlab(NULL) +
  scale_y_discrete(limits=rev)+
  theme(panel.background = element_blank(),
        legend.title = element_text(size = 50, family = "Times New Roman"),
        legend.text = element_text(size = 50, family = "Times New Roman"),
        legend.position = "right",
        plot.title = element_text(size=50, family = "Times New Roman"),
        axis.title=element_text(size=50,face="bold", family = "Times New Roman"),
        axis.title.y= element_blank(),
        axis.text.y=element_text(size=50, face="bold.italic", family = "Times New Roman"),
        axis.text.x = element_text(size=50, angle = 90, vjust = 1, hjust=1,face="bold", family = "Times New Roman")) +
  labs(fill = "β coefficient")+
  guides(fill = guide_colorbar(barwidth = 1.5,
                               barheight = 10))
p1
# create left bar for categorizing species to phyla
gutcognition$category<-as.factor(gutcognition$Phylum)
p2<-ggplot(gutcognition, aes(level_x_order, level_y_order))+
  scale_y_discrete(position = "right")+
  geom_tile(aes(fill = rev(category))) +
  xlab(NULL)  +
  theme(legend.title = element_text(size = 50, family = "Times New Roman"),
        legend.text = element_text(size = 50, family = "Times New Roman"),
        axis.title.y= element_blank(),
        axis.title=element_text(size=50,face="bold", family = "Times New Roman"),
        axis.text.x = element_text(angle = 0, hjust = 1, size=10, family = "Times New Roman")) +
  labs(fill = "Phyla")

png(file="/home/jiaqi/Jiaqi_Ni/OliveMicrobiomeCognition_Microbiome/Output/microbiome_cognition_nofdr.png",width=4500,height=2500, pointsize=50)
grid.newpage()
grid.draw(cbind(ggplotGrob(p1), ggplotGrob(p2), size = "last"))
dev.off()

########################################################## PLOT p-adjusted ##################################################
overallcog_total$stars <- cut(overallcog_total$q, breaks=c(-Inf, 0.05, 0.1, 0.25, Inf), label=c("***", "**", "*", ""))

taxacog_total$stars <- cut(taxacog_total$q, breaks=c(-Inf, 0.05, 0.1, 0.25, Inf), label=c("***", "**", "*", ""))

gutcognition <- rbind(overallcog_total, taxacog_total)
gutcognition$padjust <- p.adjust(gutcognition$p, method = "BH")


# read in table of species matched to phyla

phylum <- mia::meltAssay(tse_genus,
                         add_row_data = T,
                         add_col_data = T,
                         assay.type = "clr")


phylum <- phylum[, c("FeatureID", "Phylum")]
phylum_S<- phylum[!duplicated(phylum[, c("FeatureID", "Phylum")]), ]
phylum_S$taxon<-gsub("[^[:alnum:][:blank:]+?&/\\-]", " ", phylum_S$FeatureID)
phylum_S$taxon <- sub("^Genus ", "Genus_", phylum_S$taxon)

gutcognition<-left_join(gutcognition, phylum_S, by="taxon")

gutcognition<-gutcognition[order(gutcognition$Outcome, gutcognition$Phylum, gutcognition$taxon),]


#gutcognition <- gutcognition %>% arrange(gutcognition$Outcome, gutcognition$taxon)
feature_names <- unique(gutcognition$taxon)
# create heatmap for olive oil-taxonomy associations
level_x_order <- factor(gutcognition$Outcome, level = c('Global cognitive function', 'General cognitive function', 'Executive function', 'Attention', 'Language'))
level_y_order <- factor(gutcognition$taxon, levels = feature_names)


p1<-ggplot(gutcognition, aes(level_x_order, level_y_order)) +
  geom_tile(aes(fill = betaestimate), color = "white",
            lwd = 1.5,
            linetype = 1) +
  coord_fixed(ratio = 0.5)+
  scale_fill_gradient2(low = "blue", high = "red", space = "Lab" ) +
  # scale_fill_gradient2(low = "#075AFF",
  #                      mid = "#FFFFCC",
  #                      high = "#FF0000") +
  geom_text(aes(label=stars), color="black", size=35, show.legend = TRUE,vjust = 0.7) +
  xlab(NULL) +
  scale_y_discrete(limits=rev)+
  theme(panel.background = element_blank(),
        legend.title = element_text(size = 50, family = "Times New Roman"),
        legend.text = element_text(size = 50, family = "Times New Roman"),
        legend.position = "right",
        plot.title = element_text(size=50, family = "Times New Roman"),
        axis.title=element_text(size=50,face="bold", family = "Times New Roman"),
        axis.title.y= element_blank(),
        axis.text.y=element_text(size=50, face="bold.italic", family = "Times New Roman"),
        axis.text.x = element_text(size=50, angle = 90, vjust = 1, hjust=1,face="bold", family = "Times New Roman")) +
  labs(fill = "β coefficient")+
  guides(fill = guide_colorbar(barwidth = 1.5,
                               barheight = 10))
p1
# create left bar for categorizing species to phyla
gutcognition$category<-as.factor(gutcognition$Phylum)
p2<-ggplot(gutcognition, aes(level_x_order, level_y_order))+
  scale_y_discrete(position = "right")+
  geom_tile(aes(fill = rev(category))) +
  xlab(NULL)  +
  theme(legend.title = element_text(size = 50, family = "Times New Roman"),
        legend.text = element_text(size = 50, family = "Times New Roman"),
        axis.title.y= element_blank(),
        axis.title=element_text(size=50,face="bold", family = "Times New Roman"),
        axis.text.x = element_text(angle = 0, hjust = 1, size=10, family = "Times New Roman")) +
  labs(fill = "Phyla")

png(file="/home/jiaqi/Jiaqi_Ni/OliveMicrobiomeCognition_Microbiome/Output/microbiome_cognition_adjust.png",width=4500,height=2500, pointsize=50)
grid.newpage()
grid.draw(cbind(ggplotGrob(p1), ggplotGrob(p2), size = "last"))
dev.off()


