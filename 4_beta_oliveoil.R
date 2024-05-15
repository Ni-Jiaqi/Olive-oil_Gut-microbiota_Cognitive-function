# Title: Cross-sectional associations of olive oil consumption with human gut microbiota diversity and composition
# Author: Jiaqi Ni [ORCID](https://orcid.org/0000-0001-5625-0373)
# Description: This is the analysis script for beta diversity analyses

########################################################### STARTUP ################################################################
# Set up working environment  ------------------------------------------------------------------
rm(list=ls())
getwd()
setwd("/home/jiaqi/Jiaqi_Ni/OliveMicrobiomeCognition_Microbiome")
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
ps <- readRDS("./DATA/jiaqi_Phyloseq_110324.rds") #baseline microbiota N=656

# convert phyloseq to TSE
tse <- makeTreeSummarizedExperimentFromPhyloseq(ps) 
tse

# * Aitchison distance calculation---------------------------------------------------
##CLR transform
tse <- transformAssay(x = tse, assay.type = "counts", method = "clr", 
                      pseudocount = 1, name = "clr")

##Pick the clr assay table
taxa <- assays(tse)$clr
# Calculate Aitchison distance (euclidean distance of CLR-transformed species abundances)
beta <- vegan::vegdist(t(taxa), method = "euclidean")
beta_mtx <- as.matrix(beta)

# Calculate the cumulative variance explained by the PCs
pca <- prcomp(beta_mtx, center = TRUE, scale. = TRUE)
var_explained = format(round((pca$sdev^2 / sum(pca$sdev^2))*100,1),nsmall=1)
cum_var_explained <- cumsum(var_explained) 
# Take a PCA and choose the first 2 principal components
pca <- cmdscale(beta, k = 2, eig = TRUE)

# * extract data for dataframe ---------------------------------------------------
pca <- as.data.frame(pca$points)
##get metadata as dataframe
df <- as.data.frame(colData(tse))
# Define and fix continuous and categorical variables
catVars <- c("sexo_s1",  "nodo", "area", "grupo_int_v00", "edu_level_v00", "civil_status_v00",
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

# PCA ---------------------------------------------------
# *  total olive oil  ---------------------------------------------------
# Plot based on tertiles of total olive oil consumption 
pca$ter_ajust_olivatot_v00 <- as.factor(df$ter_ajust_olivatot_v00)
names(pca) <- c("PC1", "PC2", "total_olive_oil_consumption_tertile")
group <- pca[,"total_olive_oil_consumption_tertile", drop = FALSE]
group <- group[order(group[,1]),,drop = FALSE]
pca <- pca[match(rownames(group), rownames(pca)),]

pca.total.plot <- ggplot2::ggplot(data = pca, aes(PC1, PC2, colour = pca$total_olive_oil_consumption_tertile)) + 
  geom_point(size = 10, alpha = 1) +
  scale_color_manual(values=c("#fbea96", "#f8a804", "#bc6404")) +
  stat_ellipse(aes(group = pca$total_olive_oil_consumption_tertile), type = "t", level = 0.95, lwd = 2) +
  xlab(paste0("PC1 (", var_explained[1], "%)")) +
  ylab(paste0("PC2 (", var_explained[2], "%)")) +
  theme_minimal() +
  labs(colour = paste0("TOO tertiles:")) +
  theme(text = element_text(size=25, family = "Times New Roman"),
        legend.text = element_text(size = 20),
        legend.position = c(0.92, 0.1), 
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.background = element_rect(fill="white", size=.5, linetype="solid"))
print.total.plot <- pca.total.plot + ggpubr::border()
print.total.plot

pca.total.plot <- ggMarginal(
  print.total.plot,
  type = c("boxplot"),
  margins= "both",
  size =10,
  groupColour = F,
  groupFill = T
) 
pca.total.plot
# *  extra virgin olive oil  ---------------------------------------------------
pca <- prcomp(beta_mtx, center = TRUE, scale. = TRUE)
var_explained = format(round((pca$sdev^2 / sum(pca$sdev^2))*100,1),nsmall=1)
cum_var_explained <- cumsum(var_explained) 
# Take a PCA and choose the first 2 principal components
pca <- cmdscale(beta, k = 2, eig = TRUE)
pca <- as.data.frame(pca$points)

# Plot based on tertiles of extra virgin olive oil consumption 
pca$ter_ajust_ac_olivavir_v00 <- as.factor(df$ter_ajust_ac_olivavir_v00)
names(pca) <- c("PC1", "PC2", "extra_virgin_olive_oil_consumption_tertile")
group <- pca[,"extra_virgin_olive_oil_consumption_tertile", drop = FALSE]
group <- group[order(group[,1]),,drop = FALSE]
pca <- pca[match(rownames(group), rownames(pca)),]

pca.ex.plot <- ggplot2::ggplot(data = pca, aes(PC1, PC2, colour = pca$extra_virgin_olive_oil_consumption_tertile)) + 
  geom_point(size = 10, alpha = 1) +
  scale_color_manual(values=c("#deeeba", "#99CC66", "#336600")) +
  stat_ellipse(aes(group = pca$extra_virgin_olive_oil_consumption_tertile), type = "t", level = 0.95, lwd = 2) +
  xlab(paste0("PC1 (", var_explained[1], "%)")) +
  ylab(paste0("PC2 (", var_explained[2], "%)")) +
  theme_minimal() +
  labs(colour = paste0("VOO tertiles:")) +
  theme(text = element_text(size=25, family = "Times New Roman"),
        legend.text = element_text(size = 20),
        legend.position = c(0.92, 0.1), 
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.background = element_rect(fill="white", size=.5, linetype="solid"))
print.ex.plot <- pca.ex.plot + ggpubr::border()
print.ex.plot

pca.ex.plot <- ggMarginal(
  print.ex.plot,
  type = c("boxplot"),
  margins= "both",
  size =10,
  groupColour = F,
  groupFill = T
) 
pca.ex.plot

# * common olive oil  ---------------------------------------------------
pca <- prcomp(beta_mtx, center = TRUE, scale. = TRUE)
var_explained = format(round((pca$sdev^2 / sum(pca$sdev^2))*100,1),nsmall=1)
cum_var_explained <- cumsum(var_explained) 
# Take a PCA and choose the first 2 principal components
pca <- cmdscale(beta, k = 2, eig = TRUE)
pca <- as.data.frame(pca$points)

# Plot based on tertiles of common olive oil consumption 
pca$ter_ajust_re_oliva_v00 <- as.factor(df$ter_ajust_re_oliva_v00)
names(pca) <- c("PC1", "PC2", "refined_olive_oil_consumption_tertile")
group <- pca[,"refined_olive_oil_consumption_tertile", drop = FALSE]
group <- group[order(group[,1]),,drop = FALSE]
pca <- pca[match(rownames(group), rownames(pca)),]

pca.re.plot <- ggplot2::ggplot(data = pca, aes(PC1, PC2, colour = pca$refined_olive_oil_consumption_tertile)) + 
  geom_point(size = 10, alpha = 1) +
  scale_color_manual(values=c("bisque2", "tomato2", "darkred")) +
  stat_ellipse(aes(group = pca$refined_olive_oil_consumption_tertile), type = "t", level = 0.95, lwd = 2) +
  xlab(paste0("PC1 (", var_explained[1], "%)")) +
  ylab(paste0("PC2 (", var_explained[2], "%)")) +
  theme_minimal() +
  labs(colour = paste0("COO tertiles:")) +
  theme(text = element_text(size=25, family = "Times New Roman"),
        legend.position = c(0.92, 0.1), 
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.background = element_rect(fill="white", size=.5, linetype="solid"))
print.re.plot <- pca.re.plot + ggpubr::border()
print.re.plot

pca.re.plot <- ggMarginal(
  print.re.plot,
  type = c("boxplot"),
  margins= "both",
  size =10,
  groupColour = F,
  groupFill = T
) 
pca.re.plot

combined_beta <- ggarrange(pca.total.plot, pca.ex.plot, pca.re.plot, ncol=1, nrow=3)
combined_beta



## COMBINE 3 FIGURES
png(
  file = "/home/jiaqi/Jiaqi_Ni/OliveMicrobiomeCognition_Microbiome/Output/beta.png",
  width = 6000,
  height = 2000,
  pointsize = 50
)
ggarrange(pca.total.plot, pca.ex.plot, pca.re.plot, ncol=3, nrow=1)
dev.off()


###################################################################  PERMANOVA  #################################################
pca <- prcomp(beta_mtx, center = TRUE, scale. = TRUE)
# Take a PCA and choose the first 2 principal components
pca <- cmdscale(beta, k = 2, eig = TRUE)
pcoap <- as.data.frame(pca$points)
names(pcoap) <- c("PC1", "PC2")
# combine this PCoA data frame with metadata 
pcoap$SampleIDcounts<-rownames(pcoap)
pcoap<-pcoap[order(pcoap$SampleIDcounts),]
pcoap <-merge(pcoap, df, by="SampleIDcounts")
rownames(pcoap)<-pcoap$SampleIDcounts
# *  PERMANOVA plot (variance explained by exposures, outcomes,covariates)---------------------------------------------------
# Set seed for reproducibility
set.seed(1996)
x <- adonis2(beta ~ ajust_olivatot_v00, df, method = "euclidean", permutations = 999, by = "terms")
olicon<-x[!is.na(x$F), ]
x <- adonis2(beta ~ ter_ajust_olivatot_v00, df, method = "euclidean", permutations = 999, by = "terms")
olicat<-x[!is.na(x$F), ]
x <- adonis2(beta ~ ajust_ac_olivavir_v00, df, method = "euclidean", permutations = 999, by = "terms")
excon<-x[!is.na(x$F), ]
x <- adonis2(beta ~ ter_ajust_ac_olivavir_v00, df, method = "euclidean", permutations = 999, by = "terms")
excat<-x[!is.na(x$F), ]
x <- adonis2(beta ~ ajust_re_oliva_v00, df, method = "euclidean", permutations = 999, by = "terms")
recon<-x[!is.na(x$F), ]
x <- adonis2(beta ~ ter_ajust_re_oliva_v00, df, method = "euclidean", permutations = 999, by = "terms")
recat<-x[!is.na(x$F), ]

x <- adonis2(beta ~ zGCF_v00, df, method = "euclidean", permutations = 999, by = "terms")
gcf<-x[!is.na(x$F), ]
x <- adonis2(beta ~ zgenCF_v00, df, method = "euclidean", permutations = 999, by = "terms")
gencf<-x[!is.na(x$F), ]
x <- adonis2(beta ~ zExF_v00, df, method = "euclidean", permutations = 999, by = "terms")
exf<-x[!is.na(x$F), ]
x <- adonis2(beta ~ zattention_v00, df, method = "euclidean", permutations = 999, by = "terms")
att<-x[!is.na(x$F), ]
x <- adonis2(beta ~ zlanguage_v00, df, method = "euclidean", permutations = 999, by = "terms")
lan<-x[!is.na(x$F), ]
x <- adonis2(beta ~ zMMSE_v00, df, method = "euclidean", permutations = 999, by = "terms")
mmse<-x[!is.na(x$F), ]
x <- adonis2(beta ~ zGCF_2yc, df, method = "euclidean", permutations = 999, by = "terms")
gcf2y<-x[!is.na(x$F), ]
x <- adonis2(beta ~ zgenCF_2yc, df, method = "euclidean", permutations = 999, by = "terms")
gencf2y<-x[!is.na(x$F), ]
x <- adonis2(beta ~ zExF_2yc, df, method = "euclidean", permutations = 999, by = "terms")
exf2y<-x[!is.na(x$F), ]
x <- adonis2(beta ~ zattention_2yc, df, method = "euclidean", permutations = 999, by = "terms")
att2y<-x[!is.na(x$F), ]
x <- adonis2(beta ~ zlanguage_2yc, df, method = "euclidean", permutations = 999, by = "terms")
lan2y<-x[!is.na(x$F), ]
x <- adonis2(beta ~ zMMSE_2yc, df, method = "euclidean", permutations = 999, by = "terms")
mmse2y<-x[!is.na(x$F), ]

x <- adonis2(beta ~ edad_s1, df, method = "euclidean", permutations = 999, by = "terms")
age<-x[!is.na(x$F), ]
x <- adonis2(beta ~ sexo_s1, df, method = "euclidean", permutations = 999, by = "terms")
sex<-x[!is.na(x$F), ]
x <- adonis2(beta ~ edu_level_v00, df, method = "euclidean", permutations = 999, by = "terms")
edu<-x[!is.na(x$F), ]
x <- adonis2(beta ~ civil_status_v00, df, method = "euclidean", permutations = 999, by = "terms")
civ<-x[!is.na(x$F), ]
x <- adonis2(beta ~ area, df, method = "euclidean", permutations = 999, by = "terms")
center<-x[!is.na(x$F), ]
x <- adonis2(beta ~ imc_v00, df, method = "euclidean", permutations = 999, by = "terms")
bmi<-x[!is.na(x$F), ]
x <- adonis2(beta ~ exercise_day_v00, df, method = "euclidean", permutations = 999, by = "terms")
phy<-x[!is.na(x$F), ]
x <- adonis2(beta ~ smoke_status_v00, df, method = "euclidean", permutations = 999, by = "terms")
smoke<-x[!is.na(x$F), ]
x <- adonis2(beta ~ alcoholg_v00, df, method = "euclidean", permutations = 999, by = "terms")
alcohol<-x[!is.na(x$F), ]
x <- adonis2(beta ~ qalcoholg_v00, df, method = "euclidean", permutations = 999, by = "terms")
qalcohol<-x[!is.na(x$F), ]
x <- adonis2(beta ~ depression_v00, df, method = "euclidean", permutations = 999, by = "terms")
depre<-x[!is.na(x$F), ]
x <- adonis2(beta ~ diab_prev_s1, df, method = "euclidean", permutations = 999, by = "terms")
diab<-x[!is.na(x$F), ]
x <- adonis2(beta ~ hta_v00, df, method = "euclidean", permutations = 999, by = "terms")
hta<-x[!is.na(x$F), ]
x <- adonis2(beta ~ colest_v00, df, method = "euclidean", permutations = 999, by = "terms")
choles<-x[!is.na(x$F), ]
x <- adonis2(beta ~ antidepressant_v00, df, method = "euclidean", permutations = 999, by = "terms")
antidep<-x[!is.na(x$F), ]
x <- adonis2(beta ~ hypolipidemic_v00, df, method = "euclidean", permutations = 999, by = "terms")
hypolip<-x[!is.na(x$F), ]
x <- adonis2(beta ~ antidiabetic_v00, df, method = "euclidean", permutations = 999, by = "terms")
antidiab<-x[!is.na(x$F), ]
x <- adonis2(beta ~ antihypertensive_v00, df, method = "euclidean", permutations = 999, by = "terms")
antihyper<-x[!is.na(x$F), ]
x <- adonis2(beta ~ MED12score_v00, df, method = "euclidean", permutations = 999, by = "terms")
med12<-x[!is.na(x$F), ]
# merge PERMANOVA results
all_taxonomy.1 <-rbind(olicon,excon, recon, 
                       gcf, gencf, exf, att, lan, 
                       age, sex, edu, civ, center, bmi, phy, smoke, alcohol, depre, diab, hta, choles, med12)

# label for variables
component.1<- c('Total olive oil',
                'Virgin olive oil',
                'Common olive oil',
                'Global cognitive function',
                'General cognitive function',
                'Executive function',
                'Attention',
                'Language',
                'Age',
                'Sex',
                'Education level',
                'Civil status',
                'Geographical area',
                'BMI',
                'Physical activity',
                'Smoking status',
                'Alcohol consumption',
                'Depressive symptomatology',
                'Diabetes prevalence',
                'Hypertension prevalence',
                'Hypercholesterolemia prevalence',
                'Mediterranean diet adherence')

all_taxonomy.1<-cbind(all_taxonomy.1, component.1)
all_taxonomy_oil_1<-all_taxonomy.1[1:3,]
all_taxonomy_oil_1$cat<-"Types of olive oil"
all_taxonomy_cogn_1<-all_taxonomy.1[4:8,]
all_taxonomy_cogn_1$cat<-"Baseline cognitive performance"
all_taxonomy_cov_1<-all_taxonomy.1[9:22,]
all_taxonomy_cov_1$cat<-"Covariables"
all_taxonomy_tax_1<-rbind(all_taxonomy_oil_1, all_taxonomy_cogn_1, all_taxonomy_cov_1)
all_taxonomy_tax_1$cat<-factor(all_taxonomy_tax_1$cat, levels = c("Types of olive oil","Baseline cognitive performance",  "Covariables"))
# FDR adjustment for p-values
all_taxonomy_tax_1$fdr <-
  p.adjust(
    all_taxonomy_tax_1$`Pr(>F)`,
    method = "BH",
    n = length(all_taxonomy_tax_1$`Pr(>F)`)
  )

# create figure to show PERMANOVA results
all_taxonomy_tax_1$stars <- cut(all_taxonomy_tax_1$fdr, breaks=c(-Inf, 0.005, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))
level_y_order <- factor(all_taxonomy_tax_1$component, level = rev(all_taxonomy_tax_1$component))
color_code<-c("Types of olive oil"="#E69F00",
              "Baseline cognitive performance"="#56B4E9",
              "Covariables"="#999999")

permanova.plot1<- ggplot(all_taxonomy_tax_1, aes(level_y_order, R2, fill=cat)) +
  geom_bar(stat="identity", width = 0.6)+
  theme_classic()+
  scale_y_continuous(labels=scales::percent, limits = c(0, 0.005))+
  geom_text(aes(label=stars), color="black", size=5) +
  scale_fill_manual(values = color_code)+ coord_flip()+
  theme(axis.line = element_line(colour = "black", 
                                 linewidth = 0.3, linetype = "solid"),
        legend.position = "right",
        legend.text = element_text(size = 15,color="black"),
        legend.title = element_blank(),
        plot.title = element_blank(),
        axis.title.x=element_blank(),
        axis.title.y= element_blank(),
        axis.text.y=element_text(size=15,color="black"),
        axis.text.x = element_text(size=15,color="black",),
        text = element_text(family = "Times New Roman"))
# save plot
png(file="/home/jiaqi/Jiaqi_Ni/OliveMicrobiomeCognition_Microbiome/Output/permanova_variables.png",width=1200,height=800, pointsize=20)
permanova.plot1 + theme(plot.margin = margin(1, 1, 1, 1, "cm")) 
dev.off()

# *  PERMANOVA tests (variance across tertiles)---------------------------------------------------

# * *  total olive oil ------------------------------------------------------------- 
permanova.1.1 <- vegan::adonis2(beta ~  ter_ajust_olivatot_v00 +  edad_s1 + sexo_s1 , df, method = "euclidean", permutations = 999, by = "margin")
permanova.1.1
permanova.1.2 <- vegan::adonis2(beta ~  ter_ajust_olivatot_v00 +  edad_s1 + sexo_s1 + area + edu_level_v00 + civil_status_v00, df, method = "euclidean", permutations = 999, by = "margin")
permanova.1.2
permanova.1.3 <- vegan::adonis2(beta ~  ter_ajust_olivatot_v00 +  edad_s1 + sexo_s1 + area + edu_level_v00 + civil_status_v00 + 
                                  imc_v00 + exercise_day_v00 + smoke_status_v00 + alcoholg_v00 + qalcoholg_v00 + depression_v00 + diab_prev_s1 + hta_v00 + colest_v00 +
                                  MED12score_v00, df, method = "euclidean", permutations = 999, by = "margin")
permanova.1.3
# * *  extra virgin olive oil ------------------------------------------------------------- 
permanova.2.1 <- vegan::adonis2(beta ~  ter_ajust_ac_olivavir_v00 +  edad_s1 + sexo_s1 , df, method = "euclidean", permutations = 999, by = "margin")
permanova.2.1
permanova.2.2 <- vegan::adonis2(beta ~  ter_ajust_ac_olivavir_v00 +  edad_s1 + sexo_s1 + area + edu_level_v00 + civil_status_v00 , df, method = "euclidean", permutations = 999, by = "margin")
permanova.2.2
permanova.2.3 <- vegan::adonis2(beta ~  ter_ajust_ac_olivavir_v00 +  edad_s1 + sexo_s1 + area + edu_level_v00 + civil_status_v00
                                + imc_v00 + exercise_day_v00 + smoke_status_v00 + alcoholg_v00 + qalcoholg_v00 + depression_v00 + diab_prev_s1 + hta_v00 + colest_v00 +
                                  MED12score_v00, df, method = "euclidean", permutations = 999, by = "margin")
permanova.2.3
# * *  common olive oil ------------------------------------------------------------- 
permanova.3.1 <- vegan::adonis2(beta ~  ter_ajust_re_oliva_v00 +  edad_s1 + sexo_s1 , df, method = "euclidean", permutations = 999, by = "margin")
permanova.3.1
permanova.3.2 <- vegan::adonis2(beta ~  ter_ajust_re_oliva_v00 +  edad_s1 + sexo_s1 + area + edu_level_v00 + civil_status_v00 , df, method = "euclidean", permutations = 999, by = "margin")
permanova.3.2
permanova.3.3 <- vegan::adonis2(beta ~  ter_ajust_re_oliva_v00 +  edad_s1 + sexo_s1 + area + edu_level_v00 + civil_status_v00
                                + imc_v00 + exercise_day_v00 + smoke_status_v00 + alcoholg_v00 + qalcoholg_v00 + depression_v00 + diab_prev_s1 + hta_v00 + colest_v00 +
                                  MED12score_v00, df, method = "euclidean", permutations = 999, by = "margin")
permanova.3.3

# P-adjustment
basic.adj.beta <- data.frame()
basic.adj.beta[1,"analysis"] <- "beta too basic"
basic.adj.beta[1,"pvalue"] <- permanova.1.1$"Pr(>F)"[1]
basic.adj.beta[2,"analysis"] <- "beta too socio"
basic.adj.beta[2,"pvalue"] <- permanova.1.2$"Pr(>F)"[1]
basic.adj.beta[3,"analysis"] <- "beta too full"
basic.adj.beta[3,"pvalue"] <- permanova.1.3$"Pr(>F)"[1]
basic.adj.beta[4,"analysis"] <- "beta voo basic"
basic.adj.beta[4,"pvalue"] <- permanova.2.1$"Pr(>F)"[1]
basic.adj.beta[5,"analysis"] <- "beta voo socio"
basic.adj.beta[5,"pvalue"] <- permanova.2.2$"Pr(>F)"[1]
basic.adj.beta[6,"analysis"] <- "beta voo full"
basic.adj.beta[6,"pvalue"] <- permanova.2.3$"Pr(>F)"[1]
basic.adj.beta[7,"analysis"] <- "beta coo basic"
basic.adj.beta[7,"pvalue"] <- permanova.3.1$"Pr(>F)"[1]
basic.adj.beta[8,"analysis"] <- "beta coo socio"
basic.adj.beta[8,"pvalue"] <- permanova.3.2$"Pr(>F)"[1]
basic.adj.beta[9,"analysis"] <- "beta coo full"
basic.adj.beta[9,"pvalue"] <- permanova.3.3$"Pr(>F)"[1]
basic.adj.beta$qvalue <- p.adjust(basic.adj.beta$pvalue, method = "BH")
################################################################## REGRESSION ANALYSIS (PC1&PC2) ###################################################
# Basic model ------------------------------------------------------------- 

# baseline age, sex
## Set variables for Loop
betaindex <-c("PC1", "PC2")
### * total olive oil  (cat.) -------------------------------------------------------------
baselinebasicresults <- list()
for (i in betaindex){
  formula <- as.formula(paste(i, "~ ter_ajust_olivatot_v00 +  edad_s1 + sexo_s1"))
  baselinebasicmodels <- lm ( formula, data = pcoap)
  baselinebasicresults[[i]] <- summary(baselinebasicmodels)
  
}

for (i in betaindex){
  print(paste("Regression results for", i))
  print(baselinebasicresults[[i]])
  
}

modelsummary(baselinebasicresults, 
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "p.value")
### total olive oil  (p-trend)
baselinebasicresults <- list()
for (i in betaindex){
  formula <- as.formula(paste(i, "~ ptrend_ter_ajust_olivatot_v00 +  edad_s1 + sexo_s1"))
  baselinebasicmodels <- lm ( formula, data = pcoap)
  baselinebasicresults[[i]] <- summary(baselinebasicmodels)
  
}

for (i in betaindex){
  print(paste("Regression results for", i))
  print(baselinebasicresults[[i]])
  
}

modelsummary(baselinebasicresults, 
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "p.value")
### * total olive oil  (cont.) -------------------------------------------------------------
baselinebasicresults <- list()
for (i in betaindex){
  formula <- as.formula(paste(i, "~ ajust_tspolivatot_v00 +  edad_s1 + sexo_s1"))
  baselinebasicmodels <- lm ( formula, data = pcoap)
  baselinebasicresults[[i]] <- summary(baselinebasicmodels)
  
}

for (i in betaindex){
  print(paste("Regression results for", i))
  print(baselinebasicresults[[i]])
  
}

modelsummary(baselinebasicresults, 
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "p.value")
### * extra virgin olive oil  (cat.) -------------------------------------------------------------
baselinebasicresults <- list()
for (i in betaindex){
  formula <- as.formula(paste(i, "~ ter_ajust_ac_olivavir_v00 +  edad_s1 + sexo_s1"))
  baselinebasicmodels <- lm ( formula, data = pcoap)
  baselinebasicresults[[i]] <- summary(baselinebasicmodels)
  
}

for (i in betaindex){
  print(paste("Regression results for", i))
  print(baselinebasicresults[[i]])
  
}

modelsummary(baselinebasicresults, 
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "p.value")
### extra virgin olive oil  (p-trend)
baselinebasicresults <- list()
for (i in betaindex){
  formula <- as.formula(paste(i, "~ ptrend_ter_ajust_ac_olivavir_v00 +  edad_s1 + sexo_s1"))
  baselinebasicmodels <- lm ( formula, data = pcoap)
  baselinebasicresults[[i]] <- summary(baselinebasicmodels)
  
}

for (i in betaindex){
  print(paste("Regression results for", i))
  print(baselinebasicresults[[i]])
  
}

modelsummary(baselinebasicresults, 
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "p.value")
### * extra virgin olive oil  (cont.) -------------------------------------------------------------
baselinebasicresults <- list()
for (i in betaindex){
  formula <- as.formula(paste(i, "~ ajust_tspac_olivavir_v00 +  edad_s1 + sexo_s1"))
  baselinebasicmodels <- lm ( formula, data = pcoap)
  baselinebasicresults[[i]] <- summary(baselinebasicmodels)
  
}

for (i in betaindex){
  print(paste("Regression results for", i))
  print(baselinebasicresults[[i]])
  
}

modelsummary(baselinebasicresults, 
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "p.value")
### * common olive oil  (cat.) -------------------------------------------------------------
baselinebasicresults <- list()
for (i in betaindex){
  formula <- as.formula(paste(i, "~ ter_ajust_re_oliva_v00 +  edad_s1 + sexo_s1"))
  baselinebasicmodels <- lm ( formula, data = pcoap)
  baselinebasicresults[[i]] <- summary(baselinebasicmodels)
  
}

for (i in betaindex){
  print(paste("Regression results for", i))
  print(baselinebasicresults[[i]])
  
}

modelsummary(baselinebasicresults, 
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "p.value")
### refined olive oil  (p-trend)
baselinebasicresults <- list()
for (i in betaindex){
  formula <- as.formula(paste(i, "~ ptrend_ter_ajust_re_oliva_v00 +  edad_s1 + sexo_s1"))
  baselinebasicmodels <- lm ( formula, data = pcoap)
  baselinebasicresults[[i]] <- summary(baselinebasicmodels)
  
}

for (i in betaindex){
  print(paste("Regression results for", i))
  print(baselinebasicresults[[i]])
  
}

modelsummary(baselinebasicresults, 
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "p.value")
### * common olive oil  (cont.) -------------------------------------------------------------
baselinebasicresults <- list()
for (i in betaindex){
  formula <- as.formula(paste(i, "~ ajust_tspre_oliva_v00 +  edad_s1 + sexo_s1"))
  baselinebasicmodels <- lm ( formula, data = pcoap)
  baselinebasicresults[[i]] <- summary(baselinebasicmodels)
  
}

for (i in betaindex){
  print(paste("Regression results for", i))
  print(baselinebasicresults[[i]])
  
}

modelsummary(baselinebasicresults, 
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "p.value")
# Socio-demographically adjusted model ------------------------------------------------------------- 
# baseline age, sex,  geographical area, educational level, civil status
### * total olive oil  (cat.) -------------------------------------------------------------
baselinebasicresults <- list()
for (i in betaindex){
  formula <- as.formula(paste(i, "~ ter_ajust_olivatot_v00 +  edad_s1 + sexo_s1 + area + edu_level_v00 + civil_status_v00"))
  baselinebasicmodels <- lm ( formula, data = pcoap)
  baselinebasicresults[[i]] <- summary(baselinebasicmodels)
  
}

for (i in betaindex){
  print(paste("Regression results for", i))
  print(baselinebasicresults[[i]])
  
}

modelsummary(baselinebasicresults, 
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "p.value")
### total olive oil  (p-trend)
baselinebasicresults <- list()
for (i in betaindex){
  formula <- as.formula(paste(i, "~ ptrend_ter_ajust_olivatot_v00 +  edad_s1 + sexo_s1 + area + edu_level_v00 + civil_status_v00"))
  baselinebasicmodels <- lm ( formula, data = pcoap)
  baselinebasicresults[[i]] <- summary(baselinebasicmodels)
  
}

for (i in betaindex){
  print(paste("Regression results for", i))
  print(baselinebasicresults[[i]])
  
}

modelsummary(baselinebasicresults, 
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "p.value")

### * total olive oil  (cont.) -------------------------------------------------------------
baselinebasicresults <- list()
for (i in betaindex){
  formula <- as.formula(paste(i, "~ ajust_tspolivatot_v00 +  edad_s1 + sexo_s1 + area + edu_level_v00 + civil_status_v00"))
  baselinebasicmodels <- lm ( formula, data = pcoap)
  baselinebasicresults[[i]] <- summary(baselinebasicmodels)
  
}

for (i in betaindex){
  print(paste("Regression results for", i))
  print(baselinebasicresults[[i]])
  
}

modelsummary(baselinebasicresults, 
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "p.value")
### * extra virgin olive oil  (cat.) -------------------------------------------------------------
baselinebasicresults <- list()
for (i in betaindex){
  formula <- as.formula(paste(i, "~ ter_ajust_ac_olivavir_v00 +  edad_s1 + sexo_s1 + area + edu_level_v00 + civil_status_v00"))
  baselinebasicmodels <- lm ( formula, data = pcoap)
  baselinebasicresults[[i]] <- summary(baselinebasicmodels)
  
}

for (i in betaindex){
  print(paste("Regression results for", i))
  print(baselinebasicresults[[i]])
  
}

modelsummary(baselinebasicresults, 
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "p.value")
### extra virgin olive oil  (p-trend)
baselinebasicresults <- list()
for (i in betaindex){
  formula <- as.formula(paste(i, "~ ptrend_ter_ajust_ac_olivavir_v00 +  edad_s1 + sexo_s1 + area + edu_level_v00 + civil_status_v00"))
  baselinebasicmodels <- lm ( formula, data = pcoap)
  baselinebasicresults[[i]] <- summary(baselinebasicmodels)
  
}

for (i in betaindex){
  print(paste("Regression results for", i))
  print(baselinebasicresults[[i]])
  
}

modelsummary(baselinebasicresults, 
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "p.value")
### * extra virgin olive oil  (cont.) -------------------------------------------------------------
baselinebasicresults <- list()
for (i in betaindex){
  formula <- as.formula(paste(i, "~ ajust_tspac_olivavir_v00 +  edad_s1 + sexo_s1 + area + edu_level_v00 + civil_status_v00"))
  baselinebasicmodels <- lm ( formula, data = pcoap)
  baselinebasicresults[[i]] <- summary(baselinebasicmodels)
  
}

for (i in betaindex){
  print(paste("Regression results for", i))
  print(baselinebasicresults[[i]])
  
}

modelsummary(baselinebasicresults, 
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "p.value")
### * common olive oil  (cat.) -------------------------------------------------------------
baselinebasicresults <- list()
for (i in betaindex){
  formula <- as.formula(paste(i, "~ ter_ajust_re_oliva_v00 +  edad_s1 + sexo_s1+ area + edu_level_v00 + civil_status_v00"))
  baselinebasicmodels <- lm ( formula, data = pcoap)
  baselinebasicresults[[i]] <- summary(baselinebasicmodels)
  
}

for (i in betaindex){
  print(paste("Regression results for", i))
  print(baselinebasicresults[[i]])
  
}

modelsummary(baselinebasicresults, 
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "p.value")
### refined olive oil  (p-trend)
baselinebasicresults <- list()
for (i in betaindex){
  formula <- as.formula(paste(i, "~ ptrend_ter_ajust_re_oliva_v00 +  edad_s1 + sexo_s1 + area + edu_level_v00 + civil_status_v00"))
  baselinebasicmodels <- lm ( formula, data = pcoap)
  baselinebasicresults[[i]] <- summary(baselinebasicmodels)
  
}

for (i in betaindex){
  print(paste("Regression results for", i))
  print(baselinebasicresults[[i]])
  
}

modelsummary(baselinebasicresults, 
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "p.value")
### * common olive oil  (cont.) -------------------------------------------------------------
baselinebasicresults <- list()
for (i in betaindex){
  formula <- as.formula(paste(i, "~ ajust_tspre_oliva_v00 +  edad_s1 + sexo_s1 + area + edu_level_v00 + civil_status_v00"))
  baselinebasicmodels <- lm ( formula, data = pcoap)
  baselinebasicresults[[i]] <- summary(baselinebasicmodels)
  
}

for (i in betaindex){
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
for (i in betaindex){
  formula <- as.formula(paste(i, "~ ter_ajust_olivatot_v00 +  edad_s1 + sexo_s1 + area + edu_level_v00 + civil_status_v00 + 
                              imc_v00 + exercise_day_v00 + smoke_status_v00 + alcoholg_v00 + qalcoholg_v00 + 
                              depression_v00 + diab_prev_s1 + hta_v00 + colest_v00 + MED12score_v00"))
  baselinebasicmodels <- lm (formula, data = pcoap)
  baselinebasicresults[[i]] <- summary(baselinebasicmodels)
  
}

for (i in betaindex){
  print(paste("Regression results for", i))
  print(baselinebasicresults[[i]])
  
}

modelsummary(baselinebasicresults, 
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "p.value")

### total olive oil  (p-trend)
baselinebasicresults <- list()
for (i in betaindex){
  formula <- as.formula(paste(i, "~ ptrend_ter_ajust_olivatot_v00 +  edad_s1 + sexo_s1 + area + edu_level_v00 + civil_status_v00 + 
                              imc_v00 + exercise_day_v00 + smoke_status_v00 + alcoholg_v00 + qalcoholg_v00 + 
                              depression_v00 + diab_prev_s1 + hta_v00 + colest_v00 + MED12score_v00"))
  baselinebasicmodels <- lm ( formula, data = pcoap)
  baselinebasicresults[[i]] <- summary(baselinebasicmodels)
  
}

for (i in betaindex){
  print(paste("Regression results for", i))
  print(baselinebasicresults[[i]])
  
}

modelsummary(baselinebasicresults, 
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "p.value")

### * total olive oil  (cont.) -------------------------------------------------------------
baselinebasicresults <- list()
for (i in betaindex){
  formula <- as.formula(paste(i, "~ ajust_tspolivatot_v00 +  edad_s1 + sexo_s1 + area + edu_level_v00 + civil_status_v00 + 
                              imc_v00 + exercise_day_v00 + smoke_status_v00 + alcoholg_v00 + qalcoholg_v00 + 
                              depression_v00 + diab_prev_s1 + hta_v00 + colest_v00 + MED12score_v00"))
  baselinebasicmodels <- lm ( formula, data = pcoap)
  baselinebasicresults[[i]] <- summary(baselinebasicmodels)
  
}

for (i in betaindex){
  print(paste("Regression results for", i))
  print(baselinebasicresults[[i]])
  
}

modelsummary(baselinebasicresults, 
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "p.value")
### * extra virgin olive oil  (cat.) -------------------------------------------------------------
baselinebasicresults <- list()
for (i in betaindex){
  formula <- as.formula(paste(i, "~ ter_ajust_ac_olivavir_v00 +  edad_s1 + sexo_s1 + area + edu_level_v00 + civil_status_v00 + 
                              imc_v00 + exercise_day_v00 + smoke_status_v00 + alcoholg_v00 + qalcoholg_v00 + 
                              depression_v00 + diab_prev_s1 + hta_v00 + colest_v00 + MED12score_v00"))
  baselinebasicmodels <- lm ( formula, data = pcoap)
  baselinebasicresults[[i]] <- summary(baselinebasicmodels)
  
}

for (i in betaindex){
  print(paste("Regression results for", i))
  print(baselinebasicresults[[i]])
  
}

modelsummary(baselinebasicresults, 
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "p.value")
### extra virgin olive oil  (p-trend)
baselinebasicresults <- list()
for (i in betaindex){
  formula <- as.formula(paste(i, "~ ptrend_ter_ajust_ac_olivavir_v00 + edad_s1 + sexo_s1 + area + edu_level_v00 + civil_status_v00 + 
                              imc_v00 + exercise_day_v00 + smoke_status_v00 + alcoholg_v00 + qalcoholg_v00 + 
                              depression_v00 + diab_prev_s1 + hta_v00 + colest_v00 + MED12score_v00"))
  baselinebasicmodels <- lm ( formula, data = pcoap)
  baselinebasicresults[[i]] <- summary(baselinebasicmodels)
  
}

for (i in betaindex){
  print(paste("Regression results for", i))
  print(baselinebasicresults[[i]])
  
}

modelsummary(baselinebasicresults, 
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "p.value")
### * extra virgin olive oil  (cont.) -------------------------------------------------------------
baselinebasicresults <- list()
for (i in betaindex){
  formula <- as.formula(paste(i, "~ ajust_tspac_olivavir_v00 +  edad_s1 + sexo_s1 + area + edu_level_v00 + civil_status_v00 + 
                              imc_v00 + exercise_day_v00 + smoke_status_v00 + alcoholg_v00 + qalcoholg_v00 + 
                              depression_v00 + diab_prev_s1 + hta_v00 + colest_v00 + MED12score_v00"))
  baselinebasicmodels <- lm ( formula, data = pcoap)
  baselinebasicresults[[i]] <- summary(baselinebasicmodels)
  
}

for (i in betaindex){
  print(paste("Regression results for", i))
  print(baselinebasicresults[[i]])
  
}

modelsummary(baselinebasicresults, 
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "p.value")

### * common olive oil  (cat.) -------------------------------------------------------------
baselinebasicresults <- list()
for (i in betaindex){
  formula <- as.formula(paste(i, "~ ter_ajust_re_oliva_v00 +  edad_s1 + sexo_s1 + area + edu_level_v00 + civil_status_v00 + 
                              imc_v00 + exercise_day_v00 + smoke_status_v00 + alcoholg_v00 + qalcoholg_v00 + 
                              depression_v00 + diab_prev_s1 + hta_v00 + colest_v00 + MED12score_v00"))
  baselinebasicmodels <- lm ( formula, data = pcoap)
  baselinebasicresults[[i]] <- summary(baselinebasicmodels)
  
}

for (i in betaindex){
  print(paste("Regression results for", i))
  print(baselinebasicresults[[i]])
  
}

modelsummary(baselinebasicresults, 
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "p.value")
### refined olive oil  (p-trend)
baselinebasicresults <- list()
for (i in betaindex){
  formula <- as.formula(paste(i, "~ ptrend_ter_ajust_re_oliva_v00 +  edad_s1 + sexo_s1 + area + edu_level_v00 + civil_status_v00 + 
                              imc_v00 + exercise_day_v00 + smoke_status_v00 + alcoholg_v00 + qalcoholg_v00 + 
                              depression_v00 + diab_prev_s1 + hta_v00 + colest_v00 + MED12score_v00"))
  baselinebasicmodels <- lm ( formula, data = pcoap)
  baselinebasicresults[[i]] <- summary(baselinebasicmodels)
  
}

for (i in betaindex){
  print(paste("Regression results for", i))
  print(baselinebasicresults[[i]])
  
}

modelsummary(baselinebasicresults, 
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "p.value")
### * common olive oil  (cont.) -------------------------------------------------------------
baselinebasicresults <- list()
for (i in betaindex){
  formula <- as.formula(paste(i, "~ ajust_tspre_oliva_v00 +  edad_s1 + sexo_s1 + area + edu_level_v00 + civil_status_v00 + 
                              imc_v00 + exercise_day_v00 + smoke_status_v00 + alcoholg_v00 + qalcoholg_v00 + 
                              depression_v00 + diab_prev_s1 + hta_v00 + colest_v00 + MED12score_v00"))
  baselinebasicmodels <- lm ( formula, data = pcoap)
  baselinebasicresults[[i]] <- summary(baselinebasicmodels)
  
}

for (i in betaindex){
  print(paste("Regression results for", i))
  print(baselinebasicresults[[i]])
  
}

modelsummary(baselinebasicresults, 
             estimate = "{estimate} [{conf.low}, {conf.high}]",
             statistic = "p.value")

