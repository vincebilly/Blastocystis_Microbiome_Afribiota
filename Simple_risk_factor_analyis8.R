####################################################################################################################################
####################################################################################################################################
######                                                                                                                        ######
######                              RISK FACTOR ANALYSIS WITH qPCR  (AFRIBIOTA BLASTOCYSTIS)                                  ######
######                                                                                                                        ######
####################################################################################################################################
####################################################################################################################################

rm(list=ls())

##################################################################
###                      IMPORT LIBRARY                        ###
##################################################################

library("tidyverse")
library("reshape2")    # Load the reshape2 package for converting between long and wide format data
library("stringr")     # Load the stringr package for improved filtering of text
library("ape")         # Load the ape package for reading and modifying phylogenetic trees
library("phyloseq")    # Load the phyloseq package for microbial community analysis
library("data.table")  # Load the data.table package for better metadata manipulation
library("viridis")     # Load the viridis package for colour palettes for continuous data
library("qualpalr")    # Load the qualpalr package for colour palettes for qualitative data
library("ggplot2")     # load the ggplot2 package for visualization of data
library("vegan")       # 
library("gplots")      #
library("ggpubr")      #
library("ggforce")
library("plyr")
library("lattice")
library("Rmisc")
library("funrar")
library("permute")
library("readr")
library("forcats")
library("RColorBrewer")
library("ade4")
library("cluster")
library("clue")
library("dada2")
library("permute")
library("broom")
library("viridisLite")
library("seqinr")
library("ShortRead")
library("Biostrings")
library("VennDiagram")
library("DESeq2")
library("ggtree")
library("tidyr")
library("ggimage")
library("lubridate")    #Add vertical line in time serie plot
library("onewaytests")
library("cowplot")
library("qpdf")
library("ggVennDiagram")
library("ggpmisc")
library("vcd")
library("ggmosaic")
library("gtools")
library("forestplot")
library("ggalt")
library("dplyr")
library("ggrepel")
library("aplot")       #
library("gridExtra")   #
library("ggeffects")   #
library("rstatix")
library("emmeans")
library("samplesizeCMH")
require("FactoMineR")
library("corrplot")
library("DataExplorer")
library("pheatmap")
library("ggbiplot")
library("pls")
library("factoextra")
library("MASS")
library("ca")
library("homals")
library("cluster")
library("graphics")
require("grDevices")
library("CCA")
library("labdsv")
library("caTools")
library("randomForest")
library("rstatix")
library("ggvenn")
library("chemodiv")
library("rasterdiv")
library("indicspecies")

##################################################################

#################################################################################
#################################################################################
###                     IMPORT and CLEAN/FILTER METADATA                      ###
#################################################################################
#################################################################################

##################################################################
###             IMPORT METADATA and FILTER/CLEAN               ###
##################################################################

#load Pascale's filtered metadata
metadata1 <- read_delim("/Users/vincebilly/Desktop/vince/phd/assignment/afribiota_16s/Pascale_2023/forRachelleriskfactors_filtered.csv", delim = ",", trim_ws = TRUE)
dim(metadata1)

#load in calprotectin data, combine with other metadata
calprotectin_data <- read_delim("/Users/vincebilly/Desktop/vince/phd/assignment/afribiota_16s/Pascale_2023/CALPROTECTINE_data.txt", delim = "\t", trim_ws = TRUE)
dim(calprotectin_data)

#Checking if samples are similar between the 2 dataset
intersect(metadata1$id,calprotectin_data$id)
setdiff(metadata1$id,calprotectin_data$id)
setdiff(calprotectin_data$id,metadata1$id)

#Merge Metadata with Calprotectin data
metadata <- merge(metadata1, calprotectin_data, by = "id", all=FALSE);dim(metadata)

##################################################################

##################################################################
###                   Exclusion criteria                       ###
##################################################################

#Convert integer into numerical values
metadata$haz_cont <- as.numeric(metadata$haz_cont)
metadata$whz_cont <- as.numeric(metadata$whz_cont)
metadata$CALPROTECTINEggdePS <- as.numeric(metadata$CALPROTECTINEggdePS)
metadata$AATmggdePS <- as.numeric(metadata$AATmggdePS)

#Remove samples with haz > 2 ans <-2
#metadata <- subset(metadata, haz_cont >-2);dim(metadata)

##################################################################

#################################################################################
#################################################################################
###                        IMPORT and CLEAN/FILTER 16S                        ###
#################################################################################
#################################################################################

##################################################################
###                           New 16s                          ###
##################################################################

##Import ASV table
##################################################################

asv_16s1_spl.tbl <- as.data.frame(fread("/Users/vincebilly/Desktop/vince/phd/assignment/afribiota_16s/Pascale_2023/ASV_table16S_Pascale_2023.csv", 
                              sep=",", header = T, stringsAsFactors = FALSE));str(asv_16s1_spl.tbl)
row.names(asv_16s1_spl.tbl) <- asv_16s1_spl.tbl$V1
asv_16s1_spl.tbl$V1 <- NULL
colnames(asv_16s1_spl.tbl)

asv_16s1_spl.tbl <- t(asv_16s1_spl.tbl)

#Remove prefix "_SXXX"
row.names(asv_16s1_spl.tbl) <- gsub("_.*","",row.names(asv_16s1_spl.tbl))
row.names(asv_16s1_spl.tbl) <- gsub("ER.","SE.",row.names(asv_16s1_spl.tbl))


#asv_16s1_spl.tbl$SampleID <- row.names(asv_16s1_spl.tbl)

spl_code <- as.data.frame(fread("/Users/vincebilly/Desktop/vince/phd/assignment/afribiota_16s/Pascale_2023/sentforsequencingMI2.csv", 
                                        sep=",", header = T, stringsAsFactors = FALSE));str(spl_code)

intersect(row.names(asv_16s1_spl.tbl), spl_code$samplename_metagenomics)


row.names(asv_16s1_spl.tbl) <- gsub("ER","",row.names(asv_16s1_spl.tbl))
row.names(asv_16s1_spl.tbl) <- gsub("\\.1|\\.2|\\.3","",row.names(asv_16s1_spl.tbl))

#Select only gastric aspiration ("AG.")
asv_16s1_spl_AG.tbl <- asv_16s1_spl.tbl[grep("AG.",row.names(asv_16s1_spl.tbl)),];dim(asv_16s1_spl_AG.tbl)
#196 samples
#Select only duodenal aspiration ("AD.")
asv_16s1_spl_AD.tbl <- asv_16s1_spl.tbl[grep("AD.",row.names(asv_16s1_spl.tbl)),];dim(asv_16s1_spl_AD.tbl)
# 104 samples

#Select only samples ("S0.")
asv_16s1_spl_S0.tbl <- asv_16s1_spl.tbl[grep("S0.",row.names(asv_16s1_spl.tbl)),];dim(asv_16s1_spl_S0.tbl)
# 535 samples
intersect(row.names(asv_16s1_spl_S0.tbl), spl_code$samplename_metagenomics)
asv_16s1_spl_S0_code <- merge(spl_code,asv_16s1_spl_S0.tbl, by.x="samplename_metagenomics",by.y="row.names")
asv_16s1_spl_S0_code[,c(1,3,4)] <- NULL
asv_16s1_spl_S0_code <- asv_16s1_spl_S0_code[grep("SE.",asv_16s1_spl_S0_code$Sample_Name),];dim(asv_16s1_spl_S0_code)
row.names(asv_16s1_spl_S0_code) <- asv_16s1_spl_S0_code$Sample_Name
asv_16s1_spl_S0_code$Sample_Name <- NULL


#Select only fecal samples ("SE.")
asv_16s1_spl_SE.tbl <- asv_16s1_spl.tbl[grep("SE.",row.names(asv_16s1_spl.tbl)),];dim(asv_16s1_spl_SE.tbl)
#334 samples

#Rbind S0 and SE
dim(asv_16s1_spl_S0_code);dim(asv_16s1_spl_SE.tbl)
asv_16s1_spl_SE_all <- rbind(asv_16s1_spl_S0_code,asv_16s1_spl_SE.tbl)
row.names(asv_16s1_spl_SE_all) <- gsub("SE.","",row.names(asv_16s1_spl_SE_all))
row.names(asv_16s1_spl_SE_all) <- gsub("SE-","",row.names(asv_16s1_spl_SE_all))

#Sum of samples
334+196+104+535


##################################################################

##Import taxonomy table
##################################################################

asv_16s1_tax.tbl <- as.data.frame(fread("/Users/vincebilly/Desktop/vince/phd/assignment/afribiota_16s/Pascale_2023/tax_table16S_Pascale_2023.csv", 
                                        sep=",", header = T, stringsAsFactors = FALSE));str(asv_16s1_tax.tbl)
row.names(asv_16s1_tax.tbl) <- asv_16s1_tax.tbl$V1
asv_16s1_tax.tbl$V1 <- NULL
colnames(asv_16s1_tax.tbl)

#Remove column with ACTG sequence
asv_16s1_tax.tbl$V10<- NULL

#Extract Taxonomy table
asv_tax.tbl <- asv_16s1_tax.tbl[,c("Rank1","Rank2","Rank3","Rank4","Rank5","Rank6","Rank7")]
asv_tax.tbl[is.na(asv_tax.tbl)] <- "Un"

#create a function to paste "Un" to the name of the previous rank
mypaste = function(s) {
  s = gsub("_Un","",s)
  return(paste0(s,"_Un"))
}

#apply the function to all ranks
for (i in 2:8){
  
  if (i == 8)
  {
    colName = "Rank7"
    colNamePrevious = "Rank6"
  } else
  {
    colName = paste0("Rank", i)
    colNamePrevious = paste0("Rank", i-1)
  }
  
  w = which(asv_tax.tbl[[colName]] == "Un")
  asv_tax.tbl[[colName]][w] = mypaste(asv_tax.tbl[[colNamePrevious]][w])
}
asv_tax.mtx <- as.matrix(asv_tax.tbl);dim(asv_tax.mtx)
##################################################################

##Create Phyloseq object and Filter
##################################################################

# Create Phyloseq Object
otu = otu_table(asv_16s1_spl_SE_all, taxa_are_rows = FALSE)
tax = tax_table(asv_tax.mtx)
project_data <- phyloseq(otu,tax);project_data

#Remove ASVs that are in groups that you wish to exclude (e.g. Chloroplast and Mitochondria), making sure that NAs within those Ranks are not removed in the process
project_data <- project_data %>%
  subset_taxa(Rank4 != "Chloroplast") %>% 
  subset_taxa(Rank5 != "Mitochondria") %>% 
  subset_taxa(Rank6 != "Blastocystis")  %>%
  subset_taxa(Rank2 != "Archaeplastida") %>% #plants and algae
  subset_taxa(Rank4 != "Porifera") %>% #sponges
  subset_taxa(Rank4 != "Arthropoda") %>% #arthropods
  subset_taxa(Rank5 != "Arthropoda") %>% #arthropods
  subset_taxa(Rank6 != "Arthropoda") %>% #arthropods
  subset_taxa(Rank4 != "Mollusca") #mussels;
project_data

#Convert ASVs with counts of 2 or less to 0
otu_conversion <- as.data.frame(otu_table(project_data))
otu_table(project_data)[otu_conversion <= 2] <- 0;project_data
#Prune ASV that have less than X total reads
project_data <- prune_taxa(taxa_sums(project_data) >= 250, project_data);project_data

#rarefy
sort(sample_sums(project_data))
project_data.rarefied_initial <- rarefy_even_depth(project_data, sample.size = 5000);project_data.rarefied_initial

##################################################################

#Calculate Alpha diversity
##################################################################

otu <- otu_table(project_data.rarefied_initial)
Shannon16s_new <- diversity(otu, index = "shannon", MARGIN = 1, base = exp(1))
Simpson16s_new <- diversity(otu, index = "simpson", MARGIN = 1, base = exp(1))
otu[otu > 0] <- 1
Number_ASVs16s_new <- rowSums(otu)
meta_shannon16s_new <- cbind(Shannon16s_new, Number_ASVs16s_new,Simpson16s_new);dim(meta_shannon16s_new)

##################################################################

##################################################################
###                      IMPORT qPCR DATA                      ###
##################################################################

# Import qPCR table
##################################################################

qPCR_table1A <- read.csv2("/Users/vincebilly/Desktop/vince/phd/qpcr_blastocystis/run_and_results/2022-09-07/WellResult2.csv", sep=",")
qPCR_table1B <- read.csv2("/Users/vincebilly/Desktop/vince/phd/qpcr_blastocystis/run_and_results/2022-09-08/WellResult2.csv", sep=",")
qPCR_table2C <- read.csv2("/Users/vincebilly/Desktop/vince/phd/qpcr_blastocystis/run_and_results/2022-09-09/WellResult2.csv", sep=",")
qPCR_table2D <- read.csv2("/Users/vincebilly/Desktop/vince/phd/qpcr_blastocystis/run_and_results/2022-09-14/WellResult2.csv", sep=",")
qPCR_table3E <- read.csv2("/Users/vincebilly/Desktop/vince/phd/qpcr_blastocystis/run_and_results/2022-09-15/WellResult2.csv", sep=",")
qPCR_table3F <- read.csv2("/Users/vincebilly/Desktop/vince/phd/qpcr_blastocystis/run_and_results/2022-09-16/WellResult2.csv", sep=",")
qPCR_table7M <- read.csv2("/Users/vincebilly/Desktop/vince/phd/qpcr_blastocystis/run_and_results/2022-11-04-MEG3/WellResult2.csv", sep=",")
qPCR_table7N <- read.csv2("/Users/vincebilly/Desktop/vince/phd/qpcr_blastocystis/run_and_results/2022-11-04-MEG4/WellResult2.csv", sep=",")
qPCR_table9Q <- read.csv2("/Users/vincebilly/Desktop/vince/phd/qpcr_blastocystis/run_and_results/2022-09-27/WellResult2.csv", sep=",")
qPCR_table9R <- read.csv2("/Users/vincebilly/Desktop/vince/phd/qpcr_blastocystis/run_and_results/2022-09-29/WellResult2.csv", sep=",")
qPCR_table10S <- read.csv2("/Users/vincebilly/Desktop/vince/phd/qpcr_blastocystis/run_and_results/2022-10-21-MEG4/WellResult2.csv", sep=",")
qPCR_table10T <- read.csv2("/Users/vincebilly/Desktop/vince/phd/qpcr_blastocystis/run_and_results/2022-10-21-MEG3/WellResult2.csv", sep=",")
qPCR_table11U <- read.csv2("/Users/vincebilly/Desktop/vince/phd/qpcr_blastocystis/run_and_results/2022-10-04/WellResult2.csv", sep=",")
qPCR_table11V <- read.csv2("/Users/vincebilly/Desktop/vince/phd/qpcr_blastocystis/run_and_results/2022-10-08/WellResult2.csv", sep=",")
qPCR_table12W <- read.csv2("/Users/vincebilly/Desktop/vince/phd/qpcr_blastocystis/run_and_results/2022-10-11/WellResult2.csv", sep=",")
qPCR_table12X <- read.csv2("/Users/vincebilly/Desktop/vince/phd/qpcr_blastocystis/run_and_results/2022-10-13/WellResult2.csv", sep=",")
qPCR_table8O <- read.csv2("/Users/vincebilly/Desktop/vince/phd/qpcr_blastocystis/run_and_results/2022-11-28-MEG3/WellResult2.csv", sep=",")
qPCR_table8P <- read.csv2("/Users/vincebilly/Desktop/vince/phd/qpcr_blastocystis/run_and_results/2022-11-28-MEG4/WellResult2.csv", sep=",")
#qPCR_table2C_duplicate <- read.csv2("/Users/vincebilly/Desktop/vince/phd/qpcr_blastocystis/run_and_results/2022-12-09-MEG3/WellResult2.csv", sep=",")
#qPCR_table2D_duplicate <- read.csv2("/Users/vincebilly/Desktop/vince/phd/qpcr_blastocystis/run_and_results/2022-12-09-MEG4/WellResult3.csv", sep=",")
qPCR_table4G <- read.csv2("/Users/vincebilly/Desktop/vince/phd/qpcr_blastocystis/run_and_results/2022-12-12-MEG4/WellResult2.csv", sep=",")
qPCR_table4H <- read.csv2("/Users/vincebilly/Desktop/vince/phd/qpcr_blastocystis/run_and_results/2022-12-12-MEG3/WellResult2.csv", sep=",")


#Combine table
qPCR_table <- rbind(qPCR_table1A,qPCR_table1B,
                    qPCR_table2C,qPCR_table2D,
                    qPCR_table3E,qPCR_table3F,
                    qPCR_table7M,qPCR_table7N,
                    qPCR_table9Q,qPCR_table9R,
                    qPCR_table10S,qPCR_table10T,
                    qPCR_table11U,qPCR_table11V,
                    qPCR_table12W,qPCR_table12X,
                    qPCR_table8O,qPCR_table8P,
                    #qPCR_table2C_duplicate,qPCR_table2D_duplicate,
                    qPCR_table4G,qPCR_table4H)

qPCR_table$Sample_eppendorf_tube <- qPCR_table$Sample


#Remove suffix and Prefix from sample ID (SE, ADN, ADN3, ADN4, SE3, SE4, -, 4, SE9)
qPCR_table$Sample <- gsub("SE3-ADN4-","",qPCR_table$Sample)
qPCR_table$Sample <- gsub("SE1-ADN1","",qPCR_table$Sample)
qPCR_table$Sample <- gsub("SE3-ADN","",qPCR_table$Sample)
qPCR_table$Sample <- gsub("ADN-SE-","",qPCR_table$Sample)
qPCR_table$Sample <- gsub("-ADN2","",qPCR_table$Sample)
qPCR_table$Sample <- gsub("-ADN3","",qPCR_table$Sample)
qPCR_table$Sample <- gsub("-ADN4","",qPCR_table$Sample)
qPCR_table$Sample <- gsub("-ADN","",qPCR_table$Sample)
qPCR_table$Sample <- gsub("ADN-","",qPCR_table$Sample)
qPCR_table$Sample <- gsub("ADN","",qPCR_table$Sample)
qPCR_table$Sample <- gsub("DNA-","",qPCR_table$Sample)
qPCR_table$Sample <- gsub("-SED","",qPCR_table$Sample)
qPCR_table$Sample <- gsub("-SE34","",qPCR_table$Sample)
qPCR_table$Sample <- gsub("-SE3","",qPCR_table$Sample)
qPCR_table$Sample <- gsub("-SE4","",qPCR_table$Sample)
qPCR_table$Sample <- gsub("-SE9","",qPCR_table$Sample)
qPCR_table$Sample <- gsub("SE3-","",qPCR_table$Sample)
qPCR_table$Sample <- gsub("SE1-","",qPCR_table$Sample)
qPCR_table$Sample <- gsub("SE-","",qPCR_table$Sample)
qPCR_table$Sample <- gsub("-SE","",qPCR_table$Sample)
qPCR_table$Sample <- gsub("-4","",qPCR_table$Sample)
qPCR_table$Sample <- gsub("4-","",qPCR_table$Sample)
qPCR_table$Sample <- gsub("-","",qPCR_table$Sample)
qPCR_table$Sample <- gsub("ADB4","",qPCR_table$Sample)
qPCR_table$Sample <- gsub("AD","",qPCR_table$Sample)

#Check if samples match with MISEQ data
match_spl_metadata_qPCR <- setdiff(qPCR_table$Sample,metadata$id);match_spl_metadata_qPCR
sample_to_be_relabbled <- match_spl_metadata_qPCR[grep("B", match_spl_metadata_qPCR)]
sample_to_be_relabbled <- sample_to_be_relabbled[-grep("CPB",sample_to_be_relabbled)]

#
qPCR_table$Sample[which(qPCR_table$Sample %in% sample_to_be_relabbled)] <- gsub("B","",qPCR_table$Sample[which(qPCR_table$Sample %in% sample_to_be_relabbled)])
match_spl_metadata_qPCR <- setdiff(qPCR_table$Sample,metadata$id);match_spl_metadata_qPCR

#Correct HM178 to HMET178 (HM178 was written that way on eppendorf tube)
qPCR_table$Sample[which(qPCR_table$Sample %in% "HM178")] <- gsub("HM178","HMET178",qPCR_table$Sample[which(qPCR_table$Sample %in% "HM178")])
#Correct HME008 to HMET008 (HME008 was written that way on eppendorf tube)
qPCR_table$Sample[which(qPCR_table$Sample %in% "HME008")] <- gsub("HME008","HMET008",qPCR_table$Sample[which(qPCR_table$Sample %in% "HME008")])
qPCR_table$Sample[which(qPCR_table$Sample %in% "CPB60")] <- gsub("CPB60","CPB060",qPCR_table$Sample[which(qPCR_table$Sample %in% "CPB60")])
qPCR_table$Sample[which(qPCR_table$Sample %in% "CPB73")] <- gsub("CPB73","CPB073",qPCR_table$Sample[which(qPCR_table$Sample %in% "CPB73")])

match_spl_metadata_qPCR <- setdiff(qPCR_table$Sample,metadata$id);match_spl_metadata_qPCR

#
length(intersect(qPCR_table$Sample,metadata$id))
length(setdiff(metadata$id,qPCR_table$Sample))
length(setdiff(qPCR_table$Sample,metadata$id)) - 8 # Remove 8 samples which are the standard (1.10^4 - 1.10^11)


#Convert Undetermined Cq Value to numeric and 
qPCR_table$Cq <- gsub("Undetermined", "0", qPCR_table$Cq)
qPCR_table <- subset(qPCR_table, Cq!= "Issue")
qPCR_table$Cq <- as.numeric(qPCR_table$Cq);dim(qPCR_table)

#
qPCR_table$Standard <- gsub("Fourh","Fourth",qPCR_table$Standard)

##################################################################

## qPCR filtering steps
##################################################################

#Remove IPC and Standard
no_ipc <- subset(qPCR_table, Target == "BLASTO");dim(no_ipc)
no_ipc_no_std <- no_ipc[-grep("Standard", no_ipc$Sample),];dim(no_ipc_no_std)

# Caluclate mean of Cq value iwithin duplicate
plate1 <- no_ipc_no_std[,c("Sample","Cq")]
plate2_sup0 <- subset(plate1, Cq > 0)
plate2_sup0 <- aggregate(plate2_sup0, by = list(plate2_sup0$Sample ), FUN = mean)
plate2_inf0 <- subset(plate1, Cq == 0)
plate2_inf0 <- aggregate(plate2_inf0, by = list(plate2_inf0$Sample),FUN = mean) #Aggregate at Rank8
to_keep <-  setdiff(plate2_inf0$Group.1, plate2_sup0$Group.1)
plate2_inf0 <- subset(plate2_inf0, Group.1 %in% to_keep)
plate3 <- rbind(plate2_sup0,plate2_inf0)
plate3$Sample <- NULL
colnames(plate3) <- c("Sample","Cq");dim(plate3)

##################################################################

##################################################################
###        IMPORT/CLEAN/FILTER  MISEQ MICROBIOME INSIGHT       ###
##################################################################

#Import 18s otu table
##################################################################
asv_mi_tbl <- read.table("/Users/vincebilly/Desktop/vince/phd/assignment/afribiota/sequence_table.18s_MI.all_years_all_sample_types.corrected_sampleIDs.correct_taxIDs.txt", header=T);dim(asv_mi_tbl)
asv_mi_tbl_agg <- asv_mi_tbl[grep("SE",asv_mi_tbl$row_names),]
asv_mi_tbl_agg$row_names <- gsub("SE-","",asv_mi_tbl_agg$row_names)
row.names(asv_mi_tbl_agg) <- asv_mi_tbl_agg[,1]
asv_mi_tbl_agg[,1] <- NULL
asv_mi_tbl_agg <- as.data.frame(t(asv_mi_tbl_agg))
head(asv_mi_tbl_agg)

#Import taxonomy
asv_tax.tbl <- read.table("/Users/vincebilly/Desktop/vince/phd/assignment/afribiota/tax_midal.txt", 
                          sep="\t", header = T, quote = "", stringsAsFactors = FALSE);dim(asv_tax.tbl)
asv_tax.tbl <- asv_tax.tbl[grep("MI",asv_tax.tbl$MIrow_names),]
asv_tax.tbl$MIrow_names <- gsub("MI","",asv_tax.tbl$MIrow_names)
row.names(asv_tax.tbl) <- asv_tax.tbl$MIrow_names
asv_tax.tbl$MIrow_names <- NULL
#Convert "NA" value into "Un" (for Unknown)
asv_tax.tbl[is.na(asv_tax.tbl)] <- "Un"

#create a function to paste "Un" to the name of the previous rank
mypaste = function(s) {
  s = gsub("_Un","",s)
  return(paste0(s,"_Un"))
}
#apply the function to all ranks
for (i in 2:9){
  if (i == 9)
  {
    colName = "Accession"
    colNamePrevious = "Rank8"
  } else
  {
    colName = paste0("Rank", i)
    colNamePrevious = paste0("Rank", i-1)
  }
  w = which(asv_tax.tbl[[colName]] == "Un")
  asv_tax.tbl[[colName]][w] = mypaste(asv_tax.tbl[[colNamePrevious]][w])
}

asv_tax_MI.tbl <- asv_tax.tbl

#merge OTU TABLE and TAXONOMY TABLE
asv_tax_otu_MI <- merge(asv_tax_MI.tbl, asv_mi_tbl_agg, by="row.names", all.y=FALSE)
row.names(asv_tax_otu_MI) <- asv_tax_otu_MI[,1]
asv_tax_otu_MI[,1] <-  NULL

##################################################################

##Filtering step
##################################################################

# Create Phyloseq Object
otu_euk = otu_table(asv_mi_tbl_agg, taxa_are_rows = TRUE)
tax_euk = tax_table(as.matrix(asv_tax.tbl))
project_data_euk <- phyloseq(otu_euk,tax_euk)

#rarefy
sort(sample_sums(project_data_euk))
project_data_euk <- rarefy_even_depth(project_data_euk, sample.size = 1000);project_data_euk


#Remove ASVs that are in groups that you wish to exclude (e.g. Chloroplast and Mitochondria), making sure that NAs within those Ranks are not removed in the process
project_data_euk <- project_data_euk %>%
  subset_taxa(Rank5 != "Mitochondria")

#Convert ASVs with counts of 2 or less to 0
otu_conversion <- as.data.frame(otu_table(project_data_euk))
otu_table(project_data_euk)[otu_conversion <= 2] <- 0;project_data_euk

#Keep only intestinal euakaryotes
euk_to_keep <- c("Litostomatea_Un","Entamoeba","Trichinella","Ascaris",
                 "Trichuris","Pentatrichomonas","Hexamita","Enteromonas",
                 "Giardia","Candida","Cryptosporidium","Iodamoeba",
                 "Enterocytozoon","Sorex","Saccharomyces","Plasmodium",
                 "Entodinium","Toxocara","Unassigned_Litostomatea","Tritrichomonas",
                 "Endolimax","Malassezia","Cephalobus","Cystoisospora",
                 "Pristionchus","Enterobius","Acrostichus","Pelodera");length(euk_to_keep)
project_data_euk <- project_data_euk %>%
  subset_taxa(Rank7 %in% euk_to_keep)


##################################################################

#Calculate Alpha diversity
##################################################################

otu <- as.data.frame(t(otu_table(project_data_euk)))
Shannon_euk <- diversity(otu, index = "shannon", MARGIN = 1, base = exp(1))
Simpson_euk <- diversity(otu, "simpson", MARGIN = 1, base = exp(1))
otu[otu > 0] <- 1
Number_ASVs_euk <- rowSums(otu)
meta_shannon_euk <- cbind(Shannon_euk, Number_ASVs_euk,Simpson_euk);dim(meta_shannon)

##################################################################

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
###                                                                                                                              ###
###                                        COMPARING DATASETS TO IDENTIFY MISSING SAMPLES                                        ###
###                                                                                                                              ###
####################################################################################################################################
####################################################################################################################################
###################################################################################################################################

# Comparing data set from with all combinaison
#########################################

#Present only in Metadata
samples_diff <- setdiff(setdiff(setdiff(metadata$id,qPCR_table$Sample),row.names(asv_16s1_spl_SE_all)),colnames(asv_mi_tbl_agg))
samples_diff;length(samples_diff)

#Present only in 18s MI
samples_diff <- setdiff(setdiff(setdiff(colnames(asv_mi_tbl_agg),qPCR_table$Sample),row.names(asv_16s1_spl_SE_all)),metadata$id)
samples_diff;length(samples_diff)

#Present only in 16s MI
samples_diff <- setdiff(setdiff(setdiff(row.names(asv_16s1_spl_SE_all),qPCR_table$Sample),colnames(asv_mi_tbl_agg)),metadata$id)
samples_diff;length(samples_diff)

#Present only in qPCR
samples_diff <- setdiff(setdiff(setdiff(qPCR_table$Sample,row.names(asv_16s1_spl_SE_all)),colnames(asv_mi_tbl_agg)),metadata$id)
samples_diff;length(samples_diff)

#Present only in Metadata and 18s MI, BUT not qPCR and 16s
samples_diff <- setdiff(setdiff(intersect(colnames(asv_mi_tbl_agg),metadata$id),qPCR_table$Sample),row.names(asv_16s1_spl_SE_all))
samples_diff;length(samples_diff)

# Present in 16s, qPCR, Metadata BUT no 18s
samples_diff <- setdiff(intersect(qPCR_table$Sample,intersect(row.names(asv_16s1_spl_SE_all),metadata$id)),colnames(asv_mi_tbl_agg))
samples_diff;length(samples_diff)

# Present in 16s, 18s, Metadata BUT no qPCR
samples_diff <- setdiff(intersect(colnames(asv_mi_tbl_agg),intersect(row.names(asv_16s1_spl_SE_all),metadata$id)),qPCR_table$Sample)
samples_diff;length(samples_diff)

# Present in qPCR and Metadata BUT no 16s and 18s
samples_diff <- setdiff(setdiff(intersect(qPCR_table$Sample,metadata$id),colnames(asv_mi_tbl_agg)),row.names(asv_16s1_spl_SE_all))
samples_diff;length(samples_diff)

# Present in 16s and Metadata BUT no qPCR and 18s
samples_diff <- setdiff(setdiff(intersect(row.names(asv_16s1_spl_SE_all),metadata$id),colnames(asv_mi_tbl_agg)),qPCR_table$Sample)
samples_diff;length(samples_diff)

# Present in 18s and 16s BUT no qPCR and Metadata
samples_diff <- setdiff(setdiff(intersect(row.names(asv_16s1_spl_SE_all),colnames(asv_mi_tbl_agg)), metadata$id),qPCR_table$Sample)
samples_diff;length(samples_diff)

# Present in 18s and 16s and qPCR BUT no  Metadata
samples_diff <- intersect(setdiff(intersect(row.names(asv_16s1_spl_SE_all),colnames(asv_mi_tbl_agg)), metadata$id),qPCR_table$Sample)
samples_diff;length(samples_diff)

# Present in 18s and qPCR BUT no 16s and Metadata
samples_diff <- setdiff(setdiff(intersect(qPCR_table$Sample,colnames(asv_mi_tbl_agg)), metadata$id),row.names(asv_16s1_spl_SE_all))
samples_diff;length(samples_diff)

# Present in 18s and Metadata and qPCR BUT no 16s
samples_diff <- intersect(setdiff(intersect(metadata$id,colnames(asv_mi_tbl_agg)),row.names(asv_16s1_spl_SE_all)),qPCR_table$Sample)
samples_diff;length(samples_diff)

#Comparing 18s and metadata
unique(colnames(asv_mi_tbl_agg))
setdiff(colnames(asv_mi_tbl_agg),metadata$id)
intersect(colnames(asv_mi_tbl_agg),metadata$id)
setdiff(metadata$id, colnames(asv_mi_tbl_agg))

#Comparing 18s and qPCR
unique(colnames(asv_mi_tbl_agg))
setdiff(colnames(asv_mi_tbl_agg),plate3$Sample)
intersect(colnames(asv_mi_tbl_agg),plate3$Sample)
setdiff(plate3$Sample, colnames(asv_mi_tbl_agg))

#Comparing 18s and 16s
unique(colnames(asv_mi_tbl_agg))
setdiff(colnames(asv_mi_tbl_agg),row.names(asv_16s1_spl_SE_all))
intersect(colnames(asv_mi_tbl_agg),row.names(asv_16s1_spl_SE_all))
setdiff(row.names(asv_16s1_spl_SE_all), colnames(asv_mi_tbl_agg))

#Comparing 16s and metadata
unique(row.names(asv_16s1_spl_SE_all))
setdiff(row.names(asv_16s1_spl_SE_all),metadata$id)
intersect(row.names(asv_16s1_spl_SE_all),metadata$id)
setdiff(metadata$id, row.names(asv_16s1_spl_SE_all))

#Comparing 16s and qPCR
unique(colnames(asv_mi_tbl_agg))
setdiff(row.names(asv_16s1_spl_SE_all),metadata$id)
intersect(row.names(asv_16s1_spl_SE_all),metadata$id)
setdiff(metadata$id, row.names(asv_16s1_spl_SE_all))

#Comparing qPCR and metadata
unique(plate3$Sample)
setdiff(plate3$Sample,metadata$id)
intersect(plate3$Sample,metadata$id)
setdiff(metadata$id, plate3$Sample)

#########################################

# Final Venn Diagram
#########################################

x1 = list(Euk_18S_MI = unique(colnames(asv_mi_tbl_agg)), 
          Metadata = metadata$id,
          Bact_16s_MI = unique(row.names(asv_16s1_spl_SE_all)),
          Blasto_qPCR = unique(qPCR_table$Sample))


pdf("/Users/vincebilly/Desktop/venn.diagram_comparison_dataset.pdf")
ggVennDiagram(x1, label_alpha = 0,
              stroke_size = 0.1, set_name_size = 1) + 
  scale_fill_gradient(low="white",high = "red")
dev.off()
#########################################

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
###                                                                                                                              ###
###                                    INFLUENCE OF FACTOR ON BLASTOCYSTIS PRESENCE AND ABUNDANCE                                ###
###                                                                                                                              ###
####################################################################################################################################
####################################################################################################################################
###################################################################################################################################

##################################################################
###               COMBINE qPCR and METADATA                    ###
##################################################################

#Check if samples match with MISEQ data
qPCR_Metadata1 <- merge(plate3,metadata, by.x = "Sample", by.y = "id", all= FALSE);dim(qPCR_Metadata1)
qPCR_Metadata <- merge(qPCR_Metadata1, meta_shannon16s_new, by.x = "Sample", by.y = "row.names");dim(qPCR_Metadata1)
#qPCR_Metadata <- merge(qPCR_Metadata1, meta_shannon_euk, by.x = "Sample", by.y = "row.names", all= FALSE);dim(qPCR_Metadata)


qPCR_Metadata$Cq <- as.numeric(qPCR_Metadata$Cq)
qPCR_Metadata$Cq_presence_absence <- ifelse(qPCR_Metadata$Cq<=35 & qPCR_Metadata$Cq>10,"Presence_qPCR", "Absence_qPCR")

#Combined qPCR and metadata
qPCR_Metadata$Cq_presence_absence <- as.factor(qPCR_Metadata$Cq_presence_absence)
qPCR_Metadata <- qPCR_Metadata %>%
  reorder_levels(Cq_presence_absence, order = c("Presence_qPCR", "Absence_qPCR"))
levels(qPCR_Metadata$Cq_presence_absence)

#Reanotate variable
qPCR_Metadata$stunted <- gsub("Oui","Stunted",qPCR_Metadata$stunted)
qPCR_Metadata$stunted <- gsub("Non","Non-Stunted",qPCR_Metadata$stunted)
qPCR_Metadata$stunted <- as.factor(qPCR_Metadata$stunted)
qPCR_Metadata <- qPCR_Metadata %>%
  reorder_levels(stunted, order = c("Stunted", "Non-Stunted"))
levels(qPCR_Metadata$stunted)



##################################################################

##################################################################
## Association pays stunted   ###
##################################################################


table <- table(qPCR_Metadata$pays,qPCR_Metadata$stunted)
fisher.test(table)

table <- as.data.frame(table)
colnames(table) <- c("Country","Stunted","Number_samples")

colours <- c("darkgreen","red")

ggplot(table, aes(x=Country, y =Number_samples , fill=Stunted, color = Stunted )) +
  geom_bar(stat="identity", width=1, color="black", linetype = "solid") +  #height of the column equal the value, stacked
  scale_fill_manual(values=colours) +
  theme_bw() 

ggplot(table, aes(x=qPCR, y=Freq, fill=Miseq, color = colours )) +
  geom_bar(stat="identity", width=1, color="black", linetype = "solid") +  #height of the column equal the value, stacked
  scale_fill_manual(values=colours) +
  theme_bw() 

ggplot(table, aes(x = ))


table <- table(qPCR_Metadata$pays,qPCR_Metadata$Cq_presence_absence)
fisher.test(table)

table <- table(qPCR_Metadata$stunted,qPCR_Metadata$Cq_presence_absence)
fisher.test(table)

##################################################################

##################################################################
###                      GLOBAL PREVALENCE                     ###
##################################################################

#qPCR_Metadata2 <- qPCR_Metadata
qPCR_Metadata$Cq_trans <- replace(qPCR_Metadata$Cq, qPCR_Metadata$Cq<1,45)
plot <- ggplot(qPCR_Metadata,aes(x=Cq_trans)) +
  geom_histogram(binwidth = 1, boundary = 0, color ="White") +
  scale_x_continuous(breaks = c(0,seq(0, 55, 5))) +
  geom_vline(xintercept = 40, linetype="dotted", color = "blue", linewidth=1) +
  geom_vline(xintercept = 35, linetype="dotted", color = "red", linewidth=1) +
  theme_light() +
  theme(axis.title=element_text(size=15, face="bold"),
        axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        plot.title=element_text(size=15, face="bold")) +
  labs(y = "Number of fecal DNA samples", 
       x = "Ct values of qPCR reactions")
print(plot)


#Shapiro test
shapiro.test(qPCR_Metadata$Cq)
#Not normal distributed

data_prevalence_threshold_list <- list()
cq_threshold <- seq(20,40,0.1)
for(i in cq_threshold) {
  qPCR_table2_prev <- qPCR_Metadata
  qPCR_table2_prev$Cq_presence_absence <- ifelse(qPCR_table2_prev$Cq<=i & qPCR_Metadata$Cq>10,"Presence_qPCR", "Absence_qPCR")
  
  prevalence <- table(qPCR_table2_prev$Cq_presence_absence)[2]/sum(table(qPCR_table2_prev$Cq_presence_absence))*100
  nb_infected <- table(qPCR_table2_prev$Cq_presence_absence)[2]
  summary_prevalence_nb_infected.df <- data.frame(Cq_threshold = i,
                                                  Prevalence=prevalence,
                                                  Nb_infected=nb_infected)
  
  data_prevalence_threshold_list[[length(data_prevalence_threshold_list)+1]] <- summary_prevalence_nb_infected.df
}


# Generate table that concatenate all value
summary_stat <- data.frame()
for (i in 1:length(data_prevalence_threshold_list)){
  summary_stat <- rbind(summary_stat,data_prevalence_threshold_list[[i]])
}

scatter_plot <- ggplot(summary_stat, aes(x=Cq_threshold)) +
  geom_point( aes(y=Prevalence,color = "Prevalence")) + 
  geom_point( aes(y=Nb_infected/5,color = "Nb_infected")) + # Divide by 10 to get the same range than the temperature
  scale_color_manual(values = c("Prevalence" = "gray40", "Nb_infected" = "gray1")) +
  scale_y_continuous(
    # Features of the first axis
    name = "Prevalence",
    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*5, name="Number of infected children"),
    breaks = c(0,seq(0, 100, 5))) +
  geom_vline(xintercept = 40, linetype="dotted", color = "blue", linewidth=1) +
  geom_vline(xintercept = 35, linetype="dotted", color = "red", linewidth=1) +
  theme_light() +
  theme(axis.title.y = element_text(color = "gray40", size=15, face="bold"),
        axis.text.y = element_text(color = "gray40", size=15),
        axis.title.y.right = element_text(color = "black", size=15, face="bold"),
        axis.text.y.right = element_text(color = "black", size=15),
        axis.text.x=element_text(size=15),
        axis.title.x=element_text(size=15, face="bold"),
        legend.position = "none") +
  labs(x = "Ct values cutoff for presence of Blastocystis")

print(scatter_plot)



plot_combined <- ggarrange(plot,scatter_plot,labels = c("A","B"), ncol = 2, nrow=1);plot_combined
ggsave2(plot_combined, filename = "/Users/vincebilly/Desktop/vince/phd/statistic/afribiota_vince/Simple_risk_factor/figure/Figure1_prevalence_distribution_Cq.pdf", height = 4, width = 10)


##################################################################

##################################################################
### NORMALITY
##################################################################



##################################################################

##################################################################
###          ANOVA          ###
##################################################################


qPCR_Metadata <- qPCR_Metadata %>%
  reorder_levels(stunted, order = c("Stunted", "Non-Stunted"))
qPCR_Metadata <- qPCR_Metadata %>%
  reorder_levels(pays, order = c("Madagascar", "RCA"))
qPCR_Metadata <- qPCR_Metadata %>%
  reorder_levels(Cq_presence_absence, order = c("Presence_qPCR", "Absence_qPCR"))


##################################################################


##################################################################
## ONE WAY ANOVA SHANNON stunted
##################################################################

qPCR_Metadata %>% sample_n_by(stunted, size = 1)
levels(qPCR_Metadata$stunted)
qPCR_Metadata <- qPCR_Metadata %>%
  reorder_levels(stunted, order = c("Non", "Oui"))
qPCR_Metadata %>%
  group_by(stunted) %>%
  get_summary_stats(Shannon, type = "mean_sd")
ggboxplot(qPCR_Metadata, x = "stunted", y = "Shannon")
qPCR_Metadata %>% 
  group_by(stunted) %>%
  identify_outliers(Shannon)
# Build the linear model
model  <- lm(Shannon ~ stunted, data = qPCR_Metadata)
# Create a QQ plot of residuals
ggqqplot(residuals(model))

# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))

qPCR_Metadata %>%
  group_by(stunted) %>%
  shapiro_test(Shannon)

ggqqplot(qPCR_Metadata, "Shannon", facet.by = "stunted")


plot(model, 1)


qPCR_Metadata %>% levene_test(Shannon ~ stunted)


res.aov <- qPCR_Metadata %>% anova_test(Shannon ~ stunted)
res.aov


# Pairwise comparisons
pwc <- qPCR_Metadata %>% tukey_hsd(Shannon ~ stunted)
pwc


# Visualization: box plots with p-values
pwc <- pwc %>% add_xy_position(x = "stunted")
ggboxplot(qPCR_Metadata, x = "stunted", y = "Shannon",
          color = "stunted", palette = "jco") +
  stat_pvalue_manual(pwc, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(res.aov, detailed = TRUE),
    caption = get_pwc_label(pwc)
  )


# Welch One way ANOVA test
res.aov2 <- qPCR_Metadata %>% welch_anova_test(Shannon ~ stunted)
# Pairwise comparisons (Games-Howell)
pwc2 <- qPCR_Metadata %>% games_howell_test(Shannon ~ stunted)
# Visualization: box plots with p-values
pwc2 <- pwc2 %>% add_xy_position(x = "stunted", step.increase = 1)
ggboxplot(qPCR_Metadata, x = "stunted", y = "Shannon") +
  stat_pvalue_manual(pwc2, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(res.aov2, detailed = TRUE),
    caption = get_pwc_label(pwc2)
  )


pwc3 <- qPCR_Metadata %>% 
  pairwise_t_test(
    Shannon ~ stunted, pool.sd = FALSE,
    p.adjust.method = "bonferroni"
  )
pwc3
##################################################################

##################################################################
## ONE WAY ANOVA Number_ASVs stunted
##################################################################

qPCR_Metadata %>% sample_n_by(stunted, size = 1)
levels(qPCR_Metadata$stunted)
qPCR_Metadata <- qPCR_Metadata %>%
  reorder_levels(stunted, order = c("Non", "Oui"))
qPCR_Metadata %>%
  group_by(stunted) %>%
  get_summary_stats(Number_ASVs, type = "mean_sd")
ggboxplot(qPCR_Metadata, x = "stunted", y = "Number_ASVs")
qPCR_Metadata %>% 
  group_by(stunted) %>%
  identify_outliers(Number_ASVs)
# Build the linear model
model  <- lm(Number_ASVs ~ stunted, data = qPCR_Metadata)
# Create a QQ plot of residuals
ggqqplot(residuals(model))

# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))

qPCR_Metadata %>%
  group_by(stunted) %>%
  shapiro_test(Number_ASVs)

ggqqplot(qPCR_Metadata, "Number_ASVs", facet.by = "stunted")


plot(model, 1)


qPCR_Metadata %>% levene_test(Number_ASVs ~ stunted)


res.aov <- qPCR_Metadata %>% anova_test(Number_ASVs ~ stunted)
res.aov


# Pairwise comparisons
pwc <- qPCR_Metadata %>% tukey_hsd(Number_ASVs ~ stunted)
pwc


# Visualization: box plots with p-values
pwc <- pwc %>% add_xy_position(x = "stunted")
ggboxplot(qPCR_Metadata, x = "stunted", y = "Number_ASVs",
          color = "stunted", palette = "jco") +
  stat_pvalue_manual(pwc, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(res.aov, detailed = TRUE),
    caption = get_pwc_label(pwc)
  )


# Welch One way ANOVA test
res.aov2 <- qPCR_Metadata %>% welch_anova_test(Number_ASVs ~ stunted)
# Pairwise comparisons (Games-Howell)
pwc2 <- qPCR_Metadata %>% games_howell_test(Number_ASVs ~ stunted)
# Visualization: box plots with p-values
pwc2 <- pwc2 %>% add_xy_position(x = "stunted", step.increase = 1)
ggboxplot(qPCR_Metadata, x = "stunted", y = "Number_ASVs") +
  stat_pvalue_manual(pwc2, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(res.aov2, detailed = TRUE),
    caption = get_pwc_label(pwc2)
  )


pwc3 <- qPCR_Metadata %>% 
  pairwise_t_test(
    Number_ASVs ~ stunted, pool.sd = FALSE,
    p.adjust.method = "bonferroni"
  )
pwc3
##################################################################

##################################################################
## ONE WAY ANOVA SHANNON
##################################################################

qPCR_Metadata %>% sample_n_by(Cq_presence_absence, size = 1)
levels(qPCR_Metadata$Cq_presence_absence)
qPCR_Metadata <- qPCR_Metadata %>%
  reorder_levels(Cq_presence_absence, order = c("Presence_qPCR", "Absence_qPCR"))
qPCR_Metadata %>%
  group_by(Cq_presence_absence) %>%
  get_summary_stats(Shannon, type = "mean_sd")
ggboxplot(qPCR_Metadata, x = "Cq_presence_absence", y = "Shannon")
qPCR_Metadata %>% 
  group_by(Cq_presence_absence) %>%
  identify_outliers(Shannon)
# Build the linear model
model  <- lm(Shannon ~ Cq_presence_absence, data = qPCR_Metadata)

# Create a QQ plot of residuals
ggqqplot(residuals(model))

# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))

qPCR_Metadata %>%
  group_by(Cq_presence_absence) %>%
  shapiro_test(Shannon)

ggqqplot(qPCR_Metadata, "Shannon", facet.by = "Cq_presence_absence")


plot(model, 1)


qPCR_Metadata %>% levene_test(Shannon ~ Cq_presence_absence)


res.aov <- qPCR_Metadata %>% anova_test(Shannon ~ Cq_presence_absence)
res.aov


# Pairwise comparisons
pwc <- qPCR_Metadata %>% tukey_hsd(Shannon ~ Cq_presence_absence)
pwc


# Visualization: box plots with p-values
pwc <- pwc %>% add_xy_position(x = "Cq_presence_absence")
ggboxplot(qPCR_Metadata, x = "Cq_presence_absence", y = "Shannon",
          color = "Cq_presence_absence", palette = "jco") +
  stat_pvalue_manual(pwc, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(res.aov, detailed = TRUE),
    caption = get_pwc_label(pwc)
  )


# Welch One way ANOVA test
res.aov2 <- qPCR_Metadata %>% welch_anova_test(Shannon ~ Cq_presence_absence)
# Pairwise comparisons (Games-Howell)
pwc2 <- qPCR_Metadata %>% games_howell_test(Shannon ~ Cq_presence_absence)
# Visualization: box plots with p-values
pwc2 <- pwc2 %>% add_xy_position(x = "Cq_presence_absence", step.increase = 1)
ggboxplot(qPCR_Metadata, x = "Cq_presence_absence", y = "Shannon") +
  stat_pvalue_manual(pwc2, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(res.aov2, detailed = TRUE),
    caption = get_pwc_label(pwc2)
  )


pwc3 <- qPCR_Metadata %>% 
  pairwise_t_test(
    Shannon ~ Cq_presence_absence, pool.sd = FALSE,
    p.adjust.method = "bonferroni"
  )
pwc3
##################################################################

##################################################################
## ONE WAY ANOVA Number_ASVs
##################################################################

qPCR_Metadata %>% sample_n_by(Cq_presence_absence, size = 1)
levels(qPCR_Metadata$Cq_presence_absence)
qPCR_Metadata <- qPCR_Metadata %>%
  reorder_levels(Cq_presence_absence, order = c("Presence_qPCR", "Absence_qPCR"))
qPCR_Metadata %>%
  group_by(Cq_presence_absence) %>%
  get_summary_stats(Number_ASVs, type = "mean_sd")
ggboxplot(qPCR_Metadata, x = "Cq_presence_absence", y = "Number_ASVs")
qPCR_Metadata %>% 
  group_by(Cq_presence_absence) %>%
  identify_outliers(Number_ASVs)
# Build the linear model
model  <- lm(Number_ASVs ~ Cq_presence_absence, data = qPCR_Metadata)
# Create a QQ plot of residuals
ggqqplot(residuals(model))

# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))

qPCR_Metadata %>%
  group_by(Cq_presence_absence) %>%
  shapiro_test(Number_ASVs)

ggqqplot(qPCR_Metadata, "Number_ASVs", facet.by = "Cq_presence_absence")


plot(model, 1)


qPCR_Metadata %>% levene_test(Number_ASVs ~ Cq_presence_absence)


res.aov <- qPCR_Metadata %>% anova_test(Number_ASVs ~ Cq_presence_absence)
res.aov


# Pairwise comparisons
pwc <- qPCR_Metadata %>% tukey_hsd(Number_ASVs ~ Cq_presence_absence)
pwc


# Visualization: box plots with p-values
pwc <- pwc %>% add_xy_position(x = "Cq_presence_absence")
ggboxplot(qPCR_Metadata, x = "Cq_presence_absence", y = "Number_ASVs",
          color = "Cq_presence_absence", palette = "jco") +
  stat_pvalue_manual(pwc, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(res.aov, detailed = TRUE),
    caption = get_pwc_label(pwc)
  )


# Welch One way ANOVA test
res.aov2 <- qPCR_Metadata %>% welch_anova_test(Number_ASVs ~ Cq_presence_absence)
# Pairwise comparisons (Games-Howell)
pwc2 <- qPCR_Metadata %>% games_howell_test(Number_ASVs ~ Cq_presence_absence)
# Visualization: box plots with p-values
pwc2 <- pwc2 %>% add_xy_position(x = "Cq_presence_absence", step.increase = 1)
ggboxplot(qPCR_Metadata, x = "Cq_presence_absence", y = "Number_ASVs") +
  stat_pvalue_manual(pwc2, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(res.aov2, detailed = TRUE),
    caption = get_pwc_label(pwc2)
  )


pwc3 <- qPCR_Metadata %>% 
  pairwise_t_test(
    Number_ASVs ~ Cq_presence_absence, pool.sd = FALSE,
    p.adjust.method = "bonferroni"
  )
pwc3
##################################################################

##################################################################
### 2 WAYS ANOVA Country*stunted ~ shannon
##################################################################

qPCR_Metadata %>% sample_n_by(stunted, pays, size = 1)
qPCR_Metadata$stunted <- as.factor(qPCR_Metadata$stunted)
levels(qPCR_Metadata$stunted)
qPCR_Metadata <- qPCR_Metadata %>%
  reorder_levels(stunted, order = c("Non-Stunted", "Non-Stunted"))
qPCR_Metadata <- qPCR_Metadata %>%
  reorder_levels(pays, order = c("Madagascar", "RCA"))

qPCR_Metadata %>%
  group_by(pays, stunted) %>%
  get_summary_stats(Shannon16s_new.x, type = "mean_sd")


bxp <- ggboxplot(
  qPCR_Metadata, x = "pays", y = "Shannon16s_new.x",
  color = "stunted", palette = "jco"
)
bxp


qPCR_Metadata %>%
  group_by(stunted, pays) %>%
  identify_outliers(Shannon16s_new.x)


# Build the linear model
model  <- lm(Shannon16s_new.x ~ stunted*pays,
             data = qPCR_Metadata)
# Create a QQ plot of residuals
ggqqplot(residuals(model))


# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))


qPCR_Metadata %>%
  group_by(stunted, pays) %>%
  shapiro_test(Shannon16s_new.x)


ggqqplot(qPCR_Metadata, "Shannon16s_new.x", ggtheme = theme_bw()) +
  facet_grid(stunted ~ pays)


qPCR_Metadata %>% levene_test(Shannon16s_new.x ~ stunted*pays)


res.aov <- qPCR_Metadata %>% anova_test(Shannon16s_new.x ~ stunted * pays)
res.aov


# Group the data by stunted and fit  anova
model <- lm(Shannon16s_new.x ~ stunted * pays, data = qPCR_Metadata)
qPCR_Metadata %>%
  group_by(pays) %>%
  anova_test(Shannon16s_new.x ~ stunted , error = model)


# pairwise comparisons
library(emmeans)
pwc <- qPCR_Metadata %>% 
  group_by(pays ) %>%
  emmeans_test(Shannon16s_new.x ~ stunted, p.adjust.method = "bonferroni") 
pwc


res.aov


qPCR_Metadata %>%
  pairwise_t_test(
    Shannon16s_new.x ~ stunted, 
    p.adjust.method = "bonferroni"
  )



model <- lm(Shannon16s_new.x ~ stunted * pays, data = qPCR_Metadata)
qPCR_Metadata %>% 
  emmeans_test(
    Shannon16s_new.x ~ stunted, p.adjust.method = "bonferroni",
    model = model
  )



# Visualization: box plots with p-values
pwc <- pwc %>% add_xy_position(x = "pays")
bxp +
  stat_pvalue_manual(pwc) +
  labs(
    subtitle = get_test_label(res.aov, detailed = TRUE),
    caption = get_pwc_label(pwc)
  )

##################################################################

##################################################################
### 2 WAYS ANOVA Country*stunted ~ Number ASVS
##################################################################

qPCR_Metadata %>% sample_n_by(stunted, pays, size = 1)
levels(qPCR_Metadata$stunted)
qPCR_Metadata <- qPCR_Metadata %>%
  reorder_levels(stunted, order = c("Non", "Oui"))
qPCR_Metadata <- qPCR_Metadata %>%
  reorder_levels(pays, order = c("Madagascar", "RCA"))

qPCR_Metadata %>%
  group_by(pays, stunted) %>%
  get_summary_stats(Number_ASVs, type = "mean_sd")


bxp <- ggboxplot(
  qPCR_Metadata, x = "pays", y = "Number_ASVs",
  color = "stunted", palette = "jco"
)
bxp


qPCR_Metadata %>%
  group_by(stunted, pays) %>%
  identify_outliers(Number_ASVs)


# Build the linear model
model  <- lm(Number_ASVs ~ stunted*pays,
             data = qPCR_Metadata)
# Create a QQ plot of residuals
ggqqplot(residuals(model))


# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))


qPCR_Metadata %>%
  group_by(stunted, pays) %>%
  shapiro_test(Number_ASVs)


ggqqplot(qPCR_Metadata, "Number_ASVs", ggtheme = theme_bw()) +
  facet_grid(stunted ~ pays)


qPCR_Metadata %>% levene_test(Number_ASVs ~ stunted*pays)


res.aov <- qPCR_Metadata %>% anova_test(Number_ASVs ~ stunted * pays)
res.aov


# Group the data by stunted and fit  anova
model <- lm(Number_ASVs ~ stunted * pays, data = qPCR_Metadata)
qPCR_Metadata %>%
  group_by(pays) %>%
  anova_test(Number_ASVs ~ stunted , error = model)


# pairwise comparisons
library(emmeans)
pwc <- qPCR_Metadata %>% 
  group_by(pays ) %>%
  emmeans_test(Number_ASVs ~ stunted, p.adjust.method = "bonferroni") 
pwc


res.aov


qPCR_Metadata %>%
  pairwise_t_test(
    Number_ASVs ~ stunted, 
    p.adjust.method = "bonferroni"
  )



model <- lm(Number_ASVs ~ stunted * pays, data = qPCR_Metadata)
qPCR_Metadata %>% 
  emmeans_test(
    Number_ASVs ~ stunted, p.adjust.method = "bonferroni",
    model = model
  )



# Visualization: box plots with p-values
pwc <- pwc %>% add_xy_position(x = "pays")
bxp +
  stat_pvalue_manual(pwc) +
  labs(
    subtitle = get_test_label(res.aov, detailed = TRUE),
    caption = get_pwc_label(pwc)
  )

##################################################################

##################################################################
### 2 WAYS ANOVA Country Number_ASVs
##################################################################

qPCR_Metadata %>% sample_n_by(Cq_presence_absence, stunted, size = 1)


qPCR_Metadata %>%
  group_by(pays, stunted) %>%
  get_summary_stats(Number_ASVs, type = "mean_sd")


bxp <- ggboxplot(
  qPCR_Metadata, x = "pays", y = "Number_ASVs",
  color = "Cq_presence_absence", palette = "jco"
)
bxp


qPCR_Metadata %>%
  group_by(Cq_presence_absence, pays) %>%
  identify_outliers(Number_ASVs)


# Build the linear model
model  <- lm(Number_ASVs ~ Cq_presence_absence*pays,
             data = qPCR_Metadata)
# Create a QQ plot of residuals
ggqqplot(residuals(model))


# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))


qPCR_Metadata %>%
  group_by(Cq_presence_absence, pays) %>%
  shapiro_test(Number_ASVs)


ggqqplot(qPCR_Metadata, "Number_ASVs", ggtheme = theme_bw()) +
  facet_grid(Cq_presence_absence ~ pays)


qPCR_Metadata %>% levene_test(Number_ASVs ~ Cq_presence_absence*pays)


res.aov <- qPCR_Metadata %>% anova_test(Number_ASVs ~ Cq_presence_absence * pays)
res.aov


# Group the data by Cq_presence_absence and fit  anova
model <- lm(Number_ASVs ~ Cq_presence_absence * pays, data = qPCR_Metadata)
qPCR_Metadata %>%
  group_by(pays) %>%
  anova_test(Number_ASVs ~ Cq_presence_absence , error = model)


# pairwise comparisons
library(emmeans)
pwc <- qPCR_Metadata %>% 
  group_by(pays ) %>%
  emmeans_test(Number_ASVs ~ Cq_presence_absence, p.adjust.method = "bonferroni") 

pwc


res.aov


qPCR_Metadata %>%
  pairwise_t_test(
    Number_ASVs ~ Cq_presence_absence, 
    p.adjust.method = "bonferroni"
  )



model <- lm(Number_ASVs ~ Cq_presence_absence * pays, data = qPCR_Metadata)
qPCR_Metadata %>% 
  emmeans_test(
    Number_ASVs ~ Cq_presence_absence, p.adjust.method = "bonferroni",
    model = model
  )



# Visualization: box plots with p-values
pwc <- pwc %>% add_xy_position(x = "pays")
bxp +
  stat_pvalue_manual(pwc) +
  labs(
    subtitle = get_test_label(res.aov, detailed = TRUE),
    caption = get_pwc_label(pwc)
  )

##################################################################

##################################################################
### 2 WAYS ANOVA Country shannon
##################################################################

qPCR_Metadata %>% sample_n_by(Cq_presence_absence, stunted, size = 1)


qPCR_Metadata %>%
  group_by(pays, stunted) %>%
  get_summary_stats(Shannon, type = "mean_sd")


bxp <- ggboxplot(
  qPCR_Metadata, x = "pays", y = "Shannon",
  color = "Cq_presence_absence", palette = "jco"
)
bxp


qPCR_Metadata %>%
  group_by(Cq_presence_absence, pays) %>%
  identify_outliers(Shannon)


# Build the linear model
model  <- lm(Shannon ~ Cq_presence_absence*pays,
             data = qPCR_Metadata)
# Create a QQ plot of residuals
ggqqplot(residuals(model))


# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))


qPCR_Metadata %>%
  group_by(Cq_presence_absence, pays) %>%
  shapiro_test(Shannon)


ggqqplot(qPCR_Metadata, "Shannon", ggtheme = theme_bw()) +
  facet_grid(Cq_presence_absence ~ pays)


qPCR_Metadata %>% levene_test(Shannon ~ Cq_presence_absence*pays)


res.aov <- qPCR_Metadata %>% anova_test(Shannon ~ Cq_presence_absence * pays)
res.aov


# Group the data by Cq_presence_absence and fit  anova
model <- lm(Shannon ~ Cq_presence_absence * pays, data = qPCR_Metadata)
qPCR_Metadata %>%
  group_by(pays) %>%
  anova_test(Shannon ~ Cq_presence_absence , error = model)


# pairwise comparisons
library(emmeans)
pwc <- qPCR_Metadata %>% 
  group_by(pays ) %>%
  emmeans_test(Shannon ~ Cq_presence_absence, p.adjust.method = "bonferroni") 
pwc


res.aov


qPCR_Metadata %>%
  pairwise_t_test(
    Shannon ~ Cq_presence_absence, 
    p.adjust.method = "bonferroni"
  )



model <- lm(Shannon ~ Cq_presence_absence * pays, data = qPCR_Metadata)
qPCR_Metadata %>% 
  emmeans_test(
    Shannon ~ Cq_presence_absence, p.adjust.method = "bonferroni",
    model = model
  )



# Visualization: box plots with p-values
pwc <- pwc %>% add_xy_position(x = "pays")
bxp +
  stat_pvalue_manual(pwc) +
  labs(
    subtitle = get_test_label(res.aov, detailed = TRUE),
    caption = get_pwc_label(pwc)
  )

##################################################################

##################################################################
### 2 WAYS ANOVA Country Number_ASVs
##################################################################

qPCR_Metadata %>% sample_n_by(Cq_presence_absence, stunted, size = 1)


qPCR_Metadata %>%
  group_by(pays, stunted) %>%
  get_summary_stats(Number_ASVs, type = "mean_sd")


bxp <- ggboxplot(
  qPCR_Metadata, x = "pays", y = "Number_ASVs",
  color = "Cq_presence_absence", palette = "jco"
)
bxp


qPCR_Metadata %>%
  group_by(Cq_presence_absence, pays) %>%
  identify_outliers(Number_ASVs)


# Build the linear model
model  <- lm(Number_ASVs ~ Cq_presence_absence*pays,
             data = qPCR_Metadata)
# Create a QQ plot of residuals
ggqqplot(residuals(model))


# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))


qPCR_Metadata %>%
  group_by(Cq_presence_absence, pays) %>%
  shapiro_test(Number_ASVs)


ggqqplot(qPCR_Metadata, "Number_ASVs", ggtheme = theme_bw()) +
  facet_grid(Cq_presence_absence ~ pays)


qPCR_Metadata %>% levene_test(Number_ASVs ~ Cq_presence_absence*pays)


res.aov <- qPCR_Metadata %>% anova_test(Number_ASVs ~ Cq_presence_absence * pays)
res.aov


# Group the data by Cq_presence_absence and fit  anova
model <- lm(Number_ASVs ~ Cq_presence_absence * pays, data = qPCR_Metadata)
qPCR_Metadata %>%
  group_by(pays) %>%
  anova_test(Number_ASVs ~ Cq_presence_absence , error = model)


# pairwise comparisons
library(emmeans)
pwc <- qPCR_Metadata %>% 
  group_by(pays ) %>%
  emmeans_test(Number_ASVs ~ Cq_presence_absence, p.adjust.method = "bonferroni") 

pwc


res.aov


qPCR_Metadata %>%
  pairwise_t_test(
    Number_ASVs ~ Cq_presence_absence, 
    p.adjust.method = "bonferroni"
  )



model <- lm(Number_ASVs ~ Cq_presence_absence * pays, data = qPCR_Metadata)
qPCR_Metadata %>% 
  emmeans_test(
    Number_ASVs ~ Cq_presence_absence, p.adjust.method = "bonferroni",
    model = model
  )



# Visualization: box plots with p-values
pwc <- pwc %>% add_xy_position(x = "pays")
bxp +
  stat_pvalue_manual(pwc) +
  labs(
    subtitle = get_test_label(res.aov, detailed = TRUE),
    caption = get_pwc_label(pwc)
  )

##################################################################

##################################################################
### 2 WAYS ANOVA Stunted Shannon
##################################################################

qPCR_Metadata %>% sample_n_by(Cq_presence_absence, stunted, size = 1)
levels(qPCR_Metadata$stunted)
qPCR_Metadata <- qPCR_Metadata %>%
  reorder_levels(stunted, order = c("Non", "Oui"))


qPCR_Metadata %>%
  group_by(pays, stunted) %>%
  get_summary_stats(Shannon, type = "mean_sd")


bxp <- ggboxplot(
  qPCR_Metadata, x = "stunted", y = "Shannon",
  color = "Cq_presence_absence", palette = "jco"
)
bxp


qPCR_Metadata %>%
  group_by(Cq_presence_absence, stunted) %>%
  identify_outliers(Shannon)


# Build the linear model
model  <- lm(Shannon ~ Cq_presence_absence*stunted,
             data = qPCR_Metadata)
# Create a QQ plot of residuals
ggqqplot(residuals(model))


# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))


qPCR_Metadata %>%
  group_by(Cq_presence_absence, stunted) %>%
  shapiro_test(Shannon)


ggqqplot(qPCR_Metadata, "Shannon", ggtheme = theme_bw()) +
  facet_grid(Cq_presence_absence ~ stunted)


qPCR_Metadata %>% levene_test(Shannon ~ Cq_presence_absence*stunted)


res.aov <- qPCR_Metadata %>% anova_test(Shannon ~ Cq_presence_absence * stunted)
res.aov


# Group the data by Cq_presence_absence and fit  anova
model <- lm(Shannon ~ Cq_presence_absence * stunted, data = qPCR_Metadata)
qPCR_Metadata %>%
  group_by(stunted) %>%
  anova_test(Shannon ~ Cq_presence_absence , error = model)


# pairwise comparisons
library(emmeans)
pwc <- qPCR_Metadata %>% 
  group_by(stunted ) %>%
  emmeans_test(Shannon ~ Cq_presence_absence, p.adjust.method = "bonferroni") 
pwc


res.aov


qPCR_Metadata %>%
  pairwise_t_test(
    Shannon ~ Cq_presence_absence, 
    p.adjust.method = "bonferroni"
  )



model <- lm(Shannon ~ Cq_presence_absence * stunted, data = qPCR_Metadata)
qPCR_Metadata %>% 
  emmeans_test(
    Shannon ~ Cq_presence_absence, p.adjust.method = "bonferroni",
    model = model
  )



# Visualization: box plots with p-values
pwc <- pwc %>% add_xy_position(x = "stunted")
bxp +
  stat_pvalue_manual(pwc) +
  labs(
    subtitle = get_test_label(res.aov, detailed = TRUE),
    caption = get_pwc_label(pwc)
  )

##################################################################

##################################################################
### 2 WAYS ANOVA Stunted Number of ASVs
##################################################################

qPCR_Metadata %>% sample_n_by(Cq_presence_absence, stunted, size = 1)
levels(qPCR_Metadata$stunted)
qPCR_Metadata <- qPCR_Metadata %>%
  reorder_levels(stunted, order = c("Non", "Oui"))


qPCR_Metadata %>%
  group_by(pays, stunted) %>%
  get_summary_stats(Number_ASVs, type = "mean_sd")


bxp <- ggboxplot(
  qPCR_Metadata, x = "stunted", y = "Number_ASVs",
  color = "Cq_presence_absence", palette = "jco"
)
bxp


qPCR_Metadata %>%
  group_by(Cq_presence_absence, stunted) %>%
  identify_outliers(Number_ASVs)


# Build the linear model
model  <- lm(Number_ASVs ~ Cq_presence_absence*stunted,
             data = qPCR_Metadata)
# Create a QQ plot of residuals
ggqqplot(residuals(model))


# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))


qPCR_Metadata %>%
  group_by(Cq_presence_absence, stunted) %>%
  shapiro_test(Number_ASVs)


ggqqplot(qPCR_Metadata, "Number_ASVs", ggtheme = theme_bw()) +
  facet_grid(Cq_presence_absence ~ stunted)


qPCR_Metadata %>% levene_test(Number_ASVs ~ Cq_presence_absence*stunted)


res.aov <- qPCR_Metadata %>% anova_test(Number_ASVs ~ Cq_presence_absence * stunted)
res.aov


# Group the data by Cq_presence_absence and fit  anova
model <- lm(Number_ASVs ~ Cq_presence_absence * stunted, data = qPCR_Metadata)
qPCR_Metadata %>%
  group_by(stunted) %>%
  anova_test(Number_ASVs ~ Cq_presence_absence , error = model)


# pairwise comparisons
library(emmeans)
pwc <- qPCR_Metadata %>% 
  group_by(stunted ) %>%
  emmeans_test(Number_ASVs ~ Cq_presence_absence, p.adjust.method = "bonferroni") 
pwc


res.aov


qPCR_Metadata %>%
  pairwise_t_test(
    Number_ASVs ~ Cq_presence_absence, 
    p.adjust.method = "bonferroni"
  )



model <- lm(Number_ASVs ~ Cq_presence_absence * stunted, data = qPCR_Metadata)
qPCR_Metadata %>% 
  emmeans_test(
    Number_ASVs ~ Cq_presence_absence, p.adjust.method = "bonferroni",
    model = model
  )



# Visualization: box plots with p-values
pwc <- pwc %>% add_xy_position(x = "stunted")
bxp +
  stat_pvalue_manual(pwc) +
  labs(
    subtitle = get_test_label(res.aov, detailed = TRUE),
    caption = get_pwc_label(pwc)
  )

##################################################################

##################################################################
###   3-WAY ANOVA with LOOP for ALL VARIABLE
##################################################################

var1_vector <- c("pays","sexe","diarrhee","voie_accou","stunted","eau_fontaine","dermatoses","mere_emploi","pere_emploi","pere_rev","murs","electricite","eau","cuisine")
var2_vector <- var1_vector
LIST_3WAY_ANOVA <- list()
for (vr1 in seq_along(var1_vector)) {
  for (vr2 in seq_along(var2_vector)) {
    if(var1_vector[vr1] !=var2_vector[vr2]) {
      pwc <- qPCR_Metadata %>%
        group_by(!!parse_expr(var1_vector[vr1]),!!parse_expr(var2_vector[vr2])) %>%
        emmeans_test(Shannon16s_new ~ Cq_presence_absence, p.adjust.method = "BH")
      x <- data.frame(var1 = var1_vector[vr1],
                      var2 = var2_vector[vr2],
                      s1 = pwc$statistic[1],
                      p1 = pwc$p.adj[1],
                      s2 = pwc$statistic[2],
                      p2 = pwc$p.adj[2],
                      s3 = pwc$statistic[3],
                      p3 = pwc$p.adj[3],
                      s4 = pwc$statistic[4],
                      p4 = pwc$p.adj[4])
      LIST_3WAY_ANOVA[[length(LIST_3WAY_ANOVA) +1]] <- x
    }
  }
}

#Generate table that summarize all alpha and beta statistic
summary_stat <- data.frame()
for (i in 1:length(LIST_3WAY_ANOVA)){
  summary_stat <- rbind(summary_stat,LIST_3WAY_ANOVA[[i]])
}

#Create column that give combinason of significant between 2 level
summary_stat <- summary_stat %>% mutate(Level1_situation = case_when(s1 < 0 & p1 <0.05 ~ "Neg1Sig",
                                                                     s1 < 0 & p1 >0.05 ~ "Neg1NS",
                                                                     s1 > 0 & p1 <0.05 ~ "Pos1Sig",
                                                                     s1 > 0 & p1 >0.05 ~ "Pos1NS")
)

summary_stat <- summary_stat %>% mutate(Level2_situation = case_when(s2 < 0 & p2 <0.05 ~ "Neg2Sig",
                                                                     s2 < 0 & p2 >0.05 ~ "Neg2NS",
                                                                     s2 > 0 & p2 <0.05 ~ "Pos2Sig",
                                                                     s2 > 0 & p2 >0.05 ~ "Pos2NS")
)

summary_stat <- summary_stat %>% mutate(Level3_situation = case_when(s3 < 0 & p3 <0.05 ~ "Neg3Sig",
                                                                     s3 < 0 & p3 >0.05 ~ "Neg3NS",
                                                                     s3 > 0 & p3 <0.05 ~ "Pos3Sig",
                                                                     s3 > 0 & p3 >0.05 ~ "Pos3NS")
)

summary_stat <- summary_stat %>% mutate(Level4_situation = case_when(s4 < 0 & p4 <0.05 ~ "Neg4Sig",
                                                                     s4 < 0 & p4 >0.05 ~ "Neg4NS",
                                                                     s4 > 0 & p4 <0.05 ~ "Pos4Sig",
                                                                     s4 > 0 & p4 >0.05 ~ "Pos4NS")
)

summary_stat$Significance <- paste0(summary_stat$Level1_situation,summary_stat$Level2_situation,summary_stat$Level3_situation,summary_stat$Level4_situation)

ggplot(summary_stat, aes(x=var2,y=var1, fill=Significance)) +
  geom_tile()

pwc <- qPCR_Metadata %>%
  group_by(sexe,eau_fontaine) %>%
  emmeans_test(Shannon16s_new ~ Cq_presence_absence, p.adjust.method = "BH")

bxp <- ggboxplot(
  qPCR_Metadata, x = "eau_fontaine", y = "Shannon16s_new", 
  color = "Cq_presence_absence", palette = c("blue","orange"), facet.by = "sexe",add = "jitter",legend = "R")
bxp


##################################################################

##################################################################
### 3 WAYS ANOVA Shannon_16s   ###
##################################################################


qPCR_Metadata %>% sample_n_by(pays, Cq_presence_absence, stunted, size = 1)
qPCR_Metadata %>%
  group_by(pays, Cq_presence_absence, stunted) %>%
  get_summary_stats(Shannon16s_new, type = "mean_sd")


bxp <- ggboxplot(
  qPCR_Metadata, x = "stunted", y = "Shannon16s_new", 
  color = "Cq_presence_absence", palette = c("blue","orange"), facet.by = "pays",add = "jitter",legend = "R")
bxp


qPCR_Metadata %>%
  group_by(pays, Cq_presence_absence, stunted) %>%
  identify_outliers(Shannon16s_new)


model  <- lm(Shannon16s_new ~ pays*Cq_presence_absence*stunted, data = qPCR_Metadata)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))



qPCR_Metadata %>%
  group_by(pays, Cq_presence_absence, stunted) %>%
  shapiro_test(Shannon16s_new)
qPCR_Metadata %>%
  group_by(pays, Cq_presence_absence) %>%
  shapiro_test(Shannon16s_new)



ggqqplot(qPCR_Metadata, "Shannon16s_new", ggtheme = theme_bw()) +
  facet_grid(pays + Cq_presence_absence ~ stunted, labeller = "label_both")


qPCR_Metadata %>% levene_test(Shannon16s_new ~ pays*Cq_presence_absence*stunted)


res.aov <- qPCR_Metadata %>% anova_test(Shannon16s_new ~ pays*Cq_presence_absence*stunted)
res.aov


# Group the data by pays and 
# fit simple two-way interaction 
model  <- lm(Shannon16s_new ~ pays*Cq_presence_absence*stunted, data = qPCR_Metadata)
qPCR_Metadata %>%
  group_by(pays) %>%
  anova_test(Shannon16s_new ~ Cq_presence_absence*stunted, error = model)

# Group the data by pays and Cq_presence_absence, and fit  anova
stunted.effect <- qPCR_Metadata %>%
  group_by(pays, stunted) %>%
  anova_test(Shannon16s_new ~ Cq_presence_absence, error = model)
stunted.effect

# Pairwise comparisons
library(emmeans)
pwc <- qPCR_Metadata %>%
  group_by(pays,stunted) %>%
  emmeans_test(Shannon16s_new ~ Cq_presence_absence, p.adjust.method = "BH")
# Show comparison results for male at high Cq_presence_absence
pwc %>% filter(stunted == "Stunted")

# Estimated marginal means (i.e. adjusted means) 
# with 95% confidence interval
get_emmeans(pwc) %>% filter(stunted == "Stunted")


# Visualization: box plots with p-values
pwc <- pwc %>% add_xy_position(x = "stunted")
pwc.filtered <- pwc
ANOVA_3WAY_Shannon_Blastocystis_presence <- bxp +
  stat_pvalue_manual(pwc.filtered,
                     tip.length = 0, 
                     step.increase = 0.1, 
                     step.group.by = "pays") +
  labs(subtitle = get_test_label(res.aov, detailed = TRUE),
       caption = get_pwc_label(pwc))
ANOVA_3WAY_Shannon_Blastocystis_presence 




##################################################################

##################################################################
### 3 WAYS ANOVA Number_ASVs_16s   ###
##################################################################


qPCR_Metadata %>% sample_n_by(pays, Cq_presence_absence, stunted, size = 1)


qPCR_Metadata %>%
  group_by(pays, Cq_presence_absence, stunted) %>%
  get_summary_stats(Number_ASVs16s_new, type = "mean_sd")


bxp <- ggboxplot(
  qPCR_Metadata, x = "stunted", y = "Number_ASVs16s_new", 
  color = "Cq_presence_absence", palette = c("blue","orange"), facet.by = "pays",add = "jitter",legend = "R")
bxp


qPCR_Metadata %>%
  group_by(pays, Cq_presence_absence, stunted) %>%
  identify_outliers(Number_ASVs16s_new)


model  <- lm(Number_ASVs16s_new ~ pays*Cq_presence_absence*stunted, data = qPCR_Metadata)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))



qPCR_Metadata %>%
  group_by(pays, Cq_presence_absence, stunted) %>%
  shapiro_test(Number_ASVs16s_new)



ggqqplot(qPCR_Metadata, "Number_ASVs16s_new", ggtheme = theme_bw()) +
  facet_grid(pays + Cq_presence_absence ~ stunted, labeller = "label_both")


qPCR_Metadata %>% levene_test(Number_ASVs16s_new ~ pays*Cq_presence_absence*stunted)


res.aov <- qPCR_Metadata %>% anova_test(Number_ASVs16s_new ~ pays*Cq_presence_absence*stunted)
res.aov


# Group the data by pays and 
# fit simple two-way interaction 
model  <- lm(Number_ASVs16s_new ~ pays*Cq_presence_absence*stunted, data = qPCR_Metadata)
qPCR_Metadata %>%
  group_by(pays) %>%
  anova_test(Number_ASVs16s_new ~ Cq_presence_absence*stunted, error = model)

# Group the data by pays and Cq_presence_absence, and fit  anova
stunted.effect <- qPCR_Metadata %>%
  group_by(pays, stunted) %>%
  anova_test(Number_ASVs16s_new ~ Cq_presence_absence, error = model)
#stunted.effect %>% filter(pays == "Madagascar")

# Pairwise comparisons
library(emmeans)
pwc <- qPCR_Metadata %>%
  group_by(pays,stunted) %>%
  emmeans_test(Number_ASVs16s_new ~ Cq_presence_absence, p.adjust.method = "BH")
# Show comparison results for male at high Cq_presence_absence
pwc %>% filter(stunted == "Stunted")

# Estimated marginal means (i.e. adjusted means) 
# with 95% confidence interval
get_emmeans(pwc) %>% filter(stunted == "Stunted")


# Visualization: box plots with p-values
pwc <- pwc %>% add_xy_position(x = "stunted")
pwc.filtered <- pwc
ANOVA_3WAY_NumberASV_Blastocystis_presence <- bxp +
  stat_pvalue_manual(
    pwc.filtered, hide.ns = TRUE,
    tip.length = 0, step.increase = 0.1, step.group.by = "pays") +
  labs(
    subtitle = get_test_label(res.aov, detailed = TRUE),
    caption = get_pwc_label(pwc))


plot_combined <- ggarrange(ANOVA_3WAY_Shannon_Blastocystis_presence,ANOVA_3WAY_NumberASV_Blastocystis_presence,labels = c("A","B"), ncol = 2, nrow=1);plot_combined
ggsave2(plot_combined, filename = "/Users/vincebilly/Desktop/vince/phd/statistic/afribiota_vince/Simple_risk_factor/figure/3WAY_ANOVA_Shannon_NumberASV_Blastocystis.pdf", height = 6, width = 12)



##################################################################

##################################################################
### 3 WAYS ANOVA Shannon_16s with cheese  ###
##################################################################


qPCR_Metadata %>% sample_n_by(pays, Cq_presence_absence, fromage, size = 1)
qPCR_Metadata %>%
  group_by(pays, Cq_presence_absence, fromage) %>%
  get_summary_stats(Shannon16s_new, type = "mean_sd")


bxp <- ggboxplot(
  qPCR_Metadata, x = "fromage", y = "Shannon16s_new", 
  color = "Cq_presence_absence", palette = c("blue","orange"), facet.by = "pays",add = "jitter",legend = "R")
bxp


qPCR_Metadata %>%
  group_by(pays, Cq_presence_absence, fromage) %>%
  identify_outliers(Shannon16s_new)


model  <- lm(Shannon16s_new ~ pays*Cq_presence_absence*fromage, data = qPCR_Metadata)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))



qPCR_Metadata %>%
  group_by(pays, Cq_presence_absence, fromage) %>%
  shapiro_test(Shannon16s_new)



ggqqplot(qPCR_Metadata, "Shannon16s_new", ggtheme = theme_bw()) +
  facet_grid(pays + Cq_presence_absence ~ fromage, labeller = "label_both")


qPCR_Metadata %>% levene_test(Shannon16s_new ~ pays*Cq_presence_absence*fromage)


res.aov <- qPCR_Metadata %>% anova_test(Shannon16s_new ~ pays*Cq_presence_absence*fromage)
res.aov


# Group the data by pays and 
# fit simple two-way interaction 
model  <- lm(Shannon16s_new ~ pays*Cq_presence_absence*fromage, data = qPCR_Metadata)
qPCR_Metadata %>%
  group_by(pays) %>%
  anova_test(Shannon16s_new ~ Cq_presence_absence*fromage, error = model)

# Group the data by pays and Cq_presence_absence, and fit  anova
stunted.effect <- qPCR_Metadata %>%
  group_by(pays, fromage) %>%
  anova_test(Shannon16s_new ~ Cq_presence_absence, error = model)
stunted.effect

# Pairwise comparisons
library(emmeans)
pwc <- qPCR_Metadata %>%
  group_by(pays,fromage) %>%
  emmeans_test(Shannon16s_new ~ Cq_presence_absence, p.adjust.method = "BH")
# Show comparison results for male at high Cq_presence_absence
pwc %>% filter(fromage == "Non")

# Estimated marginal means (i.e. adjusted means) 
# with 95% confidence interval
get_emmeans(pwc) %>% filter(fromage == "Non")


# Visualization: box plots with p-values
pwc <- pwc %>% add_xy_position(x = "fromage")
pwc.filtered <- pwc %>% filter(fromage == "Non")
pwc.filtered$y.position <- c(5,5,5,5)
ANOVA_3WAY_Shannon_Blastocystis_presence <- bxp +
  stat_pvalue_manual(pwc.filtered, 
                     hide.ns = TRUE,
                     tip.length = 0, 
                     step.increase = 0.1, 
                     step.group.by = "pays") +
  labs(subtitle = get_test_label(res.aov, detailed = TRUE),
       caption = get_pwc_label(pwc))
ANOVA_3WAY_Shannon_Blastocystis_presence 




##################################################################

##################################################################
### 3 WAYS ANOVA Shannon_euk   ###
##################################################################


qPCR_Metadata %>% sample_n_by(pays, Cq_presence_absence, stunted, size = 1)


qPCR_Metadata %>%
  group_by(pays, Cq_presence_absence, stunted) %>%
  get_summary_stats(Shannon_euk, type = "mean_sd")


bxp <- ggboxplot(
  qPCR_Metadata, x = "stunted", y = "Shannon_euk", 
  color = "Cq_presence_absence", palette = "jco", facet.by = "pays"
)
bxp


qPCR_Metadata %>%
  group_by(pays, Cq_presence_absence, stunted) %>%
  identify_outliers(Shannon_euk)


model  <- lm(Shannon_euk ~ pays*Cq_presence_absence*stunted, data = qPCR_Metadata)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))



qPCR_Metadata %>%
  group_by(pays, Cq_presence_absence, stunted) %>%
  shapiro_test(Shannon_euk)



ggqqplot(qPCR_Metadata, "Shannon_euk", ggtheme = theme_bw()) +
  facet_grid(pays + Cq_presence_absence ~ stunted, labeller = "label_both")


qPCR_Metadata %>% levene_test(Shannon_euk ~ pays*Cq_presence_absence*stunted)


res.aov <- qPCR_Metadata %>% anova_test(Shannon_euk ~ pays*Cq_presence_absence*stunted)
res.aov


# Group the data by pays and 
# fit simple two-way interaction 
model  <- lm(Shannon_euk ~ pays*Cq_presence_absence*stunted, data = qPCR_Metadata)
qPCR_Metadata %>%
  group_by(pays) %>%
  anova_test(Shannon_euk ~ Cq_presence_absence*stunted, error = model)

# Group the data by pays and Cq_presence_absence, and fit  anova
stunted.effect <- qPCR_Metadata %>%
  group_by(pays, stunted) %>%
  anova_test(Shannon_euk ~ Cq_presence_absence, error = model)
#stunted.effect %>% filter(pays == "Madagascar")

stunted.effect <- qPCR_Metadata %>%
  group_by(pays, stunted ) %>%
  anova_test(Shannon_euk ~ Cq_presence_absence, error = model)

# Pairwise comparisons
library(emmeans)
pwc <- qPCR_Metadata %>%
  group_by(pays,stunted) %>%
  emmeans_test(Shannon_euk ~ Cq_presence_absence, p.adjust.method = "bonferroni") %>%
  select(-df, -statistic, -p) # Remove details
# Show comparison results for male at high Cq_presence_absence
pwc %>% filter(stunted == "Oui")

# Estimated marginal means (i.e. adjusted means) 
# with 95% confidence interval
get_emmeans(pwc) %>% filter(stunted == "Oui")


# Visualization: box plots with p-values
pwc <- pwc %>% add_xy_position(x = "stunted")
pwc.filtered <- pwc %>% filter(stunted == "Oui")
bxp +
  stat_pvalue_manual(
    pwc.filtered, color = "Cq_presence_absence", linetype = "Cq_presence_absence", hide.ns = TRUE,
    tip.length = 0, step.increase = 0.1, step.group.by = "pays"
  ) +
  labs(
    subtitle = get_test_label(res.aov, detailed = TRUE),
    caption = get_pwc_label(pwc)
  )

bxp +
  labs(
    subtitle = get_test_label(res.aov, detailed = TRUE),
    caption = get_pwc_label(pwc)
  )




##################################################################

##################################################################
### 3 WAYS ANOVA Shannon ~stunted + pays+ Blasto  ###
##################################################################




qPCR_Metadata %>%
  group_by(pays, Cq_presence_absence, stunted) %>%
  get_summary_stats(Shannon, type = "mean_sd")


bxp <- ggboxplot(
  qPCR_Metadata, x = "Cq_presence_absence", y = "Shannon", 
  color = "stunted", palette = "jco", facet.by = "pays"
)
bxp


qPCR_Metadata %>%
  group_by(pays, Cq_presence_absence, stunted) %>%
  identify_outliers(Shannon)


model  <- lm(Shannon ~ pays*Cq_presence_absence*stunted, data = qPCR_Metadata);summary(model)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))



qPCR_Metadata %>%
  group_by(pays, Cq_presence_absence, stunted) %>%
  shapiro_test(Shannon)



ggqqplot(qPCR_Metadata, "Shannon", ggtheme = theme_bw()) +
  facet_grid(pays + Cq_presence_absence ~ stunted, labeller = "label_both")


qPCR_Metadata %>% levene_test(Shannon ~ pays*Cq_presence_absence*stunted)


res.aov <- qPCR_Metadata %>% anova_test(Shannon ~ pays*Cq_presence_absence*stunted)
res.aov


# Group the data by pays and 
# fit simple two-way interaction 
model  <- lm(Shannon ~ pays*Cq_presence_absence*stunted, data = qPCR_Metadata)
qPCR_Metadata %>%
  group_by(pays) %>%
  anova_test(Shannon ~ Cq_presence_absence*stunted, error = model)

# Group the data by pays and Cq_presence_absence, and fit  anova
stunted.effect <- qPCR_Metadata %>%
  group_by(pays, Cq_presence_absence) %>%
  anova_test(Shannon ~ stunted, error = model)
#stunted.effect %>% filter(pays == "Madagascar")

stunted.effect <- qPCR_Metadata %>%
  group_by(pays, Cq_presence_absence  ) %>%
  anova_test(Shannon ~ stunted, error = model)

# Pairwise comparisons
library(emmeans)
pwc <- qPCR_Metadata %>%
  group_by(pays, Cq_presence_absence) %>%
  emmeans_test(Shannon ~ stunted, p.adjust.method = "bonferroni") %>%
  select(-df, -statistic, -p) # Remove details
# Show comparison results for male at high Cq_presence_absence
pwc %>% filter(stunted == "Oui")

# Estimated marginal means (i.e. adjusted means) 
# with 95% confidence interval
get_emmeans(pwc) %>% filter(stunted == "Oui")


# Visualization: box plots with p-values
pwc <- pwc %>% add_xy_position(x = "stunted")
pwc.filtered <- pwc %>% filter(stunted == "Oui")
bxp +
  stat_pvalue_manual(
    pwc.filtered, color = "Cq_presence_absence", linetype = "Cq_presence_absence", hide.ns = TRUE,
    tip.length = 0, step.increase = 0.1, step.group.by = "pays"
  ) +
  labs(
    subtitle = get_test_label(res.aov, detailed = TRUE),
    caption = get_pwc_label(pwc)
  )

bxp +
  labs(
    subtitle = get_test_label(res.aov, detailed = TRUE),
    caption = get_pwc_label(pwc)
  )




##################################################################

##################################################################
### 3 WAYS ANOVA Number_ASVs   ###
##################################################################


qPCR_Metadata %>% sample_n_by(pays, Cq_presence_absence, stunted, size = 1)


qPCR_Metadata %>%
  group_by(pays, Cq_presence_absence, stunted) %>%
  get_summary_stats(Number_ASVs_euk, type = "mean_sd")


bxp <- ggboxplot(
  qPCR_Metadata, x = "stunted", y = "Number_ASVs_euk", 
  color = "Cq_presence_absence", palette = "jco", facet.by = "pays"
)
bxp


qPCR_Metadata %>%
  group_by(pays, Cq_presence_absence, stunted) %>%
  identify_outliers(Number_ASVs_euk)


model  <- lm(Number_ASVs_euk ~ pays*Cq_presence_absence*stunted, data = qPCR_Metadata)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))



qPCR_Metadata %>%
  group_by(pays, Cq_presence_absence, stunted) %>%
  shapiro_test(Number_ASVs_euk)



ggqqplot(qPCR_Metadata, "Number_ASVs_euk", ggtheme = theme_bw()) +
  facet_grid(pays + Cq_presence_absence ~ stunted, labeller = "label_both")


qPCR_Metadata %>% levene_test(Number_ASVs_euk ~ pays*Cq_presence_absence*stunted)


res.aov <- qPCR_Metadata %>% anova_test(Number_ASVs_euk ~ pays*Cq_presence_absence*stunted)
res.aov


# Group the data by pays and 
# fit simple two-way interaction 
model  <- lm(Number_ASVs_euk ~ pays*Cq_presence_absence*stunted, data = qPCR_Metadata)
qPCR_Metadata %>%
  group_by(pays) %>%
  anova_test(Number_ASVs_euk ~ Cq_presence_absence*stunted, error = model)

# Group the data by pays and Cq_presence_absence, and fit  anova
stunted.effect <- qPCR_Metadata %>%
  group_by(pays, stunted) %>%
  anova_test(Number_ASVs_euk ~ Cq_presence_absence, error = model)
#stunted.effect %>% filter(pays == "Madagascar")

stunted.effect <- qPCR_Metadata %>%
  group_by(pays, stunted ) %>%
  anova_test(Number_ASVs_euk ~ Cq_presence_absence, error = model)

# Pairwise comparisons
library(emmeans)
pwc <- qPCR_Metadata %>%
  group_by(pays,stunted) %>%
  emmeans_test(Number_ASVs ~ Cq_presence_absence, p.adjust.method = "bonferroni") %>%
  select(-df, -statistic, -p) # Remove details
# Show comparison results for male at high Cq_presence_absence
pwc %>% filter(stunted == "Oui")

# Estimated marginal means (i.e. adjusted means) 
# with 95% confidence interval
get_emmeans(pwc) %>% filter(stunted == "Oui")


# Visualization: box plots with p-values
pwc <- pwc %>% add_xy_position(x = "stunted")
pwc.filtered <- pwc %>% filter(stunted == "Oui")
bxp +
  stat_pvalue_manual(
    pwc.filtered, color = "Cq_presence_absence", linetype = "Cq_presence_absence", hide.ns = TRUE,
    tip.length = 0, step.increase = 0.1, step.group.by = "pays"
  ) +
  labs(
    subtitle = get_test_label(res.aov, detailed = TRUE),
    caption = get_pwc_label(pwc)
  )

bxp +
  labs(
    subtitle = get_test_label(res.aov, detailed = TRUE),
    caption = get_pwc_label(pwc)
  )




##################################################################

##################################################################
### 3 WAYS ANOVA Calprotectine   ###
##################################################################

qPCR_Metadata$CALPROTECTINEggdePS

qPCR_Metadata %>% sample_n_by(pays, Cq_presence_absence, stunted, size = 1)
qPCR_Metadata %>%
  group_by(pays, Cq_presence_absence, stunted) %>%
  get_summary_stats(CALPROTECTINEggdePS, type = "mean_sd")


bxp <- ggboxplot(
  qPCR_Metadata, x = "stunted", y = "CALPROTECTINEggdePS", 
  color = "Cq_presence_absence", palette = c("blue","orange"), facet.by = "pays",add = "jitter",legend = "R")
bxp


qPCR_Metadata %>%
  group_by(pays, Cq_presence_absence, stunted) %>%
  identify_outliers(CALPROTECTINEggdePS)


model  <- lm(CALPROTECTINEggdePS ~ pays*Cq_presence_absence*stunted, data = qPCR_Metadata)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))



qPCR_Metadata %>%
  group_by(pays, Cq_presence_absence, stunted) %>%
  shapiro_test(CALPROTECTINEggdePS)



ggqqplot(qPCR_Metadata, "CALPROTECTINEggdePS", ggtheme = theme_bw()) +
  facet_grid(pays + Cq_presence_absence ~ stunted, labeller = "label_both")


qPCR_Metadata %>% levene_test(CALPROTECTINEggdePS ~ pays*Cq_presence_absence*stunted)


res.aov <- qPCR_Metadata %>% anova_test(CALPROTECTINEggdePS ~ pays*Cq_presence_absence*stunted)
res.aov


# Group the data by pays and 
# fit simple two-way interaction 
model  <- lm(CALPROTECTINEggdePS ~ pays*Cq_presence_absence*stunted, data = qPCR_Metadata)
qPCR_Metadata %>%
  group_by(pays) %>%
  anova_test(CALPROTECTINEggdePS ~ Cq_presence_absence*stunted, error = model)

# Group the data by pays and Cq_presence_absence, and fit  anova
stunted.effect <- qPCR_Metadata %>%
  group_by(pays, stunted) %>%
  anova_test(CALPROTECTINEggdePS ~ Cq_presence_absence, error = model)
stunted.effect

# Pairwise comparisons
library(emmeans)
pwc <- qPCR_Metadata %>%
  group_by(pays,stunted) %>%
  emmeans_test(CALPROTECTINEggdePS ~ Cq_presence_absence, p.adjust.method = "BH")
# Show comparison results for male at high Cq_presence_absence
pwc %>% filter(stunted == "Stunted")

# Estimated marginal means (i.e. adjusted means) 
# with 95% confidence interval
get_emmeans(pwc) %>% filter(stunted == "Stunted")


# Visualization: box plots with p-values
pwc <- pwc %>% add_xy_position(x = "stunted")
pwc.filtered <- pwc #%>% filter(stunted == "Stunted")
ANOVA_3WAY_Calprotectine_Blastocystis_presence <- bxp +
  stat_pvalue_manual(pwc.filtered,
                     tip.length = 0, 
                     step.increase = 0.1, 
                     step.group.by = "pays") +
  labs(subtitle = get_test_label(res.aov, detailed = TRUE),
       caption = get_pwc_label(pwc))
ANOVA_3WAY_Calprotectine_Blastocystis_presence 




##################################################################

##################################################################
### 3 WAYS ANOVA AATmggdePS   ###
##################################################################

qPCR_Metadata %>% sample_n_by(pays, Cq_presence_absence, stunted, size = 1)
qPCR_Metadata %>%
  group_by(pays, Cq_presence_absence, stunted) %>%
  get_summary_stats(AATmggdePS, type = "mean_sd")


bxp <- ggboxplot(
  qPCR_Metadata, x = "stunted", y = "AATmggdePS", 
  color = "Cq_presence_absence", palette = c("blue","orange"), facet.by = "pays",add = "jitter",legend = "R")
bxp


qPCR_Metadata %>%
  group_by(pays, Cq_presence_absence, stunted) %>%
  identify_outliers(AATmggdePS)


model  <- lm(AATmggdePS ~ pays*Cq_presence_absence*stunted, data = qPCR_Metadata)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))



qPCR_Metadata %>%
  group_by(pays, Cq_presence_absence, stunted) %>%
  shapiro_test(AATmggdePS)



ggqqplot(qPCR_Metadata, "AATmggdePS", ggtheme = theme_bw()) +
  facet_grid(pays + Cq_presence_absence ~ stunted, labeller = "label_both")


qPCR_Metadata %>% levene_test(AATmggdePS ~ pays*Cq_presence_absence*stunted)


res.aov <- qPCR_Metadata %>% anova_test(AATmggdePS ~ pays*Cq_presence_absence*stunted)
res.aov


# Group the data by pays and 
# fit simple two-way interaction 
model  <- lm(AATmggdePS ~ pays*Cq_presence_absence*stunted, data = qPCR_Metadata)
qPCR_Metadata %>%
  group_by(pays) %>%
  anova_test(AATmggdePS ~ Cq_presence_absence*stunted, error = model)

# Group the data by pays and Cq_presence_absence, and fit  anova
stunted.effect <- qPCR_Metadata %>%
  group_by(pays, stunted) %>%
  anova_test(AATmggdePS ~ Cq_presence_absence, error = model)
stunted.effect

# Pairwise comparisons
library(emmeans)
pwc <- qPCR_Metadata %>%
  group_by(pays,stunted) %>%
  emmeans_test(AATmggdePS ~ Cq_presence_absence, p.adjust.method = "BH")
# Show comparison results for male at high Cq_presence_absence
pwc %>% filter(stunted == "Stunted")

# Estimated marginal means (i.e. adjusted means) 
# with 95% confidence interval
get_emmeans(pwc) %>% filter(pays == "Madagascar")


# Visualization: box plots with p-values
pwc <- pwc %>% add_xy_position(x = "stunted")
pwc.filtered <- pwc #%>% filter(stunted == "Stunted")
ANOVA_3WAY_AAT_Blastocystis_presence <- bxp +
  stat_pvalue_manual(pwc.filtered,
                     tip.length = 0, 
                     step.increase = 0.1, 
                     step.group.by = "pays") +
  labs(subtitle = get_test_label(res.aov, detailed = TRUE),
       caption = get_pwc_label(pwc))
ANOVA_3WAY_AAT_Blastocystis_presence 




##################################################################

##################################################################
### 3 WAYS ANOVA haz_cont   ###
##################################################################


qPCR_Metadata %>% sample_n_by(pays, Cq_presence_absence, stunted, size = 1)
qPCR_Metadata %>%
  group_by(pays, Cq_presence_absence, stunted) %>%
  get_summary_stats(haz_cont, type = "mean_sd")


bxp <- ggboxplot(
  qPCR_Metadata, x = "stunted", y = "haz_cont", 
  color = "Cq_presence_absence", palette = c("blue","orange"), facet.by = "pays",add = "jitter",legend = "R")
bxp


qPCR_Metadata %>%
  group_by(pays, Cq_presence_absence, stunted) %>%
  identify_outliers(haz_cont)


model  <- lm(haz_cont ~ pays*Cq_presence_absence*stunted, data = qPCR_Metadata)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))



qPCR_Metadata %>%
  group_by(pays, Cq_presence_absence, stunted) %>%
  shapiro_test(haz_cont)



ggqqplot(qPCR_Metadata, "haz_cont", ggtheme = theme_bw()) +
  facet_grid(pays + Cq_presence_absence ~ stunted, labeller = "label_both")


qPCR_Metadata %>% levene_test(haz_cont ~ pays*Cq_presence_absence*stunted)


res.aov <- qPCR_Metadata %>% anova_test(haz_cont ~ pays*Cq_presence_absence*stunted)
res.aov


# Group the data by pays and 
# fit simple two-way interaction 
model  <- lm(haz_cont ~ pays*Cq_presence_absence*stunted, data = qPCR_Metadata)
qPCR_Metadata %>%
  group_by(pays) %>%
  anova_test(haz_cont ~ Cq_presence_absence*stunted, error = model)

# Group the data by pays and Cq_presence_absence, and fit  anova
stunted.effect <- qPCR_Metadata %>%
  group_by(pays, stunted) %>%
  anova_test(haz_cont ~ Cq_presence_absence, error = model)
stunted.effect

# Pairwise comparisons
library(emmeans)
pwc <- qPCR_Metadata %>%
  group_by(pays,stunted) %>%
  emmeans_test(haz_cont ~ Cq_presence_absence, p.adjust.method = "BH")
# Show comparison results for male at high Cq_presence_absence
pwc %>% filter(stunted == "Stunted")

# Estimated marginal means (i.e. adjusted means) 
# with 95% confidence interval
get_emmeans(pwc) %>% filter(pays == "Madagascar")


# Visualization: box plots with p-values
pwc <- pwc %>% add_xy_position(x = "stunted")
pwc.filtered <- pwc #%>% filter(stunted == "Stunted")
pwc.filtered$y.position <- c(1,1,1,1)
ANOVA_3WAY_haz_Blastocystis_presence <- bxp +
  stat_pvalue_manual(pwc.filtered,
                     tip.length = 0, 
                     step.increase = 0.1, 
                     step.group.by = "pays") +
  labs(subtitle = get_test_label(res.aov, detailed = TRUE),
       caption = get_pwc_label(pwc))
ANOVA_3WAY_haz_Blastocystis_presence 




##################################################################

##################################################################
### 3 WAYS ANOVA whz_cont   ###
##################################################################



qPCR_Metadata %>% sample_n_by(pays, Cq_presence_absence, stunted, size = 1)
qPCR_Metadata %>%
  group_by(pays, Cq_presence_absence, stunted) %>%
  get_summary_stats(whz_cont, type = "mean_sd")


bxp <- ggboxplot(
  qPCR_Metadata, x = "stunted", y = "whz_cont", 
  color = "Cq_presence_absence", palette = c("blue","orange"), facet.by = "pays",add = "jitter",legend = "R")
bxp


qPCR_Metadata %>%
  group_by(pays, Cq_presence_absence, stunted) %>%
  identify_outliers(whz_cont)


model  <- lm(whz_cont ~ pays*Cq_presence_absence*stunted, data = qPCR_Metadata)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))



qPCR_Metadata %>%
  group_by(pays, Cq_presence_absence, stunted) %>%
  shapiro_test(whz_cont)



ggqqplot(qPCR_Metadata, "whz_cont", ggtheme = theme_bw()) +
  facet_grid(pays + Cq_presence_absence ~ stunted, labeller = "label_both")


qPCR_Metadata %>% levene_test(whz_cont ~ pays*Cq_presence_absence*stunted)


res.aov <- qPCR_Metadata %>% anova_test(whz_cont ~ pays*Cq_presence_absence*stunted)
res.aov


# Group the data by pays and 
# fit simple two-way interaction 
model  <- lm(whz_cont ~ pays*Cq_presence_absence*stunted, data = qPCR_Metadata)
qPCR_Metadata %>%
  group_by(pays) %>%
  anova_test(whz_cont ~ Cq_presence_absence*stunted, error = model)

# Group the data by pays and Cq_presence_absence, and fit  anova
stunted.effect <- qPCR_Metadata %>%
  group_by(pays, stunted) %>%
  anova_test(whz_cont ~ Cq_presence_absence, error = model)
stunted.effect

# Pairwise comparisons
library(emmeans)
pwc <- qPCR_Metadata %>%
  group_by(pays,stunted) %>%
  emmeans_test(whz_cont ~ Cq_presence_absence, p.adjust.method = "BH")
# Show comparison results for male at high Cq_presence_absence
pwc %>% filter(stunted == "Stunted")

# Estimated marginal means (i.e. adjusted means) 
# with 95% confidence interval
get_emmeans(pwc) %>% filter(pays == "Madagascar")


# Visualization: box plots with p-values
pwc <- pwc %>% add_xy_position(x = "stunted")
pwc.filtered <- pwc #%>% filter(stunted == "Stunted")
ANOVA_3WAY_whz_Blastocystis_presence <- bxp +
  stat_pvalue_manual(pwc.filtered,
                     tip.length = 0, 
                     step.increase = 0.1, 
                     step.group.by = "pays") +
  labs(subtitle = get_test_label(res.aov, detailed = TRUE),
       caption = get_pwc_label(pwc))
ANOVA_3WAY_whz_Blastocystis_presence 




##################################################################

##################################################################
### 3 WAYS ANOVA hemoglobine   ###
##################################################################

qPCR_Metadata$hemoglobine <- as.numeric(qPCR_Metadata$hemoglobine)

qPCR_Metadata %>% sample_n_by(pays, Cq_presence_absence, stunted, size = 1)
qPCR_Metadata %>%
  group_by(pays, Cq_presence_absence, stunted) %>%
  get_summary_stats(hemoglobine, type = "mean_sd")


bxp <- ggboxplot(
  qPCR_Metadata, x = "stunted", y = "hemoglobine", 
  color = "Cq_presence_absence", palette = c("blue","orange"), facet.by = "pays",add = "jitter",legend = "R")
bxp


qPCR_Metadata %>%
  group_by(pays, Cq_presence_absence, stunted) %>%
  identify_outliers(hemoglobine)


model  <- lm(hemoglobine ~ pays*Cq_presence_absence*stunted, data = qPCR_Metadata)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))



qPCR_Metadata %>%
  group_by(pays, Cq_presence_absence, stunted) %>%
  shapiro_test(hemoglobine)



ggqqplot(qPCR_Metadata, "hemoglobine", ggtheme = theme_bw()) +
  facet_grid(pays + Cq_presence_absence ~ stunted, labeller = "label_both")


qPCR_Metadata %>% levene_test(hemoglobine ~ pays*Cq_presence_absence*stunted)


res.aov <- qPCR_Metadata %>% anova_test(hemoglobine ~ pays*Cq_presence_absence*stunted)
res.aov


# Group the data by pays and 
# fit simple two-way interaction 
model  <- lm(hemoglobine ~ pays*Cq_presence_absence*stunted, data = qPCR_Metadata)
qPCR_Metadata %>%
  group_by(pays) %>%
  anova_test(hemoglobine ~ Cq_presence_absence*stunted, error = model)

# Group the data by pays and Cq_presence_absence, and fit  anova
stunted.effect <- qPCR_Metadata %>%
  group_by(pays, stunted) %>%
  anova_test(hemoglobine ~ Cq_presence_absence, error = model)
stunted.effect

# Pairwise comparisons
library(emmeans)
pwc <- qPCR_Metadata %>%
  group_by(pays,stunted) %>%
  emmeans_test(hemoglobine ~ Cq_presence_absence, p.adjust.method = "BH")
# Show comparison results for male at high Cq_presence_absence
pwc %>% filter(stunted == "Stunted")

# Estimated marginal means (i.e. adjusted means) 
# with 95% confidence interval
get_emmeans(pwc) %>% filter(pays == "Madagascar")


# Visualization: box plots with p-values
pwc <- pwc %>% add_xy_position(x = "stunted")
pwc.filtered <- pwc #%>% filter(stunted == "Stunted")
ANOVA_3WAY_Hemoglobine_Blastocystis_presence <- bxp +
  stat_pvalue_manual(pwc.filtered,
                     tip.length = 0, 
                     step.increase = 0.1, 
                     step.group.by = "pays") +
  labs(subtitle = get_test_label(res.aov, detailed = TRUE),
       caption = get_pwc_label(pwc))
ANOVA_3WAY_Hemoglobine_Blastocystis_presence 




##################################################################

##################################################################
### 3 WAYS ANOVA ferritine_non_corr   ###
##################################################################

qPCR_Metadata$ferritine_non_corr <- as.numeric(qPCR_Metadata$ferritine_non_corr)

qPCR_Metadata %>% sample_n_by(pays, Cq_presence_absence, stunted, size = 1)
qPCR_Metadata %>%
  group_by(pays, Cq_presence_absence, stunted) %>%
  get_summary_stats(ferritine_non_corr, type = "mean_sd")


bxp <- ggboxplot(
  qPCR_Metadata, x = "stunted", y = "ferritine_non_corr", 
  color = "Cq_presence_absence", palette = c("blue","orange"), facet.by = "pays",add = "jitter",legend = "R")
bxp


qPCR_Metadata %>%
  group_by(pays, Cq_presence_absence, stunted) %>%
  identify_outliers(ferritine_non_corr)


model  <- lm(ferritine_non_corr ~ pays*Cq_presence_absence*stunted, data = qPCR_Metadata)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))



qPCR_Metadata %>%
  group_by(pays, Cq_presence_absence, stunted) %>%
  shapiro_test(ferritine_non_corr)



ggqqplot(qPCR_Metadata, "ferritine_non_corr", ggtheme = theme_bw()) +
  facet_grid(pays + Cq_presence_absence ~ stunted, labeller = "label_both")


qPCR_Metadata %>% levene_test(ferritine_non_corr ~ pays*Cq_presence_absence*stunted)


res.aov <- qPCR_Metadata %>% anova_test(ferritine_non_corr ~ pays*Cq_presence_absence*stunted)
res.aov


# Group the data by pays and 
# fit simple two-way interaction 
model  <- lm(ferritine_non_corr ~ pays*Cq_presence_absence*stunted, data = qPCR_Metadata)
qPCR_Metadata %>%
  group_by(pays) %>%
  anova_test(ferritine_non_corr ~ Cq_presence_absence*stunted, error = model)

# Group the data by pays and Cq_presence_absence, and fit  anova
stunted.effect <- qPCR_Metadata %>%
  group_by(pays, stunted) %>%
  anova_test(ferritine_non_corr ~ Cq_presence_absence, error = model)
stunted.effect

# Pairwise comparisons
library(emmeans)
pwc <- qPCR_Metadata %>%
  group_by(pays,stunted) %>%
  emmeans_test(ferritine_non_corr ~ Cq_presence_absence, p.adjust.method = "BH")
# Show comparison results for male at high Cq_presence_absence
pwc %>% filter(stunted == "Stunted")

# Estimated marginal means (i.e. adjusted means) 
# with 95% confidence interval
get_emmeans(pwc) %>% filter(pays == "RCA")


# Visualization: box plots with p-values
pwc <- pwc %>% add_xy_position(x = "stunted")
pwc.filtered <- pwc #%>% filter(stunted == "Stunted")
ANOVA_3WAY_ferritine_non_corr_Blastocystis_presence <- bxp +
  stat_pvalue_manual(pwc.filtered,
                     tip.length = 0, 
                     step.increase = 0.1, 
                     step.group.by = "pays") +
  labs(subtitle = get_test_label(res.aov, detailed = TRUE),
       caption = get_pwc_label(pwc))
ANOVA_3WAY_ferritine_non_corr_Blastocystis_presence 




##################################################################


plot_combined <- ggarrange(ANOVA_3WAY_AAT_Blastocystis_presence,
                           ANOVA_3WAY_ferritine_non_corr_Blastocystis_presence,
                           ANOVA_3WAY_Calprotectine_Blastocystis_presence,
                           ANOVA_3WAY_Hemoglobine_Blastocystis_presence,
                           ANOVA_3WAY_haz_Blastocystis_presence,
                           ANOVA_3WAY_whz_Blastocystis_presence,
                           labels = c("A","B"), ncol = 2, nrow=3);plot_combined
ggsave2(plot_combined, filename = "/Users/vincebilly/Desktop/vince/phd/statistic/afribiota_vince/Simple_risk_factor/figure/3WAY_ANOVA_clinical_Blastocystis.pdf", height = 17, width = 12)

##################################################################
### 3 WAYS ANOVA diet score - all   ###
##################################################################


qPCR_Metadata_veg <- qPCR_Metadata[,c("haricot","jus","pain","citrouille","tubercule","leg_vert","fr_leg_aut","cereale_racine_tbrc","legum_noix")]
for(i in colnames(qPCR_Metadata_veg)){
  qPCR_Metadata_veg[,i] <- ifelse(qPCR_Metadata_veg[,i] == "Oui",1,0)
  qPCR_Metadata_veg[,i] <- as.numeric(qPCR_Metadata_veg[,i])
}

qPCR_Metadata_meat <- qPCR_Metadata[,c("lait_hier","fromage","abat","viande","oeuf","poisson","pr_laitier","pr_carne","pr_carne_poisson","huile_viande")]
for(i in colnames(qPCR_Metadata_meat)){
  qPCR_Metadata_meat[,i] <- ifelse(qPCR_Metadata_meat[,i] == "Oui",-1,0)
  qPCR_Metadata_meat[,i] <- as.numeric(qPCR_Metadata_meat[,i])
}

qPCR_Metadata_diet <- merge(qPCR_Metadata_meat,qPCR_Metadata_veg, by="row.names")
qPCR_Metadata_diet$Row.names <- NULL

qPCR_Metadata$diet_score <- rowSums(qPCR_Metadata_diet)


qPCR_Metadata %>% sample_n_by(pays, Cq_presence_absence, stunted, size = 1)
qPCR_Metadata %>%
  group_by(pays, Cq_presence_absence, stunted) %>%
  get_summary_stats(diet_score, type = "mean_sd")


bxp <- ggboxplot(
  qPCR_Metadata, x = "stunted", y = "diet_score", 
  color = "Cq_presence_absence", palette = c("blue","orange"), facet.by = "pays",add = "jitter",legend = "R")
bxp


qPCR_Metadata %>%
  group_by(pays, Cq_presence_absence, stunted) %>%
  identify_outliers(diet_score)


model  <- lm(diet_score ~ pays*Cq_presence_absence*stunted, data = qPCR_Metadata)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))



qPCR_Metadata %>%
  group_by(pays, Cq_presence_absence, stunted) %>%
  shapiro_test(diet_score)



ggqqplot(qPCR_Metadata, "diet_score", ggtheme = theme_bw()) +
  facet_grid(pays + Cq_presence_absence ~ stunted, labeller = "label_both")


qPCR_Metadata %>% levene_test(diet_score ~ pays*Cq_presence_absence*stunted)


res.aov <- qPCR_Metadata %>% anova_test(diet_score ~ pays*Cq_presence_absence*stunted)
res.aov


# Group the data by pays and 
# fit simple two-way interaction 
model  <- lm(diet_score ~ pays*Cq_presence_absence*stunted, data = qPCR_Metadata)
qPCR_Metadata %>%
  group_by(pays) %>%
  anova_test(diet_score ~ Cq_presence_absence*stunted, error = model)

# Group the data by pays and Cq_presence_absence, and fit  anova
stunted.effect <- qPCR_Metadata %>%
  group_by(pays, stunted) %>%
  anova_test(diet_score ~ Cq_presence_absence, error = model)
stunted.effect

# Pairwise comparisons
library(emmeans)
pwc <- qPCR_Metadata %>%
  group_by(pays,stunted) %>%
  emmeans_test(diet_score ~ Cq_presence_absence, p.adjust.method = "BH")
# Show comparison results for male at high Cq_presence_absence
pwc %>% filter(stunted == "Stunted")

# Estimated marginal means (i.e. adjusted means) 
# with 95% confidence interval
get_emmeans(pwc) %>% filter(stunted == "Stunted")


# Visualization: box plots with p-values
pwc <- pwc %>% add_xy_position(x = "stunted")
pwc.filtered <- pwc #%>% filter(stunted == "Stunted")
pwc.filtered$y.position <- c(8,8,8,8)
ANOVA_3WAY_DietScore_Blastocystis_presence <- bxp +
  stat_pvalue_manual(pwc.filtered, 
                     tip.length = 0, 
                     step.increase = 0.1, 
                     step.group.by = "pays") +
  labs(subtitle = get_test_label(res.aov, detailed = TRUE),
       caption = get_pwc_label(pwc))
ANOVA_3WAY_DietScore_Blastocystis_presence 

##################################################################

##################################################################
### 3 WAYS ANOVA diet score - vegetarian   ###
##################################################################

qPCR_Metadata_veg <- qPCR_Metadata[,c("haricot","jus","pain","citrouille","tubercule","leg_vert","fr_leg_aut","cereale_racine_tbrc","legum_noix")]
for(i in colnames(qPCR_Metadata_veg)){
  qPCR_Metadata_veg[,i] <- ifelse(qPCR_Metadata_veg[,i] == "Oui",1,0)
  qPCR_Metadata_veg[,i] <- as.numeric(qPCR_Metadata_veg[,i])
}


qPCR_Metadata$diet_score_veg <- rowSums(qPCR_Metadata_veg)


qPCR_Metadata %>% sample_n_by(pays, Cq_presence_absence, stunted, size = 1)
qPCR_Metadata %>%
  group_by(pays, Cq_presence_absence, stunted) %>%
  get_summary_stats(diet_score_veg, type = "mean_sd")


bxp <- ggboxplot(
  qPCR_Metadata, x = "stunted", y = "diet_score_veg", 
  color = "Cq_presence_absence", palette = c("blue","orange"), facet.by = "pays",add = "jitter",legend = "R")
bxp


qPCR_Metadata %>%
  group_by(pays, Cq_presence_absence, stunted) %>%
  identify_outliers(diet_score_veg)


model  <- lm(diet_score_veg ~ pays*Cq_presence_absence*stunted, data = qPCR_Metadata)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))



qPCR_Metadata %>%
  group_by(pays, Cq_presence_absence, stunted) %>%
  shapiro_test(diet_score_veg)



ggqqplot(qPCR_Metadata, "diet_score_veg", ggtheme = theme_bw()) +
  facet_grid(pays + Cq_presence_absence ~ stunted, labeller = "label_both")


qPCR_Metadata %>% levene_test(diet_score_veg ~ pays*Cq_presence_absence*stunted)


res.aov <- qPCR_Metadata %>% anova_test(diet_score_veg ~ pays*Cq_presence_absence*stunted)
res.aov


# Group the data by pays and 
# fit simple two-way interaction 
model  <- lm(diet_score_veg ~ pays*Cq_presence_absence*stunted, data = qPCR_Metadata)
qPCR_Metadata %>%
  group_by(pays) %>%
  anova_test(diet_score_veg ~ Cq_presence_absence*stunted, error = model)

# Group the data by pays and Cq_presence_absence, and fit  anova
stunted.effect <- qPCR_Metadata %>%
  group_by(pays, stunted) %>%
  anova_test(diet_score_veg ~ Cq_presence_absence, error = model)
stunted.effect

# Pairwise comparisons
library(emmeans)
pwc <- qPCR_Metadata %>%
  group_by(pays,stunted) %>%
  emmeans_test(diet_score_veg ~ Cq_presence_absence, p.adjust.method = "BH")
# Show comparison results for male at high Cq_presence_absence
pwc %>% filter(stunted == "Stunted")

# Estimated marginal means (i.e. adjusted means) 
# with 95% confidence interval
get_emmeans(pwc) %>% filter(stunted == "Stunted")


# Visualization: box plots with p-values
pwc <- pwc %>% add_xy_position(x = "stunted")
pwc.filtered <- pwc #%>% filter(stunted == "Stunted")
pwc.filtered$y.position <- c(8,8,8,8)
ANOVA_3WAY_vegetarian_Blastocystis_presence <- bxp +
  stat_pvalue_manual(pwc.filtered, 
                     tip.length = 0, 
                     step.increase = 0.1, 
                     step.group.by = "pays") +
  labs(subtitle = get_test_label(res.aov, detailed = TRUE),
       caption = get_pwc_label(pwc))
ANOVA_3WAY_vegetarian_Blastocystis_presence 


##################################################################

##################################################################
### 3 WAYS ANOVA diet score - meat fat   ###
##################################################################

qPCR_Metadata_meat_fat <- qPCR_Metadata[,c("lait_hier","fromage","abat","viande","oeuf","poisson","pr_laitier","pr_carne","pr_carne_poisson","huile_viande")]
for(i in colnames(qPCR_Metadata_meat_fat)){
  qPCR_Metadata_meat_fat[,i] <- ifelse(qPCR_Metadata_meat_fat[,i] == "Oui",1,0)
  qPCR_Metadata_meat_fat[,i] <- as.numeric(qPCR_Metadata_meat_fat[,i])
}


qPCR_Metadata$diet_score_meat_fat <- rowSums(qPCR_Metadata_meat_fat)


qPCR_Metadata %>% sample_n_by(pays, Cq_presence_absence, stunted, size = 1)
qPCR_Metadata %>%
  group_by(pays, Cq_presence_absence, stunted) %>%
  get_summary_stats(diet_score_meat_fat, type = "mean_sd")


bxp <- ggboxplot(
  qPCR_Metadata, x = "stunted", y = "diet_score_meat_fat", 
  color = "Cq_presence_absence", palette = c("blue","orange"), facet.by = "pays",add = "jitter",legend = "R")
bxp


qPCR_Metadata %>%
  group_by(pays, Cq_presence_absence, stunted) %>%
  identify_outliers(diet_score_meat_fat)


model  <- lm(diet_score_meat_fat ~ pays*Cq_presence_absence*stunted, data = qPCR_Metadata)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))



qPCR_Metadata %>%
  group_by(pays, Cq_presence_absence, stunted) %>%
  shapiro_test(diet_score_meat_fat)



ggqqplot(qPCR_Metadata, "diet_score_meat_fat", ggtheme = theme_bw()) +
  facet_grid(pays + Cq_presence_absence ~ stunted, labeller = "label_both")


qPCR_Metadata %>% levene_test(diet_score_meat_fat ~ pays*Cq_presence_absence*stunted)


res.aov <- qPCR_Metadata %>% anova_test(diet_score_meat_fat ~ pays*Cq_presence_absence*stunted)
res.aov


# Group the data by pays and 
# fit simple two-way interaction 
model  <- lm(diet_score_meat_fat ~ pays*Cq_presence_absence*stunted, data = qPCR_Metadata)
qPCR_Metadata %>%
  group_by(pays) %>%
  anova_test(diet_score_meat_fat ~ Cq_presence_absence*stunted, error = model)

# Group the data by pays and Cq_presence_absence, and fit  anova
stunted.effect <- qPCR_Metadata %>%
  group_by(pays, stunted) %>%
  anova_test(diet_score_meat_fat ~ Cq_presence_absence, error = model)
stunted.effect

# Pairwise comparisons
library(emmeans)
pwc <- qPCR_Metadata %>%
  group_by(pays,stunted) %>%
  emmeans_test(diet_score_meat_fat ~ Cq_presence_absence, p.adjust.method = "BH")
# Show comparison results for male at high Cq_presence_absence
pwc %>% filter(stunted == "Stunted")

# Estimated marginal means (i.e. adjusted means) 
# with 95% confidence interval
get_emmeans(pwc) %>% filter(stunted == "Stunted")


# Visualization: box plots with p-values
pwc <- pwc %>% add_xy_position(x = "stunted")
pwc.filtered <- pwc #%>% filter(stunted == "Stunted")
pwc.filtered$y.position <- c(8,8,8,8)
ANOVA_3WAY_meat_fat_Blastocystis_presence <- bxp +
  stat_pvalue_manual(pwc.filtered, 
                     tip.length = 0, 
                     step.increase = 0.1, 
                     step.group.by = "pays") +
  labs(subtitle = get_test_label(res.aov, detailed = TRUE),
       caption = get_pwc_label(pwc))
ANOVA_3WAY_meat_fat_Blastocystis_presence 




##################################################################

##################################################################
### 3 WAYS ANOVA diet score - meat    ###
##################################################################

qPCR_Metadata_meat <- qPCR_Metadata[,c("abat","viande","oeuf","poisson","pr_carne","pr_carne_poisson")]
for(i in colnames(qPCR_Metadata_meat)){
  qPCR_Metadata_meat[,i] <- ifelse(qPCR_Metadata_meat[,i] == "Oui",1,0)
  qPCR_Metadata_meat[,i] <- as.numeric(qPCR_Metadata_meat[,i])
}


qPCR_Metadata$diet_score_meat <- rowSums(qPCR_Metadata_meat)


qPCR_Metadata %>% sample_n_by(pays, Cq_presence_absence, stunted, size = 1)
qPCR_Metadata %>%
  group_by(pays, Cq_presence_absence, stunted) %>%
  get_summary_stats(diet_score_meat, type = "mean_sd")


bxp <- ggboxplot(
  qPCR_Metadata, x = "stunted", y = "diet_score_meat", 
  color = "Cq_presence_absence", palette = c("blue","orange"), facet.by = "pays",add = "jitter",legend = "R")
bxp


qPCR_Metadata %>%
  group_by(pays, Cq_presence_absence, stunted) %>%
  identify_outliers(diet_score_meat)


model  <- lm(diet_score_meat ~ pays*Cq_presence_absence*stunted, data = qPCR_Metadata)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))



qPCR_Metadata %>%
  group_by(pays, Cq_presence_absence, stunted) %>%
  shapiro_test(diet_score_meat)



ggqqplot(qPCR_Metadata, "diet_score_meat", ggtheme = theme_bw()) +
  facet_grid(pays + Cq_presence_absence ~ stunted, labeller = "label_both")


qPCR_Metadata %>% levene_test(diet_score_meat ~ pays*Cq_presence_absence*stunted)


res.aov <- qPCR_Metadata %>% anova_test(diet_score_meat ~ pays*Cq_presence_absence*stunted)
res.aov


# Group the data by pays and 
# fit simple two-way interaction 
model  <- lm(diet_score_meat ~ pays*Cq_presence_absence*stunted, data = qPCR_Metadata)
qPCR_Metadata %>%
  group_by(pays) %>%
  anova_test(diet_score_meat ~ Cq_presence_absence*stunted, error = model)

# Group the data by pays and Cq_presence_absence, and fit  anova
stunted.effect <- qPCR_Metadata %>%
  group_by(pays, stunted) %>%
  anova_test(diet_score_meat ~ Cq_presence_absence, error = model)
stunted.effect

# Pairwise comparisons
library(emmeans)
pwc <- qPCR_Metadata %>%
  group_by(pays,stunted) %>%
  emmeans_test(diet_score_meat ~ Cq_presence_absence, p.adjust.method = "BH")
# Show comparison results for male at high Cq_presence_absence
pwc %>% filter(stunted == "Stunted")

# Estimated marginal means (i.e. adjusted means) 
# with 95% confidence interval
get_emmeans(pwc) %>% filter(stunted == "Stunted")


# Visualization: box plots with p-values
pwc <- pwc %>% add_xy_position(x = "stunted")
pwc.filtered <- pwc #%>% filter(stunted == "Stunted")
pwc.filtered$y.position <- c(8,8,8,8)
ANOVA_3WAY_meat_Blastocystis_presence <- bxp +
  stat_pvalue_manual(pwc.filtered, 
                     tip.length = 0, 
                     step.increase = 0.1, 
                     step.group.by = "pays") +
  labs(subtitle = get_test_label(res.aov, detailed = TRUE),
       caption = get_pwc_label(pwc))
ANOVA_3WAY_meat_Blastocystis_presence 

##################################################################

##################################################################
### 3 WAYS ANOVA diet score - cheese/milk    ###
##################################################################

qPCR_Metadata_milk <- qPCR_Metadata[,c("lait_hier","fromage","pr_laitier")]
for(i in colnames(qPCR_Metadata_milk)){
  qPCR_Metadata_milk[,i] <- ifelse(qPCR_Metadata_milk[,i] == "Oui",1,0)
  qPCR_Metadata_milk[,i] <- as.numeric(qPCR_Metadata_milk[,i])
}


qPCR_Metadata$diet_score_milk_cheese <- rowSums(qPCR_Metadata_milk)


qPCR_Metadata %>% sample_n_by(pays, Cq_presence_absence, stunted, size = 1)
qPCR_Metadata %>%
  group_by(pays, Cq_presence_absence, stunted) %>%
  get_summary_stats(diet_score_milk_cheese, type = "mean_sd")


bxp <- ggboxplot(
  qPCR_Metadata, x = "stunted", y = "diet_score_milk_cheese", 
  color = "Cq_presence_absence", palette = c("blue","orange"), facet.by = "pays",add = "jitter",legend = "R")
bxp


qPCR_Metadata %>%
  group_by(pays, Cq_presence_absence, stunted) %>%
  identify_outliers(diet_score_milk_cheese)


model  <- lm(diet_score_milk_cheese ~ pays*Cq_presence_absence*stunted, data = qPCR_Metadata)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))



qPCR_Metadata %>%
  group_by(pays, Cq_presence_absence, stunted) %>%
  shapiro_test(diet_score_milk_cheese)



ggqqplot(qPCR_Metadata, "diet_score_milk_cheese", ggtheme = theme_bw()) +
  facet_grid(pays + Cq_presence_absence ~ stunted, labeller = "label_both")


qPCR_Metadata %>% levene_test(diet_score_milk_cheese ~ pays*Cq_presence_absence*stunted)


res.aov <- qPCR_Metadata %>% anova_test(diet_score_milk_cheese ~ pays*Cq_presence_absence*stunted)
res.aov


# Group the data by pays and 
# fit simple two-way interaction 
model  <- lm(diet_score_milk_cheese ~ pays*Cq_presence_absence*stunted, data = qPCR_Metadata)
qPCR_Metadata %>%
  group_by(pays) %>%
  anova_test(diet_score_milk_cheese ~ Cq_presence_absence*stunted, error = model)

# Group the data by pays and Cq_presence_absence, and fit  anova
stunted.effect <- qPCR_Metadata %>%
  group_by(pays, stunted) %>%
  anova_test(diet_score_milk_cheese ~ Cq_presence_absence, error = model)
stunted.effect

# Pairwise comparisons
library(emmeans)
pwc <- qPCR_Metadata %>%
  group_by(pays,stunted) %>%
  emmeans_test(diet_score_milk_cheese ~ Cq_presence_absence, p.adjust.method = "BH")
# Show comparison results for male at high Cq_presence_absence
pwc %>% filter(stunted == "Stunted")

# Estimated marginal means (i.e. adjusted means) 
# with 95% confidence interval
get_emmeans(pwc) %>% filter(stunted == "Stunted")


# Visualization: box plots with p-values
pwc <- pwc %>% add_xy_position(x = "stunted")
pwc.filtered <- pwc #%>% filter(stunted == "Stunted")
pwc.filtered$y.position <- c(8,8,8,8)
ANOVA_3WAY_milk_Blastocystis_presence <- bxp +
  stat_pvalue_manual(pwc.filtered, 
                     tip.length = 0, 
                     step.increase = 0.1, 
                     step.group.by = "pays") +
  labs(subtitle = get_test_label(res.aov, detailed = TRUE),
       caption = get_pwc_label(pwc))
ANOVA_3WAY_milk_Blastocystis_presence 




##################################################################

##################################################################
### 3 WAYS ANOVA diet score - cheese/milk - oil    ###
##################################################################

qPCR_Metadata_milk_fat <- qPCR_Metadata[,c("lait_hier","fromage","pr_laitier","huile_viande")]
for(i in colnames(qPCR_Metadata_milk_fat)){
  qPCR_Metadata_milk_fat[,i] <- ifelse(qPCR_Metadata_milk_fat[,i] == "Oui",1,0)
  qPCR_Metadata_milk_fat[,i] <- as.numeric(qPCR_Metadata_milk_fat[,i])
}


qPCR_Metadata$diet_score_milk_cheese_fat <- rowSums(qPCR_Metadata_milk_fat)


qPCR_Metadata %>% sample_n_by(pays, Cq_presence_absence, stunted, size = 1)
qPCR_Metadata %>%
  group_by(pays, Cq_presence_absence, stunted) %>%
  get_summary_stats(diet_score_milk_cheese_fat, type = "mean_sd")


bxp <- ggboxplot(
  qPCR_Metadata, x = "stunted", y = "diet_score_milk_cheese_fat", 
  color = "Cq_presence_absence", palette = c("blue","orange"), facet.by = "pays",add = "jitter",legend = "R")
bxp


qPCR_Metadata %>%
  group_by(pays, Cq_presence_absence, stunted) %>%
  identify_outliers(diet_score_milk_cheese_fat)


model  <- lm(diet_score_milk_cheese_fat ~ pays*Cq_presence_absence*stunted, data = qPCR_Metadata)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))



qPCR_Metadata %>%
  group_by(pays, Cq_presence_absence, stunted) %>%
  shapiro_test(diet_score_milk_cheese_fat)



ggqqplot(qPCR_Metadata, "diet_score_milk_cheese_fat", ggtheme = theme_bw()) +
  facet_grid(pays + Cq_presence_absence ~ stunted, labeller = "label_both")


qPCR_Metadata %>% levene_test(diet_score_milk_cheese_fat ~ pays*Cq_presence_absence*stunted)


res.aov <- qPCR_Metadata %>% anova_test(diet_score_milk_cheese_fat ~ pays*Cq_presence_absence*stunted)
res.aov


# Group the data by pays and 
# fit simple two-way interaction 
model  <- lm(diet_score_milk_cheese_fat ~ pays*Cq_presence_absence*stunted, data = qPCR_Metadata)
qPCR_Metadata %>%
  group_by(pays) %>%
  anova_test(diet_score_milk_cheese_fat ~ Cq_presence_absence*stunted, error = model)

# Group the data by pays and Cq_presence_absence, and fit  anova
stunted.effect <- qPCR_Metadata %>%
  group_by(pays, stunted) %>%
  anova_test(diet_score_milk_cheese_fat ~ Cq_presence_absence, error = model)
stunted.effect

# Pairwise comparisons
library(emmeans)
pwc <- qPCR_Metadata %>%
  group_by(pays,stunted) %>%
  emmeans_test(diet_score_milk_cheese_fat ~ Cq_presence_absence, p.adjust.method = "BH")
# Show comparison results for male at high Cq_presence_absence
pwc %>% filter(stunted == "Stunted")

# Estimated marginal means (i.e. adjusted means) 
# with 95% confidence interval
get_emmeans(pwc) %>% filter(stunted == "Stunted")


# Visualization: box plots with p-values
pwc <- pwc %>% add_xy_position(x = "stunted")
pwc.filtered <- pwc #%>% filter(stunted == "Stunted")
pwc.filtered$y.position <- c(6,6,6,6)
ANOVA_3WAY_milk_cheese_fat_Blastocystis_presence <- bxp +
  stat_pvalue_manual(pwc.filtered, 
                     tip.length = 0, 
                     step.increase = 0.1, 
                     step.group.by = "pays") +
  labs(subtitle = get_test_label(res.aov, detailed = TRUE),
       caption = get_pwc_label(pwc))
ANOVA_3WAY_milk_cheese_fat_Blastocystis_presence 




##################################################################

plot_combined <- ggarrange(ANOVA_3WAY_DietScore_Blastocystis_presence,
                           ANOVA_3WAY_vegetarian_Blastocystis_presence,
                           ANOVA_3WAY_meat_fat_Blastocystis_presence,
                           ANOVA_3WAY_meat_Blastocystis_presence,
                           ANOVA_3WAY_milk_Blastocystis_presence,
                           ANOVA_3WAY_milk_cheese_fat_Blastocystis_presence,
                           labels = c("A","B"), ncol = 2, nrow=3);plot_combined
ggsave2(plot_combined, filename = "/Users/vincebilly/Desktop/vince/phd/statistic/afribiota_vince/Simple_risk_factor/figure/3WAY_ANOVA_diet_score_Blastocystis.pdf", height = 17, width = 12)


##################################################################
## Linear regression Abundance Blastocystis and Shannon and NUmber of ASV
##################################################################

qPCR_Metadata_presence <- subset(qPCR_Metadata, Cq_presence_absence =="Presence_qPCR")

#Shapiro test

shapiro.test(qPCR_Metadata_presence$Cq)
hist(qPCR_Metadata_presence$Cq)
#Not normal distributed

ggplot(data = qPCR_Metadata_presence, aes(x = Cq, y = Number_ASVs16s_new)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("adj.R2" ,"p", "n"))) +
  geom_point(size=1) +
  theme_classic() +
  scale_x_reverse()

ggplot(data = qPCR_Metadata_presence, aes(x = Cq, y = Shannon16s_new)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("adj.R2" ,"p", "n"))) +
  geom_point(size=1) +
  theme_classic() +
  scale_x_reverse()

ggplot(data = qPCR_Metadata_presence, aes(x = Cq, y = Number_ASVs_euk)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("adj.R2" ,"p", "n"))) +
  geom_point(size=1) +
  theme_classic() +
  scale_x_reverse()

ggplot(data = qPCR_Metadata_presence, aes(x = Cq, y = Shannon_euk)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("adj.R2" ,"p", "n"))) +
  geom_point(size=1) +
  theme_classic() +
  scale_x_reverse()

person_test <- cor.test(qPCR_Metadata_presence$Cq,qPCR_Metadata_presence$Shannon_euk, method = c("pearson", "kendall", "spearman"))
summary(person_test)


ggplot(data = qPCR_Metadata_presence, aes(x = Cq, y = Number_ASVs)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("adj.R2" ,"p", "n"))) +
  geom_point(size=1) +
  facet_wrap(~stunted, ncol=4) + 
  theme_classic() +
  scale_x_reverse()

ggplot(data = qPCR_Metadata_presence, aes(x = Cq, y = Shannon)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("adj.R2" ,"p", "n"))) +
  geom_point(size=1) +
  facet_wrap(~stunted, ncol=4) + 
  theme_classic() +
  scale_x_reverse()

ggplot(data = qPCR_Metadata_presence, aes(x = Cq, y = Number_ASVs_euk)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("adj.R2" ,"p", "n"))) +
  geom_point(size=1) +
  facet_wrap(~stunted, ncol=4) + 
  theme_classic() +
  scale_x_reverse()

ggplot(data = qPCR_Metadata_presence, aes(x = Cq, y = Shannon_euk)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("adj.R2" ,"p", "n"))) +
  geom_point(size=1) +
  facet_wrap(~stunted, ncol=4) + 
  theme_classic() +
  scale_x_reverse()



ggplot(data = qPCR_Metadata_presence, aes(x = Cq, y = Shannon)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("adj.R2" ,"p", "n"))) +
  geom_point(size=1) +
  facet_wrap(~pays, ncol=4) + 
  theme_classic() +
  scale_x_reverse()

ggplot(data = qPCR_Metadata_presence, aes(x = Cq, y = Number_ASVs)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("adj.R2" ,"p", "n"))) +
  geom_point(size=1) +
  facet_wrap(~pays, ncol=4) + 
  theme_classic() +
  scale_x_reverse()

ggplot(data = qPCR_Metadata_presence, aes(x = Cq, y = Shannon_euk)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("adj.R2" ,"p", "n"))) +
  geom_point(size=1) +
  facet_wrap(~pays, ncol=4) + 
  theme_classic() +
  scale_x_reverse()




ggplot(data = qPCR_Metadata_presence, aes(x = Cq, y = Number_ASVs_euk)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("adj.R2" ,"p", "n"))) +
  geom_point(size=1) +
  facet_wrap(~pays, ncol=4) + 
  theme_classic() +
  scale_x_reverse()

Reg_NumberASVs_Cq <- ggplot(data = qPCR_Metadata_presence, aes(x = Cq, y = Number_ASVs16s_new)) +
  stat_poly_line(color="blue") +
  stat_poly_eq(use_label(c("adj.R2" ,"p", "n"))) +
  geom_point(size=1, color="blue") +
  facet_wrap(~pays+stunted, ncol=4) + 
  theme_classic() +
  scale_x_reverse() +
  theme(axis.title=element_text(size=15),
        axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        plot.title=element_text(size=15))

Reg_Shannon_Cq <- ggplot(data = qPCR_Metadata_presence, aes(x = Cq, y = Shannon16s_new)) +
  stat_poly_line(color="blue") +
  stat_poly_eq(use_label(c("adj.R2" ,"p", "n"))) +
  geom_point(size=1, color="blue") +
  facet_wrap(~pays+stunted, ncol=4) + 
  theme_classic() +
  scale_x_reverse() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_text(size=15),
        axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        plot.title=element_text(size=15))

plot_combined <- ggarrange(Reg_Shannon_Cq,Reg_NumberASVs_Cq,labels = c("A","B"), ncol = 1, nrow=2);plot_combined
ggsave2(plot_combined, filename = "/Users/vincebilly/Desktop/vince/phd/statistic/afribiota_vince/Simple_risk_factor/figure/Reg_Shannon_NumberASVs_Cq.pdf", height = 6, width = 10)

Reg_NumberASVs_Cq <- ggplot(data = qPCR_Metadata, aes(x = Cq_trans, y = Number_ASVs)) +
  stat_poly_line(color="darkgreen") +
  stat_poly_eq(use_label(c("adj.R2" ,"p", "n"))) +
  geom_point(size=1, color="darkgreen") +
  facet_wrap(~pays+stunted+Cq_presence_absence, ncol=8) + 
  theme_classic() +
  scale_x_reverse()

Reg_Shannon_Cq <- ggplot(data = qPCR_Metadata, aes(x = Cq, y = Shannon)) +
  stat_poly_line(color="darkgreen") +
  stat_poly_eq(use_label(c("adj.R2" ,"p", "n"))) +
  geom_point(size=1, color="darkgreen") +
  facet_wrap(~pays+stunted, ncol=4) + 
  theme_classic() +
  scale_x_reverse()



ggplot(data = qPCR_Metadata_presence, aes(x = Cq, y = Number_ASVs_euk)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("adj.R2" ,"p", "n"))) +
  geom_point(size=1) +
  facet_wrap(~pays+stunted, ncol=4) + 
  theme_classic() +
  scale_x_reverse()



ggplot(data = qPCR_Metadata_presence, aes(x = Cq, y = Shannon_euk)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("adj.R2" ,"p", "n"))) +
  geom_point(size=1) +
  facet_wrap(~pays+stunted, ncol=4) + 
  theme_classic() +
  scale_x_reverse()

##################################################################

##################################################################
## Linear regression Abundance Blastocystis and haz
##################################################################

qPCR_Metadata_presence <- subset(qPCR_Metadata, Cq_presence_absence =="Presence_qPCR")

#Shapiro test

shapiro.test(qPCR_Metadata_presence$Cq)
hist(qPCR_Metadata_presence$Cq)
#Not normal distributed

qPCR_Metadata_presence$haz_cont <- as.numeric(qPCR_Metadata_presence$haz_cont)
ggplot(data = qPCR_Metadata_presence, aes(x = Cq, y = haz_cont)) +
  stat_poly_line(color="darkgreen") +
  stat_poly_eq(use_label(c("adj.R2" ,"p", "n"))) +
  geom_point(size=1, color="darkgreen") +
  facet_wrap(~pays+stunted, ncol=8) + 
  theme_classic() +
  scale_x_reverse()


##################################################################


##################################################################
## Madagascar Stunting
##################################################################

qPCR_Metadata_Mada_stunted <- subset(qPCR_Metadata, pays=="Madagascar" & stunted == "Stunted")


ggplot(qPCR_Metadata, aes(x=Cq_presence_absence, y=calprotectinelevel)) +
  geom_boxplot() +
  geom_jitter()

p <- ggboxplot(qPCR_Metadata_Mada_stunted, x = "Cq_presence_absence", y = "calprotectinelevel",
               color = "Cq_presence_absence", palette = "jco",
               add = "jitter")
#  Add p-value
p + stat_compare_means()

p <- ggboxplot(qPCR_Metadata_Mada_stunted, x = "Cq_presence_absence", y = "CALPROTECTINEggdePS",
               color = "Cq_presence_absence", palette = "jco",
               add = "jitter")
#  Add p-value
p + stat_compare_means()

qPCR_Metadata_Mada_stunted$age <- as.numeric(qPCR_Metadata_Mada_stunted$age)
p <- ggboxplot(qPCR_Metadata_Mada_stunted, x = "Cq_presence_absence", y = "age",
               color = "Cq_presence_absence", palette = "jco",
               add = "jitter")
#  Add p-value
p + stat_compare_means()

p <- ggboxplot(qPCR_Metadata_Mada_stunted, x = "Cq_presence_absence", y = "alphaantitrypsinlevel",
               color = "Cq_presence_absence", palette = "jco",
               add = "jitter")
#  Add p-value
p + stat_compare_means()

p <- ggboxplot(qPCR_Metadata_Mada_stunted, x = "Cq_presence_absence", y = "AATmggdePS",
               color = "Cq_presence_absence", palette = "jco",
               add = "jitter")
#  Add p-value
p + stat_compare_means()


p <- ggboxplot(qPCR_Metadata_Mada_stunted, x = "Cq_presence_absence", y = "Number_ASVs",
               color = "Cq_presence_absence", palette = "jco",
               add = "jitter")
#  Add p-value
p + stat_compare_means()

p <- ggboxplot(qPCR_Metadata_Mada_stunted, x = "Cq_presence_absence", y = "Shannon",
               color = "Cq_presence_absence", palette = "jco",
               add = "jitter")
#  Add p-value
p + stat_compare_means()

p <- ggboxplot(qPCR_Metadata_Mada_stunted, x = "Cq_presence_absence", y = "Simpson",
               color = "Cq_presence_absence", palette = "jco",
               add = "jitter")
#  Add p-value
p + stat_compare_means()

qPCR_Metadata_Mada_stunted$haz_cont <- as.numeric(qPCR_Metadata_Mada_stunted$haz_cont)
p <- ggboxplot(qPCR_Metadata_Mada_stunted, x = "Cq_presence_absence", y = "haz_cont",
               color = "Cq_presence_absence", palette = "jco",
               add = "jitter")
#  Add p-value
p + stat_compare_means()

qPCR_Metadata_Mada_stunted$whz_cont <- as.numeric(qPCR_Metadata_Mada_stunted$whz_cont)
p <- ggboxplot(qPCR_Metadata_Mada_stunted, x = "Cq_presence_absence", y = "whz_cont",
               color = "Cq_presence_absence", palette = "jco",
               add = "jitter")
#  Add p-value
p + stat_compare_means()

qPCR_Metadata_Mada_stunted$ferritine_non_corr <- as.numeric(qPCR_Metadata_Mada_stunted$ferritine_non_corr)
qPCR_Metadata_Mada_stunted_ferr <- subset(qPCR_Metadata_Mada_stunted, ferritine_non_corr<200)
p <- ggboxplot(qPCR_Metadata_Mada_stunted_ferr, x = "Cq_presence_absence", y = "ferritine_non_corr",
               color = "Cq_presence_absence", palette = "jco",
               add = "jitter")
#  Add p-value
p + stat_compare_means()

qPCR_Metadata_Mada_stunted$hemoglobine <- as.numeric(qPCR_Metadata_Mada_stunted$hemoglobine)
p <- ggboxplot(qPCR_Metadata_Mada_stunted, x = "Cq_presence_absence", y = "hemoglobine",
               color = "Cq_presence_absence", palette = "jco",
               add = "jitter")
#  Add p-value
p + stat_compare_means()


qPCR_Metadata_Mada_stunted$age_allaite <- as.numeric(qPCR_Metadata_Mada_stunted$age_allaite)
p <- ggboxplot(qPCR_Metadata_Mada_stunted, x = "Cq_presence_absence", y = "age_allaite",
               color = "Cq_presence_absence", palette = "jco",
               add = "jitter")
#  Add p-value
p + stat_compare_means()

qPCR_Metadata_Mada_stunted$age_alim <- as.numeric(qPCR_Metadata_Mada_stunted$age_alim)
p <- ggboxplot(qPCR_Metadata_Mada_stunted, x = "Cq_presence_absence", y = "age_alim",
               color = "Cq_presence_absence", palette = "jco",
               add = "jitter")
#  Add p-value
p + stat_compare_means()

p <- ggboxplot(qPCR_Metadata_Mada_stunted, x = "Cq_presence_absence", y = "Shannon_euk",
               color = "Cq_presence_absence", palette = "jco",
               add = "jitter")
#  Add p-value
p + stat_compare_means()

p <- ggboxplot(qPCR_Metadata, x = "Cq_presence_absence", y = "Number_ASVs_euk",
               color = "Cq_presence_absence", palette = "jco",
               add = "jitter")
#  Add p-value
p + stat_compare_means()

p <- ggboxplot(qPCR_Metadata, x = "Cq_presence_absence", y = "Shannon_euk",
               color = "Cq_presence_absence", palette = "jco",
               add = "jitter")
#  Add p-value
p + stat_compare_means()
table <- table(qPCR_Metadata_Mada_stunted$Cq_presence_absence,qPCR_Metadata_Mada_stunted$entamoeba)
fisher.test(table)

table <- as.data.frame(table)
colnames(table) <- c("Country","Stunted","Number_samples")

colours <- c("darkgreen","red","purple")

ggplot(table, aes(x=Country, y =Number_samples , fill=Stunted, color = Stunted )) +
  geom_bar(stat="identity", width=1, color="black", linetype = "solid") +  #height of the column equal the value, stacked
  scale_fill_manual(values=colours) +
  theme_bw()

##################################################################

##################################################################
###                Microbiome analysis 16s                   ###
##################################################################

otu_filtered <- as.data.frame(otu_table(project_data.rarefied_initial));rowSums(otu_filtered)
tax_filtered <- tax_table(project_data.rarefied_initial)

# Create Phyloseq Object
otu_filtered_phylo = otu_table(otu_filtered, taxa_are_rows = FALSE)
tax_filtered_phylo = tax_table(tax_filtered)
row.names(qPCR_Metadata) <- qPCR_Metadata$Sample
metadata_filtered_phylo <- sample_data(qPCR_Metadata)
subset_project_data <- phyloseq(otu_filtered_phylo,tax_filtered_phylo,metadata_filtered_phylo)


NMDS.bray <- ordinate(physeq = subset_project_data, method = "NMDS", distance = "bray")
NMDS <- as.data.frame(sample_data(subset_project_data))
bray <- as.data.frame(NMDS.bray$points)
row.names(bray) == row.names(NMDS) #sanity check #tests as true
NMDS$NMDS.bray1 <- bray$MDS1
NMDS$NMDS.bray2 <- bray$MDS2
#OPTIONAL: sort data for better looking easier to read plots
NMDS.sort <- NMDS[order(NMDS$pays, NMDS$stunted),] #if you do create a sorted plotting object, remember to modify the plots below so that you are using it as the input, not the unsorted tabl
## Create Colour Palettes ###
#n <- length(unique(NMDS$gorilla_name))
#qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
#col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))[-16]
#plain NMDS plot colored by "FACTOR_1" and shaped by "FACTOR_2"

NMDS_pays_bray <- ggplot(NMDS, aes(x=NMDS.bray1, y=NMDS.bray2, color = Cq_presence_absence)) + 
  geom_point(size=1) +
  #ggtitle(as.character(title)) +
  #geom_text(data=NMDS, aes(label=gorilla_name), vjust = 2) +
  scale_color_manual(values=c("blue","orange")) +
  scale_shape_manual(values=c(1,15)) +
  facet_wrap(~pays+stunted, ncol = 4) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")) +
  theme_classic()
print(NMDS_pays_bray)

NMDS.bray <- ordinate(physeq = subset_project_data, method = "NMDS", distance = "jaccard")
NMDS <- as.data.frame(sample_data(subset_project_data))
bray <- as.data.frame(NMDS.bray$points)
row.names(bray) == row.names(NMDS) #sanity check #tests as true
NMDS$NMDS.bray1 <- bray$MDS1
NMDS$NMDS.bray2 <- bray$MDS2
#OPTIONAL: sort data for better looking easier to read plots
NMDS.sort <- NMDS[order(NMDS$pays, NMDS$stunted),] #if you do create a sorted plotting object, remember to modify the plots below so that you are using it as the input, not the unsorted tabl
## Create Colour Palettes ###
#n <- length(unique(NMDS$gorilla_name))
#qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
#col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))[-16]
#plain NMDS plot colored by "FACTOR_1" and shaped by "FACTOR_2"

NMDS_pays_jaccard <- ggplot(NMDS, aes(x=NMDS.bray1, y=NMDS.bray2, color = Cq_presence_absence)) + 
  geom_point(size=1) +
  #ggtitle(as.character(title)) +
  #geom_text(data=NMDS, aes(label=gorilla_name), vjust = 2) +
  scale_color_manual(values=c("dodgerblue","orange")) +
  scale_shape_manual(values=c(1,15)) +
  facet_wrap(~pays+stunted, ncol = 4) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")) +
  theme_classic()
print(NMDS_pays_jaccard)

NMDS_partitioned <- ggarrange(NMDS_pays_bray,NMDS_pays_jaccard,ncol=1, nrow=2)
ggsave2(NMDS_partitioned, filename = "/Users/vincebilly/Desktop/vince/phd/statistic/afribiota_vince/Simple_risk_factor/figure/nmds_partitioned_presence.pdf", height = 6, width = 8)




project_data.bray <- phyloseq::distance(subset_project_data, method = "bray")
sample_df <- as.data.frame(as.matrix(sample_data(subset_project_data)))


res.adonis.rarefied.bray <- adonis2(project_data.bray ~ stunted, data =sample_df, method="bray");res.adonis.rarefied.bray
res.adonis.rarefied.bray <- adonis2(project_data.bray ~ pays, data =sample_df, method="bray");res.adonis.rarefied.bray
res.adonis.rarefied.bray <- adonis2(project_data.bray ~ Cq_presence_absence, data =sample_df, method="bray");res.adonis.rarefied.bray


res.adonis.rarefied.bray <- adonis2(project_data.bray ~ stunted*pays, data =sample_df, method="bray");res.adonis.rarefied.bray
res.adonis.rarefied.bray <- adonis2(project_data.bray ~ pays*stunted, data =sample_df, method="bray");res.adonis.rarefied.bray

res.adonis.rarefied.bray <- adonis2(project_data.bray ~ pays*Cq_presence_absence, data =sample_df, method="bray");res.adonis.rarefied.bray
res.adonis.rarefied.bray <- adonis2(project_data.bray ~ Cq_presence_absence*pays*, data =sample_df, method="bray");res.adonis.rarefied.bray

res.adonis.rarefied.bray <- adonis2(project_data.bray ~ stunted*Cq_presence_absence, data =sample_df, method="bray");res.adonis.rarefied.bray
res.adonis.rarefied.bray <- adonis2(project_data.bray ~ Cq_presence_absence*stunted*, data =sample_df, method="bray");res.adonis.rarefied.bray


res.adonis.rarefied.bray <- adonis2(project_data.bray ~ stunted*pays*Cq_presence_absence, data =sample_df, method="bray");res.adonis.rarefied.bray
res.adonis.rarefied.bray <- adonis2(project_data.bray ~ stunted*Cq_presence_absence*pays, data =sample_df, method="bray");res.adonis.rarefied.bray
res.adonis.rarefied.bray <- adonis2(project_data.bray ~ Cq_presence_absence*stunted*pays, data =sample_df, method="bray");res.adonis.rarefied.bray
res.adonis.rarefied.bray <- adonis2(project_data.bray ~ Cq_presence_absence*pays*stunted, data =sample_df, method="bray");res.adonis.rarefied.bray
res.adonis.rarefied.bray <- adonis2(project_data.bray ~ stunted*Cq_presence_absence*pays, data =sample_df, method="bray");res.adonis.rarefied.bray

res.adonis.rarefied.bray <- adonis2(project_data.bray ~ pays|Cq_presence_absence, data =sample_df, method="bray");res.adonis.rarefied.bray





Pvalue = round(res.adonis.rarefied.bray$aov.tab$`Pr(>F)`[1],3)
Statistic = round(res.adonis.rarefied.bray$aov.tab$R2[1],3)

NMDS.bray <- ordinate(physeq = project_data2, method = "NMDS", distance = "bray")
NMDS <- as.data.frame(sample_data(project_data2))
bray <- as.data.frame(NMDS.bray$points)
row.names(bray) == row.names(NMDS) #sanity check #tests as true
NMDS$NMDS.bray1 <- bray$MDS1
NMDS$NMDS.bray2 <- bray$MDS2
#OPTIONAL: sort data for better looking easier to read plots
NMDS.sort <- NMDS[order(NMDS$pays, NMDS$stunted),] #if you do create a sorted plotting object, remember to modify the plots below so that you are using it as the input, not the unsorted tabl
## Create Colour Palettes ###
#n <- length(unique(NMDS$gorilla_name))
#qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
#col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))[-16]
#plain NMDS plot colored by "FACTOR_1" and shaped by "FACTOR_2"

NMDS_pays_bray <- ggplot(NMDS, aes(x=NMDS.bray1, y=NMDS.bray2, color = pays)) + 
  geom_point(size=1) +
  #ggtitle(as.character(title)) +
  #geom_text(data=NMDS, aes(label=gorilla_name), vjust = 2) +
  scale_color_manual(values=c("purple","pink")) +
  scale_shape_manual(values=c(1,15)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")) +
  theme_classic()

NMDS_stunted_bray <- ggplot(NMDS, aes(x=NMDS.bray1, y=NMDS.bray2, color = stunted)) + 
  geom_point(size=1) +
  #ggtitle(as.character(title)) +
  #geom_text(data=NMDS, aes(label=gorilla_name), vjust = 2) +
  scale_color_manual(values=c("darkgreen","red")) +
  scale_shape_manual(values=c(1,15)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")) +
  theme_classic()


NMDS_Cq_presence_absence_bray <- ggplot(NMDS, aes(x=NMDS.bray1, y=NMDS.bray2, color = Cq_presence_absence)) + 
  geom_point(size=1) +
  #ggtitle(as.character(title)) +
  #geom_text(data=NMDS, aes(label=gorilla_name), vjust = 2) +
  scale_color_manual(values=c("dodgerblue","orange")) +
  scale_shape_manual(values=c(1,15)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")) +
  theme_classic()

ggarrange(NMDS_pays_bray,NMDS_stunted_bray,NMDS_Cq_presence_absence_bray, ncol=3)



NMDS_site_bray <- ggplot(NMDS, aes(x=NMDS.bray1, y=NMDS.bray2, color = Cq_presence_absence, shape = stunted)) + 
  geom_point(size=1) +
  #ggtitle(as.character(title)) +
  #geom_text(data=NMDS, aes(label=gorilla_name), vjust = 2) +
  facet_wrap(~ pays + stunted, ncol = 4) +
  scale_color_manual(values=c("dodgerblue","orange")) +
  scale_shape_manual(values=c(1,15)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")) +
  theme_classic()

#PERMANOVA with presence and Absence
list_NMDS <- list()
for(cntr in c("Madagascar","RCA")){
  subset_project_data1 <- subset_samples(subset_project_data, pays == cntr)
  for (stnd in c("Stunted","Non-Stunted")){
    subset_project_data2 <- subset_samples(subset_project_data1, stunted == stnd)
    sample_df <- as.data.frame(as.matrix(sample_data(subset_project_data2)))
    project_data.bray <- phyloseq::distance(subset_project_data2, method = "bray")
    res.adonis.rarefied.bray <- adonis2(project_data.bray ~ Cq_presence_absence, data =sample_df, method="bray")
    project_data.Jaccard  <- phyloseq::distance(subset_project_data2, method = "jaccard")
    res.adonis.rarefied.Jaccard <- adonis2(project_data.Jaccard  ~ Cq_presence_absence, data =sample_df, method="jaccard")
    table_summary <- data_frame(Country =cntr,
                                Stunting = stnd,
                                Presence = table(sample_df$Cq_presence_absence)[2],
                                Absence = table(sample_df$Cq_presence_absence)[1],
                                Df_Bray = res.adonis.rarefied.bray$Df[1],
                                F_Bray = round(res.adonis.rarefied.bray$F[1],3),
                                RSquared_Bray = round(res.adonis.rarefied.bray$R2[1],3),
                                Pval_Bray = res.adonis.rarefied.bray$`Pr(>F)`[1],
                                Df_Jaccard  = res.adonis.rarefied.Jaccard$Df[1],
                                F_Jaccard  = round(res.adonis.rarefied.Jaccard$F[1],3),
                                RSquared_Jaccard  = round(res.adonis.rarefied.Jaccard$R2[1],3),
                                Pval_Jaccard  = res.adonis.rarefied.Jaccard$`Pr(>F)`[1])
    list_NMDS[[length(list_NMDS)+1]] <- table_summary
  }
}

summary_table <- data.frame()
for(i in 1:length(list_NMDS)){
  summary_table <- rbind(summary_table,list_NMDS[[i]])
}

write_csv(summary_table,"/Users/vincebilly/Desktop/vince/phd/statistic/afribiota_vince/Simple_risk_factor/figure/PERMANOVA.csv")

#PERMANOVA with abundance
subset_project_data_presence <-  subset_samples(subset_project_data, Cq_presence_absence =="Presence_qPCR")
list_NMDS <- list()
for(cntr in c("Madagascar","RCA")){
  subset_project_data1 <- subset_samples(subset_project_data_presence, pays == cntr)
  for (stnd in c("Stunted","Non-Stunted")){
    subset_project_data2 <- subset_samples(subset_project_data1, stunted == stnd)
    sample_df <- as.data.frame(as.matrix(sample_data(subset_project_data2)))
    sample_df$Cq <- as.numeric(sample_df$Cq)
    project_data.bray <- phyloseq::distance(subset_project_data2, method = "bray")
    res.adonis.rarefied.bray <- adonis2(project_data.bray ~ Cq, data =sample_df, method="bray")
    project_data.Jaccard  <- phyloseq::distance(subset_project_data2, method = "jaccard")
    res.adonis.rarefied.Jaccard <- adonis2(project_data.Jaccard  ~ Cq, data =sample_df, method="jaccard")
    table_summary <- data_frame(Country = cntr,
                                Stunting = stnd,
                                Presence = length(sample_df$Cq_presence_absence),
                                Absence = 0,
                                Df_Bray = res.adonis.rarefied.bray$Df[1],
                                F_Bray = round(res.adonis.rarefied.bray$F[1],3),
                                RSquared_Bray = round(res.adonis.rarefied.bray$R2[1],3),
                                Pval_Bray = round(res.adonis.rarefied.bray$`Pr(>F)`[1],3),
                                Df_Jaccard  = res.adonis.rarefied.Jaccard$Df[1],
                                F_Jaccard  = round(res.adonis.rarefied.Jaccard$F[1],3),
                                RSquared_Jaccard  = round(res.adonis.rarefied.Jaccard$R2[1],3),
                                Pval_Jaccard  = round(res.adonis.rarefied.Jaccard$`Pr(>F)`[1],3))
    list_NMDS[[length(list_NMDS)+1]] <- table_summary
  }
}

summary_table_NMDS_abundance <- data.frame()
for(i in 1:length(list_NMDS)){
  summary_table_NMDS_abundance <- rbind(summary_table_NMDS_abundance,list_NMDS[[i]])
}

write_csv(summary_table_NMDS_abundance,"/Users/vincebilly/Desktop/vince/phd/statistic/afribiota_vince/Simple_risk_factor/figure/PERMANOVA_abundance.csv")


##################################################################

##################################################################
###          Microbiome analysis 16s abundance                 ###
##################################################################

qPCR_Metadata_presence <- subset(qPCR_Metadata, Cq_presence_absence =="Presence_qPCR")
qPCR_Metadata_presence <- subset(qPCR_Metadata_presence, Cq>10)

otu_filtered <- as.data.frame(t(otu_table(project_data.rarefied_initial)));rowSums(otu_filtered)
tax_filtered <- tax_table(project_data.rarefied_initial)

row.names(otu_filtered)==row.names(tax_filtered)
otu_tax_filtered <- merge(otu_filtered)

# Create Phyloseq Object
otu_filtered_phylo = otu_table(otu_filtered, taxa_are_rows = FALSE)
tax_filtered_phylo = tax_table(tax_filtered)
row.names(qPCR_Metadata_presence) <- qPCR_Metadata_presence$Sample
metadata_filtered_phylo <- sample_data(qPCR_Metadata_presence)
subset_project_data <- phyloseq(otu_filtered_phylo,tax_filtered_phylo,metadata_filtered_phylo)

otu <- otu_table

##################################################################



##################################################################
###                Microbiome analysis 18s                   ###
##################################################################

otu_filtered <- as.data.frame(otu_table(project_data_euk));rowSums(otu_filtered)
tax_filtered <- tax_table(project_data_euk)


# Create Phyloseq Object
otu_filtered_phylo = otu_table(otu_filtered, taxa_are_rows = TRUE)
tax_filtered_phylo = tax_table(tax_filtered)
row.names(qPCR_Metadata) <- qPCR_Metadata$Sample
metadata_filtered_phylo <- sample_data(qPCR_Metadata)
subset_project_data <- phyloseq(otu_filtered_phylo,tax_filtered_phylo,metadata_filtered_phylo)


project_data.bray <- phyloseq::distance(subset_project_data, method = "bray")
sample_df <- as.data.frame(as.matrix(sample_data(subset_project_data)))


res.adonis.rarefied.bray <- adonis2(project_data.bray ~ stunted, data =sample_df, method="bray");res.adonis.rarefied.bray
res.adonis.rarefied.bray <- adonis2(project_data.bray ~ pays, data =sample_df, method="bray");res.adonis.rarefied.bray
res.adonis.rarefied.bray <- adonis2(project_data.bray ~ Cq_presence_absence, data =sample_df, method="bray");res.adonis.rarefied.bray


res.adonis.rarefied.bray <- adonis(project_data.bray ~ stunted*pays, data =sample_df, method="bray");res.adonis.rarefied.bray
res.adonis.rarefied.bray <- adonis(project_data.bray ~ pays*stunted, data =sample_df, method="bray");res.adonis.rarefied.bray

res.adonis.rarefied.bray <- adonis(project_data.bray ~ pays*Cq_presence_absence, data =sample_df, method="bray");res.adonis.rarefied.bray
res.adonis.rarefied.bray <- adonis(project_data.bray ~ Cq_presence_absence*pays*, data =sample_df, method="bray");res.adonis.rarefied.bray

res.adonis.rarefied.bray <- adonis(project_data.bray ~ stunted*Cq_presence_absence, data =sample_df, method="bray");res.adonis.rarefied.bray
res.adonis.rarefied.bray <- adonis(project_data.bray ~ Cq_presence_absence*stunted*, data =sample_df, method="bray");res.adonis.rarefied.bray


res.adonis.rarefied.bray <- adonis(project_data.bray ~ stunted*pays*Cq_presence_absence, data =sample_df, method="bray");res.adonis.rarefied.bray
res.adonis.rarefied.bray <- adonis(project_data.bray ~ stunted*Cq_presence_absence*pays, data =sample_df, method="bray");res.adonis.rarefied.bray
res.adonis.rarefied.bray <- adonis(project_data.bray ~ Cq_presence_absence*stunted*pays, data =sample_df, method="bray");res.adonis.rarefied.bray
res.adonis.rarefied.bray <- adonis(project_data.bray ~ Cq_presence_absence*pays*stunted, data =sample_df, method="bray");res.adonis.rarefied.bray
res.adonis.rarefied.bray <- adonis(project_data.bray ~ stunted*Cq_presence_absence*pays, data =sample_df, method="bray");res.adonis.rarefied.bray
res.adonis.rarefied.bray <- adonis(project_data.bray ~ stunted*pays*Cq_presence_absence, data =sample_df, method="bray");res.adonis.rarefied.bray





Pvalue = round(res.adonis.rarefied.bray$aov.tab$`Pr(>F)`[1],3)
Statistic = round(res.adonis.rarefied.bray$aov.tab$R2[1],3)

NMDS.bray <- ordinate(physeq = subset_project_data, method = "NMDS", distance = "bray")
NMDS <- as.data.frame(sample_data(project_data2))
bray <- as.data.frame(NMDS.bray$points)
row.names(bray) == row.names(NMDS) #sanity check #tests as true
NMDS$NMDS.bray1 <- bray$MDS1
NMDS$NMDS.bray2 <- bray$MDS2
#OPTIONAL: sort data for better looking easier to read plots
NMDS.sort <- NMDS[order(NMDS$pays, NMDS$stunted),] #if you do create a sorted plotting object, remember to modify the plots below so that you are using it as the input, not the unsorted tabl
## Create Colour Palettes ###
#n <- length(unique(NMDS$gorilla_name))
#qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
#col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))[-16]
#plain NMDS plot colored by "FACTOR_1" and shaped by "FACTOR_2"

NMDS_pays_bray <- ggplot(NMDS, aes(x=NMDS.bray1, y=NMDS.bray2, color = pays)) + 
  geom_point(size=1) +
  #ggtitle(as.character(title)) +
  #geom_text(data=NMDS, aes(label=gorilla_name), vjust = 2) +
  scale_color_manual(values=c("purple","pink")) +
  scale_shape_manual(values=c(1,15)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")) +
  theme_classic()

NMDS_stunted_bray <- ggplot(NMDS, aes(x=NMDS.bray1, y=NMDS.bray2, color = stunted)) + 
  geom_point(size=1) +
  #ggtitle(as.character(title)) +
  #geom_text(data=NMDS, aes(label=gorilla_name), vjust = 2) +
  scale_color_manual(values=c("darkgreen","red")) +
  scale_shape_manual(values=c(1,15)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")) +
  theme_classic()


NMDS_Cq_presence_absence_bray <- ggplot(NMDS, aes(x=NMDS.bray1, y=NMDS.bray2, color = Cq_presence_absence)) + 
  geom_point(size=1) +
  #ggtitle(as.character(title)) +
  #geom_text(data=NMDS, aes(label=gorilla_name), vjust = 2) +
  scale_color_manual(values=c("dodgerblue","orange")) +
  scale_shape_manual(values=c(1,15)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")) +
  theme_classic()

ggarrange(NMDS_pays_bray,NMDS_stunted_bray,NMDS_Cq_presence_absence_bray, ncol=3)



NMDS_site_bray <- ggplot(NMDS, aes(x=NMDS.bray1, y=NMDS.bray2, color = Cq_presence_absence, shape = stunted)) + 
  geom_point(size=1) +
  #ggtitle(as.character(title)) +
  #geom_text(data=NMDS, aes(label=gorilla_name), vjust = 2) +
  facet_wrap(~ pays + stunted, ncol = 4) +
  scale_color_manual(values=c("dodgerblue","orange")) +
  scale_shape_manual(values=c(1,15)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")) +
  theme_classic()


subset_project_data2 <- subset_project_data %>%
  subset_samples(pays == "Madagascar") %>%
  subset_samples(stunted == "Non")
project_data.bray <- phyloseq::distance(subset_project_data2, method = "bray")
sample_df <- as.data.frame(as.matrix(sample_data(subset_project_data2)))
res.adonis.rarefied.bray <- adonis(project_data.bray ~ Cq_presence_absence, data =sample_df, method="bray");res.adonis.rarefied.bray

subset_project_data2 <- subset_project_data %>%
  subset_samples(pays == "Madagascar") %>%
  subset_samples(stunted == "Oui")
project_data.bray <- phyloseq::distance(subset_project_data2, method = "bray")
sample_df <- as.data.frame(as.matrix(sample_data(subset_project_data2)))
res.adonis.rarefied.bray <- adonis(project_data.bray ~ Cq_presence_absence, data =sample_df, method="bray");res.adonis.rarefied.bray

subset_project_data2 <- subset_project_data %>%
  subset_samples(pays == "RCA") %>%
  subset_samples(stunted == "Non")
project_data.bray <- phyloseq::distance(subset_project_data2, method = "bray")
sample_df <- as.data.frame(as.matrix(sample_data(subset_project_data2)))
res.adonis.rarefied.bray <- adonis(project_data.bray ~ Cq_presence_absence, data =sample_df, method="bray");res.adonis.rarefied.bray

subset_project_data2 <- subset_project_data %>%
  subset_samples(pays == "RCA") %>%
  subset_samples(stunted == "Oui")
project_data.bray <- phyloseq::distance(subset_project_data2, method = "bray")
sample_df <- as.data.frame(as.matrix(sample_data(subset_project_data2)))
res.adonis.rarefied.bray <- adonis(project_data.bray ~ Cq_presence_absence, data =sample_df, method="bray");res.adonis.rarefied.bray


##################################################################

####################################################################################################################################
####################################################################################################################################
###                                                       DESeq analysis                                                         ###
####################################################################################################################################
####################################################################################################################################

#######################################################
### Run DESeq in all countries and stunted status   ###
#######################################################

pays_vector <- c("Both_Country","Madagascar","RCA")
stunted_vector <- c("Both_clinical","Stunted","Non-Stunted")
summary_table <- data.frame(Cplt_taxo=c(NA))
for(p in 1:length(pays_vector)){
  if (pays_vector[p]=="Both_Country"){
    project_data_1 <- subset_project_data
  } else {
    project_data_1 <- subset_project_data %>%
      subset_samples(pays == pays_vector[p])
  }
  for(s in 1:length(stunted_vector)){
    if (stunted_vector[s]=="Both_clinical"){
      project_data_2 <- project_data_1
    } else {
      project_data_2 <- project_data_1 %>%
        subset_samples(stunted == stunted_vector[s])
    }
    deseq_table <- data.frame()
    lefse_table <- data.frame()                                                     #Create dataframe
    table_bubble_plot <- data.frame()
    ############################################################
    ### Aggregation at each rank individually (from 1 to 7)  ###
    ############################################################
    rank_vector <- c("Rank2","Rank3","Rank4","Rank5","Rank6","Rank7","Accession")
    for(rnk in 1:length(rank_vector)){
      if (rank_vector[rnk]=="Accession"){
        project_data_1_rank <- project_data_2
      } else {
        project_data_1_rank <- tax_glom(project_data_2, taxrank=rank_vector[rnk])
      }
      ######################################
      ###       DESEQ2 using WALD        ###
      ######################################
      project_data_1_rank_deseq <- project_data_1_rank
      otu_table(project_data_1_rank_deseq) <- otu_table(project_data_1_rank_deseq) + 1
      sample_data(project_data_1_rank_deseq)$Cq_presence_absence <- factor(sample_data(project_data_1_rank_deseq)$Cq_presence_absence, levels=c("Presence_qPCR","Absence_qPCR")) #Controling level order here
      dds.var3 <- phyloseq_to_deseq2(project_data_1_rank_deseq, design = ~ Cq_presence_absence)
      dds.var3 <- DESeq(dds.var3, test = "Wald", fitType = "parametric")
      alpha <- 0.05 #Set significance threshold for multiple test corrected p-values
      res.var3 <- results(dds.var3, cooksCutoff = FALSE, alpha = alpha, pAdjustMethod = "BH", contrast=c("Cq_presence_absence","Presence_qPCR", "Absence_qPCR"), altHypothesis = "greaterAbs") #here "level2" is experimental and "level1" is control or baseline #here we can use "contrast" to specify the comparison we are interested in
      summary(res.var3, alpha = alpha) # view a summary of the results table with a padj value < 0.05
      ##filtering the results table #
      res_p_ordered <- res.var3[order(res.var3$padj, na.last = NA), ]           #reorder the results table by adjusted p-value and remove any "NA" entries
      res_p_ordered_filt <- res_p_ordered[which(res_p_ordered$padj < alpha), ]  #filter out any results which have a adjust p-value less than alpha (0.01)
      res_p_ordered_wtaxa.var3 <- cbind(as(res_p_ordered_filt, "data.frame"), as(tax_table(project_data_1_rank_deseq)[rownames(res_p_ordered_filt), ], "matrix")) #show results preserving taxa table
      res_p_ordered_wtaxa.var3$Cplt_taxo <- paste0(res_p_ordered_wtaxa.var3$Rank1,"|",
                                               res_p_ordered_wtaxa.var3$Rank2,"|",
                                               res_p_ordered_wtaxa.var3$Rank3,"|",
                                               res_p_ordered_wtaxa.var3$Rank4,"|",
                                               res_p_ordered_wtaxa.var3$Rank5,"|",
                                               res_p_ordered_wtaxa.var3$Rank6,"|",
                                               res_p_ordered_wtaxa.var3$Rank7,"|",
                                               row.names(res_p_ordered_wtaxa.var3),"|",
                                                         rank_vector[rnk])
      table_to_keep <- res_p_ordered_wtaxa.var3[,c("Cplt_taxo","log2FoldChange")]
      deseq_table <- rbind(deseq_table,table_to_keep)
      ######################################
      ###       Build LEfSe table        ###
      ######################################
      otu <- as.data.frame(t(otu_table(project_data_1_rank)))                       #Extract otu table
      tax <- as.data.frame(as(tax_table(project_data_1_rank), "matrix"))            #Extract tax table
      tax_otu <- merge(tax, otu, by ="row.names", all = T)                           #Merge tax and otu table
      lefse_table <- rbind(lefse_table, tax_otu)                                     #Concatenate result from each Rank
      table_full_bubble_plot <- lefse_table
      table_full_bubble_plot$Cplt_taxo <- paste0(table_full_bubble_plot$Rank1,"|",
                                                 table_full_bubble_plot$Rank2,"|",
                                                 table_full_bubble_plot$Rank3,"|",
                                                 table_full_bubble_plot$Rank4,"|",
                                                 table_full_bubble_plot$Rank5,"|",
                                                 table_full_bubble_plot$Rank6,"|",
                                                 table_full_bubble_plot$Rank7,"|",
                                                 table_full_bubble_plot$Row.names,"|",
                                                   rank_vector[rnk])
      table_bubble_plot <- rbind(table_bubble_plot,table_full_bubble_plot)
    }
    ######################################
    ###       DESEQ2 using WALD        ###
    ######################################
    deseq_table$Cplt_taxo <- gsub("[|]NA","",deseq_table$Cplt_taxo)
    summary_table <- merge(summary_table,deseq_table,by="Cplt_taxo", all=TRUE)
    colnames(summary_table)[dim(summary_table)[2]] <- paste0(pays_vector[p],"_",stunted_vector[s])
    ######################################
    ###       Build LEfSe table        ###
    ######################################
    #Combine taxa name in one column separated with pipe "|"
    lefse_table$Row.names <- paste0(lefse_table$Rank1,"|",lefse_table$Rank2,"|",lefse_table$Rank3,"|",lefse_table$Rank4,"|",lefse_table$Rank5,"|",lefse_table$Rank6,"|",lefse_table$Rank7,"|",lefse_table$Rank8)
    #Remove NA in the same column
    lefse_table$Row.names <- gsub("[|]NA","",lefse_table$Row.names)
    #Remove all "Rank" columns
    lefse_table <- lefse_table[,setdiff(colnames(lefse_table),colnames(tax_table(project_data_2)))]
    #Reformat table for LEfSe analysis input
    colnames(lefse_table)[1] <- "id"
    id <- colnames(lefse_table)
    sample_df <- as.data.frame(as.matrix(sample_data(project_data_2)))
    Treatment <- c("Treatment", as.character(sample_df$Cq_presence_absence))
    lefse_table <- rbind(id,lefse_table)
    colnames(lefse_table) <- Treatment
    #Save table
    setwd("/Users/vincebilly/Desktop/vince/phd/statistic/afribiota_vince/Simple_risk_factor/figure/lefse_input/")
    write.table(lefse_table, paste0("lefse_",pays_vector[p],"_",stunted_vector[s],".txt"), sep="\t", row.names = FALSE, quote=FALSE)
    ################################################################
    ###       Export full abundance table for bubble plot        ###
    ################################################################
    setwd("/Users/vincebilly/Desktop/vince/phd/statistic/afribiota_vince/Simple_risk_factor/figure/")
    write.table(table_bubble_plot, paste0("Full_abundance_Table_",pays_vector[p],"_",stunted_vector[s],".txt"), sep="\t", row.names = FALSE, quote=FALSE)
  }
}

#######################################################

#######################################################
###             Plot All DESeq analysis             ###
#######################################################

log_melt <- summary_table
  
#Create a another colum to indicuate what is significant or not
log_melt$Mada_pos_neg_stunted <- log_melt$Madagascar_Stunted
#Convert as significant value below or equal to 0.05
log_melt$Mada_pos_neg_stunted[log_melt$Mada_pos_neg_stunted>0]  <- "Increase"
log_melt$Mada_pos_neg_stunted[log_melt$Mada_pos_neg_stunted<0]  <- "Decrease"

#Create a another colum to indicuate what is significant or not
log_melt$Mada_pos_neg_non_stunted <- log_melt$`Madagascar_Non-Stunted`
#Convert as significant value below or equal to 0.05
log_melt$Mada_pos_neg_non_stunted[log_melt$Mada_pos_neg_non_stunted>0]  <- "Increase"
log_melt$Mada_pos_neg_non_stunted[log_melt$Mada_pos_neg_non_stunted<0]  <- "Decrease"

#Create a another colum to indicuate what is significant or not
log_melt$CAR_pos_neg_stunted <- log_melt$RCA_Stunted
#Convert as significant value below or equal to 0.05
log_melt$CAR_pos_neg_stunted[log_melt$CAR_pos_neg_stunted>0]  <- "Increase"
log_melt$CAR_pos_neg_stunted[log_melt$CAR_pos_neg_stunted<0]  <- "Decrease"

#Create a another colum to indicuate what is significant or not
log_melt$CAR_pos_neg_non_stunted <- log_melt$`RCA_Non-Stunted`
#Convert as significant value below or equal to 0.05
log_melt$CAR_pos_neg_non_stunted[log_melt$CAR_pos_neg_non_stunted>0]  <- "Increase"
log_melt$CAR_pos_neg_non_stunted[log_melt$CAR_pos_neg_non_stunted<0]  <- "Decrease"




log_melt3 <- log_melt
log_melt3$order <- paste0(log_melt3$Mada_pos_neg_stunted,"_",log_melt3$Mada_pos_neg_non_stunted,
                          "_",log_melt3$CAR_pos_neg_stunted,"_",log_melt3$CAR_pos_neg_non_stunted)
head(log_melt3)
log_melt3_melt <- melt(log_melt3, id.vars=c("Cplt_taxo","order"),
                    measure.vars = c("Both_Country_Both_clinical","Madagascar_Both_clinical","RCA_Both_clinical","Both_Country_Stunted","Both_Country_Non-Stunted","Madagascar_Stunted","Madagascar_Non-Stunted","RCA_Stunted","RCA_Non-Stunted"));head(log_melt3_mel)

log_melt3_mel <- log_melt3_melt
log_melt3_mel$direction <- ifelse(log_melt3_mel$value>0,"Increase","Decrease")
log_melt3_mel$value <- abs(log_melt3_mel$value)
log_melt3_mel$Cplt_taxo <- as.character(log_melt3$Cplt_taxo)

log_melt3_mel$order <- factor(log_melt3_mel$order, levels = c("Increase_Increase_Increase_Increase",
                                                              "Decrease_Decrease_Decrease_Decrease",
                                                              "Increase_NA_Increase_NA",
                                                              "Decrease_NA_Decrease_NA",
                                                              "NA_Increase_NA_Increase",
                                                              "NA_Decrease_NA_Decrease",
                                                              "Increase_Increase_NA_NA",
                                                              "Decrease_Decrease_NA_NA",
                                                              "NA_NA_Increase_Increase",
                                                              "NA_NA_Decrease_Decrease",
                                                              "Increase_Increase_Increase_NA",
                                                              "Decrease_Decrease_Decrease_NA",
                                                              "NA_Increase_Increase_Increase",
                                                              "NA_Decrease_Decrease_Decrease",
                                                              "Increase_Increase_NA_Increase",
                                                              "Decrease_Decrease_NA_Decrease",
                                                              "Increase_NA_Increase_Increase",
                                                              "Increase_NA_NA_Increase",
                                                              "Decrease_NA_NA_Decrease",
                                                              "NA_Increase_Increase_NA","NA_Decrease_Decrease_NA",
                                                              "NA_Increase_NA_NA","NA_Decrease_NA_NA",
                                                              "NA_NA_NA_Increase","NA_NA_NA_Decrease",
                                                              "NA_NA_Increase_NA","NA_NA_Decrease_NA",
                                                              "Increase_NA_NA_NA","Decrease_NA_NA_NA",
                                                              "NA_Decrease_NA_Increase","Increase_Decrease_NA_NA","NA_NA_Increase_Decrease","Increase_Decrease_Increase_NA","Decrease_NA_NA_Increase","NA_NA_Decrease_Increase","Increase_NA_Decrease_Decrease",
                                                              "Increase_NA_Increase_Decrease","NA_Decrease_Decrease_Increase","Increase_Decrease_Decrease_NA",
                                                              "NA_Increase_Increase_Decrease","Increase_NA_NA_Decrease","Increase_NA_Decrease_NA",
                                                              "Increase_Increase_NA_Decrease","Increase_Increase_Decrease_NA","NA_Decrease_Increase_NA",
                                                              "Decrease_Decrease_NA_Increase","NA_Increase_NA_Decrease",
                                                              "Decrease_NA_Increase_NA","Decrease_Increase_Increase_NA","Decrease_Increase_NA_NA",
                                                              "NA_NA_NA_NA","Increase_Decrease_Decrease_Decrease","Decrease_NA_Decrease_Increase" )  
)

which(log_melt3_mel$order=="Decrease_Decrease_Decrease_Decrease")
which(log_melt3_mel$order=="Decrease_NA_Decrease_NA")
which(log_melt3_mel$order=="NA_Decrease_NA_Decrease")
which(log_melt3_mel$order=="Decrease_Increase_Decrease_Increase")
which(log_melt3_mel$order=="Increase_Decrease_Increase_Decrease")
which(log_melt3_mel$order=="Decrease_Decrease_Decrease_NA")
which(log_melt3_mel$order=="Decrease_NA_Decrease_Decrease")
which(log_melt3_mel$order=="Decrease_Decrease_NA_Decrease")
which(log_melt3_mel$order=="NA_Decrease_Decrease_Decrease")
which(log_melt3_mel$order=="NA_Increase_Increase_Increase")
which(log_melt3_mel$order=="Decrease_NA_NA_Decrease")

pdf("/Users/vincebilly/Desktop/vince/phd/statistic/afribiota_vince/Simple_risk_factor/figure/DESeq_summary.pdf", height = 45, width = 20)
ggplot(log_melt3_mel, aes(x=variable, y=reorder(Cplt_taxo,-as.integer(factor(order)), FUN=min), size= value, color=direction)) +
  geom_point(alpha=0.5) +
  scale_color_manual(values=c("orange","blue"), name="Change in abundance") +
  theme(axis.text.x = element_text(angle = -45, hjust = 0))
dev.off()

#PLot only the most important taxa
to_keep <- c("Increase_Increase_Increase_Increase",
             "Decrease_Decrease_Decrease_Decrease",
             "Increase_NA_Increase_NA",
             "Decrease_NA_Decrease_NA",
             "NA_Increase_NA_Increase",
             "NA_Decrease_NA_Decrease")
log_melt4_mel <- subset(log_melt3_mel, order %in% to_keep)
pdf("/Users/vincebilly/Desktop/vince/phd/statistic/afribiota_vince/Simple_risk_factor/figure/DESeq_summary_most_important.pdf", height = 10, width = 20)
ggplot(log_melt4_mel, aes(x=variable, y=reorder(Cplt_taxo,-as.integer(factor(order)), FUN=min), size= value, color=direction)) +
  geom_point(alpha=0.5) +
  scale_color_manual(values=c("orange","blue"), name="Change in abundance") +
  facet_grid(order~., scales = "free_y", space = "free") +
  theme_light() +
  theme(axis.text.x = element_text(size = 15, angle = 45, hjust = 1, face="bold"),
        axis.text.y = element_text(size = 15, hjust = 0, face="bold"),
        strip.text.y = element_blank())
dev.off()

#PLot with only 4 dataset
to_keep <- c("Madagascar_Stunted",
             "Madagascar_Non-Stunted",
             "RCA_Stunted",
             "RCA_Non-Stunted")
log_melt5_mel <- subset(log_melt4_mel, variable %in% to_keep)
pdf("/Users/vincebilly/Desktop/vince/phd/statistic/afribiota_vince/Simple_risk_factor/figure/DESeq_summary_most_important_4_datasets.pdf", height = 10, width = 20)
ggplot(log_melt5_mel, aes(x=variable, y=reorder(Cplt_taxo,-as.integer(factor(order)), FUN=min), size= value, color=direction)) +
  geom_point(alpha=0.5) +
  scale_color_manual(values=c("orange","blue"), name="Change in abundance") +
  facet_grid(order~., scales = "free_y", space = "free") +
  theme_light() +
  theme(axis.text.x = element_text(size = 15, angle = 45, hjust = 1, face="bold"),
        axis.text.y = element_text(size = 15, hjust = 0, face="bold"),
        strip.text.y = element_blank(),
        legend.position="bottom") +
  guides(fill=guide_legend(nrow=1,byrow=TRUE))
dev.off()
#######################################################

####################################################################
###  Plot with only stunted and non stunted (countries combined) ###
####################################################################

log_melt <- summary_table[,c("Cplt_taxo","Both_Country_Stunted","Both_Country_Non-Stunted")]

#Create a another column to indicuate what is significant or not
log_melt$pos_neg_stunted <- log_melt$Both_Country_Stunted
#Convert as significant value below or equal to 0.05
log_melt$pos_neg_stunted[log_melt$pos_neg_stunted>0]  <- "Increase"
log_melt$pos_neg_stunted[log_melt$pos_neg_stunted<0]  <- "Decrease"

#Create a another colum to indicuate what is significant or not
log_melt$pos_neg_non_stunted <- log_melt$`Both_Country_Non-Stunted`
#Convert as significant value below or equal to 0.05
log_melt$pos_neg_non_stunted[log_melt$pos_neg_non_stunted>0]  <- "Increase"
log_melt$pos_neg_non_stunted[log_melt$pos_neg_non_stunted<0]  <- "Decrease"

log_melt3 <- log_melt
log_melt3$order <- paste0(log_melt3$pos_neg_stunted,"_",log_melt3$pos_neg_non_stunted)
head(log_melt3)
log_melt3_melt <- melt(log_melt3, id.vars=c("Cplt_taxo","order"),
                       measure.vars = c("Both_Country_Stunted","Both_Country_Non-Stunted"));head(log_melt3_mel)


log_melt3_mel <- log_melt3_melt
log_melt3_mel$direction <- ifelse(log_melt3_mel$value>0,"Increase","Decrease")
log_melt3_mel$value_in <- log_melt3_mel$value
log_melt3_mel$value <- abs(log_melt3_mel$value)
log_melt3_mel$Cplt_taxo <- as.character(log_melt3$Cplt_taxo)

log_melt3_mel$order <- factor(log_melt3_mel$order, levels = c("Increase_Increase",
                                                              "Decrease_Decrease",
                                                              "Increase_NA",
                                                              "Decrease_NA",
                                                              "NA_Increase",
                                                              "NA_Decrease",
                                                              "Decrease_Increase",
                                                              "Increase_Decrease",
                                                              "NA_NA")  
)

pdf("/Users/vincebilly/Desktop/vince/phd/statistic/afribiota_vince/Simple_risk_factor/figure/DESeq_summary_country_combined.pdf", height = 45, width = 20)
ggplot(log_melt3_mel, aes(x=variable, y=reorder(Cplt_taxo,-as.integer(factor(order)), FUN=min), size= value, color=direction)) +
  geom_point(alpha=0.5) +
  scale_color_manual(values=c("orange","blue"), name="Change in abundance") +
  theme(axis.text.x = element_text(angle = -45, hjust = 0))
dev.off()

to_keep <- c("Increase_Increase","Decrease_Decrease")
log_melt4_mel <- subset(log_melt3_mel,order %in% to_keep)
log_melt4_mel$Cplt_taxo2 <- gsub("ASV.*","",log_melt4_mel$Cplt_taxo)
log_melt4_mel$variable <- gsub("Both_Country_Non-Stunted","Non-Stunted",log_melt4_mel$variable)
log_melt4_mel$variable <- gsub("Both_Country_Stunted","Stunted",log_melt4_mel$variable)

#Create unique code for each taxa
code_table <- expand.grid(rep(list(letters[1:26]), 3))
code_taxo <- paste0("id",code_table$Var1,code_table$Var2,code_table$Var3)[1:length(log_melt4_mel$Cplt_taxo2)]
log_melt4_mel$code_taxo <- code_taxo

#Create Rank column
log_melt4_mel$Rank <- sub(".*\\|", "", log_melt4_mel$Cplt_taxo)

#Keep last name (remove everything before last pipe)
log_melt4_mel$Cplt_taxo2 <- gsub('.{1}$','',log_melt4_mel$Cplt_taxo2)
log_melt4_mel$Cplt_taxo2 <- sub(".*\\|", "", log_melt4_mel$Cplt_taxo2)

#Check for duplicate
length(log_melt4_mel$taxo)
length(unique(log_melt4_mel$taxo))

#
log_melt4_mel$taxo <- paste0(log_melt4_mel$Cplt_taxo2,"_",log_melt4_mel$code_taxo," (",log_melt4_mel$Rank,")")

pdf("/Users/vincebilly/Desktop/vince/phd/statistic/afribiota_vince/Simple_risk_factor/figure/DESeq_summary_country_combined_similar.pdf", height = 20, width = 20)
plot_deseq_sig_all <- ggplot(log_melt4_mel, aes(x=variable, y=reorder(taxo,value_in), size= value, color=direction)) +
  geom_point(alpha=0.5) +
  theme_light() +
  scale_color_manual(values=c("orange","blue"), name="Change in abundance") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=20, face="bold"),
        axis.text.y = element_text(size=15, face="bold"))
print(plot_deseq_sig_all)
dev.off()



final_table <- data.frame()
rank_vector <- c("Rank2","Rank3","Rank4","Rank5","Rank6","Rank7","Accession")
for(rnk in 1:length(rank_vector)){
  test1_rnk <- log_melt3_mel[grep(rank_vector[rnk],log_melt3_mel$Cplt_taxo),]
  test <- as.data.frame(table(test1_rnk$order))
  test$Sample <- rank_vector[rnk]
  final_table <- rbind(final_table,test)
}

final_table$Var1 <- factor(final_table$Var1, levels = c("Decrease_Increase",
                                          "Increase_Decrease",
                                          "Decrease_NA",
                                          "Decrease_Decrease",
                                          "NA_Decrease",
                                          "Increase_NA",
                                          "Increase_Increase",
                                          "NA_Increase",
                                          "NA_NA")  
)

final_table$Sample <- factor(final_table$Sample, levels = c("Rank2","Rank3","Rank4","Rank5","Rank6","Rank7","Accession"))
myCustomPalette <- c("pink","purple","orange","red","yellow","dodgerblue","blue","lightblue","grey90")
pdf("/Users/vincebilly/Desktop/vince/phd/statistic/afribiota_vince/Simple_risk_factor/figure/proportion_ASV_deseq_allrank.pdf", height = 3, width = 10)
ggplot(final_table, aes(x = Sample, y = Freq, fill = Var1)) + 
  theme_bw() +
  theme(strip.background = element_rect(fill="white")) + #these "theme" settings determine how the facet grid looks on the plot
  theme(strip.text = element_text(colour = 'black')) +
  geom_bar(stat = "identity", width = 1.0)  +
  scale_fill_manual(values = myCustomPalette) +
  coord_flip() +
  theme(legend.position="top",
        axis.text.y = element_text(size=20, face = "bold"),
        axis.text.x = element_text(size=20, face = "bold"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank())
dev.off()

####################################################################

####################################################################
###   Total Proportion of Taxa and reads that differ in DESeq   ###
####################################################################

full_abundance <- fread("/Users/vincebilly/Desktop/vince/phd/statistic/afribiota_vince/Simple_risk_factor/figure/Full_abundance_Table_Both_Country_Both_clinical.txt")
head(full_abundance)
full_abundance <- full_abundance[,c("Cplt_taxo")]
full_abundance <- unique(full_abundance)
full_abundance$Cplt_taxo <- gsub("[|]NA","",full_abundance$Cplt_taxo)

full_abundance1 <- full_abundance[grep("Rank",full_abundance$Cplt_taxo),]
full_abundance1$Cplt_taxo <- gsub("ASV.*[|]","",full_abundance1$Cplt_taxo)
full_abundance2 <- full_abundance[grep("Accession",full_abundance$Cplt_taxo),]
full_abundance <- rbind(full_abundance1,full_abundance2)


log_melt3_mel1 <- log_melt3_mel[grep("Rank",log_melt3_mel$Cplt_taxo),]
log_melt3_mel1$Cplt_taxo <- gsub("ASV.*[|]","",log_melt3_mel1$Cplt_taxo)
log_melt3_mel2 <- log_melt3_mel[grep("Accession",log_melt3_mel$Cplt_taxo),]
log_melt3_mel <- rbind(log_melt3_mel1,log_melt3_mel2)

intersect(full_abundance$Cplt_taxo, log_melt3_mel$Cplt_taxo)
setdiff(full_abundance$Cplt_taxo, log_melt3_mel$Cplt_taxo)
setdiff(log_melt3_mel$Cplt_taxo,full_abundance$Cplt_taxo)


data <- merge(full_abundance,log_melt3_mel, by="Cplt_taxo", all=TRUE)


final_table <- data.frame()
rank_vector <- c("Rank2","Rank3","Rank4","Rank5","Rank6","Rank7","Accession")
for(rnk in 1:length(rank_vector)){
  test1_rnk <- data[grep(rank_vector[rnk],data$Cplt_taxo),]
  test <- as.data.frame(table(test1_rnk$order))
  test$Sample <- rank_vector[rnk]
  final_table <- rbind(final_table,test)
}

final_table$Var1 <- factor(final_table$Var1, levels = c("Decrease_Increase",
                                                        "Increase_Decrease",
                                                        "Decrease_NA",
                                                        "Decrease_Decrease",
                                                        "NA_Decrease",
                                                        "Increase_NA",
                                                        "Increase_Increase",
                                                        "NA_Increase",
                                                        "NA_NA")  
)

final_table$Sample <- factor(final_table$Sample, levels = c("Rank2","Rank3","Rank4","Rank5","Rank6","Rank7","Accession"))
myCustomPalette <- c("pink","purple","orange","red","yellow","dodgerblue","blue","lightblue","grey90")
pdf("/Users/vincebilly/Desktop/vince/phd/statistic/afribiota_vince/Simple_risk_factor/figure/proportion_ASV_deseq_allrank.pdf", height = 3, width = 10)
ggplot(final_table, aes(x = Sample, y = Freq, fill = Var1)) + 
  theme_bw() +
  theme(strip.background = element_rect(fill="white")) + #these "theme" settings determine how the facet grid looks on the plot
  theme(strip.text = element_text(colour = 'black')) +
  geom_bar(stat = "identity", width = 1.0)  +
  scale_fill_manual(values = myCustomPalette) +
  coord_flip() +
  theme(legend.position="top",
        axis.text.y = element_text(size=20, face = "bold"),
        axis.text.x = element_text(size=20, face = "bold"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank())
dev.off()

####################################################################

####################################################################
###                  Bubble plot with DESeq results              ###
####################################################################

#
full_abundance <- fread("/Users/vincebilly/Desktop/vince/phd/statistic/afribiota_vince/Simple_risk_factor/figure/lefse_input/lefse_Both_Country_Both_clinical.txt")
#full_abundance <- fread("/Users/vincebilly/Desktop/vince/phd/statistic/afribiota_vince/Simple_risk_factor/figure/Full_abundance_Table_Both_Country_Both_clinical1.txt")
colnames(full_abundance) <- as.character(full_abundance[1,])
full_abundance <- full_abundance[-1,]

#Select only the ones that are significant
to_keep <- intersect(full_abundance$id,log_melt4_mel$Cplt_taxo2)
full_abundance <- subset(full_abundance, id %in% to_keep)
full_abundance[full_abundance==1] <- 0 
  
#Melt
full_abundance_melt <- melt(full_abundance, id.vars=c("id"));head(full_abundance_melt)


#Merge with metadata
qPCR_Metadata_sub <- qPCR_Metadata[,c("Sample","stunted","Cq_presence_absence")]
full_abundance_melt <- merge(qPCR_Metadata_sub, full_abundance_melt, by.x ="Sample", by.y = "variable", all=FALSE)

#Formating for plotting
full_abundance_melt$value <- as.numeric(full_abundance_melt$value)
full_abundance_melt$value2 <- ifelse(full_abundance_melt$value==0,NA,full_abundance_melt$value)
#ggplot
bubble_plot <- ggplot(full_abundance_melt, aes(x=Sample, y=id, size=value2, color=Cq_presence_absence)) +
  geom_point() +
  scale_color_manual(values=c("dodgerblue","gold")) +
  facet_grid(.~stunted + Cq_presence_absence, scales = "free", space = "free") +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size=3),
        strip.text.x = element_text(size=17))

#
heatmap_plot <- ggplot(full_abundance_melt, aes(x=Sample, y=id, fill=value2)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="black",trans = 'log') +
  facet_grid(.~stunted + Cq_presence_absence, scales = "free", space = "free") +
  theme_light() +
  theme(axis.text.y = element_blank())

full_abundance_melt <- full_abundance_melt %>% mutate(value3=case_when(Cq_presence_absence=="Presence_qPCR" ~ value,
                                                                        Cq_presence_absence=="Absence_qPCR"~ -value))
heatmap_plot <- ggplot(full_abundance_melt, aes(x=Sample, y=id, fill=value)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="dodgerblue",trans = 'log') +
  facet_grid(.~stunted + Cq_presence_absence, scales = "free", space = "free") +
  theme_light() +
  theme(axis.text.y = element_blank())

#Combine plot
pdf("/Users/vincebilly/Desktop/vince/phd/statistic/afribiota_vince/Simple_risk_factor/figure/DESeq_and_bubble.pdf", height = 15, width = 30)
plot_deseq_sig_all %>% insert_right(bubble_plot, width = 10)
dev.off()

pdf("/Users/vincebilly/Desktop/vince/phd/statistic/afribiota_vince/Simple_risk_factor/figure/DESeq_and_heatmap.pdf", height = 20, width = 30)
plot_deseq_sig_all %>% insert_right(heatmap_plot, width = 5)
dev.off()


####################################################################

getwd();dir()
lefse_Madagascar_Stunted <- read.csv("/Users/vincebilly/Desktop/lefse.csv", sep = "\t")
lefse_Madagascar_Non_Stunted <- read.csv("/Users/vincebilly/Desktop/lefse.csv", sep = "\t")
lefse_RCA_Stunted <- read.csv("/Users/vincebilly/Desktop/lefse.csv", sep = "\t")
lefse_RCA_Non_Stunted <- read.csv("/Users/vincebilly/Desktop/lefse.csv", sep = "\t")

#Convert into numeric values
lefse_Madagascar_Stunted$Pval <- as.numeric(lefse_Madagascar_Stunted$Pval)
lefse_Madagascar_Non_Stunted$Pval <- as.numeric(lefse_Madagascar_Non_Stunted$Pval)
lefse_RCA_Stunted$Pval <- as.numeric(lefse_RCA_Stunted$Pval)
lefse_RCA_Non_Stunted$Pval <- as.numeric(lefse_RCA_Non_Stunted$Pval)

#Select only significant results
lefse_Madagascar_Stunted <- subset(lefse_Madagascar_Stunted,direction %in% c("Stunted","Non-Stunted"))
lefse_Madagascar_Non_Stunted <- subset(lefse_Madagascar_Non_Stunted,direction %in% c("Stunted","Non-Stunted"))
lefse_RCA_Stunted <- subset(lefse_RCA_Stunted,direction %in% c("Stunted","Non-Stunted"))
lefse_RCA_Non_Stunted <- subset(lefse_RCA_Non_Stunted,direction %in% c("Stunted","Non-Stunted"))

#Convert into negative when depleted (stunted)
lefse_Madagascar_Stunted$value <- ifelse(lefse_Madagascar_Stunted$direction=="Non-Stunted",lefse_Madagascar_Stunted$value,-lefse_Madagascar_Stunted$value)
lefse_Madagascar_Non_Stunted$value <- ifelse(lefse_Madagascar_Non_Stunted$direction=="Non-Stunted",lefse_Madagascar_Non_Stunted$value,-lefse_Madagascar_Non_Stunted$value)
lefse_RCA_Stunted$value <- ifelse(lefse_RCA_Stunted$direction=="Non-Stunted",lefse_RCA_Stunted$value,-lefse_RCA_Stunted$value)
lefse_RCA_Non_Stunted$value <- ifelse(lefse_RCA_Non_Stunted$direction=="Non-Stunted",lefse_RCA_Non_Stunted$value,-lefse_RCA_Non_Stunted$value)

#Keep only column of taxa and log change
lefse_Madagascar_Stunted <- lefse_Madagascar_Stunted[,c(1,2)]
lefse_Madagascar_Non_Stunted <- lefse_Madagascar_Non_Stunted[,c(1,2)]
lefse_RCA_Stunted <- lefse_RCA_Stunted[,c(1,2)]
lefse_RCA_Non_Stunted <- lefse_RCA_Non_Stunted[,c(1,2)]

#Rename column names with stutned and country status
colnames(lefse_Madagascar_Stunted) <- c("Cplt_taxo","lefse_Madagascar_Stunted")
colnames(lefse_Madagascar_Non_Stunted) <- c("Cplt_taxo","lefse_Madagascar_Non_Stunted")
colnames(lefse_RCA_Stunted) <- c("Cplt_taxo","lefse_RCA_Stunted")
colnames(lefse_RCA_Non_Stunted) <- c("Cplt_taxo","lefse_RCA_Non_Stunted")

#Merge table
lefse_table_sig <- merge(lefse_Madagascar_Stunted,lefse_Madagascar_Non_Stunted, by="Cplt_taxo", all=TRUE)
lefse_table_sig <- merge(lefse_table_sig,lefse_RCA_Stunted, by="Cplt_taxo", all=TRUE)
lefse_table_sig <- merge(lefse_table_sig,lefse_RCA_Non_Stunted, by="Cplt_taxo", all=TRUE)

#
lefse_table_sig_melt <- melt(lefse_table_sig, id.vars=c("Cplt_taxo"),
                             measure.vars = c("lefse_Madagascar_Stunted","lefse_Madagascar_Non_Stunted","lefse_RCA_Stunted","lefse_RCA_Non_Stunted"));head(lefse_table_sig_melt)
lefse_table_sig_melt <- lefse_table_sig_melt %>% mutate(Cplt_taxo = str_replace_all(Cplt_taxo, '[|]$', ''))
lefse_table_sig_melt <- aggregate(value ~ Cplt_taxo + variable, data = lefse_table_sig_melt, FUN = mean, na.rm = TRUE)

#Add column with order
test <- merge(lefse_table_sig_melt,log_melt3_melt[,c(1,2)],by="Cplt_taxo", all.x=TRUE, all.y=FALSE)[,c(1,4,2,3)]

#Add column with test name
test$Test <- "LEfSe"
log_melt3_melt$Test <- "DESeq"
test_melt2 <- rbind(test, log_melt3_melt)
test_melt2$direction <- ifelse(test_melt2$value>0,"Increase","Deacrease")
test_melt2$value <- abs(test_melt2$value)

test_melt2$order <- as.character(test_melt2$order)

#test_melt2 <- subset(test_melt2, Subset !="NA")
test_melt2 <- test_melt2 %>%
  mutate_at(c('order'), ~replace_na(.,"NA"))
test_melt2$order <- factor(test_melt2$order, levels = c("Increase_Increase_Increase_Increase",
                                                              "Decrease_Decrease_Decrease_Decrease",
                                                              "Increase_NA_Increase_NA",
                                                              "Decrease_NA_Decrease_NA",
                                                              "NA_Increase_NA_Increase",
                                                              "NA_Decrease_NA_Decrease",
                                                              "Increase_Increase_NA_NA",
                                                              "Decrease_Decrease_NA_NA",
                                                              "NA_NA_Increase_Increase",
                                                              "NA_NA_Decrease_Decrease",
                                                              "Increase_Increase_Increase_NA",
                                                              "Decrease_Decrease_Decrease_NA",
                                                              "NA_Increase_Increase_Increase",
                                                              "NA_Decrease_Decrease_Decrease",
                                                              "Increase_Increase_NA_Increase",
                                                              "Increase_NA_Increase_Increase",
                                                              "Increase_NA_NA_Increase",
                                                              "Decrease_NA_NA_Decrease",
                                                              "NA_Increase_Increase_NA","NA_Decrease_Decrease_NA",
                                                              "NA_Increase_NA_NA","NA_Decrease_NA_NA",
                                                              "NA_NA_NA_Increase","NA_NA_NA_Decrease",
                                                              "NA_NA_Increase_NA","NA_NA_Decrease_NA",
                                                              "Increase_NA_NA_NA","Decrease_NA_NA_NA",
                                                              "NA_Decrease_NA_Increase","Increase_Decrease_NA_NA","NA_NA_Increase_Decrease","Increase_Decrease_Increase_NA","Decrease_NA_NA_Increase","NA_NA_Decrease_Increase","Increase_NA_Decrease_Decrease",
                                                              "Increase_NA_Increase_Decrease","NA_Decrease_Decrease_Increase","Increase_Decrease_Decrease_NA",
                                                              "NA_Increase_Increase_Decrease","Increase_NA_NA_Decrease","Increase_NA_Decrease_NA",
                                                              "Increase_Increase_NA_Decrease","Increase_Increase_Decrease_NA","NA_Decrease_Increase_NA",
                                                              "Decrease_Decrease_NA_Increase","NA_Increase_NA_Decrease",
                                                              "Decrease_NA_Increase_NA","Decrease_Decrease_NA_Decrease","Decrease_Increase_Increase_NA","Decrease_Increase_NA_NA",
                                                              "NA_NA_NA_NA","Increase_Decrease_Decrease_Decrease","Decrease_NA_Decrease_Increase","NA")  
)

levels(test_melt2$order)

pdf("/Users/vincebilly/Desktop/vince/phd/statistic/afribiota_vince/Simple_risk_factor/figure/DESeq_LEfSE_summary.pdf", height = 60, width = 20)
ggplot(test_melt2, aes(x=variable, y=reorder(Cplt_taxo,-as.integer(factor(order)), FUN=min), size= value, color=direction)) +
  geom_point(alpha=0.5) +
  scale_color_manual(values=c("orange","blue"), name="Change in abundance") +
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
  facet_grid( ~ Test, scales="free")
dev.off()







log_melt3$ASV <- as.character(row.names(log_melt3))




tax_tbl <- as.data.frame(tax_table(project_data))
tax_tbl$ASV <- row.names(tax_tbl)

merge_deseq <- as.data.frame(merge(log_melt3,tax_tbl, by = "ASV"))
merge_deseq$ASV_Rank6 <- paste0(merge_deseq$ASV,"_", merge_deseq$Rank6)
merge_deseq$variable <- "Presence vs Control"
merge_deseq <- merge_deseq[order(merge_deseq$order),]

log_deseq_bubble <- ggplot(merge_deseq, aes(x=variable, y=reorder(Rank6, -order), size=control_vs_ethanol, color = pos_neg)) +
  geom_point(alpha=0.5) +
  scale_color_manual(values=c("orange","purple"), name="Diversity change") +
  scale_size(range = c(.1, 7), name="Fold change in diversity") +
  theme(axis.text=element_text(size=8),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 10),
        axis.text.y = element_text(size = 15))+
  xlab("Deseq pair-wise comparison ") +
  ylab("Rank6")



head(log_melt)
names_column <- c("stunted.x","stunted.y")
for(i in 1:length(names_column)) {
  log_melt %>% mutate(Situation =case_when(names_column[i] > 0 ~ "Inc",
                                           names_column[i] < 0 ~ "Dea",
                                           names_column[i] == NA ~ "No"))
  colnames(log_melt)[dim(log_melt)[2]] <- paste0("Situation_",names_column)
}

#######################################################

####################################################################################################################################
####################################################################################################################################
###                                                       IndVal analysis                                                        ###
####################################################################################################################################
####################################################################################################################################

#######################################################
###         According to Siobhan script             ### 
#######################################################

theme_set(theme_bw()+
            theme(strip.background = element_rect(fill="white"),
                  axis.text.y = element_text(colour = "black", size = 12, face = "bold"),
                  axis.text.x = element_text(colour = "black", face = "bold", size = 12),
                  legend.text = element_text(size = 8, face ="bold", colour ="black"),
                  legend.position = "right", axis.title.y = element_text(face = "bold", size = 14),
                  axis.title.x = element_text(face = "bold", size = 14, colour = "black"),
                  legend.title = element_text(size = 14, colour = "black", face = "bold"),
                  legend.key=element_blank(),))
                  # axis.ticks = element_blank()

dephyloseq = function(subset_project_data){
  
  ## get the metadata
  meta = as.data.frame(as.matrix(subset_project_data@sam_data))
  
  ## how many metadata columns you have 
  metacols = ncol(meta)+1
  
  ## get out the otu table 
  otu = as.data.frame(as.matrix(subset_project_data@otu_table))
  
  ## merge the metadata and otu table by the rownames (sample ids from the Illumina sequencing data)
  mo = merge(meta, otu, by=0)
  
  ## get out the taxonomy file 
  tax = as.data.frame(subset_project_data@tax_table)
  
  ## get the ASV ID out. This the matches the placeholder ASV ID in the OTU table
  tax = tax %>% rownames_to_column(var="asv_id")
  
  ## pivot longer to be able to match the ASVs in the OTU table to the taxonomy table 
  mo = mo %>% pivot_longer(cols = -c(1:metacols), names_to = "asv_id", values_to="asv_abundance")
  
  ## Join the metadata and otu table with the taoxnomy table 
  mot = full_join(mo, tax)
  
  ## Specify the output for the dephyloseq funciton 
  output = mot
}

## set working directory and filepath 
setwd("/Users/vincebilly/Desktop")
indvalpath = "/Users/vincebilly/Desktop"
path = "/Users/vincebilly/Desktop"

## read in data
sf = subset_project_data

metadata = read.csv("metadata_2023_06_21.csv")
metadata = metadata %>% column_to_rownames(var="illumina_id")

nonrare = phyloseq(sample_data(metadata), otu_table(sf@otu_table), tax_table(sf@tax_table))

subset_project_data@sam_data$read_depth_filtered = sample_sums(nonrare)

field = subset_samples(nonrare, lab_or_field_sampled=="field")

# INDVAL SUBSTRATE CALCULATIONS at ASV LEVEL 
## get the metadata and otu data together for indval
metadata = as.data.frame(as.matrix(field@sam_data)) %>% rownames_to_column(var="rownames")
otu = as.data.frame(t(as.matrix(field@otu_table))) %>% rownames_to_column(var="rownames")
tax = as.data.frame(field@tax_table) %>% rownames_to_column(var="ASV")

## merge metadata and otu 
metaasv = left_join(metadata, otu)

## only keep salinity over or equal to 25 to find the core 
metaasv$mean.conduct = as.numeric(metaasv$mean.conduct)
saldf = subset(metaasv, metaasv$mean.conduct>=20)

## get metadata size
medatada.cols = ncol(metadata)

## Get otu table where samples are rows and ASVs are columns 
otu.sal = saldf[,-c(1:medatada.cols)]

## indval comparison 
substrate <- as.character(saldf$updated_sample_type) 

## run indval
indval <- multipatt(otu.sal, substrate, duleg = TRUE, control = how(nperm=999))

# Get indval statistic
indval.16.str <- as.data.frame(indval$str)
indval.16.str$rn <- rownames(indval.16.str)

# get p-value
indval.stat <- as.data.frame(indval$sign) #get dataframe of indval statistic
indval.stat$rn <- rownames(indval.stat) # make column of ASVs

# Prevalence as dataframe
indval.prev <- as.data.frame(indval$A)
# extract rownames into column
setDT(indval.prev, keep.rownames = TRUE)[] 
## rename columns 
colnames(indval.prev) <- paste0("prev.", colnames(indval.prev))
names(indval.prev)[names(indval.prev) == 'prev.rn'] <- 'rn'

# Fidelity as dataframe
indval.fid <- as.data.frame(indval$B) 
# extract rownames into column
setDT(indval.fid, keep.rownames = TRUE)[] 
## rename columns 
colnames(indval.fid) <- paste0("fid.", colnames(indval.fid))
names(indval.fid)[names(indval.fid) == 'fid.rn'] <- 'rn'

# Join statistics together (you could do multi_join but this might crash your computer)
str.and.stat = full_join(indval.16.str, indval.stat,
                         by="rn")
prev.and.fid = full_join(indval.prev, indval.fid,
                         by="rn")
indval_table = full_join(str.and.stat, prev.and.fid,
                         by="rn")

## rename columns to join with taxonomy
names(indval_table)[names(indval_table) == 'rn'] <- 'ASV'

## merge with taxonomy
indval_table= inner_join(indval_table, tax)

coresig = subset(indval_table, indval_table$index != "NaN" & indval_table$p.value<=0.05)




# FILTER INDVAL OUTPUT TO ONLY KEEP MERISTEM ASSOCIATED BY 0.8 ###

## filter indval output
sig = subset(coresig, coresig$s.meristem=="1" & 
               coresig$stat>=0.7)



# FORMAT THE SINGIFICANT ASSOCITIONS FOR ANALYSIS/PLOTTING ###

## get asv list by threshold
#m09 = c(unique(m09$asv_sequence))
m08 = c(unique(sig$asv_sequence))

## make dataframe for entire salinity dataset
meri = subset_samples(field, updated_sample_type=="meristem")
meri = dephyloseq(meri)

## use significant indicator lists to subset the meristem dataframe
#meri09 = subset(meri, meri$asv_sequence %in% c(m09))
meri08 = subset(meri, meri$asv_sequence %in% c(m08))

## get the dataframe of the core for the paper 
coredf = ddply(meri08, c("asv_id", "Order", "Genus", "asv_sequence"),
               summarise,
               ncore = length(year))
write.csv(coredf, paste0(indvalpath, "over20_core_asv.csv"))

## make presence/abscence column to sum when summarizing 
meri08$presabs = ifelse(meri08$asv_abundance>0, "1", "0")

# dotplot of EACH CORE ASV  by salinity ###
meriybasv = ddply(meri08, c("Row.names","mean.conduct", "mean.temp", "read_depth_filtered", "Genus", "asv_id"),
                  summarise,
                  corecount = sum(as.numeric(presabs)),
                  coresum = sum(asv_abundance)) %>%
  mutate(corera = as.numeric(coresum)/as.numeric(read_depth_filtered))


coregenus=ggplot(meriybasv, aes(x=as.numeric(mean.conduct), y = corera))+
  geom_point(cex=2, alpha=0.5)+
  geom_smooth(se=F, linewidth=1, method="lm")+
  labs(x="Salintiy", y = "")+
  #annotate("text", x =14, y = 0.2, label = "r = 0.2290, p < 0.001", size=4, color="black")+
  ylim(0,1)+
  facet_wrap(.~Genus)
coregenus

# CORR TEST FOR EACH ASV ###
ver = subset(meriybasv, meriybasv$asv_id=="ASV1")
cor.test(as.numeric(ver$mean.conduct), ver$corera)

rgam = subset(meriybasv, meriybasv$asv_id=="ASV15")
cor.test(as.numeric(gam$mean.conduct), gam$corera)

cau = subset(meriybasv, meriybasv$asv_id=="ASV3")
cor.test(as.numeric(cau$mean.conduct), cau$corera)

mar = subset(meriybasv, meriybasv$asv_id=="ASV4")
cor.test(as.numeric(mar$mean.conduct), mar$corera)

coc = subset(meriybasv, meriybasv$asv_id=="ASV5")
cor.test(as.numeric(coc$mean.conduct), coc$corera)

rob = subset(meriybasv, meriybasv$asv_id=="ASV7")
cor.test(as.numeric(rob$mean.conduct), rob$corera)

# summarise the data to count the relative abundance of core per ###

meri08 = ddply(meri08, c("Row.names","mean.conduct", "mean.temp", "read_depth_filtered"),
               summarise,
               corecount = sum(as.numeric(presabs)),
               coresum = sum(asv_abundance)) %>%
  mutate(corera = as.numeric(coresum)/as.numeric(read_depth_filtered))


# CORRELATION TEST OF CORE BY SALINITY ##
cor.test(as.numeric(meri08$mean.conduct), meri08$corera)


# dotplot of ALL CORE relative abundance by salinity ##
corera=ggplot(meri08, aes(x=as.numeric(mean.conduct), y = corera))+
  geom_point(cex=2, alpha=0.5)+
  geom_smooth(se=F, linewidth=1, method="lm")+
  labs(x="Salintiy", y = "Relative abundance of core ASV
as percent of total reads")+
  annotate("text", x =14, y = 0.8, label = "r = 0.2290, p < 0.001", size=4, color="black")+
  ylim(0,1)





# ARRANGE PLOTS ##
ggarrange(corera,coregenus, 
          labels=c('A', "B"),
          widths=c(0.7, 1))

ggsave("field_core_ra_count_over20.pdf", path=path, width = 17.9, height=5.9, units="in")

#######################################################




####################################################################################################################################
####################################################################################################################################
### LefSe
####################################################################################################################################
####################################################################################################################################

## Stunted- no-stunted all samples
##################################################################

subset_project_data2 <- subset_project_data


lefse_table <- data.frame()                                                     #Create dataframe
rank <- unique(as.character(colnames(tax_table(subset_project_data2))))              #Create vector for all ranks
sample_df <- as.data.frame(as.matrix(sample_data(subset_project_data2)))

#Select for rank
for (r in 1:length(rank)) {
  project_data_r <- tax_glom(subset_project_data2, taxrank=rank[r])                  #Aggregation at each rank individually (from 1 to 7)
  otu <- as.data.frame(t(otu_table(project_data_r)))                       #Extract otu table
  tax <- as.data.frame(as(tax_table(project_data_r), "matrix"))                #Extract tax table
  tax_otu <- merge(tax,otu, by ="row.names", all = T)                           #Merge tax and otu table
  lefse_table <- rbind(lefse_table,tax_otu)                                     #Concatenate result from each Rank
}
#Combine taxa name in one column separated with pipe "|"
lefse_table$Row.names <- paste0(lefse_table$Rank1,"|",lefse_table$Rank2,"|",lefse_table$Rank3,"|",lefse_table$Rank4,"|",lefse_table$Rank5,"|",lefse_table$Rank6,"|",lefse_table$Rank7,"|",lefse_table$Rank8)
#Remove NA in the same column
lefse_table$Row.names <- gsub("[|]NA","",lefse_table$Row.names)
#Remove all "Rank" columns
lefse_table <- lefse_table[,setdiff(colnames(lefse_table),colnames(tax_table(subset_project_data2)))]
#Reformate table for lefse analysis input
colnames(lefse_table)[1] <- "id"
id <- colnames(lefse_table)
Treatment <- c("Treatment", as.character(sample_df$stunted))
lefse_table <- rbind(id,lefse_table)
colnames(lefse_table) <- Treatment
#Save table
setwd("/Users/vincebilly/Desktop/vince/phd/statistic/afribiota_vince/Simple_risk_factor/figure/lefse_input/")
write.table(lefse_table, paste0("lefse_stunted_nonstunted_all",".txt"), sep="\t", row.names = FALSE, quote=FALSE)

##################################################################


## all
##################################################################

subset_project_data2 <- subset_project_data


lefse_table <- data.frame()                                                     #Create dataframe
rank <- unique(as.character(colnames(tax_table(subset_project_data2))))              #Create vector for all ranks
sample_df <- as.data.frame(as.matrix(sample_data(subset_project_data2)))

#Select for rank
for (r in 1:length(rank)) {
  project_data_r <- tax_glom(subset_project_data2, taxrank=rank[r])                  #Aggregation at each rank individually (from 1 to 7)
  otu <- as.data.frame(t(otu_table(project_data_r)))                       #Extract otu table
  tax <- as.data.frame(as(tax_table(project_data_r), "matrix"))                #Extract tax table
  tax_otu <- merge(tax,otu, by ="row.names", all = T)                           #Merge tax and otu table
  lefse_table <- rbind(lefse_table,tax_otu)                                     #Concatenate result from each Rank
}
#Combine taxa name in one column separated with pipe "|"
lefse_table$Row.names <- paste0(lefse_table$Rank1,"|",lefse_table$Rank2,"|",lefse_table$Rank3,"|",lefse_table$Rank4,"|",lefse_table$Rank5,"|",lefse_table$Rank6,"|",lefse_table$Rank7,"|",lefse_table$Rank8)
#Remove NA in the same column
lefse_table$Row.names <- gsub("[|]NA","",lefse_table$Row.names)
#Remove all "Rank" columns
lefse_table <- lefse_table[,setdiff(colnames(lefse_table),colnames(tax_table(subset_project_data2)))]
#Reformate table for lefse analysis input
colnames(lefse_table)[1] <- "id"
id <- colnames(lefse_table)
Treatment <- c("Treatment", as.character(sample_df$Cq_presence_absence))
lefse_table <- rbind(id,lefse_table)
colnames(lefse_table) <- Treatment
#Save table
setwd("/Users/vincebilly/Desktop/vince/phd/statistic/afribiota_vince/Simple_risk_factor/figure/lefse_input/")
write.table(lefse_table, paste0("lefse_Cq_presence_absence_all",".txt"), sep="\t", row.names = FALSE, quote=FALSE)

##################################################################


## Madagascar
##################################################################

pays_vector <- c("Madagascar","RCA","")
subset_project_data2 <- subset_project_data %>%
  subset_samples(pays == pays_vector[i])


lefse_table <- data.frame()                                                     #Create dataframe
rank <- unique(as.character(colnames(tax_table(subset_project_data2))))              #Create vector for all ranks
sample_df <- as.data.frame(as.matrix(sample_data(subset_project_data2)))

#Select for rank
for (r in 1:length(rank)) {
  project_data_r <- tax_glom(subset_project_data2, taxrank=rank[r])                  #Aggregation at each rank individually (from 1 to 7)
  otu <- as.data.frame(t(otu_table(project_data_r)))                       #Extract otu table
  tax <- as.data.frame(as(tax_table(project_data_r), "matrix"))                #Extract tax table
  tax_otu <- merge(tax,otu, by ="row.names", all = T)                           #Merge tax and otu table
  lefse_table <- rbind(lefse_table,tax_otu)                                     #Concatenate result from each Rank
}
#Combine taxa name in one column separated with pipe "|"
lefse_table$Row.names <- paste0(lefse_table$Rank1,"|",lefse_table$Rank2,"|",lefse_table$Rank3,"|",lefse_table$Rank4,"|",lefse_table$Rank5,"|",lefse_table$Rank6,"|",lefse_table$Rank7,"|",lefse_table$Rank8)
#Remove NA in the same column
lefse_table$Row.names <- gsub("[|]NA","",lefse_table$Row.names)
#Remove all "Rank" columns
lefse_table <- lefse_table[,setdiff(colnames(lefse_table),colnames(tax_table(subset_project_data2)))]
#Reformate table for lefse analysis input
colnames(lefse_table)[1] <- "id"
id <- colnames(lefse_table)
Treatment <- c("Treatment", as.character(sample_df$Cq_presence_absence))
lefse_table <- rbind(id,lefse_table)
colnames(lefse_table) <- Treatment
#Save table
setwd("/Users/vincebilly/Desktop/vince/phd/statistic/afribiota_vince/Simple_risk_factor/figure/lefse_input/")
write.table(lefse_table, paste0("lefse_Cq_presence_absence_Mada",".txt"), sep="\t", row.names = FALSE, quote=FALSE)

##################################################################


#Mada non stunted
##################################################################

subset_project_data2 <- subset_project_data %>%
  subset_samples(pays == "Madagascar") %>%
  subset_samples(stunted == "Non-Stunted")


lefse_table <- data.frame()                                                     #Create dataframe
rank <- unique(as.character(colnames(tax_table(subset_project_data2))))              #Create vector for all ranks
sample_df <- as.data.frame(as.matrix(sample_data(subset_project_data2)))

#Select for rank
for (r in 1:length(rank)) {
  project_data_r <- tax_glom(subset_project_data2, taxrank=rank[r])                  #Aggregation at each rank individually (from 1 to 7)
  otu <- as.data.frame(t(otu_table(project_data_r)))                       #Extract otu table
  tax <- as.data.frame(as(tax_table(project_data_r), "matrix"))                #Extract tax table
  tax_otu <- merge(tax,otu, by ="row.names", all = T)                           #Merge tax and otu table
  lefse_table <- rbind(lefse_table,tax_otu)                                     #Concatenate result from each Rank
}
#Combine taxa name in one column separated with pipe "|"
lefse_table$Row.names <- paste0(lefse_table$Rank1,"|",lefse_table$Rank2,"|",lefse_table$Rank3,"|",lefse_table$Rank4,"|",lefse_table$Rank5,"|",lefse_table$Rank6,"|",lefse_table$Rank7,"|",lefse_table$Rank8)
#Remove NA in the same column
lefse_table$Row.names <- gsub("[|]NA","",lefse_table$Row.names)
#Remove all "Rank" columns
lefse_table <- lefse_table[,setdiff(colnames(lefse_table),colnames(tax_table(subset_project_data2)))]
#Reformate table for lefse analysis input
colnames(lefse_table)[1] <- "id"
id <- colnames(lefse_table)
Treatment <- c("Treatment", as.character(sample_df$Cq_presence_absence))
lefse_table <- rbind(id,lefse_table)
colnames(lefse_table) <- Treatment
#Save table
setwd("/Users/vincebilly/Desktop/vince/phd/statistic/afribiota_vince/Simple_risk_factor/figure/lefse_input/")
write.table(lefse_table, paste0("lefse_Cq_presence_absence_Mada_Non-stunted",".txt"), sep="\t", row.names = FALSE, quote=FALSE)

##################################################################

#Mada stunted
##################################################################

subset_project_data2 <- subset_project_data %>%
  subset_samples(pays == "Madagascar") %>%
  subset_samples(stunted == "Stunted")


lefse_table <- data.frame()                                                     #Create dataframe
rank <- unique(as.character(colnames(tax_table(subset_project_data2))))              #Create vector for all ranks
sample_df <- as.data.frame(as.matrix(sample_data(subset_project_data2)))

#Select for rank
for (r in 1:length(rank)) {
  project_data_r <- tax_glom(subset_project_data2, taxrank=rank[r])                  #Aggregation at each rank individually (from 1 to 7)
  otu <- as.data.frame(t(otu_table(project_data_r)))                       #Extract otu table
  tax <- as.data.frame(as(tax_table(project_data_r), "matrix"))                #Extract tax table
  tax_otu <- merge(tax,otu, by ="row.names", all = T)                           #Merge tax and otu table
  lefse_table <- rbind(lefse_table,tax_otu)                                     #Concatenate result from each Rank
}
#Combine taxa name in one column separated with pipe "|"
lefse_table$Row.names <- paste0(lefse_table$Rank1,"|",lefse_table$Rank2,"|",lefse_table$Rank3,"|",lefse_table$Rank4,"|",lefse_table$Rank5,"|",lefse_table$Rank6,"|",lefse_table$Rank7,"|",lefse_table$Rank8)
#Remove NA in the same column
lefse_table$Row.names <- gsub("[|]NA","",lefse_table$Row.names)
#Remove all "Rank" columns
lefse_table <- lefse_table[,setdiff(colnames(lefse_table),colnames(tax_table(subset_project_data2)))]
#Reformate table for lefse analysis input
colnames(lefse_table)[1] <- "id"
id <- colnames(lefse_table)
Treatment <- c("Treatment", as.character(sample_df$Cq_presence_absence))
lefse_table <- rbind(id,lefse_table)
colnames(lefse_table) <- Treatment
#Save table
setwd("/Users/vincebilly/Desktop/vince/phd/statistic/afribiota_vince/Simple_risk_factor/figure/lefse_input/")
write.table(lefse_table, paste0("lefse_Cq_presence_absence_Mada_stunted",".txt"), sep="\t", row.names = FALSE, quote=FALSE)

##################################################################

#RCA non stunted
##################################################################

subset_project_data2 <- subset_project_data %>%
  subset_samples(pays == "RCA") %>%
  subset_samples(stunted == "Non-Stunted")


lefse_table <- data.frame()                                                     #Create dataframe
rank <- unique(as.character(colnames(tax_table(subset_project_data2))))              #Create vector for all ranks
sample_df <- as.data.frame(as.matrix(sample_data(subset_project_data2)))

#Select for rank
for (r in 1:length(rank)) {
  project_data_r <- tax_glom(subset_project_data2, taxrank=rank[r])                  #Aggregation at each rank individually (from 1 to 7)
  otu <- as.data.frame(t(otu_table(project_data_r)))                       #Extract otu table
  tax <- as.data.frame(as(tax_table(project_data_r), "matrix"))                #Extract tax table
  tax_otu <- merge(tax,otu, by ="row.names", all = T)                           #Merge tax and otu table
  lefse_table <- rbind(lefse_table,tax_otu)                                     #Concatenate result from each Rank
}
#Combine taxa name in one column separated with pipe "|"
lefse_table$Row.names <- paste0(lefse_table$Rank1,"|",lefse_table$Rank2,"|",lefse_table$Rank3,"|",lefse_table$Rank4,"|",lefse_table$Rank5,"|",lefse_table$Rank6,"|",lefse_table$Rank7,"|",lefse_table$Rank8)
#Remove NA in the same column
lefse_table$Row.names <- gsub("[|]NA","",lefse_table$Row.names)
#Remove all "Rank" columns
lefse_table <- lefse_table[,setdiff(colnames(lefse_table),colnames(tax_table(subset_project_data2)))]
#Reformate table for lefse analysis input
colnames(lefse_table)[1] <- "id"
id <- colnames(lefse_table)
Treatment <- c("Treatment", as.character(sample_df$Cq_presence_absence))
lefse_table <- rbind(id,lefse_table)
colnames(lefse_table) <- Treatment
#Save table
setwd("/Users/vincebilly/Desktop/vince/phd/statistic/afribiota_vince/Simple_risk_factor/figure/lefse_input/")
write.table(lefse_table, paste0("lefse_Cq_presence_absence_RCA_Non-stunted",".txt"), sep="\t", row.names = FALSE, quote=FALSE)

##################################################################

#RCA stunted
##################################################################

subset_project_data2 <- subset_project_data %>%
  subset_samples(pays == "RCA") %>%
  subset_samples(stunted == "Stunted")

lefse_table <- data.frame()                                                     #Create dataframe
rank <- unique(as.character(colnames(tax_table(subset_project_data2))))              #Create vector for all ranks
sample_df <- as.data.frame(as.matrix(sample_data(subset_project_data2)))

#Select for rank
for (r in 1:length(rank)) {
  project_data_r <- tax_glom(subset_project_data2, taxrank=rank[r])                  #Aggregation at each rank individually (from 1 to 7)
  otu <- as.data.frame(t(otu_table(project_data_r)))                       #Extract otu table
  tax <- as.data.frame(as(tax_table(project_data_r), "matrix"))                #Extract tax table
  tax_otu <- merge(tax,otu, by ="row.names", all = T)                           #Merge tax and otu table
  lefse_table <- rbind(lefse_table,tax_otu)                                     #Concatenate result from each Rank
}
#Combine taxa name in one column separated with pipe "|"
lefse_table$Row.names <- paste0(lefse_table$Rank1,"|",lefse_table$Rank2,"|",lefse_table$Rank3,"|",lefse_table$Rank4,"|",lefse_table$Rank5,"|",lefse_table$Rank6,"|",lefse_table$Rank7,"|",lefse_table$Rank8)
#Remove NA in the same column
lefse_table$Row.names <- gsub("[|]NA","",lefse_table$Row.names)
#Remove all "Rank" columns
lefse_table <- lefse_table[,setdiff(colnames(lefse_table),colnames(tax_table(subset_project_data2)))]
#Reformate table for lefse analysis input
colnames(lefse_table)[1] <- "id"
id <- colnames(lefse_table)
Treatment <- c("Treatment", as.character(sample_df$Cq_presence_absence))
lefse_table <- rbind(id,lefse_table)
colnames(lefse_table) <- Treatment
#Save table
setwd("/Users/vincebilly/Desktop/vince/phd/statistic/afribiota_vince/Simple_risk_factor/figure/lefse_input/")
write.table(lefse_table, paste0("lefse_Cq_presence_absence_RCA_stunted",".txt"), sep="\t", row.names = FALSE, quote=FALSE)

##################################################################


####################################################################################################################################
####################################################################################################################################
### LefSe
####################################################################################################################################
####################################################################################################################################

## Madagascar
##################################################################

subset_project_data2 <- subset_project_data


lefse_table <- data.frame()                                                     #Create dataframe
rank <- unique(as.character(colnames(tax_table(subset_project_data2))))              #Create vector for all ranks
sample_df <- as.data.frame(as.matrix(sample_data(subset_project_data2)))

#Select for rank
for (r in 1:length(rank)) {
  project_data_r <- tax_glom(subset_project_data2, taxrank=rank[r])                  #Aggregation at each rank individually (from 1 to 7)
  otu <- as.data.frame(otu_table(project_data_r))                       #Extract otu table
  tax <- as.data.frame(as(tax_table(project_data_r), "matrix"))                #Extract tax table
  tax_otu <- merge(tax,otu, by ="row.names", all = T)                           #Merge tax and otu table
  lefse_table <- rbind(lefse_table,tax_otu)                                     #Concatenate result from each Rank
}
#Combine taxa name in one column separated with pipe "|"
lefse_table$Row.names <- paste0(lefse_table$Rank1,"|",lefse_table$Rank2,"|",lefse_table$Rank3,"|",lefse_table$Rank4,"|",lefse_table$Rank5,"|",lefse_table$Rank6,"|",lefse_table$Rank7,"|",lefse_table$Rank8)
#Remove NA in the same column
lefse_table$Row.names <- gsub("[|]NA","",lefse_table$Row.names)
#Remove all "Rank" columns
lefse_table <- lefse_table[,setdiff(colnames(lefse_table),colnames(tax_table(subset_project_data2)))]
#Reformate table for lefse analysis input
colnames(lefse_table)[1] <- "id"
id <- colnames(lefse_table)
Treatment <- c("Treatment", as.character(sample_df$Cq_presence_absence))
lefse_table <- rbind(id,lefse_table)
colnames(lefse_table) <- Treatment
#Save table
setwd("/Users/vincebilly/Desktop/vince/phd/statistic/afribiota_vince/Simple_risk_factor/figure/lefse_input/18s")
write.table(lefse_table, paste0("lefse_Cq_presence_absence",".txt"), sep="\t", row.names = FALSE, quote=FALSE)

##################################################################




## Madagascar
##################################################################

subset_project_data2 <- subset_project_data %>%
  subset_samples(pays == "Madagascar")


lefse_table <- data.frame()                                                     #Create dataframe
rank <- unique(as.character(colnames(tax_table(subset_project_data2))))              #Create vector for all ranks
sample_df <- as.data.frame(as.matrix(sample_data(subset_project_data2)))

#Select for rank
for (r in 1:length(rank)) {
  project_data_r <- tax_glom(subset_project_data2, taxrank=rank[r])                  #Aggregation at each rank individually (from 1 to 7)
  otu <- as.data.frame(t(otu_table(project_data_r)))                       #Extract otu table
  tax <- as.data.frame(as(tax_table(project_data_r), "matrix"))                #Extract tax table
  tax_otu <- merge(tax,otu, by ="row.names", all = T)                           #Merge tax and otu table
  lefse_table <- rbind(lefse_table,tax_otu)                                     #Concatenate result from each Rank
}
#Combine taxa name in one column separated with pipe "|"
lefse_table$Row.names <- paste0(lefse_table$Rank1,"|",lefse_table$Rank2,"|",lefse_table$Rank3,"|",lefse_table$Rank4,"|",lefse_table$Rank5,"|",lefse_table$Rank6,"|",lefse_table$Rank7,"|",lefse_table$Rank8)
#Remove NA in the same column
lefse_table$Row.names <- gsub("[|]NA","",lefse_table$Row.names)
#Remove all "Rank" columns
lefse_table <- lefse_table[,setdiff(colnames(lefse_table),colnames(tax_table(subset_project_data2)))]
#Reformate table for lefse analysis input
colnames(lefse_table)[1] <- "id"
id <- colnames(lefse_table)
Treatment <- c("Treatment", as.character(sample_df$Cq_presence_absence))
lefse_table <- rbind(id,lefse_table)
colnames(lefse_table) <- Treatment
#Save table
setwd("/Users/vincebilly/Desktop/vince/phd/statistic/afribiota_vince/Simple_risk_factor/figure/lefse_input/")
write.table(lefse_table, paste0("lefse_Cq_presence_absence_Mada",".txt"), sep="\t", row.names = FALSE, quote=FALSE)

##################################################################


#Mada non stunted
##################################################################

subset_project_data2 <- subset_project_data %>%
  subset_samples(pays == "Madagascar") %>%
  subset_samples(stunted == "Non")


lefse_table <- data.frame()                                                     #Create dataframe
rank <- unique(as.character(colnames(tax_table(subset_project_data2))))              #Create vector for all ranks
sample_df <- as.data.frame(as.matrix(sample_data(subset_project_data2)))

#Select for rank
for (r in 1:length(rank)) {
  project_data_r <- tax_glom(subset_project_data2, taxrank=rank[r])                  #Aggregation at each rank individually (from 1 to 7)
  otu <- as.data.frame(t(otu_table(project_data_r)))                       #Extract otu table
  tax <- as.data.frame(as(tax_table(project_data_r), "matrix"))                #Extract tax table
  tax_otu <- merge(tax,otu, by ="row.names", all = T)                           #Merge tax and otu table
  lefse_table <- rbind(lefse_table,tax_otu)                                     #Concatenate result from each Rank
}
#Combine taxa name in one column separated with pipe "|"
lefse_table$Row.names <- paste0(lefse_table$Rank1,"|",lefse_table$Rank2,"|",lefse_table$Rank3,"|",lefse_table$Rank4,"|",lefse_table$Rank5,"|",lefse_table$Rank6,"|",lefse_table$Rank7,"|",lefse_table$Rank8)
#Remove NA in the same column
lefse_table$Row.names <- gsub("[|]NA","",lefse_table$Row.names)
#Remove all "Rank" columns
lefse_table <- lefse_table[,setdiff(colnames(lefse_table),colnames(tax_table(subset_project_data2)))]
#Reformate table for lefse analysis input
colnames(lefse_table)[1] <- "id"
id <- colnames(lefse_table)
Treatment <- c("Treatment", as.character(sample_df$Cq_presence_absence))
lefse_table <- rbind(id,lefse_table)
colnames(lefse_table) <- Treatment
#Save table
setwd("/Users/vincebilly/Desktop/vince/phd/statistic/afribiota_vince/Simple_risk_factor/figure/lefse_input/")
write.table(lefse_table, paste0("lefse_Cq_presence_absence_Mada_Non-stunted",".txt"), sep="\t", row.names = FALSE, quote=FALSE)

##################################################################

#Mada stunted
##################################################################

subset_project_data2 <- subset_project_data %>%
  subset_samples(pays == "Madagascar") %>%
  subset_samples(stunted == "Oui")


lefse_table <- data.frame()                                                     #Create dataframe
rank <- unique(as.character(colnames(tax_table(subset_project_data2))))              #Create vector for all ranks
sample_df <- as.data.frame(as.matrix(sample_data(subset_project_data2)))

#Select for rank
for (r in 1:length(rank)) {
  project_data_r <- tax_glom(subset_project_data2, taxrank=rank[r])                  #Aggregation at each rank individually (from 1 to 7)
  otu <- as.data.frame(t(otu_table(project_data_r)))                       #Extract otu table
  tax <- as.data.frame(as(tax_table(project_data_r), "matrix"))                #Extract tax table
  tax_otu <- merge(tax,otu, by ="row.names", all = T)                           #Merge tax and otu table
  lefse_table <- rbind(lefse_table,tax_otu)                                     #Concatenate result from each Rank
}
#Combine taxa name in one column separated with pipe "|"
lefse_table$Row.names <- paste0(lefse_table$Rank1,"|",lefse_table$Rank2,"|",lefse_table$Rank3,"|",lefse_table$Rank4,"|",lefse_table$Rank5,"|",lefse_table$Rank6,"|",lefse_table$Rank7,"|",lefse_table$Rank8)
#Remove NA in the same column
lefse_table$Row.names <- gsub("[|]NA","",lefse_table$Row.names)
#Remove all "Rank" columns
lefse_table <- lefse_table[,setdiff(colnames(lefse_table),colnames(tax_table(subset_project_data2)))]
#Reformate table for lefse analysis input
colnames(lefse_table)[1] <- "id"
id <- colnames(lefse_table)
Treatment <- c("Treatment", as.character(sample_df$Cq_presence_absence))
lefse_table <- rbind(id,lefse_table)
colnames(lefse_table) <- Treatment
#Save table
setwd("/Users/vincebilly/Desktop/vince/phd/statistic/afribiota_vince/Simple_risk_factor/figure/lefse_input/")
write.table(lefse_table, paste0("lefse_Cq_presence_absence_Mada_stunted",".txt"), sep="\t", row.names = FALSE, quote=FALSE)

##################################################################

#RCA non stunted
##################################################################

subset_project_data2 <- subset_project_data %>%
  subset_samples(pays == "RCA") %>%
  subset_samples(stunted == "Non")


lefse_table <- data.frame()                                                     #Create dataframe
rank <- unique(as.character(colnames(tax_table(subset_project_data2))))              #Create vector for all ranks
sample_df <- as.data.frame(as.matrix(sample_data(subset_project_data2)))

#Select for rank
for (r in 1:length(rank)) {
  project_data_r <- tax_glom(subset_project_data2, taxrank=rank[r])                  #Aggregation at each rank individually (from 1 to 7)
  otu <- as.data.frame(t(otu_table(project_data_r)))                       #Extract otu table
  tax <- as.data.frame(as(tax_table(project_data_r), "matrix"))                #Extract tax table
  tax_otu <- merge(tax,otu, by ="row.names", all = T)                           #Merge tax and otu table
  lefse_table <- rbind(lefse_table,tax_otu)                                     #Concatenate result from each Rank
}
#Combine taxa name in one column separated with pipe "|"
lefse_table$Row.names <- paste0(lefse_table$Rank1,"|",lefse_table$Rank2,"|",lefse_table$Rank3,"|",lefse_table$Rank4,"|",lefse_table$Rank5,"|",lefse_table$Rank6,"|",lefse_table$Rank7,"|",lefse_table$Rank8)
#Remove NA in the same column
lefse_table$Row.names <- gsub("[|]NA","",lefse_table$Row.names)
#Remove all "Rank" columns
lefse_table <- lefse_table[,setdiff(colnames(lefse_table),colnames(tax_table(subset_project_data2)))]
#Reformate table for lefse analysis input
colnames(lefse_table)[1] <- "id"
id <- colnames(lefse_table)
Treatment <- c("Treatment", as.character(sample_df$Cq_presence_absence))
lefse_table <- rbind(id,lefse_table)
colnames(lefse_table) <- Treatment
#Save table
setwd("/Users/vincebilly/Desktop/vince/phd/statistic/afribiota_vince/Simple_risk_factor/figure/lefse_input/")
write.table(lefse_table, paste0("lefse_Cq_presence_absence_RCA_Non-stunted",".txt"), sep="\t", row.names = FALSE, quote=FALSE)

##################################################################

#RCA stunted
##################################################################

subset_project_data2 <- subset_project_data %>%
  subset_samples(pays == "RCA") %>%
  subset_samples(stunted == "Oui")


lefse_table <- data.frame()                                                     #Create dataframe
rank <- unique(as.character(colnames(tax_table(subset_project_data2))))              #Create vector for all ranks
sample_df <- as.data.frame(as.matrix(sample_data(subset_project_data2)))

#Select for rank
for (r in 1:length(rank)) {
  project_data_r <- tax_glom(subset_project_data2, taxrank=rank[r])                  #Aggregation at each rank individually (from 1 to 7)
  otu <- as.data.frame(t(otu_table(project_data_r)))                       #Extract otu table
  tax <- as.data.frame(as(tax_table(project_data_r), "matrix"))                #Extract tax table
  tax_otu <- merge(tax,otu, by ="row.names", all = T)                           #Merge tax and otu table
  lefse_table <- rbind(lefse_table,tax_otu)                                     #Concatenate result from each Rank
}
#Combine taxa name in one column separated with pipe "|"
lefse_table$Row.names <- paste0(lefse_table$Rank1,"|",lefse_table$Rank2,"|",lefse_table$Rank3,"|",lefse_table$Rank4,"|",lefse_table$Rank5,"|",lefse_table$Rank6,"|",lefse_table$Rank7,"|",lefse_table$Rank8)
#Remove NA in the same column
lefse_table$Row.names <- gsub("[|]NA","",lefse_table$Row.names)
#Remove all "Rank" columns
lefse_table <- lefse_table[,setdiff(colnames(lefse_table),colnames(tax_table(subset_project_data2)))]
#Reformate table for lefse analysis input
colnames(lefse_table)[1] <- "id"
id <- colnames(lefse_table)
Treatment <- c("Treatment", as.character(sample_df$Cq_presence_absence))
lefse_table <- rbind(id,lefse_table)
colnames(lefse_table) <- Treatment
#Save table
setwd("/Users/vincebilly/Desktop/vince/phd/statistic/afribiota_vince/Simple_risk_factor/figure/lefse_input/")
write.table(lefse_table, paste0("lefse_Cq_presence_absence_RCA_stunted",".txt"), sep="\t", row.names = FALSE, quote=FALSE)

##################################################################

####################################################################################################################################
####################################################################################################################################
###                                         Comparison Blastocystis and Metadata                                                 ###
####################################################################################################################################
####################################################################################################################################

##################################################################
###               COMBINE qPCR and METADATA                    ###
##################################################################

#Check if samples match with MISEQ data
qPCR_Metadata <- merge(plate3,metadata, by.x = "Sample", by.y = "id", all= FALSE);dim(qPCR_Metadata)


qPCR_Metadata$Cq <- as.numeric(qPCR_Metadata$Cq)
qPCR_Metadata$Cq_presence_absence <- ifelse(qPCR_Metadata$Cq<=35 & qPCR_Metadata$Cq>10,"Presence_qPCR", "Absence_qPCR")

#Combined qPCR and metadata
qPCR_Metadata$Cq_presence_absence <- as.factor(qPCR_Metadata$Cq_presence_absence)

#Reanotate variable
qPCR_Metadata$stunted <- gsub("Oui","Stunted",qPCR_Metadata$stunted)
qPCR_Metadata$stunted <- gsub("Non","Non-Stunted",qPCR_Metadata$stunted)


##################################################################


##################################################################
## Linear regression Abundance Blastocystis and haz
##################################################################

qPCR_Metadata_presence <- subset(qPCR_Metadata, Cq_presence_absence =="Presence_qPCR")

#Shapiro test

shapiro.test(qPCR_Metadata_presence$Cq)
hist(qPCR_Metadata_presence$Cq)
#Not normal distributed

qPCR_Metadata_presence$haz_cont <- as.numeric(qPCR_Metadata_presence$haz_cont)
ggplot(data = qPCR_Metadata_presence, aes(x = Cq, y = haz_cont)) +
  stat_poly_line(color="darkgreen") +
  stat_poly_eq(use_label(c("adj.R2" ,"p", "n"))) +
  geom_point(size=1, color="darkgreen") +
  facet_wrap(~pays+stunted, ncol=8) + 
  theme_classic() +
  scale_x_reverse()

ggplot(data = qPCR_Metadata_presence, aes(x = Cq, y = haz_cont)) +
  stat_poly_line(color="darkgreen") +
  stat_poly_eq(use_label(c("adj.R2" ,"p", "n"))) +
  geom_point(size=1, color="darkgreen") +
  facet_wrap(~pays, ncol=8) + 
  theme_classic() +
  scale_x_reverse()

ggplot(data = qPCR_Metadata_presence, aes(x = Cq, y = haz_cont)) +
  stat_poly_line(color="darkgreen") +
  stat_poly_eq(use_label(c("adj.R2" ,"p", "n"))) +
  geom_point(size=1, color="darkgreen") +
  facet_wrap(~stunted, ncol=8) + 
  theme_classic() +
  scale_x_reverse()

ggplot(data = qPCR_Metadata_presence, aes(x = Cq, y = haz_cont)) +
  stat_poly_line(color="darkgreen") +
  stat_poly_eq(use_label(c("adj.R2" ,"p", "n"))) +
  geom_point(size=1, color="darkgreen") +
  theme_classic() +
  scale_x_reverse()

lm1 <- lm(haz_cont ~ Cq, data = qPCR_Metadata_presence)
summary(lm1)




##################################################################


##################################################################
###                     2 WAY ANCOVA                           ###
##################################################################

ggscatter(
  qPCR_Metadata_presence, x = "Cq", y = "haz_cont",
  facet.by  = c("stunted", "pays"), 
  short.panel.labs = FALSE
)+
  stat_smooth(method = "loess", span = 0.9)

qPCR_Metadata_presence %>%
  anova_test(
    haz_cont ~ Cq + stunted + pays + 
      stunted*pays + Cq*stunted +
      Cq*pays + Cq*pays*stunted
  )

qPCR_Metadata_presence %>%
  unite(col = "group", stunted, pays) %>%
  anova_test(haz_cont ~ group*Cq)

# Fit the model, the covariate goes first
model <- lm(haz_cont ~ Cq + stunted*pays, data = qPCR_Metadata_presence)
# Inspect the model diagnostic metrics
model.metrics <- augment(model)
head(model.metrics, 3)

# Assess normality of residuals using shapiro wilk test
shapiro_test(model.metrics$.resid)

levene_test(.resid ~ stunted*pays, data = model.metrics)

res.aov <- qPCR_Metadata_presence %>% 
  anova_test(haz_cont ~ Cq + stunted*pays)
get_anova_table(res.aov)

# Effect of treatment at each level of exercise
qPCR_Metadata_presence %>%
  group_by(pays) %>%
  anova_test(haz_cont ~ Cq + stunted)

qPCR_Metadata_presence %>%
  group_by(stunted) %>%
  anova_test(haz_cont ~ Cq + pays)


##################################################################

####################################################################################################################################
####################################################################################################################################
###                                         Comparison 16s and Metadata                                                 ###
####################################################################################################################################
####################################################################################################################################

##################################################################
###               COMBINE qPCR and METADATA                    ###
##################################################################

#Check if samples match with MISEQ data
qPCR_Metadata <- merge(meta_shannon,metadata, by.x = "row.names", by.y = "id", all= FALSE);dim(qPCR_Metadata)


#Reanotate variable
qPCR_Metadata$stunted <- gsub("Oui","Stunted",qPCR_Metadata$stunted)
qPCR_Metadata$stunted <- gsub("Non","Non-Stunted",qPCR_Metadata$stunted)


##################################################################


##################################################################
## Linear regression 16s ASV and haz
##################################################################


#Shapiro test

shapiro.test(qPCR_Metadata$Number_ASVs)
hist(qPCR_Metadata$Number_ASVs)
#Not normal distributed

qPCR_Metadata$haz_cont <- as.numeric(qPCR_Metadata$haz_cont)
ggplot(data = qPCR_Metadata, aes(x = Number_ASVs, y = haz_cont)) +
  stat_poly_line(color="darkgreen") +
  stat_poly_eq(use_label(c("adj.R2" ,"p", "n"))) +
  geom_point(size=1, color="darkgreen") +
  facet_wrap(~pays+stunted, ncol=8) + 
  theme_classic() +
  scale_x_reverse()

ggplot(data = qPCR_Metadata, aes(x = Number_ASVs, y = haz_cont)) +
  stat_poly_line(color="darkgreen") +
  stat_poly_eq(use_label(c("adj.R2" ,"p", "n"))) +
  geom_point(size=1, color="darkgreen") +
  facet_wrap(~pays, ncol=8) + 
  theme_classic() +
  scale_x_reverse()

ggplot(data = qPCR_Metadata, aes(x = Number_ASVs, y = haz_cont)) +
  stat_poly_line(color="darkgreen") +
  stat_poly_eq(use_label(c("adj.R2" ,"p", "n"))) +
  geom_point(size=1, color="darkgreen") +
  facet_wrap(~stunted, ncol=8) + 
  theme_classic() +
  scale_x_reverse()

ggplot(data = qPCR_Metadata, aes(x = Number_ASVs, y = haz_cont)) +
  stat_poly_line(color="darkgreen") +
  stat_poly_eq(use_label(c("adj.R2" ,"p", "n"))) +
  geom_point(size=1, color="darkgreen") +
  theme_classic() +
  scale_x_reverse()

lm1 <- lm(haz_cont ~ Number_ASVs, data = qPCR_Metadata)
summary(lm1)




##################################################################


##################################################################
###                     2 WAY ANCOVA                           ###
##################################################################

ggscatter(
  qPCR_Metadata, x = "Number_ASVs", y = "haz_cont",
  facet.by  = c("stunted", "pays"), 
  short.panel.labs = FALSE
)+
  stat_smooth(method = "loess", span = 0.9)

qPCR_Metadata %>%
  anova_test(
    haz_cont ~ Number_ASVs + stunted + pays + 
      stunted*pays + Number_ASVs*stunted +
      Number_ASVs*pays + Number_ASVs*pays*stunted
  )

qPCR_Metadata %>%
  unite(col = "group", stunted, pays) %>%
  anova_test(haz_cont ~ group*Number_ASVs)

# Fit the model, the covariate goes first
model <- lm(haz_cont ~ Number_ASVs + stunted*pays, data = qPCR_Metadata)
# Inspect the model diagnostic metrics
model.metrics <- augment(model)
head(model.metrics, 3)

# Assess normality of residuals using shapiro wilk test
shapiro_test(model.metrics$.resid)

levene_test(.resid ~ stunted*pays, data = model.metrics)

res.aov <- qPCR_Metadata %>% 
  anova_test(haz_cont ~ Number_ASVs + stunted*pays)
get_anova_table(res.aov)

# Effect of treatment at each level of exercise
qPCR_Metadata %>%
  group_by(pays) %>%
  anova_test(haz_cont ~ Number_ASVs + stunted)

qPCR_Metadata %>%
  group_by(stunted) %>%
  anova_test(haz_cont ~ Number_ASVs + pays)


##################################################################


####################################################################################################################################
####################################################################################################################################
###                                         Comparison Blastocystis and 18s                                                ###
####################################################################################################################################
####################################################################################################################################


##################################################################
###               COMBINE qPCR and METADATA                    ###
##################################################################

#Check if samples match with MISEQ data
qPCR_Metadata <- merge(plate3,metadata, by.x = "Sample", by.y = "id", all= FALSE);dim(qPCR_Metadata)


qPCR_Metadata$Cq <- as.numeric(qPCR_Metadata$Cq)
qPCR_Metadata$Cq_presence_absence <- ifelse(qPCR_Metadata$Cq<=35 & qPCR_Metadata$Cq>10,"Presence_qPCR", "Absence_qPCR")

#Combined qPCR and metadata
qPCR_Metadata$Cq_presence_absence <- as.factor(qPCR_Metadata$Cq_presence_absence)

#Reanotate variable
qPCR_Metadata$stunted <- gsub("Oui","Stunted",qPCR_Metadata$stunted)
qPCR_Metadata$stunted <- gsub("Non","Non-Stunted",qPCR_Metadata$stunted)


##################################################################
