#Protocol for taking a phyloseq object and doing DESEQ2 tests#
#Note that all work is done on absolute numbers#
#Followed protocol in the following document from phyloseq#
#https://www.bioconductor.org/packages/devel/bioc/vignettes/phyloseq/inst/doc/phyloseq-mixture-models.html#

#Note on Muck Subcategory: Biosample = Manuscript: Muck = 3 muck characteristics (3MC), Mucky = 2 muck characteristics (2MC), Muckish = 1 muck characteristics (1MC), Not = 0 muck characteristics (0MC)#

#Note on TOM/Cu aka LOI.Cu Subcategory: Biosample = Manuscript: High.High = high total organic matter (TOM)/low copper (Cu) (HiHi); High-Low = high TOM/low Cu (HiLo); Low-Low = low TOM/low Cu (LoLo); Low-High = low TOM/high Cu (LoHi)#

#Load libraries#
library(phyloseq)
library(indicspecies)
library(vegan)
library(ggplot2)
library(DESeq2)
library(dplyr)
library(plyr)
library(tidyverse)
library(reshape)
library(FSA)
library(factoextra)
library(cluster)
library(arsenal)
library(rcompanion)

#Load the following function for DESeq2 Analysis#
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

##ASV##

getwd() 
setwd("C:/Users/bradshawd/Documents/Journal_Articles/New_Indicators/ASV")
getwd() #"C:/Users/bradshawd/Documents/Journal_Articles/New_Indicators/ASV"

###UPLOADING AND GENERAL EDITS###
#Import input files from QIIME2 if snakefile worked to phyloseq.biom#
BIOM <- import_biom(file.choose()) #phyloseq.biom
TREE =  read_tree(file.choose()) #file = tree -> tree.nwk
META <- import_qiime_sample_data(file.choose()) #map.txt

#Merge data#
dataA <- merge_phyloseq (BIOM,TREE,META)

#Alternative Uploading if snakefile did not work to phyloseq.biom#
asv_table = read.delim(file.choose(), row.names=1, header = T) #file = "feature-table.tsv    #open file in Excel/TextEdit and remove 1st row in file before import
taxa_table = read.delim(file.choose(), row.names = 1) #file= "taxonomy -> taxonomy2.tsv"   #need to delete Confidence column and separate Taxonomy column using "Text to Columns" (Data tab) in Excel by "semicolon" and fill column headers with taxonomic ranks
taxa_table = as.matrix(taxa_table)
DottedMETA = read.delim(file.choose(), row.names=1) #mapping file
dotted_meta = sample_data(DottedMETA) 
ASV = otu_table(asv_table, taxa_are_rows = TRUE)
TAX = tax_table(taxa_table)
TREE =  read_tree(file.choose()) #file = tree -> tree.nwk

dataA <- merge_phyloseq(ASV, TAX, dotted_meta, TREE)

#Change taxa names if needed#
colnames(tax_table(dataA))
colnames(tax_table(dataA))=c("Kingdom","Phylum","Class","Order","Family","Genus", "Species")
colnames(tax_table(dataA))

###Create at more filtered set of LWSS Data###
ntaxa(dataA) #166944#
nsamples(dataA)#483

LWSDataA <- subset_samples(dataA, Survey=="Lagoon Wide")
nsamples(LWSDataA)#324

LWSWDataA <- subset_samples(LWSDataA, Medium=="Water")
nsamples(LWSWDataA)#132, 11 sites in triplicate (33) x 4 sampling periods (sps)

LWSSDataA1 <- subset_samples(dataA, Survey=="Lagoon Wide-Post Hurricane")
nsamples(LWSSDataA1)#12

LWSSDataA2 <- subset_samples(LWSDataA, Medium=="Sediment")
nsamples(LWSSDataA2)#192

LWSSDataA <- merge_phyloseq(LWSSDataA1, LWSSDataA2)
nsamples(LWSSDataA)#204

##Create phyloseqs for determining and testing indicators for Sediment
#Split samples by sampling period (SP)
LWSS_W16_DataA <- subset_samples(LWSSDataA, AbbSeasonAbbYear=="W16")
nsamples(LWSS_W16_DataA)#45, 15 sites in triplicate

LWSS_D17_DataA <- subset_samples(LWSSDataA, AbbSeasonAbbYear=="D17")
nsamples(LWSS_D17_DataA)#45,15 sites in triplicate

LWSS_W17_DataA <- subset_samples(LWSSDataA, AbbSeasonAbbYear=="W17")
nsamples(LWSS_W17_DataA)#57, 19 sites in triplicate

LWSS_D18_DataA <- subset_samples(LWSSDataA, AbbSeasonAbbYear=="D18")
nsamples(LWSS_D18_DataA)#57, 19 sites in triplicate

#Combine sampling periods into first and second years
LWSS_Y1_DataA <- merge_phyloseq(LWSS_W16_DataA, LWSS_D17_DataA)
nsamples(LWSS_Y1_DataA)#90

LWSS_Y2_DataA <- merge_phyloseq(LWSS_W17_DataA, LWSS_D18_DataA)
nsamples(LWSS_Y2_DataA)#114

##Create phyloseqs for determining and testing indicators for Water
#Split samples by sampling period (SP)
LWSW_W16_DataA <- subset_samples(LWSWDataA, AbbSeasonAbbYear=="W16")
nsamples(LWSW_W16_DataA)#33, 11 sites in triplicate

LWSW_D17_DataA <- subset_samples(LWSWDataA, AbbSeasonAbbYear=="D17")
nsamples(LWSW_D17_DataA)#33,11 sites in triplicate

LWSW_W17_DataA <- subset_samples(LWSWDataA, AbbSeasonAbbYear=="W17")
nsamples(LWSW_W17_DataA)#33, 11 sites in triplicate

LWSW_D18_DataA <- subset_samples(LWSWDataA, AbbSeasonAbbYear=="D18")
nsamples(LWSW_D18_DataA)#33, 11 sites in triplicate

#Combine sampling periods into first and second years
LWSW_Y1_DataA <- merge_phyloseq(LWSW_W16_DataA, LWSW_D17_DataA)
nsamples(LWSW_Y1_DataA)#66

LWSW_Y2_DataA <- merge_phyloseq(LWSW_W17_DataA, LWSW_D18_DataA)
nsamples(LWSW_Y2_DataA)#66

##Filter yearly sample sets of all zero and then low abundance asvs
FLWSS_Y1_DataA <- filter_taxa(LWSS_Y1_DataA, function(x) sum(x) >0, TRUE)
ntaxa(FLWSS_Y1_DataA) #54839#
RFLWSS_Y1_DataAasv <- filter_taxa(FLWSS_Y1_DataA, function(x) sum(x) >20, TRUE)
ntaxa(RFLWSS_Y1_DataAasv) #7472#

FLWSS_Y2_DataA <- filter_taxa(LWSS_Y2_DataA, function(x) sum(x) >0, TRUE)
ntaxa(FLWSS_Y2_DataA) #78386#
RFLWSS_Y2_DataAasv <- filter_taxa(FLWSS_Y2_DataA, function(x) sum(x) >20, TRUE)
ntaxa(RFLWSS_Y2_DataAasv) #4169#

FLWSW_Y1_DataA <- filter_taxa(LWSW_Y1_DataA, function(x) sum(x) >0, TRUE)
ntaxa(FLWSW_Y1_DataA) #8189#
RFLWSW_Y1_DataAasv <- filter_taxa(FLWSW_Y1_DataA, function(x) sum(x) >20, TRUE)
ntaxa(RFLWSW_Y1_DataAasv) #1672#

FLWSW_Y2_DataA <- filter_taxa(LWSW_Y2_DataA, function(x) sum(x) >0, TRUE)
ntaxa(FLWSW_Y2_DataA) #17357#
RFLWSW_Y2_DataAasv <- filter_taxa(FLWSW_Y2_DataA, function(x) sum(x) >20, TRUE)
ntaxa(RFLWSW_Y2_DataAasv) #1471#





#Determine the number of samples in each category and subcategory#
length(which(sample_data(RFLWSS_Y1_DataAasv)$Estuary == "SLE")) #24
length(which(sample_data(RFLWSS_Y1_DataAasv)$Estuary == "IRL")) #66

length(which(sample_data(RFLWSS_Y1_DataAasv)$Muck == "Not")) #60
length(which(sample_data(RFLWSS_Y1_DataAasv)$Muck == "Muck")) #22
length(which(sample_data(RFLWSS_Y1_DataAasv)$Muck == "Mucky")) #4
length(which(sample_data(RFLWSS_Y1_DataAasv)$Muck == "Muckish")) #4

length(which(sample_data(RFLWSS_Y1_DataAasv)$LOI.Cu == "Low-Low")) #61
length(which(sample_data(RFLWSS_Y1_DataAasv)$LOI.Cu == "High-Low")) #18
length(which(sample_data(RFLWSS_Y1_DataAasv)$LOI.Cu == "High-High")) #11
length(which(sample_data(RFLWSS_Y1_DataAasv)$LOI.Cu == "Low-High")) #0



##NOTE##
#Use res to determine what is + or -#
#muck res: Not (+) vs Muck (-)#
#tomcu res: HL (+) vs HH (-)#
#estuary res: SLE (+) vs IRL (-)#

###RUN DESEQ2 AT asv LEVEL ON ESTUARY SUBCATEGORIES###
#Convert phyloseq to deseq2 object centered around a irlvssle factor#
estuary_cudds_asv = phyloseq_to_deseq2(RFLWSS_Y1_DataAasv,  ~ Estuary)

#Calculate geometric means prior to estimate size factors#
estuary_geoMeans_asv = apply(counts(estuary_cudds_asv), 1, gm_mean)
estuary_cudds_asv = estimateSizeFactors(estuary_cudds_asv, geoMeans = estuary_geoMeans_asv)

#Conduct DESEQ2 test#
estuary_cudds_asv = DESeq(estuary_cudds_asv, fitType="local")

#Explore the results#
estuary_DESeq2_res_asv = results(estuary_cudds_asv)
estuary_DESeq2_res_asv = estuary_DESeq2_res_asv[order(estuary_DESeq2_res_asv$padj, na.last=NA), ]
alpha = 0.05
estuary_DESeq2_sig_res_asv = estuary_DESeq2_res_asv[(estuary_DESeq2_res_asv$padj < alpha), ]

#Make dataframe with taxanomy added in#
estuary_DESeq2_sig_res_taxo_asv = cbind(as(estuary_DESeq2_sig_res_asv, "data.frame"), as(tax_table(RFLWSS_Y1_DataAasv)[rownames(estuary_DESeq2_sig_res_asv), ], "matrix"))

#Make dataframe with asvs added in from all sites#
estuary_DESeq2_sig_res_taxo_seqs_asv = cbind(as(estuary_DESeq2_sig_res_taxo_asv, "data.frame"), as(otu_table(RFLWSS_Y1_DataAasv)[rownames(estuary_DESeq2_sig_res_taxo_asv), ], "matrix"))

#Make rownames an actual column and remove old rownames#
estuary_DESeq2_sig_res_taxo_seqs_asv <- cbind(ESV.ID = rownames(estuary_DESeq2_sig_res_taxo_seqs_asv), estuary_DESeq2_sig_res_taxo_seqs_asv)
rownames(estuary_DESeq2_sig_res_taxo_seqs_asv) <-NULL

#Make a txt file of results for all sites#
write.csv(as.data.frame(estuary_DESeq2_sig_res_taxo_seqs_asv), file="estuary_DESeq2_sig_res_taxo_seqs_asv.csv")

#Detemine which subcategory is negative or positive#
estuary_DESeq2_res_asv
#estuary_DESeq2_res_asv: SLE (+) vs IRL (-)#

#Determine number of Indicators with padj >0.05 in each subcagetory#
#Number of SLE indicators# 
length(which(estuary_DESeq2_sig_res_taxo_seqs_asv$log2FoldChange > 0)) #771
#Number of IRL indicators#
length(which(estuary_DESeq2_sig_res_taxo_seqs_asv$log2FoldChange < 0)) #1506

###RUN DESEQ2 AT asv LEVEL ON MUCK EXTREME SUBCATEGORIES###

##Prepare for running DESeq2##
#Make phyloseq focusing on Muck and Not subcategories#  
RFLWSS_Y1_Muck_DataAasv = subset_samples(RFLWSS_Y1_DataAasv, Muck=="Muck")
RFLWSS_Y1_Not_DataAasv = subset_samples(RFLWSS_Y1_DataAasv, Muck=="Not")
RFLWSS_Y1_MuckvsNot_DataAasv = merge_phyloseq(RFLWSS_Y1_Muck_DataAasv, RFLWSS_Y1_Not_DataAasv)

##Run DESEQ2 and make a txt table of results##
#Convert phyloseq to deseq2 object centered around a muck factor#
muck_cudds_asv = phyloseq_to_deseq2(RFLWSS_Y1_MuckvsNot_DataAasv,  ~ Muck)

#Calculate geometric means prior to estimate size factors#
muck_geoMeans_asv = apply(counts(muck_cudds_asv), 1, gm_mean)
muck_cudds_asv = estimateSizeFactors(muck_cudds_asv, geoMeans = muck_geoMeans_asv)

#Conduct DESEQ2 test#
muck_cudds_asv = DESeq(muck_cudds_asv, fitType="local")

#Explore the results#
muck_DESeq2_res_asv = results(muck_cudds_asv)
muck_DESeq2_res_asv = muck_DESeq2_res_asv[order(muck_DESeq2_res_asv$padj, na.last=NA), ]
alpha = 0.05
muck_DESeq2_sig_res_asv = muck_DESeq2_res_asv[(muck_DESeq2_res_asv$padj < alpha), ]

#Make dataframe with taxanomy added in#
muck_DESeq2_sig_res_taxo_asv = cbind(as(muck_DESeq2_sig_res_asv, "data.frame"), as(tax_table(RFLWSS_Y1_MuckvsNot_DataAasv)[rownames(muck_DESeq2_sig_res_asv), ], "matrix"))

#Make dataframe with asvs added in from all sites#
muck_DESeq2_sig_res_taxo_seqs_asv = cbind(as(muck_DESeq2_sig_res_taxo_asv, "data.frame"), as(otu_table(RFLWSS_Y1_DataAasv)[rownames(muck_DESeq2_sig_res_taxo_asv), ], "matrix"))

#Make rownames an actual column and remove old rownames#
muck_DESeq2_sig_res_taxo_seqs_asv <- cbind(ESV.ID = rownames(muck_DESeq2_sig_res_taxo_seqs_asv), muck_DESeq2_sig_res_taxo_seqs_asv)
rownames(muck_DESeq2_sig_res_taxo_seqs_asv) <-NULL

#Make a txt file of results for all sites#
write.csv(as.data.frame(muck_DESeq2_sig_res_taxo_seqs_asv), file="muck_DESeq2_sig_res_taxo_seqs_asv.csv")

#Detemine which subcategory is negative or positive#
muck_DESeq2_res_asv
#muck_DESeq2_res_asv: Not (+) vs Muck (-)#
#Determine number of Indicators with padj >0.05 in each subcagetory#
#Number of Not indicators# 
length(which(muck_DESeq2_sig_res_taxo_seqs_asv$log2FoldChange > 0)) #670
#Number of Muck indicators#
length(which(muck_DESeq2_sig_res_taxo_seqs_asv$log2FoldChange < 0)) #149

###RUN DESEQ2 AT asv LEVEL ON TOM/Cu EXTREME SUBCATEGORIES###

##Prepare for running DESeq2##
#Make phyloseq focusing on High-High and High-Low subcategories#  
RFLWSS_Y1_HH_DataAasv = subset_samples(RFLWSS_Y1_DataAasv, LOI.Cu=="High-High")
RFLWSS_Y1_HL_DataAasv = subset_samples(RFLWSS_Y1_DataAasv, LOI.Cu=="High-Low")
RFLWSS_Y1_HHvsHL_DataAasv = merge_phyloseq(RFLWSS_Y1_HH_DataAasv, RFLWSS_Y1_HL_DataAasv)

##Run DESEQ2 and make a txt table of results##
#Convert phyloseq to deseq2 object centered around a loicu factor#
tomcu_cudds_asv = phyloseq_to_deseq2(RFLWSS_Y1_HHvsHL_DataAasv,  ~ LOI.Cu)

#Calculate geometric means prior to estimate size factors#
tomcu_geoMeans_asv = apply(counts(tomcu_cudds_asv), 1, gm_mean)
tomcu_cudds_asv = estimateSizeFactors(tomcu_cudds_asv, geoMeans = tomcu_geoMeans_asv)

#Conduct DESEQ2 test#
tomcu_cudds_asv = DESeq(tomcu_cudds_asv, fitType="local")

#Explore the results#
tomcu_DESeq2_res_asv = results(tomcu_cudds_asv)
tomcu_DESeq2_res_asv = tomcu_DESeq2_res_asv[order(tomcu_DESeq2_res_asv$padj, na.last=NA), ]
alpha = 0.05
tomcu_DESeq2_sig_res_asv = tomcu_DESeq2_res_asv[(tomcu_DESeq2_res_asv$padj < alpha), ]

#Make dataframe with taxanomy added in#
tomcu_DESeq2_sig_res_taxo_asv = cbind(as(tomcu_DESeq2_sig_res_asv, "data.frame"), as(tax_table(RFLWSS_Y1_HHvsHL_DataAasv)[rownames(tomcu_DESeq2_sig_res_asv), ], "matrix"))

#Make dataframe with asvs added in from all sites#
tomcu_DESeq2_sig_res_taxo_seqs_asv = cbind(as(tomcu_DESeq2_sig_res_taxo_asv, "data.frame"), as(otu_table(RFLWSS_Y1_DataAasv)[rownames(tomcu_DESeq2_sig_res_taxo_asv), ], "matrix"))

#Make rownames an actual column and remove old rownames#
tomcu_DESeq2_sig_res_taxo_seqs_asv <- cbind(ESV.ID = rownames(tomcu_DESeq2_sig_res_taxo_seqs_asv), tomcu_DESeq2_sig_res_taxo_seqs_asv)
rownames(tomcu_DESeq2_sig_res_taxo_seqs_asv) <-NULL

#Make a txt file of results for all sites#
write.csv(as.data.frame(tomcu_DESeq2_sig_res_taxo_seqs_asv), file="tomcu_DESeq2_sig_res_taxo_seqs_asv.csv")

#Detemine which subcategory is negative or positive#
tomcu_DESeq2_res_asv
#tomcu_DESeq2_res_asv: High.Low (+) vs High.High (-)#
#Determine number of Indicators with padj >0.05 in each subcagetory#
#Number of High.Low indicators# 
length(which(tomcu_DESeq2_sig_res_taxo_seqs_asv$log2FoldChange > 0)) #417
#Number of High.High indicators#
length(which(tomcu_DESeq2_sig_res_taxo_seqs_asv$log2FoldChange < 0)) #331

#Order rows by log2FoldChange#
muck_DESeq2_sig_res_taxo_seqs_asv <- arrange(muck_DESeq2_sig_res_taxo_seqs_asv, log2FoldChange)
tomcu_DESeq2_sig_res_taxo_seqs_asv <- arrange(tomcu_DESeq2_sig_res_taxo_seqs_asv, log2FoldChange)
estuary_DESeq2_sig_res_taxo_seqs_asv <- arrange(estuary_DESeq2_sig_res_taxo_seqs_asv, log2FoldChange)

#Create lists of just significant asvs for each category#
unfiltered_DESeq2_estuary_asvs <- subset(estuary_DESeq2_sig_res_taxo_seqs_asv, select=c(ESV.ID))
unfiltered_DESeq2_muck_asvs <- subset(muck_DESeq2_sig_res_taxo_seqs_asv, select=c(ESV.ID))
unfiltered_DESeq2_tomcu_asvs <- subset(tomcu_DESeq2_sig_res_taxo_seqs_asv, select=c(ESV.ID))



##Create indicator tables filtered of overlapping indicators with other two metadata categories##

#Filter out rows in the IRLSLE DESEQ object that are also in the LOI.Cu and Muck DESEQs#
filtered_estuary_DESeq2_sig_res_taxo_seqs_asv <- anti_join(estuary_DESeq2_sig_res_taxo_seqs_asv, unfiltered_DESeq2_tomcu_asvs, by = "ESV.ID")
filtered_estuary_DESeq2_sig_res_taxo_seqs_asv <- anti_join(filtered_estuary_DESeq2_sig_res_taxo_seqs_asv, unfiltered_DESeq2_muck_asvs, by = "ESV.ID")

#Filter out rows in the Muck DESEQ object that are also in the IRL.SLE and LOI.Cu DESEQs#
filtered_muck_DESeq2_sig_res_taxo_seqs_asv <- anti_join(muck_DESeq2_sig_res_taxo_seqs_asv, unfiltered_DESeq2_estuary_asvs, by = "ESV.ID")
filtered_muck_DESeq2_sig_res_taxo_seqs_asv <- anti_join(filtered_muck_DESeq2_sig_res_taxo_seqs_asv, unfiltered_DESeq2_tomcu_asvs, by = "ESV.ID")

#Filter out rows in the LOI.CU DESEQ object that are also in the IRL.SLE and Muck DESEQs#
filtered_tomcu_DESeq2_sig_res_taxo_seqs_asv <- anti_join(tomcu_DESeq2_sig_res_taxo_seqs_asv, unfiltered_DESeq2_estuary_asvs, by = "ESV.ID")
filtered_tomcu_DESeq2_sig_res_taxo_seqs_asv <- anti_join(filtered_tomcu_DESeq2_sig_res_taxo_seqs_asv, unfiltered_DESeq2_muck_asvs, by = "ESV.ID")



#Order rows by log2FoldChange#
filtered_muck_DESeq2_sig_res_taxo_seqs_asv <- arrange(filtered_muck_DESeq2_sig_res_taxo_seqs_asv, log2FoldChange)
filtered_tomcu_DESeq2_sig_res_taxo_seqs_asv <- arrange(filtered_tomcu_DESeq2_sig_res_taxo_seqs_asv, log2FoldChange)
filtered_estuary_DESeq2_sig_res_taxo_seqs_asv <- arrange(filtered_estuary_DESeq2_sig_res_taxo_seqs_asv, log2FoldChange)

#Make a csv file for each of the filtered tables#
write.csv(as.data.frame(filtered_estuary_DESeq2_sig_res_taxo_seqs_asv), file="filtered_estuary_DESeq2_sig_res_taxo_seqs_asv.csv")
write.csv(as.data.frame(filtered_muck_DESeq2_sig_res_taxo_seqs_asv), file="filtered_muck_DESeq2_sig_res_taxo_seqs_asv.csv")
write.csv(as.data.frame(filtered_tomcu_DESeq2_sig_res_taxo_seqs_asv), file="filtered_tomcu_DESeq2_sig_res_taxo_seqs_asv.csv")



##Determine number of filtered indicators##

#Number of Not indicators# 
length(which(filtered_muck_DESeq2_sig_res_taxo_seqs_asv$log2FoldChange > 0)) #16
#Number of Muck indicators#
length(which(filtered_muck_DESeq2_sig_res_taxo_seqs_asv$log2FoldChange < 0)) #56
#Number of High.Low indicators# 
length(which(filtered_tomcu_DESeq2_sig_res_taxo_seqs_asv$log2FoldChange > 0)) #39
#Number of High.High indicators#
length(which(filtered_tomcu_DESeq2_sig_res_taxo_seqs_asv$log2FoldChange < 0)) #142
#Number of SLE indicators# 
length(which(filtered_estuary_DESeq2_sig_res_taxo_seqs_asv$log2FoldChange > 0)) #379
#Number of IRL indicators#
length(which(filtered_estuary_DESeq2_sig_res_taxo_seqs_asv$log2FoldChange < 0)) #674

#Create lists of filtered significant asvs for each category#
filtered_DESeq2_estuary_asvs <- subset(filtered_estuary_DESeq2_sig_res_taxo_seqs_asv, select=c(ESV.ID))
filtered_DESeq2_muck_asvs <- subset(filtered_muck_DESeq2_sig_res_taxo_seqs_asv, select=c(ESV.ID))
filtered_DESeq2_tomcu_asvs <- subset(filtered_tomcu_DESeq2_sig_res_taxo_seqs_asv, select=c(ESV.ID))

##Create phyloseq objects for DESeq2 indicators##

#Create phyloseq objects of the unfiltered indicators using character string version of significant asvs##
charac_unfiltered_DESeq2_estuary_asvs <- unfiltered_DESeq2_estuary_asvs[,"ESV.ID"]
charac_unfiltered_DESeq2_estuary_asvs <- as.character(charac_unfiltered_DESeq2_estuary_asvs)
UFDestuaryDataAasv <- prune_taxa(charac_unfiltered_DESeq2_estuary_asvs, RFLWSS_Y1_DataAasv)

charac_unfiltered_DESeq2_muck_asvs <- unfiltered_DESeq2_muck_asvs[,"ESV.ID"]
charac_unfiltered_DESeq2_muck_asvs <- as.character(charac_unfiltered_DESeq2_muck_asvs)
UFDmuckDataAasv <- prune_taxa(charac_unfiltered_DESeq2_muck_asvs, RFLWSS_Y1_DataAasv)

charac_unfiltered_DESeq2_tomcu_asvs <- unfiltered_DESeq2_tomcu_asvs[,"ESV.ID"]
charac_unfiltered_DESeq2_tomcu_asvs <- as.character(charac_unfiltered_DESeq2_tomcu_asvs)
UFDtomcuDataAasv <- prune_taxa(charac_unfiltered_DESeq2_tomcu_asvs, RFLWSS_Y1_DataAasv)

#Create phyloseq objects of the filtered indicators using character string version of significant asvs##
charac_filtered_DESeq2_estuary_asvs <- filtered_DESeq2_estuary_asvs[,"ESV.ID"]
charac_filtered_DESeq2_estuary_asvs <- as.character(charac_filtered_DESeq2_estuary_asvs)
FDestuaryDataAasv <- prune_taxa(charac_filtered_DESeq2_estuary_asvs, RFLWSS_Y1_DataAasv)

charac_filtered_DESeq2_muck_asvs <- filtered_DESeq2_muck_asvs[,"ESV.ID"]
charac_filtered_DESeq2_muck_asvs <- as.character(charac_filtered_DESeq2_muck_asvs)
FDmuckDataAasv <- prune_taxa(charac_filtered_DESeq2_muck_asvs, RFLWSS_Y1_DataAasv)

charac_filtered_DESeq2_tomcu_asvs <- filtered_DESeq2_tomcu_asvs[,"ESV.ID"]
charac_filtered_DESeq2_tomcu_asvs <- as.character(charac_filtered_DESeq2_tomcu_asvs)
FDtomcuDataAasv <- prune_taxa(charac_filtered_DESeq2_tomcu_asvs, RFLWSS_Y1_DataAasv)



##Create heatmaps to see overall trends, may take a while to load##

unfiltered_DESeq2_estuary_asv_heatmap <- plot_heatmap(UFDestuaryDataAasv, taxa.order = charac_unfiltered_DESeq2_estuary_asvs, sample.order= "Estuary", sample.label = "Estuary")
unfiltered_DESeq2_estuary_asv_heatmap

unfiltered_DESeq2_muck_asv_heatmap <- plot_heatmap(UFDmuckDataAasv, taxa.order = charac_unfiltered_DESeq2_muck_asvs, sample.order= "Muck", sample.label = "Muck")
unfiltered_DESeq2_muck_asv_heatmap

unfiltered_DESeq2_tomcu_asv_heatmap <- plot_heatmap(UFDtomcuDataAasv, taxa.order = charac_unfiltered_DESeq2_tomcu_asvs, sample.order= "LOI.Cu", sample.label = "LOI.Cu")
unfiltered_DESeq2_tomcu_asv_heatmap

filtered_DESeq2_estuary_asv_heatmap <- plot_heatmap(FDestuaryDataAasv, taxa.order = charac_filtered_DESeq2_estuary_asvs, sample.order= "Estuary", sample.label = "Estuary")
filtered_DESeq2_estuary_asv_heatmap

filtered_DESeq2_muck_asv_heatmap <- plot_heatmap(FDmuckDataAasv, taxa.order = charac_filtered_DESeq2_muck_asvs, sample.order= "Muck", sample.label = "Muck")
filtered_DESeq2_muck_asv_heatmap

filtered_DESeq2_tomcu_asv_heatmap <- plot_heatmap(FDtomcuDataAasv, taxa.order = charac_filtered_DESeq2_tomcu_asvs, sample.order= "LOI.Cu", sample.label = "LOI.Cu")
filtered_DESeq2_tomcu_asv_heatmap


###CONDUCT INDICSPECIES ANALYSIS AT asv LEVEL###

#Increase max print so that all of results can be printed#
options(max.print=100000000)


##Conduct asv indicspecies analysis for Estuary Category#

#Take out OTU aka asv table from metadata focused phyloseq object#
estuary_seqs_asv = as(otu_table(RFLWSS_Y1_DataAasv), "matrix")

#Make it a data frame, transpose it, make it a data frame again#
estuary_seqs_asv <- as.data.frame(estuary_seqs_asv)
estuary_seqs_asv<- t(estuary_seqs_asv)
estuary_seqs_asv <- as.data.frame(estuary_seqs_asv)

#Take out metadata table metadata focused phyloseq object#
estuary_meta = as(sample_data(RFLWSS_Y1_DataAasv), "matrix")

#Keep only the sampleids and metadata category focusing on#
estuary_meta <- subset(estuary_meta, select=c(sampleid, Estuary))

#Reset names of two factors#
estuary_meta[estuary_meta == "IRL"] <- 1
estuary_meta[estuary_meta == "SLE"] <- 2

#Convert it to data frame, then create class object based on metadatacategory#
estuary_meta = as.data.frame(estuary_meta)
estuary_meta_group <- as.character(estuary_meta$Estuary)

#Run indicspecies#
estuary_indval_asv = multipatt(estuary_seqs_asv, estuary_meta_group, control = how(nperm=9999))

#Create a blank csv, then sink all indicspecies results into it#
sink("estuary_IS_res_asv.csv")
estuary_sig_indval_asv <- summary(estuary_indval_asv, indvalcomp=TRUE, alpha=1)
sink()

##Conduct asv indicspecies analysis for Muck Category#

#Take out OTU aka asv table from metadata focused phyloseq object#
muck_seqs_asv = as(otu_table(RFLWSS_Y1_MuckvsNot_DataAasv), "matrix")

#Make it a data frame, transpose it, make it a data frame again#
muck_seqs_asv <- as.data.frame(muck_seqs_asv)
muck_seqs_asv<- t(muck_seqs_asv)
muck_seqs_asv <- as.data.frame(muck_seqs_asv)

#Take out metadata table metadata focused phyloseq object#
muck_meta = as(sample_data(RFLWSS_Y1_MuckvsNot_DataAasv), "matrix")

#Keep only the sampleids and metadata category focusing on#
muck_meta <- subset(muck_meta, select=c(sampleid, Muck))

#Reset names of two factors#
muck_meta[muck_meta == "Muck"] <- 1
muck_meta[muck_meta == "Not"] <- 2

#Convert it to data frame, then create class object based on metadatacategory#
muck_meta = as.data.frame(muck_meta)
muck_meta_group <- as.character(muck_meta$Muck)

#Run indicspecies#
muck_indval_asv = multipatt(muck_seqs_asv, muck_meta_group, control = how(nperm=9999))

#Create a blank csv, then sink all indicspecies results into it#
sink("muck_IS_res_asv.csv")
muck_sig_indval_asv <- summary(muck_indval_asv, indvalcomp=TRUE, alpha=1)
sink()


##Conduct asv indicspecies analysis for TOM/Cu Category#

#Take out OTU aka asv table from metadata focused phyloseq object#
tomcu_seqs_asv = as(otu_table(RFLWSS_Y1_HHvsHL_DataAasv), "matrix")

#Make it a data frame, transpose it, make it a data frame again#
tomcu_seqs_asv <- as.data.frame(tomcu_seqs_asv)
tomcu_seqs_asv<- t(tomcu_seqs_asv)
tomcu_seqs_asv <- as.data.frame(tomcu_seqs_asv)

#Take out metadata table metadata focused phyloseq object#
tomcu_meta = as(sample_data(RFLWSS_Y1_HHvsHL_DataAasv), "matrix")

#Keep only the sampleids and metadata category focusing on#
tomcu_meta <- subset(tomcu_meta, select=c(sampleid, LOI.Cu))

#Reset names of two factors#
tomcu_meta[tomcu_meta == "High-High"] <- 1
tomcu_meta[tomcu_meta == "High-Low"] <- 2

#Convert it to data frame, then create class object based on metadatacategory#
tomcu_meta = as.data.frame(tomcu_meta)
tomcu_meta_group <- as.character(tomcu_meta$LOI.Cu)

#Run indicspecies#
tomcu_indval_asv = multipatt(tomcu_seqs_asv, tomcu_meta_group, control = how(nperm=9999))

#Create a blank csv, then sink all indicspecies results into it#
sink("tomcu_IS_res_asv.csv")
tomcu_sig_indval_asv <- summary(tomcu_indval_asv, indvalcomp=TRUE, alpha=1)
sink()

#Reset max print to default#
options(max.print=99999)

##Look at csv file you sinked into and create a table
##of just ASV.ID, A, B, stat, p.value, and group##

##Continue asv indicspecies analysis for Estuary Category##

#Upload new table#
estuary_IS_res_asv <- read.csv(file.choose())
#Make a new p value adjusted column#
estuary_IS_res_asv$padj <- p.adjust(estuary_IS_res_asv$p.value, method = "BH")

#Keep table of adjusted p values less than 0.05#
estuary_IS_sig_res_asv <- filter(estuary_IS_res_asv, padj < 0.05)

#Make dataframe with taxanomy added in by making ESV.ID the rownames then adding taxonomy table#
estuary_IS_sig_res_taxo_asv <- as.data.frame(estuary_IS_sig_res_asv)
rownames(estuary_IS_sig_res_taxo_asv) <- estuary_IS_sig_res_taxo_asv[,1]
estuary_IS_sig_res_taxo_asv = cbind(as(estuary_IS_sig_res_taxo_asv, "data.frame"), as(tax_table(RFLWSS_Y1_DataAasv)[rownames(estuary_IS_sig_res_taxo_asv), ], "matrix"))

#Make dataframe with asvs added in from all sites#
estuary_IS_sig_res_taxo_seqs_asv = cbind(as(estuary_IS_sig_res_taxo_asv, "data.frame"), as(otu_table(RFLWSS_Y1_DataAasv)[rownames(estuary_IS_sig_res_taxo_asv), ], "matrix"))

#Remove old rownames#
rownames(estuary_IS_sig_res_taxo_seqs_asv) <-NULL

#Make a txt file of results for all sites#
write.csv(as.data.frame(estuary_IS_sig_res_taxo_seqs_asv), file="estuary_IS_sig_res_taxo_seqs_asv.csv")

#Check the number of indicators in each group after padj#
length(which(sample_data(estuary_IS_sig_res_asv)$group == "1")) #1486 IRL
length(which(sample_data(estuary_IS_sig_res_asv)$group == "2")) #2239 SLE


##Continue asv indicspecies analysis for Muck Category##

#Upload new table#
muck_IS_res_asv <- read.csv(file.choose())

#Make a new p value adjusted column#
muck_IS_res_asv$padj <- p.adjust(muck_IS_res_asv$p.value, method = "BH")

#Keep table of adjusted p values less than 0.05#
muck_IS_sig_res_asv <- filter(muck_IS_res_asv, padj < 0.05)

#Make dataframe with taxanomy added in by making ESV.ID the rownames then adding taxonomy table#
muck_IS_sig_res_taxo_asv <- as.data.frame(muck_IS_sig_res_asv)
rownames(muck_IS_sig_res_taxo_asv) <- muck_IS_sig_res_taxo_asv[,1]
muck_IS_sig_res_taxo_asv = cbind(as(muck_IS_sig_res_taxo_asv, "data.frame"), as(tax_table(RFLWSS_Y1_MuckvsNot_DataAasv)[rownames(muck_IS_sig_res_taxo_asv), ], "matrix"))

#Make dataframe with asvs added in from all sites#
muck_IS_sig_res_taxo_seqs_asv = cbind(as(muck_IS_sig_res_taxo_asv, "data.frame"), as(otu_table(RFLWSS_Y1_DataAasv)[rownames(muck_IS_sig_res_taxo_asv), ], "matrix"))

#Remove old rownames#
rownames(muck_IS_sig_res_taxo_seqs_asv) <-NULL

#Make a txt file of results for all sites#
write.csv(as.data.frame(muck_IS_sig_res_taxo_seqs_asv), file="muck_IS_sig_res_taxo_seqs_asv.csv")

#Check the number of indicators in each group after padj#
length(which(sample_data(muck_IS_sig_res_asv)$group == "1")) #1319 Muck
length(which(sample_data(muck_IS_sig_res_asv)$group == "2")) #575 Not

##Continue asv indicspecies analysis for TOM/Cu Category##

#Upload new table#
tomcu_IS_res_asv <- read.csv(file.choose())
#Make a new p value adjusted column#
tomcu_IS_res_asv$padj <- p.adjust(tomcu_IS_res_asv$p.value, method = "BH")

#Keep table of adjusted p values less than 0.05#
tomcu_IS_sig_res_asv <- filter(tomcu_IS_res_asv, padj < 0.05)

#Make dataframe with taxanomy added in by making ESV.ID the rownames then adding taxonomy table#
tomcu_IS_sig_res_taxo_asv <- as.data.frame(tomcu_IS_sig_res_asv)
rownames(tomcu_IS_sig_res_taxo_asv) <- tomcu_IS_sig_res_taxo_asv[,1]
tomcu_IS_sig_res_taxo_asv = cbind(as(tomcu_IS_sig_res_taxo_asv, "data.frame"), as(tax_table(RFLWSS_Y1_MuckvsNot_DataAasv)[rownames(tomcu_IS_sig_res_taxo_asv), ], "matrix"))

#Make dataframe with asvs added in from all sites#
tomcu_IS_sig_res_taxo_seqs_asv = cbind(as(tomcu_IS_sig_res_taxo_asv, "data.frame"), as(otu_table(RFLWSS_Y1_DataAasv)[rownames(tomcu_IS_sig_res_taxo_asv), ], "matrix"))

#Remove old rownames#
rownames(tomcu_IS_sig_res_taxo_seqs_asv) <-NULL

#Make a txt file of results for all sites#
write.csv(as.data.frame(tomcu_IS_sig_res_taxo_seqs_asv), file="tomcu_IS_sig_res_taxo_seqs_asv.csv")

#Check the number of indicators in each group after padj#
length(which(sample_data(tomcu_IS_sig_res_asv)$group == "1")) #296 High/High
length(which(sample_data(tomcu_IS_sig_res_asv)$group == "2")) #273 High/Low

#Order rows by group#
estuary_IS_sig_res_taxo_seqs_asv <- arrange(estuary_IS_sig_res_taxo_seqs_asv, group)
muck_IS_sig_res_taxo_seqs_asv <- arrange(muck_IS_sig_res_taxo_seqs_asv, group)
tomcu_IS_sig_res_taxo_seqs_asv <- arrange(tomcu_IS_sig_res_taxo_seqs_asv, group)


#Create lists of just significant asvs for each category#
unfiltered_IS_estuary_asvs <- subset(estuary_IS_sig_res_taxo_seqs_asv, select=c(ESV.ID))
unfiltered_IS_muck_asvs <- subset(muck_IS_sig_res_taxo_seqs_asv, select=c(ESV.ID))
unfiltered_IS_tomcu_asvs <- subset(tomcu_IS_sig_res_taxo_seqs_asv, select=c(ESV.ID))


##Create indicator tables filtered of overlapping indicators with other two metadata categories##

#Filter out rows in the estuary IS object that are also in the tomcu and muck IS objects#
filtered_estuary_IS_sig_res_taxo_seqs_asv <- anti_join(estuary_IS_sig_res_taxo_seqs_asv, unfiltered_IS_tomcu_asvs, by = "ESV.ID")
filtered_estuary_IS_sig_res_taxo_seqs_asv <- anti_join(filtered_estuary_IS_sig_res_taxo_seqs_asv, unfiltered_IS_muck_asvs, by = "ESV.ID")

#Filter out rows in the muck IS object that are also in the estuary and tomcu IS objects#
filtered_muck_IS_sig_res_taxo_seqs_asv<- anti_join(muck_IS_sig_res_taxo_seqs_asv, unfiltered_IS_estuary_asvs, by = "ESV.ID")
filtered_muck_IS_sig_res_taxo_seqs_asv <- anti_join(filtered_muck_IS_sig_res_taxo_seqs_asv, unfiltered_IS_tomcu_asvs, by = "ESV.ID")

#Filter out rows in the tomcu IS object that are also in the estuary and muck IS objects#
filtered_tomcu_IS_sig_res_taxo_seqs_asv <- anti_join(tomcu_IS_sig_res_taxo_seqs_asv, unfiltered_IS_estuary_asvs, by = "ESV.ID")
filtered_tomcu_IS_sig_res_taxo_seqs_asv <- anti_join(filtered_tomcu_IS_sig_res_taxo_seqs_asv, unfiltered_IS_muck_asvs, by = "ESV.ID")

#Order rows by group#
filtered_estuary_IS_sig_res_taxo_seqs_asv <- arrange(filtered_estuary_IS_sig_res_taxo_seqs_asv, group)
filtered_muck_IS_sig_res_taxo_seqs_asv <- arrange(filtered_muck_IS_sig_res_taxo_seqs_asv, group)
filtered_tomcu_IS_sig_res_taxo_seqs_asv <- arrange(filtered_tomcu_IS_sig_res_taxo_seqs_asv, group)


#Make a csv file for each of the filtered tables#
write.csv(as.data.frame(filtered_estuary_IS_sig_res_taxo_seqs_asv), file="filtered_estuary_IS_sig_res_taxo_seqs_asv.csv")
write.csv(as.data.frame(filtered_muck_IS_sig_res_taxo_seqs_asv), file="filtered_muck_IS_sig_res_taxo_seqs_asv.csv")
write.csv(as.data.frame(filtered_tomcu_IS_sig_res_taxo_seqs_asv), file="filtered_tomcu_IS_sig_res_taxo_seqs_asv.csv")

#Determine the number of indicators in each category#
length(which(sample_data(filtered_estuary_IS_sig_res_taxo_seqs_asv)$group == "1")) #656 IRL
length(which(sample_data(filtered_estuary_IS_sig_res_taxo_seqs_asv)$group == "2")) #1542 SLE
length(which(sample_data(filtered_muck_IS_sig_res_taxo_seqs_asv)$group == "1")) #490 Muck
length(which(sample_data(filtered_muck_IS_sig_res_taxo_seqs_asv)$group == "2")) #44 Not
length(which(sample_data(filtered_tomcu_IS_sig_res_taxo_seqs_asv)$group == "1")) #112 High/High
length(which(sample_data(filtered_tomcu_IS_sig_res_taxo_seqs_asv)$group == "2")) #5 High/Low


#Create lists of filtered significant asvs for each category#
filtered_IS_estuary_asvs <- subset(filtered_estuary_IS_sig_res_taxo_seqs_asv, select=c(ESV.ID))
filtered_IS_muck_asvs <- subset(filtered_muck_IS_sig_res_taxo_seqs_asv, select=c(ESV.ID))
filtered_IS_tomcu_asvs <- subset(filtered_tomcu_IS_sig_res_taxo_seqs_asv, select=c(ESV.ID))

#Create phyloseq objects of the unfiltered indicators using character string version of significant asvs##
charac_unfiltered_IS_estuary_asvs <- unfiltered_IS_estuary_asvs[,"ESV.ID"]
charac_unfiltered_IS_estuary_asvs <- as.character(charac_unfiltered_IS_estuary_asvs)
UFISestuaryDataAasv <- prune_taxa(charac_unfiltered_IS_estuary_asvs, RFLWSS_Y1_DataAasv)

charac_unfiltered_IS_muck_asvs <- unfiltered_IS_muck_asvs[,"ESV.ID"]
charac_unfiltered_IS_muck_asvs <- as.character(charac_unfiltered_IS_muck_asvs)
UFISmuckDataAasv <- prune_taxa(charac_unfiltered_IS_muck_asvs, RFLWSS_Y1_DataAasv)

charac_unfiltered_IS_tomcu_asvs <- unfiltered_IS_tomcu_asvs[,"ESV.ID"]
charac_unfiltered_IS_tomcu_asvs <- as.character(charac_unfiltered_IS_tomcu_asvs)
UFIStomcuDataAasv <- prune_taxa(charac_unfiltered_IS_tomcu_asvs, RFLWSS_Y1_DataAasv)

#Create phyloseq objects of the filtered indicators using character string version of significant asvs##
charac_filtered_IS_estuary_asvs <- filtered_IS_estuary_asvs[,"ESV.ID"]
charac_filtered_IS_estuary_asvs <- as.character(charac_filtered_IS_estuary_asvs)
FISestuaryDataAasv <- prune_taxa(charac_filtered_IS_estuary_asvs, RFLWSS_Y1_DataAasv)

charac_filtered_IS_muck_asvs <- filtered_IS_muck_asvs[,"ESV.ID"]
charac_filtered_IS_muck_asvs <- as.character(charac_filtered_IS_muck_asvs)
FISmuckDataAasv <- prune_taxa(charac_filtered_IS_muck_asvs, RFLWSS_Y1_DataAasv)

charac_filtered_IS_tomcu_asvs <- filtered_IS_tomcu_asvs[,"ESV.ID"]
charac_filtered_IS_tomcu_asvs <- as.character(charac_filtered_IS_tomcu_asvs)
FIStomcuDataAasv <- prune_taxa(charac_filtered_IS_tomcu_asvs, RFLWSS_Y1_DataAasv)

##Create heatmaps to see overall trends, may take a while to load##

unfiltered_IS_estuary_asv_heatmap <- plot_heatmap(UFDestuaryDataAasv, taxa.order = charac_unfiltered_IS_estuary_asvs, sample.order= "Estuary", sample.label = "Estuary")
unfiltered_IS_estuary_asv_heatmap

unfiltered_IS_muck_asv_heatmap <- plot_heatmap(UFDmuckDataAasv, taxa.order = charac_unfiltered_IS_muck_asvs, sample.order= "Muck", sample.label = "Muck")
unfiltered_IS_muck_asv_heatmap

unfiltered_IS_tomcu_asv_heatmap <- plot_heatmap(UFDtomcuDataAasv, taxa.order = charac_unfiltered_IS_tomcu_asvs, sample.order= "LOI.Cu", sample.label = "LOI.Cu")
unfiltered_IS_tomcu_asv_heatmap

filtered_IS_estuary_asv_heatmap <- plot_heatmap(FDestuaryDataAasv, taxa.order = charac_filtered_IS_estuary_asvs, sample.order= "Estuary", sample.label = "Estuary")
filtered_IS_estuary_asv_heatmap

filtered_IS_muck_asv_heatmap <- plot_heatmap(FDmuckDataAasv, taxa.order = charac_filtered_IS_muck_asvs, sample.order= "Muck", sample.label = "Muck")
filtered_IS_muck_asv_heatmap

filtered_IS_tomcu_asv_heatmap <- plot_heatmap(FDtomcuDataAasv, taxa.order = charac_filtered_IS_tomcu_asvs, sample.order= "LOI.Cu", sample.label = "LOI.Cu")
filtered_IS_tomcu_asv_heatmap

###MAKE COMBO INDICATOR LISTS###

##Make "Combo" indicator lists from overlapping unfiltered DESeq2 and IS indicators##

#Keep rows in the estuary DESEQ object that are also in the IRLSLE IS object#
unfiltered_estuary_combo_sig_res_taxo_seqs_asv <- merge(estuary_IS_sig_res_asv,estuary_DESeq2_sig_res_taxo_seqs_asv, by="ESV.ID")

#Keep rows in the muck DESEQ object that are also in the Muck IS object#
unfiltered_muck_combo_sig_res_taxo_seqs_asv <- merge(muck_IS_sig_res_asv,muck_DESeq2_sig_res_taxo_seqs_asv, by="ESV.ID")

#Keep rows in the tomcu DESEQ object that are also in the LOICu IS object#
unfiltered_tomcu_combo_sig_res_taxo_seqs_asv <- merge(tomcu_IS_sig_res_asv,tomcu_DESeq2_sig_res_taxo_seqs_asv, by="ESV.ID")

##Make "Combo" indicator lists from overlapping filtered DESeq2 and IS indicators##

#Create focused filtered IS tables without taxonomy or sequences#
filtered_estuary_IS_sig_res_asv <- subset(filtered_estuary_IS_sig_res_taxo_seqs_asv, select=c(ESV.ID,A,B,stat,p.value,group,padj))

filtered_muck_IS_sig_res_asv <- subset(filtered_muck_IS_sig_res_taxo_seqs_asv, select=c(ESV.ID,A,B,stat,p.value,group,padj))

filtered_tomcu_IS_sig_res_asv <- subset(filtered_tomcu_IS_sig_res_taxo_seqs_asv, select=c(ESV.ID,A,B,stat,p.value,group,padj))

#Keep rows in the estuary DESEQ object that are also in the IRLSLE IS object#
filtered_estuary_combo_sig_res_taxo_seqs_asv <- merge(filtered_estuary_IS_sig_res_asv,filtered_estuary_DESeq2_sig_res_taxo_seqs_asv, by="ESV.ID")

#Keep rows in the muck DESEQ object that are also in the Muck IS object#
filtered_muck_combo_sig_res_taxo_seqs_asv <- merge(filtered_muck_IS_sig_res_asv,filtered_muck_DESeq2_sig_res_taxo_seqs_asv, by="ESV.ID")

#Keep rows in the tomcu DESEQ object that are also in the LOICu IS object#
filtered_tomcu_combo_sig_res_taxo_seqs_asv <- merge(filtered_tomcu_IS_sig_res_asv,filtered_tomcu_DESeq2_sig_res_taxo_seqs_asv, by="ESV.ID")

#Determine the number of indicators in each subcategory of the unfiltered Combo tables#
length(which(sample_data(unfiltered_estuary_combo_sig_res_taxo_seqs_asv)$group == "1")) #1297 IRL
length(which(sample_data(unfiltered_estuary_combo_sig_res_taxo_seqs_asv)$group == "2")) #764 SLE
length(which(sample_data(unfiltered_muck_combo_sig_res_taxo_seqs_asv)$group == "1")) #141 Muck
length(which(sample_data(unfiltered_muck_combo_sig_res_taxo_seqs_asv)$group == "2")) #488 Not
length(which(sample_data(unfiltered_tomcu_combo_sig_res_taxo_seqs_asv)$group == "1")) #304 High/High
length(which(sample_data(unfiltered_tomcu_combo_sig_res_taxo_seqs_asv)$group == "2")) #253 High/Low

#Determine the number of indicators in each subcategory of the filtered Combo tables#
length(which(sample_data(filtered_estuary_combo_sig_res_taxo_seqs_asv)$group == "1")) #415 IRL
length(which(sample_data(filtered_estuary_combo_sig_res_taxo_seqs_asv)$group == "2")) #290 SLE
length(which(sample_data(filtered_muck_combo_sig_res_taxo_seqs_asv)$group == "1")) #30 Muck
length(which(sample_data(filtered_muck_combo_sig_res_taxo_seqs_asv)$group == "2")) #14 Not
length(which(sample_data(filtered_tomcu_combo_sig_res_taxo_seqs_asv)$group == "1")) #29 High/High
length(which(sample_data(filtered_tomcu_combo_sig_res_taxo_seqs_asv)$group == "2")) #5 High/Low


#Order unfiltered Combo table rows by group#
unfiltered_estuary_combo_sig_res_taxo_seqs_asv <- arrange(unfiltered_estuary_combo_sig_res_taxo_seqs_asv, group)
unfiltered_muck_combo_sig_res_taxo_seqs_asv <- arrange(unfiltered_muck_combo_sig_res_taxo_seqs_asv, group)
unfiltered_tomcu_combo_sig_res_taxo_seqs_asv <- arrange(unfiltered_tomcu_combo_sig_res_taxo_seqs_asv, group)

#Order filtered Combo table rows by group#
filtered_estuary_combo_sig_res_taxo_seqs_asv <- arrange(filtered_estuary_combo_sig_res_taxo_seqs_asv, group)
filtered_muck_combo_sig_res_taxo_seqs_asv <- arrange(filtered_muck_combo_sig_res_taxo_seqs_asv, group)
filtered_tomcu_combo_sig_res_taxo_seqs_asv <- arrange(filtered_tomcu_combo_sig_res_taxo_seqs_asv, group)

#Make a csv file for each of the unfiltered Combo tables#
write.csv(as.data.frame(unfiltered_estuary_combo_sig_res_taxo_seqs_asv), file="unfiltered_estuary_combo_sig_res_taxo_seqs_asv.csv")
write.csv(as.data.frame(unfiltered_muck_combo_sig_res_taxo_seqs_asv), file="unfiltered_muck_combo_sig_res_taxo_seqs_asv.csv")
write.csv(as.data.frame(unfiltered_tomcu_combo_sig_res_taxo_seqs_asv), file="unfiltered_tomcu_combo_sig_res_taxo_seqs_asv.csv")

#Make a csv file for each of the filtered Combo tables#
write.csv(as.data.frame(filtered_estuary_combo_sig_res_taxo_seqs_asv), file="filtered_estuary_combo_sig_res_taxo_seqs_asv.csv")
write.csv(as.data.frame(filtered_muck_combo_sig_res_taxo_seqs_asv), file="filtered_muck_combo_sig_res_taxo_seqs_asv.csv")
write.csv(as.data.frame(filtered_tomcu_combo_sig_res_taxo_seqs_asv), file="filtered_tomcu_combo_sig_res_taxo_seqs_asv.csv")

#Create lists of unfiltered significant asvs for each category#
unfiltered_Combo_estuary_asvs <- subset(unfiltered_estuary_combo_sig_res_taxo_seqs_asv, select=c(ESV.ID))
unfiltered_Combo_muck_asvs <- subset(unfiltered_muck_combo_sig_res_taxo_seqs_asv, select=c(ESV.ID))
unfiltered_Combo_tomcu_asvs <- subset(unfiltered_tomcu_combo_sig_res_taxo_seqs_asv, select=c(ESV.ID))

#Create lists of filtered significant asvs for each category#
filtered_Combo_estuary_asvs <- subset(filtered_estuary_combo_sig_res_taxo_seqs_asv, select=c(ESV.ID))
filtered_Combo_muck_asvs <- subset(filtered_muck_combo_sig_res_taxo_seqs_asv, select=c(ESV.ID))
filtered_Combo_tomcu_asvs <- subset(filtered_tomcu_combo_sig_res_taxo_seqs_asv, select=c(ESV.ID))

##Create phyloseq objects for Combo indicators##

#Create phyloseq objects of the unfiltered indicators using character string version of significant asvs##
charac_unfiltered_Combo_estuary_asvs <- unfiltered_Combo_estuary_asvs[,"ESV.ID"]
charac_unfiltered_Combo_estuary_asvs <- as.character(charac_unfiltered_Combo_estuary_asvs)
UFCestuaryDataAasv <- prune_taxa(charac_unfiltered_Combo_estuary_asvs, RFLWSS_Y1_DataAasv)

charac_unfiltered_Combo_muck_asvs <- unfiltered_Combo_muck_asvs[,"ESV.ID"]
charac_unfiltered_Combo_muck_asvs <- as.character(charac_unfiltered_Combo_muck_asvs)
UFCmuckDataAasv <- prune_taxa(charac_unfiltered_Combo_muck_asvs, RFLWSS_Y1_DataAasv)

charac_unfiltered_Combo_tomcu_asvs <- unfiltered_Combo_tomcu_asvs[,"ESV.ID"]
charac_unfiltered_Combo_tomcu_asvs <- as.character(charac_unfiltered_Combo_tomcu_asvs)
UFCtomcuDataAasv <- prune_taxa(charac_unfiltered_Combo_tomcu_asvs, RFLWSS_Y1_DataAasv)

#Create phyloseq objects of the filtered indicators using character string version of significant asvs##
charac_filtered_Combo_estuary_asvs <- filtered_Combo_estuary_asvs[,"ESV.ID"]
charac_filtered_Combo_estuary_asvs <- as.character(charac_filtered_Combo_estuary_asvs)
FCestuaryDataAasv <- prune_taxa(charac_filtered_Combo_estuary_asvs, RFLWSS_Y1_DataAasv)

charac_filtered_Combo_muck_asvs <- filtered_Combo_muck_asvs[,"ESV.ID"]
charac_filtered_Combo_muck_asvs <- as.character(charac_filtered_Combo_muck_asvs)
FCmuckDataAasv <- prune_taxa(charac_filtered_Combo_muck_asvs, RFLWSS_Y1_DataAasv)

charac_filtered_Combo_tomcu_asvs <- filtered_Combo_tomcu_asvs[,"ESV.ID"]
charac_filtered_Combo_tomcu_asvs <- as.character(charac_filtered_Combo_tomcu_asvs)
FCtomcuDataAasv <- prune_taxa(charac_filtered_Combo_tomcu_asvs, RFLWSS_Y1_DataAasv)

##Create heatmaps to see overall trends, may take a while to load##

unfiltered_Combo_estuary_asv_heatmap <- plot_heatmap(UFDestuaryDataAasv, taxa.order = charac_unfiltered_Combo_estuary_asvs, sample.order= "Estuary", sample.label = "Estuary")
unfiltered_Combo_estuary_asv_heatmap

unfiltered_Combo_muck_asv_heatmap <- plot_heatmap(UFDmuckDataAasv, taxa.order = charac_unfiltered_Combo_muck_asvs, sample.order= "Muck", sample.label = "Muck")
unfiltered_Combo_muck_asv_heatmap

unfiltered_Combo_tomcu_asv_heatmap <- plot_heatmap(UFDtomcuDataAasv, taxa.order = charac_unfiltered_Combo_tomcu_asvs, sample.order= "LOI.Cu", sample.label = "LOI.Cu")
unfiltered_Combo_tomcu_asv_heatmap

filtered_Combo_estuary_asv_heatmap <- plot_heatmap(FDestuaryDataAasv, taxa.order = charac_filtered_Combo_estuary_asvs, sample.order= "Estuary", sample.label = "Estuary")
filtered_Combo_estuary_asv_heatmap

filtered_Combo_muck_asv_heatmap <- plot_heatmap(FDmuckDataAasv, taxa.order = charac_filtered_Combo_muck_asvs, sample.order= "Muck", sample.label = "Muck")
filtered_Combo_muck_asv_heatmap

filtered_Combo_tomcu_asv_heatmap <- plot_heatmap(FDtomcuDataAasv, taxa.order = charac_filtered_Combo_tomcu_asvs, sample.order= "LOI.Cu", sample.label = "LOI.Cu")
filtered_Combo_tomcu_asv_heatmap

###INDICATOR EFFECTIVENESS METHOD###

##NOTES ON METHOD##
#"affected' = subcategory that is more affected by the stressor you are      concerned about. For example SLE for the Estuary category becasue it has     more freshwater discharges, 3MC for the Muck category because it has more    muck characteristics, HiHi for the TOM/Cu category becasue it has more       copper#

#"non-affected" = IRL, 0MC, HiLo#

#microbially-predicted "affected" sample = samples in the partitioning around medoids (PAM) cluster with the most metadata-defined "affected" samples#

#microbially-predicted "non-affected" sample = samples in the PAM cluster with the most metadata-defined "non-affected" samples#

#Overall Idea: Use PAM clustering to split samples between the "affected" and "non-affected" clusters based upon the indicators (microbially-predicted) and see how well this classification matches the metadata-defined classification by using the product of specificity and sensitivity#

#Copy results to an Excel Sheet for further statistical testing, each combo of factors should have four percentage types (Product is Sensitivity x Specificity), for example:#
#Indicator_Test	Taxonomic_Level	Filtering_Status	Metadata_Tested	Percent_Type	Percentage#
#Indicspecies   ASV             Unfiltered        Estuary         Total         98.78#
#Indicspecies   ASV             Unfiltered        Estuary         Sensitivity   96.23#
#Indicspecies   ASV             Unfiltered        Estuary         Specificity   99.15#
#Indicspecies   ASV             Unfiltered        Estuary         Product       95.41#


###ORIGINAL DATASETS
###Split samples of orginal datasets into two clusters to see how it will compare to splitting based upon indicators###

##Test on Muck and Estuary categories at asv level##

#Remove any samples with no asvs from target phyloseq#
nsamples(RFLWSS_Y1_DataAasv) #90
RFLWSS_Y1_DataAasv = prune_samples(sample_sums(RFLWSS_Y1_DataAasv)>=1, RFLWSS_Y1_DataAasv)
nsamples(RFLWSS_Y1_DataAasv) #90

#Make a matrix with the sample_data and relabel subcategories to make metadata-defined column#
Original_Est_Muck_asv <- as.matrix(sample_data(RFLWSS_Y1_DataAasv))
Original_Est_Muck_asv[Original_Est_Muck_asv ==  "Muck"] <- 2
Original_Est_Muck_asv[Original_Est_Muck_asv ==  "Not"] <- 1
Original_Est_Muck_asv[Original_Est_Muck_asv ==  "Mucky"] <- 2
Original_Est_Muck_asv[Original_Est_Muck_asv ==  "Muckish"] <- 2
Original_Est_Muck_asv[Original_Est_Muck_asv ==  "IRL"] <- 1
Original_Est_Muck_asv[Original_Est_Muck_asv ==  "SLE"] <- 2
Original_Est_Muck_asv<-as.data.frame(Original_Est_Muck_asv)


#Create Bray-Curtis distance matrix from asv table#
bc <- vegdist(t(otu_table(RFLWSS_Y1_DataAasv)))

#Conduct pam clustering based upon 2 clusters#
pam.res <- pam(bc,2, diss=TRUE)

#Add a microbially-predicted column (Orginal_asv) to matrix above for comparison to metadata-defined column#
Original_Est_Muck_asv <- cbind(Original_Est_Muck_asv, Orginal_asv = pam.res$cluster)

##Test Estuary category##

#Test  to see how well the metadata-defined and microbially-predicted columns match one another aka Total Efficiency; if the number is below 0.5, then switch the numbers in the metadata-defiend column#

#Test for Estuary Effeciency# 0.9888889
sum(Original_Est_Muck_asv$Estuary == Original_Est_Muck_asv$Orginal_asv)/nrow(Original_Est_Muck_asv)

##Test indicator sensitivity aka the true postives (TP)/ (TP + false negatives (FN))##
#TP = metadata-defined "affected" sample correctly placed in the microbiall   y-predicted "affected" cluster#
#FN = metadata-defined "affected" sample incorrectly placed in the           microbially-predicted "non-affected" cluster#
#TP+FN = all metadata-defined "affected" samples

#Test for SLE Effeciency aka Sensitvity# 1
Orginal_SLE_asv <- filter(Original_Est_Muck_asv, Estuary == "2")
sum(Orginal_SLE_asv$Estuary == Orginal_SLE_asv$Orginal_asv)/nrow(Orginal_SLE_asv)

##Test indicator specificity aka the true negatives (TN)/ (TN + false positives (FP))##
#TP = metadata-defined "affected" sample correctly placed in the microbiall   y-predicted "affected" cluster#
#FN = metadata-defined "affected" sample incorrectly placed in the           microbially-predicted "non-affected" cluster#
#TP+FN = all metadata-defined "affected" samples

#Test for IRL Effeciency aka Specificity# 0.9848485
Orginal_IRL_asv <- filter(Original_Est_Muck_asv, Estuary == "1")
sum(Orginal_IRL_asv$Estuary == Orginal_IRL_asv$Orginal_asv)/nrow(Orginal_IRL_asv)

#PCoA with Color by Estuary#
TotalPCoA <- ordinate(RFLWSS_Y1_DataAasv,"PCoA")
p = plot_ordination(RFLWSS_Y1_DataAasv, TotalPCoA, color="Estuary") + 
  ggtitle("Total asvs separated by Estuary Category PCoA with Bray-Curtis distance") +
  geom_text(aes(label=Original_Est_Muck_asv$Orginal_asv),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p

##Test Muck category##

#Test for Muck Total Efficiency# 0.7
sum(Original_Est_Muck_asv$Muck == Original_Est_Muck_asv$Orginal_asv)/nrow(Original_Est_Muck_asv)

#Test for 3MC Effeciency aka Sensitivity# 0.4666667
Orginal_3MC_asv <- filter(Original_Est_Muck_asv, Muck == "2")
sum(Orginal_3MC_asv$Muck == Orginal_3MC_asv$Orginal_asv)/nrow(Orginal_3MC_asv)

#Test for 0MC Effeciency aka Specificity# 0.8166667
Orginal_0MC_asv <- filter(Original_Est_Muck_asv, Muck == "1")
sum(Orginal_0MC_asv$Muck == Orginal_0MC_asv$Orginal_asv)/nrow(Orginal_0MC_asv)

#PCoA with Color by Muck#
TotalPCoA <- ordinate(RFLWSS_Y1_DataAasv,"PCoA")
p = plot_ordination(RFLWSS_Y1_DataAasv, TotalPCoA, color="Muck") + 
  ggtitle("Total asvs separated by Muck Category PCoA with Bray-Curtis distance") +
  geom_text(aes(label=Original_Est_Muck_asv$Orginal_asv),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p + annotate("text", x = 0.3, y = 0.3, label = c("SN = 0.44, SP=0.84"))+ scale_color_discrete(name="Muck Characteristics", limits=c("Muck", "Mucky", "Muckish", "Not"), labels=c("3", "2", "1", "0"))

##Test TOM/Cu category##

#Remove any samples with no asvs from target phyloseq#
nsamples(RFLWSS_Y1_HHvsHL_DataAasv) #29
RFLWSS_Y1_HHvsHL_DataAasv = prune_samples(sample_sums(RFLWSS_Y1_HHvsHL_DataAasv)>=1, RFLWSS_Y1_HHvsHL_DataAasv)
nsamples(RFLWSS_Y1_HHvsHL_DataAasv) #29

#Make a matrix with the sample_data and relabel subcategories to make metadata-defined column#
Original_TOMCu_asv <- as.matrix(sample_data(RFLWSS_Y1_HHvsHL_DataAasv))
Original_TOMCu_asv[Original_TOMCu_asv ==  "High-High"] <- 1
Original_TOMCu_asv[Original_TOMCu_asv ==  "High-Low"] <- 2
Original_TOMCu_asv<-as.data.frame(Original_TOMCu_asv)

#Create Bray-Curtis distance matrix from asv table#
bc <- vegdist(t(otu_table(RFLWSS_Y1_HHvsHL_DataAasv)))

#Conduct PAM clustering based upon 2 clusters#
pam.res <- pam(bc,2, diss=TRUE)

#Add a microbially-predicted column to matrix above for comparison to metadata-defined column#
Original_TOMCu_asv <- cbind(Original_TOMCu_asv, Orginal_asv = pam.res$cluster)

#Determine Total TOM/Cu Efficiency# 0.8275862
sum(Original_TOMCu_asv$LOI.Cu == Original_TOMCu_asv$Orginal_asv)/nrow(Original_TOMCu_asv)

#Determine HiHi Efficiency  aka Sensitivity# 1
Original_TOMCu_asvHH <- filter(Original_TOMCu_asv, LOI.Cu == "1") 
sum(Original_TOMCu_asvHH$LOI.Cu == Original_TOMCu_asvHH$Orginal_asv)/nrow(Original_TOMCu_asvHH)

#Determine HiLo Efficiency aka Specificity# 0.6578947
Original_TOMCu_asvHL <- filter(Original_TOMCu_asv, LOI.Cu == "2") 
sum(Original_TOMCu_asvHL$LOI.Cu == Original_TOMCu_asvHL$Orginal_asv)/nrow(Original_TOMCu_asvHL)

#PCoA with Color by TOM/Cu#
TotalPCoA <- ordinate(RFLWSS_Y1_HHvsHL_DataAasv,"PCoA")
p = plot_ordination(RFLWSS_Y1_HHvsHL_DataAasv, TotalPCoA, color="LOI.Cu") + 
  ggtitle("Total asvs separated by LOI-Cu Category PCoA with Bray-Curtis distance") +
  geom_text(aes(label=Original_TOMCu_asv$Orginal_asv),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p

###TEST INDICSPECIES ESTUARY INDICATORS###

##Test unfiltered Estuary asvs##

#Remove any samples with no asvs from target phyloseq#
nsamples(UFISestuaryDataAasv) #90
UFISestuaryDataAasv = prune_samples(sample_sums(UFISestuaryDataAasv)>=1, UFISestuaryDataAasv)
nsamples(UFISestuaryDataAasv) #90

#Make a matrix with the sample_data and relabel subcategories to make metadata-defined column#
UF_IS_Est_asv <- as.matrix(sample_data(UFISestuaryDataAasv))
UF_IS_Est_asv[UF_IS_Est_asv ==  "SLE"] <- 2
UF_IS_Est_asv[UF_IS_Est_asv ==  "IRL"] <- 1
UF_IS_Est_asv<-as.data.frame(UF_IS_Est_asv)

#Create Bray-Curtis distance matrix from asv table#
bc <- vegdist(t(otu_table(UFISestuaryDataAasv)))

#Conduct PAM clustering based upon 2 clusters#
pam.res <- pam(bc,2, diss=TRUE)

#Add a microbially-predicted column to matrix above for comparison to metadata-defined column#
UF_IS_Est_asv <- cbind(UF_IS_Est_asv, MP_UFISEstasv = pam.res$cluster)

#Test Total Estuary Efficiency# 1
sum(UF_IS_Est_asv$Estuary == UF_IS_Est_asv$MP_UFISEstasv)/nrow(UF_IS_Est_asv)

#Test SLE Efficiency aka Sensitvity# 1
UF_IS_SLE_asv <- filter(UF_IS_Est_asv, Estuary == "2")
sum(UF_IS_SLE_asv$Estuary == UF_IS_SLE_asv$MP_UFISEstasv)/nrow(UF_IS_SLE_asv)

#Test IRL Efficiency aka Specificity# 1
UF_IS_IRL_asv <- filter(UF_IS_Est_asv, Estuary == "1")
sum(UF_IS_IRL_asv$Estuary == UF_IS_IRL_asv$MP_UFISEstasv)/nrow(UF_IS_IRL_asv)

TotalPCoA <- ordinate(UFISestuaryDataAasv,"PCoA")
p = plot_ordination(UFISestuaryDataAasv, TotalPCoA, color="Estuary") + 
  ggtitle("Unfiltered IS Estuary asv Indicators PCoA with Bray-Curtis distance") +
  geom_text(aes(label=UF_IS_Est_asv$MP_UFISEstasv),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p


##Test filtered Estuary asvs##

#Remove any samples with no asvs from target phyloseq#
nsamples(FISestuaryDataAasv) #90
FISestuaryDataAasv = prune_samples(sample_sums(FISestuaryDataAasv)>=1, FISestuaryDataAasv)
nsamples(FISestuaryDataAasv) #90

#Make a matrix with the sample_data and relabel subcategories to make metadata-defined column#
F_IS_Est_asv <- as.matrix(sample_data(FISestuaryDataAasv))
F_IS_Est_asv[F_IS_Est_asv ==  "SLE"] <- 2
F_IS_Est_asv[F_IS_Est_asv ==  "IRL"] <- 1
F_IS_Est_asv<-as.data.frame(F_IS_Est_asv)


#Create Bray-Curtis distance matrix from asv table#
bc <- vegdist(t(otu_table(FISestuaryDataAasv)))

#Conduct PAM clustering based upon 2 clusters#
pam.res <- pam(bc,2, diss=TRUE)

#Add a microbially-predicted column to matrix above for comparison to metadata-defined column#
F_IS_Est_asv <- cbind(F_IS_Est_asv, MP_FISEstasv = pam.res$cluster)

#Test Total Estuary Efficiency# 0.9666667
sum(F_IS_Est_asv$Estuary == F_IS_Est_asv$MP_FISEstasv)/nrow(F_IS_Est_asv)
#Test SLE Efficiency aka Sensitivity# 0.875
F_IS_SLE_asv <- filter(F_IS_Est_asv, Estuary == "2")
sum(F_IS_SLE_asv$Estuary == F_IS_SLE_asv$MP_FISEstasv)/nrow(F_IS_SLE_asv)
#Test IRL Efficiency aka Specificity# 0.9807692
F_IS_IRL_asv <- filter(F_IS_Est_asv, Estuary == "1")
sum(F_IS_IRL_asv$Estuary == F_IS_IRL_asv$MP_FISEstasv)/nrow(F_IS_IRL_asv)

TotalPCoA <- ordinate(FISestuaryDataAasv,"PCoA")
p = plot_ordination(FISestuaryDataAasv, TotalPCoA, color="Estuary") + 
  ggtitle("Filtered IS Estuary asv Indicators PCoA with Bray-Curtis distance") +
  geom_text(aes(label=F_IS_Est_asv$MP_FISEstasv),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p


###TEST DESEQ2 ESTUARY INDICATORS###

##Test unfiltered Estuary asvs##

#Remove any samples with no asvs from target phyloseq#
nsamples(UFDestuaryDataAasv) #90
UFDestuaryDataAasv = prune_samples(sample_sums(UFDestuaryDataAasv)>=1, UFDestuaryDataAasv)
nsamples(UFDestuaryDataAasv) #90

#Make a matrix with the sample_data and relabel subcategories to make metadata-defined column#
UF_D_Est_asv <- as.matrix(sample_data(UFDestuaryDataAasv))
UF_D_Est_asv[UF_D_Est_asv ==  "SLE"] <- 2
UF_D_Est_asv[UF_D_Est_asv ==  "IRL"] <- 1
UF_D_Est_asv<-as.data.frame(UF_D_Est_asv)


#Create Bray-Curtis distance matrix from asv table#
bc <- vegdist(t(otu_table(UFDestuaryDataAasv)))

#Conduct PAM clustering based upon 2 clusters#
pam.res <- pam(bc,2, diss=TRUE)

#Add a microbially-predicted column to matrix above for comparison to metadata-defined column#
UF_D_Est_asv <- cbind(UF_D_Est_asv, MP_UFDEstasv = pam.res$cluster)

#Test Total Estuary Efficiency# 1
sum(UF_D_Est_asv$Estuary == UF_D_Est_asv$MP_UFDEstasv)/nrow(UF_D_Est_asv)

#Test SLE Efficiency aka Sensitivity# 1
UF_D_SLE_asv <- filter(UF_D_Est_asv, Estuary == "2")
sum(UF_D_SLE_asv$Estuary == UF_D_SLE_asv$MP_UFDEstasv)/nrow(UF_D_SLE_asv)

#Test IRL Efficiency aka Specificity# 1
UF_D_IRL_asv <- filter(UF_D_Est_asv, Estuary == "1")
sum(UF_D_IRL_asv$Estuary == UF_D_IRL_asv$MP_UFDEstasv)/nrow(UF_D_IRL_asv)


TotalPCoA <- ordinate(UFDestuaryDataAasv,"PCoA")
p = plot_ordination(UFDestuaryDataAasv, TotalPCoA, color="Estuary") + 
  ggtitle("Unfiltered DESeq2 Estuary asv Indicators PCoA with Bray-Curtis distance") +
  geom_text(aes(label=UF_D_Est_asv$MP_UFDEstasv),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p


#Test filtered DESeq2 Estuary asvs#

#Remove any samples with no asvs from target phyloseq#
nsamples(FDestuaryDataAasv) #90
FDestuaryDataAasv = prune_samples(sample_sums(FDestuaryDataAasv)>=1, FDestuaryDataAasv)
nsamples(FDestuaryDataAasv) #90

#Make a matrix with the sample_data and relabel subcategories to make metadata-defined column#
F_D_Est_asv <- as.matrix(sample_data(FDestuaryDataAasv))
F_D_Est_asv[F_D_Est_asv ==  "SLE"] <- 2
F_D_Est_asv[F_D_Est_asv ==  "IRL"] <- 1
F_D_Est_asv<-as.data.frame(F_D_Est_asv)


#Create Bray-Curtis distance matrix from asv table#
bc <- vegdist(t(otu_table(FDestuaryDataAasv)))

#Conduct PAM clustering based upon 2 clusters#
pam.res <- pam(bc,2, diss=TRUE)

#Add a microbially-predicted column to matrix above for comparison to metadata-defined column#
F_D_Est_asv <- cbind(F_D_Est_asv, MP_FDEstasv = pam.res$cluster)

#Test Total Estuary Efficiency# 0.9666667
sum(F_D_Est_asv$Estuary == F_D_Est_asv$MP_FDEstasv)/nrow(F_D_Est_asv)
#Test SLE Efficiency aka Sensitivity# 0.875
F_D_SLE_asv <- filter(F_D_Est_asv, Estuary == "2")
sum(F_D_SLE_asv$Estuary == F_D_SLE_asv$MP_FDEstasv)/nrow(F_D_SLE_asv)
#Test IRL Efficiency aka Specificity# 0.9679487
F_D_IRL_asv <- filter(F_D_Est_asv, Estuary == "1")
sum(F_D_IRL_asv$Estuary == F_D_IRL_asv$MP_FDEstasv)/nrow(F_D_IRL_asv)

TotalPCoA <- ordinate(FDestuaryDataAasv,"PCoA")
p = plot_ordination(FDestuaryDataAasv, TotalPCoA, color="Estuary") + 
  ggtitle("Filtered DESeq2 Estuary asv Indicators PCoA with Bray-Curtis distance") +
  geom_text(aes(label=F_D_Est_asv$MP_FDEstasv),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p

###TEST COMBO ESTUARY INDICATORS###

##Test unfiltered Estuary asvs##

#Remove any samples with no asvs from target phyloseq#
nsamples(UFCestuaryDataAasv) #90
UFCestuaryDataAasv = prune_samples(sample_sums(UFCestuaryDataAasv)>=1, UFCestuaryDataAasv)
nsamples(UFCestuaryDataAasv) #90

#Make a matrix with the sample_data and relabel subcategories to make metadata-defined column#
UF_C_Est_asv <- as.matrix(sample_data(UFCestuaryDataAasv))
UF_C_Est_asv[UF_C_Est_asv ==  "SLE"] <- 2
UF_C_Est_asv[UF_C_Est_asv ==  "IRL"] <- 1
UF_C_Est_asv<-as.data.frame(UF_C_Est_asv)


#Create Bray-Curtis distance matrix from asv table#
bc <- vegdist(t(otu_table(UFCestuaryDataAasv)))

#Conduct PAM clustering based upon 2 clusters#
pam.res <- pam(bc,2, diss=TRUE)

#Add a microbially-predicted column to matrix above for comparison to metadata-defined column#
UF_C_Est_asv <- cbind(UF_C_Est_asv, MP_UFCEstasv = pam.res$cluster)

#Test Total Estuary Efficiency# 1
sum(UF_C_Est_asv$Estuary == UF_C_Est_asv$MP_UFCEstasv)/nrow(UF_C_Est_asv)

#Test SLE Efficiency aka Sensitivity# 1
UF_C_SLE_asv <- filter(UF_C_Est_asv, Estuary == "2")
sum(UF_C_SLE_asv$Estuary == UF_C_SLE_asv$MP_UFCEstasv)/nrow(UF_C_SLE_asv)

#Test IRL Efficiency aka Specificity# 1
UF_C_IRL_asv <- filter(UF_C_Est_asv, Estuary == "1")
sum(UF_C_IRL_asv$Estuary == UF_C_IRL_asv$MP_UFCEstasv)/nrow(UF_C_IRL_asv)

TotalPCoA <- ordinate(UFCestuaryDataAasv,"PCoA")
p = plot_ordination(UFCestuaryDataAasv, TotalPCoA, color="Estuary") + 
  ggtitle("Unfiltered Final Estuary asv Indicators PCoA with Bray-Curtis distance") +
  geom_text(aes(label=UF_C_Est_asv$MP_UFCEstasv),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p


##Test filtered Estuary asvs##

#Remove any samples with no asvs from target phyloseq#
nsamples(FCestuaryDataAasv) #90
FCestuaryDataAasv = prune_samples(sample_sums(FCestuaryDataAasv)>=1, FCestuaryDataAasv)
nsamples(FCestuaryDataAasv) #90

#Make a matrix with the sample_data and relabel subcategories to make metadata-defined column#
F_C_Est_asv <- as.matrix(sample_data(FCestuaryDataAasv))
F_C_Est_asv[F_C_Est_asv ==  "SLE"] <- 2
F_C_Est_asv[F_C_Est_asv ==  "IRL"] <- 1
F_C_Est_asv<-as.data.frame(F_C_Est_asv)


#Create Bray-Curtis distance matrix from asv table#
bc <- vegdist(t(otu_table(FCestuaryDataAasv)))

#Conduct PAM clustering based upon 2 clusters#
pam.res <- pam(bc,2, diss=TRUE)

#Add a microbially-predicted column to matrix above for comparison to metadata-defined column#
F_C_Est_asv <- cbind(F_C_Est_asv, MP_FCEstasv = pam.res$cluster)

#Test Total Estuary Efficiency# 0.9777778
sum(F_C_Est_asv$Estuary == F_C_Est_asv$MP_FCEstasv)/nrow(F_C_Est_asv)

#Test SLE Efficiency aka Sensitivity# 0.9166667
F_C_SLE_asv <- filter(F_C_Est_asv, Estuary == "2")
sum(F_C_SLE_asv$Estuary == F_C_SLE_asv$MP_FCEstasv)/nrow(F_C_SLE_asv)

#Test IRL Efficiency aka Sensitivity# 0.9807692
F_C_IRL_asv <- filter(F_C_Est_asv, Estuary == "1")
sum(F_C_IRL_asv$Estuary == F_C_IRL_asv$MP_FCEstasv)/nrow(F_C_IRL_asv)

TotalPCoA <- ordinate(FCestuaryDataAasv,"PCoA")
p = plot_ordination(FCestuaryDataAasv, TotalPCoA, color="Estuary") + 
  ggtitle("Filtered Final Estuary asv Indicators PCoA with Bray-Curtis distance") +
  geom_text(aes(label=F_C_Est_asv$MP_FCEstasv),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p


###TEST INDICSPECIES MUCK INDICATORS###

##Test unfiltered Muck asvs##

#Remove any samples with no asvs from target phyloseq#
nsamples(UFISmuckDataAasv) #90
UFISmuckDataAasv = prune_samples(sample_sums(UFISmuckDataAasv)>=1, UFISmuckDataAasv)
nsamples(UFISmuckDataAasv) #90

#Make a matrix with the sample_data and relabel subcategories to make metadata-defined column#
UF_IS_Muck_asv <- as.matrix(sample_data(UFISmuckDataAasv))
UF_IS_Muck_asv[UF_IS_Muck_asv ==  "Muck"] <- 2
UF_IS_Muck_asv[UF_IS_Muck_asv ==  "Not"] <- 1
UF_IS_Muck_asv[UF_IS_Muck_asv ==  "Mucky"] <- 2
UF_IS_Muck_asv[UF_IS_Muck_asv ==  "Muckish"] <- 2
UF_IS_Muck_asv<-as.data.frame(UF_IS_Muck_asv)


#Create Bray-Curtis distance matrix from asv table#
bc <- vegdist(t(otu_table(UFISmuckDataAasv)))

#Conduct PAM clustering based upon 2 clusters#
pam.res <- pam(bc,2, diss=TRUE)

#Add a microbially-predicted column to matrix above for comparison to metadata-defined column#
UF_IS_Muck_asv <- cbind(UF_IS_Muck_asv, MP_UFISMuckasv = pam.res$cluster)

#Test Total Muck Efficiency# 0.8
sum(UF_IS_Muck_asv$Muck == UF_IS_Muck_asv$MP_UFISMuckasv)/nrow(UF_IS_Muck_asv)

#Test 3MC Efficiency aka Sensitivity# 0.8
UF_IS_3MC_asv <- filter(UF_IS_Muck_asv, Muck == "2")
sum(UF_IS_3MC_asv$Muck == UF_IS_3MC_asv$MP_UFISMuckasv)/nrow(UF_IS_3MC_asv)

#Test 0MC Efficiency aka Specificity# 0.8
UF_IS_0MC_asv <- filter(UF_IS_Muck_asv, Muck == "1")
sum(UF_IS_0MC_asv$Muck == UF_IS_0MC_asv$MP_UFISMuckasv)/nrow(UF_IS_0MC_asv)


TotalPCoA <- ordinate(UFISmuckDataAasv,"PCoA")
p = plot_ordination(UFISmuckDataAasv, TotalPCoA, color="Muck") + 
  ggtitle("Unfiltered IS Muck asv Indicators PCoA with Bray-Curtis distance") +
  geom_text(aes(label=UF_IS_Muck_asv$MP_UFISMuckasv),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p


##Test filtered Muck asvs##

#Remove any samples with no asvs from target phyloseq#
nsamples(FISmuckDataAasv) #90
FISmuckDataAasv = prune_samples(sample_sums(FISmuckDataAasv)>=1, FISmuckDataAasv)
nsamples(FISmuckDataAasv) #90

#Make a matrix with the sample_data and relabel subcategories to make metadata-defined column#
F_IS_Muck_asv <- as.matrix(sample_data(FISmuckDataAasv))
F_IS_Muck_asv[F_IS_Muck_asv ==  "Muck"] <- 2
F_IS_Muck_asv[F_IS_Muck_asv ==  "Not"] <- 1
F_IS_Muck_asv[F_IS_Muck_asv ==  "Mucky"] <- 2
F_IS_Muck_asv[F_IS_Muck_asv ==  "Muckish"] <- 2
F_IS_Muck_asv<-as.data.frame(F_IS_Muck_asv)


#Create Bray-Curtis distance matrix from asv table#
bc <- vegdist(t(otu_table(FISmuckDataAasv)))

#Conduct PAM clustering based upon 2 clusters#
pam.res <- pam(bc,2, diss=TRUE)

#Add a microbially-predicted column to matrix above for comparison to metadata-defined column#
F_IS_Muck_asv <- cbind(F_IS_Muck_asv, MP_FISMuckasv = pam.res$cluster)

#Test Total Muck Efficiency# 0.9111111
sum(F_IS_Muck_asv$Muck == F_IS_Muck_asv$MP_FISMuckasv)/nrow(F_IS_Muck_asv)
#Test 3MC Efficiency aka  Sensitivity# 1
F_IS_3MC_asv <- filter(F_IS_Muck_asv, Muck == "2")
sum(F_IS_3MC_asv$Muck == F_IS_3MC_asv$MP_FISMuckasv)/nrow(F_IS_3MC_asv)
#Test 0MC Efficiency aka Specificity# 0.8666667
F_IS_0MC_asv <- filter(F_IS_Muck_asv, Muck == "1")
sum(F_IS_0MC_asv$Muck == F_IS_0MC_asv$MP_FISMuckasv)/nrow(F_IS_0MC_asv)

TotalPCoA <- ordinate(FISmuckDataAasv,"PCoA")
p = plot_ordination(FISmuckDataAasv, TotalPCoA, color="Muck") + 
  ggtitle("Filtered IS Muck asv Indicators PCoA with Bray-Curtis distance") +
  geom_text(aes(label=F_IS_Muck_asv$MP_FISMuckasv),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p

###TEST DESEQ2 MUCK INDICATORS###

#Test unfiltered Muck asvs##

#Remove any samples with no asvs from target phyloseq#
nsamples(UFDmuckDataAasv) #90
UFDmuckDataAasv = prune_samples(sample_sums(UFDmuckDataAasv)>=1, UFDmuckDataAasv)
nsamples(UFDmuckDataAasv) #90

#Make a matrix with the sample_data and relabel subcategories to make metadata-defined column#
UF_D_Muck_asv <- as.matrix(sample_data(UFDmuckDataAasv))
UF_D_Muck_asv[UF_D_Muck_asv ==  "Muck"] <- 2
UF_D_Muck_asv[UF_D_Muck_asv ==  "Not"] <- 1
UF_D_Muck_asv[UF_D_Muck_asv ==  "Mucky"] <- 2
UF_D_Muck_asv[UF_D_Muck_asv ==  "Muckish"] <- 2
UF_D_Muck_asv<-as.data.frame(UF_D_Muck_asv)


#Create Bray-Curtis distance matrix from asv table#
bc <- vegdist(t(otu_table(UFDmuckDataAasv)))

#Conduct PAM clustering based upon 2 clusters#
pam.res <- pam(bc,2, diss=TRUE)

#Add a microbially-predicted column to matrix above for comparison to metadata-defined column#
UF_D_Muck_asv <- cbind(UF_D_Muck_asv, MP_UFDMuckasv = pam.res$cluster)

#Test Total Muck Efficiency# 0.7333333
sum(UF_D_Muck_asv$Muck == UF_D_Muck_asv$MP_UFDMuckasv)/nrow(UF_D_Muck_asv)

#Test 3MC Efficiency aka Sensitivity# 0.6
UF_D_3MC_asv <- filter(UF_D_Muck_asv, Muck == "2")
sum(UF_D_3MC_asv$Muck == UF_D_3MC_asv$MP_UFDMuckasv)/nrow(UF_D_3MC_asv)

#Test 0MC Efficiency aka Sensitivity# 0.8
UF_D_0MC_asv <- filter(UF_D_Muck_asv, Muck == "1")
sum(UF_D_0MC_asv$Muck == UF_D_0MC_asv$MP_UFDMuckasv)/nrow(UF_D_0MC_asv)

TotalPCoA <- ordinate(UFDmuckDataAasv,"PCoA")
p = plot_ordination(UFDmuckDataAasv, TotalPCoA, color="Muck") + 
  ggtitle("Unfiltered DESeq2 Muck asv Indicators PCoA with Bray-Curtis distance") +
  geom_text(aes(label=UF_D_Muck_asv$MP_UFDMuckasv),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p


##Test filtered Muck asvs##

#Remove any samples with no asvs from target phyloseq#
nsamples(FDmuckDataAasv) #90
FDmuckDataAasv = prune_samples(sample_sums(FDmuckDataAasv)>=1, FDmuckDataAasv)
nsamples(FDmuckDataAasv) #89

#Make a matrix with the sample_data and relabel subcategories to make metadata-defined column#
F_D_Muck_asv <- as.matrix(sample_data(FDmuckDataAasv))
F_D_Muck_asv[F_D_Muck_asv ==  "Muck"] <- 2
F_D_Muck_asv[F_D_Muck_asv ==  "Not"] <- 1
F_D_Muck_asv[F_D_Muck_asv ==  "Mucky"] <- 2
F_D_Muck_asv[F_D_Muck_asv ==  "Muckish"] <- 2
F_D_Muck_asv<-as.data.frame(F_D_Muck_asv)


#Create Bray-Curtis distance matrix from asv table#
bc <- vegdist(t(otu_table(FDmuckDataAasv)))

#Conduct PAM clustering based upon 2 clusters#
pam.res <- pam(bc,2, diss=TRUE)

#Add a microbially-predicted column to matrix above for comparison to metadata-defined column#
F_D_Muck_asv <- cbind(F_D_Muck_asv, MP_FDMuckasv = pam.res$cluster)

#Test Total Muck Efficiency# 0.8764045
sum(F_D_Muck_asv$Muck == F_D_Muck_asv$MP_FDMuckasv)/nrow(F_D_Muck_asv)

#Test 3MC Efficiency aka Sensitivy# 1
F_D_3MC_asv <- filter(F_D_Muck_asv, Muck == "2")
sum(F_D_3MC_asv$Muck == F_D_3MC_asv$MP_FDMuckasv)/nrow(F_D_3MC_asv)

#Test 0MC Efficiency aka Specificity# 0.8135593
F_D_0MC_asv <- filter(F_D_Muck_asv, Muck == "1")
sum(F_D_0MC_asv$Muck == F_D_0MC_asv$MP_FDMuckasv)/nrow(F_D_0MC_asv)

TotalPCoA <- ordinate(FDmuckDataAasv,"PCoA")
p = plot_ordination(FDmuckDataAasv, TotalPCoA, color="Muck") + 
  ggtitle("Filtered DESeq2 Muck asv Indicators PCoA with Bray-Curtis distance") +
  geom_text(aes(label=F_D_Muck_asv$MP_FDMuckasv),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p

###TEST COMBO MUCK INDICATORS###

##Test unfiltered Muck asvs##

#Remove any samples with no asvs from target phyloseq#
nsamples(UFCmuckDataAasv) #90
UFCmuckDataAasv = prune_samples(sample_sums(UFCmuckDataAasv)>=1, UFCmuckDataAasv)
nsamples(UFCmuckDataAasv) #90

#Make a matrix with the sample_data and relabel subcategories to make metadata-defined column#
UF_C_Muck_asvs <- as.matrix(sample_data(UFCmuckDataAasv))
UF_C_Muck_asvs[UF_C_Muck_asvs ==  "Muck"] <- 2
UF_C_Muck_asvs[UF_C_Muck_asvs ==  "Not"] <- 1
UF_C_Muck_asvs[UF_C_Muck_asvs ==  "Mucky"] <- 2
UF_C_Muck_asvs[UF_C_Muck_asvs ==  "Muckish"] <- 2
UF_C_Muck_asvs<-as.data.frame(UF_C_Muck_asvs)


#Create Bray-Curtis distance matrix from asv table#
bc <- vegdist(t(otu_table(UFCmuckDataAasv)))

#Conduct PAM clustering based upon 2 clusters#
pam.res <- pam(bc,2, diss=TRUE)

#Add a microbially-predicted column to matrix above for comparison to metadata-defined column#
UF_C_Muck_asvs <- cbind(UF_C_Muck_asvs, MP_UFCMuckasv = pam.res$cluster)

#Test Total Muck Efficiency# 0.8111111
sum(UF_C_Muck_asvs$Muck == UF_C_Muck_asvs$MP_UFCMuckasv)/nrow(UF_C_Muck_asvs)

#Test Muck Efficiency aka Sensitivity# 0.8333333
F_C_3MC_asvs <- filter(UF_C_Muck_asvs, Muck == "2")
sum(F_C_3MC_asvs$Muck == F_C_3MC_asvs$MP_UFCMuckasv)/nrow(F_C_3MC_asvs)

#Test Muck Efficiency aka Specificity# 0.8
F_C_0MC_asvs <- filter(UF_C_Muck_asvs, Muck == "1")
sum(F_C_0MC_asvs$Muck == F_C_0MC_asvs$MP_UFCMuckasv)/nrow(F_C_0MC_asvs)

TotalPCoA <- ordinate(UFCmuckDataAasv,"PCoA")
p = plot_ordination(UFCmuckDataAasv, TotalPCoA, color="Muck") + 
  ggtitle("Unfiltered Final Muck asv Indicators PCoA with Bray-Curtis distance") +
  geom_text(aes(label=UF_C_Muck_asvs$MP_UFCMuckasv),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p

##Test filtered Muck asvs##

#Remove any samples with no asvs from target phyloseq#
nsamples(FCmuckDataAasv) #90
FCmuckDataAasv = prune_samples(sample_sums(FCmuckDataAasv)>=1, FCmuckDataAasv)
nsamples(FCmuckDataAasv) #88

#Make a matrix with the sample_data and relabel subcategories to make metadata-defined column#
F_C_Muck_asvs <- as.matrix(sample_data(FCmuckDataAasv))
F_C_Muck_asvs[F_C_Muck_asvs ==  "Muck"] <- 2
F_C_Muck_asvs[F_C_Muck_asvs ==  "Not"] <- 1
F_C_Muck_asvs[F_C_Muck_asvs ==  "Mucky"] <- 2
F_C_Muck_asvs[F_C_Muck_asvs ==  "Muckish"] <- 2
F_C_Muck_asvs<-as.data.frame(F_C_Muck_asvs)


#Create Bray-Curtis distance matrix from asv table#
bc <- vegdist(t(otu_table(FCmuckDataAasv)))

#Conduct PAM clustering based upon 2 clusters#
pam.res <- pam(bc,2, diss=TRUE)

#Add a microbially-predicted column to matrix above for comparison to metadata-defined column#
F_C_Muck_asvs <- cbind(F_C_Muck_asvs, MP_FCMuckasv = pam.res$cluster)

#Test Total Muck Efficiency# 0.5113636
sum(F_C_Muck_asvs$Muck == F_C_Muck_asvs$MP_FCMuckasv)/nrow(F_C_Muck_asvs)

#Test Muck Efficiency aka Sensitivity# 0.6
F_C_3MC_asv <- filter(F_C_Muck_asvs, Muck == "1")
sum(F_C_3MC_asv$Muck == F_C_3MC_asv$MP_FCMuckasv)/nrow(F_C_3MC_asv)

#Test Not Efficiency aka Specificity# 0.4
F_C_0MC_asv <- filter(F_C_Muck_asvs, Muck == "2")
sum(F_C_0MC_asv$Muck == F_C_0MC_asv$MP_FCMuckasv)/nrow(F_C_0MC_asv)

TotalPCoA <- ordinate(FCmuckDataAasv,"PCoA")
p = plot_ordination(FCmuckDataAasv, TotalPCoA, color="Muck") + 
  ggtitle("Filtered Final Muck asv Indicators PCoA with Bray-Curtis distance") +
  geom_text(aes(label=F_C_Muck_asvs$MP_FCMuckasv),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p

###TEST INDICSPECIES TOM/Cu INDICATORS###

##Test unfiltered TOM/Cu asvs##

#Create phyloseq object focused on HiHi and HiLo samples#
UFIShihiDataAasv = subset_samples(UFIStomcuDataAasv, LOI.Cu=="High-High")
UFIShiloDataAasv = subset_samples(UFIStomcuDataAasv, LOI.Cu=="High-Low")
UFISfoctomcuDataAasv = merge_phyloseq(UFIShihiDataAasv, UFIShiloDataAasv)

#Remove any samples with no asvs from target phyloseq#
nsamples(UFISfoctomcuDataAasv) #29
UFISfoctomcuDataAasv = prune_samples(sample_sums(UFISfoctomcuDataAasv)>=1, UFISfoctomcuDataAasv)
nsamples(UFISfoctomcuDataAasv) #29

#Make a matrix with the sample_data and relabel subcategories to make metadata-defined column#
UF_IS_TOMCu_asv <- as.matrix(sample_data(UFISfoctomcuDataAasv))
UF_IS_TOMCu_asv[UF_IS_TOMCu_asv ==  "High-High"] <- 1
UF_IS_TOMCu_asv[UF_IS_TOMCu_asv ==  "High-Low"] <- 2
UF_IS_TOMCu_asv<-as.data.frame(UF_IS_TOMCu_asv)

#Create Bray-Curtis distance matrix from asv table#
bc <- vegdist(t(otu_table(UFISfoctomcuDataAasv)))

#Conduct PAM clustering based upon 2 clusters#
pam.res <- pam(bc,2, diss=TRUE)

#Add a microbially-predicted column to matrix above for comparison to metadata-defined column#
UF_IS_TOMCu_asv <- cbind(UF_IS_TOMCu_asv, MP_UFISTOMCuasv = pam.res$cluster)

#Test Total TOM/Cu Efficiency# 0.8275862
sum(UF_IS_TOMCu_asv$LOI.Cu == UF_IS_TOMCu_asv$MP_UFISTOMCuasv)/nrow(UF_IS_TOMCu_asv)

#Test High-High Indicators aka Sensitivity# 1
UF_IS_HiHi_asv <- filter(UF_IS_TOMCu_asv, LOI.Cu == "1")
sum(UF_IS_HiHi_asv$LOI.Cu == UF_IS_HiHi_asv$MP_UFISTOMCuasv)/nrow(UF_IS_HiHi_asv)

#Test High-Low Indicators aka Specificity# 0.7222222
UF_IS_HiLo_asv <- filter(UF_IS_TOMCu_asv, LOI.Cu == "2")
sum(UF_IS_HiLo_asv$LOI.Cu == UF_IS_HiLo_asv$MP_UFISTOMCuasv)/nrow(UF_IS_HiLo_asv)

TotalPCoA <- ordinate(UFISfoctomcuDataAasv,"PCoA")
p = plot_ordination(UFISfoctomcuDataAasv, TotalPCoA, color="LOI.Cu") + 
  ggtitle("Unfiltered IS LOICu asv Indicators PCoA with Bray-Curtis distance") +
  geom_text(aes(label=UF_IS_TOMCu_asv$MP_UFISTOMCuasv),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p


##Test filtered TOM/Cu asvs##

#Create phyloseq object focused on HiHi and HiLo samples#
FIShihiDataAasv = subset_samples(FIStomcuDataAasv, LOI.Cu=="High-High")
FIShiloDataAasv = subset_samples(FIStomcuDataAasv, LOI.Cu=="High-Low")
FISfoctomcuDataAasv = merge_phyloseq(FIShihiDataAasv, FIShiloDataAasv)

#Remove any samples with no asvs from target phyloseq#
nsamples(FISfoctomcuDataAasv) #29
FISfoctomcuDataAasv = prune_samples(sample_sums(FISfoctomcuDataAasv)>=1, FISfoctomcuDataAasv)
nsamples(FISfoctomcuDataAasv) #29

#Make a matrix with the sample_data and relabel subcategories to make metadata-defined column#
F_IS_TOMCu_asv <- as.matrix(sample_data(FISfoctomcuDataAasv))
F_IS_TOMCu_asv[F_IS_TOMCu_asv ==  "High-High"] <- 1
F_IS_TOMCu_asv[F_IS_TOMCu_asv ==  "High-Low"] <- 2
F_IS_TOMCu_asv<-as.data.frame(F_IS_TOMCu_asv)


#Create Bray-Curtis distance matrix from asv table#
bc <- vegdist(t(otu_table(FISfoctomcuDataAasv)))

#Conduct PAM clustering based upon 2 clusters#
pam.res <- pam(bc,2, diss=TRUE)

#Add a microbially-predicted column to matrix above for comparison to metadata-defined column#
F_IS_TOMCu_asv <- cbind(F_IS_TOMCu_asv, MP_FISTOMCuasv = pam.res$cluster)

#Test Total TOM/Cu Efficiency# 0.6206897
sum(F_IS_TOMCu_asv$LOI.Cu == F_IS_TOMCu_asv$MP_FISTOMCuasv)/nrow(F_IS_TOMCu_asv)

#Test HiHi Efficiency aka Sensitivity# 1
F_IS_HiHi_asv <- filter(F_IS_TOMCu_asv, LOI.Cu == "1")
sum(F_IS_HiHi_asv$LOI.Cu == F_IS_HiHi_asv$MP_FISTOMCuasv)/nrow(F_IS_HiHi_asv)

#Test HiLo Efficiency aka Sensitivity# 0.3888889
F_IS_HiLo_asv <- filter(F_IS_TOMCu_asv, LOI.Cu == "2") 
sum(F_IS_HiLo_asv$LOI.Cu == F_IS_HiLo_asv$MP_FISTOMCuasv)/nrow(F_IS_HiLo_asv)

TotalPCoA <- ordinate(FISfoctomcuDataAasv,"PCoA")
p = plot_ordination(FISfoctomcuDataAasv, TotalPCoA, color="LOI.Cu") + 
  ggtitle("Filtered IS LOICu asv Indicators PCoA with Bray-Curtis distance") +
  geom_text(aes(label=F_IS_TOMCu_asv$MP_FISTOMCuasv),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p

##TEST DESEQ2 LOICu INDICATORS##

##Test unfiltered TOM/Cu asvs##

#Create phyloseq object focused on HiHi and HiLo samples#
UFDhihiDataAasv = subset_samples(UFDtomcuDataAasv, LOI.Cu=="High-High")
UFDhiloDataAasv = subset_samples(UFDtomcuDataAasv, LOI.Cu=="High-Low")
UFDfoctomcuDataAasv = merge_phyloseq(UFDhihiDataAasv, UFDhiloDataAasv)

#Remove any samples with no asvs from target phyloseq#
nsamples(UFDfoctomcuDataAasv) #29
UFDfoctomcuDataAasv = prune_samples(sample_sums(UFDfoctomcuDataAasv)>=1, UFDfoctomcuDataAasv)
nsamples(UFDfoctomcuDataAasv) #29

#Make a matrix with the sample_data and relabel subcategories to make metadata-defined column#
UF_D_TOMCu_asv <- as.matrix(sample_data(UFDfoctomcuDataAasv))
UF_D_TOMCu_asv[UF_D_TOMCu_asv ==  "High-High"] <- 1
UF_D_TOMCu_asv[UF_D_TOMCu_asv ==  "High-Low"] <- 2
UF_D_TOMCu_asv<-as.data.frame(UF_D_TOMCu_asv)


#Create Bray-Curtis distance matrix from asv table#
bc <- vegdist(t(otu_table(UFDfoctomcuDataAasv)))

#Conduct PAM clustering based upon 2 clusters#
pam.res <- pam(bc,2, diss=TRUE)

#Add a microbially-predicted column to matrix above for comparison to metadata-defined column#
UF_D_TOMCu_asv <- cbind(UF_D_TOMCu_asv, MP_UFDTOMCuasv = pam.res$cluster)

#Test Total TOM/Cu Efficiency# 0.8275862
sum(UF_D_TOMCu_asv$LOI.Cu == UF_D_TOMCu_asv$MP_UFDTOMCuasv)/nrow(UF_D_TOMCu_asv)

#Test HiHi Efficiency aka Sensitivity# 1
UF_D_HiHi_asv <- filter(UF_D_TOMCu_asv, LOI.Cu == "1")
sum(UF_D_HiHi_asv$LOI.Cu == UF_D_HiHi_asv$MP_UFDTOMCuasv)/nrow(UF_D_HiHi_asv)

#Test HiLo Efficiency aka Specificity# 0.7222222
UF_D_HiLo_asv <- filter(UF_D_TOMCu_asv, LOI.Cu == "2")
sum(UF_D_HiLo_asv$LOI.Cu == UF_D_HiLo_asv$MP_UFDTOMCuasv)/nrow(UF_D_HiLo_asv)

TotalPCoA <- ordinate(UFDfoctomcuDataAasv,"PCoA")
p = plot_ordination(UFDfoctomcuDataAasv, TotalPCoA, color="LOI.Cu") + 
  ggtitle("Unfiltered DESeq2 LOICu asv Indicators PCoA with Bray-Curtis distance") +
  geom_text(aes(label=UF_D_TOMCu_asv$MP_UFDTOMCuasv),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p


##Test filtered TOM/Cu asvs##

#Create phyloseq object focused on HiHi and HiLo samples#
FDhihiDataAasv = subset_samples(FDtomcuDataAasv, LOI.Cu=="High-High")
FDhiloDataAasv = subset_samples(FDtomcuDataAasv, LOI.Cu=="High-Low")
FDfoctomcuDataAasv = merge_phyloseq(FDhihiDataAasv, FDhiloDataAasv)

#Remove any samples with no asvs from target phyloseq#
nsamples(FDfoctomcuDataAasv) #29
FDfoctomcuDataAasv = prune_samples(sample_sums(FDfoctomcuDataAasv)>=1, FDfoctomcuDataAasv)
nsamples(FDfoctomcuDataAasv) #29


#Make a matrix with the sample_data and relabel subcategories to make metadata-defined column#
F_D_TOMCu_asv <- as.matrix(sample_data(FDfoctomcuDataAasv))
F_D_TOMCu_asv[F_D_TOMCu_asv ==  "High-High"] <- 1
F_D_TOMCu_asv[F_D_TOMCu_asv ==  "High-Low"] <- 2
F_D_TOMCu_asv<-as.data.frame(F_D_TOMCu_asv)


#Create Bray-Curtis distance matrix from asv table#
bc <- vegdist(t(otu_table(FDfoctomcuDataAasv)))

#Conduct PAM clustering based upon 2 clusters#
pam.res <- pam(bc,2, diss=TRUE)

#Add a microbially-predicted column to matrix above for comparison to metadata-defined column#
F_D_TOMCu_asv <- cbind(F_D_TOMCu_asv, MP_FDTOMCuasv = pam.res$cluster)

#Test Total TOM/Cu Efficiency# 0.862069
sum(F_D_TOMCu_asv$LOI.Cu == F_D_TOMCu_asv$MP_FDTOMCuasv)/nrow(F_D_TOMCu_asv)

#Test HiHi Efficiency aka Sensitivity# 0.6363636
F_D_HiHi_asv <- filter(F_D_TOMCu_asv, LOI.Cu == "1")
sum(F_D_HiHi_asv$LOI.Cu == F_D_HiHi_asv$MP_FDTOMCuasv)/nrow(F_D_HiHi_asv)

#Test HiLo Efficiency aka Specificity# 1
F_D_HiLo_asv <- filter(F_D_TOMCu_asv, LOI.Cu == "2")
sum(F_D_HiLo_asv$LOI.Cu == F_D_HiLo_asv$MP_FDTOMCuasv)/nrow(F_D_HiLo_asv)

TotalPCoA <- ordinate(FDfoctomcuDataAasv,"PCoA")
p = plot_ordination(FDfoctomcuDataAasv, TotalPCoA, color="LOI.Cu") + 
  ggtitle("Filtered DESeq2 LOICu asv Indicators PCoA with Bray-Curtis distance") +
  geom_text(aes(label=F_D_TOMCu_asv$MP_FDTOMCuasv),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p

##TEST COMBO TOM/CU INDICATORS##

##Test unfiltered TOM/Cu asvs##

#Create phyloseq object focused on HiHi and HiLo samples#
UFChihiDataAasv = subset_samples(UFCtomcuDataAasv, LOI.Cu=="High-High")
UFChiloDataAasv = subset_samples(UFCtomcuDataAasv, LOI.Cu=="High-Low")
UFCfoctomcuDataAasv = merge_phyloseq(UFChihiDataAasv, UFChiloDataAasv)

#Remove any samples with no asvs from target phyloseq#
nsamples(UFCfoctomcuDataAasv) #29
UFCfoctomcuDataAasv = prune_samples(sample_sums(UFCfoctomcuDataAasv)>=1, UFCfoctomcuDataAasv)
nsamples(UFCfoctomcuDataAasv) #29

#Make a matrix with the sample_data and relabel subcategories to make metadata-defined column#
UF_C_TOMCu_asv <- as.matrix(sample_data(UFCfoctomcuDataAasv))
UF_C_TOMCu_asv[UF_C_TOMCu_asv ==  "High-High"] <- 1
UF_C_TOMCu_asv[UF_C_TOMCu_asv ==  "High-Low"] <- 2
UF_C_TOMCu_asv<-as.data.frame(UF_C_TOMCu_asv)


#Create Bray-Curtis distance matrix from asv table#
bc <- vegdist(t(otu_table(UFCfoctomcuDataAasv)))

#Conduct PAM clustering based upon 2 clusters#
pam.res <- pam(bc,2, diss=TRUE)

#Add a microbially-predicted column to matrix above for comparison to metadata-defined column#
UF_C_TOMCu_asv <- cbind(UF_C_TOMCu_asv, MP_UFCTOMCuasv = pam.res$cluster)

#Test Total TOM/Cu Efficiency# 0.8275862
sum(UF_C_TOMCu_asv$LOI.Cu == UF_C_TOMCu_asv$MP_UFCTOMCuasv)/nrow(UF_C_TOMCu_asv)

#Test HiHi Efficiency aka Sensitivity# 1
UF_C_HiHi_asv <- filter(UF_C_TOMCu_asv, LOI.Cu == "1")
sum(UF_C_HiHi_asv$LOI.Cu == UF_C_HiHi_asv$MP_UFCTOMCuasv)/nrow(UF_C_HiHi_asv)

#Test HiLo Efficiency aka Specificity# 0.7222222
UF_C_HiLo_asv <- filter(UF_C_TOMCu_asv, LOI.Cu == "2")
sum(UF_C_HiLo_asv$LOI.Cu == UF_C_HiLo_asv$MP_UFCTOMCuasv)/nrow(UF_C_HiLo_asv)

TotalPCoA <- ordinate(UFCfoctomcuDataAasv,"PCoA")
p = plot_ordination(UFCfoctomcuDataAasv, TotalPCoA, color="LOI.Cu") + 
  ggtitle("Unfiltered Final LOICu asv Indicators PCoA with Bray-Curtis distance") +
  geom_text(aes(label=UF_C_TOMCu_asv$MP_UFCTOMCuasv),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p


##Test filtered TOM/Cu asvs##

#Create phyloseq object focused on HiHi and HiLo samples#
FChihiDataAasv = subset_samples(FCtomcuDataAasv, LOI.Cu=="High-High")
FChiloDataAasv = subset_samples(FCtomcuDataAasv, LOI.Cu=="High-Low")
FCfoctomcuDataAasv = merge_phyloseq(FChihiDataAasv, FChiloDataAasv)

#Remove any samples with no asvs from target phyloseq#
nsamples(FCfoctomcuDataAasv) #29
FCfoctomcuDataAasv = prune_samples(sample_sums(FCfoctomcuDataAasv)>=1, FCfoctomcuDataAasv)
nsamples(FCfoctomcuDataAasv) #29

#Make a matrix with the sample_data and relabel subcategories to make metadata-defined column#
F_C_TOMCu_asv <- as.matrix(sample_data(FCfoctomcuDataAasv))
F_C_TOMCu_asv[F_C_TOMCu_asv ==  "High-High"] <- 1
F_C_TOMCu_asv[F_C_TOMCu_asv ==  "High-Low"] <- 2
F_C_TOMCu_asv<-as.data.frame(F_C_TOMCu_asv)


#Create Bray-Curtis distance matrix from asv table#
bc <- vegdist(t(otu_table(FCfoctomcuDataAasv)))

#Conduct PAM clustering based upon 2 clusters#
pam.res <- pam(bc,2, diss=TRUE)

#Add a microbially-predicted column to matrix above for comparison to metadata-defined column#
F_C_TOMCu_asv <- cbind(F_C_TOMCu_asv, MP_FCTOMCuasv = pam.res$cluster)

#Test Total TOM/Cu Efficiency# 0.862069
sum(F_C_TOMCu_asv$LOI.Cu == F_C_TOMCu_asv$MP_FCTOMCuasv)/nrow(F_C_TOMCu_asv)

#Test HiHi Efficiency aka Sensitivity# 0.6363636
F_C_HiHi_asv <- filter(F_C_TOMCu_asv, LOI.Cu == "1")
sum(F_C_HiHi_asv$LOI.Cu == F_C_HiHi_asv$MP_FCTOMCuasv)/nrow(F_C_HiHi_asv)

#Test HiLo Efficiency aka Specificity# 0.9032258
F_C_HiLo_asv <- filter(F_C_TOMCu_asv, LOI.Cu == "2")
sum(F_C_HiLo_asv$LOI.Cu == F_C_HiLo_asv$MP_FCTOMCuasv)/nrow(F_C_HiLo_asv)


TotalPCoA <- ordinate(FCfoctomcuDataAasv,"PCoA")
p = plot_ordination(FCfoctomcuDataAasv, TotalPCoA, color="LOI.Cu") + 
  ggtitle("Filtered Final LOICu asv Indicators PCoA with Bray-Curtis distance") +
  geom_text(aes(label=F_C_TOMCu_asv$MP_FCTOMCuasv),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p

###RUNNING ABOVE SCRIPTS ON OTHER TAXONOMIC LEVELS###

##Upload QIIME2 information into phyloseq for taxnomic level you want to test##
species_aa_table = read.csv(file.choose(), row.names=1)
species_taxa_table = read.csv(file.choose(), row.names = 1)
species_taxa_table = as.matrix(species_taxa_table)

#Can use the same metadata table as ASV level#
DottedMETA = read.delim(file.choose(), row.names=1) 

dotted_meta = sample_data(DottedMETA)
Species_OTU = otu_table(species_aa_table, taxa_are_rows = TRUE)
Species_TAX = tax_table(species_taxa_table)

LWSSSpeciesA <- merge_phyloseq(Species_OTU, Species_TAX, dotted_meta)

#Filter rare organisms#
ntaxa(LWSSSpeciesA) #4120#
FLWSSSpeciesA <- filter_taxa(LWSSSpeciesA, function(x) sum(x) >0, TRUE)
ntaxa(FLWSSSpeciesA) #3428#
RLWSSSpeciesA <- filter_taxa(FLWSSSpeciesA, function(x) sum(x) >20, TRUE)
ntaxa(RLWSSSpeciesA) #1805#

###To run above scripts at any other level, just replace all lowercase "asv" with the level you are testing at to differentiate from asv level by (Ctrl + F all lowercase "asv" with "species"). ESV.ID is programmed into the results of DESeq2, thus it is used throughout despite different taxonomic levels###

##Example##
#Convert phyloseq to deseq2 object centered around a irlvssle factor#
estuary_cudds_asv = phyloseq_to_deseq2(RFLWSS_Y1_DataAasv,  ~ Estuary)

#Becomes...#

#Convert phyloseq to deseq2 object centered around a irlvssle factor#
estuary_cudds_species = phyloseq_to_deseq2(RLWSSSpeciesA,  ~ Estuary)


###TESTING PRODUCT TRENDS###

getwd() 
setwd("C:/Users/bradshawd/Documents/Journal_Articles/New_Indicators")
getwd() #"C:/Users/bradshawd/Documents/Journal_Articles/New_Indicators"

##Note: upload a table with all the total, sensitivity, specificity, and product percentages for each set of factors you want to test. In the manuscript there were 72 different combos: 3 indicator tests, 3 metadata categories, 4 taxonomic levels, 2 filtering statuses. There were also 12 additional combos associated with the orginal datasets that did not have any indicator tests done on them: 3 metadata categories and 4 taxonomic levels## 
Adjusted_Summary  <-read.csv(file.choose())

#Filter out original values to focus on tested ones#
Tested_Adjusted_Summary <-filter(Adjusted_Summary, Indicator_Test != "Original")

##Make box plots to visulize differences amoung percentage types by various factors##

#Percentage by Indicator Test#
ggplot(Adjusted_Summary, aes(x=Percent_Type, y=Percentage, color=Indicator_Test)) + 
  geom_boxplot()

#Percentage by Taxonomic Level#
ggplot(Tested_Adjusted_Summary, aes(x=Percent_Type, y=Percentage, color=Taxonomic_Level)) +   geom_boxplot()

#Percentage by Metadata Category#
ggplot(Tested_Adjusted_Summary, aes(x=Percent_Type, y=Percentage, color=Metadata_Tested)) +   geom_boxplot()

#Percentage by Filtering Status#
ggplot(Tested_Adjusted_Summary, aes(x=Percent_Type, y=Percentage, color=Filtering_Status)) +   geom_boxplot()

##Test the trends just associated with the Product##
#Filter out all metrics besides Product#
Product_Adjusted_Summary <- filter(Adjusted_Summary, Percent_Type == "Product")

#Filter out original values to focus on tested ones#
Product_Tested_Adjusted_Summary <- filter(Tested_Adjusted_Summary, Percent_Type == "Product")

#Test statistically main and pairwise#
kruskal.test(Percentage ~ Indicator_Test, data = Product_Adjusted_Summary)
dunnTest(Percentage ~ Indicator_Test, data = Product_Adjusted_Summary, method="bh")

kruskal.test(Percentage ~ Taxonomic_Level, data = Product_Tested_Adjusted_Summary)
dunnTest(Percentage ~ Taxonomic_Level, data = Product_Tested_Adjusted_Summary, method="bh")

kruskal.test(Percentage ~ Metadata_Tested, data = Product_Tested_Adjusted_Summary)
dunnTest(Percentage ~ Metadata_Tested, data = Product_Tested_Adjusted_Summary, method="bh")

kruskal.test(Percentage ~ Filtering_Status, data = Product_Tested_Adjusted_Summary)

#Make box plots above with just Product to visualize differences#
ggplot(Product_Adjusted_Summary, aes(x=Indicator_Test, y=Percentage, color=Indicator_Test)) +   geom_boxplot()+ 
  ylab("Product (%)") +
  scale_color_manual(name="Indicator Test", values = c("red", "forestgreen", "blue", "purple"))+
  #annotate("text", x = 1:4, y = 105, label = c("a", "b", "a", "b"))+
  xlab("Indicator Test")
#theme(axis.title.x=element_blank())

ggplot(Product_Tested_Adjusted_Summary, aes(x=Taxonomic_Level, y=Percentage, color=Taxonomic_Level)) + 
  geom_boxplot()+ 
  geom_boxplot() +
  xlab("Taxonomic Level") +
  ylab("Product") +
  scale_color_discrete(name="Taxonomic Level")+
  #annotate("text", x = 1:4, y = 105, label = c("a", "a", "a", "a"))+
  scale_x_discrete(limits=c("ESV", "Species", "Genus", "Family"), labels=c("ASV", "Species", "Genus", "Family"))

ggplot(Product_Tested_Adjusted_Summary, aes(x=Metadata_Tested, y=Percentage, color=Metadata_Tested)) + 
  geom_boxplot() +
  ylab("Product (%)") +
  scale_color_discrete(name="Metadata Tested", labels = c("Estuary", "TOM/Cu", "Muck"))+
  #annotate("text", x = 1:3, y = 105, label = c("a", "b", "b")) +
  #theme(axis.title.x=element_blank())+
  scale_x_discrete(name="Metadata Tested", labels = c("Estuary", "TOM/Cu", "Muck"))

ggplot(Product_Tested_Adjusted_Summary, aes(x=Indicator_Test, y=Percentage, color=Filtering_Status)) + 
  geom_boxplot()+ 
  ylab("Product (%)") +
  scale_color_discrete(name="Indicator Test")
  #annotate("text", x = 1:2, y = 105, label = c("a", "a")) +xlab("FIltering Status")
#theme(axis.title.x=element_blank())

##Conclusion: from these resutls it was decided to use the Unfiltered, "Combo", Genus level indicators for each metadata category moving forward. This reduced the number of indicator lists from 72 to 3##

###Affected Indicator Percentage (AIP)###
#A metric was designed to summarize all of the indicators associated with a category without taking into account taxonomy to allow uncultured/unknonwn/undefined indicators to still be used to assess changes in the entire indicator microbiome#
#AIP = sum of sequences associated with the affected indicators divided by the sum of the total indicator sequences#
#Estuary AIP = SLE sequences sum / (SLE sequences sum + IRL sequences sum)#
#Muck AIP = 3MC sequences sum / (3MC sequences sum + 0MC sequences sum)#
#TOM/Cu AIP = HiHi sequences sum / (HiHi sequences sum + HiLo sequences sum)#

###Create a table with all metadata, environmental varaibles, and AIP information###

##First, create a table with all metadata and environmental varibles information##

#Extract metadata table from the genus phyloseq object#
LWSS_Meta_V2 <- as(sample_data(RFLWSS_Y1_GenusAgenus), "data.frame")


##Create new columns in data table related to the Estuary AIP##

#Determine the sequence sums from Unfiltered, Combo, Genus Level IRL indicators and upload as a new column#
LWSS_Meta_V2$UF_Combo_Genus_IRL <- prune_taxa(as.character(unfiltered_estuary_combo_sig_res_taxo_seqs_genus$ESV.ID[unfiltered_estuary_combo_sig_res_taxo_seqs_genus$log2FoldChange<0]), RFLWSS_Y1_GenusAgenus) %>% sample_sums() 

#Determine the sequence sums from Unfiltered, Combo, Genus Level SLE indicators and upload as a new column#
LWSS_Meta_V2$UF_Combo_Genus_SLE <- prune_taxa(as.character(unfiltered_estuary_combo_sig_res_taxo_seqs_genus$ESV.ID[unfiltered_estuary_combo_sig_res_taxo_seqs_genus$log2FoldChange>0]), RFLWSS_Y1_GenusAgenus) %>% sample_sums() 

#Calculate Estuary AIP and add in new column#
LWSS_Meta_V2$UF_Combo_Genus_Estuary <- with(LWSS_Meta_V2, 100*UF_Combo_Genus_SLE/(UF_Combo_Genus_SLE+UF_Combo_Genus_IRL))

##Create new columns in data table related to the Muck AIP##

#Determine the sequence sums from Unfiltered, Combo, Genus Level 3MC indicators and upload as a new column#
LWSS_Meta_V2$UF_Combo_Genus_3MC <- prune_taxa(as.character(unfiltered_muck_combo_sig_res_taxo_seqs_genus$ESV.ID[unfiltered_muck_combo_sig_res_taxo_seqs_genus$log2FoldChange<0]), RFLWSS_Y1_GenusAgenus) %>% sample_sums() 

#Determine the sequence sums from Unfiltered, Combo, Genus Level 0MC indicators and upload as a new column#
LWSS_Meta_V2$UF_Combo_Genus_0MC <- prune_taxa(as.character(unfiltered_muck_combo_sig_res_taxo_seqs_genus$ESV.ID[unfiltered_muck_combo_sig_res_taxo_seqs_genus$log2FoldChange>0]), RFLWSS_Y1_GenusAgenus) %>% sample_sums() 

#Calculate Muck AIP and add in new column#
LWSS_Meta_V2$UF_Combo_Genus_Muck <- with(LWSS_Meta_V2, 100*UF_Combo_Genus_3MC/(UF_Combo_Genus_3MC+UF_Combo_Genus_0MC))

##Create new columns in data table related to the TOM/Cu AIP##

#Determine the sequence sums from Unfiltered, Combo, Genus Level HiHi indicators and upload as a new column#
LWSS_Meta_V2$UF_Combo_Genus_HiHi <- prune_taxa(as.character(unfiltered_tomcu_combo_sig_res_taxo_seqs_genus$ESV.ID[unfiltered_tomcu_combo_sig_res_taxo_seqs_genus$log2FoldChange<0]), RFLWSS_Y1_GenusAgenus) %>% sample_sums() 

#Determine the sequence sums from Unfiltered, Combo, Genus Level HiLo indicators and upload as a new column#
LWSS_Meta_V2$UF_Combo_Genus_HiLo <- prune_taxa(as.character(unfiltered_tomcu_combo_sig_res_taxo_seqs_genus$ESV.ID[unfiltered_tomcu_combo_sig_res_taxo_seqs_genus$log2FoldChange>0]), RFLWSS_Y1_GenusAgenus) %>% sample_sums() 

#Calculate TOM/Cu AIP and add in new column#
LWSS_Meta_V2$UF_Combo_Genus_TOMCu <- with(LWSS_Meta_V2, 100*UF_Combo_Genus_HiHi/(UF_Combo_Genus_HiHi+UF_Combo_Genus_HiLo))

#Upload csv with all of the environmental variables#
LWSS_Env_Y1 <- read.csv(file.choose())

colnames(LWSS_Env_Y1)

#Make the row names the Sample IDs
LWSS_Env_Y1 <- LWSS_Env_Y1 %>% remove_rownames %>% column_to_rownames(var="ï..Sample.Name")

#Combine the metadata and environmental data tables#
LWSS_Total_V2 <- merge(LWSS_Env_Y1, LWSS_Meta_V2, by="row.names")

#Make the rownames a separate column#
rownames(LWSS_Total_V2) <- LWSS_Total_V2[,1]

#Convert all occurences of NaN to NA#
LWSS_Total_V2[LWSS_Total_V2 ==  "NaN"] <- NA


###VISUALIZE AND TEST METRICS###

##Statistically test AIP differences between subcategories##

#Estuary AIP Stats#
wilcox.test(UF_Combo_Genus_Estuary ~ Estuary, data = LWSS_Total_V2)

#Muck AIP Stats#
kruskal.test(UF_Combo_Genus_Muck ~ Muck, data = LWSS_Total_V2)
dunnTest(UF_Combo_Genus_Muck ~ Muck, data = LWSS_Total_V2, method="bh")

#TOM/Cu AIP Stats#
kruskal.test(UF_Combo_Genus_TOMCu ~ LOI.Cu, data = LWSS_Total_V2)
dunnTest(UF_Combo_Genus_TOMCu ~ LOI.Cu, data = LWSS_Total_V2, method="bh")

##Create Box Plots to visualize differences between subcategories##

#Estuary AIP BP#
ggplot(data = LWSS_Total_V2, aes( x = Estuary, y = UF_Combo_Genus_Estuary, color=Estuary)) +
  geom_boxplot() +
  ylab("Affected Indicator Percentage") +
  ggtitle("Affected Indicator Percentage") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  annotate("text", x = 1:2, y = 105, label = c("a", "b")) +
  theme(axis.title.x=element_blank()) +
  #theme(legend.position="bottom")
  theme(legend.position = "none")

#Muck AIP BP#
ggplot(data = LWSS_Total_V2, aes( x = Muck, y = UF_Combo_Genus_Muck, color=Muck)) +
  geom_boxplot() +
  ylab("Affected Indicator Percentage") +
  xlab("Muck Characteristics") +
  ggtitle("Affected Indicator Percentage") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_discrete(name="Muck Characteristics", limits=c("Muck", "Mucky", "Muckish", "Not"), labels=c("3", "2", "1", "0")) +
  scale_x_discrete(limits=c("Muck", "Mucky", "Muckish", "Not"), labels=c("3", "2", "1", "0")) +
  annotate("text", x = 1:4, y = 105, label = c("a", "a", "a", "b")) +
  theme(legend.position="none")

##Create Metric Box Plots TOM/Cu Category to visualize stats##
#TOM/Cu AIP BP#
ggplot(data = LWSS_Total_V2, aes( x = LOI.Cu, y = UF_Combo_Genus_TOMCu, color=LOI.Cu)) +
  geom_boxplot() +
  ylab("Affected Indicator Percentage") +
  xlab("TOM/Cu Category") +
  #ggtitle("Affected Indicator Percentage") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_discrete(name="TOM/Cu", limits=c("High.High", "High.Low", "Low.High", "Low.Low"), labels=c("HiHi", "HiLo", "LoHi", "LoLo")) +
  scale_x_discrete(limits=c("High.High", "High.Low", "Low.High", "Low.Low"), labels=c("HiHi", "HiLo", "LoHi", "LoLo")) +
  annotate("text", x = 1:4, y = 105, label = c("a", "b", "ac", "c")) +
  #theme(axis.title.x=element_blank()) +
  #  theme(legend.position="bottom")
  theme(legend.position = "none")


##Create Box Plots to visualize differences between sites##

#Estuary AIP by Site#
ggplot(data = LWSS_Total_V2, aes( x = Site, y = UF_Combo_Genus_Estuary)) +
  geom_boxplot() +
  theme(axis.text.x  = element_text(angle=45, hjust=1))+
  ylab("Affected Indicator Percentage") +
  xlab("Site") +
  #ggtitle("Overall Sediment Samples Estuary AIP")+
  #theme(legend.position = "none") +
  #scale_x_discrete(labels=c("Barber Bridge", "Fort Pierce", "Harbor Branch Channel", "Harbortown Marina", "Hobe Sound", "Jensen Beach", "Jupiter Narrows", "Linkport", "Manatee Pocket", "Melbourne Causeway", "Merritt Island Causeway", "Middle Estuary", "North Fork", "Round Island", "Sebastian Inlet", "South Fork", "South Fork 2", "Vero Beach", "Vero Beach Marina"))+
  theme(plot.title = element_text(hjust = 0.5))

#Estuary AIP by Site by Season#
ggplot(data = LWSS_Total_V2, aes( x = Site, y = UF_Combo_Genus_Estuary, color=Season)) +
  geom_boxplot() +
  theme(axis.text.x  = element_text(angle=45, hjust=1))+
  ylab("Affected Indicator Percentage") +
  xlab("Site") +
  #scale_x_discrete(labels=c("Barber Bridge", "Fort Pierce", "Harbor Branch Channel", "Harbortown Marina", "Hobe Sound", "Jensen Beach", "Jupiter Narrows", "Linkport", "Manatee Pocket", "Melbourne Causeway", "Merritt Island Causeway", "Middle Estuary", "North Fork", "Round Island", "Sebastian Inlet", "South Fork", "South Fork 2", "Vero Beach", "Vero Beach Marina"))+
  scale_color_manual(values = c("darkgray", "black"))

#Muck AIP by Site#
ggplot(data = LWSS_Total_V2, aes( x = Site, y = UF_Combo_Genus_Muck)) +
  geom_boxplot()+
  ylab("Affected Indicator Percentage") +
  xlab("Site") +
  theme(axis.text.x  = element_text(angle=45, hjust=1)) +
  #annotate("text", x=c(1, 2, 5, 6, 7, 11, 13, 14, 15, 18, 19), y=105, label=c("0/6", "0/12","0/12", "0/12", "0/12",  "0/12", "0/12", "0/6", "0/12", "0/12", "0/6"), size=3) +
  theme(axis.text.x  = element_text(angle=45, hjust=1)) +  #scale_x_discrete(labels=c("Barber Bridge", "Fort Pierce", "Harbor Branch Channel", "Harbortown Marina", "Hobe Sound", "Jensen Beach", "Jupiter Narrows", "Linkport", "Manatee Pocket", "Melbourne Causeway", "Merritt Island Causeway", "Middle Estuary", "North Fork", "Round Island", "Sebastian Inlet", "South Fork", "South Fork 2", "Vero Beach", "Vero Beach Marina"), )+
  annotate("text", x=c(3, 4, 8, 9, 10, 12, 16, 17), y=105, label=c("17/18", "6/6", "2/12", "11/12", "9/12", "12/12","12/12", "1/12"), size=3, fontface=2)+
  scale_color_manual (name="Site", values=c("midnightblue", "blue", "red", "purple", "yellow", "violet", "maroon", "darkgreen", "dodgerblue", "gold",  "moccasin",  "darkgray", "powderblue", "green", "lightgray", "black", "slateblue", "greenyellow", "lightseagreen"), labels=c("Barber Bridge", "Fort Pierce", "Harbor Branch Channel", "Harbortown Marina", "Hobe Sound", "Jensen Beach", "Jupiter Narrows", "Linkport", "Manatee Pocket", "Melbourne Causeway", "Merritt Island Causeway", "Middle Estuary", "North Fork", "Round Island", "Sebastian Inlet", "South Fork", "South Fork 2", "Vero Beach", "Vero Beach Marina"))

#TOM/Cu AIP by Sites only with samples with >10% TOM#
ggplot(data = subset(LWSS_Total_V2, Total.Organic.Matter.... > 10), aes( x = Site, y = UF_Combo_Genus_TOMCu, color=LOI.Cu)) +
  geom_boxplot()+
  ylab("Affected Indicator Percentage") +
  xlab("Site") +
  theme(axis.text.x  = element_text(angle=45, hjust=1))+
  scale_color_manual(name="TOM/Cu", labels=c("HiHi", "HiLo"), values = c("darkgray", "black"))

##Adjust for multiple testing across all KW and MWU main tests##
indicator_stats_names2=c("IT", "TL", "MT", "FS", "LWS_M_AIP", "LWS_M_AR", "LWS_M_AS", "LWS_M_ATP", "LWS_L_AIP", "LWS_L_AR", "LWS_L_AS", "LWS_L_ATP", "LWS_E_AIP", "LWS_E_AR", "LWS_E_AS", "LWS_E_ATP", "DGS_M_AIP", "DGS_L_AIP", "ER_M_AIP", "ER_L", "CR_M", "CR_L", "WR_M", "WR_L", "NL_M", "NL_L", "SL_M", "SL_L", "HCS_M_AIP", "FP_M", "HT_M", "LP_M", "HB_M", "HCS_E_Survey", "HCS_E_Site", "FP_E", "HT_E", "LP_E", "HB_E", "DGS_E_Site", "HCS_L_AIP", "FP_L", "HB_L", "HT_L", "LP_L")
indicator_stats_names2=data.frame(indicator_stats_names2)

indicator_stats_p.value2=c(0.001321, 0.8389, 1.331e-06, 0.1952, 2.2e-16, 2.2e-16, 2.634e-15, 2.2771e-16, 1.517e-14, 1.517e-14, 6.653e-09, 1.163e-15, 2.2e-16, 2.2e-16, 2.2e-16, 2.2e-16, 8.665e-08, 3.913e-05, 0.228, 0.1718, 0.04987, 0.1574, 0.3761, 0.1349, 0.3611, 0.2001, 0.0396, 0.1129, 0.06879, 0.06081, 0.07939, 0.02732, 0.3012, 3.17e-06, 0.78, 0.02732, 0.02732, 0.2521, 0.05091, 9.196e-09, 0.1454, 0.06081, 0.2881, 0.02732, 0.4911)
indicator_stats_p.value2=data.frame(indicator_stats_p.value2)

indicator_stats2 <- cbind(indicator_stats_names2, indicator_stats_p.value2)

indicator_stats2$padj <- p.adjust(indicator_stats2$indicator_stats_p.value2, method = "BH")

##Adjust for multiple testing across all Kruskal-Wallis main tests. Here we are adjusting for all Kruskal-Wallis tests conducted in the manuscript##

#Create a list of all the names associated with each of the p-values you are adusting and make it a data frame#
indicator_stats_names=c("IT", "TL", "MT", "FS", "LWS_M_AIP", "LWS_M_AR", "LWS_M_AS", "LWS_M_ATP", "LWS_L_AIP", "LWS_L_AR", "LWS_L_AS", "LWS_L_ATP", "LWS_E_AIP", "LWS_E_AR", "LWS_E_AS", "LWS_E_ATP")
indicator_stats_names=data.frame(indicator_stats_names)

#Create a list of all the associated p-values and make it a dataframe#
indicator_stats_p.value=c(0.001321, 0.8389, 1.331e-06, 0.1952, 2.2e-16, 2.2e-16, 2.634e-15, 2.2771e-16, 1.517e-14, 1.517e-14, 6.653e-09, 1.163e-15, 2.2e-16, 2.2e-16, 2.2e-16, 2.2e-16)
indicator_stats_p.value=data.frame(indicator_stats_p.value)

#Bind together the dataframes by columns, make sure each dataframe is in the correct order#
indicator_stats <- cbind(indicator_stats_names, indicator_stats_p.value)

#Adjust the p-values using Benjamini-Hochberg#
indicator_stats$padj <- p.adjust(indicator_stats$indicator_stats_p.value, method = "BH")
indicator_stats


###OTHER USEFUL SCRIPTS FOR ANALZYING INDICATOR TRENDS###

##Breakdown of number of unknowns script##
#First step is to create datataxa which is copy of table with all the DESeq, IS, taxonomic, and otu tables so you can manipulate it without worring about loss of data#

#Second step is to remove all D_x__ labels from taxonomic columns to allow comparison between columns with lapply#

#Last step is to subtract two elements from one another to determine the number of undefineds at a particular taxonomic level#
#Element One: datataxa%>%filter(group==1)%>%nrow#
#Determines the total number of rows      associated with group 1#

#Element Two: datataxa%>%filter(group==1)%>% filter(!grepl(pattern = "uncultured", x =       Order.x)) %>%filter(!grepl(pattern = "unidentified", x = Order.x))%>% filter(!grepl          (pattern = "metagenome", x = Order.x))%>%filter(!grepl(pattern = "Unknown", x = Order.x))%    >%filter(as.character(Class.x) != as.character(Order.x))%>%filter(!grepl(pattern ="__", x    =Order.x))%>%nrow#
#Determines the number of rows that does NOT include orders with uncultured, unidentified,    metagenome, Unknonwn in its name or is a repeat of the previous taxonomic level#

#filter(!grepl(...)) removes rows with that pattern from the table#

##Get counts for unfiltered, combo, genus level Estuary indicators##

#Get abundance table for Order Level for Estuary Indicators#
Unfiltered_Combo_Genus_IRL_Orders <- as.data.frame(table(unfiltered_estuary_combo_sig_res_taxo_seqs_genus$Order.x, unfiltered_estuary_combo_sig_res_taxo_seqs_genus$log2FoldChange<0))

#Get abundance table for Family Level for Estuary Indicators#
Unfiltered_Combo_Genus_IRL_Familys <- as.data.frame(table(unfiltered_estuary_combo_sig_res_taxo_seqs_genus$Family.x, unfiltered_estuary_combo_sig_res_taxo_seqs_genus$log2FoldChange<0))

#Make data an alterable dataframe and remove all taxonomy headers for filtering #
datataxa <- unfiltered_estuary_combo_sig_res_taxo_seqs_genus
datataxa <- data.frame(lapply(datataxa, function(x) {gsub("D_0__", "", x)}))
datataxa <- data.frame(lapply(datataxa, function(x) {gsub("D_1__", "", x)}))
datataxa <- data.frame(lapply(datataxa, function(x) {gsub("D_2__", "", x)}))
datataxa <- data.frame(lapply(datataxa, function(x) {gsub("D_3__", "", x)}))
datataxa <- data.frame(lapply(datataxa, function(x) {gsub("D_4__", "", x)}))
datataxa <- data.frame(lapply(datataxa, function(x) {gsub("D_5__", "", x)}))
datataxa <- data.frame(lapply(datataxa, function(x) {gsub("D_6__", "", x)}))

##Determine Order level unknowns
#Determine number of unknowns in IRL indicators#
datataxa%>%filter(group==1)%>%nrow - datataxa%>%filter(group==1)%>% filter(!grepl(pattern = "uncultured", x = Order.x)) %>%filter(!grepl(pattern = "unidentified", x = Order.x))%>% filter(!grepl(pattern = "metagenome", x = Order.x))%>%filter(!grepl(pattern = "Unknown", x = Order.x))%>%filter(as.character(Class.x) != as.character(Order.x))%>%filter(!grepl(pattern ="__", x=Order.x))%>%nrow #15

#Determine number of unknowns in SLE indicators#
datataxa%>%filter(group==2)%>%nrow - datataxa%>%filter(group==2)%>% filter(!grepl(pattern = "uncultured", x = Order.x)) %>%filter(!grepl(pattern = "unidentified", x = Order.x))%>% filter(!grepl(pattern = "metagenome", x = Order.x))%>%filter(!grepl(pattern = "Unknown", x = Order.x))%>%filter(as.character(Class.x) != as.character(Order.x))%>%filter(!grepl(pattern ="__", x=Order.x))%>%nrow #30

##Determine Family level unknowns
#Determine number of unknowns in IRL indicators#
datataxa%>%filter(group==1)%>%nrow - datataxa%>%filter(group==1)%>% filter(!grepl(pattern = "uncultured", x = Family.x)) %>%filter(!grepl(pattern = "unidentified", x = Family.x))%>% filter(!grepl(pattern = "metagenome", x = Family.x))%>%filter(!grepl(pattern = "Unknown", x = Family.x))%>%filter(as.character(Order.x) != as.character(Family.x))%>%filter(!grepl(pattern ="__", x=Family.x))%>%nrow #34

#Determine number of unknowns in SLE indicators#
datataxa%>%filter(group==2)%>%nrow - datataxa%>%filter(group==2)%>% filter(!grepl(pattern = "uncultured", x = Family.x)) %>%filter(!grepl(pattern = "unidentified", x = Family.x))%>% filter(!grepl(pattern = "metagenome", x = Family.x))%>%filter(!grepl(pattern = "Unknown", x = Family.x))%>%filter(as.character(Order.x) != as.character(Family.x))%>%filter(!grepl(pattern ="__", x=Family.x))%>%nrow #76

##Determine Genus level unknowns
#Determine number of unknowns in IRL indicators#
datataxa%>%filter(group==1)%>%nrow - datataxa%>%filter(group==1)%>% filter(!grepl(pattern = "uncultured", x = Genus.x)) %>%filter(!grepl(pattern = "unidentified", x = Genus.x))%>% filter(!grepl(pattern = "metagenome", x = Genus.x))%>%filter(!grepl(pattern = "Unknown", x = Genus.x))%>%filter(as.character(Family.x) != as.character(Genus.x))%>%filter(!grepl(pattern ="__", x=Genus.x))%>%nrow #42

#Determine number of unknowns in SLE indicators#
datataxa%>%filter(group==2)%>%nrow - datataxa%>%filter(group==2)%>% filter(!grepl(pattern = "uncultured", x = Genus.x)) %>%filter(!grepl(pattern = "unidentified", x = Genus.x))%>% filter(!grepl(pattern = "metagenome", x = Genus.x))%>%filter(!grepl(pattern = "Unknown", x = Genus.x))%>%filter(as.character(Family.x) != as.character(Genus.x))%>%filter(!grepl(pattern ="__", x=Genus.x))%>%nrow #133


##Get counts for unfiltered, combo, genus level Muck indicators##
#Get abundance table for Order Level for Muck Indicators#
Unfiltered_Combo_Genus_Muck_Orders <- as.data.frame(table(unfiltered_muck_combo_sig_res_taxo_seqs_genus$Order.x, unfiltered_muck_combo_sig_res_taxo_seqs_genus$log2FoldChange<0))

#Get abundance table for Family Level for Muck Indicators#
Unfiltered_Combo_Genus_Muck_Familys <- as.data.frame(table(unfiltered_muck_combo_sig_res_taxo_seqs_genus$Family.x, unfiltered_muck_combo_sig_res_taxo_seqs_genus$log2FoldChange<0))

#Make data an alterable dataframe and remove all taxonomy headers for filtering #
datataxa <- unfiltered_muck_combo_sig_res_taxo_seqs_genus
datataxa <- data.frame(lapply(datataxa, function(x) {gsub("D_0__", "", x)}))
datataxa <- data.frame(lapply(datataxa, function(x) {gsub("D_1__", "", x)}))
datataxa <- data.frame(lapply(datataxa, function(x) {gsub("D_2__", "", x)}))
datataxa <- data.frame(lapply(datataxa, function(x) {gsub("D_3__", "", x)}))
datataxa <- data.frame(lapply(datataxa, function(x) {gsub("D_4__", "", x)}))
datataxa <- data.frame(lapply(datataxa, function(x) {gsub("D_5__", "", x)}))
datataxa <- data.frame(lapply(datataxa, function(x) {gsub("D_6__", "", x)}))

##Determine Order level unknowns
#Determine number of unknowns in 3MC indicators#
datataxa%>%filter(group==1)%>%nrow - datataxa%>%filter(group==1)%>% filter(!grepl(pattern = "uncultured", x = Order.x)) %>%filter(!grepl(pattern = "unidentified", x = Order.x))%>% filter(!grepl(pattern = "metagenome", x = Order.x))%>%filter(!grepl(pattern = "Unknown", x = Order.x))%>%filter(as.character(Class.x) != as.character(Order.x))%>%filter(!grepl(pattern ="__", x=Order.x))%>%nrow #13

#Determine number of unknowns in 0MC indicators#
datataxa2 <- datataxa%>%filter(group==2)%>%nrow - datataxa%>%filter(group==2)%>% filter(!grepl(pattern = "uncultured", x = Order.x)) %>%filter(!grepl(pattern = "unidentified", x = Order.x))%>% filter(!grepl(pattern = "metagenome", x = Order.x))%>%filter(!grepl(pattern = "Unknown", x = Order.x))%>%filter(as.character(Class.x) != as.character(Order.x))%>%filter(!grepl(pattern ="__", x=Order.x))%>%nrow #2

##Determine Family level unknowns
#Determine number of unknowns in 3MC indicators#
datataxa%>%filter(group==1)%>%nrow - datataxa%>%filter(group==1)%>% filter(!grepl(pattern = "uncultured", x = Family.x)) %>%filter(!grepl(pattern = "unidentified", x = Family.x))%>% filter(!grepl(pattern = "metagenome", x = Family.x))%>%filter(!grepl(pattern = "Unknown", x = Family.x))%>%filter(as.character(Order.x) != as.character(Family.x))%>%filter(!grepl(pattern ="__", x=Family.x))%>%nrow #29

#Determine number of unknowns in 0MC indicators#
datataxa%>%filter(group==2)%>%nrow - datataxa%>%filter(group==2)%>% filter(!grepl(pattern = "uncultured", x = Family.x)) %>%filter(!grepl(pattern = "unidentified", x = Family.x))%>% filter(!grepl(pattern = "metagenome", x = Family.x))%>%filter(!grepl(pattern = "Unknown", x = Family.x))%>%filter(as.character(Order.x) != as.character(Family.x))%>%filter(!grepl(pattern ="__", x=Family.x))%>%nrow #11

##Determine Genus level unknowns
#Determine number of unknowns in 3MC indicators#
datataxa%>%filter(group==1)%>%nrow - datataxa%>%filter(group==1)%>% filter(!grepl(pattern = "uncultured", x = Genus.x)) %>%filter(!grepl(pattern = "unidentified", x = Genus.x))%>% filter(!grepl(pattern = "metagenome", x = Genus.x))%>%filter(!grepl(pattern = "Unknown", x = Genus.x))%>%filter(as.character(Family.x) != as.character(Genus.x))%>%filter(!grepl(pattern ="__", x=Genus.x))%>%nrow #39

#Determine number of unknowns in 0MC indicators#
datataxa%>%filter(group==2)%>%nrow - datataxa%>%filter(group==2)%>% filter(!grepl(pattern = "uncultured", x = Genus.x)) %>%filter(!grepl(pattern = "unidentified", x = Genus.x))%>% filter(!grepl(pattern = "metagenome", x = Genus.x))%>%filter(!grepl(pattern = "Unknown", x = Genus.x))%>%filter(as.character(Family.x) != as.character(Genus.x))%>%filter(!grepl(pattern ="__", x=Genus.x))%>%nrow #24


##Get counts for unfiltered, combo, genus level TOM/Cu indicators##

#Get abundance table for Order Level for TOM/Cu Indicators#
Unfiltered_Combo_Genus_LOICu_Orders <- as.data.frame(table(unfiltered_tomcu_combo_sig_res_taxo_seqs_genus$Order.x, unfiltered_tomcu_combo_sig_res_taxo_seqs_genus$log2FoldChange<0))

#Get abundance table for Family Level for TOM/Cu Indicators#
Unfiltered_Combo_Genus_LOICu_Familys <- as.data.frame(table(unfiltered_tomcu_combo_sig_res_taxo_seqs_genus$Family.x, unfiltered_tomcu_combo_sig_res_taxo_seqs_genus$log2FoldChange<0))

#Make data an alterable dataframe and remove all taxonomy headers for filtering #
datataxa <- unfiltered_tomcu_combo_sig_res_taxo_seqs_genus
datataxa <- data.frame(lapply(datataxa, function(x) {gsub("D_0__", "", x)}))
datataxa <- data.frame(lapply(datataxa, function(x) {gsub("D_1__", "", x)}))
datataxa <- data.frame(lapply(datataxa, function(x) {gsub("D_2__", "", x)}))
datataxa <- data.frame(lapply(datataxa, function(x) {gsub("D_3__", "", x)}))
datataxa <- data.frame(lapply(datataxa, function(x) {gsub("D_4__", "", x)}))
datataxa <- data.frame(lapply(datataxa, function(x) {gsub("D_5__", "", x)}))
datataxa <- data.frame(lapply(datataxa, function(x) {gsub("D_6__", "", x)}))

##Determine Order level unknowns
#Determine number of unknowns in HiHi indicators#
datataxa%>%filter(group==1)%>%nrow - datataxa%>%filter(group==1)%>% filter(!grepl(pattern = "uncultured", x = Order.x)) %>%filter(!grepl(pattern = "unidentified", x = Order.x))%>% filter(!grepl(pattern = "metagenome", x = Order.x))%>%filter(!grepl(pattern = "Unknown", x = Order.x))%>%filter(as.character(Class.x) != as.character(Order.x))%>%filter(!grepl(pattern ="__", x=Order.x))%>%nrow #11

#Determine number of unknowns in HiLo indicators#
datataxa%>%filter(group==2)%>%nrow - datataxa%>%filter(group==2)%>% filter(!grepl(pattern = "uncultured", x = Order.x)) %>%filter(!grepl(pattern = "unidentified", x = Order.x))%>% filter(!grepl(pattern = "metagenome", x = Order.x))%>%filter(!grepl(pattern = "Unknown", x = Order.x))%>%filter(as.character(Class.x) != as.character(Order.x))%>%filter(!grepl(pattern ="__", x=Order.x))%>%nrow #14

##Determine Family level unknowns
#Determine number of unknowns in HiHi indicators#
datataxa%>%filter(group==1)%>%nrow - datataxa%>%filter(group==1)%>% filter(!grepl(pattern = "uncultured", x = Family.x)) %>%filter(!grepl(pattern = "unidentified", x = Family.x))%>% filter(!grepl(pattern = "metagenome", x = Family.x))%>%filter(!grepl(pattern = "Unknown", x = Family.x))%>%filter(as.character(Order.x) != as.character(Family.x))%>%filter(!grepl(pattern ="__", x=Family.x))%>%nrow #24

#Determine number of unknowns in HiLo indicators#
datataxa%>%filter(group==2)%>%nrow - datataxa%>%filter(group==2)%>% filter(!grepl(pattern = "uncultured", x = Family.x)) %>%filter(!grepl(pattern = "unidentified", x = Family.x))%>% filter(!grepl(pattern = "metagenome", x = Family.x))%>%filter(!grepl(pattern = "Unknown", x = Family.x))%>%filter(as.character(Order.x) != as.character(Family.x))%>%filter(!grepl(pattern ="__", x=Family.x))%>%nrow #39

##Determine Genus level unknowns
#Determine number of unknowns in HiHi indicators#
datataxa%>%filter(group==1)%>%nrow - datataxa%>%filter(group==1)%>% filter(!grepl(pattern = "uncultured", x = Genus.x)) %>%filter(!grepl(pattern = "unidentified", x = Genus.x))%>% filter(!grepl(pattern = "metagenome", x = Genus.x))%>%filter(!grepl(pattern = "Unknown", x = Genus.x))%>%filter(as.character(Family.x) != as.character(Genus.x))%>%filter(!grepl(pattern ="__", x=Genus.x))%>%nrow #33

#Determine number of unknowns in HiLo indicators#
datataxa%>%filter(group==2)%>%nrow - datataxa%>%filter(group==2)%>% filter(!grepl(pattern = "uncultured", x = Genus.x)) %>%filter(!grepl(pattern = "unidentified", x = Genus.x))%>% filter(!grepl(pattern = "metagenome", x = Genus.x))%>%filter(!grepl(pattern = "Unknown", x = Genus.x))%>%filter(as.character(Family.x) != as.character(Genus.x))%>%filter(!grepl(pattern ="__", x=Genus.x))%>%nrow #69

#Use the following script to summarize the mean, standard deviation (sd), median, and interquartile range (IQR) for a particular value (UF_Combo_Genus_Estuary) across a particular classification (AbbSiteAbbSeasonAbbYear)
ddply(F_HCS_Sediment, .(AbbSiteAbbSeasonAbbYear), summarise, mean=mean(UF_Combo_Genus_Estuary), sd=sd(UF_Combo_Genus_Estuary), median=median(UF_Combo_Genus_Estuary), IQR=IQR(UF_Combo_Genus_Estuary))





















#Extract metadata table from the genus phyloseq object#
LWSS_Meta_Y1 <- as(sample_data(RFLWSS_Y1_GenusAgenus), "data.frame")


##Create new columns in data table related to the Estuary AIP##

#Determine the sequence sums from Unfiltered, Combo, Genus Level IRL indicators and upload as a new column#
LWSS_Meta_Y1$UF_Combo_Genus_IRL <- prune_taxa(as.character(unfiltered_estuary_combo_sig_res_taxo_seqs_genus$ESV.ID[unfiltered_estuary_combo_sig_res_taxo_seqs_genus$log2FoldChange<0]), RFLWSS_Y1_GenusAgenus) %>% sample_sums() 

#Determine the sequence sums from Unfiltered, Combo, Genus Level SLE indicators and upload as a new column#
LWSS_Meta_Y1$UF_Combo_Genus_SLE <- prune_taxa(as.character(unfiltered_estuary_combo_sig_res_taxo_seqs_genus$ESV.ID[unfiltered_estuary_combo_sig_res_taxo_seqs_genus$log2FoldChange>0]), RFLWSS_Y1_GenusAgenus) %>% sample_sums() 

#Calculate Estuary AIP and add in new column#
LWSS_Meta_Y1$UF_Combo_Genus_Estuary <- with(LWSS_Meta_Y1, 100*UF_Combo_Genus_SLE/(UF_Combo_Genus_SLE+UF_Combo_Genus_IRL))

##Create new columns in data table related to the Muck AIP##

#Determine the sequence sums from Unfiltered, Combo, Genus Level 3MC indicators and upload as a new column#
LWSS_Meta_Y1$UF_Combo_Genus_3MC <- prune_taxa(as.character(unfiltered_muck_combo_sig_res_taxo_seqs_genus$ESV.ID[unfiltered_muck_combo_sig_res_taxo_seqs_genus$log2FoldChange<0]), RFLWSS_Y1_GenusAgenus) %>% sample_sums() 

#Determine the sequence sums from Unfiltered, Combo, Genus Level 0MC indicators and upload as a new column#
LWSS_Meta_Y1$UF_Combo_Genus_0MC <- prune_taxa(as.character(unfiltered_muck_combo_sig_res_taxo_seqs_genus$ESV.ID[unfiltered_muck_combo_sig_res_taxo_seqs_genus$log2FoldChange>0]), RFLWSS_Y1_GenusAgenus) %>% sample_sums() 

#Calculate Muck AIP and add in new column#
LWSS_Meta_Y1$UF_Combo_Genus_Muck <- with(LWSS_Meta_Y1, 100*UF_Combo_Genus_3MC/(UF_Combo_Genus_3MC+UF_Combo_Genus_0MC))

##Create new columns in data table related to the TOM/Cu AIP##

#Determine the sequence sums from Unfiltered, Combo, Genus Level HiHi indicators and upload as a new column#
LWSS_Meta_Y1$UF_Combo_Genus_HiHi <- prune_taxa(as.character(unfiltered_tomcu_combo_sig_res_taxo_seqs_genus$ESV.ID[unfiltered_tomcu_combo_sig_res_taxo_seqs_genus$log2FoldChange<0]), RFLWSS_Y1_GenusAgenus) %>% sample_sums() 

#Determine the sequence sums from Unfiltered, Combo, Genus Level HiLo indicators and upload as a new column#
LWSS_Meta_Y1$UF_Combo_Genus_HiLo <- prune_taxa(as.character(unfiltered_tomcu_combo_sig_res_taxo_seqs_genus$ESV.ID[unfiltered_tomcu_combo_sig_res_taxo_seqs_genus$log2FoldChange>0]), RFLWSS_Y1_GenusAgenus) %>% sample_sums() 

#Calculate TOM/Cu AIP and add in new column#
LWSS_Meta_Y1$UF_Combo_Genus_TOMCu <- with(LWSS_Meta_Y1, 100*UF_Combo_Genus_HiHi/(UF_Combo_Genus_HiHi+UF_Combo_Genus_HiLo))

#Extract metadata table from the genus phyloseq object#
LWSS_Meta_Y2 <- as(sample_data(RFLWSS_Y2_GenusAgenus), "data.frame")


##Create new columns in data table related to the Estuary AIP##

#Determine the sequence sums from Unfiltered, Combo, Genus Level IRL indicators and upload as a new column#
LWSS_Meta_Y2$UF_Combo_Genus_IRL <- prune_taxa(as.character(unfiltered_estuary_combo_sig_res_taxo_seqs_genus$ESV.ID[unfiltered_estuary_combo_sig_res_taxo_seqs_genus$log2FoldChange<0]), RFLWSS_Y2_GenusAgenus) %>% sample_sums() 

#Determine the sequence sums from Unfiltered, Combo, Genus Level SLE indicators and upload as a new column#
LWSS_Meta_Y2$UF_Combo_Genus_SLE <- prune_taxa(as.character(unfiltered_estuary_combo_sig_res_taxo_seqs_genus$ESV.ID[unfiltered_estuary_combo_sig_res_taxo_seqs_genus$log2FoldChange>0]), RFLWSS_Y2_GenusAgenus) %>% sample_sums() 

#Calculate Estuary AIP and add in new column#
LWSS_Meta_Y2$UF_Combo_Genus_Estuary <- with(LWSS_Meta_Y2, 100*UF_Combo_Genus_SLE/(UF_Combo_Genus_SLE+UF_Combo_Genus_IRL))

##Create new columns in data table related to the Muck AIP##

#Determine the sequence sums from Unfiltered, Combo, Genus Level 3MC indicators and upload as a new column#
LWSS_Meta_Y2$UF_Combo_Genus_3MC <- prune_taxa(as.character(unfiltered_muck_combo_sig_res_taxo_seqs_genus$ESV.ID[unfiltered_muck_combo_sig_res_taxo_seqs_genus$log2FoldChange<0]), RFLWSS_Y2_GenusAgenus) %>% sample_sums() 

#Determine the sequence sums from Unfiltered, Combo, Genus Level 0MC indicators and upload as a new column#
LWSS_Meta_Y2$UF_Combo_Genus_0MC <- prune_taxa(as.character(unfiltered_muck_combo_sig_res_taxo_seqs_genus$ESV.ID[unfiltered_muck_combo_sig_res_taxo_seqs_genus$log2FoldChange>0]), RFLWSS_Y2_GenusAgenus) %>% sample_sums() 

#Calculate Muck AIP and add in new column#
LWSS_Meta_Y2$UF_Combo_Genus_Muck <- with(LWSS_Meta_Y2, 100*UF_Combo_Genus_3MC/(UF_Combo_Genus_3MC+UF_Combo_Genus_0MC))

##Create new columns in data table related to the TOM/Cu AIP##

#Determine the sequence sums from Unfiltered, Combo, Genus Level HiHi indicators and upload as a new column#
LWSS_Meta_Y2$UF_Combo_Genus_HiHi <- prune_taxa(as.character(unfiltered_tomcu_combo_sig_res_taxo_seqs_genus$ESV.ID[unfiltered_tomcu_combo_sig_res_taxo_seqs_genus$log2FoldChange<0]), RFLWSS_Y2_GenusAgenus) %>% sample_sums() 

#Determine the sequence sums from Unfiltered, Combo, Genus Level HiLo indicators and upload as a new column#
LWSS_Meta_Y2$UF_Combo_Genus_HiLo <- prune_taxa(as.character(unfiltered_tomcu_combo_sig_res_taxo_seqs_genus$ESV.ID[unfiltered_tomcu_combo_sig_res_taxo_seqs_genus$log2FoldChange>0]), RFLWSS_Y2_GenusAgenus) %>% sample_sums() 

#Calculate TOM/Cu AIP and add in new column#
LWSS_Meta_Y2$UF_Combo_Genus_TOMCu <- with(LWSS_Meta_Y2, 100*UF_Combo_Genus_HiHi/(UF_Combo_Genus_HiHi+UF_Combo_Genus_HiLo))

LWSS_Meta_Total <- rbind(LWSS_Meta_Y1, LWSS_Meta_Y2)


#Upload csv with all of the environmental variables#
LWSS_Env_Total <- read.csv(file.choose())

colnames(LWSS_Env_Total)

#Make the row names the Sample IDs
LWSS_Env_Total <- LWSS_Env_Total %>% remove_rownames %>% column_to_rownames(var="ï..Sample.Name")

#Combine the metadata and environmental data tables#
LWSS_Meta_Env_Total <- merge(LWSS_Env_Total, LWSS_Meta_Total, by="row.names")

#Make the rownames a separate column#
rownames(LWSS_Meta_Env_Total) <- LWSS_Meta_Env_Total[,1]

#Convert all occurences of NaN to NA#
LWSS_Meta_Env_Total[LWSS_Meta_Env_Total ==  "NaN"] <- NA

#Add new columns combiing sampling year and target characteristics

LWSS_Meta_Env_Total$Estuary.Sampling.Year <- paste(LWSS_Meta_Env_Total$Estuary,LWSS_Meta_Env_Total$Sampling.Year)
LWSS_Meta_Env_Total$Muck.Sampling.Year <- paste(LWSS_Meta_Env_Total$Muck,LWSS_Meta_Env_Total$Sampling.Year)
LWSS_Meta_Env_Total$LOI.Cu.Sampling.Year <- paste(LWSS_Meta_Env_Total$LOI.Cu,LWSS_Meta_Env_Total$Sampling.Year)

LWSS_Meta_Env_Total$Season.Sampling.Year <- paste(LWSS_Meta_Env_Total$Season,LWSS_Meta_Env_Total$Sampling.Year)



###VISUALIZE AND TEST METRICS###

##Statistically test AIP differences between subcategories##

#Estuary AIP Stats#
kruskal.test(UF_Combo_Genus_Estuary ~ Estuary.Sampling.Year, data = LWSS_Meta_Env_Total)

dunnTest(UF_Combo_Genus_Estuary ~ Estuary.Sampling.Year, data = LWSS_Meta_Env_Total, method="bh")

Dunn <- dunnTest(UF_Combo_Genus_Estuary ~ Estuary.Sampling.Year, data = LWSS_Meta_Env_Total, method="bh")

#Get compact letter display for this dataset
Dunn = Dunn$res
cldList(comparison = Dunn$Comparison,
        p.value    = Dunn$P.adj,
        threshold  = 0.05)

# Group Letter           MonoLetter
# 1 IRLYear1      a        a  
# 2 IRLYear2      b         b 
# 3 SLEYear1      c          c
# 4 SLEYear2      c          c

#Muck AIP Stats#
kruskal.test(UF_Combo_Genus_Muck ~ Muck.Sampling.Year, data = LWSS_Meta_Env_Total)
dunnTest(UF_Combo_Genus_Muck ~ Muck.Sampling.Year, data = LWSS_Meta_Env_Total, method="bh")
Dunn = dunnTest(UF_Combo_Genus_Muck ~ Muck.Sampling.Year, data = LWSS_Meta_Env_Total, method="bh")

#Get compact letter display for this dataset
Dunn = Dunn$res
cldList(comparison = Dunn$Comparison,
        p.value    = Dunn$P.adj,
        threshold  = 0.05)

# Group Letter MonoLetter
# 1    MuckYear1      a        a  
# 2    MuckYear2      a        a  
# 3 MuckishYear1     ab        ab 
# 4 MuckishYear2      a        a  
# 5   MuckyYear1      a        a  
# 6   MuckyYear2      a        a  
# 7     NotYear1      c          c
# 8     NotYear2     bc         bc

#TOM/Cu AIP Stats#
kruskal.test(UF_Combo_Genus_TOMCu ~ LOI.Cu.Sampling.Year, data = subset(LWSS_Meta_Env_Total, TOM > 10))
dunnTest(UF_Combo_Genus_TOMCu ~ LOI.Cu.Sampling.Year, data = subset(LWSS_Meta_Env_Total, TOM > 10), method="bh")
Dunn = dunnTest(UF_Combo_Genus_TOMCu ~ LOI.Cu.Sampling.Year, data = subset(LWSS_Meta_Env_Total, TOM > 10), method="bh")

Dunn = Dunn$res
cldList(comparison = Dunn$Comparison,
        p.value    = Dunn$P.adj,
        threshold  = 0.05)

# 1 High.HighYear1      a         a 
# 2 High.HighYear2      a         a 
# 3  High.LowYear1      b          b
# 4  High.LowYear2      b          b

##Create Box Plots to visualize differences between subcategories##

#Estuary AIP BP#
ggplot(data = LWSS_Meta_Env_Total, aes( x = Estuary.Sampling.Year, y = UF_Combo_Genus_Estuary, color=Sampling.Year)) +
  geom_boxplot() +
  ylab("Affected Indicator Percentage") +
  xlab("Estuary by Sampling Year") +
  ggtitle("Affected Indicator Percentage") +
  theme(axis.text.x  = element_text(angle=45, hjust=1))+
  theme(plot.title = element_text(hjust = 0.5)) + 
  annotate("text", x = 1:4, y = c(40, 87, 105, 105), label = c("a", "b", "c", "c")) +
  #theme(legend.position="bottom")
  theme(legend.position = "none")

ggplot(data = LWSS_Meta_Env_Total, aes( x = Estuary.Sampling.Year, y = UF_Combo_Genus_Estuary, color=Sampling.Year)) +
  geom_boxplot() +
  ylab("Affected Indicator Percentage") +
  xlab("Estuary by Sampling Year") +
  theme(axis.text.x  = element_text(angle=45, hjust=1))+
  ggtitle("Affected Indicator Percentage") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  annotate("text", x = 1:4, y = c(40, 105, 87, 105), label = c("a", "c", "b", "c")) +
  scale_x_discrete(limits=c("IRL Year 1", "SLE Year 1", "IRL Year 2", "SLE Year 2")) +
  #theme(legend.position="bottom")
  theme(legend.position = "none")

#Muck AIP BP#
ggplot(data = LWSS_Meta_Env_Total, aes( x = Muck.Sampling.Year, y = UF_Combo_Genus_Muck, color=Sampling.Year)) +
  geom_boxplot() +
  ylab("Affected Indicator Percentage") +
  xlab("Muck Characteristics by Sampling Year") +
  ggtitle("Affected Indicator Percentage") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x  = element_text(angle=45, hjust=1))+
  scale_x_discrete(limits=c("Muck Year 1", "Mucky Year 1", "Muckish Year 1", "Not Year 1", "Muck Year 2", "Mucky Year 2", "Muckish Year 2", "Not Year 2"), labels=c("3 Year 1", "2 Year 1", "1 Year 1", "0 Year 1", "3 Year 2", "2 Year 2", "1 Year 2", "0 Year 2")) +
  annotate("text", x = 1:8, y = 105, label = c("a", "a", "ab", "c", "a", "a", "a", "bc")) +
  theme(legend.position="none")

ggplot(data = LWSS_Meta_Env_Total, aes( x = Muck.Sampling.Year, y = UF_Combo_Genus_Muck, color=Sampling.Year)) +
  geom_boxplot() +
  ylab("Affected Indicator Percentage") +
  xlab("Muck Characteristics by Sampling Year") +
  ggtitle("Affected Indicator Percentage") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x  = element_text(angle=45, hjust=1))+
  scale_x_discrete(limits=c("Muck Year 1", "Muck Year 2", "Mucky Year 1", "Mucky Year 2", "Muckish Year 1", "Muckish Year 2", "Not Year 1", "Not Year 2"), labels=c("3 Year 1", "3 Year 2", "2 Year 1", "2 Year 2", "1 Year 1", "1 Year 2", "0 Year 1", "0 Year 2")) +
  annotate("text", x = 1:8, y = 105, label = c("a", "a", "a", "a", "ab", "a", "c", "bc")) +
  theme(legend.position="none")

##Create Metric Box Plots TOM/Cu Category to visualize stats##
#TOM/Cu AIP BP#
ggplot(data = subset(LWSS_Meta_Env_Total, TOM > 10), aes( x = LOI.Cu.Sampling.Year, y = UF_Combo_Genus_TOMCu, color=Sampling.Year)) +
  geom_boxplot() +
  theme(axis.text.x  = element_text(angle=45, hjust=1))+
  ylab("Affected Indicator Percentage") +
  xlab("TOM/Cu Category by Year") +
  #ggtitle("Affected Indicator Percentage") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_x_discrete(limits=c("High.High Year 1", "High.Low Year 1", "High.High Year 2", "High.Low Year 2"), labels=c("HiHi Year 1", "HiLo Year 1", "HiHi Year 2", "HiLo Year 2")) +
  annotate("text", x = 1:4, y = 105, label = c("a", "b", "a", "b")) +
  #theme(axis.title.x=element_blank()) +
  #  theme(legend.position="bottom")
  theme(legend.position = "none")

ggplot(data = subset(LWSS_Meta_Env_Total, TOM > 10), aes( x = LOI.Cu.Sampling.Year, y = UF_Combo_Genus_TOMCu, color=Sampling.Year)) +
  geom_boxplot() +
  theme(axis.text.x  = element_text(angle=45, hjust=1))+
  ylab("Affected Indicator Percentage") +
  xlab("TOM/Cu Category by Year") +
  #ggtitle("Affected Indicator Percentage") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_x_discrete(labels=c("HiHi Year 1", "HiHi Year 2", "HiLo Year 1", "HiLo Year 2")) +
  annotate("text", x = 1:4, y = 105, label = c("a", "a", "b", "b")) +
  #theme(axis.title.x=element_blank()) +
  #  theme(legend.position="bottom")
  theme(legend.position = "none")

##Create Box Plots to visualize differences between sites##

#Estuary AIP by Site#
ggplot(data = LWSS_Meta_Env_Total, aes( x = Site, y = UF_Combo_Genus_Estuary)) +
  geom_boxplot() +
  theme(axis.text.x  = element_text(angle=45, hjust=1))+
  ylab("Affected Indicator Percentage") +
  xlab("Site") +
  #ggtitle("Overall Sediment Samples Estuary AIP")+
  #theme(legend.position = "none") +
  scale_x_discrete(labels=c("Barber Bridge", "Fort Pierce", "Harbor Branch Channel", "Harbortown Marina", "Hobe Sound", "Jensen Beach", "Jupiter Narrows", "Linkport", "Manatee Pocket", "Melbourne Causeway", "Merritt Island Causeway", "Middle Estuary", "North Fork", "Round Island", "Sebastian Inlet", "South Fork", "South Fork 2", "Vero Beach", "Vero Beach Marina"))+
  theme(plot.title = element_text(hjust = 0.5))

#Estuary AIP by Site by Season#
ggplot(data = LWSS_Meta_Env_Total, aes( x = Site, y = UF_Combo_Genus_Estuary, color=Season)) +
  geom_boxplot() +
  theme(axis.text.x  = element_text(angle=45, hjust=1))+
  ylab("Affected Indicator Percentage") +
  xlab("Site") +
  scale_x_discrete(labels=c("Barber Bridge", "Fort Pierce", "Harbor Branch Channel", "Harbortown Marina", "Hobe Sound", "Jensen Beach", "Jupiter Narrows", "Linkport", "Manatee Pocket", "Melbourne Causeway", "Merritt Island Causeway", "Middle Estuary", "North Fork", "Round Island", "Sebastian Inlet", "South Fork", "South Fork 2", "Vero Beach", "Vero Beach Marina"))+
  scale_color_manual(values = c("darkgray", "black"))

#Muck AIP by Site#
ggplot(data = LWSS_Meta_Env_Total, aes( x = Site, y = UF_Combo_Genus_Muck)) +
  geom_boxplot()+
  ylab("Affected Indicator Percentage") +
  xlab("Site") +
  theme(axis.text.x  = element_text(angle=45, hjust=1)) +
  #annotate("text", x=c(1, 2, 5, 6, 7, 11, 13, 14, 15, 18, 19), y=105, label=c("0/6", "0/12","0/12", "0/12", "0/12",  "0/12", "0/12", "0/6", "0/12", "0/12", "0/6"), size=3) +
  theme(axis.text.x  = element_text(angle=45, hjust=1)) +  scale_x_discrete(labels=c("Barber Bridge", "Fort Pierce", "Harbor Branch Channel", "Harbortown Marina", "Hobe Sound", "Jensen Beach", "Jupiter Narrows", "Linkport", "Manatee Pocket", "Melbourne Causeway", "Merritt Island Causeway", "Middle Estuary", "North Fork", "Round Island", "Sebastian Inlet", "South Fork", "South Fork 2", "Vero Beach", "Vero Beach Marina"), )+
  annotate("text", x=c(3, 4, 8, 9, 10, 12, 16, 17), y=105, label=c("17/18", "6/6", "2/12", "11/12", "9/12", "12/12","12/12", "1/12"), size=3, fontface=2)+
  scale_color_manual (name="Site", values=c("midnightblue", "blue", "red", "purple", "yellow", "violet", "maroon", "darkgreen", "dodgerblue", "gold",  "moccasin",  "darkgray", "powderblue", "green", "lightgray", "black", "slateblue", "greenyellow", "lightseagreen"), labels=c("Barber Bridge", "Fort Pierce", "Harbor Branch Channel", "Harbortown Marina", "Hobe Sound", "Jensen Beach", "Jupiter Narrows", "Linkport", "Manatee Pocket", "Melbourne Causeway", "Merritt Island Causeway", "Middle Estuary", "North Fork", "Round Island", "Sebastian Inlet", "South Fork", "South Fork 2", "Vero Beach", "Vero Beach Marina"))

#TOM/Cu AIP by Sites only with samples with >10% TOM#
ggplot(data = subset(LWSS_Meta_Env_Total, TOM > 10), aes( x = Site, y = UF_Combo_Genus_TOMCu, color=LOI.Cu)) +
  geom_boxplot()+
  ylab("Affected Indicator Percentage") +
  xlab("Site") +
  theme(axis.text.x  = element_text(angle=45, hjust=1))+
  scale_color_manual(name="TOM/Cu", labels=c("HiHi", "HiLo"), values = c("darkgray", "black"))
