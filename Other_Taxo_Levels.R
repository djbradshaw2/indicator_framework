##SPECIES##

getwd() #"C:/Users/bradshawd/Documents/Journal_Articles/New_Indicators"
setwd("C:/Users/bradshawd/Documents/Journal_Articles/New_Indicators/Species")
getwd() #"C:/Users/bradshawd/Documents/Journal_Articles/New_Indicators/Species"

##Upload QIIME2 information into phyloseq for taxnomic level you want to test##
species_aa_table = read.csv(file.choose(), row.names=1)
species_taxa_table = read.csv(file.choose(), row.names = 1)
species_taxa_table = as.matrix(species_taxa_table)

#Can use the same metadata table as asv level#
DottedMETA = read.delim(file.choose(), row.names=1) 

dotted_meta = sample_data(DottedMETA)
Species_OTU = otu_table(species_aa_table, taxa_are_rows = TRUE)
Species_TAX = tax_table(species_taxa_table)

SpeciesA <- merge_phyloseq(Species_OTU, Species_TAX, dotted_meta)

colnames(tax_table(SpeciesA))

###Create at more filtered set of LWSS Data###
ntaxa(SpeciesA) #4120#
nsamples(SpeciesA)#483

LWSSpeciesA <- subset_samples(SpeciesA, Survey=="Lagoon Wide")
nsamples(LWSSpeciesA)#324

LWSWSpeciesA <- subset_samples(LWSSpeciesA, Medium=="Water")
nsamples(LWSWSpeciesA)#132, 11 sites in triplicate (33) x 4 sampling periods (sps)

LWSSSpeciesA1 <- subset_samples(SpeciesA, Survey=="Lagoon Wide.Post Hurricane")
nsamples(LWSSSpeciesA1)#12

LWSSSpeciesA2 <- subset_samples(LWSSpeciesA, Medium=="Sediment")
nsamples(LWSSSpeciesA2)#192

LWSSSpeciesA <- merge_phyloseq(LWSSSpeciesA1, LWSSSpeciesA2)
nsamples(LWSSSpeciesA)#204

##Create phyloseqs for determining and testing indicators for Sediment
#Split samples by sampling period (SP)
LWSS_W16_SpeciesA <- subset_samples(LWSSSpeciesA, AbbSeasonAbbYear=="W16")
nsamples(LWSS_W16_SpeciesA)#45, 15 sites in triplicate

LWSS_D17_SpeciesA <- subset_samples(LWSSSpeciesA, AbbSeasonAbbYear=="D17")
nsamples(LWSS_D17_SpeciesA)#45,15 sites in triplicate

LWSS_W17_SpeciesA <- subset_samples(LWSSSpeciesA, AbbSeasonAbbYear=="W17")
nsamples(LWSS_W17_SpeciesA)#57, 19 sites in triplicate

LWSS_D18_SpeciesA <- subset_samples(LWSSSpeciesA, AbbSeasonAbbYear=="D18")
nsamples(LWSS_D18_SpeciesA)#57, 19 sites in triplicate

#Combine sampling periods into first and second years
LWSS_Y1_SpeciesA <- merge_phyloseq(LWSS_W16_SpeciesA, LWSS_D17_SpeciesA)
nsamples(LWSS_Y1_SpeciesA)#90

LWSS_Y2_SpeciesA <- merge_phyloseq(LWSS_W17_SpeciesA, LWSS_D18_SpeciesA)
nsamples(LWSS_Y2_SpeciesA)#114

##Create phyloseqs for determining and testing indicators for Water
#Split samples by sampling period (SP)
LWSW_W16_SpeciesA <- subset_samples(LWSWSpeciesA, AbbSeasonAbbYear=="W16")
nsamples(LWSW_W16_SpeciesA)#33, 11 sites in triplicate

LWSW_D17_SpeciesA <- subset_samples(LWSWSpeciesA, AbbSeasonAbbYear=="D17")
nsamples(LWSW_D17_SpeciesA)#33,11 sites in triplicate

LWSW_W17_SpeciesA <- subset_samples(LWSWSpeciesA, AbbSeasonAbbYear=="W17")
nsamples(LWSW_W17_SpeciesA)#33, 11 sites in triplicate

LWSW_D18_SpeciesA <- subset_samples(LWSWSpeciesA, AbbSeasonAbbYear=="D18")
nsamples(LWSW_D18_SpeciesA)#33, 11 sites in triplicate

#Combine sampling periods into first and second years
LWSW_Y1_SpeciesA <- merge_phyloseq(LWSW_W16_SpeciesA, LWSW_D17_SpeciesA)
nsamples(LWSW_Y1_SpeciesA)#66

LWSW_Y2_SpeciesA <- merge_phyloseq(LWSW_W17_SpeciesA, LWSW_D18_SpeciesA)
nsamples(LWSW_Y2_SpeciesA)#66

##Filter yearly sample sets of all zero and then low abundance speciess
FLWSS_Y1_SpeciesA <- filter_taxa(LWSS_Y1_SpeciesA, function(x) sum(x) >0, TRUE)
ntaxa(FLWSS_Y1_SpeciesA) #2760#
RFLWSS_Y1_SpeciesAspecies <- filter_taxa(FLWSS_Y1_SpeciesA, function(x) sum(x) >20, TRUE)
ntaxa(RFLWSS_Y1_SpeciesAspecies) #1531#

FLWSS_Y2_SpeciesA <- filter_taxa(LWSS_Y2_SpeciesA, function(x) sum(x) >0, TRUE)
ntaxa(FLWSS_Y2_SpeciesA) #3039#
RFLWSS_Y2_SpeciesAspecies <- filter_taxa(FLWSS_Y2_SpeciesA, function(x) sum(x) >20, TRUE)
ntaxa(RFLWSS_Y2_SpeciesAspecies) #1367#

FLWSW_Y1_SpeciesA <- filter_taxa(LWSW_Y1_SpeciesA, function(x) sum(x) >0, TRUE)
ntaxa(FLWSW_Y1_SpeciesA) #1511#
RFLWSW_Y1_SpeciesAspecies <- filter_taxa(FLWSW_Y1_SpeciesA, function(x) sum(x) >20, TRUE)
ntaxa(RFLWSW_Y1_SpeciesAspecies) #682#

FLWSW_Y2_SpeciesA <- filter_taxa(LWSW_Y2_SpeciesA, function(x) sum(x) >0, TRUE)
ntaxa(FLWSW_Y2_SpeciesA) #2087#
RFLWSW_Y2_SpeciesAspecies <- filter_taxa(FLWSW_Y2_SpeciesA, function(x) sum(x) >20, TRUE)
ntaxa(RFLWSW_Y2_SpeciesAspecies) #828#





#Determine the number of samples in each category and subcategory#
length(which(sample_data(RFLWSS_Y1_SpeciesAspecies)$Estuary == "SLE")) #24
length(which(sample_data(RFLWSS_Y1_SpeciesAspecies)$Estuary == "IRL")) #66

length(which(sample_data(RFLWSS_Y1_SpeciesAspecies)$Muck == "Not")) #60
length(which(sample_data(RFLWSS_Y1_SpeciesAspecies)$Muck == "Muck")) #22
length(which(sample_data(RFLWSS_Y1_SpeciesAspecies)$Muck == "Mucky")) #4
length(which(sample_data(RFLWSS_Y1_SpeciesAspecies)$Muck == "Muckish")) #4

length(which(sample_data(RFLWSS_Y1_SpeciesAspecies)$LOI.Cu == "Low.Low")) #61
length(which(sample_data(RFLWSS_Y1_SpeciesAspecies)$LOI.Cu == "High.Low")) #18
length(which(sample_data(RFLWSS_Y1_SpeciesAspecies)$LOI.Cu == "High.High")) #11
length(which(sample_data(RFLWSS_Y1_SpeciesAspecies)$LOI.Cu == "Low.High")) #0

##NOTE##
#Use res to determine what is + or .#
#muck res: Not (+) vs Muck (.)#
#tomcu res: HL (+) vs HH (.)#
#estuary res: SLE (+) vs IRL (.)#

###RUN DESEQ2 AT species LEVEL ON ESTUARY SUBCATEGORIES###
#Convert phyloseq to deseq2 object centered around a irlvssle factor#
estuary_cudds_species = phyloseq_to_deseq2(RFLWSS_Y1_SpeciesAspecies,  ~ Estuary)

#Calculate geometric means prior to estimate size factors#
estuary_geoMeans_species = apply(counts(estuary_cudds_species), 1, gm_mean)
estuary_cudds_species = estimateSizeFactors(estuary_cudds_species, geoMeans = estuary_geoMeans_species)

#Conduct DESEQ2 test#
estuary_cudds_species = DESeq(estuary_cudds_species, fitType="local")

#Explore the results#
estuary_DESeq2_res_species = results(estuary_cudds_species)
estuary_DESeq2_res_species = estuary_DESeq2_res_species[order(estuary_DESeq2_res_species$padj, na.last=NA), ]
alpha = 0.05
estuary_DESeq2_sig_res_species = estuary_DESeq2_res_species[(estuary_DESeq2_res_species$padj < alpha), ]

#Make dataframe with taxanomy added in#
estuary_DESeq2_sig_res_taxo_species = cbind(as(estuary_DESeq2_sig_res_species, "data.frame"), as(tax_table(RFLWSS_Y1_SpeciesAspecies)[rownames(estuary_DESeq2_sig_res_species), ], "matrix"))

#Make dataframe with speciess added in from all sites#
estuary_DESeq2_sig_res_taxo_seqs_species = cbind(as(estuary_DESeq2_sig_res_taxo_species, "data.frame"), as(otu_table(RFLWSS_Y1_SpeciesAspecies)[rownames(estuary_DESeq2_sig_res_taxo_species), ], "matrix"))

#Make rownames an actual column and remove old rownames#
estuary_DESeq2_sig_res_taxo_seqs_species <- cbind(ESV.ID = rownames(estuary_DESeq2_sig_res_taxo_seqs_species), estuary_DESeq2_sig_res_taxo_seqs_species)
rownames(estuary_DESeq2_sig_res_taxo_seqs_species) <-NULL

#Make a txt file of results for all sites#
write.csv(as.data.frame(estuary_DESeq2_sig_res_taxo_seqs_species), file="estuary_DESeq2_sig_res_taxo_seqs_species.csv")

#Detemine which subcategory is negative or positive#
estuary_DESeq2_res_species
#estuary_DESeq2_res_species: SLE (+) vs IRL (.)#

#Determine number of Indicators with padj >0.05 in each subcagetory#
#Number of SLE indicators# 
length(which(estuary_DESeq2_sig_res_taxo_seqs_species$log2FoldChange > 0)) #259
#Number of IRL indicators#
length(which(estuary_DESeq2_sig_res_taxo_seqs_species$log2FoldChange < 0)) #346

###RUN DESEQ2 AT species LEVEL ON MUCK EXTREME SUBCATEGORIES###

##Prepare for running DESeq2##
#Make phyloseq focusing on Muck and Not subcategories#  
RFLWSS_Y1_Muck_SpeciesAspecies = subset_samples(RFLWSS_Y1_SpeciesAspecies, Muck=="Muck")
RFLWSS_Y1_Not_SpeciesAspecies = subset_samples(RFLWSS_Y1_SpeciesAspecies, Muck=="Not")
RFLWSS_Y1_MuckvsNot_SpeciesAspecies = merge_phyloseq(RFLWSS_Y1_Muck_SpeciesAspecies, RFLWSS_Y1_Not_SpeciesAspecies)

##Run DESEQ2 and make a txt table of results##
#Convert phyloseq to deseq2 object centered around a muck factor#
muck_cudds_species = phyloseq_to_deseq2(RFLWSS_Y1_MuckvsNot_SpeciesAspecies,  ~ Muck)

#Calculate geometric means prior to estimate size factors#
muck_geoMeans_species = apply(counts(muck_cudds_species), 1, gm_mean)
muck_cudds_species = estimateSizeFactors(muck_cudds_species, geoMeans = muck_geoMeans_species)

#Conduct DESEQ2 test#
muck_cudds_species = DESeq(muck_cudds_species, fitType="local")

#Explore the results#
muck_DESeq2_res_species = results(muck_cudds_species)
muck_DESeq2_res_species = muck_DESeq2_res_species[order(muck_DESeq2_res_species$padj, na.last=NA), ]
alpha = 0.05
muck_DESeq2_sig_res_species = muck_DESeq2_res_species[(muck_DESeq2_res_species$padj < alpha), ]

#Make dataframe with taxanomy added in#
muck_DESeq2_sig_res_taxo_species = cbind(as(muck_DESeq2_sig_res_species, "data.frame"), as(tax_table(RFLWSS_Y1_MuckvsNot_SpeciesAspecies)[rownames(muck_DESeq2_sig_res_species), ], "matrix"))

#Make dataframe with speciess added in from all sites#
muck_DESeq2_sig_res_taxo_seqs_species = cbind(as(muck_DESeq2_sig_res_taxo_species, "data.frame"), as(otu_table(RFLWSS_Y1_SpeciesAspecies)[rownames(muck_DESeq2_sig_res_taxo_species), ], "matrix"))

#Make rownames an actual column and remove old rownames#
muck_DESeq2_sig_res_taxo_seqs_species <- cbind(ESV.ID = rownames(muck_DESeq2_sig_res_taxo_seqs_species), muck_DESeq2_sig_res_taxo_seqs_species)
rownames(muck_DESeq2_sig_res_taxo_seqs_species) <-NULL

#Make a txt file of results for all sites#
write.csv(as.data.frame(muck_DESeq2_sig_res_taxo_seqs_species), file="muck_DESeq2_sig_res_taxo_seqs_species.csv")

#Detemine which subcategory is negative or positive#
muck_DESeq2_res_species
#muck_DESeq2_res_species: Not (+) vs Muck (.)#
#Determine number of Indicators with padj >0.05 in each subcagetory#
#Number of Not indicators# 
length(which(muck_DESeq2_sig_res_taxo_seqs_species$log2FoldChange > 0)) #204
#Number of Muck indicators#
length(which(muck_DESeq2_sig_res_taxo_seqs_species$log2FoldChange < 0)) #114

###RUN DESEQ2 AT species LEVEL ON TOM/Cu EXTREME SUBCATEGORIES###

##Prepare for running DESeq2##
#Make phyloseq focusing on High.High and High.Low subcategories#  
RFLWSS_Y1_HH_SpeciesAspecies = subset_samples(RFLWSS_Y1_SpeciesAspecies, LOI.Cu=="High.High")
RFLWSS_Y1_HL_SpeciesAspecies = subset_samples(RFLWSS_Y1_SpeciesAspecies, LOI.Cu=="High.Low")
RFLWSS_Y1_HHvsHL_SpeciesAspecies = merge_phyloseq(RFLWSS_Y1_HH_SpeciesAspecies, RFLWSS_Y1_HL_SpeciesAspecies)

##Run DESEQ2 and make a txt table of results##
#Convert phyloseq to deseq2 object centered around a loicu factor#
tomcu_cudds_species = phyloseq_to_deseq2(RFLWSS_Y1_HHvsHL_SpeciesAspecies,  ~ LOI.Cu)

#Calculate geometric means prior to estimate size factors#
tomcu_geoMeans_species = apply(counts(tomcu_cudds_species), 1, gm_mean)
tomcu_cudds_species = estimateSizeFactors(tomcu_cudds_species, geoMeans = tomcu_geoMeans_species)

#Conduct DESEQ2 test#
tomcu_cudds_species = DESeq(tomcu_cudds_species, fitType="local")

#Explore the results#
tomcu_DESeq2_res_species = results(tomcu_cudds_species)
tomcu_DESeq2_res_species = tomcu_DESeq2_res_species[order(tomcu_DESeq2_res_species$padj, na.last=NA), ]
alpha = 0.05
tomcu_DESeq2_sig_res_species = tomcu_DESeq2_res_species[(tomcu_DESeq2_res_species$padj < alpha), ]

#Make dataframe with taxanomy added in#
tomcu_DESeq2_sig_res_taxo_species = cbind(as(tomcu_DESeq2_sig_res_species, "data.frame"), as(tax_table(RFLWSS_Y1_HHvsHL_SpeciesAspecies)[rownames(tomcu_DESeq2_sig_res_species), ], "matrix"))

#Make dataframe with speciess added in from all sites#
tomcu_DESeq2_sig_res_taxo_seqs_species = cbind(as(tomcu_DESeq2_sig_res_taxo_species, "data.frame"), as(otu_table(RFLWSS_Y1_SpeciesAspecies)[rownames(tomcu_DESeq2_sig_res_taxo_species), ], "matrix"))

#Make rownames an actual column and remove old rownames#
tomcu_DESeq2_sig_res_taxo_seqs_species <- cbind(ESV.ID = rownames(tomcu_DESeq2_sig_res_taxo_seqs_species), tomcu_DESeq2_sig_res_taxo_seqs_species)
rownames(tomcu_DESeq2_sig_res_taxo_seqs_species) <-NULL

#Make a txt file of results for all sites#
write.csv(as.data.frame(tomcu_DESeq2_sig_res_taxo_seqs_species), file="tomcu_DESeq2_sig_res_taxo_seqs_species.csv")

#Detemine which subcategory is negative or positive#
tomcu_DESeq2_res_species
#tomcu_DESeq2_res_species: High.Low (+) vs High.High (.)#
#Determine number of Indicators with padj >0.05 in each subcagetory#
#Number of High.Low indicators# 
length(which(tomcu_DESeq2_sig_res_taxo_seqs_species$log2FoldChange > 0)) #168
#Number of High.High indicators#
length(which(tomcu_DESeq2_sig_res_taxo_seqs_species$log2FoldChange < 0)) #101

#Order rows by log2FoldChange#
muck_DESeq2_sig_res_taxo_seqs_species <- arrange(muck_DESeq2_sig_res_taxo_seqs_species, log2FoldChange)
tomcu_DESeq2_sig_res_taxo_seqs_species <- arrange(tomcu_DESeq2_sig_res_taxo_seqs_species, log2FoldChange)
estuary_DESeq2_sig_res_taxo_seqs_species <- arrange(estuary_DESeq2_sig_res_taxo_seqs_species, log2FoldChange)

#Create lists of just significant speciess for each category#
unfiltered_DESeq2_estuary_speciess <- subset(estuary_DESeq2_sig_res_taxo_seqs_species, select=c(ESV.ID))
unfiltered_DESeq2_muck_speciess <- subset(muck_DESeq2_sig_res_taxo_seqs_species, select=c(ESV.ID))
unfiltered_DESeq2_tomcu_speciess <- subset(tomcu_DESeq2_sig_res_taxo_seqs_species, select=c(ESV.ID))



##Create indicator tables filtered of overlapping indicators with other two metadata categories##

#Filter out rows in the IRLSLE DESEQ object that are also in the LOI.Cu and Muck DESEQs#
filtered_estuary_DESeq2_sig_res_taxo_seqs_species <- anti_join(estuary_DESeq2_sig_res_taxo_seqs_species, unfiltered_DESeq2_tomcu_speciess, by = "ESV.ID")
filtered_estuary_DESeq2_sig_res_taxo_seqs_species <- anti_join(filtered_estuary_DESeq2_sig_res_taxo_seqs_species, unfiltered_DESeq2_muck_speciess, by = "ESV.ID")

#Filter out rows in the Muck DESEQ object that are also in the IRL.SLE and LOI.Cu DESEQs#
filtered_muck_DESeq2_sig_res_taxo_seqs_species <- anti_join(muck_DESeq2_sig_res_taxo_seqs_species, unfiltered_DESeq2_estuary_speciess, by = "ESV.ID")
filtered_muck_DESeq2_sig_res_taxo_seqs_species <- anti_join(filtered_muck_DESeq2_sig_res_taxo_seqs_species, unfiltered_DESeq2_tomcu_speciess, by = "ESV.ID")

#Filter out rows in the LOI.CU DESEQ object that are also in the IRL.SLE and Muck DESEQs#
filtered_tomcu_DESeq2_sig_res_taxo_seqs_species <- anti_join(tomcu_DESeq2_sig_res_taxo_seqs_species, unfiltered_DESeq2_estuary_speciess, by = "ESV.ID")
filtered_tomcu_DESeq2_sig_res_taxo_seqs_species <- anti_join(filtered_tomcu_DESeq2_sig_res_taxo_seqs_species, unfiltered_DESeq2_muck_speciess, by = "ESV.ID")



#Order rows by log2FoldChange#
filtered_muck_DESeq2_sig_res_taxo_seqs_species <- arrange(filtered_muck_DESeq2_sig_res_taxo_seqs_species, log2FoldChange)
filtered_tomcu_DESeq2_sig_res_taxo_seqs_species <- arrange(filtered_tomcu_DESeq2_sig_res_taxo_seqs_species, log2FoldChange)
filtered_estuary_DESeq2_sig_res_taxo_seqs_species <- arrange(filtered_estuary_DESeq2_sig_res_taxo_seqs_species, log2FoldChange)

#Make a csv file for each of the filtered tables#
write.csv(as.data.frame(filtered_estuary_DESeq2_sig_res_taxo_seqs_species), file="filtered_estuary_DESeq2_sig_res_taxo_seqs_species.csv")
write.csv(as.data.frame(filtered_muck_DESeq2_sig_res_taxo_seqs_species), file="filtered_muck_DESeq2_sig_res_taxo_seqs_species.csv")
write.csv(as.data.frame(filtered_tomcu_DESeq2_sig_res_taxo_seqs_species), file="filtered_tomcu_DESeq2_sig_res_taxo_seqs_species.csv")



##Determine number of filtered indicators##

#Number of Not indicators# 
length(which(filtered_muck_DESeq2_sig_res_taxo_seqs_species$log2FoldChange > 0)) #44
#Number of Muck indicators#
length(which(filtered_muck_DESeq2_sig_res_taxo_seqs_species$log2FoldChange < 0)) #49
#Number of High.Low indicators# 
length(which(filtered_tomcu_DESeq2_sig_res_taxo_seqs_species$log2FoldChange > 0)) #18
#Number of High.High indicators#
length(which(filtered_tomcu_DESeq2_sig_res_taxo_seqs_species$log2FoldChange < 0)) #17
#Number of SLE indicators# 
length(which(filtered_estuary_DESeq2_sig_res_taxo_seqs_species$log2FoldChange > 0)) #103
#Number of IRL indicators#
length(which(filtered_estuary_DESeq2_sig_res_taxo_seqs_species$log2FoldChange < 0)) #137

#Create lists of filtered significant speciess for each category#
filtered_DESeq2_estuary_speciess <- subset(filtered_estuary_DESeq2_sig_res_taxo_seqs_species, select=c(ESV.ID))
filtered_DESeq2_muck_speciess <- subset(filtered_muck_DESeq2_sig_res_taxo_seqs_species, select=c(ESV.ID))
filtered_DESeq2_tomcu_speciess <- subset(filtered_tomcu_DESeq2_sig_res_taxo_seqs_species, select=c(ESV.ID))

##Create phyloseq objects for DESeq2 indicators##

#Create phyloseq objects of the unfiltered indicators using character string version of significant speciess##
charac_unfiltered_DESeq2_estuary_speciess <- unfiltered_DESeq2_estuary_speciess[,"ESV.ID"]
charac_unfiltered_DESeq2_estuary_speciess <- as.character(charac_unfiltered_DESeq2_estuary_speciess)
UFDestuarySpeciesAspecies <- prune_taxa(charac_unfiltered_DESeq2_estuary_speciess, RFLWSS_Y1_SpeciesAspecies)

charac_unfiltered_DESeq2_muck_speciess <- unfiltered_DESeq2_muck_speciess[,"ESV.ID"]
charac_unfiltered_DESeq2_muck_speciess <- as.character(charac_unfiltered_DESeq2_muck_speciess)
UFDmuckSpeciesAspecies <- prune_taxa(charac_unfiltered_DESeq2_muck_speciess, RFLWSS_Y1_SpeciesAspecies)

charac_unfiltered_DESeq2_tomcu_speciess <- unfiltered_DESeq2_tomcu_speciess[,"ESV.ID"]
charac_unfiltered_DESeq2_tomcu_speciess <- as.character(charac_unfiltered_DESeq2_tomcu_speciess)
UFDtomcuSpeciesAspecies <- prune_taxa(charac_unfiltered_DESeq2_tomcu_speciess, RFLWSS_Y1_SpeciesAspecies)

#Create phyloseq objects of the filtered indicators using character string version of significant speciess##
charac_filtered_DESeq2_estuary_speciess <- filtered_DESeq2_estuary_speciess[,"ESV.ID"]
charac_filtered_DESeq2_estuary_speciess <- as.character(charac_filtered_DESeq2_estuary_speciess)
FDestuarySpeciesAspecies <- prune_taxa(charac_filtered_DESeq2_estuary_speciess, RFLWSS_Y1_SpeciesAspecies)

charac_filtered_DESeq2_muck_speciess <- filtered_DESeq2_muck_speciess[,"ESV.ID"]
charac_filtered_DESeq2_muck_speciess <- as.character(charac_filtered_DESeq2_muck_speciess)
FDmuckSpeciesAspecies <- prune_taxa(charac_filtered_DESeq2_muck_speciess, RFLWSS_Y1_SpeciesAspecies)

charac_filtered_DESeq2_tomcu_speciess <- filtered_DESeq2_tomcu_speciess[,"ESV.ID"]
charac_filtered_DESeq2_tomcu_speciess <- as.character(charac_filtered_DESeq2_tomcu_speciess)
FDtomcuSpeciesAspecies <- prune_taxa(charac_filtered_DESeq2_tomcu_speciess, RFLWSS_Y1_SpeciesAspecies)



##Create heatmaps to see overall trends, may take a while to load##

unfiltered_DESeq2_estuary_species_heatmap <- plot_heatmap(UFDestuarySpeciesAspecies, taxa.order = charac_unfiltered_DESeq2_estuary_speciess, sample.order= "Estuary", sample.label = "Estuary")
unfiltered_DESeq2_estuary_species_heatmap

unfiltered_DESeq2_muck_species_heatmap <- plot_heatmap(UFDmuckSpeciesAspecies, taxa.order = charac_unfiltered_DESeq2_muck_speciess, sample.order= "Muck", sample.label = "Muck")
unfiltered_DESeq2_muck_species_heatmap

unfiltered_DESeq2_tomcu_species_heatmap <- plot_heatmap(UFDtomcuSpeciesAspecies, taxa.order = charac_unfiltered_DESeq2_tomcu_speciess, sample.order= "LOI.Cu", sample.label = "LOI.Cu")
unfiltered_DESeq2_tomcu_species_heatmap

filtered_DESeq2_estuary_species_heatmap <- plot_heatmap(FDestuarySpeciesAspecies, taxa.order = charac_filtered_DESeq2_estuary_speciess, sample.order= "Estuary", sample.label = "Estuary")
filtered_DESeq2_estuary_species_heatmap

filtered_DESeq2_muck_species_heatmap <- plot_heatmap(FDmuckSpeciesAspecies, taxa.order = charac_filtered_DESeq2_muck_speciess, sample.order= "Muck", sample.label = "Muck")
filtered_DESeq2_muck_species_heatmap

filtered_DESeq2_tomcu_species_heatmap <- plot_heatmap(FDtomcuSpeciesAspecies, taxa.order = charac_filtered_DESeq2_tomcu_speciess, sample.order= "LOI.Cu", sample.label = "LOI.Cu")
filtered_DESeq2_tomcu_species_heatmap


###CONDUCT INDICSPECIES ANALYSIS AT species LEVEL###

#Increase max print so that all of results can be printed#
options(max.print=100000000)


##Conduct species indicspecies analysis for Estuary Category#

#Take out OTU aka species table from metadata focused phyloseq object#
estuary_seqs_species = as(otu_table(RFLWSS_Y1_SpeciesAspecies), "matrix")

#Make it a data frame, transpose it, make it a data frame again#
estuary_seqs_species <- as.data.frame(estuary_seqs_species)
estuary_seqs_species<- t(estuary_seqs_species)
estuary_seqs_species <- as.data.frame(estuary_seqs_species)

#Run indicspecies with same estuary_meta_group as in asv section#
estuary_indval_species = multipatt(estuary_seqs_species, estuary_meta_group, control = how(nperm=9999))

#Create a blank csv, then sink all indicspecies results into it#
sink("estuary_IS_res_species.csv")
estuary_sig_indval_species <- summary(estuary_indval_species, indvalcomp=TRUE, alpha=1)
sink()

##Conduct species indicspecies analysis for Muck Category#

#Take out OTU aka species table from metadata focused phyloseq object#
muck_seqs_species = as(otu_table(RFLWSS_Y1_MuckvsNot_SpeciesAspecies), "matrix")

#Make it a data frame, transpose it, make it a data frame again#
muck_seqs_species <- as.data.frame(muck_seqs_species)
muck_seqs_species<- t(muck_seqs_species)
muck_seqs_species <- as.data.frame(muck_seqs_species)

#Run indicspecies with same muck_meta_group as asv section#
muck_indval_species = multipatt(muck_seqs_species, muck_meta_group, control = how(nperm=9999))

#Create a blank csv, then sink all indicspecies results into it#
sink("muck_IS_res_species.csv")
muck_sig_indval_species <- summary(muck_indval_species, indvalcomp=TRUE, alpha=1)
sink()


##Conduct species indicspecies analysis for TOM/Cu Category#

#Take out OTU aka species table from metadata focused phyloseq object#
tomcu_seqs_species = as(otu_table(RFLWSS_Y1_HHvsHL_SpeciesAspecies), "matrix")

#Make it a data frame, transpose it, make it a data frame again#
tomcu_seqs_species <- as.data.frame(tomcu_seqs_species)
tomcu_seqs_species<- t(tomcu_seqs_species)
tomcu_seqs_species <- as.data.frame(tomcu_seqs_species)

#Run indicspecies with same tomcu_meta_group as asv section#
tomcu_indval_species = multipatt(tomcu_seqs_species, tomcu_meta_group, control = how(nperm=9999))

#Create a blank csv, then sink all indicspecies results into it#
sink("tomcu_IS_res_species.csv")
tomcu_sig_indval_species <- summary(tomcu_indval_species, indvalcomp=TRUE, alpha=1)
sink()

#Reset max print to default#
options(max.print=99999)

##Look at csv file you sinked into and create a table
##of just species.ID, A, B, stat, p.value, and group##

##Continue species indicspecies analysis for Estuary Category##

#Upload new table#
estuary_IS_res_species <- read.csv(file.choose())
#Make a new p value adjusted column#
estuary_IS_res_species$padj <- p.adjust(estuary_IS_res_species$p.value, method = "BH")

#Keep table of adjusted p values less than 0.05#
estuary_IS_sig_res_species <- filter(estuary_IS_res_species, padj < 0.05)

#Make dataframe with taxanomy added in by making ESV.ID the rownames then adding taxonomy table#
estuary_IS_sig_res_taxo_species <- as.data.frame(estuary_IS_sig_res_species)
rownames(estuary_IS_sig_res_taxo_species) <- estuary_IS_sig_res_taxo_species[,1]
estuary_IS_sig_res_taxo_species = cbind(as(estuary_IS_sig_res_taxo_species, "data.frame"), as(tax_table(RFLWSS_Y1_SpeciesAspecies)[rownames(estuary_IS_sig_res_taxo_species), ], "matrix"))

#Make dataframe with speciess added in from all sites#
estuary_IS_sig_res_taxo_seqs_species = cbind(as(estuary_IS_sig_res_taxo_species, "data.frame"), as(otu_table(RFLWSS_Y1_SpeciesAspecies)[rownames(estuary_IS_sig_res_taxo_species), ], "matrix"))

#Remove old rownames#
rownames(estuary_IS_sig_res_taxo_seqs_species) <-NULL

#Make a txt file of results for all sites#
write.csv(as.data.frame(estuary_IS_sig_res_taxo_seqs_species), file="estuary_IS_sig_res_taxo_seqs_species.csv")

#Check the number of indicators in each group after padj#
length(which(sample_data(estuary_IS_sig_res_species)$group == "1")) #202 IRL
length(which(sample_data(estuary_IS_sig_res_species)$group == "2")) #431 SLE


##Continue species indicspecies analysis for Muck Category##

#Upload new table#
muck_IS_res_species <- read.csv(file.choose())

#Make a new p value adjusted column#
muck_IS_res_species$padj <- p.adjust(muck_IS_res_species$p.value, method = "BH")

#Keep table of adjusted p values less than 0.05#
muck_IS_sig_res_species <- filter(muck_IS_res_species, padj < 0.05)

#Make dataframe with taxanomy added in by making ESV.ID the rownames then adding taxonomy table#
muck_IS_sig_res_taxo_species <- as.data.frame(muck_IS_sig_res_species)
rownames(muck_IS_sig_res_taxo_species) <- muck_IS_sig_res_taxo_species[,1]
muck_IS_sig_res_taxo_species = cbind(as(muck_IS_sig_res_taxo_species, "data.frame"), as(tax_table(RFLWSS_Y1_MuckvsNot_SpeciesAspecies)[rownames(muck_IS_sig_res_taxo_species), ], "matrix"))

#Make dataframe with speciess added in from all sites#
muck_IS_sig_res_taxo_seqs_species = cbind(as(muck_IS_sig_res_taxo_species, "data.frame"), as(otu_table(RFLWSS_Y1_SpeciesAspecies)[rownames(muck_IS_sig_res_taxo_species), ], "matrix"))

#Remove old rownames#
rownames(muck_IS_sig_res_taxo_seqs_species) <-NULL

#Make a txt file of results for all sites#
write.csv(as.data.frame(muck_IS_sig_res_taxo_seqs_species), file="muck_IS_sig_res_taxo_seqs_species.csv")

#Check the number of indicators in each group after padj#
length(which(sample_data(muck_IS_sig_res_species)$group == "1")) #171 Muck
length(which(sample_data(muck_IS_sig_res_species)$group == "2")) #575 Not

##Continue species indicspecies analysis for TOM/Cu Category##

#Upload new table#
tomcu_IS_res_species <- read.csv(file.choose())
#Make a new p value adjusted column#
tomcu_IS_res_species$padj <- p.adjust(tomcu_IS_res_species$p.value, method = "BH")

#Keep table of adjusted p values less than 0.05#
tomcu_IS_sig_res_species <- filter(tomcu_IS_res_species, padj < 0.05)

#Make dataframe with taxanomy added in by making ESV.ID the rownames then adding taxonomy table#
tomcu_IS_sig_res_taxo_species <- as.data.frame(tomcu_IS_sig_res_species)
rownames(tomcu_IS_sig_res_taxo_species) <- tomcu_IS_sig_res_taxo_species[,1]
tomcu_IS_sig_res_taxo_species = cbind(as(tomcu_IS_sig_res_taxo_species, "data.frame"), as(tax_table(RFLWSS_Y1_MuckvsNot_SpeciesAspecies)[rownames(tomcu_IS_sig_res_taxo_species), ], "matrix"))

#Make dataframe with speciess added in from all sites#
tomcu_IS_sig_res_taxo_seqs_species = cbind(as(tomcu_IS_sig_res_taxo_species, "data.frame"), as(otu_table(RFLWSS_Y1_SpeciesAspecies)[rownames(tomcu_IS_sig_res_taxo_species), ], "matrix"))

#Remove old rownames#
rownames(tomcu_IS_sig_res_taxo_seqs_species) <-NULL

#Make a txt file of results for all sites#
write.csv(as.data.frame(tomcu_IS_sig_res_taxo_seqs_species), file="tomcu_IS_sig_res_taxo_seqs_species.csv")

#Check the number of indicators in each group after padj#
length(which(sample_data(tomcu_IS_sig_res_species)$group == "1")) #92 High/High
length(which(sample_data(tomcu_IS_sig_res_species)$group == "2")) #111 High/Low

#Order rows by group#
estuary_IS_sig_res_taxo_seqs_species <- arrange(estuary_IS_sig_res_taxo_seqs_species, group)
muck_IS_sig_res_taxo_seqs_species <- arrange(muck_IS_sig_res_taxo_seqs_species, group)
tomcu_IS_sig_res_taxo_seqs_species <- arrange(tomcu_IS_sig_res_taxo_seqs_species, group)


#Create lists of just significant speciess for each category#
unfiltered_IS_estuary_speciess <- subset(estuary_IS_sig_res_taxo_seqs_species, select=c(ESV.ID))
unfiltered_IS_muck_speciess <- subset(muck_IS_sig_res_taxo_seqs_species, select=c(ESV.ID))
unfiltered_IS_tomcu_speciess <- subset(tomcu_IS_sig_res_taxo_seqs_species, select=c(ESV.ID))


##Create indicator tables filtered of overlapping indicators with other two metadata categories##

#Filter out rows in the estuary IS object that are also in the tomcu and muck IS objects#
filtered_estuary_IS_sig_res_taxo_seqs_species <- anti_join(estuary_IS_sig_res_taxo_seqs_species, unfiltered_IS_tomcu_speciess, by = "ESV.ID")
filtered_estuary_IS_sig_res_taxo_seqs_species <- anti_join(filtered_estuary_IS_sig_res_taxo_seqs_species, unfiltered_IS_muck_speciess, by = "ESV.ID")

#Filter out rows in the muck IS object that are also in the estuary and tomcu IS objects#
filtered_muck_IS_sig_res_taxo_seqs_species<- anti_join(muck_IS_sig_res_taxo_seqs_species, unfiltered_IS_estuary_speciess, by = "ESV.ID")
filtered_muck_IS_sig_res_taxo_seqs_species <- anti_join(filtered_muck_IS_sig_res_taxo_seqs_species, unfiltered_IS_tomcu_speciess, by = "ESV.ID")

#Filter out rows in the tomcu IS object that are also in the estuary and muck IS objects#
filtered_tomcu_IS_sig_res_taxo_seqs_species <- anti_join(tomcu_IS_sig_res_taxo_seqs_species, unfiltered_IS_estuary_speciess, by = "ESV.ID")
filtered_tomcu_IS_sig_res_taxo_seqs_species <- anti_join(filtered_tomcu_IS_sig_res_taxo_seqs_species, unfiltered_IS_muck_speciess, by = "ESV.ID")

#Order rows by group#
filtered_estuary_IS_sig_res_taxo_seqs_species <- arrange(filtered_estuary_IS_sig_res_taxo_seqs_species, group)
filtered_muck_IS_sig_res_taxo_seqs_species <- arrange(filtered_muck_IS_sig_res_taxo_seqs_species, group)
filtered_tomcu_IS_sig_res_taxo_seqs_species <- arrange(filtered_tomcu_IS_sig_res_taxo_seqs_species, group)


#Make a csv file for each of the filtered tables#
write.csv(as.data.frame(filtered_estuary_IS_sig_res_taxo_seqs_species), file="filtered_estuary_IS_sig_res_taxo_seqs_species.csv")
write.csv(as.data.frame(filtered_muck_IS_sig_res_taxo_seqs_species), file="filtered_muck_IS_sig_res_taxo_seqs_species.csv")
write.csv(as.data.frame(filtered_tomcu_IS_sig_res_taxo_seqs_species), file="filtered_tomcu_IS_sig_res_taxo_seqs_species.csv")

#Determine the number of indicators in each category#
length(which(sample_data(filtered_estuary_IS_sig_res_taxo_seqs_species)$group == "1")) #106 IRL
length(which(sample_data(filtered_estuary_IS_sig_res_taxo_seqs_species)$group == "2")) #297 SLE
length(which(sample_data(filtered_muck_IS_sig_res_taxo_seqs_species)$group == "1")) #94 Muck
length(which(sample_data(filtered_muck_IS_sig_res_taxo_seqs_species)$group == "2")) #44 Not
length(which(sample_data(filtered_tomcu_IS_sig_res_taxo_seqs_species)$group == "1")) #13 High/High
length(which(sample_data(filtered_tomcu_IS_sig_res_taxo_seqs_species)$group == "2")) #11 High/Low


#Create lists of filtered significant speciess for each category#
filtered_IS_estuary_speciess <- subset(filtered_estuary_IS_sig_res_taxo_seqs_species, select=c(ESV.ID))
filtered_IS_muck_speciess <- subset(filtered_muck_IS_sig_res_taxo_seqs_species, select=c(ESV.ID))
filtered_IS_tomcu_speciess <- subset(filtered_tomcu_IS_sig_res_taxo_seqs_species, select=c(ESV.ID))

#Create phyloseq objects of the unfiltered indicators using character string version of significant speciess##
charac_unfiltered_IS_estuary_speciess <- unfiltered_IS_estuary_speciess[,"ESV.ID"]
charac_unfiltered_IS_estuary_speciess <- as.character(charac_unfiltered_IS_estuary_speciess)
UFISestuarySpeciesAspecies <- prune_taxa(charac_unfiltered_IS_estuary_speciess, RFLWSS_Y1_SpeciesAspecies)

charac_unfiltered_IS_muck_speciess <- unfiltered_IS_muck_speciess[,"ESV.ID"]
charac_unfiltered_IS_muck_speciess <- as.character(charac_unfiltered_IS_muck_speciess)
UFISmuckSpeciesAspecies <- prune_taxa(charac_unfiltered_IS_muck_speciess, RFLWSS_Y1_SpeciesAspecies)

charac_unfiltered_IS_tomcu_speciess <- unfiltered_IS_tomcu_speciess[,"ESV.ID"]
charac_unfiltered_IS_tomcu_speciess <- as.character(charac_unfiltered_IS_tomcu_speciess)
UFIStomcuSpeciesAspecies <- prune_taxa(charac_unfiltered_IS_tomcu_speciess, RFLWSS_Y1_SpeciesAspecies)

#Create phyloseq objects of the filtered indicators using character string version of significant speciess##
charac_filtered_IS_estuary_speciess <- filtered_IS_estuary_speciess[,"ESV.ID"]
charac_filtered_IS_estuary_speciess <- as.character(charac_filtered_IS_estuary_speciess)
FISestuarySpeciesAspecies <- prune_taxa(charac_filtered_IS_estuary_speciess, RFLWSS_Y1_SpeciesAspecies)

charac_filtered_IS_muck_speciess <- filtered_IS_muck_speciess[,"ESV.ID"]
charac_filtered_IS_muck_speciess <- as.character(charac_filtered_IS_muck_speciess)
FISmuckSpeciesAspecies <- prune_taxa(charac_filtered_IS_muck_speciess, RFLWSS_Y1_SpeciesAspecies)

charac_filtered_IS_tomcu_speciess <- filtered_IS_tomcu_speciess[,"ESV.ID"]
charac_filtered_IS_tomcu_speciess <- as.character(charac_filtered_IS_tomcu_speciess)
FIStomcuSpeciesAspecies <- prune_taxa(charac_filtered_IS_tomcu_speciess, RFLWSS_Y1_SpeciesAspecies)

##Create heatmaps to see overall trends, may take a while to load##

unfiltered_IS_estuary_species_heatmap <- plot_heatmap(UFDestuarySpeciesAspecies, taxa.order = charac_unfiltered_IS_estuary_speciess, sample.order= "Estuary", sample.label = "Estuary")
unfiltered_IS_estuary_species_heatmap

unfiltered_IS_muck_species_heatmap <- plot_heatmap(UFDmuckSpeciesAspecies, taxa.order = charac_unfiltered_IS_muck_speciess, sample.order= "Muck", sample.label = "Muck")
unfiltered_IS_muck_species_heatmap

unfiltered_IS_tomcu_species_heatmap <- plot_heatmap(UFDtomcuSpeciesAspecies, taxa.order = charac_unfiltered_IS_tomcu_speciess, sample.order= "LOI.Cu", sample.label = "LOI.Cu")
unfiltered_IS_tomcu_species_heatmap

filtered_IS_estuary_species_heatmap <- plot_heatmap(FDestuarySpeciesAspecies, taxa.order = charac_filtered_IS_estuary_speciess, sample.order= "Estuary", sample.label = "Estuary")
filtered_IS_estuary_species_heatmap

filtered_IS_muck_species_heatmap <- plot_heatmap(FDmuckSpeciesAspecies, taxa.order = charac_filtered_IS_muck_speciess, sample.order= "Muck", sample.label = "Muck")
filtered_IS_muck_species_heatmap

filtered_IS_tomcu_species_heatmap <- plot_heatmap(FDtomcuSpeciesAspecies, taxa.order = charac_filtered_IS_tomcu_speciess, sample.order= "LOI.Cu", sample.label = "LOI.Cu")
filtered_IS_tomcu_species_heatmap

###MAKE COMBO INDICATOR LISTS###

##Make "Combo" indicator lists from overlapping unfiltered DESeq2 and IS indicators##

#Keep rows in the estuary DESEQ object that are also in the IRLSLE IS object#
unfiltered_estuary_combo_sig_res_taxo_seqs_species <- merge(estuary_IS_sig_res_species,estuary_DESeq2_sig_res_taxo_seqs_species, by="ESV.ID")

#Keep rows in the muck DESEQ object that are also in the Muck IS object#
unfiltered_muck_combo_sig_res_taxo_seqs_species <- merge(muck_IS_sig_res_species,muck_DESeq2_sig_res_taxo_seqs_species, by="ESV.ID")

#Keep rows in the tomcu DESEQ object that are also in the LOICu IS object#
unfiltered_tomcu_combo_sig_res_taxo_seqs_species <- merge(tomcu_IS_sig_res_species,tomcu_DESeq2_sig_res_taxo_seqs_species, by="ESV.ID")

##Make "Combo" indicator lists from overlapping filtered DESeq2 and IS indicators##

#Create focused filtered IS tables without taxonomy or sequences#
filtered_estuary_IS_sig_res_species <- subset(filtered_estuary_IS_sig_res_taxo_seqs_species, select=c(ESV.ID,A,B,stat,p.value,group,padj))

filtered_muck_IS_sig_res_species <- subset(filtered_muck_IS_sig_res_taxo_seqs_species, select=c(ESV.ID,A,B,stat,p.value,group,padj))

filtered_tomcu_IS_sig_res_species <- subset(filtered_tomcu_IS_sig_res_taxo_seqs_species, select=c(ESV.ID,A,B,stat,p.value,group,padj))

#Keep rows in the estuary DESEQ object that are also in the IRLSLE IS object#
filtered_estuary_combo_sig_res_taxo_seqs_species <- merge(filtered_estuary_IS_sig_res_species,filtered_estuary_DESeq2_sig_res_taxo_seqs_species, by="ESV.ID")

#Keep rows in the muck DESEQ object that are also in the Muck IS object#
filtered_muck_combo_sig_res_taxo_seqs_species <- merge(filtered_muck_IS_sig_res_species,filtered_muck_DESeq2_sig_res_taxo_seqs_species, by="ESV.ID")

#Keep rows in the tomcu DESEQ object that are also in the LOICu IS object#
filtered_tomcu_combo_sig_res_taxo_seqs_species <- merge(filtered_tomcu_IS_sig_res_species,filtered_tomcu_DESeq2_sig_res_taxo_seqs_species, by="ESV.ID")

#Determine the number of indicators in each subcategory of the unfiltered Combo tables#
length(which(sample_data(unfiltered_estuary_combo_sig_res_taxo_seqs_species)$group == "1")) #176 IRL
length(which(sample_data(unfiltered_estuary_combo_sig_res_taxo_seqs_species)$group == "2")) #221 SLE
length(which(sample_data(unfiltered_muck_combo_sig_res_taxo_seqs_species)$group == "1")) #55 Muck
length(which(sample_data(unfiltered_muck_combo_sig_res_taxo_seqs_species)$group == "2")) #89 Not
length(which(sample_data(unfiltered_tomcu_combo_sig_res_taxo_seqs_species)$group == "1")) #56 High/High
length(which(sample_data(unfiltered_tomcu_combo_sig_res_taxo_seqs_species)$group == "2")) #109 High/Low

#Determine the number of indicators in each subcategory of the filtered Combo tables#
length(which(sample_data(filtered_estuary_combo_sig_res_taxo_seqs_species)$group == "1")) #65 IRL
length(which(sample_data(filtered_estuary_combo_sig_res_taxo_seqs_species)$group == "2")) #80 SLE
length(which(sample_data(filtered_muck_combo_sig_res_taxo_seqs_species)$group == "1")) #22 Muck
length(which(sample_data(filtered_muck_combo_sig_res_taxo_seqs_species)$group == "2")) #24 Not
length(which(sample_data(filtered_tomcu_combo_sig_res_taxo_seqs_species)$group == "1")) #2 High/High
length(which(sample_data(filtered_tomcu_combo_sig_res_taxo_seqs_species)$group == "2")) #8 High/Low


#Order unfiltered Combo table rows by group#
unfiltered_estuary_combo_sig_res_taxo_seqs_species <- arrange(unfiltered_estuary_combo_sig_res_taxo_seqs_species, group)
unfiltered_muck_combo_sig_res_taxo_seqs_species <- arrange(unfiltered_muck_combo_sig_res_taxo_seqs_species, group)
unfiltered_tomcu_combo_sig_res_taxo_seqs_species <- arrange(unfiltered_tomcu_combo_sig_res_taxo_seqs_species, group)

#Order filtered Combo table rows by group#
filtered_estuary_combo_sig_res_taxo_seqs_species <- arrange(filtered_estuary_combo_sig_res_taxo_seqs_species, group)
filtered_muck_combo_sig_res_taxo_seqs_species <- arrange(filtered_muck_combo_sig_res_taxo_seqs_species, group)
filtered_tomcu_combo_sig_res_taxo_seqs_species <- arrange(filtered_tomcu_combo_sig_res_taxo_seqs_species, group)

#Make a csv file for each of the unfiltered Combo tables#
write.csv(as.data.frame(unfiltered_estuary_combo_sig_res_taxo_seqs_species), file="unfiltered_estuary_combo_sig_res_taxo_seqs_species.csv")
write.csv(as.data.frame(unfiltered_muck_combo_sig_res_taxo_seqs_species), file="unfiltered_muck_combo_sig_res_taxo_seqs_species.csv")
write.csv(as.data.frame(unfiltered_tomcu_combo_sig_res_taxo_seqs_species), file="unfiltered_tomcu_combo_sig_res_taxo_seqs_species.csv")

#Make a csv file for each of the filtered Combo tables#
write.csv(as.data.frame(filtered_estuary_combo_sig_res_taxo_seqs_species), file="filtered_estuary_combo_sig_res_taxo_seqs_species.csv")
write.csv(as.data.frame(filtered_muck_combo_sig_res_taxo_seqs_species), file="filtered_muck_combo_sig_res_taxo_seqs_species.csv")
write.csv(as.data.frame(filtered_tomcu_combo_sig_res_taxo_seqs_species), file="filtered_tomcu_combo_sig_res_taxo_seqs_species.csv")

#Create lists of unfiltered significant speciess for each category#
unfiltered_Combo_estuary_speciess <- subset(unfiltered_estuary_combo_sig_res_taxo_seqs_species, select=c(ESV.ID))
unfiltered_Combo_muck_speciess <- subset(unfiltered_muck_combo_sig_res_taxo_seqs_species, select=c(ESV.ID))
unfiltered_Combo_tomcu_speciess <- subset(unfiltered_tomcu_combo_sig_res_taxo_seqs_species, select=c(ESV.ID))

#Create lists of filtered significant speciess for each category#
filtered_Combo_estuary_speciess <- subset(filtered_estuary_combo_sig_res_taxo_seqs_species, select=c(ESV.ID))
filtered_Combo_muck_speciess <- subset(filtered_muck_combo_sig_res_taxo_seqs_species, select=c(ESV.ID))
filtered_Combo_tomcu_speciess <- subset(filtered_tomcu_combo_sig_res_taxo_seqs_species, select=c(ESV.ID))

##Create phyloseq objects for Combo indicators##

#Create phyloseq objects of the unfiltered indicators using character string version of significant speciess##
charac_unfiltered_Combo_estuary_speciess <- unfiltered_Combo_estuary_speciess[,"ESV.ID"]
charac_unfiltered_Combo_estuary_speciess <- as.character(charac_unfiltered_Combo_estuary_speciess)
UFCestuarySpeciesAspecies <- prune_taxa(charac_unfiltered_Combo_estuary_speciess, RFLWSS_Y1_SpeciesAspecies)

charac_unfiltered_Combo_muck_speciess <- unfiltered_Combo_muck_speciess[,"ESV.ID"]
charac_unfiltered_Combo_muck_speciess <- as.character(charac_unfiltered_Combo_muck_speciess)
UFCmuckSpeciesAspecies <- prune_taxa(charac_unfiltered_Combo_muck_speciess, RFLWSS_Y1_SpeciesAspecies)

charac_unfiltered_Combo_tomcu_speciess <- unfiltered_Combo_tomcu_speciess[,"ESV.ID"]
charac_unfiltered_Combo_tomcu_speciess <- as.character(charac_unfiltered_Combo_tomcu_speciess)
UFCtomcuSpeciesAspecies <- prune_taxa(charac_unfiltered_Combo_tomcu_speciess, RFLWSS_Y1_SpeciesAspecies)

#Create phyloseq objects of the filtered indicators using character string version of significant speciess##
charac_filtered_Combo_estuary_speciess <- filtered_Combo_estuary_speciess[,"ESV.ID"]
charac_filtered_Combo_estuary_speciess <- as.character(charac_filtered_Combo_estuary_speciess)
FCestuarySpeciesAspecies <- prune_taxa(charac_filtered_Combo_estuary_speciess, RFLWSS_Y1_SpeciesAspecies)

charac_filtered_Combo_muck_speciess <- filtered_Combo_muck_speciess[,"ESV.ID"]
charac_filtered_Combo_muck_speciess <- as.character(charac_filtered_Combo_muck_speciess)
FCmuckSpeciesAspecies <- prune_taxa(charac_filtered_Combo_muck_speciess, RFLWSS_Y1_SpeciesAspecies)

charac_filtered_Combo_tomcu_speciess <- filtered_Combo_tomcu_speciess[,"ESV.ID"]
charac_filtered_Combo_tomcu_speciess <- as.character(charac_filtered_Combo_tomcu_speciess)
FCtomcuSpeciesAspecies <- prune_taxa(charac_filtered_Combo_tomcu_speciess, RFLWSS_Y1_SpeciesAspecies)

##Create heatmaps to see overall trends, may take a while to load##

unfiltered_Combo_estuary_species_heatmap <- plot_heatmap(UFDestuarySpeciesAspecies, taxa.order = charac_unfiltered_Combo_estuary_speciess, sample.order= "Estuary", sample.label = "Estuary")
unfiltered_Combo_estuary_species_heatmap

unfiltered_Combo_muck_species_heatmap <- plot_heatmap(UFDmuckSpeciesAspecies, taxa.order = charac_unfiltered_Combo_muck_speciess, sample.order= "Muck", sample.label = "Muck")
unfiltered_Combo_muck_species_heatmap

unfiltered_Combo_tomcu_species_heatmap <- plot_heatmap(UFDtomcuSpeciesAspecies, taxa.order = charac_unfiltered_Combo_tomcu_speciess, sample.order= "LOI.Cu", sample.label = "LOI.Cu")
unfiltered_Combo_tomcu_species_heatmap

filtered_Combo_estuary_species_heatmap <- plot_heatmap(FDestuarySpeciesAspecies, taxa.order = charac_filtered_Combo_estuary_speciess, sample.order= "Estuary", sample.label = "Estuary")
filtered_Combo_estuary_species_heatmap

filtered_Combo_muck_species_heatmap <- plot_heatmap(FDmuckSpeciesAspecies, taxa.order = charac_filtered_Combo_muck_speciess, sample.order= "Muck", sample.label = "Muck")
filtered_Combo_muck_species_heatmap

filtered_Combo_tomcu_species_heatmap <- plot_heatmap(FDtomcuSpeciesAspecies, taxa.order = charac_filtered_Combo_tomcu_speciess, sample.order= "LOI.Cu", sample.label = "LOI.Cu")
filtered_Combo_tomcu_species_heatmap

###INDICATOR EFFECTIVENESS METHOD###

##NOTES ON METHOD##
#"affected' = subcategory that is more affected by the stressor you are      concerned about. For example SLE for the Estuary category becasue it has     more freshwater discharges, 3MC for the Muck category because it has more    muck characteristics, HiHi for the TOM/Cu category becasue it has more       copper#

#"non.affected" = IRL, 0MC, HiLo#

#microbially.predicted "affected" sample = samples in the partitioning around medoids (PAM) cluster with the most metadata.defined "affected" samples#

#microbially.predicted "non.affected" sample = samples in the PAM cluster with the most metadata.defined "non.affected" samples#

#Overall Idea: Use PAM clustering to split samples between the "affected" and "non.affected" clusters based upon the indicators (microbially.predicted) and see how well this classification matches the metadata.defined classification by using the product of specificity and sensitivity#

#Copy results to an Excel Sheet for further statistical testing, each combo of factors should have four percentage types (Product is Sensitivity x Specificity), for example:#
#Indicator_Test	Taxonomic_Level	Filtering_Status	Metadata_Tested	Percent_Type	Percentage#
#Indicspecies   species             Unfiltered        Estuary         Total         98.78#
#Indicspecies   species             Unfiltered        Estuary         Sensitivity   96.23#
#Indicspecies   species             Unfiltered        Estuary         Specificity   99.15#
#Indicspecies   species             Unfiltered        Estuary         Product       95.41#


###ORIGINAL DATASETS
###Split samples of orginal datasets into two clusters to see how it will compare to splitting based upon indicators###

##Test on Muck and Estuary categories at species level##

#Remove any samples with no speciess from target phyloseq#
nsamples(RFLWSS_Y1_SpeciesAspecies) #90
RFLWSS_Y1_SpeciesAspecies = prune_samples(sample_sums(RFLWSS_Y1_SpeciesAspecies)>=1, RFLWSS_Y1_SpeciesAspecies)
nsamples(RFLWSS_Y1_SpeciesAspecies) #90

#Make a matrix with the sample_data and relabel subcategories to make metadata.defined column#
Original_Est_Muck_species <- as.matrix(sample_data(RFLWSS_Y1_SpeciesAspecies))
Original_Est_Muck_species[Original_Est_Muck_species ==  "Muck"] <- 2
Original_Est_Muck_species[Original_Est_Muck_species ==  "Not"] <- 1
Original_Est_Muck_species[Original_Est_Muck_species ==  "Mucky"] <- 2
Original_Est_Muck_species[Original_Est_Muck_species ==  "Muckish"] <- 2
Original_Est_Muck_species[Original_Est_Muck_species ==  "IRL"] <- 1
Original_Est_Muck_species[Original_Est_Muck_species ==  "SLE"] <- 2
Original_Est_Muck_species<-as.data.frame(Original_Est_Muck_species)


#Create Bray.Curtis distance matrix from species table#
bc <- vegdist(t(otu_table(RFLWSS_Y1_SpeciesAspecies)))

#Conduct pam clustering based upon 2 clusters#
pam.res <- pam(bc,2, diss=TRUE)

#Add a microbially.predicted column (Orginal_species) to matrix above for comparison to metadata.defined column#
Original_Est_Muck_species <- cbind(Original_Est_Muck_species, Orginal_species = pam.res$cluster)

##Test Estuary category##

#Test  to see how well the metadata.defined and microbially.predicted columns match one another aka Total Efficiency; if the number is below 0.5, then switch the numbers in the metadata.defiend column#

#Test for Estuary Effeciency# 0.9555556
sum(Original_Est_Muck_species$Estuary == Original_Est_Muck_species$Orginal_species)/nrow(Original_Est_Muck_species)

##Test indicator sensitivity aka the true postives (TP)/ (TP + false negatives (FN))##
#TP = metadata.defined "affected" sample correctly placed in the microbiall   y.predicted "affected" cluster#
#FN = metadata.defined "affected" sample incorrectly placed in the           microbially.predicted "non.affected" cluster#
#TP+FN = all metadata.defined "affected" samples

#Test for SLE Effeciency aka Sensitvity# 1
Orginal_SLE_species <- filter(Original_Est_Muck_species, Estuary == "2")
sum(Orginal_SLE_species$Estuary == Orginal_SLE_species$Orginal_species)/nrow(Orginal_SLE_species)

##Test indicator specificity aka the true negatives (TN)/ (TN + false positives (FP))##
#TP = metadata.defined "affected" sample correctly placed in the microbiall   y.predicted "affected" cluster#
#FN = metadata.defined "affected" sample incorrectly placed in the           microbially.predicted "non.affected" cluster#
#TP+FN = all metadata.defined "affected" samples

#Test for IRL Effeciency aka Specificity# 0.9393939
Orginal_IRL_species <- filter(Original_Est_Muck_species, Estuary == "1")
sum(Orginal_IRL_species$Estuary == Orginal_IRL_species$Orginal_species)/nrow(Orginal_IRL_species)

#PCoA with Color by Estuary#
TotalPCoA <- ordinate(RFLWSS_Y1_SpeciesAspecies,"PCoA")
p = plot_ordination(RFLWSS_Y1_SpeciesAspecies, TotalPCoA, color="Estuary") + 
  ggtitle("Total speciess separated by Estuary Category PCoA with Bray.Curtis distance") +
  geom_text(aes(label=Original_Est_Muck_species$Orginal_species),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p

##Test Muck category##

#Test for Muck Total Efficiency# 0.7333333
sum(Original_Est_Muck_species$Muck == Original_Est_Muck_species$Orginal_species)/nrow(Original_Est_Muck_species)

#Test for 3MC Effeciency aka Sensitivity# 0.5666667
Orginal_3MC_species <- filter(Original_Est_Muck_species, Muck == "2")
sum(Orginal_3MC_species$Muck == Orginal_3MC_species$Orginal_species)/nrow(Orginal_3MC_species)

#Test for 0MC Effeciency aka Specificity# 0.8166667
Orginal_0MC_species <- filter(Original_Est_Muck_species, Muck == "1")
sum(Orginal_0MC_species$Muck == Orginal_0MC_species$Orginal_species)/nrow(Orginal_0MC_species)

#PCoA with Color by Muck#
TotalPCoA <- ordinate(RFLWSS_Y1_SpeciesAspecies,"PCoA")
p = plot_ordination(RFLWSS_Y1_SpeciesAspecies, TotalPCoA, color="Muck") + 
  ggtitle("Total speciess separated by Muck Category PCoA with Bray.Curtis distance") +
  geom_text(aes(label=Original_Est_Muck_species$Orginal_species),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p + annotate("text", x = 0.3, y = 0.3, label = c("SN = 0.44, SP=0.84"))+ scale_color_discrete(name="Muck Characteristics", limits=c("Muck", "Mucky", "Muckish", "Not"), labels=c("3", "2", "1", "0"))

##Test TOM/Cu category##

#Remove any samples with no speciess from target phyloseq#
nsamples(RFLWSS_Y1_HHvsHL_SpeciesAspecies) #29
RFLWSS_Y1_HHvsHL_SpeciesAspecies = prune_samples(sample_sums(RFLWSS_Y1_HHvsHL_SpeciesAspecies)>=1, RFLWSS_Y1_HHvsHL_SpeciesAspecies)
nsamples(RFLWSS_Y1_HHvsHL_SpeciesAspecies) #29

#Make a matrix with the sample_data and relabel subcategories to make metadata.defined column#
Original_TOMCu_species <- as.matrix(sample_data(RFLWSS_Y1_HHvsHL_SpeciesAspecies))
Original_TOMCu_species[Original_TOMCu_species ==  "High.High"] <- 1
Original_TOMCu_species[Original_TOMCu_species ==  "High.Low"] <- 2
Original_TOMCu_species<-as.data.frame(Original_TOMCu_species)

#Create Bray.Curtis distance matrix from species table#
bc <- vegdist(t(otu_table(RFLWSS_Y1_HHvsHL_SpeciesAspecies)))

#Conduct PAM clustering based upon 2 clusters#
pam.res <- pam(bc,2, diss=TRUE)

#Add a microbially.predicted column to matrix above for comparison to metadata.defined column#
Original_TOMCu_species <- cbind(Original_TOMCu_species, Orginal_species = pam.res$cluster)

#Determine Total TOM/Cu Efficiency# 0.6551724
sum(Original_TOMCu_species$LOI.Cu == Original_TOMCu_species$Orginal_species)/nrow(Original_TOMCu_species)

#Determine HiHi Efficiency  aka Sensitivity# 0.5454545
Original_TOMCu_speciesHH <- filter(Original_TOMCu_species, LOI.Cu == "1") 
sum(Original_TOMCu_speciesHH$LOI.Cu == Original_TOMCu_speciesHH$Orginal_species)/nrow(Original_TOMCu_speciesHH)

#Determine HiLo Efficiency aka Specificity# 0.7222222
Original_TOMCu_speciesHL <- filter(Original_TOMCu_species, LOI.Cu == "2") 
sum(Original_TOMCu_speciesHL$LOI.Cu == Original_TOMCu_speciesHL$Orginal_species)/nrow(Original_TOMCu_speciesHL)

#PCoA with Color by TOM/Cu#
TotalPCoA <- ordinate(RFLWSS_Y1_HHvsHL_SpeciesAspecies,"PCoA")
p = plot_ordination(RFLWSS_Y1_HHvsHL_SpeciesAspecies, TotalPCoA, color="LOI.Cu") + 
  ggtitle("Total speciess separated by LOI.Cu Category PCoA with Bray.Curtis distance") +
  geom_text(aes(label=Original_TOMCu_species$Orginal_species),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p

###TEST INDICSPECIES ESTUARY INDICATORS###

##Test unfiltered Estuary speciess##

#Remove any samples with no speciess from target phyloseq#
nsamples(UFISestuarySpeciesAspecies) #90
UFISestuarySpeciesAspecies = prune_samples(sample_sums(UFISestuarySpeciesAspecies)>=1, UFISestuarySpeciesAspecies)
nsamples(UFISestuarySpeciesAspecies) #90

#Make a matrix with the sample_data and relabel subcategories to make metadata.defined column#
UF_IS_Est_species <- as.matrix(sample_data(UFISestuarySpeciesAspecies))
UF_IS_Est_species[UF_IS_Est_species ==  "SLE"] <- 2
UF_IS_Est_species[UF_IS_Est_species ==  "IRL"] <- 1
UF_IS_Est_species<-as.data.frame(UF_IS_Est_species)

#Create Bray.Curtis distance matrix from species table#
bc <- vegdist(t(otu_table(UFISestuarySpeciesAspecies)))

#Conduct PAM clustering based upon 2 clusters#
pam.res <- pam(bc,2, diss=TRUE)

#Add a microbially.predicted column to matrix above for comparison to metadata.defined column#
UF_IS_Est_species <- cbind(UF_IS_Est_species, MP_UFISEstspecies = pam.res$cluster)

#Test Total Estuary Efficiency# 1
sum(UF_IS_Est_species$Estuary == UF_IS_Est_species$MP_UFISEstspecies)/nrow(UF_IS_Est_species)

#Test SLE Efficiency aka Sensitvity# 1
UF_IS_SLE_species <- filter(UF_IS_Est_species, Estuary == "2")
sum(UF_IS_SLE_species$Estuary == UF_IS_SLE_species$MP_UFISEstspecies)/nrow(UF_IS_SLE_species)

#Test IRL Efficiency aka Specificity# 1
UF_IS_IRL_species <- filter(UF_IS_Est_species, Estuary == "1")
sum(UF_IS_IRL_species$Estuary == UF_IS_IRL_species$MP_UFISEstspecies)/nrow(UF_IS_IRL_species)

TotalPCoA <- ordinate(UFISestuarySpeciesAspecies,"PCoA")
p = plot_ordination(UFISestuarySpeciesAspecies, TotalPCoA, color="Estuary") + 
  ggtitle("Unfiltered IS Estuary species Indicators PCoA with Bray.Curtis distance") +
  geom_text(aes(label=UF_IS_Est_species$MP_UFISEstspecies),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p


##Test filtered Estuary speciess##

#Remove any samples with no speciess from target phyloseq#
nsamples(FISestuarySpeciesAspecies) #90
FISestuarySpeciesAspecies = prune_samples(sample_sums(FISestuarySpeciesAspecies)>=1, FISestuarySpeciesAspecies)
nsamples(FISestuarySpeciesAspecies) #90

#Make a matrix with the sample_data and relabel subcategories to make metadata.defined column#
F_IS_Est_species <- as.matrix(sample_data(FISestuarySpeciesAspecies))
F_IS_Est_species[F_IS_Est_species ==  "SLE"] <- 2
F_IS_Est_species[F_IS_Est_species ==  "IRL"] <- 1
F_IS_Est_species<-as.data.frame(F_IS_Est_species)


#Create Bray.Curtis distance matrix from species table#
bc <- vegdist(t(otu_table(FISestuarySpeciesAspecies)))

#Conduct PAM clustering based upon 2 clusters#
pam.res <- pam(bc,2, diss=TRUE)

#Add a microbially.predicted column to matrix above for comparison to metadata.defined column#
F_IS_Est_species <- cbind(F_IS_Est_species, MP_FISEstspecies = pam.res$cluster)

#Test Total Estuary Efficiency# 0.9777778
sum(F_IS_Est_species$Estuary == F_IS_Est_species$MP_FISEstspecies)/nrow(F_IS_Est_species)
#Test SLE Efficiency aka Sensitivity# 0.9166667
F_IS_SLE_species <- filter(F_IS_Est_species, Estuary == "2")
sum(F_IS_SLE_species$Estuary == F_IS_SLE_species$MP_FISEstspecies)/nrow(F_IS_SLE_species)
#Test IRL Efficiency aka Specificity# 1
F_IS_IRL_species <- filter(F_IS_Est_species, Estuary == "1")
sum(F_IS_IRL_species$Estuary == F_IS_IRL_species$MP_FISEstspecies)/nrow(F_IS_IRL_species)

TotalPCoA <- ordinate(FISestuarySpeciesAspecies,"PCoA")
p = plot_ordination(FISestuarySpeciesAspecies, TotalPCoA, color="Estuary") + 
  ggtitle("Filtered IS Estuary species Indicators PCoA with Bray.Curtis distance") +
  geom_text(aes(label=F_IS_Est_species$MP_FISEstspecies),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p


###TEST DESEQ2 ESTUARY INDICATORS###

##Test unfiltered Estuary speciess##

#Remove any samples with no speciess from target phyloseq#
nsamples(UFDestuarySpeciesAspecies) #90
UFDestuarySpeciesAspecies = prune_samples(sample_sums(UFDestuarySpeciesAspecies)>=1, UFDestuarySpeciesAspecies)
nsamples(UFDestuarySpeciesAspecies) #90

#Make a matrix with the sample_data and relabel subcategories to make metadata.defined column#
UF_D_Est_species <- as.matrix(sample_data(UFDestuarySpeciesAspecies))
UF_D_Est_species[UF_D_Est_species ==  "SLE"] <- 2
UF_D_Est_species[UF_D_Est_species ==  "IRL"] <- 1
UF_D_Est_species<-as.data.frame(UF_D_Est_species)


#Create Bray.Curtis distance matrix from species table#
bc <- vegdist(t(otu_table(UFDestuarySpeciesAspecies)))

#Conduct PAM clustering based upon 2 clusters#
pam.res <- pam(bc,2, diss=TRUE)

#Add a microbially.predicted column to matrix above for comparison to metadata.defined column#
UF_D_Est_species <- cbind(UF_D_Est_species, MP_UFDEstspecies = pam.res$cluster)

#Test Total Estuary Efficiency# 0.9777778
sum(UF_D_Est_species$Estuary == UF_D_Est_species$MP_UFDEstspecies)/nrow(UF_D_Est_species)

#Test SLE Efficiency aka Sensitivity# 1
UF_D_SLE_species <- filter(UF_D_Est_species, Estuary == "2")
sum(UF_D_SLE_species$Estuary == UF_D_SLE_species$MP_UFDEstspecies)/nrow(UF_D_SLE_species)

#Test IRL Efficiency aka Specificity# 0.969697
UF_D_IRL_species <- filter(UF_D_Est_species, Estuary == "1")
sum(UF_D_IRL_species$Estuary == UF_D_IRL_species$MP_UFDEstspecies)/nrow(UF_D_IRL_species)


TotalPCoA <- ordinate(UFDestuarySpeciesAspecies,"PCoA")
p = plot_ordination(UFDestuarySpeciesAspecies, TotalPCoA, color="Estuary") + 
  ggtitle("Unfiltered DESeq2 Estuary species Indicators PCoA with Bray.Curtis distance") +
  geom_text(aes(label=UF_D_Est_species$MP_UFDEstspecies),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p


#Test filtered DESeq2 Estuary speciess#

#Remove any samples with no speciess from target phyloseq#
nsamples(FDestuarySpeciesAspecies) #90
FDestuarySpeciesAspecies = prune_samples(sample_sums(FDestuarySpeciesAspecies)>=1, FDestuarySpeciesAspecies)
nsamples(FDestuarySpeciesAspecies) #90

#Make a matrix with the sample_data and relabel subcategories to make metadata.defined column#
F_D_Est_species <- as.matrix(sample_data(FDestuarySpeciesAspecies))
F_D_Est_species[F_D_Est_species ==  "SLE"] <- 2
F_D_Est_species[F_D_Est_species ==  "IRL"] <- 1
F_D_Est_species<-as.data.frame(F_D_Est_species)


#Create Bray.Curtis distance matrix from species table#
bc <- vegdist(t(otu_table(FDestuarySpeciesAspecies)))

#Conduct PAM clustering based upon 2 clusters#
pam.res <- pam(bc,2, diss=TRUE)

#Add a microbially.predicted column to matrix above for comparison to metadata.defined column#
F_D_Est_species <- cbind(F_D_Est_species, MP_FDEstspecies = pam.res$cluster)

#Test Total Estuary Efficiency# 0.6222222
sum(F_D_Est_species$Estuary == F_D_Est_species$MP_FDEstspecies)/nrow(F_D_Est_species)
#Test SLE Efficiency aka Sensitivity# 0.8333333
F_D_SLE_species <- filter(F_D_Est_species, Estuary == "2")
sum(F_D_SLE_species$Estuary == F_D_SLE_species$MP_FDEstspecies)/nrow(F_D_SLE_species)
#Test IRL Efficiency aka Specificity# 0.5454545
F_D_IRL_species <- filter(F_D_Est_species, Estuary == "1")
sum(F_D_IRL_species$Estuary == F_D_IRL_species$MP_FDEstspecies)/nrow(F_D_IRL_species)

TotalPCoA <- ordinate(FDestuarySpeciesAspecies,"PCoA")
p = plot_ordination(FDestuarySpeciesAspecies, TotalPCoA, color="Estuary") + 
  ggtitle("Filtered DESeq2 Estuary species Indicators PCoA with Bray.Curtis distance") +
  geom_text(aes(label=F_D_Est_species$MP_FDEstspecies),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p

###TEST COMBO ESTUARY INDICATORS###

##Test unfiltered Estuary speciess##

#Remove any samples with no speciess from target phyloseq#
nsamples(UFCestuarySpeciesAspecies) #90
UFCestuarySpeciesAspecies = prune_samples(sample_sums(UFCestuarySpeciesAspecies)>=1, UFCestuarySpeciesAspecies)
nsamples(UFCestuarySpeciesAspecies) #90

#Make a matrix with the sample_data and relabel subcategories to make metadata.defined column#
UF_C_Est_species <- as.matrix(sample_data(UFCestuarySpeciesAspecies))
UF_C_Est_species[UF_C_Est_species ==  "SLE"] <- 2
UF_C_Est_species[UF_C_Est_species ==  "IRL"] <- 1
UF_C_Est_species<-as.data.frame(UF_C_Est_species)


#Create Bray.Curtis distance matrix from species table#
bc <- vegdist(t(otu_table(UFCestuarySpeciesAspecies)))

#Conduct PAM clustering based upon 2 clusters#
pam.res <- pam(bc,2, diss=TRUE)

#Add a microbially.predicted column to matrix above for comparison to metadata.defined column#
UF_C_Est_species <- cbind(UF_C_Est_species, MP_UFCEstspecies = pam.res$cluster)

#Test Total Estuary Efficiency# 1
sum(UF_C_Est_species$Estuary == UF_C_Est_species$MP_UFCEstspecies)/nrow(UF_C_Est_species)

#Test SLE Efficiency aka Sensitivity# 1
UF_C_SLE_species <- filter(UF_C_Est_species, Estuary == "2")
sum(UF_C_SLE_species$Estuary == UF_C_SLE_species$MP_UFCEstspecies)/nrow(UF_C_SLE_species)

#Test IRL Efficiency aka Specificity# 1
UF_C_IRL_species <- filter(UF_C_Est_species, Estuary == "1")
sum(UF_C_IRL_species$Estuary == UF_C_IRL_species$MP_UFCEstspecies)/nrow(UF_C_IRL_species)

TotalPCoA <- ordinate(UFCestuarySpeciesAspecies,"PCoA")
p = plot_ordination(UFCestuarySpeciesAspecies, TotalPCoA, color="Estuary") + 
  ggtitle("Unfiltered Final Estuary species Indicators PCoA with Bray.Curtis distance") +
  geom_text(aes(label=UF_C_Est_species$MP_UFCEstspecies),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p


##Test filtered Estuary speciess##

#Remove any samples with no speciess from target phyloseq#
nsamples(FCestuarySpeciesAspecies) #90
FCestuarySpeciesAspecies = prune_samples(sample_sums(FCestuarySpeciesAspecies)>=1, FCestuarySpeciesAspecies)
nsamples(FCestuarySpeciesAspecies) #90

#Make a matrix with the sample_data and relabel subcategories to make metadata.defined column#
F_C_Est_species <- as.matrix(sample_data(FCestuarySpeciesAspecies))
F_C_Est_species[F_C_Est_species ==  "SLE"] <- 2
F_C_Est_species[F_C_Est_species ==  "IRL"] <- 1
F_C_Est_species<-as.data.frame(F_C_Est_species)


#Create Bray.Curtis distance matrix from species table#
bc <- vegdist(t(otu_table(FCestuarySpeciesAspecies)))

#Conduct PAM clustering based upon 2 clusters#
pam.res <- pam(bc,2, diss=TRUE)

#Add a microbially.predicted column to matrix above for comparison to metadata.defined column#
F_C_Est_species <- cbind(F_C_Est_species, MP_FCEstspecies = pam.res$cluster)

#Test Total Estuary Efficiency# 0.9666667
sum(F_C_Est_species$Estuary == F_C_Est_species$MP_FCEstspecies)/nrow(F_C_Est_species)

#Test SLE Efficiency aka Sensitivity# 0.875
F_C_SLE_species <- filter(F_C_Est_species, Estuary == "2")
sum(F_C_SLE_species$Estuary == F_C_SLE_species$MP_FCEstspecies)/nrow(F_C_SLE_species)

#Test IRL Efficiency aka Sensitivity# 1
F_C_IRL_species <- filter(F_C_Est_species, Estuary == "1")
sum(F_C_IRL_species$Estuary == F_C_IRL_species$MP_FCEstspecies)/nrow(F_C_IRL_species)

TotalPCoA <- ordinate(FCestuarySpeciesAspecies,"PCoA")
p = plot_ordination(FCestuarySpeciesAspecies, TotalPCoA, color="Estuary") + 
  ggtitle("Filtered Final Estuary species Indicators PCoA with Bray.Curtis distance") +
  geom_text(aes(label=F_C_Est_species$MP_FCEstspecies),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p


###TEST INDICSPECIES MUCK INDICATORS###

##Test unfiltered Muck speciess##

#Remove any samples with no speciess from target phyloseq#
nsamples(UFISmuckSpeciesAspecies) #90
UFISmuckSpeciesAspecies = prune_samples(sample_sums(UFISmuckSpeciesAspecies)>=1, UFISmuckSpeciesAspecies)
nsamples(UFISmuckSpeciesAspecies) #90

#Make a matrix with the sample_data and relabel subcategories to make metadata.defined column#
UF_IS_Muck_species <- as.matrix(sample_data(UFISmuckSpeciesAspecies))
UF_IS_Muck_species[UF_IS_Muck_species ==  "Muck"] <- 2
UF_IS_Muck_species[UF_IS_Muck_species ==  "Not"] <- 1
UF_IS_Muck_species[UF_IS_Muck_species ==  "Mucky"] <- 2
UF_IS_Muck_species[UF_IS_Muck_species ==  "Muckish"] <- 2
UF_IS_Muck_species<-as.data.frame(UF_IS_Muck_species)


#Create Bray.Curtis distance matrix from species table#
bc <- vegdist(t(otu_table(UFISmuckSpeciesAspecies)))

#Conduct PAM clustering based upon 2 clusters#
pam.res <- pam(bc,2, diss=TRUE)

#Add a microbially.predicted column to matrix above for comparison to metadata.defined column#
UF_IS_Muck_species <- cbind(UF_IS_Muck_species, MP_UFISMuckspecies = pam.res$cluster)

#Test Total Muck Efficiency# 0.8444444
sum(UF_IS_Muck_species$Muck == UF_IS_Muck_species$MP_UFISMuckspecies)/nrow(UF_IS_Muck_species)

#Test 3MC Efficiency aka Sensitivity# 0.9333333
UF_IS_3MC_species <- filter(UF_IS_Muck_species, Muck == "2")
sum(UF_IS_3MC_species$Muck == UF_IS_3MC_species$MP_UFISMuckspecies)/nrow(UF_IS_3MC_species)

#Test 0MC Efficiency aka Specificity# 0.8
UF_IS_0MC_species <- filter(UF_IS_Muck_species, Muck == "1")
sum(UF_IS_0MC_species$Muck == UF_IS_0MC_species$MP_UFISMuckspecies)/nrow(UF_IS_0MC_species)


TotalPCoA <- ordinate(UFISmuckSpeciesAspecies,"PCoA")
p = plot_ordination(UFISmuckSpeciesAspecies, TotalPCoA, color="Muck") + 
  ggtitle("Unfiltered IS Muck species Indicators PCoA with Bray.Curtis distance") +
  geom_text(aes(label=UF_IS_Muck_species$MP_UFISMuckspecies),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p


##Test filtered Muck speciess##

#Remove any samples with no speciess from target phyloseq#
nsamples(FISmuckSpeciesAspecies) #90
FISmuckSpeciesAspecies = prune_samples(sample_sums(FISmuckSpeciesAspecies)>=1, FISmuckSpeciesAspecies)
nsamples(FISmuckSpeciesAspecies) #90

#Make a matrix with the sample_data and relabel subcategories to make metadata.defined column#
F_IS_Muck_species <- as.matrix(sample_data(FISmuckSpeciesAspecies))
F_IS_Muck_species[F_IS_Muck_species ==  "Muck"] <- 2
F_IS_Muck_species[F_IS_Muck_species ==  "Not"] <- 1
F_IS_Muck_species[F_IS_Muck_species ==  "Mucky"] <- 2
F_IS_Muck_species[F_IS_Muck_species ==  "Muckish"] <- 2
F_IS_Muck_species<-as.data.frame(F_IS_Muck_species)


#Create Bray.Curtis distance matrix from species table#
bc <- vegdist(t(otu_table(FISmuckSpeciesAspecies)))

#Conduct PAM clustering based upon 2 clusters#
pam.res <- pam(bc,2, diss=TRUE)

#Add a microbially.predicted column to matrix above for comparison to metadata.defined column#
F_IS_Muck_species <- cbind(F_IS_Muck_species, MP_FISMuckspecies = pam.res$cluster)

#Test Total Muck Efficiency# 0.9222222
sum(F_IS_Muck_species$Muck == F_IS_Muck_species$MP_FISMuckspecies)/nrow(F_IS_Muck_species)
#Test 3MC Efficiency aka  Sensitivity# 0.8
F_IS_3MC_species <- filter(F_IS_Muck_species, Muck == "2")
sum(F_IS_3MC_species$Muck == F_IS_3MC_species$MP_FISMuckspecies)/nrow(F_IS_3MC_species)
#Test 0MC Efficiency aka Specificity# 0.9333333
F_IS_0MC_species <- filter(F_IS_Muck_species, Muck == "1")
sum(F_IS_0MC_species$Muck == F_IS_0MC_species$MP_FISMuckspecies)/nrow(F_IS_0MC_species)

TotalPCoA <- ordinate(FISmuckSpeciesAspecies,"PCoA")
p = plot_ordination(FISmuckSpeciesAspecies, TotalPCoA, color="Muck") + 
  ggtitle("Filtered IS Muck species Indicators PCoA with Bray.Curtis distance") +
  geom_text(aes(label=F_IS_Muck_species$MP_FISMuckspecies),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p

###TEST DESEQ2 MUCK INDICATORS###

#Test unfiltered Muck speciess##

#Remove any samples with no speciess from target phyloseq#
nsamples(UFDmuckSpeciesAspecies) #90
UFDmuckSpeciesAspecies = prune_samples(sample_sums(UFDmuckSpeciesAspecies)>=1, UFDmuckSpeciesAspecies)
nsamples(UFDmuckSpeciesAspecies) #90

#Make a matrix with the sample_data and relabel subcategories to make metadata.defined column#
UF_D_Muck_species <- as.matrix(sample_data(UFDmuckSpeciesAspecies))
UF_D_Muck_species[UF_D_Muck_species ==  "Muck"] <- 2
UF_D_Muck_species[UF_D_Muck_species ==  "Not"] <- 1
UF_D_Muck_species[UF_D_Muck_species ==  "Mucky"] <- 2
UF_D_Muck_species[UF_D_Muck_species ==  "Muckish"] <- 2
UF_D_Muck_species<-as.data.frame(UF_D_Muck_species)


#Create Bray.Curtis distance matrix from species table#
bc <- vegdist(t(otu_table(UFDmuckSpeciesAspecies)))

#Conduct PAM clustering based upon 2 clusters#
pam.res <- pam(bc,2, diss=TRUE)

#Add a microbially.predicted column to matrix above for comparison to metadata.defined column#
UF_D_Muck_species <- cbind(UF_D_Muck_species, MP_UFDMuckspecies = pam.res$cluster)

#Test Total Muck Efficiency# 0.8555556
sum(UF_D_Muck_species$Muck == UF_D_Muck_species$MP_UFDMuckspecies)/nrow(UF_D_Muck_species)

#Test 3MC Efficiency aka Sensitivity# 0.7666667
UF_D_3MC_species <- filter(UF_D_Muck_species, Muck == "2")
sum(UF_D_3MC_species$Muck == UF_D_3MC_species$MP_UFDMuckspecies)/nrow(UF_D_3MC_species)

#Test 0MC Efficiency aka Sensitivity# 0.9
UF_D_0MC_species <- filter(UF_D_Muck_species, Muck == "1")
sum(UF_D_0MC_species$Muck == UF_D_0MC_species$MP_UFDMuckspecies)/nrow(UF_D_0MC_species)

TotalPCoA <- ordinate(UFDmuckSpeciesAspecies,"PCoA")
p = plot_ordination(UFDmuckSpeciesAspecies, TotalPCoA, color="Muck") + 
  ggtitle("Unfiltered DESeq2 Muck species Indicators PCoA with Bray.Curtis distance") +
  geom_text(aes(label=UF_D_Muck_species$MP_UFDMuckspecies),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p


##Test filtered Muck speciess##

#Remove any samples with no speciess from target phyloseq#
nsamples(FDmuckSpeciesAspecies) #90
FDmuckSpeciesAspecies = prune_samples(sample_sums(FDmuckSpeciesAspecies)>=1, FDmuckSpeciesAspecies)
nsamples(FDmuckSpeciesAspecies) #89

#Make a matrix with the sample_data and relabel subcategories to make metadata.defined column#
F_D_Muck_species <- as.matrix(sample_data(FDmuckSpeciesAspecies))
F_D_Muck_species[F_D_Muck_species ==  "Muck"] <- 1
F_D_Muck_species[F_D_Muck_species ==  "Not"] <- 2
F_D_Muck_species[F_D_Muck_species ==  "Mucky"] <- 1
F_D_Muck_species[F_D_Muck_species ==  "Muckish"] <- 1
F_D_Muck_species<-as.data.frame(F_D_Muck_species)


#Create Bray.Curtis distance matrix from species table#
bc <- vegdist(t(otu_table(FDmuckSpeciesAspecies)))

#Conduct PAM clustering based upon 2 clusters#
pam.res <- pam(bc,2, diss=TRUE)

#Add a microbially.predicted column to matrix above for comparison to metadata.defined column#
F_D_Muck_species <- cbind(F_D_Muck_species, MP_FDMuckspecies = pam.res$cluster)

#Test Total Muck Efficiency# 0.8777778
sum(F_D_Muck_species$Muck == F_D_Muck_species$MP_FDMuckspecies)/nrow(F_D_Muck_species)

#Test 3MC Efficiency aka Sensitivy# 0.8666667
F_D_3MC_species <- filter(F_D_Muck_species, Muck == "1")
sum(F_D_3MC_species$Muck == F_D_3MC_species$MP_FDMuckspecies)/nrow(F_D_3MC_species)

#Test 0MC Efficiency aka Specificity# 0.8833333
F_D_0MC_species <- filter(F_D_Muck_species, Muck == "2")
sum(F_D_0MC_species$Muck == F_D_0MC_species$MP_FDMuckspecies)/nrow(F_D_0MC_species)

TotalPCoA <- ordinate(FDmuckSpeciesAspecies,"PCoA")
p = plot_ordination(FDmuckSpeciesAspecies, TotalPCoA, color="Muck") + 
  ggtitle("Filtered DESeq2 Muck species Indicators PCoA with Bray.Curtis distance") +
  geom_text(aes(label=F_D_Muck_species$MP_FDMuckspecies),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p

###TEST COMBO MUCK INDICATORS###

##Test unfiltered Muck speciess##

#Remove any samples with no speciess from target phyloseq#
nsamples(UFCmuckSpeciesAspecies) #90
UFCmuckSpeciesAspecies = prune_samples(sample_sums(UFCmuckSpeciesAspecies)>=1, UFCmuckSpeciesAspecies)
nsamples(UFCmuckSpeciesAspecies) #90

#Make a matrix with the sample_data and relabel subcategories to make metadata.defined column#
UF_C_Muck_speciess <- as.matrix(sample_data(UFCmuckSpeciesAspecies))
UF_C_Muck_speciess[UF_C_Muck_speciess ==  "Muck"] <- 2
UF_C_Muck_speciess[UF_C_Muck_speciess ==  "Not"] <- 1
UF_C_Muck_speciess[UF_C_Muck_speciess ==  "Mucky"] <- 2
UF_C_Muck_speciess[UF_C_Muck_speciess ==  "Muckish"] <- 2
UF_C_Muck_speciess<-as.data.frame(UF_C_Muck_speciess)


#Create Bray.Curtis distance matrix from species table#
bc <- vegdist(t(otu_table(UFCmuckSpeciesAspecies)))

#Conduct PAM clustering based upon 2 clusters#
pam.res <- pam(bc,2, diss=TRUE)

#Add a microbially.predicted column to matrix above for comparison to metadata.defined column#
UF_C_Muck_speciess <- cbind(UF_C_Muck_speciess, MP_UFCMuckspecies = pam.res$cluster)

#Test Total Muck Efficiency# 0.8888889
sum(UF_C_Muck_speciess$Muck == UF_C_Muck_speciess$MP_UFCMuckspecies)/nrow(UF_C_Muck_speciess)

#Test Muck Efficiency aka Sensitivity# 0.9333333
F_C_3MC_speciess <- filter(UF_C_Muck_speciess, Muck == "2")
sum(F_C_3MC_speciess$Muck == F_C_3MC_speciess$MP_UFCMuckspecies)/nrow(F_C_3MC_speciess)

#Test Muck Efficiency aka Specificity# 0.8666667
F_C_0MC_speciess <- filter(UF_C_Muck_speciess, Muck == "1")
sum(F_C_0MC_speciess$Muck == F_C_0MC_speciess$MP_UFCMuckspecies)/nrow(F_C_0MC_speciess)

TotalPCoA <- ordinate(UFCmuckSpeciesAspecies,"PCoA")
p = plot_ordination(UFCmuckSpeciesAspecies, TotalPCoA, color="Muck") + 
  ggtitle("Unfiltered Final Muck species Indicators PCoA with Bray.Curtis distance") +
  geom_text(aes(label=UF_C_Muck_speciess$MP_UFCMuckspecies),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p

##Test filtered Muck speciess##

#Remove any samples with no speciess from target phyloseq#
nsamples(FCmuckSpeciesAspecies) #90
FCmuckSpeciesAspecies = prune_samples(sample_sums(FCmuckSpeciesAspecies)>=1, FCmuckSpeciesAspecies)
nsamples(FCmuckSpeciesAspecies) #88

#Make a matrix with the sample_data and relabel subcategories to make metadata.defined column#
F_C_Muck_speciess <- as.matrix(sample_data(FCmuckSpeciesAspecies))
F_C_Muck_speciess[F_C_Muck_speciess ==  "Muck"] <- 2
F_C_Muck_speciess[F_C_Muck_speciess ==  "Not"] <- 1
F_C_Muck_speciess[F_C_Muck_speciess ==  "Mucky"] <- 2
F_C_Muck_speciess[F_C_Muck_speciess ==  "Muckish"] <- 2
F_C_Muck_speciess<-as.data.frame(F_C_Muck_speciess)


#Create Bray.Curtis distance matrix from species table#
bc <- vegdist(t(otu_table(FCmuckSpeciesAspecies)))

#Conduct PAM clustering based upon 2 clusters#
pam.res <- pam(bc,2, diss=TRUE)

#Add a microbially.predicted column to matrix above for comparison to metadata.defined column#
F_C_Muck_speciess <- cbind(F_C_Muck_speciess, MP_FCMuckspecies = pam.res$cluster)

#Test Total Muck Efficiency# 0.8222222
sum(F_C_Muck_speciess$Muck == F_C_Muck_speciess$MP_FCMuckspecies)/nrow(F_C_Muck_speciess)

#Test Muck Efficiency aka Sensitivity# 0.7333333
F_C_3MC_species <- filter(F_C_Muck_speciess, Muck == "2")
sum(F_C_3MC_species$Muck == F_C_3MC_species$MP_FCMuckspecies)/nrow(F_C_3MC_species)

#Test Not Efficiency aka Specificity# 0.7333333
F_C_0MC_species <- filter(F_C_Muck_speciess, Muck == "1")
sum(F_C_0MC_species$Muck == F_C_0MC_species$MP_FCMuckspecies)/nrow(F_C_0MC_species)

TotalPCoA <- ordinate(FCmuckSpeciesAspecies,"PCoA")
p = plot_ordination(FCmuckSpeciesAspecies, TotalPCoA, color="Muck") + 
  ggtitle("Filtered Final Muck species Indicators PCoA with Bray.Curtis distance") +
  geom_text(aes(label=F_C_Muck_speciess$MP_FCMuckspecies),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p

###TEST INDICSPECIES TOM/Cu INDICATORS###

##Test unfiltered TOM/Cu speciess##

#Create phyloseq object focused on HiHi and HiLo samples#
UFIShihiSpeciesAspecies = subset_samples(UFIStomcuSpeciesAspecies, LOI.Cu=="High.High")
UFIShiloSpeciesAspecies = subset_samples(UFIStomcuSpeciesAspecies, LOI.Cu=="High.Low")
UFISfoctomcuSpeciesAspecies = merge_phyloseq(UFIShihiSpeciesAspecies, UFIShiloSpeciesAspecies)

#Remove any samples with no speciess from target phyloseq#
nsamples(UFISfoctomcuSpeciesAspecies) #29
UFISfoctomcuSpeciesAspecies = prune_samples(sample_sums(UFISfoctomcuSpeciesAspecies)>=1, UFISfoctomcuSpeciesAspecies)
nsamples(UFISfoctomcuSpeciesAspecies) #29

#Make a matrix with the sample_data and relabel subcategories to make metadata.defined column#
UF_IS_TOMCu_species <- as.matrix(sample_data(UFISfoctomcuSpeciesAspecies))
UF_IS_TOMCu_species[UF_IS_TOMCu_species ==  "High.High"] <- 1
UF_IS_TOMCu_species[UF_IS_TOMCu_species ==  "High.Low"] <- 2
UF_IS_TOMCu_species<-as.data.frame(UF_IS_TOMCu_species)

#Create Bray.Curtis distance matrix from species table#
bc <- vegdist(t(otu_table(UFISfoctomcuSpeciesAspecies)))

#Conduct PAM clustering based upon 2 clusters#
pam.res <- pam(bc,2, diss=TRUE)

#Add a microbially.predicted column to matrix above for comparison to metadata.defined column#
UF_IS_TOMCu_species <- cbind(UF_IS_TOMCu_species, MP_UFISTOMCuspecies = pam.res$cluster)

#Test Total TOM/Cu Efficiency# 0.8275862
sum(UF_IS_TOMCu_species$LOI.Cu == UF_IS_TOMCu_species$MP_UFISTOMCuspecies)/nrow(UF_IS_TOMCu_species)

#Test High.High Indicators aka Sensitivity# 1
UF_IS_HiHi_species <- filter(UF_IS_TOMCu_species, LOI.Cu == "1")
sum(UF_IS_HiHi_species$LOI.Cu == UF_IS_HiHi_species$MP_UFISTOMCuspecies)/nrow(UF_IS_HiHi_species)

#Test High.Low Indicators aka Specificity# 0.7222222
UF_IS_HiLo_species <- filter(UF_IS_TOMCu_species, LOI.Cu == "2")
sum(UF_IS_HiLo_species$LOI.Cu == UF_IS_HiLo_species$MP_UFISTOMCuspecies)/nrow(UF_IS_HiLo_species)

TotalPCoA <- ordinate(UFISfoctomcuSpeciesAspecies,"PCoA")
p = plot_ordination(UFISfoctomcuSpeciesAspecies, TotalPCoA, color="LOI.Cu") + 
  ggtitle("Unfiltered IS LOICu species Indicators PCoA with Bray.Curtis distance") +
  geom_text(aes(label=UF_IS_TOMCu_species$MP_UFISTOMCuspecies),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p


##Test filtered TOM/Cu speciess##

#Create phyloseq object focused on HiHi and HiLo samples#
FIShihiSpeciesAspecies = subset_samples(FIStomcuSpeciesAspecies, LOI.Cu=="High.High")
FIShiloSpeciesAspecies = subset_samples(FIStomcuSpeciesAspecies, LOI.Cu=="High.Low")
FISfoctomcuSpeciesAspecies = merge_phyloseq(FIShihiSpeciesAspecies, FIShiloSpeciesAspecies)

#Remove any samples with no speciess from target phyloseq#
nsamples(FISfoctomcuSpeciesAspecies) #29
FISfoctomcuSpeciesAspecies = prune_samples(sample_sums(FISfoctomcuSpeciesAspecies)>=1, FISfoctomcuSpeciesAspecies)
nsamples(FISfoctomcuSpeciesAspecies) #29

#Make a matrix with the sample_data and relabel subcategories to make metadata.defined column#
F_IS_TOMCu_species <- as.matrix(sample_data(FISfoctomcuSpeciesAspecies))
F_IS_TOMCu_species[F_IS_TOMCu_species ==  "High.High"] <- 1
F_IS_TOMCu_species[F_IS_TOMCu_species ==  "High.Low"] <- 2
F_IS_TOMCu_species<-as.data.frame(F_IS_TOMCu_species)


#Create Bray.Curtis distance matrix from species table#
bc <- vegdist(t(otu_table(FISfoctomcuSpeciesAspecies)))

#Conduct PAM clustering based upon 2 clusters#
pam.res <- pam(bc,2, diss=TRUE)

#Add a microbially.predicted column to matrix above for comparison to metadata.defined column#
F_IS_TOMCu_species <- cbind(F_IS_TOMCu_species, MP_FISTOMCuspecies = pam.res$cluster)

#Test Total TOM/Cu Efficiency# 1
sum(F_IS_TOMCu_species$LOI.Cu == F_IS_TOMCu_species$MP_FISTOMCuspecies)/nrow(F_IS_TOMCu_species)

#Test HiHi Efficiency aka Sensitivity# 1
F_IS_HiHi_species <- filter(F_IS_TOMCu_species, LOI.Cu == "1")
sum(F_IS_HiHi_species$LOI.Cu == F_IS_HiHi_species$MP_FISTOMCuspecies)/nrow(F_IS_HiHi_species)

#Test HiLo Efficiency aka Sensitivity# 0.3888889
F_IS_HiLo_species <- filter(F_IS_TOMCu_species, LOI.Cu == "2") 
sum(F_IS_HiLo_species$LOI.Cu == F_IS_HiLo_species$MP_FISTOMCuspecies)/nrow(F_IS_HiLo_species)

TotalPCoA <- ordinate(FISfoctomcuSpeciesAspecies,"PCoA")
p = plot_ordination(FISfoctomcuSpeciesAspecies, TotalPCoA, color="LOI.Cu") + 
  ggtitle("Filtered IS LOICu species Indicators PCoA with Bray.Curtis distance") +
  geom_text(aes(label=F_IS_TOMCu_species$MP_FISTOMCuspecies),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p

##TEST DESEQ2 LOICu INDICATORS##

##Test unfiltered TOM/Cu speciess##

#Create phyloseq object focused on HiHi and HiLo samples#
UFDhihiSpeciesAspecies = subset_samples(UFDtomcuSpeciesAspecies, LOI.Cu=="High.High")
UFDhiloSpeciesAspecies = subset_samples(UFDtomcuSpeciesAspecies, LOI.Cu=="High.Low")
UFDfoctomcuSpeciesAspecies = merge_phyloseq(UFDhihiSpeciesAspecies, UFDhiloSpeciesAspecies)

#Remove any samples with no speciess from target phyloseq#
nsamples(UFDfoctomcuSpeciesAspecies) #29
UFDfoctomcuSpeciesAspecies = prune_samples(sample_sums(UFDfoctomcuSpeciesAspecies)>=1, UFDfoctomcuSpeciesAspecies)
nsamples(UFDfoctomcuSpeciesAspecies) #29

#Make a matrix with the sample_data and relabel subcategories to make metadata.defined column#
UF_D_TOMCu_species <- as.matrix(sample_data(UFDfoctomcuSpeciesAspecies))
UF_D_TOMCu_species[UF_D_TOMCu_species ==  "High.High"] <- 1
UF_D_TOMCu_species[UF_D_TOMCu_species ==  "High.Low"] <- 2
UF_D_TOMCu_species<-as.data.frame(UF_D_TOMCu_species)


#Create Bray.Curtis distance matrix from species table#
bc <- vegdist(t(otu_table(UFDfoctomcuSpeciesAspecies)))

#Conduct PAM clustering based upon 2 clusters#
pam.res <- pam(bc,2, diss=TRUE)

#Add a microbially.predicted column to matrix above for comparison to metadata.defined column#
UF_D_TOMCu_species <- cbind(UF_D_TOMCu_species, MP_UFDTOMCuspecies = pam.res$cluster)

#Test Total TOM/Cu Efficiency# 0.8275862
sum(UF_D_TOMCu_species$LOI.Cu == UF_D_TOMCu_species$MP_UFDTOMCuspecies)/nrow(UF_D_TOMCu_species)

#Test HiHi Efficiency aka Sensitivity# 1
UF_D_HiHi_species <- filter(UF_D_TOMCu_species, LOI.Cu == "1")
sum(UF_D_HiHi_species$LOI.Cu == UF_D_HiHi_species$MP_UFDTOMCuspecies)/nrow(UF_D_HiHi_species)

#Test HiLo Efficiency aka Specificity# 0.7222222
UF_D_HiLo_species <- filter(UF_D_TOMCu_species, LOI.Cu == "2")
sum(UF_D_HiLo_species$LOI.Cu == UF_D_HiLo_species$MP_UFDTOMCuspecies)/nrow(UF_D_HiLo_species)

TotalPCoA <- ordinate(UFDfoctomcuSpeciesAspecies,"PCoA")
p = plot_ordination(UFDfoctomcuSpeciesAspecies, TotalPCoA, color="LOI.Cu") + 
  ggtitle("Unfiltered DESeq2 LOICu species Indicators PCoA with Bray.Curtis distance") +
  geom_text(aes(label=UF_D_TOMCu_species$MP_UFDTOMCuspecies),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p


##Test filtered TOM/Cu speciess##

#Create phyloseq object focused on HiHi and HiLo samples#
FDhihiSpeciesAspecies = subset_samples(FDtomcuSpeciesAspecies, LOI.Cu=="High.High")
FDhiloSpeciesAspecies = subset_samples(FDtomcuSpeciesAspecies, LOI.Cu=="High.Low")
FDfoctomcuSpeciesAspecies = merge_phyloseq(FDhihiSpeciesAspecies, FDhiloSpeciesAspecies)

#Remove any samples with no speciess from target phyloseq#
nsamples(FDfoctomcuSpeciesAspecies) #29
FDfoctomcuSpeciesAspecies = prune_samples(sample_sums(FDfoctomcuSpeciesAspecies)>=1, FDfoctomcuSpeciesAspecies)
nsamples(FDfoctomcuSpeciesAspecies) #29


#Make a matrix with the sample_data and relabel subcategories to make metadata.defined column#
F_D_TOMCu_species <- as.matrix(sample_data(FDfoctomcuSpeciesAspecies))
F_D_TOMCu_species[F_D_TOMCu_species ==  "High.High"] <- 1
F_D_TOMCu_species[F_D_TOMCu_species ==  "High.Low"] <- 2
F_D_TOMCu_species<-as.data.frame(F_D_TOMCu_species)


#Create Bray.Curtis distance matrix from species table#
bc <- vegdist(t(otu_table(FDfoctomcuSpeciesAspecies)))

#Conduct PAM clustering based upon 2 clusters#
pam.res <- pam(bc,2, diss=TRUE)

#Add a microbially.predicted column to matrix above for comparison to metadata.defined column#
F_D_TOMCu_species <- cbind(F_D_TOMCu_species, MP_FDTOMCuspecies = pam.res$cluster)

#Test Total TOM/Cu Efficiency# 0.5862069
sum(F_D_TOMCu_species$LOI.Cu == F_D_TOMCu_species$MP_FDTOMCuspecies)/nrow(F_D_TOMCu_species)

#Test HiHi Efficiency aka Sensitivity# 1
F_D_HiHi_species <- filter(F_D_TOMCu_species, LOI.Cu == "1")
sum(F_D_HiHi_species$LOI.Cu == F_D_HiHi_species$MP_FDTOMCuspecies)/nrow(F_D_HiHi_species)

#Test HiLo Efficiency aka Specificity# 0.3333333
F_D_HiLo_species <- filter(F_D_TOMCu_species, LOI.Cu == "2")
sum(F_D_HiLo_species$LOI.Cu == F_D_HiLo_species$MP_FDTOMCuspecies)/nrow(F_D_HiLo_species)

TotalPCoA <- ordinate(FDfoctomcuSpeciesAspecies,"PCoA")
p = plot_ordination(FDfoctomcuSpeciesAspecies, TotalPCoA, color="LOI.Cu") + 
  ggtitle("Filtered DESeq2 LOICu species Indicators PCoA with Bray.Curtis distance") +
  geom_text(aes(label=F_D_TOMCu_species$MP_FDTOMCuspecies),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p

##TEST COMBO TOM/CU INDICATORS##

##Test unfiltered TOM/Cu speciess##

#Create phyloseq object focused on HiHi and HiLo samples#
UFChihiSpeciesAspecies = subset_samples(UFCtomcuSpeciesAspecies, LOI.Cu=="High.High")
UFChiloSpeciesAspecies = subset_samples(UFCtomcuSpeciesAspecies, LOI.Cu=="High.Low")
UFCfoctomcuSpeciesAspecies = merge_phyloseq(UFChihiSpeciesAspecies, UFChiloSpeciesAspecies)

#Remove any samples with no speciess from target phyloseq#
nsamples(UFCfoctomcuSpeciesAspecies) #29
UFCfoctomcuSpeciesAspecies = prune_samples(sample_sums(UFCfoctomcuSpeciesAspecies)>=1, UFCfoctomcuSpeciesAspecies)
nsamples(UFCfoctomcuSpeciesAspecies) #29

#Make a matrix with the sample_data and relabel subcategories to make metadata.defined column#
UF_C_TOMCu_species <- as.matrix(sample_data(UFCfoctomcuSpeciesAspecies))
UF_C_TOMCu_species[UF_C_TOMCu_species ==  "High.High"] <- 1
UF_C_TOMCu_species[UF_C_TOMCu_species ==  "High.Low"] <- 2
UF_C_TOMCu_species<-as.data.frame(UF_C_TOMCu_species)


#Create Bray.Curtis distance matrix from species table#
bc <- vegdist(t(otu_table(UFCfoctomcuSpeciesAspecies)))

#Conduct PAM clustering based upon 2 clusters#
pam.res <- pam(bc,2, diss=TRUE)

#Add a microbially.predicted column to matrix above for comparison to metadata.defined column#
UF_C_TOMCu_species <- cbind(UF_C_TOMCu_species, MP_UFCTOMCuspecies = pam.res$cluster)

#Test Total TOM/Cu Efficiency# 0.8275862
sum(UF_C_TOMCu_species$LOI.Cu == UF_C_TOMCu_species$MP_UFCTOMCuspecies)/nrow(UF_C_TOMCu_species)

#Test HiHi Efficiency aka Sensitivity# 1
UF_C_HiHi_species <- filter(UF_C_TOMCu_species, LOI.Cu == "1")
sum(UF_C_HiHi_species$LOI.Cu == UF_C_HiHi_species$MP_UFCTOMCuspecies)/nrow(UF_C_HiHi_species)

#Test HiLo Efficiency aka Specificity# 0.7222222
UF_C_HiLo_species <- filter(UF_C_TOMCu_species, LOI.Cu == "2")
sum(UF_C_HiLo_species$LOI.Cu == UF_C_HiLo_species$MP_UFCTOMCuspecies)/nrow(UF_C_HiLo_species)

TotalPCoA <- ordinate(UFCfoctomcuSpeciesAspecies,"PCoA")
p = plot_ordination(UFCfoctomcuSpeciesAspecies, TotalPCoA, color="LOI.Cu") + 
  ggtitle("Unfiltered Final LOICu species Indicators PCoA with Bray.Curtis distance") +
  geom_text(aes(label=UF_C_TOMCu_species$MP_UFCTOMCuspecies),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p


##Test filtered TOM/Cu speciess##

#Create phyloseq object focused on HiHi and HiLo samples#
FChihiSpeciesAspecies = subset_samples(FCtomcuSpeciesAspecies, LOI.Cu=="High.High")
FChiloSpeciesAspecies = subset_samples(FCtomcuSpeciesAspecies, LOI.Cu=="High.Low")
FCfoctomcuSpeciesAspecies = merge_phyloseq(FChihiSpeciesAspecies, FChiloSpeciesAspecies)

#Remove any samples with no speciess from target phyloseq#
nsamples(FCfoctomcuSpeciesAspecies) #29
FCfoctomcuSpeciesAspecies = prune_samples(sample_sums(FCfoctomcuSpeciesAspecies)>=1, FCfoctomcuSpeciesAspecies)
nsamples(FCfoctomcuSpeciesAspecies) #29

#Make a matrix with the sample_data and relabel subcategories to make metadata.defined column#
F_C_TOMCu_species <- as.matrix(sample_data(FCfoctomcuSpeciesAspecies))
F_C_TOMCu_species[F_C_TOMCu_species ==  "High.High"] <- 1
F_C_TOMCu_species[F_C_TOMCu_species ==  "High.Low"] <- 2
F_C_TOMCu_species<-as.data.frame(F_C_TOMCu_species)


#Create Bray.Curtis distance matrix from species table#
bc <- vegdist(t(otu_table(FCfoctomcuSpeciesAspecies)))

#Conduct PAM clustering based upon 2 clusters#
pam.res <- pam(bc,2, diss=TRUE)

#Add a microbially.predicted column to matrix above for comparison to metadata.defined column#
F_C_TOMCu_species <- cbind(F_C_TOMCu_species, MP_FCTOMCuspecies = pam.res$cluster)

#Test Total TOM/Cu Efficiency# 0.9655172
sum(F_C_TOMCu_species$LOI.Cu == F_C_TOMCu_species$MP_FCTOMCuspecies)/nrow(F_C_TOMCu_species)

#Test HiHi Efficiency aka Sensitivity# 0.9090909
F_C_HiHi_species <- filter(F_C_TOMCu_species, LOI.Cu == "1")
sum(F_C_HiHi_species$LOI.Cu == F_C_HiHi_species$MP_FCTOMCuspecies)/nrow(F_C_HiHi_species)

#Test HiLo Efficiency aka Specificity# 0.9032258
F_C_HiLo_species <- filter(F_C_TOMCu_species, LOI.Cu == "2")
sum(F_C_HiLo_species$LOI.Cu == F_C_HiLo_species$MP_FCTOMCuspecies)/nrow(F_C_HiLo_species)


TotalPCoA <- ordinate(FCfoctomcuSpeciesAspecies,"PCoA")
p = plot_ordination(FCfoctomcuSpeciesAspecies, TotalPCoA, color="LOI.Cu") + 
  ggtitle("Filtered Final LOICu species Indicators PCoA with Bray.Curtis distance") +
  geom_text(aes(label=F_C_TOMCu_species$MP_FCTOMCuspecies),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p


###GENUS###

getwd() #"C:/Users/bradshawd/Documents/Journal_Articles/New_Indicators"
setwd("C:/Users/bradshawd/Documents/Journal_Articles/New_Indicators/Genus")
getwd() #"C:/Users/bradshawd/Documents/Journal_Articles/New_Indicators/Genus"

##Upload QIIME2 information into phyloseq for taxnomic level you want to test##
genus_aa_table = read.csv(file.choose(), row.names=1)
genus_taxa_table = read.csv(file.choose(), row.names = 1)
genus_taxa_table = as.matrix(genus_taxa_table)

#Can use the same metadata table as species level#
Genus_OTU = otu_table(genus_aa_table, taxa_are_rows = TRUE)
Genus_TAX = tax_table(genus_taxa_table)

GenusA <- merge_phyloseq(Genus_OTU, Genus_TAX, dotted_meta)

colnames(tax_table(GenusA))

###Create at more filtered set of LWSS Data###
ntaxa(GenusA) #2883#
nsamples(GenusA)#483

LWSGenusA <- subset_samples(GenusA, Survey=="Lagoon Wide")
nsamples(LWSGenusA)#324

LWSWGenusA <- subset_samples(LWSGenusA, Medium=="Water")
nsamples(LWSWGenusA)#132, 11 sites in triplicate (33) x 4 sampling periods (sps)

LWSSGenusA1 <- subset_samples(GenusA, Survey=="Lagoon Wide.Post Hurricane")
nsamples(LWSSGenusA1)#12

LWSSGenusA2 <- subset_samples(LWSGenusA, Medium=="Sediment")
nsamples(LWSSGenusA2)#192

LWSSGenusA <- merge_phyloseq(LWSSGenusA1, LWSSGenusA2)
nsamples(LWSSGenusA)#204

##Create phyloseqs for determining and testing indicators for Sediment
#Split samples by sampling period (SP)
LWSS_W16_GenusA <- subset_samples(LWSSGenusA, AbbSeasonAbbYear=="W16")
nsamples(LWSS_W16_GenusA)#45, 15 sites in triplicate

LWSS_D17_GenusA <- subset_samples(LWSSGenusA, AbbSeasonAbbYear=="D17")
nsamples(LWSS_D17_GenusA)#45,15 sites in triplicate

LWSS_W17_GenusA <- subset_samples(LWSSGenusA, AbbSeasonAbbYear=="W17")
nsamples(LWSS_W17_GenusA)#57, 19 sites in triplicate

LWSS_D18_GenusA <- subset_samples(LWSSGenusA, AbbSeasonAbbYear=="D18")
nsamples(LWSS_D18_GenusA)#57, 19 sites in triplicate

#Combine sampling periods into first and second years
LWSS_Y1_GenusA <- merge_phyloseq(LWSS_W16_GenusA, LWSS_D17_GenusA)
nsamples(LWSS_Y1_GenusA)#90

LWSS_Y2_GenusA <- merge_phyloseq(LWSS_W17_GenusA, LWSS_D18_GenusA)
nsamples(LWSS_Y2_GenusA)#114

##Create phyloseqs for determining and testing indicators for Water
#Split samples by sampling period (SP)
LWSW_W16_GenusA <- subset_samples(LWSWGenusA, AbbSeasonAbbYear=="W16")
nsamples(LWSW_W16_GenusA)#33, 11 sites in triplicate

LWSW_D17_GenusA <- subset_samples(LWSWGenusA, AbbSeasonAbbYear=="D17")
nsamples(LWSW_D17_GenusA)#33,11 sites in triplicate

LWSW_W17_GenusA <- subset_samples(LWSWGenusA, AbbSeasonAbbYear=="W17")
nsamples(LWSW_W17_GenusA)#33, 11 sites in triplicate

LWSW_D18_GenusA <- subset_samples(LWSWGenusA, AbbSeasonAbbYear=="D18")
nsamples(LWSW_D18_GenusA)#33, 11 sites in triplicate

#Combine sampling periods into first and second years
LWSW_Y1_GenusA <- merge_phyloseq(LWSW_W16_GenusA, LWSW_D17_GenusA)
nsamples(LWSW_Y1_GenusA)#66

LWSW_Y2_GenusA <- merge_phyloseq(LWSW_W17_GenusA, LWSW_D18_GenusA)
nsamples(LWSW_Y2_GenusA)#66

##Filter yearly sample sets of all zero and then low abundance genuss
FLWSS_Y1_GenusA <- filter_taxa(LWSS_Y1_GenusA, function(x) sum(x) >0, TRUE)
ntaxa(FLWSS_Y1_GenusA) #2087#
RFLWSS_Y1_GenusAgenus <- filter_taxa(FLWSS_Y1_GenusA, function(x) sum(x) >20, TRUE)
ntaxa(RFLWSS_Y1_GenusAgenus) #1233#

FLWSS_Y2_GenusA <- filter_taxa(LWSS_Y2_GenusA, function(x) sum(x) >0, TRUE)
ntaxa(FLWSS_Y2_GenusA) #2201#
RFLWSS_Y2_GenusAgenus <- filter_taxa(FLWSS_Y2_GenusA, function(x) sum(x) >20, TRUE)
ntaxa(RFLWSS_Y2_GenusAgenus) #1098#

FLWSW_Y1_GenusA <- filter_taxa(LWSW_Y1_GenusA, function(x) sum(x) >0, TRUE)
ntaxa(FLWSW_Y1_GenusA) #1143#
RFLWSW_Y1_GenusAgenus <- filter_taxa(FLWSW_Y1_GenusA, function(x) sum(x) >20, TRUE)
ntaxa(RFLWSW_Y1_GenusAgenus) #566#

FLWSW_Y2_GenusA <- filter_taxa(LWSW_Y2_GenusA, function(x) sum(x) >0, TRUE)
ntaxa(FLWSW_Y2_GenusA) #1555#
RFLWSW_Y2_GenusAgenus <- filter_taxa(FLWSW_Y2_GenusA, function(x) sum(x) >20, TRUE)
ntaxa(RFLWSW_Y2_GenusAgenus) #690#





#Determine the number of samples in each category and subcategory#
length(which(sample_data(RFLWSS_Y1_GenusAgenus)$Estuary == "SLE")) #24
length(which(sample_data(RFLWSS_Y1_GenusAgenus)$Estuary == "IRL")) #66

length(which(sample_data(RFLWSS_Y1_GenusAgenus)$Muck == "Not")) #60
length(which(sample_data(RFLWSS_Y1_GenusAgenus)$Muck == "Muck")) #22
length(which(sample_data(RFLWSS_Y1_GenusAgenus)$Muck == "Mucky")) #4
length(which(sample_data(RFLWSS_Y1_GenusAgenus)$Muck == "Muckish")) #4

length(which(sample_data(RFLWSS_Y1_GenusAgenus)$LOI.Cu == "Low.Low")) #61
length(which(sample_data(RFLWSS_Y1_GenusAgenus)$LOI.Cu == "High.Low")) #18
length(which(sample_data(RFLWSS_Y1_GenusAgenus)$LOI.Cu == "High.High")) #11
length(which(sample_data(RFLWSS_Y1_GenusAgenus)$LOI.Cu == "Low.High")) #0


###To run above scripts at any other level, just replace all lowercase "genus" with the level you are testing at to differentiate from genus level by (Ctrl + F all lowercase "asv" with "genus"). Also replace "DataA" with "taxolevelA". ESV.ID is programmed into the results of DESeq2, thus it is used throughout despite different taxonomic levels###

##Example##
#Convert phyloseq to deseq2 object centered around a irlvssle factor#
estuary_cudds_asv = phyloseq_to_deseq2(RFLWSS_Y1_DataAasv,  ~ Estuary)

#Becomes...#

#Convert phyloseq to deseq2 object centered around a irlvssle factor#
estuary_cudds_genus = phyloseq_to_deseq2(RLWSSGenusAgenus,  ~ Estuary)


##NOTE##
#Use res to determine what is + or .#
#muck res: Not (+) vs Muck (.)#
#tomcu res: HL (+) vs HH (.)#
#estuary res: SLE (+) vs IRL (.)#

###RUN DESEQ2 AT genus LEVEL ON ESTUARY SUBCATEGORIES###
#Convert phyloseq to deseq2 object centered around a irlvssle factor#
estuary_cudds_genus = phyloseq_to_deseq2(RFLWSS_Y1_GenusAgenus,  ~ Estuary)

#Calculate geometric means prior to estimate size factors#
estuary_geoMeans_genus = apply(counts(estuary_cudds_genus), 1, gm_mean)
estuary_cudds_genus = estimateSizeFactors(estuary_cudds_genus, geoMeans = estuary_geoMeans_genus)

#Conduct DESEQ2 test#
estuary_cudds_genus = DESeq(estuary_cudds_genus, fitType="local")

#Explore the results#
estuary_DESeq2_res_genus = results(estuary_cudds_genus)
estuary_DESeq2_res_genus = estuary_DESeq2_res_genus[order(estuary_DESeq2_res_genus$padj, na.last=NA), ]
alpha = 0.05
estuary_DESeq2_sig_res_genus = estuary_DESeq2_res_genus[(estuary_DESeq2_res_genus$padj < alpha), ]

#Make dataframe with taxanomy added in#
estuary_DESeq2_sig_res_taxo_genus = cbind(as(estuary_DESeq2_sig_res_genus, "data.frame"), as(tax_table(RFLWSS_Y1_GenusAgenus)[rownames(estuary_DESeq2_sig_res_genus), ], "matrix"))

#Make dataframe with genuss added in from all sites#
estuary_DESeq2_sig_res_taxo_seqs_genus = cbind(as(estuary_DESeq2_sig_res_taxo_genus, "data.frame"), as(otu_table(RFLWSS_Y1_GenusAgenus)[rownames(estuary_DESeq2_sig_res_taxo_genus), ], "matrix"))

#Make rownames an actual column and remove old rownames#
estuary_DESeq2_sig_res_taxo_seqs_genus <- cbind(ESV.ID = rownames(estuary_DESeq2_sig_res_taxo_seqs_genus), estuary_DESeq2_sig_res_taxo_seqs_genus)
rownames(estuary_DESeq2_sig_res_taxo_seqs_genus) <-NULL

#Make a txt file of results for all sites#
write.csv(as.data.frame(estuary_DESeq2_sig_res_taxo_seqs_genus), file="estuary_DESeq2_sig_res_taxo_seqs_genus.csv")

#Detemine which subcategory is negative or positive#
estuary_DESeq2_res_genus
#estuary_DESeq2_res_genus: SLE (+) vs IRL (.)#

#Determine number of Indicators with padj >0.05 in each subcagetory#
#Number of SLE indicators# 
length(which(estuary_DESeq2_sig_res_taxo_seqs_genus$log2FoldChange > 0)) #210
#Number of IRL indicators#
length(which(estuary_DESeq2_sig_res_taxo_seqs_genus$log2FoldChange < 0)) #279

###RUN DESEQ2 AT genus LEVEL ON MUCK EXTREME SUBCATEGORIES###

##Prepare for running DESeq2##
#Make phyloseq focusing on Muck and Not subcategories#  
RFLWSS_Y1_Muck_GenusAgenus = subset_samples(RFLWSS_Y1_GenusAgenus, Muck=="Muck")
RFLWSS_Y1_Not_GenusAgenus = subset_samples(RFLWSS_Y1_GenusAgenus, Muck=="Not")
RFLWSS_Y1_MuckvsNot_GenusAgenus = merge_phyloseq(RFLWSS_Y1_Muck_GenusAgenus, RFLWSS_Y1_Not_GenusAgenus)

##Run DESEQ2 and make a txt table of results##
#Convert phyloseq to deseq2 object centered around a muck factor#
muck_cudds_genus = phyloseq_to_deseq2(RFLWSS_Y1_MuckvsNot_GenusAgenus,  ~ Muck)

#Calculate geometric means prior to estimate size factors#
muck_geoMeans_genus = apply(counts(muck_cudds_genus), 1, gm_mean)
muck_cudds_genus = estimateSizeFactors(muck_cudds_genus, geoMeans = muck_geoMeans_genus)

#Conduct DESEQ2 test#
muck_cudds_genus = DESeq(muck_cudds_genus, fitType="local")

#Explore the results#
muck_DESeq2_res_genus = results(muck_cudds_genus)
muck_DESeq2_res_genus = muck_DESeq2_res_genus[order(muck_DESeq2_res_genus$padj, na.last=NA), ]
alpha = 0.05
muck_DESeq2_sig_res_genus = muck_DESeq2_res_genus[(muck_DESeq2_res_genus$padj < alpha), ]

#Make dataframe with taxanomy added in#
muck_DESeq2_sig_res_taxo_genus = cbind(as(muck_DESeq2_sig_res_genus, "data.frame"), as(tax_table(RFLWSS_Y1_MuckvsNot_GenusAgenus)[rownames(muck_DESeq2_sig_res_genus), ], "matrix"))

#Make dataframe with genuss added in from all sites#
muck_DESeq2_sig_res_taxo_seqs_genus = cbind(as(muck_DESeq2_sig_res_taxo_genus, "data.frame"), as(otu_table(RFLWSS_Y1_GenusAgenus)[rownames(muck_DESeq2_sig_res_taxo_genus), ], "matrix"))

#Make rownames an actual column and remove old rownames#
muck_DESeq2_sig_res_taxo_seqs_genus <- cbind(ESV.ID = rownames(muck_DESeq2_sig_res_taxo_seqs_genus), muck_DESeq2_sig_res_taxo_seqs_genus)
rownames(muck_DESeq2_sig_res_taxo_seqs_genus) <-NULL

#Make a txt file of results for all sites#
write.csv(as.data.frame(muck_DESeq2_sig_res_taxo_seqs_genus), file="muck_DESeq2_sig_res_taxo_seqs_genus.csv")

#Detemine which subcategory is negative or positive#
muck_DESeq2_res_genus
#muck_DESeq2_res_genus: Not (+) vs Muck (.)#
#Determine number of Indicators with padj >0.05 in each subcagetory#
#Number of Not indicators# 
length(which(muck_DESeq2_sig_res_taxo_seqs_genus$log2FoldChange > 0)) #164
#Number of Muck indicators#
length(which(muck_DESeq2_sig_res_taxo_seqs_genus$log2FoldChange < 0)) #94

###RUN DESEQ2 AT genus LEVEL ON TOM/Cu EXTREME SUBCATEGORIES###

##Prepare for running DESeq2##
#Make phyloseq focusing on High.High and High.Low subcategories#  
RFLWSS_Y1_HH_GenusAgenus = subset_samples(RFLWSS_Y1_GenusAgenus, LOI.Cu=="High.High")
RFLWSS_Y1_HL_GenusAgenus = subset_samples(RFLWSS_Y1_GenusAgenus, LOI.Cu=="High.Low")
RFLWSS_Y1_HHvsHL_GenusAgenus = merge_phyloseq(RFLWSS_Y1_HH_GenusAgenus, RFLWSS_Y1_HL_GenusAgenus)

##Run DESEQ2 and make a txt table of results##
#Convert phyloseq to deseq2 object centered around a loicu factor#
tomcu_cudds_genus = phyloseq_to_deseq2(RFLWSS_Y1_HHvsHL_GenusAgenus,  ~ LOI.Cu)

#Calculate geometric means prior to estimate size factors#
tomcu_geoMeans_genus = apply(counts(tomcu_cudds_genus), 1, gm_mean)
tomcu_cudds_genus = estimateSizeFactors(tomcu_cudds_genus, geoMeans = tomcu_geoMeans_genus)

#Conduct DESEQ2 test#
tomcu_cudds_genus = DESeq(tomcu_cudds_genus, fitType="local")

#Explore the results#
tomcu_DESeq2_res_genus = results(tomcu_cudds_genus)
tomcu_DESeq2_res_genus = tomcu_DESeq2_res_genus[order(tomcu_DESeq2_res_genus$padj, na.last=NA), ]
alpha = 0.05
tomcu_DESeq2_sig_res_genus = tomcu_DESeq2_res_genus[(tomcu_DESeq2_res_genus$padj < alpha), ]

#Make dataframe with taxanomy added in#
tomcu_DESeq2_sig_res_taxo_genus = cbind(as(tomcu_DESeq2_sig_res_genus, "data.frame"), as(tax_table(RFLWSS_Y1_HHvsHL_GenusAgenus)[rownames(tomcu_DESeq2_sig_res_genus), ], "matrix"))

#Make dataframe with genuss added in from all sites#
tomcu_DESeq2_sig_res_taxo_seqs_genus = cbind(as(tomcu_DESeq2_sig_res_taxo_genus, "data.frame"), as(otu_table(RFLWSS_Y1_GenusAgenus)[rownames(tomcu_DESeq2_sig_res_taxo_genus), ], "matrix"))

#Make rownames an actual column and remove old rownames#
tomcu_DESeq2_sig_res_taxo_seqs_genus <- cbind(ESV.ID = rownames(tomcu_DESeq2_sig_res_taxo_seqs_genus), tomcu_DESeq2_sig_res_taxo_seqs_genus)
rownames(tomcu_DESeq2_sig_res_taxo_seqs_genus) <-NULL

#Make a txt file of results for all sites#
write.csv(as.data.frame(tomcu_DESeq2_sig_res_taxo_seqs_genus), file="tomcu_DESeq2_sig_res_taxo_seqs_genus.csv")

#Detemine which subcategory is negative or positive#
tomcu_DESeq2_res_genus
#tomcu_DESeq2_res_genus: High.Low (+) vs High.High (.)#
#Determine number of Indicators with padj >0.05 in each subcagetory#
#Number of High.Low indicators# 
length(which(tomcu_DESeq2_sig_res_taxo_seqs_genus$log2FoldChange > 0)) #137
#Number of High.High indicators#
length(which(tomcu_DESeq2_sig_res_taxo_seqs_genus$log2FoldChange < 0)) #80

#Order rows by log2FoldChange#
muck_DESeq2_sig_res_taxo_seqs_genus <- arrange(muck_DESeq2_sig_res_taxo_seqs_genus, log2FoldChange)
tomcu_DESeq2_sig_res_taxo_seqs_genus <- arrange(tomcu_DESeq2_sig_res_taxo_seqs_genus, log2FoldChange)
estuary_DESeq2_sig_res_taxo_seqs_genus <- arrange(estuary_DESeq2_sig_res_taxo_seqs_genus, log2FoldChange)

#Create lists of just significant genuss for each category#
unfiltered_DESeq2_estuary_genuss <- subset(estuary_DESeq2_sig_res_taxo_seqs_genus, select=c(ESV.ID))
unfiltered_DESeq2_muck_genuss <- subset(muck_DESeq2_sig_res_taxo_seqs_genus, select=c(ESV.ID))
unfiltered_DESeq2_tomcu_genuss <- subset(tomcu_DESeq2_sig_res_taxo_seqs_genus, select=c(ESV.ID))



##Create indicator tables filtered of overlapping indicators with other two metadata categories##

#Filter out rows in the IRLSLE DESEQ object that are also in the LOI.Cu and Muck DESEQs#
filtered_estuary_DESeq2_sig_res_taxo_seqs_genus <- anti_join(estuary_DESeq2_sig_res_taxo_seqs_genus, unfiltered_DESeq2_tomcu_genuss, by = "ESV.ID")
filtered_estuary_DESeq2_sig_res_taxo_seqs_genus <- anti_join(filtered_estuary_DESeq2_sig_res_taxo_seqs_genus, unfiltered_DESeq2_muck_genuss, by = "ESV.ID")

#Filter out rows in the Muck DESEQ object that are also in the IRL.SLE and LOI.Cu DESEQs#
filtered_muck_DESeq2_sig_res_taxo_seqs_genus <- anti_join(muck_DESeq2_sig_res_taxo_seqs_genus, unfiltered_DESeq2_estuary_genuss, by = "ESV.ID")
filtered_muck_DESeq2_sig_res_taxo_seqs_genus <- anti_join(filtered_muck_DESeq2_sig_res_taxo_seqs_genus, unfiltered_DESeq2_tomcu_genuss, by = "ESV.ID")

#Filter out rows in the LOI.CU DESEQ object that are also in the IRL.SLE and Muck DESEQs#
filtered_tomcu_DESeq2_sig_res_taxo_seqs_genus <- anti_join(tomcu_DESeq2_sig_res_taxo_seqs_genus, unfiltered_DESeq2_estuary_genuss, by = "ESV.ID")
filtered_tomcu_DESeq2_sig_res_taxo_seqs_genus <- anti_join(filtered_tomcu_DESeq2_sig_res_taxo_seqs_genus, unfiltered_DESeq2_muck_genuss, by = "ESV.ID")



#Order rows by log2FoldChange#
filtered_muck_DESeq2_sig_res_taxo_seqs_genus <- arrange(filtered_muck_DESeq2_sig_res_taxo_seqs_genus, log2FoldChange)
filtered_tomcu_DESeq2_sig_res_taxo_seqs_genus <- arrange(filtered_tomcu_DESeq2_sig_res_taxo_seqs_genus, log2FoldChange)
filtered_estuary_DESeq2_sig_res_taxo_seqs_genus <- arrange(filtered_estuary_DESeq2_sig_res_taxo_seqs_genus, log2FoldChange)

#Make a csv file for each of the filtered tables#
write.csv(as.data.frame(filtered_estuary_DESeq2_sig_res_taxo_seqs_genus), file="filtered_estuary_DESeq2_sig_res_taxo_seqs_genus.csv")
write.csv(as.data.frame(filtered_muck_DESeq2_sig_res_taxo_seqs_genus), file="filtered_muck_DESeq2_sig_res_taxo_seqs_genus.csv")
write.csv(as.data.frame(filtered_tomcu_DESeq2_sig_res_taxo_seqs_genus), file="filtered_tomcu_DESeq2_sig_res_taxo_seqs_genus.csv")



##Determine number of filtered indicators##

#Number of Not indicators# 
length(which(filtered_muck_DESeq2_sig_res_taxo_seqs_genus$log2FoldChange > 0)) #33
#Number of Muck indicators#
length(which(filtered_muck_DESeq2_sig_res_taxo_seqs_genus$log2FoldChange < 0)) #46
#Number of High.Low indicators# 
length(which(filtered_tomcu_DESeq2_sig_res_taxo_seqs_genus$log2FoldChange > 0)) #14
#Number of High.High indicators#
length(which(filtered_tomcu_DESeq2_sig_res_taxo_seqs_genus$log2FoldChange < 0)) #14
#Number of SLE indicators# 
length(which(filtered_estuary_DESeq2_sig_res_taxo_seqs_genus$log2FoldChange > 0)) #82
#Number of IRL indicators#
length(which(filtered_estuary_DESeq2_sig_res_taxo_seqs_genus$log2FoldChange < 0)) #108

#Create lists of filtered significant genuss for each category#
filtered_DESeq2_estuary_genuss <- subset(filtered_estuary_DESeq2_sig_res_taxo_seqs_genus, select=c(ESV.ID))
filtered_DESeq2_muck_genuss <- subset(filtered_muck_DESeq2_sig_res_taxo_seqs_genus, select=c(ESV.ID))
filtered_DESeq2_tomcu_genuss <- subset(filtered_tomcu_DESeq2_sig_res_taxo_seqs_genus, select=c(ESV.ID))

##Create phyloseq objects for DESeq2 indicators##

#Create phyloseq objects of the unfiltered indicators using character string version of significant genuss##
charac_unfiltered_DESeq2_estuary_genuss <- unfiltered_DESeq2_estuary_genuss[,"ESV.ID"]
charac_unfiltered_DESeq2_estuary_genuss <- as.character(charac_unfiltered_DESeq2_estuary_genuss)
UFDestuaryGenusAgenus <- prune_taxa(charac_unfiltered_DESeq2_estuary_genuss, RFLWSS_Y1_GenusAgenus)

charac_unfiltered_DESeq2_muck_genuss <- unfiltered_DESeq2_muck_genuss[,"ESV.ID"]
charac_unfiltered_DESeq2_muck_genuss <- as.character(charac_unfiltered_DESeq2_muck_genuss)
UFDmuckGenusAgenus <- prune_taxa(charac_unfiltered_DESeq2_muck_genuss, RFLWSS_Y1_GenusAgenus)

charac_unfiltered_DESeq2_tomcu_genuss <- unfiltered_DESeq2_tomcu_genuss[,"ESV.ID"]
charac_unfiltered_DESeq2_tomcu_genuss <- as.character(charac_unfiltered_DESeq2_tomcu_genuss)
UFDtomcuGenusAgenus <- prune_taxa(charac_unfiltered_DESeq2_tomcu_genuss, RFLWSS_Y1_GenusAgenus)

#Create phyloseq objects of the filtered indicators using character string version of significant genuss##
charac_filtered_DESeq2_estuary_genuss <- filtered_DESeq2_estuary_genuss[,"ESV.ID"]
charac_filtered_DESeq2_estuary_genuss <- as.character(charac_filtered_DESeq2_estuary_genuss)
FDestuaryGenusAgenus <- prune_taxa(charac_filtered_DESeq2_estuary_genuss, RFLWSS_Y1_GenusAgenus)

charac_filtered_DESeq2_muck_genuss <- filtered_DESeq2_muck_genuss[,"ESV.ID"]
charac_filtered_DESeq2_muck_genuss <- as.character(charac_filtered_DESeq2_muck_genuss)
FDmuckGenusAgenus <- prune_taxa(charac_filtered_DESeq2_muck_genuss, RFLWSS_Y1_GenusAgenus)

charac_filtered_DESeq2_tomcu_genuss <- filtered_DESeq2_tomcu_genuss[,"ESV.ID"]
charac_filtered_DESeq2_tomcu_genuss <- as.character(charac_filtered_DESeq2_tomcu_genuss)
FDtomcuGenusAgenus <- prune_taxa(charac_filtered_DESeq2_tomcu_genuss, RFLWSS_Y1_GenusAgenus)



##Create heatmaps to see overall trends, may take a while to load##

unfiltered_DESeq2_estuary_genus_heatmap <- plot_heatmap(UFDestuaryGenusAgenus, taxa.order = charac_unfiltered_DESeq2_estuary_genuss, sample.order= "Estuary", sample.label = "Estuary")
unfiltered_DESeq2_estuary_genus_heatmap

unfiltered_DESeq2_muck_genus_heatmap <- plot_heatmap(UFDmuckGenusAgenus, taxa.order = charac_unfiltered_DESeq2_muck_genuss, sample.order= "Muck", sample.label = "Muck")
unfiltered_DESeq2_muck_genus_heatmap

unfiltered_DESeq2_tomcu_genus_heatmap <- plot_heatmap(UFDtomcuGenusAgenus, taxa.order = charac_unfiltered_DESeq2_tomcu_genuss, sample.order= "LOI.Cu", sample.label = "LOI.Cu")
unfiltered_DESeq2_tomcu_genus_heatmap

filtered_DESeq2_estuary_genus_heatmap <- plot_heatmap(FDestuaryGenusAgenus, taxa.order = charac_filtered_DESeq2_estuary_genuss, sample.order= "Estuary", sample.label = "Estuary")
filtered_DESeq2_estuary_genus_heatmap

filtered_DESeq2_muck_genus_heatmap <- plot_heatmap(FDmuckGenusAgenus, taxa.order = charac_filtered_DESeq2_muck_genuss, sample.order= "Muck", sample.label = "Muck")
filtered_DESeq2_muck_genus_heatmap

filtered_DESeq2_tomcu_genus_heatmap <- plot_heatmap(FDtomcuGenusAgenus, taxa.order = charac_filtered_DESeq2_tomcu_genuss, sample.order= "LOI.Cu", sample.label = "LOI.Cu")
filtered_DESeq2_tomcu_genus_heatmap


###CONDUCT INDICSPECIES ANALYSIS AT genus LEVEL###

#Increase max print so that all of results can be printed#
options(max.print=100000000)


##Conduct genus indicspecies analysis for Estuary Category#

#Take out OTU aka genus table from metadata focused phyloseq object#
estuary_seqs_genus = as(otu_table(RFLWSS_Y1_GenusAgenus), "matrix")

#Make it a data frame, transpose it, make it a data frame again#
estuary_seqs_genus <- as.data.frame(estuary_seqs_genus)
estuary_seqs_genus<- t(estuary_seqs_genus)
estuary_seqs_genus <- as.data.frame(estuary_seqs_genus)

#Run indicspecies with same estuary_meta_group as in asv section#
estuary_indval_genus = multipatt(estuary_seqs_genus, estuary_meta_group, control = how(nperm=9999))

#Create a blank csv, then sink all indicspecies results into it#
sink("estuary_IS_res_genus.csv")
estuary_sig_indval_genus <- summary(estuary_indval_genus, indvalcomp=TRUE, alpha=1)
sink()

##Conduct genus indicspecies analysis for Muck Category#

#Take out OTU aka genus table from metadata focused phyloseq object#
muck_seqs_genus = as(otu_table(RFLWSS_Y1_MuckvsNot_GenusAgenus), "matrix")

#Make it a data frame, transpose it, make it a data frame again#
muck_seqs_genus <- as.data.frame(muck_seqs_genus)
muck_seqs_genus<- t(muck_seqs_genus)
muck_seqs_genus <- as.data.frame(muck_seqs_genus)

#Run indicspecies with same muck_meta_group as asv section#
muck_indval_genus = multipatt(muck_seqs_genus, muck_meta_group, control = how(nperm=9999))

#Create a blank csv, then sink all indicspecies results into it#
sink("muck_IS_res_genus.csv")
muck_sig_indval_genus <- summary(muck_indval_genus, indvalcomp=TRUE, alpha=1)
sink()


##Conduct genus indicspecies analysis for TOM/Cu Category#

#Take out OTU aka genus table from metadata focused phyloseq object#
tomcu_seqs_genus = as(otu_table(RFLWSS_Y1_HHvsHL_GenusAgenus), "matrix")

#Make it a data frame, transpose it, make it a data frame again#
tomcu_seqs_genus <- as.data.frame(tomcu_seqs_genus)
tomcu_seqs_genus<- t(tomcu_seqs_genus)
tomcu_seqs_genus <- as.data.frame(tomcu_seqs_genus)

#Run indicspecies with same tomcu_meta_group as asv section#
tomcu_indval_genus = multipatt(tomcu_seqs_genus, tomcu_meta_group, control = how(nperm=9999))

#Create a blank csv, then sink all indicspecies results into it#
sink("tomcu_IS_res_genus.csv")
tomcu_sig_indval_genus <- summary(tomcu_indval_genus, indvalcomp=TRUE, alpha=1)
sink()

#Reset max print to default#
options(max.print=99999)

##Look at csv file you sinked into and create a table
##of just genus.ID, A, B, stat, p.value, and group##

##Continue genus indicspecies analysis for Estuary Category##

#Upload new table#
estuary_IS_res_genus <- read.csv(file.choose())
#Make a new p value adjusted column#
estuary_IS_res_genus$padj <- p.adjust(estuary_IS_res_genus$p.value, method = "BH")

#Keep table of adjusted p values less than 0.05#
estuary_IS_sig_res_genus <- filter(estuary_IS_res_genus, padj < 0.05)

#Make dataframe with taxanomy added in by making ESV.ID the rownames then adding taxonomy table#
estuary_IS_sig_res_taxo_genus <- as.data.frame(estuary_IS_sig_res_genus)
rownames(estuary_IS_sig_res_taxo_genus) <- estuary_IS_sig_res_taxo_genus[,1]
estuary_IS_sig_res_taxo_genus = cbind(as(estuary_IS_sig_res_taxo_genus, "data.frame"), as(tax_table(RFLWSS_Y1_GenusAgenus)[rownames(estuary_IS_sig_res_taxo_genus), ], "matrix"))

#Make dataframe with genuss added in from all sites#
estuary_IS_sig_res_taxo_seqs_genus = cbind(as(estuary_IS_sig_res_taxo_genus, "data.frame"), as(otu_table(RFLWSS_Y1_GenusAgenus)[rownames(estuary_IS_sig_res_taxo_genus), ], "matrix"))

#Remove old rownames#
rownames(estuary_IS_sig_res_taxo_seqs_genus) <-NULL

#Make a txt file of results for all sites#
write.csv(as.data.frame(estuary_IS_sig_res_taxo_seqs_genus), file="estuary_IS_sig_res_taxo_seqs_genus.csv")

#Check the number of indicators in each group after padj#
length(which(sample_data(estuary_IS_sig_res_genus)$group == "1")) #146 IRL
length(which(sample_data(estuary_IS_sig_res_genus)$group == "2")) #356 SLE


##Continue genus indicspecies analysis for Muck Category##

#Upload new table#
muck_IS_res_genus <- read.csv(file.choose())

#Make a new p value adjusted column#
muck_IS_res_genus$padj <- p.adjust(muck_IS_res_genus$p.value, method = "BH")

#Keep table of adjusted p values less than 0.05#
muck_IS_sig_res_genus <- filter(muck_IS_res_genus, padj < 0.05)

#Make dataframe with taxanomy added in by making ESV.ID the rownames then adding taxonomy table#
muck_IS_sig_res_taxo_genus <- as.data.frame(muck_IS_sig_res_genus)
rownames(muck_IS_sig_res_taxo_genus) <- muck_IS_sig_res_taxo_genus[,1]
muck_IS_sig_res_taxo_genus = cbind(as(muck_IS_sig_res_taxo_genus, "data.frame"), as(tax_table(RFLWSS_Y1_MuckvsNot_GenusAgenus)[rownames(muck_IS_sig_res_taxo_genus), ], "matrix"))

#Make dataframe with genuss added in from all sites#
muck_IS_sig_res_taxo_seqs_genus = cbind(as(muck_IS_sig_res_taxo_genus, "data.frame"), as(otu_table(RFLWSS_Y1_GenusAgenus)[rownames(muck_IS_sig_res_taxo_genus), ], "matrix"))

#Remove old rownames#
rownames(muck_IS_sig_res_taxo_seqs_genus) <-NULL

#Make a txt file of results for all sites#
write.csv(as.data.frame(muck_IS_sig_res_taxo_seqs_genus), file="muck_IS_sig_res_taxo_seqs_genus.csv")

#Check the number of indicators in each group after padj#
length(which(sample_data(muck_IS_sig_res_genus)$group == "1")) #148 Muck
length(which(sample_data(muck_IS_sig_res_genus)$group == "2")) #62 Not

##Continue genus indicspecies analysis for TOM/Cu Category##

#Upload new table#
tomcu_IS_res_genus <- read.csv(file.choose())
#Make a new p value adjusted column#
tomcu_IS_res_genus$padj <- p.adjust(tomcu_IS_res_genus$p.value, method = "BH")

#Keep table of adjusted p values less than 0.05#
tomcu_IS_sig_res_genus <- filter(tomcu_IS_res_genus, padj < 0.05)

#Make dataframe with taxanomy added in by making ESV.ID the rownames then adding taxonomy table#
tomcu_IS_sig_res_taxo_genus <- as.data.frame(tomcu_IS_sig_res_genus)
rownames(tomcu_IS_sig_res_taxo_genus) <- tomcu_IS_sig_res_taxo_genus[,1]
tomcu_IS_sig_res_taxo_genus = cbind(as(tomcu_IS_sig_res_taxo_genus, "data.frame"), as(tax_table(RFLWSS_Y1_MuckvsNot_GenusAgenus)[rownames(tomcu_IS_sig_res_taxo_genus), ], "matrix"))

#Make dataframe with genuss added in from all sites#
tomcu_IS_sig_res_taxo_seqs_genus = cbind(as(tomcu_IS_sig_res_taxo_genus, "data.frame"), as(otu_table(RFLWSS_Y1_GenusAgenus)[rownames(tomcu_IS_sig_res_taxo_genus), ], "matrix"))

#Remove old rownames#
rownames(tomcu_IS_sig_res_taxo_seqs_genus) <-NULL

#Make a txt file of results for all sites#
write.csv(as.data.frame(tomcu_IS_sig_res_taxo_seqs_genus), file="tomcu_IS_sig_res_taxo_seqs_genus.csv")

#Check the number of indicators in each group after padj#
length(which(sample_data(tomcu_IS_sig_res_genus)$group == "1")) #84 High/High
length(which(sample_data(tomcu_IS_sig_res_genus)$group == "2")) #92 High/Low

#Order rows by group#
estuary_IS_sig_res_taxo_seqs_genus <- arrange(estuary_IS_sig_res_taxo_seqs_genus, group)
muck_IS_sig_res_taxo_seqs_genus <- arrange(muck_IS_sig_res_taxo_seqs_genus, group)
tomcu_IS_sig_res_taxo_seqs_genus <- arrange(tomcu_IS_sig_res_taxo_seqs_genus, group)


#Create lists of just significant genuss for each category#
unfiltered_IS_estuary_genuss <- subset(estuary_IS_sig_res_taxo_seqs_genus, select=c(ESV.ID))
unfiltered_IS_muck_genuss <- subset(muck_IS_sig_res_taxo_seqs_genus, select=c(ESV.ID))
unfiltered_IS_tomcu_genuss <- subset(tomcu_IS_sig_res_taxo_seqs_genus, select=c(ESV.ID))


##Create indicator tables filtered of overlapping indicators with other two metadata categories##

#Filter out rows in the estuary IS object that are also in the tomcu and muck IS objects#
filtered_estuary_IS_sig_res_taxo_seqs_genus <- anti_join(estuary_IS_sig_res_taxo_seqs_genus, unfiltered_IS_tomcu_genuss, by = "ESV.ID")
filtered_estuary_IS_sig_res_taxo_seqs_genus <- anti_join(filtered_estuary_IS_sig_res_taxo_seqs_genus, unfiltered_IS_muck_genuss, by = "ESV.ID")

#Filter out rows in the muck IS object that are also in the estuary and tomcu IS objects#
filtered_muck_IS_sig_res_taxo_seqs_genus<- anti_join(muck_IS_sig_res_taxo_seqs_genus, unfiltered_IS_estuary_genuss, by = "ESV.ID")
filtered_muck_IS_sig_res_taxo_seqs_genus <- anti_join(filtered_muck_IS_sig_res_taxo_seqs_genus, unfiltered_IS_tomcu_genuss, by = "ESV.ID")

#Filter out rows in the tomcu IS object that are also in the estuary and muck IS objects#
filtered_tomcu_IS_sig_res_taxo_seqs_genus <- anti_join(tomcu_IS_sig_res_taxo_seqs_genus, unfiltered_IS_estuary_genuss, by = "ESV.ID")
filtered_tomcu_IS_sig_res_taxo_seqs_genus <- anti_join(filtered_tomcu_IS_sig_res_taxo_seqs_genus, unfiltered_IS_muck_genuss, by = "ESV.ID")

#Order rows by group#
filtered_estuary_IS_sig_res_taxo_seqs_genus <- arrange(filtered_estuary_IS_sig_res_taxo_seqs_genus, group)
filtered_muck_IS_sig_res_taxo_seqs_genus <- arrange(filtered_muck_IS_sig_res_taxo_seqs_genus, group)
filtered_tomcu_IS_sig_res_taxo_seqs_genus <- arrange(filtered_tomcu_IS_sig_res_taxo_seqs_genus, group)


#Make a csv file for each of the filtered tables#
write.csv(as.data.frame(filtered_estuary_IS_sig_res_taxo_seqs_genus), file="filtered_estuary_IS_sig_res_taxo_seqs_genus.csv")
write.csv(as.data.frame(filtered_muck_IS_sig_res_taxo_seqs_genus), file="filtered_muck_IS_sig_res_taxo_seqs_genus.csv")
write.csv(as.data.frame(filtered_tomcu_IS_sig_res_taxo_seqs_genus), file="filtered_tomcu_IS_sig_res_taxo_seqs_genus.csv")

#Determine the number of indicators in each category#
length(which(sample_data(filtered_estuary_IS_sig_res_taxo_seqs_genus)$group == "1")) #75 IRL
length(which(sample_data(filtered_estuary_IS_sig_res_taxo_seqs_genus)$group == "2")) #240 SLE
length(which(sample_data(filtered_muck_IS_sig_res_taxo_seqs_genus)$group == "1")) #78 Muck
length(which(sample_data(filtered_muck_IS_sig_res_taxo_seqs_genus)$group == "2")) #28 Not
length(which(sample_data(filtered_tomcu_IS_sig_res_taxo_seqs_genus)$group == "1")) #17 High/High
length(which(sample_data(filtered_tomcu_IS_sig_res_taxo_seqs_genus)$group == "2")) #10 High/Low


#Create lists of filtered significant genuss for each category#
filtered_IS_estuary_genuss <- subset(filtered_estuary_IS_sig_res_taxo_seqs_genus, select=c(ESV.ID))
filtered_IS_muck_genuss <- subset(filtered_muck_IS_sig_res_taxo_seqs_genus, select=c(ESV.ID))
filtered_IS_tomcu_genuss <- subset(filtered_tomcu_IS_sig_res_taxo_seqs_genus, select=c(ESV.ID))

#Create phyloseq objects of the unfiltered indicators using character string version of significant genuss##
charac_unfiltered_IS_estuary_genuss <- unfiltered_IS_estuary_genuss[,"ESV.ID"]
charac_unfiltered_IS_estuary_genuss <- as.character(charac_unfiltered_IS_estuary_genuss)
UFISestuaryGenusAgenus <- prune_taxa(charac_unfiltered_IS_estuary_genuss, RFLWSS_Y1_GenusAgenus)

charac_unfiltered_IS_muck_genuss <- unfiltered_IS_muck_genuss[,"ESV.ID"]
charac_unfiltered_IS_muck_genuss <- as.character(charac_unfiltered_IS_muck_genuss)
UFISmuckGenusAgenus <- prune_taxa(charac_unfiltered_IS_muck_genuss, RFLWSS_Y1_GenusAgenus)

charac_unfiltered_IS_tomcu_genuss <- unfiltered_IS_tomcu_genuss[,"ESV.ID"]
charac_unfiltered_IS_tomcu_genuss <- as.character(charac_unfiltered_IS_tomcu_genuss)
UFIStomcuGenusAgenus <- prune_taxa(charac_unfiltered_IS_tomcu_genuss, RFLWSS_Y1_GenusAgenus)

#Create phyloseq objects of the filtered indicators using character string version of significant genuss##
charac_filtered_IS_estuary_genuss <- filtered_IS_estuary_genuss[,"ESV.ID"]
charac_filtered_IS_estuary_genuss <- as.character(charac_filtered_IS_estuary_genuss)
FISestuaryGenusAgenus <- prune_taxa(charac_filtered_IS_estuary_genuss, RFLWSS_Y1_GenusAgenus)

charac_filtered_IS_muck_genuss <- filtered_IS_muck_genuss[,"ESV.ID"]
charac_filtered_IS_muck_genuss <- as.character(charac_filtered_IS_muck_genuss)
FISmuckGenusAgenus <- prune_taxa(charac_filtered_IS_muck_genuss, RFLWSS_Y1_GenusAgenus)

charac_filtered_IS_tomcu_genuss <- filtered_IS_tomcu_genuss[,"ESV.ID"]
charac_filtered_IS_tomcu_genuss <- as.character(charac_filtered_IS_tomcu_genuss)
FIStomcuGenusAgenus <- prune_taxa(charac_filtered_IS_tomcu_genuss, RFLWSS_Y1_GenusAgenus)

##Create heatmaps to see overall trends, may take a while to load##

unfiltered_IS_estuary_genus_heatmap <- plot_heatmap(UFDestuaryGenusAgenus, taxa.order = charac_unfiltered_IS_estuary_genuss, sample.order= "Estuary", sample.label = "Estuary")
unfiltered_IS_estuary_genus_heatmap

unfiltered_IS_muck_genus_heatmap <- plot_heatmap(UFDmuckGenusAgenus, taxa.order = charac_unfiltered_IS_muck_genuss, sample.order= "Muck", sample.label = "Muck")
unfiltered_IS_muck_genus_heatmap

unfiltered_IS_tomcu_genus_heatmap <- plot_heatmap(UFDtomcuGenusAgenus, taxa.order = charac_unfiltered_IS_tomcu_genuss, sample.order= "LOI.Cu", sample.label = "LOI.Cu")
unfiltered_IS_tomcu_genus_heatmap

filtered_IS_estuary_genus_heatmap <- plot_heatmap(FDestuaryGenusAgenus, taxa.order = charac_filtered_IS_estuary_genuss, sample.order= "Estuary", sample.label = "Estuary")
filtered_IS_estuary_genus_heatmap

filtered_IS_muck_genus_heatmap <- plot_heatmap(FDmuckGenusAgenus, taxa.order = charac_filtered_IS_muck_genuss, sample.order= "Muck", sample.label = "Muck")
filtered_IS_muck_genus_heatmap

filtered_IS_tomcu_genus_heatmap <- plot_heatmap(FDtomcuGenusAgenus, taxa.order = charac_filtered_IS_tomcu_genuss, sample.order= "LOI.Cu", sample.label = "LOI.Cu")
filtered_IS_tomcu_genus_heatmap

###MAKE COMBO INDICATOR LISTS###

##Make "Combo" indicator lists from overlapping unfiltered DESeq2 and IS indicators##

#Keep rows in the estuary DESEQ object that are also in the IRLSLE IS object#
unfiltered_estuary_combo_sig_res_taxo_seqs_genus <- merge(estuary_IS_sig_res_genus,estuary_DESeq2_sig_res_taxo_seqs_genus, by="ESV.ID")

#Keep rows in the muck DESEQ object that are also in the Muck IS object#
unfiltered_muck_combo_sig_res_taxo_seqs_genus <- merge(muck_IS_sig_res_genus,muck_DESeq2_sig_res_taxo_seqs_genus, by="ESV.ID")

#Keep rows in the tomcu DESEQ object that are also in the LOICu IS object#
unfiltered_tomcu_combo_sig_res_taxo_seqs_genus <- merge(tomcu_IS_sig_res_genus,tomcu_DESeq2_sig_res_taxo_seqs_genus, by="ESV.ID")

##Make "Combo" indicator lists from overlapping filtered DESeq2 and IS indicators##

#Create focused filtered IS tables without taxonomy or sequences#
filtered_estuary_IS_sig_res_genus <- subset(filtered_estuary_IS_sig_res_taxo_seqs_genus, select=c(ESV.ID,A,B,stat,p.value,group,padj))

filtered_muck_IS_sig_res_genus <- subset(filtered_muck_IS_sig_res_taxo_seqs_genus, select=c(ESV.ID,A,B,stat,p.value,group,padj))

filtered_tomcu_IS_sig_res_genus <- subset(filtered_tomcu_IS_sig_res_taxo_seqs_genus, select=c(ESV.ID,A,B,stat,p.value,group,padj))

#Keep rows in the estuary DESEQ object that are also in the IRLSLE IS object#
filtered_estuary_combo_sig_res_taxo_seqs_genus <- merge(filtered_estuary_IS_sig_res_genus,filtered_estuary_DESeq2_sig_res_taxo_seqs_genus, by="ESV.ID")

#Keep rows in the muck DESEQ object that are also in the Muck IS object#
filtered_muck_combo_sig_res_taxo_seqs_genus <- merge(filtered_muck_IS_sig_res_genus,filtered_muck_DESeq2_sig_res_taxo_seqs_genus, by="ESV.ID")

#Keep rows in the tomcu DESEQ object that are also in the LOICu IS object#
filtered_tomcu_combo_sig_res_taxo_seqs_genus <- merge(filtered_tomcu_IS_sig_res_genus,filtered_tomcu_DESeq2_sig_res_taxo_seqs_genus, by="ESV.ID")

#Determine the number of indicators in each subcategory of the unfiltered Combo tables#
length(which(sample_data(unfiltered_estuary_combo_sig_res_taxo_seqs_genus)$group == "1")) #128 IRL
length(which(sample_data(unfiltered_estuary_combo_sig_res_taxo_seqs_genus)$group == "2")) #179 SLE
length(which(sample_data(unfiltered_muck_combo_sig_res_taxo_seqs_genus)$group == "1")) #43 Muck
length(which(sample_data(unfiltered_muck_combo_sig_res_taxo_seqs_genus)$group == "2")) #59 Not
length(which(sample_data(unfiltered_tomcu_combo_sig_res_taxo_seqs_genus)$group == "1")) #47 High/High
length(which(sample_data(unfiltered_tomcu_combo_sig_res_taxo_seqs_genus)$group == "2")) #91 High/Low

#Determine the number of indicators in each subcategory of the filtered Combo tables#
length(which(sample_data(filtered_estuary_combo_sig_res_taxo_seqs_genus)$group == "1")) #43 IRL
length(which(sample_data(filtered_estuary_combo_sig_res_taxo_seqs_genus)$group == "2")) #62 SLE
length(which(sample_data(filtered_muck_combo_sig_res_taxo_seqs_genus)$group == "1")) #18 Muck
length(which(sample_data(filtered_muck_combo_sig_res_taxo_seqs_genus)$group == "2")) #15 Not
length(which(sample_data(filtered_tomcu_combo_sig_res_taxo_seqs_genus)$group == "1")) #3 High/High
length(which(sample_data(filtered_tomcu_combo_sig_res_taxo_seqs_genus)$group == "2")) #7 High/Low


#Order unfiltered Combo table rows by group#
unfiltered_estuary_combo_sig_res_taxo_seqs_genus <- arrange(unfiltered_estuary_combo_sig_res_taxo_seqs_genus, group)
unfiltered_muck_combo_sig_res_taxo_seqs_genus <- arrange(unfiltered_muck_combo_sig_res_taxo_seqs_genus, group)
unfiltered_tomcu_combo_sig_res_taxo_seqs_genus <- arrange(unfiltered_tomcu_combo_sig_res_taxo_seqs_genus, group)

#Order filtered Combo table rows by group#
filtered_estuary_combo_sig_res_taxo_seqs_genus <- arrange(filtered_estuary_combo_sig_res_taxo_seqs_genus, group)
filtered_muck_combo_sig_res_taxo_seqs_genus <- arrange(filtered_muck_combo_sig_res_taxo_seqs_genus, group)
filtered_tomcu_combo_sig_res_taxo_seqs_genus <- arrange(filtered_tomcu_combo_sig_res_taxo_seqs_genus, group)

#Make a csv file for each of the unfiltered Combo tables#
write.csv(as.data.frame(unfiltered_estuary_combo_sig_res_taxo_seqs_genus), file="unfiltered_estuary_combo_sig_res_taxo_seqs_genus.csv")
write.csv(as.data.frame(unfiltered_muck_combo_sig_res_taxo_seqs_genus), file="unfiltered_muck_combo_sig_res_taxo_seqs_genus.csv")
write.csv(as.data.frame(unfiltered_tomcu_combo_sig_res_taxo_seqs_genus), file="unfiltered_tomcu_combo_sig_res_taxo_seqs_genus.csv")

#Make a csv file for each of the filtered Combo tables#
write.csv(as.data.frame(filtered_estuary_combo_sig_res_taxo_seqs_genus), file="filtered_estuary_combo_sig_res_taxo_seqs_genus.csv")
write.csv(as.data.frame(filtered_muck_combo_sig_res_taxo_seqs_genus), file="filtered_muck_combo_sig_res_taxo_seqs_genus.csv")
write.csv(as.data.frame(filtered_tomcu_combo_sig_res_taxo_seqs_genus), file="filtered_tomcu_combo_sig_res_taxo_seqs_genus.csv")

#Create lists of unfiltered significant genuss for each category#
unfiltered_Combo_estuary_genuss <- subset(unfiltered_estuary_combo_sig_res_taxo_seqs_genus, select=c(ESV.ID))
unfiltered_Combo_muck_genuss <- subset(unfiltered_muck_combo_sig_res_taxo_seqs_genus, select=c(ESV.ID))
unfiltered_Combo_tomcu_genuss <- subset(unfiltered_tomcu_combo_sig_res_taxo_seqs_genus, select=c(ESV.ID))

#Create lists of filtered significant genuss for each category#
filtered_Combo_estuary_genuss <- subset(filtered_estuary_combo_sig_res_taxo_seqs_genus, select=c(ESV.ID))
filtered_Combo_muck_genuss <- subset(filtered_muck_combo_sig_res_taxo_seqs_genus, select=c(ESV.ID))
filtered_Combo_tomcu_genuss <- subset(filtered_tomcu_combo_sig_res_taxo_seqs_genus, select=c(ESV.ID))

##Create phyloseq objects for Combo indicators##

#Create phyloseq objects of the unfiltered indicators using character string version of significant genuss##
charac_unfiltered_Combo_estuary_genuss <- unfiltered_Combo_estuary_genuss[,"ESV.ID"]
charac_unfiltered_Combo_estuary_genuss <- as.character(charac_unfiltered_Combo_estuary_genuss)
UFCestuaryGenusAgenus <- prune_taxa(charac_unfiltered_Combo_estuary_genuss, RFLWSS_Y1_GenusAgenus)

charac_unfiltered_Combo_muck_genuss <- unfiltered_Combo_muck_genuss[,"ESV.ID"]
charac_unfiltered_Combo_muck_genuss <- as.character(charac_unfiltered_Combo_muck_genuss)
UFCmuckGenusAgenus <- prune_taxa(charac_unfiltered_Combo_muck_genuss, RFLWSS_Y1_GenusAgenus)

charac_unfiltered_Combo_tomcu_genuss <- unfiltered_Combo_tomcu_genuss[,"ESV.ID"]
charac_unfiltered_Combo_tomcu_genuss <- as.character(charac_unfiltered_Combo_tomcu_genuss)
UFCtomcuGenusAgenus <- prune_taxa(charac_unfiltered_Combo_tomcu_genuss, RFLWSS_Y1_GenusAgenus)

#Create phyloseq objects of the filtered indicators using character string version of significant genuss##
charac_filtered_Combo_estuary_genuss <- filtered_Combo_estuary_genuss[,"ESV.ID"]
charac_filtered_Combo_estuary_genuss <- as.character(charac_filtered_Combo_estuary_genuss)
FCestuaryGenusAgenus <- prune_taxa(charac_filtered_Combo_estuary_genuss, RFLWSS_Y1_GenusAgenus)

charac_filtered_Combo_muck_genuss <- filtered_Combo_muck_genuss[,"ESV.ID"]
charac_filtered_Combo_muck_genuss <- as.character(charac_filtered_Combo_muck_genuss)
FCmuckGenusAgenus <- prune_taxa(charac_filtered_Combo_muck_genuss, RFLWSS_Y1_GenusAgenus)

charac_filtered_Combo_tomcu_genuss <- filtered_Combo_tomcu_genuss[,"ESV.ID"]
charac_filtered_Combo_tomcu_genuss <- as.character(charac_filtered_Combo_tomcu_genuss)
FCtomcuGenusAgenus <- prune_taxa(charac_filtered_Combo_tomcu_genuss, RFLWSS_Y1_GenusAgenus)

##Create heatmaps to see overall trends, may take a while to load##

unfiltered_Combo_estuary_genus_heatmap <- plot_heatmap(UFDestuaryGenusAgenus, taxa.order = charac_unfiltered_Combo_estuary_genuss, sample.order= "Estuary", sample.label = "Estuary")
unfiltered_Combo_estuary_genus_heatmap

unfiltered_Combo_muck_genus_heatmap <- plot_heatmap(UFDmuckGenusAgenus, taxa.order = charac_unfiltered_Combo_muck_genuss, sample.order= "Muck", sample.label = "Muck")
unfiltered_Combo_muck_genus_heatmap

unfiltered_Combo_tomcu_genus_heatmap <- plot_heatmap(UFDtomcuGenusAgenus, taxa.order = charac_unfiltered_Combo_tomcu_genuss, sample.order= "LOI.Cu", sample.label = "LOI.Cu")
unfiltered_Combo_tomcu_genus_heatmap

filtered_Combo_estuary_genus_heatmap <- plot_heatmap(FDestuaryGenusAgenus, taxa.order = charac_filtered_Combo_estuary_genuss, sample.order= "Estuary", sample.label = "Estuary")
filtered_Combo_estuary_genus_heatmap

filtered_Combo_muck_genus_heatmap <- plot_heatmap(FDmuckGenusAgenus, taxa.order = charac_filtered_Combo_muck_genuss, sample.order= "Muck", sample.label = "Muck")
filtered_Combo_muck_genus_heatmap

filtered_Combo_tomcu_genus_heatmap <- plot_heatmap(FDtomcuGenusAgenus, taxa.order = charac_filtered_Combo_tomcu_genuss, sample.order= "LOI.Cu", sample.label = "LOI.Cu")
filtered_Combo_tomcu_genus_heatmap

###INDICATOR EFFECTIVENESS METHOD###

##NOTES ON METHOD##
#"affected' = subcategory that is more affected by the stressor you are      concerned about. For example SLE for the Estuary category becasue it has     more freshwater discharges, 3MC for the Muck category because it has more    muck characteristics, HiHi for the TOM/Cu category becasue it has more       copper#

#"non.affected" = IRL, 0MC, HiLo#

#microbially.predicted "affected" sample = samples in the partitioning around medoids (PAM) cluster with the most metadata.defined "affected" samples#

#microbially.predicted "non.affected" sample = samples in the PAM cluster with the most metadata.defined "non.affected" samples#

#Overall Idea: Use PAM clustering to split samples between the "affected" and "non.affected" clusters based upon the indicators (microbially.predicted) and see how well this classification matches the metadata.defined classification by using the product of specificity and sensitivity#

#Copy results to an Excel Sheet for further statistical testing, each combo of factors should have four percentage types (Product is Sensitivity x Specificity), for example:#
#Indicator_Test	Taxonomic_Level	Filtering_Status	Metadata_Tested	Percent_Type	Percentage#
#Indicgenus   genus             Unfiltered        Estuary         Total         98.78#
#Indicgenus   genus             Unfiltered        Estuary         Sensitivity   96.23#
#Indicgenus   genus             Unfiltered        Estuary         Specificity   99.15#
#Indicgenus   genus             Unfiltered        Estuary         Product       95.41#


###ORIGINAL DATASETS
###Split samples of orginal datasets into two clusters to see how it will compare to splitting based upon indicators###

##Test on Muck and Estuary categories at genus level##

#Remove any samples with no genuss from target phyloseq#
nsamples(RFLWSS_Y1_GenusAgenus) #90
RFLWSS_Y1_GenusAgenus = prune_samples(sample_sums(RFLWSS_Y1_GenusAgenus)>=1, RFLWSS_Y1_GenusAgenus)
nsamples(RFLWSS_Y1_GenusAgenus) #90

#Make a matrix with the sample_data and relabel subcategories to make metadata.defined column#
Original_Est_Muck_genus <- as.matrix(sample_data(RFLWSS_Y1_GenusAgenus))
Original_Est_Muck_genus[Original_Est_Muck_genus ==  "Muck"] <- 2
Original_Est_Muck_genus[Original_Est_Muck_genus ==  "Not"] <- 1
Original_Est_Muck_genus[Original_Est_Muck_genus ==  "Mucky"] <- 2
Original_Est_Muck_genus[Original_Est_Muck_genus ==  "Muckish"] <- 2
Original_Est_Muck_genus[Original_Est_Muck_genus ==  "IRL"] <- 1
Original_Est_Muck_genus[Original_Est_Muck_genus ==  "SLE"] <- 2
Original_Est_Muck_genus<-as.data.frame(Original_Est_Muck_genus)


#Create Bray.Curtis distance matrix from genus table#
bc <- vegdist(t(otu_table(RFLWSS_Y1_GenusAgenus)))

#Conduct pam clustering based upon 2 clusters#
pam.res <- pam(bc,2, diss=TRUE)

#Add a microbially.predicted column (Orginal_genus) to matrix above for comparison to metadata.defined column#
Original_Est_Muck_genus <- cbind(Original_Est_Muck_genus, Orginal_genus = pam.res$cluster)

##Test Estuary category##

#Test  to see how well the metadata.defined and microbially.predicted columns match one another aka Total Efficiency; if the number is below 0.5, then switch the numbers in the metadata.defiend column#

#Test for Estuary Effeciency# 0.9444444
sum(Original_Est_Muck_genus$Estuary == Original_Est_Muck_genus$Orginal_genus)/nrow(Original_Est_Muck_genus)

##Test indicator sensitivity aka the true postives (TP)/ (TP + false negatives (FN))##
#TP = metadata.defined "affected" sample correctly placed in the microbiall   y.predicted "affected" cluster#
#FN = metadata.defined "affected" sample incorrectly placed in the           microbially.predicted "non.affected" cluster#
#TP+FN = all metadata.defined "affected" samples

#Test for SLE Effeciency aka Sensitvity# 1
Orginal_SLE_genus <- filter(Original_Est_Muck_genus, Estuary == "2")
sum(Orginal_SLE_genus$Estuary == Orginal_SLE_genus$Orginal_genus)/nrow(Orginal_SLE_genus)

##Test indicator specificity aka the true negatives (TN)/ (TN + false positives (FP))##
#TP = metadata.defined "affected" sample correctly placed in the microbiall   y.predicted "affected" cluster#
#FN = metadata.defined "affected" sample incorrectly placed in the           microbially.predicted "non.affected" cluster#
#TP+FN = all metadata.defined "affected" samples

#Test for IRL Effeciency aka Specificity# 0.9242424
Orginal_IRL_genus <- filter(Original_Est_Muck_genus, Estuary == "1")
sum(Orginal_IRL_genus$Estuary == Orginal_IRL_genus$Orginal_genus)/nrow(Orginal_IRL_genus)

#PCoA with Color by Estuary#
TotalPCoA <- ordinate(RFLWSS_Y1_GenusAgenus,"PCoA")
p = plot_ordination(RFLWSS_Y1_GenusAgenus, TotalPCoA, color="Estuary") + 
  ggtitle("Total genuss separated by Estuary Category PCoA with Bray.Curtis distance") +
  geom_text(aes(label=Original_Est_Muck_genus$Orginal_genus),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p

##Test Muck category##

#Test for Muck Total Efficiency# 0.7222222
sum(Original_Est_Muck_genus$Muck == Original_Est_Muck_genus$Orginal_genus)/nrow(Original_Est_Muck_genus)

#Test for 3MC Effeciency aka Sensitivity# 0.5666667
Orginal_3MC_genus <- filter(Original_Est_Muck_genus, Muck == "2")
sum(Orginal_3MC_genus$Muck == Orginal_3MC_genus$Orginal_genus)/nrow(Orginal_3MC_genus)

#Test for 0MC Effeciency aka Specificity# 0.8
Orginal_0MC_genus <- filter(Original_Est_Muck_genus, Muck == "1")
sum(Orginal_0MC_genus$Muck == Orginal_0MC_genus$Orginal_genus)/nrow(Orginal_0MC_genus)

#PCoA with Color by Muck#
TotalPCoA <- ordinate(RFLWSS_Y1_GenusAgenus,"PCoA")
p = plot_ordination(RFLWSS_Y1_GenusAgenus, TotalPCoA, color="Muck") + 
  ggtitle("Total genuss separated by Muck Category PCoA with Bray.Curtis distance") +
  geom_text(aes(label=Original_Est_Muck_genus$Orginal_genus),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p + annotate("text", x = 0.3, y = 0.3, label = c("SN = 0.44, SP=0.84"))+ scale_color_discrete(name="Muck Characteristics", limits=c("Muck", "Mucky", "Muckish", "Not"), labels=c("3", "2", "1", "0"))

##Test TOM/Cu category##

#Remove any samples with no genuss from target phyloseq#
nsamples(RFLWSS_Y1_HHvsHL_GenusAgenus) #29
RFLWSS_Y1_HHvsHL_GenusAgenus = prune_samples(sample_sums(RFLWSS_Y1_HHvsHL_GenusAgenus)>=1, RFLWSS_Y1_HHvsHL_GenusAgenus)
nsamples(RFLWSS_Y1_HHvsHL_GenusAgenus) #29

#Make a matrix with the sample_data and relabel subcategories to make metadata.defined column#
Original_TOMCu_genus <- as.matrix(sample_data(RFLWSS_Y1_HHvsHL_GenusAgenus))
Original_TOMCu_genus[Original_TOMCu_genus ==  "High.High"] <- 1
Original_TOMCu_genus[Original_TOMCu_genus ==  "High.Low"] <- 2
Original_TOMCu_genus<-as.data.frame(Original_TOMCu_genus)

#Create Bray.Curtis distance matrix from genus table#
bc <- vegdist(t(otu_table(RFLWSS_Y1_HHvsHL_GenusAgenus)))

#Conduct PAM clustering based upon 2 clusters#
pam.res <- pam(bc,2, diss=TRUE)

#Add a microbially.predicted column to matrix above for comparison to metadata.defined column#
Original_TOMCu_genus <- cbind(Original_TOMCu_genus, Orginal_genus = pam.res$cluster)

#Determine Total TOM/Cu Efficiency# 0.6551724
sum(Original_TOMCu_genus$LOI.Cu == Original_TOMCu_genus$Orginal_genus)/nrow(Original_TOMCu_genus)

#Determine HiHi Efficiency  aka Sensitivity# 0.5454545
Original_TOMCu_genusHH <- filter(Original_TOMCu_genus, LOI.Cu == "1") 
sum(Original_TOMCu_genusHH$LOI.Cu == Original_TOMCu_genusHH$Orginal_genus)/nrow(Original_TOMCu_genusHH)

#Determine HiLo Efficiency aka Specificity# 0.7222222
Original_TOMCu_genusHL <- filter(Original_TOMCu_genus, LOI.Cu == "2") 
sum(Original_TOMCu_genusHL$LOI.Cu == Original_TOMCu_genusHL$Orginal_genus)/nrow(Original_TOMCu_genusHL)

#PCoA with Color by TOM/Cu#
TotalPCoA <- ordinate(RFLWSS_Y1_HHvsHL_GenusAgenus,"PCoA")
p = plot_ordination(RFLWSS_Y1_HHvsHL_GenusAgenus, TotalPCoA, color="LOI.Cu") + 
  ggtitle("Total genuss separated by LOI.Cu Category PCoA with Bray.Curtis distance") +
  geom_text(aes(label=Original_TOMCu_genus$Orginal_genus),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p

###TEST INDICSPECIES ESTUARY INDICATORS###

##Test unfiltered Estuary genuss##

#Remove any samples with no genuss from target phyloseq#
nsamples(UFISestuaryGenusAgenus) #90
UFISestuaryGenusAgenus = prune_samples(sample_sums(UFISestuaryGenusAgenus)>=1, UFISestuaryGenusAgenus)
nsamples(UFISestuaryGenusAgenus) #90

#Make a matrix with the sample_data and relabel subcategories to make metadata.defined column#
UF_IS_Est_genus <- as.matrix(sample_data(UFISestuaryGenusAgenus))
UF_IS_Est_genus[UF_IS_Est_genus ==  "SLE"] <- 2
UF_IS_Est_genus[UF_IS_Est_genus ==  "IRL"] <- 1
UF_IS_Est_genus<-as.data.frame(UF_IS_Est_genus)

#Create Bray.Curtis distance matrix from genus table#
bc <- vegdist(t(otu_table(UFISestuaryGenusAgenus)))

#Conduct PAM clustering based upon 2 clusters#
pam.res <- pam(bc,2, diss=TRUE)

#Add a microbially.predicted column to matrix above for comparison to metadata.defined column#
UF_IS_Est_genus <- cbind(UF_IS_Est_genus, MP_UFISEstgenus = pam.res$cluster)

#Test Total Estuary Efficiency# 1
sum(UF_IS_Est_genus$Estuary == UF_IS_Est_genus$MP_UFISEstgenus)/nrow(UF_IS_Est_genus)

#Test SLE Efficiency aka Sensitvity# 1
UF_IS_SLE_genus <- filter(UF_IS_Est_genus, Estuary == "2")
sum(UF_IS_SLE_genus$Estuary == UF_IS_SLE_genus$MP_UFISEstgenus)/nrow(UF_IS_SLE_genus)

#Test IRL Efficiency aka Specificity# 1
UF_IS_IRL_genus <- filter(UF_IS_Est_genus, Estuary == "1")
sum(UF_IS_IRL_genus$Estuary == UF_IS_IRL_genus$MP_UFISEstgenus)/nrow(UF_IS_IRL_genus)

TotalPCoA <- ordinate(UFISestuaryGenusAgenus,"PCoA")
p = plot_ordination(UFISestuaryGenusAgenus, TotalPCoA, color="Estuary") + 
  ggtitle("Unfiltered IS Estuary genus Indicators PCoA with Bray.Curtis distance") +
  geom_text(aes(label=UF_IS_Est_genus$MP_UFISEstgenus),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p


##Test filtered Estuary genuss##

#Remove any samples with no genuss from target phyloseq#
nsamples(FISestuaryGenusAgenus) #90
FISestuaryGenusAgenus = prune_samples(sample_sums(FISestuaryGenusAgenus)>=1, FISestuaryGenusAgenus)
nsamples(FISestuaryGenusAgenus) #90

#Make a matrix with the sample_data and relabel subcategories to make metadata.defined column#
F_IS_Est_genus <- as.matrix(sample_data(FISestuaryGenusAgenus))
F_IS_Est_genus[F_IS_Est_genus ==  "SLE"] <- 2
F_IS_Est_genus[F_IS_Est_genus ==  "IRL"] <- 1
F_IS_Est_genus<-as.data.frame(F_IS_Est_genus)


#Create Bray.Curtis distance matrix from genus table#
bc <- vegdist(t(otu_table(FISestuaryGenusAgenus)))

#Conduct PAM clustering based upon 2 clusters#
pam.res <- pam(bc,2, diss=TRUE)

#Add a microbially.predicted column to matrix above for comparison to metadata.defined column#
F_IS_Est_genus <- cbind(F_IS_Est_genus, MP_FISEstgenus = pam.res$cluster)

#Test Total Estuary Efficiency# 0.9777778
sum(F_IS_Est_genus$Estuary == F_IS_Est_genus$MP_FISEstgenus)/nrow(F_IS_Est_genus)
#Test SLE Efficiency aka Sensitivity# 0.9166667
F_IS_SLE_genus <- filter(F_IS_Est_genus, Estuary == "2")
sum(F_IS_SLE_genus$Estuary == F_IS_SLE_genus$MP_FISEstgenus)/nrow(F_IS_SLE_genus)
#Test IRL Efficiency aka Specificity# 1
F_IS_IRL_genus <- filter(F_IS_Est_genus, Estuary == "1")
sum(F_IS_IRL_genus$Estuary == F_IS_IRL_genus$MP_FISEstgenus)/nrow(F_IS_IRL_genus)

TotalPCoA <- ordinate(FISestuaryGenusAgenus,"PCoA")
p = plot_ordination(FISestuaryGenusAgenus, TotalPCoA, color="Estuary") + 
  ggtitle("Filtered IS Estuary genus Indicators PCoA with Bray.Curtis distance") +
  geom_text(aes(label=F_IS_Est_genus$MP_FISEstgenus),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p


###TEST DESEQ2 ESTUARY INDICATORS###

##Test unfiltered Estuary genuss##

#Remove any samples with no genuss from target phyloseq#
nsamples(UFDestuaryGenusAgenus) #90
UFDestuaryGenusAgenus = prune_samples(sample_sums(UFDestuaryGenusAgenus)>=1, UFDestuaryGenusAgenus)
nsamples(UFDestuaryGenusAgenus) #90

#Make a matrix with the sample_data and relabel subcategories to make metadata.defined column#
UF_D_Est_genus <- as.matrix(sample_data(UFDestuaryGenusAgenus))
UF_D_Est_genus[UF_D_Est_genus ==  "SLE"] <- 2
UF_D_Est_genus[UF_D_Est_genus ==  "IRL"] <- 1
UF_D_Est_genus<-as.data.frame(UF_D_Est_genus)


#Create Bray.Curtis distance matrix from genus table#
bc <- vegdist(t(otu_table(UFDestuaryGenusAgenus)))

#Conduct PAM clustering based upon 2 clusters#
pam.res <- pam(bc,2, diss=TRUE)

#Add a microbially.predicted column to matrix above for comparison to metadata.defined column#
UF_D_Est_genus <- cbind(UF_D_Est_genus, MP_UFDEstgenus = pam.res$cluster)

#Test Total Estuary Efficiency# 0.9666667
sum(UF_D_Est_genus$Estuary == UF_D_Est_genus$MP_UFDEstgenus)/nrow(UF_D_Est_genus)

#Test SLE Efficiency aka Sensitivity# 0.9583333
UF_D_SLE_genus <- filter(UF_D_Est_genus, Estuary == "2")
sum(UF_D_SLE_genus$Estuary == UF_D_SLE_genus$MP_UFDEstgenus)/nrow(UF_D_SLE_genus)

#Test IRL Efficiency aka Specificity# 0.969697
UF_D_IRL_genus <- filter(UF_D_Est_genus, Estuary == "1")
sum(UF_D_IRL_genus$Estuary == UF_D_IRL_genus$MP_UFDEstgenus)/nrow(UF_D_IRL_genus)


TotalPCoA <- ordinate(UFDestuaryGenusAgenus,"PCoA")
p = plot_ordination(UFDestuaryGenusAgenus, TotalPCoA, color="Estuary") + 
  ggtitle("Unfiltered DESeq2 Estuary genus Indicators PCoA with Bray.Curtis distance") +
  geom_text(aes(label=UF_D_Est_genus$MP_UFDEstgenus),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p


#Test filtered DESeq2 Estuary genuss#

#Remove any samples with no genuss from target phyloseq#
nsamples(FDestuaryGenusAgenus) #90
FDestuaryGenusAgenus = prune_samples(sample_sums(FDestuaryGenusAgenus)>=1, FDestuaryGenusAgenus)
nsamples(FDestuaryGenusAgenus) #90

#Make a matrix with the sample_data and relabel subcategories to make metadata.defined column#
F_D_Est_genus <- as.matrix(sample_data(FDestuaryGenusAgenus))
F_D_Est_genus[F_D_Est_genus ==  "SLE"] <- 2
F_D_Est_genus[F_D_Est_genus ==  "IRL"] <- 1
F_D_Est_genus<-as.data.frame(F_D_Est_genus)


#Create Bray.Curtis distance matrix from genus table#
bc <- vegdist(t(otu_table(FDestuaryGenusAgenus)))

#Conduct PAM clustering based upon 2 clusters#
pam.res <- pam(bc,2, diss=TRUE)

#Add a microbially.predicted column to matrix above for comparison to metadata.defined column#
F_D_Est_genus <- cbind(F_D_Est_genus, MP_FDEstgenus = pam.res$cluster)

#Test Total Estuary Efficiency# 0.6222222
sum(F_D_Est_genus$Estuary == F_D_Est_genus$MP_FDEstgenus)/nrow(F_D_Est_genus)
#Test SLE Efficiency aka Sensitivity# 0.8333333
F_D_SLE_genus <- filter(F_D_Est_genus, Estuary == "2")
sum(F_D_SLE_genus$Estuary == F_D_SLE_genus$MP_FDEstgenus)/nrow(F_D_SLE_genus)
#Test IRL Efficiency aka Specificity# 0.5454545
F_D_IRL_genus <- filter(F_D_Est_genus, Estuary == "1")
sum(F_D_IRL_genus$Estuary == F_D_IRL_genus$MP_FDEstgenus)/nrow(F_D_IRL_genus)

TotalPCoA <- ordinate(FDestuaryGenusAgenus,"PCoA")
p = plot_ordination(FDestuaryGenusAgenus, TotalPCoA, color="Estuary") + 
  ggtitle("Filtered DESeq2 Estuary genus Indicators PCoA with Bray.Curtis distance") +
  geom_text(aes(label=F_D_Est_genus$MP_FDEstgenus),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p

###TEST COMBO ESTUARY INDICATORS###

##Test unfiltered Estuary genuss##

#Remove any samples with no genuss from target phyloseq#
nsamples(UFCestuaryGenusAgenus) #90
UFCestuaryGenusAgenus = prune_samples(sample_sums(UFCestuaryGenusAgenus)>=1, UFCestuaryGenusAgenus)
nsamples(UFCestuaryGenusAgenus) #90

#Make a matrix with the sample_data and relabel subcategories to make metadata.defined column#
UF_C_Est_genus <- as.matrix(sample_data(UFCestuaryGenusAgenus))
UF_C_Est_genus[UF_C_Est_genus ==  "SLE"] <- 2
UF_C_Est_genus[UF_C_Est_genus ==  "IRL"] <- 1
UF_C_Est_genus<-as.data.frame(UF_C_Est_genus)


#Create Bray.Curtis distance matrix from genus table#
bc <- vegdist(t(otu_table(UFCestuaryGenusAgenus)))

#Conduct PAM clustering based upon 2 clusters#
pam.res <- pam(bc,2, diss=TRUE)

#Add a microbially.predicted column to matrix above for comparison to metadata.defined column#
UF_C_Est_genus <- cbind(UF_C_Est_genus, MP_UFCEstgenus = pam.res$cluster)

#Test Total Estuary Efficiency# 1
sum(UF_C_Est_genus$Estuary == UF_C_Est_genus$MP_UFCEstgenus)/nrow(UF_C_Est_genus)

#Test SLE Efficiency aka Sensitivity# 1
UF_C_SLE_genus <- filter(UF_C_Est_genus, Estuary == "2")
sum(UF_C_SLE_genus$Estuary == UF_C_SLE_genus$MP_UFCEstgenus)/nrow(UF_C_SLE_genus)

#Test IRL Efficiency aka Specificity# 1
UF_C_IRL_genus <- filter(UF_C_Est_genus, Estuary == "1")
sum(UF_C_IRL_genus$Estuary == UF_C_IRL_genus$MP_UFCEstgenus)/nrow(UF_C_IRL_genus)

TotalPCoA <- ordinate(UFCestuaryGenusAgenus,"PCoA")
p = plot_ordination(UFCestuaryGenusAgenus, TotalPCoA, color="Estuary") + 
  ggtitle("Unfiltered Final Estuary genus Indicators PCoA with Bray.Curtis distance") +
  geom_text(aes(label=UF_C_Est_genus$MP_UFCEstgenus),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p


##Test filtered Estuary genuss##

#Remove any samples with no genuss from target phyloseq#
nsamples(FCestuaryGenusAgenus) #90
FCestuaryGenusAgenus = prune_samples(sample_sums(FCestuaryGenusAgenus)>=1, FCestuaryGenusAgenus)
nsamples(FCestuaryGenusAgenus) #90

#Make a matrix with the sample_data and relabel subcategories to make metadata.defined column#
F_C_Est_genus <- as.matrix(sample_data(FCestuaryGenusAgenus))
F_C_Est_genus[F_C_Est_genus ==  "SLE"] <- 2
F_C_Est_genus[F_C_Est_genus ==  "IRL"] <- 1
F_C_Est_genus<-as.data.frame(F_C_Est_genus)


#Create Bray.Curtis distance matrix from genus table#
bc <- vegdist(t(otu_table(FCestuaryGenusAgenus)))

#Conduct PAM clustering based upon 2 clusters#
pam.res <- pam(bc,2, diss=TRUE)

#Add a microbially.predicted column to matrix above for comparison to metadata.defined column#
F_C_Est_genus <- cbind(F_C_Est_genus, MP_FCEstgenus = pam.res$cluster)

#Test Total Estuary Efficiency# 0.9666667
sum(F_C_Est_genus$Estuary == F_C_Est_genus$MP_FCEstgenus)/nrow(F_C_Est_genus)

#Test SLE Efficiency aka Sensitivity# 0.875
F_C_SLE_genus <- filter(F_C_Est_genus, Estuary == "2")
sum(F_C_SLE_genus$Estuary == F_C_SLE_genus$MP_FCEstgenus)/nrow(F_C_SLE_genus)

#Test IRL Efficiency aka Sensitivity# 1
F_C_IRL_genus <- filter(F_C_Est_genus, Estuary == "1")
sum(F_C_IRL_genus$Estuary == F_C_IRL_genus$MP_FCEstgenus)/nrow(F_C_IRL_genus)

TotalPCoA <- ordinate(FCestuaryGenusAgenus,"PCoA")
p = plot_ordination(FCestuaryGenusAgenus, TotalPCoA, color="Estuary") + 
  ggtitle("Filtered Final Estuary genus Indicators PCoA with Bray.Curtis distance") +
  geom_text(aes(label=F_C_Est_genus$MP_FCEstgenus),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p


###TEST INDICSPECIES MUCK INDICATORS###

##Test unfiltered Muck genuss##

#Remove any samples with no genuss from target phyloseq#
nsamples(UFISmuckGenusAgenus) #90
UFISmuckGenusAgenus = prune_samples(sample_sums(UFISmuckGenusAgenus)>=1, UFISmuckGenusAgenus)
nsamples(UFISmuckGenusAgenus) #90

#Make a matrix with the sample_data and relabel subcategories to make metadata.defined column#
UF_IS_Muck_genus <- as.matrix(sample_data(UFISmuckGenusAgenus))
UF_IS_Muck_genus[UF_IS_Muck_genus ==  "Muck"] <- 2
UF_IS_Muck_genus[UF_IS_Muck_genus ==  "Not"] <- 1
UF_IS_Muck_genus[UF_IS_Muck_genus ==  "Mucky"] <- 2
UF_IS_Muck_genus[UF_IS_Muck_genus ==  "Muckish"] <- 2
UF_IS_Muck_genus<-as.data.frame(UF_IS_Muck_genus)


#Create Bray.Curtis distance matrix from genus table#
bc <- vegdist(t(otu_table(UFISmuckGenusAgenus)))

#Conduct PAM clustering based upon 2 clusters#
pam.res <- pam(bc,2, diss=TRUE)

#Add a microbially.predicted column to matrix above for comparison to metadata.defined column#
UF_IS_Muck_genus <- cbind(UF_IS_Muck_genus, MP_UFISMuckgenus = pam.res$cluster)

#Test Total Muck Efficiency# 0.8555556
sum(UF_IS_Muck_genus$Muck == UF_IS_Muck_genus$MP_UFISMuckgenus)/nrow(UF_IS_Muck_genus)

#Test 3MC Efficiency aka Sensitivity# 0.9333333
UF_IS_3MC_genus <- filter(UF_IS_Muck_genus, Muck == "2")
sum(UF_IS_3MC_genus$Muck == UF_IS_3MC_genus$MP_UFISMuckgenus)/nrow(UF_IS_3MC_genus)

#Test 0MC Efficiency aka Specificity# 0.8166667
UF_IS_0MC_genus <- filter(UF_IS_Muck_genus, Muck == "1")
sum(UF_IS_0MC_genus$Muck == UF_IS_0MC_genus$MP_UFISMuckgenus)/nrow(UF_IS_0MC_genus)


TotalPCoA <- ordinate(UFISmuckGenusAgenus,"PCoA")
p = plot_ordination(UFISmuckGenusAgenus, TotalPCoA, color="Muck") + 
  ggtitle("Unfiltered IS Muck genus Indicators PCoA with Bray.Curtis distance") +
  geom_text(aes(label=UF_IS_Muck_genus$MP_UFISMuckgenus),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p


##Test filtered Muck genuss##

#Remove any samples with no genuss from target phyloseq#
nsamples(FISmuckGenusAgenus) #90
FISmuckGenusAgenus = prune_samples(sample_sums(FISmuckGenusAgenus)>=1, FISmuckGenusAgenus)
nsamples(FISmuckGenusAgenus) #90

#Make a matrix with the sample_data and relabel subcategories to make metadata.defined column#
F_IS_Muck_genus <- as.matrix(sample_data(FISmuckGenusAgenus))
F_IS_Muck_genus[F_IS_Muck_genus ==  "Muck"] <- 2
F_IS_Muck_genus[F_IS_Muck_genus ==  "Not"] <- 1
F_IS_Muck_genus[F_IS_Muck_genus ==  "Mucky"] <- 2
F_IS_Muck_genus[F_IS_Muck_genus ==  "Muckish"] <- 2
F_IS_Muck_genus<-as.data.frame(F_IS_Muck_genus)


#Create Bray.Curtis distance matrix from genus table#
bc <- vegdist(t(otu_table(FISmuckGenusAgenus)))

#Conduct PAM clustering based upon 2 clusters#
pam.res <- pam(bc,2, diss=TRUE)

#Add a microbially.predicted column to matrix above for comparison to metadata.defined column#
F_IS_Muck_genus <- cbind(F_IS_Muck_genus, MP_FISMuckgenus = pam.res$cluster)

#Test Total Muck Efficiency# 0.9111111
sum(F_IS_Muck_genus$Muck == F_IS_Muck_genus$MP_FISMuckgenus)/nrow(F_IS_Muck_genus)
#Test 3MC Efficiency aka  Sensitivity# 1
F_IS_3MC_genus <- filter(F_IS_Muck_genus, Muck == "2")
sum(F_IS_3MC_genus$Muck == F_IS_3MC_genus$MP_FISMuckgenus)/nrow(F_IS_3MC_genus)
#Test 0MC Efficiency aka Specificity# 0.8666667
F_IS_0MC_genus <- filter(F_IS_Muck_genus, Muck == "1")
sum(F_IS_0MC_genus$Muck == F_IS_0MC_genus$MP_FISMuckgenus)/nrow(F_IS_0MC_genus)

TotalPCoA <- ordinate(FISmuckGenusAgenus,"PCoA")
p = plot_ordination(FISmuckGenusAgenus, TotalPCoA, color="Muck") + 
  ggtitle("Filtered IS Muck genus Indicators PCoA with Bray.Curtis distance") +
  geom_text(aes(label=F_IS_Muck_genus$MP_FISMuckgenus),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p

###TEST DESEQ2 MUCK INDICATORS###

#Test unfiltered Muck genuss##

#Remove any samples with no genuss from target phyloseq#
nsamples(UFDmuckGenusAgenus) #90
UFDmuckGenusAgenus = prune_samples(sample_sums(UFDmuckGenusAgenus)>=1, UFDmuckGenusAgenus)
nsamples(UFDmuckGenusAgenus) #90

#Make a matrix with the sample_data and relabel subcategories to make metadata.defined column#
UF_D_Muck_genus <- as.matrix(sample_data(UFDmuckGenusAgenus))
UF_D_Muck_genus[UF_D_Muck_genus ==  "Muck"] <- 2
UF_D_Muck_genus[UF_D_Muck_genus ==  "Not"] <- 1
UF_D_Muck_genus[UF_D_Muck_genus ==  "Mucky"] <- 2
UF_D_Muck_genus[UF_D_Muck_genus ==  "Muckish"] <- 2
UF_D_Muck_genus<-as.data.frame(UF_D_Muck_genus)


#Create Bray.Curtis distance matrix from genus table#
bc <- vegdist(t(otu_table(UFDmuckGenusAgenus)))

#Conduct PAM clustering based upon 2 clusters#
pam.res <- pam(bc,2, diss=TRUE)

#Add a microbially.predicted column to matrix above for comparison to metadata.defined column#
UF_D_Muck_genus <- cbind(UF_D_Muck_genus, MP_UFDMuckgenus = pam.res$cluster)

#Test Total Muck Efficiency# 0.8555556
sum(UF_D_Muck_genus$Muck == UF_D_Muck_genus$MP_UFDMuckgenus)/nrow(UF_D_Muck_genus)

#Test 3MC Efficiency aka Sensitivity# 0.7666667
UF_D_3MC_genus <- filter(UF_D_Muck_genus, Muck == "2")
sum(UF_D_3MC_genus$Muck == UF_D_3MC_genus$MP_UFDMuckgenus)/nrow(UF_D_3MC_genus)

#Test 0MC Efficiency aka Sensitivity# 0.9
UF_D_0MC_genus <- filter(UF_D_Muck_genus, Muck == "1")
sum(UF_D_0MC_genus$Muck == UF_D_0MC_genus$MP_UFDMuckgenus)/nrow(UF_D_0MC_genus)

TotalPCoA <- ordinate(UFDmuckGenusAgenus,"PCoA")
p = plot_ordination(UFDmuckGenusAgenus, TotalPCoA, color="Muck") + 
  ggtitle("Unfiltered DESeq2 Muck genus Indicators PCoA with Bray.Curtis distance") +
  geom_text(aes(label=UF_D_Muck_genus$MP_UFDMuckgenus),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p


##Test filtered Muck genuss##

#Remove any samples with no genuss from target phyloseq#
nsamples(FDmuckGenusAgenus) #90
FDmuckGenusAgenus = prune_samples(sample_sums(FDmuckGenusAgenus)>=1, FDmuckGenusAgenus)
nsamples(FDmuckGenusAgenus) #90

#Make a matrix with the sample_data and relabel subcategories to make metadata.defined column#
F_D_Muck_genus <- as.matrix(sample_data(FDmuckGenusAgenus))
F_D_Muck_genus[F_D_Muck_genus ==  "Muck"] <- 1
F_D_Muck_genus[F_D_Muck_genus ==  "Not"] <- 2
F_D_Muck_genus[F_D_Muck_genus ==  "Mucky"] <- 1
F_D_Muck_genus[F_D_Muck_genus ==  "Muckish"] <- 1
F_D_Muck_genus<-as.data.frame(F_D_Muck_genus)


#Create Bray.Curtis distance matrix from genus table#
bc <- vegdist(t(otu_table(FDmuckGenusAgenus)))

#Conduct PAM clustering based upon 2 clusters#
pam.res <- pam(bc,2, diss=TRUE)

#Add a microbially.predicted column to matrix above for comparison to metadata.defined column#
F_D_Muck_genus <- cbind(F_D_Muck_genus, MP_FDMuckgenus = pam.res$cluster)

#Test Total Muck Efficiency# 0.6777778
sum(F_D_Muck_genus$Muck == F_D_Muck_genus$MP_FDMuckgenus)/nrow(F_D_Muck_genus)

#Test 3MC Efficiency aka Sensitivy# 1
F_D_3MC_genus <- filter(F_D_Muck_genus, Muck == "1")
sum(F_D_3MC_genus$Muck == F_D_3MC_genus$MP_FDMuckgenus)/nrow(F_D_3MC_genus)

#Test 0MC Efficiency aka Specificity# 0.5166667
F_D_0MC_genus <- filter(F_D_Muck_genus, Muck == "2")
sum(F_D_0MC_genus$Muck == F_D_0MC_genus$MP_FDMuckgenus)/nrow(F_D_0MC_genus)

TotalPCoA <- ordinate(FDmuckGenusAgenus,"PCoA")
p = plot_ordination(FDmuckGenusAgenus, TotalPCoA, color="Muck") + 
  ggtitle("Filtered DESeq2 Muck genus Indicators PCoA with Bray.Curtis distance") +
  geom_text(aes(label=F_D_Muck_genus$MP_FDMuckgenus),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p

###TEST COMBO MUCK INDICATORS###

##Test unfiltered Muck genuss##

#Remove any samples with no genuss from target phyloseq#
nsamples(UFCmuckGenusAgenus) #90
UFCmuckGenusAgenus = prune_samples(sample_sums(UFCmuckGenusAgenus)>=1, UFCmuckGenusAgenus)
nsamples(UFCmuckGenusAgenus) #90

#Make a matrix with the sample_data and relabel subcategories to make metadata.defined column#
UF_C_Muck_genuss <- as.matrix(sample_data(UFCmuckGenusAgenus))
UF_C_Muck_genuss[UF_C_Muck_genuss ==  "Muck"] <- 2
UF_C_Muck_genuss[UF_C_Muck_genuss ==  "Not"] <- 1
UF_C_Muck_genuss[UF_C_Muck_genuss ==  "Mucky"] <- 2
UF_C_Muck_genuss[UF_C_Muck_genuss ==  "Muckish"] <- 2
UF_C_Muck_genuss<-as.data.frame(UF_C_Muck_genuss)


#Create Bray.Curtis distance matrix from genus table#
bc <- vegdist(t(otu_table(UFCmuckGenusAgenus)))

#Conduct PAM clustering based upon 2 clusters#
pam.res <- pam(bc,2, diss=TRUE)

#Add a microbially.predicted column to matrix above for comparison to metadata.defined column#
UF_C_Muck_genuss <- cbind(UF_C_Muck_genuss, MP_UFCMuckgenus = pam.res$cluster)

#Test Total Muck Efficiency# 0.8666667
sum(UF_C_Muck_genuss$Muck == UF_C_Muck_genuss$MP_UFCMuckgenus)/nrow(UF_C_Muck_genuss)

#Test Muck Efficiency aka Sensitivity# 0.9333333
F_C_3MC_genuss <- filter(UF_C_Muck_genuss, Muck == "2")
sum(F_C_3MC_genuss$Muck == F_C_3MC_genuss$MP_UFCMuckgenus)/nrow(F_C_3MC_genuss)

#Test Muck Efficiency aka Specificity# 0.8333333
F_C_0MC_genuss <- filter(UF_C_Muck_genuss, Muck == "1")
sum(F_C_0MC_genuss$Muck == F_C_0MC_genuss$MP_UFCMuckgenus)/nrow(F_C_0MC_genuss)

TotalPCoA <- ordinate(UFCmuckGenusAgenus,"PCoA")
p = plot_ordination(UFCmuckGenusAgenus, TotalPCoA, color="Muck") + 
  ggtitle("Unfiltered Final Muck genus Indicators PCoA with Bray.Curtis distance") +
  geom_text(aes(label=UF_C_Muck_genuss$MP_UFCMuckgenus),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p

##Test filtered Muck genuss##

#Remove any samples with no genuss from target phyloseq#
nsamples(FCmuckGenusAgenus) #90
FCmuckGenusAgenus = prune_samples(sample_sums(FCmuckGenusAgenus)>=1, FCmuckGenusAgenus)
nsamples(FCmuckGenusAgenus) #90

#Make a matrix with the sample_data and relabel subcategories to make metadata.defined column#
F_C_Muck_genuss <- as.matrix(sample_data(FCmuckGenusAgenus))
F_C_Muck_genuss[F_C_Muck_genuss ==  "Muck"] <- 2
F_C_Muck_genuss[F_C_Muck_genuss ==  "Not"] <- 1
F_C_Muck_genuss[F_C_Muck_genuss ==  "Mucky"] <- 2
F_C_Muck_genuss[F_C_Muck_genuss ==  "Muckish"] <- 2
F_C_Muck_genuss<-as.data.frame(F_C_Muck_genuss)


#Create Bray.Curtis distance matrix from genus table#
bc <- vegdist(t(otu_table(FCmuckGenusAgenus)))

#Conduct PAM clustering based upon 2 clusters#
pam.res <- pam(bc,2, diss=TRUE)

#Add a microbially.predicted column to matrix above for comparison to metadata.defined column#
F_C_Muck_genuss <- cbind(F_C_Muck_genuss, MP_FCMuckgenus = pam.res$cluster)

#Test Total Muck Efficiency# 0.8222222
sum(F_C_Muck_genuss$Muck == F_C_Muck_genuss$MP_FCMuckgenus)/nrow(F_C_Muck_genuss)

#Test Muck Efficiency aka Sensitivity# 1
F_C_3MC_genus <- filter(F_C_Muck_genuss, Muck == "2")
sum(F_C_3MC_genus$Muck == F_C_3MC_genus$MP_FCMuckgenus)/nrow(F_C_3MC_genus)

#Test Not Efficiency aka Specificity# 0.7333333
F_C_0MC_genus <- filter(F_C_Muck_genuss, Muck == "1")
sum(F_C_0MC_genus$Muck == F_C_0MC_genus$MP_FCMuckgenus)/nrow(F_C_0MC_genus)

TotalPCoA <- ordinate(FCmuckGenusAgenus,"PCoA")
p = plot_ordination(FCmuckGenusAgenus, TotalPCoA, color="Muck") + 
  ggtitle("Filtered Final Muck genus Indicators PCoA with Bray.Curtis distance") +
  geom_text(aes(label=F_C_Muck_genuss$MP_FCMuckgenus),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p

###TEST INDICSPECIES TOM/Cu INDICATORS###

##Test unfiltered TOM/Cu genuss##

#Create phyloseq object focused on HiHi and HiLo samples#
UFIShihiGenusAgenus = subset_samples(UFIStomcuGenusAgenus, LOI.Cu=="High.High")
UFIShiloGenusAgenus = subset_samples(UFIStomcuGenusAgenus, LOI.Cu=="High.Low")
UFISfoctomcuGenusAgenus = merge_phyloseq(UFIShihiGenusAgenus, UFIShiloGenusAgenus)

#Remove any samples with no genuss from target phyloseq#
nsamples(UFISfoctomcuGenusAgenus) #29
UFISfoctomcuGenusAgenus = prune_samples(sample_sums(UFISfoctomcuGenusAgenus)>=1, UFISfoctomcuGenusAgenus)
nsamples(UFISfoctomcuGenusAgenus) #29

#Make a matrix with the sample_data and relabel subcategories to make metadata.defined column#
UF_IS_TOMCu_genus <- as.matrix(sample_data(UFISfoctomcuGenusAgenus))
UF_IS_TOMCu_genus[UF_IS_TOMCu_genus ==  "High.High"] <- 1
UF_IS_TOMCu_genus[UF_IS_TOMCu_genus ==  "High.Low"] <- 2
UF_IS_TOMCu_genus<-as.data.frame(UF_IS_TOMCu_genus)

#Create Bray.Curtis distance matrix from genus table#
bc <- vegdist(t(otu_table(UFISfoctomcuGenusAgenus)))

#Conduct PAM clustering based upon 2 clusters#
pam.res <- pam(bc,2, diss=TRUE)

#Add a microbially.predicted column to matrix above for comparison to metadata.defined column#
UF_IS_TOMCu_genus <- cbind(UF_IS_TOMCu_genus, MP_UFISTOMCugenus = pam.res$cluster)

#Test Total TOM/Cu Efficiency# 0.8275862
sum(UF_IS_TOMCu_genus$LOI.Cu == UF_IS_TOMCu_genus$MP_UFISTOMCugenus)/nrow(UF_IS_TOMCu_genus)

#Test High.High Indicators aka Sensitivity# 1
UF_IS_HiHi_genus <- filter(UF_IS_TOMCu_genus, LOI.Cu == "1")
sum(UF_IS_HiHi_genus$LOI.Cu == UF_IS_HiHi_genus$MP_UFISTOMCugenus)/nrow(UF_IS_HiHi_genus)

#Test High.Low Indicators aka Specificity# 0.7222222
UF_IS_HiLo_genus <- filter(UF_IS_TOMCu_genus, LOI.Cu == "2")
sum(UF_IS_HiLo_genus$LOI.Cu == UF_IS_HiLo_genus$MP_UFISTOMCugenus)/nrow(UF_IS_HiLo_genus)

TotalPCoA <- ordinate(UFISfoctomcuGenusAgenus,"PCoA")
p = plot_ordination(UFISfoctomcuGenusAgenus, TotalPCoA, color="LOI.Cu") + 
  ggtitle("Unfiltered IS LOICu genus Indicators PCoA with Bray.Curtis distance") +
  geom_text(aes(label=UF_IS_TOMCu_genus$MP_UFISTOMCugenus),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p


##Test filtered TOM/Cu genuss##

#Create phyloseq object focused on HiHi and HiLo samples#
FIShihiGenusAgenus = subset_samples(FIStomcuGenusAgenus, LOI.Cu=="High.High")
FIShiloGenusAgenus = subset_samples(FIStomcuGenusAgenus, LOI.Cu=="High.Low")
FISfoctomcuGenusAgenus = merge_phyloseq(FIShihiGenusAgenus, FIShiloGenusAgenus)

#Remove any samples with no genuss from target phyloseq#
nsamples(FISfoctomcuGenusAgenus) #29
FISfoctomcuGenusAgenus = prune_samples(sample_sums(FISfoctomcuGenusAgenus)>=1, FISfoctomcuGenusAgenus)
nsamples(FISfoctomcuGenusAgenus) #29

#Make a matrix with the sample_data and relabel subcategories to make metadata.defined column#
F_IS_TOMCu_genus <- as.matrix(sample_data(FISfoctomcuGenusAgenus))
F_IS_TOMCu_genus[F_IS_TOMCu_genus ==  "High.High"] <- 1
F_IS_TOMCu_genus[F_IS_TOMCu_genus ==  "High.Low"] <- 2
F_IS_TOMCu_genus<-as.data.frame(F_IS_TOMCu_genus)


#Create Bray.Curtis distance matrix from genus table#
bc <- vegdist(t(otu_table(FISfoctomcuGenusAgenus)))

#Conduct PAM clustering based upon 2 clusters#
pam.res <- pam(bc,2, diss=TRUE)

#Add a microbially.predicted column to matrix above for comparison to metadata.defined column#
F_IS_TOMCu_genus <- cbind(F_IS_TOMCu_genus, MP_FISTOMCugenus = pam.res$cluster)

#Test Total TOM/Cu Efficiency# 1
sum(F_IS_TOMCu_genus$LOI.Cu == F_IS_TOMCu_genus$MP_FISTOMCugenus)/nrow(F_IS_TOMCu_genus)

#Test HiHi Efficiency aka Sensitivity# 1
F_IS_HiHi_genus <- filter(F_IS_TOMCu_genus, LOI.Cu == "1")
sum(F_IS_HiHi_genus$LOI.Cu == F_IS_HiHi_genus$MP_FISTOMCugenus)/nrow(F_IS_HiHi_genus)

#Test HiLo Efficiency aka Sensitivity# 1
F_IS_HiLo_genus <- filter(F_IS_TOMCu_genus, LOI.Cu == "2") 
sum(F_IS_HiLo_genus$LOI.Cu == F_IS_HiLo_genus$MP_FISTOMCugenus)/nrow(F_IS_HiLo_genus)

TotalPCoA <- ordinate(FISfoctomcuGenusAgenus,"PCoA")
p = plot_ordination(FISfoctomcuGenusAgenus, TotalPCoA, color="LOI.Cu") + 
  ggtitle("Filtered IS LOICu genus Indicators PCoA with Bray.Curtis distance") +
  geom_text(aes(label=F_IS_TOMCu_genus$MP_FISTOMCugenus),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p

##TEST DESEQ2 LOICu INDICATORS##

##Test unfiltered TOM/Cu genuss##

#Create phyloseq object focused on HiHi and HiLo samples#
UFDhihiGenusAgenus = subset_samples(UFDtomcuGenusAgenus, LOI.Cu=="High.High")
UFDhiloGenusAgenus = subset_samples(UFDtomcuGenusAgenus, LOI.Cu=="High.Low")
UFDfoctomcuGenusAgenus = merge_phyloseq(UFDhihiGenusAgenus, UFDhiloGenusAgenus)

#Remove any samples with no genuss from target phyloseq#
nsamples(UFDfoctomcuGenusAgenus) #29
UFDfoctomcuGenusAgenus = prune_samples(sample_sums(UFDfoctomcuGenusAgenus)>=1, UFDfoctomcuGenusAgenus)
nsamples(UFDfoctomcuGenusAgenus) #29

#Make a matrix with the sample_data and relabel subcategories to make metadata.defined column#
UF_D_TOMCu_genus <- as.matrix(sample_data(UFDfoctomcuGenusAgenus))
UF_D_TOMCu_genus[UF_D_TOMCu_genus ==  "High.High"] <- 1
UF_D_TOMCu_genus[UF_D_TOMCu_genus ==  "High.Low"] <- 2
UF_D_TOMCu_genus<-as.data.frame(UF_D_TOMCu_genus)


#Create Bray.Curtis distance matrix from genus table#
bc <- vegdist(t(otu_table(UFDfoctomcuGenusAgenus)))

#Conduct PAM clustering based upon 2 clusters#
pam.res <- pam(bc,2, diss=TRUE)

#Add a microbially.predicted column to matrix above for comparison to metadata.defined column#
UF_D_TOMCu_genus <- cbind(UF_D_TOMCu_genus, MP_UFDTOMCugenus = pam.res$cluster)

#Test Total TOM/Cu Efficiency# 0.8275862
sum(UF_D_TOMCu_genus$LOI.Cu == UF_D_TOMCu_genus$MP_UFDTOMCugenus)/nrow(UF_D_TOMCu_genus)

#Test HiHi Efficiency aka Sensitivity# 1
UF_D_HiHi_genus <- filter(UF_D_TOMCu_genus, LOI.Cu == "1")
sum(UF_D_HiHi_genus$LOI.Cu == UF_D_HiHi_genus$MP_UFDTOMCugenus)/nrow(UF_D_HiHi_genus)

#Test HiLo Efficiency aka Specificity# 0.7222222
UF_D_HiLo_genus <- filter(UF_D_TOMCu_genus, LOI.Cu == "2")
sum(UF_D_HiLo_genus$LOI.Cu == UF_D_HiLo_genus$MP_UFDTOMCugenus)/nrow(UF_D_HiLo_genus)

TotalPCoA <- ordinate(UFDfoctomcuGenusAgenus,"PCoA")
p = plot_ordination(UFDfoctomcuGenusAgenus, TotalPCoA, color="LOI.Cu") + 
  ggtitle("Unfiltered DESeq2 LOICu genus Indicators PCoA with Bray.Curtis distance") +
  geom_text(aes(label=UF_D_TOMCu_genus$MP_UFDTOMCugenus),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p


##Test filtered TOM/Cu genuss##

#Create phyloseq object focused on HiHi and HiLo samples#
FDhihiGenusAgenus = subset_samples(FDtomcuGenusAgenus, LOI.Cu=="High.High")
FDhiloGenusAgenus = subset_samples(FDtomcuGenusAgenus, LOI.Cu=="High.Low")
FDfoctomcuGenusAgenus = merge_phyloseq(FDhihiGenusAgenus, FDhiloGenusAgenus)

#Remove any samples with no genuss from target phyloseq#
nsamples(FDfoctomcuGenusAgenus) #29
FDfoctomcuGenusAgenus = prune_samples(sample_sums(FDfoctomcuGenusAgenus)>=1, FDfoctomcuGenusAgenus)
nsamples(FDfoctomcuGenusAgenus) #29


#Make a matrix with the sample_data and relabel subcategories to make metadata.defined column#
F_D_TOMCu_genus <- as.matrix(sample_data(FDfoctomcuGenusAgenus))
F_D_TOMCu_genus[F_D_TOMCu_genus ==  "High.High"] <- 1
F_D_TOMCu_genus[F_D_TOMCu_genus ==  "High.Low"] <- 2
F_D_TOMCu_genus<-as.data.frame(F_D_TOMCu_genus)


#Create Bray.Curtis distance matrix from genus table#
bc <- vegdist(t(otu_table(FDfoctomcuGenusAgenus)))

#Conduct PAM clustering based upon 2 clusters#
pam.res <- pam(bc,2, diss=TRUE)

#Add a microbially.predicted column to matrix above for comparison to metadata.defined column#
F_D_TOMCu_genus <- cbind(F_D_TOMCu_genus, MP_FDTOMCugenus = pam.res$cluster)

#Test Total TOM/Cu Efficiency# 0.8275862
sum(F_D_TOMCu_genus$LOI.Cu == F_D_TOMCu_genus$MP_FDTOMCugenus)/nrow(F_D_TOMCu_genus)

#Test HiHi Efficiency aka Sensitivity# 0.5454545
F_D_HiHi_genus <- filter(F_D_TOMCu_genus, LOI.Cu == "1")
sum(F_D_HiHi_genus$LOI.Cu == F_D_HiHi_genus$MP_FDTOMCugenus)/nrow(F_D_HiHi_genus)

#Test HiLo Efficiency aka Specificity# 1
F_D_HiLo_genus <- filter(F_D_TOMCu_genus, LOI.Cu == "2")
sum(F_D_HiLo_genus$LOI.Cu == F_D_HiLo_genus$MP_FDTOMCugenus)/nrow(F_D_HiLo_genus)

TotalPCoA <- ordinate(FDfoctomcuGenusAgenus,"PCoA")
p = plot_ordination(FDfoctomcuGenusAgenus, TotalPCoA, color="LOI.Cu") + 
  ggtitle("Filtered DESeq2 LOICu genus Indicators PCoA with Bray.Curtis distance") +
  geom_text(aes(label=F_D_TOMCu_genus$MP_FDTOMCugenus),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p

##TEST COMBO TOM/CU INDICATORS##

##Test unfiltered TOM/Cu genuss##

#Create phyloseq object focused on HiHi and HiLo samples#
UFChihiGenusAgenus = subset_samples(UFCtomcuGenusAgenus, LOI.Cu=="High.High")
UFChiloGenusAgenus = subset_samples(UFCtomcuGenusAgenus, LOI.Cu=="High.Low")
UFCfoctomcuGenusAgenus = merge_phyloseq(UFChihiGenusAgenus, UFChiloGenusAgenus)

#Remove any samples with no genuss from target phyloseq#
nsamples(UFCfoctomcuGenusAgenus) #29
UFCfoctomcuGenusAgenus = prune_samples(sample_sums(UFCfoctomcuGenusAgenus)>=1, UFCfoctomcuGenusAgenus)
nsamples(UFCfoctomcuGenusAgenus) #29

#Make a matrix with the sample_data and relabel subcategories to make metadata.defined column#
UF_C_TOMCu_genus <- as.matrix(sample_data(UFCfoctomcuGenusAgenus))
UF_C_TOMCu_genus[UF_C_TOMCu_genus ==  "High.High"] <- 1
UF_C_TOMCu_genus[UF_C_TOMCu_genus ==  "High.Low"] <- 2
UF_C_TOMCu_genus<-as.data.frame(UF_C_TOMCu_genus)


#Create Bray.Curtis distance matrix from genus table#
bc <- vegdist(t(otu_table(UFCfoctomcuGenusAgenus)))

#Conduct PAM clustering based upon 2 clusters#
pam.res <- pam(bc,2, diss=TRUE)

#Add a microbially.predicted column to matrix above for comparison to metadata.defined column#
UF_C_TOMCu_genus <- cbind(UF_C_TOMCu_genus, MP_UFCTOMCugenus = pam.res$cluster)

#Test Total TOM/Cu Efficiency# 0.8275862
sum(UF_C_TOMCu_genus$LOI.Cu == UF_C_TOMCu_genus$MP_UFCTOMCugenus)/nrow(UF_C_TOMCu_genus)

#Test HiHi Efficiency aka Sensitivity# 1
UF_C_HiHi_genus <- filter(UF_C_TOMCu_genus, LOI.Cu == "1")
sum(UF_C_HiHi_genus$LOI.Cu == UF_C_HiHi_genus$MP_UFCTOMCugenus)/nrow(UF_C_HiHi_genus)

#Test HiLo Efficiency aka Specificity# 0.7222222
UF_C_HiLo_genus <- filter(UF_C_TOMCu_genus, LOI.Cu == "2")
sum(UF_C_HiLo_genus$LOI.Cu == UF_C_HiLo_genus$MP_UFCTOMCugenus)/nrow(UF_C_HiLo_genus)

TotalPCoA <- ordinate(UFCfoctomcuGenusAgenus,"PCoA")
p = plot_ordination(UFCfoctomcuGenusAgenus, TotalPCoA, color="LOI.Cu") + 
  ggtitle("Unfiltered Final LOICu genus Indicators PCoA with Bray.Curtis distance") +
  geom_text(aes(label=UF_C_TOMCu_genus$MP_UFCTOMCugenus),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p


##Test filtered TOM/Cu genuss##

#Create phyloseq object focused on HiHi and HiLo samples#
FChihiGenusAgenus = subset_samples(FCtomcuGenusAgenus, LOI.Cu=="High.High")
FChiloGenusAgenus = subset_samples(FCtomcuGenusAgenus, LOI.Cu=="High.Low")
FCfoctomcuGenusAgenus = merge_phyloseq(FChihiGenusAgenus, FChiloGenusAgenus)

#Remove any samples with no genuss from target phyloseq#
nsamples(FCfoctomcuGenusAgenus) #29
FCfoctomcuGenusAgenus = prune_samples(sample_sums(FCfoctomcuGenusAgenus)>=1, FCfoctomcuGenusAgenus)
nsamples(FCfoctomcuGenusAgenus) #29

#Make a matrix with the sample_data and relabel subcategories to make metadata.defined column#
F_C_TOMCu_genus <- as.matrix(sample_data(FCfoctomcuGenusAgenus))
F_C_TOMCu_genus[F_C_TOMCu_genus ==  "High.High"] <- 1
F_C_TOMCu_genus[F_C_TOMCu_genus ==  "High.Low"] <- 2
F_C_TOMCu_genus<-as.data.frame(F_C_TOMCu_genus)


#Create Bray.Curtis distance matrix from genus table#
bc <- vegdist(t(otu_table(FCfoctomcuGenusAgenus)))

#Conduct PAM clustering based upon 2 clusters#
pam.res <- pam(bc,2, diss=TRUE)

#Add a microbially.predicted column to matrix above for comparison to metadata.defined column#
F_C_TOMCu_genus <- cbind(F_C_TOMCu_genus, MP_FCTOMCugenus = pam.res$cluster)

#Test Total TOM/Cu Efficiency# 1
sum(F_C_TOMCu_genus$LOI.Cu == F_C_TOMCu_genus$MP_FCTOMCugenus)/nrow(F_C_TOMCu_genus)

#Test HiHi Efficiency aka Sensitivity# 1
F_C_HiHi_genus <- filter(F_C_TOMCu_genus, LOI.Cu == "1")
sum(F_C_HiHi_genus$LOI.Cu == F_C_HiHi_genus$MP_FCTOMCugenus)/nrow(F_C_HiHi_genus)

#Test HiLo Efficiency aka Specificity# 1
F_C_HiLo_genus <- filter(F_C_TOMCu_genus, LOI.Cu == "2")
sum(F_C_HiLo_genus$LOI.Cu == F_C_HiLo_genus$MP_FCTOMCugenus)/nrow(F_C_HiLo_genus)


TotalPCoA <- ordinate(FCfoctomcuGenusAgenus,"PCoA")
p = plot_ordination(FCfoctomcuGenusAgenus, TotalPCoA, color="LOI.Cu") + 
  ggtitle("Filtered Final LOICu genus Indicators PCoA with Bray.Curtis distance") +
  geom_text(aes(label=F_C_TOMCu_genus$MP_FCTOMCugenus),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p

##FAMILY##
getwd() #"C:/Users/bradshawd/Documents/Journal_Articles/New_Indicators"
setwd("C:/Users/bradshawd/Documents/Journal_Articles/New_Indicators/Family")
getwd() #"C:/Users/bradshawd/Documents/Journal_Articles/New_Indicators/Family"

##Upload QIIME2 information into phyloseq for taxnomic level you want to test##
family_aa_table = read.csv(file.choose(), row.names=1)
family_taxa_table = read.csv(file.choose(), row.names = 1)
family_taxa_table = as.matrix(family_taxa_table)

#Can use the same metadata table as species level#
Family_OTU = otu_table(family_aa_table, taxa_are_rows = TRUE)
Family_TAX = tax_table(family_taxa_table)

FamilyA <- merge_phyloseq(Family_OTU, Family_TAX, dotted_meta)

colnames(tax_table(FamilyA))

###Create at more filtered set of LWSS Data###
ntaxa(FamilyA) #1557#
nsamples(FamilyA)#483

LWSFamilyA <- subset_samples(FamilyA, Survey=="Lagoon Wide")
nsamples(LWSFamilyA)#324

LWSWFamilyA <- subset_samples(LWSFamilyA, Medium=="Water")
nsamples(LWSWFamilyA)#132, 11 sites in triplicate (33) x 4 sampling periods (sps)

LWSSFamilyA1 <- subset_samples(FamilyA, Survey=="Lagoon Wide.Post Hurricane")
nsamples(LWSSFamilyA1)#12

LWSSFamilyA2 <- subset_samples(LWSFamilyA, Medium=="Sediment")
nsamples(LWSSFamilyA2)#192

LWSSFamilyA <- merge_phyloseq(LWSSFamilyA1, LWSSFamilyA2)
nsamples(LWSSFamilyA)#204

##Create phyloseqs for determining and testing indicators for Sediment
#Split samples by sampling period (SP)
LWSS_W16_FamilyA <- subset_samples(LWSSFamilyA, AbbSeasonAbbYear=="W16")
nsamples(LWSS_W16_FamilyA)#45, 15 sites in triplicate

LWSS_D17_FamilyA <- subset_samples(LWSSFamilyA, AbbSeasonAbbYear=="D17")
nsamples(LWSS_D17_FamilyA)#45,15 sites in triplicate

LWSS_W17_FamilyA <- subset_samples(LWSSFamilyA, AbbSeasonAbbYear=="W17")
nsamples(LWSS_W17_FamilyA)#57, 19 sites in triplicate

LWSS_D18_FamilyA <- subset_samples(LWSSFamilyA, AbbSeasonAbbYear=="D18")
nsamples(LWSS_D18_FamilyA)#57, 19 sites in triplicate

#Combine sampling periods into first and second years
LWSS_Y1_FamilyA <- merge_phyloseq(LWSS_W16_FamilyA, LWSS_D17_FamilyA)
nsamples(LWSS_Y1_FamilyA)#90

LWSS_Y2_FamilyA <- merge_phyloseq(LWSS_W17_FamilyA, LWSS_D18_FamilyA)
nsamples(LWSS_Y2_FamilyA)#114

##Create phyloseqs for determining and testing indicators for Water
#Split samples by sampling period (SP)
LWSW_W16_FamilyA <- subset_samples(LWSWFamilyA, AbbSeasonAbbYear=="W16")
nsamples(LWSW_W16_FamilyA)#33, 11 sites in triplicate

LWSW_D17_FamilyA <- subset_samples(LWSWFamilyA, AbbSeasonAbbYear=="D17")
nsamples(LWSW_D17_FamilyA)#33,11 sites in triplicate

LWSW_W17_FamilyA <- subset_samples(LWSWFamilyA, AbbSeasonAbbYear=="W17")
nsamples(LWSW_W17_FamilyA)#33, 11 sites in triplicate

LWSW_D18_FamilyA <- subset_samples(LWSWFamilyA, AbbSeasonAbbYear=="D18")
nsamples(LWSW_D18_FamilyA)#33, 11 sites in triplicate

#Combine sampling periods into first and second years
LWSW_Y1_FamilyA <- merge_phyloseq(LWSW_W16_FamilyA, LWSW_D17_FamilyA)
nsamples(LWSW_Y1_FamilyA)#66

LWSW_Y2_FamilyA <- merge_phyloseq(LWSW_W17_FamilyA, LWSW_D18_FamilyA)
nsamples(LWSW_Y2_FamilyA)#66

##Filter yearly sample sets of all zero and then low abundance familys
FLWSS_Y1_FamilyA <- filter_taxa(LWSS_Y1_FamilyA, function(x) sum(x) >0, TRUE)
ntaxa(FLWSS_Y1_FamilyA) #1232#
RFLWSS_Y1_FamilyAfamily <- filter_taxa(FLWSS_Y1_FamilyA, function(x) sum(x) >20, TRUE)
ntaxa(RFLWSS_Y1_FamilyAfamily) #776#

FLWSS_Y2_FamilyA <- filter_taxa(LWSS_Y2_FamilyA, function(x) sum(x) >0, TRUE)
ntaxa(FLWSS_Y2_FamilyA) #1268#
RFLWSS_Y2_FamilyAfamily <- filter_taxa(FLWSS_Y2_FamilyA, function(x) sum(x) >20, TRUE)
ntaxa(RFLWSS_Y2_FamilyAfamily) #698#

FLWSW_Y1_FamilyA <- filter_taxa(LWSW_Y1_FamilyA, function(x) sum(x) >0, TRUE)
ntaxa(FLWSW_Y1_FamilyA) #620#
RFLWSW_Y1_FamilyAfamily <- filter_taxa(FLWSW_Y1_FamilyA, function(x) sum(x) >20, TRUE)
ntaxa(RFLWSW_Y1_FamilyAfamily) #332#

FLWSW_Y2_FamilyA <- filter_taxa(LWSW_Y2_FamilyA, function(x) sum(x) >0, TRUE)
ntaxa(FLWSW_Y2_FamilyA) #900#
RFLWSW_Y2_FamilyAfamily <- filter_taxa(FLWSW_Y2_FamilyA, function(x) sum(x) >20, TRUE)
ntaxa(RFLWSW_Y2_FamilyAfamily) #444#





#Determine the number of samples in each category and subcategory#
length(which(sample_data(RFLWSS_Y1_FamilyAfamily)$Estuary == "SLE")) #24
length(which(sample_data(RFLWSS_Y1_FamilyAfamily)$Estuary == "IRL")) #66

length(which(sample_data(RFLWSS_Y1_FamilyAfamily)$Muck == "Not")) #60
length(which(sample_data(RFLWSS_Y1_FamilyAfamily)$Muck == "Muck")) #22
length(which(sample_data(RFLWSS_Y1_FamilyAfamily)$Muck == "Mucky")) #4
length(which(sample_data(RFLWSS_Y1_FamilyAfamily)$Muck == "Muckish")) #4

length(which(sample_data(RFLWSS_Y1_FamilyAfamily)$LOI.Cu == "Low.Low")) #61
length(which(sample_data(RFLWSS_Y1_FamilyAfamily)$LOI.Cu == "High.Low")) #18
length(which(sample_data(RFLWSS_Y1_FamilyAfamily)$LOI.Cu == "High.High")) #11
length(which(sample_data(RFLWSS_Y1_FamilyAfamily)$LOI.Cu == "Low.High")) #0

##NOTE##
#Use res to determine what is + or .#
#muck res: Not (+) vs Muck (.)#
#tomcu res: HL (+) vs HH (.)#
#estuary res: SLE (+) vs IRL (.)#

###RUN DESEQ2 AT family LEVEL ON ESTUARY SUBCATEGORIES###
#Convert phyloseq to deseq2 object centered around a irlvssle factor#
estuary_cudds_family = phyloseq_to_deseq2(RFLWSS_Y1_FamilyAfamily,  ~ Estuary)

#Calculate geometric means prior to estimate size factors#
estuary_geoMeans_family = apply(counts(estuary_cudds_family), 1, gm_mean)
estuary_cudds_family = estimateSizeFactors(estuary_cudds_family, geoMeans = estuary_geoMeans_family)

#Conduct DESEQ2 test#
estuary_cudds_family = DESeq(estuary_cudds_family, fitType="local")

#Explore the results#
estuary_DESeq2_res_family = results(estuary_cudds_family)
estuary_DESeq2_res_family = estuary_DESeq2_res_family[order(estuary_DESeq2_res_family$padj, na.last=NA), ]
alpha = 0.05
estuary_DESeq2_sig_res_family = estuary_DESeq2_res_family[(estuary_DESeq2_res_family$padj < alpha), ]

#Make dataframe with taxanomy added in#
estuary_DESeq2_sig_res_taxo_family = cbind(as(estuary_DESeq2_sig_res_family, "data.frame"), as(tax_table(RFLWSS_Y1_FamilyAfamily)[rownames(estuary_DESeq2_sig_res_family), ], "matrix"))

#Make dataframe with familys added in from all sites#
estuary_DESeq2_sig_res_taxo_seqs_family = cbind(as(estuary_DESeq2_sig_res_taxo_family, "data.frame"), as(otu_table(RFLWSS_Y1_FamilyAfamily)[rownames(estuary_DESeq2_sig_res_taxo_family), ], "matrix"))

#Make rownames an actual column and remove old rownames#
estuary_DESeq2_sig_res_taxo_seqs_family <- cbind(ESV.ID = rownames(estuary_DESeq2_sig_res_taxo_seqs_family), estuary_DESeq2_sig_res_taxo_seqs_family)
rownames(estuary_DESeq2_sig_res_taxo_seqs_family) <-NULL

#Make a txt file of results for all sites#
write.csv(as.data.frame(estuary_DESeq2_sig_res_taxo_seqs_family), file="estuary_DESeq2_sig_res_taxo_seqs_family.csv")

#Detemine which subcategory is negative or positive#
estuary_DESeq2_res_family
#estuary_DESeq2_res_family: SLE (+) vs IRL (.)#

#Determine number of Indicators with padj >0.05 in each subcagetory#
#Number of SLE indicators# 
length(which(estuary_DESeq2_sig_res_taxo_seqs_family$log2FoldChange > 0)) #140
#Number of IRL indicators#
length(which(estuary_DESeq2_sig_res_taxo_seqs_family$log2FoldChange < 0)) #177

###RUN DESEQ2 AT family LEVEL ON MUCK EXTREME SUBCATEGORIES###

##Prepare for running DESeq2##
#Make phyloseq focusing on Muck and Not subcategories#  
RFLWSS_Y1_Muck_FamilyAfamily = subset_samples(RFLWSS_Y1_FamilyAfamily, Muck=="Muck")
RFLWSS_Y1_Not_FamilyAfamily = subset_samples(RFLWSS_Y1_FamilyAfamily, Muck=="Not")
RFLWSS_Y1_MuckvsNot_FamilyAfamily = merge_phyloseq(RFLWSS_Y1_Muck_FamilyAfamily, RFLWSS_Y1_Not_FamilyAfamily)

##Run DESEQ2 and make a txt table of results##
#Convert phyloseq to deseq2 object centered around a muck factor#
muck_cudds_family = phyloseq_to_deseq2(RFLWSS_Y1_MuckvsNot_FamilyAfamily,  ~ Muck)

#Calculate geometric means prior to estimate size factors#
muck_geoMeans_family = apply(counts(muck_cudds_family), 1, gm_mean)
muck_cudds_family = estimateSizeFactors(muck_cudds_family, geoMeans = muck_geoMeans_family)

#Conduct DESEQ2 test#
muck_cudds_family = DESeq(muck_cudds_family, fitType="local")

#Explore the results#
muck_DESeq2_res_family = results(muck_cudds_family)
muck_DESeq2_res_family = muck_DESeq2_res_family[order(muck_DESeq2_res_family$padj, na.last=NA), ]
alpha = 0.05
muck_DESeq2_sig_res_family = muck_DESeq2_res_family[(muck_DESeq2_res_family$padj < alpha), ]

#Make dataframe with taxanomy added in#
muck_DESeq2_sig_res_taxo_family = cbind(as(muck_DESeq2_sig_res_family, "data.frame"), as(tax_table(RFLWSS_Y1_MuckvsNot_FamilyAfamily)[rownames(muck_DESeq2_sig_res_family), ], "matrix"))

#Make dataframe with familys added in from all sites#
muck_DESeq2_sig_res_taxo_seqs_family = cbind(as(muck_DESeq2_sig_res_taxo_family, "data.frame"), as(otu_table(RFLWSS_Y1_FamilyAfamily)[rownames(muck_DESeq2_sig_res_taxo_family), ], "matrix"))

#Make rownames an actual column and remove old rownames#
muck_DESeq2_sig_res_taxo_seqs_family <- cbind(ESV.ID = rownames(muck_DESeq2_sig_res_taxo_seqs_family), muck_DESeq2_sig_res_taxo_seqs_family)
rownames(muck_DESeq2_sig_res_taxo_seqs_family) <-NULL

#Make a txt file of results for all sites#
write.csv(as.data.frame(muck_DESeq2_sig_res_taxo_seqs_family), file="muck_DESeq2_sig_res_taxo_seqs_family.csv")

#Detemine which subcategory is negative or positive#
muck_DESeq2_res_family
#muck_DESeq2_res_family: Not (+) vs Muck (.)#
#Determine number of Indicators with padj >0.05 in each subcagetory#
#Number of Not indicators# 
length(which(muck_DESeq2_sig_res_taxo_seqs_family$log2FoldChange > 0)) #105
#Number of Muck indicators#
length(which(muck_DESeq2_sig_res_taxo_seqs_family$log2FoldChange < 0)) #73

###RUN DESEQ2 AT family LEVEL ON TOM/Cu EXTREME SUBCATEGORIES###

##Prepare for running DESeq2##
#Make phyloseq focusing on High.High and High.Low subcategories#  
RFLWSS_Y1_HH_FamilyAfamily = subset_samples(RFLWSS_Y1_FamilyAfamily, LOI.Cu=="High.High")
RFLWSS_Y1_HL_FamilyAfamily = subset_samples(RFLWSS_Y1_FamilyAfamily, LOI.Cu=="High.Low")
RFLWSS_Y1_HHvsHL_FamilyAfamily = merge_phyloseq(RFLWSS_Y1_HH_FamilyAfamily, RFLWSS_Y1_HL_FamilyAfamily)

##Run DESEQ2 and make a txt table of results##
#Convert phyloseq to deseq2 object centered around a loicu factor#
tomcu_cudds_family = phyloseq_to_deseq2(RFLWSS_Y1_HHvsHL_FamilyAfamily,  ~ LOI.Cu)

#Calculate geometric means prior to estimate size factors#
tomcu_geoMeans_family = apply(counts(tomcu_cudds_family), 1, gm_mean)
tomcu_cudds_family = estimateSizeFactors(tomcu_cudds_family, geoMeans = tomcu_geoMeans_family)

#Conduct DESEQ2 test#
tomcu_cudds_family = DESeq(tomcu_cudds_family, fitType="local")

#Explore the results#
tomcu_DESeq2_res_family = results(tomcu_cudds_family)
tomcu_DESeq2_res_family = tomcu_DESeq2_res_family[order(tomcu_DESeq2_res_family$padj, na.last=NA), ]
alpha = 0.05
tomcu_DESeq2_sig_res_family = tomcu_DESeq2_res_family[(tomcu_DESeq2_res_family$padj < alpha), ]

#Make dataframe with taxanomy added in#
tomcu_DESeq2_sig_res_taxo_family = cbind(as(tomcu_DESeq2_sig_res_family, "data.frame"), as(tax_table(RFLWSS_Y1_HHvsHL_FamilyAfamily)[rownames(tomcu_DESeq2_sig_res_family), ], "matrix"))

#Make dataframe with familys added in from all sites#
tomcu_DESeq2_sig_res_taxo_seqs_family = cbind(as(tomcu_DESeq2_sig_res_taxo_family, "data.frame"), as(otu_table(RFLWSS_Y1_FamilyAfamily)[rownames(tomcu_DESeq2_sig_res_taxo_family), ], "matrix"))

#Make rownames an actual column and remove old rownames#
tomcu_DESeq2_sig_res_taxo_seqs_family <- cbind(ESV.ID = rownames(tomcu_DESeq2_sig_res_taxo_seqs_family), tomcu_DESeq2_sig_res_taxo_seqs_family)
rownames(tomcu_DESeq2_sig_res_taxo_seqs_family) <-NULL

#Make a txt file of results for all sites#
write.csv(as.data.frame(tomcu_DESeq2_sig_res_taxo_seqs_family), file="tomcu_DESeq2_sig_res_taxo_seqs_family.csv")

#Detemine which subcategory is negative or positive#
tomcu_DESeq2_res_family
#tomcu_DESeq2_res_family: High.Low (+) vs High.High (.)#
#Determine number of Indicators with padj >0.05 in each subcagetory#
#Number of High.Low indicators# 
length(which(tomcu_DESeq2_sig_res_taxo_seqs_family$log2FoldChange > 0)) #93
#Number of High.High indicators#
length(which(tomcu_DESeq2_sig_res_taxo_seqs_family$log2FoldChange < 0)) #66

#Order rows by log2FoldChange#
muck_DESeq2_sig_res_taxo_seqs_family <- arrange(muck_DESeq2_sig_res_taxo_seqs_family, log2FoldChange)
tomcu_DESeq2_sig_res_taxo_seqs_family <- arrange(tomcu_DESeq2_sig_res_taxo_seqs_family, log2FoldChange)
estuary_DESeq2_sig_res_taxo_seqs_family <- arrange(estuary_DESeq2_sig_res_taxo_seqs_family, log2FoldChange)

#Create lists of just significant familys for each category#
unfiltered_DESeq2_estuary_familys <- subset(estuary_DESeq2_sig_res_taxo_seqs_family, select=c(ESV.ID))
unfiltered_DESeq2_muck_familys <- subset(muck_DESeq2_sig_res_taxo_seqs_family, select=c(ESV.ID))
unfiltered_DESeq2_tomcu_familys <- subset(tomcu_DESeq2_sig_res_taxo_seqs_family, select=c(ESV.ID))



##Create indicator tables filtered of overlapping indicators with other two metadata categories##

#Filter out rows in the IRLSLE DESEQ object that are also in the LOI.Cu and Muck DESEQs#
filtered_estuary_DESeq2_sig_res_taxo_seqs_family <- anti_join(estuary_DESeq2_sig_res_taxo_seqs_family, unfiltered_DESeq2_tomcu_familys, by = "ESV.ID")
filtered_estuary_DESeq2_sig_res_taxo_seqs_family <- anti_join(filtered_estuary_DESeq2_sig_res_taxo_seqs_family, unfiltered_DESeq2_muck_familys, by = "ESV.ID")

#Filter out rows in the Muck DESEQ object that are also in the IRL.SLE and LOI.Cu DESEQs#
filtered_muck_DESeq2_sig_res_taxo_seqs_family <- anti_join(muck_DESeq2_sig_res_taxo_seqs_family, unfiltered_DESeq2_estuary_familys, by = "ESV.ID")
filtered_muck_DESeq2_sig_res_taxo_seqs_family <- anti_join(filtered_muck_DESeq2_sig_res_taxo_seqs_family, unfiltered_DESeq2_tomcu_familys, by = "ESV.ID")

#Filter out rows in the LOI.CU DESEQ object that are also in the IRL.SLE and Muck DESEQs#
filtered_tomcu_DESeq2_sig_res_taxo_seqs_family <- anti_join(tomcu_DESeq2_sig_res_taxo_seqs_family, unfiltered_DESeq2_estuary_familys, by = "ESV.ID")
filtered_tomcu_DESeq2_sig_res_taxo_seqs_family <- anti_join(filtered_tomcu_DESeq2_sig_res_taxo_seqs_family, unfiltered_DESeq2_muck_familys, by = "ESV.ID")



#Order rows by log2FoldChange#
filtered_muck_DESeq2_sig_res_taxo_seqs_family <- arrange(filtered_muck_DESeq2_sig_res_taxo_seqs_family, log2FoldChange)
filtered_tomcu_DESeq2_sig_res_taxo_seqs_family <- arrange(filtered_tomcu_DESeq2_sig_res_taxo_seqs_family, log2FoldChange)
filtered_estuary_DESeq2_sig_res_taxo_seqs_family <- arrange(filtered_estuary_DESeq2_sig_res_taxo_seqs_family, log2FoldChange)

#Make a csv file for each of the filtered tables#
write.csv(as.data.frame(filtered_estuary_DESeq2_sig_res_taxo_seqs_family), file="filtered_estuary_DESeq2_sig_res_taxo_seqs_family.csv")
write.csv(as.data.frame(filtered_muck_DESeq2_sig_res_taxo_seqs_family), file="filtered_muck_DESeq2_sig_res_taxo_seqs_family.csv")
write.csv(as.data.frame(filtered_tomcu_DESeq2_sig_res_taxo_seqs_family), file="filtered_tomcu_DESeq2_sig_res_taxo_seqs_family.csv")



##Determine number of filtered indicators##

#Number of Not indicators# 
length(which(filtered_muck_DESeq2_sig_res_taxo_seqs_family$log2FoldChange > 0)) #26
#Number of Muck indicators#
length(which(filtered_muck_DESeq2_sig_res_taxo_seqs_family$log2FoldChange < 0)) #25
#Number of High.Low indicators# 
length(which(filtered_tomcu_DESeq2_sig_res_taxo_seqs_family$log2FoldChange > 0)) #11
#Number of High.High indicators#
length(which(filtered_tomcu_DESeq2_sig_res_taxo_seqs_family$log2FoldChange < 0)) #9
#Number of SLE indicators# 
length(which(filtered_estuary_DESeq2_sig_res_taxo_seqs_family$log2FoldChange > 0)) #57
#Number of IRL indicators#
length(which(filtered_estuary_DESeq2_sig_res_taxo_seqs_family$log2FoldChange < 0)) #66

#Create lists of filtered significant familys for each category#
filtered_DESeq2_estuary_familys <- subset(filtered_estuary_DESeq2_sig_res_taxo_seqs_family, select=c(ESV.ID))
filtered_DESeq2_muck_familys <- subset(filtered_muck_DESeq2_sig_res_taxo_seqs_family, select=c(ESV.ID))
filtered_DESeq2_tomcu_familys <- subset(filtered_tomcu_DESeq2_sig_res_taxo_seqs_family, select=c(ESV.ID))

##Create phyloseq objects for DESeq2 indicators##

#Create phyloseq objects of the unfiltered indicators using character string version of significant familys##
charac_unfiltered_DESeq2_estuary_familys <- unfiltered_DESeq2_estuary_familys[,"ESV.ID"]
charac_unfiltered_DESeq2_estuary_familys <- as.character(charac_unfiltered_DESeq2_estuary_familys)
UFDestuaryFamilyAfamily <- prune_taxa(charac_unfiltered_DESeq2_estuary_familys, RFLWSS_Y1_FamilyAfamily)

charac_unfiltered_DESeq2_muck_familys <- unfiltered_DESeq2_muck_familys[,"ESV.ID"]
charac_unfiltered_DESeq2_muck_familys <- as.character(charac_unfiltered_DESeq2_muck_familys)
UFDmuckFamilyAfamily <- prune_taxa(charac_unfiltered_DESeq2_muck_familys, RFLWSS_Y1_FamilyAfamily)

charac_unfiltered_DESeq2_tomcu_familys <- unfiltered_DESeq2_tomcu_familys[,"ESV.ID"]
charac_unfiltered_DESeq2_tomcu_familys <- as.character(charac_unfiltered_DESeq2_tomcu_familys)
UFDtomcuFamilyAfamily <- prune_taxa(charac_unfiltered_DESeq2_tomcu_familys, RFLWSS_Y1_FamilyAfamily)

#Create phyloseq objects of the filtered indicators using character string version of significant familys##
charac_filtered_DESeq2_estuary_familys <- filtered_DESeq2_estuary_familys[,"ESV.ID"]
charac_filtered_DESeq2_estuary_familys <- as.character(charac_filtered_DESeq2_estuary_familys)
FDestuaryFamilyAfamily <- prune_taxa(charac_filtered_DESeq2_estuary_familys, RFLWSS_Y1_FamilyAfamily)

charac_filtered_DESeq2_muck_familys <- filtered_DESeq2_muck_familys[,"ESV.ID"]
charac_filtered_DESeq2_muck_familys <- as.character(charac_filtered_DESeq2_muck_familys)
FDmuckFamilyAfamily <- prune_taxa(charac_filtered_DESeq2_muck_familys, RFLWSS_Y1_FamilyAfamily)

charac_filtered_DESeq2_tomcu_familys <- filtered_DESeq2_tomcu_familys[,"ESV.ID"]
charac_filtered_DESeq2_tomcu_familys <- as.character(charac_filtered_DESeq2_tomcu_familys)
FDtomcuFamilyAfamily <- prune_taxa(charac_filtered_DESeq2_tomcu_familys, RFLWSS_Y1_FamilyAfamily)



##Create heatmaps to see overall trends, may take a while to load##

unfiltered_DESeq2_estuary_family_heatmap <- plot_heatmap(UFDestuaryFamilyAfamily, taxa.order = charac_unfiltered_DESeq2_estuary_familys, sample.order= "Estuary", sample.label = "Estuary")
unfiltered_DESeq2_estuary_family_heatmap

unfiltered_DESeq2_muck_family_heatmap <- plot_heatmap(UFDmuckFamilyAfamily, taxa.order = charac_unfiltered_DESeq2_muck_familys, sample.order= "Muck", sample.label = "Muck")
unfiltered_DESeq2_muck_family_heatmap

unfiltered_DESeq2_tomcu_family_heatmap <- plot_heatmap(UFDtomcuFamilyAfamily, taxa.order = charac_unfiltered_DESeq2_tomcu_familys, sample.order= "LOI.Cu", sample.label = "LOI.Cu")
unfiltered_DESeq2_tomcu_family_heatmap

filtered_DESeq2_estuary_family_heatmap <- plot_heatmap(FDestuaryFamilyAfamily, taxa.order = charac_filtered_DESeq2_estuary_familys, sample.order= "Estuary", sample.label = "Estuary")
filtered_DESeq2_estuary_family_heatmap

filtered_DESeq2_muck_family_heatmap <- plot_heatmap(FDmuckFamilyAfamily, taxa.order = charac_filtered_DESeq2_muck_familys, sample.order= "Muck", sample.label = "Muck")
filtered_DESeq2_muck_family_heatmap

filtered_DESeq2_tomcu_family_heatmap <- plot_heatmap(FDtomcuFamilyAfamily, taxa.order = charac_filtered_DESeq2_tomcu_familys, sample.order= "LOI.Cu", sample.label = "LOI.Cu")
filtered_DESeq2_tomcu_family_heatmap


###CONDUCT INDICSPECIES ANALYSIS AT family LEVEL###

#Increase max print so that all of results can be printed#
options(max.print=100000000)


##Conduct family indicspecies analysis for Estuary Category#

#Take out OTU aka family table from metadata focused phyloseq object#
estuary_seqs_family = as(otu_table(RFLWSS_Y1_FamilyAfamily), "matrix")

#Make it a data frame, transpose it, make it a data frame again#
estuary_seqs_family <- as.data.frame(estuary_seqs_family)
estuary_seqs_family<- t(estuary_seqs_family)
estuary_seqs_family <- as.data.frame(estuary_seqs_family)

#Run indicspecies with same estuary_meta_group as in asv section#
estuary_indval_family = multipatt(estuary_seqs_family, estuary_meta_group, control = how(nperm=9999))

#Create a blank csv, then sink all indicspecies results into it#
sink("estuary_IS_res_family.csv")
estuary_sig_indval_family <- summary(estuary_indval_family, indvalcomp=TRUE, alpha=1)
sink()

##Conduct family indicspecies analysis for Muck Category#

#Take out OTU aka family table from metadata focused phyloseq object#
muck_seqs_family = as(otu_table(RFLWSS_Y1_MuckvsNot_FamilyAfamily), "matrix")

#Make it a data frame, transpose it, make it a data frame again#
muck_seqs_family <- as.data.frame(muck_seqs_family)
muck_seqs_family<- t(muck_seqs_family)
muck_seqs_family <- as.data.frame(muck_seqs_family)

#Run indicspecies with same muck_meta_group as asv section#
muck_indval_family = multipatt(muck_seqs_family, muck_meta_group, control = how(nperm=9999))

#Create a blank csv, then sink all indicspecies results into it#
sink("muck_IS_res_family.csv")
muck_sig_indval_family <- summary(muck_indval_family, indvalcomp=TRUE, alpha=1)
sink()


##Conduct family indicspecies analysis for TOM/Cu Category#

#Take out OTU aka family table from metadata focused phyloseq object#
tomcu_seqs_family = as(otu_table(RFLWSS_Y1_HHvsHL_FamilyAfamily), "matrix")

#Make it a data frame, transpose it, make it a data frame again#
tomcu_seqs_family <- as.data.frame(tomcu_seqs_family)
tomcu_seqs_family<- t(tomcu_seqs_family)
tomcu_seqs_family <- as.data.frame(tomcu_seqs_family)

#Run indicspecies with same tomcu_meta_group as asv section#
tomcu_indval_family = multipatt(tomcu_seqs_family, tomcu_meta_group, control = how(nperm=9999))

#Create a blank csv, then sink all indicspecies results into it#
sink("tomcu_IS_res_family.csv")
tomcu_sig_indval_family <- summary(tomcu_indval_family, indvalcomp=TRUE, alpha=1)
sink()

#Reset max print to default#
options(max.print=99999)

##Look at csv file you sinked into and create a table
##of just family.ID, A, B, stat, p.value, and group##

##Continue family indicspecies analysis for Estuary Category##

#Upload new table#
estuary_IS_res_family <- read.csv(file.choose())
#Make a new p value adjusted column#
estuary_IS_res_family$padj <- p.adjust(estuary_IS_res_family$p.value, method = "BH")

#Keep table of adjusted p values less than 0.05#
estuary_IS_sig_res_family <- filter(estuary_IS_res_family, padj < 0.05)

#Make dataframe with taxanomy added in by making ESV.ID the rownames then adding taxonomy table#
estuary_IS_sig_res_taxo_family <- as.data.frame(estuary_IS_sig_res_family)
rownames(estuary_IS_sig_res_taxo_family) <- estuary_IS_sig_res_taxo_family[,1]
estuary_IS_sig_res_taxo_family = cbind(as(estuary_IS_sig_res_taxo_family, "data.frame"), as(tax_table(RFLWSS_Y1_FamilyAfamily)[rownames(estuary_IS_sig_res_taxo_family), ], "matrix"))

#Make dataframe with familys added in from all sites#
estuary_IS_sig_res_taxo_seqs_family = cbind(as(estuary_IS_sig_res_taxo_family, "data.frame"), as(otu_table(RFLWSS_Y1_FamilyAfamily)[rownames(estuary_IS_sig_res_taxo_family), ], "matrix"))

#Remove old rownames#
rownames(estuary_IS_sig_res_taxo_seqs_family) <-NULL

#Make a txt file of results for all sites#
write.csv(as.data.frame(estuary_IS_sig_res_taxo_seqs_family), file="estuary_IS_sig_res_taxo_seqs_family.csv")

#Check the number of indicators in each group after padj#
length(which(sample_data(estuary_IS_sig_res_family)$group == "1")) #74 IRL
length(which(sample_data(estuary_IS_sig_res_family)$group == "2")) #200 SLE


##Continue family indicspecies analysis for Muck Category##

#Upload new table#
muck_IS_res_family <- read.csv(file.choose())

#Make a new p value adjusted column#
muck_IS_res_family$padj <- p.adjust(muck_IS_res_family$p.value, method = "BH")

#Keep table of adjusted p values less than 0.05#
muck_IS_sig_res_family <- filter(muck_IS_res_family, padj < 0.05)

#Make dataframe with taxanomy added in by making ESV.ID the rownames then adding taxonomy table#
muck_IS_sig_res_taxo_family <- as.data.frame(muck_IS_sig_res_family)
rownames(muck_IS_sig_res_taxo_family) <- muck_IS_sig_res_taxo_family[,1]
muck_IS_sig_res_taxo_family = cbind(as(muck_IS_sig_res_taxo_family, "data.frame"), as(tax_table(RFLWSS_Y1_MuckvsNot_FamilyAfamily)[rownames(muck_IS_sig_res_taxo_family), ], "matrix"))

#Make dataframe with familys added in from all sites#
muck_IS_sig_res_taxo_seqs_family = cbind(as(muck_IS_sig_res_taxo_family, "data.frame"), as(otu_table(RFLWSS_Y1_FamilyAfamily)[rownames(muck_IS_sig_res_taxo_family), ], "matrix"))

#Remove old rownames#
rownames(muck_IS_sig_res_taxo_seqs_family) <-NULL

#Make a txt file of results for all sites#
write.csv(as.data.frame(muck_IS_sig_res_taxo_seqs_family), file="muck_IS_sig_res_taxo_seqs_family.csv")

#Check the number of indicators in each group after padj#
length(which(sample_data(muck_IS_sig_res_family)$group == "1")) #112 Muck
length(which(sample_data(muck_IS_sig_res_family)$group == "2")) #36 Not

##Continue family indicspecies analysis for TOM/Cu Category##

#Upload new table#
tomcu_IS_res_family <- read.csv(file.choose())
#Make a new p value adjusted column#
tomcu_IS_res_family$padj <- p.adjust(tomcu_IS_res_family$p.value, method = "BH")

#Keep table of adjusted p values less than 0.05#
tomcu_IS_sig_res_family <- filter(tomcu_IS_res_family, padj < 0.05)

#Make dataframe with taxanomy added in by making ESV.ID the rownames then adding taxonomy table#
tomcu_IS_sig_res_taxo_family <- as.data.frame(tomcu_IS_sig_res_family)
rownames(tomcu_IS_sig_res_taxo_family) <- tomcu_IS_sig_res_taxo_family[,1]
tomcu_IS_sig_res_taxo_family = cbind(as(tomcu_IS_sig_res_taxo_family, "data.frame"), as(tax_table(RFLWSS_Y1_MuckvsNot_FamilyAfamily)[rownames(tomcu_IS_sig_res_taxo_family), ], "matrix"))

#Make dataframe with familys added in from all sites#
tomcu_IS_sig_res_taxo_seqs_family = cbind(as(tomcu_IS_sig_res_taxo_family, "data.frame"), as(otu_table(RFLWSS_Y1_FamilyAfamily)[rownames(tomcu_IS_sig_res_taxo_family), ], "matrix"))

#Remove old rownames#
rownames(tomcu_IS_sig_res_taxo_seqs_family) <-NULL

#Make a txt file of results for all sites#
write.csv(as.data.frame(tomcu_IS_sig_res_taxo_seqs_family), file="tomcu_IS_sig_res_taxo_seqs_family.csv")

#Check the number of indicators in each group after padj#
length(which(sample_data(tomcu_IS_sig_res_family)$group == "1")) #53 High/High
length(which(sample_data(tomcu_IS_sig_res_family)$group == "2")) #55 High/Low

#Order rows by group#
estuary_IS_sig_res_taxo_seqs_family <- arrange(estuary_IS_sig_res_taxo_seqs_family, group)
muck_IS_sig_res_taxo_seqs_family <- arrange(muck_IS_sig_res_taxo_seqs_family, group)
tomcu_IS_sig_res_taxo_seqs_family <- arrange(tomcu_IS_sig_res_taxo_seqs_family, group)


#Create lists of just significant familys for each category#
unfiltered_IS_estuary_familys <- subset(estuary_IS_sig_res_taxo_seqs_family, select=c(ESV.ID))
unfiltered_IS_muck_familys <- subset(muck_IS_sig_res_taxo_seqs_family, select=c(ESV.ID))
unfiltered_IS_tomcu_familys <- subset(tomcu_IS_sig_res_taxo_seqs_family, select=c(ESV.ID))


##Create indicator tables filtered of overlapping indicators with other two metadata categories##

#Filter out rows in the estuary IS object that are also in the tomcu and muck IS objects#
filtered_estuary_IS_sig_res_taxo_seqs_family <- anti_join(estuary_IS_sig_res_taxo_seqs_family, unfiltered_IS_tomcu_familys, by = "ESV.ID")
filtered_estuary_IS_sig_res_taxo_seqs_family <- anti_join(filtered_estuary_IS_sig_res_taxo_seqs_family, unfiltered_IS_muck_familys, by = "ESV.ID")

#Filter out rows in the muck IS object that are also in the estuary and tomcu IS objects#
filtered_muck_IS_sig_res_taxo_seqs_family<- anti_join(muck_IS_sig_res_taxo_seqs_family, unfiltered_IS_estuary_familys, by = "ESV.ID")
filtered_muck_IS_sig_res_taxo_seqs_family <- anti_join(filtered_muck_IS_sig_res_taxo_seqs_family, unfiltered_IS_tomcu_familys, by = "ESV.ID")

#Filter out rows in the tomcu IS object that are also in the estuary and muck IS objects#
filtered_tomcu_IS_sig_res_taxo_seqs_family <- anti_join(tomcu_IS_sig_res_taxo_seqs_family, unfiltered_IS_estuary_familys, by = "ESV.ID")
filtered_tomcu_IS_sig_res_taxo_seqs_family <- anti_join(filtered_tomcu_IS_sig_res_taxo_seqs_family, unfiltered_IS_muck_familys, by = "ESV.ID")

#Order rows by group#
filtered_estuary_IS_sig_res_taxo_seqs_family <- arrange(filtered_estuary_IS_sig_res_taxo_seqs_family, group)
filtered_muck_IS_sig_res_taxo_seqs_family <- arrange(filtered_muck_IS_sig_res_taxo_seqs_family, group)
filtered_tomcu_IS_sig_res_taxo_seqs_family <- arrange(filtered_tomcu_IS_sig_res_taxo_seqs_family, group)


#Make a csv file for each of the filtered tables#
write.csv(as.data.frame(filtered_estuary_IS_sig_res_taxo_seqs_family), file="filtered_estuary_IS_sig_res_taxo_seqs_family.csv")
write.csv(as.data.frame(filtered_muck_IS_sig_res_taxo_seqs_family), file="filtered_muck_IS_sig_res_taxo_seqs_family.csv")
write.csv(as.data.frame(filtered_tomcu_IS_sig_res_taxo_seqs_family), file="filtered_tomcu_IS_sig_res_taxo_seqs_family.csv")

#Determine the number of indicators in each category#
length(which(sample_data(filtered_estuary_IS_sig_res_taxo_seqs_family)$group == "1")) #37 IRL
length(which(sample_data(filtered_estuary_IS_sig_res_taxo_seqs_family)$group == "2")) #122 SLE
length(which(sample_data(filtered_muck_IS_sig_res_taxo_seqs_family)$group == "1")) #55 Muck
length(which(sample_data(filtered_muck_IS_sig_res_taxo_seqs_family)$group == "2")) #20 Not
length(which(sample_data(filtered_tomcu_IS_sig_res_taxo_seqs_family)$group == "1")) #12 High/High
length(which(sample_data(filtered_tomcu_IS_sig_res_taxo_seqs_family)$group == "2")) #5 High/Low


#Create lists of filtered significant familys for each category#
filtered_IS_estuary_familys <- subset(filtered_estuary_IS_sig_res_taxo_seqs_family, select=c(ESV.ID))
filtered_IS_muck_familys <- subset(filtered_muck_IS_sig_res_taxo_seqs_family, select=c(ESV.ID))
filtered_IS_tomcu_familys <- subset(filtered_tomcu_IS_sig_res_taxo_seqs_family, select=c(ESV.ID))

#Create phyloseq objects of the unfiltered indicators using character string version of significant familys##
charac_unfiltered_IS_estuary_familys <- unfiltered_IS_estuary_familys[,"ESV.ID"]
charac_unfiltered_IS_estuary_familys <- as.character(charac_unfiltered_IS_estuary_familys)
UFISestuaryFamilyAfamily <- prune_taxa(charac_unfiltered_IS_estuary_familys, RFLWSS_Y1_FamilyAfamily)

charac_unfiltered_IS_muck_familys <- unfiltered_IS_muck_familys[,"ESV.ID"]
charac_unfiltered_IS_muck_familys <- as.character(charac_unfiltered_IS_muck_familys)
UFISmuckFamilyAfamily <- prune_taxa(charac_unfiltered_IS_muck_familys, RFLWSS_Y1_FamilyAfamily)

charac_unfiltered_IS_tomcu_familys <- unfiltered_IS_tomcu_familys[,"ESV.ID"]
charac_unfiltered_IS_tomcu_familys <- as.character(charac_unfiltered_IS_tomcu_familys)
UFIStomcuFamilyAfamily <- prune_taxa(charac_unfiltered_IS_tomcu_familys, RFLWSS_Y1_FamilyAfamily)

#Create phyloseq objects of the filtered indicators using character string version of significant familys##
charac_filtered_IS_estuary_familys <- filtered_IS_estuary_familys[,"ESV.ID"]
charac_filtered_IS_estuary_familys <- as.character(charac_filtered_IS_estuary_familys)
FISestuaryFamilyAfamily <- prune_taxa(charac_filtered_IS_estuary_familys, RFLWSS_Y1_FamilyAfamily)

charac_filtered_IS_muck_familys <- filtered_IS_muck_familys[,"ESV.ID"]
charac_filtered_IS_muck_familys <- as.character(charac_filtered_IS_muck_familys)
FISmuckFamilyAfamily <- prune_taxa(charac_filtered_IS_muck_familys, RFLWSS_Y1_FamilyAfamily)

charac_filtered_IS_tomcu_familys <- filtered_IS_tomcu_familys[,"ESV.ID"]
charac_filtered_IS_tomcu_familys <- as.character(charac_filtered_IS_tomcu_familys)
FIStomcuFamilyAfamily <- prune_taxa(charac_filtered_IS_tomcu_familys, RFLWSS_Y1_FamilyAfamily)

##Create heatmaps to see overall trends, may take a while to load##

unfiltered_IS_estuary_family_heatmap <- plot_heatmap(UFDestuaryFamilyAfamily, taxa.order = charac_unfiltered_IS_estuary_familys, sample.order= "Estuary", sample.label = "Estuary")
unfiltered_IS_estuary_family_heatmap

unfiltered_IS_muck_family_heatmap <- plot_heatmap(UFDmuckFamilyAfamily, taxa.order = charac_unfiltered_IS_muck_familys, sample.order= "Muck", sample.label = "Muck")
unfiltered_IS_muck_family_heatmap

unfiltered_IS_tomcu_family_heatmap <- plot_heatmap(UFDtomcuFamilyAfamily, taxa.order = charac_unfiltered_IS_tomcu_familys, sample.order= "LOI.Cu", sample.label = "LOI.Cu")
unfiltered_IS_tomcu_family_heatmap

filtered_IS_estuary_family_heatmap <- plot_heatmap(FDestuaryFamilyAfamily, taxa.order = charac_filtered_IS_estuary_familys, sample.order= "Estuary", sample.label = "Estuary")
filtered_IS_estuary_family_heatmap

filtered_IS_muck_family_heatmap <- plot_heatmap(FDmuckFamilyAfamily, taxa.order = charac_filtered_IS_muck_familys, sample.order= "Muck", sample.label = "Muck")
filtered_IS_muck_family_heatmap

filtered_IS_tomcu_family_heatmap <- plot_heatmap(FDtomcuFamilyAfamily, taxa.order = charac_filtered_IS_tomcu_familys, sample.order= "LOI.Cu", sample.label = "LOI.Cu")
filtered_IS_tomcu_family_heatmap

###MAKE COMBO INDICATOR LISTS###

##Make "Combo" indicator lists from overlapping unfiltered DESeq2 and IS indicators##

#Keep rows in the estuary DESEQ object that are also in the IRLSLE IS object#
unfiltered_estuary_combo_sig_res_taxo_seqs_family <- merge(estuary_IS_sig_res_family,estuary_DESeq2_sig_res_taxo_seqs_family, by="ESV.ID")

#Keep rows in the muck DESEQ object that are also in the Muck IS object#
unfiltered_muck_combo_sig_res_taxo_seqs_family <- merge(muck_IS_sig_res_family,muck_DESeq2_sig_res_taxo_seqs_family, by="ESV.ID")

#Keep rows in the tomcu DESEQ object that are also in the LOICu IS object#
unfiltered_tomcu_combo_sig_res_taxo_seqs_family <- merge(tomcu_IS_sig_res_family,tomcu_DESeq2_sig_res_taxo_seqs_family, by="ESV.ID")

##Make "Combo" indicator lists from overlapping filtered DESeq2 and IS indicators##

#Create focused filtered IS tables without taxonomy or sequences#
filtered_estuary_IS_sig_res_family <- subset(filtered_estuary_IS_sig_res_taxo_seqs_family, select=c(ESV.ID,A,B,stat,p.value,group,padj))

filtered_muck_IS_sig_res_family <- subset(filtered_muck_IS_sig_res_taxo_seqs_family, select=c(ESV.ID,A,B,stat,p.value,group,padj))

filtered_tomcu_IS_sig_res_family <- subset(filtered_tomcu_IS_sig_res_taxo_seqs_family, select=c(ESV.ID,A,B,stat,p.value,group,padj))

#Keep rows in the estuary DESEQ object that are also in the IRLSLE IS object#
filtered_estuary_combo_sig_res_taxo_seqs_family <- merge(filtered_estuary_IS_sig_res_family,filtered_estuary_DESeq2_sig_res_taxo_seqs_family, by="ESV.ID")

#Keep rows in the muck DESEQ object that are also in the Muck IS object#
filtered_muck_combo_sig_res_taxo_seqs_family <- merge(filtered_muck_IS_sig_res_family,filtered_muck_DESeq2_sig_res_taxo_seqs_family, by="ESV.ID")

#Keep rows in the tomcu DESEQ object that are also in the LOICu IS object#
filtered_tomcu_combo_sig_res_taxo_seqs_family <- merge(filtered_tomcu_IS_sig_res_family,filtered_tomcu_DESeq2_sig_res_taxo_seqs_family, by="ESV.ID")

#Determine the number of indicators in each subcategory of the unfiltered Combo tables#
length(which(sample_data(unfiltered_estuary_combo_sig_res_taxo_seqs_family)$group == "1")) #66 IRL
length(which(sample_data(unfiltered_estuary_combo_sig_res_taxo_seqs_family)$group == "2")) #109 SLE
length(which(sample_data(unfiltered_muck_combo_sig_res_taxo_seqs_family)$group == "1")) #30 Muck
length(which(sample_data(unfiltered_muck_combo_sig_res_taxo_seqs_family)$group == "2")) #32 Not
length(which(sample_data(unfiltered_tomcu_combo_sig_res_taxo_seqs_family)$group == "1")) #28 High/High
length(which(sample_data(unfiltered_tomcu_combo_sig_res_taxo_seqs_family)$group == "2")) #54 High/Low

#Determine the number of indicators in each subcategory of the filtered Combo tables#
length(which(sample_data(filtered_estuary_combo_sig_res_taxo_seqs_family)$group == "1")) #20 IRL
length(which(sample_data(filtered_estuary_combo_sig_res_taxo_seqs_family)$group == "2")) #37 SLE
length(which(sample_data(filtered_muck_combo_sig_res_taxo_seqs_family)$group == "1")) #9 Muck
length(which(sample_data(filtered_muck_combo_sig_res_taxo_seqs_family)$group == "2")) #9 Not
length(which(sample_data(filtered_tomcu_combo_sig_res_taxo_seqs_family)$group == "1")) #2 High/High
length(which(sample_data(filtered_tomcu_combo_sig_res_taxo_seqs_family)$group == "2")) #3 High/Low


#Order unfiltered Combo table rows by group#
unfiltered_estuary_combo_sig_res_taxo_seqs_family <- arrange(unfiltered_estuary_combo_sig_res_taxo_seqs_family, group)
unfiltered_muck_combo_sig_res_taxo_seqs_family <- arrange(unfiltered_muck_combo_sig_res_taxo_seqs_family, group)
unfiltered_tomcu_combo_sig_res_taxo_seqs_family <- arrange(unfiltered_tomcu_combo_sig_res_taxo_seqs_family, group)

#Order filtered Combo table rows by group#
filtered_estuary_combo_sig_res_taxo_seqs_family <- arrange(filtered_estuary_combo_sig_res_taxo_seqs_family, group)
filtered_muck_combo_sig_res_taxo_seqs_family <- arrange(filtered_muck_combo_sig_res_taxo_seqs_family, group)
filtered_tomcu_combo_sig_res_taxo_seqs_family <- arrange(filtered_tomcu_combo_sig_res_taxo_seqs_family, group)

#Make a csv file for each of the unfiltered Combo tables#
write.csv(as.data.frame(unfiltered_estuary_combo_sig_res_taxo_seqs_family), file="unfiltered_estuary_combo_sig_res_taxo_seqs_family.csv")
write.csv(as.data.frame(unfiltered_muck_combo_sig_res_taxo_seqs_family), file="unfiltered_muck_combo_sig_res_taxo_seqs_family.csv")
write.csv(as.data.frame(unfiltered_tomcu_combo_sig_res_taxo_seqs_family), file="unfiltered_tomcu_combo_sig_res_taxo_seqs_family.csv")

#Make a csv file for each of the filtered Combo tables#
write.csv(as.data.frame(filtered_estuary_combo_sig_res_taxo_seqs_family), file="filtered_estuary_combo_sig_res_taxo_seqs_family.csv")
write.csv(as.data.frame(filtered_muck_combo_sig_res_taxo_seqs_family), file="filtered_muck_combo_sig_res_taxo_seqs_family.csv")
write.csv(as.data.frame(filtered_tomcu_combo_sig_res_taxo_seqs_family), file="filtered_tomcu_combo_sig_res_taxo_seqs_family.csv")

#Create lists of unfiltered significant familys for each category#
unfiltered_Combo_estuary_familys <- subset(unfiltered_estuary_combo_sig_res_taxo_seqs_family, select=c(ESV.ID))
unfiltered_Combo_muck_familys <- subset(unfiltered_muck_combo_sig_res_taxo_seqs_family, select=c(ESV.ID))
unfiltered_Combo_tomcu_familys <- subset(unfiltered_tomcu_combo_sig_res_taxo_seqs_family, select=c(ESV.ID))

#Create lists of filtered significant familys for each category#
filtered_Combo_estuary_familys <- subset(filtered_estuary_combo_sig_res_taxo_seqs_family, select=c(ESV.ID))
filtered_Combo_muck_familys <- subset(filtered_muck_combo_sig_res_taxo_seqs_family, select=c(ESV.ID))
filtered_Combo_tomcu_familys <- subset(filtered_tomcu_combo_sig_res_taxo_seqs_family, select=c(ESV.ID))

##Create phyloseq objects for Combo indicators##

#Create phyloseq objects of the unfiltered indicators using character string version of significant familys##
charac_unfiltered_Combo_estuary_familys <- unfiltered_Combo_estuary_familys[,"ESV.ID"]
charac_unfiltered_Combo_estuary_familys <- as.character(charac_unfiltered_Combo_estuary_familys)
UFCestuaryFamilyAfamily <- prune_taxa(charac_unfiltered_Combo_estuary_familys, RFLWSS_Y1_FamilyAfamily)

charac_unfiltered_Combo_muck_familys <- unfiltered_Combo_muck_familys[,"ESV.ID"]
charac_unfiltered_Combo_muck_familys <- as.character(charac_unfiltered_Combo_muck_familys)
UFCmuckFamilyAfamily <- prune_taxa(charac_unfiltered_Combo_muck_familys, RFLWSS_Y1_FamilyAfamily)

charac_unfiltered_Combo_tomcu_familys <- unfiltered_Combo_tomcu_familys[,"ESV.ID"]
charac_unfiltered_Combo_tomcu_familys <- as.character(charac_unfiltered_Combo_tomcu_familys)
UFCtomcuFamilyAfamily <- prune_taxa(charac_unfiltered_Combo_tomcu_familys, RFLWSS_Y1_FamilyAfamily)

#Create phyloseq objects of the filtered indicators using character string version of significant familys##
charac_filtered_Combo_estuary_familys <- filtered_Combo_estuary_familys[,"ESV.ID"]
charac_filtered_Combo_estuary_familys <- as.character(charac_filtered_Combo_estuary_familys)
FCestuaryFamilyAfamily <- prune_taxa(charac_filtered_Combo_estuary_familys, RFLWSS_Y1_FamilyAfamily)

charac_filtered_Combo_muck_familys <- filtered_Combo_muck_familys[,"ESV.ID"]
charac_filtered_Combo_muck_familys <- as.character(charac_filtered_Combo_muck_familys)
FCmuckFamilyAfamily <- prune_taxa(charac_filtered_Combo_muck_familys, RFLWSS_Y1_FamilyAfamily)

charac_filtered_Combo_tomcu_familys <- filtered_Combo_tomcu_familys[,"ESV.ID"]
charac_filtered_Combo_tomcu_familys <- as.character(charac_filtered_Combo_tomcu_familys)
FCtomcuFamilyAfamily <- prune_taxa(charac_filtered_Combo_tomcu_familys, RFLWSS_Y1_FamilyAfamily)

##Create heatmaps to see overall trends, may take a while to load##

unfiltered_Combo_estuary_family_heatmap <- plot_heatmap(UFDestuaryFamilyAfamily, taxa.order = charac_unfiltered_Combo_estuary_familys, sample.order= "Estuary", sample.label = "Estuary")
unfiltered_Combo_estuary_family_heatmap

unfiltered_Combo_muck_family_heatmap <- plot_heatmap(UFDmuckFamilyAfamily, taxa.order = charac_unfiltered_Combo_muck_familys, sample.order= "Muck", sample.label = "Muck")
unfiltered_Combo_muck_family_heatmap

unfiltered_Combo_tomcu_family_heatmap <- plot_heatmap(UFDtomcuFamilyAfamily, taxa.order = charac_unfiltered_Combo_tomcu_familys, sample.order= "LOI.Cu", sample.label = "LOI.Cu")
unfiltered_Combo_tomcu_family_heatmap

filtered_Combo_estuary_family_heatmap <- plot_heatmap(FDestuaryFamilyAfamily, taxa.order = charac_filtered_Combo_estuary_familys, sample.order= "Estuary", sample.label = "Estuary")
filtered_Combo_estuary_family_heatmap

filtered_Combo_muck_family_heatmap <- plot_heatmap(FDmuckFamilyAfamily, taxa.order = charac_filtered_Combo_muck_familys, sample.order= "Muck", sample.label = "Muck")
filtered_Combo_muck_family_heatmap

filtered_Combo_tomcu_family_heatmap <- plot_heatmap(FDtomcuFamilyAfamily, taxa.order = charac_filtered_Combo_tomcu_familys, sample.order= "LOI.Cu", sample.label = "LOI.Cu")
filtered_Combo_tomcu_family_heatmap

###INDICATOR EFFECTIVENESS METHOD###

##NOTES ON METHOD##
#"affected' = subcategory that is more affected by the stressor you are      concerned about. For example SLE for the Estuary category becasue it has     more freshwater discharges, 3MC for the Muck category because it has more    muck characteristics, HiHi for the TOM/Cu category becasue it has more       copper#

#"non.affected" = IRL, 0MC, HiLo#

#microbially.predicted "affected" sample = samples in the partitioning around medoids (PAM) cluster with the most metadata.defined "affected" samples#

#microbially.predicted "non.affected" sample = samples in the PAM cluster with the most metadata.defined "non.affected" samples#

#Overall Idea: Use PAM clustering to split samples between the "affected" and "non.affected" clusters based upon the indicators (microbially.predicted) and see how well this classification matches the metadata.defined classification by using the product of specificity and sensitivity#

#Copy results to an Excel Sheet for further statistical testing, each combo of factors should have four percentage types (Product is Sensitivity x Specificity), for example:#
#Indicator_Test	Taxonomic_Level	Filtering_Status	Metadata_Tested	Percent_Type	Percentage#
#Indicfamily   family             Unfiltered        Estuary         Total         98.78#
#Indicfamily   family             Unfiltered        Estuary         Sensitivity   96.23#
#Indicfamily   family             Unfiltered        Estuary         Specificity   99.15#
#Indicfamily   family             Unfiltered        Estuary         Product       95.41#


###ORIGINAL DATASETS
###Split samples of orginal datasets into two clusters to see how it will compare to splitting based upon indicators###

##Test on Muck and Estuary categories at family level##

#Remove any samples with no familys from target phyloseq#
nsamples(RFLWSS_Y1_FamilyAfamily) #90
RFLWSS_Y1_FamilyAfamily = prune_samples(sample_sums(RFLWSS_Y1_FamilyAfamily)>=1, RFLWSS_Y1_FamilyAfamily)
nsamples(RFLWSS_Y1_FamilyAfamily) #90

#Make a matrix with the sample_data and relabel subcategories to make metadata.defined column#
Original_Est_Muck_family <- as.matrix(sample_data(RFLWSS_Y1_FamilyAfamily))
Original_Est_Muck_family[Original_Est_Muck_family ==  "Muck"] <- 2
Original_Est_Muck_family[Original_Est_Muck_family ==  "Not"] <- 1
Original_Est_Muck_family[Original_Est_Muck_family ==  "Mucky"] <- 2
Original_Est_Muck_family[Original_Est_Muck_family ==  "Muckish"] <- 2
Original_Est_Muck_family[Original_Est_Muck_family ==  "IRL"] <- 1
Original_Est_Muck_family[Original_Est_Muck_family ==  "SLE"] <- 2
Original_Est_Muck_family<-as.data.frame(Original_Est_Muck_family)


#Create Bray.Curtis distance matrix from family table#
bc <- vegdist(t(otu_table(RFLWSS_Y1_FamilyAfamily)))

#Conduct pam clustering based upon 2 clusters#
pam.res <- pam(bc,2, diss=TRUE)

#Add a microbially.predicted column (Orginal_family) to matrix above for comparison to metadata.defined column#
Original_Est_Muck_family <- cbind(Original_Est_Muck_family, Orginal_family = pam.res$cluster)

##Test Estuary category##

#Test  to see how well the metadata.defined and microbially.predicted columns match one another aka Total Efficiency; if the number is below 0.5, then switch the numbers in the metadata.defiend column#

#Test for Estuary Effeciency# 0.9222222
sum(Original_Est_Muck_family$Estuary == Original_Est_Muck_family$Orginal_family)/nrow(Original_Est_Muck_family)

##Test indicator sensitivity aka the true postives (TP)/ (TP + false negatives (FN))##
#TP = metadata.defined "affected" sample correctly placed in the microbiall   y.predicted "affected" cluster#
#FN = metadata.defined "affected" sample incorrectly placed in the           microbially.predicted "non.affected" cluster#
#TP+FN = all metadata.defined "affected" samples

#Test for SLE Effeciency aka Sensitvity# 0.9583333
Orginal_SLE_family <- filter(Original_Est_Muck_family, Estuary == "2")
sum(Orginal_SLE_family$Estuary == Orginal_SLE_family$Orginal_family)/nrow(Orginal_SLE_family)

##Test indicator specificity aka the true negatives (TN)/ (TN + false positives (FP))##
#TP = metadata.defined "affected" sample correctly placed in the microbiall   y.predicted "affected" cluster#
#FN = metadata.defined "affected" sample incorrectly placed in the           microbially.predicted "non.affected" cluster#
#TP+FN = all metadata.defined "affected" samples

#Test for IRL Effeciency aka Specificity# 0.9090909
Orginal_IRL_family <- filter(Original_Est_Muck_family, Estuary == "1")
sum(Orginal_IRL_family$Estuary == Orginal_IRL_family$Orginal_family)/nrow(Orginal_IRL_family)

#PCoA with Color by Estuary#
TotalPCoA <- ordinate(RFLWSS_Y1_FamilyAfamily,"PCoA")
p = plot_ordination(RFLWSS_Y1_FamilyAfamily, TotalPCoA, color="Estuary") + 
  ggtitle("Total familys separated by Estuary Category PCoA with Bray.Curtis distance") +
  geom_text(aes(label=Original_Est_Muck_family$Orginal_family),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p

##Test Muck category##

#Test for Muck Total Efficiency# 0.7444444
sum(Original_Est_Muck_family$Muck == Original_Est_Muck_family$Orginal_family)/nrow(Original_Est_Muck_family)

#Test for 3MC Effeciency aka Sensitivity# 0.6
Orginal_3MC_family <- filter(Original_Est_Muck_family, Muck == "2")
sum(Orginal_3MC_family$Muck == Orginal_3MC_family$Orginal_family)/nrow(Orginal_3MC_family)

#Test for 0MC Effeciency aka Specificity# 0.8166667
Orginal_0MC_family <- filter(Original_Est_Muck_family, Muck == "1")
sum(Orginal_0MC_family$Muck == Orginal_0MC_family$Orginal_family)/nrow(Orginal_0MC_family)

#PCoA with Color by Muck#
TotalPCoA <- ordinate(RFLWSS_Y1_FamilyAfamily,"PCoA")
p = plot_ordination(RFLWSS_Y1_FamilyAfamily, TotalPCoA, color="Muck") + 
  ggtitle("Total familys separated by Muck Category PCoA with Bray.Curtis distance") +
  geom_text(aes(label=Original_Est_Muck_family$Orginal_family),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p + annotate("text", x = 0.3, y = 0.3, label = c("SN = 0.44, SP=0.84"))+ scale_color_discrete(name="Muck Characteristics", limits=c("Muck", "Mucky", "Muckish", "Not"), labels=c("3", "2", "1", "0"))

##Test TOM/Cu category##

#Remove any samples with no familys from target phyloseq#
nsamples(RFLWSS_Y1_HHvsHL_FamilyAfamily) #29
RFLWSS_Y1_HHvsHL_FamilyAfamily = prune_samples(sample_sums(RFLWSS_Y1_HHvsHL_FamilyAfamily)>=1, RFLWSS_Y1_HHvsHL_FamilyAfamily)
nsamples(RFLWSS_Y1_HHvsHL_FamilyAfamily) #29

#Make a matrix with the sample_data and relabel subcategories to make metadata.defined column#
Original_TOMCu_family <- as.matrix(sample_data(RFLWSS_Y1_HHvsHL_FamilyAfamily))
Original_TOMCu_family[Original_TOMCu_family ==  "High.High"] <- 1
Original_TOMCu_family[Original_TOMCu_family ==  "High.Low"] <- 2
Original_TOMCu_family<-as.data.frame(Original_TOMCu_family)

#Create Bray.Curtis distance matrix from family table#
bc <- vegdist(t(otu_table(RFLWSS_Y1_HHvsHL_FamilyAfamily)))

#Conduct PAM clustering based upon 2 clusters#
pam.res <- pam(bc,2, diss=TRUE)

#Add a microbially.predicted column to matrix above for comparison to metadata.defined column#
Original_TOMCu_family <- cbind(Original_TOMCu_family, Orginal_family = pam.res$cluster)

#Determine Total TOM/Cu Efficiency# 0.6551724
sum(Original_TOMCu_family$LOI.Cu == Original_TOMCu_family$Orginal_family)/nrow(Original_TOMCu_family)

#Determine HiHi Efficiency  aka Sensitivity# 0.5454545
Original_TOMCu_familyHH <- filter(Original_TOMCu_family, LOI.Cu == "1") 
sum(Original_TOMCu_familyHH$LOI.Cu == Original_TOMCu_familyHH$Orginal_family)/nrow(Original_TOMCu_familyHH)

#Determine HiLo Efficiency aka Specificity# 0.7222222
Original_TOMCu_familyHL <- filter(Original_TOMCu_family, LOI.Cu == "2") 
sum(Original_TOMCu_familyHL$LOI.Cu == Original_TOMCu_familyHL$Orginal_family)/nrow(Original_TOMCu_familyHL)

#PCoA with Color by TOM/Cu#
TotalPCoA <- ordinate(RFLWSS_Y1_HHvsHL_FamilyAfamily,"PCoA")
p = plot_ordination(RFLWSS_Y1_HHvsHL_FamilyAfamily, TotalPCoA, color="LOI.Cu") + 
  ggtitle("Total familys separated by LOI.Cu Category PCoA with Bray.Curtis distance") +
  geom_text(aes(label=Original_TOMCu_family$Orginal_family),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p

###TEST INDICSPECIES ESTUARY INDICATORS###

##Test unfiltered Estuary familys##

#Remove any samples with no familys from target phyloseq#
nsamples(UFISestuaryFamilyAfamily) #90
UFISestuaryFamilyAfamily = prune_samples(sample_sums(UFISestuaryFamilyAfamily)>=1, UFISestuaryFamilyAfamily)
nsamples(UFISestuaryFamilyAfamily) #90

#Make a matrix with the sample_data and relabel subcategories to make metadata.defined column#
UF_IS_Est_family <- as.matrix(sample_data(UFISestuaryFamilyAfamily))
UF_IS_Est_family[UF_IS_Est_family ==  "SLE"] <- 2
UF_IS_Est_family[UF_IS_Est_family ==  "IRL"] <- 1
UF_IS_Est_family<-as.data.frame(UF_IS_Est_family)

#Create Bray.Curtis distance matrix from family table#
bc <- vegdist(t(otu_table(UFISestuaryFamilyAfamily)))

#Conduct PAM clustering based upon 2 clusters#
pam.res <- pam(bc,2, diss=TRUE)

#Add a microbially.predicted column to matrix above for comparison to metadata.defined column#
UF_IS_Est_family <- cbind(UF_IS_Est_family, MP_UFISEstfamily = pam.res$cluster)

#Test Total Estuary Efficiency# 1
sum(UF_IS_Est_family$Estuary == UF_IS_Est_family$MP_UFISEstfamily)/nrow(UF_IS_Est_family)

#Test SLE Efficiency aka Sensitvity# 1
UF_IS_SLE_family <- filter(UF_IS_Est_family, Estuary == "2")
sum(UF_IS_SLE_family$Estuary == UF_IS_SLE_family$MP_UFISEstfamily)/nrow(UF_IS_SLE_family)

#Test IRL Efficiency aka Specificity# 1
UF_IS_IRL_family <- filter(UF_IS_Est_family, Estuary == "1")
sum(UF_IS_IRL_family$Estuary == UF_IS_IRL_family$MP_UFISEstfamily)/nrow(UF_IS_IRL_family)

TotalPCoA <- ordinate(UFISestuaryFamilyAfamily,"PCoA")
p = plot_ordination(UFISestuaryFamilyAfamily, TotalPCoA, color="Estuary") + 
  ggtitle("Unfiltered IS Estuary family Indicators PCoA with Bray.Curtis distance") +
  geom_text(aes(label=UF_IS_Est_family$MP_UFISEstfamily),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p


##Test filtered Estuary familys##

#Remove any samples with no familys from target phyloseq#
nsamples(FISestuaryFamilyAfamily) #90
FISestuaryFamilyAfamily = prune_samples(sample_sums(FISestuaryFamilyAfamily)>=1, FISestuaryFamilyAfamily)
nsamples(FISestuaryFamilyAfamily) #90

#Make a matrix with the sample_data and relabel subcategories to make metadata.defined column#
F_IS_Est_family <- as.matrix(sample_data(FISestuaryFamilyAfamily))
F_IS_Est_family[F_IS_Est_family ==  "SLE"] <- 2
F_IS_Est_family[F_IS_Est_family ==  "IRL"] <- 1
F_IS_Est_family<-as.data.frame(F_IS_Est_family)


#Create Bray.Curtis distance matrix from family table#
bc <- vegdist(t(otu_table(FISestuaryFamilyAfamily)))

#Conduct PAM clustering based upon 2 clusters#
pam.res <- pam(bc,2, diss=TRUE)

#Add a microbially.predicted column to matrix above for comparison to metadata.defined column#
F_IS_Est_family <- cbind(F_IS_Est_family, MP_FISEstfamily = pam.res$cluster)

#Test Total Estuary Efficiency# 0.9777778
sum(F_IS_Est_family$Estuary == F_IS_Est_family$MP_FISEstfamily)/nrow(F_IS_Est_family)
#Test SLE Efficiency aka Sensitivity# 1
F_IS_SLE_family <- filter(F_IS_Est_family, Estuary == "2")
sum(F_IS_SLE_family$Estuary == F_IS_SLE_family$MP_FISEstfamily)/nrow(F_IS_SLE_family)
#Test IRL Efficiency aka Specificity# 0.969697
F_IS_IRL_family <- filter(F_IS_Est_family, Estuary == "1")
sum(F_IS_IRL_family$Estuary == F_IS_IRL_family$MP_FISEstfamily)/nrow(F_IS_IRL_family)

TotalPCoA <- ordinate(FISestuaryFamilyAfamily,"PCoA")
p = plot_ordination(FISestuaryFamilyAfamily, TotalPCoA, color="Estuary") + 
  ggtitle("Filtered IS Estuary family Indicators PCoA with Bray.Curtis distance") +
  geom_text(aes(label=F_IS_Est_family$MP_FISEstfamily),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p


###TEST DESEQ2 ESTUARY INDICATORS###

##Test unfiltered Estuary familys##

#Remove any samples with no familys from target phyloseq#
nsamples(UFDestuaryFamilyAfamily) #90
UFDestuaryFamilyAfamily = prune_samples(sample_sums(UFDestuaryFamilyAfamily)>=1, UFDestuaryFamilyAfamily)
nsamples(UFDestuaryFamilyAfamily) #90

#Make a matrix with the sample_data and relabel subcategories to make metadata.defined column#
UF_D_Est_family <- as.matrix(sample_data(UFDestuaryFamilyAfamily))
UF_D_Est_family[UF_D_Est_family ==  "SLE"] <- 2
UF_D_Est_family[UF_D_Est_family ==  "IRL"] <- 1
UF_D_Est_family<-as.data.frame(UF_D_Est_family)


#Create Bray.Curtis distance matrix from family table#
bc <- vegdist(t(otu_table(UFDestuaryFamilyAfamily)))

#Conduct PAM clustering based upon 2 clusters#
pam.res <- pam(bc,2, diss=TRUE)

#Add a microbially.predicted column to matrix above for comparison to metadata.defined column#
UF_D_Est_family <- cbind(UF_D_Est_family, MP_UFDEstfamily = pam.res$cluster)

#Test Total Estuary Efficiency# 0.9666667
sum(UF_D_Est_family$Estuary == UF_D_Est_family$MP_UFDEstfamily)/nrow(UF_D_Est_family)

#Test SLE Efficiency aka Sensitivity# 0.9583333
UF_D_SLE_family <- filter(UF_D_Est_family, Estuary == "2")
sum(UF_D_SLE_family$Estuary == UF_D_SLE_family$MP_UFDEstfamily)/nrow(UF_D_SLE_family)

#Test IRL Efficiency aka Specificity# 0.969697
UF_D_IRL_family <- filter(UF_D_Est_family, Estuary == "1")
sum(UF_D_IRL_family$Estuary == UF_D_IRL_family$MP_UFDEstfamily)/nrow(UF_D_IRL_family)


TotalPCoA <- ordinate(UFDestuaryFamilyAfamily,"PCoA")
p = plot_ordination(UFDestuaryFamilyAfamily, TotalPCoA, color="Estuary") + 
  ggtitle("Unfiltered DESeq2 Estuary family Indicators PCoA with Bray.Curtis distance") +
  geom_text(aes(label=UF_D_Est_family$MP_UFDEstfamily),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p


#Test filtered DESeq2 Estuary familys#

#Remove any samples with no familys from target phyloseq#
nsamples(FDestuaryFamilyAfamily) #90
FDestuaryFamilyAfamily = prune_samples(sample_sums(FDestuaryFamilyAfamily)>=1, FDestuaryFamilyAfamily)
nsamples(FDestuaryFamilyAfamily) #90

#Make a matrix with the sample_data and relabel subcategories to make metadata.defined column#
F_D_Est_family <- as.matrix(sample_data(FDestuaryFamilyAfamily))
F_D_Est_family[F_D_Est_family ==  "SLE"] <- 2
F_D_Est_family[F_D_Est_family ==  "IRL"] <- 1
F_D_Est_family<-as.data.frame(F_D_Est_family)


#Create Bray.Curtis distance matrix from family table#
bc <- vegdist(t(otu_table(FDestuaryFamilyAfamily)))

#Conduct PAM clustering based upon 2 clusters#
pam.res <- pam(bc,2, diss=TRUE)

#Add a microbially.predicted column to matrix above for comparison to metadata.defined column#
F_D_Est_family <- cbind(F_D_Est_family, MP_FDEstfamily = pam.res$cluster)

#Test Total Estuary Efficiency# 0.6
sum(F_D_Est_family$Estuary == F_D_Est_family$MP_FDEstfamily)/nrow(F_D_Est_family)
#Test SLE Efficiency aka Sensitivity# 0.7916667
F_D_SLE_family <- filter(F_D_Est_family, Estuary == "2")
sum(F_D_SLE_family$Estuary == F_D_SLE_family$MP_FDEstfamily)/nrow(F_D_SLE_family)
#Test IRL Efficiency aka Specificity# 0.530303
F_D_IRL_family <- filter(F_D_Est_family, Estuary == "1")
sum(F_D_IRL_family$Estuary == F_D_IRL_family$MP_FDEstfamily)/nrow(F_D_IRL_family)

TotalPCoA <- ordinate(FDestuaryFamilyAfamily,"PCoA")
p = plot_ordination(FDestuaryFamilyAfamily, TotalPCoA, color="Estuary") + 
  ggtitle("Filtered DESeq2 Estuary family Indicators PCoA with Bray.Curtis distance") +
  geom_text(aes(label=F_D_Est_family$MP_FDEstfamily),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p

###TEST COMBO ESTUARY INDICATORS###

##Test unfiltered Estuary familys##

#Remove any samples with no familys from target phyloseq#
nsamples(UFCestuaryFamilyAfamily) #90
UFCestuaryFamilyAfamily = prune_samples(sample_sums(UFCestuaryFamilyAfamily)>=1, UFCestuaryFamilyAfamily)
nsamples(UFCestuaryFamilyAfamily) #90

#Make a matrix with the sample_data and relabel subcategories to make metadata.defined column#
UF_C_Est_family <- as.matrix(sample_data(UFCestuaryFamilyAfamily))
UF_C_Est_family[UF_C_Est_family ==  "SLE"] <- 2
UF_C_Est_family[UF_C_Est_family ==  "IRL"] <- 1
UF_C_Est_family<-as.data.frame(UF_C_Est_family)


#Create Bray.Curtis distance matrix from family table#
bc <- vegdist(t(otu_table(UFCestuaryFamilyAfamily)))

#Conduct PAM clustering based upon 2 clusters#
pam.res <- pam(bc,2, diss=TRUE)

#Add a microbially.predicted column to matrix above for comparison to metadata.defined column#
UF_C_Est_family <- cbind(UF_C_Est_family, MP_UFCEstfamily = pam.res$cluster)

#Test Total Estuary Efficiency# 1
sum(UF_C_Est_family$Estuary == UF_C_Est_family$MP_UFCEstfamily)/nrow(UF_C_Est_family)

#Test SLE Efficiency aka Sensitivity# 1
UF_C_SLE_family <- filter(UF_C_Est_family, Estuary == "2")
sum(UF_C_SLE_family$Estuary == UF_C_SLE_family$MP_UFCEstfamily)/nrow(UF_C_SLE_family)

#Test IRL Efficiency aka Specificity# 1
UF_C_IRL_family <- filter(UF_C_Est_family, Estuary == "1")
sum(UF_C_IRL_family$Estuary == UF_C_IRL_family$MP_UFCEstfamily)/nrow(UF_C_IRL_family)

TotalPCoA <- ordinate(UFCestuaryFamilyAfamily,"PCoA")
p = plot_ordination(UFCestuaryFamilyAfamily, TotalPCoA, color="Estuary") + 
  ggtitle("Unfiltered Final Estuary family Indicators PCoA with Bray.Curtis distance") +
  geom_text(aes(label=UF_C_Est_family$MP_UFCEstfamily),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p


##Test filtered Estuary familys##

#Remove any samples with no familys from target phyloseq#
nsamples(FCestuaryFamilyAfamily) #90
FCestuaryFamilyAfamily = prune_samples(sample_sums(FCestuaryFamilyAfamily)>=1, FCestuaryFamilyAfamily)
nsamples(FCestuaryFamilyAfamily) #90

#Make a matrix with the sample_data and relabel subcategories to make metadata.defined column#
F_C_Est_family <- as.matrix(sample_data(FCestuaryFamilyAfamily))
F_C_Est_family[F_C_Est_family ==  "SLE"] <- 2
F_C_Est_family[F_C_Est_family ==  "IRL"] <- 1
F_C_Est_family<-as.data.frame(F_C_Est_family)


#Create Bray.Curtis distance matrix from family table#
bc <- vegdist(t(otu_table(FCestuaryFamilyAfamily)))

#Conduct PAM clustering based upon 2 clusters#
pam.res <- pam(bc,2, diss=TRUE)

#Add a microbially.predicted column to matrix above for comparison to metadata.defined column#
F_C_Est_family <- cbind(F_C_Est_family, MP_FCEstfamily = pam.res$cluster)

#Test Total Estuary Efficiency# 0.9666667
sum(F_C_Est_family$Estuary == F_C_Est_family$MP_FCEstfamily)/nrow(F_C_Est_family)

#Test SLE Efficiency aka Sensitivity# 0.875
F_C_SLE_family <- filter(F_C_Est_family, Estuary == "2")
sum(F_C_SLE_family$Estuary == F_C_SLE_family$MP_FCEstfamily)/nrow(F_C_SLE_family)

#Test IRL Efficiency aka Sensitivity# 1
F_C_IRL_family <- filter(F_C_Est_family, Estuary == "1")
sum(F_C_IRL_family$Estuary == F_C_IRL_family$MP_FCEstfamily)/nrow(F_C_IRL_family)

TotalPCoA <- ordinate(FCestuaryFamilyAfamily,"PCoA")
p = plot_ordination(FCestuaryFamilyAfamily, TotalPCoA, color="Estuary") + 
  ggtitle("Filtered Final Estuary family Indicators PCoA with Bray.Curtis distance") +
  geom_text(aes(label=F_C_Est_family$MP_FCEstfamily),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p


###TEST INDICSPECIES MUCK INDICATORS###

##Test unfiltered Muck familys##

#Remove any samples with no familys from target phyloseq#
nsamples(UFISmuckFamilyAfamily) #90
UFISmuckFamilyAfamily = prune_samples(sample_sums(UFISmuckFamilyAfamily)>=1, UFISmuckFamilyAfamily)
nsamples(UFISmuckFamilyAfamily) #90

#Make a matrix with the sample_data and relabel subcategories to make metadata.defined column#
UF_IS_Muck_family <- as.matrix(sample_data(UFISmuckFamilyAfamily))
UF_IS_Muck_family[UF_IS_Muck_family ==  "Muck"] <- 2
UF_IS_Muck_family[UF_IS_Muck_family ==  "Not"] <- 1
UF_IS_Muck_family[UF_IS_Muck_family ==  "Mucky"] <- 2
UF_IS_Muck_family[UF_IS_Muck_family ==  "Muckish"] <- 2
UF_IS_Muck_family<-as.data.frame(UF_IS_Muck_family)


#Create Bray.Curtis distance matrix from family table#
bc <- vegdist(t(otu_table(UFISmuckFamilyAfamily)))

#Conduct PAM clustering based upon 2 clusters#
pam.res <- pam(bc,2, diss=TRUE)

#Add a microbially.predicted column to matrix above for comparison to metadata.defined column#
UF_IS_Muck_family <- cbind(UF_IS_Muck_family, MP_UFISMuckfamily = pam.res$cluster)

#Test Total Muck Efficiency# 0.8333333
sum(UF_IS_Muck_family$Muck == UF_IS_Muck_family$MP_UFISMuckfamily)/nrow(UF_IS_Muck_family)

#Test 3MC Efficiency aka Sensitivity# 0.8333333
UF_IS_3MC_family <- filter(UF_IS_Muck_family, Muck == "2")
sum(UF_IS_3MC_family$Muck == UF_IS_3MC_family$MP_UFISMuckfamily)/nrow(UF_IS_3MC_family)

#Test 0MC Efficiency aka Specificity# 0.8333333
UF_IS_0MC_family <- filter(UF_IS_Muck_family, Muck == "1")
sum(UF_IS_0MC_family$Muck == UF_IS_0MC_family$MP_UFISMuckfamily)/nrow(UF_IS_0MC_family)


TotalPCoA <- ordinate(UFISmuckFamilyAfamily,"PCoA")
p = plot_ordination(UFISmuckFamilyAfamily, TotalPCoA, color="Muck") + 
  ggtitle("Unfiltered IS Muck family Indicators PCoA with Bray.Curtis distance") +
  geom_text(aes(label=UF_IS_Muck_family$MP_UFISMuckfamily),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p


##Test filtered Muck familys##

#Remove any samples with no familys from target phyloseq#
nsamples(FISmuckFamilyAfamily) #90
FISmuckFamilyAfamily = prune_samples(sample_sums(FISmuckFamilyAfamily)>=1, FISmuckFamilyAfamily)
nsamples(FISmuckFamilyAfamily) #90

#Make a matrix with the sample_data and relabel subcategories to make metadata.defined column#
F_IS_Muck_family <- as.matrix(sample_data(FISmuckFamilyAfamily))
F_IS_Muck_family[F_IS_Muck_family ==  "Muck"] <- 2
F_IS_Muck_family[F_IS_Muck_family ==  "Not"] <- 1
F_IS_Muck_family[F_IS_Muck_family ==  "Mucky"] <- 2
F_IS_Muck_family[F_IS_Muck_family ==  "Muckish"] <- 2
F_IS_Muck_family<-as.data.frame(F_IS_Muck_family)


#Create Bray.Curtis distance matrix from family table#
bc <- vegdist(t(otu_table(FISmuckFamilyAfamily)))

#Conduct PAM clustering based upon 2 clusters#
pam.res <- pam(bc,2, diss=TRUE)

#Add a microbially.predicted column to matrix above for comparison to metadata.defined column#
F_IS_Muck_family <- cbind(F_IS_Muck_family, MP_FISMuckfamily = pam.res$cluster)

#Test Total Muck Efficiency# 0.9111111
sum(F_IS_Muck_family$Muck == F_IS_Muck_family$MP_FISMuckfamily)/nrow(F_IS_Muck_family)
#Test 3MC Efficiency aka  Sensitivity# 0.8666667
F_IS_3MC_family <- filter(F_IS_Muck_family, Muck == "2")
sum(F_IS_3MC_family$Muck == F_IS_3MC_family$MP_FISMuckfamily)/nrow(F_IS_3MC_family)
#Test 0MC Efficiency aka Specificity# 0.9333333
F_IS_0MC_family <- filter(F_IS_Muck_family, Muck == "1")
sum(F_IS_0MC_family$Muck == F_IS_0MC_family$MP_FISMuckfamily)/nrow(F_IS_0MC_family)

TotalPCoA <- ordinate(FISmuckFamilyAfamily,"PCoA")
p = plot_ordination(FISmuckFamilyAfamily, TotalPCoA, color="Muck") + 
  ggtitle("Filtered IS Muck family Indicators PCoA with Bray.Curtis distance") +
  geom_text(aes(label=F_IS_Muck_family$MP_FISMuckfamily),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p

###TEST DESEQ2 MUCK INDICATORS###

#Test unfiltered Muck familys##

#Remove any samples with no familys from target phyloseq#
nsamples(UFDmuckFamilyAfamily) #90
UFDmuckFamilyAfamily = prune_samples(sample_sums(UFDmuckFamilyAfamily)>=1, UFDmuckFamilyAfamily)
nsamples(UFDmuckFamilyAfamily) #90

#Make a matrix with the sample_data and relabel subcategories to make metadata.defined column#
UF_D_Muck_family <- as.matrix(sample_data(UFDmuckFamilyAfamily))
UF_D_Muck_family[UF_D_Muck_family ==  "Muck"] <- 2
UF_D_Muck_family[UF_D_Muck_family ==  "Not"] <- 1
UF_D_Muck_family[UF_D_Muck_family ==  "Mucky"] <- 2
UF_D_Muck_family[UF_D_Muck_family ==  "Muckish"] <- 2
UF_D_Muck_family<-as.data.frame(UF_D_Muck_family)


#Create Bray.Curtis distance matrix from family table#
bc <- vegdist(t(otu_table(UFDmuckFamilyAfamily)))

#Conduct PAM clustering based upon 2 clusters#
pam.res <- pam(bc,2, diss=TRUE)

#Add a microbially.predicted column to matrix above for comparison to metadata.defined column#
UF_D_Muck_family <- cbind(UF_D_Muck_family, MP_UFDMuckfamily = pam.res$cluster)

#Test Total Muck Efficiency# 0.8555556
sum(UF_D_Muck_family$Muck == UF_D_Muck_family$MP_UFDMuckfamily)/nrow(UF_D_Muck_family)

#Test 3MC Efficiency aka Sensitivity# 0.7666667
UF_D_3MC_family <- filter(UF_D_Muck_family, Muck == "2")
sum(UF_D_3MC_family$Muck == UF_D_3MC_family$MP_UFDMuckfamily)/nrow(UF_D_3MC_family)

#Test 0MC Efficiency aka Sensitivity# 0.9
UF_D_0MC_family <- filter(UF_D_Muck_family, Muck == "1")
sum(UF_D_0MC_family$Muck == UF_D_0MC_family$MP_UFDMuckfamily)/nrow(UF_D_0MC_family)

TotalPCoA <- ordinate(UFDmuckFamilyAfamily,"PCoA")
p = plot_ordination(UFDmuckFamilyAfamily, TotalPCoA, color="Muck") + 
  ggtitle("Unfiltered DESeq2 Muck family Indicators PCoA with Bray.Curtis distance") +
  geom_text(aes(label=UF_D_Muck_family$MP_UFDMuckfamily),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p


##Test filtered Muck familys##

#Remove any samples with no familys from target phyloseq#
nsamples(FDmuckFamilyAfamily) #90
FDmuckFamilyAfamily = prune_samples(sample_sums(FDmuckFamilyAfamily)>=1, FDmuckFamilyAfamily)
nsamples(FDmuckFamilyAfamily) #90

#Make a matrix with the sample_data and relabel subcategories to make metadata.defined column#
F_D_Muck_family <- as.matrix(sample_data(FDmuckFamilyAfamily))
F_D_Muck_family[F_D_Muck_family ==  "Muck"] <- 2
F_D_Muck_family[F_D_Muck_family ==  "Not"] <- 1
F_D_Muck_family[F_D_Muck_family ==  "Mucky"] <- 2
F_D_Muck_family[F_D_Muck_family ==  "Muckish"] <- 2
F_D_Muck_family<-as.data.frame(F_D_Muck_family)


#Create Bray.Curtis distance matrix from family table#
bc <- vegdist(t(otu_table(FDmuckFamilyAfamily)))

#Conduct PAM clustering based upon 2 clusters#
pam.res <- pam(bc,2, diss=TRUE)

#Add a microbially.predicted column to matrix above for comparison to metadata.defined column#
F_D_Muck_family <- cbind(F_D_Muck_family, MP_FDMuckfamily = pam.res$cluster)

#Test Total Muck Efficiency# 0.9222222
sum(F_D_Muck_family$Muck == F_D_Muck_family$MP_FDMuckfamily)/nrow(F_D_Muck_family)

#Test 3MC Efficiency aka Sensitivy# 0.9166667
F_D_3MC_family <- filter(F_D_Muck_family, Muck == "1")
sum(F_D_3MC_family$Muck == F_D_3MC_family$MP_FDMuckfamily)/nrow(F_D_3MC_family)

#Test 0MC Efficiency aka Specificity# 0.9333333
F_D_0MC_family <- filter(F_D_Muck_family, Muck == "2")
sum(F_D_0MC_family$Muck == F_D_0MC_family$MP_FDMuckfamily)/nrow(F_D_0MC_family)

TotalPCoA <- ordinate(FDmuckFamilyAfamily,"PCoA")
p = plot_ordination(FDmuckFamilyAfamily, TotalPCoA, color="Muck") + 
  ggtitle("Filtered DESeq2 Muck family Indicators PCoA with Bray.Curtis distance") +
  geom_text(aes(label=F_D_Muck_family$MP_FDMuckfamily),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p

###TEST COMBO MUCK INDICATORS###

##Test unfiltered Muck familys##

#Remove any samples with no familys from target phyloseq#
nsamples(UFCmuckFamilyAfamily) #90
UFCmuckFamilyAfamily = prune_samples(sample_sums(UFCmuckFamilyAfamily)>=1, UFCmuckFamilyAfamily)
nsamples(UFCmuckFamilyAfamily) #90

#Make a matrix with the sample_data and relabel subcategories to make metadata.defined column#
UF_C_Muck_familys <- as.matrix(sample_data(UFCmuckFamilyAfamily))
UF_C_Muck_familys[UF_C_Muck_familys ==  "Muck"] <- 2
UF_C_Muck_familys[UF_C_Muck_familys ==  "Not"] <- 1
UF_C_Muck_familys[UF_C_Muck_familys ==  "Mucky"] <- 2
UF_C_Muck_familys[UF_C_Muck_familys ==  "Muckish"] <- 2
UF_C_Muck_familys<-as.data.frame(UF_C_Muck_familys)


#Create Bray.Curtis distance matrix from family table#
bc <- vegdist(t(otu_table(UFCmuckFamilyAfamily)))

#Conduct PAM clustering based upon 2 clusters#
pam.res <- pam(bc,2, diss=TRUE)

#Add a microbially.predicted column to matrix above for comparison to metadata.defined column#
UF_C_Muck_familys <- cbind(UF_C_Muck_familys, MP_UFCMuckfamily = pam.res$cluster)

#Test Total Muck Efficiency# 0.9
sum(UF_C_Muck_familys$Muck == UF_C_Muck_familys$MP_UFCMuckfamily)/nrow(UF_C_Muck_familys)

#Test Muck Efficiency aka Sensitivity# 0.9
F_C_3MC_familys <- filter(UF_C_Muck_familys, Muck == "2")
sum(F_C_3MC_familys$Muck == F_C_3MC_familys$MP_UFCMuckfamily)/nrow(F_C_3MC_familys)

#Test Muck Efficiency aka Specificity# 0.9
F_C_0MC_familys <- filter(UF_C_Muck_familys, Muck == "1")
sum(F_C_0MC_familys$Muck == F_C_0MC_familys$MP_UFCMuckfamily)/nrow(F_C_0MC_familys)

TotalPCoA <- ordinate(UFCmuckFamilyAfamily,"PCoA")
p = plot_ordination(UFCmuckFamilyAfamily, TotalPCoA, color="Muck") + 
  ggtitle("Unfiltered Final Muck family Indicators PCoA with Bray.Curtis distance") +
  geom_text(aes(label=UF_C_Muck_familys$MP_UFCMuckfamily),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p

##Test filtered Muck familys##

#Remove any samples with no familys from target phyloseq#
nsamples(FCmuckFamilyAfamily) #90
FCmuckFamilyAfamily = prune_samples(sample_sums(FCmuckFamilyAfamily)>=1, FCmuckFamilyAfamily)
nsamples(FCmuckFamilyAfamily) #90

#Make a matrix with the sample_data and relabel subcategories to make metadata.defined column#
F_C_Muck_familys <- as.matrix(sample_data(FCmuckFamilyAfamily))
F_C_Muck_familys[F_C_Muck_familys ==  "Muck"] <- 2
F_C_Muck_familys[F_C_Muck_familys ==  "Not"] <- 1
F_C_Muck_familys[F_C_Muck_familys ==  "Mucky"] <- 2
F_C_Muck_familys[F_C_Muck_familys ==  "Muckish"] <- 2
F_C_Muck_familys<-as.data.frame(F_C_Muck_familys)


#Create Bray.Curtis distance matrix from family table#
bc <- vegdist(t(otu_table(FCmuckFamilyAfamily)))

#Conduct PAM clustering based upon 2 clusters#
pam.res <- pam(bc,2, diss=TRUE)

#Add a microbially.predicted column to matrix above for comparison to metadata.defined column#
F_C_Muck_familys <- cbind(F_C_Muck_familys, MP_FCMuckfamily = pam.res$cluster)

#Test Total Muck Efficiency# 0.8555556
sum(F_C_Muck_familys$Muck == F_C_Muck_familys$MP_FCMuckfamily)/nrow(F_C_Muck_familys)

#Test Muck Efficiency aka Sensitivity# 0.9
F_C_3MC_family <- filter(F_C_Muck_familys, Muck == "2")
sum(F_C_3MC_family$Muck == F_C_3MC_family$MP_FCMuckfamily)/nrow(F_C_3MC_family)

#Test Not Efficiency aka Specificity# 0.8333333
F_C_0MC_family <- filter(F_C_Muck_familys, Muck == "1")
sum(F_C_0MC_family$Muck == F_C_0MC_family$MP_FCMuckfamily)/nrow(F_C_0MC_family)

TotalPCoA <- ordinate(FCmuckFamilyAfamily,"PCoA")
p = plot_ordination(FCmuckFamilyAfamily, TotalPCoA, color="Muck") + 
  ggtitle("Filtered Final Muck family Indicators PCoA with Bray.Curtis distance") +
  geom_text(aes(label=F_C_Muck_familys$MP_FCMuckfamily),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p

###TEST INDICSPECIES TOM/Cu INDICATORS###

##Test unfiltered TOM/Cu familys##

#Create phyloseq object focused on HiHi and HiLo samples#
UFIShihiFamilyAfamily = subset_samples(UFIStomcuFamilyAfamily, LOI.Cu=="High.High")
UFIShiloFamilyAfamily = subset_samples(UFIStomcuFamilyAfamily, LOI.Cu=="High.Low")
UFISfoctomcuFamilyAfamily = merge_phyloseq(UFIShihiFamilyAfamily, UFIShiloFamilyAfamily)

#Remove any samples with no familys from target phyloseq#
nsamples(UFISfoctomcuFamilyAfamily) #29
UFISfoctomcuFamilyAfamily = prune_samples(sample_sums(UFISfoctomcuFamilyAfamily)>=1, UFISfoctomcuFamilyAfamily)
nsamples(UFISfoctomcuFamilyAfamily) #29

#Make a matrix with the sample_data and relabel subcategories to make metadata.defined column#
UF_IS_TOMCu_family <- as.matrix(sample_data(UFISfoctomcuFamilyAfamily))
UF_IS_TOMCu_family[UF_IS_TOMCu_family ==  "High.High"] <- 1
UF_IS_TOMCu_family[UF_IS_TOMCu_family ==  "High.Low"] <- 2
UF_IS_TOMCu_family<-as.data.frame(UF_IS_TOMCu_family)

#Create Bray.Curtis distance matrix from family table#
bc <- vegdist(t(otu_table(UFISfoctomcuFamilyAfamily)))

#Conduct PAM clustering based upon 2 clusters#
pam.res <- pam(bc,2, diss=TRUE)

#Add a microbially.predicted column to matrix above for comparison to metadata.defined column#
UF_IS_TOMCu_family <- cbind(UF_IS_TOMCu_family, MP_UFISTOMCufamily = pam.res$cluster)

#Test Total TOM/Cu Efficiency# 0.8275862
sum(UF_IS_TOMCu_family$LOI.Cu == UF_IS_TOMCu_family$MP_UFISTOMCufamily)/nrow(UF_IS_TOMCu_family)

#Test High.High Indicators aka Sensitivity# 1
UF_IS_HiHi_family <- filter(UF_IS_TOMCu_family, LOI.Cu == "1")
sum(UF_IS_HiHi_family$LOI.Cu == UF_IS_HiHi_family$MP_UFISTOMCufamily)/nrow(UF_IS_HiHi_family)

#Test High.Low Indicators aka Specificity# 0.7222222
UF_IS_HiLo_family <- filter(UF_IS_TOMCu_family, LOI.Cu == "2")
sum(UF_IS_HiLo_family$LOI.Cu == UF_IS_HiLo_family$MP_UFISTOMCufamily)/nrow(UF_IS_HiLo_family)

TotalPCoA <- ordinate(UFISfoctomcuFamilyAfamily,"PCoA")
p = plot_ordination(UFISfoctomcuFamilyAfamily, TotalPCoA, color="LOI.Cu") + 
  ggtitle("Unfiltered IS LOICu family Indicators PCoA with Bray.Curtis distance") +
  geom_text(aes(label=UF_IS_TOMCu_family$MP_UFISTOMCufamily),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p


##Test filtered TOM/Cu familys##

#Create phyloseq object focused on HiHi and HiLo samples#
FIShihiFamilyAfamily = subset_samples(FIStomcuFamilyAfamily, LOI.Cu=="High.High")
FIShiloFamilyAfamily = subset_samples(FIStomcuFamilyAfamily, LOI.Cu=="High.Low")
FISfoctomcuFamilyAfamily = merge_phyloseq(FIShihiFamilyAfamily, FIShiloFamilyAfamily)

#Remove any samples with no familys from target phyloseq#
nsamples(FISfoctomcuFamilyAfamily) #29
FISfoctomcuFamilyAfamily = prune_samples(sample_sums(FISfoctomcuFamilyAfamily)>=1, FISfoctomcuFamilyAfamily)
nsamples(FISfoctomcuFamilyAfamily) #29

#Make a matrix with the sample_data and relabel subcategories to make metadata.defined column#
F_IS_TOMCu_family <- as.matrix(sample_data(FISfoctomcuFamilyAfamily))
F_IS_TOMCu_family[F_IS_TOMCu_family ==  "High.High"] <- 1
F_IS_TOMCu_family[F_IS_TOMCu_family ==  "High.Low"] <- 2
F_IS_TOMCu_family<-as.data.frame(F_IS_TOMCu_family)


#Create Bray.Curtis distance matrix from family table#
bc <- vegdist(t(otu_table(FISfoctomcuFamilyAfamily)))

#Conduct PAM clustering based upon 2 clusters#
pam.res <- pam(bc,2, diss=TRUE)

#Add a microbially.predicted column to matrix above for comparison to metadata.defined column#
F_IS_TOMCu_family <- cbind(F_IS_TOMCu_family, MP_FISTOMCufamily = pam.res$cluster)

#Test Total TOM/Cu Efficiency# 1
sum(F_IS_TOMCu_family$LOI.Cu == F_IS_TOMCu_family$MP_FISTOMCufamily)/nrow(F_IS_TOMCu_family)

#Test HiHi Efficiency aka Sensitivity# 1
F_IS_HiHi_family <- filter(F_IS_TOMCu_family, LOI.Cu == "1")
sum(F_IS_HiHi_family$LOI.Cu == F_IS_HiHi_family$MP_FISTOMCufamily)/nrow(F_IS_HiHi_family)

#Test HiLo Efficiency aka Sensitivity# 1
F_IS_HiLo_family <- filter(F_IS_TOMCu_family, LOI.Cu == "2") 
sum(F_IS_HiLo_family$LOI.Cu == F_IS_HiLo_family$MP_FISTOMCufamily)/nrow(F_IS_HiLo_family)

TotalPCoA <- ordinate(FISfoctomcuFamilyAfamily,"PCoA")
p = plot_ordination(FISfoctomcuFamilyAfamily, TotalPCoA, color="LOI.Cu") + 
  ggtitle("Filtered IS LOICu family Indicators PCoA with Bray.Curtis distance") +
  geom_text(aes(label=F_IS_TOMCu_family$MP_FISTOMCufamily),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p

##TEST DESEQ2 LOICu INDICATORS##

##Test unfiltered TOM/Cu familys##

#Create phyloseq object focused on HiHi and HiLo samples#
UFDhihiFamilyAfamily = subset_samples(UFDtomcuFamilyAfamily, LOI.Cu=="High.High")
UFDhiloFamilyAfamily = subset_samples(UFDtomcuFamilyAfamily, LOI.Cu=="High.Low")
UFDfoctomcuFamilyAfamily = merge_phyloseq(UFDhihiFamilyAfamily, UFDhiloFamilyAfamily)

#Remove any samples with no familys from target phyloseq#
nsamples(UFDfoctomcuFamilyAfamily) #29
UFDfoctomcuFamilyAfamily = prune_samples(sample_sums(UFDfoctomcuFamilyAfamily)>=1, UFDfoctomcuFamilyAfamily)
nsamples(UFDfoctomcuFamilyAfamily) #29

#Make a matrix with the sample_data and relabel subcategories to make metadata.defined column#
UF_D_TOMCu_family <- as.matrix(sample_data(UFDfoctomcuFamilyAfamily))
UF_D_TOMCu_family[UF_D_TOMCu_family ==  "High.High"] <- 1
UF_D_TOMCu_family[UF_D_TOMCu_family ==  "High.Low"] <- 2
UF_D_TOMCu_family<-as.data.frame(UF_D_TOMCu_family)


#Create Bray.Curtis distance matrix from family table#
bc <- vegdist(t(otu_table(UFDfoctomcuFamilyAfamily)))

#Conduct PAM clustering based upon 2 clusters#
pam.res <- pam(bc,2, diss=TRUE)

#Add a microbially.predicted column to matrix above for comparison to metadata.defined column#
UF_D_TOMCu_family <- cbind(UF_D_TOMCu_family, MP_UFDTOMCufamily = pam.res$cluster)

#Test Total TOM/Cu Efficiency# 0.7586207
sum(UF_D_TOMCu_family$LOI.Cu == UF_D_TOMCu_family$MP_UFDTOMCufamily)/nrow(UF_D_TOMCu_family)

#Test HiHi Efficiency aka Sensitivity# 0.8181818
UF_D_HiHi_family <- filter(UF_D_TOMCu_family, LOI.Cu == "1")
sum(UF_D_HiHi_family$LOI.Cu == UF_D_HiHi_family$MP_UFDTOMCufamily)/nrow(UF_D_HiHi_family)

#Test HiLo Efficiency aka Specificity# 0.7222222
UF_D_HiLo_family <- filter(UF_D_TOMCu_family, LOI.Cu == "2")
sum(UF_D_HiLo_family$LOI.Cu == UF_D_HiLo_family$MP_UFDTOMCufamily)/nrow(UF_D_HiLo_family)

TotalPCoA <- ordinate(UFDfoctomcuFamilyAfamily,"PCoA")
p = plot_ordination(UFDfoctomcuFamilyAfamily, TotalPCoA, color="LOI.Cu") + 
  ggtitle("Unfiltered DESeq2 LOICu family Indicators PCoA with Bray.Curtis distance") +
  geom_text(aes(label=UF_D_TOMCu_family$MP_UFDTOMCufamily),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p


##Test filtered TOM/Cu familys##

#Create phyloseq object focused on HiHi and HiLo samples#
FDhihiFamilyAfamily = subset_samples(FDtomcuFamilyAfamily, LOI.Cu=="High.High")
FDhiloFamilyAfamily = subset_samples(FDtomcuFamilyAfamily, LOI.Cu=="High.Low")
FDfoctomcuFamilyAfamily = merge_phyloseq(FDhihiFamilyAfamily, FDhiloFamilyAfamily)

#Remove any samples with no familys from target phyloseq#
nsamples(FDfoctomcuFamilyAfamily) #29
FDfoctomcuFamilyAfamily = prune_samples(sample_sums(FDfoctomcuFamilyAfamily)>=1, FDfoctomcuFamilyAfamily)
nsamples(FDfoctomcuFamilyAfamily) #29


#Make a matrix with the sample_data and relabel subcategories to make metadata.defined column#
F_D_TOMCu_family <- as.matrix(sample_data(FDfoctomcuFamilyAfamily))
F_D_TOMCu_family[F_D_TOMCu_family ==  "High.High"] <- 2
F_D_TOMCu_family[F_D_TOMCu_family ==  "High.Low"] <- 1
F_D_TOMCu_family<-as.data.frame(F_D_TOMCu_family)


#Create Bray.Curtis distance matrix from family table#
bc <- vegdist(t(otu_table(FDfoctomcuFamilyAfamily)))

#Conduct PAM clustering based upon 2 clusters#
pam.res <- pam(bc,2, diss=TRUE)

#Add a microbially.predicted column to matrix above for comparison to metadata.defined column#
F_D_TOMCu_family <- cbind(F_D_TOMCu_family, MP_FDTOMCufamily = pam.res$cluster)

#Test Total TOM/Cu Efficiency# 0.7931034
sum(F_D_TOMCu_family$LOI.Cu == F_D_TOMCu_family$MP_FDTOMCufamily)/nrow(F_D_TOMCu_family)

#Test HiHi Efficiency aka Sensitivity# 0.4545455
F_D_HiHi_family <- filter(F_D_TOMCu_family, LOI.Cu == "2")
sum(F_D_HiHi_family$LOI.Cu == F_D_HiHi_family$MP_FDTOMCufamily)/nrow(F_D_HiHi_family)

#Test HiLo Efficiency aka Specificity# 1
F_D_HiLo_family <- filter(F_D_TOMCu_family, LOI.Cu == "1")
sum(F_D_HiLo_family$LOI.Cu == F_D_HiLo_family$MP_FDTOMCufamily)/nrow(F_D_HiLo_family)

TotalPCoA <- ordinate(FDfoctomcuFamilyAfamily,"PCoA")
p = plot_ordination(FDfoctomcuFamilyAfamily, TotalPCoA, color="LOI.Cu") + 
  ggtitle("Filtered DESeq2 LOICu family Indicators PCoA with Bray.Curtis distance") +
  geom_text(aes(label=F_D_TOMCu_family$MP_FDTOMCufamily),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p

##TEST COMBO TOM/CU INDICATORS##

##Test unfiltered TOM/Cu familys##

#Create phyloseq object focused on HiHi and HiLo samples#
UFChihiFamilyAfamily = subset_samples(UFCtomcuFamilyAfamily, LOI.Cu=="High.High")
UFChiloFamilyAfamily = subset_samples(UFCtomcuFamilyAfamily, LOI.Cu=="High.Low")
UFCfoctomcuFamilyAfamily = merge_phyloseq(UFChihiFamilyAfamily, UFChiloFamilyAfamily)

#Remove any samples with no familys from target phyloseq#
nsamples(UFCfoctomcuFamilyAfamily) #29
UFCfoctomcuFamilyAfamily = prune_samples(sample_sums(UFCfoctomcuFamilyAfamily)>=1, UFCfoctomcuFamilyAfamily)
nsamples(UFCfoctomcuFamilyAfamily) #29

#Make a matrix with the sample_data and relabel subcategories to make metadata.defined column#
UF_C_TOMCu_family <- as.matrix(sample_data(UFCfoctomcuFamilyAfamily))
UF_C_TOMCu_family[UF_C_TOMCu_family ==  "High.High"] <- 1
UF_C_TOMCu_family[UF_C_TOMCu_family ==  "High.Low"] <- 2
UF_C_TOMCu_family<-as.data.frame(UF_C_TOMCu_family)


#Create Bray.Curtis distance matrix from family table#
bc <- vegdist(t(otu_table(UFCfoctomcuFamilyAfamily)))

#Conduct PAM clustering based upon 2 clusters#
pam.res <- pam(bc,2, diss=TRUE)

#Add a microbially.predicted column to matrix above for comparison to metadata.defined column#
UF_C_TOMCu_family <- cbind(UF_C_TOMCu_family, MP_UFCTOMCufamily = pam.res$cluster)

#Test Total TOM/Cu Efficiency# 0.8275862
sum(UF_C_TOMCu_family$LOI.Cu == UF_C_TOMCu_family$MP_UFCTOMCufamily)/nrow(UF_C_TOMCu_family)

#Test HiHi Efficiency aka Sensitivity# 1
UF_C_HiHi_family <- filter(UF_C_TOMCu_family, LOI.Cu == "1")
sum(UF_C_HiHi_family$LOI.Cu == UF_C_HiHi_family$MP_UFCTOMCufamily)/nrow(UF_C_HiHi_family)

#Test HiLo Efficiency aka Specificity# 0.7222222
UF_C_HiLo_family <- filter(UF_C_TOMCu_family, LOI.Cu == "2")
sum(UF_C_HiLo_family$LOI.Cu == UF_C_HiLo_family$MP_UFCTOMCufamily)/nrow(UF_C_HiLo_family)

TotalPCoA <- ordinate(UFCfoctomcuFamilyAfamily,"PCoA")
p = plot_ordination(UFCfoctomcuFamilyAfamily, TotalPCoA, color="LOI.Cu") + 
  ggtitle("Unfiltered Final LOICu family Indicators PCoA with Bray.Curtis distance") +
  geom_text(aes(label=UF_C_TOMCu_family$MP_UFCTOMCufamily),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p


##Test filtered TOM/Cu familys##

#Create phyloseq object focused on HiHi and HiLo samples#
FChihiFamilyAfamily = subset_samples(FCtomcuFamilyAfamily, LOI.Cu=="High.High")
FChiloFamilyAfamily = subset_samples(FCtomcuFamilyAfamily, LOI.Cu=="High.Low")
FCfoctomcuFamilyAfamily = merge_phyloseq(FChihiFamilyAfamily, FChiloFamilyAfamily)

#Remove any samples with no familys from target phyloseq#
nsamples(FCfoctomcuFamilyAfamily) #29
FCfoctomcuFamilyAfamily = prune_samples(sample_sums(FCfoctomcuFamilyAfamily)>=1, FCfoctomcuFamilyAfamily)
nsamples(FCfoctomcuFamilyAfamily) #29

#Make a matrix with the sample_data and relabel subcategories to make metadata.defined column#
F_C_TOMCu_family <- as.matrix(sample_data(FCfoctomcuFamilyAfamily))
F_C_TOMCu_family[F_C_TOMCu_family ==  "High.High"] <- 2
F_C_TOMCu_family[F_C_TOMCu_family ==  "High.Low"] <- 1
F_C_TOMCu_family<-as.data.frame(F_C_TOMCu_family)


#Create Bray.Curtis distance matrix from family table#
bc <- vegdist(t(otu_table(FCfoctomcuFamilyAfamily)))

#Conduct PAM clustering based upon 2 clusters#
pam.res <- pam(bc,2, diss=TRUE)

#Add a microbially.predicted column to matrix above for comparison to metadata.defined column#
F_C_TOMCu_family <- cbind(F_C_TOMCu_family, MP_FCTOMCufamily = pam.res$cluster)

#Test Total TOM/Cu Efficiency# 0.862069
sum(F_C_TOMCu_family$LOI.Cu == F_C_TOMCu_family$MP_FCTOMCufamily)/nrow(F_C_TOMCu_family)

#Test HiHi Efficiency aka Sensitivity# 0.8181818
F_C_HiHi_family <- filter(F_C_TOMCu_family, LOI.Cu == "2")
sum(F_C_HiHi_family$LOI.Cu == F_C_HiHi_family$MP_FCTOMCufamily)/nrow(F_C_HiHi_family)

#Test HiLo Efficiency aka Specificity# 0.8888889
F_C_HiLo_family <- filter(F_C_TOMCu_family, LOI.Cu == "1")
sum(F_C_HiLo_family$LOI.Cu == F_C_HiLo_family$MP_FCTOMCufamily)/nrow(F_C_HiLo_family)


TotalPCoA <- ordinate(FCfoctomcuFamilyAfamily,"PCoA")
p = plot_ordination(FCfoctomcuFamilyAfamily, TotalPCoA, color="LOI.Cu") + 
  ggtitle("Filtered Final LOICu family Indicators PCoA with Bray.Curtis distance") +
  geom_text(aes(label=F_C_TOMCu_family$MP_FCTOMCufamily),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p