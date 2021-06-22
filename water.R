###Getting your indicator lists###

#Create a phyloseq object of your target data. Here we are using the water dataset from a 11 sites in the Indian River Lagoon


#Filter out zero and low abundance asvs



#Determine the number of samples in each category and subcategory you are testing#
length(which(sample_data(RFLWSW_Y1_DataAasv)$Estuary == "SLE")) #24
length(which(sample_data(RFLWSW_Y1_DataAasv)$Estuary == "IRL")) #42

###RUN DESEQ2 AT asv LEVEL ON ESTUARY SUBCATEGORIES###

#Convert phyloseq to deseq2 object centered around a irlvssle factor#
#W_estuary = metadata category testing
#asv = taxonomic level
#cudds = DESeq2
#W 



W_estuary_cudds_asv = phyloseq_to_deseq2(RFLWSW_Y1_DataAasv,  ~ Estuary)

#Calculate geometric means prior to estimate size factors#
W_estuary_geoMeans_asv = apply(counts(W_estuary_cudds_asv), 1, gm_mean)
W_estuary_cudds_asv = estimateSizeFactors(W_estuary_cudds_asv, geoMeans = W_estuary_geoMeans_asv)

#Conduct DESEQ2 test#
W_estuary_cudds_asv = DESeq(W_estuary_cudds_asv, fitType="local")

#Explore the results#
W_estuary_DESeq2_res_asv = results(W_estuary_cudds_asv)
W_estuary_DESeq2_res_asv = W_estuary_DESeq2_res_asv[order(W_estuary_DESeq2_res_asv$padj, na.last=NA), ]
alpha = 0.05
W_estuary_DESeq2_sig_res_asv = W_estuary_DESeq2_res_asv[(W_estuary_DESeq2_res_asv$padj < alpha), ]

#Make dataframe with taxanomy added in#
W_estuary_DESeq2_sig_res_taxo_asv = cbind(as(W_estuary_DESeq2_sig_res_asv, "data.frame"), as(tax_table(RFLWSW_Y1_DataAasv)[rownames(W_estuary_DESeq2_sig_res_asv), ], "matrix"))

#Make dataframe with asvs added in from all sites#
W_estuary_DESeq2_sig_res_taxo_seqs_asv = cbind(as(W_estuary_DESeq2_sig_res_taxo_asv, "data.frame"), as(otu_table(RFLWSW_Y1_DataAasv)[rownames(W_estuary_DESeq2_sig_res_taxo_asv), ], "matrix"))

#Make rownames an actual column and remove old rownames#
W_estuary_DESeq2_sig_res_taxo_seqs_asv <- cbind(ESV.ID = rownames(W_estuary_DESeq2_sig_res_taxo_seqs_asv), W_estuary_DESeq2_sig_res_taxo_seqs_asv)
rownames(W_estuary_DESeq2_sig_res_taxo_seqs_asv) <-NULL

#Make a txt file of results for all sites#
write.csv(as.data.frame(W_estuary_DESeq2_sig_res_taxo_seqs_asv), file="W_estuary_DESeq2_sig_res_taxo_seqs_asv.csv")

#Detemine which subcategory is negative or positive#
W_estuary_DESeq2_res_asv
#W_estuary_DESeq2_res_asv: SLE (+) vs IRL (-)#

#Determine number of Indicators with padj >0.05 in each subcagetory#
#Number of SLE indicators# 
length(which(W_estuary_DESeq2_sig_res_taxo_seqs_asv$log2FoldChange > 0)) #359
#Number of IRL indicators#
length(which(W_estuary_DESeq2_sig_res_taxo_seqs_asv$log2FoldChange < 0)) #301

#Order rows by log2FoldChange#
W_estuary_DESeq2_sig_res_taxo_seqs_asv <- arrange(W_estuary_DESeq2_sig_res_taxo_seqs_asv, log2FoldChange)

#Create lists of just significant asvs for each category#
unfiltered_DESeq2_W_estuary_asvs <- subset(W_estuary_DESeq2_sig_res_taxo_seqs_asv, select=c(ESV.ID))

##Create phyloseq objects for DESeq2 indicators##

#Create phyloseq objects of the unfiltered indicators using character string version of significant asvs##
charac_unfiltered_DESeq2_W_estuary_asvs <- unfiltered_DESeq2_W_estuary_asvs[,"ESV.ID"]
charac_unfiltered_DESeq2_W_estuary_asvs <- as.character(charac_unfiltered_DESeq2_W_estuary_asvs)
UFDW_estuaryDataAasv <- prune_taxa(charac_unfiltered_DESeq2_W_estuary_asvs, RFLWSW_Y1_DataAasv)

unfiltered_DESeq2_W_estuary_asv_heatmap <- plot_heatmap(UFDW_estuaryDataAasv, taxa.order = charac_unfiltered_DESeq2_W_estuary_asvs, sample.order= "Estuary", sample.label = "Estuary")
unfiltered_DESeq2_W_estuary_asv_heatmap

###CONDUCT INDICSPECIES ANALYSIS AT asv LEVEL###

#Increase max print so that all of results can be printed#
options(max.print=100000000)


##Conduct asv indicspecies analysis for Estuary Category#

#Take out OTU aka asv table from metadata focused phyloseq object#
W_estuary_seqs_asv = as(otu_table(RFLWSW_Y1_DataAasv), "matrix")

#Make it a data frame, transpose it, make it a data frame again#
W_estuary_seqs_asv <- as.data.frame(W_estuary_seqs_asv)
W_estuary_seqs_asv<- t(W_estuary_seqs_asv)
W_estuary_seqs_asv <- as.data.frame(W_estuary_seqs_asv)

#Take out metadata table metadata focused phyloseq object#
W_estuary_meta = as(sample_data(RFLWSW_Y1_DataAasv), "matrix")

#Keep only the sampleids and metadata category focusing on#
W_estuary_meta <- subset(W_estuary_meta, select=c(sampleid, Estuary))

#Reset names of two factors#
W_estuary_meta[W_estuary_meta == "IRL"] <- 1
W_estuary_meta[W_estuary_meta == "SLE"] <- 2

#Convert it to data frame, then create class object based on metadatacategory#
W_estuary_meta = as.data.frame(W_estuary_meta)
W_estuary_meta_group <- as.character(W_estuary_meta$Estuary)

#Run indicspecies#
W_estuary_indval_asv = multipatt(W_estuary_seqs_asv, W_estuary_meta_group, control = how(nperm=9999))

#Create a blank csv, then sink all indicspecies results into it#
sink("W_estuary_IS_res_asv.csv")
W_estuary_sig_indval_asv <- summary(W_estuary_indval_asv, indvalcomp=TRUE, alpha=1)
sink()




#Reset max print to default#
options(max.print=99999)

##Look at csv file you sinked into and create a table
##of just ASV.ID, A, B, stat, p.value, and group##

##Continue asv indicspecies analysis for Estuary Category##

#Upload new table#
W_estuary_IS_res_asv <- read.csv(file.choose())
#Make a new p value adjusted column#
W_estuary_IS_res_asv$padj <- p.adjust(W_estuary_IS_res_asv$p.value, method = "BH")

#Keep table of adjusted p values less than 0.05#
W_estuary_IS_sig_res_asv <- filter(W_estuary_IS_res_asv, padj < 0.05)

#Make dataframe with taxanomy added in by making ESV.ID the rownames then adding taxonomy table#
W_estuary_IS_sig_res_taxo_asv <- as.data.frame(W_estuary_IS_sig_res_asv)
rownames(W_estuary_IS_sig_res_taxo_asv) <- W_estuary_IS_sig_res_taxo_asv[,1]
W_estuary_IS_sig_res_taxo_asv = cbind(as(W_estuary_IS_sig_res_taxo_asv, "data.frame"), as(tax_table(RFLWSW_Y1_DataAasv)[rownames(W_estuary_IS_sig_res_taxo_asv), ], "matrix"))

#Make dataframe with asvs added in from all sites#
W_estuary_IS_sig_res_taxo_seqs_asv = cbind(as(W_estuary_IS_sig_res_taxo_asv, "data.frame"), as(otu_table(RFLWSW_Y1_DataAasv)[rownames(W_estuary_IS_sig_res_taxo_asv), ], "matrix"))

#Remove old rownames#
rownames(W_estuary_IS_sig_res_taxo_seqs_asv) <-NULL

#Make a txt file of results for all sites#
write.csv(as.data.frame(W_estuary_IS_sig_res_taxo_seqs_asv), file="W_estuary_IS_sig_res_taxo_seqs_asv.csv")

#Check the number of indicators in each group after padj#
length(which(sample_data(W_estuary_IS_sig_res_asv)$group == "1")) #336 IRL
length(which(sample_data(W_estuary_IS_sig_res_asv)$group == "2")) #796 SLE

#Order rows by group#
W_estuary_IS_sig_res_taxo_seqs_asv <- arrange(W_estuary_IS_sig_res_taxo_seqs_asv, group)

#Create lists of just significant asvs for each category#
unfiltered_IS_W_estuary_asvs <- subset(W_estuary_IS_sig_res_taxo_seqs_asv, select=c(ESV.ID))

#Create phyloseq objects of the unfiltered indicators using character string version of significant asvs##
charac_unfiltered_IS_W_estuary_asvs <- unfiltered_IS_W_estuary_asvs[,"ESV.ID"]
charac_unfiltered_IS_W_estuary_asvs <- as.character(charac_unfiltered_IS_W_estuary_asvs)
UFISW_estuaryDataAasv <- prune_taxa(charac_unfiltered_IS_W_estuary_asvs, RFLWSW_Y1_DataAasv)

##Create heatmaps to see overall trends, may take a while to load##

unfiltered_IS_W_estuary_asv_heatmap <- plot_heatmap(UFDW_estuaryDataAasv, taxa.order = charac_unfiltered_IS_W_estuary_asvs, sample.order= "Estuary", sample.label = "Estuary")
unfiltered_IS_W_estuary_asv_heatmap

###MAKE COMBO INDICATOR LISTS###

##Make "Combo" indicator lists from overlapping unfiltered DESeq2 and IS indicators##

#Keep rows in the W_estuary DESEQ object that are also in the IRLSLE IS object#
unfiltered_W_estuary_combo_sig_res_taxo_seqs_asv <- merge(W_estuary_IS_sig_res_asv,W_estuary_DESeq2_sig_res_taxo_seqs_asv, by="ESV.ID")

#Determine the number of indicators in each subcategory of the unfiltered Combo tables#
length(which(sample_data(unfiltered_W_estuary_combo_sig_res_taxo_seqs_asv)$group == "1")) #242 IRL
length(which(sample_data(unfiltered_W_estuary_combo_sig_res_taxo_seqs_asv)$group == "2")) #356 SLE

#Order unfiltered Combo table rows by group#
unfiltered_W_estuary_combo_sig_res_taxo_seqs_asv <- arrange(unfiltered_W_estuary_combo_sig_res_taxo_seqs_asv, group)


#Make a csv file for each of the unfiltered Combo tables#
write.csv(as.data.frame(unfiltered_W_estuary_combo_sig_res_taxo_seqs_asv), file="unfiltered_W_estuary_combo_sig_res_taxo_seqs_asv.csv")

#Create lists of unfiltered significant asvs for each category#
unfiltered_Combo_W_estuary_asvs <- subset(unfiltered_W_estuary_combo_sig_res_taxo_seqs_asv, select=c(ESV.ID))

##Create phyloseq objects for Combo indicators##

#Create phyloseq objects of the unfiltered indicators using character string version of significant asvs##
charac_unfiltered_Combo_W_estuary_asvs <- unfiltered_Combo_W_estuary_asvs[,"ESV.ID"]
charac_unfiltered_Combo_W_estuary_asvs <- as.character(charac_unfiltered_Combo_W_estuary_asvs)
UFCW_estuaryDataAasv <- prune_taxa(charac_unfiltered_Combo_W_estuary_asvs, RFLWSW_Y1_DataAasv)

##Create heatmaps to see overall trends, may take a while to load##

unfiltered_Combo_W_estuary_asv_heatmap <- plot_heatmap(UFDW_estuaryDataAasv, taxa.order = charac_unfiltered_Combo_W_estuary_asvs, sample.order= "Estuary", sample.label = "Estuary")
unfiltered_Combo_W_estuary_asv_heatmap


##Species##

#Determine the number of samples in each category and subcategory you are testing#
length(which(sample_data(RFLWSW_Y1_DataAasv)$Estuary == "SLE")) #24
length(which(sample_data(RFLWSW_Y1_DataAasv)$Estuary == "IRL")) #42

###RUN DESEQ2 AT asv LEVEL ON ESTUARY SUBCATEGORIES###

#Convert phyloseq to deseq2 object centered around a irlvssle factor#
#W_estuary = metadata category testing
#asv = taxonomic level
#cudds = DESeq2
#W 



W_estuary_cudds_asv = phyloseq_to_deseq2(RFLWSW_Y1_DataAasv,  ~ Estuary)

#Calculate geometric means prior to estimate size factors#
W_estuary_geoMeans_asv = apply(counts(W_estuary_cudds_asv), 1, gm_mean)
W_estuary_cudds_asv = estimateSizeFactors(W_estuary_cudds_asv, geoMeans = W_estuary_geoMeans_asv)

#Conduct DESEQ2 test#
W_estuary_cudds_asv = DESeq(W_estuary_cudds_asv, fitType="local")

#Explore the results#
W_estuary_DESeq2_res_asv = results(W_estuary_cudds_asv)
W_estuary_DESeq2_res_asv = W_estuary_DESeq2_res_asv[order(W_estuary_DESeq2_res_asv$padj, na.last=NA), ]
alpha = 0.05
W_estuary_DESeq2_sig_res_asv = W_estuary_DESeq2_res_asv[(W_estuary_DESeq2_res_asv$padj < alpha), ]

#Make dataframe with taxanomy added in#
W_estuary_DESeq2_sig_res_taxo_asv = cbind(as(W_estuary_DESeq2_sig_res_asv, "data.frame"), as(tax_table(RFLWSW_Y1_DataAasv)[rownames(W_estuary_DESeq2_sig_res_asv), ], "matrix"))

#Make dataframe with asvs added in from all sites#
W_estuary_DESeq2_sig_res_taxo_seqs_asv = cbind(as(W_estuary_DESeq2_sig_res_taxo_asv, "data.frame"), as(otu_table(RFLWSW_Y1_DataAasv)[rownames(W_estuary_DESeq2_sig_res_taxo_asv), ], "matrix"))

#Make rownames an actual column and remove old rownames#
W_estuary_DESeq2_sig_res_taxo_seqs_asv <- cbind(ESV.ID = rownames(W_estuary_DESeq2_sig_res_taxo_seqs_asv), W_estuary_DESeq2_sig_res_taxo_seqs_asv)
rownames(W_estuary_DESeq2_sig_res_taxo_seqs_asv) <-NULL

#Make a txt file of results for all sites#
write.csv(as.data.frame(W_estuary_DESeq2_sig_res_taxo_seqs_asv), file="W_estuary_DESeq2_sig_res_taxo_seqs_asv.csv")

#Detemine which subcategory is negative or positive#
W_estuary_DESeq2_res_asv
#W_estuary_DESeq2_res_asv: SLE (+) vs IRL (-)#

#Determine number of Indicators with padj >0.05 in each subcagetory#
#Number of SLE indicators# 
length(which(W_estuary_DESeq2_sig_res_taxo_seqs_asv$log2FoldChange > 0)) #359
#Number of IRL indicators#
length(which(W_estuary_DESeq2_sig_res_taxo_seqs_asv$log2FoldChange < 0)) #301

#Order rows by log2FoldChange#
W_estuary_DESeq2_sig_res_taxo_seqs_asv <- arrange(W_estuary_DESeq2_sig_res_taxo_seqs_asv, log2FoldChange)

#Create lists of just significant asvs for each category#
unfiltered_DESeq2_W_estuary_asvs <- subset(W_estuary_DESeq2_sig_res_taxo_seqs_asv, select=c(ESV.ID))

##Create phyloseq objects for DESeq2 indicators##

#Create phyloseq objects of the unfiltered indicators using character string version of significant asvs##
charac_unfiltered_DESeq2_W_estuary_asvs <- unfiltered_DESeq2_W_estuary_asvs[,"ESV.ID"]
charac_unfiltered_DESeq2_W_estuary_asvs <- as.character(charac_unfiltered_DESeq2_W_estuary_asvs)
UFDW_estuaryDataAasv <- prune_taxa(charac_unfiltered_DESeq2_W_estuary_asvs, RFLWSW_Y1_DataAasv)

unfiltered_DESeq2_W_estuary_asv_heatmap <- plot_heatmap(UFDW_estuaryDataAasv, taxa.order = charac_unfiltered_DESeq2_W_estuary_asvs, sample.order= "Estuary", sample.label = "Estuary")
unfiltered_DESeq2_W_estuary_asv_heatmap

###CONDUCT INDICSPECIES ANALYSIS AT asv LEVEL###

#Increase max print so that all of results can be printed#
options(max.print=100000000)


##Conduct asv indicspecies analysis for Estuary Category#

#Take out OTU aka asv table from metadata focused phyloseq object#
W_estuary_seqs_asv = as(otu_table(RFLWSW_Y1_DataAasv), "matrix")

#Make it a data frame, transpose it, make it a data frame again#
W_estuary_seqs_asv <- as.data.frame(W_estuary_seqs_asv)
W_estuary_seqs_asv<- t(W_estuary_seqs_asv)
W_estuary_seqs_asv <- as.data.frame(W_estuary_seqs_asv)

#Take out metadata table metadata focused phyloseq object#
W_estuary_meta = as(sample_data(RFLWSW_Y1_DataAasv), "matrix")

#Keep only the sampleids and metadata category focusing on#
W_estuary_meta <- subset(W_estuary_meta, select=c(sampleid, Estuary))

#Reset names of two factors#
W_estuary_meta[W_estuary_meta == "IRL"] <- 1
W_estuary_meta[W_estuary_meta == "SLE"] <- 2

#Convert it to data frame, then create class object based on metadatacategory#
W_estuary_meta = as.data.frame(W_estuary_meta)
W_estuary_meta_group <- as.character(W_estuary_meta$Estuary)

#Run indicspecies#
W_estuary_indval_asv = multipatt(W_estuary_seqs_asv, W_estuary_meta_group, control = how(nperm=9999))

#Create a blank csv, then sink all indicspecies results into it#
sink("W_estuary_IS_res_asv.csv")
W_estuary_sig_indval_asv <- summary(W_estuary_indval_asv, indvalcomp=TRUE, alpha=1)
sink()




#Reset max print to default#
options(max.print=99999)

##Look at csv file you sinked into and create a table
##of just ASV.ID, A, B, stat, p.value, and group##

##Continue asv indicspecies analysis for Estuary Category##

#Upload new table#
W_estuary_IS_res_asv <- read.csv(file.choose())
#Make a new p value adjusted column#
W_estuary_IS_res_asv$padj <- p.adjust(W_estuary_IS_res_asv$p.value, method = "BH")

#Keep table of adjusted p values less than 0.05#
W_estuary_IS_sig_res_asv <- filter(W_estuary_IS_res_asv, padj < 0.05)

#Make dataframe with taxanomy added in by making ESV.ID the rownames then adding taxonomy table#
W_estuary_IS_sig_res_taxo_asv <- as.data.frame(W_estuary_IS_sig_res_asv)
rownames(W_estuary_IS_sig_res_taxo_asv) <- W_estuary_IS_sig_res_taxo_asv[,1]
W_estuary_IS_sig_res_taxo_asv = cbind(as(W_estuary_IS_sig_res_taxo_asv, "data.frame"), as(tax_table(RFLWSW_Y1_DataAasv)[rownames(W_estuary_IS_sig_res_taxo_asv), ], "matrix"))

#Make dataframe with asvs added in from all sites#
W_estuary_IS_sig_res_taxo_seqs_asv = cbind(as(W_estuary_IS_sig_res_taxo_asv, "data.frame"), as(otu_table(RFLWSW_Y1_DataAasv)[rownames(W_estuary_IS_sig_res_taxo_asv), ], "matrix"))

#Remove old rownames#
rownames(W_estuary_IS_sig_res_taxo_seqs_asv) <-NULL

#Make a txt file of results for all sites#
write.csv(as.data.frame(W_estuary_IS_sig_res_taxo_seqs_asv), file="W_estuary_IS_sig_res_taxo_seqs_asv.csv")

#Check the number of indicators in each group after padj#
length(which(sample_data(W_estuary_IS_sig_res_asv)$group == "1")) #336 IRL
length(which(sample_data(W_estuary_IS_sig_res_asv)$group == "2")) #796 SLE

#Order rows by group#
W_estuary_IS_sig_res_taxo_seqs_asv <- arrange(W_estuary_IS_sig_res_taxo_seqs_asv, group)

#Create lists of just significant asvs for each category#
unfiltered_IS_W_estuary_asvs <- subset(W_estuary_IS_sig_res_taxo_seqs_asv, select=c(ESV.ID))

#Create phyloseq objects of the unfiltered indicators using character string version of significant asvs##
charac_unfiltered_IS_W_estuary_asvs <- unfiltered_IS_W_estuary_asvs[,"ESV.ID"]
charac_unfiltered_IS_W_estuary_asvs <- as.character(charac_unfiltered_IS_W_estuary_asvs)
UFISW_estuaryDataAasv <- prune_taxa(charac_unfiltered_IS_W_estuary_asvs, RFLWSW_Y1_DataAasv)

##Create heatmaps to see overall trends, may take a while to load##

unfiltered_IS_W_estuary_asv_heatmap <- plot_heatmap(UFDW_estuaryDataAasv, taxa.order = charac_unfiltered_IS_W_estuary_asvs, sample.order= "Estuary", sample.label = "Estuary")
unfiltered_IS_W_estuary_asv_heatmap

###MAKE COMBO INDICATOR LISTS###

##Make "Combo" indicator lists from overlapping unfiltered DESeq2 and IS indicators##

#Keep rows in the W_estuary DESEQ object that are also in the IRLSLE IS object#
unfiltered_W_estuary_combo_sig_res_taxo_seqs_asv <- merge(W_estuary_IS_sig_res_asv,W_estuary_DESeq2_sig_res_taxo_seqs_asv, by="ESV.ID")

#Determine the number of indicators in each subcategory of the unfiltered Combo tables#
length(which(sample_data(unfiltered_W_estuary_combo_sig_res_taxo_seqs_asv)$group == "1")) #242 IRL
length(which(sample_data(unfiltered_W_estuary_combo_sig_res_taxo_seqs_asv)$group == "2")) #356 SLE

#Order unfiltered Combo table rows by group#
unfiltered_W_estuary_combo_sig_res_taxo_seqs_asv <- arrange(unfiltered_W_estuary_combo_sig_res_taxo_seqs_asv, group)


#Make a csv file for each of the unfiltered Combo tables#
write.csv(as.data.frame(unfiltered_W_estuary_combo_sig_res_taxo_seqs_asv), file="unfiltered_W_estuary_combo_sig_res_taxo_seqs_asv.csv")

#Create lists of unfiltered significant asvs for each category#
unfiltered_Combo_W_estuary_asvs <- subset(unfiltered_W_estuary_combo_sig_res_taxo_seqs_asv, select=c(ESV.ID))

##Create phyloseq objects for Combo indicators##

#Create phyloseq objects of the unfiltered indicators using character string version of significant asvs##
charac_unfiltered_Combo_W_estuary_asvs <- unfiltered_Combo_W_estuary_asvs[,"ESV.ID"]
charac_unfiltered_Combo_W_estuary_asvs <- as.character(charac_unfiltered_Combo_W_estuary_asvs)
UFCW_estuaryDataAasv <- prune_taxa(charac_unfiltered_Combo_W_estuary_asvs, RFLWSW_Y1_DataAasv)

##Create heatmaps to see overall trends, may take a while to load##

unfiltered_Combo_W_estuary_asv_heatmap <- plot_heatmap(UFDW_estuaryDataAasv, taxa.order = charac_unfiltered_Combo_W_estuary_asvs, sample.order= "Estuary", sample.label = "Estuary")
unfiltered_Combo_W_estuary_asv_heatmap

###SPECIES###

#Determine the number of samples in each category and subcategory you are testing#
length(which(sample_data(RFLWSW_Y1_SpeciesAspecies)$Estuary == "SLE")) #24
length(which(sample_data(RFLWSW_Y1_SpeciesAspecies)$Estuary == "IRL")) #42

###RUN DESEQ2 AT species LEVEL ON ESTUARY SUBCATEGORIES###

#Convert phyloseq to deseq2 object centered around a irlvssle factor#
#W_estuary = metadata category testing
#species = taxonomic level
#cudds = DESeq2
#W 



W_estuary_cudds_species = phyloseq_to_deseq2(RFLWSW_Y1_SpeciesAspecies,  ~ Estuary)

#Calculate geometric means prior to estimate size factors#
W_estuary_geoMeans_species = apply(counts(W_estuary_cudds_species), 1, gm_mean)
W_estuary_cudds_species = estimateSizeFactors(W_estuary_cudds_species, geoMeans = W_estuary_geoMeans_species)

#Conduct DESEQ2 test#
W_estuary_cudds_species = DESeq(W_estuary_cudds_species, fitType="local")

#Explore the results#
W_estuary_DESeq2_res_species = results(W_estuary_cudds_species)
W_estuary_DESeq2_res_species = W_estuary_DESeq2_res_species[order(W_estuary_DESeq2_res_species$padj, na.last=NA), ]
alpha = 0.05
W_estuary_DESeq2_sig_res_species = W_estuary_DESeq2_res_species[(W_estuary_DESeq2_res_species$padj < alpha), ]

#Make dataframe with taxanomy added in#
W_estuary_DESeq2_sig_res_taxo_species = cbind(as(W_estuary_DESeq2_sig_res_species, "data.frame"), as(tax_table(RFLWSW_Y1_SpeciesAspecies)[rownames(W_estuary_DESeq2_sig_res_species), ], "matrix"))

#Make dataframe with speciess added in from all sites#
W_estuary_DESeq2_sig_res_taxo_seqs_species = cbind(as(W_estuary_DESeq2_sig_res_taxo_species, "data.frame"), as(otu_table(RFLWSW_Y1_SpeciesAspecies)[rownames(W_estuary_DESeq2_sig_res_taxo_species), ], "matrix"))

#Make rownames an actual column and remove old rownames#
W_estuary_DESeq2_sig_res_taxo_seqs_species <- cbind(ESV.ID = rownames(W_estuary_DESeq2_sig_res_taxo_seqs_species), W_estuary_DESeq2_sig_res_taxo_seqs_species)
rownames(W_estuary_DESeq2_sig_res_taxo_seqs_species) <-NULL

#Make a txt file of results for all sites#
write.csv(as.data.frame(W_estuary_DESeq2_sig_res_taxo_seqs_species), file="W_estuary_DESeq2_sig_res_taxo_seqs_species.csv")

#Detemine which subcategory is negative or positive#
W_estuary_DESeq2_res_species
#W_estuary_DESeq2_res_species: SLE (+) vs IRL (-)#

#Determine number of Indicators with padj >0.05 in each subcagetory#
#Number of SLE indicators# 
length(which(W_estuary_DESeq2_sig_res_taxo_seqs_species$log2FoldChange > 0)) #177
#Number of IRL indicators#
length(which(W_estuary_DESeq2_sig_res_taxo_seqs_species$log2FoldChange < 0)) #168

#Order rows by log2FoldChange#
W_estuary_DESeq2_sig_res_taxo_seqs_species <- arrange(W_estuary_DESeq2_sig_res_taxo_seqs_species, log2FoldChange)

#Create lists of just significant speciess for each category#
unfiltered_DESeq2_W_estuary_speciess <- subset(W_estuary_DESeq2_sig_res_taxo_seqs_species, select=c(ESV.ID))

##Create phyloseq objects for DESeq2 indicators##

#Create phyloseq objects of the unfiltered indicators using character string version of significant speciess##
charac_unfiltered_DESeq2_W_estuary_speciess <- unfiltered_DESeq2_W_estuary_speciess[,"ESV.ID"]
charac_unfiltered_DESeq2_W_estuary_speciess <- as.character(charac_unfiltered_DESeq2_W_estuary_speciess)
UFDW_estuaryDataAspecies <- prune_taxa(charac_unfiltered_DESeq2_W_estuary_speciess, RFLWSW_Y1_SpeciesAspecies)

unfiltered_DESeq2_W_estuary_species_heatmap <- plot_heatmap(UFDW_estuaryDataAspecies, taxa.order = charac_unfiltered_DESeq2_W_estuary_speciess, sample.order= "Estuary", sample.label = "Estuary")
unfiltered_DESeq2_W_estuary_species_heatmap

###CONDUCT INDICSPECIES ANALYSIS AT species LEVEL###

#Increase max print so that all of results can be printed#
options(max.print=100000000)


##Conduct species indicspecies analysis for Estuary Category#

#Take out OTU aka species table from metadata focused phyloseq object#
W_estuary_seqs_species = as(otu_table(RFLWSW_Y1_SpeciesAspecies), "matrix")

#Make it a data frame, transpose it, make it a data frame again#
W_estuary_seqs_species <- as.data.frame(W_estuary_seqs_species)
W_estuary_seqs_species<- t(W_estuary_seqs_species)
W_estuary_seqs_species <- as.data.frame(W_estuary_seqs_species)

#Run indicspecies#
W_estuary_indval_species = multipatt(W_estuary_seqs_species, W_estuary_meta_group, control = how(nperm=9999))

#Create a blank csv, then sink all indicspecies results into it#
sink("W_estuary_IS_res_species.csv")
W_estuary_sig_indval_species <- summary(W_estuary_indval_species, indvalcomp=TRUE, alpha=1)
sink()




#Reset max print to default#
options(max.print=99999)

##Look at csv file you sinked into and create a table
##of just ASV.ID, A, B, stat, p.value, and group##

##Continue species indicspecies analysis for Estuary Category##

#Upload new table#
W_estuary_IS_res_species <- read.csv(file.choose())
#Make a new p value adjusted column#
W_estuary_IS_res_species$padj <- p.adjust(W_estuary_IS_res_species$p.value, method = "BH")

#Keep table of adjusted p values less than 0.05#
W_estuary_IS_sig_res_species <- filter(W_estuary_IS_res_species, padj < 0.05)

#Make dataframe with taxanomy added in by making ESV.ID the rownames then adding taxonomy table#
W_estuary_IS_sig_res_taxo_species <- as.data.frame(W_estuary_IS_sig_res_species)
rownames(W_estuary_IS_sig_res_taxo_species) <- W_estuary_IS_sig_res_taxo_species[,1]
W_estuary_IS_sig_res_taxo_species = cbind(as(W_estuary_IS_sig_res_taxo_species, "data.frame"), as(tax_table(RFLWSW_Y1_SpeciesAspecies)[rownames(W_estuary_IS_sig_res_taxo_species), ], "matrix"))

#Make dataframe with speciess added in from all sites#
W_estuary_IS_sig_res_taxo_seqs_species = cbind(as(W_estuary_IS_sig_res_taxo_species, "data.frame"), as(otu_table(RFLWSW_Y1_SpeciesAspecies)[rownames(W_estuary_IS_sig_res_taxo_species), ], "matrix"))

#Remove old rownames#
rownames(W_estuary_IS_sig_res_taxo_seqs_species) <-NULL

#Make a txt file of results for all sites#
write.csv(as.data.frame(W_estuary_IS_sig_res_taxo_seqs_species), file="W_estuary_IS_sig_res_taxo_seqs_species.csv")

#Check the number of indicators in each group after padj#
length(which(sample_data(W_estuary_IS_sig_res_species)$group == "1")) #124 IRL
length(which(sample_data(W_estuary_IS_sig_res_species)$group == "2")) #297 SLE

#Order rows by group#
W_estuary_IS_sig_res_taxo_seqs_species <- arrange(W_estuary_IS_sig_res_taxo_seqs_species, group)

#Create lists of just significant speciess for each category#
unfiltered_IS_W_estuary_speciess <- subset(W_estuary_IS_sig_res_taxo_seqs_species, select=c(ESV.ID))

#Create phyloseq objects of the unfiltered indicators using character string version of significant speciess##
charac_unfiltered_IS_W_estuary_speciess <- unfiltered_IS_W_estuary_speciess[,"ESV.ID"]
charac_unfiltered_IS_W_estuary_speciess <- as.character(charac_unfiltered_IS_W_estuary_speciess)
UFISW_estuaryDataAspecies <- prune_taxa(charac_unfiltered_IS_W_estuary_speciess, RFLWSW_Y1_SpeciesAspecies)

##Create heatmaps to see overall trends, may take a while to load##

unfiltered_IS_W_estuary_species_heatmap <- plot_heatmap(UFDW_estuaryDataAspecies, taxa.order = charac_unfiltered_IS_W_estuary_speciess, sample.order= "Estuary", sample.label = "Estuary")
unfiltered_IS_W_estuary_species_heatmap

###MAKE COMBO INDICATOR LISTS###

##Make "Combo" indicator lists from overlapping unfiltered DESeq2 and IS indicators##

#Keep rows in the W_estuary DESEQ object that are also in the IRLSLE IS object#
unfiltered_W_estuary_combo_sig_res_taxo_seqs_species <- merge(W_estuary_IS_sig_res_species,W_estuary_DESeq2_sig_res_taxo_seqs_species, by="ESV.ID")

#Determine the number of indicators in each subcategory of the unfiltered Combo tables#
length(which(sample_data(unfiltered_W_estuary_combo_sig_res_taxo_seqs_species)$group == "1")) #107 IRL
length(which(sample_data(unfiltered_W_estuary_combo_sig_res_taxo_seqs_species)$group == "2")) #168 SLE

#Order unfiltered Combo table rows by group#
unfiltered_W_estuary_combo_sig_res_taxo_seqs_species <- arrange(unfiltered_W_estuary_combo_sig_res_taxo_seqs_species, group)


#Make a csv file for each of the unfiltered Combo tables#
write.csv(as.data.frame(unfiltered_W_estuary_combo_sig_res_taxo_seqs_species), file="unfiltered_W_estuary_combo_sig_res_taxo_seqs_species.csv")

#Create lists of unfiltered significant speciess for each category#
unfiltered_Combo_W_estuary_speciess <- subset(unfiltered_W_estuary_combo_sig_res_taxo_seqs_species, select=c(ESV.ID))

##Create phyloseq objects for Combo indicators##

#Create phyloseq objects of the unfiltered indicators using character string version of significant speciess##
charac_unfiltered_Combo_W_estuary_speciess <- unfiltered_Combo_W_estuary_speciess[,"ESV.ID"]
charac_unfiltered_Combo_W_estuary_speciess <- as.character(charac_unfiltered_Combo_W_estuary_speciess)
UFCW_estuaryDataAspecies <- prune_taxa(charac_unfiltered_Combo_W_estuary_speciess, RFLWSW_Y1_SpeciesAspecies)

##Create heatmaps to see overall trends, may take a while to load##

unfiltered_Combo_W_estuary_species_heatmap <- plot_heatmap(UFDW_estuaryDataAspecies, taxa.order = charac_unfiltered_Combo_W_estuary_speciess, sample.order= "Estuary", sample.label = "Estuary")
unfiltered_Combo_W_estuary_species_heatmap

###GENUS###

#Determine the number of samples in each category and subcategory you are testing#
length(which(sample_data(RFLWSW_Y1_GenusAgenus)$Estuary == "SLE")) #24
length(which(sample_data(RFLWSW_Y1_GenusAgenus)$Estuary == "IRL")) #42

###RUN DESEQ2 AT genus LEVEL ON ESTUARY SUBCATEGORIES###

#Convert phyloseq to deseq2 object centered around a irlvssle factor#
#W_estuary = metadata category testing
#genus = taxonomic level
#cudds = DESeq2
#W 



W_estuary_cudds_genus = phyloseq_to_deseq2(RFLWSW_Y1_GenusAgenus,  ~ Estuary)

#Calculate geometric means prior to estimate size factors#
W_estuary_geoMeans_genus = apply(counts(W_estuary_cudds_genus), 1, gm_mean)
W_estuary_cudds_genus = estimateSizeFactors(W_estuary_cudds_genus, geoMeans = W_estuary_geoMeans_genus)

#Conduct DESEQ2 test#
W_estuary_cudds_genus = DESeq(W_estuary_cudds_genus, fitType="local")

#Explore the results#
W_estuary_DESeq2_res_genus = results(W_estuary_cudds_genus)
W_estuary_DESeq2_res_genus = W_estuary_DESeq2_res_genus[order(W_estuary_DESeq2_res_genus$padj, na.last=NA), ]
alpha = 0.05
W_estuary_DESeq2_sig_res_genus = W_estuary_DESeq2_res_genus[(W_estuary_DESeq2_res_genus$padj < alpha), ]

#Make dataframe with taxanomy added in#
W_estuary_DESeq2_sig_res_taxo_genus = cbind(as(W_estuary_DESeq2_sig_res_genus, "data.frame"), as(tax_table(RFLWSW_Y1_GenusAgenus)[rownames(W_estuary_DESeq2_sig_res_genus), ], "matrix"))

#Make dataframe with genuss added in from all sites#
W_estuary_DESeq2_sig_res_taxo_seqs_genus = cbind(as(W_estuary_DESeq2_sig_res_taxo_genus, "data.frame"), as(otu_table(RFLWSW_Y1_GenusAgenus)[rownames(W_estuary_DESeq2_sig_res_taxo_genus), ], "matrix"))

#Make rownames an actual column and remove old rownames#
W_estuary_DESeq2_sig_res_taxo_seqs_genus <- cbind(ESV.ID = rownames(W_estuary_DESeq2_sig_res_taxo_seqs_genus), W_estuary_DESeq2_sig_res_taxo_seqs_genus)
rownames(W_estuary_DESeq2_sig_res_taxo_seqs_genus) <-NULL

#Make a txt file of results for all sites#
write.csv(as.data.frame(W_estuary_DESeq2_sig_res_taxo_seqs_genus), file="W_estuary_DESeq2_sig_res_taxo_seqs_genus.csv")

#Detemine which subcategory is negative or positive#
W_estuary_DESeq2_res_genus
#W_estuary_DESeq2_res_genus: SLE (+) vs IRL (-)#

#Determine number of Indicators with padj >0.05 in each subcagetory#
#Number of SLE indicators# 
length(which(W_estuary_DESeq2_sig_res_taxo_seqs_genus$log2FoldChange > 0)) #156
#Number of IRL indicators#
length(which(W_estuary_DESeq2_sig_res_taxo_seqs_genus$log2FoldChange < 0)) #143

#Order rows by log2FoldChange#
W_estuary_DESeq2_sig_res_taxo_seqs_genus <- arrange(W_estuary_DESeq2_sig_res_taxo_seqs_genus, log2FoldChange)

#Create lists of just significant genuss for each category#
unfiltered_DESeq2_W_estuary_genuss <- subset(W_estuary_DESeq2_sig_res_taxo_seqs_genus, select=c(ESV.ID))

##Create phyloseq objects for DESeq2 indicators##

#Create phyloseq objects of the unfiltered indicators using character string version of significant genuss##
charac_unfiltered_DESeq2_W_estuary_genuss <- unfiltered_DESeq2_W_estuary_genuss[,"ESV.ID"]
charac_unfiltered_DESeq2_W_estuary_genuss <- as.character(charac_unfiltered_DESeq2_W_estuary_genuss)
UFDW_estuaryDataAgenus <- prune_taxa(charac_unfiltered_DESeq2_W_estuary_genuss, RFLWSW_Y1_GenusAgenus)

unfiltered_DESeq2_W_estuary_genus_heatmap <- plot_heatmap(UFDW_estuaryDataAgenus, taxa.order = charac_unfiltered_DESeq2_W_estuary_genuss, sample.order= "Estuary", sample.label = "Estuary")
unfiltered_DESeq2_W_estuary_genus_heatmap

###CONDUCT INDICSPECIES ANALYSIS AT genus LEVEL###

#Increase max print so that all of results can be printed#
options(max.print=100000000)


##Conduct genus indicspecies analysis for Estuary Category#

#Take out OTU aka genus table from metadata focused phyloseq object#
W_estuary_seqs_genus = as(otu_table(RFLWSW_Y1_GenusAgenus), "matrix")

#Make it a data frame, transpose it, make it a data frame again#
W_estuary_seqs_genus <- as.data.frame(W_estuary_seqs_genus)
W_estuary_seqs_genus<- t(W_estuary_seqs_genus)
W_estuary_seqs_genus <- as.data.frame(W_estuary_seqs_genus)

#Run indicspecies#
W_estuary_indval_genus = multipatt(W_estuary_seqs_genus, W_estuary_meta_group, control = how(nperm=9999))

#Create a blank csv, then sink all indicspecies results into it#
sink("W_estuary_IS_res_genus.csv")
W_estuary_sig_indval_genus <- summary(W_estuary_indval_genus, indvalcomp=TRUE, alpha=1)
sink()




#Reset max print to default#
options(max.print=99999)

##Look at csv file you sinked into and create a table
##of just ASV.ID, A, B, stat, p.value, and group##

##Continue genus indicspecies analysis for Estuary Category##

#Upload new table#
W_estuary_IS_res_genus <- read.csv(file.choose())
#Make a new p value adjusted column#
W_estuary_IS_res_genus$padj <- p.adjust(W_estuary_IS_res_genus$p.value, method = "BH")

#Keep table of adjusted p values less than 0.05#
W_estuary_IS_sig_res_genus <- filter(W_estuary_IS_res_genus, padj < 0.05)

#Make dataframe with taxanomy added in by making ESV.ID the rownames then adding taxonomy table#
W_estuary_IS_sig_res_taxo_genus <- as.data.frame(W_estuary_IS_sig_res_genus)
rownames(W_estuary_IS_sig_res_taxo_genus) <- W_estuary_IS_sig_res_taxo_genus[,1]
W_estuary_IS_sig_res_taxo_genus = cbind(as(W_estuary_IS_sig_res_taxo_genus, "data.frame"), as(tax_table(RFLWSW_Y1_GenusAgenus)[rownames(W_estuary_IS_sig_res_taxo_genus), ], "matrix"))

#Make dataframe with genuss added in from all sites#
W_estuary_IS_sig_res_taxo_seqs_genus = cbind(as(W_estuary_IS_sig_res_taxo_genus, "data.frame"), as(otu_table(RFLWSW_Y1_GenusAgenus)[rownames(W_estuary_IS_sig_res_taxo_genus), ], "matrix"))

#Remove old rownames#
rownames(W_estuary_IS_sig_res_taxo_seqs_genus) <-NULL

#Make a txt file of results for all sites#
write.csv(as.data.frame(W_estuary_IS_sig_res_taxo_seqs_genus), file="W_estuary_IS_sig_res_taxo_seqs_genus.csv")

#Check the number of indicators in each group after padj#
length(which(sample_data(W_estuary_IS_sig_res_genus)$group == "1")) #97 IRL
length(which(sample_data(W_estuary_IS_sig_res_genus)$group == "2")) #253 SLE

#Order rows by group#
W_estuary_IS_sig_res_taxo_seqs_genus <- arrange(W_estuary_IS_sig_res_taxo_seqs_genus, group)

#Create lists of just significant genuss for each category#
unfiltered_IS_W_estuary_genuss <- subset(W_estuary_IS_sig_res_taxo_seqs_genus, select=c(ESV.ID))

#Create phyloseq objects of the unfiltered indicators using character string version of significant genuss##
charac_unfiltered_IS_W_estuary_genuss <- unfiltered_IS_W_estuary_genuss[,"ESV.ID"]
charac_unfiltered_IS_W_estuary_genuss <- as.character(charac_unfiltered_IS_W_estuary_genuss)
UFISW_estuaryDataAgenus <- prune_taxa(charac_unfiltered_IS_W_estuary_genuss, RFLWSW_Y1_GenusAgenus)

##Create heatmaps to see overall trends, may take a while to load##

unfiltered_IS_W_estuary_genus_heatmap <- plot_heatmap(UFDW_estuaryDataAgenus, taxa.order = charac_unfiltered_IS_W_estuary_genuss, sample.order= "Estuary", sample.label = "Estuary")
unfiltered_IS_W_estuary_genus_heatmap

###MAKE COMBO INDICATOR LISTS###

##Make "Combo" indicator lists from overlapping unfiltered DESeq2 and IS indicators##

#Keep rows in the W_estuary DESEQ object that are also in the IRLSLE IS object#
unfiltered_W_estuary_combo_sig_res_taxo_seqs_genus <- merge(W_estuary_IS_sig_res_genus,W_estuary_DESeq2_sig_res_taxo_seqs_genus, by="ESV.ID")

#Determine the number of indicators in each subcategory of the unfiltered Combo tables#
length(which(sample_data(unfiltered_W_estuary_combo_sig_res_taxo_seqs_genus)$group == "1")) #87 IRL
length(which(sample_data(unfiltered_W_estuary_combo_sig_res_taxo_seqs_genus)$group == "2")) #149 SLE

#Order unfiltered Combo table rows by group#
unfiltered_W_estuary_combo_sig_res_taxo_seqs_genus <- arrange(unfiltered_W_estuary_combo_sig_res_taxo_seqs_genus, group)


#Make a csv file for each of the unfiltered Combo tables#
write.csv(as.data.frame(unfiltered_W_estuary_combo_sig_res_taxo_seqs_genus), file="unfiltered_W_estuary_combo_sig_res_taxo_seqs_genus.csv")

#Create lists of unfiltered significant genuss for each category#
unfiltered_Combo_W_estuary_genuss <- subset(unfiltered_W_estuary_combo_sig_res_taxo_seqs_genus, select=c(ESV.ID))

##Create phyloseq objects for Combo indicators##

#Create phyloseq objects of the unfiltered indicators using character string version of significant genuss##
charac_unfiltered_Combo_W_estuary_genuss <- unfiltered_Combo_W_estuary_genuss[,"ESV.ID"]
charac_unfiltered_Combo_W_estuary_genuss <- as.character(charac_unfiltered_Combo_W_estuary_genuss)
UFCW_estuaryDataAgenus <- prune_taxa(charac_unfiltered_Combo_W_estuary_genuss, RFLWSW_Y1_GenusAgenus)

##Create heatmaps to see overall trends, may take a while to load##

unfiltered_Combo_W_estuary_genus_heatmap <- plot_heatmap(UFDW_estuaryDataAgenus, taxa.order = charac_unfiltered_Combo_W_estuary_genuss, sample.order= "Estuary", sample.label = "Estuary")
unfiltered_Combo_W_estuary_genus_heatmap


###FAMILY###

#Determine the number of samples in each category and subcategory you are testing#
length(which(sample_data(RFLWSW_Y1_FamilyAfamily)$Estuary == "SLE")) #24
length(which(sample_data(RFLWSW_Y1_FamilyAfamily)$Estuary == "IRL")) #42

###RUN DESEQ2 AT family LEVEL ON ESTUARY SUBCATEGORIES###

#Convert phyloseq to deseq2 object centered around a irlvssle factor#
#W_estuary = metadata category testing
#family = taxonomic level
#cudds = DESeq2
#W 



W_estuary_cudds_family = phyloseq_to_deseq2(RFLWSW_Y1_FamilyAfamily,  ~ Estuary)

#Calculate geometric means prior to estimate size factors#
W_estuary_geoMeans_family = apply(counts(W_estuary_cudds_family), 1, gm_mean)
W_estuary_cudds_family = estimateSizeFactors(W_estuary_cudds_family, geoMeans = W_estuary_geoMeans_family)

#Conduct DESEQ2 test#
W_estuary_cudds_family = DESeq(W_estuary_cudds_family, fitType="local")

#Explore the results#
W_estuary_DESeq2_res_family = results(W_estuary_cudds_family)
W_estuary_DESeq2_res_family = W_estuary_DESeq2_res_family[order(W_estuary_DESeq2_res_family$padj, na.last=NA), ]
alpha = 0.05
W_estuary_DESeq2_sig_res_family = W_estuary_DESeq2_res_family[(W_estuary_DESeq2_res_family$padj < alpha), ]

#Make dataframe with taxanomy added in#
W_estuary_DESeq2_sig_res_taxo_family = cbind(as(W_estuary_DESeq2_sig_res_family, "data.frame"), as(tax_table(RFLWSW_Y1_FamilyAfamily)[rownames(W_estuary_DESeq2_sig_res_family), ], "matrix"))

#Make dataframe with familys added in from all sites#
W_estuary_DESeq2_sig_res_taxo_seqs_family = cbind(as(W_estuary_DESeq2_sig_res_taxo_family, "data.frame"), as(otu_table(RFLWSW_Y1_FamilyAfamily)[rownames(W_estuary_DESeq2_sig_res_taxo_family), ], "matrix"))

#Make rownames an actual column and remove old rownames#
W_estuary_DESeq2_sig_res_taxo_seqs_family <- cbind(ESV.ID = rownames(W_estuary_DESeq2_sig_res_taxo_seqs_family), W_estuary_DESeq2_sig_res_taxo_seqs_family)
rownames(W_estuary_DESeq2_sig_res_taxo_seqs_family) <-NULL

#Make a txt file of results for all sites#
write.csv(as.data.frame(W_estuary_DESeq2_sig_res_taxo_seqs_family), file="W_estuary_DESeq2_sig_res_taxo_seqs_family.csv")

#Detemine which subcategory is negative or positive#
W_estuary_DESeq2_res_family
#W_estuary_DESeq2_res_family: SLE (+) vs IRL (-)#

#Determine number of Indicators with padj >0.05 in each subcagetory#
#Number of SLE indicators# 
length(which(W_estuary_DESeq2_sig_res_taxo_seqs_family$log2FoldChange > 0)) #98
#Number of IRL indicators#
length(which(W_estuary_DESeq2_sig_res_taxo_seqs_family$log2FoldChange < 0)) #90

#Order rows by log2FoldChange#
W_estuary_DESeq2_sig_res_taxo_seqs_family <- arrange(W_estuary_DESeq2_sig_res_taxo_seqs_family, log2FoldChange)

#Create lists of just significant familys for each category#
unfiltered_DESeq2_W_estuary_familys <- subset(W_estuary_DESeq2_sig_res_taxo_seqs_family, select=c(ESV.ID))

##Create phyloseq objects for DESeq2 indicators##

#Create phyloseq objects of the unfiltered indicators using character string version of significant familys##
charac_unfiltered_DESeq2_W_estuary_familys <- unfiltered_DESeq2_W_estuary_familys[,"ESV.ID"]
charac_unfiltered_DESeq2_W_estuary_familys <- as.character(charac_unfiltered_DESeq2_W_estuary_familys)
UFDW_estuaryDataAfamily <- prune_taxa(charac_unfiltered_DESeq2_W_estuary_familys, RFLWSW_Y1_FamilyAfamily)

unfiltered_DESeq2_W_estuary_family_heatmap <- plot_heatmap(UFDW_estuaryDataAfamily, taxa.order = charac_unfiltered_DESeq2_W_estuary_familys, sample.order= "Estuary", sample.label = "Estuary")
unfiltered_DESeq2_W_estuary_family_heatmap

###CONDUCT INDICSPECIES ANALYSIS AT family LEVEL###

#Increase max print so that all of results can be printed#
options(max.print=100000000)


##Conduct family indicspecies analysis for Estuary Category#

#Take out OTU aka family table from metadata focused phyloseq object#
W_estuary_seqs_family = as(otu_table(RFLWSW_Y1_FamilyAfamily), "matrix")

#Make it a data frame, transpose it, make it a data frame again#
W_estuary_seqs_family <- as.data.frame(W_estuary_seqs_family)
W_estuary_seqs_family<- t(W_estuary_seqs_family)
W_estuary_seqs_family <- as.data.frame(W_estuary_seqs_family)

#Run indicspecies#
W_estuary_indval_family = multipatt(W_estuary_seqs_family, W_estuary_meta_group, control = how(nperm=9999))

#Create a blank csv, then sink all indicspecies results into it#
sink("W_estuary_IS_res_family.csv")
W_estuary_sig_indval_family <- summary(W_estuary_indval_family, indvalcomp=TRUE, alpha=1)
sink()




#Reset max print to default#
options(max.print=99999)

##Look at csv file you sinked into and create a table
##of just ASV.ID, A, B, stat, p.value, and group##

##Continue family indicspecies analysis for Estuary Category##

#Upload new table#
W_estuary_IS_res_family <- read.csv(file.choose())
#Make a new p value adjusted column#
W_estuary_IS_res_family$padj <- p.adjust(W_estuary_IS_res_family$p.value, method = "BH")

#Keep table of adjusted p values less than 0.05#
W_estuary_IS_sig_res_family <- filter(W_estuary_IS_res_family, padj < 0.05)

#Make dataframe with taxanomy added in by making ESV.ID the rownames then adding taxonomy table#
W_estuary_IS_sig_res_taxo_family <- as.data.frame(W_estuary_IS_sig_res_family)
rownames(W_estuary_IS_sig_res_taxo_family) <- W_estuary_IS_sig_res_taxo_family[,1]
W_estuary_IS_sig_res_taxo_family = cbind(as(W_estuary_IS_sig_res_taxo_family, "data.frame"), as(tax_table(RFLWSW_Y1_FamilyAfamily)[rownames(W_estuary_IS_sig_res_taxo_family), ], "matrix"))

#Make dataframe with familys added in from all sites#
W_estuary_IS_sig_res_taxo_seqs_family = cbind(as(W_estuary_IS_sig_res_taxo_family, "data.frame"), as(otu_table(RFLWSW_Y1_FamilyAfamily)[rownames(W_estuary_IS_sig_res_taxo_family), ], "matrix"))

#Remove old rownames#
rownames(W_estuary_IS_sig_res_taxo_seqs_family) <-NULL

#Make a txt file of results for all sites#
write.csv(as.data.frame(W_estuary_IS_sig_res_taxo_seqs_family), file="W_estuary_IS_sig_res_taxo_seqs_family.csv")

#Check the number of indicators in each group after padj#
length(which(sample_data(W_estuary_IS_sig_res_family)$group == "1")) #57 IRL
length(which(sample_data(W_estuary_IS_sig_res_family)$group == "2")) #137 SLE

#Order rows by group#
W_estuary_IS_sig_res_taxo_seqs_family <- arrange(W_estuary_IS_sig_res_taxo_seqs_family, group)

#Create lists of just significant familys for each category#
unfiltered_IS_W_estuary_familys <- subset(W_estuary_IS_sig_res_taxo_seqs_family, select=c(ESV.ID))

#Create phyloseq objects of the unfiltered indicators using character string version of significant familys##
charac_unfiltered_IS_W_estuary_familys <- unfiltered_IS_W_estuary_familys[,"ESV.ID"]
charac_unfiltered_IS_W_estuary_familys <- as.character(charac_unfiltered_IS_W_estuary_familys)
UFISW_estuaryDataAfamily <- prune_taxa(charac_unfiltered_IS_W_estuary_familys, RFLWSW_Y1_FamilyAfamily)

##Create heatmaps to see overall trends, may take a while to load##

unfiltered_IS_W_estuary_family_heatmap <- plot_heatmap(UFDW_estuaryDataAfamily, taxa.order = charac_unfiltered_IS_W_estuary_familys, sample.order= "Estuary", sample.label = "Estuary")
unfiltered_IS_W_estuary_family_heatmap

###MAKE COMBO INDICATOR LISTS###

##Make "Combo" indicator lists from overlapping unfiltered DESeq2 and IS indicators##

#Keep rows in the W_estuary DESEQ object that are also in the IRLSLE IS object#
unfiltered_W_estuary_combo_sig_res_taxo_seqs_family <- merge(W_estuary_IS_sig_res_family,W_estuary_DESeq2_sig_res_taxo_seqs_family, by="ESV.ID")

#Determine the number of indicators in each subcategory of the unfiltered Combo tables#
length(which(sample_data(unfiltered_W_estuary_combo_sig_res_taxo_seqs_family)$group == "1")) #52 IRL
length(which(sample_data(unfiltered_W_estuary_combo_sig_res_taxo_seqs_family)$group == "2")) #93 SLE

#Order unfiltered Combo table rows by group#
unfiltered_W_estuary_combo_sig_res_taxo_seqs_family <- arrange(unfiltered_W_estuary_combo_sig_res_taxo_seqs_family, group)


#Make a csv file for each of the unfiltered Combo tables#
write.csv(as.data.frame(unfiltered_W_estuary_combo_sig_res_taxo_seqs_family), file="unfiltered_W_estuary_combo_sig_res_taxo_seqs_family.csv")

#Create lists of unfiltered significant familys for each category#
unfiltered_Combo_W_estuary_familys <- subset(unfiltered_W_estuary_combo_sig_res_taxo_seqs_family, select=c(ESV.ID))

##Create phyloseq objects for Combo indicators##

#Create phyloseq objects of the unfiltered indicators using character string version of significant familys##
charac_unfiltered_Combo_W_estuary_familys <- unfiltered_Combo_W_estuary_familys[,"ESV.ID"]
charac_unfiltered_Combo_W_estuary_familys <- as.character(charac_unfiltered_Combo_W_estuary_familys)
UFCW_estuaryDataAfamily <- prune_taxa(charac_unfiltered_Combo_W_estuary_familys, RFLWSW_Y1_FamilyAfamily)

##Create heatmaps to see overall trends, may take a while to load##

unfiltered_Combo_W_estuary_family_heatmap <- plot_heatmap(UFDW_estuaryDataAfamily, taxa.order = charac_unfiltered_Combo_W_estuary_familys, sample.order= "Estuary", sample.label = "Estuary")
unfiltered_Combo_W_estuary_family_heatmap


###INDICATOR EFFECTIVENESS METHOD###

##NOTES ON METHOD##
#"affected' = subcategory that is more affected by the stressor you are      concerned about. For example SLE for the Estuary category becasue it has     more freshwater discharges, 3MC for the Muck category because it has more    muck characteristics, HiHi for the TOM/Cu category becasue it has more       copper#

#"non-affected" = IRL, 0MC, HiLo#

#microbially-predicted "affected" sample = samples in the partitioning around medoids (PAM) cluster with the most metadata-defined "affected" samples#

#microbially-predicted "non-affected" sample = samples in the PAM cluster with the most metadata-defined "non-affected" samples#

#Overall Idea: Use PAM clustering to split samples between the "affected" and "non-affected" clusters based upon the indicators (microbially-predicted) and see how well this classification matches the metadata-defined classification by using the product of specificity and sensitivity#

#Copy results to an Excel Sheet for further statistical testing, each combo of factors should have four percentage types (Product is Sensitivity x Specificity), for example:#
#Indicator_Test	Taxonomic_Level	Filtering_Status	Metadata_Tested	Percent_Type	Percentage#
#Indicspecies   ASV             Unfiltered        W_Estuary         Total         98.78#
#Indicspecies   ASV             Unfiltered        W_Estuary         Sensitivity   96.23#
#Indicspecies   ASV             Unfiltered        W_Estuary         Specificity   99.15#
#Indicspecies   ASV             Unfiltered        W_Estuary         Product       95.41#


###ORIGINAL DATASETS
###Split samples of orginal datasets into two clusters to see how it will compare to splitting based upon indicators###

##Test on Muck and W_Estuary categories at asv level##

#Remove any samples with no asvs from target phyloseq#
nsamples(RFLWSW_Y1_DataAasv) #66
RFLWSW_Y1_DataAasv = prune_samples(sample_sums(RFLWSW_Y1_DataAasv)>=1, RFLWSW_Y1_DataAasv)
nsamples(RFLWSW_Y1_DataAasv) #66

#Make a matrix with the sample_data and relabel subcategories to make metadata-defined column#
OriginalW_Est_asv <- as.matrix(sample_data(RFLWSW_Y1_DataAasv))
OriginalW_Est_asv[OriginalW_Est_asv ==  "IRL"] <- 1
OriginalW_Est_asv[OriginalW_Est_asv ==  "SLE"] <- 2
OriginalW_Est_asv<-as.data.frame(OriginalW_Est_asv)


#Create Bray-Curtis distance matrix from asv table#
bc <- vegdist(t(otu_table(RFLWSW_Y1_DataAasv)))

#Conduct pam clustering based upon 2 clusters#
pam.res <- pam(bc,2, diss=TRUE)

#Add a microbially-predicted column (Orginal_asv) to matrix above for comparison to metadata-defined column#
OriginalW_Est_asv <- cbind(OriginalW_Est_asv, Orginal_asv = pam.res$cluster)

##Test Estuary category##

#Test  to see how well the metadata-defined and microbially-predicted columns match one another aka Total Efficiency; if the number is below 0.5, then switch the numbers in the metadata-defiend column#

#Test for Estuary Effeciency# 0.8636364
sum(OriginalW_Est_asv$Estuary == OriginalW_Est_asv$Orginal_asv)/nrow(OriginalW_Est_asv)

##Test indicator sensitivity aka the true postives (TP)/ (TP + false negatives (FN))##
#TP = metadata-defined "affected" sample correctly placed in the microbiall   y-predicted "affected" cluster#
#FN = metadata-defined "affected" sample incorrectly placed in the           microbially-predicted "non-affected" cluster#
#TP+FN = all metadata-defined "affected" samples

#Test for SLE Effeciency aka Sensitvity# 0.625
Orginal_SLE_asv <- filter(OriginalW_Est_asv, Estuary == "2")
sum(Orginal_SLE_asv$Estuary == Orginal_SLE_asv$Orginal_asv)/nrow(Orginal_SLE_asv)

##Test indicator specificity aka the true negatives (TN)/ (TN + false positives (FP))##
#TP = metadata-defined "affected" sample correctly placed in the microbiall   y-predicted "affected" cluster#
#FN = metadata-defined "affected" sample incorrectly placed in the           microbially-predicted "non-affected" cluster#
#TP+FN = all metadata-defined "affected" samples

#Test for IRL Effeciency aka Specificity# 1
Orginal_IRL_asv <- filter(OriginalW_Est_asv, Estuary == "1")
sum(Orginal_IRL_asv$Estuary == Orginal_IRL_asv$Orginal_asv)/nrow(Orginal_IRL_asv)

#PCoA with Color by W_Estuary#
TotalPCoA <- ordinate(RFLWSW_Y1_DataAasv,"PCoA")
p = plot_ordination(RFLWSW_Y1_DataAasv, TotalPCoA, color="Estuary") + 
  ggtitle("Total asvs separated by W_Estuary Category PCoA with Bray-Curtis distance") +
  geom_text(aes(label=OriginalW_Est_asv$Orginal_asv),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p

###TEST INDICSPECIES ESTUARY INDICATORS###

##Test unfiltered Estuary asvs##

#Remove any samples with no asvs from target phyloseq#
nsamples(UFISW_estuaryDataAasv) #66
UFISW_estuaryDataAasv = prune_samples(sample_sums(UFISW_estuaryDataAasv)>=1, UFISW_estuaryDataAasv)
nsamples(UFISW_estuaryDataAasv) #66

#Make a matrix with the sample_data and relabel subcategories to make metadata-defined column#
UF_ISW_Est_asv <- as.matrix(sample_data(UFISW_estuaryDataAasv))
UF_ISW_Est_asv[UF_ISW_Est_asv ==  "SLE"] <- 2
UF_ISW_Est_asv[UF_ISW_Est_asv ==  "IRL"] <- 1
UF_ISW_Est_asv<-as.data.frame(UF_ISW_Est_asv)

#Create Bray-Curtis distance matrix from asv table#
bc <- vegdist(t(otu_table(UFISW_estuaryDataAasv)))

#Conduct PAM clustering based upon 2 clusters#
pam.res <- pam(bc,2, diss=TRUE)

#Add a microbially-predicted column to matrix above for comparison to metadata-defined column#
UF_ISW_Est_asv <- cbind(UF_ISW_Est_asv, MP_UFISW_Estasv = pam.res$cluster)

#Test Total W_Estuary Efficiency# 0.8636364
sum(UF_ISW_Est_asv$Estuary == UF_ISW_Est_asv$MP_UFISW_Estasv)/nrow(UF_ISW_Est_asv)

#Test SLE Efficiency aka Sensitvity# 0.625
UF_IS_SLE_asv <- filter(UF_ISW_Est_asv, Estuary == "2")
sum(UF_IS_SLE_asv$Estuary == UF_IS_SLE_asv$MP_UFISW_Estasv)/nrow(UF_IS_SLE_asv)

#Test IRL Efficiency aka Specificity# 1
UF_IS_IRL_asv <- filter(UF_ISW_Est_asv, Estuary == "1")
sum(UF_IS_IRL_asv$Estuary == UF_IS_IRL_asv$MP_UFISW_Estasv)/nrow(UF_IS_IRL_asv)

TotalPCoA <- ordinate(UFISW_estuaryDataAasv,"PCoA")
p = plot_ordination(UFISW_estuaryDataAasv, TotalPCoA, color="Estuary") + 
  ggtitle("Unfiltered IS W_Estuary asv Indicators PCoA with Bray-Curtis distance") +
  geom_text(aes(label=UF_ISW_Est_asv$MP_UFISW_Estasv),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p

###TEST DESEQ2 ESTUARY INDICATORS###

##Test unfiltered Estuary asvs##

#Remove any samples with no asvs from target phyloseq#
nsamples(UFDW_estuaryDataAasv) #66
UFDW_estuaryDataAasv = prune_samples(sample_sums(UFDW_estuaryDataAasv)>=1, UFDW_estuaryDataAasv)
nsamples(UFDW_estuaryDataAasv) #66

#Make a matrix with the sample_data and relabel subcategories to make metadata-defined column#
UF_DW_Est_asv <- as.matrix(sample_data(UFDW_estuaryDataAasv))
UF_DW_Est_asv[UF_DW_Est_asv ==  "SLE"] <- 2
UF_DW_Est_asv[UF_DW_Est_asv ==  "IRL"] <- 1
UF_DW_Est_asv<-as.data.frame(UF_DW_Est_asv)


#Create Bray-Curtis distance matrix from asv table#
bc <- vegdist(t(otu_table(UFDW_estuaryDataAasv)))

#Conduct PAM clustering based upon 2 clusters#
pam.res <- pam(bc,2, diss=TRUE)

#Add a microbially-predicted column to matrix above for comparison to metadata-defined column#
UF_DW_Est_asv <- cbind(UF_DW_Est_asv, MP_UFDW_Estasv = pam.res$cluster)

#Test Total W_Estuary Efficiency# 0.8636364
sum(UF_DW_Est_asv$Estuary == UF_DW_Est_asv$MP_UFDW_Estasv)/nrow(UF_DW_Est_asv)

#Test SLE Efficiency aka Sensitivity# 0.625
UF_D_SLE_asv <- filter(UF_DW_Est_asv, Estuary == "2")
sum(UF_D_SLE_asv$Estuary == UF_D_SLE_asv$MP_UFDW_Estasv)/nrow(UF_D_SLE_asv)

#Test IRL Efficiency aka Specificity# 1
UF_D_IRL_asv <- filter(UF_DW_Est_asv, Estuary == "1")
sum(UF_D_IRL_asv$Estuary == UF_D_IRL_asv$MP_UFDW_Estasv)/nrow(UF_D_IRL_asv)


TotalPCoA <- ordinate(UFDW_estuaryDataAasv,"PCoA")
p = plot_ordination(UFDW_estuaryDataAasv, TotalPCoA, color="Estuary") + 
  ggtitle("Unfiltered DESeq2 W_Estuary asv Indicators PCoA with Bray-Curtis distance") +
  geom_text(aes(label=UF_DW_Est_asv$MP_UFDW_Estasv),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p

###TEST COMBO ESTUARY INDICATORS###

##Test unfiltered Estuary asvs##

#Remove any samples with no asvs from target phyloseq#
nsamples(UFCW_estuaryDataAasv) #66
UFCW_estuaryDataAasv = prune_samples(sample_sums(UFCW_estuaryDataAasv)>=1, UFCW_estuaryDataAasv)
nsamples(UFCW_estuaryDataAasv) #66

#Make a matrix with the sample_data and relabel subcategories to make metadata-defined column#
UF_CW_Est_asv <- as.matrix(sample_data(UFCW_estuaryDataAasv))
UF_CW_Est_asv[UF_CW_Est_asv ==  "SLE"] <- 2
UF_CW_Est_asv[UF_CW_Est_asv ==  "IRL"] <- 1
UF_CW_Est_asv<-as.data.frame(UF_CW_Est_asv)


#Create Bray-Curtis distance matrix from asv table#
bc <- vegdist(t(otu_table(UFCW_estuaryDataAasv)))

#Conduct PAM clustering based upon 2 clusters#
pam.res <- pam(bc,2, diss=TRUE)

#Add a microbially-predicted column to matrix above for comparison to metadata-defined column#
UF_CW_Est_asv <- cbind(UF_CW_Est_asv, MP_UFCW_Estasv = pam.res$cluster)

#Test Total W_Estuary Efficiency# 0.8636364
sum(UF_CW_Est_asv$Estuary == UF_CW_Est_asv$MP_UFCW_Estasv)/nrow(UF_CW_Est_asv)

#Test SLE Efficiency aka Sensitivity# 0.625
UF_C_SLE_asv <- filter(UF_CW_Est_asv, Estuary == "2")
sum(UF_C_SLE_asv$Estuary == UF_C_SLE_asv$MP_UFCW_Estasv)/nrow(UF_C_SLE_asv)

#Test IRL Efficiency aka Specificity# 1
UF_C_IRL_asv <- filter(UF_CW_Est_asv, Estuary == "1")
sum(UF_C_IRL_asv$Estuary == UF_C_IRL_asv$MP_UFCW_Estasv)/nrow(UF_C_IRL_asv)

TotalPCoA <- ordinate(UFCW_estuaryDataAasv,"PCoA")
p = plot_ordination(UFCW_estuaryDataAasv, TotalPCoA, color="Estuary") + 
  ggtitle("Unfiltered Final W_Estuary asv Indicators PCoA with Bray-Curtis distance") +
  geom_text(aes(label=UF_CW_Est_asv$MP_UFCW_Estasv),hjust=0, vjust=0)
p$layers <- p$layers[-1]
p







#Extract metadata table from the genus phyloseq object#
LWSW_Meta_Y1 <- as(sample_data(RFLWSW_Y1_GenusAgenus), "data.frame")


##Create new columns in data table related to the Estuary AIP##

#Determine the sequence sums from Unfiltered, Combo, Genus Level IRL indicators and upload as a new column#
LWSW_Meta_Y1$UF_Combo_Genus_IRL <- prune_taxa(as.character(unfiltered_W_estuary_combo_sig_res_taxo_seqs_genus$ESV.ID[unfiltered_W_estuary_combo_sig_res_taxo_seqs_genus$log2FoldChange<0]), RFLWSW_Y1_GenusAgenus) %>% sample_sums() 

#Determine the sequence sums from Unfiltered, Combo, Genus Level SLE indicators and upload as a new column#
LWSW_Meta_Y1$UF_Combo_Genus_SLE <- prune_taxa(as.character(unfiltered_W_estuary_combo_sig_res_taxo_seqs_genus$ESV.ID[unfiltered_W_estuary_combo_sig_res_taxo_seqs_genus$log2FoldChange>0]), RFLWSW_Y1_GenusAgenus) %>% sample_sums() 

#Calculate Estuary AIP and add in new column#
LWSW_Meta_Y1$UF_Combo_Genus_Estuary <- with(LWSW_Meta_Y1, 100*UF_Combo_Genus_SLE/(UF_Combo_Genus_SLE+UF_Combo_Genus_IRL))

#Extract metadata table from the genus phyloseq object#
LWSW_Meta_Y2 <- as(sample_data(RFLWSW_Y2_GenusAgenus), "data.frame")


##Create new columns in data table related to the Estuary AIP##

#Determine the sequence sums from Unfiltered, Combo, Genus Level IRL indicators and upload as a new column#
LWSW_Meta_Y2$UF_Combo_Genus_IRL <- prune_taxa(as.character(unfiltered_W_estuary_combo_sig_res_taxo_seqs_genus$ESV.ID[unfiltered_W_estuary_combo_sig_res_taxo_seqs_genus$log2FoldChange<0]), RFLWSW_Y2_GenusAgenus) %>% sample_sums() 

#Determine the sequence sums from Unfiltered, Combo, Genus Level SLE indicators and upload as a new column#
LWSW_Meta_Y2$UF_Combo_Genus_SLE <- prune_taxa(as.character(unfiltered_W_estuary_combo_sig_res_taxo_seqs_genus$ESV.ID[unfiltered_W_estuary_combo_sig_res_taxo_seqs_genus$log2FoldChange>0]), RFLWSW_Y2_GenusAgenus) %>% sample_sums() 

#Calculate Estuary AIP and add in new column#
LWSW_Meta_Y2$UF_Combo_Genus_Estuary <- with(LWSW_Meta_Y2, 100*UF_Combo_Genus_SLE/(UF_Combo_Genus_SLE+UF_Combo_Genus_IRL))

LWSW_Meta_Total <- rbind(LWSW_Meta_Y1, LWSW_Meta_Y2)


#Add new columns combiing sampling year and target characteristics

LWSW_Meta_Total$Sampling.Year <- paste(LWSW_Meta_Total$AbbSeasonAbbYear)

LWSW_Meta_Total$Sampling.Year[LWSW_Meta_Total$Sampling.Year ==  "W16"] <- 1
LWSW_Meta_Total$Sampling.Year[LWSW_Meta_Total$Sampling.Year ==  "D17"] <- 1
LWSW_Meta_Total$Sampling.Year[LWSW_Meta_Total$Sampling.Year ==  "W17"] <- 2
LWSW_Meta_Total$Sampling.Year[LWSW_Meta_Total$Sampling.Year ==  "D18"] <- 2



LWSW_Meta_Total$Estuary.Sampling.Year <- paste(LWSW_Meta_Total$Estuary,LWSW_Meta_Total$Sampling.Year)

LWSW_Meta_Total$Season.Sampling.Year <- paste(LWSW_Meta_Total$Season,LWSW_Meta_Total$Sampling.Year)



###VISUALIZE AND TEST METRICS###

##Statistically test AIP differences between subcategories##

#Estuary AIP Stats#
kruskal.test(UF_Combo_Genus_Estuary ~ Estuary.Sampling.Year, data = LWSW_Meta_Total)

dunnTest(UF_Combo_Genus_Estuary ~ Estuary.Sampling.Year, data = LWSW_Meta_Total, method="bh")

Dunn <- dunnTest(UF_Combo_Genus_Estuary ~ Estuary.Sampling.Year, data = LWSW_Meta_Total, method="bh")

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

##Create Box Plots to visualize differences between subcategories##

#Estuary AIP BP#
ggplot(data = LWSW_Meta_Total, aes( x = Estuary.Sampling.Year, y = UF_Combo_Genus_Estuary, color=Season)) +
  geom_boxplot() +
  ylab("Affected Indicator Percentage") +
  xlab("Estuary by Sampling Year") +
  theme(axis.text.x  = element_text(angle=45, hjust=1))+
  ggtitle("Affected Indicator Percentage") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  annotate("text", x = 1:4, y = c(40, 105, 87, 105), label = c("a", "c", "b", "c")) +
  scale_x_discrete(limits=c("IRL 1", "SLE 1", "IRL 2", "SLE 2"), labels=c("IRL Year 1", "SLE Year 1", "IRL Year 2", "SLE Year 2")) 
  #theme(legend.position="bottom")
  #theme(legend.position = "none")

#Estuary AIP by Site#
ggplot(data = LWSW_Meta_Total, aes( x = Site, y = UF_Combo_Genus_Estuary)) +
  geom_boxplot() +
  theme(axis.text.x  = element_text(angle=45, hjust=1))+
  ylab("Affected Indicator Percentage") +
  xlab("Site") +
  #ggtitle("Overall Sediment Samples Estuary AIP")+
  #theme(legend.position = "none") +
  #scale_x_discrete(labels=c("Barber Bridge", "Fort Pierce", "Harbor Branch Channel", "Harbortown Marina", "Hobe Sound", "Jensen Beach", "Jupiter Narrows", "Linkport", "Manatee Pocket", "Melbourne Causeway", "Merritt Island Causeway", "Middle Estuary", "North Fork", "Round Island", "Sebastian Inlet", "South Fork", "South Fork 2", "Vero Beach", "Vero Beach Marina"))+
  theme(plot.title = element_text(hjust = 0.5))

#Estuary AIP by Site by Season#
ggplot(data = LWSW_Meta_Total, aes( x = Site, y = UF_Combo_Genus_Estuary, color=Season.Sampling.Year)) +
  geom_boxplot() +
  theme(axis.text.x  = element_text(angle=45, hjust=1))+
  ylab("Affected Indicator Percentage") +
  xlab("Site") 
#scale_x_discrete(labels=c("Barber Bridge", "Fort Pierce", "Harbor Branch Channel", "Harbortown Marina", "Hobe Sound", "Jensen Beach", "Jupiter Narrows", "Linkport", "Manatee Pocket", "Melbourne Causeway", "Merritt Island Causeway", "Middle Estuary", "North Fork", "Round Island", "Sebastian Inlet", "South Fork", "South Fork 2", "Vero Beach", "Vero Beach Marina"))+
scale_color_manual(values = c("darkgray", "black"))

