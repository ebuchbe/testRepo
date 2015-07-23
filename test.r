#############Script for searching the up-/downregulated genes in the set of genes that show the GATA enriched factor

results_findPnrMotiv <- read.table("testfindPnrMotifFinal.txt", header= TRUE, sep = "\t", fill = TRUE, quote = "")
PnrMotifGenes <- unique(results_findPnrMotiv[,"Name"]) #this are the 1675 genes that have the GATA motif I found
write.table (PnrMotifGenes, file = "PnrMotifGenes.txt", row.names = FALSE, col.names = FALSE, quote = FALSE) ##now convert the gene symbols into FlyBaseIds
conversion_table <- read.table("PnrMotifGenesConversionTable.txt", header= TRUE, sep = "\t", fill = TRUE, quote = "")
length(unique(conversion_table[,"submitted_id"])) #just to check
newConversionTable <- conversion_table[which(conversion_table[,"current_symbol"] %in% PnrMotifGenes),]
flybaseSymbols <- newConversionTable[,"converted_id"] ####this are the 1675 genes that have the GATA motif "high confidence list"
write.table(flybaseSymbols, file = "PnrMotifGenesFlyBase.txt", row.names= FALSE, col.names = FALSE, quote = FALSE)

DESeq2 <- read.table("DESeq2_OreR_96hvs120h_results.txt", sep = "\t")
DESeq2sign<-subset(DESeq2, padj<=0.05)
upregulated <- rownames(subset(DESeq2sign,log2FoldChange >= 0))
downregulated <- rownames(subset(DESeq2sign,log2FoldChange <= 0))
length(upregulated)
length(downregulated)
upregTargets <- subset(flybaseSymbols, flybaseSymbols %in% upregulated)
downregTargets <- subset(flybaseSymbols, flybaseSymbols %in% downregulated)
length(upregTargets)
length(downregTargets)

##### do a GO analysis with the genes that are up and downregulated
write.table(upregTargets, file = "upregTargetsORE72vs96h.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(downregTargets, file = "downregTargetsORE72vs96h.txt", row.names = FALSE, col.names = FALSE, quote = FALSE) ###convert the IDs 
up <- read.table("upregTargetsORE72vs96hconversion.txt", header= FALSE, sep = "\t", fill = TRUE, quote = "")
upSymbol <- up[,"V4"]
down <- read.table("downregTargetsORE72vs96hconversion.txt", header= FALSE, sep = "\t", fill = TRUE, quote = "")
downSymbol <- down[,"V4"]
write.table(upSymbol, file = "upSymbolOre72vs96h.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(downSymbol, file = "downSymbolOre72vs96h.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
##########I did the GO analysis with HOMER, the output directory includes a .txt file that has all the biological processes 

bp <- read.table("biological_process.txt", header= TRUE, sep = "\t", fill = TRUE, quote = "")
signBP<-subset(bp, Enrichment<=0.05)
