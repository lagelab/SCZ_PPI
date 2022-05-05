##########################################################################################
## Overlap enrichment analysis between SCZ PPI networks vs. cell-type-specific DEGs
## found in SCZ patients vs. controls (using scRNA-seq data from Ruzicka et al.)
##
## Author: Yu-Han Hsu
##########################################################################################

rm(list=ls())

library(readxl)
library(genoppi)
library(ggplot2)


# ----------------------------------------------------------------------------------------
# PPI networks to test

networks <- c('all_combined','wk2_combined','wk3_combined','wk4_combined','wk7_combined',
	'CACNA1C_combined','CUL3_wk7','HCN1_combined',
	'RIMS1_combined','SYNGAP1_combined','TCF4_wk1')
names <- c('All combined','Week 2','Week 3','Week 4','Week 7',
	'CACNA1C','CUL3','HCN1','RIMS1','SYNGAP1','TCF4')

# read in interactors and non-interactors for each network
intDf <- read.table('../data/SCZ_MasterInteractorTable_stringent.txt',
	header=T,sep='\t',stringsAsFactors=F)
intDf <- subset(intDf, ListName %in% networks)[,c('ListName','GeneName','IsInteractor')]
colnames(intDf) <- c('ListName','gene','significant')
dim(intDf)

# network size (int count)
sizeDf <- subset(intDf, significant)
sizeDf <- data.frame(table(sizeDf$ListName))
colnames(sizeDf) <- c('ListName','IntCount')
for (i in 1:length(networks)){
	sizeDf$Network[as.character(sizeDf$ListName)==networks[i]] <- names[i]
}


# ----------------------------------------------------------------------------------------
# read in SCZ DEGs in each cell type and perform overlap enrichment test for each network

sheets <- excel_sheets('../data/Ruzicka2020_media-5.xlsx')
overlapDf <- NULL
upOverlapDf <- NULL
dnOverlapDf <- NULL

# iterate through cell types
for (s in sheets) {
	print(s)
	
	# up-regulated DEGs
	tempDf <- read_excel('../data/Ruzicka2020_media-5.xlsx',sheet=s)
	upDf <- data.frame(gene=tempDf$Gene, significant=TRUE)

	# down-regulated DEGs
	tempDf <- read_excel('../data/Ruzicka2020_media-6.xlsx',sheet=s)
	dnDf <- data.frame(gene=tempDf$Gene, significant=TRUE)	
	
	# all DEGs
	degDf <- rbind(upDf,dnDf)
	
	# iterate through networks, perform hypergeometric tests
	for (i in 1:length(networks)){
		netDf <- subset(intDf,ListName==networks[i])[,c('gene','significant')]

		# up-regulated DEGs
		results <- calc_hyper(netDf,upDf,data.frame(intersectN=F))
		resDf <- data.frame(CellType=paste(s,':',nrow(upDf),sep=''),
			Network=names[i],results$statistics[,c('successInSample_count',
			'sample_count','success_count','population_count','pvalue')])
		colnames(resDf) <- c('CellType','Network','OverlapCount',
			'DegCount','IntCount','PopCount','Pvalue')
		resDf$OverlapGenes <- paste(results$genes[[1]]$successInSample_genes,collapse=',')
		upOverlapDf <- rbind(upOverlapDf,resDf)
		
		# down-regulated DEGs
		results <- calc_hyper(netDf,dnDf,data.frame(intersectN=F))
		resDf <- data.frame(CellType=paste(s,':',nrow(dnDf),sep=''),
			Network=names[i],results$statistics[,c('successInSample_count',
			'sample_count','success_count','population_count','pvalue')])
		colnames(resDf) <- c('CellType','Network','OverlapCount',
			'DegCount','IntCount','PopCount','Pvalue')
		resDf$OverlapGenes <- paste(results$genes[[1]]$successInSample_genes,collapse=',')
		dnOverlapDf <- rbind(dnOverlapDf,resDf)

		# all DEGs
		results <- calc_hyper(netDf,degDf,data.frame(intersectN=F))
		resDf <- data.frame(CellType=paste(s,':',nrow(degDf),sep=''),
			Network=names[i],results$statistics[,c('successInSample_count',
			'sample_count','success_count','population_count','pvalue')])
		colnames(resDf) <- c('CellType','Network',
			'OverlapCount','DegCount','IntCount','PopCount','Pvalue')
		resDf$OverlapGenes <- paste(results$genes[[1]]$successInSample_genes,collapse=',')
		overlapDf <- rbind(overlapDf,resDf)
	}
}


# output enirchment results tables
outDf <- data.frame(DegType=c(rep('All DEGs',nrow(overlapDf)),
	rep('Up-regulated DEGs',nrow(upOverlapDf)),
	rep('Down-regulated DEGs',nrow(dnOverlapDf))),
	rbind(overlapDf,upOverlapDf,dnOverlapDf))
write.table(outDf,'../output/SCZ_scDegEnrichment.txt',sep='\t',row.names=F,quote=F)


# ----------------------------------------------------------------------------------------
# plot enrichment results in heat maps

# set p-value label based on nominal or Bonferroni significance
bonfN <- length(sheets)*length(networks) # 20 * 11
outDf$Sig <- ifelse(outDf$Pvalue < 0.05,"*","")
outDf$Sig[outDf$Pvalue<0.05/bonfN] <- "**"

# show overlap count for nominally significant results
outDf$SigText <- ifelse(outDf$Pvalue < 0.05,paste(outDf$OverlapCount,outDf$Sig,sep=''),"")
	
# set network label and plotting order
outDf$IntCount <- sizeDf$IntCount[match(outDf$Network,sizeDf$Network)]
outDf$Network <- paste(outDf$Network,' (',outDf$IntCount,')',sep='')
netOrder <- paste(rev(names),
	' (',sizeDf$IntCount[match(rev(names),sizeDf$Network)],')',sep='')
outDf$Network <- factor(outDf$Network, levels=netOrder)

# set cell type label and plotting order
outDf$CellType <- paste(gsub(':','\n(',outDf$CellType,),')',sep='')
outDf$CellType <- factor(outDf$CellType,levels=unique(outDf$CellType))


# heat maps with colors based on -log10(P)
pdf("../output/SCZ_scDegEnrichment_HeatMaps.pdf",height=4,width=8)

for (dType in unique(outDf$DegType)) {

	p <- ggplot(subset(outDf,DegType==dType),
			aes(x=CellType,y=Network,fill=-log10(Pvalue))) +
		geom_tile() + geom_text(aes(label=SigText),size=3) +
		scale_fill_gradient(name=expression(paste(-log[10],"(P-value)",sep="")),
			low="white", high="red") +
		xlab("Cell type") + ylab("Network") + ggtitle(dType) +
		theme_bw() + 
		theme(legend.position="bottom",legend.key.size=unit(0.9,"line"),
			legend.title=element_text(size=9),legend.text=element_text(size=8),
			legend.margin=margin(0,0,0,0),legend.box.margin=margin(-5,-5,-5,-5),
			axis.title=element_text(size=9),axis.text.y=element_text(size=8),
			axis.text.x=element_text(size=8,angle=90,vjust=0.5,hjust=1))

	print(p)

}

dev.off()

