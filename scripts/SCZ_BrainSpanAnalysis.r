##########################################################################################
## BrainSpan expression profiles for various SCZ gene sets
##
## Author: Yu-Han Hsu
##########################################################################################

rm(list=ls())

library(dplyr)
library(reshape2)
library(ggplot2)


# ----------------------------------------------------------------------------------------
# BrainSpan data (exon microarray summarized to genes)

# read in expression data (rows = genes, columns = samples)
expDf <- read.csv("../data/BrainSpan/array_expression_matrix.processed.csv", header=T, stringsAsFactor=F)

# subset to data from samples in the four frontal cortex regions
expDf <- expDf %>% select(contains("gene"),
	contains("MFC"),contains("DFC"),contains("OFC"),contains("VFC"))

# read in gene and sample meta data
geneDf <- read.csv("../data/BrainSpan/rows_metadata.filtered.csv", stringsAsFactor=F)
sampDf <- read.csv("../data/BrainSpan/columns_metadata.csv", stringsAsFactors=F)
# only keep samples from cortex regions
sampDf <- subset(sampDf,Identifier %in% colnames(expDf)) 

# read in age -> developmental stage/time point mapping
# then assign each sample to a time point
ageDf <- read.csv("../data/BrainSpan/age.mapping.csv", stringsAsFactors=F)
sampDf$TimePoint <- plyr::mapvalues(sampDf$age, ageDf$Age, ageDf$Developmental.Period)

# store number of samples for each time point
sampCounts <- as.data.frame(table(sampDf$TimePoint))
colnames(sampCounts) <- c('TimePoint','NumSamples')


# ----------------------------------------------------------------------------------------
# read in gene sets

# SCHEMA: FDR < 3.7e-3, 0.05, 0.25, 0.5 genes
schemaDf <- read.csv('../data/SCHEMA2020_TableS5.csv')
schemaDf <- subset(schemaDf,!is.na(P.meta))


# Bait selection: Sets 1-3 genes
setDf <- read.table('../data/SCZ_BaitSelectionGeneList.txt',
	header=T,sep='\t',stringsAsFactors=F)


# PPI: full combined network genes
intDf <- read.table('../data/SCZ_MasterInteractorTable_full.txt',
	header=T,sep='\t',stringsAsFactors=F)
intDf <- subset(intDf,ListName=='all_combined' & IsInteractor)

# non-interactors
nonIntDf <- read.table('../data/SCZ_MasterInteractorTable_full.txt',
	header=T,sep='\t',stringsAsFactors=F)
nonIntDf <- subset(nonIntDf,ListName=='all_combined' & !IsInteractor)


# PGC3 SCZ GWAS genes in LD loci (medRxiv manuscript Supp Table 3)
# 270 regions from combined meta-analysis (r2 > 0.1 SNPs +- 50kb)
pgcDf <- read.table('../data/PGC3_SCZ_combined_270loci_GenePosPvalue.txt',
	header=T,sep='\t',stringsAsFactors=F)
# remove MHC region genes (except SYNGAP1)
pgcDf <- subset(pgcDf,SNP!='rs13195636' | Gene=='SYNGAP1')
# remove all baits
pgcDf <- subset(pgcDf, !Gene %in% c('CACNA1C','CUL3','HCN1','RIMS1','SYNGAP1','TCF4'))

# PGC3 FINEMAP/SMR genes (medRxiv manuscript Supp Table 20)
# only keep protein-coding genes included in PGC3 loci above
fineDf <- read.table('../data/PGC3_SCZ_PrioritizedGenes.txt',
	header=T,sep='\t',stringsAsFactors=F)
fineGenes <- subset(fineDf,Symbol.ID %in% pgcDf$Gene)$Symbol.ID

# FINEMAP genes
mapGenes <- subset(fineDf,Symbol.ID %in% pgcDf$Gene & 
	(nonsynPP0.10==1 | UTRPP0.10==1 | k3.5singleGene==1))$Symbol.ID

# SMR genes
smrGenes <- subset(fineDf,Symbol.ID %in% pgcDf$Gene & 
	(SMRmap==1 | SMRsingleGene==1 | HI.C.SMR==1))$Symbol.ID

# PGC3 Network genes (genes in PPI network)
netGenes <- intersect(pgcDf$Gene,intDf$GeneName)

# PGC3 Non-interactor genes
nonGenes <- intersect(pgcDf$Gene,nonIntDf$GeneName)

# PGC3 Overlap genes (FINEMAP/SMR & Network)
overGenes <- intersect(fineGenes,netGenes)

# random genes (placeholder, to be added below)
randGenes <- c()


# named list of lists containing each gene set
geneSets <- list(
	'SCHEMA EWS' = schemaDf$Gene.Symbol[schemaDf$P.meta < 2.14e-6],
	'SCHEMA 5% FDR' = schemaDf$Gene.Symbol[schemaDf$Q.meta < 0.05],
	'SCHEMA 25% FDR' = schemaDf$Gene.Symbol[schemaDf$Q.meta < 0.25],
	'SCHEMA 50% FDR' = schemaDf$Gene.Symbol[schemaDf$Q.meta < 0.5],
	
	'Set 1' = setDf$Gene[setDf$Step1],
	'Set 2' = setDf$Gene[setDf$Step2],
	'Set 3' = setDf$Gene[setDf$Step3],
	'Set Random' = randGenes,
	
	'PPI Index' = c('CACNA1C','CUL3','HCN1','RIMS1','SYNGAP1','TCF4'),
	'PPI Network' = intDf$GeneName,
	'PPI Random' = randGenes,
	
	'PGC3' = pgcDf$Gene,
	'PGC3 FINEMAP' = mapGenes,
	'PGC3 SMR' = smrGenes,
	'PGC3 Network' = netGenes,
	'PGC3 Non-interactor' = nonGenes,
	'PGC3 Overlap' = overGenes,
	'PGC3 Random' = randGenes,
	
	# individual SCZ genes
	# SCHEMA EWS genes (SETD1A not found in BrainSpan data)
	'CUL1'=c('CUL1'),'XPO7'=c('XPO7'),'TRIO' = c('TRIO'),
	'CACNA1G'=c('CACNA1G'),'SP4'=c('SP4'),'GRIA3' = c('GRIA3'),
	'GRIN2A' = c('GRIN2A'),'HERC1'=c('HERC1'),'RB1CC1'=c('RB1CC1'),

	# SCZ index genes
	'CACNA1C'=c('CACNA1C'),'CUL3'=c('CUL3'),'HCN1'=c('HCN1'),
	'RIMS1'=c('RIMS1'),'SYNGAP1'=c('SYNGAP1'),'TCF4'=c('TCF4'),

	# Set3 genes not included as index genes or SCHEMA genes
	'CACNB2'=c('CACNB2'),'CSMD1'=c('CSMD1'),'SATB2'=c('SATB2'),'ZNF804A'=c('ZNF804A'))


# filter gene sets containing only genes with available BrainSpan data
filteredSets <- lapply(geneSets, function(x) {x[x %in% geneDf$gene_symbol]})


# random genes (matching # of genes in in Set 1, PPI Network, or PGC3)
set.seed(1234)
filteredSets$'Set Random' <- sample(geneDf$gene_symbol,length(filteredSets$'Set 1'))

set.seed(1234)
filteredSets$'PPI Random' <- sample(geneDf$gene_symbol,length(filteredSets$'PPI Network'))

set.seed(1234)
filteredSets$'PGC3 Random' <- sample(geneDf$gene_symbol,length(filteredSets$PGC3))


# check # of genes in each gene set pre- and post-filter
geneCounts <- merge(plyr::ldply(geneSets, .fun=length, .id='Set'),
	plyr::ldply(filteredSets, .fun=length, .id='Set'),by='Set')
colnames(geneCounts) <- c('Set','PreGenes','NumGenes')
geneCounts <- geneCounts[order(geneCounts$Set),]
rownames(geneCounts) <- NULL

geneCounts


# ----------------------------------------------------------------------------------------
# extract expression values for each gene set at each time point

# iterate through each gene set
meltDf <- NULL
for (s in names(filteredSets)) {
	tempDf <- expDf[geneDf$row_num[geneDf$gene_symbol %in% filteredSets[[s]]],]
	tempDf <- melt(tempDf,id.vars='gene')
	tempDf$Set <- s
	
	meltDf <- rbind(meltDf,tempDf)
}

meltDf$TimePoint <- plyr::mapvalues(meltDf$variable, sampDf$Identifier, sampDf$TimePoint)


# ----------------------------------------------------------------------------------------
# calculate median, SD, and SE expression of each gene set at each time point

# median
tpMedian <- melt(dcast(meltDf, Set~TimePoint, value.var='value', fun.aggregate=median), 
	id.vars='Set', variable.name='TimePoint', value.name='Median')

# standard deviation
tpSd <- melt(dcast(meltDf, Set~TimePoint, value.var='value', fun.aggregate=sd),
	id.vars='Set', variable.name='TimePoint', value.name='SD')

# number of data points (# genes x # samples) used at each time point
tpCounts <- melt(dcast(meltDf, Set~TimePoint, value.var='value', fun.aggregate=length),
	id.vars='Set', variable.name='TimePoint', value.name='N')

# merge stats into one data frame
mergeDf <- merge(merge(tpMedian, tpSd, by=c('Set','TimePoint')), tpCounts, by=c('Set','TimePoint'))

# calculate SE: standard deviation divided by sqrt(number of data points)
mergeDf$SE <- mergeDf$SD/sqrt(mergeDf$N)

# add number of genes and samples as columns
mergeDf$nSamples <- plyr::mapvalues(mergeDf$TimePoint,sampCounts$TimePoint,sampCounts$NumSamples)
mergeDf$nGenes <- plyr::mapvalues(mergeDf$Set,geneCounts$Set,geneCounts$NumGenes)


# order table by time points (chronologically) and gene sets
mergeDf$Set <- factor(mergeDf$Set, levels=names(geneSets))
mergeDf$TimePoint <- factor(mergeDf$TimePoint, levels=unique(ageDf$Developmental.Period))
mergeDf <- mergeDf[order(mergeDf$Set,mergeDf$TimePoint),]


# OUTPUT RESULTS TABLE
write.table(mergeDf,'../output/SCZ_BrainSpanStats.txt',row.names=F,quote=F,sep='\t')


##########################################################################################
# generate BrainSpan plot for Sets 1-3

# data to for line plot
lineSets <- c('Set Random','Set 1','Set 2','Set 3','SCHEMA EWS')
lineDf <- subset(mergeDf,Set %in% lineSets)
lineDf$Set <- factor(lineDf$Set, levels=lineSets)

# data for shading
fdrSets <- c('SCHEMA EWS','SCHEMA 5% FDR','SCHEMA 25% FDR','SCHEMA 50% FDR')
fdrDf <- subset(mergeDf,Set %in% fdrSets)
fdrDf$Set <- factor(fdrDf$Set, levels=fdrSets)

fdrDf$fdr <- 1 # FDR grouping for color scale
fdrDf$fdr[fdrDf$Set=='SCHEMA 5% FDR'] <- 2
fdrDf$fdr[fdrDf$Set=='SCHEMA 25% FDR'] <- 3
fdrDf$fdr[fdrDf$Set=='SCHEMA 50% FDR'] <- 4

fdrDf$ribbon.min <- 6
fdrDf$ribbon.min[fdrDf$Set=='SCHEMA EWS'] <- fdrDf$Median[fdrDf$Set=='SCHEMA 5% FDR']
fdrDf$ribbon.min[fdrDf$Set=='SCHEMA 5% FDR'] <- fdrDf$Median[fdrDf$Set=='SCHEMA 25% FDR']
fdrDf$ribbon.min[fdrDf$Set=='SCHEMA 25% FDR'] <- fdrDf$Median[fdrDf$Set=='SCHEMA 50% FDR']
# code above requires table is sorted so time points for each gene set are in same order


pdf("../output/SCZ_Sets1-3_BrainSpanPlot.pdf",height=4,width=5)

ggplot(lineDf,aes(x=TimePoint,y=Median)) + 

# plot SCHEMA genes with GWS, FDR < 0.05, 0.25, and 0.5 as shaded regions
geom_ribbon(data=fdrDf,mapping=aes(x=TimePoint,ymax=Median,ymin=ribbon.min,
	fill=Set,group=Set),linetype='blank') +
scale_fill_manual(name='SCHEMA FDR',
	labels=c('3.7e-3','0.05','0.25','0.5'),
	values=c('grey55','grey70','grey85','grey100')) +
                      
# plot median +- SE of gene sets across time points 
geom_line(aes(color=Set,group=Set),size=0.5) + geom_point(aes(color=Set,group=Set)) +
geom_errorbar(aes(ymin=Median-SE,ymax=Median+SE,color=Set,group=Set),size=0.5,width=0.2) +
scale_color_manual(name='Gene set',
	labels=c('Random','Set 1','Set 2','Set 3','SCHEMA'),
	values=c('dodgerblue','goldenrod1','sienna2','red2','black')) +

# pre vs. postnatal line
geom_vline(xintercept=4.5,linetype='dashed') +
annotate("text",x=4.6,y=10.2,hjust=0,size=4,label='Birth') +

xlab("Developmental stage") +
ylab("Median expression (± SE)") +
theme_classic() +
theme(legend.position='bottom',legend.box='vertical',legend.margin=margin(),
	legend.title=element_text(size=8),legend.text=element_text(size=8),
	axis.text.x=element_text(angle=45,hjust=1)) +
guides(fill=guide_legend(order=2,override.aes=list(linetype=1,color='black')),
	color=guide_legend(order=1))

dev.off()


##########################################################################################
# generate BrainSpan plot for PPI network genes

# data to for line plot
lineSets <- c('PPI Random','PPI Index','PPI Network','SCHEMA EWS')
lineDf <- subset(mergeDf,Set %in% lineSets)
lineDf$Set <- factor(lineDf$Set, levels=lineSets)

# data for shading
fdrSets <- c('SCHEMA EWS','SCHEMA 5% FDR','SCHEMA 25% FDR','SCHEMA 50% FDR')
fdrDf <- subset(mergeDf,Set %in% fdrSets)
fdrDf$Set <- factor(fdrDf$Set, levels=fdrSets)

fdrDf$fdr <- 1 # FDR grouping for color scale
fdrDf$fdr[fdrDf$Set=='SCHEMA 5% FDR'] <- 2
fdrDf$fdr[fdrDf$Set=='SCHEMA 25% FDR'] <- 3
fdrDf$fdr[fdrDf$Set=='SCHEMA 50% FDR'] <- 4

fdrDf$ribbon.min <- 6
fdrDf$ribbon.min[fdrDf$Set=='SCHEMA EWS'] <- fdrDf$Median[fdrDf$Set=='SCHEMA 5% FDR']
fdrDf$ribbon.min[fdrDf$Set=='SCHEMA 5% FDR'] <- fdrDf$Median[fdrDf$Set=='SCHEMA 25% FDR']
fdrDf$ribbon.min[fdrDf$Set=='SCHEMA 25% FDR'] <- fdrDf$Median[fdrDf$Set=='SCHEMA 50% FDR']
# code above requires table is sorted so time points for each gene set are in same order


pdf("../output/SCZ_Network_BrainSpanPlot.pdf",height=4,width=5)

ggplot(lineDf,aes(x=TimePoint,y=Median)) + 

# plot SCHEMA genes with GWS, FDR < 0.05, 0.25, and 0.5 as shaded regions
geom_ribbon(data=fdrDf,mapping=aes(x=TimePoint,ymax=Median,ymin=ribbon.min,
	fill=Set,group=Set),linetype='blank') +
scale_fill_manual(name='SCHEMA FDR',
	labels=c('3.7e-3','0.05','0.25','0.5'),
	values=c('grey55','grey70','grey85','grey100')) +
                      
# plot median +- SE of gene sets across time points 
geom_line(aes(color=Set,group=Set),size=0.5) + geom_point(aes(color=Set,group=Set)) +
geom_errorbar(aes(ymin=Median-SE,ymax=Median+SE,color=Set,group=Set),size=0.5,width=0.2) +
scale_color_manual(name='Gene set',
	labels=c('Random','Index','Network','SCHEMA'),
	values=c('dodgerblue','red2','purple','black')) +

# pre vs. postnatal line
geom_vline(xintercept=4.5,linetype='dashed') +
annotate("text",x=4.6,y=10.6,hjust=0,size=4,label='Birth') +

xlab("Developmental stage") +
ylab("Median expression (± SE)") +
theme_classic() +
theme(legend.position='bottom',legend.box='vertical',legend.margin=margin(),
	legend.title=element_text(size=8),legend.text=element_text(size=8),
	axis.text.x=element_text(angle=45,hjust=1)) +
guides(fill=guide_legend(order=2,override.aes=list(linetype=1,color='black')),
	color=guide_legend(order=1))

dev.off()


##########################################################################################
# generate BrainSpan plot for PGC3 gene sets

# data to for line plot
lineSets <- c('PGC3 Random','PGC3','PGC3 FINEMAP','PGC3 SMR','PGC3 Network','PGC3 Overlap','SCHEMA EWS')
lineDf <- subset(mergeDf,Set %in% lineSets)
lineDf$Set <- factor(lineDf$Set, levels=lineSets)

# data for shading
fdrSets <- c('SCHEMA EWS','SCHEMA 5% FDR','SCHEMA 25% FDR','SCHEMA 50% FDR')
fdrDf <- subset(mergeDf,Set %in% fdrSets)
fdrDf$Set <- factor(fdrDf$Set, levels=fdrSets)

fdrDf$fdr <- 1 # FDR grouping for cUnionolor scale
fdrDf$fdr[fdrDf$Set=='SCHEMA 5% FDR'] <- 2
fdrDf$fdr[fdrDf$Set=='SCHEMA 25% FDR'] <- 3
fdrDf$fdr[fdrDf$Set=='SCHEMA 50% FDR'] <- 4

fdrDf$ribbon.min <- 6
fdrDf$ribbon.min[fdrDf$Set=='SCHEMA EWS'] <- fdrDf$Median[fdrDf$Set=='SCHEMA 5% FDR']
fdrDf$ribbon.min[fdrDf$Set=='SCHEMA 5% FDR'] <- fdrDf$Median[fdrDf$Set=='SCHEMA 25% FDR']
fdrDf$ribbon.min[fdrDf$Set=='SCHEMA 25% FDR'] <- fdrDf$Median[fdrDf$Set=='SCHEMA 50% FDR']
# code above requires table is sorted so time points for each gene set are in same order


pdf("../output/SCZ_PGC3_BrainSpanPlot.pdf",height=4.5,width=5)

ggplot(lineDf,aes(x=TimePoint,y=Median)) + 

# plot SCHEMA genes with GWS, FDR < 0.05, 0.25, and 0.5 as shaded regions
geom_ribbon(data=fdrDf,mapping=aes(x=TimePoint,ymax=Median,ymin=ribbon.min,
	fill=Set,group=Set),linetype='blank') +
scale_fill_manual(name='SCHEMA FDR',
	labels=c('3.7e-3','0.05','0.25','0.5'),
	values=c('grey55','grey70','grey85','grey100')) +
                      
# plot median +- SE of gene sets across time points 
geom_line(aes(color=Set,group=Set),size=0.5) + geom_point(aes(color=Set,group=Set)) +
geom_errorbar(aes(ymin=Median-SE,ymax=Median+SE,color=Set,group=Set),size=0.5,width=0.2) +
scale_color_manual(name='Gene set',
	labels=c('Random','PGC3','FINEMAP','SMR','Network','Overlap','SCHEMA'),
	values=c('dodgerblue','goldenrod1','sienna2','darkgreen',
		'darkorchid4','magenta','black')) +

# pre vs. postnatal line
geom_vline(xintercept=4.5,linetype='dashed') +
annotate("text",x=4.6,y=10,hjust=0,size=4,label='Birth') +

xlab("Developmental stage") +
ylab("Median expression (± SE)") +
theme_classic() +
theme(legend.position='bottom',legend.box='vertical',legend.margin=margin(),
	legend.title=element_text(size=8),legend.text=element_text(size=8),
	axis.text.x=element_text(angle=45,hjust=1)) +
guides(fill=guide_legend(order=2,override.aes=list(linetype=1,color='black')),
	color=guide_legend(order=1,byrow=T))

dev.off()


##########################################################################################
# generate BrainSpan plots for individual SCZ genes

pdf("../output/SCZ_IndividualGenes_BrainSpanPlots.pdf",height=4,width=5)

# ----------------------------------------------------------------------------------------
# (1) SCHEMA genes

# data to for line plot
lineSets <- c('CUL1','XPO7','TRIO','CACNA1G','SP4','GRIA3','GRIN2A','HERC1','RB1CC1')
lineDf <- subset(mergeDf,Set %in% lineSets)
lineDf$Gene <- factor(lineDf$Set, levels=lineSets)


ggplot(lineDf,aes(x=TimePoint,y=Median)) + 

# plot SCHEMA genes with GWS, FDR < 0.05, 0.25, and 0.5 as shaded regions
geom_ribbon(data=fdrDf,mapping=aes(x=TimePoint,ymax=Median,ymin=ribbon.min,
	fill=Set,group=Set),linetype='blank') +
scale_fill_manual(name='SCHEMA FDR',
	labels=c('3.7e-3','0.05','0.25','0.5'),
	values=c('grey55','grey70','grey85','grey100')) +
                   
# plot median +- SE for each gene across time points 
geom_line(aes(color=Gene,group=Gene),size=0.5) + geom_point(aes(color=Gene,group=Gene)) +
geom_errorbar(aes(ymin=Median-SE,ymax=Median+SE,color=Gene,group=Gene),
	size=0.5,width=0.2) +

# pre vs. postnatal line
geom_vline(xintercept=4.5,linetype='dashed') +
annotate("text",x=4.6,y=11,hjust=0,size=4,label='Birth') +

xlab("Developmental stage") +
ylab("Median expression (± SE)") +
theme_classic() +
theme(legend.position='bottom',legend.box='vertical',legend.margin=margin(),
	legend.title=element_text(size=8),legend.text=element_text(size=8),
	axis.text.x=element_text(angle=45,hjust=1)) +
guides(fill=guide_legend(order=2,override.aes=list(linetype=1,color='black')),
	color=guide_legend(order=1))


# ----------------------------------------------------------------------------------------
# (2) SCZ index genes

# data to for line plot
lineSets <- c('CACNA1C','CUL3','HCN1','RIMS1','SYNGAP1','TCF4')
lineDf <- subset(mergeDf,Set %in% lineSets)
lineDf$Gene <- factor(lineDf$Set, levels=lineSets)
	

ggplot(lineDf,aes(x=TimePoint,y=Median)) + 

# plot SCHEMA genes with GWS, FDR < 0.05, 0.25, and 0.5 as shaded regions
geom_ribbon(data=fdrDf,mapping=aes(x=TimePoint,ymax=Median,ymin=ribbon.min,
	fill=Set,group=Set),linetype='blank') +
scale_fill_manual(name='SCHEMA FDR',
	labels=c('3.7e-3','0.05','0.25','0.5'),
	values=c('grey55','grey70','grey85','grey100')) +
                   
# plot median +- SE for Random, Steps 1-3, and top SCHEMA GWS genes across time points 
geom_line(aes(color=Gene,group=Gene),size=0.5) + geom_point(aes(color=Gene,group=Gene)) +
geom_errorbar(aes(ymin=Median-SE,ymax=Median+SE,color=Gene,group=Gene),
	size=0.5,width=0.2) +

# pre vs. postnatal line
geom_vline(xintercept=4.5,linetype='dashed') +
annotate("text",x=4.6,y=12.5,hjust=0,size=4,label='Birth') +

xlab("Developmental stage") +
ylab("Median expression (± SE)") +
theme_classic() +
theme(legend.position='bottom',legend.box='vertical',legend.margin=margin(),
	legend.title=element_text(size=8),legend.text=element_text(size=8),
	axis.text.x=element_text(angle=45,hjust=1)) +
guides(fill=guide_legend(order=2,override.aes=list(linetype=1,color='black')),
	color=guide_legend(order=1))


# ----------------------------------------------------------------------------------------
# (3) Set 3 genes not included as index genes in final network

# data to for line plot
lineSets <- c('CACNB2','CSMD1','GRIN2A','SATB2','ZNF804A')
lineDf <- subset(mergeDf,Set %in% lineSets)
lineDf$Gene <- factor(lineDf$Set, levels=lineSets)

ggplot(lineDf,aes(x=TimePoint,y=Median)) + 

# plot SCHEMA genes with GWS, FDR < 0.05, 0.25, and 0.5 as shaded regions
geom_ribbon(data=fdrDf,mapping=aes(x=TimePoint,ymax=Median,ymin=ribbon.min,
	fill=Set,group=Set),
	linetype='blank') +
scale_fill_manual(name='SCHEMA FDR',
	labels=c('3.7e-3','0.05','0.25','0.5'),
	values=c('grey55','grey70','grey85','grey100')) +
                   
# plot median +- SE for Random, Steps 1-3, and top SCHEMA GWS genes across time points 
geom_line(aes(color=Gene,group=Gene),size=0.5) + geom_point(aes(color=Gene,group=Gene)) +
geom_errorbar(aes(ymin=Median-SE,ymax=Median+SE,color=Gene,group=Gene),
	size=0.5,width=0.2) +

# pre vs. postnatal line
geom_vline(xintercept=4.5,linetype='dashed') +
annotate("text",x=4.6,y=13,hjust=0,size=4,label='Birth') +

xlab("Developmental stage") +
ylab("Median expression (± SE)") +
theme_classic() +
theme(legend.position='bottom',legend.box='vertical',legend.margin=margin(),
	legend.title=element_text(size=8),legend.text=element_text(size=8),
	axis.text.x=element_text(angle=45,hjust=1)) +
guides(fill=guide_legend(order=2,override.aes=list(linetype=1,color='black')),
	color=guide_legend(order=1))


dev.off()
