##########################################################################################
## pLI score enrichment analysis for various SCZ gene sets using KS test
##
## Author: Yu-Han Hsu
##########################################################################################

rm(list=ls())

library(genoppi) # for pLI scores in gnomad_table
library(ggplot2)
library(plyr)


# ----------------------------------------------------------------------------------------
# gnomAD pLI scores
pliDf <- unique(subset(gnomad_table,!is.na(gnomad_table$pLI)))
scoreDf <- data.frame(variable='gnomAD pLI',gene=pliDf$gene,score=pliDf$pLI)
	
dim(scoreDf)


# ----------------------------------------------------------------------------------------
# gene sets to test

# SCHEMA genes
schemaDf <- read.csv('../data/SCHEMA2020_TableS5.csv') 
schemaDf <- subset(schemaDf,!is.na(P.meta))


# Sets 1-3 from GWAS data
setDf <- read.table('../data/SCZ_BaitSelectionGeneList.txt',
	header=T,sep='\t',stringsAsFactors=F)


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


# FINEMAP genes and loci containing FINEMAP genes
mapGenes <- subset(fineDf,Symbol.ID %in% pgcDf$Gene & 
	(nonsynPP0.10==1 | UTRPP0.10==1 | k3.5singleGene==1))$Symbol.ID
mapSnps <- unique(subset(pgcDf,Gene %in% mapGenes)$SNP)
mapLociGenes <- subset(pgcDf,SNP %in% mapSnps)$Gene


# SMR genes and loci containing SMR genes
smrGenes <- subset(fineDf,Symbol.ID %in% pgcDf$Gene & 
	(SMRmap==1 | SMRsingleGene==1 | HI.C.SMR==1))$Symbol.ID
smrSnps <- unique(subset(pgcDf,Gene %in% smrGenes)$SNP)
smrLociGenes <- subset(pgcDf,SNP %in% smrSnps)$Gene


# PGC3 Network genes (full PPI network, different from stringent network read in above)
netDf <- read.table('../data/SCZ_MasterInteractorTable_full.txt',
	header=T,sep='\t',stringsAsFactors=F)
netGenes <- intersect(pgcDf$Gene,
	subset(netDf,ListName=='all_combined' & IsInteractor)$GeneName)


# genes in loci containing Network genes
netSnps <- unique(subset(pgcDf,Gene %in% netGenes)$SNP)
netLociGenes <- subset(pgcDf,SNP %in% netSnps)$Gene


# PGC3 Non-interactor genes
nonGenes <- intersect(pgcDf$Gene,
	subset(netDf,ListName=='all_combined' & !IsInteractor)$GeneName)


# genes in loci containing Non-interactor genes
nonSnps <- unique(subset(pgcDf,Gene %in% nonGenes)$SNP)
nonLociGenes <- subset(pgcDf,SNP %in% nonSnps)$Gene


# PGC3 Overlap genes (FINEMAP/SMR & Network)
overGenes <- intersect(fineGenes,netGenes)


# genes in loci containing Overlap genes
overSnps <- unique(subset(pgcDf,Gene %in% overGenes)$SNP)
overLociGenes <- subset(pgcDf,SNP %in% overSnps)$Gene


# named list of lists containing each gene set
geneSets <- list(
	'SCHEMA EWS' = schemaDf$Gene.Symbol[schemaDf$P.meta < 2.14e-6],
	'SCHEMA 5% FDR' = schemaDf$Gene.Symbol[schemaDf$Q.meta < 0.05],
	
	'Set 1' = setDf$Gene[setDf$Step1],
	'Set 2' = setDf$Gene[setDf$Step2],
	'Set 3' = setDf$Gene[setDf$Step3],
	
	'PGC3' = pgcDf$Gene,
	'FINEMAP' = mapGenes,
	'SMR' = smrGenes,
	'Network' = netGenes,
	'Non-interactor' = nonGenes,
	'Overlap' = overGenes,
	
	'FINEMAP loci' = mapLociGenes,
	'SMR loci' = smrLociGenes,
	'Network loci' = netLociGenes,
	'Non-interactor loci' = nonLociGenes,
	'Overlap loci' = overLociGenes)

# check # of genes in each gene set
plyr::ldply(geneSets, .fun=length, .id='Set')


# ----------------------------------------------------------------------------------------
# table containing pairs of gene sets to use for CONDITONAL tests

# gene sets that are subset and superset of each other
subSuperSets <- data.frame(
	Subset=c('Set 2','Set 3','Set 3',
		'FINEMAP','SMR','Network','Non-interactor','Overlap',
		'FINEMAP','SMR','Network','Non-interactor','Overlap',
		'Overlap','Overlap','Overlap'),
	Superset=c('Set 1','Set 1','Set 2',
		'PGC3','PGC3','PGC3','PGC3','PGC3',
		'FINEMAP loci','SMR loci','Network loci','Non-interactor loci','Overlap loci',
		'FINEMAP','SMR','Network'),
	stringsAsFactors=F)
	

# other gene set pairs (for testing if they're different from each other)
compareSets <- data.frame(
	X=c('Set 3','Set 3','FINEMAP',
		'Network','Network','Network',
		'Overlap','Overlap'),
	Y=c('SCHEMA EWS','SCHEMA 5% FDR','SMR',
		'FINEMAP','SMR','Non-interactor',
		'SCHEMA EWS','SCHEMA 5% FDR'),
	stringsAsFactors=F)


# ----------------------------------------------------------------------------------------
# enrichment analysis

resultsTable <- NULL # enrichment results
plotDf <- NULL # store scores for plotting cumulative density

for (v in unique(scoreDf$variable)) {

	geneDf <- subset(scoreDf,variable==v)
	
	# KS test direction
	testDir <- NULL
	if (v=='gnomAD pLI') {testDir <- 'less'}
	else {testDir <- 'greater'}
	# 'less': score CDF of x lies BELOW that of y (so x are > y)
	# 'greater': score CDF of x lies ABOVE that of y (so x are < y)
	

	# ------------------------------------------------------------------------------------
	# calculate GLOBAL enrichment for gene sets (vs. genomic background)

	for (s in names(geneSets)) {

		intScores <- geneDf$score[geneDf$gene %in% geneSets[[s]]]
		genomeScores <- geneDf$score[!geneDf$gene %in% geneSets[[s]]]
	
		# KS test for genes in set vs. genomic background
		res <- ks.test(intScores,genomeScores,alternative=testDir)

		resultsTable <- rbind(resultsTable, data.frame(Variable=v,Test='Global',Set=s,
			IntCount=length(intScores),nonIntCount=length(genomeScores),pDir=testDir,
			Stat=res$statistic,Pvalue=res$p.value))
		
		# df for cumulative density plots
		plotDf <- rbind(plotDf, data.frame(Variable=v,Set=s,
			Group=c(rep('GenesInSet',length(intScores)),
			rep('Genome',length(genomeScores))),
			Score=c(intScores,genomeScores)))
	}


	# ------------------------------------------------------------------------------------
	# calculate CONDITONAL enrichment for gene sets
	# (subset vs. non-overlapping genes in superset)

	for (i in 1:nrow(subSuperSets)) {
	
		subGenes <- geneSets[[subSuperSets[i,'Subset']]]
		superGenes <- geneSets[[subSuperSets[i,'Superset']]]
		
		intScores <- geneDf$score[geneDf$gene %in% subGenes]
		nonIntScores <- geneDf$score[geneDf$gene %in% setdiff(superGenes,subGenes)]

		# KS test for genes in subset vs. superset
		res <- ks.test(intScores,nonIntScores,alternative=testDir)

		resultsTable <- rbind(resultsTable, data.frame(Variable=v,Test='Conditional',
			Set=paste(subSuperSets[i,'Subset'],'vs',subSuperSets[i,'Superset']),
			IntCount=length(intScores),nonIntCount=length(nonIntScores),pDir=testDir,
			Stat=res$statistic,Pvalue=res$p.value))
	
	}

	
	# ------------------------------------------------------------------------------------
	# compare gene sets that are not sub/supersets of each other
	# using two-sided KS test
	
	for (i in 1:nrow(compareSets)) {
	
		xGenes <- geneSets[[compareSets[i,'X']]]
		yGenes <- geneSets[[compareSets[i,'Y']]]
		
		intScores <- geneDf$score[geneDf$gene %in% setdiff(xGenes,yGenes)]
		nonIntScores <- geneDf$score[geneDf$gene %in% setdiff(yGenes,xGenes)]

		# KS test for genes in set vs. genomic background
		res <- ks.test(intScores,nonIntScores,alternative='two.sided')

		resultsTable <- rbind(resultsTable, data.frame(Variable=v,Test='Conditional',
			Set=paste(compareSets[i,'X'],'vs',compareSets[i,'Y']),
			IntCount=length(intScores),nonIntCount=length(nonIntScores),pDir='two.sided',
			Stat=res$statistic,Pvalue=res$p.value))
	
	}
	
}


# output enrichment results table
write.table(resultsTable,"../output/SCZ_pLiEnrichment.txt",
	quote=F,row.names=F,sep="\t")


# ----------------------------------------------------------------------------------------
# plot cumulative distributions of pLI scores for SCHEMA genes + Sets 1-3

plotSets <- c('Set 1','Set 2','Set 3','SCHEMA 5% FDR','SCHEMA EWS')
df <- subset(plotDf,
	(Variable=='gnomAD pLI' & Group=='GenesInSet' & Set %in% plotSets) |
	(Variable=='gnomAD pLI' & Group=='Genome' & Set=='Set 1'))

df$Set[df$Group=='Genome'] <- 'Genome'
df$Set <- factor(df$Set,levels=c('Genome',plotSets))


pdf("../output/SCZ_Sets1-3_pLiPlot.pdf",width=5,height=3)

ggplot(df,aes(x=Score,color=Set)) + stat_ecdf(geom='step',size=1) +
scale_color_manual(name='Gene set',
	labels=c('Genome','Set 1','Set 2','Set 3','SCHEMA 5% FDR','SCHEMA EWS'),
	values = c('dodgerblue','goldenrod1','sienna2','red2','grey55','black')) +

xlab('gnomAD pLI score') + ylab('Cumulative density') +
theme_classic()

dev.off()


# ----------------------------------------------------------------------------------------
# plot cumulative distributions of pLI scores for SCHEMA + PGC3 genes

plotSets <- c('PGC3','FINEMAP','SMR','Network','Overlap','SCHEMA 5% FDR','SCHEMA EWS')
df <- subset(plotDf,
	(Variable=='gnomAD pLI' & Group=='GenesInSet' & Set %in% plotSets) |
	(Variable=='gnomAD pLI' & Group=='Genome' & Set=='PGC3'))

df$Set[df$Group=='Genome'] <- 'Genome'
df$Set <- factor(df$Set,levels=c('Genome',plotSets))


pdf("../output/SCZ_PGC3_pLiPlot.pdf",width=5,height=3)

ggplot(df,aes(x=Score,color=Set)) + stat_ecdf(geom='step',size=1) +
scale_color_manual(name='Gene set',
	labels=c('Genome','PGC3','FINEMAP','SMR','Network','Overlap','SCHEMA 5% FDR','SCHEMA EWS'),
	values = c('dodgerblue','goldenrod1','sienna2','darkgreen',
		'darkorchid4','magenta','grey55','black')) +
	
xlab('gnomAD pLI score') + ylab('Cumulative density') +
theme_classic()

dev.off()

