##########################################################################################
## Rare variant risk enrichment analysis of SCZ PPI networks
## using gene-based scores and KS test
##
## Author: Yu-Han Hsu
##########################################################################################

rm(list=ls())

library(genoppi) # for pLI scores in gnomad_table
library(ggplot2)


# ----------------------------------------------------------------------------------------
# gene-based scores to use in KS tests

# gnomAD pLI scores
pliDf <- unique(subset(gnomad_table,!is.na(gnomad_table$pLI)))
scoreDf <- data.frame(variable='gnomAD pLI',gene=pliDf$gene,score=pliDf$pLI)

# SCZ/SCHEMA p-values (Singh et al.; columns to use: Gene.Symbol, P.meta, Q.meta)
schemaDf <- read.csv('../data/SCHEMA2020_TableS5.csv') 
schemaDf <- subset(schemaDf,!is.na(P.meta))
scoreDf <- rbind(scoreDf,data.frame(variable='SCHEMA P-value',
	gene=schemaDf$Gene.Symbol,score=schemaDf$P.meta))

# DD p-values (Kaplanis et al.; columns to use: symbol, denovoWEST_p_full)
ddDf <- read.table('../data/Kaplanis2019_TableS2.txt',header=T,sep='\t')
ddDf <- subset(ddDf,!is.na(denovoWEST_p_full))
scoreDf <- rbind(scoreDf,data.frame(variable='DD P-value',
	gene=ddDf$symbol,score=ddDf$denovoWEST_p_full))

# ASD q-values (Satterstrom et al.; columns to use: hugoGene, qval_dnccPTV)
asdDf <- read.table('../data/Satterstrom2020_TableS2_Autosomal.txt',header=T,sep='\t')
scoreDf <- rbind(scoreDf,data.frame(variable='ASD Q-value',
	gene=asdDf$hugoGene,score=asdDf$qval_dnccPTV))
	
dim(scoreDf)


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
# repeat anlaysis for pLI, SCHEMA, DD, and ASD scores

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
	# calculate CONDITIONAL enrichment for networks (vs. non-interactors)
	
	for (i in 1:length(networks)) {

		ints <- subset(intDf,ListName==networks[i] & significant)$gene
		nonInts <- subset(intDf,ListName==networks[i] & !significant)$gene
	
		intScores <- geneDf$score[geneDf$gene %in% ints]
		nonIntScores <- geneDf$score[geneDf$gene %in% nonInts]
		
		# KS test for ints vs. non-ints
		res <- ks.test(intScores,nonIntScores,alternative=testDir)
		
		resultsTable <- rbind(resultsTable,
			data.frame(Variable=v,Test='Conditional',Set=names[i],
			IntCount=length(intScores),nonIntCount=length(nonIntScores),pDir=testDir,
			Stat=res$statistic,Pvalue=res$p.value))
	}
	
}


# output enrichment results table
write.table(resultsTable,"../output/SCZ_RareVariantEnrichment.txt",
	quote=F,row.names=F,sep="\t")


# ----------------------------------------------------------------------------------------
# plot enrichment results

figDf <- subset(resultsTable,Set %in% names & Test=='Conditional')

# significance text labels
figDf$Sig <- ifelse(figDf$Pvalue < 0.05,"*","")
figDf$Sig[figDf$Pvalue<0.05/11] <- "**"

figDf$SigText <- ifelse(figDf$Pvalue < 0.05,
	paste(formatC(figDf$Stat,format="f",digits=2),figDf$Sig,sep=''),"")


# set network label and plotting order
figDf$IntCount <- sizeDf$IntCount[match(figDf$Set,sizeDf$Network)]
figDf$Set <- paste(figDf$Set,'\n(',figDf$IntCount,')',sep='')
netOrder <- paste(rev(names),
	'\n(',sizeDf$IntCount[match(rev(names),sizeDf$Network)],')',sep='')
figDf$Set <- factor(figDf$Set, levels=netOrder)


# set trait order
figDf$Score[figDf$Variable=='gnomAD pLI'] <- 'pLI'
figDf$Score[figDf$Variable=='SCHEMA P-value'] <- 'SCZ'
figDf$Score[figDf$Variable=='DD P-value'] <- 'DD'
figDf$Score[figDf$Variable=='ASD Q-value'] <- 'ASD'
figDf$Score <- factor(figDf$Score, levels=c('SCZ','DD','ASD','pLI'))


pdf("../output/SCZ_RareVariantEnrichment_HeatMap.pdf",height=3.5,width=4)

ggplot(figDf,aes(x=Score,y=Set,fill=-log10(Pvalue))) + 
geom_tile() + geom_text(aes(label=SigText),size=3) +
scale_fill_gradient(name=expression(paste(-log[10],"(P-value)",sep="")),low="white", high="red") +
xlab("Score") + ylab("Network") +
theme_bw() + 
theme(legend.key.size=unit(0.7,"line"),legend.title=element_text(size=9),
	legend.text=element_text(size=8),
	axis.title=element_text(size=9),axis.text=element_text(size=8))

dev.off()

