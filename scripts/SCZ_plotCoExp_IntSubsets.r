##########################################################################################
## Compare and plot co-expression z-scores of bait-interactor pairs subsetted by
## different variables, using co-expression data derived from 4 transcriptomic datasets:
## Stickels2021, Maynard2021, Velmeshev2019, BrainSpan
##
## Author: Yu-Han Hsu
##########################################################################################

rm(list=ls())

library(ggpubr)


# ----------------------------------------------------------------------------------------
# read in bait-interactor pairs with subsetting info

intDf <- read.table('../data/SCZ_IntSubsetTable.txt',header=T,sep='\t',stringsAsFactors=F)

# bin by time point
intDf$TimePoint <- factor(intDf$TimePoint,levels=sort(unique(intDf$TimePoint)))

# bin by bait log2 FC
intDf$Bait_logFC_Cutoff <- ifelse(intDf$Bait_logFC <= 1,'0','1')
intDf$Bait_logFC_Cutoff[intDf$Bait_logFC > 1.5] <- '1.5'
intDf$Bait_logFC_Cutoff <- factor(intDf$Bait_logFC_Cutoff,levels=c('0','1','1.5'))

# bin by bait FDR
intDf$Bait_FDR_Cutoff <- ifelse(intDf$Bait_FDR <= 1e-4,'1e-4','1e-3')
intDf$Bait_FDR_Cutoff[intDf$Bait_FDR > 1e-3] <- '0.01'
intDf$Bait_FDR_Cutoff[intDf$Bait_FDR > 0.01] <- '0.1'
intDf$Bait_FDR_Cutoff <- factor(intDf$Bait_FDR_Cutoff,levels=c('1e-4','1e-3','0.01','0.1'))

# bin by interactor log2 FC
intDf$Int_logFC_Cutoff <- ifelse(intDf$Int_logFC <= 0.5,'0','0.5')
intDf$Int_logFC_Cutoff[intDf$Int_logFC > 1] <- '1'
intDf$Int_logFC_Cutoff[intDf$Int_logFC > 1.5] <- '1.5'
intDf$Int_logFC_Cutoff <- factor(intDf$Int_logFC_Cutoff,levels=c('0','0.5','1','1.5'))

# bin by interactor FDR
intDf$Int_FDR_Cutoff <- ifelse(intDf$Int_FDR <= 1e-4,'1e-4','1e-3')
intDf$Int_FDR_Cutoff[intDf$Int_FDR > 1e-3] <- '0.01'
intDf$Int_FDR_Cutoff[intDf$Int_FDR > 0.01] <- '0.1'
intDf$Int_FDR_Cutoff <- factor(intDf$Int_FDR_Cutoff,levels=c('1e-4','1e-3','0.01','0.1'))

# bin by recurrence (i.e. number of IPs detecting the interaction)
for (i in 1:nrow(intDf)) {
	intDf$Recurrence[i] <- sum(intDf$Bait==intDf$Bait[i] & 
		intDf$Interactor==intDf$Interactor[i])
}
intDf$Recurrence <- factor(intDf$Recurrence,levels=sort(unique(intDf$Recurrence)))


# ----------------------------------------------------------------------------------------
# read in co-expression data

coexpFiles <- c('../data/CoExp/Stickels2021','../data/CoExp/Maynard2021',
	'../data/CoExp/Velmeshev2019','../data/CoExp/BrainSpan')
metrics <- c('FisherNegLogP','FisherNegLogP','PropRho','PropRho')
titles <- c('Stickels 2021','Maynard 2021','Velmeshev 2019','BrainSpan')

plotDf <- NULL

for (i in 1:length(titles)) {
	print(titles[i])

	inFile <- gzfile(paste(coexpFiles[i],'_SCZ-Baits_CoExp',metrics[i],'.txt.gz',sep=''),'rt')
	
	metricDf <- read.table(inFile,header=T,sep='\t',
		stringsAsFactors=F,row.names=1,check.names=F)
		
	close(inFile)

	if (metrics[i]=='FisherNegLogP') {
		# fix any infinite numbers outputted by Python
		rowNames <- rownames(metricDf)
		metricDf[metricDf=='inf'] <- Inf
		metricDf <- data.frame(apply(metricDf,2,function(x) as.numeric(x)))
		rownames(metricDf) <- rowNames

		# set Inf values to just above max of non-Inf value (per bait/row)
		metricDf <- data.frame(t(apply(metricDf,1,function(x) {
			ifelse(is.infinite(x),ceiling(max(x[!is.infinite(x)])),x)})))
	}

	# rank-based inverse normal transformation to get z-scores (per bait/row)
	metricDf <- data.frame(t(apply(metricDf,1,function(x) {
		qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))})))

	# save z-scores for bait-gene pairs
	iDf <- subset(intDf, Bait %in% rownames(metricDf) & Interactor %in% names(metricDf))
	matchRows <- match(iDf$Bait,rownames(metricDf))
	matchCols <- match(iDf$Interactor,names(metricDf))
	for (k in 1:length(matchRows)) {
		iDf$Metric[k] <- metricDf[matchRows[k],matchCols[k]]
	}

	# filter out bait-bait pairs
	iDf <- subset(iDf,Bait!=Interactor)
	plotDf <- rbind(plotDf,data.frame(iDf,DataSource=titles[i]))

}


# ----------------------------------------------------------------------------------------
# plot and compare co-expression z-scores of ints subsetted by different variables

pdf('../output/SCZ_CoExp_IntSubsets_ViolinPlots.pdf',height=2.1,width=7.5)

# set plotting order
plotDf$DataSource <- factor(plotDf$DataSource,levels=titles)

variables <- c('TimePoint','Bait_logFC_Cutoff','Bait_FDR_Cutoff',
	'Int_logFC_Cutoff','Int_FDR_Cutoff','Recurrence')
varNames <- c('Differentiation time','Bait log2 FC','Bait FDR',
	'Interactor log2 FC','Interactor FDR','Number of IPs')

for (k in 1:length(variables)) {
	
	# only plot unique bait-int pairs if looking at recurrence
	if (variables[k]=='Recurrence') { 
		plotDf <- unique(plotDf[,c('Bait','Interactor','DataSource','Metric','Recurrence')])
	}

	# set Wilcoxon test parameters
	varPairs <- combn(levels(plotDf[,variables[k]]),2)
	testList <- NULL
	for (col in 1:ncol(varPairs)) { testList[[col]] <- varPairs[,col] }	
	pCutoffs <- c(0,0.05/length(testList),0.05,1)
	pSymbols <- c('**','*','ns')

	# plot score distributions
	plot <- ggviolin(plotDf,x=variables[k],y='Metric',color=variables[k]) +
		geom_boxplot(aes_string(color=variables[k]),width=0.2,outlier.size=0.2) +
		
		facet_wrap(~DataSource,ncol=4,scales='free_y') +

		# pairwise Wilcoxon tests
		stat_compare_means(comparisons=testList,method='wilcox.test',
		symnum.args=list(cutpoints=pCutoffs,symbols=pSymbols),hide.ns=T,size=3,vjust=0.5) +

		# label size of each group (!!sym for passing string variable as symbol in fn)
		stat_summary(aes(label=..y..,color=!!sym(variables[k]),
			y=stage(Metric,after_stat=min(iDf$Metric))),
			fun=length,geom='text',size=3,vjust=1.2,position=position_dodge(width=0.75)) +

		scale_color_brewer(palette='Set1') +
		scale_y_continuous(expand=expansion(mult=c(0.1,0.05),add=c(0,0))) +
		xlab(varNames[k]) + ylab(bquote("Co-expression"~italic(Z)*"-score")) +
		theme_bw() + theme(legend.position='none',axis.text=element_text(size=8))

	print(plot)
}

dev.off()

