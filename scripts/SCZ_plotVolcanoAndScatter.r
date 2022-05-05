##########################################################################################
## Generate volcano and scatter plots for SCZ IP-MS datasets
##
## Author: Yu-Han Hsu
##########################################################################################

rm(list=ls())

library(readxl)
library(ggplot2)
library(ggrepel)
library(genoppi)


# ----------------------------------------------------------------------------------------
# read in Genoppi analysis results each IP

df <- NULL
sheetList <- excel_sheets("../data/SCZ_GenoppiResults.xlsx")
summaryDf <- NULL

for (i in 2:length(sheetList)) {
	tempDf <- read_excel("../data/SCZ_GenoppiResults.xlsx",sheet=sheetList[i])
	tempDf$IP <- sheetList[i]

	# index protein	
	bait <- strsplit(sheetList[i],'_')[[1]][1]
	tempDf$bait <- bait
	
	# define significant proteins
	tempDf$significant = tempDf$logFC>0 & tempDf$FDR<=0.1
	
	# InWeb interactors
	tempInwebDf <- get_inweb_list(bait)
	tempDf$IsInWeb <- tempDf$gene %in% subset(tempInwebDf,significant)$gene
	
	# save data to master tables for plotting
	df <- rbind(df,tempDf)
}


# ----------------------------------------------------------------------------------------
# Scatter and volcano plot for CACNA1C_wk3_1
# InWeb + calcium channel subunits overlay

subDf <- subset(df,IP=='CACNA1C_wk3_1')

pdf("../output/SCZ_CACNA1C_ScatterVolcano.pdf",height=2.5,width=2.5)

# VOLCANO
set.seed(42)

ggplot(subDf,aes(x=logFC,y=-log10(pvalue),color=significant)) +
xlab(bquote(log[2]*"(Fold change)")) + ylab(bquote(-log[10]*"(P-value)")) +
geom_vline(xintercept=0,linetype='dashed') +

# plot all proteins (green = significant, blue = not significant)
geom_point(size=1.5,stroke=0.3,
	color=ifelse(subDf$significant,"springgreen3","dodgerblue")) +

# label InWeb interactors (yellow = significant, empty black circle = not significant)
geom_point(subset(subDf,IsInWeb & significant),
	mapping=aes(x=logFC,y=-log10(pvalue)),size=1.5,stroke=0.3,color="yellow") + 
	
# label calcium channel subunits
geom_point(subset(subDf,grepl('CACN',subDf$gene)),
	mapping=aes(x=logFC,y=-log10(pvalue)),size=1.5,stroke=0.3,color="sienna1") +

geom_point(subset(subDf,IsInWeb | grepl('CACN',subDf$gene)),
	mapping=aes(x=logFC,y=-log10(pvalue)),size=1.5,stroke=0.3,color="black",shape=1) +

# label bait (red = signficant, orange = not significant)
geom_point(subset(subDf,gene==bait),
	mapping=aes(x=logFC,y=-log10(pvalue)),size=2,stroke=0.3,color="red") + 	
geom_point(subset(subDf,gene==bait),
	mapping=aes(x=logFC,y=-log10(pvalue)),size=2,stroke=0.3,color="black",shape=1) +	

geom_text_repel(subset(subDf,gene==bait | grepl('CACN',subDf$gene)),
	mapping=aes(label=gene),arrow=arrow(length=unit(0.015,'npc')),
	box.padding=unit(0.1,"lines"),point.padding=unit(0.15,"lines"),color="black",size=2) +

ggtitle('CACNA1C, week 3') + theme_classic() +
theme(plot.title=element_text(size=8),axis.title=element_text(size=8),axis.text=element_text(size=8))


# SCATTER
set.seed(42)

ggplot(subDf,aes(x=rep1,y=rep2,color=significant)) +
xlab('Replicate 1') + ylab('Replicate 2') +

# plot all proteins (green = significant, blue = not significant)
geom_point(size=1.5,stroke=0.3,
	color=ifelse(subDf$significant,"springgreen3","dodgerblue")) +

# label InWeb interactors (yellow = significant, empty black circle = not significant)
geom_point(subset(subDf,IsInWeb & significant),
	mapping=aes(x=rep1,y=rep2),size=1.5,stroke=0.3,color="yellow") + 
	
# label calcium channel subunits
geom_point(subset(subDf,grepl('CACN',subDf$gene)),
	mapping=aes(x=rep1,y=rep2),size=1.5,stroke=0.3,color="sienna1") +

geom_point(subset(subDf,IsInWeb | grepl('CACN',subDf$gene)),
	mapping=aes(x=rep1,y=rep2),size=1.5,stroke=0.3,color="black",shape=1) +

# label bait (red = signficant, orange = not significant)
geom_point(subset(subDf,gene==bait),
	mapping=aes(x=rep1,y=rep2),size=2,stroke=0.3,color="red") + 	
geom_point(subset(subDf,gene==bait),
	mapping=aes(x=rep1,y=rep2),size=2,stroke=0.3,color="black",shape=1) +	

geom_text_repel(subset(subDf,gene==bait | grepl('CACN',subDf$gene)),
	mapping=aes(label=gene),arrow=arrow(length=unit(0.015,'npc')),
	box.padding=unit(0.1,"lines"),point.padding=unit(0.15,"lines"),color="black",size=2) +

ggtitle('CACNA1C, week 3') + theme_classic() + 
theme(plot.title=element_text(size=8),axis.title=element_text(size=8),axis.text=element_text(size=8))


dev.off()


# ----------------------------------------------------------------------------------------
# Volcano plots for 6 representative IPs (one for each index protein)

ipList = c("CACNA1C_wk4_2","CUL3_wk7","HCN1_wk7_1",
	"RIMS1_wk4_1","SYNGAP1_wk4_1","TCF4_wk1")

ipNames = c("CACNA1C, week 4","CUL3, week 7","HCN1, week 7",
	"RIMS1, week 4","SYNGAP1, week 4","TCF4, week 1")

pdf("../output/SCZ_InWebVolcanos_6IPs.pdf",height=2.5,width=2.5)

for (i in 1:length(ipList)) {

	subDf <- subset(df,IP==ipList[i])
	
	# print 2 versions with different seeds to get good ggrepel labeling
	set.seed(42)
 
	p <- ggplot(subDf,aes(x=logFC,y=-log10(pvalue),color=significant)) +
		xlab(bquote(log[2]*"(Fold change)")) + ylab(bquote(-log[10]*"(P-value)")) +
		geom_vline(xintercept=0,linetype='dashed') +
		
		# plot all proteins (green = significant, blue = not significant)
		geom_point(size=1.5,stroke=0.3,
			color=ifelse(subDf$significant,"springgreen3","dodgerblue")) +
		
		# label InWeb interactors
		# (yellow = significant, empty black circle = not significant)
		geom_point(subset(subDf,IsInWeb & significant),
			mapping=aes(x=logFC,y=-log10(pvalue)),size=1.5,stroke=0.3,color="yellow") + 
		geom_point(subset(subDf,IsInWeb),mapping=aes(x=logFC,y=-log10(pvalue)),
			size=1.5,stroke=0.3,color="black",shape=1) +
		
		# label bait (red = signficant, orange = not significant)
		geom_point(subset(subDf,gene==bait),
			mapping=aes(x=logFC,y=-log10(pvalue)),size=2,stroke=0.3,color="red") + 	
		geom_point(subset(subDf,gene==bait),
			mapping=aes(x=logFC,y=-log10(pvalue)),
			size=2,stroke=0.3,color="black",shape=1) +	
		geom_text_repel(subset(subDf,gene==bait),
			mapping=aes(label=gene),arrow=arrow(length=unit(0.015,'npc')),
			box.padding=unit(0.1,"lines"),point.padding=unit(0.15,"lines"),
			color="black",size=2) +

		ggtitle(ipNames[i]) + theme_classic() +
		theme(plot.title=element_text(size=8),
			axis.title=element_text(size=8),axis.text=element_text(size=8))

	print(p)
	
	set.seed(123)
 
	p <- ggplot(subDf,aes(x=logFC,y=-log10(pvalue),color=significant)) +
		xlab(bquote(log[2]*"(Fold change)")) + ylab(bquote(-log[10]*"(P-value)")) +
		geom_vline(xintercept=0,linetype='dashed') +
		
		# plot all proteins (green = significant, blue = not significant)
		geom_point(size=1.5,stroke=0.3,
			color=ifelse(subDf$significant,"springgreen3","dodgerblue")) +
		
		# label InWeb interactors
		# (yellow = significant, empty black circle = not significant)
		geom_point(subset(subDf,IsInWeb & significant),
			mapping=aes(x=logFC,y=-log10(pvalue)),
			size=1.5,stroke=0.3,color="yellow") + 
		geom_point(subset(subDf,IsInWeb),
			mapping=aes(x=logFC,y=-log10(pvalue)),
			size=1.5,stroke=0.3,color="black",shape=1) +
		
		# label bait (red = signficant, orange = not significant)
		geom_point(subset(subDf,gene==bait),
			mapping=aes(x=logFC,y=-log10(pvalue)),size=2,stroke=0.3,color="red") + 	
		geom_point(subset(subDf,gene==bait),
			mapping=aes(x=logFC,y=-log10(pvalue)),
			size=2,stroke=0.3,color="black",shape=1) +	
		geom_text_repel(subset(subDf,gene==bait),
			mapping=aes(label=gene),arrow=arrow(length=unit(0.015,'npc')),
			box.padding=unit(0.1,"lines"),point.padding=unit(0.15,"lines"),
			color="black",size=2) +

		ggtitle(ipNames[i]) + theme_classic() +
		theme(plot.title=element_text(size=8),
			axis.title=element_text(size=8),axis.text=element_text(size=8))

	print(p)
	
}

dev.off()

