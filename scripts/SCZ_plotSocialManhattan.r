##########################################################################################
## Generate SCZ social Manhattan plot using PGC3 GWAS data
##
## Author: Yu-Han Hsu
##########################################################################################

rm(list=ls())

library(ggplot2)
library(dplyr)
library(ggrepel)
library(ggnewscale)


# ----------------------------------------------------------------------------------------
# PGC3 FINEMAP/SMR genes (medRxiv manuscript Supp Table 20)
pgcDf <- read.table('../data/PGC3_SCZ_PrioritizedGenes.txt',
	header=T,sep='\t',stringsAsFactors=F)
pgcDf$Symbol.ID[grepl('PCDHA',pgcDf$Symbol.ID)] <- 'PCDHA@'
pgcGenes <- unique(pgcDf$Symbol.ID)


# PGC3 SCZ GWAS genes in LD loci (medRxiv manuscript Supp Table 3)
geneDf <- read.table('../data/PGC3_SCZ_combined_270loci_GenePosPvalue.txt',
	header=T,sep='\t',stringsAsFactors=F)
# remove MHC region genes (except SYNGAP1)
geneDf <- subset(geneDf,SNP!='rs13195636' | Gene=='SYNGAP1')

# collapse PCDHA@ genes
geneDf[grepl('PCDHA',geneDf$Gene),] <- geneDf[geneDf$Gene=='PCDHA1',]
geneDf <- unique(geneDf)
geneDf$Gene[geneDf$Gene=='PCDHA1'] <- 'PCDHA@'

# set minimum p-value (only affects SYNGAP1)
geneDf[geneDf$Pvalue < 1e-25,"Pvalue"] <- 1e-25 
geneDf <- geneDf[order(geneDf$Chr,geneDf$Start),]
dim(geneDf)

# set chr 1-22 coordinates on x-axis
nChr <- length(unique(geneDf$Chr))
geneDf$Coord <- NA
s <- 0
nbp <- c(249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,
	141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,
	81195210,78077248,59128983,63025520,48129895,51304566) # GRCh37 chr 1-22 lengths
chrTicks <- NULL
for (i in 1:22) {
	geneDf[geneDf$Chr==i,"Coord"] <- geneDf[geneDf$Chr==i,"Start"] + s
	chrTicks <- c(chrTicks,floor(nbp[i]/2) + s)
	s <- s + nbp[i]
}


# ----------------------------------------------------------------------------------------
# read in SCZ interactors
intDf <- read.table('../data/SCZ_MasterInteractorTable_full.txt',
	header=T,sep='\t',stringsAsFactors=F)

# pull interactors from individual IPs (to keep the HCN1-SYNGAP1 interaction,
# which was excluded in combined network since SYNGAP1 is a bait)
intDf$IsCombined <- sapply(strsplit(intDf$ListName,"_"),"[[",2)=='combined'
intDf <- subset(intDf, !IsCombined & IsInteractor)
intDf$Bait <- sapply(strsplit(intDf$ListName,"_"),"[[",1)
length(unique(intDf$GeneName)) # 1239 genes (1238 interactors + SYNGAP1)


# save all PCDHA interactors as 'PCDHA@'
intDf[grep('PCDHA',intDf$Gene),'GeneName'] <- 'PCDHA@' 
intDf <- unique(intDf[,c('Bait','GeneName')])
length(unique(intDf$GeneName)) # 1229 genes after merging PCDHA@ genes

intDf <- unique(intDf[,c('Bait','GeneName')])
names(intDf) <- c('Bait','Gene')
dim(intDf)


# ----------------------------------------------------------------------------------------
# create data frame for genes with GWAS stats + are baits or interactors in SCZ IPs
plotDf <- merge(geneDf,intDf,by.x="Gene",by.y="Gene")

plotDf$Bait_Pvalue <- NA
plotDf$Bait_Coord <- NA
baitList <- unique(plotDf$Bait)
for (bait in baitList) {
	# add bait coordinates (for plotting social links between genes and baits)
	plotDf[plotDf$Bait==bait,"Bait_Pvalue"] <- geneDf$Pvalue[geneDf$Gene==bait]
	plotDf[plotDf$Bait==bait,"Bait_Coord"] <- geneDf$Coord[geneDf$Gene==bait]
}
for (bait in baitList) {
	baitDf <- data.frame(geneDf[geneDf$Gene==bait,],Bait=NA,Bait_Pvalue=NA,Bait_Coord=NA)
	plotDf <- rbind(plotDf,baitDf,row.names=NULL)
}
dim(plotDf)


# significance group (for node color)
plotDf$Significance <- 'Locus genes'
plotDf$Significance[plotDf$Gene %in% pgcGenes] <- 'FINEMAP/SMR genes'
plotDf[plotDf$Gene %in% baitList,'Significance'] <- 'Index genes'
plotDf$Significance <- factor(plotDf$Significance,
	levels=c('Index genes','Locus genes','FINEMAP/SMR genes'))


# other plotting parameters
ymin <- abs(ceiling(log10(max(plotDf$Pvalue)))) - 0.5
ymax <- abs(floor(log10(min(plotDf$Pvalue)))) + 0.5


# unique genes that show up in plotDf (which contains unique PAIRS of genes)
uniqDf <- unique(plotDf[,c("Gene","Coord","Pvalue","Significance")])
# 6 index genes + 104 locus genes


# ----------------------------------------------------------------------------------------
# generate social Manhattan plot
# label index and locus genes
# highlight IP-WB validated interactions
# highlight overlap with FINEMAP/SMR genes


# edge color based on IP-WB validation data
wbDf <- read.table('../data/SCZ_ValidationSummary.txt',
	header=T,sep='\t',stringsAsFactors=F,row.names=1)
wbDf <- subset(wbDf,Validated)
wbDf$Interactor[wbDf$Interactor=='PCDHA2'] <- 'PCDHA@'

plotDf$Validated <- FALSE
for (i in 1:nrow(wbDf)) {
	plotDf$Validated[which(plotDf$Gene==wbDf$Interactor[i] & 
		plotDf$Bait==wbDf$Bait[i])] <- TRUE
}
plotDf <- plotDf[order(plotDf$Validated),]


# node size based on # of interactions
sizeDf <- data.frame(table(plotDf$Gene))
maxNumBaits <- max(sizeDf$Freq)
uniqDf$NodeSize <- sizeDf$Freq[match(uniqDf$Gene,sizeDf$Var1)]
uniqDf$NodeSize[uniqDf$Gene %in% baitList] <- maxNumBaits


set.seed(1234) # fix ggrepel label coordinates

pdf("../output/SCZ_SocialManhattanPlot.pdf",width=7.5,height=3.5)

ggplot(uniqDf,mapping=aes(x=Coord,y=-log10(Pvalue),size=NodeSize,
	color=as.factor(Significance))) + 

# social links
geom_segment(plotDf,mapping=aes(x=Coord,y=-log10(Pvalue),
	xend=Bait_Coord,yend=-log10(Bait_Pvalue)),size=0.25,alpha=0.7,
	color=ifelse(plotDf$Validated,'blue','grey')) +

# data point for each gene
geom_point(show.legend=F) +

# gene name label
geom_text_repel(mapping=aes(label=Gene,color=as.factor(Significance),size=NodeSize),
	fontface='bold',segment.size=0.1,segment.alpha=0.5,max.overlaps=100,show.legend=F) +

scale_size(range=c(1.5,2.5)) +
scale_color_manual(labels=c("Index genes","Locus genes","FINEMAP/SMR genes"),
	values=c("red","darkorchid4","magenta")) +

scale_x_continuous(limits=c(0,max(geneDf$Coord)),label=as.character(1:22),
	breaks=chrTicks) +
scale_y_continuous(expand=c(0,0),limits=c(5.5,25.5),breaks=seq(5.5,25.5,2.5)) +

xlab("Chromosomal position") + ylab(expression(paste(-log[10],"(P-value)",sep=""))) +
theme_classic() +
theme(axis.title=element_text(size=9),axis.text=element_text(size=8))

dev.off()

