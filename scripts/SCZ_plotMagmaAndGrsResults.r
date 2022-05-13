##########################################################################################
## Generate plots to show common variant enrichment results of SCZ PPI networks
## (1) MAGMA results
## (2) GRS results
##
## Author: Yu-Han Hsu
##########################################################################################

rm(list=ls())

library(ggplot2)
library(RColorBrewer)


# ----------------------------------------------------------------------------------------
# MAGMA PLOTS

# networks to plot
networks <- c('all_combined','wk2_combined','wk3_combined','wk4_combined','wk7_combined',
	'CACNA1C_combined','CUL3_wk7','HCN1_combined',
	'RIMS1_combined','SYNGAP1_combined','TCF4_wk1')
names <- c('All combined','Week 2','Week 3','Week 4','Week 7',
	'CACNA1C','CUL3','HCN1','RIMS1','SYNGAP1','TCF4')

# network size (interactor count in "stringent" combined networks)
sizeDf <- read.table('../data/SCZ_MasterInteractorTable_stringent.txt',
	header=T,sep='\t',stringsAsFactors=F)
sizeDf <- subset(sizeDf, ListName %in% networks & IsInteractor)
sizeDf <- data.frame(table(sizeDf$ListName))
colnames(sizeDf) <- c('ListName','IntCount')


# MAGMA results (all traits)
df <- read.table('../data/SCZ_MagmaResults.txt',header=T,sep='\t',stringsAsFactors=F)

# MAGMA meta-analysis results (SCZ + height)
metaDf <- read.table('../data/SCZ_MagmaResults_meta.txt',
	header=T,sep='\t',stringsAsFactors=F)

# set network names
for (i in 1:length(networks)){
	sizeDf$Network[as.character(sizeDf$ListName)==networks[i]] <- names[i]
	df$Network[df$Network==networks[i]] <- names[i]
	metaDf$Network[metaDf$Network==networks[i]] <- names[i]
}

# set populations
metaDf$Population <- "meta-analysis"

# merge df and metaDf, add IntCount to network labels
mergeDf <- rbind(df[,c('Trait','Population','Network','BETA','SE','P')],
	metaDf[,c('Trait','Population','Network','BETA','SE','P')])
mergeDf$IntCount <- sizeDf$IntCount[match(mergeDf$Network,sizeDf$Network)]
mergeDf$Network <- paste(mergeDf$Network,'\n(',mergeDf$IntCount,')',sep='')

# set plotting orders
netOrder <- paste(rev(names),
	'\n(',sizeDf$IntCount[match(rev(names),sizeDf$Network)],')',sep='')
mergeDf$Network <- factor(mergeDf$Network, levels=netOrder)
mergeDf$Population <- factor(mergeDf$Population, levels=c("meta-analysis","EAS","EUR"))
mergeDf$Trait <- factor(mergeDf$Trait, levels=c("SCZ","ADHD","ASD","BIP","MDD","height"))

# set p-value text color/font based on nominal or Bonferroni significance
# Bonferroni: adjust for 11 networks * 2 GWAS datasets (EUR and EAS)
mergeDf$Sig <- ifelse(mergeDf$P < 0.05,"Nominal","None")
mergeDf$Sig[mergeDf$P<0.05/22] <- "Bonferroni"
mergeDf$Sig <- factor(mergeDf$Sig, levels=c("None","Nominal","Bonferroni"))

mergeDf$pFont <- ifelse(mergeDf$Sig=="Nominal" | mergeDf$Sig=="Bonferroni",4,2)
# 4=bold italic, 2=bold


# forrest plot with p-value labels: SCZ + height
pdf("../output/SCZ_MAGMA_SczHeight_ForrestPlot.pdf",height=4,width=5)

ggplot(subset(mergeDf,Trait=="SCZ"|Trait=="height"),
	aes(x=Network,y=BETA,color=Sig,group=Population,shape=Population)) + 
facet_wrap(~Trait) +
geom_hline(yintercept=0, linetype="dashed", color="grey40",size=0.5) +
geom_point(size=1.5,position=position_dodge(width=0.7)) +
geom_errorbar(aes(ymin=BETA-SE,ymax=BETA+SE),
	width=0.3,position=position_dodge(width=0.7)) +
geom_text(aes(x=Network,label=paste("P=",formatC(P,format="e",digits=1),sep=""),
	y=BETA+SE+0.05,fontface=pFont),
	position=position_dodge(width=0.7),size=1.8,show.legend=FALSE,hjust=0) +
scale_colour_manual(name="Significance",values=c("black","sienna1","red"),guide='none') +
scale_shape_manual(name="Population",values=c(15,16,17),guide=guide_legend(reverse=T)) +
xlab('Network') + ylab("Enrichment coefficient (± SE)") +
theme_bw() + 
theme(legend.position="bottom",
	legend.margin=margin(0,0,0,0),legend.box.margin=margin(-5,-5,-5,-5),
	axis.title=element_text(size=9),axis.text=element_text(size=8)) +
ylim(min(df$BETA-df$SE),max(df$BETA+df$SE)+0.6) +
coord_flip()

dev.off()


# forrest plot: all traits
pdf("../output/SCZ_MAGMA_ForrestPlot.pdf",height=4,width=10)

ggplot(mergeDf,aes(x=BETA,y=Network,color=Sig,group=Population,shape=Population)) + 
facet_wrap(~Trait,nrow=1) +
geom_vline(xintercept=0, linetype="dashed", color="grey40",size=0.5) +
geom_point(size=1.5,position=position_dodge(width=0.7)) +
geom_errorbarh(aes(xmin=BETA-SE,xmax=BETA+SE), height=0.3,position=position_dodge(width=0.7)) +
scale_colour_manual(name="Significance",values=c("black","sienna1","red"),guide='none') +
scale_shape_manual(name="Population",values=c(15,16,17),guide=guide_legend(reverse=T)) +
xlab('Enrichment coefficient (± SE)') + ylab("Network") +
theme_bw() + 
theme(legend.position="bottom",
	legend.margin=margin(0,0,0,0),legend.box.margin=margin(-5,-5,-5,-5),
	axis.title=element_text(size=9),axis.text=element_text(size=8))

dev.off()


# heat map: psychiatric traits (minus height), only include meta-analyzed SCZ results
heatDf <- subset(mergeDf,Trait!="height" & 
	!(Trait=='SCZ' & Population=='EUR') & !(Trait=='SCZ' & Population=='EAS'))
heatDf$Trait <- as.character(heatDf$Trait)
heatDf$Trait <- factor(heatDf$Trait, levels=c("SCZ","ADHD","ASD","BIP","MDD"))

# set significance marker text for nominal for Bonferroni significance
# Bonferroni: adjust for 11 networks * 2 GWAS datasets (EUR and EAS)
heatDf$Sig <- ifelse(heatDf$P < 0.05,"*","")
heatDf$Sig[heatDf$P<0.05/22] <- "**"

heatDf$SigText <- ifelse(heatDf$P < 0.05,
	paste(formatC(heatDf$BETA,format="f",digits=2),heatDf$Sig,sep=''),"")


# heat map figure panel with colors base on -log10(P)
pdf("../output/SCZ_MAGMA_HeatMap.pdf",height=4,width=3.5)

ggplot(heatDf,aes(x=Trait,y=Network,fill=-log10(P))) +
geom_tile() + geom_text(aes(label=SigText),size=3) +
scale_fill_gradient(name=expression(paste(-log[10],"(P-value)",sep="")),
	low="white", high="red") +
scale_x_discrete("Psychiatric disorder", labels=parse(text=levels(heatDf$Trait))) +
theme_bw() + 
theme(legend.position="bottom",legend.key.size=unit(0.7,"line"),
	legend.title=element_text(size=9),legend.text=element_text(size=8),
	legend.margin=margin(0,0,0,0),legend.box.margin=margin(-5,-5,-5,-5),
	axis.title=element_text(size=9),axis.text=element_text(size=8))

dev.off()


# ----------------------------------------------------------------------------------------
# GRS PLOTS

# networks to plot
networks_short <- c('all_combined','wk2_combined','wk3_combined','wk4_combined','wk7_combined',
	'CACNA1C_combined','CUL3_wk7','HCN1_combined',
	'RIMS1_combined','SYNGAP1_combined','TCF4_wk1')
names <- c('All combined','Week 2','Week 3','Week 4','Week 7',
	'CACNA1C','CUL3','HCN1','RIMS1','SYNGAP1','TCF4')
networks <- paste(networks_short,"_int",sep="")


# GRS results
df <- read.table("../data/SCZ_GrsResults.txt",header=T,sep="\t",stringsAsFactors=F)
df <- subset(df, complete.cases(df) & VARIABLE %in% networks)
df$P <- df$pval

# GRS meta-analysis results
metaDf <- read.table('../data/SCZ_GrsResults_meta.txt',
	header=T,sep='\t',stringsAsFactors=F)
metaDf <- subset(metaDf, Network %in% networks_short)

# set network names
for (i in 1:length(networks)){
	df$Network[df$VARIABLE==networks[i]] <- names[i]
	metaDf$Network[metaDf$Network==networks_short[i]] <- names[i]
}

# set populations
df$POP <- toupper(df$POP)
metaDf$POP <- "meta-analysis"

# set trait names
df$PHE[df$PHE=="scz"] <- "SCZ"
metaDf$PHE <- metaDf$Trait

# merge df and metaDf, add IntCount to network labels
mergeDf <- rbind(df[,c('PHE','POP','Network','BETA','SE','P')],
	metaDf[,c('PHE','POP','Network','BETA','SE','P')])
mergeDf$IntCount <- sizeDf$IntCount[match(mergeDf$Network,sizeDf$Network)]
mergeDf$Network <- paste(mergeDf$Network,'\n(',mergeDf$IntCount,')',sep='')

# set plotting orders
mergeDf$Network <- factor(mergeDf$Network, levels=netOrder)
mergeDf$POP <- factor(mergeDf$POP, levels=c("meta-analysis","EAS","EUR"))
mergeDf$PHE <- factor(mergeDf$PHE, levels=c("SCZ","height"))

# set p-value text color/font based on nominal or Bonferroni significance
# Bonferroni: adjust for 11 networks * 2 GWAS datasets (EUR and EAS)
mergeDf$Sig <- ifelse(mergeDf$P < 0.05,"Nominal","None")
mergeDf$Sig[mergeDf$P<0.05/22] <- "Bonferroni"
mergeDf$Sig <- factor(mergeDf$Sig, levels=c("None","Nominal","Bonferroni"))


# forrest plot: SCZ + height
pdf("../output/SCZ_GRS_ForrestPlot.pdf",height=4,width=4)

ggplot(mergeDf,aes(x=BETA,y=Network,color=Sig,group=POP,shape=POP)) + 
facet_wrap(~PHE,nrow=1) +
geom_vline(xintercept=0, linetype="dashed", color="grey40",size=0.5) +
geom_point(size=1.5,position=position_dodge(width=0.7)) +
geom_errorbarh(aes(xmin=BETA-SE,xmax=BETA+SE),
	height=0.3,position=position_dodge(width=0.7)) +
scale_colour_manual(name="Significance",values=c("black","sienna1","red"),guide='none') +
scale_shape_manual(name="Population",values=c(15,16,17),guide=guide_legend(reverse=T)) +
xlab('Enrichment coefficient (± SE)') + ylab("Network") +
theme_bw() + 
theme(legend.position="bottom",
	legend.margin=margin(0,0,0,-10),legend.box.margin=margin(-5,-5,-5,-5),
	axis.title=element_text(size=9),axis.text=element_text(size=8))

dev.off()
