##########################################################################################
## Significant interactors overlap enrichment Heatmap for each bait IPs
##
## Author: Joshua King To Ching 
##########################################################################################

library(circlize)
library(ComplexHeatmap)
library(dplyr)
library(purrr)
library(stringr)
library(readxl)
library(GGally)
input_file <- "../data/SCZ_GenoppiResults.xlsx" 
output_plots_file <- "../output/SCZ_InteractorsOverlap_Heatmap.pdf" 

# ----------------------------------------------------------------------------------------
sheets <- list(
	CACNA1C=c("CACNA1C_wk3_1", "CACNA1C_wk3_2", "CACNA1C_wk4_1", "CACNA1C_wk4_2",
		"CACNA1C_wk7"),
	HCN1=c("HCN1_wk2", "HCN1_wk3", "HCN1_wk4", "HCN1_wk7_1","HCN1_wk7_2"),
	SYNGAP1=c("SYNGAP1_wk2", "SYNGAP1_wk3", "SYNGAP1_wk4_1", "SYNGAP1_wk4_2", 
		"SYNGAP1_wk7"))

# Read input
df <- lapply(flatten(sheets), function(sheet) {
	read_excel(input_file, sheet=sheet)})
names(df) <- flatten(sheets)

# Determine protein enrichment based on cutoff threshold
fdr_cutoff <- 0.1
logFC_cutoff <- 0
df <- lapply(df, function(d){
	d$significant <- d$logFC > logFC_cutoff & d$FDR <= fdr_cutoff
	return(d)})

# ----------------------------------------------------------------------------------------
# initialize table for storing overlap enrichment statistics
statsTable <- map(sheets, ~expand.grid(x=.,y=.))

# calculate pVal of hypergeometric test where x is the sample
getOverlapStat <- function(x_name,y_name, bait) {
	x_df <- df[[x_name]]
	y_df <- df[[y_name]]
	popGenes <- intersect(x_df$gene, y_df$gene)
	popGenes <- popGenes[popGenes != bait]
	xEnriched <- x_df %>% filter(significant, gene%in%popGenes) %>% pull(gene) %>% 
		unique()
	yEnriched <- y_df %>% filter(significant, gene%in%popGenes) %>% pull(gene) %>% 
		unique()
	notXEnriched <- setdiff(popGenes, xEnriched)
	xyEnriched <- intersect(xEnriched, yEnriched)
	phyperIn <- list(q=xyEnriched, m=xEnriched, n=notXEnriched, k=yEnriched) %>% 
		map(length)
	pVal <- stats::phyper(
		phyperIn$q-1, phyperIn$m, phyperIn$n, phyperIn$k, lower.tail=F)
	overlapPercent <- phyperIn$q/phyperIn$m
	# add list of Overlapping genes to table (comma separated)
	overlapGenes <- paste(xyEnriched, collapse=",")
	return(data.frame(q=phyperIn$q, m=phyperIn$m, n=phyperIn$n, k = phyperIn$k, 
		pop=length(popGenes), pVal=pVal, overlapPer=overlapPercent, 
		overlapGenes=overlapGenes))}

# For each IP, calculate the percentage and enrichment P-values(using one-tailed 
# hypergeometric test) for overlap of interactors between two timepoint 
statsTable <- lapply(names(statsTable), function(bait) {
	overlapStat <- statsTable[[bait]] %>%
		rowwise() %>%
		transmute(getOverlapStat(as.character(x), as.character(y), bait))
	return(bind_cols(statsTable[[bait]], overlapStat))})
statsTableOut <- bind_rows(statsTable)

# Labeling the overlap enrichment statistics output
names(statsTableOut)[names(statsTableOut) == "x"] <- "Timepoint"
names(statsTableOut)[names(statsTableOut) == "y"] <- "ReferenceTimepoint"
names(statsTableOut)[names(statsTableOut) == "q"] <- "Overlap"
names(statsTableOut)[names(statsTableOut) == "m"] <- "Interactor"
names(statsTableOut)[names(statsTableOut) == "n"] <- "NonInteractor"
names(statsTableOut)[names(statsTableOut) == "k"] <- "ReferenceInteractor"
names(statsTableOut)[names(statsTableOut) == "pop"] <- "Population"
names(statsTableOut)[names(statsTableOut) == "pVal"] <- "P"
names(statsTableOut)[names(statsTableOut) == "overlapPer"] <- "OverlapPercentage"
names(statsTableOut)[names(statsTableOut) == "overlapGenes"] <- "OverlapGenes"
statsTableOut$OverlapGenes <- as.character(statsTableOut$OverlapGenes)

# add cell type information to output stat table
iN_WA01_H1 <- c(
	"CACNA1C_wk3_1", "CACNA1C_wk3_2", "CACNA1C_wk4_1", "CACNA1C_wk7", "HCN1_wk7_1")
map_timepoint_to_celltype <- function(tp) {
	if(tp %in% iN_WA01_H1) {"WA01(H1)"} else {"iPS3"}}
statsTableOut <- statsTableOut %>%
	rowwise() %>%
	mutate(CellType=map_timepoint_to_celltype(Timepoint), 
		ReferenceCellType=map_timepoint_to_celltype(ReferenceTimepoint)) %>%
	relocate(OverlapGenes, .after=last_col())  
names(statsTable) <- names(sheets)
statsTable <- map(statsTable, ~select(.,x,y,pVal,overlapPer))

# uncomment the line below to get the statistics of the overlap
# write.table(
# 	statsTableOut, "../output/SCZ_OverlapEnrichmentStats.csv", row.names=F, sep="\t")

# ----------------------------------------------------------------------------------------
# extract the p-value from the overlap enrichment statistics table and organize into list 
# for plotting
pValMatList <- statsTable %>%
	map(~reshape(idvar="x", direction="wide", timevar="y", drop="overlapPer",data=., 
	new.row.names=unique(.$x))) %>%
	map(~data.matrix(.)) %>% 
	map(~.[,-1])    
# extract the overlap percentage from the overlap enrichment statistics table and organize 
# into list for plotting
overlapPerMatList <- statsTable %>%
	map(~reshape(idvar="x", direction="wide", timevar="y", drop="pVal",data=., 
	new.row.names=unique(.$x))) %>%
	map(~data.matrix(.)) %>% 
	map(~.[,-1])    
# ----------------------------------------------------------------------------------------
pdf(file=output_plots_file)

# define plotting variables for customizing plots aesthetic and labelling
col_fun = colorRamp2(c(0, 1), c("white", "red"))
lgd = Legend(col_fun = col_fun, title = "Overlap Percentage", direction = "horizontal", 
	legend_width = unit(4, "cm"),at=c(0,0.5,1), labels = c("0%", "50%", "100%"))
draw(lgd)

# For each IP (CACNA1C, HCN1, and SYNGAP1), plot heatmap showing the overlap percentage of 
# significant interactors between different timepoints and with hierachical clustering of 
# the overlap percentages.
lapply(names(statsTable), function(bait) {
	mat <- overlapPerMatList[[bait]]
	column_labels <- 
		structure(map(colnames(mat), ~str_extract(.,"wk\\w+")),names=colnames(mat))
	row_labels <- 
		structure(map(rownames(mat), ~str_extract(.,"wk\\w+")),names=rownames(mat))
	cell_fun = function(j,i,x,y,w,h,col) {
		grid.text(paste(round(mat[i,j]*100,0), "%\npvalue:\n", 
		format(pValMatList[[bait]][i,j], scientific=T, digits=2),sep=''), 
		x, y, gp=gpar(fontsize=10))}
	bottomAnnot = HeatmapAnnotation(
		text = anno_text(column_labels[colnames(mat)], location = 0.5, rot = 0, 
		just = "center", gp = gpar(fontsize = 12)))
	Heatmap(mat,col=col_fun,cell_fun=cell_fun,name="overlap\npercentage", 
		column_title=bait, column_labels = column_labels[colnames(mat)], 
		row_labels=row_labels[rownames(mat)], show_heatmap_legend=F, 
		bottom_annotation=bottomAnnot, show_column_names=F) 
	})
dev.off()
