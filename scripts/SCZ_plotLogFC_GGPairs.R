##########################################################################################
## logFC comparison for detected proteins across timepoints (for the same bait)
##
## Author: Joshua King To Ching
##########################################################################################
library(readxl)
library(purrr)
library(GGally)
input_file <- "../data/SCZ_GenoppiResults.xlsx" 
output_plots_file <- "../output/SCZ_LogFC_GGPairs.pdf" 

sheets <- c(
	"CACNA1C_wk3_1", "CACNA1C_wk3_2", "CACNA1C_wk4_1", "CACNA1C_wk4_2", "CACNA1C_wk7",
	"HCN1_wk2", "HCN1_wk3", "HCN1_wk4", "HCN1_wk7_1","HCN1_wk7_2",
	"SYNGAP1_wk2", "SYNGAP1_wk3", "SYNGAP1_wk4_1", "SYNGAP1_wk4_2", "SYNGAP1_wk7")

# ----------------------------------------------------------------------------------------
# Read from Genoppi Results containing logFC data
df <- lapply(sheets, function(sheet) {
	read_excel(input_file, sheet=sheet)})
names(df) <- sheets

# Determine protein enrichment based on cutoff threshold
fdr_cutoff <- 0.1
logFC_cutoff <- 0
df <- lapply(df, function(d){
	d$significant <- d$logFC > logFC_cutoff & d$FDR <= fdr_cutoff
	return(d)})

# ----------------------------------------------------------------------------------------
# function to extract logFC for the given gene from the given df
extract_genes_logFC <- function(genes, df) {
	map_dbl(genes, function(geneName) {
		logfc <- dplyr::filter(df,gene==geneName)$logFC
		# only keep genes with mapping to unique entry, 
		# i.e., discard genes with <1 or >1 mapped entries in df
		if (length(logfc)==1) {logfc} else {NA}})}

# Extract logFC for genes detected in CACNA1C IPs and organize into data frame for plotting
CACNA1C_timepoints <- c(
	"CACNA1C_wk3_1", "CACNA1C_wk3_2", "CACNA1C_wk4_1", "CACNA1C_wk4_2", "CACNA1C_wk7")
CACNA1C_genes <- df[CACNA1C_timepoints] %>% map(., ~.$gene) %>% flatten_chr() %>% unique()
CACNA1C_logfc_list <- map(df[CACNA1C_timepoints], ~extract_genes_logFC(CACNA1C_genes, .))
CACNA1C_logfc_df <- map_dfc(CACNA1C_logfc_list, ~.) %>% 
	dplyr::mutate(gene=CACNA1C_genes) %>% dplyr::relocate(gene)

# Extract logFC for genes detected in HCN1 IPs and organize into data frame for plotting
HCN1_timepoints <- c("HCN1_wk2", "HCN1_wk3", "HCN1_wk4", "HCN1_wk7_1", "HCN1_wk7_2")
HCN1_genes <- df[HCN1_timepoints] %>% map(., ~.$gene) %>% flatten_chr() %>% unique()
HCN1_logfc_list <- map(df[HCN1_timepoints], ~extract_genes_logFC(HCN1_genes, .))
HCN1_logfc_df <- map_dfc(HCN1_logfc_list, ~.) %>% dplyr::mutate(gene=HCN1_genes) %>% 
	dplyr::relocate(gene)

# Extract logFC for genes detected in SYNGAP1 IPs and organize into data frame for plotting
SYNGAP1_timepoints <- c(
	"SYNGAP1_wk2", "SYNGAP1_wk3", "SYNGAP1_wk4_1", "SYNGAP1_wk4_2", "SYNGAP1_wk7")
SYNGAP1_genes <- df[SYNGAP1_timepoints] %>% map(., ~.$gene) %>% flatten_chr() %>% unique()
SYNGAP1_logfc_list <- map(df[SYNGAP1_timepoints], ~extract_genes_logFC(SYNGAP1_genes, .))
SYNGAP1_logfc_df <- map_dfc(SYNGAP1_logfc_list, ~.) %>%
	dplyr::mutate(gene=SYNGAP1_genes) %>% dplyr::relocate(gene)

# ----------------------------------------------------------------------------------------
pdf(output_plots_file,height=10/3,width=7/2)
# Plot GGPairs plot comparing logFC between samples in experiments for CACNA1C IPs 
ggpairs(CACNA1C_logfc_df[CACNA1C_timepoints], lower=list(continuous=wrap("points",size=0.05)), 
	upper=list(continuous=wrap('cor', stars=F))) + theme_bw(base_size = 5)
# Plot GGPairs plot comparing logFC between samples in experiments for HCN1 IPs 
ggpairs(HCN1_logfc_df[HCN1_timepoints], lower=list(continuous=wrap("points",size=0.05)),
	upper=list(continuous=wrap('cor', stars=F))) + theme_bw(base_size = 5)
# Plot GGPairs plot comparing logFC between samples in experiments for SYNGAP1 IPs 
ggpairs(SYNGAP1_logfc_df[SYNGAP1_timepoints], lower=list(continuous=wrap("points",size=0.05)),
	upper=list(continuous=wrap('cor', stars=F)))  + theme_bw(base_size = 5)
dev.off()
