# Ruize Liu
# GRS analysis
# Inverse-variance weighted fixed-effect meta-analysis on interaction term across cohorts

library(metafor)
args = commandArgs(trailingOnly = TRUE)

# in_list_addr_l = c('IP_and_InWeb_list_test_p.txt', 'NC_combine_test_p.txt', 'bine_1percent_list_test_p.txt', 'bine_5percent_list_test_p.txt', 'bine_10percent_list_test_p.txt')
in_list_addr_l=args

for (in_list_addr in in_list_addr_l){

	file_list=read.table(in_list_addr, sep='\t')

	dat = c()
	cohort_list = c()
	# group_list = c()
	group_name = sub('_test_p.txt', '', in_list_addr)
	# group_list  = c(group_list, group_name)
	n_files = nrow(file_list)
	for (i in 1:n_files){
		in_tab_addr = as.character(file_list[i,1])
		cohort_name = unlist(strsplit(in_tab_addr, split ='/'))[2]
		
		each_dat = read.table(in_tab_addr, sep = '\t', header=TRUE)
		each_dat$cohort = cohort_name
		each_dat$group  = group_name
		# group_list  = c(group_list, group_name)
		dat = rbind(dat, each_dat)
		cohort_list = c(cohort_list, cohort_name)
	}

	print(paste0('loaded data: ', group_name))
	# for (each_group in group_list){
	each_group = group_name 
	dir.create(paste0('meta_analysis/', each_group, '_meta/forestplot'), recursive=TRUE)
	# group_dat = dat[dat$group == each_group, ]
	group_dat = dat
	group_res = c()
	for (each_ip in unique(group_dat$int)) {
		each_ip_dat = group_dat[group_dat$int == each_ip, ]

		res = rma(dat=each_ip_dat, yi=interaction_term_beta.norm_prs., sei=interaction_term_se.norm_prs.,
			method='FE', weights=1/interaction_term_se.norm_prs.)
		out_res = c(as.character(each_ip), as.character(each_group),res$beta, res$se, res$zval, res$pval, res$ci.lb, res$ci.ub, res$QE, res$QEp, res$H2)

		pdf(paste0('meta_analysis/', each_group, '_meta/forestplot/', each_ip, '_forestplot.pdf'))
		forest(res, slab=each_ip_dat$cohort, cex=0.7, main=paste0(as.character(each_ip),'\nQE:',res$QE, '; QEp:',res$QEp))
		dev.off()
		group_res = rbind(group_res, out_res)
	}
	# }
	colnames(group_res) = c('int', 'group', 'beta', 'se', 'zval', 'pval', 'ci.lb', 'ci.ub', 'QE', 'QEp', 'H2')
	write.table(group_res, paste0('meta_analysis/', each_group, '_meta/', each_group, '_summary.txt'), sep='\t', row.names=FALSE, quote=FALSE)
	print(paste0('done: ', group_name))
}

