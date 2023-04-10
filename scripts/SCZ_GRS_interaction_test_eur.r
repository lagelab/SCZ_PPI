# Ruize Liu
# GRS analysis 
# For EUR, fit the linear model GRS ~ Int_Non + Case_Control + Int_Non * Case_Control + PCs and test interaction term: Int_Non * Case_Control.


args = commandArgs(T)

prs_dir = args[1]    # input prs dir where prs file located
base_name = args[2]  # prefix of output files: 
print(prs_dir)
print(base_name)

setwd(prs_dir)


# library(ggplot2)
file_name = list.files()
is_profile = grepl('.profile$', file_name)
file_name_pro = file_name[is_profile]

ip_name = sub('.profile$', '', file_name_pro)

n = length(file_name_pro)
for (i in 1:n) {
  each_dat = read.table(file_name_pro[i], header = T)
  if (i==1){
    res_dat = each_dat[,c(1,2,3,6)]
  } else {
    res_dat = cbind(res_dat, each_dat[,6])
  }
  colnames(res_dat)[ncol(res_dat)] = ip_name[i]
}

pca_dat = read.table('~/PGC2/PGC2_data/prune.bfile.cobg.PGC_SCZ49.sh2.menv.mds_cov', header = T)


dir.create(paste0(base_name,'_new_fig'), recursive = TRUE)
setwd(paste0(base_name,'_new_fig'))

# base_name = 'NC_combine'
write.table(res_dat, paste0(base_name,'_prs_tab.txt'), sep='\t', row.names = FALSE, quote = FALSE)

######
# dat = read.table(paste0(base_name,'prs_tab.txt'), sep='\t', header=T)
dat = res_dat


cols=c('red','green')
k = floor((ncol(dat) - 3) / 2)

dat = merge(res_dat, pca_dat, by = c('FID', 'IID'))


dat$PHENO[dat$PHENO==1] = 'control'
dat$PHENO[dat$PHENO==2] = 'case'
phe = c(dat$PHENO, dat$PHENO)

n = nrow(dat)
n_case = sum(dat$PHENO=='case')
n_control = sum(dat$PHENO=='control')


type = c(rep('int',n), rep('non',n))

C1 = c(dat$C1, dat$C1)
C2 = c(dat$C2, dat$C2)
C3 = c(dat$C3, dat$C3)
C4 = c(dat$C4, dat$C4)
C5 = c(dat$C5, dat$C5)
C6 = c(dat$C6, dat$C6)
C7 = c(dat$C7, dat$C7)
C9 = c(dat$C9, dat$C9)
C15 = c(dat$C15, dat$C15)
C18 = c(dat$C18, dat$C18)


res = c()
for (i in c(1:k)){
  j=2*i + 2
  # print(summary(dat[,j]))
  # print(summary(dat[,j+1]))
  file_name1 = colnames(dat)[j]
  file_name2 = colnames(dat)[j+1]
  prs = c(dat[,j], dat[,j+1])
  m = mean(dat[dat$PHENO=='control', j+1])
  std = sd(dat[dat$PHENO=='control', j+1])
  prs = (prs-m) / std
  prs_tab = data.frame(phe, type, prs, C1, C2, C3, C4, C5, C6, C7, C9, C15, C18)

  fit = lm(prs ~ phe + type + phe * type + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C9 + C15 + C18, data = prs_tab)
  # ans = anova(fit)
  # p_intaction = ans$`Pr(>F)`[3]
  # res =rbind(res, c(file_name1, file_name2, int_ttest_p, non_ttest_p, int_ttest_dff, non_ttest_dff, p, p_intaction) )
  s = summary(fit)
  # interaction is the last term

  n_term = nrow(s$coefficients)
  beta_interaction = s$coefficients[n_term,1]
  beta_se = s$coefficients[n_term,2]
  p_intaction = s$coefficients[n_term,4]
  res = rbind(res, c(file_name1, file_name2, beta_interaction, beta_se, p_intaction, n_case, n_control, n))
}
colnames(res) =c('int','non', 'interaction_term_beta(norm_prs)','interaction_term_se(norm_prs)','interaction_term_p', 'n_case', 'n_control', 'sample_size')
write.table(res, paste0(base_name, '_dff_interaction_term.txt'), sep='\t',quote = FALSE,row.names = FALSE)
