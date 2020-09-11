#05_adjust_pvals.R

require(tidyverse)

load("./scratch_data/cont_analyses_ests_pre_mtc.RData")
load("./scratch_data/dich_analyses_ests_pre_mtc.RData")

ests_tab <- ests1_tab %>% 
  bind_rows(ests3_tab) %>% 
  mutate(bonf_pval = p.adjust(pvalue, method="bonferroni")) %>% 
  filter(ifelse(model=="sex_diff", sex %in% c("Male", "Female"),sex %in% c( "Female")  )) %>% 
  mutate(Sex = ifelse(model=="sex_diff", sex , "Male/Female"  ),
         fdr_pval = p.adjust(pvalue, method="fdr")) %>% 
  mutate(pheno = factor(str_sub(score, 0,3), levels = c("adh","asd","scz"),labels = c("ADHD","ASD","Schizophrenia")),
         var = factor(var, levels = c("age_1st_wlk","bin_mot_18m",
                                      "age_1st_wds", "age_1st_sns", "diff_1st_wds_sns","bin_lang_3yr",  "concern_dev_exprssd_others"), 
                      labels =c("Age at first walking","Motor delays at 18 months", 
                                "Age at first words", "Age at first sentences", "Rate (age first sentences - age first words)","Language delays at 3 years","Concern regarding development expressed by others at 3 years") )) %>% 
  select(pheno, var, Sex, std.beta = est.std, raw_pval= pvalue, fdr_pval, bonf_pval )


write.table(ests_tab,file= "./output/estimates_corrected_pvals.txt", row.names = F, quote = F, sep="\t")
