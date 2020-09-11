
load( file = './data/prepped_00_vars.RData')
#install.packages("gtools", repos = "file://tsd-evs/shared/R/cran")

library(gtools)
library(tidyverse)
library(lavaan)
require('semPlot')
library(semTools)



##
#Multiple testing correction - currently for 3 PRS * 5 main vars - applied later
pval <- 0.05

###############################################################3

#Restrict to genotyped individuals for this anaylsis as complete cases needed

alldata <- alldata %>%
  filter(sex %in% c("Male","Female")) %>% 
  drop_na(IID) %>%
  droplevels()

# Secondary anaylses: Logistic regression for binary outcomes 

xs <- c("scz_X0.001_resid", "scz_X0.01_resid", "scz_X0.05_resid", "scz_X0.1_resid", 
        "scz_X0.2_resid", "scz_X0.5_resid", "scz_X1_resid",  "asd_X0.001_resid", 
        "asd_X0.01_resid", "asd_X0.05_resid", "asd_X0.1_resid", "asd_X0.2_resid", 
        "asd_X0.5_resid", "asd_X1_resid",   
        "adh_X0.001_resid", "adh_X0.01_resid", "adh_X0.05_resid", 
        "adh_X0.1_resid", "adh_X0.2_resid", "adh_X0.5_resid", "adh_X1_resid")
ys <- c( "bin_mot_18m","bin_lang_3yr", "concern_dev_exprssd_others")

ests2 <- data.frame()

params <- expand.grid(xs,ys, stringsAsFactors = FALSE)

fitcomps <- data.frame()
prprobs <- data.frame()
rsqs <- data.frame() 

model<- 'y ~ x'

for(i in 1:dim(params)[[1]]){ 
  
  alldata[,"x"] <- alldata[,params[i,1]]
  alldata[,"y"] <- alldata[,params[i,2]]
  
  fit <- sem(model, data=alldata, ordered= c("y"),group = "sex", link="probit")
  fit1 <- sem(model, data=alldata, ordered= c("y"),group = "sex", link="probit", group.equal = c("regressions"))

  
  fitcomps <- rbind(fitcomps,lavTestLRT(fit,fit1))
  
  tmprsq_vars <- as.data.frame(inspect(fit, "r2"))%>% gather(group, estimate) %>%  mutate(model= "sex_diff") %>% 
    rbind(as.data.frame(inspect(fit1, "r2")) %>% gather(group, estimate) %>%  mutate(model= "fixed_beta")) 

  rsqs <- rbind(rsqs, tmprsq_vars)
  
  
  res <- standardizedSolution(fit, ci = TRUE, level = 0.95)
  res1 <- standardizedSolution(fit1, ci = TRUE, level = 0.95)
  #  res2 <- standardizedSolution(fit2, ci = TRUE, level = 0.95)
  
  ests2 <- ests2 %>%
    bind_rows(res %>%
                mutate(model = "sex_diff")) %>%
    bind_rows(res1 %>%
                mutate(model = "fixed_beta")) 
  
  # Compute  predicted probabilities for each model at pre-defined PRS values
  
  prprobs <-  rbind(prprobs,
                    data.frame(prsval = seq(-3,3,0.1),
                               est = pnorm( res[2,"est.std"] + res[1,"est.std"] * seq(-3,3,0.1) ),
                               uci = pnorm( res[2,"est.std"] + res[1,"ci.upper"] * seq(-3,3,0.1) ),
                               lci = pnorm( res[2,"est.std"] + res[1,"ci.lower"] * seq(-3,3,0.1) ),
                               group = 1,
                               model = c("sex_diff")),
                    data.frame(prsval = seq(-3,3,0.1),
                               est = pnorm( res[9,"est.std"] + res[8,"est.std"] * seq(-3,3,0.1) ),
                               uci = pnorm( res[9,"est.std"] + res[8,"ci.upper"] * seq(-3,3,0.1) ),
                               lci = pnorm( res[9,"est.std"] + res[8,"ci.lower"] * seq(-3,3,0.1) ),
                               group = 2,
                               model = c("sex_diff")),
                    data.frame(prsval = seq(-3,3,0.1),
                               est = pnorm( res1[2,"est.std"] + res1[1,"est.std"] * seq(-3,3,0.1) ),
                               uci = pnorm( res1[2,"est.std"] + res1[1,"ci.upper"] * seq(-3,3,0.1) ),
                               lci = pnorm( res1[2,"est.std"] + res1[1,"ci.lower"] * seq(-3,3,0.1) ),
                               group = 1,
                               model = c("fixed_beta")),
                    data.frame(prsval = seq(-3,3,0.1),
                               est = pnorm( res1[9,"est.std"] + res1[8,"est.std"] * seq(-3,3,0.1) ),
                               uci = pnorm( res1[9,"est.std"] + res1[8,"ci.upper"] * seq(-3,3,0.1) ),
                               lci = pnorm( res1[9,"est.std"] + res1[8,"ci.lower"] * seq(-3,3,0.1) ),
                               group = 2,
                               model = c("fixed_beta")))
}

ests3 <- ests2 %>%
  filter(op %in% c('~')) %>% 
  mutate(score = rep(params$Var1, each = 4),
         var = rep(params$Var2, each = 4)) %>%
  mutate(sex = recode(group, `1`="Male",`2`= "Female")) %>%
  select(score,var,everything())


# Select best-fitting models

fitcomps1 <- fitcomps %>%
  mutate(score = rep(params$Var1, each = 2),
         var = rep(params$Var2, each = 2),
         pLRT = round(`Pr(>Chisq)`,2),
         model = rep(c("sex_diff","fixed_beta"), dim(fitcomps)[1]/2 ) ) %>%
  select(-AIC,-BIC, -`Pr(>Chisq)`)

fitcomp_best <- fitcomps1  %>%
  group_by(var, score) %>%
  summarise(min = min(pLRT, na.rm = T)) %>%
  mutate(model = ifelse(min > pval, "fixed_beta", "sex_diff")) %>%
  select(score,var,model) %>%
  distinct()


fitcomps_out <- fitcomps1 %>% 
  left_join(fitcomp_best %>% 
              rename(best=model)) %>% 
  mutate(selected=ifelse(model==best,"Y",NA))


rsqs_sums <- rsqs %>% 
  mutate(score = rep(params$Var1, each = 4),
         var = rep(params$Var2, each = 4)) %>% 
  group_by(var,score,model) %>% 
  summarise(sum_rsq = sum(estimate)) 

fitcomps_rsq <- fitcomp_best %>% 
  left_join(rsqs_sums) %>% 
  mutate(pheno = str_sub(score, end=3)) %>% 
  filter(!score == "asd_X0.001_resid")  ##NB here we are manually excluding the possibility of the lowest ASD threshold being selected
                                        ## for presentation, because it differs so much from the pattern at the other thresholds for the 
                                        ## lang delay variable, so is likely untrustworthy - it will be retained for the all_thresh plots
fitcomps_rsq <- fitcomps_rsq %>%
  group_by(pheno,var) %>% 
  filter(sum_rsq == max(sum_rsq)) %>% 
  ungroup()

rsqs_bysex <- rsqs %>% 
  mutate(score = rep(params$Var1, each = 4),
         var = rep(params$Var2, each = 4)) %>% 
  right_join(fitcomps_rsq)


save(fitcomps_out,fitcomp_best,fitcomps_rsq,rsqs_bysex, file = "./scratch_data/dichotomous_analyses_fitcomps.RData")


# Table with ests from best fitting model at best threshold, and corrected p.values
# corrections = 1) Bonferroni for all possible params of interest (i.e., male + female for all comparisons)
# 2) Benjamini-Hochberg FDR

ests3_tab <- ests3 %>%
  right_join(fitcomps_rsq %>%
               rename(best=model) %>%
               select(score, var, best)) %>%
  group_by(score, var) %>% 
  filter(model == best) %>% 
  ungroup()

# Save out in this form for combination with results from 02_analyses_dichotomous
# MTC performed in 05_adjust_pvals

save(ests3_tab, file="./scratch_data/dich_analyses_ests_pre_mtc.RData")

# Filter ests on best fitting models

ests3_plt <- ests3 %>%
  left_join(fitcomp_best %>%
              rename(best=model) %>%
              select(score, var, best))%>% 
  rename(Beta = est.std,
         uci = ci.upper,
         lci = ci.lower) %>%
  mutate(pheno = factor(str_sub(score, 0,3), levels = c("adh","asd","scz"),labels = c("ADHD","ASD","Schizophrenia")),
         score = factor(str_sub(score, start = 6), levels=c("5e.08_resid","1e.06_resid", "0.0001_resid", 
                                                            "0.001_resid", "0.01_resid", "0.05_resid", "0.1_resid", 
                                                            "0.2_resid", "0.5_resid", "1_resid"), labels=c("p<5e-08","p<1e-06","p<0.0001","p<0.001", "p<0.01", 
                                                                                                           "p<0.05","p<0.1","p<0.2","p<0.5","All SNPs")),
         mysize = case_when(score == "p<5e-08" ~ 1,
                            score == "p<1e-06" ~ 2,
                            score == "p<0.0001" ~ 3,
                            score == "p<0.001" ~ 4,
                            score == "p<0.01" ~ 5,
                            score == "p<0.05" ~ 6,
                            score == "p<0.1" ~ 7,
                            score == "p<0.2" ~ 8,
                            score == "p<0.5" ~ 9,
                            score == "All SNPs" ~ 10),
         sex = ifelse(model != "sex_diff", "Male/Female", sex),
         drop = ifelse(sex=="Male/Female" & group==1, 1,0 ),
         var = factor(var, levels = c("bin_mot_18m","bin_lang_3yr",  "concern_dev_exprssd_others"), 
                      labels =c("Motor delays\nat 18 months", "Language delays\nat 3 years","Concern regarding \ndevelopment expressed \nby others at 3 years") ),
         model = factor(model, levels =c("fixed_beta", "sex_diff"), labels=c("Betas constrained\nacross sex", "Fully sex-stratified") )) %>%
  filter(drop==0) %>%
  group_by(score, var) %>%
  mutate(myalpha = case_when(best == "fixed_beta" & model %in% c("Betas constrained\nacross sex")~ 1,
                             best == "fixed_beta" & model %in% c("Fully sex-stratified")~ 0.3,
                             best == "sex_diff" & model %in% c("Fully sex-stratified")~ 1,
                             best == "sex_diff" & model %in% c("Betas constrained\nacross sex")~ 0.3)) %>%
  droplevels()%>%
  mutate(lci = ifelse(myalpha<1,NA,lci),
         uci = ifelse(myalpha<1,NA,uci),
          sex = factor(sex, levels= c("Male", "Male/Female", "Female")))



ests3_plt_maxrsq <- ests3 %>%
  left_join(fitcomps_rsq %>%
              rename(best=model) %>%
              select(score, var, best))%>%
  rename(Beta = est.std,
         uci = ci.upper,
         lci = ci.lower) %>%
  mutate(pheno = factor(str_sub(score, 0,3), levels = c("adh","asd","scz"),labels = c("ADHD","ASD","Schizophrenia")),
         score = factor(str_sub(score, start = 6), levels=c("5e.08_resid","1e.06_resid", "0.0001_resid", 
                                                            "0.001_resid", "0.01_resid", "0.05_resid", "0.1_resid", 
                                                            "0.2_resid", "0.5_resid", "1_resid"), labels=c("p<5e-08","p<1e-06","p<0.0001","p<0.001", "p<0.01", 
                                                                                                           "p<0.05","p<0.1","p<0.2","p<0.5","All SNPs")),
         mysize = case_when(score == "p<5e-08" ~ 1,
                            score == "p<1e-06" ~ 2,
                            score == "p<0.0001" ~ 3,
                            score == "p<0.001" ~ 4,
                            score == "p<0.01" ~ 5,
                            score == "p<0.05" ~ 6,
                            score == "p<0.1" ~ 7,
                            score == "p<0.2" ~ 8,
                            score == "p<0.5" ~ 9,
                            score == "All SNPs" ~ 10),
         sex = ifelse(model != "sex_diff", "Male/Female", sex),
         drop = ifelse(sex=="Male/Female" & group==1, 1,0 ),
         var = factor(var, levels = c("bin_mot_18m","bin_lang_3yr",  "concern_dev_exprssd_others"), 
                      labels =c("Motor delays\nat 18 months", "Language delays\nat 3 years","Concern regarding \ndevelopment expressed \nby others at 3 years") ),
         model = factor(model, levels =c("fixed_beta", "sex_diff"), labels=c("Betas constrained\nacross sex", "Fully sex-stratified") )) %>%
  filter(drop==0) %>%
  group_by(score, var) %>%
  mutate(myalpha = case_when(best == "fixed_beta" & model %in% c("Betas constrained\nacross sex")~ 1,
                             best == "fixed_beta" & model %in% c("Fully sex-stratified")~ 0.3,
                             best == "sex_diff" & model %in% c("Fully sex-stratified")~ 1,
                             best == "sex_diff" & model %in% c("Betas constrained\nacross sex")~ 0.3)) %>%
  droplevels()%>%
  mutate(lci = ifelse(myalpha<1,NA,lci),
         uci = ifelse(myalpha<1,NA,uci),
         sex = factor(sex, levels= c("Male", "Male/Female", "Female"))) %>% 
  drop_na(best)

save(ests3_plt,ests3_plt_maxrsq,fitcomp_best,ests3, prprobs, file = "./scratch_data/dichotomous_analyses_estimates.RData")




# Probit regression betas are not informative in and of themselves,
# so we want to plot predicted probabilities with CIs from best-fitting models

# Predicted probabilities for different PRS

prprobs2 <- prprobs %>%
  mutate(score = rep(params$Var1, each = dim(prprobs)[1]/dim(params)[1]),
         var = rep(params$Var2, each = dim(prprobs)[1]/dim(params)[1]),
         sex = factor(group, levels = c(1,2), labels = c("Females", "Males")),
         i_est = 1-est,
         i_uci = 1-uci,
         i_lci = 1-lci,
         r_prsval = prsval*-1) %>%
  right_join(fitcomp_best %>% select(score, var,model))


exemplar <- prprobs2 %>%
  filter(grepl("adh", score),
         var == "bin_lang_3yr")   %>%
  mutate(score = factor(str_sub(score, start = 6), levels=c("5e.08_resid","1e.06_resid", "0.0001_resid", 
                                                            "0.001_resid", "0.01_resid", "0.05_resid", "0.1_resid", 
                                                            "0.2_resid", "0.5_resid", "1_resid"), labels=c("p<5e-08","p<1e-06","p<0.0001","p<0.001", "p<0.01", 
                                                                                                           "p<0.05","p<0.1","p<0.2","p<0.5","All SNPs")))

tiff(paste0("figures/multigroup_secondary_analysis_predprobs_ex1.tiff"), res = 300, compression = "lzw", unit = "in",
     height = 6, width =6)

ggplot(data= exemplar,aes(x=r_prsval,y=i_est))+
  geom_line(aes(colour=score)) +
  geom_ribbon(aes(ymin=i_uci,ymax=i_lci, fill=score), alpha=0.08) +
  facet_grid(sex ~ .)+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 30, b = , l = 0)),
        text =element_text(size = 14.5),
        panel.background = element_rect(fill = "white", colour = "grey80"),
        panel.spacing = unit(0.2, "lines"),
        panel.grid.major = element_line(colour="white"),
        strip.text = element_text(size=11),
        legend.title = element_text(size =11)) +
  scale_x_continuous(name = "ADHD PRS (standardised)")+
  scale_y_continuous(name = "Model-predicted probability of \nlanguage delays at 3yrs")

dev.off()


exemplar2 <- prprobs2 %>%
  filter(grepl("scz", score),
         var == "concern_dev_exprssd_others")   %>%
  mutate(score = factor(str_sub(score, start = 6), levels=c("5e.08_resid","1e.06_resid", "0.0001_resid", 
                                                            "0.001_resid", "0.01_resid", "0.05_resid", "0.1_resid", 
                                                            "0.2_resid", "0.5_resid", "1_resid"), labels=c("p<5e-08","p<1e-06","p<0.0001","p<0.001", "p<0.01", 
                                                                                                           "p<0.05","p<0.1","p<0.2","p<0.5","All SNPs")))
tiff(paste0("figures/multigroup_secondary_analysis_predprobs_ex2.tiff"), res = 300, compression = "lzw", unit = "in",
     height = 6, width =6)

ggplot(data= exemplar2,aes(x=r_prsval,y=i_est))+
  geom_line(aes(colour=score)) +
  geom_ribbon(aes(ymin=i_uci,ymax=i_lci, fill=score), alpha=0.08) +
  facet_grid(sex ~ .)+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 30, b = , l = 0)),
        text =element_text(size = 14.5),
        panel.background = element_rect(fill = "white", colour = "grey80"),
        panel.spacing = unit(0.2, "lines"),
        panel.grid.major = element_line(colour="white"),
        strip.text = element_text(size=11),
        legend.title = element_text(size =11)) +
  scale_x_continuous(name = "SCZ PRS (standardised)")+
  scale_y_continuous(name = "Model-predicted probability of mother reporting others having \n expressed concern at child's development at 3yrs")


dev.off()

