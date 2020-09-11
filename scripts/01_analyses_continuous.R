load( file = './data/prepped_00_vars.RData')
#install.packages("gtools", repos = "file://tsd-evs/shared/R/cran")

#library(gtools)
library(tidyverse)
library(lavaan)
require('semPlot')
library(semTools)

#Do not restrict to genotyped individuals as fiml can allow us to use incomplete cases

alldata <- alldata %>%
  filter(sex %in% c("Male","Female")) %>% 
  drop_na(IID) %>%
  droplevels()

pval <- 0.05
#Multiple testing correction is now applied post hoc to critical pval thresh

# show underlying bysex associations

#Create "rate of language development" variable
alldata$diff_1st_wds_sns <- alldata$age_1st_sns - alldata$age_1st_wds


ad <- alldata %>% 
  drop_na(sex) %>% 
  group_by( sex) %>%
  select(IID, sex, age_1st_sns, age_1st_wds, age_1st_wlk_18m,age_1st_wlk_3yr,age_1st_wlk,diff_1st_wds_sns, contains("resid"))%>%
  drop_na(IID) %>%
  gather(score, val, -sex, -age_1st_sns,-age_1st_wds, -age_1st_wlk,-age_1st_wlk_18m,-age_1st_wlk_3yr, -diff_1st_wds_sns, -IID) %>%
  ungroup() %>%
  mutate(pheno = factor(str_sub(score, 0,3), levels = c("adh","asd","scz"),labels = c("ADHD","Autism","Schizophrenia")),
         score = factor(str_sub(score, start = 6), levels=c("5e.08_resid","1e.06_resid", "0.0001_resid", 
                                                            "0.001_resid", "0.01_resid", "0.05_resid", "0.1_resid", 
                                                            "0.2_resid", "0.5_resid", "1_resid"), labels=c("p<5e-08","p<1e-06","p<0.0001","p<0.001", "p<0.01", 
                                                                                                           "p<0.05","p<0.1","p<0.2","p<0.5","All SNPs"))) %>%
  filter(pheno %in% c("ADHD","Autism","Schizophrenia"),
         score %in% c( "p<0.001", "p<0.01", 
                       "p<0.05","p<0.1","p<0.2","p<0.5","All SNPs"))
RColorBrewer::brewer.pal(3, "Dark2") # Get palette to match other plots (Males green Females Purple)

gg3 <- 
  
  ggplot(ad, aes(x=val,y=age_1st_wlk, colour= sex, fill = sex)) + 
  geom_smooth(method = "lm", na.rm=T, fullrange = TRUE, alpha=0.3 ) + 
  facet_grid(pheno ~ score) +
  scale_colour_manual(values = c("#1B9E77", "#7570B3"), name="Sex" )+
  scale_fill_manual(values = c("#1B9E77", "#7570B3"), name="Sex")+
  geom_hline( aes(yintercept=mean(age_1st_wlk, na.rm = TRUE)), color = "grey50", linetype = 2)+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 30, b = , l = 0)),
        text =element_text(size = 14.5),
        panel.background = element_rect(fill = "white", colour = "grey80"),
        panel.spacing = unit(0.2, "lines"),
        panel.grid.major = element_line(colour="white"),
        strip.text = element_text(size=11, face = "bold" ),
        strip.text.y = element_text(angle = 0),
        strip.background = element_rect(fill = "grey90", colour = "white"),
        legend.title = element_text(size =11),
        legend.position = "bottom") +
  coord_cartesian(xlim=c(-3.5,3.5), ylim=c(12.5,13.5))+
  scale_x_continuous(name = "PGS (standardised)") + 
  scale_y_continuous(name = "Age at first walking (months)", breaks = c(12.6,12.8,13.0,13.2,13.4))

tiff("figures/walking_prs_bysex.tiff", res = 600, compression = "lzw", unit = "in",
     height = 6, width =6)

gg3  
dev.off()

gg4 <- ggplot(ad, aes(x=val,y=age_1st_wds, colour= sex, fill = sex)) + 
  geom_smooth(method = "lm", na.rm=T, fullrange = TRUE, alpha=0.3 ) + 
  facet_grid(pheno ~ score, switch = "y" ) +
  scale_colour_manual(values = c("#1B9E77", "#7570B3"), name="Sex" )+
  scale_fill_manual(values = c("#1B9E77", "#7570B3"), name="Sex")+
  geom_hline( aes(yintercept=mean(age_1st_wds, na.rm = TRUE)), color = "grey50", linetype = 2)+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 0, b = , l = 60)),
        text =element_text(size = 14.5),
        panel.background = element_rect(fill = "white", colour = "grey80"),
        panel.spacing = unit(0.2, "lines"),
        panel.grid.major = element_line(colour="white"),
        strip.text = element_text(size=11, face = "bold" ),
        strip.text.y = element_text(angle = 180),
        strip.background = element_rect(fill = "grey90", colour = "white"),
        legend.title = element_text(size =11),
        legend.position = "bottom") +
  coord_cartesian(xlim=c(-3.5,3.5), ylim=c(11,16))+
  scale_x_continuous(name = "PGS (standardised)") + 
  scale_y_continuous(name = "Age at first words (months)\n", breaks = c(11,12,13,14,15,16), position="right")


tiff("figures/words_prs_bysex.tiff", res = 600, compression = "lzw", unit = "in",
     height = 6, width =6)
gg4  
dev.off()


# Primary analyses - lavaan linear model with fiml for missing data


names(alldata)

model<- 'y ~ x'
#model2<- 'y ~ x + sex + x_sex_int'

# Omitting lowest PRS threshold as these contain few SNPs for ASD/ADHD

xs <- c("scz_X0.001_resid", "scz_X0.01_resid", "scz_X0.05_resid", "scz_X0.1_resid", 
        "scz_X0.2_resid", "scz_X0.5_resid", "scz_X1_resid",  "asd_X0.001_resid", 
        "asd_X0.01_resid", "asd_X0.05_resid", "asd_X0.1_resid", "asd_X0.2_resid", 
        "asd_X0.5_resid", "asd_X1_resid",   
        "adh_X0.001_resid", "adh_X0.01_resid", "adh_X0.05_resid", 
        "adh_X0.1_resid", "adh_X0.2_resid", "adh_X0.5_resid", "adh_X1_resid")
ys <- c("age_1st_wlk","age_1st_wds", "age_1st_sns","diff_1st_wds_sns" )

ests <- data.frame()
fitcomps <- data.frame()
rsqs <- data.frame() 

params <- expand.grid(xs,ys, stringsAsFactors = FALSE)

for(i in 1:dim(params)[[1]]){ 
  
  alldata[,"x"] <- alldata[,params[i,1]]
  
  alldata[,"y"] <- alldata[,params[i,2]]

  fit <- sem(model, data=alldata, missing="fiml.x", group = "sex")
  fit1 <- sem(model, data=alldata, missing="fiml.x", group = "sex", group.equal = c("regressions"))
  fit2 <- sem(model, data=alldata, missing="fiml.x", group = "sex", group.equal = c("regressions","intercepts"))
  fit3 <- sem(model, data=alldata, missing="fiml.x", group = "sex", group.equal = c("regressions","intercepts","residuals"))
 
  fitcomps <- rbind(fitcomps,lavTestLRT(fit,fit1),
                    lavTestLRT(fit,fit2)[2,], lavTestLRT(fit,fit3)[2,])
  
  tmprsq_vars <- as.data.frame(inspect(fit, "r2"))%>% gather(group, estimate) %>%  mutate(model= "sex_diff") %>% 
    rbind(as.data.frame(inspect(fit1, "r2")) %>% gather(group, estimate) %>%  mutate(model= "fixed_beta")) %>% 
    rbind(as.data.frame(inspect(fit2, "r2")) %>% gather(group, estimate) %>%  mutate(model= "fixed_beta_int")) %>% 
    rbind(as.data.frame(inspect(fit3, "r2")) %>% gather(group, estimate) %>%  mutate(model= "fixed_beta_int_var")) 
  
  rsqs <- rbind(rsqs, tmprsq_vars)
  
  ests <- ests %>%
    bind_rows(standardizedSolution(fit, ci = TRUE, level = 0.95) %>%
                mutate(model = "sex_diff")) %>%
    bind_rows(standardizedSolution(fit1, ci = TRUE, level = 0.95) %>%
                mutate(model = "fixed_beta")) %>%
    bind_rows(standardizedSolution(fit2, ci = TRUE, level = 0.95) %>%
                mutate(model = "fixed_beta_int")) %>%
    bind_rows(standardizedSolution(fit3, ci = TRUE, level = 0.95) %>%
                mutate(model = "fixed_beta_int_var")) 
}

ests1 <- ests %>%
  filter(op %in% c('~')) %>%
  mutate(score = rep(params$Var1, each = 8),
         var = rep(params$Var2, each = 8)) %>%
  mutate(sex = recode(group, `1`="Male",`2`= "Female")) %>%
  select(score,var,everything())

fitcomps1 <- fitcomps %>%
  mutate(score = rep(params$Var1, each = 4),
         var = rep(params$Var2, each = 4),
         pLRT = round(`Pr(>Chisq)`,3),
         model = rep(c("sex_diff","fixed_beta","fixed_beta_int","fixed_beta_int_var"), dim(fitcomps)[1]/4 ) )

fitcomp_best <- fitcomps1  %>%
  group_by(var, score) %>%
  summarise(min = min(pLRT, na.rm = T),
            max = max(pLRT, na.rm = T)) %>%
  mutate(best = ifelse(min > pval, "fixed_beta_int_var", "some_sex_diff")) 
  

fitcomp_best2 <- fitcomps1 %>%
  left_join(fitcomp_best) %>%
  filter(best=="fixed_beta_int_var",
         Df == 3) %>%
  select(score,var,best)

fitcomp_best3 <- fitcomps1 %>%
 left_join(fitcomp_best) %>%
  filter(best=="some_sex_diff") %>%
  mutate(best = ifelse(Df==1 & pLRT<pval, "sex_diff", NA)) %>%
  filter(best == "sex_diff") %>%
  select(score, var, best)

fitcomp_best4 <- fitcomps1 %>%
  anti_join(fitcomp_best2 %>% select(score,var)) %>%
  anti_join(fitcomp_best3 %>% select(score,var)) %>%
  mutate(best = ifelse(Df==2 & pLRT<pval, "fixed_beta", NA ))%>%
  filter(best == "fixed_beta") %>%
  select(score, var, best)
  
fitcomp_best5 <- fitcomps1 %>%
  anti_join(fitcomp_best2 %>% select(score,var)) %>%
  anti_join(fitcomp_best3 %>% select(score,var)) %>%
  anti_join(fitcomp_best4 %>% select(score,var)) %>%
  mutate(best = ifelse(Df==3 & pLRT<pval, "fixed_beta_int", NA )) %>%
  filter(best == "fixed_beta") %>%
  select(score, var, best)

fitcomps_final <- fitcomp_best2 %>%
  bind_rows(fitcomp_best3) %>%
  bind_rows(fitcomp_best4) %>%
  bind_rows(fitcomp_best5) %>%
  rename(model = best) %>%
  left_join(fitcomps1)

# Should be no LRT p-values below pval in the final df:

fitcomps_final %>% filter(pLRT < pval)


fitcomps_out <- fitcomps1 %>% 
  left_join(fitcomps_final %>% 
              select(score,var,best=model)) %>% 
  mutate(selected = ifelse(model ==best, "Y", ""))

# Incorporate rsq info

rsqs_sums <- rsqs %>% 
  mutate(score = rep(params$Var1, each = 8),
         var = rep(params$Var2, each = 8)) %>% 
  group_by(var,score,model) %>% 
  summarise(sum_rsq = sum(estimate))

fitcomps_rsq <- fitcomps_final %>% 
  left_join(rsqs_sums) %>% 
  mutate(pheno = str_sub(score, end=3)) %>% 
  group_by(pheno,var) %>% 
  filter(sum_rsq == max(sum_rsq)) %>% 
  ungroup()

rsqs_bysex <- rsqs %>% 
  mutate(score = rep(params$Var1, each = 8),
         var = rep(params$Var2, each = 8)) %>% 
  right_join(fitcomps_rsq)


save(fitcomps_rsq, fitcomps_out,rsqs_bysex, file = "./scratch_data/continuous_analyses_fitcomps.RData")

# Table with ests from best fitting model at best threshold, and corrected p.values
# corrections = 1) Bonferroni for all possible params of interest (i.e., male + female for all comparisons)
# 2) Benjamini-Hochberg FDR

ests1_tab <- ests1 %>%
  right_join(fitcomps_rsq %>%
              rename(best=model) %>%
              select(score, var, best)) %>%
  group_by(score, var) %>% 
  filter(model == best) %>% 
  ungroup()
  
# Save out in this form for combination with results from 02_analyses_dichotomous
# MTC performed in 05_adjust_pvals

save(ests1_tab, file="./scratch_data/cont_analyses_ests_pre_mtc.RData")


# Plot results

ests1_plt <- ests1 %>%
  left_join(fitcomps_final %>%
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
         x = 1,
         var = factor(var, levels = c("age_1st_wlk", "age_1st_wds", "age_1st_sns", "diff_1st_wds_sns"), 
                      labels = c("Age at first \nwalking", "Age at first \nwords", "Age at first \nsentences", "Rate (age first \nsentences - age \nfirst words)") ),   
         sex = ifelse(model != "sex_diff", "Male/Female", sex),
         drop = ifelse(sex=="Male/Female" & group==1, 1,0 ),
         model = factor(model, levels =c("fixed_beta_int_var","fixed_beta_int", "fixed_beta", "sex_diff"), labels=c("Fully constrained \nacross sex","Betas and intercepts\nconstrained across sex", "Only betas constrained\nacross sex", "Fully sex-stratified") )) %>%
  filter(drop==0) %>%
  droplevels() %>%
  group_by(score, var) %>%
  mutate(drop = case_when(best == "fixed_beta_int_var" & model %in% c("Betas and intercepts\nconstrained across sex", "Only betas constrained\nacross sex")~ 1,
                          best == "sex_diff" & model %in% c("Fully constrained \nacross sex", "Betas and intercepts\nconstrained across sex")~ 1,
                          best == "fixed_beta_int" & model %in% c("Fully constrained \nacross sex", "Only betas constrained\nacross sex")~ 1,
                          best == "fixed_beta" & model %in% c("Fully constrained \nacross sex", "Betas and intercepts\nconstrained across sex")~ 1)) %>%
  filter(is.na(drop)) %>%
  mutate(myalpha = case_when(best == "fixed_beta_int_var" & model %in% c("Fully constrained \nacross sex")~ 1,
                             best == "fixed_beta_int_var" & model %in% c("Fully sex-stratified")~ 0.3,
                             best == "sex_diff" & model %in% c("Fully sex-stratified")~ 1,
                             best == "sex_diff" & model %in% c("Only betas constrained\nacross sex")~ 0.3,
                             best == "fixed_beta" & model %in% c("Only betas constrained\nacross sex")~ 1,
                             best == "fixed_beta" & model %in% c("Fully sex-stratified")~ 0.3,
                             best == "fixed_beta_int" & model %in% c("Betas and intercepts\nconstrained across sex")~ 1,
                             best == "fixed_beta_int" & model %in% c("Fully sex-stratified")~ 0.3)) %>%
  droplevels()%>%
  mutate(lci = ifelse(myalpha<1,NA,lci),
         uci = ifelse(myalpha<1,NA,uci),
         sex = factor(sex, levels= c("Male", "Male/Female", "Female")))
  

ests1_plt_maxrsq <- ests1 %>%
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
         x = 1,
         var = factor(var, levels = c("age_1st_wlk", "age_1st_wds", "age_1st_sns", "diff_1st_wds_sns"), 
                      labels = c("Age at first \nwalking", "Age at first \nwords", "Age at first \nsentences", "Rate (age first \nsentences - age \nfirst words)") ),   
         sex = ifelse(model != "sex_diff", "Male/Female", sex),
         drop = ifelse(sex=="Male/Female" & group==1, 1,0 ),
         model = factor(model, levels =c("fixed_beta_int_var","fixed_beta_int", "fixed_beta", "sex_diff"), labels=c("Fully constrained \nacross sex","Betas and intercepts\nconstrained across sex", "Only betas constrained\nacross sex", "Fully sex-stratified") )) %>%
  filter(drop==0) %>%
  droplevels() %>%
  group_by(score, var) %>%
  mutate(drop = case_when(best == "fixed_beta_int_var" & model %in% c("Betas and intercepts\nconstrained across sex", "Only betas constrained\nacross sex")~ 1,
                          best == "sex_diff" & model %in% c("Fully constrained \nacross sex", "Betas and intercepts\nconstrained across sex")~ 1,
                          best == "fixed_beta_int" & model %in% c("Fully constrained \nacross sex", "Only betas constrained\nacross sex")~ 1,
                          best == "fixed_beta" & model %in% c("Fully constrained \nacross sex", "Betas and intercepts\nconstrained across sex")~ 1)) %>%
  filter(is.na(drop)) %>%
  mutate(myalpha = case_when(best == "fixed_beta_int_var" & model %in% c("Fully constrained \nacross sex")~ 1,
                             best == "fixed_beta_int_var" & model %in% c("Fully sex-stratified")~ 0.3,
                             best == "sex_diff" & model %in% c("Fully sex-stratified")~ 1,
                             best == "sex_diff" & model %in% c("Only betas constrained\nacross sex")~ 0.3,
                             best == "fixed_beta" & model %in% c("Only betas constrained\nacross sex")~ 1,
                             best == "fixed_beta" & model %in% c("Fully sex-stratified")~ 0.3,
                             best == "fixed_beta_int" & model %in% c("Betas and intercepts\nconstrained across sex")~ 1,
                             best == "fixed_beta_int" & model %in% c("Fully sex-stratified")~ 0.3)) %>%
  droplevels()%>%
  mutate(lci = ifelse(myalpha<1,NA,lci),
         uci = ifelse(myalpha<1,NA,uci),
         sex = factor(sex, levels= c("Male", "Male/Female", "Female"))) %>% 
  drop_na(best)

save(ests1_plt,ests1,ests1_plt_maxrsq, file = "./scratch_data/continuous_analyses_estimates.RData")



### Re-run underlying association plots selecting for max rsq

ad_wlk <- ad %>% 
  select(IID,sex,age_1st_wlk, pheno, score,val) %>% 
  right_join(fitcomps_rsq %>% 
               select(score,var) %>% 
               filter(var=="age_1st_wlk") %>% 
               mutate(pheno = factor(str_sub(score, 0,3), levels = c("adh","asd","scz"),labels = c("ADHD","Autism","Schizophrenia")),
                      score = factor(str_sub(score, start = 6), levels=c("5e.08_resid","1e.06_resid", "0.0001_resid", 
                                                                         "0.001_resid", "0.01_resid", "0.05_resid", "0.1_resid", 
                                                                         "0.2_resid", "0.5_resid", "1_resid"), labels=c("p<5e-08","p<1e-06","p<0.0001","p<0.001", "p<0.01", 
                                                                                                                        "p<0.05","p<0.1","p<0.2","p<0.5","All SNPs"))))
wlklabs <- ad_wlk %>% 
  select(pheno,score) %>% 
  distinct()

 p1 <- 
    
    ggplot(ad_wlk, aes(x=val,y=age_1st_wlk, colour= sex, fill = sex)) + 
    geom_smooth(method = "lm", na.rm=T, fullrange = TRUE, alpha=0.3 ) + 
    facet_grid(. ~ pheno) +
    geom_label(data=wlklabs,aes(x=0,y=10.2,label= score), colour="black", fill="white")+
    scale_colour_manual(values = c("#1B9E77", "#7570B3"), name="Sex" )+
    scale_fill_manual(values = c("#1B9E77", "#7570B3"), name="Sex")+
    geom_hline( aes(yintercept=mean(age_1st_wlk, na.rm = TRUE)), color = "grey50", linetype = 2)+
    theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = , l = 0)),
          axis.title.x = element_text(margin = margin(t = 10, r = 0, b = , l = 0)),
          text =element_text(size = 14.5),
          panel.background = element_rect(fill = "white", colour = "grey80"),
          panel.spacing = unit(0.2, "lines"),
          panel.grid.major = element_line(colour="white"),
          strip.text = element_text(size=11, face = "bold" ),
          strip.text.y = element_text(angle = 0),
          strip.background = element_rect(fill = "grey90", colour = "white"),
          legend.title = element_text(size =11),
          legend.position = "bottom") +
    coord_cartesian(xlim=c(-3.5,3.5), ylim=c(10,16))+
    scale_x_continuous(name = "PGS (standardised)") + 
    scale_y_continuous(name = "Age (months)", breaks = seq(10,16,1))+
   ggtitle("First unsupported walking")
  
 p1 
 
 ad_wds <- ad %>% 
   select(IID,sex,age_1st_wds, pheno, score,val) %>% 
   right_join(fitcomps_rsq %>% 
                select(score,var) %>% 
                filter(var=="age_1st_wds") %>% 
                mutate(pheno = factor(str_sub(score, 0,3), levels = c("adh","asd","scz"),labels = c("ADHD","Autism", "Schizophrenia")),
                       score = factor(str_sub(score, start = 6), levels=c("5e.08_resid","1e.06_resid", "0.0001_resid", 
                                                                          "0.001_resid", "0.01_resid", "0.05_resid", "0.1_resid", 
                                                                          "0.2_resid", "0.5_resid", "1_resid"), labels=c("p<5e-08","p<1e-06","p<0.0001","p<0.001", "p<0.01", 
                                                                                                                         "p<0.05","p<0.1","p<0.2","p<0.5","All SNPs"))))
 wdslabs <- ad_wds %>% 
   select(pheno,score) %>% 
   distinct()
 
 p2 <- 
   
   ggplot(ad_wds, aes(x=val,y=age_1st_wds, colour= sex, fill = sex)) + 
   geom_smooth(method = "lm", na.rm=T, fullrange = TRUE, alpha=0.3 ) + 
   facet_grid(. ~ pheno) +
   geom_label(data=wlklabs,aes(x=0,y=10.2,label= score), colour="black", fill="white")+
   scale_colour_manual(values = c("#1B9E77", "#7570B3"), name="Sex" )+
   scale_fill_manual(values = c("#1B9E77", "#7570B3"), name="Sex")+
   geom_hline( aes(yintercept=mean(age_1st_wds, na.rm = TRUE)), color = "grey50", linetype = 2)+
   theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = , l = 0)),
         axis.title.x = element_text(margin = margin(t = 10, r = 0, b = , l = 0)),
         text =element_text(size = 14.5),
         panel.background = element_rect(fill = "white", colour = "grey80"),
         panel.spacing = unit(0.2, "lines"),
         panel.grid.major = element_line(colour="white"),
         strip.text = element_text(size=11, face = "bold" ),
         strip.text.y = element_text(angle = 0),
         strip.background = element_rect(fill = "grey90", colour = "white"),
         legend.title = element_text(size =11),
         legend.position = "bottom") +
   coord_cartesian(xlim=c(-3.5,3.5), ylim=c(10,16))+
   scale_x_continuous(name = "PGS (standardised)") + 
   scale_y_continuous(name = "Age (months)", breaks = seq(10,16,1))+
   ggtitle("First use of words")
 
 p2 

 ad_sns <- ad %>% 
   select(IID,sex,age_1st_sns, pheno, score,val) %>% 
   right_join(fitcomps_rsq %>% 
                select(score,var) %>% 
                filter(var=="age_1st_sns") %>% 
                mutate(pheno = factor(str_sub(score, 0,3), levels = c("adh","asd","scz"),labels = c("ADHD","Autism", "Schizophrenia")),
                       score = factor(str_sub(score, start = 6), levels=c("5e.08_resid","1e.06_resid", "0.0001_resid", 
                                                                          "0.001_resid", "0.01_resid", "0.05_resid", "0.1_resid", 
                                                                          "0.2_resid", "0.5_resid", "1_resid"), labels=c("p<5e-08","p<1e-06","p<0.0001","p<0.001", "p<0.01", 
                                                                                                                         "p<0.05","p<0.1","p<0.2","p<0.5","All SNPs"))))
 snslabs <- ad_sns %>% 
   select(pheno,score) %>% 
   distinct()
 
 p3 <- 
   
   ggplot(ad_sns, aes(x=val,y=age_1st_sns, colour= sex, fill = sex)) + 
   geom_smooth(method = "lm", na.rm=T, fullrange = TRUE, alpha=0.3 ) + 
   facet_grid(. ~ pheno) +
   geom_label(data=wlklabs,aes(x=0,y=17.2,label= score), colour="black", fill="white")+
   scale_colour_manual(values = c("#1B9E77", "#7570B3"), name="Sex" )+
   scale_fill_manual(values = c("#1B9E77", "#7570B3"), name="Sex")+
   geom_hline( aes(yintercept=mean(age_1st_sns, na.rm = TRUE)), color = "grey50", linetype = 2)+
   theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = , l = 0)),
         text =element_text(size = 14.5),
         panel.background = element_rect(fill = "white", colour = "grey80"),
         panel.spacing = unit(0.2, "lines"),
         panel.grid.major = element_line(colour="white"),
         strip.text = element_text(size=11, face = "bold" ),
         strip.text.y = element_text(angle = 0),
         strip.background = element_rect(fill = "grey90", colour = "white"),
         legend.title = element_text(size =11),
         legend.position = "bottom") +
   coord_cartesian(xlim=c(-3.5,3.5), ylim=c(17,23))+
   scale_x_continuous(name = "PGS (standardised)") + 
   scale_y_continuous(name = "Age (months)", breaks = seq(17,23,1))+
   ggtitle("First use of sentences")
 
 p3 
 

 ad_rate <- ad %>% 
   select(IID,sex,diff_1st_wds_sns, pheno, score,val) %>% 
   right_join(fitcomps_rsq %>% 
                select(score,var) %>% 
                filter(var=="diff_1st_wds_sns") %>% 
                mutate(pheno = factor(str_sub(score, 0,3), levels = c("adh","asd","scz"),labels = c("ADHD","Autism", "Schizophrenia")),
                       score = factor(str_sub(score, start = 6), levels=c("5e.08_resid","1e.06_resid", "0.0001_resid", 
                                                                          "0.001_resid", "0.01_resid", "0.05_resid", "0.1_resid", 
                                                                          "0.2_resid", "0.5_resid", "1_resid"), labels=c("p<5e-08","p<1e-06","p<0.0001","p<0.001", "p<0.01", 
                                                                                                                         "p<0.05","p<0.1","p<0.2","p<0.5","All SNPs"))))
 ratelabs <- ad_rate %>% 
   select(pheno,score) %>% 
   distinct()
 
 p4 <- 
   
   ggplot(ad_rate, aes(x=val,y=diff_1st_wds_sns, colour= sex, fill = sex)) + 
   geom_smooth(method = "lm", na.rm=T, fullrange = TRUE, alpha=0.3 ) + 
   facet_grid(. ~ pheno) +
   geom_label(data=ratelabs,aes(x=0,y=4.2,label= score), colour="black", fill="white")+
   scale_colour_manual(values = c("#1B9E77", "#7570B3"), name="Sex" )+
   scale_fill_manual(values = c("#1B9E77", "#7570B3"), name="Sex")+
   geom_hline( aes(yintercept=mean(diff_1st_wds_sns, na.rm = TRUE)), color = "grey50", linetype = 2)+
   theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = , l = 0)),
         text =element_text(size = 14.5),
         panel.background = element_rect(fill = "white", colour = "grey80"),
         panel.spacing = unit(0.2, "lines"),
         panel.grid.major = element_line(colour="white"),
         strip.text = element_text(size=11, face = "bold" ),
         strip.text.y = element_text(angle = 0),
         strip.background = element_rect(fill = "grey90", colour = "white"),
         legend.title = element_text(size =11),
         legend.position = "bottom") +
   coord_cartesian(xlim=c(-3.5,3.5), ylim=c(4,10))+
   scale_x_continuous(name = "PGS (standardised)") + 
   scale_y_continuous(name = "Lag (months)", breaks = seq(4,10,1))+
   ggtitle("Rate of language development \n(age sentences - age words)")
 
 p4 

 
 require(patchwork)
 
 patchw <-p1 + p2 + p3 + p4 & theme(legend.direction="vertical" ,
                                    plot.title = element_text(size = 14, face = "bold"),
                                    strip.text = element_text(angle = 0, colour="white"),
                                    strip.background = element_rect(fill = "grey30", colour = "white"),
                                    panel.background = element_rect(fill = "white", colour = "grey30"),
                                    panel.grid.major = element_line(colour="grey90"))
 
 
 tiff("figures/patch_bysex_all_rsqbest.tiff", res = 600, compression = "lzw", unit = "in",
      height = 9, width =12)
 
 
 patchw + 
   plot_layout(guides="collect") +
   plot_annotation(tag_levels = 'A',
                   caption = "Note: Within-panel labels show PGS threshold selected for maximising model R-squared across sex; \n
                   Y-axes differ but range remains of constant size to facilitate comparison of effect sizes")
   
 
 dev.off()
 