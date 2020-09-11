#06_validation_sensitivity.R

#install.packages("powerAnalysis", repos = "file://tsd-evs/shared/R/cran")
library(powerAnalysis)
library(MASS)
library(yarrr)
library(tidyverse)
library(foreign)
library(lavaan)


load( file = './data/prepped_00_vars.RData')
q5yrs <- read.spss("N:/data/durable/data/MoBaPhenoData/PDB2306_MoBa_V12/SPSS/PDB2306_Q5yrs_v12.sav", to.data.frame = TRUE) 
q8yrs <- read.spss("N:/data/durable/data/MoBaPhenoData/PDB2306_MoBa_V12/SPSS/PDB2306_Q8yrs_v12.sav", to.data.frame = TRUE) 

## Check associations between early development measures and mat-rep NDDs/problems age 5/8

# Get mat-rep NDDs/problems
q5yrsel <- q5yrs %>%
  select(PREG_ID_2306,BARN_NR, LL76,LL77,  #Motor (2=Yes,2=Doctor confirmed)
         LL79,LL80,  #Language
         LL82, LL83, #Hyperactivity/ADHD
         LL85, LL86, #Attention
         LL88, LL89, #Autism
         LL91, LL92) #Aspergers

q8yrsel <- q8yrs %>%
  select(PREG_ID_2306,BARN_NR, NN37, NN38, NN39, #Psychomotor (Current, Ever, Specialist = 2)
                               NN41, NN42, NN43, #Language
                               NN45, NN46, NN47, #Hyperactivity
                               NN49, NN50, NN51, #Attention
                               NN53, NN54, NN55) #Autism

ndd_probs <- q5yrsel %>% 
  full_join(q8yrsel) %>% 
  mutate(motor = case_when(LL76=="Yes"|NN37==1|NN38==1 ~ 1,
                           TRUE~0),
         motor_spec = case_when(LL77=="Yes"|NN39=="Yes" ~ 1,
                             TRUE~0),
         lang = case_when(LL79=="Yes"|NN41==1|NN42==1 ~ 1,
                       TRUE~0),
         lang_spec = case_when(LL80=="Yes"|NN43=="Yes" ~ 1,
                            TRUE~0),
         adhd = case_when(LL82=="Yes"|LL85=="Yes"|NN45==1|NN46==1|NN49==1|NN50==1 ~ 1,
                       TRUE~0),
         adhd_spec = case_when(LL83=="Yes"|LL86=="Yes"|NN47=="Yes"|NN51=="Yes" ~ 1,
                            TRUE~0),
         asd = case_when(LL88=="Yes"|LL91=="Yes"|NN53==1|NN54==1 ~ 1,
                      TRUE~0),
         asd_spec = case_when(LL89=="Yes"|LL92=="Yes"|NN55=="Yes" ~ 1,
                           TRUE~0))
alldata$diff_1st_wds_sns <- alldata$age_1st_sns - alldata$age_1st_wds
valdata <- alldata %>% 
  left_join(ndd_probs) 
valid_tests <- list()
# lang
lang_langspec <- table(valdata$bin_lang_3yr, valdata$lang_spec)
ES.chisq.assoc(lang_langspec)
lang_adhdspec <- table(valdata$bin_lang_3yr, valdata$adhd_spec)
ES.chisq.assoc(lang_adhdspec)
lang_asdsepc <- table(valdata$bin_lang_3yr, valdata$asd_spec)
ES.chisq.assoc(lang_asdsepc)
# mot
mot_motspec <- table(valdata$bin_mot_18m, valdata$motor_spec)
ES.chisq.assoc(mot_motspec)
mot_adhdspec <- table(valdata$bin_mot_18m, valdata$adhd_spec)
ES.chisq.assoc(mot_adhdspec)
mot_asdsepc <- table(valdata$bin_mot_18m, valdata$asd_spec)
ES.chisq.assoc(mot_asdsepc)

#GLMs for age at variables
valdata <- valdata %>%
  filter(sex %in% c("Male","Female")) %>%
  drop_na(IID) %>%
  droplevels()
xs <- c("age_1st_wlk", "age_1st_wds", "age_1st_sns", "diff_1st_wds_sns")
ys <- c(  "motor_spec",  "lang_spec", 
          "adhd_spec", "asd_spec")

ests2 <- data.frame()

params <- expand.grid(xs,ys, stringsAsFactors = FALSE)

fitcomps <- data.frame()
prprobs <- data.frame()
rsqs <- data.frame() 

model<- 'y ~ x'

for(i in 1:dim(params)[[1]]){ 
  
  valdata[,"x"] <- valdata[,params[i,1]]
  valdata[,"y"] <- valdata[,params[i,2]]
  
  fit <- sem(model, data=valdata, ordered= c("y"),group = "sex", link="probit")
  fit1 <- sem(model, data=valdata, ordered= c("y"),group = "sex", link="probit", group.equal = c("regressions"))
  
  
  fitcomps <- rbind(fitcomps,lavTestLRT(fit,fit1))
  
  tmprsq_vars <- as.data.frame(inspect(fit, "r2"))%>% gather(group, estimate) %>%  mutate(model= "sex_diff") %>% 
    rbind(as.data.frame(inspect(fit1, "r2")) %>% gather(group, estimate) %>%  mutate(model= "fixed_beta")) 
  
  rsqs <- rbind(rsqs, tmprsq_vars)
  
  
  res <- standardizedSolution(fit, ci = TRUE, level = 0.95)
  res1 <- standardizedSolution(fit1, ci = TRUE, level = 0.95)
  #  res2 <- standardizedSolution(fit2, ci = TRUE, level = 0.95)
  
  ests2 <- ests2 %>%
    bind_rows(res %>%
                mutate(model = "Sex diffs")) %>%
    bind_rows(res1 %>%
                mutate(model = "No sex diffs")) 
  
}
ests3 <- ests2 %>%
  filter(op %in% c('~')) %>% 
  mutate(predictor = rep(params$Var1, each = 4),
         outcome = rep(params$Var2, each = 4)) %>%
  mutate(sex = recode(group, `1`="Male",`2`= "Female")) %>%
  select(predictor,outcome,everything())


fitcomps1 <- fitcomps %>%
  mutate(predictor = rep(params$Var1, each = 2),
         outcome = rep(params$Var2, each = 2),
         pLRT = round(`Pr(>Chisq)`,2),
         model = rep(c("sex_diff","fixed_beta"), dim(fitcomps)[1]/2 ) ) %>%
  select(-AIC,-BIC, -`Pr(>Chisq)`)

fitcomp_best <- fitcomps1  %>%
  group_by(outcome, predictor) %>%
  summarise(min = min(pLRT, na.rm = T)) %>%
  mutate(model = ifelse(min > 0.04, "No sex diffs", "Sex diffs")) %>%
  select(predictor,outcome,model) %>%
  distinct()

est_fra_best <-fitcomp_best %>% 
  left_join(ests3) %>% 
  select(Predictor=predictor, Outcome=outcome, model, sex, Beta=est.std, ci.lower, ci.upper,raw_Pval=pvalue)

save(est_fra_best, file="./scratch_data/sens2_glms_fitted.RData")


# Dichotomise "age at" variables and run as probit regressions

#Create "rate of language development" variable
alldata <- valdata 


alldata$diff_1st_wds_sns <- alldata$age_1st_sns - alldata$age_1st_wds

wlk <- quantile(alldata$age_1st_wlk, c(0.03,0.97), na.rm=T)
wds <- quantile(alldata$age_1st_wds, c(0.03,0.97), na.rm=T)
sns <- quantile(alldata$age_1st_sns, c(0.03,0.97), na.rm=T)
rate <- quantile(alldata$diff_1st_wds_sns, c(0.03,0.97), na.rm=T)

alldata <- alldata %>% 
  mutate(early_wlk= case_when(age_1st_wlk<=wlk[1] ~ 1,
                               age_1st_wlk>wlk[1] ~ 0),
         late_wlk= case_when( age_1st_wlk<wlk[2] ~ 0,
                               age_1st_wlk>=wlk[2] ~ 1),
         early_wds= case_when(age_1st_wds<=wds[1] ~ 1,
                              age_1st_wds>wds[1] ~ 0),
         late_wds= case_when( age_1st_wds<wds[2] ~ 0,
                              age_1st_wds>=wds[2] ~ 1),
         early_sns= case_when(age_1st_sns<=sns[1] ~ 1,
                              age_1st_sns>sns[1] ~ 0),
         late_sns= case_when( age_1st_sns<sns[2] ~ 0,
                              age_1st_sns>=sns[2] ~ 1),
         slow_rate= case_when(diff_1st_wds_sns<=rate[1] ~ 1,
                              diff_1st_wds_sns>rate[1] ~ 0),
         fast_rate= case_when(diff_1st_wds_sns<rate[2] ~ 0,
                              diff_1st_wds_sns>=rate[2] ~ 1))

latesns_lang <- table(alldata$bin_lang_3yr, alldata$late_sns)
ES.chisq.assoc(latesns_lang)
fisher.test(latesns_lang)
latewds_lang <- table(alldata$bin_lang_3yr, alldata$late_wds)
ES.chisq.assoc(latewds_lang)
fisher.test(latewds_lang)

earlysns_lang <- table(alldata$bin_lang_3yr, alldata$early_sns)
ES.chisq.assoc(earlysns_lang)
fisher.test(earlysns_lang)
earlywds_lang <- table(alldata$bin_lang_3yr, alldata$early_wds)
ES.chisq.assoc(earlywds_lang)
fisher.test(earlywds_lang)

slow_rate_lang <- table(alldata$bin_lang_3yr, alldata$slow_rate)
ES.chisq.assoc(slow_rate_lang)
fisher.test(slow_rate_lang)
fast_rate_lang <- table(alldata$bin_lang_3yr, alldata$fast_rate)
ES.chisq.assoc(fast_rate_lang)
fisher.test(fast_rate_lang)

latewlk_mot <- table(alldata$bin_mot_18m, alldata$late_wlk)
ES.chisq.assoc(latewlk_mot)
fisher.test(latewlk_mot)
earlywlk_mot <- table(alldata$bin_mot_18m, alldata$early_wlk)
ES.chisq.assoc(earlywlk_mot)
fisher.test(earlywlk_mot)

save(alldata, file = './data/sensitivity_06_vars.RData')

# Re-reun all analyses excluding individuals with mat-rep diagnoses