#install.packages("gtools", repos = "file://tsd-evs/shared/R/cran")

library(tidyverse)


load( file = './data/prepped_00_vars.RData')

seldata <- alldata %>%
  select(PREG_ID_2306, BARN_NR, sex, IID, matches("resid"),age_1st_wlk,denom_18m, denom_3yr,denom_5yr, age_1st_wds,bin_lang_3yr,
         age_1st_sns, bin_mot_18m,concern_dev_exprssd_others) %>%
  mutate(label = case_when(!is.na(IID) ~ "Geno",
                           is.na(IID) ~"noGeno")) %>%
  mutate(label2 = case_when(!is.na(age_1st_wds) ~ "asked_AAFW",
                            is.na(age_1st_wds)&(!is.na(denom_5yr)) ~"notasked_AAFW")) %>%
  filter(sex %in% c("Male","Female")) %>% 
  droplevels()


Genodat<- seldata %>% 
  filter(label=="Geno") %>% 
  select(-matches("resid"))

psych::describeBy(Genodat, group="sex", mat=TRUE)%>% 
  rownames_to_column() %>% 
  select(var= rowname, sex=group1, n, mean, sd,min, max)


b <- data.frame(table(Genodat$sex,Genodat$bin_mot_18m)) %>% 
  bind_rows(data.frame(table(Genodat$sex,Genodat$bin_lang_3yr)),
            data.frame(table(Genodat$sex,Genodat$concern_dev_exprssd_others))) %>% 
  rename(sex= Var1,yn = Var2) %>% 
  mutate(var = rep(c("bin_mot_18m","bin_lang_3yr","concern_dev_exprssd_others"), each=4)) %>%
  filter(yn==1)

noGenodat<- seldata %>% 
  filter(label=="noGeno")%>% 
  select(-matches("resid"))

psych::describe(noGenodat)



genomales <- Genodat %>% 
  filter(sex=="Male")

genofemales <- Genodat %>% 
  filter(sex=="Female")

conc_summ <- Genodat %>%
  mutate(concern_dev_exprssd_others=as.integer(concern_dev_exprssd_others)-1) %>% 
  group_by(sex) %>% 
  summarise(sum_conc=sum(concern_dev_exprssd_others, na.rm=T),
            n_conc=sum(!is.na(concern_dev_exprssd_others)))
  
prop.test(x=c(168,86),n=c(4808,4610))

# Testing selection on key vars

ttests1 <- list(lapply(split(seldata,seldata$sex),function(x)with(x, t.test(age_1st_wlk~label))),
                lapply(split(seldata,seldata$sex),function(x)with(x, t.test(age_1st_wds~label))),
                lapply(split(seldata,seldata$sex),function(x)with(x, t.test(age_1st_sns~label))))

selection_cont <- sapply(ttests1, function(x) {
  c(x$Male$estimate[1],
    x$Male$estimate[2],
    x$Male$estimate[1]-x$Male$estimate[2],
    ci.lower = x$Male$conf.int[1],
    ci.upper = x$Male$conf.int[2],
    p.value = x$Male$p.value)
}) %>% as.data.frame()  %>% 
  bind_cols(sapply(ttests1, function(x) {
    c(x$Female$estimate[1],
      x$Female$estimate[2],
      x$Female$estimate[1]-x$Female$estimate[2],
      ci.lower = x$Female$conf.int[1],
      ci.upper = x$Female$conf.int[2],
      p.value = x$Female$p.value)
  }) %>% as.data.frame()) %>% 
  t()%>% as.data.frame() %>% 
  mutate(sex=c(rep("Male",3),rep("Female",3)),
         var=rep(c("age_1st_wlk","age_1st_wds","age_1st_sns"),2)) %>%
  `colnames<-`(c("mean_geno", "mean_nogeno", "t_stat","lci","uci","pval","sex","var")) 



# Prportions

summdata <- seldata %>%
  select(PREG_ID_2306,BARN_NR,sex,bin_lang_3yr,bin_mot_18m,concern_dev_exprssd_others, label) %>%
  mutate(concern_dev_exprssd_others=as.integer(concern_dev_exprssd_others)) %>% 
  group_by(sex, label) %>% 
  summarise(sum_lang=sum(bin_lang_3yr, na.rm=T),
            n_lang=sum(!is.na(bin_lang_3yr)),
            sum_mot=sum(bin_mot_18m, na.rm=T),
            n_mot=sum(!is.na(bin_mot_18m)),
            sum_conc=sum(concern_dev_exprssd_others, na.rm=T),
            n_conc=sum(!is.na(concern_dev_exprssd_others))) %>% 
  gather(key=stat, value=value, -sex,-label) %>%
  separate(stat, into=c("type","pheno")) %>%
  unite(label_type, c("label","type")) %>% 
  spread(key=label_type, value=value) 
summdata

selection2 <- data.frame()
for(i in 1:dim(summdata)[1]){
  a <- prop.test(x=c(summdata[i,]$Geno_sum, summdata[i,]$noGeno_sum),
                 n= c(summdata[i,]$Geno_n, summdata[i,]$noGeno_n))
  selection2 <- rbind(selection2,
                      c(a$estimate, a$p.value, a$conf.int))
}
colnames(selection2)<- c("prop_Geno","prop_noGeno", "pval","lower.ci","upper.ci")

selection_proportions <- selection2  %>% 
  mutate(pheno= summdata$pheno,
         sex= summdata$sex, 
         diff_prop = prop_Geno-prop_noGeno)

# Prportions language delays per availability of age at first talking Q
q5yrs <- read.spss("N:/data/durable/data/MoBaPhenoData/PDB2306_MoBa_V12/SPSS/PDB2306_Q5yrs_v12.sav", to.data.frame = TRUE) %>% 
  select(PREG_ID_2306,BARN_NR,VERSJON_SKJEMA_5AAR_TBL1) 



summdata2 <- seldata %>%
  left_join(q5yrs) %>% 
  select(PREG_ID_2306,BARN_NR,sex,bin_lang_3yr,VERSJON_SKJEMA_5AAR_TBL1) %>%
  drop_na(VERSJON_SKJEMA_5AAR_TBL1) %>% 
  mutate(q_version = case_when(VERSJON_SKJEMA_5AAR_TBL1 == "SKJEMA_5AAR                                                                " ~ "AAFW",
                               VERSJON_SKJEMA_5AAR_TBL1 == "SKJEMA_5AARB                                                               " ~"noAAFW")) %>% 
  select(-VERSJON_SKJEMA_5AAR_TBL1) %>% 
  group_by(q_version) %>% 
  summarise(sum_lang=sum(bin_lang_3yr, na.rm=T),
            n_lang=sum(!is.na(bin_lang_3yr))) %>% 
  gather(key=stat, value=value, -q_version) %>%
  unite(q_vrsn_stat, c("q_version","stat")) %>% 
  spread(key=q_vrsn_stat, value)
summdata2

results2 <- data.frame()

for(i in 1:dim(summdata2)[1]){
  a <- prop.test(x=c(summdata2[i,]$AAFW_sum_lang, summdata2[i,]$noAAFW_sum_lang),
                 n= c(summdata2[i,]$AAFW_n_lang, summdata2[i,]$noAAFW_n_lang))
  results2 <- rbind(results2,
                    c(a$estimate, a$p.value, a$conf.int))
}
colnames(results2)<- c("prop_AAFW","prop_noAAFW", "pval","lower.ci","upper.ci")

aafw_selection <- results2  %>% 
  mutate(diff_prop = prop_AAFW-prop_noAAFW)


# Prportions language delays by sex

summdata3 <- Genodat %>%
  select(PREG_ID_2306,BARN_NR,sex,bin_lang_3yr) %>%
  drop_na(bin_lang_3yr) %>% 
  summarise(sum_lang=sum(bin_lang_3yr, na.rm=T),
            n_lang=sum(!is.na(bin_lang_3yr))) 
summdata3




a <- prop.test(x=c(281),
               n= c(9486))


colnames(results2)<- c("prop_AAFW","prop_noAAFW", "pval","lower.ci","upper.ci")

aafw_selection <- results2  %>% 
  mutate(diff_prop = prop_AAFW-prop_noAAFW)

# Testing selective attrition - analyses are not longitudinal so omitting this unless asked

## restricting only to those always invited

status <- foreign::read.spss("N:/data/durable/data/MoBaPhenoData/PDB2306_MoBa_V12/SPSS/PDB2306_status_v12.sav", to.data.frame = TRUE) %>% 
  select(PREG_ID_2306, status_skj5_v12, status_skj6_v12, status_skj5y_v12, status_skj8y_v12) %>% 
  filter_at(vars(matches("status")), all_vars(. %in% c("x", "1", "2","3")) ) %>% 
  mutate_at(vars(matches("status")), funs(factor(.))) %>% 
  droplevels() %>% 
  mutate(data_18m =  fct_recode(status_skj5_v12, yes="x", yes="1", no = "2", no="3"),
         data_3yr =  fct_recode(status_skj6_v12, yes="x", yes="1", no = "2", no="3"),
         data_5yr =  fct_recode(status_skj5y_v12, yes="x", yes="1", no = "2", no="3"),
         data_8yr =  fct_recode(status_skj8y_v12, yes="x", yes="1", no = "2", no="3")) %>% 
  right_join(seldata %>% 
               select(PREG_ID_2306, matches("resid")))

seladata <-status



ttests3 <-list(  t.test(seladata$adh_X1_resid~seladata$data_18m),
                 t.test(seladata$asd_X1_resid~seladata$data_18m),
                 t.test(seladata$scz_X1_resid~seladata$data_18m),
                 t.test(seladata$adh_X1_resid~seladata$data_3yr),
                 t.test(seladata$asd_X1_resid~seladata$data_3yr),
                 t.test(seladata$scz_X1_resid~seladata$data_3yr),
                 t.test(seladata$adh_X1_resid~seladata$data_5yr),
                 t.test(seladata$asd_X1_resid~seladata$data_5yr),
                 t.test(seladata$scz_X1_resid~seladata$data_5yr))


selective_attrition <- sapply(ttests3, function(x) {
  c(x$estimate[1],
    x$estimate[2],
    x$estimate[1]-x$estimate[2],
    ci.lower = x$conf.int[1],
    ci.upper = x$conf.int[2],
    p.value = x$p.value)
})  %>% 
  t() %>%
  as.data.frame() 

colnames(selective_attrition) <- c("mean_data", "mean_nodata", "t_stat","lci","uci","pval")


selective_attrition <- selective_attrition %>%
  mutate(wave=rep(c("18m","3yr","5yr"), each=3),
         PGS=rep(c("adhd","asd","scz"),3))  



##Save out


save(aafw_selection, selection_proportions,selection_cont,selective_attrition,file = "./output/selection.RData")

