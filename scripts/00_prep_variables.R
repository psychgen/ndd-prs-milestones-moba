# Read in and prepare data for milestones prs analysis

library(foreign)
library(tidyverse)


q4 <- read.spss("N:/data/durable/data/MoBaPhenoData/PDB2306_MoBa_V12/SPSS/PDB2306_Q4_6months_v12.sav", to.data.frame = TRUE) 
q5 <- read.spss("N:/data/durable/data/MoBaPhenoData/PDB2306_MoBa_V12/SPSS/PDB2306_Q5_18months_v12.sav", to.data.frame = TRUE) 
q6 <- read.spss("N:/data/durable/data/MoBaPhenoData/PDB2306_MoBa_V12/SPSS/PDB2306_Q6_3yrs_v12.sav", to.data.frame = TRUE) 
q5yrs <- read.spss("N:/data/durable/data/MoBaPhenoData/PDB2306_MoBa_V12/SPSS/PDB2306_Q5yrs_v12.sav", to.data.frame = TRUE) 
mbrn <- read.spss("N:/data/durable/data/MoBaPhenoData/PDB2306_MoBa_V12/SPSS/PDB2306_MBRN_541_v12.sav", to.data.frame = TRUE) 

## PRIMARY ANALYSIS VARIABLES

# Age at first walking

# Walk unaided 18months Y/N; version CDE	EE986
# -	If yes, age at which(months); version CDE	EE987
# Still not walking unaided 18months; version AB	EE401
# -	Is walking unaided, age at which (months) version AB	EE400

q5sel <- q5 %>%
  select(PREG_ID_2306,BARN_NR, EE986, EE987, EE401, EE400)

table(q5sel$EE986, q5sel$EE401, useNA = "ifany")

# Exclusions:
# Not walking at 18months (no on EE986, 1 on EE401, declared walking age >18months)
# More than one checkbox
# Implausible values? (some people responding in years? drop <6) 
# + if EE987 is missing, use EE400

q5sel <- q5sel %>%
  filter(!EE986 %in% c("No", "More than 1 check box filled in"),
         is.na(EE401)) %>%
  mutate(age_1st_wlk_18m = ifelse(!is.na(EE987), EE987, EE400)) %>%
  filter(age_1st_wlk_18m >=6) %>%
  select(PREG_ID_2306,BARN_NR, age_1st_wlk_18m) %>%
  mutate(bin_1st_wlk_15m = case_when(age_1st_wlk_18m>=15 ~ 1,
                                     age_1st_wlk_18m<15 ~ 0))


# Age at first walking (II)

# Age at first walking 36months  	GG27
# Still not walking 36months GG28

q6wlk <- q6 %>%
  select(PREG_ID_2306, BARN_NR,GG27,GG28) %>%
  mutate(age_1st_wlk_3yr = GG27) %>%
  filter(age_1st_wlk_3yr >=6) 

q5q6wlk <- q5sel %>%
  full_join(q6wlk) %>%
  mutate(age_1st_wlk = case_when(!is.na(age_1st_wlk_18m) & !is.na(age_1st_wlk_3yr) ~ round(((age_1st_wlk_18m+age_1st_wlk_3yr)/2),0),
                                 !is.na(age_1st_wlk_18m) & is.na(age_1st_wlk_3yr) ~age_1st_wlk_18m,
                                 is.na(age_1st_wlk_18m) & !is.na(age_1st_wlk_3yr) ~age_1st_wlk_3yr)) %>% 
  select(PREG_ID_2306,BARN_NR,age_1st_wlk_18m,age_1st_wlk_3yr,age_1st_wlk)

# Age at first talking

# Language development stage 36months	GG226
# Age first words (5yr) version A 	LL117
# -	Still not using words version A	LL118
# Age first words+sentences (5yr) version A 	LL119
# -	Still not using words+sentences version A	LL120
# First word <2years (5yr) version B	LL465
# Words+sentences<2.5yrs(5yr) version B	LL466
# Språk20 (LL190-212) and CCC-2 (LL213-221 &480-483) for improving imputation 


q5yrsel <- q5yrs %>%
  select(PREG_ID_2306,BARN_NR, LL117, LL118, LL119, LL120, LL465, LL466, 
         paste0("LL", seq(190,221,1)), LL480, LL481, LL482, LL483)

# Exclusions:
# Still not using words (1 on LL118 or declared first usage >60months)
# Implausible values? (some people responding in years? drop <6) 

q5yrsel <- q5yrsel %>%
  filter(is.na(LL118)) %>%
  filter(!LL465 %in% "More than 1 check box filled in") %>%
  mutate(age_1st_wds = ifelse(LL117>=6, LL117,NA) ,
         age_1st_sns = ifelse(LL119>=6, LL119,NA) ,
         bin_1st_wds_2yr = case_when(LL117<24 ~ 0,
                                     LL117>=24 ~ 1,
                                     LL465=="No" ~ 1,
                                     LL465=="Yes" ~0),
         bin_1st_sns_2yr = case_when(LL119<30 ~ 0,
                                     LL119<30 ~ 1,
                                     LL466=="No" ~ 1,
                                     LL466=="Yes" ~ 0),
         denom_5yr=1)  

table(q5yrsel$bin_1st_wds_2yr)

rvrsd <- c("LL216","LL218","LL219","LL220","LL221","LL480", "LL481", "LL482", "LL483")
ccc2 <- c(rvrsd, "LL213", "LL214", "LL215", "LL217") 
q5yrsel_scales <- q5yrsel %>%
  select(PREG_ID_2306,BARN_NR, paste0("LL", seq(190,221,1)), LL480, LL481, LL482, LL483) %>%
  gather(item,val,-PREG_ID_2306,-BARN_NR) %>%
  mutate(true_val = case_when(val %in% c("Completely wrong (1)","No","Seldom or never","Sometimes") ~ 1,
                              val %in% c("Yes","Regularily","Often/ Always","(2)")~ 2,
                              val %in% c("Both right and wrong (3)")~ 3,
                              val %in% c("(4)")~ 4,
                              val %in% c("Completely right (5)")~ 2,
                              val %in% c("Agree")~ 5),
         sscale = ifelse(item %in% ccc2, "ccc2","sprk20")) %>%
  mutate(true_val = ifelse(item %in% rvrsd, ((-1*true_val)+3),true_val))


q5yrsel_scales <- q5yrsel_scales %>% 
  group_by( sscale) %>%
  summarize(items_scale = length(unique(item))) %>%
  right_join(q5yrsel_scales) %>%
  group_by(PREG_ID_2306,BARN_NR, sscale, items_scale) %>%
  summarize(items_present = sum(!is.na(val)),
            score = mean(true_val, na.rm=T) ) %>%
  ungroup()


q5yrsel_scales <- q5yrsel_scales %>%
  mutate(sc_score = ifelse(items_present >= (items_scale/2), round(score*items_scale,0), NA))


a <- q5yrsel_scales %>%
  filter(sscale %in% 'ccc2')
table(a$sc_score)
hist(a$sc_score)

b <- q5yrsel_scales %>%
  filter(sscale %in% 'sprk20')
hist(b$sc_score)
table(b$sc_score)

q5yrtlk<- q5yrsel  %>%
  left_join(q5yrsel_scales %>%
              select(PREG_ID_2306,BARN_NR, sscale, sc_score) %>%
              spread(key=sscale, value = sc_score) ) %>%
  select(PREG_ID_2306,BARN_NR, age_1st_wds, age_1st_sns, denom_5yr) 

## Binary motor development

# 18mo

q5bin <- q5 %>%
  select(PREG_ID_2306,BARN_NR, ALDERRETUR_S5,
         EE874, EE840, EE841, EE842,
         EE401, EE986, EE406, EE407,  EE800, EE801, EE802, EE185, EE186)

q5bin <- q5bin %>%
  mutate(bin_mot_18m = ifelse(EE401==1 | EE986=="No" | EE406 == "Not Yet" | EE407 == "Not Yet", 1,0),
         mot_delay_18m = ifelse(EE800 ==1 | EE801 ==1 | EE185 =="Yes", 1,0),
         mot_delay_spec_18m = ifelse(EE802 =="Yes" | EE186 == "Yes", 1,0),
         lang_delay_18m = ifelse(EE840==1 | EE841==1, 1,0),
         lang_delay_spec_18m =ifelse(EE842=="Yes", 1,0) ) %>%
  rename(denom_18m = ALDERRETUR_S5) %>%
  select(PREG_ID_2306,BARN_NR,denom_18m,bin_mot_18m)


## Binary language development

# 36mo

q6bin <- q6 %>%
  select(PREG_ID_2306,BARN_NR,ALDERRETUR_S6,
         GG241, GG242, GG239, GG229, GG286, GG94, GG95, GG96, GG226,
         GG28, GG38, GG39, GG40)
q6bin <- q6bin %>%
  mutate(bin_lang_3yr = ifelse( GG239 == "Not Yet" | GG226 %in% c("Not yet talking", "He/she is talking, but you can not understand him/her", 
                                                                                      "Talking in one-word utterances such as «milk» or «down»", 
                                                                                      "Talking in 2 to 3 word phrases, such as «me got ball» or «give doll»"), 1,0 ),
         mot_delay_3yr = ifelse(GG28 ==1 | GG39 == 1, 1,0),
         mot_delay_spec_3yr = ifelse(GG40 == "Yes", 1,0),
         lang_delay_3yr = ifelse(GG94 == 1 | GG95 == 1, 1,0),
         lang_delay_spec_3yr =ifelse(GG96 == "Yes", 1,0) ) %>%
  rename(denom_3yr = ALDERRETUR_S6) %>%
  select(PREG_ID_2306,BARN_NR,denom_3yr,bin_lang_3yr)


# General concern about development

q6sel2 <- q6 %>%
  select(PREG_ID_2306,BARN_NR, "GG382") %>%
  mutate(concern_dev_exprssd_others = case_when(GG382 == "Yes" ~ 1,
                                                GG382 == "No" ~ 0)) 

 
# Get sex variable from Birth Registry 

sex_mbrn <- mbrn %>%
  select(PREG_ID_2306, BARN_NR, sex = KJONN, SVLEN)


  

 ############################################################
 ############################################################
## Need to restrict to unrelated sample for child-only analysis

incl <- read.table('N:/data/durable/data/genetic/qcd_genetic_data/relatedness_exclusion_flag_list.txt', header= T) %>% 
  filter(children_only_analysis==0)



# Need file that matches ID numbers and PREG_IDs

IDS <- read.table('N:/data/durable/data/Linkage files/core_IDs&covars_hrv_njl_v2.txt', header= T) %>%
  filter(IID %in% incl$IID,
         Role == "Child") %>%
  select(PREG_ID_2306, IID)


# Read in scz PRS - we only need those individuals who have PRSs so can drop_na on this

scz_prs <- read.table('./data/prs_for_analyses/PRS/off_scz_108.prs', header=T) %>% 
  select(IID, matches("resid")) %>%
  rename_at(vars(contains("resid")), funs(paste0("scz_",.)))
asd_prs <- read.table('./data/prs_for_analyses/PRS/off_asd.prs', header=T) %>% 
  select(IID, matches("resid")) %>%
  rename_at(vars(contains("resid")), funs(paste0("asd_",.)))
adhd_prs <- read.table('./data/prs_for_analyses/PRS/off_adhd.prs', header=T) %>% 
  select(IID, matches("resid")) %>%
  rename_at(vars(contains("resid")), funs(paste0("adh_",.)))
# intell_prs <- read.table('./data/prs_for_analyses/PRS/off_intell.prs', header=T) %>% 
#   select(IID, matches("resid")) %>%
#   rename_at(vars(contains("resid")), funs(paste0("iiq_",.)))

rm(q4, q5, q6, q5yrs)

# Join everything together

alldata <- sex_mbrn %>%
  full_join(q5q6wlk) %>%
  full_join(q5yrtlk) %>% 
  full_join(q5bin) %>%
  full_join(q6bin) %>%
  full_join(q6sel2) %>% 
  full_join(IDS) %>%
  full_join(scz_prs) %>%
  full_join(asd_prs) %>%
  full_join(adhd_prs) %>%
  distinct()

head(alldata)

# Assume no delays for anyone who didn't endorse but did complete the questionnaire

alldata <- alldata %>%
  mutate(bin_mot_18m = ifelse(!is.na(denom_18m)& is.na(bin_mot_18m), 0, bin_mot_18m ),
         bin_lang_3yr = ifelse(!is.na(denom_3yr)& is.na(bin_lang_3yr), 0, bin_lang_3yr),
         concern_dev_exprssd_others = ifelse(!is.na(denom_3yr)& is.na(concern_dev_exprssd_others), 0, concern_dev_exprssd_others)) 


#  Drop any hanging empty factor levels

alldata <- alldata %>%
  droplevels()  %>%
  drop_na(sex) %>%
  select(PREG_ID_2306,IID,BARN_NR,sex,denom_18m,denom_3yr,denom_5yr,
         age_1st_wlk_18m,age_1st_wlk_3yr,age_1st_wlk,
         age_1st_wds,age_1st_sns,
         bin_lang_3yr,
         bin_mot_18m,
         concern_dev_exprssd_others,
         matches("resid"))

#These vars
table(alldata$bin_mot_18m)
table(alldata$bin_lang_3yr)
table(alldata$concern_dev_exprssd_others)

names(alldata)

save(alldata, file = './data/prepped_00_vars.RData')


