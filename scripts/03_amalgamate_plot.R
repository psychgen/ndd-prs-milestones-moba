load( file = './scratch_data/continuous_analyses_estimates.RData')


load( file = './scratch_data/dichotomous_analyses_estimates.RData')


load( file = './data/prepped_00_vars.RData')

install.packages("tidyverse", repos = "file://tsd-evs/shared/R/cran")

library(tidyverse)
#library(Hmisc)

#Plots for all ests across all thresholds

all_ests <- ests1_plt %>%
  ungroup() %>% 
  bind_rows(ests3_plt %>%  ungroup()) 
all_ests <- all_ests %>%
  mutate_at(vars(score,var,pheno, model, sex), funs(factor(.))) %>%
  mutate(mtc_sig = factor(ifelse(pvalue<(0.05/21), 1, 0)),
         pheno = str_replace(pheno, "ASD", "Autism"))


# MOTOR plot

motplotdat_all <- all_ests %>%
  filter(var %in%c("Age at first \nwalking","Motor delays\nat 18 months")) %>%
  droplevels()%>%
  mutate(var = factor(var, levels=c("Age at first \nwalking","Motor delays\nat 18 months") ) ) %>% 
  as.data.frame()


# LANG plot

langplotdat_all <- all_ests %>%
  filter(var %in% c("Age at first \nwords","Age at first \nsentences", "Rate (age first \nsentences - age \nfirst words)","Language delays\nat 3 years")) %>%
  droplevels() %>%
  mutate(var = factor(var, levels=c("Age at first \nwords","Age at first \nsentences", "Rate (age first \nsentences - age \nfirst words)","Language delays\nat 3 years") ) ) %>% 
  as.data.frame()



# concern plot

concplotdat_all <- all_ests %>%
  filter(var %in%c("Concern regarding \ndevelopment expressed \nby others at 3 years")) %>%
  droplevels() %>%
  as.data.frame()



plotdats <- list(motplotdat_all,
                 langplotdat_all,
                 concplotdat_all)
modnames <- c("motor_all", 
              "language_all",
              "concern_all")


pd <- position_dodge(width = 1)

for(i in 1:length(plotdats)){
  
  mtc_sig_dat <- as.data.frame(plotdats[[i]]) %>%
    select(pheno, var, mtc_sig) %>%
    group_by(pheno,var) %>%
    distinct()
  
  tiff(paste0("figures/multigroup_analyses_", modnames[[i]],"_bestmodel.tiff"), res = 600, compression = "lzw", unit = "in",
       height = 10, width =8)
  
  print(ggplot(plotdats[[i]],aes(x= score, y = Beta, shape= mtc_sig)) +
          geom_point(aes( size=factor(mysize), ymin = lci, ymax = uci,  alpha = myalpha,  group = interaction(score,sex), colour =sex),stroke=0, position = pd) +
          geom_errorbar(aes(colour=sex, ymin = lci, ymax = uci,  alpha = myalpha,  group = interaction(score,sex)),position = pd , size = 0.6,width = 0) +
         # geom_text(data=mtc_sig_dat, aes(label=mtc_sig, x=3.7, y=-0.15), cex =8, colour= "grey20")+
          facet_grid( var ~ pheno , scales = "fixed") +
          scale_shape_manual(values=c( 19, 17), name = "Passes multiple\n testing correction", breaks=c(0,1), labels=c("No", "Yes") )+
          scale_size_discrete(range = c(1.6,5), name = "PGS at threshold:", labels=levels(motplotdat_all$score))+
          scale_alpha_continuous(range = c(0.37,1)) +
          scale_colour_brewer(type = "qual", palette = 2, direction =1, name="Sex")+
          scale_fill_brewer(type = "qual", palette = 2, direction =1, name="Sex")+
          geom_hline( aes(yintercept=0), color = "grey70", linetype = 2, size=0.8)+
          theme(axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                axis.title.x = element_blank(),
                axis.title.y = element_text(margin = margin(t = 0, r = 30, b = , l = 0)),
                text =element_text(size = 14.5),
                panel.background = element_rect(fill = "white", colour = "grey90"),
                panel.spacing = unit(0.2, "lines"),
                panel.grid.major = element_line(colour="white"),
                strip.text = element_text(size=11, face = "bold" ),
                strip.text.y = element_text(angle = 0),
                legend.title = element_text(size =11),
                strip.background = element_rect(fill = "grey90", colour = "white")) +
          coord_cartesian(ylim = c(-0.14,0.20))+
          labs(shape="Best-fitting\nmodel")+ scale_y_continuous(breaks = c(-0.10,-0.05,0,0.05,0.10,0.15,0.20)) +
          guides(alpha=FALSE))
        
    
  dev.off()

}


## Now plots for max_rsq models only


all_maxrsq_ests <- ests1_plt_maxrsq %>%
  ungroup() %>% 
  bind_rows(ests3_plt_maxrsq %>% ungroup()) %>% 
  ungroup()
all_maxrsq_ests <- all_maxrsq_ests %>%
  mutate_at(vars(score,var,pheno, model, sex), funs(factor(.))) %>%
  mutate(mtc_sig = factor(ifelse(pvalue<(0.05/21), 1, 0))) %>% 
  mutate(domain = factor(case_when(var %in% c( "Age at first \nwalking","Motor delays\nat 18 months")~"Motor",
                                   var %in% c("Age at first \nwords", "Age at first \nsentences", "Rate (age first \nsentences - age \nfirst words)","Language delays\nat 3 years")~"Language",
                                   var %in% c( "Concern regarding \ndevelopment expressed \nby others at 3 years")~"General")),
         var = factor(var,levels= c("Age at first \nwalking","Motor delays\nat 18 months",
                                "Age at first \nwords", "Age at first \nsentences", "Rate (age first \nsentences - age \nfirst words)","Language delays\nat 3 years",
                                "Concern regarding \ndevelopment expressed \nby others at 3 years") ),
         pheno = str_replace(pheno, "ASD", "Autism"))
pd <- position_dodge(width = 0.4)


all_labs<- all_maxrsq_ests %>% 
  select(pheno,domain,var,score)

p1<-ggplot(all_maxrsq_ests,aes(x= fct_rev(var), y = Beta)) +
  geom_hline( aes(yintercept=0), color = "grey70", linetype = 2, size=0.8)+
  geom_point(aes(  alpha = myalpha,  group = interaction(score,sex), colour =sex),shape= 19,size=4,stroke=0, position = pd) +
  geom_errorbar(aes(colour=sex, ymin = lci, ymax = uci,  alpha = myalpha-0.2,  group = interaction(score,sex)),position = pd , size = 1.4,width = 0) +
  #geom_label(data=all_labs,aes(x=fct_rev(var),y=0.025,label= score),size=3, colour="black", fill="white", nudge_x = -0.5)+
  facet_grid( fct_rev(domain) ~ pheno , scales = "free", space="free") +
  coord_flip()+
  scale_shape_manual(values=c( 19, 17), name = "Passes multiple\n testing correction", breaks=c(0,1), labels=c("No", "Yes") )+
  scale_alpha_continuous(range = c(0.37,1)) +
  scale_colour_brewer(type = "qual", palette = 2, direction =1, name="Sex")+
  scale_fill_brewer(type = "qual", palette = 2, direction =1, name="Sex")+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 0, b = , l = 0)),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = , l = 0)),
        text =element_text(size = 14.5),
        axis.text.x  = element_text(angle=90, vjust=0.5),
        panel.background = element_rect(fill = "white", colour = "grey30"),
        panel.spacing = unit(1, "lines"),
        panel.grid.major = element_line(colour="grey90"),
        strip.text = element_text(size=11, face = "bold", colour= "white" ),
        strip.text.y = element_text(angle = 0, colour= "white"),
        legend.title = element_text(size =11),
        strip.background = element_rect(fill = "grey30", colour = "white")) +
  labs(shape="Best-fitting\nmodel")+ 
  scale_y_continuous("Standardized beta", breaks = c(-0.10,-0.05,0,0.05,0.10,0.15,0.20)) +
  scale_x_discrete("Measure") +
  guides(alpha=FALSE, colour = guide_legend(reverse = T))


tiff(paste0("figures/multigroup_allmodels_bestrsq.tiff"), res = 600, compression = "lzw", unit = "in",
     height = 8, width =12)

p1

dev.off()
# Probit regression betas are not informative in and of themselves,
# so we want to plot predicted probabilities with CIs from best-fitting models

# Predicted probabilities for different PGS


xs <- c("scz_X0.001_resid", "scz_X0.01_resid", "scz_X0.05_resid", "scz_X0.1_resid", 
        "scz_X0.2_resid", "scz_X0.5_resid", "scz_X1_resid",  "asd_X0.001_resid", 
        "asd_X0.01_resid", "asd_X0.05_resid", "asd_X0.1_resid", "asd_X0.2_resid", 
        "asd_X0.5_resid", "asd_X1_resid",   
        "adh_X0.001_resid", "adh_X0.01_resid", "adh_X0.05_resid", 
        "adh_X0.1_resid", "adh_X0.2_resid", "adh_X0.5_resid", "adh_X1_resid")
ys <- c( "bin_mot_18m","bin_lang_3yr", "concern_dev_exprssd_others")


params <- expand.grid(xs,ys, stringsAsFactors = FALSE)

load( file = './scratch_data/dichotomous_analyses_estimates.RData')
prprobs2 <- prprobs %>%
  mutate(score = rep(params$Var1, each = dim(prprobs)[1]/dim(params)[1]),
         var = rep(params$Var2, each = dim(prprobs)[1]/dim(params)[1]),
         sex = recode(group, `1`="Male",`2`= "Female"),
         i_est = 1-est,
         i_uci = 1-uci,
         i_lci = 1-lci,
         r_prsval = prsval*-1) %>%
  right_join(fitcomp_best %>% select(score, var,model))


exemplar <- prprobs2 %>%
  filter(score =="adh_X0.01_resid",
         var == "bin_lang_3yr")                              

tiff(paste0("figures/multigroup_predprobs_lanfdelay001.tiff"), res = 600, compression = "lzw", unit = "in",
     height = 6, width =6)

ggplot(data= exemplar,aes(x=r_prsval,y=i_est))+
  geom_line(aes(colour=sex), size=0.9) +
  geom_ribbon(aes(ymin=i_uci,ymax=i_lci, fill=sex), alpha=0.15) +
  scale_color_manual("", values = c( "#7570B3", "#1B9E77" ))+
  scale_fill_manual("",values = c("#7570B3", "#1B9E77" ))+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 30, b = , l = 0)),
        text =element_text(size = 14.5),
        panel.background = element_rect(fill = "white", colour = "grey80"),
        panel.spacing = unit(0.2, "lines"),
        panel.grid.major = element_line(colour="white"),
        strip.text = element_text(size=11),
        legend.title = element_text(size =11)) +
  scale_x_continuous(name = "ADHD PGS")+
  scale_y_continuous(name = "Model-predicted probability of \nlanguage delays at 3yrs")

dev.off()


##

# Amalgamate fit tables

load( file = './scratch_data/dichotomous_analyses_fitcomps.RData')

rsqs_bysex_dich <- rsqs_bysex

load( file = './scratch_data/continuous_analyses_fitcomps.RData')


rsqs_bysex_all<-rsqs_bysex %>% 
  select(pheno,score,group,var,model,rsq=estimate) %>% 
  bind_rows(rsqs_bysex_dich %>% 
              select(pheno,score,group,var,model,rsq=estimate)) %>% 
  mutate(pheno = factor(str_sub(score, 0,3), levels = c("adh","asd","scz"),labels = c("ADHD","Autism","SCZ")),
         score = factor(str_sub(score, start = 6), levels=c("0.001_resid", "0.01_resid", "0.05_resid", "0.1_resid", 
                      "0.2_resid", "0.5_resid", "1_resid"), labels=c("p<0.001", "p<0.01", 
                                                                                                           "p<0.05","p<0.1","p<0.2","p<0.5","All SNPs")),
         var = factor(var, levels = c("bin_mot_18m","age_1st_wlk", "age_1st_wds", "age_1st_sns", "diff_1st_wds_sns","bin_lang_3yr",  "concern_dev_exprssd_others"), 
                      labels =c("Motor delays\nat 18 months", "Age at first \nwalking", "Age at first \nwords", "Age at first \nsentences", "Rate (age first \nsentences - age \nfirst words)",
                                "Language delays\nat 3 years","Concern regarding \ndevelopment expressed \nby others at 3 years") ))%>% 
  select(pheno,group,var,model,score,rsq) 




