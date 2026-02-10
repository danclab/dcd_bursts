library("lme4")
library("car")
library("lsmeans")
library("Rmisc")
library("ggplot2")

for(pc in c(6,8)) {
  all_data<- read.csv(paste0('./trial_pc',pc,'_burst_rate_windows.csv'))
  all_data$subject<-as.factor(all_data$subject)
  all_data$group<-as.factor(all_data$group)
  all_data$hemi<-as.factor(all_data$hemi)
  all_data$condition<-as.factor(all_data$condition)
  all_data$tertile<-as.factor(all_data$tertile)
  
  for(cond in c('exe fine', 'exe gross')) {
    for(hemi in c('contra','ipsi')) {
      hemi_cond_data<-all_data[all_data$hemi==hemi & all_data$condition==cond,]
      for(t in 0:2){
        print(paste0('PC',pc,': ',hemi, ' ', cond, ' Tertile ',t+1))
        tert_data<-hemi_cond_data[hemi_cond_data$tertile==t,]
        lmer_model <- lmer(rate ~ group+(1|subject), data = tert_data, control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))
        print(lmer_model)
        
        lmer_results<-Anova(lmer_model, type = 3)
        print(lmer_results)
      }
    }
  }
}