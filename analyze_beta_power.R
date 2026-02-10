library("lme4")
library("car")
library("lsmeans")
library("Rmisc")
library("ggplot2")

all_data<- read.csv(paste0('./trial_beta_power_windows.csv'))
all_data$subject<-as.factor(all_data$subject)
all_data$group<-as.factor(all_data$group)
all_data$hemi<-as.factor(all_data$hemi)
all_data$condition<-as.factor(all_data$condition)
all_data$window<-as.factor(all_data$window)

contra_exe_f_data<-all_data[all_data$hemi=='contra' & all_data$condition=='exe fine',]
lmer_model <- lmer(power ~ group*window+(1|subject), data = contra_exe_f_data, control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))
print(lmer_model)

lmer_results<-Anova(lmer_model, type = 3)
print(lmer_results)

pw1<-lsmeans(lmer_model, pairwise~group|window, adjust='tukey')
print(summary(pw1)$contrasts)
pw2<-lsmeans(lmer_model, pairwise~window|group, adjust='tukey')
print(summary(pw2)$contrasts)


ipsi_exe_f_data<-all_data[all_data$hemi=='ipsi' & all_data$condition=='exe fine',]
lmer_model <- lmer(power ~ group*window+(1|subject), data = ipsi_exe_f_data, control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))
print(lmer_model)
  
lmer_results<-Anova(lmer_model, type = 3)
print(lmer_results)
  
pw1<-lsmeans(lmer_model, pairwise~group|window, adjust='tukey')
print(summary(pw1)$contrasts)
pw2<-lsmeans(lmer_model, pairwise~window|group, adjust='tukey')
print(summary(pw2)$contrasts)


contra_exe_g_data<-all_data[all_data$hemi=='contra' & all_data$condition=='exe gross',]
lmer_model <- lmer(power ~ group*window+(1|subject), data = contra_exe_g_data, control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))
print(lmer_model)

lmer_results<-Anova(lmer_model, type = 3)
print(lmer_results)

pw1<-lsmeans(lmer_model, pairwise~window, adjust='tukey')
print(summary(pw1)$contrasts)


ipsi_exe_g_data<-all_data[all_data$hemi=='ipsi' & all_data$condition=='exe gross',]
lmer_model <- lmer(power ~ group*window+(1|subject), data = ipsi_exe_g_data, control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))
print(lmer_model)

lmer_results<-Anova(lmer_model, type = 3)
print(lmer_results)

pw1<-lsmeans(lmer_model, pairwise~group|window, adjust='tukey')
print(summary(pw1)$contrasts)
pw2<-lsmeans(lmer_model, pairwise~window|group, adjust='tukey')
print(summary(pw2)$contrasts)

