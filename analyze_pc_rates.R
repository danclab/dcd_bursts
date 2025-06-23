library("lme4")
library("car")
library("lsmeans")
library("Rmisc")
library("ggplot2")

pc<-7
all_data<- read.csv(paste0('./data/derivatives/NEARICA_behav/pc_',pc,'_rates.csv'))

all_data$subject_idx<-as.factor(all_data$subject_idx)
all_data$group<-as.factor(all_data$group)
all_data$condition<-as.factor(all_data$condition)
all_data$pc_bin<-as.factor(all_data$pc_bin)

for(q in 0:3) {
  print(q)
  q_data<-all_data[all_data$pc_bin==q,]
  lmer_model <- lmer(rate ~ group*condition+(1|subject_idx), data = q_data, control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))
  print(lmer_model)
  
  lmer_results<-Anova(lmer_model, type = 3)
  print(lmer_results)
  
  if(lmer_results$`Pr(>Chisq)`[4]<0.05) {
    pw1<-lsmeans(lmer_model, pairwise~group|condition, adjust='tukey')
    print(summary(pw1)$contrasts)
    pw2<-lsmeans(lmer_model, pairwise~condition|group, adjust='tukey')
    print(summary(pw2)$contrasts)
  } else {
    if(lmer_results$`Pr(>Chisq)`[3]<0.05) {
      pw1<-lsmeans(lmer_model, pairwise~condition, adjust='tukey')
      print(summary(pw1)$contrasts)
      conditions<-unique(all_data$condition)
      for(condition in conditions) {
        m_cond<-mean(q_data$rate[q_data$condition==condition])
        sd_cond<-sd(q_data$rate[q_data$condition==condition])
        print(paste0('M',condition,'=',m_cond,', SD=',sd_cond))
      }
    }
    if(lmer_results$`Pr(>Chisq)`[2]<0.05) {
      pw1<-lsmeans(lmer_model, pairwise~group, adjust='tukey')
      print(summary(pw1)$contrasts)
      
      m_dcd<-mean(q_data$rate[q_data$group=='dcd'])
      sd_dcd<-sd(q_data$rate[q_data$group=='dcd'])
      m_typ<-mean(q_data$rate[q_data$group=='typ'])
      sd_typ<-sd(q_data$rate[q_data$group=='typ'])
      print(paste0('Mdcd=',m_dcd,', SD=',sd_dcd,', Mtyp=',m_typ,', SD=',sd_typ))
    }
  }
}


ipsi_data<-read.csv(paste0('./data/derivatives/NEARICA_behav/pc_',pc,'_rates_ipsi.csv'))
ipsi_data$hemi<-'ipsi'
contra_data<-read.csv(paste0('./data/derivatives/NEARICA_behav/pc_',pc,'_rates_contra.csv'))
contra_data$hemi<-'contra'
all_data<-rbind(ipsi_data, contra_data)
all_data$hemi<-as.factor(all_data$hemi)
all_data$subject_idx<-as.factor(all_data$subject_idx)
all_data$group<-as.factor(all_data$group)
all_data$condition<-as.factor(all_data$condition)
all_data$pc_bin<-as.factor(all_data$pc_bin)

for(q in 0:3) {
  print(q)
  q_data<-all_data[all_data$pc_bin==q,]
  lmer_model <- lmer(rate ~ group*hemi+(1+condition|subject_idx), data = q_data, control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))
  #print(lmer_model)
  
  lmer_results<-Anova(lmer_model, type = 3)
  print(lmer_results)
  
  # Three-way
  if(lmer_results$`Pr(>Chisq)`[4]<0.05){
    pw1<-lsmeans(lmer_model, pairwise~group|hemi, adjust='tukey')
    print(summary(pw1)$contrasts)
    pw2<-lsmeans(lmer_model, pairwise~hemi|group, adjust='tukey')
    print(summary(pw2)$contrasts)
    
    for(g in c('dcd','typ')) {
      m_ipsi<-mean(q_data$rate[q_data$group==g & q_data$hemi=='ipsi'])
      sd_ipsi<-sd(q_data$rate[q_data$group==g & q_data$hemi=='ipsi'])
      m_contra<-mean(q_data$rate[q_data$group==g & q_data$hemi=='contra'])
      sd_contra<-sd(q_data$rate[q_data$group==g & q_data$hemi=='contra'])
      print(paste0(g,': Mipsi=',m_ipsi,', SD=',sd_ipsi,', Mcontra=',m_contra,', SD=',sd_contra))
    }
  } else {
    if(lmer_results$`Pr(>Chisq)`[3]<0.05) {
      pw1<-lsmeans(lmer_model, pairwise~hemi, adjust='tukey')
      print(summary(pw1)$contrasts)
      
      m_ipsi<-mean(q_data$rate[q_data$hemi=='ipsi'])
      sd_ipsi<-sd(q_data$rate[q_data$hemi=='ipsi'])
      m_contra<-mean(q_data$rate[q_data$hemi=='contra'])
      sd_contra<-sd(q_data$rate[q_data$hemi=='contra'])
      print(paste0('Mipsi=',m_ipsi,', SD=',sd_ipsi,', Mcontra=',m_contra,', SD=',sd_contra))
    }
    if(lmer_results$`Pr(>Chisq)`[2]<0.05) {
      pw1<-lsmeans(lmer_model, pairwise~group, adjust='tukey')
      print(summary(pw1)$contrasts)
      
      m_dcd<-mean(q_data$rate[q_data$group=='dcd'])
      sd_dcd<-sd(q_data$rate[q_data$group=='dcd'])
      m_typ<-mean(q_data$rate[q_data$group=='typ'])
      sd_typ<-sd(q_data$rate[q_data$group=='typ'])
      print(paste0('Mdcd=',m_dcd,', SD=',sd_dcd,', Mtyp=',m_typ,', SD=',sd_typ))
    }
  }
}