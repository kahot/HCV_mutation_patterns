#GLM/Beta regression preparation:
library(tidyverse)
library(zoo)
library(purrr)
library(MASS)
library(betareg)
library(miscTools)
library(colorspace)
source("Rscripts/baseRscript.R")

####################

#Read the Betareg formatted data
brdata2<-list()
subs<-c("1A", "1B","3A")
for (i in 1:length(subs)){
        df<-read.csv(paste0("Output_all/Betareg/brData_withGenes.",subs[i],".csv"),stringsAsFactors = F, row.names = 1)
        brdata2[[i]]<-df
        names(brdata2)[i]<-subs[i]
}

################
#run Beta Regression for each dataset
#1.  Subtype 1A
brData <-brdata2[[1]]
brData2<-brData[complete.cases(brData),]
#remove the stop sites
brData2<-brData2[brData2$Stop == 0,]

# run the best models (mod1 and mod2) from previous studies
mod1<-betareg(mean ~ t + c + g + CpG + Nonsyn +bigAAChange+ t:Nonsyn + c:Nonsyn + g:Nonsyn +
              Core +E1 +HVR1 +E2 +NS1 +NS2+NS3 +NS4A+NS5A+NS5B, data = brData2)
summary(mod1)

#The best fit model from Fitness Cost Analysis (remove two genes):  
mod2 <- update(mod1, ~. -NS3 - NS5A)
summary(mod2)
#            Estimate Std. Error  z value Pr(>|z|)    
#(Intercept) -4.03696    0.02362 -170.927  < 2e-16 ***
#t            0.14808    0.02662    5.562 2.66e-08 ***
#c           -0.52754    0.02626  -20.093  < 2e-16 ***
#g           -0.72692    0.02932  -24.790  < 2e-16 ***
#CpG         -0.08886    0.01791   -4.961 7.00e-07 ***
#Nonsyn      -0.95365    0.02945  -32.387  < 2e-16 ***
#bigAAChange -0.16704    0.01882   -8.878  < 2e-16 ***
#Core        -0.32350    0.02635  -12.279  < 2e-16 ***
#E1           0.01651    0.02382    0.693   0.4882    
#HVR1         0.27903    0.05915    4.717 2.39e-06 ***
#E2           0.09398    0.01842    5.103 3.34e-07 ***
#NS1          0.04810    0.03950    1.218   0.2233    
#NS2          0.10187    0.02208    4.613 3.96e-06 ***
#NS4A        -0.06578    0.04494   -1.464   0.1433    
#NS5B        -0.21094    0.02043  -10.326  < 2e-16 ***
#t:Nonsyn    -0.24249    0.03537   -6.855 7.13e-12 ***
#c:Nonsyn    -0.08226    0.03603   -2.283   0.0224 *  
#g:Nonsyn     0.17275    0.03714    4.651 3.30e-06 ***


#remove other ns genes
mod3 <- update(mod2, ~. -E1 -NS4A)
summary(mod3)

mod4 <- update(mod3, ~. -NS1)
summary(mod4)


AIC(mod1,mod2,mod3,mod4)
#     df       AIC
#mod1 21 -67783.90
#mod2 19 -67787.68
#mod3 17 -67788.82
#mod4 16 -67789.36

sum.g1<-summary(mod1)
sum.g2<-summary(mod2)
sum.g3<-summary(mod3)

modcoef1<-sum.g1$coefficients
coef1<-modcoef1[[1]]
write.csv(coef1,"Output_all/Betareg/BetaReg_mod1_1A.csv")
modcoef2<-sum.g2$coefficients
coef2<-modcoef2[[1]]
write.csv(coef2,"Output_all/Betareg/BetaReg_mod2_1A.csv")
modcoef3<-sum.g3$coefficients
coef3<-modcoef3[[1]]
write.csv(coef3,"Output_all/Betareg/BetaReg_mod3_1A.csv")

#######################################################################################################

# 2.  Subtype 1B
brData <-brdata2[[2]]
brData2<-brData[complete.cases(brData),] 
brData2<-brData2[brData2$Stop == 0,] #7988

# run the models
mod1<-betareg(mean ~ t + c + g + CpG + Nonsyn +bigAAChange+ t:Nonsyn + c:Nonsyn + g:Nonsyn +
                      Core +E1 +HVR1 +E2 +NS1 +NS2+NS3 +NS4A+NS5A+NS5B,data = brData2)
summary(mod1)

# mod2:  
mod2 <- update(mod1, ~. -NS3 - NS5A)
summary(mod2)

#remove other ns genes
mod3<- update(mod2, ~. -E1 -NS4A)
summary(mod3)

mod4 <- update(mod3, ~. -NS1)
summary(mod4)


AIC(mod1,mod2,mod3,mod4)
#     df       AIC
#mod1 21 -67369.70
#mod2 19 -67371.86
#mod3 17 -67374.34
#mod4 16 -67362.34

sum.g1<-summary(mod1)
sum.g2<-summary(mod2)
sum.g3<-summary(mod3)

modcoef1<-sum.g1$coefficients
coef1<-modcoef1[[1]]
write.csv(coef1,"Output_all/Betareg/BetaReg_mod1_1B.csv")
modcoef2<-sum.g2$coefficients
coef2<-modcoef2[[1]]
write.csv(coef2,"Output_all/Betareg/BetaReg_mod2_1B.csv")
modcoef3<-sum.g3$coefficients
coef3<-modcoef3[[1]]
write.csv(coef3,"Output_all/Betareg/BetaReg_mod3_1B.csv")


########################################################################


#3.Subtype 3A
brData<-brdata2[[3]]
brData2<-brData[complete.cases(brData),] 
brData2<-brData2[brData2$Stop == 0,] #8196

mod1 <- betareg(mean ~ t + c + g + CpG  + t:Nonsyn + c:Nonsyn + g:Nonsyn + Nonsyn +bigAAChange + 
                          Core +E1 +HVR1++E2 +NS1 +NS2++NS3 +NS4A+NS5A+NS5B, data = brData2)
summary(mod1)
        
#mod2
mod2 <- update(mod1, ~. -NS3 -NS5A)
summary(mod2)
        
#mod3
mod3 <- update(mod2, ~. -E1 -NS4A)
#mod4
mod4 <- update(mod3, ~. -NS1)


AIC(mod1,mod2,mod3, mod4)
#       df       AIC
#mod1 21 -68407.58
#mod2 19 -68402.88
#mod3 17 -68406.15
#mod4 16 -68404.62

sum.g1<-summary(mod1)
sum.g2<-summary(mod2)
sum.g3<-summary(mod3)

modcoef1<-sum.g1$coefficients
coef1<-modcoef1[[1]]
write.csv(coef1,"Output_all/Betareg/BetaReg_mod1_3A.csv")
modcoef2<-sum.g2$coefficients
coef2<-modcoef2[[1]]
write.csv(coef2,"Output_all/Betareg/BetaReg_mod2_3A.csv")
modcoef3<-sum.g3$coefficients
coef3<-modcoef3[[1]]
write.csv(coef3,"Output_all/Betareg/BetaReg_mod3_3A.csv")



####################
####################
#Calculate the effect sizes:
source("Rscripts/BetaEffectSize.R")

#Use mod3
effects<-list()
subs<-c("1A", "1B","3A")
for (i in 1:length(subs)){
        df<-read.csv(paste0("Output_all/Betareg/BetaReg_mod3_",subs[i],".csv"), stringsAsFactors = F, row.names = 1)
        effect<-BetaEffectSizeDF(df)
        write.csv(effect, paste0("Output_all/Betareg/Effectsize_",subs[i],".csv"))
        effects[[i]]<-effect
        names(effects)[i]<-subs[i]
}



###################
## plot the effects of each component for the 3 subotypes
ef.df<-list()
ef<-data.frame(factor=rownames(effects[[1]]))
for (i in 1:3){
        df<-effects[[i]]
        df2<-ef
        df2$Effect<-as.numeric(df$Effect)*100
        df2$Subtype=subs[i]
        ef.df[[i]]<-df2
}
BReffects<-do.call(rbind, ef.df)
BReffects<-BReffects[BReffects$factor!="(Intercept)",]
write.csv(BReffects,"Output_all/Betareg/EffectSizes_mod3_all.csv" )

colors2<-qualitative_hcl(6, palette="Dark3")
BReffects$factor<-factor(BReffects$factor, levels=BReffects$factor[1:15])

BReffects$Subtype[BReffects$Subtype=="1A"]<-"1a"
BReffects$Subtype[BReffects$Subtype=="1B"]<-"1b"
BReffects$Subtype[BReffects$Subtype=="3A"]<-"3a"




ggplot(BReffects, aes(factor(factor), Effect, color=Subtype, fill =Subtype)) +
        geom_bar(stat="identity", position="dodge")+
        scale_fill_manual(values=paste0(colors2[c(1,3,5)],"CC"))+
        scale_color_manual(values=colors2[c(1,3,5)])+
        theme_test() +
        theme(axis.text=element_text(size=13),
               axis.title=element_text(size=14,))+
        theme(axis.text.x = element_text(angle =45, hjust = 1, color=1))+
        labs(x="", y="Estimated effects (%)")+
        theme(panel.grid.major.y = element_line(color="gray90"))

ggsave("Output_all/Figures/Betareg_Effectsizes_Mod3.pdf",width = 9, height = 5.5)



####
#Test the differences in mut freq patterns between the subtypes
library(MASS)
library(DescTools)

eff1A<-effects[[1]]
estMF<-eff1A[,4:5]
colnames(estMF)[2]<-"1A"
eff1B<-effects[[2]]
eff3A<-effects[[3]]

estMF$`3A`<-as.numeric(eff3A$estMF)
estMF<-estMF[,-1]

estMF$`1A`<-as.numeric(estMF$`1A`)
DF<-estMF
#DF<-data.frame(sapply(estMF, function(x) as.numeric(as.character(x))))
#colnames(DF)<-c("1A","3A")
#DF<-lapply(estMF, function(x){ as.numeric(as.character(x))} )

#chisquare test
#1. 1A vs.3A
chisq.test(DF)
#X-squared = 0.00102, df = 15, p-value = 1

#2. 1A vs. 1B
estMF$`1B`<-as.numeric(eff1B$estMF)
DF<-estMF[,c("1A","1B")]
chisq.test(DF)
#X-squared = 0.0007473, df = 15, p-value = 1

#2. 1B vs. 3A
DF<-estMF[,c("1B","3A")]
chisq.test(DF)
#X-squared = 0.0006021, df = 15, p-value = 1

