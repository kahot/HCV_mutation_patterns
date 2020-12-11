# Compare the sites that have the same ref nucleotides and mutation types across 
# subtypes ("Conserved Sites")
library(stringr)
library(tidyverse)
library(zoo)
library(reshape2)
library(colorspace)
library(RColorBrewer)
library(sfsmisc)
source("Rscripts/baseRscript.R")

#Geonotype colors 
colors2<-qualitative_hcl(6, palette="Dark3")
# 1A=1, 1B=3, 3A=5, c("#E16A86","#50A315","#009ADE")
col2_light<-qualitative_hcl(6, palette="Set3")
div.colors<-c(colors2[1],col2_light[1],colors2[3],col2_light[3],colors2[5],col2_light[5])

subs<-c("1A","1B","3A")


###read the MF data that are the same as ref to be consistent. 
sameDF<-read.csv("Output_all/MVF/Ts.Same_Mean_3Subtypes.csv", stringsAsFactors = F, row.names = 1)
sameDF<-sameDF[sameDF$merged.pos<8649& sameDF$merged.pos>=264,] 
sameDF$gene<-factor(sameDF$gene, levels=c("5' UTR","Core","E1", "HVR1","E2","P7", "NS2","NS3","NS4A","NS4B","NS5A","NS5B"))


#attach the depth info:
for (g in 1:3){
        gty<-subs[g]
        dt1<-read.csv(paste0("Output_all/MVF/Total_Reads.",gty,".csv"), stringsAsFactors = F, row.names = 1)
        dt1$Depth<-rowSums(dt1[2:ncol(dt1)], na.rm=T)
        depth<-dt1[,c("merged.pos","Depth")]
        colnames(depth)[2]<-paste0("depth.",gty)
        sameDF<-merge(sameDF, depth, by="merged.pos")
}

########
## With depth info, create "conserved-site summary"

cons<-sameDF[sameDF$ref.1A==sameDF$ref.1B & sameDF$ref.1A==sameDF$ref.3A & sameDF$ref.1B==sameDF$ref.3A,] 
cons<-cons[cons$Type.1A==cons$Type.1B & cons$Type.1A==cons$Type.3A & cons$Type.1B==cons$Type.3A|cons$gene=="5' UTR",] 
cons<-cons[!is.na(cons$merged.pos),] #5003 sites

cons5UTR<-cons[cons$gene=="5' UTR",] 
cons<-cons[cons$gene!="5' UTR",]

#Reformat data for analysis 
mf<-data.frame()
for(g in 1:3){
        dat<-cons[,c("merged.pos","gene",paste0("mean.",subs[g]), paste0("Type.",subs[g]), paste0("depth.",subs[g]))]
        colnames(dat)[3:5]<-c("mean","Type","depth")
        dat$Subtype<-subs[g]
        mf<-rbind(mf,dat)
}
mf1<-mf[!is.na(mf$mean),]
mf1<-mf1[mf1$Type!="stop",]

genenames<-levels(sameDF$gene)
genenames<-genenames[2:12]

#Summary table
#calculate means
SumMFGenes<-aggregate(mf1$mean,by=list(mf1$Subtype,mf1$gene, mf1$Type),FUN=mean)
#calculate SE
meSE<-SumMFGenes[,1:3]
i=1
for (t in c("nonsyn", "syn")){
        for (k in 1:length(genenames)){
                for (g in 1:3){
                        df<-mf1[mf1$gene==genenames[k]&mf1$Subtype==subs[g]&mf1$Type==t,]
                        meSE$SE[i]<-sqrt(mean(df$mean)*(1-mean(df$mean))/mean(df$depth, na.rm=T))
                        i=i+1
                }
        }
}
sumGsame<-cbind(SumMFGenes, meSE$SE)
colnames(sumGsame)<-c("Subtype","Gene","Type","Mean","SE")
sumGsame$Gene<-factor(sumGsame$Gene, levels=c("Core","E1", "HVR1","E2","P7", "NS2","NS3","NS4A","NS4B","NS5A","NS5B"))


#Compare the sites that had the same type and nucleotide vs. all sites across 3 subtpyes:
sumGsame$group<-"Conserved"
sameDF2<-sameDF[sameDF$gene!="5' UTR",]
varied<-sameDF2[!(sameDF2$merged.pos %in% cons$merged.pos),] 

mf2<-data.frame()
for(g in 1:3){
        dat<-varied[,c("merged.pos","gene",paste0("mean.",subs[g]), paste0("Type.",subs[g]), paste0("depth.",subs[g]))]
        colnames(dat)[3:5]<-c("mean","Type","depth")
        dat$Subtype<-subs[g]
        mf2<-rbind(mf2,dat)
}

#rmove NA and stop sites
mf2<-mf2[!is.na(mf2$mean),]
mf2<-mf2[mf2$Type!="stop",]

SumMFGenes2<-aggregate(mf2$mean,by=list(mf2$Subtype,mf2$gene, mf2$Type),FUN=mean)

meSE2<-SumMFGenes2[,1:3]
i=1
for (t in c("nonsyn", "syn")){
        for (k in 1:length(genenames)){
                for (g in 1:3){
                        df<-mf2[mf2$gene==genenames[k]&mf2$Subtype==subs[g]&mf2$Type==t,]
                        df<-df[!is.na(df$mean),]
                        meSE2$SE[i]<-sqrt(mean(df$mean, na.rm=T)*(1-mean(df$mean, na.rm=T))/mean(df$depth, na.rm=T))
                        i=i+1
                }
        }
}


sumGvary<-cbind(SumMFGenes2, meSE2$SE)
colnames(sumGvary)<-c("Subtype","Gene","Type","Mean","SE")
sumGsame$Gene<-factor(sumGsame$Gene, levels=c("Core","E1", "HVR1","E2","P7", "NS2","NS3","NS4A","NS4B","NS5A","NS5B"))

sumGvary$group<-"non-conserved"
        
sumGAll<-rbind(sumGvary,sumGsame)

### Plot the results together
co<-rep(div.colors[c(2,4,6)], times=24)

labs<-list(
        "nonsyn"= "Nonsyn",
        'syn'= "Syn")


plot_labeller <- function(variable,value){
        return(labs[value])
}



ggplot(sumGAll, aes(x=Gene, y=Mean, shape=group, color=Subtype, fill=Subtype))+
        geom_point(aes(group=factor(Subtype), color=factor(Subtype)), position=position_dodge(width=0.8),size =2)+
        scale_color_manual(values=colors2[c(1,3,5)])+
        scale_fill_manual(values=paste0(colors2[c(1,3,5)],"66"), guide = 'none')+
        scale_shape_manual(values=c(16,25))+
        geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE,group=factor(Subtype), color=factor(Subtype)), width=.2, size=.2, position=position_dodge(width=0.8))+
        theme(axis.title.x=element_blank())+ylab("Minor variant frequency")+
        theme_bw()+
        labs(x="")+
        theme(axis.text.x = element_text(size=12),axis.title.y = element_text(size=13))+
        theme(panel.grid.major.x = element_blank())+
        geom_vline(xintercept = c(1:10)+0.5,  
                   color = "gray60", size=.4)+
        guides(
                shape = guide_legend(order = 2),
                color=guide_legend(order=1)
        )+
        labs(shape="")+
        facet_grid(Type~., scales="free",space="free", labeller=plot_labeller)
ggsave(filename="Output_all/Figures/MVF_byGene.Consvd.vs.non-consvd.pdf", width = 8, height = 6.5)

  
## stats: Test mutation frequencies at conserved vs. non-conserved sites (syn and nonsyn separately)
#mf1 == conserved, mf2= non-conserved sites

wilcoxresG.s<-data.frame(matrix(nrow=11, ncol=3))
rownames(wilcoxresG.s)<-genenames
colnames(wilcoxresG.s)<- subs                    

wilcoxresG.ns<-data.frame(matrix(nrow=11, ncol=3))
rownames(wilcoxresG.ns)<-genenames
colnames(wilcoxresG.ns)<- subs    

for (i in 1:3){
    for (g in 1:11){
        s<-wilcox.test(mf1$mean[mf1$gene==genenames[g]&mf1$Subtype==subs[i]&mf1$Type=="syn"], 
                       mf2$mean[mf2$gene==genenames[g]& mf2$Subtype==subs[i]&mf2$Type=="syn"], alternative = "less", paired = FALSE) 
        wilcoxresG.s[g,i]<-s[[3]][1]
        ns<-wilcox.test(mf1$mean[mf1$gene==genenames[g]& mf1$Subtype==subs[i]&mf1$Type=="nonsyn"], mf2$mean[mf2$gene==genenames[g]& mf2$Subtype==subs[i]&mf2$Type=="nonsyn"], alternative = "less", paired = FALSE) 
        wilcoxresG.ns[g,i]<-ns[[3]][1]
    }
}

write.csv(wilcoxresG.s, "Output_all/MVF/WilcoxTest_Conserved.vs.nonConserved.syn.csv",row.names = T)
write.csv(wilcoxresG.ns, "Output_all/MVF/WilcoxTest_Conserved.vs.nonConserved.nonsyn.csv",row.names = T)
