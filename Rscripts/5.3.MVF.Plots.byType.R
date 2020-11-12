### Plot Transition MVF separately for Nonsyn and syn
library(zoo)
library(reshape2)
library(colorspace)
#library(sfsmisc)
library(grid)
source("Rscripts/baseRscript.R")

#Subtype colors 
colors2<-qualitative_hcl(6, palette="Dark3")
subs<-c("1A","1B","3A")

##################  Use Transiton MVF with only sites that are maj==ref for plotting by mutation types  #####

sameDF<-read.csv("Output_all/MVF/Ts.Same_Mean_3subtypes.csv", stringsAsFactors = F, row.names = 1)
sameDF$gene[sameDF$gene=="NS1(P7)"]<-"NS1"

for (g in 1:3){
        sub<-subs[g]
        dt1<-read.csv(paste0("Output_all/MVF/Total_Reads.",sub,".csv"), stringsAsFactors = F, row.names = 1)
        dt1$Depth<-rowSums(dt1[2:ncol(dt1)], na.rm=T)
        depth<-dt1[,c("merged.pos","Depth")]
        colnames(depth)[2]<-paste0("depth.",sub)
        sameDF<-merge(sameDF, depth, by="merged.pos")
}

mf3<-data.frame()
for(g in 1:3){
        dat<-sameDF[,c("merged.pos","gene",paste0("mean.",subs[g]),paste0("Type.",subs[g]),paste0("depth.",subs[g]))]
        colnames(dat)[3:5]<-c("mean", "Type", "depth")
        dat$Subtype<-subs[g]
        mf3<-rbind(mf3,dat)
}

mf3<-mf3[!is.na(mf3$mean),]
#remove 5'UTR
mf3<-mf3[mf3$gene!="5' UTR",]
SumMFGenes3<-aggregate(mf3$mean,by=list(mf3$Subtype,mf3$gene, mf3$Type),FUN=mean)
genenames<-unique(SumMFGenes3$Group.2)

meSE3<-SumMFGenes3[,1:3]
i=1
for (t in c("nonsyn", "syn", "stop")){
        for (k in 1:length(genenames)){
                for (g in 1:3){
                        df<-mf3[mf3$gene==genenames[k]&mf3$Subtype==subs[g]&mf3$Type==t,]
                        df[is.na(df$mean),]
                        meSE3$SE[i]<-sqrt(mean(df$mean, na.rm=T)*(1-mean(df$mean, na.rm=T))/mean(df$depth, na.rm=T))
                        i=i+1
                }
        }
}

sumG3<-cbind(SumMFGenes3, meSE3$SE)
colnames(sumG3)<-c("Subtype","Gene","Type","Mean","SE")
sumG3$Gene<-factor(sumG3$Gene, levels=c("Core","E1", "HVR1","E2","NS1", "NS2","NS3","NS4A","NS4B","NS5A","NS5B"))
sumG3$Type<-factor(sumG3$Type, levels=c("syn","nonsyn","stop"))

write.csv(sumG3, "Output_all/MVF/MVFsummary.byGene.byType.withStop.csv")


### 
# Wilcoxin test for comparing nonsyn vs syn

sumTs<-sameDF[,c("merged.pos", "gene", "mean.1A","mean.1B", "mean.3A","Type.1A","Type.1B","Type.3A")]
Comb<-t(combn(1:3, 2))
combnames<-c("1a-1b","1a-3a","1b-3a")
genenames<-levels(sumG3$Gene)
wilcox.resS<- data.frame(matrix (nrow=22, ncol=3))
wilcox.resNS<- data.frame(matrix(nrow=22, ncol=3))
wilcox.res2<- data.frame(matrix (nrow=22, ncol=3))
wilcox.resST<- data.frame(matrix(nrow=22, ncol=3))
wilcox.NSvsST<-data.frame(matrix(nrow=22, ncol=3))
Stopsites<-data.frame(matrix(nrow=11, ncol=3))

rownames(wilcox.resS)<-c(paste0(genenames,"_greater"),paste0(genenames,"_less"))
colnames(wilcox.resS)<- combnames                 
rownames(wilcox.resNS)<-c(paste0(genenames,"_greater"),paste0(genenames,"_less"))
colnames(wilcox.resNS)<- combnames                   
rownames(wilcox.res2)<-c(paste0(genenames,"_greater"),paste0(genenames,"_less"))
colnames(wilcox.res2)<- subs                    

rownames(wilcox.resST)<-c(paste0(genenames,"_greater"),paste0(genenames,"_less"))
colnames(wilcox.resST)<- combnames            
rownames(wilcox.NSvsST)<-c(paste0(genenames,"_greater"),paste0(genenames,"_less"))
colnames(wilcox.NSvsST)<- subs               

rownames(Stopsites)<-genenames
colnames(Stopsites)<- subs   

for (i in 1:3){
        d<-sumTs[,c("merged.pos", "gene", paste0("mean.",subs[i]),paste0("Type.",subs[i]))]
        colnames(d)[3:4]<-c("mean", "Type")
        
        k=2
        n1<-Comb[i,1]+k
        n2<-Comb[i,2]+k
        
        for (g in 1:11){
                rs1<-wilcox.test(sumTs[sumTs$gene==genenames[g]& sumTs[paste0("Type.",subs[Comb[i,1]])]=="syn" ,n1],  sumTs[sumTs$gene==genenames[g] & sumTs[paste0("Type.",subs[Comb[i,2]])]=="syn",n2], alternative = "greater", paired = FALSE) 
                wilcox.resS[g,i]<-rs1[[3]][1]
                rs2<-wilcox.test(sumTs[sumTs$gene==genenames[g]& sumTs[paste0("Type.",subs[Comb[i,1]])]=="syn" ,n1],  sumTs[sumTs$gene==genenames[g] & sumTs[paste0("Type.",subs[Comb[i,2]])]=="syn",n2], alternative = "less", paired = FALSE) 
                wilcox.resS[(g+11),i]<-rs2[[3]][1]
                
                rn1<-wilcox.test(sumTs[sumTs$gene==genenames[g]& sumTs[paste0("Type.",subs[Comb[i,1]])]=="nonsyn",n1],sumTs[sumTs$gene==genenames[g] & sumTs[paste0("Type.",subs[Comb[i,2]])]=="nonsyn",n2], alternative = "greater", paired = FALSE) 
                wilcox.resNS[g,i]<-rn1[[3]][1]
                rn2<-wilcox.test(sumTs[sumTs$gene==genenames[g]& sumTs[paste0("Type.",subs[Comb[i,1]])]=="nonsyn",n1],sumTs[sumTs$gene==genenames[g] & sumTs[paste0("Type.",subs[Comb[i,2]])]=="nonsyn",n2], alternative = "less", paired = FALSE) 
                wilcox.resNS[(g+11),i]<-rn2[[3]][1]
                
                r1st<-wilcox.test(sumTs[sumTs$gene==genenames[g]&sumTs[paste0("Type.",subs[Comb[i,1]])]=="stop",n1],  sumTs[sumTs$gene==genenames[g] & sumTs[paste0("Type.",subs[Comb[i,2]])]=="stop",n2], alternative = "greater", paired = FALSE) 
                wilcox.resST[g,i]<-r1st[[3]][1]
                r2st<-wilcox.test(sumTs[sumTs$gene==genenames[g]&sumTs[paste0("Type.",subs[Comb[i,1]])]=="stop",n1],  sumTs[sumTs$gene==genenames[g] & sumTs[paste0("Type.",subs[Comb[i,2]])]=="stop",n2], alternative = "less", paired = FALSE) 
                wilcox.resST[(g+11),i]<-r2st[[3]][1]
                               
                rt1<-wilcox.test(d$mean[d$gene==genenames[g]& d$Type=="nonsyn"], d$mean[d$gene==genenames[g] & d$Type=="stop"], alternative = "greater", paired = FALSE) 
                wilcox.NSvsST[g,i]<-rt1[[3]][1]
                rt2<-wilcox.test(d$mean[d$gene==genenames[g]& d$Type=="nonsyn"], d$mean[d$gene==genenames[g] & d$Type=="stop"], alternative = "less", paired = FALSE) 
                wilcox.NSvsST[(g+11),i]<-rt2[[3]][1]
                
                r2<-wilcox.test(d$mean[d$gene==genenames[g]& d$Type=="syn"], d$mean[d$gene==genenames[g] & d$Type=="nonsyn"], alternative = "greater", paired = FALSE) 
                wilcox.res2[g,i]<-r2[[3]][1]    
                r22<-wilcox.test(d$mean[d$gene==genenames[g]& d$Type=="syn"], d$mean[d$gene==genenames[g] & d$Type=="nonsyn"], alternative = "less", paired = FALSE) 
                wilcox.res2[(g+11),i]<-r22[[3]][1]    
                
                Stopsites[g,i]<-length(d$mean[d$gene==genenames[g]& d$Type=="stop"])
                
        }
                
}


write.csv(wilcox.resS, "Output_all/MVF/WilcoxTest_byGene.Syn.csv",row.names = T)
write.csv(wilcox.resNS,"Output_all/MVF/WilcoxTest_byGene.NS.csv",row.names = T)
write.csv(wilcox.res2, "Output_all/MVF/WilcoxTest_NSvsS.csv",row.names = T)
write.csv(wilcox.resST,"Output_all/MVF/WilcoxTest_byGene.stop.csv",row.names = T)
write.csv(wilcox.NSvsST,"Output_all/MVF/WilcoxTest_NSvsStop.csv",row.names = T)

apply(wilcox.resS, 1, function(x) x=length(x[x<0.05]) ) 
#syn HVR1 1b-3a is significant
wilcox.resS$`1b-3a`[3]<-wilcox.resS$`1b-3a`[14]
apply(wilcox.resNS, 1, function(x) x=length(x[x<0.05]) ) 
apply(wilcox.res2, 1, function(x) x=length(x[x<0.05]) ) 
apply(wilcox.resST, 1, function(x) x=length(x[x<0.05]) ) 
apply(wilcox.NSvsST, 1, function(x) x=length(x[x<0.05]) ) 

Stopsites
#     1A 1B 3A
#Core 17 18 20
#E1   19 15 19
#HVR1  1  1  2
#E2   55 54 32
#NS1   7  5  7
#NS2  26 22 26
#NS3  36 35 38
#NS4A  3  3  4
#NS4B 26 26 28
#NS5A 44 52 35
#NS5B 21 21 19

# Do not use HVR1 data for stop sites. Too small data points.


Test.results<-data.frame(gene=genenames)
for (i in 1:11){
        symb1<-NA
        symb2<-NA
        symb3<-NA
        if (wilcox.resS$`1a-1b`[i]<0.05) symb1="a"
        if (wilcox.resS$`1a-3a`[i]<0.05) symb2="b"
        if (wilcox.resS$`1b-3a`[i]<0.05) symb3="c"
        
        symb<-c(symb1, symb2, symb3) 
        symb<-symb[!is.na(symb)]
        symb<-paste(symb, collapse = ',' )
        if (length(symb)==0) symb<-''
        Test.results$Syn[i]<-symb
        
        nsymb1<-NA
        nsymb2<-NA
        nsymb3<-NA
        if (wilcox.resNS$`1a-1b`[i]<0.05) nsymb1="a"
        if (wilcox.resNS$`1a-3a`[i]<0.05) nsymb2="b"
        if (wilcox.resNS$`1b-3a`[i]<0.05) nsymb3="c"
        
        nsymb<-c(nsymb1, nsymb2, nsymb3) 
        nsymb<-nsymb[!is.na(nsymb)]
        nsymb<-paste(nsymb, collapse=',')
        if (length(nsymb)==0) nsymb<-''
        Test.results$Nonsyn[i]<-nsymb
        
        Stsymb1<-NA
        Stsymb2<-NA
        Stsymb3<-NA
        if (wilcox.resST$`1a-1b`[i]<0.05) Stsymb1="a"
        if (wilcox.resST$`1a-3a`[i]<0.05) Stsymb2="b"
        if (wilcox.resST$`1b-3a`[i]<0.05) Stsymb2="c"
        
        Stsymb<-c(Stsymb1, Stsymb2, Stsymb3) 
        Stsymb<-Stsymb[!is.na(Stsymb)]
        Stsymb<-paste(Stsymb, collapse=',')
        if (length(Stsymb)==0) Stsymb<-''
        Test.results$Stop[i]<-Stsymb
}

      
ggplot(sumG3, aes(x=Gene, y=Mean, shape=Type, color=Subtype, fill=Subtype))+
        geom_point(aes(group=factor(Subtype), color=factor(Subtype)),position=position_dodge(width=0.8),size =2)+
        scale_color_manual(values=colors2[c(1,3,5)])+
        scale_fill_manual(values=paste0(colors2[c(1,3,5)],"66"), guide = 'none')+
        scale_shape_manual(values=c(16,4,23), labels=c("Syn","Nonsyn","Nonsense"))+
        geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE, group=factor(Subtype), color=factor(Subtype)), width=.2, size=.2, position=position_dodge(width=0.8))+
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
        labs(shape="Mutation type")+
        annotate(geom="text", x=c(1:2,4:11), y=-0.0002, label="***",color ='gray40', size=3)+
        annotate(geom="text", x=c(1:2,4:7,9:11), y=-0.0005, label="***",color ='gray60', size=3)+
        annotate(geom="text", x=8, y=-0.0005, label="  *",color ='gray60', size=3)+
        annotate(geom="text", x=3, y=-0.0005, label="na",color ='gray', size=3)+
        annotate(geom="text", x=1:11, y=0.021,  label=Test.results$Syn,, color="gray40", size=2.5 )+
        annotate(geom="text", x=1:11, y=0.0205,  label=Test.results$Nonsyn,, color="gray60", size=2.5 )+
        annotation_custom(textGrob(label = "Syn", hjust = 0, gp = gpar(cex = .5, col="gray40")),  xmin=11.7,xmax=11.7,ymin=0.021,ymax=0.021)+
        annotation_custom(textGrob(label = "Nonyn", hjust = 0, gp = gpar(cex = .5,col="gray60")), xmin=11.7,xmax=11.7,ymin=0.0205,ymax=0.0205)+
        #annotation_custom(textGrob(label = "S vs NS", hjust = 0, gp = gpar(cex = .5,col="gray40")), xmin=11.7,xmax=11.7,ymin=0,ymax=0)+
        annotation_custom(pointsGrob(pch=16, gp=gpar(cex = .4,col="gray60")), xmin=11.8,xmax=11.8,ymin=0,ymax=0 )+
        annotation_custom(textGrob(label ="vs", hjust = 0, gp = gpar(cex = .5,col="gray40")), xmin=11.95,xmax=11.95,ymin=0,ymax=0)+
        annotation_custom(pointsGrob(pch=4, gp=gpar(cex = .35,col="gray60")), xmin=12.3,xmax=12.3,ymin=0,ymax=0)+
        annotation_custom(pointsGrob(pch=4, gp=gpar(cex = .35,col="gray60")), xmin=11.8,xmax=11.8,ymin=-0.0006,ymax=-0.0006)+
        annotation_custom(textGrob(label ="vs", hjust = 0, gp = gpar(cex = .5,col="gray60")), xmin=11.95,xmax=11.95,ymin=-0.0006,ymax=-0.0006)+
        annotation_custom(pointsGrob(pch=23, gp=gpar(cex = .4,col="gray60")), xmin=12.3,xmax=12.3,ymin=-0.0006,ymax=-0.0006 )+
        #annotation_custom(textGrob(label = "NS vs ", hjust = 0, gp = gpar(cex = .5,col="gray60")), xmin=11.7,xmax=11.7,ymin=-0.0005,ymax=-0.0005)+
        coord_cartesian(clip = "off")

ggsave(filename="Output_all/Figures/MVF.byGene_byType_stats.pdf", width = 8, height = 5)
