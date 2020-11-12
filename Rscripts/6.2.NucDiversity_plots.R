#Use overveiw2 data to avoid losing data in the filtered overview (overview3) for Fst calculation
library(ggplot2)
library(reshape2)
#library(ggpubr)
#library(ggthemes)
#library(plotrix)
library(grid)
#library(purrr)
#library(zoo)
#library(tidyverse)
#library(gplots)
library(colorspace)

source("Rscripts/baseRscript.R")
subs<-c("1A","1B","3A")
colors2<-qualitative_hcl(6, palette="Dark3")
col2_light<-qualitative_hcl(6, palette="Set3")
div.colors<-c(colors2[1],col2_light[1],colors2[3],col2_light[3],colors2[5],col2_light[5])


## Calculate pi (nuc diversity) per gene by Subtype and plot them
### addd gene annotation info for all genes
Genes<-read.csv("Data/HCV_annotations_joined.csv", stringsAsFactors = F)
Genes$Gene[6]<-"NS1"
gene.vector<-c()
for (i in 1:12){
        gene.vector<-c(gene.vector, rep(Genes$Gene[i],times=Genes$start[i+1]-Genes$start[i]))
}

n<-data.frame(pos=1:length(gene.vector))
gene.vec<-cbind(n,gene.vector)
colnames(gene.vec)[2]<-"gene"


#Ds, paste0("Output_all/Diversity/Ds_",subs[g],".csv"))
Ds.1A<-read.csv("Output_all/Diversity/Ds_1A.csv", stringsAsFactors = F, row.names = 1)
Ds.1B<-read.csv("Output_all/Diversity/Ds_1B.csv", stringsAsFactors = F, row.names = 1)
Ds.3A<-read.csv("Output_all/Diversity/Ds_3A.csv", stringsAsFactors = F, row.names = 1)


#depth info
Summary<-read.csv("Output_all/MVF/Ts.MV_Mean_3subtypes_depth.csv",stringsAsFactors = F, row.names = 1)


## Merge metadata to D values caluclated in 6.1.
merged.meta<-read.csv("Output_all/merged.metadata.csv", stringsAsFactors = F, row.names = 1)
meta1A<-merged.meta[,c("org.pos.1A","Type.1A")]
colnames(meta1A)[1]<-"pos"
Pi1A<-merge(meta1A, Ds.1A, by="pos")
colnames(Pi1A)[2:3]<-c("Type","gene")

meta1B<-merged.meta[,c("org.pos.1B","Type.1B")]
colnames(meta1B)[1]<-"pos"
Pi1B<-merge(meta1B, Ds.1B, by="pos")
colnames(Pi1B)[2:3]<-c("Type","gene")

meta3A<-merged.meta[,c("org.pos.3A","Type.3A")]
colnames(meta3A)[1]<-"pos"
Pi3A<-merge(meta3A, Ds.3A, by="pos")
colnames(Pi3A)[2:3]<-c("Type","gene")

 
#remove 5'UTR and calculate piN/piS ratios
PiNS<-data.frame()
for (g in 1:3){
        df<-get(paste0("Pi",subs[g]))
        dfN<-df[df$gene!="5' UTR" & df$Type=="nonsyn",]
        dfS<-df[df$gene!="5' UTR" & df$Type=="syn",]
        
        piN<-data.frame(Gene=Genes$Gene[2:12])     
        piS<-data.frame(Gene=Genes$Gene[2:12])     

        for (i in 1:11){
                dfN2<-dfN[dfN$gene==Genes$Gene[(i+1)],]
                piN[i, "Mean"]<-sum(dfN2[,4:ncol(dfN2)], na.rm=T)/ sum(!is.na(dfN2[,3:ncol(dfN2)]))
                piN[i,"SE"]<-mean(std.error(dfN2[,3:ncol(dfN2)], na.rm=T), na.rm=T)
                dfS2<-dfS[dfS$gene==Genes$Gene[(i+1)],]
                piS[i, "Mean"]<-sum(dfS2[,4:ncol(dfS2)], na.rm=T)/ sum(!is.na(dfS2[,3:ncol(dfS2)]))
                piS[i,"SE"]<-mean(std.error(dfS2[,3:ncol(dfS2)], na.rm=T), na.rm=T)
        }
       
        piN$Group<-paste0(subs[g],".Nonsyn")
        piN$Type<-"Nonsyn"
        piS$Group<-paste0(subs[g],".Syn")
        piS$Type<-"Syn"
        pi<-rbind(piN,piS)
        pi$Subtype<-subs[g]
        PiNS<-rbind(PiNS, pi)
        
}

PiNS$Gene<-factor(PiNS$Gene, levels=c("Core","E1", "HVR1","E2","NS1", "NS2","NS3","NS4A","NS4B","NS5A","NS5B"))
write.csv(PiNS, "Output_all/Diversity/PiNS_Summary.csv")


## Plot Pi
#ggplot(PiNS, aes(x=Gene, y=Mean, shape=Type,color=Subtype))+
#        geom_point(position=position_dodge(width=0.8),size =2)+scale_color_manual(values=colors2[c(1,3,5)])+
#        scale_shape_manual(values=c(4,19))+
#        #geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width=.2, size=.2, position=position_dodge(width=0.8))+
#        theme(axis.title.x=element_blank())+ylab(expression(paste("Nucleotide diversity (",pi,")")))+
#        theme_bw()+
#        labs(x="")+
#        theme(axis.text.x = element_text(size=12),axis.title.y = element_text(size=13))+
#        theme(panel.grid.major.x = element_blank())+
#        geom_vline(xintercept = c(1:10)+0.5,  
#                   color = "gray60", size=.4)+
#        ylim(0.005,0.055)
#
##ggsave(filename="Output_all/Figures/Pi.byType.byGene.bySubtype_points.pdf",width = 8, height = 5)

#plot PiN/PiS
Pitb<-data.frame()
for (g in 1:3){
        pi.tb<-data.frame(Gene=Genes$Gene[2:12])
        for ( i in 1:11) { 
            pi.tb$pNpS[i]<-PiNS$Mean[PiNS$Subtype==subs[g] & PiNS$Gene==Genes$Gene[i+1] & PiNS$Type=="Nonsyn"]/PiNS$Mean[PiNS$Subtype==subs[g] & PiNS$Gene==Genes$Gene[i+1] & PiNS$Type=="Syn"]
        }
        pi.tb$Subtype<-subs[g]
        Pitb<-rbind(Pitb, pi.tb)
}
Pitb$Gene<-factor(Pitb$Gene, levels=c("Core","E1", "HVR1","E2","NS1", "NS2","NS3","NS4A","NS4B","NS5A","NS5B"))

ggplot(Pitb, aes(x=Gene, y=pNpS, color=Subtype))+
        geom_point(position=position_dodge(width=0.8),size =3)+scale_color_manual(values=colors2[c(1,3,5)])+
        scale_shape_manual(values=c(4,19))+
        theme(axis.title.x=element_blank())+ylab(expression(paste(pi[N],"/",pi[S])))+
        theme_bw()+
        labs(x="")+
        theme(axis.text.x = element_text(size=12),axis.title.y = element_text(size=13))+
        theme(panel.grid.major.x = element_blank())+
        geom_vline(xintercept = c(1:10)+0.5,  color = "gray60", size=.4)+
        geom_hline(yintercept = 1,  color = "gray40", size=.3, linetype=4)+
        labs(color="Subtype")
ggsave(filename="Output_all/Figures/PiN.over.piS.bygene_bySubtype.point.pdf",width = 8, height = 5)


# Stats --Wilcoxin test of pi values 
Comb<-t(combn(1:3, 2))
combnames<-c("1A-1B","1A-3A","1B-3A")
genenames<-as.character(unique(Pitb$Gene))

#data formatting
Pi1A$mean.1A<-rowMeans(Pi1A[4:ncol(Pi1A)], na.rm = T)
meta1a<-merged.meta[,c("merged.pos","org.pos.1A","Type.1A")]
colnames(meta1a)[2]<-"pos"
nuc.div<-merge(meta1a, Pi1A[,c(1:3,199)], by="pos", all.x=T)
nuc.div<-nuc.div[,c("merged.pos", "gene", "Type.1A","mean.1A")]

Pi1B$mean.1B<-rowMeans(Pi1B[4:ncol(Pi1B)], na.rm = T)
meta1b<-merged.meta[,c("merged.pos","org.pos.1B","Type.1B")]
colnames(meta1b)[2]<-"pos"
nuc.div2<-merge(meta1b, Pi1B[,c(1:3,ncol(Pi1B))], by="pos", all.x=T)
nuc.div2<-nuc.div2[,-c(1,4)]

nuc.div<-merge(nuc.div, nuc.div2[,c(1,2,4)], by="merged.pos")

Pi3A$mean.3A<-rowMeans(Pi3A[4:ncol(Pi3A)], na.rm = T)
meta3a<-merged.meta[,c("merged.pos","org.pos.3A","Type.3A")]
colnames(meta3a)[2]<-"pos"
nuc.div2<-merge(meta3a, Pi3A[,c(1:3,ncol(Pi3A))], by="pos", all.x=T)
nuc.div2<-nuc.div2[,-c(1,4)]

nuc.div<-merge(nuc.div, nuc.div2[,c(1,2,4)], by="merged.pos")
nuc.div<-nuc.div[,c("merged.pos","gene","mean.1A","mean.1B","mean.3A","Type.1A","Type.1B","Type.3A")]

#prepare empty data frames
wilcox.resS<- data.frame(matrix (nrow=22, ncol=3))
wilcox.resNS<- data.frame(matrix(nrow=22, ncol=3))
wilcox.res2<- data.frame(matrix (nrow=22, ncol=3))
wilcox.resST<- data.frame(matrix(nrow=22, ncol=3))
wilcox.NSvsST<-data.frame(matrix(nrow=22, ncol=3))

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

for (i in 1:3){
    d<-nuc.div[,c("merged.pos", "gene", paste0("mean.",subs[i]),paste0("Type.",subs[i]))]
    colnames(d)[3:4]<-c("mean", "Type")
    
    k=2
    n1<-Comb[i,1]+k
    n2<-Comb[i,2]+k
    
    for (g in 1:11){
        rs1<-wilcox.test(nuc.div[nuc.div$gene==genenames[g]& nuc.div[paste0("Type.",subs[Comb[i,1]])]=="syn" ,n1],  nuc.div[nuc.div$gene==genenames[g] & nuc.div[paste0("Type.",subs[Comb[i,2]])]=="syn",n2], alternative = "greater", paired = FALSE) 
        wilcox.resS[g,i]<-rs1[[3]][1]
        rs2<-wilcox.test(nuc.div[nuc.div$gene==genenames[g]& nuc.div[paste0("Type.",subs[Comb[i,1]])]=="syn" ,n1],  nuc.div[nuc.div$gene==genenames[g] & nuc.div[paste0("Type.",subs[Comb[i,2]])]=="syn",n2], alternative = "less", paired = FALSE) 
        wilcox.resS[(g+11),i]<-rs2[[3]][1]
        
        rn1<-wilcox.test(nuc.div[nuc.div$gene==genenames[g]& nuc.div[paste0("Type.",subs[Comb[i,1]])]=="nonsyn",n1],nuc.div[nuc.div$gene==genenames[g] & nuc.div[paste0("Type.",subs[Comb[i,2]])]=="nonsyn",n2], alternative = "greater", paired = FALSE) 
        wilcox.resNS[g,i]<-rn1[[3]][1]
        rn2<-wilcox.test(nuc.div[nuc.div$gene==genenames[g]& nuc.div[paste0("Type.",subs[Comb[i,1]])]=="nonsyn",n1],nuc.div[nuc.div$gene==genenames[g] & nuc.div[paste0("Type.",subs[Comb[i,2]])]=="nonsyn",n2], alternative = "less", paired = FALSE) 
        wilcox.resNS[(g+11),i]<-rn2[[3]][1]
        
        r1st<-wilcox.test(nuc.div[nuc.div$gene==genenames[g]&nuc.div[paste0("Type.",subs[Comb[i,1]])]=="stop",n1],  nuc.div[nuc.div$gene==genenames[g] & nuc.div[paste0("Type.",subs[Comb[i,2]])]=="stop",n2], alternative = "greater", paired = FALSE) 
        wilcox.resST[g,i]<-r1st[[3]][1]
        r2st<-wilcox.test(nuc.div[nuc.div$gene==genenames[g]&nuc.div[paste0("Type.",subs[Comb[i,1]])]=="stop",n1],  nuc.div[nuc.div$gene==genenames[g] & nuc.div[paste0("Type.",subs[Comb[i,2]])]=="stop",n2], alternative = "less", paired = FALSE) 
        wilcox.resST[(g+11),i]<-r2st[[3]][1]
        
        rt1<-wilcox.test(d$mean[d$gene==genenames[g]& d$Type=="nonsyn"], d$mean[d$gene==genenames[g] & d$Type=="stop"], alternative = "greater", paired = FALSE) 
        wilcox.NSvsST[g,i]<-rt1[[3]][1]
        rt2<-wilcox.test(d$mean[d$gene==genenames[g]& d$Type=="nonsyn"], d$mean[d$gene==genenames[g] & d$Type=="stop"], alternative = "less", paired = FALSE) 
        wilcox.NSvsST[(g+11),i]<-rt2[[3]][1]
        
        r2<-wilcox.test(d$mean[d$gene==genenames[g]& d$Type=="syn"], d$mean[d$gene==genenames[g] & d$Type=="nonsyn"], alternative = "greater", paired = FALSE) 
        wilcox.res2[g,i]<-r2[[3]][1]    
        r22<-wilcox.test(d$mean[d$gene==genenames[g]& d$Type=="syn"], d$mean[d$gene==genenames[g] & d$Type=="nonsyn"], alternative = "less", paired = FALSE) 
        wilcox.res2[(g+11),i]<-r22[[3]][1]    
        
    }
    
}

write.csv(wilcox.resS,"Output_all/Diversity/WilcoxTest.Pi_byGene.Syn.csv",row.names = T)
write.csv(wilcox.resNS,"Output_all/Diversity/WilcoxTest.Pi_byGene.NS.csv",row.names = T)
write.csv(wilcox.res2,"Output_all/Diversity/WilcoxTest.Pi_NSvsS.csv",row.names = T)
write.csv(wilcox.resST,"Output_all/Diversity/WilcoxTest.Pi_byGene.Stop.csv",row.names = T)
write.csv(wilcox.NSvsST,"Output_all/Diversity/WilcoxTest.Pi_NSvsStop.csv",row.names = T)
###
apply(wilcox.res2, 1, function(x) x=length(x[x<0.05]) ) 
apply(wilcox.NSvsST, 1, function(x) x=length(x[x<0.05]) ) 
apply(wilcox.resS, 1, function(x) x=length(x[x<0.05]) ) 
apply(wilcox.resNS, 1, function(x) x=length(x[x<0.05]) ) 
apply(wilcox.resST, 1, function(x) x=length(x[x<0.05]) ) 


## Plot with stat results
Test.results<-data.frame(gene=genenames)
for (i in 1:11){
        symb1<-NA
        symb2<-NA
        symb3<-NA
        if (wilcox.resS$`1A-1B`[i]<0.05) symb1="a"
        if (wilcox.resS$`1A-3A`[i]<0.05) symb2="b"
        if (wilcox.resS$`1B-3A`[i]<0.05) symb3="c"
        
        symb<-c(symb1, symb2, symb3) 
        symb<-symb[!is.na(symb)]
        symb<-paste(symb, collapse = ',' )
        if (length(symb)==0) symb<-''
        Test.results$Syn[i]<-symb
        
        nsymb1<-NA
        nsymb2<-NA
        nsymb3<-NA
        if (wilcox.resNS$`1A-1B`[i]<0.05) nsymb1="a"
        if (wilcox.resNS$`1A-3A`[i]<0.05) nsymb2="b"
        if (wilcox.resNS$`1B-3A`[i]<0.05) nsymb2="c"
        
        nsymb<-c(nsymb1, nsymb2, nsymb3) 
        nsymb<-nsymb[!is.na(nsymb)]
        nsymb<-paste(nsymb, collapse=',')
        if (length(nsymb)==0) nsymb<-''
        Test.results$Nonsyn[i]<-nsymb
        
        Stsymb1<-NA
        Stsymb2<-NA
        Stsymb3<-NA
        if (wilcox.resST$`1A-1B`[i]<0.05) Stsymb1="a"
        if (wilcox.resST$`1A-3A`[i]<0.05) Stsymb2="b"
        if (wilcox.resST$`1B-3A`[i]<0.05) Stsymb2="c"
        
        Stsymb<-c(Stsymb1, Stsymb2, Stsymb3) 
        Stsymb<-Stsymb[!is.na(Stsymb)]
        Stsymb<-paste(Stsymb, collapse=',')
        if (length(Stsymb)==0) Stsymb<-''
        Test.results$Stop[i]<-Stsymb
}

PiNS$Type<-factor(PiNS$Type,levels=c("Syn","Nonsyn"))
ggplot(PiNS, aes(x=Gene, y=Mean, shape=Type,color=Subtype))+
        geom_point(aes(group=factor(Subtype), color=factor(Subtype)), position=position_dodge(width=0.8),size =2)+scale_color_manual(values=colors2[c(1,3,5)])+
        scale_shape_manual(values=c(19,4))+
        #geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width=.2, size=.2, position=position_dodge(width=0.8))+
        theme(axis.title.x=element_blank())+ylab(expression(paste("Nucleotide diversity (",pi,")")))+
        theme_bw()+
        labs(x="")+
        theme(axis.text.x = element_text(size=12),axis.title.y = element_text(size=13))+
        theme(panel.grid.major.x = element_blank())+
        geom_vline(xintercept = c(1:10)+0.5,  
                   color = "gray60", size=.4)+
        guides(shape = guide_legend(order = 2),
                color=guide_legend(order=1))+
        ylim(0.005,0.045)+labs(shape="Mutation type", color="Subtype")+
        annotate(geom="text", x=c(1:2,4:11), y=0.005, label="***",color ='gray40', size=3)+
        annotate(geom="text", x=1:11, y=0.045,  label=Test.results$Syn,, color="gray40", size=2.5)+
        annotate(geom="text", x=1:11, y=0.0438,  label=Test.results$Nonsyn,, color="gray60", size=2.5 )+
        annotation_custom(textGrob(label = "Syn", hjust = 0, gp = gpar(cex = .5,col="gray40")), xmin=11.7,xmax=11.7,ymin=0.045,ymax=0.045)+
        annotation_custom(textGrob(label = "Nonyn", hjust = 0, gp = gpar(cex = .5,col="gray60")), xmin=11.7,xmax=11.7,ymin=0.0438,ymax=0.0438)+
        annotation_custom(pointsGrob(pch=16, gp=gpar(cex = .4,col="gray60")), xmin=11.8,xmax=11.8,ymin=0.005,ymax=0.005 )+
        annotation_custom(textGrob(label ="vs", hjust = 0, gp = gpar(cex = .5,col="gray40")), xmin=11.95,xmax=11.95,ymin=0.005,ymax=0.005)+
        annotation_custom(pointsGrob(pch=4, gp=gpar(cex = .35,col="gray60")), xmin=12.3,xmax=12.3,ymin=0.005,ymax=0.005 )+
        coord_cartesian(clip = "off")
ggsave(filename="Output_all/Figures/Pi.byType.byGene.bySubtype_points.pdf",width = 8, height = 5)
