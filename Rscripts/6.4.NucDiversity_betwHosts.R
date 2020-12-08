# Compare nucleotide diersity between hosts vs. within hosts
library(ggplot2)
library(colorspace)
library(dplyr)
library(reshape2)
source("Rscripts/baseRscript.R")
subs<-c("1A","1B","3A")

colors2<-qualitative_hcl(6, palette="Dark3")
#col2_light<-qualitative_hcl(6, palette="Set3")
#div.colors<-c(colors2[1],col2_light[1],colors2[3],col2_light[3],colors2[5],col2_light[5])


## Calculate pi (nuc diversity) per gene by genotype and plot them

mergedmeta<-read.csv("Output_all/merged.metadata.csv", stringsAsFactors = F, row.names = 1)
fastalist<-list.files("Data/", pattern="NCBI.fasta")

genes<-unique(mergedmeta$gene)
genes<-genes[-1]

for (i in 1:3){
        gtype<-subs[i]
        seq<-read.fasta(paste0("Data/", fastalist[i]))
        Seq<-do.call(rbind, seq)
        Seq<-data.frame(t(Seq), stringsAsFactors = F)
        
        if (i!=3) n=342
        if (i==3) n=340
        
        Seq$pos<-n:(n+nrow(Seq)-1)
        
        meta<-mergedmeta[,c("merged.pos","gene","codon",paste0("org.pos.",gtype),paste0("ref.",gtype),
                            paste0("Type.",gtype),paste0("makesCpG.",gtype) )]
        Seq<-merge(meta, Seq, by.x= paste0("org.pos.",gtype), by.y="pos")
        
        write.csv(Seq, paste0("Output_all/betwHost/Sequence_metadata.",gtype,".csv"))
}

for (i in 1:3){
        gtype<-subs[i]
        Seq<-read.csv(paste0("Output_all/betwHost/Sequence_metadata.",gtype,".csv"), stringsAsFactors = F, row.names = 1)
        
        nucdiv<-data.frame(gene=genes)
        
        for (j in 1:length(genes)){
                genedf<-Seq[Seq$gene==genes[j], ]
                if (i==3&j==3) {
                        genedf<-genedf[,-which(colnames(genedf)=="KY620815")]
                        genedf<-genedf[,-which(colnames(genedf)=="KY620411")]
                        genedf<-genedf[,-which(colnames(genedf)=="KY620638")]
                } 
                if (i==3&j==5) {
                        genedf<-genedf[,-which(colnames(genedf)=="KY620815")]
                        genedf<-genedf[,-which(colnames(genedf)=="KY620683")]
                        genedf<-genedf[,-which(colnames(genedf)=="KY620544")]
                } 
                geneSyn<-genedf[genedf[paste0("Type.",gtype)]=="syn",]
                geneNS<-genedf[genedf[paste0("Type.",gtype)]=="nonsyn",]
                s<-t(geneSyn[,8:ncol(genedf)])
                #replace - with N
                s[s=="-"]<-"n"
                Syn<-as.DNAbin(s)
                
                ns<-t(geneNS[,8:ncol(genedf)])
                ns[ns=="-"]<-"n"
                
                Nonsyn<-as.DNAbin(ns)
                nucdiv[j,"syn"]<-nuc.div(Syn, pairwise.deletion=T)
                nucdiv[j,"nonsyn"]<-nuc.div(Nonsyn, pairwise.deletion=T)
        }
        write.csv(nucdiv, paste0("Output_all/betwHost/NucDiv.by.gene.",gtype,".csv"))
}


### 
######
d<-data.frame()
for (i in 1:3){
        df<-read.csv(paste0("Output_all/betwHost/NucDiv.by.gene.",subs[i],".csv"))
        df$Subtype<-subs[i]
        d<-rbind(d,df )
}
        
d<-d[,-1]
write.csv(d, "Output_all/betwHost/NucDiv_summary_3Subtypes.csv")


Ndiv<-melt(d, id.vars=c("gene","Subtype"))
colnames(Ndiv)[1]<-"Gene"
colnames(Ndiv)[3:4]<-c("Type","Pi")
        
Ndiv$Gene<-factor(Ndiv$Gene, levels=c("Core","E1", "HVR1","E2","P7", "NS2","NS3","NS4A","NS4B","NS5A","NS5B"))

## nucleotide diversity plot
ggplot(Ndiv, aes(x=Gene, y=Pi, shape=Type,color=Subtype))+
        geom_point(position=position_dodge(width=0.8),size =2)+
        scale_color_manual(values=colors2[c(1,3,5)], labels=c("1a","1b",'3a'))+
        scale_shape_manual(values=c(19,4))+
        labs(color="Subtype")+
        labs(shape="Mutation type")+
        #geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width=.2, size=.2, position=position_dodge(width=0.8))+
        theme(axis.title.x=element_blank())+ylab(expression(paste("Nucleotide diversity (",pi,")")))+
        theme_bw()+
        labs(x="")+
        theme(axis.text.x = element_text(size=12),axis.title.y = element_text(size=13))+
        theme(panel.grid.major.x = element_blank())+
        geom_vline(xintercept = c(1:10)+0.5,  
                   color = "gray60", size=.4)+
        guides(shape = guide_legend(order = 2),color=guide_legend(order=1))

ggsave(filename="Output_all/Figures/betwHost_Pi.byType.byGene.pdf",width = 8, height = 5)


#plot PiN/PiS
Pitb<-data.frame()
PiNS<-d
for (g in 1:3){
        pi.tb<-data.frame(Gene=genes)
        for ( i in 1:11) { 
            pi.tb$pNpS[i]<-PiNS$nonsyn[PiNS$Subtype==subs[g] & PiNS$gene==genes[i]]/PiNS$syn[PiNS$Subtype==subs[g] & PiNS$gene==genes[i] ]
        }
        pi.tb$Subtype<-subs[g]
        Pitb<-rbind(Pitb, pi.tb)
}

Pitb$Gene<-factor(Pitb$Gene, levels=c("Core","E1", "HVR1","E2","P7", "NS2","NS3","NS4A","NS4B","NS5A","NS5B"))

#point plot
ggplot(Pitb, aes(x=Gene, y=pNpS, color=Subtype))+
        geom_point(position=position_dodge(width=0.8),size =3)+
        scale_color_manual(values=colors2[c(1,3,5)], labels=c("1a","1b","3a"))+
        scale_shape_manual(values=c(4,19))+
        theme(axis.title.x=element_blank())+ylab(expression(paste(pi[N],"/",pi[S])))+
        theme_bw()+
        labs(x="")+
        theme(axis.text.x = element_text(size=12),axis.title.y = element_text(size=13))+
        theme(panel.grid.major.x = element_blank())+
        geom_vline(xintercept = c(1:10)+0.5,  color = "gray60", size=.4)+
        geom_hline(yintercept = 1,  color = "gray40", size=.3, linetype=4)+
        labs(color="Subtype")

ggsave(filename="Output_all/betwHost/PiN.over.piS.bygene_bygenotype.pdf",width = 8, height = 5)


## Within-host pi data ##
InvivoPi<-read.csv("Output_all/Diversity/PiNS_Summary.csv", stringsAsFactors = F, row.names = 1)

#plot PiN/PiS
InvivoP<-data.frame()
for (g in 1:3){
        pi.tb<-data.frame(Gene=genes)
        for (i in 1:11) { 
                pi.tb$pNpS[i]<-InvivoPi$Mean[InvivoPi$Subtype==subs[g] & InvivoPi$Gene==genes[i] & InvivoPi$Type=="Nonsyn"]/
                        InvivoPi$Mean[InvivoPi$Subtype==subs[g] & InvivoPi$Gene==genes[i] & InvivoPi$Type=="Syn"]
        }
        pi.tb$Subtype<-subs[g]
        InvivoP<-rbind(InvivoP, pi.tb)
}

InvivoP$Level<-"Within-host"
Pitb$Level<-"Between-host"

Pitb2<-rbind(InvivoP, Pitb)
Pitb2$Level<-factor(Pitb2$Level, levels=c("Within-host", "Between-host"))
Pitb2$Gene<-factor(Pitb2$Gene, levels=c("Core","E1", "HVR1","E2","P7", "NS2","NS3","NS4A","NS4B","NS5A","NS5B"))
colnames(Pitb2)[3]<-"Subtype"

ggplot(Pitb2, aes(x=Gene, y=pNpS, shape=Level,color=Subtype, fill=Subtype))+
        geom_point(aes(group=factor(Subtype), color=factor(Subtype)),position=position_dodge(width=0.8),size =2)+
        
        #geom_point(position=position_dodge(width=0.8),size =2)+
        scale_color_manual(values=colors2[c(1,3,5)], labels=c("1a","1b","3a"))+
        scale_shape_manual(values=c(19,25))+
        scale_fill_manual(values=paste0(colors2[c(1,3,5)],"66"), guide = 'none')+
        guides(shape = guide_legend(title = NULL))+
          
        theme(axis.title.x=element_blank())+ylab(expression(paste(pi[N],"/",pi[S])))+
        theme_bw()+
        guides(shape = guide_legend(order = 2, title=NULL),color=guide_legend(order=1))+
        labs(x="")+
        theme(axis.text.x = element_text(size=12),axis.title.y = element_text(size=13))+
        theme(panel.grid.major.x = element_blank())+
        geom_vline(xintercept = c(1:10)+0.5,  
                   color = "gray60", size=.4)+
                geom_hline(yintercept = 1,  color = "gray40", size=.3, linetype=4)

ggsave(filename="Output_all/Figures/pNpS.byLevel.pdf",width = 8, height = 5)



