#Beta regression preparation:
library(tidyverse)
library(zoo)
library(purrr)
library(miscTools)
source("Rscripts/baseRscript.R")

# Read the summary data (only Maj==Ref)
subs<-c("1A", "1B","3A")
sum<-read.csv("Output_all/MVF/Ts.Same_Mean_3subtypes.csv", stringsAsFactors = F, row.names = 1)

#### For preparing Beta Regression data formatting
nucord <- c("a", "t", "c", "g")

for (k in 2:length(subs)){
        name<-subs[k]
        dat<-sum[,c("merged.pos",paste0("makesCpG.",name),paste0("bigAAChange.",name),paste0("mean.",name))]
        brform<-as.matrix(data.frame(a="",t="",c="",g="", Syn ="",Nonsyn="",Stop=""))
        nuc<-as.matrix(data.frame(a="",t="",c="",g=""))
        for (i in 1:nrow(sum)){
                atcg <- c(0,0,0,0)
                atcg[which(nucord == sum[i,paste0("ref.",name)])] <- 1
                nuc<-insertRow(nuc,i,atcg)
                nonsyn <- as.numeric(regexpr("nonsyn",sum[i,paste0("Type.",name)]) > 0)
                stop <- as.numeric(regexpr("stop",sum[i,paste0("Type.",name)]) > 0)
                syn<-as.numeric(regexpr("^syn",sum[i,paste0("Type.",name)]) > 0)
                new<-c(atcg,syn,nonsyn,stop)
                brform<-insertRow(brform,i,new)
        }
        BrData<-cbind(dat$merged.pos,brform[1:nrow(dat),])
        BrData<-cbind(BrData,dat[,2:4])
        colnames(BrData)[1]<-"pos"
        colnames(BrData)[9]<-"CpG"
        colnames(BrData)[10]<-"bigAAChange"
        colnames(BrData)[11]<-'mean'
        BrData <-BrData[!is.na(BrData$mean),]
        write.csv(BrData, paste0("Output_all/Betareg/Betareg_",name,".csv"))
}

#Attache the subtype info
brdata<-list()
subs<-c("1A", "1B","3A")
for (i in 1:length(subs)){
    df<-read.csv(paste0("Output_all/Betareg/Betareg_",subs[i],".csv"),stringsAsFactors = F, row.names = 1)
    brdata[[i]]<-df
    names(brdata)[i]<-subs[i]
}        

### addd gene annotation info for all genes
genes<-read.csv("Data/HCV_annotations_joined.csv")
genes$Gene<- factor(genes$Gene, levels=c("5' UTR","Core","E1", "HVR1","E2","P7", "NS2","NS3","NS4A","NS4B","NS5A","NS5B"))

genenames<-genes$Gene
gene<-c()
for (i in 2:12){
    gene<-c(gene, rep(i, times=genes$end[i]-genes$start[i]+1))
}

n<-data.frame(pos=342:(length(gene)+341))
g<-cbind(n,gene)
for (k in 1:length(brdata)){
    dat1<-brdata[[k]]
    brData<-merge(dat1, g, by ="pos")
    
    for (i in 2:12){
        gname<-paste(genes$Gene[i])
        n<-ncol(brData)
        brData[,n+1]<-0
        colnames(brData)[n+1]<-gname
        brData[brData$gene==i,n+1]<-1
    }
     write.csv(brData,paste0("Output_all/Betareg/brData_withGenes.",subs[k],".csv"))
}
