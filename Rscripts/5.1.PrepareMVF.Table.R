library(stringr)
library(tidyverse)
library(zoo)
library(reshape2)
library(colorspace)
library(sfsmisc)
source("Rscripts/baseRscript.R")


#If run previously, skip to ## **START** ##
#1. Create non-filtered Overview files
MergedPositions<-read.csv("Data/MergedPositionInfo.csv", stringsAsFactors = F)
sub<-c("1A","1B","3A")
    
for (i in 1:3){
    lname<-paste0("HCV",sub[i])
    lS<-list.files(paste0("Output",sub[i],"/Overview2/"), pattern="overview2.csv")
    assign(lname,lS)
}
    
for (g in 1:3){
    flist<-get(paste0("HCV",sub[g]))
    
    Positions<-MergedPositions[,c("merged.pos",paste0("org.pos.",sub[g]))]
    pos<-data.frame(merged.pos=Positions$merged.pos)
    for (i in 1:length(flist)){
        dat<-read.csv(paste0("Output",sub[g],"/Overview2/",flist[i]), stringsAsFactors = F, row.names = 1)
        colnames(dat)[1]<-paste0("org.pos.",sub[g])
        
        dat2<-merge(Positions,dat, by=paste0("org.pos.",sub[g]),all.x=T )
        dat3<-merge(pos, dat2, by="merged.pos", all.x=T)
        
        fname<-substr(paste(flist[i]),start=1,stop=7)
        write.csv(dat3,paste0("Output_all/Overview",sub[g],"/Overview2.2/", fname, "_overview2.2.csv"))
    }
}
  
    
#2. Create transition mutation summary of Overview2.2 (unfiltered) 
#2.1 Create merged metadata since some are missing from "Output_all/Ts_summary_metadata.1A.csv"

Metadata<-list()
for (i in 1:3){
        flist<-list.files(paste0("Output_all/Overview",sub[i],"/Overview2.2/"), pattern="overview2.2.csv")
        meta1<-read.csv(paste0("Output_all/",sub[i],".ref.overview.csv"), stringsAsFactors =F, row.names = 1)
        colnames(meta1)[2:8]<-paste0(colnames(meta1)[2:8],".",sub[i])
        colnames(meta1)[1]<-paste0("org.pos.", sub[i])
        meta2<-read.csv(paste0("Output_all/Overview",sub[i],"/Overview2.2/", flist[1]), stringsAsFactors =F, row.names = 1)
        meta<-merge(meta2[,c(1,2)], meta1, by=paste0("org.pos.", sub[i]), all=T)
        meta<-meta[,c(2,1,3:ncol(meta))]
        meta<-meta[order(meta$merged.pos),]
        
        Metadata[[i]]<-meta
        }

merged.meta<-do.call(cbind, Metadata)
merged.meta<-merged.meta[,-c(10,19)]

genes<-read.csv("Data/HCV_annotations_joined.csv", stringsAsFactors = F)
gene.vector<-c()
for (i in 1:(nrow(genes)-1)){
        gene.vector<-c(gene.vector, rep(genes$Gene[i],times=genes$start[i+1]-genes$start[i]))
}

genetable<-data.frame("merged.pos"=c(1:length(gene.vector)))
genetable$gene<-gene.vector
merged.meta<-merge(genetable,merged.meta, by="merged.pos")
#add codon positions
codon<-rep(1:3, times=nrow(merged.meta)/3)
codon<-c(2,3,codon)

merged.meta<-cbind(codon[1:nrow(merged.meta)],merged.meta)
merged.meta<-merged.meta[c(2,1,3:ncol(merged.meta))]
colnames(merged.meta)[2]<-"codon"

write.csv(merged.meta, "Output_all/merged.metadata.csv")

#############################################################
#merged.meta<-read.csv("Output_all/merged.metadata.csv",stringsAsFactors = F,row.names = 1)

## 1) combine all "minor variant" frequency tables regardless of Maj==Ref or not

Trans<-list()
Transv<-list()
Allmf<-list()
for (f in 1:3){
        flist<-list.files(paste0("Output_all/Overview",sub[f],"/Overview2.2/"),pattern="overview2.2.csv")
        
        MutFreq_Ts<-list()
        MutFreq_all<-list()
        MutFreq_Tvs<-list()
        
        for (i in 1:length(flist)){
                dat<-read.csv(paste0("Output_all/Overview",sub[f],"/Overview2.2/",flist[i]), stringsAsFactors = F, row.names=1)
                #filter out the sites with read depth <1000
                low_reads<-which(dat$TotalReads<1000) 
                dat[low_reads,c(7:ncol(dat))]<-NA
                
                #get the Ts. mutation frequency
                MutFreq_Ts[[i]]<-dat[,c("merged.pos","freq.Ts")] 
                MutFreq_Tvs[[i]]<-dat[,c("merged.pos","freq.transv")] 
                MutFreq_all[[i]]<-dat[,c("merged.pos","freq.mutations")] 
                
                filename<-substr(paste0(flist[i]),start=1,stop=7)
                names(MutFreq_Ts)[i]<-filename
                names(MutFreq_all)[i]<-filename
                names(MutFreq_Tvs)[i]<-filename
        }
        
        Ts<-MutFreq_Ts%>% purrr::reduce(full_join, by='merged.pos')
        Ts$mean<-rowMeans(Ts[2:ncol(Ts)],na.rm=T)
        
        Tvs<-MutFreq_Tvs%>% purrr::reduce(full_join, by='merged.pos')
        Tvs$mean<-rowMeans(Tvs[2:ncol(Tvs)],na.rm=T)
        
        all<-MutFreq_all%>% purrr::reduce(full_join, by='merged.pos')
        all$mean<-rowMeans(all[2:ncol(all)],na.rm=T)
        
        #Replace the sites that did not have at least 1/3 samples with NA
        s<-length(flist)
        Ts2<-Ts[Ts$merged.pos>=264,]
        Ts2$sum<-apply(Ts2[2:(s+1)],1,function(x) sum(!is.na(x)))
        Ts2$keep0.3<-(Ts2$sum/s)>=1/3
        Ts2$mean[ Ts2$keep0.3==F]<-NA
        
        Tvs2<-Tvs[Tvs$merged.pos>=264,]
        Tvs2$sum<-apply(Tvs2[2:(s+1)],1,function(x) sum(!is.na(x)))
        Tvs2$keep0.3 <-(Tvs2$sum/s)>=1/3
        Tvs2$mean[Tvs2$keep0.3==F]<-NA
        
        All2<-all[all$merged.pos>=264,]
        All2$sum<-apply(All2[2:(s+1)],1,function(x) sum(!is.na(x)))
        All2$keep0.3 <-(All2$sum/s)>=1/3
        All2$mean[All2$keep0.3==F]<-NA

         Ts2<- Ts2[,c("merged.pos","mean")]
        Tvs2<-Tvs2[,c("merged.pos","mean")]
        All2<-All2[,c("merged.pos","mean")]
        
        colnames( Ts2)[2]<-paste0("mean.",sub[f]) 
        colnames(Tvs2)[2]<-paste0("mean.",sub[f]) 
        colnames(All2)[2]<-paste0("mean.",sub[f]) 
        
        Trans[[f]]<-Ts2
        Transv[[f]]<-Tvs2
        Allmf[[f]]<-All2
}

#Combine the mean mut freq of 3 subtypes & meta data into 1 table.
Summary<-merge(merged.meta,Trans[[1]],by="merged.pos", all.x = T)
for (f in 2:3){
    Summary<-merge(Summary, Trans[[f]],by="merged.pos", all.x = T )
}
write.csv(Summary,"Output_all/MVF/Ts.MinorVariant_Mean_3subtypes.csv")

#Combine all transv freq
Summary<-merge(merged.meta,Transv[[1]],by="merged.pos", all.x = T)
for (f in 2:3){
    Summary<-merge(Summary, Transv[[f]],by="merged.pos", all.x = T )
}
write.csv(Summary,"Output_all/MVF/Tvs.MinorVariant_Mean_3subtypes.csv")

### Combined all MVF 3 subtypes
Summary<-merge(merged.meta,Allmf[[1]],by="merged.pos", all.x = T)
for (f in 2:3){
    Summary<-merge(Summary, Allmf[[f]],by="merged.pos", all.x = T )
}
write.csv(Summary,"Output_all/MVF/All.MinorVariant_Mean_3subtypes.csv")


#S<-read.csv("Output_all/MVF/Tvs.MinorVariant_Mean_3subtypes.csv")
#a1<-S[!is.na(S$mean.1A),]
#b1<-S[!is.na(S$mean.1B),]
#a3<-S[!is.na(S$mean.3A),]
#################################
#2) separate the sites that are Maj=Ref vs. Maj!=Ref for comparing the mutation types
#merged.meta<-read.csv("Output_all/merged.metadata.csv",stringsAsFactors = F,row.names = 1)
subs<-c("1A","1B","3A")

T.same<-list()
Transv<-list()
Allmf<-list()
for (f in 1:3){
        flist<-list.files(paste0("Output_all/Overview",subs[f],"/Overview2.2/"),pattern="overview2.2.csv")
        
        MutFreq_Ts.same<-list()
        MutFreq_Tvs.same<-list()
        MutFreq_all.same<-list()
        
        for (i in 1:length(flist)){
                dat<-read.csv(paste0("Output_all/Overview",subs[f],"/Overview2.2/",flist[i]), stringsAsFactors = F, row.names = 1)
                filename<-substr(paste0(flist[i]),start=1,stop=7)
                
                #remove the both ends of the genome without info
                dat<-dat[dat$merged.pos>261 & dat$merged.pos<8700,]
                
                #filter out the sites with read depth <1000
                low_reads<-which(dat$TotalReads<1000) 
                dat[low_reads,c(9:20)]<-NA
                
                #Use the sites with maj=ref
                remove<-which(dat$MajNt!=paste0(dat$ref))
                dat[remove, 9:20]<-NA
      
                MutFreq_Ts.same[[i]]<-dat[,c("merged.pos","freq.Ts.ref")] 
                names(MutFreq_Ts.same)[i]<-filename
                
                MutFreq_Tvs.same[[i]]<-dat[,c("merged.pos","freq.transv.ref")] 
                MutFreq_all.same[[i]]<-dat[,c("merged.pos","freq.mutations.ref")] 
                names(MutFreq_Tvs.same)[i]<-filename
                names(MutFreq_all.same)[i]<-filename
                
        }
         
        Ts.same<-MutFreq_Ts.same %>% purrr::reduce(full_join, by='merged.pos')
        Ts.same$mean<-rowMeans(Ts.same[2:ncol(Ts.same)],na.rm=T)
        
        Tvs.same<-MutFreq_Tvs.same %>% purrr::reduce(full_join, by='merged.pos')
        Tvs.same$mean<-rowMeans(Tvs.same[2:ncol(Tvs.same)],na.rm=T)
        
        All.same<-MutFreq_all.same %>% purrr::reduce(full_join, by='merged.pos')
        All.same$mean<-rowMeans(All.same[2:ncol(All.same)],na.rm=T)
        
        
        s<-length(flist)
        Ts.same$sum<-apply(Ts.same[2:(s+1)],1,function(x) sum(!is.na(x)))
        Ts.same$keep0.3<-(Ts.same$sum/s)>=1/3
        
        Tvs.same$sum<-apply(Tvs.same[2:(s+1)],1,function(x) sum(!is.na(x)))
        Tvs.same$keep0.3<-(Tvs.same$sum/s)>=1/3
        
        All.same$sum<-apply(All.same[2:(s+1)],1,function(x) sum(!is.na(x)))
        All.same$keep0.3<-(All.same$sum/s)>=1/3
        
        Ts.same$mean[Ts.same$keep0.3==F]<-NA
        Tvs.same$mean[Tvs.same$keep0.3==F]<-NA
        All.same$mean[All.same$keep0.3==F]<-NA
        
         Ts.same<- Ts.same[,c("merged.pos","mean")]
        Tvs.same<-Tvs.same[,c("merged.pos","mean")]
        All.same<-All.same[,c("merged.pos","mean")]
        
        colnames( Ts.same)[2]<-paste0("mean.",subs[f]) 
        colnames(Tvs.same)[2]<-paste0("mean.",subs[f]) 
        colnames(All.same)[2]<-paste0("mean.",subs[f]) 

        T.same[[f]]<-Ts.same
        Transv[[f]]<-Tvs.same
        Allmf[[f]]<-All.same

}

#Combine the mean mut freq of 3 subtypes into 1 table.
merged2<-merged.meta[merged.meta$merged.pos>263&merged.meta$merged.pos<8700,]
Summary<-merge(merged2,T.same[[1]],by="merged.pos", all.x = T)
for (f in 2:3){
    Summary<-merge(Summary, T.same[[f]],by="merged.pos", all.x = T )
}
write.csv(Summary,"Output_all/MVF/Ts.Same_Mean_3subtypes.csv")

Summary<-merge(merged2,Transv[[1]],by="merged.pos", all.x = T)
for (f in 2:3){
    Summary<-merge(Summary, Transv[[f]],by="merged.pos", all.x = T )
}
write.csv(Summary,"Output_all/MVF/Tvs.Same_Mean_3subtypes.csv")

Summary.All<-merge(merged2,Allmf[[1]],by="merged.pos", all.x = T)
for (f in 2:3){
    Summary.All<-merge(Summary.All,Allmf[[f]],by="merged.pos", all.x = T )
}
write.csv(Summary.All,"Output_all/MVF/All.Same_Mean_3subtypes.csv")


#############################################################
#############################################################

#Create mean and sd summary data tables
Summary<-read.csv("Output_all/MVF/Ts.MinorVariant_Mean_3subtypes.csv",stringsAsFactors = F, row.names = 1)

Sum2<-read.csv("Output_all/MVF/Tvs.MinorVariant_Mean_3subtypes.csv",stringsAsFactors = F, row.names = 1)
Sum3<-read.csv("Output_all/MVF/All.MinorVariant_Mean_3subtypes.csv",stringsAsFactors = F, row.names = 1)

Sum2<-Sum2[,c("merged.pos","mean.1A", "mean.1B", "mean.3A")]
colnames(Sum2)[2:4]<-c("mean.tvs.1A", "mean.tvs.1B", "mean.tvs.3A")
Sum3<-Sum3[,c("merged.pos","mean.1A", "mean.1B", "mean.3A")]
colnames(Sum3)[2:4]<-c("mean.all.1A", "mean.all.1B", "mean.all.3A")

SumT1<-Summary[,c("merged.pos", "mean.1A","mean.1B", "mean.3A")]
SumT<-merge(SumT1, Sum2, by="merged.pos")
SumT<-merge(SumT, Sum3, by="merged.pos")

write.csv(SumT,"Output_all/MVF/SummaryAll.csv")
###### Get the depth info for overview2.2 ######
## Total Reads
for (g in 1:3){
        sb<-sub[g] 
        SeqDt<-list.files(paste0("Output_all/Overview",sb, "/Overview2.2/"),pattern=".csv")
        Reads<-list()
        for (i in 1:length(SeqDt)){   
                id<-substr(paste(SeqDt[i]),start=1,stop=7)
                print(id)
                DF<-read.csv(paste0("Output_all/Overview",sb, "/Overview2.2/",SeqDt[i]),stringsAsFactors=FALSE, row.names = 1)
                DF$TotalReads[DF$TotalReads<1000]<-NA
                Reads[[i]]<-DF[,c("merged.pos","TotalReads")]
                names(Reads)[i]<-id
                
        }
        for (i in 1:length(Reads)) {
                colnames(Reads[[i]])<-c("merged.pos",paste0(names(Reads[i])))
        }
        
        Reads_Total<-Reads%>% purrr::reduce(full_join, by='merged.pos')
        write.csv(Reads_Total, paste0("Output_all/MVF/Total_Reads.",sb, ".csv"))
}


#attach the depth info:
positions<-Summary[,c("merged.pos", "org.pos.1A","org.pos.1B","org.pos.3A")]
for (g in 1:3){
        sb<-sub[g]
        dt1<-read.csv(paste0("Output_all/MVF/Total_Reads.",sb,".csv"), stringsAsFactors = F, row.names = 1)
        dt1$Depth<-rowSums(dt1[2:ncol(dt1)], na.rm=T)
        depth<-dt1[,c("merged.pos","Depth")]
        colnames(depth)[2]<-paste0("depth.",sb)
        Summary<-merge(Summary, depth, by="merged.pos")
}

write.csv(Summary,"Output_all/MVF/Ts.MV_Mean_3subtypes_depth.csv")


SumT<-merge(SumT, Summary[,c("merged.pos","depth.1A","depth.1B","depth.3A")], by="merged.pos")
write.csv(SumT,"Output_all/MVF/MVF.all.Ts.Tvs.depth.csv")

## create a mean/sd summary table:
tb<-data.frame(type=colnames(SumT)[2:10])
tb$Subtype<-rep(c("1A","1B","3A"),times=3)
tb$Type<-c(rep("Transition",times=3), rep("Transversion",times=3), rep("Total MV",times=3))

for (i in 1:9){
        if (i%%3==1) gty="1A"
        if (i%%3==2) gty="1B"
        if (i%%3==0) gty="3A"
        n<-i+1
        tb$Mean[i]<-mean(SumT[,n], na.rm=T)
        tb$SE[i]<-sqrt(tb$Mean[i]*(1-tb$Mean[i])/mean(SumT[,paste0("depth.",gty)], na.rm=T)) 
}
write.csv(tb, "Output_all/MVF/Summary_3subtypes.csv")


