#Calculate mutation frequencies at each site from the frequency table file created in Step 2 (2.Pileup_Rsamtools.R)
# in ='/CSV/xxxx_un.csv' and '/CSV/xxxx_me.csv', out='/SeqData/SeqData_xxxxxx.csv'

library(tidyverse)
library(plyr)
source("Rscripts/baseRscript.R")

#Specify the Subtype
subtypes<-c("1A", "1B", "3A")

for (g in 1:3) {
    sub<-subtypes[g]
    
    # read the CSV files (unmerged + merged separately) from the step2
    HCVFiles<-list.files(paste0("Output",sub,"/CSV/"),pattern="un.csv")
    MergeFiles<-list.files(paste0("Output",sub,"/CSV/"), pattern="me.csv")

    if (g==3) start<-262
    else {start<- 264}
     
    no<-data.frame("pos"=c(start:8800))
    
    for (i in 1:length(MergeFiles)){
            print(i)
            id<-substr(paste(HCVFiles[i]),start=1,stop=7)
            print(id)
            
            #Read in CSV files for merged and unmerged and combined them
            merge<-read.csv(paste("Output",sub,"/CSV/",MergeFiles[i],sep=""))
            merge<-merge[,-c(1,2,8,9)]
            colnames(merge)<-c("pos","mA","mC","mG","mT", "mDel","mIns")
            unmerge<-read.csv(paste("Output",sub,"/CSV/",HCVFiles[i],sep=""))
            unmerge<-unmerge[,-c(1,2,8,9)]
            colnames(unmerge)<-c("pos","uA","uC","uG","uT", "uDel","uIns")
    
            SeqData<-join(merge,unmerge,by="pos")
            SeqData[is.na(SeqData)] <- 0
            SeqData$A<-SeqData$mA+SeqData$uA
            SeqData$C<-SeqData$mC+SeqData$uC
            SeqData$G<-SeqData$mG+SeqData$uG
            SeqData$t<-SeqData$mT+SeqData$uT
            SeqData$Del<-SeqData$mDel+SeqData$uDel
            SeqData$Ins<-SeqData$mIns+SeqData$uIns
            SeqData<-SeqData[,c("pos","A","C","G","t","Del","Ins")]
            SeqData$TotalReads<-rowSums(SeqData[2:5],na.rm=T)
            SeqData$TotalReads_indels<-rowSums(SeqData[2:7],na.rm=T)
            colnames(SeqData)[2:7]<-c("a","c","g","t","deletion","insertion")
            
            #determine the majority nucleotide base at each site
            SeqData$MajNt<-apply(SeqData[,2:5],1,function(x) c("a","c","g","t")[which.max(x)])
            
            #read the refrence sequence:
            reference<-read.dna(paste0("Data/HCV",sub,"_Consensus.fasta"), format = "fasta",as.character=TRUE)
            ref.code<-reference[start:8800]
            SeqData<-merge(no,SeqData,by="pos",all.x=T)
            SeqData$ref<-ref.code[1:length(SeqData[,1])]
            
            #nuceotides from transition mutations
            SeqData$transition.maj[j]<-sapply(SeqData$MajNt, function(x) transition(x))        
            SeqData$transition.ref[j]<-sapply(SeqData$ref, function(x) transition(x))
        
            #rearrange the columns
            SeqData<-SeqData[,c("a","c","g","t","deletion","insertion","pos","TotalReads","TotalReads_indels","MajNt","ref","transition.maj","transition.ref")]
            
            #determine mutation frequencies of every site.
            for (k in 1:nrow(SeqData)){
                    if (is.na(SeqData$MajNt[k])) {
                            SeqData$freq.Ts[k]<-NA #transition mutations (minor varient freq)
                            SeqData$freq.Ts.ref[k]<-NA # transition mutations in realtive to the reference seq
                            
                            SeqData$freq.transv[k]<-NA #transversion mutations
                            SeqData$freq.transv.ref[k]<-NA
                            SeqData$freq.transv1[k]<-NA
                            SeqData$freq.transv2[k]<-NA
                            SeqData$freq.transv1.ref[k]<-NA
                            SeqData$freq.transv2.ref[k]<-NA
                            
                            SeqData$freq.mutations.ref[k]<-NA #all mutations
                            SeqData$freq.mutations[k]<-NA
                            SeqData$freq.mutations.indels[k]<-NA
                            SeqData$freq.mutations.indels.ref[k]<-NA
                    }
                    else {
                            MajNum <- SeqData [k,paste0(SeqData$MajNt[k])]
                            MutNum1<- SeqData [k,paste0(SeqData$transition.maj[k])]
                            WTNum <- SeqData [k,paste0(SeqData$ref[k])]
                            MutNum2<- SeqData [k,paste0(SeqData$transition.ref[k])]
                            
                            SeqData$freq.Ts[k]<-MutNum1/SeqData$TotalReads[k]
                            SeqData$freq.Ts.ref[k]<-MutNum2/SeqData$TotalReads[k]
                            
                            
                            #mutation frequencies of all ransversion mutataions
                            if (SeqData$MajNt[k]=="a"|SeqData$MajNt[k]=='g'){
                                    TrvMutNum<-SeqData[k,"c"]+SeqData[k,"t"]}
                            if (SeqData$MajNt[k]=="c"|SeqData$MajNt[k]=="t"){
                                    TrvMutNum<-SeqData[k,"a"]+SeqData[k,"g"]}
                            SeqData$freq.transv[k]<-TrvMutNum/SeqData$TotalReads[k]
                            if (SeqData$ref[k]=="a"|SeqData$ref[k]=='g'){
                                    TrvMutNum2<-SeqData[k,"c"]+SeqData[k,"t"]}
                            if (SeqData$ref[k]=="c"|SeqData$ref[k]=="t"){
                                    TrvMutNum2<-SeqData[k,"a"]+SeqData[k,"g"]}
                            SeqData$freq.transv.ref[k]<-TrvMutNum2/SeqData$TotalReads[k]
                            
                            #Frequenceis for specific transversion mutations (1 & 2)
                            Tvs1Num<-SeqData[k,paste0(transv1(SeqData$MajNt[k]))]
                            Tvs2Num<-SeqData[k,paste0(transv2(SeqData$MajNt[k]))]
                            SeqData$freq.transv1[k]<-Tvs1Num/SeqData$TotalReads[k]
                            SeqData$freq.transv2[k]<-Tvs2Num/SeqData$TotalReads[k]
                            Tvs1rNum<-SeqData[k,paste0(transv1(SeqData$ref[k]))]
                            Tvs2rNum<-SeqData[k,paste0(transv2(SeqData$ref[k]))]
                            SeqData$freq.transv1.ref[k]<-Tvs1rNum/SeqData$TotalReads[k]
                            SeqData$freq.transv2.ref[k]<-Tvs2rNum/SeqData$TotalReads[k]
                            
                            
                            #Frequencies of all SNPs (no indels)
                            AllMutNum<-SeqData$TotalReads[k]-MajNum
                            AllMutNum2<-SeqData$TotalReads[k]-WTNum
                            
                            SeqData$freq.mutations[k]<-AllMutNum/SeqData$TotalReads[k]
                            SeqData$freq.mutations.ref[k]<-AllMutNum2/SeqData$TotalReads[k]
                            
                            #with indels
                            SeqData$freq.mutations.indels[k]<-AllMutNum/SeqData$TotalReads_indels[k]
                            SeqData$freq.mutations.indels.ref[k]<-AllMutNum2/SeqData$TotalReads_indels[k]
                        }
                }
            
            write.csv(SeqData,paste0("Output",sub,"/SeqData/SeqData_",id,".csv"))
    s}
            
        
}
