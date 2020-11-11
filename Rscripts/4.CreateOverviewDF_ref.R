#Create the Overview dataframe of mutation freqeuncies, types of mutations etc.
# in ='/SeqData/SeqData_xxxxxx.csv', out = '/Overview2/xxxx_overview2.csv'


library(dplyr)
source("Rscripts/baseRscript.R")

#Specify the Subtype
subtypes<-c("1A", "1B", "3A")

for (g in 1:3) {
    sub<-subtypes[g]

	#get the file name
	HCVFiles_SeqData<-list.files(paste0("Output",sub,"/SeqData/"),pattern="SeqData")

	for (i in 1:length(HCVFiles_SeqData)){   
        id<-substr(paste(HCVFiles_SeqData[i]),start=9,stop=15)
        print(id)
        OverviewDF<-read.csv(paste0("Output",sub,"/SeqData/",HCVFiles_SeqData[i]),stringsAsFactors=FALSE, row.names = 1)
       
        #Add information whether the mutation at a particular site is Syn, nonSyn, or Stop codon.        
        TypeOfSite<-c()
        TypeOfSite.tv1<-c()
        TypeOfSite.tv2<-c()
        for (codon in 1:(nrow(OverviewDF)/3)) {#for each codon in the sequence
                positions <- c(codon*3-2,codon*3-1, codon*3)
                WTcodon <- OverviewDF$ref[positions]  
                if (is.na(WTcodon[1])|is.na(WTcodon[2])|is.na(WTcodon[3])){ 
                        WTcodon<-c('n','n','n')
                        mutant1codon<-c('n','n','n')
                        mutant2codon<-c('n','n','n')
                        mutant3codon<-c('n','n','n')}
                else{                        
                        mutant1codon <- c(transition(WTcodon[1]), WTcodon[2:3])  #If the first position has transistion mutation, it's labeld as mutatnt1codon.
                        mutant2codon <- c(WTcodon[1],transition(WTcodon[2]), WTcodon[3])
                        mutant3codon <- c(WTcodon[1:2], transition(WTcodon[3]))
                        
                        #transversion mutation to 'a' or 'c'
                        mutant1codon.tv1 <- c(transv1(WTcodon[1]), WTcodon[2:3]) 
                        mutant2codon.tv1 <- c(WTcodon[1],transv1(WTcodon[2]), WTcodon[3])
                        mutant3codon.tv1 <- c(WTcodon[1:2], transv1(WTcodon[3]))
                        #transversion mutation to 'g' or 't'
                        mutant1codon.tv2 <- c(transv2(WTcodon[1]), WTcodon[2:3])  
                        mutant2codon.tv2 <- c(WTcodon[1],transv2(WTcodon[2]), WTcodon[3])
                        mutant3codon.tv2 <- c(WTcodon[1:2], transv2(WTcodon[3]))
                }
                
                
                TypeOfSite<-c(TypeOfSite,typeofsitefunction(WTcodon,mutant1codon))
                TypeOfSite<-c(TypeOfSite,typeofsitefunction(WTcodon,mutant2codon))
                TypeOfSite<-c(TypeOfSite,typeofsitefunction(WTcodon,mutant3codon))
                
                TypeOfSite.tv1<-c(TypeOfSite.tv1,typeofsitefunction(WTcodon,mutant1codon.tv1))
                TypeOfSite.tv1<-c(TypeOfSite.tv1,typeofsitefunction(WTcodon,mutant2codon.tv1))
                TypeOfSite.tv1<-c(TypeOfSite.tv1,typeofsitefunction(WTcodon,mutant3codon.tv1))
                
                TypeOfSite.tv2<-c(TypeOfSite.tv2,typeofsitefunction(WTcodon,mutant1codon.tv2))
                TypeOfSite.tv2<-c(TypeOfSite.tv2,typeofsitefunction(WTcodon,mutant2codon.tv2))
                TypeOfSite.tv2<-c(TypeOfSite.tv2,typeofsitefunction(WTcodon,mutant3codon.tv2))                
        }
        
        OverviewDF$Type<-TypeOfSite[1:length(OverviewDF$pos)]
        OverviewDF$Type.tv1<-TypeOfSite.tv1[1:length(OverviewDF$pos)]
        OverviewDF$Type.tv2<-TypeOfSite.tv2[1:length(OverviewDF$pos)]
        
        OverviewDF<-OverviewDF[,-c(1:6)]
    
        for (k in 1:nrow(OverviewDF)){
                if (k%%3==1){
                        if (is.na(OverviewDF$MajNt[k])|is.na(OverviewDF$MajNt[k+1])|is.na(OverviewDF$MajNt[k+2])) { OverviewDF$MajAA[k]<-"NA"
                        OverviewDF$WTAA[k]<-"NA"
                        OverviewDF$MUTAA[k]<-"NA"
                        OverviewDF$TVS1_AA[k]<-"NA"
                        OverviewDF$TVS2_AA[k]<-"NA"}
                        else {  OverviewDF$MajAA[k] = seqinr::translate(OverviewDF$MajNt[c(k,k+1,k+2)])
                        OverviewDF$WTAA[k] = seqinr::translate(OverviewDF$ref[c(k,k+1,k+2)])
                        OverviewDF$MUTAA[k] = seqinr::translate(c(transition(OverviewDF$ref[k]),OverviewDF$ref[c(k+1,k+2)]))
                        OverviewDF$TVS1_AA[k] = seqinr::translate(c(transv1(OverviewDF$ref[k]),OverviewDF$ref[c(k+1,k+2)]))
                        OverviewDF$TVS2_AA[k] = seqinr::translate(c(transv2(OverviewDF$ref[k]),OverviewDF$ref[c(k+1,k+2)]))}
                } 
                if (k%%3==2){
                        if (is.na(OverviewDF$MajNt[k-1])|is.na(OverviewDF$MajNt[k])|is.na(OverviewDF$MajNt[k+1]))  {OverviewDF$MajAA[k]<-"NA"
                        OverviewDF$WTAA[k]<-"NA"
                        OverviewDF$MUTAA[k]<-"NA"
                        OverviewDF$TVS1_AA[k]<-"NA"
                        OverviewDF$TVS2_AA[k]<-"NA"}
                        else {  OverviewDF$MajAA[k] = seqinr::translate(OverviewDF$MajNt[c(k-1,k,k+1)])
                        OverviewDF$WTAA[k] = seqinr::translate(OverviewDF$ref[c(k-1,k,k+1)])
                        OverviewDF$MUTAA[k] = seqinr::translate(c(OverviewDF$ref[c(k-1)],transition(OverviewDF$ref[k]),OverviewDF$ref[c(k+1)]))
                        OverviewDF$TVS1_AA[k] = seqinr::translate(c(OverviewDF$ref[c(k-1)],transv1(OverviewDF$ref[k]),OverviewDF$ref[c(k+1)]))
                        OverviewDF$TVS2_AA[k] = seqinr::translate(c(OverviewDF$ref[c(k-1)],transv2(OverviewDF$ref[k]),OverviewDF$ref[c(k+1)]))}
                }
                if (k%%3==0){
                        if (is.na(OverviewDF$MajNt[k-2])|is.na(OverviewDF$MajNt[k-1])|is.na(OverviewDF$MajNt[k]))  {  OverviewDF$MajAA[k]<-"NA"
                        OverviewDF$WTAA[k]<-"NA"
                        OverviewDF$MUTAA[k]<-"NA"
                        OverviewDF$TVS1_AA[k]<-"NA"
                        OverviewDF$TVS2_AA[k]<-"NA"}
                        else {  OverviewDF$MajAA[k] = seqinr::translate(OverviewDF$MajNt[c(k-2,k-1,k)])
                        OverviewDF$WTAA[k] = seqinr::translate(OverviewDF$ref[c(k-2,k-1,k)])
                        OverviewDF$MUTAA[k] = seqinr::translate(c(OverviewDF$ref[c(k-2,k-1)],transition(OverviewDF$ref[k])))
                        OverviewDF$TVS1_AA[k] = seqinr::translate(c(OverviewDF$ref[c(k-2,k-1)],transv1(OverviewDF$ref[k])))
                        OverviewDF$TVS2_AA[k] = seqinr::translate(c(OverviewDF$ref[c(k-2,k-1)],transv2(OverviewDF$ref[k])))}
                }
                
        }

        #Add whether AA change is drastic & makes a CpG site
        OverviewDF$bigAAChange<-0
        OverviewDF$bigAAChange.tv1<-0
        OverviewDF$bigAAChange.tv2<-0
        OverviewDF$makesCpG <- 0
        OverviewDF$makesCpG.tvs <- 0
        OverviewDF$makesCpG.tv1 <- 0
        OverviewDF$makesCpG.tv2 <- 0
        
        for(j in 2:nrow(OverviewDF)-1){
                WT <- amCat(OverviewDF[j,'WTAA'])
                MUT <- amCat(OverviewDF[j,'MUTAA'])
                MUT1<-amCat(OverviewDF[j,'TVS1_AA'])
                MUT2<-amCat(OverviewDF[j,'TVS2_AA'])
                
                if (WT != MUT) OverviewDF$bigAAChange[j] <- 1
                if (WT != MUT1) OverviewDF$bigAAChange.tv1[j] <- 1
                if (WT != MUT2) OverviewDF$bigAAChange.tv2[j] <- 1
                
                trip <- OverviewDF$ref[c(j-1, j,j+1)]
                if (is.na(trip[1])|is.na(trip[2])|is.na(trip[3])) 
                        next
                else{
                        if (trip[1] == "c" & trip[2] == "a" ) OverviewDF$makesCpG[j] <- 1 
                        if (trip[2] == "t" & trip[3] == "g")  OverviewDF$makesCpG[j] <- 1
                        if (trip[1] == "c" & (trip[2]=="c"|trip[2]=='t')) OverviewDF$makesCpG.tvs[j] <- 1
                        if (trip[3] == "g" & (trip[2]=="a"|trip[2]=="g")) OverviewDF$makesCpG.tvs[j] <- 1
                        
                        if (trip[1] == "c" & (trip[2]=="c"|trip[2]=='t')) OverviewDF$makesCpG.tv2[j] <- 1                                
                        if (trip[3] == "g" & (trip[2]=="a"|trip[2]=="g")) OverviewDF$makesCpG.tv1[j] <- 1
                    }
        }
        
        
        write.csv(OverviewDF,paste0("Output",sub,"/Overview2/",id,"overview2.csv"))
	}        


}



#################
### Read depths for all files ###
HCVFiles_SeqData<-list.files(paste0("Output",sub,"/SeqData/"),pattern="SeqData")
ReadsSummary<-data.frame(Sample_ID=matrix(nrow=length(HCVFiles_SeqData)))
ReadsSummary$MaxDepth<-""
ReadsSummary$AveDepth<-""

for (i in 1:length(HCVFiles_SeqData)){
        print(i)
        id<-substr(paste(HCVFiles_SeqData[i]),start=9,stop=15)
        ReadsSummary$Sample_ID[i]<-id
        print(id)
        SeqData<-read.csv(paste0("Output",geno,"/SeqData/",HCVFiles_SeqData[i],sep=""))
        ReadsSummary$MaxDepth[i]<-max(SeqData$TotalReads,na.rm=T)
        ReadsSummary$AveDepth[i]<-mean(SeqData$TotalReads,na.rm=T)
}

write.csv(ReadsSummary,paste0("Output",sub,"/SeqData/ReadsSummary_",sub,".csv"))      




