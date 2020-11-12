#Create the Overview dataframe of mutation freqeuncies, types of mutations etc.
# in ='/SeqData/SeqData_xxxxxx.csv', out = '/Overview2/xxxx_overview2.csv'


library(dplyr)
source("Rscripts/baseRscript.R")

#Specify the Subtype
subtypes<-c("1A", "1B", "3A")

for (g in 1:2) {
    sub<-subtypes[g]

	#get the file name
	HCVFiles_SeqData<-list.files(paste0("Output",sub,"/SeqData/"),pattern="SeqData")
    
	df<-read.csv(paste0("Output",sub,"/SeqData/",HCVFiles_SeqData[1]),stringsAsFactors=FALSE, row.names = 1)
    df<-df[,c("pos","ref")]
    df$transition.ref<-sapply(df$ref,function(x) transition(x))
    
    #Add information whether the mutation at a particular site is Syn, nonSyn, or Stop codon.        
    TypeOfSite<-c()
    for (codon in 1:(nrow(df)/3)) {#for each codon in the sequence
            positions <- c(codon*3-2,codon*3-1, codon*3)
            WTcodon <- df$ref[positions]  
            mutant1codon <- c(transition(WTcodon[1]), WTcodon[2:3])  #If the first position has transistion mutation, it's labeld as mutatnt1codon.
            mutant2codon <- c(WTcodon[1],transition(WTcodon[2]), WTcodon[3])
            mutant3codon <- c(WTcodon[1:2], transition(WTcodon[3]))
    
            TypeOfSite<-c(TypeOfSite,typeofsitefunction(WTcodon,mutant1codon))
            TypeOfSite<-c(TypeOfSite,typeofsitefunction(WTcodon,mutant2codon))
            TypeOfSite<-c(TypeOfSite,typeofsitefunction(WTcodon,mutant3codon))
                
    }
        
    df$Type<-TypeOfSite[1:nrow(df)]
    
    for (k in 1:nrow(df)){
            if (k%%3==1){
                df$WTAA[k] = seqinr::translate(df$ref[c(k,k+1,k+2)])
                df$MUTAA[k] = seqinr::translate(c(transition(df$ref[k]),df$ref[c(k+1,k+2)]))
                        } 
            if (k%%3==2){
                    df$WTAA[k] = seqinr::translate(df$ref[c(k-1,k,k+1)])
                    df$MUTAA[k] = seqinr::translate(c(df$ref[c(k-1)],transition(df$ref[k]),df$ref[c(k+1)]))
                }
            if (k%%3==0){
                    df$WTAA[k] = seqinr::translate(df$ref[c(k-2,k-1,k)])
                    df$MUTAA[k] = seqinr::translate(c(df$ref[c(k-2,k-1)],transition(df$ref[k])))
            }
            
    }

        #Add whether AA change is drastic & makes a CpG site
        df$bigAAChange<-0
        df$makesCpG <- 0
        
        for (j in 2:(nrow(df)-1)){
                WT <- amCat(df[j,'WTAA'])
                MUT <-amCat(df[j,'MUTAA'])
                
                if (WT != MUT)  df$bigAAChange[j] <- 1
                
                trip <- df$ref[c(j-1, j,j+1)]
                if (is.na(trip[1])|is.na(trip[2])|is.na(trip[3])) 
                        next
                else{
                        if (trip[1] == "c" & trip[2] == "a" ) df$makesCpG[j] <- 1 
                        if (trip[2] == "t" & trip[3] == "g")  df$makesCpG[j] <- 1
                    }
        }
        
        
        write.csv(df,paste0("Output_all/",sub, ".ref.overview.csv"))
	}        
