#Quality filter the sam files mapped to each sample's consensus.
#Reads are mapped separately for merged (me) and unmerged (un) files
#Sam files are avilable at Figshare.

#mapping is done for 
library(stringr)
library(ape)
library(seqinr)
library(e1071)
library(msa)
library(car)
library(readtext)

suppressMessages(library(msa))

#Specify the subtype: 
subtypes<-c("1A", "1B", "3A")

for (g in 1:3){
    #specify the subtype
    sub<-subtypes[g]
    
    #list the sam file names (sam files are not availabel at github due to the large size)
    HCVsams<-list.files(paste0("Output",sub,"/sam/"),recursive = F,pattern="sam")
    
    coln<-c('QNAME','Flag','RefName','Pos','MapQ','cigar','MRNM','Mpos','isize','seq','Qual','tag1','tag2','tag3','tag4','tag5')
    
    ham.distance<-list()
    filenames<-list()
    sam_indels<-list()
    ham.distance.indel<-list()
    filtered_sam<-list()
    sum<-list()
    
    for (i in 1:length(HCVsams)){
        
        print(i)
        print(HCVsams[i])
        
        #read the sam file
        sam<-read.table(paste0("Output",sub,"/sam/",HCVsams[i]),skip=3, col.names=coln, sep = "\t",fill=T, comment.char="",quote= "")
        #remove the comment columns not needed for the analysis
        sam<-sam[,1:11]
        #select only the mapped sequences
        sam<-subset(sam, MapQ>0&MapQ<61)
        sam$seq<-as.character(sam$seq)
        print(paste("# of mapped reads: ",nrow(sam)))
        mapped.reads<-nrow(sam)
        file.name<-substr(paste(HCVsams[i]),start=1,stop=7)
        fname2<-substr(paste(HCVsams[i]),start=1,stop=10)
        
        #plot the mappling quality as a histogram 
        hist(sam$MapQ,main=paste('Histogram of Mapping Quliaty for ',fname2))
        
        #read the consesus sequence for each sample
        consensus<-read.dna(paste0("Output",sub,"/Consensus/",file.name,"_consensus.fasta"), format = "fasta",as.character=TRUE)
        
        #cheange characters in the consensus sequence to the upper case
        consensus<-as.character(sapply(consensus, function(v) {
            if (is.character(v)) return(toupper(v))
            else return(v)
        }))
        #replace ? with N
        consensus<-recode(consensus,'"?"="N"')
        
        #remove rows containing indels as it substantially increases the hamming distance
        cigars<-as.character(sam$cigar)
        indel<-grepl("I|D|N",cigars)
        sam2<-sam[!indel,]  #reads without indels
        print(paste("# of reads w/o indels: ",nrow(sam2)))
        no_indels<-nrow(sam2)
        
        #create a list with reads with indels only
        sam_indels<-sam[indel,]
        with_indels<-nrow(sam_indels)
        
        #calculate the hamming distance of each read (without indels) to its consensus
        H<-c()
        for (j in 1:nrow(sam2)) {
            reads<-sam2$seq[j]
            reads<-strsplit(reads,"")
            reads<-reads[[1]]
            if ((sam2$Pos[j]+length(reads)-1)>length(consensus)) {
                cut<-(sam2$Pos[j]+length(reads)-1)-length(consensus)
                reads<-reads[1:(length(reads)-cut)]
                ref=consensus[sam2$Pos[j]:length(consensus)]
            }       else {ref<-consensus[sam2$Pos[j]:(sam2$Pos[j]+length(reads)-1)]
            }
            m<-cbind(ref,reads)
            H[j]<-hamming.distance(m[,1],m[,2])
        }
        
        #Plot the results 
        filename<-paste0("Output",sub,"/HammingDistanceFiltering/",fname2,".pdf")
        pdf(filename, width =10, height = 5)
        par(mfrow=c(1,2))
        par(mar = c(5,4,4,2))
        plot(table(H),main=paste(fname2),xlab='Hamming distance',ylab='Counts')
        
        ham.distance[[i]]<-H
        names(ham.distance[i])<-paste(fname2)
        filenames[[i]]<-paste(fname2)
        
        #reads with indels
        Hindel<-c()
        for (k in 1:nrow(sam_indels)){
            read<-sam_indels$seq[k]
            read.s<-strsplit(read,"")
            read.s<-read.s[[1]]
            if ((sam_indels$Pos[k]+length(read.s)-1)>length(consensus)) {
                cut<-(sam_indels$Pos[k]+length(read.s)-1)-length(consensus)
                read.s<-read.s[1:(length(read.s)-cut)]
                ref=consensus[sam_indels$Pos[k]:length(consensus)]
            }        else {ref<-consensus[sam_indels$Pos[k]:(sam_indels$Pos[k]+length(read.s)-1)]
            }
            
            ref<-paste(ref, sep="", collapse="")
            dna<-DNAStringSet(c(paste(ref),paste(read)))
            names(dna)<-c(paste0('ref',k),paste0(file.name,".read",paste(k)))
            invisible(capture.output(align<-msa(dna)))
            aligned<-DNAStringSet(unlist(align))
            m2<-cbind(paste(aligned[1]),paste(aligned[2]))
            m2<-strsplit(m2,"")
            m2.m<-do.call(cbind,m2)
            
            Hindel[k]<-hamming.distance(m2.m[,1],m2.m[,2])-1
        }
        plot(table(Hindel),main=paste(fname2,"w/ indels"),xlab='Hamming distance',ylab='Counts')
        ham.distance.indel[[i]]<-Hindel
        names(ham.distance.indel[i])<-paste(fname2)
        
        
        #eliminate the reads with hamming distance>10 and save as a new sam file
        Large.ham<-(H>10)
        sam_re<-sam2[which(Large.ham==F),]
        no_indels_removed<-length(Large.ham[Large.ham==T])
        print(paste("# of removed reads w/o indels: ",no_indels_removed))
        
        Large.ham2<-(Hindel>10)
        sam_re2<-sam_indels[which(Large.ham2==F),]
        with_indels_removed<-length(Large.ham2[Large.ham2==T])
        print(paste("# of removed reads with indels: ",with_indels_removed))
        sam_RC<-rbind(sam_re,sam_re2)
        write.table(sam_RC, paste0("Output",sub,"/sam2/", fname2,"-filtered.sam"),sep="\t", quote=F,row.names=F,col.names=F)
        
        sum[[i]]<-data.frame(fname2,mapped.reads,no_indels, with_indels,no_indels_removed,with_indels_removed)
        dev.off()
    }
    
}






