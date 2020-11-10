#Create a frequency table using Rsamtools from a bam file.
# After the step1 read_check, convert the sam files to bam files, sort and index them using samtools.
# all bam files should be saved in the respective subtype(=XX)'s  "Output_XX/bam2/" directory 
# bam files are not available at github due to their large sizes
# in ='.bam' and '.bam.bai', out='.csv' (a freqeuncy table)

library(Rsamtools)
library(stringr)
source("Rscripts/pileupFreq.R")

subtypes<-c("1A", "1B", "3A")


#Adjust the max depth parameter based on data files
for (g in 1:3) {
    sub<-subtypes[g]
    
    #list of bam files to be processed
    bamfiles<-list.files(paste0("Output",sub,"/bam2/"),pattern="bam$")
    
    for (i in 1:length(bamfiles)){
        bam<-bamfiles[i]
        index<-paste0(paste0("Output",sub,"/bam2/",bam),'.bai')
        bf<-BamFile(paste0("Output",subo,"/bam2/",bam), index=index)
        
        file.name<-paste(bam)
        file.name<-substr(file.name,start=1,stop=10 )
        p_param <- PileupParam(max_depth=90000,distinguish_strands=FALSE,include_insertions=TRUE)
        result<-pileup(bf, pileupParam = p_param)
        summary<-pileupFreq(result)
        
        print(file.name)
        write.csv(summary, file=paste0("Output",sub,"/CSV/",file.name,".csv",collapse=""))
        
    }
}