#Calculate nucleotide diversity values based on Nelson et al. (2015) method.
library(purrr)
library(dplyr)
source("Rscripts/baseRscript.R")
subs<-c("1A","1B","3A")

#dir.create("Output3A/Overview_D2/")
#dir.create("Output1B/Overview_D2/")
#dir.create("Output1A/Overview_D2/")


for (f in 1:3){
        flist1<-list.files(paste0("Output",subs[f],"/Overview2/"),pattern="overview2.csv")
        flist2<-list.files(paste0("Output",subs[f],"/SeqData/"),pattern="SeqData")
        
        Overview1<-list()
        for (i in 1:length(flist1)){ 
                overview<-read.csv(paste0("Output",subs[f],"/Overview2/",flist1[i]),stringsAsFactors=FALSE, row.names=1)
                #transition info only
                overview<-overview[,c("pos","TotalReads","MajNt","ref","freq.Ts.ref","Type","WTAA","MUTAA","bigAAChange","makesCpG")]
                low_reads<-which(overview$TotalReads<1000) 
                overview[low_reads,"freq.Ts.ref"]<-NA
                #remove maj!=ref sites
                remove<-which(overview$ref!=overview$MajNt)
                overview[remove,"freq.Ts.ref"]<-NA
                filename<-substr(paste(flist1[i]),start=1,stop=6)
                
                seq<-read.csv(paste0("Output",subs[f],"/SeqData/",flist2[i]),stringsAsFactors=FALSE, row.names=1)
                seq<-seq[,c(7,1:4)]
                dat<-merge(seq, overview, by="pos", all.y=T )
                dat$D<-""
                for (k in 1:nrow(dat)){
                        if (is.na(dat$freq.Ts.ref[k]))  dat$D[k]<-NA
                        else{
                                ac<-dat[k, "a"]*dat[k,"c"]
                                ag<-dat[k, "a"]*dat[k,"g"]
                                at<-dat[k, "a"]*dat[k,"t"]
                                cg<-dat[k, "c"]*dat[k,"g"]
                                ct<-dat[k, "c"]*dat[k,"t"]
                                gt<-dat[k, "g"]*dat[k,"t"]
                                m<-(dat[k,"TotalReads"]^2-dat[k,"TotalReads"])/2
                                dat$D[k]<-(ac+ag+at+cg+ct+gt)/m
                        }
                }
                
                Overview1[[i]]<-dat
                names(Overview1)[i]<-filename
                print(paste(i,filename))
                write.csv(dat,paste0("Output",subs[f],"/Overview_D2/",filename,"_overviewD.csv"))
        }
        
        listname<-paste0("Overview1.", subs[f])
        assign(listname, Overview1)
        
        pi<-data.frame(SampleID=names(Overview1))
        for (i in 1:length(Overview1)){
                        df<-Overview1[[i]]
                        n<-nrow(df[!is.na(df$D),])
                        df$D<-as.numeric(df$D)
                        pi$Pi[i]<-(sum(df$D, na.rm=T))/n
                        
        }
        write.csv(pi, paste0("Output_all/Diversity/NucDiversity.per.sample.",subs[f],".csv"))
}


## Calculate pi (nuc diversity) per gene by genotype and plot them
### addd gene annotation info for all genes
Genes<-read.csv("Data/HCV_annotations_joined.csv", stringsAsFactors = F)
Genes$Gene<-factor(Genes$Gene, levels=c("5' UTR","Core","E1", "HVR1","E2","P7", "NS2","NS3","NS4A","NS4B","NS5A","NS5B"))

gene.vector<-c()
for (i in 1:12){
        gene.vector<-c(gene.vector, rep(paste(Genes$Gene[i]),times=Genes$start[i+1]-Genes$start[i]))
}

n<-data.frame(pos=1:length(gene.vector))
gene.vec<-cbind(n,gene.vector)
colnames(gene.vec)[2]<-"gene"
gene.vec$gene<-as.character(gene.vec$gene)

for(g in 1:3){
        flist<-list.files(paste0("Output",subs[g],"/Overview_D2/"), pattern="overviewD.csv")
        
        Dlist<-list()
        for (i in 1:length(flist)){
                dt<-read.csv(paste0("Output",subs[g],"/Overview_D2/",flist[i]), stringsAsFactors = F, row.names = 1)
                Dlist[[i]]<-dt[,c("pos","D")] 
                fname<-substr(flist[i], start = 1, stop = 6)
                names(Dlist)[i]<-fname
        }
        
        for (i in 1:length(flist)) {
                colnames(Dlist[[i]])<-c("pos",paste0(names(Dlist[i])))
        }
        
        Ds<-Dlist%>% purrr::reduce(full_join, by='pos')
        Ds<-merge(gene.vec, Ds, "pos")
        #lname<-paste0("Ds.",subs[g])
        #assign(lname, Ds)
        write.csv(Ds, paste0("Output_all/Diversity/Ds_",subs[g],".csv"))

}


SummaryPi<-data.frame()
for (g in 1:3){
        df<-read.csv(paste0("Output_all/Diversity/Ds_",subs[g],".csv"), stringsAsFactors = F, row.names = 1)
        colnames(df)[2]<-"gene"
        lname<-paste0("Ds.",subs[g])
        assign(lname, df)
        
        pi.genes<-data.frame(Gene=Genes$Gene[1:12])
        for (i in 1:nrow(pi.genes)){
               df2<-df[df$gene==Genes$Gene[i],]
               pi.genes[i, "Mean"]<-sum(df2[,3:ncol(df2)], na.rm=T)/ sum(!is.na(df2[,3:ncol(df2)]))
               pi.genes[i,"SE"]<-mean(std.error(df2[,3:ncol(df2)], na.rm=T), na.rm=T) 
        }
        pi.genes$Subtype<-subs[g]
        SummaryPi<-rbind(SummaryPi, pi.genes)
        
}

write.csv(SummaryPi, "Output_all/Diversity/Pi_byGene_bySubtype.csv")




