library(ggplot2)
library(reshape2)
library(colorspace)

source("Rscripts/baseRscript.R")
source("Rscripts/label_scientific.R")

#Geonotype colors 
colors2<-qualitative_hcl(6, palette="Dark3")
subtypes<-c("1A","1B","3A")


### create data frame of Read depth for used sites (sties > 1000 read depth) ###
ReadsSum<-list()
for (s in 1:3){
    SeqDatFiles<-list.files(paste0("Output",subtypes[s],"/SeqData/"),pattern="SeqData")
    ReadsSummary<-data.frame(Sample_ID=matrix(nrow=length(SeqDatFiles)))
    ReadsSummary$MaxDepth<-""
    ReadsSummary$AveDepth<-""
    for (i in 1:length(SeqDatFiles)){
        id<-substr(paste(SeqDatFiles[i]),start=9,stop=14)
        cat(i, id, "\n")
        ReadsSummary$Sample_ID[i]<-id
        SeqData<-read.csv(paste0("Output",subtypes[s],"/SeqData/",SeqDatFiles[i],sep=""))
        SeqData<-SeqData[SeqData$TotalReads>=1000,]
        ReadsSummary$MaxDepth[i]<-max(SeqData$TotalReads,na.rm=T)
        ReadsSummary$AveDepth[i]<-mean(SeqData$TotalReads,na.rm=T)
    }
    ReadsSummary$Subtype<-subtypes[s]
    ReadsSum[[s]]<-ReadsSummary
}

ReadDepth<-do.call(rbind, ReadsSum)

#Attach viral_load data
vl<-read.csv("Data/Viral_load.csv", stringsAsFactors = F)

SampleSummary<-merge(ReadDepth, vl[,c(1,3)], by="Sample_ID")
write.csv(SampleSummary,"Output_all/ReadDepth_Summary.csv")     

depths<-data.frame(subtype=c('all',subtypes))
for (i in 1:4){
    if (i==1){
        depths$mean[i]<-mean(SampleSummary$AveDepth, na.rm = T)
        depths$se[i]<-std.error(SampleSummary$AveDepth, na.rm = T) 
    }
    else {
        k=i-1
        depths$mean[i]<-mean(SampleSummary$AveDepth[SampleSummary$Subtype==subtypes[k]], na.rm = T)
        depths$se[i]<-std.error(SampleSummary$AveDepth[SampleSummary$Subtype==subtypes[k]], na.rm = T) 
    }
}

depths
#subtype     mean       se
#1     all 6115.048 340.7325
#2      1A 6178.412 426.7979
#3      1B 5890.685 510.3417
#4      3A 5919.040 590.0699

wilcox.test(SampleSummary$AveDepth[SampleSummary$Subtype==subtypes[1]], SampleSummary$AveDepth[SampleSummary$Subtype==subtypes[2]],alternative="greater" )
wilcox.test(SampleSummary$AveDepth[SampleSummary$Subtype==subtypes[1]], SampleSummary$AveDepth[SampleSummary$Subtype==subtypes[3]],alternative="greater" )
wilcox.test(SampleSummary$AveDepth[SampleSummary$Subtype==subtypes[3]], SampleSummary$AveDepth[SampleSummary$Subtype==subtypes[2]],alternative="greater" )
#W = 1802, p-value = 0.817
#W = 3710, p-value = 0.5952
#W = 374, p-value = 0.7101

#####################################
## correlation between diversity and viral laod/read depth 

#viral_load & depth data
vl<-read.csv("Output_all/ReadDepth_Summary.csv", row.names=1, stringsAsFactors = F)
#diversity data
Pi<-read.csv(paste0("Output_all/Diversity/NucDiversity.per.sample.1A.csv"), stringsAsFactors = F, row.names = 1)
Pi$Subtype="1A"
for(g in 2:3){
    dt<-read.csv(paste0("Output_all/Diversity/NucDiversity.per.sample.",subtypes[g],".csv"), stringsAsFactors = F, row.names = 1)
    dt$Subtype=subtypes[g]
    Pi<-rbind(Pi, dt)
}

colnames(Pi)[1]<-"Sample_ID"
div<-merge(Pi,vl[,c("Sample_ID","Viral_load","AveDepth")], by="Sample_ID")
div2<-div[!is.na(div$Viral_load),]


ggplot(data=div2, aes(x=Viral_load, y=Pi, color=Subtype))+
    geom_point()+
    scale_color_manual(values=colors2[c(1,3,5)])+
    scale_x_continuous(trans='log10',labels=label_scientific)+
    xlab("viral load")+ylab('nucleotide diversity')+
    theme_bw()
ggsave("Output_all/Diversity/Viral_load.vs.Pi.pdf", width = 5.5, height = 4)


#Test viral load vs. pi
cor.test(div2$Viral_load, div2$Pi, method="spearman")
#S = 703030, p-value = 0.9215
#rho =0.007800393 

#test within subtypes

ggplot(data=div2[div2$Subtype=="1A",], aes(x=Viral_load, y=Pi, color=Subtype))+
    geom_point()+
    scale_color_manual(values=colors2[c(1,3,5)],guide='none')+
    scale_x_continuous(trans='log10',labels=label_scientific)+
    xlab("viral load")+ylab('nucleotide diversity')+
    theme_bw()+
    ggtitle("1a")

ggsave("Output_all/Diversity/Viral_load.vs.Pi_1A.pdf", width = 4.5, height = 4)

cor.test(div2$Viral_load[div2$Subtype=="1A"], div2$Pi[div2$Subtype=="1A"], method="spearman")
#S = 342260, p-value = 0.8165
#rho = 0.02071239 

ggplot(data=div2[div2$Subtype=="1B",], aes(x=Viral_load, y=Pi, color=Subtype))+
    geom_point()+
    scale_color_manual(values=colors2[c(3,5)],guide='none')+
    scale_x_continuous(trans='log10',labels=label_scientific)+
    xlab("viral load")+ylab('nucleotide diversity')+
    theme_bw()+ggtitle("1b")
ggsave("Output_all/Diversity/Viral_load.vs.Pi_1B.pdf", width = 4.5, height = 4)

cor.test(div2$Viral_load[div2$Subtype=="1B"], div2$Pi[div2$Subtype=="1B"], method="spearman")
#S = 304, p-value = 0.8517
#rho = -0.06293706 

ggplot(data=div2[div2$Subtype=="3A",], aes(x=Viral_load, y=Pi, color=Subtype))+
    geom_point()+
    scale_color_manual(values=colors2[c(5)], guide='none')+
    scale_x_continuous(trans='log10',labels=label_scientific)+
    xlab("viral load")+ylab('nucleotide diversity')+
    theme_bw()+ ggtitle("3a")
ggsave("Output_all/Diversity/Viral_load.vs.Pi_3A.pdf", width = 4.5, height = 4)

cor.test(div2$Viral_load[div2$Subtype=='3A'], div2$Pi[div2$Subtype=='3A'], method="spearman")
#p-value = 0.7284
#rho =0.0785089


### read depth vs. diversity

ggplot(data=div, aes(x=AveDepth, y=Pi, color=Subtype))+
    geom_point()+
    scale_color_manual(values=colors2[c(1,3,5)])+
    scale_x_continuous(trans='log10',labels=label_scientific)+
    xlab("Mean read depth")+ylab('nucleotide diversity')+
    theme_bw()
ggsave("Output_all/Diversity/ReadDepth.vs.Pi.pdf", width = 5.5, height = 4)


cor.test(div$Pi, div$AveDepth, method="spearman")
#S = 2813200, p-value = 0.7748
#rho = -0.01799299 
cor.test(div$AveDepth[div$Subtype=='1A'], div$Pi[div$Subtype=='1A'], method="spearman")
#S = 1292200, p-value = 0.5258
# rho = -0.04566347 
cor.test(div$AveDepth[div$Subtype=='1B'], div$Pi[div$Subtype=='1B'], method="spearman")
#S = 1020, p-value = 0.1346
#  rho 
#0.3376623
cor.test(div$AveDepth[div$Subtype=='3A'], div$Pi[div$Subtype=='3A'], method="spearman")
#S = 10308, p-value = 0.7929
#        rho 
#-0.04331984 


