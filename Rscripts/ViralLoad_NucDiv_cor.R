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
ReadDepth$AveDepth<-as.numeric(ReadDepth$AveDepth)

#Attach viral_load data
vl<-read.csv("Data/Viral_load.csv", stringsAsFactors = F)


SampleSummary<-merge(ReadDepth, vl[,c(1,3)], by="Sample_ID")
write.csv(SampleSummary,"Output_all/ReadDepth_Summary.csv")     

#SampleSummary<-read.csv("Output_all/ReadDepth_Summary.csv",stringsAsFactors = F, row.names = 1)   

depths<-data.frame(subtype=c('all',subtypes))
viralloads<-data.frame(subtype=c('all',subtypes))

for (i in 1:4){
    if (i==1){
        depths$mean[i]<-mean(SampleSummary$AveDepth, na.rm = T)
        depths$se[i]<-std.error(SampleSummary$AveDepth, na.rm = T) 
        viralloads$mean[i]<-mean(SampleSummary$Viral_load, na.rm = T)
        viralloads$se[i]<-std.error(SampleSummary$Viral_load, na.rm = T) 
    }
    else {
        k=i-1
        depths$mean[i]<-mean(SampleSummary$AveDepth[SampleSummary$Subtype==subtypes[k]], na.rm = T)
        depths$se[i]<-std.error(SampleSummary$AveDepth[SampleSummary$Subtype==subtypes[k]], na.rm = T) 
        viralloads$mean[i]<-mean(SampleSummary$Viral_load[SampleSummary$Subtype==subtypes[k]], na.rm = T)
        viralloads$se[i]<-std.error(SampleSummary$Viral_load[SampleSummary$Subtype==subtypes[k]], na.rm = T) 
        
    }
}

depths
#  subtype     mean       se
#1     all 6115.048 340.7325
#2      1A 6178.412 426.7979
#3      1B 5890.685 510.3417
#4      3A 5919.040 590.0699

viralloads
#  subtype     mean      se
#1     all 16846508 6176537
#2      1A 19923297 7796421
#3      1B  5460833 1592055
#4      3A  5155555 1317333


wilcox.test(SampleSummary$AveDepth[SampleSummary$Subtype==subtypes[1]], SampleSummary$AveDepth[SampleSummary$Subtype==subtypes[2]],alternative="greater" )
wilcox.test(SampleSummary$AveDepth[SampleSummary$Subtype==subtypes[1]], SampleSummary$AveDepth[SampleSummary$Subtype==subtypes[3]],alternative="greater" )
wilcox.test(SampleSummary$AveDepth[SampleSummary$Subtype==subtypes[3]], SampleSummary$AveDepth[SampleSummary$Subtype==subtypes[2]],alternative="greater" )
#W = 1802, p-value = 0.817
#W = 3710, p-value = 0.5952
#W = 374, p-value = 0.7101

#VL betw/genotypes 
wilcox.test(vl$Viral_load[vl$Subtype=="1A"], vl$Viral_load[vl$Subtype=="3A"], "greater", paired =F)
#W = 1354.5, p-value = 0.6129
wilcox.test(vl$Viral_load[vl$Subtype=="1A"], vl$Viral_load[vl$Subtype=="1B"], "greater", paired =F)
#W = 682.5, p-value = 0.739
wilcox.test(vl$Viral_load[vl$Subtype=="1B"], vl$Viral_load[vl$Subtype=="3A"], "greater", paired =F)
#W = 141, p-value = 0.3796


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
    xlab("Viral load")+ylab(expression(paste("Nucleotide diversity (",pi,")")))+
    theme_bw()
ggsave("Output_all/Diversity/Viral_load.vs.Pi.pdf", width = 5.5, height = 4)


#Test viral load vs. pi
cor.test(div2$Viral_load, div2$Pi, method="spearman")
#S = 704000, p-value = 0.9351
#rho =0.006443285  

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
#S = 343010, p-value = 0.8351
#rho = 0.01857932  

ggplot(data=div2[div2$Subtype=="1B",], aes(x=Viral_load, y=Pi, color=Subtype))+
    geom_point()+
    scale_color_manual(values=colors2[c(3,5)],guide='none')+
    scale_x_continuous(trans='log10',labels=label_scientific)+
    xlab("viral load")+ylab('nucleotide diversity')+
    theme_bw()+ggtitle("1b")
ggsave("Output_all/Diversity/Viral_load.vs.Pi_1B.pdf", width = 4.5, height = 4)

cor.test(div2$Viral_load[div2$Subtype=="1B"], div2$Pi[div2$Subtype=="1B"], method="spearman")
#S = 303.03, p-value = 0.8542
#rho = -0.05954475 

ggplot(data=div2[div2$Subtype=="3A",], aes(x=Viral_load, y=Pi, color=Subtype))+
    geom_point()+
    scale_color_manual(values=colors2[c(5)], guide='none')+
    scale_x_continuous(trans='log10',labels=label_scientific)+
    xlab("viral load")+ylab('nucleotide diversity')+
    theme_bw()+ ggtitle("3a")
ggsave("Output_all/Diversity/Viral_load.vs.Pi_3A.pdf", width = 4.5, height = 4)

cor.test(div2$Viral_load[div2$Subtype=='3A'], div2$Pi[div2$Subtype=='3A'], method="spearman")
#S = 1626, p-value = 0.7171
#rho =0.08189777 


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
cor.test(div3$AveDepth[div3$Subtype=='1A'], div3$Pi[div3$Subtype=='1A'], method="spearman")
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


### read depth vs. VL
#remove the 2 outliers in read depth for correlation analysis
div2<-div2[div2$AveDepth<=40000,]

cor.test(div2$AveDepth,div2$Viral_load,  method="spearman")
#S = 757420, p-value = 0.1679
#rho = -0.109548 
cor.test(div2$AveDepth[div2$Subtype=='1A'], div2$Viral_load[div2$Subtype=='1A'], method="spearman")
#S = 359050, p-value = 0.3913
# rho = -0.07702219 
cor.test(div2$AveDepth[div2$Subtype=='1B'], div2$Viral_load[div2$Subtype=='1B'], method="spearman")
#S = 275.98, p-value = 0.9139
#  rho 
#0.03502632 
cor.test(div2$AveDepth[div2$Subtype=='3A'], div2$Viral_load[div2$Subtype=='3A'], method="spearman")
#S = 2390.2, p-value = 0.1107
#        rho 
#-0.3496188 

#plot the results
ggplot(data=div2, aes(x=Viral_load, y=AveDepth, color=Subtype))+
    geom_point()+
    scale_color_manual(values=colors2[c(1,3,5)])+
    scale_x_continuous(trans='log10',labels=label_scientific)+
    xlab("Viral load")+ylab('Mean read depth')+
    theme_bw()
ggsave("Output_all/Diversity/ReadDepth.vs.VL.pdf", width = 5.5, height = 4)





#Plot the viral load of each subtype
ggplot(data=div2, aes(x=Subtype, y=Viral_load))+
    geom_boxplot()+
    #scale_color_manual(values=colors2[c(5)], guide='none')+
    scale_y_continuous(trans='log10',labels=label_scientific)+
    xlab("Subtype")+ylab('Viral load')+
    theme_bw()

