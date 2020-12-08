#Calculate within- and between-subtype Fst using D-values 
source("Rscripts/baseRscript.R")
library(grid)
library(purrr)
library(zoo)
library(tidyverse)
library(gplots)
library(colorspace)

subs<-c("1A","1B","3A")

# read all seq data
for (g  in 1:3){
    HCVFiles<-list.files(paste0("Output", subs[g],"/SeqData/"),pattern="SeqData")
    SeqData<-list()
    for (i in 1:length(HCVFiles)){ 
        df<-read.csv(paste0("Output", subs[g],"/SeqData/",HCVFiles[i]),stringsAsFactors=FALSE, row.names = 1)
        SeqData[[i]]<-df
        names(SeqData)[i]<-substr(paste(HCVFiles[i]),start=1,stop=7)
        listname<-paste0("SeqData",subs[g])
        assign(listname, SeqData)
    }
    
}

merged.meta<-read.csv("Output_all/merged.metadata.csv", stringsAsFactors = F, row.names = 1)
mergepo<-data.frame(merged.pos=merged.meta[264:8699,"merged.pos"])

for (f in 1:3){
        SeqData1<-get(paste0("SeqData",subs[f]))
        seqdata1<-SeqData1[[1]]
        
        #remove the sites maj!=ref
        seqdata1<-seqdata1[seqdata1$MajNt==paste0(seqdata1$ref),]
        #remove the sites with total reads <1000 
        low_reads<-which(seqdata1$TotalReads<1000) 
        seqdata1[low_reads,c("a","c","g","t")]<-NA
        
        seqdata1<-seqdata1[!is.na(seqdata1$a),]
        a<-seqdata1[,c("pos", "a")]
        c<-seqdata1[,c("pos", "c")]
        g<-seqdata1[,c("pos", "g")]
        t<-seqdata1[,c("pos", "t")]
        
        
        for (i in 2:length(SeqData1)){ 
                seqdata<-SeqData1[[i]]
                seqdata<-seqdata[seqdata$MajNt==paste0(seqdata$ref),]
                low_reads<-which(seqdata$TotalReads<1000) 
                seqdata[low_reads,c("a","c","g","t")]<-NA
                seqdata<-seqdata[!is.na(seqdata$a),]
                
                a<-merge(a, seqdata[,c("pos","a")], by="pos", all = T)
                c<-merge(c, seqdata[,c("pos","c")], by="pos", all = T)
                g<-merge(g, seqdata[,c("pos","g")], by="pos", all = T)
                t<-merge(t, seqdata[,c("pos","t")], by="pos", all = T)
        }
        
        s<-length(SeqData1)
        a$aSum<-rowSums(a[2:(s+1)], na.rm=T)
        c$cSum<-rowSums(c[2:(s+1)], na.rm=T)
        g$gSum<-rowSums(g[2:(s+1)], na.rm=T)
        t$tSum<-rowSums(t[2:(s+1)], na.rm=T)
        
        dat<-merge(a[,c("pos","aSum")],c[,c("pos","cSum")], by="pos", all=T)
        dat<-merge(dat,g[,c("pos","gSum")], by="pos", all=T)
        dat<-merge(dat,t[,c("pos","tSum")], by="pos", all=T)
        colnames(dat)<-c("pos","a","c","g","t")
        dat$Dt<-""
        
        for (k in 1:nrow(dat)){
                if (rowSums(dat[k,2:5])<=3000)  dat$Dt[k]<-NA
                else{
                    ac<-dat[k, "a"]*dat[k,"c"]
                    ag<-dat[k, "a"]*dat[k,"g"]
                    at<-dat[k, "a"]*dat[k,"t"]
                    cg<-dat[k, "c"]*dat[k,"g"]
                    ct<-dat[k, "c"]*dat[k,"t"]
                    gt<-dat[k, "g"]*dat[k,"t"]
                    
                    m<-(rowSums(dat[k,2:5])^2-rowSums(dat[k,2:5]))/2
                    dat$Dt[k]<-(ac+ag+at+cg+ct+gt)/m
                }
        }
        colnames(dat)[1]<-paste0("org.pos.",subs[f])
        dt<-merge(dat, merged.meta, by=paste0("org.pos.",subs[f]), all.y=T )
        dt<-dt[order(dt$merged.pos),]
        write.csv(dt, paste0("Output_all/Diversity/Dt_",subs[f],".csv"))
}        

f=3
dt<-read.csv(paste0("Output_all/Diversity/Dt_",subs[f],".csv"))
levels(dt$gene)
dt$gene<-as.character((dt$gene))
dt$gene[dt$merged.pos>=2598 &dt$merged.pos<=2786]<-"P7"
write.csv(dt, paste0("Output_all/Diversity/Dt_",subs[f],".csv"))

#data trimming (remove extra rows -run only once)
for (i in 1:3){
        dt<-read.csv(paste0("Output_all/Diversity/Dt_", subs[i],".csv"), stringsAsFactors = F, row.names = 1)
        dt<-dt[,c(7,2:6,1,8:31)]

        if (i==1) dt<-dt[dt$merged.pos>263&dt$merged.pos<8679,]
        if (i==2) dt<-dt[dt$merged.pos>263&dt$merged.pos<8639,]
        if (i==3) dt<-dt[dt$merged.pos>263&dt$merged.pos<8642,]
        
        write.csv(dt,paste0("Output_all/Diversity/Dt2_",subs[i],".csv") )
}

       
###### Calculate Fst between subtypes
subs<-c("1A","1B","3A")
CombG<-t(combn(subs,2))
combnames<-c("1A-1B","1A-3A","1B-3A")
FstList<-list()

#normalize the uneven reads between subtypes
for (f in 1:3){
        g1<-CombG[f,1]
        g2<-CombG[f,2]
        fname<-paste0(g1,"-",g2)
        print(fname)
        
        df1<-read.csv(paste0("Output_all/Diversity/Dt2_", g1,".csv"), stringsAsFactors = F, row.names = 1)
        df2<-read.csv(paste0("Output_all/Diversity/Dt2_", g2,".csv"), stringsAsFactors = F, row.names = 1)
        n<-min(nrow(df1),nrow(df2))
        df1<-df1[1:n,]
        df2<-df2[1:n,]
        
        for (i in 2:5){
                df1[,i]<-as.numeric(df1[,i])
                df2[,i]<-as.numeric(df2[,i])
        }
        
        FstDF<-data.frame(merged.pos=df1$merged.pos)
        FstDF$DtT<-""
        
        #divider (m) =(m^2-m)/2   For m=2*m1, divider(m)= 2m1^2-m1
        for (k in 1:nrow(FstDF)){
                if (is.na(df1$Dt[k])|is.na(df2$Dt[k])) FstDF$DtT[k]<-NA
                else{
                        #normalize the reads between the subtypes
                        adj<-sum(df1[k,c("a","g","c","t")])/sum(df2[k,c("a","g","c","t")])

                        ac<-(df1[k, "a"]+df2[k, "a"]*adj)*(df1[k,"c"]+df2[k,"c"]*adj)
                        ag<-(df1[k, "a"]+df2[k, "a"]*adj)*(df1[k,"g"]+df2[k,"g"]*adj)
                        at<-(df1[k, "a"]+df2[k, "a"]*adj)*(df1[k,"t"]+df2[k,"t"]*adj)
                        cg<-(df1[k, "c"]+df2[k, "c"]*adj)*(df1[k,"g"]+df2[k,"g"]*adj)
                        ct<-(df1[k, "c"]+df2[k, "c"]*adj)*(df1[k,"t"]+df2[k,"t"]*adj)
                        gt<-(df1[k, "g"]+df2[k, "g"]*adj)*(df1[k,"t"]+df2[k,"t"]*adj)
                        
                        N<-sum(df1[k,2:5])+adj*sum(df2[k,2:5])
                        m<-(N^2-N)/2
                        FstDF$DtT[k]<-(ac+ag+at+cg+ct+gt)/m
                }
        }
        
        
        FstDF$DtT<-as.numeric(FstDF$DtT)
        
        for (k in 1:nrow(FstDF)){
                FstDF$Dt.ave[k]<-mean(c(df1$Dt[k],df2$Dt[k]))
                FstDF$Fst[k]<-(FstDF$DtT[k]- FstDF$Dt.ave[k])/ FstDF$DtT[k]
        }
        
        tname<-paste0("Fst.",fname)
        assign(tname,FstDF)
        
        FstList[[f]]<-FstDF
        names(FstList)[f]<-tname
        
        write.csv(FstDF,paste0("Output_all/Diversity/",tname,".csv"))
}

### Calculate average Pi for for each subtype and pairs of subtypes
Pi<-data.frame(Sample=c("1A","1B", "3A", combnames ))
for (i in 1:3){
    g1<-CombG[i,1]
    g2<-CombG[i,2]
    fname<-combnames[i]
    #D-values for pairs of subtypes
    dT<-read.csv(paste0("Output_all/Diversity/Fst.", fname,".csv"), stringsAsFactors = F, row.names = 1)
    Pi$pi[Pi$Sample==fname]<-mean(dT$DtT, na.rm=T)
    
    #D-values for a single subtype
    df1<-read.csv(paste0("Output_all/Diversity/Dt2_", subs[i],".csv"), stringsAsFactors = F, row.names = 1)
    Pi$pi[Pi$Sample==subs[i]]<-mean(df1$Dt, na.rm=T)
}

#write.csv("Output_all/Diversity/Pi.csv")


##### 
### Plot Fst & Pi
FstList<-list()
fstfiles<-list.files("Output_all/Diversity/", pattern=glob2rx("*^Fst*.csv*"))
for (i in 1:3){
        FstList[[i]]<-read.csv(paste0("Output_all/Diversity/",fstfiles[i]),stringsAsFactors = F, row.names = 1)
        names(FstList)[i]<-substr(fstfiles[i], start=1, stop = 9)
}

fst.subs<-data.frame(pop=combnames)
for (i in 1:3){
    dat<-FstList[[i]]
    fst.subs$Fst[i]<-mean(dat$Fst, na.rm = T)
}


#colors2<-qualitative_hcl(6, palette="Dark3")
#col2_2<-paste0(colors2,"99")

#Reformat Pi for creating a heatmap 

Pi$Sample1<-substr(Pi$Sample, 1,2)
Pi$Sample2[1:3]<-Pi$Sample1[1:3]
Pi$Sample2[4:6]<-substr(Pi$Sample[4:6], 4,5)
Pi<-Pi[,-1]
Pi<-Pi[,c("Sample1","Sample2","pi")]
Pi[7,1:2]<-c("1B","1A")
Pi[8,1:2]<-c("3A","1A")
Pi[9,1:2]<-c("3A","1B")
Pi$pi[c(7:9)]<-''
Pi$Fst<-''
Pi$Fst[c(1:6)]<-''
Pi$Fst[c(7:9)]<-fst.subs$Fst[1:3]
colnames(Pi)[3]<-"Pi"
Pi$Pi<-as.numeric(Pi$Pi)
Pi$Fst<-as.numeric(Pi$Fst)

#write.csv(Pi,"Output_all/Diversity/Pi.Fst_for_heatmap.csv")
#Pi<-read.csv("Output_all/Diversity/Pi.Fst_for_heatmap.csv", row.names = 1)

colorRampBlue <- colorRampPalette(c("white", "steelblue1", "blue3"))
colorRampRed <- colorRampPalette(c("white", "maroon3"))

ggplot(data=Pi,aes(Sample1,Sample2,fill=Fst))+
    geom_tile(color="gray60", size=.3)+
    scale_fill_gradientn(colors=colorRampBlue(64),limit=c(0,0.4), na.value = 'white')+
    theme_bw()+
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(), 
          panel.grid.major = element_blank(),panel.border = element_blank())+
    #geom_rect(aes(xmin = 1 - 0.5, xmax = 3 + 0.5, ymin = 1 - 0.5, ymax = 3 + 0.5),
    #        fill = "transparent", color = "gray40", size = .3)+
    geom_text(aes(label = round(Fst, 4)))

ggsave("Output_all/Diversity/Fst_lowerTriangle.pdf", height = 5, width = 5.5)

ggplot()+
    geom_tile(data=Pi, aes(Sample1,Sample2, fill=Pi),color="gray60", size=.3)+ 
    #scale_fill_gradientn(colors=colorRampBlue(64))+
    theme_bw()+
    #geom_tile(data=subset(Div, as.integer(factor(Sample1))> as.integer(factor(Sample2))),aes(Sample1,Sample2,fill=Pi))+
    scale_fill_gradientn(colors=colorRampRed(64),limit=c(0,0.2),na.value = "white")+
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(), 
          panel.grid.major = element_blank(),panel.border = element_blank())+
    #geom_rect(aes(xmin = 1 - 0.5, xmax = 3 + 0.5, ymin = 1 - 0.5, ymax = 3 + 0.5),
    #        fill = "transparent", color = "gray40", size = .3)+
    geom_text(data=Pi,aes(Sample1,Sample2, fill=Pi,label = round(Pi, digits=3)))

ggsave("Output_all/Diversity/Pi_upperTriangle.pdf", height = 5, width = 5.5)







