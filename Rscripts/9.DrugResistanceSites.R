# Look at RAVs (frequencies and proportions) and create a plot
library(tidyverse)
library(colorspace)
source("Rscripts/baseRscript.R")

cold<-qualitative_hcl(6, palette="Cold")
colors<-cold[c(1:3)]
#Dark shades of the above colors
colorsD<-c("#8E86C4", "#6990C1","#4198B6")

#specific nucleotide positions with known drug resistant mutations 
drsites<-read.csv("Data/HCV_RAVsites.csv",stringsAsFactors = F)
meta<-read.csv("Output_all/merged.metadata.csv",stringsAsFactors = F, row.names = 1)

#attach the original position info for all genotypes
for (i in 1:nrow(drsites)){
        drsites$pos.1A[i]<-meta$org.pos.1A[which(meta$merged.pos==drsites$merged.pos[i])]
        drsites$pos.1B[i]<-meta$org.pos.1B[which(meta$merged.pos==drsites$merged.pos[i])]
        drsites$pos.3A[i]<-meta$org.pos.3A[which(meta$merged.pos==drsites$merged.pos[i])]
}


#########

#1. extract the mut freq at drug resistant sites.
subs<-c("1A","1B","3A")
for (j in 1:length(subs)){
        HCVFiles<-list.files(paste0("Output",subs[j],"/Overview2/"), pattern="overview2.csv")
        
        dr.sites<-list()
        diff<-list()
        diff.count<-list()
        DR_mutfreq<-list()
        
        #first deal with the DRV with one mutation 
        DR<-drsites[drsites$genotype==subs[j],]
        
        #create an id column
        for (i in 1:nrow(DR)){
                if (DR$Need_both[i]=="y") DR$ID[i]<- paste0(DR$Name[i],'.',DR$merged.pos[i])
                else DR$ID[i]<-DR$Name[i]
        }
        
        for (i in 1:length(HCVFiles)){ 
                df<-read.csv(paste0("Output",subs[j],"/Overview2/",HCVFiles[i]),stringsAsFactors=FALSE, row.names = 1)
                dname<-substr(paste(HCVFiles[i]),start=1,stop=7)
                dr<-DR
                cname<-paste0("pos.",subs[j])
                DRsites<-df[df$pos %in% dr[,cname],]

                #count the number of samples fixed with the RAVs
                dr$obs<-0
                for (k in 1:nrow(dr)){
                        pos<-dr[k,cname]
                        if (is.na(DRsites$MajNt[DRsites$pos==pos])|is.na(DRsites$ref[DRsites$pos==pos])) next
                        if (DRsites$MajNt[DRsites$pos==pos]!=DRsites$ref[DRsites$pos==pos]){
                                if (dr$Type[k]=="Ts") {mutnt<-transition(DRsites$ref[DRsites$pos==pos]) 
                                        if (DRsites$MajNt[DRsites$pos==pos]==mutnt) dr$obs[k]<-1} 
                                if (dr$Type[k]=='Tv1') {mutnt<-transv1(DRsites$ref[DRsites$pos==pos])
                                        if (DRsites$MajNt[DRsites$pos==pos]==mutnt) dr$obs[k]<-1}
                                if (dr$Type[k]=='Tv2') {mutnt<-transv2(DRsites$ref[DRsites$pos==pos])
                                        if (DRsites$MajNt[DRsites$pos==pos]==mutnt) dr$obs[k]<-1} 
                                if (dr$Type[k]=='Tv') {mutnt1<-transv1(DRsites$ref[DRsites$pos==pos])
                                        mutnt2<-transv2(DRsites$ref[DRsites$pos==pos])
                                        if (DRsites$MajNt[DRsites$pos==pos]==mutnt1|DRsites$MajNt[DRsites$pos==pos]==mutnt2) dr$obs[k]<-1  }    
                                if (dr$Type[k]=='Ts,Tv1') {mutnt1<-transition(DRsites$ref[DRsites$pos==pos])
                                        mutnt2<-transv1(DRsites$ref[DRsites$pos==pos])    
                                        if (DRsites$MajNt[DRsites$pos==pos]==mutnt1|DRsites$MajNt[DRsites$pos==pos]==mutnt2) dr$obs[k]<-1 }  
                                if (dr$Type[k]=='Ts,Tv2') {mutnt1<-transition(DRsites$ref[DRsites$pos==pos])
                                        mutnt2<-transv2(DRsites$ref[DRsites$pos==pos])    
                                        if (DRsites$MajNt[DRsites$pos==pos]==mutnt1|DRsites$MajNt[DRsites$pos==pos]==mutnt2) dr$obs[k]<-1 }     
                                if (dr$Type[k]=='Ts,Tv') dr$obs[k]<-1
                        }
                }
                
                #store the observation number in a list
                diff.count[i]<-sum(dr$obs==1)
                names(diff.count)[i]<-dname
                
                #obtain the mutation freq.
                dr$freq<-""
                for (k in 1:nrow(dr)){
                        if (dr$Type[k]=='Tv1')    dr$freq[k]<-DRsites$freq.transv1.ref[DRsites$pos==dr[k,cname]]
                        if (dr$Type[k]=='Tv2')    dr$freq[k]<-DRsites$freq.transv2.ref[DRsites$pos==dr[k,cname]]
                        if (dr$Type[k]=='Tv')     dr$freq[k]<-DRsites$freq.transv.ref [DRsites$pos==dr[k,cname]]
                        if (dr$Type[k]=='Ts')     dr$freq[k]<-DRsites$freq.Ts.ref     [DRsites$pos==dr[k,cname]]
                        if (dr$Type[k]=='Ts,Tv1') dr$freq[k]<-(DRsites$freq.Ts.ref    [DRsites$pos==dr[k,cname]]+DRsites$freq.transv1.ref[DRsites$pos==dr[k,cname]])
                        if (dr$Type[k]=='Ts,Tv2') dr$freq[k]<-(DRsites$freq.Ts.ref    [DRsites$pos==dr[k,cname]]+DRsites$freq.transv2.ref[DRsites$pos==dr[k,cname]])
                        if (dr$Type[k]=='Ts,Tv')  dr$freq[k]<-DRsites$freq.mutations.ref[DRsites$pos==dr[k,cname]]
                }
                
                if (i==1) write.csv(dr, paste0("Output_all/DR/DRsites.",dname,".csv"))
                dr$freq<-as.numeric(dr$freq)
                # fixed to RAV 
                diff[[i]]<-dr[,c("ID","obs")]
                names(diff)[i]<-dname
                #mut.freq
                DR_mutfreq[[i]]<-dr[,c("ID","freq")]
                names(DR_mutfreq)[i]<-dname
        }
        
        for (i in 1:length(DR_mutfreq)) {
                colnames(DR_mutfreq[[i]])<-c("ID",paste0(names(DR_mutfreq[i])))
                colnames(diff[[i]])<-c("ID",paste0(names(diff[i])))
        }

        DR.mutated.counts<-do.call(rbind,diff.count)
        
        DR_MutFreq<-DR_mutfreq%>% purrr::reduce(full_join, by='ID')
        DR_diff<-diff %>% purrr::reduce(full_join, by='ID')
        write.csv(DR_MutFreq, paste0("Output_all/DR/RAV.MutationFreq_summary.",   subs[j],".csv"))
        write.csv(DR_diff,    paste0("Output_all/DR/RAV.counts.MutFreq_summary.", subs[j],".csv"))
}
 
#create a plot
subs<-c("1A","1B","3A")
fixed<-data.frame(Subtype=subs)
for (j in 1:length(subs)){
        HCVFiles<-list.files(paste0("Output",subs[j],"/Overview2/"), pattern="overview2.csv")
        s<-length(HCVFiles)
        
        DR_MutFreq<-read.csv(paste0("Output_all/DR/RAV.MutationFreq_summary.",  subs[j],".csv"), stringsAsFactors = F, row.name=1)
        DR_diff   <-read.csv(paste0("Output_all/DR/RAV.counts.MutFreq_summary.",subs[j],".csv"), stringsAsFactors = F, row.name=1)

        #ns3n <-nrow(drsites[drsites$Gene=="NS3" &drsites$genotype==geno[j],])
        #ns5an<-nrow(drsites[drsites$Gene=="NS5A"&drsites$genotype==geno[j],])
        #ns5bn<-nrow(drsites[drsites$Gene=="NS5B"&drsites$genotype==geno[j],])
        
        if (j==1) dr<-read.csv("Output_all/DR/DRsites.D75002-.csv", stringsAsFactors = F, row.names = 1)
        if (j==2) dr<-read.csv("Output_all/DR/DRsites.D75046-.csv", stringsAsFactors = F, row.names = 1)
        if (j==3) dr<-read.csv("Output_all/DR/DRsites.D75003-.csv", stringsAsFactors = F, row.names = 1)
        
        ns3n <-nrow(dr[dr$Gene=="NS3" ,])
        ns5an<-nrow(dr[dr$Gene=="NS5A",])
        ns5bn<-nrow(dr[dr$Gene=="NS5B",])
        
        #count the number of patients with fixed RAVs
        DR_diff$total<-apply(DR_diff[2:(s+1)],1,sum, na.rm=T)
        #% of patients with fixed RAVs
        DR_diff$Percent_all<-format(round(DR_diff$total/s*100, 1), nsmall=1)
        DR_diff$Percent_all<-as.numeric(DR_diff$Percent_all)
        DR_diff$percent<-as.numeric(DR_diff$Percent_all)
        #% of patients with fixed RAVs for the sites that had at least 1 fixed sample
        DR_diff$percent[DR_diff$percent==0]<-NA
        
        #samples that had fixed RAVs
        fixed$no_fixed_RAVs[j]<-nrow(DR_diff[!is.na(DR_diff$percent),])
        #proportion of RAVsites that had samples with fixed RAVs
        fixed$prop_fixed_RAVs[j]<-fixed$no_fixed_RAVs[j]/nrow(DR_MutFreq)
        #proportion of samples with fixed RAVs
        fixed$prop_fixed_samples[j]<-fixed$no_fixed_RAVs[j]/s
        #average % fixed RAVs
        fixed$mean.percent_all[j]<-mean(DR_diff$Percent_all, na.rm=T)
        #Average % of patients that had fixed RAVs for those sites that had at least 1 sample with fixed RAV
        fixed$mean.percent_fixedOnly[j]<-mean(DR_diff$percent, na.rm=T)
        
        #create a figure
        pdf(paste0("Output_all/DR/DrugResi.",subs[j],".pdf"), height = 8, width = 15)
        layout(matrix(c(1,2), byrow = TRUE), heights=c(1,3))
        par(mar=c(0, 4, 4, .5))
        
        plot(DR_diff$percent, type = "n", xlim = c(1, nrow(dr)), ylim=c(0,55),axes = FALSE, ylab = "% fixed to RAVs", xlab = "")
        axis(side=2, col="gray30")
        abline(v=c((0:nrow(dr))+0.5), col="gray60", lty=1, lwd=.1)
        abline(v=ns3n+.5,col='gray50', lwd=3)
        abline(v=ns3n+ns5an+.5,col='gray50', lwd=3 )
        
        #Plots of proportion of fixed RAVs
        for (i in 1:nrow(dr)){
                if (dr$Type[i]=='Ts') segments(i,0,i,DR_diff$percent[i],col="#87AEDF",lwd=10,lend=2)
                #if (dr$Type[i]=='Tv1'|dr$Type[i]=='Tv2'|dr$Type[i]=='Tv') barplot(DR_diff$percent[i],col="purple",add=T)
                if (dr$Type[i]=='Tv1'|dr$Type[i]=='Tv2'|dr$Type[i]=='Tv') segments(i,0,i,DR_diff$percent[i],,lwd=10, col="purple",lend=2)
                }
        
        rect(-1,42,5,50, col="white",border ="white")
        rect(-1,48 ,0, 50,density = NULL, col="#87AEDF",border ="#87AEDF")
        text(1,48,"Transition",col="black", cex=.8, adj=0)
        rect(-1,42 ,0, 44,density = NULL, col="purple",border ="purple")
        text(1,42,"Transversion",col="black", cex=.8, adj=0)

        
        # mf plots
        par(mar=c(6, 4, 0.5, .5))
        ymin <- -5
        plot(0, type = "n", xlim = c(1, nrow(DR_MutFreq)), ylim = c(ymin, 0), axes = FALSE, ylab = "Frequency of RAVs", 
             xlab = "")
        axis(side = 2, at = seq(0, ymin, by=-1), labels = expression(10^0, 10^-1, 10^-2, 10^-3, 10^-4,10^-5), las = 2)
        n<-seq(1, by = 2, len = (nrow(DR_MutFreq)/2))

        genesDF<-data.frame("name"= c("NS3","NS5A","NS5B"), "Begin"= c(1,ns3n+1,(ns3n+ns5an+1)),"End"= c(ns3n,(ns3n+ns5an),nrow(DR_MutFreq)))
        for (i in 1:3){
                nvec<-seq(genesDF$Begin[i],genesDF$End[i],2)
                nvec2<-seq(genesDF$Begin[i]+1,genesDF$End[i],2)
                
                if (nvec[1]%%2==0) {v2<-nvec; v1<-nvec2}
                if (nvec[1]%%2==1) {v1<-nvec; v2<-nvec2}
                        
                abline(v=v1, col=paste0(colors[i],"66"),lty=1, lwd=16)
                abline(v=v2, col=paste0(colors[i],"1A"),lty=1, lwd=16)
        }
        
        for (i in 1:nrow(DR_MutFreq)){
                if (i<=genesDF$End[1]){
                        xjit <- rnorm(s, 0, .1)
                        points(rep(i,s)+xjit, log10(DR_MutFreq[i,2:(s+1)]),pch=16,cex=0.4,col=colorsD[1])}
                if (i<=genesDF$End[2]& i>=genesDF$Begin[2]){
                        xjit <- rnorm(s, 0, .1)
                        points(rep(i,s)+xjit, log10(DR_MutFreq[i,2:(s+1)]),pch=16,cex=0.4,col=colorsD[2])}
                if (i<=genesDF$End[3]& i>=genesDF$Begin[3]){
                        xjit <- rnorm(s, 0, .1)
                        points(rep(i,s)+xjit, log10(DR_MutFreq[i,2:(s+1)]),pch=16,cex=0.3, col=colorsD[3])}
        }
        
        #mtext(side = 3, at = 1:nrow(DR_MutFreq), text = paste(DR_diff$total), cex = .8)
        mtext(side = 1, at = 1:nrow(DR_MutFreq), text = paste(dr$ID), las=2,padj=0, cex = .8)
        
        abline(v=ns3n+.5,col='gray50', lwd=3)
        abline(v=ns3n+ns5an+.5,col='gray50', lwd=3 )
        
        rect(0.5,-5,ns3n+.5,-5.19,density = NULL, angle = 45,col="white",border =colors[1])
        text(ns3n/2+.5,-5.08,paste0(genesDF$name[1]),col="black", cex=.8)
        rect(ns3n+.5,-5, ns3n+ns5an+.5 ,-5.19,density = NULL, angle = 45,col="white",border =colors[2])
        text(ns3n+ns5an-ns5an/2+.5,-5.08,paste0(genesDF$name[2]),col="black", cex=.8)
        rect(ns3n+ns5an+.5,-5,nrow(dr)+.5,-5.19,density = NULL, angle = 45,col="white",border =colors[3])
        text(nrow(dr)-(nrow(dr)-ns3n-ns5an)/2+.5,-5.08,paste0(genesDF$name[3]),col="black", cex=.8)
        dev.off()
}

write.csv(fixed,"Output_all/DR/Fixed_sites_summary.csv")


subs<-c("1A","1B","3A")
ravF<-data.frame(Subtype=subs)
ravF_gene<-list()
Fixedrav<-list()
for (j in 1:length(subs)){
        HCVFiles<-list.files(paste0("Output",subs[j],"/Overview2/"), pattern="overview2.csv")
        s<-length(HCVFiles)
        
        DR_MutFreq<-read.csv(paste0("Output_all/DR/RAV.MutationFreq_summary.", subs[j],".csv"),  stringsAsFactors = F, row.name=1)
        df<-DR_MutFreq
        df[df>0.5]<-NA
        
        DR_MutFreq$Mean<-rowMeans(DR_MutFreq[,2:ncol(DR_MutFreq)], na.rm=T)
        df$Mean<-rowMeans(df[,2:ncol(df)], na.rm = T)
        
        ravF$MeanFreq[j]<-mean(DR_MutFreq$Mean, na.rm=T)
        ravF$MeanFreq_notFixed[j]<-mean(df$Mean)
        
        #a number of patients that showed rav
        DR_MutFreq$No.samples.withRAVs<-rowSums(DR_MutFreq[,2:(s+1)]!=0, na.rm=T)
        #ravF$MeanCount[j]<-mean(DR_MutFreq$Counts, na.rm=T)
        ravF$Percent.samples.withRAVs[j]<-mean(DR_MutFreq$No.samples.withRAVs, na.rm=T)/s*100
        
        ravsites<- drsites[drsites$genotype==subs[j],]
        
        #look at mut freq for differnet genes and types
        DR_MutFreq$Gene<-ravsites$Gene
        DR_MutFreq$multi.change<-ravsites$Need_both
        df$Gene<-ravsites$Gene
        df$multi.change<-ravsites$Need_both
        rav_gene<-aggregate(DR_MutFreq$No.samples.withRAVs,by=list(DR_MutFreq$Gene), mean) 
        colnames(rav_gene)<-c("Gene","Counts")
        rav_gene$Proportion<-rav_gene$Counts/s
        
        #non-fixed RAV freq.
        mf<-aggregate(df$Mean,by=list(df$Gene), mean, na.rm=T)
        rav_gene$meanMF<-mf$x
        #compare mut freq between rav that need two changes 
        mf2<-aggregate(df$Mean, by=list(df$multi.change, df$Gene), mean)
        rav_gene$meanMF.y<-mf2$x[mf2$Group.1=="y"]
        rav_gene$meanMF.n<-mf2$x[mf2$Group.1=="n"]
        
        
        
        ravF_gene[[j]]<-rav_gene
        names(ravF_gene)[j]<-subs[j]
        
        
        DR_diff<-read.csv(paste0("Output_all/DR/RAV.counts.MutFreq_summary.",subs[j],".csv"), stringsAsFactors = F, row.name=1)
        
        DR_diff$Gene<-ravsites$Gene
        DR_diff$Need_both<-ravsites$Need_both
        #count the number of patients with fixed RAVs
        DR_diff$total<-apply(DR_diff[2:(s+1)],1,sum, na.rm=T)
        
        #Look at the sites with fixed RAVs 
        drd<-DR_diff[DR_diff$total>0,]
        fix.ct<-aggregate(DR_diff$total, by=list(DR_diff$Need_both, DR_diff$Gene), mean)
        colnames(fix.ct)<-c("Need_both","Gene","Ave.no.of.fixed.samples")
        fix.ct$Ave.percent.fixed.samples.<-fix.ct$Ave.no.of.fixed.samples/s*100
        fix.ct2<-aggregate(drd$total, by=list(drd$Need_both, drd$Gene), mean)
        if (j==3){
            fix.ct$Ave.percent.fixed.samples_withFixed<-fix.ct2$x/s*100
        }
        else{
            fix.ct$Ave.percent.fixed.samples_withFixed<-c(fix.ct2$x/s*100, 0, 0)
        }
        #without q80k for 1a
        if (j==1){
            dr2<-DR_diff[DR_diff$ID!="Q80K",]
            fix1<-aggregate(dr2$total, by=list(dr2$Need_both, dr2$Gene), mean)
            colnames(fix1)<-c("Need_both","Gene","Ave.no.of.fixed.samples")
            fix1$Ave.percent.fixed.samples.<-fix1$Ave.no.of.fixed.samples/s*100
            drd2<-dr2[dr2$total>0,]
            fix2<-aggregate(drd2$total, by=list(drd2$Need_both, drd2$Gene), mean)
            fix1$Ave.percent.fixed.samples_withFixed<-c(fix2$x/s*100, 0, 0)
            
        }
        
        wilcox.test(DR_diff$total[DR_diff$Need_both=="y"],DR_diff$total[DR_diff$Need_both=="n"],alternative = "greater" )
        #W = 746.5, p-value = 0.257 (3a)
        #W = 846, p-value = 0.03871 (1a)
        #W = 839, p-value = 0.2541
        
        
        ct<-data.frame(table(DR_diff$Need_both, DR_diff$Gene))
        fix.ct$number<-ct$Freq
        Fixedrav[[j]]<-fix.ct
        names(Fixedrav)[j]<-subs[j]
        
        
}
write.csv(ravF, "Output_all/DR/RAV_freq.csv")

ravf<-data.frame()
fixr<-data.frame()
for (i in 1:3){
    df<-ravF_gene[[i]]
    df$subtype<-subs[i]
    ravf<-rbind(ravf,df)
    df2<-Fixedrav[[i]]
    df2$subtype<-subs[i]
    fixr<-rbind(fixr,df2)
    
}

write.csv(ravf, "Output_all/DR/RAV.freq.gene.summary.csv")
write.csv(fixr, "Output_all/DR/Fixed_rav.summary.csv")
