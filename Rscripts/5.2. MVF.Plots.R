library(reshape2)
library(colorspace)
library(ggplot2)
source("Rscripts/baseRscript.R")
source("Rscripts/label_scientific.R")

#Subtype colors 
colors2<-qualitative_hcl(6, palette="Dark3")
#1A=1, 1B=3, 3A=5, c("#E16A86","#50A315","#009ADE")

sub<-c("1A","1B","3A")


## All MVF figures ##
Summary<-read.csv("Output_all/MVF/All.MinorVariant_Mean_3subtypes.csv",stringsAsFactors = F, row.names = 1)

## Plot across genomes all together (Figure 1)
Summary2<-Summary[Summary$merged.pos<=8500 & Summary$merged.pos>=264,]

genes<-read.csv("Data/HCV_annotations_joined.csv", stringsAsFactors = F)
genes$Gene[genes$Gene=="NS1(P7)"]<-"NS1"
genenames<-genes$Gene[1:12]

Summary2$roll1001a<-c(rep(NA, times=99), rollmean(Summary2$mean.1A, k=100, na.rm=T))
Summary2$roll1001b<-c(rep(NA, times=99), rollmean(Summary2$mean.1B, k=100, na.rm=T))
Summary2$roll1003a<-c(rep(NA, times=99), rollmean(Summary2$mean.3A, k=100, na.rm=T))


#craete a plot
pdf("Output_all/Figures/All.MVF_AcrossGenome.pdf",width=12,height=4.5)
plot(mean.1A~merged.pos, data=Summary2,t="n",log='y',yaxt='n',xlab='Genome position', ylab="Minor variant frequency",
     ylim=c(10^-4,0.4),xlim=c(264,8500))
eaxis(side = 2, at = 10^((-0):(-(4))), cex=2)
for(k in 1:4){abline(h = 1:10 * 10^(-k), col = "gray90",lty=1,lwd=.5)}

points(mean.1B~merged.pos,pch=16, data=Summary2, col=paste0(colors2[3],"B3"), cex=0.2)
points(mean.3A~merged.pos,pch=16, data=Summary2, col=paste0(colors2[5],"B3"), cex=0.2)
points(mean.1A~merged.pos,pch=16, data=Summary2, col=paste0(colors2[1],"B3"), cex=0.2)

lines(roll1001b~merged.pos, data=Summary2, col=colors2[3],lwd=1)
lines(roll1003a~merged.pos, data=Summary2, col=colors2[5],lwd=1)
lines(roll1001a~merged.pos, data=Summary2, col=colors2[1],lwd=1)

legend("topright",legend=c("1a","1b","3a"), col=colors2[c(1,3,5)], pch=16, bty = "n", cex=1)

genes2<-genes[1:12,]
#genes2$Gene[6]<-"NS1"

abline(v=genes2$end, col="gray80", lwd=.5)

ylow<-0.0001;yhigh<-.4
for (j in 1:nrow(genes2)){
        xleft<- genes2$start[j]
        xright<-genes2$start[j+1]
        if (j==1){
                
                rect(-200,ylow,xright,ylow*1.6,density = NULL, angle = 45,col="white",border ="gray40",lwd=1.5)
                text(150,1.27*ylow,paste0(genes2$Gene[j]),col="black", cex=0.7)
        }
        else if (j==6){
                rect(xleft,ylow,xright,ylow*1.6,density = NULL, angle = 45,col="white",border ="gray40",lwd=1.5)
                text(xleft+100,1.27*ylow,paste0(genes2$Gene[j]),col="black", cex=0.7)
                
        }
        else if (j==4|j==9){
                rect(xleft,ylow,xright,1.6*ylow,density = NULL, angle = 45,col="white",border ="gray40",lwd=1.5)
                text(xleft+50,2.1*ylow,paste0(genes2$Gene[j]),col="black", cex=0.7)
        }
        else if (j==12){
                rect(xleft,ylow,genes2$end[j],1.6*ylow,density = NULL, angle = 45,col="white",border ="gray40",lwd=1.5)
                text(xleft+600,1.27*ylow,paste0(genes2$Gene[j]),col="black", cex=0.7)
        }
        else{rect(xleft,ylow,xright,1.6*ylow,density = NULL, angle = 45,col="white",border ="gray40",lwd=1.5)
                text(xright-(xright-xleft)/2,1.27*ylow,paste0(genes2$Gene[j]),col="black", cex=0.7)}
}
box()
dev.off()


####
#plot the means with boxplots

tb<-read.csv ("Output_all/MVF/Summary_3subtypes.csv", stringsAsFactors = F, row.names = 1)
Sum<-read.csv("Output_all/MVF/MVF.all.Ts.Tvs.depth.csv", row.names = 1,stringsAsFactors = F)
Sum<-Sum[1:8800,]
tb2<-Sum[,c(2:10)]
tb2m<-melt(tb2)
tb2m$Subtype<-c(rep("1A", times=8800),rep("1B", times=8800),rep("3A", times=8800),rep("1A", times=8800),rep("1B", times=8800),rep("3A", times=8800),rep("1A", times=8800),rep("1B", times=8800),rep("3A", times=8800))
tb2m$Type<-c(rep("Transition", times=8800*3), rep("Transversion", times=8800*3), rep("Total MV", times=8800*3))
tb2m<-tb2m[!is.na(tb2m$value),]

SumT<-read.csv("Output_all/MVF/SummaryAll.csv", stringsAsFactors = F, row.names = 1)

col80<- paste0(colors2,"CC")
col80<-col80[c(1,3,5)]

#Wilcoxon rank-sum test
Comb<-t(combn(1:3, 2))
combnames<-c("1a-1b","1a-3a","1b-3a")

results.Ts<-list()
results.Tvs<-list()
results.all<-list()

wilcox.res<-data.frame(matrix(nrow=3, ncol=3))
rownames(wilcox.res)<-c("Ts","Tvs","All")
colnames(wilcox.res)<-combnames                       
for (i in 1:3){
    k=1
    n1<-Comb[i,1]+k
    n2<-Comb[i,2]+k
    r1<-wilcox.test(SumT[,n1], SumT[,n2], alternative = "greater", paired = FALSE) 
    wilcox.res[1,i]<-r1[[3]][1]
    
    k=4
    n1<-Comb[i,1]+k
    n2<-Comb[i,2]+k
    r2<-wilcox.test(SumT[,n1], SumT[,n2], alternative = "greater", paired = FALSE) 
    wilcox.res[2,i]<-r2[[3]][1]
    
    k=7
    n1<-Comb[i,1]+k
    n2<-Comb[i,2]+k
    r3<-wilcox.test(SumT[,n1], SumT[,n2], alternative = "greater", paired = FALSE) 
    wilcox.res[3,i]<-r2[[3]][1]
    
}
write.csv(wilcox.res,"Output_all/MVF/WilcoxTest_MVF.bySubtype.csv",row.names = T)

Test.results<-data.frame(mutation=c("Ts","Tvs","All"))
for (i in 1:3){
    symb1<-NA
    symb2<-NA
    symb3<-NA
    if (wilcox.res$`1a-1b`[i]<0.05) symb1="a"
    if (wilcox.res$`1a-3a`[i]<0.05) symb2="b"
    if (wilcox.res$`1b-3a`[i]<0.05) symb3="c"
    
    symb<-c(symb1, symb2, symb3) 
    symb<-symb[!is.na(symb)]
    symb<-paste(symb, collapse = ',' )
    if (length(symb)==0) symb<-''
    Test.results$Sig[i]<-symb
}



ggplot()+
        geom_boxplot(data=tb2m,aes(x=Type, y=value,fill=Subtype, middle=mean(value), 
                                   color=Subtype),outlier.alpha = 0.2)+ 
        labs(x="")+ylab("Minor variant frequency")+
        scale_y_continuous(trans = 'log10', labels=label_scientific)+
        scale_fill_manual(values=paste0(colors2[c(1,3,5)],"CC"), labels=c("1a","1b","3a"))+
        geom_point(data=tb,aes(x=Type,y=Mean,color=Subtype, fill=Subtype), position=position_dodge(.75), shape=21, color="gray30", bg=rep(colors2[c(1,3,5)], times=3))+
        geom_errorbar(data=tb,aes(x=Type,y=Mean, ymin=Mean-SE, ymax=Mean+SE, fill=Subtype, color=Subtype),position=position_dodge(.75), width=.2, color="gray40")+
        scale_color_manual(values=colors2[c(1,3,5)])+
        theme_bw()+
        theme(axis.text.x = element_text(size=13),axis.title.y = element_text(size=15), axis.text.y=element_text(size=12))+
        geom_vline(xintercept = c(1:2)+0.5,  
                   color = "gray70", size=.3)+
        theme(panel.grid.major.x = element_blank())+
        theme(legend.position = "none")+
        annotate(geom="text", x=1:3, y=0.2,  label=Test.results$Sig[c(3,1,2)], color="red", size=2.5 )

ggsave(filename="Output_all/Figures/AveMVF.bytype.bysubtype_boxplot.pdf",width = 5.5, height = 6)


########
## Plot means and SEs for each gene

#attach the gene info to Sum
genes<-read.csv("Data/HCV_annotations_joined.csv", stringsAsFactors = F)
gene.vector<-c()
for (i in 1:(nrow(genes)-1)) gene.vector<-c(gene.vector, rep(genes$Gene[i],times=genes$start[i+1]-genes$start[i]))
genetable<-data.frame("merged.pos"=c(1:length(gene.vector)))
genetable$Gene<-gene.vector

Sum<-merge(genetable,Sum, by="merged.pos")

mf<-data.frame()
subs<-c("1A","1B","3A")
for (g in 1:3){
        dat<-Sum[,c("merged.pos","Gene",paste0("mean.all.",subs[g]),paste0("depth.",subs[g]))]
        colnames(dat)[3:4]<-c("mean","depth")
        dat$Subtype<-subs[g]
        mf<-rbind(mf,dat)
}

mf1<-mf[!is.na(mf$mean),]
SumMFGenes<-aggregate(mf1$mean,by=list(mf1$Subtype,mf1$Gene),FUN=mean)
meSE<-data.frame(Subtype=rep(subs, times=3), Gene=rep(genenames, each=3))
i=1
for (k in 1:length(genenames)){
        for (g in 1:3){
                df<-mf1[mf1$Gene==genenames[k]&mf1$Subtype==subs[g],]
                meSE$SE[i]<-sqrt(mean(df$mean)*(1-mean(df$mean))/mean(df$depth, na.rm=T))
                i=i+1
        }
}
sumG<-cbind(SumMFGenes, meSE$SE)
colnames(sumG)<-c("Subtype","Gene","Mean","SE")
write.csv(sumG,"Output_all/MVF/Summary.byGene.csv")

#which subtype has the highest mean?
for (i in 1:12){
        df<-sumG[sumG$Gene==genenames[i],]
        m<-which.max(df$Mean)
        print(paste(genenames[i], subs[m]))
}

#which subtype has the lowest mean?
for (i in 1:12){
        df<-sumG[sumG$Gene==genenames[i],]
        m<-which.min(df$Mean)
        print(paste(genenames[i], subs[m]))
}

sumG$Gene<-factor(sumG$Gene, levels=c("5' UTR","Core","E1", "HVR1","E2","NS1", "NS2","NS3","NS4A","NS4B","NS5A","NS5B"))
mf1$Gene<-factor(mf1$Gene, levels=c("5' UTR","Core","E1", "HVR1","E2","NS1", "NS2","NS3","NS4A","NS4B","NS5A","NS5B"))

sumG$Subtype[sumG$Subtype=="1A"]<-"1a"
sumG$Subtype[sumG$Subtype=="1B"]<-"1b"
sumG$Subtype[sumG$Subtype=="3A"]<-"3a"
mf1$Subtype[mf1$Subtype=="1A"]<-"1a"
mf1$Subtype[mf1$Subtype=="1B"]<-"1b"
mf1$Subtype[mf1$Subtype=="3A"]<-"3a"
colnames(mf1)[c(2,5)]<-c("Gene","Subtype")

### Wilcoxon Test 
sums<-Sum[,c("merged.pos", "Gene", "mean.all.1A","mean.all.1B", "mean.all.3A")]
genenames<-unique(genes$Gene)
Comb<-t(combn(1:3, 2))
combnames<-c("1a-1b","1a-3a","1b-3a")

wilcox.res<-data.frame(matrix(nrow=12, ncol=3))
rownames(wilcox.res)<-genenames[1:12]
colnames(wilcox.res)<-combnames                       
for (i in 1:3){
        k=2
        n1<-Comb[i,1]+k
        n2<-Comb[i,2]+k
        for (g in 1:12){
                r1<-wilcox.test(sums[sums$Gene==genenames[g],n1], sums[sums$Gene==genenames[g],n2], alternative = "greater", paired = FALSE) 
                wilcox.res[g,i]<-r1[[3]][1]
        }
}

write.csv(wilcox.res,"Output_all/MVF/WilcoxTest_all.byGene.csv",row.names = T)


Test.results<-data.frame(gene=genenames[1:12])
for (i in 1:12){
        symb1<-NA
        symb2<-NA
        symb3<-NA
        if (wilcox.res$`1a-1b`[i]<0.05) symb1="a"
        if (wilcox.res$`1a-3a`[i]<0.05) symb2="b"
        if (wilcox.res$`1b-3a`[i]<0.05) symb3="c"
        
        symb<-c(symb1, symb2, symb3) 
        symb<-symb[!is.na(symb)]
        symb<-paste(symb, collapse = ',' )
        if (length(symb)==0) symb<-''
        Test.results$Sig[i]<-symb
}

## plot with significance info
sumG$Gene<-factor(sumG$Gene, levels=c("5' UTR","Core","E1", "HVR1","E2","NS1", "NS2","NS3","NS4A","NS4B","NS5A","NS5B"))
mf1$Gene<-factor(mf1$Gene, levels=c("5' UTR","Core","E1", "HVR1","E2","NS1", "NS2","NS3","NS4A","NS4B","NS5A","NS5B"))

ggplot()+
        geom_boxplot(data=mf1, aes(x=Gene, y=mean, fill=Subtype, color=Subtype),outlier.size=0.3,outlier.alpha = 0.2 )+
        scale_fill_manual(values=paste0(colors2[c(1,3,5)], "CC"))+
        scale_color_manual(values=colors2[c(1,3,5)])+
        scale_y_continuous(trans = 'log10', labels=label_scientific)+
        geom_errorbar(data=sumG, aes(x=Gene, y=Mean, ymin=Mean-SE, ymax=Mean+SE, group=Subtype), position=position_dodge(width=0.75),  width=.3)+
        geom_point(data=sumG, aes(x=Gene, y=Mean, group=Subtype), position=position_dodge(width=0.75), pch=21, color="gray20", bg=rep(colors2[c(1,3,5)], times=12))+
        theme_bw()+theme(axis.text=element_text(size=11), axis.title=element_text(size=13))+ylab("Minor variant frequency")+
        theme(panel.grid.major.x=element_blank(),axis.title.y = element_text(size=13))+
        geom_vline(xintercept = c(1:11)+0.5,  color = "gray70", size=.4)+
        theme(axis.title.x=element_blank())+
        annotate(geom="text", x=1:12, y=0.2,  label=Test.results$Sig,, color="red", size=2.5 )
ggsave(filename="Output_all/Figures/MVF.byGene.bySubtype_box_stats.pdf", width = 8.5, height = 5)
