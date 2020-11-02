
#load relevant libraries and read consensusfasta file

library(ape)
library(seqinr)
library(pegas)
library(sfsmisc)
library(ggplot2)
library(scales)
library(plotrix)
library(RColorBrewer)


#* Transition function*
transition<-function(nuc){
    if (is.na(nuc)|nuc=="-"){ 
                return (NA)
                next}
    else if (nuc=="a") return("g")
    else if (nuc=="g") return("a")
    else if (nuc=="c") return("t")
    else if (nuc=="t") return("c")
    else if (nuc=='n') return('n')
}

transv1<-function(nuc){
        if (is.na(nuc) ==T ){ 
                return (NA)
                next}
        else if (nuc=="a") return("c")
        else if (nuc=="g") return("c")
        else if (nuc=="c") return("a")
        else if (nuc=="t") return("a")
        else if (nuc=='n') return('n')
}

transv2<-function(nuc){
        if (is.na(nuc) ==T ){ 
                return ("NA")
                next}
        else if (nuc=="a") return("t")
        else if (nuc=="g") return("t")
        else if (nuc=="c") return("g")
        else if (nuc=="t") return("g")
        else if (nuc=='n') return('n')
}






typeofsitefunction<-function(WTcodon, mutantcodon){
        WTAA<-seqinr::translate(WTcodon)
        MUTAA<-seqinr::translate(mutantcodon)
        if (WTAA == "X"){ 
                return ("NA")
                next}
        else if (WTAA == MUTAA) return ("syn")
        else if (MUTAA == "*") return ("stop")
        else return ("nonsyn")
}


EstimatedS <- function(mu, meanfreq){
        if (is.na(mu)|is.na(meanfreq) ==T){ 
                return ("NA")
                next}
         else if (meanfreq == 0) return (1)
         else return (min(c(mu/meanfreq,1)))
}

#Amino acid changes
pos <- "R|H|K"
neg <- "D|E"
unc <- "S|T|N|Q"
spe <- "C|U|G|P"
hyd <- "A|I|L|F|M|W|Y|V"
amCat <- function(AA){
    if(regexpr(pos, AA) > 0){ return(0) }
    if(regexpr(neg, AA) > 0){ return(1) }
    if(regexpr(unc, AA) > 0){ return(2) }
    if(regexpr(spe, AA) > 0){ return(3) }
    if(regexpr(hyd, AA) > 0){ return(4) }
    return(5)
}

#Colors
#3 colors for Genotypes)
cols2<-c("#66CCEE" ,"#EE6677", "#228833")

