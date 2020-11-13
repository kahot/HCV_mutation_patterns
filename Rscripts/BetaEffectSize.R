#result=glm or betareg result (the output of betareg())
BetaEffectSize<-function(result){
        res<-summary(result)
        model<-data.frame(res$coefficients[[1]])
        model$Effect<-''
        for (i in 1:length(row.names(model)) ){
                if (i==1) model$Effect[1]<- exp(model[1,i])
                
                else{
                        model$Effect[i]<- (((exp(model[1,1] + model$Estimate[i]) - exp(model[1,1])) /exp(model[1,1]))) }
        }
        return(model)
}

BetaEffectSizeDF<-function(model){
        model$estMF<-''
        model$Effect<-''
        for (i in 1:length(row.names(model)) ){
                if (i==1) {model$estMF[1]<- exp(model[1,i])
                           model$Effect[1]<-0}
                
                else{   model$estMF[i]<-exp(model[1,1] + model$Estimate[i])
                        model$Effect[i]<- (((exp(model[1,1] + model$Estimate[i]) - exp(model[1,1])) /exp(model[1,1]))) }
        }
        return(model)
}
