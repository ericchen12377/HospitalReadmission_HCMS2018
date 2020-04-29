#This is the package for COX model and C-index
library(survival)
#This is the package for Parametric Hazard model 
library(parfm)
#These two packages below are for C-index
library(Hmisc)
library(survcomp)


# This function is used to get C-index for each different model, using uncensored data 
# There are three packages can solve for C-index, with slight difference. 
# The results of this function will return three values of C-index for each model.
GetCindex.uncensored = function(predict,times){
  require(survcomp)
  require(Hmisc)
  require(survival)
  status = rep(1,length(predict))
  c.index = concordance.index(x=predict, surv.time=times, surv.event=status,method="noether")$c.index
  cindex.orig=1-rcorr.cens(predict,Surv(times))
  surv.c = survConcordance(Surv(times) ~ predict)$concordance
  p = c()
  p[1] = c.index
  p[2] = cindex.orig[1]
  p[3] = surv.c
  p = t(as.matrix(p))
  colnames(p)<-c('No tie risk','rcorrcens','Tie risk')
  return(p)
}

#Model 1: COX model with random effects in patients
#This is the function to generate risk score for COX model, both train and test data
#The result returns a list 
coxphrisk_p = function(fit,train_data,test_data,dc){
  risk = list()
  train_data = train_data[with(train_data, order(train_data$VisitLink)), ]
  test_data = test_data[with(test_data, order(test_data$VisitLink)), ]
  train_ID = as.factor(train_data$VisitLink)
  train_ID = as.integer(train_ID)
  train_data = cbind(train_ID,train_data)
  train_frail = c()
  test_frail = c()
  for( i in 1:length(train_data[,1])){
    train_frail[i] = fit$frail[train_ID[i]]
    for(j in 1:length(test_data[,1])){
      if(test_data$VisitLink[j]==train_data$VisitLink[i]){
        test_frail[j]=train_frail[i]
      }
    }
  }
  tfrail = test_frail
  tfrail[is.na(tfrail)]=0
  test_frail[is.na(test_frail)]= mean(tfrail)
  risk1 = exp(fit$linear.predictors+train_frail)
  lp = c()
  risk2 = c()
  fit$coefficients[is.na(fit$coefficients)]=0
  for( k in 1:length(test_data[,1])){
    lp[k] = as.matrix(test_data[k,dc])%*%as.matrix(fit$coefficients)
    risk2[k] = exp(lp[k]+test_frail[k])
  }
  risk[[1]]<-risk1
  risk[[2]]<-risk2
  return(risk)
}

#This function is used to generate C-index and risk scores for both train and test data, 
#and estimation of random effect variables, which is frailty. 
#The results are returned in a list
gammaFrailty.up = function(train_data,test_data){
  result.cox.frailty <- coxph(Surv(train_data$DtoA) ~.+frailty(train_data$VisitLink,distribution = "gamma"),data = train_data[,-c(1:5)])
  print(summary(result.cox.frailty))
  risk = coxphrisk_p(result.cox.frailty,train_data,test_data,dc = c(6:40))
  cindex = list()
  cindex[[1]] = GetCindex.uncensored(risk[[1]], train_data$DtoA)
  cindex[[2]] = GetCindex.uncensored(risk[[2]], test_data$DtoA)
  cindex[[3]] = result.cox.frailty$frail
  cindex[[4]] = result.cox.frailty
  cindex[[5]] = risk
  return(cindex)
}

#Model 2: Parametric regression model with exponential baseline function, no random effect
#This is the function for model 2, with results of c-index, risk scores for train and test
#The results are returned as a list
#################Parametric regression###########################
expreg.u = function(train_data,test_data){
  ######weibull#####
  surv<-survreg(Surv(train_data$DtoA)~.,data=train_data[,-c(1:5)],dist="exponential")
  print(summary(surv))
  cindex = list()
  #####c-index for train######
  risk = exp(predict(surv, train_data, type="link"))  
  cindex[[1]] = GetCindex.uncensored(risk, train_data$DtoA)
  #####c-index for test######
  risk = exp(predict(surv, test_data, type="link"))   
  cindex[[2]] = GetCindex.uncensored(risk, test_data$DtoA)
  cindex[[3]] = surv 
  return(cindex)
}



