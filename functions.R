#all the libraries used for survival models
library(survival)
library(coxme)
library(frailtypack)
library(Hmisc)
library(survcomp)
library(bujar)
library(glmnet)
library(randomForest)
library(gbm)
library(tree)
library(parfm)
library(gdata)
library(e1071)
library(neuralnet)



#C-index function# combining three different packages library(Hmisc), library(survcomp), library(survival)
GetCindex = function(predict,times,status){
  require(survcomp)
  require(Hmisc)
  require(survival)
  c.index = concordance.index(x=predict, surv.time=times, surv.event=status, method="noether")$c.index
  cindex.orig=1-rcorr.cens(predict,Surv(times, status))
  surv.c = survConcordance(Surv(times, status) ~ predict)$concordance
  p = c()
  p[1] = c.index
  p[2] = cindex.orig[1]
  p[3] = surv.c
  p = t(as.matrix(p))
  colnames(p)<-c('No tie risk','rcorrcens','Tie risk')
return(p)
}
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


#coxme risk prediction###########
coxmerisk = function(fit,data,value,dc){
  num = length(fit$frail)
  f = list()
  rname = list()
  rnamevalue = list()
  frail<-list()
  tfrail<-rep(0,length(data[,1]))
  for ( i in 1:num){
    f[[i]] = as.matrix(fit$frail[[i]])
    rname[[i]] = row.names(f[[i]])
    rvalue<-c()
    for( j in 1:length(rname[[i]])){
      rvalue[j] = eval(parse(text = rname[[i]][j]))
    }
    rnamevalue[[i]] = rvalue
    rfrail = c()
    for( k in 1:length(data[,1])){
      for(j in 1:length(rname[[i]])){
        if(value[[i]][k] == rvalue[j]){
          rfrail[k] = f[[i]][j]
        }
      }
    }
    rrfrail = rfrail
    rrfrail[is.na(rrfrail)]=0
    rfrail[is.na(rfrail)]= mean(rrfrail)
    frail[[i]] = rfrail
    tfrail = tfrail+frail[[i]]
  }
  lp<-c()
  risk<-c()
  fit$coefficients[is.na(fit$coefficients)]=0
  for( k in 1:length(data[,1])){
    lp[k] = as.matrix(data[k,dc])%*%as.matrix(fit$coefficients)
    risk[k] = exp(lp[k]+tfrail[k])
  }
return(risk)
}
coxmerisk_nf = function(fit,data,value,dc){
  lp<-c()
  risk<-c()
  fit$coefficients[is.na(fit$coefficients)]=0
  for( k in 1:length(data[,1])){
    lp[k] = as.matrix(data[k,dc])%*%as.matrix(fit$coefficients)
    risk[k] = exp(lp[k])
  }
  return(risk)
}
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
coxphrisk_r = function(fit,train_data,test_data,dc){
  risk = list()
  train_data = train_data[with(train_data, order(train_data$PL_NCHS2006)), ]
  test_data = test_data[with(test_data, order(test_data$PL_NCHS2006)), ]
  train_ID = as.factor(train_data$PL_NCHS2006)
  train_ID = as.integer(train_ID)
  train_data = cbind(train_ID,train_data)
  train_frail = c()
  test_frail = c()
  for( i in 1:length(train_data[,1])){
    train_frail[i] = fit$frail[train_ID[i]]
    for(j in 1:length(test_data[,1])){
      if(test_data$PL_NCHS2006[j]==train_data$PL_NCHS2006[i]){
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
coxphrisk_h = function(fit,train_data,test_data,dc){
  risk = list()
  train_data = train_data[with(train_data, order(train_data$DSHOSPID)), ]
  test_data = test_data[with(test_data, order(test_data$DSHOSPID)), ]
  train_ID = as.factor(train_data$DSHOSPID)
  train_ID = as.integer(train_ID)
  train_data = cbind(train_ID,train_data)
  train_frail = c()
  test_frail = c()
  for( i in 1:length(train_data[,1])){
    train_frail[i] = fit$frail[train_ID[i]]
    for(j in 1:length(test_data[,1])){
      if(test_data$DSHOSPID[j]==train_data$DSHOSPID[i]){
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
#############Frailty model - patient###########
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

gaussianFrailty.up = function(train_data,test_data){
  result.cox.frailty <- coxph(Surv(train_data$DtoA) ~.+frailty(train_data$VisitLink,distribution = "gaussian"),data = train_data[,-c(1:5)])
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

tFrailty.up = function(train_data,test_data){
  result.cox.frailty <- coxph(Surv(train_data$DtoA) ~.+frailty(train_data$VisitLink,distribution = "t"),data = train_data[,-c(1:5)])
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

#############Frailty model - region###########
gammaFrailty.ur = function(train_data,test_data){
  result.cox.frailty <- coxph(Surv(train_data$DtoA) ~.+frailty(train_data$PL_NCHS2006,distribution = "gamma"),data = train_data[,-c(1:5)])
  print(summary(result.cox.frailty))
  risk = coxphrisk_r(result.cox.frailty,train_data,test_data,dc = c(6:40))
  cindex = list()
  cindex[[1]] = GetCindex.uncensored(risk[[1]], train_data$DtoA)
  cindex[[2]] = GetCindex.uncensored(risk[[2]], test_data$DtoA)
  cindex[[3]] = result.cox.frailty$frail
  cindex[[4]] = result.cox.frailty
  cindex[[5]] = risk
  return(cindex)
}

gaussianFrailty.ur = function(train_data,test_data){
  result.cox.frailty <- coxph(Surv(train_data$DtoA) ~.+frailty(train_data$PL_NCHS2006,distribution = "gaussian"),data = train_data[,-c(1:5)])
  print(summary(result.cox.frailty))
  risk = coxphrisk_r(result.cox.frailty,train_data,test_data,dc = c(6:40))
  cindex = list()
  cindex[[1]] = GetCindex.uncensored(risk[[1]], train_data$DtoA)
  cindex[[2]] = GetCindex.uncensored(risk[[2]], test_data$DtoA)
  cindex[[3]] = result.cox.frailty$frail
  cindex[[4]] = result.cox.frailty
  cindex[[5]] = risk
  return(cindex)
}

tFrailty.ur = function(train_data,test_data){
  result.cox.frailty <- coxph(Surv(train_data$DtoA) ~.+frailty(train_data$PL_NCHS2006,distribution = "t"),data = train_data[,-c(1:5)])
  print(summary(result.cox.frailty))
  risk = coxphrisk_r(result.cox.frailty,train_data,test_data,dc = c(6:40))
  cindex = list()
  cindex[[1]] = GetCindex.uncensored(risk[[1]], train_data$DtoA)
  cindex[[2]] = GetCindex.uncensored(risk[[2]], test_data$DtoA)
  cindex[[3]] = result.cox.frailty$frail
  cindex[[4]] = result.cox.frailty
  cindex[[5]] = risk
  return(cindex)
}

#############Frailty model - hosptial###########
gammaFrailty.uh = function(train_data,test_data){
  result.cox.frailty <- coxph(Surv(train_data$DtoA) ~.+frailty(train_data$DSHOSPID,distribution = "gamma"),data = train_data[,-c(1:5)])
  print(summary(result.cox.frailty))
  risk = coxphrisk_h(result.cox.frailty,train_data,test_data,dc = c(6:40))
  cindex = list()
  cindex[[1]] = GetCindex.uncensored(risk[[1]], train_data$DtoA)
  cindex[[2]] = GetCindex.uncensored(risk[[2]], test_data$DtoA)
  cindex[[3]] = result.cox.frailty$frail
  cindex[[4]] = result.cox.frailty
  cindex[[5]] = risk
  return(cindex)
}

gaussianFrailty.uh = function(train_data,test_data){
  result.cox.frailty <- coxph(Surv(train_data$DtoA) ~.+frailty(train_data$DSHOSPID,distribution = "gaussian"),data = train_data[,-c(1:5)])
  print(summary(result.cox.frailty))
  risk = coxphrisk_h(result.cox.frailty,train_data,test_data,dc = c(6:40))
  cindex = list()
  cindex[[1]] = GetCindex.uncensored(risk[[1]], train_data$DtoA)
  cindex[[2]] = GetCindex.uncensored(risk[[2]], test_data$DtoA)
  cindex[[3]] = result.cox.frailty$frail
  cindex[[4]] = result.cox.frailty
  cindex[[5]] = risk
  return(cindex)
}

tFrailty.uh = function(train_data,test_data){
  result.cox.frailty <- coxph(Surv(train_data$DtoA) ~.+frailty(train_data$DSHOSPID,distribution = "t"),data = train_data[,-c(1:5)])
  print(summary(result.cox.frailty))
  risk = coxphrisk_h(result.cox.frailty,train_data,test_data,dc = c(6:40))
  cindex = list()
  cindex[[1]] = GetCindex.uncensored(risk[[1]], train_data$DtoA)
  cindex[[2]] = GetCindex.uncensored(risk[[2]], test_data$DtoA)
  cindex[[3]] = result.cox.frailty$frail
  cindex[[4]] = result.cox.frailty
  cindex[[5]] = risk
  return(cindex)
}
gamma_no_Frailty.u = function(train_data,test_data){
  result.cox.frailty <- coxph(Surv(train_data$DtoA) ~.+frailty(train_data$VisitLink,distribution = "gamma"),data = train_data[,-c(1:5)])
  print(summary(result.cox.frailty))
  risk1 = predict(result.cox.frailty, train_data,type=c("risk"))
  risk2 = predict(result.cox.frailty, test_data,type=c("risk"))
  cindex_nf = list()
  cindex_nf[[1]] = GetCindex.uncensored(risk1, train_data$DtoA)
  cindex_nf[[2]] = GetCindex.uncensored(risk2, test_data$DtoA)
  return(cindex_nf)
}
gaussian_no_Frailty.u = function(train_data,test_data){
  result.cox.frailty <- coxph(Surv(train_data$DtoA) ~.+frailty(train_data$VisitLink,distribution = "gaussian"),data = train_data[,-c(1:5)])
  print(summary(result.cox.frailty))
  risk1 = predict(result.cox.frailty, train_data,type=c("risk"))
  risk2 = predict(result.cox.frailty, test_data,type=c("risk"))
  cindex_nf = list()
  cindex_nf[[1]] = GetCindex.uncensored(risk1, train_data$DtoA)
  cindex_nf[[2]] = GetCindex.uncensored(risk2, test_data$DtoA)
  return(cindex_nf)
}
t_no_Frailty.u = function(train_data,test_data){
  result.cox.frailty <- coxph(Surv(train_data$DtoA) ~.+frailty(train_data$VisitLink,distribution = "t"),data = train_data[,-c(1:5)])
  print(summary(result.cox.frailty))
  risk1 = predict(result.cox.frailty, train_data,type=c("risk"))
  risk2 = predict(result.cox.frailty, test_data,type=c("risk"))
  cindex_nf = list()
  cindex_nf[[1]] = GetCindex.uncensored(risk1, train_data$DtoA)
  cindex_nf[[2]] = GetCindex.uncensored(risk2, test_data$DtoA)
  return(cindex_nf)
}

######################COX,LASSO-COX, EN-LASSO-COX##############
COX.u = function(train_data,test_data){
  result.cox <- coxph(Surv(train_data$DtoA) ~.,data = train_data[,-c(1:5)])
  print(summary(result.cox))
  cindex = list()
  #####c-index for train######
  risk = exp(result.cox$linear.predictors)
  cindex[[1]] = GetCindex.uncensored(risk, train_data$DtoA)
  #####c-index for test######
  risk = predict(result.cox,test_data,type = 'risk')
  cindex[[2]] = GetCindex.uncensored(risk, test_data$DtoA)
  cindex[[3]] = result.cox
  return(cindex)
}
LASSO_COX.u = function(train_data,test_data){
  x_train <- as.matrix(train_data[,-c(1:5)])
  x_test <- as.matrix(test_data[,-c(1:5)])
  y <- Surv(train_data$DtoA)
  cv.fit <- cv.glmnet(x_train, y, family="cox", alpha=1)
  print(coef(cv.fit))
  cindex = list()
  #####c-index for train######
  risk = exp(predict(cv.fit, newx=x_train, type="link"))  
  cindex[[1]] = GetCindex.uncensored(risk, train_data$DtoA)
  #####c-index for test######
  risk = exp(predict(cv.fit, newx=x_test, type="link")) 
  cindex[[2]] = GetCindex.uncensored(risk, test_data$DtoA)
  cindex[[3]] = coef(cv.fit)
  return(cindex)
}
EN_LASSO_COX.u = function(train_data,test_data){
  x_train <- as.matrix(train_data[,-c(1:5)])
  x_test <- as.matrix(test_data[,-c(1:5)])
  y <- Surv(train_data$DtoA)
  cv.fit <- cv.glmnet(x_train, y, family="cox", alpha=0.1)
  print(coef(cv.fit))
  cindex = list()
  #####c-index for train######
  risk = exp(predict(cv.fit, newx=x_train, type="link"))  
  cindex[[1]] = GetCindex.uncensored(risk, train_data$DtoA)
  #####c-index for test######
  risk = exp(predict(cv.fit, newx=x_test, type="link")) 
  cindex[[2]] = GetCindex.uncensored(risk, test_data$DtoA)
  cindex[[3]] = coef(cv.fit) 
  return(cindex)
}
###############################################################

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
weibullreg.u = function(train_data,test_data){
  ######weibull#####
  surv<-survreg(Surv(train_data$DtoA)~.,data=train_data[,-c(1:5)],dist="weibull")
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
logisticreg.u = function(train_data,test_data){

  surv<-survreg(Surv(train_data$DtoA)~.,data=train_data[,-c(1:5)],dist="logistic")
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
loglogisticreg.u = function(train_data,test_data){

  surv<-survreg(Surv(train_data$DtoA)~.,data=train_data[,-c(1:5)],dist="loglogistic")
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
lognormalreg.u = function(train_data,test_data){

  surv<-survreg(Surv(train_data$DtoA)~.,data=train_data[,-c(1:5)],dist="lognormal")
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
##########################linear models###############
########linear regression
linearreg.u = function(train_data,test_data){
  lm<-lm(train_data$DtoA~., data=train_data[,-(1:5)])
  print(summary(lm))
  cindex = list()
  #####c-index for train######
  risk = predict(lm, train_data, type="response")
  cindex[[1]] = GetCindex.uncensored(risk, train_data$DtoA)
  #####c-index for test######
  risk = predict(lm, test_data, type="response")   
  cindex[[2]] = GetCindex.uncensored(risk, test_data$DtoA)
  cindex[[3]] = lm 
  return(cindex)
}
tobit.u = function(train_data,test_data){
  surv<-survreg(Surv(train_data$DtoA)~.,data=train_data[,-c(1:5)],dist="gaussian")
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
###############other uncensored methods################
LASSO.u = function(train_data,test_data){
  x_train <- as.matrix(train_data[,-c(1:5)])
  x_test <- as.matrix(test_data[,-c(1:5)])
  y <- train_data$DtoA
  cv.fit <- cv.glmnet(x_train, y,alpha=1)
  print(coef(cv.fit))
  cindex = list()
  #####c-index for train######
  risk = exp(predict(cv.fit, newx=x_train, type="link"))  
  cindex[[1]] = GetCindex.uncensored(risk, train_data$DtoA)
  #####c-index for test######
  risk = exp(predict(cv.fit, newx=x_test, type="link")) 
  cindex[[2]] = GetCindex.uncensored(risk, test_data$DtoA)
  cindex[[3]] = coef(cv.fit) 
  return(cindex)
}
Ridge.u = function(train_data,test_data){
  x_train <- as.matrix(train_data[,-c(1:5)])
  x_test <- as.matrix(test_data[,-c(1:5)])
  y <- train_data$DtoA
  cv.fit <- cv.glmnet(x_train, y,alpha=0)
  print(coef(cv.fit))
  cindex = list()
  #####c-index for train######
  risk = exp(predict(cv.fit, newx=x_train, type="link"))  
  cindex[[1]] = GetCindex.uncensored(risk, train_data$DtoA)
  #####c-index for test######
  risk = exp(predict(cv.fit, newx=x_test, type="link")) 
  cindex[[2]] = GetCindex.uncensored(risk, test_data$DtoA)
  cindex[[3]] = coef(cv.fit)  
  return(cindex)
}
####################classification methods########
rforest.u = function(train_data,test_data,mtry,ntree){
  rforestmodel=randomForest(train_data$DtoA~.,data=train_data[,-c(1:5)],mtry = mtry, ntree=ntree,importance = TRUE )
  cindex = list()
  #####c-index for train######
  risk = rforestmodel$predicted 
  cindex[[1]] = GetCindex.uncensored(risk, train_data$DtoA)
  #####c-index for test######
  risk = predict(rforestmodel, test_data, type="response") 
  cindex[[2]] = GetCindex.uncensored(risk, test_data$DtoA)
  cindex[[3]] = rforestmodel 
  return(cindex)
}
boostingtree.u = function(train_data,test_data,dist,depth,shrinkage){
  boostmodel=gbm(train_data$DtoA~., data=train_data[,-c(1:5)], distribution = "gaussian", interaction.depth = 3, shrinkage = 0.15)
  
  cindex = list()
  #####c-index for train######
  risk = boostmodel$fit  
  cindex[[1]] = GetCindex.uncensored(risk, train_data$DtoA)
  #####c-index for test######
  risk = predict(boostmodel, test_data, n.trees =100,type="link") 
  cindex[[2]] = GetCindex.uncensored(risk, test_data$DtoA)
  cindex[[3]] = boostmodel  
  return(cindex)
}
tree.u = function(train_data,test_data){
  tree.model<-tree(DtoA~.,data=train_data[,-c(1:5)])
  cindex = list()
  risk = predict(tree.model,train_data) 
  cindex[[1]] = GetCindex.uncensored(risk, train_data$DtoA)
  #####c-index for test######
  risk = predict(tree.model, test_data) 
  cindex[[2]] = GetCindex.uncensored(risk, test_data$DtoA)
  cindex[[3]] = tree.model  
  return(cindex)
}

svm.u = function(train_data,test_data){
  svm <- svm(train_data$DtoA ~ ., data = train_data[,-c(1:5)], type='eps-regression', kernel="linear",scale = FALSE)
  cindex = list()
  risk = predict(svm,train_data) 
  cindex[[1]] = GetCindex.uncensored(risk, train_data$DtoA)
  #####c-index for test######
  risk = predict(svm, test_data) 
  cindex[[2]] = GetCindex.uncensored(risk, test_data$DtoA)
  cindex[[3]] = svm  
  return(cindex)
}


##################Cox mixed effect###########
##################Parametric Frailty model############################
#####################################Censored methods############################################################
#############Frailty model###########
gammaFrailty= function(train_data,test_data){
  result.cox.frailty <- coxph(Surv(train_data$DtoA,train_data$status) ~.+frailty(train_data$VisitLink,distribution = "gamma"),data = train_data[,-c(1:5)])
  print(summary(result.cox.frailty))
  cindex = list()
  risk = coxphrisk(result.cox.frailty,train_data,test_data,dc = c(6:40))
  #####c-index for train######
  cindex[[1]] = GetCindex(risk[[1]], train_data$DtoA,train_data$status)
  #####c-index for test######
  cindex[[2]] = GetCindex(risk[[2]], test_data$DtoA,train_data$status)
   
  return(cindex)
}
gaussianFrailty = function(train_data,test_data){
  result.cox.frailty <- coxph(Surv(train_data$DtoA,train_data$status) ~.+frailty(train_data$VisitLink,distribution = "gaussian"),data = train_data[,-c(1:5)])
  print(summary(result.cox.frailty))
  cindex = list()
  risk = coxphrisk(result.cox.frailty,train_data,test_data,dc = c(6:40))
  #####c-index for train######
  cindex[[1]] = GetCindex(risk[[1]], train_data$DtoA,train_data$status)
  #####c-index for test######
  cindex[[2]] = GetCindex(risk[[2]], test_data$DtoA,train_data$status)
   
  return(cindex)
}
tFrailty = function(train_data,test_data){
  result.cox.frailty <- coxph(Surv(train_data$DtoA,train_data$status) ~.+frailty(train_data$VisitLink,distribution = "t"),data = train_data[,-c(1:5)])
  print(summary(result.cox.frailty))
  cindex = list()
  risk = coxphrisk(result.cox.frailty,train_data,test_data,dc = c(6:40))
  #####c-index for train######
  cindex[[1]] = GetCindex(risk[[1]], train_data$DtoA,train_data$status)
  #####c-index for test######
  cindex[[2]] = GetCindex(risk[[2]], test_data$DtoA,train_data$status)
   
  return(cindex)
}
######################COX,LASSO-COX, EN-LASSO-COX##############
COX = function(train_data,test_data){
  result.cox <- coxph(Surv(train_data$DtoA,train_data$status) ~.,data = train_data[,-c(1:5)])
  print(summary(result.cox))
  cindex = list()
  #####c-index for train######
  risk = exp(result.cox$linear.predictors)
  cindex[[1]] = GetCindex(risk[[1]], train_data$DtoA,train_data$status)
  #####c-index for test######
  risk = predict(result.cox,test_data,type = 'risk')
  cindex[[2]] = GetCindex(risk[[2]], test_data$DtoA,train_data$status)
  cindex[[3]] = result.cox
  return(cindex)
}
LASSO_COX = function(train_data,test_data){
  x_train <- as.matrix(train_data[,-c(1:5)])
  x_test <- as.matrix(test_data[,-c(1:5)])
  y <- Surv(train_data$DtoA)
  cv.fit <- cv.glmnet(x_train, y, family="cox", alpha=1)
  print(coef(cv.fit))
  cindex = list()
  #####c-index for train######
  risk = exp(predict(cv.fit, newx=x_train, type="link"))  
  cindex[[1]] = GetCindex(risk[[1]], train_data$DtoA,train_data$status)
  #####c-index for test######
  risk = exp(predict(cv.fit, newx=x_test, type="link")) 
  cindex[[2]] = GetCindex(risk[[2]], test_data$DtoA,train_data$status)
  cindex[[3]] = cv.fit
  return(cindex)
}
EN_LASSO_COX = function(train_data,test_data){
  x_train <- as.matrix(train_data[,-c(1:5)])
  x_test <- as.matrix(test_data[,-c(1:5)])
  y <- Surv(train_data$DtoA)
  cv.fit <- cv.glmnet(x_train, y, family="cox", alpha=0.1)
  print(coef(cv.fit))
  cindex = list()
  #####c-index for train######
  risk = exp(predict(cv.fit, newx=x_train, type="link"))  
  cindex[[1]] = GetCindex(risk[[1]], train_data$DtoA,train_data$status)
  #####c-index for test######
  risk = exp(predict(cv.fit, newx=x_test, type="link")) 
  cindex[[2]] = GetCindex(risk[[2]], test_data$DtoA,train_data$status)
  cindex[[3]] = cv.fit 
  return(cindex)
}
###############################################################

#################Parametric regression###########################
weibullreg = function(train_data,test_data){
  ######weibull#####
  surv<-survreg(Surv(train_data$DtoA,train_data$status)~.,data=train_data[,-c(1:5)],dist="weibull")
  print(summary(surv))
  cindex = list()
  #####c-index for train######
  risk = exp(predict(surv, train_data, type="link"))  
  cindex[[1]] = GetCindex(risk[[1]], train_data$DtoA,train_data$status)
  #####c-index for test######
  risk = exp(predict(surv, test_data, type="link"))   
  cindex[[2]] = GetCindex(risk[[2]], test_data$DtoA,train_data$status)
  cindex[[3]] = surv 
  return(cindex)
}
logisticreg = function(train_data,test_data){
  
  surv<-survreg(Surv(train_data$DtoA,train_data$status)~.,data=train_data[,-c(1:5)],dist="logistic")
  print(summary(surv))
  cindex = list()
  #####c-index for train######
  risk = exp(predict(surv, train_data, type="link"))  
  cindex[[1]] = GetCindex(risk[[1]], train_data$DtoA,train_data$status)
  #####c-index for test######
  risk = exp(predict(surv, test_data, type="link"))   
  cindex[[2]] = GetCindex(risk[[2]], test_data$DtoA,train_data$status)
  cindex[[3]] = surv  
  return(cindex)
}
loglogisticreg = function(train_data,test_data){
  
  surv<-survreg(Surv(train_data$DtoA,train_data$status)~.,data=train_data[,-c(1:5)],dist="loglogistic")
  print(summary(surv))
  cindex = list()
  #####c-index for train######
  risk = exp(predict(surv, train_data, type="link"))  
  cindex[[1]] = GetCindex(risk[[1]], train_data$DtoA,train_data$status)
  #####c-index for test######
  risk = exp(predict(surv, test_data, type="link"))   
  cindex[[2]] = GetCindex(risk[[2]], test_data$DtoA,train_data$status)
  cindex[[3]] = surv  
  return(cindex)
}
lognormalreg = function(train_data,test_data){
  
  surv<-survreg(Surv(train_data$DtoA,train_data$status)~.,data=train_data[,-c(1:5)],dist="lognormal")
  print(summary(surv))
  cindex = list()
  #####c-index for train######
  risk = exp(predict(surv, train_data, type="link"))  
  cindex[[1]] = GetCindex(risk[[1]], train_data$DtoA,train_data$status)
  #####c-index for test######
  risk = exp(predict(surv, test_data, type="link"))   
  cindex[[2]] = GetCindex(risk[[2]], test_data$DtoA,train_data$status)
  cindex[[3]] = surv  
  return(cindex)
}





##################################Parametric frailty model########################################
#################################                         ########################################
base_exp = function(data,mod){
  mod[2,1]*rep(1,length(data[,1]))
}
base_weibull = function(data,mod){
  mod[2,1]*mod[3,1]*(data$DtoA^(mod[2,1]-1))
}
base_gompertz = function(data,mod){
  mod[3,1]*exp(mod[2,1]*data$DtoA)
}
base_lognormal = function(data,mod){
  z = (log(data$DtoA)-mod[2,1])/mod[3,1]
  base = dnorm(z)/(mod[3,1]*data$DtoA*(1-pnorm(z)))
  return(base)
}
base_loglogistic = function(data,mod){
  num = exp(mod[2,1])*mod[3,1]*(data$DtoA^(mod[3,1]-1))
  den = 1+exp(mod[2,1])*(data$DtoA^mod[3,1])
  base = num/den
  return(base)
}


parfmrisk = function(fit,train_data,test_data,train_cluster,test_cluster,base,numbase){
  risk = list()
  train_ID = as.factor(train_cluster)
  train_ID = as.integer(train_ID)
  train_data = cbind(train_ID,train_data)
  train_frail = c()
  test_frail = c()
  fit_frail = predict(fit)
  for( i in 1:length(train_data[,1])){
    train_frail[i] = fit_frail[train_ID[i]]
    for(j in 1:length(test_data[,1])){
      if(test_cluster[j]==train_cluster[i]){
        test_frail[j]=train_frail[i]
      }
    }
  }
  tfrail = test_frail
  tfrail[is.na(tfrail)]=0
  test_frail[is.na(test_frail)]= mean(tfrail)
  coefficient = fit[,1][-(1:(1+numbase))]
  coefficient[is.na(coefficient)]=0
  risk1 = base[[1]]*exp(as.matrix(train_data[,-c(1:6)])%*%coefficient)*train_frail
  risk2 = base[[2]]*exp(as.matrix(test_data[,-c(1:5)])%*%coefficient)*test_frail
  risk[[1]]<-risk1
  risk[[2]]<-risk2
  risk[[3]]<-train_frail
  risk[[4]]<-test_frail
  return(risk)
}


##################################Hosptial####################################################
parfm_exp_gamma_hosptial = function(train_data,test_data){
  train_data = train_data[with(train_data, order(DSHOSPID)), ]
  test_data = test_data[with(test_data, order(DSHOSPID)), ]
mod = parfm(Surv(DtoA,status) ~ 
                LOS
              +AGE
              +ASOURCE2
              +ASOURCE3
              +ASOURCE4
              +ASOURCE5
              +ASCHED1
              +DISPUNIFORM2
              +DISPUNIFORM5
              +DISPUNIFORM6
              +DISPUNIFORM7
              +DISPUNIFORM20
              +DISPUNIFORM99
              +BMONTH
              +FEMALE1
              +HISPANIC_X2
              +HISPANIC_X3
              +HISPANIC_X4
              +RACE2
              +RACE3
              +RACE4
              +RACE5
              +RACE6
              +RACE7
              +AWEEKEND1
              +AMONTH
              +DMONTH
              +DNR1
              +DQTR
              +NCHRONIC
              +NDX
              +NECODE
              +NPR
              +TVisits
              +NVisit, cluster = "DSHOSPID",
              data = train_data, dist = "exponential", frailty = "gamma")
  base = list()
  base[[1]] = base_exp(train_data,mod)
  base[[2]] = base_exp(test_data,mod)
  train_cluster = train_data$DSHOSPID
  test_cluster = test_data$DSHOSPID
  risk = parfmrisk(mod,train_data,test_data,train_cluster,test_cluster,base,numbase = 1)
  cindex = list()
  cindex[[1]] = GetCindex.uncensored(risk[[1]], train_data$DtoA)
  cindex[[2]] = GetCindex.uncensored(risk[[2]], test_data$DtoA)
  cindex[[3]] = risk[[3]]
  cindex[[4]] = risk[[4]]
  cindex[[5]] = mod
  return(cindex)
}

parfm_weibull_gamma_hosptial = function(train_data,test_data){
  train_data = train_data[with(train_data, order(DSHOSPID)), ]
  test_data = test_data[with(test_data, order(DSHOSPID)), ]
  mod = parfm(Surv(DtoA,status) ~ 
                LOS
              +AGE
              +ASOURCE2
              +ASOURCE3
              +ASOURCE4
              +ASOURCE5
              +ASCHED1
              +DISPUNIFORM2
              +DISPUNIFORM5
              +DISPUNIFORM6
              +DISPUNIFORM7
              +DISPUNIFORM20
              +DISPUNIFORM99
              +BMONTH
              +FEMALE1
              +HISPANIC_X2
              +HISPANIC_X3
              +HISPANIC_X4
              +RACE2
              +RACE3
              +RACE4
              +RACE5
              +RACE6
              +RACE7
              +AWEEKEND1
              +AMONTH
              +DMONTH
              +DNR1
              +DQTR
              +NCHRONIC
              +NDX
              +NECODE
              +NPR
              +TVisits
              +NVisit, cluster = "DSHOSPID",
              data = train_data, dist = "weibull", frailty = "gamma")
  base = list()
  base[[1]] = base_weibull(train_data,mod)
  base[[2]] = base_weibull(test_data,mod)
  train_cluster = train_data$DSHOSPID
  test_cluster = test_data$DSHOSPID
  risk = parfmrisk(mod,train_data,test_data,train_cluster,test_cluster,base,numbase = 2)
  cindex = list()
  cindex[[1]] = GetCindex.uncensored(risk[[1]], train_data$DtoA)
  cindex[[2]] = GetCindex.uncensored(risk[[2]], test_data$DtoA)
  cindex[[3]] = risk[[3]]
  cindex[[4]] = risk[[4]]
  cindex[[5]] = mod
  return(cindex)
}




parfm_lognormal_gamma_hosptial = function(train_data,test_data){
  train_data = train_data[with(train_data, order(DSHOSPID)), ]
  test_data = test_data[with(test_data, order(DSHOSPID)), ]
  mod = parfm(Surv(DtoA,status) ~ 
                LOS
              +AGE
              +ASOURCE2
              +ASOURCE3
              +ASOURCE4
              +ASOURCE5
              +ASCHED1
              +DISPUNIFORM2
              +DISPUNIFORM5
              +DISPUNIFORM6
              +DISPUNIFORM7
              +DISPUNIFORM20
              +DISPUNIFORM99
              +BMONTH
              +FEMALE1
              +HISPANIC_X2
              +HISPANIC_X3
              +HISPANIC_X4
              +RACE2
              +RACE3
              +RACE4
              +RACE5
              +RACE6
              +RACE7
              +AWEEKEND1
              +AMONTH
              +DMONTH
              +DNR1
              +DQTR
              +NCHRONIC
              +NDX
              +NECODE
              +NPR
              +TVisits
              +NVisit, cluster = "DSHOSPID",
              data = train_data, dist = "lognormal", frailty = "gamma")
  base = list()
  base[[1]] = base_lognormal(train_data,mod)
  base[[2]] = base_lognormal(test_data,mod)
  train_cluster = train_data$DSHOSPID
  test_cluster = test_data$DSHOSPID
  risk = parfmrisk(mod,train_data,test_data,train_cluster,test_cluster,base,numbase = 2)
  cindex = list()
  cindex[[1]] = GetCindex.uncensored(risk[[1]], train_data$DtoA)
  cindex[[2]] = GetCindex.uncensored(risk[[2]], test_data$DtoA)
  cindex[[3]] = risk[[3]]
  cindex[[4]] = risk[[4]]
  cindex[[5]] = mod
  return(cindex)
}


parfm_loglogistic_gamma_hosptial = function(train_data,test_data){
  train_data = train_data[with(train_data, order(DSHOSPID)), ]
  test_data = test_data[with(test_data, order(DSHOSPID)), ]
  mod = parfm(Surv(DtoA,status) ~ 
                LOS
              +AGE
              +ASOURCE2
              +ASOURCE3
              +ASOURCE4
              +ASOURCE5
              +ASCHED1
              +DISPUNIFORM2
              +DISPUNIFORM5
              +DISPUNIFORM6
              +DISPUNIFORM7
              +DISPUNIFORM20
              +DISPUNIFORM99
              +BMONTH
              +FEMALE1
              +HISPANIC_X2
              +HISPANIC_X3
              +HISPANIC_X4
              +RACE2
              +RACE3
              +RACE4
              +RACE5
              +RACE6
              +RACE7
              +AWEEKEND1
              +AMONTH
              +DMONTH
              +DNR1
              +DQTR
              +NCHRONIC
              +NDX
              +NECODE
              +NPR
              +TVisits
              +NVisit, cluster = "DSHOSPID",
              data = train_data, dist = "loglogistic", frailty = "gamma")
  base = list()
  base[[1]] = base_loglogistic(train_data,mod)
  base[[2]] = base_loglogistic(test_data,mod)
  train_cluster = train_data$DSHOSPID
  test_cluster = test_data$DSHOSPID
  risk = parfmrisk(mod,train_data,test_data,train_cluster,test_cluster,base,numbase = 2)
  cindex = list()
  cindex[[1]] = GetCindex.uncensored(risk[[1]], train_data$DtoA)
  cindex[[2]] = GetCindex.uncensored(risk[[2]], test_data$DtoA)
  cindex[[3]] = risk[[3]]
  cindex[[4]] = risk[[4]]
  cindex[[5]] = mod
  return(cindex)
}

##############################Patient###################################################
parfm_exp_gamma_patient = function(train_data,test_data){
  train_data = train_data[with(train_data, order(VisitLink)), ]
  test_data = test_data[with(test_data, order(VisitLink)), ]
  mod = parfm(Surv(DtoA,status) ~ 
                LOS
              +AGE
              +ASOURCE2
              +ASOURCE3
              +ASOURCE4
              +ASOURCE5
              +ASCHED1
              +DISPUNIFORM2
              +DISPUNIFORM5
              +DISPUNIFORM6
              +DISPUNIFORM7
              +DISPUNIFORM20
              +DISPUNIFORM99
              +BMONTH
              +FEMALE1
              +HISPANIC_X2
              +HISPANIC_X3
              +HISPANIC_X4
              +RACE2
              +RACE3
              +RACE4
              +RACE5
              +RACE6
              +RACE7
              +AWEEKEND1
              +AMONTH
              +DMONTH
              +DNR1
              +DQTR
              +NCHRONIC
              +NDX
              +NECODE
              +NPR
              +TVisits
              +NVisit, cluster = "VisitLink",
              data = train_data, dist = "exponential", frailty = "gamma")
  base = list()
  base[[1]] = base_exp(train_data,mod)
  base[[2]] = base_exp(test_data,mod)
  train_cluster = train_data$VisitLink
  test_cluster = test_data$VisitLink
  risk = parfmrisk(mod,train_data,test_data,train_cluster,test_cluster,base,numbase = 1)
  cindex = list()
  cindex[[1]] = GetCindex.uncensored(risk[[1]], train_data$DtoA)
  cindex[[2]] = GetCindex.uncensored(risk[[2]], test_data$DtoA)
  cindex[[3]] = risk[[3]]
  cindex[[4]] = risk[[4]]
  cindex[[5]] = mod
  return(cindex)
}


parfm_weibull_gamma_patient = function(train_data,test_data){
  train_data = train_data[with(train_data, order(VisitLink)), ]
  test_data = test_data[with(test_data, order(VisitLink)), ]
  mod = parfm(Surv(DtoA,status) ~ 
                LOS
              +AGE
              +ASOURCE2
              +ASOURCE3
              +ASOURCE4
              +ASOURCE5
              +ASCHED1
              +DISPUNIFORM2
              +DISPUNIFORM5
              +DISPUNIFORM6
              +DISPUNIFORM7
              +DISPUNIFORM20
              +DISPUNIFORM99
              +BMONTH
              +FEMALE1
              +HISPANIC_X2
              +HISPANIC_X3
              +HISPANIC_X4
              +RACE2
              +RACE3
              +RACE4
              +RACE5
              +RACE6
              +RACE7
              +AWEEKEND1
              +AMONTH
              +DMONTH
              +DNR1
              +DQTR
              +NCHRONIC
              +NDX
              +NECODE
              +NPR
              +TVisits
              +NVisit, cluster = "VisitLink",
              data = train_data, dist = "weibull", frailty = "gamma")
  base = list()
  base[[1]] = base_weibull(train_data,mod)
  base[[2]] = base_weibull(test_data,mod)
  train_cluster = train_data$VisitLink
  test_cluster = test_data$VisitLink
  risk = parfmrisk(mod,train_data,test_data,train_cluster,test_cluster,base,numbase = 2)
  cindex = list()
  cindex[[1]] = GetCindex.uncensored(risk[[1]], train_data$DtoA)
  cindex[[2]] = GetCindex.uncensored(risk[[2]], test_data$DtoA)
  cindex[[3]] = risk[[3]]
  cindex[[4]] = risk[[4]]
  cindex[[5]] = mod
  return(cindex)
}


parfm_lognormal_gamma_patient = function(train_data,test_data){
  train_data = train_data[with(train_data, order(VisitLink)), ]
  test_data = test_data[with(test_data, order(VisitLink)), ]
  mod = parfm(Surv(DtoA,status) ~ 
                LOS
              +AGE
              +ASOURCE2
              +ASOURCE3
              +ASOURCE4
              +ASOURCE5
              +ASCHED1
              +DISPUNIFORM2
              +DISPUNIFORM5
              +DISPUNIFORM6
              +DISPUNIFORM7
              +DISPUNIFORM20
              +DISPUNIFORM99
              +BMONTH
              +FEMALE1
              +HISPANIC_X2
              +HISPANIC_X3
              +HISPANIC_X4
              +RACE2
              +RACE3
              +RACE4
              +RACE5
              +RACE6
              +RACE7
              +AWEEKEND1
              +AMONTH
              +DMONTH
              +DNR1
              +DQTR
              +NCHRONIC
              +NDX
              +NECODE
              +NPR
              +TVisits
              +NVisit, cluster = "VisitLink",
              data = train_data, dist = "lognormal", frailty = "gamma")
  base = list()
  base[[1]] = base_lognormal(train_data,mod)
  base[[2]] = base_lognormal(test_data,mod)
  train_cluster = train_data$VisitLink
  test_cluster = test_data$VisitLink
  risk = parfmrisk(mod,train_data,test_data,train_cluster,test_cluster,base,numbase = 2)
  cindex = list()
  cindex[[1]] = GetCindex.uncensored(risk[[1]], train_data$DtoA)
  cindex[[2]] = GetCindex.uncensored(risk[[2]], test_data$DtoA)
  cindex[[3]] = risk[[3]]
  cindex[[4]] = risk[[4]]
  cindex[[5]] = mod
  return(cindex)
}

parfm_loglogistic_gamma_patient = function(train_data,test_data){
  train_data = train_data[with(train_data, order(VisitLink)), ]
  test_data = test_data[with(test_data, order(VisitLink)), ]
  mod = parfm(Surv(DtoA,status) ~ 
                LOS
              +AGE
              +ASOURCE2
              +ASOURCE3
              +ASOURCE4
              +ASOURCE5
              +ASCHED1
              +DISPUNIFORM2
              +DISPUNIFORM5
              +DISPUNIFORM6
              +DISPUNIFORM7
              +DISPUNIFORM20
              +DISPUNIFORM99
              +BMONTH
              +FEMALE1
              +HISPANIC_X2
              +HISPANIC_X3
              +HISPANIC_X4
              +RACE2
              +RACE3
              +RACE4
              +RACE5
              +RACE6
              +RACE7
              +AWEEKEND1
              +AMONTH
              +DMONTH
              +DNR1
              +DQTR
              +NCHRONIC
              +NDX
              +NECODE
              +NPR
              +TVisits
              +NVisit, cluster = "VisitLink",
              data = train_data, dist = "loglogistic", frailty = "gamma")
  base = list()
  base[[1]] = base_loglogistic(train_data,mod)
  base[[2]] = base_loglogistic(test_data,mod)
  train_cluster = train_data$VisitLink
  test_cluster = test_data$VisitLink
  risk = parfmrisk(mod,train_data,test_data,train_cluster,test_cluster,base,numbase = 2)
  cindex = list()
  cindex[[1]] = GetCindex.uncensored(risk[[1]], train_data$DtoA)
  cindex[[2]] = GetCindex.uncensored(risk[[2]], test_data$DtoA)
  cindex[[3]] = risk[[3]]
  cindex[[4]] = risk[[4]]
  cindex[[5]] = mod
  return(cindex)
}

#########################################Region######################################
parfm_exp_gamma_region = function(train_data,test_data){
  train_data = train_data[with(train_data, order(PL_NCHS2006)), ]
  test_data = test_data[with(test_data, order(PL_NCHS2006)), ]
  mod = parfm(Surv(DtoA,status) ~ 
                LOS
              +AGE
              +ASOURCE2
              +ASOURCE3
              +ASOURCE4
              +ASOURCE5
              +ASCHED1
              +DISPUNIFORM2
              +DISPUNIFORM5
              +DISPUNIFORM6
              +DISPUNIFORM7
              +DISPUNIFORM20
              +DISPUNIFORM99
              +BMONTH
              +FEMALE1
              +HISPANIC_X2
              +HISPANIC_X3
              +HISPANIC_X4
              +RACE2
              +RACE3
              +RACE4
              +RACE5
              +RACE6
              +RACE7
              +AWEEKEND1
              +AMONTH
              +DMONTH
              +DNR1
              +DQTR
              +NCHRONIC
              +NDX
              +NECODE
              +NPR
              +TVisits
              +NVisit, cluster = "PL_NCHS2006",
              data = train_data, dist = "exponential", frailty = "gamma")
  base = list()
  base[[1]] = base_exp(train_data,mod)
  base[[2]] = base_exp(test_data,mod)
  train_cluster = train_data$PL_NCHS2006
  test_cluster = test_data$PL_NCHS2006
  risk = parfmrisk(mod,train_data,test_data,train_cluster,test_cluster,base,numbase = 1)
  cindex = list()
  cindex[[1]] = GetCindex.uncensored(risk[[1]], train_data$DtoA)
  cindex[[2]] = GetCindex.uncensored(risk[[2]], test_data$DtoA)
  cindex[[3]] = risk[[3]]
  cindex[[4]] = risk[[4]]
  cindex[[5]] = mod
  return(cindex)
}

parfm_weibull_gamma_region = function(train_data,test_data){
  train_data = train_data[with(train_data, order(PL_NCHS2006)), ]
  test_data = test_data[with(test_data, order(PL_NCHS2006)), ]
  mod = parfm(Surv(DtoA,status) ~ 
                LOS
              +AGE
              +ASOURCE2
              +ASOURCE3
              +ASOURCE4
              +ASOURCE5
              +ASCHED1
              +DISPUNIFORM2
              +DISPUNIFORM5
              +DISPUNIFORM6
              +DISPUNIFORM7
              +DISPUNIFORM20
              +DISPUNIFORM99
              +BMONTH
              +FEMALE1
              +HISPANIC_X2
              +HISPANIC_X3
              +HISPANIC_X4
              +RACE2
              +RACE3
              +RACE4
              +RACE5
              +RACE6
              +RACE7
              +AWEEKEND1
              +AMONTH
              +DMONTH
              +DNR1
              +DQTR
              +NCHRONIC
              +NDX
              +NECODE
              +NPR
              +TVisits
              +NVisit, cluster = "PL_NCHS2006",
              data = train_data, dist = "weibull", frailty = "gamma")
  base = list()
  base[[1]] = base_weibull(train_data,mod)
  base[[2]] = base_weibull(test_data,mod)
  train_cluster = train_data$PL_NCHS2006
  test_cluster = test_data$PL_NCHS2006
  risk = parfmrisk(mod,train_data,test_data,train_cluster,test_cluster,base,numbase = 2)
  cindex = list()
  cindex[[1]] = GetCindex.uncensored(risk[[1]], train_data$DtoA)
  cindex[[2]] = GetCindex.uncensored(risk[[2]], test_data$DtoA)
  cindex[[3]] = risk[[3]]
  cindex[[4]] = risk[[4]]
  cindex[[5]] = mod
  return(cindex)
}




parfm_lognormal_gamma_region = function(train_data,test_data){
  train_data = train_data[with(train_data, order(PL_NCHS2006)), ]
  test_data = test_data[with(test_data, order(PL_NCHS2006)), ]
  mod = parfm(Surv(DtoA,status) ~ 
                LOS
              +AGE
              +ASOURCE2
              +ASOURCE3
              +ASOURCE4
              +ASOURCE5
              +ASCHED1
              +DISPUNIFORM2
              +DISPUNIFORM5
              +DISPUNIFORM6
              +DISPUNIFORM7
              +DISPUNIFORM20
              +DISPUNIFORM99
              +BMONTH
              +FEMALE1
              +HISPANIC_X2
              +HISPANIC_X3
              +HISPANIC_X4
              +RACE2
              +RACE3
              +RACE4
              +RACE5
              +RACE6
              +RACE7
              +AWEEKEND1
              +AMONTH
              +DMONTH
              +DNR1
              +DQTR
              +NCHRONIC
              +NDX
              +NECODE
              +NPR
              +TVisits
              +NVisit, cluster = "PL_NCHS2006",
              data = train_data, dist = "lognormal", frailty = "gamma")
  base = list()
  base[[1]] = base_lognormal(train_data,mod)
  base[[2]] = base_lognormal(test_data,mod)
  train_cluster = train_data$PL_NCHS2006
  test_cluster = test_data$PL_NCHS2006
  risk = parfmrisk(mod,train_data,test_data,train_cluster,test_cluster,base,numbase = 2)
  cindex = list()
  cindex[[1]] = GetCindex.uncensored(risk[[1]], train_data$DtoA)
  cindex[[2]] = GetCindex.uncensored(risk[[2]], test_data$DtoA)
  cindex[[3]] = risk[[3]]
  cindex[[4]] = risk[[4]]
  cindex[[5]] = mod
  return(cindex)
}


parfm_loglogistic_gamma_region = function(train_data,test_data){
  train_data = train_data[with(train_data, order(PL_NCHS2006)), ]
  test_data = test_data[with(test_data, order(PL_NCHS2006)), ]
  mod = parfm(Surv(DtoA,status) ~ 
                LOS
              +AGE
              +ASOURCE2
              +ASOURCE3
              +ASOURCE4
              +ASOURCE5
              +ASCHED1
              +DISPUNIFORM2
              +DISPUNIFORM5
              +DISPUNIFORM6
              +DISPUNIFORM7
              +DISPUNIFORM20
              +DISPUNIFORM99
              +BMONTH
              +FEMALE1
              +HISPANIC_X2
              +HISPANIC_X3
              +HISPANIC_X4
              +RACE2
              +RACE3
              +RACE4
              +RACE5
              +RACE6
              +RACE7
              +AWEEKEND1
              +AMONTH
              +DMONTH
              +DNR1
              +DQTR
              +NCHRONIC
              +NDX
              +NECODE
              +NPR
              +TVisits
              +NVisit, cluster = "PL_NCHS2006",
              data = train_data, dist = "loglogistic", frailty = "gamma")
  base = list()
  base[[1]] = base_loglogistic(train_data,mod)
  base[[2]] = base_loglogistic(test_data,mod)
  train_cluster = train_data$PL_NCHS2006
  test_cluster = test_data$PL_NCHS2006
  risk = parfmrisk(mod,train_data,test_data,train_cluster,test_cluster,base,numbase = 2)
  cindex = list()
  cindex[[1]] = GetCindex.uncensored(risk[[1]], train_data$DtoA)
  cindex[[2]] = GetCindex.uncensored(risk[[2]], test_data$DtoA)
  cindex[[3]] = risk[[3]]
  cindex[[4]] = risk[[4]]
  cindex[[5]] = mod
  return(cindex)
}

#######################################################################
cumhaz_exp_t = function(t,mod){
  mod[2,1]*t
}
cumhaz_weibull_t = function(t,mod){
  mod[3,1]*(t^(mod[2,1]))
}

cumhaz_lognormal_t = function(t,mod){
  z = (log(t)-mod[2,1])/mod[3,1]
  cum = -log(1-pnorm(z))
  return(cum)
}
cumhaz_loglogistic_t = function(t,mod){
  den = 1+exp(mod[2,1])*(t^mod[3,1])
  cum = log(den)
  return(cum)
}

cumhaz_exp = function(data,mod){
  mod[2,1]*data$DtoA
}
cumhaz_weibull = function(data,mod){
  mod[3,1]*(data$DtoA^(mod[2,1]))
}
cumhaz_gompertz = function(data,mod){
  (mod[3,1]/mod[2,1])*exp(mod[2,1]*data$DtoA-1)
}
cumhaz_lognormal = function(data,mod){
  z = (log(data$DtoA)-mod[2,1])/mod[3,1]
  cum = -log(1-pnorm(z))
  return(cum)
}
cumhaz_loglogistic = function(data,mod){
  den = 1+exp(mod[2,1])*(data$DtoA^mod[3,1])
  cum = log(den)
  return(cum)
}

