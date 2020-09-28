#rm(list=ls())
library(tidyverse)
library(gridExtra)
library(Rsolnp)
library(data.table)
library(effsize)


#----------------------------------------------------------#
# Q learning
#----------------------------------------------------------#

func_qlearning = function(Model, param, Est, c, r, Sim, r1, r2){
  
  T <- length(c)
  
  beta <- param[1]
  alpha <- param[2]
  
  if(grepl("F",Model) == 1){
    phi <- param[3]
    Npara <- 3
    mu <- 0.5
  } else {
    phi <- 0
    Npara <-2
    mu <- 0.5
  }
  
  
  p1 <- numeric(T) # choice probability of option1 
  p1[1] <- 0.5 # initial choice probability
  
  RPE <- numeric(T)
  cQ <- numeric(T) # Q of the chosen choice in each trial
  cQ[1] <- mu
  Q <- matrix(numeric(2*T), nrow=2)
  Q[,1] <- mu
  
  # initialize log-likelihood
  
  ll <- 0
  
  for (t in 1:T ){
    
    p1[t] <- 1/(1+exp(-beta * (Q[1,t] - Q[2,t]))) 
    
    #--- USE for simulation ---#
    if(Sim == 1) {
      if(p1[t] > runif(1)){ 
        c[t] <- 1
        r[t] <- r1[t]
      } else {
        c[t] <- 2
        r[t] <- r2[t]
      }
    }
    #--------------------------#
    
    ll <- ll + (c[t]==1) * log(p1[t]) +  (c[t]==2) * log(1-p1[t]) 
    
    
    ##---   update value   ---##
    
    cQ[t] <- Q[c[t],t]
    RPE[t] <- r[t] - Q[c[t],t]
    if (t < T) {
      Q[c[t],t+1] <- Q[c[t],t] + alpha*RPE[t] 
      
      Q[3-c[t],t+1] <- Q[3-c[t],t] + phi*(mu - Q[3-c[t],t])
      
    }
  }
  
  return(list(negll = -ll, r = r, c = c, p1 = p1, Q = Q, cQ = cQ, RPE = RPE, Npara = Npara))
  
}

#------------------------------------#
# Data Generation
#------------------------------------#

func_datageneration = function(Model, trueparam, RHtype, N){
  
  Sim = 1
  DF <- data.frame()
  
  for (i in 1:N){
    
    r <- numeric(nTrial)
    c <- numeric(nTrial)
    
    if(is.numeric(RHtype)){
      
      r1 <- as.numeric()
      r2 <- as.numeric()
      m1 <- as.numeric()
      m2 <- as.numeric()
      for(ii in 0:RHtype){
        if(ii%%2 == 0){
          r1 <- c(r1,sample(rbetter))
          r2 <- c(r2,sample(rworse))
          m1 <- c(m1,rep(0.8,length(rbetter)))
          m2 <- c(m2,rep(0.2,length(rworse)))
        } else {
          r1 <- c(r1,sample(rworse))
          r2 <- c(r2,sample(rbetter))
          m1 <- c(m1,rep(0.2,length(rworse)))
          m2 <- c(m2,rep(0.8,length(rbetter)))
        }
      }
      
    }
    
    result <- func_qlearning(Model, trueparam, Est, c, r, Sim, r1, r2)
    
    # fMRI signal
    y <- result$RPE + rnorm(nTrial, mean = 0, sd = 0.5) 
    
    temp <- data.frame(RHtype = RHtype,
                       id = i,
                       trial = seq(1,length(result$c),1),
                       r1 = r1,
                       r2 = r2,
                       m1 = m1,
                       m2 = m2,
                       Q1 = result$Q[1,], 
                       Q2 = result$Q[2,],
                       cQ = result$cQ,
                       RPE = result$RPE,
                       y = y, # as brain activity
                       c = result$c, 
                       r = result$r, 
                       p1 =result$p1)
    DF <- rbind(DF,temp)
  }
  DF <<- DF
}

#------------------------------------#
# Model Fitting
#------------------------------------#


#----- func_minimize -----#

func_minimize <- function(Model, param, Est, c, r) {
  Sim=0
  result_fit<-NULL
  ret <- func_qlearning(Model, param, Est, c, r, Sim, NULL, NULL)
  result_fit <<- ret
  return(ret$negll) 
} 

#----- ParameterFit -----# 

func_ParameterFit = function(Est, Model, DF){
  
  DF_param <- data.frame()
  DF_fit <- data.frame()
  
  for (CurrentSub in 1:max(DF$id)) {
    
    if(length(DF$id[DF$id == CurrentSub])==0){ next }
    
    cdata <- subset(DF,id == CurrentSub)
    
    fvalmin = Inf;
    
    for (idx in 1:3) {
      
      names <- c("beta","alpha")
      
      if(grepl("F",Model) == 1){
        names <- c(names,"alphaF")
      }
      
      initparam <- runif(length(names), 0, 1) 
      lblist <- rep(0,length(names))
      ublist <- rep(1,length(names))
      ublist[1] <- 20
      
      ##------------------------------------##
      ##      Optimization with solnp()     ##
      ##------------------------------------##
      sol<-solnp(pars=initparam, fun=func_minimize, LB=lblist, UB=ublist,
                 control=list(trace=0),
                 c = cdata$c,
                 r = cdata$r,
                 Model=Model,
                 Est=Est)
      
      nll<-sol$values[length(sol$values)]
      print(nll)# show negative LL
      print(sol$pars)  # show parameters
      
      if(nll < fvalmin){# when the LL is minimum then it will be recorded
        fvalmin<-nll
        paramest<-sol$pars
        sol_best<-sol
      }
    }
    
    ##--- For regression analyses ---##
    temp_fit <- cbind(id = CurrentSub,
                      trial = seq(1,length(result_fit$c),1),
                      Q1_fit = result_fit$Q[1,], 
                      Q2_fit = result_fit$Q[2,],
                      cQ_fit = result_fit$cQ,
                      RPE_fit = result_fit$RPE)
    DF_fit <- rbind(DF_fit,temp_fit) 
    ##-------------------------------##
    
    aic <- 2*fvalmin + 2*length(names) 
    bic <- 2*fvalmin + length(names)*log(nrow(cdata))
    negll <- fvalmin 
    
    df1 <- data.frame(sub = CurrentSub, negll = negll, aic =aic, bic = bic,
                      param1 = paramest[1],
                      param2 = paramest[2],
                      param3 = paramest[3],
                      param4 = paramest[4],
                      param5 = paramest[5])
    colnames(df1)<-c("sub","negll","aic","bic",names)
    print(round(df1,2))
    
    ##---  estimated parameters   ---##
    temp_param <- cbind(
      Est = Est,
      Model = Model,
      id = CurrentSub,
      negLL = negll,
      beta = paramest[1],
      alphaL = paramest[2]
    )
    if(grepl("F",Model) == 1){temp_param <- cbind(temp_param, alphaF = paramest[3]) }
    
    DF_param <- rbind(DF_param,temp_param) 
    ##-------------------------------##
  }
  
  DF_fit <<- DF_fit
  DF_param <<- DF_param
  
}

#------------------------------------#
# Regression
#------------------------------------#
# GLM1: simple regression on fit RPE only (beta.delta0)
# GLM2: multiple regression on fit RPE and reward (beta.delta, beta.r)

func_regression <- function(D){
  dfstat <- data.table()
  for(i in 1:max(D$id)){
    dftmp <- D %>% filter(id == i)
    
    dftmp <- dftmp %>% mutate(y=y-mean(y),
                              RPE=RPE-mean(RPE),
                              RPE_fit=RPE_fit-mean(RPE_fit),
                              nQ_fit = -cQ_fit + mean(cQ_fit)) 
    
    fit.GLM1 <- glm(y ~ RPE_fit, data = dftmp)
    fit.GLM2 <- glm(y ~ r + nQ_fit, data = dftmp)
    fit.GLM2p <- glm(y ~ r + RPE_fit, data = dftmp)
    
    dftmp2 <- data.table(group = unique(dftmp$group), 
                         subject = i,
                         beta.GLM1.delta = fit.GLM1$coefficients["RPE_fit"],
                         beta.GLM2.r = fit.GLM2$coefficients["r"],
                         beta.GLM2.nQ = fit.GLM2$coefficients["nQ_fit"],
                         beta.GLM2p.r = fit.GLM2p$coefficients["r"],
                         beta.GLM2p.delta = fit.GLM2p$coefficients["RPE_fit"],
                         AIC.GLM1 = AIC(fit.GLM1),
                         AIC.GLM2 = AIC(fit.GLM2),
                         AIC.GLM2p = AIC(fit.GLM2p)
    )
    dfstat <- dplyr::bind_rows(dfstat, dftmp2)
  }
  dfstat <<- dfstat
}

#------------------------------------#
# t test
#------------------------------------#

ttest <- function(x1,x2,pair){
  t <- t.test(x1,x2, paired=pair,var.equal = T)
  tvalue <- round(t$statistic,2)
  df <- t$parameter
  pvalue <- round(t$p.value,4)
  cohenD <- cohen.d(x1,x2, paired=pair)$estimate
  #cat(paste("t(",df,") = ",tvalue, " (p = ", pvalue, "; d = ", cohenD,")", sep=""))
  return(list(df=df,tvalue=tvalue,pvalue=pvalue,cohenD=cohenD))
}

#------------------------------------#
# FUNCTIONS to change Reward History
#------------------------------------#
func_ParameterFit_Reverse <- function(Est,Data,Model){
  
  DF_param_all <- data.frame()
  DF_fit_all <- data.frame()
  
  for (i in 1:7){
    tempRHtype <- c(0,1,2,5,8,17,35)[i]
    func_ParameterFit(Est, Model, Data %>% filter(RHtype == tempRHtype))
    
    temp1 <- DF_param %>%  mutate(RHtype = tempRHtype)
    DF_param_all <- rbind(DF_param_all,temp1)
    
    temp2 <- DF_fit %>%  mutate(RHtype = tempRHtype)
    DF_fit_all <- rbind(DF_fit_all,temp2)
  }
  DF_param_all <<- DF_param_all
  DF_fit_all <<- DF_fit_all
}

#------------------------------------#
## Regression for Reverse conditions
#------------------------------------#

func_Regression_Reverse <- function(Data){
  D <- Data %>% filter(group=="FHnoF")
  func_regression(D)
  Reg_FHnoF <- tbl_df(dfstat)
  
  D <- Data %>% filter(group=="FLnoF")
  func_regression(D)
  Reg_FLnoF <- tbl_df(dfstat)
  
  D <- Data %>% filter(group=="FHF")
  func_regression(D)
  Reg_FHF <- tbl_df(dfstat)
  
  D <- Data %>% filter(group=="FLF")
  func_regression(D)
  Reg_FLF <- tbl_df(dfstat)
  
  temp_Reg_All <<- dplyr::bind_rows(Reg_FHnoF,Reg_FLnoF,Reg_FHF,Reg_FLF)
}

#------------------------------------#
## Effect size
#------------------------------------#

func_cohenD <- function(Reg_All){
  
  #--------
  # (A) Fit without forgetting 
  #--------
  
  # GLM1 delta 
  y <- Reg_All %>% filter(group == "FHnoF") %>% select(beta.GLM1.delta)
  x <- Reg_All %>% filter(group == "FLnoF") %>% select(beta.GLM1.delta)
  ttest(x[[1]],y[[1]],0) -> ttest.noF_GLM1d
  
  
  # GLM2 r
  y <- Reg_All %>% filter(group == "FHnoF") %>% select(beta.GLM2.r)
  x <- Reg_All %>% filter(group == "FLnoF") %>% select(beta.GLM2.r)
  ttest(x[[1]],y[[1]],0) -> ttest.noF_GLM2r
  
  # GLM2 nV 
  y <- Reg_All %>% filter(group == "FHnoF") %>% select(beta.GLM2.nQ)
  x <- Reg_All %>% filter(group == "FLnoF") %>% select(beta.GLM2.nQ)
  ttest(x[[1]],y[[1]],0) -> ttest.noF_GLM2nV
  
  #--------
  # (B) Fit with forgetting 
  #--------
  
  # GLM1 delta 
  y <- Reg_All %>% filter(group == "FHF") %>% select(beta.GLM1.delta)
  x <- Reg_All %>% filter(group == "FLF") %>% select(beta.GLM1.delta)
  ttest(x[[1]],y[[1]],0) -> ttest.F_GLM1d
  
  # GLM2 r
  y <- Reg_All %>% filter(group == "FHF") %>% select(beta.GLM2.r)
  x <- Reg_All %>% filter(group == "FLF") %>% select(beta.GLM2.r)
  ttest(x[[1]],y[[1]],0) -> ttest.F_GLM2r
  
  # GLM2 nV
  y <- Reg_All %>% filter(group == "FHF") %>% select(beta.GLM2.nQ)
  x <- Reg_All %>% filter(group == "FLF") %>% select(beta.GLM2.nQ)
  ttest(x[[1]],y[[1]],0) -> ttest.F_GLM2nV
  
  temp_cohenD_All <<- data.frame(
    noF_GLM1d = ttest.noF_GLM1d$cohenD,
    noF_GLM2r = ttest.noF_GLM2r$cohenD,
    noF_GLM2nV = ttest.noF_GLM2nV$cohenD,
    F_GLM1d = ttest.F_GLM1d$cohenD,
    F_GLM2r = ttest.F_GLM2r$cohenD,
    F_GLM2nV = ttest.F_GLM2nV$cohenD,
    p.noF_GLM1d = ttest.noF_GLM1d$pvalue,
    p.noF_GLM2r = ttest.noF_GLM2r$pvalue,
    p.noF_GLM2nV = ttest.noF_GLM2nV$pvalue,
    p.F_GLM1d = ttest.F_GLM1d$pvalue,
    p.F_GLM2r = ttest.F_GLM2r$pvalue,
    p.F_GLM2nV = ttest.F_GLM2nV$pvalue) 
}