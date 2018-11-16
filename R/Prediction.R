# Copyright 2018 Observational Health Data Sciences and Informatics
#
# This file is part of DistributedRegressionEval
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


expit = function(x){exp(x)/(1+exp(x))}

#likelihood function for logistic regression, the input X is a n*d matrix where
#each patient has d covariates stored in each row.
Lik = function(beta,X,Y){
  design = cbind(1,X)
  sum(Y*(design%*%t(t(beta)))-log(1+exp(design%*%t(t(beta)))))/length(Y)
}

#Y is the binary vector with disease status per patient (1 indicates case and 0 indicates control)

#first order gradient
Lgradient = function(beta,X,Y){
  design = cbind(1,X)
  t(Y-expit(design%*%t(t(beta))))%*%design/length(Y)
}

#first-order surogate likelihood, suppose the local data are stored in Xlocal,Ylocal
SL = function(beta, Xlocal,Ylocal, L){
  -Lik(beta,Xlocal,Ylocal) - L%*%beta
}

#second-order gradient
Lgradient2 =function(beta,X){
  design = cbind(1,X)
  Z=expit(design%*%beta)
  t(c(-1*Z*(1-Z))*design)%*%design/nrow(X)

}

#second-order surogate likelihood
SL2 = function(beta, beta0, Xlocal, Ylocal, L, L2){
  beta_temp = beta - beta0
  -Lik(beta,Xlocal,Ylocal) - L%*%beta - 0.5*t(beta_temp)%*%L2%*%t(t(beta_temp))
}


##Meat of the sandwich estimator of variance ####
Lgradient_meat = function(beta,X,Y){
  design = cbind(1,X)
  Z = expit(design%*%beta)
  t(c(Y-Z)*design)%*%(c(Y-Z)*design)/nrow(X)
}

#Sandwich estimator for variance
#N is the total sample size###
Sandwich = function(beta,X,Y,N){
  mat_L1 = Lgradient_meat(beta,X,Y)
  mat_L2 = Lgradient2(beta,X)
  inv_L2 = solve.default(mat_L2)
  out = inv_L2%*%mat_L1%*%inv_L2/N
  return(out)
}


evaluateOdal <- function(studyFolder, outcomeId) {
  # outcomeId <- 4
  data <- readRDS(file.path(studyFolder, sprintf("data_o%s.rds", outcomeId)))
  localDb <- "ccae"
  otherDbs <- unique(data$database[data$database != localDb])
  Ylocal <- data$y[data$database == localDb]
  Xlocal <- as.matrix(data[data$database == localDb, -which(colnames(data) %in% c("y", "database"))])
  Yall <- data$y
  Xall <- as.matrix(data[, -which(colnames(data) %in% c("y", "database"))])
  K <- length(otherDbs)
  p <- ncol(Xlocal) + 1
  nlocal <- length(Ylocal)

  ################################################
  # PART 3.. Initialization     #
  ###############################################

  #Local estimator is obtained for initialization
  fit0 = summary(glm(Ylocal~Xlocal, family = "binomial"(link = "logit")))
  beta0 = fit0$coefficients[,1]
  #Then beta0 is passed to each site
  v_beta0 = fit0$cov.scaled#covariance matrix of beta0


  ##########################################
  # PART 4.. Calculation in each site     #
  #########################################
  #calculate the gradient in each site
  L = matrix(0,nrow = K,ncol = p) # Store the first order gradient (each is a p-dimensional vector)
  L2 = matrix(0,nrow = K,ncol = p^2)# Store the second order gradient (each is a p*p matrix, we expand it into a vector)


  nsite = rep(0,K)  #sample size in each site
  for (i in 1:K) {
    Ysite <- data$y[data$database == otherDbs[i]]
    Xsite <- as.matrix(data[data$database == otherDbs[i], -which(colnames(data) %in% c("y", "database"))])
    nsite[i] = length(Ysite)
    L[i,] = Lgradient(beta0,Xsite,Ysite)
    L2[i,] = as.vector(Lgradient2(beta0,Xsite))
  }
  #Each L[i,] is transfered to the local site for ODAL1
  # L[i,] and L2[i,] are transfered to the local site for ODAL2


  ##########################################
  # PART 5.. Local process    #
  #########################################


  #Calculate the global first order gradient L
  L = apply(diag(nsite)%*%L,2,sum)/(sum(nsite)+nlocal)

  #Calculate the global first order gradient L
  L2_local = Lgradient2(beta0,Xlocal)
  #Calculate the global second order gradient L2
  L2 = apply(diag(c(nlocal,nsite))%*%rbind(as.vector(L2_local),L2),2,sum)/(sum(nsite)+nlocal)
  L2 = matrix(L2,ncol = p, nrow = p)-L2_local

  ######################################
  # PART 6.. Estimation and Inference  #
  #####################################
  #Point estimation
  ODAL1 = optim(beta0,SL,control = list(maxit = 10000,reltol = 1e-16), Xlocal = Xlocal, Ylocal = Ylocal, L = L)$par
  ODAL2 = optim(beta0,SL2,control = list(maxit = 10000,reltol = 1e-16), Xlocal = Xlocal, Ylocal = Ylocal, L = L, L2 = L2, beta0 = beta0)$par

  ####Variance####
  var1 = tryCatch(Sandwich(ODAL1,Xlocal,Ylocal,(sum(nsite)+nlocal)),error=function(err) matrix(data=NA,nrow=length(ODAL1),ncol=length(ODAL1)))
  var2 = Sandwich(ODAL2,Xlocal,Ylocal,(sum(nsite)+nlocal))
  var_ODAL1 = diag(var1)
  var_ODAL2 = diag(var2)


  ######################################
  # PART 7.. Pooled estimator  #
  #####################################
  fitall = summary(glm(Yall~Xall, family = "binomial"(link = "logit")))
  betaall = fitall$coefficients[,1]
  v_betaall = fitall$cov.scaled



  ######################################
  # PART 8.. Meta estimator  #
  #####################################
  Beta = matrix(0,nrow = K,ncol = p) # Store the estimator in each site
  VBeta = matrix(0,nrow = K,ncol = p)# Store the covariance matrix from each site (each is a p*p matrix, we expand it into a vector)

  #fit logistic regression in each site
  for (i in 1:K){
    Ysite <- data$y[data$database == otherDbs[i]]
    Xsite <- as.matrix(data[data$database == otherDbs[i], -which(colnames(data) %in% c("y", "database"))])
    Beta[i,] = tryCatch(glm(Ysite~Xsite, family = "binomial"(link = "logit"))$coefficients[,1],error=function(err) rep(NA,length(betaall)))
    if(sum(is.na(Beta[i,]))!=0){
      Beta[i,] = rep(NA,length(betaall))
      VBeta[i,] = rep(NA,length(betaall))
    }
    else{ VBeta[i,] = summary(glm(Ysite~Xsite, family = "binomial"(link = "logit")))$coefficients[,2]^2}
  }

  Beta = rbind(Beta,beta0)
  VBeta = rbind(VBeta,diag(v_beta0))

  #estimate from meta-analysis
  betameta = apply(Beta/VBeta,2,function(x){sum(x, na.rm = T)})/apply(1/VBeta,2,function(x){sum(x, na.rm = T)})
  vmeta = 1/apply(1/VBeta,2,function(x){sum(x, na.rm = T)})


  ######################################
  # PART 9.. output  #
  #####################################
  output = list(beta0,v_beta0,ODAL1,var1,ODAL2,var2,betaall,v_betaall,betameta,vmeta)

  pathToCsv <- system.file("settings", "CohortsToCreate.csv", package = "DistributedRegressionEval")
  cohortsToCreate <- read.csv(pathToCsv)
  outcomeName <- cohortsToCreate$name[cohortsToCreate$cohortId == outcomeId]

  saveRDS(output, file.path(studyFolder, sprintf("output_%s.rds", outcomeName)))

  pretty <- rbind(beta0, ODAL1, ODAL2, betaall, betameta)
  write.csv(pretty, file.path(studyFolder, sprintf("output_%s.csv", outcomeName)))
}
