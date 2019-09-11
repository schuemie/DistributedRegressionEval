##################################################
###Distributed algorithm for logistic regression##
##################################################
###############Updated on 07/16/2019##############

expit = function(x){exp(x)/(1+exp(x))}
library(MASS)
###########################################################################################
#The function "compare.methods()" returns the point estimation and hessian matrix for the following methods
#1. Logistic regression using only local data
#2. Logistic regression using all data from multiple sites (pooled together).
#3. Distributed algorithm ODAL1 **mean** using local estiator (from #1) as initial value.
#4. Distributed algorithm ODAL1 **median** using local estiator (from #1) as initial value.
#5. Distributed algorithm ODAL2 **mean** using local estiator (from #1) as initial value.
#6. Distributed algorithm ODAL2 **median** using local estiator (from #1) as initial value.
#7. Fix effect Meta-analysis
#8. Distributed algorithm ODAL1 **mean** using meta estiator (from #5) as initial value.
#9. Distributed algorithm ODAL1 **median** using meta estiator (from #5) as initial value.
#10. Distributed algorithm ODAL2 **mean** using meta estiator (from #5) as initial value.
#11. Distributed algorithm ODAL2 **median** using meta estiator (from #5) as initial value.
###########################################################################################



###########################################################################################
#The input of the function is
#1. Xall: p predictors for all N patients as a N*p matrix (N is the total sample size of all sites)
#2. Yall: the corresponding N binary outcomes as a vector.
#3. site: Index of the clinical sites. Local site should be labeled as 0, and other K sites should be labeled from 1 to K.


compare.methods = function(Xall, Yall, site){

  Xlocal = Xall[which(site==0),] # extract the local site X, which indicator "site = 0"
  Ylocal = Yall[which(site==0)] # extract the local site Y, which indicator "site = 0"
  K = length(unique(site))-1 # K (number of sites excludes local) = (the number of sites (includes local) - 1)
  p = dim(Xall)[2]+1 # number of variables include the intercept
  nlocal = length(Ylocal) # number of patients in local site
  ################################################
  #      Functions              #
  ###############################################

  #likelihood function for logistic regression, the input X is a n*d matrix where
  #each patient has d covariates stored in each row.
  Lik = function(beta,X,Y){
    design = cbind(1,X)
    sum(Y*(design%*%t(t(beta)))-log(1+exp(design%*%t(t(beta)))))/length(Y)
  }

  #Y is the binary vector with disease status per patient (1 indicates case and 0 indicates control)

  ####### first order gradient ######
  Lgradient = function(beta,X,Y){
    design = cbind(1,X)
    t(Y-expit(design%*%t(t(beta))))%*%design/length(Y)
  }

  ###### first-order surogate likelihood, **MEAN** ######
  ###### suppose the local data are stored in Xlocal,Ylocal ######
  SL = function(beta){
    -Lik(beta,Xlocal,Ylocal) - L%*%beta
  }

  ###### first-order surogate likelihood, **MEDIAN** ######
  SL_median = function(beta){
    -Lik(beta,Xlocal,Ylocal) - L_median%*%beta
  }

  ####### second-order gradient ######
  Lgradient2 =function(beta,X){
    design = cbind(1,X)
    Z=expit(design%*%beta)
    t(c(-Z*(1-Z))*design)%*%design/nrow(X)
  }

  ####### second-order surogate likelihood ######
  SL2 = function(beta){
    beta_temp = beta - betabar
    -Lik(beta,Xlocal,Ylocal) - L%*%beta - 0.5*t(beta_temp)%*%L2%*%t(t(beta_temp))
  }

  ####### second-order surogate likelihood, **MEDIAN** ######
  SL2_median = function(beta){
    beta_temp = beta - betabar
    -Lik(beta,Xlocal,Ylocal) - L_median%*%beta - 0.5*t(beta_temp)%*%L2_median%*%t(t(beta_temp))
  }

  # function to calcualte the hessian matrix for each site
  hessian_mat <- function(beta){
    L2 = matrix(0,nrow = K,ncol = p^2) # Store the second order gradient (each is a p*p matrix, we expand it into a vector)
    L2_local = Lgradient2(beta,Xlocal)

    nsite = rep(0,K)  #sample size in each site
    for (i in 1:K){
      Xsite = Xall[which(site==i),] #Predictors in the ith site
      Ysite = Yall[which(site==i)]                #Outcomes in the ith site
      nsite[i] = length(Ysite)
      L2[i,] = as.vector(Lgradient2(beta,Xsite))
    }
    return(rbind(as.vector(L2_local),L2))
  }


  ################################################
  # PART 1.. Initialization     #
  ###############################################

  #Local estimator is obtained for initialization
  fit0 = summary(glm(Ylocal~Xlocal, family = "binomial"(link = "logit")))
  beta0 = fit0$coefficients[,1]
  nam = c("Intercept", colnames(Xall))
  names(beta0) = nam
  #Then beta0 is passed to each site
  v_beta0 = fit0$cov.scaled#covariance matrix of beta0
  colnames(v_beta0) =nam
  rownames(v_beta0) =nam

  ##########################################
  # PART 2.. Calculation in each site     #
  #########################################
  #calculate the gradient in each site
  L = matrix(0,nrow = K,ncol = p) # Store the first order gradient (each is a p-dimensional vector)
  L2 = matrix(0,nrow = K,ncol = p^2) # Store the second order gradient (each is a p*p matrix, we expand it into a vector)
  betabar = beta0

  ##########Calculate the global second order gradient L2
  L2_local = Lgradient2(betabar,Xlocal)

  nsite = rep(0,K)  #sample size in each site
  for (i in 1:K){
    Xsite = Xall[which(site==i),] #Predictors in the ith site
    Ysite = Yall[which(site==i)]                #Outcomes in the ith site
    nsite[i] = length(Ysite)
    L[i,] = Lgradient(betabar,Xsite,Ysite)
    L2[i,] = as.vector(Lgradient2(betabar,Xsite))
  }
  hessian_1 = rbind(as.vector(L2_local),L2)
  #Each L[i,] is transfered to the local site for ODAL1
  # L[i,] and L2[i,] are transfered to the local site for ODAL2


  ##########################################
  # PART 3.. Local process    #
  #########################################
  ##########Calculate the global first order gradient L
  L_local = Lgradient(betabar,Xlocal,Ylocal)

  #median (median should be calcualted first, because "mean" updates L)
  L_all_median = apply(simplify2array(rbind(L,L_local)), 2, median)
  L_median = L_all_median - L_local

  #mean
  L_all = apply(rbind(diag(nsite)%*%L,L_local*nlocal),2,sum)/(sum(nsite)+nlocal)
  L =L_all-L_local


  #median (median should be calcualted first, because "mean" updates L)
  L2_all_median = apply(simplify2array(rbind(as.vector(L2_local),L2)), 2, median)
  L2_all_median = matrix(L2_all_median,ncol = p, nrow = p)
  L2_median = L2_all_median-L2_local

  #mean
  L2_all = apply(diag(c(nlocal,nsite))%*%rbind(as.vector(L2_local),L2),2,sum)/(sum(nsite)+nlocal)
  L2_all = matrix(L2_all,ncol = p, nrow = p)
  L2 = L2_all-L2_local



  ######################################
  # PART 4.. Estimation and Inference  #
  #####################################
  #Point estimation ODAL1
  ODAL1_mean = optim(betabar,SL,control = list(maxit = 10000,reltol = 1e-16))$par
  ODAL1_median = optim(betabar,SL_median,control = list(maxit = 10000,reltol = 1e-16))$par

  #Point estimation ODAL2
  ODAL2_mean = optim(betabar,SL2,control = list(maxit = 10000,reltol = 1e-16))$par
  ODAL2_median = optim(betabar,SL2_median,control = list(maxit = 10000,reltol = 1e-16))$par

  ####Variance: return the heissian matrix####
  hessian_ODAL1_local_mean = hessian_mat(ODAL1_mean)
  hessian_ODAL1_local_median = hessian_mat(ODAL1_median)
  hessian_ODAL2_local_mean = hessian_mat(ODAL2_mean)
  hessian_ODAL2_local_median = hessian_mat(ODAL2_median)

  ######################################
  # PART 5.. Pooled estimator  #
  #####################################
  fitall = summary(glm(Yall~Xall, family = "binomial"(link = "logit")))
  betaall = fitall$coefficients[,1]
  v_betaall = fitall$cov.scaled

  names(betaall) = nam
  colnames(v_betaall) =nam
  rownames(v_betaall) =nam


  ######################################
  # PART 6.. Meta estimator  #
  #####################################
  Beta = matrix(0,nrow = K,ncol = p) # Store the estimator in each site
  VBeta = matrix(0,nrow = K,ncol = p)# Store the covariance matrix from each site (each is a p*p matrix, we expand it into a vector)

  #fit logistic regression in each site
  for (i in 1:K){
    Xsite = Xall[which(site==i),] #Predictors in the ith site
    Ysite = Yall[which(site==i)]                #Outcomes in the ith site
    Beta[i,] = tryCatch(glm(Ysite~Xsite, family = "binomial"(link = "logit"))$coefficients,error=function(err) rep(NA,length(betaall)))
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
  # PART 7.. ODAL + Meta  #
  #####################################

  L = matrix(0,nrow = K,ncol = p) # Store the first order gradient (each is a p-dimensional vector)
  L2 = matrix(0,nrow = K,ncol = p^2)# Store the second order gradient (each is a p*p matrix, we expand it into a vector)
  betabar = betameta

  ##########Calculate the global second order gradient L2
  L2_local = Lgradient2(betabar,Xlocal)

  nsite = rep(0,K)  #sample size in each site
  for (i in 1:K){
    Xsite = Xall[which(site==i),] #Predictors in the ith site
    Ysite = Yall[which(site==i)]                #Outcomes in the ith site          #Outcomes in the ith site
    nsite[i] = length(Ysite)
    L[i,] = Lgradient(betabar,Xsite,Ysite)
    L2[i,] = as.vector(Lgradient2(betabar,Xsite))
  }
  hessian_2 = rbind(as.vector(L2_local),L2)
  #Each L[i,] is transfered to the local site for ODAL1
  # L[i,] and L2[i,] are transfered to the local site for ODAL2

  ##########Calculate the global first order gradient L
  L_local = Lgradient(betabar,Xlocal,Ylocal)

  #median (here median should be calcualted first, because "mean" updates L)
  L_all_median = apply(simplify2array(rbind(L,L_local)), 2, median)
  L_median = L_all_median - L_local

  #mean
  L_all = apply(rbind(diag(nsite)%*%L,L_local*nlocal),2,sum)/(sum(nsite)+nlocal)
  L =L_all-L_local


  #median (here median should be calcualted first, because "mean" updates L)
  L2_all_median = apply(simplify2array(rbind(as.vector(L2_local),L2)), 2, median)
  L2_all_median = matrix(L2_all_median,ncol = p, nrow = p)
  L2_median = L2_all_median-L2_local

  #mean
  L2_all = apply(diag(c(nlocal,nsite))%*%rbind(as.vector(L2_local),L2),2,sum)/(sum(nsite)+nlocal)
  L2_all = matrix(L2_all,ncol = p, nrow = p)
  L2 = L2_all-L2_local


  #Point estimation ODAL1
  meta_ODAL1_mean = optim(betabar,SL,control = list(maxit = 10000,reltol = 1e-16))$par
  meta_ODAL1_median = optim(betabar,SL_median,control = list(maxit = 10000,reltol = 1e-16))$par

  #Point estimation ODAL2
  meta_ODAL2_mean = optim(betabar,SL2,control = list(maxit = 10000,reltol = 1e-16))$par
  meta_ODAL2_median = optim(betabar,SL2_median,control = list(maxit = 10000,reltol = 1e-16))$par

  ####Variance: return the heissian matrix####
  hessian_ODAL1_meta_mean = hessian_mat(meta_ODAL1_mean)
  hessian_ODAL1_meta_median = hessian_mat(meta_ODAL1_median)
  hessian_ODAL2_meta_mean = hessian_mat(meta_ODAL2_mean)
  hessian_ODAL2_meta_median = hessian_mat(meta_ODAL2_median)

  ######################################
  # PART 8.. output  #
  #####################################
  output = list(beta_local = beta0,
                beta_pooled = betaall,
                beta_meta = betameta,
                ODAL1_local_mean = ODAL1_mean,
                ODAL1_local_median = ODAL1_median,
                ODAL2_local_mean = ODAL2_mean,
                ODAL2_local_median = ODAL2_median,
                ODAL1_meta_mean = meta_ODAL1_mean,
                ODAL1_meta_median = meta_ODAL1_median,
                ODAL2_meta_mean = meta_ODAL2_mean,
                ODAL2_meta_median = meta_ODAL2_median,
                var_matrix_local = as.vector(v_beta0),
                var_matrix_pooled = as.vector(v_betaall),
                var_meta = vmeta,
                hessian_1_betalocal = hessian_1,
                hessian_2_betameta = hessian_2,
                hessian_ODAL1_local_mean = hessian_ODAL1_local_mean,
                hessian_ODAL1_local_median = hessian_ODAL1_local_median,
                hessian_ODAL2_local_mean = hessian_ODAL2_local_mean,
                hessian_ODAL2_local_median = hessian_ODAL2_local_median,
                hessian_ODAL1_meta_mean = hessian_ODAL1_meta_mean,
                hessian_ODAL1_meta_median = hessian_ODAL1_meta_median,
                hessian_ODAL2_meta_mean = hessian_ODAL2_meta_mean,
                hessian_ODAL2_meta_median = hessian_ODAL2_meta_median)
  return(output)
}




evaluateOdal <- function(studyFolder, outcomeId, skipJmdc = FALSE, splitMdcr = FALSE) {
  writeLines(paste("Evaluating outcome", outcomeId))
  skipJmdc <- TRUE
  skipMdcr <- TRUE # Note: can only skip MDCR if also skipping JMDC
  outcomeId <- 3 # 5 = AMI, 3 = stroke
  data <- readRDS(file.path(studyFolder, sprintf("data_o%s.rds", outcomeId)))

  if (skipJmdc) {
    writeLines("Skipping JMDC")
    data <- data[data$database != "Jmdc", ]
  }
  if (skipMdcr) {
    writeLines("Skipping MDCR")
    data <- data[data$database != "mdcr", ]
  }

  data$age_in_years <- NULL
  data$`gender_=_FEMALE` <- NULL
  data$time <- NULL
  # data$Obesity <- NULL

  Yall <- data$y
  Xall <- as.matrix(data[, -which(colnames(data) %in% c("y", "database"))])
  site <- rep(0, nrow(data))
  # site[data$database == "Panther"] <- 1

  site[data$database == "mdcd"] <- 1
  site[data$database == "optum"] <- 2
  site[data$database == "mdcr"] <- 3
  site[data$database == "Jmdc"] <- 4
  # out = compare.methods(Xall,Yall,site)

  # pathToCsv <- system.file("settings", "CohortsToCreate.csv", package = "DistributedRegressionEval")
  # cohortsToCreate <- read.csv(pathToCsv)
  # outcomeName <- cohortsToCreate$name[cohortsToCreate$cohortId == outcomeId]

  if (outcomeId == 5) {
    outcomeName <- "AMI"
  } else if (outcomeId == 3) {
    outcomeName <- "stroke"
  }
  # siteName <- "ccae"
  # saveRDS(out, file.path(studyFolder, sprintf("output_ODAL_median_%s_%s%s.rds", outcomeName, siteName)))

  len_site = length(unique(site))
  i <- 0
  for (i in 0:(len_site - 1)) {
    if (i == 0) {
      siteName <- "ccae"
    } else if (i == 1) {
      siteName <- "mdcd"
    } else if (i == 2) {
      siteName <- "optum"
    } else if (i == 3) {
      siteName <- "mdcr"
    } else if (i == 4) {
      siteName <- "jmdc"
    }
    print(paste("I am working on ",i,"-th site as local site.",sep = ""))
    newSite = (site - i) %% len_site # change site number (0,1,2,3,...,K) to (K-i+1,..,K,0,1,2,...,K-i),
    out <- compare.methods(Xall,Yall,newSite)
    if (!skipJmdc & !skipMdcr) {
    postFix <- ""
    } else if (skipJmdc & !skipMdcr) {
      postFix <- "_noJmdc"
    } else if (skipJmdc & skipMdcr) {
      postFix <- "_noJmdcNorMdcr"
    }
    saveRDS(out, file.path(studyFolder, sprintf("output_ODAL_median_%s_%s%s.rds", outcomeName, siteName, postFix)))
  }
}



######################################################################################
##### use each site as the local site to run the function ############################
######################################################################################
# len_site = length(unique(site))
# output = list(list())
# for (i in 0:(len_site-1)){
#   print(paste("I am working on ",i,"-th site as local site.",sep = ""))
#   site = (site - i)%%len_site # change site number (0,1,2,3,...,K) to (K-i+1,..,K,0,1,2,...,K-i),
#   output[[i+1]] = compare.methods(Xall,Yall,site)
# }





