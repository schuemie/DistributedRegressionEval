
##################################################
###Distributed algorithm for logistic regression##
##################################################
###############Updated on 12/21/2018##############

expit = function(x){exp(x)/(1+exp(x))}
library(MASS)
###########################################################################################
#The function "compare.methods()" returns the point estimation and variance estimation for the following methods
#1. Logistic regression using only local data
#2. Logistic regression using all data from multiple sites (pooled together).
#3. Distributed algorithm ODAL1 using local estiator (from #1) as initial value.
#4. Distributed algorithm ODAL2 using local estiator (from #1) as initial value.
#5. Fix effect Meta-analysis
#6. Distributed algorithm ODAL1 using meta estiator (from #5) as initial value.
#7. Distributed algorithm ODAL2 using meta estiator (from #5) as initial value.
###########################################################################################
compare.methods = function(Xall, Yall, site){
  Xlocal = Xall[which(site==0),]
  Ylocal = Yall[which(site==0)]
  K = length(unique(site))-1
  p = dim(Xall)[2]+1
  nlocal = length(Ylocal)
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

  #first order gradient
  Lgradient = function(beta,X,Y){
    design = cbind(1,X)
    t(Y-expit(design%*%t(t(beta))))%*%design/length(Y)
  }

  #first-order surogate likelihood, suppose the local data are stored in Xlocal,Ylocal
  SL = function(beta){
    -Lik(beta,Xlocal,Ylocal) - L%*%beta
  }

  #second-order gradient
  Lgradient2 =function(beta,X){
    design = cbind(1,X)
    Z=expit(design%*%beta)
    t(c(-Z*(1-Z))*design)%*%design/nrow(X)

  }

  #second-order surogate likelihood
  SL2 = function(beta){
    beta_temp = beta - betabar
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
    inv_L2 = ginv(mat_L2)
    out = inv_L2%*%mat_L1%*%inv_L2/N
    return(out)
  }

  ################################################
  # PART 1.. Initialization     #
  ###############################################

  #Local estimator is obtained for initialization
  fit0 = summary(glm(Ylocal~Xlocal, family = "binomial"(link = "logit")))
  beta0 = fit0$coefficients[,1]
  nam = c("Intercept", colnames(Xall))
  #nam = c("Intercept", "X1",    "X2",    "X3",    "X4")
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
  L2 = matrix(0,nrow = K,ncol = p^2)# Store the second order gradient (each is a p*p matrix, we expand it into a vector)
  betabar = beta0

  nsite = rep(0,K)  #sample size in each site
  for (i in 1:K){
    Xsite = Xall[which(site==i),] #Predictors in the ith site
    Ysite = Yall[which(site==i)]                #Outcomes in the ith site
    nsite[i] = length(Ysite)
    L[i,] = Lgradient(betabar,Xsite,Ysite)
    L2[i,] = as.vector(Lgradient2(betabar,Xsite))
  }
  #Each L[i,] is transfered to the local site for ODAL1
  # L[i,] and L2[i,] are transfered to the local site for ODAL2


  ##########################################
  # PART 3.. Local process    #
  #########################################


  #Calculate the global first order gradient L
  L_local = Lgradient(betabar,Xlocal,Ylocal)

  L_all = apply(rbind(diag(nsite)%*%L,L_local*nlocal),2,sum)/(sum(nsite)+nlocal)
  L =L_all-L_local
  #Calculate the global first order gradient L
  L2_local = Lgradient2(betabar,Xlocal)
  #Calculate the global second order gradient L2
  L2_all = apply(diag(c(nlocal,nsite))%*%rbind(as.vector(L2_local),L2),2,sum)/(sum(nsite)+nlocal)
  L2_all = matrix(L2_all,ncol = p, nrow = p)
  L2 = L2_all-L2_local




  ##########################################
  # PART 4.. Local + Newton #
  #########################################
  #LN = beta0 - ginv(L2_all)%*%t(L)


  ######################################
  # PART 4.. Estimation and Inference  #
  #####################################
  #Point estimation
  o1 <- optim(betabar,SL,control = list(maxit = 10000,reltol = 1e-16))
  ODAL1 = o1$par
  o2 <- optim(betabar,SL2,control = list(maxit = 10000,reltol = 1e-16))
  ODAL2 = o2$par

  ####Variance####
  #var1 = Sandwich(ODAL1,Xlocal,Ylocal,(sum(nsite)+nlocal))
  #var2 = Sandwich(ODAL2,Xlocal,Ylocal,(sum(nsite)+nlocal))
  #var_ODAL1 = diag(var1)
  #var_ODAL2 = diag(var2)


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

  nsite = rep(0,K)  #sample size in each site
  for (i in 1:K){
    Xsite = Xall[which(site==i),] #Predictors in the ith site
    Ysite = Yall[which(site==i)]                #Outcomes in the ith site          #Outcomes in the ith site
    nsite[i] = length(Ysite)
    L[i,] = Lgradient(betabar,Xsite,Ysite)
    L2[i,] = as.vector(Lgradient2(betabar,Xsite))
  }
  #Each L[i,] is transfered to the local site for ODAL1
  # L[i,] and L2[i,] are transfered to the local site for ODAL2

  #Calculate the global first order gradient L
  L_local = Lgradient(betabar,Xlocal,Ylocal)
  L_all = apply(rbind(diag(nsite)%*%L,L_local*nlocal),2,sum)/(sum(nsite)+nlocal)
  L = L_all-L_local

  #Calculate the global first order gradient L
  L2_local = Lgradient2(betabar,Xlocal)
  #Calculate the global second order gradient L2
  L2 = apply(diag(c(nlocal,nsite))%*%rbind(as.vector(L2_local),L2),2,sum)/(sum(nsite)+nlocal)
  L2_all = matrix(L2,ncol = p, nrow = p)
  L2 = L2_all-L2_local

  #Point estimation
  ODAL_meta1 = optim(betabar,SL,control = list(maxit = 10000,reltol = 1e-16))$par
  ODAL_meta2 = optim(betabar,SL2,control = list(maxit = 10000,reltol = 1e-16))$par

  ####Variance####
  #var1_odalmeta = Sandwich(ODAL_meta1,Xlocal,Ylocal,(sum(nsite)+nlocal))
  #var2_odalmeta = Sandwich(ODAL_meta2,Xlocal,Ylocal,(sum(nsite)+nlocal))

  #LN_meta = betameta - ginv(L2_all)%*%t(L)



  ######################################
  # PART 8.. output  #
  #####################################
  #output = list(estimation_local = beta0,var_local = v_beta0,
  # estimation_odal1 =ODAL1,var_odal1 = var1,
  #estimation_odal2 = ODAL2,var_odal2 = 0,
  # estimation_pooled = betaall,var_pooled = v_betaall,
  # estimation_meta = betameta,var_meta = vmeta,
  # estimation_metaodal1 = ODAL_meta1,var_metaodal1 = var1_odalmeta
  #,
  #estimation_metaodal2 = ODAL_meta2,var_metaodal2 = 0
  # )

  # out = c(beta0,ODAL1,ODAL2, betaall,betameta, ODAL_meta1,ODAL_meta2,vmeta,as.vector(v_betaall),as.vector(v_beta0))
  out <- rbind(beta0, ODAL1, ODAL2, betaall, betameta, ODAL_meta1, ODAL_meta2)
  return(out)
}


######################################################################################
evaluateOdal <- function(studyFolder, outcomeId, skipJmdc = FALSE, splitCcae = FALSE) {
  writeLines(paste("Evaluating outcome", outcomeId))
  # outcomeId <- 5
  data <- readRDS(file.path(studyFolder, sprintf("data_o%s.rds", outcomeId)))
  # normalize(data)
  if (skipJmdc) {
    writeLines("Skipping JMDC")
    data <- data[data$database != "Jmdc", ]
  }
  if (splitCcae) {
    writeLines("Splitting CCAE")
    data <- data[data$database == "ccae", ]
    set.seed(123)
    data$subsetId <- sample.int(5, nrow(data), replace = TRUE)
    idx <- data$subsetId != 1
    data$database[idx] <- paste(data$database[idx], data$subsetId[idx])
    data$subsetId <- NULL
  }
  # data$Alcohol_dependence <- NULL
  # data$Major_depressive_disorder <- NULL
  data$age_in_years <- NULL
  data$`gender_=_FEMALE` <- NULL
  data$time <- NULL
  # data$Obesity <- NULL

  Yall <- data$y
  Xall <- as.matrix(data[, -which(colnames(data) %in% c("y", "database"))])
  site <- rep(0, nrow(data))
  # site[data$database == "Panther"] <- 1
  site[data$database == "Jmdc"] <- 1
  site[data$database == "mdcd"] <- 2
  site[data$database == "optum"] <- 3
  site[data$database == "mdcr"] <- 4
  out = compare.methods(Xall,Yall,site)

  pathToCsv <- system.file("settings", "CohortsToCreate.csv", package = "DistributedRegressionEval")
  cohortsToCreate <- read.csv(pathToCsv)
  outcomeName <- cohortsToCreate$name[cohortsToCreate$cohortId == outcomeId]

  postfix <- ""
  if (skipJmdc) {
    postfix <- "_noJmdc"
  }
  if (splitCcae) {
    postfix <- "_splitCcae"
  }
  write.csv(out, file.path(studyFolder, sprintf("output_%s%s.csv", outcomeName, postfix)))

  # library(fmsb)
  # vizData <- as.data.frame(out)
  # vizData <- vizData[rev(c("betaall", "ODAL1", "ODAL2", "betameta")), ]
  # # vizData <- vizData[rev(c("betaall", "ODAL2", "betameta")), ]
  # vizData <- vizData[, 2:ncol(vizData)]
  # # colors_border = rev(c(rgb(0,0,0),  rgb(0.2,0.5,0.5,0.8), rgb(0.8,0.2,0.5,0.8) ))
  # # colors_in = rev(c(rgb(0,0,0, 0), rgb(0.2,0.5,0.5,0.4), rgb(0.8,0.2,0.5,0.4)))
  # colors_border = rev(c(rgb(0,0,0),  rgb(0.2,0.5,0.5,0.7), rgb(0.8,0.2,0.5,0.7) , rgb(0.7,0.5,0.1,0.7) ))
  # colors_in = rev(c(rgb(0,0,0, 0), rgb(0.2,0.5,0.5,0.25), rgb(0.8,0.2,0.5,0.25) , rgb(0.7,0.5,0.1,0.25) ))
  # # vizData <- rbind(apply(vizData, 2, max), apply(vizData, 2, min), vizData)
  # vizData <- rbind(rep(max(vizData), ncol(vizData)), rep(min(vizData), ncol(vizData)), vizData)
  #
  # radarchart( vizData  , axistype=1 ,
  #             #custom polygon
  #             pcol=colors_border , pfcol=colors_in , plwd=4 , plty=c(1,1,1,2),
  #             #custom the grid
  #             cglcol="grey", cglty=1, axislabcol="grey", caxislabels=format(seq(min(vizData),max(vizData), length.out = 5), digits = 2), cglwd=0.8,
  #             #custom labels
  #             vlcex=0.8
  # )
  # legend(x=1, y=1, legend = rownames(vizData[-c(1,2),]), bty = "n", pch=20 , col=colors_border , text.col = "grey", cex=1.2, pt.cex=3)

}


normalize <- function(data) {
  for (i in 1:ncol(data)) {
    if (is.numeric(data[, i]) && max(data[, i]) != 0) {
      data[, i] <- data[, i] / max(data[, i])
    }
  }
}


