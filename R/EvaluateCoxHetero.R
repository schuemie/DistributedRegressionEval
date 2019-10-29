## this file is for running ODAC assuming heterogeous baseline hazard among sites
## please update ODACO package to v2.0,
##    DistCox() now accept input data as data frame
##    please make the first 3 columns of the input data as: site id, time, status
##    and the rest numeric covariates
## I'm still using the example of AMI / Stroke
## Db.names <- c('ccae',	'Jmdc',	'mdcd',	'optum')
## please also specify the var names for (time, status,  X1, X2) in the script below
## output are .csv files of beta and var est, and also a pdf of cumulative hazards plot


# use all 5 sites
require(KernSmooth)
Db.names <- c('ccae', 'Jmdc',		'mdcd',	'optum') # 'mdcr',
evaluateCoxHetero <- function(studyFolder, data, outcomeName){
  # the first 3 columns of the input data frame are: site id, time, status
  # and the rest are covariates (numeric)
  n.total <- nrow(data)
  Db.names <- unique(data$database)

  colnames(data) <- gsub("=", "", colnames(data))
  predictors <- colnames(data)
  predictors <- predictors[!(predictors %in% c("time", "y", "database"))]
  px <- length(predictors)
  form <- as.formula(paste("Surv(time, y) ~", paste(predictors, collapse = " + ")))

  pooledCox <- coxph(form, data=data)
  ## stratified Cox will be used as golden standard
  stratForm <- as.formula(paste("Surv(time, y) ~", paste(predictors, collapse = " + "), "+ strata(database)"))
  stratCox <- coxph(stratForm, data=data)

  # get localCox
  local.all <- list()
  time.max <- H0.max <- h0.max <- 0
  sum_inv_var <- matrix(0, px, px)
  sum_inv_var_b <- rep(0, px)

  for(idb in 1:length(Db.names)){
    localDb <- Db.names[idb]
    cat(localDb, '...\n')
    # Now don't need to transfer data frame to list
    localCox <- coxph(form, data[data$database == localDb, ])
    # prepare weighted est
    sum_inv_var <- sum_inv_var + solve(localCox$var)
    sum_inv_var_b <- sum_inv_var_b + solve(localCox$var, localCox$coef)

    local.all[[localDb]]$localCox <- localCox
    local.all[[localDb]]$n.localDb <- sum(data$database == localDb)
    # cumulative baseline haz
    H0 <- basehaz(localCox, centered=FALSE)
    local.all[[localDb]]$H0.localCox <- H0
    time.max <- max(time.max, max(H0$time))
    H0.max <- max(H0.max, max(H0$hazard))

    # set a common beta to calculate basehaz at each site, must be $coefficients, not $coef
    localCox$coefficients <- stratCox$coef
    H0 <- basehaz(localCox, centered=FALSE)
    local.all[[localDb]]$H0.stratCox <- H0

    # get baseline haz, keep h0 only at event times (h0=0 at censored times)
    h0 <- diff(H0$hazard)
    time0 <- H0$time
    time0 <- time0[-length(time0)]
    time0 = time0[h0!=0]
    h0 = h0[h0!=0]

    local.all[[localDb]]$h0 <- list(y=h0, x=time0)
    # kernal smoothing h0
    local.all[[localDb]]$h0.ks.20 <- ksmooth(time0, h0, bandwidth = 20)
    local.all[[localDb]]$h0.ks.50 <- ksmooth(time0, h0, bandwidth = 50)
    local.all[[localDb]]$h0.ks.normal <- ksmooth(time0, h0, bandwidth = 50, kernel = 'normal')
    h0.max <- max(h0.max, max(h0))

    # ss <- survfit(localCox)
  }

  # get avgCox: inverse-variance (i.e. hessian) weighted est
  var_avgCox <- solve(sum_inv_var)
  beta_avgCox <- var_avgCox %*% (sum_inv_var_b)


  # plot the cumulative hazard curve of each site
  # pdf(file.path(studyFolder, sprintf('cumhaz_all_sites_localCox_%s.pdf', outcomeName)))
  # for(idb in 1:length(Db.names)){
  #   localDb <- Db.names[idb]
  #   if(localDb==Db.names[1])
  #     plot(hazard ~ time, data=local.all[[localDb]]$H0.localCox, type='l',
  #          xlim=c(0, time.max), ylim=c(0, H0.max),
  #          xlab='time', ylab='cumhaz', main=paste0(outcomeName, ', beta_localCox') )
  #   else
  #     lines(hazard ~ time, data=local.all[[localDb]]$H0.localCox, col=idb, lty=idb)
  #
  # }
  # legend('topleft', legend=Db.names,  # x=0.1*time.max, y=0.9*H0.max
  #        cex=.8, col=1:length(Db.names), lty=1:length(Db.names))
  # dev.off()

  pdf(file.path(studyFolder, sprintf('basehaz_stratCox_%s.pdf', outcomeName)))
  for(idb in 1:length(Db.names)){
    localDb <- Db.names[idb]
    if(localDb==Db.names[1])
      # plot(hazard ~ time, data=local.all[[localDb]]$H0.stratCox, type='l',
      #      xlim=c(0, time.max), ylim=c(0, H0.max),
      #      xlab='time', ylab='cumhaz', main=paste0(outcomeName, ', beta_stratCox'))
      plot(local.all[[localDb]]$h0$y ~ local.all[[localDb]]$h0$x, type='s',
           xlim=c(0, time.max), ylim=c(0, h0.max),
           xlab='time', ylab='baseline hazard', main=paste0(outcomeName, ', beta_stratCox'))
    else
      lines(local.all[[localDb]]$h0$y ~ local.all[[localDb]]$h0$x, type='s',
            col=idb, lty=idb)
  }
  legend('topleft', legend=Db.names,
         cex=.8, col=1:length(Db.names), lty=1:length(Db.names))
  dev.off()

  ## ksmooth(bandwidth = 20), box kernel
  pdf(file.path(studyFolder, sprintf('basehaz_ks20_stratCox_%s.pdf', outcomeName)))
  for(idb in 1:length(Db.names)){
    localDb <- Db.names[idb]
    if(localDb==Db.names[1])
      plot(local.all[[localDb]]$h0.ks.20$y ~ local.all[[localDb]]$h0.ks.20$x, type='l',
           xlim=c(0, time.max), ylim=c(0, h0.max),
           xlab='time', ylab='baseline hazard', main=paste0(outcomeName, ', beta_stratCox (ksmooth, bw=20)'))
    else
      lines(local.all[[localDb]]$h0.ks.20$y ~ local.all[[localDb]]$h0.ks.20$x, type='l',
            col=idb, lty=idb)
  }
  legend('topleft', legend=Db.names,
         cex=.8, col=1:length(Db.names), lty=1:length(Db.names))
  dev.off()

  ## ksmooth(bandwidth = 50)
  pdf(file.path(studyFolder, sprintf('basehaz_ks50_stratCox_%s.pdf', outcomeName)))
  for(idb in 1:length(Db.names)){
    localDb <- Db.names[idb]
    if(localDb==Db.names[1])
      plot(local.all[[localDb]]$h0.ks.50$y ~ local.all[[localDb]]$h0.ks.50$x, type='l',
           xlim=c(0, time.max), ylim=c(0, h0.max),
           xlab='time', ylab='baseline hazard', main=paste0(outcomeName, ', beta_stratCox (ksmooth, bw=50)'))
    else
      lines(local.all[[localDb]]$h0.ks.50$y ~ local.all[[localDb]]$h0.ks.50$x, type='l',
            col=idb, lty=idb) #
  }
  legend('topleft', legend=Db.names,
         cex=.8, col=1:length(Db.names), lty=1:length(Db.names))
  dev.off()

  ## ksmooth(bandwidth = 50) using normal kernel
  pdf(file.path(studyFolder, sprintf('basehaz_ks50_normal_stratCox_%s.pdf', outcomeName)))
  for(idb in 1:length(Db.names)){
    localDb <- Db.names[idb]
    if(localDb==Db.names[1])
      plot(local.all[[localDb]]$h0.ks.normal$y ~ local.all[[localDb]]$h0.ks.normal$x, type='l',
           xlim=c(0, time.max), ylim=c(0, h0.max),
           xlab='time', ylab='baseline hazard', main=paste0(outcomeName, ', beta_stratCox (ksmooth, bw=50, normal)'))
    else
      lines(local.all[[localDb]]$h0.ks.normal$y ~ local.all[[localDb]]$h0.ks.normal$x, type='l',
            col=idb, lty=idb) #
  }
  legend('topleft', legend=Db.names,
         cex=.8, col=1:length(Db.names), lty=1:length(Db.names))
  dev.off()



  ## using each Db as Local and run DistCox
  ## now assume heterogeneous baseline hazards
  # Reformat data:
  sites <- unique(data$database)
  sites <- data.frame(database = sites,
                      id = 1:length(sites),
                      stringsAsFactors = FALSE)
  reformData <- merge(data, sites)
  reformData <- reformData[, c("id", "time", "y", predictors)]
  for (localId in sites$id) {
    localDb <- sites$database[sites$id == localId]
    cat(localDb, '...\n')

    # proposed DistCox, assume common baseline hazards (this is already studied...)
    distCoxAvgInit <- ODACO::DistCox(mydata = reformData,
                                     id.local = localId,
                                     init_est = beta_avgCox,
                                     output.ODACO1 = T,
                                     strat = F)

    # proposed DistCox, assume heterogeneous baseline hazards
    distCoxAvgInit.H <- ODACO::DistCox(mydata = reformData,
                                       id.local = localId,
                                       init_est = beta_avgCox,
                                       output.ODACO1 = T,
                                       strat = T)
    # assign(paste0('distCoxAvgInit.H.', localDb), distCoxAvgInit.H)

    tmp = rbind(local.all[[localDb]]$localCox$coef,
                c(beta_avgCox),
                c(distCoxAvgInit$beta_tilde1),
                c(distCoxAvgInit$beta_tilde),
                c(distCoxAvgInit.H$beta_tilde1),
                c(distCoxAvgInit.H$beta_tilde),
                pooledCox$coef,
                stratCox$coef,
                local.all[[localDb]]$n.localDb)
    results <- data.frame(description = c("Local",
                                          "Meta",
                                          "ODAC1",
                                          "ODAC2",
                                          "ODACH1",
                                          "ODACH2",
                                          "pooledCox",
                                          "stratCox",
                                          "sampleSize"), tmp)
    write.csv(results, file.path(studyFolder, sprintf("%s_beta_%s.csv", localDb, outcomeName)), row.names = FALSE)

    res.var <- rbind(local.all[[localDb]]$localCox$var,
                     var_avgCox,
                     distCoxAvgInit$sol1$hessian,
                     distCoxAvgInit$sol$hessian,
                     distCoxAvgInit.H$sol1$hessian,
                     distCoxAvgInit.H$sol$hessian,
                     pooledCox$var,
                     stratCox$var )

    write.csv(res.var, file.path(studyFolder, sprintf("%s_var_%s.csv", localDb, outcomeName)), row.names = FALSE)
  }
}


evaluateCoxHeteroAllOutcomes <- function(){
  library(survival)
  studyFolder <- "r:/DistributedCoxHeteroEval"
  pathToCsv <- system.file("settings", "KnownPredictors.csv", package = "DistributedRegressionEval")
  outcomes <- read.csv(pathToCsv)
  outcomes <- outcomes[!duplicated(outcomes$outcomeId), c("outcomeId", "outcomeName")]

  for (outcomeId in outcomes$outcomeId) {
    # outcomeId <- 5 # 3 = stroke, 5 = AMI
    outcomeName <- outcomes$outcomeName[outcomes$outcomeId == outcomeId]
    writeLines(paste("Evaluating outcome", outcomeName))

    data <- readRDS(file.path(studyFolder, sprintf("data_o%s.rds", outcomeId)))

    data$age_in_years <- (data$age_in_years - 50) / 20

    if (outcomeId == 3) {
      data <- data[data$database != "Jmdc", ]
    }

    # Remove for stroke model:
    # if (outcomeId == 3) {
    #   data$age_in_years <- NULL
    #   data$`gender_=_FEMALE` <- NULL
    #   data$Major_depressive_disorder <- NULL
    # }

    evaluateCoxHetero(studyFolder, data, outcomeName)
  }
}


