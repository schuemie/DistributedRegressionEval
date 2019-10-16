## this file is for running ODAC assuming heterogeous baseline hazard among sites
## please update ODACO package to v2.0,
##    DistCox() now accept input data as data frame
##    please make the first 3 columns of the input data as: site id, time, status
##    and the rest numeric covariates
## I'm still using the example of AMI / Stroke
## Db.names <- c('ccae',	'Jmdc',	'mdcd',	'optum')
## please also specify the var names for (time, status,  X1, X2) in the script below
## output are .csv files of beta and var est, and also a pdf of cumulative hazards plot

# Jdbc has convergence issue, can you use mdcr?
Db.names <- c('ccae',	'mdcr',	'mdcd',	'optum')
evaluateCoxHetero <- function(studyFolder, data, outcomeName){
  # the first 3 columns of the input data frame are: site id, time, status
  # and the rest are covariates (numeric)
  n.total <- nrow(data)
  Db.names <- unique(data$database)

  colnames(data) <- gsub("=", "", colnames(data))
  predictors <- colnames(data)
  predictors <- predictors[!(predictors %in% c("time", "y", "database"))]
  form <- as.formula(paste("Surv(time, y) ~", paste(predictors, collapse = " + ")))

  pooledCox <- coxph(form, data=data)
  ## stratified Cox will be used as golden standard
  stratForm <- as.formula(paste("Surv(time, y) ~", paste(predictors, collapse = " + "), "+ strata(database)"))
  stratCox <- coxph(stratForm, data=data)

  # get localCox
  # also plot the cumulative hazard curve of each site
  pdf(file.path(studyFolder, sprintf('cumhaz_all_sites_%s.pdf', outcomeName)))
  for(idb in 1:length(Db.names)){
    localDb <- Db.names[idb]
    cat(localDb, '...\n')
    # Now don't need to transfer data frame to list
    localCox <- coxph(form, data[data$database == localDb, ])

    assign(paste0('localCox.', localDb), localCox)
    assign(paste0("n.", localDb), sum(data$database == localDb))

    ss <- survfit(localCox)
    if(localDb==Db.names[1])
      plot(ss$cumhaz ~ ss$time, type='l', xlab='time', ylab='cumhaz', main=outcomeName)
    else
      lines(ss$cumhaz~ ss$time, col=idb, lty=idb)
    # plot(ss, fun='cumhaz', main=paste0(localDb, ", cumulative hazard") )
  }
  legend(x=0.1*max(ss$time), y=0.9*max(ss$cumhaz), Db.names,
         cex=.8, col=1:length(Db.names), lty=1:length(Db.names))
  dev.off()


  # get avgCox: inverse-variance (i.e. hessian) weighted est
  var_avgCox <- solve(solve(localCox.ccae$var)+ solve(localCox.mdcr$var) +
                        solve(localCox.mdcd$var) + solve(localCox.optum$var))
  beta_avgCox <- var_avgCox %*% (solve(localCox.ccae$var, localCox.ccae$coef)+ solve(localCox.mdcr$var, localCox.mdcr$coef) +
                                   solve(localCox.mdcd$var,localCox.mdcd$coef) + solve(localCox.optum$var,localCox.optum$coef))

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

    tmp = rbind(get(paste0('localCox.', localDb))$coef,
                c(beta_avgCox),
                c(distCoxAvgInit$beta_tilde1),
                c(distCoxAvgInit$beta_tilde),
                c(distCoxAvgInit.H$beta_tilde1),
                c(distCoxAvgInit.H$beta_tilde),
                pooledCox$coef,
                stratCox$coef,
                get(paste0("n.", localDb)))
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

    res.var <- rbind(get(paste0('localCox.', localDb))$var,
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

  outcomeId <- 3 # 3 = stroke, 5 = AMI
  writeLines(paste("Evaluating outcome", outcomeId))

  data <- readRDS(file.path(studyFolder, sprintf("data_o%s.rds", outcomeId)))

  # Remove for stroke model:
  # if (outcomeId == 3) {
  #   data$age_in_years <- NULL
  #   data$`gender_=_FEMALE` <- NULL
  #   data$Major_depressive_disorder <- NULL
  # }

  if (outcomeId == 5) {
    outcomeName <- "AMI"
  } else if (outcomeId == 3) {
    outcomeName <- "stroke"
  }
  evaluateCoxHetero(studyFolder, data, outcomeName)

}


