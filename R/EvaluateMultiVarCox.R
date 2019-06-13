evaluateCox <- function(studyFolder, outcomeId) {
  writeLines(paste("Evaluating outcome", outcomeId))
  # outcomeId <- 5 # 3 = stroke, 5 = AMI
  data <- readRDS(file.path(studyFolder, sprintf("data_o%s.rds", outcomeId)))

  # Remove for stroke model:
  if (outcomeId == 3) {
    data$age_in_years <- NULL
    data$`gender_=_FEMALE` <- NULL
    data$Major_depressive_disorder <- NULL
  }

  # Reformat for ODACO:
  data$status <- data$y != 0
  predictors <- colnames(data)[1:(which(colnames(data) == "y")-1)]

  data <- data[, c("time", "status", predictors, "database")]
  n.total <- nrow(data)
  # localDb <- "panther"

  Db.names <- c('ccae',	'Jmdc',	'mdcd',	'optum',	'mdcr')
  for(localDb in Db.names){
    cat(localDb, '...\n')
    # coxph using local machine only
    otherDbs <- unique(data$database[data$database != localDb])

    dataLocal <- data[data$database == localDb, ]
    dataLocal$database <- NULL
    assign(paste0('n.', localDb), nrow(dataLocal))

    dataRemote <- data[data$database != localDb, ]
    database <- dataRemote$database
    dataRemote$database <- NULL
    dataRemote <- split(dataRemote, database)

    # transfer data frame to list
    dataCombined <- ODACO::combine_data(dataLocal, dataRemote, col_time = 1, col_event = 2)

    localCox <- ODACO::my_coxph(dataCombined$local_data)
    assign(paste0('localCox.', localDb), localCox)
  }
  # inverse-variance (i.e. hessian) weighted est
  # beta_avgCox <- (localCox.ccae$par*localCox.ccae$hessian +
  #                   localCox.Jmdc$par*localCox.Jmdc$hessian +
  #                   localCox.mdcd$par*localCox.mdcd$hessian +
  #                   localCox.optum$par*localCox.optum$hessian +
  #                   localCox.mdcr$par*localCox.mdcr$hessian) / (localCox.ccae$hessian+ localCox.Jmdc$hessian + localCox.mdcd$hessian + localCox.optum$hessian +localCox.mdcr$hessian)
  beta_avgCox <- solve((localCox.ccae$hessian+ localCox.Jmdc$hessian + localCox.mdcd$hessian + localCox.optum$hessian +     localCox.mdcr$hessian),
                       t(localCox.ccae$par%*%localCox.ccae$hessian +
                          localCox.Jmdc$par %*%localCox.Jmdc$hessian +
                          localCox.mdcd$par%*%localCox.mdcd$hessian +
                          localCox.optum$par%*%localCox.optum$hessian +
                          localCox.mdcr$par%*%localCox.mdcr$hessian) )


  # use pooled data, coxph() in pkg 'survival'
  pooledCox <- survival::coxph(Surv(dataCombined$all_data$t_surv, dataCombined$all_data$ind_event) ~ dataCombined$all_data$X)

  ## try using each Db as Local and run DistCox
  Big.Db.names <- c('ccae',	'mdcd',	'optum',	'mdcr')
  for(localDb in Big.Db.names){
  # localDb <- "ccae"
    cat(localDb, '...\n')
    otherDbs <- unique(data$database[data$database != localDb])

    dataLocal <- data[data$database == localDb, ]
    dataLocal$database <- NULL
    assign(paste0('n.', localDb), nrow(dataLocal))

    dataRemote <- data[data$database != localDb, ]
    database <- dataRemote$database
    dataRemote$database <- NULL
    dataRemote <- split(dataRemote, database)

    # transfer data frame to list
    dataCombined <- ODACO::combine_data(dataLocal, dataRemote, col_time = 1, col_event = 2)

    # proposed DistCox, use likelihood from local and gradient from remote
    distCoxLocalInit <- ODACO::DistCox(local_data = dataCombined$local_data,
                                       all_data = dataCombined$all_data,
                                       init_est = 'local')
    assign(paste0('distCoxLocalInit.', localDb), distCoxLocalInit)

    ## this was a mistake: init_est = 'all' use pooledCox instead of Avg as init, plz ignore this option
    # distCoxAvgInit <- ODACO::DistCox(local_data = dataCombined$local_data,
    #                                  all_data = dataCombined$all_data,
    #                                  init_est = 'all')

    ## Use inverse-variance weighted localCox est as init
    distCoxAvgInit <- ODACO::DistCox(local_data = dataCombined$local_data,
                                     all_data = dataCombined$all_data,
                                     init_est = beta_avgCox)
    assign(paste0('distCoxAvgInit.', localDb), distCoxAvgInit)

  }




  # comparison: DistCox (beta_tilde) obtains coef est better than using local only, closer to use pooled data
  for (localDb in Big.Db.names){
    results <- rbind(c(get(paste0("localCox.", localDb))$par),
                     c(beta_avgCox),
                     c(get(paste0("distCoxLocalInit.", localDb))$beta_tilde),
                     c(get(paste0("distCoxAvgInit.", localDb))$beta_tilde),
                     c(get(paste0("distCoxLocalInit.", localDb))$beta_tilde2),
                     c(get(paste0("distCoxAvgInit.", localDb))$beta_tilde2),
                     c(pooledCox$coef))
    row.names(results) <- c("Local beta",
                            "Average beta",
                            "ODACO local init beta_tilde",
                            "ODACO average init beta_tilde",
                            "ODACO2 local init beta_tilde",
                            "ODACO2 average init beta_tilde",
                            "pooled beta_pkg")
    write.csv(results, file.path(studyFolder, sprintf("output_%s.csv", localDb)))
  }

  sd_avgCox <- sqrt(diag(solve((localCox.ccae$hessian*n.ccae+ localCox.Jmdc$hessian*n.Jmdc + localCox.mdcd$hessian*n.mdcd
                                + localCox.optum$hessian*n.optum + localCox.mdcr$hessian*n.mdcr))))

  for (localDb in Big.Db.names){
    Cox.sd <- rbind(c(sqrt(diag(solve(get(paste0("localCox.", localDb))$hessian))/ get(paste0('n.', localDb)))),
                    c(sd_avgCox),
                    c(get(paste0("distCoxLocalInit.", localDb))$sigma_tilde),
                    c(get(paste0("distCoxAvgInit.", localDb))$sigma_tilde),
                    c(get(paste0("distCoxLocalInit.", localDb))$sigma_tilde2),
                    c(get(paste0("distCoxAvgInit.", localDb))$sigma_tilde2),
                    c(summary(pooledCox)$coef[,3]))
    row.names(Cox.sd) <- c("Local sd",
                           "Average sd",
                           "ODACO local init sd_tilde",
                           "ODACO average init sd_tilde",
                           "ODACO2 local init sd_tilde",
                           "ODACO2 average init sd_tilde",
                           "pooled sd_pkg")
    write.csv(Cox.sd, file.path(studyFolder, sprintf("sd_output_%s.csv", localDb)))
  }

  for (localDb in Big.Db.names){
    local.hess <- get(paste0("distCoxLocalInit.", localDb))$sol2$hessian
    avg.hess  <- get(paste0("distCoxAvgInit.", localDb))$sol2$hessian
    write.csv(rbind(local.hess, avg.hess),
              file.path(studyFolder, sprintf("hessian_output_%s.csv", localDb)))

  }

}
