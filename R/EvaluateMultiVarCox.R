evaluateCox <- function(studyFolder, outcomeId) {
  writeLines(paste("Evaluating outcome", outcomeId))
  # outcomeId <- 5
  data <- readRDS(file.path(studyFolder, sprintf("data_o%s.rds", outcomeId)))
  data$age_in_years <- NULL
  data$`gender_=_FEMALE` <- NULL
  data$Major_depressive_disorder <- NULL

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
  # for(localDb in Db.names){
  localDb <- "ccae"
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

  # }




  # comparison: DistCox (beta_tilde) obtains coef est better than using local only, closer to use pooled data
  # results <- data.frame(description = c("Local beta",
  #                                       "ODACO local init beta_tilde",
  #                                       "ODACO local init beta_N",
  #                                       "ODACO average init beta_tilde",
  #                                       "ODACO average init beta_N",
  #                                       "pooled beta"),
  #                       beta = c(localCox$par,
  #                                distCoxLocalInit$beta_tilde,
  #                                distCoxLocalInit$beta_N,
  #                                distCoxAvgInit$beta_tilde,
  #                                distCoxAvgInit$beta_N,
  #                                pooledCox$coef))

  tmp <- rbind(get(paste0('localCox.', 'ccae'))$par,
               c(beta_avgCox),
  c(distCoxLocalInit$beta_tilde),
  c(distCoxAvgInit$beta_tilde),
  c(distCoxLocalInit$beta_tilde2),
  c(distCoxAvgInit$beta_tilde2),
  # distCoxAvgInit.ccae$beta_N,
  pooledCox$coef)
  row.names(tmp) <- c("Local beta",
                      "Average beta",
                      "ODACO local init beta_tilde",
                      "ODACO average init beta_tilde",
                      "ODACO2 local init beta_tilde",
                      "ODACO2 average init beta_tilde",
                      # "pooled beta_N",
                      "pooled beta_pkg")
  write.csv(tmp, file.path(studyFolder, "output.csv"))

  return(NULL)
  # Stopping here ------------------------

  results <- data.frame(description = c("Local beta",
                                        "Average beta",
                                        "ODACO local init beta_tilde",
                                        "ODACO average init beta_tilde",
                                        "ODACO2 local init beta_tilde",
                                        "ODACO2 average init beta_tilde",
                                        # "pooled beta_N",
                                        "pooled beta_pkg"),
                        tmp)


  write.csv(results, file.path(studyFolder, "results.csv"), row.names = FALSE)

  # compare s.d. estimates
  # cbind(distCox$sigma_tilde,
  #       summary(pooledCox)$coef[,3] )
  Cox.sd <- data.frame(description = c("Local sd",
                                       "ODACO local init sd_tilde",
                                       "ODACO average init sd_tilde",
                                       "pooled sd_N",
                                       "pooled sd_pkg"),
                       ccae = c(sqrt(diag(solve(localCox.ccae$hessian))/ n.ccae),
                                distCoxLocalInit.ccae$sigma_tilde,
                                distCoxAvgInit.ccae$sigma_tilde,
                                sqrt(diag(solve(distCoxAvgInit.ccae$sol_N$hessian))/ n.total),
                                summary(pooledCox)$coef[,3]),
                       Jmdc = c(sqrt(diag(solve(localCox.Jmdc$hessian))/ n.Jmdc),
                                distCoxLocalInit.Jmdc$sigma_tilde,
                                distCoxAvgInit.Jmdc$sigma_tilde,
                                sqrt(diag(solve(distCoxAvgInit.Jmdc$sol_N$hessian))/ n.total),
                                summary(pooledCox)$coef[,3]),
                       mdcd = c(sqrt(diag(solve(localCox.mdcd$hessian))/ n.mdcd),
                                distCoxLocalInit.mdcd$sigma_tilde,
                                distCoxAvgInit.mdcd$sigma_tilde,
                                sqrt(diag(solve(distCoxAvgInit.mdcd$sol_N$hessian))/ n.total),
                                summary(pooledCox)$coef[,3]),
                       optum = c(sqrt(diag(solve(localCox.optum$hessian))/ n.optum),
                                 distCoxLocalInit.optum$sigma_tilde,
                                 distCoxAvgInit.optum$sigma_tilde,
                                 sqrt(diag(solve(distCoxAvgInit.optum$sol_N$hessian))/ n.total),
                                 summary(pooledCox)$coef[,3]),
                       panther = c(sqrt(diag(solve(localCox.panther$hessian))/ n.panther),
                                   distCoxLocalInit.panther$sigma_tilde,
                                   distCoxAvgInit.panther$sigma_tilde,
                                   sqrt(diag(solve(distCoxAvgInit.panther$sol_N$hessian))/ n.total),
                                   summary(pooledCox)$coef[,3]))
  write.csv(Cox.sd, file.path(studyFolder, "sd-results.csv"), row.names = FALSE)
}
