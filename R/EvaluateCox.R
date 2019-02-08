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



evaluateCox <- function(studyFolder) {
  data <- readRDS(file.path(studyFolder, "coxData.rds"))

  # Reformat for ODACO:
  data$status <- data$outcomeCount != 0
  data <- data[, c("survivalTime", "status", "treatment" , "database")]

  # localDb <- "panther"

  Db.names <- c('ccae',	'Jmdc',	'mdcd',	'optum',	'panther')
  for(localDb in Db.names){
    cat(localDb, '...\n')
    # coxph using local machine only
    dataLocal <- data[data$database == localDb, ]
    dataLocal$database <- NULL
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
  beta_avgCox <- (localCox.ccae$par*localCox.ccae$hessian +
                    localCox.Jmdc$par*localCox.Jmdc$hessian +
                    localCox.mdcd$par*localCox.mdcd$hessian +
                    localCox.optum$par*localCox.optum$hessian +
                    localCox.panther$par*localCox.panther$hessian) / (localCox.ccae$hessian+ localCox.Jmdc$hessian + localCox.mdcd$hessian + localCox.optum$hessian +localCox.panther$hessian)

  ## try using each Db as Local and run DistCox
  for(localDb in Db.names){
    cat(localDb, '...\n')
    otherDbs <- unique(data$database[data$database != localDb])

    dataLocal <- data[data$database == localDb, ]
    dataLocal$database <- NULL

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


  # use pooled data, coxph() in pkg 'survival'
  require(survival)
  pooledCox <- survival::coxph(Surv(dataCombined$all_data$t_surv, dataCombined$all_data$ind_event) ~ dataCombined$all_data$X)

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

  results <- data.frame(description = c("Local beta",
                                        "ODACO local init beta_tilde",
                                        "ODACO average init beta_tilde",
                                        "pooled beta_N",
                                        "pooled beta_pkg"),
                        ccae = c(localCox.ccae$par,
                                 distCoxLocalInit.ccae$beta_tilde,
                                 distCoxAvgInit.ccae$beta_tilde,
                                 distCoxAvgInit.ccae$beta_N,
                                 pooledCox$coef),
                        Jmdc = c(localCox.Jmdc$par,
                                 distCoxLocalInit.Jmdc$beta_tilde,
                                 distCoxAvgInit.Jmdc$beta_tilde,
                                 distCoxAvgInit.Jmdc$beta_N,
                                 pooledCox$coef),
                        mdcd = c(localCox.mdcd$par,
                                 distCoxLocalInit.mdcd$beta_tilde,
                                 distCoxAvgInit.mdcd$beta_tilde,
                                 distCoxAvgInit.mdcd$beta_N,
                                 pooledCox$coef),
                        optum = c(localCox.optum$par,
                                  distCoxLocalInit.optum$beta_tilde,
                                  distCoxAvgInit.optum$beta_tilde,
                                  distCoxAvgInit.optum$beta_N,
                                  pooledCox$coef),
                        panther = c(localCox.panther$par,
                                    distCoxLocalInit.panther$beta_tilde,
                                    distCoxAvgInit.panther$beta_tilde,
                                    distCoxAvgInit.panther$beta_N,
                                    pooledCox$coef))


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
                                sqrt(diag(solve(distCoxAvgInit.ccae$sol_N$hessian))/ n.ccae),
                                summary(pooledCox)$coef[,3]),
                       Jmdc = c(sqrt(diag(solve(localCox.Jmdc$hessian))/ n.Jmdc),
                                distCoxLocalInit.Jmdc$sigma_tilde,
                                distCoxAvgInit.Jmdc$sigma_tilde,
                                sqrt(diag(solve(distCoxAvgInit.Jmdc$sol_N$hessian))/ n.Jmdc),
                                summary(pooledCox)$coef[,3]),
                       mdcd = c(sqrt(diag(solve(localCox.mdcd$hessian))/ n.mdcd),
                                distCoxLocalInit.mdcd$sigma_tilde,
                                distCoxAvgInit.mdcd$sigma_tilde,
                                sqrt(diag(solve(distCoxAvgInit.mdcd$sol_N$hessian))/ n.mdcd),
                                summary(pooledCox)$coef[,3]),
                       optum = c(sqrt(diag(solve(localCox.optum$hessian))/ n.optum),
                                 distCoxLocalInit.optum$sigma_tilde,
                                 distCoxAvgInit.optum$sigma_tilde,
                                 sqrt(diag(solve(distCoxAvgInit.optum$sol_N$hessian))/ n.optum),
                                 summary(pooledCox)$coef[,3]),
                       panther = c(sqrt(diag(solve(localCox.panther$hessian))/ n.panther),
                                   distCoxLocalInit.panther$sigma_tilde,
                                   distCoxAvgInit.panther$sigma_tilde,
                                   sqrt(diag(solve(distCoxAvgInit.panther$sol_N$hessian))/ n.panther),
                                   summary(pooledCox)$coef[,3]))
  write.csv(Cox.sd, file.path(studyFolder, "sd-results.csv"), row.names = FALSE)

}
