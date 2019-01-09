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
  localDb <- "panther"
  otherDbs <- unique(data$database[data$database != localDb])

  # Reformat for ODACO:
  data$status <- data$outcomeCount != 0
  data <- data[, c("survivalTime", "status", "treatment" , "database")]
  dataLocal <- data[data$database == localDb, ]
  dataLocal$database <- NULL
  dataRemote <- data[data$database != localDb, ]
  database <- dataRemote$database
  dataRemote$database <- NULL
  dataRemote <- split(dataRemote, database)

  # transfer data frame to list
  dataCombined <- ODACO::combine_data(dataLocal, dataRemote, col_time = 1, col_event = 2)

  # coxph using local machine only
  localCox <- ODACO::my_coxph(dataCombined$local_data)

  # proposed DistCox, use likelihood from local and gradient from remote
  distCoxLocalInit <- ODACO::DistCox(local_data = dataCombined$local_data,
                            all_data = dataCombined$all_data,
                            init_est = 'local')

  distCoxAvgInit <- ODACO::DistCox(local_data = dataCombined$local_data,
                            all_data = dataCombined$all_data,
                            init_est = 'all')

  # use pooled data, coxph() in pkg 'survival'
  pooledCox <- survival::coxph(Surv(dataCombined$all_data$t_surv, dataCombined$all_data$ind_event) ~ dataCombined$all_data$X)

  # comparison: DistCox (beta_tilde) obtains coef est better than using local only, closer to use pooled data
  results <- data.frame(description = c("Local beta",
                                        "ODACO local init beta_tilde",
                                        "ODACO local init beta_N",
                                        "ODACO average init beta_tilde",
                                        "ODACO average init beta_N",
                                        "pooled beta"),
                        beta = c(localCox$par,
                                 distCoxLocalInit$beta_tilde,
                                 distCoxLocalInit$beta_N,
                                 distCoxAvgInit$beta_tilde,
                                 distCoxAvgInit$beta_N,
                                 pooledCox$coef))
  write.csv(results, file.path(studyFolder, "results.csv"), row.names = FALSE)

  # compare s.d. estimates
  cbind(distCox$sigma_tilde,
        summary(pooledCox)$coef[,3] )
}
