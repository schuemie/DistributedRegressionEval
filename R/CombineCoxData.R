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

#' @export
combineCoxData <- function(studyFolder, sampleSize = 1000000) {
  # sampleSize = 10000
  folders <- list.files(studyFolder, include.dirs = TRUE)
  data <- data.frame()
  for (folder in folders) {
    if (dir.exists(file.path(studyFolder, folder)) && file.exists(file.path(studyFolder, folder, "matchedPop.rds"))) {
      file <- file.path(studyFolder, folder, "matchedPop.rds")
      dataDb <- readRDS(file)
      if (nrow(dataDb) > sampleSize) {
        dataDb <- dataDb[sample.int(nrow(dataDb), sampleSize), ]
      }
      dataDb$database <- folder
      data <- plyr::rbind.fill(data, dataDb)
    }
  }
  saveRDS(data, file.path(studyFolder, "coxData.rds"))
}

createSummary <- function(studyFolder) {
  data <- readRDS(file.path(studyFolder, "coxData.rds"))
  createSummaryStats <- function(database) {
    dataDb <- data[data$database == database, ]
    stats <- data.frame(database = database,
                        treated = sum(dataDb$treatment == 1),
                        comparator = sum(dataDb$treatment == 0),
                        treatedEvents = sum(dataDb$treatment == 1 & dataDb$outcomeCount != 0),
                        comparatorEvents = sum(dataDb$treatment == 0 & dataDb$outcomeCount != 0))
    return(stats)
  }
  summaryStats <- lapply(unique(data$database), createSummaryStats)
  summaryStats <- do.call("rbind", summaryStats)

  write.csv(t(summaryStats), file.path(studyFolder, "Summary.csv"))
}
