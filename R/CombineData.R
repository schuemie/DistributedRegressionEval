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

combineData <- function(studyFolder, sampleSize = 100000) {
  # sampleSize = 10000
  folders <- list.files(studyFolder, include.dirs = TRUE)
  outcomeIds <- c(3, 4, 5, 6)
  for (outcomeId in outcomeIds) {
    data <- data.frame()
    columnsToDrop <- c()
    for (folder in folders) {
      if (dir.exists(file.path(studyFolder, folder))) {
        file <- file.path(studyFolder, folder, sprintf("data_o%s.rds", outcomeId))
        dataDb <- readRDS(file)
        if (nrow(dataDb) > sampleSize) {
          dataDb <- dataDb[sample.int(nrow(dataDb), sampleSize), ]
        }
        for (i in 1:ncol(dataDb)) {
          if (sum(dataDb[, i]) == 0 && colnames(dataDb)[i] != "y") {
            columnsToDrop <- c(columnsToDrop, colnames(dataDb)[i])
          }
        }
        dataDb$database <- folder
        data <- plyr::rbind.fill(data, dataDb)

      }
    }
    for (i in 1:ncol(data)) {
       if (any(is.na(data[, i]))) {
         columnsToDrop <- c(columnsToDrop, colnames(data)[i])
       }
    }
    if (length(columnsToDrop) > 0) {
      data <- data[, -which(colnames(data) %in% columnsToDrop)]
    }
    saveRDS(data, file.path(studyFolder, sprintf("data_o%s.rds", outcomeId)))
  }
}

createSummary <- function(studyFolder) {
  outcomeIds <- c(3, 4, 5, 6)
  pathToCsv <- system.file("settings", "CohortsToCreate.csv", package = "DistributedRegressionEval")
  cohortsToCreate <- read.csv(pathToCsv)
  first <- TRUE
  summaryStats <- NULL
  for (outcomeId in outcomeIds) {
    data <- readRDS(file.path(studyFolder, sprintf("data_o%s.rds", outcomeId)))
    if (first) {
      summaryStats <- aggregate(y ~ database, data, length)
      colnames(summaryStats)[2] <- "subjects"
      first <- FALSE
    }
    outcomeCount <- aggregate(y ~ database, data, sum)
    colnames(outcomeCount)[2] <- as.character(cohortsToCreate$name[cohortsToCreate$cohortId == outcomeId])
    summaryStats <- merge(summaryStats, outcomeCount)
  }
  write.csv(summaryStats, file.path(studyFolder, "summaryStats.csv"), row.names = FALSE)
}