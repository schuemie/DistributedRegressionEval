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
combineData <- function(studyFolder, dropColumns = FALSE, sampleSize = 100000) {
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
      if (dropColumns) {
        writeLines(paste("Dropping columns:", paste(columnsToDrop, collapse = ", ")))
        data <- data[, -which(colnames(data) %in% columnsToDrop)]
      } else {
        writeLines(paste("Columns", paste(columnsToDrop, collapse = ", "), "are all zero in one or more DBs"))
        for (i in 1:ncol(data)) {
          data[is.na(data[, i]), i] <- 0
        }
      }
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
  require(dplyr)
  for (outcomeId in outcomeIds) {
    data <- readRDS(file.path(studyFolder, sprintf("data_o%s.rds", outcomeId)))
    means <- data %>% select(-y) %>% group_by(database) %>% summarize_all(mean)
    counts <- data %>% group_by(database) %>% summarize(n = n(), outcomes = sum(y))
    counts <- rename(counts, !!paste0("o_", outcomeId) := outcomes)
    if (first) {
      summaryStats <- inner_join(means, counts)
      first <- FALSE
    } else {
      newCols <- names(means)[!(names(means) %in% names(summaryStats))]
      summaryStats <- inner_join(summaryStats, means %>% select(database, newCols))
      summaryStats <- inner_join(summaryStats, counts %>% select(-n))
    }
  }
  summaryStats <- summaryStats %>% select(database, n, o_3, o_4, o_5, o_6, everything())
  write.csv(t(summaryStats), file.path(studyFolder, "Summary.csv"))

}
