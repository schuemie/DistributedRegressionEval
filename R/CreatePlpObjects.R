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


#' Create the PLP objects
#'
#' @details
#' This function will create the PLP objects.
#'
#' @param connectionDetails    An object of type \code{connectionDetails} as created using the
#'                             \code{\link[DatabaseConnector]{createConnectionDetails}} function in the
#'                             DatabaseConnector package.
#' @param cdmDatabaseSchema    Schema name where your patient-level data in OMOP CDM format resides.
#'                             Note that for SQL Server, this should include both the database and
#'                             schema name, for example 'cdm_data.dbo'.
#' @param cohortDatabaseSchema Schema name where intermediate data can be stored. You will need to have
#'                             write priviliges in this schema. Note that for SQL Server, this should
#'                             include both the database and schema name, for example 'cdm_data.dbo'.
#' @param cohortTable          The name of the table that will be created in the work database schema.
#'                             This table will hold the exposure and outcome cohorts used in this
#'                             study.
#' @param oracleTempSchema     Should be used in Oracle to specify a schema where the user has write
#'                             priviliges for storing temporary tables.
#' @param outputFolder         Name of local folder to place results; make sure to use forward slashes
#'                             (/)
#'
#' @export
createPlpObjects <- function(connectionDetails,
                             cdmDatabaseSchema,
                             cohortDatabaseSchema,
                             cohortTable = "cohort",
                             oracleTempSchema,
                             outputFolder) {
  pathToCsv <- system.file("settings", "KnownPredictors.csv", package = "DistributedRegressionEval")
  knownPredictors <- read.csv(pathToCsv)
  covariateSettings <- FeatureExtraction::createCovariateSettings(useDemographicsAgeGroup = TRUE,
                                                                  useDemographicsGender = TRUE,
                                                                  useConditionGroupEraLongTerm = TRUE,
                                                                  longTermStartDays = -365,
                                                                  endDays = 0,
                                                                  includedCovariateConceptIds = unique(knownPredictors$conceptId))
  plpData <- PatientLevelPrediction::getPlpData(connectionDetails = connectionDetails,
                                                cdmDatabaseSchema = cdmDatabaseSchema,
                                                oracleTempSchema = oracleTempSchema,
                                                cohortDatabaseSchema = cohortDatabaseSchema,
                                                cohortTable = cohortTable,
                                                cohortId = 102,
                                                washoutPeriod = 365,
                                                covariateSettings = covariateSettings,
                                                outcomeDatabaseSchema = cohortDatabaseSchema,
                                                outcomeTable = cohortTable,
                                                outcomeIds = c(3, 4, 5, 6),
                                                firstExposureOnly = TRUE)

  PatientLevelPrediction::savePlpData(plpData, file.path(outputFolder, "plpData"))

  # plpData <- PatientLevelPrediction::loadPlpData(file.path(outputFolder, "plpData"))
  covariateData <- FeatureExtraction::tidyCovariateData(covariates = plpData$covariates,
                                                        covariateRef = plpData$covariateRef,
                                                        populationSize = plpData$metaData$populationSize,
                                                        minFraction = 0,
                                                        normalize = TRUE,
                                                        removeRedundancy = TRUE)
  covariates <- ff::as.ram(ff::as.ram(covariateData$covariates))
  covariateRef <- ff::as.ram(plpData$covariateRef)
  s <- summary(plpData)
  outcomeCounts <- s$outcomeCounts
  outcomeCounts <- outcomeCounts[order(-outcomeCounts$personCount), ]
  safeOutcomeId <- outcomeCounts$outcomeId[1]

  for (outcomeId in s$metaData$outcomeIds) {
    population <- PatientLevelPrediction::createStudyPopulation(plpData,
                                                                outcomeId = outcomeId,
                                                                includeAllOutcomes = TRUE,
                                                                requireTimeAtRisk = TRUE ,
                                                                minTimeAtRisk =  365,
                                                                riskWindowStart = 1,
                                                                addExposureDaysToStart = FALSE,
                                                                riskWindowEnd = 366,
                                                                addExposureDaysToEnd = FALSE,
                                                                removeSubjectsWithPriorOutcome = TRUE)
    if (is.null(population)) {
      # population is set to null if there's no one with the outcome. Create study population with 'safe' outcome, then set
      # otucomeCount to 0:
      population <- PatientLevelPrediction::createStudyPopulation(plpData,
                                                                  outcomeId = safeOutcomeId,
                                                                  includeAllOutcomes = FALSE,
                                                                  requireTimeAtRisk = TRUE ,
                                                                  minTimeAtRisk =  365,
                                                                  riskWindowStart = 1,
                                                                  addExposureDaysToStart = FALSE,
                                                                  riskWindowEnd = 366,
                                                                  addExposureDaysToEnd = FALSE,
                                                                  removeSubjectsWithPriorOutcome = FALSE)
      population$outcomeCount <- 0
    }
    covariateIds <- covariateRef$covariateId[covariateRef$conceptId %in% knownPredictors$conceptId[knownPredictors$outcomeId == outcomeId]]
    covariateSubset <- covariates[covariates$covariateId %in% covariateIds & covariates$rowId %in% population$rowId, ]

    # Sparse to dense:
    ncovars <- length(covariateIds)
    nrows <- nrow(population)
    m <- matrix(0, nrows, ncovars)
    rowIs <- match(covariateSubset$rowId, population$rowId)
    columnIs <- match(covariateSubset$covariateId, covariateIds)
    for (i in 1:nrow(covariateSubset)) {
      m[rowIs[i], columnIs[i]] <- 1
    }
    data <- as.data.frame(m)
    columnNames <- as.character(covariateRef$covariateName[match(covariateIds, covariateRef$covariateId)])
    columnNames <- gsub(" ", "_", gsub(".*: ", "", columnNames))
    colnames(data) <- columnNames
    data$y <- as.integer(population$outcomeCount != 0)

    fileName <-  file.path(outputFolder, paste0("data_o", outcomeId, ".rds"))
    saveRDS(data, fileName)
  }
}
