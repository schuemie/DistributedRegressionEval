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


#' Create the Cox objects
#'
#' @details
#' This function will create the Cox objects.
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
createCoxObjects <- function(connectionDetails,
                             cdmDatabaseSchema,
                             cohortDatabaseSchema,
                             cohortTable = "cohort",
                             oracleTempSchema,
                             outputFolder) {

  covariateSettings <- FeatureExtraction::createDefaultCovariateSettings(excludedCovariateConceptIds = 21603933,
                                                                         addDescendantsToExclude = TRUE)
  cmData <- CohortMethod::getDbCohortMethodData(connectionDetails = connectionDetails,
                                                cdmDatabaseSchema = cdmDatabaseSchema,
                                                oracleTempSchema = oracleTempSchema,
                                                outcomeDatabaseSchema = cohortDatabaseSchema,
                                                outcomeTable = cohortTable,
                                                covariateSettings = covariateSettings,
                                                targetId = 1118084,
                                                comparatorId = 1124300,
                                                outcomeIds = 7,
                                                excludeDrugsFromCovariates = FALSE,
                                                removeDuplicateSubjects = TRUE,
                                                restrictToCommonPeriod = TRUE,
                                                washoutPeriod = 365,
                                                firstExposureOnly = TRUE,
                                                maxCohortSize = 100000)

  CohortMethod::saveCohortMethodData(cmData, file.path(outputFolder, "cmData"), compress = TRUE)

  # cmData <- CohortMethod::loadCohortMethodData(file.path(outputFolder, "cmData"))

  studyPop <- CohortMethod::createStudyPopulation(cohortMethodData = cmData,
                                                  outcomeId = 7,
                                                  removeSubjectsWithPriorOutcome = TRUE,
                                                  minDaysAtRisk = 1,
                                                  riskWindowStart = 0,
                                                  addExposureDaysToStart = FALSE,
                                                  riskWindowEnd = 30,
                                                  addExposureDaysToEnd = TRUE)

  ps <- CohortMethod::createPs(cohortMethodData = cmData,
                               population = studyPop,
                               prior = Cyclops::createPrior("laplace", exclude = c(0), useCrossValidation = TRUE),
                               control = Cyclops::createControl(cvType = "auto",
                                                                startingVariance = 0.01,
                                                                noiseLevel = "quiet",
                                                                tolerance = 2e-07,
                                                                cvRepetitions = 1,
                                                                threads = 10))
  saveRDS(ps, file.path(outputFolder, "ps.rds"))

  matchedPop <- CohortMethod::matchOnPs(ps, maxRatio = 1)

  # strataSizes <- aggregate(stratumId ~ rowId, matchedPop, length)

  saveRDS(matchedPop, file.path(outputFolder, "matchedPop.rds"))

  # matchedPop <- readRDS(file.path(outputFolder, "matchedPop.rds"))
}
