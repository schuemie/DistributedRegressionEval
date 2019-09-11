library(DistributedRegressionEval)
library(survival)
# install.packages("c:/temp/ODACO_2.0.tar.gz", repos = NULL, type = "source")

options(fftempdir = "r:/FFtemp")

maxCores <- 30
studyFolder <- "r:/DistributedCoxEval"
oracleTempSchema <- NULL
user <- Sys.getenv("redShiftUser")
pw <- Sys.getenv("redShiftPassword")

# JMDC settings (RedShift) -------------------------------------------------
cdmDatabaseSchema <- "cdm"
cohortDatabaseSchema <- "scratch_mschuemi"
cohortTable <- "informed_priors_jmdc"
databaseName <- "Jmdc"
outputFolder <- file.path(studyFolder, "Jmdc")
connectionString <- Sys.getenv("jmdcRedShiftConnectionString")
connectionDetails <- DatabaseConnector::createConnectionDetails(dbms = "redshift",
                                                                connectionString = connectionString,
                                                                user = user,
                                                                password = pw)

# PanTher settings (RedShift) -------------------------------------------------
cdmDatabaseSchema <- "cdm"
cohortDatabaseSchema <- "scratch_mschuemi"
cohortTable <- "informed_priors"
databaseName <- "panther"
outputFolder <- file.path(studyFolder, "panther")
connectionString <- Sys.getenv("pantherRedShiftConnectionString")
connectionDetails <- DatabaseConnector::createConnectionDetails(dbms = "redshift",
                                                                connectionString = connectionString,
                                                                user = user,
                                                                password = pw)

# CCAE settings (RedShift) -------------------------------------------------
cdmDatabaseSchema <- "cdm"
cohortDatabaseSchema <- "scratch_mschuemi"
cohortTable <- "informed_priors"
databaseName <- "ccae"
outputFolder <- file.path(studyFolder, "ccae")
connectionString <- Sys.getenv("ccaeRedShiftConnectionString")
connectionDetails <- DatabaseConnector::createConnectionDetails(dbms = "redshift",
                                                                connectionString = connectionString,
                                                                user = user,
                                                                password = pw)

# MDCD settings (RedShift) -------------------------------------------------
cdmDatabaseSchema <- "cdm"
cohortDatabaseSchema <- "scratch_mschuemi"
cohortTable <- "informed_priors"
databaseName <- "mdcd"
outputFolder <- file.path(studyFolder, "mdcd")
connectionString <- Sys.getenv("mdcdRedShiftConnectionString")
connectionDetails <- DatabaseConnector::createConnectionDetails(dbms = "redshift",
                                                                connectionString = connectionString,
                                                                user = user,
                                                                password = pw)

# Optum settings (RedShift) -------------------------------------------------
cdmDatabaseSchema <- "cdm"
cohortDatabaseSchema <- "scratch_mschuemi"
cohortTable <- "informed_priors"
databaseName <- "optum"
outputFolder <- file.path(studyFolder, "optum")
connectionString <- Sys.getenv("optumDodRedShiftConnectionString")
connectionDetails <- DatabaseConnector::createConnectionDetails(dbms = "redshift",
                                                                connectionString = connectionString,
                                                                user = user,
                                                                password = pw)


createCohorts(connectionDetails = connectionDetails,
              cdmDatabaseSchema = cdmDatabaseSchema,
              cohortDatabaseSchema = cohortDatabaseSchema,
              cohortTable = cohortTable,
              oracleTempSchema = oracleTempSchema,
              outputFolder = outputFolder)

createCoxObjects(connectionDetails = connectionDetails,
                 cdmDatabaseSchema = cdmDatabaseSchema,
                 cohortDatabaseSchema = cohortDatabaseSchema,
                 cohortTable = cohortTable,
                 oracleTempSchema = oracleTempSchema,
                 outputFolder = outputFolder)

combineCoxData(studyFolder = studyFolder, sampleSize = 1e9)

library(survival)
evaluateCox.prop(studyFolder = studyFolder)

