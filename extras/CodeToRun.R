library(DistributedRegressionEval)
options(fftempdir = "r:/FFtemp")

maxCores <- 30
studyFolder <- "r:/DistributedRegressionEval"
dbms <- "pdw"
user <- NULL
pw <- NULL
server <- Sys.getenv("PDW_SERVER")
port <- Sys.getenv("PDW_PORT")
oracleTempSchema <- NULL
connectionDetails <- DatabaseConnector::createConnectionDetails(dbms = dbms,
                                                                server = server,
                                                                user = user,
                                                                password = pw,
                                                                port = port)

# CCAE settings ----------------------------------------------------------------
cdmDatabaseSchema <- "cdm_truven_ccae_v778.dbo"
cohortDatabaseSchema <- "scratch.dbo"
cohortTable <- "mschuemi_informed_priors_ccae"
databaseName <- "CCAE"
outputFolder <- file.path(studyFolder, "ccae")

# MDCD settings ----------------------------------------------------------------
cdmDatabaseSchema <- "cdm_truven_mdcd_v780.dbo"
cohortDatabaseSchema <- "scratch.dbo"
cohortTable <- "mschuemi_informed_priors_mdcd"
databaseName <- "MDCD"
outputFolder <- file.path(studyFolder, "mdcd")

# MDCR settings ----------------------------------------------------------------
cdmDatabaseSchema <- "cdm_truven_mdcr_v779.dbo"
cohortDatabaseSchema <- "scratch.dbo"
cohortTable <- "mschuemi_informed_priors_mdcr"
databaseName <- "MDCR"
outputFolder <- file.path(studyFolder, "mdcr")

# Optum settings ----------------------------------------------------------------
cdmDatabaseSchema <- "cdm_optum_extended_dod_v774.dbo"
cohortDatabaseSchema <- "scratch.dbo"
cohortTable <- "mschuemi_informed_priors_optum"
databaseName <- "Optum"
outputFolder <- file.path(studyFolder, "optum")

# JMDC settings ----------------------------------------------------------------
cdmDatabaseSchema <- "cdm_jmdc_v773.dbo"
cohortDatabaseSchema <- "scratch.dbo"
cohortTable <- "mschuemi_informed_priors_jmdc"
databaseName <- "Jmdc"
outputFolder <- file.path(studyFolder, "Jmdc")



createCohorts(connectionDetails = connectionDetails,
              cdmDatabaseSchema = cdmDatabaseSchema,
              cohortDatabaseSchema = cohortDatabaseSchema,
              cohortTable = cohortTable,
              oracleTempSchema = oracleTempSchema,
              outputFolder = outputFolder)

createPlpObjects(connectionDetails = connectionDetails,
                 cdmDatabaseSchema = cdmDatabaseSchema,
                 cohortDatabaseSchema = cohortDatabaseSchema,
                 cohortTable = cohortTable,
                 oracleTempSchema = oracleTempSchema,
                 outputFolder = outputFolder)



