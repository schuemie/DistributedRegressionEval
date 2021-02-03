###
#
# Main function: outputs meta and mega-analysis estimates given a data frame.
# The options are for specifying the columns for ID, Group, and Response. The
# Group column will be used to separate the dataframe into a list with each
# element being data from one group.
#
# Output:
#   x.size: number of subjects in each site
#   meta.list: list of meta-analysis estimates in format c(intercept, effects)
#   meta.a, mega.a: intercept estimates
#   meta.b, mega.b: effect size estimates
#   ...var: corresponding var-covar matrix
#   iter: first K elements are iterations for each site, K+1st element is for mega
#   n.a, n.b: count of sites where Newton's method converges
#
###


meta.mega <- function(df, id = "id", group = "group", response = "response") {
  # set parameters
  n = ncol(df) - 3 # number of variables
  KK <-  length(unique(getElement(df, group)))

  x <- df[, !(names(df) %in% c(id, group, response))] # subset to desired variables
  xnames <- names(x)  #get covariates names
  x <- split(x, as.factor(getElement(df, group)), drop = TRUE) # split into list
  x <- lapply(x, function(a) as.matrix(a))

  x.size <- unlist(lapply(x, function(a) nrow(a))) # get number of subjects in each site

  y <- df[, (names(df) %in% c(response))]
  y <- split(y, as.factor(getElement(df, group)), drop = TRUE)
  y <- lapply(y, function(a) as.matrix(a))



  # Meta (inverse-variance weighted estimator)

  # regression formula
  myformula <- as.formula( paste(response,
                                 paste(xnames,collapse = "+"),
                                 sep="~"))

  meta.tmp <- rep(0,length(xnames)+1)
  var.tmp <- array(0,dim=c(length(xnames)+1,length(xnames)+1))
  for (k in 1:KK) {
    tmpdat <- data.frame(y[[k]], x[[k]])
    colnames(tmpdat)[1] <- response
    tmpglm <- glm(myformula,data = tmpdat, family = "binomial")
    param <- coef(summary(tmpglm))
    parvcov <- vcov(tmpglm)
    std <- param[,2]
    meta.list <- param[,1]

    if(anyNA(std)) {
      print(sprintf("Study %s not included in meta-analysis estimate", k))
    } else {
      meta.tmp <- meta.tmp + meta.list %*% solve(parvcov)
      var.tmp <- var.tmp + solve(parvcov)
    }
  }

  meta.b <-  c(meta.tmp%*% solve(var.tmp))
  meta.se <- c(diag(sqrt(solve(var.tmp))))

  meta.p <- sapply(meta.b/meta.se, pnorm, lower.tail = FALSE)*2

  meta.out <- data.frame(Estimate = meta.b,
                         Std.Error = meta.se, Pvalue = meta.p,
                         ci.lwr = meta.b- 1.96*meta.se,
                         cl.upr = meta.b+ 1.96*meta.se)


  # Mega
  megaglm <- glm(myformula,data = df, family = "binomial")
  param <- coef(summary(megaglm))
  parvcov <- vcov(megaglm)
  mega.b <- param[,1]
  mega.se <- diag(sqrt(vcov(megaglm)))

  mega.out <- data.frame(Estimate = param[,1],
                         Std.Error = param[,2], Pvalue = param[,4],
                         ci.lwr = mega.b- 1.96*mega.se,
                         cl.upr = mega.b+ 1.96*mega.se)

  return(list(
    x.size = x.size,
    meta.out = meta.out,
    mega.out = mega.out
  ))
}

evaluateMega <- function(studyFolder, skipJmdc = FALSE, splitCcae = FALSE) {
  studyFolder <- "r:/DistributedRegressionEval"
  # AMI ---------------------------------------------------------
  outcomeId = 5
  pathToCsv <- system.file("settings", "CohortsToCreate.csv", package = "DistributedRegressionEval")
  cohortsToCreate <- read.csv(pathToCsv)
  outcomeName <- cohortsToCreate$name[cohortsToCreate$cohortId == outcomeId]

  writeLines(paste("Evaluating outcome", outcomeName))
  data <- readRDS(file.path(studyFolder, sprintf("data_o%s.rds", outcomeId)))
  data$age_in_years <- NULL
  data$`gender_=_FEMALE` <- NULL
  data$time <- NULL

  data$id <- seq(1, nrow(data))
  colnames(data)[colnames(data) == "database"] <- "group"
  colnames(data)[colnames(data) == "y"] <- "response"

  results <- meta.mega(df = data)
  saveRDS(results, file.path(studyFolder, "AMI_meta_mega_results.rds"))

    # Heart failure ---------------------------
  outcomeId = 4
  pathToCsv <- system.file("settings", "CohortsToCreate.csv", package = "DistributedRegressionEval")
  cohortsToCreate <- read.csv(pathToCsv)
  outcomeName <- cohortsToCreate$name[cohortsToCreate$cohortId == outcomeId]

  writeLines(paste("Evaluating outcome", outcomeName))
  data <- readRDS(file.path(studyFolder, sprintf("data_o%s.rds", outcomeId)))
  data$time <- NULL

  data$id <- seq(1, nrow(data))
  colnames(data)[colnames(data) == "gender_=_FEMALE"] <- "Female"
  colnames(data)[colnames(data) == "database"] <- "group"
  colnames(data)[colnames(data) == "y"] <- "response"

  results <- meta.mega(df = data)
  saveRDS(results, file.path(studyFolder, "HF_meta_mega_results.rds"))
}
