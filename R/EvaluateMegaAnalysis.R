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
  astart <- rep(0, KK)
  bstart <- matrix(0, nrow = KK, ncol = n, byrow = TRUE)
  iter <- rep(0, KK+1)

  x <- df[, !(names(df) %in% c(id, group, response))] # subset to desired variables
  x <- split(x, as.factor(getElement(df, group)), drop = TRUE) # split into list
  x <- lapply(x, function(a) as.matrix(a))

  x.size <- unlist(lapply(x, function(a) nrow(a))) # get number of subjects in each site

  y <- df[, (names(df) %in% c(response))]
  y <- split(y, as.factor(getElement(df, group)), drop = TRUE)
  y <- lapply(y, function(a) as.matrix(a))

  # Meta (inverse-variance weighted estimator)
  tmp.ab <- rep(0, n+1)
  tmp.hess <- array(0, dim = c(n+1, n+1))
  n.a <- 0
  n.b <- 0
  meta.list <- list()

  for (k in c(1:KK)) {
    param <- newton.logistic(x[[k]], y[[k]], astart[k], bstart[k,], return.fail = TRUE)
    a.est <- param[[1]]
    b.est <- param[[2]]
    iter[k] <- param[[3]]
    hess <- param[[4]]
    meta.list[[k]] <- c(a.est, b.est)

    if(anyNA(hess)) {print(sprintf("Study %s not included in meta-analysis estimate", k))} else {
      tmp.ab <- tmp.ab + hess %*% c(a.est, b.est)
      tmp.hess <- tmp.hess + hess

      if (!param[[5]]) {
        n.a = n.a + 1
        n.b = n.b + 1
      } else {print(sprintf("Hessian not invertible for study %s", k))}
    }
  }

  meta.var <- try(chol2inv(chol(tmp.hess)))

  if(class(meta.var) == "matrix"){
    meta.ab <- meta.var %*% tmp.ab

    meta.a <- meta.ab[1]
    meta.b <- meta.ab[2:(n+1)]
    meta.a.var <- meta.var[1,1]
    meta.b.var <- meta.var[2:(n+1), 2:(n+1)]
  } else {
    meta.a <- NA
    meta.b <- NA
    meta.a.var <- NA
    meta.b.var <- NA
  }

  # Mega
  param <- glore.list(x, y, astart[1], bstart[1,], return.fail = TRUE)
  mega.a <- param[[1]]
  mega.b <- param[[2]]
  mega.a.var <- param[[3]]
  mega.b.var <- param[[4]]
  iter[KK + 1] <- param[[5]] # keeps track of total iterations

  return(list(
    x.size = x.size,
    meta.list = meta.list,
    meta.a = meta.a, mega.a = mega.a,
    meta.b = meta.b, mega.b = mega.b,
    meta.a.var = meta.a.var, mega.a.var = mega.a.var,
    meta.b.var = meta.b.var, mega.b.var = mega.b.var,
    iter = iter,
    n.a = n.a, n.b = n.b
  ))
}

pr <- function(x, a, b) {
  c(exp(a + x%*%b) / (1+exp(a + x%*%b)))
}

glore.list <- function(x, y, a0, b0, niter = 100, nstudy, n = ncol(x[[1]]), return.fail = FALSE, prec = 10^-5, step = 1, damped = FALSE, damping = 0.5, epsilon = 10^-4, debug = FALSE) { # for all studies, a0 and b0 are starting values, niter is max number of iterations
  fail = FALSE # boolean for tracking failure
  step.orig <- step

  for (iter in 1:niter) {
    total.iter <- iter

    grad.sum = rep(0, n + 1)
    hess.sum = matrix(rep(0, (n+1)^2), nrow = n + 1, ncol = n + 1)

    x.bind <- x[[1]] # for computing logl in the damping
    y.bind <- y[[1]]

    for (k in 1:nstudy) {
      x.fix <- cbind(rep(1, nrow(x[[k]])), x[[k]])
      p <- pr(x[[k]], a0, b0)
      pp <- matrix(rep(p * (1-p), n + 1), n + 1, nrow(x.fix), byrow = TRUE)

      grad.sum = grad.sum + c(t(x.fix) %*% (y[[k]] - p))
      hess.sum = hess.sum + (t(x.fix) * pp) %*% x.fix
      # hess.sum = hess.sum + t(x.fix) %*% diag(p * (1-p)) %*% x.fix

      if (k >= 2) {
        x.bind <- rbind(x.bind, x[[k]])
        y.bind <- c(y.bind, y[[k]])
      }
    }

    bf <- try(chol2inv(chol(hess.sum)))

    if(class(bf) == "matrix"){
      newton.step = bf %*% grad.sum
      step <- step.orig
      if (damped) {
        logl.post <- logl(c(a0, b0) + step * newton.step, x.bind, y.bind)
        logl.pre <- logl(c(a0, b0), x.bind, y.bind)
        damp.step <- epsilon * step * t(grad.sum) %*% newton.step

        if (!is.na(logl.post) & !is.na(logl.pre + damp.step)) {
          while (logl.post < logl.pre + damp.step) {
            step = damping * step
            # message(step) # for testing

            logl.post <- logl(c(a0, b0) + step * newton.step, x.bind, y.bind)
            logl.pre <- logl(c(a0, b0), x.bind, y.bind)
            damp.step <- epsilon * step * t(grad.sum) %*% newton.step

            if (is.na(logl.post) | is.na(logl.pre + damp.step) | step < 10^-50) {break}
          }
        }
      }

      if (debug) {
        message(sprintf("%s: %.2f, %.2f, %.2f, %.2f, %.2f", iter, a0, b0, newton.step[1], newton.step[2], logl(c(a0, b0), x.bind, y.bind)))
      }

      a1 = a0 + step * newton.step[1]
      b1 = b0 + step * newton.step[2:(n+1)]

      if (return.fail) {
        grad.sum1 = rep(0, n + 1)
        hess.sum1 = matrix(rep(0, (n+1)^2), nrow = n + 1, ncol = n + 1, byrow = TRUE)

        for (k in 1:nstudy) {
          x.fix <- cbind(rep(1, nrow(x[[k]])), x[[k]])
          p <- pr(x[[k]], a1, b1)
          pp <- matrix(rep(p * (1-p), n + 1), n + 1, nrow(x.fix))

          grad.sum1 = grad.sum1 + c(t(x.fix) %*% (y[[k]] - p))
          hess.sum1 = hess.sum1 + (t(x.fix) * pp) %*% x.fix
          # hess.sum1 = hess.sum1 + t(x.fix) %*% diag(p * (1-p)) %*% x.fix

          if (k >= 2) {
            x.bind <- rbind(x.bind, x[[k]])
            y.bind <- c(y.bind, y[[k]])
          }
        }

        if (is.null(tryCatch(chol2inv(chol(hess.sum1)), error = function(e) NULL))) {
          if (debug) {
            message("Returning last iteration before failure")
          }
          fail = TRUE

          a.var <- bf[1,1]
          b.var <- bf[2:(n+1), 2:(n+1)]
          return(list(c(a0), c(b0), a.var, b.var, total.iter-1, fail))
        }
      }

    } else {
      a1 <- NA
      b1 <- NA
      a.var <- NA
      b.var <- NA
    }

    if (is.na(a1) || is.na(b1)) {
      if (return.fail) {
        fail = TRUE
        return(list(c(a0), c(b0), a.var, b.var, total.iter-1,fail))
      } else {
        return(list(c(a1), c(b1), a.var, b.var, total.iter,fail))
      }

    }
    else if (abs(a1 - a0) > prec && abs(b1 - b0) > prec) {
      a0 <- a1
      b0 <- b1
    }
    else {
      grad.sum1 = rep(0, n + 1)
      hess.sum1 = matrix(rep(0, (n+1)^2), nrow = n + 1, ncol = n + 1)

      for (k in 1:nstudy) {
        x.fix <- cbind(rep(1, nrow(x[[k]])), x[[k]])
        p <- pr(x[[k]], a1, b1)
        pp <- matrix(rep(p * (1-p), n + 1), n + 1, nrow(x.fix))

        grad.sum1 = grad.sum1 + c(t(x.fix) %*% (y[[k]] - p))
        hess.sum1 = hess.sum1 + (t(x.fix) * pp) %*% x.fix

        if (k >= 2) {
          x.bind <- rbind(x.bind, x[[k]])
          y.bind <- c(y.bind, y[[k]])
        }
      }

      bf1 <- try(chol2inv(chol(hess.sum1)))

      a.var <- bf1[1,1]
      b.var <- bf1[2:(n+1), 2:(n+1)]

      return(list(c(a1), c(b1), a.var, b.var, total.iter,fail))
    }
  }

  a.var <- bf[1,1]
  b.var <- bf[2:(n+1), 2:(n+1)]
  return(list(c(a1), c(b1), a.var, b.var, total.iter,fail))
  print("GLORE did not converge")

}

newton.logistic <- function(x, y, a0, b0, niter = 100, return.fail = FALSE, prec = 1e-6, step = 1, damped = FALSE, damping = 0.5, epsilon = 10^-4, debug = FALSE) { # for single study, a0 and b0 are starting values
  fail = FALSE # boolean for tracking failure
  step.orig <- step

  x.fix <- cbind(rep(1, nrow(x)), x)

  for (iter in 1:niter) {
    total.iter <- iter

    p <- pr(x, a0, b0)
    pp <- matrix(rep(p * (1-p), n + 1), n + 1, nrow(x), byrow = TRUE) # faster way to implement t(x.fix) %*% diag(p * (1-p))
    # hess <- t(x.fix) %*% diag(p * (1-p)) %*% x.fix
    hess <- (t(x.fix) * pp) %*% x.fix
    bf <- try(chol2inv(chol(hess)))
    bg <- c(t(x.fix) %*% (y - p))

    # if (debug) {
    #   message(logl(c(a0, b0), x, y))
    #   message(a0)
    #   message(b0)
    #   message(bf)
    # }

    if(class(bf) == "matrix"){
      # implements damped Newton's method: while f(x + step * inv.hess%*%grad) > f(x) + alpha * step * t(grad)%*%inv.hess%*%grad, step = damping * step
      newton.step <- bf %*% bg
      step <- step.orig
      if (damped) {
        logl.post <- logl(c(a0, b0) + step * newton.step, x, y)
        logl.pre <- logl(c(a0, b0), x, y)
        damp.step <- epsilon * step * t(bg) %*% c(newton.step)

        if (!is.na(logl.post) & !is.na(logl.pre + damp.step)) {
          while(logl.post < logl.pre + damp.step) {
            step = damping * step
            # message(step) # for testing

            logl.post <- logl(c(a0, b0) + step * newton.step, x, y)
            logl.pre <- logl(c(a0, b0), x, y)
            damp.step <- epsilon * step * t(bg) %*% c(newton.step)

            if (is.na(logl.post) | is.na(logl.pre + damp.step) | step < 10^-50) {break}
          }
        }
      }

      if (debug) {
        message(sprintf("%s: %.2f, %.2f, %.2f, %.2f, %.2f", iter, a0, b0, newton.step[1], newton.step[2], logl(c(a0, b0), x, y)))
      }

      a1 = a0 + step * newton.step[1]
      b1 = b0 + step * newton.step[2:(n+1)]

      if (return.fail) {
        p1 <- pr(x, a1, b1)
        pp1 <- matrix(rep(p1 * (1-p1), n + 1), n + 1, nrow(x), byrow = TRUE) # faster way to implement t(x.fix) %*% diag(p * (1-p))
        hess1 <- (t(x.fix) * pp1) %*% x.fix
        # hess1 <- t(x.fix) %*% diag(p1 * (1-p1)) %*% x.fix
        bf1 <- try(chol2inv(chol(hess1)))

        if (class(bf1) != "matrix") {
          if (debug) {
            message("Returning last iteration before failure")
          }
          fail = TRUE

          return(list(c(a0), c(b0), total.iter-1, hess, fail))
        }
      }

    } else {
      a1 <- NA
      b1 <- NA
      hess <- NA
    }

    if (is.na(a1) || is.na(b1)) {
      if (return.fail) {
        fail = TRUE
        return(list(c(a0),c(b0), total.iter-1, hess, fail))
      } else {
        return(list(c(a1),c(b1), total.iter, hess,fail))
      }

    } else if (abs(a1 - a0) > prec && abs(b1 - b0) > prec) {
      a0 <- a1
      b0 <- b1
    } else {
      p1 <- pr(x, a1, b1)
      pp1 <- matrix(rep(p1 * (1-p1), n + 1), n + 1, nrow(x), byrow = TRUE)

      hess1 <- (t(x.fix) * pp1) %*% x.fix
      return(list(c(a1),c(b1), total.iter, hess1,fail))
    }
  }

  return(list(c(a1),c(b1), total.iter, hess, fail))
  print("Newton logistic did not converge")
}


evaluateMega <- function(studyFolder, skipJmdc = FALSE, splitCcae = FALSE) {
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

  dbs <- unique(data$database)
  x <- lapply(dbs, function(db) as.matrix(data[data$database == db, -which(colnames(data) %in% c("y", "database"))]))
  y <- lapply(dbs, function(db) data$y[data$database == db])

  a0 <- 0
  b0 <- rep(0, ncol(x[[1]]))

  out1 <- glore.list(x = x,
                    y = y,
                    a0 = a0,
                    b0 = b0,
                    nstudy = length(x),
                    debug = FALSE)


  pretty1 <- t(c(out1[[1]], out1[[2]]))

  # USing meta-analysis estimate as initial value:
  a0 <- -4.821312871
  b0 <- c(-0.069581107,	0.552218731,	0.863125647,	-0.346934496,	0.61961744,	0.181170585)
  out2 <- glore.list(x = x,
                    y = y,
                    a0 = a0,
                    b0 = b0,
                    nstudy = length(x),
                    debug = FALSE)

  pretty2 <- t(c(out2[[1]], out2[[2]]))
  pretty <- rbind(pretty1, pretty2)
  pretty <- cbind(pretty, iterations = c(out1[[5]], out2[[5]]))
  write.csv(pretty, file.path(studyFolder, sprintf("output_mega_%s.csv", outcomeName)))

  # Heart failure ---------------------------
  outcomeId = 4
  pathToCsv <- system.file("settings", "CohortsToCreate.csv", package = "DistributedRegressionEval")
  cohortsToCreate <- read.csv(pathToCsv)
  outcomeName <- cohortsToCreate$name[cohortsToCreate$cohortId == outcomeId]

  writeLines(paste("Evaluating outcome", outcomeName))
  data <- readRDS(file.path(studyFolder, sprintf("data_o%s.rds", outcomeId)))
  data$time <- NULL

  dbs <- unique(data$database)
  x <- lapply(dbs, function(db) as.matrix(data[data$database == db, -which(colnames(data) %in% c("y", "database"))]))
  y <- lapply(dbs, function(db) data$y[data$database == db])

  a0 <- 0
  b0 <- rep(0, ncol(x[[1]]))

  out1 <- glore.list(x = x,
                    y = y,
                    a0 = a0,
                    b0 = b0,
                    nstudy = length(x),
                    debug = FALSE)


  pretty1 <- t(c(out1[[1]], out1[[2]]))

  # USing meta-analysis estimate as initial value:
  a0 <- -7.925287993
  b0 <- c(-0.12358508,	0.318519908,	0.437623879,	0.380657722,	0.389520888,	0.498572917,	0.651213266,	0.375646035,	0.579801763,	0.457754935,	0.233975656,	0.205752966,	0.058501978)
  out2 <- glore.list(x = x,
                    y = y,
                    a0 = a0,
                    b0 = b0,
                    nstudy = length(x),
                    debug = FALSE)
  pretty2 <- t(c(out2[[1]], out2[[2]]))
  pretty <- rbind(pretty1, pretty2)
  pretty <- cbind(pretty, iterations = c(out1[[5]], out2[[5]]))
  write.csv(pretty, file.path(studyFolder, sprintf("output_mega_%s.csv", outcomeName)))

}
