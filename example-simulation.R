# This script recreates a single run from the simulation setup described in the paper.
# The main function is covVimLyl. 
#
# covVimLyl(time,
#           status,
#           treatment,
#           confounders,
#           l,
#           nuisance,
#           tauLearner,
#           folds)
#           
# time: vector of event times
# status: vector of event types (1: event of interest, 2: competing event, 0: censoring)
# treatment: vector of treatment indicators (factor vector with 2 levels)
# confounders: data.frame of confounders. 
# l: vector of indices of confounders for which variable importance is to be estimated.
# nuisance: list of speficications for nuisance estimators on the following form:
#                                 nuisance = list(surv = list(method, ...),
#                                                 comp = list(method, ...),
#                                                 cens = list(method, ...),
#                                                 prop = list(method, ...),
#                                                 pseudoReg = list(method, ...))
#                           Here, method can be different for each of the nuisance parameters (as well as ...). 
#                           Available methods for surv, comp and cens are: "coxph", "gam", "rfsrc", "ranger", "rfsrcTune", "glmnet".
#                             - see functions/survML.R for details.
#                           Available methods for prop are: "glm", "gam", "ranger", "rfsrc".
#                             - see functions/propML.R for details.
#                           Available methods for pseudoReg are: "gam", "ranger", "rfsrc".
# tauLearner: character ("S" or "T"). Defines meta-learner for tau-estimation.
# folds: integer. Defines the number of folds used for cross-fitting.

library(survival)
library(mgcv)
library(riskRegression)

# set path equal to path of current script
path <- getwd()

f1_tmp <- list.files(paste0(path, "/functions"), full.names = TRUE)
f1 <- f1_tmp[grepl(".R", f1_tmp)]
lapply(f1, source)


# Cox-Weibull distribution function and its invers
inv_weibull <- function(u, theta, scale, shape){
  t <- (-log(u) / (scale * exp(theta))) ^ (1 / shape)
  return(t)
}  

weibull <- function(t, theta, scale, shape){
  p <- exp(-scale*exp(theta)*t^shape)
  return(p)
}


expit <- function(x) exp(x)/(1+exp(x))



# Set parameters 
times <- 30                      # event horizon
shape <- 2                      # shape parameter for weibull distributed event times 
scale <- 0.0025                    # scale parameter for weibull distributed event times 
n <- 500                       # number of observations

x1 <- runif(n, min = -1, max = 1)
x2 <- runif(n, min = -1, max = 1)
x3 <- runif(n, min = -1, max = 1)
x4 <- runif(n, min = -1, max = 1)
xx <- data.frame(x1,x2,x3,x4)
A <- rbinom(n, 1, expit(0.5*x1 + 0.5*x2))
survTime <- inv_weibull(runif(n), theta = -x1 - x2 - 0.2*x3 + A*(-2 + 0.5*x1 - 0.3*x2), scale = scale, shape = shape)
compTime <- inv_weibull(runif(n), theta = -x1 - x2 - 0.2*x3 + A, scale = 0.1*scale, shape = shape)
eventTime <- pmin(survTime, compTime)
censTime <- inv_weibull(runif(n), theta = -0.5*x1, scale = 0.1*scale, shape = shape)
obsTime <- pmin(eventTime, censTime)
delta <- 0*(censTime == obsTime) + 1*(survTime == obsTime) + 2*(compTime == obsTime)
A_fac <- factor(A)
###########################

correctGam <- covVimLyl(time = obsTime, status = delta, treatment = A_fac, confounders = xx,
                         l = 1:4,
                         evaltimes = times,
                         nuisance = list(surv = list(method = "coxph", form = "x1*A+x2*A+x3"),
                                         comp = list(method = "coxph", form = "x1+x2+x3+A"),
                                         cens = list(method = "coxph", form = "x1"),
                                         prop = list(method = "glm"),
                                         pseudoRegTau = list(method = "gam"),
                                         pseudoRegX = list(method = "gam")),
                         tauLearner = "S",
                         folds = 1)

correctGamCF <- covVimLyl(time = obsTime, status = delta, treatment = A_fac, confounders = xx,
                           l = 1:4,
                           evaltimes = times,
                           nuisance = list(surv = list(method = "coxph", form = "x1*A+x2*A+x3"),
                                           comp = list(method = "coxph", form = "x1+x2+x3+A"),
                                           cens = list(method = "coxph", form = "x1"),
                                           prop = list(method = "glm"),
                                           pseudoRegTau = list(method = "gam"),
                                           pseudoRegX = list(method = "gam")),
                           tauLearner = "S",
                           folds = 10)


RFRF <- covVimLyl(time = obsTime, status = delta, treatment = A_fac, confounders = xx,
                   l = 1:4,
                   evaltimes = times,
                   nuisance = list(surv = list(method = "rfsrc"),
                                   comp = list(method = "rfsrc"),
                                   cens = list(method = "rfsrc"),
                                   prop = list(method = "rfsrc"),
                                   pseudoRegTau = list(method = "rfsrc"),
                                   pseudoRegX = list(method = "rfsrc")),
                   tauLearner = "S",
                   folds = 1)

RFRFCF <- covVimLyl(time = obsTime, status = delta, treatment = A_fac, confounders = xx,
                     l = 1:4,
                     evaltimes = times,
                     nuisance = list(surv = list(method = "rfsrc"),
                                     comp = list(method = "rfsrc"),
                                     cens = list(method = "rfsrc"),
                                     prop = list(method = "rfsrc"),
                                     pseudoRegTau = list(method = "rfsrc"),
                                     pseudoRegX = list(method = "rfsrc")),
                     tauLearner = "S",
                     folds = 10)


# ATE estimates corresponding to section 5.1.
resAte <- data.frame(method = c("corGam", "corGamCF", "RFRF", "RFRFCF"),
                     rbind(correctGam$ate, correctGamCF$ate,
                           RFRF$ate, RFRFCF$ate))

# Best partially linear projection estimates corresponding to section 5.2.
resMod <- data.frame(method = c(rep("corGam", 4), rep("corGamCF", 4),
                                rep("RFRF", 4), rep("RFRFCF", 4)) ,
                     rbind(correctGam$modifiers, correctGamCF$modifiers,
                           RFRF$modifiers, RFRFCF$modifiers))


