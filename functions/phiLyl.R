phiLyl <- function(surv, comp, cens, treat, newtimes, newstatus, newX, newA, predtimes){
  
  nobs <- length(newtimes)
  time.jump1 <- unique(sort(surv$time[surv$status == 1]))
  time.jump2 <- unique(sort(comp$time[comp$status == 1]))
  time.event <- unique(sort(c(time.jump1, time.jump2)))
  tauInd1 <- max(which(time.jump1 <= predtimes))
  tauInd2 <- max(which(time.jump2 <= predtimes))
  tauIndEvent <- max(which(time.event <= predtimes))
  time.zero <- c(0, time.event)
  tauIndZero <- max(which(time.zero <= predtimes))
  pmhat <- predict(treat, newX = newX)
  # pmhat[pmhat == 1] <- 1 - 1e-06  # "incomment" to prevent crash due to propensity positivity violation. 
  # pmhat[pmhat == 0] <- 1e-06 # "incomment" to prevent crash due to propensity positivity violation. 
  
  # M1 integral
  
  Lambda1.te <- predict(surv, newtimes = time.zero, newX = newX, newA = newA)$cumhaz
  Lambda2.te <- predict(comp, newtimes = time.zero, newX = newX, newA = newA)$cumhaz
  
  sHats.te <- exp(-Lambda1.te - Lambda2.te)
  gHats.te <- predict(cens, newtimes = time.zero, newX = newX, newA = newA)$surv
  
  dL1.te <- t(diff(t(Lambda1.te)))
  dL2.te <- t(diff(t(Lambda2.te)))
  
  dF1.te <- sHats.te[,-length(time.zero)]*dL1.te
  dF2.te <- sHats.te[,-length(time.zero)]*dL2.te
  F1.te <- cbind(0, rowCumSum(dF1.te))
  F2.te <- cbind(0, rowCumSum(dF2.te))
  
  dN1 <- sapply(predtimes, function(t0){
    1*(newtimes <= t0 & newstatus == 1)
  })
  
  dN2 <- sapply(predtimes, function(t0){
    1*(newtimes <= t0 & newstatus == 2)
  })
  
  mgEvalIndex <- do.call(cbind,lapply(predtimes, function(tau){
    prodlim::sindex(jump.times = time.zero, eval.times = pmin(newtimes,tau))
  }))
  
  if(!(predtimes %in% time.zero)){
    time.zero[tauIndZero + 1] <- predtimes
  }
  
  H11 <- t((predtimes - time.zero[1:tauIndZero])*t(1 + F1.te[, 1:tauIndZero]/sHats.te[, 1:tauIndZero])) -
    (1/sHats.te[, 1:tauIndZero]) * rowCumSum(t(diff(time.zero[1:(tauIndZero + 1)]) * t(F1.te[, 1:tauIndZero]))[, tauIndZero:1])[, tauIndZero:1]
  
  H21 <- t((predtimes - time.zero[1:tauIndZero])*t(F1.te[, 1:tauIndZero]/sHats.te[, 1:tauIndZero])) -
    (1/sHats.te[, 1:tauIndZero]) * rowCumSum(t(diff(time.zero[1:(tauIndZero + 1)]) * t(F1.te[, 1:tauIndZero]))[, tauIndZero:1])[, tauIndZero:1]
  
  H22 <- t((predtimes - time.zero[1:tauIndZero])*t(1 + F2.te[, 1:tauIndZero]/sHats.te[, 1:tauIndZero])) -
    (1/sHats.te[, 1:tauIndZero]) * rowCumSum(t(diff(time.zero[1:(tauIndZero + 1)]) * t(F2.te[, 1:tauIndZero]))[, tauIndZero:1])[, tauIndZero:1]
  
  H12 <- t((predtimes - time.zero[1:tauIndZero])*t(F2.te[, 1:tauIndZero]/sHats.te[, 1:tauIndZero])) -
    (1/sHats.te[, 1:tauIndZero]) * rowCumSum(t(diff(time.zero[1:(tauIndZero + 1)]) * t(F2.te[, 1:tauIndZero]))[, tauIndZero:1])[, tauIndZero:1]
  
  L1.1_tmp <- rowCumSum((H11[, 1:tauIndZero]/gHats.te[, 1:tauIndZero])*dL1.te[, 1:tauIndZero])
  L1.2_tmp <- rowCumSum((H12[, 1:tauIndZero]/gHats.te[, 1:tauIndZero])*dL1.te[, 1:tauIndZero])
  
  L2.1_tmp <- rowCumSum((H21[, 1:tauIndZero]/gHats.te[, 1:tauIndZero])*dL2.te[, 1:tauIndZero])
  L2.2_tmp <- rowCumSum((H22[, 1:tauIndZero]/gHats.te[, 1:tauIndZero])*dL2.te[, 1:tauIndZero])
  
  L1.1 <- do.call(cbind, lapply(seq_along(predtimes), function(x){
    tmp <- matrix(0, nrow = nobs, ncol = 1)
    tmp[mgEvalIndex[, x] > 0] <- L1.1_tmp[cbind(1:nobs, mgEvalIndex[, x])]
    tmp
  }))
  
  L1.2 <- do.call(cbind, lapply(seq_along(predtimes), function(x){
    tmp <- matrix(0, nrow = nobs, ncol = 1)
    tmp[mgEvalIndex[, x] > 0] <- L1.2_tmp[cbind(1:nobs, mgEvalIndex[, x])]
    tmp
  }))
  
  L2.2 <- do.call(cbind, lapply(seq_along(predtimes), function(x){
    tmp <- matrix(0, nrow = nobs, ncol = 1)
    tmp[mgEvalIndex[, x] > 0] <- L2.2_tmp[cbind(1:nobs, mgEvalIndex[, x])]
    tmp
  }))
  
  L2.1 <- do.call(cbind, lapply(seq_along(predtimes), function(x){
    tmp <- matrix(0, nrow = nobs, ncol = 1)
    tmp[mgEvalIndex[, x] > 0] <- L2.1_tmp[cbind(1:nobs, mgEvalIndex[, x])]
    tmp
  }))
  
  N1.1 <- do.call(cbind, lapply(seq_along(predtimes), function(x){
    tmp <- matrix(0, nrow = nobs, ncol = 1)
    tmp[mgEvalIndex[, x] > 0] <- (H11[, 1:tauIndZero]/gHats.te[, 1:tauIndZero])[cbind(1:nobs, mgEvalIndex[, x])] * dN1[, x][mgEvalIndex[, x] > 0]
    tmp
  }))
  
  N1.2 <- do.call(cbind, lapply(seq_along(predtimes), function(x){
    tmp <- matrix(0, nrow = nobs, ncol = 1)
    tmp[mgEvalIndex[, x] > 0] <- (H12[, 1:tauIndZero]/gHats.te[, 1:tauIndZero])[cbind(1:nobs, mgEvalIndex[, x])] * dN1[, x][mgEvalIndex[, x] > 0]
    tmp
  }))
  
  N2.2 <- do.call(cbind, lapply(seq_along(predtimes), function(x){
    tmp <- matrix(0, nrow = nobs, ncol = 1)
    tmp[mgEvalIndex[, x] > 0] <- (H22[, 1:tauIndZero]/gHats.te[, 1:tauIndZero])[cbind(1:nobs, mgEvalIndex[, x])] * dN2[, x][mgEvalIndex[, x] > 0]
    tmp
  }))
  
  N2.1 <- do.call(cbind, lapply(seq_along(predtimes), function(x){
    tmp <- matrix(0, nrow = nobs, ncol = 1)
    tmp[mgEvalIndex[, x] > 0] <- (H21[, 1:tauIndZero]/gHats.te[, 1:tauIndZero])[cbind(1:nobs, mgEvalIndex[, x])] * dN2[, x][mgEvalIndex[, x] > 0]
    tmp
  }))
  
  M1.1 <- N1.1 - L1.1
  M2.1 <- N2.1 - L2.1
  
  M1.2 <- N1.2 - L1.2
  M2.2 <- N2.2 - L2.2
  
  # calculate counterfactual L_j's
  times.tmp <- time.event
  if(!(predtimes %in% time.event)){
    times.tmp[tauIndEvent + 1] <- predtimes
  }
  times.tmp <- times.tmp[1:(tauIndEvent + 1)]
  
  Lambda1.te.a1 <- predict(surv, newtimes = time.event, newX = newX, newA = 1)$cumhaz
  Lambda1.te.a0 <- predict(surv, newtimes = time.event, newX = newX, newA = 0)$cumhaz
  Lambda2.te.a1 <- predict(comp, newtimes = time.event, newX = newX, newA = 1)$cumhaz
  Lambda2.te.a0 <- predict(comp, newtimes = time.event, newX = newX, newA = 0)$cumhaz
  
  dL1.te.a1 <- cbind(Lambda1.te.a1[, 1], t(diff(t(Lambda1.te.a1))))
  dL2.te.a1 <- cbind(Lambda2.te.a1[, 1], t(diff(t(Lambda2.te.a1))))
  dL1.te.a0 <- cbind(Lambda1.te.a0[, 1], t(diff(t(Lambda1.te.a0))))
  dL2.te.a0 <- cbind(Lambda2.te.a0[, 1], t(diff(t(Lambda2.te.a0))))
  
  F1.te.a1 <- rowCumSum(exp(-(Lambda1.te.a1 + Lambda2.te.a1)) * dL1.te.a1)
  F1.te.a0 <- rowCumSum(exp(-(Lambda1.te.a0 + Lambda2.te.a0)) * dL1.te.a0)
  
  F2.te.a1 <- rowCumSum(exp(-(Lambda1.te.a1 + Lambda2.te.a1)) * dL2.te.a1)
  F2.te.a0 <- rowCumSum(exp(-(Lambda1.te.a0 + Lambda2.te.a0)) * dL2.te.a0)
  
  L_1.1 <- F1.te.a1[,1:tauIndEvent] %*% diff(times.tmp)
  L_1.0 <- F1.te.a0[,1:tauIndEvent] %*% diff(times.tmp)
  
  L_2.1 <- F2.te.a1[,1:tauIndEvent] %*% diff(times.tmp)
  L_2.0 <- F2.te.a0[,1:tauIndEvent] %*% diff(times.tmp)
  
  # CF phi's
  phi1.1 <- L_1.1 + (1*(newA == 1)/pmhat)*(M1.1 + M2.1)
  phi1.0 <- L_1.0 + (1*(newA == 0)/(1 - pmhat))*(M1.1 + M2.1)
  
  phi2.1 <- L_2.1 + (1*(newA == 1)/pmhat)*(M1.2 + M2.2)
  phi2.0 <- L_2.0 + (1*(newA == 0)/(1 - pmhat))*(M1.2 + M2.2)
  
  phi1 <- phi1.1 - phi1.0
  phi2 <- phi2.1 - phi2.0
  
 
  list(phi1 = phi1, phi2 = phi2, phi1.1 = phi1.1, phi1.0 = phi1.0, phi2.1 = phi2.1, phi2.0 = phi2.0,
       L1 = L_1.1 - L_1.0, L2 = L_2.1 - L_2.0, L_1.1 = L_1.1, L_1.0 = L_1.0, L_2.1 = L_2.1, L_2.0 = L_2.0)
}
























