phiLyl_old <- function(surv, comp, cens, treat, newtimes, newstatus, newX, newA, predtimes){
  
  browser()
  
  nobs <- length(newtimes)
  time.jump1 <- unique(sort(surv$time[surv$status == 1]))
  time.jump2 <- unique(sort(comp$time[comp$status == 1]))
  time.event <- unique(sort(c(time.jump1, time.jump2)))
  tauInd1 <- max(which(time.jump1 <= predtimes))
  tauInd2 <- max(which(time.jump2 <= predtimes))
  tauIndEvent <- max(which(time.event <= predtimes))
  
  # pmhat
  pmhat <- predict(treat, newX = newX)
  pmhat[pmhat == 1] <- 1 - 1e-06
  pmhat[pmhat == 0] <- 1e-06
  
  # M1 integral
  
  Lambda1.te <- predict(surv, newtimes = c(0, time.event), newX = newX, newA = newA)$cumhaz
  Lambda2.te <- predict(comp, newtimes = c(0, time.event), newX = newX, newA = newA)$cumhaz
  
  # Lambda1.obs.t1 <- predict(surv, newtimes = time.jump1, newX = newX, newA = newA)$cumhaz
  # Lambda2.obs.t1 <- predict(comp, newtimes = time.jump1, newX = newX, newA = newA)$cumhaz
  # 
  # Lambda1.obs.t2 <- predict(surv, newtimes = time.jump2, newX = newX, newA = newA)$cumhaz
  # Lambda2.obs.t2 <- predict(comp, newtimes = time.jump2, newX = newX, newA = newA)$cumhaz

  sHats.te <- exp(-Lambda1.te - Lambda2.te)
  gHats.te <- predict(cens, newtimes = c(0, time.event), newX = newX, newA = newA)$surv
  
  # sHats.t1 <- exp(-Lambda1.obs.t1 - Lambda2.obs.t1)
  # gHats.t1 <- predict(cens, newtimes = time.jump1, newX = newX, newA = newA)$surv
  # 
  # sHats.t2 <- exp(-Lambda1.obs.t2 - Lambda2.obs.t2)
  # gHats.t2 <- predict(cens, newtimes = time.jump2, newX = newX, newA = newA)$surv
  
  # dL1.te <- cbind(Lambda1.te[, 1], t(diff(t(Lambda1.te))))
  # dL2.te <- cbind(Lambda2.te[, 1], t(diff(t(Lambda2.te))))
  dL1.te <- t(diff(t(Lambda1.te)))
  dL2.te <- t(diff(t(Lambda2.te)))
  
  dF1.te <- sHats.te[,-length(time.event)]*dL1.te
  dF2.te <- sHats.te[,-length(time.event)]*dL2.te
  F1.te <- rowCumSum(dF1.te)
  F2.te <- rowCumSum(dF2.te)
  
  #dL1 <- cbind(Lambda2.obs.t1[, 1], t(diff(t(Lambda1.obs.t1))))
  dN1 <- sapply(predtimes, function(t0){
    1*(newtimes <= t0 & newstatus == 1)
  })
  
  dN2 <- sapply(predtimes, function(t0){
    1*(newtimes <= t0 & newstatus == 2)
  })
  
  mgEvalIndex <- do.call(cbind,lapply(predtimes, function(tau){
    prodlim::sindex(jump.times = time.event, eval.times = pmin(newtimes,tau))
  }))
  
  
  H11 <- do.call(cbind, lapply(c(0, time.event[1:tauIndEvent]), function(x){
    Hij(Fj = F1.te, surv = sHats.te, iequalj = TRUE, time.jump = c(0, time.event), smin = x, smax = predtimes)
  }))
  
  H21 <- do.call(cbind, lapply(c(0, time.event[1:tauIndEvent]), function(x){
    Hij(Fj = F1.te, surv = sHats.te, iequalj = FALSE, time.jump = c(0, time.event), smin = x, smax = predtimes)
  }))
  
  H22 <- do.call(cbind, lapply(c(0, time.event[1:tauIndEvent]), function(x){
    Hij(Fj = F2.te, surv = sHats.te, iequalj = TRUE, time.jump = c(0, time.event), smin = x, smax = predtimes)
  }))
  
  H12 <- do.call(cbind, lapply(c(0, time.event[1:tauIndEvent]), function(x){
    Hij(Fj = F2.te, surv = sHats.te, iequalj = FALSE, time.jump = c(0, time.event), smin = x, smax = predtimes)
  }))
  
  L1.1_tmp <- rowCumSum((H11[, 1:tauIndEvent]/gHats.te[, 1:tauIndEvent])*dL1.te[, 1:tauIndEvent])
  L1.2_tmp <- rowCumSum((H12[, 1:tauIndEvent]/gHats.te[, 1:tauIndEvent])*dL1.te[, 1:tauIndEvent])
  
  L2.1_tmp <- rowCumSum((H21[, 1:tauIndEvent]/gHats.te[, 1:tauIndEvent])*dL2.te[, 1:tauIndEvent])
  L2.2_tmp <- rowCumSum((H22[, 1:tauIndEvent]/gHats.te[, 1:tauIndEvent])*dL2.te[, 1:tauIndEvent])
  
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
    tmp[mgEvalIndex[, x] > 0] <- (H11[, 1:tauIndEvent]/gHats.te[, 1:tauIndEvent])[cbind(1:nobs, mgEvalIndex[, x])] * dN1[, x][mgEvalIndex[, x] > 0]
    tmp
  }))
  
  N1.2 <- do.call(cbind, lapply(seq_along(predtimes), function(x){
    tmp <- matrix(0, nrow = nobs, ncol = 1)
    tmp[mgEvalIndex[, x] > 0] <- (H12[, 1:tauIndEvent]/gHats.te[, 1:tauIndEvent])[cbind(1:nobs, mgEvalIndex[, x])] * dN1[, x][mgEvalIndex[, x] > 0]
    tmp
  }))
  
  N2.2 <- do.call(cbind, lapply(seq_along(predtimes), function(x){
    tmp <- matrix(0, nrow = nobs, ncol = 1)
    tmp[mgEvalIndex[, x] > 0] <- (H22[, 1:tauIndEvent]/gHats.te[, 1:tauIndEvent])[cbind(1:nobs, mgEvalIndex[, x])] * dN2[, x][mgEvalIndex[, x] > 0]
    tmp
  }))
  
  N2.1 <- do.call(cbind, lapply(seq_along(predtimes), function(x){
    tmp <- matrix(0, nrow = nobs, ncol = 1)
    tmp[mgEvalIndex[, x] > 0] <- (H21[, 1:tauIndEvent]/gHats.te[, 1:tauIndEvent])[cbind(1:nobs, mgEvalIndex[, x])] * dN2[, x][mgEvalIndex[, x] > 0]
    tmp
  }))
  
  M1.1 <- N1.1 - L1.1
  M2.1 <- N2.1 - L2.1
  
  M1.2 <- N1.2 - L1.2
  M2.2 <- N2.2 - L2.2
  
  # calculate CF L_j's
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
  
  #### tau = F1CF ####
  
  L1.1_tmp_AR <- rowCumSum(((F1.te[, 1:tauIndEvent] - F1.te[, tauIndEvent]) /(gHats.te[, 1:tauIndEvent]*sHats.te[, 1:tauIndEvent]) + 1/gHats.te[, 1:tauIndEvent])*dL1.te[, 1:tauIndEvent])
  L1.2_tmp_AR <- rowCumSum(((F2.te[, 1:tauIndEvent] - F2.te[, tauIndEvent]) /(gHats.te[, 1:tauIndEvent]*sHats.te[, 1:tauIndEvent]))*dL1.te[, 1:tauIndEvent])
  
  L2.2_tmp_AR <- rowCumSum(((F2.te[, 1:tauIndEvent] - F2.te[, tauIndEvent]) /(gHats.te[, 1:tauIndEvent]*sHats.te[, 1:tauIndEvent]) + 1/gHats.te[, 1:tauIndEvent])*dL2.te[, 1:tauIndEvent])
  L2.1_tmp_AR <- rowCumSum(((F1.te[, 1:tauIndEvent] - F1.te[, tauIndEvent]) /(gHats.te[, 1:tauIndEvent]*sHats.te[, 1:tauIndEvent]))*dL2.te[, 1:tauIndEvent])
  
  L1.1_AR <- do.call(cbind, lapply(seq_along(predtimes), function(x){
    tmp <- matrix(0, nrow = nobs, ncol = 1)
    tmp[mgEvalIndex[, x] > 0] <- L1.1_tmp_AR[cbind(1:nobs, mgEvalIndex[, x])]
    tmp
  }))
  
  L1.2_AR <- do.call(cbind, lapply(seq_along(predtimes), function(x){
    tmp <- matrix(0, nrow = nobs, ncol = 1)
    tmp[mgEvalIndex[, x] > 0] <- L1.2_tmp_AR[cbind(1:nobs, mgEvalIndex[, x])]
    tmp
  }))
  
  L2.1_AR <- do.call(cbind, lapply(seq_along(predtimes), function(x){
    tmp <- matrix(0, nrow = nobs, ncol = 1)
    tmp[mgEvalIndex[, x] > 0] <- L2.1_tmp_AR[cbind(1:nobs, mgEvalIndex[, x])]
    tmp
  }))
  
  L2.2_AR <- do.call(cbind, lapply(seq_along(predtimes), function(x){
    tmp <- matrix(0, nrow = nobs, ncol = 1)
    tmp[mgEvalIndex[, x] > 0] <- L2.2_tmp_AR[cbind(1:nobs, mgEvalIndex[, x])]
    tmp
  }))
  
  N1.1_AR <- do.call(cbind, lapply(seq_along(predtimes), function(x){
    tmp <- matrix(0, nrow = nobs, ncol = 1)
    tmp[mgEvalIndex[, x] > 0] <- ((F1.te[, 1:tauIndEvent] - F1.te[, tauIndEvent]) /(gHats.te[, 1:tauIndEvent]*sHats.te[, 1:tauIndEvent]) + 1/gHats.te[, 1:tauIndEvent])[cbind(1:nobs, mgEvalIndex[, x])] * dN1[, x][mgEvalIndex[, x] > 0]
    tmp
  }))
  
  N1.2_AR <- do.call(cbind, lapply(seq_along(predtimes), function(x){
    tmp <- matrix(0, nrow = nobs, ncol = 1)
    tmp[mgEvalIndex[, x] > 0] <- ((F2.te[, 1:tauIndEvent] - F2.te[, tauIndEvent]) /(gHats.te[, 1:tauIndEvent]*sHats.te[, 1:tauIndEvent]))[cbind(1:nobs, mgEvalIndex[, x])] * dN1[, x][mgEvalIndex[, x] > 0]
    tmp
  }))
  
  N2.2_AR <- do.call(cbind, lapply(seq_along(predtimes), function(x){
    tmp <- matrix(0, nrow = nobs, ncol = 1)
    tmp[mgEvalIndex[, x] > 0] <- ((F2.te[, 1:tauIndEvent] - F2.te[, tauIndEvent]) /(gHats.te[, 1:tauIndEvent]*sHats.te[, 1:tauIndEvent]) + 1/gHats.te[, 1:tauIndEvent])[cbind(1:nobs, mgEvalIndex[, x])] * dN2[, x][mgEvalIndex[, x] > 0]
    tmp
  }))
  
  N2.1_AR <- do.call(cbind, lapply(seq_along(predtimes), function(x){
    tmp <- matrix(0, nrow = nobs, ncol = 1)
    tmp[mgEvalIndex[, x] > 0] <- ((F1.te[, 1:tauIndEvent] - F1.te[, tauIndEvent]) /(gHats.te[, 1:tauIndEvent]*sHats.te[, 1:tauIndEvent]))[cbind(1:nobs, mgEvalIndex[, x])] * dN2[, x][mgEvalIndex[, x] > 0]
    tmp
  }))
  
  M1.1_AR <- N1.1_AR - L1.1_AR
  M2.1_AR <- N2.1_AR - L2.1_AR
  
  M1.2_AR <- N1.2_AR - L1.2_AR
  M2.2_AR <- N2.2_AR - L2.2_AR
  
  phi1.1_AR <- F1.te.a1[, tauIndEvent] + (1*(newA == 1)/pmhat)*(M1.1_AR + M2.1_AR)
  phi1.0_AR <- F1.te.a0[, tauIndEvent]  + (1*(newA == 0)/(1 - pmhat))*(M1.1_AR + M2.1_AR)
  
  phi2.1_AR <- F2.te.a1[, tauIndEvent] + (1*(newA == 1)/pmhat)*(M1.2_AR + M2.2_AR)
  phi2.0_AR <- F2.te.a0[, tauIndEvent] + (1*(newA == 0)/(1 - pmhat))*(M1.2_AR + M2.2_AR)
  
  phi1_AR <- phi1.1_AR - phi1.0_AR
  phi2_AR <- phi2.1_AR - phi2.0_AR
  
  ###################
  
  # list(phi1 = phi1, phi2 = phi2, phi1.1 = phi1.1, phi1.0 = phi1.0, phi2.1 = phi2.1, phi2.0 = phi2.0,
  #      phi1_AR = phi1_AR, phi2_AR = phi2_AR, phi1.1_AR = phi1.1_AR, phi1.0_AR = phi1.0_AR, phi2.1_AR = phi2.1_AR, phi2.0_AR = phi2.0_AR,
  #      L1 = L_1.1 - L_1.0, L2 = L_2.1 - L_2.0, L_1.1 = L_1.1, L_1.0 = L_1.0, L_2.1 = L_2.1, L_2.0 = L_2.0,
  #      F1 = F1.te[,tauIndEvent], F1CFdiff = F1.te.a1[,tauIndEvent] - F1.te.a0[,tauIndEvent])
  
  list(phi1 = phi1, phi2 = phi2, phi1.1 = phi1.1, phi1.0 = phi1.0, phi2.1 = phi2.1, phi2.0 = phi2.0,
       L1 = L_1.1 - L_1.0, L2 = L_2.1 - L_2.0, L_1.1 = L_1.1, L_1.0 = L_1.0, L_2.1 = L_2.1, L_2.0 = L_2.0)
}































phiLyl_tmp <- function(surv, comp, cens, treat, newtimes, newstatus, newX, newA, predtimes){
  #browser()
  nobs <- length(newtimes)
  time.jump1 <- unique(sort(surv$time[surv$status == 1]))
  time.jump2 <- unique(sort(comp$time[comp$status == 1]))
  time.event <- unique(sort(c(time.jump1, time.jump2)))
  tauInd1 <- max(which(time.jump1 <= predtimes))
  tauInd2 <- max(which(time.jump2 <= predtimes))
  tauIndEvent <- max(which(time.event <= predtimes))
  
  # pmhat
  pmhat <- predict(treat, newX = newX)
  
  # M1 integral
  
  # Lambda1.te <- predict(surv, newtimes = time.event, newX = newX, newA = newA)$cumhaz
  # Lambda2.te <- predict(comp, newtimes = time.event, newX = newX, newA = newA)$cumhaz
  
  Lambda1.obs.t1 <- predict(surv, newtimes = time.jump1, newX = newX, newA = newA)$cumhaz
  Lambda2.obs.t1 <- predict(comp, newtimes = time.jump1, newX = newX, newA = newA)$cumhaz

  Lambda1.obs.t2 <- predict(surv, newtimes = time.jump2, newX = newX, newA = newA)$cumhaz
  Lambda2.obs.t2 <- predict(comp, newtimes = time.jump2, newX = newX, newA = newA)$cumhaz
  
  # sHats.te <- exp(-Lambda1.te - Lambda2.te)
  # gHats.te <- predict(cens, newtimes = time.event, newX = newX, newA = newA)$surv
  
  sHats.t1 <- exp(-Lambda1.obs.t1 - Lambda2.obs.t1)
  gHats.t1 <- predict(cens, newtimes = time.jump1, newX = newX, newA = newA)$surv

  sHats.t2 <- exp(-Lambda1.obs.t2 - Lambda2.obs.t2)
  gHats.t2 <- predict(cens, newtimes = time.jump2, newX = newX, newA = newA)$surv
  
  # dL1.te <- cbind(Lambda1.te[, 1], t(diff(t(Lambda1.te))))
  # dL2.te <- cbind(Lambda2.te[, 1], t(diff(t(Lambda2.te))))
  dL1.t1 <- cbind(Lambda1.obs.t1[, 1], t(diff(t(Lambda1.obs.t1))))
  dL2.t2 <- cbind(Lambda1.obs.t1[, 1], t(diff(t(Lambda1.obs.t1))))
  
  dF1.te <- sHats.te*dL1.te
  dF2.te <- sHats.te*dL2.te
  F1.te <- rowCumSum(dF1.te)
  F2.te <- rowCumSum(dF2.te)
  
  #dL1 <- cbind(Lambda2.obs.t1[, 1], t(diff(t(Lambda1.obs.t1))))
  dN1 <- sapply(predtimes, function(t0){
    1*(newtimes <= t0 & newstatus == 1)
  })
  
  dN2 <- sapply(predtimes, function(t0){
    1*(newtimes <= t0 & newstatus == 2)
  })
  
  mgEvalIndex <- do.call(cbind,lapply(predtimes, function(tau){
    prodlim::sindex(jump.times = time.event, eval.times = pmin(newtimes,tau))
  }))
  
  times.tmp <- time.event
  if(!(predtimes %in% time.event)){
    times.tmp[tauIndEvent + 1] <- predtimes
  }
  times.tmp <- times.tmp[1:(tauIndEvent + 1)]
  
  F1.int.tmp <- t(diff(times.tmp)*t(F1.te[,1:tauIndEvent]))
  F1.int <- rowCumSum(F1.int.tmp[, rev(1:tauIndEvent)])[, rev(1:tauIndEvent)]
  
  F2.int.tmp <- t(diff(times.tmp)*t(F2.te[,1:tauIndEvent]))
  F2.int <- rowCumSum(F2.int.tmp[, rev(1:tauIndEvent)])[, rev(1:tauIndEvent)]
  
  H11 <- (predtimes - time.event[1:tauIndEvent]) +
    t((predtimes - time.event[1:tauIndEvent])*t((F1.te/sHats.te)[, 1:tauIndEvent])) + 
    F1.int/sHats.te[, 1:tauIndEvent]
  
  H21 <- t((predtimes - time.event[1:tauIndEvent])*t((F1.te/sHats.te)[, 1:tauIndEvent])) + 
    F1.int/sHats.te[, 1:tauIndEvent]
  
  H22 <- (predtimes - time.event[1:tauIndEvent]) +
    t((predtimes - time.event[1:tauIndEvent])*t((F2.te/sHats.te)[, 1:tauIndEvent])) + 
    F2.int/sHats.te[, 1:tauIndEvent]
  
  H12 <- t((predtimes - time.event[1:tauIndEvent])*t((F2.te/sHats.te)[, 1:tauIndEvent])) + 
    F2.int/sHats.te[, 1:tauIndEvent]
  
  
  # mgEvalIndex1 <- do.call(cbind,lapply(predtimes, function(tau){
  #   prodlim::sindex(jump.times = time.jump1, eval.times = pmin(newtime,tau))
  # }))
  # 
  # mgEvalIndex2 <- do.call(cbind,lapply(predtimes, function(tau){
  #   prodlim::sindex(jump.times = time.jump2, eval.times = pmin(newtime,tau))
  # }))
  
  L1.1_tmp <- rowCumSum((H11/gHats.te[, 1:tauIndEvent])*dL1.te[, 1:tauIndEvent])
  L1.2_tmp <- rowCumSum((H12/gHats.te[, 1:tauIndEvent])*dL1.te[, 1:tauIndEvent])
  
  L2.1_tmp <- rowCumSum((H21/gHats.te[, 1:tauIndEvent])*dL2.te[, 1:tauIndEvent])
  L2.2_tmp <- rowCumSum((H22/gHats.te[, 1:tauIndEvent])*dL2.te[, 1:tauIndEvent])
  
  # integrand1 <- (exp(-Lambda1.obs.t1 - Lambda2.astar.t1)[,1:tauInd1] - V1obs)/(sHats1[,1:tauInd1]*gHats1[,1:tauInd1])
  # L1_tmp <- rowCumSum(integrand1*dL1[,1:tauInd1])
  
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
    tmp[mgEvalIndex[, x] > 0] <- (H11/gHats.te[, 1:tauIndEvent])[cbind(1:nobs, mgEvalIndex[, x])] * dN1[, x][mgEvalIndex[, x] > 0]
    tmp
  }))
  
  N1.2 <- do.call(cbind, lapply(seq_along(predtimes), function(x){
    tmp <- matrix(0, nrow = nobs, ncol = 1)
    tmp[mgEvalIndex[, x] > 0] <- (H12/gHats.te[, 1:tauIndEvent])[cbind(1:nobs, mgEvalIndex[, x])] * dN1[, x][mgEvalIndex[, x] > 0]
    tmp
  }))
  
  N2.2 <- do.call(cbind, lapply(seq_along(predtimes), function(x){
    tmp <- matrix(0, nrow = nobs, ncol = 1)
    tmp[mgEvalIndex[, x] > 0] <- (H22/gHats.te[, 1:tauIndEvent])[cbind(1:nobs, mgEvalIndex[, x])] * dN2[, x][mgEvalIndex[, x] > 0]
    tmp
  }))
  
  N2.1 <- do.call(cbind, lapply(seq_along(predtimes), function(x){
    tmp <- matrix(0, nrow = nobs, ncol = 1)
    tmp[mgEvalIndex[, x] > 0] <- (H21/gHats.te[, 1:tauIndEvent])[cbind(1:nobs, mgEvalIndex[, x])] * dN2[, x][mgEvalIndex[, x] > 0]
    tmp
  }))
  
  M1.1 <- N1.1 - L1.1
  M2.1 <- N2.1 - L2.1
  
  M1.2 <- N1.2 - L1.2
  M2.2 <- N2.2 - L2.2
  
  
  
  # calculate CF L_j's
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
  phi1.0 <- L_1.0 + (1*(newA == 1)/(1 - pmhat))*(M1.1 + M2.1)
  
  phi2.1 <- L_2.1 + (1*(newA == 1)/pmhat)*(M1.2 + M2.2)
  phi2.0 <- L_2.0 + (1*(newA == 1)/(1 - pmhat))*(M1.2 + M2.2)
  
  phi1 <- phi1.1 - phi1.0
  phi2 <- phi2.1 - phi2.0
  
  #### tau = F1CF ####
  
  L1.1_tmp_AR <- rowCumSum(((F1.te[, 1:tauIndEvent] - F1.te[, tauIndEvent]) /(gHats.te[, 1:tauIndEvent]*sHats.te[, 1:tauIndEvent]) + 1/gHats.te[, 1:tauIndEvent])*dL1.te[, 1:tauIndEvent])
  L1.2_tmp_AR <- rowCumSum(((F1.te[, 1:tauIndEvent] - F1.te[, tauIndEvent]) /(gHats.te[, 1:tauIndEvent]*sHats.te[, 1:tauIndEvent]))*dL1.te[, 1:tauIndEvent])
  
  L2.2_tmp_AR <- rowCumSum(((F2.te[, 1:tauIndEvent] - F2.te[, tauIndEvent]) /(gHats.te[, 1:tauIndEvent]*sHats.te[, 1:tauIndEvent]) + 1/gHats.te[, 1:tauIndEvent])*dL2.te[, 1:tauIndEvent])
  L2.1_tmp_AR <- rowCumSum(((F2.te[, 1:tauIndEvent] - F2.te[, tauIndEvent]) /(gHats.te[, 1:tauIndEvent]*sHats.te[, 1:tauIndEvent]))*dL2.te[, 1:tauIndEvent])
  
  L1.1_AR <- do.call(cbind, lapply(seq_along(predtimes), function(x){
    tmp <- matrix(0, nrow = nobs, ncol = 1)
    tmp[mgEvalIndex[, x] > 0] <- L1.1_tmp_AR[cbind(1:nobs, mgEvalIndex[, x])]
    tmp
  }))
  
  L1.2_AR <- do.call(cbind, lapply(seq_along(predtimes), function(x){
    tmp <- matrix(0, nrow = nobs, ncol = 1)
    tmp[mgEvalIndex[, x] > 0] <- L1.1_tmp_AR[cbind(1:nobs, mgEvalIndex[, x])]
    tmp
  }))
  
  L2.1_AR <- do.call(cbind, lapply(seq_along(predtimes), function(x){
    tmp <- matrix(0, nrow = nobs, ncol = 1)
    tmp[mgEvalIndex[, x] > 0] <- L2.1_tmp_AR[cbind(1:nobs, mgEvalIndex[, x])]
    tmp
  }))
  
  L2.2_AR <- do.call(cbind, lapply(seq_along(predtimes), function(x){
    tmp <- matrix(0, nrow = nobs, ncol = 1)
    tmp[mgEvalIndex[, x] > 0] <- L2.2_tmp_AR[cbind(1:nobs, mgEvalIndex[, x])]
    tmp
  }))
  
  N1.1_AR <- do.call(cbind, lapply(seq_along(predtimes), function(x){
    tmp <- matrix(0, nrow = nobs, ncol = 1)
    tmp[mgEvalIndex[, x] > 0] <- ((F1.te[, 1:tauIndEvent] - F1.te[, tauIndEvent]) /(gHats.te[, 1:tauIndEvent]*sHats.te[, 1:tauIndEvent]) + 1/gHats.te[, 1:tauIndEvent])[cbind(1:nobs, mgEvalIndex[, x])] * dN1[, x][mgEvalIndex[, x] > 0]
    tmp
  }))
  
  N1.2_AR <- do.call(cbind, lapply(seq_along(predtimes), function(x){
    tmp <- matrix(0, nrow = nobs, ncol = 1)
    tmp[mgEvalIndex[, x] > 0] <- ((F1.te[, 1:tauIndEvent] - F1.te[, tauIndEvent]) /(gHats.te[, 1:tauIndEvent]*sHats.te[, 1:tauIndEvent]))[cbind(1:nobs, mgEvalIndex[, x])] * dN1[, x][mgEvalIndex[, x] > 0]
    tmp
  }))
  
  N2.2_AR <- do.call(cbind, lapply(seq_along(predtimes), function(x){
    tmp <- matrix(0, nrow = nobs, ncol = 1)
    tmp[mgEvalIndex[, x] > 0] <- ((F2.te[, 1:tauIndEvent] - F2.te[, tauIndEvent]) /(gHats.te[, 1:tauIndEvent]*sHats.te[, 1:tauIndEvent]) + 1/gHats.te[, 1:tauIndEvent])[cbind(1:nobs, mgEvalIndex[, x])] * dN2[, x][mgEvalIndex[, x] > 0]
    tmp
  }))
  
  N2.1_AR <- do.call(cbind, lapply(seq_along(predtimes), function(x){
    tmp <- matrix(0, nrow = nobs, ncol = 1)
    tmp[mgEvalIndex[, x] > 0] <- ((F2.te[, 1:tauIndEvent] - F2.te[, tauIndEvent]) /(gHats.te[, 1:tauIndEvent]*sHats.te[, 1:tauIndEvent]))[cbind(1:nobs, mgEvalIndex[, x])] * dN2[, x][mgEvalIndex[, x] > 0]
    tmp
  }))
  
  browser()
  
  M1.1_AR <- N1.1_AR - L1.1_AR
  M2.1_AR <- N2.1_AR - L2.1_AR
  
  M1.2_AR <- N1.2_AR - L1.2_AR
  M2.2_AR <- N2.2_AR - L2.2_AR
  
  phi1.1_AR <- F1.te.a1[, tauIndEvent] + (1*(newA == 1)/pmhat)*(M1.1_AR + M2.1_AR)
  phi1.0_AR <- F1.te.a0[, tauIndEvent]  + (1*(newA == 1)/(1 - pmhat))*(M1.1_AR + M2.1_AR)
  
  phi2.1_AR <- F2.te.a1[, tauIndEvent] + (1*(newA == 1)/pmhat)*(M1.2_AR + M2.2_AR)
  phi2.0_AR <- F2.te.a0[, tauIndEvent] + (1*(newA == 1)/(1 - pmhat))*(M1.2_AR + M2.2_AR)
  
  phi1_AR <- phi1.1_AR - phi1.0_AR
  phi2_AR <- phi2.1_AR - phi2.0_AR
  
  ###################
  
  
  list(phi1 = phi1, phi2 = phi2, phi1.1 = phi1.1, phi1.0 = phi1.0, phi2.1 = phi2.1, phi2.0 = phi2.0,
       phi1_AR = phi1_AR, phi2_AR = phi2_AR, phi1.1_AR = phi1.1_AR, phi1.0_AR = phi1.0_AR, phi2.1_AR = phi2.1_AR, phi2.0_AR = phi2.0_AR,
       L1 = L_1.1 - L_1.0, L2 = L_2.1 - L_2.0, L_1.1 = L_1.1, L_1.0 = L_1.0, L_2.1 = L_2.1, L_2.0 = L_2.0,
       F1 = F1.te[,tauIndEvent], F1CFdiff = F1.te.a1[,tauIndEvent] - F1.te.a0[,tauIndEvent])
}






























phiLyl_tmp2 <- function(surv, comp, cens, treat, newtimes, newstatus, newX, newA, predtimes){
  
  #browser()
  
  nobs <- length(newtimes)
  time.jump1 <- unique(sort(surv$time[surv$status == 1]))
  time.jump2 <- unique(sort(comp$time[comp$status == 1]))
  time.event <- unique(sort(c(time.jump1, time.jump2)))
  tauInd1 <- max(which(time.jump1 <= predtimes))
  tauInd2 <- max(which(time.jump2 <= predtimes))
  tauIndEvent <- max(which(time.event <= predtimes))
  time.zero <- c(0, time.event)
  
  # pmhat
  pmhat <- predict(treat, newX = newX)
  pmhat[pmhat == 1] <- 1 - 1e-06
  pmhat[pmhat == 0] <- 1e-06
  
  # M1 integral
  
  Lambda1.te <- predict(surv, newtimes = c(0, time.event), newX = newX, newA = newA)$cumhaz
  Lambda2.te <- predict(comp, newtimes = c(0, time.event), newX = newX, newA = newA)$cumhaz
  
  # Lambda1.obs.t1 <- predict(surv, newtimes = time.jump1, newX = newX, newA = newA)$cumhaz
  # Lambda2.obs.t1 <- predict(comp, newtimes = time.jump1, newX = newX, newA = newA)$cumhaz
  # 
  # Lambda1.obs.t2 <- predict(surv, newtimes = time.jump2, newX = newX, newA = newA)$cumhaz
  # Lambda2.obs.t2 <- predict(comp, newtimes = time.jump2, newX = newX, newA = newA)$cumhaz
  
  sHats.te <- exp(-Lambda1.te - Lambda2.te)
  gHats.te <- predict(cens, newtimes = c(0, time.event), newX = newX, newA = newA)$surv
  
  # sHats.t1 <- exp(-Lambda1.obs.t1 - Lambda2.obs.t1)
  # gHats.t1 <- predict(cens, newtimes = time.jump1, newX = newX, newA = newA)$surv
  # 
  # sHats.t2 <- exp(-Lambda1.obs.t2 - Lambda2.obs.t2)
  # gHats.t2 <- predict(cens, newtimes = time.jump2, newX = newX, newA = newA)$surv
  
  # dL1.te <- cbind(Lambda1.te[, 1], t(diff(t(Lambda1.te))))
  # dL2.te <- cbind(Lambda2.te[, 1], t(diff(t(Lambda2.te))))
  dL1.te <- t(diff(t(Lambda1.te)))
  dL2.te <- t(diff(t(Lambda2.te)))
  
  dF1.te <- sHats.te[,-length(time.event)]*dL1.te
  dF2.te <- sHats.te[,-length(time.event)]*dL2.te
  F1.te <- rowCumSum(dF1.te)
  F2.te <- rowCumSum(dF2.te)
  
  #dL1 <- cbind(Lambda2.obs.t1[, 1], t(diff(t(Lambda1.obs.t1))))
  dN1 <- sapply(predtimes, function(t0){
    1*(newtimes <= t0 & newstatus == 1)
  })
  
  dN2 <- sapply(predtimes, function(t0){
    1*(newtimes <= t0 & newstatus == 2)
  })
  
  mgEvalIndex <- do.call(cbind,lapply(predtimes, function(tau){
    prodlim::sindex(jump.times = time.event, eval.times = pmin(newtimes,tau))
  }))
  
  # if(!(predtimes %in% time.event)){
  #   time.event[tauIndEvent + 1] <- predtimes
  # }
  
  H11 <- t((predtimes - time.event[1:tauIndEvent])*t(1 + F1.te[, 1:tauIndEvent]/sHats.te[, 1:tauIndEvent])) -
    (1/sHats.te[, 1:tauIndEvent]) * rowCumSum(t(diff(time.event[1:(tauIndEvent + 1)]) * t(F1.te[, 1:tauIndEvent]))[, tauIndEvent:1])[, tauIndEvent:1]
  
  H21 <- t((predtimes - time.event[1:tauIndEvent])*t(F1.te[, 1:tauIndEvent]/sHats.te[, 1:tauIndEvent])) -
    (1/sHats.te[, 1:tauIndEvent]) * rowCumSum(t(diff(time.event[1:(tauIndEvent + 1)]) * t(F1.te[, 1:tauIndEvent]))[, tauIndEvent:1])[, tauIndEvent:1]
  
  H22 <- t((predtimes - time.event[1:tauIndEvent])*t(1 + F2.te[, 1:tauIndEvent]/sHats.te[, 1:tauIndEvent])) -
    (1/sHats.te[, 1:tauIndEvent]) * rowCumSum(t(diff(time.event[1:(tauIndEvent + 1)]) * t(F2.te[, 1:tauIndEvent]))[, tauIndEvent:1])[, tauIndEvent:1]
  
  H12 <- t((predtimes - time.event[1:tauIndEvent])*t(F2.te[, 1:tauIndEvent]/sHats.te[, 1:tauIndEvent])) -
    (1/sHats.te[, 1:tauIndEvent]) * rowCumSum(t(diff(time.event[1:(tauIndEvent + 1)]) * t(F2.te[, 1:tauIndEvent]))[, tauIndEvent:1])[, tauIndEvent:1]
  
  
  H11_old <- do.call(cbind, lapply(c(0, time.event[1:tauIndEvent]), function(x){
    Hij(Fj = F1.te, surv = sHats.te, iequalj = TRUE, time.jump = c(0, time.event), smin = x, smax = predtimes)
  }))
  
  # H21 <- do.call(cbind, lapply(c(0, time.event[1:tauIndEvent]), function(x){
  #   Hij(Fj = F1.te, surv = sHats.te, iequalj = FALSE, time.jump = c(0, time.event), smin = x, smax = predtimes)
  # }))
  # 
  # H22 <- do.call(cbind, lapply(c(0, time.event[1:tauIndEvent]), function(x){
  #   Hij(Fj = F2.te, surv = sHats.te, iequalj = TRUE, time.jump = c(0, time.event), smin = x, smax = predtimes)
  # }))
  # 
  # H12 <- do.call(cbind, lapply(c(0, time.event[1:tauIndEvent]), function(x){
  #   Hij(Fj = F2.te, surv = sHats.te, iequalj = FALSE, time.jump = c(0, time.event), smin = x, smax = predtimes)
  # }))
  
  L1.1_tmp <- rowCumSum((H11[, 1:tauIndEvent]/gHats.te[, 1:tauIndEvent])*dL1.te[, 1:tauIndEvent])
  L1.2_tmp <- rowCumSum((H12[, 1:tauIndEvent]/gHats.te[, 1:tauIndEvent])*dL1.te[, 1:tauIndEvent])
  
  L2.1_tmp <- rowCumSum((H21[, 1:tauIndEvent]/gHats.te[, 1:tauIndEvent])*dL2.te[, 1:tauIndEvent])
  L2.2_tmp <- rowCumSum((H22[, 1:tauIndEvent]/gHats.te[, 1:tauIndEvent])*dL2.te[, 1:tauIndEvent])
  
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
    tmp[mgEvalIndex[, x] > 0] <- (H11[, 1:tauIndEvent]/gHats.te[, 1:tauIndEvent])[cbind(1:nobs, mgEvalIndex[, x])] * dN1[, x][mgEvalIndex[, x] > 0]
    tmp
  }))
  
  N1.2 <- do.call(cbind, lapply(seq_along(predtimes), function(x){
    tmp <- matrix(0, nrow = nobs, ncol = 1)
    tmp[mgEvalIndex[, x] > 0] <- (H12[, 1:tauIndEvent]/gHats.te[, 1:tauIndEvent])[cbind(1:nobs, mgEvalIndex[, x])] * dN1[, x][mgEvalIndex[, x] > 0]
    tmp
  }))
  
  N2.2 <- do.call(cbind, lapply(seq_along(predtimes), function(x){
    tmp <- matrix(0, nrow = nobs, ncol = 1)
    tmp[mgEvalIndex[, x] > 0] <- (H22[, 1:tauIndEvent]/gHats.te[, 1:tauIndEvent])[cbind(1:nobs, mgEvalIndex[, x])] * dN2[, x][mgEvalIndex[, x] > 0]
    tmp
  }))
  
  N2.1 <- do.call(cbind, lapply(seq_along(predtimes), function(x){
    tmp <- matrix(0, nrow = nobs, ncol = 1)
    tmp[mgEvalIndex[, x] > 0] <- (H21[, 1:tauIndEvent]/gHats.te[, 1:tauIndEvent])[cbind(1:nobs, mgEvalIndex[, x])] * dN2[, x][mgEvalIndex[, x] > 0]
    tmp
  }))
  
  M1.1 <- N1.1 - L1.1
  M2.1 <- N2.1 - L2.1
  
  M1.2 <- N1.2 - L1.2
  M2.2 <- N2.2 - L2.2
  
  # calculate CF L_j's
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
  
  #### tau = F1CF ####
  
  L1.1_tmp_AR <- rowCumSum(((F1.te[, 1:tauIndEvent] - F1.te[, tauIndEvent]) /(gHats.te[, 1:tauIndEvent]*sHats.te[, 1:tauIndEvent]) + 1/gHats.te[, 1:tauIndEvent])*dL1.te[, 1:tauIndEvent])
  L1.2_tmp_AR <- rowCumSum(((F2.te[, 1:tauIndEvent] - F2.te[, tauIndEvent]) /(gHats.te[, 1:tauIndEvent]*sHats.te[, 1:tauIndEvent]))*dL1.te[, 1:tauIndEvent])
  
  L2.2_tmp_AR <- rowCumSum(((F2.te[, 1:tauIndEvent] - F2.te[, tauIndEvent]) /(gHats.te[, 1:tauIndEvent]*sHats.te[, 1:tauIndEvent]) + 1/gHats.te[, 1:tauIndEvent])*dL2.te[, 1:tauIndEvent])
  L2.1_tmp_AR <- rowCumSum(((F1.te[, 1:tauIndEvent] - F1.te[, tauIndEvent]) /(gHats.te[, 1:tauIndEvent]*sHats.te[, 1:tauIndEvent]))*dL2.te[, 1:tauIndEvent])
  
  L1.1_AR <- do.call(cbind, lapply(seq_along(predtimes), function(x){
    tmp <- matrix(0, nrow = nobs, ncol = 1)
    tmp[mgEvalIndex[, x] > 0] <- L1.1_tmp_AR[cbind(1:nobs, mgEvalIndex[, x])]
    tmp
  }))
  
  L1.2_AR <- do.call(cbind, lapply(seq_along(predtimes), function(x){
    tmp <- matrix(0, nrow = nobs, ncol = 1)
    tmp[mgEvalIndex[, x] > 0] <- L1.2_tmp_AR[cbind(1:nobs, mgEvalIndex[, x])]
    tmp
  }))
  
  L2.1_AR <- do.call(cbind, lapply(seq_along(predtimes), function(x){
    tmp <- matrix(0, nrow = nobs, ncol = 1)
    tmp[mgEvalIndex[, x] > 0] <- L2.1_tmp_AR[cbind(1:nobs, mgEvalIndex[, x])]
    tmp
  }))
  
  L2.2_AR <- do.call(cbind, lapply(seq_along(predtimes), function(x){
    tmp <- matrix(0, nrow = nobs, ncol = 1)
    tmp[mgEvalIndex[, x] > 0] <- L2.2_tmp_AR[cbind(1:nobs, mgEvalIndex[, x])]
    tmp
  }))
  
  N1.1_AR <- do.call(cbind, lapply(seq_along(predtimes), function(x){
    tmp <- matrix(0, nrow = nobs, ncol = 1)
    tmp[mgEvalIndex[, x] > 0] <- ((F1.te[, 1:tauIndEvent] - F1.te[, tauIndEvent]) /(gHats.te[, 1:tauIndEvent]*sHats.te[, 1:tauIndEvent]) + 1/gHats.te[, 1:tauIndEvent])[cbind(1:nobs, mgEvalIndex[, x])] * dN1[, x][mgEvalIndex[, x] > 0]
    tmp
  }))
  
  N1.2_AR <- do.call(cbind, lapply(seq_along(predtimes), function(x){
    tmp <- matrix(0, nrow = nobs, ncol = 1)
    tmp[mgEvalIndex[, x] > 0] <- ((F2.te[, 1:tauIndEvent] - F2.te[, tauIndEvent]) /(gHats.te[, 1:tauIndEvent]*sHats.te[, 1:tauIndEvent]))[cbind(1:nobs, mgEvalIndex[, x])] * dN1[, x][mgEvalIndex[, x] > 0]
    tmp
  }))
  
  N2.2_AR <- do.call(cbind, lapply(seq_along(predtimes), function(x){
    tmp <- matrix(0, nrow = nobs, ncol = 1)
    tmp[mgEvalIndex[, x] > 0] <- ((F2.te[, 1:tauIndEvent] - F2.te[, tauIndEvent]) /(gHats.te[, 1:tauIndEvent]*sHats.te[, 1:tauIndEvent]) + 1/gHats.te[, 1:tauIndEvent])[cbind(1:nobs, mgEvalIndex[, x])] * dN2[, x][mgEvalIndex[, x] > 0]
    tmp
  }))
  
  N2.1_AR <- do.call(cbind, lapply(seq_along(predtimes), function(x){
    tmp <- matrix(0, nrow = nobs, ncol = 1)
    tmp[mgEvalIndex[, x] > 0] <- ((F1.te[, 1:tauIndEvent] - F1.te[, tauIndEvent]) /(gHats.te[, 1:tauIndEvent]*sHats.te[, 1:tauIndEvent]))[cbind(1:nobs, mgEvalIndex[, x])] * dN2[, x][mgEvalIndex[, x] > 0]
    tmp
  }))
  
  M1.1_AR <- N1.1_AR - L1.1_AR
  M2.1_AR <- N2.1_AR - L2.1_AR
  
  M1.2_AR <- N1.2_AR - L1.2_AR
  M2.2_AR <- N2.2_AR - L2.2_AR
  
  phi1.1_AR <- F1.te.a1[, tauIndEvent] + (1*(newA == 1)/pmhat)*(M1.1_AR + M2.1_AR)
  phi1.0_AR <- F1.te.a0[, tauIndEvent]  + (1*(newA == 0)/(1 - pmhat))*(M1.1_AR + M2.1_AR)
  
  phi2.1_AR <- F2.te.a1[, tauIndEvent] + (1*(newA == 1)/pmhat)*(M1.2_AR + M2.2_AR)
  phi2.0_AR <- F2.te.a0[, tauIndEvent] + (1*(newA == 0)/(1 - pmhat))*(M1.2_AR + M2.2_AR)
  
  phi1_AR <- phi1.1_AR - phi1.0_AR
  phi2_AR <- phi2.1_AR - phi2.0_AR
  
  ###################
  
  # list(phi1 = phi1, phi2 = phi2, phi1.1 = phi1.1, phi1.0 = phi1.0, phi2.1 = phi2.1, phi2.0 = phi2.0,
  #      phi1_AR = phi1_AR, phi2_AR = phi2_AR, phi1.1_AR = phi1.1_AR, phi1.0_AR = phi1.0_AR, phi2.1_AR = phi2.1_AR, phi2.0_AR = phi2.0_AR,
  #      L1 = L_1.1 - L_1.0, L2 = L_2.1 - L_2.0, L_1.1 = L_1.1, L_1.0 = L_1.0, L_2.1 = L_2.1, L_2.0 = L_2.0,
  #      F1 = F1.te[,tauIndEvent], F1CFdiff = F1.te.a1[,tauIndEvent] - F1.te.a0[,tauIndEvent])
  
  list(phi1 = phi1, phi2 = phi2, phi1.1 = phi1.1, phi1.0 = phi1.0, phi2.1 = phi2.1, phi2.0 = phi2.0,
       L1 = L_1.1 - L_1.0, L2 = L_2.1 - L_2.0, L_1.1 = L_1.1, L_1.0 = L_1.0, L_2.1 = L_2.1, L_2.0 = L_2.0)
}




























phiLyl <- function(surv, comp, cens, treat, newtimes, newstatus, newX, newA, predtimes){
  
  #browser()
  
  nobs <- length(newtimes)
  time.jump1 <- unique(sort(surv$time[surv$status == 1]))
  time.jump2 <- unique(sort(comp$time[comp$status == 1]))
  time.event <- unique(sort(c(time.jump1, time.jump2)))
  tauInd1 <- max(which(time.jump1 <= predtimes))
  tauInd2 <- max(which(time.jump2 <= predtimes))
  tauIndEvent <- max(which(time.event <= predtimes))
  time.zero <- c(0, time.event)
  tauIndZero <- max(which(time.zero <= predtimes))
  # pmhat
  pmhat <- predict(treat, newX = newX)
  pmhat[pmhat == 1] <- 1 - 1e-06
  pmhat[pmhat == 0] <- 1e-06
  
  # M1 integral
  
  Lambda1.te <- predict(surv, newtimes = time.zero, newX = newX, newA = newA)$cumhaz
  Lambda2.te <- predict(comp, newtimes = time.zero, newX = newX, newA = newA)$cumhaz
  
  # Lambda1.obs.t1 <- predict(surv, newtimes = time.jump1, newX = newX, newA = newA)$cumhaz
  # Lambda2.obs.t1 <- predict(comp, newtimes = time.jump1, newX = newX, newA = newA)$cumhaz
  # 
  # Lambda1.obs.t2 <- predict(surv, newtimes = time.jump2, newX = newX, newA = newA)$cumhaz
  # Lambda2.obs.t2 <- predict(comp, newtimes = time.jump2, newX = newX, newA = newA)$cumhaz
  
  sHats.te <- exp(-Lambda1.te - Lambda2.te)
  gHats.te <- predict(cens, newtimes = time.zero, newX = newX, newA = newA)$surv
  
  # sHats.t1 <- exp(-Lambda1.obs.t1 - Lambda2.obs.t1)
  # gHats.t1 <- predict(cens, newtimes = time.jump1, newX = newX, newA = newA)$surv
  # 
  # sHats.t2 <- exp(-Lambda1.obs.t2 - Lambda2.obs.t2)
  # gHats.t2 <- predict(cens, newtimes = time.jump2, newX = newX, newA = newA)$surv
  
  # dL1.te <- cbind(Lambda1.te[, 1], t(diff(t(Lambda1.te))))
  # dL2.te <- cbind(Lambda2.te[, 1], t(diff(t(Lambda2.te))))
  dL1.te <- t(diff(t(Lambda1.te)))
  dL2.te <- t(diff(t(Lambda2.te)))
  
  dF1.te <- sHats.te[,-length(time.zero)]*dL1.te
  dF2.te <- sHats.te[,-length(time.zero)]*dL2.te
  F1.te <- cbind(0, rowCumSum(dF1.te))
  F2.te <- cbind(0, rowCumSum(dF2.te))
  
  #dL1 <- cbind(Lambda2.obs.t1[, 1], t(diff(t(Lambda1.obs.t1))))
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
  
  
  H11_old <- do.call(cbind, lapply(c(0, time.event[1:tauIndEvent]), function(x){
    Hij(Fj = F1.te, surv = sHats.te, iequalj = TRUE, time.jump = c(0, time.event), smin = x, smax = predtimes)
  }))
  
  # H21 <- do.call(cbind, lapply(c(0, time.event[1:tauIndEvent]), function(x){
  #   Hij(Fj = F1.te, surv = sHats.te, iequalj = FALSE, time.jump = c(0, time.event), smin = x, smax = predtimes)
  # }))
  # 
  # H22 <- do.call(cbind, lapply(c(0, time.event[1:tauIndEvent]), function(x){
  #   Hij(Fj = F2.te, surv = sHats.te, iequalj = TRUE, time.jump = c(0, time.event), smin = x, smax = predtimes)
  # }))
  # 
  # H12 <- do.call(cbind, lapply(c(0, time.event[1:tauIndEvent]), function(x){
  #   Hij(Fj = F2.te, surv = sHats.te, iequalj = FALSE, time.jump = c(0, time.event), smin = x, smax = predtimes)
  # }))
  
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
  
  # calculate CF L_j's
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
  
  #### tau = F1CF ####
  
  # L1.1_tmp_AR <- rowCumSum(((F1.te[, 1:tauIndEvent] - F1.te[, tauIndEvent]) /(gHats.te[, 1:tauIndEvent]*sHats.te[, 1:tauIndEvent]) + 1/gHats.te[, 1:tauIndEvent])*dL1.te[, 1:tauIndEvent])
  # L1.2_tmp_AR <- rowCumSum(((F2.te[, 1:tauIndEvent] - F2.te[, tauIndEvent]) /(gHats.te[, 1:tauIndEvent]*sHats.te[, 1:tauIndEvent]))*dL1.te[, 1:tauIndEvent])
  # 
  # L2.2_tmp_AR <- rowCumSum(((F2.te[, 1:tauIndEvent] - F2.te[, tauIndEvent]) /(gHats.te[, 1:tauIndEvent]*sHats.te[, 1:tauIndEvent]) + 1/gHats.te[, 1:tauIndEvent])*dL2.te[, 1:tauIndEvent])
  # L2.1_tmp_AR <- rowCumSum(((F1.te[, 1:tauIndEvent] - F1.te[, tauIndEvent]) /(gHats.te[, 1:tauIndEvent]*sHats.te[, 1:tauIndEvent]))*dL2.te[, 1:tauIndEvent])
  # 
  # L1.1_AR <- do.call(cbind, lapply(seq_along(predtimes), function(x){
  #   tmp <- matrix(0, nrow = nobs, ncol = 1)
  #   tmp[mgEvalIndex[, x] > 0] <- L1.1_tmp_AR[cbind(1:nobs, mgEvalIndex[, x])]
  #   tmp
  # }))
  # 
  # L1.2_AR <- do.call(cbind, lapply(seq_along(predtimes), function(x){
  #   tmp <- matrix(0, nrow = nobs, ncol = 1)
  #   tmp[mgEvalIndex[, x] > 0] <- L1.2_tmp_AR[cbind(1:nobs, mgEvalIndex[, x])]
  #   tmp
  # }))
  # 
  # L2.1_AR <- do.call(cbind, lapply(seq_along(predtimes), function(x){
  #   tmp <- matrix(0, nrow = nobs, ncol = 1)
  #   tmp[mgEvalIndex[, x] > 0] <- L2.1_tmp_AR[cbind(1:nobs, mgEvalIndex[, x])]
  #   tmp
  # }))
  # 
  # L2.2_AR <- do.call(cbind, lapply(seq_along(predtimes), function(x){
  #   tmp <- matrix(0, nrow = nobs, ncol = 1)
  #   tmp[mgEvalIndex[, x] > 0] <- L2.2_tmp_AR[cbind(1:nobs, mgEvalIndex[, x])]
  #   tmp
  # }))
  # 
  # N1.1_AR <- do.call(cbind, lapply(seq_along(predtimes), function(x){
  #   tmp <- matrix(0, nrow = nobs, ncol = 1)
  #   tmp[mgEvalIndex[, x] > 0] <- ((F1.te[, 1:tauIndEvent] - F1.te[, tauIndEvent]) /(gHats.te[, 1:tauIndEvent]*sHats.te[, 1:tauIndEvent]) + 1/gHats.te[, 1:tauIndEvent])[cbind(1:nobs, mgEvalIndex[, x])] * dN1[, x][mgEvalIndex[, x] > 0]
  #   tmp
  # }))
  # 
  # N1.2_AR <- do.call(cbind, lapply(seq_along(predtimes), function(x){
  #   tmp <- matrix(0, nrow = nobs, ncol = 1)
  #   tmp[mgEvalIndex[, x] > 0] <- ((F2.te[, 1:tauIndEvent] - F2.te[, tauIndEvent]) /(gHats.te[, 1:tauIndEvent]*sHats.te[, 1:tauIndEvent]))[cbind(1:nobs, mgEvalIndex[, x])] * dN1[, x][mgEvalIndex[, x] > 0]
  #   tmp
  # }))
  # 
  # N2.2_AR <- do.call(cbind, lapply(seq_along(predtimes), function(x){
  #   tmp <- matrix(0, nrow = nobs, ncol = 1)
  #   tmp[mgEvalIndex[, x] > 0] <- ((F2.te[, 1:tauIndEvent] - F2.te[, tauIndEvent]) /(gHats.te[, 1:tauIndEvent]*sHats.te[, 1:tauIndEvent]) + 1/gHats.te[, 1:tauIndEvent])[cbind(1:nobs, mgEvalIndex[, x])] * dN2[, x][mgEvalIndex[, x] > 0]
  #   tmp
  # }))
  # 
  # N2.1_AR <- do.call(cbind, lapply(seq_along(predtimes), function(x){
  #   tmp <- matrix(0, nrow = nobs, ncol = 1)
  #   tmp[mgEvalIndex[, x] > 0] <- ((F1.te[, 1:tauIndEvent] - F1.te[, tauIndEvent]) /(gHats.te[, 1:tauIndEvent]*sHats.te[, 1:tauIndEvent]))[cbind(1:nobs, mgEvalIndex[, x])] * dN2[, x][mgEvalIndex[, x] > 0]
  #   tmp
  # }))
  # 
  # M1.1_AR <- N1.1_AR - L1.1_AR
  # M2.1_AR <- N2.1_AR - L2.1_AR
  # 
  # M1.2_AR <- N1.2_AR - L1.2_AR
  # M2.2_AR <- N2.2_AR - L2.2_AR
  # 
  # phi1.1_AR <- F1.te.a1[, tauIndEvent] + (1*(newA == 1)/pmhat)*(M1.1_AR + M2.1_AR)
  # phi1.0_AR <- F1.te.a0[, tauIndEvent]  + (1*(newA == 0)/(1 - pmhat))*(M1.1_AR + M2.1_AR)
  # 
  # phi2.1_AR <- F2.te.a1[, tauIndEvent] + (1*(newA == 1)/pmhat)*(M1.2_AR + M2.2_AR)
  # phi2.0_AR <- F2.te.a0[, tauIndEvent] + (1*(newA == 0)/(1 - pmhat))*(M1.2_AR + M2.2_AR)
  # 
  # phi1_AR <- phi1.1_AR - phi1.0_AR
  # phi2_AR <- phi2.1_AR - phi2.0_AR
  
  ###################
  
  # list(phi1 = phi1, phi2 = phi2, phi1.1 = phi1.1, phi1.0 = phi1.0, phi2.1 = phi2.1, phi2.0 = phi2.0,
  #      phi1_AR = phi1_AR, phi2_AR = phi2_AR, phi1.1_AR = phi1.1_AR, phi1.0_AR = phi1.0_AR, phi2.1_AR = phi2.1_AR, phi2.0_AR = phi2.0_AR,
  #      L1 = L_1.1 - L_1.0, L2 = L_2.1 - L_2.0, L_1.1 = L_1.1, L_1.0 = L_1.0, L_2.1 = L_2.1, L_2.0 = L_2.0,
  #      F1 = F1.te[,tauIndEvent], F1CFdiff = F1.te.a1[,tauIndEvent] - F1.te.a0[,tauIndEvent])
  
  list(phi1 = phi1, phi2 = phi2, phi1.1 = phi1.1, phi1.0 = phi1.0, phi2.1 = phi2.1, phi2.0 = phi2.0,
       L1 = L_1.1 - L_1.0, L2 = L_2.1 - L_2.0, L_1.1 = L_1.1, L_1.0 = L_1.0, L_2.1 = L_2.1, L_2.0 = L_2.0)
}



































phiLyl_backup <- function(surv, comp, cens, treat, newtimes, newstatus, newX, newA, predtimes){
  
  browser()
  
  nobs <- length(newtimes)
  time.jump1 <- unique(sort(surv$time[surv$status == 1]))
  time.jump2 <- unique(sort(comp$time[comp$status == 1]))
  time.event <- unique(sort(c(time.jump1, time.jump2)))
  tauInd1 <- max(which(time.jump1 <= predtimes))
  tauInd2 <- max(which(time.jump2 <= predtimes))
  tauIndEvent <- max(which(time.event <= predtimes))
  time.zero <- c(0, time.event)
  
  # pmhat
  pmhat <- predict(treat, newX = newX)
  pmhat[pmhat == 1] <- 1 - 1e-06
  pmhat[pmhat == 0] <- 1e-06
  
  # M1 integral
  
  Lambda1.te <- predict(surv, newtimes = c(0, time.event), newX = newX, newA = newA)$cumhaz
  Lambda2.te <- predict(comp, newtimes = c(0, time.event), newX = newX, newA = newA)$cumhaz
  
  # Lambda1.obs.t1 <- predict(surv, newtimes = time.jump1, newX = newX, newA = newA)$cumhaz
  # Lambda2.obs.t1 <- predict(comp, newtimes = time.jump1, newX = newX, newA = newA)$cumhaz
  # 
  # Lambda1.obs.t2 <- predict(surv, newtimes = time.jump2, newX = newX, newA = newA)$cumhaz
  # Lambda2.obs.t2 <- predict(comp, newtimes = time.jump2, newX = newX, newA = newA)$cumhaz
  
  sHats.te <- exp(-Lambda1.te - Lambda2.te)
  gHats.te <- predict(cens, newtimes = c(0, time.event), newX = newX, newA = newA)$surv
  
  # sHats.t1 <- exp(-Lambda1.obs.t1 - Lambda2.obs.t1)
  # gHats.t1 <- predict(cens, newtimes = time.jump1, newX = newX, newA = newA)$surv
  # 
  # sHats.t2 <- exp(-Lambda1.obs.t2 - Lambda2.obs.t2)
  # gHats.t2 <- predict(cens, newtimes = time.jump2, newX = newX, newA = newA)$surv
  
  # dL1.te <- cbind(Lambda1.te[, 1], t(diff(t(Lambda1.te))))
  # dL2.te <- cbind(Lambda2.te[, 1], t(diff(t(Lambda2.te))))
  dL1.te <- t(diff(t(Lambda1.te)))
  dL2.te <- t(diff(t(Lambda2.te)))
  
  dF1.te <- sHats.te[,-length(time.event)]*dL1.te
  dF2.te <- sHats.te[,-length(time.event)]*dL2.te
  F1.te <- rowCumSum(dF1.te)
  F2.te <- rowCumSum(dF2.te)
  
  #dL1 <- cbind(Lambda2.obs.t1[, 1], t(diff(t(Lambda1.obs.t1))))
  dN1 <- sapply(predtimes, function(t0){
    1*(newtimes <= t0 & newstatus == 1)
  })
  
  dN2 <- sapply(predtimes, function(t0){
    1*(newtimes <= t0 & newstatus == 2)
  })
  
  mgEvalIndex <- do.call(cbind,lapply(predtimes, function(tau){
    prodlim::sindex(jump.times = time.event, eval.times = pmin(newtimes,tau))
  }))
  
  # if(!(predtimes %in% time.event)){
  #   time.event[tauIndEvent + 1] <- predtimes
  # }
  
  H11 <- t((predtimes - time.event[1:tauIndEvent])*t(1 + F1.te[, 1:tauIndEvent]/sHats.te[, 1:tauIndEvent])) -
    (1/sHats.te[, 1:tauIndEvent]) * rowCumSum(t(diff(time.event[1:(tauIndEvent + 1)]) * t(F1.te[, 1:tauIndEvent]))[, tauIndEvent:1])[, tauIndEvent:1]
  
  H21 <- t((predtimes - time.event[1:tauIndEvent])*t(F1.te[, 1:tauIndEvent]/sHats.te[, 1:tauIndEvent])) -
    (1/sHats.te[, 1:tauIndEvent]) * rowCumSum(t(diff(time.event[1:(tauIndEvent + 1)]) * t(F1.te[, 1:tauIndEvent]))[, tauIndEvent:1])[, tauIndEvent:1]
  
  H22 <- t((predtimes - time.event[1:tauIndEvent])*t(1 + F2.te[, 1:tauIndEvent]/sHats.te[, 1:tauIndEvent])) -
    (1/sHats.te[, 1:tauIndEvent]) * rowCumSum(t(diff(time.event[1:(tauIndEvent + 1)]) * t(F2.te[, 1:tauIndEvent]))[, tauIndEvent:1])[, tauIndEvent:1]
  
  H12 <- t((predtimes - time.event[1:tauIndEvent])*t(F2.te[, 1:tauIndEvent]/sHats.te[, 1:tauIndEvent])) -
    (1/sHats.te[, 1:tauIndEvent]) * rowCumSum(t(diff(time.event[1:(tauIndEvent + 1)]) * t(F2.te[, 1:tauIndEvent]))[, tauIndEvent:1])[, tauIndEvent:1]
  
  
  H11_old <- do.call(cbind, lapply(c(0, time.event[1:tauIndEvent]), function(x){
    Hij(Fj = F1.te, surv = sHats.te, iequalj = TRUE, time.jump = c(0, time.event), smin = x, smax = predtimes)
  }))
  
  # H21 <- do.call(cbind, lapply(c(0, time.event[1:tauIndEvent]), function(x){
  #   Hij(Fj = F1.te, surv = sHats.te, iequalj = FALSE, time.jump = c(0, time.event), smin = x, smax = predtimes)
  # }))
  # 
  # H22 <- do.call(cbind, lapply(c(0, time.event[1:tauIndEvent]), function(x){
  #   Hij(Fj = F2.te, surv = sHats.te, iequalj = TRUE, time.jump = c(0, time.event), smin = x, smax = predtimes)
  # }))
  # 
  # H12 <- do.call(cbind, lapply(c(0, time.event[1:tauIndEvent]), function(x){
  #   Hij(Fj = F2.te, surv = sHats.te, iequalj = FALSE, time.jump = c(0, time.event), smin = x, smax = predtimes)
  # }))
  
  L1.1_tmp <- rowCumSum((H11[, 1:tauIndEvent]/gHats.te[, 1:tauIndEvent])*dL1.te[, 1:tauIndEvent])
  L1.2_tmp <- rowCumSum((H12[, 1:tauIndEvent]/gHats.te[, 1:tauIndEvent])*dL1.te[, 1:tauIndEvent])
  
  L2.1_tmp <- rowCumSum((H21[, 1:tauIndEvent]/gHats.te[, 1:tauIndEvent])*dL2.te[, 1:tauIndEvent])
  L2.2_tmp <- rowCumSum((H22[, 1:tauIndEvent]/gHats.te[, 1:tauIndEvent])*dL2.te[, 1:tauIndEvent])
  
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
    tmp[mgEvalIndex[, x] > 0] <- (H11[, 1:tauIndEvent]/gHats.te[, 1:tauIndEvent])[cbind(1:nobs, mgEvalIndex[, x])] * dN1[, x][mgEvalIndex[, x] > 0]
    tmp
  }))
  
  N1.2 <- do.call(cbind, lapply(seq_along(predtimes), function(x){
    tmp <- matrix(0, nrow = nobs, ncol = 1)
    tmp[mgEvalIndex[, x] > 0] <- (H12[, 1:tauIndEvent]/gHats.te[, 1:tauIndEvent])[cbind(1:nobs, mgEvalIndex[, x])] * dN1[, x][mgEvalIndex[, x] > 0]
    tmp
  }))
  
  N2.2 <- do.call(cbind, lapply(seq_along(predtimes), function(x){
    tmp <- matrix(0, nrow = nobs, ncol = 1)
    tmp[mgEvalIndex[, x] > 0] <- (H22[, 1:tauIndEvent]/gHats.te[, 1:tauIndEvent])[cbind(1:nobs, mgEvalIndex[, x])] * dN2[, x][mgEvalIndex[, x] > 0]
    tmp
  }))
  
  N2.1 <- do.call(cbind, lapply(seq_along(predtimes), function(x){
    tmp <- matrix(0, nrow = nobs, ncol = 1)
    tmp[mgEvalIndex[, x] > 0] <- (H21[, 1:tauIndEvent]/gHats.te[, 1:tauIndEvent])[cbind(1:nobs, mgEvalIndex[, x])] * dN2[, x][mgEvalIndex[, x] > 0]
    tmp
  }))
  
  M1.1 <- N1.1 - L1.1
  M2.1 <- N2.1 - L2.1
  
  M1.2 <- N1.2 - L1.2
  M2.2 <- N2.2 - L2.2
  
  # calculate CF L_j's
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
  
  #### tau = F1CF ####
  
  L1.1_tmp_AR <- rowCumSum(((F1.te[, 1:tauIndEvent] - F1.te[, tauIndEvent]) /(gHats.te[, 1:tauIndEvent]*sHats.te[, 1:tauIndEvent]) + 1/gHats.te[, 1:tauIndEvent])*dL1.te[, 1:tauIndEvent])
  L1.2_tmp_AR <- rowCumSum(((F2.te[, 1:tauIndEvent] - F2.te[, tauIndEvent]) /(gHats.te[, 1:tauIndEvent]*sHats.te[, 1:tauIndEvent]))*dL1.te[, 1:tauIndEvent])
  
  L2.2_tmp_AR <- rowCumSum(((F2.te[, 1:tauIndEvent] - F2.te[, tauIndEvent]) /(gHats.te[, 1:tauIndEvent]*sHats.te[, 1:tauIndEvent]) + 1/gHats.te[, 1:tauIndEvent])*dL2.te[, 1:tauIndEvent])
  L2.1_tmp_AR <- rowCumSum(((F1.te[, 1:tauIndEvent] - F1.te[, tauIndEvent]) /(gHats.te[, 1:tauIndEvent]*sHats.te[, 1:tauIndEvent]))*dL2.te[, 1:tauIndEvent])
  
  L1.1_AR <- do.call(cbind, lapply(seq_along(predtimes), function(x){
    tmp <- matrix(0, nrow = nobs, ncol = 1)
    tmp[mgEvalIndex[, x] > 0] <- L1.1_tmp_AR[cbind(1:nobs, mgEvalIndex[, x])]
    tmp
  }))
  
  L1.2_AR <- do.call(cbind, lapply(seq_along(predtimes), function(x){
    tmp <- matrix(0, nrow = nobs, ncol = 1)
    tmp[mgEvalIndex[, x] > 0] <- L1.2_tmp_AR[cbind(1:nobs, mgEvalIndex[, x])]
    tmp
  }))
  
  L2.1_AR <- do.call(cbind, lapply(seq_along(predtimes), function(x){
    tmp <- matrix(0, nrow = nobs, ncol = 1)
    tmp[mgEvalIndex[, x] > 0] <- L2.1_tmp_AR[cbind(1:nobs, mgEvalIndex[, x])]
    tmp
  }))
  
  L2.2_AR <- do.call(cbind, lapply(seq_along(predtimes), function(x){
    tmp <- matrix(0, nrow = nobs, ncol = 1)
    tmp[mgEvalIndex[, x] > 0] <- L2.2_tmp_AR[cbind(1:nobs, mgEvalIndex[, x])]
    tmp
  }))
  
  N1.1_AR <- do.call(cbind, lapply(seq_along(predtimes), function(x){
    tmp <- matrix(0, nrow = nobs, ncol = 1)
    tmp[mgEvalIndex[, x] > 0] <- ((F1.te[, 1:tauIndEvent] - F1.te[, tauIndEvent]) /(gHats.te[, 1:tauIndEvent]*sHats.te[, 1:tauIndEvent]) + 1/gHats.te[, 1:tauIndEvent])[cbind(1:nobs, mgEvalIndex[, x])] * dN1[, x][mgEvalIndex[, x] > 0]
    tmp
  }))
  
  N1.2_AR <- do.call(cbind, lapply(seq_along(predtimes), function(x){
    tmp <- matrix(0, nrow = nobs, ncol = 1)
    tmp[mgEvalIndex[, x] > 0] <- ((F2.te[, 1:tauIndEvent] - F2.te[, tauIndEvent]) /(gHats.te[, 1:tauIndEvent]*sHats.te[, 1:tauIndEvent]))[cbind(1:nobs, mgEvalIndex[, x])] * dN1[, x][mgEvalIndex[, x] > 0]
    tmp
  }))
  
  N2.2_AR <- do.call(cbind, lapply(seq_along(predtimes), function(x){
    tmp <- matrix(0, nrow = nobs, ncol = 1)
    tmp[mgEvalIndex[, x] > 0] <- ((F2.te[, 1:tauIndEvent] - F2.te[, tauIndEvent]) /(gHats.te[, 1:tauIndEvent]*sHats.te[, 1:tauIndEvent]) + 1/gHats.te[, 1:tauIndEvent])[cbind(1:nobs, mgEvalIndex[, x])] * dN2[, x][mgEvalIndex[, x] > 0]
    tmp
  }))
  
  N2.1_AR <- do.call(cbind, lapply(seq_along(predtimes), function(x){
    tmp <- matrix(0, nrow = nobs, ncol = 1)
    tmp[mgEvalIndex[, x] > 0] <- ((F1.te[, 1:tauIndEvent] - F1.te[, tauIndEvent]) /(gHats.te[, 1:tauIndEvent]*sHats.te[, 1:tauIndEvent]))[cbind(1:nobs, mgEvalIndex[, x])] * dN2[, x][mgEvalIndex[, x] > 0]
    tmp
  }))
  
  M1.1_AR <- N1.1_AR - L1.1_AR
  M2.1_AR <- N2.1_AR - L2.1_AR
  
  M1.2_AR <- N1.2_AR - L1.2_AR
  M2.2_AR <- N2.2_AR - L2.2_AR
  
  phi1.1_AR <- F1.te.a1[, tauIndEvent] + (1*(newA == 1)/pmhat)*(M1.1_AR + M2.1_AR)
  phi1.0_AR <- F1.te.a0[, tauIndEvent]  + (1*(newA == 0)/(1 - pmhat))*(M1.1_AR + M2.1_AR)
  
  phi2.1_AR <- F2.te.a1[, tauIndEvent] + (1*(newA == 1)/pmhat)*(M1.2_AR + M2.2_AR)
  phi2.0_AR <- F2.te.a0[, tauIndEvent] + (1*(newA == 0)/(1 - pmhat))*(M1.2_AR + M2.2_AR)
  
  phi1_AR <- phi1.1_AR - phi1.0_AR
  phi2_AR <- phi2.1_AR - phi2.0_AR
  
  ###################
  
  # list(phi1 = phi1, phi2 = phi2, phi1.1 = phi1.1, phi1.0 = phi1.0, phi2.1 = phi2.1, phi2.0 = phi2.0,
  #      phi1_AR = phi1_AR, phi2_AR = phi2_AR, phi1.1_AR = phi1.1_AR, phi1.0_AR = phi1.0_AR, phi2.1_AR = phi2.1_AR, phi2.0_AR = phi2.0_AR,
  #      L1 = L_1.1 - L_1.0, L2 = L_2.1 - L_2.0, L_1.1 = L_1.1, L_1.0 = L_1.0, L_2.1 = L_2.1, L_2.0 = L_2.0,
  #      F1 = F1.te[,tauIndEvent], F1CFdiff = F1.te.a1[,tauIndEvent] - F1.te.a0[,tauIndEvent])
  
  list(phi1 = phi1, phi2 = phi2, phi1.1 = phi1.1, phi1.0 = phi1.0, phi2.1 = phi2.1, phi2.0 = phi2.0,
       L1 = L_1.1 - L_1.0, L2 = L_2.1 - L_2.0, L_1.1 = L_1.1, L_1.0 = L_1.0, L_2.1 = L_2.1, L_2.0 = L_2.0)
}
