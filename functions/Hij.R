Hij <- function(Fj, surv, iequalj = TRUE, time.jump, smin, smax){
  maxInd <- max(which(time.jump <= smax))
  minInd <- min(which(time.jump >= smin))
  times.tmp <- time.jump
  if(!(smax %in% time.jump)){
    times.tmp[maxInd + 1] <- smax
  }
  times.tmp <- times.tmp[minInd:(maxInd + 1)]
  evalInd <- which(time.jump >= smin & time.jump <= smax)
  tmp <- 1*iequalj + (Fj[, evalInd[1]] - Fj[, evalInd])/surv[, evalInd[1]]
  res <- rowSums(t(diff(times.tmp)*t(tmp)))
  res
}