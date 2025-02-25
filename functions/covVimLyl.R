covVimLyl <- function(time,
                             status,
                             treatment,
                             confounders,
                             l,
                             evaltimes,
                             nuisance = list(),
                             tauLearner = "S",
                             folds = 1){
  #browser()
  require(riskRegression)
  A <- treatment
  treatLabs <- levels(treatment)
  treatment <- 1*(treatment == treatLabs[2])
  confounders <- as.data.frame(confounders)
  xl_numeric <- list()
  for(i in seq_along(l)){
    if(class(confounders[, i]) == "numeric"){
      xl_numeric[[i]] <- confounders[, i]
    } else if(class(confounders[, i]) == "factor"){
      xl_numeric_labs <- levels(confounders[, i])
      xl_numeric[[i]] <- 1*(confounders[, i] == xl_numeric_labs[2])
    }
  }
  
  N <- length(time)
  lminus <- lapply(l, function(x){setdiff(1:dim(confounders)[2], x)})
  k <- sample(1:folds, length(status), replace = TRUE)
  if(folds == 1){
    pi <- do.call(propML, c(nuisance$prop, list(X = confounders, A = A)))
    if(tauLearner == "S"){
      S1 <- do.call(survML, c(nuisance$surv, list(time = time, status = 1*(status == 1), X = confounders, A = treatment)))
      S2 <- do.call(survML, c(nuisance$comp, list(time = time, status = 1*(status == 2), X = confounders, A = treatment)))
      C <- do.call(survML, c(nuisance$cens, list(time = time, status = 1*(status == 0), X = confounders, A = treatment)))
      phi_tmp <- phiLyl(surv = S1,
                     comp = S2,
                        cens = C,
                        treat = pi,
                        newtimes = time,
                        newstatus = status,
                        newX = confounders,
                        newA = treatment,
                        predtimes = evaltimes)
      phi <- phi_tmp$phi1
      eif_L1.1_uncen <- phi_tmp$phi1.1
      eif_L1.0_uncen <- phi_tmp$phi1.0
      tau <- phi_tmp$L1
      
      taul_tmp <- lapply(seq_along(l), function(x){
        do.call(pseudoReg, c(nuisance$pseudoRegTau, list(pseudo = tau[,1], X = confounders[, lminus[[x]], drop = FALSE])))
      })
      taul <- lapply(seq_along(l), function(x){
        predict(taul_tmp[[x]], newX = confounders[, lminus[[x]], drop = FALSE])
      })
      
      EXl_Xlminus_tmp <- lapply(seq_along(l), function(x){
        do.call(pseudoReg, c(nuisance$pseudoRegX, list(pseudo = confounders[, l[x]], X = confounders[, lminus[[x]], drop = FALSE])))
      })
      EXl_Xlminus <- lapply(seq_along(l), function(x){
        predict(EXl_Xlminus_tmp[[x]], newX = confounders[, lminus[[x]], drop = FALSE])
      })
      
    } else if(tauLearner == "T"){
      S1.1 <- do.call(survML, c(nuisance$surv, list(time = time[treatment == 1], status = 1*(status[treatment == 1] == 1), X = confounders[treatment == 1, ])))
      S1.0 <- do.call(survML, c(nuisance$surv, list(time = time[treatment == 0], status = 1*(status[treatment == 0] == 1), X = confounders[treatment == 0, ])))
      S2.1 <- do.call(survML, c(nuisance$surv, list(time = time[treatment == 1], status = 1*(status[treatment == 1] == 2), X = confounders[treatment == 1, ])))
      S2.0 <- do.call(survML, c(nuisance$surv, list(time = time[treatment == 0], status = 1*(status[treatment == 0] == 2), X = confounders[treatment == 0, ])))
      C.1 <- do.call(survML, c(nuisance$cens, list(time = time[treatment == 1], status = 1*(status[treatment == 1]==0), X = confounders[treatment == 1, ])))
      C.0 <- do.call(survML, c(nuisance$cens, list(time = time[treatment == 0], status = 1*(status[treatment == 0]==0), X = confounders[treatment == 0, ])))
      phi1 <- phiLyl(surv = S1.1,
                      comp = S2.1,
                         cens = C.1,
                         treat = pi,
                         newtimes = time,
                         newstatus = status,
                         newX = confounders,
                         newA = treatment,
                         predtimes = evaltimes)
      phi0 <- phiLyl(surv = S1.0,
                      comp = S2.0,
                         cens = C.0,
                         treat = pi,
                         newtimes = time,
                         newstatus = status,
                         newX = confounders,
                         newA = treatment,
                         predtimes = evaltimes)
      phi <- phi1$phi1.1 - phi0$phi1.0
      eif_L1.1_uncen <- phi1$phi1.1
      eif_L1.0_uncen <- phi0$phi1.0
      tau <- phi1$L_1.1 - phi0$L_1.0
      taul_tmp <- lapply(seq_along(l), function(x){
        do.call(pseudoReg, c(nuisance$pseudoRegTau, list(pseudo = tau[,1], X = confounders[, lminus[[x]], drop = FALSE])))
      })
      taul <- lapply(seq_along(l), function(x){
        predict(taul_tmp[[x]], newX = confounders[, lminus[[x]], drop = FALSE])
      })
      
      EXl_Xlminus_tmp <- lapply(seq_along(l), function(x){
        do.call(pseudoReg, c(nuisance$pseudoRegX, list(pseudo = confounders[, l[x]], X = confounders[, lminus[[x]], drop = FALSE])))
      })
      EXl_Xlminus <- lapply(seq_along(l), function(x){
        predict(EXl_Xlminus_tmp[[x]], newX = confounders[, lminus[[x]], drop = FALSE])
      })
    }
    
    res_tmp <- lapply(seq_along(l), function(x){
      eif_Gamma_j_uncen <- (phi - taul[[x]])*(xl_numeric[[x]] - EXl_Xlminus[[x]])
      Gamma_j <- mean(eif_Gamma_j_uncen)
      eif_Gamma_j <- eif_Gamma_j_uncen - Gamma_j
      var_Gamma_j <- (1/N)*mean(eif_Gamma_j^2)
      se_Gamma_j <- sqrt(var_Gamma_j)
      pval_Gamma_j <- 2*(1 - pnorm(abs(Gamma_j)/sqrt(var_Gamma_j)))
      
      eif_chi_j_uncen <- (xl_numeric[[x]] - EXl_Xlminus[[x]])^2
      chi_j <- mean(eif_chi_j_uncen)
      eif_chi_j <- eif_chi_j_uncen - chi_j
      var_chi_j <- (1/N)*mean(eif_chi_j^2)
      se_chi_j <- sqrt(var_chi_j)
      
      omega_j <- Gamma_j/chi_j
      eif_omega_j <- (eif_Gamma_j - omega_j*eif_chi_j)/chi_j
      var_omega_j <- (1/N)*mean(eif_omega_j^2)
      se_omega_j <- sqrt(var_omega_j)
      pval_omega_j <- 2*(1 - pnorm(abs(omega_j)/sqrt(var_omega_j)))
      
      data.frame(l = names(confounders)[l[x]],
                 Gamma_j = Gamma_j, se_Gamma_j = se_Gamma_j, pval_Gamma_j = pval_Gamma_j,
                 chi_j = chi_j, se_chi_j = se_chi_j,
                 omega_j = omega_j, se_omega_j = se_omega_j, pval_omega_j = pval_omega_j)
    })
    res <- do.call(rbind, res_tmp)
    
    eif_taup_uncen <- phi
    taup <- mean(phi)
    eif_taup <- eif_taup_uncen - taup
    var_taup <- (1/N)*mean(eif_taup^2)
    
    eif_L1.1_uncen <- phi_tmp$phi1.1
    EL1.1 <- mean(eif_L1.1_uncen)
    eif_L1.1 <- eif_L1.1_uncen - EL1.1
    var_EL1.1 <- (1/N)*mean(eif_L1.1^2)
    
    eif_L1.0_uncen <- phi_tmp$phi1.0
    EL1.0 <- mean(eif_L1.0_uncen)
    eif_L1.0 <- eif_L1.0_uncen - EL1.0
    var_EL1.0 <- (1/N)*mean(eif_L1.0^2)
    
  } else{
    vals <- lapply(1:folds, function(i){
      pi <- do.call(propML, c(nuisance$prop, list(X = confounders[k != i,], A = A[k != i])))
      if(tauLearner == "S"){
        S1 <- do.call(survML, c(nuisance$surv, list(time = time[k != i], status = 1*(status[k != i] == 1), X = confounders[k != i,], A = treatment[k != i])))
        S2 <- do.call(survML, c(nuisance$comp, list(time = time[k != i], status = 1*(status[k != i] == 2), X = confounders[k != i,], A = treatment[k != i])))
        C <- do.call(survML, c(nuisance$cens, list(time = time[k != i], status = 1*(status[k != i] == 0), X = confounders[k != i,], A = treatment[k != i])))
        phi_tmp <- phiLyl(surv = S1,
                       comp = S2,
                          cens = C,
                          treat = pi,
                          newtimes = time,
                          newstatus = status,
                          newX = confounders,
                          newA = treatment,
                          predtimes = evaltimes)
        phi <- phi_tmp$phi1
        eif_L1.1_uncen <- phi_tmp$phi1.1
        eif_L1.0_uncen <- phi_tmp$phi1.0
        tau_tmp <- phi_tmp$L1
  
        taul_tmp <- lapply(seq_along(l), function(x){
          do.call(pseudoReg, c(nuisance$pseudoRegTau, list(pseudo = tau_tmp[k != i, 1], X = confounders[k != i, lminus[[x]], drop = FALSE])))
        })
        taul <- lapply(seq_along(l), function(x){
          predict(taul_tmp[[x]], newX = confounders[k == i, lminus[[x]], drop = FALSE])
        })
        
        EXl_Xlminus_tmp <- lapply(seq_along(l), function(x){
          do.call(pseudoReg, c(nuisance$pseudoRegX, list(pseudo = confounders[k != i, l[x]], X = confounders[k != i, lminus[[x]], drop = FALSE])))
        })
        EXl_Xlminus <- lapply(seq_along(l), function(x){
          predict(EXl_Xlminus_tmp[[x]], newX = confounders[k == i, lminus[[x]], drop = FALSE])
        })
      } else if(tauLearner == "T"){
        S1.1 <- do.call(survML, c(nuisance$surv, list(time = time[k != i & treatment == 1], status = 1*(status[k != i & treatment == 1] == 1), X = confounders[k != i & treatment == 1,])))
        S1.0 <- do.call(survML, c(nuisance$surv, list(time = time[k != i & treatment == 0], status = 1*(status[k != i & treatment == 0] == 1), X = confounders[k != i & treatment == 0,])))
        S2.1 <- do.call(survML, c(nuisance$surv, list(time = time[k != i & treatment == 1], status = 1*(status[k != i & treatment == 1] == 2), X = confounders[k != i & treatment == 1,])))
        S2.0 <- do.call(survML, c(nuisance$surv, list(time = time[k != i & treatment == 0], status = 1*(status[k != i & treatment == 0] == 2), X = confounders[k != i & treatment == 0,])))
        C.1 <- do.call(survML, c(nuisance$cens, list(time = time[k != i & treatment == 1], status = 1*(status[k != i & treatment == 1] == 0), X = confounders[k != i & treatment == 1,])))
        C.0 <- do.call(survML, c(nuisance$cens, list(time = time[k != i & treatment == 0], status = 1*(status[k != i & treatment == 0] == 0), X = confounders[k != i & treatment == 0,])))
        phi1 <- phiLyl(surv = S1.1,
                        comp = S2.1,
                           cens = C.1,
                           treat = pi,
                           newtimes = time,
                           newstatus = status,
                           newX = confounders,
                           newA = treatment,
                           predtimes = evaltimes)
        phi0 <- phiLyl(surv = S1.0,
                         comp = S2.0,
                           cens = C.0,
                           treat = pi,
                           newtimes = time,
                           newstatus = status,
                           newX = confounders,
                           newA = treatment,
                           predtimes = evaltimes)
        phi <- phi1$phi1.1 - phi0$phi1.0
        eif_L1.1_uncen <- phi1$phi1.1
        eif_L1.0_uncen <- phi0$phi1.0
        tau_tmp <- phi1$L_1.1 - phi0$L_1.0
        taul_tmp <- lapply(seq_along(l), function(x){
          do.call(pseudoReg, c(nuisance$pseudoRegTau, list(pseudo = tau_tmp[k != i, 1], X = confounders[k != i, lminus[[x]], drop = FALSE])))
        })
        taul <- lapply(seq_along(l), function(x){
          predict(taul_tmp[[x]], newX = confounders[k == i, lminus[[x]], drop = FALSE])
        })
        
        EXl_Xlminus_tmp <- lapply(seq_along(l), function(x){
          do.call(pseudoReg, c(nuisance$pseudoRegX, list(pseudo = confounders[k != i, l[x]], X = confounders[k != i, lminus[[x]], drop = FALSE])))
        })
        EXl_Xlminus <- lapply(seq_along(l), function(x){
          predict(EXl_Xlminus_tmp[[x]], newX = confounders[k == i, lminus[[x]], drop = FALSE])
        })
      }
      
      eif_Gamma_j_uncen <- lapply(seq_along(l), function(x){
        (phi[k == i, ] - taul[[x]])*(xl_numeric[[x]][k == i] - EXl_Xlminus[[x]])
      })
      Gamma_j <- lapply(eif_Gamma_j_uncen, mean)
      
      eif_chi_j_uncen <- lapply(seq_along(l), function(x){
        (xl_numeric[[x]][k == i] - EXl_Xlminus[[x]])^2
      })  
      chi_j <- lapply(eif_chi_j_uncen, mean)
      
      eif_taup_uncen <- phi[k == i, ]
      taup <- mean(eif_taup_uncen)
      
      EL1.1 <- mean(eif_L1.1_uncen[k == i, ])
      EL1.0 <- mean(eif_L1.0_uncen[k == i, ])
      
      vals <- list(nk = sum(k == i), eif_Gamma_j_uncen = eif_Gamma_j_uncen,
                   eif_taup_uncen = eif_taup_uncen, eif_chi_j_uncen = eif_chi_j_uncen,
                   Gamma_j = Gamma_j, chi_j = chi_j,
                   taup = taup, EL1.1 = EL1.1, EL1.0 = EL1.0,
                   eif_L1.1_uncen  = eif_L1.1_uncen[k == i, ], eif_L1.0_uncen = eif_L1.0_uncen[k == i, ])
      vals
    })
    
    res_tmp <- lapply(seq_along(l), function(j){
      Gamma_j <- sum(sapply(vals, function(x) x$nk/N * x$Gamma_j[[j]]))
      chi_j <- sum(sapply(vals, function(x) x$nk/N * x$chi_j[[j]]))
      omega_j <- Gamma_j/chi_j
      
      var_Gamma_j <- sum(sapply(vals, function(x) x$nk/N * mean((x$eif_Gamma_j_uncen[[j]] - Gamma_j)^2)))/N
      var_chi_j <- sum(sapply(vals, function(x) x$nk/N * mean((x$eif_chi_j_uncen[[j]] - chi_j)^2)))/N
      var_omega_j <- sum(sapply(vals, function(x) x$nk/N * mean(((x$eif_Gamma_j_uncen[[j]] - Gamma_j - omega_j*(x$eif_chi_j_uncen[[j]] - chi_j))/chi_j)^2)))/N
      
      pval_Gamma_j <- 2*(1 - pnorm(abs(Gamma_j)/sqrt(var_Gamma_j)))
      pval_omega_j <- 2*(1 - pnorm(abs(omega_j)/sqrt(var_omega_j)))
      data.frame(l = names(confounders)[l[j]],
                 Gamma_j = Gamma_j, se_Gamma_j = sqrt(var_Gamma_j), pval_Gamma_j = pval_Gamma_j,
                 chi_j = chi_j, se_chi_j = sqrt(var_chi_j),
                 omega_j, se_omega_j = sqrt(var_omega_j), pval_omega_j = pval_omega_j)
    })
    res <- do.call(rbind, res_tmp)
    
    taup <- sum(sapply(vals, function(x) x$nk/N * x$taup))
    var_taup <- sum(sapply(vals, function(x) x$nk/N * mean((x$eif_taup_uncen - taup)^2)))/N
    
    EL1.1 <- sum(sapply(vals, function(x) x$nk/N * x$EL1.1))
    var_EL1.1 <- sum(sapply(vals, function(x) x$nk/N * mean((x$eif_L1.1_uncen - EL1.1)^2)))/N
    
    EL1.0 <- sum(sapply(vals, function(x) x$nk/N * x$EL1.0))
    var_EL1.0 <- sum(sapply(vals, function(x) x$nk/N * mean((x$eif_L1.0_uncen - EL1.0)^2)))/N
    
    
  }
  
  list(modifiers = res,
       ate = data.frame(ATE = taup, seAte = sqrt(var_taup),
                        EL1.1 = EL1.1, seL1.1 = sqrt(var_EL1.1),
                        EL1.0 = EL1.0, seL1.0 = sqrt(var_EL1.0)),
       treatment_diff = paste0(rev(treatLabs), collapse = " - "))
}

































