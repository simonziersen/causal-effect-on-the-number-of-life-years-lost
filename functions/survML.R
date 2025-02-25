survML <- function(time, status, X, A = NULL, method = "gam", ...){
  require(mgcv)
  require(randomForestSRC)
  obj <- list()
  obj$time <- time
  obj$status <- status
  obj$A <- A
  
  obj$cts <- which(!sapply(X, is.factor)) # based on factor vs. non-factor variables
  if(!is.null(A)){
    obj$x <- data.frame(cbind(X, A = A))
  } else{
    obj$x <- data.frame(X)
  }
  class(obj) <- paste0("survML.", method)
  fit <- survMLfit(obj, ...)
  fit
}

survMLfit <- function(obj, ...){
  UseMethod("survMLfit")
}

survMLfit.survML.coxph <- function(obj, interactions = FALSE, form = NULL){
  require(survival)
  time <- obj$time
  status <- obj$status
  if(!is.null(form)){
    form <- as.formula(paste0("Surv(time,status==1)~", form))
  } else{
    if(!("data.frame" %in% class(obj$x))){
      obj$x <- data.frame(obj$x)
    }
    p <- dim(obj$x)[2]
    if(length(setdiff(1:p, obj$cts)) == 0){
      form <- as.formula(paste0("Surv(time,status==1)~",
                                paste0(colnames(obj$x[, obj$cts, drop = FALSE]), collapse = "+")))
    } else{
      if(interactions){
        form <- as.formula(paste0("Surv(time,status==1)~",
                                  paste0(colnames(obj$x[, obj$cts, drop = FALSE]), collapse = "*"), "*",
                                  paste0(colnames(obj$x[, setdiff(1:p, obj$cts), drop = FALSE]))))
      } else{
        form <- as.formula(paste0("Surv(time,status==1)~",
                                  paste0(colnames(obj$x[, obj$cts, drop = FALSE]), collapse = "+"), "+",
                                  paste0(colnames(obj$x[, setdiff(1:p, obj$cts), drop = FALSE]), collapse = "+")))
      }
    }
  }
  
  
  m <- coxph(form, data = obj$x, x = TRUE)
  obj$m <- m
  obj
}

survMLfit.survML.glmnet <- function(obj, alpha = 0.5){
  require(glmnet)
  time <- obj$time
  status <- obj$status
  m <- cv.glmnet(y = Surv(time,status), x = as.matrix(obj$x), alpha = alpha, family = "cox")
  obj$m <- m
  obj
}

survMLfit.survML.gam <- function(obj, interactions = FALSE, formula = NULL){
  #browser()
  time <- obj$time
  if(!("data.frame" %in% class(obj$x))){
    obj$x <- data.frame(obj$x)
  }
  p <- dim(obj$x)[2]
  dis <- setdiff(1:p, obj$cts)
  if(!is.null(formula)){
    form <- as.formula(paste0("Surv(time,status==1)~", formula))
  } else{
    if(interactions){
      tmp <- NULL
      if(length(dis) == 0){
        form <- paste0("s(", paste0(colnames(obj$x[, obj$cts, drop = FALSE]), collapse = ","), ")")
      } else{
        for(i in 1:length(dis)){
          tmp[i] <- paste0("s(", paste0(colnames(obj$x[, obj$cts, drop = FALSE]), collapse = ","), ",by=", colnames(obj$x[, dis, drop = FALSE]), ")")
        }
        form <- as.formula(paste0(paste0("time~", tmp, collapse = "+"), "+", paste0(colnames(obj$x[, dis, drop = FALSE]))))
      }
      
    } else{
      if(length(dis) == 0){
        form <- as.formula(paste0("time~",
                                  paste0("s(", colnames(obj$x[, obj$cts, drop = FALSE]), ")", collapse = "+")))
      } else{
        form <- as.formula(paste0("time~",
                                  paste0("s(", colnames(obj$x[, obj$cts, drop = FALSE]), ")", collapse = "+"), "+",
                                  paste0(colnames(obj$x[, dis, drop = FALSE]), collapse = "+")))
      }
      
    }
  }
  
  m <- mgcv::gam(form, family = cox.ph(), data = obj$x, weights = obj$status)
  obj$m <- m
  obj
}

survMLfit.survML.ranger <- function(obj, ...){
  require(ranger)
  time <- obj$time
  data_tmp <- data.frame(time = obj$time, status = obj$status, obj$x)
  m <- ranger(Surv(time, status)~., data = data_tmp)
  obj$m <- m
  obj
}

survMLfit.survML.rfsrc <- function(obj, ...){
  require(randomForestSRC)
  time <- obj$time
  data_tmp <- data.frame(time = obj$time, status = obj$status, obj$x)
  m <- rfsrc(Surv(time, status)~., data = data_tmp, ...)
  obj$m <- m
  obj
}

survMLfit.survML.rfsrcTune <- function(obj, ...){
  require(randomForestSRC)
  time <- obj$time
  data_tmp <- data.frame(time = obj$time, status = obj$status, obj$x)
  m_tmp <- tune(Surv(time, status)~., data = data_tmp, doBest = TRUE, ...)
  m <- rfsrc(Surv(time, status)~., data = data_tmp, nodesize = m_tmp$optimal["nodesize"], mtry = m_tmp$optimal["nodesize"], ...)
  obj$m <- m
  obj
}

predict.survML.coxph <- function(obj, newtimes, newX, newA = NULL){
  #browser()
  N <- dim(newX)[1]
  tl <- length(newtimes)
  if(!is.null(obj$A) & is.null(newA)){
    stop("model fitted with A, but no new A is provided. newA should be a binary vector of 0's and 1's")
  } else if(!is.null(newA)){
    #x <- cbind(newX, A = as.factor(newA))
    x <- cbind(newX, A = newA)
  } else{
    x <- newX
  }
  surv <- predictCox(obj$m, newdata = x, type = "survival", times = newtimes)$survival
  cumhaz <- -log(surv)
  list(surv = surv, cumhaz = cumhaz)
}

predict.survML.gam <- function(obj, newtimes, newX, newA = NULL){
  N <- dim(newX)[1]
  tl <- length(newtimes)
  if(!is.null(obj$A) & is.null(newA)){
    stop("model fitted with A, but no new A is provided. newA should be a binary vector of 0's and 1's")
  } else if(!is.null(newA)){
    x <- cbind(newX, A = newA)
  } else{
    x <- newX
  }
  predData <- data.frame(time = rep(newtimes, each = N), x[rep(1:N, tl),])
  tmp <- predict(obj$m, newdata = predData, type = "response", se = FALSE)
  surv <- matrix(tmp, nrow = N, ncol = tl)
  cumhaz <- -log(surv)
  list(surv = surv, cumhaz = cumhaz)
}

predict.survML.ranger <- function(obj, newtimes, newX, newA = NULL){
  
  N <- dim(newX)[1]
  tl <- length(newtimes)
  if(!is.null(obj$A) & is.null(newA)){
    stop("model fitted with A, but no new A is provided. newA should be a binary vector of 0's and 1's")
  } else if(!is.null(newA)){
    x <- cbind(newX, A = newA)
  } else{
    x <- newX
  }
  tmp <- predict(obj$m, data = x)
  minTime <- min(obj$m$unique.death.times)
  surv <- sapply(newtimes, function(t0){
    if(t0 < minTime){
      return(rep(1, N))
    } 
    ind <- max(which(t0 >= obj$m$unique.death.times))
    tmp$survival[, ind]
  })
  cumhaz <- -log(surv)
  list(surv = surv, cumhaz = cumhaz)
}

predict.survML.rfsrc <- function(obj, newtimes, newX, newA = NULL){
  N <- dim(newX)[1]
  tl <- length(newtimes)
  
  if(!is.null(obj$A) & is.null(newA)){
    stop("model fitted with A, but no new A is provided. newA should be a binary vector of 0's and 1's")
  } else if(!is.null(newA)){
    x <- cbind(newX, A = newA)
  } else{
    x <- newX
  }
  tmp <- predict(obj$m, newdata = x)
  minTime <- min(obj$m$time.interest)
  surv <- sapply(newtimes, function(t0){
    if(t0 < minTime){
      return(rep(1, N))
    } 
    ind <- max(which(t0 >= sort(obj$m$time.interest)))
    tmp$survival[, ind]
  })
  cumhaz <- -log(surv)
  list(surv = surv, cumhaz = cumhaz)
}

predict.survML.rfsrcTune <- function(obj, newtimes, newX, newA = NULL){
  
  N <- dim(newX)[1]
  tl <- length(newtimes)
  
  if(!is.null(obj$A) & is.null(newA)){
    stop("model fitted with A, but no new A is provided. newA should be a binary vector of 0's and 1's")
  } else if(!is.null(newA)){
    x <- cbind(newX, A = newA)
  } else{
    x <- newX
  }
  tmp <- predict(obj$m, newdata = x)
  minTime <- min(obj$m$time.interest)
  surv <- sapply(newtimes, function(t0){
    if(t0 < minTime){
      return(rep(1, N))
    } 
    ind <- max(which(t0 >= sort(obj$m$time.interest)))
    tmp$survival[, ind]
  })
  cumhaz <- -log(surv)
  list(surv = surv, cumhaz = cumhaz)
}

predict.survML.glmnet <- function(obj, newtimes, newX, newA = NULL){
  N <- dim(newX)[1]
  tl <- length(newtimes)
  if(!is.null(obj$A) & is.null(newA)){
    stop("model fitted with A, but no new A is provided. newA should be a binary vector of 0's and 1's")
  } else if(!is.null(newA)){
    x <- as.matrix(cbind(newX, A = newA))
  } else{
    x <- as.matrix(newX)
  }
  tmp <- survfit(obj$m, s = obj$m$lambda.min, x = as.matrix(obj$x), y = Surv(obj$time, obj$status), newx = x)
  surv <- sapply(newtimes, function(t0){
    ind <- max(which(t0 >= sort(unique(obj$time))))
    t(tmp$surv)[, ind]
  })
  cumhaz <- -log(surv)
  list(surv = surv, cumhaz = cumhaz)
}
  
  
  