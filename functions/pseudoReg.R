pseudoReg <- function(pseudo, X, method = "gam", ...){
  obj <- list()
  obj$x <- X
  obj$pseudo <- pseudo
  #obj$cts <- 1:dim(X)[2] # just for example - has to be rewritten to work on generic covariates
  obj$cts <- which(!sapply(X, is.factor)) # based on factor vs. non-factor variables
  class(obj) <- paste0("pseudoReg.", method)
  fit <- pseudoRegfit(obj, ...)
  fit
}

pseudoRegfit <- function(obj, ...){
  UseMethod("pseudoRegfit")
}

pseudoRegfit.pseudoReg.gam <- function(obj, interactions = FALSE){
  require(mgcv)
  p <- dim(obj$x)[2]
  dis <- setdiff(1:p, obj$cts)
  pseudo <- obj$pseudo
  if(interactions){
    tmp <- NULL
    if(length(dis) == 0){
      tmp <- paste0("s(", paste0(colnames(obj$x[, obj$cts, drop = FALSE]), collapse = ","), ")")
    } else{
      for(i in 1:length(dis)){
        tmp[i] <- paste0("s(", paste0(colnames(obj$x[, obj$cts, drop = FALSE]), collapse = ","), ",by=", colnames(obj$x[, dis, drop = FALSE]), ")")
      }
    }
    form <- as.formula(paste0("pseudo~", tmp, collapse = "+"))
  } else{
    if(length(dis) == 0){
      form <- as.formula(paste0("pseudo~",
                                paste0("s(", colnames(obj$x[, obj$cts, drop = FALSE]), ")", collapse = "+")))
    } else{
      form <- as.formula(paste0("pseudo~",
                                paste0("s(", colnames(obj$x[, obj$cts, drop = FALSE]), ")", collapse = "+"),
                                "+",
                                paste0(colnames(obj$x[, dis, drop = FALSE]), collapse = "+")))
    }
  }
  m <- gam(form, data = obj$x)
  obj$m <- m
  obj
}

pseudoRegfit.pseudoReg.ranger <- function(obj, ...){
  require(ranger)
  tmp <- as.numeric(obj$pseudo)
  data_tmp <- cbind(tmp, obj$x)
  m <- ranger(tmp ~ ., data = data_tmp)
  obj$m <- m
  obj
}

pseudoRegfit.pseudoReg.rfsrc <- function(obj, ...){
  require(randomForestSRC)
  tmp <- obj$pseudo
  data_tmp <- cbind(tmp, obj$x)
  m <- rfsrc(tmp ~ ., data = data_tmp)
  obj$m <- m
  obj
}

predict.pseudoReg.gam <- function(obj, newX){
  as.numeric(predict(obj$m, newdata = newX, se = FALSE))
}

predict.pseudoReg.ranger <- function(obj, newX){
  as.numeric(predict(obj$m, data = newX)$predictions)
}

predict.pseudoReg.rfsrc <- function(obj, newX){
  if(is.numeric(obj$pseudo)){
    return(as.numeric(predict(obj$m, newdata = newX)$predicted))
  } else if(is.factor(obj$pseudo)){
    return(as.numeric(predict(obj$m, newdata = newX)$predicted[,2]))
  }
}





