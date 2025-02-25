propML <- function(A, X, method = "gam", ...){
  obj <- list()
  obj$A <- A
  obj$x <- X
  #obj$cts <- 1:dim(X)[2] # just for example - has to be rewritten to work on generic covariates
  obj$cts <- which(!sapply(X, is.factor)) # based on factor vs. non-factor variables
  class(obj) <- paste0("propML.", method)
  fit <- propMLfit(obj, ...)
  fit
}

propMLfit <- function(obj, ...){
  UseMethod("propMLfit")
}

propMLfit.propML.glm <- function(obj, interactions = FALSE){
  A <- obj$A
  if(interactions){
    form <- as.formula(paste0("A~", paste0(colnames(obj$x), collapse = "*")))
  } else{
    form <- as.formula(paste0("A~", paste0(colnames(obj$x), collapse = "+")))
  }
  m <- glm(form, family = "binomial", data = obj$x)
  obj$m <- m
  obj
}

propMLfit.propML.gam <- function(obj, interactions = FALSE){
  p <- dim(obj$x)[2]
  dis <- setdiff(1:p, obj$cts)
  A <- obj$A
  if(interactions){
    tmp <- NULL
    if(length(dis) == 0){
      tmp <- paste0("s(", paste0(colnames(obj$x[, obj$cts, drop = FALSE]), collapse = ","), ")")
    }
    for(i in 1:length(dis)){
      tmp[i] <- paste0("s(", paste0(colnames(obj$x[, obj$cts, drop = FALSE]), collapse = ","), ",by=", colnames(obj$x[, dis, drop = FALSE]), ")")
    }
    form <- as.formula(paste0("A~", tmp, collapse = "+"))
  } else{
    if(length(dis) == 0){
      form <- as.formula(paste0("A~",
                                paste0("s(", colnames(obj$x[, obj$cts, drop = FALSE]), ")", collapse = "+")))
    } else{
      form <- as.formula(paste0("A~",
                                paste0("s(", colnames(obj$x[, obj$cts, drop = FALSE]), ")", collapse = "+"), "+",
                                paste0(colnames(obj$x[, dis, drop = FALSE]))))
    }
  }
  m <- gam(form, family = "binomial", data = obj$x)
  obj$m <- m
  obj
}

propMLfit.propML.ranger <- function(obj){
  A <- as.numeric(obj$A)
  data_tmp <- cbind(A, obj$x)
  m <- ranger(A ~ ., data = data_tmp)
  obj$m <- m
  obj
}

propMLfit.propML.rfsrc <- function(obj, ...){
  require(randomForestSRC)
  if(!is.factor(obj$A))
    stop("treatment must be a factor")
  data_tmp <- cbind(A = obj$A, obj$x)
  m <- rfsrc(A ~ ., data = data_tmp, ...)
  obj$m <- m
  obj
}

predict.propML.glm <- function(obj, newX){
  as.numeric(predict(obj$m, newdata = newX, type = "response", se = FALSE))
}

predict.propML.gam <- function(obj, newX){
  as.numeric(predict(obj$m, newdata = newX, type = "response", se = FALSE))
}

predict.propML.ranger <- function(obj, newX){
  as.numeric(predict(obj$m, data = newX)$predictions)
}

predict.propML.rfsrc <- function(obj, newX){
  as.numeric(predict(obj$m, newdata = newX)$predicted[,2])
}