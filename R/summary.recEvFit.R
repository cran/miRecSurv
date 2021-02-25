summary.recEvFit <- function(object, ...)
{
  if (class(object)!="recEvFit") stop("Object should be of class recEcvFit")
  call   <- attr(object, "Call")
  m      <- attr(object, "m")
  X      <- attr(object, "X")
  coefs  <- do.call(rbind, object$coeff)
  aux1   <- lapply(object$vcov, diag)
  aux2   <- lapply(aux1, sqrt)
  sds    <- do.call(rbind, aux2)
  gcoef  <- colMeans(coefs)
  expco  <- exp(gcoef)
  wvar   <- (1/m)*colSums(sds^2)
  bvar   <- (1/(m-1))*colSums((coefs-gcoef)^2)
  gsd    <- sqrt(wvar+(1+1/m)*bvar)
  aux3   <- lapply(object$fit, summary)
  aux4   <- lapply(aux3, coefficients)
  p_name <- ifelse(attr(object, "frailty")==TRUE, "p", "Pr(>|z|)")
  aux5   <- lapply(aux4, "[", , p_name)
  pvals  <- do.call(rbind, aux5)
  pvalsf <- colMedians(pvals)
  aux6   <- ifelse(attr(object, "frailty")==FALSE, lapply(aux4, "[", , "z"), lapply(aux4, "[", , "Chisq"))
  u_aux6 <- unlist(aux6)
  if (attr(object, "frailty")==TRUE) 
  {
    gcoef  <- c(gcoef, NA)
    expco  <- c(expco, NA)
    gsd    <- c(gsd, NA)
  }
  res    <- cbind(gcoef, expco, gsd, u_aux6, pvalsf)
  colnames(res) <- c("coef", "exp(coef)", "se(coef)", colnames(coefficients(summary(object$fit[[1]])))[4], "Pr(>|z|)")
  if (attr(object, "frailty")==TRUE) rownames(res) <- c(colnames(X)[2:length(colnames(attr(object, "X")))], attr(object, "oldInd"), "frailty")
  if (attr(object, "frailty")==FALSE) rownames(res) <- c(colnames(X)[2:length(colnames(attr(object, "X")))], attr(object, "oldInd"))
  aux7   <- do.call(rbind, object$AIC)
  aic <- mean(aux7)
  imp_table <- summary(object$CMP)$DF
  rownames(imp_table) <- c(colnames(X), "S:(Intercept)")
  imp_table <- data.matrix(cbind(imp_table[, 1:3], as.numeric(imp_table[, 4])))
  colnames(imp_table)[4] <- "Pr(>|z|)"
  ans <- list(coefficients = res, aic=aic, call=call, CMP=imp_table)
  
  class(ans) <- "summary.recEvFit"
  return(ans)
}