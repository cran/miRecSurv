recEvFit <- function(formula, data, id, prevEp, riskBef, oldInd,
                     frailty=FALSE, m=5, seed=NA)
{
  complete.eprev <- function(data, id_n, id, var){
    Daux <- data[data[, id]==id_n,]
    n <- dim(Daux)[1]
    eval(parse(text=paste0("aux <- seq_len(n)+(", "Daux[1,", var, "]-1)")))
    return(aux)
  }
  
  if (!is.na(seed)) set.seed(seed)
  Call <- match.call()
  if (missing(data))
    data <- environment(formula)
  mf <- match.call()
  mm <- match(c("formula", "data"), names(mf), 0)
  mf <- mf[c(1, mm)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  X <- model.matrix(mt, mf)
  Y <- model.response(mf, "numeric")
  
  terms <- as.character(formula)
  sta <- ifelse(grepl("-", terms[2])==FALSE, stri_extract_all_words(terms[2])[[1]][2],
                stri_extract_all_words(terms[2])[[1]][3])
  sto <- ifelse(grepl("-", terms[2])==FALSE, stri_extract_all_words(terms[2])[[1]][3],
                stri_extract_all_words(terms[2])[[1]][2])
  data$ft <- data[, sto]-data[, sta]
  status <- stri_extract_all_words(terms[2])[[1]][4]
  sample <- data
  covars <- stri_split_coll(terms[3][[1]], "+")[[1]]
  for (i in 1:length(covars))
  {
    covars[i] <- gsub(" ", "", covars[i])
  }
  data[, oldInd][is.na(data[, oldInd])] <- 0
  data$oldInd2 <- ifelse(data[, oldInd] < 0, (-1)*data[, oldInd], data[, oldInd])
  accum_sample <- accum.sample(sample, id, status, covars, riskBef, oldInd)
  mod <- glm.cmp(formula.lambda=accum_sample$obs.ep.accum~as.matrix(accum_sample[,covars])+offset(log(accum_sample$time.accum)), data=accum_sample)  
  beta <- unlist(coef(mod))
  rv <- t(chol(vcov(mod)))
  if (!is.na(seed)) set.seed(seed)
  b.star <- beta + rv %*% rnorm(ncol(rv))
  if (b.star[length(b.star)] < 0.05) b.star[length(b.star)] <- 0.05
  if (b.star[length(b.star)] > 500)  b.star[length(b.star)] <- 500
  lambda <- ifelse(data[, riskBef]==TRUE, exp(X %*% b.star[1:(length(b.star)-1)])*data$oldInd2, 0)
  
  for (i in 1:m)
  {
    for (j in 1:nrow(data))
    {
      if(data[j, riskBef]==TRUE)
      {
        eval(parse(text=paste0("data$EprevCOMPoiss",i,"[",j,"] <- rcom(1, lambda[",j,"], nu=exp(b.star[", length(b.star), "]))")))
      }else{
        eval(parse(text=paste0("data$EprevCOMPoiss",i,"[",j,"] <- 0")))
      }
    }
  }
  id_2 <- which(names(data)==id)
  rr <- unique(data[, id_2])
  for (i in 1:m)
  {
    col  <- which(colnames(data)==paste0("EprevCOMPoiss", i))
    k <- lapply(unique(data[, id_2]), 
                function(k) complete.eprev(data, k, id, col))
    eval(parse(text=paste0("data$EprevCOMPoissDef", i, "<- unlist(k)"))) 
  }
  for (i in 1:m)
  {
    eval(parse(text=paste0("data$EprevCOMPoiss", i, "<-NULL")))
  }
  if (frailty==TRUE) d <- dim(X)[2]+1
  if (frailty==FALSE) d <- dim(X)[2]
  p_vals  <- matrix(nrow=d, ncol=m)
  coefs   <- matrix(nrow=d, ncol=m)
  expcoef <- matrix(nrow=d, ncol=m)
  sds     <- matrix(nrow=d, ncol=m)
  
  if (frailty==TRUE)
  {
    for (i in 1:m)
    {
      eval(parse(text=paste0("modPWP <-
      coxph(", deparse(formula),
                "+strata(as.factor(100000*", riskBef, "+EprevCOMPoissDef", i, "))+oldInd2+
                 frailty(", id, "), data=data)")))
      coefs[, i]   <- summary(modPWP)$coefficients[, 1]
      expcoef[, i] <- exp(summary(modPWP)$coefficients[, 1])
      sds[, i]     <- summary(modPWP)$coefficients[, 2]
      p_vals[, i]  <- summary(modPWP)$coefficients[, 6]
    }
  }else{
    for (i in 1:m)
    {
      eval(parse(text=paste0("modPWP <-
      coxph(", deparse(formula),
                  "+strata(as.factor(100000*", riskBef, "+EprevCOMPoissDef", i, "))+oldInd2, data=data)")))
      coefs[, i]   <- summary(modPWP)$coefficients[, 1]
      expcoef[, i] <- summary(modPWP)$coefficients[, 2]
      sds[, i]     <- summary(modPWP)$coefficients[, 3]
      p_vals[, i]  <- summary(modPWP)$coefficients[, 5]
    }
  }
  data$oldInd2 <- NULL
  data$ft <- NULL
  gcoef <- rowMeans(coefs)
  expco <- exp(gcoef)
  wvar  <- (1/m)*rowSums(sds^2)
  bvar  <- (1/(m-1))*rowSums((coefs-gcoef)^2)
  gsd   <- sqrt(wvar+(1+1/m)*bvar)
  pvalsf<- rowMedians(p_vals)
  res <- list()
  res[[1]] <- data
  res[[2]] <- cbind(gcoef, expco, gsd, pvalsf)
  colnames(res[[2]]) <- c("coef", "exp(coef)", "se(coef)", "Pr(>|z|)")
  if (frailty==TRUE) rownames(res[[2]]) <- c(colnames(X)[2:length(colnames(X))], oldInd, "frailty")
  if (frailty==FALSE) rownames(res[[2]]) <- c(colnames(X)[2:length(colnames(X))], oldInd)
  class(res) <- "recEvFit"
  return(res)
}