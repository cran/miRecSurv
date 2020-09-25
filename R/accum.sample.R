accum.sample <- function (data, id, status, covars, riskBef, oldInd) 
{
  data2 <- split(data, as.factor(data[, id]))
  for (i in 1:length(data2)) {
    if (dim(data2[[i]])[1] != 0) {
      data2[[i]]$obs.ep.accum[1] <- sum(data2[[i]][, status])
      data2[[i]]$time.accum[1] <- sum(na.omit(data2[[i]]$ft))
    }
    data2[[i]] <- data2[[i]][1:1, ]
  }
  data2 <- do.call(rbind, data2)
  data2 <- data2[, c(id, covars, "obs.ep.accum", "time.accum", riskBef, oldInd)]
  return(data2)
}