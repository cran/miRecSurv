rcom <- function (n, lambda, nu, log.z = NULL)
{
  if (lambda < 0 || nu < 0)
    stop("Invalid arguments, only defined for lambda >= 0, nu >= 0")
  if (is.null(log.z))
    log.z = com.compute.log.z(lambda, nu)
  r = NULL
  for (i in 1:n) {
    log.prob = log(runif(1))
    j = 0
    while (1) {
      new.log.prob = com.log.difference(log.prob, com.log.density(j,
                                                                  lambda, nu, log.z))
      if (is.nan(new.log.prob))
        break
      log.prob = new.log.prob
      j = j + 1
    }
    r = c(r, j)
  }
  return(r)
}
