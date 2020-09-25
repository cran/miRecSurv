com.log.density <- function (x, lambda, nu, log.z = NULL)
{
  if (lambda < 0 || nu < 0)
    stop("Invalid arguments, only defined for lambda >= 0, nu >= 0")
  if (x < 0 || x != floor(x))
    return(0)
  if (is.null(log.z)) {
    log.z = com.compute.log.z(lambda, nu)
  }
  return((x * log(lambda) - nu * com.log.factorial(x)) - log.z)
}
