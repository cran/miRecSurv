com.compute.log.z <- function (lambda, nu, log.error = 0.001)
{
  if (lambda < 0 || nu < 0)
    stop("Invalid arguments, only defined for lambda >= 0, nu >= 0")
  z = -Inf
  z.last = 0
  j = 0
  while (abs(z - z.last) > log.error) {
    z.last = z
    z = com.log.sum(z, j * log(lambda) - nu * com.log.factorial(j))
    j = j + 1
  }
  return(z)
}
