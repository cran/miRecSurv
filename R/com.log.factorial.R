com.log.factorial <- function (x)
{
  if (is.vector(x) && length(x) > 1) {
    for (i in 1:length(x)) x[i] = com.log.factorial(x[i])
    return(x)
  }
  else if (is.numeric(x)) {
    if (x == 0) {
      x = 1
    }
    return(sum(log(seq(from = 1, to = x, by = 1))))
  }
  else {
    stop("x must be a vector or number.")
  }
}
