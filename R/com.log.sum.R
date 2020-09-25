com.log.sum <- function (x, y)
{
  if (x == -Inf) {
    return(y)
  }
  else if (y == -Inf) {
    return(x)
  }
  else if (x > y) {
    return(x + log(1 + exp(y - x)))
  }
  else {
    return(y + log(1 + exp(x - y)))
  }
}
