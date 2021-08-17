com.log.difference <- function (x, y)
{
  if (x == -Inf) {
    return(NaN)
  }
  else if (y == -Inf) {
    return(x)
  }
  else if (x > y) {
    return(x + log(1 - exp(y - x)))
  }
  else {
    return(NaN)
  }
}
