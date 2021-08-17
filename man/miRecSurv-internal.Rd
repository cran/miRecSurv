\encoding{utf8}
\name{miRecSurv-internal}
\alias{miRecSurv-internal}
\alias{com.compute.log.z}
\alias{com.log.density}
\alias{com.log.factorial}
\alias{com.log.difference}
\alias{com.log.sum}
\alias{rcom}
\alias{accum.sample}
\alias{summary.recEvFit}
\alias{print.summary.recEvFit}

\title{Internal miRecSurv functions}
\description{Internal miRecSurv functions}
\usage{
com.compute.log.z(lambda, nu, log.error = 0.001)
com.log.density(x, lambda, nu, log.z = NULL)
com.log.difference(x, y)
com.log.factorial(x)
com.log.sum(x, y)
rcom(n, lambda, nu, log.z = NULL)
accum.sample(data, id, status, covars, riskBef, oldInd) 
\method{summary}{recEvFit}(object, ...)
\method{print.summary}{recEvFit}(x, ...)
}

\details{These functions are not to be called by the user.}
\seealso{
\code{\link{miRecSurv-package}}, \code{\link{recEvFit}}
}
\keyword{internal}
