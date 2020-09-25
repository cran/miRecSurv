\encoding{utf8}
\name{miRecSurv-package}
\alias{miRecSurv-package}
\alias{miRecSurv}
\docType{package}
\title{
Left-Censored Recurrent Events Survival Models
}
\description{
Fitting recurrent events survival models for left-censored data with multiple imputation of the number of previous episodes.
    See Hernández-Herrera, G, Moriña, D and Navarro, A (2020) <arXiv:2007.15031>.
}
\details{
\tabular{ll}{
Package: \tab miRecSurv\cr
Type: \tab Package\cr
Version: \tab 0.0.1\cr
Date: \tab 2020-9-23\cr
License: \tab GPL version 2 or newer\cr
LazyLoad: \tab yes\cr
}

}
\author{
David Moriña (University of Barcelona), Gilma Hernández-Herrera (Universidad de Antioquía), Albert Navarro (Universitat Autònoma de Barcelona)

Mantainer: David Moriña <dmorina@ub.edu>
}
\keyword{ package }
\seealso{
\code{\link{recEvFit}}
}
\examples{
data(sim.data)
fit <- recEvFit(Surv(start2, stop2, status)~x+x.1+x.2, data=sim.data,
                id="nid", prevEp = "obs.episode",
                riskBef = "risk.bef", oldInd = "old", frailty=FALSE, m=5)
fit[[2]]  ### Fit summary table
head(fit[[1]]) ### Data with imputed number of previous episodes
}
