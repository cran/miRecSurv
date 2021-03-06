\encoding{utf8}
\name{recEvFit}
\alias{recEvFit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Left-censored recurrent events survival models
}
\description{
The function allows the user to fit recurrent events survival models.
}
\usage{
recEvFit(formula, data, id, prevEp, riskBef, oldInd,
         frailty=FALSE, m=5, seed=NA, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
a formula object, with the response on the left of a \code{~} operator, and the terms on the right. The response must be a survival object as returned by the Surv function.
}
  \item{data}{
a data.frame in which to interpret the variables named in the formula.
}
  \item{id}{
subject identifier.
}
  \item{prevEp}{
known previous episodes.
}
  \item{riskBef}{
indicator for new individual in the cohort (\code{riskBef==FALSE}) or subject who was at risk before the start of follow-up (\code{riskBef==TRUE}).
}
  \item{oldInd}{
time an individual has been at risk prior to the follow-up.
}
  \item{frailty}{
should the model include a frailty term. Defaults to \code{FALSE}.
}
  \item{m}{
number of multiple imputations. The default is \code{m=5.}
}
  \item{seed}{
an integer that is used as argument by the \code{set.seed} function for offsetting the random number generator. Default is to leave the random number generator alone.
}
\item{...}{
extra arguments to pass to \code{coxph}.
}
}

\value{
A list with seven elements:
 \item{fit}{
a list with all the \code{coxph} objects fitted for each imputed dataset.
}
 \item{coeff}{
a list with the vectors of coefficients from the models fitted to each imputed dataset
}
 \item{loglik}{
a list with the loglikelihood for each model fitted.
}
 \item{vcov}{
a list with the variance-covariance matrices for the parameters fitted for each of the imputed datasets.
}
 \item{AIC}{
a list with the AIC of each of the models fitted.
}
\item{CMP}{
summary tables of the fitted COMPoisson models used for imputing missing values
}
 \item{data.impute}{
the original dataset with the multiple imputed variables as final columns.
}
}

\author{
David Moriña (University of Barcelona), Gilma Hernández-Herrera (Universidad de Antioquía), Albert Navarro (Universitat Autònoma de Barcelona)

Mantainer: David Moriña <dmorina@ub.edu>
}

\seealso{
\code{\link{miRecSurv-package}}
}
\examples{
data(sim.data)
fit <- recEvFit(Surv(start2, stop2, status)~x+x.1+x.2, data=sim.data,
                id="nid", prevEp = "obs.episode",
                riskBef = "risk.bef", oldInd = "old", frailty=FALSE, m=5)
summary(fit)
}

\keyword{ radir }
