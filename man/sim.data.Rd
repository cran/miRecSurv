\name{sim.data}
\docType{data}
\alias{sim.data}
\title{Simulated data set}
\description{
  This data corresponds to a recurrent events simulated cohort using the \code{survsim} package.
}
\usage{sim.data}
\format{A \code{data.frame} with 668 rows and 17 columns, including:
\enumerate{
    \item{nid}{an integer number that identifies the subject.}
    \item{real.episode }{number of the episode corresponding to the real history of the individual.}
    \item{obs.episode }{number of the episode corresponding to the follow-up time of the individual.}
    \item{time }{time until the corresponding event happens (or time to subject drop-out), regarding the beginning of the follow-up time.}
    \item{status }{logical value indicating if the episode corresponds to an event or a drop-out.}
    \item{start }{time at which an episode starts, taking the beginning of follow-up as the origin of the time scale.}
    \item{stop }{time at which an episode ends, taking the beginning of follow-up as the origin of the time scale.}
    \item{time2 }{time until the corresponding event happens (or time to subject drop-out), in calendar time.}
    \item{start2 }{time at which an episode starts, where the time scale is calendar time.}
    \item{stop2 }{time at which an episode ends, where the time scale is calendar time.}
    \item{old }{real value indicating the time that the individual was at risk before the beginning of follow-up.}
    \item{risk.bef }{factor that indicates if an individual was at risk before the beginning of follow-up or not.}
    \item{long }{time not at risk immediately after an episode.}
    \item{z}{Individual heterogeneity.}
    \item{x}{binomial covariate.}
    \item{x.1}{binomial covariate.}
    \item{x.2}{binomial covariate.}
 }
}
\keyword{datasets}
