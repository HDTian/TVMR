\name{getX}
\alias{getX}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
simulate a data for the instrument and exposure data with time-varying variable information
}
\description{
simulate a data for the instrument and exposure data with time-varying variable information
}
\usage{
getX(  N=10000,
       J=30,
       ZXmodel='A',
       MGX_used=NA,
       nSparse=10
       )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{
a integer indicates the sample size
}
 \item{J}{
a integer indicates the number of genetic instrument
}
\item{ZXmodel}{
a character indicates the form of the instrument-exposure model. Can be \code{'A'}, \code{'B'}, \code{'C'}, or \code{'D'}
}
\item{MGX_used}{
a matrix indicates the user-defined instrumental effect on the exposure. The default value is \code{NA}, corresponding to the case that the genetic effect is randomly valued.
}
\item{nSparse}{
a integer indicates the number of the random measured timepoints for each individual.
}

}

\value{
\code{getX} returns a list, which contains the information used for reproducing simulation, and the data matrix.
}
\examples{
RES<-getX(J=30,ZXmodel='D')
}

