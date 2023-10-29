\name{getDat}
\alias{getDat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
simulate a dat with time-varying variable information
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getDat(N=10000, #sample size
       J=30,# number of SNPs (three-valued SNP)
       ZXmodel='A',
       MGX_used=NA, #if user define MGX, then will use this MGX; this is important for two-sample case where certain SNP should have the same effect
       nSparse=10
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
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\code{getDat} returns a list, which contains the information used for reproducing simulation, and the data matrix.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
RES<-getX(J=30,ZXmodel='D')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
