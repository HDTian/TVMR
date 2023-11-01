\name{IS}
\alias{IS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Instrumental strength with the conditional F statistics
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
IS(J, K, timepoints, datafull)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{J}{
a integer indicates the number of instruments
}
  \item{K}{
a integer indicates the number of exposures
}
  \item{timepoints}{
a vector indicates the index of the exposures used for calculating the conditional F statistics
}
  \item{datafull}{
a data frame with the colums corresponding in order to the instruments, the exposures and the outcome
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
It returns a table with the colunms of the coefficient of determination, F statistics, and the conditional F statistics for each exposures selected
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
