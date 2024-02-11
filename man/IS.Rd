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
a data frame with the columns corresponding in order to the instruments, the exposures and the outcome
}
}

\value{
It returns a table with the colunms of the coefficient of determination, F statistics, and the conditional F statistics for each exposures selected
}
