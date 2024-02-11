\name{plotEifun}
\alias{plotEifun}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
draw the eigenfunction plot based on a FPCA result
}

\usage{
plotEifun(res)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{res}{
a list result of FPCA. \code{res} is usually derived by \code{FPCA} from the package fdapace. \code{res} must contains the following elements:
\code{res$lambda}: the vector correpsonding the fraction of variance explained by each principal component (i.e. eigenvalues).
\code{res$workGrid}: a vector indicates the working time Grid (such time points are mainly used for visulization purpose)
\code{res$phi}: a matrix, where each column indicates the values of one eigenfunction over \code{res$workGrid}
}
}

\value{
a gg-plot where each curve corresponds to one eigenfunction with the corresponding the eigenvalue (i.e. the fraction of variance explained by this eigenfunction).

}



