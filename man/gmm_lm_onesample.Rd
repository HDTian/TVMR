\name{gmm_lm_onesample}
\alias{gmm_lm_onesample}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multivariable GMM and Kleibergen's Lagrange multplier (LM) statistics with one-sample individual-level data
}

\usage{
gmm_lm_onesample(X,
                 Y,
                 Z,
                 beta0=NA
                 )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{Z}{n x J instrument matrix}
\item{X}{n x K exposure matrix}
\item{Y}{n x 1 outcome vector}
\item{beta0}{the tested null of the causal parameter value (for LM test only)}
}

\value{
a result list, containing
\item{gmm_est}{K vector of causal effect estimates using GMM}
\item{gmm_se}{K vector of standard errors corresponding to \code{gmm_est}}
\item{variance_matrix}{K x K variance matrix corresponding to \code{gmm_est}}
\item{gmm_pval}{K vector of p-values corresponding to gmm_est}
\item{Q_stat}{overidentification test statistic}
\item{Q_pval}{overidentification test p-value}
\item{lm_stat}{the LM statistic value}
\item{lm_pval}{p-value of the LM test of the null hypothesis H0: beta=beta0}

}

\author{
Haodong Tian, Ashish Patel
}


