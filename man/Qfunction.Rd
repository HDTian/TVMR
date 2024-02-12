\name{Qfunction}
\alias{Qfunction}
\title{
Q function
}
\description{
obtain the Q statistic value and alos return the inference results based on Q statistic (e.g. weak-IV-robust estimation), based on an inputted effect parameter.
}
\usage{
Qfunction(v,
          by,
          byse,
          B,
          BX,
          Sigma,
          Gam
          )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{v}{a vector indicates the initial value of the effect parameter; length is L.}
\item{by}{a vector indicates the genetic associations with the outcome; length is J.}
\item{byse}{a vector indicates the standard error of \code{by}; length is J.}
\item{B}{a transforming matrix; dimision is K x L. If using the full eigenfunction as the basisfunction, B is just the identity matrix I.}
\item{BX}{a matrix corresponds to the genetic association with the (transformed) exposure; dimension is J x K. }
\item{Sigma}{a array (K*K*J), each surface of which is the covariance matrix of the association estimators \code{BX[j,]}. }
\item{Gam}{a Gamma matrix, each row of which indicates the covariance of the estimated genetic association with the outcome and the estimated genetic assocation with the (transformed) exposure; dimension is J x K.}
}

\value{
return a list, containing
\item{original_input_v}{the original inputted value of the effect parameter.}
\item{inference_results}{the Q-based (e.g. weak-IV-robust) estimation results, containing the parameter estimate and its standard errors.}
\item{Est}{the estimated effext parameter.}
\item{Estvar}{the variance matrix of \code{Est}}
\item{Qvalue}{the Q satistic value}
}



