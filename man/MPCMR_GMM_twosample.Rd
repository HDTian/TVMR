\name{MPCMR_GMM_twosample}
\alias{MPCMR_GMM_twosample}
\title{
Two-sample time-varying MR fitting
}
\description{
Multiple-principal-component Mendelian randomization fitting with GMM methods. This function is desinged for the two-sample setting.
}
\usage{
MPCMR_GMM_twosample(Gmatrix,
                    res,
                    by_used,
                    sy_used,
                    ny_used,
                    nPC=NA,
                    nL=NA,
                    eigenfit=TRUE,
                    polyfit=TRUE,
                    LMCI=TRUE,
                    LMCI2=TRUE,
                    nLM=20,
                    Parallel=TRUE,
                    cores_used=NA,
                    XYmodel=NA
                    )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Gmatrix}{
a matrix indicates the instrument information
}
  \item{res}{
a list result of FPCA. \code{res} is usually derived by \code{FPCA} from the package fdapace. \code{res} must contains the following elements:
\code{res$cumFVE}: a vector corresponds to the cumulative fraction-of-variance-explained (FVE).
\code{res$xiEst}: a matrix contains the principal components values for each individual, where the rows correspond to the individuals and columns correspond to the peincipal compoemnts.
\code{res$workGrid}: a vector indicates the working time Grid (such time points are mainly used for visulization purpose)
\code{res$phi}: a matrix, where each column indicates the values of one eigenfunction over \code{res$workGrid}
}
  \item{by_used}{
a vector indicates the estimated genetic association with the outcome. The order of the genetic variants in \code{by_used} should be consistent with that order in \code{Gmatrix}
}
  \item{sy_used}{
a vector indicates the standard errors of the estimated genetic association with the outcome
}
  \item{ny_used}{
a integer indicates the smaple size of the outcome data
}
  \item{nPC}{
a integer indicates the number of principal components used for MPCMR fitting. The default value is the number of principal compoments that just explain more than 95 percent variations.
}
  \item{nL}{
a integer indicates the degree of the polynomial (the number of polynomial basisfunctions)
}
  \item{LMCI}{
logic. Whether to calculate the CI with LM statistic where the basisfunction is the eigenfunction.
}
  \item{LMCI2}{
logic. Whether to calculate the CI with LM statistic where the basisfunction is the polynomial.
}
  \item{nLM}{
a integer indicates the number of increasing point for each dimension when calculating the CI with LM statistic.
}
  \item{Parallel}{
logic. Whether to use parallel computing. The default value is \code{TRUE} and the cores used are the maximal cores munus one.
}
  \item{XYmodel}{
a character indicates the XY model. It should be only used for simulation purpose.
}

}


\details{
Note that you should have individual-level data containing the genetic variants (i.e. genotype) information and the longitudinal information of the exposure of interest. Both information should be contained simultaneously for each individual.

The longitudinal information must contain both the exposure level and its corresponding measured time point (age). It allows for the exposure to be measured at different time points (ages) for every individual, and each individual can have a sparse measurement.

If your individual outcome data is in two-sample with the exposure data or you just wish to treat your data as the two-sample case (e.g. your overlapping sample contains only a small fraction of identical individuals), then obtain the summary statistics from the individual outcome data, and then fit the MPCMR with summary outcome data.
}

\value{
\code{MPCMR_GMM_twosample} retuns a list, consisting of various results, inclduing fitted parameters (and their standard errors), the weak IV assessment results, the IV validity assessment results, the fitted curve.

\item{nPC_used}{how many principal components were used in MPCMR.}
\item{L}{the number of polynomial basis function used for fitting MPCMR.}
\item{K}{the number of eigenfunction basis function used for fitting MPCMR.}
\item{ISres}{the table results of the instrument strength. The columns are coefficient of determination, F value, conditional F, Q statistic value, degree-of-freedom of Q, and the p-value, respectively.}
\item{scatterp}{the MR scatter plot corresponding to the genetic association with the first and second principal compoments.}
\item{one_sample_size}{the sample size of the data finally used for MPCMR fitting.}
\item{IV_validity_test}{the IV validity test results, where the three values are the Q statistic, the degree-of-freedom and the p-value, respectively.}
\item{MPCMRest}{the fitted parameters for the eigenfunction basis set.}
\item{MPCMRvar}{the corresponding variance matrix of \code{MPCMRest}.}
\item{p1}{the fitted curve with eigenfunctin as the bsis function.}
\item{ggdata1}{the dataframe used for producing \code{p1} via ggplot.}
\item{p2}{the fitted curve with polynomial as the bsis function.}
\item{ggdata2}{the dataframe used for producing \code{p2} via ggplot.}
\item{IV_validity_and_basisfunction_test}{the IV validity test results considering the parameteric (polynomial) basis function,  where the three values are the Q statistic, the degree-of-freedom and the p-value, respectively.}

\code{SE},\code{MSE},\code{Co},\code{Coverage_rate},\code{Co_LM},\code{Coverage_rate_LM},\code{sig_points},\code{sig_points_LM} gives some useful information when the true effect function is known (given by the argument \code{XYmodel}). They are used for simulation design.

The result name with the complementary symbol \code{_p} represents the results when the basis function are polynomial functions. For example, \code{MPCMRest_p} is the the fitted parameters for the polynmomial functions.

}

\author{
Haodong Tian
}

\examples{
###see README.md file and TVMR/sim_real_illustration/MPCMR_illustration.R from GitHub
}
