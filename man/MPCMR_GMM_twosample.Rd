\name{MPCMR_GMM_twosample}
\alias{MPCMR_GMM_twosample}
%- Also NEED an '\alias' for EACH other topic documented here.
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
                    nPC = NA,
                    nL = NA,
                    LMCI = TRUE,
                    LMCI2 = TRUE,
                    nLM = 20,
                    Parallel = TRUE,
                    XYmodel = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Gmatrix}{
a matrix indicates the instrument information
}
  \item{res}{
a list result contains the principal component matrix and the eigenfunctions. \code{res} is usually derived by \code{FPCA}.
}
  \item{by_used}{
a vector indicates the estimated genetic association with the outcome
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
\code{MPCMR_GMM_twosample} retuns a list, consisting of the fitted parameters (and their standard errors), the weak IV assessment results and the fitted curve.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
set.seed(seed)
RES<-getX(N=N_used,J=30,ZXmodel=ZXmodel_used)
res <- FPCA(RES$Ly_sim, RES$Lt_sim,list(dataType='Sparse', error=TRUE, verbose=TRUE))

DAT<-getY(RES, XYmodel=XYmodel_used)

#get another samples for GY data
set.seed(seed+10000)
RES2<-getX(N=N_used,J=30,ZXmodel=ZXmodel_used, MGX_used=RES$details$MGX) #使用之前同一个original data中的genetic effect matrix
DAT2<-getY(RES2, XYmodel=XYmodel_used)


#and prepare the summary information of the GY data
by<-rep(NA,RES2$details$J); sy<-rep(NA,RES2$details$J)
for(j in 1:RES2$details$J){
   GYfit<-lm(  DAT2$Y ~  DAT2[,j] )
   by[j]<-as.numeric(  summary( GYfit )$coef[2,1]  )
   sy[j]<-as.numeric(  summary( GYfit )$coef[2,2]  )
  }
ny<-nrow(DAT2)

#MPCMR fitting
MPCMRres<-MPCMR_GMM_twosample(  Gmatrix=RES$DAT[,1:RES$details$J],
                                 res=res,
                                 by_used=by,
                                 sy_used=sy,
                                 ny_used=ny,
                                 nPC=NA,
                                 nL=NA,
                                 LMCI=TRUE,
                                 LMCI2=TRUE,
                                 nLM=20,
                                 Parallel=TRUE,
                                 XYmodel=XYmodel_used
                                )

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
