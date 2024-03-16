# TVMR

[Estimating time-varying exposure effects through continuous-time
modelling in Mendelian randomization](https://arxiv.org/abs/2403.05336)

If you are interested to learn how the exposure affects the outcome over
time (e.g. whether the early-age exposure has more impact on the outcome
than the later-age exposure), the time-varying effect function will be
useful to you.

The TVMR method relies on principal components to estimate effect
function, therefore also called multi-principal-component MR (MPCMR)

## Load the package

Install the TVMR package from Github with:

``` r
devtools::install_github("HDTian/TVMR")
```

``` r
library(TVMR)
```

## Prepare the data

Make sure that you have *individual-level* data containing the genetic
variants (i.e. genotype) information and the *longitudinal* information
of the exposure of interest. Both information should be contained
simultaneously for each individual. Denote the data by `Dat`.

The *longitudinal* information must contain both the exposure level and
its corresponding measured time point (age)[^1]. It allows for the
exposure to be measured at different time points (ages) for every
individual, and each individual can have a sparse measurement.

[^1]: the longitudinal information for the measured exposure and the
    measured time point may be stored as a string in one column of `Dat`

You will also have the outcome data, which can be *summary-data* or
*individual-data*; *one-sample* or *two-sample* as the exposure data.
Depending on the specific data setting, you can use different functions
of TVMR. See MPCMR fitting section below.

## Data cleaning

This section reminds you of cleaning data in your study. You may remove
correlated genetic variants, rescale the variables, limit your data to
certain subgroups, remove or impute the missing values of your data, or
choose a specific time region for the exposure.

## Functional dimension reduction

You can use functional principal components analysis (FPCA) to achieve
dimension reduction[^2]. Any packages for FPCA can be used, and here we
suggest using `FPCA` function in the package *fdapace*[^3]:

[^2]: You can use any dimensional reduction to express the individual
    exposure trajectory by the linear additive sum with the summary
    non-functional components and their corresponding functional
    variables. If you choose to use your own components (and the
    corresponding functions), you can omit the FPCA and directly fit it
    by MPCMR regarding your components as the 'principal components' and
    the corresponding functions as the 'eigenfunction' (note that you
    must use parametric basis function like polynomial rather than
    eigenfunction in the later estimation)

[^3]: <https://cran.r-project.org/web/packages/fdapace/index.html>

``` r
install.packages('fdapace'); library(fdapace)
```

First prepare the list `List_exposure`, each entry containing the vector
of individual measured exposure, and the corresponding list `List_time`,
each containing the corresponding measured timepoints for the exposure.

Note that `List_exposure` and `List_time` should be based on `Dat` and
keep the order consistent

Then run the FPCA

``` r
my_res <- FPCA(List_exposure, List_time, list(dataType='Sparse', error=TRUE, verbose=TRUE)  )
```

`my_res` is a list containing the results of FPCA and will be used for
MPCMR fitting.

## MPCMR fitting

Recall that you have the individual-level data for the genotype and the
longitudinal exposure, denoted by `Dat`.

### One-sample individual data

If your individual outcome data is in one sample with `Dat`, prepare the
following variables (based on `Dat`): 
* the matrix of genetic variants:
`my_Gmatrix`
* the vector of outcome: `my_Yvector`

Then run MPCMR:

``` r
MPCMRres<-MPCMR_GMM( Gmatrix=my_Gmatrix, res=my_res ,Yvector=my_Yvector)
```

### Overlapping individual data

If your individual outcome data is in overlap-sample with `Dat` and
assume the ID vector for the exposure and outcome data are `ID_X` and
`ID_Y`, prepare the following variables: 
* the matrix of genetic
variants for the exposure: `my_Gmatrix`;
* the matrix of genetic
variants for the outcome: `my_Gymatrix`;
* the vector of outcome:
`my_Yvector`;
* the ID vector indicating the possible overlapping
samples of the exposure and outcome data:
`myIDmatch <- match( ID_X, ID_Y)`

Then run MPCMR:

``` r
MPCMRres<-MPCMR_GMM( Gmatrix=my_Gmatrix, res=my_res, Yvector=my_Yvector, Gymatrix=my_Gymatrix, IDmatch=my_IDmatch)
```

### Two-sample individual data

If your individual outcome data is in two-sample with `Dat` or you just
wish to treat your data as the two-sample case (e.g. your overlapping
sample contains only a small fraction of identical individuals), then
obtain the summary statistics from the individual outcome data, and
refer to the summary data fit scenario below.

### Summary data

If your outcome data is only summary information, prepare the following
variables: 
* the matrix of genetic variants for the exposure:
`my_Gmatrix`;
* the vector of the genetic association with the outcome:
`my_by_used` (note that the order of the genetic variants in
`my_by_used` should be consistent with that order in `my_Gmatrix`);
* the corresponding vector of standard error for `my_by_used`:
`my_sy_used`;
* the sample size of the outcome data from which the
summary information was obtained: `my_ny_used`

Then run MCPMR:

``` r
MPCMRres<-MPCMR_GMM_twosample( Gmatrix=my_Gmatrix, res=my_res, by_used=my_by_used, sy_used=my_sy_used, ny_used=my_ny_used )
```

## Other arguments in MPCMR fitting

We highlight some arguments in MPCMR fitting.

### User-defined basis function

MPCMR provides two choices of basis function: the eigenfunction and the
polynomial set. The default arguments will run both basis function sets,
and the degree-of-freedom of the both basis function sets is equal to
the number of principal components that just explain 95% variation in
FPCA.

You can use your specific degree-of-freedom for basis function in
time-varying analysis. Assume your exposure data supports four principal
components, and you only wish to use the first 2 eigenfunctions as the
basis function; run

``` r
MPCMRres<-MPCMR_GMM( ..., nPC=2, polyfit=FALSE)
```

(`polyfit=FALSE` will deprecate the MPCMR fitting with polynomial basis
function, so you only fit MRCMR with eigenfunction basis set)

If you only wish to use the linear basis function (i.e. the first two
polynomial functions), run

``` r
MPCMRres<-MPCMR_GMM( ..., nL=2, eigenfit=FALSE)
```

(similarly, `polyfit=FALSE` will deprecate the MPCMR fitting with
eigenfunction basis set, so you only fit MRCMR with polynomial basis
function)

### Identification-robust inference

MPCMR also supports the identification-robust inference based on
Kleibergen's Lagrange multiplier (LM) statistic.

The default setting will calculate the LM confidence intervals and will
occupy most of your laptop cores (if you use an 8-core machine, 7 cores
will be used). The calculation of LM confidence intervals could be very
time-consuming and also occupy much RAM. One 8-core machine could take
5\~10 hours to finish the calculation of LM CI. If you are impatient, we
suggest running MPCMR with HPC. Also be careful about the RAM in HPC, as
insufficient memory could easily cause an out-of-memory problem. If you
use HPC, we suggest claiming more total RAM or Max RAM per CPU (e.g. via
SLURM) and also controlling the cores used.

If you do not wish to calculate the LM CI, simply run

``` r
MPCMRres<-MPCMR_GMM( ..., LMCI=FALSE, LMCI2=FALSE)
```

## Results

`MPCMRres` contains many results you may wish to take. Try `?MPCMR_GMM`
to see more details.

To see the instrument strength information:

``` r
MPCMRres$ISres
```

``` r
#            RR         F       cF   Qvalue df pvalue
#PC1 0.02214637  7.525909 6.553622 174.5873 29      0
#PC2 0.02960446 10.137684 8.827976 213.5132 29      0
```

(here columns are coefficient of determination, F value, conditional F,
Q statistic value, degree-of-freedom of Q, and the p-value,
respectively)

To draw the fitted effect function:

``` r
MPCMRres$p1  #MPCMRres$p2 returns the polynomial fitting result
```

<img src="https://github.com/HDTian/TVMR/assets/127906571/2b6b7616-b4e8-4720-ba8c-eecea77fbbca" alt="fit_example" width="765"/>

(here the black solid curve is the fitted effect function, the black
dashed curves are CI via GMM, and the grey dashed curves are LM CI. The
blue curve is the true effect function.)
