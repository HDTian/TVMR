# TVMR
Time-varying Mendelian randomization with continuous-time modeling and return the time-varying effect estimate (a functional objective)

If you are interested to learn how the exposure affects the outcome over time (e.g. whether the early-age exposure has more impact on the outcome than the later-age exposure), the time-varying effect function will be useful to you. 

The TVMR method relies on principal components to estimate effect function, therefore also called multi-principal-component MR (MPCMR) 

## Load the package
Install the TVMR package from Github with:
```R
devtools::install_github("HDTian/TVMR")
```
```R
library(TVMR)
```

## Prepare the data
Make sure that you have *individual-level* data containing the genetic variants (i.e. genotype) information and the *longitudinal* information of the exposure of interest. Both information should be contained simultaneously for each individual.

The *longitudinal* information must contain both the exposure level and its corresponding measured time point (age). It allows for the exposure to be measured at different time points (ages) for every individual, and each individual can have a sparse measurement.

You will also have the outcome data, which can be *summary-data* or *individual-data*; *one-sample* or *two-sample* as the exposure data. Depending on the specific data setting, you can use different functions embedded by TVMR. 


## Functional dimension reduction
You can use functional principal components (FPCA) to achieve dimension reduction[^1]. Any packages for FPCA can be used, and here we suggest using `FPCA` function in the package *fdapace*:
```R
install.packages('fdapace'); library(fdapace)
```

First prepare the list `List_exposure`, each entry containing the vector of individual measured exposure, and the corresponding list `List_time`, each containing the corresponding measured timepoints for the exposure.

Then run the FPCA
```R
res <- FPCA(List_exposure, List_time,list(dataType='Sparse', error=TRUE, verbose=TRUE))
```
`res` is a list containing the results of FPCA and will be used for MPCMR fitting.


[^1]: You can use any dimensional reduction to express the individual exposure trajectory by the linear additive sum with the summary non-functional components and their corresponding functional variables. If you choose to use your own components (and the corresponding functions), you can omit the FPCA and directly fit it by MPCMR regarding your components as the 'principal components' and the corresponding functions as the 'eigenfunction' (note that you must use parametric basis function like polynomial rather than eigenfunction in the later estimation) 

## MPCMR fitting



## Results 


