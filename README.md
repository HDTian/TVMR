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
Make sure that you have *individual-level* data containing the genetic variants (i.e. genotype) information and the *longitudinal* information of the exposure of interest. Both information should be contained simultaneously for each individual. Denote the data by `Dat1`.

The *longitudinal* information must contain both the exposure level and its corresponding measured time point (age)[^1]. It allows for the exposure to be measured at different time points (ages) for every individual, and each individual can have a sparse measurement.

[^1]: the longitudinal information for the measured exposure and the measured time point may be stored as a string in one column of `Dat1`

You will also have the outcome data, which can be *summary-data* or *individual-data*; *one-sample* or *two-sample* as the exposure data. Depending on the specific data setting, you can use different functions of TVMR. See MPCMR fitting section below. 

## Data cleaning
This section reminds you of cleaning data in your study. You may remove correlated genetic variants, limit your data to certain subgroups, remove or permute the missing values of your data, or choose a specific time region for the exposure.


## Functional dimension reduction
You can use functional principal components (FPCA) to achieve dimension reduction[^2]. Any packages for FPCA can be used, and here we suggest using `FPCA` function in the package *fdapace*[^3]:
```R
install.packages('fdapace'); library(fdapace)
```

First prepare the list `List_exposure`, each entry containing the vector of individual measured exposure, and the corresponding list `List_time`, each containing the corresponding measured timepoints for the exposure. `List_exposure` and `List_time` should be based on `Dat1`

Then run the FPCA
```R
res <- FPCA(List_exposure, List_time, list(dataType='Sparse', error=TRUE, verbose=TRUE)  )
```
`res` is a list containing the results of FPCA and will be used for MPCMR fitting.


[^2]: You can use any dimensional reduction to express the individual exposure trajectory by the linear additive sum with the summary non-functional components and their corresponding functional variables. If you choose to use your own components (and the corresponding functions), you can omit the FPCA and directly fit it by MPCMR regarding your components as the 'principal components' and the corresponding functions as the 'eigenfunction' (note that you must use parametric basis function like polynomial rather than eigenfunction in the later estimation) 

[^3]: https://cran.r-project.org/web/packages/fdapace/index.html

## MPCMR fitting
Recall that you have the individual-level data for the genotype and the longitudinal exposure, denoted by `Dat`. 

Prepare the following  


## Results 


