
<!-- README.md is generated from README.Rmd. Please edit that file -->
Installing the CCCBiclust Package
---------------------------------

``` r
install.packages("devtools") # If not yet installed on your R Version
devtools::install_github("hadley/devtools") # Only run this if your currently installed 
                                            # devtools version is <= 1.12 (recursive dependencies bug)

devtools::install_github("ewouddt/CCCBiclust")
```

Should the installation of `CCCBiclust` or `devtools::install_github("hadley/devtools")` throw an error, please install the dependencies manually, then try:

``` r
install.packages(c("flexclust","biclust"))
devtools::install_github("ewouddt/CCCBiclust")
```
