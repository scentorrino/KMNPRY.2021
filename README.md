# KMNPRY.2021

This is the R package `KMNPRY.2021' written and maintained by Samuele Centorrino (scentorrino@proton.me), which aims at reproducing the results of the paper by Kahn, Mohaddes, Ng, Pesaran, Raissi, and Yang, 2021, Energy Economics (https://www.sciencedirect.com/science/article/pii/S0140988321004898)

## Installation

You can install the package from this github repository

```r
install.packages('devtools')
devtools::install_github('scentorrino/KMNPRY.2021',build_vignettes=FALSE)
```

The vignette is computationally intensive at this stage. If you wish to run the vignette, change
the option build_vignettes to TRUE in the devtools::install_github() call.

## See the vignette and the code 

To check the main vignette coming with the package, you can type 

```r
vignette("kmnpry_2021_replication", package = "KMNPRY.2021")
```

If you wish to extract the R code used to produce the vignette, you can run the following

```r
library("knitr")
vigSrc <- list.files(pattern = "Rmd$",
                    system.file("doc", package = "KMNPRY.2021"),
                    full.names = TRUE)
purl(vigSrc)
```
or directly copy/paste the .R file from Github to your local machine. 