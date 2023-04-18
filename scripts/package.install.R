
# Rscript scripts/package.install.R
rm(list=ls())

## packages installed by package remotes ---
if (!require("remotes", quietly = TRUE))
    install.packages("remotes")

pkgs <- c("cgwtools","corrplot",
          "digest","doMC","dplyr","DT",
          "forestplot","glmnet","glue","Hmisc",
          "htmltools","matrixStats","mfp",
          "pec","plyr","rlang","rmarkdown",
          "rms", "superpc","survival",
          "survivalROC","survminer",
          "yaml","tibble","tools")
newpkgs <- pkgs[!(pkgs %in% installed.packages()[, "Package"])]

##
repos <- "http://cran.us.r-project.org"
for(x in newpkgs) {
    if(!require(x, quietly = TRUE)) { remotes::install_version(x, upgrade="never", repos = repos) }
}

## packages installed by package BiocManager ---
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if(!require("impute", quietly = TRUE)) { BiocManager::install("impute") }

print("Required R packages are installed successfully!")



