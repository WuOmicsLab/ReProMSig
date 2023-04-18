#
# @Copyright: Peking University Cancer Hospital, All Rights Reserved.
# @Author: Lihua Cao
# @Date: 2022-01
# @LastEditTime: 2023-04-16
# @LastEditors: Lihua Cao
# @Description: Run the installation script to install the required packages automatically before using the pipeline.
#

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
if(!require("riskRegression", quietly = TRUE)) { remotes::install_version("", upgrade="never", repos = repos, version = "1.4.3") }
for(x in newpkgs) {
    if(!require(x, quietly = TRUE)) { remotes::install_version(x, upgrade="never", repos = repos) }
}

## packages installed by package BiocManager ---
if(!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if(!require("impute", quietly = TRUE)) { BiocManager::install("impute") }

print("Required R packages are installed successfully!")




