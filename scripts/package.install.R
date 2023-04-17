# Rscript scripts/package.install.R
## packages installed by package remotes ------
if (!require("remotes", quietly = TRUE))
    install.packages("remotes")

packages <- c("cgwtools","openxlsx","doMC","DT",
              "yaml","Hmisc","RSQLite","DBI","corrplot",
              "forestplot","glmnet","superpc","mfp",
              "pec","rms","survminer","survival","survivalROC",
              "dplyr","plyr","purrr","rlang",
              "glue","matrixStats","magrittr","digest",
              "htmltools","rmarkdown","mailR")

newPackages <- packages[!(packages %in% installed.packages()[,"Package"])]
repos <- "http://cran.us.r-project.org"
for(x in newPackages) { 
    if(!require(x, quietly = TRUE)) { remotes::install_version(x, upgrade="never", repos = repos) }
}

## packages installed by package BiocManager ------
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
packages <- c("impute", "plotly", "sva")
newPackages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(newPackages) > 0) BiocManager::install(newPackages)

print("Required R packages are installed successfully!")
