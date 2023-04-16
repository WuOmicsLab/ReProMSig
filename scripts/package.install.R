#

packages <- c("cgwtools","corrplot","DBI","digest","doMC","dplyr","DT",
              "forestplot","glmnet","glue","Hmisc","htmltools","impute",
              "magrittr","mailR","matrixStats","mfp","openxlsx","pec",
              "plyr","purrr","rlang","rmarkdown","rms","RSQLite","superpc",
              "survival","survivalROC","survminer","sva","tibble","yaml")

newPackages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(newPackages) > 0) install.packages(newPackages)

head(installed.packages())


