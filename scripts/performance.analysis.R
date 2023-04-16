#
# @Copyright: Peking University Cancer Hospital, All Rights Reserved.
# @Author: Lihua Cao
# @Date: 2021-05-19 17:55:54
# @LastEditTime: 2021-12-17 17:11:02
# @LastEditors: Lihua Cao
# @Description: File content
# Rscript /opt/shiny-server/apps/repromsig/scripts/performance.analysis.R /opt/shiny-server/apps/repromsig/ColoGuide_Stage_II_local/output/sig.ini


# HEADER ------------------------------------------------------------------
rm(list=ls())
options(stringsAsFactors = FALSE)
DEBUG_MODE <- FALSE
ARGS_MODE <- TRUE

# Input/Output/Lib/Config/Params ------------------------------------------
# 1) Parameters
if(ARGS_MODE) {
	args<-commandArgs(TRUE)
	if(length(args)!=1) {
		stop("Usage: Rscript performance.analysis.R [user.config.ini.file]\n")
	}
	user.config.ini.file <- args[1]
} else {
	user.config.ini.file <- "/opt/shiny-server/apps/repromsig/ColoGuide_Stage_II_local/output/sig.ini"
}

# 2) Library
library(dplyr)
library(tibble)

# 3) derive files from conf ini
get_value.ccb <- function(config_file, key = NULL) {
    if(!file.exists(config_file)) {
        warning("Config ini file was not found!")
    } else {
        ini_df <- read.table(config_file, fill = TRUE, header = FALSE,sep='\t',stringsAsFactors=FALSE)
        ini_df <- ini_df[grep(pattern = "^\\[", x = ini_df[, 1], invert = T), ]
        ini_df <- do.call(rbind.data.frame,strsplit(ini_df,'='))
        colnames(ini_df) <- c('key','path')
        ini_df$key <- gsub(" ", "", ini_df$key, fixed = TRUE)
        ini_df$path <- gsub(" ", "", ini_df$path, fixed = TRUE)
        if(any(duplicated(ini_df$key))) {
            stop("exist duplicated keys in your ini file, please check")
        }
        rownames(ini_df) <- ini_df$key
        if(is.null(key)) {
            stop("key is not specified.")
        } else {
            return(as.character(ini_df[key,'path']))
        }
    }
}

script.dir <- get_value.ccb(config_file = user.config.ini.file,  key = 'script_dir')[[1]]
model_analysis_rdata.file <- get_value.ccb(config_file = user.config.ini.file,  key = 'model_analysis_rdata')[[1]]
independence_analysis_rdata.file <- get_value.ccb(config_file = user.config.ini.file,  key = 'independence_analysis_rdata')[[1]]

# 4) define/source functions
source(paste0(script.dir,'ccb.helper.R'))

# readin and preprocess -------------------------------------------
load(model_analysis_rdata.file)

if(nrow(coef1) > 0) {
    load(independence_analysis_rdata.file)

    #lsdata(model_analysis_rdata.file)
    mysetwd(paste0(outpath,'/performance/'))
    rdata.filename <- paste0(b.name.sig, '.performance.analysis.RData')

    train.roc <- train_multi_sig
    debug(train.roc)
    ##
    if(b.sigtype.roc == "score") {
        sams <- intersect(risk_score$SampleID, rownames(train.roc))
        train.roc <- train.roc[sams, ]
        id0 <- which(colnames(train.roc) %in% c("RS", "Signature_score", "Predictive_group"))
        train.roc[,id0] <- as.character(train.roc[,id0])
        train.roc[,id0] <- risk_score[risk_score[,1] %in% sams, "RS"]
    } else {
        colnames(train.roc)[which(colnames(train.roc) == "Signature_score")] <- "Signature_group"
    }

    # ROC evaluation of risk score and other CP features -----------------------------------
    train.roc <- train.roc[order(rownames(train.roc)), ]
    roc_variables <- colnames(train.roc)
    roc_variables <- roc_variables[!roc_variables %in% c("SampleID", "time", "status")]
    
    ##
    models <- list()
    predicts <- list()
    z=1
    key_df <- data.frame(SampleID = rownames(train.roc)) 
    for(x in roc_variables) {
        df <- train.roc[, c('time', 'status', x)]
        df <- df[complete.cases(df),]
        cox0 <- survival::coxph(survival::Surv(time, status == 1)~., data = df, x=TRUE, y=TRUE)
        pred0 <- data.frame(SampleID = rownames(df), RS = predict(cox0, type="lp"))
        predicts[[z]] <- merge.data.frame(key_df, pred0, by = "SampleID", all.x = TRUE) %>%
                                          arrange(SampleID) %>% select(RS) %>%
                                          data.frame(check.names = F)

        ## remove cox models with NA coefficient                                    
        if(length(which(is.na(coefficients(cox0)))) == 0) {
            models[[z]] <- cox0
            names(models)[z] <- x
        }
        z=z+1
    }

    ## the combined model
    if(length(roc_variables) > 1) {
        df <- train.roc
        df <- df[complete.cases(df),]
        cox0 <- survival::coxph(survival::Surv(time, status == 1)~., data = df, x=T, y=T)
        pred0 <- data.frame(SampleID = rownames(df), RS = predict(cox0, type="lp"))
        predicts[[z]] <- merge.data.frame(key_df, pred0, by = "SampleID", all.x = T) %>%
                                          arrange(SampleID) %>% select(RS) %>%
                                          data.frame(check.names = F)

        ## remove cox models with NA coefficient                                    
        if(length(which(is.na(coefficients(cox0)))) == 0) {
            models[[z]] <- cox0
            names(models)[z] <- "Combined"
        }
    }
    ##
    LP <- data.frame(do.call(cbind, predicts), check.names = F)
    colnames(LP)[1:length(roc_variables)] <- roc_variables
    if(length(roc_variables) > 1) { colnames(LP)[ncol(LP)] <- "Combined" }
    LP <- data.frame(train.roc[,c('time','status')], LP, check.names = F)
    
    ##
    roc_variables <- colnames(LP)[-(1:2)]
    key_id <- intersect(c("RS", "Signature_score", "Signature_group", "Predictive_group"), roc_variables)
    if(length(roc_variables) > 1) {
        roc_variables <- c(roc_variables[!roc_variables %in% c("Combined", key_id)], key_id, "Combined")
    }
    LP <- LP[,colnames(LP) %in% c("time", "status", roc_variables)]
    models <- models[intersect(names(models), roc_variables)]
    save(train.roc, LP, models, roc_variables, file = rdata.filename)

    ##
    section <- 'Prediction'
    key <- 'performance_analysis_rdata'
    value <- paste0(paste0(outpath, '/performance/'), rdata.filename)
    system(paste("python3", paste0(script.dir,"/","write_ini.py"), user.config.ini.file, section, key, value))
    message_nosignature <- "Successful!"
} else {
    message_nosignature <- "There was no signature detected!"
}
print(message_nosignature)
