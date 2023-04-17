
# Rscript scripts/independence.analysis.R ColoGuide_Stage_II_local/output/sig.ini
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
		stop("Usage: Rscript independence.analysis.R [user.config.ini.file]\n")
	}
	user.config.ini.file <- args[1]
} else {
	user.config.ini.file <- "ColoGuide_Stage_II_local/output/sig.ini"
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

# 4) define/source functions
source(paste0(script.dir,'ccb.helper.R'))

# readin and preprocess -------------------------------------------
load(model_analysis_rdata.file)
# cgwtools::lsdata(model_analysis_rdata.file)
mysetwd(paste0(outpath, '/independence/'))
rdata.filename <- paste0(b.name.sig, '.independence.analysis.RData')

# survival data preprocess ---------------------------------------------------------
if(nrow(coef1) > 0) {
    if(b.signature.type == "Predictive") {
        risk_score0 <- train_stratify %>% select(c("SampleID", "Type", "Treatment_type"))
        colnames(risk_score0)[2] <- "Predictive_group"
        if(only_treatment_pts_train) {
            train_surv0 <- rbind(train_surv, train_surv_sg)
            b.upload.ann0 <- rbind(b.upload.ann, b.upload.ann.sg)
        } else {
            train_surv0 <- train_surv
            b.upload.ann0 <- b.upload.ann
        }
        risk_score0$Predictive_group <- relevel(as.factor(risk_score0$Predictive_group), ref = "No-benefit")
        risk_score0$Treatment_type <- relevel(as.factor(risk_score0$Treatment_type), ref = "Control")
    } else {
        risk_score0 <- train_stratify %>% select(c("SampleID", "Type"))
        colnames(risk_score0)[2] <- "RS"
        risk_score0$RS <- relevel(as.factor(risk_score0$RS), ref = "Low risk")
        train_surv0 <- train_surv
        b.upload.ann0 <- b.upload.ann
    }
    str(risk_score0)
    train <- merge.data.frame(x = risk_score0, y = train_surv0, by='SampleID', all.x=T)

    ##
    if(!is.null(b.clinical.variables) | !is.null(b.interaction.test)) {
        if(b.signature.type == "Predictive") { b.clinical.variables = b.interaction.test }
        variables <- intersect(colnames(b.upload.ann0), b.clinical.variables)
        variables <- variables[!variables %in% c("time", "status")]
        train <- merge.data.frame(x = train, y = b.upload.ann0[, c("SampleID", variables)], by='SampleID', all.x=T)
    }
    rownames(train) <- train$SampleID

    ## preprocess train, leves with < 3 values were setted to NA
    ids <- colnames(train)
    ids <- ids[!ids %in% c("SampleID", "time", "status")]
    for(x in ids) {
        y <- train[,x]
        if(!is.numeric(y)) {
             y <- table(y)
             y <- names(y[y < 3])
             if(length(y) > 0) {
                 train[train[,x] %in% y, x] <- NA
             }
        }
    }

    ##
    variables0 <- colnames(train)
    variables0 <- unlist(lapply(variables0, function(x) {if(n_distinct(train[!is.na(train[,x]),x]) > 1) x }))
    train <- train[,variables0]
    debug(train)
    str(train)

    # Standard cox regression -------------------------------------------------
    variables <- colnames(train)
    variables <- variables[!(variables %in% c("SampleID", "time", "status"))]
    variables <- intersect(unique(c("Predictive_group", "RS", sort(variables))), variables)

    if(!is.null(b.interaction.test)) {
        variables <- intersect(variables, c("Predictive_group", "Treatment_type", b.interaction.test))
    } else if(is.null(b.interaction.test) & b.signature.type == "Predictive") {
        variables <- c("Predictive_group", "Treatment_type")
    }
 
    ## univariate and multivariate cox ------------------------------------------------
    clinic_ann <- train[ ,colnames(train) %in% c("SampleID", variables)]
    pts_surv <- train[, c('SampleID','time','status')]
    
    ##
    uni_sig_multi = FALSE
    if(b.signature.type == "Predictive") {
        uni_cox_out <- uni_cox.preditive(pts_surv, clinic_ann, treatment_columnID = "Treatment_type")
        multi_cox_out <- multi_cox.preditive(pts_surv, clinic_ann, treatment_columnID = "Treatment_type")
    } else {
        uni_cox_out <- uni_cox.prognosis(pts_surv, clinic_ann)
        if(uni_sig_multi) {
            df <- uni_cox_out %>% filter(!is.na(P) & P != "")
            mode(df$P) <- "numeric"
            variables0 <- unique(df[df$P < 0.1, "Type"])
            clinic_ann <- clinic_ann[, c("SampleID", variables0)]
        }
        multi_cox_out <- multi_cox.prognosis(pts_surv, clinic_ann)
    }

    ## combine uni and multivariate into one table ---
    sig_variable = c()
    uni_multi_cox_report = data.frame()

    if(is.null(b.clinical.variables)) { multi_cox_out <- uni_cox_out }

    res <- merge.data.frame(uni_cox_out, multi_cox_out, by = c("Type", "Variable"),
                            all.x = T, sort = F,
                            suffixes = c("<br />Univariate", "<br />Multivariate"))

    web_variables <- c("Type", "Variable", "Events/N<br />Univariate",
                       "HR for events (95% CI)<br />Univariate",
                       "P<br />Univariate", "P_interaction<br />Univariate",
                       "Events/N<br />Multivariate",
                       "HR for events (95% CI)<br />Multivariate",
                       "P<br />Multivariate",
                       "P_interaction<br />Multivariate")
    res <- res[, colnames(res) %in% web_variables]

    if(b.signature.type == "Predictive") {
        colnames(res) <- c("Variable", "Value", "No. of Events<br />Univariate analysis",
                           "Hazard Ratio (95% CI)<br />Univariate analysis",
                           "P<br />Univariate analysis", "P interaction<br />Univariate analysis",
                           "No. of Events<br />Multivariate analysis",
                           "Hazard Ratio (95% CI)<br />Multivariate analysis",
                           "P<br />Multivariate analysis",
                           "P interaction<br />Multivariate analysis")
    } else {
        colnames(res) <- c("Variable", "Value", "No. of Events<br />Univariate analysis",
                           "Hazard Ratio (95% CI)<br />Univariate analysis",
                           "P<br />Univariate analysis", "No. of Events<br />Multivariate analysis",
                           "Hazard Ratio (95% CI)<br />Multivariate analysis",
                           "P<br />Multivariate analysis")
    }
    
    res[which(res$"No. of Events<br />Univariate analysis" != ""),
              c("No. of Events", "No. of Patients")] <- do.call(rbind.data.frame,
                                                            strsplit(res$"No. of Events<br />Univariate analysis", "/")
                                                        )

    ##
    if(b.signature.type == "Predictive") {
        res <- res %>% select("Variable", "Value", "No. of Patients", "No. of Events",
                              "Hazard Ratio (95% CI)<br />Univariate analysis",
                              "P<br />Univariate analysis", "P interaction<br />Univariate analysis",
                              "No. of Events<br />Multivariate analysis",
                              "Hazard Ratio (95% CI)<br />Multivariate analysis",
                              "P<br />Multivariate analysis",
                              "P interaction<br />Multivariate analysis")
    } else {
        res <- res %>% select("Variable", "Value", "No. of Patients", "No. of Events",
                            "Hazard Ratio (95% CI)<br />Univariate analysis",
                            "P<br />Univariate analysis", "No. of Events<br />Multivariate analysis",
                            "Hazard Ratio (95% CI)<br />Multivariate analysis",
                            "P<br />Multivariate analysis")
    }
    res[is.na(res)] <- ""
    uni_multi_cox_report <- res

    ## multivariate sig variables  ------------------------------------------------
    if(is.null(b.clinical.variables)) { multi_cox_out = data.frame() }
    mode(multi_cox_out$P) <- 'numeric'
    sig_variable <- sort(unique(multi_cox_out[!is.na(multi_cox_out$P) &
                                multi_cox_out$P < 0.05 &
                                multi_cox_out$P != "", "Type"]))
    sig_variable <- unique(c(sig_variable, "RS", "Predictive_group"))

    ##
    train_multi_sig <- train[, intersect(c('time', 'status', sig_variable), colnames(train))]

    if(!is.null(b.nomogram.variables) & b.nomogram.option == "user") {
        nom_vars <- intersect(colnames(b.upload.ann0), b.nomogram.variables)
        nom_vars0 <- nom_vars[!nom_vars %in% colnames(train)]
        if(length(nom_vars0) > 0) {
            train_nom <- merge.data.frame(x = train, y = b.upload.ann0[, c("SampleID", nom_vars0)], by='SampleID', all.x=T)
        } else {
            train_nom = train
        }
        sig_variable <- unique(c(nom_vars, "RS", "Predictive_group"))
        train_multi_sig <- train_nom[, intersect(c('time', 'status', sig_variable), colnames(train_nom))]
        rownames(train_multi_sig) <- train_nom$SampleID
    } else if(is.null(sig_variable) | (is.null(b.nomogram.variables) & b.nomogram.option == "user" & b.have.profiles == "Yes")) {
        sig_variable <- c("RS", "Predictive_group")
        train_multi_sig <- train[, intersect(c('time', 'status', sig_variable), colnames(train))]
    }

    ##
    uni_cox <- uni_cox_out[,-1]
    uni_cox <- uni_cox[!colnames(uni_cox) %in% "Coefficient"]

    multi_cox <- multi_cox_out[,-1]
    multi_cox <- multi_cox[!colnames(multi_cox) %in% "Coefficient"]

    ##
    ifelse(b.signature.type == "Predictive",
           labeltext <- c("Value", "Events/N\n", "HR", "lower95", "upper95", "Hazard Ratio (95% CI)\n", "P\n", "P_interaction\n"), 
           labeltext <- c("Value", "Events/N", "HR", "lower95", "upper95", "Hazard Ratio (95% CI)", "P"))
    colnames(uni_cox) <- labeltext
    if(!is.null(b.clinical.variables)) {
        colnames(multi_cox) <- labeltext
    }

    ##
    colnames(train_multi_sig)[which(colnames(train_multi_sig) == "RS")] <- "Signature_score"
    uni_cox[uni_cox=="RS"] <- "Signature score"
    multi_cox[multi_cox=="RS"] <- "Signature score"
    uni_multi_cox_report[uni_multi_cox_report=="RS"] <- "Signature score"

    ##
    if(b.signature.type == "Predictive") {
        colnames(uni_multi_cox_report)[2] <- "Value<br />(Control, Treatment)"
        colnames(uni_cox)[1] <- "Value\n(Control, Treatment)"
        if(!is.null(b.clinical.variables)) { colnames(multi_cox)[1] <- "Value\n(Control, Treatment)" }
    }

    uni_multi_cox_report <- uni_multi_cox_report %>% select(-c("No. of Events<br />Multivariate analysis"))
    uni_multi_cox_report <- data.frame("Raw Rank" = 1:nrow(uni_multi_cox_report), uni_multi_cox_report, check.names = F)
    
    ##
    train_multi_sig_nom_t = data.frame()
    train_multi_sig_nom_c = data.frame()

    if(b.have.profiles == "No") {
        nom_ids <- sort(unique(coef1$Type))
        nom_ids <- do.call(rbind, lapply(nom_ids, function(x) {
            if((!class(Train[,x]) %in% c("character", "factor")) & b.predictor.selection == "Yes") {
                y <- coef1[coef1$Type == x, "Variable"]
                c(paste0(x, "_mfp"), paste0(x, "_mfp_", y))
            } else { c(x, x) }
        }))
        sig_variables <- c("time", "status", nom_ids[,2])

        ## all patients
        train_multi_sig_nom <- train_model[, sig_variables]
        colnames(train_multi_sig_nom) <- c("time", "status", nom_ids[,1])

        ## nomogram data for control and treatment patients for predictive signature
        if(b.signature.type == "Predictive") {
            chemo_sams0 <- intersect(chemo_sams, rownames(traintopred))
            sg_sams0 <- intersect(sg_sams, rownames(traintopred))
            train_multi_sig_nom_t <- traintopred[chemo_sams0, sig_variables]
            colnames(train_multi_sig_nom_t) <- c("time", "status", nom_ids[,1])

            train_multi_sig_nom_c <- traintopred[sg_sams0, sig_variables]
            colnames(train_multi_sig_nom_c) <- c("time", "status", nom_ids[,1])
        }
    } else {
        train_multi_sig_nom <- train_multi_sig
        ## nomogram data for control and treatment patients for predictive signature
        if(b.signature.type == "Predictive") {
            chemo_sams0 <- intersect(chemo_sams, rownames(train_multi_sig))
            sg_sams0 <- intersect(sg_sams, rownames(train_multi_sig))
            train_multi_sig_nom_t <- train_multi_sig[chemo_sams0,]
            train_multi_sig_nom_c <- train_multi_sig[sg_sams0, ]
        }
    }

    ## 
    if(is.null(b.clinical.variables)) {
        uni_multi_cox_report <- uni_multi_cox_report[, -grep("Multivariate analysis", colnames(uni_multi_cox_report))]
    }
    save(train_multi_sig, train_multi_sig_nom,
         train_multi_sig_nom_t,
         train_multi_sig_nom_c,
         uni_cox, multi_cox,
         uni_multi_cox_report,
         file = rdata.filename)

    section <- 'Building_signature'
    key <- 'independence_analysis_rdata'
    value <- paste0(paste0(outpath,'/independence/'), rdata.filename)
    system(paste("python3", paste0(script.dir,"/","write_ini.py"), user.config.ini.file, section, key, value))
    message_nosignature <- "Successful!"
} else {
    message_nosignature <- "There was no signature detected!"
}
print(message_nosignature)

