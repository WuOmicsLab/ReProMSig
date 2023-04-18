#
# @Copyright: Peking University Cancer Hospital, All Rights Reserved.
# @Author: Lihua Cao
# @Date: 2022-01
# @LastEditTime: 2023-04-16
# @LastEditors: Lihua Cao
# @Description: Export the RData file that could be uploaded to the ReProSig website, for displaying and sharing.
#

# HEADER ------------------------------------------------------------------
rm(list=ls())
options(stringsAsFactors = FALSE)
DEBUG_MODE = FALSE
ARGS_MODE = TRUE

# Input/Output/Lib/Config/Params ------------------------------------------
# 1) Parameters
if(ARGS_MODE) {
  args<-commandArgs(TRUE)
  if(length(args)!=2) {
    stop("Usage: Rscript local.rdata.generate.R [user.config.ini.file] [tripod.ini.file]\n")
  }
  user.config.ini.file <- args[1]
  tripod.ini.file <- args[2]
} else {
  user.config.ini.file <- "ColoGuide_Stage_II_local/output/sig.ini"
  tripod.ini.file <- "ColoGuide_Stage_II_local/output/tripod.ini"
}

# 2) Library
library(dplyr)

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
user_filtered_datasets.rdata <- get_value.ccb(config_file = user.config.ini.file,  key = 'user_filtered_datasets_rdata')[[1]]
user_parameters.rdata <- get_value.ccb(config_file = user.config.ini.file,  key = 'user_parameters_rdata')[[1]]

model_analysis_rdata <- get_value.ccb(config_file = user.config.ini.file, key = 'model_analysis_rdata')[[1]]
independence_analysis_rdata <- get_value.ccb(config_file = user.config.ini.file, key = 'independence_analysis_rdata')[[1]]
performance_analysis_rdata <- get_value.ccb(config_file = user.config.ini.file, key = 'performance_analysis_rdata')[[1]]
external_evaluate_analysis_rdata <- get_value.ccb(config_file = user.config.ini.file, key = 'external_evaluate_analysis_rdata')[[1]]
add_info_tripod_rdata <- get_value.ccb(config_file = user.config.ini.file, key = 'add_info_tripod_rdata')[[1]]

table_cp_stat_list_rdata <- get_value.ccb(config_file = tripod.ini.file, key = 'table_cp_stat_list_rdata')[[1]]
table_surv_prob_rdata <- get_value.ccb(config_file = tripod.ini.file, key = 'table_surv_prob_rdata')[[1]]
user_tripod_rdata <- get_value.ccb(config_file = tripod.ini.file, key = 'user_tripod_rdata')[[1]]

# 4) define/source functions
source(paste0(script.dir,'ccb.helper.R'))


# readin and preprocess  -------------------------------------------------------------------------
load(user_filtered_datasets.rdata)
load(user_parameters.rdata)

load(model_analysis_rdata)
load(independence_analysis_rdata)
load(performance_analysis_rdata)
load(external_evaluate_analysis_rdata)
load(add_info_tripod_rdata)

load(table_cp_stat_list_rdata)
load(table_surv_prob_rdata)
load(user_tripod_rdata)

## parameters ------------
outpath <- b.user.sig.path
mysetwd(outpath)

dir.create(paste0(b.user.sig.path, "/upload/"))
rdata.filename <- paste0(outpath, "/upload/", b.name.sig.id, ".", Sys.Date(), ".local.RData")

signature_result_introduction_text <- paste0("<ul>
						<p><li><b>Summary</b> tab shows summary of the developed signature, clinicopathological statistics, molecular profiles distribution for used training dataset and validation dataset(s)</li></p>
                                             <p><li><b>Signature details</b> tab introduces the performance of the developed signature from multiple perspectives tab shows the details of the developed prediction model, including predictors, regression coefficients, hazard ratios (HRs) accompanied by 95% confidence interval, and <I>P</I> values. For expression signature, the Pearson correlation between each pair of the signature predictors is plotted for visual check if their expression/methylation levels are independent.</li></p>
                                             <p><li><b>Signature performance</b> tab introduces the performance of the developed signature from multiple perspectives</li></p>
                                             <ul>			
                                             <p><li>For prognostic molecular signature, independence test is used to evaluate the independence of the developed molecular signature, compared with and the clinicopathological variables in the training dataset.</li></p>
                                             <p><li>For predictive molecular signature, univariate and multivariate Cox regression analyses were performed to test for interaction between treatment response and the predicted treatment-benefit groups in terms of survival, in the training dataset.</li></p>
                                             <p><li>Kaplan-Meier method is applied to assess the difference of survival between different patient groups.</li></p>
                                             <p><li>The ",paste0(b.times,collapse=", "),"-months time-dependent receiver operating characteristic (ROC) is generated to quantify and compare the discrimination ability of the molecular signature and clinicopathological predictors in the training dataset.</li></p>
                                             <p><li>The calibration plots at ",paste0(b.times,collapse=", "),"-months are generated to assess the calibration performance of the prognostic/predictive molecular signatures by comparing predicted probabilities and observed outcome.</li></p>
                                             <p><li>The prediction error (PE) curve is generated to evaluate overall performance of the developed model in the training dataset.</li></p>
                                             </ul>
                                             <p><li><b>Clinicopathological association</b> tab: Association analysis is applied to evaluate the correlation between the developed prognosis model and clinicopathological variables in the training and validation datasets respectively.</li></p>
                                             <p><li><b>Nomogram</b> tab shows the details of tools helping clinical decision, based on the developed model.</li></p>
                                             <ul>
                                             <p><li>Nomogram is be generated for predicting ",paste0(b.times,collapse=", "),"-months ",b.endpoint," probabilities of a new case, which may serve as a tool facilitating clinical decision.</li></p>
                                             <p><li>The calibration plots at ",paste0(b.times,collapse=", "),"-months are generated to assess the calibration performance of the Nomogram by comparing predicted probabilities and observed outcome.</li></p>
                                             </ul>
                                             <p><li><b>TRIPOD report</b> tab shows a single structured HTML file following TRIPOD guideline check items, which can be edited, downloaded and shared to public.</li></p>
                                             </ul>")

## parameters selected in step1 of web version --
local_msig <- TRUE
b.min.age <- NA
b.T <- c("T1", "T2", "T3", "T4", "NA")
b.N <- c("N0", "N1", "N2", "N3", "NA")
b.M <- c("M0", "M1", "NA")
b.Gender <- c("Male", "Female", "NA")
b.ethnicity <- c("Asian", "Black or African American", "White", "NA")
b.stage.edition <- c("AJCC 6th", "AJCC 7th", "AJCC 8th", "NA")
b.stage <- c("I", "II", "III", "IV", "NA")
b.is.primary.treatment <- "Yes"

## 
b.predictive.treatment <- unique(c(b.control.treatment.type, b.treatment.treatment.type))
b.predictive.treatment.setting <- unique(c(b.control.treatment.setting,b.treatment.treatment.setting))
b.predictive.regimen <- unique(c(b.control.regimen,b.treatment.regimen))

##
Train <- as.matrix(Train)
debug(Train)
save(local_msig,
      # user-uploaded datas
      user_filtered_datasets_info, 
      b.upload.predictors,
      b.clinicopathological.predictors,

      # step1 parameters
      b.signature.type, b.is.primary.treatment,
      b.disease.type, b.sample.type, b.primary.site,
      b.have.profiles, b.profiling, b.endpoint, b.min.age,
      b.stage.edition, b.stage, b.T, b.N, b.M, b.Gender, b.ethnicity,
      b.prognostic.treatment, b.prognostic.treatment.setting,
      b.cp.variables, predictor_type, 

      # for predictive signature
      b.predictive.treatment, b.predictive.treatment.setting, b.predictive.regimen,
      b.control.treatment.type, b.control.treatment.setting, b.control.regimen,
      b.treatment.treatment.type, b.treatment.treatment.setting, b.treatment.regimen,
      
      # step3 parameters
      b.user.sig.path, b.name.sig, b.name.sig.id,  
      b.times, b.group, b.group.cutoff, b.cutpoint.method,
      b.high2, b.low2, b.high3, b.low3, b.moderate3, 
      b.clinical.variables, b.interaction.test,
      b.combine.variables, b.combine.age, b.combine.stage1, b.combine.stage2,
      b.combine.n1, b.combine.n2, b.combine.t1, b.combine.t2,  
      b.predictor.selection, b.predictor.selection.methods,
      b.bootstrap.frequency, b.bootstrap.iterations,

      # step2 parameters
      select_training_dataset, used_validation_dataset,
      
      # parameters generated by signature building and evaluation analysis
      patient_characteristic_legend, patient_characteristic_summary.df,
      train_valid_df_list_names, 
      signature_result_introduction_text, signature_summary_text,
      signature_datasets_text, coef1, coef_report, risk_score_formula, 
      b.profile_data.list, Train, train_stratify, roc_variables,

      uni_multi_cox_report, pts_ann_sb, testsets.name, test_pts_ann_sb.list,
      

      train_surv, train_surv_sg, test_surv.list,
      continuous_annvar, geneset_name,
      valid_warning.list,

      ## report paras
      generator, email, organization, signature_description, 
      local_sig_df, association_text, 
      local_molecular_profiling, local_min_age, local_max_age,
      local_model_summary, local_model_text,

      file = rdata.filename
)

##
if(nrow(coef1) > 0) {
    cgwtools::resave(
      risk_score, train_multi_sig, test_rs_plot.list, 
      CP_Stat_report, surv_prob_report, 

      uni_cox, multi_cox, LP, models, train.roc,
      train_multi_sig_nom, train_multi_sig_nom.list,
      train_surv, train_surv_sg, km_input.list,
      figure_km_validsets_exist, valid_KM_prex.list,
      file = rdata.filename
  )
}