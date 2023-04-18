# Rscript scripts/tripod.report.html.R ColoGuide_Stage_II_local/input/reporting.yaml ColoGuide_Stage_II_local/output/sig.ini ColoGuide_Stage_II_local/output/tripod.ini


# HEADER ------------------------------------------------------------------
rm(list=ls())
options(stringsAsFactors = FALSE)
DEBUG_MODE = FALSE
ARGS_MODE = TRUE

# Input/Output/Lib/Config/Params ------------------------------------------
# 1) Parameters
if(ARGS_MODE) {
  args<-commandArgs(TRUE)
  if(length(args)!=3) {
    stop("Usage: Rscript tripod.report.html.R [reporting.yaml.file] [user.config.ini.file] [tripod.ini.file]\n")
  }
  reporting.yaml.file <- args[1]
  user.config.ini.file <- args[2]
  tripod.ini.file <- args[3]
} else {
  reporting.yaml.file <- "ColoGuide_Stage_II_local/input/reporting.yaml"
  user.config.ini.file <- "ColoGuide_Stage_II_local/output/sig.ini"
  tripod.ini.file <- "ColoGuide_Stage_II_local/output/tripod.ini"
}


# 2) Library
library(yaml)
library(cgwtools)
library(DT)
library(dplyr)
library(digest)
library(htmltools)

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

getValidationValues <- function(yaml_key1,yaml_key2) {
  value <- unique(c(conf[[yaml_key1]][[yaml_key2]]$`validation_cohort-1`,conf[[yaml_key1]][[yaml_key2]]$`validation_cohort-2`,conf[[yaml_key1]][[yaml_key2]]$`validation_cohort-3`))
  if(!is.null(value)){
    value <- gsub("\n"," <br />", value) 
    value <- paste0(value,collapse = "<br /><br /> ")
  } 
  return(value)
}

script.dir <- get_value.ccb(config_file = user.config.ini.file, key = 'script_dir')[[1]]
repromsig.dir <- get_value.ccb(config_file = user.config.ini.file, key = 'repromsig_dir')[[1]]
user_parameters.rdata <- get_value.ccb(config_file = user.config.ini.file, key = 'user_parameters_rdata')[[1]]
user_filtered_datasets.rdata <- get_value.ccb(config_file = user.config.ini.file, key = 'user_filtered_datasets_rdata')[[1]]
user_uploaded_anno.rdata <- get_value.ccb(config_file = user.config.ini.file, key = 'user_uploaded_anno_rdata')[[1]]
b.user.sig.path <- get_value.ccb(config_file = user.config.ini.file, key = 'b_user_sig_path')[[1]]

# readin and preprocess ---------------------------
conf <- read_yaml(reporting.yaml.file)

# source functions ---
source(paste0(script.dir, "/ccb.helper.R"))
js_dir = paste0(repromsig.dir, "/config/")


# user input ---------------------
generator <- conf$Lead_contact
email <- conf$Email_address
organization <- conf$Organization
signature_description <- conf$Signature_summary$Signature_description

item4a_training_source <- conf$Item4_Source_of_data$Item4a_training_Study_design
item4a_validation_source <- getValidationValues(yaml_key1 = "Item4_Source_of_data",yaml_key2 = "Item4a_validation_Study_design")
item4a_training_study_date <- conf$Item4_Source_of_data$Item4b_training_Key_study_dates
item4a_validation_study_date <- getValidationValues(yaml_key1 = "Item4_Source_of_data",yaml_key2 = "Item4b_validation_Key_study_dates")

item5a_training_setting <- conf$Item5_Participants$Item5a_training_Key_study_setting
item5a_validation_setting <- getValidationValues(yaml_key1 = "Item5_Participants",yaml_key2 = "Item5a_validation_Key_study_setting")
item5b_training_exclusion_user <- conf$Item5_Participants$Item5b_training_Exclusion_criteria_Clinical_anlaysis
item5b_training_inclusion_user <- conf$Item5_Participants$Item5b_training_Inclusion_criteria_Clinical_anlaysis
item5b_validation_exclusion_user <- getValidationValues(yaml_key1 = "Item5_Participants",yaml_key2 = "Item5b_validation_Exclusion_criteria_Clinical_anlaysis")
item5b_validation_inclusion_user <- getValidationValues(yaml_key1 = "Item5_Participants",yaml_key2 = "Item5b_validation_Inclusion_criteria_Clinical_anlaysis")
item5b_validation_similar <- conf$Item5_Participants$Item5b_Is_it_similar_to_the_eligibility_criteria_used_in_training_dataset
item5c_training <- conf$Item5_Participants$Item5c_training_Received_treatments
item5c_validation <- getValidationValues(yaml_key1 = "Item5_Participants",yaml_key2 = "Item5c_validation_Received_treatments")

item6a_training_outcome <- conf$Item6_Outcome$Item6a_outcome
item6a_training_methods <- conf$Item6_Outcome$Item6a_How_and_when_assessed
item6b_training <- conf$Item6_Outcome$Item6b_Actions_for_blind_assessment_of_the_outcome

item7a_training <- conf$Item7_Predictors$Item7a_Predictors_used_in_developing_the_multivariable_prediction_model
item7b_training <- conf$Item7_Predictors$Item7b_Actions_for_blind_assessment_of_predictors

item8_training <- conf$Item8_Sample_Size$Item8_How_the_study_size_arrived_at

item9_training_methods_user <- conf$Item9_Missing_Data$Item9_training_Missing_value_handling
item9_validation_methods_user <- conf$Item9_Missing_Data$Item9_validation_Missing_value_handling
item10e_updating <- conf$Item10_Statistical_Analysis_Methods$Item10e_Model_updating_arising_from_the_validation

item12_setting <- conf$Item12_Development_vs_Validation$Item12_Setting
item12_eligibility <- conf$Item12_Development_vs_Validation$Item12_Eligibility_criteria
item12_outcome <- conf$Item12_Development_vs_Validation$Item12_Outcome
item12_predictors <- conf$Item12_Development_vs_Validation$Item12_Predictors

user_upload_diagram_path <- conf$Item13_Participants$Item13a_Study_flow_diagram

# web ---------------------
item10a_handling_molecular = item10a_implausible = item10a_handling_clinical = item10b_before_genes = item10b_before = item10b_model = item10c_signature = item10c_nomogram = item10c_online = item10d_independence = item10d_interaction = item10d_roc = item10d_km = item10d_pe = item10d_calibration = item10d_km_pro = item10d_association = item11_groups = ""


# derive files from conf ini ------------
user.parameters.rdata <- user_parameters.rdata
user.filtered.datasets.rdata <- user_filtered_datasets.rdata

penalized_cox_analysis_rdata.file <- get_value.ccb(config_file = user.config.ini.file, key = 'model_analysis_rdata')[[1]]
independence_analysis_rdata.file <- get_value.ccb(config_file = user.config.ini.file, key = 'independence_analysis_rdata')[[1]]
performance_analysis_rdata.file <- get_value.ccb(config_file = user.config.ini.file, key = 'performance_analysis_rdata')[[1]]
external_evaluate_analysis_rdata.file <- get_value.ccb(config_file = user.config.ini.file, key = 'external_evaluate_analysis_rdata')[[1]]

table_cp_stat_list_rdata <- get_value_tripod.ccb(config_file = tripod.ini.file, key = 'table_cp_stat_list_rdata', space = TRUE)[[1]]
table_cp_stat_list_rdata <- gsub(" ","",table_cp_stat_list_rdata)
table_uni_multi_cox_rdata <- get_value_tripod.ccb(config_file = tripod.ini.file, key = 'table_uni_multi_cox')[[1]]
table_uni_multi_cox_rdata <- gsub(" ","",table_uni_multi_cox_rdata)
table_surv_prob_rdata <- get_value_tripod.ccb(config_file = tripod.ini.file, key = 'table_surv_prob_rdata')[[1]]
table_surv_prob_rdata <- gsub(" ","",table_surv_prob_rdata)
user_tripod_rdata <- get_value_tripod.ccb(config_file = tripod.ini.file, key = 'user_tripod_rdata')[[1]]

user_tripod_rdata <- gsub(" ","", user_tripod_rdata)
add_info_tripod_rdata <- paste0(b.user.sig.path,"/rda/add_info_tripod.RData")

# load data ----
load(user.filtered.datasets.rdata)
load(user.parameters.rdata)
load(penalized_cox_analysis_rdata.file)

if(nrow(coef1) > 0) {
  load(independence_analysis_rdata.file)
  load(performance_analysis_rdata.file)
  load(external_evaluate_analysis_rdata.file)
  
  load(table_surv_prob_rdata)
  load(table_uni_multi_cox_rdata)
  load(table_cp_stat_list_rdata)
  load(user_tripod_rdata)
  
  # Generate the items that tripod needs ----
  current.date <- date()
  date1 <- gsub(' ','_',current.date)
  date2 <- gsub(':','_',date1)
  date2 <- gsub('__','_',date2)
  date2_md5 <- digest(date2,algo="crc32")
  b.name.sig.id = paste0(b.name.sig,"_",date2_md5)
  b.name.sig.id = gsub("__","_",b.name.sig.id)
  
  noused_validation_dataset <- names(valid_warning.list)
  used_validation_dataset <- names(test_rs_plot.list)
  used_validation_datasets <- paste0(used_validation_dataset, collapse = ", ")
  used_validation_dataset_pas <- paste0(used_validation_dataset,collapse = ", ")
  
  train_valid_excluded_cp_list <- user_excluded_cp_list[c(select_training_dataset,used_validation_dataset)]
  names(train_valid_excluded_cp_list) <- c("Training dataset",used_validation_dataset)
  
  ##
  train_df <- train.cp.raw
  local_molecular_profiling <- sort(unique(train_df$Molecular_profiling))
  local_min_age <- min(train_df$Age)
  local_max_age <- max(train_df$Age)
  
  # parameters for web
  b.profiling <- sort(unique(user_filtered_datasets_info$molecular_profiling))  
  b.disease.type <- user_filtered_datasets_info %>% 
    select(disease_type) %>% unlist() %>% as.character() %>%
    gsub(" \\(\\d.*\\)", "", .) %>% unique() %>% sort()
  b.sample.type <- sort(unique(user_filtered_datasets_info$sample_type))
  b.primary.site <- sort(unique(user_filtered_datasets_info$primary_site))
  
  ##
  ids <- c("treatment_type", "treatment_setting", "regimen")
  res <- lapply(ids, function(x) {
    user_filtered_datasets_info %>% 
      select(x) %>% unlist() %>% as.character() %>% 
      strsplit(., "<br />") %>% unlist() %>% 
      gsub(" \\(\\d.*\\)", "", .) %>% unique() %>% sort()
  })
  b.prognostic.treatment <- res[[1]]
  b.prognostic.treatment.setting <- res[[2]]
  b.prognostic.regimen <- res[[3]]
  
  
  ##
  train_valid_df_list <- c(list(train_df), valid.cp.list[used_validation_dataset])
  names(train_valid_df_list) <- c("Training dataset", used_validation_dataset)
  train_valid_df_list_names <- names(train_valid_df_list)
  
  m.b.endpoint <- paste0(b.endpoint,"_months")
  b.endpoint.status <- paste0(b.endpoint,"_status")
  if(b.cutpoint.method == "X-tile") {
    b.low2 <- b.low2_xtile
    b.high2 <- b.high2_xtile
    resave(b.low2, b.high2, file = user.parameters.rdata)
  }
  
  ##
  geneset_name = "user provided"
  if(geneset_name == "user provided") {
    geneset_name_text <- "(user provided)"
  } else if (geneset_name == "ColoGuideStageII_genes") {
    geneset_name_text <- "(ColoGuide_Stage_II-Msig candidate genes)"
  } else if (geneset_name == "ColoGuideStageIII_genes") {
    geneset_name_text <- "(ColoGuide_Stage_III-Msig candidate genes)"
  }	else {
    geneset_name_text <- paste0("from MsigDB ",geneset_name," gene set (V7.4)")
  }
  
  ##
  train.profiles_type <- NULL
  predictors_num <- nrow(b.upload.predictors)
  if(b.have.profiles == "Yes") {
    # Signature type
    b.signature.type.pro <- paste0(b.signature.type," molecular signature")
    # Profiling type for training and validation
    train_info <- subset(user_filtered_datasets_info,dataset_name == select_training_dataset)
    train.profiling <- train_info$molecular_profiling
    train.profiles_type <- strsplit(train.profiling," profiling ")[[1]][1]
    train.profiling_type <- tolower(strsplit(train.profiling," profiling ")[[1]][2])
    if(length(used_validation_dataset) > 0){
      valid_info <- subset(user_filtered_datasets_info,dataset_name %in% used_validation_dataset)
      valid.profiling <- valid_info$molecular_profiling
      valid.profiling.list <- list(valid.profiling) 
      names(valid.profiling.list[[1]]) <- used_validation_dataset
    }
  } else {
    # Signature type
    b.signature.type.pro <- paste0(b.signature.type," traditional model")
  }
  sigtype.roc <- "signature group"
  
  
  ##
  if(!is.null(train.profiles_type)){
    if(train.profiles_type == "Protein"){
      predictor_type <- "proteins"
    } else if(train.profiles_type == "Non-coding RNA") {
      predictor_type <- "ncRNA"
    } else {
      predictor_type <- "genes"
    }
  } else {
    predictor_type <- "predictors"
  }
  
  
  if(b.have.profiles == "Yes") {
    candidates_text <- paste0(predictors_num, " candidate ", predictor_type, " ", geneset_name_text, " were used in analysis")
  } else {
    candidates_text <- paste0("Started from ", length(b.clinicopathological.predictors), " candidate clinicopathological variables")
  }
  
  ##
  b.group.cutoff <- "2 Groups"
  if(b.cutpoint.method == "Percentile") { b.group.cutoff <- gsub("groups", "Groups", b.group) }
  
  ##
  predictor_selection_text <- NULL
  if(b.predictor.selection == "Yes") {
    if(b.have.profiles == "Yes") {
      if(b.predictor.selection.methods == "SPCA") {
        predictor_selection_method_text <- "supervised principal components approach (SPCA)"
      } else {						
        predictor_selection_method_text <- "penalized LASSO Cox regression"
      }
    } else {
      predictor_selection_method_text <- "multivariable fractional polynomial (MFP) analysis"
    }
    if(is.na(b.bootstrap.frequency)) {b.bootstrap.frequency = 0}
    predictor_selection_text <- paste0("Selection of robust candidates by a bootstrapping approach (",
                                       b.bootstrap.iterations, " iterations, occurrence > ",
                                       b.bootstrap.frequency, "% frequency) using ",
                                       predictor_selection_method_text
    )
    predictor_selection_text <- paste0("<li><span style = 'color:#34621d'>", predictor_selection_text, "</span></li>")
  }
  
  developed_model_text <- "The final prediction model built by multivariate Cox regression"
  
  signature_generation_text <- paste0("<p><b>Signature generation:&nbsp;</b><ul><li><span style = 'color:#34621d'>",
                                      candidates_text, "</span></li>",
                                      predictor_selection_text,
                                      "<li><span style = 'color:#34621d'>",
                                      developed_model_text, "</span></li></ul></p>"
  )
  
  signature_summary_text <- paste0("<p><b>Signature type:&nbsp;</b><span style = 'color:#34621d'>", 
                                   b.signature.type.pro, "</span></p>", 
                                   signature_generation_text,
                                   "<p><b>Patient risk stratification:&nbsp;</b><span style = 'color:#34621d'>",
                                   b.group.cutoff, "</span></p>"
  )
  
  ##
  if(length(valid_warning.list) == 0) {
    if (length(used_validation_dataset) == 0) {
      signature_datasets_text <- glue::glue("<p><b>Training dataset:&nbsp;</b><span style = 'color:#34621d'>{select_training_dataset}</span></p>")
    } else if (length(used_validation_dataset) == 1) {
      signature_datasets_text <- glue::glue("<p><b>Training dataset:&nbsp;</b><span style = 'color:#34621d'>{select_training_dataset}</span></p><p><b>Validation dataset:&nbsp;</b><span style = 'color:#34621d'>{used_validation_dataset_pas}</span></p>")
    } else {
      signature_datasets_text <- glue::glue("<p><b>Training dataset:&nbsp;</b><span style = 'color:#34621d'>{select_training_dataset}</span></p><p><b>Validation datasets:&nbsp;</b><span style = 'color:#34621d'>{used_validation_dataset_pas}</span></p>")
    }
  } else {
    noused_valid_warning <- paste0("Not all signature variables were exist or levels of some signature variables were not same as training dataset for ",
                                   paste0(noused_validation_dataset, collapse = ", "))
    if(length(used_validation_dataset) == 0) {
      signature_datasets_text <- glue::glue("<p><b>Training dataset:&nbsp;</b><span style = 'color:#34621d'>{select_training_dataset}</span></p><p><b>Validation dataset:&nbsp;</b><span style = 'color:#34621d'>No validation dataset</span></p><p><b>Warning:&nbsp;</b><span style = 'color:#34621d'>{noused_valid_warning}</span></p>")
    } else {
      signature_datasets_text <- glue::glue("<p><b>Training dataset:&nbsp;</b><span style = 'color:#34621d'>{select_training_dataset}</span></p><p><b>Validation datasets:&nbsp;</b><span style = 'color:#34621d'>{used_validation_dataset_pas}</span></p><p><b>Warning:&nbsp;</b><span style = 'color:#34621d'>{noused_valid_warning}</span></p>")
    }
  }
  
  # I0. TRIPOD title -------------------------------------
  tripod_generator <- paste0("<b>Lead contact:</b> ",generator)
  tripod_email <- paste0("<b>Email address:</b> ",email)
  tripod_organization <- paste0("<b>Organization:</b> ",organization)
  create_date <- Sys.Date()
  tripod_create_date <- paste0("Date created: ",create_date)
  
  # I. Signature summary -------------------------------------
  signature_summary_texts <- paste0("<p><b>Signature name:&nbsp;</b><span style = 'color:#34621d'>", b.name.sig,"</span></p>", signature_summary_text)
  
  # Dataset summary
  train_valid_df_list_names0 <- gsub("Training dataset", select_training_dataset, train_valid_df_list_names)
  
  user_filtered_datasets_info0 <- user_filtered_datasets_info %>%
    filter(dataset_name %in% c(select_training_dataset, used_validation_dataset)) %>%
    mutate(numberpts = do.call(rbind, strsplit(numberpts, " / "))[,1] %>% as.numeric()) %>% 
    mutate(endpoint = do.call(rbind, strsplit(endpoint, "\\("))[,2] %>%
             gsub(")", ", ", .) %>%
             paste0(b.endpoint, " (", ., numberpts, " patients)"))
  if(b.have.profiles == "Yes") {
    user_filtered_datasets_info0 <- user_filtered_datasets_info0 %>%
      mutate(molecular_profiling = tolower(gsub("profiling", "profiles", molecular_profiling))) %>%
      #select(dataset_id, geo_url, dataset_name,
      select(dataset_id, dataset_name,
             numberpts, primary_site, disease_type,
             molecular_profiling, endpoint)
  } else {
    user_filtered_datasets_info0 <- user_filtered_datasets_info0 %>%
      #select(dataset_id, geo_url,
      select(dataset_id,
             dataset_name, numberpts,
             primary_site, disease_type,
             endpoint)
  }
  
  ##
  tripod_dataset_infos <- data.frame()
  ids <- names(conf$Signature_summary)[-which(names(conf$Signature_summary) %in% "Signature_description")]
  nullToNA <- function(x) { 
    x[is.null(x)] <- NA
    return(x) 
  } 
  
  for(x in ids) {
    list0 <- conf$Signature_summary[[x]]
    name0 <- list0$dataset_name
    ##
    if(!is.null(name0)) {
      meta_infor <- data.frame("dataset_name" = name0, "pubmed_id" = nullToNA(list0$pubmed_id), 
                               "accession_number" = nullToNA(list0$accession_number), "accession_url" = nullToNA(list0$accession_url), 
                               "clinical_number" = nullToNA(list0$clinical_number), "clinical_url" = nullToNA(list0$clinical_url))
      tripod_dataset_infos <- rbind(tripod_dataset_infos, meta_infor)
    }
  } 
  
  ##
  tripod_datasets_summary_df <- get_tripod_dataset_summary(tripod_dataset_infos = tripod_dataset_infos,
                                                           user_filtered_datasets_info0 = user_filtered_datasets_info0,
                                                           b.have.profiles = b.have.profiles,
                                                           select_training_dataset = select_training_dataset,
                                                           used_validation_dataset = used_validation_dataset)
  
  
  save(create_date, tripod_create_date, train_valid_df_list,
       signature_summary_texts, signature_summary_text, train.profiles_type,
       predictor_type, generator, email, organization, signature_description,
       b.signature.type.pro, b.group.cutoff, km_input.list, valid_KM_prex.list,
       b.name.sig.id, b.name.sig, user_filtered_datasets_info0, tripod_dataset_infos,
       train_valid_df_list_names, train_valid_df_list_names0, 
       tripod_datasets_summary_df,
       used_validation_dataset,
       signature_datasets_text, 
       file = add_info_tripod_rdata)
  
  
  # for pearson correlation ----
  profile_type <- ""
  if(!is.null(train.profiles_type)) {
    if(train.profiles_type == "Methylation") {
      profile_type <- "methylation"
    } else if(train.profiles_type %in% c("Somatic mutations", "Germline mutations")) {
      profile_type <- ""
    } else {
      profile_type <- "expression"
    }
  }
  
  # II. Items of methods -------------------------------------
  # Initial values for items of methods
  # Item5b (web)
  user_excluded_info_list <- lapply(names(train_valid_excluded_cp_list), function(data_name) {
    df <- train_valid_excluded_cp_list[[data_name]]
    sample_num <- nrow(df)
    
    if(nrow(df) == 0){
      exclusion_web <- paste0(data_name,": No patients were excluded by ReProMSig.")
    }else {
      excluded.primary.site <- excluded_sample_info(df=df,column_name="Primary_site",input=b.primary.site)
      excluded.disease.type <- excluded_sample_info(df=df,column_name="Disease_type",input=b.disease.type)
      excluded.sample.type <- excluded_sample_info(df=df,column_name="Sample_type",input=b.sample.type)
      excluded.stage.edition <- excluded_sample_info(df=df,column_name="Stage_edition",input=b.stage.edition)
      excluded.stage <- excluded_sample_info(df=df,column_name="Stage",input=b.stage)
      excluded.T <- excluded_sample_info(df=df,column_name="T",input=b.T)
      excluded.N <- excluded_sample_info(df=df,column_name="N",input=b.N)
      excluded.M <- excluded_sample_info(df=df,column_name="M",input=b.M)
      excluded.gender <- excluded_sample_info(df=df,column_name="Gender",input=b.Gender)
      excluded.ethnicity <- excluded_sample_info(df=df,column_name="Ethnicity",input=b.ethnicity)
      
      if(b.have.profiles == "Yes"){
        excluded.profiling <- excluded_sample_info(df=df,column_name="Molecular_profiling",input=b.profiling)
      } else{
        excluded.profiling <- ""
      }
      excluded.is.primary.treatment <- excluded_sample_info(df=df,column_name="Is_primary_treatment",input=b.is.primary.treatment)
      if(b.signature.type=="Prognostic"){
        excluded.treatment.type <- excluded_sample_info(df=df,column_name="Treatment_type",input=b.prognostic.treatment)
        excluded.treatment.setting <- excluded_sample_info(df=df,column_name="Treatment_setting",input=b.prognostic.treatment.setting)
        excluded.regimen <- ""  
      } else {
        excluded.treatment.type <- excluded_sample_info(df=df,column_name="Treatment_type",input=b.predictive.treatment)
        excluded.treatment.setting <- excluded_sample_info(df=df,column_name="Treatment_setting",input=b.predictive.treatment.setting)
        excluded.regimen <- excluded_sample_info(df=df,column_name="Regimen",input=b.predictive.regimen)
      }
      
      # age
      if(!is.na(b.min.age)){
        df$Age <- as.numeric(df$Age)
        excluded.max.age <- max(df$Age[!is.na(df$Age)])
        if(excluded.max.age > b.min.age){
          excluded.min.age = paste0("age > ",b.min.age)
        } else {
          excluded.min.age =""
        }
      } else{
        excluded.min.age =""
      }
      
      # outcome
      df[,m.b.endpoint] <- as.numeric(df[, m.b.endpoint])
      df[,b.endpoint.status] <- as.numeric(df[, b.endpoint.status])
      missing_outcome_sample <- c(df[df[,m.b.endpoint][is.na(df[,m.b.endpoint])], "Sample_ID"],
                                  df[df[,b.endpoint.status][is.na(df[,b.endpoint.status])], "Sample_ID"])
      missing_outcome_sample <- unique(missing_outcome_sample)
      if(length(missing_outcome_sample) == 0){
        missing_outcome <- ""
      } else {
        missing_outcome <- "- missing outcome"
      }
      
      excluded.info <- c(excluded.primary.site, excluded.disease.type,
                         excluded.sample.type,
                         excluded.stage.edition, excluded.stage,
                         excluded.T, excluded.N, excluded.M,
                         excluded.gender, excluded.ethnicity,
                         excluded.is.primary.treatment,
                         excluded.treatment.type,
                         excluded.treatment.setting, excluded.regimen,
                         excluded.min.age, missing_outcome)
      
      excluded.infos <- excluded.info[!excluded.info %in% ""]
      excluded.infos <- paste0(excluded.infos,collapse="<br/ >")
      exclusion_web <- paste0(data_name,": ", sample_num, 
                              " patients meeting the following criteria were excluded from the ReProMSig.<br/ >",
                              excluded.infos)
    }
    return(exclusion_web)
  })
  
  names(user_excluded_info_list) <- names(train_valid_excluded_cp_list)
  
  user_included_info_list <- lapply(names(train_valid_df_list), function(data_name) {
    df <- train_valid_df_list[[data_name]]
    excluded_df <- train_valid_excluded_cp_list[[data_name]]
    sample_num <- nrow(df)
    
    if(nrow(excluded_df) == 0){
      inclusion_web <- paste0(data_name,": All patients were included in ReProMSig.")
    } else {
      included.primary.site <- included_sample_info(df=df,column_name="Primary_site",input=b.primary.site)
      included.disease.type <- included_sample_info(df=df,column_name="Disease_type",input=b.disease.type)
      included.sample.type <- included_sample_info(df=df,column_name="Sample_type",input=b.sample.type)
      included.stage.edition <- included_sample_info(df=df,column_name="Stage_edition",input=b.stage.edition)
      included.stage <- included_sample_info(df=df,column_name="Stage",input=b.stage)
      included.T <- included_sample_info(df=df,column_name="T",input=b.T)
      included.N <- included_sample_info(df=df,column_name="N",input=b.N)
      included.M <- included_sample_info(df=df,column_name="M",input=b.M)
      included.gender <- included_sample_info(df=df,column_name="Gender",input=b.Gender)
      included.ethnicity <- included_sample_info(df=df,column_name="Ethnicity",input=b.ethnicity)
      if(b.have.profiles == "Yes"){
        included.profiling <- included_sample_info(df=df,column_name="Molecular_profiling",input=b.profiling)
      } else{
        included.profiling <- ""
      }
      included.is.primary.treatment <- included_sample_info(df=df,column_name="Is_primary_treatment",input=b.is.primary.treatment)
      
      if(b.signature.type=="Prognostic"){
        included.treatment.type <- included_sample_info(df=df,column_name="Treatment_type",input=b.prognostic.treatment)
        included.treatment.setting <- included_sample_info(df=df,column_name="Treatment_setting",input=b.prognostic.treatment.setting)
        included.regimen <- ""
      } else {
        included.treatment.type <- included_sample_info(df=df,column_name="Treatment_type",input=b.predictive.treatment)
        included.treatment.setting <- included_sample_info(df=df,column_name="Treatment_setting",input=b.predictive.treatment.setting)
        included.regimen <- included_sample_info(df=df,column_name="Regimen",input=b.predictive.regimen)
      }
      
      # age
      if(!is.na(b.min.age)){
        included.min.age = paste0("age >=",b.min.age)
      } else{
        included.min.age = ""
      }
      
      included.info <- c(included.primary.site,included.disease.type,included.sample.type,
                         included.stage.edition,included.stage,included.T,included.N,included.M,
                         included.gender,included.ethnicity,included.is.primary.treatment,
                         included.treatment.type,included.treatment.setting,included.regimen,
                         included.min.age)
      
      included.infos <- included.info[!included.info %in% ""]
      included.infos <- paste0(included.infos,collapse="<br/ >")
      inclusion_web <- paste0(data_name,": ",sample_num,
                              " patients meeting the following criteria were included in ReProMSig.<br/ >",
                              included.infos)
    }
    return(inclusion_web)
  })
  names(user_included_info_list) <- names(train_valid_df_list)
  
  # Items for training dataset
  item5b_training_exclusion_web <- user_excluded_info_list[["Training dataset"]]
  item5b_training_exclusion_web <- gsub("Training dataset: ","",item5b_training_exclusion_web)
  item5b_training_inclusion_web <- user_included_info_list[["Training dataset"]]
  item5b_training_inclusion_web <- gsub("Training dataset: ","",item5b_training_inclusion_web)
  
  # Items for validation
  if(length(used_validation_dataset) > 0) {
    valid_excluded_info_list <-user_excluded_info_list[used_validation_dataset]
    item5b_validation_exclusion_web <- do.call(paste,c(valid_excluded_info_list,sep = "<br /><br /> "))
    valid_included_info_list <-user_included_info_list[used_validation_dataset]
    item5b_validation_inclusion_web <- do.call(paste,c(valid_included_info_list,sep = "<br /><br /> "))
    
    if(b.have.profiles == "Yes") {  
      item9_validation_methods <- paste0("Assuming occurrence of missing values in ", 
                                         predictor_type,
                                         " profiles were random events, KNN method (impute.knn function in R package impute version 1.60.0) was applied to impute missing values by ReProMSig.<br />Patients with missing outcome or with >=10% missing clinicopathological variables, were excluded from ReProMSig analysis."
      )  
    } else {
      item9_validation_methods <- "Patients with missing outcome or with >=10% missing clinicopathological variables, were excluded from ReProMSig analysis."
    }
  } else {
    item4a_validation_source <- ""
    item4a_validation_study_date <- "" 
    item5a_validation_setting <- ""
    item5b_validation_exclusion_user <- ""
    item5b_validation_exclusion_web <- ""
    item5b_validation_inclusion_user <- ""
    item5b_validation_inclusion_web <- ""
    item5c_validation <- ""
    item9_validation_methods <- ""
  }
  
  if(uni_cox_variables == "" || length(uni_cox_variables) == 0) {
    uni_variables_pas <- paste0("The molecular signature prediction (risk groups) was included in the independence test.")
  } else {
    if(length(uni_cox_variables) == 1) {
      uni_variables_pas <- paste0("Both the molecular signature prediction (risk groups) and clinicopathological variable (", 
                                  paste0(uni_cox_variables,collapse=", "), 
                                  ") were included in the independence test.")
    } else {
      uni_variables_pas <- paste0("Both the molecular signature prediction (risk groups) and clinicopathological variables (",
                                  paste0(uni_cox_variables,collapse=", "),
                                  ") were included in the independence test.")
    }
  }
  
  if(multi_cox_variables == "" || length(multi_cox_variables) == 0) {
    multi_variables_pas <- ""
  } else {
    if(length(uni_cox_variables) == 1) {
      multi_variables_pas <- paste0("In a multivariate analysis, the interaction between predicted treatment-benefit groups and treatment were adjusted for other clinicopathological variable (", paste0(multi_cox_variables,collapse=", "),"), demonstrating persistence of the strength of the interaction. We also examined the interactions between clinicopathological variable and treatment.")
    } else {
      multi_variables_pas <- paste0("In a multivariate analysis, the interaction between predicted treatment-benefit groups and treatment were adjusted for other clinicopathological variables (",paste0(multi_cox_variables,collapse=", "),"), demonstrating persistence of the strength of the interaction. We also examined the interactions between clinicopathological variables and treatment.")
    }
  }
  
  ##
  num0 <- ncol(train_multi_sig_nom)
  if(num0 >= 4){
    nomogram_variables_pas <- paste0(nomogram_variables,collapse=", ")
  } else {
    nomogram_variables_pas <- ""
  }
  
  # molecular / traditional
  if(b.have.profiles == "Yes") {
    if(ROC_variables == "" || length(ROC_variables) == 0) {
      if(ncol(train_multi_sig_nom) >= 4) {
        roc_variables_pas <- paste0("ROC curves for molecular signature and the combined model were plotted.")
      } else {
        roc_variables_pas <- paste0("ROC curves for molecular signature was plotted.")
      }
    } else {
      if(length(ROC_variables) == 1) {
        if(ncol(train_multi_sig_nom) >= 4) {
          roc_variables_pas <- paste0("ROC curves for molecular signature, clinicopathological variable (", 
                                      paste0(ROC_variables,collapse=", "), 
                                      ") and the combined model were plotted.")
        } else {
          roc_variables_pas <- paste0("ROC curves for molecular signature and clinicopathological variable (", 
                                      paste0(ROC_variables,collapse=", "),
                                      ") were plotted.")
        }
      } else {
        if(ncol(train_multi_sig_nom) >= 4) {
          roc_variables_pas <- paste0("ROC curves for molecular signature, clinicopathological variables (",
                                      paste0(ROC_variables,collapse=", "),
                                      ") and the combined model were plotted.")
        } else {
          roc_variables_pas <- paste0("ROC curves for molecular signature and clinicopathological variables (",
                                      paste0(ROC_variables,collapse=", "),
                                      ") were plotted.")
        }
      }
    }
    
    if(PE_variables == "" || length(PE_variables) == 0) {
      if(ncol(train_multi_sig_nom) >= 4) {
        pe_variables_pas <- paste0("PE curves for molecular signature and the combined model were plotted.")
      } else {
        pe_variables_pas <- paste0("PE curves for molecular signature was plotted.")
      }
    } else {
      if(length(PE_variables) == 1) {
        if(ncol(train_multi_sig_nom) >= 4) {
          pe_variables_pas <- paste0("PE curves for molecular signature, clinicopathological variable (",
                                     paste0(PE_variables,collapse=", "),
                                     ") and the combined model were plotted.")
        } else {
          pe_variables_pas <- paste0("PE curves for molecular signature and clinicopathological variable (", 
                                     paste0(PE_variables,collapse=", "),
                                     ") were plotted.")
        }
      } else {
        if(ncol(train_multi_sig_nom) >= 4) {
          pe_variables_pas <- paste0("PE curves for molecular signature, clinicopathological variables (",
                                     paste0(PE_variables,collapse=", "),
                                     ") and the combined model were plotted.")
        } else {
          pe_variables_pas <- paste0("PE curves for molecular signature, clinicopathological variables (",
                                     paste0(PE_variables,collapse=", "),
                                     ") were plotted.")
        }
      }
    }
    
    # add feature num
    item9_training_methods <- paste0("Assuming occurrence of missing values in ", 
                                     predictor_type,
                                     " were random events, KNN method (impute.knn function in R package impute version 1.60.0) was applied to impute missing values by ReProMSig.\nPatients with missing outcome or with >=10% missing clinicopathological variables, were excluded from ReProMSig analysis.")
    
    item10a_handling_molecular0 <- paste0(train.profiles_type," values were converted to log2 transformed.")
    if(!is.null(b.batch.correction)){
      item10a_handling_molecular1 <- "Using training dataset gene expression values as the reference, ComBat function in R package \"sva\" (version 3.34.0) was applied to reduce the likelihood of batch effects from nonbiological technical biases for each validation dataset profiles."
    } else {
      item10a_handling_molecular1 <- ""
    }
    item10a_handling_molecular <- paste(item10a_handling_molecular0,item10a_handling_molecular1)
    item10a_implausible <- paste0("Extreme values (e.g., outliers) in ", 
                                  tolower(train.profiles_type),
                                  " profiles were regarded as NA, and imputed by KNN method, using impute.knn function in R package \"impute \" (version 1.60.0).")
    
    predictors_num <- nrow(b.upload.predictors)
    predictors <- as.character(b.upload.predictors[,1])
    item10b_before_genes0 <- ""
    
    if(train.profiles_type == "Expression") {
      if(geneset_name == "user provided") {
        item10b_before_genes0 <- paste0("The following genes (n=",predictors_num,") from XXX analysis were used for predictor selection and modelling.")
      } else if (geneset_name == "ColoGuideStageII_genes") {
        item10b_before_genes0 <- paste0("The following genes (n=",predictors_num,
                                        ") from ColoGuide_Stage_II-Msig candidate genes were used for predictor selection and modelling.") 
      } else if(geneset_name == "ColoGuideStageIII_genes") {
        item10b_before_genes0 <- paste0("The following genes (n=",predictors_num,
                                        ") from ColoGuide_Stage_III-Msig candidate genes were used for predictor selection and modelling.")
      } else {
        item10b_before_genes0 <- paste0("The following genes (n=",predictors_num,") from ",
                                        paste0("MsigDB ",geneset_name," gene set (V7.4)"),
                                        " were used for predictor selection and modelling.")
      }
    } else {
      if(geneset_name == "user provided") {
        item10b_before_genes0 <- paste0("The following predictors (n=",predictors_num,") from XXX analysis were used for predictor selection and modelling.")
      } else {
        item10b_before_genes0 <- paste0("The following predictors (n=",predictors_num,") from ",
                                        paste0("MsigDB ",geneset_name,
                                               " gene set (V7.4)")," were used for predictor selection and modelling.")
      }
    }
    
    item10b_before_genes1 <- paste0(predictors,collapse = " ")
    item10b_before_genes <- paste0(item10b_before_genes0,"\n",item10b_before_genes1)
    
    if(b.predictor.selection == "Yes") {
      if(b.predictor.selection.methods == "SPCA") {
        item10b_before.method <- "supervised principal components approach (SPCA, R package \"SPCA\")"
      } else {						
        item10b_before.method <- "LASSO penalized Cox regression (R package \"glmnet\")"
      }
      item10b_before <- paste0("A bootstrapping procedure was applied to select molecular predictors associated with outcome from the list of candidate predictors shown above. For each bootstrap resample (samples were drawn with replacement keeping the same size as the original training dataset), ",item10b_before.method," was applied. The entire process was repeated ", b.bootstrap.iterations," times, and predictors with frequency of more than ",b.bootstrap.frequency,"% occurrence were taken as potential signature predictors for prognosis.")
    }
    
    ##
    if(LASSO_cox) {
      item10b_model <- "LASSO penalized Cox regression (R package \"glmnet\") was generated using the selected molecular predictors, which optimum penalty parameter \"lambda\" was chosen using ten-fold cross-validation to minimise the mean cross-validation partial likelihood error rate."
      item10c_signature <- "Individual signature score was calculated by a weighted sum of the predictors in the generated LASSO penalized Cox regression model, in which weights are the corresponding regression coefficients."
    } 
    if(standard_cox) {
      item10b_model <- paste0("Univariate Cox regression model (R package \"survival\") was ﬁtted to the expression data (training dataset) of the selected molecular predictors, to identify robust molecular signatures (using 0.05 as the significance level) associated with ",b.endpoint,", which were used to generate the multivariable prediction signature (R package \"survival\").")
      item10c_signature <- "Individual signature score was calculated by a weighted sum of the predictors in the generated multivariate Cox regression model, in which weights are the corresponding regression coefficients (i.e., log hazard ratios). The signature score for a patient is the log relative hazard compared with a hypothetical patient whose signature score is zero."
    }
    
    if(b.nomogram.option == "user"){
      nomogram.option = "user provided"
    } else {
      nomogram.option = "identified by independence test"
    }
    
    if(ncol(train_multi_sig_nom) >= 4){
      clinical_factors0 <- colnames(train_multi_sig_nom)[!colnames(train_multi_sig_nom) %in% c("time","status","Signature_score","Signature_group","Predictive_group")]
      clinical_factors <- paste0(clinical_factors0,collapse=", ")
      
      item10c_nomogram <- paste0("Multivariate Cox regression model integrating the molecular signature prediction (",sigtype.roc,") and clinicopathological factors (",clinical_factors,") (",nomogram.option,"), was developed for nomograms. The nomogram was then created using patients in training dataset, which could be used to estimate the ",b.endpoint," probabilities at ",paste0(b.times,collapse=", "),"-months for single patient.")
      
      item10c_online <- paste0("An online research tool for single patient prediction of risk score (i.e., the exponential of signature score) and ",b.endpoint," probabilities at ",paste0(b.times,collapse=", "),"-months using the molecular signature is available from https://omics.bjcancer.org/prognosis/.")
      
      item10d_calibration <- paste0("The calibration plots at ",paste0(b.times,collapse=", "),"-months were used to assess the consistency between the actual and predicted ",b.endpoint," probabilities from the molecular signature and the combined model in training dataset.")
    } else {
      item10c_nomogram <- "Nomogram was not constructed as no clinicopathological variable was selected in the development form."
      
      item10c_online <- paste0("An online research tool for single patient prediction of risk score (i.e., the exponential of signature score) and ",b.endpoint," probabilities at ",paste0(b.times,collapse=", "),"-months using the molecular signature is available from https://omics.bjcancer.org/prognosis/.")
      
      item10d_calibration <- paste0("The calibration plots at ",paste0(b.times,collapse=", "),"-months were used to assess the consistency between the actual and predicted ",b.endpoint," probabilities from the molecular signature in training dataset.")
    }
    
    item10d_roc <- paste0(paste0(b.times,collapse=", "),"-months time-dependent ROC analysis was performed to examine the prognostic accuracy of the developed model in the training dataset. ",roc_variables_pas," An area under the ROC curve (AUC) of 0.5 indicates no discrimination, whereas an AUC of 1.0 indicates perfect discrimination.")
    
    item10d_pe <- paste0("PE curve analysis was applied in the training dataset to examine the prognosis prediction error rate, and ten-fold cross-validation cumulative prediction error was computed using Kaplan-Meier estimation as reference. Models with smaller area under the curve indicates a relatively lower error rate. ",pe_variables_pas)
    
  } else {
    if(ROC_variables == "" || length(ROC_variables) == 0) {
      roc_variables_pas <- paste0("ROC curve for the prediction model was plotted.")
    } else {
      if(length(ROC_variables) == 1) {
        roc_variables_pas <- paste0("ROC curves for the prediction model, clinicopathological variable (",paste0(ROC_variables,collapse=", "),") were plotted.")
      } else {
        roc_variables_pas <- paste0("ROC curves for the prediction model, clinicopathological variables (",paste0(ROC_variables,collapse=", "),") were plotted.")
      }
    }
    
    if(PE_variables == "" || length(PE_variables) == 0) {
      pe_variables_pas <- paste0("PE curve for the prediction model was plotted.")
    } else {
      if(length(PE_variables) == 1) {
        pe_variables_pas <- paste0("PE curves for the prediction model, clinicopathological variable (",paste0(PE_variables,collapse=", "),") were plotted.")
      } else {
        pe_variables_pas <- paste0("PE curves for the prediction model, clinicopathological variables (",paste0(PE_variables,collapse=", "),") were plotted.")
      }
    }
    
    item9_training_methods <- "Patients with missing outcome or with >=10% missing clinicopathological variables, were excluded from ReProMSig analysis."
    
    if("Stage" %in% b.combine.variables){
      stage1 <- paste0(b.combine.stage1,collapse ="+")
      stage2 <- paste0(b.combine.stage2,collapse ="+")
      item10a.stage <- paste0("- Stage was merged into two categories: ",stage1," and ",stage2)
    } else {
      item10a.stage <- ""
    }
    
    if("T" %in% b.combine.variables){
      t1 <- paste0(b.combine.t1,collapse ="+")
      t2 <- paste0(b.combine.t2,collapse ="+")
      item10a.t <- paste0("- T stage was merged into two categories: ",t1," and ",t2)
    } else {
      item10a.t <- ""
    }
    if("N" %in% b.combine.variables){
      n1 <- paste0(b.combine.n1,collapse ="+")
      n2 <- paste0(b.combine.n2,collapse ="+")
      item10a.n <- paste0("- N stage was merged into two categories: ",n1," and ",n2)
    } else {
      item10a.n <- ""
    }
    if("Age" %in% b.combine.variables){
      cutoff = as.numeric(b.combine.age)
      item10a.age <- paste0("- Age was dichotomized into two groups using",cutoff,"as the cutoff.")
    } else {
      item10a.age <- ""
    }
    mfp_continuous_variables_pas <- paste0(mfp_continuous_variables, collapse = ", ")
    mfp_transform <- paste0("After assessment of nonlinearity for continuous clinicopathological variables, multivariable fractional polynomial (MFP) approach was applied to (",mfp_continuous_variables_pas,"), to achieve a good approximation of linear relationship with outcome after MFP transformation.")
    
    item10a.combine <- c(item10a.stage,item10a.t,item10a.n,item10a.age,mfp_transform)
    item10a.combine <- item10a.combine[which(!item10a.combine %in% "")]
    item10a_handling_clinical <- paste0(item10a.combine,collapse = "<br />")
    
    item10b_before <- paste0("A bootstrapping procedure was applied to select clinicopathological predictors associated with outcome. Multivariable fractional polynomial (MFP) approach with backward elimination at the 5% significance level was used to selection of features, along with exploring the presence and transformation of nonlinear relationships with outcome for the continuous predictors. For each bootstrap resample (samples were drawn with replacement keeping the same size as the original training dataset), MFP model were created. The entire process was repeated ",b.bootstrap.iterations," times, and predictors with frequency of more than ",b.bootstrap.frequency," occurrence were taken as potential signature predictors for prognosis.")
    
    item10b_model <- paste0("Univariate Cox regression model (R package \"survival\") was ﬁtted to the training dataset using the selected clinicopathological predictors, to identify robust predictors (using 0.05 as the significance level) associated with ",b.endpoint,", which were used to generate the final multivariable Cox regression model (R package \"survival\").")
    
    item10c_signature <- "Individual signature score was calculated by a weighted sum of the predictors in the generated multivariate Cox regression model, in which weights are the corresponding regression coefficients (i.e., log hazard ratios). The signature score for a patient is the log relative hazard compared with a hypothetical patient whose signature score is zero."
    
    clinical_factors0 <- colnames(train_multi_sig_nom)[!colnames(train_multi_sig_nom) %in% c("time","status","Signature_score","Signature_group","Predictive_group")]
    clinical_factors <- paste0(clinical_factors0,collapse=", ")
    
    item10c_nomogram <- paste0("Multivariate Cox regression model integrating predictors (",clinical_factors,") of the prediction model was developed for nomogram. The nomogram was then created using patients in training dataset, which could be used to estimate the ",b.endpoint," probabilities ",paste0(b.times,collapse=", "),"-months for single patient.")
    
    item10c_online <- paste0("An online research tool for single patient prediction of risk score (i.e., the exponential of signature score) and ",b.endpoint," probabilities at ",paste0(b.times,collapse=", "),"-months using the prediction model is available from https://omics.bjcancer.org/prognosis/")
    
    item10d_roc <- paste0(paste0(b.times,collapse=", "),"-months time-dependent ROC analysis was performed to examine the prognostic accuracy of the developed model in the training dataset. ",roc_variables_pas," An area under the ROC curve (AUC) of 0.5 indicates no discrimination, whereas an AUC of 1.0 indicates perfect discrimination.")
    
    item10d_pe <- paste0("PE curve analysis was applied in the training dataset to examine the prognosis prediction error rate, and ten-fold cross-validation cumulative prediction error was computed using Kaplan-Meier estimation as reference. Models with smaller area under the curve indicates a relatively lower error rate. ",pe_variables_pas)
    
    item10d_calibration <- paste0("The calibration plots at ",paste0(b.times,collapse=", "),"-months were used to assess the consistency between the actual and predicted ",b.endpoint," probabilities from the prognosis model in training dataset.")
  }
  
  # prognostic / predictive
  if(length(used_validation_dataset) == 0) {
    km_dataset_text <- "in the training dataset"
  } else if(length(used_validation_dataset) == 1) {
    km_dataset_text <- "in the training and validation dataset respectively"
  } else {
    km_dataset_text <- "in the training and validation datasets respectively"
  }
  if(b.signature.type == "Prognostic") {
    if(multi_cox_variables == "" || length(multi_cox_variables) == 0) {
      item10d_independence <- paste0("Univariate Cox regression analysis was performed to test the association between the signature and ",b.endpoint,".")
    } else {
      item10d_independence <- paste0("Univariate and multivariate Cox regression analyses were performed to test whether individual prognosis factor is an independent factor in predicting patients ",b.endpoint,". ",uni_variables_pas)
    }
    
    item10d_interaction <- "We did not examine interaction terms but relied on the main effects of the selected."
    if(b.group == "2 groups") {
      item10d_km <- paste0("Discrimination is visually inspected from the spread of Kaplan-Meier curves for each predicted risk group (high-risk, and low-risk), ",km_dataset_text,". Differences in the probability of ",b.endpoint," between risk groups (high risk vs low risk) were tested by the two-sided log-rank test. In addition, Kaplan-Meier plots also presented the total number of patients, the number of events (outcome) for each risk group.")
    } else {
      item10d_km <- paste0("Discrimination is visually inspected from the spread of Kaplan-Meier curves for each predicted risk group (high-risk, intermediate risk and low-risk), ",km_dataset_text,". Differences in the probability of ",b.endpoint," between risk groups (intermediate risk vs low risk, and high risk vs low risk) were tested by the two-sided log-rank test. In addition, Kaplan-Meier plots also presented the total number of patients, the number of events (outcome) for each risk group.")
    }
    
    item10d_km_pro <- paste0(b.endpoint," with 95% confidence intervals at ",paste0(b.times,collapse=", "),"-months ",km_dataset_text,", were calculated by Kaplan-Meier method for patients in each risk group.")
    
    if(!is.null(b.cp.variables)) {
      item10d_association <- "To assess whether the developed prediction model is correlated with clinicopathological variables, two-sided Fisher's exact test (for categorical variables) and two-sided t test (for continuous variables) were applied to measure association or difference between predicted risk groups."
    } else {
      item10d_association <- ""
    }
  } else {
    if(multi_cox_variables == "" || length(multi_cox_variables) == 0) {
      item10d_interaction <- paste0("Using the training cohort, univariate Cox regression analysis was performed to test the association between the signature and ",b.endpoint)
    } else {
      item10d_interaction <- paste0("Using the training cohort, univariate and multivariate Cox regression analyses were performed to test for interaction between treatment and the predicted treatment-benefit groups (treatment benefit and no-benefit). ",multi_variables_pas," Furthermore, we conducted a subgroup analysis, for each subgroup (e.g., treatment-benefit group), the difference in outcome between patients with and without treatment was tested by log-rank test, to examine if patients of this group received significant treatment benefit compared to those who did not underwent treatment. P value less than 0.05 by the likelihood ratio test was considered signiﬁcant.")
    }
    
    if(b.group == "2 groups") {
      item10d_km <- paste0("Discrimination is visually inspected from the spread of Kaplan-Meier curves for each predicted treatment-benefit group (benefit, and no-benefit; treatment and control in benefit group; treatment and control in no-benefit group), ",km_dataset_text,". Differences in the probability of ",b.endpoint," between risk groups (benefit vs no-benefit; treatment vs control in benefit group; treatment vs control in no-benefit group) were tested by the two-sided log-rank test. In addition, Kaplan-Meier plots also presented the total number of patients, the number of events (outcome) for each treatment-benefit group.")
    } else {
      item10d_km <- paste0("Discrimination is visually inspected from the spread of Kaplan-Meier curves for each predicted treatment-benefit group (benefit, intermediate group, and no-benefit; treatment and control in benefit group; treatment and control in intermediate group; treatment and control in no-benefit group), ",km_dataset_text,". Differences in the probability of ",b.endpoint," between risk groups (benefit vs no-benefit; intermediate group vs no-benefit; treatment vs control in benefit group; treatment vs control in intermediate group; treatment vs control in no-benefit group) were tested by the two-sided log-rank test. In addition, Kaplan-Meier plots also presented the total number of patients, the number of events (outcome) for each treatment-benefit group.")
    }
    
    item10d_km_pro <- paste0(b.endpoint," with 95% confidence intervals at ",paste0(b.times,collapse=", "),"-months ",km_dataset_text,", were calculated by Kaplan-Meier method for patients in each treatment-benefit group.")
    
    if(!is.null(b.cp.variables)) {
      item10d_association <- "To assess whether the developed prediction model is correlated with clinicopathological variables, two-sided Fisher's exact test (for categorical variables) and two-sided t test (for continuous variables) were applied to measure association or difference between predicted treatment-benefit groups."
    } else {
      item10d_association <- ""
    }
  }
  
  
  ##
  if(b.group == "2 groups") {
    if(b.cutpoint.method == "Percentile") {
      b.high2.cutoff = b.low2.cutoff = paste0(b.high2,"th percentile")
    } else {
      b.high2.cutoff = b.low2.cutoff = round(top_train,3)
    }
    if(length(used_validation_dataset) == 0) {
      item11_groups <- paste0("Patients were stratified into two risk groups on the basis of signature score distribution: low-risk (signature score < ",b.low2.cutoff,") and high-risk (signature score >= ",b.high2.cutoff,").")
    } else {
      item11_groups <- paste0("Patients were stratified into two risk groups on the basis of signature score distribution in each dataset: low-risk (signature score < ",b.low2.cutoff,") and high-risk (signature score >= ",b.high2.cutoff,").")
    }
  } else {
    if(length(used_validation_dataset) == 0) {
      item11_groups <- paste0("Patients were stratified into three risk groups on the basis of signature score distribution: low-risk (signature score < ",b.low3,"th percentile), intermediate-risk (signature score >= ",b.low3,"th and signature score < ",b.high3,"th percentile) and high-risk (signature score >= ",b.high3,"th percentile).")
    } else {
      item11_groups <- paste0("Patients were stratified into three risk groups on the basis of signature score distribution in each dataset: low-risk (signature score < ",b.low3,"th percentile), intermediate-risk (signature score >= ",b.low3,"th and signature score < ",b.high3,"th percentile) and high-risk (signature score >= ",b.high3,"th percentile).")
    }
  }
  
  # Item4a,13a ----
  da_followups <- c()
  for(i in train_valid_df_list_names) {
    dacp <- train_valid_df_list[[i]]
    
    # followup
    ifelse(i == "Training dataset", i_type <- "Training dataset" , i_type <- "Validation dataset")
    if(i == "Training dataset"){
      median_followup <- user_filtered_datasets_info[user_filtered_datasets_info$dataset_name==select_training_dataset,"median_followup"]
    }else{
      median_followup <- user_filtered_datasets_info[user_filtered_datasets_info$dataset_name==i,"median_followup"]
    }
    
    median_followup <- strsplit(median_followup,": ")[[1]][2]
    dacp[,m.b.endpoint] <- as.numeric(dacp[,m.b.endpoint])
    dacp[,b.endpoint.status] <- as.numeric(dacp[,b.endpoint.status])
    followup_value <- dacp[!is.na(dacp[,b.endpoint.status]) & !is.na(dacp[,m.b.endpoint]),m.b.endpoint]
    followup_range <- paste0(min(followup_value),"-",max(followup_value))
    da_followup <- c(i_type, i, length(followup_value), nrow(dacp),
                     round(length(followup_value)/nrow(dacp)*100,0),
                     b.endpoint, median_followup,followup_range)
    da_followups <- rbind(da_followups,da_followup)
  }
  
  da_followups <- as.data.frame(da_followups)
  colnames(da_followups) <- c("type","dataset_name","followup_num","all_num","percent","endpoint","median_followup","range_followup")
  da_followups$item13a <- paste0(da_followups$dataset_name,": ",
                                 da_followups$followup_num, " out of ",
                                 da_followups$all_num, " (",
                                 da_followups$percent, "%) had follow-up of ",
                                 da_followups$endpoint,
                                 " with a median follow-up (by reverse Kaplan-Meier method) of ",
                                 da_followups$median_followup, " months (range :",
                                 da_followups$range_followup, ").")
  
  da_followups$item13a_html <- paste0(da_followups$dataset_name, ": ", 
                                      da_followups$followup_num," out of ", 
                                      da_followups$all_num, " (",
                                      da_followups$percent, "%) had follow-up of ",
                                      da_followups$endpoint, 
                                      " with a median follow-up (by reverse Kaplan-Meier method) of <b>",
                                      da_followups$median_followup, "</b> months (range: <b>",
                                      da_followups$range_followup, "</b>).")
  
  # Item of results ----------------------------------------
  # item13 patient characteristic ----
  item13a_training_text_html <- da_followups[da_followups$type=="Training dataset","item13a_html"]
  item13a_validations_text_html <- paste0(da_followups[da_followups$type=="Validation dataset","item13a_html"], collapse = "<br /><br /> ")
  patient_characteristic_summary.df <- get_summary_test(train_valid_df_list = train_valid_df_list,
                                                        signature.type = b.signature.type,
                                                        used_validation_dataset = used_validation_dataset)
  
  patient_characteristic_summary.df$Variable[which(patient_characteristic_summary.df$Variable %in% "Age")] <- "Age (year)"
  unique_variables <- unique(patient_characteristic_summary.df$Variable)
  patient_characteristic_summary.df$Characteristic[which(patient_characteristic_summary.df$Characteristic %in% unique_variables)] <- NA
  
  patient_characteristic_summary.df <- patient_characteristic_summary.df[apply(patient_characteristic_summary.df[,!colnames(patient_characteristic_summary.df) %in% c("Variable","Raw Rank")], 1, function(x) !all(is.na(x))),]
  
  patient_characteristic_summary_data <- get_format_datatable(train_valid_df_list_names = train_valid_df_list_names,
                                                              signature.type = b.signature.type, 
                                                              df = patient_characteristic_summary.df,
                                                              analysis = "Participant characteristic",
                                                              train_name = select_training_dataset,
                                                              js_dir = js_dir)
  
  if(b.signature.type == "Prognostic") {
    if(length(used_validation_dataset) > 1) {
      patient_characteristic_legend <- "Key clinical and pathological characteristics of participants were included. P denotes the P value of the characteristics distribution comparison between the training dataset and each validation dataset. Categorical variables were compared using fisher'exact test, and continuous variables using t-test."
    } else if(length(used_validation_dataset) == 1) {
      patient_characteristic_legend <- "Key clinical and pathological characteristics of participants were included. P denotes the P value of the characteristics distribution comparison between the training dataset and validation dataset. Categorical variables were compared using fisher'exact test, and continuous variables using t-test."
    } else {
      patient_characteristic_legend <- "Key clinical and pathological characteristics of participants were included in training dataset."
    }
  } else {
    if(length(used_validation_dataset) > 1) {
      patient_characteristic_legend <- "Key clinical and pathological characteristics of participants were included. P in each dataset column denotes the P value of the characteristics distribution comparison between the control group and treatment group patients. P in \"Dataset Comparison\" column denotes the P value of the characteristics distribution comparison between the training dataset and each validation dataset. Categorical variables were compared using fisher'exact test, and continuous variables using t-test."
    } else if(length(used_validation_dataset) == 1) {
      patient_characteristic_legend <- "Key clinical and pathological characteristics of participants were included. P in each dataset column denotes the P value of the characteristics distribution comparison between the control group and treatment group patients. P in \"Dataset Comparison\" column denotes the P value of the characteristics distribution comparison between the training dataset and validation dataset. Categorical variables were compared using fisher'exact test, and continuous variables using t-test."
    } else {
      patient_characteristic_legend <- "Key clinical and pathological characteristics of participants were included in training dataset. P denotes the P value of the characteristics distribution comparison between the control group and treatment group patients. Categorical variables were compared using fisher'exact test, and continuous variables using t-test."
    }
  }
  
  # Item_13b-c ----
  if(b.signature.type == "Prognostic") {
    if(length(used_validation_dataset) > 1) {
      item13b_c_title <- "<p><span style='font-size:15px;background-color:#B1BB4E;color:white;display:block;padding:5px'><b>Item 13b/c. Comparison of participant characteristics between the training and validation datasets</b></span></p>"
      item13b_c_patient_characteristic_legend <- "Key clinical and pathological characteristics of participants were included. P denotes the P value of the characteristics distribution comparison between the training dataset and each validation dataset. Categorical variables were compared using fisher'exact test, and continuous variables using t-test."
    } else if(length(used_validation_dataset) == 1) {
      item13b_c_title <- "<p><span style='font-size:15px;background-color:#B1BB4E;color:white;display:block;padding:5px'><b>Item 13b/c. Comparison of participant characteristics between the training and validation dataset</b></span></p>"
      item13b_c_patient_characteristic_legend <- "Key clinical and pathological characteristics of participants were included. P denotes the P value of the characteristics distribution comparison between the training dataset and each validation dataset. Categorical variables were compared using fisher'exact test, and continuous variables using t-test."
    } else {
      item13b_c_title <- "<p><span style='font-size:15px;background-color:#B1BB4E;color:white;display:block;padding:5px'><b>Item 13b/c. Characteristics of participants</b></span></p>"
      item13b_c_patient_characteristic_legend <- "Key clinical and pathological characteristics of participants were included in training dataset."
    }
  } else {
    if(length(used_validation_dataset) > 1) {
      item13b_c_title <- "<p><span style='font-size:15px;background-color:#B1BB4E;color:white;display:block;padding:5px'><b>Item 13b/c. Comparison of participant characteristics between the training and validation datasets</b></span></p>"
      item13b_c_patient_characteristic_legend <- "Key clinical and pathological characteristics of participants were included. P in each dataset column denotes the P value of the characteristics distribution comparison between the control group and treatment group patients. P in \"Dataset Comparison\" column denotes the P value of the characteristics distribution comparison between the training dataset and each validation dataset. Categorical variables were compared using fisher'exact test, and continuous variables using t-test."
    } else if(length(used_validation_dataset) == 1) {
      item13b_c_title <- "<p><span style='font-size:15px;background-color:#B1BB4E;color:white;display:block;padding:5px'><b>Item 13b/c. Comparison of participant characteristics between the training and validation dataset</b></span></p>"
      item13b_c_patient_characteristic_legend <- "Key clinical and pathological characteristics of participants were included. P in each dataset column denotes the P value of the characteristics distribution comparison between the control group and treatment group patients. P in \"Dataset Comparison\" column denotes the P value of the characteristics distribution comparison between the training dataset and validation dataset. Categorical variables were compared using fisher'exact test, and continuous variables using t-test."
    } else {
      item13b_c_title <- "<p><span style='font-size:15px;background-color:#B1BB4E;color:white;display:block;padding:5px'><b>Item 13b/c. Characteristics of participants</b></span></p>"
      item13b_c_patient_characteristic_legend <- "Key clinical and pathological characteristics of participants were included in training dataset. P denotes the P value of the characteristics distribution comparison between the control group and treatment group patients. Categorical variables were compared using fisher'exact test, and continuous variables using t-test."
    }
  }
  
  if(b.have.profiles == "Yes") {
    figure_exp_distr_train_validsets = get_value_tripod.ccb(config_file = tripod.ini.file, key = 'figure_exp_distr_train_validsets', space = TRUE)[[1]]
    
    item13b_c_profile_title <- paste0(train.profiles_type," profile distribution") 
    if(length(used_validation_dataset) > 0){
      if(length(used_validation_dataset) > 1) {
        item13b_c_profile_distribution_text <- glue::glue("The density plots present the distributions of the log2-transformed {tolower(train.profiles_type)} level for the training dataset {select_training_dataset}, and the validation datasets {used_validation_datasets}.")
      } else{
        item13b_c_profile_distribution_text <- glue::glue("The density plots present the distributions of the log2-transformed {tolower(train.profiles_type)} level for the training dataset {select_training_dataset}, and the validation dataset {used_validation_datasets}.")
      }
    } else {
      item13b_c_profile_distribution_text <- glue::glue("The density plots present the distributions of the log2-transformed {tolower(train.profiles_type)} level for the training dataset {select_training_dataset}.")
    }
    
  } else {
    item13b_c_profile_title <- NULL
    figure_exp_distr_train_validsets <- NULL
    item13b_c_profile_distribution_text <- NULL
  }
  
  # Item_14 ----
  if(exists("c_preds_unicox")) {
    table_uni_cox <- get_rowgroup_table(df = c_preds_unicox, coef = FALSE,
                                        have.profiles = b.have.profiles,
                                        train.profiles_type = train.profiles_type)
    if(b.have.profiles == "Yes") {
      if(nrow(c_preds_unicox) < 100) {
        item14_unicox_legend <- paste0("Unadjusted association between each candidate predictor and ",b.endpoint,".")
      } else {
        item14_unicox_legend <- paste0("Unadjusted association between top 100 candidate predictors and ",b.endpoint,".")
      }
    } else {
      item14_unicox_legend <- paste0("Unadjusted association between each candidate clinicopathological variable and ",b.endpoint,".")
    }
  } else {
    table_uni_cox <- data.frame()
    item14_unicox_legend <- ""
  }
  
  # item15a
  item15a_signature_text <- paste0(length(unique(coef_report[,2])), 
                                   " variables comprising the developed signature and their corresponding coefficients were calculated by the Cox regression model.<br />The <b>signature score</b> of a patient can be calculated using the following formula: <br />",
                                   risk_score_formula, 
                                   ".<br />A higher signature score indicates a relative poorer prognosis.")
  
  colnames(coef_report) <- gsub("\n","<br />",colnames(coef_report))
  coef_report_new <- get_rowgroup_table(df = coef_report, have.profiles = b.have.profiles, train.profiles_type = train.profiles_type)
  
  # item15b 
  table_uni_multi_cox <- get_format_datatable(train_valid_df_list_names = train_valid_df_list_names,
                                              signature.type = b.signature.type,
                                              df = uni_multi_cox_report,
                                              analysis = "uni_multi_cox", js_dir = js_dir)
  
  # item16
  if(b.signature.type == "Prognostic") {
    if(is.null(b.clinical.variables) & is.null(b.interaction.test)) {
      item14_cox_legend <- "Univariate hazard ratio of each variable with a 95% confidence interval provided in parentheses, indicate multiplicative effects on the hazard. P denotes the statistical significances in the hazard ratio between a test group relative to a reference group, which were tested by the wald test."
    } else {
      item14_cox_legend <- "Univariate hazard ratio and multivariate hazard ratio of each variable with a 95% confidence interval provided in parentheses, indicate multiplicative effects on the hazard. P denotes the statistical significances in the hazard ratio between a test group relative to a reference group, which were tested by the wald test. For multivariate Cox analysis, the relationship between the variable of interest and the outcome was evaluated after adjusting for potential confounding variables that may be related to the outcome."
    }
  } else {
    if(is.null(b.clinical.variables) & is.null(b.interaction.test)) {
      item14_cox_legend <- "Univariate hazard ratio of each variable with a 95% confidence interval provided in parentheses, indicate multiplicative effects on the hazard. P denotes P values for survival difference of subgroup of clinicopathological variables in different therapies. P_interaction denotes P values for association between the variables of interest (including predictive group) and therapy benefit. The P values were tested by Wald test."
    } else {
      item14_cox_legend <- "Univariate hazard ratio and multivariate hazard ratio of each variable with a 95% confidence interval provided in parentheses, indicate multiplicative effects on the hazard. P denotes P values for survival difference of subgroup of clinicopathological variables in different therapies. P_interaction denotes P values for association between the variables of interest (including predictive group) and therapy benefit. The P values were tested by Wald test. For multivariate Cox analysis, the relationship between the variable of interest and the outcome was evaluated after adjusting for potential confounding variables that may be related to the outcome."
    }
  }
  
  # Item_15a ----
  item15a_signature_table_title <- paste0("Signature ",predictor_type," and the corresponding coefficients")
  item15a_signature_text <- paste0(length(unique(coef_report[,2])), 
                                   " variables comprising the developed signature and their corresponding coefficients were calculated by the Cox regression model.<br />The <b>signature score</b> of a patient can be calculated using the following formula: <br />",
                                   risk_score_formula,
                                   ".<br />A higher signature score indicates a relative poorer prognosis.")
  
  if(b.have.profiles == "Yes") {
    if(nrow(coef_report) > 1) {
      figure_sig_corrplot_exist <- TRUE
      figure_sig_corrplot = get_value_tripod.ccb(config_file = tripod.ini.file, key = 'figure_sig_corrplot', space = TRUE)[[1]]
      item15a_profile_correlation_title <- paste0("Pearson correlation plot for pairwise ", profile_type, " comparison among signature ", predictor_type)
      
      item15a_profile_correlation_text <- paste0("The correlation plot shows the Pearson correlation coefficients of ", 
                                                 profile_type, 
                                                 " profiles between each pair of signature predictors in the training dataset. From the plot, the correlation pattern of signature predictors can be visually checked if their expression levels are independent. The correlation coefficient ranges from -1 to 1. A smaller absolute value implies that a lower linear dependency, indicating this pair of signature predictors maybe is independent.")
      
      item15a_profile_correlation_title <- gsub(" ", " ", item15a_profile_correlation_title)
      item15a_profile_correlation_text <- gsub(" ", " ", item15a_profile_correlation_text) 
    } else {
      figure_sig_corrplot_exist <- FALSE
      item15a_profile_correlation_title <- ""
      figure_sig_corrplot <- NULL
      item15a_profile_correlation_text <- paste0("There has no Pearson correlation plot of ",
                                                 profile_type,
                                                 " profiles in the training dataset because of one signature predictor.")
      item15a_profile_correlation_text <- gsub(" "," ",item15a_profile_correlation_text)
    }
  } else {
    figure_sig_corrplot_exist <- FALSE
    item15a_profile_correlation_title <- ""
    figure_sig_corrplot <- NULL
    item15a_profile_correlation_text <- NULL
  }
  
  # Item_15b ----
  item15b_nomogram_title <- paste0("Nomogram for prediction of ",paste0(b.times,collapse=", "),"-months ",b.endpoint," probabilities")
  item15b_nomogram_text <- paste0(item10c_nomogram, 
                                  " In the nomogram, each factor was assigned a weighted score. To use this nomogram, first draw a line straight upward to the Point's axis to determine how many points toward the probability of each variable. Then sum the points achieved for each of the predictors. Last locate the final sum on the Total Points axis and draw a line straight down to ﬁnd the patient's probability of ",
                                  b.endpoint,".")
  
  if(b.signature.type == "Prognostic") { 
    if(ncol(train_multi_sig_nom)>=4) {
      figure_nomogram_exist = get_value_tripod.ccb(config_file = tripod.ini.file, key = 'figure_nomogram_exist', space = TRUE)[[1]]
      figure_nomogram_exist <- as.logical(figure_nomogram_exist)
      figure_nomogram = get_value_tripod.ccb(config_file = tripod.ini.file, key = 'figure_nomogram', space = TRUE)[[1]]
      
      figure_calibrate_lines = get_value_tripod.ccb(config_file = tripod.ini.file, key = 'figure_calibrate_lines', space = TRUE)[[1]]
      if(length(b.times) == 1) {
        item16_calibration_text <- paste("Calibration plot at ", 
                                         paste0(b.times,collapse=", "), 
                                         "-month was generated to explore the performance of the nomogram. Nomogram-predicted probabilities and observed outcome were plotted on the X-axis and Y-axis, respectively. A calibration plot along the 45-degree line indicates perfect consistency between the actual and nomogram-predicted prognosis. The vertical bars represent 95% CIs.")
      } else {
        item16_calibration_text <- paste("Calibration plots at ", 
                                         paste0(b.times,collapse=", "),
                                         "-months were generated to explore the performance of the nomogram. Nomogram-predicted probabilities and observed outcome were plotted on the X-axis and Y-axis, respectively. A calibration plot along the 45-degree line indicates perfect consistency between the actual and nomogram-predicted prognosis. The vertical bars represent 95% CIs.")
      }
      
      if(b.have.profiles == "Yes") {
        if(b.nomogram.option == "test") {
          item16_calibration_title <- "Calibration of the combined model integrating the signature and independent clinicopathological variables (training dataset)"
        } else {
          item16_calibration_title <- "Calibration of the combined model integrating the signature and key clinicopathological variables (training dataset)"
        }
      } else {
        item16_calibration_title <- "Calibration of the nomogram developed based on the signature (training dataset)"
      }
      
    } else {
      figure_nomogram_exist = FALSE
      figure_nomogram <- NULL
      figure_calibrate_lines <- NULL
      item15b_nomogram_text = "Nomogram was not constructed as no clinicopathological variable was selected in the development form."
      item16_calibration_text <- "Calibration curve was not generated as no nomogram was constructed."
      item16_calibration_title <- NULL
    }
    resave(figure_nomogram_exist, figure_nomogram,
           item15b_nomogram_title, item15b_nomogram_text,
           item16_calibration_title,
           figure_calibrate_lines, item16_calibration_text,
           file = add_info_tripod_rdata)
    
  } else {
    if(ncol(train_multi_sig_nom_c) >= 4 & ncol(train_multi_sig_nom_t) >= 4) {
      figure_nomogram_exist = get_value_tripod.ccb(config_file = tripod.ini.file, key = 'figure_nomogram_exist', space = TRUE)[[1]]
      figure_nomogram_exist <- as.logical(figure_nomogram_exist)
      
      if(length(b.times) == 1) {
        item16_calibration_text <- paste("Calibration plot at ", 
                                         paste0(b.times,collapse=", "), 
                                         "-month was generated to explore the performance of the nomogram. Nomogram-predicted probabilities and observed outcome were plotted on the X-axis and Y-axis, respectively. A calibration plot along the 45-degree line indicates perfect consistency between the actual and nomogram-predicted prognosis. The vertical bars represent 95% CIs.")
      } else {
        item16_calibration_text <- paste("Calibration plots at ",
                                         paste0(b.times,collapse=", "),
                                         "-months were generated to explore the performance of the nomogram. Nomogram-predicted probabilities and observed outcome were plotted on the X-axis and Y-axis, respectively. A calibration plot along the 45-degree line indicates perfect consistency between the actual and nomogram-predicted prognosis. The vertical bars represent 95% CIs.")
      }
      
      if(b.have.profiles == "Yes"){
        if(b.nomogram.option == "test"){
          item16_calibration_title <- "Calibration of the combined model integrating the signature and independent clinicopathological variables (training dataset)"
        } else {
          item16_calibration_title <- "Calibration of the combined model integrating the signature and key clinicopathological variables (training dataset)"
        }
      } else {
        item16_calibration_title <- "Calibration of the nomogram developed based on the signature (training dataset)"
      }
      
      figure_nomogram_all_pts = get_value_tripod.ccb(config_file = tripod.ini.file, key = 'figure_nomogram_all_pts', space = TRUE)[[1]]
      
      figure_calibrate_lines_all_pts = get_value_tripod.ccb(config_file = tripod.ini.file, key = 'figure_calibrate_lines_all_pts', space = TRUE)[[1]]
    } else {
      figure_nomogram_exist = FALSE
      figure_nomogram_all_pts <- NULL
      figure_calibrate_lines_all_pts <- NULL
      item15b_nomogram_text <- "Nomogram was not constructed as no clinicopathological variable was selected in the development form."
      item16_calibration_text <- "Calibration curve was not generated as no nomogram was constructed."
      item16_calibration_title <- NULL
    }
    
    resave(item15b_nomogram_title, item15b_nomogram_text,
           figure_nomogram_exist, figure_nomogram_all_pts,
           item16_calibration_text, item16_calibration_title, 
           figure_calibrate_lines_all_pts,
           file = add_info_tripod_rdata)
  }
  
  # Item_16 ----
  # Clinicopathological Association tab ----
  cp_variables <- unique(CP_Stat_report$Variable)
  missing_valid_dataset_text <- NULL
  
  if(!is.null(b.cp.variables)){
    train_valid_CP <- paste0("P<br />",train_valid_df_list_names)
    train_valid_CP <- train_valid_CP[train_valid_CP %in% colnames(CP_Stat_report)]
    train_valid_CP <- gsub("P<br />","",train_valid_CP)
    cp_variables_text <- ifelse(length(cp_variables) == 1,
                                paste0(cp_variables," variable."),
                                paste0(paste0(cp_variables,collapse = ", "),
                                       " variables."))
    
    missing_valid_dataset <- setdiff(c("Training dataset",used_validation_dataset),train_valid_CP)
    if(length(missing_valid_dataset) == 0) {
      missing_valid_dataset_text <- NULL
    } else if (length(missing_valid_dataset) == 1) {
      missing_valid_dataset_text <- paste0(missing_valid_dataset, " dataset was not shown due to lack of ",cp_variables_text)
    } else {
      missing_valid_dataset_text <- paste0(paste0(missing_valid_dataset,collapse=", "), " datasets were not shown due to lack of ",cp_variables_text)
    }
    
    ##
    if(b.group == "2 groups") {
      CP_Stat_report_new <- get_format_datatable(train_valid_df_list_names = train_valid_CP,
                                                 signature.type = b.signature.type,
                                                 df = CP_Stat_report,
                                                 analysis = "Association analysis",
                                                 b.group = "2 groups", js_dir = js_dir)
    } else {
      CP_Stat_report_new <- get_format_datatable(train_valid_df_list_names = train_valid_CP,
                                                 signature.type = b.signature.type,
                                                 df = CP_Stat_report,
                                                 analysis = "Association analysis",
                                                 b.group = "3 groups", js_dir = js_dir)
    }
    
    ##
    if(b.signature.type == "Predictive") {
      if(length(train_valid_CP) > 2) {
        association_text <- paste0("Association analysis was employed to investigate whether the signature is correlated with clinicopathological variables and other molecular features in both training and validation datasets. Samples with valid values (such as non-missing, non-unknown) were compared in the analysis. P denotes the P values comparing the associations between clinicopathological/molecular variables and the predictive groups in each dataset. Benefit P denotes the P values comparing the associations between clinicopathological/molecular variables and therapies in benefit group. No-benefit P denotes the P values comparing the associations between clinicopathological/molecular variables and therapies in no-benefit group. P-subgroup denotes the P values between therapies (control and treatment) and the predictive groups in each subgroup of clinicopathological/molecular variables. Fisher'exact test was used for association analysis.<br />", missing_valid_dataset_text)
      } else if(length(train_valid_CP) == 2){
        association_text <- paste0("Association analysis was employed to investigate whether the signature is correlated with clinicopathological variables and other molecular features in both training and validation dataset. Samples with valid values (such as non-missing, non-unknown) were compared in the analysis. P denotes the P values comparing the associations between clinicopathological/molecular variables and the predictive groups in each dataset. Benefit P denotes the P values comparing the associations between clinicopathological/molecular variables and therapies in benefit group. No-benefit P denotes the P values comparing the associations between clinicopathological/molecular variables and therapies in no-benefit group. P-subgroup denotes the P values between therapies (control and treatment) and the predictive groups in each subgroup of clinicopathological/molecular variables. Fisher'exact test was used for association analysis.<br />", missing_valid_dataset_text)
      } else {
        association_text <- paste0("Association analysis was employed to investigate whether the signature is correlated with clinicopathological variables and other molecular features in training dataset. Samples with valid values (such as non-missing, non-unknown) were compared in the analysis. P denotes the P values comparing the associations between clinicopathological/molecular variables and the predictive groups in training dataset. Benefit P denotes the P values comparing the associations between clinicopathological/molecular variables and therapies in benefit group. No-benefit P denotes the P values comparing the associations between clinicopathological/molecular variables and therapies in no-benefit group. P-subgroup denotes the P values between therapies (control and treatment) and the predictive groups in each subgroup of clinicopathological/molecular variables. Fisher'exact test was used for association analysis.<br />", missing_valid_dataset_text)
      }
    } else {
      if(length(train_valid_CP) > 2) {
        association_text <- paste0("Association analysis was employed to investigate whether the signature is correlated with clinicopathological variables and other molecular features in both training and validation datasets. Samples with valid values (such as non-missing, non-unknown) were compared in the analysis. P denotes the P values comparing the associations between the clinicopathological/molecular variables and the risk groups in each dataset. Fisher'exact test was used for association analysis.<br />", missing_valid_dataset_text)
      } else if(length(train_valid_CP) == 2){
        association_text <- paste0("Association analysis was employed to investigate whether the signature is correlated with clinicopathological variables and other molecular features in both training and validation dataset. Samples with valid values (such as non-missing, non-unknown) were compared in the analysis. P denotes the P values comparing the associations between the clinicopathological/molecular variables and the risk groups in each dataset. Fisher'exact test was used for association analysis.<br />", missing_valid_dataset_text)
      } else {
        association_text <- paste0("Association analysis was employed to investigate whether the signature is correlated with clinicopathological variables and other molecular features in training dataset. Samples with valid values (such as non-missing, non-unknown) were compared in the analysis. P denotes the P values comparing the associations between the clinicopathological/molecular variables and the risk groups in training dataset. Fisher'exact test was used for association analysis.<br />", missing_valid_dataset_text)
      }
    }
    item16_association_text <- association_text
  } else {
    association_text <- "<p style=\"margin-left:-2%;font-size:17px\">Association analysis was not conducted as no clinicopathological variable was selected for analysis.</p>"
    item16_association_text <- ""
  }
  
  ifelse(is.null(b.cp.variables), CP_association_exists <- FALSE, CP_association_exists <- TRUE)
  
  # km training
  figure_km_trainset = get_value_tripod.ccb(config_file = tripod.ini.file, key = 'figure_km_trainset', space = TRUE)[[1]]
  # km validation
  figure_km_validsets_exist = get_value_tripod.ccb(config_file = tripod.ini.file, key = 'figure_km_validsets_exist', space = TRUE)[[1]]
  if(is.na(figure_km_validsets_exist)){
    figure_km_validsets_exist <- FALSE
  }
  figure_km_validsets_exist <- as.logical(figure_km_validsets_exist)
  if(figure_km_validsets_exist){
    figure_km_validsets = get_value_tripod.ccb(config_file = tripod.ini.file, key = 'figure_km_validsets', space = TRUE)[[1]]
  }else{ 
    figure_km_validsets <- NULL
  }
  
  #
  if(b.signature.type == "Predictive"){
    risk_groups <- "treatment-benefit groups"
  } else {
    risk_groups <- "risk groups"
  }
  if(length(used_validation_dataset) == 0) {
    item16_KM_text <- paste0("Kaplan-Meier survival curves for ", 
                             b.endpoint,
                             " in training dataset, stratified by the prediction model (high-risk and low risk). The predictive performance of the prognostic signature was evaluated by the two‐sided log‐rank test in training dataset. P < 0.05 was considered as statistically significant, and 95% confidence intervals are presented in brackets. HR, hazard ratio.")
    item16_surv_title <- paste0(paste0(b.times,collapse="/")," months ",b.endpoint," for different signature-predicted ",risk_groups)
    item16_surv_text <- paste0(paste0(b.times,collapse="/")," months ", b.endpoint,
                               " of patients in different risk groups for training dataset were predicted by the prediction model.")
    
  } else if(length(used_validation_dataset) == 1){
    item16_KM_text <- paste0("Kaplan-Meier survival curves for ", 
                             b.endpoint,
                             " in training and validation dataset, stratified by the prediction model (high-risk and low risk). The performance of the prognostic signature was evaluated by the two‐sided log‐rank test in both training dataset and validation dataset. P < 0.05 was considered as statistically significant, and 95% confidence intervals are presented in brackets. HR, hazard ratio.")
    
    item16_surv_title <- paste0(paste0(b.times,collapse="/")," months ", b.endpoint," for different signature-predicted ",risk_groups)
    item16_surv_text <- paste0(paste0(b.times,collapse="/")," months ", b.endpoint,
                               " of patients in different risk groups for each dataset were predicted by the prediction model.")
  } else {
    item16_KM_text <- paste0("Kaplan-Meier survival curves for ",
                             b.endpoint," in training and validation datasets, stratified by the prediction model (high-risk and low risk). The performance of the prognostic signature was evaluated by the two‐sided log‐rank test in both training dataset and validation datasets. P < 0.05 was considered as statistically significant, and 95% confidence intervals are presented in brackets. HR, hazard ratio.")
    item16_surv_title <- paste0(paste0(b.times,collapse="/")," months ",b.endpoint," for different signature-predicted ",risk_groups)
    item16_surv_text <- paste0(paste0(b.times,collapse="/")," months ",b.endpoint,
                               " of patients in different risk groups for each dataset were predicted by the prediction model.")
  }
  
  
  surv_prob_report_new <- get_format_datatable(train_valid_df_list_names = train_valid_df_list_names,
                                               signature.type = b.signature.type,
                                               df = surv_prob_report,
                                               analysis = "Kaplan_Meier_analysis",
                                               b.times=b.times, b.endpoint=b.endpoint,
                                               js_dir = js_dir)
  
  # roc
  figure_roc = get_value_tripod.ccb(config_file = tripod.ini.file, key = 'figure_roc', space = TRUE)[[1]]
  item16_ROC_text <- paste0("Time-dependent receiver operating characteristic (ROC) curves show the sensitivity and specificity of different variables in prognosis prediction. It could help to quantify and compare the discrimination ability of the signature and clinicopathological prognosis factors using the ", paste0(b.times,collapse=", "), "-months time-dependent ROC analysis in the training dataset. The area under the curve (AUC) ranges from 0.5 (no discrimination) to a theoretical maximum of 1, which were texted in the legend for each prognosis model. A prognosis model with larger AUC indicates a better performance.")
  
  # pe
  figure_pe = get_value_tripod.ccb(config_file = tripod.ini.file, key = 'figure_pe', space = TRUE)[[1]]
  item16_PE_text <- paste0("Prediction error (PE) curves were used to evaluate the performance of prediction models in the training dataset, and ten-fold cross-validation cumulative prediction error were computed using Kaplan-Meier estimation as reference. Integrated Brier score (IBScore) was defined as the area under a prediction error curve, which were texted in the legend for each prognosis model. A prognosis model with smaller IBScore indicates a better performance.<br />", PE_na_sams_reason)
  
  # item16 molecular signature
  if(b.have.profiles == "Yes") {
    figure_calibrate_lines_signature = get_value_tripod.ccb(config_file = tripod.ini.file, key = 'figure_calibrate_lines_signature', space = TRUE)[[1]]
    
    if(file.exists(figure_calibrate_lines_signature)) {
      if(length(b.times) == 1) {
        item16_calibration_molecular_legend <- paste0("Calibration plot at ", paste0(b.times,collapse=", "), 
                                                      "-month was generated to explore the performance of the signature. Signature-predicted probabilities and observed outcome were plotted on the X-axis and Y-axis, respectively. A calibration plot along the 45-degree line indicates perfect consistency between the actual and signature-predicted prognosis. The vertical bars represent 95% CIs.")
      } else {
        item16_calibration_molecular_legend <- paste0("Calibration plots at ", paste0(b.times,collapse=", "),
                                                      "-months were generated to explore the performance of the signature. Signature-predicted probabilities and observed outcome were plotted on the X-axis and Y-axis, respectively. A calibration plot along the 45-degree line indicates perfect consistency between the actual and signature-predicted prognosis. The vertical bars represent 95% CIs.")
      }
      figure_calibrate_lines_signature_exist <- TRUE
      
    } else {
      figure_calibrate_lines_signature_exist <- FALSE
      figure_calibrate_lines_signature <- NULL
      item16_calibration_molecular_legend <- "There has no calibration plot for molecular signature."
      
    }
    resave(figure_calibrate_lines_signature_exist,
           figure_calibrate_lines_signature,
           item16_calibration_molecular_legend,
           file = add_info_tripod_rdata)
  }
  
  # item13 diagram
  if(!is.null(user_upload_diagram_path)) {
    if(grepl("\\.png$", user_upload_diagram_path)) {
      figure_flow_diagram_exists <- TRUE
    } else if (grepl("\\.jpg$", user_upload_diagram_path)) {
      figure_flow_diagram_exists <- TRUE
    } else {
      figure_flow_diagram_exists <- FALSE
      user_upload_diagram_path <- NULL
    }
  } else {
    figure_flow_diagram_exists <- FALSE
  }
  
  resave(item4a_training_source, item4a_validation_source,
         item4a_training_study_date, item4a_validation_study_date,
         item5a_training_setting, item5a_validation_setting,
         item5b_training_exclusion_user, item5b_training_exclusion_web,
         item5b_training_inclusion_user, item5b_training_inclusion_web,
         item5b_validation_exclusion_user, item5b_validation_exclusion_web,
         item5b_validation_inclusion_user, item5b_validation_inclusion_web,
         item5b_validation_similar,
         item5c_training, item5c_validation,
         item6a_training_outcome, item6a_training_methods, item6b_training,
         item7a_training, item7b_training, item8_training,
         item9_training_methods, item9_validation_methods,
         item9_training_methods_user, item9_validation_methods_user,
         item10a_handling_molecular, item10a_implausible, item10a_handling_clinical,
         item10b_before_genes, item10b_before, item10b_model,
         item10c_signature, item10c_nomogram, item10c_online,
         item10d_independence, item10d_interaction, item10d_km,
         item10d_roc, item10d_pe, item10d_calibration,
         item10d_km_pro, item10d_association,
         item10e_updating, item10e_updating, item11_groups,
         item12_setting, item12_eligibility, item12_outcome, item12_predictors,
         item13a_training_text_html, item13a_validations_text_html,
         user_upload_diagram_path, figure_flow_diagram_exists,
         file = add_info_tripod_rdata)
  
  # local_model_summary, local_model_text
  training_cp_df <- train.cp.raw
  model_summary <- NULL
  model_text <- NULL
  if(b.have.profiles == "No") {
    variable_id <- which(colnames(coef_report) %in% c("Type","Variable"))
    uniq_variables <- unique(coef_report[,variable_id])
    
    for(i in uniq_variables) {
      sub_coef1 <- coef_report[coef_report[,variable_id] == i,"Characteristic"]
      if(length(sub_coef1) >1) {
        sub_coef2 <- gsub(paste0(i," \\(vs\\. "),"",sub_coef1)
        sub_coef3 <- gsub("\\)","",sub_coef2[1])
        sub_coefs <- sort(c(sub_coef3,sub_coef2[2:length(sub_coef2)]))
        sub_text <- paste0("<p>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;",i," (categorical): ",paste0(sub_coefs,collapse = ", "),"</p>")
        model_text <- c(model_text,sub_text)
        sub_summary <- data.frame("Variable" = i, "Summary" = paste0(sub_coefs,collapse = ", "), "Class" = "categorical")
        model_summary <- rbind.data.frame(model_summary,sub_summary)
      } else {
        if(i == "Age") {
          min_age <- min(training_cp_df[,i])
          sub_coefs <- paste0(i,": ",i," >= ",min_age)
          sub_text <- paste0("<p>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;",i," (continuous): >= ",min_age,"</p>")
          model_text <- c(model_text,sub_text)
          sub_summary <- data.frame("Variable" = i, "Summary" = paste0(">= ",min_age), "Class" = "continuous")
          model_summary <- rbind.data.frame(model_summary,sub_summary)
        } else {
          min_i <- min(training_cp_df[,i])
          max_i <- max(training_cp_df[,i])
          sub_coefs <- paste0(i,": ",i,"(",min_i,",",max_i,")")
          sub_text <- paste0("<p>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;",i," (continuous): ",min_i," - ",max_i,"</p>")
          model_text <- c(model_text,sub_text)
          sub_summary <- data.frame("Variable" = i, "Summary" = paste0(min_i," - ",max_i), "Class" = "continuous")
          model_summary <- rbind.data.frame(model_summary,sub_summary)
        }
      }
    }
    colnames(model_summary) <- c("Variable","Summary","Class")
    rownames(model_summary) <- model_summary$Variable
    model_text <- paste0(model_text,collapse = "")
  }
  local_model_summary <- model_summary
  local_model_text <- model_text
  
  ##
  resave(da_followups, surv_prob_report, CP_association_exists, CP_Stat_report,
         surv_prob_report, patient_characteristic_summary.df,
         patient_characteristic_summary_data, CP_Stat_report_new,
         table_uni_multi_cox, coef_report_new,
         surv_prob_report_new, table_uni_cox,
         item14_unicox_legend, item13b_c_title,
         item13b_c_patient_characteristic_legend,
         figure_exp_distr_train_validsets, item13b_c_profile_title,
         item13b_c_profile_distribution_text,
         item14_cox_legend, item15a_signature_table_title, item15a_signature_text,
         figure_sig_corrplot_exist, figure_sig_corrplot,
         item15a_profile_correlation_title,item15a_profile_correlation_text,
         figure_km_trainset, figure_km_validsets_exist, figure_km_validsets,
         item16_KM_text, item16_surv_title, item16_surv_text,
         figure_roc, item16_ROC_text, figure_pe,
         item16_PE_text, item16_association_text, 
         patient_characteristic_legend,
         geneset_name, b.profiling,
         b.disease.type, b.sample.type, b.primary.site,
         association_text, local_min_age, local_max_age,
         local_molecular_profiling,
         local_model_summary, local_model_text,
         b.prognostic.treatment, b.prognostic.treatment.setting, b.prognostic.regimen,
         
         file = add_info_tripod_rdata)
  
  
  section <- "Signature_reporting"
  key <- "add_info_tripod_rdata"
  value <- add_info_tripod_rdata
  system(paste("python3", paste0(script.dir, "/", "write_ini.py"), user.config.ini.file, section, key, value))
  
  # my signatures table -------------------------------------------------------
  colnames(train_df) <- tolower(colnames(train_df))
  
  # is primary treatment
  m_is_primary_treatment <- paste0(unique(train_df$is_primary_treatment),collapse=", ")
  
  # treatment type ---
  treatment_type_id <- which(colnames(train_df) %in% c("Treatment type","treatment_type","Treatment_type","treatment type"))
  if(length(treatment_type_id) == 1){
    m_treatment_type <- paste0(unique(train_df[,treatment_type_id][which(!is.na(train_df[,treatment_type_id]))]),collapse=", ")
  }else{
    m_treatment_type <- NA
  }
  
  # treatment setting ---
  treatment_setting_id <- which(colnames(train_df) %in% c("Treatment setting","treatment_setting","Treatment_setting","treatment setting"))
  if(length(treatment_setting_id) == 1){
    m_treatment_setting <- paste0(unique(train_df[,treatment_setting_id][which(!is.na(train_df[,treatment_setting_id]))]),collapse=", ")
  }else{
    m_treatment_setting <- NA
  }
  
  # regimen ---
  regimen_id <- which(colnames(train_df) %in% c("Regimen","regimen"))
  if(length(regimen_id) == 1){
    m_regimen <- paste0(unique(train_df[,regimen_id][which(!is.na(train_df[,regimen_id]))]),collapse=",, ")
  }else{
    m_regimen <- NA
  }
  
  # gender
  m_gender <- paste0(unique(train_df$gender),collapse=", ")
  
  # ethnicity
  m_ethnicity<- paste0(unique(train_df$ethnicity),collapse=", ")
  
  # stage edition ----
  stage_edition_id <- which(colnames(train_df) %in% c("Stage_edition","stage_edition"))
  if(length(stage_edition_id) == 1){
    m_stage_edition <- sort(unique(train_df[,stage_edition_id][which(!is.na(train_df[,stage_edition_id]))]))
    m_stage_edition <- paste0(m_stage_edition,collapse=", ")
  }else{
    m_stage_edition <- NA
  }
  
  # stage ----
  stage_id <- which(colnames(train_df) %in% c("Stage","stage"))
  if(length(stage_id) == 1){
    m_stage <- sort(unique(train_df[,stage_id][which(!is.na(train_df[,stage_id]))]))
    m_stage <- paste0(m_stage,collapse=", ")
  }else{
    m_stage <- NA
  }
  
  ##
  signature_id <- b.name.sig.id
  m_development_date = as.character(Sys.time())
  m_signature_name = b.name.sig
  if(b.have.profiles == "No") {
    reference_id <- grep("\\(vs.",coef_report$Characteristic)
    coef_report <- coef_report[-reference_id,]
    m_signature_gene = paste0(coef_report$Characteristic,collapse=", ")
    coef_report$Coefficient <- round(as.numeric(coef_report$Coefficient),3)
    m_coefs = paste0(coef_report$Coefficient,collapse=",")
  }
  m_signature_gene = paste0(coef_report$Characteristic,collapse=", ")
  coef_report$Coefficient <- round(as.numeric(coef_report$Coefficient),3)
  m_coefs = paste0(coef_report$Coefficient,collapse=",")
  m_primary_site = b.primary.site
  m_disease_type <- paste0(unique(train_df$disease_type),collapse=", ")
  m_sample_type = b.sample.type
  
  if(b.have.profiles == "Yes"){
    m_platform_type = paste0(unique(train_df$molecular_profiling),collapse=", ")
  }else{
    m_platform_type = "No profiling"
  }
  
  m_signature_type = paste0(b.signature.type," signature")
  m_endpoint = b.endpoint
  
  m_training = select_training_dataset
  m_validation = paste0(used_validation_dataset,collapse="<br />")
  
  m_signature_path = b.user.sig.path
  
  m_signature_group = "User-developed signature"
  
  m_source = "Private signature"
  
  local_sig_df = data.frame("signature_id"=signature_id, "signature_name"=m_signature_name, "signature_description" = signature_description,
                            "development_date"=m_development_date, "signature_gene"=m_signature_gene,
                            "coefficients"=m_coefs,"primary_site"=m_primary_site,
                            "disease_type"=m_disease_type,"sample_type"=m_sample_type,
                            "molecular_profiling"=m_platform_type,"signature_type"=m_signature_type,
                            "is_primary_treatment"=m_is_primary_treatment,
                            "treatment_type"=m_treatment_type,"treatment_setting"=m_treatment_setting,
                            "regimen"=m_regimen,"gender"=m_gender,"ethnicity"=m_ethnicity,
                            "endpoint"=m_endpoint,"stage_edition"=m_stage_edition,"stage"=m_stage,
                            "training_dataset"=m_training,"validation_datasets"=m_validation,
                            "signature_path"=m_signature_path,"signature_group"=m_signature_group,
                            "source"=m_source)
  
  resave(local_sig_df, file = add_info_tripod_rdata)
  
  # download tripod -------------------------------------------------------
  dir.create(paste0(b.user.sig.path, "/upload/"))
  file <- paste0(b.user.sig.path, "/upload/", b.name.sig.id, ".", Sys.Date(), ".report.html")
  src <- ifelse(b.have.profiles == "Yes",
                normalizePath(paste0(script.dir,"/rmarkdown_exp.Rmd")),
                normalizePath(paste0(script.dir,"/rmarkdown_noexp.Rmd")))
  
  
  rmd.file <- ifelse(b.have.profiles == "Yes",
                     paste0(repromsig.dir, "/config/rmarkdown_exp.Rmd"),
                     paste0(repromsig.dir, "/config/rmarkdown_noexp.Rmd"))
  
  
  owd <- setwd(tempdir())
  on.exit(setwd(owd))
  file.copy(src, rmd.file, overwrite = TRUE)
  library(rmarkdown)
  out <- render(rmd.file, html_document())
  
  file.rename(out, file)
  
  if(file.exists(file)) {
    print(paste0("Output created: ", file))
  }
} else {
  message_nosignature <- "There was no signature detected!"
  print(message_nosignature)
}