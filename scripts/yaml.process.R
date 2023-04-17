# /home/pub/tools/R-3.6.2/bin/Rscript /opt/shiny-server/apps/repromsig/scripts/yaml.process.R /opt/shiny-server/apps/repromsig/ColoGuide_Stage_II_local/input/analysis.yaml


# HEADER ------------------------------------------------------------------ 
rm(list=ls())
options(stringsAsFactors = FALSE)
DEBUG_MODE = FALSE
ARGS_MODE = TRUE

# Input/Output/Lib/Config/Params ------------------------------------------
# 1) Parameters
if(ARGS_MODE) {
	args<-commandArgs(TRUE)
	if(length(args)!=1) {
		stop("Usage: Rscript yaml.process.R [analysis.yaml]\n")
	}
	analysis.yaml.file <- args[1]
} else {
    analysis.yaml.file <- "/opt/shiny-server/apps/repromsig/ColoGuide_Stage_II_local/input/analysis.yaml"
}

# 2) Library
library(dplyr)
library(yaml)
library(Hmisc)
library(survival)
library(survminer)

# readin and preprocess  ---------------------------
conf <- read_yaml(analysis.yaml.file)

# source functions ---
repromsig.dir <- conf$repromsig_dir
script.dir <- paste0(repromsig.dir, "/scripts/")
source(paste0(script.dir,'ccb.helper.R'))


## indir, outdir and paths for public datasets ---
in.dir <- conf$input_dir
b.user.sig.path <- conf$output_dir
mysetwd(b.user.sig.path)
mysetwd(in.dir)

## generate user_filtered_datasets_info (a data frame) of User_filtered_datasets.RData ------------------------------------------------------
b.upload.predictors <- read.file(conf$candidate_genes_file)

user_filtered_datasets_info <- data.frame()


# Check training and validation dataset(s) ------------------------------------------------------
if(is.null(conf$training_cohort$dataset_name)) {
    stop("The dataset name of training cohort is missing in the yaml file.")
}
ids <- names(conf)[grep("validation_cohort", names(conf))] %>% c("training_cohort", .)

# datasets selected
select_ds <- do.call(c, lapply(ids, function(x) { conf[[x]]$dataset_name }))
select_training_dataset <- select_ds[1]

## 
if(is.null(conf$training_cohort$patient_annotation_file)) {
    stop("The patient_annotation_file for training cohort is missing in yaml file.")
}

if(!is.null(conf$training_cohort$molecular_profile_file)) {
    molecluar_profiles_included = b.have.profiles = "Yes"
} else {
    molecluar_profiles_included = b.have.profiles = NULL
}

if(!is.null(conf$training_cohort$molecular_profile_file) & is.null(conf$training_cohort$molecular_profiling)) {
    stop("The molecular_profiling (Molecular Assay Platform) for training cohort is missing in yaml file.")
}


## for private datasets, generate the user_uploaded_list of User_uploaded_anno.RData
# and save the expr matrix to rdata file ------------------------------------------------------
is_primary_treatment <- "Yes"
source <- "Private dataset"
user_uploaded_list <- list()
user_filtered_datasets_info <- data.frame()
data_ids <- c()
select_ds <- c()
z = 1
sr <- "Private dataset"

for(x in ids) {
    list0 <- conf[[x]]
    name0 <- list0$dataset_name
    ##
    primary_site <- capitalize(tolower(list0$primary_site))
    profile_platform <- list0$molecular_profiling
    profile_platform <- ifelse(is.null(profile_platform), NA, profile_platform)
    meta_infor1 <- data.frame(name0, primary_site, list0$sample_type, 
                              profile_platform, list0$log_transform_type,
                              is_primary_treatment, source
                             )         
    user_filtered_datasets_info <- rbind(user_filtered_datasets_info, meta_infor1)


    # Only keep datasets with patient_annotation_file
    if((sr == "Private dataset" & !is.null(list0$patient_annotation_file))) {
        if(!is.null(list0$molecular_profile_file)) {
            if(file.exists(list0$molecular_profile_file)) {
                expr <- read.file(list0$molecular_profile_file)
            } else {
                expr <- data.frame()
                Warning(paste0("The molecular_profile_file: ", list0$molecular_profile_file, " for training cohort was not found."))
            }

            dir.create(paste0(b.user.sig.path, "/rda/"))
            save(expr, file = paste0(b.user.sig.path, "/rda/", name0, "_exp.RData"))
        }

        ## cp table 
        pts_ann <- read.file(list0$patient_annotation_file)


        user_uploaded_list[[z]] <- pts_ann
        data_ids <- c(data_ids, name0)
        z=z+1
        ##
        if(x == "training_cohort") {
            b.primary.site <- capitalize(tolower(list0$primary_site))
            b.sample.type <- capitalize(tolower(list0$sample_type))
        }
        select_ds <- c(select_ds, name0)
    }
}

names(user_uploaded_list) <- data_ids
select_training_dataset <- select_ds[1]
select_validation_dataset <- select_ds[-1]

colnames(user_filtered_datasets_info) <- c("dataset_name", "primary_site", "sample_type",
                                           "molecular_profiling", "log_transform_type",
                                           "is_primary_treatment", "source")



# preprocessing the patient_annotation_file for each valid cohort (user_uploaded_list) ------
user_uploaded_list_2 <- list()
meta_infor <- data.frame()
z=1
for(d_name in names(user_uploaded_list)) {
    print(paste0(d_name," dataset starts processing ..."))
    cp_df <- user_uploaded_list[[d_name]]

    ##
    u_disease_type <- unique(cp_df$Disease_type)
    if(is.na(u_disease_type)) {
        stop(paste0("Please provide the disease type in ", d_name, "'s patient annotation file."))
    } 

    u_primary_site <- unique(cp_df$Primary_site)
    if(is.na(u_primary_site)) {
        stop(paste0("Please provide the primary site in ", d_name, "'s patient annotation file."))
    }

    ## 
    name0 <- gsub("\\.|-", "_", colnames(cp_df))
    name1 <- tolower(name0)
    if(nrow(cp_df) == 0) { stop("The provided patient annotation file is empty.") }
    colnames(cp_df) <- tolower(gsub("\\.|-", "_", colnames(cp_df)))
    cp_df$dataset_id <- d_name

    # Is primary treatment processs ----
    id0 <- which(colnames(cp_df) %in% c("is_primary_treatment", "is primary treatment", "is.primary.treatment"))
    if(length(id0) != 1) { stop("The \"Is_primary_treatment\" column in the patient annotation file was not found.") }
    colnames(cp_df)[id0[1]] <- 'is_primary_treatment'
    ##
    cp_df[is.na(cp_df$is_primary_treatment), 'is_primary_treatment'] <- "NA"
    cp_df[cp_df$is_primary_treatment %in% c("no","No"), 'is_primary_treatment'] <- "No"
    cp_df[cp_df$is_primary_treatment %in% c("yes","Yes"), 'is_primary_treatment'] <- "Yes"
    tx <- cp_df %>% filter(!is_primary_treatment %in% c("No","Yes")) %>% select(is_primary_treatment) %>% unlist() %>% unique()
    if(length(tx) > 0) { stop("Please check the annotation file as the optional value is \"Yes\" or \"No\" for \"Is_primary_treatment\" column.") }
    ## 
    tx <- as.data.frame(table(cp_df$is_primary_treatment))
    if(nrow(tx) > 0) {
        tx$com <- paste0(tx[,1]," (",tx[,2], ")")
        is_primary_treatment <- paste0(tx$com,collapse = "<br />")
    }

    # Patient ID process ----
    pts_id <- which(colnames(cp_df) %in% c("patient_id","patient id","patient.id","patientid"))
    sam_id <- which(colnames(cp_df) %in% c("sample_id","sample id","sample.id","sampleid"))

    if(length(pts_id) >= 1) {
        colnames(cp_df)[pts_id[1]] <- 'patient_id'
        anno_p_id <- cp_df$patient_id
        na_id <- cp_df %>% filter(is.na(patient_id) | patient_id == "") %>% nrow()
        dup_id <- anno_p_id[duplicated(anno_p_id)] 
        
        if(na_id > 0) {
            stop("Empty patient ID(s) were found in the uploaded annotation file.")
        } else {
            if(length(dup_id) > 0){
                dup_id_p <- paste0(dup_id, collapse = ", ")
                dup_id_warning <- paste0("The same patient ID (", dup_id_p, ") exist(s) in the annotation file.")
                stop(dup_id_warning)
            } else {
                if(length(sam_id) == 1) {
                    colnames(cp_df)[sam_id[1]] <- 'sample_id'
                    anno_s_id <- cp_df$sample_id
                    na_samid <- cp_df %>% filter(is.na(sample_id) | sample_id == "") %>% nrow()
                    dup_s_id <- anno_s_id[duplicated(anno_s_id)] 
                    if(na_samid > 0) { stop("Empty sample ID(s) were found in the annotation file.") }

                    if(length(dup_s_id) > 0){
                        dup_s_id_p <- paste0(dup_s_id, collapse = ", ")
                        dup_s_id_warning <- paste0("Please check the annotation file as the same sample ID (", 
                                                                dup_s_id_p,") exist(s) in the file.")
                        stop(dup_s_id_warning)
                    }
                } else {
                    cp_df$sample_id <- cp_df$patient_id
                    anno_s_id <- cp_df$sample_id
                }
            }
        }
    } else {
        stop("The \"Patient_ID\" column was not found in the patient annotation file.")
    }

    cp_df[cp_df == "NA" | cp_df == "na" | cp_df == "<NA>" | cp_df == "#NULL!" | cp_df == "?" | cp_df == "" | is.na(cp_df)] <- NA

    # Profiling process
    if(all(is.na(cp_df$molecular_profiling))) {
        cp_df$molecular_profiling <- "NA"
        cp_df$platform <- "NA"
    }

    # Sample type process ----
    ids <- which(colnames(cp_df) %in% c("sample_type", "sample type", "sample.type"))
    if(length(ids) >= 1) {
        colnames(cp_df)[ids[1]] <- 'sample_type'
        cp_df$sample_type <- tolower(cp_df$sample_type)

        cp_df[!is.na(cp_df$sample_type) & cp_df$sample_type %in% c("primary tumor"), 'sample_type'] <- 'Primary tumor'
        cp_df[!is.na(cp_df$sample_type) & cp_df$sample_type %in% c("metastatic tumor"), 'sample_type'] <- 'Metastatic tumor'
        cp_df[!is.na(cp_df$sample_type) & cp_df$sample_type %in% c("recurrent tumor"), 'sample_type'] <- 'Recurrent tumor'
        cp_df[is.na(cp_df$sample_type) | cp_df$sample_type %in% c("na", "unknown","Unknown","[Discrepancy]","Discrepancy"), 'sample_type'] <- "NA"
        cp_df[!cp_df$sample_type %in% c('NA','Primary tumor','Metastatic tumor','Recurrent tumor'), 'sample_type'] <- "Other"
    } else {
        cp_df$sample_type <- "NA"
    }
    ##
    tx <- as.data.frame(table(cp_df$sample_type))
    if(nrow(tx) > 0) {
        tx$com <- paste0(tx[,1], " (", tx[,2], ")")
        sample_type <- paste0(tx$com, collapse = "<br />")
    }

    # treatment_type process ----
    ids <- which(colnames(cp_df) %in% c("treatment","treatment_type","treatment.type","treatment type"))
    if(length(ids) >= 1) {
        colnames(cp_df)[ids[1]] <- 'treatment_type'
        cp_df[is.na(cp_df$treatment_type), 'treatment_type'] <- "NA"
    } else {
        cp_df$treatment_type <- "NA"
    }

    tx <- as.data.frame(table(cp_df[,'treatment_type']))
    if(nrow(tx) > 0) {
        tx$com <- paste0(tx[,1], " (", tx[,2], ")")
        treatment_type <- paste0(tx$com, collapse = "<br />")
    }

    # treatment_setting process ----
    ids <- which(colnames(cp_df) %in% c("treatment_setting","treatment.setting","treatment setting"))
    if(length(ids) >= 1) {
        colnames(cp_df)[ids[1]] <- 'treatment_setting'
        cp_df[is.na(cp_df$treatment_setting), 'treatment_setting'] <- "NA"
        cp_df[cp_df$treatment_type=="Surgery", 'treatment_setting'] <- "Not applicable"
    } else {
        cp_df$treatment_setting <- "NA"
    }

    tx <- as.data.frame(table(cp_df$treatment_setting))
    if(nrow(tx) > 0) {
        tx$com <- paste0(tx[,1], " (", tx[,2], ")")
        treatment_setting <- paste0(tx$com, collapse = "<br />")
    }

    # regimen ----
    ids <- which(colnames(cp_df) %in% c("regimen", "regimens"))
    if(length(ids) >= 1) {
        colnames(cp_df)[ids[1]] <- 'regimen'
        cp_df[is.na(cp_df$regimen), 'regimen'] <- "NA"
        cp_df[cp_df$treatment_type == "Surgery", 'regimen'] <- "Not applicable"
    } else {
        cp_df$regimen <- "NA"
    }

    tx <- as.data.frame(table(cp_df$regimen))
    if(nrow(tx) > 0) {
        tx$com <- paste0(tx[,1], " (", tx[,2], ")")
        regimen <- paste0(tx$com, collapse = "<br />")
    } else {
        regimen <- "NA"
    }

    # Age, gender, race, t, n, m, stage process ----
    cp_df <- age_gender_race_process.func(df = cp_df, column_ids = c("age", "gender", "ethnicity"))
    cp_df <- t_n_m_process.func(df = cp_df, column_ids = c("t", "n", "m"))
    cp_df <- stage_process.func(df = cp_df, column_ids = c("stage_edition", "stage"))
    ##
    age <- age_summary(age_value = cp_df$age)
    gender <- all_sub(colname="gender", subclass=c('Female', 'Male', 'NA'), df = cp_df)
    tt <- all_sub(colname = "t", subclass = c('T1','T2','T3','T4','NA'), df = cp_df)
    nn <- all_sub(colname = "n", subclass = c('N0','N1','N2','N3','NA'), df = cp_df)
    mm <- all_sub(colname = "m", subclass=c('M0','M1','NA'), df = cp_df)
    stage <- all_sub(colname = "stage", subclass=c('I','II','III','IV','NA'), df = cp_df)

    ##
    tx <- as.data.frame(table(cp_df$stage_edition))
    if(nrow(tx) > 0) {
        tx$com <- paste0(tx[,1]," (",tx[,2],")")
        stage_edition <- paste0(tx$com, collapse = "<br />")
    }

    ##
    tx <- as.data.frame(table(cp_df$ethnicity))
    if(nrow(tx) > 0) {
        tx$com <- paste0(tx[,1]," (",tx[,2],")")
        ethnicity <- paste0(tx$com, collapse = "<br />")
    }

    # endpoint process ----
    endpoint.list <- endpoint_process.func(df=cp_df, column_ids = c("os", "rfs", "dfs", "dss", "pfs", "ttp"))
    endpoint_events <- endpoint.list$endpoint_events
    median_followups <- endpoint.list$median_followups
    zero_status <- endpoint.list$zero_status
    zero_endpoint_status <- endpoint.list$zero_endpoint_status
    missing_endpoint <- endpoint.list$missing_endpoint

    ##
    if(is.null(median_followups)) { 
        median_followups <- "NA"
    } else{
        median_followups <- paste0(median_followups, collapse = "<br />")
    }

    ##
    if(missing_endpoint < 6) {
        if(zero_status == 0) {
            if(is.null(endpoint_events)) {
                endpoint_events <- "NA"
            } else {
                endpoint_events <- paste0(endpoint_events, collapse = "<br />")
            }

            ##
            cp_df[cp_df == "na"] <- "NA"
            
            # combine fixed columns and user custom columns ----
            fixed_colname <- c("patient_id","sample_id","dataset_id","platform","molecular_profiling",
                                "primary_site","disease_type","sample_type","is_primary_treatment",
                                "treatment_type","treatment_setting","regimen",
                                "age","gender","ethnicity","t","n","m","stage_edition","stage",
                                "os_months","os_status","rfs_months","rfs_status",
                                "dfs_months","dfs_status", "dss_months","dss_status",
                                "pfs_months","pfs_status","ttp_months","ttp_status")

            custom_colname <- setdiff(colnames(cp_df), fixed_colname)
            custom_name_id <- which(name1 %in% custom_colname)
            raw_cp_custom_colnames <- name0[custom_name_id]
            fixed_col <- cp_df[, fixed_colname]

            colnames(fixed_col) <- c("Patient_ID","Sample_ID","Dataset_ID","Platform","Molecular_profiling",
                                    "Primary_site","Disease_type","Sample_type","Is_primary_treatment",
                                    "Treatment_type","Treatment_setting","Regimen",
                                    "Age","Gender","Ethnicity","T","N","M","Stage_edition","Stage",
                                    "OS_months","OS_status","RFS_months","RFS_status",
                                    "DFS_months","DFS_status", "DSS_months","DSS_status",
                                    "PFS_months","PFS_status","TTP_months","TTP_status")

            custom_col <- as.data.frame(cp_df[, custom_colname])
            colnames(custom_col) <- raw_cp_custom_colnames
            cp_df <- cbind(fixed_col, custom_col)
            sample_id <- paste0(cp_df$sample_id, collapse = ",")
            
            # custom column ----
            coms <- c()
            for (i in raw_cp_custom_colnames){
                if(!is.numeric(cp_df[,i])) {
                    cp_df[is.na(cp_df[, i]), i] <- "NA"
                    custom_df <- data.frame((table(cp_df[,i])))
                    if(nrow(custom_df)==1 & custom_df[1, "Var1"] == "NA") {} else {
                        custom_df$com <- paste0(custom_df$Var1," (",custom_df$Freq,")")
                        com <- paste0(custom_df$com, collapse = ", ")
                        com <- paste0("<b>",i,"</b>: ",com)
                        coms <- c(coms, com)
                    }
                } else {
                    i_value <- custom_col[,i][!is.na(custom_col[,i])]
                    if(length(i_value) > 0) {
                        if(length(unique(i_value)) > 5) {
                            i_value <- custom_col[,i][!is.na(custom_col[,i])]
                            mean_i <- round(mean(i_value),0)
                            mean_i_text <- paste0("Mean, ",mean_i)
                            median_i <- round(median(i_value),0)
                            median_i_text <- paste0("Median, ",median_i)
                            min_i <- round(min(i_value),0)
                            max_i <- round(max(i_value),0)
                            range_i_text <- paste0("Range, ",min_i,"-",max_i)
                            com <- paste0(c(mean_i_text,median_i_text,range_i_text),collapse = "; ")
                            com <- paste0("<b>",i,"</b>: ",com)
                            coms <- c(coms,com)
                        } else {
                            custom_df <- data.frame((table(custom_col[,i])))
                            custom_df$com <- paste0(custom_df$Var1, " (", custom_df$Freq, ")")
                            com <- paste0(custom_df$com, collapse = ", ")
                            com <- paste0("<b>", i, "</b>: ", com)
                            coms <- c(coms, com)
                        }
                    }
                }
            }
            custom <- paste0(coms, collapse = "\r\n")
            has_anno <- "yes"
        } else {
            zero_endpoint_status0 <- paste0(zero_endpoint_status,collapse = ", ")
            if(length(zero_endpoint_status) > 1) {
                zero_endpoint_status_warning <- paste0(zero_endpoint_status0," columns have no event (i.e. all 0) in the uploaded annotation file.")
            } else {
                zero_endpoint_status_warning <- paste0(zero_endpoint_status0," column has no event (i.e. all 0) in the uploaded annotation file.")
            }
            stop(zero_endpoint_status_warning)
        } # zero status
    } else {
        stop("Endpoint columns were not found.")
    }

    ##
    numberpts <- paste0(nrow(cp_df), " / ", nrow(cp_df))
    sample_id <- paste(cp_df$Sample_ID, collapse = ",")
    ##
    tx <- as.data.frame(table(cp_df$Disease_type))
    if(nrow(tx) > 0) {
        tx$com <- paste0(tx[,1]," (",tx[,2],")")
        disease_type <- paste0(tx$com, collapse = "<br />")
    }
    meta_infor <- rbind(meta_infor, data.frame(d_name, d_name, numberpts,
                                u_primary_site, disease_type, sample_type,
                                d_name, d_name, is_primary_treatment, 
                                treatment_type, treatment_setting, regimen,
                                endpoint_events, median_followups,
                                age, gender, ethnicity, tt, nn, mm,
                                stage_edition, stage, custom, sample_id))
    ##
    user_uploaded_list_2[[z]] <- cp_df
    z=z+1
}

names(user_uploaded_list_2) <- names(user_uploaded_list)
user_uploaded_list <- user_uploaded_list_2
save(user_uploaded_list, file = paste0(b.user.sig.path, "/rda/User_uploaded_anno.RData"))



##
colnames(meta_infor) <- c("dataset_id", "dataset_name", "numberpts",
                          "primary_site", "disease_type", "sample_type",
                          "clinicopathological_info_rdata_variable_name",
                          "profile_matrix_rdata_variable_name", "is_primary_treatment",
                          "treatment_type", "treatment_setting", "regimen", 
                          "endpoint", "median_followup",
                          "age", "gender", "ethnicity", "t", "n", "m", 
                          "stage_edition", "stage", "custom", "sample_id")
                          
##
meta_infor_add <- meta_infor %>% select(dataset_id, dataset_name, numberpts, disease_type,
                      clinicopathological_info_rdata_variable_name,
                      profile_matrix_rdata_variable_name,
                      treatment_type, treatment_setting, regimen, 
                      endpoint, median_followup,
                      age, gender, ethnicity, t, n, m, 
                      stage_edition, stage, custom, sample_id)

user_filtered_datasets_info <- left_join(user_filtered_datasets_info, meta_infor_add, by = "dataset_name")
save(user_filtered_datasets_info, file = paste0(b.user.sig.path, "/rda/User_filtered_datasets.RData"))

## to generate User_parameters.RData ---------------------------------------------------------------------------------
# basic parameters ---
conf_basic_paras <- conf$basic_settings
# parameter names modification 
paras <- sort(names(conf_basic_paras))
conf_basic_paras <- conf_basic_paras[paras]
signature_type = b.signature.type = "Prognostic"

id_map <- c(signature_name = "b.name.sig", 
            endpoint = "b.endpoint", time_intervals_months = "b.times",
            variables_for_independence_test = "b.clinical.variables",
            variables_for_association_analysis = "b.cp.variables",
            number_of_groups = "b.group")

if(!all(paras %in% names(id_map))) { warning(paste0(paras[which(!paras %in% names(id_map))], " was not found!")) }
if(!all(names(conf_basic_paras) == names(id_map[paras]))) { stop("Names not match!!!")}

names(conf_basic_paras) <- id_map[paras]
for(x in names(conf_basic_paras)) { assign(x, conf_basic_paras[[x]]) }


# advanced parameters ---
adv_default <- list(predictor_selection = "Yes", predictor_selection_method = "SPCA",
                    bootstrap_iterations = 200, bootstrap_frequency = 45, batch_correction = NULL,
                    signature_generation_method = "COX", method_to_stratify_patients = "Percentile",
                    number_of_groups = 2, "2groups_high_percentile" = 50, "2groups_low_percentile" = 50, 
                    "3groups_high_percentile" = 75, "3groups_moderate_percentile" = c(25, 75), 
                    "3groups_low_percentile" = 25, variables_for_subgroup_Kaplan_Meier_analysis = NULL, 
                    combine_age_cutoff = NULL, combine_stage_category1 = NULL, combine_stage_category2 = NULL,
                    combine_N_category1 = NULL, combine_N_category2 = NULL, combine_T_category1 = NULL, combine_T_category2 = NULL,
                    variables_for_nomogram = NULL, candidate_predictors_for_signature_wout_molecular_data = NULL,
                    patients_used_for_predictive_signature = NULL, variables_for_interaction_test = NULL,
                    control_regimen = NULL, control_treatment_setting = NULL, control_treatment_type = NULL, 
                    treatment_regimen = NULL, treatment_treatment_setting = NULL, treatment_treatment_type = NULL)
##
if("advanced_settings" %in% names(conf)) {
    conf_adv_paras <- conf$advanced_settings
} else {
    conf_adv_paras <- adv_default
}

conf_adv_paras <- c(conf_adv_paras, adv_default[!names(adv_default) %in% names(conf_adv_paras)])



## parameter names modification 
paras <- sort(names(conf_adv_paras))
conf_adv_paras <- conf_adv_paras[paras]

id_map <- c("b.predictor.selection","b.predictor.selection.methods",
            "b.bootstrap.iterations","b.bootstrap.frequency",
            "b.batch.correction", "signature_generation_method",
            "b.cutpoint.method","b.group",
            "b.high2","b.low2","b.high3","b.moderate3","b.low3",
            "b.subgroup.variable","b.combine.age","b.combine.stage1","b.combine.stage2",
            "b.combine.n1","b.combine.n2","b.combine.t1","b.combine.t2",
            "b.nomogram.variables", "b.clinicopathological.predictors",
            "b.predictive.modeling.data","b.interaction.test",
            "b.control.regimen","b.control.treatment.setting","b.control.treatment.type",
            "b.treatment.regimen","b.treatment.treatment.setting","b.treatment.treatment.type")

names(id_map) <- c("predictor_selection","predictor_selection_method",
                   "bootstrap_iterations","bootstrap_frequency",
                   "batch_correction", "signature_generation_method",
                   "method_to_stratify_patients","number_of_groups",
                   "2groups_high_percentile","2groups_low_percentile",
                   "3groups_high_percentile","3groups_moderate_percentile","3groups_low_percentile",
                   "variables_for_subgroup_Kaplan_Meier_analysis",
                   "combine_age_cutoff","combine_stage_category1","combine_stage_category2",
                   "combine_N_category1","combine_N_category2","combine_T_category1","combine_T_category2",
                   "variables_for_nomogram", "candidate_predictors_for_signature_wout_molecular_data",

                   "patients_used_for_predictive_signature","variables_for_interaction_test",
                   "control_regimen","control_treatment_setting","control_treatment_type",
                   "treatment_regimen","treatment_treatment_setting","treatment_treatment_type")

if(!all(paras %in% names(id_map))) { warning(paste0(paras[which(!paras %in% names(id_map))], " was not found!")) }
if(!all(names(conf_adv_paras) == names(id_map[paras]))) { stop("Names not match!!!")}
names(conf_adv_paras) <- id_map[paras]


for(x in names(conf_adv_paras)) {
    if(x == "b.cutpoint.method") {
        y <- conf_adv_paras[[x]]
        y <- gsub("ReProMSig_defined", "X-tile", y)
        assign(x, y)
    } else if(x == "b.group") {
        y <- conf_adv_paras[[x]]
        y <- paste0(y, " groups")
        assign(x, y)
    } else if(x == "b.combine.age") {
        y <- conf_adv_paras[[x]]
        y <- ifelse(is.null(y), "", as.numeric(y))
        assign(x, y)
    } else if(x == "signature_generation_method") {
        y <- conf_adv_paras[[x]]
        if(y == "COX") {
            LASSO_cox <- FALSE
            standard_cox <- TRUE
        } 
        if(y != "COX") {
            LASSO_cox <- TRUE
            standard_cox <- FALSE
        }
    } else {
        assign(x, conf_adv_paras[[x]])
    }
}

## save out
b.sigtype.roc <- "score"
b.tool.cons <- "nomogram"
b.nomogram.option <- "user"

list0 <- conf_adv_paras[grep("b.combine.", names(conf_adv_paras))]
x <- do.call(c, lapply(names(list0), function(x) { if(!is.null(list0[[x]])) { x } }))
if(!is.null(x)) { x <- gsub("b.combine.|\\d", "", names(x)) %>% unique() %>% capitalize(.) }
b.combine.variables <- x

##
save(b.user.sig.path, b.primary.site, b.sample.type,
     b.have.profiles, b.signature.type, b.name.sig,
     b.upload.predictors, select_training_dataset, select_validation_dataset,
     b.endpoint, b.times,

     b.predictor.selection, b.predictor.selection.methods,
     b.bootstrap.frequency, b.bootstrap.iterations,        
     b.batch.correction, 
     LASSO_cox, standard_cox, b.sigtype.roc,                

     b.clinical.variables, b.cp.variables,
     b.subgroup.variable, b.combine.variables, b.combine.age,
     b.combine.n1, b.combine.n2,
     b.combine.stage1, b.combine.stage2,
     b.combine.t1, b.combine.t2, 
     b.tool.cons, b.nomogram.option, b.nomogram.variables,

     b.cutpoint.method, b.group, b.high2, b.low2,
     b.high3, b.moderate3, b.low3,

     b.predictive.modeling.data, b.interaction.test,
     b.control.regimen, b.control.treatment.setting, b.control.treatment.type,
     b.treatment.regimen, b.treatment.treatment.setting, b.treatment.treatment.type,

     b.clinicopathological.predictors, # variables used to develop signature wout molecular data
     file = paste0(b.user.sig.path, "/rda/User_parameters.RData"))


##
user.config.ini.file <- paste0(b.user.sig.path, "/", "sig.ini")
section <- "Building_signature"
keys <- c("user_filtered_datasets_rdata", "user_uploaded_anno_rdata", "user_parameters_rdata",
          "script_dir", "repromsig_dir", "b_user_sig_path")
values <- c(paste0(b.user.sig.path, "/rda/User_filtered_datasets.RData"),
            paste0(b.user.sig.path, "/rda/User_uploaded_anno.RData"),
            paste0(b.user.sig.path, "/rda/User_parameters.RData"),
            script.dir, repromsig.dir, b.user.sig.path)

for(i in 1:length(keys)) {
    system(paste("python3", paste0(script.dir, "/", "write_ini.py"), user.config.ini.file, section, keys[i], values[i]))
}
