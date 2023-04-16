#
# @Copyright: Peking University Cancer Hospital, All Rights Reserved.
# @Author: Lihua Cao
# @Date: 2021-11-29 10:26:15
# @LastEditTime: 2021-12-20 17:55:03
# @LastEditors: Lihua Cao
# @Description: File content
# Rscript /opt/shiny-server/apps/repromsig/scripts/model.analysis.R /opt/shiny-server/apps/repromsig/ColoGuide_Stage_II_local/output/sig.ini

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
		stop("Usage: Rscript model.analysis.R [user.config.ini.file]\n")
	}
	user.config.ini.file <- args[1]
} else {
    user.config.ini.file <- "/opt/shiny-server/apps/repromsig/ColoGuide_Stage_II_local/output/sig.ini"
}

# 2) Library
library(dplyr)
library(tibble)
library(matrixStats)
library(mfp)
library(purrr)

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
#public_exp_path <- get_value.ccb(config_file = user.config.ini.file,  key = 'public_exp_path')[[1]]
#public_cp_sqlite_path <- get_value.ccb(config_file = user.config.ini.file,  key = 'public_cp_sqlite')[[1]]
user_filtered_datasets.rdata <- get_value.ccb(config_file = user.config.ini.file,  key = 'user_filtered_datasets_rdata')[[1]]
user_uploaded_anno.rdata <- get_value.ccb(config_file = user.config.ini.file,  key = 'user_uploaded_anno_rdata')[[1]]
user_parameters.rdata <- get_value.ccb(config_file = user.config.ini.file,  key = 'user_parameters_rdata')[[1]]

# 4) define/source functions
source(paste0(script.dir,'ccb.helper.R'))

# readin and preprocess  -------------------------------------------------------------------------
load(user_filtered_datasets.rdata)
# cgwtools::lsdata(user_filtered_datasets.rdata)
load(user_parameters.rdata)
# cgwtools::lsdata(user_parameters.rdata)
print(c(select_training_dataset, select_validation_dataset))
# user_path <- dirname(dirname(b.user.sig.path))
user_path <- b.user.sig.path

# load public exp/cp data ------------------------------------------------------------------
sr = "Private dataset"
if(sr == "Public dataset") {
    pub_cp <- DBI::dbConnect(RSQLite::SQLite(), dbname=public_cp_sqlite_path)
    if(b.primary.site == "Stomach") {
        pub_cp_df <- DBI::dbGetQuery(pub_cp, "select * from gc_cp_datasets")
    } else if(b.primary.site == "Colorectal") {
        pub_cp_df <- DBI::dbGetQuery(pub_cp, "select * from crc_cp_datasets")
    }
    DBI::dbDisconnect(pub_cp)
}

# get exp/cp subset datasets ------------------------------------------------------------------
u_cp_list = list()
u_cp_list_ids = c()
user_excluded_cp_list = list()
u_exp_list = list()
u_exp_list_ids = c()

for (i in 1:nrow(user_filtered_datasets_info)) {
    # i = 1
	df <- user_filtered_datasets_info[i, ]
	samid <- strsplit(df$sample_id, ",")[[1]]

	cp_id0 <- paste0(df$dataset_id, "_cp")
	cp_id <- df$clinicopathological_info_rdata_variable_name
    
    if(b.have.profiles == "Yes") {
        exp_id0 <- paste0(df$dataset_id, "_exp")
        exp_id <- df$profile_matrix_rdata_variable_name
    }

	if(df$source %in% "Public dataset") {
        # filter exp subset data
        if(b.have.profiles == "Yes") {
            exp.rda.file <- paste0(public_exp_path, "/", exp_id0, ".RData")
            load(exp.rda.file)
            colnames(public_exp)[1] <- "Symbol"
            expr <- public_exp %>% select("Symbol", samid)
            u_exp_list <- c(u_exp_list, list(expr))
            u_exp_list_ids <- c(u_exp_list_ids, exp_id)
            names(u_exp_list) <- u_exp_list_ids
        }

		# filter cp subset data 
		cp_sub <- subset(pub_cp_df, (pub_cp_df$dataset_id == df$dataset_id) & (pub_cp_df$sample_id %in% samid))
        cp_del <- subset(pub_cp_df, (pub_cp_df$dataset_id == df$dataset_id) & !(pub_cp_df$sample_id %in% samid))

        if(b.primary.site == "Stomach") {
            colnames(cp_sub) <- stomach_pub_cp_colnames
            colnames(cp_del) <- stomach_pub_cp_colnames
        } else if(b.primary.site == "Colorectal") {
            colnames(cp_sub) <- colon_pub_cp_colnames
            colnames(cp_del) <- colon_pub_cp_colnames
        }
        ##
		u_cp_list <- c(u_cp_list, list(cp_sub))
		u_cp_list_ids <- c(u_cp_list_ids, cp_id)
		names(u_cp_list) <- u_cp_list_ids
		user_excluded_cp_list <- c(user_excluded_cp_list, list(cp_del))
		names(user_excluded_cp_list) <- u_cp_list_ids
	} else if(df$source %in% "Private dataset") { 
		# filter exp subset data
        if(b.have.profiles == "Yes") {
            exp.rda.file <- paste0(user_path, "/rda/", exp_id0, ".RData")
            load(exp.rda.file)
            # cgwtools::lsdata(exp.rda.file)
            colnames(expr)[1] <- "Symbol"
            expr <- expr %>% select("Symbol", samid)
            u_exp_list <- c(u_exp_list, list(expr))
            u_exp_list_ids <- c(u_exp_list_ids, exp_id)
            names(u_exp_list) <- u_exp_list_ids
        }

		# filter cp subset data 
        #load(paste0(user_path, "/", "User_uploaded_anno.RData"))
        load(user_uploaded_anno.rdata)
        # cgwtools::lsdata(user_uploaded_anno.rdata)
        cp_id0 <- paste0(df$dataset_id)
		u_cp <- user_uploaded_list[[cp_id0]]
		cp_sub <- subset(u_cp, u_cp$Sample_ID %in% samid)
        cp_del <- subset(u_cp, !u_cp$Sample_ID %in% samid)

        ##
		u_cp_list <- c(u_cp_list, list(cp_sub))
		u_cp_list_ids <- c(u_cp_list_ids, cp_id)
		names(u_cp_list) <- u_cp_list_ids

        ##
        user_excluded_cp_list <- c(user_excluded_cp_list, list(cp_del))
		names(user_excluded_cp_list) <- u_cp_list_ids
	}
    print(i)
}
user_filter_cp_list <- u_cp_list


## continuous_annvar
continuous_annvar <- ""
df <- get_df_summary(signature.type=b.signature.type, cp_df=u_cp_list[[select_training_dataset]], return_data_type=TRUE)
if(nrow(df) > 0) { continuous_annvar <- df %>% filter(Type=="continuous") %>% select(Characteristic) %>% unlist() %>% as.character() } 


# generate sig.ini file ------------------------------------------------------------------
outpath <- b.user.sig.path
mysetwd(paste0(outpath, "/model/"))

user.config.ini.file <- paste0(outpath, "/", "sig.ini")
section <- "Building_signature"
key <- paste0(b.name.sig, "_signature_dir")
value <- outpath
system(paste("python3", paste0(script.dir, "/", "write_ini.py"), user.config.ini.file, section, key, value))
rdata.filename <- paste0(b.name.sig, '.model.analysis.RData')

# Load and select datasets -------------------------------------------
# get expression profiles and cp table for training set ---
cp.name <- user_filtered_datasets_info %>%
                        filter(dataset_name == select_training_dataset) %>%
                        select(clinicopathological_info_rdata_variable_name) %>% unlist()

valid_samid <- strsplit(user_filtered_datasets_info %>%
                        filter(dataset_name == select_training_dataset) %>%
                        select("sample_id") %>% unlist(), ",")[[1]]

## preprocessing the exprmat of training cohort ------------------
if(b.have.profiles == "Yes") {
    exp.name <- user_filtered_datasets_info %>%
                        filter(dataset_name == select_training_dataset) %>%
                        select(profile_matrix_rdata_variable_name) %>% unlist()

    train.exp.log <- user_filtered_datasets_info %>%
                        filter(dataset_name == select_training_dataset) %>%
                        select("log_transform_type") %>%
                        unlist() %>% as.character()

    train.exp <- u_exp_list[[exp.name]] %>% set_names("Symbol", colnames(.)[-1])
    train.exp <- train.exp %>% select(Symbol, intersect(valid_samid, colnames(train.exp)))
    debug(train.exp, 5)

    ## unique by gene
    f <- data.frame(table(train.exp$Symbol))
    r <- f %>% filter(Freq >=2) %>% select(Var1) %>% unlist() %>% as.character()
    s <- f %>% filter(Freq == 1) %>% select(Var1) %>% unlist() %>% as.character()

    if(length(r) > 0) {
        df <- aggregate(train.exp[train.exp$Symbol %in% r, -1], list(train.exp[train.exp$Symbol %in% r, 1]), mean, na.rm = TRUE)
        colnames(df)[1] <- "Symbol"
        train.exp <- rbind(train.exp[train.exp$Symbol %in% s, ], df)
    }
}

# preprocessing the cp table of training cohort ---
df <- u_cp_list[[cp.name]] %>% tibble %>% mycolumn_to_rownames("Sample_ID")
df$Therapy <- df$Treatment_type
if(b.signature.type == "Predictive") {
    df[df$Treatment_type %in% b.treatment.treatment.type & 
       df$Treatment_setting %in% b.treatment.treatment.setting & 
       df$Regimen %in% b.treatment.regimen, "Therapy"] <- "Treatment"
 
    df[df$Treatment_type %in% b.control.treatment.type & 
       df$Treatment_setting %in% b.control.treatment.setting & 
       df$Regimen %in% b.control.regimen, "Therapy"] <- "Control"
}
df <- df %>% filter(!is.na(df[, paste0(b.endpoint, "_months")]) & !is.na(df[,paste0(b.endpoint, "_status")]))

train.cp <- df[intersect(valid_samid, rownames(df)),]
train.cp <- train.cp %>% mutate(Treatment_type = Therapy) %>% select(-Therapy)
train.cp.raw <- train.cp

## for all na treatment column
if(all(is.na(train.cp$Treatment_type)) | all(train.cp$Treatment_type == "NA")) { train.cp$Treatment_type <- "notsure" }
debug(train.cp, 5)

# get validation exp and cp table
valid.exp.list <- list()
valid.log.list <- list()
valid.cp.list <- list()

if(length(select_validation_dataset) >= 1) {
	for(i in 1:length(select_validation_dataset)) {
        # i = 2
		x <- select_validation_dataset[i]
		cp.name <- user_filtered_datasets_info %>%
                         filter(dataset_name == x) %>%
                         select(clinicopathological_info_rdata_variable_name) %>% unlist()
        df <- u_cp_list[[cp.name]]  
        df$Therapy <- df$Treatment_type
        if(b.signature.type == "Predictive") {
            df[df$Treatment_type %in% b.treatment.treatment.type & 
               df$Treatment_setting %in% b.treatment.treatment.setting & 
               df$Regimen %in% b.treatment.regimen, "Therapy"] <- "Treatment"
        
            df[df$Treatment_type %in% b.control.treatment.type & 
               df$Treatment_setting %in% b.control.treatment.setting & 
               df$Regimen %in% b.control.regimen, "Therapy"] <- "Control"
        }

        df <- df %>% filter(!is.na(df[, paste0(b.endpoint, "_months")]) & !is.na(df[,paste0(b.endpoint, "_status")]))
        valid.cp.list[[i]] <- df %>% mutate(Treatment_type = Therapy) %>% select(-Therapy)
		names(valid.cp.list)[i] <- cp.name

        ##
        if(b.have.profiles == "Yes") {                         
            exp.name <- user_filtered_datasets_info %>%
                            filter(dataset_name == x) %>%
                            select(profile_matrix_rdata_variable_name) %>%
                            unlist()

            exp.log <- user_filtered_datasets_info %>%
                            filter(dataset_name == x) %>%
                            select(log_transform_type) %>%
                            unlist() %>% as.character()


            valid.exp.list[[i]] <- u_exp_list[[exp.name]]
            valid.log.list [[i]] <- exp.log
            names(valid.exp.list)[i] <- exp.name
            names(valid.log.list)[i] <- gsub("_exp$", "", exp.name)
        }
        print(i)
	}
}


# Preprocess datasets -------------------------------------------            
colnames(train.cp)[1:length(fixed_variables)] <- fixed_variables ## Note: 1:length(fixed_variables) must be the fixed_variables, note the order
train.cp[train.cp == "NA" | train.cp == "na" | train.cp == "Unknown" | train.cp == "unknown" | train.cp == "[Discrepancy]"] <- NA
train.cp <- train.cp[, which(apply(train.cp, 2, function(x) { n_distinct(unique(x[!is.na(x)]))}) >= 1)]

## remove predictors missing in >= 50% patients for predictive signature
if(b.signature.type == "Predictive") {
    variables0 <- colnames(train.cp)[!colnames(train.cp) %in% c("Treatment_type", exclude_variables)]
    variables0 <- unlist(lapply(variables0, function(x) {if(n_distinct(train.cp[!is.na(train.cp[,x]), x]) > 1) x }))
    variables_exclude <- unlist(lapply(variables0, function(x) {if((nrow(train.cp[is.na(train.cp[,x]),])/nrow(train.cp)) > 0.5) x }))
    train.cp <- train.cp[, !colnames(train.cp) %in% variables_exclude]
}
debug(train.cp)


## b.upload.ann0
x <- colnames(train.cp)
x <- x[!x %in% exclude_variables]
x <- unique(c("Treatment_type", x))
b.upload.ann0 <- train.cp[, intersect(c("SampleID", paste0(b.endpoint, "_months"), paste0(b.endpoint, "_status"), x), colnames(train.cp))]
colnames(b.upload.ann0)[2:3] <- c('time','status')
mode(b.upload.ann0$time) <- "numeric"
mode(b.upload.ann0$status) <- "numeric"


## filter out samples with unvalid follow-up infor
b.upload.ann0 <- b.upload.ann0 %>% filter((!is.na(time) & time >= 0) & (!is.na(status)))
b.upload.ann0$time <- (b.upload.ann0$time*30)
b.upload.ann0[b.upload.ann0$time == 0, "time"] <- 1
colnames(b.upload.ann0) <- gsub(" ", "_", colnames(b.upload.ann0))
debug(b.upload.ann0)

train_surv0 <- b.upload.ann0[, c('SampleID','time','status')]
debug(train_surv0)

## Age and stage process
df <- b.upload.ann0
if(!is.null(b.combine.variables)) {
    age_group <- ifelse("Age" %in% b.combine.variables, TRUE, FALSE)
    if(age_group & ("Age" %in% colnames(df))) { df$Age <- age_comb.func(df$Age, cutoff = as.numeric(b.combine.age)) }

    ##
    stage_group <- ifelse("Stage" %in% b.combine.variables, TRUE, FALSE)
    if(stage_group & ("Stage" %in% colnames(df))) { df$Stage <- stage.func(df$Stage, group1.ids = b.combine.stage1, group2.ids = b.combine.stage2) }

    ##
    TN_ids <- c("T", "N")
    TN_comb_ids <- list(T = list(b.combine.t1, b.combine.t2), N = list(b.combine.n1, b.combine.n2))
    for(x in TN_ids) {
        g1 <- TN_comb_ids[[x]][[1]]
        g2 <- TN_comb_ids[[x]][[2]]
        group <- ifelse(x %in% b.combine.variables, TRUE, FALSE)
        if(group & (x %in% colnames(df))) { df[,x] <- TNM_comb.func(df[,x], group1.ids = g1, group2.ids = g2) }
    }
}
b.upload.ann0 <- df


## patient for subgroup analysis
pts_ann_sb = data.frame()
if(!is.null(b.subgroup.variable)) { pts_ann_sb <- b.upload.ann0[, c("SampleID", b.subgroup.variable)] }

# trans train exp into log2 and trans into row-sample and column-gene format
if(b.have.profiles == "Yes") {
    train <- train.exp %>% tibble %>% mycolumn_to_rownames("Symbol")
    debug(train,5)
    dim(train)
    
    if(train.exp.log == 'nonlog') { df <- data.frame(t(log2((train[, -1] + 1))), check.names = F) }
    if(train.exp.log == 'log10') { df <- data.frame(t(log2((10^(train[, -1]) + 1))), check.names = F) }
    if(train.exp.log %in% c('log2', 'Not applicable')) { df <- data.frame(t(train[, -1]), check.names = F) }
    df <- as.numeric.df(df) %>% as.matrix()
    debug(df, 5)

    # remove genes with sd == 0
    df <- df[, (colSds(df, na.rm = T) > 0)] %>% data.frame(check.names = F) %>% t() # colMeans(log2train0, na.rm = T) > 0
    ## knn impute
    Train0 <- t(knnimpute(df))
} else {
    Train0 <- b.upload.ann0[,-1]
}
debug(Train0, 5)
dim(Train0)

# for predictive signature, split train into test (e.g., chemo) and control (e.g., surgery only) samples ------------
only_treatment_pts_train = FALSE
if(b.signature.type == "Predictive") { only_treatment_pts_train <- ifelse(b.predictive.modeling.data == "treatment", TRUE, FALSE) }

if(b.signature.type == "Predictive") {
    chemo_sams <- train.cp %>% filter(Treatment_type == "Treatment" &
                                      Treatment_setting %in% b.treatment.treatment.setting &
                                      Regimen %in% b.treatment.regimen) %>%
                                      select(SampleID) %>% unlist() %>% as.character()
                                      

    if(only_treatment_pts_train) {
        Train <- Train0[intersect(rownames(Train0), chemo_sams), ]
        train_surv <- train_surv0[intersect(rownames(train_surv0), chemo_sams),]
        b.upload.ann <- b.upload.ann0[intersect(rownames(b.upload.ann0), chemo_sams),]
    } else {
        Train <- Train0
        train_surv <- train_surv0
        b.upload.ann <- b.upload.ann0
    }

    ##
    sg_sams <- train.cp %>% filter(Treatment_type == "Control" &
                                   Treatment_setting %in% b.control.treatment.setting &
                                   Regimen %in% b.control.regimen) %>%
                                   select(SampleID) %>% unlist() %>% as.character()
                                   
    Train_sg <- Train0[intersect(rownames(Train0), sg_sams), ]
    train_surv_sg <- train_surv0[intersect(rownames(train_surv0), sg_sams), ]
    b.upload.ann.sg <- b.upload.ann0[intersect(rownames(b.upload.ann0), sg_sams), ]
} else {
    chemo_sams <- ""
    sg_sams <- ""
    Train <- Train0
    train_surv <- train_surv0
    b.upload.ann <- b.upload.ann0

    Train_sg <- data.frame()
    train_surv_sg <- data.frame()
    b.upload.ann.sg <- data.frame()
}

debug(Train, 5)
debug(train_surv)
debug(b.upload.ann, 5)

##
sams <- intersect.multi(list(rownames(Train), rownames(b.upload.ann), rownames(train_surv)))
Train <- Train[sams, ]
b.upload.ann <- b.upload.ann[sams, ]
train_surv <- train_surv[sams, ]
train_surv_mat <- as.matrix(train_surv[sams, -1])
debug(train_surv_mat)

##
if(b.signature.type == "Predictive") {
    sams <- intersect(rownames(Train_sg), rownames(train_surv_sg))
    train_surv_sg <- train_surv_sg[sams, ]
    Train_sg <- Train_sg[sams, ]
    b.upload.ann.sg <- b.upload.ann.sg[sams,]
}

## Using exprmat of trainset as the reference, ComBat function was applied to reduce the likelihood of batch effects from nonbiological technical biases for each validation dataset profiles ---------------
combat <- ifelse(!is.null(b.batch.correction), TRUE, FALSE)
if(b.have.profiles == "Yes" & combat) {
    if(ncol(Train) > 10) {
        debug(Train, 5)
        train0 <- data.frame(Symbol = colnames(Train), t(Train), check.names=F)
        debug(train0, 5)

        valid.exp.list.combat <- lapply(names(valid.exp.list), function(data_name) {
            # names(valid.exp.list)[1]->data_name
            exp_df <- valid.exp.list[[data_name]] %>%
                                set_names("Symbol", colnames(.)[-1]) %>%
                                tibble %>% mycolumn_to_rownames("Symbol")
            debug(exp_df,5)
            dim(exp_df)

            # trans into log2 for validation dataset
            logtype <- valid.log.list[[data_name]]

            if(grepl("nonlog", logtype)) {df <- log2((exp_df[, -1] + 1))}
            if(grepl("log10", logtype)) {df <- log2((10^(exp_df[, -1]) + 1))}
            if(grepl("log2|Not applicable", logtype)) {df <- exp_df[, -1]}


            # merge two datasets and get batch info ----
            test0 <- data.frame(Symbol = rownames(df), df, check.names=F)
            debug(test0,5)

            ##
            mergedata <- merge(train0, test0, by = "Symbol", all = FALSE) %>%
                                tibble %>% column_to_rownames("Symbol") %>%
                                as.matrix()
            sams <- colnames(train0)[-1]
            sams <- intersect(colnames(mergedata), c(sams, paste0(sams, ".x")))
            batch_info1 <- sams %>% data.frame(check.names = F) %>%
                                    mutate(SampleID = sams, batch = "Training dataset") %>%
                                    select(SampleID, batch)

            sams <- colnames(test0)[-1]
            sams <- intersect(colnames(mergedata), c(sams, paste0(sams, ".y")))
            batch_info2 <- sams %>% data.frame(check.names = F) %>%
                                    mutate(SampleID = sams, batch = data_name) %>% 
                                    select(SampleID, batch)
            batch_infos <- rbind(batch_info1, batch_info2) %>% tibble %>% mycolumn_to_rownames("SampleID")
            res_combat <- sva::ComBat(dat = mergedata, batch_infos$batch,
                                      mod = NULL,
                                      ref.batch = "Training dataset",
                                      par.prior = T)

            # split merged data using sampleID ----
            sams <- colnames(test0)[-1]
            sams <- intersect(colnames(res_combat), c(sams, paste0(sams, ".y")))
            test0 <- res_combat %>% data.frame(check.names = F) %>% select(sams)
            test0 <- data.frame(Symbol = rownames(test0), test0, check.names=F)       
            colnames(test0) <- gsub(".y$", "", colnames(test0))
            debug(test0, 5)
            return(test0)
        })
        names(valid.exp.list.combat) <- names(valid.exp.list)
        valid.exp.list <- valid.exp.list.combat
        print("Combat successful!")
    }
}

# extract candidate genes ---------------
if(!is.null(b.upload.predictors)) {
    Train <- Train[, intersect(unique(as.character(b.upload.predictors[,1])), colnames(Train))] %>% as.matrix()
}

##
top_variance_genes = FALSE
if(top_variance_genes & ncol(Train) > 1000) {
    var_top_genes <- names(rev(sort(apply(Train, 2, var))))[1:1000]
    Train <- Train[, var_top_genes]
}
dim(Train)
debug(Train,5)


## univariate cox for candidate genes or cp predictors ------------------------------------------
web_variables <- c("Type", "Variable", "Events/N", "HR for events (95% CI)", "P")
c_preds <- colnames(Train)
c_preds <- c_preds[!c_preds %in% c("time", "status", "Treatment_type")]
c_preds_unicox <- do.call(rbind.data.frame, lapply(c_preds, function(x) {
    pts_surv = train_surv
    clinic_ann = data.frame(SampleID = rownames(Train), Train[, x], check.names=F) %>% na.omit()
    colnames(clinic_ann)[2] <- x
    if(b.have.profiles == "Yes" | (b.have.profiles == "No" & n_distinct(clinic_ann[,2]) >= 2 & all(table(clinic_ann[,2]) >= 5))) {
        res <- uni_cox.prognosis(pts_surv, clinic_ann)
        mode(res$P) <- 'numeric'
        res <- res[, colnames(res) %in% web_variables]
        colnames(res) <- c("Variable", "Value", "No. of Events", "Hazard Ratio (95% CI)", "P")
        res[which(res$"No. of Events" != ""), c("No. of Events", "No. of Patients")] <- do.call(rbind.data.frame, strsplit(res$"No. of Events", "/"))
        res <- res %>% select("Variable", "Value", "No. of Patients", "No. of Events", "Hazard Ratio (95% CI)","P")
        res[is.na(res)] <- ""
        res
    }
}))

if(b.have.profiles == "Yes") {
    c_preds_unicox <- c_preds_unicox %>% arrange(P)
    if(nrow(c_preds_unicox) > 100) { c_preds_unicox <- c_preds_unicox[1:100, ] }
}

c_preds_unicox <- data.frame("Raw Rank" = 1:nrow(c_preds_unicox), c_preds_unicox, check.names = F)
if(b.have.profiles == "Yes") { c_preds_unicox <- c_preds_unicox %>% select("Raw Rank", "Variable", "Hazard Ratio (95% CI)","P") }


## save key variables ------------------------------------------
#b.cutpoint.method = "X-tile" , b.method, b.predictive.treatment, b.predictive.treatment
save_out = TRUE
if(save_out) {
    if(b.group == "2 groups") {
        save(outpath, select_training_dataset, b.name.sig, b.sample.type, 
             b.have.profiles, b.clinical.variables, b.upload.ann, b.tool.cons,
             b.combine.variables, b.combine.age, b.combine.stage1, b.combine.stage2,
             b.combine.t1, b.combine.t2, b.combine.n1, b.combine.n2, b.nomogram.variables, b.nomogram.option,
             b.control.treatment.type, b.treatment.treatment.type, b.sigtype.roc,
             Train, train_surv, Train_sg, train_surv_sg, b.upload.ann.sg, only_treatment_pts_train,
             train.cp.raw, c_preds_unicox, b.signature.type, valid.exp.list, valid.cp.list, valid.log.list,
             b.cutpoint.method, b.group, b.high2, b.low2, lbs_fun,b.interaction.test,
             b.subgroup.variable, b.cp.variables, b.times, pts_ann_sb, chemo_sams, sg_sams, b.predictor.selection,
             continuous_annvar, user_filter_cp_list, user_excluded_cp_list, file = rdata.filename)
    } else {
        save(outpath, select_training_dataset, b.name.sig, b.sample.type,
            b.have.profiles, b.clinical.variables, b.upload.ann, b.tool.cons,
            b.combine.variables, b.combine.age, b.combine.stage1, b.combine.stage2,
            b.combine.t1, b.combine.t2, b.combine.n1, b.combine.n2, b.nomogram.variables, b.nomogram.option,
            b.control.treatment.type, b.treatment.treatment.type, b.sigtype.roc,
            Train, train_surv, Train_sg, train_surv_sg, b.upload.ann.sg, only_treatment_pts_train,
            train.cp.raw, c_preds_unicox, b.signature.type, valid.exp.list, valid.cp.list, valid.log.list,
            b.cutpoint.method, b.group, b.high3, b.moderate3, b.low3, lbs_fun,
            b.subgroup.variable, b.cp.variables, b.times, pts_ann_sb, b.interaction.test, chemo_sams, sg_sams,
            b.predictor.selection, continuous_annvar,
            user_filter_cp_list, user_excluded_cp_list, file = rdata.filename)
    }
}

## modeling ------------------------------------------
## feature selection ---
boot_n <- b.bootstrap.iterations
# boot_n <- 10
boot_freq <- b.bootstrap.frequency/100
lambda <- ""
coef_lasso_raw <- data.frame()

## step 1. Generate a bootstrap sample, by sampling n individuals with replacement from the original sample.
if(boot_n > 1) {
    samid <- rownames(Train)
    bs_seed <- seq(1, 100000, length = boot_n)
    bs_samid.list <- lapply(bs_seed, function(x) {
        set.seed(x)
        sample(samid, length(samid), replace = TRUE)
    })
    names(bs_samid.list) <- paste0("bs", 1:boot_n)
} else {
    bs_samid.list = list(all_pts = rownames(Train))
}

##
coef1 = data.frame()
fit_cv = NULL
sigfeas = ""
mfp_feas_freq = ""

if(b.have.profiles == "Yes") {
    siggenes = colnames(Train)
    if(!b.name.sig %in% c("ColoGuidePro", "GeneExpressScore")) {
        if(b.predictor.selection == "Yes") {
            # Develop SPCA model using the bootstrap datasets and determine important genes with freq > xx%
            if(b.predictor.selection.methods == "SPCA") {
                print("Start bootstrapped SPCA analysis!")
                spca_genes <- do.call(c, lapply(names(bs_samid.list), function(x) {
                    samid0 <- unique(bs_samid.list[[x]])
                    exp <- Train[samid0, ]
                    surv0 <- train_surv[samid0, ]
                    train_spca <- t(exp)
                    
                    # Fit SPCA model to get important genes from candidate genes
                    status0 = surv0$status
                    imp_feas = ""
                    if(!(all(status0==0)) & all(table(status0) >=5)) {
                        imp_feas <- spca.msig.func(train.exp = train_spca, sample.survival = surv0, unit = "day", n.folds = 10)
                        if(nrow(imp_feas) > 0) { imp_feas = imp_feas$Name }
                    }
                    imp_feas
                }))

                ##
                spca_genes_freq <- (table(spca_genes)/boot_n)
                spca_genes_freq["MT1M"]

                siggenes = names(spca_genes_freq[spca_genes_freq >= boot_freq])
                print("Bootstrapped SPCA successful!")
            }

            # Develop LASSO penalized Cox models using the bootstrap datasets and determine important genes with freq > xx% 
            if(b.predictor.selection.methods == "LASSO") {
                print("Start bootstrapped LASSO penalized Cox analysis!")
                # Develop the LASSO penalized cox model using the bootstrap samples and determine genes with freq > 60%
                lasso_genes <- do.call(c, lapply(names(bs_samid.list), function(x) {
                    samid0 <- bs_samid.list[[x]]
                    exp <- Train[samid0, ]
                    surv0 <- train_surv_mat[samid0, ]
                    
                    # Fit LASSO model to get important genes from candidate genes
                    status0 <- surv0[,"status"]
                    imp_feas <- ""
                    if(!(all(status0==0)) & all(table(status0) >=5)) {
                        model0 <- lasso_cox.func(x = as.matrix(exp), y = surv0)                        
                        # genes with nonzero coefs were signature genes
                        cv_lasso <- model0[[1]]
                        lambda <- model0[[2]]
                        coef0 <- coef(cv_lasso, s = lambda)
                        coef0 <- data.frame(symbol = names(coef0[,1]), coefficients = coef0[,1])
                        coef1 <- coef0[coef0$coefficients != 0,]
                        imp_feas = coef1$symbol
                    } else { imp_feas = "" }
                    imp_feas
                }))
                lasso_genes_freq <- (table(lasso_genes)/boot_n)
                siggenes <- names(lasso_genes_freq[lasso_genes_freq >= boot_freq])
                print("Bootstrapped LASSO penalized Cox analysis successful!")
            }
        }
    }
    siggenes <- siggenes[siggenes != ""]

    # Fit standard Cox to select candidate signature genes -----------------------------------------
    all(rownames(Train) == rownames(train_surv_mat))
    train_model = Train

    if(b.name.sig == "ColoGuidePro") { 
        standard_cox = FALSE
        LASSO_cox = TRUE
    }

    ##
    if(n_distinct(siggenes) > 1) {
        ## standard uni Cox
        pts_surv <- train_surv
        clinic_ann <- data.frame(SampleID = rownames(Train), Train[, siggenes], check.names=F)
        colnames(clinic_ann)[-1] <- siggenes
        uni_cox <- uni_cox.prognosis(pts_surv, clinic_ann)
        multi_cox0 <- multi_cox.prognosis(pts_surv, clinic_ann)
        #
        uni_cox$Coefficient <- multi_cox0[match(uni_cox$Variable, multi_cox0$Variable), "Coefficient"]
        mode(uni_cox$P) <- 'numeric'
        coef1 <- uni_cox
        if(b.name.sig == "ColoGuidePro") { coef1$Coefficient <- ColoGuidePro_coeff }
        if(b.name.sig == "GeneExpressScore") { coef1$Coefficient <- GeneExpressScore_coeff }

        ## standard uni and multi Cox
        if(standard_cox) {
            print("Start standard Cox analysis!")
            coef1 <- uni_cox
            if(nrow(coef1)>=1) {
                # Surgery only patients: none of the seven genes were associated with survival
                if(b.signature.type == "Predictive") {
                    siggenes0 <- coef1$Variable
                    pts_surv <- train_surv_sg
                    clinic_ann <- data.frame(SampleID = rownames(Train_sg), Train_sg[, siggenes0], check.names=F)
                    colnames(clinic_ann)[-1] <- siggenes0
                    uni_cox_sg <- uni_cox.prognosis(pts_surv, clinic_ann)
                    uni_cox_sg <- uni_cox_sg[,c("Variable", "HR for events (95% CI)", "P")]
                    colnames(uni_cox_sg)[-1] <- paste0(colnames(uni_cox_sg)[-1], "\nControl_group")
                    coef1 <- merge.data.frame(coef1, uni_cox_sg, by="Variable", all.x=T)
                }

                ##
                siggenes <- coef1$Variable
                df <- data.frame(Train[, siggenes],check.names = F)
                colnames(df) <- siggenes
                train_bs_sig <- cbind(train_surv_mat, df)
                ##
                train_model <- train_bs_sig
                multi_cox <- survival::coxph(survival::Surv(time, status==1) ~ ., data = train_model)
                fit_cv <- multi_cox
            }
            print("Standard Cox analysis successful!")
        }

        ## LASSO penalized Cox
        if(LASSO_cox) {
            Train <- data.frame(Train[, siggenes],check.names = F)
            colnames(Train) <- siggenes
            train_model <- Train
            
            if(!b.name.sig %in% c("ColoGuidePro", "GeneExpressScore")) {
                # lasso model with selected genes
                model_lasso <- lasso_cox.func(x = as.matrix(train_model), y = train_surv_mat)
                # genes with nonzero coefs were signature genes
                cv_lasso <- model_lasso[[1]]
                lambda <- model_lasso[[2]]
                #lambda <- model_lasso[[1]]$lambda.1se
                coef_lasso <- coef(cv_lasso, s = lambda)
                coef_lasso <- data.frame(symbol = names(coef_lasso[,1]), coefficients = coef_lasso[,1])
                coef_lasso_raw <- coef_lasso
                coef_lasso <- coef_lasso[coef_lasso$coefficients != 0,]
                print(coef_lasso)
                fit_cv <- cv_lasso
            } else {
                if(b.name.sig == "ColoGuidePro") { coef_lasso = data.frame(row.names = siggenes, symbol = siggenes, ColoGuidePro_coeff) }
                if(b.name.sig == "GeneExpressScore") { coef_lasso = data.frame(row.names = siggenes, symbol = siggenes, GeneExpressScore_coeff) }
                coef_lasso_raw <- coef_lasso
            }
            ##
            if(nrow(coef_lasso) > 0) {
                feas <- intersect(coef1[,1], coef_lasso[,1])
                coef1 <- coef1[coef1[,1] %in% feas,]
                coef1[,"Coefficient"] <- coef_lasso[match(coef1[,1], coef_lasso[,1]), 2]
            } else {
                coef1 = data.frame()
            }
        }
    }
} else {
    # preprocess training dataset ---
    debug(Train)
    variables <- colnames(Train)[!colnames(Train) %in% c("time", "status", "Treatment_type")]

    ## only include b.clinicopathological.predictors
    variables <- intersect(b.clinicopathological.predictors, variables)
    variables <- unlist(lapply(variables, function(x) { if(n_distinct(Train[!is.na(Train[,x]),x]) > 1) x }) )
    Train <- Train[, colnames(Train) %in% c("time", "status", "Treatment_type", variables)]

    ## filter predictors with >10% NA
    na_r <- apply(Train[, variables], 2, na_ratio)
    na_ids <- names(na_r[na_r >= 0.1])
    if(length(na_ids) >0) { Train <- Train[,!colnames(Train) %in% na_ids] }
    variables <- names(na_r[na_r < 0.1])

    ##
    res <- do.call(cbind.data.frame, lapply(variables, function(x) {
        x <- Train[, x]
        ifelse(is.numeric(x) & length(x) >= 5, x <- x, x <- as.character(x))
        x
    }))
    colnames(res) <- variables
    Train <- cbind(Train[, c("time", "status")], res)
    debug(Train)

    # start bootstrapped mfp analysis ---
    categ_x <- unlist(lapply(variables, function(x) {if(is.character(Train[,x]) | is.factor(Train[,x])) {x}}))
    formula0 <- do.call(c, lapply(variables, function(x) {
        value <- Train[,x]
        if(is.numeric(value)) {
            out <- paste0("fp(",x,", df = 2, select = 0.05)")
        } else {
            out <- x
        }
        out
    }))

    formula0 <- paste(formula0, collapse = " + ")
    formula0 <- as.formula(paste("Surv(time, status) ~ ", formula0, sep = ""))

    # bootstrapped mfp model ---------------
    if(b.predictor.selection == "Yes") {
        print("Start bootstrapped MFP analysis!")
        mfp_feas <- do.call(c, lapply(names(bs_samid.list), function(x) {
            samid0 <- unique(bs_samid.list[[x]])
            train_mfp <- Train[samid0, ]
            
            # fit mfp model
            mfp_fit0 <- mfp(formula0, family = cox, data = train_mfp)
            sfit0 <- summary(mfp_fit0)
            coef_mfp <- as.data.frame(sfit0$coefficients)
            coef_mfp <- coef_mfp[!is.na(coef_mfp$'Pr(>|z|)') & coef_mfp$'Pr(>|z|)' < 0.05,]
            coef_mfp <- data.frame(rownames(coef_mfp), coef_mfp, check.names = F)
            colnames(coef_mfp) <- c('Variable', 'Coef', "HR", "se_coef", "z", "P")

            ## get id for category variables
            if(nrow(coef_mfp) > 0) {
                for(x in coef_mfp$Variable) {
                    x_rawid <- unlist(lapply(categ_x, function(y) {
                            id0 <- grep(paste0("^", y), x)
                            if(length(id0) >= 1) { y }
                    }))
                    if(!is.null(x_rawid)) {
                        coef_mfp[coef_mfp$Variable == x, "Variable"] <- x_rawid
                    }
                }
            }
            ifelse(nrow(coef_mfp) > 0, imp_feas <- unique(coef_mfp$Variable), imp_feas <- "")
            imp_feas
        }))
        mfp_feas_freq <- (table(mfp_feas)/boot_n)
        sigfeas <- names(mfp_feas_freq[mfp_feas_freq >= boot_freq])
        print("Bootstrapped MFP analysis successful!")
    } else {
        sigfeas <- variables
        mfp_feas_freq <- rep(1, length(variables))
        names(mfp_feas_freq) <- variables
    }

    ## standard cox regression
    train_model <- Train
    if(n_distinct(sigfeas) >= 1) {
        print("Start standard Cox analysis!")
        if(b.signature.type == "Predictive") {
            Train_sg <- Train_sg[, colnames(Train)]
            df <- rbind(Train, Train_sg)
        } else {
            df <- Train
        }
        attach(df)

        ##
        continuous_x <- unlist(lapply(variables, function(x) {if(is.numeric(df[, x])) { x }}))
        df_continuous <- data.frame()

        if(!is.null(continuous_x)) {
            df_continuous <- do.call(cbind.data.frame, 
                                        lapply(continuous_x, function(x) {
                                            id0 <- grep(x, sigfeas)
                                            if(length(id0) >= 1) {
                                                id1 <- match(sigfeas[id0], names(mfp_feas_freq))
                                                formula1 <- names(which.max(mfp_feas_freq[id1]))
                                                res <- data.frame(eval(parse(text = formula1)))
                                                if(b.predictor.selection == "Yes") {
                                                    colnames(res) <- paste0(x, "_mfp_", formula1)
                                                } else {
                                                    colnames(res) <- x
                                                }
                                                res
                                            }
                                        }))
            if(nrow(df_continuous) > 0) { rownames(df_continuous) <- rownames(df) }
            debug(df_continuous)
        }

        ##
        categ_x <- unlist(lapply(variables, function(x) {if(is.character(df[,x]) | is.factor(df[,x])) {x}}))
        df_categ <- data.frame()
        if(!is.null(categ_x)) {
            categ_x <- unlist(lapply(categ_x, function(x) {
                    id0 <- grep(x, sigfeas)
                    if(length(id0) >= 1) { x }
            }))
            df_categ <- data.frame(df[, categ_x], check.names = F)
            if(nrow(df_categ) > 0) {
                colnames(df_categ) <- categ_x
                rownames(df_categ) <- rownames(df)
            }
        }
        list2comb <- list(df_continuous, df_categ)

        ##
        res <- df[,c("time", "status")]
        for(i in 1:length(list2comb)) {
            if(nrow(list2comb[[i]]) > 1) {
                res <- cbind(res, list2comb[[i]])
            }
        }
        df <- res
        detach(df)

        ##
        Train <- df[rownames(Train),  ]
        Train_sg <- df[rownames(Train_sg),  ]
        debug(Train)
        ##
        train_model <- Train
        multi_cox <- survival::coxph(survival::Surv(time,status==1) ~ ., data = train_model)
        fit_cv <- multi_cox
        
        ##
        pts_surv <- data.frame(SampleID = rownames(Train), Train[,c("time", "status")])
        clinic_ann <- data.frame(SampleID = rownames(Train),
                                 Train[,!colnames(Train) %in% c("time", "status")],
                                 check.names = F)
        coef1 <- multi_cox.prognosis(pts_surv, clinic_ann)
        coef1$Type <- do.call(c, lapply(coef1$Type, function(x) { strsplit(x, "_mfp_")[[1]][1] }))
        coef1$Variable <- do.call(c, lapply(coef1$Variable, function(x) {
                                            ifelse(grepl("_mfp_",x), strsplit(x, "_mfp_")[[1]][2],
                                            strsplit(x, "_mfp_")[[1]][1])}))
        print("Standard Cox analysis successful!")
    }
}
print(coef1)
print(fit_cv)



# Evalution in trainset ---------------------------------------------------
if(nrow(coef1) > 0) {
    print("Start evalution in training dataset!")
    ## calculate risk score for chemo trainset 
    require(survival)
    ## remove pts with new level out of training set
    if(b.signature.type == "Predictive" & b.have.profiles == "No") {
        categ_x_sg <- unlist(lapply(colnames(Train_sg)[-(1:2)], function(x) {
                                        if(is.character(Train_sg[,x]) | is.factor(Train_sg[,x])) {x}
                                        }
                                    )
                            )
        for(x in categ_x_sg) {
            Train_sg <- Train_sg[Train_sg[,x] %in% unique(Train[,x]),]
        }
    }

    ##
    if(b.signature.type == "Predictive" & only_treatment_pts_train) {
        Train_sg <- Train_sg[, colnames(Train)]
        df <- rbind(Train, Train_sg)
    } else {
        df <- Train
    }
    traintopred <- df

    ##
    if(b.have.profiles == "Yes") {
        print(all(coef1[,1] %in% colnames(df)))
        if(standard_cox) {
            require(survival)
            risk_score <- data.frame(SampleID = rownames(df), RS = predict(fit_cv, type = "lp", newdata = as.data.frame(df)))
        } else if(LASSO_cox & !(b.name.sig %in% c("ColoGuidePro", "GeneExpressScore"))) {
            risk_score <- data.frame(SampleID = rownames(df), RS = as.numeric(predict(fit_cv, newx = as.matrix(df), s = lambda, type='link')))
        } else if(b.name.sig %in% c("ColoGuidePro", "GeneExpressScore")) {
            df <- df[,coef1$Variable]
            risk_score <- data.frame(SampleID = rownames(df), RS = apply(df, 1, function(x) {
                x <- as.numeric(x)
                sum(x * coef1$Coefficient)
            }))
        }
    } else {
        risk_score <- data.frame(SampleID = rownames(df), RS = predict(fit_cv, type = "lp", newdata = as.data.frame(df)))
    }


    ## remove pts with NA RS
    risk_score <- risk_score %>% filter(!is.na(RS))
    rs0 <- risk_score[,2]
    names(rs0) <- risk_score$SampleID

    # KM plot split by risk score level from user setted cutoff --
    g1 <- ifelse(b.signature.type == "Predictive", "Benefit", "High risk")
    g2 <- ifelse(b.signature.type == "Predictive", "No-benefit", "Low risk")
    g3 <- ifelse(b.signature.type == "Predictive", "Intermediate group", "Intermediate risk")

    ##
    b.high2_xtile = ""
    b.low2_xtile = ""
    ##
    if(b.cutpoint.method == "Percentile") {
        if(b.group == "2 groups") {
            high <- as.numeric(b.high2)/100
            moderate <- NULL
            low <- as.numeric(b.low2)/100
        } else {
            high <- as.numeric(b.high3)/100
            moderate <- as.numeric(b.moderate3)/100
            low <- as.numeric(b.low3)/100
        }
    } else {
        p_x <- seq(10, 90, by = 2)/100
        names(p_x) <- seq(10, 90, by = 2)
        rs_x <- data.frame(RS = rs0)

        p_xtiles = do.call(rbind.data.frame, lapply(names(p_x), function(x) {
            rs_x0 <- rs_x
            cutp <- quantile(rs0, probs = p_x[x])
            rs_x0[rs_x0$RS > cutp, 'RS_group'] <- g1
            rs_x0[rs_x0$RS <= cutp, 'RS_group'] <- g2
            rs_x0 <- rs_x0 %>% filter(!is.na(RS_group)) %>%
                               rownames_to_column(var = "SampleID") %>%
                               select(SampleID, RS_group)
            ##
            rs_x1 <- rs_x0$RS_group
            names(rs_x1) <- rs_x0$SampleID
            leg0 <- paste(paste0(names(table(rs_x1)), ' (n = '), table(rs_x1), ')',sep='')
            surv_df <- train_surv[, c("time", "status")]
            res <- data.frame(cutoff = x, survival_1VS2_diffcal.fun(surv_df, rs_x1))
            res
        }))

        p_xtiles <- p_xtiles %>% filter(HR != "Inf" & HR != "-Inf" & number1 > 10 & number2 > 10) %>%
                                mutate(p = as.numeric(p), HR = as.numeric(HR)) %>%
                                filter(HR > 1)

        cutoff <- p_xtiles[which.min(p_xtiles$p), 1] %>% as.numeric()
        b.high2_xtile <- cutoff
        b.low2_xtile <- cutoff
        ##
        high <- as.numeric(cutoff)/100
        moderate <- NULL
        low <- as.numeric(cutoff)/100
        write.file(p_xtiles, "p_xtiles.csv")
    }

    ##
    top = quantile(rs0, probs = high)
    bottom = quantile(rs0, probs = low)
    top_train = top
    bottom_train = bottom
    ##
    rs0 <- data.frame(RS = rs0)
    if((top == bottom) | is.null(moderate))  {
        rs0[rs0$RS > top, 'RS_group'] <- g1
        rs0[rs0$RS <= bottom, 'RS_group'] <- g2
    } else if(top > bottom & !is.null(moderate)) {
        rs0[rs0$RS > top, 'RS_group'] <- g1
        rs0[rs0$RS > bottom & rs0$RS <= top, 'RS_group'] <- g3
        rs0[rs0$RS <= bottom, 'RS_group'] <- g2
    }
    
    ##
    rs0 = rs0 %>% filter(!is.na(RS_group))
    rs1 = rs0$RS_group
    names(rs1) = rownames(rs0)
    stratify <- data.frame(SampleID = names(rs1), Type = rs1, signature = b.name.sig)
    
    # Clinicopathological Association  --------------------------------------
    cp_rs <- stratify %>% filter(signature == b.name.sig) %>% select(SampleID, Type)
    cp_rs <- merge.data.frame(cp_rs, b.upload.ann0, by='SampleID', all.x=T)
    cp_rs <- cp_rs[!is.na(cp_rs$Type),]
    debug(cp_rs)
    ##
    train_stratify <- cp_rs
    debug(train_stratify)
    score0 <- rs0[match(train_stratify$SampleID, rownames(rs0)),]
    write.file(data.frame(train_stratify, score = score0), file = paste0(b.name.sig, ".risk_score.csv"))

    ##
    d_id <- "Training dataset"
    train_CP_Stat_report <- data.frame()
    variables <- sort(colnames(cp_rs)[!colnames(cp_rs) %in% c("SampleID", "Type", "time", "status")])

    if(!is.null(b.cp.variables)) {
        variables <- intersect(variables, b.cp.variables)
        ifelse(b.signature.type == "Predictive",
               pts_type <- c("Training dataset", sort(unique(cp_rs$Type))), 
               pts_type <- "Training dataset")

        ##
        train_CP_Stat.list = list()
        key_df = data.frame()
        include_continous = FALSE

        for(type in pts_type) {
            ifelse(type == "Training dataset", cp_rs0 <- cp_rs, cp_rs0 <- cp_rs %>% filter(Type == type))
            clust_column <- ifelse(type == "Training dataset", "Type", "Treatment_type")

            Stat <- lapply(variables, function(x) {
                # variables[1] -> x 
                rs0 <- cp_rs0[, c(clust_column, x)]
                rs0 <- rs0[!is.na(rs0[,1]) & !is.na(rs0[,2]), ]
                rs0[is.na(rs0)] <- "Unknown"
                rs0 <- rs0 %>% na.omit() %>% arrange(rs0[,1])
                null0 <- data.frame(matrix("", nrow=1, ncol = n_distinct(rs0[,1]) + 2,
                                           dimnames = list("1", c("Variable", "Level", unique(rs0[,1])))),
                                           check.names=F)
                null0[1, 1:2] <- x

                if(n_distinct(rs0[,1]) > 1 & n_distinct(rs0[,1]) < 6) {
                    if(n_distinct(rs0[,2]) > 1 & n_distinct(rs0[,2]) < 6) {
                        rs1 <- t(as.data.frame.matrix(table(rs0)))
                        rs2 <- rbind(null0, data.frame(Variable = x, Level = rownames(rs1), rs1, check.names = F))

                        ft.mat <- rs1[!rownames(rs1) %in% 'Unknown',]
                        if(class(ft.mat) == "matrix") {
                            ft <- fisher.test(ft.mat, workspace = 2e5, simulate.p.value = TRUE, B = 1e5)
                            pval0 = ifelse(ft$p.value < 0.01,
                                           format(ft$p.value, digits=2, scientific=TRUE),
                                           round(ft$p.value, 3))
                            rs2 <- data.frame(rs2, P = c(pval0, rep('', (nrow(rs2)-1))), check.names=F)
                            rs2
                        }
                    } else if(include_continous & n_distinct(rs0[,1]) == 2 & is.numeric(rs0[,2]) & n_distinct(rs0[,2]) >= 6) {
                        rs1 <- split(rs0, rs0[,clust_column])
                        x0 <- rs1[[1]][,2]
                        y0 <- rs1[[2]][,2]
                        pval0 <- t.test(x0, y0)$p.value
                        pval0 <- ifelse(pval0<0.01,
                                        format(pval0, digits=2, scientific=TRUE),
                                        round(pval0, 3))
                        null0[,3] <- round(mean(x0, na.rm = TRUE), 3)
                        null0[,4] <- round(mean(y0, na.rm = TRUE), 3)
                        rs2 <- data.frame(null0, P = pval0, check.names=F)
                        rs2
                    }
                }
            })

            res <- do.call(rbind.data.frame, Stat)
            if(b.signature.type == "Predictive" & type != pts_type[1]) {
                colnames(res)[-(1:2)] <- paste0(type, "<br />", colnames(res)[-(1:2)])
            }
            train_CP_Stat.list <- c(train_CP_Stat.list, list(res))
            key_df <- rbind(key_df, res[, c("Variable", "Level")])
        }
        
        key_df <- unique(key_df)
        ids <- sort(unique(key_df$Variable))
        key_df <- do.call(rbind.data.frame, lapply(ids, function(x) { key_df[key_df$Variable == x, ] }))
        train_CP_Stat_report <- left_join.multi(key.df = key_df, ID = c("Variable", "Level"), list.df = train_CP_Stat.list)
        colnames(train_CP_Stat_report)[-(1:2)] <- paste0(colnames(train_CP_Stat_report)[-(1:2)], "<br />", d_id)
        
        if(b.signature.type == "Predictive") {
            n_g = ifelse(b.group == "2 groups", 2, 3)
            if(n_g == 2) {
                s_ids <- sort(c(seq(6, 13, by = 3)[1:n_g], (seq(6, 13, by = 3)+1)[1:n_g]))
            } else {
                s_ids <- sort(c(seq(7, 15, by = 3)[1:n_g], (seq(7, 15, by = 3)+1)[1:n_g]))
            }
            p_sub <- apply(train_CP_Stat_report[, s_ids], 1, function(x) {
                x <- as.numeric(x)
                if(sum(is.na(x)) == 0) {
                    p0 <- fisher.test(matrix(x, nrow = n_g))$p.value
                    p0 <- ifelse(p0 < 0.01, format(p0, digits = 2, scientific = TRUE), round(p0, 3))
                } else { "" }
            })
            train_CP_Stat_report <- data.frame(train_CP_Stat_report, p_sub, check.names = F)
            colnames(train_CP_Stat_report)[ncol(train_CP_Stat_report)] <- paste0("P-subgroup", "<br />", d_id)
            train_CP_Stat_report[train_CP_Stat_report == "Treatment_type"] <- "Therapy"
        }
        debug(train_CP_Stat_report)
    }

    # Estimating x-year survival -------------------------------------------
    surv_prob <- train_stratify %>% select(c("SampleID", "time", "status", "Type", "Treatment_type"))
    types <- sort(unique(surv_prob$Type))
    h1 <- ifelse(b.signature.type == "Predictive", "Group", "Signature group")

    surv_prob.list <- lapply(types, function(x) {
        df0 <- surv_prob[surv_prob$Type == x, ]
        if(b.signature.type == "Predictive") {
            df0[df0$Treatment_type %in% b.treatment.treatment.type & 
                df0$Treatment_setting %in% b.treatment.treatment.setting & 
                df0$Regimen %in% b.treatment.regimen, "Treatment_type"] <- "Treatment"
        
            df0[df0$Treatment_type %in% b.control.treatment.type & 
                df0$Treatment_setting %in% b.control.treatment.setting & 
                df0$Regimen %in% b.control.regimen, "Treatment_type"] <- "Control"
            treat_types <- sort(unique(df0$Treatment_type))
        }

        ## for type
        N <- nrow(df0)
        n <- nrow(df0[df0$status == 1, ])
        N <- paste0(N, " (", n, ")")
        ##
        yrs <- b.times/12
        res1 <- do.call(cbind, lapply(yrs, function(n0) {
            if(max(df0$time) > 365*n0) {
                sf0 <- summary(survfit(Surv(time, status) ~ 1, data = df0), times = 365*n0)
                surv_p = round(sf0$surv*100, 2)
                surv_ci = paste0(round(sf0$lower * 100, 2), " - ", round(sf0$upper * 100, 2))
                surv_ci = paste0(surv_p, " (", surv_ci, ")")
            } else { surv_ci = NA }
            surv_ci
        }))

        res1 <- data.frame(x, x, N, res1, check.names = F)
        colnames(res1) <- c("Key", h1, paste0("All Patients (No. of Events)<br />", d_id),
                            paste0(b.times, "-months Survival Probability (%, 95% Cl)<br />", d_id))

        ## for predictive signature
        if(b.signature.type == "Predictive") {
            for(treat_type0 in treat_types) {
                df1 <- df0[df0$Treatment_type == treat_type0, ]
                N1 = paste0(nrow(df1), " (", nrow(df1[df1$status == 1, ]), ")")
                res2 <- do.call(cbind, lapply(yrs, function(n0) {
                    if(max(df1$time) > 365*n0) {
                        sf1 <- summary(survfit(Surv(time, status) ~ 1, data = df1), times = 365.25*n0)
                        surv_ci1 = paste0(round(sf1$surv*100, 1), " (",
                                          paste0(round(sf1$lower * 100, 2), " - ",
                                          round(sf1$upper * 100, 2)), ")")
                    } else { surv_ci1 = NA }
                    surv_ci1
                }))
                res2 <- data.frame(paste0(x, "_", treat_type0), treat_type0, N1, res2, check.names = F)
                colnames(res2) <- c("Key", h1, paste0("All Patients (No. of Events)<br />", d_id),
                                     paste0(b.times, "-months Survival Probability (%, 95% Cl)<br />", d_id))
                res1 <- rbind(res1, res2)
            }
        }
        if(b.signature.type == "Predictive") { res1 <- data.frame("Predictive group" = x, res1, check.names = F) }
        res1
    })

    train_surv_prob_report <- do.call(rbind.data.frame, surv_prob.list)
    if(save_out) {
        cgwtools::resave(standard_cox, LASSO_cox, lambda, coef_lasso_raw, fit_cv,
                        traintopred, risk_score, top, bottom, top_train, bottom_train,
                        b.high2_xtile, b.low2_xtile,
                        train_stratify, train_CP_Stat_report,
                        train_surv_prob_report, file = rdata.filename)
    }
    message_nosignature <- "Successful!"
    print("Evalution in training dataset successful!")
} else {
    message_nosignature <- "There was no signature detected!"
}

## modify header
if(nrow(coef1) > 0) {
    coef_report <- coef1 %>% select(Type, Variable, Coefficient, "HR for events (95% CI)", P)
    colnames(coef_report) <- c("Variable", "Characteristic", "Coefficient", "HR (95% CI)", "P")
    coef_report <- data.frame("Raw Rank" = 1:nrow(coef_report), coef_report, check.names=F)
    if(b.have.profiles == "Yes" & b.signature.type == "Prognostic") {
        coef_report <- coef_report %>% mutate(P = ifelse(P < 0.01, format(P, digits = 2, scientific = TRUE), round(P, 3))) %>%
                                       mutate(P = as.character(P))
    }
    if(b.have.profiles == "Yes") { coef_report <- coef_report %>% select(-Variable) }
} else {
    coef_report <- data.frame()
}

##
if(save_out) {
    cgwtools::resave(train_model, coef1, sigfeas, mfp_feas_freq, coef_report, b.endpoint, message_nosignature, file = rdata.filename)
}
section = 'Building_signature'
key = 'model_analysis_rdata'
value <- paste0(paste0(outpath, '/model/'), rdata.filename)
system(paste("python3", paste0(script.dir, "/", "write_ini.py"), user.config.ini.file, section, key, value))


# Generate summary plot data for trainset and validation datasets -------------------------------------------------------------------------
cp.list <- c(list("Training dataset" = train.cp), valid.cp.list)
b.sum_data.list <- lapply(names(cp.list), function(data_name) {
    pts.ann <- cp.list[[data_name]]
    pts.ann[pts.ann == "NA" | pts.ann == "na" | pts.ann == "Unknown" | pts.ann == "unknown" | pts.ann == "[Discrepancy]"] <- NA
    pts.ann <- pts.ann[, which(apply(pts.ann, 2, function(x) { n_distinct(unique(x[!is.na(x)])) }) >= 2)]
    df <- pts.ann[,!colnames(pts.ann) %in% exclude_variables]
    debug(df)

    if(class(df) == "data.frame") {
        if(!is.null(b.combine.variables)) {
            age_group <- ifelse("Age" %in% b.combine.variables, TRUE, FALSE)
            if(age_group & ("Age" %in% colnames(df))) { df$Age <- age_comb.func(df$Age, cutoff = as.numeric(b.combine.age)) }

            ##
            stage_group <- ifelse("Stage" %in% b.combine.variables, TRUE, FALSE)
            if(stage_group & ("Stage" %in% colnames(df))) { df$Stage <- stage.func(df$Stage, group1.ids = b.combine.stage1, group2.ids = b.combine.stage2) }

            ##
            TN_ids <- c("T", "N")
            TN_comb_ids <- list(T = list(b.combine.t1, b.combine.t2), N = list(b.combine.n1, b.combine.n2))
            for(x in TN_ids) {
                g1 <- TN_comb_ids[[x]][[1]]
                g2 <- TN_comb_ids[[x]][[2]]
                group <- ifelse(x %in% b.combine.variables, TRUE, FALSE)
                if(group & (x %in% colnames(df))) { df[,x] <- TNM_comb.func(df[,x], group1.ids = g1, group2.ids = g2) }
            }
        }

        # Generate summary table and plot data 
        b.pie_plot_data <- apply(df[, !colnames(df) %in% "SampleID"], 2, function(x) {
                                    x <- data.frame(table(x))
                                    x <- data.frame(x, Fraction = paste0(round(x$Freq/sum(x$Freq), 2) * 100, '%'))
                                    colnames(x) <- c('Variable', 'Freq', 'Fraction')
                                    x
                                    }
                                )
        b.sum_data <- c(list(pts_infor = pts.ann), b.pie_plot_data)
    } else {
        df <- pts.ann[, colnames(pts.ann) %in% exclude_variables]
        b.sum_data <- c(list(pts_infor = df), list())
    }
    b.sum_data
})

names(b.sum_data.list) <- names(cp.list)
len0 <- lapply(b.sum_data.list, function(x) { length(x)})
b.sum_data.list <- b.sum_data.list[which(len0 > 0)]

# Generate log2 expression profiles data for trainset and validation datasets-------------------------------------------------------------------------
if(b.have.profiles == "Yes") {
    if(length(valid.exp.list) >= 1) {
        testsets.name <- gsub("_exp$", '', names(valid.exp.list))
        exp.list <- c(list(Train0), valid.exp.list)
        names(exp.list) <- c("Training dataset", testsets.name)
    } else {
        exp.list <- list("Training dataset" = Train0)
    }

    b.profile_data.list <- lapply(names(exp.list), function(data_name) {
        exp_df <- exp.list[[data_name]]
        debug(exp_df,5)
        if(data_name != "Training dataset") {          
            # trans into log2 for validation dataset
            logtype <- valid.log.list[[data_name]]
            if(grepl('nonlog', logtype)) {logd2exp <- log2((exp_df[,-1]+1))}
            if(grepl('log10', logtype)) {log2exp <- log2((10^(exp_df[,-1])+1))}
            if(grepl('log2|Not applicable', logtype)) { log2exp <- exp_df[,-1] }
            log2exp <- data.frame(t(log2exp), check.names = F)
            log2exp[is.na(log2exp)] <- mean(apply(log2exp, 1, mean, na.rm=T))
        } else {
            log2exp <- exp_df
        }
        debug(log2exp,5)

        # remove outliers
        x <- as.numeric(as.matrix(log2exp))
        qnt <- quantile(x, probs=c(.25, .75), na.rm = T)
        caps <- quantile(x, probs=c(.05, .95), na.rm = T)
        H <- 1.5 * IQR(x, na.rm = T)
        x[x < (qnt[1] - H)] <- NA
        x[x > (qnt[2] + H)] <- NA
        x <- x[!is.na(x)]
        sample(x, 10000, replace=T)
    })
    names(b.profile_data.list) <- names(exp.list)
} else {
    b.profile_data.list <- list()
}
if(save_out) {
    cgwtools::resave(b.sum_data.list, b.profile_data.list, file = rdata.filename)
}
# lsdata(rdata.filename)

print(message_nosignature)


