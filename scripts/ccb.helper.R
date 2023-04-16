# Read in and write out -----------------------------------------------------------
# Author: Lihua cao. Create date: Feb 2020. support xlsx, csv and plain text format.
## support xlsx, csv, maf and plain text format, if your file contain "#",please careful, and check the row and column number
read.file <- function(file,startRow=1,header=TRUE) {
    if(!file.exists(file)) {
        warning("File was not found!")
    } else {
        if(grepl("\\.xlsx$", file)) {
            require(openxlsx)
            y <- openxlsx::read.xlsx(file,colNames = TRUE,startRow = startRow,check.names=FALSE)
        } else if(grepl("\\.csv$", file)) {
            y <- read.csv(file,header = header, check.names=FALSE,fill=TRUE)
        } else if(grepl("\\.maf$", file)) {
            y <- data.frame(data.table::fread(input = file, sep = "\t", stringsAsFactors = FALSE, verbose = FALSE, 
            data.table = TRUE, showProgress = TRUE, header = TRUE, fill = TRUE, skip = "Hugo_Symbol"))
        } else {
            y<-read.delim(file,check.names=FALSE,header=header,fill=TRUE)
        }
        return(y)
    }
}


## support xlsx, csv, maf and plain text format
write.file <- function(x,file,quote=FALSE,col.names=TRUE,row.names=FALSE,sep='\t') {
    if(grepl("\\.xlsx$", file)) {
        require(openxlsx)
        openxlsx::write.xlsx(x=x, file=file, row.names = row.names, col.names = col.names)
    } else if(grepl("\\.csv$", file)) {
        write.csv(x=x,file=file,row.names=row.names)
    } else {
	write.table(x=x,file=file,quote=quote,col.names=col.names,row.names=row.names,sep=sep)
    }
}


# config ini read in Author: Yang DU. Create date: March 2020. 
## 192 web; Modified by lihua cao.
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


##
get_value_tripod.ccb <- function(config_file, key = NULL, space = FALSE) {
  if(!file.exists(config_file)) {
    warning("Config ini file was not found.")
  } else {
    ini_df <- read.table(config_file, fill = TRUE, header = FALSE,sep='\t',stringsAsFactors=FALSE)
    ini_df <- ini_df[grep(pattern = "^\\[", x = ini_df[, 1], invert = T), ]
    ini_df <- do.call(rbind.data.frame,strsplit(ini_df,'='))
    colnames(ini_df) <- c('key','value')
    ini_df$key <- gsub(" ", "", ini_df$key, fixed = TRUE)
    if (space == TRUE){
      ini_df$value <- gsub(" ", "", ini_df$value, fixed = TRUE)
    }
    if(any(duplicated(ini_df$key))) {
      stop("exist duplicated keys in your ini file, please check")
    }
    rownames(ini_df) <- ini_df$key
    if(is.null(key)) {
      stop("key is not specified.")
    } else {
      return(as.character(ini_df[key,'value']))
    }
  }
}


# Basic calculation or process ----------------------------------------------------
mysetwd <- function(dir, final=FALSE) {
  if(!file.exists(dir))
    dir.create(dir, recursive = TRUE)
  else if (final)
    stop(paste0("The ouput dir ", dir, " exists. Have you done the analysis already?"))
  setwd(dir)
}

##
h <- function(df,limit='all') {
  if(limit=='all') {
    print(head(df))
  } else {
    print(head(df[1:limit, 1:limit]))
  }
}

debug.message <- function(x) if(DEBUG_MODE) message(x)

debug <- function(x,limit='all') if(DEBUG_MODE) h(x,limit=limit)

##
na_ratio<-function(x) {
  ((length(x[is.na(x) | x=="NA" | x=="na"]))/length(x))->y
  return(y)
}

# Matrix/data.frame processing  ---------------------------------------------------
as.numeric.df <- function(df) {
  if(nrow(df)>1) {
    df.rst <- apply(df, 2, as.numeric)
  } else {
    df.rst <- data.frame(matrix(apply(df, 2, as.numeric),nc=ncol(df)))
  }
  rownames(df.rst) <-rownames(df)
  colnames(df.rst) <-colnames(df)
  df.rst
}

##
intersect.multi <- function(in.list) {
  for(i0 in 1:(length(in.list)-1)) {
    intersect(in.list[[i0]],in.list[[i0+1]])->in.list[[i0+1]]
  }
  y<-in.list[[length(in.list)]]
  return(y)
}

# Function: left_join multi data.frames (>=3)
left_join.multi <- function(key.df, ID, list.df, keep = FALSE, suffix = c(".x", ".y")) {
  require(dplyr)
  a = key.df
  y <- list()
  for(i0 in 1:length(list.df)) {
    # i0 =1
    b = list.df[[i0]] 
    ids <- names(b)
    if(any(duplicated(ids))) {
        dupids <- ids[duplicated(ids)]
        for(x in dupids) {
            x0 = ids[ids==x]
            ids[ids==x] <- paste0(x, "ccb_suffix_", 1:length(x0))
        }
    }
    names(b) <- ids
    res <- left_join(a, b, by = ID, keep = keep, suffix = suffix)
    if(i0 > 1) { res <- res %>% select(colnames(res)[!colnames(res) %in% ID]) }
    y <- c(y, list(res))
  }
  y <- do.call(cbind.data.frame, y)
  names(y) <- gsub("ccb_suffix_\\d", "", names(y))
  return(y)
}

## impute expmat with NA, gene name as rownames
knnimpute <- function(data){
  m <- data
  debug.message(paste("Before knnimpute:", nrow(data),"rows,", ncol(data), "columns"))
  h <- impute::impute.knn(as.matrix(m))
  h <- h$data
  rownames(h) <- rownames(m)
  debug.message(paste("After knnimpute:", nrow(h),"rows,", ncol(h), "columns"))
  return(h)
}

##
mycolumn_to_rownames <- function (.data, var = "rowname") {
    require(dplyr)
    require(rlang)
    stopifnot(is.data.frame(.data))
    if (has_rownames(.data)) {
        cnd_signal(error_already_has_rownames())
    }
    if (!has_name(.data, var)) {
        cnd_signal(error_unknown_column_names(var))
    }
    .data <- as.data.frame(.data)
    rownames(.data) <- .data[[var]]
    #.data[[var]] <- NULL
    .data
}

## group age
age_comb.func <- function(age, cutoff = 60) {
    mode(age) <- 'numeric'
    age[!is.na(age) & age >= cutoff] <- paste0('>=', cutoff)
    age[!is.na(age) & age != paste0('>=', cutoff)] <- paste0('<', cutoff)
    return(age)
}


## group stage
stage.func <- function(stage, group = TRUE, group1.ids = NULL, group2.ids = NULL) {
        stage <- gsub('stage ', '', stage)
        stage <- toupper(stage)
        stage[!is.na(stage) & stage %in% c("IA","IB","IC",'1A','1B','1C','1')] <- 'I'
        stage[!is.na(stage) & stage %in% c("IIA","IIB","IIC",'2A','2B','2C','2')] <- 'II'
        stage[!is.na(stage) & stage %in% c("IIIA","IIIB","IIIC",'3A','3B','3C','3')] <- 'III'
        stage[!is.na(stage) & stage %in% c("IVA","IVB","IVC",'4A','4B','4C','4')] <- 'IV'
        stage[!is.na(stage) & !stage %in% c('I','II','III','IV')] <- "Other"
        if(group & !is.null(group1.ids) & !is.null(group2.ids)) {
            stage[!is.na(stage) & stage %in% group1.ids] <- paste(group1.ids, collapse = "+")
            stage[!is.na(stage) & stage %in% group2.ids] <- paste(group2.ids, collapse = "+")
        }
    return(stage)
}

## group T/N/M
TNM_comb.func <- function(x, group = TRUE, group1.ids = NULL, group2.ids = NULL) {
  if(group & !is.null(group1.ids) & !is.null(group2.ids)) {
      x[!is.na(x) & x %in% group1.ids] <- paste(group1.ids, collapse = "+")
      x[!is.na(x) & x %in% group2.ids] <- paste(group2.ids, collapse = "+")
  } else {
    x
  }
  return(x)
}

## from tt
age_gender_race_process.func <- function(df, column_ids = c("age", "gender", "ethnicity")) {
    # df = df; column_ids = c("age", "gender", "ethnicity")
    for(x in column_ids) {
        # x = "age"
        if(x == "ethnicity") { id0 <- which(colnames(df) %in% c("race", "ethnicity")) }
        if(x != "ethnicity") { id0 <- which(colnames(df) %in% x) }

        if(length(id0) >= 1) {
            colnames(df)[id0] <- x
            if(x == "age") {
                mode(df[,x]) <- "numeric"
            } else if(x == "gender") {
                df[,x] <- toupper(df[,x])
                df[!is.na(df[,x]) & df[,x] %in% c("F", "FEMALE"), "gender"] <- "Female"
                df[!is.na(df[,x]) & df[,x] %in% c("M", "MALE"), "gender"] <- "Male"
                df[is.na(df[,x]) | !df[,x] %in% c("Female", "Male"), "gender"] <- "NA"
            } else if(x == "ethnicity") {
                df[,x] <- tolower(df[,x])
                df[!is.na(df[,x]) & df[,x] %in% c("asian"), "ethnicity"]<- "Asian"
                df[!is.na(df[,x]) & df[,x] %in% c("white"), "ethnicity"]<- "White"
                df[!is.na(df[,x]) & df[,x] %in% c("black", "black or african american", "african american"), 'ethnicity'] <- "Black or African American"
                df[is.na(df[,x]) | df[,x] %in% c("na", "unknown", "[discrepancy]", "discrepancy"), "ethnicity"] <- "NA"
                df[!df[,x] %in% c("Asian", "White", "Black or African American", "NA"), "ethnicity"] <- "Other"
            }
        } else {
            df[, x] <- "NA"
        }
    }
    res <- df
    return(res)
}

## from tt
t_n_m_process.func <- function(df, column_ids = c("t", "n", "m")) {
    for(x in column_ids) {
        id0 <- which(colnames(df) %in% c(x, paste0(x, "_stage"), paste0(x, " stage")))
        if(length(id0) >= 1) {
            colnames(df)[id0] <- x
            df[, x] <- toupper(df[,x])
            if(x == "t") {
                df[!is.na(df[, x]) & df[, x] %in%
                                c("T1","T1A","T1B","T1C","T1D","1","T1A1","T1A2","T1A3","T1A4",
                                  "T1B1","T1B2","T1B3","T1B4","T1C1","T1C2","T1C3","T1C4",
                                  "T1D1","T1D2","T1D3","T1D4"), x] <- 'T1'

                df[!is.na(df[,x]) & df[,x] %in%
                                    c("T2","T2A","T2B","T2C","T2D","2","T2A1","T2A2","T2A3","T2A4",
                                      "T2B1","T2B2","T2B3","T2B4","T2C1","T2C2","T2C3","T2C4",
                                      "T2D1","T2D2","T2D3","T2D4"), x] <- 'T2'

                df[!is.na(df[,x]) & df[,x] %in%
                                    c("T3","T3A","T3B","T3C","T3D","3","T3A1","T3A2","T3A3","T3A4",
                                      "T3B1","T3B2","T3B3","T3B4","T3C1","T3C2","T3C3","T3C4",
                                      "T3D1","T3D2","T3D3","T3D4"), x] <- 'T3'

                df[!is.na(df[,x]) & df[,x] %in%
                                    c("T4","T4A","T4B","T4C","T4D","4","T4A1","T4A2","T4A3","T4A4",
                                      "T4B1","T4B2","T4B3","T4B4","T4C1","T4C2","T4C3","T4C4",
                                      "T4D1","T4D2","T4D3","T4D4"), x] <- 'T4'
                df[is.na(df[,x]) | !df[,x] %in% c('T1','T2','T3','T4'), x] <- "NA"
            } else if(x == "n") {
                df[!is.na(df[,x]) & df[,x] %in%
                                    c("0","N0","N0A","N0B","N0C","N0D","N0A1","N0A2","N0A3","N0A4",
                                    "N0B1","N0B2","N0B3","N0B4","N0C1","N0C2","N0C3","N0C4",
                                    "N0D1","N0D2","N0D3","N0D4"), 'n'] <- 'N0'
                df[!is.na(df[,x]) & df[,x] %in%
                                    c("1","N1","N1A","N1B","N1C","N1D","N1A1","N1A2","N1A3","N1A4",
                                    "N1B1","N1B2","N1B3","N1B4","N1C1","N1C2","N1C3","N1C4",
                                    "N1D1","N1D2","N1D3","N1D4"), 'n'] <- 'N1'
                df[!is.na(df[,x]) & df[,x] %in%
                                    c("2","N2","N2A","N2B","N2C","N2D","N2A1","N2A2","N2A3","N2A4",
                                    "N2B1","N2B2","N2B3","N2B4","N2C1","N2C2","N2C3","N2C4",
                                    "N2D1","N2D2","N2D3","N2D4"), 'n'] <- 'N2'
                df[!is.na(df[,x]) & df[,x] %in%
                                    c("3","N3","N3A","N3B","N3C","N3D","N3A1","N3A2","N3A3","N3A4",
                                    "N3B1","N3B2","N3B3","N3B4","N3C1","N3C2","N3C3","N3C4",
                                    "N3D1","N3D2","N3D3","N3D4"), 'n'] <- 'N3'
                df[is.na(df[,x]) | !df[,x] %in% c('N0','N1','N2','N3'), 'n'] <- "NA"

            } else if(x == "m") {
                # m process
                df[!is.na(df[,x]) & df[,x] %in%
                                    c("0","M0","M0A","M0B","M0C","M0D","M0A1","M0A2","M0A3","M0A4",
                                        "M0B1","M0B2","M0B3","M0B4","M0C1","M0C2","M0C3","M0C4",
                                        "M0D1","M0D2","M0D3","M0D4"),'m']<-'M0'
                df[!is.na(df[,x]) & df[,x] %in%
                                    c("1","M1","M1A","M1B","M1C","M1D","M1A1","M1A2","M1A3","M1A4",
                                        "M1B1","M1B2","M1B3","M1B4","M1C1","M1C2","M1C3","M1C4",
                                        "M1D1","M1D2","M1D3","M1D4"),'m']<-'M1'
                df[is.na(df[,x]) | !df[,x] %in% c('M0','M1'),'m'] <- "NA"
            }
        } else {
            df[, x] <- "NA"
        }
    }
    res <- df
    return(res)
}

## from tt
stage_process.func <- function(df, column_ids = c("stage_edition", "stage")) {
    for(x in column_ids) {
        if(x == "stage_edition") { id0 <- which(colnames(df) %in% c('stage_edition', 'edition', 'stage_version', 'version')) }
        if(x == "stage") { id0 <- which(colnames(df) %in% c('stage')) }

        if(length(id0) >= 1) {
            colnames(cp_df)[id0] <- x

            if(x == "stage_edition") {
                df[, x] <- tolower(df[,x])
                df[!is.na(df[,x]) & df[,x] %in%
                                    c("ajcc_1st","ajcc.1st","ajcc 1st","ajcc 1","ajcc1st",'ajcc1','1st','1','one'),
                                    'stage_edition'] <- 'AJCC 1st'
                df[!is.na(df[,x]) & df[,x] %in% 
                                    c("ajcc_2nd","ajcc.2nd","ajcc 2nd","ajcc 2","ajcc2nd",'ajcc2','2nd','2','two'),
                                    'stage_edition']<-'AJCC 2nd'
                df[!is.na(df[,x]) & df[,x] %in%
                                    c("ajcc_3rd","ajcc.3rd","ajcc 3rd","ajcc 3","ajcc3rd",'ajcc3','3rd','3','three'),
                                    'stage_edition']<-'AJCC 3rd'
                df[!is.na(df[,x]) & df[,x] %in%
                                    c("ajcc_4th","ajcc.4th","ajcc 4th","ajcc 4","ajcc4th",'ajcc4','4th','4','four'),
                                    'stage_edition']<-'AJCC 4th'
                df[!is.na(df[,x]) & df[,x] %in%
                                    c("ajcc_5th","ajcc.5th","ajcc 5th","ajcc 5","ajcc5th",'ajcc5','5th','5','five'),
                                    'stage_edition']<-'AJCC 5th'
                df[!is.na(df[,x]) & df[,x] %in%
                                    c("ajcc_6th","ajcc.6th","ajcc 6th","ajcc 6","ajcc6th",'ajcc6','6th','6','six'),
                                    'stage_edition']<- 'AJCC 6th'
                df[!is.na(df[,x]) & df[,x] %in%
                                    c("ajcc_7th","ajcc.7th","ajcc 7th","ajcc 7","ajcc7th",'ajcc7','7th','7','seven'),
                                    'stage_edition']<-'AJCC 7th'
                df[!is.na(df[,x]) & df[,x] %in%
                                    c("ajcc_8th","ajcc.8th","ajcc 8th","ajcc 8","ajcc8th",'ajcc8','8th','8','eight'),
                                    'stage_edition']<-'AJCC 8th'
                df[!is.na(df[,x]) & df[,x] %in% 
                                    c("ajcc_9th","ajcc.9th","ajcc 9th","ajcc 9","ajcc9th",'ajcc9','9th','9','nine'),
                                    'stage_edition']<-'AJCC 9th'
                df[!is.na(df[,x]) & df[,x] %in%
                                    c("dukes stage","dukes","dukes_stage",'dukes.stage','duke'),
                                    'stage_edition']<-'Dukes'
                df[is.na(df[,x]) | df[,x] %in%
                                    c("na","unknown","Unknown","[Discrepancy]","Discrepancy"),
                                    'stage_edition'] <- "NA"
                df[!df[,x] %in%
                                    c('AJCC 1st','AJCC 2nd','AJCC 3rd','AJCC 4th','AJCC 5th','AJCC 6th','AJCC 7th','AJCC 8th','AJCC 9th','Dukes','NA'),
                                    'stage_edition'] <- "Other"
            } else if(x == "stage") {
                # stage process
                df[,x] <- toupper(df[,x])
                df[!is.na(df[,x]) & df[,x] %in% c("IA","IB","IC",'1A','1B','1C','1'), 'stage']<-'I'
                df[!is.na(df[,x]) & df[,x] %in% c("IIA","IIB","IIC",'2A','2B','2C','2'), 'stage']<-'II'
                df[!is.na(df[,x]) & df[,x] %in% c("IIIA","IIIB","IIIC",'3A','3B','3C','3'), 'stage']<-'III'
                df[!is.na(df[,x]) & df[,x] %in% c("IVA","IVB","IVC",'4A','4B','4C','4'), 'stage']<-'IV'
                df[is.na(df[,x]) | !df[,x] %in% c('I','II','III','IV'), 'stage'] <- "NA"
            }
        } else {
            df[, x] <- "NA"
        }
    }
    res <- df
    return(res)
}

## from tt
endpoint_process.func <- function(df, column_ids = c("os", "rfs", "dfs", "dss", "pfs", "ttp")) {
    # df = cp_df; column_ids = c("os", "rfs", "dfs", "dss", "pfs", "ttp")
    # endpoint process ----
    endpoint_events <- c()
    median_followups <- c()
    zero_status <- 0
    zero_endpoint_status <- c()
    missing_endpoint <- 0
    
    df0 <- df[, c("os_months","os_status","rfs_months","rfs_status",
       "dfs_months","dfs_status", "dss_months","dss_status",
       "pfs_months","pfs_status","ttp_months","ttp_status")]
    debug(df0)


    for(x in column_ids) {
        # x="os"
        id0 <- which(colnames(df) %in% c(paste0(x, "_days"), paste0(x, "_months")))
        id1 <- which(colnames(df) %in% paste0(x, "_status"))

        if(length(id0) == 1){
            if(colnames(df)[id0] == paste0(x, "_days")) {
                colnames(df)[id0] <- paste0(x, "_months")
                df[, paste0(x, "_months")] <- as.numeric(df[, paste0(x, "_months")])/30.41667
            } else {
                colnames(df)[id0] <- paste0(x, "_months")
            }
            
            if(length(id1) == 1) {
                mode(df[, paste0(x, "_months")]) <- "numeric"
                mode(df[, paste0(x, "_status")]) <- "numeric"        
                df_f <- df[!is.na(df[, paste0(x, "_months")]) & !is.na(df[, paste0(x, "_status")]), ]

                if(nrow(df_f) > 0) {
                    evt <- sum(df_f[, paste0(x, "_status")][!is.na(df_f[, paste0(x, "_status")])] == 1)
                    endpoint_event <- c()
                    if(evt > 0) { endpoint_event <- paste0("<b>", toupper(x), "</b>: ", nrow(df_f), " (", evt, " events)") }
                    
                    if(length(unique(df_f[, paste0(x, "_status")])) == 1) {
                        if(unique(df_f[, paste0(x, "_status")]) == 1) {
                            med_f <- median(df_f[, paste0(x, "_months")])
                            med_f <- paste0("<b>", toupper(x), "</b>: ", med_f)
                        } else if(unique(df_f[, paste0(x, "_status")]) == 0) {
                            zero_status <- zero_status + 1
                            zero_endpoint_status <- c(zero_endpoint_status, paste0(toupper(x), "_status"))
                            med_f <- c()
                        } else {
                            med_f <- c()
                        }
                    } else {
                        df_f0 <- df_f
                        colnames(df_f0)[which(colnames(df_f0) == paste0(x, "_months"))] <- "time"
                        colnames(df_f0)[which(colnames(df_f0) == paste0(x, "_status"))] <- "status"
                        fit.null <- surv_fit(Surv(time = time, status == 0) ~ 1, data = df_f0)
                        med_f <- round(surv_median(fit.null)$median, 2)
                        med_f <- paste0("<b>", toupper(x), "</b>: ", med_f)
                    }
                } else {
                    med_f <- c()
                    endpoint_event <- c()
                    missing_endpoint <- missing_endpoint + 1
                }
                median_followups <- c(median_followups, med_f)
                endpoint_events <- c(endpoint_events, endpoint_event)
            }
        } else {
            df[, paste0(x, "_months")] <- NA
            df[, paste0(x, "_status")] <- NA
            missing_endpoint <- missing_endpoint + 1
        }
    }
    res <- list(endpoint_events,  median_followups,  zero_status, zero_endpoint_status, missing_endpoint)
    names(res) <- c("endpoint_events",  "median_followups",  "zero_status", "zero_endpoint_status", "missing_endpoint")
    res
}

# Association/correlation/model --------------------------------------------------------------
# glmnet lbs_fun
lbs_fun <- function(fit,coef_non0_1,coef_non0_2,lab.cex=1,...) {
  L <- length(fit$lambda)
  x1 <- log(fit$lambda[L])
  x2 <- log(fit$lambda[1])
  
  y <- fit$beta[, L]
  y1<-y[coef_non0_1$symbol]
  y2<-y[coef_non0_2$symbol]
  
  labs1 <- names(y1)
  labs2 <- names(y2)
  
  text(x1-1, y1,cex=lab.cex, labels=labs1,adj =0, ...)
  text(x1-2.3, y2,cex=lab.cex, labels=labs2,adj =0, ...)
}

## LASSO penalized Cox model ---
lasso_cox.func <- function(x, y, alpha = 1, family='cox', nfolds=10, parallel = TRUE, cores = 8) {
    doMC::registerDoMC(cores = cores)
    set.seed(2019)
    cv_lasso <- glmnet::cv.glmnet(x = as.matrix(x), y = y, alpha = alpha,
                                  family = family, nfolds = nfolds, 
                                  parallel=parallel)
    minlambda <- cv_lasso$lambda.min
    out <- list(cv_lasso, minlambda)
    return(out)
}

# Fit SPCA model to get important genes from candidate genes ---
spca.msig.func <- function(train.exp, sample.survival, unit=c('day','month'), outprex = "spca", n.folds = 5) {
    require(superpc)
    mexp0 <- train.exp
    surv0 <- sample.survival
    colnames(surv0)[1:3] <- c("SampleID", "time", "status")
    ##
    surv0 <- surv0[!is.na(surv0$status) & !is.na(surv0$time),]
    #rownames(surv0) <- surv0$SampleID
    mode(surv0$time) <- 'numeric'
    mode(surv0$status) <- 'numeric'
    if(unit == 'month'){ surv0$time <- surv0$time * 30 }

    ##
    sams = intersect(colnames(mexp0), surv0$SampleID)
    surv0 = surv0[surv0$SampleID %in% sams, ]
    x = mexp0[, sams]
    y = surv0$time
    censoring.status = surv0$status
    featurenames = rownames(x)
    data <- list(x = x, y = y, censoring.status = censoring.status, featurenames = featurenames)
    data.test <- data

    # train the model. This step just computes the scores for each feature
    train.obj <- superpc.train(data, type = "survival")

    # note for regression (non-survival) data, we leave the component "censoring.status"
    # out of the data object, and call superpc.train with type="regression".
    # otherwise the superpc commands are all the same

    # cross-validate the model
    #set.seed(2021)
    #cv.obj <- superpc.cv(train.obj, data, n.fold = n.folds)

    # plot the cross-validation curves. From this plot we see that the 1st
    #pdf(file = paste0(outprex, '.plotcv.pdf'), width = 4, height = 4)
    #superpc.plotcv(cv.obj)
    #dev.off()

    # See pdf version of the cv plot
    # here we have the luxury of test data, so we can compute the likelihood ratio statistic
    # over the test data and plot them. We see that the threshold of xx works pretty well
    set.seed(2021)
    lrtest.obj <- superpc.lrtest.curv(train.obj, data, data.test)
    lr <- lrtest.obj$lrtest
    thresholds <- lrtest.obj$threshold
    num_feas <- lrtest.obj$num.features

    ##
    th0 <- thresholds[which(lr == max(lr))[1]]
    id0 <- which(thresholds == th0)

    if(lrtest.obj$num.features[id0] < 5) {
        ids = which(num_feas >= 5)
        lr = lr[ids]
        thresholds <- thresholds[ids]
        th0 <- thresholds[which(lr == max(lr))[1]]
    }
    print(th0)

    #pdf(file = paste0(outprex, '.plotlrtest.pdf'), width = 4, height = 4)
    #superpc.plot.lrtest(lrtest.obj)
    #dev.off()

    # See pdf version of the lrtest plot
    # now we derive the predictor of survival for the test data,
    # and then then use it
    # as the predictor in a Cox model. We see that the 1st supervised PC is
    # highly significant; the next two are not
    #set.seed(2021)
    #fit.cts <- superpc.predict(train.obj, data, data.test, threshold = th0,
                               #n.components = 3, prediction.type = "continuous")
    #superpc.fit.to.outcome(train.obj, data.test, fit.cts$v.pred)

    # sometimes a discrete (categorical) predictor is attractive.
    # Here we form two groups by cutting the predictor at its median
    # and then plot Kaplan-Meier curves for the two groups
    #fit.groups <- superpc.predict(train.obj, data, data.test, threshold = th0,
                                   #n.components = 1, prediction.type="discrete")
    #superpc.fit.to.outcome(train.obj, data.test, fit.groups$v.pred)

    #pdf(file = paste0(outprex, '.KM.median.pdf'), width = 6, height = 6)
    #plot(survfit(Surv(data.test$y, data.test$censoring.status)~fit.groups$v.pred), col = c('royalblue', 'red'), xlab="time", ylab="Prob survival", lwd =3)
    #legend("topright", legend = c("High SPCA score", "Low SPCA score"), lwd =2, col = c('red', "royalblue"), bty = 'n')
    #dev.off()

    #See pdf version of the survival plot
    # Finally, we look for a predictor of survival a small number of
    # genes (rather than all 1000 genes). We do this by computing an importance
    # score for each equal its correlation with the supervised PC predictor.
    # Then we soft threshold the importance scores, and use the shrunken
    # scores as gene weights to from a reduced predictor. Cross-validation
    # gives us an estimate of the best amount to shrink and an idea of
    # how well the shrunken predictor works.
    set.seed(2021)
    fit.red <- superpc.predict.red(train.obj, data, data.test, threshold = th0)
    #fit.redcv <- superpc.predict.red.cv(fit.red, cv.obj, data, threshold = th0)
    #pdf(file= paste0(outprex, '.plotred.lrtest.pdf'), width = 6, height = 6)
    #superpc.plotred.lrtest(fit.redcv)
    #dev.off()

    #See pdf version of this plot
    # Finally we list the significant genes, in order of decreasing importance score
    imp_feas <- superpc.listfeatures(data, train.obj, fit.red,
                        #fitred.cv = fit.redcv,
                        component.number=1,
                        #num.features=NULL
                        )
    imp_feas <- as.data.frame(imp_feas)
    return(imp_feas)
}

# survival related ----------------------------------------------
survival_1VS2_diffcal.fun <- function(sample.survival, pts.cluster) {
  require(survival)
  require(survminer)
  # get survival data for selected patients
  colnames(sample.survival) <- c("time", "status")
  overlapname<-intersect(rownames(sample.survival), names(pts.cluster))
  x <- data.frame(sample.survival[overlapname,], cluster=pts.cluster[overlapname])

  # survival analysis
  mode(x[,"time"]) <- 'numeric'
  mode(x[,"status"]) <- 'numeric'
  fit<-survfit(Surv(time, status)~cluster, data=x)
  surv_diff<-survdiff(Surv(time, status)~cluster, data=x)
  
  p.val <- (1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1))
  p <-ifelse(p.val >=0.01,format(p.val ,scientific=F,digits=3),format(p.val ,scientific=TRUE,digits=2))
  
  HR1 = round((surv_diff$obs[1]/surv_diff$exp[1])/(surv_diff$obs[2]/surv_diff$exp[2]), 2)
  HR=paste0("HR = ",HR1)
  up95 = round(exp(log(HR1) + qnorm(0.975)*sqrt(1/surv_diff$exp[2]+1/surv_diff$exp[1])),2)
  low95 = round(exp(log(HR1) - qnorm(0.975)*sqrt(1/surv_diff$exp[2]+1/surv_diff$exp[1])),2)
  CI<-paste(round(low95, 2), round(up95, 2), sep='-')
  group<-as.vector(surv_diff$n)
  #tempres<-data.frame(group=paste(group,collapse=";"), p=p,'HR'= HR1,check.names=F)
  tempres<-data.frame(t(data.frame(group)), p=p,'HR'= HR1,check.names=F)
  colnames(tempres)[1:length(group)] <- paste0("number", 1:length(group))
  return(tempres)
}


##
survival_1VS2_plot_multi3 <- function(sample.survival, pts.cluster, control_group = NULL,
                                      plot.sigOnly =TRUE, imageName="KM",
                                      palette = c("#d42a24", "#00427d"),                                      
                                      pval.coord = c(max((sample.survival[,1]))*0.1,0.2),
                                      pval.size = 3.5,
                                      risk.table.fontsize = 3.5,
                                      legend.title = "Group",
                                      legend.labs = c('High','Low'),
                                      #legend = c(0.7, 0.9),
                                      title = "", xlab = "Time (months)",
                                      ylab = "Probability of survival",
                                      
                                      risk.table = "abs_pct",
                                      #risk.table.col = "strata",
                                      risk.table.title = "Number at risk: n (%)",
                                      ncensor.plot=FALSE,
                                      
                                      surv.plot.height = 0.75, risk.table.height = 0.25,
                                      cumevents=FALSE, cumevents.title = "Cumulative number of deaths",
                                      font.legend = c(10, "plain", "black"),
                                      font.title = c(10, "bold", "black"), 
                                      font.x = c(10, "plain", "black"),
                                      font.y = c(10, "plain", "black"), 
                                      font.tickslab = c(10, "plain", "black"),
                                      font.table.main = 10, 
                                      #size=c(0.7,0.7,0.7,0.7),
                                      surv.scale = "percent",
                                      newpage = FALSE) {
  require(survival)
  require(survminer)
  # get survival data for selected patients
  colnames(sample.survival) = c("time", "status")
  overlapname = intersect(rownames(sample.survival), names(pts.cluster))
  x = data.frame(sample.survival[overlapname,], cluster=pts.cluster[overlapname])

  # survival analysis
  mode(x[,"time"]) = 'numeric'
  mode(x[,"status"]) = 'numeric'
  fit = survfit(Surv(time, status)~cluster, data=x)
  surv_diff = survdiff(Surv(time, status)~cluster, data=x)
  group = as.vector(surv_diff$n)
  p.val = (1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1))

  ##
  if(!is.null(control_group)) {
      levs = sort(unique(x$cluster))
      levs = levs[levs!=control_group]
      in.list = lapply(levs, function(lev0){
          # lev0 = levs[1]
          x %>% filter(cluster==control_group | cluster==lev0)
      })
      names(in.list) <- levs
  } else { in.list = list(all = x) }

  ##
  p_text.list = lapply(names(in.list), function(type0) {
    # type0 = names(in.list)[1]
    x = in.list[[type0]]
    surv_diff = survdiff(Surv(time, status)~cluster, data=x)
    ##
    p.val = (1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1))
    p = ifelse(p.val >=0.01,format(p.val ,scientific=F,digits=3),format(p.val ,scientific=TRUE,digits=2))
    
    ##
    HR1 = round((surv_diff$obs[1]/surv_diff$exp[1])/(surv_diff$obs[2]/surv_diff$exp[2]), 2)
    HR = paste0("HR = ",HR1)
    up95 = round(exp(log(HR1) + qnorm(0.975)*sqrt(1/surv_diff$exp[2]+1/surv_diff$exp[1])),2)
    low95 = round(exp(log(HR1) - qnorm(0.975)*sqrt(1/surv_diff$exp[2]+1/surv_diff$exp[1])),2)
    CI = paste(round(low95, 2), round(up95, 2), sep='-')
    ##
    if(type0 == "all") {
        p_text = paste0('P = ',p,'; ',paste0(HR,' (',CI,')'))
    } else {
        p_text = paste0(type0, ': P = ',p,'; ',paste0(HR,' (',CI,')'))
    }
    p_text
  })

  p_text = do.call(c,p_text.list)
  p_text = paste(p_text, collapse = "\n")

  # survival plot 
  if(!plot.sigOnly | (min(group) >= 5 & p.val <=0.05)) {
    theme <- theme(axis.line = element_line(colour = "black"),
                   panel.border = element_blank(),
                   panel.background = element_blank())
    ##
    survp <- ggsurvplot(fit, data = x,
                        pval = p_text, pval.size=pval.size, pval.coord=pval.coord,
                        break.x.by=NULL, xscale=30,
                        ##
                        font.legend = font.legend,
                        font.title = font.title, 
                        font.x = font.x,
                        font.y = font.y,
                        font.tickslab = font.tickslab,

                        ##
                        main = "", title=title,
                        xlab=xlab, ylab=ylab,
                        #legend = legend,
                        legend.title = legend.title,
                        legend.labs = legend.labs,

                        # Specific to the risk table
                        risk.table = risk.table, # Add risk table
                        #risk.table.col = risk.table.col, # Change risk table color by groups
                        risk.table.title = risk.table.title,
                        risk.table.pos = "out",
                        risk.table.fontsize = risk.table.fontsize,
                        risk.table.y.text.col = TRUE,
                        risk.table.y.text = FALSE, 
                        fontsize  = pval.size,
                        #ncensor.plot=ncensor.plot,

                        risk.table.height = risk.table.height,
                        surv.plot.height = surv.plot.height,
                        cumevents=cumevents, cumevents.title = cumevents.title,
                        ##
                        #surv.median.line = "hv", # Specify median survival
                        ggtheme = theme,
                        palette = palette,
                        #size = size,
                        #tables.theme = theme,
                        tables.theme = theme_cleantable(),
                        surv.scale = surv.scale,
                        axes.offset = TRUE)
                        
    survp$table <- survp$table + theme(plot.title = element_text(size=font.table.main))
  }
  return(survp)
}


# univariate cox regression for preditive model ---
uni_cox.preditive <- function(pts_surv, clinic_ann, treatment_columnID = NULL) {
    require(survival)
    require(plyr)
    if(is.null(treatment_columnID)) { stop("treatment columnID is not specified.") }
    ##
    train_df <- merge.data.frame(x=clinic_ann, y=pts_surv,by='SampleID',all.x=T)
    rownames(train_df) <- train_df[,1]
    train_df <- train_df[!is.na(train_df$time) & train_df$time>0,]
    
    ## univariate cox
    variables <- colnames(train_df)[which(!colnames(train_df) %in% c('time','status', treatment_columnID))[-1]]
    variables <- unlist(lapply(variables, function(x) {if(n_distinct(train_df[,x]) > 1) x }))

    uni_cox_out <- lapply(variables, function(x) {
        if(n_distinct(train_df[,x]) > 1 & n_distinct(train_df[,x]) < 6) {
            df <- train_df[,c("time", "status", x, treatment_columnID)]
            df <- df[!is.na(df[,x]), ]
            colnames(df)[(which(colnames(df) %in% x))] <- "x"
            colnames(df)[(which(colnames(df) %in% treatment_columnID))] <- "treatment"
            df$x <- as.factor(df$x)
            df$treatment <- as.factor(df$treatment)

            ## subgroup analysis
            levs <- levels(df$x)
            cox0 <- do.call(rbind.data.frame, lapply(levs, function(lev0) {
                # levs[1]->lev0
                train_df0 <- df %>% filter(x == lev0) %>% select(-x)
                train_df0[train_df0 == "NA" | train_df0 == "Unknown" | train_df0 == "[Discrepancy]"] <- NA

                if(n_distinct(train_df0$treatment, na.rm = T) >= 2 & all(table(train_df0$treatment) >= 2)) {
                    cox0 <- summary(survival::coxph(Surv(time, status)~treatment, data=train_df0, na.action = na.omit))
                    res0 <- cox0$coefficients
                    if(nrow(res0) == 1) {
                        ci <- t(data.frame(c(lower95 = round(cox0$conf.int[,"lower .95"], 3),
                                                upper95 = round(cox0$conf.int[,"upper .95"], 3),
                                                CI95 = paste0(round(res0[,"exp(coef)"],3),
                                                " (", round(cox0$conf.int[,"lower .95"], 3),
                                                "-", round(cox0$conf.int[,"upper .95"], 3), ")"))))
                    } else {
                        ci <- data.frame(lower95 = round(cox0$conf.int[,"lower .95"], 3),
                                            upper95 = round(cox0$conf.int[,"upper .95"], 3),
                                            CI95 = paste0(round(res0[,"exp(coef)"], 3),
                                            " (", round(cox0$conf.int[,"lower .95"], 3),
                                            "-", round(cox0$conf.int[,"upper .95"], 3), ")"),
                                            check.names=FALSE)
                    } 
                    res <- data.frame(res0, ci, check.names=F)
                    ##
                    N = cox0$n
                    n = cox0$nevent
                    res$'Events/N' <- paste0(n, "/", N)
                    res$'exp(coef)' <- round(res$'exp(coef)', 3)
                    res$'Pr(>|z|)' <- ifelse(res$'Pr(>|z|)' >= 0.01,
                                                format(res$'Pr(>|z|)', digits=3),
                                                format(res$'Pr(>|z|)', scientific = T, digits=2))
                } else {
                    res <- as.data.frame(matrix(NA, nrow=1, ncol=9))
                    colnames(res) <- c("coef", "exp(coef)", "se(coef)", "z", "Pr(>|z|)", "lower95", "upper95", "CI95", "Events/N")
                    N = nrow(train_df0)
                    n = nrow(subset(train_df0, status == 1))
                    res$'Events/N' <- paste0(n, "/", N)
                }
                res
            }))
            cox0 <- data.frame(Variable = levs, cox0, check.names=F, row.names=2:(nrow(cox0)+1))
            treat_n <- table(df[,c("x", "treatment")])
            treat_id <- paste(colnames(treat_n), collapse = ",")
            treat_n <- data.frame(paste(treat_n[,1], treat_n[,2],sep = ", "))
            names(treat_n) <- treat_id
            cox0$Variable <- paste0(paste(cox0$Variable, treat_n[,1], sep = " ("), ")")

            ## interaction analysis
            cox1 <- summary(survival::coxph(survival::Surv(time,status==1) ~ x*treatment, data=df), na.action=na.omit)
            cox1 <- data.frame(cox1$coefficients, check.names=F)
            p = cox1$'Pr(>|z|)'

            ## combine
            cox0 <- data.frame(cox0[, c("Variable", "Events/N", "exp(coef)", "lower95", "upper95", "CI95", "Pr(>|z|)")],
                               P_interaction = "", check.names=F)
            p <- ifelse(p[length(p)] >= 0.01,
                        format(p[length(p)], digits=3),
                        format(p[length(p)], scientific = T, digits=2))

            cox1 <- data.frame(Variable = x, t(data.frame(rep("", ncol(cox0)-2))),
                                P_interaction = p, row.names =1, check.names=F)
            colnames(cox1) <- colnames(cox0)
            uni_cox <- rbind(cox1, cox0)
            colnames(uni_cox) <- c("Variable", "Events/N", "HR", "lower95", "upper95", "HR for events (95% CI)", "P", "P_interaction")
            uni_cox <- data.frame(Type = x, uni_cox, check.names=F)
        } else {
            uni_cox <- data.frame()
        }
        uni_cox
    })
    res <- do.call(rbind.data.frame, uni_cox_out)
    res[is.na(res)] <- ""
    res
}

# multivariate cox regression for preditive model ---
multi_cox.preditive <- function(pts_surv, clinic_ann, treatment_columnID = NULL) {
    require(survival)
    require(plyr)
    if(is.null(treatment_columnID)) { stop("treatment columnID is not specified.") }
    ##
    train_df <- merge.data.frame(x=clinic_ann, y=pts_surv,by='SampleID',all.x=T)
    rownames(train_df) <- train_df[,1]
    train_df <- train_df[!is.na(train_df$time) & train_df$time>0,]

    ## multivariate cox
    variables <- colnames(train_df)[which(!colnames(train_df) %in% c('time','status', treatment_columnID))[-1]]
    variables <- unlist(lapply(variables, function(x) {if(n_distinct(train_df[!is.na(train_df[,x]),x]) > 1) x }))
    ##
    variables_exclude <- unlist(lapply(variables, function(x) {if((nrow(train_df[is.na(train_df[,x]),])/nrow(train_df)) > 0.5) x }))
    variables <- variables[!variables %in% variables_exclude]
    train_df <- train_df[, !colnames(train_df) %in% variables_exclude]

    multi_cox_out <- lapply(variables, function(x) {
        if(n_distinct(train_df[,x]) > 1 & n_distinct(train_df[,x]) < 6) {
            df <- train_df
            df <- df[!is.na(df[,x]), ]
            colnames(df)[(which(colnames(df) %in% x))] <- "x"
            colnames(df)[(which(colnames(df) %in% treatment_columnID))] <- "treatment"
            df$x <- as.factor(df$x)
            df$treatment <- as.factor(df$treatment)

            ## subgroup analysis
            levs <- levels(df$x)
            cox0 <- do.call(rbind.data.frame, lapply(levs, function(lev0) {
                train_df0 <- df[df$x == lev0, !colnames(df) %in% c("SampleID", "x")]
                train_df0 <- train_df0[, unique(c("treatment", colnames(train_df0)))]
                train_df0[train_df0 == "NA" | train_df0 == "Unknown" | train_df0 == "[Discrepancy]"] <- NA
                ids <- colnames(train_df0)[!colnames(train_df0) %in% c("treatment", "time", "status")]
                df0 <- data.frame(train_df0[, ids])
                colnames(df0) <- ids
                id0 <- which(apply(df0, 2, function(x) n_distinct(x[!is.na(x)])) > 1)
                train_df0 <- train_df0[, c("treatment", "time", "status", names(id0))]
                ##
                res <- as.data.frame(matrix(NA, nrow=1, ncol=9))
                colnames(res) <- c("coef", "exp(coef)", "se(coef)", "z", "Pr(>|z|)", "lower95", "upper95", "CI95", "Events/N")
                N = nrow(train_df0)
                n = nrow(subset(train_df0, status == 1))
                res$'Events/N' <- paste0(n, "/", N)
                ##                
                if(n_distinct(train_df0$treatment, na.rm = T) >= 2 & all(table(train_df0$treatment) >= 2)) {
                    ## remove category variables with levels less than 3 values
                    ids <- do.call(c, lapply(colnames(train_df0), function(x) {
                        if(n_distinct(train_df0[,x], na.rm = T) > 1 & n_distinct(train_df0[,x], na.rm = T) < 6) {
                            if(all(table(train_df0[,x]) >= 3)) { x }
                        } else { x }
                    }))
                    train_df0 <- train_df0[,ids]
                    ##
                    if(all(c("time", "status") %in% colnames(train_df0))) {
                        cox0 <- summary(coxph(Surv(time, status==1)~., data = train_df0, na.action = na.omit))
                        ci <- data.frame(lower95 = round(cox0$conf.int[,"lower .95"], 3),
                                         upper95 = round(cox0$conf.int[,"upper .95"], 3),
                                         CI95 = paste0(round(cox0$coefficients[,"exp(coef)"], 3),
                                                     " (", round(cox0$conf.int[,"lower .95"], 3),
                                                     "-", round(cox0$conf.int[,"upper .95"], 3), ")"),
                                         check.names = FALSE)
                                    
                        res <- data.frame(cox0$coefficients, ci, check.names=F)
                        N = cox0$n
                        n = cox0$nevent
                        res$'Events/N' <- paste0(n, "/", N)
                        res$'exp(coef)' <- round(res$'exp(coef)', 3)
                        res$'Pr(>|z|)' <- ifelse(res$'Pr(>|z|)' >= 0.01,
                                                 round(res$'Pr(>|z|)', 3),
                                                 format(res$'Pr(>|z|)', scientific = T, digits=2))
                        res <- res[1,]
                    }
                }
                res
            }))

            cox0 <- data.frame(Variable = levs, cox0, row.names = 2:(nrow(cox0)+1), check.names=F)

            treat_n <- table(df[,c("x", "treatment")])
            treat_id <- paste(colnames(treat_n), collapse = ",")
            treat_n <- data.frame(paste(treat_n[,1], treat_n[,2],sep = ", "))
            names(treat_n) <- treat_id
            cox0$Variable <- paste0(paste(cox0$Variable, treat_n[,1], sep = " ("), ")")
            
            ## interaction analysis
            formula0 <- paste(c("x*treatment", variables[variables!=x]), collapse = " + ")
            formula0 <- as.formula(paste("Surv(time,status==1)  ~ ",formula0, sep = ""))
            cox1 <- summary(coxph(formula0 , data=df), na.action = na.omit)
            cox1 <- data.frame(cox1$coefficients, check.names=F)
            cox1 <- cox1[nrow(cox1), ]
            p <- cox1$'Pr(>|z|)'

            ## combine
            cox0 <- data.frame(cox0[, c("Variable", "Events/N", "exp(coef)", "lower95", "upper95", "CI95", "Pr(>|z|)")],
                                P_interaction = "", check.names=F)
            p <- ifelse(p[length(p)] >= 0.01,
                        format(p[length(p)], digits=3),
                        format(p[length(p)], scientific = T, digits=2))

            cox1 <- data.frame(Variable = x, t(data.frame(rep("", ncol(cox0)-2))),
                                P_interaction = p, row.names =1, check.names=F)
            colnames(cox1) <- colnames(cox0)
            multi_cox <- rbind(cox1, cox0)
            colnames(multi_cox) <- c("Variable", "Events/N", "HR", "lower95",
                                     "upper95", "HR for events (95% CI)", 
                                     "P", "P_interaction")
            multi_cox <- data.frame(Type = x, multi_cox, check.names=F)
        } else {
            multi_cox <- data.frame()
        }
        multi_cox
    })

    res <- do.call(rbind.data.frame, multi_cox_out)
    res[is.na(res)] <- ""
    res
}

# univariate cox regression for prognosis model ---
uni_cox.prognosis <- function(pts_surv, clinic_ann) {
    require(survival)
    require(plyr)
    train_df <- merge.data.frame(x=clinic_ann, y=pts_surv,by='SampleID',all.x=T)
    rownames(train_df) <- train_df[,1]
    train_df <- train_df[!is.na(train_df$time) & train_df$time>0,]

    ## univariate cox
    variables <- colnames(train_df)[which(!colnames(train_df) %in% c('time','status'))[-1]]
    if(length(variables)>=1) {
        uni_cox_out <- lapply(variables, function(x) {
            # variables[3] -> x
            if(length(unique(train_df[,x])) > 1) {
                    df <- train_df[,c("time", "status", x)]
                    df <- df[!is.na(df[,x]), ]
                    colnames(df)[(which(colnames(df) %in% x))] <- "x"

                    ## unicox analysis
                    train_df0 <- df
                    cox0 <- summary(survival::coxph(Surv(time, status)~x, data=train_df0, na.action = na.omit))

                    ##
                    res0 <- cox0$coefficients
                    if(nrow(res0) == 1) {
                        ci <- t(data.frame(c(lower95 = round(cox0$conf.int[,"lower .95"], 3),
                                             upper95 = round(cox0$conf.int[,"upper .95"], 3),
                                             CI95 = paste0(round(res0[,"exp(coef)"],3),
                                             " (", round(cox0$conf.int[,"lower .95"], 3),
                                             "-", round(cox0$conf.int[,"upper .95"],3), ")"))))
                    } else {
                        ci <- data.frame(lower95 = round(cox0$conf.int[,"lower .95"], 3),
                                         upper95 = round(cox0$conf.int[,"upper .95"], 3),
                                         CI95 = paste0(round(res0[,"exp(coef)"], 3),
                                         " (", round(cox0$conf.int[,"lower .95"], 3),
                                         "-", round(cox0$conf.int[,"upper .95"], 3), ")"),
                                         check.names=FALSE)
                    }

                    res <- data.frame(res0, ci, check.names=F)
                    N = cox0$n
                    n = cox0$nevent
                    res$'Events/N' <- paste0(n, "/", N)
                    res$coef <- round(res$coef, 3)
                    res$'exp(coef)' <- round(res$'exp(coef)', 3)
                    res$'Pr(>|z|)' <- ifelse(res$'Pr(>|z|)' >= 0.01,
                                            format(res$'Pr(>|z|)', digits=3),
                                            format(res$'Pr(>|z|)', scientific = T, digits=2))
                    cox0 <- res
                    if(is.numeric(train_df[,x])) {
                        cox0 <- data.frame(Variable = x, cox0, check.names=F, row.names=2:(nrow(cox0)+1))
                    } else {
                        cox0 <- data.frame(Variable = gsub("^x", "",rownames(res)), cox0, check.names=F, row.names=2:(nrow(cox0)+1))
                    }
                    ## 
                    if(is.character(train_df[,x])) {
                        ref_lev = paste0(" (vs. ", sort(train_df[,x])[1], ")")
                    } else if(is.factor(train_df[,x])) {
                        ref_lev = paste0(" (vs. ", levels(train_df[,x])[1], ")")
                    } else {
                        ref_lev = ""
                    }

                    ## 
                    cox0 <- data.frame(cox0[, c("Variable", "coef",  "Events/N", "exp(coef)", "lower95", "upper95", "CI95", "Pr(>|z|)")], check.names=F)

                    if(is.numeric(train_df[,x])) {
                       uni_cox <- cox0
                    } else {
                        cox1 <- data.frame(Variable = paste0(x, ref_lev),
                            t(data.frame(rep("", ncol(cox0)-1))),
                            row.names = 1, check.names = F)
                        colnames(cox1) <- colnames(cox0)
                        uni_cox <- rbind(cox1, cox0)
                    }
                    colnames(uni_cox) <- c("Variable", "Coefficient", "Events/N", "HR", "lower95", "upper95", "HR for events (95% CI)", "P")

                    # recalculate events/N for character/factor variables
                    if(is.character(train_df[,x]) | is.factor(train_df[,x])) {
                        if(is.character(train_df[,x])) {
                            levs <- sort(unique(train_df[,x]))
                        } else if(is.factor(train_df[,x])) {
                            levs <- levels(train_df[,x])
                        }                        
                        for(lev0 in levs) {
                            N = nrow(train_df[!is.na(train_df[,x]) & train_df[,x]==lev0,])
                            n = nrow(train_df[(!is.na(train_df[,x]) & train_df[,x]==lev0) &
                                            (!is.na(train_df$status) & train_df[,"status"]==1),])
                            if(lev0 == levs[1]) {
                                uni_cox[1, "Events/N"] <- paste0(n, "/", N)
                            } else {
                                uni_cox[uni_cox$Variable == lev0, "Events/N"] <- paste0(n, "/", N)
                            }
                        }
                    }
                    uni_cox <- data.frame(Type = x, uni_cox, check.names=F)
                    uni_cox
            }
        })
        res <- do.call(rbind.data.frame, uni_cox_out)
    } else {
        res <- data.frame()
    }
    res
}

# multivariate cox regression for prognosis model ---
multi_cox.prognosis <- function(pts_surv, clinic_ann) {
    require(survival); require(plyr)
    train_df <- merge.data.frame(x=clinic_ann, y=pts_surv,by='SampleID',all.x = TRUE)
    rownames(train_df) <- train_df[,1]
    train_df <- train_df[!is.na(train_df$time) & train_df$time>0,]

    ## multivariate cox
    train_df <- train_df[,which(apply(train_df, 2, function(x) {length(unique(x))}) >=2)]
    variables <- colnames(train_df)[which(!colnames(train_df) %in% c('time','status', "SampleID"))]
    train_df <- train_df[,c('time','status',variables)]
    if(length(variables)>=2) {
        multi_cox_out <- lapply(variables, function(x) {
            # variables[6] -> x
            train_df0 <- train_df[, unique(c('time','status', x, variables))]
            N = nrow(train_df0[!is.na(train_df0[,x]),])
            n = nrow(train_df0[!is.na(train_df0[,x]) & train_df0$status==1,])
            
            ## multiicox analysis
            colnames(train_df0)[-(1:2)] <- paste0(colnames(train_df0)[-(1:2)],'=')
            cox0 <- summary(survival::coxph(Surv(time, status)~., data=train_df0, na.action = na.omit))
            ci <- data.frame(lower95 = round(cox0$conf.int[,"lower .95"], 3),
                             upper95 = round(cox0$conf.int[,"upper .95"], 3),
                             CI95 = paste0(round(cox0$coefficients[,"exp(coef)"],3),
                             " (", round(cox0$conf.int[,"lower .95"], 3),
                             "-", round(cox0$conf.int[,"upper .95"], 3), ")"), check.names=FALSE)
                    
            res <- data.frame(cox0$coefficients, ci, check.names=FALSE)
            
            res$coef <- round(res$coef, 3)
            res$'Events/N' <- paste0(n, "/", N)
            res$'exp(coef)' <- round(res$'exp(coef)', 3)
            res$'Pr(>|z|)' <- ifelse(res$'Pr(>|z|)' >= 0.01,
                                     round(res$'Pr(>|z|)', digits=3),
                                     format(res$'Pr(>|z|)', scientific = TRUE, digits=2))
            res <- res[grep(paste0("`", x, "="), rownames(res), fixed = TRUE),]
            cox0<- res

            if(!(is.character(train_df0[,paste0(x,'=')]) | is.factor(train_df0[,paste0(x,'=')]))) {
                levs <- gsub('`','', do.call(rbind, strsplit(rownames(res),'=`')))
            } else {
                levs <- gsub('`','', do.call(rbind, strsplit(rownames(res),'=`')))[,2]
            }
            cox0 <- data.frame(Variable = levs, cox0, check.names=F, row.names=2:(nrow(cox0)+1))

            ## 
            if(is.character(train_df[,x])) {
                ref_lev = paste0(" (vs. ", sort(train_df[,x])[1], ")")
            } else if(is.factor(train_df[,x])) {
                ref_lev = paste0(" (vs. ", levels(train_df[,x])[1], ")")
            } else { ref_lev = "" }

            ## 
            cox0 <- data.frame(cox0[, c("Variable", "coef", "Events/N", "exp(coef)", "lower95", "upper95", "CI95", "Pr(>|z|)")], check.names=F)
            if(is.numeric(train_df[,x])) {
                multi_cox <- cox0
            } else {
                cox1 <- data.frame(Variable = paste0(x, ref_lev),
                                t(data.frame(rep("", ncol(cox0)-1))),
                                row.names = 1, check.names = F)
                colnames(cox1) <- colnames(cox0)
                multi_cox <- rbind(cox1, cox0)
            }
            colnames(multi_cox) <- c("Variable", "Coefficient", "Events/N", "HR", 
                                     "lower95", "upper95", "HR for events (95% CI)", "P")
            multi_cox <- data.frame(Type = x, multi_cox, check.names=F)

            # recalculate events/N for character/factor variables
            if(is.character(train_df[,x]) | is.factor(train_df[,x])) {
                if(is.character(train_df[,x])) {
                    levs <- sort(unique(train_df[,x]))
                } else if(is.factor(train_df[,x])) {
                    levs <- levels(train_df[,x])
                }                        
                for(lev0 in levs) {
                    N = nrow(train_df[!is.na(train_df[,x]) & train_df[,x]==lev0,])
                    n = nrow(train_df[(!is.na(train_df[,x]) & train_df[,x]==lev0) &
                                    (!is.na(train_df$status) & train_df[,"status"]==1),])
                    if(lev0 == levs[1]) {
                        multi_cox[1, "Events/N"] <- paste0(n, "/", N)
                    } else {
                        multi_cox[multi_cox$Variable == lev0, "Events/N"] <- paste0(n, "/", N)
                    }
                }
            }
            multi_cox
        })
        res <- do.call(rbind.data.frame, multi_cox_out)
    } else {
        res <- data.frame()
    }
    res
}




## fixed ----------------------------------------------
stomach_pub_cp_colnames <- c("Patient_ID","Sample_ID","Dataset_ID","Platform","Molecular_profiling",
                                  "Primary_site","Disease_type","Sample_type","Is_primary_treatment",
                                  "Treatment_type","Treatment_setting","Regimen","Age","Gender","Ethnicity",
                                  "T","N","M","Stage_edition","Stage",
                                  "OS_months","OS_status","RFS_months","RFS_status",
                                  "DFS_months","DFS_status","DSS_months","DSS_status",
                                  "PFS_months","PFS_status","TTP_months","TTP_status",
                                  "Lauren","EBV","Molecular_subtype","pStage")

colon_pub_cp_colnames <- c("Patient_ID","Sample_ID","Dataset_ID","Platform","Molecular_profiling",
                                  "Primary_site","Disease_type","Sample_type","Is_primary_treatment",
                                  "Treatment_type","Treatment_setting","Regimen","Age","Gender","Ethnicity",
                                  "T","N","M","Stage_edition","Stage",
                                  "OS_months","OS_status","RFS_months","RFS_status",
                                  "DFS_months","DFS_status","DSS_months","DSS_status",
                                  "PFS_months","PFS_status","TTP_months","TTP_status",
                                  "MSI_status","Location")

fixed_variables <- c("Patient_ID","SampleID","Dataset_ID","Platform","Molecular_profiling",
                     "Primary_site","Disease_type","Sample_type","Is_primary_treatment",
                     "Treatment_type","Treatment_setting","Regimen","Age","Gender",
                     "Ethnicity","T","N","M","Stage_edition","Stage",
                     "OS_months","OS_status","RFS_months","RFS_status","DFS_months","DFS_status",
                     "DSS_months","DSS_status","PFS_months","PFS_status","TTP_months","TTP_status")
                     
exclude_variables <- c("Patient_ID", "SampleID", "Sample_ID", "Dataset_ID", "Platform", "Molecular_profiling",
                       "Primary_site", "Disease_type", "Sample_type", "Is_primary_treatment",
                       "Treatment_setting", "Regimen","Stage_edition",
                       "OS_months", "OS_status", "RFS_months", "RFS_status",
                       "DFS_months", "DFS_status", "DSS_months", "DSS_status",
                       "PFS_months", "PFS_status", "TTP_months", "TTP_status")



## ploting functions ----------------------------------------------
density.1group <- function(x, legend='', main='',xlab='', ylab='Density', 
                           border_col=rgb(192,192,0,max = 255),
                           polygon_col=rgb(192,192,0, 88,max = 255),
                           cex=1.5,cex.axis=1.5) {
  par(mar=c(4,4,4,2))
  dx<-density(x)
  ##
  xmin=(min(x) - (max(x)-min(x))/2)
  xmax=(max(x) + (max(x)-min(x))/2)
  ymin=min(dx$y)
  ymax=(max(dx$y)+(max(dx$y)-min(dx$y))/5)
  ##
  plot(dx,bty='n',xlab='',ylab='',xaxs='i',yaxs='i',main='',xlim=c(xmin,xmax),ylim=c(ymin,ymax),cex.axis=cex.axis)
  polygon(dx, col=polygon_col, border=border_col)
  legend('topright', legend=legend, ncol=1, bty='n', cex=cex)
  ##
  mtext(main,side=3,line=0.5,cex=cex)
  mtext(xlab,side=1,line=2.5,cex=cex)
  mtext(ylab,side=2,line=2.5,cex=cex)
}


##
forestplot.func <- function(signature_type, cox_table, continuous_variables = c(), labeltext = c("Value", "Events/N", "Hazard Ratio (95% CI)", "P"), cex=4) {
  
  require(forestplot)
  raw_col1 <- colnames(cox_table)[1]
  colnames(cox_table)[1] <- "Variable"
  
  if(signature_type=="Predictive"){
    # Bold sub-title
    cox_table$summ <- F
    cox_table[(cox_table[,"Events/N\n"]=="" & cox_table[,"HR"]==""),'summ'] <- T
    IS.summary <- c(T,cox_table$summ)
    #
    cox_table$summ1 <- "   "
    cox_table[(cox_table[,"Events/N\n"]=="" & cox_table[,"HR"]==""),'summ1'] <- ""
  }else{
    cox_table$summ <- F
    cox_table[(cox_table[,"HR"]==""),'summ'] <- T
    IS.summary <- c(T,cox_table$summ)
    #
    cox_table$summ1 <- "   "
    cox_table[(cox_table[,"HR"]==""),'summ1'] <- ""
  }
  #cox_table[cox_table$Variable %in% c("Age","age","RS","rs","Signature Score","Signature score","Signature_score"),'summ1'] <- ""
  cox_table[cox_table$Variable %in% c(continuous_variables, c("Age","age","RS","rs","Signature Score","Signature score","Signature_score")),'summ1'] <- ""
  cox_table$Value <- paste0(cox_table$summ1,cox_table$Variable)
  colnames(cox_table)[colnames(cox_table)=="Value"] <- raw_col1
  cox_colnames <- colnames(cox_table)
  cox_table <- rbind(cox_colnames, cox_table)
  cox_table$HR <- as.numeric(cox_table$HR)
  cox_table$lower95 <- as.numeric(cox_table$lower95)
  cox_table$upper95 <- as.numeric(cox_table$upper95)
  
  cox_table[cox_table=="Inf"] <- NA
  cox_labeltext = cox_table[,labeltext]
  
  forest_figure <- forestplot(labeltext=cox_labeltext, 
                              mean=cox_table$HR,   
                              lower=cox_table$lower95, 
                              upper=cox_table$upper95,  
                              zero=1,            
                              boxsize=0.3,
                              graphwidth=unit(40, 'mm'),
                              #col = fpColors(lines = "#990000", box = "#660000", zero = "darkblue"),
                              lwd.ci = 2,#
                              ci.vertices =T,
                              is.summary=IS.summary,
                              hrzl_lines=list("2" = gpar(lwd=2, col="black")),
                              txt_gp=fpTxtGp(
                                label=gpar(cex=1),
                                ticks=gpar(cex=1), 
                                xlab=gpar(cex=1.5), 
                                title=gpar(cex=2)),
                              graph.pos=3, cex=cex)
}



##
msig_km.plot <- function(df, b.group, b.signature.type, b.endpoint, png.file = NULL) {
    # df with columns IDs: SampleID, time, status, Risk_group, Dataset, subgroup
    ds <- unique(df$Dataset)
    kmplots <- list()
    z=1
    for(ds0 in ds) {
        # ds0 = ds[1]
        df1 <- df[df$Dataset == ds0,]
        rownames(df1) <- df1$SampleID
        rs0 <- df1$Risk_group
        names(rs0) <- df1$SampleID
        surv_df = df1[,c("time", "status")]
        leg0 = paste(paste0(names(table(rs0)), ' (n = '), table(rs0), ')',sep='')

        sbs0 <- unique(df1$subgroup)
        id0 <- paste0(ds0, " (", sbs0,")")
        title0 = ifelse(sbs0 == "all_pts", ds0, id0)

        ##
        control_group = NULL
        if("Low risk" %in% rs0) { control_group = "Low risk" }
       
        ##
        set.seed(2019)
        if(all(table(rs0)>=5) & n_distinct(rs0) > 1) {
            kmplots[[z]] <- survival_1VS2_plot_multi3(surv_df, rs0, control_group = control_group,
                                                        plot.sigOnly = FALSE,
                                                        pval.coord = c(0,0.1),
                                                        ylab = paste0("Probability of ", b.endpoint),
                                                        legend.title = "",
                                                        legend.labs = leg0, palette = get.color.km(leg0),
                                                        title = title0)
            z=z+1
        }
    }

    ##
    n0 = 1
    if(b.group == "2 groups") {
        n1 = ifelse(b.signature.type != "Predictive", 1, 3)
        width0 = ifelse(b.signature.type != "Predictive", 8, 24)
    } else {
        n1 = ifelse(b.signature.type != "Predictive", 1, 4)
        width0 = ifelse(b.signature.type != "Predictive", 8, 32)
    }
    ##
    if(!is.null(png.file)) {
        png(file=png.file, width = width0*300, height = 5*300, res=300)
    }
    arrange_ggsurvplots(kmplots, print=T, nrow = n0, ncol=n1)
    if(!is.null(png.file)) { dev.off() }
}



## get.color.km(c("Benefit (n=10)", "Non-benefit (n=20)"))
get.color.km = function(x){
    ids = list(test = c("High", "High risk", "Non-benefit", "No-benefit", "Control"),
               moderate = c("Moderate", "Intermediate risk"),
               control = c("Low", "Low risk", "Benefit", "Treatment"))

    palette = c("#d42a24", "orange",'#00427d')

    n = length(ids)
    cols = sapply(x, function(i) {
        # i = x[1]
        i <- strsplit(i, " \\(")[[1]][1]
        id0 <- grep(i, unlist(ids))
        if(length(id0) == 0){
            col=NA
        }else{
            for(j in 1:n){
                grep(i, ids[[j]])
                if(length(grep(i, ids[[j]])) > 0) {
                    col = palette[j]
                    break()
                }
            }
        }
        return(col)
    })

    cols[is.na(cols)] = palette[(n+1):(n+sum(is.na(cols)))]
    return(cols)
}




## from tt --
get_df_summary <- function(signature.type="Predictive",cp_df,dataset_name = NULL,return_data_type=FALSE){
  
  if(!is.null(cp_df)){
    treatment_type_report_id <- which(colnames(cp_df) %in% "Treatment_type_report")
	if(length(treatment_type_report_id) > 0){
		colnames(cp_df)[treatment_type_report_id] <- "Therapy"
	}
	
    if(signature.type == "Predictive"){
      fixed_char_colname <- c("Age","Gender","Ethnicity","Therapy","T","N","M","Stage")
    }else{
      fixed_char_colname <- c("Age","Gender","Ethnicity","Treatment_type","T","N","M","Stage")
    }
    
    fixed_colname <- c("Patient_ID","Sample_ID","Dataset_ID","Platform","Molecular_profiling","Primary_site","Disease_type","Sample_type","Is_primary_treatment","Treatment_type","Treatment_setting","Regimen","Age","Gender","Ethnicity","T","N","M","Stage_edition","Stage","OS_months","OS_status","RFS_months","RFS_status","DFS_months","DFS_status","DSS_months","DSS_status","PFS_months","PFS_status","TTP_months","TTP_status")
    fixed_char_colname <- intersect(fixed_char_colname,colnames(cp_df))
    custom_colname <- setdiff(colnames(cp_df),c(fixed_colname,"Therapy","pStage"))
    custom_char_colname <- setdiff(custom_colname,c("OS_months","RFS_months","DFS_months","DSS_months","PFS_months","TTP_months"))
    char_colname <- c(fixed_char_colname,custom_char_colname)
    char_colname <- unique(char_colname)
    cp_df[cp_df == "NA" | cp_df == "na" | cp_df == "Unknown" | cp_df == "unknown" | cp_df == "[Discrepancy]" | cp_df == "#NULL!" | cp_df == "?" | cp_df == ""] <- NA
    
    
    if(return_data_type){
      #data_types <- data.frame()
	  data_types <- data.frame(Characteristic= character(), Type= character(), stringsAsFactors=FALSE)
      for(i in char_colname){
		#print(i)
        if(sum(is.na(cp_df[,i])) == nrow(cp_df)){
          next
        }
		
        rows <- c()
        NA_id <- c()
		if(i == "Age") {
			type = "continuous"
			cp_df[,i] <- as.numeric(cp_df[,i])
		} else {
			i_vas <- as.numeric(cp_df[,i])
			#if(length(!is.na(i_vas)) > 0) {
			if(length(i_vas[!is.na(i_vas)]) > 0) {
				if(is.numeric(cp_df[,i]) & length(unique(cp_df[,i])) > 5){
					type = "continuous"
				}else{
					type = "categorical"
				}
			} else {
				type = "categorical"
			}
			
		}
		#print(type)
        #data_type=c(i,type)
		data_type=data.frame("Characteristic"=i,"Type"=type)
		data_types <- rbind.data.frame(data_types,data_type)
		#print(data_type)
        #data_types <- rbind(data_types,data_type)
		#colnames(data_types) <- c("Characteristic","Type")
      }
      #data_types <- as.data.frame(data_types)
      ##data_types <- data_types
      #colnames(data_types) <- c("Characteristic","Type")
      return(data_types)
      
    } else {
      
      all_rows <- c()
      
      for(i in char_colname){
        #i="Therapy"
        if(sum(is.na(cp_df[,i])) == nrow(cp_df)){
          next
        }
        rows <- c()
        NA_id <- c()
		if(i == "Age") {
			type = "continuous"
			cp_df[,i] <- as.numeric(cp_df[,i])
		} else {
			i_vas <- as.numeric(cp_df[,i])
			if(length(i_vas[!is.na(i_vas)]) > 0) {
				if(is.numeric(cp_df[,i]) & length(unique(cp_df[,i])) > 5){
					type = "continuous"
				}else{
					type = "categorical"
				}
			} else {
				type = "categorical"
			}
		}
        
        if(type == "continuous"){
          # Age or custom col ----
          num <- length(cp_df[!is.na(cp_df$Age),i])
          #if(i == "Age"){
          #  all <- c("Age (year)",num,NA,NA,NA)
          #}else{
          #  all <- c(i,num,NA,NA,NA)
          #}
          mean_row <- c("Mean",num,NA,mean(round(mean(cp_df[!is.na(cp_df[,i]) ,i]),0)),NA)  
          median_row <- c("Median",num,NA,median(round(mean(cp_df[!is.na(cp_df[,i]) ,i]),0)),NA)
          min_i <- min(cp_df[!is.na(cp_df[,i]),i])
          max_i <- max(cp_df[!is.na(cp_df[,i]),i])
          range_row <- c("Range",num,NA,paste0(min_i,"-",max_i),NA)
          #rows <- rbind(all,mean_row,median_row,range_row)
          rows <- rbind(mean_row,median_row,range_row)
          rows <- as.data.frame(rows)
          rows <- rows[,c(1,2,4,5)]
          col_names <- c("Number","All patients","P value")
          col_names <- c("Characteristic",paste0(dataset_name,"<br />",col_names[1:2]),paste0("Training dataset vs ",dataset_name,"<br />",col_names[3]))
          #col_names <- c("Number","Percent","All patients","P value")
          #col_names <- c("Characteristic",paste0(dataset_name,"<br />",col_names[1:3]),paste0("Training dataset vs ",dataset_name,"<br />",col_names[4]))
          colnames(rows) <- col_names
          rows$Key <- paste0(i,"__",rows$Characteristic)
          rows$Key1 <- i
          #rows <- rows[,-1]
          
        } else {
          
          cp_df[(is.na(cp_df[,i]) | cp_df[,i] == "na" | cp_df[,i] == "Unknown" | cp_df[,i] == "unknown" | cp_df[,i] == "[Discrepancy]" | cp_df[,i] == "#NULL!" | cp_df[,i] == "?" | cp_df[,i] == ""),i] <- "NA"
          #i_all <- c(i,NA,NA,NA,NA)
          df <- as.data.frame(table(cp_df[,i]))
          df$Var1 <- as.character(df$Var1)
          NA_id <- which(df$Var1=="NA")
          if(length(NA_id)==1) {
            df$Var1[NA_id] <- gsub("NA","zzzz",df$Var1[NA_id])
            vars <- c(df$Var1[-NA_id],df$Var1[NA_id])
            vars <- data.frame(vars)
            mer_df <- merge(vars,df,by.x = "vars",by.y = "Var1",all.x = TRUE, sort = FALSE)
          } else {
            mer_df <- df
          }
          colnames(mer_df) <- c("Characteristic","Number")
          mer_df$Percent <- round(mer_df$"Number"/sum(mer_df$"Number")*100,2)
          mer_df$"All patients" <- paste0(mer_df$"Number"," (",mer_df$Percent,")")
          mer_df$P <- NA
          #rows <- rbind(i_all,mer_df)
          rows <- mer_df
          rows <- as.data.frame(rows)
          rows <- rows[,c(1,2,4,5)]
          col_names <- c("Number","All patients","P value")
          col_names <- c("Characteristic",paste0(dataset_name,"<br />",col_names[1:2]),paste0("Training dataset vs ",dataset_name,"<br />",col_names[3]))
          #col_names <- c("Number","Percent","All patients","P value")
          #col_names <- c("Characteristic",paste0(dataset_name,"<br />",col_names[1:3]),paste0("Training dataset vs ",dataset_name,"<br />",col_names[4]))
          colnames(rows) <- col_names
          rows$Key <- paste0(i,"__",rows$Characteristic)
          rows$Key1 <- i
          
        }
        
        all_rows <- rbind(all_rows,rows)
        all_rows <- unique(all_rows)
      }
      rownames(all_rows) <- all_rows$Key
      return(all_rows)
    }
  } else {
    return(NULL)
  }
}


all_sub <- function(colname,subclass,df) {
  x <- as.data.frame(table(df[,colname]))
  if(nrow(x)==1 & x[1,"Var1"]=="NA"){
    x = "NA"
  } else {
    sub=subclass
    non_sub = sub[!sub %in% x[,1]]
    if(length(non_sub)>0){
      non_sub = as.data.frame(non_sub)
      non_sub$Freq=0
      colnames(non_sub)=c("Var1","Freq")
      x = rbind(x,non_sub)
      x=t(x)
      colnames(x)=x[1,]
      x=x[,sub]
      x=t(x)
    }else{
      x=t(x)
      colnames(x)=x[1,]
      x=x[,sub]
      x=t(x)
    }
    #x=paste0(x[,1]," ",x[,2],collapse = "\r\n")
    x = paste0(x[,2],collapse = ",")
    x = gsub(" ","",x)
  }
  return(x)
}

age_summary <- function(age_value) {
  age_value <- as.numeric(age_value)
  age_value <- age_value[!is.na(age_value)]
  if(length(age_value) > 0) {  
  	mean_age <- round(mean(age_value),0)
  	mean_age_text <- paste0("<b>Mean</b>: ",mean_age)
  	median_age <- round(median(age_value),0)
  	median_age_text <- paste0("<b>Median</b>: ",median_age)
  	min_age <- round(min(age_value),0)
  	max_age <- round(max(age_value),0)
  	range_age_text <- paste0("<b>Range</b>: ",min_age,"-",max_age)
  	age <- paste0(c(mean_age_text,median_age_text,range_age_text),collapse = "<br />")				  
  } else {
  	age <- "NA"
  }
  return(age)
} 

get_tripod_dataset_summary <- function(tripod_dataset_infos, user_filtered_datasets_info0, b.have.profiles, select_training_dataset, used_validation_dataset) {
  if(nrow(tripod_dataset_infos) > 0) {
    geo_url_id <- which(colnames(mtcars) %in% "geo_url")
    if(length(geo_url_id) > 0) {
       subset(mtcars, select = -geo_url_id)
    }
    tripod_datasets_summary.df <- merge.data.frame(user_filtered_datasets_info0,tripod_dataset_infos,by.x = "dataset_name",by.y="dataset_name",all.x=TRUE,sort=F)
    for(i in 1:nrow(tripod_datasets_summary.df)){
      #if(!is.na(tripod_datasets_summary.df[i,"pubmed_id"])){
	  if(tripod_datasets_summary.df[i,"pubmed_id"] != "" & !is.na(tripod_datasets_summary.df[i,"pubmed_id"])) {
		pubmed_ids <- strsplit(tripod_datasets_summary.df[i,"pubmed_id"],", ")[[1]]
		#pubmed_ids <- gsub(" ","",pubmed_ids)
		tripod_datasets_summary.df[i,"pubmed_url"] <- paste0(paste0("https://pubmed.ncbi.nlm.nih.gov/",pubmed_ids),collapse = ", ")
        #tripod_datasets_summary.df[i,"pubmed_url"] <- paste0("https://pubmed.ncbi.nlm.nih.gov/",tripod_datasets_summary.df[i,"pubmed_id"])
      } else {
        tripod_datasets_summary.df[i,"pubmed_url"] <- ""
      }
    }
    
  } else {
    if(b.have.profiles == "Yes"){
      user_filtered_datasets_info0$pubmed_id <- ""
      #user_filtered_datasets_info0$molecular_database <- ""
      user_filtered_datasets_info0$accession_number <- ""
      user_filtered_datasets_info0$accession_url <- ""
      #user_filtered_datasets_info0$clinical_database <- ""
      user_filtered_datasets_info0$clinical_number <- ""
      user_filtered_datasets_info0$clinical_url <- ""
      
    } else{
      user_filtered_datasets_info0$pubmed_id <- ""
      #user_filtered_datasets_info0$clinical_database <- ""
      user_filtered_datasets_info0$clinical_number <- ""
      user_filtered_datasets_info0$clinical_url <- ""
    }
    
    tripod_datasets_summary.df <- user_filtered_datasets_info0
	tripod_datasets_summary.df$pubmed_url <- ""
  }
  
  
  if(b.have.profiles == "Yes") {
    tripod_datasets_summary_df=B_train_vali_data(select_training_dataset=select_training_dataset,select_validation_dataset=used_validation_dataset,userupload_sqlite_path = FALSE, user_filtered_datasets_info = tripod_datasets_summary.df,filtered_datasets_type="B_User-filtered datasets") #%>%
      #mutate(pubmed_id = glue::glue(as.character(tags$a(href ="{pubmed_url}",target="_blank","{pubmed_id}")))) #%>%
	  tripod_datasets_summary_df = tripod_datasets_summary_href(df = tripod_datasets_summary_df,url_colname="pubmed_url",number_colname="pubmed_id")
    tripod_datasets_summary_df = tripod_datasets_summary_href(df = tripod_datasets_summary_df,url_colname="accession_url",number_colname="accession_number")
	  tripod_datasets_summary_df = tripod_datasets_summary_href(df = tripod_datasets_summary_df,url_colname="clinical_url",number_colname="clinical_number")
    #colnames(tripod_datasets_summary_df) <- c("Dataset name","Dataset ID","GEO url","Number of pts","Primary site","Disease type","Molecular profiling","Endpoint","Pubmed","Molecular dataset","Molecular dataset URL","Clinical trial Registry","Clinical trial URL","Pubmed URL","Class")
    colnames(tripod_datasets_summary_df) <- c("Dataset name","Dataset ID","Number of pts","Primary site","Disease type","Molecular profiling","Endpoint","Pubmed","Molecular dataset","Molecular dataset URL","Clinical trial Registry","Clinical trial URL","Pubmed URL","Class")
    tripod_datasets_summary_df <- tripod_datasets_summary_df[,c("Dataset name","Primary site","Disease type","Molecular dataset","Endpoint","Pubmed","Clinical trial Registry","Class")]
  } else {
    tripod_datasets_summary_df=B_train_vali_data(select_training_dataset=select_training_dataset,select_validation_dataset=used_validation_dataset,userupload_sqlite_path = FALSE, user_filtered_datasets_info = tripod_datasets_summary.df,filtered_datasets_type="B_User-filtered datasets") #%>%
      #mutate(pubmed_id = glue::glue(as.character(tags$a(href ="{pubmed_url}",target="_blank","{pubmed_id}")))) 
    tripod_datasets_summary_df = tripod_datasets_summary_href(df = tripod_datasets_summary_df,url_colname="pubmed_url",number_colname="pubmed_id")  
	  tripod_datasets_summary_df = tripod_datasets_summary_href(df = tripod_datasets_summary_df,url_colname="clinical_url",number_colname="clinical_number")
    #colnames(tripod_datasets_summary_df) <- c("Dataset name","Dataset ID","GEO url","Number of pts","Primary site","Disease type","Endpoint","Pubmed","Clinical trial Registry","Clinical trial URL","Pubmed URL","Class")
    colnames(tripod_datasets_summary_df) <- c("Dataset name","Dataset ID","Number of pts","Primary site","Disease type","Endpoint","Pubmed","Clinical trial Registry","Clinical trial URL","Pubmed URL","Class")
    tripod_datasets_summary_df <- tripod_datasets_summary_df[,c("Dataset name","Primary site","Disease type","Endpoint","Pubmed","Clinical trial Registry","Class")]
  }
  return(tripod_datasets_summary_df)
}

# tripod item
excluded_sample_info <- function(df,column_name,input,title=TRUE) {
    excluded.sample <- setdiff(unique(df[,column_name]),input)
    excluded.sample[is.na(excluded.sample)] = "NA"
    column_name <- gsub("_", " ",column_name)
    if(column_name %in% c("T","N","M")){
      column_name0 <- paste0(column_name," stage")
      column_name <- column_name0
    }else if(column_name == "Stage"){
      column_name0 <- column_name
    } else{
      column_name0 <- tolower(column_name)
    }
    
    if(length(excluded.sample) == 0){
      excluded.samples <- ""
    } else {
      if("NA" %in% excluded.sample){
        excluded.sample <- c(sort(excluded.sample[!"NA" %in% excluded.sample]),"NA")
        excluded.sample0 <- paste0("missing ",column_name0)
        excluded.sample1 <- excluded.sample[!excluded.sample %in% "NA"]
        if(column_name == "Stage"){
          excluded.samples <- paste0(c(paste0("stage ",excluded.sample1),excluded.sample0),collapse=", ")
        } else{
          excluded.samples <- paste0(c(excluded.sample1,excluded.sample0),collapse=" or ")
        }
        
      } else {
        excluded.sample <- sort(excluded.sample)
        if(column_name == "Stage"){
          excluded.samples <- paste0(paste0("stage ",excluded.sample),collapse=" or ")
        } else {
          excluded.samples <- paste0(excluded.sample,collapse=" or ")
        }
        
      }
    }
    if(!column_name %in% c("T stage","N stage","M stage","Stage")){
      excluded.samples <- tolower(excluded.samples)
    }
    if(title){
      if(excluded.samples != ""){
        excluded.samples <- paste0("&nbsp;&nbsp;&nbsp;&nbsp;-  ",column_name,": ",excluded.samples)
      }
    }
    return(excluded.samples)
  }
  
  included_sample_info <- function(df,column_name,input,title=TRUE) {
    included.sample <- intersect(unique(df[,column_name]),input)
    included.sample[is.na(included.sample)] = "NA"
    column_name <- gsub("_", " ",column_name)
    if(column_name %in% c("T","N","M")){
      column_name0 <- paste0(column_name," stage")
      column_name <- column_name0
    }else if(column_name == "Stage"){
      column_name0 <- column_name
    } else {
      column_name0 <- tolower(column_name)
    }
    
    if(length(included.sample) == 0){
      included.samples <- ""
    } else {
      if("NA" %in% included.sample){
        included.sample <- c(sort(included.sample[!"NA" %in% included.sample]),"NA")
        included.sample0 <- paste0("missing ",column_name0)
        included.sample1 <- included.sample[!included.sample %in% "NA"]
        if(column_name == "Stage"){
          included.samples <- paste0(c(paste0("stage ",included.sample1),included.sample0),collapse=", ")
        } else{
          included.samples <- paste0(c(included.sample1,included.sample0),collapse=" or ")
        }
        
      } else {
        included.sample <- sort(included.sample)
        if(column_name == "Stage"){
          included.samples <- paste0(paste0("stage ",included.sample),collapse=" or ")
        } else {
          included.samples <- paste0(included.sample,collapse=" or ")
        }
      }
    }
    if(!column_name %in% c("T stage","N stage","M stage","Stage")){
      included.samples <- tolower(included.samples)
    }
    if(title){
      if(included.samples != ""){
        included.samples <- paste0("&nbsp;&nbsp;&nbsp;&nbsp;-  ",column_name,": ",included.samples)
      }
    }
    return(included.samples)
  }


##
get_summary_test <- function(train_valid_df_list, signature.type, used_validation_dataset=used_validation_dataset) {
  train_valid_df_rows_list <- list()
  train_valid_df_rows_list_names <- c()
  key_df <- data.frame()
  all_data_types <- data.frame()
  
  for (j in names(train_valid_df_list)) {
    print(j)
    cp_df <- train_valid_df_list[[j]]
    
    if(signature.type=="Predictive"){
      control_df <- cp_df[cp_df$Treatment_type == "Control",]
      treatment_df <- cp_df[cp_df$Treatment_type == "Treatment",]
      
      df_rows <- get_df_summary(cp_df=cp_df,dataset_name=j)
      control_df_rows <- get_df_summary(cp_df=control_df,dataset_name=paste0(j,"<br />Control"))
      treatment_df_rows <- get_df_summary(cp_df=treatment_df,dataset_name=paste0(j,"<br />Treatment"))
      train_valid_df_rows_list_names <- names(train_valid_df_rows_list)
      train_valid_df_rows_list <- c(train_valid_df_rows_list,list(df_rows),list(control_df_rows),list(treatment_df_rows))
      names(train_valid_df_rows_list) <- c(train_valid_df_rows_list_names,j,paste0(j,"_control"),paste0(j,"_treatment"))
      key_df <- rbind(key_df, df_rows[, c("Key1","Key","Characteristic")])
      
      # output data type
      j_data_types <- get_df_summary(cp_df=cp_df,dataset_name=j,return_data_type=TRUE)
      all_data_types <- rbind(all_data_types,j_data_types)
      
    }else{
      df_rows <- get_df_summary(signature.type="Prognostic",cp_df=cp_df,dataset_name=j)
      train_valid_df_rows_list_names <- names(train_valid_df_rows_list)
      train_valid_df_rows_list <- c(train_valid_df_rows_list,list(df_rows))
      names(train_valid_df_rows_list) <- c(train_valid_df_rows_list_names,j)
      key_df <- rbind(key_df, df_rows[, c("Key1","Key","Characteristic")])
      
      # output data type
      j_data_types <- get_df_summary(signature.type="Prognostic",cp_df=cp_df,dataset_name=j,return_data_type=TRUE)
      all_data_types <- rbind(all_data_types,j_data_types)
    }
  }
  
  key_df <- unique(key_df)
  key_df1 <- key_df[order(key_df$Key1,key_df$Key),"Key"]
  key_df1<- as.data.frame(key_df1)
  colnames(key_df1) <- "Key"
  
  all_data_types <- as.data.frame(all_data_types)
  all_data_types <- unique(all_data_types)
  
  
  res0 <- data.frame()
  for(i in 1:length(train_valid_df_rows_list)) {
    
    i_df <- train_valid_df_rows_list[[i]]
    ids <- which(colnames(i_df) %in% c("Key1","Characteristic"))
    if(length(ids) > 0) {
      i_df <- i_df[-ids]
    }
    if(i == 1) {
      res0 <- merge.data.frame(key_df1,i_df,by = "Key",all.x = TRUE,sort = FALSE)
    } else {
      res0 <- merge.data.frame(res0,i_df,by = "Key",all.x = TRUE,sort = FALSE)
    }
    
    res0 <- res0[order(res0$Key),]
    
  }
  
  res01 <- do.call(rbind.data.frame,strsplit(res0$Key,'__'))
  colnames(res01) <- c("Level", "Characteristic")
  res0 <- cbind(res01,res0)
  
  add_rows <- data.frame()
  for(m in 1:nrow(all_data_types)) {
    if(all_data_types[m,"Characteristic"] == "Age"){
      add_row <- c("Age","Age (year)","Age__!!!!!",rep(NA,ncol(res0)-3))
    }else{
      add_row <- c(all_data_types[m,"Characteristic"],all_data_types[m,"Characteristic"],paste0(all_data_types[m,"Characteristic"],"__!!!!!"),rep(NA,ncol(res0)-3))
    }
    
    add_rows <- rbind(add_rows,add_row)
    
  }
  
  colnames(add_rows) <- colnames(res0)
  
  res00 <- rbind(add_rows,res0)
  class(res00)
  res00 <- res00[order(res00$Key),]
  
  if(signature.type=="Predictive"){
    
    add_cols <- paste0(names(train_valid_df_list),"<br />Control vs Treatment<br />P value")
    for(aa in add_cols){
      res00$aa <- NA
      colnames(res00)[colnames(res00)=="aa"] <- aa
    }
  }
  
  # p
  train_df <- train_valid_df_list[["Training dataset"]]
  
  for(n in 1:nrow(all_data_types)){
    print(n)
	#n = 1
    char <- all_data_types[n,"Characteristic"]
    type <- all_data_types[n,"Type"]
    col_id <- paste0(char,"__!!!!!")
    
    sub_res00_z <- res00[res00$Level==char,]
    sub_res00 <- sub_res00_z[sub_res00_z$Characteristic!="zzzz",]
    training_v_name <- paste0("Training dataset","<br />Number")
    
    
    if(nrow(sub_res00) > 2){
      
      # p value within control and treatment
      if(signature.type=="Predictive"){
        for(ct in 1:length(train_valid_df_list)) {
          #ct = 1
          cp_df <- train_valid_df_list[[ct]]
          control_df <- cp_df[cp_df$Therapy == "Control",]
          treatment_df <- cp_df[cp_df$Therapy == "Treatment",]
          ct_p_name <- add_cols[ct]
          control_v_name <- paste0(names(train_valid_df_list)[ct],"<br />Control<br />Number")
          treatment_v_name <- paste0(names(train_valid_df_list)[ct],"<br />Treatment<br />Number")
          
          if(type == "continuous") {
            control_v  <- control_df[!is.na(control_df[,char]),char]
            treatment_v <- treatment_df[!is.na(treatment_df[,char]),char]
            if(length(control_v) > 0 & length(treatment_v) > 0) {
				ct_v <- t.test(control_v,treatment_v)
				pval <- ifelse(ct_v$p.value<0.01,
							format(ct_v$p.value, digits=3, scientific=TRUE),
							round(ct_v$p.value, 3))
				
				res00[res00$Key==col_id,ct_p_name] = pval
            }
          } else {
            
            control_v <- sub_res00[-1,control_v_name]
            control_nona <- length(control_v[!is.na(control_v)])
            treatment_v <- sub_res00[-1,treatment_v_name]
            treatment_nona <- length(treatment_v[!is.na(treatment_v)])
            
            if(control_nona > 0 & treatment_nona > 0){
              ct_mat <- sub_res00[-1,c(control_v_name,treatment_v_name)]
              ct_mat[is.na(ct_mat)] = 0
              colnames(ct_mat) <- c("co","t")
              mode(ct_mat$co) <- "numeric"
              mode(ct_mat$t) <- "numeric"
              ct_mat <- as.matrix(ct_mat)
              
              ct_v <- fisher.test(ct_mat, workspace = 2e5, simulate.p.value = TRUE, B = 1e5)
              pval <- ifelse(ct_v$p.value<0.01,
                             format(ct_v$p.value, digits=3, scientific=TRUE),
                             round(ct_v$p.value, 3))
              
              res00[res00$Key==col_id,ct_p_name] = pval
              
            }
          }
        }
      }
      
      # p value between datasets
      for(va in used_validation_dataset){
        
        valid_df <- train_valid_df_list[[va]]
        valid_v_name <- paste0(va,"<br />Number")
        
        p_name <- paste0("Training dataset vs ",va,"<br />P value")
        
        if(type == "continuous"){
          if(char %in% colnames(train_df)) {
			training_v <- as.numeric(train_df[!is.na(as.numeric(train_df[,char])),char])
		  } else {
			training_v <- c()
		  }
		  if(char %in% colnames(valid_df)) {
			valid_v <- as.numeric(valid_df[!is.na(as.numeric(valid_df[,char])),char])
		  } else {
			valid_v <- c()
		  }
          
          if(length(training_v) > 0 & length(valid_v) > 0) {
			t_v <- t.test(training_v,valid_v)
			pval <- ifelse(t_v$p.value<0.01,
							format(t_v$p.value, digits=2, scientific=TRUE),
							round(t_v$p.value, 3))
			
			res00[res00$Key==col_id,p_name] = pval
          }
        }else{
          training_v <- sub_res00[-1,training_v_name]
          training_nona <- length(training_v[!is.na(training_v)])
          valid_v <- sub_res00[-1,valid_v_name]
          valid_nona <- length(valid_v[!is.na(valid_v)])
          num_nona <- training_nona+valid_nona
          if(training_nona > 0 & valid_nona > 0 & num_nona > 2){
            t_v_mat <- sub_res00[-1,c(training_v_name,valid_v_name)]
            t_v_mat[is.na(t_v_mat)] = 0
            colnames(t_v_mat) <- c("t","v")
            mode(t_v_mat$t) <- "numeric"
            mode(t_v_mat$v) <- "numeric"
            t_v_mat <- as.matrix(t_v_mat)
            
            t_v <- fisher.test(t_v_mat, workspace = 2e5, simulate.p.value = TRUE, B = 1e5)
			pval <- ifelse(t_v$p.value<0.01,
                           format(t_v$p.value, digits=2, scientific=TRUE),
                           round(t_v$p.value, 3))
            
            res00[res00$Key==col_id,p_name] = pval
          }
        }
      }
    }
  }
  
  res00 <- res00[order(res00$Key),]
  
  valid_keep_cols <- c()
  if(signature.type == "Predictive"){
    train_keep_cols <- c("Level","Characteristic","Training dataset<br />All patients","Training dataset<br />Control<br />All patients","Training dataset<br />Treatment<br />All patients","Training dataset<br />Control vs Treatment<br />P value")
    for(va in used_validation_dataset){
      valid_keep_col <- c(paste0(va,"<br />All patients"),paste0(va,"<br />Control<br />All patients"),paste0(va,"<br />Treatment<br />All patients"),paste0(va,"<br />Control vs Treatment<br />P value"))
      valid_keep_cols <- c(valid_keep_cols,valid_keep_col)
    }
    if(length(used_validation_dataset) > 0){
      valid_keep_col1 <- paste0("Training dataset vs ",used_validation_dataset,"<br />P value")
    }else{
      valid_keep_col1 <- NULL
    }
    train_valid_keep_cols <- c(train_keep_cols,valid_keep_cols,valid_keep_col1)
    
  }else{
    train_keep_cols <- c("Level","Characteristic","Training dataset<br />All patients")
    for(va in used_validation_dataset){
      valid_keep_col <- c(paste0(va,"<br />All patients"),paste0("Training dataset vs ",va,"<br />P value"))
      valid_keep_cols <- c(valid_keep_cols,valid_keep_col)
    }
    train_valid_keep_cols <- c(train_keep_cols,valid_keep_cols)
  }
  
  res00 <- res00[,train_valid_keep_cols]
  res00$Characteristic <- gsub("zzzz","Unknown",res00$Characteristic)
  res00$Characteristic <- gsub("Therapy","Therapy",res00$Characteristic)
  res00$Level <- gsub("Therapy","Therapy",res00$Level)
  
  if(signature.type == "Predictive"){
    train_valid_keep_cols <- gsub("Control<br />All patients","Control",train_valid_keep_cols)
    train_valid_keep_cols <- gsub("Treatment<br />All patients","Treatment",train_valid_keep_cols)
  }else{
    train_valid_keep_cols <- gsub("<br />All patients","",train_valid_keep_cols)
  }
  
  colnames(res00) <- train_valid_keep_cols
  colnames(res00)[1] <- "Variable"
  res00 <- data.frame("Raw Rank" = 1:nrow(res00), res00, check.names=F)
  res00 <- res00[,c("Variable","Raw Rank",colnames(res00)[!colnames(res00) %in% c("Variable","Raw Rank")])]
  
  return(res00)
  
}

# set color to table
get_format_datatable <- function(train_valid_df_list_names, signature.type, df, analysis = "Participant characteristic", train_name = "Training dataset", b.group = "2 groups", b.times=c(12,36,60), b.endpoint="OS", js_dir) {
  
  if(nrow(df) > 0) {
  format_colus <- c()
  format_colu_list <- c()
  format_colu_list_name <- c()
  format_colors <- c()
  format_border_list <- list()
  format_border_list_name <- c()
  merge_header_col0 <- c()
  th_value_list <- list()
  th_name_list <- list()
  
  rank_id <- which(colnames(df) %in% "Raw Rank")
  if(length(rank_id > 0)){
    df <- df %>% select(-c(rank_id))
  }
  
  if(analysis == "Participant characteristic"){
	colnames(df)[which(colnames(df) %in% "Characteristic")] <- "Value"
	colnames(df)[which(colnames(df) %in% "Training dataset")] <- train_name
  }
  
  if(analysis == "uni_multi_cox"){
	if("Hazard Ratio (95% CI)<br />Multivariate analysis" %in% colnames(df)) { 
		train_valid_df_list_names <- c("Univariate analysis","Multivariate analysis")
		id <- which(colnames(df) %in% "No. of Events<br />Multivariate analysis")
		if(length(id > 0)){
			df <- df %>% select(-c(id))
		}
	} else {
		train_valid_df_list_names <- "Univariate analysis"
	}
	colnames(df) <- gsub("Characteristic","Value",colnames(df))
	colnames(df) <- gsub("Value<br />\\(Control, Treatment)","Value\n(Control, Treatment)",colnames(df))
	
	variable_id <- which(colnames(df) %in% "Variable")
	if(length(variable_id) > 0){
      df <- df[,c("Variable",colnames(df)[2:ncol(df)])]
    }
  } else {
    used_validation_dataset <- train_valid_df_list_names[-1]
  }
  
  rep_color_time <- ceiling((length(train_valid_df_list_names)-1)/7)
  format_dt_color1 <- c('#34621d',rep(c('#d34b00','#dd8e0a','#e6c710','#9aa440','#347e86','#304e98','#816ab3'),rep_color_time))
  format_dt_color2 <- c('#96c999',rep(c('#ffd7b9','#f9cc82','#f9eda6','#d7dda6','#84c8d0','#88a1da','#ccc2e1'),rep_color_time))
  
  # set color
  for(i in 1:length(train_valid_df_list_names)){
    #i =1
    i_name <- train_valid_df_list_names[i]
    
    if(analysis == "Participant characteristic"){
      if(signature.type=="Predictive"){
        
        format_colu <- paste0(i_name,c("<br />All patients","<br />Control","<br />Treatment","<br />Control vs Treatment<br />P value"))
        rep_time <- 4
        
        # border
        format_border <- paste0(i_name,"<br />Control vs Treatment<br />P value")
        
      }else{
        if(i_name == "Training dataset"){
          rep_time <- 1
          format_colu <- train_name
          # border
          format_border <- train_name
        }else{
          rep_time <- 2
          format_colu <- c(i_name,paste0("Training dataset vs ",i_name,"<br />P value"))
          # border
          format_border <- paste0("Training dataset vs ",i_name,"<br />P value")
        }
        
      }
      
      uncol_num <- length(which(c("Variable","Value") %in% colnames(df)))
      pa_id <- i+uncol_num-1
      color_id <- format_dt_color1[i]
      merge_header_col <- paste0("$(thead).closest('thead').find('th').eq(",pa_id,").css({'backgroundColor':'",color_id,"','color':'white','text-align':'center'})") 
      merge_header_col0 <- c(merge_header_col0,merge_header_col)
      
    }else if(analysis == "Kaplan_Meier_analysis"){
      if(signature.type == "Predictive"){
        df <- df[,c("Predictive group","Group",colnames(df)[3:ncol(df)])]
	  }
      
	  times_title <- paste0(b.times,"-months Survival Probability (%, 95% Cl)<br />")
	  format_colu <- paste0(c("All Patients (No. of Events)<br />",times_title),i_name)
	  rep_time <- length(b.times)+1
	  
      # header color
      uncol_num <- length(which(c("Group","Predictive Groups","Predictive Group","Predictive groups","Predictive group","Risk Groups","Risk Group","Signature group","Signature Group") %in% colnames(df)))
      pa_id <- i+uncol_num-1
      color_id <- format_dt_color1[i]
      merge_header_col <- paste0("$(thead).closest('thead').find('th').eq(",pa_id,").css({'backgroundColor':'",color_id,"','color':'white','text-align':'center'})")  
      merge_header_col0 <- c(merge_header_col0,merge_header_col)
      
      # border
      format_border <- paste0(b.times[length(b.times)],"-months Survival Probability (%, 95% Cl)<br />",i_name)
	  
    }else if(analysis == "Association analysis"){
      
      colnames(df) <- gsub("Level","Value",colnames(df))
      df <- df[,c("Variable",colnames(df)[2:ncol(df)])]
      
      if(signature.type == "Predictive"){
        
		if(b.group == "2 groups") {
		    format_colu <- paste0(c("Benefit<br />","No-benefit<br />","P<br />","Benefit<br />Control<br />","Benefit<br />Treatment<br />","Benefit<br />P<br />","No-benefit<br />Control<br />","No-benefit<br />Treatment<br />","No-benefit<br />P<br />","P-subgroup<br />"),i_name)
			rep_time <- 10
        } else {
			format_colu <- paste0(c("Benefit<br />","Intermediate group<br />","No-benefit<br />","P<br />","Benefit<br />Control<br />","Benefit<br />Treatment<br />","Benefit<br />P<br />","Intermediate group<br />Control<br />","Intermediate group<br />Treatment<br />","Intermediate group<br />P<br />","No-benefit<br />Control<br />","No-benefit<br />Treatment<br />","No-benefit<br />P<br />","P-subgroup<br />"),i_name)
			rep_time <- 14
		}
		uncol_num <- length(which(c("Variable","Value") %in% colnames(df)))
        pa_id <- i+uncol_num-1
        color_id <- format_dt_color1[i]
        merge_header_col <- paste0("$(thead).closest('thead').find('th').eq(",pa_id,").css({'backgroundColor':'",color_id,"','color':'white','text-align':'center'})") 
        merge_header_col0 <- c(merge_header_col0,merge_header_col)
        
        # border
		if(b.group == "2 groups") {
			format_border_p <- paste0("P<br />",i_name)
			format_border_b <- paste0("Benefit<br />P<br />",i_name)
			format_border_n <- paste0("No-benefit<br />P<br />",i_name)
			format_border_ps <- paste0("P-subgroup<br />",i_name)
			
			if(!format_border_p %in% colnames(df)) {
				format_border_p <- NULL
			}	
			if(!format_border_b %in% colnames(df)) {
				format_border_b <- NULL
			}
			if(!format_border_n %in% colnames(df)) {
				format_border_n <- NULL
			}
			if(!format_border_ps %in% colnames(df)) {
				format_border_ps <- NULL
			}	
		
		} else {
			format_border_p <- paste0("P<br />",i_name)
			format_border_b <- paste0("Benefit<br />P<br />",i_name)
			format_border_i <- paste0("Intermediate group<br />P<br />",i_name)
			format_border_n <- paste0("No-benefit<br />P<br />",i_name)
			format_border_ps <- paste0("P-subgroup<br />",i_name)
			
			if(!format_border_p %in% colnames(df)) {
				format_border_p <- NULL
			}	
			if(!format_border_b %in% colnames(df)) {
				format_border_b <- NULL
			}
			if(!format_border_i %in% colnames(df)) {
				format_border_i <- NULL
			}
			if(!format_border_n %in% colnames(df)) {
				format_border_n <- NULL
			}
			if(!format_border_ps %in% colnames(df)) {
				format_border_ps <- NULL
			}	
		
		}
		
      }else{
		if(b.group == "2 groups") {
			format_colu <- paste0(c("High risk<br />","Low risk<br />","P<br />"),i_name)
			rep_time <- 3
        } else {
			format_colu <- paste0(c("High risk<br />","Intermediate risk<br />","Low risk<br />","P<br />"),i_name)
			rep_time <- 4
		}
		uncol_num <- length(which(c("Variable","Value") %in% colnames(df)))
        
        pa_id <- i+uncol_num-1
        color_id <- format_dt_color1[i]
        merge_header_col <- paste0("$(thead).closest('thead').find('th').eq(",pa_id,").css({'backgroundColor':'",color_id,"','color':'white','text-align':'center'})") 
        merge_header_col0 <- c(merge_header_col0,merge_header_col)
        
        # border
        format_border <- paste0("P<br />",i_name)
      }
      
    } else if(analysis == "uni_multi_cox") {
      if(signature.type == "Predictive"){
        format_colu <- paste0(c("Hazard Ratio (95% CI)<br />","P<br />","P interaction<br />"),i_name)
        rep_time <- 3
        
        # border
        format_border <- paste0("P interaction<br />",i_name)
        
      }else{
        format_colu <- paste0(c("Hazard Ratio (95% CI)<br />","P<br />"),i_name)
        rep_time <- 2
        
        # border
        format_border <- paste0("P<br />",i_name)
        
      }
      
      uncol_num <- length(which(c("Variable","Value\n(Control, Treatment)","Value","No. of Patients","No. of Events") %in% colnames(df)))
      
      pa_id <- i+uncol_num-1
      color_id <- format_dt_color1[i]
      merge_header_col <- paste0("$(thead).closest('thead').find('th').eq(",pa_id,").css({'backgroundColor':'",color_id,"','color':'white','text-align':'center'})") 
      merge_header_col0 <- c(merge_header_col0,merge_header_col)
      
    }
    
    # body color
    format_colus <- c(format_colus,format_colu)
    format_colu_list_name <- names(format_colu_list)
    format_colu_list <- c(format_colu_list,list(format_colu))
    names(format_colu_list) <- c(format_colu_list_name,i_name)
    
    format_color <- rep(format_dt_color2[i],rep_time)
    format_colors <- c(format_colors,format_color)
    
    if(analysis == "Association analysis" & signature.type == "Predictive"){
		if(b.group == "2 groups") {
			format_border_list <-c(format_border_list,list(format_border_p,format_border_b,format_border_n,format_border_ps))
		} else {
			format_border_list <-c(format_border_list,list(format_border_p,format_border_b,format_border_i,format_border_n,format_border_ps))
		}
    }else{
      format_border_list <- c(format_border_list,list(format_border))
    }
    
  }
  
  pc_compared_p_header_col <- c()
  compared_p_header_name <- c()
  compared_p_num <- 0
  
  rowspan_header_col0 <- c()
  ro_id <- 0
  for(ro in c("Variable","Level","Characteristic","Value","Group","Groups","Value\n(Control, Treatment)","No. of Patients","No. of Events")){
    rowspan_header_col <- paste0("$(thead).closest('thead').find('th').eq(",ro_id,").css({'backgroundColor':'#D0CECE','color':'black','text-align':'center'})") 
    rowspan_header_col0 <- c(rowspan_header_col0,rowspan_header_col)
    ro_id <- ro_id+1
  }
  
  
  if(analysis == "Participant characteristic") {
	if(length(train_valid_df_list_names) > 1){
	  valid_names <- train_valid_df_list_names[-1]
	  valid_names <- paste0("Validation dataset ",1:length(valid_names))
	} else {
	  valid_names <- NULL
	}
	train_valid_df_list_names0 <- c("Training dataset",valid_names)

    # merge column and row
    th_value_list <- list(htmltools::withTags(th(rowspan = 2, 'Variable')),
                          htmltools::withTags(th(rowspan = 2, 'Value')))
    
    # compared p header	
    if(signature.type=="Predictive"){
      for (na in train_valid_df_list_names){
        th_value_list <- c(th_value_list,list(htmltools::withTags(th(colspan = 4, na))))
      }
      
      if(length(train_valid_df_list_names) > 1){
        compared_p_num <- 1
        compared_p_col_id <- uncol_num+length(train_valid_df_list_names)
        pc_compared_p_header_col <- paste0("$(thead).closest('thead').find('th').eq(",compared_p_col_id,").css({'backgroundColor':'#34621d','color':'white','text-align':'center'})")
        
        
        predictive_dt_p_colus <- paste0("Training dataset vs ",c(train_valid_df_list_names[-1]),"<br />P value")
        format_colus <- c(format_colus,predictive_dt_p_colus)
        format_colors <- c(format_colors,rep('#96c999',length(train_valid_df_list_names)-1))
        compared_p_header_name <- paste0("P\n(Trainging dataset vs ",used_validation_dataset,")")
        th_value_list <- c(th_value_list,list(htmltools::withTags(th(colspan = length(train_valid_df_list_names), 'Dataset Comparison (All Patients)'))))
      }
      th_name_list <- lapply(c(rep(c('All Patients', 'Control',
                                     'Treatment','P'), length(train_valid_df_list_names)),compared_p_header_name), htmltools::withTags(th))
    } else {
      
      th_names <- c()
	  na_id = 0
      for (na in train_valid_df_list_names){
        na_id = na_id+1
		na0 <- train_valid_df_list_names0[na_id]
        if(na == "Training dataset"){
          th_name <- train_name
          th_value_list0 <- list(htmltools::withTags(th(colspan = 1, na0)))
          
        }else{
          th_name <- c(na,"P")
          th_value_list0 <- list(htmltools::withTags(th(colspan = 2, na0)))
        }
        th_names <- c(th_names,th_name)
        th_value_list <- c(th_value_list,th_value_list0)
      }
      
      th_name_list <- lapply(th_names, htmltools::withTags(th))
    } 
    
    rowsGroup_list <- list(which(colnames(df) %in% "Variable")-1)
    
  } else if (analysis == "Kaplan_Meier_analysis") {
    
    if(signature.type == "Predictive") {
      th_value_list <- list(htmltools::withTags(th(rowspan = 2, 'Predictive group')),
                            htmltools::withTags(th(rowspan = 2, 'Group')))
      
    } else {
      th_value_list <- list(
                            htmltools::withTags(th(rowspan = 2, 'Signature group')))
    }
    
    colspan_num <- length(b.times)+1
    for (na in train_valid_df_list_names) {
      th_value_list <- c(th_value_list,list(htmltools::withTags(th(colspan = colspan_num, na))))
    }
    
	times_title <- paste0(b.times,"-months ",b.endpoint," (%, 95% Cl)")
	th_name_list <- lapply(rep(c('All Patients (No. of Events)', times_title), length(train_valid_df_list_names)), htmltools::withTags(th))
    rowsGroup_list <- list(which(colnames(df) %in% "Predictive group")-1)
    
  } else if (analysis == "Association analysis") {
    th_value_list <- list(htmltools::withTags(th(rowspan = 2, 'Variable')),
                          htmltools::withTags(th(rowspan = 2, 'Value')))
    
    if(signature.type == "Predictive") {
      for (na in train_valid_df_list_names) {
		if(b.group == "2 groups"){
			th_value_list <- c(th_value_list,list(htmltools::withTags(th(colspan = 10, na))))
		} else {
			th_value_list <- c(th_value_list,list(htmltools::withTags(th(colspan = 14, na))))
		}
      }
      
	  if(b.group == "2 groups") {
		th_name_list <- lapply(rep(c('Benefit','No-benefit','P','Benefit\nControl','Benefit\nTreatment','Benefit\nP',
                                   'No-benefit\nControl','No-benefit\nTreatment','No-benefit\nP','P-subgroup'), 
                                 length(train_valid_df_list_names)), htmltools::withTags(th))
      } else {
		th_name_list <- lapply(rep(c('Benefit','Intermediate group','No-benefit','P',
									'Benefit\nControl','Benefit\nTreatment','Benefit\nP',
									'Intermediate group\nControl','Intermediate group\nTreatment','Intermediate group\nP',
                                    'No-benefit\nControl','No-benefit\nTreatment','No-benefit\nP','P-subgroup'), 
                                 length(train_valid_df_list_names)), htmltools::withTags(th))
	  }
      
    } else {
      
      for (na in train_valid_df_list_names){
		if(b.group == "2 groups"){
			th_value_list <- c(th_value_list,list(htmltools::withTags(th(colspan = 3, na))))
		} else {
			th_value_list <- c(th_value_list,list(htmltools::withTags(th(colspan = 4, na))))
		}
      }
      
	  if(b.group == "2 groups") {
		th_name_list <- lapply(rep(c('High risk', 'Low risk',
                                   'P'), length(train_valid_df_list_names)), htmltools::withTags(th))
	  } else {
		th_name_list <- lapply(rep(c('High risk', 'Intermediate risk', 'Low risk',
                                   'P'), length(train_valid_df_list_names)), htmltools::withTags(th))
	  }
      
    }
    
    rowsGroup_list <- list(which(colnames(df) %in% "Variable")-1)
    
    
  } else if (analysis == "uni_multi_cox") {
    Characteristic_id <- which(colnames(df) %in% c("Value\n(Control, Treatment)","Value"))
    if(length(Characteristic_id > 0)){
      Characteristic_name <- colnames(df)[Characteristic_id]
      th_value_list <- list(htmltools::withTags(th(rowspan = 2, 'Variable')),
                            htmltools::withTags(th(rowspan = 2, Characteristic_name)),
                            htmltools::withTags(th(rowspan = 2, 'No. of Patients')),
                            htmltools::withTags(th(rowspan = 2, 'No. of Events')))
    } else {
      Characteristic_name <- colnames(df)[1]
      th_value_list <- list(htmltools::withTags(th(rowspan = 2, 'Variable')),
                            htmltools::withTags(th(rowspan = 2, Characteristic_name)),
                            htmltools::withTags(th(rowspan = 2, 'No. of Patients')),
                            htmltools::withTags(th(rowspan = 2, 'No. of Events')))
    }
    for (na in train_valid_df_list_names) {
      if(signature.type == "Predictive") {
        th_value_list <- c(th_value_list,list(htmltools::withTags(th(colspan = 3, na))))
        th_name_list <- lapply(rep(c('Hazard Ratio (95% CI)', 'P', 'P interaction'), length(train_valid_df_list_names)), htmltools::withTags(th))
      } else {
        th_value_list <- c(th_value_list,list(htmltools::withTags(th(colspan = 2, na))))
        th_name_list <- lapply(rep(c('Hazard Ratio (95% CI)', 'P'), length(train_valid_df_list_names)), htmltools::withTags(th))
      }
    }
    rowsGroup_list <- list(which(colnames(df) %in% "Variable")-1)
  }
  
  # set header color
  notmerge_header_col0 <- c()
  
  for(colu in 1:length(format_colus)) {
    te_name <- format_colus[colu]
    te_id <- which(colnames(df) %in% te_name)-1
    te_id <- te_id+length(train_valid_df_list_names)+compared_p_num
    color_id <- format_colors[colu]
    notmerge_header_col <- paste0("$(thead).closest('thead').find('th').eq(",te_id,").css('backgroundColor', '",color_id,"')")
    notmerge_header_col0 <- c(notmerge_header_col0,notmerge_header_col)
  }
  
  headerCallback_texts <- paste0(c(rowspan_header_col0,merge_header_col0,notmerge_header_col0,pc_compared_p_header_col),collapse = ";")
  headerCallback_texts <- paste0("function( thead, data, start, end, display ) {",headerCallback_texts,"}")
  
  
  # merge column and row
  sketch = htmltools::withTags(table(
    class = 'display',
    thead(
      tr(
        th_value_list
      ),
      tr(
        th_name_list
      )
    )
  ))
  
  #
  if(nrow(df) > 10) {
    scrollY_value <- 500
  } else {
    scrollY_value = FALSE
  }
  
  data <- DT::datatable(df, escape=FALSE, rownames=FALSE, extensions = c('RowGroup','Buttons'), selection="single", container = sketch, #class = 'cell-border stripe',
                        options = list(scrollX = TRUE, scrollY = scrollY_value, dom = 'Bt',buttons = c('csv', 'excel', 'print'), pageLength = 100, 
                                       columnDefs = list(list(className = 'dt-center', targets = "_all")), #
                                       rowsGroup = rowsGroup_list, ordering=F, 
                                       headerCallback = JS(headerCallback_texts)
                                       
                        )) 
  rowspan_format_border <- colnames(df)[which(colnames(df) %in% c("Variable","Characteristic","Value","Level","Risk Group","Risk Groups","Group","Signature Group","Signature group","Predictive Group","Predictive Groups","Predictive group","Predictive groups","Value\n(Control, Treatment)","No. of Patients","No. of Events"))] 
  
  format_border_list <- c(format_border_list,list(rowspan_format_border))
  for(dt in 1:length(format_border_list)){
    data <- data %>%
      DT::formatStyle(format_border_list[[dt]], borderRight = '1px solid #ddd')
  }
  
  path <- js_dir
  dep <- htmltools::htmlDependency(
    "RowsGroup", "2.0.0", 
    path, script = "dataTables.rowsGroup.js")
  data$dependencies <- c(data$dependencies, list(dep))
  data
  
  
  return(data)
  
  } else {
	return(df)
  }
}

# rowgroup column 
get_rowgroup_table <- function(df, coef = TRUE, have.profiles, train.profiles_type = NULL) {
	rank_id <- which(colnames(df) %in% "Raw Rank")
	if(length(rank_id > 0)){
		df <- df %>% select(-c(rank_id))
	}
	
	if(coef) {
		
		df$Coefficient <- as.numeric(df$Coefficient) 
		df$Coefficient <- round(as.numeric(df$Coefficient),3)
		df$P <- as.numeric(df$P) 
		df$P <- ifelse(df$P<0.01, format(df$P, digits=3, scientific=TRUE), round(df$P,3))
		
		if(have.profiles == "No") {
			df <- df[,c("Variable",colnames(df)[2:ncol(df)])]
			format_border_list <- list("Variable")
			rowsGroup_list <- list(which(colnames(df) %in% "Variable")-1)
			
			if(nrow(df) > 15){
				data <- DT::datatable({df}, escape=FALSE, rownames=FALSE, extensions = c('RowGroup','Buttons'), selection="single",options = list(scrollX = TRUE, scrollY = 500, dom = 'Bt', buttons = c('csv', 'excel', 'print'), pageLength = 100, rowsGroup = rowsGroup_list, columnDefs = list(list(className = 'dt-left', targets = "_all"))))
			}else{
				data <- DT::datatable({df}, escape=FALSE, rownames=FALSE, extensions = c('RowGroup','Buttons'), selection="single", options = list(scrollX = TRUE, dom = 'Bt', buttons = c('csv', 'excel', 'print'), pageLength = 100, rowsGroup = rowsGroup_list, columnDefs = list(list(className = 'dt-left', targets = "_all"))))
			}
			
			for(dt in 1:length(format_border_list)){
				data <- data %>%
				DT::formatStyle(format_border_list[[dt]], borderRight = '1px solid #ddd')
			}	
		} else {
			
			if(!is.null(train.profiles_type)){
				if(train.profiles_type == "Protein"){
					colnames(df)[which(colnames(df) %in% "Characteristic")] <- "Protein"
				}	else if(train.profiles_type == "Non-coding") {
					colnames(df)[which(colnames(df) %in% "Characteristic")] <- "ncRNA"
				} else {
					colnames(df)[which(colnames(df) %in% "Characteristic")] <- "Gene symbol"
				}
			}
			
			if(nrow(df) > 15){
				data <- DT::datatable({df}, escape=FALSE, rownames=FALSE, extensions = 'Buttons', selection="single", options = list(scrollX = TRUE, scrollY = 500, dom = 'Bt', buttons = c('csv', 'excel', 'print'), pageLength = 100, columnDefs = list(list(className = 'dt-left', targets = "_all"))))
			}else{
				data <- DT::datatable({df}, escape=FALSE, rownames=FALSE, extensions = 'Buttons', selection="single", options = list(scrollX = TRUE, dom = 'Bt', buttons = c('csv', 'excel', 'print'), pageLength = 100, columnDefs = list(list(className = 'dt-left', targets = "_all"))))
			}
		}
	} else {
		df$P <- as.numeric(df$P)
		df$P <- ifelse(df$P<0.01, format(df$P, digits=3, scientific=TRUE), round(df$P,3))
		if(have.profiles == "Yes") {
			if(!is.null(train.profiles_type)){
				if(train.profiles_type == "Protein"){
					colnames(df)[which(colnames(df) %in% "Variable")] <- "Protein"
				}	else if(train.profiles_type == "Non-coding") {
					colnames(df)[which(colnames(df) %in% "Variable")] <- "ncRNA"
				} else {
					colnames(df)[which(colnames(df) %in% "Variable")] <- "Gene symbol"
				}
			}
			 
			data <- DT::datatable({df}, escape=FALSE, rownames=FALSE, extensions = 'Buttons', selection="single", options = list(scrollX = TRUE, dom = 'Blfrtip', buttons = c('csv', 'excel', 'print'), pageLength = 10, columnDefs = list(list(className = 'dt-left', targets = "_all"))))
		} else {
			format_border_list <- list("Variable")
			rowsGroup_list <- list(which(colnames(df) %in% "Variable")-1)
			
			data <- DT::datatable({df}, escape=FALSE, rownames=FALSE, extensions = c('RowGroup','Buttons'), selection="single",options = list(scrollX = TRUE, dom = 'Blfrtip', buttons = c('csv', 'excel', 'print'), pageLength = 10, rowsGroup = rowsGroup_list, columnDefs = list(list(className = 'dt-left', targets = "_all"))))
			
			for(dt in 1:length(format_border_list)){
				data <- data %>%
				DT::formatStyle(format_border_list[[dt]], borderRight = '1px solid #ddd')
			}	
		}
		
	}

    return(data)
}
		

##

# query selected training/validation datasets (User-filtered datasets) from user upload sqlite ----
B_train_vali_data <- function(select_training_dataset,select_validation_dataset,userupload_sqlite_path,user_filtered_datasets_info = FALSE,filtered_datasets_type){
  if(userupload_sqlite_path == FALSE) {
    user_filtered_sqlite <- user_filtered_datasets_info
  } else {
    mydb_user <- DBI::dbConnect(RSQLite::SQLite(),dbname=userupload_sqlite_path)
    user_filtered_sqlite <- DBI::dbGetQuery(mydb_user,paste0("select * from userupload_datasets where dataset_group='",filtered_datasets_type,"'")) 
    dbDisconnect(mydb_user)
  }
  train_data = user_filtered_sqlite[user_filtered_sqlite$dataset_name==select_training_dataset,]
  if(nrow(train_data) > 0) {
    train_data$class="Training dataset"
  }
  vali_datas=c()
  if(length(select_validation_dataset) > 0){ 
    for(vali in select_validation_dataset) {
      vali_data = user_filtered_sqlite[user_filtered_sqlite$dataset_name==vali,]
      vali_datas = rbind(vali_datas,vali_data)
    }
    if(nrow(vali_datas) > 0) {
      vali_datas$class="Validation dataset"
    }
  }
  train_vali_data=rbind(train_data,vali_datas)
  return(train_vali_data)
}


##
  
tripod_datasets_summary_href <- function(df, url_colname, number_colname, urlsep = ", ", idsep = ", "){
  if(nrow(df) > 0) {
    for(i in 1:nrow(df)) {
      #i = 4
      not_ref <- FALSE
      sub_df = df[i,]
      if(!is.na(sub_df[,url_colname]) & !is.null(sub_df[,number_colname])) {
        urls <- strsplit(sub_df[,url_colname],urlsep)[[1]]
        
        numbers <- strsplit(sub_df[,number_colname],idsep)[[1]]
        
        if(length(urls) == 1 || length(numbers) == 1) {
          
          if(length(urls) == 1) {
            if(length(numbers) == 0){
              numbers0 <- numbers
            } else if(length(numbers) == 1) {
              numbers0 <- strsplit(numbers," ")[[1]]
            } else {
              numbers0 <- strsplit(numbers[1]," ")[[1]]
            }
          } 
          
          if(length(numbers) == 1) {
            if(length(urls) == 0) {
              not_ref <- TRUE
            } else if (length(urls) >= 1) {
              urls <- urls[1]
            }
          }
          
          href_number <- numbers0[grep("^GSE",numbers0)]
          if(length(href_number)==1) {
            other_text <- paste0(numbers0[!numbers0 %in% href_number],collapse = " ")
          } else {
            href_number = numbers0[1]
            other_text <- paste0(numbers0[!numbers0 %in% href_number],collapse = " ")
          }
          if(not_ref) {
            df[i,number_colname] = df[i,number_colname]
            
          } else {
            if(length(other_text) > 0) {
              df[i,number_colname] = glue::glue(as.character(tags$a(href ="{urls}",target="_blank","{href_number}"))," ",{other_text})
            } else {
              df[i,number_colname] = glue::glue(as.character(tags$a(href ="{urls}",target="_blank","{href_number}")))
            }
          }
          
          
        } else {
          hrefs <- sapply(numbers,function(x) {
            #x = "GSE14333 (218 in total)"
            id = which(numbers %in% x)
            url = urls[id]
            url <- url[!is.na(url)]
            x0 <- strsplit(x," ")[[1]]
            href_x <- x0[grep("^GSE",x0)]
            if(length(href_x)==1) {
              x_other_text <- paste(x0[!x0 %in% href_x],collapse = " ")
            } else {
              href_x = x0[1]
              x_other_text <- paste(x0[!x0 %in% href_x],collapse = " ")
            }
            if(length(url) >= 1) {
              url <- url[1]
              if(length(x_other_text) > 0) {
                href = glue::glue(as.character(tags$a(href ="{url}",target="_blank","{href_x}"))," ",{x_other_text})
              } else {
                href = glue::glue(as.character(tags$a(href ="{url}",target="_blank","{href_x}")))
              }
              
            } else {
              url <- NA
              href = x
            }
            
            return(href)
          })
          
          df[i,number_colname] = paste0(hrefs,collapse = ", ")
        }
      }
    }
  }  
  return(df)
}


# change colors of dataset names in TRIPOD report  
style_datasetname <- function(dataset_name,dataset_name0,dataset_style,text){
  text0 <- strsplit(text, " ")[[1]]
  change_text = FALSE
  if(dataset_name0[1] %in% text0){
    dataset_name1 = dataset_name0
    name_col = "name0"
    combined_col = "combined0"
    change_text = TRUE
  } 
  #if(dataset_name[1] %in% text0){
  #  dataset_name1 = dataset_name
  #  name_col = "name"
  #  combined_col = "combined"
  #  change_text = TRUE
  #}
  
  if(change_text){
    inum=0
    for(i in dataset_name1){
      inum = inum+1
      i0 <- paste0("^",i,"$")
      new_i <- dataset_style[dataset_style[,name_col]==i,combined_col]
      text <- strsplit(text, " ")[[1]]
      text <- paste(gsub(i0, new_i, text), collapse = " ")
    }
  }
  
  return(text)
}