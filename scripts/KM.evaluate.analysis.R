# Rscript /opt/shiny-server/apps/repromsig/scripts/KM.evaluate.analysis.R /opt/shiny-server/apps/repromsig/ColoGuide_Stage_II_local/output/sig.ini


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
	user.config.ini.file <- "/opt/shiny-server/apps/repromsig/ColoGuide_Stage_II_local/output/sig.ini"
}

# 2) Library
library(survival)
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

##
script.dir <- get_value.ccb(config_file = user.config.ini.file,  key = 'script_dir')[[1]]
model_analysis_rdata.file <- get_value.ccb(config_file = user.config.ini.file, key = 'model_analysis_rdata')[[1]]


# 4) define/source functions
source(paste0(script.dir,'ccb.helper.R'))

# readin and preprocess -------------------------------------------
load(model_analysis_rdata.file)
# cgwtools::lsdata(model_analysis_rdata.file)
mysetwd(paste0(outpath,'/external_evaluate/'))
rdata.filename <- paste0(b.name.sig, '.external_evaluate_analysis_rdata.RData')

if(nrow(coef1) > 0) {
    model_coef <- coef1
    model_coef <- model_coef %>% filter(!is.na(P) & P != "")

    # preprocess valid sets -------------------------------------------
    if((b.have.profiles == "Yes" & length(valid.exp.list) >= 1) |
       (b.have.profiles == "No" & length(valid.cp.list) >= 1)) {

        ifelse(b.have.profiles == "Yes",
               valid.list <- valid.exp.list,
               valid.list <- valid.cp.list
               )
        testsets.name <- gsub("_exp$", '', names(valid.list))

        ##
        test_rs_plot.list = list()
        test_surv_prob_report.list = list()
        surv_prob_key_df.list = list()
        test_pts_ann_sb.list = list()
        test_surv.list = list()

        CP_Stat.list <- list()
        CP_Stat_key_df.list <- list()
        valid_warning.list <- list()
        z=1 
        z0=1

        for(I in 1:length(testsets.name)) {
            ## cp table
            d_id <- testsets.name[I]
            test_cp <- valid.cp.list[[d_id]]
            colnames(test_cp)[1:length(fixed_variables)] <- fixed_variables
            rownames(test_cp) <- test_cp$SampleID
            
            test_cp[test_cp == "NA" | test_cp == "na" | test_cp == "<NA>" | test_cp == "Unknown" | test_cp == "unknown" | test_cp == "[Discrepancy]"] <- NA
            test_cp <- test_cp[, which(apply(test_cp, 2, function(x) { n_distinct(unique(x[!is.na(x)])) }) >= 2)]
            debug(test_cp)

            ##
            variables0 <- colnames(test_cp)[!colnames(test_cp) %in% exclude_variables]
            variables0 <- intersect(c("SampleID", paste0(b.endpoint, "_months"),
                                    paste0(b.endpoint, "_status"), variables0),
                                    colnames(test_cp))

            ## Age and stage process
            df <- test_cp
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
            test_cp <- df
            debug(test_cp)

            ##
            pts_ann_filter0 <- test_cp[, variables0]    
            colnames(pts_ann_filter0)[1:3] <- c('SampleID', 'time', 'status')
            mode(pts_ann_filter0$time) <- "numeric"
            mode(pts_ann_filter0$status) <- "numeric"
            pts_ann_filter0$time <- (pts_ann_filter0$time*30)
            colnames(pts_ann_filter0) <- gsub(" ", "_", colnames(pts_ann_filter0))
            debug(pts_ann_filter0)

            ## patient for subgroup analysis
            test_pts_ann_sb <- data.frame()
            if(!is.null(b.subgroup.variable)) {
                if(b.subgroup.variable %in% colnames(pts_ann_filter0)) {
                    test_pts_ann_sb <- unique(pts_ann_filter0[, c("SampleID", b.subgroup.variable)])
                }
            }
            test_pts_ann_sb.list[[z]] <- test_pts_ann_sb

            ##
            test_surv0 <- pts_ann_filter0[, c('SampleID', 'time', 'status')] %>%
                                mutate(status = as.integer(status)) %>%
                                filter(status %in% c(0,1))
            debug(test_surv0)
            test_surv.list[[z]] <- test_surv0[,-1]

            ## exp, pts.surv and pts ann
            if(b.have.profiles == "Yes") {
                Test0 <- valid.list[[d_id]]
                logtype <- valid.log.list[[d_id]]
                colnames(Test0)[1] <- 'Symbol'

                ## unique by gene
                f <- data.frame(table(Test0$Symbol))
                r <- f %>% filter(Freq >=2) %>% select(Var1) %>% unlist() %>% as.character()
                s <- f %>% filter(Freq == 1) %>% select(Var1) %>% unlist() %>% as.character()

                if(length(r) > 0) {
                    exp_rec <- aggregate(Test0[Test0$Symbol %in% r, -1],list(Test0[Test0$Symbol %in% r, 1]), mean)
                    colnames(exp_rec)[1] <- "Symbol"
                    Test0 <- rbind(Test0[Test0$Symbol %in% s, ], exp_rec)
                }
                Test0 <- Test0 %>% tibble %>% mycolumn_to_rownames("Symbol")
                debug(Test0, 5)

                ## 
                if(LASSO_cox) {
                    features <- coef_lasso_raw[,1]
                } else {
                    features <- model_coef$Variable
                } 
                siggenes_overlap <- intersect(features, rownames(Test0))
                siggenes_exist <- ifelse(length(siggenes_overlap) >= 1, TRUE, FALSE)
            } else {
                Test0 <- pts_ann_filter0
                variable0 <- unique(model_coef$Type)
                variables <- colnames(Test0)

                siggenes_overlap <- do.call(c, lapply(variable0, function(x) {
                    if((is.character(Train[,x]) | is.factor(Train[,x])) &  (x %in% colnames(Test0))) {
                        levs0 <- sort(unique(Train[,x]))
                        levs1 <- sort(unique(Test0[,x]))
                        if(all(levs0 %in% levs1) & all(levs1 %in% levs0)) {x}
                    } else if(!(is.character(Train[,x]) | is.factor(Train[,x]))) {
                        x
                    }
                }))
                siggenes_exist <- all(variable0 %in% siggenes_overlap)
            }
            debug(Test0,3)
            siggenes_exist
            

            if(siggenes_exist) {
                # trans into log2 for validation dataset
                if(b.have.profiles == "Yes") {
                    if(grepl('nonlog', logtype)) {Test0 <- log2((Test0[,-1]+1))}
                    if(grepl('log10', logtype)) {Test0 <- log2((10^(Test0[,-1])+1))}
                    if(grepl('log2|Not applicable', logtype)) {Test0 <- Test0[,-1]}

                    Test0 <- data.frame(t(Test0), check.names = F)
                    Test0[is.na(Test0)] <- mean(apply(Test0,1,mean, na.rm=T))
                    debug(Test0,5)

                    ##
                    Test0.mean <- apply(Test0, 1, mean)
                    df <- data.frame(Test0[,siggenes_overlap], check.names = F)
                    colnames(df) <- siggenes_overlap
                    rownames(df) <- rownames(Test0)

                    if(LASSO_cox) {
                        na.marker <- coef_lasso_raw[!coef_lasso_raw[,1] %in% colnames(df), 1]
                    } else {
                        na.marker <- model_coef[!model_coef$Variable %in% colnames(df), 'Variable']
                    }
                    
                    if(length(na.marker)>=1) {
                        na.marker.df <- data.frame(matrix(rep(Test0.mean,
                                                   length(na.marker)),
                                                   ncol=length(na.marker),
                                                   nrow=nrow(df)),
                                                   check.names = F)
                        colnames(na.marker.df) <- na.marker
                        df <- data.frame(df, na.marker.df, check.names = F)
                    }
                    Test0 <- df
                    debug(Test0)
                } else {
                    Test0 <- Test0[, siggenes_overlap]
                    variables <- colnames(Test0)
                    res <- do.call(cbind.data.frame, lapply(variables, function(x) {
                                x <- Test0[,x]
                                ifelse(is.numeric(x) & length(x) >= 5, x <- x, x <- as.character(x))
                                x
                            }))
                    colnames(res) <- variables
                    rownames(res) <- rownames(Test0)
                    Test0 <- res

                    ##
                    attach(Test0)
                    continuous_x <- unlist(lapply(variables, function(x) {if(is.numeric(Test0[,x])) {x}}))
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
                        rownames(df_continuous) <- rownames(Test0)
                        debug(df_continuous)
                    }

                    ##
                    categ_x <- do.call(c, lapply(variables, function(x) { if(is.character(Test0[,x]) | is.factor(df[,x])) {x}}))
                    df_categ <- data.frame()
                    if(!is.null(categ_x)) {
                        categ_x <- do.call(c, lapply(categ_x, function(x) {
                                id0 <- grep(x, model_coef$Type)
                                if(length(id0) >= 1) { x }
                        }))
                        df_categ <- data.frame(Test0[, categ_x], check.names = F)
                        colnames(df_categ) <- categ_x
                        rownames(df_categ) <- rownames(Test0)
                    }

                    list2comb <- list(df_continuous, df_categ)
                    res <- data.frame(SampleID = rownames(Test0))
                    for(i in 1:length(list2comb)) {
                        if(nrow(list2comb[[i]]) > 1) {
                            res <- cbind(res, list2comb[[i]])
                        }
                    }
                    Test0 <- res %>% tibble %>% column_to_rownames("SampleID")
                    detach(Test0)
                }
                debug(Test0)

                ##
                sams <- intersect.multi(list(rownames(Test0), rownames(pts_ann_filter0), rownames(test_surv0)))
                test_surv0 <- test_surv0[sams, ]
                df <- data.frame(Test0[sams, ], check.names = F)
                rownames(df) <- sams
                colnames(df) <- colnames(Test0)
                Test0 <- df
                pts_ann_filter0 <- pts_ann_filter0[sams, ]
                
                # calculate risk score using coef for validation dataset -------------------------------------------
                if(b.have.profiles == "Yes") {
                    if(standard_cox) {
                        require(survival)
                        risk_score_test <- data.frame(SampleID = rownames(Test0), 
                                                      RS = predict(fit_cv, newdata = Test0, type="lp")) %>%
                                                      filter(!is.na(RS))
                    } else if(LASSO_cox & !(b.name.sig %in% c("ColoGuidePro", "GeneExpressScore"))) {
                        require(glmnet)
                        risk_score_test <- data.frame(SampleID = rownames(Test0), 
                                                      RS = as.numeric(predict(fit_cv, newx = as.matrix(Test0), s = lambda, type = "link"))) %>%
                                                      filter(!is.na(RS))
                    } else if(b.name.sig %in% c("ColoGuidePro", "GeneExpressScore")) {
                        Test0 <- Test0[, coef1$Variable]
                        risk_score_test <- data.frame(SampleID = rownames(Test0), RS = apply(Test0, 1, function(x) {
                            x <- as.numeric(x)
                            sum(x * coef1$Coefficient)
                        })) %>%
                        filter(!is.na(RS))
                    }
                } else {
                        require(survival)
                        risk_score_test <- data.frame(SampleID = rownames(Test0), 
                                                      RS = predict(fit_cv, newdata = Test0, type="lp")) %>%
                                                      filter(!is.na(RS))
                }

                ##
                rs0 <- risk_score_test$RS
                names(rs0) <- risk_score_test$SampleID
                ##
                g1 <- ifelse(b.signature.type == "Predictive", "Benefit", "High risk")
                g2 <- ifelse(b.signature.type == "Predictive", "No-benefit", "Low risk")
                g3 <- ifelse(b.signature.type == "Predictive", "Intermediate group", "Intermediate risk")

                ## KM plot split by multi-cutoffs risk score
                if(b.cutpoint.method == "Percentile") {
                    if(b.group == "2 groups") {
                        high <- as.numeric(b.high2)/100
                        moderate <- NULL
                        low <- as.numeric(b.low2)/100
                    } else {
                        high <- as.numeric(b.high3)/100
                        moderate <- b.moderate3
                        low <- as.numeric(b.low3)/100
                    }
                    top = quantile(rs0, probs = high)
                    bottom = quantile(rs0, probs = low)
                } else {
                    top = top_train
                    bottom = bottom_train
                    moderate <- NULL
                }

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
                rs0 <- rs0[!is.na(rs0$RS_group),]
                rs1 <- rs0$RS_group
                names(rs1) <- rownames(rs0)

                stratify <- data.frame(SampleID = names(rs1), Type = rs1)
                debug(stratify)
                print(table(stratify[,2]))

                ##
                if(b.signature.type == "Predictive") {
                    cr <- pts_ann_filter0[,c("SampleID", "Treatment_type")]
                    test_stratify <- merge.data.frame(stratify, cr, by='SampleID', all.x = T)
                } else {
                    test_stratify <- stratify
                }
                print(table(test_stratify[,2]))
                test_rs_plot.list[[z]] <- test_stratify

                # Estimating x-year survival -------------------------------------------
                test_surv_prob <- merge.data.frame(test_surv0, test_stratify, by = "SampleID", all = F)
                types <- sort(unique(test_surv_prob$Type))
                h1 <- ifelse(b.signature.type == "Predictive", "Group", "Signature group")
                test_surv_prob.list = list()
                key_df = data.frame()
                
                for(x in types) {
                    df0 <- test_surv_prob[test_surv_prob$Type == x, ]
                    if(b.signature.type == "Predictive") { treat_types <- sort(unique(df0$Treatment_type)) }

                    ## for type
                    yrs <- b.times/12
                    res1 <- do.call(cbind, lapply(yrs, function(n0) {
                        tmp = df0 %>% filter(time > 365*n0)
                        if(max(df0$time) > 365*n0 & nrow(tmp) > 1) {
                            sf0 <- summary(survfit(Surv(time, status) ~ 1, data = df0), times = 365*n0)
                            surv_p = round(sf0$surv*100, 2)
                            surv_ci = paste0(round(sf0$lower * 100, 2), " - ", round(sf0$upper * 100, 2))
                            surv_ci = paste0(surv_p, " (", surv_ci, ")")
                        } else {
                            surv_ci = NA
                        }
                        surv_ci
                    }))
                    N = paste0(nrow(df0), " (", nrow(df0[df0$status == 1, ]), ")")
                    res1 <- data.frame(x, x, N, res1, check.names = F)
                    colnames(res1) <- c("Key", h1,
                                         paste0("All Patients (No. of Events)<br />", d_id),
                                         paste0(b.times, "-months Survival Probability (%, 95% Cl)<br />", d_id)
                                        )

                    ## for predictive signature
                    if(b.signature.type == "Predictive") {
                        for(treat_type0 in treat_types) {
                            df1 <- df0[df0$Treatment_type == treat_type0, ]
                            N1 <- paste0(nrow(df1), " (", nrow(df1[df1$status == 1, ]), ")")
                            res2 <- do.call(cbind, lapply(yrs, function(n0) {
                                tmp <- df1 %>% filter(time > 365*n0)
                                if(max(df1$time) > 365*n0 & nrow(tmp) > 1) {
                                    sf1 <- summary(survfit(Surv(time, status) ~ 1, data = df1), times = 365*n0)
                                    surv_ci1 <- paste0(round(sf1$surv*100, 2), " (",
                                                        paste0(round(sf1$lower * 100, 2), " - ",
                                                        round(sf1$upper * 100, 2)), ")")
                                } else {
                                    surv_ci1 = NA
                                }
                                surv_ci1
                            }))
                            res2 <- data.frame(paste0(x, "_", treat_type0), treat_type0, N1, res2, check.names = F)
                            colnames(res2) <- c("Key", h1,
                                                paste0("All Patients (No. of Events)<br />", d_id),
                                                paste0(b.times, "-months Survival Probability (%, 95% Cl)<br />", d_id)
                                               )
                            res1 <- rbind(res1, res2)
                        }
                    }

                    if(b.signature.type == "Predictive") {
                        res1 <- data.frame("Predictive group" = x, res1, check.names = F)
                        key_df <- rbind(key_df, res1[, 1:3])
                    } else {
                        key_df <- rbind(key_df, res1[, 1:2])
                    }
                    test_surv_prob.list <- c(test_surv_prob.list, list(res1))
                    print(x)
                }

                test_surv_prob_report <- do.call(rbind.data.frame, test_surv_prob.list)
                test_surv_prob_report.list[[z]] <- test_surv_prob_report
                surv_prob_key_df.list[[z]] <- key_df

                # Clinicopathological Association  -------------------------------------------   
                cp_rs <- test_stratify[, c('SampleID', 'Type')]
                cp_rs <- merge.data.frame(cp_rs, pts_ann_filter0, by='SampleID', all.x=T) %>% filter(!is.na(Type))
                debug(cp_rs)

                ##
                variables <- sort(colnames(cp_rs)[!colnames(cp_rs) %in% c("SampleID", "Type", "time", "status")])
                test_CP_Stat_report <- data.frame()
                key_df <- data.frame()

                if(!is.null(b.cp.variables)) {
                    variables = intersect(variables, b.cp.variables)
                    ifelse(b.signature.type == "Predictive", pts_type <- c(d_id, sort(unique(cp_rs$Type))), pts_type <- d_id)
                    
                    ##
                    include_continous = FALSE
                    test_CP_Stat.list = list()
                    key_df = data.frame()
                    for(type in pts_type) {
                        ifelse(type == d_id, cp_rs0 <- cp_rs, cp_rs0 <- cp_rs[cp_rs$Type == type, ])
                        clust_column <- ifelse(type == d_id, "Type", "Treatment_type")

                        Stat <- lapply(variables, function(x) {
                            rs0 <- cp_rs0[, c(clust_column, x)]
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
                                        pval0 <- ifelse(ft$p.value<0.01,
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
                                                    format(pval0, digits=2, scientific=T),
                                                    round(pval0, 3))
                                    null0[,3] <- round(mean(x0, na.rm = T), 3)
                                    null0[,4] <- round(mean(y0, na.rm = T), 3)
                                    rs2 <- data.frame(null0, P = pval0, check.names=F)
                                    rs2
                                }
                            }
                        })
                        res <- do.call(rbind.data.frame, Stat)

                        ##
                        if(nrow(res) > 0) {
                            if(b.signature.type == "Predictive" & type != pts_type[1]) {
                                colnames(res)[-(1:2)] <- paste0(type, "<br />", colnames(res)[-(1:2)])
                            }
                            test_CP_Stat.list <- c(test_CP_Stat.list, list(res))
                            key_df <- rbind(key_df, res[, c("Variable", "Level")])
                        }
                    }

                    ##
                    if(nrow(key_df) > 0) {
                        key_df <- unique(key_df)
                        ids <- sort(unique(key_df$Variable))
                        key_df <- do.call(rbind.data.frame, lapply(ids, function(x) {key_df[key_df$Variable == x, ] })) 
                        test_CP_Stat_report <- left_join.multi(key.df = key_df, ID = colnames(key_df), list.df = test_CP_Stat.list)
                        colnames(test_CP_Stat_report)[-(1:2)] <- paste0(colnames(test_CP_Stat_report)[-(1:2)], "<br />", d_id)
                    
                        ## subgroup p
                        if(b.signature.type == "Predictive") {
                            n_g = ifelse(b.group == "2 groups", 2, 3)
                            if(n_g == 2) {
                                s_ids <- sort(c(seq(6, 13, by = 3)[1:n_g], (seq(6, 13, by = 3)+1)[1:n_g]))
                            } else {
                                s_ids <- sort(c(seq(7, 15, by = 3)[1:n_g], (seq(7, 15, by = 3)+1)[1:n_g]))
                            }

                            if(ncol(test_CP_Stat_report) >= 10) {
                                p_sub <- apply(test_CP_Stat_report[, s_ids], 1, function(x) {
                                    x <- as.numeric(x)
                                    if(sum(is.na(x)) == 0) {
                                        p0 <- fisher.test(matrix(x, nrow = n_g))$p.value
                                        p0 <- ifelse(p0 < 0.01, format(p0, digits = 2, scientific = TRUE), round(p0, 3))
                                    } else { "" }
                                })
                                test_CP_Stat_report <- data.frame(test_CP_Stat_report, p_sub, check.names = F)
                                colnames(test_CP_Stat_report)[ncol(test_CP_Stat_report)] <- paste0("P-subgroup", "<br />", d_id)
                            }
                            test_CP_Stat_report[test_CP_Stat_report == "Treatment_type"] <- "Therapy"
                        }
                    }
                }

                CP_Stat.list[[z]] <- test_CP_Stat_report
                CP_Stat_key_df.list[[z]] <- key_df

                ##
                names(test_rs_plot.list)[z] <- d_id
                names(test_surv.list)[z] <- d_id
                names(test_pts_ann_sb.list)[z] <- d_id
                names(test_surv_prob_report.list)[z] <- d_id
                names(surv_prob_key_df.list)[z] <- d_id
                names(CP_Stat.list)[z] <- d_id
                z=z+1
            } else {
                valid_warning.list[[z0]] <- paste0("Not all signature variables were present or \n levels of some signature variables were not same as training dataset for ", d_id)
                names(valid_warning.list)[z0] <- d_id
                z0=z0+1
            }
            print(I)
        }

        ## merge surv prob tables
        if(length(surv_prob_key_df.list) > 0) {
            key_df <- unique(do.call(rbind, surv_prob_key_df.list))
            surv_p_key_df <- unique(rbind(train_surv_prob_report[, colnames(key_df)], key_df))
            surv_prob_report <- left_join.multi(key.df = surv_p_key_df,
                                                ID = colnames(surv_p_key_df),
                                                list.df = c(list(train_surv_prob_report),
                                                test_surv_prob_report.list))
        } else { 
            surv_prob_report <- train_surv_prob_report
        }


        ## merge CP stat tables
        CP_Stat_report = data.frame()
        if(!is.null(b.cp.variables)) {
            CP_Stat.list <- lapply(CP_Stat.list, function(x) {if(nrow(x) == 0) {NULL} else {x} })
            CP_Stat.list <- CP_Stat.list[!sapply(CP_Stat.list, is.null)]
            CP_Stat_report <- train_CP_Stat_report            
            if(length(CP_Stat_key_df.list) > 0) {
                CP_Stat_key_df <- unique(do.call(rbind, CP_Stat_key_df.list))
                if(nrow(CP_Stat_key_df) > 0) {
                    if(b.signature.type == "Predictive") {
                        CP_Stat_key_df[CP_Stat_key_df == "Treatment_type"] <- "Therapy"
                    }
                    cp_key_df <- unique(rbind(train_CP_Stat_report[,c("Variable", "Level")], CP_Stat_key_df))
                    ids <- sort(unique(cp_key_df$Variable))
                    cp_key_df <- do.call(rbind.data.frame, lapply(ids, function(x) { cp_key_df[cp_key_df$Variable == x, ] }))
                    CP_Stat_report <- left_join.multi(key.df = cp_key_df,
                                                      ID = colnames(cp_key_df),
                                                      list.df = c(list(train_CP_Stat_report),
                                                      CP_Stat.list)
                                                     )
                }
            }
            CP_Stat_report <- data.frame("Raw Rank" = 1:nrow(CP_Stat_report), CP_Stat_report, check.names=F)
        }
    } else {
            testsets.name = c()
            test_surv.list = list()
            test_pts_ann_sb.list = list()
            surv_prob_report = train_surv_prob_report
            CP_Stat_report = data.frame()
            test_rs_plot.list = list()
            valid_warning.list <- list()
            if(!is.null(b.cp.variables)) {
                CP_Stat_report <- train_CP_Stat_report            
                CP_Stat_report <- data.frame("Raw Rank" = 1:nrow(CP_Stat_report), CP_Stat_report, check.names=F)
            }
    }
    surv_prob_report <- surv_prob_report %>% select(-Key)
    surv_prob_report <- data.frame("Raw Rank" = 1:nrow(surv_prob_report), surv_prob_report, check.names=F)
    save(testsets.name, test_surv.list, test_pts_ann_sb.list,
         surv_prob_report, CP_Stat_report, test_rs_plot.list,
         valid_warning.list, file = rdata.filename)

    ##
    section <- 'Prediction'
    key <- 'external_evaluate_analysis_rdata'
    value <- paste0(paste0(outpath,'/external_evaluate/'), rdata.filename)
    system(paste("python3", paste0(script.dir,"/","write_ini.py"), user.config.ini.file, section, key, value))
    message_nosignature <- "Successful!"
} else {
    message_nosignature <- "There was no signature detected!"
}
print(message_nosignature)
# lsdata(rdata.filename)



