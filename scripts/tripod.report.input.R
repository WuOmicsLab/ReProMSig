#
# @Copyright: Peking University Cancer Hospital, All Rights Reserved.
# @Author: Lihua Cao
# @Date: 2022-01
# @LastEditTime: 2023-04-16
# @LastEditors: Lihua Cao
# @Description: Generate the variables, tables and graphs for a signature that will be shown in the reporting file.
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
	if(length(args)!=1) {
		stop("Usage: Rscript tripod.report.input.R [user.config.ini.file]\n")
	}
	user.config.ini.file <- args[1]
} else {
	user.config.ini.file <- "ColoGuide_Stage_II_local/output/sig.ini"
}

# 2) Library
library(dplyr)
library(tibble)
library(survival)
library(rms)
library(pec)

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
model_analysis_rdata.file <- get_value.ccb(config_file = user.config.ini.file, key = 'model_analysis_rdata')[[1]]
independence_analysis_rdata.file <- get_value.ccb(config_file = user.config.ini.file, key = 'independence_analysis_rdata')[[1]]
performance_analysis_rdata.file <- get_value.ccb(config_file = user.config.ini.file, key = 'performance_analysis_rdata')[[1]]
external_evaluate_analysis_rdata.file <- get_value.ccb(config_file = user.config.ini.file, key = 'external_evaluate_analysis_rdata')[[1]]

# 4) define/source functions
source(paste0(script.dir,'ccb.helper.R'))


# readin and preprocess  -------------------------------------------------------------------------
load(model_analysis_rdata.file)
# cgwtools::lsdata(model_analysis_rdata.file)
tripod.ini.file = paste0(outpath,'/','tripod.ini')
mysetwd(paste0(outpath,'/tripod/'))
coef1
rdata.filename = paste0(b.name.sig, '.tripod.RData')

## plot 
if(nrow(coef1) > 0) {
    load(independence_analysis_rdata.file)
    load(performance_analysis_rdata.file)
    load(external_evaluate_analysis_rdata.file)
    time_n <- length(b.times)

    # Items of results  ----------------------------------------------------------------------------------------------
    # Item_13b ---------
    section <- 'Item_13b'
    if(b.have.profiles != "No") {
        png(file=paste0("figure_exp_distr_train_validsets.png"), width = 12*300, height=6*300, res=300)
        par(mfrow=c(ceiling((length(b.profile_data.list)-1)/3)+1,3))
        names(b.profile_data.list)[1] <- "Training dataset"
        sapply(names(b.profile_data.list), function(x) {
            exp0 <- b.profile_data.list[[x]]
            if(x=="Training dataset") {
                density.1group(exp0, main=names(b.profile_data.list)[1], xlab='Log2 abundance', ylab='Density', cex=0.8,cex.axis=1)
                plot(1:10,1:10,type='n',xaxt='n',yaxt='n',xlab='',ylab='',bty='n')
                plot(1:10,1:10,type='n',xaxt='n',yaxt='n',xlab='',ylab='',bty='n')
            } else {
                density.1group(exp0, border_col=rgb(59,125,196,max = 255),
                               polygon_col=rgb(59,125,196, 88,max = 255),
                               main=x, xlab='Log2 abundance', ylab='Density',
                               cex=0.8, cex.axis=1)
            }
        })
        dev.off()

        ##
        key <- "figure_exp_distr_train_validsets"
        value <- paste0(outpath, '/tripod/figure_exp_distr_train_validsets.png')
        system(paste("python3", paste0(script.dir,"/","write_ini.py"), tripod.ini.file, section, key, value))
    }

    # Item_14b ---------
    if(nrow(uni_multi_cox_report) > 0) {
        if(b.signature.type == "Predictive") {
            colnames(uni_multi_cox_report)[which(colnames(uni_multi_cox_report) == "Value<br />(Control, Treatment)")] <- "Value"
        }
        ##
        variables <- unique(uni_multi_cox_report$Variable)
        uni_multi_cox_report <- do.call(rbind, lapply(variables, function(x) {
            df1 <- uni_multi_cox_report %>% filter(Variable == x, Value != "Unknown")
            df2 <- uni_multi_cox_report %>% filter(Variable == x, Value == "Unknown")
            rbind(df1, df2)
        }))
        uni_multi_cox_report$'Raw Rank' <- (1:nrow(uni_multi_cox_report))

        if(b.signature.type == "Predictive") {
            colnames(uni_multi_cox_report)[which(colnames(uni_multi_cox_report) == "Value")] <- "Value<br />(Control, Treatment)"
        }
        ##
        if(b.signature.type != "Predictive") {
            ids <- grep("P<br />", colnames(uni_multi_cox_report))
            ids <- do.call(c, lapply(1:nrow(uni_multi_cox_report), function(i) {
                x0 <- as.character(uni_multi_cox_report[i, ids])
                x0 <- x0[!is.na(x0) & x0!=""]
                if(length(x0) == 0) { i }
            }))

            ##
            uni_multi_cox_report[ids, "Value"] <- do.call(c, lapply(ids, function(i) {
                ref0 <- strsplit(uni_multi_cox_report[i, "Value"], "\\(vs. ")[[1]][2]
                gsub(")$", " (ref.)", ref0)
            }))
        }
    } else {
        uni_multi_cox_report <- uni_cox
    }

    ##
    uni_multi_cox_report[uni_multi_cox_report == "Signature score"] <- "Signature group"
    save(uni_multi_cox_report, file = paste0(b.name.sig, '.uni_multi_cox.RData'))
    write.file(uni_multi_cox_report, file = paste0(b.name.sig, '.uni_multi_cox.csv'))

    ##
    section <- 'Item_14b'
    keys <- c("figure_uni_cox", "figure_multi_cox", "table_uni_multi_cox")
    he <- ifelse(nrow(uni_cox) > 12, ceiling(nrow(uni_cox)*0.5), ceiling(nrow(uni_cox)*0.7))
    
    ifelse(b.signature.type == "Predictive",
           labeltext <- c("Value\n(Control, Treatment)", "Events/N\n", "Hazard Ratio (95% CI)\n", "P\n", "P_interaction\n"), 
           labeltext <- c("Value", "Events/N", "Hazard Ratio (95% CI)", "P"))

    ##
    uni_cox[is.na(uni_cox)] <- ""
    uni_cox[,1] <- gsub("^RS\\ \\(", "Signature_group\\ \\(", uni_cox[,1])
    ##
    png(file='figure_uni_cox_forest.png', width = 12*300, height = he*300, res=300)
    par(mar=c(1,1,1,1))
    forestplot.func(signature_type = b.signature.type,
                    cox_table = uni_cox, labeltext = labeltext,
                    continuous_variables = continuous_annvar)
    dev.off()

    ##
    if(nrow(multi_cox) > 0) {
        multi_cox[is.na(multi_cox)] <- ""
        multi_cox[,1] <- gsub("^RS\\ \\(", "Signature_group\\ \\(", multi_cox[,1])
        png(file='figure_multi_cox_forest.png', width = 12*300, height = he*300, res=300)
        par(mar=c(1,1,1,1))
        forestplot.func(signature_type = b.signature.type,
                        cox_table = multi_cox, labeltext = labeltext,
                        continuous_variables = continuous_annvar)
        dev.off()
    }

    ##
    values <- c(paste0(outpath, '/tripod/figure_uni_cox_forest.png'), 
                paste0(outpath, '/tripod/figure_multi_cox_forest.png'),
                paste0(outpath, "/tripod/", b.name.sig, '.uni_multi_cox.RData'))
               
    for(i in 1:length(keys)) {
        system(paste("python3",paste0(script.dir,"/","write_ini.py"),tripod.ini.file,section,keys[i],values[i]))
    }

    ## variables
    if(nrow(uni_multi_cox_report) > 0) {
        uni_cox_variables <- sort(unique(uni_multi_cox_report[which(uni_multi_cox_report$'P<br />Univariate analysis' != ""), "Variable"]))
        uni_cox_variables <- uni_cox_variables[!uni_cox_variables %in%
                                                c("Signature score", "Signature group", "Signature_score", "Signature_group", "Predictive_group")
                                              ]
        multi_cox_variables <- sort(unique(uni_multi_cox_report[which(uni_multi_cox_report$'P<br />Multivariate analysis' != ""), "Variable"]))
        multi_cox_variables <- multi_cox_variables[!multi_cox_variables %in%
                                                    c("Signature score", "Signature group", "Signature_score", "Signature_group", "Predictive_group")
                                                  ]
    } else {
        uni_cox_variables <- c()
        multi_cox_variables <- c()
    }
    save(uni_cox_variables, multi_cox_variables, file = rdata.filename)

    # Item_15a ---------
    section <- 'Item_15a'
    write.file(coef_report, file="coef.txt")
    
    ## risk score formula
    coef_report0 <- coef_report
    if(length(grep("\\(vs\\.", coef_report$Characteristic)) > 0) {
        coef_report0 <- coef_report[-grep("\\(vs\\.", coef_report$Characteristic), ]
    }
    mode(coef_report0$Coefficient) <- "numeric"

    formula0 <- paste0(round(coef_report0$Coefficient, 3), "*<b>", coef_report0$Characteristic, "</b>")
    formula0 <- paste(formula0, collapse = " + ")
    formula0 <- gsub(" \\+ -", " - ", formula0)
    risk_score_formula <- paste0("Signature score = ", formula0)
    
    ## mfp transformed variables
    if(b.have.profiles == "No") {
        x <- table(coef_report$Variable)
        mfp_continuous_variables <- sort(unique(names(x[x==1])))
    } else {
        mfp_continuous_variables <- c()
    }
    
    # signature exp corplot 
    if(b.have.profiles == "Yes" & nrow(coef1) > 1) {
        signature.corr <- cor(Train[,coef1[,1]])
        res1 <- corrplot::cor.mtest(Train[,coef1[,1]], conf.level = 0.95)
        col3 <- colorRampPalette(c('DodgerBlue3','white', "OrangeRed"))

        ##
        addCoef.col = "black"
        number.cex=0.8
        tl.cex=1
        cl.cex = 1
        cl.offset = 0.5
        if(ncol(signature.corr) >= 20){
            addCoef.col = NULL
            number.cex = 0.5
            tl.cex = 0.8
            cl.cex = 1
            cl.offset = 1
        }
        ##
        size0 <- ifelse(ncol(signature.corr) >= 50, 10, 8)
        png(file='figure_sig_corrplot.png', width = size0*300, height=size0*300, res=300)
        corrplot::corrplot(signature.corr,method="color",type='lower', col=col3(20), cl.pos = "b", bg=col3(20),
                           number.cex=number.cex, tl.cex=tl.cex, cl.cex=cl.cex, tl.srt=45, addCoef.col = addCoef.col, tl.col = "black",
                           number.font=1, cl.lim=c(-1,1), cl.length=length(seq(-1,1,by=0.2)), cl.offset = cl.offset, mar=c(0,2,2,5)
                           )
        dev.off()
    }
    
    if(b.have.profiles == "Yes" & nrow(coef1) > 1) {
        keys <- c("table_coef", "figure_sig_corrplot")
        values <- c(paste0(outpath, "/tripod/", "coef.txt"),
                    paste0(outpath, "/tripod/", "figure_sig_corrplot.png"))
    } else {
        keys <- "table_coef"
        values <- paste0(outpath, "/tripod/", "coef.txt")
    }

    for(i in 1:length(keys)) {
        system(paste("python3", paste0(script.dir,"/","write_ini.py"), tripod.ini.file, section, keys[i], values[i]))
    }


    # Item_15b ---------
    # calibrate lines for signature model ---
    section = 'Item_15b'
    if(b.have.profiles == "Yes") {
        figure_calibrate_lines_signature <- paste0(outpath, "/tripod/", 'figure_calibrate_lines_signature.png')
        ##
        df <- train_multi_sig_nom
        sams <- intersect(risk_score$SampleID, rownames(df))
        df <- df[sams, ]
        id0 <- which(colnames(df) %in% c("RS", "Signature_score", "Predictive_group"))
        df[,id0] <- as.character(df[,id0])
        df[,id0] <- risk_score[risk_score[,1] %in% sams, "RS"]
        df <- df[, colnames(df) %in% c("time", "status", "Signature_score", "Predictive_group")]
        ##
        variable0 <- colnames(df)[3]
        df0 <- data.frame(df[,variable0], check.names = F)
        colnames(df0)[1] <- variable0
        dd <<- datadist(df0)
        options(datadist="dd")
        formula0 <- paste("Surv(time, status==1)", "~", paste(variable0, collapse="+"))

        ##
        u.list <- b.times * 30.41667
        names(u.list) <- paste0(b.times, "-months")
        
        ##
        png(file='figure_calibrate_lines_signature.png', width = 4*time_n*300, height = 4*300, res=300)
        par(mfrow=c(1, time_n))
        sapply(names(u.list), function(x0) {
            u <- as.numeric(u.list[x0])
            cph_cbl <- cph(as.formula(formula0), data = df, x = TRUE, y = TRUE, surv = TRUE, na.action=na.delete, time.inc = u)

            if(length(cph_cbl) > 1) {
                m <- ifelse(nrow(df) > 200, floor(nrow(df)/5), floor(nrow(df)/3))
                if(nrow(df) < 50) { m <- floor(nrow(df)/2) }
                
                if(nrow(df[df$time < u,]) > 0) {
                    cal0 <- calibrate(cph_cbl, cmethod='KM', method="boot", u = u, m = m, B = 1000)
                    par(mar=c(5,5,2,2))
                    plot(cal0, lwd=3, lty=1,
                         errbar.col = "#0076C0",
                         xlim=c(0,1), ylim=c(0,1), 
                         xlab="", ylab="", col = "#C06253",
                         subtitles=F, cex.axis=1)
                    lines(cal0[, c('mean.predicted', "KM")], type = 'b', lwd = 3, pch = 16, col = "#2166AC")
                    mtext(paste0("Predicted ", x0, " ", b.endpoint, " (%)"),side=1,line=3,cex=1)
                    mtext(paste0("Observed ", x0, " ", b.endpoint, " (%)"), side=2, line=3, cex=1)
                    abline(0, 1, lty=3, lwd=3, col = "#0076C0")
                    box(lwd = 1)
                }
            }
        })
        dev.off()
        
        ##
        key <- "figure_calibrate_lines_signature"
        value <- figure_calibrate_lines_signature
        system(paste("python3", paste0(script.dir,"/","write_ini.py"), tripod.ini.file, section, key, value))
    }
    
    # Constructed nomgram combing significant variables ------
    variables <- colnames(train_multi_sig_nom)
    variables <- variables[!variables %in% c("SampleID", "time", "status")]
    ##
    nomogram_variables <- variables[!variables %in% c("Signature_score", "Predictive_group")]
    nomogram_variables <- sort(unique(nomogram_variables))

    ##
    variables <- c(variables[variables %in% c("Signature_score", "Signature_group", "Predictive_group")],
                   variables[!variables %in% c("Signature_score", "Signature_group", "Predictive_group")])

    ##
    if(b.signature.type == "Predictive") {
        train_multi_sig_nom.list <- list(all_pts = train_multi_sig_nom, untreated = train_multi_sig_nom_c, treated = train_multi_sig_nom_t)
    } else {
        train_multi_sig_nom.list <- list(all_pts = train_multi_sig_nom)
    }
    
    ##
    if(length(variables) >= 2) {
        lapply(names(train_multi_sig_nom.list), function(x) {
            df <- train_multi_sig_nom.list[[x]]
            if(b.signature.type == "Predictive" & b.have.profiles == "Yes") {
                df <- df %>% mutate(Signature_group = as.character(Predictive_group)) %>% select(-Predictive_group)
                df[df=="No-benefit"] <- "Low risk"
                df[df=="Benefit"] <- "High risk"
            }
            colnames(df)[colnames(df) == "Signature_score"] <- "Signature_group"
            variables[variables %in% c("Signature_score", "Predictive_group")] <- "Signature_group"             

            ## 1)data
            dd <<- datadist(df[,variables])
            options(datadist="dd")

            ## 2) Cox modelling
            formula0 <- paste("Surv(time, status==1)", "~", paste(variables, collapse = "+"))
            cph0 <- cph(as.formula(formula0), data =df, x=TRUE, y=TRUE, surv=TRUE, na.action=na.delete)

            ## 3) Nomogram generation
            if("nomogram" %in% b.tool.cons) {
                file_nom <- ifelse(b.signature.type != "Predictive", "figure_nomogram.png", paste0("figure_nomogram_", x,".png"))
                width0 <- ifelse(length(variables) >= 10, 12*1.5, 12*1.2)
                height0 <- ifelse(length(variables) >= 10, 8*1.5, 8*1.2)
                png(file=file_nom, width = width0*300, height = height0*300, res=300)
                par(mar=c(1,1,1,1))
                ##
                if(length(cph0) > 1) {
                    med <- Quantile(cph0) # median survival time
                    surv <- Survival(cph0) # model survial model
                    ## 
                    funlabel <- paste0(b.times, "-months", " ", b.endpoint)
                    b.times_day <- b.times * 30.41667
                    fun <- lapply(b.times_day, function(n0) {
                        function(x, day = n0) surv(day, x)
                    })
                    ##
                    nom <- nomogram(cph0, fun = fun, funlabel = funlabel, lp = FALSE, maxscale = 100)
                    plot(nom,
                        #lplabel="Linear score",
                        xfrac=0.5, varname.label=TRUE, varname.label.sep="=", ia.space=0.2, 
                        tck=NA, tcl=-0.20, lmgp=0.3,
                        points.label='Points', total.points.label='Total Points',
                        total.sep.page=FALSE, cap.labels=FALSE,
                        cex.var = 1.6, cex.axis = 1.05, lwd=5,
                        label.every = 1, col.grid = gray(c(0.8, 0.95)))
                    }
                dev.off()
                
                ## 1/3/5 years calibrate lines of nomogram ---
                u.list <- b.times * 30.41667
                names(u.list) <- paste0(b.times, "-months")
                ##
                file_cal <- ifelse(b.signature.type != "Predictive", "figure_calibrate_lines.png", paste0("figure_calibrate_lines_", x,".png"))
                png(file=file_cal, width = 4*time_n*300, height = 4*300, res=300)
                par(mfrow=c(1, time_n))
                sapply(names(u.list), function(x0) {
                    u <- as.numeric(u.list[x0])
                    cph_cbl <- cph(as.formula(formula0), data = df, x = TRUE, y = TRUE, surv = TRUE, na.action=na.delete, time.inc = u)
                    if(length(cph_cbl) > 1) {
                        m <- ifelse(nrow(df) > 200, floor(nrow(df)/5), floor(nrow(df)/3))
                        if(nrow(df) < 50) { m <- floor(nrow(df)/2) }

                        if(nrow(df[df$time < u,]) > 0) {
                            cal0 <- calibrate(cph_cbl, cmethod="KM", method="boot", u = u, m = m, B = 1000) 
                            par(mar=c(5,5,2,2))
                            plot(cal0, lwd=3, lty=1, errbar.col = "#0076C0", xlim=c(0,1), ylim=c(0,1), 
                                 xlab="", ylab="", col = "#C06253", subtitles=F, cex.axis=1)

                            lines(cal0[, c('mean.predicted', "KM")], type = 'b', lwd = 3, pch = 16, col = "#2166AC")
                            mtext(paste0("Predicted ", x0, " ", b.endpoint, " (%)"),side=1,line=3,cex=1)
                            mtext(paste0("Observed ", x0, " ", b.endpoint, " (%)"), side=2, line=3, cex=1)
                            abline(0, 1, lty=3, lwd=3, col = "#0076C0")
                            box(lwd = 1)
                        }
                        print(x0)
                    }
                })
                dev.off()

                ##
                figure_nomogram <- paste0(outpath, "/tripod/", file_nom)
                figure_calibrate_lines <- paste0(outpath, "/tripod/", file_cal)
            } else {
                figure_nomogram=""
                figure_calibrate_lines=""
            }

            ## 
            if(b.signature.type == "Predictive") {
                keys <- c("figure_nomogram_exist",
                          paste0("figure_nomogram_", x),
                          paste0("figure_calibrate_lines_", x)
                         )
            } else {
                keys <- c("figure_nomogram_exist", "figure_nomogram", "figure_calibrate_lines")
            }
            ##
            values <- c(ifelse("nomogram" %in% b.tool.cons, TRUE, FALSE),
                        figure_nomogram,
                        figure_calibrate_lines
                        )
            for(i in 1:length(keys)) {
                system(paste("python3",paste0(script.dir,"/","write_ini.py"),tripod.ini.file,section,keys[i],values[i]))
            }
        })
    }



    # Item_16 --------------------------------------
    section = 'Item_16'
    value1 = "xxxxxx"
    # ROC method KM using maker RS ---
    color1 <- c("#34621d","#d34b00","#dd8e0a","#e6c710","#9aa440","#347e86","#304e98","#816ab3",
                "#C0C000", "#3B7DC4", "#F0C430", '#0072B5FF', 'blueviolet', "#AFABAB", "#709770",
                "#6FAC47","#71E945","#CCE744","green4","#DCA8A1", "orange", "#6BDFDC", "#23a09d",
                "#6BDFDC", "#D4D4D6", "#4C72B0", "#C44E52", "#DE6591","#E4CF4D", "#DEE2AF", "#D66644",
                "#835DD5","#C9ADD4","#7EE5B3", "#78B9DA","#A433EA","#D4AB67","#B6DB79","#D28ADA",
                "#E253CB","#7288D5")
    color2 <- 'darkred'
    color <- rep('x', (ncol(LP)-2))
    names(color) <- colnames(LP)[-(1:2)]

    id0 <- match(colnames(LP)[-(1:2)], c("RS", "Signature_score", "Signature_group", "Predictive_group"))
    color[which(!is.na(id0))] <- color2
    color[which(is.na(id0))] <- color1[1:(ncol(LP)-3)]
    
    ##
    ROC_variables <- colnames(LP)[-(1:2)]
    ROC_variables <- sort(unique(ROC_variables[!ROC_variables %in% c("time", "status", "Combined", "Signature_score", "Signature_group")]))

    ##
    width0 = 6*time_n
    png(file = paste0("figure_roc.png"), width = width0*300, height = 6*300, res=300)
    par(mfrow=c(1,time_n))
    ##
    cutoffs <- b.times * 30.41667
    names(cutoffs) <- paste0(b.times, "-months")

    for(type0 in names(cutoffs)) {
        cutoff <- cutoffs[type0]
        ID <- type0
        ##
        ROC <- list()
        z=1
        for(x in colnames(LP)[-(1:2)]) {
            LP0 <- LP[,c("time", "status", x)]
            LP0 <- LP0[complete.cases(LP0),]
            ROC[[z]] <- survivalROC::survivalROC(Stime=LP0$time,
                                                status=LP0$status, marker = LP0[,x],
                                                predict.time = cutoff,
                                                method="KM")
            z=z+1
        }
        names(ROC) <- colnames(LP)[-(1:2)]

        ## order by AUC
        auc <- unlist(lapply(names(ROC), function(X) {
        round(ROC[[X]]$AUC,3)
        }))
        names(auc) <- names(ROC)
        auc <- rev(sort(auc))

        ##
        par(mar=c(5,5,3,2))
        i0=1
        plot(ROC[[i0]]$FP, ROC[[i0]]$TP, type="l",
             col=color[i0], xlim=c(0,1), ylim=c(0,1),
             xlab="", ylab="", main="", sub=' ', lwd=3, lty=1,
             cex.axis=1.5, bty='n', cex.main=1.5)
        abline(0, 1, col="gray",lty=3)
        main0 <- paste0("Training dataset", " (", type0,")")
        mtext('1 - Specificity',side=1,line= 3,cex=1.5)
        mtext('Sensitivity',side=2,line=3,cex=1.5)
        mtext(main0, side=3,line=1,cex=1.5)
        if(length(ROC)>=2) {
            for(i0 in 2:length(ROC)) {
                lines(ROC[[i0]]$FP, ROC[[i0]]$TP, type="l",col=color[i0],xlim=c(0,1), ylim=c(0,1),lwd=3,lty=1)
            }
        }
        leg <- legend('bottomright', rep('', (ncol(LP)-2)), x.intersp=1, y.intersp=0.8,
                    col = rep('white',(ncol(LP)-2)), bty = "n", seg.len=1,cex=1, 
                    lty = rep(1,(ncol(LP)-2)), lwd=3)

        leg <- legend(x = unique(leg$text$x)-0.5,y = leg$text$y[1]+0.05, legend = names(auc), 
                    x.intersp=1, y.intersp=0.8, lty = rep(1,(ncol(LP)-2)), lwd=3,
                    col=color[names(auc)], bty = "n", seg.len=1,cex=1)
        ##
        x0 = unique(leg$text$x)
        y0 = leg$text$y
        text(x=x0-0.15,y=y0,labels=auc,cex=1)
        text(x=x0 - 0.15,y=y0[1]+0.05,labels='AUC',cex=1)
        text(x=x0 + 0.07,y=y0[1]+0.05,labels="Characteristic",cex=1)
        segments(x0=x0-0.2, y0 = y0[1]+0.025, x1=x0+0.22, y1 = y0[1]+0.025, lwd = 1.5)
    }
    dev.off()


    # Prediction error curves ---
    Model <- models
    data <- train.roc[,c('time', 'status', roc_variables[!roc_variables %in% "Combined"])]

    if("Combined" %in% colnames(LP)) {
        sams <- intersect(rownames(data), rownames(LP))
        data <- data.frame(data[sams, ], Combined = LP[sams, "Combined"], check.names = F)
    }

    # remove samples with NA for pec modeling analysis
    PE_na_sams_reason <- ""
    PE_na_sams_num <- nrow(data[!complete.cases(data),])
    data <- data[complete.cases(data),]
    perror <- pec(Model, data=data,
                  formula = Surv(time,status==1)~1,
                  keep.matrix, splitMethod='cvK',
                  na.action = na.pass, verbose=T)
    # 
    crps0 <- crps(perror, times = max(b.times * 30.41667), start=0)[,1] # max(data$time)
    crps0 <- round(c(sort(crps0[!names(crps0) %in% "Reference"]), crps0[names(crps0) %in% "Reference"]), 3)
    color_pe <- c(Reference = "black", color)

    ##
    png(file='figure_pe.png', width = 6*300, height = 6*300, res=300)
    par(mar=c(5,5,2,2))
    plot(perror, smooth=T, add.refline=T,
         col=color_pe, lwd=3, legend=F,
         xlab='', ylab='', cex.axis=1.5,
         lty=c(3, rep(1,(ncol(LP)-2))),
         xlim=c(0, max(b.times * 30.41667))
         )

    mtext('Time (days)', side=1, line=2.5, cex=1.5)
    mtext('Predition error', side=2, line=2.5, cex=1.5)

    leg <- legend('bottomright', legend = names(crps0),
                  x.intersp=1, y.intersp=0.8,
                  lty=c(rep(1,(ncol(LP)-2)), 3),
                  col=color_pe[names(crps0)],
                  bty="n", seg.len=1, cex=1, lwd=3)


    ##
    x0 <- unique(leg$text$x)
    y0 <- leg$text$y
    text(x=x0-x0/7, y=y0+0.005, labels=crps0,cex=1, adj=c(1,1))
    text(x=x0-x0/11.5, y=y0[1]+0.025, labels='IBScore',cex=1, adj=c(1,1))
    text(x=x0+x0/3.06, y=y0[1]+0.025, labels="Characteristic",cex=1, adj=c(1,1))
    segments(x0=x0-x0/3.5, y0=y0[1]+0.01,
             x1=x0+x0/3.3, y1=y0[1]+0.01,
             lwd=1.5)
    dev.off()
    
    ##
    PE_variables <- sort(unique(colnames(data)[!colnames(data) %in%c("time", "status", "Combined", "Signature_score", "Signature_group")]))
    x0 <- ifelse(length(PE_variables) > 1, "predictors", "predictor")
    x1 <- paste(PE_variables, collapse = ", ")

    if(PE_na_sams_num > 0) {
        PE_na_sams_reason <- paste0(PE_na_sams_num, " samples with missing values in ", x0, " (", x1, ") were excluded in the PE analysis")
    }

    # KM plot split by risk score on trainset ---
    if(b.signature.type == "Predictive") {
        lev <- sort(unique(train_stratify$Type))
        train_lev_rs_plot.list <- lapply(lev, function(x) {
            df <- train_stratify[train_stratify$Type == x, ]
            res <- df$Treatment_type
            names(res) <- df$SampleID
            res
        })
        names(train_lev_rs_plot.list) <- lev
        x <- train_stratify$Treatment_type
        names(x) <- train_stratify$SampleID
        train_rs_plot.list <- c(list("Training dataset" = x), train_lev_rs_plot.list)
    } else {
        x <- train_stratify$Type
        names(x) <- train_stratify$SampleID
        train_rs_plot.list <- list("Training dataset" = x)
    }
    surv.train <- unique(rbind(train_surv, train_surv_sg))[,-1]

    # input matrix for KM plots split by risk score on trainset and validation datasets ------------    
    train_km_input <- do.call(rbind.data.frame,
                        lapply(names(train_rs_plot.list), function(x) {
       
        rs0 <- train_rs_plot.list[[x]]
        rs0 <- data.frame(SampleID = names(rs0), Risk_group = rs0, Dataset = x, subgroup = "all_pts")

        surv0 <- surv.train %>% rownames_to_column(var = "SampleID")
        out <- merge.data.frame(surv0, rs0, by = "SampleID") %>% 
                        filter(!is.na(time) & !is.na(status) & !is.na(Risk_group))

        # for each subgroup---
        if(nrow(pts_ann_sb) > 1) {
            sub_levs <- sort(unique(pts_ann_sb[, b.subgroup.variable]))   
            for(lev0 in sub_levs) {
                sb_sams <- unique(pts_ann_sb[pts_ann_sb[,2] == lev0, "SampleID"])
                if(n_distinct(sb_sams) >= 3) {
                    sams <- intersect(out$SampleID, sb_sams)
                    out <- rbind(out, out[out$SampleID %in% sams,] %>%
                                    mutate(subgroup = paste0(b.subgroup.variable, "_", lev0)))
                    out$subgroup <- gsub("/", "_",out$subgroup)
                }
            }
        }
        out
    }))

    train_km_input.list <- split(train_km_input, train_km_input$subgroup)
    names(train_km_input.list) <- paste0("Training_dataset__msig__", names(train_km_input.list))
    names(train_km_input.list)[grep("all_pts", names(train_km_input.list))] <- "Training_dataset"
    km_input.list <- list("Training_dataset" = train_km_input.list)

    # input matrix for KM plots split by risk score on each validset
    figure_km_validsets_exist <- FALSE
    if(length(test_rs_plot.list) >= 1) {
        valid_km_input.list = list()
        z=1
        for(I in 1:length(test_rs_plot.list)) {
            surv0 <- test_surv.list[[I]] %>% rownames_to_column(var = "SampleID")
            df <- test_rs_plot.list[[I]]
            table(df[,2])

            if(b.signature.type == "Predictive") {
                test_lev_rs_plot.list <- split(df, df$Type)
                test_rs_plot0.list <- c(list(df), test_lev_rs_plot.list)
            } else {
                test_rs_plot0.list <- list(df)
            }
            d_id <- testsets.name[I]
            names(test_rs_plot0.list)[1] <- d_id
            test_pts_ann_sb <- test_pts_ann_sb.list[[d_id]]

            # input matrix for KM plots for each validset
            valid_km_input <- do.call(rbind.data.frame, lapply(names(test_rs_plot0.list), function(x) {
                df0 <- test_rs_plot0.list[[x]]
                clust_column <- ifelse(b.signature.type == "Predictive", "Treatment_type", "Type")
                rs0 <- df0[,clust_column]
                names(rs0) <- df0$SampleID
                rs0 <- data.frame(SampleID = names(rs0), Risk_group = rs0, Dataset = x, subgroup = "all_pts")
                out <- merge.data.frame(surv0, rs0, by = "SampleID") %>% filter(!is.na(time) & !is.na(status) & !is.na(Risk_group))

                # for each subgroup---
                if(nrow(test_pts_ann_sb) > 1) {
                    sub_levs <- sort(unique(test_pts_ann_sb[, b.subgroup.variable]))   
                    for(lev0 in sub_levs) {
                        sb_sams <- unique(test_pts_ann_sb[test_pts_ann_sb[,2] == lev0, "SampleID"])
                        if(n_distinct(sb_sams) >= 3) {
                            sams <- intersect(out$SampleID, sb_sams)
                            out <- rbind(out, out[out$SampleID %in% sams,] %>% mutate(subgroup = paste0(b.subgroup.variable, "_", lev0)))
                            out$subgroup <- gsub("/", "_",out$subgroup)
                        }
                    }
                }
                out
            }))

            ##
            valid_km_input.list0 <- split(valid_km_input, valid_km_input$subgroup)
            names(valid_km_input.list0) <- paste0("Validation_dataset_", I, "__msig__", names(valid_km_input.list0))
            names(valid_km_input.list0)[grep("all_pts", names(valid_km_input.list0))] <- paste0("Validation_dataset_", I)
            valid_km_input.list[[z]] <- valid_km_input.list0
            z=z+1
        }
        names(valid_km_input.list) <- paste0("Validation_dataset_", 1:length((test_rs_plot.list)))
        km_input.list <- c(km_input.list, valid_km_input.list)
        figure_km_validsets_exist <- TRUE
    }

    ##
    ids <- do.call(c, lapply(names(km_input.list), function(x) {
        df <- km_input.list[[x]][[1]]
        if(all(table(df$Risk_group) >= 5) & n_distinct(df$Risk_group) > 1) { x }
    }))
    km_input.list <- km_input.list[names(km_input.list) %in% ids]

    # KM ploting ---
    keys_km = c()
    values_km = c()
    for(x in names(km_input.list)) {
        list0 <- km_input.list[[x]]
        sbs <- names(list0)
        ##
        for(y in sbs) {
            df0 <- list0[[y]]
            if(all(table(df0$Risk_group) >= 5) & n_distinct(df0$Risk_group) >= 2) {
                y0 <- strsplit(y, "__msig__")[[1]][2]
                y0 <- gsub("/", "_", y0)
                file0 <- paste0("figure_km_", x, ".", y0, ".png")
                file0 <- gsub("\\.NA\\.", "\\.", file0)
                print(n_distinct(df0$Risk_group))
                msig_km.plot(df = df0, b.group = b.group,
                             b.signature.type = b.signature.type,
                             b.endpoint = b.endpoint,
                             png.file = file0)

                key0 <- gsub(".png", "", file0)
                keys_km <- c(keys_km, gsub("\\.", "_", key0))
                values_km <- c(values_km, paste0(outpath, "/tripod/", file0))
            }
        }
    }

    list.files(getwd(), pattern = "figure_km_")
    valid_KM_prex.list <- as.list(paste0("figure_km_", names(km_input.list)))
    names(valid_KM_prex.list) <- names(km_input.list)

    # estimated x-year survival ---
    surv_prob_report[is.na(surv_prob_report)] <- "-"
    save(surv_prob_report, file = paste0(b.name.sig, '.surv_prob.RData'))
    write.file(surv_prob_report, file = paste0(b.name.sig, '.surv_prob.csv'))
    
    # Clinicopathological association table ---
    if(!is.null(b.cp.variables)) {
        variables <- unique(CP_Stat_report$Variable)
        CP_Stat_report <- do.call(rbind, lapply(variables, function(x) {
            df1 <- CP_Stat_report %>% filter(Variable == x, Level != "Unknown")
            df2 <- CP_Stat_report %>% filter(Variable == x, Level == "Unknown")
            rbind(df1, df2)
        }))
        CP_Stat_report <- CP_Stat_report %>% filter(Level != "Unknown")
        CP_Stat_report$'Raw Rank' <- (1:nrow(CP_Stat_report))

        ##
        ids <- grep("P<br />", colnames(CP_Stat_report))
        ids <- do.call(c, lapply(1:nrow(CP_Stat_report), function(i) {
            x0 <- as.character(CP_Stat_report[i, ids])
            x0 <- x0[!is.na(x0) & x0!=""]
            if(length(x0) > 0) { i }
        }))
        CP_Stat_report[ids, "Level"] <- ""
    }
    save(CP_Stat_report, file = paste0(b.name.sig, '.cp_stat.RData'))
    write.file(CP_Stat_report, file = paste0(b.name.sig, '.cp_stat.csv'))

    
    ## write into tripod.ini.file ---
    keys = c("figure_roc", "figure_pe",
             "figure_km_validsets_exist",
              keys_km,
             "table_surv_prob_rdata", "table_cp_stat_list_rdata")

    values <- c(paste0(outpath, "/tripod/", 'figure_roc.png'), 
                paste0(outpath, "/tripod/", 'figure_pe.png'), 
                figure_km_validsets_exist, values_km,
                paste0(outpath, "/tripod/", b.name.sig, '.surv_prob.RData'),
                paste0(outpath, "/tripod/", b.name.sig, '.cp_stat.RData')
                )
                
    for(i in 1:length(keys)) {
        system(paste("python3",paste0(script.dir,"/","write_ini.py"), tripod.ini.file,section,keys[i],values[i]))
    }

    ## save key data in tripod Rdata
    cgwtools::resave(risk_score_formula, mfp_continuous_variables, nomogram_variables, 
                     ROC_variables, PE_variables, PE_na_sams_num, PE_na_sams_reason,
                     km_input.list, valid_KM_prex.list, 
                     train_multi_sig_nom.list,
                     file = rdata.filename)

    # write tripod.RData into ini file --------------------------------------
    section <- 'variables'
    key <- "user_tripod_rdata"
    value <- paste0(outpath, "/tripod/", b.name.sig, ".tripod.RData")
    system(paste("python3",paste0(script.dir,"/","write_ini.py"),tripod.ini.file, section, key, value))
    message_nosignature <- "Successful!"
} else {
    message_nosignature <- "There was no signature detected!"
}
print(message_nosignature)
# lsdata(rdata.filename)

