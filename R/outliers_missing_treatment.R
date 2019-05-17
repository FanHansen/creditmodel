#'  Missing Treatment
#'
#' \code{process_nas_var} is for missing value analysis and treatment using knn imputation, central impulation and random imputation.
#' \code{process_nas} is a simpler wrapper for \code{process_nas_var}.
#'
#' @param dat A data.frame with independent variables.
#' @param x_list Names of independent variables.
#' @param x  The name of variable to process.
#' @param ex_cols A list of excluded variables. Regular expressions can also be used to match variable names. Default is NULL.
#' @param nas_rate A list contains nas rate of each variable.
#' @param mat_nas_shadow A shadow matrix of variables which contain nas.
#' @param dt_nas_random A data.frame with random nas imputation.
#' @param missing_type Type of missing, genereted by code{\link{analysis_nas}}
#' @param method  The methods of imputation by knn."median" is knn imputation by k neighbors median.
#' @param class_var Logical, nas analysis of the nominal variables. Default is TRUE.
#' @param default_miss Logical. If TRUE, assigning the missing values to -1 or "Missing", otherwise ,processing the missing values according to the results of missing analysis.
#' @param parallel Logical, parallel computing. Default is FALSE.
#' @param note Logical, outputs info. Default is TRUE.
#' @param save_data Logical. If TRUE, save missing analysis to \code{dir_path}
#' @param dir_path The path for periodically saved missing analysis file. Default is "./variable".
#' @param file_name The file name for periodically saved missing analysis file. Default is NULL.
#' @param ... Other parameters.
#' @return A dat frame with no NAs.
#'
#' @examples
#' dat_na = process_nas(dat = UCICreditCard[1:1000,], default_miss = FALSE,
#' target = "default.payment.next.month",
#' parallel = FALSE,ex_cols = "ID$" ,method = "median")
#'
#' @export

process_nas <- function(dat, x_list = NULL, default_miss = TRUE, class_var = FALSE,
                        parallel = FALSE, ex_cols = NULL, method = "median", note = FALSE,
                        save_data = TRUE, file_name = NULL, dir_path = tempdir(), ...) {
    if (note) cat("[NOTE]", "process nas ...\n")
    if (save_data) {
        dir_path = ifelse(!is.character(dir_path), tempdir(), dir_path)
        if (!dir.exists(dir_path)) dir.create(dir_path)
        if (!is.character(file_name)) file_name = NULL
        }
    mat_nas_shadow = nas_rate = na_vars = dt_nas_random = missing_type = NULL
    x_list = get_x_list(x_list = x_list, dat_train = dat, dat_test = NULL, ex_cols = ex_cols)
    mat_nas_shadow <- get_shadow_nas(dat[, x_list])
    if (length(mat_nas_shadow) > 0) {
        nas_rate = colSums(mat_nas_shadow) / nrow(mat_nas_shadow)
        na_vars = names(mat_nas_shadow)
        if (default_miss) {
            dt_nas_random = NULL
            missing_type = NULL
        } else {
            dt_nas_random = get_nas_random(dat)
            missing_type = analysis_nas(dat = dat, class_var = class_var, nas_rate = nas_rate,
                                        na_vars = na_vars, mat_nas_shadow = mat_nas_shadow,
                                        dt_nas_random = dt_nas_random)
        }

        dat[, na_vars] <- loop_function(func = process_nas_var, x_list = na_vars,
                                        args = list(dat = dat, nas_rate = nas_rate,
                                                    mat_nas_shadow = mat_nas_shadow,
                                                    dt_nas_random = dt_nas_random,
                                                    missing_type = missing_type,
                                                    method = method, default_miss = default_miss,
                                                    note = note, save_data = save_data,
                                                    file_name = file_name, dir_path = dir_path),
                                        bind = "cbind", as_list = FALSE, parallel = parallel)
    }
    return(dat)
}

#' Missing Analysis
#'
#' #' \code{analysis_nas} is for understanding the reason for missing data and understand distribution of missing data so we can categorise it as: 
#' \itemize{
#'   \item Missing completely at random(MCAR)
#'   \item MMissing at random(MAR), or
#'   \item Missing not at random, also known as IM. 
#' }
#' @param dat A data.frame with independent variables and target variable.
#' @param class_var Logical, nas analysis of the nominal variables. Default is TRUE.
#' @param na_vars  Names of variables which contain nas.
#' @param nas_rate A list contains nas rate of each variable.
#' @param mat_nas_shadow A shadow matrix of variables which contain nas.
#' @param dt_nas_random A data.frame with random nas imputation.
#' @param ... Other parameters.
#' @return A data.frame with outliers analysis for each variable.
#' @export

analysis_nas <- function(dat, class_var = FALSE, nas_rate = NULL, na_vars = NULL, mat_nas_shadow = NULL, dt_nas_random = NULL, ...) {
    if (is.null(mat_nas_shadow)) { mat_nas_shadow = get_shadow_nas(dat) }
    if (is.null(nas_rate)) { nas_rate = colSums(mat_nas_shadow) / nrow(mat_nas_shadow) }
    if (is.null(na_vars)) { na_vars = names(mat_nas_shadow) }
    if (is.null(dt_nas_random)) { dt_nas_random = get_nas_random(dat) }
    missing_type_var = NULL
    if (length(mat_nas_shadow) > 0) {
        char_x_list = get_names(dat = dt_nas_random, types = c('character', 'factor'),
                                ex_cols = NULL, get_ex = FALSE)
        num_x_list = get_names(dat = dt_nas_random, types = c('numeric', 'integer', 'double'),
                               ex_cols = NULL, get_ex = FALSE)
        corr_random_num = corr_random_char = NULL
        if (length(num_x_list) > 1) {
            corr_random_num = cor(dt_nas_random[num_x_list], method = "spearman")
        }
        if (class_var && length(char_x_list) > 3) {
            corr_random_char = char_cor(dat = dt_nas_random, x_list = char_x_list,
                                        parallel = FALSE, note = FALSE)
        }
        mising_type_var = vapply(na_vars, function(x) {
            nas = which(is.na(dat[, x]))
            nas_rate_x = nas_rate[x]
            if (nas_rate_x == 0 | !is.element(x, names(mat_nas_shadow))) {
                missing_type = "No_Nas"
            } else {
                if (any(c('numeric', 'integer', 'double') == class(dat[, x])[1]) &&
                    !is.null(corr_random_num) && is.element(x, names(corr_random_num))) {
                    #correlation of nas
                    #correlation of nas and other variables .
                    corr_random_x = corr_random_num[, x]
                    corr_random_x_1 = corr_random_x[!is.na(corr_random_x)]
                    if (length(corr_random_x_1) > 5) {
                        max_cor = order(abs(corr_random_x_1), decreasing = TRUE)[1:5]
                    } else {
                        max_cor = order(abs(corr_random_x_1), decreasing = TRUE)[1:length(abs(corr_random_x_1))]
                    }
                    cor_names = names(corr_random_x_1[max_cor])
                    corr_random_others = cor(dt_nas_random[, cor_names], y = mat_nas_shadow[, x],
                    method = "spearman", use = "pairwise.complete.obs")

                    cor_self = corr_random_others[which(rownames(corr_random_others) == x)]
                    if (is.na(cor_self)) {
                        cor_self = 0
                    }
                    cor_others = corr_random_others[which(rownames(corr_random_others) != x)]
                    if (length(!is.na(cor_others)) > 5) {
                        cor_others_mean = mean(abs(cor_others[order(abs(cor_others), decreasing = TRUE)[1:5]]), na.rm = TRUE)
                    } else {
                        cor_others_mean = mean(abs(cor_others[order(abs(cor_others),
                        decreasing = TRUE)[1:length(!is.na(abs(cor_others)))]]), na.rm = TRUE)
                    }

                    # misssing type of variale.
                    missing_type = ifelse(cor_self > 0.05 & cor_others_mean > 0.05, "IM",
                                          ifelse(cor_others_mean > 0.05, "MAR", "MCAR"))
                    rm(cor_self, corr_random_others, cor_others, corr_random_x_1, cor_others_mean)
                } else {
                    if (any(c('factor', 'character') == class(dat[, x])[1]) &&
                        !is.null(corr_random_char) && is.element(x, colnames(corr_random_char))) {
                        corr_random_x = corr_random_char[, x]
                        corr_random_x_1 = corr_random_x[!is.na(corr_random_x)]
                        if (length(corr_random_x_1) > 6) {
                            max_cor = order(abs(corr_random_x_1), decreasing = TRUE)[1:6]
                        } else {
                            max_cor = order(abs(corr_random_x_1), decreasing = TRUE)[1:length(abs(corr_random_x_1))]
                        }
                        cor_names = names(corr_random_x_1[max_cor])
                        corr_random_others = char_cor_vars(dt_nas_random[cor_names], mat_nas_shadow[x])
                        cor_others = corr_random_others[-1]
                        cor_self = corr_random_others[1]
                        if (length(!is.na(cor_others)) > 5) {
                            cor_others_mean = mean(abs(cor_others[1:5]), na.rm = TRUE)
                        } else {
                            cor_others_mean = mean(abs(cor_others[1:length(!is.na(abs(cor_others)))]), na.rm = TRUE)
                        }
                        if (is.na(cor_self)) {
                            cor_self = 0
                        }
                        # misssing type of variale.
                        missing_type = ifelse(cor_self > 0.03 & cor_others_mean > 0.03, "IM",
                        ifelse(cor_others_mean > 0.03, "MAR", "MCAR"))
                        rm(cor_self, corr_random_others, cor_others, corr_random_x_1, cor_others_mean)
                    } else {
                        missing_type = "IM"
                    }
                }
            }
            missing_type
        }, FUN.VALUE = character("1"))

    }
    return(mising_type_var)
}


#' @rdname process_nas
#' @export

process_nas_var <- function(dat = dat, x, default_miss = TRUE, nas_rate = NULL,
                            mat_nas_shadow = NULL, dt_nas_random = NULL, missing_type = NULL,
                            method = "median", note = FALSE, save_data = TRUE,
                            file_name = NULL, dir_path = tempdir(), ...) {
    nas <- which(is.na(dat[, x]))
    if (is.null(nas_rate)) {
        nas_rate_x = length(nas) / length(dat[, x])
    } else {
        nas_rate_x = nas_rate[x]
    }
    len_num = length(unique(dat[, x]))
    if (default_miss) {
        if (nas_rate_x > 0 & len_num > 0) {
            if (all(c('numeric', 'integer', 'double', 'Date') != class(dat[, x])[1])) {
                dat[nas, x] = "Missing"
            } else {
                dat[nas, x] = -1
            }
        }
    } else {
        if (is.null(missing_type)) {
            missing_type = analysis_nas(dat, mat_nas_shadow = mat_nas_shadow, dt_nas_random = dt_nas_random)
        }
        missing_type_x = "No_NAs"
        if (length(missing_type) > 0 && is.element(x, names(missing_type))) {
            missing_type_x <- missing_type[x]
        }
        nas_analysis = data.frame(Feature = x, miss_rate = as_percent(nas_rate_x, digits = 6),
                                  miss_type = missing_type_x)
        if (save_data) {
            save_dt(nas_analysis, dir_path = dir_path,
                  file_name = ifelse(is.null(file_name), "missing.analysis", paste(file_name, "missing.analysis", sep = ".")),
                  append = TRUE, note = FALSE)
        }
        if (note) cat(paste0(paste(unlist(nas_analysis), collapse = "\t"), "\n"))
        is_not_na <- unlist(dat[, x][!is.na(dat[, x])])
        if (nas_rate_x > 0 & len_num > 0) {
            if (!is.na(missing_type_x) && missing_type_x != "No_NAs") {
                if (!is.element(class(dat[, x])[1], c('numeric', 'integer', 'double', 'Date'))) {
                    if ((missing_type_x == "IM" & (length(nas) > 30 | nas_rate_x > 0.01)) |
                        nas_rate_x > 0.1) {
                        dat[nas, x] = "Missing"
                    } else {
                        dat[, x] = knn_nas_imp(dat = dat, x, mat_nas_shadow = mat_nas_shadow,
                                               dt_nas_random = dt_nas_random)
                    }
                } else {
                    if ((missing_type_x == "IM" & (length(nas) > 30 | nas_rate_x > 0.01)) | nas_rate_x > 0.1) {
                        dat[nas, x] = -1
                    } else {
                        if ((missing_type_x == "MCAR" | length(nas) < 10 | nas_rate_x < 0.001) &
                            length(unique(is_not_na)) > 10) {
                            set.seed(46)
                            dat[nas, x] = sample(is_not_na, size = length(nas), replace = TRUE)
                        } else {
                            dat[, x] = knn_nas_imp(dat = dat, x, mat_nas_shadow = mat_nas_shadow,
                            dt_nas_random = dt_nas_random, method = method)
                        }
                    }
                }
            }
        }
    }
    return(dat[x])
}




#'  Imputate nas using KNN
#'
#' This function is not intended to be used by end user.
#' @param dat A data.frame with independent variables.
#' @param x  The name of variable to process.
#' @param nas_rate A list contains nas rate of each variable.
#' @param mat_nas_shadow A shadow matrix of variables which contain nas.
#' @param dt_nas_random A data.frame with random nas imputation.
#' @param k  Number of neighbors of each obs which x is missing.
#' @param scale Logical.Standardization of variable.
#' @param method The methods of imputation by knn."median" is knn imputation by k neighbors median.
#' @export

knn_nas_imp <- function(dat, x, nas_rate = NULL, mat_nas_shadow = NULL,
dt_nas_random = NULL, k = 10, scale = FALSE, method = 'median') {
    dat = checking_data(dat)
    miss_x <- which(!complete.cases(dat[, x]))
    if (length(miss_x) > 0) {
        n_row <- nrow(dat)
        if (any(is.element(class(dat[, x]), c("integer", "numeric", "double")))) {
            if (is.null(mat_nas_shadow)) {
                mat_nas_shadow = get_shadow_nas(dat)
            }
            if (is.null(dt_nas_random)) {
                dt_nas_random <- get_nas_random(dat)
            }
            mat_cor = cor_names = nums_names = NULL
            nums_names = get_names(dat = dt_nas_random, types = c('numeric', 'integer', 'double'),
            get_ex = FALSE, ex_cols = c(x))
            mat_cor = cor(dt_nas_random[, nums_names], y = mat_nas_shadow[, x], method = "spearman")
            if (length(mat_cor) > 1) {
                mat_cor_sub = mat_cor[which(!is.na(mat_cor)),]
                if (length(mat_cor_sub) > 5) {
                    max_cor = order(abs(mat_cor_sub), decreasing = TRUE)[1:5]
                } else {
                    max_cor = order(abs(mat_cor_sub), decreasing = TRUE)[1:length(!is.na(abs(mat_cor_sub)))]
                }
                cor_names = names(mat_cor_sub[max_cor])
                dm = data.frame(dat[x], dt_nas_random[, cor_names])
                dm = as.matrix(dm)
                no_ms <- which(!is.na(dm[, x]))
                set.seed(46)
                no_miss_x <- sample(no_ms, round(n_row / 3))
                no_miss_obs <- dm[no_miss_x,]
                na_distance = na_dist = NULL
                if (nrow(no_miss_obs) > k) {
                    dat[miss_x, x] <- vapply(miss_x, function(i) {
                        na_distance <- scale(no_miss_obs[, - 1], dm[i, - 1], FALSE)
                        na_distance <- ifelse(na_distance > 0, 1, na_distance)
                        na_dist <- sqrt(drop(na_distance ^ 2 %*% rep(1, ncol(na_distance))))
                        k_neighbors <- order(na_dist)[seq(k)]
                        if (method == "median") {
                            dat[i, x] = get_median(dat[base::setdiff(1:n_row, miss_x), x][k_neighbors])
                        } else {
                            dat[i, x] = get_median(dat[base::setdiff(1:n_row, miss_x), x][k_neighbors],
                                                   exp(-na_dist[k_neighbors]))
                        }
                    }, FUN.VALUE = numeric(1)
                    )
                } else {
                    dat[miss_x, x] = -1
                }
            } else {
                dat[miss_x, x] = -1
            }
        } else {
            dat[miss_x, x] = get_median(dat[base::setdiff(1:n_row, miss_x), x])
        }
    }
    rm(na_dist, na_distance, dm, mat_cor, mat_cor_sub, no_miss_obs)
    return(dat[x])
}


#' get_shadow_nas
#'
#' This function is not intended to be used by end user.
#'
#' @param  dat  A data.frame contained only predict variables.
#' @export
get_shadow_nas <- function(dat) {
    nas_shadow <- as.data.frame(abs(is.na(dat)))
    mat_nas_shadow = list()
    if (length(nas_shadow) > 0) {
        low_variance = vapply(nas_shadow,
        function(x) max(table(x, useNA = "always")) / nrow(nas_shadow),
        FUN.VALUE = numeric(1))
        mat_nas_shadow = nas_shadow[which(low_variance < 1)]
    }
    return(mat_nas_shadow)
}



#' get_nas_random
#'
#' This function is not intended to be used by end user.
#'
#' @param  dat  A data.frame contained only predict variables.
#' @export

get_nas_random <- function(dat) {
    dt_nas_random = c()
    if (length(dat) > 0) {
        dt_nas_random <- quick_as_df(lapply(dat, function(x) {
            nas <- which(is.na(x))
            is_not_na <- unlist(x[!is.na(x)])
            if (length(unique(is_not_na)) <= 1) {
                if (length(unique(is_not_na)) == 1 && unique(is_not_na) == 0) {
                    x[nas] = 1
                } else {
                    x[nas] = -1
                }
            } else {
                x[nas] = sample(is_not_na, size = length(nas), replace = TRUE)
            }
            x
        }))
        low_variance = vapply(dt_nas_random,
        function(x) max(table(x, useNA = "always")) / nrow(dt_nas_random),
        FUN.VALUE = numeric(1))
        dt_nas_random = dt_nas_random[which(low_variance < 1)]
    }
    return(dt_nas_random)
}



#' Outliers Treatment
#'
#'
#' \code{outliers_kmeans_lof} is for outliers detection and treatment using Kmeans and Local Outlier Factor (lof)
#' \code{process_outliers} is a simpler wrapper for \code{outliers_kmeans_lof}.
#' @param dat Dataset with independent variables and target variable.
#' @param target The name of target variable.
#' @param x_list Names of independent variables.
#' @param x  The name of variable to process.
#' @param ex_cols A list of excluded variables. Regular expressions can also be used to match variable names. Default is NULL.
#' @param kc Number of  clustering centers for Kmeans
#' @param kn Number of neighbors for LOF.
#' @param parallel Logical, parallel computing.
#' @param note Logical, outputs info. Default is TRUE.
#' @param process Logical, process outliers, not just analysis.
#' @param save_data Logical. If TRUE, save outliers analysis file to the specified folder at \code{dir_path}
#' @param dir_path The path for periodically saved outliers analysis file. Default is "./variable".
#' @param file_name The file name for periodically saved outliers analysis file. Default is NULL.
#'
#' @return A data frame with outliers process to all the variables.
#'
#' @examples
#' dat_out = process_outliers(UCICreditCard[1:10000,],
#'                           target = "default.payment.next.month",
#'                           ex_cols = "date$", kc = 3, kn = 10, parallel = FALSE)
#' @export




process_outliers <- function(dat, target, ex_cols = NULL, kc = 3, kn = 5, x_list = NULL,
                             parallel = FALSE,note = FALSE, process = TRUE,
                             save_data = TRUE, file_name = NULL, dir_path = tempdir()) {
    if (note) cat("[NOTE]", "process outliers using kmeans and lof...\n")
    if (save_data) {
        dir_path = ifelse(!is.character(dir_path), tempdir(), dir_path)
        if (!dir.exists(dir_path)) dir.create(dir_path)
        if (!is.character(file_name)) file_name = NULL
        }
    x_list = get_x_list(x_list = x_list, dat_train = dat, dat_test = NULL, ex_cols = ex_cols)
    if (length(x_list) > 0) {
        num_x_list = get_names(dat = dat[x_list], types = c('numeric', 'integer', 'double'),
                           ex_cols = c(target, ex_cols), get_ex = FALSE)
    } else {
        stop("No variable in the x_list or ex_col excludes all variables.\n ")
    }

    if (length(num_x_list) > 0) {
        dat[ ,num_x_list] = loop_function(func = outliers_kmeans_lof, x_list = num_x_list,
                                          args = list(dat = dat, target = target, kc = kc, kn = kn,
                                                      note = note, process = process,
                                                      save_data = save_data, file_name = file_name,
                                                      dir_path = dir_path),
                                          bind = "cbind", as_list = FALSE, parallel = parallel)
    }
    return(dat)
}


#' @rdname process_outliers
#' @export

outliers_kmeans_lof <- function(dat, x, target = NULL, kc = 3, kn = 5,
                                note = FALSE, process = TRUE,
                                save_data = TRUE, file_name = NULL, dir_path = tempdir()) {
    lof <- outliers_detection(dat = dat, x)
    len_num <- length(unique(dat[, x]))
    if (length(lof) > 0) {
        outlier_type = analysis_outliers(dat = dat, x, target = target, lof = lof)
    } else {
        outlier_type ="no_outlier"
    }
    out_analysis = data.frame(Feature = x,
                              out_rate = as_percent(length(lof) / nrow(dat), digits = 6),
                              out_type = outlier_type)

    if(save_data){
      save_dt(out_analysis, dir_path = dir_path,
              file_name = ifelse(is.null(file_name), "outliers.analysis", paste(file_name, "outliers.analysis", sep = ".")),
              append = TRUE, note = FALSE)
    }
    if (note) cat(paste(unlist(out_analysis), collapse = "\t"), "\n")
    if (process) {
        if (len_num > 0) {
            if (outlier_type =="no_outlier") {
                dat[, x] <- dat[, x]
            } else {
                if (outlier_type == "random") {
                    dat[lof, x] <- sample(dat[, x][!is.na(dat[, x])], size = 1, replace = TRUE)
                } else {
                    if (outlier_type == "NAs") {
                        dat[lof, x] <- NA
                    } else {
                        top_rate = length(lof) / length(dat[, x])
                        if (outlier_type == "floor") {
                            p1 <- ifelse(any(c('numeric', 'integer', 'double') == class(dat[, x])),
                                         quantile(dat[, x], 0.005, na.rm = TRUE, type = 3),
                                         get_median(dat[, x]))
                            dat[lof, x] <- p1
                        } else {
                            if (outlier_type == "top") {
                                p99 <- ifelse(any(c('numeric', 'integer', 'double') == class(dat[, x])),
                                              quantile(dat[, x], 0.995, na.rm = TRUE, type = 3),
                                              get_median(dat[, x]))
                                dat[lof, x] <- p99
                            }
                        }
                    }
                }
            }
        }
    }
    return(dat[x])
}



#' Outliers Analysis
#'
#' #' \code{analysis_outliers} is the function for outliers analysis.
#' @param dat A data.frame with independent variables and target variable.
#' @param target The name of target variable.
#' @param x  The name of variable to process.
#' @param lof  Outliers of each variable detected by \code{outliers_detection}.
#' @return A data.frame with outliers analysis for each variable.
#' @export

analysis_outliers <- function(dat, target, x, lof = NULL) {
    if (is.null(lof)) {
        lof = outliers_detection(dat = dat, x)
    }
    dm_x = dat[[x]]
    outlier_type ="no_outlier"
    if (length(lof) > 10) {
        corr_lof_y = corr_lof_na = corr_lof_x = corr_xy = 1
        is_lof <- abs(dm_x %in% dm_x[lof])
        if (any(is.element(class(dm_x) ,  c('Date', 'numeric', 'integer', 'double')))) {
            corr_lof_x = cor(dm_x, y = is_lof, method = "spearman", use = "pairwise.complete.obs")
        } else {
            table_x_lof = table(dm_x, is_lof)
            corr_lof_x = sqrt(try(chisq.test(table_x_lof, correct = T,
                                             simulate.p.value = TRUE)$statistic,
                                  silent = TRUE) / sum(table_x_lof))
        }
        ctb_na = ctb_y = matrix(0, ncol = 2, nrow = 2)
        if (!is.null(target)) {
            corr_lof_y = cor(dat[, target], y = is_lof, method = "spearman")
            if (any(is.element(class(dm_x), c('Date', 'numeric', 'integer', 'double')))) {
                corr_xy = cor(dat[, target],
                              y = dm_x,
                              method = "spearman",
                              use = "pairwise.complete.obs")
            } else {
                table_xy = table(dm_x, dat[, target])
                corr_xy = sqrt(try(chisq.test(table_xy, correct = TRUE,
                                              simulate.p.value = TRUE)$statistic,
                                   silent = TRUE) / sum(table_xy))
            }
            ctb_y <- as.matrix(table(is_lof, dat[, target]))
        }
        NAs <- which(is.na(dm_x))
        if (length(NAs) > 30) {
            is_na <- abs(is.na(dm_x))
            corr_lof_na = cor(is_na, y = is_lof, method = "spearman")
            ctb_na <- as.matrix(table(is_lof, is_na))
        }
        ctb_list = list(ctb_y, ctb_na)
        effect_sz = odds_ratio = dif_1_rate = c()
        for (j in 1:2) {
            dt_gb = ctb_list[[j]]
            if (any(dt_gb[1, ] <= 10, na.rm = TRUE) | any(dt_gb[2, ] <= 10, na.rm = TRUE)) {
                effect_sz[j] = 0
                odds_ratio[j] = 1
                dif_1_rate[j] = 0
            } else {
                effect_sz[j] = sqrt(try(chisq.test(dt_gb, correct = TRUE,
                                                   simulate.p.value = TRUE)$statistic,
                                        silent = TRUE) / sum(dt_gb))
                odds_ratio[j] = (dt_gb[1, 1] * dt_gb[2, 2]) / (dt_gb[2, 1] * dt_gb[1, 2])
                dif_1_rate[j] = dt_gb[1, 2] / (dt_gb[1, 1] + dt_gb[1, 2]) - dt_gb[2, 2] / (dt_gb[2, 1] + dt_gb[2, 2])
            }
        }
        if (abs(corr_lof_y) < 0.01 & abs(corr_lof_x) < 0.01 & effect_sz[1] < 0.05 &
            abs(odds_ratio[1] - 1) < 0.1 & abs(dif_1_rate[1]) < 0.01) {
            outlier_type = "random"
        } else {
            if (abs(corr_lof_na) < 0.01 & effect_sz[2] < 0.1 & abs(odds_ratio[2] - 1) < 0.05 &
                abs(dif_1_rate[2]) < 0.01) {
                outlier_type = "NAs"
            } else {
                if ((corr_lof_y < 0 & corr_xy < 0) | (corr_lof_y > 0 & corr_xy > 0)) {
                    outlier_type = "top"
                } else {
                    if ((corr_lof_y < 0 & corr_xy > 0) | (corr_lof_y > 0 & corr_xy < 0)) {
                        outlier_type = "floor"
                    }
                }
            }
        }
    } else {
        if (length(lof) > 0) {
            outlier_type = "random"
        } else {
            outlier_type ="no_outlier"
        }
    }
    return(outlier_type)
}

#' Outliers Detection
#' \code{outliers_detection} is for outliers detecting using Kmeans and Local Outlier Factor (lof)
#' @param dat A data.frame with independent variables.
#' @param x  The name of variable to process.
#' @param kc Number of  clustering centers for Kmeans
#' @param kn Number of neighbors for LOF.
#' @return  Outliers of each variable.
#' @export

outliers_detection <- function(dat, x, kc = 3, kn = 5) {
    dm_x = dat[, x]
    len = length(unique(dm_x))
    lof = NULL
    if (len > 50 & any(c('Date', 'numeric', 'integer', 'double') == class(dm_x)[1])) {
        dm_s = scale(dm_x)
        dm_s[is.na(dm_s)] <- get_median(as.numeric(dm_s))
        #  packages("mclust")
        set.seed(46)
        #G = if (len > 100) { c(1:min(ceiling(len / 50), 5)) } else { c(1:2) }
        #d_clust = Mclust(as.matrix(dm_s), G = G, verbose = FALSE)
        kc = ifelse(len > 100, kc, 2)
        km = kmeans(dm_s, centers = kc)
        #distance of samples
        center <- list()
        distance <- list()
        for (i in 1:kc) {
            center[[i]] = matrix(km$centers[i,], nrow = length(dm_s), ncol = 1, byrow = T)
            distance[[i]] = sqrt(rowSums((dm_s - center[[i]]) ^ 2))
        }
        distance_matrix = Reduce("cbind", distance)
        min_dist = apply(distance_matrix, 1, mean)
        lof_sub = which(min_dist > quantile(min_dist, 0.999, na.rm = TRUE, type = 3))
        max_out = c()
        for (i in 1:length(km$centers)) { max_out[i] = length(which(km$cluster == i)) }
        lof = base :: intersect(which(km$cluster %in% km$cluster[lof_sub]),
                                which(km$cluster %in% which(max_out / sum(max_out) < 0.03)))
        rm(lof_sub, dm_s, km, center, distance_matrix, min_dist, max_out)
        if (length(lof) / length(dm_x) > 0.05) {
            lof = NULL
        } else {
            if (length(lof) / length(dm_x) <= 0.001 | length(lof)< kn + 3 ) {
                lof = lof
            } else {
                outlier_scores = local_outlier_factor(dat = dm_x[lof], k = kn )
                if (length(lof[which(outlier_scores > 1.2)]) > 1) {
                    lof = lof[which(outlier_scores > 1.2)]
                    rm(outlier_scores)
                } else {
                    lof = NULL
                }
            }
        }
    } else {
        if (len > 2 & len < 50 & any(c('factor', 'character') == class(dm_x)[1])) {
            dm_x <- as.character(dm_x)
            lof = which(dm_x %in% names(which(table(dm_x, useNA = "no") / length(dm_x) < 0.005)))
            rm(dm_x)
        }
    }
    return(lof)
}
