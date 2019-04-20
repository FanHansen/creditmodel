#'  Process NAs
#'
#' \code{process_nas_var} is the function for missing value analysis and process missing value by knn imputation or central impulation or random imputation.
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
#' @param default_miss Logical. If TRUE, assigning the missing values to -1 or "Unknown", otherwise ,processing the missing values according to the results of missing analysis.
#' @param parallel Logical, parallel computing. Default is FALSE.
#' @param note Logical, outputs info. Default is TRUE.
#' @param dir_path The path for periodically saved missing analysis file. Default is "./variable".
#' @param ... Other parameters.
#' @return A dat frame with no NAs.
#'
#' @examples
#' \dontrun{
#' dat_na = process_nas(dat = UCICreditCard, default_miss = FALSE, 
#' target = "default.payment.next.month", 
#' parallel = FALSE,ex_cols = "ID$" method = "median")
#' }
#'
#' @export

process_nas <- function(dat, x_list = NULL, default_miss = TRUE, class_var = FALSE,
parallel = FALSE, ex_cols = NULL, method = "median", note = TRUE, dir_path = "./variable",...) {
    if (note) cat("[NOTE]", "process nas ...\n")
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
            na_vars = na_vars, mat_nas_shadow = mat_nas_shadow, dt_nas_random = dt_nas_random)
        }

        dat[, na_vars] <- loop_function(func = process_nas_var, x_list = na_vars,
        args = list(dat = dat, nas_rate = nas_rate, mat_nas_shadow = mat_nas_shadow,
        dt_nas_random = dt_nas_random, missing_type = missing_type,
        method = method, default_miss = default_miss, note = note,
        dir_path = dir_path), bind = "cbind", as_list = FALSE, parallel = parallel)
    }
    return(dat)
}

#' Missing Analysis
#'
#' #' \code{analysis_nas} is the function for missing value analysis
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
        char_x_list = get_names(dat = dt_nas_random, types = c('character', 'factor'), ex_cols = NULL, get_ex = FALSE)
        num_x_list = get_names(dat = dt_nas_random, types = c('numeric', 'integer', 'double'), ex_cols = NULL, get_ex = FALSE)
        corr_random_num = corr_random_char = NULL
        if (length(num_x_list)>1) {
            corr_random_num = cor(dt_nas_random[num_x_list], method = "spearman")
        }
        if (class_var && length(char_x_list) > 3) {
            corr_random_char = char_cor(dat = dt_nas_random, x_list = char_x_list, parallel = FALSE, note = FALSE)
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
                    missing_type = ifelse(cor_self > 0.05 & cor_others_mean > 0.05, "IM", ifelse(cor_others_mean > 0.05, "MAR", "MCAR"))
                    rm(cor_self, corr_random_others, cor_others, corr_random_x_1, cor_others_mean)
                } else {
                    if (any(c('factor', 'character') == class(dat[, x])[1]) && !is.null(corr_random_char) &&
                    is.element(x, colnames(corr_random_char))) {
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
method = "median", note = TRUE, dir_path = "./variable",  ...) {
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
                dat[nas, x] = "Unknown"
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
        nas_analysis = data.frame(Feature = x, miss_rate = as_percent(nas_rate_x, digits = 6), miss_type = missing_type_x)
        save_dt(nas_analysis, dir_path = dir_path, file_name = "missing_analysis", append = TRUE, note = FALSE)
        if (note) cat(paste0(paste(unlist(nas_analysis), collapse = "\t"), "\n"))
        is_not_na <- unlist(dat[, x][!is.na(dat[, x])])
        if (nas_rate_x > 0 & len_num > 0) {
            if (!is.na(missing_type_x) && missing_type_x != "No_NAs") {
                if (all(c('numeric', 'integer', 'double', 'Date') != class(dat[, x])[1])) {
                    if ((missing_type_x == "IM" & (length(nas) > 30 | nas_rate_x > 0.01)) | nas_rate_x > 0.1) {
                        dat[nas, x] = "Unknown"
                    } else {
                        dat[, x] = knn_nas_imp(dat = dat, x, mat_nas_shadow = mat_nas_shadow, dt_nas_random = dt_nas_random)
                    }
                } else {
                    if ((missing_type_x == "IM" & (length(nas) > 30 | nas_rate_x > 0.01)) | nas_rate_x > 0.1) {
                        dat[nas, x] = -1
                    } else {
                        if ((missing_type_x == "MCAR" | length(nas) < 10 | nas_rate_x < 0.001) & length(unique(is_not_na)) > 10) {
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
        if (any(c("integer", "numeric", "double") == class(dat[, x][1]))) {
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
            if (length(mat_cor) >1) {
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
                            dat[i, x] =  get_median(dat[base::setdiff(1:n_row, miss_x), x][k_neighbors], exp(-na_dist[k_neighbors]))
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



#' get central value.
#'
#' This function is not intended to be used by end user. 
#'
#' @param  x  A vector or list.
#' @param  weight_avg  avarage weigth to calculate means.
#' @export
#' @importFrom stats aggregate

get_median <- function(x, weight_avg = NULL) {
    if (any(c('numeric', 'integer', 'double') == class(x)[1])) {
        if (is.null(weight_avg)) {
            central_value <- median(x, na.rm = T)
        } else {
            central_value <- ifelse(sum(weight_avg) > 0, sum(x * (weight_avg / sum(weight_avg))), NA)
        }
    } else {
        x <- as.factor(x)
        if (is.null(weight_avg)) {
            central_value <- levels(x)[which.max(table(x))]
        } else {
            central_value <- levels(x)[which.max(aggregate(weight_avg, list(x), sum)[, 2])]
        }
    }
    return(central_value)
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
                if (length(unique(is_not_na)) ==1 &&  unique(is_not_na) == 0) {
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
