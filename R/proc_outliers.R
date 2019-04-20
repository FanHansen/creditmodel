#' Process Outliers
#'
#'
#' \code{outliers_kmeans_lof} Outliers processing using Kmeans and Local Outlier Factor (lof)
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
#' @param dir_path The path for periodically saved outliers analysis file. Default is "./variable".
#'
#' @return A data frame with outliers process to all the variables.
#'
#' @examples
#' \dontrun{
#' # load germancredit dat
#' dat(UCI_Credit_Card)
#' dat_out = process_outliers(dat_clean, 
#' target = "default.payment.next.month",
#' ex_cols = "date$", kc = 3, kn = 10, parallel = FALSE)
#' }
#'
#' @export



process_outliers <- function(dat, target, ex_cols = NULL, kc = 3, kn = 5, x_list = NULL, parallel = FALSE,
note = TRUE, process = TRUE, dir_path = "./variable") {
    if (note) cat("[NOTE]", "process outliers using kmeans and lof...\n")
    x_list = get_x_list(x_list = x_list, dat_train = dat, dat_test = NULL, ex_cols = ex_cols)
    num_x_list = get_names(dat = dat[x_list], types = c('numeric', 'integer', 'double'), ex_cols = c(target, ex_cols), get_ex = FALSE)
    if (length(num_x_list) > 0) {
        dat[, num_x_list] = loop_function(func = outliers_kmeans_lof, x_list = num_x_list, args = list(dat = dat, target = target, kc = kc, kn = kn, note = note, process = process, dir_path = dir_path), bind = "cbind", as_list = FALSE, parallel = parallel)
    }
    return(dat)
}


#' @rdname process_outliers
#' @export

outliers_kmeans_lof <- function(dat, x, target = NULL, kc = 3, kn = 5, note = TRUE, process = TRUE, dir_path = "./variable") {
    lof <- outliers_detection(dat = dat, x)
    len_num <- length(unique(dat[, x]))
    if (length(lof) > 0) {
        outlier_type = analysis_outliers(dat = dat, x, target = target, lof = lof)
    } else {
        outlier_type ="no_outlier"
    }
    out_analysis = data.frame(Feature = x, out_rate = as_percent(length(lof) / nrow(dat), digits = 6), out_type = outlier_type)
    save_dt(out_analysis, dir_path = dir_path , file_name = "outliers_analysis", append = TRUE, note = FALSE)
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
                            p1 <- ifelse(any(c('numeric', 'integer', 'double') == class(dat[, x])), quantile(dat[, x], 0.005, na.rm = TRUE, type = 3), get_median(dat[, x]))
                            dat[lof, x] <- p1
                        } else {
                            if (outlier_type == "top") {
                                p99 <- ifelse(any(c('numeric', 'integer', 'double') == class(dat[, x])), quantile(dat[, x], 0.995, na.rm = TRUE, type = 3), get_median(dat[, x]))
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
        if (any(c('Date', 'numeric', 'integer', 'double') == class(dm_x)[1])) {
            corr_lof_x = cor(dm_x, y = is_lof, method = "spearman", use = "pairwise.complete.obs")
        } else {
            table_x_lof = table(dm_x, is_lof)
            corr_lof_x = sqrt(try(chisq.test(table_x_lof, correct = T, simulate.p.value = TRUE)$statistic, silent = TRUE) / sum(table_x_lof))
        }
        ctb_na = ctb_y = matrix(0, ncol = 2, nrow = 2)
        if (!is.null(target)) {
            corr_lof_y = cor(dat[, target], y = is_lof, method = "spearman")
            if (any(c('Date', 'numeric', 'integer', 'double') == class(dm_x)[1])) {
                corr_xy = cor(dat[, target], y = dm_x, method = "spearman", use = "pairwise.complete.obs")
            } else {
                table_xy = table(dm_x, dat[, target])
                corr_xy = sqrt(try(chisq.test(table_xy, correct = T, simulate.p.value = TRUE)$statistic, silent = TRUE) / sum(table_xy))
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
        effect_sz = odds_ratio = dif_bad_rate = c()
        for (j in 1:2) {
            dt_gb = ctb_list[[j]]
            if (any(dt_gb[1, ] <= 10, na.rm = TRUE) | any(dt_gb[2, ] <= 10, na.rm = TRUE)) {
                effect_sz[j] = 0
                odds_ratio[j] = 1
                dif_bad_rate[j] = 0
            } else {
                effect_sz[j] = sqrt(try(chisq.test(dt_gb, correct = T, simulate.p.value = TRUE)$statistic, silent = TRUE) / sum(dt_gb))
                odds_ratio[j] = (dt_gb[1, 1] * dt_gb[2, 2]) / (dt_gb[2, 1] * dt_gb[1, 2])
                dif_bad_rate[j] = dt_gb[1, 2] / (dt_gb[1, 1] + dt_gb[1, 2]) - dt_gb[2, 2] / (dt_gb[2, 1] + dt_gb[2, 2])
            }
        }
        if (abs(corr_lof_y) < 0.01 & abs(corr_lof_x) < 0.01 & effect_sz[1] < 0.05 & abs(odds_ratio[1] - 1) < 0.1 & abs(dif_bad_rate[1]) < 0.01) {
            outlier_type = "random"
        } else {
            if (abs(corr_lof_na) < 0.01 & effect_sz[2] < 0.1 & abs(odds_ratio[2] - 1) < 0.05 & abs(dif_bad_rate[2]) < 0.01) {
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
        set.seed(2019)
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
        lof = base :: intersect(which(km$cluster %in% km$cluster[lof_sub]), which(km$cluster %in% which(max_out / sum(max_out) < 0.03)))
        rm(lof_sub, dm_s, km, center, distance_matrix, min_dist, max_out)
        if (length(lof) / length(dm_x) > 0.05) {
            lof = NULL
        } else {
            if (length(lof) / length(dm_x) <= 0.001) {
                lof = lof
            } else {
                outlier_scores = local_outlier_factor(dm_x[lof], k = kn)
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

#' local_outlier_factor
#' \code{local_outlier_factor}  is function for calculating the lof factor for a data set using knn
#' This function is not intended to be used by end user. 
#'
#' @param  dat  A data.frame contained only predict variables.
#' @param k Number of neighbors for LOF.
#' @export

local_outlier_factor <- function(dat, k) {
    dat <- as.matrix(dat)
    distance_data <- k_distance(dat, k)
    p <- dim(distance_data)[2L]
    local_reach_density_data <- local_reachability_density(distance_data, k)
    lof <- rep(0, p)
    # computer the local outlier factor of each observation in dat
    lof <- sapply(1:p, function(i) {
        n_neighbor <- distance_data[2, i] - distance_data[1, i] + 1
        j <- seq(0, (n_neighbor - 1))
        local_factor <- sum(local_reach_density_data[distance_data[3 + j, i]] / local_reach_density_data[i]) / n_neighbor
        lof[i] <- local_factor
    }
    )
    lof
}

#' k_distance
#' \code{k_distance} is for returning an object in which each column
#' This function is not intended to be used by end user. 
#'
#' @param  dataset  A data.frame contained only predict variables.
#' @param  neighbors  A data.frame contained neighbors of each obs.
#' @export

k_distance <- function(dataset, neighbors) {
    row_num <- dim(dataset)[1L]
    matrix_nrow <- neighbors * 2 + 2
    knn_distance <- matrix(0, nrow = matrix_nrow, ncol = row_num)
    for (i in 1:row_num) {
        neighdist <- k_distance_neighbor(dataset[i,], dataset, neighbors)
        x <- length(neighdist)
        if (x > matrix_nrow) {
            knn_distance <- rbind(knn_distance, matrix(rep(0, (x - matrix_nrow) * row_num), ncol = row_num))
            matrix_nrow <- x
        }
        knn_distance[1:x, i] <- neighdist
    }
    return(knn_distance[1:matrix_nrow,])
}

#' k_distance_neighbor
#'
#' \code{k_distance_neighbor}  is used for returning knn distance neighbors.
#' This function is not intended to be used by end user. 
#'
#' @param  dat  A data.frame contained only predict variables.
#' @param x  The name of variable to process.
#' @param k Number of neighbors for LOF.
#' @export
k_distance_neighbor <- function(x, dat, k) {
    data_temp <- as.matrix(dat)
    row_num <- dim(dat)[1L]
    dimnames(data_temp) <- NULL
    difference <- scale(data_temp, x, FALSE)
    dist_temp <- drop(difference ^ 2 %*% rep(1, ncol(dat)))
    dist_temp <- sqrt(dist_temp)
    distance_order <- order(dist_temp)
    nn_distance <- dist_temp[distance_order]
    knn_distance <- nn_distance[k + 1]
    neibs <- drop(nn_distance[nn_distance <= knn_distance])
    neibs <- neibs[-1]
    neibr_num <- length(neibs)
    neighbor_index <- distance_order[1:neibr_num + 1]
    neihbor_num1 <- neibr_num + 3
    neihbor_num2 <- neibr_num + neibr_num + 2
    return(c(neihbor_num1, neihbor_num2, neighbor_index, neibs))
}



#' local_reachability_density
#'
#' \code{local_reachability_density} is used for calculating local reachability density.
#'
#' This function is not intended to be used by end user. 
#'
#' @param distance_data A matrix with distance of each observation.
#' @param k Number of neighbors for LOF.
#' @export
local_reachability_density <- function(distance_data, k) {
    p <- dim(distance_data)[2]
    local_reach_density <- rep(0, p)
    local_reach_density <- sapply(1:p, function(i) {
        j <- seq(3, 3 + (distance_data[2, i] - distance_data[1, i]))
        neibr_num <- distance_data[2, i] - distance_data[1, i] + 1
        data_temp <- rbind(diag(distance_data[distance_data[2, distance_data[j, i]], distance_data[j, i]]), distance_data[j + neibr_num, i])
        #calculate reachability
        reach <- 1 / (sum(apply(data_temp, 2, function(x) max(x))) / neibr_num)
        local_reach_density[i] <- reach
    })
    local_reach_density
}
