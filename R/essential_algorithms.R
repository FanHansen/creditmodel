#' local_outlier_factor
#' \code{local_outlier_factor}  is function for calculating the lof factor for a data set using knn
#' This function is not intended to be used by end user.
#'
#' @param dat A data.frame contained only predict variables.
#' @param k Number of neighbors for LOF.Default is 10.
#' @export

local_outlier_factor <- function(dat, k = 10) {
    dat = as.matrix(dat)
    row_num = dim(dat)[1L]
    matrix_nrow = k * 2 + 2
    dist_dat = matrix(0, nrow = matrix_nrow, ncol = row_num)
    for (i in 1:row_num) {
        #get k distance neighbors
        data_temp = as.matrix(dat)
        dimnames(data_temp) = NULL
        difference = scale(data_temp, dat[i,], FALSE)
        dist_temp = drop(difference ^ 2 %*% rep(1, ncol(dat)))
        dist_temp = sqrt(dist_temp)
        distance_order = order(dist_temp)
        nn_distance = dist_temp[distance_order]
        knn_distance = nn_distance[k + 1]
        neibs = drop(nn_distance[nn_distance <= knn_distance])
        neibs = neibs[-1]
        neibr_num = length(neibs)
        neighbor_index = distance_order[1:neibr_num + 1]
        neihbor_num1 = neibr_num + 3
        neihbor_num2 = neibr_num + neibr_num + 2
        neighdist = c(neihbor_num1, neihbor_num2, neighbor_index, neibs)
        neighdist[is.na(neighdist)] = 0
        x <- length(neighdist)
        if (x > matrix_nrow) {
            dist_dat = rbind(dist_dat, matrix(rep(0, (x - matrix_nrow) * row_num), ncol = row_num))
            matrix_nrow = x
        }
        dist_dat[1:x, i] = neighdist
        dist_dat[1:matrix_nrow, i]
    }

    #calculating local reachability density
    p = dim(dist_dat)[2]
    lrd_dat = rep(0, p)
    lrd_dat = sapply(1:p, function(i) {
        j = seq(3, 3 + (dist_dat[2, i] - dist_dat[1, i]))
        neibr_num = dist_dat[2, i] - dist_dat[1, i] + 1
        dt_temp = rbind(diag(dist_dat[dist_dat[2, dist_dat[j, i]], dist_dat[j, i]]),
                        dist_dat[j + neibr_num, i])
        #calculate reachability
        reach = 1 / (sum(apply(dt_temp, 2, function(x) max(x))) / neibr_num)
        reach
    })
    lof = rep(0, p)
    # computer the local outlier factor of each observation in dat
    lof = sapply(1:p, function(i) {
        n_neighbor = dist_dat[2, i] - dist_dat[1, i] + 1
        j = seq(0, (n_neighbor - 1))
        local_factor = sum(lrd_dat[dist_dat[3 + j, i]] / lrd_dat[i]) / n_neighbor
        lof[i] = local_factor
    }
    )
    lof
}

#' Fuzzy Cluster means.
#'
#' This function is used for Fuzzy Clustering.
#'
#' @param dat  A data.frame contained only predict variables.
#' @param kc The number of cluster center (default is 2),
#' @param sf Default is 2.
#' @param nstart The number of random groups (default is 1),
#' @param max_iter Max iteration number(default is 100) .
#' @param epsm Default is 1e-06.
#' @references
#' Bezdek, James C. "FCM: The fuzzy c-means clustering algorithm".
#' Computers & Geosciences (0098-3004),\url{https://doi.org/10.1016/0098-3004(84)90020-7}
#' @export
fuzzy_cluster_means <- function(dat, kc = 2, sf = 2, nstart = 1, max_iter = 100, epsm = 1e-06) {
    dat = as.matrix(dat)
    set.seed(46)
    init_centers_id = sample(1:nrow(dat), kc)
    init_centers = as.matrix(dat[init_centers_id,])
    re_best = fuzzy_cluster(dat = dat, kc = kc, init_centers = init_centers, sf = sf, max_iter = max_iter, epsm = epsm)
    if (nstart > 1) {
        i = 2
        while (i <= nstart) {
            centers_id = sample(1:nrow(dat), kc)
            centers = as.matrix(dat[centers_id,])
            rand_best = fuzzy_cluster(dat, kc = kc, init_centers = centers, sf = sf, max_iter = max_iter, epsm = epsm)
            if (rand_best$obj_fun <= re_best$obj_fun) {
                re_best = rand_best
            }
            i = i + 1
        }
    }
    return(re_best)
}


#' @param init_centers Initial centers of obs.
#' @rdname fuzzy_cluster_means
#' @export

fuzzy_cluster <- function(dat, kc = 2, init_centers, sf = 3, max_iter = 100, epsm = 1e-06) {

    dist = sapply(1:kc, function(j) rowSums(scale(dat, init_centers[j,], FALSE) ^ 2))
    history_obj = c()
    obj_fun_old = Inf
    stop_con = TRUE
    iter = 1
    while (stop_con) {
        s = (1 / (dist + epsm)) ^ (1 / (sf - 1))
        md = s / (s %*% matrix(1, kc, kc))
        t1 = t(md ^ sf) %*% dat
        t2 = t(md ^ sf) %*% matrix(1, nrow(dat), ncol(dat))
        centers = t1 / t2
        dist = sapply(1:kc, function(j) rowSums(scale(dat, centers[j,], FALSE) ^ 2))
        obj_fun = sum(md ^ sf * dist)
        stop_con = abs(obj_fun - obj_fun_old) > 0.0001 && (iter < max_iter)
        obj_fun_old = obj_fun
        history_obj = c(history_obj, obj_fun)
        iter = iter + 1
    }
    cluster = apply(md, 1, which.max)
    fcm_res = list(md, obj_fun, history_obj, init_centers, cluster)
    names(fcm_res) = c("membership_degree", "obj_fun", "history_obj", "init_centers", "cluster")
    return(fcm_res)
}


#' euclid_dist
#'
#' This function is not intended to be used by end user.
#'
#' @param  x  A list
#' @param  y  A list
#' @param margin  rows or cols
#' @export
euclid_dist <- function(x, y, margin = 1) {
    x = as.matrix(x)
    y = as.matrix(y)
    if (margin == 1) {
        nr = nrow(y)
        dist = sapply(1:nr, function(j) rowSums(scale(x, y[j,], FALSE) ^ 2))
    } else {
        nc = ncol(y)
        dist = sapply(1:nc, function(j) colSums(scale(x, y[, j], FALSE) ^ 2))
    }
    return(dist)
}

#' cos_sim
#'
#' This function is not intended to be used by end user.
#'
#' @param  x  A list
#' @param  y  A list
#' @param margin  rows or cols
#' @export

cos_sim <- function(x, y, margin = 1) {
    x = as.matrix(x)
    y = as.matrix(y)
    if (margin == 1) {
        nr = nrow(y)
        dist = sapply(1:nr, function(j) colSums(t(x) * y[j,]) / sqrt(rowSums(x ^ 2) * sum(y[j,] ^ 2)))
        colnames(dist) = rownames(y)
    } else {
        nc = ncol(y)
        dist = sapply(1:nc, function(j) colSums(x * y[, j]) / sqrt(colSums(x ^ 2) * sum(y[, j] ^ 2)))
        colnames(dist) = colnames(y)
    }
    return(dist)
}

#' Cramer's V matrix between categorical variables.
#'
#' \code{char_cor_vars} is function for calculating Cramer's V matrix between categorical variables.
#' \code{char_cor} is function for calculating the correlation coefficient between variables by cremers 'V
#' @param dat A data frame.
#' @param x  The name of variable to process.
#' @param x_list Names of independent variables.
#' @param ex_cols A list of excluded variables. Regular expressions can also be used to match variable names. Default is NULL.
#' @param parallel Logical, parallel computing. Default is FALSE.
#' @param note  Logical. Outputs info. Default is TRUE.
#' @return  A list contains correlation index of x with other variables in dat.
#' @examples
#' \dontrun{
#' char_x_list = get_names(dat = UCICreditCard,
#' types = c('factor', 'character'),
#' ex_cols = "ID$|date$|default.payment.next.month$", get_ex = FALSE)
#'  char_cor(dat = UCICreditCard[char_x_list])
#' }
#' @export
char_cor_vars <- function(dat, x) {
    dat = as.data.frame(dat) %>% merge_category() %>% mutate_if(is.factor, as.character)

    vapply(seq_along(dat), function(j) {
        if (length(x) > 1 | length(unlist(x)) > 1) {
            cross_table <- table(unlist(x), dat[, j])
        } else {

            cross_table <- table(dat[, x], dat[, j])

        }
        sqrt(chisq.test(cross_table, correct = TRUE,
        simulate.p.value = TRUE)$statistic / (sum(cross_table, na.rm = TRUE) * min(ncol(cross_table) - 1,
                                                                     nrow(cross_table) - 1, na.rm = TRUE)))

    }, FUN.VALUE = numeric(1))
}

#' @rdname char_cor_vars
#' @export

char_cor <- function(dat, x_list = NULL, ex_cols = "date$", parallel = FALSE, note = FALSE) {
    if (note) {
        cat("[NOTE] Computing the correlation matrix of factor or character variables.\n")
    }
    if (is.null(x_list)) {
        #Obtaining factor or character variables
        x_list = get_names(dat = dat, types = c('factor', 'character'),
                           ex_cols = ex_cols, get_ex = FALSE)
    }
    if (length(x_list) > 0) {
        #calculating the correlation coefficient between variables-cremers'V
        character_cor = loop_function(func = char_cor_vars, x_list = x_list,
                                      args = list(dat = dat[x_list]),
                                      bind = "cbind", as_list = FALSE,
                                      parallel = parallel)
        character_cor = as.matrix(character_cor)
        colnames(character_cor) = rownames(character_cor) = x_list
    } else {
        character_cor = NULL
    }
    return(character_cor)
}

#' auc_value
#' \code{auc_value} is for get best lambda required in lasso_filter. This function required in \code{lasso_filter}
#' @param prob A list of redict probability or score.
#' @param target Vector of target.
#' @return Lanmbda value
#' @export


auc_value = function(target, prob) {
    prob.rank = rank(prob)
    cnt_1 = sum(target)
    cnt_0 = length(target) - cnt_1
    prob_1 = prob.rank[target == 1]
    u = sum(prob_1) - cnt_1 * (cnt_1 + 1) / 2
    exp(log(u) - log(cnt_1) - log(cnt_0))
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



