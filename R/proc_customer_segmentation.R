#' Customer Segmentation
#'
#' \code{customer_segmentation} is  a function for clustering and find the best segment variable.
#' @param  dat  A data.frame contained only predict variables.
#' @param x_list A list of x variables.
#' @param ex_cols A list of excluded variables. Default is NULL.
#' @param  cluster_control  A list controls cluster. kc is the number of cluster center (default is 2), nstart is the number of random groups (default is 1), max_iter max iteration number(default is 100) .
#' \itemize{
#'   \item \code{meth} Method of clustering. Provides two mehods,"Kmeans" and "FCM(Fuzzy Cluster Means)"(default is "Kmeans").
#'   \item \code{kc}  Number of cluster center (default is 2).
#'   \item \code{nstart} Number of random groups (default is 1).
#'   \item \code{max_iter}  Max iteration number(default is 100).
#' }
#' @param  tree_control  A list of controls for desison tree to find the best segment variable.
#' \itemize{
#'   \item \code{cv_folds}  Number of cross-validations(default is 5).
#'   \item \code{maxdepth} Maximum depth of a tree(default is kc +1).
#'   \item \code{minbucket}  Minimum percent of observations in any terminal <leaf> node (default is nrow(dat) / (kc + 1)).
#' }
#' @param file_name  The name for periodically saved segmentation file. Default is "customer_seg".
#' @param dir_path The path for periodically saved segmentation file. Default is "./segmentation".
#' @return  A "data.frame" object contains cluster results.
#' @references
#' Bezdek, James C. "FCM: The fuzzy c-means clustering algorithm".
#' Computers & Geosciences (0098-3004),\url{https://doi.org/10.1016/0098-3004(84)90020-7}
#' @examples
#' \dontrun{ 
#' library(creditmodel)
#' clust <- customer_segmentation(dat = lendingclub,
#' x_list = NULL, ex_cols = "id$|loan_status", 
#' cluster_control = list(meth = "FCM", kc = 2), 
#' tree_control = list(minbucket = round(nrow(lendingclub) / 10)), 
#' file_name = "cus_seg", dir_path = "./seg")
#' dat_sub1 <- subset(lendingclub, purpose == "debt_consolidation")
#' dat_sub2 <- subset(lendingclub, purpose != "debt_consolidation")
#' }
#' @importFrom rpart rpart rpart.control
#' @export

customer_segmentation <- function(dat, x_list = NULL, ex_cols = NULL,
                                  cluster_control = list(meth = "Kmeans", kc = 2,
    nstart = 1, epsm = 1e-06, sf = 2, max_iter = 100),
    tree_control = list(cv_folds = 5, maxdepth = kc + 1, minbucket = nrow(dat) / (kc + 1)),
    file_name = "customer_seg", dir_path = "./segmentation") {
    opt = options(scipen = 100, stringsAsFactors = FALSE, digits = 6) #
    #----------------------
    x_list = get_x_list(x_list = x_list, dat_train = dat, ex_cols = ex_cols)
    dir_path = ifelse(is.null(dir_path) | !is.character(dir_path) || !grepl(".|/", dir_path), "./segmentation", dir_path)
    file_name = ifelse(is.null(file_name) | !is.character(file_name), "customer_seg", file_name)
    if (!dir.exists(dir_path)) dir.create(dir_path)
    if (any(is.na(dat[x_list]))) {
        dat = process_nas(dat = dat, x_list = x_list, default_miss = TRUE, ex_cols =  ex_cols, parallel = FALSE, method = "median",note = FALSE)
    }
    dat = one_hot_encoding(dat)
    dat = low_variance_filter(dat, lvp = 0.9)
    dat_s <- apply(dat, 2, min_max_norm)
    meth = ifelse(!is.null(cluster_control[["meth"]]), cluster_control[["meth"]], "FCM")
    nstart = ifelse(!is.null(cluster_control[["nstart"]]), cluster_control[["nstart"]], 1)
    epsm = ifelse(!is.null(cluster_control[["epsm"]]), cluster_control[["epsm"]], 1e-06)
    sf = ifelse(!is.null(cluster_control[["sf"]]), cluster_control[["sf"]], 2)
    max_iter = ifelse(!is.null(cluster_control[["max_iter"]]), cluster_control[["max_iter"]], 100)
    kc = ifelse(!is.null(cluster_control[["kc"]]), cluster_control[["kc"]], 2)
    if (meth == "FCM") {
        cluster_res = fuzzy_cluster_means(dat = dat_s, kc = kc, nstart = nstart, max_iter = max_iter, sf = sf, epsm = epsm) #Clustering
    } else {
        cluster_res = kmeans(x = dat_s, centers = kc, nstart = nstart, iter.max = max_iter)
    }

    #The result of cluster analysis
    dt_cluster_res = data.frame(dat, cluster_id = cluster_res$cluster)
    save_dt(dt_cluster_res, file_name = file_name, dir_path = dir_path)
    dt_cluster_res$cluster_id = as.factor(as.character(dt_cluster_res$cluster_id))
    #fund the best Clustering variable by desision tree.
    tree_formula <- as.formula(paste("cluster_id", paste(names(dt_cluster_res[, - length(dt_cluster_res)]), collapse = "+"), sep = ' ~ '))
    cv_folds = ifelse(!is.null(tree_control[["cv_folds"]]), tree_control[["cv_folds"]], 5)
    maxdepth = ifelse(!is.null(tree_control[["maxdepth"]]), tree_control[["maxdepth"]], kc + 1)
    minbucket = max(min(table(dt_cluster_res$cluster_id)) / (kc + 1), tree_control[["minbucket"]])
    trcontrol <- rpart.control(minbucket = minbucket, xval = cv_folds, maxdepth = maxdepth)
    set.seed(46)
    fit <- rpart(data = dt_cluster_res, formula = tree_formula
                 , control = trcontrol
                 , parms = list(split = "information"))
    print(fit, minlength = 1, spaces = 20, digits = getOption("digits"))
    plot(fit, uniform = TRUE, branch = 0.4, margin = ifelse(kc > 3, 0.1, 0.2), main = paste(meth,"Customer Segmentation"), compress = TRUE)
    text(fit, use.n = ifelse(kc > 3, FALSE, TRUE), cex = ifelse(kc > 3, 0.6, 0.8), col = 2)
    dev.print(png, file = paste0(dir_path, "/", file_name, "_tree_result_plot.png"), width = 1024, height = 768)
    summary(fit, digits = getOption("digits"), file = paste0(dir_path, "/", file_name, "_tree_result.txt")) 
    options(opt) # reset
    return(cluster_res)
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
        dist = sapply(1:nr, function(j) rowSums(scale(x, y[j, ], FALSE) ^ 2))
    } else {
        nc = ncol(y)
        dist = sapply(1:nc, function(j) colSums(scale(x, y[ ,j], FALSE) ^ 2))
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
