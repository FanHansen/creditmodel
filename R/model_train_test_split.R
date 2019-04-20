#' Train-Test-Split
#'
#' \code{train_test_split} Functions for partition of data.
#' @param dat A data.frame with independent variables and target variable.
#' @param prop The percentage of train data samples after the partition.
#' @param split_type  Methods for partition. 
#' \itemize{
#'   \item "Random" is to split train & test set randomly. 
#'   \item "OOT" is to split by time for observation over time test.
#'   \item "byRow" is to split by rownumbers.
#' }
#' @param occur_time The name of the variable that represents the time at which each observation takes place. It is used for "OOT" split.
#' @param cut_date Time points for spliting data sets, e.g. : spliting Actual and Expected data sets.
#' @param start_date The earliest occurrence time of observations.
#' @param save_data Logical, save results in locally specified folder. Default is FALSE.
#' @param file_name The name for periodically saved data file. Default is "dat".
#' @param dir_path The path for periodically saved data file. Default is "./data".
#' @param seed  Random number seed. Default is 46.
#' @param note Logical. Outputs info. Default is TRUE.
#' @return A list of indices (train-test)
#' @examples
#' \dontrun{
#' train_test <- train_test_split(lendingclub, 
#' split_type = "OOT", prop = 0.7, 
#' occur_time = "issue_d", seed = 12, save_data = FALSE)
#' dat_train = train_test$train
#' dat_test = train_test$test
#' }
#' @importFrom stats quantile ecdf
#' @export




train_test_split <- function(dat, prop = 0.7, split_type = c("Random", "OOT", "byRow"),
occur_time = NULL, cut_date = NULL, start_date = NULL,  save_data = FALSE,
dir_path = "./data/", file_name = "dat", note = TRUE, seed = 43) {
    if (prop > 1 || !is.numeric(prop)) {
        warning("[Invalid]  prop is not a numeric or more than 1,  reset to 0.7.\n")
        prop = 0.7
    }
    if (!is.element(split_type, c("OOT", "Random", "byRow"))) {
        stop("split_type must be either 'OOT' or 'Random' or 'byRow'.\n")
    }
    if (length(split_type)> 1) {
        warning("your split_type is more than one and only the first one is selected.\n")
    }
    if (length(split_type) == 0) {
        warning("split_type is missing,  set 'Random' by default.\n")
        split_type = "Random"
    }
    if (split_type[1] == "OOT" & !is.null(occur_time) && any(names(dat) == occur_time)) {
        dat =  time_transfer(dat, date_cols = occur_time)
        if (is_date(dat[, occur_time])) {
            if (is.null(cut_date)) {
                cut_date = date_cut(dat_time = dat[, occur_time], pct = prop)
            }
            if (is.null(start_date)) {
                start_date = date_cut(dat_time = dat[, occur_time], pct = 0)
            }
            dat[, occur_time] = as.Date(dat[, occur_time])
            test = dat[which(dat[,occur_time] >= cut_date),]
            train = dat[which(dat[, occur_time] >= start_date & dat[, occur_time] < cut_date),]
            if (note) cat(paste("[NOTE]", "total:", nrow(dat), "--->test:",
            nrow(test), " train:", nrow(train), ".\n", sep = "", collapse = "\n"))
            } else {
                if (!is.null(seed)) set.seed(seed) else set.seed(46)
                sub = sample(1:nrow(dat), round(nrow(dat) * prop))
                train = dat[sub,]
                test = dat[-sub,]
                if (note) cat(paste("[NOTE]", "total:", nrow(dat), "--->test:",
                nrow(test), " train:", nrow(train), ".\n", sep = "", collapse = "\n"))
                warning(paste(occur_time, "is  not date or time, unable to use OOT , split random.\n"))
            }

        } else {
            if (!is.null(seed)) set.seed(seed) else set.seed(46)
            sub = sample(1:nrow(dat), round(nrow(dat) * prop))
            train = dat[sub,]
            test = dat[-sub,]
            if (note) cat(paste("[NOTE]", "total:", nrow(dat), "--->test:",
            nrow(test), " train:", nrow(train), ".\n", sep = "", collapse = "\n"))
    }
    if (split_type[1] == "Random") {
        if (!is.null(seed)) set.seed(seed) else set.seed(46)
        sub = sample(1:nrow(dat), round(nrow(dat) * prop))
        train = dat[sub,]
        test = dat[-sub,]
        if (note) cat(paste("[NOTE]", "total:", nrow(dat), "--->test:",
        nrow(test), " train:", nrow(train), ".\n", sep = "", collapse = "\n"))
    }
    if (split_type[1] == "byRow") {   
        sub = 1:round(nrow(dat) * prop)
        train = dat[sub,]
        test = dat[-sub,]
        if (note) (paste("[NOTE]", "total:", nrow(dat), "--->test:",
        nrow(test), " train:", nrow(train), ".\n", sep = "", collapse = "\n"))
    }
    if (save_data) {
        if (!dir.exists(dir_path)) dir.create(dir_path)
        save_dt(train, file_name = paste("train", file_name, sep = "_"), dir_path = dir_path)
        save_dt(test, file_name = paste("test", file_name, sep = "_"), dir_path = dir_path)
    }
    return(list(test = test, train = train))
}


#' Stratified Folds
#'
#' this function creates stratified folds for cross validation.
#'
#' @param dat  A data.frame.
#' @param k  k is an integer specifying the number of folds.
#' @param occur_time time variable for creating OOT folds. Default is NULL.
#' @param seed A seed. Default is 46.
#' @return a list of indices
#' @examples
#' \dontrun{
#' CV_folds = cv_split(dat = lendingclub, k=5)
#' }
#' @importFrom stats quantile ecdf
#' @export
cv_split <- function(dat,k= 5, occur_time = NULL,seed = 46) {
    cv_list = list()
    dat = checking_data(dat = dat, occur_time = occur_time)
    if (!is.null(seed)) set.seed(seed) else set.seed(46)
    if (!is.null(occur_time) && is.element(occur_time, names(dat)) && is_date(dat[, occur_time])) {
        date_n = quantile(ecdf(dat[, occur_time]), seq(0, 1, by = 0.01))
        date_q = as.double(sub("%", "", names(date_n))) / 100
        prop = round(1 / k, 2)
        date_temp = date_n[which(date_q == prop)]
        if (nchar(date_temp)<= 7) {
            for (i in 1:k) {
                cv_list[[i]] = which(dat[, occur_time] >= min(as.Date(date_n[which(date_q >= prop * (k - i))], origin = "1970-01-01")) &
                dat[, occur_time] < min(as.Date(date_n[which(date_q >= prop * (k - i + 1))], origin = "1970-01-01")))
            }
        } else {
            for (i in 1:k) {
                cv_list[[i]] = which(dat[, occur_time] >= min(as.Date.POSIXct(date_n[which(date_q >= prop * (k - i))], origin = "1970-01-01")) &
                dat[, occur_time] < min(as.Date.POSIXct(date_n[which(date_q >= prop * (k - i + 1))], origin = "1970-01-01")))
            }
        }
        null_cv = which(sapply(cv_list, function(i) length(i) == 0))
        if (length(null_cv) > 0) {
            not_in_cv = which(!1:nrow(dat) %in% unlist(cv_list))
            for (i in null_cv[which(null_cv < k)]) {
                cv_list[[i]] = sample(not_in_cv, ceiling(length(not_in_cv) / length(null_cv)))
            }
            cv_list[[k]] = which(!1:nrow(dat) %in% unlist(cv_list))
        }

        cv_list = cv_list[which(sapply(cv_list, function(i) length(i) > 0))]
        length(unlist(cv_list))
    } else {
        nr = nrow(dat)
        chaos_n = sample(rep(1:k, ceiling(nr / k))[1:nr], nr) 
        dat_seq = 1:nr
        cv_list = lapply(1:k, function(x) dat_seq[chaos_n == x]) 
    }  
    return(cv_list)
}