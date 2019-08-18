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
#' train_test <- train_test_split(lendingclub,
#' split_type = "OOT", prop = 0.7,
#' occur_time = "issue_d", seed = 12, save_data = FALSE)
#' dat_train = train_test$train
#' dat_test = train_test$test
#' @importFrom stats quantile ecdf
#' @export

train_test_split <- function(dat, prop = 0.7, split_type = c("Random", "OOT", "byRow"),
occur_time = NULL, cut_date = NULL, start_date = NULL, save_data = FALSE,
dir_path = tempdir(), file_name = NULL, note = FALSE, seed = 43) {

    if (prop > 1 || !is.numeric(prop)) {
        warning("[Invalid]  prop is not a numeric or more than 1,  reset to 0.7.\n")
        prop = 0.7
    }
    if (!is.element(split_type, c("OOT", "Random", "byRow"))) {
        stop("split_type must be either 'OOT' or 'Random' or 'byRow'.\n")
    }
    if (length(split_type) > 1) {
        warning("your split_type is more than one and only the first one is selected.\n")
    }
    if (length(split_type) == 0) {
        warning("split_type is missing,  set 'Random' by default.\n")
        split_type = "Random"
    }
    if (split_type[1] == "OOT" & !is.null(occur_time) && any(names(dat) == occur_time)) {
        dat = time_transfer(dat, date_cols = occur_time)
        if (is_date(dat[, occur_time])) {
            if (is.null(cut_date)) {
                cut_date = date_cut(dat_time = dat[, occur_time], pct = prop)
            }
            if (is.null(start_date)) {
                start_date = date_cut(dat_time = dat[, occur_time], pct = 0)
            }
            dat[, occur_time] = as.Date(dat[, occur_time])
            test = dat[which(dat[, occur_time] >= cut_date),]
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
        dir_path = ifelse(!is.character(dir_path),
                      tempdir(), dir_path)
        if (!dir.exists(dir_path)) dir.create(dir_path)
        if (!is.character(file_name)) file_name = NULL
        save_dt(train, file_name = ifelse(is.null(file_name), "dat.train", paste(file_name, "dat.train", sep = ".")), dir_path = dir_path)
        save_dt(test, file_name = ifelse(is.null(file_name), "dat.test", paste(file_name, "dat.test", sep = ".")), dir_path = dir_path)
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
#' sub = cv_split(UCICreditCard, k = 30)[[1]]
#' dat = UCICreditCard[sub,]
#' @importFrom stats quantile ecdf
#' @export
cv_split <- function(dat, k = 5, occur_time = NULL, seed = 46) {
    cv_list = list()
    dat = checking_data(dat = dat, occur_time = occur_time)
    if (!is.null(seed)) set.seed(seed) else set.seed(46)
    if (!is.null(occur_time) && is.element(occur_time, names(dat)) &&
            is_date(dat[, occur_time])) {
        date_n = quantile(ecdf(dat[, occur_time]), seq(0, 1, by = 0.01))
        date_q = as.double(sub("%", "", names(date_n))) / 100
        prop = round(1 / k,  2)
        date_temp = date_n[which(date_q == prop)]
        if (length(date_temp) >0 && nchar(date_temp) <= 7) {
            for (i in 1:k) {
                cv_list[[i]] = which(dat[, occur_time] >= min(as.Date(date_n[which(date_q >= prop * (k - i))],
                                                                      origin = "1970-01-01")) &
                                       dat[, occur_time] < min(as.Date(date_n[which(date_q >= prop * (k - i + 1))],
                                                                       origin = "1970-01-01")))
            }
        } else {
            for (i in 1:k) {
                cv_list[[i]] = which(dat[, occur_time] >= min(as.Date.POSIXct(date_n[which(date_q >= prop * (k - i))],
                                                                              origin = "1970-01-01")) &
                                       dat[, occur_time] < min(as.Date.POSIXct(date_n[which(date_q >= prop * (k - i + 1))],
                                                                               origin = "1970-01-01")))
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

    } else {
        nr = nrow(dat)
        k = ceiling(k)
        if (!is.null(seed)) set.seed(seed) else set.seed(46)
        chaos_n = sample(rep(1:k, ceiling(nr / k))[1:nr], nr)
        dat_seq = 1:nr
        cv_list = lapply(1:k, function(x) dat_seq[chaos_n == x])
    }
    return(cv_list)
}



#' Packages required and intallment
#'
#' \code{require_packages} is function for librarying required packages and  installing missing packages if needed.
#' @param pkg A list or vector of names of required packages.
#' @return  packages installed and library.
#' @examples
#' \dontrun{
#' require_packages(c("dplyr","ggplot2"))
#' }
#' @export
require_packages <- function(pkg) {
    opt = options("warn" = -1) # suppress warnings
    new_pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new_pkg) > 0) {
        install.packages(new_pkg, dependencies = TRUE)
        cat("[NOTE] Installs missing packages if needed")
        cat(paste("[NOTE]", "packages", new_pkg, " are installed!"))
    }
    sapply(pkg, require, ch = TRUE)
    options(opt) # reset warnings
}


#' List as data.frame quickly
#'
#' \code{quick_as_df} is function for fast dat frame  transfromation.
#' @param df_list A list of data.
#' @return  packages installed and library,
#'
#' @examples
#'
#' UCICreditCard = quick_as_df(UCICreditCard)
#'
#' @export
quick_as_df <- function(df_list) {
    class(df_list) <- "data.frame"
    attr(df_list, "row.names") <- .set_row_names(length(df_list[[1]]))
    df_list
}


#' Save data
#'
#' \code{save_dt} is for saving a data.frame or a list fast.
#' @param dt A data frame or list.
#' @param file_name A string. The file name to save breaks_list.
#' @param dir_path  A string. The dir path to save breaks_list.
#' @param note  Logical. Outputs info.Default is TRUE.
#' @param as_list  Logical. List format or data.frame format to save. Default is FALSE.
#' @param row_names  Logical,retain rownames.
#' @param append Logical, append newdata to old.
#' @examples
#' save_dt(UCICreditCard,"UCICreditCard", tempdir())
#' @importFrom data.table fwrite fread dcast melt
#' @export

save_dt <- function(dt, file_name = "dat", dir_path = getwd(),
                    note = FALSE, as_list = FALSE, row_names = FALSE, append = FALSE) {

    if (!dir.exists(dir_path)) dir.create(dir_path)
    if (dir.exists(paste0(dir_path, '/', file_name, ".csv")))file.remove(list.files(paste0(dir_path, '/', file_name, ".csv"), recursive = TRUE, full.names = TRUE))
    if (note) {
        cat(paste0("[NOTE]", "Saved to ", "'",paste0(dir_path, '/', file_name, ".csv"),"'\n"))
    }
    if (as_list) {
        fwrite(list(dt), paste0(dir_path, '/', file_name, ".csv"),
               append = append,col.names = FALSE)
    } else {
        fwrite(as.data.frame(dt), paste0(dir_path, '/', file_name, ".csv"),
               append = append, row.names = row_names)
    }
}

#' Number of digits
#'
#' \code{digits_num} is for caculating optimal digits number for numeric variables.
#' @param dat_x  A numeric variable.
#' @return  A number of digits
#' @examples
#' \dontrun{
#' digits_num(lendingclub[,"dti"])
#' # 7
#' }
#' @export


digits_num =function(dat_x) {
    options(scipen = 200, stringsAsFactors = FALSE)
    digits1 = digits2 = 16
    dat_x = unique(unlist(dat_x[!is.na(dat_x)]))
    if (length(dat_x) > 0 && any(is.element(class(dat_x), c("integer", "numeric",
        "double")))) {
        digits1 = vapply(dat_x, function(num) {
            char_num = as.character(gsub("-", "",
                num))
            n_num = as.numeric(char_num) %% 1
            if (!is.null(n_num) && !is.na(n_num) && is.numeric(n_num)) {
                if (n_num == 0) {
                    t_lens = nchar(char_num)
                    left_comma = t_lens
                    right_comma = t_lens - left_comma
                }
                else {
                    comma_p = gregexpr("[.]", char_num)[[1]][1]
                    t_lens = nchar(char_num)
                    left_comma = comma_p - 1
                    right_comma = t_lens - 1 - left_comma
                }
                right_comma
            }
        }, FUN.VALUE = numeric(1))
        digits2 = max(digits1)
    }
	digits2 = ifelse(digits2 > 16, 16, digits2)
	return(digits2)
}




#' is_date
#'
#' \code{is_date} is a small function for distinguishing time formats
#' @param x  list or vectors
#' @return  A Date.
#' @examples
#' is_date(lendingclub$issue_d)
#' @export
is_date = function(x) {
    any(class(x) %in% c("Date", "POSIXlt", "POSIXct", "POSIXt"))
}

#' Date Time Cut Point
#'
#' \code{date_cut} is  a small function to get date point.
#' @param dat_time  time vectors.
#' @param pct  the percent of cutting. Default: 0.7.
#' @return  A Date.
#' @examples
#' date_cut(dat_time = lendingclub$issue_d, pct = 0.8)
#' #"2018-08-01"
#' @importFrom stats quantile ecdf
#' @export

date_cut = function(dat_time, pct = 0.7) {
    dat_time = as.Date(dat_time)
    if (is_date(dat_time)) {
        date_n = quantile(ecdf(dat_time), seq(0, 1, by = 0.01))
        date_q = as.double(sub("%", "", names(date_n))) / 100
        date_temp = date_n[which(date_q == pct)]
        if (nchar(date_temp) > 7) {
            cut_date= min(as.Date.POSIXct(date_n[which(date_q >= pct)], origin = "1970-01-01"))
        } else {
            cut_date = min(as.Date(date_n[which(date_q >= pct)], origin = "1970-01-01"))
        }
        return(cut_date)
    } else {
        stop(paste("Not Date or Time.\n"))
    }
}


#' Percent Format
#'
#' \code{as_percent} is  a small function for making percent format..
#' @param x  A numeric vector or  list.
#' @param digits  Number of digits.Default: 2.
#' @return  x with percent format.
#' @examples
#' as_percent(0.2363, digits = 2)
#' as_percent(1)
#' @export
as_percent <- function(x, digits = 2) {
    x = as.numeric(x)
    pct = round(x, digits) * 100
    x_pct = paste0(pct, ifelse(is.finite(x), "%", ""))
    return(x_pct)
}

#' %islike%
#' Fuzzy String matching
#'
#' @param x  A string.
#' @param y  A string.
#' @return  Logical.
#' @examples
#'  "xyz"  %islike% "yz$"
#' @export

'%islike%' <- function(x, y) {
    grx = FALSE
    x = gsub("\\$|\\*|\\+|\\?|\\[|\\^|\\{|\\}|\\\\|\\(|\\)|\\|\\)|\\]", "", x)
    y = gsub("\\{|\\}", "", y)

    if (any(x != '') & any(y != '') & any(!is.null(x)) & any(!is.null(y)) & any(!is.na(x)) & any(!is.na(y))) {
        y = unique(unlist(y))
        y = y[which(y != '')]
        if (length(y[!grepl("\\$|\\|", y)]) > 0) {
            grx = Reduce(function(u, v) paste(u, v, sep = '|'), paste0("^", y[!grepl("\\|", y)], sep = "$"))
        }
        if (length(y[grepl("\\$|\\|", y)]) > 0) {
            grx = paste(grx, y[grepl("\\|", y)], sep = '|')
        }
    }
    grepl(grx, x)
}


#' %alike%
#' Fuzzy String matching
#' @param x  A string.
#' @param y  A string.
#' @return  Logical.
#' @examples
#' "xyz"  %alike% "xy"
#' @export

'%alike%' <- function(x, y) {
    x = gsub("\\$|\\*|\\+|\\?|\\[|\\^|\\{|\\}|\\\\|\\(|\\)|\\|\\)|\\]", "", x)
    if (any(x != '') & any(y != '') & any(!is.null(x)) & any(!is.null(y)) & any(!is.na(x)) & any(!is.na(y))) {
        y = unlist(y)
        y = y[which(y != '')]
        x = unlist(x)
        x[which(x != '')]
        grx1 = Reduce(function(u, v) paste(u, v, sep = '|'), y)
        grx2 = Reduce(function(u, v) paste(u, v, sep = '|'), x)
        if (any(grepl(grx1, x))) {
            grpl2 = grepl(grx1, x)
        } else {
            grpl2 = grepl(grx2, y)
        }
    } else {
        grpl2 = grepl(FALSE, x)
    }
    return(grpl2)
}


#'  Rename
#'
#' \code{re_name} is  for renaming variables.
#' @param dat A data frame with vairables to rename.
#' @param newname  New names of vairables.
#' @param oldname  Old names of vairables.
#' @return  data with new variable names.
#' @examples
#' dt = re_name(dat = UCICreditCard, "default.payment.next.month" , "target")
#' names(dt['target'])
#' @export


re_name = function(dat, oldname = c(), newname = c()) {
    ind <- which(names(dat) %in% oldname)
    names(dat)[ind] = newname
    dat
}

#' Min Max Normalization
#'
#' \code{min_max_norm} is for normalizing each column vector of matrix 'x' using min_max normalization
#' @param x  Vector
#' @return Normalized vector
#' @examples
#' dat_s = apply(UCICreditCard[,12:14], 2, min_max_norm)
#' @export

min_max_norm <- function(x) {
    ((x - min(x , na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)))
}
#' Functions for vector operation.
#'
#' @param x  A data.frame or Matrix.
#' @param na.rm  Logical, remove NAs.
#' @return  A data.frame or Matrix.
#' @examples
#' #any row has missing values
#' row_amy =  rowAny(UCICreditCard[8:10])
#' #rows which is all missing values
#' row_na =  rowAllnas(UCICreditCard[8:10])
#' #cols which is all missing values
#' col_na =  colAllnas(UCICreditCard[8:10])
#' #cols which is all zeros
#' row_zero =  colAllzeros(UCICreditCard[8:10])
#' #sum all numbers of a row
#' row_all =  rowAll(UCICreditCard[8:10])
#' #caculate cv of a row
#' row_cv =  rowCVs(UCICreditCard[8:10])
#' #caculate sd of a row
#' row_sd =  rowSds(UCICreditCard[8:10])
#' #caculate sd of a column
#' col_sd =  colSds(UCICreditCard[8:10])
#' @export

rowAny <- function(x) rowSums(x, na.rm = TRUE) > 0

#' @rdname rowAny
#' @export
rowAllnas <- function(x) rowSums(is.na(x)) == length(x)

#' @rdname rowAny
#' @export
colAllnas <- function(x) colSums(is.na(x)) == nrow(x)


#' @rdname rowAny
#' @export
colAllzeros <- function(x) colSums(x) == 0

#' @rdname rowAny
#' @export
rowAll <- function(x) rowSums(x, na.rm = TRUE) == ncol(x)


#' @rdname rowAny
#' @export

rowCVs <- function(x, na.rm = FALSE) {
    ifelse(rowAllnas(x), NA,
            ifelse(rowAny(x) > 0,
            sqrt(rowSums((x - rowMeans(x, na.rm = na.rm)) ^ 2, na.rm = na.rm) / length(x)) / rowMeans(x, na.rm = na.rm), 0))
}

#' @rdname rowAny
#' @export
rowSds <- function(x, na.rm = FALSE) sqrt(rowSums((x - rowMeans(x, na.rm = na.rm)) ^ 2, na.rm = na.rm) / length(x))

colSds <- function(x, na.rm = FALSE) sqrt(colSums((x - colMeans(x, na.rm = na.rm)) ^ 2, na.rm = na.rm) / nrow(x))
#' @rdname rowAny
#' @export

colSds <- function(x, na.rm = TRUE) {
    lapply(x, function(i)sd(i,na.rm = TRUE))
}

#' @rdname rowAny
#' @export
rowMaxs <- function(x, na.rm = FALSE) {
    maxs <- apply(x, 1, function(i) ifelse(length(i) > 1, i[which.max(i)], i))
    as.numeric(maxs)
}

#' @rdname rowAny
#' @export
rowMins <- function(x, na.rm = FALSE) {
    mins <- apply(x, 1, function(i) i[which.min(i)])
    as.numeric(mins)
}

#' @rdname rowAny
#' @export
rowMaxMins <- function(x, na.rm = FALSE) {
    max_mins <- apply(x, 1, function(i) i[which.max(i)] - i[which.min(i)])
    as.numeric(max_mins)
}

#' @rdname rowAny
#' @export
colMaxMins <- function(x, na.rm = FALSE) {
    max_mins <- apply(x, 2, function(i) i[which.max(i)] - i[which.min(i)])
    as.numeric(max_mins)
}


#'  Get Variable Names
#'
#' \code{get_names} is  for getting names of particular classes of variables
#' @param dat A data.frame with independent variables and target variable.
#' @param types  The class or types of variables which names to get. Default: c('numeric', 'integer', 'double')
#' @param ex_cols A list of excluded variables. Regular expressions can also be used to match variable names. Default is NULL.
#' @param get_ex  Logical ,if TRUE, return a list contains names of excluded variables.
#' @return  A list contains names of variables
#' @seealso \code{\link{get_x_list}}
#' @examples
#' x_list = get_names(dat = UCICreditCard, types = c('factor', 'character'),
#' ex_cols = c("default.payment.next.month","ID$|_date$"), get_ex = FALSE)
#' x_list = get_names(dat = UCICreditCard, types = c('numeric', 'character', "integer"),
 #' ex_cols = c("default.payment.next.month", "ID$|SEX "), get_ex = FALSE)
#' @export



get_names <- function(dat, types = c('logical', 'factor', 'character', 'numeric',
                                        'integer', 'double', "Date", "POSIXlt", "POSIXct", "POSIXt"), ex_cols = NULL, get_ex = FALSE) {
    if (is.null(types)) {
        stop("types is missing!")
    }
    if (is.null(ex_cols)) {
        sel_names = names(dat)[sapply(dat, function(x) any(is.element(class(x), types)))]
        ex_names = names(dat)[!(sapply(dat, function(x) any(is.element(class(x), types))))]
    } else {
        var_names = names(dat)[sapply(dat, function(x) any(is.element(class(x), types)))]
        if (length(ex_cols) > 1 || !grepl("\\$|\\*|\\+|\\?|\\[|\\^|\\{|\\}|\\\\|\\|\\)|\\]", ex_cols) ) {
            ex_vars = names(dat)[colnames(dat) %in% ex_cols]
        } else {
            ex_vars = names(dat)[colnames(dat) %alike% ex_cols]
        }
        ex_types = names(dat)[!(sapply(dat, function(x) any(is.element(class(x), types))))]
        ex_names = unique(c(ex_vars, ex_types))
        sel_names = setdiff(var_names, ex_names)
    }
    if (get_ex) {
        dat = dat[ex_names]
    } else {
        dat = dat[sel_names]
    }
    var_names <- names(dat)
    return(var_names)
}


#' Get X List.
#'
#' \code{get_x_list} is  for getting intersect names of x_list, train and test.
#' @param dat_train  A data.frame with independent variables.
#' @param dat_test  Another data.frame.
#' @param x_list Names of independent variables.
#' @param ex_cols A list of excluded variables. Regular expressions can also be used to match variable names. Default is NULL.
#' @return  A list contains names of variables
#' @seealso \code{\link{get_names}}
#' @examples
#' x_list = get_x_list(x_list = NULL,dat_train = UCICreditCard,
#' ex_cols = c("default.payment.next.month","ID$|_date$"))
#' @export


get_x_list <- function(x_list = NULL, dat_train = NULL, dat_test = NULL, ex_cols = NULL) {
    if (!is.null(dat_train)) {
        if (is.null(x_list) | length(x_list) <1 ) {
            if (is.null(dat_test)) {
                x_list_retain = get_names(dat = dat_train,
                                          types = c('character', 'factor', 'numeric', 'integer', 'double'),
                                          ex_cols = ex_cols, get_ex = FALSE)
            } else {
                x_list_t = get_names(dat = dat_train,
                                     types = c('character', 'factor', 'numeric', 'integer', 'double'),
                                     ex_cols = ex_cols, get_ex = FALSE)
                x_list_s = get_names(dat = dat_test,
                                     types = c('character', 'factor', 'numeric', 'integer', 'double'),
                                     ex_cols = ex_cols, get_ex = FALSE)
                x_list_retain = intersect(x_list_t, x_list_s)
            }
        } else {
            if (is.null(dat_test)) {
                x_list_ts = get_names(dat = dat_train,
                                      types = c('character', 'factor', 'numeric', 'integer', 'double'),
                                      ex_cols = ex_cols, get_ex = FALSE)
                x_input = x_list %in% x_list_ts
            } else {
                x_list_t = get_names(dat = dat_train,
                                     types = c('character', 'factor', 'numeric', 'integer', 'double'),
                                     ex_cols = ex_cols, get_ex = FALSE)
                x_list_s = get_names(dat = dat_test,
                                     types = c('character', 'factor', 'numeric', 'integer', 'double'),
                                     ex_cols = ex_cols, get_ex = FALSE)
                x_list_ts = intersect(x_list_t, x_list_s)
                x_excluded_ts = setdiff(x_list_t, x_list_s)
                if (length(x_excluded_ts) > 0) {
                    cat(paste("[EXCLUDED]", "variables which are not both in train & test :\n",
                              paste(x_excluded_ts, collapse = " \n"), ".\n",sep = ""))
                }
                x_input = x_list %in% x_list_ts
            }
            if (!all(x_input)) {
                x_retain = x_list[which(x_input)]
                x_excluded = x_list[which(!x_input)]
                cat(paste("[EXCLUDED]", "which  are not in x_list: \n",
                          paste(x_excluded, collapse = " \n"), ".\n",sep = ""))
            } else {
                x_retain = x_list
            }
            x_type = sapply(dat_train[, x_retain],
                            class) %in% c('character', 'factor', 'numeric', 'integer', 'double')

            if (!all(x_type)) {
                x_list_retain = x_retain[which(x_type)]
                x_excluded = x_retain[which(!x_type)]
                cat(paste("[EXCLUDED]", "which  are not numeric or character or factor:\n",
                          paste(x_excluded, collapse = " \n"), ".\n", sep = ""))
            } else {
                x_list_retain = x_retain
            }
        }
    }
    return(x_list_retain)
}


#'Checking Data
#'
#'
#' \code{checking_data}  cheking dat before processing.
#' @param dat A data.frame with independent variables and target variable.
#' @param target The name of target variable. Default is NULL.
#' @param occur_time The name of the variable that represents the time at which each observation takes place.
#' @param pos_flag The value of positive class of target variable, default: "1".
#' @param note Logical.Outputs info.Default is TRUE.
#' @return data.frame
#' @examples
#' dat = checking_data(dat = UCICreditCard, target = "default.payment.next.month")
#' @export




checking_data <- function(dat = NULL, target = NULL, occur_time = NULL,
                          note = FALSE, pos_flag = NULL) {
    if (note) {
        (cat("[NOTE] checking dat and target format.\n"))
    }
    if (is.null(dat)) {
        warning("dat is null.\n")
    } else {
        if (!(class(dat)[1] == "data.frame")) {
            if (any(is.element(class(dat), c("data.table", "list", "tbl_df", "tbl", "matrix"))) && length(dim(dat)) == 2) {
                dat = as.data.frame(dat)
                cat(paste("[NOTE]", "convert", class(dat)[1], "to data.frame.\n"))
            } else {
                warning("Dat is not two-dimensional.\n")
            }
        }
    }
    if (!is.null(target)) {
        if (!is.character(target) || length(target) > 1) {
            warning(paste("target is not a string or a name.\n", sep = "\t"))
        } else {
            if (length(unique(dat[, target])) < 2) {
                warning(paste("Unique values of target is only one.\n", sep = "\t"))
            } else {
                if (length(unique(dat[, target])) == 2) {

                    if (is.null(pos_flag)) {
                        pos_flag = list("1", "bad", 1, "Bad", "positive", "pos", "Positive", "Pos")
                    }
                    if (length(which(dat[, target] %in% pos_flag)) != 0) {
                        if (!is.numeric(dat[, target]) || is.numeric(dat[, target]) && !all(sort(unique(dat[, target])) == c(0, 1))) {
                            pos <- unique(dat[, target])[which(unique(dat[, target]) %in% pos_flag)]
                            dat[, target] = ifelse(dat[, target] %in% pos_flag, 1, 0)
                            warning(paste("The  values in of target has been encoded", pos, "=1 and others = 0", " \n"))
                        }
                    } else {
                        warning(paste("Positive values of", target, "is not in pos_flag:", paste(pos_flag, collapse = ","), "\nplease set pos_flag. \n", sep = "\t"))
                    }
                }
            }
        }
    }
    if (!is.null(occur_time)) {
        if (is.element(occur_time, names(dat))) {
            dat <- time_transfer(dat, date_cols = occur_time)
            if (!is_date(dat[, occur_time])) {
                warning(paste("occur_time:", occur_time, "is not time or date.\n", sep = "\t"))
            }
        } else {
            warning(paste("occur_time:", occur_time, "is not in data.\n", sep = "\t"))
        }
    }
    return(dat)
}

#' Parallel computing and export variables to global Env.
#'
#' This function  is not intended to be used by end user.
#' @param parallel  A logical, default is TRUE.
#'
#' @return  parallel works.
#' @importFrom parallel detectCores  clusterExport clusterCall makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar% %do%  registerDoSEQ
#' @export
start_parallel_computing <- function(parallel = TRUE) {

    parallelType <- if (.Platform$OS.type == "windows")
        "snow" else "multicore"

    numCores <- parallel::detectCores()
    attr(parallel, "type") <- parallelType
    attr(parallel, "cores") <- numCores

    if (parallel) {
        if (parallelType == "snow") {

            cl <- parallel::makeCluster(numCores, type = "PSOCK")
            attr(parallel, "cluster") <- cl

            varlist <- ls(envir = parent.frame(), all.names = TRUE)
            varlist <- varlist[varlist != "..."]
            parallel::clusterExport(cl, varlist = varlist, envir = parent.frame())
            parallel::clusterExport(cl, varlist = ls(envir = globalenv(), all.names = TRUE), envir = globalenv())

            pkgs <- .packages()
            lapply(pkgs, function(pkg)
                parallel::clusterCall(cl, library, package = pkg,
                                     character.only = TRUE))
            doParallel::registerDoParallel(cl, cores = numCores)
        }
        else if (parallelType == "multicore") {
            cl <- parallel::makeCluster(numCores, type = "FORK")

            doParallel::registerDoParallel(cl, cores = numCores[1])

            attr(parallel, "cluster") <- cl
        }
        else { stop("Only 'snow' and 'multicore' clusters allowed!") }
        }
    return(parallel)
}

#' Stop parallel computing
#'
#' This function  is not intended to be used by end user.
#' @param cluster  Parallel works.
#'
#' @return  stop clusters.
#'
#' @importFrom parallel  stopCluster
#' @importFrom foreach  registerDoSEQ
#' @export
stop_parallel_computing <- function(cluster) {
    parallel::stopCluster(cluster)
    foreach::registerDoSEQ()
    invisible()
}



#'  Loop Function.
#' #' \code{loop_function} is an iterator to loop through
#' @param func  A function.
#' @param args  A list of argauments required by function.
#' @param x_list  Names of objects to loop through.
#' @param bind  Complie results, "rbind" & "cbind" are available.
#' @param parallel  Logical, parallel computing.
#' @param as_list  Logical, whether outputs  to be a list.
#'
#' @return   A data.frame or list
#'
#' @examples
#' dat = UCICreditCard[24:26]
#' num_x_list = get_names(dat = dat, types = c('numeric', 'integer', 'double'),
#'                       ex_cols = NULL, get_ex = FALSE)
#' dat[ ,num_x_list] = loop_function(func = outliers_kmeans_lof, x_list = num_x_list,
#'                                    args = list(dat = dat),
#'                                    bind = "cbind", as_list = FALSE,
#'                                  parallel = FALSE)
#' @importFrom foreach foreach %dopar% %do%
#' @export

loop_function <- function(func = NULL, args = list(data = NULL), x_list = NULL,
                          bind = "rbind", parallel = TRUE, as_list = FALSE) {
    opt = options(scipen = 200, stringsAsFactors = FALSE, "warn" = -1) # suppress warnings
    df_list = df_tbl = NULL
    if (parallel) {
        parallel <- start_parallel_computing(parallel)
        stopCluster <- TRUE
    } else {
        parallel <- stopCluster <- FALSE
    }
    on.exit(if (parallel & stopCluster) stop_parallel_computing(attr(parallel, "cluster")))

    if (is.null(func)) {
        stop("The function is missing. Please input a split_type function.\n")
    }
    if (is.null(x_list)) {
        x_list = get_names(data, types = c("character", "factor", "numeric", "integer", "double"))
        warning(paste("x_list is NULL, use all variables.\n"))
    }
    i. = NULL
    if (!parallel) {
        funct = function(i.) {
            try(do.call(func, c(args, x = i.)), silent = FALSE)
        }
        df_list <- lapply(x_list, funct)
        if (as_list) {
            df_tbl <- df_list
            names(df_tbl) <- x_list
        } else {
            df_tbl <- as.data.frame(Reduce(bind, df_list))
        }
    } else {
        df_list <- foreach(i. = x_list, .errorhandling = c('pass')) %dopar% {
            try(do.call(func, args = c(args, x = i.)), silent = FALSE)
        }
        if (as_list) {
            df_tbl <- df_list
            names(df_tbl) <- x_list
        } else {
            df_tbl <- as.data.frame(Reduce(bind, df_list))
        }
    }
    options(opt) # reset warnings
    return(df_tbl)
}
