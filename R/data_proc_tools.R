#' require_packages
#'
#' \code{require_packages} is function for librarying required packages and  installing missing packages if needed.
#' @param pkg A list or vector of names of required packages.
#' @return  packages installed and library.
#' @examples
#' \dontrun{
#' require_packages(c("data.table","ggplot2"))
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


#' quick_as_df
#'
#' \code{quick_as_df} is function for fast dat frame  transfromation.
#' @param df_list A list of data.
#' @return  packages installed and library,
#'
#' @examples
#' \dontrun{
#' UCICreditCard = quick_as_df(UCICreditCard)
#' }
#' @export
quick_as_df <- function(df_list) {
    class(df_list) <- "data.frame"
    attr(df_list, "row.names") <- .set_row_names(length(df_list[[1]]))
    df_list
}


#' Save Data.frame or List
#' 
#' \code{save_dt} is for saving a data.frame or a list fast.
#' @param dt A data frame or list.
#' @param file_name A string. The file name to save breaks_list.
#' @param dir_path  A string. The dir path to save breaks_list. 
#' @param note  Logical.Outputs info.Default is TRUE.
#' @param as_list  Logical. List format or data.frame format to save.Default is FALSE.
#' @param row_names  Logical,retain rownames.
#' @param append Logical, append newdata to old.
#' @examples
#' \dontrun{
#' save_dt(dat = UCICreditCard,"UCICreditCard",dir_path = "./data")
#' }
#' @importFrom data.table fwrite fread dcast melt
#' @export

save_dt <- function(dt, file_name = NULL, dir_path = "./data", note = TRUE, as_list = FALSE,
row_names = FALSE, append = FALSE) {

    if (!dir.exists(dir_path)) dir.create(dir_path)
    if (dir.exists(paste0(dir_path, '/', file_name, ".csv"))) file.remove(list.files(paste0(dir_path, '/', file_name, ".csv"), recursive = TRUE, full.names = TRUE))
    if (note) {
        cat(paste0("[NOTE]", "Saved to ", "'",paste0(dir_path, '/', file_name, ".csv"),"'\n"))
    } 
    if (as_list) {
        fwrite(list(dt), paste0(dir_path, '/', file_name, ".csv"), append = append,col.names = FALSE)
    } else {
        fwrite(as.data.frame(dt), paste0(dir_path, '/', file_name, ".csv"), append = append, row.names = row_names)
    }
}

#' digits_num
#' 
#' \code{digits_num} is for caculating optimal digits number for numeric variables.
#' @param num A sample of numeric variable.
#' @return  A number of digits
#' @examples
#' \dontrun{
#' digits_num(3.1415296)
#' }
#' @export

digits_num <- function(num) {
    options(scipen = 100)
    char_num <- as.character(gsub("-", "", num))
    digits1 = 0;
    digits2 = 0;

    n_num = as.numeric(char_num) %% 1

    if (!is.null(n_num) && !is.na(n_num) && is.numeric(n_num)) {
        if (n_num == 0) {
            t_lens = nchar(char_num)
            left_comma <- t_lens
            right_comma <- t_lens - left_comma
        } else {
            comma_p <- gregexpr("[.]", char_num)[[1]][1]
            t_lens <- nchar(char_num)
            left_comma <- comma_p - 1
            right_comma <- t_lens - 1 - left_comma
        }
        if (left_comma >= 3) {
            digits1 = 2 - left_comma
        } else {
            if (left_comma == 2 & right_comma >= 1) {
                digits1 = 1
            } else {
                if (right_comma >= 1) {
                    digits1 = right_comma
                }
            }
        }
    } else {
        digits1 = 2
    }
    digits2 = min(digits1, getOption("digits"))
    digits2
}


#' one-hot encoding
#'
#' \code{one_hot_encoding} is for converting the factor or character variables into multiple columns
#' @param dat A dat frame.
#' @param cat_vars The name or Column index list to be one_hot encoded. 
#' @param merge_cat Logical. If TRUE, to merge categories greater than 8, default is TRUE.
#' @param ex_cols  Variables to be  excluded, use regular expression matching
#' @param na_act Logical,If true, the missing value is processed, if FALSE missing value is omitted .
#' @param note Logical.Outputs info.Default is TRUE.
#' @return A dat frame with the one hot encoding applied to all the variables with type as factor or character.
#' @seealso \code{\link{de_one_hot_encoding}}
#' @examples
#' dat1 = one_hot_encoding(dat = UCICreditCard, 
#' cat_vars = c("SEX", "MARRIAGE"), 
#' merge_cat = TRUE, na_act = TRUE)
#' dat2 = de_one_hot_encoding(dat_one_hot = dat1, 
#' cat_vars = c("SEX","MARRIAGE"), na_act = FALSE)
#' 
#' @export
one_hot_encoding = function(dat, cat_vars = NULL, ex_cols =NULL, merge_cat = TRUE,na_act = TRUE,note = TRUE) {
   if(note) cat("[NOTE]  one-hot encoding for charactor or factor.\n")
    if (class(dat)[1] != "data.frame") {
        dat <- as.data.frame(dat)
    }
    if (is.null(cat_vars)) {
        cat_vars <- get_names(dat = dat, types = c("character", "factor"), ex_cols = ex_cols)
    }
    if (length(cat_vars) > 0) {
        if (na_act) {
            dat[, cat_vars] = process_nas(dat[cat_vars], note = FALSE)
        }
        if (merge_cat) {
            dat[, cat_vars] = merge_category(dat[cat_vars] ,note = FALSE)
        }
        for (i in cat_vars) {
            if (is.factor(dat[, i]) || is.character(dat[, i])) {
                col_name = i
                dat[, i] = sapply(dat[, i], function(x) gsub(" |\"|\\$|\\*|\\+|\\?|\\[|\\^|\\{|\\}|\\\\|\\(|\\)|\\|\\)|\\]|\\.|\\-", "", x))
                cat_list <- unique(dat[, i])
                encode_cols <- length(cat_list)
                #Create individual column for every unique value in the variable
                for (j in 1:encode_cols) {
                    one_hot_name <- (paste(col_name, ".", cat_list[j], ".", sep = ""))
                    dat[, one_hot_name] <- ifelse(dat[, i] == cat_list[j] & !is.na(dat[, i]), 1, 0)
                }
            }
        }
        dat = dat[, - which(names(dat) %in% cat_vars)]
    } 
    return(dat)
}


#' recovery one-hot encoding
#'
#' \code{de_one_hot_encoding} is for one-hot encoding recovery processing
#' @param dat_one_hot A dat frame with the one hot encoding variables
#' @param cat_vars  variables to be recovery processed, default is null, if null, find these variables through regular expressions . 
#' @param na_act Logical,If true, the missing value is  assigned as "Unknown", if FALSE missing value is omitted, the default is TRUE.
#' @param note Logical.Outputs info.Default is TRUE.
#' @return A dat frame with the one hot encoding recorery character variables
#' @seealso \code{\link{one_hot_encoding}}
#' @examples
#' \dontrun{
#' dat1 = one_hot_encoding(dat = UCICreditCard, 
#' cat_vars = c("SEX", "MARRIAGE"),
#' merge_cat = TRUE, na_act = TRUE)
#' dat2 = de_one_hot_encoding(dat_one_hot = dat1, 
#' cat_vars = c("SEX","MARRIAGE"), 
#' na_act = FALSE)
#' }
#' @export

de_one_hot_encoding = function(dat_one_hot, cat_vars = NULL, na_act = TRUE,note = TRUE) {
    if(note)cat("[NOTE] recovery one-hot encoding for charactor or factor.\n")
    if (class(dat_one_hot)[1] != "data.frame") {
        dat_one_hot <- as.data.frame(dat_one_hot)
    }

    if (is.null(cat_vars)) {
        char_names = one_hot_names = one_hot_names = c()
        for (i in 1:length(dat_one_hot)) {
            char_names[i] <- sub(paste0("\\.$"), "", colnames(dat_one_hot)[i])
            if (!is.null(char_names[i]) && !is.na(char_names[i]) && char_names[i] == colnames(dat_one_hot)[i]) {
                char_names[i] <- NA
            }
            one_hot_names[i] <- try(strsplit(char_names[i], "[.]")[[1]][1], silent = TRUE)
        }
        cat_vars <- unique(one_hot_names[!is.na(one_hot_names)])
    }

    one_hot_vars <- unlist(sapply(cat_vars, function(x) grep(paste0(x, "\\.", "\\S{1,100}", "\\."), paste(colnames(dat_one_hot)))))

    de_cat_vars = intersect(cat_vars, unique(gsub("\\d{1}$", "", names(one_hot_vars))))

    if (length(de_cat_vars) > 0) {
        dat_one_hot[, de_cat_vars] <- lapply(de_cat_vars, function(x) {
            grx = cv_cols = names_1 = re_code = NULL
            grx = paste0(x, "\\.", "\\S{1,100}", "\\.$")
            cv_cols =  grep(grx, paste(colnames(dat_one_hot)))
            names_1 = colnames(dat_one_hot)[cv_cols]
            if (na_act) {
                re_code =  rep("other", nrow(dat_one_hot))
            } else {
                re_code  =  rep(NA, nrow(dat_one_hot))
            }
            for (i in 1:(length(names_1))) {
                re_code[which(dat_one_hot[cv_cols][i] == 1)] = strsplit(names_1[i], "[.]")[[1]][2]
            }
            return(re_code)
        }
    )
        names(dat_one_hot[, de_cat_vars]) =  de_cat_vars
        dat_one_hot = data.frame(dat_one_hot, stringsAsFactors = FALSE)[, - one_hot_vars]
    }
    return(dat_one_hot)
}


#' Time Format Transfering
#' 
#' \code{time_transfer} is for transfering time variables to time format.
#' @param dat A data frame
#' @param date_cols  Names of time variable or regular expressions for finding time variables. Default is  "DATE$|time$|date$|timestamp$|stamp$".
#' @param ex_cols Names of excluded variables. Regular expressions can also be used to match variable names. Default is NULL.
#' @param note   Logical, outputs info. Default is TRUE.
#' @return  A data.frame with transfermed time variables.
#' @examples
#' \dontrun{
#' #transfer a variable.
#' dat = time_transfer(dat = lendingclub,date_cols = "issue_d")
#' #transfer a group of variables with similar name.
#' dat = time_transfer(dat = lendingclub,date_cols = "_d$")
#' #transfer all the time variable.
#' dat = time_transfer(dat = lendingclub,date_cols = NULL)
#' }
#' @export

time_transfer <- function(dat, date_cols = "DATE$|time$|date$|timestamp$|stamp$", ex_cols = NULL, note = FALSE) {
    dat <- checking_data(dat)
    if (note) {
        cat("[NOTE] format time variables.\n")
    }
    x_list = get_x_list(x_list = NULL, dat_train = dat, dat_test = NULL, ex_cols = ex_cols)
    date_cols1 =  NULL
    if (!is.null(date_cols)) {
        date_cols1 <- names(dat[x_list])[colnames(dat[x_list]) %islike% date_cols]
    } else {
        date_cols1 = names(dat[x_list])
    }
    df_date = dat[date_cols1]
    df_date <- df_date[!colAllnas(df_date)]
    df_date = df_date[!sapply(df_date, is_date)]
    if (dim(df_date)[2] != 0) {
        df_date_cols <- names(df_date)
        t_sample <- list()
        t_len <- list()
        tryCatch({
        for (x in 1:ncol(df_date)) {
            t_sample[[x]] = min(unlist(lapply(as.character(sample(na.omit(df_date[[x]]), 1)), function(i) {
                if (nchar(i) >= 8) { nchar(i) } else { 0 }
                })), na.rm = T)
           t_len[[x]] = unlist(as.character(sample(na.omit(df_date[[x]]), 1)))
        }
    }, error = function(e) { cat("ERROR :", conditionMessage(e), "\n") }, warning = function(w) { "" })
        date_cols2 = which(t_sample != 0)
        for (x in date_cols2) {
            if (t_sample[[x]] >= 8 & t_sample[[x]] <= 10 & grepl(pattern = "^[2]{1}[0]{1}[0-3]{1}[0-9]{1}-[0-9]{1}-[0-9]{1,2}$|^[2]{1}[0]{1}[0-3]{1}[0-9]{1}-[0-9]{1,2}-[0-3]{1}[0-9]{1}$|^[1]{1}[9]{1}[0-9]{2}-[0-9]{1,2}-[0-3]{1}[0-9]{1}$", x = gsub(" ", "", substr(t_len[[x]], 1, 10)))) {
                df_date[[x]] = as.Date(as.character(df_date[[x]]), "%Y-%m-%d")
            }
            if (t_sample[[x]] > 10 & grepl(pattern = "^[2]{1}[0]{1}[0-3]{1}[0-9]{1}-[0-9]{1}-[0-3]{1}[0-9]{1,2}$|^[2]{1}[0]{1}[0-3]{1}[0-9]{1}-[0-9]{1,2}-[0-3]{1}[0-9]{1}$|^[1]{1}[9]{1}[0-9]{2}-[0-9]{1,2}-[0-3]{1}[0-9]{1}$", x = gsub(" ", "", substr(t_len[[x]], 1, 10))) & grepl(pattern = "^[0-1]{1}[0-9]{1}:[0-9]{2}", gsub(" ", "", substr(t_len[[x]], 12, nchar(t_len[[x]]))))) {
                df_date[[x]] = as.Date(as.character(df_date[[x]]))
            }
            if (t_sample[[x]] >= 7 & t_sample[[x]] <= 9 & grepl(pattern = "^[2]{1}[0]{1}[0-5]{1}[0-9]{1}[0-9]{1}[0-3]{1}[0-9]{1,2}$|^[2]{1}[0]{1}[0-5]{1}[0-9]{1}[0-9]{1,2}[0-3]{1}[0-9]{1}$|^[1]{1}[9]{1}[0-9]{2}[0-9]{1,2}[0-3]{1}[0-9]{1}$", x =t_len[[x]])) {
                df_date[[x]] = as.Date(as.character(df_date[[x]]), "%Y%m%d")
            }
            if (t_sample[[x]] > 10 & grepl(pattern = "^[2]{1}[0]{1}[0-5]{1}[0-9]{1}[0-9]{1}[0-3]{1}[0-9]{1,2}$|^[2]{1}[0]{1}[0-5]{1}[0-9]{1}[0-9]{1,2}[0-3]{1}[0-9]{1}$|^[1]{1}[9]{1}[0-9]{2}[0-9]{1,2}[0-3]{1}[0-9]{1}$", x = gsub(" ", "", substr(t_len[[x]], 1, 10))) & grepl(pattern = "^[0-9]{1,2}:[0-9]{2}", gsub(" ", "", substr(t_len[[x]], 10, nchar(t_len[[x]]))))) {
                df_date[[x]] = as.Date( as.character(df_date[[x]]))
            }
            if (t_sample[[x]] >= 8 & t_sample[[x]] <= 10 & grepl(pattern = "^[2]{1}[0]{1}[0-3]{1}[0-9]{1}/[0-9]{1,2}/[0-9]{1,2}$|^[1]{1}[9]{1}[0-9]{2}/[0-9]{1,2}/[0-9]{1,2}$", x = gsub(" ", "", substr(t_len[[x]], 1, 10)))) {
                df_date[[x]] = as.Date(as.character(df_date[[x]]), "%Y/%m/%d")
            }
            if (t_sample[[x]] > 10 & grepl(pattern = "^[2]{1}[0]{1}[0-3]{1}[0-9]{1}/[0-9]{1,2}/[0-9]{1,2}$|^[2]{1}[0]{1}[0-3]{1}[0-9]{1}/[0-9]{1}/[0-9]{1}$|^[2]{1}[0]{1}[0-3]{1}[0-9]{1}/[0-9]{1,2}/[0-3]{1}[0-9]{1}$|^[2]{1}[0]{1}[0-3]{1}[0-9]{1}/[0-9]{1}/[0-3]{1}[0-9]{1}$|^[1]{1}[9]{1}[0-9]{2}/[0-9]{1,2}/[0-3]{1}[0-9]{1}$", x = gsub(" ", "", substr(t_len[[x]], 1, 10))) & grepl(pattern = "^[0-9]{1,2}:[0-9]{2}", gsub(" ", "", substr(t_len[[x]], 11, nchar(t_len[[x]]))))) {
                df_date[[x]] = as.POSIXct(as.character(df_date[[x]]))
            }
            if (t_sample[[x]] == 10 & grepl(pattern = "^[1]{1}[0-9]{9}", x =t_len[[x]])) {
                df_date[[x]] <- as.POSIXct(as.numeric(df_date[[x]]), origin = "1970-01-01")
            }
            if (t_sample[[x]] == 13 & grepl(pattern = "^[1]{1}[0-9]{12}", x =t_len[[x]])) {
                df_date[[x]] <- as.POSIXct(as.numeric(df_date[[x]]) / 1000, origin = "1970-01-01")
            }
        }
        dat[df_date_cols] <- df_date
        rm(df_date)
    } else {
        dat = dat
    }
    return(dat)
}

#' is_date
#' 
#' \code{is_date} is a small function for distinguishing time formats 
#' @param x  list or vectors
#' @return  A Date.
#' @examples
#' \dontrun{
#' is_date(lendingclub$issue_d)
#' }
#' @export
is_date = function(x) {
    any(class(x) %in% c("Date", "POSIXlt", "POSIXct", "POSIXt"))
}

#' date cut point
#' 
#' \code{date_cut} is  a small function to get date point.
#' @param dat_time  time vectors.
#' @param pct  the percent of cutting. Default: 0.7.
#' @return  A Date.
#' @examples
#' \dontrun{
#' date_cut(dat_time = lendingclub$issue_d, pct = 0.8)
#' }
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


#' Percent Format.
#' 
#' \code{as_percent} is  a small function for making percent format..
#' @param x  A numeric vector or  list.
#' @param digits  Number of digits.Default: 2.
#' @return  x with percent format.
#' @examples
#' \dontrun{
#' as_percent(1)
#' }
#' @export
as_percent <- function(x, digits = 2) {
    x = as.numeric(x)
    pct = round(x, digits) * 100
    x_pct = paste0(pct, ifelse(is.finite(x), "%", ""))
    return(x_pct)
}

#' Recovery Percent Format.
#' 
#' \code{de_percent} is  a small function for recoverying percent format..
#' @param x  Character with percent formant.
#' @param digits  Number of digits.Default: 2.
#' @return  x without percent format.
#' @examples
#' \dontrun{
#' de_percent("100%")
#' }
#' @export

de_percent <- function(x, digits = 2) {
    x = as.character(x)
    round(as.numeric(gsub("%", "", x)) / 100, digits = digits)
}


#' Merge Category
#' 
#' \code{merge_category} is  for merging   category of nominal variables which number of categories is more than m or percent of samples in any categories is less than p.
#' @param dat A data frame with x and target.
#' @param ex_cols A list of excluded variables. Default is NULL.
#' @param p The minimum percent of samples in a category to merge.
#' @param m The minimum number of categories.
#' @param note Logical, outputs info. Default is TRUE.
#' @return  A data.frame with merged category variables.
#' @examples
#' \dontrun{
#' dat =  merge_category(lendingclub,ex_cols = "id$|_d$")
#' }
#' @export

merge_category <- function(dat, ex_cols = "date$|id$|time$|DATA$|ID$|TIME$", p = 0.01, m = 10, note = FALSE) {
    opt = options("warn" = -1) # suppress warnings
    if (note) {
        (cat("[NOTE] merge categories which percent is less than 0.001 or  obs number is less than 10.\n"))
    }
    char_list = get_names(dat = dat, types = c('factor', 'character'), ex_cols = ex_cols, get_ex = FALSE)

    for (x in char_list) {
        dt_x = table(as.character(dat[, x]), useNA = "no")
        merge_cat = which(dt_x < nrow(dat) * p)
        over_vars = order(abs(dt_x), decreasing = TRUE)[m:length(dt_x)]
        char_num = tryCatch({ as.numeric(names(dt_x)) }, error = function(e) { cat("ERROR :", conditionMessage(e), "\n") }, warning = function(w) { as.numeric(names(dt_x)) })
        char_num_ind = which(!is.na(char_num))
        if ((length(merge_cat) > 0 & length(over_vars) > 0) && round(length(char_num_ind) / length(dt_x),2) < 0.8) {
            max_class = unique(c(over_vars, merge_cat))
            dat[which(dat[, x] %in% names(dt_x[max_class])), x] = "other"
        } 
    }
    options(opt) # reset warnings
    return(dat)
}

#' character to number
#' 
#' \code{char_to_num} is  for transfering character variables which are actually numerical to numeric.
#' @param dat A data frame
#' @param ex_cols A list of excluded variables. Regular expressions can also be used to match variable names. Default is NULL.
#' @param note Logical, outputs info. Default is TRUE.
#' @return  A data.frame
#' @examples
#' \dontrun{
#' char_to_num(lendingclub,ex_cols = "id$|_d$")
#' }
#' @export

char_to_num <- function(dat, note = FALSE, ex_cols = "date$|id$|time$|DATA$|ID$|TIME$") {
    opt = options("warn" = -1) # suppress warnings
    if (note) {
        cat("[NOTE] transfer character variables which are actually numerical to numeric.\n")
    }
    char_list = get_names(dat = dat, types = c('factor', 'character'), ex_cols = ex_cols, get_ex = FALSE)

    for (x in char_list) {
        dt_x = table(as.character(dat[, x]), useNA = "no")
        char_num = tryCatch({ as.numeric(names(dt_x)) }, error = function(e) { cat("ERROR :", conditionMessage(e), "\n") }, warning = function(w) { as.numeric(names(dt_x)) })
        char_num_ind = which(!is.na(char_num))
        if (length(char_num_ind) > 0 && length(dt_x) > 5 && round(length(char_num_ind) / length(dt_x),2) >= 0.8) {
            dat[, x] = as.numeric(as.character(dat[, x]))
        }
    }
    options(opt) # reset warnings
    return(dat)
}





#' %islike%
#' Fuzzy String matching
#' 
#' @param x  A string.
#' @param y  A string.
#' @return  Logical.
#' @examples
#' \dontrun{
#'  "xyz"  %islike% "xy"
#' }
#' @export


'%islike%' <- function(x, y) {
    grx = FALSE
    x = gsub("\\$|\\*|\\+|\\?|\\[|\\^|\\{|\\}|\\\\|\\(|\\)|\\|\\)|\\]", "", x)
    y = gsub("\\{|\\}", "", y)

    if (any(x != '') & any(y != '') & any(!is.null(x)) & any(!is.null(y)) & any(!is.na(x)) & any(!is.na(y))) {
        y = unique(unlist(y))
        y = y[which(y != '')]
        if (length(y[!grepl("\\|", y)]) > 0) {
            grx = Reduce(function(u, v) paste(u, v, sep = '|'), paste0("^", y[!grepl("\\|", y)], sep = "$"))
        } 
        if (length(y[grepl("\\|", y)]) > 0) {
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
#' \dontrun{
#'  "xyz"  %alike% "xy"
#' }
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
#' \dontrun{
#'  dt = re_name(dat = UCICreditCard, "default.payment.next.month" , "target")
#' }
#' @export


re_name <- function(dat, oldname = c() , newname = c()) {
    ind <- which(names(dat) == oldname)
    names(dat)[ind] = newname
    dat
}


#' min_max_norm
#' 
#' \code{min_max_norm} Functions for vector operation.
#' @param x  A data.frame or Matrix.
#' @return  A list
#' @examples
#' \dontrun{
#' dat = one_hot_encoding(dat = UCICreditCard, ex_col = "ID$|_date$",merge_cat = TRUE, na_act = TRUE)
#' dat = low_variance_filter(dat, lvp = 0.9)
#' dat_s = apply(dat, 2, min_max_norm)
#' }
#' @export

min_max_norm <- function(x) {
    ((x - min(x)) / (max(x) - min(x)))
}
#' Functions for vector operation.
#' 
#' @param x  A data.frame or Matrix.
#' @param na.rm  Logical, remove NAs.
#' @return  A list
#' @examples
#' \dontrun{
#'  rowAny(UCICreditCard[8:10])
#'  rowAllnas(UCICreditCard)
#' }
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


#' @rdname rowAny
#' @export
colSds <- function(x, na.rm = FALSE) sqrt(colSums((x - colMeans(x, na.rm = na.rm)) ^ 2, na.rm = na.rm) / nrow(x))

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
#' \dontrun{
#' x_list = get_names(dat = UCICreditCard, types = c('factor', 'character'), 
#' ex_cols = c("default.payment.next.month","ID$|_date$"), get_ex = FALSE)
#' x_list = get_names(dat = UCICreditCard, 
#' ex_cols = c("default.payment.next.month","ID$|_date$"), get_ex = FALSE)
#' }
#' @export

get_names <- function(dat, types = c('numeric', 'integer', 'double'), ex_cols = NULL, get_ex = FALSE) {
    if (is.null(types)) {
        stop("types is missing!")
    }
    if (is.null(ex_cols)) {
        sel_names = names(dat)[sapply(dat, class) %islike% types]
        ex_names = names(dat)[!(sapply(dat, class) %islike% types)]
    } else {
        var_names = names(dat)[sapply(dat, class) %islike% types]
        ex_vars = names(dat)[colnames(dat) %islike% ex_cols]
        ex_types = names(dat)[!(sapply(dat, class) %islike% types)]
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
#' \dontrun{
#' x_list = get_x_list(x_list = NULL,dat_train = UCICreditCard, 
#' ex_cols = c("default.payment.next.month","ID$|_date$"))
#' }
#' @export

get_x_list <- function(x_list = NULL, dat_train = NULL, dat_test = NULL, ex_cols = NULL) {
    if (!is.null(dat_train)) {
        if (is.null(x_list) | length(x_list) <1 ) {
            if (is.null(dat_test)) {
                x_list_retain = get_names(dat = dat_train, types = c('character', 'factor', 'numeric', 'integer', 'double'), ex_cols = ex_cols, get_ex = FALSE)
            } else {
                x_list_t = get_names(dat = dat_train, types = c('character', 'factor', 'numeric', 'integer', 'double'), ex_cols = ex_cols, get_ex = FALSE)
                x_list_s = get_names(dat = dat_test, types = c('character', 'factor', 'numeric', 'integer', 'double'), ex_cols = ex_cols, get_ex = FALSE)
                x_list_retain = intersect(x_list_t, x_list_s)
            }
        } else {
            if (is.null(dat_test)) {
                x_list_ts = get_names(dat = dat_train, types = c('character', 'factor', 'numeric', 'integer', 'double'), ex_cols = ex_cols, get_ex = FALSE)
                x_input = x_list %in% x_list_ts
            } else {
                x_list_t = get_names(dat = dat_train, types = c('character', 'factor', 'numeric', 'integer', 'double'), ex_cols = ex_cols, get_ex = FALSE)
                x_list_s = get_names(dat = dat_test, types = c('character', 'factor', 'numeric', 'integer', 'double'), ex_cols = ex_cols, get_ex = FALSE)
                x_list_ts = intersect(x_list_t, x_list_s)
                x_excluded_ts = setdiff(x_list_t, x_list_s)
                if (length(x_excluded_ts) > 0) {
                    cat(paste("[EXCLUDED]", "variables which  are not both in train & test :\n", paste(x_excluded_ts, collapse = " \n"), ".\n",sep = ""))
                }
                x_input = x_list %in% x_list_ts
            }
            if (!all(x_input)) {
                x_retain = x_list[which(x_input)]
                x_excluded = x_list[which(!x_input)]
                cat(paste("[EXCLUDED]", "which  are not in x_list: \n", paste(x_excluded, collapse = " \n"), ".\n",sep = ""))
            } else {
                x_retain = x_list
            }
            x_type = sapply(dat_train[, x_retain], class) %in% c('character', 'factor', 'numeric', 'integer', 'double')

            if (!all(x_type)) {
                x_list_retain = x_retain[which(x_type)]
                x_excluded = x_retain[which(!x_type)]
                cat(paste("[EXCLUDED]", "which  are not numeric or character or factor:\n", paste(x_excluded, collapse = " \n"), ".\n", sep = ""))
            } else {
                x_list_retain = x_retain
            }
        }
    } 
    return(x_list_retain)
}
