#' Data Cleaning
#'
#'
#'The \code{data_cleansing} function is a simpler wrapper for data cleaning functions, such as delete variables that values are all NAs;checking dat and target format.;delete low variance variables.;replace null or NULL or blank with NA; encode variables which NAs &  miss value rate is more than 95% as 1,0 ;encode variables which unique value  rate is  more than 95% as 1,0; merge categories of character variables that  is more than 8; transfer time variables to dateformation; remove duplicated observations;process outliers;process NAs.
#' @param dat A data frame with x and target.
#' @param target The name of target variable.
#' @param miss_values  Other extreme value might be used to represent missing values, e.g: -9999, -9998. These miss_values will be encoded to -1 or "Missing".
#' @param x_list  A list of x variables.
#' @param ex_cols  A list of excluded variables. Default is NULL.
#' @param obs_id  The name of ID of observations.Default is NULL.
#' @param occur_time  The name of occur time of observations.Default is NULL.
#' @param pos_flag The value of positive class of target variable, default: "1".
#' @param outlier_proc  Logical, process outliers or not. Default is TRUE.
#' @param missing_proc  Logical, process nas or not. Default is TRUE.
#' @param low_var  Logical, delete low variance variables or not. Default is TRUE.
#' @param one_hot  Logical. If TRUE, one-hot_encoding  of category variables. Default is FASLE.
#' @param parallel  Logical, parallel computing or not. Default is FALSE.
#' @param note Logical. Outputs info. Default is TRUE.
#' @param save_data  Logical, save the result or not. Default is FALSE.
#' @param file_name  The name for periodically saved data file. Default is NULL.
#' @param dir_path The path for periodically saved data file. Default is tempdir().
#'
#' @return A preprocessed data.frame
#' @seealso
#' \code{\link{remove_duplicated}},
#' \code{\link{null_blank_na}},
#' \code{\link{entry_rate_na}},
#' \code{\link{low_variance_filter}},
#' \code{\link{process_nas}},
#' \code{\link{process_outliers}}
#' @examples
#' #data cleaning
#' dat_cl <- data_cleansing(dat = UCICreditCard[1:2000,],
#'                        target = "default.payment.next.month",
#'                        x_list = NULL,
#'                        obs_id = "ID",
#'                        occur_time = "apply_date",
#'                        ex_cols = c("PAY_6|BILL_"),
#'                        outlier_proc = TRUE,
#'                        missing_proc = TRUE,
#'                       one_hot = FALSE,
#'                        low_var = TRUE,
#'                        save_data = FALSE)
#'
#' @importFrom dplyr group_by mutate summarize  summarise n  count %>% filter mutate_if
#' @importFrom data.table fwrite melt fread dcast
#' @export



data_cleansing <- function(dat, target = NULL, x_list = NULL, obs_id = NULL, occur_time = NULL,
                          pos_flag = NULL, miss_values = NULL, ex_cols = NULL,
                          outlier_proc = TRUE, missing_proc = TRUE,
                          low_var = TRUE, one_hot = FALSE,
                          parallel = FALSE, note = FALSE,
                          save_data = FALSE, file_name = NULL, dir_path = tempdir()) {
    #delete variables that values are all nas.
    opt = options(scipen = 200, "warn" = -1, stringsAsFactors = FALSE, digits = 10) # suppress warnings
    dat = checking_data(dat = dat, target = target, pos_flag = pos_flag,
                            occur_time = occur_time, note = note)
    if (save_data) {
        dir_path = ifelse(!is.character(dir_path), tempdir(), dir_path)
        if (!dir.exists(dir_path)) dir.create(dir_path)
        if (!is.character(file_name)) file_name = NULL
        }
    if (note) cat("[NOTE] data cleansing.\n")

    x_list = get_names(dat = dat,
                              types = c('logical', 'factor', 'character', 'numeric',
                                        'integer', 'double', "Date", "POSIXlt", "POSIXct", "POSIXt"),
                              ex_cols = ex_cols, get_ex = FALSE)
    ex_x_cols = ex_x_cols2 = NULL
    ex_x_cols = get_names(dat = dat,
                              types = c('logical', 'factor', 'character', 'numeric',
                                        'integer', 'double', "Date", "POSIXlt", "POSIXct", "POSIXt"),
                              ex_cols = ex_cols, get_ex = TRUE)
    if (length(dat) > 0) {
        if (!is.null(x_list)) {
            dat = dat[unique(c(obs_id, target, occur_time, x_list))]
        }

        dat =  null_blank_na(dat = dat, miss_values = miss_values, note = note)

        #delete variables that all nas.
        if (low_var) {
            #delecte low vaiance variables
            dat = dat[!colAllnas(dat)] %>%
                      low_variance_filter(lvp = 0.98, note = note)
            x_list_2 = get_names(dat = dat,
                              types = c('logical', 'factor', 'character', 'numeric',
                                        'integer', 'double', "Date", "POSIXlt", "POSIXct", "POSIXt"),
                              ex_cols = ex_cols, get_ex = FALSE)
            ex_x_cols2 = setdiff(x_list, x_list_2)
            ex_x_cols = unique(c(ex_x_cols, ex_x_cols2))
        }
        if (length(ex_x_cols) > 0) {
            cat(paste("[NOTE] exclued variables : ", paste(ex_x_cols, collapse = " , "), ".\n",
                      sep = "", collapse = "\n"))
        }
        #transfer time variables to date formation
        dat = time_transfer(dat = dat, date_cols = NULL,
        ex_cols = c(obs_id, target, occur_time, ex_cols),
                        note = note) %>%
          char_to_num(ex_cols = c(obs_id, target, occur_time, ex_cols), note = note) %>%
        remove_duplicated(target = target, obs_id = obs_id, occur_time = occur_time, note = note) %>% entry_rate_na(nr = 0.95, note = note) %>%
        #merge categories of character variables that  is more than 8
        merge_category(ex_cols = c(obs_id, target, occur_time, ex_cols),
                         p = 0.001, m = 20, note = note)
        char_x_list = get_names(dat = dat,
                                types = c('factor', 'character'),
                                ex_cols = c(obs_id, target, occur_time, ex_cols),
                                get_ex = FALSE)
        num_x_list = get_names(dat = dat,
                               types = c('numeric', 'integer', 'double'),
                               ex_cols = c(obs_id, target, occur_time, ex_cols),
                               get_ex = FALSE)
        date_x_list = get_names(dat = dat,
                                types = c("Date", "POSIXlt", "POSIXct", "POSIXt"),
                                ex_cols = c(obs_id, target, occur_time, ex_cols),
                                get_ex = FALSE)
        flag_list = get_names(dat = dat,
                              types = c('logical', 'factor', 'character', 'numeric',
                                        'integer', 'double', "Date", "POSIXlt", "POSIXct", "POSIXt"),
                              ex_cols = c(obs_id, occur_time, target),
                              get_ex = TRUE)
        dat = dat[, c(flag_list, date_x_list, char_x_list, num_x_list)]
       re_x_list = gsub("-|\\$|\\*|\\+|\\?|\\[|\\^|\\{|\\}|\\\\|\\(|\\)|\\|\\)|\\]", "", c(char_x_list, num_x_list))
        if (one_hot & length(char_x_list) > 0) {
            dat = one_hot_encoding(dat = dat, cat_vars = char_x_list,
                                  na_act = FALSE, note = note)
        }
        if (outlier_proc) {
            dat = process_outliers(dat = dat, target = target,
                                   x_list = num_x_list, ex_cols = c(obs_id, occur_time, target),
                                   parallel = parallel, note = note, save_data = save_data, file_name = file_name,dir_path = dir_path)
        }
        if (missing_proc) {
            dat = process_nas(dat = dat, class_var = FALSE, x_list = c(num_x_list, char_x_list),
                              ex_cols = c(obs_id, occur_time, target),
                              default_miss = FALSE, parallel = parallel,
                              method = "median", note = note, save_data = save_data, file_name = file_name, dir_path = dir_path)
        } else {
            dat = process_nas(dat = dat, class_var = FALSE, x_list = c(num_x_list, char_x_list),
                              ex_cols = c(obs_id, occur_time, target),
                              default_miss = TRUE, parallel = parallel,
                              method = "median", note = note, save_data = save_data, file_name = file_name, dir_path = dir_path)
        }
        re_x_list = gsub("-|\\$|\\*|\\+|\\?|\\[|\\^|\\{|\\}|\\\\|\\(|\\)|\\|\\)|\\]", "_", c(char_x_list, num_x_list))
        dat = re_name(dat, oldname = c(char_x_list, num_x_list), newname = re_x_list)
        if (save_data) {
            save_dt(dat, dir_path = dir_path, file_name = ifelse(is.null(file_name), "dat.cl", paste(file_name, "dat.cl", sep = ".")), append = FALSE, note = note)
        }
    }
    options(opt) # reset warnings
    return(dat)
}


#' Remove Duplicated Observations
#'
#' \code{remove_duplicated} is the function to remove duplicated observations
#'
#' @param dat A data frame with x and target.
#' @param target The name of target variable.
#' @param obs_id The name of ID of observations. Default is NULL.
#' @param occur_time The name of occur time of observations.Default is NULL.
#' @param note   Logical.Outputs info.Default is TRUE.
#'
#' @return  A data.frame
#'
#' @examples
#' datss = remove_duplicated(dat = UCICreditCard,
#' target = "default.payment.next.month",
#' obs_id = "ID", occur_time =  "apply_date")
#' @importFrom dplyr group_by mutate summarize  summarise n  count %>% filter mutate_if distinct ungroup
#' @export


remove_duplicated <- function(dat = dat, obs_id = NULL, occur_time = NULL,
                              target = NULL , note = FALSE) {
    if(note)cat("[NOTE]", "remove duplicated observations...\n")
    y_1 = id_1 = apply_date_1 = NULL
    if (!is.null(obs_id)) {
        if (!is.null(target)) {
            if (!is.null(occur_time)) {
                dat[c("id_1", "apply_date_1", "y_1")] = dat[c(obs_id, occur_time, target)]
                if (is.numeric(dat$y_1)) {
                    dat <- dat %>% dplyr::group_by(id_1) %>%
                    dplyr::filter( y_1 == max(y_1) & apply_date_1 == max(apply_date_1)) %>%
                    dplyr::distinct(id_1, .keep_all = TRUE) %>%
                    dplyr::ungroup()
                } else {
                    dat <- dat %>% dplyr::group_by(id_1) %>% base :: as.factor(y_1) %>%
                        dplyr::filter(y_1 == levels(y_1)[which.min(table(y_1))] & apply_date_1 == max(apply_date_1)) %>%
                        dplyr::distinct(id_1, .keep_all = TRUE) %>%
                        dplyr::ungroup()
                }
            } else {
                dat[c("id_1", "y_1")] = dat[c(obs_id, target)]
                if (is.numeric(dat$y_1)) {
                    dat <- dat %>% dplyr::group_by(id_1) %>%
                    dplyr::filter(y_1 == max(y_1) )%>%
                    dplyr::distinct(id_1, .keep_all = TRUE) %>%
                    dplyr::ungroup()
                } else {
                    dat <- dat %>% dplyr::group_by(id_1) %>% as.factor(y_1) %>%
                        dplyr::filter(y_1 == levels(y_1)[which.min(table(y_1))]) %>%
                        dplyr::distinct(id_1, .keep_all = TRUE) %>%
                        dplyr::ungroup()
                }
            }
        } else {
            if (!is.null(occur_time)) {
                dat[c("id_1", "apply_date_1")] = dat[c(obs_id, occur_time)]
                dat <- dat %>% dplyr::group_by(id_1) %>%
                dplyr::filter(apply_date_1 == max(apply_date_1)) %>%
                dplyr::distinct(id_1, .keep_all = TRUE) %>% dplyr::ungroup()
            } else {
                dat[c("id_1")] = dat[c(obs_id)]
                dat <- dat %>% dplyr::group_by(id_1) %>%
                dplyr::distinct(id_1, .keep_all = TRUE) %>% dplyr::ungroup()
            }
        }
        dat[c("id_1", "apply_date_1", "y_1")] = NULL
    }
    quick_as_df(dat)
}



#' Encode NAs
#'
#' \code{null_blank_na} is the function to  replace null ,NULL, blank or other missing vaules with NA.
#'
#' @param dat A data frame with x and target.
#' @param miss_values  Other extreme value might be used to represent missing values, e.g: -9999, -9998. These miss_values will be encoded to -1 or "Missing".
#' @param note   Logical.Outputs info.Default is TRUE.
#'
#' @return  A data.frame
#'
#' @examples
#' datss = null_blank_na(dat = UCICreditCard[1:1000, ], miss_values =list(-1,-2))
#' @export
#' @importFrom dplyr group_by mutate summarize  summarise n  count %>% filter mutate_if
null_blank_na <- function(dat, miss_values = NULL, note = FALSE) {
    if(note)cat("[NOTE]", "replace null or blank or miss_values, with NA\n")
    dat %>%
      mutate_if(is.character, function(x) {
        ifelse(is.null(x) | x == "" | x == "null" | x == "NULL" | x == " " | x %in% miss_values, NA, x)}) %>%
      mutate_if(is.numeric, function(x) {
        ifelse(x %in% miss_values, NA, x)
        })
}


#' Max Percent of Missing Value
#'
#' \code{entry_rate_na} is the function to recode variables which NAs or  miss value rate is more than 95%.
#'
#' @param dat A data frame with x and target.
#' @param nr  The maximum percent of NAs.
#' @param note   Logical.Outputs info.Default is TRUE.
#'
#' @return  A data.frame
#'
#' @examples
#' datss = entry_rate_na(dat = lendingclub[1:1000, ], nr = 0.98)
#' @export


entry_rate_na <- function(dat, nr = 0.95, note = FALSE) {
    if(note)cat("[NOTE]","process NAs & special value rate is more than 95%\n")
    n_row <- nrow(dat)
    NAs_rate <- vapply(dat, function(x) {
        sum(is.na(x)) / n_row
    }, FUN.VALUE = numeric(1))
    dat[which(NAs_rate > nr)] <- lapply(dat[which(NAs_rate > nr)], function(x) {
        if (is.element(class(x), c('numeric', 'integer', 'double'))) {
            ifelse( is.na(x), -1,  0)
        } else {
                x
        }
    })
    return(dat)
}



#'  Filtering Low Variance Variables
#'
#' \code{low_variance_filter} is the function for filtering low variance variables.
#'
#' @param dat A data frame with x and target.
#' @param lvp  The maximum percent of unique values (including NAs).
#' @param note   Logical.Outputs info.Default is TRUE.
#' @return  A data.frame
#'
#' @examples
#' dat = low_variance_filter(lendingclub[1:1000, ], lvp = 0.9)
#'
#' @export
low_variance_filter <- function(dat, lvp = 0.97, note = FALSE) {
    if (note) {
        (cat("[NOTE] delete low variance variables.\n"))
    }
    dat <- dat[which(colSums(is.na(dat)) != nrow(dat))]
    is_na <- as.data.frame(abs(is.na(dat)))
    dat <- dat[which(apply(is_na, 2, function(x) sum(is.na(x)) / nrow(is_na)) < lvp)]
    low_variance <- vapply(dat, function(x) max(table(x, useNA = "always")) / nrow(dat),
                           FUN.VALUE = numeric(1))
    dat <- dat[which(low_variance < lvp)]
    return(dat)
}
