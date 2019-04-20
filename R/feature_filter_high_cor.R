#' high_cor_filter
#'
#'
#' \code{fast_high_cor_filter} In a highly correlated variable group, select the  variable with the highest IV.
#' \code{high_cor_filter} In a highly correlated variable group, select the  variable with the highest IV.
#' @param dat A data.frame with independent variables.
#' @param p  Threshold of correlation between features. Default is 0.7.
#' @param x_list Names of independent variables.
#' @param com_list   A data.frame with important values of each variable. eg : IV_list
#' @param ex_cols A list of excluded variables. Regular expressions can also be used to match variable names. Default is NULL.
#' @param cor_class  Culculate catagery variables's correlation matrix. Default is FALSE.
#' @param parallel Logical, parallel computing. Default is FALSE.
#' @param onehot one-hot-encoding independent variables.
#' @param note  Logical. Outputs info. Default is TRUE.
#' @param save_data Logical, save results in locally specified folder. Default is FALSE
#' @param file_name  The name for periodically saved results files. Default is "Feature_selected_COR".
#' @param dir_path The path for periodically saved results files. Default is "./variable".
#' @param ...  Additional parameters.
#' @return  A list of selected variables.
#' @seealso \code{\link{get_correlation_group}}, \code{\link{reduce_high_cor}}, \code{\link{char_cor_vars}}
#' @examples
#' \dontrun{
#' # calculate iv for each variable.
#' iv_list = feature_select_wrapper(dat_train = UCICreditCard, dat_test = NULL, 
#' target = "default.payment.next.month", 
#' occur_time = "apply_date", 
#' filter = c("IV"), cv_folds = 1, iv_cp = 0.01, 
#' psi_cp = 0.1, xgb_cp = 0, cor_cp = 0.98, 
#' ex_cols = "ID$|date$|default.payment.next.month$", 
#' save_data = FALSE, vars_name = FALSE)
#' fast_high_cor_filter(dat = UCICreditCard, 
#' com_list = iv_list, save_data = TRUE, 
#' ex_cols = "ID$|date$|default.payment.next.month$",
#' p = 0.6, cor_class = FALSE ,var_name = FALSE)
#' }
#' @export

fast_high_cor_filter <- function(dat, p = 0.7, x_list = NULL, com_list = NULL,
ex_cols = NULL, save_data = FALSE, cor_class = TRUE, parallel = FALSE,note = TRUE, 
file_name ="Feature_selected_COR" ,dir_path = "./variable", ...) {
    if(note)cat("[NOTE] Fast dimension reduction for highly correlated variables. \n")
    dat = checking_data(dat)
    dat  =  time_transfer(dat)  
    if (!is.null(x_list)) {
        dat = dat[, unique(c(x_list))]
    }
    dat <- low_variance_filter(dat,lvp = 1)
    ex_x_list = get_names(dat = dat, types = c('factor', 'character', 'numeric', 'integer', 'double'),
                          ex_cols = ex_cols, get_ex = TRUE)

    if (is.null(com_list)) {
        stop("The comparison list is empty. For comparisons between variables, IV or PSI or other index must be used.")
    }
    char_x_list = num_x_list = NULL
    if (cor_class) {
        char_x_list = get_names(dat = dat, types = c('factor', 'character'),
        ex_cols = ex_cols, get_ex = FALSE)
        num_x_list = get_names(dat = dat, types = c('numeric', 'integer', 'double'),
        ex_cols = ex_cols, get_ex = FALSE)
        if (length(num_x_list) > 2) {
            cor_mat_num = cor(dat[num_x_list], method = "spearman", use = "complete.obs")
            cor_nums <- reduce_high_cor(cor_mat = cor_mat_num, p = p, com_list = com_list, x_list = num_x_list)
        } else {
            cor_nums = num_x_list
        }
        if (length(num_x_list) > 2) {
            cor_mat_char = char_cor(dat = dat, x_list = char_x_list, parallel = parallel)
            cor_chars <- reduce_high_cor(cor_mat = cor_mat_char, p = p, com_list = com_list, x_list = char_x_list)
        } else {
            cor_chars = num_x_list
        }
        cor_vars <- unique(c(cor_chars, cor_nums))
        var_list = cor_vars    
    } else {
        char_x_list = get_names(dat = dat, types = c('factor', 'character'), ex_cols = ex_cols, get_ex = FALSE)
        num_x_list = get_names(dat = dat, types = c('numeric', 'integer', 'double'), ex_cols = ex_cols, get_ex = FALSE)
        if (length(num_x_list) > 2) {
            cor_mat_num = cor(dat[num_x_list], method = "spearman", use = "complete.obs")
            cor_vars = reduce_high_cor(cor_mat = cor_mat_num, p = p, com_list = com_list, x_list = num_x_list)
        } else {
            cor_vars = num_x_list
        }
       
        var_list = c(char_x_list, cor_vars)
    }
       
    if (save_data) {
        save_dt(var_list, file_name = file_name, dir_path = dir_path, note = note, as_list = TRUE)
    }
    return(var_list)
}



#' @rdname fast_high_cor_filter
#' @export

high_cor_filter <- function(dat, com_list = NULL, x_list = NULL, ex_cols = NULL,
onehot = TRUE, parallel = TRUE, p = 0.7,
dir_path = "./vars",  save_data = TRUE, note = TRUE, ...) {
    if(note)cat("[NOTE] Dimension reduction for highly correlated variables. \n")
    dat = checking_data(dat = dat)
    dat <- time_transfer(dat)
    dat <- merge_category(dat, note = FALSE)

    if (!is.null(x_list)) {
        dat = dat[, unique(c(x_list))]
    }
    if (onehot) {
        #if one-hot of charactor of factor variables.
        dat <- one_hot_encoding(dat)
    }

    #obtain the exclueded variables.
    ex_list <- get_names(dat = dat, types = c('factor', 'character', 'numeric', 'integer', 'double'), ex_cols = ex_cols, get_ex = TRUE)
    #obtain the numeric variables.
    num_x_list = get_names(dat = dat, types = c('numeric', 'integer', 'double'), ex_cols = ex_cols, get_ex = FALSE)
    #obtain the character or factor variables.
    char_x_list = get_names(dat = dat, types = c('factor', 'character'), ex_cols = ex_cols, get_ex = FALSE)

    if (note) cat("[NOTE] Calculate the correlation matrix of numeric variables. \n")
    cor_mat_num = cor(dat[num_x_list], method = "spearman")
    #calculate the correlation matrix of character or factor variables.
    cor_mat_char = char_cor(dat = dat, x_list = char_x_list, parallel = parallel)
    # obtain highly correlated variable groups.
    group_vars <- c(get_correlation_group(cor_mat_num, p = p), get_correlation_group(cor_mat_char, p = round(p / 1.5, 1)))
    group_len <- sapply(group_vars, function(x) length(x))
    single_group_vars <- unlist(group_vars[group_len == 1])
    multi_group_vars <- group_vars[group_len > 1]
    cat("[NOTE] Selecting the variable with the highest IV in a highly correlated variable group . \n")
    x = multi_group_vars[[1]]
    sel_vars <- vapply(multi_group_vars, function(x) {

        #In a highly correlated variable group, the variable with the highest IV  was selected.
        if (!is.null(com_list) & all(x%in% as.character(com_list[,1]))) {
            x_group = com_list[which(as.character(com_list[,1])%in% x),]
            goup_max = x_group[which.max(x_group[, 2]), 1]
        } else {
            #If any variable in a group is not in the comparison list, or the comparison list is missing, the variable with the smallest average correlation coefficient was selected.
                if (any(x %in% num_x_list)) {
                    min_cor <- colMeans(cor_mat_num)[which(colnames(cor_mat_num) %in% x)]
                    goup_max = names(which.min(min_cor))
                } else {
                    min_cor <- colMeans(cor_mat_char)[which(colnames(cor_mat_char) %in% x)]
                    goup_max = names(which.min(min_cor))
                }           
        }
        return(goup_max)
    }, FUN.VALUE = character(1))
    cor_vars = c(single_group_vars, sel_vars)
    dat <- dat[cor_vars]
    #return to the original form of one-hot encoding variables
    dat <- de_one_hot_encoding(dat)
    var_list = colnames(dat)
    if (save_data) {
        save_dt(var_list, file_name = "vars_high_cor_filter", dir_path = dir_path, note = note, as_list = TRUE)
        save_dt(group_vars, file_name = "vars_cor_group", dir_path = dir_path, note = note, as_list = TRUE)
    }
    return(var_list)

}


#' Compare the two highly correlated variables
#'
#' \code{reduce_high_cor} is function for comparing the two highly correlated variables, select a variable with the largest IV value.
#'
#' @param cor_mat A correlation matrix.
#' @param p  The threshold of high correlation.
#' @param x_list Names of independent variables.
#' @param com_list  A data.frame with important values of each variable. eg : IV_list.
#' @param retain Logical, output selected variables, if FALSE, output filtered variables.
#' @return  A list of selected variables.
#' @export
reduce_high_cor <- function(cor_mat, p = 0.90, x_list = NULL, com_list = NULL, retain = TRUE) {
    cols = NULL
    if (!is.null(cor_mat) & !is.null(com_list)) {
        x_com_list = com_list[, "Feature"]
        x_cor_mat = rownames(cor_mat)
        x_list = intersect(intersect(com_list[, "Feature"], rownames(cor_mat)), x_list)
        cor_mat <- cor_mat[x_list, x_list]
        vars_num <- dim(cor_mat)[1]
        if (length(vars_num) > 0 && vars_num > 2) {
            if (!isTRUE(all.equal(cor_mat, t(cor_mat)))) stop("correlation matrix is not symmetric")
            cor_mat <- abs(cor_mat)
            delete_cols <- rep(FALSE, vars_num)
            cor_mat2 <- cor_mat
            diag(cor_mat2) <- NA
            IV_t <- t(com_list)
            colnames(IV_t) <- IV_t[1,]
            IV_cor <- IV_t[2,]

            for (i in 1:(vars_num - 1)) {
                if (!any(cor_mat2[!is.na(cor_mat2)] > p)) {
                    break
                }
                if (delete_cols[i]) next
                for (j in (i + 1):vars_num) {
                    if (!delete_cols[i] & !delete_cols[j]) {
                        if (cor_mat[i, j] > p) {
                            iv1 <- as.numeric(IV_cor[colnames(cor_mat)[i]])
                            iv2 <- as.numeric(IV_cor[colnames(cor_mat)[j]])
                            if (!is.na(iv1) & !is.na(iv2) & length(iv1) > 0 & length(iv2) > 0) {
                                if (iv1 <= iv2) {
                                    delete_cols[i] <- TRUE
                                    cor_mat2[i,] <- NA
                                    cor_mat2[, i] <- NA
                                }
                                else {
                                    delete_cols[j] <- TRUE
                                    cor_mat2[j,] <- NA
                                    cor_mat2[, j] <- NA
                                }
                            }
                        }
                    }
                }
            }
            if (retain) {
                cols = colnames(cor_mat2[, which(!delete_cols)])
            } else {
                cols = colnames(cor_mat2[, which(delete_cols)])
            }
        } else {
            if (retain) {
                cols = colnames(cor_mat)
            }
        }
    }
    cols
}


#' get_correlation_group
#'
#'
#' \code{get_correlation_group} is funtion for  obtaining highly correlated variable groups.
#' \code{select_cor_group} is funtion for selecting highly correlated variable group.
#' \code{select_cor_list} is funtion for selecting highly correlated variable list.
#' @param cor_mat  A correlation matrix of independent variables.
#' @param p  Threshold of correlation between features. Default is 0.7.
#' @return  A list of selected variables.
#' @examples
#' \dontrun{
#' cor_mat = cor(UCICreditCard[8:20], 
#' use = "complete.obs", method = "spearman")
#' get_correlation_group(cor_mat, p = 0.6 )
#' }
#' @export

get_correlation_group <- function(cor_mat, p = 0.8) {
    cat("[NOTE] Getting highly correlated groups of variables. \n")
    vars_num <- dim(cor_mat)[1]
    cor_vars_list = correlation_sub = cor_vars_list_final = cor_arr = NULL

    if (length(vars_num) > 0 && vars_num > 1) {
        diag(cor_mat) <- NA
        if (!any(abs(cor_mat)[!is.na(abs(cor_mat))] > p)) {
            cor_vars_list = colnames(cor_mat)
        } else {
            correlation_sub <- data.frame(which(abs(cor_mat) > p, arr.ind = TRUE))
            correlation_sub <- subset(correlation_sub, col != row)
            cor_arr <- list()
            for (i in unique(correlation_sub$col)) {
                cor_arr[[i]] <- sort(unique(unlist(correlation_sub[which(correlation_sub$col == i),])))
            }
            cor_vars = unique(cor_arr[!sapply(cor_arr, function(x) is.null(x))])
            cor_vars_list_final = select_cor_group(cor_vars)
            cor_vars_list <- list()
            cor_vars_list <- lapply(1:length(cor_vars_list_final),
                              function(x) colnames(cor_mat[, unlist(cor_vars_list_final[x])]))
                              cor_vars_list = append(cor_vars_list,
            colnames(cor_mat)[which(!(colnames(cor_mat) %in% unlist(cor_vars_list)))])
        }
    } else {
        cor_vars_list = colnames(cor_mat)
    }
    return(cor_vars_list)
}



#' @param cor_vars  Correlated variables.
#' @rdname get_correlation_group
#' @export

select_cor_group <- function(cor_vars) {
    cor_vars_group = list()
    cor_vars_group_final = cor_vars
    for (i in 1:length(cor_vars)) {
        cor_vars_group[[i]] = select_cor_list(cor_vars_group_final)
        if (length(cor_vars_group[[i]]) == length(cor_vars_group_final)) break
        cor_vars_group_final = cor_vars_group[[i]]
    }
    return(cor_vars_group_final)
}


#' @param cor_vars_list  List of correlated variable
#' @rdname get_correlation_group
#' @export
select_cor_list <- function(cor_vars_list) {
    cor_vars_list2 <- list()
    for (i in 1:length(cor_vars_list)) {
        cor_vars_list2[[i]] <- lapply(cor_vars_list, function(x) base :: setdiff(unlist(x), unlist(cor_vars_list[i])))
    }
    n_list = length(cor_vars_list2)
    cor_vars_list3 <- unique(lapply(1:n_list,
                            function(i) {
                                cor_vars_sub = cor_vars_list2[[i]]
                                n_vars = length(cor_vars_sub)
                                ind = sapply(1:n_vars, function(i) length(cor_vars_sub[[i]]) != length(cor_vars_list[[i]]))
                                unique(sort(unlist(cor_vars_list[ind])))
                            }))
    return(cor_vars_list3)
}


#' Categery variables' correlation matrix by Cremers'V.
#'
#' \code{char_cor_vars} is function for calculating the correlation coefficent between a Variable and other variables by Cremers'V.
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
char_cor_vars <- function(dat, x ) {
    vapply(seq_along(dat), function(j) {
        if (length(x) > 1| length(unlist(x)) >1) {
            cross_table <- table(unlist(x), dat[, j])         
        } else {
            cross_table <- table(dat[, x], dat[, j])
        }    
        sqrt(chisq.test(cross_table, correct = T,
        simulate.p.value = TRUE)$statistic / (sum(cross_table) * min(ncol(cross_table) - 1, nrow(cross_table) - 1)))
    }, FUN.VALUE = numeric(1))
}

#' @rdname char_cor_vars
#' @export

char_cor <- function(dat,  x_list = NULL, ex_cols = "date$", parallel = FALSE, note = FALSE) {
    if (note) {
        cat("[NOTE] Computing the correlation matrix of factor or character variables.\n")
    }
    if (is.null(x_list)) {
        #Obtaining factor or character variables
        x_list = get_names(dat = dat, types = c('factor', 'character'), ex_cols = ex_cols, get_ex = FALSE)
    }
    if (length(x_list) > 0) {
        #calculating the correlation coefficient between variables-cremers'V
        character_cor = loop_function(func = char_cor_vars, x_list = x_list,
        args = list(dat = dat[x_list]), bind = "cbind", as_list = FALSE, parallel = parallel)
        character_cor <- as.matrix(character_cor)
        colnames(character_cor) = rownames(character_cor) = x_list
    } else {
        character_cor = NULL
    }
    return(character_cor)
}

