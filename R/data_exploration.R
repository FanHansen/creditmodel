#' Data Exploration
#'
#' #'The \code{data_exploration} includes both univariate and bivariate analysis and ranges from univariate statistics and frequency distributions, to correlations, cross-tabulation and characteristic analysis. 
#' @param dat A data.frame with x and target.
#' @param target The name of target variable.
#' @param save_data Logical. If TRUE, save  files to the specified folder at \code{dir_path}
#' @param dir_path The path for periodically saved outliers analysis file. Default is tempdir().
#' @param file_name The file name for periodically saved outliers analysis file. Default is NULL.
#' @return A list contains both categrory and numeric variable analysis.
#' @examples
#' data_ex = data_exploration(dat = lendingclub[1:1000,],target = "loan_status")
#' @importFrom dplyr group_by mutate summarize  summarise n  count %>% filter mutate_if
#' @importFrom data.table fwrite melt fread dcast
#' @export 

data_exploration <- function(dat, target = NULL, save_data = FALSE, file_name = NULL, dir_path = tempdir()) {
    opt = options("warn" = -1, scipen = 100, stringsAsFactors = FALSE, digits = 2) # suppress warnings
    view_char = view_char1 = view_num1 = view_num = Feature =  NULL
    dat1 = dat %>% checking_data(target = target) %>% null_blank_na() %>% time_transfer() %>% mutate_if(is.character, as.factor)
    char_x_list = get_names(dat = dat1,
                                types = c('factor', 'character', "ordered"),
                                ex_cols = target,
                                get_ex = FALSE)
    num_x_list = get_names(dat = dat1,
                               types = c('numeric', 'integer', 'double', "Date",
                               "POSIXlt", "POSIXct", "POSIXt"),
                               ex_cols = target,
                               get_ex = FALSE)
    cat(paste(paste("Observations : ", nrow(dat1)), paste("Numeric variable: ", length(num_x_list)), paste("Category variable: ", length(char_x_list)), sep = "\n"))
    #numeric varaibles
    if (length(num_x_list) > 0) {
        #summary
        dat_num = summary(dat1[num_x_list])
        dat_num = as.data.frame(dat_num)
        view_num = cbind(as.character(dat_num[, "Var2"]), t(as.data.frame(strsplit(dat_num[, "Freq"], ":"))))
        view_num = apply(view_num, 2, function(i) gsub(" ", "", i))
        colnames(view_num) = c("Feature", "Summary", "View")

        view_num[, "Summary"] = ifelse(is.na(view_num[, "Summary"]) | view_num[, "Summary"] %alike% c("NA's", "Other"), "NMiss", view_num[, "Summary"])
        #reshape data
        view_num = dcast(setDT(as.data.frame(view_num)), Feature ~ Summary, value.var = "View")
        if (!is.element("NMiss", names(view_num))) {
            view_num$NMiss = 0
        }
        view_num[is.na(view_num)] = 0
        num_var_list = gsub(" ", "", as.character(unlist(view_num[, "Feature"])))
        view_num$Std = colSds(dat1[num_var_list], na.rm = TRUE)
        if (!is.null(target)) {
            view_num$Corr = cor(dat1[num_var_list], y = dat1[, target],
                    method = "spearman", use = "pairwise.complete.obs")
        } else {
            view_num$Corr = rep(NA, length(num_var_list))
        }
        names(view_num) = c("Feature", "25%", "75%", "Max", "Mean", "Median", "Min", "NMiss", "Std", "Corr")
        view_num1 = as.data.frame(view_num[, c("Feature", "NMiss", "Corr", "Max", "Min", "Mean", "Median", "Std", "25%", "75%")])
        if (save_data) {
            dir_path = ifelse(!is.character(dir_path), tempdir(), dir_path)
            if (!dir.exists(dir_path)) dir.create(dir_path)
            if (!is.character(file_name)) { file_name = NULL }
            save_dt(view_num1, file_name = ifelse(is.null(file_name), "feature.exploration.num", paste(file_name, "feature.exploration.num", sep = ".")), dir_path = dir_path, note = TRUE, as_list = TRUE)
        }
    }
    #category varaibles
    if (length(char_x_list) > 0) {

        dat_char = summary(dat1[char_x_list])
        dat_char = as.data.frame(dat_char)
        dat_char = apply(dat_char, 2, function(i) gsub(" ", "", i))
        view_char = cbind(as.character(dat_char[, "Var2"]), as.data.frame(gsub(":", " : ", dat_char[, "Freq"])))
        view_char = as.data.frame(view_char, row.names = NULL)
        names(view_char) = c("Feature", "Value")

        n_len = view_char %>% dplyr::group_by(Feature) %>% dplyr::count() %>% dplyr::ungroup()
        v_name = strsplit(paste0("Value", 1:max(n_len$n), sep = "", collapse = ","), ",")
        view_char$ids = rep(unlist(v_name), length(char_x_list))
        #reshape data
        view_char1 = data.table::dcast(setDT(as.data.frame(view_char)), Feature ~ ids, value.var = c("Value"))
        #correlation & missing count.
        if (length(unique(dat1[, target])) > 10) {
            dat1[, target] = split_bins(dat = dat1, x = target, bins_no = FALSE)
        }
        dat1[char_x_list] = merge_category(dat1[char_x_list], m = 10)
        Corr_dt = data.frame(Feature = char_x_list, NMiss = unlist(lapply(char_x_list, function(x) sum(is.na(dat1[, x])))))
        view_char1 = merge(view_char1, Corr_dt)
        view_char1[is.na(view_char1)] = "--"
        view_char1 = as.data.frame(view_char1)
        view_char1 = view_char1[, c("Feature", "NMiss", unlist(v_name))]
        if (save_data) {
            dir_path = ifelse(!is.character(dir_path), tempdir(), dir_path)
            if (!dir.exists(dir_path)) dir.create(dir_path)
            if (!is.character(file_name)) { file_name = NULL }
            save_dt(view_char1, file_name = ifelse(is.null(file_name), "feature.exploration.char", paste(file_name, "feature.exploration.char", sep = ".")), dir_path = dir_path, note = TRUE, as_list = TRUE)
        }
    }

    return(list(sum_num = view_num1, sum_char = view_char1))
    options(opt)
}


