#' split_bins
#'
#' \code{split_bins} is  for binning using breaks.
#' @param dat A data.frame with independent variables.
#' @param x  The name of an independent variable.
#' @param breaks  Breaks for binning.
#' @param bins_no Number the generated bins. Default is TRUE.
#' @return  A data.frame with Bined x.
#' @examples
#' bins = split_bins(dat = UCICreditCard,
#' x = "PAY_AMT1", breaks = NULL, bins_no = TRUE)
#' @export


split_bins <- function(dat, x, breaks = NULL,  bins_no = TRUE) {
    digits_x = digits_num(dat[ ,x])
    opt = options(stringsAsFactors = FALSE, digits = digits_x) #
    if (length(breaks) < 1) {
        breaks = get_breaks(dat, x, target = NULL, best = FALSE, equal_bins = TRUE, g = 5,note = FALSE)
    }
    sp_value_num = sp_value_char = NULL
    if (any(c("integer", "numeric", "double") == class(dat[, x])[1])) {
        breaks = sort(unlist(unique(c(-Inf, breaks, Inf))))
        bins_1 = cut(round(dat[, x], digits = digits_x),breaks = unique(breaks), dig.lab = digits_x,
                     ordered = TRUE, include.lowest = FALSE, right = TRUE)
        if (bins_no) {
            bins_0 = paste("0", as.numeric(bins_1), sep = "")
            bins = paste(bins_0, bins_1, sep = ".")
            bins[which(as.numeric(bins_1) >= 10)] =
              paste(as.numeric(bins_1[which(as.numeric(bins_1) >= 10)]),
                    bins_1[which(as.numeric(bins_1) >= 10)], sep = ".")
        } else {
            bins = as.character(bins_1)
        }
    } else {
        breaks = unique(breaks)
        if (any(grepl("\\|", breaks))) {
            breaks_s = strsplit(breaks, "\\|")
        } else {
            breaks_s = breaks
        }
        dat[which(!(dat[, x] %in% unlist(breaks_s))), x] = get_median(dat[, x])
        if (length(breaks_s) > 0) {
            for (i in 1:length(breaks_s)) {
                if (length(which(dat[, x] %in% unlist(breaks_s[[i]]))) > 1) {
                    if (i < 10) {
                        dat[which(dat[, x] %in% unlist(breaks_s[[i]])), x] =
                          paste(ifelse(bins_no, paste0("0", i),""), paste(breaks_s[[i]], collapse = ";"), sep = ".")
                    } else {
                        dat[which(dat[, x] %in% unlist(breaks_s[[i]])), x] =
                          paste(ifelse(bins_no, paste0("0", i), ""), paste(breaks_s[[i]], collapse = ";"), sep = ".")
                    }
                }
            }
        }
        bins = dat[, x]
    }
    options(opt) # reset
    return(bins)
}


#' Generates Best Breaks for Binning
#'
#' \code{get_breaks} is for generating optimal binning for numerical and nominal variables.
#' The \code{get_breaks_all}  is a simpler wrapper for \code{get_breaks}.
#' @param dat A data frame with x and target.
#' @param target The name of target variable.
#' @param sp_values A list of missing values.
#' @param x_list A list of x variables.
#' @param ex_cols A list of excluded variables. Default is NULL.
#' @param pos_flag The value of positive class of target variable, default: "1".
#' @param occur_time The name of the variable that represents the time at which each observation takes place.
#' @param oot_pct  Percentage of observations retained for overtime test (especially to calculate PSI). Defualt is 0.7
#' @param best  Logical, if TRUE, merge initial breaks to get optimal breaks for binning.
#' @param equal_bins  Logical, if TRUE, equal sample size initial breaks generates.If FALSE , tree breaks generates using desison tree.
#' @param g  Integer, number of initial bins for equal_bins.
#' @param tree_control the list of tree parameters.
#' \itemize{
#'   \item \code{p} the minimum percent of observations in any terminal <leaf> node. 0 < p< 1; 0.01 to 0.1 usually work.
#'   \item \code{cp} complexity parameter. the larger, the more conservative the algorithm will be. 0 < cp< 1 ; 0.0001 to 0.0000001 usually work.
#'   \item \code{xval} number of cross-validations.Default: 5
#'   \item \code{max_depth} maximum depth of a tree. Default: 10
#' }
#' @param bins_control the list of parameters.
#' \itemize{
#'   \item \code{bins_num} The maximum number of bins. 5 to 10 usually work. Default: 10
#'   \item \code{bins_pct} The minimum percent of observations in any bins. 0 < bins_pct < 1 , 0.01 to 0.1 usually work. Default: 0.02
#'   \item \code{b_chi} The minimum threshold of chi-square merge. 0 < b_chi< 1; 0.01 to 0.1 usually work. Default: 0.02
#'   \item \code{b_odds} The minimum threshold of  odds merge. 0 < b_odds < 1; 0.05 to 0.2 usually work. Default: 0.1
#'   \item \code{b_psi} The maximum threshold of PSI in any bins. 0 < b_psi < 1 ; 0 to 0.1 usually work. Default: 0.05
#'   \item \code{b_or} The maximum threshold of G/B index in any bins.  0 < b_or < 1 ; 0.05 to 0.3 usually work. Default: 0.15
#'   \item \code{odds_psi} The maximum threshold of Training and Testing G/B index PSI in any bins. 0 < odds_psi < 1 ; 0.01 to 0.3 usually work. Default: 0.1
#'   \item \code{mono} Monotonicity of all bins, the larger, the more nonmonotonic the bins will be.  0 < mono < 0.5 ; 0.2 to 0.4 usually work. Default: 0.2
#'   \item \code{kc} number of cross-validations. 1 to 5 usually work. Default: 1
#' }
#' @param parallel Logical, parallel computing or not. Default is FALSE.
#' @param save_data Logical, save results in locally specified folder. Default is TRUE
#' @param file_name File name that save results in locally specified folder. Default is "breaks_list".
#' @param dir_path Path to save results. Default is "./variable"
#' @param note Logical.Outputs info.Default is TRUE.
#' @param ... Additional parameters.
#'
#' @return  A table containing a list of splitting points for each independent variable.
#' @seealso \code{\link{get_tree_breaks}}, \code{\link{cut_equal}}, \code{\link{select_best_class}}, \code{\link{select_best_breaks}}
#' @examples
#' #controls
#' tree_control = list(p = 0.02, cp = 0.000001, xval = 5, maxdepth = 10)
#' bins_control = list(bins_num = 10, bins_pct = 0.02, b_chi = 0.02, b_odds = 0.1,
#'                    b_psi = 0.05, b_or = 15, mono = 0.2, odds_psi = 0.1, kc = 5)
#' # get categrory variable breaks
#' b <-  get_breaks(dat = UCICreditCard[1:1000,], x = "MARRIAGE",
#'                 target = "default.payment.next.month",
#'                 occur_time = "apply_date",
#'                 sp_values = list(-1, "Missing"),
#'                 tree_control = tree_control, bins_control = bins_control)
#' # get numeric variable breaks
#' b2 <-  get_breaks(dat = UCICreditCard[1:1000,], x = "PAY_2",
#'                  target = "default.payment.next.month",
#'                  occur_time = "apply_date",
#'                  sp_values = list(-1, "Missing"),
#'                  tree_control = tree_control, bins_control = bins_control)
#' # get breaks of all predictive variables
#' b3 <-  get_breaks_all(dat = UCICreditCard[1:1000,], target = "default.payment.next.month",
#'                      x_list = c("MARRIAGE","PAY_2"),
#'                      occur_time = "apply_date", ex_cols = "ID",
#'                      sp_values = list(-1, "Missing"),
#'                     tree_control = tree_control, bins_control = bins_control,
#'                      save_data = TRUE)
#'
#' @export

get_breaks_all <- function(dat, target = NULL, x_list = NULL, ex_cols = NULL,
                           pos_flag = NULL, occur_time = NULL, oot_pct = 0.7,
                           best = TRUE, equal_bins = FALSE, g = 10, sp_values = NULL,
                           tree_control = list(p = 0.05, cp = 0.000001, xval = 5, maxdepth = 10),
                           bins_control = list(bins_num = 10, bins_pct = 0.05,
                                               b_chi = 0.05, b_odds = 0.1,
                                               b_psi = 0.05, b_or = 0.15, mono = 0.3,
                                               odds_psi = 0.2, kc = 1),
                           parallel = FALSE, note = FALSE, save_data = TRUE,
                           file_name = NULL, dir_path = tempdir(), ...) {
    opt = options(stringsAsFactors = FALSE)
    dat <- checking_data(dat = dat, target = target, pos_flag = pos_flag)
    if (note) {
        cat(paste("[NOTE]", "start cutting breaks ...."), "\n")
    }
    if (is.null(x_list)) {
        x_list = get_names(dat = dat, types = c('factor', 'character', 'numeric', 'integer', 'double'),
        ex_cols = c(target, occur_time, ex_cols), get_ex = FALSE)
    }
    breaks_list = loop_function(func = get_breaks,
                                x_list = x_list,
                                args = list(dat = dat, target = target,
                                            occur_time = occur_time, oot_pct = oot_pct,
                                            pos_flag = pos_flag, best = best,
                                            tree_control = tree_control,
                                            sp_values = sp_values, equal_bins = equal_bins,
                                            g = g, bins_control = bins_control, note = note),
                                as_list = TRUE, bind = "cbind", parallel = parallel)
    breaks_list1 = lapply(1:length(breaks_list), function(i) data.frame(Feature = names(breaks_list)[i],
    cuts = t(data.frame(t(breaks_list[[i]]))), row.names = NULL))
    breaks_list2 = as.data.frame(rbindlist(breaks_list1))
    if (save_data) {
        dir_path = ifelse(!is.character(dir_path),
                      tempdir(), dir_path)
        if (!dir.exists(dir_path)) dir.create(dir_path)
        if (!is.character(file_name)) file_name = NULL
        save_dt(breaks_list2, file_name = ifelse(is.null(file_name), "breaks.list", paste(file_name, "breaks.list", sep = ".")),
                dir_path = dir_path, note = note, as_list = FALSE)
    }
    options(opt) # reset
    return(breaks_list2)
}


#' @param x  The Name of an independent variable.
#' @rdname get_breaks_all
#' @export

get_breaks <- function(dat, x, target = NULL, pos_flag = NULL,
                       best = TRUE, equal_bins = FALSE, g = 10,
                       sp_values = NULL, occur_time = NULL, oot_pct = 0.7,
                       tree_control = NULL, bins_control = NULL, note = FALSE, ...) {
    dat = checking_data(dat = dat, target = target, pos_flag = pos_flag, occur_time = occur_time)
    if (equal_bins | is.null(target)) {
        tree_breaks = cut_equal(dat_x = dat[, x], g = g, sp_values = sp_values)
    } else {
        tree_breaks <- get_tree_breaks(dat = dat, x = x, target = target, pos_flag = pos_flag,
            tree_control = tree_control, sp_values = sp_values)
    }
    if (!is.null(target) && best) {
        if (any(c("integer", "numeric", "double") == class(dat[, x]))) {
            breaks <- select_best_breaks(dat = dat, x = x, target = target,
                                             occur_time = occur_time, oot_pct = oot_pct,
                                             breaks = tree_breaks, pos_flag = pos_flag,
                                             sp_values = sp_values,
                                             bins_control = bins_control)
        } else {
            breaks <- select_best_class(dat = dat, x = x, target = target,
                                            occur_time = occur_time, oot_pct = oot_pct,
                                            breaks = tree_breaks, pos_flag = pos_flag,
                                            sp_values = sp_values,
                                            bins_control = bins_control)
        }
    } else {
        breaks = tree_breaks
    }
    if (note) cat(paste0(x, "   :   ", paste(breaks, sep = ",", collapse = ","), sep = "\t\t", collapse = "\t"), "\n")
    return(breaks)
}


#' Getting the breaks for terminal nodes from decision tree
#'
#' \code{get_tree_breaks} is for generating initial braks by decision tree for a numerical or nominal variable.
#' The \code{get_breaks} function is a simpler wrapper for \code{get_tree_breaks}.
#' @param dat A data frame with x and target.
#' @param x  name of variable to cut breaks by tree.
#' @param target The name of target variable.
#' @param sp_values A list of special value. Default: NULL.
#' @param pos_flag The value of positive class of target variable, default: "1".
#' @param tree_control the list of parameters to control cutting initial breaks by decision tree.
#' \itemize{
#'   \item \code{p} the minimum percent of observations in any terminal <leaf> node. 0 < p< 1; 0.01 to 0.1 usually work.
#'   \item \code{cp} complexity parameter. the larger, the more conservative the algorithm will be. 0 < cp< 1 ; 0.0001 to 0.0000001 usually work.
#'   \item \code{xval} number of cross-validations.Default: 5
#'   \item \code{max_depth} maximum depth of a tree. Default: 10
#' }
#' @seealso \code{\link{get_breaks}}, \code{\link{get_breaks_all}}
#' @examples
#' #tree breaks
#' tree_control = list(p = 0.02, cp = 0.000001, xval = 5, maxdepth = 10)
#' tree_breaks = get_tree_breaks(dat = UCICreditCard, x = "MARRIAGE",
#' target = "default.payment.next.month", tree_control = tree_control)
#' @importFrom rpart rpart rpart.control path.rpart
#' @export

get_tree_breaks <- function(dat, x, target, pos_flag = NULL,
                            tree_control = list(p = 0.02, cp = 0.000001, xval = 5, maxdepth = 10),
                            sp_values = NULL) {

    dat = checking_data(dat = dat, target = target, pos_flag = pos_flag)
    opt = options(stringsAsFactors = FALSE, digits = 10) #
    sp_value_char = sp_value_num = tree_breaks = NULL
    x_miss = any(dat[, x] %in% sp_values)
    if (!is.null(sp_values) && x_miss) {
        sp_value_char = unlist(sp_values[sapply(sp_values, is.character)])
        sp_value_num = unlist(sp_values[sapply(sp_values, is.numeric)])
        dat1 = dat[!(dat[, x] %in% sp_values),]
    } else {
        dat1 = dat
    }
    if (length(unique(dat1[, x])) > 1) {
        digits_x = digits_num(dat[, x])
        cp = ifelse(!is.null(tree_control[["cp"]]), tree_control[["cp"]], 0.00001)
        xval = ifelse(!is.null(tree_control[["xval"]]), tree_control[["xval"]], 5)
        maxdepth = ifelse(!is.null(tree_control[["maxdepth"]]), tree_control[["maxdepth"]], 10)
        p = ifelse(!is.null(tree_control[["p"]]), tree_control[["p"]], 0.02)
        trcontrol <- rpart.control(minbucket = round(nrow(dat1) * p),
                               cp = cp,
                               xval = xval,
                               maxdepth = maxdepth)

        Formula <- as.formula(paste(target, names(dat1[x]), sep = ' ~ '))
        set.seed(46)
        fit <- rpart(data = dat1, formula = Formula
                 , control = trcontrol
                 , parms = list(split = "information"))
        if (any(is.element(class(dat1[, x]), c("integer", "numeric", "double")))) {
            if (is.null(fit$splits[, 4])) {
                tree_breaks = cut_equal(dat1[, x], g = 20, sp_values = sp_values)
            } else {
                tree_breaks = sort(fit$splits[, 4])

                tree_breaks = floor(tree_breaks * 10 ^ (digits_x)) / 10 ^ (digits_x)
            }
            if (!is.null(sp_values) && x_miss && length(sp_value_num) > 0) {
                tree_breaks = sort(unique(append(c(tree_breaks, Inf), sp_value_num, 0)))
            }
        } else {
            if (is.null(fit$splits[, 4])) {
                tree_breaks = cut_equal(dat1[, x], g = 10, sp_values = sp_values)
            } else {
                rpart.rules <- path.rpart(fit,
                                  rownames(fit$frame)[fit$frame$var == "<leaf>"],
                                  print.it = FALSE)
                tree_breaks = list()
                for (i in 1:length(rpart.rules)) {
                    tree_breaks[i] = strsplit(gsub(paste0(x, "="), "",
                                         rpart.rules[[i]][length(rpart.rules[[i]])]),
                                    split = ",")
                }
            }
            if (!is.null(sp_values) && x_miss && length(sp_value_char) > 0) {
                tree_breaks = unique(append(sp_value_char, tree_breaks, 0))
            }
        }
    } else {
        tree_breaks = as.list(unique(dat[, x]))
    }
    options(opt) # reset
    rm(dat1)
    return(unique(c(tree_breaks)))
}

#' Generating Initial Equal Size Sample Bins
#'
#' \code{cut_equal} is used to generate initial breaks for equal frequency binning.
#' @param dat_x  A vector of an variable x.
#' @param g numeric, number of initial bins for equal_bins.
#' @param sp_values a list of special value. Default: list(-1, "Missing")
#' @seealso \code{\link{get_breaks}}, \code{\link{get_breaks_all}},\code{\link{get_tree_breaks}}
#' @examples
#' #equal sample size breaks
#' equ_breaks = cut_equal(dat = UCICreditCard[, "PAY_AMT2"], g = 10)
#'
#' @export
#' @importFrom stats aggregate approx

cut_equal <- function(dat_x, g = 10, sp_values = list(-1, "Missing")) {
    digits_x = digits_num(dat_x)
    opt = options(stringsAsFactors = FALSE, digits = digits_x + 1) #
    sp_value_char = sp_value_num = NULL
    dat_x = unlist(dat_x)
    x_miss = any(dat_x %in% sp_values)
    if (!is.null(sp_values) && x_miss) {
        sp_value_char = unlist(sp_values[sapply(sp_values, is.character)])
        sp_value_num = unlist(sp_values[sapply(sp_values, is.numeric)])
        dat_x = dat_x[!(dat_x %in% sp_values)]
    }
    if (any(is.element(class(dat_x), c("integer", "numeric", "double")))) {
        dat_x = round(dat_x, digits_x)
        if (length(unique(dat_x)) < 2) {
            cuts = unique(dat_x)
        } else {
            none_na_num <- sum(!is.na(dat_x))
            if (g < 1) stop("g must be >=1")
            tbl_x = table(dat_x)
            x_unique_value = as.double(names(tbl_x))
            cum_sum = sort(cumsum(tbl_x))
            cuts_sum = approx(cum_sum, x_unique_value, xout = (1:g) * round(none_na_num / g),
            method = "constant", ties = "ordered", rule = 2, f = 1)$y
            n_cuts = table(cuts_sum)
            max_n_cuts = as.numeric(names(n_cuts)[which.max(n_cuts)])
            cuts_unique = unique(cuts_sum)
            cuts = cuts_unique[-which.max(cuts_unique)]
            if (length(cuts_unique) <= 2 & length(x_unique_value) > 1) {
                x_unique_ncuts = x_unique_value[which(x_unique_value != max_n_cuts)]
                min_x_unique_value = x_unique_ncuts[which.min(x_unique_ncuts)]
                cuts = append(cuts, min_x_unique_value, 1)
            }
        }
        cuts = round(cuts, digits = digits_x)
        if (!is.null(sp_values) && x_miss && length(sp_value_num) > 0) {
            cuts = sort(unique(append(c(cuts, Inf), sp_value_num, 0)))
        } else {
            cuts = sort(unique(c(cuts, Inf)))
        }
    } else {
        if (length(unique(dat_x)) < 2) {
            cuts = unique(dat_x)
        } else {
            cuts = as.list(unique(dat_x))
        }
        if (!is.null(sp_values) && x_miss && length(sp_value_char) > 0) {
            cuts = unique(append(sp_value_char, cuts, 0))
        } else {
            cuts = unique(cuts)
        }
    }
    rm(dat_x)
    options(opt) # reset
    return(unique(c(cuts)))
}

#' Generates Best Binning Breaks
#'
#'
#' \code{select_best_class} & \code{select_best_breaks} are  for merging initial breaks of variables using chi-square, odds-ratio,PSI,G/B index and so on.
#' The \code{get_breaks}  is a simpler wrapper for \code{select_best_class} & \code{select_best_class}.
#' @param dat A data frame with x and target.
#' @param target The name of target variable.
#' @param breaks Splitting points for an independent variable. Default is NULL.
#' @param sp_values A list of special value.
#' @param x  The name of variable to process.
#' @param pos_flag The value of positive class of target variable, default: "1".
#' @param occur_time The name of the variable that represents the time at which each observation takes place.
#' @param oot_pct  The percentage of Actual and Expected set for PSI calculating.
#' @param bins_control the list of parameters.
#' \itemize{
#'   \item \code{bins_num} The maximum number of bins. 5 to 10 usually work. Default: 10
#'   \item \code{bins_pct} The minimum percent of observations in any bins. 0 < bins_pct < 1 , 0.01 to 0.1 usually work. Default: 0.02.
#'   \item \code{b_chi} The minimum threshold of chi-square merge. 0 < b_chi< 1; 0.01 to 0.1 usually work. Default: 0.02.
#'   \item \code{b_odds} The minimum threshold of  odds merge. 0 < b_odds < 1; 0.05 to 0.2 usually work. Default: 0.1.
#'   \item \code{b_psi} The maximum threshold of PSI in any bins. 0 < b_psi < 1 ; 0 to 0.1 usually work. Default: 0.05.
#'   \item \code{b_or} The maximum threshold of G/B index in any bins.  0 < b_or < 1 ; 0.05 to 0.3 usually work. Default: 0.15.
#'   \item \code{odds_psi} The maximum threshold of Training and Testing G/B index PSI in any bins. 0 < odds_psi < 1 ; 0.01 to 0.3 usually work. Default: 0.1.
#'   \item \code{mono} Monotonicity of all bins, the larger, the more nonmonotonic the bins will be.  0 < mono < 0.5 ; 0.2 to 0.4 usually work. Default: 0.2.
#'   \item \code{kc} number of cross-validations. 1 to 5 usually work. Default: 1.
#' }
#' @param ... Other parameters.
#' @return A list of breaks for x.
#' @details
#'  The folloiwing is the list of Reference Principles
#' \itemize{
#'      \item 1.The increasing or decreasing trend of variables is consistent with the actual business experience.(The percent of Non-monotonic intervals of which are not head or tail is less than 0.35)
#'      \item 2.Maximum 10 intervals for a single variable.
#'      \item 3.Each interval should cover more than 2% of the model development samples.
#'      \item 4. Each interval needs at least 30 or 1%  1 samples. .
#'      \item 5.Combining the values of blank, missing or other special value into the same interval called Missing.
#'      \item 6.The difference of Chi effect size between intervals should be at least 0.03 or more.
#'      \item 7.The difference of absolute odds ratio between intervals should be at least 0.1 or more.
#'      \item 8.The difference of 1 rate between intervals should be at least 1/10 of the total 1 rate..
#'      \item 9.The difference of G/B index between intervals should be at least 15 or more.
#'      \item 10.The PSI of each interval should be less than 0.05.
#'   }
#' @seealso
#' \code{\link{get_tree_breaks}},
#' \code{\link{cut_equal}},
#' \code{\link{get_breaks}}
#' @examples
#' #equal sample size breaks
#' equ_breaks = cut_equal(dat = UCICreditCard[, "PAY_AMT2"], g = 10)
#'
#' # select best bins
#' bins_control = list(bins_num = 10, bins_pct = 0.02, b_chi = 0.02,
#' b_odds = 0.1, b_psi = 0.05, b_or = 0.15, mono = 0.3, odds_psi = 0.1, kc = 1)
#' select_best_breaks(dat = UCICreditCard, x = "PAY_AMT2", breaks = equ_breaks,
#' target = "default.payment.next.month", occur_time = "apply_date",
#' sp_values = NULL, bins_control = bins_control)
#' @export


select_best_class <- function(dat, x, target, breaks = NULL, occur_time = NULL, oot_pct = 0.7,
                              pos_flag = NULL, bins_control = NULL, sp_values = NULL, ...) {
    dat = checking_data(dat = dat, target = target, pos_flag = pos_flag)
    if (is.null(breaks) || any(is.na(breaks)) || length(breaks) < 1) {
        stop("breaks is missing.\n")
    }
    if (!is.character(dat[, x])) {
        stop(paste(x, "must be a character.\n"))
    }
    if (length(breaks) > 2) {
        break_class = sp_value_char = NULL
        x_miss = any(dat[, x] %in% sp_values)
        if (!is.null(sp_values) && x_miss) {
            sp_value_char = unlist(sp_values[sapply(sp_values, is.character)])
            miss_class = unlist(breaks[sapply(breaks, function(x) any(sp_value_char %in% x))])
            non_miss_class = breaks[!sapply(breaks, function(x) any(sp_value_char %in% x))]
            breaks = unique(non_miss_class)
            dat = dat[dat[, x] %in% unlist(breaks),]
        } else {
            breaks = unique(breaks)
        }
        b_chi = ifelse(!is.null(bins_control[["b_chi"]]), bins_control[["b_chi"]], 0.02) / 2
        b_odds = ifelse(!is.null(bins_control[["b_odds"]]), bins_control[["b_odds"]], 0.1) / 2
        bins_num = ifelse(!is.null(bins_control[["bins_num"]]), bins_control[["bins_num"]], 10)
        bins_pct = ifelse(!is.null(bins_control[["bins_pct"]]), bins_control[["bins_pct"]], 0.02)
        b_psi = ifelse(!is.null(bins_control[["b_psi"]]), bins_control[["b_psi"]], 0.05) * 1.5
        b_or = ifelse(!is.null(bins_control[["b_or"]]), bins_control[["b_or"]], 0.1) / 2
        odds_psi = ifelse(!is.null(bins_control[["odds_psi"]]), bins_control[["odds_psi"]], 0.1) * 1.5
        kc = ifelse(!is.null(bins_control[["kc"]]), bins_control[["kc"]], 5)
        if (!is.null(kc) && kc > 1) {
            cv_list = cv_split(dat, k = kc, occur_time = occur_time, seed = 46)
        } else {
            cv_list = cv_split(dat, k = 1 / (1 - oot_pct), occur_time = occur_time, seed = 46)
            kc = 1
        }
        breaks_cv = iv_list_sub = psi_list_sub = list()
        for (k in 1:kc) {
            dat_train = dat[-cv_list[[k]],]
            dat_test = dat[cv_list[[k]],]
            break_class = breaks
            dt_bins = NULL
            while (TRUE) {
                if (length(unique(break_class)) <= 1 | length(dat_train) < 1 | length(dat_test) < 1) break
                dt_bins = get_psi_iv(dat = dat_train, x = x, dat_test = dat_test, target = target,
                pos_flag = pos_flag, breaks = break_class, breaks_list = NULL, occur_time = occur_time,
                oot_pct = oot_pct, bins_total = FALSE, note = FALSE, bins_no = TRUE)
                gb = dt_bins[, c("expected_0", "expected_1")]
                cut_psi = dt_bins[, "PSIi"]
                bins_odds_ratio_s = dt_bins[, "odds_ratio_s"]
                gb_index = dt_bins[, "odds_ratio"]
                gb_percent = dt_bins[, "%expected"]
                dif_gb = c()
                for (brk in 1:(dim(gb)[1] - 1)) {
                    cross_table = rbind(gb[brk,], gb[brk + 1,])
                    a = cross_table[1, 1]
                    b = cross_table[1, 2]
                    c = cross_table[2, 1]
                    d = cross_table[2, 2]
                    if (any(is.na(cross_table)) | any(cross_table < 30) |
                        any(cross_table[, 2] / (cross_table[, 1] + cross_table[, 2]) < 0.01)) {
                        dif_gb[brk] = 0
                    } else {
                        dif_gb[brk] = gb_index[brk] - gb_index[brk + 1]
                    }
                }
                if (length(break_class) <= 2) {
                    break
                } else {
                    if (any(abs(dif_gb) < b_or)) {
                        if (which.min(abs(dif_gb)) < length(break_class)) {
                            break_class[[which.min(abs(dif_gb)) + 1]] = append(unlist(break_class[which.min(abs(dif_gb))]),
                            unlist(break_class[which.min(abs(dif_gb)) + 1]), 0)
                        } else {
                            break_class[[which.min(abs(dif_gb)) - 1]] = append(unlist(break_class[which.min(abs(dif_gb))]),
                            unlist(break_class[which.min(abs(dif_gb)) - 1]), 0)
                        }
                        break_class = break_class[-which.min(abs(dif_gb))]
                    } else {
                        if (length(break_class) > bins_num) {
                            if (which.min(abs(dif_gb)) < length(break_class)) {
                                break_class[[which.min(abs(dif_gb)) + 1]] = append(unlist(break_class[which.min(abs(dif_gb))]),
                                unlist(break_class[which.min(abs(dif_gb)) + 1]), 0)
                            } else {
                                break_class[[which.min(abs(dif_gb)) - 1]] = append(unlist(break_class[which.min(abs(dif_gb))]),
                                unlist(break_class[which.min(abs(dif_gb)) - 1]), 0)
                            }
                            break_class = break_class[-which.min(abs(dif_gb))]
                        } else {
                            if (length(gb_percent) > 2 & any(gb_percent < bins_pct)) {
                                bins_pct_cut = dt_bins[, "cuts"][[which.min(gb_percent)]]
                                max_psi_bin = dt_bins[, "bins"][[which.min(gb_percent)]]
                                min_pct_m = which(dt_bins[, "bins"] %islike% max_psi_bin)

                                if (length(min_pct_m) > 0 && min_pct_m < length(break_class)) {
                                    break_class[[min_pct_m + 1]] = append(bins_pct_cut,
                                                                          break_class[[min_pct_m + 1]], 0)
                                } else {
                                    break_class[[min_pct_m - 1]] = append(bins_pct_cut,
                                                                          break_class[[min_pct_m - 1]], 0)
                                }
                                break_class = break_class[-min_pct_m]
                            } else {
                                if (length(cut_psi) > 2 & any(cut_psi > b_psi)) {
                                    max_psi_cut = dt_bins[, "cuts"][[which.max(cut_psi)]]
                                    max_psi_bin = dt_bins[, "bins"][[which.max(cut_psi)]]
                                    max_psi_m = which(dt_bins[, "bins"] %islike% max_psi_bin)
                                    if (length(max_psi_m) > 0 && max_psi_m < length(break_class)) {
                                        break_class[[max_psi_m + 1]] = append(max_psi_cut,
                                                                              break_class[[max_psi_m + 1]], 0)
                                    } else {
                                        break_class[[max_psi_m - 1]] = append(max_psi_cut,
                                                                              break_class[[max_psi_m - 1]], 0)
                                    }
                                    break_class = break_class[-max_psi_m]
                                } else {
                                    if (length(bins_odds_ratio_s) > 2 & any(bins_odds_ratio_s > odds_psi)) {
                                        max_odds_psi_cuts = dt_bins[, "cuts"][[which.max(bins_odds_ratio_s)]]
                                        max_odds_psi_bin = dt_bins[, "bins"][[which.max(bins_odds_ratio_s)]]
                                        max_odds_psi_m = which(dt_bins[, "bins"] %islike% max_odds_psi_bin)
                                        if (length(max_odds_psi_m) > 0 && max_odds_psi_m < length(break_class)) {
                                            break_class[[max_odds_psi_m + 1]] = append(max_odds_psi_cuts,
                                                                                     break_class[[max_odds_psi_m + 1]], 0)
                                        } else {
                                            break_class[[max_odds_psi_m - 1]] = append(max_odds_psi_cuts,
                                                                                     break_class[[max_odds_psi_m - 1]], 0)
                                        }
                                        break_class = break_class[-max_odds_psi_m]
                                    } else {
                                        break
                                    }
                                }
                            }
                        }
                    }
                }
            }
            dt_iv = get_iv(dat = dat_test, x = x, target = target, pos_flag = pos_flag,
                           breaks = break_class, note = FALSE)
            iv_list_sub[[k]] = unlist(dt_iv[, "IV"])
            if (is.null(dt_bins)) {
                dt_bins = get_psi_iv(dat = dat_train, x = x, dat_test = dat_test,
                                     target = target, pos_flag = pos_flag,
                                     breaks = break_class, breaks_list = NULL,
                                     occur_time = occur_time, oot_pct = oot_pct,
                                     bins_total = FALSE, note = FALSE, bins_no = TRUE)
            }
            psi_list_sub[[k]] = sum(dt_bins[, "PSIi"])
            breaks_cv[[k]] = break_class
        }
        max_iv = unlist(iv_list_sub) / ifelse(max(unlist(iv_list_sub)) == 0, 1, max(unlist(iv_list_sub)))
        min_psi = (1 - unlist(psi_list_sub)) / (1 - min(unlist(psi_list_sub)))
        best_ind = which.max(max_iv * 0.8 + min_psi * 0.2)
        best_class = unique(breaks_cv[[best_ind]])
        if (length(best_ind) > 0) {
            best_class = unique(breaks_cv[[best_ind]])
        } else {
            best_class = unique(unlist(breaks_cv))
        }
        if (!is.null(sp_values) && x_miss && length(miss_class) > 0) {
            best_class = unique(append(c(best_class), miss_class, 0))
        }
    } else {
        best_class = breaks
    }
    return(best_class)
}


#' @rdname select_best_class
#' @export

select_best_breaks <- function(dat, x, target, breaks = NULL, pos_flag = NULL,
                               sp_values = NULL, occur_time = NULL, oot_pct = 0.7,
                               bins_control = NULL, ...) {
    dat = checking_data(dat = dat, target = target, pos_flag = pos_flag)
    if (is.null(breaks) || any(is.na(breaks)) || length(breaks) < 1) {
        stop("breaks is missing")
    }
    if (!any(c("integer", "numeric", "double") == class(dat[, x]))) {
        stop("x must be numeric.")
    }
    if (length(breaks) > 2) {
        break_points = sp_value_num = NULL
        x_miss = any(dat[, x] %in% sp_values)
        if (any(is.na(dat[, x])) | !is.null(sp_values) && x_miss) {
            sp_value_num = unlist(sp_values[sapply(sp_values, is.numeric)])
            miss_num = unlist(breaks[sapply(breaks, function(x) any(sp_value_num %in% x))])
            breaks = breaks[!sapply(breaks, function(x) any(sp_value_num %in% x))]
            dat = dat[!(dat[, x] %in% miss_num | is.na(dat[, x])),]
        }
        b_chi = ifelse(!is.null(bins_control[["b_chi"]]), bins_control[["b_chi"]], 0.01)
        b_odds = ifelse(!is.null(bins_control[["b_odds"]]), bins_control[["b_odds"]], 0.1)
        bins_num = ifelse(!is.null(bins_control[["bins_num"]]), bins_control[["bins_num"]], 10)
        bins_pct = ifelse(!is.null(bins_control[["bins_pct"]]), bins_control[["bins_pct"]], 0.02)
        b_psi = ifelse(!is.null(bins_control[["b_psi"]]), bins_control[["b_psi"]], 0.05)
        b_or = ifelse(!is.null(bins_control[["b_or"]]), bins_control[["b_or"]], 0.15)
        odds_psi = ifelse(!is.null(bins_control[["odds_psi"]]), bins_control[["odds_psi"]], 0.1)
        kc = ifelse(!is.null(bins_control[["kc"]]), bins_control[["kc"]], 5)
        mono = ifelse(!is.null(bins_control[["mono"]]), bins_control[["mono"]], 0.3)
        # kc = 1
        if (!is.null(kc) && kc > 1) {
            cv_list = cv_split(dat, k = kc, occur_time = occur_time, seed = 46)
        } else {
            cv_list = cv_split(dat, k = 1 / (1 - oot_pct), occur_time = occur_time, seed = 46)
            kc = 1
        }
        breaks_cv = iv_list_sub = psi_list_sub = list()

        for (k in 1:kc) {
            dat_train = dat[-cv_list[[k]],]
            dat_test = dat[cv_list[[k]],]
            break_points = unique(unlist(c(breaks, Inf)))
            dt_bins = NULL
            while (TRUE) {
                if (length(unique(break_points)) <= 2 | length(dat_train) < 1 | length(dat_test) < 1) break
                dt_bins = get_psi_iv(dat = dat_train, dat_test = dat_test, x = x, target = target,
                                     pos_flag = pos_flag, breaks = break_points, breaks_list = NULL,
                                     occur_time = occur_time, oot_pct = oot_pct,
                                     bins_total = FALSE, note = FALSE, bins_no = TRUE)
                gb = dt_bins[, c("expected_0", "expected_1")]
                cut_psi = dt_bins[, "PSIi"]
                bins_odds_ratio_s = dt_bins[, "odds_ratio_s"]
                gb_index = dt_bins[, "odds_ratio"]
                gb_percent = dt_bins[, "%expected"]
                effect_sz = odds_ratio = dif_gb = c()
                for (brk in 1:(dim(gb)[1] - 1)) {
                    cross_table = rbind(gb[brk,], gb[brk + 1,])
                    a = cross_table[1, 1]
                    b = cross_table[1, 2]
                    c = cross_table[2, 1]
                    d = cross_table[2, 2]
                    if (any(is.na(cross_table)) | any(cross_table < 30) | any(cross_table[, 2] / (cross_table[, 1] + cross_table[, 2]) < 0.01)) {
                        effect_sz[brk] = 0
                        odds_ratio[brk] = 1
                        dif_gb[brk] = 0
                    } else {
                        effect_sz[brk] = ((a * d - b * c) / 10000) / ((sqrt((a + b)) / 10000) * sqrt((c + d)) * sqrt((a + c)) * sqrt((b + d)))
                        odds_ratio[brk] = (a * d / 10000) / (b * c / 10000)
                        dif_gb[brk] = gb_index[brk] - gb_index[brk + 1]
                    }
                }
                gb_single = non_mono_break = gb_single_break = gb_single_pct = c()
                if (length(gb_index) > 3) {
                    for (i in 2:(length(gb_index) - 1)) {
                        gb_single[i] = (gb_index[i] > gb_index[i - 1] &
                                          gb_index[i] > gb_index[i + 1]) | (gb_index[i] < gb_index[i - 1] &
                                                                              gb_index[i] < gb_index[i + 1])
                   
                    }

                    gb_single_break = which(gb_single)                 
                    gb_single_pct = length(gb_single_break) / (length(gb_index))
                    if (length(gb_single_break)>0 && gb_single_pct > mono) {
                        non_mono_break = sort(gb_single_break)
                        # non_mono_break = gb_single_break
                    }
                }
                if (length(break_points) <= 2) {
                    break
                } else {
                    if (any(abs(effect_sz) < b_chi)) {
                        min_chi_bin = which.min(abs(effect_sz))
                        break_points = break_points[-c(min_chi_bin)]
                    } else {
                        if (any(abs(odds_ratio - 1) < b_odds)) {
                            min_odds_bin = which.min(abs(odds_ratio - 1))
                            break_points = break_points[-c(min_odds_bin)]
                        } else {
                            if (any(abs(dif_gb) < b_or)) {
                                min_gb_bin = which.min(abs(dif_gb))
                                break_points = break_points[-min_gb_bin]
                            } else {
                                if (length(non_mono_break) > 0) {
                                    min_bins_gb = min(abs(dif_gb[non_mono_break]))[1]
                                    min_gb_bin = which(abs(dif_gb) == min_bins_gb)
                                    break_points = break_points[-min_gb_bin]
                                } else {
                                    if (length(break_points) > bins_num) {
                                        min_gb_bin = which.min(abs(dif_gb))
                                        break_points = break_points[-min_gb_bin]
                                    } else {
                                        if (length(gb_percent) > 2 & any(gb_percent < bins_pct)) {
                                            bins_pct_bin = which.min(gb_percent)
                                            if (length(bins_pct_bin) > 0 &&
                                                bins_pct_bin == length(gb_percent)) {
                                                break_points = break_points[-c(bins_pct_bin - 1)]
                                            } else {
                                                break_points = break_points[-bins_pct_bin]
                                            }
                                        } else {
                                            if (length(cut_psi) > 2 & any(cut_psi > b_psi)) {
                                                max_psi_bin = which.max(cut_psi)
                                                if (length(max_psi_bin) > 0 &&
                                                    max_psi_bin == length(cut_psi)) {
                                                    break_points = break_points[-c(max_psi_bin - 1)]
                                                } else {
                                                    break_points = break_points[-max_psi_bin]
                                                }
                                            } else {
                                                if (length(bins_odds_ratio_s) > 2 & any(bins_odds_ratio_s > odds_psi)) {
                                                    max_odds_psi_bin = which.max(bins_odds_ratio_s)
                                                    if (length(max_odds_psi_bin) > 0 &&
                                                        max_odds_psi_bin == length(bins_odds_ratio_s)) {
                                                        break_points = break_points[-c(max_odds_psi_bin - 1)]
                                                    } else {
                                                        break_points = break_points[-max_odds_psi_bin]
                                                    }
                                                } else {
                                                    break
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            dt_iv = get_iv(dat = dat_test, x = x, target = target, pos_flag = pos_flag,
                           breaks = break_points, note = FALSE)
            iv_list_sub[[k]] = unlist(dt_iv[, "IV"])
            if (is.null(dt_bins)) {
                dt_bins = get_psi_iv(dat = dat_train, dat_test = dat_test, x = x,
                                     target = target, pos_flag = pos_flag,
                                     breaks = break_points, breaks_list = NULL,
                                     occur_time = occur_time, oot_pct = oot_pct,
                                     bins_total = FALSE, note = FALSE, bins_no = TRUE)
            }
            psi_list_sub[[k]] = sum(dt_bins[, "PSIi"])
            breaks_cv[[k]] = break_points
        }
        max_iv = unlist(iv_list_sub) / ifelse(max(unlist(iv_list_sub)) == 0, 1, max(unlist(iv_list_sub)))
        min_psi = (1 - unlist(psi_list_sub)) / (1 - min(unlist(psi_list_sub)))
        best_ind = which.max(max_iv * 0.8 + min_psi * 0.2)
        best_class = unique(breaks_cv[[best_ind]])
        if (length(best_ind) > 0) {
            best_breaks = sort(unlist(unique(c(breaks_cv[[best_ind]], Inf))))
        } else {
            best_breaks = sort(unique(unlist(unique(c(breaks_cv, Inf)))))
        }

        if (!is.null(sp_values) && x_miss && length(miss_num) > 0) {
            best_breaks = sort(unlist(unique(append(best_breaks, miss_num, 0))))
        }
    } else {
        best_breaks = sort(unique(unlist(c(breaks, Inf))))
    }
    return(best_breaks)
}
