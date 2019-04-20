#' Calculate IV & PSI
#'
#'
#' \code{get_iv_psi}  is used to calculate Information Value (IV)  and Population Stability Index (PSI) of an independent variable.
#' \code{get_iv_psi_all} can loop through IV & PSI for all specified independent variables.
#' @param dat A data.frame with independent variables and target variable.
#' @param dat_test  A data.frame of test data. Default is NULL.
#' @param target The name of target variable.
#' @param x_list Names of independent variables.
#' @param x  The name of an independent variable.
#' @param ex_cols A list of excluded variables. Regular expressions can also be used to match variable names. Default is NULL.
#' @param pos_flag The value of positive class of target variable, default: "1".
#' @param breaks_list A table containing a list of splitting points for each independent variable. Default is NULL.
#' @param breaks Splitting points for an independent variable. Default is NULL.
#' @param bins_total Logical, total sum for each variable.
#' @param occur_time The name of the variable that represents the time at which each observation takes place.
#' @param oot_pct  Percentage of observations retained for overtime test (especially to calculate PSI). Defualt is 0.7
#' @param best  Logical, merge initial breaks to get optimal breaks for binning.
#' @param equal_bins  Logical, generates initial breaks for equal frequency binning.
#' @param g  Number of initial breakpoints for equal frequency binning.
#' @param tree_control  Parameters of using Decision Tree to segment initial breaks. See detials: \code{\link{get_tree_breaks}}
#' @param bins_control  Parameters  used to control binning.  See detials: \code{\link{select_best_class}}, \code{\link{select_best_breaks}}
#' @param as_table Logical, output results in a table. Default is TRUE.
#' @param bins_no Logical, add serial numbers to bins. Default is FALSE.
#' @param parallel Logical, parallel computing. Default is FALSE.
#' @param note   Logical, outputs info. Default is TRUE.
#' @seealso \code{\link{get_iv}},\code{\link{get_iv_all}},\code{\link{get_psi}},\code{\link{get_psi_all}}
#' @examples
#' \dontrun{
#' iv_list = get_psi_iv_all(dat = UCICreditCard,
#' x_list = names(UCICreditCard)[3:10],
#' target = "default.payment.next.month", ex_cols = "ID|apply_date")
#' get_psi_iv(UCICreditCard, x = "PAY_3",
#' target = "default.payment.next.month",bins_total =FALSE)
#' }
#' @importFrom data.table fwrite melt fread dcast
#' @export



get_psi_iv_all <- function(dat, dat_test = NULL, x_list = NULL, target, ex_cols = NULL, pos_flag = NULL,
breaks_list = NULL, occur_time = NULL, oot_pct = 0.7, equal_bins = FALSE, tree_control = NULL, bins_control = NULL,
bins_total = TRUE, best = TRUE, g = 10, as_table = TRUE, note = TRUE, parallel = FALSE, bins_no = FALSE) {
    dat <- checking_data(dat = dat, target = target, pos_flag = pos_flag)
    opt = options(stringsAsFactors = FALSE) # 
    if (note) {
        cat(paste("[NOTE] start calculating IV & PSI...."), "\n")
    }
    x_list = get_x_list(x_list = x_list, dat_train = dat, dat_test = dat_test, ex_cols = c(target, occur_time, ex_cols))
    psi_iv_list = loop_function(func = get_psi_iv, x_list = x_list, args = list(dat = dat, dat_test = dat_test, breaks = NULL,
    breaks_list = breaks_list, target = target, pos_flag = pos_flag, best = best, equal_bins = equal_bins, tree_control = tree_control,
    oot_pct = oot_pct, occur_time = occur_time, bins_control = bins_control, g = g, note = note, bins_total = bins_total,
    as_table = as_table, bins_no = bins_no), bind = "rbind", parallel = parallel)
    options(opt) # reset
    return(psi_iv_list)
}

#' @rdname get_psi_iv_all
#' @export

get_psi_iv <- function(dat, dat_test = NULL, x, target, pos_flag = NULL, breaks = NULL, breaks_list = NULL,
occur_time = NULL, oot_pct = 0.7, equal_bins = FALSE, tree_control = NULL, bins_control = NULL,
bins_total = TRUE, best = TRUE, g = 10, as_table = TRUE, note = TRUE,bins_no = FALSE) {
    if (is.null(target)) {
        stop("target is missing!")
    }
    if (is.null(breaks)) {
        if (!is.null(breaks_list)) {
            breaks = breaks_list[which(as.character(breaks_list[,"Feature"]) == names(dat[x])), "cuts"]
        }
        if (is.null(breaks)) {
            breaks = get_breaks(dat = dat, x = x, target = target, pos_flag = pos_flag, bins_control = bins_control,
            occur_time = occur_time, oot_pct = oot_pct, equal_bins = equal_bins,
            best = best, tree_control = tree_control, g = g, note = FALSE)
        }
    }

    if (is.null(dat_test)) {
        df_bins <- split_bins(dat = dat, x = x, breaks = breaks, bins_no = bins_no)
        df_bins = as.data.frame(cbind(dat[occur_time], bins = df_bins, target = dat[, target]))
        train_test = train_test_split(dat = df_bins, prop = oot_pct, split_type = "OOT",
        occur_time = occur_time, save_data = FALSE, note = FALSE)
        dfe <- train_test$train
        dfa <- train_test$test
        dfa$ae = "actual"
        dfe$ae = "expected"
        dfe_unique = unique(dfe[, "bins"])
        dfa_unique = unique(dfa[, "bins"])
        if (any(dfa_unique != dfe_unique, na.rm = TRUE)) {
            dfa = subset(dfa, dfa$bins %in% dfe_unique)
        }
        if (any(dfe_unique != dfa_unique, na.rm = TRUE)) {
            dfe = subset(dfe, dfe$bins %in% dfa_unique)
        }
        df_ae = rbind(dfa, dfe)
    } else {
        dfe_bins <- split_bins(dat = dat, x = x, breaks = breaks,bins_no = bins_no)
        dfe = as.data.frame(cbind(bins = dfe_bins, target = dat[, target]))
        dfa_bins <- split_bins(dat = dat_test, x = x, breaks = breaks, bins_no = bins_no )
        dfa = as.data.frame(cbind(bins = dfa_bins, target = dat_test[, target]))
        dfe_unique = unique(dfe[, "bins"])
        dfa_unique = unique(dfa[, "bins"])

        if (any(dfa_unique != dfe_unique, na.rm = TRUE)) {
            dfa = subset(dfa, dfa$bins %in% dfe_unique)
        }
        if (any(dfe_unique != dfa_unique, na.rm = TRUE)) {
            dfe = subset(dfe, dfe$bins %in% dfa_unique)
        }
        dfa$ae = "actual"
        dfe$ae = "expected"
        df_ae = rbind(dfa, dfe)
    }

    bins_psi_iv <- data.table::dcast(df_ae, bins ~ ae + target, fun.aggregate = length, value.var = "ae")

    bins_psi_iv$actual_0[which(bins_psi_iv$actual_0 == 0 | is.na(bins_psi_iv$actual_0))] = 1
    bins_psi_iv$actual_1[which(bins_psi_iv$actual_1 == 0 | is.na(bins_psi_iv$actual_1))] = 1
    bins_psi_iv$expected_0[which(bins_psi_iv$expected_0 == 0 | is.na(bins_psi_iv$expected_0))] = 1
    bins_psi_iv$expected_1[which(bins_psi_iv$expected_1 == 0 | is.na(bins_psi_iv$expected_1))] = 1

    bins_psi_iv <- within(bins_psi_iv, {
       cuts = breaks[1:nrow(bins_psi_iv)]
        actual_0 = as.numeric(actual_0)
        expected_0 = as.numeric(expected_0)
        actual_1 = as.numeric(actual_1)
        expected_1 = as.numeric(expected_1)
        Feature = names(dat[x])
        `#total` = actual_0 + expected_0 + actual_1 + expected_1
        `%total` = round((actual_0 + expected_0 + actual_1 + expected_1) / sum(`#total`), 2)
        `%totalB` = round((actual_1 + expected_1) / `#total`, 2)
        good_pct = (actual_0 + expected_0) / sum(actual_0 + expected_0)
        bad_pct = (actual_1 + expected_1) / sum(actual_1 + expected_1)
        `#actual` = actual_1 + actual_0
        `%actual` = round((actual_1 + actual_0) / sum(`#actual`), 2)
        `#expected` = expected_1 + expected_0
        `%expected` = round((expected_1 + expected_0) / sum(`#expected`), 2)
        actual_pct_1 = (actual_1) / sum(actual_1)
        expected_pct_1 = (expected_1) / sum(expected_1)
        `%expectedB` = round(expected_1 / (expected_1 + expected_0), 2)
        `%actualB` = round(actual_1 / (actual_1 + actual_0), 2)
        `G/B_index` = round((expected_0 / expected_1) / (sum(expected_0) / sum(expected_1)), 2)
        `ac_G/B_index` = round((actual_0 / actual_1) / (sum(actual_0) / sum(actual_1)), 2)
        GB_psi_i = round(((expected_0 / expected_1) / (sum(expected_0 / expected_1)) -
                         (actual_0 / actual_1) / (sum(actual_0/ actual_1))) * log(((expected_0 / expected_1) / (sum(expected_0 / expected_1))) / ((actual_0 / actual_1) / (sum(actual_0/ actual_1)))), 3)
        GB_psi = sum(GB_psi_i)
        PSIi = round((`#actual` / sum(`#actual`) - `#expected` / sum(`#expected`)) * log((`#actual` / sum(`#actual`)) / (`#expected` / sum(`#expected`))), 3)
        PSI = sum(PSIi)
        IVi = round((good_pct - bad_pct) * log(good_pct / bad_pct), 3)
        IV = sum(IVi)
    })
    if (as_table) {
        df_psi_iv = bins_psi_iv[c("Feature", "bins", "cuts", "#total", "#expected", "expected_0", "expected_1",
        "#actual", "actual_0", "actual_1", "%total", "%expected", "%actual", "%totalB",
        "%expectedB", "%actualB", "G/B_index","ac_G/B_index","GB_psi_i","GB_psi","PSIi", "PSI", "IVi", "IV")]
        if (bins_total) {
            sums = c()
            for (i in 1:length(df_psi_iv)) {
                sums[i] = paste(rep("--", max(nchar(names(df_psi_iv)[i]),
            max(sapply(df_psi_iv[i], function(j) nchar(j))),na.rm = TRUE)), collapse = " ")
            }
            df_psi_iv <- rbind(df_psi_iv, sums)
        }
    } else {
        df_psi_iv = data.frame(Feature = x, IV = as.numeric(sum(bins_psi_iv$IVi)), PSI = as.numeric(sum(bins_psi_iv$PSIi)))
    }
    if (note) {
        cat(paste(x, " IV :", as.numeric(sum(bins_psi_iv$IVi)), "PSI: ", as.numeric(sum(bins_psi_iv$PSIi)), sep = "   "), "\n")
    }
    return(df_psi_iv)
}


#' Calculate Information Value (IV)
#' \code{get_iv}  is used to calculate Information Value (IV) of an independent variable.
#' \code{get_iv_all} can loop through IV for all specified independent variables.
#' @param dat A data.frame with independent variables and target variable.
#' @param target The name of target variable.
#' @param x_list Names of independent variables.
#' @param x  The name of an independent variable.
#' @param ex_cols A list of excluded variables. Regular expressions can also be used to match variable names. Default is NULL.
#' @param pos_flag Value of positive class, Default is "1".
#' @param breaks_list A table containing a list of splitting points for each independent variable. Default is NULL.
#' @param best  Logical, merge initial breaks to get optimal breaks for binning.
#' @param equal_bins  Logical, generates initial breaks for equal frequency binning.
#' @param g  Number of initial breakpoints for equal frequency binning.
#' @param tree_control  Parameters of using Decision Tree to segment initial breaks. See detials: \code{\link{get_tree_breaks}}
#' @param bins_control  Parameters  used to control binning.  See detials: \code{\link{select_best_class}}, \code{\link{select_best_breaks}}
#' @param parallel Logical, parallel computing. Default is FALSE.
#' @param note  Logical, outputs info. Default is TRUE.
#'
#' @seealso \code{\link{get_iv}},\code{\link{get_iv_all}},\code{\link{get_psi}},\code{\link{get_psi_all}}
#' @references Information Value Statistic:Bruce Lund, Magnify Analytics Solutions, a Division of Marketing Associates, Detroit, MI(Paper AA - 14 - 2013)
#' @details
#' IV Rules of Thumb for evaluating the strength a predictor
#' Less than 0.02:unpredictive
#' 0.02 to 0.1:weak
#' 0.1 to 0.3:medium
#' 0.3 + :strong
#' @examples
#' \dontrun{
#' iv_list = get_iv_all(dat = UCICreditCard, 
#' x_list = names(UCICreditCard)[3:10],
#' target = "default.payment.next.month", 
#' ex_cols = "ID|apply_date")
#' get_iv(UCICreditCard, x = "PAY_3", 
#' target = "default.payment.next.month",
#' bins_total = FALSE)
#' }
#' @export


get_iv_all <- function(dat, x_list = NULL, ex_cols = NULL, breaks_list = NULL,
                       target = NULL, pos_flag = NULL, best = TRUE,
equal_bins = FALSE, tree_control = NULL, bins_control = NULL,
g = 10, parallel = FALSE, note = TRUE) {
    dat = checking_data(dat = dat, target = target, pos_flag = pos_flag)
    opt = options(stringsAsFactors = FALSE) # 
    if (note) {
        cat(paste("[NOTE] start calculating IV...."), "\n")
    }
    x_list = get_x_list(x_list = x_list, dat_train = dat, dat_test =NULL, ex_cols = c(target, ex_cols))
    iv_list = loop_function(func = get_iv, x_list = x_list, args = list(dat = dat, breaks = NULL,
    breaks_list = breaks_list, target = target, pos_flag = pos_flag, best = best,
    equal_bins = equal_bins, tree_control = tree_control, bins_control = bins_control,
    g = g, note = note), bind = "rbind", parallel = parallel)
    options(opt) # reset
    return(iv_list)
}


#' @param breaks Splitting points for an independent variable. Default is NULL.
#' @rdname get_iv_all
#' @export


get_iv <- function(dat, x, target = NULL, pos_flag = NULL,breaks= NULL,  breaks_list = NULL,
 best = TRUE, equal_bins = FALSE, tree_control = NULL, bins_control = NULL, g = 10, note = TRUE) {

    IV = good =  bad = NULL

    if (is.null(target)) {
        stop("target is missing!")
    }
    if (is.null(breaks)) {
        if (!is.null(breaks_list)) {
            breaks = breaks_list[which(as.character(breaks_list[, "Feature"]) == names(dat[x])), "cuts"]
        }
        if (is.null(breaks)) {
            breaks = get_breaks(dat = dat, x = x, target = target, pos_flag = pos_flag,
            bins_control = bins_control, equal_bins = equal_bins, best = best,
            tree_control = tree_control, g = g, note = FALSE)
        }
    }
    best_bins <- split_bins(dat = dat, x = x, breaks = breaks, bins_no = TRUE)
    dt_bins <- table(best_bins, dat[, target])
    rm(best_bins)
    dt_bins[which(dt_bins == 0)] = 1
    dt_bins[which(is.na(dt_bins))] = 1
    bins_df = data.frame(unclass(dt_bins))
    rm(dt_bins)
    if (all(names(bins_df) == c("X0", "X1"))) {
        names(bins_df) =  c("good", "bad")
    } else {
        if (all(names(bins_df) == c("X1", "X0"))) {
            names(bins_df) = c("bad", "good")
        } else {
            stop(paste(target, "is neither 1 nor 0./n"))
        }
    }
    #IV
    bins_df = within(bins_df, {
        `%totalG` = good / sum(good)
        `%totalB` = bad / sum(bad)
        IVi = round((`%totalG` - `%totalB`) * log(`%totalG` / `%totalB`), 3)
    })
    iv_df =  data.frame(Feature = x, IV = as.numeric(sum(bins_df$IVi)))
    rm(bins_df)
    iv_df = within(iv_df, {
        strength = "Suspicious"
        strength[IV <= 0.01] <- "Unpredictive"
        strength[IV > 0.01 & IV <= 0.02] <- "Very Weak"
        strength[IV > 0.02 & IV <= 0.05] <- "Weak"
        strength[IV > 0.05 & IV <= 0.1] <- "Medium"
        strength[IV > 0.1 & IV <= 0.3] <- "Strong"
        strength[IV <= 3 & IV > 0.3] <- "Very Strong"
    })
    if (note) {
        cat(paste(x, " IV :", iv_df$IV, iv_df$strength, sep = "   "), "\n")
    }

    iv_df
}


#' Calculate Population Stability Index (PSI)
#' \code{get_psi} is used to calculate Population Stability Index (PSI)  of an independent variable.
#' \code{get_psi_all} can loop through PSI for all specified independent variables.
#' @param dat A data.frame with independent variables and target variable.
#' @param dat_test  A data.frame of test data. Default is NULL.
#' @param x_list Names of independent variables.
#' @param x  The name of an independent variable.
#' @param ex_cols Names of excluded variables. Regular expressions can also be used to match variable names. Default is NULL.
#' @param pos_flag Value of positive class, Default is "1".
#' @param breaks_list A table containing a list of splitting points for each independent variable. Default is NULL.
#' @param breaks Splitting points for an independent variable. Default is NULL.
#' @param g  Number of initial breakpoints for equal frequency binning.
#' @param occur_time The name of the variable that represents the time at which each observation takes place.
#' @param start_date The earliest occurrence time of observations.
#' @param cut_date Time points for spliting data sets, e.g. : spliting Actual and Expected data sets.
#' @param oot_pct  Percentage of observations retained for overtime test (especially to calculate PSI). Defualt is 0.7
#' @param as_table Logical, output results in a table. Default is TRUE.
#' @param bins_no Logical, add serial numbers to bins. Default is TRUE.
#' @param parallel Logical, parallel computing. Default is FALSE.
#' @param note   Logical, outputs info. Default is TRUE.
#' @seealso \code{\link{get_iv}},\code{\link{get_iv_all}},\code{\link{get_psi}},\code{\link{get_psi_all}}
#' @details
#' PSI Rules for evaluating the stability of a predictor
#' Less than 0.02: Very stable
#' 0.02 to 0.1: Stable
#' 0.1 to 0.2: Unstable
#' 0.2 to 0.5] : Change
#' more than 0.5: Great change
#' @examples
#' \dontrun{
#' #  dat_test is null
#' get_psi(dat = UCICreditCard, x = "PAY_3", occur_time = "apply_date")
#' # dat_test is not all 
#' # train_test split
#' train_test = train_test_split(dat = UCICreditCard, prop = 0.7, split_type = "OOT",
#'                              occur_time = "apply_date", start_date = NULL, cut_date = NULL,
#'                             save_data = FALSE, note = FALSE)
#' dat_ex = train_test$train
#' dat_ac = train_test$test
#' # generate psi table
#' get_psi(dat = dat_ex, dat_test = dat_ac, x = "PAY_3",
#'        occur_time = "apply_date", bins_no = TRUE)
#' }
#' @export



get_psi_all <- function(dat, x_list = NULL, dat_test = NULL, breaks_list = NULL, occur_time = NULL,
start_date = NULL, cut_date = NULL, oot_pct = 0.7, pos_flag = NULL,
parallel = FALSE, ex_cols = NULL, as_table = FALSE, g = 10, bins_no = TRUE, note = TRUE) {
    if (note) {
        cat(paste("[NOTE] start calculating PSI...."), "\n")
    }
    opt = options(stringsAsFactors = FALSE) # 
    if (is.null(x_list)) {
        if (!is.null(breaks_list)) {
            x_list = unique(as.character(breaks_list[, "Feature"]))
        } else {
            x_list = get_x_list(x_list = x_list, dat_train = dat, dat_test = dat_test, ex_cols = c(occur_time, ex_cols))
        }
    }
    if (is.null(dat_test) && !is.null(occur_time) && any(names(dat) == occur_time)) {
        if (!is_date(dat[, occur_time])) {
            dat =  time_transfer(dat, date_cols = occur_time, note = FALSE)
        }
        if (is_date(dat[, occur_time])) {
            if (is.null(cut_date)) {
                cut_date = date_cut(dat_time = dat[, occur_time], pct = oot_pct)
            }
            if (is.null(start_date)) {
                start_date = date_cut(dat_time = dat[, occur_time], pct = 0)
            }
        } else {
            stop(paste(occur_time, "is not Date or Time"))
        }
    } 
    psi_list = loop_function(func = get_psi, x_list = x_list, args = list(dat = dat, dat_test = dat_test, 
    breaks = NULL, breaks_list = breaks_list, occur_time = occur_time, start_date = start_date, cut_date = cut_date,
    oot_pct = oot_pct, pos_flag = pos_flag, as_table = as_table, g = g, note = note, bins_no = bins_no), bind = "rbind", parallel = parallel)
    options(opt) # reset
    return(psi_list)
}


#' @rdname get_psi_all
#' @export

get_psi <- function(dat, x, dat_test = NULL, occur_time = NULL, start_date = NULL, cut_date = NULL,
 pos_flag = NULL, breaks = NULL, breaks_list = NULL, oot_pct = 0.7,  g = 10,
 as_table = TRUE, note = TRUE, bins_no = TRUE) {
    bins = PSI = actual = expected = NULL
    if (!is.null(breaks_list)) {
        breaks = breaks_list[which(as.character(breaks_list[, "Feature"]) == names(dat[x])), "cuts"]
    }
    if (is.null(breaks)) {
        breaks = get_breaks(dat = dat, x = x,  pos_flag = pos_flag, equal_bins = TRUE, best = FALSE, g = g, note = FALSE)
    }

    if (is.null(dat_test)) {
        dat$bins = split_bins(dat = dat, x = x, breaks = breaks, bins_no = bins_no)
        df_ae = train_test_split(dat = dat, prop = oot_pct, split_type = "OOT", occur_time = occur_time,
        start_date = start_date, cut_date = cut_date, save_data = FALSE, note = FALSE)
        dfe = df_ae$train
        dfa = df_ae$test
        dfe_unique = unique(dfe[, "bins"])
        dfa_unique = unique(dfa[, "bins"])


        if (any(dfa_unique != dfe_unique, na.rm = TRUE)) {
            dfa = subset(dfa, dfa$bins %in% dfe_unique)
        }
        if (any(dfe_unique != dfa_unique, na.rm = TRUE)) {
            dfe = subset(dfe, dfe$bins %in% dfa_unique)
        }
        dfa$ae = "actual"
        dfe$ae = "expected"
        df_ae = rbind(dfa, dfe)
    } else {
        dat$ae = "expected"
        dat_test$ae = "actual"
        dat$bins = split_bins(dat = dat, x = x, breaks = breaks, bins_no = bins_no)
        dat_test$bins = split_bins(dat = dat_test, x = x, breaks = breaks, bins_no = bins_no)
        dfe = dat[, c("ae", "bins")]
        dfa = dat_test[, c("ae", "bins")]
        dfe_unique = unique(dfe[, "bins"])
        dfa_unique = unique(dfa[, "bins"])

        if (any(dfa_unique != dfe_unique, na.rm = TRUE)) {
            dfa = subset(dfa, dfa$bins %in% dfe_unique)
        }
        if (any(dfe_unique != dfa_unique, na.rm = TRUE)) {
            dfe = subset(dfe, dfe$bins %in% dfa_unique)
        }
        df_ae = rbind(dfa, dfe)
    }
    df_psi = data.table::dcast(df_ae, bins ~ ae, fun.aggregate = length, value.var = "ae")

    df_psi$actual[which(df_psi$actual == 0 | is.na(df_psi$actual))] = 1
    df_psi$expected[which(df_psi$expected == 0 | is.na(df_psi$expected))] = 1

    df_psi =  within(df_psi, {
        Ac_pct = actual / sum(actual)
        Ex_pct = expected / sum(expected)
        PSI_i = round((Ac_pct - Ex_pct) * log(Ac_pct / Ex_pct), 3)
    })

    if (as_table) {
        dat_psi = data.frame(Feature = x, Bins = df_psi$bins, actual = df_psi$actual, expected = df_psi$expected,
                        Ac_pct = as_percent(df_psi$Ac_pct, digits = 3), Ex_pct = as_percent(df_psi$Ex_pct, digits = 3),
                        PSI_i = df_psi$PSI_i, PSI = as.numeric(sum(df_psi$PSI_i)))
    } else {
        dat_psi = data.frame(Feature = x, PSI = as.numeric(sum(df_psi$PSI_i)))
    }
    dt_psi <- within(dat_psi, {
        stability = "Great change"
        stability[PSI <= 0.02] <- "Very stable"
        stability[PSI > 0.02 & PSI <= 0.1] <- "Stable"
        stability[PSI > 0.1 & PSI <= 0.2] <- "Unstable"
        stability[PSI > 0.2 & PSI <= 0.5] <- "Change"
        stability[PSI > 0.5] <- "Great change"
        })
    if (note) {
        cat(paste(x, " PSI :", as.numeric(dt_psi$PSI[1]), dt_psi$stability[1], sep = "  "), "\n")
    }
    rm(df_psi,dt_psi)
    return(dat_psi)
}