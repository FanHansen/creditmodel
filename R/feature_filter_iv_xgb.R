#' Select Features using IV & PSI
#'
#'
#' \code{psi_iv_filter}  is for selecting important and stable features using IV & PSI.
#' @param dat A data.frame with independent variables and target variable.
#' @param dat_test  A data.frame of test data. Default is NULL.
#' @param target The name of target variable.
#' @param x_list Names of independent variables.
#' @param ex_cols A list of excluded variables. Regular expressions can also be used to match variable names. Default is NULL.
#' @param pos_flag The value of positive class of target variable, default: "1".
#' @param occur_time The name of the variable that represents the time at which each observation takes place.
#' @param oot_pct  Percentage of observations retained for overtime test (especially to calculate PSI). Defualt is 0.7
#' @param breaks_list A table containing a list of splitting points for each independent variable. Default is NULL.
#' @param iv_i The minimum threshold of IV. 0 < iv_i ; 0.01 to 0.1 usually work. Default: 0.01
#' @param psi_i The maximum threshold of PSI.  0 <= psi_i <=1; 0.05 to 0.2 usually work. Default: 0.1
#' @param vars_name Logical, output a list of filtered variables or table with detailed IV and PSI value of each variable. Default is FALSE.
#' @param parallel Logical, parallel computing. Default is FALSE.
#' @param note Logical, outputs info. Default is TRUE.
#' @param save_data Logical, save results in locally specified folder. Default is TRUE
#' @param file_name The name for periodically saved results files.  Default is "Featrue_importance_IV_PSI".
#' @param dir_path The path for periodically saved results files.  Default is "./variable".
#' @param ... Other parameters.
#' @return 
#' A list with the following elements:
#' \itemize{
#'   \item \code{Feature} Selected variables.
#'   \item \code{IV} IV of variables.
#'   \item \code{ PSI} PSI of variables.
#' }
#' @seealso \code{\link{xgb_filter}}, \code{\link{gbm_filter}}, \code{\link{feature_select_wrapper}}
#' @examples
#' \dontrun{
#' psi_iv_filter(dat= UCICreditCard, target = "default.payment.next.month", 
#' occur_time = "apply_date", pos_flag = list(1), parallel = FALSE,
#' ex_cols = "ID$|date$|default.payment.next.month$")
#' }

#' @export

psi_iv_filter <- function(dat, dat_test = NULL, target, x_list = NULL, breaks_list = NULL,
pos_flag = NULL, ex_cols = NULL, occur_time = NULL, oot_pct = 0.7, psi_i = 0.1, iv_i = 0.01,
vars_name = FALSE, save_data = TRUE, note = TRUE, parallel = TRUE,
file_name = "Featrue_importance_IV_PSI", dir_path = "./variable", ...) {

    IV = equal_bins = best = NULL

    if (note) {
        cat(paste("[NOTE]", "Feature filtering by IV & PSI .\n"))
    }
    dat = checking_data(dat = dat, target = target, occur_time = occur_time, pos_flag = pos_flag)
    if (is.null(dat_test)) {
        train_test = train_test_split(dat, split_type = "OOT", prop = 0.7,
        occur_time = occur_time, seed = 46, save_data = FALSE)
        dat_train = train_test$train
        dat_test = train_test$test
    } else {
        dat_train = dat
        dat_test = checking_data(dat = dat_test, target = target, occur_time = occur_time, pos_flag = pos_flag)
    }
    x_list = get_x_list(x_list = x_list, dat_train = dat, dat_test = dat_test, ex_cols = c(target, occur_time, ex_cols))
    if (any(is.na(dat_train[x_list]))) {
        dat_train = low_variance_filter(dat = dat_train, lvp = 1, note = TRUE)
        dat_train = process_nas(dat = dat, x_list = x_list, default_miss = TRUE,
        ex_cols = c(occur_time, target, ex_cols), parallel = parallel, method = "median")
    }
    if (any(is.na(dat_test[x_list]))) {
        dat_test = low_variance_filter(dat = dat_test, lvp = 1, note = TRUE)
        dat_test = process_nas(dat = dat_test, x_list = x_list, default_miss = TRUE,
        ex_cols = c(occur_time, target, ex_cols), parallel = parallel, method = "median")
    }

    x_list = get_x_list(x_list = x_list, dat_train = dat_train, dat_test = dat_test, ex_cols = c(target, occur_time, ex_cols))
    psi_sel = iv_sel = vars_psi_sel = iv_list = psi_list = iv_psi_sel = vars_sel = NULL
    psi_list = get_psi_all(dat = dat_train, dat_test = dat_test, x_list = x_list,
                           breaks_list = breaks_list, g = 5, parallel = parallel, note = note, as_table = FALSE)
    psi_sel = psi_list[psi_list$PSI <= psi_i,][1:2]
    select_vars_psi = as.character(psi_sel[, "Feature"])

    if (length(select_vars_psi) <= 1) {
        psi_sel = psi_list[psi_list$PSI <= 0.3,][1:2]
        select_vars_psi = as.character(psi_sel[, "Feature"])
    }

    if (!is.null(target) && is.element(target, colnames(dat_train))) {
        iv_list_train = get_iv_all(dat = dat_train, target = target, x_list = x_list,
        parallel = parallel, breaks_list = breaks_list, pos_flag = pos_flag, g = 50,  note = note)
        iv_sel = subset(iv_list_train, iv_list_train$IV > iv_i & iv_list_train$IV < 2)
        if (any(iv_list_train$IV > 2)) {
            cat(paste(paste(iv_list_train[which(iv_list_train$IV > 2), "Feature"], collapse = ","),  "IV  is too high to be doubted.\n"))
        }
        select_vars_iv = as.character(iv_sel[, "Feature"])
        if (length(select_vars_iv) <= 1) {
            iv_sel = subset(iv_list_train, iv_list_train$IV > 0)
            select_vars_iv = as.character(iv_sel[, "Feature"])
        }
    }
    if (length(select_vars_psi) > 0 & length(select_vars_iv) > 0) {
        iv_psi_sel = merge(iv_sel[1:2], psi_sel[1:2], by = "Feature")
        iv_psi_sel = iv_psi_sel[order(iv_psi_sel$IV, decreasing = TRUE),]
        vars_sel <- as.character(iv_psi_sel[, "Feature"])
    } else {
        if (length(select_vars_iv) > 0) {
            iv_psi_sel = iv_sel[1:2]
            iv_psi_sel = iv_psi_sel[order(iv_psi_sel$IV, decreasing = TRUE),]
            vars_sel = as.character(iv_psi_sel[, "Feature"])
        }
    }
    if (length(vars_sel) < 0) {
        vars_sel = x_list
      warning("No feature satisfies the criteria for IV & PSI feature selection, use the previous x_list.\n")
    }
    if (save_data) {
        if (!dir.exists(dir_path)) dir.create(dir_path)
        save_dt(iv_psi_sel, file_name = file_name, dir_path = dir_path, note = note)
        save_dt(psi_list, file_name = "PSI_list", dir_path = dir_path, note = note)
        save_dt(iv_list_train, file_name = "IV_list", dir_path = dir_path, note = note)
    }
    if (vars_name) {
        return(c(vars_sel))
    } else {
        return(iv_psi_sel)
    }
}


#' Select Features using XGB
#'
#'
#' \code{xgb_filter} is for selecting important features using xgboost.
#' @param dat_train A data.frame with independent variables and target variable.
#' @param dat_test  A data.frame of test data. Default is NULL.
#' @param target The name of target variable.
#' @param occur_time The name of the variable that represents the time at which each observation takes place.
#' @param x_list Names of independent variables.
#' @param ex_cols A list of excluded variables. Regular expressions can also be used to match variable names. Default is NULL.
#' @param pos_flag The value of positive class of target variable, default: "1".
#' @param xgb_params Parameters of xgboost.The complete list of parameters is available at: \url{ http://xgboost.readthedocs.io/en/latest/parameter.html}. 
#' @param cv_folds Number of cross-validations. Default: 5.
#' @param cp Threshold of XGB feature's Gain. Default is 1/number of independent variables.
#' @param seed  Random number seed. Default is 46.
#' @param vars_name Logical, output a list of filtered variables or table with detailed IV and PSI value of each variable. Default is FALSE.
#' @param note Logical, outputs info. Default is TRUE.
#' @param save_data Logical, save results results in locally specified folder. Default is TRUE
#' @param file_name The name for periodically saved results files.  Default is "Featrue_importance_XGB".
#' @param dir_path The path for periodically saved results files.  Default is "./variable".
#' @param ... Other parameters to pass to xgb_params.
#' @return Selected variables.
#' @seealso \code{\link{psi_iv_filter}}, \code{\link{gbm_filter}}, \code{\link{feature_select_wrapper}}
#' @examples
#' \dontrun{
#' xgb_filter(dat_train = UCICreditCard, dat_test = NULL, 
#' target = "default.payment.next.month", occur_time = "apply_date", 
#' ex_cols = "ID$|date$|default.payment.next.month$", vars_name = FALSE)
#' }
#' @importFrom xgboost xgb.importance xgb.train xgb.DMatrix
#' @export

xgb_filter <- function(dat_train, dat_test = NULL, target = "flag", pos_flag = NULL, x_list = NULL, occur_time = NULL, ex_cols = NULL,
xgb_params = list(nrounds = 1000, max.depth = 6, eta = 0.1, min_child_weight = 1, subsample = 1,
colsample_bytree = 1, gamma = 0, max_delta_step = 0, early_stopping_rounds = 100, eval_metric = "auc", objective = "binary:logistic"),
cv_folds = 3, cp = NULL, seed = 46, vars_name = TRUE, note = TRUE, save_data = FALSE,
file_name = "Featrue_importance_XGB", dir_path = "./variable", ...) {

    if (note) {
        cat(paste("[NOTE]", "Feature filtering by xgboost.\n"))
    }
    if (!dir.exists(dir_path)) dir.create(dir_path)
    #get parameters
    nrounds = ifelse(!is.null(xgb_params[["nrounds"]]), xgb_params[["nrounds"]], 2000)
    max.depth = ifelse(!is.null(xgb_params[["max.depth"]]), xgb_params[["max.depth"]], 6)
    eta = ifelse(!is.null(xgb_params[["eta"]]), xgb_params[["eta"]], 0.1)
    min_child_weight = ifelse(!is.null(xgb_params[["min_child_weight"]]), xgb_params[["min_child_weight"]], 1)
    subsample = ifelse(!is.null(xgb_params[["subsample"]]), xgb_params[["subsample"]], 1)
    colsample_bytree = ifelse(!is.null(xgb_params[["colsample_bytree"]]), xgb_params[["colsample_bytree"]], 1)
    gamma = ifelse(!is.null(xgb_params[["gamma"]]), xgb_params[["gamma"]], 0)
    max_delta_step = ifelse(!is.null(xgb_params[["max_delta_step"]]), xgb_params[["max_delta_step"]], 0)
    early_stopping_rounds = ifelse(!is.null(xgb_params[["early_stopping_rounds"]]), xgb_params[["early_stopping_rounds"]], 100)
    eval_metric = ifelse(!is.null(xgb_params[["eval_metric"]]), xgb_params[["eval_metric"]], "auc")
    objective = ifelse(!is.null(xgb_params[["objective"]]), xgb_params[["objective"]], "binary:logistic")
    cp = ifelse(is.null(cp), 0, cp)


    dat_train = checking_data(dat = dat_train, target = target, pos_flag = pos_flag)
    if (!is.null(dat_test)) {
        dat_test = checking_data(dat = dat_test, target = target, pos_flag = pos_flag)
    } else {
        train_test = train_test_split(dat_train, split_type = "OOT", prop = 0.7, occur_time = occur_time, seed = 46, save_data = FALSE)
        dat_train = train_test$train
        dat_test = train_test$test
    }
    x_list = get_x_list(x_list = x_list, dat_train = dat_train, dat_test = dat_test, ex_cols = c(target, occur_time, ex_cols))

    com_list = unique(c(target, occur_time, x_list))

    dat_train = dat_train[, com_list]
    dat_test = dat_test[, com_list]
    dat_ts = rbind(dat_train, dat_test)
    dat_ts = low_variance_filter(dat = dat_ts, lvp = 1, note = FALSE)
    char_x_list = get_names(dat = dat_ts, types = c('character', 'factor'),
                            ex_cols = c(target, occur_time, ex_cols), get_ex = FALSE)
    if (length(char_x_list) > 0) {
        dat_ts = one_hot_encoding(dat = dat_ts, cat_vars = char_x_list, na_act = FALSE, note = FALSE)
    }

    xg_list = get_names(dat = dat_ts, types = c('numeric', 'integer', 'double'),
                        ex_cols = c(target, occur_time, ex_cols), get_ex = FALSE)

    if (!is.null(cv_folds) && cv_folds > 1) {
        cv_list = cv_split(dat_ts, k = cv_folds, occur_time = occur_time, seed = 46)
        k = cv_folds
        if(note)cat(paste("[NOTE]", k, "folds validation to select best variables.\n"))
    } else {
        nr = nrow(dat_train)
        train_test = train_test_split(dat_ts, split_type = "byRow", prop = nr / nrow(dat_ts),
        occur_time = occur_time, seed = 46, save_data = FALSE)
        dat_train = train_test$train
        dat_test = train_test$test
        k = 1
    }
    dt_imp_XGB = best_score= list()
    for (i in 1:k) {
        if (k > 1) {
            train_sub = dat_ts[-cv_list[[i]],]
            test_sub = dat_ts[cv_list[[i]],]
        } else {
            train_sub = dat_train
            test_sub = dat_test
        }
            x_train = as.matrix(train_sub[, xg_list])
            y_train = as.numeric(as.character(train_sub[, target]))
            xgb_train = list(data = x_train, label = y_train)
            dtrain = xgb.DMatrix(data = xgb_train$data, label = xgb_train$label)

            x_test = as.matrix(test_sub[, xg_list])
            y_test = as.numeric(as.character(test_sub[, target]))
            xgb_test = list(data = x_test, label = y_test)
            dtest <- xgb.DMatrix(data = xgb_test$data, label = xgb_test$label)
            watchlist <- list(train = dtrain, eval = dtest)
        # Train a model
        if (!is.null(seed)) set.seed(seed) else set.seed(46)
        xgb_model_new = xgb.train(data = dtrain,
                                  watchlist = watchlist,
                             nrounds = nrounds,
                             max.depth = max.depth,
                             eta = eta,
                             min_child_weight = min_child_weight,
                             subsample = subsample,
                             colsample_bytree = colsample_bytree,
                             gamma = gamma,
                             max_delta_step = max_delta_step,
                             early_stopping_rounds = early_stopping_rounds,
                             eval_metric = eval_metric,
                             objective = objective,
                                  verbose = ifelse(note,1,0),
                                  maximize = TRUE)
        # feature importance

        dat_names = dimnames(x_train)[[2]]
        imp_XGB = xgb.importance(dat_names, model = xgb_model_new)
        imp_XGB = data.frame(Feature = imp_XGB[, "Feature"],
                             Importance = round(imp_XGB[, 'Gain'], 5),
                             stringsAsFactors = FALSE)
        names(imp_XGB) = c("Feature", paste("Importance.cv", i, sep = "_"))
        dt_imp_XGB[[i]] = imp_XGB
        best_score[[i]] = xgb_model_new$best_score
    }

    merge_all_xy = function(x, y) {
        merge(x, y, by.x = "Feature", by.y = "Feature", all = FALSE)
    }
    dt_imp_var <- Reduce("merge_all_xy", dt_imp_XGB)
    dt_imp_var[is.na(dt_imp_var)] = 0
    xgb_auc = unlist(best_score) - 0.5
    dt_imp_var = transform(dt_imp_var, Imp_Means_XGB = rowSums(dt_imp_var[1:k + 1] * (xgb_auc / sum(xgb_auc))))

    imp_xgb = dt_imp_var[which(dt_imp_var$Imp_Means_XGB > cp),]
    if (length(imp_xgb) <= 1) {
        imp_xgb = dt_imp_var[which(dt_imp_var$Imp_Means_XGB > 0),]
    }
    imp_xgb = imp_xgb[order(imp_xgb$Imp_Means_XGB, decreasing = TRUE),]
    imp_vars_xgb = imp_xgb[,"Feature"]
    dat_ts = de_one_hot_encoding(dat_one_hot = dat_ts[imp_vars_xgb],
                                 cat_vars = char_x_list, na_act = TRUE, note = FALSE)
    imp_vars = get_names(dat = dat_ts, types = c('character', 'factor', 'numeric', 'integer', 'double'),
                         ex_cols = c(target, occur_time, ex_cols), get_ex = FALSE)
    if (save_data) {
        save_dt(imp_vars, file_name = paste("select_vars_XGB", sep = "_"), dir_path = dir_path, note = note, as_list = TRUE)
        save_dt(dt_imp_var, file_name = file_name, dir_path = dir_path, note = note)
    }
    if (vars_name) {
        return(c(imp_vars))
    } else {
        return(imp_xgb[c("Feature", "Imp_Means_XGB")])
    }
}


#' Select Features using GBM
#'
#'
#' @description \code{gbm_filter}  is for selecting important features using GBM.
#' @param dat A data.frame with independent variables and target variable.
#' @param target The name of target variable.
#' @param x_list Names of independent variables.
#' @param ex_cols A list of excluded variables. Regular expressions can also be used to match variable names. Default is NULL.
#' @param pos_flag The value of positive class of target variable, default: "1".
#' @param gbdt_params Parameters of GBM.The complete list of parameters is available at: \code{\link{gbm}}. 
#' @param cores_num The number of CPU cores to use. 
#' @param seed  Random number seed. Default is 46.
#' @param vars_name Logical, output a list of filtered variables or table with detailed IV and PSI value of each variable. Default is TRUE.
#' @param note Logical, outputs info. Default is TRUE.
#' @param save_data Logical, save results results in locally specified folder. Default is TRUE
#' @param file_name The name for periodically saved results files.  Default is "Featrue_importance_GBDT".
#' @param dir_path The path for periodically saved results files.  Default is "./variable".
#' @param ... Other parameters to pass to gbdt_params.
#' @return Selected variables.
#' @seealso \code{\link{psi_iv_filter}}, \code{\link{xgb_filter}}, \code{\link{feature_select_wrapper}}
#' @examples
#' \dontrun{
#' gbm_filter(dat = UCICreditCard, target = "default.payment.next.month", 
#' occur_time = "apply_date", 
#' ex_cols = "ID$|date$|default.payment.next.month$", 
#' gbdt_params = list(n.trees = 100, interaction.depth = 6, shrinkage = 0.01, 
#' bag.fraction = 0.5, train.fraction = 0.7, n.minobsinnode = 30, 
#' cv.folds = 5, best_iter = TRUE, method = "cv"),
#' vars_name = FALSE)
#' }
#' @importFrom gbm gbm gbm.perf
#' @export


gbm_filter <- function(dat, target = "flag", x_list = NULL, ex_cols = NULL, pos_flag = NULL,
gbdt_params = list(n.trees = 1000, interaction.depth = 6, shrinkage = 0.01,
bag.fraction = 0.5, train.fraction = 0.7, n.minobsinnode = 30, cv.folds = 5, best_iter = TRUE, method = "cv"),
cores_num = 2, vars_name = TRUE, note = TRUE, save_data = FALSE,
file_name = "Featrue_importance_GBDT", dir_path = "./variable", seed = 46, ...) {
    dat = checking_data(dat = dat, target = target, pos_flag = pos_flag)
    if (note) {
        cat(paste("[NOTE]", "Feature filtering by gbm.\n"))
    }
    if (!dir.exists(dir_path)) dir.create(dir_path)
    #get parameters
    n.trees = ifelse(!is.null(gbdt_params[["n.trees"]]), gbdt_params[["n.trees"]], 100)
    interaction.depth = ifelse(!is.null(gbdt_params[["interaction.depth"]]), gbdt_params[["interaction.depth"]], 6)
    shrinkage = ifelse(!is.null(gbdt_params[["shrinkage"]]), gbdt_params[["shrinkage"]], 0.01)
    n.minobsinnode = ifelse(!is.null(gbdt_params[["n.minobsinnode"]]), gbdt_params[["n.minobsinnode"]], 30)
    bag.fraction = ifelse(!is.null(gbdt_params[["bag.fraction"]]), gbdt_params[["bag.fraction"]], 0.5)
    train.fraction = ifelse(!is.null(gbdt_params[["train.fraction"]]), gbdt_params[["train.fraction"]], 1)
    cv.folds = ifelse(!is.null(gbdt_params[["cv.folds"]]), gbdt_params[["cv.folds"]], 5)
    best_iter = ifelse(!is.null(gbdt_params[["best_iter"]]), gbdt_params[["best_iter"]], TRUE)
    method = ifelse(!is.null(gbdt_params[["method"]]), gbdt_params[["method"]], "cv")
    x_list = get_x_list(x_list = x_list, dat_train = dat, dat_test = NULL, ex_cols = c(target, ex_cols))
    dat = dat[, c(target, x_list)]
    char_x_list = get_names(dat = dat, types = c('character', 'factor'), ex_cols = c(target, ex_cols), get_ex = FALSE)
    if (length(char_x_list) > 0) {
        dat = one_hot_encoding(dat = dat, cat_vars = char_x_list, na_act = FALSE, note = FALSE)
    }
    gbm_list = get_names(dat = dat, types = c('numeric', 'integer', 'double'), ex_cols = c(target, ex_cols), get_ex = FALSE)
    Formula = as.formula(paste(target, paste(gbm_list, collapse = ' + '), sep = ' ~ '))
    if (!is.null(seed)) set.seed(seed) else set.seed(46)
    gbdt_model_new = gbm(
                Formula,
                data = dat, # include variables for use only.
                distribution = "bernoulli", # a loss function
                n.trees = n.trees, # the number of iterations
                shrinkage = shrinkage, # shrinkage or learning rate, It is important to know that smaller values of shrinkage (almost) always give improved predictive performance. 0.001 to 0.1 usually work
                interaction.depth = interaction.depth, # the depth of each tree
                bag.fraction = bag.fraction, # subsampling rate, (0.5 is recommended)
                train.fraction = train.fraction, # fraction of data for training,  first train.fraction*N used for training
                n.minobsinnode = n.minobsinnode, # minimum number of obs in terminal node
                cv.folds = cv.folds, # do 3-fold cross-validation
                class.stratify.cv = TRUE,
                keep.data = FALSE, # keep a copy of the dataset with the object
                verbose = TRUE, # don't print out progress
                n.cores = cores_num
                )

    best.iter = gbm.perf(gbdt_model_new, method = "cv", plot.it = TRUE, oobag.curve = FALSE)
    dt_gbm = as.data.frame(summary(gbdt_model_new, best.iter))
    dt_imp_gbm = data.frame(Feature = dt_gbm[, "var"], Imp_GBM = round(dt_gbm[, 'rel.inf'], 5), stringsAsFactors = FALSE)
    imp_gbm = subset(dt_imp_gbm, dt_imp_gbm$Imp_GBM > 0)
    imp_gbm_vars = imp_gbm[, "Feature"]
    dat = de_one_hot_encoding(dat[imp_gbm_vars], cat_vars = char_x_list, na_act = TRUE, note = FALSE)
    imp_vars = get_names(dat = dat, types = c('character', 'factor', 'numeric', 'integer', 'double'),
                         ex_cols = c(target, ex_cols), get_ex = FALSE)
    if (save_data) {
        save_dt(imp_vars, file_name = paste("select_vars_XGB", sep = "_"), dir_path = dir_path, note = FALSE, as_list = TRUE)
        save_dt(dt_imp_gbm, file_name = file_name, dir_path = dir_path, note = FALSE)
    }
    if (vars_name) {
        return(c(imp_vars))
    } else {
        return(imp_gbm)
    }
}

#' Feature Selection Wrapper
#' @description \code{feature_select_wrapper} This function uses four different methods (IV, PSI, correlation, xgboost) in order to select important features.The correlation algorithm must be used with IV.
#' @param dat_train A data.frame with independent variables and target variable.
#' @param dat_test  A data.frame of test data. Default is NULL.
#' @param target The name of target variable.
#' @param x_list Names of independent variables.
#' @param occur_time The name of the variable that represents the time at which each observation takes place.
#' @param ex_cols A list of excluded variables. Regular expressions can also be used to match variable names. Default is NULL.
#' @param pos_flag The value of positive class of target variable, default: "1".
#' @param filter The methods for selecting important and stable variables. 
#' @param breaks_list A table containing a list of splitting points for each independent variable. Default is NULL.
#' @param cv_folds Number of cross-validations. Default: 5.
#' @param hopper Logical.Filtering screening. Default is FALSE.
#' @param vars_name Logical, output a list of filtered variables or table with detailed IV and PSI value of each variable. Default is FALSE.
#' @param iv_cp The minimum threshold of IV. 0 < iv_i ; 0.01 to 0.1 usually work. Default: 0.02
#' @param psi_cp The maximum threshold of PSI.  0 <= psi_i <=1; 0.05 to 0.2 usually work. Default: 0.1
#' @param cor_cp Threshold of correlation between features. 0 <= cor_cp <=1; 0.7 to 0.98 usually work. Default is 0.98.
#' @param xgb_cp Threshold of XGB feature's Gain. 0 <= xgb_cp <=1. Default is 1/number of independent variables.
#' @param seed  Random number seed. Default is 46.
#' @param parallel Logical, parallel computing. Default is FALSE.
#' @param note Logical.Outputs info.Default is TRUE.
#' @param save_data Logical, save results in locally specified folder. Default is TRUE.
#' @param file_name The name for periodically saved results files. Default is "select_vars".
#' @param dir_path  The path for periodically saved results files. Default is "./variable"
#' @param ... Other parameters.
#' @return  A list of selected features
#' @seealso \code{\link{psi_iv_filter}}, \code{\link{xgb_filter}}, \code{\link{gbm_filter}}
#' @examples
#' \dontrun{
#' feature_select_wrapper(dat_train = UCICreditCard,
#' dat_test = NULL, target = "default.payment.next.month", 
#' occur_time = "apply_date", filter = c("IV", "PSI", "COR", "XGB"), 
#' cv_folds = 1, iv_cp = 0.01, psi_cp = 0.1, xgb_cp = 0, cor_cp = 0.98,
#' ex_cols = "ID$|date$|default.payment.next.month$",vars_name = FALSE)
#' }
#' @importFrom xgboost xgb.importance xgb.train xgb.DMatrix
#' @importFrom dplyr %>% group_by summarise
#' @export



feature_select_wrapper = function(dat_train, dat_test = NULL, x_list = NULL, target = NULL,
pos_flag = NULL, occur_time = NULL, ex_cols = NULL,
filter = c("IV", "PSI", "XGB"), cv_folds = 1, iv_cp = 0.01, psi_cp = 0.1, xgb_cp = 0, cor_cp = 0.98,
breaks_list = NULL, hopper = FALSE, vars_name = TRUE, parallel = FALSE, note = TRUE,
seed = 46, save_data = FALSE, file_name=  "select_vars", dir_path = "./variable", ...) {
    if (!is.null(filter) && any(is.element(filter, c("IV", "PSI", "XGB")))) {
        filter = filter
    } else {
        filter = c("IV", "PSI", "XGB")
    }
    cv_folds = ifelse(!is.null(cv_folds)&& is.numeric(cv_folds), cv_folds, 1)
    iv_cp = ifelse(!is.null(iv_cp) && is.numeric(iv_cp), iv_cp, 0.01)
    psi_cp = ifelse(!is.null(psi_cp) && is.numeric(psi_cp), psi_cp, 0.1)
    xgb_cp = ifelse(!is.null(xgb_cp) && is.numeric(xgb_cp), xgb_cp, 0)
   cor_cp = ifelse(!is.null(cor_cp) && is.numeric(cor_cp), cor_cp, 0.98)
    cat(paste("[NOTE] Feature filtering by", paste(filter, collapse = " & "), ".\n"))
    dat_train = checking_data(dat = dat_train, target = target, pos_flag = pos_flag)
    if (!is.null(dat_test)) {
        dat_test = checking_data(dat = dat_test, target = target, pos_flag = pos_flag)
        x_list = get_x_list(x_list = x_list, dat_train = dat_train, dat_test = dat_test,
        ex_cols = c(target, occur_time, ex_cols))
        com_list = unique(c(target, occur_time, x_list))
        dat_train = dat_train[, com_list]
        dat_test = dat_test[, com_list]
        dat_ts = rbind(dat_train, dat_test)
        dat_ts = low_variance_filter(dat = dat_ts, lvp = 1, note = FALSE)
        char_x_list = get_names(dat = dat_ts, types = c('character', 'factor'),
        ex_cols = c(target, occur_time, ex_cols), get_ex = FALSE)
        if (length(char_x_list) > 0) {
            dat_ts = one_hot_encoding(dat = dat_ts, cat_vars = char_x_list, na_act = FALSE, note = FALSE)
        }
        nr = nrow(dat_train)
        train_test = train_test_split(dat_ts, split_type = "byRow", prop = nr / nrow(dat_ts),
          occur_time = occur_time, seed = 46, save_data = FALSE,note = FALSE)
        dat_train = train_test$train
        dat_test = train_test$test
    } else {
        dat_train = low_variance_filter(dat = dat_train, lvp = 1, note = FALSE)
        char_x_list = get_names(dat = dat_train, types = c('character', 'factor'),
        ex_cols = c(target, occur_time, ex_cols), get_ex = FALSE)
        if (length(char_x_list) > 0) {
            dat_train = one_hot_encoding(dat = dat_train, cat_vars = char_x_list, na_act = FALSE, note = FALSE)
        }
        train_test = train_test_split(dat = dat_train, split_type = "OOT", prop = 0.7,
        occur_time = occur_time, seed = 46, save_data = FALSE)
        dat_train = train_test$train
        dat_test = train_test$test
    }

    imp_list = get_names(dat = dat_ts, types = c('numeric', 'integer', 'double'),
                        ex_cols = c(target, occur_time, ex_cols), get_ex = FALSE)
    dt_imp = PSI = IV = Feature = Imp_Means_XGB = NULL
    select_vars_iv = select_vars_psi = select_vars_xgb = select_vars_cor = imp_list
    if (any(filter == "PSI")) {
        psi_list_train <- get_psi_all(dat = dat_train, dat_test = dat_test, x_list = imp_list,
        breaks_list = breaks_list, ex_cols = ex_cols, g = 5, parallel = parallel, note = FALSE, as_table = FALSE)
        psi_list_t = subset(psi_list_train, PSI <= psi_cp)
        select_vars_psi = as.character(psi_list_t[, "Feature"])
        if (length(select_vars_psi) <= 1) {
            select_vars_psi = imp_list
            warning("There is no variable that meets the threshold of PSI selection.\n")
        }
        psi_list_t$Feature = gsub("\\.\\S{1,100}\\.", "", psi_list_t$Feature)
        psi_list_t = psi_list_t %>% dplyr::group_by(Feature) %>% dplyr::summarise(PSI = sum(PSI))
        if (save_data) {
            save_dt(psi_list_t, as_list = FALSE, row_names = FALSE, file_name = "Feature_Filter_PSI", dir_path = dir_path)
        }
        if (hopper) {
            imp_list = select_vars_psi
        }
      
        dt_imp = psi_list_t

    }
    if (any(filter == "IV") & !is.null(target) && is.element(target, colnames(dat_train)) ) {
        iv_list = get_iv_all(dat = dat_train, target = target, x_list = imp_list, ex_cols = ex_cols,
        pos_flag = pos_flag, equal_bins = TRUE, best = FALSE, breaks_list = breaks_list, g = 50, note = FALSE, parallel = parallel)
        iv_list_t = subset(iv_list, IV > iv_cp & IV <= 2)[, c("Feature", "IV")]
        if (any(iv_list$IV > 2)) {
            cat(paste(paste(iv_list[which(iv_list$IV > 2), "Feature"], collapse = ","), ": IV  is too high to be doubted.\n"))
        }
        select_vars_iv = as.character(iv_list_t[, "Feature"])
        if (length(select_vars_iv) <= 1) {
            select_vars_iv = imp_list
            warning("There is no variable that meets the threshold of IV selection.\n")
        }
        if (hopper) {
            imp_list = select_vars_iv
        }
        if (any(filter == "COR")) {
            select_vars_cor = fast_high_cor_filter(dat = dat_train, x_list = imp_list, com_list = iv_list,
           ex_cols = ex_cols, cor_class = FALSE, p = cor_cp, save_data = save_data,note = note,
            file_name = "Feature_filter_COR", dir_path = dir_path)
            select_vars_cor = as.character(select_vars_cor)
            if (length(select_vars_cor) <= 1) {
                select_vars_cor = x_list
                warning("There is no variable that meets the threshold of COR selection.\n")
            }
            if (hopper) {
                imp_list = select_vars_cor
            }
        }
        iv_list_t$Feature = gsub("\\.\\S{1,100}\\.", "", iv_list_t$Feature)
        iv_list_t = iv_list_t %>% dplyr::group_by(Feature) %>% dplyr::summarise(IV = sum(IV))
        if (save_data) {
            save_dt(iv_list_t, as_list = FALSE, row_names = FALSE, file_name = "Feature_Filter_IV", dir_path = dir_path)
        }

        if (is.null(dt_imp)) {
            dt_imp = iv_list_t
        } else {
            dt_imp = merge(iv_list_t, dt_imp)
        }

    }
    if (any(filter == "XGB") & !is.null(target) && is.element(target, colnames(dat_train))) {
        xgb_list = xgb_filter(dat_train = dat_train, dat_test = dat_test, target = target, x_list = imp_list,
            occur_time = occur_time, ex_cols = ex_cols, xgb_params = list(nrounds = 1000,
        max.depth = 6, eta = 0.1, min_child_weight = 1, subsample = 1, colsample_bytree = 1,
        gamma = 0, max_delta_step = 0, early_stopping_rounds = 50, eval_metric = "auc", objective = "binary:logistic"),
        cv_folds = cv_folds, cp = xgb_cp, seed = seed, note = note, vars_name = FALSE,
        save_data = FALSE)
        select_vars_xgb = xgb_list[,"Feature"]
        if (length(select_vars_xgb) <= 1) {
            select_vars_xgb = imp_list
            warning("There is no variable that meets the threshold of XGB selection.\n")
        }
        if (hopper) {
            imp_list = select_vars_xgb
        }
        xgb_list$Feature = gsub("\\.\\S{1,100}\\.", "", xgb_list$Feature)
        xgb_list = xgb_list %>% dplyr::group_by(Feature) %>% dplyr::summarise(Imp_Means_XGB = sum(Imp_Means_XGB))
        if (save_data) {
            save_dt(xgb_list, as_list = FALSE, row_names = FALSE, file_name = "Feature_Filter_XGB", dir_path = dir_path)
        }

        if (is.null(dt_imp)) {
            dt_imp = xgb_list
        } else {
            dt_imp = merge(xgb_list, dt_imp)
        }
    }

    imp_select_vars = unique(Reduce("intersect", list(select_vars_iv, select_vars_psi, select_vars_cor, select_vars_xgb)))
 
    if (length(imp_select_vars) < 0) {
        imp_select_vars = unique(Reduce("union", list(select_vars_iv, select_vars_psi, select_vars_cor, select_vars_xgb)))
        dat = de_one_hot_encoding(dat_train[imp_select_vars], cat_vars = char_x_list, na_act = TRUE, note = FALSE)
        imp_vars = get_names(dat = dat, types = c('character', 'factor', 'numeric', 'integer', 'double'),
                         ex_cols = c(target, ex_cols), get_ex = FALSE)
    } else {
        dt_imp = dt_imp[order(dt_imp$IV, decreasing = TRUE),]
        imp_vars = as.character( dt_imp[ ,"Feature"])
    }
  

    if (save_data) {
        save_dt(imp_vars, as_list = TRUE, row_names = FALSE, file_name = file_name, dir_path = dir_path)
        save_dt(dt_imp, as_list = FALSE, row_names = FALSE, file_name = paste0(file_name,"dt"), dir_path = dir_path)
    }
    if (vars_name) {
        return(imp_vars)
    } else {
        return(dt_imp)
    }
}
