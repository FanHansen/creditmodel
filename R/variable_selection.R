#'  Variable reduction based on Information Value & Population Stability Index filter
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
#' @param save_data Logical, save results in locally specified folder. Default is FALSE.
#' @param file_name The name for periodically saved results files.  Default is "Featrue_importance_IV_PSI".
#' @param dir_path The path for periodically saved results files.  Default is tempdir().
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
#' psi_iv_filter(dat= UCICreditCard[1:1000,c(2,4,8:9,26)],
#'              target = "default.payment.next.month",
#'              occur_time = "apply_date",
#'              parallel = FALSE)
#' @export

psi_iv_filter <- function(dat, dat_test = NULL, target, x_list = NULL,
                          breaks_list = NULL,pos_flag = NULL,
                          ex_cols = NULL, occur_time = NULL,
                          oot_pct = 0.7, psi_i = 0.1, iv_i = 0.01,
                          vars_name = FALSE, note = TRUE, parallel = FALSE,
                          save_data = FALSE, file_name = NULL,
                          dir_path = tempdir(), ...) {

  IV = equal_bins = best = NULL
  if (note) {
    cat(paste("[NOTE]", "Feature filtering by IV & PSI .\n"))
  }
  dat = checking_data(dat = dat, target = target, occur_time = occur_time, pos_flag = pos_flag)
  if (is.null(dat_test)) {
    train_test = train_test_split(dat, split_type = "OOT", prop = 0.7,
                                  occur_time = occur_time, seed = 46, save_data = FALSE, note = FALSE)
    dat_train = train_test$train
    dat_test = train_test$test
  } else {
    dat_train = dat
    dat_test = checking_data(dat = dat_test, target = target,
                             occur_time = occur_time, pos_flag = pos_flag)
  }
  x_list = get_x_list(x_list = x_list, dat_train = dat, dat_test = dat_test,
                      ex_cols = c(target, occur_time, ex_cols))
  if (any(is.na(dat_train[x_list]))) {
    dat_train = low_variance_filter(dat = dat_train, lvp = 1, note = FALSE)
    dat_train = process_nas(dat = dat, x_list = x_list, default_miss = TRUE,
                            ex_cols = c(occur_time, target, ex_cols), parallel = parallel, method = "median")
  }
  if (any(is.na(dat_test[x_list]))) {
    dat_test = low_variance_filter(dat = dat_test, lvp = 1, note = FALSE)
    dat_test = process_nas(dat = dat_test, x_list = x_list, default_miss = TRUE,
                           ex_cols = c(occur_time, target, ex_cols), parallel = parallel, method = "median")
  }

  x_list = get_x_list(x_list = x_list, dat_train = dat_train, dat_test = dat_test,
                      ex_cols = c(target, occur_time, ex_cols))
  psi_sel = iv_sel = vars_psi_sel = iv_list = psi_list = iv_psi_sel = vars_sel = NULL
  psi_list = get_psi_all(dat = dat_train, dat_test = dat_test,
                         x_list = x_list,
                         breaks_list = breaks_list, g = 5,
                         parallel = parallel, note = note,
                         as_table = FALSE)
  psi_sel = psi_list[psi_list$PSI <= psi_i,][1:2]
  select_vars_psi = as.character(psi_sel[, "Feature"])

  if (length(select_vars_psi) <= 1) {
    psi_sel = psi_list[psi_list$PSI <= 0.3,][1:2]
    select_vars_psi = as.character(psi_sel[, "Feature"])
  }

  if (!is.null(target) && is.element(target, colnames(dat_train))) {
    iv_list_train = get_iv_all(dat = dat_train, target = target, x_list = x_list,
                               parallel = parallel, breaks_list = breaks_list,
                               pos_flag = pos_flag, g = 50,  note = note)
    iv_sel = subset(iv_list_train, iv_list_train$IV > iv_i & iv_list_train$IV < 2)
    if (any(iv_list_train$IV > 2)) {
      cat(paste(paste(iv_list_train[which(iv_list_train$IV > 2), "Feature"], collapse = ","),
                "IV  is too high to be doubted.\n"))
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
    dir_path = ifelse(!is.character(dir_path),  tempdir(), dir_path)
    if (!dir.exists(dir_path)) dir.create(dir_path)
    if (!is.character(file_name)) { file_name = NULL }
    save_dt(iv_psi_sel, file_name = ifelse(is.null(file_name), "feature.IV_PSI", paste(file_name, "feature.COR.group", sep = ".")), dir_path = dir_path, note = note)
    save_dt(psi_list, file_name = ifelse(is.null(file_name), "feature.PSI.table", paste(file_name, "feature.PSI.table", sep = ".")), dir_path = dir_path, note = note)
    save_dt(iv_list_train, file_name = ifelse(is.null(file_name), "feature.IV.table", paste(file_name, "feature.IV.table", sep = ".")), dir_path = dir_path, note = note)
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
#' @param save_data Logical, save results results in locally specified folder. Default is FALSE.
#' @param file_name The name for periodically saved results files.  Default is "Featrue_importance_XGB".
#' @param dir_path The path for periodically saved results files.  Default is "./variable".
#' @param ... Other parameters to pass to xgb_params.
#' @return Selected variables.
#' @seealso \code{\link{psi_iv_filter}}, \code{\link{gbm_filter}}, \code{\link{feature_select_wrapper}}
#' @examples
#' mm <- xgb_filter(dat_train = UCICreditCard[1:1000,c(2,4,8:9,26)], dat_test = NULL,
#' target = "default.payment.next.month", occur_time = "apply_date",cv_folds = 1,
#' ex_cols = "ID$|date$|default.payment.next.month$", vars_name = FALSE)
#' @importFrom xgboost xgb.importance xgb.train xgb.DMatrix
#' @export

xgb_filter <- function(dat_train, dat_test = NULL, target = "flag", pos_flag = NULL,
                       x_list = NULL, occur_time = NULL, ex_cols = NULL,
                       xgb_params = list(nrounds = 1000, max.depth = 6, eta = 0.1,
                                         min_child_weight = 1, subsample = 1,
                                         colsample_bytree = 1, gamma = 0, max_delta_step = 0,
                                         early_stopping_rounds = 100,
                                         eval_metric = "auc", objective = "binary:logistic"),
                       cv_folds = 3, cp = NULL, seed = 46, vars_name = TRUE,
                       note = TRUE, save_data = FALSE,
                       file_name = NULL, dir_path = tempdir(), ...) {
  if (note) cat(paste("[NOTE]", "Feature filtering by XGB.\n"))
  #get parameters
  nrounds = ifelse(!is.null(xgb_params[["nrounds"]]),
                   xgb_params[["nrounds"]], 2000)
  max.depth = ifelse(!is.null(xgb_params[["max.depth"]]),
                     xgb_params[["max.depth"]], 6)
  eta = ifelse(!is.null(xgb_params[["eta"]]),
               xgb_params[["eta"]], 0.1)
  min_child_weight = ifelse(!is.null(xgb_params[["min_child_weight"]]),
                            xgb_params[["min_child_weight"]], 1)
  subsample = ifelse(!is.null(xgb_params[["subsample"]]),
                     xgb_params[["subsample"]], 1)
  colsample_bytree = ifelse(!is.null(xgb_params[["colsample_bytree"]]),
                            xgb_params[["colsample_bytree"]], 1)
  gamma = ifelse(!is.null(xgb_params[["gamma"]]),
                 xgb_params[["gamma"]], 0)
  max_delta_step = ifelse(!is.null(xgb_params[["max_delta_step"]]),
                          xgb_params[["max_delta_step"]], 0)
  early_stopping_rounds = ifelse(!is.null(xgb_params[["early_stopping_rounds"]]),
                                 xgb_params[["early_stopping_rounds"]], 100)
  eval_metric = ifelse(!is.null(xgb_params[["eval_metric"]]),
                       xgb_params[["eval_metric"]], "auc")
  objective = ifelse(!is.null(xgb_params[["objective"]]),
                     xgb_params[["objective"]], "binary:logistic")
  cp = ifelse(is.null(cp), 0, cp)

  dat_train = checking_data(dat = dat_train, target = target, pos_flag = pos_flag)
  if (!is.null(dat_test)) {
    dat_test = checking_data(dat = dat_test, target = target, pos_flag = pos_flag)
  } else {
    train_test = train_test_split(dat_train, split_type = "OOT", prop = 0.7,note = FALSE,
                                  occur_time = occur_time, seed = 46, save_data = FALSE)
    dat_train = train_test$train
    dat_test = train_test$test
  }
  x_list = get_x_list(x_list = x_list, dat_train = dat_train, dat_test = dat_test,
                      ex_cols = c(target, occur_time, ex_cols))

  com_list = unique(c(target, occur_time, x_list))

  dat_train = dat_train[, com_list]
  dat_test = dat_test[, com_list]
  dat_ts = rbind(dat_train, dat_test)
  dat_ts = low_variance_filter(dat = dat_ts, lvp = 1, note = FALSE)
  char_x_list = get_names(dat = dat_ts, types = c('character', 'factor'),
                          ex_cols = c(target, occur_time, ex_cols),
                          get_ex = FALSE)
  if (length(char_x_list) > 0) {
    dat_ts = one_hot_encoding(dat = dat_ts, cat_vars = char_x_list,
                              na_act = FALSE, note = FALSE)
  }

  xg_list = get_names(dat = dat_ts, types = c('numeric', 'integer', 'double'),
                      ex_cols = c(target, occur_time, ex_cols), get_ex = FALSE)

  if (!is.null(cv_folds) && cv_folds > 1) {
    cv_list = cv_split(dat_ts, k = cv_folds, occur_time = occur_time, seed = 46)
    k = cv_folds
  } else {
    nr = nrow(dat_train)
    train_test = train_test_split(dat_ts, split_type = "byRow", prop = nr / nrow(dat_ts),
                                  occur_time = occur_time, seed = 46, note = FALSE, save_data = FALSE)
    dat_train = train_test$train
    dat_test = train_test$test
    k = 1
  }
  dt_imp_XGB = bst_score= list()
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
                              verbose = 0,
                              maximize = TRUE)
    # feature importance

    dat_names = dimnames(x_train)[[2]]
    imp_XGB = xgb.importance(dat_names, model = xgb_model_new)
    imp_XGB = data.frame(Feature = imp_XGB[, "Feature"],
                         Importance = round(imp_XGB[, 'Gain'], 5),
                         stringsAsFactors = FALSE)
    names(imp_XGB) = c("Feature", paste("Importance.cv", i, sep = "_"))
    dt_imp_XGB[[i]] = imp_XGB
    bst_score[[i]] = xgb_model_new$best_score
  }

  merge_all_xy = function(x, y) {
    merge(x, y, by.x = "Feature", by.y = "Feature", all = FALSE)
  }
  dt_imp_var <- Reduce("merge_all_xy", dt_imp_XGB)
  dt_imp_var[is.na(dt_imp_var)] = 0
  xgb_auc = unlist(bst_score) - 0.5
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
    dir_path = ifelse(!is.character(dir_path), tempdir(), dir_path)
    if (!dir.exists(dir_path)) dir.create(dir_path)
    if (!is.character(file_name)) { file_name = NULL }
    save_dt(imp_vars, file_name = ifelse(is.null(file_name), "feature.XGB", paste(file_name, "feature.XGB", sep = ".")),
            dir_path = dir_path, note = note, as_list = TRUE)
    save_dt(dt_imp_var, file_name = ifelse(is.null(file_name), "feature.XGB.table", paste(file_name, "feature.XGB.table", sep = ".")),
            dir_path = dir_path, note = note)
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
#' @param GBM.params Parameters of GBM.The complete list of parameters is available at: \code{\link{gbm}}.
#' @param cores_num The number of CPU cores to use.
#' @param seed  Random number seed. Default is 46.
#' @param vars_name Logical, output a list of filtered variables or table with detailed IV and PSI value of each variable. Default is TRUE.
#' @param note Logical, outputs info. Default is TRUE.
#' @param save_data Logical, save results results in locally specified folder. Default is FALSE.
#' @param file_name The name for periodically saved results files.  Default is "Featrue_importance_GBDT".
#' @param dir_path The path for periodically saved results files.  Default is "./variable".
#' @param ... Other parameters to pass to gbdt_params.
#' @return Selected variables.
#' @seealso \code{\link{psi_iv_filter}}, \code{\link{xgb_filter}}, \code{\link{feature_select_wrapper}}
#' @examples
#' GBM.params = gbm_params(n.trees = 2, interaction.depth = 2, shrinkage = 0.1,
#'                        bag.fraction = 1, train.fraction = 1,
#'                        n.minobsinnode = 30,
#'                      cv.folds = 2)
#' \dontrun{
#'  features <- gbm_filter(dat = UCICreditCard[1:1000, c(8:12, 26)],
#'          target = "default.payment.next.month",
#'       occur_time = "apply_date",
#'      GBM.params = GBM.params
#'        , vars_name = FALSE)
#'}
#' @importFrom gbm gbm gbm.perf
#' @export


gbm_filter <- function(dat, target = NULL, x_list = NULL, ex_cols = NULL, pos_flag = NULL,
                       GBM.params = gbm_params(),
                       cores_num = 2, vars_name = TRUE, note = TRUE, save_data = FALSE,
                       file_name = NULL,
                       dir_path = tempdir(), seed = 46, ...) {
  dat = checking_data(dat = dat, target = target, pos_flag = pos_flag)
  if (note) {
    cat(paste("[NOTE]", "Feature filtering by GBM.\n"))
  }

  #get parameters
  n.trees = ifelse(!is.null(GBM.params[["n.trees"]]), GBM.params[["n.trees"]], 100)
  interaction.depth = ifelse(!is.null(GBM.params[["interaction.depth"]]), GBM.params[["interaction.depth"]], 6)
  shrinkage = ifelse(!is.null(GBM.params[["shrinkage"]]), GBM.params[["shrinkage"]], 0.1)
  n.minobsinnode = ifelse(!is.null(GBM.params[["n.minobsinnode"]]), GBM.params[["n.minobsinnode"]], 30)
  bag.fraction = ifelse(!is.null(GBM.params[["bag.fraction"]]), GBM.params[["bag.fraction"]], 0.5)
  train.fraction = ifelse(!is.null(GBM.params[["train.fraction"]]), GBM.params[["train.fraction"]], 1)
  cv.folds = ifelse(!is.null(GBM.params[["cv.folds"]]), GBM.params[["cv.folds"]], 0)
  x_list = get_x_list(x_list = x_list, dat_train = dat, dat_test = NULL, ex_cols = c(target, ex_cols))
  dat = dat[, c(target, x_list)]
  char_x_list = get_names(dat = dat, types = c('character', 'factor'),
                          ex_cols = c(target, ex_cols), get_ex = FALSE)
  if (length(char_x_list) > 0) {
    dat = one_hot_encoding(dat = dat, cat_vars = char_x_list, na_act = FALSE, note = FALSE)
  }
  gbm_list = get_names(dat = dat, types = c('numeric', 'integer', 'double'),
                       ex_cols = c(target, ex_cols), get_ex = FALSE)
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
    verbose = FALSE, # don't print out progress
    n.cores = cores_num
  )

  best.iter = gbm.perf(gbdt_model_new, method = "cv", plot.it = FALSE, oobag.curve = FALSE)
  dt_gbm = as.data.frame(summary(gbdt_model_new, best.iter, plotit = FALSE ))
  dt_imp_gbm = data.frame(Feature = dt_gbm[, "var"],
                          Imp_GBM = round(dt_gbm[, 'rel.inf'], 5),
                          stringsAsFactors = FALSE)
  imp_gbm = subset(dt_imp_gbm, dt_imp_gbm$Imp_GBM > 0)
  imp_gbm_vars = imp_gbm[, "Feature"]
  dat = de_one_hot_encoding(dat[imp_gbm_vars], cat_vars = char_x_list,
                            na_act = TRUE, note = FALSE)
  imp_vars = get_names(dat = dat,
                       types = c('character', 'factor', 'numeric', 'integer', 'double'),
                       ex_cols = c(target, ex_cols),
                       get_ex = FALSE)
  if (save_data) {
    dir_path = ifelse(!is.character(dir_path), tempdir(), dir_path)
    if (!dir.exists(dir_path)) dir.create(dir_path)
    if (!is.character(file_name)) { file_name = NULL }
    save_dt(imp_vars, file_name = ifelse(is.null(file_name), "feature.GBM",
                                         paste(file_name, "feature.GBM", sep = ".")), dir_path = dir_path, note = FALSE, as_list = TRUE)
    save_dt(dt_imp_gbm, file_name = ifelse(is.null(file_name), "feature.GBM.table",
                                           paste(file_name, "feature.GBM.table", sep = ".")), dir_path = dir_path, note = FALSE)
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
#' @param note Logical.Outputs info. Default is TRUE.
#' @param save_data Logical, save results in locally specified folder. Default is FALSE.
#' @param file_name The name for periodically saved results files. Default is "select_vars".
#' @param dir_path  The path for periodically saved results files. Default is "./variable"
#' @param ... Other parameters.
#' @return  A list of selected features
#' @seealso \code{\link{psi_iv_filter}}, \code{\link{xgb_filter}}, \code{\link{gbm_filter}}
#' @examples
#' feature_select_wrapper(dat_train = UCICreditCard[1:1000,c(8:12,26)],
#'                       dat_test = NULL, target = "default.payment.next.month",
#'                       occur_time = "apply_date", filter = c("IV", "PSI"),
#'                       cv_folds = 1, iv_cp = 0.01, psi_cp = 0.1, xgb_cp = 0, cor_cp = 0.98,
#'                       vars_name = FALSE,note = FALSE)
#' @importFrom xgboost xgb.importance xgb.train xgb.DMatrix
#' @importFrom dplyr %>% group_by summarise
#' @export




feature_select_wrapper = function(dat_train, dat_test = NULL, x_list = NULL, target = NULL,
                                  pos_flag = NULL, occur_time = NULL, ex_cols = NULL,
                                  filter = c("IV", "PSI", "XGB","COR"),
                                  cv_folds = 1,
                                  iv_cp = 0.01, psi_cp = 0.1, xgb_cp = 0, cor_cp = 0.98,
                                  breaks_list = NULL, hopper = FALSE, vars_name = TRUE,
                                  parallel = FALSE, note = TRUE, seed = 46,
                                  save_data = FALSE, file_name= NULL,
                                  dir_path = tempdir(), ...) {

  if (!is.null(filter) && any(is.element(filter, c("IV", "PSI", "XGB","COR")))) {
    filter = filter
  } else {
    filter = c("IV","COR", "PSI", "XGB")
  }
  cv_folds = ifelse(!is.null(cv_folds)&& is.numeric(cv_folds), cv_folds, 1)
  iv_cp = ifelse(!is.null(iv_cp) && is.numeric(iv_cp), iv_cp, 0.01)
  psi_cp = ifelse(!is.null(psi_cp) && is.numeric(psi_cp), psi_cp, 0.1)
  xgb_cp = ifelse(!is.null(xgb_cp) && is.numeric(xgb_cp), xgb_cp, 0)
  cor_cp = ifelse(!is.null(cor_cp) && is.numeric(cor_cp), cor_cp, 0.98)
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
      dat_ts = one_hot_encoding(dat = dat_ts, cat_vars = char_x_list,
                                na_act = FALSE, note = FALSE)
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
      dat_train = one_hot_encoding(dat = dat_train, cat_vars = char_x_list,
                                   na_act = FALSE, note = FALSE)
    }
    train_test = train_test_split(dat = dat_train, split_type = "OOT", prop = 0.7,
                                  occur_time = occur_time, seed = 46, save_data = FALSE, note = FALSE)
    dat_train = train_test$train
    dat_test = train_test$test
  }

  imp_list = get_names(dat = dat_train, types = c('numeric', 'integer', 'double'),
                       ex_cols = c(target, occur_time, ex_cols), get_ex = FALSE)
  dt_imp = PSI = IV = Feature = Imp_Means_XGB = NULL
  select_vars_iv = select_vars_psi = select_vars_xgb = select_vars_cor = imp_list
  if (save_data) {
    dir_path = ifelse(!is.character(dir_path), tempdir(), dir_path)
    if (!dir.exists(dir_path)) dir.create(dir_path)
    if (!is.character(file_name)) file_name = NULL
  }
  if (any(filter == "PSI")) {
    if (note) cat(paste("[NOTE]", "Feature filtering by PSI.\n"))
    psi_list_train <- get_psi_all(dat = dat_train, dat_test = dat_test, x_list = imp_list,
                                  breaks_list = breaks_list, ex_cols = ex_cols, g = 5,
                                  parallel = parallel, note = FALSE, as_table = FALSE)
    psi_list_t = subset(psi_list_train, PSI <= psi_cp)
    select_vars_psi = as.character(psi_list_t[, "Feature"])
    if (length(select_vars_psi) <= 1) {
      select_vars_psi = imp_list
      warning("There is no variable that meets the threshold of PSI selection.\n")
    }
    psi_list_t$Feature = gsub("\\.\\S{1,100}\\.", "", psi_list_t$Feature)
    psi_list_t = psi_list_t %>% dplyr::group_by(Feature) %>% dplyr::summarise(PSI = sum(PSI))
    if (save_data) {
      save_dt(psi_list_t, as_list = FALSE, row_names = FALSE, note = note,
              file_name = ifelse(is.null(file_name), "feature.filter.PSI", paste(file_name, "feature.filter.PSI", sep = ".")), dir_path = dir_path)
    }
    if (hopper) {
      imp_list = select_vars_psi
    }

    dt_imp = psi_list_t

  }
  if (any(filter == "IV" | filter == "COR") & !is.null(target) &&
      is.element(target, colnames(dat_train))) {
    if (note) cat(paste("[NOTE]", "Feature filtering by IV.\n"))
    iv_list = get_iv_all(dat = dat_train, target = target,
                         x_list = imp_list, ex_cols = ex_cols,
                         pos_flag = pos_flag,
                         equal_bins = TRUE, best = FALSE,
                         breaks_list = breaks_list, g = 50,
                         note = FALSE, parallel = parallel)
    iv_list_t = subset(iv_list, IV > iv_cp & IV <= 2)[, c("Feature", "IV")]
    if (any(iv_list$IV > 2)) {
      cat(paste(paste(iv_list[which(iv_list$IV > 2), "Feature"], collapse = ","),
                ": IV  is too high to be doubted.\n"))
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
      if(note) cat(paste("[NOTE]", "Feature filtering by COR.\n"))
      select_vars_cor = fast_high_cor_filter(dat = dat_train,
                                             x_list = imp_list, com_list = iv_list,
                                             ex_cols = ex_cols, cor_class = FALSE,
                                             p = cor_cp, note = FALSE,
                                             save_data = FALSE)
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
      save_dt(iv_list_t, as_list = FALSE, row_names = FALSE,
              file_name = ifelse(is.null(file_name), "feature.filter.IV", paste(file_name, "feature.filter.IV", sep = ".")), dir_path = dir_path)
    }

    if (is.null(dt_imp)) {
      dt_imp = iv_list_t
    } else {
      dt_imp = merge(iv_list_t, dt_imp)
    }

  }
  if (any(filter == "XGB") & !is.null(target) && is.element(target, colnames(dat_train))) {
    xgb_list = xgb_filter(dat_train = dat_train, dat_test = dat_test, target = target,
                          x_list = imp_list,occur_time = occur_time, ex_cols = ex_cols,
                          xgb_params = list(nrounds = 1000,max.depth = 6, eta = 0.1,
                                            min_child_weight = 1, subsample = 1,
                                            colsample_bytree = 1,gamma = 0,
                                            max_delta_step = 0, early_stopping_rounds = 100,
                                            eval_metric = "auc",
                                            objective = "binary:logistic"),
                          cv_folds = cv_folds, cp = xgb_cp, seed = seed, note = note,
                          vars_name = FALSE, save_data = TRUE)
    select_vars_xgb = xgb_list[,"Feature"]
    if (length(select_vars_xgb) <= 1) {
      select_vars_xgb = imp_list
      warning("There is no variable that meets the threshold of XGB selection.\n")
    }
    if (hopper) {
      imp_list = select_vars_xgb
    }
    xgb_list$Feature = gsub("\\.\\S{1,100}\\.", "", xgb_list$Feature)
    xgb_list = xgb_list %>% dplyr::group_by(Feature) %>%
      dplyr::summarise(Imp_Means_XGB = sum(Imp_Means_XGB))
    if (save_data) {
      save_dt(xgb_list, as_list = FALSE, row_names = FALSE, note = note,
              file_name = ifelse(is.null(file_name), "feature.filter.XGB", paste(file_name, "feature.filter.XGB", sep = ".")), dir_path = dir_path)
    }
    if (is.null(dt_imp)) {
      dt_imp = xgb_list
    } else {
      dt_imp = merge(xgb_list, dt_imp)
    }
  }

  imp_select_vars = unique(Reduce("intersect", list(select_vars_iv, select_vars_psi,
                                                    select_vars_cor, select_vars_xgb)))

  if (length(imp_select_vars) < 0) {
    imp_select_vars = unique(Reduce("union", list(select_vars_iv, select_vars_psi,
                                                  select_vars_cor, select_vars_xgb)))
    dat = de_one_hot_encoding(dat_train[imp_select_vars], cat_vars = char_x_list,
                              na_act = TRUE, note = FALSE)
    imp_vars = get_names(dat = dat,
                         types = c('character', 'factor', 'numeric', 'integer', 'double'),
                         ex_cols = c(target, ex_cols),
                         get_ex = FALSE)
  } else {
    dt_imp = dt_imp[order(dt_imp$IV, decreasing = TRUE),]
    imp_vars = as.character( dt_imp[ ,"Feature"])
  }


  if (save_data) {
    save_dt(imp_vars, as_list = TRUE, row_names = FALSE, note = note,
            file_name = ifelse(is.null(file_name), "feature.filter", paste(file_name, "feature.filter", sep = ".")), dir_path = dir_path)
    save_dt(dt_imp, as_list = FALSE, row_names = FALSE, note =  note,
            file_name = ifelse(is.null(file_name), "feature.filter.table", paste(file_name, "feature.filter.table", sep = ".")), dir_path = dir_path)
  }
  if (vars_name) {
    return(imp_vars)
  } else {
    return(dt_imp)
  }
}

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
#' @param save_data Logical, save results in locally specified folder. Default is FALSE.
#' @param file_name  The name for periodically saved results files. Default is "Feature_selected_COR".
#' @param dir_path The path for periodically saved results files. Default is "./variable".
#' @param ...  Additional parameters.
#' @return  A list of selected variables.
#' @seealso \code{\link{get_correlation_group}}, \code{\link{reduce_high_cor}}, \code{\link{char_cor_vars}}
#' @examples
#' # calculate iv for each variable.
#' iv_list = feature_select_wrapper(dat_train = UCICreditCard[1:1000,], dat_test = NULL,
#' target = "default.payment.next.month",
#' occur_time = "apply_date",
#' filter = c("IV"), cv_folds = 1, iv_cp = 0.01,
#' ex_cols = "ID$|date$|default.payment.next.month$",
#' save_data = FALSE, vars_name = FALSE)
#' fast_high_cor_filter(dat = UCICreditCard[1:1000,],
#' com_list = iv_list, save_data = FALSE,
#' ex_cols = "ID$|date$|default.payment.next.month$",
#' p = 0.9, cor_class = FALSE ,var_name = FALSE)
#' @export

fast_high_cor_filter <- function(dat, p = 0.7, x_list = NULL, com_list = NULL,
                                 ex_cols = NULL, save_data = FALSE, cor_class = TRUE,
                                 parallel = FALSE, note = FALSE,
                                 file_name = NULL, dir_path = tempdir(), ...) {
  if (note) cat("[NOTE] Fast dimension reduction for highly correlated variables. \n")
  dat = checking_data(dat)
  dat = time_transfer(dat)
  if (!is.null(x_list)) {
    dat = dat[, unique(c(x_list))]
  }
  dat <- low_variance_filter(dat, lvp = 1)
  ex_x_list = get_names(dat = dat,
                        types = c('factor', 'character', 'numeric', 'integer', 'double'),
                        ex_cols = ex_cols, get_ex = TRUE)

  if (is.null(com_list)) {
    stop("The comparison list is empty. \nFor comparisons between variables, IV or PSI or other index must be used.")
  } else {
    com_list = as.data.frame(com_list)
  }
  char_x_list = num_x_list = NULL
  if (cor_class) {
    char_x_list = get_names(dat = dat, types = c('factor', 'character'),
                            ex_cols = ex_cols, get_ex = FALSE)
    num_x_list = get_names(dat = dat, types = c('numeric', 'integer', 'double'),
                           ex_cols = ex_cols, get_ex = FALSE)
    if (length(num_x_list) > 2) {
      cor_mat_num = cor(dat[num_x_list], method = "spearman", use = "complete.obs")
      cor_nums <- reduce_high_cor(cor_mat = cor_mat_num, p = p,
                                  com_list = com_list, x_list = num_x_list)
    } else {
      cor_nums = num_x_list
    }
    if (length(num_x_list) > 2) {
      cor_mat_char = char_cor(dat = dat, x_list = char_x_list, parallel = parallel)
      cor_chars <- reduce_high_cor(cor_mat = cor_mat_char, p = p,
                                   com_list = com_list, x_list = char_x_list)
    } else {
      cor_chars = num_x_list
    }
    cor_vars <- unique(c(cor_chars, cor_nums))
    var_list = cor_vars
  } else {
    char_x_list = get_names(dat = dat, types = c('factor', 'character'),
                            ex_cols = ex_cols, get_ex = FALSE)
    num_x_list = get_names(dat = dat, types = c('numeric', 'integer', 'double'),
                           ex_cols = ex_cols, get_ex = FALSE)
    if (length(num_x_list) > 2) {
      cor_mat_num = cor(dat[num_x_list], method = "spearman", use = "complete.obs")
      cor_vars = reduce_high_cor(cor_mat = cor_mat_num, p = p,
                                 com_list = com_list, x_list = num_x_list)
    } else {
      cor_vars = num_x_list
    }
    var_list = c(char_x_list, cor_vars)
  }

  if (save_data) {
    dir_path = ifelse(!is.character(dir_path), tempdir(), dir_path)
    if (!dir.exists(dir_path)) dir.create(dir_path)
    if (!is.character(file_name)) { file_name = NULL }
    save_dt(var_list, file_name = ifelse(is.null(file_name), "feature.COR", paste(file_name, "feature.COR", sep = ".")), dir_path = dir_path, note = note, as_list = TRUE)
  }
  return(var_list)
}

#' @rdname fast_high_cor_filter
#' @export

high_cor_filter <- function(dat, com_list = NULL, x_list = NULL, ex_cols = NULL,
                            onehot = TRUE, parallel = TRUE, p = 0.7, file_name = NULL,
                            dir_path = tempdir(), save_data = FALSE, note = FALSE, ...) {
  if (note) cat("[NOTE] Dimension reduction for highly correlated variables. \n")
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
  ex_list <- get_names(dat = dat,
                       types = c('factor', 'character', 'numeric', 'integer', 'double'),
                       ex_cols = ex_cols, get_ex = TRUE)
  #obtain the numeric variables.
  num_x_list = get_names(dat = dat,
                         types = c('numeric', 'integer', 'double'),
                         ex_cols = ex_cols, get_ex = FALSE)
  #obtain the character or factor variables.
  char_x_list = get_names(dat = dat, types = c('factor', 'character'),
                          ex_cols = ex_cols, get_ex = FALSE)

  if (note) cat("[NOTE] Calculate the correlation matrix of numeric variables. \n")
  cor_mat_num = cor(dat[num_x_list], method = "spearman")
  #calculate the correlation matrix of character or factor variables.
  cor_mat_char = char_cor(dat = dat, x_list = char_x_list, parallel = parallel)
  # obtain highly correlated variable groups.
  group_vars <- c(get_correlation_group(cor_mat_num, p = p),
                  get_correlation_group(cor_mat_char, p = round(p / 1.5, 1)))
  group_len <- sapply(group_vars, function(x) length(x))
  single_group_vars <- unlist(group_vars[group_len == 1])
  multi_group_vars <- group_vars[group_len > 1]
  cat("[NOTE] Selecting the variable with the highest IV in a highly correlated variable group . \n")
  x = multi_group_vars[[1]]
  sel_vars <- vapply(multi_group_vars, function(x) {

    #In a highly correlated variable group, the variable with the highest IV  was selected.
    if (!is.null(com_list) & all(x %in% as.character(com_list[, 1]))) {
      x_group = com_list[which(as.character(com_list[, 1]) %in% x),]
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
    dir_path = ifelse(!is.character(dir_path), tempdir(), dir_path)
    if (!dir.exists(dir_path)) dir.create(dir_path)
    if (!is.character(file_name)) { file_name = NULL }
    save_dt(var_list, file_name = ifelse(is.null(file_name), "feature.COR.list", paste(file_name, "feature.COR.list", sep = ".")),
            dir_path = dir_path, note = note, as_list = TRUE)
    save_dt(group_vars, file_name = ifelse(is.null(file_name), "feature.COR.group", paste(file_name, "feature.COR.group", sep = ".")),
            dir_path = dir_path, note = note, as_list = TRUE)
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
    cor_vars_list2[[i]] <- lapply(cor_vars_list, function(x) base::setdiff(unlist(x), unlist(cor_vars_list[i])))
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

#' Variable selection by LASSO
#'
#' @description  \code{lasso_filter} filter variables by lasso.
#'
#' @param dat_train  A data.frame with independent variables and target variable.
#' @param dat_test  A data.frame of test data. Default is NULL.
#' @param target The name of target variable.
#' @param x_list Names of independent variables.
#' @param ex_cols A list of excluded variables. Regular expressions can also be used to match variable names. Default is NULL.
#' @param pos_flag The value of positive class of target variable, default: "1".
#' @param sim_sign The coefficients of all variables should be all negetive or positive, after turning to woe. Default is "negetive" for pos_flag is "1".
#' @param best_lambda  Metheds of best lambda stardards using to filter variables by LASSO. There are 3 methods: ("lambda.auc", "lambda.ks", "lambda.sim_sign") . Default is  "lambda.auc".
#' @param plot.it Logical, shrinkage plot. Default is TRUE.
#' @param seed  Random number seed. Default is 46.
#' @param save_data Logical, save results in locally specified folder. Default is FALSE
#' @param file_name  The name for periodically saved results files. Default is "Feature_selected_LASSO".
#' @param dir_path The path for periodically saved results files. Default is "./variable".
#' @return A list of filtered x variables by lasso.
#' @examples
#'  sub = cv_split(UCICreditCard, k = 40)[[1]]
#'  dat = UCICreditCard[sub,]
#'  dat = re_name(dat, "default.payment.next.month", "target")
#'  dat = data_cleansing(dat, target = "target", obs_id = "ID", occur_time = "apply_date",
#'   miss_values = list("", -1))
#'
#'  train_test <- train_test_split(dat, split_type = "OOT", prop = 0.7,
#'                                   occur_time = "apply_date")
#'  dat_train = train_test$train
#'  dat_test = train_test$test
#'  #get breaks of all predictive variables
#'  x_list = c("PAY_0", "LIMIT_BAL", "PAY_AMT5", "EDUCATION", "PAY_3", "PAY_2")
#'  breaks_list <- get_breaks_all(dat = dat_train, target = "target",
#'                                 x_list = x_list, occur_time = "apply_date", ex_cols = "ID",
#'   save_data = FALSE, note = FALSE)
#'  #woe transform
#'  train_woe = woe_trans_all(dat = dat_train,
#'                             target = "target",
#'                             breaks_list = breaks_list,
#'                             woe_name = FALSE)
#'  test_woe = woe_trans_all(dat = dat_test,
#'                          target = "target",
#'                            breaks_list = breaks_list,
#'                            note = FALSE)
#'  lasso_filter(dat_train = train_woe, dat_test = test_woe,
#'          target = "target", x_list = x_list,
#'        save_data = FALSE, plot.it = TRUE)
#' @importFrom glmnet cv.glmnet glmnet
#' @import ggplot2
#' @importFrom gridExtra arrangeGrob
#' @export


lasso_filter <- function(dat_train, dat_test = NULL, target = NULL,
                         x_list = NULL, pos_flag = NULL,
                         ex_cols = NULL, sim_sign = "negtive",
                         best_lambda = "lambda.auc",
                         save_data = FALSE, plot.it = TRUE, seed = 46,
                         file_name = NULL, dir_path = tempdir()) {
  cat("[NOTE] variable selection by LASSO. \n")
  opt = options(scipen = 200, "warn" = -1, stringsAsFactors = FALSE, digits = 10) # suppress warnings
  dat_train = checking_data(dat = dat_train, target = target, pos_flag = pos_flag)
  x_list = get_x_list(x_list = x_list, dat_train = dat_train, dat_test = dat_test,
                      ex_cols = c(target, ex_cols))
  var_list = unique(c(target, x_list))
  if (!is.null(dat_test)) {
    dat_ts = rbind(dat_train[, var_list], dat_test[, var_list])
  } else {
    dat_ts = dat_train
  }

  dat_ts = low_variance_filter(dat = dat_ts, lvp = 1, note = FALSE)
  if (any(is.na(dat_ts[x_list]))) {
    dat_ts = process_nas(dat = dat_ts, x_list = var_list,
                         ex_cols = c(target, ex_cols),
                         method = "median")
  }

  char_x_list = get_names(dat = dat_ts, types = c('character', 'factor'),
                          ex_cols = c(target, ex_cols), get_ex = FALSE)
  if (length(char_x_list) > 0) {
    dat_ts = one_hot_encoding(dat = dat_ts, cat_vars = char_x_list, na_act = FALSE, note = FALSE)
  }
  num_x_list = get_names(dat = dat_ts, types = c('numeric', 'integer', 'double'),
                         ex_cols = c(target, ex_cols), get_ex = FALSE)

  if (is.null(dat_test)) {
    train_test = train_test_split(dat_ts, split_type = "Random", prop = 0.3,
                                  seed = 46, note = FALSE, save_data = FALSE)
  } else {
    nr = nrow(dat_train)
    train_test = train_test_split(dat_ts, split_type = "byRow", prop = nr / nrow(dat_ts),
                                  seed = 46, note = FALSE, save_data = FALSE)
  }
  dat_train = train_test$train
  dat_test = train_test$test
  #cross validation can also be used to select lambda.
  x_train = as.matrix(dat_train[, num_x_list])
  y_train = as.numeric(as.character(dat_train[, target]))
  x_test = as.matrix(dat_test[, num_x_list])
  y_test = as.numeric(as.character(dat_test[, target]))

  # Train a model
  if (!is.null(seed)) set.seed(seed) else set.seed(46)
  lasso_model <- glmnet(y = y_train,
                        x = x_train, #X must be a matrix.
                        family = "binomial", #For non-negative count dependent variable.
                        alpha = 1) #1for lasso
  bst_lambda = bst_lambda_id = best_ks =  bst_coefficients = bst_vars = shrink_vars = c()
  lambda_ks_auc = get_auc_ks_lambda(lasso_model, x_test = x_test, y_test = y_test, save_data = save_data, plot_show = plot.it, file_name = file_name, dir_path = dir_path)

  lasso_lambda_sim = get_sim_sign_lambda(lasso_model = lasso_model,
                                         sim_sign = sim_sign)

  if (best_lambda == "lambda.auc" | !is.element(best_lambda, c("lambda.auc", "lambda.ks", "lambda.sim_sign"))) {
    b_lambda = lambda_ks_auc$lambda.auc
  }
  if (best_lambda == "lambda.ks") {
    b_lambda = lambda_ks_auc$lambda.ks
  }

  if (best_lambda == "lambda.sim_sign") {
    b_lambda = max(lambda_ks_auc[["lambda.auc"]],
                   lambda_ks_auc[["lambda.ks"]],
                   lasso_lambda_sim)
  }
  coefficients = coef(lasso_model, s = b_lambda)
  #variable coefficients  which coefficient are not zero
  bst_vars = which(coefficients[-1] != 0)
  if (length(bst_vars) < 2) {
    b_lambda = max(lambda_ks_auc[["lambda.auc"]],
                   lambda_ks_auc[["lambda.ks"]])
    coefficients = coef(lasso_model, s = b_lambda)
    bst_vars = which(coefficients[-1] != 0)
  }
  if (length(bst_vars) < 2) {
    b_lambda = lambda_ks_auc[["lambda.auc"]]
    coefficients = coef(lasso_model, s = b_lambda)
    bst_vars = which(coefficients[-1] != 0)
  }
  shrinkage_vars = num_x_list[bst_vars]
  dat_lasso = de_one_hot_encoding(dat_one_hot = dat_train[shrinkage_vars],
                                  cat_vars = char_x_list,
                                  na_act = TRUE, note = FALSE)
  lasso_vars = get_x_list(x_list = x_list, dat_train = dat_lasso,
                          ex_cols = c(target, ex_cols))
  if (save_data) {
    dir_path = ifelse(!is.character(dir_path), tempdir(), dir_path)
    if (!dir.exists(dir_path)) dir.create(dir_path)
    if (!is.character(file_name)) { file_name = NULL }
    save_dt(lasso_vars, file_name = ifelse(is.null(file_name), "feature.filter.LASSO", paste(file_name, "feature.filter.LASSO", sep = ".")), dir_path = dir_path, as_list = TRUE)
  }
  options(opt) # reset warnings
  return(lasso_vars)
}


#' get_auc_ks_lambda
#' \code{get_auc_ks_lambda} is for get best lambda required in lasso_filter. This function required in \code{lasso_filter}
#' @param lasso_model A lasso model genereted by glmnet.
#' @param x_test A matrix of test dataset with x.
#' @param y_test A matrix of y test dataset with y.
#' @param plot_show Logical, if TRUE plot the results. Default is TRUE.
#' @param save_data Logical, save results in locally specified folder. Default is FALSE
#' @param file_name  The name for periodically saved results files. Default is NULL.
#' @param dir_path The path for periodically saved results files.
#' @return Lanmbda values with max K-S and AUC.
#' @seealso \code{\link{lasso_filter}}, \code{\link{get_sim_sign_lambda}}
#' @importFrom glmnet cv.glmnet glmnet
#' @import ggplot2
#' @importFrom gridExtra arrangeGrob
#' @export
get_auc_ks_lambda <- function(lasso_model, x_test, y_test, save_data = FALSE, plot_show = TRUE, file_name = NULL, dir_path = tempdir()) {
  test_pre = predict(lasso_model, newx = x_test, s = lasso_model$lambda, type = "response") # make predictions
  KS = apply(test_pre, 2, function(x) round(ks_value(score = x, target = y_test, g = 20), 4))
  AUC = apply(test_pre, 2, function(x) round(auc_value(target = y_test, prob = x), 4))

  best_ks = max(KS, na.rm = TRUE)[1]
  bst_ks_lambda = lasso_model$lambda[which(KS == best_ks)[1]]
  KS_lambda = data.frame(vars_num = lasso_model$df,
                         KS = KS,
                         lambda = lasso_model$lambda)
  ks_lanm = KS_lambda[-1,]

  best_auc = max(AUC, na.rm = TRUE)[1]
  bst_auc_lambda = lasso_model$lambda[which(AUC == best_auc)[1]]

  AUC_lambda = data.frame(vars_num = lasso_model$df,
                          AUC = AUC,
                          lambda = lasso_model$lambda)
  auc_lanm = AUC_lambda[-1,]

  plot_ks = ggplot(data = ks_lanm, aes(x = -log(ks_lanm$lambda), y = ks_lanm$KS)) +
    geom_line(aes(color = "KS"),
              position = position_dodge(width = 0.5), #color = love_color("sky_blue"),
              size = 1) +
    geom_point(position = position_dodge(width = 0.5),
               fill = 'white', color = love_color("deep_red"),
               size = 1, shape = 21) +
    geom_line(aes(y = -ks_lanm$lambda[seq(1, by = ceiling(length(ks_lanm$vars_num) / 20),
                                          length.out = length(ks_lanm$vars_num))] * max(ks_lanm$KS) * 10 + max(ks_lanm$KS),
                  x = -log(ks_lanm$lambda)[seq(1, by = ceiling(length(ks_lanm$vars_num) / 20),
                                               length.out = length(ks_lanm$vars_num))],
                  color = "Number of Variables"),
              position = position_dodge(width = 0.5), size = 1) +
    geom_point(aes(y = -ks_lanm$lambda[seq(1, by = ceiling(length(ks_lanm$vars_num) / 20),
                                           length.out = length(ks_lanm$vars_num))] * max(KS) * 10 + max(KS),
                   x = -log(ks_lanm$lambda)[seq(1, by = ceiling(length(ks_lanm$vars_num) / 20),
                                                length.out = length(ks_lanm$vars_num))]),
               position = position_dodge(width = 0.5),
               fill = 'white', color = love_color("light_grey"), size = 1, shape = 21) +
    geom_text(aes(y = -ks_lanm$lambda[seq(1, by = ceiling(length(ks_lanm$vars_num) / 20),
                                          length.out = length(ks_lanm$vars_num))] * max(KS) * 10 + max(ks_lanm$KS),
                  x = -log(ks_lanm$lambda)[seq(1, by = ceiling(length(ks_lanm$vars_num) / 20),
                                               length.out = length(ks_lanm$vars_num))],
                  label = paste(ks_lanm$vars_num[seq(1, by = ceiling(length(ks_lanm$vars_num) / 20),
                                                     length.out = length(ks_lanm$vars_num))])),
              position = position_dodge(width = 0.5),
              colour = 'black', size = 2.5, vjust = -1) +
    annotate(geom = 'text',
             x = -log(bst_ks_lambda),
             y = best_ks * 1.1,
             label = paste("max K-S", best_ks, sep = "\n")) +
    scale_color_manual(values = c('KS' = love_color("sky_blue"),
                                  'Number of Variables' = love_color("light_yellow"))) +
    ylab("K-S") +
    ggtitle("LASSO: K-S with different lambda") +
    plot_theme(legend.position = "top", angle = 0)

  plot_auc = ggplot(data = auc_lanm, aes(x = -log(auc_lanm$lambda), y = auc_lanm$AUC)) +
    geom_line(aes(color = "AUC"),
              position = position_dodge(width = 0.5), #color = love_color("sky_blue"),
              size = 1) +
    geom_point(position = position_dodge(width = 0.5),
               fill = 'white', color = love_color("deep_red"),
               size = 1, shape = 21) +
    geom_line(aes(y = -auc_lanm$lambda[seq(1, by = ceiling(length(auc_lanm$vars_num) / 20),
                                           length.out = length(auc_lanm$vars_num))] * max(auc_lanm$AUC) * 10 + max(auc_lanm$AUC),
                  x = -log(auc_lanm$lambda)[seq(1, by = ceiling(length(auc_lanm$vars_num) / 20),
                                                length.out = length(auc_lanm$vars_num))],
                  color = "Number of Variables"),
              position = position_dodge(width = 0.5), size = 1) +
    geom_point(aes(y = -auc_lanm$lambda[seq(1, by = ceiling(length(auc_lanm$vars_num) / 20),
                                            length.out = length(auc_lanm$vars_num))] * max(AUC) * 10 + max(AUC),
                   x = -log(auc_lanm$lambda)[seq(1, by = ceiling(length(auc_lanm$vars_num) / 20),
                                                 length.out = length(auc_lanm$vars_num))]),
               position = position_dodge(width = 0.5),
               fill = 'white', color = love_color("light_grey"), size = 1, shape = 21) +
    geom_text(aes(y = -auc_lanm$lambda[seq(1, by = ceiling(length(auc_lanm$vars_num) / 20),
                                           length.out = length(auc_lanm$vars_num))] * max(AUC) * 10 + max(auc_lanm$AUC),
                  x = -log(auc_lanm$lambda)[seq(1, by = ceiling(length(auc_lanm$vars_num) / 20),
                                                length.out = length(auc_lanm$vars_num))],
                  label = paste(auc_lanm$vars_num[seq(1, by = ceiling(length(auc_lanm$vars_num) / 20),
                                                      length.out = length(auc_lanm$vars_num))])),
              position = position_dodge(width = 0.5),
              colour = 'black', size = 2.5, vjust = -1) +
    annotate(geom = 'text',
             x = -log(bst_auc_lambda),
             y = best_auc * 1.1,
             label = paste("max AUC", best_auc, sep = " \n")) +
    scale_color_manual(values = c('AUC' = love_color("deep_purple"),
                                  'Number of Variables' = love_color("deep_green"))) +
    ylab("AUC") +
    ggtitle("LASSO: AUC with different lambda") +
    plot_theme(legend.position = "top", angle = 0)
  if (save_data) {
    dir_path = ifelse(!is.character(dir_path), tempdir(), dir_path)
    if (!dir.exists(dir_path)) dir.create(dir_path)
    if (!is.character(file_name)) { file_name = NULL }
    save_dt(KS_lambda, file_name = ifelse(is.null(file_name), "LASSO.lambda.KS",
                                          paste(file_name, "LASSO.lambda.KS", sep = ".")), dir_path = dir_path)
    save_dt(AUC_lambda, file_name = ifelse(is.null(file_name), "LASSO.lambda.AUC", paste(file_name, "LASSO.lambda.AUC", sep = ".")), dir_path = dir_path)
    ggsave(device = "png",
           filename = paste0(dir_path, "/", ifelse(is.null(file_name), "LASSO.lambda.KS_AUC",
                                                   paste(file_name, "LASSO.lambda.KS_AUC", sep = ".")), ".png"),
           plot = arrangeGrob(grobs = list(plot_ks, plot_auc), ncol = 2, nrow = 1), dpi = "retina", width = 8)
  }
  if (plot_show) {
    plot(arrangeGrob(grobs = list(plot_ks, plot_auc), ncol = 2, nrow = 1))
  }
  return(list(lambda.ks = bst_ks_lambda, lambda.auc = bst_auc_lambda))
}

#' get_sim_sign_lambda
#' \code{get_sim_sign_lambda} is for get Best lambda required in lasso_filter. This function required in \code{lasso_filter}
#' @param lasso_model A lasso model genereted by glmnet.
#' @param sim_sign  Default is "negtive". This is related to pos_plag. If pos_flag equals 1 or 1, the value must be set to negetive. If pos_flag equals 0 or 0, the value must be set to positive.
#' @return Lanmbda value
#' @details
#' lambda.sim_sign give the model with the same positive or negetive coefficients of all variables.
#'
#' @export
get_sim_sign_lambda <- function(lasso_model, sim_sign = "negtive") {
  lambda = lasso_model$lambda

  coefs <- c()
  if (sim_sign == "negtive") {
    for (i in 1:length(lambda)) {
      coefs[i] = coef(lasso_model, s = lambda[i])[-1] > 0
    }
  } else {
    for (i in 1:length(lambda)) {
      coefs[i] = coef(lasso_model, s = lambda[i])[-1] < 0
    }
  }
  if (length(coefs) > 0) {
    lambda.sim_sign = lambda[sum(!coefs)]
  } else {
    stop(paste("no variables coeficeient is positive or negetive"))
  }
  lambda.sim_sign
}

#' PCA Dimension Reduction
#'
#' \code{PCA_reduce} is used for PCA reduction of high demension data .
#' @param train A data.frame with independent variables and target variable.
#' @param test  A data.frame of test data.
#' @param mc   Threshold of cumulative imp.
#' @examples
#' \dontrun{
#' num_x_list = get_names(dat = UCICreditCard, types = c('numeric'),
#' ex_cols = "ID$|date$|default.payment.next.month$", get_ex = FALSE)
#'  PCA_dat = PCA_reduce(train = UCICreditCard[num_x_list])
#' }
#' @export
PCA_reduce <- function(train = train, test = NULL, mc = 0.9){
  train = train[which(complete.cases(train)),]
  train_mean = apply(train, 2, mean)
  train_std = apply(train, 2, sd)
  dat = as.matrix((train - train_mean) / train_std)
  cov_data = cov(dat)
  eigen_data = eigen(cov_data)
  eigen_value = eigen_data$values
  eigen_vector = eigen_data$vectors
  order_value = order(eigen_value, decreasing = T)
  values = eigen_value[order_value]
  value_sum = sum(values)
  cum_var = cumsum(values) / value_sum
  k_order_value = 1
  for (i in order_value) {
    if (cum_var[i] >= mc) break
    k_order_value = which(order_value <= i + 1)
  }
  order_vector = eigen_vector[, k_order_value]
  if (is.null(test)) {
    principal = dat %*% order_vector
  } else {
    test = test[which(complete.cases(test))]
    test_data = as.matrix((test - train_mean) / train_std)
    principal = test_data %*% order_vector
  }
  PCA_data = data.frame(PCA = principal)
  mypaste = function(x, y) paste(x, y, sep = ".")
  sp = ""
  names_strsplit = lapply(names(train), function(x) strsplit(x, sp)[[1]])
  for (i in 1:length(PCA_data)) {
    names(PCA_data)[i] = paste(Reduce("mypaste", Reduce("union", names_strsplit)),
                               names(PCA_data)[i], sep = "_")
  }
  return(PCA_data)
}
