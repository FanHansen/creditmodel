#'Training model
#'
#' \code{training_model} Model builder
#' @param model_name  A string, name of the project. Default is NULL..
#' @param dat_train A data.frame with independent variables and target variable.
#' @param dat_test  A data.frame of test data. Default is NULL.
#' @param target The name of target variable. 
#' @param x_list Names of independent variables. Default is NULL.
#' @param occur_time The name of the variable that represents the time at which each observation takes place.Default is NULL.
#' @param pos_flag The value of positive class of target variable, default: "1".
#' @param ex_cols Names of excluded variables. Regular expressions can also be used to match variable names. Default is NULL.
#' @param prop Percentage of train-data after the partition. Default: 0.7.
#' @param obs_id  The name of ID of observations or key variable of data. Default is NULL.
#' @param miss_values  Other extreme value might be used to represent missing values, e.g: -9999, -9998. These miss_values will be encoded to -1 or "Unknown".
#' @param preproc Logical. Preprocess data. Default is TRUE
#' @param outlier_proc Logical. If TRUE,  Outliers processing using Kmeans and Local Outlier Factor. Default is TRUE
#' @param missing_proc Logical. If TRUE, missing value analysis and process missing value by knn imputation or central impulation or random imputation. Default is TRUE
#' @param default_miss Logical. If TRUE, assigning the missing values to -1 or "Unknown", otherwise, processing the missing values according to the results of missing analysis. See details at: \code{\link{process_nas}}
#' @param feature_filter  Parameters for selecting important and stable features.See details at: \code{\link{feature_select_wrapper}}
#' @param algorithm  Algorithms for training a model. list("LR", "XGB", "GBDT", "RF") are available.
#' @param LR.params  Parameters of logistic regression & scorecard. See details at :  \code{\link{lr_params}}. 
#' \itemize{
#'   \item \code{tree_control} the list of parameters to control cutting initial breaks by decision tree. See details at: \code{\link{get_tree_breaks}}
#'   \item \code{bins_control} the list of parameters to control merging initial breaks. See details at: \code{\link{select_best_breaks}},\code{\link{select_best_class}}
#'   \item \code{best_lambda} Metheds of best lanmbda stardards using to filter variables by LASSO.There are four methods: ("lambda.min", "lambda.1se", "lambda.05se" , "lambda.sim_sign") . Default is  "lambda.sim_sign". See details at: \code{\link{get_best_lambda}}
#'   \item \code{obsweight} An optional vector of  'prior weights' to be used in the fitting process. Should be NULL or a numeric vector. If you oversample or cluster diffrent datasets to training the LR model, you need to set this parameter to ensure that the probability of logistic regression output is the same as that before oversampling or segmentation. e.g.:There are 10,000 good obs and 500 bad obs before oversampling or under-sampling, 5,000 good obs and 3,000 bad obs after oversampling. Then this parameter should be set to c(10000/5000, 500/3000). Default is NULL..
#'   \item \code{forced_in}Names of forced input variables. Default is NULL.
#'   \item \code{sp_values} Vaules will be in separate bins.e.g. list(-1, "Unknown")  means that -1 & Unknown as special values.Default is NULL.
#'   \item \code{step_wise} Logical, stepwise method. Default is TRUE.
#'   \item \code{score_card} Logical, transfer woe to a standard scorecard. If TRUE, Output scorecard, and score prediction, otherwise output probability. Default is TRUE.
#'   \item \code{cor_p} The maximum threshold of correlation.0 <= cor_p <=1; 0.5 to 0.8 usually work.  Default: 0.7.
#'   \item \code{iv_i} The minimum threshold of IV. 0 < iv_i ; 0.01 to 0.1 usually work. Default: 0.01
#'   \item \code{psi_i} The maximum threshold of PSI.  0 <= psi_i <=1; 0.05 to 0.2 usually work. Default: 0.1
#' }
#' @param XGB.params Parameters of xgboost. See details at :  \code{\link{xgb_params}}. 
#' @param GBM.params Parameters of GBM. See details at :  \code{\link{gbm_params}}. 
#' @param RF.params  Parameters of Random Forest. See details at :  \code{\link{rf_params}}. 
#' @param breaks_list A table containing a list of splitting points for each independent variable. Default is NULL.
#' @param parallel  Default is FALSE
#' @param cores_num The number of CPU cores to use. 
#' @param save_pmml Logical, save model in PMML format. Default is TRUE.
#' @param plot_show  Logical, show model performance in current graphic device. Default is FALSE.
#' @param ...  Other parameters.
#' @param seed  Random number seed. Default is 46.
#' @return A list containing Model Objects.
#' @seealso   \code{\link{train_test_split}},\code{\link{cleaning_data}}, \code{\link{feature_select_wrapper}},   \code{\link{lr_params}}, \code{\link{xgb_params}}, \code{\link{gbm_params}}, \code{\link{rf_params}},\code{\link{fast_high_cor_filter}},\code{\link{get_breaks_all}},\code{\link{lasso_filter}}, \code{\link{woe_trans_all}}, \code{\link{get_logistic_coef}}, \code{\link{score_transfer}},\code{\link{get_score_card}}, \code{\link{model_key_index}},\code{\link{ks_psi_plot}},\code{\link{get_plots}},\code{\link{ks_table_plot}}
#' @examples
#' \dontrun{
#' B_model = training_model(dat_train = lendingclub,
#' model_name = "lendingclub",target = "loan_status",
#' occur_time = "issue_d", obs_id = "id", prop = 0.7, 
#' algorithm = list("LR", "XGB"),seed = 46)
#' }
#' @import ggplot2 
#' @importFrom gridExtra arrangeGrob
#' @importFrom foreach foreach
#' @importFrom gbm gbm gbm.perf
#' @importFrom xgboost xgb.importance xgb.train xgb.DMatrix xgb.dump
#' @importFrom dplyr group_by mutate summarize  summarise n  count %>% filter  left_join
#' @importFrom data.table fwrite fread dcast melt
#' @importFrom randomForest randomForest tuneRF combine
#' @importFrom XML saveXML
#' @importFrom pmml pmml
#' @importFrom glmnet cv.glmnet glmnet
#' @importFrom grDevices dev.print png rgb
#' @importFrom graphics par plot text
#' @importFrom stats  ar as.formula binomial chisq.test coef complete.cases cor cov glm  kmeans median  na.omit  na.pass predict reorder runif sd ts
#' @importFrom utils  capture.output  data  install.packages  installed.packages
#' @export


training_model <- function(model_name = "model", dat_train, dat_test = NULL,
                           target, pos_flag = NULL, x_list = NULL, prop = 0.7,
                           occur_time = NULL, obs_id = NULL, miss_values = NULL, ex_cols = NULL, 
                           preproc = TRUE, outlier_proc = TRUE, missing_proc = TRUE, default_miss = FALSE,
                           feature_filter = list(  filter = c("IV", "PSI", "COR", "XGB"),
                                                 cv_folds = 1, iv_cp = 0.02, psi_cp = 0.1, xgb_cp = 0, hopper = FALSE),
                           algorithm = list("LR", "XGB"),
                          LR.params = lr_params(), XGB.params = xgb_params(),
                          GBM.params = gbm_params(), RF.params = rf_params(),
                          breaks_list = NULL, parallel = FALSE,
                         cores_num = NULL, save_pmml = TRUE, plot_show = FALSE,  seed = 46, ...) {

    opt = options(scipen = 100, stringsAsFactors = FALSE, digits = 10) 
    if (!dir.exists("./model")) dir.create("./model")
    if (!dir.exists("./data")) dir.create("./data")
    if (!dir.exists("./variable")) dir.create("./variable")
    if (!dir.exists("./performance")) dir.create("./performance")
    if (!dir.exists("./predict")) dir.create("./predict")

    if (is.null(algorithm)) {
        stop(paste("algorithm is missing.\n"))
    }
    if (!any(is.element(algorithm, c("LR", "XGB", "GBDT", "RF")))) {
        stop("In algorithm, only  LR, XGB, GBDT, RF are supported.\n")
    }
    if (length(x_list) > 0 && any(x_list == target)) {
        stop(paste("x_list  contains", target, ".\n"))
    }
    if (is.null(dat_test)) {
        train_test <- train_test_split(dat_train, split_type = "OOT", prop = prop, occur_time = occur_time, seed = seed, save_data = TRUE, dir_path = "./data")
        dat_train = train_test$train
        dat_test = train_test$test
    }
    x_list = get_x_list(x_list = x_list, dat_train = dat_train, dat_test = dat_test, ex_cols = c(target, obs_id, occur_time, ex_cols))
    if (preproc) {
        dat_train = cleaning_data(dat = dat_train, target = target, x_list = x_list, obs_id = obs_id, occur_time = occur_time,
        pos_flag = pos_flag, miss_values = miss_values, ex_cols = ex_cols, outlier_proc = outlier_proc,
        missing_proc = missing_proc, default_miss = default_miss, low_var = TRUE, parallel = parallel, save_data = TRUE,
        dir_path = "./data", file_name = "dat_train_cl")
        dat_test = cleaning_data(dat = dat_test, target = target, x_list = x_list, obs_id = obs_id, occur_time = occur_time,
        pos_flag = pos_flag, miss_values = miss_values, ex_cols = ex_cols, outlier_proc = FALSE,note = FALSE,
        missing_proc = TRUE, default_miss = TRUE, low_var = FALSE, parallel = parallel, save_data = TRUE,
        dir_path = "./data", file_name = "dat_test_cl")
        x_list = get_x_list(x_list = x_list, dat_train = dat_train, dat_test = dat_test, ex_cols = c(target, obs_id, occur_time, ex_cols))
    }
    #prepare parallel computing
    if (parallel) {
        parallel <- start_parallel_computing(parallel)
        stopCluster <- TRUE
    } else {
        parallel <- stopCluster <- FALSE
    }
    if (is.null(cores_num)) {
        cores_num = parallel::detectCores()
    }
    on.exit(if (parallel & stopCluster) stop_parallel_computing(attr(parallel, "cluster")))

    if (!is.null(feature_filter)) {
        sel_list = NULL
        if (!is.null(feature_filter[["filter"]]) && any(is.element(feature_filter[["filter"]], c("IV", "PSI", "XGB","COR")))) {
            filter = feature_filter[["filter"]]
        } else {
            filter = c("IV", "PSI", "XGB","COR")
        }
        cv_folds = ifelse(!is.null(feature_filter[["cv_folds"]]), feature_filter[["cv_folds"]], 1)
        iv_cp = ifelse(!is.null(feature_filter[["iv_cp"]]), feature_filter[["iv_cp"]], 0.01)
        psi_cp = ifelse(!is.null(feature_filter[["psi_cp"]]), feature_filter[["psi_cp"]], 0.1)
        xgb_cp = ifelse(!is.null(feature_filter[["xgb_cp"]]), feature_filter[["xgb_cp"]], 0)
        cor_cp = ifelse(!is.null(feature_filter[["cor_cp"]]), feature_filter[["cor_cp"]], 0.98)
        hopper = ifelse(!is.null(feature_filter[["hopper"]]), feature_filter[["hopper"]], TRUE)
        sel_list = feature_select_wrapper(dat_train = dat_train, dat_test = dat_test, target = target,
        pos_flag = pos_flag, x_list = x_list, occur_time = occur_time, ex_cols = c(obs_id, occur_time, target, ex_cols),
        filter = filter, cv_folds = cv_folds, iv_cp = iv_cp, psi_cp = psi_cp, cor_cp = cor_cp, xgb_cp = xgb_cp, hopper = hopper,
        vars_name = TRUE, save_data = TRUE, parallel = parallel, cores_num = cores_num, seed = seed, file_name = paste0(model_name, ".Feature_Select_", paste(filter, collapse = "_")), dir_path = "./variable")
        if (length(sel_list) > 0) {
            x_list = get_x_list(x_list = sel_list, dat_train = dat_train, dat_test = dat_test, ex_cols = c(target, obs_id, occur_time, ex_cols))
        } else {
            warning("No feature satisfies the criteria for feature selection, use the previous x_list.\n")
        }
    }
    model_new = lr_model_new = xgb_model_new = gbdt_model_new = rf_model_new = NULL
    # Train Models
    if (length(unique(dat_train[, target])) == 2) {
        if (any(algorithm == "LR")) {
            cat("[NOTE] Training Logistic Regression Model.\n")
            LR_dir_path = "LR"
            LR_model_dir_path = paste0("./model/", LR_dir_path)
            LR_var_dir_path = paste0("./variable/", LR_dir_path)
            LR_perf_dir_path = paste0("./performance/", LR_dir_path)
            LR_pre_dir_path = paste0("./predict/", LR_dir_path)
            LR_data_dir_path = paste0("./data/", LR_dir_path)
       
            if (!dir.exists(LR_model_dir_path)) dir.create(LR_model_dir_path)
            if (!dir.exists(LR_var_dir_path)) dir.create(LR_var_dir_path)
            if (!dir.exists(LR_perf_dir_path)) dir.create(LR_perf_dir_path)
            if (!dir.exists(LR_pre_dir_path)) dir.create(LR_pre_dir_path)
            if (!dir.exists(LR_data_dir_path)) dir.create(LR_data_dir_path)
            if (is.null(LR.params)) {
                LR.params = lr_params()
            }
            obsweight = ifelse(!is.null(LR.params[["obsweight"]]), LR.params[["obsweight"]], 1)
            step_wise = ifelse(!is.null(LR.params[["step_wise"]]), LR.params[["step_wise"]], TRUE)
            cor_p = ifelse(!is.null(LR.params[["cor_p"]]), LR.params[["cor_p"]], 0.7)
            best_lambda = ifelse(!is.null(LR.params[["best_lambda"]]), LR.params[["best_lambda"]], "lambda.sim_sign")
            sp_values =LR.params[["sp_values"]]
            forced_in =LR.params[["forced_in"]]
            iv_i = ifelse(!is.null(LR.params[["iv_i"]]), LR.params[["iv_i"]], 0.01)
            psi_i = ifelse(!is.null(LR.params[["psi_i"]]), LR.params[["psi_i"]], 0.2)
            score_card = ifelse(!is.null(LR.params[["score_card"]]), LR.params[["score_card"]], TRUE)
            tree_control = if (!is.null(LR.params["tree_control"])) {
                LR.params[["tree_control"]]
            } else {
                list(p = 0.02, cp = 0.00000001, xval = 5, maxdepth = 10)
            }

          bins_control = if (!is.null(LR.params["bins_control"])) {
                LR.params[["bins_control"]]
            } else {
                list(bins_num = 10, bins_pct = 0.05, b_chi = 0.02, b_odds = 0.1, b_psi = 0.03, b_gb = 0.15, mono = 0.2,gb_psi = 0.15,kc = 1)
            }
            if (length(obsweight) > 1 && is.vector(obsweight, mode = "numeric")) {
                obsweighted <- ifelse(dat_train[, target] == 0, obsweight[1], obsweight[2])
            } else {
                obsweighted = NULL
            }
           bins_params = bins_control

            #bins_table
            if (is.null(breaks_list) && length(breaks_list) < 2) {
                breaks_list = get_breaks_all(dat = dat_train, occur_time = occur_time, oot_pct = prop, x_list = x_list,
                ex_cols = c(obs_id, occur_time, target, ex_cols), target = target, pos_flag = pos_flag,
                tree_control = tree_control, bins_control = bins_control, sp_values = sp_values,
                parallel = parallel, note = TRUE, file_name = paste0(model_name, ".Breaks_List"), dir_path = LR_var_dir_path)
            }
            #------------------------------------------
            #psi_iv_filter          
            iv_psi_list = psi_iv_filter(dat = dat_train, dat_test = dat_test, x_list = x_list, target = target,
            occur_time = occur_time, oot_pct = prop, pos_flag = pos_flag, parallel = parallel,
            ex_cols = c(obs_id, occur_time, target, ex_cols), breaks_list = breaks_list, psi_i = psi_i, iv_i = iv_i,
            save_data = TRUE, file_name = paste0(model_name, ".Feature_Imp_PSI_IV"), dir_path = LR_var_dir_path)
            if (length(iv_psi_list) < 1) {
                stop(paste("No variable satisfies the psi & iv condition."))
            }
            iv_psi_select_vars = as.character(iv_psi_list[, "Feature"])

            save_dt(iv_psi_select_vars, as_list = TRUE, row_names = FALSE, file_name = paste0(model_name, ".Feature_select_IV_PSI"), dir_path = LR_var_dir_path)
    
            #--------------------------------------------------------
            #woe transform
            train_woe = woe_trans_all(dat = dat_train, x_list = unlist(iv_psi_select_vars), target = target,
            ex_cols = c(target, occur_time, obs_id), bins_table = NULL, breaks_list = breaks_list,
            file_name = "train_woe", dir_path = LR_data_dir_path, note = TRUE, save_data = TRUE, woe_name = FALSE, parallel = parallel)
            test_woe = woe_trans_all(dat = dat_test, x_list = unlist(iv_psi_select_vars), target = target,
            ex_cols = c(target, occur_time, obs_id), bins_table = NULL, breaks_list = breaks_list,
            file_name = "test_woe", dir_path = LR_data_dir_path, note = FALSE, woe_name = FALSE, save_data = TRUE, parallel = parallel)

            cor_select_vars = fast_high_cor_filter(dat = train_woe, x_list = iv_psi_select_vars, com_list = iv_psi_list,
            save_data = TRUE, ex_cols = c(target, occur_time, obs_id, ex_cols), p = cor_p,
            file_name = paste0(model_name, ".Feature_select_COR"), dir_path = LR_var_dir_path)
            #lasso woe filter
            lasso_vars = lasso_filter(dat_train = train_woe, dat_test = test_woe, x_list = cor_select_vars,
            target = target, ex_cols = c(target, occur_time, obs_id, ex_cols), best_lambda = best_lambda,
            sim_sign = "negtive", plot.it = plot_show, save_data = FALSE, seed = seed,
            parallel = parallel, dir_path = LR_var_dir_path)

            save_dt(lasso_vars, as_list = TRUE, row_names = FALSE,
            file_name = paste0(model_name, ".Feature_select_LASSO"), dir_path = LR_var_dir_path)
         
            Formula = as.formula(paste(target, paste(unique(c(forced_in, lasso_vars)), collapse = ' + '), sep = ' ~ '))
            if (!is.null(seed)) set.seed(seed) else set.seed(46)
            lr_model = glm(Formula, data = train_woe[, c(target, unique(c(forced_in, lasso_vars)))], family = binomial(logit), weights = obsweighted)
            dt_coef = data.frame(summary(lr_model)$coefficients)
            lg_coef = subset(dt_coef, abs(dt_coef$Estimate) > 0)
            glm_vars = row.names(lg_coef)[-1]
            Formula = as.formula(paste(target, paste(unique(c(forced_in, glm_vars)), collapse = ' + '), sep = ' ~ '))
            lr_model_new = glm(Formula, data = train_woe[, c(target, unique(c(forced_in, glm_vars)))], family = binomial(logit), weights = obsweighted)
            if (step_wise) {
                lr_model_step = stats::step(lr_model_new, dir_pathection = "both", trace = TRUE)
                dt_step_coef = data.frame(summary(lr_model_step)$coefficients)
                step_vars = row.names(dt_step_coef)[-1]
                Formula = as.formula(paste(target, paste(unique(c(forced_in, step_vars)), collapse = ' + '), sep = ' ~ '))
                lr_model_new = glm(Formula, data = train_woe[, c(target, unique(c(forced_in, step_vars)))], family = binomial(logit), weights = obsweighted)
            }

            dt_imp_LR = get_logistic_coef(lg_model = lr_model_new, save_data = FALSE)
            lr_vars = dt_imp_LR[-1, "Feature"]
            save_dt(lr_vars, as_list = TRUE, file_name = paste0(model_name, ".Model_Features"), dir_path = LR_var_dir_path, note = TRUE)
            LR_iv_psi <- subset(iv_psi_list, iv_psi_list$Feature %in% lr_vars)[1:3]
            LR_iv_psi <- rbind(c("(Intercept)", 0, 0), LR_iv_psi)
            dt_imp_LR = merge(dt_imp_LR, LR_iv_psi)
            imp_vars = dt_imp_LR[-1, "Feature"]
            save_dt(dt_imp_LR, file_name = paste0(model_name, ".Model_Coef"), dir_path = LR_perf_dir_path, note = TRUE)
    
            if (length(imp_vars) > 1) {
                cor_plot(dat = train_woe, x_list = imp_vars, dir_path = LR_perf_dir_path, gtitle = paste0(model_name, ".LR"))
            }
            bins_table = get_bins_table_all(dat = dat_train, target = target,
              ex_cols = c(target, occur_time, obs_id, ex_cols), x_list = lr_vars, breaks_list = breaks_list,
            note = FALSE, save_data = TRUE, file_name = paste0(model_name, ".Feature_Bins"), dir_path = LR_perf_dir_path)
            if (score_card) {
                # standard socre card
                LR_score_card <- get_score_card(lg_model = lr_model_new, bins_table, target = target, dir_path = LR_perf_dir_path)
                #socre transforming
                train_pred = dat_train[, c(obs_id, occur_time, target)]
                test_pred = dat_test[, c(obs_id, occur_time, target)]
                train_pred$pred_LR = score_transfer(model = lr_model_new,  tbl_woe = train_woe, save_data = FALSE)[, "score"]
                test_pred$pred_LR = score_transfer(model = lr_model_new,  tbl_woe = test_woe, save_data = FALSE)[, "score"]
            } else {
                train_pred = dat_train[, c(obs_id, occur_time, target)]
                test_pred = dat_test[, c(obs_id, occur_time, target)]
                train_pred$pred_LR = round(predict(lr_model_new, train_woe[, imp_vars], type = "response"), 5)
                test_pred$pred_LR = round(predict(lr_model_new, test_woe[, imp_vars], type = "response"),5)
            }
            save_dt(train_pred, file_name = "train_pred_LR", dir_path = LR_pre_dir_path, note = FALSE)
            save_dt(test_pred, file_name = "test_pred_LR", dir_path = LR_pre_dir_path, note = FALSE)

            if (isTRUE(all.equal(names(train_pred), names(test_pred))) & all(lr_vars %in% names(dat_test))) {
                ks_psi_plot(train_pred = train_pred, test_pred = test_pred, score = "pred_LR",
                target = target, plot_show = plot_show, g_width = 10, dir_path = LR_perf_dir_path, gtitle = paste0(model_name, ".LR"))
                tb_pred = ks_table_plot(train_pred = train_pred, test_pred = test_pred, score = "pred_LR",
                target = target, g = 20, g_width = 13, plot_show = FALSE,
                dir_path = LR_perf_dir_path, gtitle = paste0(model_name, ".LR"))
                get_plots(dat_train = dat_train, dat_test = dat_test, x_list = lr_vars, target = target, ex_cols = ex_cols,
                breaks_list = breaks_list, pos_flag = pos_flag, occur_time = occur_time, parallel = parallel,
                g_width = 8,  dir_path = LR_perf_dir_path)
                key_index = model_key_index(tb_pred)
                params_key_index = data.frame(c(bins_params), key_index)
                save_dt(params_key_index, file_name = paste0(model_name, "LR_Params"), dir_path = LR_perf_dir_path,
                append = TRUE, note = FALSE)
            }
            if (save_pmml) {
                cat(paste("[NOTE] convert Logistic Regression model to pmml.\n"))
                model_pmml <- pmml(lr_model_new,
                model.name = "Logistic_Regression_Model",
                 description = "Logistic Regression Model",
                 copyright = NULL,
                transforms = NULL,
                unknownValue = NULL,
                 weights = NULL)
                saveXML(model_pmml, file = paste0(LR_model_dir_path,"/lr_model.pmml"))
            }

            save(lr_model_new, file = paste0(LR_model_dir_path, "/lg_model.RData")) #save model
            model_new = list(lr_model = lr_model_new)
        }

        if (any(algorithm == "XGB")) {
            cat("[NOTE] Training XGboost Model.\n")
            XGB_dir_path ="XGB"
            XGB_model_dir_path = paste0("./model/", XGB_dir_path)
            XGB_var_dir_path = paste0("./variable/", XGB_dir_path)
            XGB_perf_dir_path = paste0("./performance/", XGB_dir_path)
            XGB_pred_dir_path = paste0("./predict/", XGB_dir_path)
            XGB_data_dir_path = paste0("./data/", XGB_dir_path)
            if (!dir.exists(XGB_model_dir_path)) dir.create(XGB_model_dir_path)
            if (!dir.exists(XGB_var_dir_path)) dir.create(XGB_var_dir_path)
            if (!dir.exists(XGB_perf_dir_path)) dir.create(XGB_perf_dir_path)
            if (!dir.exists(XGB_pred_dir_path)) dir.create(XGB_pred_dir_path)
            if (!dir.exists(XGB_data_dir_path)) dir.create(XGB_data_dir_path)

            #get parameters
            if (is.null(XGB.params)) {
                XGB.params = xgb_params()
            }
            nrounds = ifelse(!is.null(XGB.params[["nrounds"]]), XGB.params[["nrounds"]], 5000)
            if (!is.null(XGB.params[["params"]])) {
                params = XGB.params[["params"]]
            } else {
                params = list(max.depth = 6, eta = 0.1,
                min_child_weight = 1, subsample = 1,
                colsample_bytree = 1, gamma = 0, max_delta_step = 0,
                eval_metric = "auc", objective = "binary:logistic" )
            }
            early_stopping_rounds = ifelse(!is.null(XGB.params[["early_stopping_rounds"]]), XGB.params[["early_stopping_rounds"]], 300)

            char_x_list = get_names(dat = dat_train[, x_list], types = c('character', 'factor'),
            ex_cols = c(target, obs_id, occur_time, ex_cols), get_ex = FALSE)

            if (length(char_x_list) > 0) {
                nr = nrow(dat_train)
                var_list = unique(c(target, obs_id, occur_time, x_list))
                dat_ts = rbind(dat_train[, var_list], dat_test[, var_list])
                dat_ts = one_hot_encoding(dat = dat_ts, cat_vars = char_x_list, na_act = FALSE)
                dat_ts = low_variance_filter(dat = dat_ts, lvp = 1, note = FALSE)
                dat_train = dat_ts[1:nr,]
                dat_test = dat_ts[-c(1:nr),]
                x_list = get_names(dat = dat_train, types = c('numeric', 'integer', 'double'),
                ex_cols = c(target, obs_id, occur_time, ex_cols), get_ex = FALSE)
                rm(dat_ts)
            }

            # Generate XGBoost feature map
            # mpg.fmap = genFMap(dat_train[, x_list])
            # Generate XGBoost DMatrix

            x_train = as.matrix(dat_train[, x_list])
            y_train = as.numeric(as.character(dat_train[, target]))
            xgb_train = list(data = x_train, label = y_train)
            dtrain = xgb.DMatrix(data = xgb_train$data, label = xgb_train$label)
            x_test = as.matrix(dat_test[, x_list])
            y_test = as.numeric(as.character(dat_test[, target]))
            xgb_test = list(data = x_test, label = y_test)
            dtest <- xgb.DMatrix(data = xgb_test$data, label = xgb_test$label)
            watchlist <- list(train = dtrain, eval = dtest)
            # Train a model
            if (!is.null(seed)) set.seed(seed) else set.seed(46)
            xgb_model_new = xgb.train(data = dtrain,
                                  watchlist = watchlist,
                             nrounds = nrounds,
                            params = params,
                             early_stopping_rounds = early_stopping_rounds,
                                  verbose = 1,
                                  maximize = TRUE)
            # feature importance
            dat_names <- dimnames(x_train)[[2]]
            imp_XGB <- xgb.importance(dat_names, model = xgb_model_new)
            imp_XGB = as.data.frame(imp_XGB)
            dt_imp_XGB = data.frame(Feature = imp_XGB[, "Feature"],
                             Importance = round(imp_XGB[, 'Gain'], 5),
                            stringsAsFactors = FALSE)

            save_dt(dt_imp_XGB, file_name = paste0(model_name, ".Featrue_imp_XGB"), dir_path = XGB_var_dir_path, note = FALSE)

            train_pred = dat_train[, c(obs_id, occur_time, target)]
            test_pred = dat_test[, c(obs_id, occur_time, target)]
            train_pred$pred_XGB = round(predict(xgb_model_new, as.matrix(dat_train[, x_list]), type = "response"), 5)
            test_pred$pred_XGB = round(predict(xgb_model_new, as.matrix(dat_test[, x_list]), type = "response"), 5)
            save_dt(train_pred, file_name ="train_pred_XGB", dir_path = XGB_pred_dir_path, note = FALSE)
            save_dt(test_pred, file_name = "test_pred_XGB", dir_path = XGB_pred_dir_path, note = FALSE)
            if (isTRUE(all.equal(names(dat_test), names(dat_train)))) {
                ks_psi_plot(train_pred = train_pred, test_pred = test_pred, gtitle = paste0(model_name, ".XGB"), score = "pred_XGB",
                plot_show = plot_show, target = target, dir_path = XGB_perf_dir_path)
                tb_pred = ks_table_plot(train_pred = train_pred, test_pred = test_pred, target = target, score = "pred_XGB",
                gtitle = paste0(model_name, ".XGB"), g = 20, plot_show = FALSE, dir_path = XGB_perf_dir_path)
                key_index = model_key_index(tb_pred)
                params_key_index = data.frame(c(params), key_index)
                save_dt(params_key_index, file_name = paste0(model_name, ".XGB_Params"), dir_path = XGB_perf_dir_path, append = TRUE, note = FALSE)
            }
            if (save_pmml) {
     
                cat(paste("[NOTE] convert XGboost model to pmml.\n"))
                # save the tree information file
                xgb.dump(xgb_model_new, paste0(XGB_model_dir_path, "/xgb_model.dumped.trees"))
                # Export the model to PMML.                 
                xgb_model_pmml = pmml(xgb_model_new,
                inputFeatureNames = x_list,
                outputLabelName = target,
                outputCategories = c(0, 1),
                xgbDumpFile = paste0(XGB_model_dir_path, "/xgb_model.dumped.trees"),
                model.name = "xgboost_Model",
                app.name = "R-PMML",
                description = "Extreme Gradient Boosting Model",
                copyright = NULL, transforms = NULL,
                unknownValue = NULL,
                parentInvalidValueTreatment = "returnInvalid",
                childInvalidValueTreatment = "asIs"
                )
                saveXML(xgb_model_pmml, file = paste0(XGB_model_dir_path, "/xgb_modgbmel.pmml"))
            }
            save(xgb_model_new, file = paste0(XGB_model_dir_path, "/xgb_model.RData"))#save model
            model_new = append(model_new, list(xgb_model = xgb_model_new), 1)
            rm(x_train, y_train, xgb_model_new )
        }
        if (any(algorithm == "GBDT")) {
            cat("[NOTE] Training GBDT Model.\n")

            GBDT_dir_path ="GBDT"
            GBDT_model_dir_path = paste0("./model/", GBDT_dir_path)
            GBDT_var_dir_path = paste0("./variable/", GBDT_dir_path)
            GBDT_perf_dir_path = paste0("./performance/", GBDT_dir_path)
            GBDT_pred_dir_path = paste0("./predict/", GBDT_dir_path)
            GBDT_data_dir_path = paste0("./data/", GBDT_dir_path)
            if (!dir.exists(GBDT_model_dir_path)) dir.create(GBDT_model_dir_path)
            if (!dir.exists(GBDT_var_dir_path)) dir.create(GBDT_var_dir_path)
            if (!dir.exists(GBDT_perf_dir_path)) dir.create(GBDT_perf_dir_path)
            if (!dir.exists(GBDT_pred_dir_path)) dir.create(GBDT_pred_dir_path)
            if (!dir.exists(GBDT_data_dir_path)) dir.create(GBDT_data_dir_path)
            if (is.null(GBM.params)) {
                GBM.params = gbm_params()
            }
            n.trees = ifelse(!is.null(GBM.params[["n.trees"]]), GBM.params[["n.trees"]], 100)
            interaction.depth = ifelse(!is.null(GBM.params[["interaction.depth"]]), GBM.params[["interaction.depth"]], 6)
            shrinkage = ifelse(!is.null(GBM.params[["shrinkage"]]), GBM.params[["shrinkage"]], 0.01)
            n.minobsinnode = ifelse(!is.null(GBM.params[["n.minobsinnode"]]), GBM.params[["n.minobsinnode"]], 30)
            bag.fraction = ifelse(!is.null(GBM.params[["bag.fraction"]]), GBM.params[["bag.fraction"]], 0.5)
            train.fraction = ifelse(!is.null(GBM.params[["train.fraction"]]), GBM.params[["train.fraction"]], 1)
            cv.folds = ifelse(!is.null(GBM.params[["cv.folds"]]), GBM.params[["cv.folds"]], 5)


            char_x_list = get_names(dat = dat_train[, x_list], types = c('character', 'factor'),
            ex_cols = c(target, obs_id, occur_time, ex_cols), get_ex = FALSE)
            if (length(char_x_list) > 0) {
                nr = nrow(dat_train)
                var_list = unique(c(target, obs_id, occur_time, x_list))
                dat_ts = rbind(dat_train[, var_list], dat_test[, var_list])
                dat_ts = one_hot_encoding(dat = dat_ts, cat_vars = char_x_list, na_act = FALSE)
                dat_train = dat_ts[1:nr,]
                dat_test = dat_ts[-c(1:nr),]
                x_list = get_x_list(x_list = NULL, dat_train = dat_train, dat_test = dat_test,
                ex_cols = c(target, obs_id, occur_time, ex_cols))
            }

            Formula = as.formula(paste(target, paste(x_list, collapse = ' + '), sep = ' ~ '))
            if (!is.null(seed)) set.seed(seed) else set.seed(46)
            gbdt_model_new = gbm(
                Formula,
                data = dat_train, # include variables for use only.
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
            # check performance using cross-validation.
            best.iter = gbm.perf(gbdt_model_new, method = "cv", plot.it = TRUE, oobag.curve = FALSE)
            imp_gbdt = as.data.frame(summary(gbdt_model_new, best.iter))
            dt_imp_GBDT = data.frame(Feature = imp_gbdt[, "var"],
                             Importance = round(imp_gbdt[, 'rel.inf'], 5),
                            stringsAsFactors = FALSE)
            imp_vars_gbdt = subset(dt_imp_GBDT, dt_imp_GBDT$Importance > 0)[, "Feature"]
            save_dt(dt_imp_GBDT, file_name = paste0(model_name, "Featrue_Imp_GBDT"), dir_path = GBDT_var_dir_path, note = FALSE)

            train_pred = dat_train[, c(obs_id, occur_time, target)]
            test_pred = dat_test[, c(obs_id, occur_time, target)]
            train_pred$pred_GBDT = round(predict(gbdt_model_new, dat_train[, x_list], best.iter, type = "response"), 5)
            test_pred$pred_GBDT = round(predict(gbdt_model_new, dat_test[, x_list], best.iter, type = "response"), 5)
            save_dt(train_pred, file_name = "train_pred_GBDT", dir_path = GBDT_pred_dir_path, note = FALSE)
            save_dt(test_pred, file_name = "test_pred_GBDT", dir_path = GBDT_pred_dir_path, note = FALSE)

            if (isTRUE(all.equal(names(dat_test), names(dat_train)))) {
                ks_psi_plot(train_pred = train_pred, test_pred = test_pred, gtitle = GBDT_dir_path,
                target = target, score = "pred_GBDT", plot_show = plot_show, dir_path = GBDT_perf_dir_path)
                tb_pred = ks_table_plot(train_pred = train_pred, test_pred = test_pred, target = target,
                score = "pred_GBDT", gtitle = GBDT_dir_path, g = 20, plot_show = FALSE, dir_path = GBDT_perf_dir_path)
                key_index = model_key_index(tb_pred)
                params_key_index = data.frame(c(GBM.params), key_index)
                save_dt(params_key_index, file_name = "params_GBDT", dir_path = GBDT_perf_dir_path, append = TRUE, note = FALSE)
            }
            if (save_pmml) {
                cat(paste("[NOTE] convert GBDT model to pmml.\n"))

                gbdt_model_pmml = pmml(gbdt_model_new,
                model.name = "GBDT_Model",
                app.name = "R-PMML",
                description = "Gradient Boosting Decision Tree Model",
                copyright = NULL,
                transforms = NULL,
                unknownValue = NULL,
                parentInvalidValueTreatment = "returnInvalid",
                childInvalidValueTreatment = "asIs"
                )
                saveXML(gbdt_model_pmml, file = paste0(GBDT_model_dir_path, "/gbdt_model.pmml"))
            }
            save(gbdt_model_new, file = paste0(GBDT_model_dir_path, "/gbdt_model.RData")) #save model     
            model_new = append(model_new, list(gbdt_model = gbdt_model_new), 1)
            rm(gbdt_model_new)
        }

        if (any(algorithm == "RF")) {
            cat("[NOTE] RandomForest model is building.\n")
       
            RF_dir_path = "RF"
            RF_model_dir_path = paste0("./model/", RF_dir_path)
            RF_var_dir_path = paste0("./variable/", RF_dir_path)
            RF_perf_dir_path = paste0("./performance/", RF_dir_path)
            RF_pred_dir_path = paste0("./predict/", RF_dir_path)
            RF_data_dir_path = paste0("./data/", RF_dir_path)
            if (!dir.exists(RF_model_dir_path)) dir.create(RF_model_dir_path)
            if (!dir.exists(RF_var_dir_path)) dir.create(RF_var_dir_path)
            if (!dir.exists(RF_perf_dir_path)) dir.create(RF_perf_dir_path)
            if (!dir.exists(RF_pred_dir_path)) dir.create(RF_pred_dir_path)
            if (!dir.exists(RF_data_dir_path)) dir.create(RF_data_dir_path) 
            `%DO%` = if (parallel) `%dopar%` else `%do%`
            if (is.null(RF.params)) {
               RF.params = rf_params()
            }
            ntree = ifelse(!is.null(RF.params[["ntree"]]), RF.params[["ntree"]], 100)
            nodesize = ifelse(!is.null(RF.params[["nodesize"]]), RF.params[["nodesize"]], 30)
            samp_rate = ifelse(!is.null(RF.params[["samp_rate"]]), RF.params[["samp_rate"]], 0.1)
            tune_rf = ifelse(!is.null(RF.params[["tune_rf"]]), RF.params[["tune_rf"]], TRUE)

            char_x_list = get_names(dat = dat_train[, x_list], types = c('character', 'factor'),
            ex_cols = c(target, obs_id, occur_time, ex_cols), get_ex = FALSE)
            if (length(char_x_list) > 0) {
                nr = nrow(dat_train)
                var_list = unique(c(target, obs_id, occur_time, x_list))
                dat_ts = rbind(dat_train[, var_list], dat_test[, var_list])
                dat_ts = one_hot_encoding(dat = dat_ts, cat_vars = char_x_list, na_act = FALSE)
                dat_train = dat_ts[1:nr,]
                dat_test = dat_ts[-c(1:nr),]
                x_list = get_names(dat = dat_train, types = c('numeric', 'integer', 'double'),
                ex_cols = c(target, obs_id, occur_time, ex_cols), get_ex = FALSE)
            }
           
            dat_train[, target] = as.factor(as.character(dat_train[, target]))
            if (!is.null(seed)) set.seed(seed) else set.seed(46)
            if (tune_rf) {
                #Tune Random Forest Model
                tRF = foreach(n_tree = rep(round(10 / cores_num), cores_num),
                .combine = randomForest::combine, .packages = "randomForest") %DO% {
                    tuneRF(
                         x = as.matrix(dat_train[, x_list]),
                         y = dat_train[, target],
                         stepFactor = 0.5,
                         plot = FALSE,
                         ntreeTry = n_tree,
                         trace = TRUE,
                         improve = 0.05,
                         doBest = TRUE,
                         strata = dat_train[, target],
                          sampsize = c(1, 1) * nrow(dat_train) * samp_rate,
                           nodesize = nodesize,
                           importance = TRUE,
                           proximity = FALSE
                         )
                }

                n_tree = tRF$ntree
                mtry = tRF$mtry
                imp_rf = as.data.frame(randomForest::importance(tRF))
               vars_rf = rownames(imp_rf[which(imp_rf$MeanDecreaseAccuracy > 0),])          
            } else {
                n_tree = ntree
                mtry = floor(sqrt(length(x_list)))
                vars_rf = x_list
            }

            #Fit the Random Forest Model After Tuning 
            Formula = as.formula(paste(target, paste(vars_rf, collapse = ' + '), sep = ' ~ '))
            if (!is.null(seed)) set.seed(seed) else set.seed(46)
            rf_model_new <- foreach(n_tree = rep(round(ntree / cores_num), cores_num),
            .combine = randomForest::combine, .packages = "randomForest") %DO% {
                randomForest(
                        Formula,
                           data = dat_train,
                          strata = dat_train[, target],
                           sampsize = c(1, 1) * nrow(dat_train) * samp_rate,
                           ntree = n_tree,
                    nodesize = nodesize,
                           importance = TRUE,
                           proximity = TRUE,
                           mtry = mtry)
            }
            imp_rf = randomForest::importance(rf_model_new)
            dt_imp_RF = data.frame(Feature = row.names(imp_rf),
                             Importance = round(imp_rf[, 'MeanDecreaseAccuracy'], 5),
                            stringsAsFactors = FALSE)
            imp_vars_rf = subset(dt_imp_RF, dt_imp_RF$Importance > 0)[, "Feature"]

            save_dt(dt_imp_RF, file_name = paste0(model_name, "Feature_Imp_RF"), dir_path = RF_var_dir_path, note = FALSE)

            train_pred = dat_train[, c(obs_id, occur_time, target)]
            test_pred = dat_test[, c(obs_id, occur_time, target)]
            train_pred$pred_RF= round(predict(rf_model_new, dat_train[, vars_rf], type = c("prob"))[, 2], 5)
            test_pred$pred_RF = round(predict(rf_model_new, dat_test[, vars_rf], type = c("prob"))[, 2], 5)
            save_dt(train_pred, file_name = "train_pred_RF", dir_path = RF_pred_dir_path, note = FALSE)
            save_dt(test_pred, file_name = "test_pred_RF", dir_path = RF_pred_dir_path, note = FALSE)


            if (isTRUE(all.equal(names(dat_test), names(dat_train)))) {
                ks_psi_plot(train_pred = train_pred, test_pred = test_pred, gtitle = RF_dir_path, target = target,
                score = "pred_RF", plot_show = plot_show, dir_path = RF_perf_dir_path)
                tb_pred = ks_table_plot(train_pred = train_pred, test_pred = test_pred, target = target,
                score = "pred_RF", gtitle = RF_dir_path, g = 20, plot_show = FALSE, dir_path = RF_perf_dir_path)
                key_index = model_key_index(tb_pred)
                params_key_index = data.frame(c(RF.params), key_index)
                save_dt(params_key_index, file_name = "params_RF", dir_path = RF_perf_dir_path, append = TRUE, note = FALSE)
            }
            if (save_pmml) {
                cat(paste("[NOTE] convert Random Forest model to pmml.\n"))

                rf_model_pmml = pmml(tRF,
                model.name = "randomForest_Model",
                app.name = "R-PMML",
                description = "Random Forest Tree Model",
                copyright = NULL,
                transforms = NULL,
                unknownValue = NULL,
                parentInvalidValueTreatment = "returnInvalid",
                childInvalidValueTreatment = "asIs")
                saveXML(rf_model_pmml, file = paste0(RF_model_dir_path, "/rf_model.pmml"))
            }
            save(rf_model_new, file = paste0(RF_model_dir_path, "/rf_model.RData")) #save model
            model_new = append(model_new, list(rf_model = rf_model_new), 1)
            rm(rf_model_new)
        }
    } else {
        stop(paste("target must be binomial.\n"))
    }
    return(model_new)
    options(opt) # reset
}

#'Logistic Regression & Scorecard Parameters
#'
#' \code{lr_params} is the list of parameters to train a LR model or Scorecard using in  \code{\link{training_model}}.
#' @param tree_control the list of parameters to control cutting initial breaks by decision tree. See details at: \code{\link{get_tree_breaks}}
#' @param bins_control  the list of parameters to control merging initial breaks. See details at: \code{\link{select_best_breaks}},\code{\link{select_best_class}}
#' @param best_lambda  Metheds of best lanmbda stardards using to filter variables by LASSO.There are four methods: ("lambda.min", "lambda.1se", "lambda.05se" , "lambda.sim_sign") . Default is  "lambda.sim_sign". See details at: \code{\link{get_best_lambda}}
#' @param obsweight An optional vector of  'prior weights' to be used in the fitting process. Should be NULL or a numeric vector. If you oversample or cluster diffrent datasets to training the LR model, you need to set this parameter to ensure that the probability of logistic regression output is the same as that before oversampling or segmentation. e.g.:There are 10,000 good obs and 500 bad obs before oversampling or under-sampling, 5,000 good obs and 3,000 bad obs after oversampling. Then this parameter should be set to c(10000/5000, 500/3000). Default is NULL..
#' @param forced_in Names of forced input variables. Default is NULL.
#' @param sp_values  Vaules will be in separate bins.e.g. list(-1, "Unknown")  means that -1 & Unknown as special values.Default is NULL.
#' @param step_wise  Logical, stepwise method. Default is TRUE.
#' @param score_card  Logical, transfer woe to a standard scorecard. If TRUE, Output scorecard, and score prediction, otherwise output probability. Default is TRUE.
#' @param cor_p  The maximum threshold of correlation. Default: 0.8.
#' @param iv_i The minimum threshold of IV. 0.01 to 0.1 usually work. Default: 0.02
#' @param psi_i The maximum threshold of PSI. 0.1 to 0.3 usually work. Default: 0.1.
#' @param ... Other parameters
#' @return A list of parameters.
#' @seealso  \code{\link{training_model}}, \code{\link{xgb_params}}, \code{\link{gbm_params}}, \code{\link{rf_params}}
#' @export

lr_params = function(tree_control = list(p = 0.02, cp = 0.00000001, xval = 5, maxdepth = 10),
                     bins_control = list(bins_num = 10, bins_pct = 0.05, b_chi = 0.02, b_odds = 0.1,
                     b_psi = 0.03, b_gb = 0.15, mono = 0.2, gb_psi = 0.15, kc = 1),
                                           best_lambda = "lambda.sim_sign", sp_values = NULL,
                                           forced_in = NULL, obsweight = c(1, 1),
                                           step_wise = TRUE, score_card = TRUE,
                                           cor_p = 0.8, iv_i = 0.02, psi_i = 0.1, ...) {
    structure(list(tree_control = tree_control, bins_control = bins_control,
                                           best_lambda = best_lambda, sp_values = sp_values,
                                           forced_in = forced_in, obsweight = obsweight,
                                           step_wise = step_wise, score_card = score_card,
                                           cor_p = cor_p, iv_i = iv_i, psi_i = psi_i))
}

#'Logistic Regression & Scorecard Parameters
#'
#' \code{xgb_params} is the list of parameters to train a LR model or Scorecard using in \code{\link{training_model}}.
#' @param nrounds Max number of boosting iterations.
#' @param params  A list contains parameters of xgboost.The complete list of parameters is available at: \url{ http://xgboost.readthedocs.io/en/latest/parameter.html}
#' @param early_stopping_rounds  If NULL, the early stopping function is not triggered. If set to an integer k, training with a validation set will stop if the performance doesn't improve for k rounds. 
#' @param ... Other parameters
#' @return A list of parameters.
#' @seealso \code{\link{training_model}}, \code{\link{lr_params}},\code{\link{gbm_params}}, \code{\link{rf_params}}
#' @export

xgb_params = function(nrounds = 1000, params = list(max.depth = 6, eta = 0.1, min_child_weight = 1, subsample = 1,
colsample_bytree = 1, gamma = 0, max_delta_step = 0, eval_metric = "auc", objective = "binary:logistic"), early_stopping_rounds = 100, ...) {
    structure(list(nrounds = nrounds, params = params, early_stopping_rounds = early_stopping_rounds))
}

#' GBM Parameters
#'
#' \code{gbm_params} is the list of parameters to train a GBM using in  \code{\link{training_model}}.
#' @param n.trees Integer specifying the total number of trees to fit. This is equivalent to the number of iterations and the number of basis functions in the additive expansion. Default is 100.
#' @param interaction.depth Integer specifying the maximum depth of each tree(i.e., the highest level of variable interactions allowed) . A value of 1 implies an additive model, a value of 2 implies a model with up to 2 - way interactions, etc. Default is 1.
#' @param n.minobsinnode  Integer specifying the minimum number of observations in the terminal nodes of the trees. Note that this is the actual number of observations, not the total weight.
#' @param shrinkage a shrinkage parameter applied to each tree in the expansion. Also known as the learning rate or step - size reduction; 0.001 to 0.1 usually work, but a smaller learning rate typically requires more trees. Default is 0.1 .
#' @param bag.fraction  the fraction of the training set observations randomly selected to propose the next tree in the expansion. This introduces randomnesses into the model fit. If bag.fraction < 1 then running the same model twice will result in similar but different fits. gbm uses the R random number generator so set.seed can ensure that the model can be reconstructed. Preferably, the user can save the returned gbm.object using save. Default is 0.5 .
#' @param train.fraction The first train.fraction * nrows(data) observations are used to fit the gbm and the remainder are used for computing out-of-sample estimates of the loss function.
#' @param cv.folds  Number of cross - validation folds to perform. If cv.folds > 1 then gbm, in addition to the usual fit, will perform a cross - validation, calculate an estimate of generalization error returned in cv.error.
#' @param ... Other parameters
#' @return A list of parameters.
#' @details See details at: \code{gbm}
#' @seealso \code{\link{training_model}}, \code{\link{lr_params}}, \code{\link{xgb_params}}, \code{\link{rf_params}}
#' @export


gbm_params = function(n.trees = 1000, interaction.depth = 6, shrinkage = 0.01,
bag.fraction = 0.5, train.fraction = 0.7, n.minobsinnode = 30, cv.folds = 5, ...) {
    structure(list(n.trees = n.trees, interaction.depth = interaction.depth, shrinkage = shrinkage,
    bag.fraction = bag.fraction, train.fraction = train.fraction, n.minobsinnode = n.minobsinnode, cv.folds = cv.folds))
}

#' Random Forest Parameters
#'
#' \code{rf_params} is the list of parameters to train a Random Forest using in  \code{\link{training_model}}.
#' @param ntree Number of trees to grow. This should not be set to too small a number, to ensure that every input row gets predicted at least a few times. 
#' @param nodesize  Minimum size of terminal nodes. Setting this number larger causes smaller trees to be grown (and thus take less time). Note that the default values are different for classification (1) and regression (5).
#' @param samp_rate   Percentage of sample to draw. Default is 0.2.
#' @param tune_rf A logical.If TRUE, then tune Random Forest model.Default is FALSE.
#' @param ... Other parameters
#' @return A list of parameters.
#' @details See details at : \url{https://www.stat.berkeley.edu/~breiman/Using_random_forests_V3.1.pdf}
#' @seealso  \code{\link{training_model}}, \code{\link{lr_params}}, \code{\link{gbm_params}}, \code{\link{xgb_params}}
#' @export
rf_params = function(ntree = 100, nodesize = 30, samp_rate = 0.5,
tune_rf = FALSE, ...) {
    structure(list(ntree = ntree, nodesize = nodesize, samp_rate = samp_rate,
    tune_rf = tune_rf))
}
