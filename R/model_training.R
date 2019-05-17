#'Training model
#'
#' \code{training_model} Model builder
#' @param model_name  A string, name of the project. Default is "mymodel"
#' @param model_path The path for periodically saved data file. Default is \code{tempdir()}.
#' @param dat A data.frame with independent variables and target variable.
#' @param dat_test  A data.frame of test data. Default is NULL.
#' @param target The name of target variable.
#' @param x_list Names of independent variables. Default is NULL.
#' @param occur_time The name of the variable that represents the time at which each observation takes place.Default is NULL.
#' @param pos_flag The value of positive class of target variable, default: "1".
#' @param ex_cols Names of excluded variables. Regular expressions can also be used to match variable names. Default is NULL.
#' @param prop Percentage of train-data after the partition. Default: 0.7.
#' @param obs_id  The name of ID of observations or key variable of data. Default is NULL.
#' @param miss_values  Other extreme value might be used to represent missing values, e.g: -9999, -9998. These miss_values will be encoded to -1 or "Missing".
#' @param preproc Logical. Preprocess data. Default is TRUE
#' @param outlier_proc Logical. If TRUE,  Outliers processing using Kmeans and Local Outlier Factor. Default is TRUE
#' @param missing_proc Logical. If TRUE, missing value analysis and process missing value by knn imputation or central impulation or random imputation. Default is TRUE
#' @param one_hot  Logical. If TRUE, one-hot_encoding  of category variables. Default is FASLE.
#' @param feature_filter  Parameters for selecting important and stable features.See details at: \code{\link{feature_select_wrapper}}
#' @param algorithm  Algorithms for training a model. list("LR", "XGB", "GBDT", "RF") are available.
#' @param LR.params  Parameters of logistic regression & scorecard. See details at :  \code{\link{lr_params}}.
#' \itemize{
#'   \item \code{tree_control} the list of parameters to control cutting initial breaks by decision tree. See details at: \code{\link{get_tree_breaks}}
#'   \item \code{bins_control} the list of parameters to control merging initial breaks. See details at: \code{\link{select_best_breaks}},\code{\link{select_best_class}}
#'   \item \code{best_lambda} Metheds of best lanmbda stardards using to filter variables by LASSO. There are 3 methods: ("lambda.auc", "lambda.ks", "lambda.sim_sign") . Default is  "lambda.sim_sign".
#'   \item \code{obsweight} An optional vector of  'prior weights' to be used in the fitting process. Should be NULL or a numeric vector. If you oversample or cluster diffrent datasets to training the LR model, you need to set this parameter to ensure that the probability of logistic regression output is the same as that before oversampling or segmentation. e.g.:There are 10,000 0 obs and 500 1 obs before oversampling or under-sampling, 5,000 0 obs and 3,000 1 obs after oversampling. Then this parameter should be set to c(10000/5000, 500/3000). Default is NULL..
#'   \item \code{forced_in}Names of forced input variables. Default is NULL.
#'   \item \code{sp_values} Vaules will be in separate bins.e.g. list(-1, "Missing")  means that -1 & Missing as special values.Default is NULL.
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
#' @seealso   \code{\link{train_test_split}},\code{\link{data_cleansing}}, \code{\link{feature_select_wrapper}},   \code{\link{lr_params}}, \code{\link{xgb_params}}, \code{\link{gbm_params}}, \code{\link{rf_params}},\code{\link{fast_high_cor_filter}},\code{\link{get_breaks_all}},\code{\link{lasso_filter}}, \code{\link{woe_trans_all}}, \code{\link{get_logistic_coef}}, \code{\link{score_transfer}},\code{\link{get_score_card}}, \code{\link{model_key_index}},\code{\link{ks_psi_plot}},\code{\link{get_plots}},\code{\link{ks_table_plot}}
#' @examples
#' sub = cv_split(UCICreditCard, k = 30)[[1]]
#' dat = UCICreditCard[sub,]
#' dat = re_name(dat, "default.payment.next.month", "target")
#' dat = data_cleansing(dat, target = "target", obs_id = "ID", 
#' occur_time = "apply_date", miss_values = list("", -1, -2))
#' train_test <- train_test_split(dat, split_type = "OOT", prop = 0.7,
#'                                 occur_time = "apply_date")
#' dat_train = train_test$train
#' dat_test = train_test$test
#' x_list = c("PAY_0", "LIMIT_BAL", "PAY_AMT5", "PAY_3", "PAY_2")
#' B_model = training_model(dat = dat_train,
#' model_name = "UCICreditCard", target = "target", x_list = x_list,
#' occur_time = "apply_date", obs_id = "ID", dat_test = dat_test,
#'                            preproc = FALSE,
#'                            feature_filter = NULL,
#'                            algorithm = list("LR"),
#'                            LR.params = lr_params(lasso = FALSE, 
#'                            step_wise = FALSE, vars_plot = FALSE),
#'                            XGB.params = xgb_params(),
#'                            breaks_list = NULL,
#'                            parallel = FALSE, cores_num = NULL,
#'                            save_pmml = FALSE, plot_show = FALSE,
#'                            model_path = tempdir(),
#'                            seed = 46)
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



training_model <- function(model_name = "mymodel", dat, dat_test = NULL,
                           target = NULL,occur_time = NULL, obs_id = NULL,
                           x_list = NULL, ex_cols = NULL,
                           pos_flag = NULL, prop = 0.7,
                           preproc = TRUE, one_hot = FALSE, miss_values = NULL,
                           outlier_proc = TRUE, missing_proc = TRUE, 
                           feature_filter = list(filter = c("IV", "PSI", "COR", "XGB"),
                                                 iv_cp = 0.02, psi_cp = 0.1, xgb_cp = 0,
                                                 cv_folds = 1, hopper = FALSE),
                           algorithm = list("LR", "XGB","GBM","RF"),
                           LR.params = lr_params(),
                           XGB.params = xgb_params(),
                           GBM.params = gbm_params(),
                           RF.params = rf_params(),
                           breaks_list = NULL,
                           parallel = FALSE,cores_num = NULL,
                           save_pmml = FALSE, plot_show = FALSE,
                           model_path = getwd(),
                           seed = 46, ...) {

    opt = options(stringsAsFactors = FALSE)
    if (is.null(algorithm)) {
        stop(paste("algorithm is missing.\n"))
    }
    if (!any(is.element(algorithm, c("LR", "XGB", "GBM", "RF")))) {
        stop("In algorithm, only  LR, XGB, GBM, RF are supported.\n")
    }
    if (length(x_list) > 0 && any(x_list == target)) {
        stop(paste("x_list  contains", target, ".\n"))
    }
    if (!dir.exists(model_path)) dir.create(model_path)
    if (!is.character(model_name)) model_name = "my_model"
    model_path = ifelse(!is.character(model_path) ,
                      paste(tempdir(), model_name, sep = "/"), paste(model_path, model_name, sep = "/"))
    if (!dir.exists(model_path)) dir.create(model_path)
    model_dir_path = paste(model_path, "model", sep = "/")
    data_dir_path = paste(model_path, "data", sep = "/")
    var_dir_path = paste(model_path, "variable", sep = "/")
    perf_dir_path = paste(model_path, "performance", sep = "/")
   pred_dir_path = paste(model_path, "predict", sep = "/")
    if (!dir.exists(model_dir_path)) dir.create(model_dir_path)
    if (!dir.exists(data_dir_path)) dir.create(data_dir_path)
    if (!dir.exists(var_dir_path)) dir.create(var_dir_path)
    if (!dir.exists(perf_dir_path)) dir.create(perf_dir_path)
    if (!dir.exists(pred_dir_path)) dir.create(pred_dir_path)
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

    
    if (!is.null(dat_test)) {
        dat = checking_data(dat = dat, target = target, pos_flag = pos_flag)
        dat_test = checking_data(dat = dat_test, target = target, pos_flag = pos_flag)
        x_list = get_x_list(x_list = x_list,
                        dat_train = dat, dat_test = dat_test,
                        ex_cols = c(target, obs_id, occur_time, ex_cols))
        nr = nrow(dat)
        com_list = unique(c(obs_id, occur_time, target, x_list))
        dat = dat[, com_list]
        dat_test = dat_test[, com_list]
        dat = rbind(dat, dat_test)
    } else {
        dat = checking_data(dat = dat, target = target, pos_flag = pos_flag)
        x_list = get_x_list(x_list = x_list,
                        dat_train = dat, dat_test = dat_test,
                        ex_cols = c(target, obs_id, occur_time, ex_cols))
    }

    #Data cleansing & preperation
    if (preproc) {
        dat = data_cleansing(dat = dat, x_list = x_list,
                                  target = target, obs_id = obs_id, occur_time = occur_time,
                                  ex_cols = ex_cols,
                                  outlier_proc = outlier_proc,
                                  missing_proc = missing_proc,
                                  miss_values = miss_values, one_hot = one_hot,
                                  low_var = TRUE, parallel = parallel, note = TRUE,
                                  save_data = TRUE, dir_path = data_dir_path, file_name = model_name)

        x_list = get_x_list(x_list = x_list,
                            dat_train = dat, ex_cols = c(target, obs_id, occur_time, ex_cols))
    }
    #train test spliting
    if (is.null(dat_test)) {
        train_test <- train_test_split(dat, split_type = "OOT",
                                       prop = prop, occur_time = occur_time, note = TRUE,
                                      save_data = TRUE, dir_path = data_dir_path,
                                      file_name = model_name,
                                      seed = seed)
        dat_train = train_test$train
        dat_test = train_test$test
    } else {
        train_test = train_test_split(dat, split_type = "byRow", prop = nr / nrow(dat),
        occur_time = occur_time, seed = seed, note = FALSE,
        save_data = TRUE, dir_path = data_dir_path, file_name = model_name)
        dat_train = train_test$train
        dat_test = train_test$test
    }
  
    if (!is.null(feature_filter)) {
        sel_list = NULL
        if (!is.null(feature_filter[["filter"]]) &&
            any(is.element(feature_filter[["filter"]], c("IV", "PSI", "XGB","COR")))) {
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
        sel_list = feature_select_wrapper(dat_train = dat_train, dat_test = dat_test,
                                          target = target,
                                          x_list = x_list, occur_time = occur_time,
                                          ex_cols = c(obs_id, occur_time, target, ex_cols),
                                          filter = filter, cv_folds = cv_folds,
                                          iv_cp = iv_cp, psi_cp = psi_cp,
                                          cor_cp = cor_cp, xgb_cp = xgb_cp, hopper = hopper,
                                          vars_name = TRUE, save_data = TRUE, parallel = parallel,
                                          cores_num = cores_num, seed = seed,
                                          file_name = model_name,note = TRUE,
                                          dir_path = var_dir_path)
        if (length(sel_list) > 0) {
            x_list = get_x_list(x_list = sel_list,
                                dat_train = dat_train, dat_test = dat_test,
                                ex_cols = c(target, obs_id, occur_time, ex_cols))
        } else {
            warning("No feature satisfies the criteria for feature selection, use the previous x_list.\n")
        }
    }
    model_new = lr_model_new = xgb_model_new = gbdt_model_new = rf_model_new = NULL
    # Train Models
    if (length(unique(dat_train[, target])) == 2) {
        if (any(algorithm == "LR")) {
            cat("[NOTE] Training Logistic Regression Model.\n")

            LR_model_dir_path = paste(model_dir_path, "LR", sep = "/")
            LR_var_dir_path = paste(var_dir_path, "LR", sep = "/")
            LR_perf_dir_path = paste(perf_dir_path, "LR", sep = "/")
            LR_pred_dir_path = paste(pred_dir_path, "LR", sep = "/")
            LR_data_dir_path = paste(data_dir_path, "LR", sep = "/")

            if (!dir.exists(LR_model_dir_path)) dir.create(LR_model_dir_path)
            if (!dir.exists(LR_var_dir_path)) dir.create(LR_var_dir_path)
            if (!dir.exists(LR_perf_dir_path)) dir.create(LR_perf_dir_path)
            if (!dir.exists(LR_pred_dir_path)) dir.create(LR_pred_dir_path)
            if (!dir.exists(LR_data_dir_path)) dir.create(LR_data_dir_path)
            if (is.null(LR.params)) {
                LR.params = lr_params()
            }
            obsweight = ifelse(!is.null(LR.params[["obsweight"]]), LR.params[["obsweight"]], 1)
            step_wise = ifelse(!is.null(LR.params[["step_wise"]]), LR.params[["step_wise"]], TRUE)
            lasso = ifelse(!is.null(LR.params[["lasso"]]), LR.params[["lasso"]], TRUE)
            vars_plot = ifelse(!is.null(LR.params[["vars_plot"]]), LR.params[["vars_plot"]], TRUE)
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
                list(bins_num = 10, bins_pct = 0.05,
                     b_chi = 0.02, b_odds = 0.1, b_psi = 0.03,
                     b_or = 0.15, mono = 0.2,odds_psi = 0.15,kc = 1)
            }
            if (length(obsweight) > 1 && is.vector(obsweight, mode = "numeric")) {
                obsweighted <- ifelse(dat_train[, target] == 0, obsweight[1], obsweight[2])
            } else {
                obsweighted = NULL
            }
           bins_params = bins_control

            #bins breaks

            if (is.null(breaks_list) || length(breaks_list) < 1) {
                breaks_list = get_breaks_all(dat = dat_train,
                                             x_list = x_list,
                                             ex_cols = c(obs_id, occur_time, target, ex_cols),
                                             occur_time = occur_time, oot_pct = prop,
                                             target = target,
                                             tree_control = tree_control,
                                             bins_control = bins_control,
                                             best = TRUE,
                                             sp_values = sp_values, parallel = parallel,
                                             note = TRUE,
                                             save_data = TRUE,
                                             file_name = paste0(model_name, ".Breaks_List"),
                                             dir_path = LR_var_dir_path)
            }

            #psi_iv_filter
            iv_psi_list = psi_iv_filter(dat = dat_train, dat_test = dat_test,
                                        x_list = x_list, target = target,
                                        occur_time = occur_time, oot_pct = prop,
                                        parallel = parallel,
                                        ex_cols = c(obs_id, occur_time, target, ex_cols),
                                        breaks_list = breaks_list,
                                        psi_i = psi_i, iv_i = iv_i,
                                        note = TRUE,
                                        save_data = TRUE,
                                        file_name = model_name,
                                        dir_path = LR_var_dir_path)
            if (length(iv_psi_list) < 1) {
                stop(paste("No variable satisfies the psi & iv condition."))
            }

            select_vars = as.character(iv_psi_list[, "Feature"])
            save_dt(select_vars, as_list = TRUE, row_names = FALSE, note = TRUE,
                    file_name = paste0(model_name, ".feature.IV_PSI"),
                    dir_path = LR_var_dir_path)

            #woe transform
            train_woe = woe_trans_all(dat = dat_train,
                                      x_list = unlist(select_vars),
                                      target = target,
                                      ex_cols = c(target, occur_time, obs_id),
                                      bins_table = NULL,
                                      breaks_list = breaks_list,
                                      woe_name = FALSE,note = TRUE, parallel = parallel,
                                      save_data = TRUE, file_name = paste0(model_name, ".train"),
                                      dir_path = LR_data_dir_path)
            test_woe = woe_trans_all(dat = dat_test,
                                     x_list = unlist(select_vars),
                                     target = target,
                                     ex_cols = c(target, occur_time, obs_id),
                                     bins_table = NULL,
                                     breaks_list = breaks_list,
                                     note = TRUE, woe_name = FALSE, parallel = parallel,
                                     save_data = TRUE, file_name = paste0(model_name, ".test"),
                                     dir_path = LR_data_dir_path)
            #high correlation filter
            select_vars = fast_high_cor_filter(dat = train_woe,
                                                   x_list = select_vars,
                                                   com_list = iv_psi_list,
                                                   ex_cols = c(target, occur_time, obs_id, ex_cols),
                                                   p = cor_p,
                                                   note = TRUE,
                                                   save_data = TRUE,
                                                   file_name = model_name,
                                                   dir_path = LR_var_dir_path)
            #lasso woe filter
            if (lasso) {
                select_vars = lasso_filter(dat_train = train_woe, dat_test = test_woe,
                                      x_list = select_vars,
                                      target = target,
                                      ex_cols = c(target, occur_time, obs_id, ex_cols),
                                      sim_sign ="negtive",
                                      best_lambda = best_lambda, 
                                      plot.it = plot_show, seed = seed, 
                                      save_data = TRUE, file_name = model_name, dir_path = LR_var_dir_path)
            }
            save_dt(select_vars, as_list = TRUE, row_names = FALSE, note = TRUE,
                    file_name = paste0(model_name, ".feature.premodel"), dir_path = LR_var_dir_path)
            #training lr model
            Formula = as.formula(paste(target, paste(unique(c(forced_in, select_vars)), collapse = ' + '),
                                       sep = ' ~ '))
            if (!is.null(seed)) set.seed(seed) else set.seed(46)
            lr_model = glm(Formula,
                           data = train_woe[, c(target, unique(c(forced_in, select_vars)))],
                           family = binomial(logit),
                           weights = obsweighted)
            dt_coef = data.frame(summary(lr_model)$coefficients)
            lg_coef = subset(dt_coef, abs(dt_coef$Estimate) > 0)
            glm_vars = row.names(lg_coef)[-1]
            Formula = as.formula(paste(target, paste(unique(c(forced_in, glm_vars)), collapse = ' + '), sep = ' ~ '))
            lr_model_new = glm(Formula,
                               data = train_woe[, c(target, unique(c(forced_in, glm_vars)))],
                               family = binomial(logit),
                               weights = obsweighted)
            #step wise
            if (step_wise) {
                lr_model_step = stats::step(lr_model_new, dir_pathection = "both", trace = TRUE)
                dt_step_coef = data.frame(summary(lr_model_step)$coefficients)
                step_vars = row.names(dt_step_coef)[-1]
                Formula = as.formula(paste(target, paste(unique(c(forced_in, step_vars)), collapse = ' + '), sep = ' ~ '))
                lr_model_new = glm(Formula,
                                   data = train_woe[, c(target, unique(c(forced_in, step_vars)))],
                                   family = binomial(logit),
                                   weights = obsweighted)
            }
            # get lr coef
            dt_imp_LR = get_logistic_coef(lg_model = lr_model_new, file_name = model_name,
                                                dir_path = LR_perf_dir_path, save_data = TRUE)
            lr_vars = dt_imp_LR[-1, "Feature"]
            save_dt(lr_vars, as_list = TRUE,
                    file_name = paste0(model_name, ".feature.model"),
                    dir_path = LR_var_dir_path, note = TRUE)
            LR_iv_psi <- subset(iv_psi_list, iv_psi_list$Feature %in% lr_vars)[1:3]
            LR_iv_psi <- rbind(c("(Intercept)", 0, 0), LR_iv_psi)
            dt_imp_LR = merge(dt_imp_LR, LR_iv_psi)
            imp_vars = dt_imp_LR[-1, "Feature"]
            save_dt(dt_imp_LR,
                    file_name =model_name,
                    dir_path = LR_perf_dir_path, note = TRUE)
            #correlation matrix plot of input variables
            if (vars_plot & length(imp_vars) > 1) {
                cor_plot(dat = train_woe, x_list = imp_vars,
                         dir_path = LR_perf_dir_path,save_data = TRUE,
                         gtitle = paste0(model_name, ".LR"))
                get_plots(dat_train = dat_train, dat_test = dat_test,
                          x_list = lr_vars, target = target,
                          ex_cols = ex_cols,
                          breaks_list = breaks_list,
                          save_data = TRUE,
                          occur_time = occur_time, parallel = parallel, plot_show = plot_show,
                          g_width = 8, dir_path = LR_perf_dir_path)
            }
            bins_table = get_bins_table_all(dat = dat_train, target = target,
                                            ex_cols = c(target, occur_time, obs_id, ex_cols),
                                            x_list = lr_vars, breaks_list = breaks_list,
                                            occur_time = occur_time,
                                            oot_pct = prop,
                                            note = TRUE, save_data = TRUE,
                                            file_name =model_name,
                                            dir_path = LR_perf_dir_path)
            if (score_card) {
                # standard socre card
                LR_score_card <- get_score_card(lg_model = lr_model_new,
                                                bins_table, target = target,
                                                file_name = model_name,
                                                dir_path = LR_perf_dir_path,
                save_data = TRUE)
                #socre transforming
                train_pred = dat_train[ c(obs_id, occur_time, target)]
                test_pred = dat_test[c(obs_id, occur_time, target)]
                train_pred$pred_LR = score_transfer(model = lr_model_new,
                                                    tbl_woe = train_woe,
                                                    save_data = TRUE)[, "score"]
                test_pred$pred_LR = score_transfer(model = lr_model_new,
                                                   tbl_woe = test_woe,                                               
                                                   save_data = TRUE)[, "score"]
            } else {
                train_pred = dat_train[c(obs_id, occur_time, target)]
                test_pred = dat_test[c(obs_id, occur_time, target)]
                train_pred$pred_LR = round(predict(lr_model_new, train_woe[, imp_vars], type = "response"), 5)
                test_pred$pred_LR = round(predict(lr_model_new, test_woe[, imp_vars], type = "response"),5)
            }
            save_dt(train_pred, file_name =paste(model_name, "LR.pred.train", sep = "."),
                    dir_path = LR_pred_dir_path, note = TRUE)
            save_dt(test_pred, file_name = paste(model_name, "LR.pred.test", sep = "."),
                    dir_path = LR_pred_dir_path, note = TRUE)
            #plot the model results
            if (isTRUE(all.equal(names(train_pred), names(test_pred))) &
                all(lr_vars %in% names(dat_test))) {
                ks_psi_plot(train_pred = train_pred, test_pred = test_pred,
                            score = "pred_LR",target = target,
                            plot_show = plot_show, g_width = 10,
                           save_data = TRUE,
                            dir_path = LR_perf_dir_path, gtitle = paste0(model_name, ".LR"))
                tb_pred = ks_table_plot(train_pred = train_pred, test_pred = test_pred,
                                        score = "pred_LR", target = target,
                                        g = 20, g_width = 13, plot_show = FALSE, save_data = TRUE,
                                        dir_path = LR_perf_dir_path,
                                        gtitle = paste0(model_name, ".LR"))         
                key_index = model_key_index(tb_pred)
                params_key_index = data.frame(c(bins_params), key_index)
                save_dt(params_key_index,
                        file_name = paste0(model_name, "LR_Params"),
                        dir_path = LR_perf_dir_path,
                append = TRUE, note = TRUE)
            }
            #transfer to pmml
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
            XGB_model_dir_path = paste(model_dir_path, "XGB", sep = "/")
            XGB_var_dir_path = paste(var_dir_path, "XGB", sep = "/")
            XGB_perf_dir_path = paste(perf_dir_path, "XGB", sep = "/")
            XGB_pred_dir_path = paste(pred_dir_path, "XGB", sep = "/")
            XGB_data_dir_path = paste(data_dir_path, "XGB", sep = "/")
            if (!dir.exists(XGB_model_dir_path)) dir.create(XGB_model_dir_path)
            if (!dir.exists(XGB_var_dir_path)) dir.create(XGB_var_dir_path)
            if (!dir.exists(XGB_perf_dir_path)) dir.create(XGB_perf_dir_path)
            if (!dir.exists(XGB_pred_dir_path)) dir.create(XGB_pred_dir_path)
            if (!dir.exists(XGB_data_dir_path)) dir.create(XGB_data_dir_path)

            #get parameters
            if (is.null(XGB.params)) {
                XGB.params = xgb_params()
            }
            nrounds = ifelse(!is.null(XGB.params[["nrounds"]]),
                             XGB.params[["nrounds"]], 5000)
            if (!is.null(XGB.params[["params"]])) {
                params = XGB.params[["params"]]
            } else {
                params = list(max.depth = 6, eta = 0.1,
                min_child_weight = 1, subsample = 1,
                colsample_bytree = 1, gamma = 0, max_delta_step = 0,
                eval_metric = "auc", objective = "binary:logistic" )
            }
            early_stopping_rounds = ifelse(!is.null(XGB.params[["early_stopping_rounds"]]),
                                           XGB.params[["early_stopping_rounds"]], 300)

            char_x_list = get_names(dat = dat_train[, x_list],
                                    types = c('character', 'factor'),
                                    ex_cols = c(target, obs_id, occur_time, ex_cols),
                                    get_ex = FALSE)

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

            save_dt(dt_imp_XGB, file_name = paste0(model_name, ".Featrue_imp_XGB"),
                    dir_path = XGB_var_dir_path, note = FALSE)

            train_pred = dat_train[c(obs_id, occur_time, target)]
            test_pred = dat_test[ c(obs_id, occur_time, target)]
            train_pred$pred_XGB = round(predict(xgb_model_new,
                                                x_train, type = "response"), 5)
            test_pred$pred_XGB = round(predict(xgb_model_new,
                                               x_test, type = "response"), 5)
            save_dt(train_pred, file_name ="train_pred_XGB", dir_path = XGB_pred_dir_path, note = FALSE)
            save_dt(test_pred, file_name = "test_pred_XGB", dir_path = XGB_pred_dir_path, note = FALSE)
            if (isTRUE(all.equal(names(train_pred), names(test_pred)))) {
                ks_psi_plot(train_pred = train_pred, test_pred = test_pred,
                            score = "pred_XGB", target = target,
                            plot_show = plot_show, g_width = 10,
                           save_data = TRUE,
                            dir_path = XGB_perf_dir_path, gtitle = paste0(model_name, ".XGB"))
                tb_pred = ks_table_plot(train_pred = train_pred, test_pred = test_pred,
                                        score = "pred_XGB", target = target,
                                        g = 20, g_width = 13, plot_show = FALSE, save_data = TRUE,
                                        dir_path = XGB_perf_dir_path,
                                        gtitle = paste0(model_name, ".XGB"))
                key_index = model_key_index(tb_pred)
                params_key_index = data.frame(c(params), key_index)
                save_dt(params_key_index,
                        file_name = paste0(model_name, ".XGB_Params"),
                        dir_path = XGB_perf_dir_path, append = TRUE, note = TRUE)
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
        if (any(algorithm == "GBM")) {
            cat("[NOTE] Training GBM Model.\n")

            GBM_model_dir_path = paste(model_dir_path, "GBM", sep = "/")
            GBM_var_dir_path = paste(var_dir_path, "GBM", sep = "/")
            GBM_perf_dir_path = paste(perf_dir_path, "GBM", sep = "/")
            GBM_pred_dir_path = paste(pred_dir_path, "GBM", sep = "/")
            GBM_data_dir_path = paste(data_dir_path, "GBM", sep = "/")
            if (!dir.exists(GBM_model_dir_path)) dir.create(GBM_model_dir_path)
            if (!dir.exists(GBM_var_dir_path)) dir.create(GBM_var_dir_path)
            if (!dir.exists(GBM_perf_dir_path)) dir.create(GBM_perf_dir_path)
            if (!dir.exists(GBM_pred_dir_path)) dir.create(GBM_pred_dir_path)
            if (!dir.exists(GBM_data_dir_path)) dir.create(GBM_data_dir_path)
            if (is.null(GBM.params)) {
                GBM.params = gbm_params()
            }
            n.trees = ifelse(!is.null(GBM.params[["n.trees"]]),
                             GBM.params[["n.trees"]], 100)
            interaction.depth = ifelse(!is.null(GBM.params[["interaction.depth"]]),
                                       GBM.params[["interaction.depth"]], 6)
            shrinkage = ifelse(!is.null(GBM.params[["shrinkage"]]),
                               GBM.params[["shrinkage"]], 0.01)
            n.minobsinnode = ifelse(!is.null(GBM.params[["n.minobsinnode"]]),
                                    GBM.params[["n.minobsinnode"]], 30)
            bag.fraction = ifelse(!is.null(GBM.params[["bag.fraction"]]),
                                  GBM.params[["bag.fraction"]], 0.5)
            train.fraction = ifelse(!is.null(GBM.params[["train.fraction"]]),
                                    GBM.params[["train.fraction"]], 1)
            cv.folds = ifelse(!is.null(GBM.params[["cv.folds"]]),
                              GBM.params[["cv.folds"]], 5)


            char_x_list = get_names(dat = dat_train[, x_list],
                                    types = c('character', 'factor'),
                                    ex_cols = c(target, obs_id, occur_time, ex_cols),
                                    get_ex = FALSE)
            if (length(char_x_list) > 0) {
                nr = nrow(dat_train)
                var_list = unique(c(target, obs_id, occur_time, x_list))
                dat_ts = rbind(dat_train[, var_list], dat_test[, var_list])
                dat_ts = one_hot_encoding(dat = dat_ts, cat_vars = char_x_list, na_act = FALSE)
                dat_train = dat_ts[1:nr,]
                dat_test = dat_ts[-c(1:nr),]
                x_list = get_x_list(x_list = NULL,
                                    dat_train = dat_train, dat_test = dat_test,
                                    ex_cols = c(target, obs_id, occur_time, ex_cols))
            }

            Formula = as.formula(paste(target, paste(x_list, collapse = ' + '), sep = ' ~ '))
            if (!is.null(seed)) set.seed(seed) else set.seed(46)
            gbm_model_new = gbm(
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
            best.iter = gbm.perf(gbm_model_new, method = "cv", plot.it = TRUE, oobag.curve = FALSE)
            imp_gbm = as.data.frame(summary(gbm_model_new, best.iter))
            dt_imp_GBM = data.frame(Feature = imp_gbm[, "var"],
                             Importance = round(imp_gbm[, 'rel.inf'], 5),
                            stringsAsFactors = FALSE)
            imp_vars_gbm = subset(dt_imp_GBM, dt_imp_GBM$Importance > 0)[, "Feature"]
            save_dt(dt_imp_GBM, file_name = paste0(model_name, "Featrue_Imp_GBM"),
                    dir_path = GBM_var_dir_path, note = FALSE)

            train_pred = dat_train[ c(obs_id, occur_time, target)]
            test_pred = dat_test[c(obs_id, occur_time, target)]
            train_pred$pred_GBM = round(predict(gbm_model_new, dat_train[, x_list],
                                                 best.iter, type = "response"), 5)
            test_pred$pred_GBM = round(predict(gbm_model_new, dat_test[, x_list],
                                                best.iter, type = "response"), 5)
            save_dt(train_pred, file_name = "train_pred_GBM",
                    dir_path = GBM_pred_dir_path, note = FALSE)
            save_dt(test_pred, file_name = "test_pred_GBM",
                    dir_path = GBM_pred_dir_path, note = FALSE)

            if (isTRUE(all.equal(names(train_pred), names(test_pred)))) {
                ks_psi_plot(train_pred = train_pred, test_pred = test_pred,
                            score = "pred_GBM", target = target,
                            plot_show = plot_show, g_width = 10,
                           save_data = TRUE,
                            dir_path = GBM_perf_dir_path, gtitle = paste0(model_name, ".GBM"))
                tb_pred = ks_table_plot(train_pred = train_pred, test_pred = test_pred,
                                        score = "pred_GBM", target = target,
                                        g = 20, g_width = 13, plot_show = FALSE, save_data = TRUE,
                                        dir_path = GBM_perf_dir_path,
                                        gtitle = paste0(model_name, ".GBM"))

                key_index = model_key_index(tb_pred)
                params_key_index = data.frame(c(GBM.params), key_index)
                save_dt(params_key_index, file_name = "params_GBM",
                        dir_path = GBM_perf_dir_path, append = TRUE, note = FALSE)
            }
            if (save_pmml) {
                cat(paste("[NOTE] convert GBM model to pmml.\n"))

                gbm_model_pmml = pmml(gbm_model_new,
                model.name = "GBM_Model",
                app.name = "R-PMML",
                description = "Gradient Boosting Decision Tree Model",
                copyright = NULL,
                transforms = NULL,
                unknownValue = NULL,
                parentInvalidValueTreatment = "returnInvalid",
                childInvalidValueTreatment = "asIs"
                )
                saveXML(gbm_model_pmml, file = paste0(GBM_model_dir_path, "/gbm_model.pmml"))
            }
            #save model
            save(gbm_model_new, file = paste0(GBM_model_dir_path, "/gbm_model.RData"))
            model_new = append(model_new, list(gbm_model = gbm_model_new), 1)
            rm(gbm_model_new)
        }

        if (any(algorithm == "RF")) {
            cat("[NOTE] RandomForest model is building.\n")

            RF_model_dir_path = paste(model_dir_path, "RF", sep = "/")
            RF_var_dir_path = paste(var_dir_path, "RF", sep = "/")
            RF_perf_dir_path = paste(perf_dir_path, "RF", sep = "/")
            RF_pred_dir_path = paste(pred_dir_path, "RF", sep = "/")
            RF_data_dir_path = paste(data_dir_path, "RF", sep = "/")
            if (!dir.exists(RF_model_dir_path)) dir.create(RF_model_dir_path)
            if (!dir.exists(RF_var_dir_path)) dir.create(RF_var_dir_path)
            if (!dir.exists(RF_perf_dir_path)) dir.create(RF_perf_dir_path)
            if (!dir.exists(RF_pred_dir_path)) dir.create(RF_pred_dir_path)
            if (!dir.exists(RF_data_dir_path)) dir.create(RF_data_dir_path)
            `%DO%` = if (parallel) `%dopar%` else `%do%`
            if (is.null(RF.params)) {
               RF.params = rf_params()
            }
            ntree = ifelse(!is.null(RF.params[["ntree"]]),
                           RF.params[["ntree"]], 100)
            nodesize = ifelse(!is.null(RF.params[["nodesize"]]),
                              RF.params[["nodesize"]], 30)
            samp_rate = ifelse(!is.null(RF.params[["samp_rate"]]),
                               RF.params[["samp_rate"]], 0.1)
            tune_rf = ifelse(!is.null(RF.params[["tune_rf"]]),
                             RF.params[["tune_rf"]], TRUE)

            char_x_list = get_names(dat = dat_train[, x_list],
                                    types = c('character', 'factor'),
                                    ex_cols = c(target, obs_id, occur_time, ex_cols),
                                    get_ex = FALSE)
            if (length(char_x_list) > 0) {
                nr = nrow(dat_train)
                var_list = unique(c(target, obs_id, occur_time, x_list))
                dat_ts = rbind(dat_train[, var_list], dat_test[, var_list])
                dat_ts = one_hot_encoding(dat = dat_ts, cat_vars = char_x_list, na_act = FALSE)
                dat_train = dat_ts[1:nr,]
                dat_test = dat_ts[-c(1:nr),]
                x_list = get_names(dat = dat_train,
                                   types = c('numeric', 'integer', 'double'),
                                   ex_cols = c(target, obs_id, occur_time, ex_cols),
                                   get_ex = FALSE)
            }

            dat_train[, target] = as.factor(as.character(dat_train[, target]))
            if (!is.null(seed)) set.seed(seed) else set.seed(46)
            if (tune_rf) {
                #Tune Random Forest Model
                tRF = foreach(n_tree = rep(round(10 / cores_num), cores_num),
                .combine = randomForest::combine,
                .packages = "randomForest") %DO% {
                    tuneRF(x = as.matrix(dat_train[, x_list]),
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
                           proximity = FALSE)
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
                                    .combine = randomForest::combine,
                                    .packages = "randomForest") %DO% {
                                      randomForest(Formula,
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

            save_dt(dt_imp_RF, file_name = paste0(model_name, "Feature_Imp_RF"),
                    dir_path = RF_var_dir_path, note = FALSE)

            train_pred = dat_train[c(obs_id, occur_time, target)]
            test_pred = dat_test[ c(obs_id, occur_time, target)]
            train_pred$pred_RF= round(predict(rf_model_new,
                                              dat_train[, vars_rf],
                                              type = c("prob"))[, 2], 5)
            test_pred$pred_RF = round(predict(rf_model_new,
                                              dat_test[, vars_rf],
                                              type = c("prob"))[, 2], 5)
            save_dt(train_pred, file_name = "train_pred_RF",
                    dir_path = RF_pred_dir_path, note = FALSE)
            save_dt(test_pred, file_name = "test_pred_RF",
                    dir_path = RF_pred_dir_path, note = FALSE)

            if (isTRUE(all.equal(names(train_pred), names(test_pred)))) {
                ks_psi_plot(train_pred = train_pred, test_pred = test_pred,
                            gtitle = paste0(model_name, ".RF"), target = target,
                            score = "pred_RF", plot_show = plot_show,
                            dir_path = RF_perf_dir_path, save_data = TRUE )
                tb_pred = ks_table_plot(train_pred = train_pred, test_pred = test_pred,
                                        target = target,score = "pred_RF",
                                        gtitle = paste0(model_name, ".RF"), g = 20, plot_show = FALSE,
                                        dir_path = RF_perf_dir_path, save_data = TRUE)
                key_index = model_key_index(tb_pred)
                params_key_index = data.frame(c(RF.params), key_index)
                save_dt(params_key_index, file_name = "params_RF",
                        dir_path = RF_perf_dir_path,
                        append = TRUE, note = FALSE)
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
#' @param best_lambda  Metheds of best lanmbda stardards using to filter variables by LASSO. There are 3 methods: ("lambda.auc", "lambda.ks", "lambda.sim_sign") . Default is  "lambda.auc".
#' @param obsweight An optional vector of  'prior weights' to be used in the fitting process. Should be NULL or a numeric vector. If you oversample or cluster diffrent datasets to training the LR model, you need to set this parameter to ensure that the probability of logistic regression output is the same as that before oversampling or segmentation. e.g.:There are 10,000 0 obs and 500 1 obs before oversampling or under-sampling, 5,000 0 obs and 3,000 1 obs after oversampling. Then this parameter should be set to c(10000/5000, 500/3000). Default is NULL..
#' @param forced_in Names of forced input variables. Default is NULL.
#' @param sp_values  Vaules will be in separate bins.e.g. list(-1, "Missing")  means that -1 & Missing as special values.Default is NULL.
#' @param lasso  Logical, if TRUE, variables filtering by LASSO. Default is TRUE.
#' @param vars_plot  Logical, if TRUE, plot distribution and correlation of input variables . Default is TRUE.
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
                     bins_control = list(bins_num = 10, bins_pct = 0.05, b_chi = 0.02,
                                         b_odds = 0.1, b_psi = 0.03, b_or = 0.15,
                                         mono = 0.2, odds_psi = 0.15, kc = 1),
                     best_lambda = "lambda.sim_sign",
                     sp_values = NULL,
                     forced_in = NULL, obsweight = c(1, 1), lasso = TRUE, vars_plot = TRUE,
                     step_wise = TRUE, score_card = TRUE,
                     cor_p = 0.8, iv_i = 0.02, psi_i = 0.1, ...) {
    structure(list(tree_control = tree_control, bins_control = bins_control,
              best_lambda = best_lambda,
              sp_values = sp_values,
                   forced_in = forced_in, obsweight = obsweight,lasso = lasso,
                   step_wise = step_wise, score_card = score_card, vars_plot = vars_plot,
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

xgb_params = function(nrounds = 1000,
                      params = list(max.depth = 6, eta = 0.1, min_child_weight = 1,
                                    subsample = 1,colsample_bytree = 1, gamma = 0,
                                    max_delta_step = 0, eval_metric = "auc",
                                    objective = "binary:logistic"),
                      early_stopping_rounds = 100, ...) {
    structure(list(nrounds = nrounds, params = params,
                   early_stopping_rounds = early_stopping_rounds))
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
                      bag.fraction = 0.5, train.fraction = 0.7, n.minobsinnode = 30,
                      cv.folds = 5, ...) {
    structure(list(n.trees = n.trees, interaction.depth = interaction.depth,
                   shrinkage = shrinkage, bag.fraction = bag.fraction,
                   train.fraction = train.fraction, n.minobsinnode = n.minobsinnode,
                   cv.folds = cv.folds))
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
rf_params = function(ntree = 100, nodesize = 30, samp_rate = 0.5,tune_rf = FALSE, ...) {
    structure(list(ntree = ntree, nodesize = nodesize,
                   samp_rate = samp_rate, tune_rf = tune_rf))
}


#' Score Transformation
#'
#' \code{score_transfer} is  for transfer woe to score.
#' @param model A data frame with x and target.
#' @param tbl_woe a data.frame with woe variables.
#' @param a  Base line of score.
#' @param b  Numeric.Increased scores from doubling Odds.
#' @param save_data Logical, save results in locally specified folder. Default is TRUE
#' @param file_name The name for periodically saved score file. Default is "dat_score".
#' @param dir_path  The path for periodically saved score file.  Default is "./data"
#' @return  A data.frame with variables which values transfered to score.
#' @examples
#' # dataset spliting
#' sub = cv_split(UCICreditCard, k = 30)[[1]]
#' dat = UCICreditCard[sub,]
#' #rename the target variable
#' dat = re_name(dat, "default.payment.next.month", "target")
#' dat = data_cleansing(dat, target = "target", obs_id = "ID", 
#' occur_time = "apply_date", miss_values =  list("", -1))
#' #train_ test pliting
#' train_test <- train_test_split(dat, split_type = "OOT", prop = 0.7,
#'                                 occur_time = "apply_date")
#' dat_train = train_test$train
#' dat_test = train_test$test
#' #get breaks of all predictive variables
#' x_list = c("PAY_0", "LIMIT_BAL", "PAY_AMT5", "EDUCATION", "PAY_3", "PAY_2")
#' breaks_list <- get_breaks_all(dat = dat_train, target = "target",
#'                               x_list = x_list, occur_time = "apply_date", ex_cols = "ID", 
#' save_data = TRUE, note = FALSE)
#' #woe transforming
#' train_woe = woe_trans_all(dat = dat_train,
#'                           target = "target",
#'                           breaks_list = breaks_list,
#'                           woe_name = FALSE)
#' test_woe = woe_trans_all(dat = dat_test,
#'                        target = "target",
#'                          breaks_list = breaks_list,
#'                          note = FALSE)
#' Formula = as.formula(paste("target", paste(x_list, collapse = ' + '), sep = ' ~ '))
#' set.seed(46)
#' lr_model = glm(Formula, data = train_woe[, c("target", x_list)], family = binomial(logit))
#' #get LR coefficient
#' dt_imp_LR = get_logistic_coef(lg_model = lr_model, save_data = TRUE)
#' bins_table = get_bins_table_all(dat = dat_train, target = "target",
#'                                occur_time = "apply_date",
#'                                 x_list = x_list,
#'                                breaks_list = breaks_list, note = FALSE)
#' #score card
#' LR_score_card <- get_score_card(lg_model = lr_model, bins_table, target = "target")
#' #scoring
#' train_pred = dat_train[, c("ID", "apply_date", "target")]
#' test_pred = dat_test[, c("ID", "apply_date", "target")]
#' train_pred$pred_LR = score_transfer(model = lr_model,
#'                                                     tbl_woe = train_woe,
#'                                                     save_data = TRUE)[, "score"]
#' 
#' test_pred$pred_LR = score_transfer(model = lr_model,
#' tbl_woe = test_woe, save_data = TRUE)[, "score"]
#' @export

score_transfer <- function(model, tbl_woe, a = 600, b = 50,
                           file_name = NULL, dir_path = tempdir(),
                           save_data = TRUE) {
    coef = model$coefficients
    glm_vars <- names(coef)[-1]
    A = a
    B = b / log(2)
    base_score = A - B * coef[1]

    tbl_woe = tbl_woe[c(glm_vars)]
    score_name = c()
    for (i in glm_vars) {
        tbl_woe[i] = (-1) * B * tbl_woe[i] * coef[i]
        score_name[i] = gsub("_woe", "_score", i)
    }
    names(tbl_woe) = score_name
    tbl_woe$score <- apply(tbl_woe[1:length(tbl_woe)], MARGIN = 1, function(x) sum(x))
    tbl_woe$score <- round(tbl_woe$score + base_score, 2)
    if (save_data) {
        dir_path = ifelse(!is.character(dir_path),
                      tempdir(), dir_path)
        if (!dir.exists(dir_path)) dir.create(dir_path)
        if (!is.character(file_name)) file_name = NULL

        save_dt(tbl_woe, file_name = ifelse(is.null(file_name), "dat.score", paste(file_name, "dat.score", sep = ".")), dir_path = dir_path, note = FALSE)
    }
    tbl_woe
}

#' Score Card
#'
#' \code{get_score_card} is  for generating a stardard scorecard
#' @param lg_model An object of glm model.
#' @param target The name of target variable.
#' @param bins_table a data.frame generated by \code{\link{get_bins_table}}
#' @param a  Base line of score.
#' @param b  Numeric.Increased scores from doubling Odds.
#' @param save_data Logical, save results in locally specified folder. Default is TRUE
#' @param file_name  The name for periodically saved scorecard file. Default is "LR_Score_Card".
#' @param dir_path  The path for periodically saved scorecard file. Default is "./model"
#' @return  scorecard
#' @export
#' @examples
#' # dataset spliting
#' sub = cv_split(UCICreditCard, k = 30)[[1]]
#' dat = UCICreditCard[sub,]
#' #rename the target variable
#' dat = re_name(dat, "default.payment.next.month", "target")
#' dat = data_cleansing(dat, target = "target", obs_id = "ID", 
#' occur_time = "apply_date", miss_values =  list("", -1))
#' #train_ test pliting
#' train_test <- train_test_split(dat, split_type = "OOT", prop = 0.7,
#'                                 occur_time = "apply_date")
#' dat_train = train_test$train
#' dat_test = train_test$test
#' #get breaks of all predictive variables
#' x_list = c("PAY_0", "LIMIT_BAL", "PAY_AMT5", "EDUCATION", "PAY_3", "PAY_2")
#' breaks_list <- get_breaks_all(dat = dat_train, target = "target",
#'                               x_list = x_list, occur_time = "apply_date", ex_cols = "ID", 
#' save_data = TRUE, note = FALSE)
#' #woe transforming
#' train_woe = woe_trans_all(dat = dat_train,
#'                           target = "target",
#'                           breaks_list = breaks_list,
#'                           woe_name = FALSE)
#' test_woe = woe_trans_all(dat = dat_test,
#'                        target = "target",
#'                          breaks_list = breaks_list,
#'                          note = FALSE)
#' Formula = as.formula(paste("target", paste(x_list, collapse = ' + '), sep = ' ~ '))
#' set.seed(46)
#' lr_model = glm(Formula, data = train_woe[, c("target", x_list)], family = binomial(logit))
#' #get LR coefficient
#' dt_imp_LR = get_logistic_coef(lg_model = lr_model, save_data = TRUE)
#' bins_table = get_bins_table_all(dat = dat_train, target = "target",
#'                                occur_time = "apply_date",
#'                                 x_list = x_list,
#'                                breaks_list = breaks_list, note = FALSE)
#' #score card
#' LR_score_card <- get_score_card(lg_model = lr_model, bins_table, target = "target")
#' #scoring
#' train_pred = dat_train[, c("ID", "apply_date", "target")]
#' test_pred = dat_test[, c("ID", "apply_date", "target")]
#' train_pred$pred_LR = score_transfer(model = lr_model,
#'                                                     tbl_woe = train_woe,
#'                                                     save_data = TRUE)[, "score"]
#' 
#' test_pred$pred_LR = score_transfer(model = lr_model,
#' tbl_woe = test_woe, save_data = TRUE)[, "score"]

get_score_card <- function(lg_model, target, bins_table, a = 600, b = 50,
                           file_name = NULL, dir_path = tempdir(),
                           save_data = TRUE) {
    coef = lg_model$coefficients
    glm_vars <- gsub("_woe", "", names(coef))
    names(coef) <- glm_vars
    A = a
    B = b / log(2)
    base_score = A - B * coef[1]
    dt_score_card <- bins_table[which(as.character(bins_table[, "Feature"]) %in% glm_vars),
                                c("Feature", "cuts", "bins", "woe")]

    for (i in glm_vars) {
        dt_score_card[which(as.character(dt_score_card[, "Feature"]) == i), "coefficient"] = round(coef[i], 5)
    }

    for (i in glm_vars) {
        dt_score_card[which(as.character(dt_score_card[, "Feature"]) == i), "score"] =
        round((-1) * B * as.numeric(dt_score_card[which(as.character(dt_score_card[, "Feature"]) == i), "woe"]) * coef[i], 2)
    }
    Intercept = c("Intercept", "", "", "", round(coef[1], 5), paste("Base:", round(base_score)))
    dt_score_card = rbind(Intercept, dt_score_card)
    if (save_data) {
        dir_path = ifelse(!is.character(dir_path),
                      tempdir(), dir_path)
        if (!dir.exists(dir_path)) dir.create(dir_path)
        if (!is.character(file_name)) file_name = NULL
        save_dt(dt_score_card, file_name = ifelse(is.null(file_name), "scorecard", paste(file_name, "scorecard", sep = ".")), dir_path = dir_path, note = FALSE)
    }
    dt_score_card
}



#' get logistic coef
#'
#' \code{get_logistic_coef} is  for geting logistic coefficient.
#' @param lg_model  An object of logistic model.
#' @param save_data Logical, save the result or not. Default is TRUE.
#' @param file_name  The name for periodically saved coefficient file.  Default is "LR_coef".
#' @param dir_path  The Path for periodically saved coefficient file. Default is "./model".
#' @return  A data.frame with logistic coefficients.
#' @examples
#' # dataset spliting
#' sub = cv_split(UCICreditCard, k = 30)[[1]]
#' dat = UCICreditCard[sub,]
#' #rename the target variable
#' dat = re_name(dat, "default.payment.next.month", "target")
#' dat = data_cleansing(dat, target = "target", obs_id = "ID", 
#' occur_time = "apply_date", miss_values =  list("", -1))
#' #train_ test pliting
#' train_test <- train_test_split(dat, split_type = "OOT", prop = 0.7,
#'                                 occur_time = "apply_date")
#' dat_train = train_test$train
#' dat_test = train_test$test
#' #get breaks of all predictive variables
#' x_list = c("PAY_0", "LIMIT_BAL", "PAY_AMT5", "EDUCATION", "PAY_3", "PAY_2")
#' breaks_list <- get_breaks_all(dat = dat_train, target = "target",
#'                               x_list = x_list, occur_time = "apply_date", ex_cols = "ID", 
#' save_data = TRUE, note = FALSE)
#' #woe transforming
#' train_woe = woe_trans_all(dat = dat_train,
#'                           target = "target",
#'                           breaks_list = breaks_list,
#'                           woe_name = FALSE)
#' test_woe = woe_trans_all(dat = dat_test,
#'                        target = "target",
#'                          breaks_list = breaks_list,
#'                          note = FALSE)
#' Formula = as.formula(paste("target", paste(x_list, collapse = ' + '), sep = ' ~ '))
#' set.seed(46)
#' lr_model = glm(Formula, data = train_woe[, c("target", x_list)], family = binomial(logit))
#' #get LR coefficient
#' dt_imp_LR = get_logistic_coef(lg_model = lr_model, save_data = TRUE)
#' bins_table = get_bins_table_all(dat = dat_train, target = "target",
#'                                occur_time = "apply_date",
#'                                 x_list = x_list,
#'                                breaks_list = breaks_list, note = FALSE)
#' #score card
#' LR_score_card <- get_score_card(lg_model = lr_model, bins_table, target = "target")
#' #scoring
#' train_pred = dat_train[, c("ID", "apply_date", "target")]
#' test_pred = dat_test[, c("ID", "apply_date", "target")]
#' train_pred$pred_LR = score_transfer(model = lr_model,
#'                                                     tbl_woe = train_woe,
#'                                                     save_data = TRUE)[, "score"]
#' 
#' test_pred$pred_LR = score_transfer(model = lr_model,
#' tbl_woe = test_woe, save_data = TRUE)[, "score"]
#' @importFrom car vif
#' @export
get_logistic_coef = function(lg_model, file_name = NULL,
                             dir_path = tempdir(), save_data = TRUE) {
    lg_coef = data.frame(summary(lg_model)$coefficients)
    lg_coef[4] = round(lg_coef[4], 5)
    lg_coef[, "Feature"] = row.names(lg_coef)
    if (length(row.names(lg_coef)) > 2) {
        lg_coef[-1, "vif"] = car::vif(lg_model)
    } else {
        lg_coef[-1, "vif"] = 0
    }
    names(lg_coef) <- c("estimate", "std.error", "Z_value", "P_value", "Feature", "vif")
    lg_coef = lg_coef[c("Feature", "estimate", "std.error", "Z_value", "P_value", "vif")]
    if (save_data) {
        dir_path = ifelse(!is.character(dir_path), tempdir(), dir_path)
        if (!dir.exists(dir_path)) dir.create(dir_path)
        if (!is.character(file_name)) file_name = NULL
        save_dt(lg_coef, file_name = ifelse(is.null(file_name), "logistic.coef", paste(file_name, "logistic.coef", sep = ".")), dir_path = dir_path, note = FALSE)
    }
    return(lg_coef)
}
