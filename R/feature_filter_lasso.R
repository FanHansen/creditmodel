#' Selected Variables by LASSO
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
#' @param best_lambda  Best lanmbda stardards. one of ("lambda.min", "lambda.1se", "lambda.05se" , "lambda.sim_sign"). Default is  "lambda.min". 
#' @param plot.it Logical, shrinkage plot. Default is TRUE.
#' @param parallel Logical, parallel computing. Default is FALSE.
#' @param seed  Random number seed. Default is 46.
#' @param save_data Logical, save results in locally specified folder. Default is TRUE
#' @param file_name  The name for periodically saved results files. Default is "Feature_selected_LASSO".
#' @param dir_path The path for periodically saved results files. Default is "./variable".
#' @return A list of filtered x variables by lasso.
#' @examples
#' \dontrun{
#' lasso_filter(dat_train = UCICreditCard, 
#' target = "default.payment.next.month", 
#' best_lambda = "lambda.min")
#' }
#'
#' @importFrom glmnet cv.glmnet glmnet
#' @import ggplot2
#' @importFrom gridExtra arrangeGrob
#' @export

lasso_filter <- function(dat_train, dat_test = NULL, target = NULL, x_list = NULL, pos_flag = list(1, "1", "bad"),
ex_cols = NULL, best_lambda = "lambda.min", sim_sign = "negtive", save_data = TRUE, parallel = FALSE,
plot.it = TRUE, seed = 46, file_name = "Feature_selected_LASSO", dir_path = "./variable") {
    cat("[NOTE] Dimension reduction of variables with Lasso. \n")
    opt = options("warn" = -1, scipen = 100, stringsAsFactors = FALSE, digits = 10) # suppress warnings
    if (parallel) {
        parallel <- start_parallel_computing(parallel)
        stopCluster <- TRUE
    } else {
        parallel <- stopCluster <- FALSE
    }
    on.exit(if (parallel & stopCluster) stop_parallel_computing(attr(parallel, "cluster")))

    dat_train = checking_data(dat = dat_train, target = target)
    x_list = get_x_list(x_list = x_list, dat_train = dat_train, dat_test = dat_test, ex_cols = c(target,  ex_cols))
    var_list = unique(c(target, x_list))
    dat_ts = rbind(dat_train[, var_list], dat_test[, var_list])
    dat_ts = low_variance_filter(dat = dat_ts, lvp = 1, note = FALSE)
    if (any(is.na(dat_ts[x_list]))) {
        dat_ts = process_nas(dat = dat_ts, x_list = var_list, ex_cols = c(target, ex_cols), parallel = parallel, method = "median")
    }

    char_x_list = get_names(dat = dat_ts, types = c('character', 'factor'), ex_cols = c(target,ex_cols), get_ex = FALSE)
    if (length(char_x_list) > 0) {
        dat_ts = one_hot_encoding(dat = dat_ts, cat_vars = char_x_list, na_act = FALSE, note = FALSE)
    }
    num_x_list = get_names(dat = dat_ts, types = c('numeric', 'integer', 'double'), ex_cols = c(target, ex_cols), get_ex = FALSE)

    if (is.null(dat_test)) {
        train_test = train_test_split(dat_ts, split_type = "Random", prop = 0.3, seed = 46, note = FALSE, save_data = FALSE)
    } else {
        nr = nrow(dat_train)
        train_test = train_test_split(dat_ts, split_type = "byRow", prop = nr / nrow(dat_ts), seed = 46, note = FALSE, save_data = FALSE)
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

    cat(paste("[NOTE] select best lanmbda by max K-S. \n"))
    bst_lanmbda = bst_lanmbda_id = best_ks = best_score = bst_coefficients = bst_vars = shrink_vars = c()
    test_pre = predict(lasso_model, newx = x_test, s = lasso_model$lambda, type = "response") # make predictions
    KS = apply(test_pre, 2, function(x) round(ks_value(score = x, target = y_test, g = 20), 2))
    best_ks = max(KS, na.rm = TRUE)[1]
    bst_ks_lanmbda = lasso_model$lambda[which(KS == best_ks)[1]]
    #cross validation can also be used to select lambda.
    cat(paste("[NOTE] cross validation to maximize AUC.\n"))
    if (!is.null(seed)) set.seed(seed) else set.seed(46)
    lasso_model_cv = cv.glmnet(y = y_train,
                       x = x_train,
                       family = "binomial",
                       nfolds = 5,
                       type.measure = "auc",
                       parallel = parallel)

    lasso_lambda_cv = get_best_lambda(lasso_model = lasso_model,
                                      lasso_model_cv = lasso_model_cv, type.measure = "auc", sim_sign = sim_sign)
    b_lambda = max(lasso_lambda_cv[[best_lambda]], lasso_lambda_cv[["lambda.min"]], bst_ks_lanmbda)
    coefficients = coef(lasso_model, s = b_lambda)
    #variable coefficients  which coefficient are not zero  
    bst_vars = which(coefficients[-1] != 0)
    shrinkage_vars = num_x_list[bst_vars]
    dat_lasso = de_one_hot_encoding(dat_one_hot = dat_train[shrinkage_vars], cat_vars = char_x_list, na_act = TRUE, note = FALSE)
    lasso_vars = get_x_list(x_list = x_list, dat_train = dat_lasso, ex_cols = c(target, ex_cols))
    KS_lanmbda = data.frame(vars_num = lasso_model$df, KS = KS, lanmbda = lasso_model$lambda)
    save_dt(KS_lanmbda, file_name = "KS_lanmbda", dir_path = dir_path)
    #dev.print(png, file = paste0(dir_path, "/", file_name, "_coef&auc.png"), width = 800, height = 500)
    ks_lanm = KS_lanmbda[-1,]
    plot_ks = ggplot(data = ks_lanm, aes(x = -log(ks_lanm$lanmbda), y = ks_lanm$KS)) +
            geom_line(aes(color = "KS"),
              position = position_dodge(width = 0.5), #color = love_color("sky_blue"),
              size = 1) + geom_point(
               position = position_dodge(width = 0.5),
               fill = 'white', color = love_color("deep_red"),
              size = 1, shape = 21) +
            geom_line(aes(y = -ks_lanm$lanmbda[seq(1, by = ceiling(length(ks_lanm$vars_num) / 20),
            length.out = length(ks_lanm$vars_num))] * max(ks_lanm$KS) * 10 + max(ks_lanm$KS),
        x = -log(ks_lanm$lanmbda)[seq(1, by = ceiling(length(ks_lanm$vars_num) / 20),
        length.out = length(ks_lanm$vars_num))], color = "Number of Variables"),
               position = position_dodge(width = 0.5), # color = love_color("light_yellow"),
               size = 1) + geom_point(aes(y = -ks_lanm$lanmbda[seq(1, by = ceiling(length(ks_lanm$vars_num) / 20),
               length.out = length(ks_lanm$vars_num))] * max(KS) * 10 + max(KS),
        x = -log(ks_lanm$lanmbda)[seq(1, by = ceiling(length(ks_lanm$vars_num) / 20),
        length.out = length(ks_lanm$vars_num))]),
               position = position_dodge(width = 0.5),
               fill = 'white', color = love_color("light_grey"), size = 1, shape = 21) +
            geom_text(aes(y = -ks_lanm$lanmbda[seq(1, by = ceiling(length(ks_lanm$vars_num) / 20),
            length.out = length(ks_lanm$vars_num))] * max(KS) * 10 + max(ks_lanm$KS),
        x = -log(ks_lanm$lanmbda)[seq(1, by = ceiling(length(ks_lanm$vars_num) / 20),
        length.out = length(ks_lanm$vars_num))],
        label = paste(ks_lanm$vars_num[seq(1, by = ceiling(length(ks_lanm$vars_num) / 20),
        length.out = length(ks_lanm$vars_num))])),
              position = position_dodge(width = 0.5),
              colour = 'black', size = 2.5, vjust = -1) +
            annotate(geom = 'text', x = -log(bst_ks_lanmbda),
          y = best_ks + 0.01,
              label = paste("max K-S:", best_ks, sep = " ")) +
            scale_color_manual(values = c('KS' = love_color("sky_blue"), 'Number of Variables' = love_color("light_yellow"))) +
            ylab("K-S") +
            ggtitle("LASSO: K-S with different lanmbda") +
            plot_theme(legend.position = "top", angle = 0)
    ggsave(device = "png", filename =  paste0(dir_path, "/", "lasso_max_ks", ".png"), plot = plot_ks, dpi = "retina", width = 8)
    if (plot.it) {
        par(mar = c(5, 3, 3, 1) + 0.1, mfrow = c(1, 2))
        plot(lasso_model, xvar = "lambda")
        plot(lasso_model_cv)  
        plot(plot_ks)
    }
    if (save_data) {
        save_dt(lasso_vars, file_name = file_name, dir_path = dir_path, as_list = TRUE)
    }
    options(opt) # reset warnings
    return(lasso_vars)
}



#' get_best_lambda
#' \code{plot_theme} is for get Best lambda required in lasso_filter. This function required in \code{lasso_filter}
#' @param lasso_model A lasso model genereted by glmnet.
#' @param lasso_model_cv A cv lasso model genereted by cv.glmnet
#' @param type.measure  Default is  "auc".
#' @param sim_sign  Default is "negtive". This is related to pos_plag. If pos_flag equals 1 or bad, the value must be set to negetive. If pos_flag equals 0 or good, the value must be set to positive.
#' @return Four lanmbda values with different thresholds.
#' @details
#' lambda.min give a model with the best performance but the least number of independent variable.
#' lambda.lse give the simplest model in the range of a standard deviation
#' lambda.05se give the simplest model in the range of a half standard deviation
#' lambda.sim_sign give the model with the same positive or negetive coefficients of all variables.
#' lambda.sim_sign give the model with the same positive or negetive coefficients of all variables.
#' 
#' @export
get_best_lambda <- function(lasso_model, lasso_model_cv, type.measure = "auc", sim_sign = "negtive") {
    lambda = lasso_model_cv$glmnet.fit$lambda
    if (type.measure == "auc") {
        cvm = -lasso_model_cv$cvm
    } else {
        cvm = lasso_model_cv$cvm
    }
    cvsd = lasso_model_cv$cvsd
    # min lambda
    cvmin = min(cvm, na.rm = TRUE)
    idmin = cvm <= cvmin
    lambda.min = max(lambda[idmin], na.rm = TRUE)
    #1 standard cv_folds lambda
    id_min = match(lambda.min, lambda)
    semin = (cvm + cvsd)[id_min]
    idmin = cvm <= semin
    lambda.1se = max(lambda[idmin], na.rm = TRUE)
    #0.5 standard cv_folds lambda
    se05min = (cvm + 0.5 * cvsd)[id_min]
    idmin = cvm <= se05min
    lambda.05se = max(lambda[idmin], na.rm = TRUE)
    #same sign  lambda
    beta_x = lasso_model$beta@x
    beta_i = lasso_model$beta@i
    coefs <- c()
    if (sim_sign == "negtive") {
        for (i in 1:length(lambda)) {
            coefs[i] = any(coef(lasso_model, s = lambda[i])[-1] > 0)
        }
    } else {
        for (i in 1:length(lambda)) {
            coefs[i] = any(coef(lasso_model, s = lambda[i])[-1] < 0)
        }
    }
    if (length(coefs) > 0) {
        lambda.sim_sign = lambda[sum(!coefs)]
    } else {
        stop(paste("no variables coeficeient is positive or negetive"))
    }
    list(lambda.min = lambda.min, lambda.1se = lambda.1se, lambda.05se = lambda.05se, lambda.sim_sign = lambda.sim_sign)
}