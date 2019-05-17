#' Plot Independent Variables
#'
#' \code{plot_vars} is used for independent variables with target variable ploting.
#' \code{get_plots} can loop through plots for all specified independent variables.
#' @param dat_train A data.frame with independent variables and target variable.
#' @param dat_test  A data.frame of test data. Default is NULL.
#' @param target The name of target variable.
#' @param x_list Names of independent variables.
#' @param x  The name of an independent variable.
#' @param ex_cols A list of excluded variables. Regular expressions can also be used to match variable names. Default is NULL.
#' @param pos_flag Value of positive class, Default is "1".
#' @param breaks_list A table containing a list of splitting points for each independent variable. Default is NULL.
#' @param occur_time The name of the variable that represents the time at which each observation takes place.
#' @param g_width  The width of graphs.
#' @param best  Logical, merge initial breaks to get optimal breaks for binning.
#' @param equal_bins  Logical, generates initial breaks for equal frequency binning.
#' @param g  Number of initial breakpoints for equal frequency binning.
#' @param tree_control  Parameters of using Decision Tree to segment initial breaks. See detials: \code{\link{get_tree_breaks}}
#' @param bins_control  Parameters  used to control binning.  See detials: \code{\link{select_best_class}}, \code{\link{select_best_breaks}}
#' @param parallel Logical, parallel computing. Default is FALSE.
#' @param plot_show Logical, show model performance in current graphic device. Default is FALSE.
#' @param save_data Logical, save results in locally specified folder. Default is TRUE
#' @param dir_path The path for periodically saved graphic files.
#' @param file_name  The name for periodically saved data file. Default is NULL.
#' @examples
#' train_test <- train_test_split(UCICreditCard[1:1000,], split_type = "Random",
#'  prop = 0.8, save_data = TRUE)
#' dat_train = train_test$train
#' dat_test = train_test$test
#' get_plots(dat_train[, c(8, 26)], dat_test = dat_test[, c(8, 26)],
#' target = "default.payment.next.month")
#' @import ggplot2
#' @importFrom gridExtra arrangeGrob tableGrob
#' @importFrom dplyr group_by mutate summarize  summarise n  count %>% filter mutate_if distinct ungroup
#' @importFrom data.table melt
#' @export


get_plots <- function(dat_train, dat_test = NULL, x_list = NULL,
                      target = NULL, ex_cols = NULL, breaks_list = NULL,
                      pos_flag = NULL, occur_time = NULL, equal_bins = FALSE, best = TRUE,
                      g = 20, tree_control = NULL, bins_control = NULL, plot_show = TRUE,
                      save_data = TRUE, file_name = NULL,
                      parallel = FALSE, g_width = 8, dir_path = tempdir()) {

   cat(paste("[NOTE] Plot input variables ...."), "\n")
    opt = options(stringsAsFactors = FALSE) #
    if (save_data) {
        dir_path = paste0(dir_path, "/variable_plot/")
        if (!dir.exists(dir_path)) dir.create(dir_path)
        if (dir.exists(dir_path)) { file.remove(list.files(dir_path, recursive = TRUE, full.names = TRUE)) }
        }
    dat_train = checking_data(dat = dat_train, target = target, pos_flag = pos_flag)
    if (is.null(x_list)) {
        if (is.null(x_list)) {
            if (!is.null(breaks_list)) {
                x_list = unique(as.character(breaks_list[, "Feature"]))
            } else {
                x_list = get_x_list(x_list = x_list,
                                    dat_train = dat_train,
                                    dat_test = dat_test,
                                    ex_cols = c(target, occur_time, ex_cols))
            }
        }
    } else {
       x_list <- gsub("_woe$|_pred$","", x_list)
    }
    if (!is.null(dat_test)) {
        dat_test = checking_data(dat = dat_test, target = target, pos_flag = pos_flag)
        x_list = get_x_list(x_list = x_list, dat_train = dat_train, dat_test = dat_test,
                            ex_cols = c(target, occur_time, ex_cols))
        com_list = unique(c(target, occur_time, x_list))
        dat_train = dat_train[, com_list]
        dat_test = dat_test[, com_list]
        dat_ts = rbind(dat_train, dat_test)
        if (all(unique(dat_ts[, target]) != c("0", "1"))) {
            if (!is.null(pos_flag)) {
                dat_ts$target = ifelse(dat_ts[, target] %in% pos_flag, "1", "0")
            } else {
                pos_flag = list("1", 1, "bad", "positive")
                dat_ts$target = ifelse(dat_ts[, target] %in% pos_flag, "1", "0")
            }
            if (length(unique(dat_ts$target)) == 1) {
                stop("pos_flag is missing.\n")
            }
        } else {
            dat_ts$target = as.character(dat_ts[, target])
        }
        nr = nrow(dat_train)
        train_test = train_test_split(dat_ts, split_type = "byRow", prop = nr / nrow(dat_ts),
          occur_time = occur_time, seed = 46, save_data = TRUE, note = FALSE)
        dat_train = train_test$train
        dat_test = train_test$test
    } else {
        if (all(unique(dat_train[, target]) != c("0", "1"))) {
            if (!is.null(pos_flag)) {
                dat_train$target = ifelse(dat_train[, target] %in% pos_flag, "1", "0")
            } else {
                pos_flag = list("1", 1, "bad", "positive")
                dat_train$target = ifelse(dat_train[, target] %in% pos_flag, "1", "0")
            }
            if (length(unique(dat_train$target)) == 1) {
                stop("pos_flag is missing.\n")
            }
        } else {
            dat_train$target = as.character(dat_train[, target])
        }
        train_test = train_test_split(dat = dat_train, split_type = "OOT", prop = 0.7,
                                      occur_time = occur_time, seed = 46, save_data = TRUE)
        dat_train = train_test$train
        dat_test = train_test$test
    }

    df_ae_list = loop_function(func = plot_vars, x_list = x_list,
                               args = list(dat_train = dat_train, dat_test = dat_test,
                                           target = target, breaks_list = breaks_list,
                                           pos_flag = pos_flag,
                                           occur_time = occur_time,
                                           equal_bins = equal_bins, best = best, g = g,
                                           tree_control = tree_control,
                                           bins_control = bins_control,
                                           plot_show = plot_show,
                                           g_width = g_width, dir_path = dir_path, save_data = save_data),
                               bind = "rbind", parallel = parallel, as_list = FALSE)
    options(opt) # reset
    if (save_data) {
        save_dt(df_ae_list, dir_path = dir_path, file_name = ifelse(is.null(file_name), "plot.table", paste(file_name, "plot.table", sep = ".")), append = FALSE, note = FALSE)
    }

    return(df_ae_list)
}



#' @rdname get_plots
#' @export

plot_vars <- function(dat_train, x, dat_test = NULL, target = "target",
                      g_width = 8,breaks_list = NULL,
                      pos_flag = list("1", 1, "bad", "positive"),
                      occur_time = NULL,
                      equal_bins = FALSE, best = TRUE,
                      g = 20, tree_control = NULL,
                      bins_control = NULL,
                      plot_show = TRUE,
                      save_data = TRUE,
                      dir_path = tempdir()) {

    digits_x = digits_num(dat_train[, x])
    opt = options(stringsAsFactors = FALSE, digits = digits_x + 1) #

    df_ae <- get_psi_iv(dat = dat_train, dat_test = dat_test,
                        x = x, target = target, pos_flag = pos_flag,
                        breaks_list = breaks_list,
                        occur_time = occur_time,
                        equal_bins = equal_bins,
                        tree_control = tree_control,
                        bins_control = bins_control,
                        bins_total = FALSE,
                        best = best, g = g, as_table = TRUE,
                        note = FALSE, bins_no = TRUE)
    ae_total <- data.table::melt(df_ae[c("bins", "%actual", "%expected")],
                                 id.vars = c("bins"),
                                 variable.name = "actual_expected",
                                 value.name = "value")
    ae_1 <- data.table::melt(df_ae[c("bins", "%actual_1", "%expected_1")],
                               id.vars = c("bins"),
                               variable.name = "actual_expected",
                               value.name = "value")
    xn = NULL
    if (class(dat_train[, x]) %in% c("numeric", "double","integer") && length(unique(dat_train[,x]))>10) {
        med <- dat_train %>%
          dplyr::mutate(xn = dat_train[, x]) %>%
          dplyr::group_by(target) %>%
          dplyr::summarise(grp.mean = quantile(xn, 0.5, na.rm = TRUE, type = 3))
        none_na_num <- sum(!is.na(dat_train[, x]))
        tbl_x <- table(dat_train[, x])
        x_unique_value <- as.double(names(tbl_x))
        cum_sum <- cumsum(tbl_x)
        cuts_sum <- approx(cum_sum, x_unique_value, xout = (1:100) * none_na_num / 100,
                           method = "constant", rule = 2, f = 1)$y
        if (length(unique(cuts_sum) )> 10) {
            x_cuts <- cuts_sum[length(cuts_sum) - 1]
            dat_train <- subset(dat_train, dat_train[, x] < x_cuts)
        }
        #desity plot
        plot_1 <- ggplot(dat_train, aes(x = dat_train[, x])) +
          geom_density(aes(fill = dat_train$target), alpha = 0.3) +
          stat_density(geom = "line", position = "identity", size = 0.6,
                       aes(color = dat_train$target)) +
          scale_fill_manual(values = c('1' = love_color("shallow_red"),
                                       '0' = love_color("sky_blue"))) +
          scale_color_manual(values = c('1' = love_color("shallow_red"),
                                        '0' = love_color("sky_blue"))) +
          geom_vline(data = med,linetype = "dashed", size = 0.7,
                     aes(xintercept = med$grp.mean, color = med$target)) +
          xlab(x) +
          ggtitle(paste("Density of", x)) +
          plot_theme(legend.position = c(.9, .9), title_size = 9,
                     axis_title_size = 8)

    } else {
        #relative frequency histogram
        data1 <- dat_train %>%
          dplyr::mutate(xn = dat_train[, x]) %>%
          dplyr::filter(target %in% c("0", "1")) %>%
          dplyr::group_by(xn) %>% dplyr::count(xn, target) %>%
          dplyr::mutate(percent = n / sum(n))

        plot_1 <- ggplot(data1,aes(x = data1$xn, y = data1$percent, fill = reorder(data1$target, n))) +
          geom_bar(stat = "identity", position = position_stack()) +
          geom_text(aes(label = paste(as_percent(data1$percent, digits = 3))),
                    size = ifelse(nrow(ae_total) > 10, 2.1,
                                ifelse(nrow(ae_total) > 5, 2.5,
                                       ifelse(nrow(ae_total) > 3, 3, 3.3))), vjust = 1, colour = 'white', position = position_stack()) +
          guides(fill = guide_legend(reverse = F)) +
          ggtitle(paste("Relative Frequency of", x)) +
          ylab("Percent") + xlab(x) +
          scale_fill_manual(values = c('0' = love_color("green_cyan"),
                                       '1' = love_color("deep_grey"))) +
          plot_theme(legend.position = "top", title_size = 9,
                                       axis_title_size = 8)
    }

    plot_2 <- ggplot(ae_total, aes(x = ae_total$bins,
                                   y = ae_total$value,
                                   fill = ae_total$actual_expected)) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.7)) +
      geom_text(aes(y = ae_total$value,
                    label = paste(as_percent(ae_total$value, digits = 3))),
                position = position_dodge(width = 0.7),
                size = ifelse(nrow(ae_total) > 10, 2.6,
                                ifelse(nrow(ae_total) > 5, 2.8,
                                       ifelse(nrow(ae_total) > 3, 3, 3.2))), vjust = 1, hjust = 0.3, colour = "white") +
      geom_line(aes(x = factor(ae_1[[1]]),
                    y = as.numeric(ae_1$value) * max(ae_total$value) * 4,
                    color = ae_1$actual_expected,
                    linetype = ae_1$actual_expected,
                    group = ae_1$actual_expected),
                position = position_dodge(width = 0.5),size = 1) +
      geom_point(aes(y = as.numeric(ae_1$value) * max(ae_total$value) * 4,
                     color = ae_1$actual_expected,
                     group = ae_1$actual_expected),
                 position = position_dodge(width = 0.5),
                 fill = 'white', color = love_color("deep_red"), size = 2, shape = 21) +
      geom_text(aes( y = as.numeric(ae_1$value) * max(ae_total$value) * 4,
                     label = paste(as_percent(ae_1$value, digits = 3))),
                position = position_dodge(width = 0.5),
                colour = 'black',
                size = ifelse(nrow(ae_total) > 10, 2.6,
                                ifelse(nrow(ae_total) > 5, 2.8,
                                       ifelse(nrow(ae_total) > 3, 3, 3.2))), vjust = -0.1) +
      annotate(geom = 'text',
               x = dim(ae_total)[1] / 3,
               y = max(c(ae_total$value, as.numeric(ae_1$value) * max(ae_total$value) * 4)) + 0.09,
               label = paste(paste("IV:", sum(df_ae$IVi)), paste('PSI:',sum(df_ae$PSIi)), sep = "   ")) +
      scale_fill_manual(values = c('%actual' = love_color("deep_grey"),
                                   '%expected' = love_color("light_yellow"))) +
      scale_color_manual(values = c('%actual_1' = love_color("shallow_red"),
                                    '%expected_1' = love_color("sky_blue"))) +
      ylim(c(-0.01, max(c(ae_total$value, as.numeric(ae_1$value) * max(ae_total$value) * 4)) + 0.1)) +
      guides(fill = guide_legend(reverse = TRUE)) +
      xlab(x) + ylab("Total Percent") +
      ggtitle("Distribution of Train/Expected and Test/Actual") +
      plot_theme(legend.position = "top",
                 title_size = 9,
                 axis_title_size = 8,
                 angle = ifelse(nrow(ae_total) > 10, 60,
                                ifelse(nrow(ae_total) > 5, 40,
                                       ifelse(nrow(ae_total) > 3, 20, 10))),
                 axis_size_x = ifelse(max(nchar(ae_total$bins)) > 30, 5,
               ifelse(max(nchar(ae_total$bins)) > 20, 6,
                      ifelse(max(nchar(ae_total$bins)) > 10, 7 ,8))))

    if (save_data) {
        ggsave(paste0(dir_path, paste(x, "png", sep = '.')),
           plot = arrangeGrob(grobs = list(plot_1, plot_2), ncol = 2, nrow = 1),
           width = g_width, height = g_width / 2, dpi = "retina")
    }
    if (plot_show) {
        plot(arrangeGrob(grobs = list(plot_1, plot_2), ncol = 2, nrow = 1))
    }
    return(df_ae)
    options(opt) # reset
}





#' ks_value
#'
#' \code{ks_value} is for get plots for a  variable.
#' @param dat_pred A data frame with predict prob or score.
#' @param target The name of target variable.
#' @param score The name of prob or score variable.
#' @param g Number of breaks for prob or score.
#' @param breaks Splitting points of prob or score.
#' @examples
#' sub = cv_split(UCICreditCard, k = 30)[[1]]
#' dat = UCICreditCard[sub,]
#' dat = re_name(dat, "default.payment.next.month", "target")
#' dat = data_cleansing(dat, target = "target", obs_id = "ID", 
#' occur_time = "apply_date", miss_values = list("", -1))
#' 
#' train_test <- train_test_split(dat, split_type = "OOT", prop = 0.7,
#'                                 occur_time = "apply_date")
#' dat_train = train_test$train
#' dat_test = train_test$test
#' x_list = c("PAY_0", "LIMIT_BAL", "PAY_AMT5", "PAY_3", "PAY_2")
#' Formula = as.formula(paste("target", paste(x_list, collapse = ' + '), sep = ' ~ '))
#' set.seed(46)
#' lr_model = glm(Formula, data = dat_train[, c("target", x_list)], family = binomial(logit))
#' 
#' dat_train$pred_LR = round(predict(lr_model, dat_train[, x_list], type = "response"), 5)
#' dat_test$pred_LR = round(predict(lr_model, dat_test[, x_list], type = "response"), 5)
#' ks_value(dat_pred = dat_train,
#' target = "target", score ="pred_LR", g = 20)
#' @import ggplot2
#' @importFrom gridExtra arrangeGrob tableGrob
#' @importFrom dplyr group_by mutate summarize  summarise n  count %>% filter mutate_if distinct ungroup
#' @importFrom data.table melt dcast fread fwrite
#' @export

ks_value <- function(dat_pred = NULL, target =NULL, score =NULL, g = 20, breaks = NULL) {

    good = bins = bad = NULL
    KS = 0 
    if (length(unlist(target)) > 1 & length(unlist(score)) > 1 && is.numeric(target) & is.numeric(score)) {
        digits_x = digits_num(score)
        opt = options(stringsAsFactors = FALSE, digits = digits_x + 1) #
        dat_pred = data.frame(target = target, score = score)
        dat_pred = checking_data(dat = dat_pred)
        if (is.null(breaks)) {
            breaks = get_breaks(dat = dat_pred, x = "score", target = "target", equal_bins = TRUE, best = FALSE, g = g, note = FALSE)
        }
        dat_pred$target = ifelse(dat_pred[, "target"] %in% list("0", "good", 0), 0, 1)
    } else {
        digits_x = digits_num(dat_pred[,score])
        opt = options(stringsAsFactors = FALSE, digits = digits_x + 1) #
        if (is.null(breaks)) {
            breaks = get_breaks(dat = dat_pred, x = score, target = target, equal_bins = TRUE, best = FALSE, g = g, note = FALSE)
        }
        dat_pred$target = ifelse(dat_pred[, target] %in% list("0", "good", 0), 0, 1)
        dat_pred$score = dat_pred[, score]
    }

    dat_pred$bins = split_bins(dat = dat_pred, x = "score", breaks = breaks, bins_no = TRUE)

    dat_sum <- dat_pred %>% filter(target %in% c(1, 0)) %>% group_by(bins) %>%
         dplyr::count(bins, target) %>% mutate(percent = n / sum(n)) %>% as.data.frame()

    dat_sum <- data.table :: dcast(dat_sum[c("bins", "target", "n")], bins ~ target, value.var = "n")
    dat_sum[is.na(dat_sum)] <- 0
    dat_sum = data.frame(unclass(dat_sum))
    if (all(names(dat_sum) == c("bins","X0", "X1"))) {
        names(dat_sum) = c("bins","good", "bad")
    } else {
        if (all(names(dat_sum) == c("bins","X1", "X0"))) {
            names(dat_sum) = c("bins","bad", "good")
        } else {
            stop(paste(target, "is neither 1 nor 0./n"))
        }
    }
    dt_ks <- transform(dat_sum,
                    cumsumG = round((cumsum(bad) / sum(bad)),4),
                     cumsumB = round((cumsum(good) / sum(good)),4)    
                     )
    KS = max(abs(round(dt_ks$cumsumB - dt_ks$cumsumG, 4)), na.rm = TRUE)
    options(opt) # reset
    return(KS)
}



#' ks_table & plot
#'
#' \code{ks_table} is for generating a model performance table.
#' \code{ks_table_plot} is for ploting the table generated by \code{ks_table}
#' \code{ks_psi_plot} is for K-S & PSI distrbution ploting.
#' @param train_pred A data frame of training with predicted prob or score.
#' @param test_pred A data frame of validation with predict prob or score.
#' @param target The name of target variable.
#' @param score The name of prob or score variable.
#' @param g Number of breaks for prob or score.
#' @param g_width Width of graphs.
#' @param plot_show Logical, show model performance in current graphic device. Default is FALSE.
#' @param gtitle The title of the graph & The name for periodically saved graphic file. Default is "_ks_psi_table".
#' @param breaks Splitting points of prob or score.
#' @param save_data Logical, save results in locally specified folder. Default is TRUE
#' @param dir_path The path for periodically saved graphic files.
#' @param file_name  The name for periodically saved data file. Default is NULL.
#' @param pos_flag The value of positive class of target variable, default: "1".
#' @examples
#' sub = cv_split(UCICreditCard, k = 30)[[1]]
#' dat = UCICreditCard[sub,]
#' dat = re_name(dat, "default.payment.next.month", "target")
#' dat = data_cleansing(dat, target = "target", obs_id = "ID", 
#' occur_time = "apply_date", miss_values = list("", -1))
#' 
#' train_test <- train_test_split(dat, split_type = "OOT", prop = 0.7,
#'                                 occur_time = "apply_date")
#' dat_train = train_test$train
#' dat_test = train_test$test
#' x_list = c("PAY_0", "LIMIT_BAL", "PAY_AMT5", "PAY_3", "PAY_2")
#' Formula = as.formula(paste("target", paste(x_list, collapse = ' + '), sep = ' ~ '))
#' set.seed(46)
#' lr_model = glm(Formula, data = dat_train[, c("target", x_list)], family = binomial(logit))
#' 
#' dat_train$pred_LR = round(predict(lr_model, dat_train[, x_list], type = "response"), 5)
#' dat_test$pred_LR = round(predict(lr_model, dat_test[, x_list], type = "response"), 5)
#' # model evaluation
#' ks_psi_plot(train_pred = dat_train, test_pred = dat_test,
#'                             score = "pred_LR", target = "target",
#'                             plot_show = TRUE)
#' tb_pred <- ks_table_plot(train_pred = dat_train, test_pred = dat_test,
#'                                         score = "pred_LR", target = "target",
#'                                      g = 10, g_width = 13, plot_show = FALSE)
#' key_index = model_key_index(tb_pred)
#' @import ggplot2
#' @importFrom gridExtra arrangeGrob tableGrob ttheme_default
#' @importFrom dplyr group_by mutate summarize  summarise n  count %>% filter
#' @importFrom data.table melt dcast
#' @export


ks_table <- function(train_pred, test_pred, target =NULL, score = "score",
                     g = 10, breaks = NULL, pos_flag = NULL) {
    opt = options(stringsAsFactors = FALSE, digits = 6) #
    `train_GB_index` = `test_GB_index` = `G` = `bins` = `B` = `%train` = `%test` = `#train` = `#test` = NULL
    if (is.null(target)) {
        stop("target is missing!\n")
    }
    if (is.null(breaks)) {
        breaks = get_breaks(dat = train_pred, x = score,
                            target = target, equal_bins = TRUE,
                            best = FALSE, g = g, note = FALSE)
    }
    train_pred$bins = split_bins(dat = train_pred, x = score, breaks = breaks, bins_no = TRUE)
    test_pred$bins = split_bins(dat = test_pred, x = score, breaks = breaks, bins_no = TRUE)
    if (!is.null(target)) {
        if (length(unique(train_pred[, target])) > 1) {
            if (length(unique(test_pred[, target])) > 1) {
                if (is.null(pos_flag)) {
                    train_pred$target = ifelse(train_pred[, target] %in% list("1", "1", "Bad", 1), "B", "G")
                    test_pred$target = ifelse(test_pred[, target] %in% list("1", "1", "Bad", 1), "B", "G")
                } else {
                    train_pred$target = ifelse(train_pred[, target] %in% pos_flag, "B", "G")
                    test_pred$target = ifelse(test_pred[, target] %in% pos_flag, "B", "G")
                    if (length(unique(train_pred$target)) == 1) {
                        stop(paste("The value in pos_flag is not one of the value of train_pred's target.\n"))
                    }
                    if (length(unique(test_pred$target)) == 1) {
                        stop(paste("The value in pos_flag is not one of the value of test_pred's target.\n"))
                    }
                }
            } else {
                stop(paste("The value of test_pred's target is unique.\n"))
            }

        } else {
            stop(paste("The value of train_pred's target is unique.\n"))
        }
    } else {
        stop(paste("The target variable is missing.\n"))
    }

    train_sum <- train_pred %>%
      dplyr::filter(train_pred$target %in% c("B", "G")) %>%
      dplyr::group_by(bins) %>%
      dplyr::count(bins, target) %>%
      dplyr::mutate(percent = n / sum(n)) %>%
      as.data.frame()
    test_sum <- test_pred %>%
      dplyr::filter(test_pred$target %in% c("B", "G")) %>%
      dplyr::group_by(bins) %>%
      dplyr::count(bins, target) %>%
      dplyr::mutate(percent = n / sum(n)) %>%
      as.data.frame()

    train_sum <- data.table :: dcast(train_sum[c("bins", "target", "n")],
                                     bins ~ target, value.var = "n")
    train_sum[is.na(train_sum)] <- 0
    test_sum <- data.table :: dcast(test_sum[c("bins", "target", "n")],
                                    bins ~ target, value.var = "n")
    test_sum[is.na(test_sum)] <- 0
    train_ks <- transform(train_sum,
                     train_total = G + B,
                     `%train_total` = round((G+B) / sum(G + B), 2),
                     `%train_B` = round(B / (G + B), 3),
                     `%train_cumG` = round(cumsum(G) / sum(G), 2),
                     `%train_cumB` = round(cumsum(B) / sum(B), 2),
                     `train_K-S` = abs(round((cumsum(B) / sum(B)) - (cumsum(G) / sum(G)), 2)))
    test_ks <- transform(test_sum,
                     test_total = G + B,
                     `%test_total` = round((G + B) / (sum(G + B)), 2),
                     `%test_B` = round(B / (G + B), 3),
                     `%test_cumG` = round(cumsum(G) / sum(G), 2),
                     `%test_cumB` = round(cumsum(B) / sum(B), 2),
                     `test_K-S` = abs(round((cumsum(B) / sum(B)) - (cumsum(G) / sum(G)), 2))
                     )

    dt_ks = merge(train_ks[c(1, 4:9)], test_ks[c(1, 4:9)])
    names(dt_ks) = c("bins", "#train", "%train", "%train_B",
                     "%train_cumG", "%train_cumB","train_K-S",
                     "#test", "%test", "%test_B", "%test_cumG",
                     "%test_cumB",  "test_K-S")
    dt_ks = dt_ks %>% dplyr::mutate(PSI = round((`%train` - `%test`) * log(`%train` / `%test`), 3),
                                    `#total` = `#train` + `#test` )
    dt_ks =  dt_ks[c("bins", "#total", "#train", "#test", "%train", "%test", "%train_B",
                     "%test_B", "%train_cumG", "%train_cumB", "%test_cumG", "%test_cumB",
                     "train_K-S", "test_K-S", "PSI")]
    options(opt) # reset
    return(dt_ks)
}
#' ks_table_plot
#'
#' @rdname ks_table
#' @export



ks_table_plot <- function(train_pred, test_pred, target = "target", score = "score",
                          g = 10,plot_show = TRUE, g_width = 12,file_name = NULL, save_data = TRUE,
                          dir_path = tempdir(), gtitle = NULL) {

    opt = options(stringsAsFactors = FALSE, digits = 6) #


    ` %train_cumG` = `%train_cumB` = `%train_B` = `%train` = `%test_cumG` = `%test_cumB` = `%test_B` =
    `%test` =  `%train_cumG` = NULL

    tb_pred = ks_table(train_pred=train_pred,test_pred = test_pred, target = target, score = score, g = g)
    total <- c("Total",
               sum(tb_pred$`#total`,na.rm = TRUE),
               sum(tb_pred$`#train`, na.rm = TRUE),
               sum(tb_pred$`#test`, na.rm = TRUE),
               as_percent(sum(tb_pred$`#train`, na.rm = TRUE) / sum(tb_pred$`#total`, na.rm = TRUE), 2),
               as_percent(sum(tb_pred$`#test`, na.rm = TRUE) / sum(tb_pred$`#total` , na.rm = TRUE), 2),
               as_percent(sum(tb_pred$`%train_B` * tb_pred$`#train`, na.rm = TRUE) / sum(tb_pred$`#train`,  na.rm = TRUE), 2),
               as_percent(sum(tb_pred$`%test_B` * tb_pred$`#test`,na.rm = TRUE) / sum(tb_pred$`#test`, na.rm = TRUE), 2), "100%", "100%",
             "100%", "100%",
             max(tb_pred$`train_K-S`,na.rm = TRUE),
             max(tb_pred$`test_K-S`, na.rm = TRUE),
             sum(tb_pred$psi, na.rm = TRUE))

    dt_pred <- tb_pred[c("bins", "#total", "#train", "#test", "%train", "%test", "%train_B",
                         "%test_B", "%train_cumG", "%train_cumB", "%test_cumG",
                         "%test_cumB", "train_K-S", "test_K-S", "PSI")]
    dt_pred <- transform(dt_pred,
                         `%train` = as_percent(`%train`, digits = 2),
                         `%test` = as_percent(`%test`, digits = 2),
                         `%train_B` = as_percent(`%train_B`, digits = 3),
                         `%test_B` = as_percent(`%test_B`, digits = 3),
                         `%train_cumG` = as_percent(`%train_cumG`, digits = 2),
                         `%train_cumB` = as_percent(`%train_cumB`, digits = 2),
                         `%test_cumG` = as_percent(`%test_cumG`, digits = 2),
                         `%test_cumB` = as_percent(`%test_cumB`, digits = 2)
                         )

   dt_pred <- rbind(dt_pred, total)
    names(dt_pred) = c("bins", "#total", "#train", "#test", "%train", "%test", "%train_B",
                       "%test_B", "%train_cumG", "%train_cumB", "%test_cumG",
                       "%test_cumB", "train_K-S", "test_K-S", "PSI")
    if (save_data) {
        dir_path = ifelse( !is.character(dir_path),
                      tempdir(), dir_path)
        if (!dir.exists(dir_path)) dir.create(dir_path)
        if (!is.character(gtitle)) { gtitle = paste0("mymodel") }
        tb_ks = tableGrob(dt_pred, rows = NULL,
                     cols = c("bins", "#total", "#train", "#test",
                                  "%train", "%test", "%train_B", "%test_B",
                                  "%train_cumG", "%train_cumB", "%test_cumG",
                                  "%test_cumB", "train_K-S", "test_K-S", "PSI"),
                     theme = ttheme_default(base_size = 8, base_colour = "black",
                                            base_family = "", parse = FALSE,
                                            padding = unit(c(3, 3), "mm")
                                            )
                     )
        ggsave(paste0(dir_path, "/", paste(paste(gtitle, "_ks_psi_table"), "png", sep = '.')),
            arrangeGrob(grobs = list(tb_ks), ncol = 1, nrow = 1), dpi = "retina", width = g_width)
        save_dt(dt_pred, file_name = paste(gtitle, "_ks_psi_table"), dir_path = dir_path)
        }
    return(dt_pred)
    options(opt) # reset
}


#' ks_table_plot
#'
#' @rdname ks_table
#' @export

ks_psi_plot <- function(train_pred, test_pred, target = "target", score = "score",
                        gtitle = NULL, plot_show = TRUE, g_width = 12, save_data = TRUE,
                        breaks = NULL, g = 10, dir_path = tempdir()) {
    opt = options(stringsAsFactors = FALSE, digits = 6) #
   `value` = `train_test` = `bins` =   `%train_cumG` =   `%train_cumB` =  `%test_cumG` =  `%test_cumB` = NULL
    if (!dir.exists(dir_path)) dir.create(dir_path)
    if (is.null(gtitle)) { gtitle = paste0("Model") }
    tb_ks = ks_table(train_pred = train_pred, test_pred = test_pred,
                     target = target, score = score, g = g, breaks = breaks)

    ks <- rbind(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0,0), tb_ks)
    ks_plot <- ggplot(ks, aes(x = reorder(as.factor(as_percent(round(cumsum(ks$`%train`), 1))), cumsum(ks$`%train`)))) +
      geom_line(aes(y = `%train_cumG`, group = 1, color = "%train_cumG"), size = 1) +
      geom_point(aes(x = which.max(ks$`train_K-S`),
                     y = as.numeric(ks$`%train_cumG`[which.max(ks$`train_K-S`)])),
                 size = 2, shape = 21, fill = 'white', color = "#085A9C") +
      geom_segment(aes(x = which.max(ks$`train_K-S`),
                       y = as.numeric(ks$`%train_cumG`[which.max(ks$`train_K-S`)]) + 0.01,
                       xend = which.max(ks$`train_K-S`),
                       yend = as.numeric(ks$`%train_cumB`[which.max(ks$`train_K-S`)]) - 0.01),
                   colour = love_color("deep_grey"), linetype = "dashed",
                   arrow = arrow(ends = "both", length = unit(0.2, "cm"))) +
      geom_line(aes(y = `%train_cumB`, color = '%train_cumB', group = 1), size = 1) +
      geom_point(aes(x = which.max(ks$`train_K-S`),
                     y = as.numeric(ks$`%train_cumG`[which.max(ks$`train_K-S`)])),
                 size = 2, shape = 21, fill = 'white', color = '#ca3e1c') +
      geom_point(aes(x = which.max(ks$`train_K-S`),
                     y = as.numeric(ks$`%train_cumB`[which.max(ks$`train_K-S`)])),
                 size = 2, shape = 21, fill = 'white', color = '#ca3e1c') +
      annotate(geom = 'text', x = 7, y = 0.1,
               label = paste('train K-S : ', as_percent(max(round(ks$`train_K-S`, 2)), digits = 2)),
                vjust = 1.5) +
      geom_line(aes(y = `%test_cumG`, group = 1, color = "%test_cumG"),
                linetype = "dashed", size = 1) +
      geom_point(aes(x = which.max(ks$`test_K-S`),
                     y = as.numeric(ks$`%test_cumG`[which.max(ks$`test_K-S`)])),
                 size = 2, shape = 21, fill = 'white', color = "#085A9C") +
      geom_segment(aes(x = which.max(ks$`test_K-S`),
                       y = as.numeric(ks$`%test_cumG`[which.max(ks$`test_K-S`)]) + 0.01,
                       xend = which.max(ks$`test_K-S`),
                       yend = as.numeric(ks$`%test_cumB`[which.max(ks$`test_K-S`)]) - 0.01),
                   colour = love_color("deep_grey"), linetype = "dashed",
                   arrow = arrow(ends = "both", length = unit(0.2, "cm"))) +
      geom_line(aes(y = `%test_cumB`, color = '%test_cumB', group = 1),
                linetype = "dashed", size = 1) +
      geom_point(aes(x = which.max(ks$`test_K-S`),
                     y = as.numeric(ks$`%test_cumG`[which.max(ks$`test_K-S`)])),
                 size = 2, shape = 21, fill = 'white', color = '#ca3e1c') +
      geom_point(aes(x = which.max(ks$`test_K-S`),
                     y = as.numeric(ks$`%test_cumB`[which.max(ks$`test_K-S`)])),
                 size = 2, shape = 21, fill = 'white', color = '#ca3e1c') +
      annotate(geom = 'text', x = 7, y = 0.15,vjust = 1.5,
               label = paste('test K-S : ', as_percent(max(round(ks$`test_K-S`, 2)), digits = 2))) +
      annotate("text", vjust = -1, label = "1", size = 4, colour = "black",
               x = which.max(ks$`train_K-S`) - 1,
               y = as.numeric(ks$`%train_cumB`[which.max(ks$`train_K-S`)])) +
      annotate("text", vjust = 1, label = "0", size = 4, colour = "black",
               x = which.max(ks$`test_K-S`) + 1,
               y = as.numeric(ks$`%test_cumG`[which.max(ks$`test_K-S`)])) +
      scale_colour_manual(values = c("%train_cumG" = love_color("dark_blue"),
                                     "%test_cumG" = love_color("dark_green"),
                                     "%train_cumB" = love_color("dark_red2"),
                                     "%test_cumB" = love_color("dark_purple"))) +
      labs(x = "% of Total", y = "% of CumSum G/B",
           title = paste(gtitle,"K-S : Train vs. Test" ))+
      plot_theme(legend.position = "top", angle = 0)

    ts_total <- data.table::melt(tb_ks[c("bins", "%train", "%test")], id.vars = c("bins"),
                                 variable.name = "train_test", value.name = "value")
    ts_1 <- data.table::melt(tb_ks[c("bins", "%train_B", "%test_B")], id.vars = c("bins"),
                               variable.name = "train_test", value.name = "value")

    psi_plot <- ggplot(ts_total, aes(x = bins, y = value, fill = train_test)) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.7)) +
      geom_text(aes(y = value, label = paste(as_percent(value, digits = 3))),
                position = position_dodge(width = 0.7),
                size = 3, vjust = 1, hjust = 0.3, colour = "white") +
      geom_line(aes(x = factor(ts_1[[1]]),
                    y = as.numeric(ts_1$value) * max(ts_total$value) * 4,
                    color = ts_1$train_test,
                    linetype = ts_1$train_test,
                    group = ts_1$train_test),
                position = position_dodge(width = 0.5),size = 1) +
      geom_point(aes(y = as.numeric(ts_1$value) * max(ts_total$value) * 4,
                     color = ts_1$train_test,
                     group = ts_1$train_test),
                 position = position_dodge(width = 0.5),
                 fill = 'white', size = 2, shape = 21) +
      geom_text(aes(y = as.numeric(ts_1$value) * max(ts_total$value) * 4,
                    label = paste(as_percent(ts_1$value, digits = 3))),
                position = position_dodge(width = 0.5),
                colour = 'black', size = 3, vjust = -0.1) +
      annotate(geom = 'text',
             x = dim(ts_total)[1] / 3,
             y = max(c(ts_total$value, as.numeric(ts_1$value) * max(ts_total$value) * 4)) + 0.09,
             label = paste('PSI', sum(tb_ks$psi), sep = " : ")) +
      scale_fill_manual(values = c('%train' = love_color("deep_grey"),
                                   '%test' = love_color("light_yellow"))) +
      scale_color_manual(values = c('%train_B' = love_color("shallow_red"),
                                    '%test_B' = love_color("sky_blue"))) +
      ylim(c(-0.01, max(c(ts_total$value, as.numeric(ts_1$value) * max(ts_total$value) * 4)) + 0.1)) +
      guides(fill = guide_legend(reverse = TRUE)) +
      xlab(score) +
      ylab("Total Percent") +
      ggtitle(paste(gtitle,"Train and Test Distribution"))+
      plot_theme(legend.position = "top",
                 angle = ifelse(nrow(ts_total) > 10, 50,
                                ifelse(nrow(ts_total) > 5, 30, 0)))
    if (save_data) {
        dir_path = ifelse(is.null(dir_path),
                      tempdir(), dir_path)
        if (!dir.exists(dir_path)) dir.create(dir_path)
        if ( !is.character(gtitle)) { gtitle = paste0("mymodel") }
        ggsave(filename = paste0(dir_path, "/", paste(paste0(gtitle, "_ks_psi_plot"), "png", sep = '.')),
           device = "png",
           plot = arrangeGrob(grobs = list(ks_plot, psi_plot), ncol = 2, nrow = 1),
           dpi = "retina", width = g_width, height = g_width / 2)
    }
    if (plot_show) {
        ks_psi_plot = arrangeGrob(grobs = list(ks_plot, psi_plot), ncol = 2, nrow = 1)
        return(plot(ks_psi_plot))
    }
    options(opt) # reset
}

#' Correlation Plot
#'
#' \code{cor_plot} is for ploting correlation matrix
#' @param dat A data.frame with independent variables and target variable.
#' @param x_list Names of independent variables.
#' @param plot_show  Logical, show graph in current graphic device.
#' @param save_data Logical, save results in locally specified folder. Default is TRUE
#' @param gtitle The title of the graph & The name for periodically saved graphic file. Default is "_correlation_of_variables".
#' @param dir_path The path for periodically saved graphic files. Default is "./model/LR"
#' @examples
#' train_test <- train_test_split(UCICreditCard,
#' split_type = "Random", prop = 0.8,save_data = TRUE)
#' dat_train = train_test$train
#' dat_test = train_test$test
#' cor_plot(dat_train[,8:12],plot_show = TRUE)
#' @import ggplot2
#' @importFrom gridExtra arrangeGrob tableGrob
#' @importFrom ggcorrplot ggcorrplot
#' @export

cor_plot = function(dat, dir_path = tempdir(), x_list = NULL,
                    gtitle = NULL, save_data = TRUE, plot_show = FALSE) {
    if (!is.null(x_list)) {
        dat = dat[, c(x_list)]
    }
    num_x_list = get_names(dat = dat,
                           types = c('numeric', 'integer', 'double'),
                           ex_cols = "",
                           get_ex = FALSE)
    cor_mat <- cor(dat[, num_x_list],use = "complete.obs")
    if (save_data) {
        save_dt(cor_mat, file_name = paste(gtitle, "correlation_matrix"),
                dir_path = dir_path, note = FALSE, row_names = TRUE)
    }
    cor_p = ggcorrplot(cor_mat, hc.order = TRUE, lab = TRUE,
                       outline.color = "white",
                       ggtheme = plot_theme(),
                       show.legend = FALSE,
                       title = paste(gtitle, "Correlation  of  Variables"),
                       colors = c(love_color("deep_red"), "white", love_color("dark_green")),
                       lab_col = "black",
                       lab_size = ifelse(nrow(cor_mat) < 10, 4,
                                         ifelse(nrow(cor_mat) > 20, 2, 3)),
                       pch = ifelse(nrow(cor_mat) < 10, 4,
                                    ifelse(nrow(cor_mat) > 20, 2, 3)), pch.col = "black",
                       pch.cex = ifelse(nrow(cor_mat) < 10, 4,
                                        ifelse(nrow(cor_mat) > 20, 2, 3)),
                       tl.cex = ifelse(nrow(cor_mat) < 10, 12,
                                       ifelse(nrow(cor_mat) > 20, 8, 10)), tl.col = "black",
                       tl.srt = ifelse(nrow(cor_mat) < 10, 20,
                                       ifelse(nrow(cor_mat) > 20, 45, 30)))
    if (save_data) {
        dir_path = ifelse(!is.character(dir_path) ,
                      tempdir(), dir_path)
        if (!dir.exists(dir_path)) dir.create(dir_path)
        if ( !is.character(gtitle)) { gtitle = paste0("mymodel") }
        ggsave(paste0(dir_path, "/", paste(paste0(gtitle, "_correlation_of_variables"), "png", sep = '.')), cor_p, dpi = "retina", width = 8)
    }
    if (plot_show) {
        plot(cor_p)
    }
}
#' model_key_index
#' \code{model_key_index} is for get plots for a  variable.
#' @param tb_pred  A table generated by code{\link{ks_table}}
#' @rdname ks_table
#' @export
model_key_index <- function(tb_pred) {
    key_index = NULL
    if (all(names(tb_pred) %in% c("train_K-S", "test_K-S", "PSI"))) {
        b_psi = as.numeric(tb_pred[nrow(tb_pred), "PSI"])
        train_KS = as.numeric(tb_pred[nrow(tb_pred), "train_K-S"])
        test_KS = as.numeric(tb_pred[nrow(tb_pred), "test_K-S"])
        key_index = data.frame(train_KS = train_KS, test_KS = test_KS, PSI = b_psi)
    } else {
        key_index = data.frame(train_KS = NA, test_KS = NA, PSI = NA)
    }
    return(key_index)
}


#' love_color
#'
#' \code{love_color} is for get plots for a  variable.
#' @param  color The name of colors.
#' @examples
#' love_color("dark_cyan")
#' @import ggplot2
#' @importFrom gridExtra arrangeGrob tableGrob
#' @export

love_color <- function(color = "dark_cyan") {
    color_board <- c(
    dark_cyan = rgb(20, 78, 100, maxColorValue = 255),
    deep_cyan = rgb(28, 82, 110, maxColorValue = 255),
    light_cyan = rgb(0, 93, 125, maxColorValue = 255),
    shallow_cyan = rgb(34, 116, 135, maxColorValue = 255),
     green_cyan = rgb(0, 139, 197, maxColorValue = 255),
    dark_green = rgb(43, 149, 136, maxColorValue = 255),
    deep_green = rgb(84, 139, 84, maxColorValue = 255),
    light_green = rgb(50, 150, 78, maxColorValue = 255),
    shallow_green = rgb(95, 160, 78, maxColorValue = 255),
    pale_green = rgb(180, 202, 198, maxColorValue = 255),
    dark_grey = rgb(102, 102, 102, maxColorValue = 255),
    deep_grey = rgb(128, 129, 128, maxColorValue = 255),
    light_grey = rgb(191, 192, 191, maxColorValue = 255),
    shallow_grey = rgb(169, 169, 169, maxColorValue = 255),
    pale_grey = rgb(240, 240, 240, maxColorValue = 255),
    dark_red = rgb(154, 42, 42, maxColorValue = 255),
    deep_red = rgb(174, 61, 63, maxColorValue = 255),
    light_red = rgb(202, 62, 28, maxColorValue = 255),
    shallow_red = rgb(252, 74, 42, maxColorValue = 255),
    dark_blue = rgb(33, 71, 117, maxColorValue = 255),
    deep_blue = rgb(32, 81, 139, maxColorValue = 255),
    light_blue = rgb(22, 83, 161, maxColorValue = 255),
    pale_blue = rgb(144, 190, 216, maxColorValue = 255),
    sky_blue = rgb(0, 142, 213, maxColorValue = 255),
    water_blue = rgb(60, 140, 190, maxColorValue = 255),
    grey_blue = rgb(80, 99, 139, maxColorValue = 255),
    dark_purple = rgb(120, 78, 100, maxColorValue = 255),
    deep_purple = rgb(170, 44, 105, maxColorValue = 255),
    light_purple = rgb(212, 137, 168, maxColorValue = 255),
    deep_orange = rgb(255, 140, 0, maxColorValue = 255),
    gold_yellow = rgb(255, 215, 0, maxColorValue = 255),
    light_yellow = rgb(207, 177, 81, maxColorValue = 255),
    dark_purple = rgb(71, 0, 123, maxColorValue = 255),
    light_orange = rgb(241, 156, 0, maxColorValue = 255),
   dark_red2 = rgb(181, 0, 0, maxColorValue = 255),
  shallow_orange = rgb(241, 121, 0, maxColorValue = 255),
  lake_blue = rgb(0, 91, 181, maxColorValue = 255)
    )
    return(color_board[[color]])
}
#' plot_theme
#'
#' \code{plot_theme} is a simper wrapper of theme for ggplot2.
#' @param legend.position see details at: code{legend.position}
#' @param angle see details at:  code{axis.text.x}
#' @param legend_size  see details at:  code{legend.text}
#' @param axis_size_x see details at:  code{axis.text.x}
#' @param axis_size_y see details at:  code{axis.text.y}
#' @param axis_title_size see details at:  code{axis.title.x}
#' @param title_size see details at:  code{plot.title}
#' @param title_vjust see details at:  code{plot.title}
#' @param title_hjust see details at:  code{plot.title}
#' @param linetype see details at:  code{panel.grid.major}
#' @param face see details at:  code{axis.title.x}
#' @details see details at: code{theme}
#' @import ggplot2
#' @importFrom gridExtra arrangeGrob
#' @export

plot_theme <- function(legend.position = "top", angle = 30,
                       legend_size = 7, axis_size_y = 8,
                       axis_size_x = 8, axis_title_size = 10,
                       title_size = 11, title_vjust = 0, title_hjust = 0,
                       linetype = "dotted", face = "bold") {
    plot_theme <- theme_bw() + theme(legend.position = legend.position,
                              panel.border = element_blank(),
                              panel.grid.major = element_line(linetype = linetype),
                              panel.grid.minor = element_blank(),
                              plot.title = element_text(face = "bold", size = title_size,
                                                        vjust = title_vjust, hjust = title_hjust),
                              legend.title = element_blank(),
                              legend.text = element_text(size = legend_size),
                              legend.key = element_blank(),
                              axis.text.x = element_text(size = axis_size_x,
                                                         vjust = 0.5, angle = angle),
                              axis.text.y = element_text(size = axis_size_y),
                              axis.title.x = element_text(size = axis_title_size, face = face),
                              axis.title.y = element_text(size = axis_title_size),
                              strip.text = element_text(size = 10, face = face),
                              strip.background = element_blank())
    return(plot_theme)
}

