#'Checking Data
#'
#'
#' \code{checking_data}  cheking dat before processing.
#' @param dat A data.frame with independent variables and target variable.
#' @param target The name of target variable. Default is NULL.
#' @param occur_time The name of the variable that represents the time at which each observation takes place.
#' @param pos_flag The value of positive class of target variable, default: "1".
#' @param note Logical.Outputs info.Default is TRUE.
#' @return data.frame
#' @examples
#' \dontrun{
#' # load germancredit dat
#' data(UCI_Credit_Card)
#' dat = checking_data(dat = UCI_Credit_Card)
#' }
#' @export



checking_data <- function(dat = NULL, target = NULL, occur_time = NULL, note = FALSE, pos_flag = NULL) {
    if (note) {
        (cat("[NOTE] checking dat and target format.\n"))
    }
    if (is.null(dat)) {
        warning("dat is null.\n")
    } else {
        if (!(class(dat)[1] == "data.frame")) {
            if (class(dat)[1] %in% c("data.table", "list", "tbl_df", "tbl") && length(dim(dat)) == 2) {
                dat = as.data.frame(dat)
                cat(paste("[NOTE]", "convert", class(dat)[1], "to data.frame.\n"))
            } else {
                warning("dat is neither a dat frame nor a dat table.\n")
            }
        }
    }
    if (!is.null(target)) {
        if (length(unique(dat[, target])) < 2) {
            warning(paste("Unique values of ", target, "is", "only one.\n", sep = "\t"))
        } else {
            if (length(unique(dat[, target])) == 2) {

                if (is.null(pos_flag)) {
                    pos_flag = list("1", "bad", 1, "B", "positive", "pos", "Positive", "Pos")
                }
                if (length(which(dat[, target] %in% pos_flag)) != 0) {
                    if (!is.numeric(dat[, target]) || is.numeric(dat[, target]) && !all(sort(unique(dat[, target])) == c(0, 1))) {
                        pos <- unique(dat[, target])[which(unique(dat[, target]) %in% pos_flag)]
                        dat[, target] = ifelse(dat[, target] %in% pos_flag, 1, 0)
                        warning(paste("The  values in of target has been encoded", pos, "=1 and others = 0", " \n"))
                    }
                } else {
                    warning(paste("Positive values of", target, "is not in pos_flag:", paste(pos_flag, collapse = ","), "\nplease set pos_flag. \n", sep = "\t"))
                }
            }
        }
    }    
    if (!is.null(occur_time)) {
        if (is.element(occur_time, names(dat))) {
            dat <- time_transfer(dat, date_cols = occur_time)
            if (!is_date(dat[, occur_time])) {
                warning(paste("occur_time:", occur_time, "is not time or date.\n", sep = "\t"))
            }
        } else {
            warning(paste("occur_time:", occur_time, "is not in data.\n", sep = "\t"))
        }
    }
    return(dat)
}
