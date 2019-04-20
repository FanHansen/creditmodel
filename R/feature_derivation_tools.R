#' derived_interval
#'
#' This function is not intended to be used by end user. 
#'
#' @param  dat_s  A data.frame contained only predict variables.
#' @param interval_type  Available of c("cnt_interval", "time_interval")
#' @export
#' @importFrom data.table first
derived_interval <- function(dat_s, interval_type = c("cnt_interval", "time_interval")) {
    interval_list <- apply(dat_s, 1, function(m) {
        if (interval_type == "time_interval") {
            cnt_ind = inter_ind = which(!is.na(m) | m != 0)
        } else {
            cnt_ind = which(m >= 0)
            inter_ind = unlist(m, use.names = FALSE)[c(cnt_ind)]
        }
        interval <- rep(NA, length(m))
        if (length(cnt_ind) > 1) {
            interval[cnt_ind] <- vapply(1:(length(inter_ind)), function(i) {

                ifelse(i <= length(inter_ind), abs(inter_ind[i] - inter_ind[i +1]), NA)

            }, FUN.VALUE = numeric(1))
        }
        interval = c(abs(1 - data.table::first(inter_ind)), interval[-length(interval)])
        interval
    })
    interval_list = as.data.frame(t(interval_list))
    interval_list
}


#' derived_pct
#'
#' This function is not intended to be used by end user. 
#'
#' @param  dat_s  A data.frame contained only predict variables.
#' @param pct_type  Available of "total_pct"
#' @export
derived_pct <- function(dat_s, pct_type = "total_pct") {
    dat_s[is.na(dat_s)] <- 0
    if (pct_type == "total_pct") {
        pct_list =  dat_s / rowSums(dat_s, na.rm = TRUE)    
    } else {
        cnt_pct_list= dat_s / rowSums(dat_s, na.rm = TRUE)
        pct_list = apply(cnt_pct_list, 1, function(x) cumsum(x))
        pct_list = as.data.frame(t(pct_list))
    }
   
    pct_list
}


#' derived_partial_acf
#'
#' This function is not intended to be used by end user. 
#'
#' @param  dat_s  A data.frame
#' @export

derived_partial_acf <- function(dat_s) {
    dat_s[is.na(dat_s)] <- 0
    p_acf <- apply(dat_s, 1, function(x) ifelse(length(unique(x)) > 2, mean(abs(ar(ts(x), FALSE,
    length(unique(x)) - 1, na.action = na.pass)$partialacf)), NA))
    p_acf
}

#' sim_str
#'
#' This function is not intended to be used by end user. 
#'
#' @param a A string
#' @param b  A string
#' @param sep Seprater of strings. Default is "_|[.]|[A-Z]".
#' @export
sim_str <- function(a, b, sep = "_|[.]|[A-Z]") {
    intersect(strsplit(a, sep)[[1]], strsplit(b, sep)[[1]])
}


