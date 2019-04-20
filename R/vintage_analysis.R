#' vintage_function
#' \code{vintage_function} is for vintage analysis.
#'
#' This function is not intended to be used by end user. 
#'
#' @param  dat  A data.frame contained id, occur_time, mob, status ...
#' @param obs_id  The name of ID of observations or key variable of data. Default is NULL.
#' @param occur_time The name of the variable that represents the time at which each observation takes place.
#' @param MOB  Mobility of book
#' @param period Period of event to analysiss. Default is "monthly"
#' @param status Status of observations 
#' @param amount  The name of variable representing amount. Default is NULL.
#' @param by_out  Output: amount (amt) or count (cnt)
#' @param start_date The earliest occurrence time of observations.
#' @param end_date The latest occurrence time of observations.
#' @param dead_status Status of dead observations.
#' @importFrom dplyr group_by mutate summarize  summarise n  count %>% filter left_join
#' @importFrom data.table dcast melt fread fwrite
#' @export

vintage_function <- function(dat, obs_id = NULL, occur_time = NULL, MOB = NULL,
period = "monthly", status = NULL, amount = NULL, by_out = "cnt",
start_date = NULL, end_date = NULL, dead_status = 30) {
    dat = checking_data(dat = dat, occur_time = occur_time)
    dat = time_transfer(dat, date_cols = occur_time)
    ID = aging = cohort = sta = amt =  NULL
    if (!is.null(occur_time)&&is.null(start_date)) {
        start_date = date_cut(dat_time = dat[, occur_time], pct = 0)
   }
    if (!is.null(occur_time) && is.null(end_date)) {
        end_date = date_cut(dat_time = dat[, occur_time], pct = 1)
    }
    dat = dat[which(dat[, occur_time] >= start_date & dat[, occur_time] <= end_date),]
    if (!is.null(MOB)) {
        dat$aging = as.numeric(dat[, MOB])
    } else {
        if (!is.null(obs_id)) {
            dat$ID = dat[, obs_id]
            dat = dat %>% dplyr::group_by(ID) %>% mutate(aging = seq(1, length(ID)) - 1) %>% ungroup()
        } else {
            stop("MOB & obs_id  are both missing.\n")
        }     
    }
    if (!is.null(status)){
        dat$sta = dat[, status]
    }
    dat <- dat %>% filter(aging != 'NA')
    if (period == "weekly") { 
        dat$cohort <- cut(dat[, occur_time], breaks = "week")
    } else {
        dat$cohort = cut(dat[, occur_time], breaks = "month")     
    } 
    if (by_out == "amt" && !is.null(amount)) {
        #amount
        dat$amt = as.numeric(dat[, amount])
        dat1 = dat %>% dplyr::group_by(cohort, aging) %>% dplyr::summarise(n = round(sum(amt, na.rm = TRUE), 0))
        dat1 = data.table::dcast(dat1, cohort ~ aging, value.var = "n")
        if (length(unique(dat$sta)) > 1) {
            if (!is.null(dead_status)) {
                if (is.numeric(dead_status)) {
                    dat$sta = as.numeric(dat$sta)
                    dat2 = dat %>% dplyr::filter(sta > dead_status) %>% dplyr::group_by(cohort, aging) %>% dplyr::summarise(n = round(sum(amt, na.rm = TRUE), 0))
                } else {
                    if (is.character(dead_status)) {
                        if (is.element(dead_status, unique(dat$sta))) {
                            dat2 = dat %>% dplyr::filter(status == dead_status) %>% dplyr::group_by(cohort, aging) %>% dplyr::summarise(n = round(sum(amt, na.rm = TRUE), 0))
                        } else {
                            stop("dead_status is not one of status.\n")
                        }
                    } else {
                        stop("dead_status is neither numeric nor character.\n")
                    }
                }
            } else {
                stop("dead_status is missing.\n")
            }
        } else {
            stop("the unique value of status is less than one.\n")
        }
     
        dat2 = data.table::dcast(dat2, cohort ~ aging, value.var = "n")
    } else {
        #count
        dat1 = dat %>% dplyr:: group_by(cohort, aging) %>% dplyr::count(cohort, aging)
        dat1 = data.table::dcast(dat1, cohort ~ aging, value.var = "n")
        if (length(unique(dat$sta)) > 1) {
            if (!is.null(dead_status)) {
                if (is.numeric(dead_status)) {
                    dat$sta = as.numeric(dat$sta)
                    dat2 = dat %>% dplyr::filter(sta > dead_status) %>% dplyr::group_by(cohort, aging) %>% dplyr::count(cohort, aging)
                } else {
                    if (is.character(dead_status)) {
                        if (is.element(dead_status, unique(dat$sta))) {
                            dat2 = dat %>% dplyr::filter(sta == dead_status) %>% dplyr::group_by(cohort, aging) %>% dplyr::count(cohort, aging)
                        } else {
                            stop("dead_status is not one of status.\n")
                        }                   
                    } else {
                        stop("dead_status is neither numeric nor character.\n")
                    }
            }
            } else {
                stop("dead_status is missing.\n")
            }
        } else {
            stop("the unique value of status is less than one.\n")
        }     
        dat2 = data.table::dcast(dat2, cohort ~ aging, value.var = "n")
    }
    dat3 = dplyr::left_join(dat1, dat2, by = "cohort")
    M = colnames(dat1[3:length(dat1)])
    N = colnames(dat2[2:length(dat2)])
    A = matrix(dat3[3:length(dat1)])
    A = sapply(A, function(x) gsub("^0$", NA, x))    
    B = matrix(dat3[(length(dat1) + 1):length(dat3)])
    rownames(B) = N
    MN <- matrix(lapply(1:length(setdiff(M, N)), function(i) runif(dim(dat3)[1], 0, 0)))
    rownames(MN) = setdiff(M, N)
    B = rbind(MN, B)
    B = as.matrix(B[order(as.numeric(rownames(B))),])
    B = sapply(B, function(x) ifelse(is.na(x),0,x))
    dat4 = data.frame(B / as.numeric(A))
    dat4[dim(dat4)[1], 1] = ""
    dat4 = data.frame(lapply(dat4, function(x) as_percent(x, 4)), stringsAsFactors = FALSE)
    dat5 = data.frame(dat3[1:2], dat4)
    dat5 = data.frame(lapply(dat5, function(x) sub("NA|NaN|Inf", "", as.character(x))), stringsAsFactors = FALSE)
    dat5[is.na(dat5)] = ""
    if (by_out == "amt") {
        for (i in 1:(length(dat5) - 2)) {
            names(dat5) = c("cohort", "total_amt", 1:i)
        }
    } else {
        for (i in 1:(length(dat5) - 2)) {
            names(dat5) = c("cohort", "total_cnt", 1:i)
        }
    }
    return(dat5)
}



