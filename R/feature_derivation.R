
#' Derivation of Behavioral Variables
#'
#' This function is used for derivating behavioral variables . 
#'
#' @param  dat  A data.frame contained only predict variables.
#' @param  grx  Regular expressions used to match variable names.
#' @param  grx_x  Regular expression used to match a group of variable names.
#' @param  td  Number of variables to derivate.
#' @param  der  Variables to derivate
#' @param parallel Logical, parallel computing. Default is FALSE.
#' @importFrom data.table setDT :=  rbindlist 
#' @export

derived_ts_vars <- function(dat, grx, td = 12,
der = c("cvs", "sums", "means", "maxs", "max_mins", "time_intervals", "cnt_intervals", "total_pcts", "cum_pcts", "partial_acfs"),
parallel = TRUE) {
    cat(paste("derived variables of", paste(der),". \n"))
    if (parallel) {
        parallel <- start_parallel_computing(parallel)
        stopCluster <- TRUE
    } else {
        parallel <- stopCluster <- FALSE
    }
    on.exit(if (parallel & stopCluster) stop_parallel_computing(attr(parallel, "cluster")))
    i. = NULL
    if (!parallel) {
        df_cv_list <- lapply(unlist(grx), function(grx) derived_ts(dat, grx, td = td, der = der))
        df_cv_list <- as.data.frame(Reduce("cbind", df_cv_list))
    } else {
        df_cv_list <- foreach(i. = unlist(grx), .combine = "c", .errorhandling = c('pass')) %dopar% {
            try(do.call(derived_ts, args = list(dat = dat, grx_x = i., td = td,
            der = der)), silent = TRUE)
        }
        df_cv_list <- as.data.frame(df_cv_list)
    }
    return(df_cv_list)
}

#' @rdname derived_ts_vars
#' @importFrom stringr str_extract
#' @export


derived_ts <- function(dat = dat, grx_x = NULL, td = 12, der = c("cvs", "sums", "means", "maxs", "max_mins","time_intervals","cnt_intervals","total_pcts","cum_pcts","partial_acfs")) {
    setDT(dat)
    cv_cols <- grep(grx_x, paste(colnames(dat)))[1:td]
    cv_cols <- cv_cols[!is.na(cv_cols)]
    #cv_folds
    if (length(cv_cols) > 0) {
        name_n = orignal_nam = sim_nam =  str_num = c()
        orignal_nam <- names(dat[, cv_cols, with = FALSE])
        str_num = as.numeric(str_extract(orignal_nam, "\\d+"))
        if (!any(is.na(str_num)) && length(str_num) == td) {
            name_n = paste(min(str_num), max(str_num), sep = "to")
        }
        sim_nam = paste(unique(lapply(1:(length(orignal_nam) - 1),
        function(x) sim_str(orignal_nam[x], orignal_nam[x + 1])))[[1]], collapse = "_")
        if (any(der == "cvs")) {
            dat = dat[, paste(sim_nam, name_n, "_cvs", sep = "") := ifelse(rowAllnas(dat[, cv_cols, with = FALSE]), NA,
            rowCVs(dat[, cv_cols, with = FALSE], na.rm = TRUE))]
        }
        if (any(der == "sums")) {
            dat = dat[, paste(sim_nam, name_n, "sums", sep = "_") := ifelse(rowAllnas(dat[, cv_cols, with = FALSE]), NA,
            rowSums(dat[, cv_cols, with = FALSE], na.rm = TRUE))]
        }
        if (any(der == "means")) {
            dat = dat[, paste(sim_nam, name_n, "means", sep = "_") := ifelse(rowAllnas(dat[, cv_cols, with = FALSE]), NA,
            rowMeans(dat[, cv_cols, with = FALSE], na.rm = TRUE))]
        }
        if (any(der == "maxs")) {
            dat = dat[, paste(sim_nam, name_n, "maxs", sep = "_") := ifelse(rowAllnas(dat[, cv_cols, with = FALSE]), NA,
            rowMaxs(dat[, cv_cols, with = FALSE]))]
        }

        if (any(der == "max_mins")) {
            dat = dat[, paste(sim_nam, name_n, "max_mins", sep = "_") := ifelse(rowAllnas(dat[, cv_cols, with = FALSE]), NA,
            rowMaxMins(dat[, cv_cols, with = FALSE], na.rm = TRUE))]
        }
        if (any(der == "partial_acfs")) {
            dat = dat[, paste(sim_nam, name_n, "partial_acfs", sep = "_") := derived_partial_acf(dat[, cv_cols, with = FALSE])]
        }
        if (any(der == "time_intervals")) {
            dat = dat[, paste(orignal_nam, "time_intervals", sep = "_") := derived_interval(dat[, cv_cols, with = FALSE],
            interval_type = "time_interval")]
        }
        if (any(der == "cnt_intervals")) {
            dat = dat[, paste(orignal_nam, "cnt_intervals", sep = "_") := derived_interval(dat[, cv_cols, with = FALSE],
            interval_type = "cnt_interval")]
        }
        if (any(der == "total_pcts")) {
            dat = dat[, paste(orignal_nam, "total_pcts", sep = "_") := derived_pct(dat[, cv_cols, with = FALSE],
            pct_type = "total_pct")]
        }
        if (any(der == "cum_pcts")) {
            dat = dat[, paste(orignal_nam, "cum_pcts", sep = "_") := derived_pct(dat[, cv_cols, with = FALSE],
            pct_type = "cum_pct")]
        }   
    }
    return(dat)
}



#' Processing of Time or Date Variables
#'
#' This function is not intended to be used by end user. 
#'
#' @param  df_tm  A data.frame 
#' @param  x  Time variable.
#' @param  enddate  End time.
#' @export
#' @importFrom data.table  hour setnames
time_vars_process <- function(df_tm = df_tm, x, enddate = "occur_time") {
    if ((class(df_tm[[x]])[1] == 'POSIXct' | class(df_tm[[x]])[1] == 'Date') & x != enddate) {
        mydata <- within(df_tm, {
            new = as.numeric(df_tm[[enddate]] - as.Date(df_tm[[x]]))
            new2 = ifelse(is.na(df_tm[[x]]) == TRUE, 0, 1)
            new3 = hour(round(as.POSIXct(df_tm[[x]])))
        })
        new_name <- c(paste(x, "_", enddate, "_duration", sep = ""), paste(x, "IF", sep = "_"), paste(x, "hours", sep = "_"))
        setnames(mydata, c("new", "new2", "new3"), c(new_name))
        return(mydata[, new_name])
    }
}


#' time_varieble
#'
#' This function is not intended to be used by end user. 
#'
#' @param  dat  A data.frame.
#' @param  date_cols  Time variables.
#' @param  enddate  End time.
#' @export
time_varieble <- function(dat, date_cols = NULL, enddate = NULL) {
  
   dat = checking_data(dat = dat)
    date_cols1 = NULL
    if (!is.null(date_cols)) {
        date_cols1 <- names(dat)[colnames(dat) %islike% c(enddate,date_cols)]
    } else {
        date_cols1 = names(dat)
    }
    df_date = dat[date_cols1]
    df_date = time_transfer(dat = df_date, date_cols = c(enddate, date_cols))
    df_date <- df_date[!colAllnas(df_date)]
    df_tm = df_date[sapply(df_date, is_date)]

    time_vars_list <- lapply(date_cols1, function(x) time_vars_process(df_tm = df_tm,  x, enddate = enddate))
    index <- 0;
    j <- 1
    for (i in 1:length(time_vars_list)) {
        if (is.null(time_vars_list[[i]])) {
            index[j] <- i
            j <- j + 1
        }
    }
    tm_vars_tbl <- as.data.frame(Reduce("cbind", time_vars_list[-index]) )

    return(tm_vars_tbl)
}


#' Processing of Address Variables
#'
#' This function is not intended to be used by end user. 
#'
#' @param df_city A data.frame.
#' @param x Variables of city,
#' @param city_class  Class or levels of cities.
#' @export
city_varieble_process <- function(df_city, x, city_class) {
    if (class(df_city)[1] != "data.frame") {
        df_city <- as.data.frame(df_city)
    }
    df_city <- within(df_city, {
        city_level <- NA
        city_level[df_city[[x]] %alike% city_class[1]] <- 1
        city_level[df_city[[x]] %alike% city_class[2]] <- 2
        city_level[df_city[[x]] %alike% city_class[3]] <- 3
        city_level[df_city[[x]] %alike% city_class[4]] <- 4
        city_level[df_city[[x]] %alike% city_class[5]] <- 5
        city_level[df_city[[x]] %alike% city_class[6]] <- 6
        city_level[is.null(df_city[[x]]) == TRUE | df_city[[x]] == "NULL" | df_city[[x]] == "" |
            df_city[[x]] == "Unknown" | city_level == "NA" | df_city[[x]] == "NA"] <- -1
        city_level[is.na(city_level)] <- -1
    })
    NAsRate <- length(which(df_city$city_level == -1)) / nrow(df_city)
    if (NAsRate >= 0.3 & NAsRate < 0.6) {
        df_city2 <- data.frame()
        df_city2 <- within(df_city, {
            city_level <- NA
            city_level[df_city[[x]] %alike% city_class[1]] <- 1
            city_level[df_city[[x]] %alike% city_class[2]] <- 2
            city_level[df_city[[x]] %alike% city_class[3]] <- 3
            city_level[df_city[[x]] %alike% city_class[4]] <- 4
            city_level[df_city[[x]] %alike% city_class[5]] <- 4
            city_level[df_city[[x]] %alike% city_class[6]] <- 4
            city_level[is.null(df_city[[x]]) == TRUE | df_city[[x]] == "NULL" | df_city[[x]] == "" |
                df_city[[x]] == "Unknown" | city_level == "NA" | df_city[[x]] == "NA"] <- -1
            city_level[is.na(city_level)] <- -1
        })
    }
    if (NAsRate >= 0.6) {
        df_city3 <- data.frame()
        df_city3 <- within(df_city, {
            city_level <- NULL
            city_level[df_city[[x]] %alike% city_class[1]] <- 1
            city_level[df_city[[x]] %alike% city_class[2]] <- 1
            city_level[df_city[[x]] %alike% city_class[3]] <- 1
            city_level[df_city[[x]] %alike% city_class[4]] <- 1
            city_level[df_city[[x]] %alike% city_class[5]] <- 1
            city_level[df_city[[x]] %alike% city_class[6]] <- 1
            city_level[is.null(df_city[[x]]) == TRUE | df_city[[x]] == "NULL" | df_city[[x]] == "" |
                df_city[[x]] == "Unknown" | city_level == "NA" | df_city[[x]] == "NA"] <- -1
            city_level[is.na(city_level)] <- -1
        })
    }
    city_level_name <- paste(x, "city_level", sep = "_")
    df_city <- re_name(df_city, city_level ,city_level_name)
    return(df_city[city_level_name])
}

#' city_varieble
#'
#' This function is used for city variables derivation. 
#'
#' @param df  A data.frame.
#' @param city_cols Variables of city,
#' @param city_pattern  Regular expressions, used to match city variable names. Default is "city$".
#' @param city_class  Class or levels of cities.
#' @param parallel Logical, parallel computing. Default is TRUE.
#' @importFrom dplyr group_by mutate summarize  summarise n  count %>% filter left_join
#' @importFrom parallel detectCores  clusterExport clusterCall makeCluster stopCluster
#' @importFrom doParallel registerDoParallel  
#' @importFrom foreach foreach %dopar% %do%  registerDoSEQ
#' @export
city_varieble <- function(df = df, city_cols = NULL, city_pattern = "city$", city_class = city_class, parallel = TRUE) {
    if (class(df)[1] != "data.frame") {
        df <- as.data.frame(df)
    }
    if (is.null(city_cols)) {
        city_index <- grepl(city_pattern, paste(colnames(df)))
        city_cols <- names(df[city_index])
    } else {
        city_cols <- names(df[city_cols])
    }
    df_city = df[, city_cols]
    if (parallel) {
        parallel <- start_parallel_computing(parallel)
        stopCluster <- TRUE
    } else {
        parallel <- stopCluster <- FALSE
    }
    on.exit(if (parallel & stopCluster) stop_parallel_computing(attr(parallel, "cluster")))
    i. = NULL
    df_city_list = list()
    if (!parallel) {
        df_city_list <- lapply(city_cols, function(x) city_varieble_process(df_city, x, city_class))
        df_city_tbl <- Reduce("cbind", df_city_list) %>% as.data.frame()
    } else {
        df_city_list <- foreach(i. = city_cols, .combine = "c") %dopar% {
            try(do.call(city_varieble_process, args = list(df_city = df_city, x = i., city_class = city_class)), silent = TRUE)
        }
        df_city_tbl <- as.data.frame(df_city_list)
    }
    return(df_city_tbl)
}

#' add_variable_process
#'
#' This function is not intended to be used by end user. 
#'
#' @param  add  A data.frame contained address variables.
#' @export
add_variable_process <- function(add) {
    # acquire a sets of addresses 
    add1 = as.data.frame(add)
    sim1 = colname1 =list()
    for (i in 1:ncol(add1)) {
        if (i >= ncol(add1)) break
        sim1[[i]] <- apply(add1[, i:ncol(add1)], 2,
                       function(x) {
                           ifelse( add1[, i] %alike% x, 1, 0)})
        colname1[[i]] = lapply(names(add1)[i:(ncol(add1))], function(n) paste(names(add1)[i], n, sep = '_WITH_'))
    }
    sim1 = data.frame(t(unlist(sim1)), stringsAsFactors = FALSE)
    names(sim1) = unlist(colname1)
    # find the variables which are computing similarity with themselves
    splitvar <- strsplit(names(sim1), "_WITH_")
    vars <- c()
    for (i in 1:(length(sim1))) {
        if (splitvar[[i]][1] == splitvar[[i]][2]) {
            vars[[i]] <- names(sim1)[i]
        } else {
            vars[[i]] <- NA
        }
    }
    # get the final results
    sim = sim1[is.na(vars)]
    simm <- as.vector(sim)
    return(simm)
}


#' address_varieble
#'
#' This function is not intended to be used by end user. 
#'
#' @param df  A data.frame.
#' @param address_cols Variables of address,
#' @param address_pattern  Regular expressions, used to match address variable names.
#' @param parallel Logical, parallel computing. Default is TRUE.
#' @importFrom dplyr group_by mutate summarize  summarise n  count %>% filter left_join
#' @importFrom parallel detectCores  clusterExport clusterCall makeCluster stopCluster
#' @importFrom doParallel registerDoParallel  
#' @importFrom foreach foreach %dopar% %do%  registerDoSEQ
#' @export
address_varieble <- function(df, address_cols = NULL, address_pattern = NULL, parallel = TRUE) {
    if (class(df)[1] != "data.frame") {
        df <- as.data.frame(df)
    }
    if (is.null(address_cols)) {
        address_cols <- grepl(address_pattern, paste(colnames(df)))
        address_vars = names(df)[address_cols]
    } else {
        address_vars <- names(df[address_cols])
    }
    df_add = df[address_vars]
    if (parallel) {
        parallel <- start_parallel_computing(parallel)
        stopCluster <- TRUE
    } else {
        parallel <- stopCluster <- FALSE
    }
    on.exit(if (parallel & stopCluster) stop_parallel_computing(attr(parallel, "cluster")))
    i. = NULL
    df_add_list = list()
    if (!parallel) {
        df_add_list <- lapply(1:nrow(df_add), function(i.) add_variable_process(add = df_add[i.,]))
        df_add_tbl <- Reduce("cbind", df_add_list) %>% as.data.frame()
    } else {
        df_add_list <- foreach(i. = 1:nrow(df_add), .combine = "c") %dopar% {
            try(do.call(add_variable_process, args = list(add = df_add[i.,])), silent = TRUE)
        }
        df_add_tbl <- as.data.frame(df_add_list)
    }
    return(df_add_tbl)
}



#' variable_process
#'
#' This function is not intended to be used by end user. 
#'
#' @param  add  A data.frame
#' @importFrom data.table :=
#' @export
variable_process <- function(add) {
    td = new3= new2= grx_x= colname1=NULL

    # acquire a sets of addresses 
    cv_cols <- grep(grx_x, paste(colnames(dat)))[1:td]
    cv_cols <- cv_cols[!is.na(cv_cols)]
    #cv_folds
    if (length(cv_cols) > 0) {
        dat = dat[, new2 := ifelse(rowAllnas(dat[, cv_cols, with = FALSE]), NA,
        rowSums(dat[, cv_cols, with = FALSE], na.rm = TRUE))]
        dat = dat[, new3 := ifelse(rowAllnas(dat[, cv_cols, with = FALSE]), NA,
        rowMeans(dat[, cv_cols, with = FALSE], na.rm = TRUE))]
    }
    sim1 = data.frame(t(unlist(sim1)), stringsAsFactors = FALSE)
    names(sim1) = unlist(colname1)
    # find the variables which are computing similarity with themselves
    splitvar <- strsplit(names(sim1), "_WITH_")
    vars <- c()
    for (i in 1:(length(sim1))) {
        if (splitvar[[i]][1] == splitvar[[i]][2]) {
            vars[[i]] <- names(sim1)[i]
        } else {
            vars[[i]] <- NA
        }
    }
    # get the final results
    sim = sim1[is.na(vars)]
    simm <- as.vector(sim)
    return(simm)
}
