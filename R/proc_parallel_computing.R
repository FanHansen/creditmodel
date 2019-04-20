#' Parallel computing and export variables to global Env.
#'
#' @param parallel  A logical, default is TRUE.
#'
#' @return  parallel works.
#'
#' @examples
#' \dontrun{
#' parallel = start_parallel_computing(parallel)
#' }
#' @importFrom parallel detectCores  clusterExport clusterCall makeCluster stopCluster
#' @importFrom doParallel registerDoParallel  
#' @importFrom foreach foreach %dopar% %do%  registerDoSEQ
#' @export
start_parallel_computing <- function(parallel = TRUE) {

    parallelType <- if (.Platform$OS.type == "windows")
        "snow" else "multicore"

    numCores <- parallel::detectCores()
    attr(parallel, "type") <- parallelType
    attr(parallel, "cores") <- numCores

    if (parallel) {
        if (parallelType == "snow") {

            cl <- parallel::makeCluster(numCores, type = "PSOCK")
            attr(parallel, "cluster") <- cl
  
            varlist <- ls(envir = parent.frame(), all.names = TRUE)
            varlist <- varlist[varlist != "..."]
            parallel::clusterExport(cl, varlist = varlist, envir = parent.frame())
            parallel::clusterExport(cl, varlist = ls(envir = globalenv(), all.names = TRUE), envir = globalenv())

            pkgs <- .packages()
            lapply(pkgs, function(pkg)
                parallel::clusterCall(cl, library, package = pkg,
                                     character.only = TRUE))
            doParallel::registerDoParallel(cl, cores = numCores)
        }
        else if (parallelType == "multicore") {
            cl <- parallel::makeCluster(numCores, type = "FORK")

            doParallel::registerDoParallel(cl, cores = numCores[1])

            attr(parallel, "cluster") <- cl
        }
        else { stop("Only 'snow' and 'multicore' clusters allowed!") }
        }
    return(parallel)
}

#' Stop parallel computing
#'
#' @param cluster  Parallel works.
#'
#' @return  stop clusters.
#'
#' @examples
#' \dontrun{
#' parallel = start_parallel_computing(parallel)
#' stop_parallel_computing(attr(parallel, "cluster"))
#' }
#' @importFrom parallel  stopCluster
#' @importFrom foreach  registerDoSEQ
#' @export
stop_parallel_computing <- function(cluster) {
    parallel::stopCluster(cluster)
    foreach::registerDoSEQ()
    invisible()
}



#'  Loop Function. 
#' #' \code{loop_function} is an iterator to loop through 
#' @param func  A function.
#' @param args  A list of argauments required by function.
#' @param x_list  Names of objects to loop through.
#' @param bind  Complie results, "rbind" & "cbind" are available.
#' @param parallel  Logical, parallel computing.
#' @param as_list  Logical, whether outputs  to be a list.
#'
#' @return   A data.frame or list
#'
#' @examples
#' \dontrun{
#' parallel = start_parallel_computing(parallel)
#' stop_parallel_computing(attr(parallel, "cluster")))
#' }
#' @importFrom foreach foreach %dopar% %do% 
#' @export

loop_function <- function(func = NULL, args = list(data = NULL), x_list = NULL, bind = "rbind", parallel = TRUE, as_list = FALSE) {
    opt = options("warn" = -1) # suppress warnings
    df_list = df_tbl = NULL
    if (parallel) {
        parallel <- start_parallel_computing(parallel)
        stopCluster <- TRUE
    } else {
        parallel <- stopCluster <- FALSE
    }
    on.exit(if (parallel & stopCluster) stop_parallel_computing(attr(parallel, "cluster")))
 
    if (is.null(func)) {
        stop("The function is missing. Please input a split_type function.\n")
    }
    if (is.null(x_list)) {
        x_list = get_names(data, types = c("character","factor","numeric","integer","double"))
        warning(paste("x_list is NULL, use all variables.\n" ))
    }
    i. = NULL
    if (!parallel) {
        funct = function(i.) {
            try(do.call(func, c(args,x = i.)), silent = FALSE)
        }
        df_list <- lapply(x_list, funct)
        if (as_list) {
            df_tbl <- df_list
            names(df_tbl) <- x_list
        } else {
            df_tbl <- as.data.frame(Reduce(bind, df_list))
        }
    } else {
        df_list <- foreach(i. = x_list, .errorhandling = c('pass')) %dopar% {
            try(do.call(func, args = c(args,x = i.)), silent = FALSE)
        }
        if (as_list) {
            df_tbl <- df_list
            names(df_tbl) <- x_list
        } else {
            df_tbl <- as.data.frame(Reduce(bind, df_list))
        }
    }
    options(opt) # reset warnings
    return(df_tbl)
}

