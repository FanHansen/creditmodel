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

PCA_reduce <- function(train = train, test = NULL, mc = 0.9) {
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
        names(PCA_data)[i] = paste(Reduce("mypaste", Reduce("union", names_strsplit)), names(PCA_data)[i], sep = "_")
    }
    return(PCA_data)
}
