#' UCI Credit Card data
#' 
#'  This research aimed at the case of customers's default payments in Taiwan and compares the predictive accuracy of probability of default among six data mining methods.
#'  This research employed a binary variable, default payment (Yes = 1, No = 0), as the response variable. This study reviewed the literature and used the following 24 variables as explanatory variables
#' 
#' \itemize{
#'   \item ID: Customer id
#'   \item apply_date: This is a fake occur time. 
#'   \item LIMIT_BAL: Amount of the given credit (NT dollar): it includes both the individual consumer credit and his/her family (supplementary) credit. 
#'   \item SEX:  Gender (male; female). 
#'   \item EDUCATION: Education (1 = graduate school; 2 = university; 3 = high school; 4 = others).
#'   \item MARRIAGE: Marital status (1 = married; 2 = single; 3 = others). 
#'   \item AGE: Age (year)
#' History of past payment. We tracked the past monthly payment records (from April to September, 2005) as follows:
#'   \item PAY_0: the repayment status in September
#'   \item PAY_2: the repayment status in August
#'   \item PAY_3: ...
#'   \item PAY_4: ...
#'   \item PAY_5: ...
#'   \item PAY_6: the repayment status in April
#' The measurement scale for the repayment status is: -1 = pay duly; 1 = payment delay for one month; 2 = payment delay for two months;...;8 = payment delay for eight months; 9 = payment delay for nine months and above. 
#' Amount of bill statement (NT dollar)
#'   \item BILL_AMT1: amount of bill statement in September
#'   \item BILL_AMT2: mount of bill statement in August
#'   \item BILL_AMT3: ...
#'   \item BILL_AMT4: ...
#'   \item BILL_AMT5: ...
#'   \item BILL_AMT6: amount of bill statement in April
#'  Amount of previous payment (NT dollar)
#'   \item PAY_AMT1: amount paid in September
#'   \item PAY_AMT2: amount paid in August
#'   \item PAY_AMT3: ....
#'   \item PAY_AMT4: ...
#'   \item PAY_AMT5: ...
#'   \item PAY_AMT6: amount paid in April
#'   \item default.payment.next.month: default payment (Yes = 1, No = 0), as the response variable
#' }
#' 
#' @docType data
#' @keywords datasets
#' @format A data frame with 30000 rows and 26 variables.
#' @name UCICreditCard
#' @source \url{http://archive.ics.uci.edu/ml/datasets/default+of+credit+card+clients}
#' @seealso \code{\link{lendingclub}}
NULL