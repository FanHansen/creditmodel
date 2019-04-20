#' split_bins
#' 
#' \code{split_bins} is  for binning using breaks.
#' @param dat A data.frame with independent variables.
#' @param x  The name of an independent variable.
#' @param breaks  Breaks for binning.
#' @param bins_no Number the generated bins. Default is TRUE.
#' @return  A data.frame with Bined x.
#' @examples
#' \dontrun{
#' breaks = get_breaks(dat = UCICreditCard, x = "PAY_AMT1", 
#' target = "default.payment.next.month", occur_time = "apply_date", 
#' sp_values = list(-1, "Unknown"), 
#' tree_control = NULL, bins_control = NULL)
#' bins = split_bins(dat = UCICreditCard, 
#' x = "PAY_AMT1", breaks = breaks, bins_no = TRUE)
#' }
#' @export


split_bins <- function(dat, x, breaks = NULL,  bins_no = TRUE) {
    opt = options(stringsAsFactors = FALSE) # 
    if (length(breaks) < 1) {
        breaks = get_breaks(dat, x, target = NULL, best = FALSE, equal_bins = TRUE, g = 5,note = FALSE)
    }
    miss_value_num = miss_value_char = NULL
    if (any(c("integer", "numeric", "double") == class(dat[, x]))) {
        breaks = sort(unlist(unique(c(-Inf, breaks, Inf))))
        bins_1 = cut(dat[, x], breaks = unique(breaks), dig.lab = 10, ordered = TRUE, include.lowest = FALSE, right = FALSE)
        if (bins_no) {
            bins_0 = paste("0", as.numeric(bins_1), sep = "")
            bins = paste(bins_0, bins_1, sep = ".")
            bins[which(as.numeric(bins_1) >= 10)] = paste(as.numeric(bins_1[which(as.numeric(bins_1) >= 10)]), bins_1[which(as.numeric(bins_1) >= 10)], sep = ".")
        } else {
            bins = as.character(bins_1)
        }
    } else {
        breaks = unique(breaks)
        if (any(grepl("\\|", breaks))) {
            breaks_s = strsplit(breaks, "\\|")
        } else {
            breaks_s = breaks
        }
        dat[which(!(dat[, x] %in% unlist(breaks_s))), x] = get_median(dat[, x])
        if (length(breaks_s) > 0) {
            for (i in 1:length(breaks_s)) {
                if (length(which(dat[, x] %in% unlist(breaks_s[[i]]))) > 1) {
                    if (i < 10) {
                        dat[which(dat[, x] %in% unlist(breaks_s[[i]])), x] = paste(paste0("0", i), paste(breaks_s[[i]], collapse = ";"), sep = ".")
                    } else {
                        dat[which(dat[, x] %in% unlist(breaks_s[[i]])), x] = paste(paste0( i), paste(breaks_s[[i]], collapse = ";"), sep = ".")
                    }            
                }
            }
        }
        bins = dat[, x]
    }
    options(opt) # reset
    return(bins)
}
