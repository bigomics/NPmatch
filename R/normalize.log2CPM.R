#' Normalize Data
#'
#' @param X raw data matrix (i.e., gene expression dataset): rows are features/genes; columns are samples.
#' @param total total count within each sample
#' @param prior offset integer added prior log transformation
#'
#' @return CPM normalized data matrix
#' @export
#' @import stats
#' @import utils
#'
#' @examples# load the example dataset provided and run:
#' nX <- normalize.log2CPM(X)
#'
normalize.log2CPM <- function(X,
                              total = 1e6,
                              prior = 1) {

    ## If total counts per sample is < 1e6 uses average count, else use 1e6.
    total0 <- mean(Matrix::colSums(X, na.rm = TRUE))
    total <- ifelse(total0 < 1e6, total0, 1e6)

    message("[log2CPM] setting column sums to = ", round(total, 2))
    totcounts <- Matrix::colSums(X, na.rm = TRUE)
    cpm <- sweep(X, 2, totcounts, FUN = "/") * total
    nX <- log2(prior + cpm)
    return(nX)

}
