#' Run the regression tree analysis with specified number of splits
#'
#' \code{get.klderror.null} This function generates an area code for longline catch allocation
#'
#' @export

get.klderror.null <- function (response.mat)
{
  sumerror <- 0
  nlngth.ints <- ncol(response.mat)
  prop.bar <- apply(response.mat, 2, mean)
  for (j in 1:nlngth.ints) {
    props.tmp <- response.mat[response.mat[, j] > 0, j]
    sumerror <- sumerror + sum(props.tmp * log(props.tmp/prop.bar[j]))
  }
  return(sumerror)
}
