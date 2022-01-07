#' Run the regression tree analysis with specified number of splits
#'
#' \code{get.klderror.null} This function generates an area code for longline catch allocation
#'
#' @export

get.klderror.null <- function (response.mat, weight)
{
  sumerror <- 0
  nlngth.ints <- ncol(response.mat)
  prop.bar <- apply(response.mat, 2, mean)
  for (j in 1:nlngth.ints) {
    props.tmp <- response.mat[response.mat[, j] > 0, j]
    weight.tmp <- weight[response.mat[, j] > 0]
    # weight.tmp <- 1 # change in the future
    sumerror <- sumerror + sum(props.tmp * log(props.tmp/prop.bar[j]) * weight.tmp)
  }
  return(sumerror)
}
