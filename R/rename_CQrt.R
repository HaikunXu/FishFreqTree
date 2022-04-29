#' Run the regression tree analysis with specified number of splits
#'
#' \code{rename_CQrt} This function finds the best split based on the improvement
#'
#' @export

rename_CQrt <- function(split) {

  split$Value[which((split$Key=="CQrt")&(split$Value==1))] <- "14;23"
  split$Value[which((split$Key=="CQrt")&(split$Value==2))] <- "124;3"
  split$Value[which((split$Key=="CQrt")&(split$Value==3))] <- "134;2"

  # split$Value[which((split$Key=="Qrt")&(split$Value==1))] <- "1;234"
  # split$Value[which((split$Key=="Qrt")&(split$Value==2))] <- "12;34"
  # split$Value[which((split$Key=="Qrt")&(split$Value==3))] <- "123;4"

  return(split)
}
