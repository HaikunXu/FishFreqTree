#' Aggregate length counts data
#'
#' \code{lf.demean} aggregates counts data by year, quarter, lat, lon, and user-sepcified length bins
#'
#' @param LF The length frequency data frame input; must include four columns: lat, lon, year, and quarter
#' @param fcol The first column in the data frame with length frequency info
#' @param lcol The last column in the data frame with length frequency info
#' @param bins Names of all bins included in LF
#'
#' @export

lf.demean <- function(LF,fcol,lcol,bins) {
  # check data first
  if((lcol-fcol+1)!= length(bins)) stop("Error! The number of bins does not match the number of LF columns specified")
  else names(LF)[fcol:lcol] <- bins

  # if(is.null(LF[["weight"]])==TRUE) LF$weight <- 1
  # else print("Weight is used to compute the mean length!")

  LF_select <- LF[,c("ID","year","quarter","lat","lon",paste0(bins))]
  LF_long <- data.frame(tidyr::gather(LF_select,fcol:lcol,key = "Length",value = "LF"))

  LF_mean <- dplyr::mutate(dplyr::group_by(LF_long, year, quarter, Length), LF_mean=sum(LF))
  LF_mean$LF_demean <- LF_mean$LF / LF_mean$LF_mean

  LF_demean <- tidyr::spread(dplyr::select(LF_mean,ID,year,quarter,lat,lon,Length,LF_demean),Length,LF_demean,fill = 0)

  return(LF_demean)
}
