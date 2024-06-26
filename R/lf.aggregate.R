#' Aggregate length counts data
#'
#' \code{lf.aggregate} aggregates counts data by year, quarter, lat, lon, and user-sepcified length bins
#'
#' @param LC The length frequency data frame input; must include four columns: lat, lon, year, and quarter
#' @param fcol The first column in the data frame with length frequency info
#' @param lcol The last column in the data frame with length frequency info
#' @param bins Names of all bins included in LF
#' @param new_bins Names of all bins included in LF
#' @param LengthOnly Whether aggregate data by new length bins only
#' @param minN minimal sample size required for a year-quarter-lat-lon stratum
#'
#' @export

lf.aggregate <- function(LC,fcol,lcol,bins,new_bins,LengthOnly=FALSE,minN=0) {
  # check data first
  if((lcol-fcol+1)!= length(bins)) stop("Error! The number of bins does not match the number of LF columns specified")
  else names(LC)[fcol:lcol] <- bins

  if(is.null(LC[["weight"]])==TRUE) LC$weight <- 1
  else print("Weight is used to aggregate observations!")

  if(LengthOnly==FALSE) { # aggregate by lat, lon, year, quarter, and length

    LC_select <- LC[,c("year","quarter","lat","lon",paste0(bins),"weight")]
    LC_long <- data.frame(tidyr::gather(LC_select,fcol:lcol,key = "length",value = "number"))

    # user-specified bins
    LC_long$Length <- cut(as.numeric(LC_long$length), breaks = c(new_bins, Inf), right=F, labels = new_bins)

    # remove Length = NA
    LC_long <- dplyr::filter(LC_long,is.na(Length)==FALSE)

    # total counts per year, quarter, lat, lon, and length (L)
    LC_total <- dplyr::summarise(dplyr::group_by(LC_long, year, lat, lon, quarter, Length), total_n=sum(number*weight))

    # remove the year-quarter-lat-lon strata with less than minN samples
    LC_total <- dplyr::mutate(dplyr::group_by(LC_total, year, lat, lon, quarter), total_N=sum(total_n))
    LC_total <- dplyr::filter(LC_total,total_N>minN)

    # aggregate counts to lf
    LF <- dplyr::mutate(dplyr::group_by(LC_total, year, lat, lon, quarter), lf=total_n/total_N)

    # spread the lf
    LF_aggregate <- tidyr::spread(dplyr::select(LF,year,quarter,lat,lon,Length,lf),Length,lf,fill = 0)
  }
  else { # aggregate by length

    LC_select <- LC[,c("year","quarter","lat","lon",paste0(bins),"weight")]
    LC_select$ID <- seq(1,nrow(LC_select)) # each row has a unique ID
    LC_long <- data.frame(tidyr::gather(LC_select,fcol:lcol,key = "length",value = "number"))

    # user-specified bins
    LC_long$Length <- cut(as.numeric(LC_long$length), breaks = c(new_bins, Inf), right=F, labels = new_bins)

    # remove Length = NA
    LC_long <- dplyr::filter(LC_long,is.na(Length)==FALSE)

    # total counts per length (L)
    LC_total <- dplyr::summarise(dplyr::group_by(LC_long, ID, year, quarter, lat, lon, Length), total_n=sum(number))

    # remove the year-quarter-lat-lon strata with less than minN samples
    LC_total <- dplyr::mutate(dplyr::group_by(LC_total, year, lat, lon, quarter), total_N=sum(total_n))
    LC_total <- dplyr::filter(LC_total,total_N>=minN)

    # aggregate counts to lf
    LF <- dplyr::mutate(dplyr::group_by(LC_total, ID), lf=total_n/sum(total_n))

    # spread the lf
    LF_aggregate <- tidyr::spread(dplyr::select(LF,ID,year,quarter,lat,lon,Length,lf),Length,lf,fill = 0)

  }
  return(LF_aggregate)
}
