#' Find the best split
#'
#' \code{find_split} This function finds the best split based on the improvement
#'
#' @param lf_prop The length frequency data frame input; must include columns lat, lon, year, and quarter
#' @param frstcol.lf The first column in the data frame with length frequency info
#' @param lstcol.lf The first column in the data frame with length frequency info
#' @param quarter Whether to consider quarter as a splitting dimention; default = TRUE
#'
#' @export

find_split <- function(lf_prop,frstcol.lf,lstcol.lf,lat.min,lon.min,year.min,quarter,year) {
  # run the code to get improvement by every split
  tmpcl.kld <- simult.tree.kld.FINAL(lf_prop,frstcol.lf,lstcol.lf,lat.min=lat.min,lon.min=lon.min,year.min = year.min)

  # combine the result tables into a single data frame
  tmpcl.kld.lat <- data.frame("Improvement"=tmpcl.kld$lf.lat[,1],
                              "Key"="Lat",
                              "Value"=tmpcl.kld$lf.lat[,2])
  tmpcl.kld.lon <- data.frame("Improvement"=tmpcl.kld$lf.lon[,1],
                              "Key"="Lon",
                              "Value"=tmpcl.kld$lf.lon[,2])

  if(quarter==TRUE) {

  if((is.null(tmpcl.kld$lf.cyclic.qrtr)==TRUE)&(is.null(tmpcl.kld$lf.qrtr)==FALSE)) {
    tmpcl.kld.qrt <- data.frame("Improvement"=tmpcl.kld$lf.qrtr[,1],
                                "Key"="Qrt",
                                "Value"=tmpcl.kld$lf.qrtr[,2])
    tmpcl.kld.df <- rbind(tmpcl.kld.lat,tmpcl.kld.lon,tmpcl.kld.qrt)
  }
  if((is.null(tmpcl.kld$lf.cyclic.qrtr)==FALSE)&(is.null(tmpcl.kld$lf.qrtr)==TRUE)) {
    tmpcl.kld.cqrt <- data.frame("Improvement"=tmpcl.kld$lf.cyclic.qrtr[,1],
                                 "Key"="CQrt",
                                 "Value"=as.numeric(tmpcl.kld$lf.cyclic.qrtr[,2]))
    tmpcl.kld.df <- rbind(tmpcl.kld.lat,tmpcl.kld.lon,tmpcl.cyclic.kld.qrt)
  }
  if((is.null(tmpcl.kld$lf.cyclic.qrtr)==TRUE)&(is.null(tmpcl.kld$lf.qrtr)==TRUE)) {
    tmpcl.kld.df <- rbind(tmpcl.kld.lat,tmpcl.kld.lon)
  }
  if((is.null(tmpcl.kld$lf.cyclic.qrtr)==FALSE)&(is.null(tmpcl.kld$lf.qrtr)==FALSE)) {
    tmpcl.kld.cqrt <- data.frame("Improvement"=tmpcl.kld$lf.cyclic.qrtr[,1],
                                 "Key"="CQrt",
                                 "Value"=as.numeric(tmpcl.kld$lf.cyclic.qrtr[,2]))
    tmpcl.kld.qrt <- data.frame("Improvement"=tmpcl.kld$lf.qrtr[,1],
                                "Key"="Qrt",
                                "Value"=tmpcl.kld$lf.qrtr[,2])
    tmpcl.kld.df <- rbind(tmpcl.kld.lat,tmpcl.kld.lon,tmpcl.kld.qrt,tmpcl.kld.cqrt)
  }
  }
  else {tmpcl.kld.df <- rbind(tmpcl.kld.lat,tmpcl.kld.lon)}

  if (year==TRUE) {
    tmpcl.kld.year <- data.frame("Improvement"=tmpcl.kld$lf.year[,1],
                                "Key"="Year",
                                "Value"=tmpcl.kld$lf.year[,2])
    tmpcl.kld.df <- rbind(tmpcl.kld.df,tmpcl.kld.year)
  }
  # sort the improvement and the 1st one is the best split

  tmpcl.kld.df <- tmpcl.kld.df[order(tmpcl.kld.df$Improve,decreasing = TRUE),]

  return(tmpcl.kld.df)
  }
