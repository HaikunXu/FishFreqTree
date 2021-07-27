#' Run the regression tree analysis with specified number of splits
#'
#' \code{find_split} This function finds the best split based on the improvement
#'
#' @param lf_prop The length frequency data frame input; must include columns lat, lon, year, and quarter
#' @param frstcol.lf The first column in the data frame with length frequency info
#' @param lstcol.lf The first column in the data frame with length frequency info
#'
#' @export

find_split <- function(lf_prop,frstcol.lf,lstcol.lf) {
  # run the code to get improvement by every split
  tmpcl.kld <- simult.tree.kld.FINAL(lf_prop,frstcol.lf,lstcol.lf)

  # combine the result tables into a single data frame
  tmpcl.kld.lat <- data.frame("Improve"=tmpcl.kld$lf.lat[,1],
                              "Key"="Lat",
                              "Value"=tmpcl.kld$lf.lat[,2])
  tmpcl.kld.lon <- data.frame("Improve"=tmpcl.kld$lf.lon[,1],
                              "Key"="Lon",
                              "Value"=tmpcl.kld$lf.lon[,2])
  tmpcl.kld.qrt <- data.frame("Improve"=tmpcl.kld$lf.qrtr[,1],
                              "Key"="Qrt",
                              "Value"=tmpcl.kld$lf.qrtr[,2])
  if(is.null(tmpcl.kld$lf.cyclic.qrtr)) {
    tmpcl.kld.df <- rbind(tmpcl.kld.lat,tmpcl.kld.lon,tmpcl.kld.qrt)
  }
  else {
    tmpcl.kld.cqrt <- data.frame("Improve"=tmpcl.kld$lf.cyclic.qrtr[,1],
                                 "Key"="CQrt",
                                 "Value"=tmpcl.kld$lf.cyclic.qrtr[,2])
    tmpcl.kld.df <- rbind(tmpcl.kld.lat,tmpcl.kld.lon,
                          tmpcl.kld.qrt,tmpcl.kld.cqrt)
  }

  # sort the improvement and the 1st one is the best split
  tmpcl.kld.df <- arrange(tmpcl.kld.df,desc(Improve))

  # print(tmpcl.kld.df[1:3,])

  return(tmpcl.kld.df)
  }
