#' Run the regression tree analysis with specified number of splits
#'
#' \code{simult.tree.kld.FINAL} This function finds the best split based on the improvement
#'
#' @param lfinput.frm The length frequency data frame input; must include columns lat, lon, year, and quarter
#' @param frstcol.lf The first column in the data frame with length frequency info
#' @param lstcol.lf The first column in the data frame with length frequency info
#'
#' @export

simult.tree.kld.FINAL <- function(lfinput.frm,frstcol.lf,lstcol.lf,lat.min=lat.min,lon.min=lon.min)
{
  # 10.9.2012: added line to skip quarter if only one unique quarter value
  # calling function for kld part of simultaneous tree
  # note: 0 (1) is left (right) or <=x (>x)
  # lfinput.frm is a data frame with columns for latitude, longitude,
  #     quarter and the length bins (rows are samples, columns predictor
  #     variables and length bins).
  # it is assumed that length bin columns are consecutive from frstcol.lf
  #     to lstcol.lf.
  #

  # LATITUDE
  unq.lats<-sort(unique(lfinput.frm$lat))
  nunqlats<-length(unq.lats)
  #
  if(nunqlats>lat.min) {
    lfimp.lat<-rep(0,(nunqlats-1))

  for(i in lat.min:(nunqlats-lat.min)){
    # make left-right split flag
    lftrght.splitflg<-rep(0,length(lfinput.frm$lat))
    lftrght.splitflg[lfinput.frm$lat>unq.lats[i]]<-1
    # compute kld contribution for this split
    lfimp.lat[i]<-imp.kld.FINAL.R(as.matrix(lfinput.frm[,frstcol.lf:lstcol.lf]),lftrght.splitflg)
  }
  }
  else {
    nunqlats <- 2
    lfimp.lat <- -1
  }
  #
  # LONGITUDE
  unq.lons<-sort(unique(lfinput.frm$lon))
  nunqlons<-length(unq.lons)
  #
  if(nunqlons>lon.min) {
    lfimp.lon<-rep(0,(nunqlons-1))

  for(i in lon.min:(nunqlons-lon.min)){
    # make left-right split flag
    lftrght.splitflg<-rep(0,length(lfinput.frm$lon))
    lftrght.splitflg[lfinput.frm$lon>unq.lons[i]]<-1
    # compute kld contribution for this split
    lfimp.lon[i]<-imp.kld.FINAL.R(as.matrix(lfinput.frm[,frstcol.lf:lstcol.lf]),lftrght.splitflg)
  }
  }
  else {
    nunqlons <- 2
    lfimp.lon <- -1
  }
  #
  # QUARTER
  # first the splits that arise from treating quarter as a continuous variable
  unq.qrtrs<-sort(unique(lfinput.frm$quarter))
  nunqqrtrs<-length(unq.qrtrs)
  # 10.9.12: added if() statement below
  if(nunqqrtrs>1) {
    lfimp.qrtr<-rep(NA,(nunqqrtrs-1))
    #
    for(i in 1:(nunqqrtrs-1)){
      # make left-right split flag
      lftrght.splitflg<-rep(0,length(lfinput.frm$quarter))
      lftrght.splitflg[lfinput.frm$quarter>unq.qrtrs[i]]<-1
      # compute kld contribution for this split
      lfimp.qrtr[i]<-imp.kld.FINAL.R(as.matrix(lfinput.frm[,frstcol.lf:lstcol.lf]),lftrght.splitflg)
    }
    # now the cyclic quarter splits (only if there are >= 3 quarters in input data)
    # for 3 quarters this will be one of: 2;1,3 or 3;2,4 or 3;1,4 or 2;1,4
    # for 4 quarters there a three cases: 1,4;2,3 and 1,2,4;3 and 1,3,4;2
    if(nunqqrtrs==3) {
      lftrght.splitflg<-rep(0,length(lfinput.frm$quarter))
      lftrght.splitflg[lfinput.frm$quarter==unq.qrtrs[1] | lfinput.frm$quarter==unq.qrtrs[3]]<-1
      acyclic.qrtr<-paste(unq.qrtrs[2],paste(unq.qrtrs[1],unq.qrtrs[3],sep=","),sep=";")
      # compute kld contribution for this split
      lfimp.cyclic.qrtr<-imp.kld.FINAL.R(as.matrix(lfinput.frm[,frstcol.lf:lstcol.lf]),lftrght.splitflg)
    }
    #
    if(nunqqrtrs==4) {
      lfimp.cyclic.qrtr<-rep(NA,3)
      acyclic.qrtr<-rep(NA,3)
      #
      lftrght.splitflg<-rep(0,length(lfinput.frm$quarter))
      lftrght.splitflg[lfinput.frm$quarter==2 | lfinput.frm$quarter==3]<-1
      acyclic.qrtr[1]<-paste(paste("1","4",sep=""),paste("2","3",sep=""),sep=";")
      # compute kld contribution for this split
      lfimp.cyclic.qrtr[1]<-imp.kld.FINAL.R(as.matrix(lfinput.frm[,frstcol.lf:lstcol.lf]),lftrght.splitflg)
      #
      lftrght.splitflg<-rep(0,length(lfinput.frm$quarter))
      lftrght.splitflg[lfinput.frm$quarter==3]<-1
      acyclic.qrtr[2]<-paste("1","2","4",";","3",sep="")
      # compute kld contribution for this split
      lfimp.cyclic.qrtr[2]<-imp.kld.FINAL.R(as.matrix(lfinput.frm[,frstcol.lf:lstcol.lf]),lftrght.splitflg)
      #
      lftrght.splitflg<-rep(0,length(lfinput.frm$quarter))
      lftrght.splitflg[lfinput.frm$quarter==2]<-1
      acyclic.qrtr[3]<-paste("1","3","4",";","2",sep="")
      # compute kld contribution for this split
      lfimp.cyclic.qrtr[3]<-imp.kld.FINAL.R(as.matrix(lfinput.frm[,frstcol.lf:lstcol.lf]),lftrght.splitflg)
    }
  }
  #
  # make list of results (modified 10.9.2012)
  # no quarter splits
  if(nunqqrtrs==1) {
    output.frm<-list(lf.lat=data.frame(lfimp.lat,unq.lats[1:(nunqlats-1)]),lf.lon=data.frame(lfimp.lon,unq.lons[1:(nunqlons-1)]))
  }
  # quarter splits but no cyclic quarters
  if(nunqqrtrs==2){
    output.frm<-list(lf.lat=data.frame(lfimp.lat,unq.lats[1:(nunqlats-1)]),lf.lon=data.frame(lfimp.lon,unq.lons[1:(nunqlons-1)]),lf.qrtr=data.frame(lfimp.qrtr,unq.qrtrs[1:(nunqqrtrs-1)]))
  }
  # have quarter splits with cyclic quarters
  if(nunqqrtrs>2){
    output.frm<-list(lf.lat=data.frame(lfimp.lat,unq.lats[1:(nunqlats-1)]),lf.lon=data.frame(lfimp.lon,unq.lons[1:(nunqlons-1)]),lf.qrtr=data.frame(lfimp.qrtr,unq.qrtrs[1:(nunqqrtrs-1)]), lf.cyclic.qrtr=data.frame(lfimp.cyclic.qrtr,acyclic.qrtr))
  }
  #
  return(output.frm)
}
