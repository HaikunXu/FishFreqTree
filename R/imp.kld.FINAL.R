#' Run the regression tree analysis with specified number of splits
#'
#' \code{imp.kld.FINAL.R} This function generates an area code for longline catch allocation
#'
#' @export


imp.kld.FINAL.R <- function(lf.y,lf.splitflg)
{
  # computes l-f contribution to combined improvement for simultaneous tree.
  # lf.y is matrix of l-f proportions (rows are samples, columns are length-bins).
  # lf.splitflg gives 0 (1) for left (right) at a given split value of a specific predictor.
  # lf.splitflg is a vector with 0 or 1 for each sample in lf.y matrix.
  # this function is called by simult.tree.kld.FINAL.f
  # improvement is based on the Kullback-Leilber impurity of Lennert-Cody et al 2010.
  # improvement is computed from entropy.
  # the improvement for this split, improve.kld, is eq. (2) of Lennert-Cody et al. 2010.
  # reference: Fisheries Research 102(2010):323-326.
  #
  # first get overall, left and right mean proportions for each length-bin
  #  (means are computed over samples in parent and daughter groups)
  lf.pbar.all<-apply(lf.y,2,mean)
  if(!is.vector(lf.y[lf.splitflg==0,])){
    lf.pbar.lft<-apply(lf.y[lf.splitflg==0,],2,mean)
  } else {
    lf.pbar.lft<-lf.y[lf.splitflg==0,]
  }
  if(!is.vector(lf.y[lf.splitflg==1,])) {
    lf.pbar.rght<-apply(lf.y[lf.splitflg==1,],2,mean)
  } else {
    lf.pbar.rght<-lf.y[lf.splitflg==1,]
  }
  #
  # improvement for this split
  # m is number of length bins
  # the various n are number of samples in parent and daughter groups
  m<-ncol(lf.y)
  nlft<-length(lf.splitflg[lf.splitflg==0])
  nrght<-length(lf.splitflg[lf.splitflg==1])
  nall<-length(lf.splitflg)
  # compute entropy for parent and daughter groups
  plogp.all<-0
  plogp.lft<-0
  plogp.rght<-0
  # looping over length bins
  for(j in 1:m){
    if(lf.pbar.all[j]>0) {
      plogp.all<-plogp.all+(lf.pbar.all[j]*log(lf.pbar.all[j]))
    }
    #
    if(lf.pbar.lft[j]>0) {
      plogp.lft<-plogp.lft+(lf.pbar.lft[j]*log(lf.pbar.lft[j]))
    }
    #
    if(lf.pbar.rght[j]>0) {
      plogp.rght<-plogp.rght+(lf.pbar.rght[j]*log(lf.pbar.rght[j]))
    }
  }
  #
  improve.kld<-(-1*nall*plogp.all)-(-1*nlft*plogp.lft-nrght*plogp.rght)
  #
  return(improve.kld)
}
