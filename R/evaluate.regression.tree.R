#' Evaluate a pre-specified regression tree
#'
#' \code{evaluate.regression.tree} This function is the main regression tree function
#'
#' @param LF The length frequency data frame input; must include four columns: lat, lon, year, and quarter
#' @param fcol The first column in the data frame with length frequency info
#' @param lcol The last column in the data frame with length frequency info
#' @param Flagcol The column for pre-defined fishery definitions
#' @param bins Names of all bins included in LF
#' @param save_dir The directory where results will be saved
#' @param folder_name The name of the folder where results will be saved
#' @param pdf Whether to save figures in pdf - default is in png
#' @param lf_scale scale mean lf plot
#'
#' @export

evaluate_regression_tree <- function(LF,fcol,lcol,Flagcol,bins,save_dir,folder_name,pdf=FALSE,lf_scale=1.5) {

  # make sure bins are consistent with LF input
  if((lcol-fcol+1)!= length(bins)) stop("Error! The number of bins does not match the number of LF columns specified")
  else names(LF)[fcol:lcol] <- bins

  LF$dummy <- FALSE

  # make sure LF sums to 1 for each row
  row_sum <- apply(LF[which(LF$dummy==FALSE),fcol:lcol],1,sum)
  if(sum(abs(row_sum-1)>0.05)>0) {
    plot(row_sum)
    # print(which(abs(row_sum-1)>0.1))
    stop("Error! LF does not sum to 1 for at least one row.")
  }

  if(is.null(LF[["weight"]])==TRUE) LF$weight <- 1
  else print("Weight is used to compute the percentage of varaibility explained!")

  unlink(paste0(save_dir,folder_name), recursive = TRUE)
  dir.create(paste0(save_dir,folder_name))

  e0 <- get.klderror.null(as.matrix(LF[LF$dummy==FALSE,fcol:lcol]),LF$weight[LF$dummy==FALSE]) # null (no stratification)

  e <- 0
  for (i in unique(LF[[Flagcol]])) {
    e <- e + get.klderror.null(as.matrix(LF[LF[[Flagcol]]==i&LF$dummy==FALSE,fcol:lcol]),LF$weight[LF[[Flagcol]]==i&LF$dummy==FALSE])
  }

  # print to the screen

  cat("\n")
  cat(paste0("******  All results are saved in folder ",save_dir,folder_name,"  ******\n\n"))
  cat(paste0((e0-e)/e0*100,"% variance explained\n\n"))

  # plot result
  Flag <- LF[[Flagcol]]

  # Compare the LF among cells
  if(pdf==FALSE) png(paste0(save_dir,folder_name,"/split(lf).png"),width = 800,height = 800)
  else pdf(paste0(save_dir,folder_name,"/split(lf).pdf"),width = 5,height = 5)

  for (j in sort(unique(LF[[Flagcol]]))) {
    LF_plot <- LF[which(Flag==j),]
    LF_mean <- apply(LF_plot,2,mean)
    Length <- bins
    LF_mean <- LF_mean[fcol:lcol]

    if(j==min(sort(unique(LF[[Flagcol]])))) {
      plot(x=Length,y=LF_mean,pch=toString(j),
           ylim = c(0,max(LF_mean)*lf_scale),
           main = paste0(round((e0-e)/e0*100,2),"% variance explained"),
           cex=2)
      lines(x=Length,y=LF_mean)
    }
    else {
      points(x=Length,y=LF_mean,pch=toString(j),cex=2)
      lines(x=Length,y=LF_mean)
    }
  }

  dev.off()

}
