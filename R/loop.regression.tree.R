#' Loop the regression tree analysis with specified number of splits
#'
#' \code{loop.regression.tree} This function loops the run_regression_tree function
#'
#' @param LF The length frequency data frame input; must include four columns: lat, lon, year, and quarter
#' @param fcol The first column in the data frame with length frequency info
#' @param lcol The last column in the data frame with length frequency info
#' @param bins Names of all bins included in LF
#' @param Nsplit The number of splits
#' @param save_dir The directory where results will be saved
#' @param select_matrix User-specified split matrix to be explored
#' @param quarter Whether to consider quarter as a splitting dimension; default = TRUE
#' @param year Whether to consider year as a splitting dimension; default = FALSE
#' @param include_dummy Whether to include dummy data; default = FALSE
#'
#' @return return the input LF with cell number and the percentage of total LF variance explained by each split
#'
#' @export

loop_regression_tree <- function(LF,fcol,lcol,bins,Nsplit,save_dir,select_matrix,lat.min=1,lon.min=1,year.min=1,quarter=TRUE,year=FALSE,include_dummy=FALSE) {

  for (i in 1:nrow(select_matrix)) {
    if(i==1) {
      selecti <- select_matrix[1,]
      LF_loop <- run_regression_tree(LF,fcol,lcol,bins,Nsplit,save_dir,manual = TRUE,select=selecti,lat.min,lon.min,year.min,quarter,year,include_dummy)
      Imp_DF <- c(selecti,LF_loop$Record$Var_explained[Nsplit])
    }
    else {
      selecti <- select_matrix[i,]
      LF_loop <- run_regression_tree(LF,fcol,lcol,bins,Nsplit,save_dir,manual = TRUE,select=selecti,lat.min,lon.min,year.min,quarter,year,include_dummy)
      Imp_DF <- rbind(Imp_DF,c(selecti,LF_loop$Record$Var_explained[Nsplit]))
    }
  }

  Imp_DF <- data.frame(Imp_DF)
  names(Imp_DF) <- c(paste0("select",1:Nsplit),"Var_explained")
  row.names(Imp_DF) <- 1:nrow(Imp_DF)

  Imp_DF_sorted <- Imp_DF[order(Imp_DF$Var,decreasing = TRUE),]

  select <- as.numeric(Imp_DF_sorted[1,1:Nsplit]) # the first row has the highest % variance explained

  LF_Tree <- run_regression_tree(LF,fcol,lcol,bins,Nsplit,save_dir,manual = TRUE, select, quarter=quarter, include_dummy=include_dummy)

  cat("\n\n")
  cat("***************************************************************************************************************************************************************\n")
  cat("Below shows the loop summary results (table saved in loop.csv)\n\n")

  cat("Unsorted:\n")
  print(Imp_DF, row.names = FALSE)

  cat("\nSorted:\n")
  print(Imp_DF_sorted, row.names = FALSE)

  write.csv(Imp_DF_sorted,file=paste0(save_dir,"loop.csv"),row.names = FALSE)

  return(list("LF_Tree"=LF_Tree,"Imp_DF"=Imp_DF,"Imp_DF_sorted"=Imp_DF_sorted))
}
