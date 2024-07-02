#' Visualize the input LF data
#'
#' \code{make.lf.cell} This function plot mean LF by cell
#'
#' @param LF The length frequency data frame input; must include four columns: lat, lon, year, and quarter
#' @param fcol The first column in the data frame with length frequency info
#' @param lcol The last column in the data frame with length frequency info
#' @param bins Names of all bins included in LF
#' @param save_dir The directory where results will be saved
#' @param plot_name The name of the plot to be saved
#' @param plot_format The format (e.g., png and pdf) of the plot to be saved
#' @param width The width of the plot to be saved
#' @param height The height of the plot to be saved
#' @param lwd The width of lf bars
#'
#' @export

make.lf.cell <- function(LF,fcol,lcol,bins,save_dir,plot_name="LF_map",plot_format="png",lwd=1, width = 8, height = 5) {
  # check data first
  if((lcol-fcol+1)!= length(bins)) stop("Error! The number of bins does not match the number of LF columns specified")
  else names(LF)[fcol:lcol] <- bins

  LF_plot <- LF[,c("year","quarter","lat","lon",paste0(bins),"Flag")]
  LF_long <- data.frame(tidyr::gather(LF_plot,fcol:lcol,key = "length",value = "lf"))
  LF_long$length <- as.numeric(LF_long$length)
  LF_long$Flag <- as.factor(LF_long$Flag)

  LF_mean <- dplyr::summarise(dplyr::group_by(LF_long, Flag, length),lf_mean=sum(lf))
  LF_mean <- dplyr::mutate(dplyr::group_by(LF_mean, Flag),lf_mean=lf_mean/sum(lf_mean))

  lf.map <- ggplot2::ggplot(data=LF_mean) +
    ggplot2::geom_line(ggplot2::aes(x=length,y=lf_mean,color=Flag),size=lwd) +
    ggplot2::geom_text(ggplot2::aes(x=length,y=lf_mean,label=Flag,size=2)) +
    ggplot2::theme_bw() +
    ggplot2::xlab("Length (cm)") +
    ggplot2::ylab("Length frquency")

  ggplot2::ggsave(lf.map,filename = paste0(save_dir,plot_name,".",plot_format),width = width, height = height)

  return(lf.map)
  }
