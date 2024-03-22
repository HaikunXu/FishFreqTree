#' Visualize the input LF data
#'
#' \code{make.lf.map} This function plot the input LF data as maps
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

make.lf.map <- function(LF,fcol,lcol,bins,save_dir,plot_name="LF_map",plot_format="png",width=length(unique(LF$lon)),height=length(unique(LF$lat)),lwd=1) {
  # check data first
  if((lcol-fcol+1)!= length(bins)) stop("Error! The number of bins does not match the number of LF columns specified")
  else names(LF)[fcol:lcol] <- bins

  if(is.null(LF[["weight"]])==TRUE) LF$weight <- 1
  else print("Weight is used to compute the mean length freqeuncy!")

  LF_plot <- LF[,c("year","quarter","lat","lon",paste0(bins),"weight")]
  LF_long <- data.frame(tidyr::gather(LF,fcol:lcol,key = "length",value = "lf"))
  LF_long$length <- as.numeric(LF_long$length)

  # Reverse the lat grid so that negative lat are below positive lat
  LF_long$lat <- factor(LF_long$lat, levels = rev(levels(factor(LF_long$lat))))

  LF_mean <- dplyr::summarise(dplyr::group_by(LF_long, lat, lon, length),lf_mean=sum(lf*weight)/sum(weight))

  lf.map <- ggplot2::ggplot(data=LF_mean) +
    ggplot2::geom_linerange(ggplot2::aes(x=length,ymax=lf_mean,ymin=0),size=lwd) +
    ggplot2::facet_grid(lat~lon) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.spacing = ggplot2::unit(0, "lines")) +
    ggplot2::xlab("Length (cm)") +
    ggplot2::ylab("Length frquency")

  ggplot2::ggsave(lf.map,filename = paste0(save_dir,plot_name,".",plot_format),width = width, height=height)

  return(lf.map)
  }
