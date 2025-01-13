#' Plot mean length map based on the input LF data
#'
#' \code{make.meanl.map} This function plot mean length as maps
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
#' @param s The size of squares
#'
#' @export

make.meanl.map <- function(LF,fcol,lcol,bins,save_dir,plot_name="MeanL_map",plot_format="png",width=length(unique(LF$lon))/2,height=length(unique(LF$lat))/2) {
  # check data first
  if((lcol-fcol+1)!= length(bins)) stop("Error! The number of bins does not match the number of LF columns specified")
  else names(LF)[fcol:lcol] <- bins

  if(is.null(LF[["weight"]])==TRUE) LF$weight <- 1
  else print("Weight is used to compute the mean length!")

  LF_plot <- LF[,c("year","quarter","lat","lon",paste0(bins),"weight")]
  LF_long <- data.frame(tidyr::gather(LF,fcol:lcol,key = "length",value = "lf"))
  LF_long$length <- as.numeric(LF_long$length)

  L_mean <- dplyr::summarise(dplyr::group_by(LF_long, lat, lon),mean_length=sum(lf*length*weight)/sum(lf*weight))

  write.csv(L_mean, file = paste0(save_dir,"MeanL.csv"), row.names = FALSE)

  wmap <- ggplot2::map_data("world")

  meanl.map <- ggplot2::ggplot(data=L_mean) +
    ggplot2::geom_tile(ggplot2::aes(x=lon,y=lat,fill=mean_length),color="black") +
    ggplot2::scale_fill_distiller(palette = "Spectral") +
    ggplot2::theme_bw() +
    ggplot2::xlab("Lon") +
    ggplot2::ylab("Lat") +
    ggplot2::geom_polygon(data=wmap,ggplot2::aes(long, lat, group = group),fill = "black",colour = "white",alpha = 1,lwd=0.5) +
    ggplot2::coord_quickmap(ylim = c(min(L_mean$lat) - 5, max(L_mean$lat) + 5),
                            xlim = c(min(L_mean$lon) - 5, max(L_mean$lon) + 5))

  ggplot2::ggsave(meanl.map,filename = paste0(save_dir,plot_name,".",plot_format),width = width, height=height)

  return(meanl.map)
  }
