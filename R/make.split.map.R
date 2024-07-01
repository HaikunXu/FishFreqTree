#' Map the split outcome
#'
#' \code{make.split.map} This function color-codes cells based on the regression tree outcome
#'
#' @param LF The length frequency data frame input; must include four columns: lat, lon, year, and quarter
#' @param Nsplit The number of splits
#' @param save_dir The directory where results will be saved
#' @param plot_name The name of the plot to be saved
#' @param plot_format The format (e.g., png and pdf) of the plot to be saved
#' @param width The width of the plot to be saved
#' @param height The height of the plot to be saved
#' @param s The size of squares
#'
#' @export

make.split.map <- function(LF,Nsplit,save_dir,plot_name="Split_map",plot_format="png",width=length(unique(LF$lon))/2,height=length(unique(LF$lat))/2,text=FALSE,s=10) {
  wmap <- ggplot2::map_data("world")

  LF$Cell <- factor(LF[[paste0("Flag", Nsplit)]])

  if (text == FALSE) {
    split.map <- ggplot2::ggplot(data = LF) +
      ggplot2::geom_tile(ggplot2::aes(x = lon, y = lat, fill = Cell), color = "black") +
      ggplot2::geom_polygon(
        data = wmap,
        ggplot2::aes(long, lat, group = group),
        fill = "black",
        colour = "white",
        alpha = 1,
        lwd = 0.5
      ) +
      ggplot2::theme_bw() +
      ggplot2::xlab("Lon") +
      ggplot2::ylab("Lat") +
      ggplot2::coord_quickmap(ylim = c(min(LF$lat), max(LF$lat)), xlim = c(min(LF$lon), max(LF$lon)))
  }
  else {
    split.map <- ggplot2::ggplot(data = LF) +
      ggplot2::geom_polygon(
        data = wmap,
        ggplot2::aes(long, lat, group = group),
        fill = "black",
        colour = "white",
        alpha = 1,
        lwd = 0.5
      ) +
      ggplot2::geom_text(ggplot2::aes(x = lon, y = lat, label = Cell, color = Cell),
                         size = s) +
      ggplot2::theme_bw() +
      ggplot2::xlab("Lon") +
      ggplot2::ylab("Lat") +
      ggplot2::coord_quickmap(ylim = c(min(LF$lat), max(LF$lat)), xlim = c(min(LF$lon), max(LF$lon)))
  }

  ggplot2::ggsave(
    split.map,
    filename = paste0(save_dir, plot_name, ".", plot_format),
    width = width,
    height = height
  )

  return(split.map)
  }
