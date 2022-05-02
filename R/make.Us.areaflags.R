#' Run the regression tree analysis with specified number of splits
#'
#' \code{make.Us.areaflags.f} This function generates an area code for longline catch allocation
#'
#' @export

make.Us.areaflags.f <- function(lf.input,key,value,split_num,area_num)
{
  # area flags for input data frame (lf.input) (assuming 5 deg lat x 10 deg lon resolution)
  # input data can be either l-f or cpue, if variable names work
  # strat.flgs<-matrix(NA,ncol=7,nrow=nrow(lf.input))

  #
  if(key=="Lat") Flag <- ifelse(lf.input$lat<=value,1,2)
  if(key=="Lon") Flag <- ifelse(lf.input$lon<=value,1,2)
  if(key=="Qrt") Flag <- ifelse(lf.input$quarter<=value,1,2)
  if(key=="CQrt") {
    if(value==1) Flag <- ifelse((lf.input$quarter==2)|(lf.input$quarter==3),1,2)
    if(value==2) Flag <- ifelse(lf.input$quarter==3,1,2)
    if(value==3) Flag <- ifelse(lf.input$quarter==2,1,2)
  }
  if(key=="Year") Flag <- ifelse(lf.input$year<=value,1,2)

  if(split_num==1) lf.input[["Flag1"]] <- Flag
  else {
    Flag_new <- lf.input[[paste0("Flag",split_num-1)]]

    Flag_new[which(Flag_new==area_num)] <- paste0(
      Flag_new[which(Flag_new==area_num)],Flag[which(Flag_new==area_num)])

    # print(Flag_new)

    lf.input[[paste0("Flag",split_num)]] <- as.numeric(factor(Flag_new))
  }
  #
  return(lf.input)
}
