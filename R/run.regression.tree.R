#' Run the regression tree analysis with specified number of splits
#'
#' \code{run_regression_tree} This function generates an area code for longline catch allocation
#'
#' @param LF The length frequency data frame input; must include four columns: lat, lon, year, and quarter
#' @param fcol The first column in the data frame with length frequency info
#' @param lcol The last column in the data frame with length frequency info
#' @param bins Names of all bins included in LF
#' @param Nsplit The number of splits
#' @param save_dir The directory where results will be saved
#' @param manual Whether to use user-specified splits; default = FALSE
#' @param select User-specified splits; default = NA
#' @param lat.min Minimal number of lat grids allowed for a cell, which could not be split beyond that
#' @param lon.min Minimal number of lon grids allowed for a cell, which could not be split beyond that
#' @param year.min Minimal number of years allowed for a cell, which could not be split beyond that
#' @param quarter Whether to consider quarter as a splitting dimension; default = TRUE
#' @param year Whether to consider year as a splitting dimension; default = FALSE
#' @param include_dummy Whether to include dummy data; default = FALSE
#'
#' @return return the input LF with cell number and the percentage of total LF variance explained by each split
#'
#' @export

run_regression_tree <- function(LF,fcol,lcol,bins,Nsplit,save_dir,manual = FALSE,select=NA,lat.min=1,lon.min=1,year.min=1,quarter=TRUE,year=FALSE,include_dummy=FALSE) {

  if(manual==FALSE) select <- rep(1,Nsplit) # every split choose the one with the max improvement
  if(include_dummy==FALSE) LF$dummy <- FALSE
  if(manual==TRUE&(sum(select>1)>1)) print("Warning!!! You selection is hierarchical.")

  # make sure bins are consistent with LF input
  if((lcol-fcol+1)!= length(bins)) stop("Error! The number of bins does not match the number of LF columns specified")
  else names(LF)[fcol:lcol] <- bins

  # make sure LF sums to 1 for each row
  # row_sum <- apply(LF[,fcol:lcol],1,sum)
  # print(row_sum)
  # if(sum(row_sum>1)>0) stop("Error! LF does not sum to 1 for at least one row.")

  var_exp <- rep(NA,Nsplit) # variance explained
  folder_name <- paste0(gsub(", ","",toString(select)))
  select_name <- paste0(gsub(", ","",toString(select)),"/")
  unlink(paste0(save_dir,folder_name), recursive = TRUE)
  dir.create(paste0(save_dir,select_name))

  plot_year <- FALSE

  for (i in 1:Nsplit) {

    if(i==1) {
      # whole ALB area
      split <- find_split(LF[LF$dummy==FALSE,],fcol,lcol,lat.min,lon.min,year.min,quarter,year)
      # save result as a csv.file
      write.csv(rename_CQrt(split),file=paste0(save_dir,select_name,"split",i,".csv"),row.names = FALSE)

      # if((manual==FALSE)|(i %in% user_split$Number == FALSE))
      LF <- make.Us.areaflags.f(LF, as.character(split$Key[select[i]]), as.numeric(split$Value[select[i]]),1,0)

      # LF <- make.Us.areaflags.f(LF, user_split$Key[which(user_split$Number==i)], user_split$Value[which(user_split$Number==i)],1,0)

      e0 <- get.klderror.null(as.matrix(LF[LF$dummy==FALSE,fcol:lcol])) # null (no stratification)
      e1 <- get.klderror.null(as.matrix(LF[LF$Flag1==1&LF$dummy==FALSE,fcol:lcol]))
      e2 <- get.klderror.null(as.matrix(LF[LF$Flag1==2&LF$dummy==FALSE,fcol:lcol]))

      var_exp[1] <- (e0-e1-e2)/e0

      # print to the screen
      cat("\n\n")
      cat("***Note***: below shows the best splits in order, please check the save figures under the directory save_dir to better understand the meaning of each split\n\n")
      cat(paste0("******  Results are saved in folder ",save_dir,select_name,"  ******\n\n"))

      if(select[i]==1)
        print(paste0("Best 1st split: ",split[1,2],"<=",split[1,3]))
      else
        print(paste0("Manual 1st split: ",split[select[i],2],"<=",split[select[i],3]))
      cat(paste0((e0-e1-e2)/e0*100,"% variance explained\n\n"))

      # plot result
      Flag <- LF[[paste0("Flag",i)]]
      Quarter <- LF$quarter
      xlim <- c(min(LF$lon),max(LF$lon))
      ylim <- c(min(LF$lat),max(LF$lat))
      tlim <- c(min(LF$year),max(LF$year))

      # check whether are year splits; if yes, plot by year blocks
      year_block <- data.frame(Cell = 1:(i+1), ymin = 0, ymax = 0, id = 0)

      for (j in 1:(i+1)) {
        LF_plot <- LF[which(LF$Flag1==j),]
        year_block[j,2:3] <- c(min(LF_plot$year), max(LF_plot$year))
      }

      year_block_id <- data.frame(block=sort(unique(year_block$ymin+year_block$ymax)),
                                  id=rank(unique(year_block$ymin+year_block$ymax)))

      for (j in 1:(i+1)) {
        year_block$id[j] = year_block_id$id[which(year_block_id$block==(year_block$ymin[j]+year_block$ymax[j]))]
      }

      if(split[select[i],2]=="Year") plot_year <- TRUE
      else year_block$id <- 1

      for (y in 1:length(unique(year_block$id))) {
        # plot the spatial distribution of cells
        yearmin <- min(year_block$ymin[which(year_block$id==y)])
        yearmax <- max(year_block$ymax[which(year_block$id==y)])

        png(paste0(save_dir,select_name,"split",i,"(map for years ",yearmin,"-",yearmax,").png"),width = 1000,height = 1000)
        par(mfrow=c(2,2))

        for (q in 1:4) {
          for (j in which(year_block$id==y)) {
            LF_plot <- LF[which(Flag==j&Quarter==q),]
            if(j==min(which(year_block$id==y))) {
              plot(x=LF_plot$lon,y=LF_plot$lat,pch=toString(j),
                   xlim = xlim, ylim = ylim, xlab = "Lon", ylab = "Lat",
                   main = paste0("(Quarter ",q,")"," Split#",i,": ",split[select[i],2],"<=",split[select[i],3]))
            }
            else {
              points(x=LF_plot$lon,y=LF_plot$lat,pch=toString(j))
            }
          }
        }

        dev.off()
      }

      # Compare the LF among cells
      png(paste0(save_dir,select_name,"split",i,"(lf).png"),width = 500,height = 500)

      for (j in 1:(i+1)) {
        LF_plot <- LF[which(Flag==j),]
        LF_mean <- apply(LF_plot,2,mean)
        Length <- bins
        LF_mean <- LF_mean[fcol:lcol]

        if(j==1) {
          plot(x=Length,y=LF_mean,pch=toString(j),
               ylim = c(0,max(LF_mean)),
               main = paste0(round((e0-e1-e2)/e0*100,2),"% variance explained"))
          lines(x=Length,y=LF_mean)
        }
        else {
          points(x=Length,y=LF_mean,pch=toString(j))
          lines(x=Length,y=LF_mean)
        }
      }

      dev.off()

      # Compare Improvement against lat and lon
      png(paste0(save_dir,select_name,"split",i,"(latlon).png"),width = 1000,height = 500)
      par(mfrow=c(1,2))

      ymax <- max(split$Improve[which(split$Key=="Lat"|split$Key=="Lon")])

      split_plot <- split[which(split$Key=="Lat"),]
      plot(x=split_plot$Value,y=split_plot$Improve,
           xlim = ylim, ylim = c(0, ymax), xlab = "Lat", ylab = "Improvement",
           main = paste0(" Split#",i))

      split_plot <- split[which(split$Key=="Lon"),]
      plot(x=split_plot$Value,y=split_plot$Improve,
           xlim = xlim, ylim = c(0, ymax), xlab = "Lon", ylab = "Improvement",
           main = paste0(" Split#",i))

      dev.off()

      # Compare Improvement against year
      if(year==TRUE) {
        png(paste0(save_dir,select_name,"split",i,"(year).png"),width = 500,height = 500)

        split_plot <- split[which(split$Key=="Year"),]
        plot(x=split_plot$Value,y=split_plot$Improve,
             xlim = tlim, xlab = "Year", ylab = "Improvement",
             main = paste0(" Split#",i))

        dev.off()
      }

    }

    else {
      # imp <- rep(NA,i)
      # loop for every possible split
      for (ii in 1:i) {
        LF_raw <- LF[LF$dummy==FALSE,]
        j <- which(LF_raw[[paste0("Flag",i-1)]] == ii)
        LF_data <- LF_raw[j,]

        split <- find_split(LF_data,fcol,lcol,lat.min,lon.min,year.min,quarter,year)
        split$Cell <- ii
        if(ii==1) split_raw <- split
        else split_raw <- rbind(split_raw,split)

        LF_raw <- make.Us.areaflags.f(LF_raw,
                                      as.character(split$Key[1]),
                                      as.numeric(split$Value[1]),i,ii)

        for (k in 1:(i+1)) {
          if(k==1) e <- get.klderror.null(as.matrix(LF_raw[LF_raw[[paste0("Flag",i)]]==1,fcol:lcol]))
          else e <- e + get.klderror.null(as.matrix(LF_raw[LF_raw[[paste0("Flag",i)]]==k,fcol:lcol]))
        }

        # imp[ii] <- (e0-e)/e0

      }

      # j <- which(LF[[paste0("Flag",i-1)]] == ii)
      # LF_data <- LF[j,]
      # split <- find_split(LF_data,fcol,lcol,lat.min,lon.min)

      # save result as a csv.file

      split_raw <- split_raw[order(split_raw$Improve,decreasing = TRUE),]
      write.csv(rename_CQrt(split_raw),file=paste0(save_dir,select_name,"split",i,".csv"),row.names = FALSE)

      # the Cell with the most improvement
      ii <- split_raw$Cell[select[i]]

      # if((manual==FALSE)|(i %in% user_split$Number == FALSE))

      LF <- make.Us.areaflags.f(LF, as.character(split_raw$Key[select[i]]), as.numeric(split_raw$Value[select[i]]),i,ii)

      # else
      # LF <- make.Us.areaflags.f(LF, user_split$Key[which(user_split$Number==i)], user_split$Value[which(user_split$Number==i)],i,ii)

      for (k in 1:(i+1)) {
        if(k==1) e <- get.klderror.null(as.matrix(LF[LF[[paste0("Flag",i)]]==1&LF$dummy==FALSE,fcol:lcol]))
        else e <- e + get.klderror.null(as.matrix(LF[LF[[paste0("Flag",i)]]==k&LF$dummy==FALSE,fcol:lcol]))
      }

      var_exp[i] <- (e0-e)/e0

      # print result to screen
      if(select[i]==1) {
        if(i==2)
          print(paste0("Conditional best 2nd split is for cell ",ii," in split",i-1,".png: ",split_raw[1,2],"<=",split_raw[1,3]))
        if(i==3)
          print(paste0("Conditional best 3rd split is for cell ",ii," in split",i-1,".png: ",split_raw[1,2],"<=",split_raw[1,3]))
        if(i>3)
          print(paste0("Conditional best ",i,"th"," split is for cell ",ii," in split",i-1,".png: ",split_raw[1,2],"<=",split_raw[1,3]))
      }
      else {
        if(i==2)
          print(paste0("Manual 2nd split is for cell ",ii," in split",i-1,".png: ",split_raw[select[i],2],"<=",split_raw[select[i],3]))
        if(i==3)
          print(paste0("Manual 3rd split is for cell ",ii," in split",i-1,".png: ",split_raw[select[i],2],"<=",split_raw[select[i],3]))
        if(i>3)
          print(paste0("Manual ",i,"th"," split is for cell ",ii," in split",i-1,".png: ",split_raw[select[i],2],"<=",split_raw[select[i],3]))
      }

      cat(paste0((e0-e)/e0*100,"% variance explained\n\n"))

      # plot result
      Flag <- LF[[paste0("Flag",i)]]
      Quarter <- LF$quarter

      # png(paste0(save_dir,select_name,"split",i,".png"),width = 1000,height = 500)
      # par(mfrow=c(1,2))

      # plot the spatial distribution of cells

      # check whether are year splits; if yes, plot by year blocks
      year_block <- data.frame(Cell = 1:(i+1), ymin = 0, ymax = 0, id = 0)

      for (j in 1:(i+1)) {
        LF_plot <- LF[which(LF[[paste0("Flag",i)]]==j),]
        year_block[j,2:3] <- c(min(LF_plot$year), max(LF_plot$year))
      }

      year_block_id <- data.frame(block=sort(unique(year_block$ymin+year_block$ymax)),
                                  id=rank(unique(year_block$ymin+year_block$ymax)))

      for (j in 1:(i+1)) {
        year_block$id[j] = year_block_id$id[which(year_block_id$block==(year_block$ymin[j]+year_block$ymax[j]))]
      }

      if(split_raw$Key[select[i]]=="Year") plot_year <- TRUE
      if(plot_year==FALSE) year_block$id <- 1

      for (y in 1:length(unique(year_block$id))) {

        # plot the spatial distribution of cells
        yearmin <- min(year_block$ymin[which(year_block$id==y)])
        yearmax <- max(year_block$ymax[which(year_block$id==y)])

        png(paste0(save_dir,select_name,"split",i,"(map for years ",yearmin,"-",yearmax,").png"),width = 1000,height = 1000)
        par(mfrow=c(2,2))

        for (q in 1:4) {
          for (j in which(year_block$id==y)) {
            LF_plot <- LF[which(Flag==j&Quarter==q),]
            if(j==min(which(year_block$id==y))) {
              plot(x=LF_plot$lon,y=LF_plot$lat,pch=toString(j),
                   xlim = xlim, ylim = ylim, xlab = "Lon", ylab = "Lat",
                   main = paste0("(Quarter ",q,")"," Split#",i,": ",split[select[i],2],"<=",split[select[i],3]))
            }
            else {
              points(x=LF_plot$lon,y=LF_plot$lat,pch=toString(j))
            }
          }
        }

        dev.off()
      }

      # Compare the LF among cells
      png(paste0(save_dir,select_name,"split",i,"(lf).png"),width = 500,height = 500)

      for (j in 1:(i+1)) {
        LF_plot <- LF[which(Flag==j),]
        LF_mean <- apply(LF_plot,2,mean)
        Length <- bins
        LF_mean <- LF_mean[fcol:lcol]

        if(j==1) {
          plot(x=Length,y=LF_mean,pch=toString(j),
               ylim = c(0,max(LF_mean)),
               main = paste0(round((e0-e)/e0*100,2),"% variance explained"))
          lines(x=Length,y=LF_mean)
        }
        else {
          points(x=Length,y=LF_mean,pch=toString(j))
          lines(x=Length,y=LF_mean)
        }
      }

      dev.off()

      # Compare Improvement against lat and lon
      png(paste0(save_dir,select_name,"split",i,"(latlon).png"),width = 1000,height = 500)
      par(mfrow=c(1,2))

      ymax <- max(split_raw$Improve[which(split_raw$Key=="Lat"|split_raw$Key=="Lon")])

      for (j in 1:i) {

        split_plot <- split_raw[which(split_raw$Key=="Lat"&split_raw$Cell==j),]
        if(j==1) {
          plot(x=split_plot$Value,y=split_plot$Improve,pch=toString(j),
               xlim = ylim, ylim = c(0,ymax), xlab = "Lat", ylab = "Improvement",
               main = paste0(" Split#",i))
        }
        else points(x=split_plot$Value,y=split_plot$Improve,pch=toString(j))
      }
      for (j in 1:i) {
        split_plot <- split_raw[which(split_raw$Key=="Lon"&split_raw$Cell==j),]
        if(j==1) {
          plot(x=split_plot$Value,y=split_plot$Improve,pch=toString(j),
               xlim = xlim, ylim = c(0,ymax), xlab = "Lon", ylab = "Improvement",
               main = paste0(" Split#",i))
        }
        else points(x=split_plot$Value,y=split_plot$Improve,pch=toString(j))
      }
      dev.off()

      # Compare Improvement against year
      if(year==TRUE) {
        png(paste0(save_dir,select_name,"split",i,"(year).png"),width = 500,height = 500)

        for (j in 1:i) {

          split_plot <- split_raw[which(split_raw$Key=="Year"&split_raw$Cell==j),]
          if(j==1) {
            plot(x=split_plot$Value,y=split_plot$Improve,pch=toString(j),
                 xlim = tlim, xlab = "Year", ylab = "Improvement",
                 main = paste0(" Split#",i))
          }
          else points(x=split_plot$Value,y=split_plot$Improve,pch=toString(j))
        }

        dev.off()

      }
    }

  }

  return(list("LF"=LF,"Var"=var_exp))
}
