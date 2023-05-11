#' Run the regression tree analysis with specified number of splits
#'
#' \code{run.regression.tree} This function is the main regression tree function
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
#' @param pdf Whether to save figures in pdf - default is in png
#'
#' @return return the input LF with cell number and the percentage of total LF variance explained by each split
#'
#' @export

run_regression_tree <- function(LF,fcol,lcol,bins,Nsplit,save_dir,manual = FALSE,select=NA,lat.min=1,lon.min=1,year.min=1,quarter=TRUE,year=FALSE,include_dummy=FALSE,pdf=FALSE,lf_scale=1.5) {

  print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
  print("WARNING: lat and lon in the input data should be for grid centers")
  print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")

  if(manual==FALSE) select <- rep(1,Nsplit) # every split choose the one with the max improvement
  if(include_dummy==FALSE) LF$dummy <- FALSE
  if(manual==TRUE&(sum(select>1)>1)) print("Warning!!! You selection is hierarchical.")

  # make sure bins are consistent with LF input
  if((lcol-fcol+1)!= length(bins)) stop("Error! The number of bins does not match the number of LF columns specified")
  else names(LF)[fcol:lcol] <- bins

  # make sure LF sums to 1 for each row
  row_sum <- apply(LF[which(LF$dummy==FALSE),fcol:lcol],1,sum)
  if(sum(abs(row_sum-1)>0.05)>0) {
    plot(row_sum)
    # print(which(abs(row_sum-1)>0.1))
    stop("Error! LF does not sum to 1 for at least one row.")
  }

  if(is.null(LF[["weight"]])==TRUE) LF$weight <- 1
  else print("Weight is used to compute the percentage of varaibility explained!")

  Record <- data.frame(Key=rep(NA,Nsplit),Value=rep(NA,Nsplit),Cell=rep(NA,Nsplit),Var_explained=rep(NA,Nsplit))
  row.names(Record) <- paste0("Split",1:Nsplit)
  folder_name <- paste0(gsub(", ","",toString(select)))
  select_name <- paste0(gsub(", ","",toString(select)),"/")
  unlink(paste0(save_dir,folder_name), recursive = TRUE)
  dir.create(paste0(save_dir,select_name))

  plot_year <- FALSE

  for (i in 1:Nsplit) {

    if(i==1) {
      # whole area
      split <- find_split(LF[LF$dummy==FALSE,],fcol,lcol,lat.min,lon.min,year.min,quarter,year)
      split$Var_explained <- 0

      e0 <- get.klderror.null(as.matrix(LF[LF$dummy==FALSE,fcol:lcol]),LF$weight[LF$dummy==FALSE]) # null (no stratification)

      for (sp in 1:nrow(split)) {
        # print(split)
        if(split$Improvement[sp]>0) {
          LF_raw <- make.Us.areaflags.f(LF[LF$dummy==FALSE,],
                                        as.character(split$Key[sp]),
                                        as.numeric(split$Value[sp]),i,ii)

          e <- get.klderror.null(as.matrix(LF_raw[LF_raw[[paste0("Flag",i)]]==1,fcol:lcol]),LF_raw$weight[LF_raw[[paste0("Flag",i)]]==1]) +
            get.klderror.null(as.matrix(LF_raw[LF_raw[[paste0("Flag",i)]]==2,fcol:lcol]),LF_raw$weight[LF_raw[[paste0("Flag",i)]]==2])

          split$Var_explained[sp] <- (e0-e)/e0
        }
      }

      split_save <- rename_CQrt(split)
      split_save$Rank <- seq(1,nrow(split))
      # save result as a csv.file
      write.csv(split_save[,2:5],file=paste0(save_dir,select_name,"split",i,".csv"),row.names = FALSE)
      write.csv(split_save[,c(2,3,1)],file=paste0(save_dir,select_name,"improvement-split",i,".csv"),row.names = FALSE)

      # if((manual==FALSE)|(i %in% user_split$Number == FALSE))

      LF <- make.Us.areaflags.f(LF, as.character(split$Key[select[i]]), as.numeric(split$Value[select[i]]),1,0)

      # LF <- make.Us.areaflags.f(LF, user_split$Key[which(user_split$Number==i)], user_split$Value[which(user_split$Number==i)],1,0)

      e1 <- get.klderror.null(as.matrix(LF[LF$Flag1==1&LF$dummy==FALSE,fcol:lcol]),LF$weight[LF$Flag1==1&LF$dummy==FALSE])
      e2 <- get.klderror.null(as.matrix(LF[LF$Flag1==2&LF$dummy==FALSE,fcol:lcol]),LF$weight[LF$Flag1==2&LF$dummy==FALSE])

      Record[1,4] <- (e0-e1-e2)/e0
      Record[1,2] <- split[select[i],3]
      Record[1,1] <- as.character(split[select[i],2])
      # Record[1,4] <- as.character(split[select[i],2])

      # print to the screen
      cat("\n")
      cat(paste0("***Note***: Below shows the best splits in order, please check saved figures to better understand the meaning of each split\n"))
      cat(paste0("******  All results are saved in folder ",save_dir,select_name,"  ******\n\n"))

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

        # quarterly maps
        if(quarter==TRUE) {
        if(pdf==FALSE) png(paste0(save_dir,select_name,"split",i,"(quarterly maps for years ",yearmin,"-",yearmax,").png"),width = 1500,height = 1500)
        else pdf(paste0(save_dir,select_name,"split",i,"(quarterly maps for years ",yearmin,"-",yearmax,").pdf"),width = 10,height = 10)

        par(mfrow=c(2,2))

        for (q in 1:4) {
          for (j in which(year_block$id==y)) {
            LF_plot <- LF[which(Flag==j&Quarter==q),]
            if(j==min(which(year_block$id==y))) {
              plot(x=LF_plot$lon,y=LF_plot$lat,pch=toString(j),cex=2,
                   xlim = xlim, ylim = ylim, xlab = "Lon", ylab = "Lat",
                   main = paste0("(Quarter ",q,")"," Split#",i,": ",split[select[i],2],"<=",split[select[i],3]))
            }
            else {
              points(x=LF_plot$lon,y=LF_plot$lat,pch=toString(j),cex=2)
            }
          }
        }

        dev.off()
        }

        # annual maps
        if(pdf==FALSE) png(paste0(save_dir,select_name,"split",i,"(annual maps for years ",yearmin,"-",yearmax,").png"),width = 800,height = 800)
        else pdf(paste0(save_dir,select_name,"split",i,"(annual maps for years ",yearmin,"-",yearmax,").pdf"),width = 5,height = 5)

          for (j in which(year_block$id==y)) {
            LF_plot <- LF[which(Flag==j),]
            if(j==min(which(year_block$id==y))) {
              plot(x=LF_plot$lon,y=LF_plot$lat,pch=toString(j),cex=2,
                   xlim = xlim, ylim = ylim, xlab = "Lon", ylab = "Lat",
                   main = paste0("Split#",i,": ",split[select[i],2],"<=",split[select[i],3]))
            }
            else {
              points(x=LF_plot$lon,y=LF_plot$lat,pch=toString(j),cex=2)
            }
          }

        dev.off()
      }

      # Compare the LF among cells
      if(pdf==FALSE) png(paste0(save_dir,select_name,"split",i,"(lf).png"),width = 800,height = 800)
      else pdf(paste0(save_dir,select_name,"split",i,"(lf).pdf"),width = 5,height = 5)

      for (j in 1:(i+1)) {
        LF_plot <- LF[which(Flag==j),]
        LF_mean <- apply(LF_plot,2,mean)
        Length <- bins
        LF_mean <- LF_mean[fcol:lcol]

        if(j==1) {
          plot(x=Length,y=LF_mean,pch=toString(j),
               ylim = c(0,max(LF_mean)*lf_scale),
               main = paste0(round((e0-e1-e2)/e0*100,2),"% variance explained"),
               cex=2)
          lines(x=Length,y=LF_mean)
        }
        else {
          points(x=Length,y=LF_mean,pch=toString(j),cex=2)
          lines(x=Length,y=LF_mean)
        }
      }

      dev.off()

      # Compare Improvement against lat and lon
      if(pdf==FALSE) png(paste0(save_dir,select_name,"split",i,"(latlon).png"),width = 1000,height = 500)
      else pdf(paste0(save_dir,select_name,"split",i,"(latlon).pdf"),width = 10,height = 5)

      par(mfrow=c(1,2))

      # ymax <- max(split$Improve[which(split$Key=="Lat"|split$Key=="Lon")])

      split_plot <- split[which(split$Key=="Lat"),]
      plot(x=split_plot$Value,y=split_plot$Improve/max(split_plot$Improve),
           xlim = ylim, ylim = c(0, 1), xlab = "Lat", ylab = "Imp/max(Imp)",
           main = paste0(" Split#",i))

      split_plot <- split[which(split$Key=="Lon"),]
      plot(x=split_plot$Value,y=split_plot$Improve/max(split_plot$Improve),
           xlim = xlim, ylim = c(0, 1), xlab = "Lon", ylab = "Imp/max(Imp)",
           main = paste0(" Split#",i))

      dev.off()

      # Compare Improvement against year
      if(year==TRUE) {
        if(pdf==FALSE) png(paste0(save_dir,select_name,"split",i,"(year).png"),width = 800,height = 800)
        else pdf(paste0(save_dir,select_name,"split",i,"(year).pdf"),width = 5,height = 5)

        split_plot <- split[which(split$Key=="Year"),]
        plot(x=split_plot$Value,y=split_plot$Improve/max(split_plot$Improve),
             xlim = tlim, xlab = "Year", ylab = "Imp/max(Imp)",
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
        split$Var_explained <- rep(0,nrow(split))

        for (sp in 1:nrow(split)) {
          if(split$Improvement[sp]>0) {
            LF_raw <- make.Us.areaflags.f(LF_raw,
                                          as.character(split$Key[sp]),
                                          as.numeric(split$Value[sp]),i,ii)

            for (k in 1:(i+1)) {
              if(k==1) e <- get.klderror.null(as.matrix(LF_raw[LF_raw[[paste0("Flag",i)]]==1,fcol:lcol]),LF_raw$weight[LF_raw[[paste0("Flag",i)]]==1])
              else e <- e + get.klderror.null(as.matrix(LF_raw[LF_raw[[paste0("Flag",i)]]==k,fcol:lcol]),LF_raw$weight[LF_raw[[paste0("Flag",i)]]==k])
            }

            split$Var_explained[sp] <- (e0-e)/e0
          }
        }

        if(ii==1) split_raw <- split
        else split_raw <- rbind(split_raw,split)
      }

      # j <- which(LF[[paste0("Flag",i-1)]] == ii)
      # LF_data <- LF[j,]
      # split <- find_split(LF_data,fcol,lcol,lat.min,lon.min)

      # save result as a csv.file

      write.csv(rename_CQrt(split_raw[,c(4,2,3,1)]),file=paste0(save_dir,select_name,"improvement-split",i,".csv"),row.names = FALSE)

      split_raw <- split_raw[order(split_raw$Var,decreasing = TRUE),]
      split_raw$Rank <- seq(1,nrow(split_raw))
      write.csv(rename_CQrt(split_raw[,c(4,2,3,5,6)]),file=paste0(save_dir,select_name,"split",i,".csv"),row.names = FALSE)

      # the Cell with the most improvement
      ii <- split_raw$Cell[select[i]]

      # if((manual==FALSE)|(i %in% user_split$Number == FALSE))

      LF <- make.Us.areaflags.f(LF, as.character(split_raw$Key[select[i]]), as.numeric(split_raw$Value[select[i]]),i,ii)

      # else
      # LF <- make.Us.areaflags.f(LF, user_split$Key[which(user_split$Number==i)], user_split$Value[which(user_split$Number==i)],i,ii)

      for (k in 1:(i+1)) {
        if(k==1) e <- get.klderror.null(as.matrix(LF[LF[[paste0("Flag",i)]]==1&LF$dummy==FALSE,fcol:lcol]),LF$weight[LF[[paste0("Flag",i)]]==1&LF$dummy==FALSE])
        else e <- e + get.klderror.null(as.matrix(LF[LF[[paste0("Flag",i)]]==k&LF$dummy==FALSE,fcol:lcol]),LF$weight[LF[[paste0("Flag",i)]]==k&LF$dummy==FALSE])
      }

      Record[i,4] <- (e0-e)/e0
      Record[i,2] <- split_raw[select[i],3]
      Record[i,1] <- as.character(split_raw[select[i],2])
      Record[i,3] <- ii

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

        # quarter maps
        if(quarter==TRUE) {
        if(pdf==FALSE) png(paste0(save_dir,select_name,"split",i,"(quarterly map for years ",yearmin,"-",yearmax,").png"),width = 1500,height = 1500)
        else pdf(paste0(save_dir,select_name,"split",i,"(quarterly map for years ",yearmin,"-",yearmax,").pdf"),width = 10,height = 10)

        par(mfrow=c(2,2))

        for (q in 1:4) {
          for (j in which(year_block$id==y)) {
            LF_plot <- LF[which(Flag==j&Quarter==q),]
            if(j==min(which(year_block$id==y))) {
              plot(x=LF_plot$lon,y=LF_plot$lat,pch=toString(j),cex=2,
                   xlim = xlim, ylim = ylim, xlab = "Lon", ylab = "Lat",
                   main = paste0("(Quarter ",q,")"," Split#",i,": ",split_raw[select[i],2],"<=",split_raw[select[i],3]))
            }
            else {
              points(x=LF_plot$lon,y=LF_plot$lat,pch=toString(j),cex=2)
            }
          }
        }

        dev.off()
        }

        # annual maps
        if(pdf==FALSE) png(paste0(save_dir,select_name,"split",i,"(annual map for years ",yearmin,"-",yearmax,").png"),width = 800,height = 800)
        else pdf(paste0(save_dir,select_name,"split",i,"(annual map for years ",yearmin,"-",yearmax,").pdf"),width = 5,height = 5)

        for (j in which(year_block$id==y)) {
            LF_plot <- LF[which(Flag==j),]
            if(j==min(which(year_block$id==y))) {
              plot(x=LF_plot$lon,y=LF_plot$lat,pch=toString(j),cex=2,
                   xlim = xlim, ylim = ylim, xlab = "Lon", ylab = "Lat",
                   main = paste0("Split#",i,": ",split_raw[select[i],2],"<=",split_raw[select[i],3]))
            }
            else {
              points(x=LF_plot$lon,y=LF_plot$lat,pch=toString(j),cex=2)
            }
      }
      dev.off()
      }

      # Compare the LF among cells
      if(pdf==FALSE) png(paste0(save_dir,select_name,"split",i,"(lf).png"),width = 800,height = 800)
      else pdf(paste0(save_dir,select_name,"split",i,"(lf).pdf"),width = 5,height = 5)

      for (j in 1:(i+1)) {
        LF_plot <- LF[which(Flag==j),]
        LF_mean <- apply(LF_plot,2,mean)
        Length <- bins
        LF_mean <- LF_mean[fcol:lcol]

        if(j==1) {
          plot(x=Length,y=LF_mean,pch=toString(j),
               ylim = c(0,max(LF_mean)*lf_scale),
               main = paste0(round((e0-e)/e0*100,2),"% variance explained"),
               cex = 2)
          lines(x=Length,y=LF_mean)
        }
        else {
          points(x=Length,y=LF_mean,pch=toString(j),cex=2)
          lines(x=Length,y=LF_mean)
        }
      }

      dev.off()

      # Compare Improvement against lat and lon
      if(pdf==FALSE) png(paste0(save_dir,select_name,"split",i,"(latlon).png"),width = 1000,height = 500)
      else pdf(paste0(save_dir,select_name,"split",i,"(latlon).pdf"),width = 10,height = 5)

      par(mfrow=c(1,2))

      # ymax <- max(split_raw$Improve[which(split_raw$Key=="Lat"|split_raw$Key=="Lon")])

      for (j in 1:i) {

        split_plot <- split_raw[which(split_raw$Key=="Lat"&split_raw$Cell==j),]
        if(j==1) {
          plot(x=split_plot$Value,y=split_plot$Improve/max(split_plot$Improve),pch=toString(j),
               xlim = ylim, ylim = c(0,1), xlab = "Lat", ylab = "Imp/max(Imp)",
               main = paste0(" Split#",i))
        }
        else points(x=split_plot$Value,y=split_plot$Improve/max(split_plot$Improve),pch=toString(j))
      }
      for (j in 1:i) {
        split_plot <- split_raw[which(split_raw$Key=="Lon"&split_raw$Cell==j),]
        if(j==1) {
          plot(x=split_plot$Value,y=split_plot$Improve/max(split_plot$Improve),pch=toString(j),
               xlim = xlim, ylim = c(0,1), xlab = "Lon", ylab = "Imp/max(Imp)",
               main = paste0(" Split#",i))
        }
        else points(x=split_plot$Value,y=split_plot$Improve/max(split_plot$Improve),pch=toString(j))
      }
      dev.off()

      # Compare Improvement against year
      if(year==TRUE) {
        if(pdf==FALSE) png(paste0(save_dir,select_name,"split",i,"(year).png"),width = 800,height = 800)
        else pdf(paste0(save_dir,select_name,"split",i,"(year).pdf"),width = 5,height = 5)

        for (j in 1:i) {

          split_plot <- split_raw[which(split_raw$Key=="Year"&split_raw$Cell==j),]
          if(j==1) {
            plot(x=split_plot$Value,y=split_plot$Improve/max(split_plot$Improve),pch=toString(j),
                 xlim = tlim, xlab = "Year", ylab = "Imp/max(Imp)",
                 main = paste0(" Split#",i))
          }
          else points(x=split_plot$Value,y=split_plot$Improve/max(split_plot$Improve),pch=toString(j))
        }

        dev.off()

      }
    }

  }

  # fix cyclic quarter's output names
  Record <- rename_CQrt(Record)

  # adjust the lat and lon split values
  lat_all <- sort(as.numeric(unique(LF$lat[LF$dummy==FALSE])))
  if(length(unique(lat_all[2:length(lat_all)]-lat_all[1:(length(lat_all)-1)]))==1) lat_adjust <- (lat_all[2] - lat_all[1]) / 2
  else {
    lat_adjust <- 0
    print("Irregular latitudinal grid! No latitudinal adjustment is applied to Record.csv!")
  }
  Record$Value[which(Record$Key=="Lat")] <- as.numeric(Record$Value[which(Record$Key=="Lat")]) + lat_adjust

  lon_all <- sort(as.numeric(unique(LF$lon[LF$dummy==FALSE])))
  if(length(unique(lon_all[2:length(lon_all)]-lon_all[1:(length(lon_all)-1)]))==1) lon_adjust <- (lon_all[2] - lon_all[1]) / 2
  else {
    lon_adjust <- 0
    print("Irregular longitudinal grid! No longitudinal adjustment is applied to Record.csv!")
  }
  Record$Value[which(Record$Key=="Lon")] <- as.numeric(Record$Value[which(Record$Key=="Lon")]) + lon_adjust

  # save the split results
  write.csv(Record,file=paste0(save_dir,select_name,"Record.csv"))

  return(list("LF"=LF, "Record"=Record))
}
