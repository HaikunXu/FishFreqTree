#' Run the regression tree analysis with specified number of splits
#'
#' \code{run_regression_tree} This function generates an area code for longline catch allocation
#'
#' @param LF The length frequency data frame input; must include columns lat, lon, year, and quarter
#' @param fcol The first column in the data frame with length frequency info
#' @param lcol The first column in the data frame with length frequency info
#' @param Nsplit The number of splits
#' @param save_dir The directory where results will be saved
#' @param manual Whether to use user-specified splits; default = FALSE
#' @param select User-specified splits; default = NA
#' @param quarter Whether to consider quarter as a splitting dimention; default = TRUE
#'
#'
#' @export

run_regression_tree <- function(LF,fcol,lcol,Nsplit,save_dir,manual=FALSE,select=NA,lat.min=1,lon.min=1,quarter=TRUE) {

if(manual==FALSE) select <- rep(1,Nsplit) # every split choose the one with the max improvement

var_exp <- rep(NA,Nsplit) # variance explained
folder_name <- paste0(gsub(", ","",toString(select)))
select_name <- paste0(gsub(", ","",toString(select)),"/")
unlink(paste0(save_dir,folder_name), recursive = TRUE)
dir.create(paste0(save_dir,select_name))

for (i in 1:Nsplit) {

  png(paste0(save_dir,select_name,"split",i,".png"),width = 1000,height = 500)
  par(mfrow=c(1,2))

  if(i==1) {
    # whole ALB area
    split <- find_split(LF,fcol,lcol,lat.min,lon.min,quarter)
    # save result as a csv.file
    write.csv(split,file=paste0(save_dir,select_name,"split",i,".csv"),row.names = FALSE)

    # if((manual==FALSE)|(i %in% user_split$Number == FALSE))
      LF <- make.Us.areaflags.f(LF, as.character(split$Key[select[i]]), as.numeric(split$Value[select[i]]),1,0)
    # else
      # LF <- make.Us.areaflags.f(LF, user_split$Key[which(user_split$Number==i)], user_split$Value[which(user_split$Number==i)],1,0)

    e0 <- get.klderror.null(as.matrix(LF[,fcol:lcol])) # null (no stratification)
    e1 <- get.klderror.null(as.matrix(LF[LF$Flag1==1,fcol:lcol]))
    e2 <- get.klderror.null(as.matrix(LF[LF$Flag1==2,fcol:lcol]))

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
    xlim <- c(min(LF$lon),max(LF$lon))
    ylim <- c(min(LF$lat),max(LF$lat))

    # plot the spatial distribution of cells
    for (j in 1:(i+1)) {
      LF_plot <- LF[which(Flag==j),]
      if(j==1) {
        plot(x=LF_plot$lon,y=LF_plot$lat,pch=toString(j),
             xlim = xlim, ylim = ylim, xlab = "Lon", ylab = "Lat",
             main = paste0("Split#",i,": ",split[select[i],2],"<=",split[select[i],3]))
      }
      else {
        points(x=LF_plot$lon,y=LF_plot$lat,pch=toString(j))
      }
    }

    # Compare the LF among cells
    for (j in 1:(i+1)) {
      LF_plot <- LF[which(Flag==j),]
      LF_mean <- apply(LF_plot,2,mean)
      Length <- as.numeric(names(LF_mean)[fcol:lcol])
      LF_mean <- LF_mean[fcol:lcol]

      if(j==1) {
        plot(x=Length,y=LF_mean,pch=toString(j),
             ylim = c(0,0.3),
             main = paste0(round((e0-e1-e2)/e0*100,2),"% variance explained"))
        lines(x=Length,y=LF_mean)
      }
      else {
        points(x=Length,y=LF_mean,pch=toString(j))
        lines(x=Length,y=LF_mean)
      }
    }

    dev.off()
  }

  else {
    imp <- rep(NA,i)
    # for loop for every possible split
    for (ii in 1:i) {
      LF_raw <- LF
      j <- which(LF[[paste0("Flag",i-1)]] == ii)
      LF_data <- LF[j,]
      split <- find_split(LF_data,fcol,lcol,lat.min,lon.min,quarter)
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

      imp[ii] <- (e0-e)/e0

    }

    # the one with the most improvement
    ii <- which(imp==max(imp))

    # j <- which(LF[[paste0("Flag",i-1)]] == ii)
    # LF_data <- LF[j,]
    # split <- find_split(LF_data,fcol,lcol,lat.min,lon.min)

    # save result as a csv.file
    split_raw <- split_raw[order(split_raw$Improve,decreasing = TRUE),]
    write.csv(split_raw,file=paste0(save_dir,select_name,"split",i,".csv"),row.names = FALSE)

    # if((manual==FALSE)|(i %in% user_split$Number == FALSE))
      LF <- make.Us.areaflags.f(LF, as.character(split_raw$Key[select[i]]), as.numeric(split_raw$Value[select[i]]),i,ii)
    # else
      # LF <- make.Us.areaflags.f(LF, user_split$Key[which(user_split$Number==i)], user_split$Value[which(user_split$Number==i)],i,ii)

    for (k in 1:(i+1)) {
      if(k==1) e <- get.klderror.null(as.matrix(LF[LF[[paste0("Flag",i)]]==1,fcol:lcol]))
      else e <- e + get.klderror.null(as.matrix(LF[LF[[paste0("Flag",i)]]==k,fcol:lcol]))
    }

    var_exp[i] <- (e0-e)/e0

    # print result to screen
    if(select[i]==1) {
    if(i==2)
      print(paste0("Best 2nd split is for cell ",ii," in split",i-1,".png: ",split_raw[1,2],"<=",split_raw[1,3]))
    if(i==3)
      print(paste0("Best 3rd split is for cell ",ii," in split",i-1,".png: ",split_raw[1,2],"<=",split_raw[1,3]))
    if(i>3)
      print(paste0("Best ",i,"th"," split is for cell ",ii," in split",i-1,".png: ",split_raw[1,2],"<=",split_raw[1,3]))
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

      # png(paste0(save_dir,select_name,"split",i,".png"),width = 1000,height = 500)
      # par(mfrow=c(1,2))

      # plot the spatial distribution of cells
      for (j in 1:(i+1)) {
        LF_plot <- LF[which(Flag==j),]
        if(j==1) {
          plot(x=LF_plot$lon,y=LF_plot$lat,pch=toString(j),
               xlim = xlim, ylim = ylim, xlab = "Lon", ylab = "Lat",
               main = paste0("Split#",i,": ",split_raw[select[i],2],"<=",split_raw[select[i],3]))
        }
        else {
          points(x=LF_plot$lon,y=LF_plot$lat,pch=toString(j))
        }
      }

      # Compare the LF among cells
      for (j in 1:(i+1)) {
        LF_plot <- LF[which(Flag==j),]
        LF_mean <- apply(LF_plot,2,mean)
        Length <- as.numeric(names(LF_mean)[fcol:lcol])
        LF_mean <- LF_mean[fcol:lcol]

        if(j==1) {
          plot(x=Length,y=LF_mean,pch=toString(j),
               ylim = c(0,0.3),
               main = paste0(round((e0-e)/e0*100,2),"% variance explained"))
          lines(x=Length,y=LF_mean)
        }
        else {
          points(x=Length,y=LF_mean,pch=toString(j))
          lines(x=Length,y=LF_mean)
        }
      }

      dev.off()
    }
}

return(list("LF"=LF,"Var"=var_exp))
}
