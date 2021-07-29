#' Run the regression tree analysis with specified number of splits
#'
#' \code{run_regression_tree} This function generates an area code for longline catch allocation
#'
#' @param LF The length frequency data frame input; must include columns lat, lon, year, and quarter
#' @param fcol The first column in the data frame with length frequency info
#' @param lcol The first column in the data frame with length frequency info
#' @param Nsplit The number of splits
#' @param save_dir The directory where figures will be saved
#'
#' @export


run_regression_tree <- function(LF,fcol,lcol,Nsplit,save_dir) {

var_exp <- rep(NA,Nsplit) # variance explained

for (i in 1:Nsplit) {

  png(paste0(save_dir,"split",i,".png"),width = 1000,height = 500)
  par(mfrow=c(1,2))

  if(i==1) {
    # whole ALB area
    split <- find_split(LF,fcol,lcol)

    LF <- make.Us.areaflags.f(LF, as.character(split$Key[1]), as.numeric(split$Value[1]),1,0)

    e0 <- get.klderror.null(as.matrix(LF[,fcol:lcol])) # null (no stratification)
    e1 <- get.klderror.null(as.matrix(LF[LF$Flag1==1,fcol:lcol]))
    e2 <- get.klderror.null(as.matrix(LF[LF$Flag1==2,fcol:lcol]))

    var_exp[1] <- (e0-e1-e2)/e0

    # print to the screen
    cat("\n")
    cat("Note: below shows the best splits in order, please check the save figures under save_dir to better understand the meaning of each split\n\n")
    print(paste0("Best 1st split: ",split[1,2],"<=",split[1,3]))
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
             main = paste0("Split#",i,": ",split[1,2],"<=",split[1,3]))
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
      split <- find_split(LF_data,fcol,lcol)

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

    j <- which(LF[[paste0("Flag",i-1)]] == ii)
    LF_data <- LF[j,]
    split <- find_split(LF_data,fcol,lcol)

    LF <- make.Us.areaflags.f(LF,
                                  as.character(split$Key[1]),
                              as.numeric(split$Value[1]),i,ii)

    for (k in 1:(i+1)) {
      if(k==1) e <- get.klderror.null(as.matrix(LF[LF[[paste0("Flag",i)]]==1,fcol:lcol]))
      else e <- e + get.klderror.null(as.matrix(LF[LF[[paste0("Flag",i)]]==k,fcol:lcol]))
    }

    var_exp[i] <- (e0-e)/e0

    # print result to screen
    if(i==2) print(paste0("Best 2nd split: ",split[1,2],"<=",split[1,3]))
    if(i==3) print(paste0("Best 3rd split: ",split[1,2],"<=",split[1,3]))
    if(i>3) print(paste0("Best ",i,"th"," split: ",split[1,2],"<=",split[1,3]))

      cat(paste0((e0-e)/e0*100,"% variance explained\n\n"))

      # plot result
      Flag <- LF[[paste0("Flag",i)]]

      png(paste0(save_dir,"split",i,".png"),width = 1000,height = 500)
      par(mfrow=c(1,2))

      # plot the spatial distribution of cells
      for (j in 1:(i+1)) {
        LF_plot <- LF[which(Flag==j),]
        if(j==1) {
          plot(x=LF_plot$lon,y=LF_plot$lat,pch=toString(j),
               xlim = xlim, ylim = ylim, xlab = "Lon", ylab = "Lat",
               main = paste0("Split#",i,": ",split[1,2],"<=",split[1,3]))
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

return(LF)
}
