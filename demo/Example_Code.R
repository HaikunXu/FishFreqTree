# devtools::install_github('HaikunXu/RegressionTree',ref='main')
library(RegressionTree)

load(file="demo/LF.RData")

# column names are length bins (in cm); must also include four columns in the data frame:
# lat, lon, year, quarter

head(LF)

fcol <- 5 # the first column with LF_Tree info
lcol <- 17 # the last column with LF_Tree info
Nsplit <- 3 # the number of splits (the number of cells - 1)
save_dir <- "demo/"

# run the regression tree with the best three splits
# results are saved in the folder 111 under save_dir
LF_Tree <- run_regression_tree(LF,fcol,lcol,Nsplit,save_dir,year = TRUE)

head(LF_Tree$LF)
# The last few columns with names Flag are the cell numbers under each split

# # add a dummy data at lat=-20, lon=60
# LF$dummy <- FALSE
# LF2 <- rbind(LF,c(31,1,-20,60,rep(0,13),TRUE))

LF_Tree <- run_regression_tree(LF,fcol,lcol,Nsplit,save_dir)

# run the regression tree again with the second best 2th split
# results are saved in the folder 121 under save_dir
LF_Tree <- run_regression_tree(LF,fcol,lcol,Nsplit,save_dir,manual = TRUE, select = c(1,2,1))

# This setting leads to even better improvement (16.4% vs. 16.2% variance explained)
# than the default run

# loop the regression tree for various combinations of splits
loop_dir <- paste0(save_dir,"loop/")
dir.create(loop_dir)
LF_Tree_Loop <- loop_regression_tree(LF,fcol,lcol,Nsplit,save_dir=loop_dir,max_select = 2)

head(LF_Tree_Loop$LF)
