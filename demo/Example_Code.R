# devtools::install_github('HaikunXu/RegressionTree',ref='main')
library(FishFreqTree)

load(file="manual/LF.RData")

# column names are length bins (in cm); must also include four columns in the data frame:
# lat, lon, year, quarter
LF$weight <- 1

head(LF)

fcol <- 5 # the first column with LF_Tree info
lcol <- 17 # the last column with LF_Tree info
bins <- seq(30,150,10)
Nsplit <- 3 # the number of splits (the number of cells - 1)
save_dir <- "demo/"

# plot lf data as maps
# make.lf.map(LF,fcol,lcol,bins,save_dir)
# plot mean length as maps
# make.meanl.map(LF,fcol,lcol,bins,save_dir)

# run the regression tree with the best three splits
# results are saved in the folder 111 under save_dir
LF_Tree <- run_regression_tree(LF,fcol,lcol,bins,Nsplit,save_dir)

head(LF_Tree$LF)
make.split.map(LF_Tree$LF,Nsplit,save_dir)
# The last few columns with names Flag are the cell numbers under each split

# add a dummy data at lat=-20, lon=60
LF$dummy <- FALSE


LF_Tree <- run_regression_tree(LF,fcol,lcol,bins,Nsplit,save_dir,include_dummy = TRUE)


LF_Tree <- run_regression_tree(LF,fcol,lcol,Nsplit,save_dir)

# run the regression tree again with the second best 2th split
# results are saved in the folder 121 under save_dir
LF_Tree <- run_regression_tree(LF,fcol,lcol,Nsplit,save_dir,manual = TRUE, select = c(1,2,1))

# This setting leads to even better improvement (16.4% vs. 16.2% variance explained)
# than the default run

# loop the regression tree for various combinations of splits
loop_dir <- paste0(save_dir,"loop/")
dir.create(loop_dir)
LF_Tree_Loop <- loop_regression_tree(LF,fcol,lcol,bins,Nsplit=4,save_dir=loop_dir,max_select = 3)

head(LF_Tree_Loop$LF)
