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
LF_Tree <- run_regression_tree(LF,fcol,lcol,Nsplit,save_dir)

head(LF_Tree$LF)
# The last few columns with names Flag are the cell numbers under each split

# run the regression tree again with the second best 2th split
# results are saved in the folder 121 under save_dir
LF_Tree <- run_regression_tree(LF,fcol,lcol,Nsplit,save_dir,manual = TRUE, select = c(1,2,1))

# This setting leads to even better improvement (16.4% vs. 16.2% variance explained)
# than the default run

# loop the regression tree for various combinations of splits
Var <- loop_regression_tree(LF,fcol,lcol,Nsplit,save_dir,max_select = 3)

# column Var is the % variance explained by each split combination
Var
# the best splits is 1, 2, 1

# use the best combination of splits found above through the loop
select <- as.numeric(Var[1,1:Nsplit]) # the first row has the highest % variance explained

LF_Tree <- run_regression_tree(LF,fcol,lcol,Nsplit,save_dir,manual = TRUE, select)
head(LF_Tree$LF)
