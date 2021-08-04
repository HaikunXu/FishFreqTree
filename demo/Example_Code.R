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

# run the regression tree
LF_Tree <- run_regression_tree(LF,fcol,lcol,Nsplit,save_dir)

head(LF_Tree$LF)
# The last few columns with names Flag are the cell numbers under each split

# run the regression tree
LF_Tree <- run_regression_tree(LF,fcol,lcol,Nsplit,save_dir,lat.min=1,lon.min=2)

# user-specified regression tree
select <- c(1,2,1) # use the best 1st, 2nd, and 4th splits and the second best 3rd split

LF_Tree <- run_regression_tree(LF,fcol,lcol,Nsplit,save_dir,manual = TRUE, select)

# loop the regression tree
Var <- loop_regression_tree(LF,fcol,lcol,Nsplit,save_dir,max_select = 3)

Var

# backup code
# library(tidyverse)
#
# # plot the spatial distribution of each cell
# for (i in 1:Nsplit) {
#   LF_Tree$Flag <- LF_Tree[[paste0("Flag",i)]]
#   f <- ggplot(data=LF_Tree) +
#     geom_point(aes(x=lon,y=lat,color=as.factor(Flag))) +
#     theme_bw()
#
#   ggsave(f,file=paste0(i,"Flag.png"),width=6,height=6)
# }
#
# # plot the mean LF of each cell
# for (i in 1:Nsplit) {
#   LF_Tree$Flag <- LF_Tree[[paste0("Flag",i)]]
#
#   LF_Tree_final <- LF_Tree %>% gather(fcol:lcol,key="Length",value = "LF_Tree2") %>%
#     mutate(Length=as.numeric(Length)) %>%
#     group_by(Flag,Length) %>%
#     summarise(LF_Tree3=mean(LF_Tree2))
#
#   f <- ggplot(data=LF_Tree_final) +
#     geom_line(aes(x=Length,y=LF_Tree3,color=as.factor(Flag))) +
#     theme_bw()
#     # facet_wrap(~Flag)
#
#   ggsave(f,file=paste0(i,"LF.png"),width=6,height=6)
# }
