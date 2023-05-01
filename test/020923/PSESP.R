library(FishFreqTree)

save_dir<- "test/020923/"
fcol <- 6 # the first column with LF_Tree info
lcol <- 19 # the last column with LF_Tree info
Nsplit <- 3 # the number of splits (the number of cells - 1)
bins<-c(20, 30,  40,  50,  60,  70,  80,  90, 100, 110, 120, 130, 140, 150)
LF_Tree <- run_regression_tree(LF=LF,fcol=fcol,
                               lcol=lcol,bins=bins,Nsplit=Nsplit,save_dir=save_dir,year=FALSE,
                               manual = TRUE, select=c(2,1,1))
