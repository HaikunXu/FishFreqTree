# The IATTC regression tree R package for length frequency data

-   This is the GitHub repository for [IATTC](https://www.iattc.org/HomeENG.htm)'s regression tree algorithm on length frequency data.

-   The R codes for regression tree analysis were originally developed by Cleridy Lennart-Cody and then modified by Haikun Xu to make it automatic as a R package.

-   Please contact Haikun ([hkxu\@iattc.org](mailto:hkxu@iattc.org)) for any questions related to the package

-   Reference: Fisheries Research 102(2010):323-326

-   Please use the R version 3.5 for this package; the latest version 4.x may not work for this package

-   How to install the package: devtools::install_github('HaikunXu/RegressionTree',ref='main')

-   User manual can be found at: <https://github.com/HaikunXu/FishFreqTree/blob/main/manual/Manual.pdf>

------------------------------------------------------------------------

## **Note:**

### Data Format

The input data frame should include at least four columns named exactly as "lat", "lon", "year", and "quarter". The input data frame should also include various columns that record **length frequency** information with column names = length bin. **This regression tree package works with length frequency data so please make sure the input values sum to 1 across length bins.** An example of the input data can be found [here](https://github.com/HaikunXu/RegressionTree/blob/main/demo/LF.RData).

### Model description

This package finds the best multi-cell combination for a length frequency data based on the proportion of variance explained. The variables that are current considered in the code include latitude, longitude, quarter/cyclic quarter, and year (can be turned on by using year=TRUE). For those who don't consider quarter as a splitting dimension (e.g., your model has a time step of one year), please still add a column named "quarter" to the input data with values = 1. In the main functions this package provides (run_regression_tree and loop_regression_tree), you can manually turn off the quarter dimension by adding "quarter = FALSE" as a function argument.

### Main functions

run_regression_tree (type ?run_regression_tree on the console for more info): run the regression tree

loop_regression_tree (type ?loop_regression_tree on the console for more info): loop the regression tree

### Code description

For the nth best split, the code first loops over all existing n cells that are defined by the previous n-1 splits, to find the best split (the one that leads to the maximum variance explained) for every cell. Then those best cell-specific splits are compared to find the split that results in the maximum variance explained. This split is the nth best split. This process is iterated until reaching the maximum number of splits specified by the user.

Users should combine the output figures with output tables to understand the best splits in order. Also, the advanced feature (see the example code for more details) in the package allows users to manually specify some or all splits.
