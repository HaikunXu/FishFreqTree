# The IATTC regression tree package for length frequency data

-   This is the GitHub repository for IATTC's (<https://www.iattc.org/HomeENG.htm>) regression tree algorithm on length frequency data.

-   The R codes for regression tree analysis were originally developed by Cleridy Lennart-Cody and then modified by Haikun Xu to make it automatic as a R package.

-   Please contact Haikun ([hkxu\@iattc.org](mailto:hkxu@iattc.org)) for any questions related to the package

-   Reference: Fisheries Research 102(2010):323-326

-   A example code is provided to show how to use the package: <https://github.com/HaikunXu/RegressionTree/blob/main/demo/Example_Code.R>

# 

**Note: For the nth best split, the code first loops over all existing n cells that are defined by the previous n-1 splits, to find the best split (the one that leads to the maximum variance explained) for every cell. Then those best cell-specific splits are compared to find the split that results in the maximum variance explained. This split is the nth best split. This process is iterated until reaching the maximum number of splits specified by the user.**

**Users can combine the output figures with output tables to understand the best splits in order. Also, the advanced feature (see the example code for more details) in the package allows users to manually specify some or all splits for more flexible splits.**
