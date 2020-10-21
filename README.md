# spatial-stats-optimisation
Code associated with a paper by Lee and Meeks on a graph optimisation method for analysis of spatial data

The main function is iterative_opt, which takes the following two arguments:
1.	A base binary K * K neighbourhood (adjacency) matrix based on a purely geographic rule such as border sharing.
2.	A K * 1 vector containing the residual (after covariate adjustment) spatial structure in the data on the linear predictor scale.

