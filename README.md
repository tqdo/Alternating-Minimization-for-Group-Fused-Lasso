# Group Fused Lasso in C++
This repository contains two efficient algorithms to solve the Group Fused Lasso (Bleakley & Vert, 2011) problem. 
The first one solves the primal problem using Alternating Minimization algorithm.
The second one solves the dual problem using Block coordinate descent.

Both algorithms are coded in C++ using Armadillo library. The download link for Armadillo is http://arma.sourceforge.net/download.html.

Users only need to specify the text file, which contains the data, and how many lambdas along the path they want to 
get the solution.
The data should be n*p, where n is number of observations and p is the dimension of the signal.
