#######################################
####                               ####
####            MOIF v0.1          ####
####                               ####
#######################################


1) Description of the library
This library implements the MultiObjective Implicit Filtering (MOIF) algorithm proposed in [1] for multiobjective optimization problems with box-constraints. 
                                                        
                                                
   min f_1(x), f_2(x), ..., f_m(x)
   s.t. l <= x <= u

The performed optimization can return either an approximation of the entire Pareto front or a single approximate Pareto point.

Reference 
[1]: G. Cocchi, G. Liuzzi, A. Papini, M. Sciandrone.
An implicit filtering algorithm for derivative-free multiobjective optimization with box constraints. Computational Optimization and Applications, doi:10.1007/s10589-017-9953-2 (2017)

2) Usage of the library
"moif" is the core function of the MOIF library. The user must provide the matlab function func_F which requires one input argument (the vector of variables x) and returns the m-dimentional vector of objective function values

F(x) = [f_1(x); f_2(x); ... ; f_m(x)]

3) Parameters of the library
For setting up the algorithm parameters, please modify the parameters_moif.m file. It contains all the implicit filtering parameters, line search parameters and stopping rules.

4) Example usage of the library
An example usage of the library on a test problem is provided in the file example.m