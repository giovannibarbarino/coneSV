# Computing cone-constrained singular values of matrices

The methods in the repository solve the problem

$$\min u^TAv \quad : \quad u \in P, \quad v \in Q,\quad  \\\|u\\\|=\\\|v\\\|=1$$   

where

$$P = \\\{u | u = Gx,\\, x\leq 0\\\}, \qquad Q = \\\{u | u = Hx,\\, x \geq 0\\\}$$


## Running the Algorithms

The main algorithms and their respective codes are in the folder "algorithms" with names

 -  Sequential Regularized Partial Linearization (SRPL) :  AlternatingOptimization.m
 -  Extrapolated Alternating Optimization (E-AO) :  srpl_poly.m
 -  Brute-Force Active-Set method (BFAS) :  bfas_timestamps_test.m
 -  Gurobi (Gur) :  generators_gurobi_uv_test.m
 
 To produce the results of Table 7, the following additional algorithm is used
 
 -  Sequential Regularized Partial Linearization (SRPL) for cones of matrices : ToRun_SDP_N_stat.m
 
 Full documentation on how to use the codes is found inside the respective files.
 
## Running the Experiments

Reference: Giovanni Barbarino, Nicolas Gillis, David Sossa, <br>
''Computing cone-constrained singular values of matrices'', 2025. 
 
 The results in the Figure and Tables of the article are generated by the respective codes found in the folder "Experiments". <br>
 The external data used to test the algorithms and produce Table 3 can be found in the folder "Benchmarks".
 
 
 
 

 
 
