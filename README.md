## NYU HPC HW4.
### Mu-Hua Chien 
### mhc431@nyu.edu

### Generalize 2-D Multigrid Implementation

1. MultiGrid_v_cycle: serial verion of MultiGrid 2D

2. MultiGrid_v_cycle_omp: OpenMP parallelization version of MultiGrid 2D.

Here, we test Multigrid through a v-cycle, even though the result doesn't quite converge. 

For the pre-smoother uses 5 iteration, and post-smoother uses 10x(nlevel - level+ 1) iteration,  

e.g, if nlevel = 10, the 5th level do 60 iteration at post-smoother. 

Scaling test 

Strong scaling tests with 12 level, computation domain = 4096x4096 .

Initial residual = 8.191000e+03. Final residula = 3.204935e+02.  

(For 4096x4094 domain, pure 1000 Jocobi iteration residual = 8.191000e+03 -> 8.150787e+03.)

| N threads  | 1 | 2 | 4 | 8 | 16 |
|---|---|---|---|---|---|
|time | 50.4842 | 27.3417 | 17.1256 | 15.3324 | 15.0685 |

Weak scaling tests with 9 level computation damian = 2560x2560x n per threads. 

Initial residual = 1.023900e+04. Final residual = 6.756190e+03. (N thread = 16)

Initial residual = 2.559000e+03. Final residual = 1.053838e+02. (N thread = 1) 

| N threads  | 1 | 2 | 4 | 8 | 16 |
|---|---|---|---|---|---|
|time | 3.7630 | 4.0533 | 5.1864 | 9.1620 | 18.6332 |

We can observe that about N = 17M data the data reaches the bottle neck of memory access. 

Hence for number of threads great than 8, the result doesn't scale.  

3. MultiGrid_v_cycle_mpi: MPI parallelization version of MultiGrid 2D. 


Weak Scaling 

| N nodes  | 1 | 4 | 16 | 64 | 256 |
|---|---|---|---|---|---|
|time | 2.6645 | 3.6800 | 7.2071 | 7.2439 | - |

Strong Scaling

| N nodes  | 1 | 4 | 16 | 64 | 256 |
|---|---|---|---|---|---|
|time | 86.7919 | 32.5180 | 15.6227 | 3.9531 | - |
