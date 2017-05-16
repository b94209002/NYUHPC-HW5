## NYU HPC HW4.
### Mu-Hua Chien 
### mhc431@nyu.edu

### Generalize 2-D Multigrid Implementation

1. MultiGrid_v_cycle: serial verion of MultiGrid 2D

2. MultiGrid_v_cycle_omp: OpenMP parallelization version of MultiGrid 2D.

Scaling test 

Strong scaling tests with 12 level, computation domain = 4096x4096 .

| N threads  | 1 | 2 | 4 | 8 | 16 |
|---|---|---|---|---|---|
|time | 50.3665 | 27.3401 | 17.1039 | 15.2863 | 11.3946 |

Weak scaling tests with 9 level computation damian = 2560x2560x n threads. 

| N threads  | 1 | 2 | 4 | 8 | 16 |
|---|---|---|---|---|---|
|time | 3.7501 | 4.0582 | 5.1845 | 9.1484 | 18.4537 |

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
