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



3. MultiGrid_v_cycle_mpi: MPI parallelization version of MultiGrid 2D. 






