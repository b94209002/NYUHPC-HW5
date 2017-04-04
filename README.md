## NYU HPC HW4.
### Mu-Hua Chien 
### mhc431@nyu.edu

mpi_bug1: 

1. For the rank 2, the recv before the send function, so the rank cannot send the information before it receive. Hence it hangs.   

2. The pair of send and recv function do not have the same tag, so the cannot receive the information.   


mpi_bug2:
 
The alpha and beta does not share the same data type, so the datatype in mpi_irecv, MPI_float, is modified. 

mpi_bug3: 

This file does not include mpi_initial and mpi_finalize. 

mpi_bug4: 

The master rank does not see the mpi_reduce. Add the end line before the master rank print the answer. 

