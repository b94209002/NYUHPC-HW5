## NYU HPC HW4.
### Mu-Hua Chien 
### mhc431@nyu.edu

mpi_bug1: 

1. For the rank 2, the recv before the send function, so the rank cannot send the information before it receive. Hence it hangs.   

2. The pair of send and recv function do not have the same tag, so the cannot receive the information.   


mpi_bug2:
 

