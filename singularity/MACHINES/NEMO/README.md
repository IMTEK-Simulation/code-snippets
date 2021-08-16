## AMD nodes: 

These containers might not be compatible with the AMD nodes. I experienced MPI I/O Problems and segmentation faults (even in serial) that didn't arrise on the normal containers.

To avoid AMD nodes even in the express queue, add the option `-l feature=intel`

