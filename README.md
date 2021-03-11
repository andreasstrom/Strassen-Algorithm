# A Strassen Algroithm for matrix multiplication
An implementation of the Strassen algorithm in C. This implementation stores the matrices as 1D arrays, 
utilizes a stride-1 naive algorithm for multiplying "small enough" matrices, and is parallelized for 1, 7 or 49 threads using POSIX threads.
