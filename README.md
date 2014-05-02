rand_approx
===========

Randomized svd calculator

To compile:

edit makefile to point to your compiler and BLAS and LAPACK libraries

./make

The main executable is in bin/ and its usage is given below

./main matrix-file matrix-transpose-file p q output-file numthreads

Notes:
1. The sparse matrix files must be in Matrix Market[1] format. 
2. p is the number of samples and will be the number of approximate singular
   values computed.
3. q is the number of power iterations used.
4. The singular values are written to the output-file in csv format. 
5. numthreads is the number of threads to be used for the parallel matvec.

[1] http://math.nist.gov/MatrixMarket/ 
