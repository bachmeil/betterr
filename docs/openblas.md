# Working With OpenBLAS

You can check out the example file titled testblas.d to see how you can
use OpenBLAS for low-level matrix computations. Note a couple of things,
however.

**It's probably not going to speed things up in most cases.** It seems
natural to think that using OpenBLAS will lead to better performance for
matrix multiplication and other such operations. It will *if you're
reusing a matrix/vector*. If you already have all vectors and matrices
allocated, these functions will be much faster. If you need to allocate
a new matrix/vector to hold the result, there's nothing to gain, since
the matrix operations already call into BLAS to do the calculations.
The R functions will, in most cases, allocate a new matrix/vector to hold
the results. You can avoid that by using BLAS.

**Not everything has been demonstrated in that file.** I just showed how
to use OpenBLAS. I generally don't have a need to go this low-level, so
I have not created a full library.
