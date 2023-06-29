# Calling the Matrix Package

The Matrix package is installed by default with R. The main reason to be
interested in it is because it exposes low-level (i.e., C) functions
for linear algebra. You should read the testcran.d file for examples and
a way of using it with a unique pointer.

You need to link to Matrix.so. You can find the location by running this
command in R:

```
find.package("Matrix")
```

and then appending /libs/Matrix.so to the end of the output. On my Ubuntu,
I'm linking to -L/usr/lib/R/library/Matrix/libs/Matrix.so.

The functionality currently exposed is:

```
Robj R_geMatrix_as_matrix(Robj from, Robj ndense);
Robj denseLU_determinant(Robj obj, Robj logarithm);
Robj dgeMatrix_trf_(Robj obj, int warn);
Robj dgeMatrix_trf(Robj obj, Robj warn);
Robj dgeMatrix_norm(Robj obj, Robj type);
Robj dgeMatrix_rcond(Robj obj, Robj type);
Robj dgeMatrix_determinant(Robj obj, Robj logarithm);
Robj dgeMatrix_solve(Robj a);
Robj dgeMatrix_matrix_solve(Robj a, Robj b);
Robj dgeMatrix_crossprod(Robj x, Robj trans);
Robj geMatrix_crossprod(Robj x, Robj trans);
Robj dgeMatrix_dgeMatrix_crossprod(Robj x, Robj y, Robj trans);
Robj geMatrix_dgeMatrix_crossprod(Robj x, Robj y, Robj trans);
Robj dgeMatrix_matrix_crossprod(Robj x, Robj y, Robj trans);
Robj geMatrix_matrix_crossprod(Robj x, Robj y, Robj trans);
Robj dgeMatrix_matrix_mm(Robj a, Robj b, Robj right);
Robj geMatrix_matrix_mm(Robj a, Robj b, Robj right);
Robj dgeMatrix_Schur(Robj x, Robj vectors, Robj isDGE);
Robj dgeMatrix_svd(Robj x, Robj nu, Robj nv);
Robj dgeMatrix_exp(Robj x);
Robj R_dense_colSums(Robj obj, Robj narm, Robj mean);
Robj R_dense_rowSums(Robj obj, Robj narm, Robj mean);
Robj R_matrix_as_dense(Robj from, Robj code, Robj uplo, Robj diag);
Robj lsq_dense_Chol(Robj x, Robj y);
Robj lsq_dense_QR(Robj x, Robj y);
Robj lapack_qr(Robj X, Robj tol);
Robj unpackedMatrix_transpose(Robj from);
Robj unpackedMatrix_diag_get(Robj obj, Robj nms);
Robj unpackedMatrix_diag_set(Robj obj, Robj val);
Robj unpackedMatrix_symmpart(Robj from);
Robj matrix_symmpart(Robj from);
Robj unpackedMatrix_skewpart(Robj from);
Robj matrix_skewpart(Robj from);
```

This should do the bulk of the linear algebra you'd want to do in C. At
this time, I don't have a lot of documentation; you may need to figure
some of it out for yourself.

I may improve all of this, including the documentation, at some point. It's
something of a distraction right now, since it's not exactly the purpose
of this project. More likely, I'd fork the Matrix package and add functionality
that makes it easy to work with the usual matrix type in R rather than the
special Matrix type (which is similar, but distinct, and requires conversion).
