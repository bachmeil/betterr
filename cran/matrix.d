module cran.matrix;

import betterr.r, betterr.rdata;
import std.conv, std.exception, std.stdio;
/* Wraps things in the Matrix package. Allows a lot of matrix functionality
 * for dense and sparse matrices. */

void init() {
  evalRQ("library(Matrix)");
}

extern (C) {
	/*From dense.c line 961: kind A `char` flag, one of `'.'`, `'d'`, `'l'`, and `'n'`, 
  indicating the "kind" of `.geMatrix` desired.  A dot `'.'` means 
  to preserve the "kind" of `from`. 
  * 
  * new: see https://github.com/cran/Matrix/blob/ff84c47d2cb6bca859a841b8d1358ab33ac2a71a/src/dense.c#L946
  * but usually want it to be 0 
  * transpose_if_vector: Set to 1 if you want to drop dimensions like R does
  * 0 otherwise */
  // Coerce geMatrix to R matrix
  // I believe ndense is true if it's numerical and dense
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
	
  /* Code is an R string equal to "ge", "tr", "sy", "tp", or "sp"
   * uplo is 'U' or 'L'
   * diag is 'N' or 'U' */
  Robj R_matrix_as_dense(Robj from, Robj code, Robj uplo, Robj diag);
  Robj lsq_dense_Chol(Robj x, Robj y);
  
  /* Only works if x and y both have one column. */
  Robj lsq_dense_QR(Robj x, Robj y);
  
  /* X is the matrix. tol is the tolerance. R uses default tolerance
   * 0.0000001. */
  Robj lapack_qr(Robj X, Robj tol);

  Robj unpackedMatrix_transpose(Robj from);
	Robj unpackedMatrix_diag_get(Robj obj, Robj nms);
	Robj unpackedMatrix_diag_set(Robj obj, Robj val);
	Robj unpackedMatrix_symmpart(Robj from);
	Robj matrix_symmpart(Robj from);
	Robj unpackedMatrix_skewpart(Robj from);
	Robj matrix_skewpart(Robj from);
  Robj Mmatrix(Robj args);

  // .Internal(matrix(data, nrow, ncol, byrow, dimnames, missing(nrow), missing(ncol)))
  //~ Robj do_matrix(Robj call, Robj op, Robj args, Robj rho);
  //~ void unsafe_matprod(double *x, int nrx, int ncx, double *y, int nry, int ncy, double *z);
}

struct DenseMatrix {
  long rows;
  long cols;
  RData data;
  double * ptr;
  alias data this;

  /* What's different about this versus this(string, long, long) is that
   * it's less efficient but you don't need to know the dimensions at
   * compile time. */ 
  this(string code) {
		data.data = new RDataStorage;
    import std.datetime;
    auto ct = Clock.currTime().toISOString;
    data.data.name = "x" ~ ct[0..15] ~ ct[16..$];
    /* Necessary to do this manually because it forces it to be a symmetric matrix otherwise. */
    data.data.x = evalR(data.name ~ ` <- ` ~ code ~ `;class(` ~ data.name ~ `) <- 'dgeMatrix';` ~ data.name);
    data.data.refcount = 1;

    int * dims = INTEGER(R_do_slot(data.x, RSymbol("Dim")));
    rows = dims[0];
    cols = dims[1];
    ptr = REAL(R_do_slot(data.x, RSymbol("x")));
  }

  this(long r, long c) {
    this("Matrix(as.numeric(NA), nrow=" ~ std.conv.to!string(r) ~ ", ncol=" ~ std.conv.to!string(c) ~ ", sparse=FALSE)");
  }

  this(string code, long r, long c) {
		data.data = new RDataStorage;
    import std.datetime;
    auto ct = Clock.currTime().toISOString;
    data.data.name = "x" ~ ct[0..15] ~ ct[16..$];
    data.data.x = evalR(data.name ~ ` <- new("dgeMatrix", x=` ~ code ~ `, Dim=c(` ~ std.conv.to!string(r) ~ `L, ` ~ std.conv.to!string(c) ~ `L))`);
    data.data.refcount = 1;
    rows = r;
    cols = c;
    ptr = REAL(R_do_slot(data.x, RSymbol("x"))); 
  }

  double opIndex(long r, long c) {
    enforce(r < this.rows, "First index exceeds the number of rows");
    enforce(c < this.cols, "Second index exceeds the number of columns");
    return ptr[c*this.rows+r];
  }

  void opIndexAssign(double v, long r, long c) {
    ptr[c*rows+r] = v;
  }
  
  string toString() {
    printR(data.x);
    return "";
  }
  
  double det() {
    return dgeMatrix_determinant(data.x, RFalse).scalar;
  }
  
  /* These are pure C functions for when you need speed. R is not involved
   * at all. You're on your own for memory management. You need to 
   * manually unprotect, or perhaps set up reference counting,
   * but you're on your own however you deal with it. This isn't as bad
   * as it might seem - the inside of a loop is really the only time the
   * speed will matter.
   * 
   * See the example in testcran.d that shows how to use std.typecons.Unique
   * to manage the memory inside a loop. */
  //~ Robj solve() {
    //~ return Rf_protect(dgeMatrix_solve(data.x));
  //~ }
  
  //~ Robj matmul(DenseMatrix m) {
    //~ return Rf_protect(dgeMatrix_matrix_mm(data.x, m.data.x, RFalse));
  //~ }
  
  //~ Robj transpose() {
    //~ return Rf_protect(unpackedMatrix_transpose(data.x));
  //~ }
}

