module cblas.interface;
import betterr.matrix, betterr.r, betterr.vector;
import std.conv;

/* This module contains functions that provide a convenient interface to
 * BLAS and LAPACK functionality as you'd find in OpenBLAS. If you
 * include this module in your program, you need to link to OpenBlas or
 * another BLAS. */

// Efficient copy
void copy_fast(T)(double * source, T n, double * target) {
	cblas_dcopy(n.to!int, source, 1, target, 1);
}

void copy_fast(double[] source, double[] target) {
	enforce(source.length == target.length, "Dimensions do not match");
	copy_fast(source.ptr, source.length, target.ptr);
}

void copy_fast(Matrix source, Matrix target) {
	enforce(source.rows == target.rows, "Rows do not match");
	enforce(source.cols == target.cols, "Columns do not match");
	copy_fast(source.ptr, source.rows*source.cols, target.ptr);
}

void copy_fast(Vector source, Vector target) {
	enforce(source.rows == target.rows, "Rows do not match");
	copy_fast(source.ptr, source.rows, target.ptr);
}

// Efficient swapping
void swap_fast(T)(double * x, double * y, T n) {
	cblas_dswap(n.to!int, x, 1, y, 1);
}

void swap_fast(double[] x, double[] y) {
	enforce(x.length == y.length, "Dimensions do not match");
	cblas_dswap(x.ptr, y.ptr, x.length);
}

void swap_fast(Matrix x, Matrix y) {
	enforce(x.rows == y.rows, "Rows do not match");
	enforce(x.cols == y.cols, "Columns do not match");
	cblas_dswap(x.ptr, y.ptr, x.rows*x.cols);
}

void swap_fast(Vector x, Vector y) {
	enforce(x.length == y.length, "Dimensions do not match");
	cblas_dswap(x.ptr, y.ptr, x.length);
}

// Efficient scaling
// Note: x is overwritten, so send a copy if you need the original
void scale_fast(T)(double * x, T n, double alpha) {
	cblas_dscal(n.to!int, alpha, x, 1);
}

void scale_fast(double[] x, double alpha) {
	cblas_dscal(x.ptr, x.length, alpha);
}

void scale_fast(Matrix x, double alpha) {
	cblas_dscal(x.ptr, x.rows*x.cols, alpha);
}

void scale_fast(Vector x, double alpha) {
	cblas_dscal(x.ptr, x.length, alpha);
}

// Summation: Note that BLAS does not provide a sum of vector elements.

// Efficient dot product (sum of vector elements)
// For calculating predictions from a regression model
double dot_fast(T)(double * x, double * y, T n) {
	return cblas_ddot(n.to!int, x, 1, y, 1);
}

double dot_fast(Vector x, Vector y) {
	enforce(x.rows == y.rows, "Dimensions do not match");
	return dot_fast(x.ptr, y.ptr, x.rows);
}

double dot_fast(double[] x, double[] y) {
	enforce(x.length == y.length, "Dimensions do not match");
	return dot_fast(x.ptr, y.ptr, x.length);
}

// Matrix multiplication
void matmul_fast(Matrix x, Matrix y, Matrix result) {
	cblas_dgemm(Order.ColMajor, Transpose.NoTrans, Transpose.NoTrans,
		x.rows.to!int, y.cols.to!int, x.cols.to!int, 1.0, x.ptr, x.rows.to!int,
		y.ptr, y.rows.to!int, 0.0, result.ptr, result.rows.to!int);
}
