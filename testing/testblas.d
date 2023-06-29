import cblas.cblas;
import betterr.matrix, betterr.r, betterr.vector;
import std.conv, std.stdio;

void main() {
  startR();
  auto m1 = Matrix(2,2);
  Column(m1, 0) = [1.1, 2.2];
  Column(m1, 1) = [3.3, 4.4];
  
  auto m2 = Matrix(2,2);
  Column(m2, 0) = [0.1, 0.2];
  Column(m2, 1) = [0.3, 0.4];
  
  auto m3 = Matrix(2,2);
  
  gemm(Order.ColMajor, Transpose.NoTrans, Transpose.NoTrans,
    2, 2, 2, 1.0, m1.ptr, 2, m2.ptr, 2, 0.0, m3.ptr, 2);
  m3.print("Matrix product of m1 and m2 calculated by BLAS");
  
  evalRQ(`m1 <- matrix(c(1.1, 2.2, 3.3, 4.4), nrow=2)`);
  evalRQ(`m2 <- matrix(c(0.1, 0.2, 0.3, 0.4), nrow=2)`);
  evalRQ(`print("Matrix product calculated by R:");print(m1 %*% m2)`);
  
  // Example using the low-level interface
  // Does no checks, does not handle NA values
  auto m4 = Matrix(3,2);
  Column(m4, 0) = [1.1, 2.2, 3.3];
  Column(m4, 1) = [4.4, 5.5, 6.6];
  
  auto m5 = Matrix(2,3);
  Row(m5, 0) = [0.1, 0.2, 0.3];
  Row(m5, 1) = [0.4, 0.5, 0.6];
  
  auto m6 = Matrix(3,3);
  matmul_fast(m4, m5, m6);
  m6.print("Matrix product of m4 and m5");
  
  auto v1 = Vector([2.4, 4.8, 7.2]);
  auto v2 = Vector([2.1, 4.5, 6.9]);
  swap(v1.length.to!int, v1.ptr, 1, v2.ptr, 1);
  v1.print("v1 after the swap");
  v2.print("v2 after the swap");
  
  /* Using the fast function
   * Swap the elements back */
  swap_fast(v1, v2);
  v1.print("v1 after the swap");
  v2.print("v2 after the swap");
  
  /* Swap works with Matrix types too */
  m1.print("m1 before the swap");
  m2.print("m2 before the swap");
  swap_fast(m1, m2);
  m1.print("m1 after the swap");
  m2.print("m2 after the swap");
  
  /* Scaling */
  v1.print("v1 before scaling");
  scal_fast(v1, 1.2);
  v1.print("v1 after scaling");

  m1.print("m1 before scaling");
  scal_fast(m1, 1.2);
  m1.print("m1 after scaling");
  
  /* Copying */
  auto v3 = Vector(v1.length);
  copy_fast(v1, v3);
  v1.print("v1");
  v3.print("Copy of v1");
  
  auto m7 = Matrix(m1.rows, m2.rows);
  copy_fast(m1, m7);
  m1.print("m1");
  m7.print("Copy of m1");
  
  /* ax + y */
  axpy_fast(m7, m7, 1.5);
  m7.print("m7");
  m7.print("1.5*m7 + m7");
  
  /* Dot product */
  v1.print("v1");
  v2.print("v2");
  writeln("Dot product of v1 and v2: ", dot_fast(v1, v2));
  
  /* Norm of vector */
  writeln("Norm of v1: ", dnrm2_fast(v1));
  
  /* Sum of absolute values */
  writeln("Sum of abs([1.1, -1.1, -2.2]): ", asum_fast([1.1, -1.1, -2.2]));
  
  /* Index of largest absolute value */
  writeln("Index of largest of [1.1, -3.6, 2.4, 3.1, -2.8]: ", idamax_fast([1.1, -3.6, 2.4, 3.1, -2.8]));

	/* Matrix-vector multiply */
	m1.print("m1");
	auto d1 = Vector([0.2, 0.4]);
	auto res = Vector(m1.rows);
	writeln("res: ", res);
	dgemv_fast(m1, d1, res);
	res.print("Product of m1 and d1");
	
  closeR();
}

/* These are examples of functions you might want to write. You might
 * want to customize them for your use case. */
void matmul_fast(Matrix x, Matrix y, Matrix result) {
  cblas_dgemm(Order.ColMajor, Transpose.NoTrans, Transpose.NoTrans,
    x.rows.to!int, y.cols.to!int, x.cols.to!int, 1.0, x.ptr, x.rows.to!int,
    y.ptr, y.rows.to!int, 0.0, result.ptr, result.rows.to!int);
}

/* Works with either Matrix or Vector
 * Matrix needs to be square or the dimensions won't make sense */
void swap_fast(T)(T x, T y) {
  cblas_dswap(x.length.to!int, x.ptr, 1, y.ptr, 1);
}

void scal_fast(T)(T x, double a) {
  cblas_dscal(x.length.to!int, a, x.ptr, 1);
}

void copy_fast(T)(T x, T y) {
  cblas_dcopy(x.length.to!int, x.ptr, 1, y.ptr, 1);
}

void axpy_fast(T)(T x, T y, double a) {
  cblas_daxpy(x.length.to!int, a, x.ptr, 1, y.ptr, 1);
}

double dot_fast(Vector x, Vector y) {
  return cblas_ddot(x.length.to!int, x.ptr, 1, y.ptr, 1);
}

double dnrm2_fast(Vector x) {
  return cblas_dnrm2(x.length.to!int, x.ptr, 1);
}

double asum_fast(T)(T x) {
	return cblas_dasum(x.length.to!int, x.ptr, 1);
}

ulong idamax_fast(T)(T x) {
	return cblas_idamax(x.length.to!int, x.ptr, 1);
}

void dgemv_fast(T)(Matrix m, T v, T result, double alpha=1.0, double beta=0.0) {
	cblas_dgemv(Order.ColMajor, Transpose.NoTrans, m.rows.to!int, m.cols.to!int,
		alpha, m.ptr, m.rows.to!int, v.ptr, 1, beta, result.ptr, 1);
}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
