import betterr.matrix, betterr.r;
import cran.matrix;
import std.conv, std.exception, std.stdio, std.typecons, std.utf;

void main() {
	startR();
	cran.matrix.init();
	
	auto m1 = evalR("tmp <- Matrix(c(1.1, 2.2, 3.3, 4.4), ncol=2);tmp");
	write("m1: ");
	printR(m1);
	writeln();
	
	auto m2= dgeMatrix_matrix_mm(m1, m1, RFalse);
	write("m2: ");
	printR(m2);
	writeln();

	auto m3= geMatrix_matrix_mm(m1, m1, RFalse);
	write("m3: ");
	printR(m3);
	writeln();

	/* These operations demonstrate the use of the functions without any
	 * memory management technique. Since the functions protect the Robj,
	 * that memory will never be freed. */
  write("Determinant of m1: ");
  printR(dgeMatrix_determinant(m1, RFalse));
  writeln();
  
  write("Inverse of m1: ");
  printR(dgeMatrix_solve(m1));
  writeln();
  
  /* Low-level access to the memory. Normally indexing will work better,
   * but you can do it this way if you want complete control. */
  double * ptr = REAL(R_do_slot(m1, RSymbol("x")));
  writeln("The individual array elements of m1:");
  writeln(ptr[0]);
  writeln(ptr[1]);
  writeln(ptr[2]);
  writeln(ptr[3]);
  writeln();
  
  ptr[2] = -1.6;
	write("m1 after changing the [0,1] element: ");
  printR(m1);
  writeln();
  
  //~ printR(Rf_ScalarReal(1.889));
  //~ printR(RTrue);
  //~ printR(RFalse);
  
  /* Can convert to a regular R matrix if you want to do something that
   * can't be done with the Matrix package. */
  writeln("Converting m1 to a regular R matrix:");
  printR(R_geMatrix_as_matrix(m1, RTrue));
  writeln();
  
  writeln("Did that convert m1 to a regular R matrix? ", isMatrix(R_geMatrix_as_matrix(m1, RTrue)));
  writeln("Is m2 a regular R matrix? ", isMatrix(m2));
  writeln();
  
  write("Transpose of m1: ");
  printR(unpackedMatrix_transpose(m1));
  writeln();
  
  auto m4 = DenseMatrix(4, 4);
  write("Value of m4 on creation: ");
  writeln(m4);
  writeln();
  
  /* Use indexing to change the elements */
  foreach(col; 0..4) {
    foreach(row; 0..4) {
      m4[row, col] = col*0.25+row;
    }
  }
  write("m4 after updating the values: ");
  writeln(m4);
  writeln();
  
  /* As noted above, those functions above require you to do the
   * unprotection of the Robj manually. In some of the calls, since the
   * Robj was not saved, there is no way to ever free that memory
   * without killing the R session using closeR() and starting over.
   * 
   * One alternative, that is demonstrated here, and should work most of
   * the time, is to use std.typecons.Unique. As you can verify, both
   * m5 and m6 are freed. In practice, you would not want the writeln
   * statement inside ~this. */
  writeln("Let's demonstrate the use of std.typecons.Unique.");
  void foo() { 
		auto mat1 = DenseMatrix(2, 2);
		mat1[0,0] = 1.0;
		mat1[1,0] = 7.5;
		mat1[0,1] = 2.4;
		mat1[1,1] = 4.0;
		writeln("mat1: ", mat1);
		writeln();
		
		//~ Unique!RPointer m5 = new RPointer(mat1.solve());
		//~ write("The inverse of mat1: ");
		//~ print(m5);
		//~ writeln();
		
		//~ Unique!RPointer m6 = new RPointer(mat1.transpose());
		//~ writeln("The transpose of mat1: ");
		//~ print(m6);
		//~ writeln();
    
    /* Use MatrixPointer for convenience. */
    
	}
	foo();
  writeln("Done demonstrating. All freeing should now be done.");
  
  /* Least squares */
  Robj lhs = evalR("x27450023 <- matrix(rnorm(100), ncol=1);x27450023");
  Robj rhs = evalR("x27450024 <- matrix(rnorm(300), ncol=3);x27450024");
  Robj reg = lsq_dense_Chol(rhs, lhs);
  writeln("Regression results");
  printR(reg);
  writeln();
  
  //~ Robj reg2 = lsq_dense_QR(rhs, lhs);
  //~ writeln("Second regression:");
  //~ printR(reg2);
  //~ writeln();
  
  writeln("QR decomposition of X:");
  printR(lapack_qr(rhs, 0.0000001.robj));
  writeln();
  
  Robj f = evalR("z <- function(x,y,z=2) { return(x^2) };z");
  writeln("Function f:");
  printR(f);
  writeln("Formals of f:");
  printR(FORMALS(f));
  writeln("FORMALS of f is type: ", TYPEOF(FORMALS(f)));
  writeln("Body of f:");
  printR(BODY(f));
  writeln("BODY of f is type: ", TYPEOF(BODY(f)));
  writeln("Cloenv of f:");
  printR(CLOENV(f));
  writeln("CLOENV of f is type: ", TYPEOF(CLOENV(f)));
  
  Robj mm1 = Rf_protect(evalR("matrix(seq(1.1, 4.4, 1.1), ncol=2)"));
  Robj mm2 = Rf_protect(evalR("matrix(seq(5.5, 8.8, 1.1), ncol=2)"));
  Robj mm3 = Rf_protect(evalR("matrix(0, nrow=2, ncol=2)"));
  //~ unsafe_matprod(REAL(mm1), 2, 2, REAL(mm2), 2, 2, REAL(mm3));
  printR(mm1);
  printR(mm2);
  printR(mm3);
  Rf_unprotect(3);
  auto dm1 = DenseMatrix("seq(1.4, 6.3, 0.7)", 4, 2);
  auto dm2 = DenseMatrix("seq(1.4, 6.3, 0.7)", 2, 4);
  Unique!MatrixPointer mp3 = dm1.matmul(dm2);
  mp3.print("The product of dm1 and dm2");
	closeR();
  auto mp4 = mp3.matmul(mp3);
  mp4.print("The square of mp3");
}

struct RPointer {
	Robj obj;
  
  this(Robj x) {
    obj = x;
  }

	~this() {
		writeln("Unprotecting this:");
		printR(obj);
		Rf_unprotect_ptr(obj);
	}
}

void print(ref Unique!RPointer urp) {
	printR(urp.obj);
}

struct MatrixPointer {
	Robj obj;
  
  this(Robj x) {
    obj = x;
  }
  
  void print(string msg="") {
    if (msg.length > 0) {
      writeln(msg ~ ":");
    }
    printR(obj);
  }
  
	~this() {
		writeln("Unprotecting ", obj);
		printR(obj);
		Rf_unprotect_ptr(obj);
	}

  Unique!MatrixPointer matmul(ref Unique!MatrixPointer m) {
    Unique!MatrixPointer result = new MatrixPointer(Rf_protect(dgeMatrix_matrix_mm(obj, m.obj, RFalse)));
    return result;
  }

  Unique!MatrixPointer matmul(DenseMatrix m) {
    Unique!MatrixPointer result = new MatrixPointer(Rf_protect(dgeMatrix_matrix_mm(obj, m.x, RFalse)));
    return result;
  }
}

Unique!MatrixPointer matmul(DenseMatrix m1, DenseMatrix m2) {
  Unique!MatrixPointer result = new MatrixPointer(Rf_protect(dgeMatrix_matrix_mm(m1.x, m2.x, RFalse)));
  return result;
}

Unique!MatrixPointer matmul(DenseMatrix m1, ref Unique!MatrixPointer m2) {
  Unique!MatrixPointer result = new MatrixPointer(Rf_protect(dgeMatrix_matrix_mm(m1.x, m2.obj, RFalse)));
  return result;
}
