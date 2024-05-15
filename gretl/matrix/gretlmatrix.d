/++
A D programming language library for matrix operations and matrix algebra.
+/
module gretl.matrix;
import std.conv, std.stdio, std.string, std.typecons, std.utf;
import gretltypes, matrix;

private struct DocPage {}

/++
 + Some advantages of this library:
 + 
 + $(LIST
 + * Small size
 + * Can make changes to suit your needs
 + * No need to keep up with the version of libgretl installed on your machine
 + * No need to deal with shared libraries)
 +/
DocPage Advantages;

/++
A struct that wraps a pointer to a `gretl_matrix` struct allocated by 
libgretl. 
 
You will not normally work directly with RawGretlMatrix. You will instead
work with the [RCMatrix] type, which is a reference counted wrapper struct
over RawGretlMatrix.

Examples:
---
auto m = RCMatrix(2,2);
m[0, 0..$] = [1.1, 3.3];
m[1, 0..$] = [4.4, 1.1];
m.print("This is a new matrix");

writeln("The determinant of m is: ", m.det);
writeln("The inverse of m is: ", m.inv);
---

+/
struct RawGretlMatrix {
  gretl_matrix * gm;
  
  /++
   +/
  this(gretl_matrix * _gm) {
    gm = _gm;
  }
  
  /++
   +/
  this(long rows, long cols) {
    gm = gretl_matrix_alloc(rows.to!int, cols.to!int);
  }
  
  this(double[] v, long rows, long cols) {
    assert(v.length == rows*cols, "Dimensions don't match");
    gm = gretl_matrix_alloc(rows.to!int, cols.to!int);
    foreach(ii, a; v) {
      gm.val[ii.to!int] = a;
    }
  }
  
  /++
   + This is used by SafeRefCounted to free the underlying struct allocated by libgretl. Never called by the user.
   +/
  ~this() {
    gretl_matrix_free(gm);
  }
  
  /++
   + Duplicates the matrix. 
   + 
   + Returns:
   + A newly allocated [RCMatrix] with all data copied into it.
   + 
   + Examples:
   + ---
   + m.dup;
   + ---
   +/
  RCMatrix dup() {
    auto result = RCMatrix(gm.rows, gm.cols);
    gretl_matrix_copy_data(result.gm, this.gm);
    return result;
  }

  /++
   + Returns the value of an element from the underlying data array of the matrix.
   + 
   + 
   + Returns the value of an element from the underlying data array of the matrix. Data is stored by column. Element (0,0) is the first element of the data array, element (1,0) is the second element, and so on.
   + 
   + $(TIP This is occasionally useful, with one example being the copying of data from one matrix to another manually, but most of the time you'll want to use the conventional two-dimension syntax.)
   + 
   + $(PITFALL Unlike some languages used for scientific computing, including R, Matlab, and Fortran, this library follows the traditional practice in the C community of zero-based indexing. The first element of the data array underlying a (5x3) matrix has index 0 and the last has index 14.)
   + 
   + Example:
   + ---
   + auto m = RCMatrix(4,4);
   + 
   + // Pull element 14 from the data array. 
   + // Same as m[2,3].
   + double el = m[14];
   + ---
   +/
  double opIndex(long ii) {
    assert(ii >= 0, "Array index cannot be negative");
    assert(ii < gm.rows*gm.cols, "Array index exceeds the number of elements");
    return gm.val[ii.to!int];
  }
  
  /++
   + Get the value of a matrix element using the conventional two-dimension syntax, where the first element is the row and the second element is the column.
   + 
   + $(TIP Unlike some languages used for scientific computing, including R, Matlab, and Fortran, this library follows the traditional practice in the C community of zero-based indexing. A (4x4) matrix has upper left element (0,0) and bottom right element (3,3).)
   + 
   + $(PITFALL Although your program will stop with an error message if you attempt to get the value of element (4,4) from a (4x4) matrix, it's only an error because you've specified an index that's too large. Requesting element (1,1) is a valid operation with a (4x4) matrix, so if you're using one-based indexing, the fact that your program runs does not imply that it's correct.)
   + 
   + Example:
   + ---
   + auto m = RCMatrix(4,4);
   + double el = m[2,2];
   + 
   + // Your program will die with an out-of-bounds error if you do this
   + double el = m[4,4];
   + ---
   +/
  double opIndex(long rr, long cc) {
		assert(rr >= 0, "Row index cannot be negative");
		assert(cc >= 0, "Column index cannot be negative");
		assert(rr.to!int < gm.rows, "Row index out of bounds");
		assert(cc.to!int < gm.cols, "Column index out of bounds");
		return gm.val[to!int(cc*gm.rows+rr)];
	}
	
	double opIndex(Tuple!(int, int) ind) {
		return opIndex(ind.expand);
	}
	
	double opIndex(Element el) {
		return opIndex(el.index.expand);
	}
	
	struct SliceIndex {
		long start;
		long end;
	}
	
	MatrixRef opIndex(SliceIndex ind0, SliceIndex ind1) {
		return MatrixRef(this, ind0.start, ind0.end, ind1.start, ind1.end);
	}
	
	MatrixRef opIndex(long ind0, SliceIndex ind1) {
		return MatrixRef(this, ind0, ind0, ind1.start, ind1.end);
	}
	
	MatrixRef opIndex(SliceIndex ind0, long ind1) {
		return MatrixRef(this, ind0.start, ind0.end, ind1, ind1);
	}

  void opIndexAssign(double a, long ii) {
    assert(ii >= 0, "Array index cannot be negative");
    assert(ii < gm.rows*gm.cols, "Array index exceeds the number of elements");
    gm.val[ii.to!int] = a;
  }

  void opIndexAssign(double a, long rr, long cc) {
    assert(rr >= 0, "Array index cannot be negative");
    assert(cc >= 0, "Array index cannot be negative");
    assert(rr < gm.rows, "Row index out of bounds");
    assert(cc < gm.cols, "Column index out of bounds");
    gm.val[to!int(cc*gm.rows+rr)] = a;
  }
  
	void opIndexAssign(T)(T obj, SliceIndex ind0, SliceIndex ind1) {
		auto tmp = MatrixRef(this, ind0.start, ind0.end, ind1.start, ind1.end);
		tmp = obj;
	}
	
	void opIndexAssign(T)(T obj, long ind0, SliceIndex ind1) {
		auto tmp = MatrixRef(this, ind0, ind0, ind1.start, ind1.end);
		tmp = obj;
	}
	
	void opIndexAssign(T)(T obj, SliceIndex ind0, long ind1) {
		auto tmp = MatrixRef(this, ind0.start, ind0.end, ind1, ind1);
		tmp = obj;
	}

  int opDollar(size_t dim)() {
		static if(dim == 0) {
			return gm.rows;
		} else static if(dim == 1) {
			return gm.cols;
		} else {
			assert(false, "Cannot have a dimension greater than two on a matrix.");
		}
	}
  
  //~ Why doesn't this throw an error at compile time?
  //~ It throws an error at runtime saying opSlice has not been defined
	//~ Tuple!(long, long) opSlice(size_t dim)(long s, long e) {
		//~ return SliceIndex(s, e);
	//~ }

	SliceIndex opSlice(size_t dim)(long s, long e) {
		return SliceIndex(s, e);
	}

  RCMatrix opBinary(string op: "+")(RCMatrix m) {
    assert(this.rows == m.rows, "Number of rows does not match");
    assert(this.cols == m.cols, "Number of columns does not match");
    auto result = RCMatrix(this.rows, this.cols);
    gretl_matrix_add(this.gm, m.gm, result.gm);
    return result;
  }
  
  RCMatrix opBinary(string op: "-")(RCMatrix m) {
    assert(this.rows == m.rows, "Number of rows does not match");
    assert(this.cols == m.cols, "Number of columns does not match");
    auto result = RCMatrix(this.rows, this.cols);
    gretl_matrix_subtract(this.gm, m.gm, result.gm);
    return result;
  }
  
  RCMatrix opBinary(string op: "*")(RCMatrix m) {
    int err;
    return RCMatrix(gretl_matrix_multiply_new(this.gm, m.gm, &err));
  }
  
  RCMatrix opBinary(string op: "/")(RCMatrix m) {
    int err;
    return RCMatrix(gretl_matrix_divide(this.gm, m.gm, GRETL_MOD_NONE, &err));
  }
  
  double[] array() {
		return gm.val[0..gm.rows*gm.cols];
	}
  
  void val(double a) {
    gm.val[0..gm.rows*gm.cols] = a;
  }
  
  void val(double[] v) {
    assert(v.length == gm.rows*gm.cols, "Non-conformable array");
    foreach(ii, a; v) {
      gm.val[ii.to!int] = a;
    }
  }
  
  double * ptr() {
    return gm.val;
  }
  
  int is_complex() {
    return gm.is_complex;
  }
  
  double[] opSlice() {
    return gm.val[0..gm.rows*gm.cols];
  }
  
  int rows() {
    return gm.rows;
  }
  
  int cols() {
    return gm.cols;
  }
  
  int length() {
    return gm.rows*gm.cols;
  }
  
  matrix_info * info() {
		return gm.info;
	}
	
	void diag(double[] v) {
	}
	
	MatrixElements el() {
		return MatrixElements(gm);
	}
	
	auto range() {
		return MatrixElements(gm).range();
	}
	
	void t1(long val) {
		gretl_matrix_set_t1(this.gm, val.to!int);
	}
	
	void t2(long val) {
		gretl_matrix_set_t2(this.gm, val.to!int);
	}
	
	int t1() {
		return gretl_matrix_get_t1(this.gm);
	}
	
	int t2() {
		return gretl_matrix_get_t2(this.gm);
	}
	
	bool isDated() {
		return gretl_matrix_is_dated(this.gm).to!bool;
	}
  
  bool isIdentity() {
    return gretl_is_identity_matrix(this.gm).to!bool;
  }
  
  bool isZero() {
    return gretl_is_zero_matrix(this.gm).to!bool;
  }
  
  RCMatrix finiteElements() {
    int err;
    return RCMatrix(gretl_matrix_isfinite(this.gm, &err));
  }
  
  // Ugly but supposedly faster this way
  bool opEquals()(auto ref const RCMatrix m) const {
    int err;
    return gretl_matrices_are_equal(this.gm, m.gm, double.epsilon, &err).to!bool;
  }
  
  RCMatrix cov() {
    int err;
    return RCMatrix(gretl_covariance_matrix(this.gm, false.to!bool, 1, &err));
  } 
  
  RCMatrix cor() {
    int err;
    return RCMatrix(gretl_covariance_matrix(this.gm, true.to!bool, 1, &err));
  } 
  
  RCMatrix shape(long rr, long cc) {
    int err;
    return RCMatrix(gretl_matrix_shape(this.gm, rr.to!int, cc.to!int, &err));
  }
  
  RCMatrix trimRows(long trimTop, long trimBottom) {
    int err;
    return RCMatrix(gretl_matrix_trim_rows(this.gm, trimTop.to!int,
      trimBottom.to!int, &err));
  }
  
  RCMatrix rowMin() {
    int err;
    return RCMatrix(gretl_matrix_minmax(this.gm, 0, 0, 0, true, &err));
  }
  
  RCMatrix rowMax() {
    int err;
    return RCMatrix(gretl_matrix_minmax(this.gm, 1, 0, 0, true, &err));
  }
  
  RCMatrix colMin() {
    int err;
    return RCMatrix(gretl_matrix_minmax(this.gm, 0, 1, 0, true, &err));
  }
  
  RCMatrix colMax() {
    int err;
    return RCMatrix(gretl_matrix_minmax(this.gm, 1, 1, 0, true, &err));
  }
  
  RCMatrix rowMinIndex() {
    int err;
    return RCMatrix(gretl_matrix_minmax(this.gm, 0, 0, 1, true, &err));
  }
  
  RCMatrix rowMaxIndex() {
    int err;
    return RCMatrix(gretl_matrix_minmax(this.gm, 1, 0, 1, true, &err));
  }
  
  RCMatrix colMinIndex() {
    int err;
    return RCMatrix(gretl_matrix_minmax(this.gm, 0, 1, 1, true, &err));
  }
  
  RCMatrix colMaxIndex() {
    int err;
    return RCMatrix(gretl_matrix_minmax(this.gm, 1, 1, 1, true, &err));
  }
  
  double min() {
    int err;
    return gretl_matrix_global_minmax(this.gm, 0, &err);
  }
  
  double max() {
    int err;
    return gretl_matrix_global_minmax(this.gm, 1, &err);
  }
  
  double sum() {
    int err;
    return gretl_matrix_global_sum(this.gm, &err);
  }
  
  /++
   + Returns a newly allocated matrix holding the first p principal 
   + components of matrix m.
   +/
  RCMatrix pca(RCMatrix m, long p) {
    int err;
    return RCMatrix(gretl_matrix_pca(m.gm, p.to!int, OPT_NONE, &err));
  }
  
  RCMatrix selectRows(RCMatrix m, RCMatrix sel) {
    int err;
    return RCMatrix(gretl_matrix_bool_sel(m.gm, sel.gm, 1, &err));
  }
  
  RCMatrix selectColumns(RCMatrix m, RCMatrix sel) {
    int err;
    return RCMatrix(gretl_matrix_bool_sel(m.gm, sel.gm, 0, &err));
  }
  
  RCMatrix sort(RCMatrix m, long col) {
    int err;
    return RCMatrix(gretl_matrix_sort_by_column(m.gm, col.to!int, &err));
  }
  
  RCMatrix vectorSort(RCMatrix m, bool descending=false) {
    int err;
    return RCMatrix(gretl_vector_sort(m.gm, descending.to!int, &err));
  }

	void colnames(char** s) {
		gretl_matrix_set_colnames(this.gm, s);
	}

	void colnames(string[] names) {
		assert(names.length == this.cols, "Number of names needs to match number of columns");
		auto cstrings = new char*[names.length];
		foreach(ii, name; names) {
			cstrings[ii] = cast(char*) name.toStringz();
		}
		colnames(cstrings.ptr);
	}

	void rownames(char** s) {
		gretl_matrix_set_rownames(this.gm, s);
	}	

	void rownames(string[] names) {
		assert(names.length == this.rows, "Number of names needs to match number of rows");
		auto cstrings = new char*[names.length];
		foreach(ii, name; names) {
			cstrings[ii] = cast(char*) name.toStringz();
		}
		rownames(cstrings.ptr);
	}

	string[] colnames() {
		char** names = cast(char**) gretl_matrix_get_colnames(this.gm);
		string[] result;
		foreach(ii; 0..this.cols) {
			result ~= names[ii].to!string;
		}
		return result;
	}

	string[] rownames() {
		char** names = cast(char**) gretl_matrix_get_rownames(this.gm);
		string[] result;
		foreach(ii; 0..this.rows) {
			result ~= names[ii].to!string;
		}
		return result;
	}
}
/++
 + A reference counted struct that wraps [RawGretlMatrix]. Automatically
 + frees a `gretl_matrix *` that was allocated by libgretl when the
 + reference count hits zero.
 +/
alias RCMatrix = SafeRefCounted!RawGretlMatrix;

void print(RCMatrix m, string msg="") {
	gretl_matrix_print(m.gm, toUTFz!(char*)(msg));
}

RCMatrix inv(RCMatrix m) {
	auto result = m.dup;
	gretl_invert_matrix(result.gm);
	return result;
}

double det(RCMatrix m) {
	assert(m.rows == m.cols, "det requires a square matrix");
	int err;
	return gretl_matrix_determinant(m.gm, &err);
}

double logdet(RCMatrix m) {
	assert(m.rows == m.cols, "logdet requires a square matrix");
	int err;
	return gretl_matrix_log_determinant(m.gm, &err);
}

double logabsdet(RCMatrix m) {
	assert(m.rows == m.cols, "logabsdet requires a square matrix");
	int err;
	return gretl_matrix_log_abs_determinant(m.gm, &err);
}

RCMatrix solve(RCMatrix m1, RCMatrix m2) {
	auto a = m1.dup;
	auto b = m2.dup;
	gretl_matrix_solve(a.gm, b.gm);
	return b;
}

RCMatrix multiplyMod(string mod1, string mod2)(RCMatrix m1, RCMatrix m2) {
	static if (mod1 == "none") {
		int nrow = m1.rows;
		GretlMatrixMod _mod1 = GRETL_MOD_NONE;
	}
	static if (mod1 == "t") {
		int nrow = m1.cols;
		GretlMatrixMod _mod1 = GRETL_MOD_TRANSPOSE;
	}
	static if (mod2 == "none") {
		int ncol = m2.cols;
		GretlMatrixMod _mod2 = GRETL_MOD_NONE;
	}
	static if (mod2 == "t") {
		int ncol = m2.rows;
		GretlMatrixMod _mod2 = GRETL_MOD_TRANSPOSE;
	}
	GretlMatrixMod _mod3 = GRETL_MOD_NONE;
	auto result = RCMatrix(nrow, ncol);
	gretl_matrix_multiply_mod(m1.gm, _mod1, m2.gm, _mod2, result.gm, _mod3);
	return result;
}

void multiplyMod(string mod1, string mod2, string mod3)(RCMatrix m1, RCMatrix m2, RCMatrix m3) {
	static if (mod1 == "none") {
		int nrow = m1.rows;
		GretlMatrixMod _mod1 = GRETL_MOD_NONE;
	}
	static if (mod1 == "t") {
		int nrow = m1.cols;
		GretlMatrixMod _mod1 = GRETL_MOD_TRANSPOSE;
	}
	static if (mod2 == "none") {
		int ncol = m2.cols;
		GretlMatrixMod _mod2 = GRETL_MOD_NONE;
	}
	static if (mod2 == "t") {
		int ncol = m2.rows;
		GretlMatrixMod _mod2 = GRETL_MOD_TRANSPOSE;
	}
	static if (mod3 == "none") {
		GretlMatrixMod _mod3 = GRETL_MOD_NONE;
	}
	static if (mod3 == "cumulate") {
		GretlMatrixMod _mod3 = GRETL_MOD_CUMULATE;
	}
	static if (mod3 == "decrement") {
		GretlMatrixMod _mod3 = GRETL_MOD_DECREMENT;
	}
	gretl_matrix_multiply_mod(m1.gm, _mod1, m2.gm, _mod2, m3.gm, _mod3);
}

void eyeKron(long p, RCMatrix m, RCMatrix result) {
	gretl_matrix_I_kronecker(p.to!int, m.gm, result.gm);
}

RCMatrix eyeKron(long p, RCMatrix m) {
	int err;
	return RCMatrix(gretl_matrix_I_kronecker_new(p.to!int, m.gm, &err));
}

void kronEye(long p, RCMatrix m, RCMatrix result) {
	gretl_matrix_kronecker_I(m.gm, p.to!int, result.gm);
}

RCMatrix kronEye(long p, RCMatrix m) {
	int err;
	return RCMatrix(gretl_matrix_kronecker_I_new(m.gm, p.to!int, &err));
}

RCMatrix kron(RCMatrix m1, RCMatrix m2) {
	int err;
	return RCMatrix(gretl_matrix_kronecker_product_new(m1.gm, m2.gm, &err));
}

void kron(RCMatrix m1, RCMatrix m2, RCMatrix result) {
	gretl_matrix_kronecker_product(m1.gm, m2.gm, result.gm);
}

RCMatrix hdproduct(RCMatrix m1, RCMatrix m2) {
	int err;
	return RCMatrix(gretl_matrix_hdproduct_new(m1.gm, m2.gm, &err));
}

RCMatrix pow(RCMatrix m, double k) {
	int err;
	return RCMatrix(gretl_matrix_pow(m.gm, k, &err));
}

double dotProduct(string _mod1="none", string _mod2="none")(RCMatrix m1, RCMatrix m2) {
	int err;
	static if (_mod1 == "none") {
		GretlMatrixMod mod1 = GRETL_MOD_NONE;
	} else if (_mod1 == "t") {
		GretlMatrixMod mod1 = GRETL_MOD_TRANSPOSE;
	} else {
		static assert(0);
	}
	static if (_mod2 == "none") {
		GretlMatrixMod mod2 = GRETL_MOD_NONE;
	} else if (_mod2 == "t") {
		GretlMatrixMod mod2 = GRETL_MOD_TRANSPOSE;
	} else {
		static assert(0);
	}
	return gretl_matrix_dot_product(m1.gm, mod1, m2.gm, mod2, &err);
}

// fun is "+", "*", or "mean"
// dim is "row" or "col"
RCMatrix vectorStat(string fun, string dim)(RCMatrix m) {
	static if (fun == "+") {
		GretlVecStat vs = V_SUM;
	} else static if (fun == "*") {
		GretlVecStat vs = V_PROD;
	} else static if (fun == "mean") {
		GretlVecStat vs = V_MEAN;
	} else {
		static assert(false, `In vectorStat: fun has to be "+", "*", or "mean"`);
	}
	static if ( (dim == "row") || (dim == "rows") ) {
		int rowwise = true;
	} else static if ( (dim == "col") || (dim == "cols") ) {
		int rowwise = false;
	} else {
		static assert(false, `In vectorStat: dim has to be "row" or "col"`);
	}
	int err;
	return RCMatrix(gretl_rmatrix_vector_stat(m.gm, vs, rowwise, true, &err));
}

RCMatrix rowSums(RCMatrix m) {
	return vectorStat!("+", "rows")(m);
}
	
RCMatrix rowProducts(RCMatrix m) {
	return vectorStat!("*", "rows")(m);
}
	
RCMatrix rowMeans(RCMatrix m) {
	return vectorStat!("mean", "rows")(m);
}
	
RCMatrix colSums(RCMatrix m) {
	return vectorStat!("+", "cols")(m);
}
	
RCMatrix colProducts(RCMatrix m) {
	return vectorStat!("*", "cols")(m);
}
	
RCMatrix colMeans(RCMatrix m) {
	return vectorStat!("mean", "cols")(m);
}

RCMatrix colStdDev(RCMatrix m, long df=0) {
	int err;
	int _df = df.to!int;
	return RCMatrix(gretl_matrix_column_sd (m.gm, _df, true, &err));
}

RCMatrix rowDemean(RCMatrix m) {
	auto result = m.dup;
	gretl_matrix_demean_by_row(result.gm);
	return result;
}

RCMatrix centerColumns(RCMatrix m) {
	auto result = m.dup;
	gretl_matrix_center(result.gm, true);
	return result;
}

RCMatrix standardize(RCMatrix m, long df=0) {
	int _df = df.to!int;
	auto result = m.dup;
	gretl_matrix_standardize(m.gm, _df, true);
	return result;
}

RCMatrix leftDivide(RCMatrix m1, RCMatrix m2) {
	int err;
	return RCMatrix(gretl_matrix_divide(m1.gm, m2.gm, GRETL_MOD_NONE, &err));
}

RCMatrix rightDivide(RCMatrix m1, RCMatrix m2) {
	int err;
	return RCMatrix(gretl_matrix_divide(m1.gm, m2.gm, GRETL_MOD_TRANSPOSE, &err));
}

double rcond(RCMatrix m) {
	int err;
	return gretl_general_matrix_rcond(m.gm, &err);
}

double condIndex(RCMatrix m) {
  int err;
  return gretl_matrix_cond_index(m.gm, &err);
}

RCMatrix chol(RCMatrix m) {
  auto result = m.dup;
  gretl_matrix_cholesky_decomp(result.gm);
  return result;
}

RCMatrix psdRoot(RCMatrix m) {
  auto result = m.dup;
  int err = gretl_matrix_psd_root(result.gm, false);
  return result;
}

int qrRank(RCMatrix m) {
  int err;
  return gretl_check_QR_rank(m.gm, &err, null);
}

int rank(RCMatrix m) {
  int err;
  return gretl_matrix_rank(m.gm, double.nan, &err);
}

struct QR {
  RCMatrix Q;
  RCMatrix R;
}

QR qr(RCMatrix m) {
  auto q = m.dup;
  auto r = RCMatrix(m.rows, m.cols);
  gretl_matrix_QR_decomp(q.gm, r.gm);
  return QR(q, r);
}

struct OLSFit {
  RCMatrix y;
  RCMatrix x;
  RCMatrix coeff;
  RCMatrix vcv;
  RCMatrix residuals;
  double s2;
}

OLSFit ols(RCMatrix y, RCMatrix x) {
  auto b = RCMatrix(x.cols, 1);
  auto vcv = RCMatrix(x.cols, x.cols);
  auto uhat = RCMatrix(y.rows, 1);
  double s2;
  gretl_matrix_ols(y.gm, x.gm, b.gm, vcv.gm, uhat.gm, &s2);
  return OLSFit(y, x, b, vcv, uhat, s2);
}

struct MultiOLSFit {
  RCMatrix y;
  RCMatrix x;
  RCMatrix coeff;
  RCMatrix residuals;
  RCMatrix xtxi;
}

MultiOLSFit olsMulti(RCMatrix y, RCMatrix x) {
  auto b = RCMatrix(x.cols, y.cols);
  auto e = RCMatrix(y.rows, y.cols);
  gretl_matrix * xtxi;
  gretl_matrix_multi_ols(y.gm, x.gm, b.gm, e.gm, &xtxi);
  return MultiOLSFit(y, x, b, e, RCMatrix(xtxi));
}

OLSFit olsRestricted(RCMatrix y, RCMatrix x, RCMatrix R, RCMatrix q) {
  auto b = RCMatrix(x.cols, 1);
  auto vcv = RCMatrix(x.cols, x.cols);
  auto uhat = RCMatrix(y.rows, 1);
  double s2;
  gretl_matrix_restricted_ols(y.gm, x.gm, R.gm, q.gm, b.gm, vcv.gm, uhat.gm, &s2);
  return OLSFit(y, x, b, vcv, uhat, s2);
}

MultiOLSFit olsMultiRestricted(RCMatrix y, RCMatrix x, RCMatrix R, RCMatrix q) {
  auto b = RCMatrix(x.cols, y.cols);
  auto e = RCMatrix(y.rows, y.cols);
  gretl_matrix * xtxi;
  gretl_matrix_restricted_multi_ols(y.gm, x.gm, R.gm, q.gm, b.gm, e.gm, &xtxi);
  return MultiOLSFit(y, x, b, e, RCMatrix(xtxi));
}

double rsq(RCMatrix y, RCMatrix x, RCMatrix b) {
  int err;
  return gretl_matrix_r_squared(y.gm, x.gm, b.gm, &err);
}

RCMatrix columnwiseProduct(RCMatrix m1, RCMatrix m2, RCMatrix m3) {
  auto result = RCMatrix(m1.rows, m1.cols*m2.cols);
  gretl_matrix_columnwise_product(m1.gm, m2.gm, m3.gm, result.gm);
  return result;
}

RCMatrix mirror(char uplo)(RCMatrix m) {
	auto result = m.dup;
	static if( (uplo != 'U') && (uplo != 'L') ) {
		static assert(false, "uplo needs to be 'U' to copy the upper triangle to the lower triangle, or 'L' for the converse");
	}
	gretl_matrix_mirror(result.gm, uplo);
	return result;
}

RCMatrix eigenval(RCMatrix m) {
	int err;
	return RCMatrix(gretl_general_matrix_eigenvals(m.gm, &err));
}

double minEigenval(RCMatrix m) {
	int err;
	return gretl_symm_matrix_lambda_min(m.gm, &err);
}

double maxEigenval(RCMatrix m) {
	int err;
	return gretl_symm_matrix_lambda_max(m.gm, &err);
}

struct Eigen {
	RCMatrix values;
	RCMatrix vectors;
}

Eigen eigen(RCMatrix m) {
	auto vecs = m.dup;
	int err;
	Eigen result;
	result.values = RCMatrix(gretl_symmetric_matrix_eigenvals(m.gm, true, &err));
	result.vectors = vecs;
	return result;
}

struct SVD {
	RCMatrix u;
	RCMatrix s;
	RCMatrix vt;
}

SVD svd(RCMatrix m) {
	gretl_matrix * u;
	gretl_matrix * s;
	gretl_matrix * vt;
	gretl_matrix_SVD(m.gm, &u, &s, &vt, true);
	return SVD(RCMatrix(u), RCMatrix(s), RCMatrix(vt));
}

RCMatrix rightNullspace(RCMatrix m) {
	int err;
	return RCMatrix(gretl_matrix_right_nullspace(m.gm, &err));
}

RCMatrix leftNullspace(RCMatrix m) {
	return leftNullspace!"none"(m);
}

RCMatrix leftNullspace(string mod)(RCMatrix m) {
	int err;
	static if (mod == "none") {
		return RCMatrix(gretl_matrix_left_nullspace(m.gm, GRETL_MOD_NONE, &err));
	} else static if (mod == "t") {
		return RCMatrix(gretl_matrix_left_nullspace(m.gm, GRETL_MOD_TRANSPOSE, &err));
	} else {
		static assert(false, "In leftNullspace: mod has to be \"none\" or \"t\"");
	}
}

RCMatrix rbind(RCMatrix m1, RCMatrix m2) {
	int err;
	return RCMatrix(gretl_matrix_row_concat(m1.gm, m2.gm, &err));
}

RCMatrix cbind(RCMatrix m1, RCMatrix m2) {
	int err;
	return RCMatrix(gretl_matrix_col_concat(m1.gm, m2.gm, &err));
}

RCMatrix directSum(RCMatrix m1, RCMatrix m2) {
	int err;
	return RCMatrix(gretl_matrix_direct_sum(m1.gm, m2.gm, &err));
}

RCMatrix cumsum(RCMatrix m) {
	int err;
	return RCMatrix(gretl_matrix_cumcol(m.gm, &err));
}

RCMatrix diff(RCMatrix m) {
	int err;
	return RCMatrix(gretl_matrix_diffcol(m.gm, double.nan, &err));
}

/* k is a vector of per column lag orders */
RCMatrix lag(RCMatrix m, RCMatrix k) {
	return lag!true(m, k);
}

RCMatrix lag(bool byvar)(RCMatrix m, RCMatrix k) {
	static if (byvar) {
		return RCMatrix(gretl_matrix_lag(m.gm, k.gm, OPT_NONE, double.nan));
	} else {
		return RCMatrix(gretl_matrix_lag(m.gm, k.gm, OPT_L, double.nan));
	}
}
RCMatrix olsVarCov(RCMatrix V, double s2) {
	auto result = V.dup;
	get_ols_vcv(result.gm, &s2);
	return result;
}

RCMatrix uhat(RCMatrix y, RCMatrix x, RCMatrix b) {
	auto result = RCMatrix(y.rows, 1);
	get_ols_uhat(y.gm, x.gm, b.gm, result.gm);
	return result;
}

RCMatrix moorePenrose(RCMatrix m) {
	auto result = m.dup;
	gretl_matrix_moore_penrose(m.gm, double.nan);
	return result;
}

void qform(string mod1, string mod2)(RCMatrix m1, RCMatrix m2, RCMatrix result) {
	static if (mod1 == "none") {
		GretlMatrixMod amod = GRETL_MOD_NONE;
	} else static if (mod1 == "t") {
		GretlMatrixMod amod = GRETL_MOD_TRANSPOSE;
	} else {
		static assert(false, `In qform: mod1 has to be "none" or "t"`);
	}
	static if (mod2 == "none") {
		GretlMatrixMod cmod = GRETL_MOD_NONE;
	} else static if (mod2 == "cumulate") {
		GretlMatrixMod cmod = GRETL_MOD_CUMULATE;
	} else static if (mod2 == "decrement") {
		GretlMatrixMod cmod = GRETL_MOD_DECREMENT;
	} else {
		static assert(false, `In qform: mod2 has to be "none", "cumulate", or "decrement"`);
	}
	gretl_matrix_qform(m1.gm, amod, m2.gm, result.gm, cmod);
}

RCMatrix qform(RCMatrix m1, RCMatrix m2) {
	auto result = RCMatrix(m1.rows, m1.rows);
	qform!("none", "none")(m1, m2, result);
	return result;
}

/* m1 needs to be a vector, for when I make these changes in the future */
double qformScalar(RCMatrix m1, RCMatrix m2) {
  int err;
  return gretl_scalar_qform(m1.gm, m2.gm, &err);
}

/* m1 is a vector */
RCMatrix diagSandwich(RCMatrix m1, RCMatrix m2) {
  auto result = RCMatrix(m2.rows, m2.rows);
  gretl_matrix_diagonal_sandwich(m1.gm, m2.gm, result.gm);
  return result;
}

/* Weights is a vector */
RCMatrix weightedCovariagram(RCMatrix m1, RCMatrix m2, long p, RCMatrix weights) {
  int err;
  return RCMatrix(gretl_matrix_covariogram(m1.gm, m2.gm, weights.gm, p.to!int, &err));
}

RCMatrix covariagram(RCMatrix m1, RCMatrix m2, long p) {
  int err;
  return RCMatrix(gretl_matrix_covariogram(m1.gm, m2.gm, null, p.to!int, &err));
}

RCMatrix gginv(RCMatrix m) {
	int err;
	return RCMatrix(gretl_matrix_GG_inverse(m.gm, &err));
}


struct MatrixElements {
	int rows;
	int cols;
	double[] val;
	
	this(gretl_matrix * m) {
		rows = m.rows;
		cols = m.cols;
		val = m.val[0..m.rows*m.cols];
	}
	
	double[] array() {
		return this.val;
	}
	
	double opIndex(int rr, int cc) {
		return this.val[cc*rows + rr];
	}
	
	RCMatrix opBinary(string op, T)(T m) 
	if (is(T==RCMatrix) || is(T==MatrixElements)) {
		auto result = RCMatrix(m.rows, m.cols);
		static if(op == "+") {
			result.array[] = this.array[] + m.array[];
		}
		static if(op == "-") {
			result.array[] = this.array[] - m.array[];
		}
		static if(op == "*") {
			result.array[] = this.array[] * m.array[];
		}
		static if(op == "/") {
			result.array[] = this.array[] / m.array[];
		}
		return result;
	}
	
	RCMatrix opBinaryRight(string op)(RCMatrix m) {
		auto result = RCMatrix(m.rows, m.cols);
		static if(op == "+") {
			result.array[] = m.array[] + this.array[];
		}
		static if(op == "-") {
			result.array[] = m.array[] - this.array[];
		}
		static if(op == "*") {
			result.array[] = m.array[] * this.array[];
		}
		static if(op == "/") {
			result.array[] = m.array[] / this.array[];
		}
		return result;
	}
	
	auto diag() {
		struct ElementRange {
			int currentRow = 0;
			int currentColumn = 0;
			int rows;
			int cols;
			double[] val;
			
			this(int rr, int cc, double[] v) {
				rows = rr;
				cols = cc;
				val = v;
			}
				
			bool empty() {
				return currentColumn >= cols;
			}
			
			Element front() {
				Element result;
				result.index = tuple(currentRow, currentColumn);
				result.value = val[currentColumn*rows+currentRow];
				return result;
			}
			
			void popFront() {
				currentRow += 1;
				currentColumn += 1;
			}
		}
		
		return ElementRange(this.rows, this.cols, this.val);		
	}
	
	/* Returns a range of all elements */
	auto range() {
		struct ElementRange {
			int currentRow = 0;
			int currentColumn = 0;
			int rows;
			int cols;
			double[] val;
			
			this(int rr, int cc, double[] v) {
				rows = rr;
				cols = cc;
				val = v;
			}
				
			bool empty() {
				return currentColumn >= cols;
			}
			
			Element front() {
				Element result;
				result.index = tuple(currentRow, currentColumn);
				result.value = val[currentColumn*rows+currentRow];
				return result;
			}
			
			void popFront() {
				currentRow += 1;
				if (currentRow >= rows) {
					currentRow = 0;
					currentColumn += 1;
				}
			}
		}
		
		return ElementRange(this.rows, this.cols, this.val);		
	}
}

struct Element {
	import std.typecons;
	Tuple!(int, int) index;
	double value;
}

struct MatrixRef {
	long startRow;
	// Excluded
	long endRow;
	long startColumn;
	// Excluded
	long endColumn;
	long sourceRows;
	double[] sourceData;
	
	alias dup this;
	
	this(RCMatrix m, long sr, long er, long sc, long ec) {
		startRow = sr;
		endRow = er;
		startColumn = sc;
		endColumn = ec;
		sourceRows = m.rows;
		sourceData = m.array;
	}
	
	this(gretl_matrix * gm, long sr, long er, long sc, long ec) {
		startRow = sr;
		endRow = er;
		startColumn = sc;
		endColumn = ec;
		sourceRows = gm.rows;
		sourceData = gm.val[0..gm.rows*gm.cols];
	}
	
	// Need ref to prevent double free
	this(ref RawGretlMatrix m, long sr, long er, long sc, long ec) {
		this(m.gm, sr, er, sc, ec);
	}		
			
	this(RCMatrix m) {
		startRow = 0;
		endRow = m.rows;
		startColumn = 0;
		endColumn = m.cols;
		sourceRows = m.rows;
		sourceData = m.array;
	}
	
	long rows() {
		return endRow-startRow;
	}
	
	long cols() {
		return endColumn-startColumn;
	}
	
	long length() {
		return rows()*cols();
	}
	
	long sourceIndex(long rr, long cc) {
		return cc*sourceRows + rr;
	}
	
	long trueIndex(long rr, long cc) {
		return (cc+startColumn)*sourceRows + rr+startRow;
	}

	double opIndex(long r, long c) {
		return sourceData[trueIndex(r, c)];
	}
	
	void opIndexAssign(double val, long r, long c) {
		sourceData[trueIndex(r, c)] = val;
	}
	
	void opAssign(double a) {
		foreach(col; startColumn..endColumn) {
			foreach(row; startRow..endRow) {
				sourceData[sourceIndex(row, col)] = a;
			}
		}
	}
	
	void opAssign(RCMatrix m) {
		assert(this.rows == m.rows, "Rows don't match");
		assert(this.cols == m.cols, "Columns don't match");
		sourceData[] = m.array[];
	}
	
	void opAssign(MatrixRef mr) {
		assert(this.rows == mr.rows, "Rows don't match");
		assert(this.cols == mr.cols, "Columns don't match");
		foreach(col; 0..this.cols) {
			foreach(row; 0..this.rows) {
				this[row, col] = mr[row, col];
			}
		}
	}
	
	RCMatrix dup() {
		auto result = RCMatrix(this.rows, this.cols);
		foreach(col; 0..this.cols) {
			foreach(row; 0..this.rows) {
				result[row, col] = this[row, col];
			}
		}
		return result;
	}
}

struct SubmatrixElements {
	/* The row and column elements can individually be specified as a
	 * scalar, sequence, or array. 
	 * 
	 * Iteration produces all elements.*/
}

/* .ref returns a reference to the matrix. That means there's no copying
 * of data if you take a slice. While it's faster, there are downsides:
 * 
 * - Access is invalid if the original matrix is freed.
 * - You have to remember that you're changing the original matrix
 *   as well as the reference.
 * - If you pass a reference to a matrix to other functions, things
 *   get complicated. You essentially have global data if you do that.
 * - You don't have full access to everything. For instance, there's no
 *   direct access granted to the underlying gretl_matrix.
 * 
 * Use references locally and only if the performance is really
 * needed.
 */

/* .el says to treat the matrix "element-by-element". Multiplication
 * becomes element-by-element multiplication and so on.
 * 
 * You can also use a foreach to iterate over the elements. Examples:
 * 
 * - foreach(index; m.el) will iterate over all of the elements of the
 *   matrix.
 * - foreach(index; diag(m.el)) iterates over all elements of the
 *   diagonal.
 * - foreach(index; diagAndBelow(m.el)) iterates over all elements on
 *   the diagonal and below.
 * - foreach(index; column(m.el, 2)) iterates over all elements of
 *   column 2.
 * - foreach(index; m.el[1..4, 2..6]) iterates over all elements of
 *   that submatrix.
 * 
 * In these examples, index is a tuple holding the two index values.
 */
