module betterr.matrix;

import betterr.list, betterr.rdata, betterr.vector;
import betterr.r;
import std.algorithm, std.array, std.conv, std.exception;
import std.range, std.stdio, std.traits;
struct Matrix {
  long rows;
  long cols;
  RData data;
  double * ptr;
  alias data this;
  
  this(string code) {
    data = RData(code);
    rows = Rf_nrows(data.x);
    cols = Rf_ncols(data.x);
    ptr = REAL(data.x);
  }

  this(long r, long c) {
    this("matrix(0.0, nrow=" ~ r.to!string ~ ", ncol=" ~ c.to!string ~ ")");
  }
  
  this(Matrix m) {
    this(m.data.name);
    /* Force a copy */
    this.data.x = evalR(this.name ~ "[]");
  }
 
  this(Submatrix sm) {
    data = RData("matrix(0.0, nrow=" ~ sm.rows.to!string ~ ", ncol=" ~ sm.cols.to!string ~ ")");
    rows = sm.rows;
    cols = sm.cols;
    ptr = REAL(data.x);
    foreach(jj; 0..cols) {
			foreach(ii; 0..rows) {
				this[ii, jj] = sm[ii, jj];
			}
		}
  }

  this(Vector v) {
    this("as.matrix(" ~ v.name ~ ")");
  }
  
  this(Vector v, long r, long c) {
    enforce(v.length == r*c, "Number of elements does not match");
    this("matrix(" ~ v.name ~ ", nrow=" ~ r.to!string ~ ", ncol=" ~ c.to!string ~ ")");
  }
  
  this(double[] v, long r, long c) {
    enforce(v.length == r*c, "Number of elements does not match");
    this(r, c);
    foreach(ii; 0..v.length) {
      this.ptr[ii] = v[ii];
    }
  }
  
  this(RData rd) {
		if (Rf_isMatrix(rd.x)) {
			this(rd.name ~ "[]");
		} else {
			this(rd.name);
		}
  }
  
  Matrix dup() {
    auto result = Matrix(this.rows, this.cols);
    foreach(ii; 0..this.rows*this.cols) {
      result.ptr[ii] = this.ptr[ii];
    }
    return result;
  }

  double opIndex(long r, long c) {
    enforce(r < this.rows, "First index exceeds the number of rows");
    enforce(c < this.cols, "Second index exceeds the number of columns");
    return ptr[c*this.rows+r];
  }

  void opIndexAssign(double v, long r, long c) {
    ptr[c*rows+r] = v;
  }
  
  long opDollar(int dim)() {
    static if (dim == 0) {
      return this.rows;
    } else if (dim == 1) {
      return this.cols;
    }
  }
  
  double opIndex(Element el) {
    return opIndex(el.row, el.column);
  }

  void opIndexAssign(double v, Element el) {
    opIndexAssign(v, el.row, el.column);
  }
	
	string indexString(Vector v) {
		return "as.integer(" ~ v.name ~ ")";
	}
  
  string indexString(long[] ind) {
    return "c(" ~ ind.map!(x => to!string(x+1)).join(", ") ~ ")";
	}

  string indexString(int[] ind) {
    return "c(" ~ ind.map!(x => to!string(x+1)).join(", ") ~ ")";
	}

  string indexString(long ind) {
		return to!string(ind+1);
	}
	
	string indexString(AllElements ae) {
		return "";
	}
	
  string indexString(long begin, long end) {
    return to!string(begin+1) ~ ":" ~ to!string(end);
  }
  
  string indexString(SliceIndex si) {
    return indexString(si.start, si.end);
  }
  
  /* Multidimensional indexing 
   * Returning a string allows long[2] to be used to specify which
   * columns/rows you want.
   * Note that m[a,b] <- x operations can lead to a reallocation!
   * Much easier, and more efficient, to do all this with a Submatrix. */
  SliceIndex opSlice(long dim)(long begin, long end) {
    return SliceIndex([begin, end]);
  }
  
  void fillWith(T)(T v) {
		enforce(this.rows * this.cols == v.length, "fillWith requires the same number of elements on both sides of =.");
		foreach(ii; 0..v.length) {
			this.ptr[ii] = v[ii];
		}
	}
  
  Matrix opIndex(T1, T2)(T1 rr, T2 cc) {
		return Matrix(this.name ~ "[" ~ indexString(rr) ~ ", " ~ indexString(cc) ~ ", drop=FALSE]");
  }
  
  void opIndexAssign(T1, T2)(double a, T1 rr, T2 cc) {
    auto sm = this.reference();
    sm[rr, cc] = a;
  }
  
  void opIndexAssign(T0, T1, T2)(T0 m, T1 rr, T2 cc) 
    if (is(T0 == Matrix) | is(T0 == Submatrix)) {
      auto sm = this.reference();
      sm[rr, cc] = m;
    }
  
	void opIndexAssign(T0, T1)(T0 v, T1 rr, long col) 
    if (is(T0 == Vector) | is(T0 == Column)
				| is (T0 == Row) | is(T0 == double[])) {
			auto sm = this.reference();
			sm[rr, col] = v;
		}

	void opIndexAssign(T0, T1)(T0 v, long row, T1 cc) 
    if (is(T0 == Vector) | is(T0 == Column)
				| is (T0 == Row) | is(T0 == double[])) {
			auto sm = this.reference();
			sm[row, cc] = v;
		}

  auto opBinary(string op, T)(T x) {
    static if (op == "~") {
      return cbind(this, x);
    }
  }

  Submatrix reference() {
    return Submatrix(this, [0, 0], [rows, cols]);
  }
  
  Submatrix sub() {
    return Submatrix(this, [0, 0], [rows, cols]);
  }
  
  Vector rowSums() {
    return Vector("rowSums(" ~ this.data.name ~ ")");
  }
  
  Vector colSums() {
    return Vector("colSums(" ~ this.data.name ~ ")");
  }
  
  Vector rowMeans() {
    return Vector("rowMeans(" ~ this.data.name ~ ")");
  }
  
  Vector colMeans() {
    return Vector("colMeans(" ~ this.data.name ~ ")");
  }
  
  Vector row(long ii) {
		return Vector(this.data.name ~ "[" ~ to!string(ii+1) ~ ",]");
	}
  
  Vector column(long ii) {
		return Vector(this.data.name ~ "[," ~ to!string(ii+1) ~ "]");
	}
	
  Vector column(string name) {
		return Vector(this.data.name ~ "[,'" ~ name ~ "']");
	}
	
	Matrix columns(long c0, long c1) {
		return Matrix(this.data.name ~ "[," ~ to!string(c0+1) ~ ":" ~ to!string(c1) ~ "]");
	}
	
	Vector lastrow() {
		return row(rows-1);
	}
	
	Vector lastcolumn() {
		return column(cols-1);
	}
  
  long length() {
    return rows*cols;
  }
}

/* f should be able to be evaluated at any positive number
 * If it can't be evaluated at 1.0 for some reason, you can
 * specify val.
 * 
 * f should take a Row as an argument.
 * 
 * Does not currently work with lambdas. */
auto mapRows(alias f)(Matrix m) {
  alias T = double[];
  static if (is(ReturnType!f == double)) {
    auto result = Vector(m.rows);
    foreach(ii; 0..m.rows) {
      result[ii] = f(Row(m, ii));
    }
    return result;
  }
  static if (is(ReturnType!f == long) | is(ReturnType!f == int)) {
    auto result = IntVector(m.rows);
    foreach(ii; 0..m.rows) {
      result[ii] = f(Row(m, ii));
    }
    return result;
  }
  static if (is(ReturnType!f == bool)) {
    auto result = BoolVector(m.rows);
    foreach(ii; 0..m.rows) {
      result[ii] = f(Row(m, ii));
    }
    return result;
  }
  static if (is(ReturnType!f == T)) {
    auto result = Matrix(m.rows, f(Row(m,0)).length);
    foreach(ii; 0..m.rows) {
      auto tmp = f(Row(m, ii));
      foreach(jj; 0..tmp.length) {
        result[ii, jj] = tmp[jj];
      }
    }
    return result;
  }
}

/* These should *not* take a Submatrix, Column, or Row as an argument.
 * Those are reference types. There needs to be a name inside R, which
 * would required *more* copying if they're used. 
 * 
 * A Vector should work though. A Vector will always be treated as a
 * matrix with one column. To treat it as a row, use v.rowMatrix. */
Matrix matmul(Matrix x, Matrix y) {
  return Matrix(x.name ~ "%*%" ~ y.name);
}

Matrix matmul(Vector v, Matrix y) {
	return matmul(v.matrix, y);
}

Matrix matmul(Matrix x, Vector v) {
	return matmul(x, v.matrix);
}

Matrix elmul(Matrix x, Matrix y) {
  return Matrix(x.name ~ "*" ~ y.name);
}

Matrix plus(Matrix x, Matrix y) {
  return Matrix(x.name ~ "+" ~ y.name);
}

Matrix minus(Matrix x, Matrix y) {
  return Matrix(x.name ~ "-" ~ y.name);
}

Matrix div(Matrix x, Matrix y) {
  return Matrix(x.name ~ "/" ~ y.name);
}

Matrix mul(Matrix x, double y) {
  return Matrix(x.name ~ "*" ~ y.to!string);
}

Matrix plus(Matrix x, double y) {
  return Matrix(x.name ~ "+" ~ y.to!string);
}

Matrix minus(Matrix x, double y) {
  return Matrix(x.name ~ "-" ~ y.to!string);
}

Matrix div(Matrix x, double y) {
  return Matrix(x.name ~ "/" ~ y.to!string);
}

Matrix mul(double y, Matrix x) {
  return Matrix(x.name ~ "*" ~ y.to!string);
}

Matrix plus(double y, Matrix x) {
  return Matrix(x.name ~ "+" ~ y.to!string);
}

Matrix minus(double y, Matrix x) {
  return Matrix(y.to!string ~ "-" ~ x.name);
}

Matrix div(double y, Matrix x) {
  return Matrix(y.to!string ~ "/" ~ x.name);
}

Matrix t(Matrix x) {
  return Matrix("t(" ~ x.name ~ ")");
}

Matrix solve(Matrix x) {
  return Matrix("solve(" ~ x.name ~ ")");
}

Matrix solve(Matrix x, Matrix y) {
	enforce(x.rows == y.rows, "The first argument to solve has to be a square matrix");
	return Matrix("solve(" ~ x.name ~ ", " ~ y.name ~ ")");
}

Vector solve(Matrix x, Vector y) {
	enforce(x.rows == y.rows, "The first argument to solve has to be a square matrix");
	return Vector("solve(" ~ x.name ~ ", " ~ y.name ~ ")");
}

Matrix inv(Matrix x) {
  return solve(x);
}

Vector diag(Matrix x) {
  return Vector("diag(" ~ x.name ~ ")");
}

Matrix kronecker(Matrix x, Matrix y) {
  return Matrix("kronecker(" ~ x.name ~ "," ~ y.name ~ ")");
}

Matrix crossprod(Matrix x, Matrix y) {
  return Matrix("crossprod(" ~ x.name ~ ", " ~ y.name ~ ")");
}

Matrix tcrossprod(Matrix x, Matrix y) {
  return Matrix("tcrossprod(" ~ x.name ~ ", " ~ y.name ~ ")");
}

Matrix crossprod(Matrix x) {
  return Matrix("crossprod(" ~ x.name ~ ")");
}

Matrix tcrossprod(Matrix x) {
  return Matrix("tcrossprod(" ~ x.name ~ ")");
}

double det(Matrix x) {
  return evalR("det(" ~ x.name ~ ")").scalar;
}

Matrix diag(long ii) {
  return Matrix("diag(" ~ ii.to!string ~ ")");
}

Matrix eye(long ii) {
  return Matrix("diag(" ~ ii.to!string ~ ")");
}

Matrix cbind(T1, T2)(T1 m, T2 v) {
  return Matrix("cbind(" ~ m.name ~ ", " ~ v.name ~ ")");
}

Matrix rbind(T1, T2)(T1 m, T2 v) {
  return Matrix("rbind(" ~ m.name ~ ", " ~ v.name ~ ")");
}

Matrix cbind(Matrix m, Vector v) {
  return Matrix("cbind(" ~ m.name ~ ", " ~ v.name ~ ")");
}

Matrix rbind(Matrix m, Vector v) {
  return Matrix("rbind(" ~ m.name ~ ", " ~ v.name ~ ")");
}

Matrix cbind(Matrix m, double[] arr) {
  return Matrix("cbind(" ~ m.name ~ ", " ~ Vector(arr).name ~ ")");
}

Matrix rbind(Matrix m, double[] arr) {
  return Matrix("rbind(" ~ m.name ~ ", " ~ Vector(arr).name ~ ")");
}

struct SVD {
  Vector d;
  Matrix u;
  Matrix v;
  
  this(Matrix m) {
    auto tmp = List("svd(" ~ m.name ~ ")");
    d = Vector(tmp["d"]);
    u = Matrix(tmp["u"]);
    v = Matrix(tmp["v"]);
  }
}

struct Eigen {
  Vector values;
  Matrix vectors;
  
  this(Matrix m) {
    auto tmp = List("eigen(" ~ m.name ~ ")");
    values = Vector(tmp["values"]);
    vectors = Matrix(tmp["vectors"]);
  }
}

/* Note: This is a reference type that exists solely for efficiency reasons. If the underlying
 matrix Robj changes, this will no longer work.  Otherwise you could just use a Vector. */
struct Column {
  long rows;
  double * ptr;
  private long currentElement = 0;
  
  this(Matrix m, long col) {
    enforce(col < m.cols, "Column number is out of range, the matrix has only " ~ m.cols.to!string ~ " columns.");
    rows = m.rows;
    ptr = &m.ptr[col*m.rows];
  }
  
  double opIndex(long r) {
    return ptr[r];
  }
  
  void opIndexAssign(double v, long r) {
    ptr[r] = v;
  }
  
  long opDollar() {
		return rows;
	}
  
  void opAssign(T)(T v) {
    enforce(rows == v.length, "Attempting to copy into Column an object with the wrong number of elements");
    foreach(ii; 0..rows) {
      ptr[ii] = v[ii];
    }
  }
 
  void opAssign(double x) {
    foreach(ii; 0..rows) { 
      ptr[ii] = x;
    }
  }
  
  void opOpAssign(string op)(double a) {
		static if (op == "+") {
			ptr[0..rows] += a;
		}
		static if (op == "-") {
			ptr[0..rows] -= a;
		}
		static if (op == "*") {
			ptr[0..rows] *= a;
		}
		static if (op == "/") {
			ptr[0..rows] /= a;
		}
	}
	
	double[] reference() {
		return ptr[0..rows];
	}

	long length() {
		return rows;
	}

  bool empty() { 
    return currentElement == rows; 
  }
  
  double front() { 
    return ptr[currentElement]; 
  }
  
  void popFront() {
    currentElement += 1;
  }
}

struct Row {
  long rows;
  long cols;
  double * ptr;
  private long currentElement = 0;
  
  this(Matrix m, long row) {
    enforce(row < m.rows, "Column number is out of range, the matrix has only " ~ m.cols.to!string ~ " columns.");
    cols = m.cols;
    rows = m.rows;
    ptr = &m.ptr[row];
  }
  
  double opIndex(long c) {
    return ptr[c*rows];
  }
  
  void opIndexAssign(double v, long c) {
    ptr[c*rows] = v;
  }
  
  void opAssign(T)(T v) {
    enforce(cols == v.length, "Attempting to copy into Row an object with the wrong number of elements");
    foreach(ii; 0..cols) {
      this[ii] = v[ii];
    }
  }
 
  void opAssign(double x) {
    foreach(ii; 0..cols) { 
      this[ii] = x;
    }
  }
  
  void opOpAssign(string op)(double a) {
		static if (op == "+") {
			foreach(ii; 0..cols) {
				this[ii] = this[ii]+a;
			}
		}
		static if (op == "-") {
			foreach(ii; 0..cols) {
				this[ii] = this[ii]-a;
			}
		}
		static if (op == "*") {
			foreach(ii; 0..cols) {
				this[ii] = this[ii]*a;
			}
		}
		static if (op == "/") {
			foreach(ii; 0..cols) {
				this[ii] = this[ii]/a;
			}
		}
	}
	
	long opDollar() {
		return cols;
	}
	
	long length() {
		return cols;
	}

  bool empty() { 
    return currentElement == cols; 
  }
  
  double front() { 
    return this[currentElement]; 
  }
  
  void popFront() {
    currentElement += 1;
  }
}

/* Submatrix holds a reference to a block of a matrix. This allows
 * access/manipulation without copying.
 * 
 * rows and cols refer to the Submatrix, so that it can be substituted
 * for a Matrix.
 * 
 * matrows is the number of rows in the underlying Matrix. We need that
 * for indexing.
 * 
 * offset is the [0,0] element of the Submatrix.
 *   */
struct Submatrix {
  long rows;
  long cols;
  long matrows;
  long offset;
  double * ptr;
  
  this(Matrix m, long[2] element0, long[2] element1) {
    auto row0 = element0[0];
    auto col0 = element0[1];
    auto row1 = element1[0];
    auto col1 = element1[1];
    enforce(row0 < row1, "row0 has to be less than row1");
    enforce(row1 <= m.rows, "Row number is out of range, the matrix has only " ~ m.rows.to!string ~ " rows.");
    enforce(col0 < col1, "col0 has to be less than col1");
    enforce(col1 <= m.cols, "Column number is out of range, the matrix has only " ~ m.cols.to!string ~ " columns.");
    matrows = m.rows;
    rows = row1-row0;
    cols = col1-col0;
    offset = col0*m.rows+row0;
    ptr = m.ptr;
  }
  
  this(Submatrix sm, long[2] element0, long[2] element1) {
    auto row0 = element0[0];
    auto col0 = element0[1];
    auto row1 = element1[0];
    auto col1 = element1[1];
    enforce(row0 < row1, "row0 has to be less than row1");
    enforce(row1 <= sm.rows, "Row number is out of range, the matrix has only " ~ sm.rows.to!string ~ " rows.");
    enforce(col0 < col1, "col0 has to be less than col1");
    enforce(col1 <= sm.cols, "Column number is out of range, the matrix has only " ~ sm.cols.to!string ~ " columns.");
    matrows = sm.matrows;
    rows = row1-row0;
    cols = col1-col0;
    offset = offset + col0*sm.matrows+row0;
    ptr = sm.ptr;
  }
  
  double opIndex(long rr, long cc) {
    /* cc*rows is the start of column cc in the submatrix */
    return ptr[offset+cc*matrows+rr];
  }
  
  void opIndexAssign(double v, long rr, long cc) {
    ptr[offset+cc*matrows+rr] = v;
  }

  double opIndex(Element el) {
    return opIndex(el.row, el.column);
  }

  void opIndexAssign(double v, Element el) {
    opIndexAssign(v, el.row, el.column);
  }
  
  long opDollar(int dim)() {
    static if (dim == 0) {
      return this.rows;
    } else if (dim == 1) {
      return this.cols;
    }
  }
  
  Submatrix reference() {
		Submatrix result;
		result.rows = rows;
		result.cols = cols;
		result.matrows = matrows;
		result.offset = offset;
		result.ptr = ptr;
    return result;
  } 

	long[2] opSlice(long dim)(long begin, long end) {
		return [begin, end];
	}
	
	void opIndexAssign(T1, T2)(double a, T1 rr, T2 cc) {
    foreach(col; SliceIndex(cc, cols)) {
      foreach(row; SliceIndex(rr, rows)) {
        this[row, col] = a;
      }
    }
	}
  
	void opIndexAssign(T0, T1, T2)(T0 m, T1 rr, T2 cc) 
    if (is(T0 == Matrix) | is(T0 == Submatrix)) {
      auto csi = SliceIndex(cc, cols);
      auto rsi = SliceIndex(rr, rows);
      enforce(rsi.length == m.rows, "Rows do not match");
      enforce(csi.length == m.cols, "Columns do not match");
      foreach(jj, col; csi.enumerate) {
        foreach(ii, row; rsi.enumerate) {
          this[row, col] = m[ii, jj];
        }
      }
		}
  
	void opIndexAssign(T0, T1)(T0 v, T1 rr, long col) 
    if (is(T0 == Vector) | is(T0 == Column)
				| is (T0 == Row) | is(T0 == double[])) {
      auto rsi = SliceIndex(rr, rows);
      enforce(rsi.length == v.length, "Length does not match");
      foreach(ii, row; rsi.enumerate) {
        this[row, col] = v[ii];
      }
		}
	
	void opIndexAssign(T0, T1)(T0 v, long row, T1 cc) 
    if (is(T0 == Vector) | is(T0 == Column)
				| is (T0 == Row) | is(T0 == double[])) {
      auto csi = SliceIndex(cc, cols);
      enforce(csi.length == v.length, "Length does not match");
      foreach(ii, col; csi.enumerate) {
        this[row, col] = v[ii];
      }
		}
	
/* This will work for a Matrix or Submatrix */
  void opAssign(T)(T m) {
    enforce(rows == m.rows, "Attempting to copy into Submatrix an object with the wrong number of rows");
    enforce(cols == m.cols, "Attempting to copy into Submatrix an object with the wrong number of columns");
    foreach(el; Elements(rows, cols)) {
      this[el] = m[el];
    }
  }
 
  void opAssign(double x) {
    foreach(el; Elements(rows, cols)) { 
      this[el] = x;
    }
  }
  
  void opOpAssign(string op)(double a) {
		static if (op == "+") {
			foreach(el; Elements(rows, cols)) {
				this[el] = this[el]+a;
			}
		}
		static if (op == "-") {
			foreach(el; Elements(rows, cols)) {
				this[el] = this[el]-a;
			}
		}
		static if (op == "*") {
			foreach(el; Elements(rows, cols)) {
				this[el] = this[el]*a;
			}
		}
		static if (op == "/") {
			foreach(el; Elements(rows, cols)) {
				this[el] = this[el]/a;
			}
		}
	}
}

struct Element {
  long row;
  long column;
}

/* Returns all elements of a matrix */
struct Elements {
  long rows;
  long cols;
  long currentRow=0;
  long currentCol=0;

  bool empty() {
    return currentCol >= cols;
  }

  Element front() {
    return Element(currentRow, currentCol);
  }

  void popFront() {
    if (currentRow >= rows-1) {
      currentCol += 1;
      currentRow = 0;
    } else {
      currentRow += 1;
    }
  }
}

/* Common parts of a matrix to fill. By default, fills by column, but
 * a template parameter will allow to fill by row. This is done for
 * efficiency and to prevent mistakes. */
struct ColumnFill {
  long rows;
  long cols;
  double * ptr;
  long column;
  long length;
  private long fillPointer = 0;
  
	this(Matrix m, long col) {
		rows = m.rows;
		cols = m.cols;
		ptr = m.ptr;
		length = rows;
		column = col;
	}
	
	void opOpAssign(string op)(double a) if (op == "~") {
		enforce(!full(), "ColumnFill: No capacity to add more elements.");
    ptr[column*rows+fillPointer] = a;
		fillPointer += 1;
	}
	
	void opOpAssign(string op, T)(T v) if (op == "~") {
		enforce(this.capacity >= v.length, "ColumnFill: capacity " ~ this.capacity.to!string ~ " not large enough for " ~ v.length.to!string ~ " elements.");
		foreach(ii; 0..v.length) {
			opOpAssign!"~"(v[ii]);
		}
	}
	
	bool full() {
		return this.capacity <= 0;
	}
	
	long capacity() {
		return length - fillPointer;
	}
}
	
struct RowFill {
  long rows;
  long cols;
  double * ptr;
  long column;
  long length;
  private long fillPointer = 0;
  
	this(Matrix m, long row) {
		rows = m.rows;
		cols = m.cols;
		ptr = &m.ptr[row];
		length = cols;
	}
	
	void opOpAssign(string op)(double a) if (op == "~") {
		enforce(!full(), "RowFill: No capacity to add more elements.");
    ptr[rows*fillPointer] = a;
		fillPointer += 1;
	}
	
	void opOpAssign(string op, T)(T v) if (op == "~") {
		enforce(this.capacity >= v.length, "RowFill: capacity " ~ this.capacity.to!string ~ " not large enough for " ~ v.length.to!string ~ " elements.");
		foreach(ii; 0..v.length) {
			opOpAssign!"~"(v[ii]);
		}
	}
	
	bool full() {
		return this.capacity <= 0;
	}
	
	long capacity() {
		return length - fillPointer;
	}
}

struct DiagonalFill {
	long rows;
	long cols;
	double * ptr;
	long length;
	private long fillPointer = 0;
	
	this(Matrix m) {
		enforce(m.rows == m.cols, "Matrix must be diagonal to use DiagonalFill");
		rows = m.rows;
		cols = m.cols;
		ptr = m.ptr;
		length = m.rows;
	}
	
	void opOpAssign(string op)(double a) if (op == "~") {
		enforce(!full(), "DiagonalFill: No capacity to add more elements");
		ptr[fillPointer*rows+fillPointer] = a;
		fillPointer += 1;
	}
	
	void opOpAssign(string op, T)(T v) if (op == "~") {
		enforce(this.capacity >= v.length, "DiagonalFill: capacity " ~ this.capacity.to!string ~ " not large enough for " ~ v.length.to!string ~ " elements.");
		foreach(ii; 0..v.length) {
			opOpAssign!"~"(v[ii]);
		}
	}
	
	bool full() {
		return this.capacity <= 0;
	}
	
	long capacity() {
		return length - fillPointer;
	}
}

struct AboveDiagonalFill {
  long rows;
  long cols;
  double * ptr;
  private Element[] elements;

  this(Matrix m) {
    enforce(m.rows == m.cols, "Matrix must be diagonal to use AboveDiagonalFill");
    rows = m.rows;
    cols = m.cols;
    ptr = m.ptr;
    if (cols > 1) {
			foreach(col; 1..cols) {
				foreach(row; 0..col) {
					elements ~= Element(row, col);
				}
			}
		}
  }

  bool full() {
    return elements.length == 0;
  }

	long capacity() {
		return elements.length;
	}

  void opOpAssign(string op)(double a) if (op == "~") {
    enforce(elements.length > 0, "No more room to add elements");
    ptr[elements[0].column*rows+elements[0].row] = a;
    elements.popFrontN(1);
  }

  void opOpAssign(string op, T)(T v) if (op == "~") {
		enforce(capacity >= v.length, "Not enough room to insert " ~ v.length.to!string ~ " elements. Capacity is only " ~ capacity.to!string ~ ".");
		foreach(ii; 0..v.length) {
			ptr[elements[0].column*rows+elements[0].row] = v[ii];
			elements.popFrontN(1);
		}
  }
}

struct BelowDiagonalFill {
  long rows;
  long cols;
  double * ptr;
  private Element[] elements;

  this(Matrix m) {
    enforce(m.rows == m.cols, "Matrix must be diagonal to use AboveDiagonalFill");
    rows = m.rows;
    cols = m.cols;
    ptr = m.ptr;
		foreach(col; 0..cols-1) {
			foreach(row; col+1..rows) {
				elements ~= Element(row, col);
			}
		}
  }

  bool full() {
    return elements.length == 0;
  }

	long capacity() {
		return elements.length;
	}

  void opOpAssign(string op)(double a) if (op == "~") {
    enforce(elements.length > 0, "No more room to add elements");
    ptr[elements[0].column*rows+elements[0].row] = a;
    elements.popFrontN(1);
  }

  void opOpAssign(string op, T)(T v) if (op == "~") {
		enforce(capacity >= v.length, "Not enough room to insert " ~ v.length.to!string ~ " elements. Capacity is only " ~ capacity.to!string ~ ".");
		foreach(ii; 0..v.length) {
			ptr[elements[0].column*rows+elements[0].row] = v[ii];
			elements.popFrontN(1);
		}
  }
}

struct MatrixFill {
	Submatrix sm;
	private Element[] elements;
	alias sm this;
	
	this(Matrix m) {
		sm = Submatrix(m, [0,0], [m.rows,m.cols]);
		foreach(col; 0..sm.cols) {
			foreach(row; 0..sm.rows) {
				elements ~= Element(row, col);
			}
		}
	}
	
	this(Submatrix _sm) {
		sm.rows = _sm.rows;
		sm.cols = _sm.cols;
		sm.matrows = _sm.matrows;
		sm.offset = _sm.offset;
		sm.ptr = _sm.ptr;
		foreach(col; 0..sm.cols) {
			foreach(row; 0..sm.rows) {
				elements ~= Element(row, col);
			}
		}
	}

  bool full() {
    return elements.length == 0;
  }

	long capacity() {
		return elements.length;
	}

  void opOpAssign(string op)(double a) if (op == "~") {
    enforce(elements.length > 0, "No more room to add elements");
    ptr[elements[0].column*rows+elements[0].row] = a;
    elements.popFrontN(1);
  }

  void opOpAssign(string op, T)(T v) if (op == "~") {
		enforce(capacity >= v.length, "Not enough room to insert " ~ v.length.to!string ~ " elements. Capacity is only " ~ capacity.to!string ~ ".");
		foreach(ii; 0..v.length) {
			ptr[elements[0].column*rows+elements[0].row] = v[ii];
			elements.popFrontN(1);
		}
  }
}

struct BlockFill {
	Submatrix sm;
	private Element[] elements;
	
	this(Matrix m, long[2] element0, long[2] element1) {
		sm = Submatrix(m, element0, element1);
		foreach(col; 0..sm.cols) {
			foreach(row; 0..sm.rows) {
				elements ~= Element(row, col);
			}
		}
	}
	
	this(Submatrix _sm) {
		sm.rows = _sm.rows;
		sm.cols = _sm.cols;
		sm.matrows = _sm.matrows;
		sm.offset = _sm.offset;
		sm.ptr = _sm.ptr;
		foreach(col; 0..sm.cols) {
			foreach(row; 0..sm.rows) {
				elements ~= Element(row, col);
			}
		}
	}

  bool full() {
    return elements.length == 0;
  }

	long capacity() {
		return elements.length;
	}

  void opOpAssign(string op)(double a) if (op == "~") {
    enforce(elements.length > 0, "No more room to add elements");
    ptr[elements[0].column*rows+elements[0].row] = a;
    elements.popFrontN(1);
  }

  void opOpAssign(string op, T)(T v) if (op == "~") {
		enforce(capacity >= v.length, "Not enough room to insert " ~ v.length.to!string ~ " elements. Capacity is only " ~ capacity.to!string ~ ".");
		foreach(ii; 0..v.length) {
			ptr[elements[0].column*rows+elements[0].row] = v[ii];
			elements.popFrontN(1);
		}
  }
}

struct SliceIndex {
	long start;
	long end;
	private long offset = 0;
	
	this(AllElements ae, long _end) {
		start = 0;
		end = _end;
	}
	
	this(long index, long unused) {
		start = index;
		end = index+1;
	}
	
	this(long[2] ind, long unused) {
		start = ind[0];
		end = ind[1];
	}
	
	this(SliceIndex si, long unused) {
		start = si.start;
		end = si.end;
	}
	
	this(long[2] ind) {
		start = ind[0];
		end = ind[1];
	}
	
	long length() {
		return end - start;
	}
	
	bool empty() {
		return (offset >= end - start);
	}

	long front() {
		return start + offset;
	}

	void popFront() {
		offset += 1;
	}
}

struct IntMatrix {
  long rows;
  long cols;
  RData data;
  int * ptr;
  alias data this;

  this(string code) {
    data = RData(code);
    rows = Rf_nrows(data.x);
    cols = Rf_ncols(data.x);
    ptr = INTEGER(data.x);
  }

  this(long r, long c) {
    this("matrix(0L, nrow=" ~ r.to!string ~ ", ncol=" ~ c.to!string ~ ")");
  }

  int opIndex(long r, long c) {
    enforce(r < this.rows, "First index exceeds the number of rows");
    enforce(c < this.cols, "Second index exceeds the number of columns");
    return ptr[c*this.rows+r];
  }

  void opIndexAssign(long v, long r, long c) {
    ptr[c*rows+r] = v.to!int;
  }
}

struct BoolMatrix {
  long rows;
  long cols;
  RData data;
  int * ptr;
  alias data this;

  this(string code) {
    data = RData(code);
    rows = Rf_nrows(data.x);
    cols = Rf_ncols(data.x);
    ptr = INTEGER(data.x);
  }

  this(long r, long c) {
    this("matrix(FALSE, nrow=" ~ r.to!string ~ ", ncol=" ~ c.to!string ~ ")");
  }

  bool opIndex(long r, long c) {
    enforce(r < this.rows, "First index exceeds the number of rows");
    enforce(c < this.cols, "Second index exceeds the number of columns");
    return ptr[c*this.rows+r].to!bool;
  }

  void opIndexAssign(bool v, long r, long c) {
    ptr[c*rows+r] = v.to!int;
  }
}
