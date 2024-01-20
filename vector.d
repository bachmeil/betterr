module betterr.vector;
import betterr.matrix, betterr.rdata;
import betterr.r;
import std.conv, std.exception, std.stdio, std.traits, std.utf;

struct Vector {
  long rows;
  RData data;
  double * ptr;
  alias data this;
  
  this(string code) {
    data = RData(code);
    rows = Rf_length(data.x);
    ptr = REAL(data.x);
  }

  this(long r) {
    this("double(" ~ r.to!string ~ ")");
  }
  
  /* Used to allow assignment to a zero-length Vector. */
  void initialize(long r) {
    rows = r;
    data = RData("double(" ~ r.to!string ~ ")");
    ptr = REAL(data.x);
  }
  
  this(Vector v) {
    this(v.data.name ~ "[]");
  }
  
  this(RData rd) {
		if (Rf_isVector(rd.x)) {
			this(rd.name ~ "[]");
		} else {
			this(rd.name);
		}
  }
  
  this(double[] v) {
    this(v.length);
    foreach(ii; 0..to!int(v.length)) {
      ptr[ii] = v[ii];
    }
  }

  this(long[] v) {
    this(v.length);
    foreach(ii; 0..to!int(v.length)) {
      ptr[ii] = double(v[ii]);
    }
  }

  this(int[] v) {
    this(v.length);
    foreach(ii; 0..to!int(v.length)) {
      ptr[ii] = double(v[ii]);
    }
  }

  double opIndex(long r) {
    enforce(r < rows, "Index out of range: Index on Vector is too large");
    return ptr[r];
  }
  
  Vector opIndex(long[] obs) {
		auto result = Vector(obs.length);
		foreach(ii; 0..obs.length) {
			result[ii] = this[obs[ii]];
		}
		return result;
	}

  void opIndexAssign(double v, long r) {
    enforce(r < rows, "Index out of range: index on RVector is too large");
    ptr[r] = v;
  }

  void opAssign(Vector v) {
    /* This allows the case where you need to assign to a newly constructed
     * Vector for which the length is zero. */
    if (rows == 0) {
      initialize(v.length);
    }
    enforce(v.length == rows, "Vectors have different lengths");
    foreach(ii; 0..v.length) {
      this[ii] = v[ii];
    }
  }
  
  void opAssign(double[] v) {
    if (rows == 0) {
      initialize(v.length);
    }
    enforce(v.length == rows, "Vectors have different lengths");
    foreach(ii; 0..v.length) {
      this[ii] = v[ii];
    }
  }
  
  void opAssign(double a) {
    ptr[0..rows] = a;
  }
  
  Vector opSlice(long i, long j) {
    enforce(j <= rows, "Index out of range: index on RVector slice is too large. index=" ~ to!string(j) ~ " # rows=" ~ to!string(rows));
    enforce(i < j, "First index has to be less than second index");
    return Vector(data.name ~ "[" ~ to!string(i+1) ~ ":" ~ to!string(j) ~ "]");
  }
  
  Vector opSlice() {
    return opSlice(0, rows);
  }
  
  void opSliceAssign(double a) {
		this.ptr[0..this.rows] = a;
	}
	
	void opSliceAssign(double[] v) {
		enforce(v.length == this.rows, "Wrong number of elements");
		this.ptr[0..this.rows] = v[];
	}
	
	void opSliceAssign(double a, long ind0, long ind1) {
		enforce(ind0 >= 0, "First index cannot be negative");
		enforce(ind1 <= this.rows, "Second index out of bounds");
		enforce(ind0 <= ind1, "Second index cannot be less than the first");
		this.ptr[ind0..ind1] = a;
	}
  
	void opSliceAssign(double[] v, long ind0, long ind1) {
		enforce(ind0 >= 0, "First index cannot be negative");
		enforce(ind1 <= this.rows, "Second index out of bounds");
		enforce(ind0 <= ind1, "Second index cannot be less than the first");
		enforce(ind1-ind0 == v.length, "Wrong number of elements");
		this.ptr[ind0..ind1] = v[];
	}
  
  double last() {
		return this[rows-1];
	}
  
  long opDollar() {
		return rows;
	}

  void print(string msg="") {
    if (msg != "") { writeln(msg, ":"); }
    printR(data.x);
  }

  long length() {
    return rows;
  }
	
  bool empty() {
    return rows == 0;
  }

  double front() {
    return this[0];
  }

  void popFront() {
    ptr = &ptr[1];
    rows -= 1;
  }
  
  Matrix matrix() {
		return Matrix("as.matrix(" ~ this.name ~ ")");
	}
	
	Matrix rowMatrix() {
		return Matrix("t(as.matrix(" ~ this.name ~ "))");
	}
  
  double[] opCast(T: double[])() {
    double[] result;
    result.reserve(rows);
    foreach(val; this) {
      result ~= val;
    }
    return result;
  }
  
  auto opBinary(string op, T)(T x) {
    static if (op == "~") {
      return cbind(this, x);
    }
  }
  
  auto opBinary(string op)(Vector v) {
		static if (op == "+") {
			return Vector(this.name ~ " + " ~ v.name);
		}
		
		static if (op == "-") {
			return Vector(this.name ~ " - " ~ v.name);
		}
		
		static if (op == "*") {
			return Vector(this.name ~ " * " ~ v.name);
		}
		
		static if (op == "/") {
			return Vector(this.name ~ " / " ~ v.name);
		}

    static if (op == "~") {
      return cbind(this, v);
    }
	}
}

Vector head(Vector v, long n=6) {
	return Vector(`head(` ~ v.name ~ `, ` ~ n.to!string ~ `)`);
}

Vector tail(Vector v, long n=6) {
	return Vector(`tail(` ~ v.name ~ `, ` ~ n.to!string ~ `)`);
}

/* f should be able to be evaluated at any positive number
 * If it can't be evaluated at 1.0 for some reason, you can
 * specify val. */
auto mapRows(alias f, double val=1.0)(Vector v) {
  alias T = double[];
  static if (is(typeof(f(val)) == double)) {
    auto result = Vector(v.rows);
    foreach(ii; 0..v.rows) {
      result[ii] = f(v[ii]);
    }
    return result;
  }
  static if (is(typeof(f(val)) == long) | is(typeof(f(val)) == int)) {
    auto result = IntVector(v.rows);
    foreach(ii; 0..v.rows) {
      result[ii] = f(v[ii]);
    }
    return result;
  }
  static if (is(typeof(f(val)) == bool)) {
    auto result = BoolVector(v.rows);
    foreach(ii; 0..v.rows) {
      result[ii] = f(v[ii]);
    }
    return result;
  }
  static if (is(typeof(f(val)) == T)) {
    auto result = Matrix(v.rows, f(val).length);
    foreach(ii; 0..v.rows) {
      auto tmp = f(v[ii]);
      foreach(jj; 0..tmp.length) {
        result[ii, jj] = tmp[jj];
      }
    }
    return result;
  }
}

/* This is for filling the elements and nothing else. */
struct Fill {
	double * ptr;
	long length;
	private long fillPointer = 0;
	
	this(Vector v) {
		ptr = v.ptr;
		length = v.length;
	}
	
	/* User code should check that all elements have been filled.
	 * Especially important if reusing the Vector. */
	bool full() {
		return this.capacity <= 0;
	}
	
	long capacity() {
		return length - fillPointer;
	}
	
	void opOpAssign(string op)(double a) if (op == "~") {
		enforce(!full(), "Fill: Cannot append an additional element because it is full.");
		ptr[fillPointer] = a;
		fillPointer += 1;
	}
	
	void opOpAssign(string op, T)(T v) if (op == "~") {
		enforce(this.capacity >= v.length, "Fill: capacity " ~ this.capacity.to!string ~ " not large enough for " ~ v.length.to!string ~ " elements.");
		foreach(ii; 0..v.length) {
			opOpAssign!"~"(v[ii]);
		}
	}
	
	void reset() {
		fillPointer = 0;
	}
}

/* This could be combined with Fill, but this is less error-prone for
 * the user. */
struct PartialFill {
	double * ptr;
	long length;
	private long fillPointer = 0;
	
	/* Following D practices for slices, end is not included. */
	this(Vector v, long start, long end) {
		ptr = &v.ptr[start];
		length = end - start;
	}
	
	/* User code should check that all elements have been filled.
	 * Especially important if reusing the Vector. */
	bool full() {
		return this.capacity <= 0;
	}
	
	long capacity() {
		return length - fillPointer;
	}
	
	void opOpAssign(string op)(double a) if (op == "~") {
		enforce(!full(), "PartialFill: Cannot append an additional element because it is full.");
		ptr[fillPointer] = a;
		fillPointer += 1;
	}
	
	void opOpAssign(string op, T)(T v) if (op == "~") {
		enforce(this.capacity >= v.length, "PartialFill: capacity " ~ this.capacity.to!string ~ " not large enough for " ~ v.length.to!string ~ " elements.");
		foreach(ii; 0..v.length) {
			opOpAssign!"~"(v[ii]);
		}
	}
	
	void reset() {
		fillPointer = 0;
	}
}
			
struct IntVector {
  long rows;
  RData data;
  int * ptr;
  alias data this;
  
  this(string code) {
    data = RData(code);
    rows = Rf_length(data.x);
    ptr = INTEGER(data.x);
  }

  this(long r) {
    this("integer(" ~ r.to!string ~ ")");
  }
  
  /* Used to allow assignment to a zero-length Vector. */
  void initialize(long r) {
    rows = r;
    data = RData("integer(" ~ r.to!string ~ ")");
    ptr = INTEGER(data.x);
  }
  
  this(IntVector v) {
    this(v.data.name ~ "[]");
  }
  
  /* Tempting to try to do this more efficiently, but what happens
   * if the underlying pointer to the data changes, if for instance
   * you change the data? There is no guarantee that the pointer will
   * stay the same. */
  this(RData rd) {
		if (Rf_isVector(rd.x)) {
			this(rd.name ~ "[]");
		} else {
			this(rd.name);
		}
  }
  
  this(long[] v) {
    this(v.length);
    foreach(ii; 0..to!int(v.length)) {
      ptr[ii] = v[ii].to!int;
    }
  }

  this(int[] v) {
    this(v.length);
    foreach(ii; 0..to!int(v.length)) {
      ptr[ii] = double(v[ii]);
    }
  }

  int opIndex(long r) {
    enforce(r < rows, "Index out of range: Index on IntVector is too large");
    return ptr[r];
  }
  
  IntVector opIndex(long[] obs) {
		auto result = IntVector(obs.length);
		foreach(ii; 0..obs.length) {
			result[ii] = this[obs[ii]];
		}
		return result;
	}

  void opIndexAssign(long v, long r) {
    enforce(r < rows, "Index out of range: index on IntVector is too large");
    ptr[r] = v.to!int;
  }

  void opIndexAssign(int v, long r) {
    enforce(r < rows, "Index out of range: index on IntVector is too large");
    ptr[r] = v;
  }

  void opAssign(IntVector v) {
    if (rows == 0) {
      initialize(v.length);
    }
    enforce(v.length == rows, "Vectors have different lengths");
    foreach(ii; 0..v.length) {
      this[ii] = v[ii];
    }
  }
  
  void opAssign(long[] v) {
    if (rows == 0) {
      initialize(v.length);
    }
    enforce(v.length == rows, "Vectors have different lengths");
    foreach(ii; 0..v.length) {
      this[ii] = v[ii].to!int;
    }
  }
  
  void opAssign(int[] v) {
    if (rows == 0) {
      initialize(v.length);
    }
    enforce(v.length == rows, "Vectors have different lengths");
    foreach(ii; 0..v.length) {
      this[ii] = v[ii];
    }
  }
  
  IntVector opSlice(long i, long j) {
    enforce(j <= rows, "Index out of range: index on RVector slice is too large. index=" ~ to!string(j) ~ " # rows=" ~ to!string(rows));
    enforce(i < j, "First index has to be less than second index");
    return IntVector(data.name ~ "[" ~ to!string(i+1) ~ ":" ~ to!string(j) ~ "]");
  }
  
  IntVector opSlice() {
    return opSlice(0, rows);
  }
  
  int last() {
		return this[rows-1];
	}
  
  long opDollar() {
		return rows;
	}
  
  string rcode() {
    string result = "c(";
    foreach(ii; 0..this.length) {
      if (ii > 0) {
        result ~= ", " ~ this[ii].to!string;
      } else {
        result ~= this[ii].to!string;
      }
    }
    return result ~ ")";
  }

  void print(string msg="") {
    if (msg != "") { writeln(msg, ":"); }
    printR(data.x);
  }

  long length() {
    return rows;
  }
	
  bool empty() {
    return rows == 0;
  }

  int front() {
    return this[0];
  }

  void popFront() {
    ptr = &ptr[1];
    rows -= 1;
  }
}

struct BoolVector {
  long rows;
  RData data;
  int * ptr;
  alias data this;
  
  this(string code) {
    data = RData(code);
    rows = Rf_length(data.x);
    ptr = LOGICAL(data.x);
  }

  this(long r) {
    this("logical(" ~ r.to!string ~ ")");
  }
  
  this(BoolVector v) {
    this(v.data.name ~ "[]");
  }
  
  /* Tempting to try to do this more efficiently, but what happens
   * if the underlying pointer to the data changes, if for instance
   * you change the data? There is no guarantee that the pointer will
   * stay the same. */
  this(RData rd) {
		if (Rf_isVector(rd.x)) {
			this(rd.name ~ "[]");
		} else {
			this(rd.name);
		}
  }
  
  this(bool[] v) {
    this(v.length);
    foreach(ii; 0..to!int(v.length)) {
      ptr[ii] = v[ii].to!int;
    }
  }

  this(int[] v) {
    this(v.length);
    foreach(ii; 0..to!int(v.length)) {
      ptr[ii] = v[ii];
    }
  }

  bool opIndex(long r) {
    enforce(r < rows, "Index out of range: Index on IntVector is too large");
    return ptr[r].to!bool;
  }
  
  BoolVector opIndex(long[] obs) {
		auto result = BoolVector(obs.length);
		foreach(ii; 0..obs.length) {
			result[ii] = opIndex(obs[ii]);
		}
		return result;
	}

  void opIndexAssign(bool v, long r) {
    enforce(r < rows, "Index out of range: index on IntVector is too large");
    ptr[r] = v.to!int;
  }

  void opIndexAssign(int v, long r) {
    enforce(r < rows, "Index out of range: index on IntVector is too large");
    ptr[r] = v;
  }

  void opAssign(BoolVector v) {
    enforce(v.length == rows, "Vectors have different lengths");
    foreach(ii; 0..v.length) {
      ptr[ii] = v.ptr[ii];
    }
  }
  
  void opAssign(bool[] v) {
    enforce(v.length == rows, "Vectors have different lengths");
    foreach(ii; 0..v.length) {
      this[ii] = v[ii].to!int;
    }
  }
  
  void opAssign(int[] v) {
    enforce(v.length == rows, "Vectors have different lengths");
    foreach(ii; 0..v.length) {
      this[ii] = v[ii];
    }
  }
  
  void opAssign(long[] v) {
    enforce(v.length == rows, "Vectors have different lengths");
    foreach(ii; 0..v.length) {
      this[ii] = v[ii].to!int;
    }
  }
  
  BoolVector opSlice(long i, long j) {
    enforce(j <= rows, "Index out of range: index on RVector slice is too large. index=" ~ to!string(j) ~ " # rows=" ~ to!string(rows));
    enforce(i < j, "First index has to be less than second index");
    return BoolVector(data.name ~ "[" ~ to!string(i+1) ~ ":" ~ to!string(j) ~ "]");
  }
  
  BoolVector opSlice() {
    return opSlice(0, rows);
  }
  
  bool last() {
		return this[rows-1];
	}
  
  long opDollar() {
		return rows;
	}

  void print(string msg="") {
    if (msg != "") { writeln(msg, ":"); }
    printR(data.x);
  }

  long length() {
    return rows;
  }
	
  bool empty() {
    return rows == 0;
  }

  bool front() {
    return this[0];
  }

  void popFront() {
    ptr = &ptr[1];
    rows -= 1;
  }
}

struct StringVector {
  long rows;
  RData data;
  alias data this;
  
  this(string code) {
    data = RData(code);
    rows = Rf_length(data.x);
  }

  this(long r) {
    this("character(" ~ r.to!string ~ ")");
  }
  
  this(StringVector v) {
    this(v.data.name ~ "[]");
  }
  
  this(RData rd) {
		if (Rf_isVector(rd.x)) {
			this(rd.name ~ "[]");
		} else {
			this(rd.name);
		}
  }
  
  this(string[] v) {
    this(v.length);
    foreach(ii; 0..to!int(v.length)) {
      SET_STRING_ELT(data.x, ii, Rf_mkChar(toUTFz!(char*)(v[ii])));
    }
  }

  string opIndex(long r) {
    enforce(r < rows, "Index out of range: Index on StringVector is too large");
    return to!string(R_CHAR(STRING_ELT(data.x, r.to!int)));
  }
  
  StringVector opIndex(long[] obs) {
		auto result = StringVector(obs.length);
		foreach(ii; 0..obs.length) {
      SET_STRING_ELT(result.x, ii.to!int, STRING_ELT(data.x, ii.to!int));
		}
		return result;
	}

  void opIndexAssign(string s, long r) {
    enforce(r < rows, "Index out of range: index on StringVector is too large");
    SET_STRING_ELT(data.x, r.to!int, Rf_mkChar(toUTFz!(char*)(s)));
  }

  void opAssign(StringVector v) {
    enforce(v.length == rows, "StringVectors have different lengths");
    foreach(ii; 0..v.length) {
      SET_STRING_ELT(data.x, ii.to!int, (STRING_ELT(v.x, ii.to!int)));
    }
  }
  
  void opAssign(string[] v) {
    enforce(v.length == rows, "Vectors have different lengths");
    foreach(ii; 0..v.length) {
      this[ii] = v[ii];
    }
  }
  
  void opAssign(string s) {
    foreach(ii; 0..this.length) {
      this[ii] = s;
    }
  }

  StringVector opSlice(long i, long j) {
    enforce(j <= rows, "Index out of range: index on StringVector slice is too large. index=" ~ to!string(j) ~ " # rows=" ~ to!string(rows));
    enforce(i < j, "First index has to be less than second index");
    return StringVector(data.name ~ "[" ~ to!string(i+1) ~ ":" ~ to!string(j) ~ "]");
  }
  
  StringVector opSlice() {
    return opSlice(0, rows);
  }
  
  string[] opCast(T: string[])() {
    string[] result;
    foreach(ii; 0..this.length) {
      result[ii] = this[ii];
    }
    return result;
  }
  
  string last() {
		return this[rows-1];
	}
  
  long opDollar() {
		return rows;
	}

  void print(string msg="") {
    if (msg != "") { writeln(msg, ":"); }
    printR(data.x);
  }

  long length() {
    return rows;
  }
}
