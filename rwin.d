module betterr.rwin;

import std.algorithm, std.array, std.conv, std.datetime, std.math;
import std.meta, std.range, std.stdio, std.string, std.sumtype;
import std.traits, std.utf;
import std.exception;

struct AllElements {};
public AllElements _all;

version(Windows) {
  alias vvptr = extern(C) void function();
  vvptr startR;
  vvptr closeR;

  alias vrobjptr = extern(C) void function(Robj);
  vrobjptr Rf_unprotect_ptr;
  vrobjptr Rf_PrintValue;

  alias robjintptr = extern(C) Robj function(int);
  robjintptr Rf_ScalarInteger;

  alias v_r_cp = extern(C) void function(Robj, char *);
  v_r_cp passToR;

  alias r_cp = extern(C) Robj function(char *);
  r_cp evalInR;

  alias v_cp = extern(C) void function(char *);
  v_cp evalQuietlyInR;

  // rpath: Path of R.dll. Use forward slashes rather than backslashes
  // rinsidepath: Path of libRInside.dll. Use forward slashes.
  //   Easiest to find this directory with `find.package("RInside")` in R.
  // Set errors to true if you want to see error messages for the
  // loading.
  void loadRFunctions(string rpath, string rinsidepath, bool errors=false) {
    import std.process, std.utf;
    import core.sys.windows.windows;
    import std.windows.syserror;

    // libRInside.dll relies on R.dll, and this is how it knows where to look
    environment["PATH"] = rpath;

    // Error 126 means *something* can't be found, including a dependency.
    HMODULE rinside = LoadLibraryW(toUTF16z(rinsidepath ~ "/libRInside.dll"));
    if (errors) {
      writeln(sysErrorString(GetLastError()));
    }
    
    // These are already defined in rwin.d. If you declare them auto or some
    // other type, you won't be assigning the value to the right function pointer.
    // It will be local to this function, and thus useless.
    startR = cast(vvptr) GetProcAddress(rinside, "setupRinC");
    closeR = cast(vvptr) GetProcAddress(rinside, "teardownRinC");
    passToR = cast(v_r_cp) GetProcAddress(rinside, "passToR");
    evalInR = cast(r_cp) GetProcAddress(rinside, "evalInR");
    evalQuietlyInR = cast(v_cp) GetProcAddress(rinside, "evalQuietlyInR");

    HMODULE rdll = LoadLibraryW(toUTF16z(rpath ~ "/R.dll"));
    if (errors) {
      writeln(sysErrorString(GetLastError()));
    }

    Rf_unprotect_ptr = cast(vrobjptr) GetProcAddress(rdll, "Rf_unprotect_ptr");
    Rf_PrintValue = cast(vrobjptr) GetProcAddress(rdll, "Rf_PrintValue");
    Rf_ScalarInteger = cast(robjintptr) GetProcAddress(rdll, "Rf_ScalarInteger");
  }
}

// version(gretl) {
//   import gretl.matrix;
// }

struct sexprec {}
alias Robj = sexprec*;

// This enables reference counting of objects allocated by R
// Handles unprotecting for you
// DOES NOT DO THE *PROTECTING* OF AN R OBJECT
// It only stores a protected object and unprotects it when there are no further references to it
// You need to create the ProtectedRObject when you allocate
// It is assumed that you will not touch the .robj directly
// unprotect is needed because only some Robj will need to be unprotected
struct RObjectStorage {
  Robj ptr;
  bool unprotect;
  int refcount;
}

struct ProtectedRObject {
  RObjectStorage * data;
  alias data this;
	
  // x should already be protected
  // ProtectedRObject is for holding an Robj, not for allocating it
  this(Robj x, bool u=false) {
    data = new RObjectStorage();
    data.ptr = x;
    data.refcount = 1;
    data.unprotect = u;
  }
	
  this(this) {
    if (data !is null) {
      if (data.unprotect) {
        enforce(data !is null, "data should never be null inside an ProtectedRObject. You must have created an ProtectedRObject without using the constructor.");
        data.refcount += 1;
      }
    }
  }
  
  /* I don't entirely understand this. Sometimes data is null if unprotect is
   * false. Maybe the compiler knows unprotect is false so it'll never
   * need to execute this part of the code? */
  ~this() {
    if (data !is null) {
      if (data.unprotect) {
        enforce(data !is null, "Calling the destructor on a ProtectedRObject when data is null. You must have created an ProtectedRObject without using the constructor.");
        data.refcount -= 1;
        if (data.refcount == 0) {
          //~ writeln("unprotecting");
          Rf_unprotect_ptr(data.ptr);
        }
      }
    } else {
      writeln(&this, " data is null");
    }
  }
	
  Robj robj() {
    return data.ptr;
  }
}

extern (C) {
	//void passToR(Robj x, char * name);
	//Robj evalInR(char * cmd);
  // This needs to be @nogc so I can call it from the destructor
	//@nogc void evalQuietlyInR(char * cmd);
	void setupRinC();
	void teardownRinC();
}

void toR(T)(T x, string name) {
	passToR(x.robj, toUTFz!(char*)(name));
}

void toR(Robj x, string name) {
	passToR(x, toUTFz!(char*)(name));
}

// void toR(string[] s, string name) {
// 	passToR(s.robj, toUTFz!(char*)(name));
// }

Robj evalR(string cmd) {
  writeln("inside evalR");
  writeln("evalInR: ", evalInR);
	return evalInR(toUTFz!(char*)(cmd));
}

void evalRQ(string cmd) {
	evalQuietlyInR(toUTFz!(char*)(cmd));
}

void evalRQ(string[] cmds) {
	foreach(cmd; cmds) {
		evalQuietlyInR(toUTFz!(char*)(cmd));
	}
}

// void assertR(bool test, string msg) {
//   if (!test) {
//     Rf_error( toUTFz!(char*)("Error in D code: " ~ msg) );
//   }
// }

void printR(Robj x) {
  Rf_PrintValue(x);
}

// void printR(ProtectedRObject x) {
//   Rf_PrintValue(x.robj);
// }

// void printR(string s) {
// 	evalRQ(`print(` ~ s ~ `)`);
// }

// void source(string s) {
// 	evalRQ(`source("` ~ s ~ `")`);
// }

// int length(Robj x) {
//   return Rf_length(x);
// }

// bool isVector(Robj x) {
//   return to!bool(Rf_isVector(x));
// }

// bool isMatrix(Robj x) {
//   return to!bool(Rf_isMatrix(x));
// }

// bool isNumeric(Robj x) {
//   return to!bool(Rf_isNumeric(x));
// }

// bool isInteger(Robj x) {
//   return to!bool(Rf_isInteger(x));
// }

// /* In case you have to pass an R symbol name to a C function */
// Robj RSymbol(string sym) {
//   return Rf_install(toUTFz!(char*)(sym));
// }

// Robj RTrue() {
//   return Rf_ScalarLogical(1);
// }

// Robj RFalse() {
//   return Rf_ScalarLogical(0);
// }

// // RList is for passing data from D to R in a list
// // It's the only way to pass multiple values back to R
// struct RList {
//   ProtectedRObject data;
//   int length; // Length of the underlying Robj, which can never change
//   string[] names;
//   int fillPointer = 0;
//   private int counter = 0; // Used for foreach

//   this(int n) {
//     Robj temp;
//     Rf_protect(temp = Rf_allocVector(19, n));
//     data = ProtectedRObject(temp, true);
//     length = n;
//     names = new string[n];
//   }

//   // For an existing list - by default, assumes the list is already protected
//   // This list is full by construction
//   this(Robj v, bool u=false) {
//     enforce(to!bool(Rf_isVectorList(v)), "Cannot pass a non-list to the constructor for an RList");
//     data = ProtectedRObject(v, u);
//     length = v.length;
//     names = v.names;
//     fillPointer = v.length;
//   }
	
//   version(standalone) {
//     this(string name) {
//       this(evalR(name));
//     }
//   }
	
//   Robj opIndex(int ii) {
//     enforce(ii < length, "RList index has to be less than the number of elements");
//     return VECTOR_ELT(data.robj, ii);
//   }
  
// 	Robj opIndex(string name) {
//     auto ind = countUntil!"a == b"(names, name);
//     if (ind == -1) { enforce(false, "No element in the list with the name " ~ name); }
//     return opIndex(ind.to!int);
// 	}
  
//   void unsafePut(Robj x, int ii) {
//     enforce(ii < length, "RList index has to be less than the number of elements. Index " ~ to!string(ii) ~ " >= length " ~ to!string(length));
//     SET_VECTOR_ELT(data.robj, ii, x);
//   }
	
//   void put(Robj x, string name) {
//     enforce(fillPointer < length, "RList is full - cannot add more elements");
//     SET_VECTOR_ELT(data.robj, fillPointer, x);
//     names[fillPointer] = name;
//     fillPointer += 1;
//   }
	
//   void put(Robj x) {
//     put(x, "");
//   }
	
//   void opIndexAssign(Robj x, string name) {
//     put(x, name);
//   }

//   // No opIndexAssign(Robj, int): Use unsafePut instead
//   // Not clear what is going on if we allow rl[3] = x notation
//   // Should not usually want to put an element into a specific index
  
//   void opIndexAssign(RMatrix rm, string name) {
//     put(rm.robj, name);
//   }

//   void opIndexAssign(RVector rv, string name) {
//     put(rv.robj, name);
//   }
  
//   void opIndexAssign(RString rs, string name) {
//     put(rs.robj, name);
//   }
  
//   void opIndexAssign(string s, string name) {
//     put(s.robj, name);
//   }
  
//   void opIndexAssign(string[] sv, string name) {
//     put(RStringArray(sv).robj, name);
//   }
  
//   void opIndexAssign(double v, string name) {
//     put(v.robj, name);
//   }
	
//   void opIndexAssign(double[] vec, string name) {
//     put(vec.robj, name);
//   }
	
//   void opIndexAssign(int v, string name) {
//     put(v.robj, name);
//   }
  
//   bool empty() {
//     return counter == length;
//   }

//   Robj front() {
//     return this[counter];
//   }

//   void popFront() {
//     counter += 1;
//   }
  
//   Robj robj() {
//     setAttrib(data.robj, "names", names.robj);
//     return data.robj;
//   }
// }

// string toString(Robj cstr) {
//   return to!string(R_CHAR(cstr));
// }

// string toString(Robj sv, int ii) {
//   return to!string(R_CHAR(STRING_ELT(sv, ii)));
// }

// string[] stringArray(Robj sv) {
//   string[] result;
//   foreach(ii; 0..Rf_length(sv)) {
//     result ~= toString(sv, ii);
//   }
//   return result;
// }

// string[] stringArray(string name) {
// 	Robj sv = evalR(name);
// 	string[] result;
// 	foreach(ii; 0..Rf_length(sv)) {
// 		result ~= toString(sv, ii);
// 	}
// 	return result;
// }

// struct RString {
//   ProtectedRObject data;
  
//   this(string str) {
//     Robj temp;
//     Rf_protect(temp = Rf_allocVector(16, 1));
//     data = ProtectedRObject(temp, true);
//     SET_STRING_ELT(data.ptr, 0, Rf_mkChar(toUTFz!(char*)(str)));
//   }

//   Robj robj() {
//     return data.ptr;
//   }
// }

// Robj getAttrib(Robj x, string attr) {
//   return Rf_getAttrib(x, RString(attr).robj);
// }

// Robj getAttrib(ProtectedRObject x, string attr) {
//   return Rf_getAttrib(x.ptr, RString(attr).robj);
// }

// Robj getAttrib(Robj x, RString attr) {
//   return Rf_getAttrib(x, attr.robj);
// }

// Robj getAttrib(ProtectedRObject x, RString attr) {
//   return Rf_getAttrib(x.ptr, attr.robj);
// }

// string[] names(Robj x) {
//   return stringArray(getAttrib(x, "names"));
// }

// void setAttrib(Robj x, string attr, ProtectedRObject val) {
//   Rf_setAttrib(x, RString(attr).robj, val.robj);
// }

// void setAttrib(Robj x, RString attr, ProtectedRObject val) {
//   Rf_setAttrib(x, attr.robj, val.robj);
// }

// void setAttrib(Robj x, string attr, Robj val) {
//   Rf_setAttrib(x, RString(attr).robj, val);
// }

// void setAttrib(Robj x, RString attr, Robj val) {
//   Rf_setAttrib(x, attr.robj, val);
// }

// Robj robj(double x) {
//   return Rf_ScalarReal(x);
// }

// // Copies
// Robj robj(double[] v) {
//   return RVector(v).robj;
// }

Robj robj(int x) {
  return Rf_ScalarInteger(x);
}

// Robj robj(string s) {
//   return RString(s).robj;
// }

// Robj robj(string[] sv) {
//   return RStringArray(sv).robj;
// }

// // Copies
// double[] array(Robj rv) {
//   enforce(isVector(rv), "In array(Robj rv): Cannot convert non-vector R object to double[]");
//   enforce(isNumeric(rv), "In array(Robj rv): Cannot convert non-numeric R object to double[]");
//   double[] result;
//   result.reserve(rv.length);
//   double * ptr = REAL(rv);
//   foreach(ii; 0..rv.length.to!int) {
// 		result ~= ptr[ii];
// 	}
// 	return result;
// }

// ProtectedRObject RStringArray(string[] sv) {
//   Robj temp;
//   Rf_protect(temp = Rf_allocVector(16, to!int(sv.length)));
//   auto result = ProtectedRObject(temp, true);
//   foreach(ii; 0..to!int(sv.length)) {
//     SET_STRING_ELT(result.robj, ii, Rf_mkChar(toUTFz!(char*)(sv[ii])));
//   }
//   return result;
// }

// ulong[3] tsp(Robj rv) {
//   auto tsprop = RVector(getAttrib(rv, "tsp"));
//   ulong[3] result;
//   result[0] = lround(tsprop[0]*tsprop[2])+1;
//   result[1] = lround(tsprop[1]*tsprop[2])+1;
//   result[2] = lround(tsprop[2]);
//   return result;
// }

// double scalar(Robj rx) {
//   return Rf_asReal(rx); 
// }

// double scalar(T: double)(Robj rx) {
//   return Rf_asReal(rx);
// }

// int scalar(T: int)(Robj rx) { 
//   return Rf_asInteger(rx); 
// }

// long scalar(T: long)(Robj rx) { 
//   return to!long(rx.scalar!int); 
// }

// ulong scalar(T: ulong)(Robj rx) { 
//   return to!ulong(rx.scalar!int); 
// }

// string scalar(T: string)(Robj rx) { 
//   return to!string(R_CHAR(STRING_ELT(rx,0))); 
// }

// double scalar(string name) {
// 	return Rf_asReal(evalR(name)); 
// }

// double scalar(T: double)(string name) {
//   return Rf_asReal(evalR(name)); 
// }

// int scalar(T: int)(string name) { 
//   return Rf_asInteger(evalR(name)); 
// }

// long scalar(T: long)(string name) { 
//   return to!long(evalR(name).scalar!int); 
// }

// ulong scalar(T: ulong)(string name) { 
//   return to!ulong(evalR(name).scalar!int); 
// }

// string scalar(T: string)(string name) { 
//   return to!string(R_CHAR(STRING_ELT(evalR(name),0))); 
// }

// struct RMatrix {
//   ProtectedRObject data;
//   int rows;
//   int cols;
//   double * ptr;
  
//   this(int r, int c) {
//     Robj temp;
//     Rf_protect(temp = Rf_allocMatrix(14, r, c));
//     data = ProtectedRObject(temp, true);
//     ptr = REAL(temp);
//     rows = r;
//     cols = c;
//   }
  
//   version(gretl) {
//     this(T)(T m) if (is(T == DoubleMatrix) || is(T == GretlMatrix)) {
//       Robj temp;
//       Rf_protect(temp = Rf_allocMatrix(14, m.rows, m.cols));
//       data = ProtectedRObject(temp, true);
//       ptr = REAL(temp);
//       rows = m.rows;
//       cols = m.cols;
//       ptr[0..m.rows*m.cols] = m.ptr[0..m.rows*m.cols];
//     }
//   }

//   version(gretl) {
//     GretlMatrix mat() {
//       GretlMatrix result;
//       result.rows = this.rows;
//       result.cols = this.cols;
//       result.ptr = this.ptr;
//       return result;
//     }
  	
//     alias mat this;
//   }

//   /* Normally this will be a matrix allocated inside R, and as such, it will already be protected. Nonetheless you have the option to protect by setting the second argument to false. */
//   this(Robj rm, bool u=false) {
//     enforce(isMatrix(rm), "Constructing RMatrix from something not a matrix"); 
//     enforce(isNumeric(rm), "Constructing RMatrix from something that is not numeric");
//     data = ProtectedRObject(rm, u);
//     ptr = REAL(rm);
//     rows = Rf_nrows(rm);
//     cols = Rf_ncols(rm);
//   }
  
//   version(standalone) {
//     this(string name) {
//       this(evalR(name));
//     }
//   }
	
//   // Use this only with objects that don't need protection
//   // For "normal" use that's not an issue
//   this(ProtectedRObject rm) {
//     this(rm.ptr);
//   }
	
//   this(RVector v) {
//     data = v.data;
//     rows = v.rows;
//     cols = 1;
//     ptr = v.ptr;
//   }

//   double opIndex(int r, int c) {
//     enforce(r < this.rows, "First index exceeds the number of rows");
//     enforce(c < this.cols, "Second index exceeds the number of columns");
//     return ptr[c*this.rows+r];
//   }

//   void opIndexAssign(double v, int r, int c) {
//     ptr[c*rows+r] = v;
//   }

//   void opAssign(double val) {
//     ptr[0..this.rows*this.cols] = val;
//   }
  
//   void opAssign(RMatrix m) {
//     Robj temp;
//     Rf_protect(temp = Rf_allocMatrix(14, m.rows, m.cols));
//     data = ProtectedRObject(temp, true);
//     ptr = REAL(temp);
//     rows = m.rows;
//     cols = m.cols;
//     ptr[0..m.rows*m.cols] = m.ptr[0..m.rows*m.cols];
//   }

//   version(gretl) {
//     void opAssign(T)(T m) if (is(T == DoubleMatrix) || is(T == GretlMatrix)) {
//       enforce(rows == m.rows, "Number of rows in source (" ~ to!string(m.rows) ~ ") is different from number of rows in destination (" ~ rows ~ ").");
//       enforce(cols == m.cols, "Number of columns in source (" ~ to!string(m.rows) ~ ") is different from number of columns in destination (" ~ rows ~ ").");
//       ptr[0..m.rows*m.cols] = m.ptr[0..m.rows*m.cols];
//     }
//   }

//   RMatrix opBinary(string op)(double a) {
//     static if(op == "+") {
//       return matrixAddition(this, a);
//     }
//     static if(op == "-") {
//       return matrixSubtraction(this, a);
//     }
//     static if(op == "*") {
//       return matrixMultiplication(this, a);
//     }
//     static if(op == "/") {
//       return matrixDivision(this, a);
//     }
//   }

//   RMatrix opBinaryRight(string op)(double a) {
//     static if(op == "+") {
//       return matrixAddition(this, a);
//     }
//     static if(op == "-") {
//       return matrixSubtraction(a, this);
//     }
//     static if(op == "*") {
//       return matrixMultiplication(this, a);
//     }
//     static if(op == "/") {
//       return matrixDivision(a, this);
//     }
//   }

// 	RMatrix matrixAddition(RMatrix m, double a) {
// 		auto result = RMatrix(m.rows, m.cols);
// 		foreach(ii; 0..m.rows*m.cols) {
// 			result.ptr[ii] = m.ptr[ii] + a;
// 		}
// 		return result;
// 	}

// 	RMatrix matrixSubtraction(RMatrix m, double a) {
// 		auto result = RMatrix(m.rows, m.cols);
// 		foreach(ii; 0..m.rows*m.cols) {
// 			result.ptr[ii] = m.ptr[ii] - a;
// 		}
// 		return result;
// 	}

// 	RMatrix matrixSubtraction(double a, RMatrix m) {
// 		auto result = RMatrix(m.rows, m.cols);
// 		foreach(ii; 0..m.rows*m.cols) {
// 			result.ptr[ii] = a - m.ptr[ii];
// 		}
// 		return result;
// 	}

// 	RMatrix matrixMultiplication(RMatrix m, double a) {
// 		auto result = RMatrix(m.rows, m.cols);
// 		foreach(ii; 0..m.rows*m.cols) {
// 			result.ptr[ii] = a*m.ptr[ii];
// 		}
// 		return result;
// 	}

// 	RMatrix matrixDivision(RMatrix m, double a) {
// 		auto result = RMatrix(m.rows, m.cols);
// 		foreach(ii; 0..m.rows*m.cols) {
// 			result.ptr[ii] = m.ptr[ii]/a;
// 		}
// 		return result;
// 	}

// 	RMatrix matrixDivision(double a, RMatrix m) {
// 		auto result = RMatrix(m.rows, m.cols);
// 		foreach(ii; 0..m.rows*m.cols) {
// 			result.ptr[ii] = a/m.ptr[ii];
// 		}
// 		return result;
// 	}
  
//   RMatrix sameSize() {
//     return RMatrix(rows, cols);
//   }

//   Robj robj() {
//     return data.robj;
//   }
// }

// void print(RMatrix m, string msg="") {
//   writeln(msg);
//   foreach(row; 0..m.rows) {
//     foreach(col; 0..m.cols) {
//       write(m[row,col], " ");
//     }
//     writeln("");
//   }
// }

// // Copies
// RMatrix dup(RMatrix rm) { 
//   RMatrix result = RMatrix(Rf_protect(Rf_duplicate(rm.robj)), true);
//   return result;
// }

// struct MatrixIndex {
//   int rows;
//   int cols;
//   int currentRow=0;
//   int currentCol=0;
  
//   this(RMatrix rm) {
//     rows = rm.rows;
//     cols = rm.cols;
//   }

//   bool empty() {
//     return currentCol >= cols;
//   }

//   int[2] front() {
//     return [currentRow, currentCol];
//   }

//   void popFront() {
//     if (currentRow >= rows-1) {
//       currentCol += 1;
//       currentRow = 0;
//     } else {
//       currentRow += 1;
//     }
//   }
// }  

// struct TransposeIndex {
//   int rows;
//   int cols;
//   int currentRow=0;
//   int currentCol=0;
  
//   this(RMatrix rm) {
//     rows = rm.rows;
//     cols = rm.cols;
//   }

//   bool empty() {
//     return currentCol >= cols;
//   }

//   int[2] front() {
//     return [currentCol, currentRow];
//   }

//   void popFront() {
//     if (currentRow >= rows-1) {
//       currentCol += 1;
//       currentRow = 0;
//     } else {
//       currentRow += 1;
//     }
//   }
// }  

// struct DiagonalIndex {
//   int rows;
//   int cols;
//   int currentRow=0;
//   int currentCol=0;
  
//   this(RMatrix rm) {
//     rows = rm.rows;
//     cols = rm.cols;
//   }

//   bool empty() {
//     return (currentCol >= cols) | (currentRow >= rows);
//   }

//   int[2] front() {
//     return [currentRow, currentCol];
//   }

//   void popFront() {
//     currentCol += 1;
//     currentRow += 1;
//   }
// }  

// struct BelowDiagonalIndex {
//   int rows;
//   int cols;
//   int currentRow=1;
//   int currentCol=0;
  
//   this(RMatrix rm) {
//     rows = rm.rows;
//     cols = rm.cols;
//   }

//   bool empty() {
//     return currentCol >= cols;
//   }

//   int[2] front() {
//     return [currentRow, currentCol];
//   }

//   void popFront() {
//     if (currentRow >= rows) {
//       currentCol += 1;
//       currentRow = currentCol+1;
//     } else {
//       currentRow += 1;
//     }
//   }
// }  

// struct AboveDiagonalIndex {
//   int rows;
//   int cols;
//   int currentRow=0;
//   int currentCol=1;
  
//   this(RMatrix rm) {
//     rows = rm.rows;
//     cols = rm.cols;
//   }

//   bool empty() {
//     return currentCol >= cols;
//   }

//   int[2] front() {
//     return [currentRow, currentCol];
//   }

//   void popFront() {
//     if (currentRow >= currentCol-1) {
//       currentCol += 1;
//       currentRow = 0;
//     } else {
//       currentRow += 1;
//     }
//   }
// }  

// struct RVector {
//   int rows;
//   double * ptr;
//   ProtectedRObject data;
  
//   version(gretl) {
//     GretlMatrix mat() {
//       GretlMatrix result;
//       result.rows = this.rows;
//       result.cols = 1;
//       result.ptr = this.ptr;
//       return result;
//     }
		
//     alias mat this;
//   }
  
//   this(int r) {
//     Robj temp;
//     Rf_protect(temp = Rf_allocVector(14,r));
//     data = ProtectedRObject(temp, true);
//     rows = r;
//     ptr = REAL(temp);
//   }

//   this(Robj rv, bool u=false) {
//     enforce(isVector(rv), "In RVector constructor: Cannot convert non-vector R object to RVector");
//     enforce(isNumeric(rv), "In RVector constructor: Cannot convert non-numeric R object to RVector");
//     data = ProtectedRObject(rv, u);
//     rows = rv.length;
//     ptr = REAL(rv);
//   }
  
// 	this(string name) {
// 		this(evalR(name));
// 	}
  
//   this(ProtectedRObject rv, bool u=false) {
//     this(rv.robj, u);
//   }

//   this(T)(T v) {
//     Robj temp;
//     Rf_protect(temp = Rf_allocVector(14, to!int(v.length)));
//     data = ProtectedRObject(temp, true);
//     rows = to!int(v.length);
//     ptr = REAL(temp);
//     foreach(ii; 0..to!int(v.length)) {
//       ptr[ii] = v[ii];
//     }
//   }

//   double opIndex(int r) {
//     enforce(r < rows, "Index out of range: index on RVector is too large");
//     return ptr[r];
//   }
  
//   RVector opIndex(int[] obs) {
// 		auto result = RVector(to!int(obs.length));
// 		foreach(ii; 0..to!int(obs.length)) {
// 			result[ii] = this[obs[ii]];
// 		}
// 		return result;
// 	}

//   void opIndexAssign(double v, int r) {
//     enforce(r < rows, "Index out of range: index on RVector is too large");
//     ptr[r] = v;
//   }

//   void opAssign(T)(T x) {
//     enforce(x.length == rows, "Cannot assign to RVector from an object of a different length");
//     foreach(ii; 0..to!int(x.length)) {
//       this[ii] = x[ii];
//     }
//   }
  
//   RVector opSlice(int i, int j) {
//     enforce(j <= rows, "Index out of range: index on RVector slice is too large. index=" ~ to!string(j) ~ " # rows=" ~ to!string(rows));
//     enforce(i < j, "First index has to be less than second index");
//     RVector result = this;
//     result.rows = j-i;
//     result.ptr = &ptr[i];
//     result.data = data;
//     return result;
//   }

//   void print(string msg="") {
//     if (msg != "") { writeln(msg, ":"); }
//     foreach(val; this) {
//       writeln(val);
//     }
//   }

//   int length() {
//     return rows;
//   }
	
//   bool empty() {
//     return rows == 0;
//   }

//   double front() {
//     return this[0];
//   }

//   void popFront() {
//     ptr = &ptr[1];
//     rows -= 1;
//   }

//   double[] array() {
//     double[] result;
//     result.reserve(rows);
//     foreach(val; this) {
//       result ~= val;
//     }
//     return result;
//   }

//   Robj robj() {
//     return data.robj;
//   }
// }

// double fromLast(RVector rv, int ii) {
// 	return rv[rv.length-ii-1];
// }

// double last(RVector rv) {
// 	return rv[rv.length-1];
// }

// struct RIntVector {
//   ProtectedRObject data;
//   ulong length;
//   int * ptr;

//   this(int r) {
//     Robj temp;
//     Rf_protect(temp = Rf_allocVector(13, r));
//     data = ProtectedRObject(temp, true);
//     length = r;
//     ptr = INTEGER(temp);
//   }

//   this(int[] v) {
//     Robj temp;
//     Rf_protect(temp = Rf_allocVector(13, to!int(v.length)));
//     data = ProtectedRObject(temp);
//     length = v.length;
//     ptr = INTEGER(temp);
//     foreach(ii, val; v) {
//       this[ii.to!int] = val;
//     }
//   }

//   this(Robj rv, bool u=false) {
//     enforce(isVector(rv), "In RVector constructor: Cannot convert non-vector R object to RVector");
//     enforce(isInteger(rv), "In RVector constructor: Cannot convert non-integer R object to RVector");
//     data = ProtectedRObject(rv);
//     length = rv.length;
//     ptr = INTEGER(rv);
//   }

//   int opIndex(int obs) {
//     enforce(obs < length, "Index out of range: index on RIntVector is too large");
//     return ptr[obs];
//   }

//   void opIndexAssign(int val, int obs) {
//     enforce(obs < length, "Index out of range: index on RIntVector is too large");
//     ptr[obs] = val;
//   }

//   void opAssign(int[] v) {
//     foreach(ii, val; v) {
//       this[ii.to!int] = val;
//     }
//   }

//   RIntVector opSlice(int i, int j) {
//     enforce(j < length, "Index out of range: index on RIntVector slice is too large");
//     enforce(i < j, "First index on RIntVector slice has to be less than the second index");
//     RIntVector result;
//     result.data = data;
//     result.length = j-i;
//     result.ptr = &ptr[i];
//     return result;
//   }

//   int[] array() {
//     int[] result;
//     result.reserve(length);
//     foreach(val; this) {
//       result ~= val;
//     }
//     return result;
//   }

//   void print() {
//     foreach(val; this) {
//       writeln(val);
//     }
//   }

//   bool empty() {
//     return length == 0;
//   }

//   int front() {
//     return this[0];
//   }

//   void popFront() {
//     ptr = &ptr[1];
//     length -= 1;
//   }
  
//   Robj robj() {
//     return data.robj;
//   }
// }

// /* This is only for working with arrays that were allocated inside R. */
// struct RArrayInsideR {
// 	Robj data; // For now, this needs to be protected/freeing is up to the user - call .free
// 	double * ptr;
// 	private int k = 3; // Number of dimensions, for now only 3 is supported
// 	alias data this;
	
// 	this(Robj x) {
// 		enforce(Rf_isArray(x), "You have passed a non-array to the constructor of RArrayRO.");
// 		data = x;
// 		ptr = REAL(x);
// 	}
	
// 	void opAssign(Robj x) {
// 		enforce(Rf_isArray(x), "You have passed a non-array to the constructor of RArrayRO.");
// 		data = x;
// 		ptr = REAL(x);
// 	}
	
// 	// Can add a similar function for dimNames if desired
// 	int[] dims() {
// 		Robj d = Rf_getAttrib(data, R_DimSymbol);
// 		int len = d.length;
// 		return INTEGER(d)[0..len];
// 	}
	
// 	// Convert a 3D array into an array of RMatrix structs, each holding
// 	// a slice of the array
// 	RMatrix[] unpack() {
// 		RMatrix[] result;
// 		int[] d = this.dims();
// 		int nr = d[0]; // Number of rows in each matrix
// 		int nc = d[1]; // Number of columns in each matrix
// 		int nm = d[2]; // Number of matrices in the 3D array
// 		foreach(ii; 0..nm) {
// 			result ~= slice(nr, nc, ii);
// 		}
// 		return result;
// 	}
	
// 	// For now, only supports a slice on the first two dimensions
// 	// Copies the data into an array
// 	// In principle, this does not need to be copied, but I'm not going to
// 	// dig into it now - that requires a lot of care to get right
// 	RMatrix slice(int nr, int nc, int k) {
// 		auto result = RMatrix(nr, nc);
// 		int a1 = k*nr*nc;
// 		int a2 = (k+1)*nr*nc;
// 		result.ptr[0..nr*nc] = ptr[a1..a2];
// 		return result;
// 	}

// 	void free() {
// 		if (data !is null) {
// 			Rf_unprotect_ptr(data);
// 		}
// 	}
// }

// // Constants pulled from the R API, for compatibility
// immutable double M_E=2.718281828459045235360287471353;
// immutable double M_LOG2E=1.442695040888963407359924681002;
// immutable double M_LOG10E=0.434294481903251827651128918917;
// immutable double M_LN2=0.693147180559945309417232121458;
// immutable double M_LN10=2.302585092994045684017991454684; 
// immutable double M_PI=3.141592653589793238462643383280;
// immutable double M_2PI=6.283185307179586476925286766559; 
// immutable double M_PI_2=1.570796326794896619231321691640;
// immutable double M_PI_4=0.785398163397448309615660845820;
// immutable double M_1_PI=0.318309886183790671537767526745;
// immutable double M_2_PI=0.636619772367581343075535053490;
// immutable double M_2_SQRTPI=1.128379167095512573896158903122;
// immutable double M_SQRT2=1.414213562373095048801688724210;
// immutable double M_SQRT1_2=0.707106781186547524400844362105;
// immutable double M_SQRT_3=1.732050807568877293527446341506;
// immutable double M_SQRT_32=5.656854249492380195206754896838;
// immutable double M_LOG10_2=0.301029995663981195213738894724;
// immutable double M_SQRT_PI=1.772453850905516027298167483341;
// immutable double M_1_SQRT_2PI=0.398942280401432677939946059934;
// immutable double M_SQRT_2dPI=0.797884560802865355879892119869;
// immutable double M_LN_SQRT_PI=0.572364942924700087071713675677;
// immutable double M_LN_SQRT_2PI=0.918938533204672741780329736406;
// immutable double M_LN_SQRT_PId2=0.225791352644727432363097614947;

extern (C) {
//   double * REAL(Robj x);
//   int * INTEGER(Robj x);
//   const(char) * R_CHAR(Robj x);
//   int * LOGICAL(Robj x);
//   Robj STRING_ELT(Robj x, int i);
//   Robj VECTOR_ELT(Robj x, int i);
//   Robj SET_VECTOR_ELT(Robj x, int i, Robj v);
//   void SET_STRING_ELT(Robj x, int i, Robj v);
//   int Rf_length(Robj x);
//   int Rf_ncols(Robj x);
//   int Rf_nrows(Robj x);
//   extern __gshared Robj R_NilValue;
//   alias RNil = R_NilValue;
//   extern __gshared Robj R_DimSymbol;
//   extern __gshared Robj	R_GlobalEnv;
  
//   void Rf_PrintValue(Robj x);
//   int Rf_isArray(Robj x);
//   int Rf_isInteger(Robj x);
//   int Rf_isList(Robj x);
//   int Rf_isLogical(Robj x);
//   int Rf_isMatrix(Robj x);
//   int Rf_isNull(Robj x);
//   int Rf_isNumber(Robj x);
//   int Rf_isNumeric(Robj x);
//   int Rf_isReal(Robj x);
//   int Rf_isVector(Robj x);
//   int Rf_isVectorList(Robj x);
//   Robj Rf_protect(Robj x);
//   Robj Rf_unprotect(int n);
  // Robj Rf_unprotect_ptr(Robj x);
//   Robj Rf_listAppend(Robj x, Robj y);
//   Robj Rf_duplicate(Robj x);
//   double Rf_asReal(Robj x);
//   int Rf_asInteger(Robj x);
//   Robj Rf_ScalarReal(double x);
//   Robj Rf_ScalarInteger(int x);
//   Robj Rf_ScalarLogical(int x);
//   Robj Rf_getAttrib(Robj x, Robj attr);
//   Robj Rf_setAttrib(Robj x, Robj attr, Robj val);
//   Robj Rf_mkChar(const char * str);
//   Robj Rf_mkString(const char * str);
//   void Rf_error(const char * msg);
//   void R_CheckUserInterrupt();
  
//   Robj Rf_VectorToPairList(Robj x);
    
//   // type is 0 for NILSXP, 13 for integer, 14 for real, 19 for VECSXP
//   Robj Rf_allocVector(uint type, int n);
//   Robj Rf_allocMatrix(uint type, int rows, int cols);
//   Robj Rf_allocList(int);
  
//   // Documented in the R Extensions manual
//   // Install creates an Robj for a symbol name
//   Robj Rf_install(const char * sym);
  
//   // envir.c
//   void Rf_defineVar(Robj symbol, Robj value, Robj rho);
  
//   // GET_SLOT
//   // Use RSymbol(string) to get symbol
//   Robj R_do_slot(Robj obj, Robj symbol);
//   // SET_SLOT
//   Robj R_do_slot_assign(Robj obj, Robj symbol, Robj value);
//   // HAS_SLOT
//   Robj R_has_slot(Robj obj, Robj symbol);
  
//   /* Work with closures */
//   Robj FORMALS(Robj fn);
//   Robj BODY(Robj fn);
//   Robj CLOENV(Robj fn);
  
//   int TYPEOF(Robj obj);
//   void SET_TYPEOF(Robj, int);
  
//   Robj TAG(Robj e);
//   Robj CAR(Robj e);
//   Robj CDR(Robj e);
//   Robj SETCAR(Robj x, Robj y);
//   Robj SETCDR(Robj x, Robj y);
//   Robj SETCADR(Robj x, Robj y);
//   Robj SETCADDR(Robj x, Robj y);
//   Robj SETCADDDR(Robj x, Robj y);
//   Robj SETCAD4R(Robj e, Robj y);
//   Robj Rf_cons(Robj x, Robj y);
  
//   /* External pointers */
// 	Robj R_MakeExternalPtrFn(void function() p, Robj tag, Robj prot);
//   Robj R_MakeExternalPtr(void * p, Robj tag, Robj prot);
//   void * R_ExternalPtrAddr(Robj s);
// 	void function() R_ExternalPtrAddrFn(Robj s);
// 	Robj R_ExternalPtrTag(Robj s);
// 	Robj R_ExternalPtrProtected(Robj s);
// 	void R_ClearExternalPtr(Robj s);
// 	void R_SetExternalPtrAddr(Robj s, void *p);
// 	void R_SetExternalPtrTag(Robj s, Robj tag);
// 	void R_SetExternalPtrProtected(Robj s, Robj p);
	
//   // I don't use these, and don't know enough about them to mess with them
//   // They are documented in the R extensions manual.
//   double gammafn(double);
//   double lgammafn(double);
//   double lgammafn_sign(double, int *);
//   double digamma(double);
//   double trigamma(double);
//   double tetragamma(double);
//   double pentagamma(double);
//   double beta(double, double);
//   double lbeta(double, double);
//   double choose(double, double);
//   double lchoose(double, double);
//   double bessel_i(double, double, double);
//   double bessel_j(double, double);
//   double bessel_k(double, double, double);
//   double bessel_y(double, double);
//   double bessel_i_ex(double, double, double, double *);
//   double bessel_j_ex(double, double, double *);
//   double bessel_k_ex(double, double, double, double *);
//   double bessel_y_ex(double, double, double *);
        
        
//   /** Calculate exp(x)-1 for small x */
//   double expm1(double);
        
//   /** Calculate log(1+x) for small x */
//   double log1p(double);
        
//   /** Returns 1 for positive, 0 for zero, -1 for negative */
//   double sign(double x);
        
//   /** |x|*sign(y)
//    *  Gives x the same sign as y
//    */   
//   double fsign(double x, double y);
        
//   /** R's signif() function */
//   double fprec(double x, double digits);
        
//   /** R's round() function */
//   double fround(double x, double digits);
        
//   /** Truncate towards zero */
//   double ftrunc(double x);
        
//   /** Same arguments as the R functions */ 
//   double dnorm4(double x, double mu, double sigma, int give_log);
//   double pnorm(double x, double mu, double sigma, int lower_tail, int log_p);
//   double qnorm(double p, double mu, double sigma, int lower_tail, int log_p);
//   void pnorm_both(double x, double * cum, double * ccum, int i_tail, int log_p); /* both tails */
//   /* i_tail in {0,1,2} means: "lower", "upper", or "both" :
//      if(lower) return *cum := P[X <= x]
//      if(upper) return *ccum := P[X > x] = 1 - P[X <= x] */

//   /** Same arguments as the R functions */ 
//   double dunif(double x, double a, double b, int give_log);
//   double punif(double x, double a, double b, int lower_tail, int log_p);
//   double qunif(double p, double a, double b, int lower_tail, int log_p);

//   /** These do not allow for passing argument rate as in R 
//       Confirmed that otherwise you call them the same as in R */
//   double dgamma(double x, double shape, double scale, int give_log);
//   double pgamma(double q, double shape, double scale, int lower_tail, int log_p);
//   double qgamma(double p, double shape, double scale, int lower_tail, int log_p);
        
//   /** Unless otherwise noted from here down, if the argument
//    *  name is the same as it is in R, the argument is the same.
//    *  Some R arguments are not available in C */
//   double dbeta(double x, double shape1, double shape2, int give_log);
//   double pbeta(double q, double shape1, double shape2, int lower_tail, int log_p);
//   double qbeta(double p, double shape1, double shape2, int lower_tail, int log_p);

//   /** Use these if you want to set ncp as in R */
//   double dnbeta(double x, double shape1, double shape2, double ncp, int give_log);
//   double pnbeta(double q, double shape1, double shape2, double ncp, int lower_tail, int log_p);
//   double qnbeta(double p, double shape1, double shape2, double ncp, int lower_tail, int log_p);

//   double dlnorm(double x, double meanlog, double sdlog, int give_log);
//   double plnorm(double q, double meanlog, double sdlog, int lower_tail, int log_p);
//   double qlnorm(double p, double meanlog, double sdlog, int lower_tail, int log_p);

//   double dchisq(double x, double df, int give_log);
//   double pchisq(double q, double df, int lower_tail, int log_p);
//   double qchisq(double p, double df, int lower_tail, int log_p);

//   double dnchisq(double x, double df, double ncp, int give_log);
//   double pnchisq(double q, double df, double ncp, int lower_tail, int log_p);
//   double qnchisq(double p, double df, double ncp, int lower_tail, int log_p);

//   double df(double x, double df1, double df2, int give_log);
//   double pf(double q, double df1, double df2, int lower_tail, int log_p);
//   double qf(double p, double df1, double df2, int lower_tail, int log_p);

//   double dnf(double x, double df1, double df2, double ncp, int give_log);
//   double pnf(double q, double df1, double df2, double ncp, int lower_tail, int log_p);
//   double qnf(double p, double df1, double df2, double ncp, int lower_tail, int log_p);

//   double dt(double x, double df, int give_log);
//   double pt(double q, double df, int lower_tail, int log_p);
//   double qt(double p, double df, int lower_tail, int log_p);

//   double dnt(double x, double df, double ncp, int give_log);
//   double pnt(double q, double df, double ncp, int lower_tail, int log_p);
//   double qnt(double p, double df, double ncp, int lower_tail, int log_p);

//   double dbinom(double x, double size, double prob, int give_log);
//   double pbinom(double q, double size, double prob, int lower_tail, int log_p);
//   double qbinom(double p, double size, double prob, int lower_tail, int log_p);

//   double dcauchy(double x, double location, double scale, int give_log);
//   double pcauchy(double q, double location, double scale, int lower_tail, int log_p);
//   double qcauchy(double p, double location, double scale, int lower_tail, int log_p);
        
//   /** scale = 1/rate */
//   double dexp(double x, double scale, int give_log);
//   double pexp(double q, double scale, int lower_tail, int log_p);
//   double qexp(double p, double scale, int lower_tail, int log_p);

//   double dgeom(double x, double prob, int give_log);
//   double pgeom(double q, double prob, int lower_tail, int log_p);
//   double qgeom(double p, double prob, int lower_tail, int log_p);

//   double dhyper(double x, double m, double n, double k, int give_log);
//   double phyper(double q, double m, double n, double k, int lower_tail, int log_p);
//   double qhyper(double p, double m, double n, double k, int lower_tail, int log_p);

//   double dnbinom(double x, double size, double prob, int give_log);
//   double pnbinom(double q, double size, double prob, int lower_tail, int log_p);
//   double qnbinom(double p, double size, double prob, int lower_tail, int log_p);

//   double dnbinom_mu(double x, double size, double mu, int give_log);
//   double pnbinom_mu(double q, double size, double mu, int lower_tail, int log_p);

//   double dpois(double x, double lambda, int give_log);
//   double ppois(double x, double lambda, int lower_tail, int log_p);
//   double qpois(double p, double lambda, int lower_tail, int log_p);

//   double dweibull(double x, double shape, double scale, int give_log);
//   double pweibull(double q, double shape, double scale, int lower_tail, int log_p);
//   double qweibull(double p, double shape, double scale, int lower_tail, int log_p);

//   double dlogis(double x, double location, double scale, int give_log);
//   double plogis(double q, double location, double scale, int lower_tail, int log_p);
//   double qlogis(double p, double location, double scale, int lower_tail, int log_p);

//   double ptukey(double q, double nranges, double nmeans, double df, int lower_tail, int log_p);
//   double qtukey(double p, double nranges, double nmeans, double df, int lower_tail, int log_p);
}

// void addToR(string varname, Robj value) {
// 	Rf_defineVar(Rf_install(toUTFz!(char*)(varname)), value, R_GlobalEnv);
// }
	
	
	
	
	
	
	
	
	
	
