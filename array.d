module betterr.array;
import betterr.rdata;
import betterr.r;
import betterr.matrix, betterr.vector;
import std.conv, std.exception, std.range, std.stdio;

struct RArray {
  long dims;
  long[] lengths;
  long length;
  RData data;
  double * ptr;
  alias data this;

  this(string code) {
    data = RData(code);
    ptr = REAL(data.x);
  }

  this(long[] dd ...) {
    dims = dd.length;
    lengths = dd;
    length = 1;
    string tmp;
    foreach(ii, d; dd) {
      length *= d;
      if (ii > 0) {
        tmp ~= ", ";
      }
      tmp ~= d.to!string;
    }
    this("array(0.0, dim=c(" ~ tmp ~ "))");
  }
  
  double opIndex(long[] dd ...) {
    enforce(dd.length == dims, "Wrong number of dimensions on RArray");
    long start = 0;
    long end = length;
    double result;
    foreach(_ii; 0..dims) {
      long ii = (dims - 1) - _ii;
      if (ii == 0) {
        result = ptr[start+dd[ii]];
      } else {
        long sublength = end - start;
        end = start+(dd[ii]+1)*(sublength / lengths[ii]);
        start = start+dd[ii]*(sublength / lengths[ii]);
      }
    }
    return result;
  }
  
  Subarray sub() {
    Subarray result;
    result.rarray = &this;
    result.dims = this.dims;
    result.lengths = this.lengths;
    return result;
  }
  
  void opIndexAssign(double val, long[] dd ...) {
    enforce(dd.length == dims, "Wrong number of dimensions on RArray");
    long start = 0;
    long end = length;
    foreach(_ii; 0..dims) {
      long ii = (dims - 1) - _ii;
      if (ii == 0) {
        ptr[start+dd[ii]] = val;
      } else {
        long sublength = end - start;
        end = start+(dd[ii]+1)*(sublength / lengths[ii]);
        start = start+dd[ii]*(sublength / lengths[ii]);
      }
    }
  }
  
  void opAssign(T)(T v) {
    enforce(this.length == v.length, "Dimensions of RArray do not match what it is filled with");
    foreach(ii; 0..this.length) {
      ptr[ii] = v[ii];
    }
  }
  
  void opAssign(double a) {
    foreach(ii; 0..this.length) {
      ptr[ii] = a;
    }
  }
  
  void opAssign(Matrix m) {
    enforce(this.length == m.rows*m.cols, "Number of elements does not match");
    foreach(ii; 0..this.length) {
      ptr[ii] = m.ptr[ii];
    }
  }
  
  long opDollar(int dim)() {
    return lengths[dim];
  }
}

struct Subarray {
  RArray * rarray;
  long dims;
  long[] lengths;
  
  string opSlice(long dim)(long begin, long end) {
    return to!string(begin+1) ~ ":" ~ to!string(end);
  }

  /* I don't know of any way to do this using opIndex inside RArray without overriding
   * the opIndex above, in which case you lose the ability to work with
   * scalar elements directly - a huge efficiency problem. */
  RArray opIndex(...) {
    string a = _arguments.to!string;
    enforce(_arguments.length == rarray.dims, "Wrong number of index arguments");
    import core.vararg;
    string cmd = rarray.name ~ "[";
    foreach(dim; 0.._arguments.length) {
      if (dim > 0) {
        cmd ~= ", ";
      }
      if (_arguments[dim] == typeid(long)) {
        cmd ~= (va_arg!(long)(_argptr)+1).to!string;
      } 
      if (_arguments[dim] == typeid(int)) {
        cmd ~= (va_arg!(int)(_argptr)+1).to!string;
      } 
      else if (_arguments[dim] == typeid(string)) {
        cmd ~= va_arg!(string)(_argptr);
      }
      else if (_arguments[dim] == typeid(AllElements)) {
        cmd ~= "1:" ~ lengths[dim].to!string;
      }
      else {
        enforce(false, "Index argument needs to be type long, int, string, or AllElements");
      }
    }
    cmd ~= "]";
    writeln(cmd);
    auto result = RArray(cmd);
    auto v = IntVector("dim(" ~ rarray.name ~ ")");
    result.dims = v.length;
    foreach(ii; 0..v.length) {
      result.lengths ~= v.ptr[ii];
    }
    return result;
  }
  
  void opIndexAssign(T)(T obj, ...) {
    string argno = _arguments.length.to!string;
    enforce(_arguments.length == dims, "Wrong number of index arguments " ~ argno);
    import core.vararg;
    string cmd = rarray.name ~ "[";
    foreach(dim; 0.._arguments.length) {
      if (dim > 0) {
        cmd ~= ", ";
      }
      if (_arguments[dim] == typeid(long)) {
        cmd ~= (va_arg!(long)(_argptr)+1).to!string;
      } 
      if (_arguments[dim] == typeid(int)) {
        cmd ~= (va_arg!(int)(_argptr)+1).to!string;
      } 
      else if (_arguments[dim] == typeid(string)) {
        cmd ~= va_arg!(string)(_argptr);
      }
      else if (_arguments[dim] == typeid(AllElements)) {
        cmd ~= "1:" ~ lengths[dim].to!string;
      }
      else {
        enforce(false, "Index argument needs to be type long, int, string, or AllElements");
      }
    }
    cmd ~= "] <- " ~ obj.name;
    writeln(cmd);
    evalRQ(cmd);
    rarray.data.x = evalR(rarray.data.name);
  }
}
