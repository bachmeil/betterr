module betterr.list;

import betterr.r;
import betterr.rdata, betterr.vector;
import std.conv, std.exception, std.stdio, std.string;

struct List {
  RData data;
  alias data this;
  
  this(string code) {
    data = RData(code);
  }
  
  this(List x) {
    this(x.name ~ "[]");
  }
  
  this(long n) {
    this("vector('list', " ~ n.to!string ~ ")");
  }
  
  this(Robj obj) {
    enforce(Rf_isVectorList(obj), "Cannot construct List with anything other than an R list");
    import std.datetime;
    auto ct = Clock.currTime().toISOString;
    string name = "x" ~ ct[0..15] ~ ct[16..$];
    toR(Rf_duplicate(obj), name);
    this(name);
  }
  
  this(RData rd) {
    this(rd.x);
  }
  
  /* This has to return RData because there is no way to know the type
   * You'll have to do a conversion yourself */
  RData opIndex(string n) {
    return RData(data.name ~ "[['" ~ n ~ "']]");
  }
  
  RData opIndex(long ind) {
    return RData(data.name ~ "[[" ~ (ind+1).to!string ~ "]]");
  }
  
  List opIndex(string[] names) {
    auto result = List(names.length);
    foreach(name; names) {
      result[name] = this[name];
    }
    return result;
  }
  
  List opIndex(long[] indexes) {
    auto result = List(indexes.length);
    foreach(counter, index; indexes) {
      result[counter] = this[index];
    }
    return result;
  }    
  
  /* Note that the underlying Robj pointer can (and probably does)
   * change when you update the list. */
  void opIndexAssign(RData d, string n) {
    evalRQ(data.name ~ `[['` ~ n ~ `']] <- ` ~ d.name);
    data.update();
  }
  
  void opIndexAssign(T)(T x, string n) {
    static if (__traits(isArithmetic, T)) {
      evalRQ(data.name ~ `[['` ~ n ~ `']] <- ` ~ x.to!string);
      data.update();
    } else {
      opIndexAssign(x.data, n);
    }
  }
  
  void opIndexAssign(bool v, string n) {
    if (v) {
      evalRQ(data.name ~ `[['` ~ n ~ `']] <- TRUE`);
    } else {
      evalRQ(data.name ~ `[['` ~ n ~ `']] <- FALSE`);
    }
    data.update();
  }
  
  void opIndexAssign(string s, string n) {
    evalRQ(data.name ~ `[['` ~ n ~ `']] <- '` ~ s.replace(`'`, `\'`) ~ `'`);
    data.update();
  }
  
  void opIndexAssign(RData d, long ind) {
    evalRQ(data.name ~ `[[` ~ to!string(ind+1) ~ `]] <- ` ~ d.name);
    data.update();
  }
  
  void opIndexAssign(T)(T a, long ind)
    if (is(T == int) || is(T == double) || is(T == long)) {
      evalRQ(data.name ~ `[[` ~ to!string(ind+1) ~ `]] <- ` ~ a.to!string);
      data.update();
    }
  
  void opIndexAssign(T)(T x, long ind) {
    static if (__traits(isArithmetic, T)) {
      evalRQ(data.name ~ `[[` ~ to!string(ind+1) ~ `]] <- ` ~ x.to!string);
      data.update();
    } else {
      opIndexAssign(x.data, ind);
    }
  }

  void opIndexAssign(bool v, long ind) {
    if (v) {
      evalRQ(data.name ~ `[[` ~ to!string(ind+1) ~ `]] <- TRUE`);
    } else {
      evalRQ(data.name ~ `[[` ~ to!string(ind+1) ~ `]] <- FALSE`);
    }
    data.update();
  }
  
  void opIndexAssign(string s, long ind) {
    evalRQ(data.name ~ `[[` ~ to!string(ind+1) ~ `]] <- '` ~ s.replace(`'`, `\'`) ~ `'`);
    data.update();
  }
  
  List opSlice(long i, long j) {
    return List(data.name ~ "[" ~ to!string(i+1) ~ ":" ~ to!string(j) ~ "]");
  }
  
  string[] names() {
    return stringArray("names(" ~ data.name ~ ")");
  }
  
  int opDollar() {
    return data.x.length;
  }
  
  List dup() {
    return List(this.name ~ "[]");
  }
  
  /* Matrix, Vector, List, DataFrame, etc. */
  void opIndexAssign(T)(T x, long ind) {
    opIndexAssign(x.data, ind);
  }

  void print(string msg="") {
    if (msg != "") { writeln(msg, ":"); }
    printR(data.x);
  }
}
