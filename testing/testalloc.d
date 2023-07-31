/* Using a custom allocator
 * Work with D allocated data from inside R
 * Note: You have to allocate 80 bits to hold the sxpinfo header
 * plus three pointers and the node data.
 * That means you need to allocate 80 bits at the start of the array
 * (10 double values) that you don't use. You could end up with a segfault
 * if you don't allocate that extra memory at the beginning. 
 * 
 * https://cran.r-project.org/doc/manuals/r-release/R-ints.html#The-_0027data_0027 */
import betterr.r;
import std.conv, std.stdio;
import core.stdc.string;

extern(C) {
  __gshared void * dataptr;

  struct R_allocator_t {
    void * function(R_allocator_t, int) mem_alloc; /* malloc equivalent */
    void function(R_allocator_t, int) mem_free;  /* free equivalent */
    void * res;                /* reserved (maybe for copy) - must be NULL */
    void * data;               /* custom data for the allocator implementation */
  }

  extern(C) void * thealloc(R_allocator_t al, int n) {
    return dataptr;
  }

  extern(C) void thefreer(R_allocator_t al, int n) {}

  Robj Rf_allocVector3(int type, int length, R_allocator_t * allocator);
}

Robj asRVector(double[] v) {
  dataptr = v.ptr;
  R_allocator_t allocTest;
  allocTest.mem_alloc = &thealloc;
  allocTest.mem_free = &thefreer;
  return Rf_allocVector3(14, to!int(v.length-10), &allocTest);
}

void main() {
  startR();
  writeln(int.sizeof);
  writeln(uint.sizeof);
  Robj ry = Rf_allocVector(14, 100);
  writeln(ry);
  uint j;
  memcpy(&j, ry+5, 1);
  writeln("j: ", j);
  writeln(ry+10);
  writeln(REAL(ry));
  
  int l, tl;
  memcpy(&l, REAL(ry)-16, 4);
  writeln("l: ", l);
  memcpy(&tl, REAL(ry)-8, 4);
  writeln("tl: ", tl);
  
  R_allocator_t allocTest;
  allocTest.mem_alloc = &thealloc;
  allocTest.mem_free = &thefreer;
  double[] x;
  foreach(ii; 0..10) {
    x ~= 0.0;
  }
  foreach(ii; 0..100) {
    x ~= double(ii)*0.2;
  }
  dataptr = x.ptr;
  writeln(x);
  
  /* Now see if we can pass the data to R */
  Robj rx = Rf_allocVector3(14, 100, &allocTest);
  
  writeln("length: ", rx.length);
  
  toR(rx, "rx");
  evalRQ("print(rx)");
  evalRQ("print(mean(rx))");
  evalRQ("print(2*rx)");
  
  // Put the extra 80 bits for R at the beginning
  double[] xvec;
  foreach(ii; 0..10) {
    xvec ~= 0.0;
  }
  xvec ~= [4.7, 3.9, 3.5, 2.6];
  Robj rxvec = xvec.asRVector;
  double[] yvec;
  foreach(ii; 0..10) {
    yvec ~= 0.0;
  }
  yvec ~= [1.1, 2.2, 3.3, 4.4];
  Robj ryvec = yvec.asRVector;
  
  ryvec.toR("yvec");
  rxvec.toR("xvec");
  evalRQ("print(lm(yvec ~ xvec))");
  
  closeR();
}
