/* Using a custom allocator
 * Work with D allocated data from inside R
 * I don't understand this whole thing well enough to trust it
 * But it works */
import betterr.r;
import std.conv, std.stdio;

extern(C) {
  __gshared void * dataptr;

  struct R_allocator_t {
    void * function(R_allocator_t, int) mem_alloc; /* malloc equivalent */
    void function(R_allocator_t, int) mem_free;  /* free equivalent */
    void * res;                /* reserved (maybe for copy) - must be NULL */
    void * data;               /* custom data for the allocator implementation */
  }

  extern(C) void * thealloc(R_allocator_t al, int n) {
    // What's with the addition?
    // Look at this https://github.com/wch/r-source/blob/5434c055666c6c233c949a56f9dc38e4bc31f265/src/main/memory.c#L2824
    // R_size_t hdrsize = sizeof(SEXPREC_ALIGN);
    // custom_node_alloc(allocator, hdrsize + size * sizeof(VECREC))
    // Assumes you'll allocate an extra hdrsize before the data
    // That appears to be equal to 10*double.sizeof = 80
    // It feels like something's going to blow up, but nothing has
    return dataptr - 80;
  }

  extern(C) void thefreer(R_allocator_t al, int n) {}

  Robj Rf_allocVector3(int type, int length, R_allocator_t * allocator);
}

void main() {
  startR();
  R_allocator_t allocTest;
  allocTest.mem_alloc = &thealloc;
  allocTest.mem_free = &thefreer;
  double[] x;
  foreach(ii; 0..100) {
    x ~= double(ii)*0.2;
  }
  dataptr = x.ptr;
  writeln(x);
  
  /* Now see if we can pass the data to R */
  Robj rx = Rf_allocVector3(14, 100, &allocTest);
  toR(rx, "rx");
  evalRQ("print(rx)");
  evalRQ("print(mean(rx))");
  closeR();
}
