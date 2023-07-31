import betterr.r, betterr.matrix, betterr.random, betterr.vector;
import std.conv, std.stdio;

/* This is used for custom allocation so you can work with data 
 * allocated by D within R. Note that you need to allocate an extra 80
 * bits at the beginning so R can work with your data. If you let R
 * allocate the data, you don't need to mess with any of this. */
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

/* Here's the code that does the actual work. */
void main() {
  startR();
  
  evalRQ("library(rstan)");
  "
data {
  int<lower=0> N;
  real y[N];
} 

parameters {
  real mu;
} 

model {
  mu ~ normal(0, 10);
  y ~ normal(mu, 1); 
} 

".robj.toR("stanmodelcode");
  "normal1".robj.toR("model_name");
  double[] y;
  // The 80 bits R needs
  foreach(ii; 0..10) { y ~= 0.0; }
  Vector rand = rnorm(20); // You could send rand directly to R, but then you wouldn't be using a custom allocator
  foreach(ii; 0..20) {
    y ~= rand[ii];
  }
  Robj ry = y.asRVector;
  
  // Send the data from D to R
  ry.toR("y");
  
  // Run the R code verbatim, in practice someone familiar with RStan
  // would want to create a better interface so it feels like D code
  // I don't have the knowledge to do that
  evalRQ("rr <- stan_model(model_code = stanmodelcode, model_name = model_name, verbose = TRUE)");
  evalRQ("dat <- list(N = 20, y = y)");
  evalRQ("f <- sampling(rr, data = dat, init = 0, iter = 2012, sample_file = 'norm1.csv')");
  auto f = Matrix(`summary(f)[["summary"]]`);
  writeln(f);
  
  closeR();
}
