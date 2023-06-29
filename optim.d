module betterr.optim;

import betterr.r;
import std.conv, std.file, std.stdio;

struct OptimSolution {
  double[] sol;
  double[] init;
  bool fail;
  int fncount;
  double objFunction;
  string algorithm;

  void print(string msg="") {
    if (msg.length > 0) {
      writeln(msg);
    }
    writeln("---");
    writeln("Output of " ~ algorithm ~ " minimization");
    writeln("---");
    writeln("Solution vector: ", sol);
    writeln("Objective function at the solution: ", objFunction);
    writeln("Failure code: ", fail);
    writeln("Number of function evaluations: ", fncount);
    writeln("Initial values: ", init);
  }
}

struct NelderMead {
  optimfn fn;
  double abstol = -double.infinity;
  double intol = 0.00000001;
  double alpha = 1.0;
  double beta = 0.5;
  double gamma = 2.0;
  bool trace = false;
  int maxit = 500;

  this(optimfn f) {
    fn = f;
  }

  OptimSolution solve(double * par, double * init, int npar, void * ex = null) {
    double objFunction;
    int fail;
    int fncount;
    OptimSolution result;
    foreach(ii; 0..npar) {
      result.init ~= init[ii];
    }
    nmmin(npar, init, par, &objFunction, fn, &fail, abstol, intol, ex,
      alpha, beta, gamma, trace.to!int, &fncount, maxit);
    result.sol = par[0..npar];
    result.fail = fail.to!bool;
    result.fncount = fncount;
    result.objFunction = objFunction;
    result.algorithm = "Nelder-Mead";
    return result;
  }

  OptimSolution solve(T)(T starting, void * ex = null) {
    double[] par;
    double[] init;
    foreach(ii; 0..starting.length) {
      /* Both of these are mutated, so we need to create both. */
      par ~= starting[ii];
      init ~= starting[ii];
    }
    return solve(par.ptr, init.ptr, starting.length.to!int, ex);
  }
}

struct BFGS {
  optimfn fn;
  optimgr gr;
  double abstol = -double.infinity;
  double reltol = 0.00000001;
  bool trace = false;
  int maxit = 100;
  int report = 10;

  this(optimfn f, optimgr g) {
    fn = f;
    gr = g;
  }

  /* I don't know what mask is, but the default value works. */
  OptimSolution solve(double[] starting, void * ex = null, int[] mask = []) {
    double objFunction;
    int fail;
    int fncount;
    int grcount;
    if (mask.length == 0) {
      foreach (_; 0..starting.length) {
        mask ~= 1;
      }
    }
    double[] init = starting.dup;
    vmmin(to!int(starting.length), starting.ptr, &objFunction, fn, gr, maxit,
      to!int(trace), mask.ptr, abstol, reltol, report, ex, &fncount, &grcount, &fail);
    return OptimSolution(starting, init, to!bool(fail), fncount, objFunction, "BFGS");
  }

  OptimSolution solve(T)(T starting, void * ex = null, int[] mask = []) {
    double objFunction;
    int fail;
    int fncount;
    int grcount;
    if (mask.length == 0) {
      foreach (_; 0..starting.length) {
        mask ~= 1;
      }
    }
    double[] init;
    double[] _starting;
    foreach(ii; 0..starting.length) {
      init ~= starting[ii];
      _starting ~= starting[ii];
    }
    vmmin(to!int(starting.length), starting.ptr, &objFunction, fn, gr, maxit,
      to!int(trace), mask.ptr, abstol, reltol, report, ex, &fncount, &grcount, &fail);
    return OptimSolution(_starting, init, to!bool(fail), fncount, objFunction, "BFGS");
  }
}

struct ConjugateGradient {
  optimfn fn;
  optimgr gr;
  double abstol = -double.infinity;
  double reltol = 0.00000001;
  bool trace = false;
  int maxit = 100;
	int type = 1; // 1 is Fletcher-Reeves, 2 is Polak-Ribiere, 3 is Beale-Sorenson

  this(optimfn f, optimgr g) {
    fn = f;
    gr = g;
  }

  OptimSolution solve(double[] starting, void * ex = null) {
    double[] result = starting.dup;
    double[] init = starting.dup; // I think this is necessary. nmmin must mutate both init and result.
    double objFunction;
    int fail;
    int fncount;
    int grcount;
    cgmin(to!int(init.length), init.ptr, result.ptr, &objFunction, fn, gr, &fail, 
			abstol, reltol, ex, type, to!int(trace), &fncount, &grcount, maxit);
    return OptimSolution(result, starting, to!bool(fail), fncount, objFunction, "Conjugate Gradient");
  }

  OptimSolution solve(T)(T starting, void * ex = null) {
    double[] result = starting.dup;
    double[] init = starting.dup; // I think this is necessary. nmmin must mutate both init and result.
    double objFunction;
    int fail;
    int fncount;
    int grcount;
    cgmin(to!int(init.length), init.ptr, result.ptr, &objFunction, fn, gr, &fail, 
			abstol, reltol, ex, type, to!int(trace), &fncount, &grcount, maxit);
    return OptimSolution(result, starting, to!bool(fail), fncount, objFunction, "Conjugate Gradient");
  }
}

struct Bounded {
  optimfn fn;
  optimgr gr;
  int lmm = 5; // Maximum number of variable metric corrections
  double factr = 0.0000001;
  double pgtol = 0.0;
  bool trace = false;
  int maxit = 100;
	int type = 1; // 1 is Fletcher-Reeves, 2 is Polak-Ribiere, 3 is Beale-Sorenson
	int report = 10;

  this(optimfn f, optimgr g) {
    fn = f;
    gr = g;
  }

  // bounds: For each parameter, identify bounds (nbd in the R source) as
  //   0 for unbounded
  //   1 for lower bound
  //   2 for lower and upper bound
  //   3 for upper bound
  OptimSolution solve(double[] starting, double[] lower, double[] upper, int[] bounds, void * ex = null) {
    double[] result = starting.dup;
    double objFunction;
    int fail;
    int fncount;
    int grcount;
    char[60] msg; // This is the maximum string length in the optim function in optim.c (currently line 327)
    lbfgsb(to!int(starting.length), lmm, result.ptr, lower.ptr, upper.ptr, bounds.ptr, &objFunction, fn, gr,
 			&fail, ex, factr, pgtol, &fncount, &grcount, maxit, &msg[0], to!int(trace), report);
    return OptimSolution(result, starting, to!bool(fail), fncount, objFunction, "LBFGS-B");
  }
  
  OptimSolution solve(T)(T _starting, double[] lower, double[] upper, int[] bounds, void * ex = null) {
    double[] result;
    double[] starting;
    foreach(ii; 0.._starting.length) {
      result ~= _starting[ii];
      starting ~= _starting[ii];
    }
    double objFunction;
    int fail;
    int fncount;
    int grcount;
    char[60] msg; // This is the maximum string length in the optim function in optim.c (currently line 327)
    lbfgsb(to!int(starting.length), lmm, result.ptr, lower.ptr, upper.ptr, bounds.ptr, &objFunction, fn, gr,
 			&fail, ex, factr, pgtol, &fncount, &grcount, maxit, &msg[0], to!int(trace), report);
    return OptimSolution(result, starting, to!bool(fail), fncount, objFunction, "LBFGS-B");
  }
}

struct SA {
  optimfn fn;
  int trace = 100; // Apparently, this function treats trace like the others treat report
  int maxit = 10000;
	int tmax = 10;
	double temp = 10.0;

  this(optimfn f) {
    fn = f;
  }

  OptimSolution solve(double[] starting) {
    // Called by function genptry in src/appl/optim.c
    // Will segfault if we don't do this - as I learned after a long debugging session
    // This will make your program die: `if (!isNull(OS->R_gcall))` because OS will be null 
    struct Unused {
      Robj R_fcall;
      Robj R_gcall;
    }
    double[] init = starting.dup;
    double objFunction = 0;
    Unused unused;
    unused.R_fcall = RNil;
    unused.R_gcall = RNil;
    void * ex = cast(void*) &unused;
    samin(to!int(starting.length), starting.ptr, &objFunction, fn, maxit, tmax, temp, trace, ex);
    return OptimSolution(starting, init, false, 0, objFunction, "Simulated Annealing");
  }

  OptimSolution solve(T)(T _starting) {
    // Called by function genptry in src/appl/optim.c
    // Will segfault if we don't do this - as I learned after a long debugging session
    // This will make your program die: `if (!isNull(OS->R_gcall))` because OS will be null 
    struct Unused {
      Robj R_fcall;
      Robj R_gcall;
    }
    double[] starting;
    double[] init;
    foreach(ii; 0.._starting.length) {
      starting ~= _starting[ii];
      init ~= _starting[ii];
    }
    double objFunction = 0;
    Unused unused;
    unused.R_fcall = RNil;
    unused.R_gcall = RNil;
    void * ex = cast(void*) &unused;
    samin(to!int(starting.length), starting.ptr, &objFunction, fn, maxit, tmax, temp, trace, ex);
    return OptimSolution(starting, init, false, 0, objFunction, "Simulated Annealing");
  }
}

extern(C) {
  alias optimfn = double function(int, double*, void*);
  alias optimgr = void function(int, double*, double*, void*);

  void nmmin(int n, double *xin, double *x, double *Fmin, optimfn fn,
    int *fail, double abstol, double intol, void *ex,
	  double alpha, double beta, double gamma, int trace,
	  int *fncount, int maxit);

  void vmmin(int n, double *x, double *Fmin,
	   optimfn fn, optimgr gr, int maxit, int trace,
	   int *mask, double abstol, double reltol, int nREPORT,
	   void *ex, int *fncount, int *grcount, int *fail);
	   
	void cgmin(int n, double *xin, double *x, double *Fmin,
	   optimfn fn, optimgr gr, int *fail, double abstol,
	   double intol, void *ex, int type, int trace,
	   int *fncount, int *grcount, int maxit);
	   
	void lbfgsb(int n, int lmm, double *x, double *lower,
	    double *upper, int *nbd, double *Fmin, optimfn fn,
	    optimgr gr, int *fail, void *ex, double factr,
	    double pgtol, int *fncount, int *grcount,
	    int maxit, char *msg, int trace, int nREPORT);
	    
	void samin(int n, double *x, double *Fmin, optimfn fn, int maxit,
	   int tmax, double temp, int trace, void *ex);
}
