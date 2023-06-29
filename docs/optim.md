# Optimization

The optimization routines used by [optim](https://web.mit.edu/~r/current/lib/R/library/stats/html/optim.html) are [exposed as part of the R API](https://rstudio.github.io/r-manuals/r-exts/The-R-API.html#optimization). This module provides a convenient interface for D programs needing to solve numerical optimization problems.

## Examples

It's probably easiest in most cases to get started with this module by [studying the examples](https://github.com/bachmeil/betterr/blob/main/testing/testoptim.d).

## Objective Function and Gradient

Regardless of the algorithm, you always need to define the objective function. These routines do minimization, so if you are maximizing a function, with the most common case being maximum likelihood estimation, you need the negative of the objective function. You may also need the gradient function depending on the algorithm.

Both functions need to be marked `extern(C)` because the underlying libraries are written in C. We will minimize the function `x^2 + y^2`. The solution is easily seen to be at `x=0, y=0`. Aliases for them are defined as

```
alias optimfn = double function(int, double*, void*);
alias optimgr = void function(int, double*, double*, void*);
```

For the objective function in this example, we have

```
extern(C) {
  double f(int n, double * par, void * ex) {
    return par[0]*par[0] + par[1]*par[1];
  }

  void g(int n, double * par, double * gr, void * ex) {
    gr[0] = 2*par[0];
    gr[1] = 2*par[1];
  }
}
```

## Result

All optimization algorithms return a struct `OptimSolution` with the following data:

- double[] sol: The solution vector
- double[] init: The starting values
- bool fail: Takes the value true if the optimization failed
- int fncount: Number of function evaluations
- double objFunction: The value of the objective function at sol
- string algorithm: The algorithm used

## Nelder-Mead

This is the default algorithm used by optim. It doesn't require derivatives, and it's robust, with the main downside being that it's slow relative to the other algorithms.

You call the constructor with the objective function:

```
auto nm = NelderMead(&f);
```

To do the optimization, you call the `solve` method with a vector of starting values:

```
OptimSolution sol = nm.solve([3.5, -5.5]);
```

Alternatively, you can use the lower-level interface that is available for the rare case where you need more flexibility. The first argument is a `double *` to an array holding the solution, the second is a `double *` to an array holding the starting values, the third is an `int` holding the number of parameters, and the fourth (optional) is a `void *` that points to data used to evaluate the objective function and/or the gradient function.

```
double[] starting = [3.5, -5.5];
double[] solution = [0.0, 0.0];
sol = nm.solve(solution.ptr, starting.ptr, 2);
```

Finally, you could call the function `nmmin` directly, but since that's less convenient, more error-prone, and provides no additional functionality beyond the low-level `solve` function, I will not write about it further. You can read the betterr.optim source file and the R API documentation if you really want those details.

The NelderMead struct contains numerous parameters you can adjust, similar to the options you set when calling `optim` in R.

- double abstol = -double.infinity;
- double intol = 0.00000001;
- double alpha = 1.0;
- double beta = 0.5;
- double gamma = 2.0;
- bool trace = false;
- int maxit = 500;

See [the help file for optim](https://web.mit.edu/~r/current/lib/R/library/stats/html/optim.html) if you want more information about those parameters. `intol` is used to calculate the convergence tolerance. In [the source code](https://github.com/wch/r-source/blob/8d985707638d3e1b20df24fe48c7e47347656f8f/src/appl/optim.c#L302), the convergence tolerance is defined as `convtol = intol * (fabs(f) + intol)`, where `f` is the value of the objective function evaluated at the current parameter vector. I believe that what is called `intol` in the C source is the same as `reltol` in the control argument to `optim`, because the help says "Relative convergence tolerance. The algorithm stops if it is unable to reduce the value by a factor of reltol * (abs(val) + reltol) at a step. Defaults to sqrt(.Machine$double.eps), typically about 1e-8." That is why I have set the default value for `intol` to 0.00000001.

## BFGS

For BFGS, you use a `BFGS` struct. The options you can set are

```
optimfn fn;
optimgr gr;
double abstol = -double.infinity;
double reltol = 0.00000001;
bool trace = false;
int maxit = 100;
int report = 10;
```

`fn` is the function evaluating the objective function. `gr` is the function evaluating the gradient. The other options are similar to [those of optim](https://web.mit.edu/~r/current/lib/R/library/stats/html/optim.html). Note that the gradient function is required.

## Conjugate Gradient

For conjugate gradient, you use a `ConjugateGradient` struct. The available options are

```
optimfn fn;
optimgr gr;
double abstol = -double.infinity;
double reltol = 0.00000001;
bool trace = false;
int maxit = 100;
int type = 1; // 1 is Fletcher-Reeves, 2 is Polak-Ribiere, 3 is Beale-Sorenson
```

The gradient function is required.

## L-BFGS (Bounds constraints)

This algorithm allows you to impose bounds constraints, such that there are upper and/or lower limits on individual parameters. For this algorithm, you use the `Bounded` struct. The available options are

```
optimfn fn;
optimgr gr;
int lmm = 5; // Maximum number of variable metric corrections
double factr = 0.0000001;
double pgtol = 0.0;
bool trace = false;
int maxit = 100;
int type = 1; // 1 is Fletcher-Reeves, 2 is Polak-Ribiere, 3 is Beale-Sorenson
int report = 10;
```

The gradient function is required.

## Simulated Annealing

This algorithm is not documented in a way that allows you to use it. After a debugging session inside the R source code, I found out what was going on, but the end user of this library does not need to worry about it. You use simulated annealing through the `SA` struct. The available options are

```
optimfn fn;
int trace = 100; // Apparently, this function treats trace like the others treat report
int maxit = 10000;
int tmax = 10;
double temp = 10.0;
```


