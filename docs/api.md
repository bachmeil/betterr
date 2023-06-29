# R API

There are some useful bits pulled from the R API. They are imported as part of the betterr.r module.

## Constants

Documentation is [here](https://rstudio.github.io/r-manuals/r-exts/The-R-API.html#mathematical-constants).

```
immutable double M_E
immutable double M_LOG2E
immutable double M_LOG10E
immutable double M_LN2
immutable double M_LN10 
immutable double M_PI
immutable double M_2PI
immutable double M_PI_2
immutable double M_PI_4
immutable double M_1_PI
immutable double M_2_PI
immutable double M_2_SQRTPI
immutable double M_SQRT2
immutable double M_SQRT1_2
immutable double M_SQRT_3
immutable double M_SQRT_32
immutable double M_LOG10_2
immutable double M_SQRT_PI
immutable double M_1_SQRT_2PI
immutable double M_SQRT_2dPI
immutable double M_LN_SQRT_PI
immutable double M_LN_SQRT_2PI
immutable double M_LN_SQRT_PId2
```

## Probability distributions

These are C functions, but they work in most cases similar to the R versions. Documentation is [here](https://rstudio.github.io/r-manuals/r-exts/The-R-API.html#distribution-functions).

```
/** Same arguments as the R functions */ 
double dnorm4(double x, double mu, double sigma, int give_log);
double pnorm(double x, double mu, double sigma, int lower_tail, int log_p);
double qnorm(double p, double mu, double sigma, int lower_tail, int log_p);
void pnorm_both(double x, double * cum, double * ccum, int i_tail, int log_p); /* both tails */
/* i_tail in {0,1,2} means: "lower", "upper", or "both" :
   if(lower) return *cum := P[X <= x]
   if(upper) return *ccum := P[X > x] = 1 - P[X <= x] */

/** Same arguments as the R functions */ 
double dunif(double x, double a, double b, int give_log);
double punif(double x, double a, double b, int lower_tail, int log_p);
double qunif(double p, double a, double b, int lower_tail, int log_p);

/** These do not allow for passing argument rate as in R 
    Confirmed that otherwise you call them the same as in R */
double dgamma(double x, double shape, double scale, int give_log);
double pgamma(double q, double shape, double scale, int lower_tail, int log_p);
double qgamma(double p, double shape, double scale, int lower_tail, int log_p);
      
/** Unless otherwise noted from here down, if the argument
 *  name is the same as it is in R, the argument is the same.
 *  Some R arguments are not available in C */
double dbeta(double x, double shape1, double shape2, int give_log);
double pbeta(double q, double shape1, double shape2, int lower_tail, int log_p);
double qbeta(double p, double shape1, double shape2, int lower_tail, int log_p);

/** Use these if you want to set ncp as in R */
double dnbeta(double x, double shape1, double shape2, double ncp, int give_log);
double pnbeta(double q, double shape1, double shape2, double ncp, int lower_tail, int log_p);
double qnbeta(double p, double shape1, double shape2, double ncp, int lower_tail, int log_p);

double dlnorm(double x, double meanlog, double sdlog, int give_log);
double plnorm(double q, double meanlog, double sdlog, int lower_tail, int log_p);
double qlnorm(double p, double meanlog, double sdlog, int lower_tail, int log_p);

double dchisq(double x, double df, int give_log);
double pchisq(double q, double df, int lower_tail, int log_p);
double qchisq(double p, double df, int lower_tail, int log_p);

double dnchisq(double x, double df, double ncp, int give_log);
double pnchisq(double q, double df, double ncp, int lower_tail, int log_p);
double qnchisq(double p, double df, double ncp, int lower_tail, int log_p);

double df(double x, double df1, double df2, int give_log);
double pf(double q, double df1, double df2, int lower_tail, int log_p);
double qf(double p, double df1, double df2, int lower_tail, int log_p);

double dnf(double x, double df1, double df2, double ncp, int give_log);
double pnf(double q, double df1, double df2, double ncp, int lower_tail, int log_p);
double qnf(double p, double df1, double df2, double ncp, int lower_tail, int log_p);

double dt(double x, double df, int give_log);
double pt(double q, double df, int lower_tail, int log_p);
double qt(double p, double df, int lower_tail, int log_p);

double dnt(double x, double df, double ncp, int give_log);
double pnt(double q, double df, double ncp, int lower_tail, int log_p);
double qnt(double p, double df, double ncp, int lower_tail, int log_p);

double dbinom(double x, double size, double prob, int give_log);
double pbinom(double q, double size, double prob, int lower_tail, int log_p);
double qbinom(double p, double size, double prob, int lower_tail, int log_p);

double dcauchy(double x, double location, double scale, int give_log);
double pcauchy(double q, double location, double scale, int lower_tail, int log_p);
double qcauchy(double p, double location, double scale, int lower_tail, int log_p);
      
/** scale = 1/rate */
double dexp(double x, double scale, int give_log);
double pexp(double q, double scale, int lower_tail, int log_p);
double qexp(double p, double scale, int lower_tail, int log_p);

double dgeom(double x, double prob, int give_log);
double pgeom(double q, double prob, int lower_tail, int log_p);
double qgeom(double p, double prob, int lower_tail, int log_p);

double dhyper(double x, double m, double n, double k, int give_log);
double phyper(double q, double m, double n, double k, int lower_tail, int log_p);
double qhyper(double p, double m, double n, double k, int lower_tail, int log_p);

double dnbinom(double x, double size, double prob, int give_log);
double pnbinom(double q, double size, double prob, int lower_tail, int log_p);
double qnbinom(double p, double size, double prob, int lower_tail, int log_p);

double dnbinom_mu(double x, double size, double mu, int give_log);
double pnbinom_mu(double q, double size, double mu, int lower_tail, int log_p);

double dpois(double x, double lambda, int give_log);
double ppois(double x, double lambda, int lower_tail, int log_p);
double qpois(double p, double lambda, int lower_tail, int log_p);

double dweibull(double x, double shape, double scale, int give_log);
double pweibull(double q, double shape, double scale, int lower_tail, int log_p);
double qweibull(double p, double shape, double scale, int lower_tail, int log_p);

double dlogis(double x, double location, double scale, int give_log);
double plogis(double q, double location, double scale, int lower_tail, int log_p);
double qlogis(double p, double location, double scale, int lower_tail, int log_p);

double ptukey(double q, double nranges, double nmeans, double df, int lower_tail, int log_p);
double qtukey(double p, double nranges, double nmeans, double df, int lower_tail, int log_p);
```

## Special functions

The documentation for these functions [can be found here](https://rstudio.github.io/r-manuals/r-exts/The-R-API.html#mathematical-functions).

```
double gammafn(double);
double lgammafn(double);
double lgammafn_sign(double, int *);
double digamma(double);
double trigamma(double);
double tetragamma(double);
double pentagamma(double);
double beta(double, double);
double lbeta(double, double);
double choose(double, double);
double lchoose(double, double);
double bessel_i(double, double, double);
double bessel_j(double, double);
double bessel_k(double, double, double);
double bessel_y(double, double);
double bessel_i_ex(double, double, double, double *);
double bessel_j_ex(double, double, double *);
double bessel_k_ex(double, double, double, double *);
double bessel_y_ex(double, double, double *);
```

## A few others

```
/** Calculate exp(x)-1 for small x */
double expm1(double);
      
/** Calculate log(1+x) for small x */
double log1p(double);
      
/** Returns 1 for positive, 0 for zero, -1 for negative */
double sign(double x);
      
/** |x|*sign(y)
 *  Gives x the same sign as y
 */   
double fsign(double x, double y);
      
/** R's signif() function */
double fprec(double x, double digits);
      
/** R's round() function */
double fround(double x, double digits);
      
/** Truncate towards zero */
double ftrunc(double x);
```

## Working with Robj

The struct that holds all R data is called `Robj`, for "R object". If you
want low-level access for performance reasons, you may want some helper
functions for working with Robj structs.

### Basics

```
void Rf_PrintValue(Robj x);
int Rf_isArray(Robj x);
int Rf_isInteger(Robj x);
int Rf_isList(Robj x);
int Rf_isLogical(Robj x);
int Rf_isMatrix(Robj x);
int Rf_isNull(Robj x);
int Rf_isNumber(Robj x);
int Rf_isNumeric(Robj x);
int Rf_isReal(Robj x);
int Rf_isVector(Robj x);
int Rf_isVectorList(Robj x);
```

### Protecting objects from the R garbage collector

- `Robj Rf_protect(Robj x)`
- `Robj Rf_unprotect(int n)`: Unprotect the n most recent protected objects.
- `Robj Rf_unprotect_ptr(Robj x)`: Unprotect a specific object.

### Copying an object

- `Robj Rf_duplicate(Robj x)`

### Scalar conversions

```
double Rf_asReal(Robj x);
int Rf_asInteger(Robj x);
Robj Rf_ScalarReal(double x);
Robj Rf_ScalarInteger(int x);
Robj Rf_ScalarLogical(int x);
```

There are also these functions:

```
double scalar(Robj rx)
double scalar(T: double)(Robj rx)
int scalar(T: int)(Robj rx)
long scalar(T: long)(Robj rx)
ulong scalar(T: ulong)(Robj rx)
string scalar(T: string)(Robj rx)
double scalar(string name)
double scalar(T: double)(string name)
int scalar(T: int)(string name)
long scalar(T: long)(string name)
ulong scalar(T: ulong)(string name)
string scalar(T: string)(string name)
```

### Getting and setting attributes on R objects

```
Robj Rf_getAttrib(Robj x, Robj attr);
Robj Rf_setAttrib(Robj x, Robj attr, Robj val);
```

### Creating a string inside R

- `Robj Rf_mkString(const char * str)`

### Create an R object that represents a particular symbol

- `Robj Rf_install(const char * sym)`
- `Robj RSymbol(string sym)`: Easier than `Rf_install` because you can pass a string.

### Working with S4 slots

```
Robj R_do_slot(Robj obj, Robj symbol);
Robj R_do_slot_assign(Robj obj, Robj symbol, Robj value);
Robj R_has_slot(Robj obj, Robj symbol);
```

### Special values

- `RNil`: R NULL
- `Robj RTrue()`: Mostly used as a flag passed to C functions.
- `Robj RFalse()`
- `R_DimSymbol`: R object that refers to the "dim" symbol.
- `R_GlobalEnv`: R object that refers to the global environment.

### Pointers to the array holding data inside an Robj

```
double * REAL(Robj x);
int * INTEGER(Robj x);
const(char) * R_CHAR(Robj x);
int * LOGICAL(Robj x);
```

### Get the length of an R object and the number of rows or columns of a matrix

```
int Rf_length(Robj x);
int Rf_ncols(Robj x);
int Rf_nrows(Robj x);
```

### Convert D data structures to Robj and vice versa

```
Robj robj(double x)
Robj robj(double[] v)
Robj robj(int x)
Robj robj(string s)
Robj robj(string[] sv)
string toString(Robj cstr)
string toString(Robj sv, int ii)
string[] stringArray(Robj sv)
string[] stringArray(string name)
```

### Place objects into R

```
void toR(T)(T x, string name)
void toR(Robj x, string name)
void toR(string[] s, string name)
```

### Evaluating R code

- `Robj evalR(string cmd)`: Evaluate `cmd` inside R and return an Robj, a pointer to the data.
- `void evalRQ(string cmd)`: Evaluate `cmd` inside R but return nothing.
- `void evalRQ(string[] cmds)`: Evaluate all of the commands in `cmds` one at a time, returning nothing.

