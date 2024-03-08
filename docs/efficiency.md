# Isn't this slow?

No. It isn't slow.

One of the frustrating things about this project is the misconception that calling R means you end up with D code that runs at the same speed as equivalent R code. Using betterr is NOT the equivalent of running R code and capturing the output. Here's a brief explanation that will hopefully clarify why that's wrong and why betterr code can be very performant.

## What this library isn't doing

Let's start by destroying the notion that betterr runs R code and captures the output.

Suppose you're doing a simple simulation. For each replication, you generate random data for variables y and x, do a linear regression of y on x, and save the slope coefficient. In the alternate universe where betterr runs R code and captures the output, you'd write a loop that looks something like this:

```
double[] result;
foreach(ii; 0..1000) {
  R("y <- rnorm(200)");
  R("x <- rnorm(200)");
  R("fit <- lm(y ~ x)");
  result ~= R("coefficients(fit)[2]").output.to!double;
}
```

In the above code, `R` is a function that executes a string of code inside the R interpreter. `.output` converts the R output into a string. The output string is converted to a numeric type and saved. If this example were realistic, it would be pointless to call into R. You'd be writing strings that are snippets of R code. It would be ugly. It would not take advantage of any of D's nice features. It would be even slower than R (which TBH is pretty fast for many things these days).

## This is the code you'd actually write

In reality, your program would look like this:

```
double[] result;
foreach(ii; 0..1000) {
  auto y = rnorm(200);
  auto x = rnorm(200);
  auto fit = lm(y, x);
  result ~= fit.beta[1];
}
```

The first thing to note is that it looks like regular D code. There's not a single snippet of R code.

Further, note that this is pretty fast out of the box. The `rnorm` function calls into R, but R calls into C to generate the random numbers. Similarly, the `lm` function calls into R, but the R code calls a Fortran function, with the convenience of doing various checks to make sure things won't blow up on you. 

`fit.beta` allows you to access the C array holding the estimated coefficients. Since you're using a pointer to access the array, `fit.beta[1]` uses pointer arithmetic to return one element of the array, and it does it at the same speed it'd be done in C.

## If you want to optimize it further

The code above will be plenty fast for most purposes. If you want to optimize it, though, you do have that option. Start by creating Vectors like this:

```
auto y = Vector(200);
auto x = Vector(200);
```

Now fill in the elements using the GNU Scientific Library:

```
import gslheaders;

gsl_rng * r = gsl_rng_alloc(gsl_rng_mt19937);
scope(exit) { gsl_rng_free(r); }

foreach(ii; 0..200) {
  y[ii] = gsl_ran_ugaussian(r);
  x[ii] = gsl_ran_ugaussian(r);
  auto fit = dqrls(y, x);
  result ~= fit.beta[1];
}
```

These optimizations will be incorporated into betterr over time, so that they'll be used by default if the external libraries are available. The point of this example is to demonstrate that you can do arbitrary optimizations because you have a D program and D has full control of the data. The reason R is involved is because it adds new functionality that otherwise is not available in a D program.

## You're actually calling a C library

The functionality of the R interpreter is provided by a shared library written in C, called libR.so on Linux. When R does something, it calls functions in that library. Consider the creation of a new vector holding 10 elements:

```
double(10)
```

The vector is created by calling this C function in libR.so:

```
Rf_allocVector(14, 10);
```

You don't need R to allocate a vector. You can call the function inside libR.so directly:

```
Robj x = Rf_allocVector(14, 10);
```

That's as fast as C, and for good reason - it's 100% C. R is not involved in any way. You can get a pointer to the underlying data array by doing this:

```
double * ptr = REAL(x);
```

and you can get or set elements like this:

```
double z = ptr[2];
ptr[6] = 4.7;
```

Again, R is not involved. Everything here is D calling C functions and working with pointers. What is really nice is that you may want to call a specialized function that is only written in R. You can pass a pointer to the data to the R interpreter like this:

```
toR(x, "xx");
```

and then you can call an R function using that data:

```
evalR("plot(xx)");
```

To sum this up,

- betterr provides a wrapper over a large number of C functions.
- You can efficiently pass data between D and R, because you're only passing a single pointer when you do.
- R provides a very, very large amount of functionality that will never be written in D.

To be useful, you need a complete set of data structures that are common to both languages. Those data structures need to be convenient to use (you should not have to deal with pointers, for instance) and it should feel like you're writing an idiomtic D program. This library provides that type of access to R matrix, vector, list, array, data frame, and ts objects. Since R is garbage collected, you need to protect data from the garbage collector. The data structures of this library handle all of that for you automatically.

## But your program probably won't be slow

Here are some reasons your program is unlikely to be slow - and will almost always be faster than R.

*Access to the underlying data arrays.* As demonstrated above, you have direct access to the data. You get and set elements the same as you would the elements of a `double[]`. Accessing elements through R brings with it considerable overhead. It checks the type of the arguments, checks for valid data, and so on. You bypass that completely. This alone will lead to big speedups over R.

*A lot of R is just a thin wrapper over C functions.* If you want to multiply two large matrices, there'll be a little overhead relative to calling into BLAS directly, but the overhead is minimal relative to the work that needs to be done for the multiplication. Numerous operations like reading in data, estimating regression coefficients, and sorting data are similar, with almost all of the work being done in C.

*You can optimize as much as you want.* For linear algebra, you can use the Matrix package (installed by default) to allocate and perform operations on matrices and vectors. There are examples provided to show how to manage the memory using unique pointers, but you can use any strategy you want, including manually protecting and unprotecting. If you prefer to use OpenBLAS, which doesn't do any checks but gives the best possible performance, you can do that. testblas.d includes some examples. I've never had a need to call BLAS directly, but it's easy to do if you want.

*The R API provides functionality written in C.* You can use [the R API](https://rstudio.github.io/r-manuals/r-exts/The-R-API.html) for random distributions, a variety of functions, and numerical optimization routines.

*Internally, some of the functions use compiled functions rather than calling R.* The Mir libraries are a set of high-quality libraries providing fundamental operations such as calculating the mean, median, and quantiles. Mir is optimized D code.

*Scalar RNG is efficient.* Generators are used to make scalar RNG move quickly. If you need every last cycle, you can call the GNU Scientific Library random number generation functions. ImportC does all the work needed for interoperability. There's a GSL-compatible D implementation of one of L'Ecuyer's parallel RNG functions.

In short, there are not a lot of things you'd want to do where you're stuck with what we think of as "R performance". Additional optimizations beyond those provided by betterr are probably (though not necessarily) a waste of your time.

## There are reasons other than performance to write your program in D

One of the things that annoys me is the claim that performance is the only benefit of using D. That's so far from correct that it's not even wrong. Maybe you want to add a little statistical/numerical functionality to a D program that is otherwise non-numerical. All you need in that case is the functionality, not 100% optimized performance. D is a good general purpose programming language. You can have all of R's functionality plus the beauty of static typing and everything else that comes with writing a program in D.

Along these lines, you can read my thoughts on [scaling to zero here](scaling.md). 
