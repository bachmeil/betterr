# Isn't this slow?

The first objection people have when they see this is that calling into R means you end up with D code that runs at the same speed as equivalent R code. In that case, you should just write your program in R. This view is wrong for two reasons:

- It's not the equivalent of running R code and capturing the output. Most of the time you should end up with good performance, and some of the time you should get performance equivalent to C.
- There are good reasons to write your program in D even if it doesn't run faster than an equivalent R program.

Allow me to elaborate on each point.

## What this library isn't doing

Let's start by clarifying what is not happening by using an example of a simple simulation. For each replication, you generate random data for variables y and x, regress y on x, and save the slope coefficient. You'd set up a loop in D that looks like this:

```
double[] result;
foreach(ii; 0..1000) {
  R("y <- rnorm(200)");
  R("x <- rnorm(200)");
  R("fit <- lm(y ~ x)");
  result ~= R("coefficients(fit)[2]").output.to!double;
}
```

In the above code, `R` is a function that executes a string of code inside the R interpreter. `.output` captures the output from R into a string. The output string is converted to a numeric type and saved. If this example were realistic, it would be pointless to call into R. Your program would be slow. It would be ugly. It would not take advantage of any of D's nice features.

## The code you'd actually write

This is how your program would look in reality:

```
double[] result;
foreach(ii; 0..1000) {
  auto y = rnorm(200);
  auto x = rnorm(200);
  auto fit = lm(y, x);
  result ~= fit.beta[1];
}
```

The first thing to note is that it looks like regular D code. You're not passing around snippets of R code. The second thing to note is that it's pretty fast. The `rnorm` function calls into R, but R calls into C to generate the random numbers. Similarly, the `lm` function calls into R, but the R code calls a C function, with the convenience of doing various checks to make sure things won't blow up. `fit.beta` allows you to access the C array holding the estimated coefficients. Since you're using a pointer to access the array, `fit.beta[1]` uses pointer arithmetic to return one element of the array.

## If you want to optimize it

This code will for most purposes deliver adequate performance. If you want to optimize it, though, you have that option. You could start by using a Generator to generate `y` and `x`. You'd create Vectors like this:

```
auto y = Vector(200);
auto x = Vector(200);
```

You'd create the Generator like this:

```
Generator!"norm" normGen;
```

and you'd fill the elements like this:

```
foreach(ii; 0..200) {
  y[ii] = norm.draw();
  x[ii] = norm.draw();
}
```

This would improve performance by reducing the number of allocations - every call to `rnorm` allocates a new vector to hold the result. The `Generator` generates many draws at once, so it seldom does allocations. If even that isn't good enough, you can call a function written in D or C to take the draws. Note that since `y` and `x` were allocated inside R, the call to do the regression isn't affected:

```
auto fit = lm(y, x);
```

You can optimize the parts you want, but you don't lose access to any functionality when you do.

## What's going on

The functionality of the R interpreter is provided by a shared library written in C, called libR.so. When R does something, it calls functions in that library. Consider the creation of a new vector holding 10 elements:

```
double(10)
```

The vector is created by calling a C function in libR.so:

```
Rf_allocVector(14, 10);
```

You don't need R to allocate a vector. You can call the function inside libR.so directly:

```
Robj x = Rf_allocVector(14, 10);
```

That's as fast as C, because it's 100% C - R is not in any way involved. You can get a pointer to the underlying data array by doing this:

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

- This library provides a wrapper over a large number of C functions.
- You can efficiently pass data between D and R, because you're only passing a single pointer when you do. 

To be useful, you need a complete set of data structures that are common to both languages. Those data structures need to be convenient to use (you should not have to deal with pointers, for instance) and it should feel like you're writing an idiomtic D program. This library provides that type of access to R matrix, vector, list, array, data frame, and ts objects. Since R is garbage collected, you need to protect data from the garbage collector. The data structures of this library handle all of that for you automatically.

## Your program probably won't be slow

Here are some reasons your program is unlikely to be slow - and will almost always be faster than R.

*Access to the underlying data arrays.* As demonstrated above, you have direct access to the underlying data array. You get and set them the same as you would the elements of a `double[]`. Accessing elements in R brings considerable overhead. It checks the type of the arguments, checks for valid data, and so on. You can bypass that completely.

*A lot of R is just a thin wrapper over C functions.* If you want to multiply two large matrices, there'll be a little overhead relative to calling into BLAS directly, but the overhead will be minimal. Numerous operations like reading in data, estimating regression coefficients, and sorting data are similar, with most of the work done by C.

*You can optimize as much as you want.* For linear algebra, you can use the Matrix package (installed by default) to allocate and perform operations on matrices and vectors. There are examples provided to show how to manage the memory using unique pointers, but you can use any strategy you want, including manually protecting and unprotecting. If you prefer to use OpenBLAS, which doesn't do any checks but gives the best possible performance, you can do that. testblas.d includes some examples. I've never had a need to call BLAS directly, but it's easy to do if you want.

*The R API provides functionality written in C.* You can use [the R API](https://rstudio.github.io/r-manuals/r-exts/The-R-API.html) for random distributions, a variety of functions, and numerical optimization routines.

*Internally, some of the functions use Mir rather than calling R.* The Mir libraries are a set of high-quality libraries providing fundamental operations such as calculating the mean, median, and quantiles. Mir is optimized D code.

*Generators are provided to reduce the overhead of scalar RNG.* It's hard to completely remove R from the equation for some tasks. If you want to simulate draws from particular distributions, it's easiest to do that by calling R, especially when you're doing parallel random number generation. The Generator will remove most of the overhead associated with R's RNG functions.

In short, there are not a lot of things you'd want to do where you're stuck with what we think of as "R performance". Additional optimizations beyond those provided in this library are probably not the best use of your time.

## There are reasons other than performance to write your program in D

One of the things that annoys me is the claim that performance is the only benefit of using D. That is, in my opinion, so far from correct that it's not even wrong. Maybe you want to add a little statistical/numerical functionality to a D program that is otherwise non-numerical. You need the functionality in that case, not 100% optimized performance. But beyond that, D is a good general purpose programming language. You can have all of R's functionality plus the beauty of static typing and everything else that comes with writing a program in D.

Along these lines, you can read my thoughts on [scaling to zero here](scaling.md). 
