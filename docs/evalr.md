# Evaluating arbitrary R code

One of the nice features of this library is that it provides a set of data structures that are sufficient to work with most of the output from functions inside R. There are two ways to evaluate code. The first is with `evalR`. It runs the code you've specified and returns the output as an Robj struct (which is how ALL data in R is stored). You could do the following:

```
Robj x = evalR("rnorm(15)");
```

`x` would hold a Robj struct with a pointer to the output of that command. There are two problems with using this primitive approach:

- The storage for the output of `rnorm(15)` was allocated by R. Since R has no way to know what you're doing with the data inside your D program, it could be reclaimed by the garbage collector at any time. You could wrap the `evalR` command like this: `Rf_protect(evalR("rnorm(15)"))`. That would prevent R from ever collecting the data, so to prevent a memory leak, you'd need a corresponding call `Rf_unprotect_ptr(x)` or `Rf_unprotect(1)`. See [the R extensions manual](https://rstudio.github.io/r-manuals/r-exts/System-and-foreign-language-interfaces.html#handling-the-effects-of-garbage-collection) for further information. The bottom line is that this is messy. You're having to manually manage memory *on top of a garbage collector*.
- In addition to the above, you don't have easy access to the elements of the array. The better way to deal with both is to use the provided `Vector` struct. In other words, you can do this

```
auto x = Vector("rnorm(15)");
```

The `Vector` struct uses reference counting to handle the memory management for you, and it allows convenient access to the data, for instance:

```
double x1 = x[1];
x[2] = 4.9;
```

The other way to evaluate R code is with `evalRQ`. This is much simpler, because it does not return anything. One example is printing something to the screen:

```
evalRQ(`print("This was printed by R")`);
```

It's also useful for intermediate results. Maybe you are reading in a big dataset and you only need one variable.

```
evalRQ(`data.raw <- read.csv("file.csv")
auto x = Vector("data.raw[,2]");
```

`evalRQ` was used to read in the data, then `Vector` was used to acquire a pointer to the variable you want to work with in D.
