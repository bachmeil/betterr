# Overview

Usage:

```
import betterr.baser;
```

This module exposes some base R functionality to D.

# Source Code and Examples

- [Source Code](https://github.com/bachmeil/betterr/blob/main/baser.d)
- [Examples](https://github.com/bachmeil/betterr/blob/main/testing/testbase.d)

# Functions

These functions take a single argument, and function similar to the underlying R functions. Links are to the R documentation, which provides everything you need to know. Where it applies, you can call the functions with D scalars or arrays.

- [abs](https://web.mit.edu/r/current/lib/R/library/base/html/MathFun.html)
- [sqrt](https://web.mit.edu/r/current/lib/R/library/base/html/MathFun.html)
- [ceiling](https://web.mit.edu/r/current/lib/R/library/base/html/Round.html)
- [floor](https://web.mit.edu/r/current/lib/R/library/base/html/Round.html)
- [trunc](https://web.mit.edu/r/current/lib/R/library/base/html/Round.html)
- [cosh](https://web.mit.edu/r/current/lib/R/library/base/html/Hyperbolic.html)
- [sinh](https://web.mit.edu/r/current/lib/R/library/base/html/Hyperbolic.html)
- [tanh](https://web.mit.edu/r/current/lib/R/library/base/html/Hyperbolic.html)
- [acosh](https://web.mit.edu/r/current/lib/R/library/base/html/Hyperbolic.html)
- [asinh](https://web.mit.edu/r/current/lib/R/library/base/html/Hyperbolic.html)
- [atanh](https://web.mit.edu/r/current/lib/R/library/base/html/Hyperbolic.html)
- [log10](https://web.mit.edu/r/current/lib/R/library/base/html/Log.html)
- [log2](https://web.mit.edu/r/current/lib/R/library/base/html/Log.html)
- [log1p](https://web.mit.edu/r/current/lib/R/library/base/html/Log.html)
- [exp](https://web.mit.edu/r/current/lib/R/library/base/html/Log.html)
- [expm1](https://web.mit.edu/r/current/lib/R/library/base/html/Log.html)
- [cumsum](https://web.mit.edu/r/current/lib/R/library/base/html/cumsum.html)
- [cumprod](https://web.mit.edu/r/current/lib/R/library/base/html/cumsum.html)
- [cummax](https://web.mit.edu/r/current/lib/R/library/base/html/cumsum.html)
- [cummin](https://web.mit.edu/r/current/lib/R/library/base/html/cumsum.html)
- [rev](https://web.mit.edu/r/current/lib/R/library/base/html/rev.html)
- [summary](https://web.mit.edu/r/current/lib/R/library/base/html/summary.html)

These functions are the same as those above, but take an optional second argument. If set to true, NA values are removed before the calculation is done.

- [max](https://web.mit.edu/r/current/lib/R/library/base/html/Extremes.html)
- [min](https://web.mit.edu/r/current/lib/R/library/base/html/Extremes.html)
- [sum](https://web.mit.edu/r/current/lib/R/library/base/html/Extremes.html)

These functions take optional arguments. At this time, the optional arguments are a string (with comma separating arguments) that is passed directly to R. See the R documentation for the available optional arguments.

- [round](https://web.mit.edu/r/current/lib/R/library/base/html/Round.html)
- [signif](https://web.mit.edu/r/current/lib/R/library/base/html/Round.html)
- [log](https://web.mit.edu/r/current/lib/R/library/base/html/Log.html)
- [logb](https://web.mit.edu/r/current/lib/R/library/base/html/Log.html)
- [prod](https://web.mit.edu/r/current/lib/R/library/base/html/prod.html)
- [range](https://web.mit.edu/r/current/lib/R/library/base/html/range.html)
- [rank](https://web.mit.edu/r/current/lib/R/library/base/html/rank.html)
- [sort](https://web.mit.edu/r/current/lib/R/library/base/html/sort.html)

Finally, we have

`seq(double from, double including, double by)`

- `from`: Starting value of the sequence
- `including`: Ending value of the sequence - following R, included in the sequence
- `by`: Increment

Note that these do not need to be integers, so this function is different from R's `:` operator.

[R documentation](https://web.mit.edu/r/current/lib/R/library/base/html/seq.html)

# Fast Fundamental Operations

Although the primary goal of this project is to provide all of R's functionality to D programs, relying on R for fundamental building blocks such as summation or calculating the mean will lead to a serious performance hit. As an example, look at what happens if we use the R interpreter to calculate the mean of a `double[]`:

- Allocate a new Vector.
- Copy the elements into the Vector.
- Tell R to calculate the mean of the Vector.
- Get a pointer to the solution.
- Convert the value at the pointer to a double.

That's not a problem if you're going to calculate the mean of a `double[100]` one time. If you're going to calculate the mean of a `double[]` holding tens of thousands of elements on the inside of a loop that's executed a billion times, the above is ridiculously inefficient. To the extent possible, we want to make sure this type of operation is performant, so what actually happens under the hood is that the mean is calculated by Mir. The following operations currently have efficient implementations using Mir or Phobos:

- `double mean(double[] x)`
- `double mean(T)(T v)` for T a Vector or Matrix
- `double mean(double[] x, bool narm=false)`
- `double mean(T)(T x, bool narm=false)` for T a Vector or Matrix
- `double median(double[] x)`
- `double median(T)(T v)`
- `double median(string filtered, T)(T v)` for T a Vector or Matrix
- `double median(double[] x, bool narm=false)`
- `double median(T)(T x, bool narm=false)` for T a Vector or Matrix
- `double sum(double[] x)`
- `double sum(T)(T v)` for T a Vector or Matrix
- `double sum(double[] x, bool narm=false)`
- `double sum(T)(T x, bool narm=false)` for T a Vector or Matrix
- `double max(double[] x)`
- `double max(T)(T x)` for T a Vector or Matrix
- `double max(double[] x, bool narm=false)`
- `double max(T)(T x, bool narm=false)` for T a Vector or Matrix
- `double min(double[] x)`
- `double min(T)(T x)` for T a Vector or Matrix
- `double min(double[] x, bool narm=false)`
- `double min(T)(T x, bool narm=false)` for T a Vector or Matrix
- `double sd(double[] x)`
- `double sd(T)(T x)` for T a Vector or Matrix
- `double sd(double[] x, bool narm=false)`
- `double sd(T)(T x, bool narm=false)` for T a Vector or Matrix
- `double var(double[] x)`
- `double var(T)(T x)` for T a Vector or Matrix
- `double var(double[] x, bool narm=false)`
- `double var(T)(T x, bool narm=false)` for T a Vector or Matrix
- `double quantile(double[] x, double p)` 
- `double quantile(double[] x, double p, bool narm)`
- `double quantile(T)(T x, double p)` for T a Vector or Matrix
- `double quantile(T)(T x, double p, bool narm)` for T a Vector or Matrix
- `double[] quantile(double[] x, double[] p)`
- `double[] quantile(T)(T x, double[] p)` for T a Vector or Matrix
- `double[] quantile(double[] x, double[] p, bool narm)`
- `double[] quantile(T)(T x, double[] p, bool narm)` for T a Vector or Matrix

The following make use of R due to lack of an alternative implementation:

- `double mean(double[] x, double trim=0, bool narm=false)`
- `double mean(T)(T v, double trim=0, bool narm=false)` for T a Vector or Matrix

