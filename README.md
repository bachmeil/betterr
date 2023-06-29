# D as a Better R

I've made this repository available to share with others how I use D as a "better R". You'll find libraries, examples, documentation, and posts. You can think of the library as enabling R as a shared library. Your D program has access to all of R's data structures and functions.

The name "better R" is loosely inspired by [Better C](https://dlang.org/spec/betterc.html). The idea is that R is fully integrated into your D program, giving you access to everything in R, and making it possible to pass data between the two languages, but allowing you the option to improve upon the way R does things as much as you'd like. In some sense, R is a proper subset of D.

Feel free to use it, but things might break. There's no intention of turning this into a Dub package; I don't have time to maintain such a thing. Since it's open source, anyone that wants to is free to do so.

# Current Status (June 2023)

I've implemented support for the following data structures and types:

- Vector
- Matrix
- Data frame
- List
- Scalar double, int, bool, and string - R does not have a scalar type, but this library makes it appear that it does
- Arrays
- Time series (single)

And the following functionality:

- OLS regression
- Random number generation and sampling with or without replacement
- [dqp] functions for working with distributions
- Numerous functions such as `abs` and `acosh`.
- Basic plotting
- Reading of data files

# Documentation

There's some documentation [at the website](https://bachmeil.github.io/betterr). Probably more useful is to look at the examples in the [testing directory](https://github.com/bachmeil/betterr/tree/main/testing). For the most part, I've tried to maintain the behavior of R, so you can usually read the R documentation.

# Example

I'm running Ubuntu 22.04. I have a program `testvector.d`:

```
import callr.vector;
import callr.r;
import std.conv, std.stdio;

void main() {
  startR();
  
  auto x = Vector(5);
  x[0] = 1.5;
  x[3] = 4.1;
  x.print("Test vector");
  
  auto y = Vector("rnorm(25)");
  y.print();
  
  auto y2 = Vector(y.length);
  foreach(ii; 0..y.length) {
    y2[ii] = y[ii]*y[ii];
  }
  y2.print();
  y2[4..10].print();
  y2[].print();
  y2[[0,2,4,6,8]].print();
  
  writeln(y2.to!(double[]));
  
  closeR();
}
```

Explanations of some lines:

```
startR();
```

This initializes the R interpreter. If you don't call this, you'll see lots of stuff is `null` and get lots of segfaults. Always check that you have this line in your program if you can't figure out why you're getting strange results. Note that your program **will** compile and run without it. You just won't be getting usable output.

```
auto x = Vector(5);
x[0] = 1.5;
x[3] = 4.1;
x.print("Test vector");
```

`x` is a vector with five elements. On creation, each element is initialized to 0.0. `Vector` is a D struct, and can be passed around and manipulated like any other D struct. Internally, `x` holds information about the data held inside R, but you should never need to worry about that. From the perspective of the user of callr, you should rarely need to know anything about the internals.

Two of the elements are changed to non-zero values and `x` is printed. You can optionally include a string with a message to be printed before printing out the data.

```
auto y = Vector("rnorm(25)");
y.print();
```

You can capture the output of any R expression for use inside your D program. In this case, R generates 25 draws from a standard normal distribution. In the future, a wrapper will be written so you can do `Vector y = rnorm(25)`, or you can create a struct `RNG` with configuration options.

```
auto y2 = Vector(y.length);
foreach(ii; 0..y.length) {
	y2[ii] = y[ii]*y[ii];
}
y2.print();
```

This shows how to get and manipulate elements of a `Vector` from D. This is done efficiently. No R code is generated for this; it is using a pointer to work directly with the elements of the underlying data arrays for `y` and `y2`.

```
y2[4..10].print();
y2[].print();
y2[[0,2,4,6,8]].print();
```

You can take a slice of a vector. The first case is the usual D slice notation. Note that `Vector` indexing starts with 0 and the second index is exclusive. This is because we are adding functionality to D programs, *not* translating D code to R. In R, `y[4:10]` has seven elements, the first being element four and the last being element 10. `y[4..10]` has six elements, the first being element 5 and the last being element 10.

`y2[]` takes a slice of the entire vector `y2`. This leads to copying of the `y2` from R.

`y2[[0,2,4,6,8]].print();` prints the every other element of `y2` from the first to the ninth.

```
writeln(y2.to!(double[]));
```

You can use `std.conv.to` to convert a `Vector` to a `double[]`. This involves copying. In general, there is no reason to do that, but you might want to pass the data to an existing D function that takes a `double[]` as an argument.

```
closeR();
```

This does any needed cleanup.

# Requirements

OS: Development is done on Linux. Everything has been confirmed to work on WSL. It shouldn't take much to get it working on Mac and Windows, but since that's not my area of expertise, I can't help you (contributions welcome!)

Software: You need a relatively recent version of R and you need to install the R package RInside.

# Installation

I do my compilation with the -i flag rather than using Dub. This is what's easiest for me, but there's no reason you can't use Dub or some other approach.

## Initial Setup

1. Clone the betterr repo.
2. Change the variable REPO in the first line of the Makefile to the directory holding the betterr repo.
3. Set LINK to include the correct path of libRInside.so. If you're not sure how to find it, [instructions are here](https://bachmeil.github.io/betterr/compiling.html).

## Adding to a project

1. Copy the Makefile into your project's directory.
2. Run `make dep`.

# Compiling

Add a target to the Makefile. It will look something like this (see the testing directory for many examples):

```
ldmd2 -i program.d $(LINK)
```

# Speed

You might expect this to be slow. After all, you're running interpreted R code, which defeats the purpose of using D. The short answer is that it's unlikely to be slow and it doesn't defeat the purpose of using D. The long answer [can be found here](https://bachmeil.github.io/betterr/efficiency.html).









