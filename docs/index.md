# D as a Better R

This project is a collection of the things I've been doing since 2013 that
connect the D and R languages.[^1] There are also some libraries
that can be used in cases where efficiency is critical. The intended audience
is an academic researcher doing the type of data analysis that gets done with
R or Python, but who has a preference to write their program in D, whether for
speed, static typing, or the nice features of the language. Even though the
intended audience is academic researchers doing empirical work, it is likely
to be of interest to anyone doing statistical analysis, and to many others
doing scientific and numerical computing in software such as
Matlab.[^2] 

The emphasis is on functionality and the speed with which you
can write correct code. When this conflicts with performance, I've been
willing to let my programs run for a little bit longer. I want to minimize the
learning curve as much as possible. For instance, one design goal is
that you don't need to know anything about memory allocation,
memory management, or garbage collection. As soon as you introduce those
requirements, you've lost almost the entire community of academic
researchers, and the ones that persevere are likely to get things wrong
and waste a lot of time not doing the analysis they're supposed to be
doing. I'm not targeting C++ or Rust programmers, who will surely be
repulsed by the focus on getting work done and the absence of premature
optimization.

[^1]: I released an earlier project called embedr, and even wrote [a blog 
post about it for the D blog](https://dlang.org/blog/2020/01/27/d-for-data-science-calling-r-from-d/). The latest version of that project [can be found here](https://github.com/bachmeil/embedrv2).
While that project
only made available a fraction of the things I'd done to that point.

[^2]: Although performance is suitable for most uses, it does not and never will
provide the fastest code possible. You'll want to look at libraries such
as Mir if that's your goal, because this project doesn't have much to
offer you on that dimension.

My earlier focus was on writing D functions and calling them from R,
motivated by what Rcpp had done for C++ usage. My [embedrv2](https://github.com/bachmeil/embedrv2)
project makes this easy by using metaprogramming to write the bindings
for you. Over time I've come to the conclusion that this is the wrong
direction for interoperability. I want to write complete programs
in D, not R programs with a few bottlenecks rewritten in D. The approach
I've taken is to treat the full R language as a library called by my D
program. My D programs have full access to everything in base R, all
packages, all libraries with an R interface...literally everything you
get with R is available in D, in a fully seamless fashion that feels like
it's D all the way down.

Now, there are certainly some downsides to this. You have to write wrappers
over R. A few things are inherently slower once you involve R, so you need
to rewrite those pieces in D or find a library in another language that
does what you need. You have to link to shared libraries. I had to
learn how the R internals work in order to build the lowest-level foundation
for interoperability.

# Data Sharing

The starting point for interoperability is data sharing. I
needed a way to create, access, and manipulate R data 
structures from D. I needed to work with data frames, vectors, matrices,
lists, and so on from my D program, but they had to simultaneously
be available to R. This was not a trivial undertaking. 
Three things would be required for this to be practical:

- Allocation of data structures. This is easy. All of R's data structures
are a single C struct under the hood. You can send a string of code to R
and then capture the output, which is a pointer, or you can directly call
functions in libR.so, which do the allocation and return a pointer. The only
differences in the two approaches are that the latter is faster and the
former creates a variable inside R that can be passed as an argument to
R functions.
- Releasing the data structures to the R garbage collector when they're
no longer in use. It's easy to *protect* objects from the garbage collector.
What's more difficult is to recognize when they're no longer needed and
remove the protection exactly once. Failure to unprotect creates a memory leak,
unprotecting too soon results in segfaults, and unprotecting twice kills
your program. I implemented a reference counting approach to handle all 
of these memory management issues. To my
knowledge, all of the bugs have been worked out. It's been a long time
since I had any issues related to memory. I'm currently moving from the
reference counting system I implemented myself to `SafeRefCounted`. If
there are any remaining issues I don't know about, transitioning to
`SafeRefCounted`, which wasn't around when I wrote my reference
counting code, should fix them.
- Convenient access to the data from D. Working out the allocation
and garbage collection for matrix isn't terribly helpful if you still
need to work with individual elements using pointer arithmetic. It's 
only sustainable if you have syntax like `m[1,4] = 3.2`, `m[2..6, 0..3] = 4.0`, and
`m[1,3..$]`. Once you have that, you can build on it with things like
`Row(m,2) = [3.6, -1.1, -2.4, 4.8]`. These conveniences 
don't take much time individually, but there are many of
them, and comprehensively testing all of them and fixing bugs is a slow 
process for a side project.

# Running R Code and Calling R Functions

Once I had the data sharing under control, I needed a way to tell R what
to do. The simplest version of this is passing a string of R code to the
R shared library, having it evaluate the code, and returning a pointer to
the output (if any). In other cases, you're calling the C functions that
R calls under the hood, where you pass R data structures as arguments
and receive an R data structure as the output. Finally, you can call
R packages that provide an interface to C, C++, or Fortran code directly.
This is the same as calling C functions, but you have to find and link to
the shared library for that package, as every R package with compiled code
has its own shared library.

# Bottlenecks

Although your program is probably going to be more than fast enough out
of the box, [as I wrote about here](efficiency.html), there will be cases
in which it's worth your while to speed things up. The leading example
is a simulation that's repeated many times. Parallel random number generation
and optimized linear algebra libraries are an essential part of the
toolbox for simulation. I've done three things to facilitate this:

- Ported one of L'Ecuyer's parallel random number generators from Java to D.
- Stripped out the random number generators for most of the distributions
in the GNU Scientific Library and polish them so they can be compiled with
ImportC.
- Stripped the matrix library in Gretl into its a standalone library and
got it compiling with ImportC.

Other optimizations have been implemented. I recommend reading [the discussion here](efficiency.html)
for more on this topic. If worst comes to worst, you always have the options
to rewrite bottlenecks in D or to call any available library with a C
interface. I want to emphasize that this should not usually be necessary.
You should be able to write your program and not worry about the speed.

# Documentation

This is always tough. It takes a long time to write good documentation.

# I'm Sharing My Work

Please note that I view this as *sharing my work* rather than *releasing a library*. 
I'm showing you what has worked for me, and giving you the details of my 
workflow, but it's not polished like a well-run open source project. 
Maybe others will see the value and pitch in. If not, well, I'll 
continue to use it and update it as I do. I'm largely
indifferent on adoption, but I would be happy if someone else found it
valuable, even if I'm not going to go out of my way to make that happen.

# How Does This Compare To Project X?

There have been many previous efforts on this front. There have been
statistical libraries and optimization libraries and plotting libraries
and so on. Mir is an impressive project with some great pieces.

The problem with that approach is that there will never be a pure D solution that's remotely suitable for the day-to-day work of a data analyst. It's simply too big of a task to write everything from scratch in D, and honestly, it's pointless to do so. Reusing code written in other languages was never a problem for R, Python, or any other community, so why not do the same thing for D?

This project intends to provide a complete solution for the data analyst. One of the goals is a convenient, efficient D solution for *any* data analysis you'd want to do. Once you discard the silly objections to calling other languages that are prevalent in the D community (seriously, D even compiles C code!) it's not that big of a task. D can give you access to

- Everything written in D, C, and Fortran.
- Everything written in/for R (you can treat R as a shared library and call any of its packages).
- Everything written in Python and Julia (you can reuse the bridges that have been created to work with R).
- Anything written in C++ that has an R interface, which includes thousands of packages, and since R interfaces are just C interfaces, you can actually reuse that infrastructure to call most C++ libraries.

There are limitations, perhaps, but it would be silly to argue D is not a viable choice for the vast majority of data analysis being done in 2024.

# Better R?

Walter Bright has referred to a particular use case for D as
"better C", sometimes jokingly calling it [DasBetterC](https://dlang.org/blog/2018/06/11/dasbetterc-converting-make-c-to-d/).
Although the use case is somewhat different, the name "better R" is a
reasonable description of what I'm doing with this project:

- All of R is available to a D program (you can run arbitrary R code,
capture the output, and access it from D, most of the time without any
copying).
- You can rewrite R programs in D, sprinkling in D features where they
add convenience or speed.
- You can modify the behavior of R where you don't like it. One of the
things I wish R did differently is dropping dimensions. A row of a matrix
will silently become a vector, with a completely different interpretation.
That doesn't happen when I access the matrix from D.
- You add D's static typing to R. You can't imagine how good that is
given that R is not only a dynamic language, but that it was developed
at the same place that gave us C. (Both were created at Bell Labs. 
The first release of S was four years after the first release of C.
Both of them having single-letter names was not an accident.)
- You can use D's compile time features on top of a vanilla R program.

# Getting Started

Probably the best way to figure out what is going on is to look at the
example and then look at the tests linked in the next section. After that,
you can read the other documentation. You can also start a discussion on
the Github repo if you have questions or wonder how to do something.

In the meantime, here's a short example to give you an idea of what to expect:

```
void main() {
  // Initialization
  startR();
  // Generate 100 draws from a N(0.5, sd=0.1) distribution
  Vector x = rnorm(100, 0.5, 0.1);
  // Create a 3x2 matrix and fill the columns
  auto m = Matrix(3,2);
  Column(m,0) = [1.1, 2.2, 3.3];
  Column(m,1) = [4.4, 5.5, 6.6];
  // Calculate the inverse of the transpose
  auto m2 = solve(t(m));
  // Modify the second and third elements of the first column of m
  m[1..$,0] = [7.7, 8.8];
  // Choose x and y to minimize x^2 + y^2
  // Use Nelder-Mead with initial guesses 3.5 and -5.5
  auto nm = NelderMead(&f);
  OptimSolution sol = nm.solve([3.5, -5.5]);
  // Clean up
  closeR();
}

extern(C) {
  double f(int n, double * par, void * ex) {
    return par[0]*par[0] + par[1]*par[1];
  }
}
```

# Usage

- [Installation](installation.html)
- [Compiling](compiling.html)
- [Example](example.html)
- [Tests demonstrating most of the functionality](https://github.com/bachmeil/betterr/tree/main/testing)

# Modules

- [betterr.baser](base.html)
- [betterr.matrix](matrix.html)
- [betterr.vector](vector.html)
- [betterr.random](random.html)
- [betterr.ts](ts.html)
- [betterr.optim](optim.html)
- [betterr.plot](plot.html)
- [betterr.array](array.html)
- [betterr.dataframe](dataframe.html)
- [betterr.list](list.html)
- [betterr.lm](lm.html)
- [Pieces of the R API](api.html)
- [Quadratic programming](quadprog.html)

# Notes on Various Topics

- [Using betterr with WSL](wsl.html)
- [Evaluating arbitrary R code](evalr.html)
- [Accessing databases](databases.html)
- [Generating scalar random variables efficiently](randomscalar.html)
- [Parallel random number generation](prng.html)
- [Running Better R programs in parallel](parallelrun.html)
- [Using OpenBLAS for matrix calculations](openblas.html)
- [Calling the Matrix package from D](matrixpackage.html)
- [Passing an existing variable from D into R](setvar.html)
- [R lazy copying](lazycopy.html)
- [Understanding argument lists in the R source](arglists.html)
- [Passing function pointers around](funcptr.html)
- [Using the GSL to generate random numbers sequentially or in parallel](gslrng.html)
- [Options for numerical programming in D](numerical.html)
- [Compiling With ImportC](importc.html)

# Thoughts

- [Efficiency](efficiency.html)
- [Scaling to Zero](scaling.html)
- [Some Ways To Improve This Project](improve.html)
