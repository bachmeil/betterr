# D as a Better R

This project is a collection of things I've been doing over the last
decade to connect R and D. Although I released an earlier project called
embedr, even writing a blog post about it for the D blog, that project
only made available a fraction of the things I'd done to that point.

In addition, while the interoperability went in both directions (D calling
R functions and R calling D functions) the focus was primarily on support
for writing D functions that would be called from R. I don't have an
opinion on whether that was the right thing to do, but over time I have
come to the conclusion that I want to write my programs in D, but without
having to give up the functionality of R.

The starting point for interoperability is data sharing. That meant I
needed to have a way to create, access, and manipulate the main R data 
structures from D and R. I'd need to work with data frames, vectors, matrices,
lists, and so on from my D program, but importantly, I needed them to
be available to R at the same time. This was not a trivial undertaking. 
You need three things for this to be practical:

- Allocation of the data structures. This is fairly easy. You can send
some code to R and capture a pointer to the data, or you can call
functions in libR.so, the shared library written in C and Fortran that 
R uses internally.
- Releasing the data structures to the R garbage collector when they're
no longer in use. It's easy to *protect* objects from the garbage collector.
What's more difficult is to recognize when they're no longer needed and
remove the protection. Failure to unprotect creates a memory leak, and
unprotecting too soon results in segfaults. I've implemented a reference
counting approach to handle all of the memory management issues. To my
knowledge, all of the bugs have been worked out. It's been a long time
since I had any issues related to memory.
- Convenient access to the data from within D. Working out the allocation
and garbage collection for matrix isn't terribly helpful if you still
need to work with individual elements using pointer arithmetic. It's 
only sustainable if you have syntax like `m[1,4] = 3.2`, `m[2..6, 0..3] = 4.0`, and
`m[1,3..$]`. Once you have that, you can build on it with things like
`Row(m,2) = [3.6, -1.1, -2.4, 4.8]`. Any one of these conveniences 
doesn't take much time to implement on its own, but there are many of
them, and comprehensively testing them all is a slow process for a side
project.

Please note that I view this as *sharing my work* rather than *releasing a library*. 
I'm showing you what has worked for me, giving you the details of my 
workflow, but it's not polished like a well-run open source project. 
Maybe others will see the value and pitch in. If not, well, I'll 
continue to use it and update it as I do. I'm largely
indifferent on adoption, but I would be happy someone else found it
valuable.

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

# Thoughts

- [Efficiency](efficiency.html)
- [Scaling to Zero](scaling.html)
- [Some Ways To Improve This Project](improve.html)
