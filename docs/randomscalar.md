# Taking Scalar Random Draws Efficiently

One application where efficiency is critical involves taking many draws
from a distribution. When simulating from the posterior distribution for
Bayesian inference, for example, it's not unusual to take billions of
draws. Ideally, you will be calling D, C, or Fortran functions to take
those draws.

There are variety of readily available options to do that. Indeed, you
can directly call the C functions in libR, you can use Mir for a D
solution, or you can call a C library like the GSL. This is relatively
straightforward if you are taking draws serially on a single processor.
In 2023, you really need support for parallel random number generation,
due to the fact that all modern desktop computers have multicore 
processors. Parallel random number generation is harder to get correct.
The easiest thing to do is use the builtin parallel random number
generation functionality that comes with a default installation of R.

The downside of R for scalar random number generation is the
considerable overhead associated with taking a single draw from a
distribution. A better strategy is to take k >> 1 draws (which will be
done by a C function, store them in a vector, and treat the elements of 
the vector as draws one-by-one.

The `Generator` struct exists for this purpose. It tracks the current 
index of the storage vector, and each time you take a draw, the index is 
incremented. Once you've exhausted all draws, it refills the Vector and 
resets the index to zero. There is limited overhead associated with
using a `Generator` relative to a D/C/Fortran library, due to the
is the incrementing of the index and checking if the Vector is
empty. I believe the easy access to random number generation from a
wide variety of distributions coupled with straightforward support for
parallel generation makes up for the overhead. I may at some point
optimize this in some way, but for now I have little incentive.

The number of draws stored in the Vector is by default set to 1000,
but that can be changed. Generators are provided for all of the 
distributions with random number generators in base R, and it's 
straightforward to set up a custom generator for any distribution.

To generate 50,000 draws that modify a standard normal, you can do 
something like this (I'm aware this isn't something you'd do in
practice):

```
auto v = Vector(50_000);
Generator!"norm" normGen;
foreach(ii; 0..50_000) {
  v[ii] = normGen.draw()/3.0;
}
```

To take only 100 draws at a time, you'd adjust the compile time parameters
passed to the constructor:

```
auto v = Vector(50_000);
Generator!("norm", 100) normGen;
foreach(ii; 0..50_000) {
  v[ii] = normGen.draw()/3.0;
}
```

If you want to change the parameters of the distribution:

```
auto v = Vector(50_000);
Generator!"norm" normGen;
normGen.mean = -4.6;
normGen.sd = 2.1;
foreach(ii; 0..50_000) {
  v[ii] = normGen.draw()/3.0;
}
```

This is how you'd set up a custom Generator for the previous example:

```
auto v = Vector(50_000);
Generator!("custom", 1000, double) customGen;
customGen.cmd = "rnorm(1000, mean=-4.6, sd=2.1)";
foreach(ii; 0..50_000) {
  v[ii] = customGen.draw()/3.0;
}
```

For a custom generator, note the following:

- There is a third compile time parameter denoting the type of the
scalar returned for each draw. It can be double, int, or bool.
- You set `cmd` equal to a string holding the R call that generates
the draws each time the vector is refilled.
- The second compile time parameter is ignored. `cmd` specifies the
number of draws to take on a refill. It's a bit clumsy to pass an unused
paramter, but doing it this way keeps the library code simple.
