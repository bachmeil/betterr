# Running Programs in Parallel

The availability of multiprocessor desktop computers has led to the
heavy use of parallel programming techniques for numerical work. One
downside of embedding R is that libR is not designed to be used from a
parallel program. You can forget about using this library with 
std.parallelism. Your program will crash and burn if you try. One
reason is that R uses a global stack to track allocations. If you try to
call multiple functions at the same time, and they're all tracking the
same set of allocations, your program will die with a notification that
there's a stack imbalance.

That doesn't mean you can't use this library in a parallelized program.
It just means you have to do one of the following:

- Do all the parallel work inside R. If you're doing only random number
generation in parallel, you might be able to do it with an R function
that is called from a single D process.
- Use std.parallelism on parts that do not require libR. An example
would be generating all of the random numbers for a simulation using
multiple cores in R, and then taking that data and processing it in
parallel using D with std.parallelism.
- Using GNU parallel to run multiple D programs.

The latter is my preferred approach. Since it's easy to work with
multiple random number streams simultaneously, all you need to do is
write a D program that takes the stream number as an argument.

An example Makefile line taken from a real-world program is

```
parallel ./simulation {} $(shell grep '^core id' /proc/cpuinfo |sort -u|wc -l) ::: $(shell seq 1 "$(shell grep '^core id' /proc/cpuinfo |sort -u|wc -l)")
```

That's a bit of a mess, but the ugliness comes from automatic detection of
the number of available cores. If I want to hard code it to four cores,
I can rewrite that as

```
parallel ./simulation {} 4 ::: $(shell seq 1 4)
```

That calls `simulation` with the first argument set to 4 (the total number
of cores) and the second equal to 1, 2, 3, or 4 (this process number).
Those numbers are used within my program to set up the parallel RNG and
appropriately divide the work.

This is a simple approach to embarrassingly parallel tasks such as simulations.

## Example

Here's a quick example (omitting some imports) to show how you can run a betterr program in
parallel that uses R's built-in parallel random number generator.

```
import betterr;

void main(string[] args) {
  auto thisProcess = args[1].to!int;
  auto nprocs = args[2].to!int;
  startR();
  
  // This sets the process to the value passed at the command line
  // 500 is the seed, which is then the same for all instances of your
  // program. The seed is 1 if you don't specify it.
  prngInit(thisProcess, 500);
  long reps = ceiling(1000/nprocs);
  // Do whatever you need to do in your simulation
  // Save the output to a file or do whatever you'd usually do with
  // GNU parallel output.
  closeR();
}
```

To run this on eight cores, you'd compile it and then run it like this:

```
parallel ./example {} 8 ::: $(shell seq 1 8)
```





