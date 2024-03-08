# Generating random numbers sequentially and in parallel using the GSL

Elsewhere I discuss strategies for drawing random numbers at suitable speed using only R. ([Drawing vectors of random numbers](random.html)), ([taking scalar draws](randomscalar.html)), and ([parallel RNG](prng.html)). Those will probably be fast enough for most uses. However, in some cases, you'll want to be fully efficient. A simulation can take so many draws that you may wish to get rid of the overhead of R altogether by using the GSL for RNG.

ImportC allows D to compile C code, which allows you to call those functions directly from your D program, with no bindings or any other effort. The file gslheaders.c contains the C code needed to call GSL. (You'll still need to install GSL. On Ubuntu, it's package libgsl-dev. Maybe in the future we can eliminate the need for the dynamic library, but not yet.)

The file prng.d is a port of parallel RNG code provided by Pierre L'Ecuyer. You can see an example of how it works in the file testgslprng.d. At some point I may add documentation, but there's not much to add. The GSL functions [are documented here](https://www.gnu.org/software/gsl/doc/html/randist.html).
