# Some Options for Numerical Computing in D

One of the common criticisms of D is that it does not have the ecosystem to get your work done. In my opinion, that is inaccurate if you want to the type of computing I do, related to statistics and econometrics. This is a discussion of some of the options I've used over the years. Yes, it involves calling into C libraries, but that doesn't seem to surface as a criticism in discussions of Python, R, Julia, Matlab, and basically every statistics program in heavy use.

Calling D from R. Originally, this was what I thought was the approach to take. R is popular. It's easy to compile D code into a shared library and call it from R. My embedrv2 project makes this really easy by writing the wrapper for you.

Why did I move away from this? There are a couple of reasons. First, if I prefer writing D code, why should I use an approach that's mostly R with a little D sprinkled in? Second, you're using D in a very limited sense by merely rewriting bottlenecks for speed. This approach doesn't do anything to provide numerical libraries for D.

Calling R from D. I actually used this approach first. As in, back in the summer of 2013 when I used D for the first time in a serious capacity, I wrote a program that called into R. Over time I've come to the realization that this is the correct approach to take. Dirk was kind enough to add my small C API to RInside. That turns R into a C shared library that can be called by any language that has C interoperability. With a very minimal amount of effort, you suddenly have all of R, both the base language and every package, available to any D program. Any C, C++ or Fortran library with an R interface now has a D interface. Note that these interfaces are efficient - the R interface is actually compiled C code that you call from D. You have database interfaces, you can use all of Tensorflow, there's Stan for Bayesian inference, and on and on. While I haven't written convenient wrappers for everything, which will obviously never happen, it's quite easy to do so yourself. The low-level interoperability is all handled by betterr.

Calling Gretl. This library may not be known outside of economics, but it contains a lot of things that are of interest for data analysis. It's written in C. One of the more appealing things about Gretl is that it contains an easy-to-use matrix library that sits on top of BLAS and LAPACK. I realeased a package years ago. I'm updating it now, so it's not currently online. Gretl is a mature, well-tested library that's been in development for more than 20 years.

Calling the GNU Scientific Library (GSL). This is a well-known library in numerical computing circles. One way I use it is for parallel random number generation. I ported one of Lecuyer's functions to D, and I ported most of the RNG distributions to D, so I can trivially add PRNG to my programs. This is, of course, an important thing for simulation in 2024, as 32 or more cores on a desktop computer becomes the norm.

Calling Julia from D. I haven't actually done anything with this, but a while back Symmetry made something available. If there is something I need from Julia I will be able to call Julia from D. I simply have not had a need for it to this point.

dstats. This is a pure D library that was written many years ago by David Simcha. It includes a range of functionality, some of which has made its way into betterr.

Mir. This is an impressive body of work. I can't say I've used it a lot, but I initially had it included in betterr for speedy low-level operations.

BLAS and LAPACK. You can always fall back on these if you need speedy matrix operations. There are D wrappers available on code.dlang.org.

ggplotd. I've never used it, but it's motivated by ggplot and adds plotting capabilities.

scid. I haven't looked at this recently, but it includes numerical integrationa and differentiation, special functions, linear algebra, and nonlinear equation solvers.

tsv-utils. A little different from the other packages above, but provides tools to work with large text data files. Originally written for use at EBay.

CUDA. I've never done anything with this, but there's a project out there that is claimed to be a very nice way to work with GPUs. That's all I can say about it. The [repo is here](https://github.com/libmir/dcompute) and there were a [couple of articles](https://dlang.org/blog/2017/07/17/dcompute-gpgpu-with-native-d-for-opencl-and-cuda/) on [the D blog](https://dlang.org/blog/2017/10/30/d-compute-running-d-on-the-gpu/).
