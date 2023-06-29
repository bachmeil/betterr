# Overview

Usage:

```
import betterr.random;
```

This module exposes all of the random number generation function in base R to a D program.

# Source Code and Examples

- [Source Code](https://github.com/bachmeil/betterr/blob/main/random.d)
- [Examples](https://github.com/bachmeil/betterr/blob/main/testing/testrandom.d)

# Structs

## Sample

```
struct Sample {
	Vector x;
	long size;
	bool replace = false;
	Vector prob;
}
```

### Data

- `x`: Vector holding the data to be sampled.
- `size`: Number of elements to sample from x.
- `replace`: True if sampling from x is done with replacement.
- `prob`: A vector specifying the probabilities of drawing each of the elements of `x`.

### Methods

- `draw()`: Return a vector holding the sample drawn based on the 
configuration you have specified.

See also: [R documentation](https://web.mit.edu/r/current/lib/R/library/base/html/sample.html)

# Functions

## rnorm

`Vector rnorm(long n, double mean=0, double sd=1)`

- `n`: Number of draws
- `mean`
- `sd`: Standard deviation

See also: [rnorm](https://web.mit.edu/r/current/lib/R/library/stats/html/Normal.html)

## runif

`Vector runif(long n, double min=0, double max=1)`

- `n`: Number of draws
- `min`
- `max`

See also: [runif](https://web.mit.edu/r/current/lib/R/library/stats/html/Uniform.html)

## rgamma

`Vector rgamma(long n, double shape, double rate=double.nan, double scale=1.0)`

R has named optional parameters. D does not. Set rate to double.nan to leave it unspecified.

- `n`: Number of draws
- `shape`
- `rate`
- `scale`

See also: [rgamma](https://web.mit.edu/r/current/lib/R/library/stats/html/GammaDist.html)

## rbeta

`Vector rbeta(long n, double shape1, double shape2, double ncp=0.0)`

- `n`: Number of draws
- `shape1`
- `shape2`
- `ncp`: Non-centrality parameter

See also: [rbeta](https://web.mit.edu/r/current/lib/R/library/stats/html/Beta.html)

## rbinom

`Vector rbinom(long n, long size, double prob)`

- `n`: Number of draws
- `size`: Number of trials
- `prob`: Probability of success on each trial

See also: [rbinom](https://web.mit.edu/r/current/lib/R/library/stats/html/Binomial.html)

## rcauchy

`Vector rcauchy(long n, double location=0, double scale=1)`

- `n`: Number of draws
- `location`
- `scale`

See also: [rcauchy](https://web.mit.edu/r/current/lib/R/library/stats/html/Cauchy.html)

## rchisq

`Vector rchisq(long n, double df, double ncp)`

- `n`: Number of draws
- `df`: Degrees of freedom (>= 0)
- `ncp`: Non-centrality parameter (>= 0)

See also: [rchisq](https://web.mit.edu/r/current/lib/R/library/stats/html/Chisquare.html)

## rexp

`Vector rexp(long n, double rate=1.0)`

- `n`: Number of draws
- `rate`: Mean is 1/rate

See also: [rexp](https://web.mit.edu/r/current/lib/R/library/stats/html/Exponential.html)

## rf

`Vector rf(long n, long df1, long df2, double ncp)`  
`Vector rf(long n, long df1, long df2)`

- `n`: Number of draws
- `df1`: First degrees of freedom parameter
- `df2`: Second degrees of freedom parameter
- `ncp`: Non-centrality parameter

See also: [rf](https://web.mit.edu/r/current/lib/R/library/stats/html/Fdist.html)

## rgeom

`Vector rgeom(T)(long n, T prob)`

- `n`: Number of draws
- `prob`: A Vector of probabilities or anything that converts to a Vector, including a double[].

See also: [rgeom](https://web.mit.edu/r/current/lib/R/library/stats/html/Geometric.html)

## rhyper

`Vector rhyper(long nn, long m, long n, long k)`  
`Vector rhyper(Vector nn, long m, long n, long k)`  
`Vector rhyper(long[] nn, long m, long n, long k)`

- `nn`: Number of draws (observations). If it is Vector or long[], then the number of draws is the number of elements of `nn`. This is consistent with the behavior of the underlying R function.
- `m`: Number of white balls
- `n`: Number of black balls
- `k`: Number of balls that are drawn from the urn to generate each observation

See also: [rhyper](https://web.mit.edu/r/current/lib/R/library/stats/html/Hypergeometric.html)

## rlnorm

`Vector rlnorm(long n, double meanlog=0, double sdlog=1)`

- `n`: Number of draws
- `meanlog`: Mean of the distribution on the logscale
- `sdlog`: Standard deviation of the distribution on the logscale

See also: [rlnorm](https://web.mit.edu/r/current/lib/R/library/stats/html/Lognormal.html)

## rmultinom

`Matrix rmultinom(T)(long n, long size, T prob)`

- `n`: Number of random vectors to draw
- `size`: Total number of elements drawn
- `prob`: Vector of probabilities of each of the elements that are drawn. The number of elements in each of the `n` vectors is `size/prob.rows`. `prob` can be a Vector, or anything that converts to a Vector.

See also: [rmultinom](https://web.mit.edu/r/current/lib/R/library/stats/html/Multinom.html)

## rnbinom

`Vector rnbinom(long n, double size, double prob, double mu)`  
`Vector rnbinom(Vector n, double size, double prob, double mu)`  
`Vector rnbinom(T)(T[] n, double size, double prob, double mu)`  
`Vector rnbinom(long n, double size, double prob)`  
`Vector rnbinom(Vector n, double size, double prob)`  
`Vector rnbinom(T)(T[] n, double size, double prob)`

- `n`: Number of draws
- `size`: Target for number of successful trials
- `prob`: Probability of success in each trial
- `mu`: Alternative to specifying prob

If you want to specify `mu` rather than `prob`, set `prob` to double.nan. If `n` is a Vector or converts to a Vector, then the number of draws is `n.rows`.

See also: [rnbinom](https://web.mit.edu/r/current/lib/R/library/stats/html/NegBinomial.html)

## rpois

`Vector rpois(double n, double lambda)`

- `n`: Number of draws
- `lambda`: Mean

See also: [rpois](https://web.mit.edu/r/current/lib/R/library/stats/html/Poisson.html)

## rt

`Vector rt(long n, double df, double ncp)`  
`Vector rt(long n, double df)`  
`Vector rt(Vector n, double df, double ncp)`  
`Vector rt(Vector n, double df)`  
`Vector rt(T)(T[] n, double df, double ncp)`  
`Vector rt(T)(T[] n, double df)`

- `n`: Number of draws. If `n` is a Vector or converts to a Vector, the number of draws is `n.rows`.
- `df`: Degrees of freedom
- `ncp`: Non-centrality parameter

See also: [rt](https://web.mit.edu/r/current/lib/R/library/stats/html/TDist.html)

## rweibull

`Vector rweibull(long n, double shape, double scale=1.0)`

- `n`: Number of draws
- `shape`
- `scale`

See also: [rweibull](https://web.mit.edu/r/current/lib/R/library/stats/html/Weibull.html)
