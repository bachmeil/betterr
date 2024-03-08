/++
This module opens R's random number generation functions to D.
It exposes all functionality currently in base R.


There are two other options: D's RNG capabilities and the R API RNG
functions. This has the advantage of providing straightforward access to
parallel random number generators, since they're built into base R.
+/
module betterr.random;

import betterr.r;
import betterr.matrix, betterr.vector;
import std.conv, std.exception, std.stdio;

void setSeed(long seed) {
	evalRQ("set.seed(" ~ seed.to!string ~ ")");
}

struct Sample {
	Vector x;
	long size;
	bool replace = false;
	Vector prob;
	
	this(T)(T _x, long _size, bool _replace, T _prob) {
		x = _x.to!Vector;
		size = _size;
		replace = _replace;
		prob = _prob.to!Vector;
	}

	this(T)(T _x, long _size, bool _replace) {
		x = _x.to!Vector;
		size = _size;
		replace = _replace;
	}

	this(T)(T _x) {
		x = _x.to!Vector;
		size = x.rows;
	}

	Vector draw() {
		string cmd = "sample(";
		/* sample has a questionable implementation where a length of one
		 * for x causes x to be silently converted to 1:as.integer(x). That
		 * kind of thing is not consistent with D values. This eliminates
		 * that behavior. */
		if (x.rows == 1) {
			cmd ~= "c(" ~ x[0].to!string ~ ", " ~ x[0].to!string ~ ")";
		} else {
			cmd ~= x.name;
		}
		cmd ~= ", size=" ~ size.to!string;
		if (replace) {
			cmd ~= ", replace=TRUE";
		} else {
			enforce(x.rows >= size, "Cannot take " ~ size.to!string ~ " draws from " ~ x.rows.to!string ~ " elements of x unless you're sampling with replacement. Set replace to true if you want to do this.");
			cmd ~= ", replace=FALSE";
		}
		if (prob.rows == x.rows) {
			cmd ~= ", prob=" ~ prob.name;
		} else {
			enforce(prob.rows == 0, "prob inside Sample struct needs to be empty
or the same length as x");
		}
		cmd ~= ")";
		return Vector(cmd);
	}
}

/* dnorm, qnorm, and pnorm are already included in embedr.r
 * Those are direct calls to C code.
 * Handling rnorm is better this way, because (i) it returns a vector,
 * (ii) you can do things like set the generator when calling R.
 */
Vector rnorm(long n, double mean=0, double sd=1) {
	string cmd = "rnorm(" ~ n.to!string ~ ", mean=" ~ mean.to!string ~ ", sd=" ~ sd.to!string ~ ")";
	return Vector(cmd);
}

Vector runif(long n, double min=0, double max=1) {
	string cmd = "runif(" ~ n.to!string ~ ", min=" ~ min.to!string ~ ", max=" ~ max.to!string ~ ")";
	return Vector(cmd);
}

/* R has named optional parameters. D does not. Set rate to double.nan 
 * to leave it unspecified. */
Vector rgamma(long n, double shape, double rate=double.nan, double scale=1.0) {
	import std.math.traits: isNaN;
	string cmd;
	if (isNaN(rate)) {
		cmd = "rgamma(" ~ n.to!string ~ ", " ~ shape.to!string ~ ", scale=" ~ scale.to!string  ~ ")";
	} else {
		cmd = "rgamma(" ~ n.to!string ~ ", " ~ shape.to!string ~ ", rate=" ~ rate.to!string  ~ ")";
	}
	return Vector(cmd);
}

Vector rbeta(long n, double shape1, double shape2, double ncp=0.0) {
	string cmd = "rbeta(" ~ n.to!string ~ ", " ~ shape1.to!string ~ ", " ~ shape2.to!string ~ ", " ~ ncp.to!string ~ ")";
	return Vector(cmd);
}

IntVector rbinom(long n, long size, double prob) {
	string cmd = "rbinom(" ~ n.to!string ~ ", " ~ size.to!string ~ ", " ~ prob.to!string ~ ")";
	return IntVector(cmd);
}

Vector rcauchy(long n, double location=0, double scale=1) {
	string cmd = "rcauchy(" ~ n.to!string ~ ", " ~ location.to!string ~ ", " ~ scale.to!string ~ ")";
	return Vector(cmd);
}

/* df can be non-integer according to the R documentation */
Vector rchisq(long n, double df, double ncp) {
	enforce(df >= 0.0, "rchisq: df has to be non-negative");
	enforce(ncp >= 0.0, "rchisq: ncp has to be non-negative");
	string cmd = "rchisq(" ~ n.to!string ~ ", " ~ df.to!string ~ ", " ~ ncp.to!string ~ ")";
	return Vector(cmd);
}

Vector rexp(long n, double rate=1.0) {
	string cmd = "rexp(" ~ n.to!string ~ ", " ~ rate.to!string ~ ")";
	return Vector(cmd);
}

/* See the R documentation for details on ncp. That is why there are two versions */
Vector rf(long n, long df1, long df2, double ncp) {
	string cmd = "rf(" ~ n.to!string ~ ", " ~ df1.to!string ~ ", " ~ df2.to!string ~ ", " ~ ncp.to!string  ~ ")";
	return Vector(cmd);
}

Vector rf(long n, long df1, long df2) {
	string cmd = "rf(" ~ n.to!string ~ ", " ~ df1.to!string ~ ", " ~ df2.to!string ~ ")";
	return Vector(cmd);
}

/* prob needs to be a vector or anything that converts to a Vector */
IntVector rgeom(T)(long n, T prob) {
  auto p = Vector(prob);
  string cmd = "rgeom(" ~ n.to!string ~ ", " ~ p.name ~ ")";
  return IntVector(cmd);
}

IntVector rhyper(long nn, long m, long n, long k) {
  enforce(m+n >= k, "rhyper: Cannot draw more values (k) than the total number of values to choose from without replacement (m+n)");
  string cmd = "rhyper(" ~ nn.to!string ~ ", " ~ m.to!string ~ ", " ~ n.to!string ~ ", " ~ k.to!string ~ ")";
  return IntVector(cmd);
}

IntVector rhyper(Vector nn, long m, long n, long k) {
  return rhyper(nn.rows, m, n, k);
}

IntVector rhyper(long[] nn, long m, long n, long k) {
  return rhyper(nn.length, m, n, k);
}

Vector rlnorm(long n, double meanlog=0, double sdlog=1) {
  string cmd = "rlnorm(" ~ n.to!string ~ ", " ~ meanlog.to!string ~ ", " ~ sdlog.to!string ~ ")";
  return Vector(cmd);
}

IntMatrix rmultinom(T)(long n, long size, T prob) {
  auto p = Vector(prob);
  string cmd = "rmultinom(" ~ n.to!string ~ ", " ~ size.to!string ~ ", " ~ p.name ~ ")";
  return IntMatrix(cmd);
}

/* R's implementation allows for the specification of either prob or mu,
 * but not both. If you want mu, specify prob to be double.nan. 
 * Returns an integer matrix. 
 * 
 * For some reason this *might* return a double, but I have no idea why.
 * The help says that happens if the max integer is exceeded, but that's
 * definitely not the case. */
IntVector rnbinom(long n, double size, double prob, double mu) {
	import std.math.traits: isNaN;
	string cmd;
	if (isNaN(prob)) {
		cmd = "as.integer(rnbinom(" ~ n.to!string ~ ", " ~ size.to!string ~ ", mu=" ~ mu.to!string ~ "))";
	} else {
		cmd = "as.integer(rnbinom(" ~ n.to!string ~ ", " ~ size.to!string ~ ", prob=" ~ prob.to!string ~ "))";
	}
	return IntVector(cmd);
}

IntVector rnbinom(Vector n, double size, double prob, double mu) {
	return rnbinom(n.rows, size, prob, mu);
}

IntVector rnbinom(T)(T[] n, double size, double prob, double mu) {
	return rnbinom(n.length, size, prob, mu);
}

IntVector rnbinom(long n, double size, double prob) {
	return rnbinom(n, size, prob, 0.0);
}

IntVector rnbinom(Vector n, double size, double prob) {
	return rnbinom(n.rows, size, prob, 0.0);
}

IntVector rnbinom(T)(T[] n, double size, double prob) {
	return rnbinom(n.length, size, prob, 0.0);
}

/* R does not require either parameter to be an integer. */
IntVector rpois(double n, double lambda) {
  string cmd = "rpois(" ~ n.to!string ~ ", " ~ lambda.to!string ~ ")";
	return IntVector(cmd);
}

Vector rt(long n, double df, double ncp) {
	string cmd = "rt(" ~ n.to!string ~ ", " ~ df.to!string ~ ", " ~ ncp.to!string ~ ")";
	return Vector(cmd);
}

Vector rt(long n, double df) {
	string cmd = "rt(" ~ n.to!string ~ ", " ~ df.to!string ~ ")";
	return Vector(cmd);
}

Vector rt(Vector n, double df, double ncp) {
	return rt(n.rows, df, ncp);
}

Vector rt(Vector n, double df) {
	return rt(n.rows, df);
}

Vector rt(T)(T[] n, double df, double ncp) {
	return rt(n.length, df, ncp);
}

Vector rt(T)(T[] n, double df) {
	return rt(n.length, df);
}

Vector rweibull(long n, double shape, double scale=1.0) {
	string cmd = "rweibull(" ~ n.to!string ~ ", " ~ shape.to!string ~ ", " ~ scale.to!string ~ ")";
	return Vector(cmd);
}

/* Specialized generators for all the distributions above
 * Define a custom generator, for instance to use with a third-party
 * package that generates random numbers, by supplying cmd.
 * 
 * If you choose type custom, you'll have to add cmd manually. */
struct Generator(string type, long _length=1000, T=double) {
  long length = _length;
  /* To guarantee that refill gets called first */
  int currentIndex = _length; 
  
  static string drawfn(string type="double") {
    return type ~ ` draw() {
      if (currentIndex >= length) {
        refill();
        currentIndex = 0;
      }
      auto result = data[currentIndex];
      currentIndex += 1;
      return result;
    }`;
  }

  static if(type == "custom") {
    string cmd;
    Vector data;

    static if(is(T == double)) {
      void refill() {
        data = Vector(cmd);
        length = data.length;
      }
      mixin(drawfn("double"));
    }
    static if(is(T == int)) {
      void refill() {
        data = IntVector(cmd);
        length = data.length;
      }
      mixin(drawfn("int"));
    }
    static if(is(T == bool)) {
      void refill() {
        data = BoolVector(cmd);
        length = data.length;
      }
      mixin(drawfn("bool"));
    }
  }
  static if(type == "norm") {
    Vector data;
    double mean = 0.0;
    double sd = 1.0;
    
    void refill() {
      data = rnorm(length, mean, sd);
    }
    
    mixin(drawfn());
  }
  static if(type == "unif") {
    Vector data;
    double min = 0.0;
    double max = 1.0;
    
    void refill() {
      data = runif(length, min, max);
    }
    
    mixin(drawfn());
  }
  static if(type == "gamma") {
    Vector data;
    double shape;
    double rate = double.nan;
    double scale = 1.0;
    
    void refill() {
      data = rgamma(length, shape, rate, scale);
    }
    
    mixin(drawfn());
  }
  static if(type == "beta") {
    Vector data;
    double shape1;
    double shape2;
    double ncp = 0.0;
    
    void refill() {
      data = rbeta(length, shape1, shape2, ncp);
    }
    
    mixin(drawfn());
  }
  static if(type == "binom") {
    IntVector data;
    long size;
    double prob;
    
    void refill() {
      data = rbinom(length, size, prob);
    }
    
    mixin(drawfn("int"));
  }
}

void prngInit(long stream, long seed=1) {
  evalRQ("library(parallel)");
  evalRQ(`RNGkind("L'Ecuyer-CMRG")`);
  setSeed(seed);
  foreach(_; 0..stream) {
    evalRQ(".Random.seed <- nextRNGStream(.Random.seed)");
  }
}

/* GSL versions
 * Take a pointer and a length
 * That's a general version because these functions can be called many
 * times and I don't want to make decisions for the caller about data
 * allocation. There's an expectation of higher-level wrappers being
 * called. */
version(gsl) {
  import gslheaders;
  
  double rnorm(gsl_rng * r) {
    return gsl_ran_ugaussian(r);
  }
  
  void rnorm(gsl_rng * r, double * v, long n) {
    foreach(ii; 0..n) {
      v[ii] = gsl_ran_ugaussian(r);
    }
  }
  
  void rnorm(gsl_rng * r, double[] v) {
    rnorm(r, v.ptr, v.length);
  }
  
  double rnorm(gsl_rng * r, double sigma) {
    return gsl_ran_gaussian(r, sigma);
  }
  
  double rnorm(gsl_rng * r, double mu, double sigma) {
    return mu + gsl_ran_gaussian(r, sigma);
  }
  
  void rnorm(gsl_rng * r, double * v, long n, double sigma) {
    foreach(ii; 0..n) {
      v[ii] = rnorm(r, sigma);
    }
  }
  
  void rnorm(gsl_rng * r, double[] v, double sigma) {
    rnorm(r, v.ptr, v.length, sigma);
  }
  
  void rnorm(gsl_rng * r, double * v, long n, double mu, double sigma) {
    foreach(ii; 0..n) {
      v[ii] = rnorm(r, mu, sigma);
    }
  }
  
  void rnorm(gsl_rng * r, double[] v, double mu, double sigma) {
    rnorm(r, v.ptr, v.length, mu, sigma);
  }
  
  double[] rnorm(gsl_rng * r, long k) {
    auto result = new double[k];
    rnorm(r, result);
    return result;
  }
  
  double[] rnorm(gsl_rng * r, long k, double sigma) {
    auto result = new double[k];
    rnorm(r, result, sigma);
    return result;
  }
  
  double[] rnorm(gsl_rng * r, long k, double mu, double sigma) {
    auto result = new double[k];
    rnorm(r, result, mu, sigma);
    return result;
  }
    
    
  //~ double runif(gsl_rng *r) {
  //~ }
  
  //~ double rgamma(gsl_rng *r) {
  //~ }
}
  
  
  
  
