module betterr.baser;

import betterr.list, betterr.rdata, betterr.vector;
import betterr.r;
import std.array, std.conv, std.string;

/* Note: This file is harder to read because of the use of mixins.
 * This was the best way to handle it - better than duplicating functions
 * and better than templates. */

/* Function takes one argument, but requires different treatment
 * for Vector, double, int, and long arguments. Returns the same
 * type as the argument. */
string function0(string name, string rname) {
  return `double NAME(double x) {
  return evalR("RNAME(" ~ x.to!string ~ ")").scalar;
}

int NAME(int x) {
  return evalR("RNAME(" ~ x.to!string ~ "L)").scalar!int;
}

long NAME(long x) {
  return evalR("RNAME(" ~ x.to!string ~ "L)").scalar!long;
}

Vector NAME(Vector v) {
  return Vector("RNAME(" ~ v.name ~ ")");
}

Vector NAME(double[] arr) {
  auto v = Vector(arr);
  return Vector("RNAME(" ~ v.name ~ ")");
}
`.replace("RNAME", rname).replace("NAME", name);
}

string function0(string name) {
  return function0(name, name);
}

string function1(string name, string rname) {
  return `double ` ~ name ~ `(RData obj, bool removeNA=false) {
  if (removeNA) {
    return evalR("` ~ rname ~ `(" ~ obj.name ~ ", na.rm=TRUE)").scalar;
  } else {
    return evalR("` ~ rname ~ `(" ~ obj.name ~ ")").scalar;
  }
}
`;
}

string function1(string name) {
  return function1(name, name);
}

/* This function type takes double, long, and Vector arguments,
 * returning a double in the first two cases, and a Vector in the last.
 * There are no optional arguments. */
string function2(string name, string rname) {
  return `
double NAME(double x) {
  return evalR("RNAME(" ~ x.to!string ~ ")").scalar!double;
}

double NAME(long x) {
  return evalR("RNAME(" ~ x.to!string ~ ")").scalar!double;
}

Vector NAME(Vector v) {
  return Vector("RNAME(" ~ v.name ~ ")");
}

Vector NAME(double[] arr) {
  auto v = Vector(arr);
  return Vector("RNAME(" ~ v.name ~ ")");
}
`.replace("RNAME", rname).replace("NAME", name);
}

string function2(string name) {
  return function2(name, name);
}

/* This function type takes a Vector as an argument. It returns a
 * Vector. There are no optional arguments. */
string function3(string name, string rname) {
  return `
Vector NAME(Vector v) {
  return Vector("RNAME(" ~ v.name ~ ")");
}

Vector NAME(double[] arr) {
  auto v = Vector(arr);
  return Vector("RNAME(" ~ v.name ~ ")");
}
`.replace("RNAME", rname).replace("NAME", name);
}

string function3(string name) {
  return function3(name, name);
}

/* This function type takes double, long, and Vector arguments,
 * returning a double in the first two cases, and a Vector in the last.
 * There are optional arguments with a default value. */
string function4(string name, string rname, string _default) {
  return `
double NAME(double x, string opt="OPTIONS") {
  return evalR("RNAME(" ~ x.to!string ~ ", " ~ opt ~ ")").scalar!double;
}

double NAME(long x, string opt="OPTIONS") {
  return evalR("RNAME(" ~ x.to!string ~ ", " ~ opt ~ ")").scalar!double;
}

Vector NAME(Vector v, string opt="OPTIONS") {
  return Vector("RNAME(" ~ v.name ~ ", " ~ opt ~ ")");
}

Vector NAME(double[] arr, string opt="OPTIONS") {
  auto v = Vector(arr);
  return Vector("RNAME(" ~ v.name ~ ", " ~ opt ~ ")");
}
`.replace("RNAME", rname).replace("NAME", name).replace("OPTIONS", _default);
}

/* This function type takes double, long, and Vector arguments,
 * returning a double in all cases. There are optional arguments with 
 * a default value. 
 * 
 * Obviously a lot of code duplication, but it's possible I'll change
 * this in the future, and this makes it easier to work with. It's not
 * a problem anyway since I'm using this to generate code. */
string function5(string name, string rname, string _default) {
  return `
double NAME(double x, string opt="OPTIONS") {
  return evalR("RNAME(" ~ x.to!string ~ ", " ~ opt ~ ")").scalar!double;
}

double NAME(long x, string opt="OPTIONS") {
  return evalR("RNAME(" ~ x.to!string ~ ", " ~ opt ~ ")").scalar!double;
}

double NAME(Vector v, string opt="OPTIONS") {
  return evalR("RNAME(" ~ v.name ~ ", " ~ opt ~ ")").scalar!double;
}

double NAME(double[] arr, string opt="OPTIONS") {
  auto v = Vector(arr);
  return evalR("RNAME(" ~ v.name ~ ", " ~ opt ~ ")").scalar!double;
}
`.replace("RNAME", rname).replace("NAME", name).replace("OPTIONS", _default);
}

/* This function type takes a Vector as an argument. It returns a
 * double. There are optional arguments. */
string function6(string name, string rname, string _default) {
  return `
double NAME(Vector v, string opt="OPTIONS") {
  return evalR("RNAME(" ~ v.name ~ ", " ~ opt ~ ")").scalar!double;
}

double NAME(double[] arr, string opt="OPTIONS") {
  auto v = Vector(arr);
  return evalR("RNAME(" ~ v.name ~ ", " ~ opt ~ ")").scalar!double;
}
`.replace("RNAME", rname).replace("NAME", name).replace("OPTIONS", _default);
}

/* This function type takes a Vector as an argument. It returns a
 * Vector. There are optional arguments. */
string function7(string name, string rname, string _default) {
  return `
Vector NAME(Vector v, string opt="OPTIONS") {
  return Vector("RNAME(" ~ v.name ~ ", " ~ opt ~ ")");
}

Vector NAME(double[] arr, string opt="OPTIONS") {
  auto v = Vector(arr);
  return Vector("RNAME(" ~ v.name ~ ", " ~ opt ~ ")");
}
`.replace("RNAME", rname).replace("NAME", name).replace("OPTIONS", _default);
}

mixin(function0("abs"));
mixin(function0("sqrt"));
mixin(function0("ceiling"));
mixin(function0("floor"));
mixin(function0("trunc"));
mixin(function4("round", "round", "digits = 0"));
mixin(function4("signif", "signif", "digits = 6"));

/* removeNA: set to true to remove NA values before computing */
mixin(function1("sum"));

/* This is what we're generating with these mixins */
//~ double mean(RData obj, bool removeNA=false) {
  //~ if (removeNA) {
    //~ return evalR("mean(" ~ obj.name ~ ", na.rm=TRUE)").scalar;
  //~ } else {
    //~ return evalR("mean(" ~ obj.name ~ ")").scalar;
  //~ }
//~ }

mixin(function2("cosh"));
mixin(function2("sinh"));
mixin(function2("tanh"));
mixin(function2("acosh"));
mixin(function2("asinh"));
mixin(function2("atanh"));

mixin(function3("cumsum"));
mixin(function3("cumprod"));
mixin(function3("cummax"));
mixin(function3("cummin"));

mixin(function4("log", "log", "base = exp(1)"));
mixin(function4("logb", "log", "base = exp(1)"));
mixin(function2("log10"));
mixin(function2("log2"));
mixin(function2("log1p"));
mixin(function2("exp"));
mixin(function2("expm1"));

mixin(function6("prod", "prod", "na.rm = FALSE"));
mixin(function7("range", "range", "na.rm = FALSE, finite = FALSE"));
mixin(function7("rank", "rank", "na.last = TRUE, ties.method = 'average'"));

mixin(function3("rev"));

Vector seq(double from, double including, double by) {
  return Vector("seq(" ~ from.to!string ~ ", " ~ including.to!string ~ ", " ~ by.to!string ~ ")");
}

mixin(function7("sort", "sort", "decreasing = FALSE, na.last = NA"));

List summary(Vector x) {
  return List("as.list(summary(" ~ x.name ~ "))");
}

List summary(double[] x) {
  return summary(Vector(x));
}

/* For a single quantile of a series, use this */
//~ double quantile(T)(T x, double q) {
  //~ return evalR("quantile(" ~ x.name ~ ", " ~ q.to!string ~ ")").scalar;
//~ }

/* For the full functionality of the quantiles function, use this struct. */
struct Quantiles {
  Vector probs;
  bool narm = false;
  long type = 7;
  
  this(Vector _probs) {
    probs = _probs;
  }
  
  this(T)(T _probs) {
    probs = Vector(_probs);
  }
  
  Vector values(T)(T x) {
    string probsPart;
    string narmPart;
    if (narm) {
      narmPart = ", na.rm=TRUE";
    } else {
      narmPart = ", na.rm=FALSE";
    }
    return Vector("quantile(" ~ x.name ~ ", probs=" ~ probs.name ~ narmPart ~ ", type=" ~ type.to!string ~ ")");
  }
  
  double value(T)(T x) {
    string narmPart;
    if (narm) {
      narmPart = ", na.rm=TRUE";
    } else {
      narmPart = ", na.rm=FALSE";
    }
    return Vector("quantile(" ~ x.name ~ ", probs=" ~ probs[0].to!string ~ narmPart ~ ", type=" ~ type.to!string ~ ")");
  }
}

/* These are fundamental building blocks, so they should be fast.
 * Speed isn't be an issue for a large vector inside R, but we need to be
 * performant for double[] as well. */
double mean(double x) {
  return x;
}

version(mir) {
	/* If you don't want to remove NaN values, so if there are any, it
	 * returns double.nan. */
	double mean(double[] x) {
		import mir.stat.descriptive.univariate;
		return mir.stat.descriptive.univariate.mean(x);
	}

	/* If you want the mean of a Vector or Matrix, can avoid calling into R
	 * Returns double.nan if any missing values. */
	double mean(T)(T v) {
		return mean(v.ptr[0..v.length]);
	}

	/* Maybe a better way to do this, but it allows me to pass the filtered
	 * values to this function without first converting to double[]. */
	double mean(string filtered, T)(T v) {
		import mir.stat.descriptive.univariate;
		return mir.stat.descriptive.univariate.mean(v);
	}

	/* This is hopefully faster than shipping it to R. */
	double mean(double[] x, bool narm=false) {
		if (!narm) {
			return mean(x);
		}
		import std.math.traits: isNaN;
		import std.algorithm.iteration: filter;
		return mean!""(x.filter!(a => !isNaN(a)));
	}

	double mean(T)(T x, bool narm=false) {
		return mean(x.ptr[0..x.length], narm);
	}

	/* If you want to trim, it's going to R. It would be a lot of work for
	 * a minor use case. The cost of copying at that point is small
	 * relative to all the other stuff that needs to be done for a trimmed
	 * mean calculation. */
	double mean(double[] x, double trim=0, bool narm=false) {
		return mean(Vector(x), trim, narm);
	}

	double mean(T)(T v, double trim=0, bool narm=false) {
		return evalR("mean(" ~ v.name ~ ", trim=" ~ trim.to!string ~ ", na.rm=" ~ boolString(narm) ~ ")").scalar;
	}

	/* Now do median */
	double median(double x) {
		return x;
	}

	/* If you don't want to remove NaN values, so if there are any, it
	 * returns double.nan. */
	double median(double[] x) {
		import mir.stat.descriptive.univariate;
		return mir.stat.descriptive.univariate.median(x);
	}

	/* If you want the mean of a Vector or Matrix, can avoid calling into R
	 * Returns double.nan if any missing values. */
	double median(T)(T v) {
		return median(v.ptr[0..v.length]);
	}

	/* Maybe a better way to do this, but it allows me to pass the filtered
	 * values to this function without first converting to double[]. */
	double median(string filtered, T)(T v) {
		import mir.stat.descriptive.univariate;
		return mir.stat.descriptive.univariate.median(v);
	}

	/* This is hopefully faster than shipping it to R. */
	double median(double[] x, bool narm=false) {
		if (!narm) {
			return mean(x);
		}
		import std.math.traits: isNaN;
		import std.algorithm.iteration: filter;
		return median!""(x.filter!(a => !isNaN(a)).array);
	}

	double median(T)(T x, bool narm=false) {
		return median(x.ptr[0..x.length], narm);
	}

	/* And sum */
	double sum(double x) {
		return x;
	}

	/* If you don't want to remove NaN values, so if there are any, it
	 * returns double.nan. */
	double sum(double[] x) {
		import mir.math.sum;
		return mir.math.sum.sum!"fast"(x);
	}

	/* If you want the sum of a Vector or Matrix, can avoid calling into R
	 * Returns double.nan if any missing values. */
	double sum(T)(T v) {
		return sum(v.ptr[0..v.length]);
	}

	/* Maybe a better way to do this, but it allows me to pass the filtered
	 * values to this function without first converting to double[]. */
	double sum(string filtered, T)(T v) {
		import mir.math.sum;
		return mir.math.sum.sum!"fast"(v);
	}

	/* This is hopefully faster than shipping it to R. */
	double sum(double[] x, bool narm=false) {
		if (!narm) {
			return sum(x);
		}
		import std.math.traits: isNaN;
		import std.algorithm.iteration: filter;
		return sum!""(x.filter!(a => !isNaN(a)));
	}

	double sum(T)(T x, bool narm=false) {
		return sum(x.ptr[0..x.length], narm);
	}


	/* Returns NaN if any are NaN */
	double max(double[] x) {
		import std.algorithm.iteration: fold;
		import std.math.traits: isNaN;
		return x.fold!((a,b) => a.isNaN || b.isNaN? real.nan: a < b? b: a);
	}

	double max(T)(T x) {
		return max(x.ptr[0..x.length]);
	}

	/* Set narm to true to remove NaN values first */
	double max(double[] x, bool narm=false) {
		if (!narm) {
			return max(x);
		}
		import std.algorithm.searching: maxElement;
		import std.algorithm.iteration: filter;
		import std.math.traits: isNaN;
		auto tmp = x.filter!(a => !isNaN(a));
		if (tmp.empty) {
			return double.nan;
		} else {
			return tmp.maxElement;
		}
	}

	double max(T)(T x, bool narm=false) {
		if (!narm) {
			return max(x);
		} else {
			return max(x.ptr[0..x.length]);
		}
	}

	/* Returns NaN if any are NaN */
	double min(double[] x) {
		import std.algorithm.iteration: fold;
		import std.math.traits: isNaN;
		return x.fold!((a,b) => a.isNaN || b.isNaN? real.nan: a < b? a: b);
	}

	double min(T)(T x) {
		return min(x.ptr[0..x.length]);
	}

	/* Set narm to true to remove NaN values first */
	double min(double[] x, bool narm=false) {
		if (!narm) {
			return max(x);
		}
		import std.algorithm.searching: minElement;
		import std.algorithm.iteration: filter;
		import std.math.traits: isNaN;
		auto tmp = x.filter!(a => !isNaN(a));
		if (tmp.empty) {
			return double.nan;
		} else {
			return tmp.minElement;
		}
	}

	double min(T)(T x, bool narm=false) {
		if (!narm) {
			return min(x);
		} else {
			return min(x.ptr[0..x.length], true);
		}
	}

	double sd(double[] x) {
		import mir.math.stat;
		return mir.math.stat.standardDeviation(x);
	}

	double sd(T)(T x) {
		return sd(x.ptr[0..x.length]);
	}

	double sd(double[] x, bool narm=false) {
		if (!narm) {
			return sd(x);
		} else {
			import std.algorithm.iteration: filter;
			import std.math.traits: isNaN;
			import mir.math.stat;
			return mir.math.stat.standardDeviation(x.filter!(a => !isNaN(a)));
		}
	}

	double sd(T)(T x, bool narm=false) {
		return sd(x.ptr[0..x.length], narm);
	}

	double var(double[] x) {
		import mir.math.stat;
		return mir.math.stat.variance(x);
	}

	double var(T)(T x) {
		return var(x.ptr[0..x.length]);
	}

	double var(double[] x, bool narm=false) {
		if (!narm) {
			return var(x);
		} else {
			import std.algorithm.iteration: filter;
			import std.math.traits: isNaN;
			import mir.math.stat;
			return mir.math.stat.variance(x.filter!(a => !isNaN(a)));
		}
	}

	double var(T)(T x, bool narm=false) {
		return var(x.ptr[0..x.length], narm);
	}

	double quantile(double[] x, double p) {
		import mir.stat.descriptive.univariate;
		return mir.stat.descriptive.univariate.quantile(x, p);
	}

	/* By default, double.nan is counted as the median of the others */
	double quantile(double[] x, double p, bool narm) {
		if (!narm) {
			return quantile(x, p);
		} else {
			import std.algorithm.iteration: filter;
			import std.math.traits: isNaN;
			return quantile(x.filter!(a => !isNaN(a)).array, p);
		}
	}

	double quantile(T)(T x, double p) {
		return quantile(x.ptr[0..x.length], p);
	}

	double quantile(T)(T x, double p, bool narm) {
		return quantile(x.ptr[0..x.length], p, narm);
	}

	double[] quantile(double[] x, double[] p) {
		import mir.stat.descriptive.univariate;
		return mir.stat.descriptive.univariate.quantile!("type7", false, true)(x, p.dup).array;
	}

	double[] quantile(T)(T x, double[] p) {
		return quantile(x.ptr[0..x.length], p.dup);
	}

	double[] quantile(double[] x, double[] p, bool narm) {
		if (!narm) {
			return quantile(x, p);
		} else {
			import mir.stat.descriptive.univariate;
			import std.algorithm.iteration: filter;
			import std.math.traits: isNaN;
			return mir.stat.descriptive.univariate.quantile!("type7", false, true)(x.filter!(a => !isNaN(a)).array, p.dup).array;
		}
	}

	double[] quantile(T)(T x, double[] p, bool narm) {
		if (!narm) {
			return quantile(x, p);
		} else {
			return quantile(x.ptr[0..x.length], p.dup, narm);
		}
	}
}

string boolString(bool x) {
	if (x) {
		return "TRUE";
	} else {
		return "FALSE";
	}
}
