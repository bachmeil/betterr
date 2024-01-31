module betterr.lm;

import std.array, std.conv, std.exception, std.stdio;
import betterr.r;
import betterr.dataframe, betterr.list, betterr.matrix, betterr.rdata, betterr.vector;

/* This is a configuration. Additional options will be added in the
 * future. The data is not stored in here. 
 * We need to handle the many ways data can be sent to lm. 
 * 
 * We use the Subset struct mainly so it can take a single integer or
 * an array of integers as the argument. */
struct LMConfig {
  string lhs;
  string[] rhs;
  bool intercept = true;
  Subset subset;
  
  struct Subset {
		long[] obs;
		alias obs this;
		
		this(long ii) {
			obs = [ii];
		}
		
		this(long[] ii) {
			obs = ii;
		}
		
		/* Following standard D practice, this is exclusive of i2. */
		this(long i1, long i2) {
			foreach(ii; i1..i2) {
				obs ~= ii;
			}
		}
		
		void opAssign(long a) {
			obs = [a];
		}
		
		void opOpAssign(string op)(long ii) {
			if (op == "~") {
				obs ~= ii;
			}
		}

		void opOpAssign(string op)(long[] ii) {
			if (op == "~") {
				obs ~= ii;
			}
		}

		void opOpAssign(string op)(long i1, long i2) {
			if (op == "~") {
				foreach(ii; i1..i2) {
					obs ~= ii;
				}
			}
		}
		
		IntVector vec() {
			return IntVector(obs);
		}			
	}
	
	string rhsFormula() {
		if (rhs.length > 0) {
			import std.string: join;
			return rhs.join("+");
		} else {
			return "";
		}
	}
}

LMFit lm(T1, T2)(T1 y, T2 x) {
	return lm(y.name, x.name, "", "", true);
}	

LMFit lm(T1, T2)(T1 y, T2 x, LMConfig conf) {
	if ( (conf.subset !is null) & (conf.subset.length > 0) ) {
		auto s = conf.subset.vec;
		LMFit result = lm(y.name, x.name, s.name, "", conf.intercept);
		s.free();
		return result;
	} else {
		return lm(y.name, x.name, "", "", conf.intercept);
	}
}

LMFit lm(DataFrame df, LMConfig conf) {
	if (conf.subset.length > 0) {
		return lm(conf.lhs, conf.rhsFormula, conf.subset.vec.name, df.name, conf.intercept);
	} else {
		return lm(conf.lhs, conf.rhsFormula, "", df.name, conf.intercept);
	}
}

/* These regressions should work down to this. Arguments are variable 
 * names as they appear inside R. */
LMFit lm(string y, string x, string subset, string data, bool intercept=true) {
	string optional;
	if (subset.length > 0) {
		optional ~= ", subset=" ~ subset;
	}
	if (data.length > 0) {
		optional ~= ", data=" ~ data;
	}
	string rhs = x;
	if (!intercept) {
		rhs ~= "-1";
	}
	string cmd = `lm(` ~ y ~ ` ~ ` ~ rhs ~ optional ~ `)`;
	
	LMFit result;
	result.fit = List(cmd);
	result.summary = List(`summary(` ~ result.fit.name ~ `)`);
	result.intercept = intercept;
	return result;
}		

/* Restructure a bit. Pass LMConfig struct as an optional argument if
 * you want to set arguments. The current approach is clumsy. */

struct LMFit {
	List fit;
  List summary;
  bool intercept;
  
  double pred(double x) {
		Vector b = beta();
		if (intercept) {
			enforce(b.length == 2, "Wrong number of predictors in call to pred");
			return b[0] + b[1]*x;
		} else {
			enforce(b.length == 1, "Wrong number of predictors in call to pred");
			return b[0]*x;
		}
	}
	
	double pred(Vector x) {
		Vector b = beta();
		double result = 0.0;
		if (intercept) {
			result += b[0];
			foreach(ii; 1..b.length) {
				result += b[ii]*x[ii-1];
			}
			return result;
		} else {
			foreach(ii; 0..b.length) {
				result += b[ii]*x[ii];
			}
			return result;
		}
	}
	
	void print(string msg="") {
		if (msg.length > 0) {
			writeln(msg ~ ":");
		}
		printR(fit.x);
	}
  
  Vector beta() {
		return Vector(fit["coefficients"]);
	}
	
	List coefficients() {
		return List(summary.name ~ "['coefficients']");
	}
	
	Vector residuals() {
		return Vector(fit["residuals"]);
	}
	
	Vector fittedValues() {
		return Vector(fit["fitted.values"]);
	}
	
	int dfResidual() {
		return fit["df.residual"].as!int;
	}
	
	List model() {
		return List(fit.name ~ "['model']");
	}
	
	double sigma() {
		return summary["sigma"].as!double;
	}
	
	double rsq() {
		return summary["r.squared"].as!double;
	}
	
	double adjrsq() {
		return summary["adj.r.squared"].as!double;
	}
	
	double fstat() {
		return summary["fstatistic"].as!double;
	}
	
	List unscaledCov() {
		return List(summary.name ~ "['cov.unscaled']");
	}
}

extern(C) {
	void dqrls_(double *x, int *n, int *p, double *y, int *ny,
		     double *tol, double *b, double *rsd, double *qty, int *k, 
		     int *jpvt, double *qraux, double *work);
}

struct RegOutput {
	double[] coef;
	double[] residuals;
	double[] effects;
	int rank;
}

RegOutput dqrls(Vector y, Matrix x, double tol=1e-07) {
	RegOutput result;
	int n = x.rows.to!int;
	int p = x.cols.to!int;
	int ny = 1;
	result.coef = new double[p];
	result.residuals = new double[n];
	result.effects = new double[n];
	double[] xreg;
	xreg.reserve(n*p);
	foreach(val; x.ptr[0..n*p]) {
		xreg ~= val;
	}
	//~ foreach(val; y.ptr[0..n]) {
		//~ yreg ~= val;
		//~ result.residuals ~= val;
		//~ result.effects ~= val;
	//~ }
	auto pivot = new int[p];
	auto qraux = new double[p];
	auto work = new double[2*p];
	dqrls_(xreg.ptr, &n, &p, y.ptr, &ny, &tol, result.coef.ptr,
		result.residuals.ptr, result.effects.ptr, &(result.rank),
		pivot.ptr, qraux.ptr, work.ptr);
	return result;
}
