module cran.vars;
import betterr.r, betterr.list;
import std.conv;
import std.stdio, std.sumtype;

struct VAR {
  long p = 1;
  string type = "const";
  long lagMax = -1;
  string ic = "AIC";
  
  List fit(T)(T x) {
    string cmd = `vars::VAR(` ~ x.name ~ `, p=` ~ p.to!string ~ `, type="` ~ type ~ `"`;
    if (lagMax > 0) {
      cmd ~= `, lag.max=` ~ lagMax.to!string;
    } 
    cmd ~= `, ic="` ~ ic ~ `")`;
    writeln(cmd);
    return VAREst(RData(cmd));
  }
}

struct VAREst {
	RData output;
	
	List varresult() {
		return List(output.name ~ "$varresult");
	}
	
	Matrix datamat() {
		return Matrix(output.name ~ "$datamat");
	}
	
	Matrix y() {
		return Matrix(output.name ~ "$y");
	}
	
	string type() {
		return evalR(output.name ~ "$type").scalar!string;
	}
	
	int p() {
		return Rf_asInteger(evalR(output.name ~ "$p"));
	}
	
	int K() {
		return Rf_asInteger(evalR(output.name ~ "$K"));
	}
	
	int obs() {
		return Rf_asInteger(evalR(output.name ~ "$obs"));
	}
	
	int totobs() {
		return Rf_asInteger(evalR(output.name ~ "$totobs"));
	}
	
	MaybeMatrix restrictions() {
		auto tmp = Rf_protect(evalR(output.name ~ "$restrictions"));
		if (Rf_isNull(tmp)) {
			Rf_unprotect(1);
			return null;
		} else {
			Rf_unprotect(1);
			return Matrix(output.name ~ "$restrictions");
		}
	}	
}

alias MaybeMatrix = SumType!(null, Matrix);

struct VAREquations {
	LMFit[] equations;
	alias this equations;
	
	this(List varfit) {
		auto tmp = Rf_protect(varfit.name ~ "$varresult");
		foreach(ii; 0..tmp.length) {
			equations ~= LMFit(tmp[ii]);
		}
	}
}
	
	
	
	
	
	
