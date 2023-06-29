import betterr.r, betterr.vector;
import std.stdio;

extern(C) {
	Robj callfunction(Robj ptr, Robj par);
}

Robj fun1(Robj v) {
	double result = 0.0;
	auto vec = REAL(v)[0..v.length];
	foreach(val; vec) {
		result += val;
	}
	return result.robj;
}

Robj fun2(Robj v) {
	auto vec = REAL(v)[0..v.length];
	return robj(100.0*(vec[1]-vec[0]*vec[0])^^2 + (1-vec[0])^^2);
} 

Robj fun3(Robj v) {
	auto result = Rf_protect(Rf_allocVector(14, 2));
	double * res = REAL(result);
	auto vec = REAL(v)[0..v.length];
	res[0] = -400.0*vec[0]*(vec[1]-vec[0]*vec[0]) - 2*(1-vec[0]);
	res[1] = 200.0*(vec[1]-vec[0]*vec[0]);
	Rf_unprotect(1);
	return result;
}

struct Data {
	double[] x;
	
	Robj sum() {
		double result = 0.0;
		foreach(val; x) {
			result += val;
		}
		return robj(result);
	}
}

/* Be sure to not omit the extern(C) or it will reverse the order of
 * the arguments! */
extern(C) Robj fun4(Robj par, Robj ptr) {
	void * tmp = R_ExternalPtrAddr(ptr);
	Data * tmp2 = cast(Data *) tmp;
	return tmp2.sum();
}

Robj toExtPtr(alias f)() {
	extern(C) void function() f2 = cast(void function()) &f;
	return R_MakeExternalPtrFn(f2, RNil, Rf_protect(RNil));
}

void main() {
	startR();
	extern(C) void function() func2 = cast(void function()) &fun1;
	Robj p = R_MakeExternalPtrFn(func2, RNil, Rf_protect(RNil));
	toR(p, "thefunction");
	evalRQ("library(funcptr)");
	writeln("Got here");
	evalRQ("print(.Call('callfunction', thefunction, c(1.4, 3.2, 1.7)))");
	toR(toExtPtr!(fun2)(), "theobjfn");
	toR(toExtPtr!(fun3)(), "thegradfn");
	evalRQ("f1 <- function(par) { .Call('callfunction', theobjfn, par) }");
	evalRQ("g1 <- function(par) { .Call('callfunction', thegradfn, par) }");
	writeln("here");
	printR("f1(c(1,2))");
	printR("g1");
	auto tmp = Vector([1.0, 2.0]);
	printR(fun3(tmp.data.x));
	printR("g1(c(1,2))");
	printR("theobjfn");
	printR(evalR("constrOptim(c(-1.2,0.9), f1, g1, ui = rbind(c(-1,0), c(0,-1)), ci = c(-1,-1))"));
	writeln("here 2");
	
	Data dd;
	dd.x = [1.1, 2.2, 3.3];
	writeln("pointer to dd is ", &dd);
	writeln("void pointer to dd is ", cast(void*) &dd);
	
	toR(toExtPtr!fun4, "function4");
	printR(R_MakeExternalPtr(cast(void *) &dd, RNil, Rf_protect(RNil)));
	writeln(R_ExternalPtrAddr(R_MakeExternalPtr(cast(void *) &dd, RNil, Rf_protect(RNil))));
	toR(R_MakeExternalPtr(cast(void *) &dd, RNil, Rf_protect(RNil)), "dataptr");
	evalRQ("print(.Call('callfunction2', function4, NULL, dataptr))");
	writeln(dd);
	closeR();
}
