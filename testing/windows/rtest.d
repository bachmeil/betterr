import bindings;
import std;

void main() {
	bind_rinside("C:/Users/lanceb/packages/RInside/libs/x64", "libRInside", true);
	bind_libr("C:/Users/lanceb/R/R-4.4.1/bin/x64", "R", true);
	writeln(startR);
	writeln(setupRinC);
	
	startR();
	Rf_PrintValue(Rf_ScalarReal(-3.467));
	Rf_PrintValue(evalR("rnorm(12)"));
	evalRQ("library(urca)");
	Rf_PrintValue(evalR("ur.df"));
	Rf_PrintValue(RNil);
	closeR();
}
