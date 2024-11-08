struct Sexprec {}
alias Robj = Sexprec*;
import std.utf;

extern(C) {
	double* function(Sexprec*) REAL;
	int* function(Sexprec*) INTEGER;
	void function(Sexprec*) Rf_PrintValue;
	Sexprec* function(double) Rf_ScalarReal;
	const(char)* function(Sexprec*) R_CHAR;
	Sexprec* function(Sexprec*, int) VECTOR_ELT;
	Sexprec* function(uint, int) Rf_allocVector;
	Sexprec* function(uint, int, int) Rf_allocMatrix;
	Sexprec* function(int) Rf_allocList;

	void function(Sexprec*, char*) passToR;
	Sexprec* function(char*) evalInR;
	void function(char*) evalQuietlyInR;
	void function() setupRinC;
	void function() teardownRinC;
}

alias startR = setupRinC;
alias closeR = teardownRinC;

void bind_libr(string dllPath, string dllName, bool errors=false) {
	import std.process, std.stdio, std.utf;
	import core.sys.windows.windows;
	import std.windows.syserror;
	import std.file, std.stdio;

	// Error 126 means *something* can't be found, including a dependency.
	HMODULE libr = LoadLibraryW(toUTF16z(dllName ~ ".dll"));
	if (errors) {
		writeln(sysErrorString(GetLastError()));
	}
	REAL = cast(typeof(REAL)) GetProcAddress(libr, "REAL");
	INTEGER = cast(typeof(INTEGER)) GetProcAddress(libr, "INTEGER");
	Rf_PrintValue = cast(typeof(Rf_PrintValue)) GetProcAddress(libr, "Rf_PrintValue");
	Rf_ScalarReal = cast(typeof(Rf_ScalarReal)) GetProcAddress(libr, "Rf_ScalarReal");
	R_CHAR = cast(typeof(R_CHAR)) GetProcAddress(libr, "R_CHAR");
	VECTOR_ELT = cast(typeof(VECTOR_ELT)) GetProcAddress(libr, "VECTOR_ELT");
	Rf_allocVector = cast(typeof(Rf_allocVector)) GetProcAddress(libr, "Rf_allocVector");
	Rf_allocMatrix = cast(typeof(Rf_allocMatrix)) GetProcAddress(libr, "Rf_allocMatrix");
	Rf_allocList = cast(typeof(Rf_allocList)) GetProcAddress(libr, "Rf_allocList");
}

void bind_rinside(string dllPath, string dllName, bool errors=false) {
	import std.process, std.stdio, std.utf;
	import core.sys.windows.windows;
	import std.windows.syserror;

	// Error 126 means *something* can't be found, including a dependency.
	HMODULE rinside = LoadLibraryW(toUTF16z(dllName ~ ".dll"));
	if (errors) {
		writeln(sysErrorString(GetLastError()));
	}
	passToR = cast(typeof(passToR)) GetProcAddress(rinside, "passToR");
	evalInR = cast(typeof(evalInR)) GetProcAddress(rinside, "evalInR");
	evalQuietlyInR = cast(typeof(evalQuietlyInR)) GetProcAddress(rinside, "evalQuietlyInR");
	setupRinC = cast(typeof(setupRinC)) GetProcAddress(rinside, "setupRinC");
	teardownRinC = cast(typeof(teardownRinC)) GetProcAddress(rinside, "teardownRinC");
}

Robj evalR(string cmd) {
	return evalInR(toUTFz!(char*)(cmd));
}

void evalRQ(string cmd) {
	evalQuietlyInR(toUTFz!(char*)(cmd));
}

Robj RNil() {
	return Rf_allocVector(0, 1);
}
