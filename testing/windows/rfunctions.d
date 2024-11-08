struct Sexprec {}
alias Robj = Sexprec*;

struct RAPIBindings {
    extern(C) {
		double * REAL(Sexprec * x);
		int * INTEGER(Sexprec * x);
		void Rf_PrintValue(Robj x);
		Robj Rf_ScalarReal(double x);
		const(char) * R_CHAR(Sexprec * x);
    Sexprec * VECTOR_ELT(Sexprec * x, int i);
		Robj Rf_allocVector(uint type, int n);
		Robj Rf_allocMatrix(uint type, int rows, int cols);
		Robj Rf_allocList(int);
		Robj Rf_install(const char * sym);
  }
}
