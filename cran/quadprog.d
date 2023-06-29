module cran.quadprog;
import betterr.matrix, betterr.vector;
import std.algorithm.comparison, std.conv, std.stdio;

struct QPSolution {
	Vector sol;
	double value;
	Vector unconstrained;
	Vector lagrangian;
	int[] iterations;
	int[] active;
	
	void print() {
		sol.print("Solution vector");
		writeln("\nValue of the objective function: ", value);
		unconstrained.print("\nUnconstrained solution");
		lagrangian.print("\nLagrangian vector");
		writeln("\nIterations: ", iterations);
		writeln("Index of active constraints: ", active);
	}
}

QPSolution qpsolve(Matrix dmat, Vector dvec, Matrix amat, Vector bvec, int meq=0, bool factorized=false) {
  Matrix dmatcopy = Matrix(dmat);
  Vector dveccopy = Vector(dvec);
  Vector bveccopy = Vector(bvec);
	int n = dmat.rows.to!int;
	int fddmat = n;
	auto sol = Vector(n);
	auto lagr = Vector(amat.cols);
	double crval = 0.0;
	int fdamat = n;
	int q = amat.cols.to!int;
	auto iact = new int[q];
	iact[] = 0;
	int nact = 0;
	int[] iter = [0, 0];
	int r = min(n, q);
	auto work = Vector(2*n + r*(r+5)/2 + 2*q + 1);
	work = 0.0;
	auto ierr = to!int(factorized);
	
	qpgen2_(dmatcopy.ptr, dveccopy.ptr, &fddmat, &n, sol.ptr, lagr.ptr, &crval,
		amat.ptr, bvec.ptr, &fdamat, &q, &meq, iact.ptr, &nact, iter.ptr, work.ptr, &ierr);

	QPSolution result;
	result.sol = sol;
	result.value = crval;
	result.unconstrained = dvec;
	result.lagrangian = lagr;
	result.iterations = iter;
	result.active = iact;
	return result;
}

extern(C) {
	void qpgen2_(double *dmat, double *dvec, int *fddmat, int *n,
      double *sol, double *lagr, double *crval,
      double *amat, double *bvec, int *fdamat, int *q,
      int *meq, int *iact, int *nact, int *iter, double *work, int *ierr);
}
