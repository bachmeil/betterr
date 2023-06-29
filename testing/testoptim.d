import betterr.optim, betterr.r, betterr.vector;
import std.stdio;

extern(C) {
  double f(int n, double * par, void * ex) {
    return par[0]*par[0] + par[1]*par[1];
  }

  void g(int n, double * par, double * gr, void * ex) {
    gr[0] = 2*par[0];
    gr[1] = 2*par[1];
  }
}

void main() {
  startR();

  auto nm = NelderMead(&f);
  OptimSolution sol = nm.solve([3.5, -5.5]);
  sol.print();
  sol = nm.solve(Vector([3.5, -5.5]));
  sol.print();
  double[] starting = [3.5, -5.5];
  double[] solution = [3.5, -5.5];
  sol = nm.solve(solution.ptr, starting.ptr, 2);
  sol.print();
  
  auto bfgs = BFGS(&f, &g);
  starting = [3.5, -5.5];
  OptimSolution sol2 = bfgs.solve(starting);
  writeln();
  sol2.print;

  auto cg = ConjugateGradient(&f, &g);
  starting = [3.5, -5.5];
  OptimSolution sol3 = cg.solve(starting);
  writeln();
  sol3.print;

  auto bounded = Bounded(&f, &g);
  starting = [3.5, -5.5];
  // Impose an upper bound on the first parameter, no bound on the second
  OptimSolution sol4 = bounded.solve(starting, [-10.0, -10.0], [-5.0, 10.0], [3, 0]);
  writeln();
  sol4.print;

  auto sa = SA(&f);
  starting = [3.5, -5.5];
  writeln("Starting sann");
  OptimSolution sol5 = sa.solve(starting);
  writeln("Finished with sann");
  sol5.print;

  auto sa2 = SA(&f);
  sa2.maxit = 100_000;
  auto starting2 = Vector([3.5, -5.5]);
  writeln("Starting sann");
  OptimSolution sol6 = sa2.solve(starting);
  writeln("Finished with sann");
  sol6.print;

  closeR();
}

