import cran.quadprog;
import betterr.matrix, betterr.vector;
import betterr.r;
import std.stdio;

/* This is the example in the manual for solve.QP */
void main() {
  startR();

  writeln("D version:");
  Matrix dmat = eye(3);
  auto dvec = Vector([0.0, 5.0, 0.0]);
  auto amat = Matrix(3,3);
  Row(amat, 0) = [-4.0, 2.0,  0.0];
  Row(amat, 1) = [-3.0, 1.0, -2.0];
  Row(amat, 2) = [ 0.0, 0.0,  1.0];
  auto bvec = Vector([-8.0, 2.0, 0.0]);
  QPSolution sol = qpsolve(dmat, dvec, amat, bvec);
  sol.print();
  
  writeln("\nR version:");
  evalRQ(["library(quadprog)",
    "Dmat       <- matrix(0,3,3)",
    "diag(Dmat) <- 1",
    "dvec       <- c(0,5,0)",
    "Amat       <- matrix(c(-4,-3,0,2,1,0,0,-2,1),3,3)",
    "bvec       <- c(-8,2,0)",
    "print(solve.QP(Dmat,dvec,Amat,bvec=bvec))"]);
  
  closeR();
}
