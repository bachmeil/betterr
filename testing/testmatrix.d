import betterr.list, betterr.matrix;
import betterr.r;
import betterr.vector;
import std.conv, std.range, std.stdio;

void main() {
  startR();
  auto m1 = Matrix(2,2);
  m1[0,0] = 4.7;
  m1[0,1] = 4.8;
  m1[1,0] = 4.9;
  m1[1,1] = 5.0;
  auto m2 = Matrix(2,1);
  m2[0,0] = 5.1;
  m2[1,0] = 5.2;
  writeln(m1);
  writeln(m2);
  writeln(m1.matmul(m2));
  auto m3 = Matrix(3,2);
  Column(m3, 0) = [1.1, 2.2, 3.3];
  Column(m3, 1) = [4.4, 5.5, 6.6];
  writeln(m3);
  auto col0 = Column(m3, 0);
  foreach(ii, val; col0.enumerate) {
    col0[ii] = val*2;
  }
  writeln(m3);
  writeln(t(m3));
  writeln(m3.plus(2));
  writeln(m3.mul(2.7));
  writeln(solve(m1));
  writeln(inv(m1));
  writeln(diag(m1));
  auto svd = SVD(m1);
  writeln(svd.d);
  writeln(svd.u);
  writeln(svd.v);
  auto eigen = Eigen(m1);
  writeln(eigen.values);
  writeln(eigen.vectors);
  writeln(crossprod(m1));
  writeln(tcrossprod(m1));
  writeln(det(m1));
  writeln(diag(4));
  writeln(eye(4));
  writeln("m1: ", m1);
  Column(m1, 1) *= 1.2;
  writeln("m1: ", m1);
  Row(m1, 1) -= 0.7;
  writeln("m1: ", m1);
  
  auto m5 = Matrix(4,4);
  Column(m5, 0) = [1.1, 2.0, 3.0, 4.0];
  Column(m5, 1) = [5.0, 6.0, 7.0, 8.0];
  Column(m5, 2) = [9.0, 10.0, 11.0, 12.0];
  Column(m5, 3) = [13.0, 14.0, 15.0, 16.0];
  m5.print("m5");
  m5.Submatrix([0, 0], [2, 2]) = m5.Submatrix([2, 2], [4, 4]);
  m5.print("New m5"); 
  foreach(el; Elements(3,3)) {
    writeln(el);
  }
  foreach(el; Elements(4,2)) {
    writeln(el);
  }
  
  /* Multidimensional indexing */
  m5[_all, _all].print("The entire m5 matrix");
  m5[0..2, 0..2].print("[0,0] to [2,2] block");
  m5[0..2, _all].print("First two rows");
  m5[_all, 0..2].print("First two columns");
  m5[[0, 2], _all].print("First and third rows");
  m5[_all, [0, 2, 3]].print("First, third, and fourth columns");
  m5[1, [0, 2]].print("First and third elements of the second row");
  
  auto m6 = m5.dup();
  m6.print("Copy of m5");
  auto sm6 = m6.reference();
  sm6[0..2, 0..3] = 1.7;
  m6.print("m6 after filling in the upper left block");
  
  m6 = m5;
  m6[1, _all] = 9.4;
  m6.print("Set second row of m6 to 9.4");
  
  m6 = m5;
  m6[1, _all] = m6[0, _all];
  m6.print("Set the second row of m6 equal to the first");
  
  m6 = m5;
  m6[2..4, 2..4] = m6[0..2, 2..4];
  m6.print("Set the lower right block of m6 equal to the upper right block");
  
  m6[1..$, 2..$].print("Outer three rows and two columns of m6");
  
  /* Fill a column in a less error-prone way. Especially useful for
   * reusing an already allocated matrix. You cannot fill extra elements
   * (which will not necessarily be caught by other means) and you can
   * check that you've filled all elements. */
  auto cf = ColumnFill(m5, 2);
  cf ~= [1.2, 1.3, 1.4, 1.5];
  writeln("Has column 2 been filled completely? ", cf.full);
  m5.print("m5 after filling column 2");
  
  /* Map each row */
  m6.print("m6");
  double calcSum(Row r) {
    double result = 0.0;
    foreach(val; r) {
      result += val;
    }
    return result;
  }
  auto m7 = m6.mapRows!calcSum;
  m7.print("Sum of each row");
  
  double[] twice(Row r) {
    double[] result;
    foreach(val; r) {
      result ~= val*2.0;
    }
    return result;
  }
  auto m8 = m6.mapRows!twice;
  m8.print("Double the value of each element");
  closeR();
}
