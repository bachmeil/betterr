import betterr.baser, betterr.list, betterr.vector, betterr.matrix;
import betterr.r;
import std.conv, std.stdio;

void main() {
  startR();
  auto rlist = List(0);
  rlist["description"] = "This R list is part of my D program.";
  rlist["description"].print("description");
  rlist["scalar double"] = 3.7;
  rlist["scalar double"].print("Scalar (double)");
  rlist["scalar int"] = 4;
  rlist["scalar int"].print("Scalar (int)");
  rlist["scalar bool"] = true;
  rlist["vector"] = Vector([-6.2, 4.8, 3.7]);
  rlist["vector"].print("Vector");
  auto m = Matrix(seq(0.5, 4.5, 0.5), 3, 3);
  rlist["double matrix"] = m;
  rlist.print("\nTest of list functionality");
  rlist[0..3].print("\nFirst three elements of this list");
  rlist.dup.print("\nCopy of the list");
  rlist[2..$].print("\nDropping the first two elements of the list");
  string s1 = rlist[0].as!string;
  writeln("s1: ", s1);
  writeln(s1.length);
  double x1 = rlist[1].as!double;
  writeln("x1: ", x1);
  int x2 = rlist[2].as!int;
  writeln("x2: ", x2);
  bool b1 = rlist[3].as!bool;
  writeln("b1: ", b1);
  auto m1 = rlist["double matrix"].as!Matrix;
  m1.print("Copy of the matrix");
  auto v1 = rlist["vector"].as!Vector;
  v1.print("Copy of the vector");
  closeR();
}
  
