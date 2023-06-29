import betterr.vector, betterr.matrix;
import betterr.r;
import std.conv, std.math, std.stdio;
import std.range: enumerate;

void main() {
  startR();
  
  auto x = Vector(5);
  x[0] = 1.5;
  x[3] = 4.1;
  x.print("Test vector");
  
  auto y = Vector("rnorm(25)");
  y.print();
  
  auto y2 = Vector(y.length);
  foreach(ii; 0..y.length) {
    y2[ii] = y[ii]*y[ii];
  }
  y2.print();
  y2[4..10].print();
  y2[].print();
  y2[[0,2,4,6,8]].print();
  
  writeln(y2.to!(double[]));
  
  /* Fill a vector in a way that is less likely to lead to error. */
  auto v1 = Vector(5);
  auto f1 = Fill(v1);
  f1 ~= 1.1;
  writeln("Is v1 full? ", f1.full);
	v1.print("v1 is not filled");
  f1 ~= 1.1;
  f1 ~= 1.1;
  f1 ~= 1.1;
  f1 ~= 1.1;
  writeln("Is v1 full? ", f1.full);
	v1.print("v1 after filling");
  // Error if you run this line
  //~ f1 ~= 1.1;
  
  //~ Now reuse v1 to prevent reallocation, but without accidentally
  //~ filling the wrong number of elements.
  f1.reset;
  //~ Error if you run this line
  //~ f1 ~= [2.2, 2.2, 2.2, 2.2, 2.2, 2.2];
  f1 ~= [2.2, 2.2, 2.2, 2.2, 2.2];
  writeln("Is v1 full again? ", f1.full);
  v1.print("v1 after refilling");
  
  /* Can also partially fill a vector with error checking. Useful if you
   * only want to update part of the data. */
  auto f2 = PartialFill(v1, 1, 4);
  f2 ~= 3.3;
  f2 ~= 3.3;
  f2 ~= 3.3;
  //~ Error if you run this line
  //~ f2 ~= 3.3;
  writeln("Have we filled all of elements 1 through 3? ", f2.full);
  v1.print("v1 after filling elements 1 through 3");
  
  f2.reset;
  //~ Error if you run this line
  //~ f2 ~= [4.4, 4.4, 4.4, 4.4];
  f2 ~= [4.4, 4.4, 4.4];
  writeln("Have we refilled all of elements 1 through 3? ", f2.full);
  v1.print("v1 after refilling elements 1 through 3");
  
  //~ Can use a Vector as well
  f2.reset;
  f2 ~= Vector([5.5, 5.5, 5.5]);
  writeln("Have we refilled all of elements 1 through 3? ", f2.full);
  v1.print("v1 after refilling elements 1 through 3");
  
  //~ And a Row or Column of a Matrix
  auto mat1 = Matrix(3, 3);
  Row(mat1, 0) = [1.2, 1.3, 1.4];
  f2.reset;
  f2 ~= Row(mat1, 0);
  writeln("Have we refilled all of elements 1 through 3? ", f2.full);
  v1.print("v1 after refilling elements 1 through 3");
  
  Column(mat1, 2) = [1.5, 1.6, 1.7];
  f2.reset;
  f2 ~= Column(mat1, 2);
  writeln("Have we refilled all of elements 1 through 3? ", f2.full);
  v1.print("v1 after refilling elements 1 through 3");
  
  auto iv = IntVector(5);
  iv = [1, 2, 3, 4, 5];
  iv.print("An integer vector");
  
  auto bv = BoolVector([true, false, false, true, true]);
  bv.print("A boolean vector");
  
  auto sv = StringVector(["Foo", "Bar", "Baz", "D is kind of cool!"]);
  sv.print("A string vector");
  writeln(sv[1]);
  
  foreach(val; Vector([1.1, 2.2, 3.3])) {
		writeln(val);
	}
	
  foreach(ii, val; Vector([1.1, 2.2, 3.3]).enumerate) {
		writeln(ii+1, " ", val);
	}
  
  /* You can use map in a generic fashion
   * If the function returns a double, it returns a Vector
   * If the function returns a double[], it returns a Matrix
   */
  auto v2 = Vector([1.1, 2.2, 3.3, 4.4]);
  auto v3 = v2.mapRows!(a => a^^2);
  v3.print("Square each element");
	
  auto m3 = v2.mapRows!(a => [a^^2, a^^3]);
  m3.print("Square and cube each element");
  
  auto v4 = v2.mapRows!(a => (a > 3.0));
  v4.print("Check if above 3");
  
  auto v5 = v2.mapRows!(a => lround(a));
  v5.print("Convert to int");
  
  closeR();
}
