import mir.stat.descriptive.univariate;
import std.stdio;

void main() {
  writeln(mean([1.1, 2.2, 3.3]));
  double[] x = [];
  writeln(mean(x));
  writeln(mean([1.1, double.nan, 2.2, 3.3]));
}
