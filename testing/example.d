import betterr.r, betterr.vector;
import std.stdio;

void main() {
  startR();
  auto v = Vector([1.1, 2.2, 3.3]);
  writeln(v);
  foreach(ii; 0..3) {
    v[ii] = v[ii]*2.5;
  }
  writeln(v);
  closeR();
}
