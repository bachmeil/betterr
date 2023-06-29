import betterr.array, betterr.random, betterr.matrix, betterr.vector;
import betterr.r;
import std.stdio;

void main() {
  startR();
  auto arr = RArray(2,3,4,2);
  foreach(ii; 0..48) {
    arr.ptr[ii] = double(ii);
  }
  printR(arr.x);
  writeln(arr[1,1,1,1]);
  writeln(arr[1,2,3,1]);
  writeln(arr[1,1,3,0]);
  arr[1,1,1,1] = 33.7;
  printR(arr.x);
  arr = 2.7;
  printR(arr.x);
  arr = rnorm(48);
  printR(arr.x);
  auto arr2 = arr.sub[0, 0, 0, _all];
  printR(arr2.x);
  auto arr3 = arr.sub[0, 0..3, 0, _all];
  printR(arr3.x);
  arr.sub[0..2, 0, 0, 0] = arr.sub[0..2, 1, 1, 1];
  printR(arr.x);
  auto m = Matrix([1.1, 2, 3, 4, 5, 6], 2, 3);
  arr.sub[0..2, 0..3, 1, 1] = m;
  writeln(arr);
  closeR();
}
