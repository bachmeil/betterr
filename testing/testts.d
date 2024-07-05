import betterr.ts, betterr.random, betterr.matrix, betterr.vector;
import betterr.r;
import std.conv, std.stdio;

void main() {
  startR();
  
  auto ts = TS!12(rnorm(120), [1990, 4]);
  auto ts2 = ts.lag(2);
  MultipleTS!12 mts;
  mts["y"] = ts;
  mts["x"] = ts2;
  auto mts2 = mts[];
  mts2.print("Both series");
  
  auto mts3 = MultipleTS!12(["y": ts, "x": ts2]);
  mts3[].print("Both series again");
  
  auto mts4 = MTS!12(["y": ts, "x": ts2]);
  mts4.print("And again");
  
  MTS!12 mts5 = ["y": ts, "x": ts2];
  mts5.print("And yet again");
  
  MultipleTS!12 mts6 = ["y": ts, "x": ts2];
  writeln("One more time: ", mts6);
  
  mts6.to!(MTS!12).print("Another time");
  
  // Don't worry about the names
  // Just create a MTS struct with the data
  MTS!12 mts7 = [ts, ts2];
  mts7.print("Without providing names");
  
  mts5["y"].print("y value");
  mts5[0].print("y value by index");
  writeln(mts5.names);
  writeln(mts6.keys);
  
  closeR();
}
