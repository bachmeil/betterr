import betterr.ts, betterr.random, betterr.matrix, betterr.vector;
import betterr.r;
import std.stdio;

void main() {
  startR();
  
  auto ts = TS(rnorm(120), [1990, 4], 12);
  writeln(ts.frequency);
  writeln(ts.fromLong(ts.start));
  
  auto ts2 = ts.lag(2);
  writeln(ts2.frequency);
  ts2.print("Second lag");
  ts.lag(3).print("Third lag");
  ts.lead(1).print("First lead");
  ts.diff(2).print("Time 2 difference");
  ts.pct(1).print("Percent change");
  
  writeln(ts2[1995,3]);
  writeln(ts2[[1995,3]..[1998,1]]);
  writeln(ts2.until([1995, 3]));
  writeln(ts2.starting([1995, 3]));
  
  closeR();
}
