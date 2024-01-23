import betterr.ts, betterr.random, betterr.matrix, betterr.vector;
import betterr.r;
import std.stdio;

void main() {
  startR();
  
  auto ts = TS!(12)(rnorm(120), [1990, 4]);
  ts.print("ts");
  writeln(ts.frequency);
  
  auto ts2 = ts.lag(2);
  ts2.print("Second lag of ts");
  
  auto tp = TimePeriod(2012, 1, 12);
  writeln(tp - 12);
  writeln(tp + 12);
  writeln(tp - 5);
  writeln(tp + 60);
  writeln(tp + 59);
  writeln(tp + 58);
  writeln(tp + 61);
  writeln(tp + 62);
  writeln(tp < [2012, 1]);
  writeln(tp > [2012, 1]);
  writeln(tp == [2012, 1]);
  //~ writeln(tp + (-1));
  //~ ts.lag(3).print("Third lag");
  //~ ts.lead(1).print("First lead");
  //~ ts.diff(2).print("Time 2 difference");
  //~ ts.pct(1).print("Percent change");
  
  //~ writeln(ts2[1995,3]);
  //~ writeln(ts2[[1995,3]..[1998,1]]);
  //~ writeln(ts2.until([1995, 3]));
  //~ writeln(ts2.starting([1995, 3]));
  
  //~ MTS multiple = tsCombine([ts, ts2]);
  //~ multiple.print("Combining ts and ts2");
  
  closeR();
}
