import betterr.dataframe, betterr.list;
import cran.vars;
import betterr.r;
import std.conv, std.stdio;

void main() {
  startR();
  auto df = DataFile("macrodata.csv");
  df.commentChar = "#";
  auto dataset = df.read();
  VAR varconf = {p: 4, type: "both"};
  auto fit = varconf.fit(dataset);
  fit.print("Estimated VAR model");
  closeR();
}
