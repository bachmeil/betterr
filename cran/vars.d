module cran.vars;
import betterr.r, betterr.list;
import std.conv;
import std.stdio;

struct VAR {
  long p = 1;
  string type = "const";
  long lagMax = -1;
  string ic = "AIC";
  
  List fit(T)(T x) {
    string cmd = `vars::VAR(` ~ x.name ~ `, p=` ~ p.to!string ~ `, type="` ~ type ~ `"`;
    if (lagMax > 0) {
      cmd ~= `, lag.max=` ~ lagMax.to!string;
    } 
    cmd ~= `, ic="` ~ ic ~ `")`;
    writeln(cmd);
    return List(cmd);
  }
}
