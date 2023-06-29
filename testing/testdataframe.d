import betterr.dataframe, betterr.vector;
import betterr.r;
import std.conv, std.stdio;

void main() {
  startR();
	auto df = DataFile("unrate.csv");
	DataFrame unrate = df.read();
	unrate.print("Unemployment rate, monthly, Jan 2022-Jan 2023");
	printR(unrate.classInfo);
	Vector u = unrate["UNRATE"];
	u.print("Just the unemployment rate");
	printR(u.classInfo);
	DataFrame u2 = unrate[["DATE", "UNRATE"]];
	u2.print("This is a data frame");
	DataFrame u3 = unrate[["UNRATE"]];
	writeln("u3.name: ", u3.name);
	writeln(u3.data.name);
	u3.print("This is also a data frame");
	writeln("Columns in u3: ", u3.ncol);
	writeln("Rows in u3: ", u3.nrow);
	writeln(u2.ncol);
	writeln(u2.length);
	writeln("Names of u2: ", u2.names);
	
	auto v = Vector([1.1, 2, 3, 4, 5, 6]);
	auto cat = Cat("testfile.txt");
	cat.meta = "This is test data";
	cat.sep = "\n";
	cat.apply(v);
	closeR();
}
