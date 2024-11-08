import makebindings;
import std;

void main() {
	writeln(makeBindings!([["RAPIBindings", "libr"], ["RInsideBindings", "rinside"]])());
}
