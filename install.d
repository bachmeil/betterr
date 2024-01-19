import std.algorithm, std.array, std.file, std.path;
import std.process, std.stdio, std.string;

void main(string[] args) {
	assert(args.length > 1, "You need to supply an installation directory");
	string rinsideLocation = executeShell(`Rscript -e "cat(find.package('RInside'))"`).output;
	assert(!rinsideLocation.startsWith("Error "), "You do not have RInside installed. Install it and then run this program.");
	writeln(executeShell(`mkdir -p ` ~ args[1] ~ `/betterr`).output);
	foreach(src; sourceFiles) {
		copy(src, args[1] ~ "/betterr/" ~ src);
	}
	std.file.write(args[1] ~ "/Makefile", "LINK=-L" ~ rinsideLocation ~ "/lib/libRInside.so -L-lR\n\napp:\n\tldmd2 -i file.d $(LINK)\n");
}

string[] sourceFiles() {
	return dirEntries(".", "*.d", SpanMode.shallow)
        .filter!(a => a.isFile)
        .filter!(a => a != "install.d")
        .map!((return a) => baseName(a.name))
        .array;
}
