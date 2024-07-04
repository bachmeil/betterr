import std.algorithm, std.array, std.file, std.path;
import std.process, std.stdio, std.string;

/* Add these options:
 * 
 * --full: Install everything
 * --link: Create a symbolic link rather than copying (updates are automatic but might break your code)
 * --dir: Source directory (use the current directory if not specified)
 * --target: Target directory (required if no dir option; optional otherwise, using the current directory if not specified)
 * --force: Don't do any checks for existing directories or files; things may possibly get overwritten
 */

void main(string[] args) {
	assert(args.length > 1, "You need to supply an installation directory");
	string rinsideLocation = executeShell(`Rscript -e "cat(find.package('RInside'))"`).output;
	assert(!rinsideLocation.startsWith("Error "), "You do not have RInside installed. Install it and then run this program.");
	writeln(executeShell(`mkdir -p ` ~ args[1] ~ `/betterr`).output);
	foreach(src; sourceFiles) {
		copy(src, args[1] ~ "/betterr/" ~ src);
	}
	// This works for now
	writeln(executeShell(`mkdir -p ` ~ args[1] ~ `/dstats`).output);
	copy("dstats/summary.d", args[1] ~ "/dstats/summary.d");
	// This works for now
	writeln(executeShell(`mkdir -p ` ~ args[1] ~ `/gretl/matrix`).output);
	writeln(executeShell(`mkdir -p ` ~ args[1] ~ `/gsl/rng`).output);
	writeln(executeShell(`cp gretl/matrix/* ` ~ args[1] ~ `/gretl/matrix/`).output);
	writeln(executeShell(`cp gsl/rng/* ` ~ args[1] ~ `/gsl/rng/`).output);
	if (!exists(args[1] ~ "/Makefile")) {
		std.file.write(args[1] ~ "/Makefile", "LINK=-L" ~ rinsideLocation ~ "/lib/libRInside.so -L-lR\n\napp:\n\tldmd2 -i file.d $(LINK)\n");
	}
}

string[] sourceFiles() {
	return dirEntries(".", "*.d", SpanMode.shallow)
        .filter!(a => a.isFile)
        .filter!(a => (a != "install.d") && (a != "update.d"))
        .map!((return a) => baseName(a.name))
        .array;
}
