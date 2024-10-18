/* This file enables installation in a directory for compilation using 
 * the -i option. It would be straightforward to create a Dub package,
 * but I prefer to keep the dependencies inside the project repo in
 * case there are future updates.
 * 
 * You use it by going to this directory (the cloned repo) and running
 * this script using rdmd. Example:
 * 
 * rdmd ./install.d /home/gunther/myrepo
 * 
 * Then you compile your source file in /home/gunther/myrepo like this:
 * 
 * ldmd2 -i file.d
 * 
 * (Note that I don't use Dub to build my programs.) To update your
 * installation, use the update.d script instead. This is for the
 * creation of a new betterr project.
 * 
 * You can optionally add the name of your source file as an argument.
 * It will use that rather than file.d as the name of your source file.
 */
import std.algorithm, std.array, std.file, std.path;
import std.process, std.stdio, std.string;

void main(string[] args) {
	assert(args.length > 1, "You need to supply an installation directory");
	string rinsideLocation = executeShell(`Rscript -e "cat(find.package('RInside'))"`).output;
	assert(!rinsideLocation.startsWith("Error "), "You do not have RInside installed. Install it and then run this program.");
  
  // Copy the betterr core source files
	writeln(executeShell(`mkdir -p ` ~ args[1] ~ `/betterr`).output);
	foreach(src; sourceFiles) {
		copy(src, args[1] ~ "/betterr/" ~ src);
	}

  // Copy the dstats source file
	writeln(executeShell(`mkdir -p ` ~ args[1] ~ `/dstats`).output);
	copy("dstats/summary.d", args[1] ~ "/dstats/summary.d");
  
	// Copy the Gretl matrix files
	writeln(executeShell(`mkdir -p ` ~ args[1] ~ `/gretl/matrix`).output);
	writeln(executeShell(`cp gretl/matrix/* ` ~ args[1] ~ `/gretl/matrix/`).output);
  
  // Copy the GSL RNG files
	writeln(executeShell(`mkdir -p ` ~ args[1] ~ `/gsl/rng`).output);
	writeln(executeShell(`cp gsl/rng/* ` ~ args[1] ~ `/gsl/rng/`).output);

  // Create the Makefile
  string appFilename = "file";
  if (args.length > 2) {
    appFilename = setExtension(args[2], "d");
  }
	if (!exists(args[1] ~ "/Makefile")) {
		std.file.write(args[1] ~ "/Makefile", "LINK=-L" ~ rinsideLocation ~ "/lib/libRInside.so -L-lR -P-Igsl/rng gsl/rng/*.c -P-Igretl/matrix gretl/matrix/*.c gretl/matrix/*.d -L-lopenblas\n\napp:\n\tldmd2 -i " ~ appFilename ~ " $(LINK)\n");
	}
  if (!exists(args[1] ~ "/" ~ appFilename)) {
    std.file.write(args[1] ~ "/" ~ appFilename, "");
  }
}

string[] sourceFiles() {
	return dirEntries(".", "*.d", SpanMode.shallow)
        .filter!(a => a.isFile)
        .filter!(a => (a != "install.d") && (a != "update.d"))
        .map!((return a) => baseName(a.name))
        .array;
}
