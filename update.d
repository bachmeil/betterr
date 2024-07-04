/* Usage:
 * rdmd update.d [path] [files (optional)]
 * 
 * Examples:
 * 
 * // Update all betterr core files
 * rdmd update.d /home/phil/test 
 * 
 * // Update only ts.d
 * rdmd update.d /home/phil/test ts
 * 
 * // Update ts.d and array.d
 * rdmd update.d /home/phil/test ts array */
import std.algorithm, std.array, std.file, std.path;
import std.process, std.stdio, std.string;

void main(string[] args) {
	assert(args.length > 1, "You need to supply an installation directory");
	if (args.length == 2) { 
		writeln("Updating all betterr core files");
		foreach(src; sourceFiles) {
			string fn = args[1] ~ "/betterr/" ~ setExtension(src, "d");
			if (exists(fn)) {
				remove(fn);
			}
			copy(src, fn);
		}
	} else {
		writeln("Updating only these core files: ", args[2..$]);
		foreach(src; args[2..$]) {
			string fn = args[1] ~ "/betterr/" ~ setExtension(src, "d");
			writeln(fn);
			if (exists(fn)) {
				remove(fn);
			}
			copy(setExtension(src, "d"), fn);
		}
	}
}

string[] sourceFiles() {
	return dirEntries(".", "*.d", SpanMode.shallow)
        .filter!(a => a.isFile)
        .filter!(a => (a != "install.d") && (a != "update.d"))
        .map!((return a) => baseName(a.name))
        .array;
}
