/* Windows has issues with the search path for DLLs.
 * And TBH I have no idea how any of it works even after researching
 * it for two hours.
 * The only thing I've found that works every time is to copy all the
 * DLLs, along with their dependencies, into the same directory as
 * the exe produced by the D compiler. */
public import rfunctions, rinsidefunctions;
import std;

string makeBindings(string[2][] info)() {
    string result = "extern(C) {\n";
    string result2;
    static foreach(fs; info) {
			result2 ~= `void bind_` ~ fs[1] ~ `(string dllPath, string dllName, bool errors=false) {
import std.process, std.stdio, std.utf;
import core.sys.windows.windows;
import std.windows.syserror;

// In some cases (RInside) I had to do this for it to work
// environment["PATH"] = dllPath;

// Error 126 means *something* can't be found, including a dependency.
HMODULE ` ~ fs[1] ~ ` = LoadLibraryW(toUTF16z(dllPath ~ "/" ~ dllName ~ ".dll"));
if (errors) {
writeln(sysErrorString(GetLastError()));
}
`;

        foreach(memberName; __traits(allMembers, mixin(fs[0]))) {
            //pragma(msg, memberName);
            //pragma(msg, ReturnType!(__traits(getMember, mixin(fs), memberName)).stringof);
            //pragma(msg, Parameters!(__traits(getMember, mixin(fs), memberName)).stringof);
            string ret = ReturnType!(__traits(getMember, mixin(fs[0]), memberName)).stringof;
            string par = Parameters!(__traits(getMember, mixin(fs[0]), memberName)).stringof;
            result ~= ret ~ " function" ~ par ~ " " ~ memberName ~ ";\n";
            result2 ~= memberName ~ " = cast(typeof(" ~ memberName ~ ")) GetProcAddress(" ~ fs[1] ~ ", \"" ~ memberName ~ "\");\n";
        }
        result2 ~= "}\n\n";
    }
    return result ~ "}\n\n" ~ result2;
}

//enum values = makeBindings!([["Bindings", "libr"], ["Bindings2", "rinside"]])();

//~ void main() {
	//~ writeln(values);
//~ }
//~ struct Sexprec {}
