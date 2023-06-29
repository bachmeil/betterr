/* This module is of limited use. It's a thin wrapper over R code.
 * It's easier to write the R code to do your plotting directly, given
 * all the options you have for plotting. */
module betterr.plot;
import betterr.r;
import std.path, std.stdio;

struct Plot {
	string name;
	string type;
	string main;
	string sub;
	string xlab;
	string ylab;
	
	this(string _name) {
		name = _name;
	}
	
	this(T)(T x) {
		name = x.name;
	}
	
	void create(string fn) {
		string cmd = `pdf("` ~ setExtension(fn, "pdf") ~ `"); plot(` ~ name;
		if (type.length > 0) {
			cmd ~= `, type='` ~ type ~ `'`;
		}
		if (main.length > 0) {
			cmd ~= `, main='` ~ main ~ `'`;
		}
		if (sub.length > 0) {
			cmd ~= `, sub='` ~ sub ~ `'`;
		}
		if (xlab.length > 0) {
			cmd ~= `, xlab='` ~ xlab ~ `'`;
		}
		if (ylab.length > 0) {
			cmd ~= `, ylab='` ~ ylab ~ `'`;
		}
		cmd ~= "); dev.off()";
		writeln(cmd);
		evalRQ(cmd);
	}
}
