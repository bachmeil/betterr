module betterr.dataframe;
import betterr.rdata, betterr.vector;
import betterr.r;
import std.array, std.conv, std.file, std.path, std.stdio;

struct DataFrame {
  RData data;
  alias data this;
  
  this(string code) {
    data = RData(code);
  }
  
  this(DataFrame df) {
    data = this(df.name);
  }
  
  this(Vector v) {
		this("data.frame(" ~ v.data.name ~ ")");
	}
  
  /* A string returns a vector */
  Vector opIndex(string n) {
    return Vector(data.name ~ "[,'" ~ n ~ "']");
  }

	/* A string array returns a data frame, even if there's only one element.
	 * This is an example of silent conversions in R that drive users nuts.
	 * The static typing of D prevents that from happening, and that's a
	 * good thing. */
  DataFrame opIndex(string[] ns) {
		if (ns.length == 1) {
			auto result = DataFrame(opIndex(ns[0]));
			result.setNames(ns[0]);
			return result;
		} else {
			return DataFrame(data.name ~ "[," ~ ns.stringArray ~ "]");
		}
  }
  
  Vector column(long ii) {
    return Vector(data.name ~ "[," ~ ii.to!string ~ "]");
  }
  
  Vector row(long ii) {
    return Vector(data.name ~ "[" ~ ii.to!string ~ ",]");
  }

  Vector column(string name) {
    return Vector(data.name ~ "[," ~ name ~ "]");
  }

	void setNames(string[] names) {
		string cmd = "names(" ~ data.name ~ ") <- " ~ names.stringArray;
		evalRQ(cmd);
		/* This is a gotcha - R might change the pointer when doing this,
		 * so we need to update data. */
		this.update;
	}
	
	long ncol() {
		return data.apply("ncol").scalar!long;
	}
	
	long nrow() {
		return data.apply("nrow").scalar!long;
	}
	
	long length() {
		return data.apply("length").scalar!long;
	}
	
	string[] names() {
		return betterr.r.stringArray(data.apply("names"));
	}
	
	void setNames(string name) {
		setNames([name]);
	}

  void print(string msg="") {
    if (msg != "") { writeln(msg, ":"); }
    printR(data.x);
  }
}

/* What follows was an early attempt at wrapping certain R functionality.
 * I eventually gave up because this wasn't any better than writing the
 * R code directly and then pulling the data into D. I've left it in
 * case I decide to return to it in the future. */
struct DataFile {
	string name;
	bool header = true;
	long skip = 0;
	string sep = ",";
	string naStrings = "NA";
	string commentChar = "";

	DataFrame read() {
    string cmd = "read.csv('" ~ name ~ "'";
    if (header) {
      cmd ~= ", header=TRUE";
    } else {
      cmd ~= ", header=FALSE";
    }
    if (skip > 0) {
      cmd ~= ", skip=" ~ skip.to!string;
    }
    if (sep != ",") {
      cmd ~= ", sep='" ~ sep ~ "'";
    }
    if (naStrings != "NA") {
      cmd ~= ", na.strings='" ~ naStrings ~ "'";
    }
    if (commentChar != "") {
      cmd ~= ", comment.char='" ~ commentChar ~ "'";
    }
    cmd ~= ")";
    writeln(cmd);
    if (extension(name) == ".csv") {
      return DataFrame(cmd);
    } else {
      return DataFrame();
    }
	}
}

string stringArray(string[] arr) {
	return `c("` ~ arr.join(`","`) ~ `")`;
}

struct Cat {
	string file;
	string sep = " ";
	string meta;
	bool append;
	
	void apply(string code) {
		string cmd = "cat(" ~ code;
		if (file.length > 0) {
			cmd ~= ", file=\"" ~ file ~ "\"";
		}
		cmd ~= ", sep=\"" ~ sep ~ "\"";
		if (append) {
			cmd ~= ", append=TRUE";
		}
		cmd ~= ")";
		evalRQ(cmd);
		if (meta.length > 0) {
			std.file.write("meta-" ~ file, meta);
		}
	}
	
	void apply(T)(T x) {
		apply(x.name);
	}
} 
	
	
	
	
	
	
	
	
	
