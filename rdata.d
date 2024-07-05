module betterr.rdata;
import betterr.r;
import betterr.matrix, betterr.vector, betterr.array, betterr.list, betterr.dataframe;
import std.conv, std.exception, std.stdio, std.string;

/* This module includes the data type that holds any SEXP struct in R.
 * Therefore, it is used to hold data for a vector, matrix, list, data
 * frame, or anything else.
 * 
 * This module will be imported by any other module that wraps an R
 * data type.
 * 
 * We hold a pointer to the data because we may, for instance, want to
 * get or modify elements of data inside R, and that would be much slower
 * if it was done using text snippets. */
struct RDataStorage {
	string name;
  Robj x;
  int refcount;
}
 
struct RData {
  RDataStorage * data;
  alias data this;
  /* This is the code used to tell R to remove the variable.
   * I previously used string concatenation to do that inside the
   * destructor, but the program died with an uninformative error
   * message. We can avoid allocation by creating this array. We know
   * the name is always 23 characters, calling rm adds 4, and it has
   * to be null terminated, so it's always 28 characters. Send release.ptr
   * to R. */
  private char[28] release;
  
  this(string code) {
		//writeln(code); // Useful for debugging, since this is where many errors occur
		data = new RDataStorage;
    import std.datetime;
    // This occasionally gives confusing errors if you have a large number of zeros
    // Doing this avoids those errors
    string ct = std.conv.to!string(Clock.currTime().toISOString.leftJustifier(23, '0'));
    data.name = "x" ~ ct[0..15] ~ ct[16..$];
    data.x = evalR(name ~ ` <- ` ~ code ~ `;` ~ name);
    data.refcount = 1;
    
    // See note above about release
    release[0..3] = ['r', 'm', '('];
    foreach(ii, char ch; name) {
      release[ii+3] = name[ii];
    }
    release[26] = ')';
    release[27] = '\0';
  }
  
  this(RData rd) {
    data = new RDataStorage;
    data.name = rd.name;
    data.x = rd.x;
    data.refcount = refcount+1;
  }
  
	Robj apply(string s) {
		return evalR(s ~ "(" ~ name ~ ")");
	}
	
  Robj classInfo() {
		return apply("class");
	}

	/* name won't change, but x might change silently. Have to update x
	 * after changes to the attributes of x. Putting it out here because
	 * it's passed by reference. */	
	void update() {
		data.x = evalR(data.name);
	}

  this(this) {
    enforce(x !is null, "x should never be null inside an RData struct. Are you sure you called the constructor?");
    refcount += 1;
  }
	
	/* Since D does not have default constructors for structs, and I don't
	 * have a desire to implement a workaround just for this, we'll only
	 * go forward with removal inside R if this is not null. If it's null,
	 * then clearly there's nothing on the R side of the world that needs
	 * to be cleared from memory.
	 * 
	 * Note that we don't need to mess with a check for data.x being null,
	 * since we don't work with the pointer, and in fact the pointer is
	 * not of any interest here. The use of pointers is more of an 
	 * optimization that holds a reference to temporary data. 
   * 
   * Note: Don't use string concatenation inside here. I did and my
   * program died with an invalid memory operation error. You don't
   * want to allocate in here. I've made this @nogc. */
  @nogc ~this() {
		if (data !is null) {
			refcount -= 1;
			if (refcount == 0) {
				evalQuietlyInR(release.ptr);
			}
		}
  }
  
  /* Only call this if you know what you're doing, which is probably never. 
   * Add a message for debugging. */
  void free(bool printmsg = false)(string msg="") {
		static if(printmsg) {
			writeln("Freeing RData manually: ", msg);
			if (data is null) {
				writeln("data is null");
			} else {
				writeln(this);
			}
		}
		evalRQ(`rm(` ~ name ~ `)`);
		data = null;
	}
  
  /* Use `as` rather than `to` because toString messes things up for strings. */
  T as(T)() {
    static if(is(T == double)) {
      return data.x.scalar;
    }
    static if(is(T == int)) {
      return data.x.scalar!int;
    }
    static if(is(T == long)) {
      return data.x.scalar!long;
    }
    static if(is(T == bool)) {
      return data.x.scalar!int.to!bool;
    }
    static if(is(T == string)) {
      return data.x.scalar!string;
    }
    static if(is(T == Matrix)) {
      return Matrix(this);
    }
    static if(is(T == Vector)) {
      return Vector(this);
    }
    static if(is(T == List)) {
      return List(this);
    }
    static if(is(T == DataFrame)) {
      return DataFrame(this);
    }
  }
  
  string toString() {
    printR(data.x);
    return "";
  }
  
  void print(string msg="") {
    if (msg != "") { writeln(msg, ":"); }
    printR(data.x);
  }
}
