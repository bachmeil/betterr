module betterr.array;
import betterr.rdata;
import betterr.r;
import betterr.matrix, betterr.vector;
import std.conv, std.exception, std.range, std.stdio;

struct TS(int freq) {
  /* The value of _start and _end depend on the frequency. I don't know
   * if we should require the frequency as a compile-time parameter,
   * but it does make some things easier. */
  RData data;
  double * ptr;
  int frequency = freq;
  alias data this;
  
  static if(freq == 1) {
		private long _start;
		private long _end;
		long start() {
			return _start;
		}
		long end() {
			return _end;
		}
		void start(long s) {
			_start = s;
		}
		void end(long e) {
			_end = e;
		}
	} else {
		private long[] _start;
		private long[] _end;
		long[] start() {
			return _start;
		}
		long[] end() {
			return _end;
		}
		void start(long[] s) {
			_start = s;
		}
		void end(long[] e) {
			_end = e;
		}
		long start(T: long)() {
			return frequency*(_start[0] - 1900) + (_start[1]-1);
		}
		long end(T: long)() {
			return frequency*(_start[0] - 1900) + (_start[1]-1);
		}
	}

  long asLong(long[] d) {
    if (frequency != 1) {
      return (d[0]-1900)*frequency+(d[1]-1);
    } else {
      return d[0];
    }
  }
  
  long[] fromLong(long d) {
    return [(d/frequency)+1900, d%frequency+1];
  }
  
  long asLong(long[] d, long f) {
    if (f != 1) {
      return (d[0]-1900)*f+(d[1]-1);
    } else {
      return d[0];
    }
  }
  
  long[] fromLong(long d, long f) {
    return [(d/f)+1900, d%f+1];
  }
  
  /* To access an existing ts object that's already inside R */
  this(string code) {
    data = RData(code);
    ptr = REAL(data.x);
    frequency = freq;
    auto tmpStart = IntVector("as.integer(start(" ~ data.name ~ "))");
    auto tmpEnd = IntVector("as.integer(end(" ~ data.name ~ "))");
    static if(freq == 1) {
			_start = tmpStart[0];
			_end = tmpEnd[0];
		} else {
			_start = [tmpStart[0], tmpStart[1]];
			_end = [tmpEnd[0], tmpEnd[1]];
		}
  }
  
  /* For an annual ts */
  //~ this(string code, long _start) {
    //~ data = RData("ts(" ~ code ~ ", start=" ~ _start.to!string ~ ", frequency=1)");
    //~ ptr = REAL(data.x);
    //~ frequency = 1;
    //~ start = this.start;
    //~ end = this.start+data.x.length;
  //~ }

  /* Monthly or quarterly */
  //~ this(string code, long[2] _start) {
    //~ data = RData("ts(" ~ code ~ ", start=c(" ~ _start[0].to!string ~ ", " 
      //~ ~ _start[1].to!string ~ "), frequency=" ~ _frequency.to!string ~ ")");
    //~ ptr = REAL(data.x);
    //~ frequency = _frequency;
    //~ start = asLong(_start);
    //~ end = start+data.x.length;
  //~ }

  /* Slicing is done using dates rather than observations, which wouldn't
   * really make a lot of sense. */
  //~ this(T)(T obj, long _start) {
    //~ data = RData("ts(" ~ obj.name ~ ", start=" ~ _start.to!int ~ ", frequency=1)");
    //~ ptr = REAL(data.x);
    //~ frequency = 1;
    //~ start = _start;
    //~ end = _start+obj.length-1;
  //~ }
  
  //~ this(T)(T obj, long[2] _start) {
    //~ data = RData("ts(" ~ obj.name ~ ", start=c(" ~ _start[0].to!string ~ ", " 
      //~ ~ _start[1].to!string ~ "), frequency=" ~ _frequency.to!string ~ ")");
    //~ ptr = REAL(data.x);
    //~ frequency = _frequency;
    //~ start = asLong(_start);
    //~ end = start+data.x.length;
  //~ }
  
  this(T)(T v, long[] s) if (__traits(hasMember, v, "data") && 
	(is(typeof(__traits(getMember, v, "data")) == RData))) {
    data = RData("ts(" ~ v.name ~ ", start=c(" ~ s[0].to!string ~ ", " 
      ~ s[1].to!string ~ "), frequency=" ~ frequency.to!string ~ ")");
    ptr = REAL(data.x);
    frequency = this.frequency;
    _start = s;
    _end = fromLong(this.start!long + data.x.length - 1);
  }
  
  static if(freq == 1) {
		double opIndex(long d) {
			return ptr[d - this._start];
		}
	}
  
  //~ double opIndex(long d0, long d1) {
    //~ enforce(frequency > 1, "Indexing with two values requires monthly or quarterly data");
    //~ return ptr[asLong([d0, d1]) - start];
  //~ }
  
  //~ /* Since it's a date, the end point is included */
  //~ TS opSlice(long s, long e) {
    //~ return TS("window(" ~ this.name ~ ", start=" ~ s.to!string ~ ", end=" ~ e.to!string ~ ")");
  //~ }
  
  //~ TS opSlice(long[] s, long[] e) {
    //~ enforce(frequency > 1, "Do not use an array to denote the year for annual data");
    //~ return TS("window(" ~ this.name ~ ", start=c(" ~ s[0].to!string ~ ", " 
      //~ ~ s[1].to!string ~ "), end=c(" ~ e[0].to!string ~ ", " 
      //~ ~ e[1].to!string ~ "))");
  //~ }
  
  //~ TS until(long e) {
    //~ return TS("window(" ~ this.name ~ ", end=" ~ e.to!string ~ ")");
  //~ }

  //~ TS until(long[2] e) {
    //~ return TS("window(" ~ this.name ~ ", end=c(" ~ e[0].to!string ~ ", " 
      //~ ~ e[1].to!string ~ "))");
  //~ }
  
  //~ TS starting(long s) {
    //~ return TS("window(" ~ this.name ~ ", start=" ~ s.to!string ~ ")");
  //~ }
  
  //~ TS starting(long[2] s) {
    //~ return TS("window(" ~ this.name ~ ", start=c(" ~ s[0].to!string ~ ", " 
      //~ ~ s[1].to!string ~ "))");
  //~ }
  
  TS lag(long k) {
    return TS!freq("lag(" ~ this.name ~ ", " ~ to!string(-k) ~ ")");
  }
  
  //~ TS lead(long k) {
    //~ return lag(-k);
  //~ }
  
  //~ TS diff(long k) {
    //~ return TS("diff(" ~ this.name ~ ", " ~ k.to!string ~ ")");
  //~ }
  
  //~ TS pct(long k) {
    //~ return TS("(diff(" ~ this.name ~ ", " ~ k.to!string ~ ")/lag(" ~ this.name ~ ", " ~ to!string(-k) ~ "))");
  //~ }
  
  void print(string msg="") {
    if (msg.length > 0) {
      writeln("--------\n", msg, "\n--------");
    }
    static if(freq == 1) {
      writeln("\nStart: ", start);
      writeln("End: ", end);
    } else {
      writeln("\nStart: ", _start);
      writeln("End: ", _end);
    }
    writeln("Frequency: ", frequency);
    writeln();
    printR(data.x);
  }
}


struct MTS {
	RData data;
  long start;
  long end;
  int frequency;
  double * ptr;
  alias data this;
  
  /* In the future, make this a TS[string] to make printing sensible */
    
  this(string code) {
    data = RData(code);
    ptr = REAL(data.x);
    frequency = INTEGER(evalR("as.integer(frequency(" ~ data.name ~ "))"))[0];
    auto tmp = IntVector("as.integer(start(" ~ data.name ~ "))");
    start = asLong([tmp[0], tmp[1]], frequency);
    tmp = IntVector("as.integer(end(" ~ data.name ~ "))");
    end = asLong([tmp[0], tmp[1]], frequency);
  }  
  
  long asLong(long[2] d) {
    if (frequency != 1) {
      return (d[0]-1900)*frequency+(d[1]-1);
    } else {
      return d[0];
    }
  }
  
  long[2] fromLong(long d) {
    return [(d/frequency)+1900, d%frequency+1];
  }
  
  long asLong(long[2] d, long f) {
    if (f != 1) {
      return (d[0]-1900)*f+(d[1]-1);
    } else {
      return d[0];
    }
  }
  
  long[2] fromLong(long d, long f) {
    return [(d/f)+1900, d%f+1];
  }
  
  void print(string msg="") {
    if (msg.length > 0) {
      writeln("--------\n", msg, "\n--------");
    }
    if (frequency == 1) {
      writeln("\nStart: ", start);
      writeln("End: ", end);
    } else {
      writeln("\nStart: ", fromLong(start));
      writeln("End: ", fromLong(end));
    }
    writeln("Frequency: ", frequency);
    writeln();
    printR(data.x);
  }
}
	

MTS tsCombine(int f)(TS!f[] series) {
	enforce(series.length > 1, "tsCombine requires multiple time series");
	MTS result;
	evalRQ(`..tmp <- ` ~ series[0].data.name);
	foreach(var; series[1..$]) {
		evalRQ(`..tmp <- cbind(..tmp, ` ~ var.data.name ~ `)`);
	}
	evalRQ(`..tmp <- na.omit(..tmp)`);
	evalRQ(`colnames(..tmp) <- paste0("V", 1:ncol(..tmp))`);
	return MTS(`..tmp`);
}
		













