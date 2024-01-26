module betterr.array;
import betterr.rdata;
import betterr.r;
import betterr.matrix, betterr.vector;
import std.conv, std.exception, std.range, std.stdio;
import std.algorithm.comparison: max, min;

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
		long longStart() {
			return _start;
		}
		long longEnd() {
			return _end;
		}
	} else {
		private TimePeriod _start;
		private TimePeriod _end;
		long[] start() {
			return _start.array;
		}
		long[] end() {
			return _end.array;
		}
		void start(long[] s) {
			_start.year = s[0];
			_start.minor = s[1];
		}
		void end(long[] e) {
			_end.year = e[0];
			_end.minor = e[1];
		}
		// Shouldn't normally use these, but they're here if you want them
		long longStart() {
			return _start.to!long;
		}
		long longEnd() {
			return _end.to!long;
		}
	}

  //~ long asLong(long[] d) {
    //~ if (frequency != 1) {
      //~ return (d[0]-1900)*frequency+(d[1]-1);
    //~ } else {
      //~ return d[0];
    //~ }
  //~ }
  
  //~ long[] fromLong(long d) {
    //~ return [(d/frequency)+1900, d%frequency+1];
  //~ }
  
  //~ long asLong(long[] d, long f) {
    //~ if (f != 1) {
      //~ return (d[0]-1900)*f+(d[1]-1);
    //~ } else {
      //~ return d[0];
    //~ }
  //~ }
  
  //~ long[] fromLong(long d, long f) {
    //~ return [(d/f)+1900, d%f+1];
  //~ }
  
  /* To access an existing ts object that's already inside R */
  this(string code) {
    data = RData(code);
    ptr = REAL(data.x);
    _start = TimePeriod(freq);
    _end = TimePeriod(freq);
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
  
  static if(freq == 1) {
		this(string code, long s) {
			data = RData("ts(" ~ code ~ ", start=" ~ _start.to!string ~ ")");
			ptr = REAL(data.x);
			_start = s;
			_end = s + data.x.length-1;
		}
    
		this(T)(T v, long s) 
		if (__traits(hasMember, v, "data") && 
		(is(typeof(__traits(getMember, v, "data")) == RData))) {
			this(v.name, s);
		}
	} else {
		this(string code, long[] s) {
			data = RData("ts(" ~ code ~ ", start=c(" ~ s[0].to!string ~ ", " 
				~ s[1].to!string ~ "), frequency=" ~ frequency.to!string ~ ")");
			ptr = REAL(data.x);
			_start = TimePeriod(freq);
			_end = TimePeriod(freq);
			_start = s;
			_end = _start + (data.x.length-1);
		}

		/* Anything with an RData member named data */
		this(T)(T v, long[] s) 
		if (__traits(hasMember, v, "data") && 
		(is(typeof(__traits(getMember, v, "data")) == RData))) {
			//~ data = RData("ts(" ~ v.name ~ ", start=c(" ~ s[0].to!string ~ ", " 
				//~ ~ s[1].to!string ~ "), frequency=" ~ frequency.to!string ~ ")");
			//~ ptr = REAL(data.x);
			//~ _start = TimePeriod(freq);
			//~ _end = TimePeriod(freq);
			//~ _start = s;
			//~ _end = _start + (data.x.length-1);
			this(v.name, s);
		}
	}

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
  
  bool notBefore(long[] x, long[] y) {
    if (x[0] > y[0]) {
      return true;
    } else if ((x[0] == y[0]) && (x[1] >= y[1])) {
      return true;
    } else {
      return false;
    }
  }
  
  bool notBefore(long[] x, TimePeriod y) {
    return notBefore(x, y.array);
  }
  
  bool notBefore(TimePeriod x, long[] y) {
    return notBefore(x.array, y);
  }
  
  bool notBefore(TimePeriod x, TimePeriod y) {
    return notBefore(x.array, y.array);
  }
  
  bool notAfter(long[] x, long[] y) {
    if (x[0] < y[0]) {
      return true;
    } else if ((x[0] == y[0]) && (x[1] <= y[1])) {
      return true;
    } else {
      return false;
    }
  }
  
  bool notAfter(long[] x, TimePeriod y) {
    return notAfter(x, y.array);
  }
  
  bool notAfter(TimePeriod x, long[] y) {
    return notAfter(x.array, y);
  }
  
  bool notAfter(TimePeriod x, TimePeriod y) {
    return notAfter(x.array, y.array);
  }
  
  static if(freq == 1) {
		double opIndex(long d) {
			return ptr[d - this._start];
		}

    /* Since it's a date, the end point is included */
    TS opSlice(long s, long e) {
      return TS("window(" ~ this.name ~ ", start=" ~ s.to!string ~ ", end=" ~ e.to!string ~ ")");
    }

    TS until(long e) {
      return TS("window(" ~ this.name ~ ", end=" ~ e.to!string ~ ")");
    }

    TS starting(long s) {
      return TS("window(" ~ this.name ~ ", start=" ~ s.to!string ~ ")");
    }
	} else {
    /* If you want a single date returned as a TS, use the slice operator
     * with the same value for start and end. */
    double opIndex(long y, long m) {
      enforce(notBefore([y, m], this.start), "Index precedes the start of the TS");
      enforce(notAfter([y, m], this.end), "Index is after the end of the TS");
      return ptr[ [y, m] - this._start ];
    }
    
    /* Since it's a date, the end point is included */
    TS opSlice(long[] s, long[] e) {
      enforce(notAfter(s, e), "Start date cannot be after the end date");
      enforce(notBefore(s, this.start), "Start date prior to start of series");
      enforce(notAfter(e, this.end), "End date after start of series");
      return TS("window(" ~ this.name ~ ", start=c(" ~ s[0].to!string ~ ", " 
        ~ s[1].to!string ~ "), end=c(" ~ e[0].to!string ~ ", " 
        ~ e[1].to!string ~ "))");
    }

    TS until(long[] e) {
      return opSlice(this.start, e);
      //~ return TS("window(" ~ this.name ~ ", end=c(" ~ e[0].to!string ~ ", " 
        //~ ~ e[1].to!string ~ "))");
    }
    
    TS starting(long[] s) {
      return opSlice(s, this.end);
      //~ return TS("window(" ~ this.name ~ ", start=c(" ~ s[0].to!string ~ ", " 
        //~ ~ s[1].to!string ~ "))");
    }
  }
  
  TS lag(long k) {
    return TS!freq("lag(" ~ this.name ~ ", " ~ to!string(-k) ~ ")");
  }
  
  TS lead(long k) {
    return lag(-k);
  }
  
  TS diff(long k) {
    return TS("diff(" ~ this.name ~ ", " ~ k.to!string ~ ")");
  }
  
  TS pctChange(long k) {
    return TS("(diff(" ~ this.name ~ ", " ~ k.to!string ~ ")/lag(" ~ this.name ~ ", " ~ to!string(-k) ~ "))");
  }
  
  double[] array() {
		return ptr[0..(_end-_start+1)];
	}
  
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

struct TimePeriod {
	long year;
	long minor;
	long frequency;
	
	this(long f) {
		frequency = f;
	}
	
	this(long y, long m, long f) {
		year = y;
		minor = m;
		frequency = f;
	}
	
	TimePeriod opBinary(string op: "+")(long k) {
		enforce(k >= 0, "Can only add a positive number to a TimePeriod");
		TimePeriod result;
		long tmp = (minor-1) + k;
		result.year = this.year + (tmp / frequency);
		result.minor = (tmp % frequency) + 1;
		result.frequency = frequency;
		return result;
	}
	
	TimePeriod opBinary(string op: "-")(long k) {
		enforce(k >= 0, "Can only subtract a positive number from a TimePeriod");
		TimePeriod result;
		result.year = this.year - k/frequency;
		result.minor = this.minor - k % frequency;
		if (result.minor < 1) {
			result.minor += frequency;
			result.year -= 1;
		}
		result.frequency = this.frequency;
		return result;
	}
	
	long opBinary(string op: "-")(TimePeriod j) {
		return this.to!long - j.to!long;
	}
	
  long asLong(long[] d) {
    return frequency*(d[0]-1900) + (d[1]-1);
  }

  long opBinaryRight(string op: "-")(long[] d) {
    return asLong(d) - asLong([year, minor]);
  }
  
	long[] array() {
		return [year, minor];
	}

	long opCast(T: long)() {
		return frequency*(year - 1900) + (minor-1);
	}
	
	void opAssign(long[] d) {
		year = d[0];
		minor = d[1];
	}
	
	void opAssign(int[] d) {
		year = d[0];
		minor = d[1];
	}
	
	void opAssign(long d) {
		year = (d/frequency) + 1900;
		minor = (d%frequency) + 1;
	}
}

struct MTS {
	RData data;
  long start;
  long end;
  long frequency;
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
  
  this(long ncol, long s, long e) {
		writeln("e: ", e);
		writeln("s: ", s);
		string code = "ts(matrix(0.0, nrow=" ~ to!string(e-s+1) ~ ", ncol=" ~ to!string(ncol) ~ "), start=" ~ to!string(s) ~ ")";
		data = RData(code);
		ptr = REAL(data.x);
		start = s;
		end = s;
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
	

MTS tsCombine(long f)(TS!f[] series) {
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
		
struct MTSTransform {
	/* This struct will be used to hold info on transformations to make
	 * to the elements of a MTS struct. The purpose is mainly efficiency,
	 * because you need only to specify the dataset, and then one matrix
	 * is allocated, and it's filled one time, using the information about
	 * the transformation. */
	 TSTransform[] data;
	 long start;
	 long end;
	 long nrow;
	 
	 MTS create() {
		 start = long.min;
		 end = long.max;
		 foreach(var; data) {
			 start = max(start, var.modStart());
			 end = min(end, var.modEnd());
		 }
		 auto result = MTS(data.length, start, end);
		 nrow = end - start + 1;
		 foreach(ii, var; data) {
			 auto tmp = result.ptr[ii*nrow..(ii+1)*nrow];
			 var.compute(tmp, start, end);
		 }
		 return result;
	 }
}

struct TSTransform {
	void delegate(ref double[], long, long) compute;
	long delegate() modStart;
	long delegate() modEnd;
}

TSTransform Lag(long f)(TS!f var, long k=1) {
	double[] source = var.array;
	long arrayStart = var.longStart;
	long arrayEnd = var.longEnd;
	
	/* Guaranteed all elements from s to e can be computed */
	void compute(ref double[] target, long s, long e) {
		/* s-k is the date of the lag of the first observation after
		 * transformation.
		 * Then we subtract arrayStart to get the index inside source.
		 * Do the same through e+1 (since e+1 is not included). */
		target[0..$] = source[(s-k-arrayStart)..(e+1-k-arrayStart)];
	}
	
	/* First non-missing observation available for the transformed series */
	long modStart() {
		return arrayStart+k;
	}
	
	/* Last non-missing observation available for the transformed series */
	long modEnd() {
		return arrayEnd+k;
	}
	
	return TSTransform(&compute, &modStart, &modEnd);
}






