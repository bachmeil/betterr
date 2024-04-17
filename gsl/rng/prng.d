module prng;

import std.conv, std.stdio, std.utf;
import gslrng;

/*
This is a D port of the file MRG31k3p.java from project SSJ. All of the
code that was ported was originally written in Java.

The original code can be found at [https://github.com/umontreal-simul/ssj/blob/master/src/main/java/umontreal/ssj/rng/MRG31k3p.java](https://github.com/umontreal-simul/ssj/blob/master/src/main/java/umontreal/ssj/rng/MRG31k3p.java).

Where noted, some supporting functions were taken from other parts of the SSJ
project and ported to D as well.

Author: Lance Bachmeier
Date: May 1, 2017
This update: March 6, 2024 Remove some clutter such as the N(0,1) generation function

The license is Apache 2.0, which is the same as the original code.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.


------------------------------


The `Stream` struct holds the state of an individual
generator. The state is defined by the six integers `x11`, `x12`, `x13`,
`x21`, `x22`, and `x23`. The initial values come from the `StreamGen`
struct that creates the `Stream`. There are also a bunch of constant 
terms, which I've made into enums since they'll never change.

When creating a multidimensional array in Java, the row index comes 
first. `m[2][6]` refers to the 7th element of the third row. We can do 
the same thing in D.

The struct has a single method, `gen`, which returns the next U(0,1)
draw and updates the state.
 
This class is for use with D functions. Use GSLStream with C functions.
*/

class Stream {
  // This is the state that changes when calling the RNG
  // stream and substream are set when the stream is created
  int x11, x12, x13, x21, x22, x23;
  
  // These are all constants
	enum long serialVersionUID = 70510L;
	enum int M1 = 2147483647;    //2^31 - 1
  enum int M2 = 2147462579;    //2^31 - 21069
  enum int MASK12 = 511;       //2^9 - 1
  enum int MASK13 = 16777215;  //2^24 - 1
  enum int MASK2 = 65535;      //2^16 - 1
  enum int MULT2 = 21069;
  enum double NORM = 4.656612873077392578125e-10;
  
  enum int[][] A1p0 =
    [ [0, 4194304, 129],
      [1, 0, 0],
      [0, 1, 0] ];
  enum int[][] A2p0 =
     [ [32768, 0, 32769],
       [1, 0, 0],
       [0, 1, 0] ];

  enum int[][] A1p72 =
    [ [1516919229, 758510237, 499121365],
      [1884998244, 1516919229, 335398200],
      [601897748, 1884998244, 358115744] ];
  enum int[][] A2p72 =
    [ [1228857673, 1496414766, 954677935],
      [1133297478, 1407477216, 1496414766],
      [2002613992, 1639496704, 1407477216] ];

  enum int[][] A1p134 =
    [ [1702500920, 1849582496, 1656874625],
      [828554832, 1702500920, 1512419905],
      [1143731069, 828554832, 102237247] ];
  enum int[][] A2p134 =
    [ [796789021, 1464208080, 607337906],
      [1241679051, 1431130166, 1464208080],
      [1401213391, 1178684362, 1431130166] ];
      
  this(int[6] xs) {
		x11 = xs[0];
		x12 = xs[1];
		x13 = xs[2];
		x21 = xs[3];
		x22 = xs[4];
		x23 = xs[5];
	}

	double gen() {
    int y1, y2;
    
    //first component
    y1 = ((x12 & MASK12) << 22) + (x12 >>> 9) + ((x13 & MASK13) << 7) + (x13 >>> 24);
    if (y1 < 0 || y1 >= M1) {
			y1 -= M1;
		}
		y1 += x13;
		if (y1 < 0 || y1 >= M1) {
			y1 -= M1;
		}
		
		x13 = x12;
		x12 = x11;
		x11 = y1;
		
		//second component
		y1 = ((x21 & MASK2) << 15) + (MULT2 * (x21 >>> 16));
		if (y1 < 0 || y1 >= M2) {
			y1 -= M2;
		}
		y2 = ((x23 & MASK2) << 15) + (MULT2 * (x23 >>> 16));
		if (y2 < 0 || y2 >= M2) {
			y2 -= M2;
		}
		y2 += x23;
		if (y2 < 0 || y2 >= M2) {
			y2 -= M2;
		}
		y2 += y1;
		if (y2 < 0 || y2 >= M2) {
			y2 -= M2;
		}
		
		x23 = x22;
		x22 = x21;
		x21 = y2;
		
		//Must never return either 0 or 1
		if (x11 <= x21) {
			return (x11 - x21 + M1) * NORM;
		} else {
			return (x11 - x21) * NORM;
		}
	}
}

/*
This is used to pass data to GSL functions. If you are only using D,
stick with `Stream`, as that is passed by reference. You need to pass
by reference, because that allows the state of the generator to be
updated. I wouldn't use structs at all if not for the need to interact
with C libraries.
*/ 

struct GSLStream {
  // This is the state that changes when calling the RNG
  // stream and substream are set when the stream is created
  int x11, x12, x13, x21, x22, x23;
  
  // These are all constants
	enum long serialVersionUID = 70510L;
	enum int M1 = 2147483647;    //2^31 - 1
  enum int M2 = 2147462579;    //2^31 - 21069
  enum int MASK12 = 511;       //2^9 - 1
  enum int MASK13 = 16777215;  //2^24 - 1
  enum int MASK2 = 65535;      //2^16 - 1
  enum int MULT2 = 21069;
  enum double NORM = 4.656612873077392578125e-10;
  
  enum int[][] A1p0 =
    [ [0, 4194304, 129],
      [1, 0, 0],
      [0, 1, 0] ];
  enum int[][] A2p0 =
     [ [32768, 0, 32769],
       [1, 0, 0],
       [0, 1, 0] ];

  enum int[][] A1p72 =
    [ [1516919229, 758510237, 499121365],
      [1884998244, 1516919229, 335398200],
      [601897748, 1884998244, 358115744] ];
  enum int[][] A2p72 =
    [ [1228857673, 1496414766, 954677935],
      [1133297478, 1407477216, 1496414766],
      [2002613992, 1639496704, 1407477216] ];

  enum int[][] A1p134 =
    [ [1702500920, 1849582496, 1656874625],
      [828554832, 1702500920, 1512419905],
      [1143731069, 828554832, 102237247] ];
  enum int[][] A2p134 =
    [ [796789021, 1464208080, 607337906],
      [1241679051, 1431130166, 1464208080],
      [1401213391, 1178684362, 1431130166] ];
      
  this(int[6] xs) {
		x11 = xs[0];
		x12 = xs[1];
		x13 = xs[2];
		x21 = xs[3];
		x22 = xs[4];
		x23 = xs[5];
	}

	double gen() {
    int y1, y2;
    
    //first component
    y1 = ((x12 & MASK12) << 22) + (x12 >>> 9) + ((x13 & MASK13) << 7) + (x13 >>> 24);
    if (y1 < 0 || y1 >= M1) {
			y1 -= M1;
		}
		y1 += x13;
		if (y1 < 0 || y1 >= M1) {
			y1 -= M1;
		}
		
		x13 = x12;
		x12 = x11;
		x11 = y1;
		
		//second component
		y1 = ((x21 & MASK2) << 15) + (MULT2 * (x21 >>> 16));
		if (y1 < 0 || y1 >= M2) {
			y1 -= M2;
		}
		y2 = ((x23 & MASK2) << 15) + (MULT2 * (x23 >>> 16));
		if (y2 < 0 || y2 >= M2) {
			y2 -= M2;
		}
		y2 += x23;
		if (y2 < 0 || y2 >= M2) {
			y2 -= M2;
		}
		y2 += y1;
		if (y2 < 0 || y2 >= M2) {
			y2 -= M2;
		}
		
		x23 = x22;
		x22 = x21;
		x21 = y2;
		
		//Must never return either 0 or 1
		if (x11 <= x21) {
			return (x11 - x21 + M1) * NORM;
		} else {
			return (x11 - x21) * NORM;
		}
	}
}

/*
The `StreamGen` struct is used to create the `Stream` structs that are
used for random number generation. It holds all the state necessary to
generate a new parallel random number generation stream. The only methods
a user should ever call are `createStream` and `createGSLStream`, and 
the latter only when planning to call GSL.
*/

struct StreamGen {
	private int[6] curr_stream = [12345, 12345, 12345, 12345, 12345, 12345];
	private int[6] stream, substream;
	private int x11, x12, x13, x21, x22, x23;
  enum int[][] A1p134 =
    [ [1702500920, 1849582496, 1656874625],
      [828554832, 1702500920, 1512419905],
      [1143731069, 828554832, 102237247] ];
  enum int[][] A2p134 =
    [ [796789021, 1464208080, 607337906],
      [1241679051, 1431130166, 1464208080],
      [1401213391, 1178684362, 1431130166] ];
	enum int M1 = 2147483647;
  enum int M2 = 2147462579;
	
	private void resetStartStream() {
		substream[] = stream[];
		resetStartSubstream();
	}
	
	private void resetStartSubstream() {
    x11 = substream[0];
    x12 = substream[1];
    x13 = substream[2];
    x21 = substream[3];
    x22 = substream[4];
    x23 = substream[5];
	}

	Stream createStream() {
		stream[] = curr_stream[];
		resetStartStream();
		multMatVect(curr_stream[], A1p134, M1, A2p134, M2);
		return new Stream([x11, x12, x13, x21, x22, x23]);
	}
}

void multMatVect(int[] v, int[][] A, int m1, int[][] B, int m2) {
  int[3] vv;
  vv[] = v[0..3];
  matVecModM(A, vv, vv[], m1);
  v[0..3] = vv[];
  vv[] = v[3..6];
  matVecModM(B, vv, vv[], m2);
  v[3..6] = vv[];
}

/*
from ssj/ArithmeticMod.java
*/
void matVecModM(int[][] A, int[] s, int[] v, int m) {
	auto x = new int[v.length];
	foreach (ii; 0..v.length) {
		x[ii] = 0;
		foreach (jj; 0..s.length) {
			x[ii] = multModM(A[ii][jj], s[jj], x[ii], m);
		}
	}
	v[] = x[];
}

/*
from ssj/ArithmeticMod.java
*/
int multModM (int a, int s, int c, int m) {
  int r = to!int( (a.to!long * s + c) % m);
  return r < 0 ? r + m : r;
}

/*
from ssj/ArithmeticMod.java

This code was pretty confusing. The Java version contains
`if ((v -= a1 * m) < 0.0)`
It is hard to understand that. Clearly the value of `v` is changed.
But what if the `if` statement is false? Is `v` still changed?

Some playing around in Java suggests that `v` is permanently changed from
that point on. Therefore the return statement in the `else` block is also
affected.

It turns out that this function is not needed at the current time.
*/

double multModM (double a, double s, double c, double m) {	
	enum double two53 = 9007199254740992.0;
	enum double two17 = 131072.0;
	int a1;
	double v = a * s + c;
	if (v >= two53 || v <= -two53 ) {
		a1 = to!int(a/two17);
		a -= a1 * two17;
		v  = a1 * s;
		a1 = to!int(v/m);
		v -= a1 * m;
		v  = v * two17 + a * s + c;
	}
	a1 = to!int(v/m);
	double v2 = v - a1*m;
	if (v2 < 0.0) {
		return v2 + m;
	}	else {
		return v2;
	}
}

// End of the parallel generation functions

/*
Generate a draw from U(0,1)
*/
double genDoubleGSL(ref GSLStream s) {
	return s.gen();
}

double genDouble(Stream s) {
	return s.gen();
}

/*
Generate one int between `i` and `j`.
*/

int genIntegerGSL(ref GSLStream s, int i, int j) {
	double u = s.gen();
	if (u <= 0.0) {
		return i;
	} else if (u >= 1.0) {
		return j;
	} else {
		return i + to!int(u*(j-i+1.0));
	}
}

int genInteger(Stream s, int i, int j) {
	double u = s.gen();
	if (u <= 0.0) {
		return i;
	} else if (u >= 1.0) {
		return j;
	} else {
		return i + to!int(u*(j-i+1.0));
	}
}

/*
Generate one draw from N(mu, sigma^2).
*/

//~ double genNormalGSL(ref GSLStream s, double mu=0.0, double sigma=1.0) {
	//~ return mu + sigma * inverseF01(s.gen());
//~ }

//~ double genNormal(Stream s, double mu=0.0, double sigma=1.0) {
	//~ return mu + sigma * inverseF01(s.gen());
//~ }

/*
This is code for interacting with GSL. The `gsl_rng_type` struct is
passed to libgsl as part of a `gsl_rng` struct. The only way I've found
to make this work is to create the `RNG` struct and then use `alias this`
to pass a pointer to the `gsl_rng` struct to the GSL functions.
`alias ptr this` inside `RNG` allows a more natural calling mechanism
where you create an `RNG` and then pass it to GSL functions, without
having to worry about pointers to data inside `RNG`.
*/

extern(C) void prngSet(void * state, ulong seed) {}

extern(C) ulong prngGet(void * state) { 
	GSLStream * s = cast(GSLStream*) state;
	return (*s).genIntegerGSL(0, int.max);
}

extern(C) double prngGetDouble(void * state) {
	GSLStream * s = cast(GSLStream*) state;
	return (*s).gen();
}

struct GslRng {
	gsl_rng gr;
	gsl_rng * ptr;
	gsl_rng_type grt;
	GSLStream gs;
	
	alias ptr this;
	
	this(Stream s) {
		gs = GSLStream([s.x11, s.x12, s.x13, s.x21, s.x22, s.x23]);
		grt = gsl_rng_type("prng".toUTFz!(char*), ulong.max, 0, 
			GSLStream.sizeof, &prngSet, &prngGet, &prngGetDouble);
		gr.type = &grt;
		gr.state = cast(void*) &gs;
		ptr = &gr;
	}
	
	void print() {
		writeln([gs.x11, gs.x12, gs.x13, gs.x21, gs.x22, gs.x23]);
		writeln("ptr: ", ptr, " ", &gr);
		writeln("state: ", *(cast(GSLStream*) gr.state));
		writeln("gen: ", prngGetDouble(gr.state));
		writeln("function: ", &prngGetDouble, " ", gr.type.get_double);
		writeln("gen again: ", gr.type.get_double(gr.state));
	}
}
