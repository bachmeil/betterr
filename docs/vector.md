# Vector

A `Vector` is a struct that holds information about a vector that has been allocated inside R. The data is:

```
long rows;
RData data;
double * ptr;
```

You can use assignment and slicing as in languages like R and Matlab. Keep in mind that
slices do not include the last element. Here's an example:

```
// Create a new vector with space for five elements
auto x = Vector(5);

// Set all five elements to 1.7
x = 1.7;

// Print x with a message preceding it
x.print("Test vector");

// You can also use the empty slice operator to set all elements to the same scalar
x[] = -6.2;
x.print("Should be -6.2");

// You can use the empty slice operator to copy in the elements of a double[]
// There is a check that the lengths match
x[] = [1.1, 2.2, 3.3, 4.4, 5.5];
x.print("Should be increasing values");

// You can set a subset of the vector to the same scalar
x[2..5] = -1.2;
x.print("Last three elements should be -1.2");

// You can also set a subset of the vector to the elements of a double[]
// There's a check that the lengths match
x[2..5] = [-0.2, -0.3, -0.4];
x.print("Last three elements changed");
```

# Functionality

## Construction

### this(string code)

You can pass a string of R code and store the result as a Vector.

```
auto x = Vector(rnorm(10));
```

### this(long r)

You can send the number of elements and a Vector with that many elements will be allocated.
The elements are initialized to zero.

```
auto x = Vector(12)
```
  
### void initialize(long r)

This is used for assigning to a zero-length Vector. Otherwise there is
no way to assign to a Vector that is part of a struct that has already
been constructed, since the bounds checking will fail.

### this(Vector v)

Makes a copy in R, but forces the allocation of a new Vector.

```
auto y3 = Vector(x);
```

### this(RData rd)
  
Makes a copy in R, but forces the allocation of a new Vector.

### this(double[] v)

Allocates a new Vector and copies the elements of v into it.

```
auto x = Vector([1.1, 2.2, 3.3]);
```

### this(long[] v)

### this(int[] v)

## Indexing

```
// Indexing
v[1]

// Grab the first, third, and fifth elements
v[[0, 2, 4]]

// Set elements
v[1] = 3.3;
```

### double opIndex(long r)

Get one element. Indexing starts at zero.
  
### Vector opIndex(long[] obs)

Return a new vector holding the values at the indexes in `obs`.

### void opIndexAssign(double v, long r)

Assign to the element at index r.

## Assignment

### void opAssign(Vector v) 

Checks that the number of elements match. The exception is if there are
zero rows, which is taken to mean this has not been initialized. Will
allocate the vector and copy the elements of `v` into it.
  
### void opAssign(double[] v)

Checks that the number of elements match. The exception is if there are
zero rows, which is taken to mean this has not been initialized. Will
allocate the vector and copy the elements of `v` into it.
 
## Slicing

### Vector opSlice(long i, long j)

`this[0..3]` returns a new Vector with three elements. `j` is not included.
  
### Vector opSlice()

`this[]` returns a new Vector with all elements of `this`.

### void opSliceAssign(double a)

This does the same thing as opAssign.

### void opSliceAssign(double[] v)

This does the same thing as opAssign.

### void opSliceAssign(double a, long ind0, long ind1) 

### void opSliceAssign(double[] v, long ind0, long ind1)

## Range support

A `Vector` is a range. You can do this type of thing:

```
auto v = Vector([1.1, 2.2, 3.3]);

foreach(val; Vector([1.1, 2.2, 3.3])) {
	writeln(val);
}

import std.range: enumerate;
foreach(ii, val; Vector([1.1, 2.2, 3.3]).enumerate) {
	writeln(ii+1, " ", val);
}
```

## Miscellaneous functions

### double last()

The last value
  
### long opDollar()

Allows you to do `v[3..$]`

### void print(string msg="")

### long length()

### bool empty()

### double front()

### void popFront()

### Matrix matrix()

Returns a Matrix holding the data in this as one column.
	
###	Matrix rowMatrix()

Returns a Matrix holding the data in this as one row.

### double[] opCast(T: double[])()

Convert this to a `double[]` using std.conv.to:

```
auto v = Vector([1.1, 2.2, 3.3]);
double[] vv = v.to!(double[]);
```
  
### Vector opBinary(string op)(Vector v)

Allows element-by-element operations like this:

```
auto v = Vector([1.1, 2.2, 3.3]);
auto v2 = Vector([4.4, 5.5, 6.6]);
v += v2;
v -= v2;
v *= v2;
v /= v2;
```

### Vector head(Vector v, long n=6)

Calls R's `head` function. Returns a new Vector of length `n`. If `n` is negative, drops `n` elements from the end, so the length is `v.length - n`.

### Vector tail(Vector v, long n=6)

Calls R's `tail` function. Returns a new Vector of length `n`.  If `n` is negative, drops `n` elements from the front, so the length is `v.length - n`.

## struct Fill

Used to fill the elements of a Vector. The main case is where you want to refill a Vector inside a loop and you'd like to avoid allocation for performance reasons. Fill is a reference to the Vector.

```
auto vraw = Vector(10);
auto v = Fill(vraw);
v ~= 1.1;
v ~= 2.2;
v ~= 3.3;
writeln(v.full()); // false, because the capacity is 7
writeln(v.capacity()); // 7
v ~= [4.4, 5.5, 6.6, 7.7, 8.8, 9.9, 11.0];
writeln(v.full()); // true, because the capacity is 0
v ~= 12.1; // Error, because there's no place to put another element
vraw.print("The vector"); // Work with vraw for anything other than filling
v.reset(); // Start over without creating garbage and avoiding allocation
```

## struct PartialFill

Same as `Fill`, but in case you only want to work with part of a Vector rather than the whole thing. The main case is where you're only updating some of the elements inside a loop because you're either not updating all elements or you're doing the update in several places.

```
auto vraw = Vector(10);
auto v = PartialFill(vraw, 3, 7); // You only want to fill vraw[3]..vraw[6]
```

## struct IntVector

`Vector` holds double precision values. `IntVector` holds integer values.

## struct BoolVector

`Vector` holds double precision values. `BoolVector` holds true/false values. Under the covers, since we're working with C, it holds `int` values. The interface is designed to feel like you're working with bool values.

```
auto bv = BoolVector(2);
bv[0] = true;
bv[1] = false;
writeln(bv[0]); // true, even though ptr[0] is 1
```

## struct StringVector

Create and work with a vector of R strings. Working with the underlying data is messy (they're C strings, after all). You should work with a `StringVector` instead of working with the R data directly.

```
auto sv = StringVector(["Foo", "Bar", "Baz", "D is kind of cool!"]);
sv.print("A string vector");
writeln(sv[1]);
```

