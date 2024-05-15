# Gretl Matrix Library

BetterR includes most of the matrix library from [Gretl](https://gretl.sourceforge.net/). I stripped out only the necessary pieces and created a standalone library that's included as part of BetterR. This code has not been ported from C to D - it's the original C code, and it's compiled with ImportC. That matters, because you have stable, well-tested code that's been in use for many years. The matrix library makes its low-level calls to [OpenBLAS](https://github.com/OpenMathLib/OpenBLAS) for speed.

There's matrix functionality elsewhere in BetterR. The Gretl matrix library is intended for use as a C library, it has less functionality than that of R, and will typically be faster than equivalent matrix operations with R functions. The R functions also call fast low-level matrix algebra libraries, but they do additional work like checking the type of the matrix.

## RCMatrix

The fundamental data type is `RCMatrix`, where the RC stands for "reference counted". `RCMatrix` is an alias for `SafeRefCounted!RawGretlMatrix`, where `RawGretlMatrix` is a struct that wraps a pointer to `gretl_matrix` (the type used by Gretl for its matrix operations). [SafeRefCounted](https://dlang.org/phobos/std_typecons.html#SafeRefCounted) handles the allocation and freeing of data by calling the appropriate Gretl functions. You should not need to worry about memory management when using this library.

## Methods

You can call the following methods on a `RCMatrix`:

### this(long rows, long cols)

Allocate a new (3x8) matrix:

```
auto m = RCMatrix(3, 8);
```

### this(double[] v, long rows, long cols)

Allocate a new (3x8) matrix and copy the data from `v` into it:

```
auto v = new double[24];
v[] = 1.3;
auto m = RCMatrix(v, 3, 8);
```

### RCMatrix dup()

Allocate a new matrix with the same dimensions and copy the data into it.

```
auto m2 = m.dup();
```

### Indexing

You can do the usual matrix indexing: `m[3, 7]`. You can also use a single
index, as that is sometimes useful, even if it should generally be avoiding
for being error-prone: `m[12]`. The index moves over columns, so if `m`
is (4x6), then `m[12]` is the same as `m[0,3]`. The same operations work 
for assignment. For example, `m[3, 7] = 1.6`.

Keep in mind that D is a zero-indexed language, so if there are six columns,
the first column has index 0 and the last has index 5.

### opDollar

You can use the `$` the same as with a D array. `$` is the length of the
array.

### Binary operations

You can do binary operations with `RCMatrix`.









