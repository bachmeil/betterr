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

You can do binary operations with `RCMatrix`:

- `+` and `-` return a new RCMatrix.
- `*` is matrix multiplication.
- `/` is matrix division in the sense of Matlab (solves a least squares problem)

### double[] array()

Returns a slice to the data in the matrix. Be aware that this is a reference,
for speed, and changes to the returned slice will affect the original matrix.
Use dup if that's not what you want.

```
auto mm = m.array();
mm[0] = 2.2; // Changes m[0,0]
auto mm2 = m.array.dup(); // New allocation plus copy of the data
mm2[0] = 2.2; // m[0,0] unaffected
```

### val

Set all elements of the matrix:

```
m.val = -3.2; // All elements are set to -3.2
m.val = [1.1, 2.2, 3.3, 4.4]; // Fills by column, m needs to have 4 elements
```

### ptr

Pointer to the underlying data array. Should rarely be used.

### rows

Returns the number of rows.

### cols

Returns the number of columns.

### length

Returns the total number of elements.

### cov

Returns the covariance matrix for the columns.

### cor

Returns the correlation matrix for the columns.

### shape

Reshapes the matrix to the given dimensions.

### RCMatrix trimRows(long trimTop, long trimBottom)

New matrix after trimming the given numbers of rows from the top and bottom.

### RCMatrix rowMin()

New matrix holding the minimum values of each row.

### RCMatrix rowMax()

New matrix holding the maximum values of each row.

### RCMatrix colMin()

New matrix holding the minimum values of each column.

### RCMatrix colMax()

New matrix holding the maximum values of each column.

### RCMatrix rowMinIndex()

New matrix holding the index of the minimum values of each row.

### RCMatrix rowMaxIndex()

New matrix holding the index of the maximum values of each row.

### RCMatrix colMinIndex()

New matrix holding the index of the minimum values of each column.

### RCMatrix colMaxIndex()

New matrix holding the index of the maximum values of each column.

### double min()

Minimum value of the entire matrix.

### double max()

Maximum value of the entire matrix.

### double sum()

Sum of all elements in the matrix.

### RCMatrix pca(RCMatrix m, long p)

The first p principal components of m.

### RCMatrix selectRows(RCMatrix m, RCMatrix sel)

New matrix holding the rows specified in sel.

### RCMatrix selectColumns(RCMatrix m, RCMatrix sel)

New matrix holding the columns specified in sel.

### RCMatrix sort(RCMatrix m, long col)

Sort m by column col.

### colnames, rownames

You can get and set the column and row names.

### void print(RCMatrix m, string msg="")

Print the matrix m with optional initial message msg.

### RCMatrix inv(RCMatrix m)

Inverse of m.

### double det(RCMatrix m)

Determinant of square matrix m.

### double logdet(RCMatrix m)

Log determinant of m.

### double logabsdet(RCMatrix m)

Log of the absolute value of the determinant.

### RCMatrix solve(RCMatrix m1, RCMatrix m2)

Solve a system of linear equations.

### 
