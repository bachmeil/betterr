# Matrix

A `Matrix` is a struct that holds information about a matrix that has been allocated inside R. The data is:

```
long rows;
long cols;
RData data;
double * ptr;
```

where `rows` and `cols` are what you'd expect, `data` is the reference to the data allocated by R, and `ptr` is a pointer to the underlying data array. Allocation of a new matrix is done by R. `ptr` allows access and setting of elements to be done without involving R, which in many cases is the bottleneck when using R for simulations or other situations that require working with individual elements. It also allows speed gains for reasons that may be surprising. Consider this code in R:

```
m[1:3, 2:4] <- m2
```

The underlying data structure might change when you assign to a `Matrix`, due to R's delayed copying. Imagine doing that in the middle of a loop that is executed millions of times. The equivalent code in D looks like this:

```
m[1..3, 2..4] = m2;
```

The dimensions of `m2` can be confirmed to match `m[1..3, 2..4]` and the corresponding elements of `m2` can be copied into the right place in `m` directly by using `ptr`. There is never a need to reallocate `m` and copy its elements for this operation.

# Functionality

## Construction

### this(long r, long c)

### this(Matrix m)

### this(Submatrix sm)

### this(Vector v)

Copies the elements of `v` into a newly allocated matrix with dimensions (v.length x 1).

### this(Vector v, long r, long c)

Copies the elements of `v` into a newly allocated matrix with dimensions (r x c). Fills by column rather than row.

### dup

Allocates a new matrix and copies the elements into it.

### this(RData rd)

Should rarely be used in user code

### this(string code)

Should rarely be used in user code

### Examples

```
auto m1 = Matrix(10, 25);
auto m2 = Matrix(m1); // R creates a copy of m1
auto m1ref = m1.reference(); // Creates a reference to avoid copying
auto m3 = m1ref[0..2, 0..2]; // The right side does not create a new matrix
auto v = Vector([1.1, 2.2, 3.3, 4.4]); // Use assignment to avoid the creation of a vector
auto m4 = Matrix(v); // Matrix with one column and four rows
auto m5 = Matrix(v, 2, 2); // Matrix with two columns and two rows
auto m6 = m5.dup;
```

## Indexing

### Overview

Multidimensional slicing allows the use of standard matrix notation.

```
m[1,4]
m[1..3, 4]
m[1..3, 4..8]
m[3, 0..$] // $ has the usual meaning; takes all elements of row 3
```

There is also a special struct type used to take all elements of a row or column:

```
m[3, _all] // All elements of the third row
m[_all, 3] // All elements of the third column
m[_all, _all] // The entire matrix
```

Finally, you might want to pull out a block with non-consecutive rows or columns. It is the same as R:

```
m[3, [1, 4, 9]] // Column 1, 4, and 9 of row 3
```

The above has the R equivalent `m[4, c(2, 5, 10)]`.

### Operations that avoid allocating a new Matrix

Accessing a single element returns a double:

```
double x = m[4,7];
```

Note that R is not involved in this operation. The above is equivalent to something like this:

```
double x = m.ptr[39];
```

Any operations on a Submatrix avoid allocating a new matrix (see the next section for details). That is done by calling `reference`:

```
Submatrix sm = m.reference();
```

### Operations that allocate a new matrix

Anything that returns more than a single element will allocate a new Matrix. In other words, this returns a new (2 x 2) matrix:

```
m[0..2, 0..2];
```

Suppose you're doing this operation:

```
m1 = m[0..2, 0..2];
```

The right side will allocate a new Matrix, the elements of that block of `m` will be copied into the new matrix, and then those elements will be copied into `m1`. Since the right side will never be used again, the new matrix will be destroyed. There is no easy solution to this. In this example

```
auto m2 = m[0..2, 0..2];
```

you want the right side to return a new Matrix. In principle, that could be solved by instead writing

```
Matrix m2 = m[0..2, 0..2];
```

but then we'd have to ban the use of `auto`. The solution I have adopted is to create a Submatrix, which is a reference to the underlying Matrix.

```
auto sm = m.reference;
m2[0..2, 0..2] = sm[0..2, 0..2]; // Obviously a Submatrix, avoids a copy
Matrix m3 = sm[0..2, 0..2]; // A new Matrix
```

I've tried it the other way around, where indexing a block of a Matrix always returns a Submatrix, but it's too complicated that way. An unnecessary allocation, as much as it might slow your code, is better than the program dying with a segmentation fault. I don't write `reference` very often in my own code. It's a simple optimization in return for never having to worry about segfaults and incorrect results.

## Other Functions

These generally work as expected, so there's not much elaboration here.

### Vector rowSums()

### Vector colSums()

### Vector rowMeans()

### Vector colMeans()

### Vector row(long ii)

Returns a newly allocated vector with the elements of row ii copied into it.

### Vector column(long ii)

Returns a newly allocated vector with the elements of column ii copied into it.

### Vector lastrow()

### Vector lastcolumn()

### Matrix matmul(Matrix x, Matrix y)

Matrix multiplication

### Matrix matmul(Vector v, Matrix y)

Convert the Vector to a Matrix, then apply matrix multiplication

### Matrix matmul(Matrix x, Vector v)

Convert the Vector to a Matrix, then apply matrix multiplication

### Matrix elmul(Matrix x, Matrix y)

Element-by-element multiplication of x and y. Explicit naming is used to avoid confusion.

### Matrix plus(Matrix x, Matrix y)

### Matrix minus(Matrix x, Matrix y)

### Matrix div(Matrix x, Matrix y) 

### Matrix mul(Matrix x, double y)

### Matrix plus(Matrix x, double y) 

### Matrix minus(Matrix x, double y)

### Matrix div(Matrix x, double y)

### Matrix mul(double y, Matrix x)

### Matrix plus(double y, Matrix x)

### Matrix minus(double y, Matrix x) 

### Matrix div(double y, Matrix x) 

### Matrix t(Matrix x) 

Transpose of x

### Matrix solve(Matrix x) 

Following R, returns the inverse of x. See [the R documentation](https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/solve).

### Matrix solve(Matrix x, Matrix y) 

Solution of a system of equations. See [the R documentation](https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/solve).

### Vector solve(Matrix x, Vector y) 

Solution of a system of equations. See [the R documentation](https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/solve).

### Matrix inv(Matrix x)

Inverse of x. Same as `solve(x)`.

### Vector diag(Matrix x) 

Returns the diagonal of x as a Vector. Requires x to be square.

### Matrix kronecker(Matrix x, Matrix y) 

### Matrix crossprod(Matrix x, Matrix y) 

### Matrix tcrossprod(Matrix x, Matrix y)

### Matrix crossprod(Matrix x)

### Matrix tcrossprod(Matrix x) 

### double det(Matrix x)

### Matrix diag(long ii)

Returns an (ii x ii) identity matrix.

### Matrix eye(long ii)

Returns an (ii x ii) identity matrix.

### Matrix cbind(Matrix m, Vector v)

Add v as a new column to m. A new Matrix is allocated. m is not affected.

### Matrix rbind(Matrix m, Vector v)

Add v as a new row to m. A new Matrix is allocated. m is not affected.

### Matrix cbind(Matrix m, double[] arr)

Add arr as a new column to m. A new Matrix is allocated. m is not affected.

### Matrix rbind(Matrix m, double[] arr)

Add arr as a new row to m. A new Matrix is allocated. m is not affected.

## Other structs

### struct SVD

A struct holding the results of a call to `svd`. See [the details here](https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/svd).

```
auto s = SVD(m);
Vector d = s.d; // See the R documentation for the interpretation of these members
Matrix u = s.u;
Matrix v = s.v;
```

### struct Eigen

A struct holding the results of a call to `eigen`. See [the details here](https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/eigen).

```
auto e = Eigen(m);
Vector values = e.values; // See the R documentation for the interpretation of these members
Matrix vectors = m.vectors;
```

### Column

Creates a reference to one column of a matrix, for convenience. You can get and set individual elements. Is a range, so you can use it with foreach.

### Row

Creates a reference to one row of a matrix, for convenience. You can get and set individual elements. Is a range, so you can use it with foreach.

### ColumnFill

Used to fill all the elements of one column of a matrix safely. See the discussion in [Fill for Vector](vector.html#struct-fill).

### RowFill

Used to fill all the elements of one column of a matrix safely. See the discussion in [Fill for Vector](vector.html#struct-fill).

### DiagonalFill

Used to fill all the elements of the diagonal of a matrix safely. See the discussion in [Fill for Vector](vector.html#struct-fill).

### AboveDiagonalFill

Used to fill all the elements above the diagonal of a matrix safely, filling by column. See the discussion in [Fill for Vector](vector.html#struct-fill).

### BelowDiagonalFill

Used to fill all the elements below the diagonal of a matrix safely, filling by column. See the discussion in [Fill for Vector](vector.html#struct-fill).

### MatrixFill

Used to fill all the elements of a matrix safely, filling by column. See the discussion in [Fill for Vector](vector.html#struct-fill).

### BlockFill

Used to fill all the elements of a block of a matrix safely, filling by column. See the discussion in [Fill for Vector](vector.html#struct-fill).

## IntMatrix

Holds a matrx of int values.

## BoolMatrix

Holds a matrix of bool values.






