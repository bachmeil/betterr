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

```
auto m = Matrix(4, 6);
```

Allocates a new (r x c) Matrix in R. It has a unique name that's stored in `data`. The pointer `ptr` points to the underlying data array of length rc. Note that `ptr` changes over time (even frequently, depending on what you're doing). It's rarely a good idea to store `ptr` anywhere else.

### this(Matrix m)

Allocates a new Matrix with the same dimensions as `m`. Makes a copy of the data in `m`.

### this(Submatrix sm)

Allocates a new Matrix with the same dimensions as `sm`. Makes a copy of the data in `sm`.

### this(Vector v)

Copies the elements of `v` into a newly allocated Matrix with dimensions (v.length x 1).

### this(Vector v, long r, long c)

Copies the elements of `v` into a newly allocated Matrix with dimensions (r x c). Fills by column, not row.

### dup

Allocates a new Matrix of the same dimension and copies the elements into it.

### this(RData rd)

Should rarely be used in user code. Creates a new Matrix and copies the data into it. If the Robj inside `rd` is a Matrix, the dimensions will be the same as that Matrix. If the Robj is a Vector, it will be a Matrix with one column.

### this(string code)

Creates a new Matrix and copies the output of evaluating `code` into it. If `code` evaluates to a Matrix, the dimensions will be the same as that Matrix. If the Robj is a Vector, it will be a Matrix with one column.

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
// Alternative syntax: auto sm = m.sub;
// Submatrix, avoids a copy
m2[0..2, 0..2] = sm[0..2, 0..2];
// New Matrix; calls the Matrix constructor
Matrix m3 = sm[0..2, 0..2]; 
```

### Shouldn't slicing return a reference?

One could make the argument that any time you take a slice, it should return a reference to the corresponding parts of the matrix, not a new matrix with those elements copied into it. That's consistent with D's array slicing, where `v[1..4]` is a reference to that part of the array, and if you want a copy, you have to use `dup`. Indeed, that was the initial design, but this fails spectacularly unless you're really careful.

Consider how slices being references can go wrong in vanilla D code:

```
import std;
void main() {
    auto z = [1.1, 2, 3];
    writeln(z.ptr);
    auto z2 = z[];
    // Same as z.ptr
    writeln(z2.ptr);
    
    z ~= 4;
    // Now they're different due to a reallocation
    writeln(z.ptr);
    writeln(z2.ptr);
}
```

If you're not careful, `z2` might not be pointing to what you think it's pointing to, and you might end up with a disastrous outcome. The good news is that you're probably not going to run into many problems writing vanilla D code with slices (at least I don't).

The same problem exists when slicing returns a reference to the Matrix, but on a bigger scale. It's quite common to do things like take a row or a column of a Matrix in numerical code, and operations such as modifying the elements of a Matrix will generally lead to a reallocation, much more so than with D's built-in arrays. Holding a copy of the pointer simply does not work because the probability of it becoming invalid is so high. It's more complicated than writing C.

I considered an alternative solution. Rather than storing a pointer to the underlying data array, I can store the name of the variable in the Submatrix. That adds considerable overhead. On *every* access, you have to request the Robj that goes with the name, and then you have to get the pointer to the underlying data array. It would be an understatement to say this is inefficient.

Something that might work is to have the Submatrix hold a pointer to the Matrix. Then on each access, grab the pointer to the data array. While I won't rule out doing this in the future, since the syntax would be convenient and it would be consistent with other slicing in D, I'm hesitant to add a second pointer. Someone wanting speed, which is really the only reason to put up with the inconvenience of a reference type, is unlikely to want an extra level of indirection.

For better or worse, the current design requires you to explicitly specify that you want a reference. That's the clearest for the reader of the code and delivers the best performance. Almost certainly something shorter than `reference` will be used. You can limit the use of references to only those cases where they're crucial for performance, and you can limit the set of opportunities to mess things up.

## Other Functions

These work as expected, so there's not much elaboration needed.

### Vector rowSums()

### Vector colSums()

### Vector rowMeans()

### Vector colMeans()

### Vector row(long ii)

Returns a newly allocated Vector with the elements of row ii copied into it.

### Vector column(long ii)

Returns a newly allocated Vector with the elements of column ii copied into it.

### Vector lastrow()

Returns a newly allocated Vector with the elements of the last row copied into it.

### Vector lastcolumn()

Returns a newly allocated Vector with the elements of the last column copied into it.

### Matrix matmul(Matrix x, Matrix y)

Returns a newly allocated Matrix holding the product of `x` and `y`. Note that this is matrix multiplication, not element-by-element multiplication. The equivalent of R's `%*%` operator.

### Matrix matmul(Vector v, Matrix y)

Converts `v` to a Matrix with one column, then does matrix multiplication.

### Matrix matmul(Matrix x, Vector v)

Converts `v` to a Matrix with one column, then does matrix multiplication.

### Matrix elmul(Matrix x, Matrix y)

Element-by-element multiplication of x and y. Explicit naming is used to avoid confusion.

### Matrix plus(Matrix x, Matrix y)

Returns a newly allocated matrix holding the sum of `x` and `y`.

### Matrix minus(Matrix x, Matrix y)

Returns a newly allocated matrix holding `x - y`.

### Matrix div(Matrix x, Matrix y) 

Returns a newly allocated matrix holding the element-by-element division `x / y`.

### Matrix mul(Matrix x, double a)

Returns a newly allocated matrix holding `ax`.

### Matrix plus(Matrix x, double a)

Returns a newly allocated matrix holding the result of adding `a` to every element of `x`.

### Matrix minus(Matrix x, double a)

Returns a newly allocated matrix holding the result of subtracting `a` from every element of `x`.

### Matrix div(Matrix x, double a)

Returns a newly allocated matrix holding the result of dividing every element of `x` by `a`.

### Matrix mul(double a, Matrix x)

Returns a newly allocated matrix holding `ax`.

### Matrix plus(double a, Matrix x)

Returns a newly allocated matrix holding the result of adding `a` to every element of `x`.

### Matrix minus(double a, Matrix x) 

Returns a newly allocated matrix holding the result of subtracting every element of `x` from `a`.

### Matrix div(double a, Matrix x) 

Returns a newly allocated matrix holding the result of dividing a by every element of `x`.

### Matrix t(Matrix x) 

Returns a newly allocated Matrix holding the transpose of `x`.

### Matrix solve(Matrix x) 

Following R, returns a newly allocated Matrix holding the inverse of `x`. See [the R documentation](https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/solve).

### Matrix solve(Matrix a, Matrix b) 

Solution of a system of equations. Solves `aX=b` for `X`. See [the R documentation](https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/solve).

### Vector solve(Matrix a, Vector b) 

Solution of a system of equations. Solves `aX=b` for `X`. See [the R documentation](https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/solve).

### Matrix inv(Matrix x)

Same as `solve(x)`.

### Vector diag(Matrix x) 

Returns a newly allocated Vector holding the elements of the diagonal of `x`. `x` is required to be square.

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






