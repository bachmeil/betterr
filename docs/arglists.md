# Understanding Argument Lists in the R source

This page provides some additional details about argument lists, beyond
the things that are described in [this section of the R extensions manual](https://rstudio.github.io/r-manuals/r-exts/System-and-foreign-language-interfaces.html#evaluating-r-expressions-from-c).
Argument lists are passed from R to C when calling functions using
.External or .Internal, for instance. Consider the matrix function
inside R, which ends with the line

```
.Internal(matrix(data, nrow, ncol, byrow, dimnames, missing(nrow), missing(ncol)))
```

The arguments are passed to the C function do_matrix as a single *pairlist*.
You can in many cases do function calls entirely in C, without any R code,
by constructing an appropriate pairlist with the arguments. Keep the
following in mind before you go "optimizing" your code with pairlists
to eliminate R from the equation.

- This is unlikely in most cases to be more efficient than the much simpler approach of
sending a string of R code to libR using evalR or evalRQ. You have to
convert all the arguments into Robj structs, and that's potentially
less efficient than passing a string of code that's parsed by R. It's
always a lot more work and more error-prone to construct a pairlist.
- It will be the only way to call C functions that don't
have a direct way to be called from R. It may also be worthwhile if you're
calling a particular C function many times in a loop, but the R version of
the function adds significant overhead by doing expensive 
transformations and type checks you don't need.

Let's look at an example using the Mmatrix function in the Matrix package. It's used to
create a Matrix (more general than the matrix in base R). [Here's a link to the function source](https://github.com/cran/Matrix/blob/6cf48cd5e2a14d81580201b3e0bf923b8f782705/src/Mutils.c#L1226).

The function signature is `SEXP Mmatrix(SEXP args)`. That's not particularly
helpful in the absence of further information about the contents of `args`.
Much of the time, `args` in the source refers to the arguments to an
.External or .Internal call, and the name of the function goes first.
The documentation for this function is written out in the comment above it.

```
External(Mmatrix,
//	     data, nrow, ncol, byrow, dimnames,
//	     missing(nrow), missing(ncol))
```

First, we need to figure out what values to send for each of the
eight arguments. Then we need to figure out how to construct a pairlist
of length 8 that 23 can pass to `Mmatrix`.

## Determining the argument values

As noted, we need to supply eight items in the pairlist of arguments to
Mmatrix. Here are the items, all of which have to be converted to Robj
structs, rather than native C types like int and double:

1\. The symbol name for the C function (Mmatrix). It is *not* a string,
it's an Robj struct of type SYMSXP. The easy way to create that Robj
is to use the function `Rf_install`, which converts a C string into an
appropriate Robj. To do this from D, I will use

```
Rf_install(toUTFz!(char*)("Matrix"))
```

2\. The data used to initialize the Matrix. That not only provides the
initial values (if any) to use on construction of the Matrix, it also
identifies the type of Matrix you're creating. You can use any numeric vector or
matrix for the data, but we'll create a double precision vector
with 0 elements, which will identify the type as double, but initialize
all elements to NA. That is accomplished with `Rf_allocVector(14,0)`.

3\. and 4. The row and column dimensions. We'll create
a (3x2) Matrix, so these are `Rf_ScalarInteger(3)` and `Rf_ScalarInteger(2)`,
respectively, or their shorter aliases `3.robj` and `2.robj`.

5\. A flag that's true if the Matrix should be filled by row rather
than by column. Since we're filling by column, the flag will be false: `RFalse`.

6\. The dimension names. That's somewhat messy to do from D, so let's leave it
empty: `RNil`. It's easier - and more efficient because there's less
copying of strings - to work with dimension names from R.

7\. and 8. Flags telling if the number of rows and number of columns,
respectively, have not been specified (and thus have to be inferred).
Since we've specified them already, these arguments will both be `RFalse`.

## Putting them into a pairlist

Now that we know the representation of all eight arguments, we need to
put them into a single pairlist. This step is a bit interesting (and possibly 
even confusing). What we'll be doing is known in the Lisp world as
constructing [a linked list](https://gigamonkeys.com/book/they-called-it-lisp-for-a-reason-list-processing.html).
In practice, we have a bunch of pairs of pointers. The first is a pointer
to an Robj that holds that element of the data, and the second is a
pointer to either another pair (if it's not the last element) or
NULL (if it's the last element of the list). As a side note, one of the 
creators of the linked list was Herbert Simon, a Nobel Prize recipient in 
economics who won a Turing Award in his spare time.

The first step is to create a vector holding the pairlist. Since we need
eight elements:

```
Robj argListPointer = Rf_protect(Rf_allocList(8));
```

You can think about the layout of argListPointer like this, where each
line represents one element:

```
[ptr0, ptr1]
[ptr1, ptr2]
[ptr2, ptr3]
[ptr3, ptr4]
[ptr4, ptr5]
[ptr5, ptr6]
[ptr6, ptr7]
[ptr7, NULL]
```

argListPointer holds a pointer to the first pair, `[ptr0, ptr1]`. The
first pointer of a pairlist is called the CAR and the second is called
the CDR. These names are borrowed from Lisp. Don't worry about their
meaning because [they are basically nonsense that wasn't given much thought way back in 1959](http://www.iwriteiam.nl/HaCAR_CDR.html).

```
SET_TYPEOF(argListPointer, 6);
```

This specifies that argListPointer points to a special type of pairlist 
that holds language objects. To traverse and fill each element of the
pairlist, we'll need a temporary Robj that we can change. It will start
with the first pair, so

```
Robj fillPointer = argListPointer;
```

Now we can go through and fill in the items outlined in the previous
section. To set the value at ptr0, we need to set the CAR equal to the
symbol name, which we do using the SETCAR function:

```
SETCAR(fillPointer, Rf_install(toUTFz!(char*)("Matrix")));
```

Now we need to set the second element to hold the data argument. As you
can see above, the pointer to that element is ptr1. The
trick is to set fillPointer equal to the second pointer in the first
element, which is also ptr1. (Hence the name "linked list".) This code 
does what we need:

```
fillPointer = CDR(fillPointer);
SETCAR(fillPointer, Rf_allocVector(14,0));
```

The other elements can be filled accordingly:

```
fillPointer = CDR(fillPointer);
SETCAR(fillPointer, 3.robj);
fillPointer = CDR(fillPointer);
SETCAR(fillPointer, 2.robj);
fillPointer = CDR(fillPointer);
SETCAR(fillPointer, RFalse);
fillPointer = CDR(fillPointer);
SETCAR(fillPointer, RNil);
fillPointer = CDR(fillPointer);
SETCAR(fillPointer, RFalse);
fillPointer = CDR(fillPointer);
SETCAR(fillPointer, RFalse);
```

Now that we've done that, all elements of argListPointer have been filled.
You can print it out to confirm:

```
printR(argListPointer);
```

The last step is to send argListPointer to Mmatrix and verify that it
creates the Matrix we're after:

```
printR(Mmatrix(argListPointer));
```

In practice, you'd want to (i) save the output of the call to Mmatrix as an
Robj, (ii) save it as a reference counted struct, or (iii) give it a name
and put it into the R global environment. Examples of all three can be found
in testing/testpair.d.
