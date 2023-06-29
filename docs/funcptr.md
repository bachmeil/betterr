# Passing function pointers around

It is possible to pass data pointers and function pointers to R. Since all
data in R has to be in the form of an Robj, you need to use functions to
do the conversions from pointer to Robj and vice versa. The primary
reason you'd want to do this is if you need to pass a function to be
evaluated to R. 

For instance, consider R's optimization routines. If you
want to use the algorithms underlying `optim`, you can access the C
functions directly. If, on the other hand, you want to call `constrOptim`,
there's no way to do that without writing the objective function and
gradient function in R. This is precisely the case where you want to use
a compiled language. The evaluation of the objective function is potentially
expensive, and it may be called many times as it iterates to convergence,
and many solutions may be needed (think about bootstrapping).

The usual approach to compiled functions in R is to create a shared
library, load the library, and call it from R using `.Call`. You can
avoid the need to create a shared library by passing pointers to the
objective function and the data into R, then having a C function that
uses those pointers to call the function using the data as an argument.

The package [funcptr](https://github.com/bachmeil/r-funcptr) handles
the messy part for you. That package can be installed using devtools. 
You can see examples of its usage in testing/testfunptr.d.

Note that this is not something I've done often, so while the example
works, there's not a lot of documentation. You're on your own if you
don't understand how C function pointers and casting them works.
