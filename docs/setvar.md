# Passing data from D to R

In some cases, it is simplest to pass data from D to R in string format. For instance, this is about as efficient as you'll get for setting a scalar integer inside R:

```
evalRQ("x <- 3L");
```

The time it takes for R to evaluate that string is negligible. That wouldn't work so well, however, if you have a 200x200 matrix. Suppose you allocate a Matrix `m` inside D that is the output of some linear algebra operations, and now that you're done with efficiently carrying out those operations using an optimized BLAS, you want to create a variable named `mat` inside R so it can be used in a machine learning library.

All you have to do is this:

```
toR(m.data.x, "mat");
```

`m.data` is where information about the data is stored inside a Matrix struct. `m.data.x` is the Robj that points to the data.
