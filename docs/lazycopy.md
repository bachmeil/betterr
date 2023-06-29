# R lazy copying

For performance reasons, R delays copying as long as possible. Consider this code:

```
x <- rnorm(10)
y <- x
```

There is no reason for R to copy the elements of `x` into `y` at this point. It is sufficient - and much faster - for `y` to hold a reference to `x`. If this is the only thing that's ever done with `y`

```
print(y[4])
print(y[7])
```

there is no reason to ever make a copy of `x`, because a reference is just fine until there's mutation of `x` or `y`. This is a problem if you're working with `y` from D. If you did the above, and then you changed the elements of `y` from inside a D function, you'd actually be changing the elements of both `x` and `y`, because there's nothing that triggers the copy.

The library handles all of this for you. If you're doing the calls to R yourself, though, you'll have to account for it. One way is to use braces when you do the copy. Here's an example for a matrix:

```
m <- matrix(1:4, ncol=2)
mm <- m[]
```

Or you can do an assignment after creating `mm`:

```
m <- matrix(1:4, ncol=2)
mm <- m
mm[1,1] <- mm[1,1]
```

As soon as you do the assignment, `mm` is no longer a reference. You can see an example of this in testing/memory.d.

## Additional reading

[Finding and setting variables](https://rstudio.github.io/r-manuals/r-exts/System-and-foreign-language-interfaces.html#finding-and-setting-variables)

[Named objects and copying](https://rstudio.github.io/r-manuals/r-exts/System-and-foreign-language-interfaces.html#named-objects-and-copying)
