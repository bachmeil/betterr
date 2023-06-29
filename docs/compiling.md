# Compilation

You need to make the betterr modules available for import, and you need to link to `libR.so` and `libRInside.so`.

What works easiest for me is to create a subdirectory named betterr and create a symbolic link in there to all the betterr modules. `libR.so` should by default be in a standard location, so the linker should be able to find it without a path using the usual gcc syntax `-lR`. In contrast, `libRInside.so` will almost certainly not be in a standard location, so you *will* need to provide the full path. This is what compilation looks like on my Ubuntu 22.04 machine:

```
ldmd2 -i program.d -L/usr/lib/R/library/RInside/lib/libRInside.so -L-lR
```

## Finding libRInside.so

If you're not sure where to find `libRInside.so`, you can run this command in the terminal:

```
R -s -e 'cat(paste0(find.package("RInside"), "/lib/libRInside.so\n"))'
```

You can copy and paste the output into the Makefile. Alternatively, you can run this command inside R:

```
find.package("RInside")
```

It's in the /lib subdirectory relative to that directory. On my system, the call to `find.package` returns

```
[1] "/usr/lib/R/library/RInside"
```

which translates into the linker directive above.

## Finding Matrix.so

If you're wanting to do low-level linear algebra operations for performance reasons (probably not the case), you'll need to link to Matrix.so. You can find its location by running this command in the terminal:

```
R -s -e 'cat(paste0(find.package("Matrix"), "/libs/Matrix.so\n"))'
```

