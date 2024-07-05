# Compiling With ImportC

As of this writing, two of the components of BetterR are compiled using 
ImportC, the GSL random number generators and the matrix library from
Gretl. This is a big win, because we can use well-tested libraries that
have been used for years, allowing us to trust that most bugs will have
been found and various edge cases have been solved.

Since C files don't have access to D's module system, we're left with
two options: we can dump the C source and header files in the project directory,
or we can provide the compiler with the information it needs to find those
files. Note that this is true even if you have a D wrapper written over
the top of the C functionality. The D wrapper needs the location of the 
C source and header files.

This is what I add to the compilation command for the GSL RNG library:

```
-P-Igsl/rng gsl/rng/*.c
```

The `-P` passes info to the gcc preprocessor. `-Igsl/rng` tells it to
look in `gsl/rng` (relative directory in this case) for header files.
`gsl/rng/*.c` says to look in that directory for C source files. In the
case that I wanted to exclude some of the C files, I'd have to specify
them manually. If you're doing that, in my opinion, you should probably
find a better way to organize your files.

I prefer to make it so that the user of the library only imports D modules.
The downside of doing that is that the user may not understand why they
have to provide the location of C files to the compiler.
