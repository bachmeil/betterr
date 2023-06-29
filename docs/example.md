# An Example

Here's a simple example to get you started. Save the following in a file called example.d.

```
import betterr.r, betterr.vector;
import std.stdio;

void main() {
  startR();
  auto v = Vector([1.1, 2.2, 3.3]);
  writeln(v);
  foreach(ii; 0..3) {
    v[ii] = v[ii]*2.5;
  }
  writeln(v);
  closeR();
}
```

Compile as described on the [compiling](compiling.html) page and run:

```
ldmd2 -i example.d -L/usr/lib/R/library/RInside/lib/libRInside.so -L-lR
./example
```

## Notes

1. Always start your program with `startR()` and end it with `closeR()`. If you don't do that, you're likely to run into segfaults for no obvious reason. The first thing you should check if getting segfaults is that you've called `startR`.
2. `betterr.r` provides basic interoperability features. `startR` and `closeR` are located in that module.
3. `betterr.vector` provides the `Vector` data structure.
4. `writeln` uses R to print the value of an R object, so the printed output will be in a familiar format if you've worked with R before.
