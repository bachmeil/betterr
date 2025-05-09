# Notes for developers

These notes are useful if you want to work on the betterr library itself,
or if you want to understand how it works under the covers.

## rdata

This is the module to start with. Data is (with one exception) allocated
and managed by R.

```
struct RDataStorage {
	string name;
  Robj x;
  int refcount;
}
```

The elements of that struct are

- `name`: The name inside R that refers to the data. You'll use this if
you evaluate a string of R code that operates on a variable. For instance,
to run a regression, the `lm` module evaluates the `lm` function inside R
with the first argument being the R name of the variable.
- `x`: This is a pointer to the data. It enables things like getting data
from a vector directly, without having R in the middle, which would be
very slow.
- `refcount`: Used for the reference counting I rolled using Adam Ruppe's
book.

`RDataStorage` is wrapped by

```
struct RData {
  RDataStorage * data;
  alias data this;
```

This is done for purposes of reference counting. If we didn't use reference
counting, we'd allocate data inside R but never release it, which would
obviously eventually lead to running out of memory. When the reference
count hits zero, `rm` is called inside R.

The most common usage of `RData` is `RData(code)` where `code` is a
string of R code. An example is `auto x = RData(i"mean($(var.name))".text);`.


