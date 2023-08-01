# List

Basic facilities for working with an R list.

## Constructing

### this(string code)

Executes `code` inside R and stores the result as a List.

### this(List x)

Puts a copy of x in a new list.

### this(long n)

Allocates a new list with n elements. Not important in practice, because
R will automatically resize the list as necessary.
  
### this(Robj obj)

Mostly for internal use.

### this(RData rd)

Mostly for internal use.

## Indexing

You can index by name or number.

```
x["integer element"] = 3;
x[2] = "Third element of the list";
```

Note that when indexing like this, the return type is unknown, so you always
get an element of type `RData`. To specify the type, as you'll usually want
to do, use the `as` method:

```
auto j = x["integer element"].as!int;
auto m = x[3].as!Matrix;
```

## Getting multiple elements

You can do the usual slicing. It returns a new list holding the elements.

```
List y = x[0..3]; // First three elements of x
List z = x[2..$]; // Dollar operator works
```

Alternatively, if you want to grab multiple nonconsecutive elements, or
multiple elements by name, you can do this:

```
List y = x[[0, 1, 2]];
List z = x[["One Element's Name", "Another Element's Name"]];
```

## Duplicating

You can use `dup` to make a copy:

```
List y = x.dup();
```

## Setting elements

Just as you can get elements by name or index number, you can also set
elements that way.

```
x["Item Name"] = rnorm(100);
y[4] = 6.2;
```

## Names

You can get a string[] holding all names in the list by calling `names`:

```
x["Random"] = rnorm(100);
x["Not Random"] = 6.2;
writeln(x.names);
```

