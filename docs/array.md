# RArray

An RArray is a struct that holds an n-dimensional array allocated inside R. You can see the functionality it provides by reading this code:

```
auto arr = RArray(2,3,4,2);
foreach(ii; 0..48) {
  arr.ptr[ii] = double(ii);
}
```

`arr` is a 4-dimensional array with dimensions 2x3x4x2, for 48 total elements. The constructor takes an group of arguments of any length.

```
printR(arr.x);
writeln(arr);
```

You can print the array using the low-level interface (`printR`) or the high-level interface (`writeln`).

```
writeln(arr[1,1,1,1]);
writeln(arr[1,2,3,1]);
writeln(arr[1,1,3,0]);
```

You can pull out individual elements using the usual notation.

```
arr[1,1,1,1] = 33.7;
```

You can set elements using the usual notation.

```
arr = 2.7;
```

You can set all elements of an array to the same value.

```
arr = rnorm(48);
```

A check will be done to be sure the dimensions match. Since `arr` has 48
elements, this operation works. This would throw an error:

```
arr = rnorm(47);
```

You cannot assign to an array incompletely; the above will not fill the first 47 elements of `arr`.

```
auto arr2 = arr.sub[0, 0, 0, _all];
```

The combination of a variable number of arguments and multidimensional slicing means you have to create a
reference to the array in order to do multidimensional slicing if you want to be able to assign to a
multidimensional slice. `.sub` returns a reference to the array.

```
auto arr3 = arr.sub[0, 0..3, 0, _all];
arr.sub[0..2, 0, 0, 0] = arr.sub[0..2, 1, 1, 1];
```

`_all` means all elements on that dimension. `0..2` means the first two elements of that dimension.

```
auto m = Matrix([1.1, 2, 3, 4, 5, 6], 2, 3);
arr.sub[0..2, 0..3, 1, 1] = m;
```

You can assign other structs, such as a Matrix, to a subset of an array, so long as the dimensions match.
