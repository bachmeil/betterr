# Warning

This is an experimental module. What's here will (probably) work correctly, but it's all subject to breaking changes as I decide what works best for my use cases, so you should expect that any code you write now will fail to compile in the future. You can avoid issues by storing a copy of this file in the project directory with the rest of your source files. Also, there's no guarantee that this documentation will be accurate or complete while I'm designing the module. In some cases I write documentation before a feature has been added.

# struct TS

This struct holds a single time series. Although the underlying data is stored as a vector, it's required to know the start date, end data, and frequency properties. If you don't have all that information, the time series is not fully defined, and there were will be some operations you can't do. For instance, taking the mean of a single time series requires only the vector of data, but taking the lag is only meaningful if you know the frequency and start date. The data array is the same for both the series and its lags. What changes is the interpretation of the dates attached to observations. Operations like summing time series only make sense if they are the same frequency and cover the same dates.

Since the frequency is required, it's a required template parameter. There is currently support for the three frequencies I use most often:

- Annual: `TS!1` or `AnnualTS`
- Quarterly: `TS!4` or `QuarterlyTS`
- Monthly: `TS!12` or `MonthlyTS`

Other frequencies will be supported in the future using the facilities provided by R. One thing to note about these frequencies is that the type of the date depends on the frequency. An annual time series has an integer date such as 2014. Quarterly and monthly time series have dates that are two integers, one for the year and one for the minor item, so that January 2014 and 2014Q1 are represented by the date [2014,1].

Indexing of time series is done by date. For an annual time series, y[2014] is the observation corresponding to 2014, not element number 2014. For a monthly time series, y[2014, 1] is the observation corresponding to January 2014. Similarly, slicing is done by dates. An annual time series would be sliced using y[2014..2023]. A quarterly time series would be sliced using y[[2014,1]..[2023,4]]. If for some reason you actually want to work with the data array, you're discarding the time series properties, and the efficient way to do that is to call y.array, which returns a double[] with no allocation or copying of data.

# Construction

## this(string code)

This is primarily used internally. `code` is passed to R and the output is captured. It is assumed that the code run inside R produces a ts object. For instance, if you had a vector x inside R and you wanted to convert it to a TS struct accessible from D, you could do something like this:

```
auto xx = TS("ts(x)");
```

## this(string code, long _start)

Converts `code` to an annual ts object, with the first observation being `_start`. If you have vector x inside R and want to convert it to a TS struct starting in 1959:

```
auto xx = TS("x", 1959);
```

## this(string code, long[2] _start, int _frequency)

If you have vector x inside R and want to convert it to a monthly TS struct starting in January 1959:

```
auto xx = TS("x", [1959, 1], 12);
```

## this(T)(T obj, long _start)

Convert the D struct obj into an annual time series starting in `_start`.

```
auto v = rnorm(20);
auto x = TS(v, 2003);
```

## this(T)(T obj, long[2] _start, int _frequency)

Convert the D struct obj into a monthly or quarterly time series starting in `_start`.

```
auto v = rnorm(120);
auto x = TS(v, [2013, 1]);
```

# Indexing

Indexing of TS elements is done by date rather than element number. The underlying data array changes routinely for TS structs, as you take lags, differences, percentage changes, etc., so you want to use dates instead. (If you really want to work with the underlying data array, you can do that as ptr[ii], since ptr is a pointer to the data array.)

## One element of an annual series

```
double v = x[1987];
```

## One element of a monthly or quarterly series

```
double v = x[2017, 4];
```

## Slicing an annual series

Note: The last date is included.

```
TS y = x[1980..2022];
```

## Slicing a monthly or quarterly series

Note: The last date is included.

```
TS y = x[ [1980,1]..[2022,6] ];
```

## until(long e)

Returns all observations of an annual series through year e.

```
TS y = x.until(2008);
```

## TS until(long[2] e)

Returns all observations of a quarterly or monthly series through e.

```
TS y = x.until(2008, 8);
```

## TS starting(long s)

Returns all observations of an annual series from year s to the end.

```
TS y = x.starting(2008);
```

## TS starting(long[2] s)

Returns all observations of a quarterly or monthly series from year s to the end.

```
TS y = x.starting(2008, 8);
```

## TS lag(long k)

Returns the TS created by taking the kth lag. Creates a new TS struct without mutating this.

Note that unlike R, this uses the standard definition of a lag, such that lag 1 is the previous observation.

```
TS y = x.lag(4); // Year earlier observations of a quarterly series
```

## TS lead(long k)

Returns the series from k time periods in the future. This matches R's non-standard definition of a lag, so that `x.lead(2)` is the equivalent of `lag(x, 2)` in R.

```
TS y = x.lead(4); // Year in the future observations of a quarterly series
```

## TS diff(long k)

Difference between the values of x at time t and t-k. Equivalent to R's diff function.

```
TS dx = x.diff(1);
```

## TS pct(long k)

Calculates the percentage change from time t-k to time t. 

## void print(string msg="")

Prints the series and metadata to the screen, optionally preceded by the message msg.

```
x.print("Monthly copper prices from January 1990 to March 2023");
```

# Transforming TS Objects

There are basic facilities for taking transformations of a TS struct
(lag, lead, difference, etc.) If you have a few series of a few thousand
observations, you can take the easy route and bind TS objects together
to get an MTS object. That may end up being too slow for a simulation
or bootstrap where you're combining many time series into an MTS struct.

To provide for greater efficiency, you can define a function that returns
a TSTransform struct:

```
struct TSTransform {
	void delegate(ref double[], long, long) compute;
	/* First non-missing observation available for the transformed series */
	long modStart;
	/* Last non-missing observation available for the transformed series */
	long modEnd;
}
```

modStart and modEnd are straightforward; you just calculate how the
start and end dates of the series change with the transformation.

The compute delegate is a little more complicated. You have to create
this delegate and return it in the TSTransform struct. The first argument
is a `double[]` that points to a column of the mts object allocated by
R that will be storing the data. The second is the date of the first
observation to do the transformation, and the third is the date of the
last observation to do the transformation.

It may be more clear looking at how this is implemented for a lag
transformation. The user specifies the order of the lag (1, 2, or
whatever). The construction function looks like this:

```
TSTransform Lag(long f)(TS!f var, long k=1) {
	double[] source = var.array;
	long arrayStart = var.longStart;
	long arrayEnd = var.longEnd;
	
	void compute(ref double[] target, long s, long e) {
		target[0..$] = source[(s-k-arrayStart)..(e+1-k-arrayStart)];
	}
	
	return TSTransform(&compute, arrayStart+k, arrayEnd+k);
}
```

The user passes a TS struct of arbitrary frequency and the lag order, k.
There are three pieces of data stored in the delegate `compute`. The
underlying data array of the TS is stored in source. The start and end
dates of the TS are stored in arrayStart and arrayEnd, respectively.

The compute function does the computation of the kth lag of TS var from
time s to time e. In this case, it fills target with the appropriate
elements from source.

The MTSTransform struct holds all of the transformations that will be in
the dataset after storing all of the transformed data series.

As an example, suppose you have a TS named `x`, and you want an MTS
holding the first and second lags of `x`. You do that with these three
lines of code:

```
MTSTransform transform;
transform.data = [Lag(x, 1), Lag(x, 2)];
MTS newts = transform.create();
```

This is a simple and efficient approach for combining many transformed
time series into one MTS struct when you can use the built-in
transformations. It's general enough that you can do *any* transformation.
For instance, you can create a construction function that takes a group
of TS series and creates a new series that is the difference between
two of them. You can access arbitrary data from inside the delegate.







