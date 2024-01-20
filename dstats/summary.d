/**Summary statistics such as mean, median, sum, variance, skewness, kurtosis.
 * Except for median and median absolute deviation, which cannot be calculated
 * online, all summary statistics have both an input range interface and an
 * output range interface.
 *
 * Notes: The put method on the structs defined in this module returns this by
 *        ref.  The use case for returning this is to enable these structs
 *        to be used with std.algorithm.reduce.  The rationale for returning
 *        by ref is that the return value usually won't be used, and the
 *        overhead of returning a large struct by value should be avoided.
 *
 * Bugs:  This whole module assumes that input will be doubles or types implicitly
 *        convertible to double.  No allowances are made for user-defined numeric
 *        types such as BigInts.  This is necessary for simplicity.  However,
 *        if you have a function that converts your data to doubles, most of
 *        these functions work with any input range, so you can simply map
 *        this function onto your range.
 *
 * Author:  David Simcha
 */
/*
 * Copyright (C) 2008-2010 David Simcha
 *
 * License:
 * Boost Software License - Version 1.0 - August 17th, 2003
 *
 * Permission is hereby granted, free of charge, to any person or organization
 * obtaining a copy of the software and accompanying documentation covered by
 * this license (the "Software") to use, reproduce, display, distribute,
 * execute, and transmit the Software, and to prepare derivative works of the
 * Software, and to permit third-parties to whom the Software is furnished to
 * do so, all subject to the following:
 *
 * The copyright notices in the Software and this entire statement, including
 * the above license grant, this restriction and the following disclaimer,
 * must be included in all copies of the Software, in whole or in part, and
 * all derivative works of the Software, unless such copies or derivative
 * works are solely in the form of machine-executable object code generated by
 * a source language processor.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
 * SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
 * FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 */


module dstats.summary;

import std.functional, std.conv, std.range, std.array,
    std.traits, std.math, std.typetuple;

import std.algorithm : reduce, min, max, swap, map, filter;

/**Finds median of an input range in O(N) time on average.  In the case of an
 * even number of elements, the mean of the two middle elements is returned.
 * This is a convenience founction designed specifically for numeric types,
 * where the averaging of the two middle elements is desired.  A more general
 * selection algorithm that can handle any type with a total ordering, as well
 * as selecting any position in the ordering, can be found at
 * dstats.sort.quickSelect() and dstats.sort.partitionK().
 * Allocates memory, does not reorder input data.*/
// This is the old code. Changing to the double[] version below eliminates
// the dependency on dstats.alloc altogether.
//~ double median(T)(T data)
//~ if(doubleInput!(T)) {
    //~ auto alloc = newRegionAllocator();
    //~ auto dataDup = alloc.array(data);
    //~ return medianPartition(dataDup);
//~ }

double median(double[] x) {
	return medianPartition(x);
}

/**Median finding as in median(), but will partition input data such that
 * elements less than the median will have smaller indices than that of the
 * median, and elements larger than the median will have larger indices than
 * that of the median. Useful both for its partititioning and to avoid
 * memory allocations.  Requires a random access range with swappable
 * elements.*/
double medianPartition(T)(T data)
if(isRandomAccessRange!(T) &&
   is(ElementType!(T) : double) &&
   hasSwappableElements!(T) &&
   hasLength!(T))
{
    if(data.length == 0) {
        return double.nan;
    }
    // Upper half of median in even length case is just the smallest element
    // with an index larger than the lower median, after the array is
    // partially sorted.
    if(data.length == 1) {
        return data[0];
    } else if(data.length & 1) {  //Is odd.
        return cast(double) partitionK(data, data.length / 2);
    } else {
        auto lower = partitionK(data, data.length / 2 - 1);
        auto upper = ElementType!(T).max;

        // Avoid requiring slicing to be supported.
        foreach(i; data.length / 2..data.length) {
            if(data[i] < upper) {
                upper = data[i];
            }
        }
        return lower * 0.5 + upper * 0.5;
    }
}

/**Plain old data holder struct for median, median absolute deviation.
 * Alias this'd to the median absolute deviation member.
 */
struct MedianAbsDev {
    double median;
    double medianAbsDev;

    alias medianAbsDev this;
}

/**Calculates the median absolute deviation of a dataset.  This is the median
 * of all absolute differences from the median of the dataset.
 *
 * Returns:  A MedianAbsDev struct that contains the median (since it is
 * computed anyhow) and the median absolute deviation.
 *
 * Notes:  No bias correction is used in this implementation, since using
 * one would require assumptions about the underlying distribution of the data.
 */
MedianAbsDev medianAbsDev(T)(T data)
if(doubleInput!(T)) {
    auto alloc = newRegionAllocator();
    auto dataDup = alloc.array(data);
    immutable med = medianPartition(dataDup);
    immutable len = dataDup.length;
    alloc.freeLast();

    double[] devs = alloc.uninitializedArray!(double[])(len);

    size_t i = 0;
    foreach(elem; data) {
        devs[i++] = abs(med - elem);
    }
    auto ret = medianPartition(devs);
    alloc.freeLast();
    return MedianAbsDev(med, ret);
}

/**Output range to calculate the mean online.  Getter for mean costs a branch to
 * check for N == 0.  This struct uses O(1) space and does *NOT* store the
 * individual elements.
 *
 * Note:  This struct can implicitly convert to the value of the mean.
 *
 * Examples:
 * ---
 * Mean summ;
 * summ.put(1);
 * summ.put(2);
 * summ.put(3);
 * summ.put(4);
 * summ.put(5);
 * assert(summ.mean == 3);
 * ---*/
struct Mean {
private:
    double result = 0;
    double k = 0;

public:
    ///// Allow implicit casting to double, by returning the current mean.
    alias mean this;

    ///
    void put(double element) pure nothrow @safe {
        result += (element - result) / ++k;
    }

    /**Adds the contents of rhs to this instance.
     *
     * Examples:
     * ---
     * Mean mean1, mean2, combined;
     * foreach(i; 0..5) {
     *     mean1.put(i);
     * }
     *
     * foreach(i; 5..10) {
     *     mean2.put(i);
     * }
     *
     * mean1.put(mean2);
     *
     * foreach(i; 0..10) {
     *     combined.put(i);
     * }
     *
     * assert(approxEqual(combined.mean, mean1.mean));
     * ---
     */
     void put(typeof(this) rhs) pure nothrow @safe {
         immutable totalN = k + rhs.k;
         result = result * (k / totalN) + rhs.result * (rhs.k / totalN);
         k = totalN;
     }

    const pure nothrow @property @safe {

        ///
        double sum() {
            return result * k;
        }

        ///
        double mean() {
            return (k == 0) ? double.nan : result;
        }

        ///
        double N() {
            return k;
        }

        /**Simply returns this.  Useful in generic programming contexts.*/
        Mean toMean() {
            return this;
        }
    }

    ///
    string toString() const {
        return to!(string)(mean);
    }
}

/**Finds the arithmetic mean of any input range whose elements are implicitly
 * convertible to double.*/
Mean mean(T)(T data)
if(doubleIterable!(T)) {

    static if(isRandomAccessRange!T && hasLength!T) {
        // This is optimized for maximum instruction level parallelism:
        // The loop is unrolled such that there are 1 / (nILP)th the data
        // dependencies of the naive algorithm.
        enum nILP = 8;

        Mean ret;
        size_t i = 0;
        if(data.length > 2 * nILP) {
            double k = 0;
            double[nILP] means = 0;
            for(; i + nILP < data.length; i += nILP) {
                immutable kNeg1 = 1 / ++k;

                foreach(j; StaticIota!nILP) {
                    means[j] += (data[i + j] - means[j]) * kNeg1;
                }
            }

            ret.k = k;
            ret.result = means[0];
            foreach(m; means[1..$]) {
                ret.put( Mean(m, k));
            }
        }

        // Handle the remainder.
        for(; i < data.length; i++) {
            ret.put(data[i]);
        }
        return ret;

    } else {
        // Just submit everything to a single Mean struct and return it.
        Mean meanCalc;

        foreach(element; data) {
            meanCalc.put(element);
        }
        return meanCalc;
    }
}

/**Finds the sum of an input range whose elements implicitly convert to double.
 * User has option of making U a different type than T to prevent overflows
 * on large array summing operations.  However, by default, return type is
 * T (same as input type).*/
U sum(T, U = Unqual!(ForeachType!(T)))(T data)
if(doubleIterable!(T)) {

    static if(isRandomAccessRange!T && hasLength!T) {
        enum nILP = 8;
        U[nILP] sum = 0;

        size_t i = 0;
        if(data.length > 2 * nILP) {

            for(; i + nILP < data.length; i += nILP) {
                foreach(j; StaticIota!nILP) {
                    sum[j] += data[i + j];
                }
            }

            foreach(j; 1..nILP) {
                sum[0] += sum[j];
            }
        }

        for(; i < data.length; i++) {
            sum[0] += data[i];
        }

        return sum[0];
    } else {
        U sum = 0;
        foreach(elem; data) {
            sum += elem;
        }

        return sum;
    }
}

/**Output range to compute mean, stdev, variance online.  Getter methods
 * for stdev, var cost a few floating point ops.  Getter for mean costs
 * a single branch to check for N == 0.  Relatively expensive floating point
 * ops, if you only need mean, try Mean.  This struct uses O(1) space and
 * does *NOT* store the individual elements.
 *
 * Note:  This struct can implicitly convert to a Mean struct.
 *
 * References: Computing Higher-Order Moments Online.
 * http://people.xiph.org/~tterribe/notes/homs.html
 *
 * Examples:
 * ---
 * MeanSD summ;
 * summ.put(1);
 * summ.put(2);
 * summ.put(3);
 * summ.put(4);
 * summ.put(5);
 * assert(summ.mean == 3);
 * assert(summ.stdev == sqrt(2.5));
 * assert(summ.var == 2.5);
 * ---*/
struct MeanSD {
private:
    double _mean = 0;
    double _var = 0;
    double _k = 0;
public:
    ///
    void put(double element) pure nothrow @safe {
        immutable kMinus1 = _k;
        immutable delta = element - _mean;
        immutable deltaN = delta / ++_k;

        _mean += deltaN;
        _var += kMinus1 * deltaN * delta;
        return;
    }

    /// Combine two MeanSD's.
    void put(typeof(this) rhs) pure nothrow @safe {
        if(_k == 0) {
            foreach(ti, elem; rhs.tupleof) {
                this.tupleof[ti] = elem;
            }

            return;
        } else if(rhs._k == 0) {
            return;
        }

        immutable totalN = _k + rhs._k;
        immutable delta = rhs._mean - _mean;
        _mean = _mean * (_k / totalN) + rhs._mean * (rhs._k / totalN);

        _var = _var + rhs._var + (_k / totalN * rhs._k * delta * delta);
        _k = totalN;
    }

    const pure nothrow @property @safe {

        ///
        double sum() {
            return _k * _mean;
        }

        ///
        double mean() {
            return (_k == 0) ? double.nan : _mean;
        }

        ///
        double stdev() {
            return sqrt(var);
        }

        ///
        double var() {
            return (_k < 2) ? double.nan : _var / (_k - 1);
        }

        /**
        Mean squared error.  In other words, a biased estimate of variance.
        */
        double mse() {
            return (_k < 1) ? double.nan : _var / _k;
        }

        ///
        double N() {
            return _k;
        }

        /**Converts this struct to a Mean struct.  Also called when an
         * implicit conversion via alias this takes place.
         */
        Mean toMean() {
            return Mean(_mean, _k);
        }

        /**Simply returns this.  Useful in generic programming contexts.*/
        MeanSD toMeanSD() const  {
            return this;
        }
    }

    alias toMean this;

    ///
    string toString() const {
        return text("N = ", cast(ulong) _k, "\nMean = ", mean, "\nVariance = ",
               var, "\nStdev = ", stdev);
    }
}

/**Puts all elements of data into a MeanSD struct,
 * then returns this struct.  This can be faster than doing this manually
 * due to ILP optimizations.*/
MeanSD meanStdev(T)(T data)
if(doubleIterable!(T)) {

    MeanSD ret;

    static if(isRandomAccessRange!T && hasLength!T) {
        // Optimize for instruction level parallelism.
        enum nILP = 6;
        double k = 0;
        double[nILP] means = 0;
        double[nILP] variances = 0;
        size_t i = 0;

        if(data.length > 2 * nILP) {
            for(; i + nILP < data.length; i += nILP) {
                immutable kMinus1 = k;
                immutable kNeg1 = 1 / ++k;

                foreach(j; StaticIota!nILP) {
                    immutable double delta = data[i + j] - means[j];
                    immutable deltaN = delta * kNeg1;

                    means[j] += deltaN;
                    variances[j] += kMinus1 * deltaN * delta;
                }
            }

            ret._mean = means[0];
            ret._var = variances[0];
            ret._k = k;

            foreach(j; 1..nILP) {
                ret.put( MeanSD(means[j], variances[j], k));
            }
        }

        // Handle remainder.
        for(; i < data.length; i++) {
            ret.put(data[i]);
        }
    } else {
        foreach(elem; data) {
            ret.put(elem);
        }
    }
    return ret;
}

/**Finds the variance of an input range with members implicitly convertible
 * to doubles.*/
double variance(T)(T data)
if(doubleIterable!(T)) {
    return meanStdev(data).var;
}

/**Calculate the standard deviation of an input range with members
 * implicitly converitble to double.*/
double stdev(T)(T data)
if(doubleIterable!(T)) {
    return meanStdev(data).stdev;
}

/* This is all that was needed from dstats.sort. */


/* Returns the index, NOT the value, of the median of the first, middle, last
 * elements of data.*/
size_t medianOf3(alias compFun, T)(T[] data) {
    alias binaryFun!(compFun) comp;
    immutable size_t mid = data.length / 2;
    immutable uint result = ((cast(uint) (comp(data[0], data[mid]))) << 2) |
                            ((cast(uint) (comp(data[0], data[$ - 1]))) << 1) |
                            (cast(uint) (comp(data[mid], data[$ - 1])));

    assert(result != 2 && result != 5 && result < 8); // Cases 2, 5 can't happen.
    switch(result) {
        case 1:  // 001
        case 6:  // 110
            return data.length - 1;
        case 3:  // 011
        case 4:  // 100
            return 0;
        case 0:  // 000
        case 7:  // 111
            return mid;
        default:
            assert(0);
    }
    assert(0);
}

/**Partitions the input data according to compFun, such that position k contains
 * the kth largest/smallest element according to compFun.  For all elements e
 * with indices < k, !compFun(data[k], e) is guaranteed to be true.  For all
 * elements e with indices > k, !compFun(e, data[k]) is guaranteed to be true.
 * For example, if compFun is "a < b", all elements with indices < k will be
 * <= data[k], and all elements with indices larger than k will be >= k.
 * Reorders any additional input arrays in lockstep.
 *
 * Examples:
 * ---
 * auto foo = [3, 1, 5, 4, 2].dup;
 * auto secondSmallest = partitionK(foo, 1);
 * assert(secondSmallest == 2);
 * foreach(elem; foo[0..1]) {
 *     assert(elem <= foo[1]);
 * }
 * foreach(elem; foo[2..$]) {
 *     assert(elem >= foo[1]);
 * }
 * ---
 *
 * Returns:  The kth element of the array.
 */
ElementType!(T[0]) partitionK(alias compFun = "a < b", T...)(T data, ptrdiff_t k)
in {
    assert(data.length > 0);
    size_t len = data[0].length;
    foreach(array; data[1..$]) {
        assert(array.length == len);
    }
} do {
    // Don't use the float-to-int trick because it's actually slower here
    // because the main part of the algorithm is O(N), not O(N log N).
    return partitionKImpl!compFun(data, k);
}

/*private*/ ElementType!(T[0]) partitionKImpl(alias compFun, T...)(T data, ptrdiff_t k) {
    alias binaryFun!(compFun) comp;

    {
        immutable size_t med3 = medianOf3!(comp)(data[0]);
        foreach(array; data) {
            auto temp = array[med3];
            array[med3] = array[$ - 1];
            array[$ - 1] = temp;
        }
    }

    ptrdiff_t lessI = -1, greaterI = data[0].length - 1;
    auto pivot = data[0][$ - 1];
    while(true) {
        while(comp(data[0][++lessI], pivot)) {}
        while(greaterI > 0 && comp(pivot, data[0][--greaterI])) {}

        if(lessI < greaterI) {
            foreach(array; data) {
                auto temp = array[lessI];
                array[lessI] = array[greaterI];
                array[greaterI] = temp;
            }
        } else break;
    }
    foreach(array; data) {
        auto temp = array[lessI];
        array[lessI] = array[$ - 1];
        array[$ - 1] = temp;
    }

    if((greaterI < k && lessI >= k) || lessI == k) {
        return data[0][k];
    } else if(lessI < k) {
        foreach(ti, array; data) {
            data[ti] = array[lessI + 1..$];
        }
        return partitionK!(compFun, T)(data, k - lessI - 1);
    } else {
        foreach(ti, array; data) {
            data[ti] = array[0..min(greaterI + 1, lessI)];
        }
        return partitionK!(compFun, T)(data, k);
    }
}


/* This is all that was needed from dstats.base */


// Returns the number of dimensions in an array T.
package template nDims(T)
{
    static if(isArray!T)
    {
        enum nDims = 1 + nDims!(typeof(T.init[0]));
    }
    else
    {
        enum nDims = 0;
    }
}

/**
Tests whether T is an input range whose elements can be implicitly
converted to doubles.*/
template doubleInput(T) {
    enum doubleInput = isInputRange!(T) && is(ElementType!(T) : double);
}

/**Tests whether T is iterable and has elements of a type implicitly
 * convertible to double.*/
template doubleIterable(T) {
    static if(!isIterable!T) {
        enum doubleIterable = false;
    } else {
        enum doubleIterable = is(ForeachType!(T) : double);
    }
}

// Used for ILP optimizations.
template StaticIota(size_t upTo) {
    static if(upTo == 0) {
        alias TypeTuple!() StaticIota;
    } else {
        alias TypeTuple!(StaticIota!(upTo - 1), upTo - 1) StaticIota;
    }
}

// Uses Gauss-Jordan elim. w/ row pivoting to invert from.  Stores the results
// in to and leaves from in an undefined state.
void invert(double[][] from, double[][] to) {
		// Normalize.
		foreach(i, row; from) {
				double absMax = 1.0 / reduce!(max)(map!(abs)(row[0..from.length]));
				row[] *= absMax;
				to[i][] = 0;
				to[i][i] = absMax;
		}

		foreach(col; 0..from.length) {
				size_t bestRow;
				double biggest = 0;
				foreach(row; col..from.length) {
						if(abs(from[row][col]) > biggest) {
								bestRow = row;
								biggest = abs(from[row][col]);
						}
				}

				swap(from[col], from[bestRow]);
				swap(to[col], to[bestRow]);
				immutable pivotFactor = from[col][col];

				foreach(row; 0..from.length) if(row != col) {
						immutable ratio = from[row][col] / pivotFactor;

						// If you're ever looking to optimize this code, good luck.  The
						// bottleneck is almost ENTIRELY this one line:
						from[row][] -= from[col][] * ratio;
						to[row][] -= to[col][] * ratio;
				}
		}

		foreach(i; 0..from.length) {
				immutable diagVal = from[i][i];
				from[i][] /= diagVal;
				to[i][] /= diagVal;
		}
}

