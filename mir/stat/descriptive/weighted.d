/++
This module contains algorithms for descriptive statistics with weights.

License: $(HTTP www.apache.org/licenses/LICENSE-2.0, Apache-2.0)

Authors: John Michael Hall

Copyright: 2022 Mir Stat Authors.

Macros:
SUBREF = $(REF_ALTTEXT $(TT $2), $2, mir, stat, $1)$(NBSP)
MATHREF = $(GREF_ALTTEXT mir-algorithm, $(TT $2), $2, mir, math, $1)$(NBSP)
NDSLICEREF = $(GREF_ALTTEXT mir-algorithm, $(TT $2), $2, mir, ndslice, $1)$(NBSP)
T2=$(TR $(TDNW $(LREF $1)) $(TD $+))
T3=$(TR $(TDNW $(LREF $1)) $(TD $2) $(TD $+))
T4=$(TR $(TDNW $(LREF $1)) $(TD $2) $(TD $3) $(TD $4))

+/

module mir.stat.descriptive.weighted;

import mir.math.sum: ResolveSummationType, Summation, Summator;

private void putter2(Slices, T, U, Summation summation1, Summation summation2)
    (scope Slices slices, ref Summator!(T, summation1) seed1, ref Summator!(U, summation2) seed2)
{
    import mir.functional: Tuple;
    static if (is(Slices == Tuple!(V1, V2), V1, V2)) {
        seed1.put(slices[0]);
        seed2.put(slices[1]);
    } else {
        import mir.ndslice.internal: frontOf2;
        do
        {
            frontOf2!(slices)[0].putter2(seed1, seed2);
            slices.popFront;
        }
        while(!slices.empty);
    }
}

/++
Assumptions used for weighted moments
+/
enum AssumeWeights : bool
{
    /++
    Primary, does not assume weights sum to one
    +/
    primary,
    
    /++
    Assumes weights sum to one
    +/
    sumToOne
}

/++
Output range for wmean.
+/
struct WMeanAccumulator(T, Summation summation, AssumeWeights assumeWeights,
                        U = T, Summation weightsSummation = summation)
{
    import mir.ndslice.slice: isConvertibleToSlice, isSlice, kindOf;
    import std.range.primitives: isInputRange;
    import std.traits: isIterable;

    ///
    Summator!(T, summation) wsummator;

    static if (!assumeWeights) {
        ///
        Summator!(U, weightsSummation) weights;
    }

    ///
    F wmean(F = T)() const @safe @property pure nothrow @nogc
    {
        static if (assumeWeights) {
            return this.wsum!F;
        } else {
            assert(this.weight!F != 0, "weight must not equal zero");
            return this.wsum!F / this.weight!F;
        }
    }

    ///
    F wsum(F = T)() const @safe @property pure nothrow @nogc
    {
        return cast(F) wsummator.sum;
    }

    ///
    F weight(F = U)() const @safe @property pure nothrow @nogc
    {
        return cast(F) weights.sum;
    }

    ///
    void put(Slice1, Slice2)(Slice1 s, Slice2 w)
        if (isSlice!Slice1 && isSlice!Slice2)
    {
        static assert (Slice1.N == Slice2.N, "s and w must have the same number of dimensions");
        static assert (kindOf!Slice1 == kindOf!Slice2, "s and w must have the same kind");

        import mir.ndslice.slice: Contiguous;
        import mir.ndslice.topology: zip, map;

        assert(s._lengths == w._lengths, "WMeanAcumulator.put: both slices must have the same lengths");

        static if (kindOf!Slice1 != Contiguous && Slice1.N > 1) {
            assert(s.strides == w.strides, "WMeanAccumulator.put: cannot put canonical and universal slices when strides do not match");
            auto combine = s.zip!true(w);
        } else {
            auto combine = s.zip!false(w);
        }

        static if (assumeWeights) {
            auto combine2 = combine.map!"a * b";
            wsummator.put(combine2);
        } else {
            auto combine2 = combine.map!("b", "a * b");
            combine2.putter2(weights, wsummator);
        }
    }

    ///
    void put(SliceLike1, SliceLike2)(SliceLike1 s, SliceLike2 w)
        if (isConvertibleToSlice!SliceLike1 && !isSlice!SliceLike1 &&
            isConvertibleToSlice!SliceLike2 && !isSlice!SliceLike2)
    {
        import mir.ndslice.slice: toSlice;
        this.put(s.toSlice, w.toSlice);
    }

    ///
    void put(Range)(Range r)
        if (isIterable!Range && !assumeWeights)
    {
        import mir.primitives: hasShape, elementCount;
        static if (hasShape!Range) {
            wsummator.put(r);
            weights.put(cast(U) r.elementCount);
        } else {
            foreach(x; r)
            {
                this.put(x);
            }
        }
    }

    ///
    void put(RangeA, RangeB)(RangeA r, RangeB w)
        if (isInputRange!RangeA && !isConvertibleToSlice!RangeA &&
            isInputRange!RangeB && !isConvertibleToSlice!RangeB)
    {
        do
        {
            assert(!(!r.empty && w.empty) && !(r.empty && !w.empty),
                   "r and w must both be empty at the same time, one cannot be empty while the other has remaining items");
            this.put(r.front, w.front);
            r.popFront;
            w.popFront;
        } while(!r.empty || !w.empty); // Using an || instead of && so that the loop does not end early. mis-matched lengths of r and w sould be caught by above assert
    }

    ///
    void put()(T x, U w)
    {
        static if (!assumeWeights) {
            weights.put(w);
        }
        wsummator.put(x * w);
    }

    ///
    void put()(T x)
        if (!assumeWeights)
    {
        weights.put(cast(U) 1);
        wsummator.put(x);
    }

    ///
    void put(F = T, G = U)(WMeanAccumulator!(F, summation, assumeWeights, G, weightsSummation) wm)
        if (!assumeWeights) // because calculating is easier. When assumeWeightsSumtoOne = true, need to divide original wsummator and wm by 2.
    {
        weights.put(cast(U) wm.weights);
        wsummator.put(cast(T) wm.wsummator);
    }
}

/// Assume weights sum to 1
version(mir_stat_test)
@safe pure nothrow
unittest
{
    import mir.math.sum: Summation;
    import mir.ndslice.slice: sliced;
    import mir.test: should;

    WMeanAccumulator!(double, Summation.pairwise, AssumeWeights.sumToOne) x;
    x.put([0.0, 1, 2, 3, 4].sliced, [0.2, 0.2, 0.2, 0.2, 0.2].sliced);
    x.wmean.should == 2;
    x.put(5, 0.0);
    x.wmean.should == 2;
}

// dynamic array test, assume weights sum to 1
version(mir_stat_test)
@safe pure nothrow
unittest
{
    import mir.math.sum: Summation;
    import mir.test: should;

    WMeanAccumulator!(double, Summation.pairwise, AssumeWeights.sumToOne) x;
    x.put([0.0, 1, 2, 3, 4], [0.2, 0.2, 0.2, 0.2, 0.2]);
    x.wmean.should == 2;
}

// static array test, assume weights sum to 1
version(mir_stat_test)
@safe pure nothrow @nogc
unittest
{
    import mir.math.sum: Summation;
    import mir.test: should;

    WMeanAccumulator!(double, Summation.pairwise, AssumeWeights.sumToOne) x;
    static immutable y = [0.0, 1, 2, 3, 4];
    static immutable w = [0.2, 0.2, 0.2, 0.2, 0.2];
    x.put(y, w);
    x.wmean.should == 2;
}

// 2-d slice test, assume weights sum to 1
version(mir_stat_test)
@safe pure
unittest
{
    import mir.math.sum: Summation;
    import mir.ndslice.fuse: fuse;
    import mir.test: shouldApprox;

    WMeanAccumulator!(double, Summation.pairwise, AssumeWeights.sumToOne) x;
    auto y = [
        [0.0, 1, 2],
        [3.0, 4, 5]
    ].fuse;
    auto w = [
        [1.0 / 21, 2.0 / 21, 3.0 / 21],
        [4.0 / 21, 5.0 / 21, 6.0 / 21]
    ].fuse;
    x.put(y, w);
    x.wmean.shouldApprox == 70.0 / 21;
}

// universal 2-d slice test, assume weights sum to 1, using map
version(mir_stat_test)
@safe pure nothrow
unittest
{
    import mir.math.sum: Summation;
    import mir.ndslice.topology: iota, map, universal;
    import mir.test: shouldApprox;

    WMeanAccumulator!(double, Summation.pairwise, AssumeWeights.sumToOne) x;
    auto y = iota([2, 3]).universal;
    auto w = iota([2, 3], 1).map!(a => a / 21.0).universal;
    x.put(y, w);
    x.wmean.shouldApprox == 70.0 / 21;
}

// 2-d canonical slice test, assume weights sum to 1, using map
version(mir_stat_test)
@safe pure nothrow
unittest
{
    import mir.math.sum: Summation;
    import mir.ndslice.topology: canonical, iota, map;
    import mir.test: shouldApprox;

    WMeanAccumulator!(double, Summation.pairwise, AssumeWeights.sumToOne) x;
    auto y = iota([2, 3]).canonical;
    auto w = iota([2, 3], 1).map!(a => a / 21.0).canonical;
    x.put(y, w);
    x.wmean.shouldApprox == 70.0 / 21;
}

/// Do not assume weights sum to 1
version(mir_stat_test)
@safe pure nothrow
unittest
{
    import mir.math.sum: Summation;
    import mir.ndslice.slice: sliced;
    import mir.test: shouldApprox;

    WMeanAccumulator!(double, Summation.pairwise, AssumeWeights.primary) x;
    x.put([0.0, 1, 2, 3, 4].sliced, [1, 2, 3, 4, 5].sliced);
    x.wmean.shouldApprox == 40.0 / 15;
    x.put(5, 6);
    x.wmean.shouldApprox == 70.0 / 21;
}

// dynamic array test, do not assume weights sum to 1
version(mir_stat_test)
@safe pure nothrow
unittest
{
    import mir.math.sum: Summation;
    import mir.test: shouldApprox;

    WMeanAccumulator!(double, Summation.pairwise, AssumeWeights.primary) x;
    x.put([0.0, 1, 2, 3, 4], [1, 2, 3, 4, 5]);
    x.wmean.shouldApprox == 40.0 / 15;
}

// static array test, do not assume weights sum to 1
version(mir_stat_test)
@safe pure nothrow @nogc
unittest
{
    import mir.math.sum: Summation;
    import mir.test: shouldApprox;

    WMeanAccumulator!(double, Summation.pairwise, AssumeWeights.primary) x;
    static immutable y = [0.0, 1, 2, 3, 4];
    static immutable w = [1, 2, 3, 4, 5];
    x.put(y, w);
    x.wmean.shouldApprox == 40.0 / 15;
}

// 2-d slice test, do not assume weights sum to 1
version(mir_stat_test)
@safe pure
unittest
{
    import mir.math.sum: Summation;
    import mir.ndslice.fuse: fuse;
    import mir.test: shouldApprox;

    WMeanAccumulator!(double, Summation.pairwise, AssumeWeights.primary) x;
    auto y = [
        [0.0, 1, 2],
        [3.0, 4, 5]
    ].fuse;
    auto w = [
        [1.0, 2, 3],
        [4.0, 5, 6]
    ].fuse;
    x.put(y, w);
    x.wmean.shouldApprox == 70.0 / 21;
}

// universal slice test, do not assume weights sum to 1
version(mir_stat_test)
@safe pure nothrow
unittest
{
    import mir.math.sum: Summation;
    import mir.ndslice.topology: iota, universal;
    import mir.test: shouldApprox;

    WMeanAccumulator!(double, Summation.pairwise, AssumeWeights.primary) x;
    auto y = iota(6).universal;
    auto w = iota([6], 1).universal;
    x.put(y, w);
    x.wmean.shouldApprox == 70.0 / 21;
}

// canonical slice test, do not assume weights sum to 1
version(mir_stat_test)
@safe pure nothrow
unittest
{
    import mir.math.sum: Summation;
    import mir.ndslice.topology: canonical, iota;
    import mir.test: shouldApprox;

    WMeanAccumulator!(double, Summation.pairwise, AssumeWeights.primary) x;
    auto y = iota(6).canonical;
    auto w = iota([6], 1).canonical;
    x.put(y, w);
    x.wmean.shouldApprox == 70.0 / 21;
}

// 2-d universal slice test, do not assume weights sum to 1
version(mir_stat_test)
@safe pure nothrow
unittest
{
    import mir.math.sum: Summation;
    import mir.ndslice.topology: iota, universal;
    import mir.test: shouldApprox;

    WMeanAccumulator!(double, Summation.pairwise, AssumeWeights.primary) x;
    auto y = iota([2, 3]).universal;
    auto w = iota([2, 3], 1).universal;
    x.put(y, w);
    x.wmean.shouldApprox == 70.0 / 21;
}

// 2-d canonical slice test, do not assume weights sum to 1
version(mir_stat_test)
@safe pure nothrow
unittest
{
    import mir.math.sum: Summation;
    import mir.ndslice.topology: canonical, iota;
    import mir.test: shouldApprox;

    WMeanAccumulator!(double, Summation.pairwise, AssumeWeights.primary) x;
    auto y = iota([2, 3]).canonical;
    auto w = iota([2, 3], 1).canonical;
    x.put(y, w);
    x.wmean.shouldApprox == 70.0 / 21;
}

/// Assume no weights, like MeanAccumulator
version(mir_stat_test)
@safe pure nothrow
unittest
{
    import mir.math.sum: Summation;
    import mir.ndslice.slice: sliced;
    import mir.test: shouldApprox;

    WMeanAccumulator!(double, Summation.pairwise, AssumeWeights.primary) x;
    x.put([0.0, 1, 2, 3, 4].sliced);
    x.wmean.shouldApprox == 2;
    x.put(5);
    x.wmean.shouldApprox == 2.5;
}

// dynamic array test, assume no weights
version(mir_stat_test)
@safe pure nothrow
unittest
{
    import mir.math.sum: Summation;
    import mir.ndslice.slice: sliced;
    import mir.test: should;

    WMeanAccumulator!(double, Summation.pairwise, AssumeWeights.primary) x;
    x.put([0.0, 1, 2, 3, 4]);
    x.wmean.should == 2;
}

// static array test, assume no weights
version(mir_stat_test)
@safe pure nothrow @nogc
unittest
{
    import mir.math.sum: Summation;
    import mir.test: should;

    WMeanAccumulator!(double, Summation.pairwise, AssumeWeights.primary) x;
    static immutable y = [0.0, 1, 2, 3, 4];
    x.put(y);
    x.wmean.should == 2;
}

// Adding WMeanAccmulators
version(mir_stat_test)
@safe pure nothrow
unittest
{
    import mir.math.sum: Summation;
    import mir.test: shouldApprox;

    double[] x = [0.0, 1.0, 1.5, 2.0, 3.5, 4.25];
    double[] y = [2.0, 7.5, 5.0, 1.0, 1.5, 0.0];
    
    WMeanAccumulator!(float, Summation.pairwise, AssumeWeights.primary) m0;
    m0.put(x);
    WMeanAccumulator!(float, Summation.pairwise, AssumeWeights.primary) m1;
    m1.put(y);
    m0.put(m1);
    m0.wmean.shouldApprox == 29.25 / 12;
}

// repeat test, assume weights sum to 1
version(mir_stat_test)
@safe pure nothrow
unittest
{
    import mir.math.sum: Summation;
    import mir.ndslice.slice: slicedField;
    import mir.ndslice.topology: iota, map, repeat;
    import mir.test: shouldApprox;

    WMeanAccumulator!(double, Summation.pairwise, AssumeWeights.primary) x;
    auto y = iota(6);
    auto w = repeat(1.0, 6).map!(a => a / 6.0).slicedField;
    x.put(y, w);
    x.wmean.shouldApprox == 15.0 / 6;
}

// repeat test, do not assume weights sum to 1
version(mir_stat_test)
@safe pure nothrow
unittest
{
    import mir.math.sum: Summation;
    import mir.ndslice.slice: slicedField;
    import mir.ndslice.topology: iota, repeat;
    import mir.test: shouldApprox;

    WMeanAccumulator!(double, Summation.pairwise, AssumeWeights.primary) x;
    auto y = iota(6);
    auto w = repeat(1.0, 6).slicedField;
    x.put(y, w);
    x.wmean.shouldApprox == 15.0 / 6;
}

// range test without shape, assume weights sum to 1
version(mir_stat_test)
@safe pure nothrow
unittest
{
    import mir.math.sum: Summation;
    import std.algorithm: map;
    import std.range: iota;
    import mir.test: shouldApprox;

    WMeanAccumulator!(double, Summation.pairwise, AssumeWeights.sumToOne) x;
    auto y = iota(6);
    auto w = iota(1, 7).map!(a => a / 21.0);
    x.put(y, w);
    x.wmean.shouldApprox == 70.0 / 21;
}

// range test without shape, do not assume weights sum to 1
version(mir_stat_test)
@safe pure nothrow
unittest
{
    import mir.math.sum: Summation;
    import mir.test: shouldApprox;
    import std.range: iota;

    WMeanAccumulator!(double, Summation.pairwise, AssumeWeights.primary) x;
    auto y = iota(6);
    auto w = iota(1, 7);
    x.put(y, w);
    x.wmean.shouldApprox == 70.0 / 21;
}

// complex test, do not assume weights sum to 1
version(mir_stat_test)
@safe pure nothrow
unittest
{
    import mir.complex;
    import mir.math.sum: Summation;
    import mir.test: should;
    alias C = Complex!double;

    WMeanAccumulator!(C, Summation.pairwise, AssumeWeights.primary, double) x;
    x.put([C(1, 3), C(2), C(3)]);
    x.wmean.should == C(2, 1);
}

/++
Computes the weighted mean of the input.

By default, if `F` is not floating point type or complex type, then the result
will have a `double` type if `F` is implicitly convertible to a floating point 
type or a type for which `isComplex!F` is true.

Params:
    F = controls type of output
    summation = algorithm for calculating sums (default: Summation.appropriate)
    assumeWeights = true if weights are assumed to add to 1 (default = AssumeWeights.primary)
    G = controls the type of weights
Returns:
    The weighted mean of all the elements in the input, must be floating point or complex type

See_also: 
    $(MATHREF sum, Summation),
    $(MATHREF stat, mean),
    $(MATHREF stat, meanType)
+/
template wmean(F, Summation summation = Summation.appropriate,
               AssumeWeights assumeWeights = AssumeWeights.primary, 
               G = F, Summation weightsSummation = Summation.appropriate)
    if (!is(F : AssumeWeights))
{
    import mir.math.common: fmamath;
    import mir.math.stat: meanType;
    import mir.ndslice.slice: isConvertibleToSlice;
    import std.traits: isIterable;

    /++
    Params:
        s = slice-like
        w = weights
    +/
    @fmamath meanType!F wmean(SliceA, SliceB)(SliceA s, SliceB w)
        if (isConvertibleToSlice!SliceA && isConvertibleToSlice!SliceB)
    {
        import core.lifetime: move;

        alias H = typeof(return);
        WMeanAccumulator!(H, ResolveSummationType!(summation, SliceA, H), assumeWeights, 
                          G, ResolveSummationType!(weightsSummation, SliceB, G)) wmean;
        wmean.put(s.move, w.move);
        return wmean.wmean;
    }

    /++
    Params:
        r = range, must be finite iterable
    +/
    @fmamath meanType!F wmean(Range)(Range r)
        if (isIterable!Range)
    {
        import core.lifetime: move;

        alias H = typeof(return);
        WMeanAccumulator!(H, ResolveSummationType!(summation, Range, H), assumeWeights, G, ResolveSummationType!(weightsSummation, Range, G)) wmean;
        wmean.put(r.move);
        return wmean.wmean;
    }
}

/// ditto
template wmean(Summation summation = Summation.appropriate,
               AssumeWeights assumeWeights = AssumeWeights.primary,
               Summation weightsSummation = Summation.appropriate)
{
    import mir.math.common: fmamath;
    import mir.math.stat: meanType;
    import mir.ndslice.slice: isConvertibleToSlice;
    import std.traits: isIterable;

    /++
    Params:
        s = slice-like
        w = weights
    +/
    @fmamath meanType!SliceA wmean(SliceA, SliceB)(SliceA s, SliceB w)
        if (isConvertibleToSlice!SliceA && isConvertibleToSlice!SliceB)
    {
        import core.lifetime: move;
        import mir.math.sum: sumType;

        alias F = typeof(return);
        return .wmean!(F, summation, assumeWeights, sumType!SliceB, weightsSummation)(s.move, w.move);
    }

    /++
    Params:
        r = range, must be finite iterable
    +/
    @fmamath meanType!Range wmean(Range)(Range r)
        if (isIterable!Range)
    {
        import core.lifetime: move;

        alias F = typeof(return);
        return .wmean!(F, summation, assumeWeights, F, weightsSummation)(r.move);
    }
}

/// ditto
template wmean(F, AssumeWeights assumeWeights, Summation summation = Summation.appropriate, 
               G = F, Summation weightsSummation = Summation.appropriate)
    if (!is(F : AssumeWeights))
{
    import mir.math.common: fmamath;
    import mir.math.stat: meanType;
    import mir.ndslice.slice: isConvertibleToSlice;
    import std.traits: isIterable;

    /++
    Params:
        s = slice-like
        w = weights
    +/
    @fmamath meanType!F wmean(SliceA, SliceB)(SliceA s, SliceB w)
        if (isConvertibleToSlice!SliceA && isConvertibleToSlice!SliceB)
    {
        import core.lifetime: move;
        import mir.math.sum: sumType;

        alias H = typeof(return);
        return .wmean!(H, summation, assumeWeights, G, weightsSummation)(s.move, w.move);
    }

    /++
    Params:
        r = range, must be finite iterable
    +/
    @fmamath meanType!Range wmean(Range)(Range r)
        if (isIterable!Range)
    {
        import core.lifetime: move;

        alias F = typeof(return);
        return .wmean!(F, summation, assumeWeights, G, weightsSummation)(r.move);
    }
}

/// ditto
template wmean(F, bool assumeWeights, string summation = "appropriate", 
               G = F, string weightsSummation = "appropriate")
    if (!is(F : AssumeWeights))
{
    mixin("alias wmean = .wmean!(F, Summation." ~ summation ~ ", cast(AssumeWeights) assumeWeights, G, Summation." ~ weightsSummation ~ ");");
}

/// ditto
template wmean(bool assumeWeights, string summation = "appropriate",
               string weightsSummation = "appropriate")
{
    mixin("alias wmean = .wmean!(Summation." ~ summation ~ ", cast(AssumeWeights) assumeWeights, Summation." ~ weightsSummation ~ ");");
}

/// ditto
template wmean(F, string summation, bool assumeWeights = false,
               G = F, string weightsSummation = "appropriate")
    if (!is(F : AssumeWeights))
{
    mixin("alias wmean = .wmean!(F, Summation." ~ summation ~ ", cast(AssumeWeights) assumeWeights, G, Summation." ~ weightsSummation ~ ");");
}

/// ditto
template wmean(string summation, bool assumeWeights = false,
               string weightsSummation = "appropriate")
{
    mixin("alias wmean = .wmean!(Summation." ~ summation ~ ", cast(AssumeWeights) assumeWeights, Summation." ~ weightsSummation ~ ");");
}

/// ditto
template wmean(F, string summation, G, string weightsSummation, bool assumeWeights)
    if (!is(F : AssumeWeights))
{
    mixin("alias wmean = .wmean!(F, Summation." ~ summation ~ ", cast(AssumeWeights) assumeWeights, G, Summation." ~ weightsSummation ~ ");");
}

/// ditto
template wmean(string summation, string weightsSummation, bool assumeWeights = false)
{
    mixin("alias wmean = .wmean!(Summation." ~ summation ~ ", cast(AssumeWeights) assumeWeights, Summation." ~ weightsSummation ~ ");");
}

///
version(mir_stat_test)
@safe pure nothrow
unittest
{
    import mir.complex;
    import mir.ndslice.slice: sliced;
    import mir.test: should, shouldApprox;
    alias C = Complex!double;

    wmean([1.0, 2, 3], [1, 2, 3]).shouldApprox == (1.0 + 4.0 + 9.0) / 6;
    wmean!true([1.0, 2, 3], [1.0 / 6, 2.0 / 6, 3.0 / 6]).shouldApprox == (1.0 + 4.0 + 9.0) / 6;
    wmean([C(1, 3), C(2), C(3)], [1, 2, 3]).should == C(14.0 / 6, 3.0 / 6);

    wmean!float([0, 1, 2, 3, 4, 5].sliced(3, 2), [1, 2, 3, 4, 5, 6].sliced(3, 2)).shouldApprox == 70.0 / 21;

    static assert(is(typeof(wmean!float([1, 2, 3], [1, 2, 3])) == float));
}

/// If weights are not provided, then behaves like mean
version(mir_stat_test)
@safe pure nothrow
unittest
{
    import mir.complex;
    import mir.ndslice.slice: sliced;
    import mir.test: should;
    alias C = Complex!double;

    wmean([1.0, 2, 3]).should == 2;
    wmean([C(1, 3), C(2), C(3)]).should == C(2, 1);

    wmean!float([0, 1, 2, 3, 4, 5].sliced(3, 2)).should == 2.5;

    static assert(is(typeof(wmean!float([1, 2, 3])) == float));
}

/// Weighted mean of vector
version(mir_stat_test)
@safe pure nothrow
unittest
{
    import mir.ndslice.slice: sliced;
    import mir.ndslice.topology: iota, map;
    import mir.test: shouldApprox;

    auto x = [0.0, 1.0, 1.5, 2.0, 3.5, 4.25,
              2.0, 7.5, 5.0, 1.0, 1.5, 0.0].sliced;
    auto w = iota([12], 1);
    auto w_SumToOne = w.map!(a => a / 78.0);

    x.wmean.shouldApprox == 29.25 / 12;
    x.wmean(w).shouldApprox == 203.0 / 78;
    x.wmean!true(w_SumToOne).shouldApprox == 203.0 / 78;
}

/// Weighted mean of matrix
version(mir_stat_test)
@safe pure
unittest
{
    import mir.ndslice.fuse: fuse;
    import mir.ndslice.topology: iota, map;
    import mir.test: shouldApprox;

    auto x = [
        [0.0, 1.0, 1.5, 2.0, 3.5, 4.25],
        [2.0, 7.5, 5.0, 1.0, 1.5, 0.0]
    ].fuse;
    auto w = iota([2, 6], 1);
    auto w_SumToOne = w.map!(a => a / 78.0);

    x.wmean.shouldApprox == 29.25 / 12;
    x.wmean(w).shouldApprox == 203.0 / 78;
    x.wmean!true(w_SumToOne).shouldApprox == 203.0 / 78;
}

/// Column mean of matrix
version(mir_stat_test)
@safe pure
unittest
{
    import mir.algorithm.iteration: all;
    import mir.math.common: approxEqual;
    import mir.ndslice.fuse: fuse;
    import mir.ndslice.topology: alongDim, byDim, iota, map, universal;

    auto x = [
        [0.0, 1.0, 1.5, 2.0, 3.5, 4.25],
        [2.0, 7.5, 5.0, 1.0, 1.5, 0.0]
    ].fuse;
    auto w = iota([2], 1).universal;
    auto result = [4.0 / 3, 16.0 / 3, 11.5 / 3, 4.0 / 3, 6.5 / 3, 4.25 / 3];

    // Use byDim or alongDim with map to compute mean of row/column.
    assert(x.byDim!1.map!(a => a.wmean(w)).all!approxEqual(result));
    assert(x.alongDim!0.map!(a => a.wmean(w)).all!approxEqual(result));

    // FIXME
    // Without using map, computes the mean of the whole slice
    // assert(x.byDim!1.wmean(w) == x.sliced.wmean);
    // assert(x.alongDim!0.wmean(w) == x.sliced.wmean);
}

/// Can also set algorithm or output type
version(mir_stat_test)
@safe pure nothrow
unittest
{
    import mir.ndslice.slice: sliced;
    import mir.ndslice.topology: repeat, universal;
    import mir.test: shouldApprox;

    //Set sum algorithm (also for weights) or output type

    auto a = [1, 1e100, 1, -1e100].sliced;

    auto x = a * 10_000;
    auto w1 = [1, 1, 1, 1].sliced;
    auto w2 = [0.25, 0.25, 0.25, 0.25].sliced;

    x.wmean!"kbn"(w1).shouldApprox == 20_000 / 4;
    x.wmean!(true, "kbn")(w2).shouldApprox == 20_000 / 4;
    x.wmean!("kbn", true)(w2).shouldApprox == 20_000 / 4;
    x.wmean!("kbn", true, "pairwise")(w2).shouldApprox == 20_000 / 4;
    x.wmean!(true, "kbn", "pairwise")(w2).shouldApprox == 20_000 / 4;
    x.wmean!"kb2"(w1).shouldApprox == 20_000 / 4;
    x.wmean!"precise"(w1).shouldApprox == 20_000 / 4;
    x.wmean!(double, "precise")(w1).shouldApprox == 20_000.0 / 4;

    auto y = uint.max.repeat(3);
    y.wmean!ulong([1, 1, 1].sliced.universal).shouldApprox == 12884901885 / 3;
}

/++
For integral slices, can pass output type as template parameter to ensure output
type is correct.
+/
version(mir_stat_test)
@safe pure nothrow
unittest
{
    import mir.math.common: approxEqual;
    import mir.ndslice.slice: sliced;
    import mir.test: shouldApprox;

    auto x = [0, 1, 1, 2, 4, 4,
              2, 7, 5, 1, 2, 0].sliced;
    auto w = [1, 2, 3,  4,  5,  6,
              7, 8, 9, 10, 11, 12].sliced;

    auto y = x.wmean(w);
    y.shouldApprox(1.0e-10) == 204.0 / 78;
    static assert(is(typeof(y) == double));

    x.wmean!float(w).shouldApprox(1.0e-10) == 204f / 78;
}

/++
Mean works for complex numbers and other user-defined types (provided they
can be converted to a floating point or complex type)
+/
version(mir_test_weighted)
@safe pure nothrow
unittest
{
    import mir.complex;
    import mir.ndslice.slice: sliced;
    import mir.test: should;
    alias C = Complex!double;

    auto x = [C(1.0, 2), C(2, 3), C(3, 4), C(4, 5)].sliced;
    auto w = [1, 2, 3, 4].sliced;
    x.wmean(w).should == C(3, 4);
}

/// Compute weighted mean tensors along specified dimention of tensors
version(mir_stat_test)
@safe pure
unittest
{
    import mir.ndslice.fuse: fuse;
    import mir.ndslice.slice: sliced;
    import mir.ndslice.topology: alongDim, as, iota, map, universal;
    /++
      [[0,1,2],
       [3,4,5]]
     +/
    auto x = [
        [0, 1, 2],
        [3, 4, 5]
    ].fuse.as!double;
    auto w = [
        [1, 2, 3],
        [4, 5, 6]
    ].fuse;
    auto w1 = [1, 2].sliced.universal;
    auto w2 = [1, 2, 3].sliced;

    assert(x.wmean(w) == (70.0 / 21));

    auto m0 = [(0.0 + 6.0) / 3, (1.0 + 8.0) / 3, (2.0 + 10.0) / 3];
    assert(x.alongDim!0.map!(a => a.wmean(w1)) == m0);
    assert(x.alongDim!(-2).map!(a => a.wmean(w1)) == m0);

    auto m1 = [(0.0 + 2.0 + 6.0) / 6, (3.0 + 8.0 + 15.0) / 6];
    assert(x.alongDim!1.map!(a => a.wmean(w2)) == m1);
    assert(x.alongDim!(-1).map!(a => a.wmean(w2)) == m1);

    assert(iota(2, 3, 4, 5).as!double.alongDim!0.map!wmean == iota([3, 4, 5], 3 * 4 * 5 / 2));
}

// test chaining
version(mir_stat_test)
@safe pure nothrow
unittest
{
    import mir.test: should;
    [1.0, 2, 3, 4].wmean.should == 2.5;
}

// additional alongDim tests
version(mir_stat_test)
@safe pure nothrow
unittest
{
    import mir.algorithm.iteration: all;
    import mir.math.stat: meanType;
    import mir.ndslice.topology: iota, alongDim, map;

    auto x = iota([2, 2], 1);
    auto w = iota([2], 2);
    auto y = x.alongDim!1.map!(a => a.wmean(w));
    static assert(is(meanType!(typeof(y)) == double));
}

// @nogc test
version(mir_stat_test)
@safe pure nothrow @nogc
unittest
{
    import mir.ndslice.slice: sliced;
    import mir.test: shouldApprox;

    static immutable x = [0.0, 1.0, 1.5, 2.0, 3.5, 4.25,
                          2.0, 7.5, 5.0, 1.0, 1.5, 0.0];
    static immutable w = [1.0, 2, 3,  4,  5,  6,
                            7, 8, 9, 10, 11, 12];

    x.wmean.shouldApprox == 29.25 / 12;
    x.wmean(w).shouldApprox == 203.0 / 78;
}

/++
Output range for wsum.
+/
struct WSummator(T, Summation summation, U = T)
{
    import mir.ndslice.slice: isConvertibleToSlice, isSlice, kindOf;
    import std.range.primitives: isInputRange;
    import std.traits: isIterable;

    ///
    Summator!(T, summation) wsummator;

    ///
    F wsum(F = T)() const @safe @property pure nothrow @nogc
    {
        return cast(F) wsummator.sum;
    }

    ///
    void put(Slice1, Slice2)(Slice1 s, Slice2 w)
        if (isSlice!Slice1 && isSlice!Slice2)
    {
        static assert (Slice1.N == Slice2.N, "s and w must have the same number of dimensions");
        static assert (kindOf!Slice1 == kindOf!Slice2, "s and w must have the same kind");

        import mir.ndslice.slice: Contiguous;
        import mir.ndslice.topology: zip, map;

        assert(s._lengths == w._lengths, "WMeanAcumulator.put: both slices must have the same lengths");

        static if (kindOf!Slice1 != Contiguous && Slice1.N > 1) {
            assert(s.strides == w.strides, "WMeanAccumulator.put: cannot put canonical and universal slices when strides do not match");
            auto combine = s.zip!true(w);
        } else {
            auto combine = s.zip!false(w);
        }

        auto combine2 = combine.map!"a * b";
        wsummator.put(combine2);
    }

    ///
    void put(SliceLike1, SliceLike2)(SliceLike1 s, SliceLike2 w)
        if (isConvertibleToSlice!SliceLike1 && !isSlice!SliceLike1 &&
            isConvertibleToSlice!SliceLike2 && !isSlice!SliceLike2)
    {
        import mir.ndslice.slice: toSlice;
        this.put(s.toSlice, w.toSlice);
    }

    ///
    void put(Range)(Range r)
        if (isIterable!Range)
    {
        wsummator.put(r);
    }

    ///
    void put(RangeA, RangeB)(RangeA r, RangeB w)
        if (isInputRange!RangeA && !isConvertibleToSlice!RangeA &&
            isInputRange!RangeB && !isConvertibleToSlice!RangeB)
    {
        do
        {
            assert(!(!r.empty && w.empty) && !(r.empty && !w.empty),
                   "r and w must both be empty at the same time, one cannot be empty while the other has remaining items");
            this.put(r.front, w.front);
            r.popFront;
            w.popFront;
        } while(!r.empty || !w.empty); // Using an || instead of && so that the loop does not end early. mis-matched lengths of r and w sould be caught by above assert
    }

    ///
    void put()(T x, U w)
    {
        wsummator.put(x * w);
    }

    ///
    void put()(T x)
    {
        wsummator.put(x);
    }

    ///
    void put(F = T, G = U)(WSummator!(F, summation, G) wm)
    {
        wsummator.put(cast(T) wm.wsummator);
    }
}

///
version(mir_stat_test)
@safe pure nothrow
unittest
{
    import mir.math.sum: Summation;
    import mir.ndslice.slice: sliced;
    import mir.test: should;

    WSummator!(double, Summation.pairwise) x;
    x.put([0.0, 1, 2, 3, 4].sliced, [1, 2, 3, 4, 5].sliced);
    x.wsum.should == 40;
    x.put(5, 6);
    x.wsum.should == 70;
}

// dynamic array test
version(mir_stat_test)
@safe pure nothrow
unittest
{
    import mir.math.sum: Summation;
    import mir.test: should;

    WSummator!(double, Summation.pairwise) x;
    x.put([0.0, 1, 2, 3, 4], [1, 2, 3, 4, 5]);
    x.wsum.should == 40;
}

// static array test
version(mir_stat_test)
@safe pure nothrow @nogc
unittest
{
    import mir.math.sum: Summation;
    import mir.test: should;

    WSummator!(double, Summation.pairwise) x;
    static immutable y = [0.0, 1, 2, 3, 4];
    static immutable w = [1, 2, 3, 4, 5];
    x.put(y, w);
    x.wsum.should == 40;
}

// 2-d slice test
version(mir_stat_test)
@safe pure
unittest
{
    import mir.math.sum: Summation;
    import mir.ndslice.fuse: fuse;
    import mir.test: should;

    WSummator!(double, Summation.pairwise) x;
    auto y = [
        [0.0, 1, 2],
        [3.0, 4, 5]
    ].fuse;
    auto w = [
        [1.0, 2, 3],
        [4.0, 5, 6]
    ].fuse;
    x.put(y, w);
    x.wsum.should == 70;
}

// universal slice test
version(mir_stat_test)
@safe pure nothrow
unittest
{
    import mir.math.sum: Summation;
    import mir.ndslice.topology: iota, universal;
    import mir.test: should;

    WSummator!(double, Summation.pairwise) x;
    auto y = iota(6).universal;
    auto w = iota([6], 1).universal;
    x.put(y, w);
    x.wsum.should == 70;
}

// canonical slice test
version(mir_stat_test)
@safe pure nothrow
unittest
{
    import mir.math.sum: Summation;
    import mir.ndslice.topology: canonical, iota;
    import mir.test: should;

    WSummator!(double, Summation.pairwise) x;
    auto y = iota(6).canonical;
    auto w = iota([6], 1).canonical;
    x.put(y, w);
    x.wsum.should == 70;
}

// 2-d universal slice test
version(mir_stat_test)
@safe pure nothrow
unittest
{
    import mir.math.sum: Summation;
    import mir.ndslice.topology: iota, universal;
    import mir.test: should;

    WSummator!(double, Summation.pairwise) x;
    auto y = iota([2, 3]).universal;
    auto w = iota([2, 3], 1).universal;
    x.put(y, w);
    x.wsum.should == 70;
}

// 2-d canonical slice test
version(mir_stat_test)
@safe pure nothrow
unittest
{
    import mir.math.sum: Summation;
    import mir.ndslice.topology: canonical, iota;
    import mir.test: should;

    WSummator!(double, Summation.pairwise) x;
    auto y = iota([2, 3]).canonical;
    auto w = iota([2, 3], 1).canonical;
    x.put(y, w);
    x.wsum.should == 70;
}

/// Assume no weights, like Summator
version(mir_stat_test)
@safe pure nothrow
unittest
{
    import mir.math.sum: Summation;
    import mir.ndslice.slice: sliced;
    import mir.test: should;

    WSummator!(double, Summation.pairwise) x;
    x.put([0.0, 1, 2, 3, 4].sliced);
    x.wsum.should == 10;
    x.put(5);
    x.wsum.should == 15;
}

// dynamic array test, assume no weights
version(mir_stat_test)
@safe pure nothrow
unittest
{
    import mir.math.sum: Summation;
    import mir.test: should;

    WSummator!(double, Summation.pairwise) x;
    x.put([0.0, 1, 2, 3, 4]);
    x.wsum.should == 10;
}

// static array test, assume no weights
version(mir_stat_test)
@safe pure nothrow @nogc
unittest
{
    import mir.math.sum: Summation;
    import mir.test: should;

    WSummator!(double, Summation.pairwise) x;
    static immutable y = [0.0, 1, 2, 3, 4];
    x.put(y);
    x.wsum.should == 10;
}

// Adding WSummators
version(mir_stat_test)
@safe pure nothrow
unittest
{
    import mir.math.sum: Summation;
    import mir.test: should;

    double[] x = [0.0, 1.0, 1.5, 2.0, 3.5, 4.25];
    double[] y = [2.0, 7.5, 5.0, 1.0, 1.5, 0.0];

    WSummator!(float, Summation.pairwise) m0;
    m0.put(x);
    WSummator!(float, Summation.pairwise) m1;
    m1.put(y);
    m0.put(m1);
    m0.wsum.should == 29.25;
}

// repeat test
version(mir_stat_test)
@safe pure nothrow
unittest
{
    import mir.math.sum: Summation;
    import mir.ndslice.slice: slicedField;
    import mir.ndslice.topology: iota, repeat;
    import mir.test: should;

    WSummator!(double, Summation.pairwise) x;
    auto y = iota(6);
    auto w = repeat(1.0, 6).slicedField;
    x.put(y, w);
    x.wsum.should == 15;
}

// range test without shape
version(mir_stat_test)
@safe pure nothrow
unittest
{
    import mir.math.sum: Summation;
    import mir.test: should;
    import std.range: iota;

    WSummator!(double, Summation.pairwise) x;
    auto y = iota(6);
    auto w = iota(1, 7);
    x.put(y, w);
    x.wsum.should == 70;
}

// complex test
version(mir_stat_test)
@safe pure nothrow
unittest
{
    import mir.math.sum: Summation;
    import mir.complex;
    import mir.test: should;
    alias C = Complex!double;

    WSummator!(C, Summation.pairwise, double) x;
    x.put([C(1, 3), C(2), C(3)]);
    x.wsum.should == C(6, 3);
}

/++
Computes the weighted sum of the input.

Params:
    F = controls type of output
    summation = algorithm for calculating sums (default: Summation.appropriate)
    G = controls the type of weights
Returns:
    The weighted sum of all the elements in the input

See_also: 
    $(MATHREF sum, Summation)
+/
template wsum(F, Summation summation = Summation.appropriate, G = F)
{
    import mir.math.common: fmamath;
    import mir.math.sum: sumType;
    import mir.ndslice.slice: isConvertibleToSlice;
    import std.traits: isIterable;

    /++
    Params:
        s = slice-like
        w = weights
    +/
    @fmamath sumType!F wsum(SliceA, SliceB)(SliceA s, SliceB w)
        if (isConvertibleToSlice!SliceA && isConvertibleToSlice!SliceB)
    {
        import core.lifetime: move;

        alias H = typeof(return);
        WSummator!(H, ResolveSummationType!(summation, SliceA, H), G) wsum;
        wsum.put(s.move, w.move);
        return wsum.wsum;
    }

    /++
    Params:
        r = range, must be finite iterable
    +/
    @fmamath sumType!F wsum(Range)(Range r)
        if (isIterable!Range)
    {
        import core.lifetime: move;

        alias H = typeof(return);
        WSummator!(H, ResolveSummationType!(summation, Range, H), G) wsum;
        wsum.put(r.move);
        return wsum.wsum;
    }
}

/// ditto
template wsum(Summation summation = Summation.appropriate)
{
    import mir.math.common: fmamath;
    import mir.math.sum: sumType;
    import mir.ndslice.slice: isConvertibleToSlice;
    import std.traits: isIterable;

    /++
    Params:
        s = slice-like
        w = weights
    +/
    @fmamath sumType!SliceA wsum(SliceA, SliceB)(SliceA s, SliceB w)
        if (isConvertibleToSlice!SliceA && isConvertibleToSlice!SliceB)
    {
        import core.lifetime: move;
        import mir.math.sum: sumType;

        alias F = typeof(return);
        return .wsum!(F, summation, sumType!SliceB)(s.move, w.move);
    }

    /++
    Params:
        r = range, must be finite iterable
    +/
    @fmamath sumType!Range wsum(Range)(Range r)
        if (isIterable!Range)
    {
        import core.lifetime: move;

        alias F = typeof(return);
        return .wsum!(F, summation, F)(r.move);
    }
}

/// ditto
template wsum(F, string summation, G = F)
{
    mixin("alias wsum = .wsum!(F, Summation." ~ summation ~ ", G);");
}

/// ditto
template wsum(string summation)
{
    mixin("alias wsum = .wsum!(Summation." ~ summation ~ ");");
}

///
version(mir_stat_test)
@safe pure nothrow
unittest
{
    import mir.complex;
    import mir.math.common: approxEqual;
    import mir.ndslice.slice: sliced;
    import mir.test: should;
    alias C = Complex!double;

    wsum([1, 2, 3], [1, 2, 3]).should == (1 + 4 + 9);
    wsum([C(1, 3), C(2), C(3)], [1, 2, 3]).should == C((1 + 4 + 9), 3);

    wsum!float([0, 1, 2, 3, 4, 5].sliced(3, 2), [1, 2, 3, 4, 5, 6].sliced(3, 2)).should == 70;

    static assert(is(typeof(wmean!float([1, 2, 3], [1, 2, 3])) == float));
}

/// If weights are not provided, then behaves like sum
version(mir_stat_test)
@safe pure nothrow
unittest
{
    import mir.complex;
    import mir.ndslice.slice: sliced;
    import mir.test: should;
    alias C = Complex!double;

    wsum([1.0, 2, 3]).should == 6;
    wsum([C(1, 3), C(2), C(3)]).should == C(6, 3);

    wsum!float([0, 1, 2, 3, 4, 5].sliced(3, 2)).should == 15;

    static assert(is(typeof(wsum!float([1, 2, 3])) == float));
}

/// Weighted sum of vector
version(mir_stat_test)
@safe pure nothrow
unittest
{
    import mir.ndslice.slice: sliced;
    import mir.ndslice.topology: iota, map;
    import mir.test: should;

    auto x = [0.0, 1.0, 1.5, 2.0, 3.5, 4.25,
              2.0, 7.5, 5.0, 1.0, 1.5, 0.0].sliced;
    auto w = iota([12], 1);

    x.wsum.should == 29.25;
    x.wsum(w).should == 203;
}

/// Weighted sum of matrix
version(mir_stat_test)
@safe pure
unittest
{
    import mir.ndslice.fuse: fuse;
    import mir.ndslice.topology: iota;
    import mir.test: should;

    auto x = [
        [0.0, 1.0, 1.5, 2.0, 3.5, 4.25],
        [2.0, 7.5, 5.0, 1.0, 1.5, 0.0]
    ].fuse;
    auto w = iota([2, 6], 1);

    x.wsum.should == 29.25;
    x.wsum(w).should == 203;
}

/// Column sum of matrix
version(mir_stat_test)
@safe pure
unittest
{
    import mir.algorithm.iteration: all;
    import mir.math.common: approxEqual;
    import mir.ndslice.fuse: fuse;
    import mir.ndslice.topology: alongDim, byDim, iota, map, universal;

    auto x = [
        [0.0, 1.0, 1.5, 2.0, 3.5, 4.25],
        [2.0, 7.5, 5.0, 1.0, 1.5, 0.0]
    ].fuse;
    auto w = iota([2], 1).universal;
    auto result = [4, 16, 11.5, 4, 6.5, 4.25];

    // Use byDim or alongDim with map to compute sum of row/column.
    assert(x.byDim!1.map!(a => a.wsum(w)).all!approxEqual(result));
    assert(x.alongDim!0.map!(a => a.wsum(w)).all!approxEqual(result));

    // FIXME
    // Without using map, computes the sum of the whole slice
    // assert(x.byDim!1.wsum(w) == x.sliced.wsum);
    // assert(x.alongDim!0.wsum(w) == x.sliced.wsum);
}

/// Can also set algorithm or output type
version(mir_stat_test)
@safe pure nothrow
unittest
{
    import mir.ndslice.slice: sliced;
    import mir.ndslice.topology: repeat, universal;
    import mir.test: should;

    //Set sum algorithm (also for weights) or output type

    auto a = [1, 1e100, 1, -1e100].sliced;

    auto x = a * 10_000;
    auto w1 = [1, 1, 1, 1].sliced;
    auto w2 = [0.25, 0.25, 0.25, 0.25].sliced;

    x.wsum!"kbn"(w1).should == 20_000;
    x.wsum!"kbn"(w2).should == 20_000 / 4;
    x.wsum!"kb2"(w1).should == 20_000;
    x.wsum!"precise"(w1).should == 20_000;
    x.wsum!(double, "precise")(w1).should == 20_000;

    auto y = uint.max.repeat(3);
    y.wsum!ulong([1, 1, 1].sliced.universal).should == 12884901885;
}

/++
wsum works for complex numbers and other user-defined types
+/
version(mir_test_weighted)
@safe pure nothrow
unittest
{
    import mir.complex;
    import mir.ndslice.slice: sliced;
    import mir.test: should;
    alias C = Complex!double;

    auto x = [C(1.0, 2), C(2, 3), C(3, 4), C(4, 5)].sliced;
    auto w = [1, 2, 3, 4].sliced;
    x.wsum(w).should == C(30, 40);
}

/// Compute weighted sum tensors along specified dimention of tensors
version(mir_stat_test)
@safe pure
unittest
{
    import mir.ndslice.fuse: fuse;
    import mir.ndslice.slice: sliced;
    import mir.ndslice.topology: alongDim, as, iota, map, universal;
    /++
      [[0,1,2],
       [3,4,5]]
     +/
    auto x = [
        [0, 1, 2],
        [3, 4, 5]
    ].fuse.as!double;
    auto w = [
        [1, 2, 3],
        [4, 5, 6]
    ].fuse;
    auto w1 = [1, 2].sliced.universal;
    auto w2 = [1, 2, 3].sliced;

    assert(x.wsum(w) == 70);

    auto m0 = [(0 + 6), (1 + 8), (2 + 10)];
    assert(x.alongDim!0.map!(a => a.wsum(w1)) == m0);
    assert(x.alongDim!(-2).map!(a => a.wsum(w1)) == m0);

    auto m1 = [(0 + 2 + 6), (3 + 8 + 15)];
    assert(x.alongDim!1.map!(a => a.wsum(w2)) == m1);
    assert(x.alongDim!(-1).map!(a => a.wsum(w2)) == m1);
}

// test chaining
version(mir_stat_test)
@safe pure nothrow
unittest
{
    import mir.test: should;
    [1.0, 2, 3, 4].wsum.should == 10;
}

// additional alongDim tests
version(mir_stat_test)
@safe pure nothrow
unittest
{
    import mir.math.sum: sumType;
    import mir.ndslice.topology: iota, alongDim, map;

    auto x = iota([2, 2], 1);
    auto w = iota([2], 2);
    auto y = x.alongDim!1.map!(a => a.wsum(w));
    static assert(is(sumType!(typeof(y)) == long));
}

// @nogc test
version(mir_stat_test)
@safe pure nothrow @nogc
unittest
{
    import mir.test: should;

    static immutable x = [0.0, 1.0, 1.5, 2.0, 3.5, 4.25,
                          2.0, 7.5, 5.0, 1.0, 1.5, 0.0];
    static immutable w = [1.0, 2, 3,  4,  5,  6,
                            7, 8, 9, 10, 11, 12];

    x.wsum.should == 29.25;
    x.wsum(w).should == 203.0;
}
