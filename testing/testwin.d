//pragma(lib, "libRInside.lib");
//pragma(lib, "R.dll");

//import betterr.baser, betterr.vector;
import betterr.rwin;
import std.stdio;
import core.runtime;

//
//startrptr startR;

//
//auto startR = GetProcAddress(handle, "startR");

void main() {
  writeln("Starting...");
  loadRFunctions("C:/Program Files/R/R-4.3.1/bin/x64",
    "C:/Users/bachm/AppData/Local/R/win-library/4.3/RInside/lib/x64");
  startR();
  // Error
  //evalInR(toUTFz!(char*)("3-'v'"));
  import std.utf: toUTFz;
  evalInR(toUTFz!(char*)("print(1.7)"));
  evalR("print(3)");
  printR(robj(4));
  closeR();
  writeln("Finished...");
}


  //writeln(6);
  //writeln(handle is null);
  //
  //writeln(startR);
  //startR();
  //printR(robj(4));
  //auto v = Vector([1.4, 1.7, 3.6, - 3.2, -2.6]);
  // v.print("Vector");
  // writeln("Absolute value of v: ", abs(v));
  
  // auto v2 = Vector([1.4, 1.7, 3.6]);
  // v2.print("Vector v2");
  // writeln("Square root of v2: ", sqrt(v2));
  // writeln("Absolute value of -3: ", abs(-3));
  // writeln("Absolute value of -3.7: ", abs(-3.7));
  // writeln("Square root of 4: ", sqrt(4));
  // writeln("Square root of 5: ", sqrt(5));
  // writeln("Square root of 5.0: ", sqrt(5.0));
  // writeln("Square root of 4.5: ", sqrt(4.5));
  
  // writeln("Ceiling of 4.2: ", ceiling(4.2));
  // writeln("Ceiling of -4.2: ", ceiling(-4.2));
  // writeln("Ceiling of 4.0: ", ceiling(4.0));
  // writeln("Ceiling of -4.0: ", ceiling(-4.0));
  // writeln("Ceiling of -4: ", ceiling(-4));
  // writeln("Ceiling of 4: ", ceiling(4));
  // writeln("Round 7.4: ", round(7.4));
  // writeln("Round 7.48723 to 3 places after the decimal: ", round(7.48723, "digits=3"));
  // writeln("Save 3 significant digits of 7.48723: ", signif(7.48723, "digits=3"));
  
  // writeln("atanh(0.1): ", atanh(0.1));
  // atanh(Vector([0.1, 0.2, 0.3])).print("atanh of a vector");
  // atanh([0.1, 0.2, 0.3]).print("atanh of a vector");
  
  // cumsum([1.2, 2.4, 3.6]).print("Cumulative sum");
  // cumprod([1.2, 2.4, 3.6]).print("Cumulative product");
  // cummax([1.2, 2.4, 3.6]).print("Cumulative max");
  // cummin([1.2, 2.4, 3.6]).print("Cumulative min");
  
  // log([1.2, 2.4, 3.6]).print("log");
  // writeln("log(3.5): ", log(3.5));
  // writeln("log(3.5, base=6): ", log(3.5, "base=6"));
  // writeln("exp(1): ", exp(1));
  
  // writeln("Product of the elements [2.1, 3, 5]: ", prod([2.1, 3, 5]));
  // range([0.5, 2.5, 1.5, 2.4]).print("Range of elements");
  // rank([0.5, 2.5, 1.5, 2.4]).print("Elements ranked");
  // rev([0.5, 2.5, 1.5, 2.4]).print("Elements reversed");
  // seq(0.5, 1.5, 0.02).print("A sequence");
  // sort([0.5, 2.5, 1.5, 2.4]).print("Sorted");
  // sort([0.5, 2.5, 1.5, 2.4], "decreasing = TRUE").print("Reverse sorted");
  // auto summ = summary([0.5, 2.5, 1.5, 2.4]);
  // writeln(summ["1st Qu."]);
  // writeln(summ);
  
  // /* Test the implementations based on Mir */
  // writeln(mean([1.1, 2.2, 3.3]));
  // writeln(mean([1.1, 2.2, double.nan, 3.3]));
  // writeln(mean(v2));
  
  // /* Also based on Mir, but remove missing values first */
  // writeln(mean([1.1, 2.2, double.nan, 3.3], false));
  // writeln(mean([1.1, 2.2, double.nan, 3.3], true));
  // Vector v3 = Vector([1.1, double.nan, 3.3, 4.4]);
  // writeln(mean(v3, false));
  // writeln(mean(v3, true));
  
  // /* Test the implementations based on Mir */
  // writeln("-- Test median --");
  // writeln(median([1.1, 2.2, 3.3]));
  // writeln(median([1.1, 2.2, double.nan, 3.3]));
  // writeln(median(v2));
  // writeln(median([1.1, 2.2, double.nan, 3.3], false));
  // writeln(median([1.1, 2.2, double.nan, 3.3], true));
  // writeln(median(v3, false));
  // writeln(median(v3, true));
  
  // /* Test the implementations based on Mir */
  // writeln("-- Test sum --");
  // writeln(sum([1.1, 2.2, 3.3]));
  // writeln(sum([1.1, 2.2, double.nan, 3.3]));
  // writeln(sum(v2));
  // writeln(sum([1.1, 2.2, double.nan, 3.3], false));
  // writeln(sum([1.1, 2.2, double.nan, 3.3], true));
  // writeln(sum(v3, false));
  // writeln(sum(v3, true));
  
  // /* Use Phobos */
  // writeln("-- Test max --");
  // writeln(max([1.1, 2.2, 3.3, 1.7]));
  // writeln(max([1.1, 2.2, 3.3, double.nan, 1.7]));
  // writeln(max(Vector([1.1, 2.2, 3.3, 1.7])));
  // writeln(max(Vector([1.1, 2.2, 3.3, double.nan, 1.7])));
  // writeln(max([1.1, 2.2, 3.3, double.nan, 1.7], true));
  // writeln(max([double.nan, double.nan], true));
  // writeln(max(Vector([1.1, 2.2, 3.3, double.nan, 1.7]), true));
  // writeln(max(Vector([double.nan, double.nan]), true));

  // writeln("-- Test min --");
  // writeln(min([1.1, 2.2, 0.2, 3.3, 1.7]));
  // writeln(min([1.1, 2.2, 0.2, 3.3, double.nan, 1.7]));
  // writeln(min(Vector([1.1, 2.2, 0.2, 3.3, 1.7])));
  // writeln(min(Vector([1.1, 2.2, 0.2, 3.3, double.nan, 1.7])));
  // writeln(min([1.1, 2.2, 0.2, 3.3, double.nan, 0.2, 1.7], true));
  // writeln(min([double.nan, double.nan], true));
  // writeln(min(Vector([1.1, 2.2, 0.2, 3.3, double.nan, 1.7]), true));
  // writeln(min(Vector([double.nan, double.nan]), true));
  
  // writeln("-- Standard deviation --");
  // writeln(sd([1.1, 2.2, 0.2, 3.3, 1.7]));
  // writeln(sd([1.1, 2.2, 0.2, 3.3, double.nan, 1.7]));
  // writeln(sd(Vector([1.1, 2.2, 0.2, 3.3, 1.7])));
  // writeln(sd(Vector([1.1, 2.2, 0.2, 3.3, double.nan, 1.7])));
  // writeln(sd([1.1, 2.2, 0.2, 3.3, double.nan, 1.7], false));
  // writeln(sd([1.1, 2.2, 0.2, 3.3, double.nan, 1.7], true));
  // writeln(sd(Vector([1.1, 2.2, 0.2, 3.3, double.nan, 1.7]), false));
  // writeln(sd(Vector([1.1, 2.2, 0.2, 3.3, double.nan, 1.7]), true));
  
  // writeln("-- Variance --");
  // writeln(var([1.1, 2.2, 0.2, 3.3, 1.7]));
  // writeln(var([1.1, 2.2, 0.2, 3.3, double.nan, 1.7]));
  // writeln(var(Vector([1.1, 2.2, 0.2, 3.3, 1.7])));
  // writeln(var(Vector([1.1, 2.2, 0.2, 3.3, double.nan, 1.7])));
  // writeln(var([1.1, 2.2, 0.2, 3.3, double.nan, 1.7], false));
  // writeln(var([1.1, 2.2, 0.2, 3.3, double.nan, 1.7], true));
  // writeln(var(Vector([1.1, 2.2, 0.2, 3.3, double.nan, 1.7]), false));
  // writeln(var(Vector([1.1, 2.2, 0.2, 3.3, double.nan, 1.7]), true));
  
  // writeln("-- Quantiles --");
  // writeln(quantile([1.1, 2.2, 3.3, 4.4, 5.5], 0.5));
  // writeln(quantile([1.1, 2.2, 3.3, 4.4, 5.5], 0.8));
  // writeln(quantile([1.1, 2.2, 3.3, 4.4, 5.5], 0.25));
  // writeln(quantile([1.1, 2.2, double.nan, 3.3, 4.4, 5.5], 0.25, false));
  // writeln(quantile([1.1, 2.2, 3.3, 3.3, 4.4, 5.5], 0.25, false));
  // writeln(quantile([1.1, 2.2, double.nan, 3.3, 4.4, 5.5], 0.25, true));
  // writeln(quantile(Vector([1.1, 2.2, double.nan, 3.3, 4.4, 5.5]), 0.25, true));
  // writeln(quantile(Vector([1.1, 2.2, 3.3, 4.4, 5.5]), 0.5));
  // writeln(quantile(Vector([1.1, 2.2, 3.3, 4.4, 5.5]), 0.8));
  // writeln(quantile(Vector([1.1, 2.2, 3.3, 4.4, 5.5]), 0.25));
  // writeln(quantile([1.1, 2.2, 3.3, 4.4, 5.5], [0.2, 0.5, 0.8]));
  // writeln(quantile([1.1, 2.2, 3.3, double.nan, 4.4, 5.5], [0.2, 0.5, 0.8], true));
  // writeln(quantile(Vector([1.1, 2.2, 3.3, 4.4, 5.5]), [0.2, 0.5, 0.8]));
  // writeln(quantile(Vector([1.1, 2.2, 3.3, double.nan, 4.4, 5.5]), [0.2, 0.5, 0.8], true));

  // //~ writeln(function4("log", "log", "base = exp(1)"));*/
  // closeR();
//}
