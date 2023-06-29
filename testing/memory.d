/* This is probably only of interest if you're familiar with the inner
 * workings of R. */
import betterr.matrix;
import betterr.r;
import std.stdio, std.utf;

void main() {
  startR();
  
  writeln("\nCreate a matrix and then create a copy of that matrix.");
	auto yval = evalR("y <- matrix(1:9, ncol=3); y");
	auto xval = evalR("x <- y; x");
	printR(yval);
	printR(xval);
	writeln("\nThey have the same pointer:");
	writeln(yval, " ", xval);
	writeln("\nChanging an element in one changes it in both:");
	INTEGER(xval)[0] = 10;
	printR(yval);
	printR(xval);
	writeln("\nCreate a second copy, but force the copy:");
	auto zval = evalR("z <- y[]");
	writeln("\nThey have the same values:");
	printR(yval);
	printR(zval);
	writeln("\nBut they have different pointers:");
	writeln(yval, " ", zval);
	writeln("\nChanging an element in one does not affect the other:");
	INTEGER(zval)[1] = 11;
	printR(yval);
	printR(zval);
	writeln("\nIt doesn't matter when using [], because there's going to be a copy made:");
	auto wval = evalR("w <- y[]");
	writeln(yval, " ", wval);
	writeln("\nWe don't have the problem with vector -> matrix copies");
	evalRQ("v <- 1:9; m <- as.matrix(v)");
	writeln(evalR("v"), " ", evalR("m"));
	printR(Rf_install(toUTFz!(char*)("tmp")));
	printR(3.robj);
	/* Add a variable to the current environment */
	addToR("tmp", 3.robj);
	printR("environment()");
	printR(evalR("tmp"));
	
	/* Create a matrix from D and insert into R */
	Robj temp;
	Rf_protect(temp = Rf_allocMatrix(14, 2, 4));
	REAL(temp)[0..8] = [1.1, 2.2, 3.3, 4.4,
											5.5, 6.6, 7.7, 8.8];
	printR(temp);
	addToR("temp", temp);
	printR("temp");
	writeln(evalR("temp"));
	writeln(temp);
	
	printR(evalR("temp"));
	printR(temp);
	Rf_unprotect_ptr(temp); // Eligible to be, but might not be, collected
	evalRQ("rm(temp)");
	printR(temp);
	writeln("-----------");
  // This will throw an error if you uncomment it.
	//printR("temp");
	
	closeR();
}
