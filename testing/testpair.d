/* I was playing around with a few things in this file. It's probably
 * not of interest to anyone else, but I have no reason to remove it
 * from the repo. */
import betterr.r;
import cran.matrix;
import std.conv, std.exception, std.stdio, std.typecons, std.utf;

void main() {
  startR();
  Robj argListPointer = Rf_protect(Rf_allocList(8));
  SET_TYPEOF(argListPointer, 6);
  Robj fillPointer = argListPointer;
  SETCAR(fillPointer, Rf_install(toUTFz!(char*)("Matrix")));
  fillPointer = CDR(fillPointer);
  SETCAR(fillPointer, Rf_allocVector(14,0));
  fillPointer = CDR(fillPointer);
  SETCAR(fillPointer, 3.robj);
  fillPointer = CDR(fillPointer);
  SETCAR(fillPointer, 2.robj);
  fillPointer = CDR(fillPointer);
  SETCAR(fillPointer, RFalse);
  fillPointer = CDR(fillPointer);
  SETCAR(fillPointer, RNil);
  fillPointer = CDR(fillPointer);
  SETCAR(fillPointer, RFalse);
  fillPointer = CDR(fillPointer);
  SETCAR(fillPointer, RFalse);
  printR(argListPointer);
  printR(Mmatrix(argListPointer));
  // Save to an Robj
  Robj mat1 = Rf_protect(Mmatrix(argListPointer));
  printR(mat1);
  // No longer need mat1, so unprotect it manually
  Rf_unprotect(1);
  // This unprotects automatically when it goes out of scope
  // As a unique pointer, you shouldn't have to worry about references
  // to it floating around.
  Unique!RPointer mat2 = new RPointer(Rf_protect(Mmatrix(argListPointer)));
  mat2.print("Unique pointer version");
  // Put it in the global environment so it's accessible from R
  // If you protect it, you'll want to unprotect it
  // It won't be GCed if it's active in R
  // It will never be GCed if it's protected
  Robj mat3 = Mmatrix(argListPointer);
  mat3.toR("matrix.version3");
  evalRQ(`cat("Printing this from R\n");print(matrix.version3)`);
  closeR();
}

struct RPointer {
	Robj obj;
  
  this(Robj x) {
    obj = x;
  }

	~this() {
		Rf_unprotect_ptr(obj);
	}
  
  void print(string msg="") {
    if (msg.length > 0) {
      writeln(msg ~ ":");
    }
    printR(obj);
  }
}
