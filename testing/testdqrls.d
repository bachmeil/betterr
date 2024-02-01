import betterr.lm, betterr.random, betterr.matrix, betterr.vector;
import betterr.r;
import std.stdio;

void main() {
	startR();
	
	setSeed(100);
	auto y = Vector(rnorm(50));
	auto xtmp = rnorm(100);
	auto x = Matrix(xtmp, 50, 2);
	x.print("x");
	auto fit = dqrls(y, x);
	y.print("y");
	x.print("x");
	writeln("coefficients: ", fit.coef);
	writeln("residuals: ", fit.residuals);
	closeR();
}
