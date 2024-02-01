import betterr.dataframe, betterr.lm;
import betterr.r;
import std.conv, std.stdio;

void main() {
  startR();
  auto df = DataFile("macrodata.csv");
  df.commentChar = "#";
  auto dataset = df.read();
  LMConfig model = {
    rhs: ["ffr", "inflation", "unrate"], 
    lhs: "employment",
    intercept: true
  };
  auto fit = lm(dataset, model);
  fit.print("Regression output");
  fit.summary.print("Summary of regression output");
  fit.beta.print("Coefficients");
  writeln("Intercept of the regression: ", fit.beta[0]);
  /* Warning: Save model.beta in its own variable if you're going to
	 * heavily use the elements. Calling model.beta() requires a call to
	 * R and the creation of a new struct. Extreme inefficiency, but no
	 * big deal for a quick use like here. */
	writeln(fit.dfResidual);
	//~ fit.model.print("Model Matrix");
	writeln(fit.sigma);
	writeln(fit.rsq);
	writeln(fit.adjrsq);
	writeln(fit.fstat);
	fit.unscaledCov.print("Unscaled covariance matrix");
  fit.nwCov.print("Newey-West covariance matrix");
  fit.nwStdErrors.print("Newey-West standard errors");
	fit.coefficients.print("Coefficients and other information");
  fit.nwCoefficients.print("Newey-West adjustment");
  fit.whiteCoefficients.print("White adjustment");
  fit.whiteCov.print("White covariance matrix");
  closeR();
}
