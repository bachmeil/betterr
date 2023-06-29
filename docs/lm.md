# lm

Functionality for basic OLS regression. More will be added.

## struct LMFit

Holds the output of a regression.

- `Vector beta()`: Vector of coefficients
- `List coefficients()`: Coefficients plus standard errors and other information, as provided by `summary`
- `Vector residuals()`
-	`Vector fittedValues()`
- `int dfResidual()`
-	`List model()`
-	`double sigma()`
-	`double rsq()`
-	`double adjrsq()`
-	`double fstat()`	
-	`List unscaledCov()`
- `void print(string msg="")`: Print the output provided by `lm` (not `summary`), with an optional message
- `double pred(double x)`: Makes a prediction with new data `x`, assuming this is a simple linear regression (one regressor)
- `double pred(Vector x)`: Makes a prediction with new data `x`, for any number of regressors. Automatically handles the intercept, so you should not include a 1 or any other number in `x` to account for the intercept.

## Estimation

- `LMFit lm(T1, T2)(T1 y, T2 x)`
- `LMFit lm(T1, T2)(T1 y, T2 x, LMConfig conf)`
- `LMFit lm(DataFrame df, LMConfig conf)`

## struct LMConfig

Configuration for a regression.

- `string lhs`
- `string[] rhs`
- `bool intercept = true`
- `Subset subset`


