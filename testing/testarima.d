import betterr.everything;
import betterr.random;
import rarima;

void main() {
  startR();

  auto ts = TS!12(rnorm(120), [1990, 4]);
  auto fit = arima(ts, 1, 1);
  fit.coef.print("Coefficients");
  fit.pred(12).print("Predictions for the next 12 months");

  closeR();
}
