/* Lance Bachmeier
 * May 9, 2025
 * 
 * This module wraps the arima functionality in base R.
 * The goal is to add ARIMA functionality to D and then move on to other
 * things without going down rabbit holes. */
import betterr.everything;

struct ArimaFit(long f) {
  List fit;
  long p;
  long q;
  
  this(string n, long[2] order) {
    this(n, order[0], order[1]);
  }
  
  this(string n, long _p, long _q) {
    string cmd = i"arima($(n), order=c($(_p), 0, $(_q)))".text;
    fit = List(cmd);
    p = _p;
    q = _q;
  }
  
  Vector coef() {
    return Vector(i"$(fit.name)$coef".text);
  }
  
  TS!f pred(long h=1) {
    return TS!f(i"predict($(fit.name), $(h))$pred".text);
  }
}

ArimaFit!f arima(long f)(TS!f x, long[2] order) {
  return ArimaFit!f(x.name, order);
}

ArimaFit!f arima(long f)(TS!f x, long p, long q) {
  return ArimaFit!f(x.name, p, q);
}
