import betterr.random, betterr.vector;
import betterr.r;
import std.stdio;

void main() {
  startR();
	auto s = Sample([2.5, 5.0, 7.5, 10.0]);
	/* prob has not been set, so it has length 0 */
	writeln(s.prob.rows);
	s.size = 25;
	s.replace = true;
	writeln(s.draw());

	auto s2 = Sample([4.7], 25, true);
	writeln(s2.draw());

	// This is an error
	//~ auto s3 = Sample([4.7], 25, false);
	//~ writeln(s3.draw());
	
	/* This should yield 5.0 75% of the time */
	auto s4 = Sample([2.5, 5.0, 7.5, 10.0], 400, true, [0.1, 0.7, 0.1, 0.1]);
	writeln(s4.draw());

	/* Set the seed to make this reproducible */
	auto s5 = Sample([2.5, 5.0, 7.5, 10.0], 10, true);
	setSeed(100);
	writeln(s5.draw());
	setSeed(100);
	writeln(s5.draw());

	rnorm(100, 0.5, 0.1).print("Vector of 100 normal draws");
	runif(100, 0.5, 1.0).print("Vector of 100 uniform draws");
	rgamma(100, 0.5).print("Vector of 100 gamma draws");
	rbeta(100, 0.5, 2.4).print("Vector of 100 beta draws");
	rbinom(10, 8000, 0.4).print("Vector of 10 binomial draws, number of successes out of 8000 for each trial");
	rcauchy(100).print("Vector of 100 draws from a Cauchy distribution");
	rchisq(100, 20, 0.0).print("Vector of 100 draws from a chi-squared distribution");
	rexp(100, 2.5).print("Vector of 100 draws from an exponential distribution");
	rf(100, 10, 30).print("Vector of 100 draws from an F distribution");
	rf(100, 10, 30, 0.0).print("Vector of 100 draws from an F distribution");
  rgeom(100, Vector([0.5, 0.3, 0.1])).print("Vector of 100 draws from a geometric distribution");
  rgeom(100, [0.5, 0.3, 0.1]).print("Vector of 100 draws from a geometric distribution");
  rhyper(100, 40, 20, 6).print("Vector of 100 draws from a hypergeometric distribution");
  rhyper([147, 22, 67], 40, 20, 6).print("Vector of 3 draws from a hypergeometric distribution");
  rhyper([147, 22, 67], 40, 20, 6).print("Vector of 3 draws from a hypergeometric distribution");
  rlnorm(100).print("Vector of 100 draws from a log normal distribution");
  rmultinom(10, 20, [1.0, 3, 6, 10]).print("10 vectors drawn from a multinomial distribution");
  rnbinom(100, 20, 0.2, 0.0).print("100 draws from a negative binomial distribution");
  rnbinom(100, 20, double.nan, 12.0).print("100 draws from a negative binomial distribution");
  rnbinom([1, 2, 3, 4, 5], 20, 0.2, 0.0).print("5 draws from a negative binomial distribution");
  rnbinom([1, 2, 3, 4, 5], 20, 0.2).print("5 draws from a negative binomial distribution");
  rpois(100, 4).print("100 draws from a Poisson distribution");
  rt(100, 2, 0.0).print("100 draws from a t distribution");
  rt(100, 2).print("100 draws from a t distribution");
  rt([1, 2, 3, 4, 5], 2).print("5 draws from a t distribution");
  rweibull(100, 1.0, 1.0).print("100 draws from a Weibull distribution");
  
  // Test the generator
  // Print out 100 standard normal draws
  // But draw only 5 at a time
  Generator!("norm", 5) norm;
  foreach(_; 0..100) {
    write(norm.draw(), " ");
  }
  writeln();

  Generator!"norm" norm2;
  foreach(_; 0..100) {
    write(norm2.draw(), " ");
  }
  writeln();
  
  Generator!"unif" unif1;
  unif1.min = 2.5;
  unif1.max = 7.5;
  foreach(_; 0..10) {
    write(unif1.draw(), " ");
  }
  writeln();
  
  Generator!"binom" binom;
  binom.size = 8000;
  binom.prob = 0.4;
  foreach(_; 0..10) {
    write(binom.draw(), " ");
  }
  writeln();
  
  auto v = Vector(50);
  Generator!("custom", 8, double) customGen;
  customGen.cmd = "rnorm(1000, mean=-4.6, sd=2.1)";
  foreach(ii; 0..50) {
    v[ii] = customGen.draw()/3.0;
  }
  v.print("50 draws from a custom generator");

  /* Check that PRNG works */
  prngInit(1);
  writeln(rnorm(10));
  writeln();
  prngInit(7);
  writeln(rnorm(10));
  writeln();
  prngInit(1);
  writeln(rnorm(10));
  writeln();
  prngInit(7);
  writeln(rnorm(10));
  writeln();

	closeR();
}
