/*
 * This module provides the following functionality:
 * 
 * - Basic GSL random number generation functions.
 * - The Mersenne-Twister generator (GSL's default) for generating random
 *   draws on a single processor.
 * - GSL's random number generators for specific distributions.
 * 
 * Other GSL generators are not supported at this time. For parallel
 * random number generation, use gsl.rng.parallel.
 * 
 * Since this uses ImportC, you need to add a command like -P-Igsl/rng gsl/rng/*.c
 * when compiling. That tells the D compiler where to find the header
 * and source files.
 */
module gsl.rng.random;
import distributions, gslrng, mt;
import core.stdc.config;

alias RNGType = gsl_rng_type*;
alias RNG = gsl_rng*;

RNGType MT19937() { return gsl_rng_mt19937; }
RNGType MT19937_1999() { return gsl_rng_mt19937_1999; }
RNGType MT19937_1998() { return gsl_rng_mt19937_1998; }

RNG allocRNG(RNGType type) { return gsl_rng_alloc(type); }
RNG allocRNG() { return gsl_rng_alloc(gsl_rng_mt19937); }
void freeRNG(RNG r) { gsl_rng_free(r); }
void setSeed(RNG r, c_ulong seed) { gsl_rng_set(r, seed); }

/* Using the R naming convention of `r` at the start of functions that
 * generate random draws. Same parameters as GSL unless otherwise noted.
 * d is for the pdf and p is for the cdf. */
double rbeta(RNG r, double a, double b) {
  return gsl_ran_beta(r, a, b);
}

double dbeta(double x, double a, double b) {
  return gsl_ran_beta_pdf(x, a, b);
}
