/* Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 James Theiler, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */
#include <math.h>
#include <rancommon.h>
__import rancommon;
__import specfun;

/* The bernoulli distribution has the form,

   prob(0) = 1-p, prob(1) = p

   */

unsigned int
gsl_ran_bernoulli (const gsl_rng * r, double p)
{
  double u = gsl_rng_uniform (r) ;

  if (u < p)
    {
      return 1 ;
    }
  else
    {
      return 0 ;
    }
}

double
gsl_ran_bernoulli_pdf (const unsigned int k, double p)
{
  if (k == 0)
    {
      return 1 - p ;
    }
  else if (k == 1)
    {
      return p ;
    }
  else
    {
      return 0 ;
    }
}

/* The beta distribution has the form

   p(x) dx = (Gamma(a + b)/(Gamma(a) Gamma(b))) x^(a-1) (1-x)^(b-1) dx

   The method used here is the one described in Knuth */

double
gsl_ran_beta (const gsl_rng * r, const double a, const double b)
{
  if ( (a <= 1.0) && (b <= 1.0) )
    {
      double U, V, X, Y;
      while (1)
        {
          U = gsl_rng_uniform_pos(r);
          V = gsl_rng_uniform_pos(r);
          X = pow(U, 1.0/a);
          Y = pow(V, 1.0/b);
          if ((X + Y ) <= 1.0)
            {
              if (X + Y > 0)
                {
                  return X/ (X + Y);
                }
              else
                {
                  double logX = log(U)/a;
                  double logY = log(V)/b;
                  double logM = logX > logY ? logX: logY;
                  logX -= logM;
                  logY -= logM;
                  return exp(logX - log(exp(logX) + exp(logY)));
                }
            }
        }
    }
  else
    {
      double x1 = gsl_ran_gamma (r, a, 1.0);
      double x2 = gsl_ran_gamma (r, b, 1.0);
      return x1 / (x1 + x2);
    }
}

double
gsl_ran_beta_pdf (const double x, const double a, const double b)
{
  if (x < 0 || x > 1)
    {
      return 0 ;
    }
  else 
    {
      double p;

      double gab = gsl_sf_lngamma (a + b);
      double ga = gsl_sf_lngamma (a);
      double gb = gsl_sf_lngamma (b);
      
      if (x == 0.0 || x == 1.0) 
        {
	  if (a > 1.0 && b > 1.0)
	    {
	      p = 0.0;
	    }
	  else
	    {
	      p = exp (gab - ga - gb) * pow (x, a - 1) * pow (1 - x, b - 1);
	    }
        }
      else
        {
          p = exp (gab - ga - gb + log(x) * (a - 1)  + log1p(-x) * (b - 1));
        }

      return p;
    }
}

/* The Bivariate Gaussian probability distribution is 

   p(x,y) dxdy = (1/(2 pi sigma_x sigma_y sqrt(c))) 
    exp(-((x/sigma_x)^2 + (y/sigma_y)^2 - 2 r (x/sigma_x)(y/sigma_y))/2c) dxdy 

   where c = 1-r^2
*/

void
gsl_ran_bivariate_gaussian (const gsl_rng * r, 
                            double sigma_x, double sigma_y, double rho,
                            double *x, double *y)
{
  double u, v, r2, scale;

  do
    {
      /* choose x,y in uniform square (-1,-1) to (+1,+1) */

      u = -1 + 2 * gsl_rng_uniform (r);
      v = -1 + 2 * gsl_rng_uniform (r);

      /* see if it is in the unit circle */
      r2 = u * u + v * v;
    }
  while (r2 > 1.0 || r2 == 0);

  scale = sqrt (-2.0 * log (r2) / r2);

  *x = sigma_x * u * scale;
  *y = sigma_y * (rho * u + sqrt(1 - rho*rho) * v) * scale;
}

double
gsl_ran_bivariate_gaussian_pdf (const double x, const double y, 
                                const double sigma_x, const double sigma_y,
                                const double rho)
{
  double u = x / sigma_x ;
  double v = y / sigma_y ;
  double c = 1 - rho*rho ;
  double p = (1 / (2 * M_PI * sigma_x * sigma_y * sqrt(c))) 
    * exp (-(u * u - 2 * rho * u * v + v * v) / (2 * c));
  return p;
}

/* The binomial distribution has the form,

   prob(k) =  n!/(k!(n-k)!) *  p^k (1-p)^(n-k) for k = 0, 1, ..., n

   This is the algorithm from Knuth */

/* Default binomial generator is now in binomial_tpe.c */

unsigned int
gsl_ran_binomial_knuth (const gsl_rng * r, double p, unsigned int n)
{
  unsigned int i, a, b, k = 0;

  while (n > 10)        /* This parameter is tunable */
    {
      double X;
      a = 1 + (n / 2);
      b = 1 + n - a;

      X = gsl_ran_beta (r, (double) a, (double) b);

      if (X >= p)
        {
          n = a - 1;
          p /= X;
        }
      else
        {
          k += a;
          n = b - 1;
          p = (p - X) / (1 - X);
        }
    }

  for (i = 0; i < n; i++)
    {
      double u = gsl_rng_uniform (r);
      if (u < p)
        k++;
    }

  return k;
}

double
gsl_ran_binomial_pdf (const unsigned int k, const double p,
                      const unsigned int n)
{
  if (k > n)
    {
      return 0;
    }
  else
    {
      double P;

      if (p == 0) 
        {
          P = (k == 0) ? 1 : 0;
        }
      else if (p == 1)
        {
          P = (k == n) ? 1 : 0;
        }
      else
        {
          double ln_Cnk = gsl_sf_lnchoose (n, k);
          P = ln_Cnk + k * log (p) + (n - k) * log1p (-p);
          P = exp (P);
        }

      return P;
    }
}

/* The Cauchy probability distribution is 

   p(x) dx = (1/(pi a)) (1 + (x/a)^2)^(-1) dx

   It is also known as the Lorentzian probability distribution */

double
gsl_ran_cauchy (const gsl_rng * r, const double a)
{
  double u;
  do
    {
      u = gsl_rng_uniform (r);
    }
  while (u == 0.5);

  return a * tan (M_PI * u);
}

double
gsl_ran_cauchy_pdf (const double x, const double a)
{
  double u = x / a;
  double p = (1 / (M_PI * a)) / (1 + u * u);
  return p;
}

/* The chisq distribution has the form

   p(x) dx = (1/(2*Gamma(nu/2))) (x/2)^(nu/2 - 1) exp(-x/2) dx

   for x = 0 ... +infty */

double
gsl_ran_chisq (const gsl_rng * r, const double nu)
{
  double chisq = 2 * gsl_ran_gamma (r, nu / 2, 1.0);
  return chisq;
}

double
gsl_ran_chisq_pdf (const double x, const double nu)
{
  if (x < 0)
    {
      return 0 ;
    }
  else
    {
      if(nu == 2.0)
        {
          return exp(-x/2.0) / 2.0;
        }
      else
        {
          double p;
          double lngamma = gsl_sf_lngamma (nu / 2);

          p = exp ((nu / 2 - 1) * log (x/2) - x/2 - lngamma) / 2;
          return p;
        }
    }
}

/* The Dirichlet probability distribution of order K-1 is 

     p(\theta_1,...,\theta_K) d\theta_1 ... d\theta_K = 
        (1/Z) \prod_i=1,K \theta_i^{alpha_i - 1} \delta(1 -\sum_i=1,K \theta_i)

   The normalization factor Z can be expressed in terms of gamma functions:

      Z = {\prod_i=1,K \Gamma(\alpha_i)} / {\Gamma( \sum_i=1,K \alpha_i)}  

   The K constants, \alpha_1,...,\alpha_K, must be positive. The K parameters, 
   \theta_1,...,\theta_K are nonnegative and sum to 1.

   The random variates are generated by sampling K values from gamma
   distributions with parameters a=\alpha_i, b=1, and renormalizing. 
   See A.M. Law, W.D. Kelton, Simulation Modeling and Analysis (1991).

   Gavin E. Crooks <gec@compbio.berkeley.edu> (2002)
*/

static void ran_dirichlet_small (const gsl_rng * r, const size_t K, const double alpha[], double theta[]);

void
gsl_ran_dirichlet (const gsl_rng * r, const size_t K,
                   const double alpha[], double theta[])
{
  size_t i;
  double norm = 0.0;

  for (i = 0; i < K; i++)
    {
      theta[i] = gsl_ran_gamma (r, alpha[i], 1.0);
    }
  
  for (i = 0; i < K; i++)
    {
      norm += theta[i];
    }

  if (norm < GSL_SQRT_DBL_MIN)  /* Handle underflow */
    {
      ran_dirichlet_small (r, K, alpha, theta);
      return;
    }

  for (i = 0; i < K; i++)
    {
      theta[i] /= norm;
    }
}


/* When the values of alpha[] are small, scale the variates to avoid
   underflow so that the result is not 0/0.  Note that the Dirichlet
   distribution is defined by a ratio of gamma functions so we can
   take out an arbitrary factor to keep the values in the range of
   double precision. */

static void 
ran_dirichlet_small (const gsl_rng * r, const size_t K,
                     const double alpha[], double theta[])
{
  size_t i;
  double norm = 0.0, umax = 0;

  for (i = 0; i < K; i++)
    {
      double u = log(gsl_rng_uniform_pos (r)) / alpha[i];
      
      theta[i] = u;

      if (u > umax || i == 0) {
        umax = u;
      }
    }
  
  for (i = 0; i < K; i++)
    {
      theta[i] = exp(theta[i] - umax);
    }
  
  for (i = 0; i < K; i++)
    {
      theta[i] = theta[i] * gsl_ran_gamma (r, alpha[i] + 1.0, 1.0);
    }

  for (i = 0; i < K; i++)
    {
      norm += theta[i];
    }

  for (i = 0; i < K; i++)
    {
      theta[i] /= norm;
    }
}

double
gsl_ran_dirichlet_pdf (const size_t K,
                       const double alpha[], const double theta[])
{
  return exp (gsl_ran_dirichlet_lnpdf (K, alpha, theta));
}

double
gsl_ran_dirichlet_lnpdf (const size_t K,
                         const double alpha[], const double theta[])
{
  /*We calculate the log of the pdf to minimize the possibility of overflow */
  size_t i;
  double log_p = 0.0;
  double sum_alpha = 0.0;

  for (i = 0; i < K; i++)
    {
      log_p += (alpha[i] - 1.0) * log (theta[i]);
    }

  for (i = 0; i < K; i++)
    {
      sum_alpha += alpha[i];
    }

  log_p += gsl_sf_lngamma (sum_alpha);

  for (i = 0; i < K; i++)
    {
      log_p -= gsl_sf_lngamma (alpha[i]);
    }

  return log_p;
}

/*
   Random Discrete Events
   
   Given K discrete events with different probabilities P[k]
   produce a value k consistent with its probability.

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License as
   published by the Free Software Foundation; either version 3 of the
   License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.  You should have received
   a copy of the GNU General Public License along with this library; if
   not, write to the Free Software Foundation, Inc., 51 Franklin Street,
   Fifth Floor, Boston, MA 02110-1301, USA.
*/

/*
 * Based on: Alastair J Walker, An efficient method for generating
 * discrete random variables with general distributions, ACM Trans
 * Math Soft 3, 253-256 (1977).  See also: D. E. Knuth, The Art of
 * Computer Programming, Volume 2 (Seminumerical algorithms), 3rd
 * edition, Addison-Wesley (1997), p120.

 * Walker's algorithm does some preprocessing, and provides two
 * arrays: floating point F[k] and integer A[k].  A value k is chosen
 * from 0..K-1 with equal likelihood, and then a uniform random number
 * u is compared to F[k].  If it is less than F[k], then k is
 * returned.  Otherwise, A[k] is returned.
   
 * Walker's original paper describes an O(K^2) algorithm for setting
 * up the F and A arrays.  I found this disturbing since I wanted to
 * use very large values of K.  I'm sure I'm not the first to realize
 * this, but in fact the preprocessing can be done in O(K) steps.

 * A figure of merit for the preprocessing is the average value for
 * the F[k]'s (that is, SUM_k F[k]/K); this corresponds to the
 * probability that k is returned, instead of A[k], thereby saving a
 * redirection.  Walker's O(K^2) preprocessing will generally improve
 * that figure of merit, compared to my cheaper O(K) method; from some
 * experiments with a perl script, I get values of around 0.6 for my
 * method and just under 0.75 for Walker's.  Knuth has pointed out
 * that finding _the_ optimum lookup tables, which maximize the
 * average F[k], is a combinatorially difficult problem.  But any
 * valid preprocessing will still provide O(1) time for the call to
 * gsl_ran_discrete().  I find that if I artificially set F[k]=1 --
 * ie, better than optimum! -- I get a speedup of maybe 20%, so that's
 * the maximum I could expect from the most expensive preprocessing.
 * Folding in the difference of 0.6 vs 0.75, I'd estimate that the
 * speedup would be less than 10%.

 * I've not implemented it here, but one compromise is to sort the
 * probabilities once, and then work from the two ends inward.  This
 * requires O(K log K), still lots cheaper than O(K^2), and from my
 * experiments with the perl script, the figure of merit is within
 * about 0.01 for K up to 1000, and no sign of diverging (in fact,
 * they seemed to be converging, but it's hard to say with just a
 * handful of runs).

 * The O(K) algorithm goes through all the p_k's and decides if they
 * are "smalls" or "bigs" according to whether they are less than or
 * greater than the mean value 1/K.  The indices to the smalls and the
 * bigs are put in separate stacks, and then we work through the
 * stacks together.  For each small, we pair it up with the next big
 * in the stack (Walker always wanted to pair up the smallest small
 * with the biggest big).  The small "borrows" from the big just
 * enough to bring the small up to mean.  This reduces the size of the
 * big, so the (smaller) big is compared again to the mean, and if it
 * is smaller, it gets "popped" from the big stack and "pushed" to the
 * small stack.  Otherwise, it stays put.  Since every time we pop a
 * small, we are able to deal with it right then and there, and we
 * never have to pop more than K smalls, then the algorithm is O(K).

 * This implementation sets up two separate stacks, and allocates K
 * elements between them.  Since neither stack ever grows, we do an
 * extra O(K) pass through the data to determine how many smalls and
 * bigs there are to begin with and allocate appropriately.  In all
 * there are 2*K*sizeof(double) transient bytes of memory that are
 * used than returned, and K*(sizeof(int)+sizeof(double)) bytes used
 * in the lookup table.
   
 * Walker spoke of using two random numbers (an integer 0..K-1, and a
 * floating point u in [0,1]), but Knuth points out that one can just
 * use the integer and fractional parts of K*u where u is in [0,1].
 * In fact, Knuth further notes that taking F'[k]=(k+F[k])/K, one can
 * directly compare u to F'[k] without having to explicitly set
 * u=K*u-int(K*u).

 * Usage:

 * Starting with an array of probabilities P, initialize and do
 * preprocessing with a call to:

 *    gsl_rng *r;
 *    gsl_ran_discrete_t *f;
 *    f = gsl_ran_discrete_preproc(K,P);
   
 * Then, whenever a random index 0..K-1 is desired, use

 *    k = gsl_ran_discrete(r,f);

 * Note that several different randevent struct's can be
 * simultaneously active.

 * Aside: A very clever alternative approach is described in
 * Abramowitz and Stegun, p 950, citing: Marsaglia, Random variables
 * and computers, Proc Third Prague Conference in Probability Theory,
 * 1962.  A more accesible reference is: G. Marsaglia, Generating
 * discrete random numbers in a computer, Comm ACM 6, 37-38 (1963).
 * If anybody is interested, I (jt) have also coded up this version as
 * part of another software package.  However, I've done some
 * comparisons, and the Walker method is both faster and more stingy
 * with memory.  So, in the end I decided not to include it with the
 * GSL package.
   
 * Written 26 Jan 1999, James Theiler, jt@lanl.gov
 * Adapted to GSL, 30 Jan 1999, jt

 */

#define DEBUG 0
#define KNUTH_CONVENTION 1      /* Saves a few steps of arithmetic
                                 * in the call to gsl_ran_discrete()
                                 */

typedef struct {                /* struct for Walker algorithm */
    size_t K;
    size_t *A;
    double *F;
} gsl_ran_discrete_t;

/*** Begin Stack (this code is used just in this file) ***/

/* Stack code converted to use unsigned indices (i.e. s->i == 0 now
   means an empty stack, instead of -1), for consistency and to give a
   bigger allowable range. BJG */

typedef struct {
    size_t N;                      /* max number of elts on stack */
    size_t *v;                     /* array of values on the stack */
    size_t i;                      /* index of top of stack */
} gsl_stack_t;

static gsl_stack_t *
new_stack(size_t N) {
    gsl_stack_t *s;
    s = (gsl_stack_t *)malloc(sizeof(gsl_stack_t));
    s->N = N;
    s->i = 0;                  /* indicates stack is empty */
    s->v = (size_t *)malloc(sizeof(size_t)*N);
    return s;
}

static int
push_stack(gsl_stack_t *s, size_t v)
{
    if ((s->i) >= (s->N)) {
      return -1; /* stack overflow (shouldn't happen) */
    }
    (s->v)[s->i] = v;
    s->i += 1;
    return 0;
}

static size_t pop_stack(gsl_stack_t *s)
{
    if ((s->i) == 0) {
      GSL_ERROR ("internal error - stack exhausted", GSL_ESANITY);
    }
    s->i -= 1;
    return ((s->v)[s->i]);
}

static inline size_t size_stack(const gsl_stack_t *s)
{
    return s->i;
}

static void free_stack(gsl_stack_t *s)
{
    free((char *)(s->v));
    free((char *)s);
}

/*** End Stack ***/


/*** Begin Walker's Algorithm ***/

gsl_ran_discrete_t *
gsl_ran_discrete_preproc(size_t Kevents, const double *ProbArray)
{
    size_t k,b,s;
    gsl_ran_discrete_t *g;
    size_t nBigs, nSmalls;
    gsl_stack_t *Bigs;
    gsl_stack_t *Smalls;
    double *E;
    double pTotal = 0.0, mean, d;
    
    if (Kevents < 1) {
      /* Could probably treat Kevents=1 as a special case */

      GSL_ERROR_VAL ("number of events must be a positive integer", 
                        GSL_EINVAL, 0);
    }

    /* Make sure elements of ProbArray[] are positive.
     * Won't enforce that sum is unity; instead will just normalize
     */

    for (k=0; k<Kevents; ++k) {
        if (ProbArray[k] < 0) {
          GSL_ERROR_VAL ("probabilities must be non-negative",
                            GSL_EINVAL, 0) ;
        }
        pTotal += ProbArray[k];
    }

    /* Begin setting up the main "object" (just a struct, no steroids) */
    g = (gsl_ran_discrete_t *)malloc(sizeof(gsl_ran_discrete_t));
    g->K = Kevents;
    g->F = (double *)malloc(sizeof(double)*Kevents);
    g->A = (size_t *)malloc(sizeof(size_t)*Kevents);

    E = (double *)malloc(sizeof(double)*Kevents);

    if (E==NULL) {
      GSL_ERROR_VAL ("Cannot allocate memory for randevent", GSL_ENOMEM, 0);
    }

    for (k=0; k<Kevents; ++k) {
        E[k] = ProbArray[k]/pTotal;
    }

    /* Now create the Bigs and the Smalls */
    mean = 1.0/Kevents;
    nSmalls=nBigs=0;
    {
      /* Temporarily use which[k] = g->A[k] to indicate small or large */
      size_t * const which = g->A;

      for (k=0; k<Kevents; ++k) {
        if (E[k] < mean) { 
          ++nSmalls; which[k] = 0;
        } else { 
          ++nBigs; which[k] = 1; 
        }
      }

      Bigs   = new_stack(nBigs);
      Smalls = new_stack(nSmalls);
      for (k=0; k<Kevents; ++k) {
        gsl_stack_t * Dest = which[k] ? Bigs : Smalls;
        int status = push_stack(Dest,k);
        if (status)
          GSL_ERROR_VAL ("failed to build stacks", GSL_EFAILED, 0);
      }
    }

    /* Now work through the smalls */
    while (size_stack(Smalls) > 0) {
        s = pop_stack(Smalls);
        if (size_stack(Bigs) == 0) {
            (g->A)[s]=s;
            (g->F)[s]=1.0;
            continue;
        }
        b = pop_stack(Bigs);
        (g->A)[s]=b;
        (g->F)[s]=Kevents*E[s];
#if DEBUG
        fprintf(stderr,"s=%2d, A=%2d, F=%.4f\n",s,(g->A)[s],(g->F)[s]);
#endif        
        d = mean - E[s];
        E[s] += d;              /* now E[s] == mean */
        E[b] -= d;
        if (E[b] < mean) {
            push_stack(Smalls,b); /* no longer big, join ranks of the small */
        }
        else if (E[b] > mean) {
            push_stack(Bigs,b); /* still big, put it back where you found it */
        }
        else {
            /* E[b]==mean implies it is finished too */
            (g->A)[b]=b;
            (g->F)[b]=1.0;
        }
    }
    while (size_stack(Bigs) > 0) {
        b = pop_stack(Bigs);
        (g->A)[b]=b;
        (g->F)[b]=1.0;
    }
    /* Stacks have been emptied, and A and F have been filled */

    if ( size_stack(Smalls) != 0) {
      GSL_ERROR_VAL ("Smalls stack has not been emptied",
                     GSL_ESANITY, 0 );
    }
    
#if 0
    /* if 1, then artificially set all F[k]'s to unity.  This will
     * give wrong answers, but you'll get them faster.  But, not
     * that much faster (I get maybe 20%); that's an upper bound
     * on what the optimal preprocessing would give.
     */
    for (k=0; k<Kevents; ++k) {
        (g->F)[k] = 1.0;
    }
#endif

#if KNUTH_CONVENTION
    /* For convenience, set F'[k]=(k+F[k])/K */
    /* This saves some arithmetic in gsl_ran_discrete(); I find that
     * it doesn't actually make much difference.
     */
    for (k=0; k<Kevents; ++k) {
        (g->F)[k] += k;
        (g->F)[k] /= Kevents;
    }
#endif    

    free_stack(Bigs);
    free_stack(Smalls);
    free((char *)E);

    return g;
}

size_t
gsl_ran_discrete(const gsl_rng *r, const gsl_ran_discrete_t *g)
{
    size_t c=0;
    double u,f;
    u = gsl_rng_uniform(r);
#if KNUTH_CONVENTION
    c = (u*(g->K));
#else
    u *= g->K;
    c = u;
    u -= c;
#endif
    f = (g->F)[c];
    /* fprintf(stderr,"c,f,u: %d %.4f %f\n",c,f,u); */
    if (f == 1.0) return c;

    if (u < f) {
        return c;
    }
    else {
        return (g->A)[c];
    }
}

void gsl_ran_discrete_free(gsl_ran_discrete_t *g)
{
    RETURN_IF_NULL (g);
    free((char *)(g->A));
    free((char *)(g->F));
    free((char *)g);
}

double
gsl_ran_discrete_pdf(size_t k, const gsl_ran_discrete_t *g)
{
    size_t i,K;
    double f,p=0;
    K= g->K;
    if (k>K) return 0;
    for (i=0; i<K; ++i) {
        f = (g->F)[i];
#if KNUTH_CONVENTION
        f = K*f-i;
#endif        
        if (i==k) {
            p += f;
        } else if (k == (g->A)[i]) {
            p += 1.0 - f;
        }
    }
    return p/K;
}

/* This is the uniform distribution in the range [a, b)

   p(x) dx = 1/(b-a) dx   if  a <= x < b
   .....   = 0            otherwise 

 */

double
gsl_ran_flat (const gsl_rng * r, const double a, const double b)
{
  double u = gsl_rng_uniform (r);

  /* A uniform distribution over [a,b) */

  return a * (1 - u) + b * u;
}

double
gsl_ran_flat_pdf (double x, const double a, const double b)
{
  if (x < b && x >= a)
    {
      return 1 / (b - a);
    }
  else
    {
      return 0;
    }
}

static double gamma_large (const gsl_rng * r, const double a);
static double gamma_frac (const gsl_rng * r, const double a);

/* The Gamma distribution of order a>0 is defined by:

   p(x) dx = {1 / \Gamma(a) b^a } x^{a-1} e^{-x/b} dx

   for x>0.  If X and Y are independent gamma-distributed random
   variables of order a1 and a2 with the same scale parameter b, then
   X+Y has gamma distribution of order a1+a2.

   The algorithms below are from Knuth, vol 2, 2nd ed, p. 129. */

double
gsl_ran_gamma_knuth (const gsl_rng * r, const double a, const double b)
{
  /* assume a > 0 */
  unsigned int na = floor (a);

  if(a >= UINT_MAX) 
    {
      return b * (gamma_large (r, floor (a)) + gamma_frac (r, a - floor (a)));
    }
  else if (a == na)
    {
      return b * gsl_ran_gamma_int (r, na);
    }
  else if (na == 0)
    {
      return b * gamma_frac (r, a);
    }
  else
    {
      return b * (gsl_ran_gamma_int (r, na) + gamma_frac (r, a - na)) ;
    }
}

double
gsl_ran_gamma_int (const gsl_rng * r, const unsigned int a)
{
  if (a < 12)
    {
      unsigned int i;
      double prod = 1;

      for (i = 0; i < a; i++)
        {
          prod *= gsl_rng_uniform_pos (r);
        }

      /* Note: for 12 iterations we are safe against underflow, since
         the smallest positive random number is O(2^-32). This means
         the smallest possible product is 2^(-12*32) = 10^-116 which
         is within the range of double precision. */

      return -log (prod);
    }
  else
    {
      return gamma_large (r, (double) a);
    }
}

static double
gamma_large (const gsl_rng * r, const double a)
{
  /* Works only if a > 1, and is most efficient if a is large

     This algorithm, reported in Knuth, is attributed to Ahrens.  A
     faster one, we are told, can be found in: J. H. Ahrens and
     U. Dieter, Computing 12 (1974) 223-246.  */

  double sqa, x, y, v;
  sqa = sqrt (2 * a - 1);
  do
    {
      do
        {
          y = tan (M_PI * gsl_rng_uniform (r));
          x = sqa * y + a - 1;
        }
      while (x <= 0);
      v = gsl_rng_uniform (r);
    }
  while (v > (1 + y * y) * exp ((a - 1) * log (x / (a - 1)) - sqa * y));

  return x;
}

static double
gamma_frac (const gsl_rng * r, const double a)
{
  /* This is exercise 16 from Knuth; see page 135, and the solution is
     on page 551.  */

  double p, q, x, u, v;

  if (a == 0) {
    return 0;
  }

  p = M_E / (a + M_E);
  do
    {
      u = gsl_rng_uniform (r);
      v = gsl_rng_uniform_pos (r);

      if (u < p)
        {
          x = exp ((1 / a) * log (v));
          q = exp (-x);
        }
      else
        {
          x = 1 - log (v);
          q = exp ((a - 1) * log (x));
        }
    }
  while (gsl_rng_uniform (r) >= q);

  return x;
}

double
gsl_ran_gamma_pdf (const double x, const double a, const double b)
{
  if (x < 0)
    {
      return 0 ;
    }
  else if (x == 0)
    {
      if (a == 1)
        return 1/b ;
      else
        return 0 ;
    }
  else if (a == 1)
    {
      return exp(-x/b)/b ;
    }
  else 
    {
      double p;
      double lngamma = gsl_sf_lngamma (a);
      p = exp ((a - 1) * log (x/b) - x/b - lngamma)/b;
      return p;
    }
}


/* New version based on Marsaglia and Tsang, "A Simple Method for
 * generating gamma variables", ACM Transactions on Mathematical
 * Software, Vol 26, No 3 (2000), p363-372.
 *
 * Implemented by J.D.Lamb@btinternet.com, minor modifications for GSL
 * by Brian Gough
 */

double
gsl_ran_gamma_mt (const gsl_rng * r, const double a, const double b)
{
  return gsl_ran_gamma (r, a, b);
}

double
gsl_ran_gamma (const gsl_rng * r, const double a, const double b)
{
  /* assume a > 0 */

  if (a < 1)
    {
      double u = gsl_rng_uniform_pos (r);
      return gsl_ran_gamma (r, 1.0 + a, b) * pow (u, 1.0 / a);
    }

  {
    double x, v, u;
    double d = a - 1.0 / 3.0;
    double c = (1.0 / 3.0) / sqrt (d);

    while (1)
      {
        do
          {
            x = gsl_ran_gaussian_ziggurat (r, 1.0);
            v = 1.0 + c * x;
          }
        while (v <= 0);

        v = v * v * v;
        u = gsl_rng_uniform_pos (r);

        if (u < 1 - 0.0331 * x * x * x * x) 
          break;

        if (log (u) < 0.5 * x * x + d * (1 - v + log (v)))
          break;
      }
    
    return b * d * v;
  }
}

/* position of right-most step */
#define PARAM_R 3.44428647676

/* tabulated values for the heigt of the Ziggurat levels */
static const double ytab[128] = {
  1, 0.963598623011, 0.936280813353, 0.913041104253,
  0.892278506696, 0.873239356919, 0.855496407634, 0.838778928349,
  0.822902083699, 0.807732738234, 0.793171045519, 0.779139726505,
  0.765577436082, 0.752434456248, 0.739669787677, 0.727249120285,
  0.715143377413, 0.703327646455, 0.691780377035, 0.68048276891,
  0.669418297233, 0.65857233912, 0.647931876189, 0.637485254896,
  0.62722199145, 0.617132611532, 0.607208517467, 0.597441877296,
  0.587825531465, 0.578352913803, 0.569017984198, 0.559815170911,
  0.550739320877, 0.541785656682, 0.532949739145, 0.524227434628,
  0.515614886373, 0.507108489253, 0.498704867478, 0.490400854812,
  0.482193476986, 0.47407993601, 0.466057596125, 0.458123971214,
  0.450276713467, 0.442513603171, 0.434832539473, 0.427231532022,
  0.419708693379, 0.41226223212, 0.404890446548, 0.397591718955,
  0.390364510382, 0.383207355816, 0.376118859788, 0.369097692334,
  0.362142585282, 0.355252328834, 0.348425768415, 0.341661801776,
  0.334959376311, 0.328317486588, 0.321735172063, 0.31521151497,
  0.308745638367, 0.302336704338, 0.29598391232, 0.289686497571,
  0.283443729739, 0.27725491156, 0.271119377649, 0.265036493387,
  0.259005653912, 0.253026283183, 0.247097833139, 0.241219782932,
  0.235391638239, 0.229612930649, 0.223883217122, 0.218202079518,
  0.212569124201, 0.206983981709, 0.201446306496, 0.195955776745,
  0.190512094256, 0.185114984406, 0.179764196185, 0.174459502324,
  0.169200699492, 0.1639876086, 0.158820075195, 0.153697969964,
  0.148621189348, 0.143589656295, 0.138603321143, 0.133662162669,
  0.128766189309, 0.123915440582, 0.119109988745, 0.114349940703,
  0.10963544023, 0.104966670533, 0.100343857232, 0.0957672718266,
  0.0912372357329, 0.0867541250127, 0.082318375932, 0.0779304915295,
  0.0735910494266, 0.0693007111742, 0.065060233529, 0.0608704821745,
  0.056732448584, 0.05264727098, 0.0486162607163, 0.0446409359769,
  0.0407230655415, 0.0368647267386, 0.0330683839378, 0.0293369977411,
  0.0256741818288, 0.0220844372634, 0.0185735200577, 0.0151490552854,
  0.0118216532614, 0.00860719483079, 0.00553245272614, 0.00265435214565
};

/* tabulated values for 2^24 times x[i]/x[i+1],
 * used to accept for U*x[i+1]<=x[i] without any floating point operations */
static const unsigned long ktab[128] = {
  0, 12590644, 14272653, 14988939,
  15384584, 15635009, 15807561, 15933577,
  16029594, 16105155, 16166147, 16216399,
  16258508, 16294295, 16325078, 16351831,
  16375291, 16396026, 16414479, 16431002,
  16445880, 16459343, 16471578, 16482744,
  16492970, 16502368, 16511031, 16519039,
  16526459, 16533352, 16539769, 16545755,
  16551348, 16556584, 16561493, 16566101,
  16570433, 16574511, 16578353, 16581977,
  16585398, 16588629, 16591685, 16594575,
  16597311, 16599901, 16602354, 16604679,
  16606881, 16608968, 16610945, 16612818,
  16614592, 16616272, 16617861, 16619363,
  16620782, 16622121, 16623383, 16624570,
  16625685, 16626730, 16627708, 16628619,
  16629465, 16630248, 16630969, 16631628,
  16632228, 16632768, 16633248, 16633671,
  16634034, 16634340, 16634586, 16634774,
  16634903, 16634972, 16634980, 16634926,
  16634810, 16634628, 16634381, 16634066,
  16633680, 16633222, 16632688, 16632075,
  16631380, 16630598, 16629726, 16628757,
  16627686, 16626507, 16625212, 16623794,
  16622243, 16620548, 16618698, 16616679,
  16614476, 16612071, 16609444, 16606571,
  16603425, 16599973, 16596178, 16591995,
  16587369, 16582237, 16576520, 16570120,
  16562917, 16554758, 16545450, 16534739,
  16522287, 16507638, 16490152, 16468907,
  16442518, 16408804, 16364095, 16301683,
  16207738, 16047994, 15704248, 15472926
};

/* tabulated values of 2^{-24}*x[i] */
static const double wtab[128] = {
  1.62318314817e-08, 2.16291505214e-08, 2.54246305087e-08, 2.84579525938e-08,
  3.10340022482e-08, 3.33011726243e-08, 3.53439060345e-08, 3.72152672658e-08,
  3.8950989572e-08, 4.05763964764e-08, 4.21101548915e-08, 4.35664624904e-08,
  4.49563968336e-08, 4.62887864029e-08, 4.75707945735e-08, 4.88083237257e-08,
  5.00063025384e-08, 5.11688950428e-08, 5.22996558616e-08, 5.34016475624e-08,
  5.44775307871e-08, 5.55296344581e-08, 5.65600111659e-08, 5.75704813695e-08,
  5.85626690412e-08, 5.95380306862e-08, 6.04978791776e-08, 6.14434034901e-08,
  6.23756851626e-08, 6.32957121259e-08, 6.42043903937e-08, 6.51025540077e-08,
  6.59909735447e-08, 6.68703634341e-08, 6.77413882848e-08, 6.8604668381e-08,
  6.94607844804e-08, 7.03102820203e-08, 7.11536748229e-08, 7.1991448372e-08,
  7.2824062723e-08, 7.36519550992e-08, 7.44755422158e-08, 7.52952223703e-08,
  7.61113773308e-08, 7.69243740467e-08, 7.77345662086e-08, 7.85422956743e-08,
  7.93478937793e-08, 8.01516825471e-08, 8.09539758128e-08, 8.17550802699e-08,
  8.25552964535e-08, 8.33549196661e-08, 8.41542408569e-08, 8.49535474601e-08,
  8.57531242006e-08, 8.65532538723e-08, 8.73542180955e-08, 8.8156298059e-08,
  8.89597752521e-08, 8.97649321908e-08, 9.05720531451e-08, 9.138142487e-08,
  9.21933373471e-08, 9.30080845407e-08, 9.38259651738e-08, 9.46472835298e-08,
  9.54723502847e-08, 9.63014833769e-08, 9.71350089201e-08, 9.79732621669e-08,
  9.88165885297e-08, 9.96653446693e-08, 1.00519899658e-07, 1.0138063623e-07,
  1.02247952126e-07, 1.03122261554e-07, 1.04003996769e-07, 1.04893609795e-07,
  1.05791574313e-07, 1.06698387725e-07, 1.07614573423e-07, 1.08540683296e-07,
  1.09477300508e-07, 1.1042504257e-07, 1.11384564771e-07, 1.12356564007e-07,
  1.13341783071e-07, 1.14341015475e-07, 1.15355110887e-07, 1.16384981291e-07,
  1.17431607977e-07, 1.18496049514e-07, 1.19579450872e-07, 1.20683053909e-07,
  1.21808209468e-07, 1.2295639141e-07, 1.24129212952e-07, 1.25328445797e-07,
  1.26556042658e-07, 1.27814163916e-07, 1.29105209375e-07, 1.30431856341e-07,
  1.31797105598e-07, 1.3320433736e-07, 1.34657379914e-07, 1.36160594606e-07,
  1.37718982103e-07, 1.39338316679e-07, 1.41025317971e-07, 1.42787873535e-07,
  1.44635331499e-07, 1.4657889173e-07, 1.48632138436e-07, 1.50811780719e-07,
  1.53138707402e-07, 1.55639532047e-07, 1.58348931426e-07, 1.61313325908e-07,
  1.64596952856e-07, 1.68292495203e-07, 1.72541128694e-07, 1.77574279496e-07,
  1.83813550477e-07, 1.92166040885e-07, 2.05295471952e-07, 2.22600839893e-07
};


double
gsl_ran_gaussian_ziggurat (const gsl_rng * r, const double sigma)
{
  unsigned long int i, j;
  int sign;
  double x, y;

  const unsigned long int range = r->type->max - r->type->min;
  const unsigned long int offset = r->type->min;

  while (1)
    {
      if (range >= 0xFFFFFFFF)
        {
          unsigned long int k = gsl_rng_get(r) - offset;
          i = (k & 0xFF);
          j = (k >> 8) & 0xFFFFFF;
        }
      else if (range >= 0x00FFFFFF)
        {
          unsigned long int k1 = gsl_rng_get(r) - offset;
          unsigned long int k2 = gsl_rng_get(r) - offset;
          i = (k1 & 0xFF);
          j = (k2 & 0x00FFFFFF);
        }
      else
        {
          i = gsl_rng_uniform_int (r, 256); /*  choose the step */
          j = gsl_rng_uniform_int (r, 16777216);  /* sample from 2^24 */
        }

      sign = (i & 0x80) ? +1 : -1;
      i &= 0x7f;

      x = j * wtab[i];

      if (j < ktab[i])
        break;

      if (i < 127)
        {
          double y0, y1, U1;
          y0 = ytab[i];
          y1 = ytab[i + 1];
          U1 = gsl_rng_uniform (r);
          y = y1 + (y0 - y1) * U1;
        }
      else
        {
          double U1, U2;
          U1 = 1.0 - gsl_rng_uniform (r);
          U2 = gsl_rng_uniform (r);
          x = PARAM_R - log (U1) / PARAM_R;
          y = exp (-PARAM_R * (x - 0.5 * PARAM_R)) * U2;
        }

      if (y < exp (-0.5 * x * x))
        break;
    }

  return sign * sigma * x;
}

/* Of the two methods provided below, I think the Polar method is more
 * efficient, but only when you are actually producing two random
 * deviates.  We don't produce two, because then we'd have to save one
 * in a static variable for the next call, and that would screws up
 * re-entrant or threaded code, so we only produce one.  This makes
 * the Ratio method suddenly more appealing.
 *
 * [Added by Charles Karney] We use Leva's implementation of the Ratio
 * method which avoids calling log() nearly all the time and makes the
 * Ratio method faster than the Polar method (when it produces just one
 * result per call).  Timing per call (gcc -O2 on 866MHz Pentium,
 * average over 10^8 calls)
 *
 *   Polar method: 660 ns
 *   Ratio method: 368 ns
 *
 */

/* Polar (Box-Mueller) method; See Knuth v2, 3rd ed, p122 */

double
gsl_ran_gaussian (const gsl_rng * r, const double sigma)
{
  double x, y, r2;

  do
    {
      /* choose x,y in uniform square (-1,-1) to (+1,+1) */
      x = -1 + 2 * gsl_rng_uniform_pos (r);
      y = -1 + 2 * gsl_rng_uniform_pos (r);

      /* see if it is in the unit circle */
      r2 = x * x + y * y;
    }
  while (r2 > 1.0 || r2 == 0);

  /* Box-Muller transform */
  return sigma * y * sqrt (-2.0 * log (r2) / r2);
}

/* Ratio method (Kinderman-Monahan); see Knuth v2, 3rd ed, p130.
 * K+M, ACM Trans Math Software 3 (1977) 257-260.
 *
 * [Added by Charles Karney] This is an implementation of Leva's
 * modifications to the original K+M method; see:
 * J. L. Leva, ACM Trans Math Software 18 (1992) 449-453 and 454-455. */

double
gsl_ran_gaussian_ratio_method (const gsl_rng * r, const double sigma)
{
  double u, v, x, y, Q;
  const double s = 0.449871;    /* Constants from Leva */
  const double t = -0.386595;
  const double a = 0.19600;
  const double b = 0.25472;
  const double r1 = 0.27597;
  const double r2 = 0.27846;

  do                            /* This loop is executed 1.369 times on average  */
    {
      /* Generate a point P = (u, v) uniform in a rectangle enclosing
         the K+M region v^2 <= - 4 u^2 log(u). */

      /* u in (0, 1] to avoid singularity at u = 0 */
      u = 1 - gsl_rng_uniform (r);

      /* v is in the asymmetric interval [-0.5, 0.5).  However v = -0.5
         is rejected in the last part of the while clause.  The
         resulting normal deviate is strictly symmetric about 0
         (provided that v is symmetric once v = -0.5 is excluded). */
      v = gsl_rng_uniform (r) - 0.5;

      /* Constant 1.7156 > sqrt(8/e) (for accuracy); but not by too
         much (for efficiency). */
      v *= 1.7156;

      /* Compute Leva's quadratic form Q */
      x = u - s;
      y = fabs (v) - t;
      Q = x * x + y * (a * y - b * x);

      /* Accept P if Q < r1 (Leva) */
      /* Reject P if Q > r2 (Leva) */
      /* Accept if v^2 <= -4 u^2 log(u) (K+M) */
      /* This final test is executed 0.012 times on average. */
    }
  while (Q >= r1 && (Q > r2 || v * v > -4 * u * u * log (u)));

  return sigma * (v / u);       /* Return slope */
}

double
gsl_ran_gaussian_pdf (const double x, const double sigma)
{
  double u = x / fabs (sigma);
  double p = (1 / (sqrt (2 * M_PI) * fabs (sigma))) * exp (-u * u / 2);
  return p;
}

double
gsl_ran_ugaussian (const gsl_rng * r)
{
  return gsl_ran_gaussian (r, 1.0);
}

double
gsl_ran_ugaussian_ratio_method (const gsl_rng * r)
{
  return gsl_ran_gaussian_ratio_method (r, 1.0);
}

double
gsl_ran_ugaussian_pdf (const double x)
{
  return gsl_ran_gaussian_pdf (x, 1.0);
}

/* The sum of N samples from an exponential distribution gives an
   Erlang distribution

   p(x) dx = x^(n-1) exp (-x/a) / ((n-1)!a^n) dx

   for x = 0 ... +infty */

double
gsl_ran_erlang (const gsl_rng * r, const double a, const double n)
{
  return gsl_ran_gamma (r, n, a);
}

double
gsl_ran_erlang_pdf (const double x, const double a, const double n)
{
  if (x <= 0) 
    {
      return 0 ;
    }
  else
    {
      double p;
      double lngamma = gsl_sf_lngamma (n);

      p = exp ((n - 1) * log (x/a) - x/a - lngamma) / a;
      return p;
    }
}

/* The exponential distribution has the form

   p(x) dx = exp(-x/mu) dx/mu

   for x = 0 ... +infty */

double
gsl_ran_exponential (const gsl_rng * r, const double mu)
{
  double u = gsl_rng_uniform (r);

  return -mu * log1p (-u);
}

double
gsl_ran_exponential_pdf (const double x, const double mu)
{
  if (x < 0)
    {
      return 0 ;
    }
  else
    {
      double p = exp (-x/mu)/mu;
      
      return p;
    }
}

/* The exponential power probability distribution is  

   p(x) dx = (1/(2 a Gamma(1+1/b))) * exp(-|x/a|^b) dx

   for -infty < x < infty. For b = 1 it reduces to the Laplace
   distribution. 

   The exponential power distribution is related to the gamma
   distribution by E = a * pow(G(1/b),1/b), where E is an exponential
   power variate and G is a gamma variate.

   We use this relation for b < 1. For b >=1 we use rejection methods
   based on the laplace and gaussian distributions which should be
   faster.  For b>4 we revert to the gamma method.

   See P. R. Tadikamalla, "Random Sampling from the Exponential Power
   Distribution", Journal of the American Statistical Association,
   September 1980, Volume 75, Number 371, pages 683-686.
   
*/

double
gsl_ran_exppow (const gsl_rng * r, const double a, const double b)
{
  if (b < 1 || b > 4)
    {
      double u = gsl_rng_uniform (r);
      double v = gsl_ran_gamma (r, 1 / b, 1.0);
      double z = a * pow (v, 1 / b);

      if (u > 0.5)
        {
          return z;
        }
      else
        {
          return -z;
        }
    }
  else if (b == 1)
    {
      /* Laplace distribution */
      return gsl_ran_laplace (r, a);
    }
  else if (b < 2)
    {
      /* Use laplace distribution for rejection method, from Tadikamalla */

      double x, h, u;

      double B = pow (1 / b, 1 / b);

      do
        {
          x = gsl_ran_laplace (r, B);
          u = gsl_rng_uniform_pos (r);
          h = -pow (fabs (x), b) + fabs (x) / B - 1 + (1 / b);
        }
      while (log (u) > h);

      return a * x;
    }
  else if (b == 2)
    {
      /* Gaussian distribution */
      return gsl_ran_gaussian (r, a / sqrt (2.0));
    }
  else
    {
      /* Use gaussian for rejection method, from Tadikamalla */

      double x, h, u;

      double B = pow (1 / b, 1 / b);

      do
        {
          x = gsl_ran_gaussian (r, B);
          u = gsl_rng_uniform_pos (r);
          h = -pow (fabs (x), b) + (x * x) / (2 * B * B) + (1 / b) - 0.5;
        }
      while (log (u) > h);

      return a * x;
    }
}

double
gsl_ran_exppow_pdf (const double x, const double a, const double b)
{
  double p;
  double lngamma = gsl_sf_lngamma (1 + 1 / b);
  p = (1 / (2 * a)) * exp (-pow (fabs (x / a), b) - lngamma);
  return p;
}

/* The two-sided exponential probability distribution is  

   p(x) dx = (1/(2 a)) * exp(-|x/a|) dx

   for -infty < x < infty. It is also known as the Laplace distribution.  */

double
gsl_ran_laplace (const gsl_rng * r, const double a)
{
  double u;
  do
    {
      u = 2 * gsl_rng_uniform (r) - 1.0;
    }
  while (u == 0.0);

  if (u < 0)
    {
      return a * log (-u);
    }
  else
    {
      return -a * log (u);
    }
}

double
gsl_ran_laplace_pdf (const double x, const double a)
{
  double p = (1/(2*a)) * exp (-fabs (x)/a);
  return p;
}

/* The F distribution has the form

   p(x) dx = (nu1^(nu1/2) nu2^(nu2/2) Gamma((nu1 + nu2)/2) /
   Gamma(nu1/2) Gamma(nu2/2)) *
   x^(nu1/2 - 1) (nu2 + nu1 * x)^(-nu1/2 -nu2/2) dx

   The method used here is the one described in Knuth */

double
gsl_ran_fdist (const gsl_rng * r, const double nu1, const double nu2)
{

  double Y1 =  gsl_ran_gamma (r, nu1 / 2, 2.0);
  double Y2 =  gsl_ran_gamma (r, nu2 / 2, 2.0);

  double f = (Y1 * nu2) / (Y2 * nu1);

  return f;
}

double
gsl_ran_fdist_pdf (const double x, const double nu1, const double nu2)
{
  if (x < 0)
    {
      return 0 ;
    }
  else
    {
      double p;
      double lglg = (nu1 / 2) * log (nu1) + (nu2 / 2) * log (nu2) ;

      double lg12 = gsl_sf_lngamma ((nu1 + nu2) / 2);
      double lg1 = gsl_sf_lngamma (nu1 / 2);
      double lg2 = gsl_sf_lngamma (nu2 / 2);

      p =
	exp (lglg + lg12 - lg1 - lg2 + (nu1 / 2 - 1) * log (x) -
	     ((nu1 + nu2) / 2) * log (nu2 + nu1 * x));

      return p;
    }
}

double
gsl_ran_gaussian_tail (const gsl_rng * r, const double a, const double sigma)
{
  /* Returns a gaussian random variable larger than a
   * This implementation does one-sided upper-tailed deviates.
   */

  double s = a / sigma;

  if (s < 1)
    {
      /* For small s, use a direct rejection method. The limit s < 1
         can be adjusted to optimise the overall efficiency */

      double x;

      do
        {
          x = gsl_ran_gaussian (r, 1.0);
        }
      while (x < s);
      return x * sigma;
    }
  else
    {
      /* Use the "supertail" deviates from the last two steps
       * of Marsaglia's rectangle-wedge-tail method, as described
       * in Knuth, v2, 3rd ed, pp 123-128.  (See also exercise 11, p139,
       * and the solution, p586.)
       */

      double u, v, x;

      do
        {
          u = gsl_rng_uniform (r);
          do
            {
              v = gsl_rng_uniform (r);
            }
          while (v == 0.0);
          x = sqrt (s * s - 2 * log (v));
        }
      while (x * u > s);
      return x * sigma;
    }
}

double
gsl_ran_gaussian_tail_pdf (const double x, const double a, const double sigma)
{
  if (x < a)
    {
      return 0;
    }
  else
    {
      double N, p;
      double u = x / sigma ;

      double f = gsl_sf_erfc (a / (sqrt (2.0) * sigma));

      N = 0.5 * f;

      p = (1 / (N * sqrt (2 * M_PI) * sigma)) * exp (-u * u / 2);

      return p;
    }
}

double
gsl_ran_ugaussian_tail (const gsl_rng * r, const double a)
{
  return gsl_ran_gaussian_tail (r, a, 1.0) ;
}

double
gsl_ran_ugaussian_tail_pdf (const double x, const double a)
{
  return gsl_ran_gaussian_tail_pdf (x, a, 1.0) ;
}

/* Geometric distribution (bernoulli trial with probability p) 

   prob(k) =  p (1 - p)^(k-1) for n = 1, 2, 3, ...

   It gives the distribution of "waiting times" for an event that
   occurs with probability p. */

unsigned int
gsl_ran_geometric (const gsl_rng * r, const double p)
{
  double u = gsl_rng_uniform_pos (r);

  unsigned int k;

  if (p == 1)
    {
      k = 1;
    }
  else
    {
      k = log (u) / log (1 - p) + 1;
    }

  return k;
}

double
gsl_ran_geometric_pdf (const unsigned int k, const double p)
{
  if (k == 0)
    {
      return 0 ;
    }
  else if (k == 1)
    {
      return p ;
    }
  else
    {
      double P = p * pow (1 - p, k - 1.0);
      return P;
    }
}

/* The Type I Gumbel distribution has the form,

   p(x) dx = a b exp(-(b exp(-ax) + ax)) dx

   and the Type II Gumbel distribution has the form,

   p(x) dx = b a x^-(a+1) exp(-b x^-a)) dx

 */

double
gsl_ran_gumbel1 (const gsl_rng * r, const double a, const double b)
{
  double x = gsl_rng_uniform_pos (r);

  double z = (log(b) - log(-log(x))) / a;

  return z;
}

double
gsl_ran_gumbel1_pdf (const double x, const double a, const double b)
{
  double p = a * b *  exp (-(b * exp(-a * x) + a * x));
  return p;
}

double
gsl_ran_gumbel2 (const gsl_rng * r, const double a, const double b)
{
  double x = gsl_rng_uniform_pos (r);

  double z = pow(-b / log(x), 1/a);

  return z;
}

double
gsl_ran_gumbel2_pdf (const double x, const double a, const double b)
{
  if (x <= 0)
    {
      return 0 ;
    }
  else
    {
      double p = b * a *  pow(x,-(a+1)) * exp (-b * pow(x, -a));
      return p;
    }
}

/* The hypergeometric distribution has the form,

   prob(k) =  choose(n1,t) choose(n2, t-k) / choose(n1+n2,t)

   where choose(a,b) = a!/(b!(a-b)!) 

   n1 + n2 is the total population (tagged plus untagged)
   n1 is the tagged population
   t is the number of samples taken (without replacement)
   k is the number of tagged samples found

*/

unsigned int
gsl_ran_hypergeometric (const gsl_rng * r, unsigned int n1, unsigned int n2, 
                        unsigned int t)
{
  const unsigned int n = n1 + n2;

  unsigned int i = 0;
  unsigned int a = n1;
  unsigned int b = n1 + n2;
  unsigned int k = 0;

  if (t > n)
    {
      t = n ;
    }

  if (t < n / 2) 
    {
      for (i = 0 ; i < t ; i++)
        {
          double u = gsl_rng_uniform(r) ;
          
          if (b * u < a)
            {
              k++ ;
              if (k == n1)
                return k ;
              a-- ;
            }
          b-- ;
        }
      return k;
    }
  else
    {
      for (i = 0 ; i < n - t ; i++)
        {
          double u = gsl_rng_uniform(r) ;
          
          if (b * u < a)
            {
              k++ ;
              if (k == n1)
                return n1 - k ;
              a-- ;
            }
          b-- ;
        }
      return n1 - k;
    }


}

double
gsl_ran_hypergeometric_pdf (const unsigned int k, 
                            const unsigned int n1, 
                            const unsigned int n2, 
                            unsigned int t)
{
  if (t > n1 + n2)
    {
      t = n1 + n2 ;
    }

  if (k > n1 || k > t)
    {
      return 0 ;
    }
  else if (t > n2 && k + n2 < t )
    {
      return 0 ;
    }
  else 
    {
      double p;
      
      double c1 = gsl_sf_lnchoose(n1,k);
      double c2 = gsl_sf_lnchoose(n2,t-k);
      double c3 = gsl_sf_lnchoose(n1+n2,t);

      p = exp(c1 + c2 - c3) ;

      return p;
    }
}

double
gsl_ran_landau_pdf(const double x)
{
  static double P1[5] =
    {
      0.4259894875E0, -0.1249762550E0, 0.3984243700E-1,
      -0.6298287635E-2, 0.1511162253E-2
    };
  static double P2[5] =
    {
      0.1788541609E0, 0.1173957403E0, 0.1488850518E-1,
      -0.1394989411E-2, 0.1283617211E-3
    };
  static double P3[5] =
    {
      0.1788544503E0, 0.9359161662E-1, 0.6325387654E-2,
      0.6611667319E-4, -0.2031049101E-5
    };
  static double P4[5] =
    {
      0.9874054407E0, 0.1186723273E3, 0.8492794360E3,
      -0.7437792444E3, 0.4270262186E3
    };
  static double P5[5] =
    {
      0.1003675074E1, 0.1675702434E3, 0.4789711289E4,
      0.2121786767E5, -0.2232494910E5
    };
  static double P6[5] =
    {
      0.1000827619E1, 0.6649143136E3, 0.6297292665E5,
      0.4755546998E6, -0.5743609109E7
    };

  static double Q1[5] =
    {
      1.0, -0.3388260629E0, 0.9594393323E-1,
      -0.1608042283E-1, 0.3778942063E-2
    };
  static double Q2[5] =
    {
      1.0, 0.7428795082E0, 0.3153932961E0,
      0.6694219548E-1, 0.8790609714E-2
    };
  static double Q3[5] =
    {
      1.0, 0.6097809921E0, 0.2560616665E0,
      0.4746722384E-1, 0.6957301675E-2
    };
  static double Q4[5] =
    {
      1.0, 0.1068615961E3, 0.3376496214E3,
      0.2016712389E4, 0.1597063511E4
    };
  static double Q5[5] =
    {
      1.0, 0.1569424537E3, 0.3745310488E4,
      0.9834698876E4, 0.6692428357E5
    };
  static double Q6[5] =
    {
      1.0, 0.6514101098E3, 0.5697473333E5,
      0.1659174725E6, -0.2815759939E7
    };

  static double A1[3] =
    {
      0.4166666667E-1, -0.1996527778E-1, 0.2709538966E-1
    };
  static double A2[2] =
    {
      -0.1845568670E1, -0.4284640743E1
    };

  double U, V, DENLAN;

  V = x;
  if (V < -5.5)
    {
      U = exp(V + 1.0);
      DENLAN = 0.3989422803 * (exp( -1 / U) / sqrt(U)) *
        (1 + (A1[0] + (A1[1] + A1[2] * U) * U) * U);
    }
  else if (V < -1)
    {
      U = exp( -V - 1);
      DENLAN = exp( -U) * sqrt(U) *
        (P1[0] + (P1[1] + (P1[2] + (P1[3] + P1[4] * V) * V) * V) * V) /
        (Q1[0] + (Q1[1] + (Q1[2] + (Q1[3] + Q1[4] * V) * V) * V) * V);
    }
  else if (V < 1)
    {
      DENLAN = (P2[0] + (P2[1] + (P2[2] + (P2[3] + P2[4] * V) * V) * V) * V) /
        (Q2[0] + (Q2[1] + (Q2[2] + (Q2[3] + Q2[4] * V) * V) * V) * V);
    }
  else if (V < 5)
    {
      DENLAN = (P3[0] + (P3[1] + (P3[2] + (P3[3] + P3[4] * V) * V) * V) * V) /
        (Q3[0] + (Q3[1] + (Q3[2] + (Q3[3] + Q3[4] * V) * V) * V) * V);
    }
  else if (V < 12)
    {
      U = 1 / V;
      DENLAN = U * U *
        (P4[0] + (P4[1] + (P4[2] + (P4[3] + P4[4] * U) * U) * U) * U) /
        (Q4[0] + (Q4[1] + (Q4[2] + (Q4[3] + Q4[4] * U) * U) * U) * U);
    }
  else if (V < 50)
    {
      U = 1 / V;
      DENLAN = U * U *
        (P5[0] + (P5[1] + (P5[2] + (P5[3] + P5[4] * U) * U) * U) * U) /
        (Q5[0] + (Q5[1] + (Q5[2] + (Q5[3] + Q5[4] * U) * U) * U) * U);
    }
  else if (V < 300)
    {
      U = 1 / V;
      DENLAN = U * U *
        (P6[0] + (P6[1] + (P6[2] + (P6[3] + P6[4] * U) * U) * U) * U) /
        (Q6[0] + (Q6[1] + (Q6[2] + (Q6[3] + Q6[4] * U) * U) * U) * U);
    }
  else
    {
      U = 1 / (V - V * log(V) / (V + 1));
      DENLAN = U * U * (1 + (A2[0] + A2[1] * U) * U);
    }

  return DENLAN;
}

#if 0 /* Not needed yet */
/* This function is a translation from the original Fortran of the
 * CERN library routine DISLAN, the integral from -inf to x of the
 * Landau p.d.f.
 */
static
double
gsl_ran_landau_dislan(const double x)
{
  static double P1[5] =
    {
      0.2514091491E0, -0.6250580444E-1,
      0.1458381230E-1, -0.2108817737E-2,
      0.7411247290E-3
    };

  static double P2[4] =
    {
      0.2868328584E0, 0.3564363231E0,
      0.1523518695E0, 0.2251304883E-1
    };

  static double P3[4] =
    {
      0.2868329066E0, 0.3003828436E0,
      0.9950951941E-1, 0.8733827185E-2
    };

  static double P4[4] =
    {
      0.1000351630E1, 0.4503592498E1,
      0.1085883880E2, 0.7536052269E1
    };

  static double P5[4] =
    {
      0.1000006517E1, 0.4909414111E2,
      0.8505544753E2, 0.1532153455E3
    };

  static double P6[4] =
    {
      0.1000000983E1, 0.1329868456E3,
      0.9162149244E3, -0.9605054274E3
    };

  static double Q1[5] =
    {
      1.0, -0.5571175625E-2,
      0.6225310236E-1, -0.3137378427E-2,
      0.1931496439E-2
    };

  static double Q2[4] =
    {
      1.0, 0.6191136137E0,
      0.1720721448E0, 0.2278594771E-1
    };

  static double Q3[4] =
    {
      1.0, 0.4237190502E0,
      0.1095631512E0, 0.8693851567E-2
    };

  static double Q4[4] =
    {
      1.0, 0.5539969678E1,
      0.1933581111E2, 0.2721321508E2
    };

  static double Q5[4] =
    {
      1.0, 0.5009928881E2,
      0.1399819104E3, 0.4200002909E3
    };

  static double Q6[4] =
    {
      1.0, 0.1339887843E3,
      0.1055990413E4, 0.5532224619E3
    };

  static double A1[3] =
    {
      -0.4583333333E0, 0.6675347222E0, -0.1641741416E1
    };

  static double A2[3] =
    {
      1.0, -0.4227843351E0, -0.2043403138E1
    };

  double U, V, DISLAN;

  V = x;
  if (V < -5.5)
    {
      U = exp(V + 1);
      DISLAN = 0.3989422803 * exp( -1 / U) * sqrt(U) *
               (1 + (A1[0] + (A1[1] + A1[2] * U) * U) * U);
    }
  else if (V < -1)
    {
      U = exp( -V - 1);
      DISLAN = (exp( -U) / sqrt(U)) *
               (P1[0] + (P1[1] + (P1[2] + (P1[3] + P1[4] * V) * V) * V) * V) /
               (Q1[0] + (Q1[1] + (Q1[2] + (Q1[3] + Q1[4] * V) * V) * V) * V);
    }
  else if (V < 1)
    {
      DISLAN = (P2[0] + (P2[1] + (P2[2] + P2[3] * V) * V) * V) /
               (Q2[0] + (Q2[1] + (Q2[2] + Q2[3] * V) * V) * V);
    }
  else if (V < 4)
    {
      DISLAN = (P3[0] + (P3[1] + (P3[2] + P3[3] * V) * V) * V) /
               (Q3[0] + (Q3[1] + (Q3[2] + Q3[3] * V) * V) * V);
    }
  else if (V < 12)
    {
      U = 1 / V;
      DISLAN = (P4[0] + (P4[1] + (P4[2] + P4[3] * U) * U) * U) /
               (Q4[0] + (Q4[1] + (Q4[2] + Q4[3] * U) * U) * U);
    }
  else if (V < 50)
    {
      U = 1 / V;
      DISLAN = (P5[0] + (P5[1] + (P5[2] + P5[3] * U) * U) * U) /
               (Q5[0] + (Q5[1] + (Q5[2] + Q5[3] * U) * U) * U);
    }
  else if (V < 300)
    {
      U = 1 / V;
      DISLAN = (P6[0] + (P6[1] + (P6[2] + P6[3] * U) * U) * U) /
               (Q6[0] + (Q6[1] + (Q6[2] + Q6[3] * U) * U) * U);
    }
  else
    {
      U = 1 / (V - V * log(V) / (V + 1));
      DISLAN = 1 - (A2[0] + (A2[1] + A2[2] * U) * U) * U;
    }

  return DISLAN;
}
#endif

double
gsl_ran_landau(const gsl_rng * r)
{
  static double F[983] =
    {
      0.0000000,   /* Add empty element [0] to account for difference 
                      between C and Fortran convention for lower bound. */
      00.000000, 00.000000, 00.000000, 00.000000, 00.000000,
      -2.244733, -2.204365, -2.168163, -2.135219, -2.104898,
      -2.076740, -2.050397, -2.025605, -2.002150, -1.979866,
      -1.958612, -1.938275, -1.918760, -1.899984, -1.881879,
      -1.864385, -1.847451, -1.831030, -1.815083, -1.799574,
      -1.784473, -1.769751, -1.755383, -1.741346, -1.727620,
      -1.714187, -1.701029, -1.688130, -1.675477, -1.663057,
      -1.650858, -1.638868, -1.627078, -1.615477, -1.604058,
      -1.592811, -1.581729, -1.570806, -1.560034, -1.549407,
      -1.538919, -1.528565, -1.518339, -1.508237, -1.498254,
      -1.488386, -1.478628, -1.468976, -1.459428, -1.449979,
      -1.440626, -1.431365, -1.422195, -1.413111, -1.404112,
      -1.395194, -1.386356, -1.377594, -1.368906, -1.360291,
      -1.351746, -1.343269, -1.334859, -1.326512, -1.318229,
      -1.310006, -1.301843, -1.293737, -1.285688, -1.277693,
      -1.269752, -1.261863, -1.254024, -1.246235, -1.238494,
      -1.230800, -1.223153, -1.215550, -1.207990, -1.200474,
      -1.192999, -1.185566, -1.178172, -1.170817, -1.163500,
      -1.156220, -1.148977, -1.141770, -1.134598, -1.127459,
      -1.120354, -1.113282, -1.106242, -1.099233, -1.092255,
      -1.085306, -1.078388, -1.071498, -1.064636, -1.057802,
      -1.050996, -1.044215, -1.037461, -1.030733, -1.024029,
      -1.017350, -1.010695, -1.004064, -0.997456, -0.990871,
      -0.984308, -0.977767, -0.971247, -0.964749, -0.958271,
      -0.951813, -0.945375, -0.938957, -0.932558, -0.926178,
      -0.919816, -0.913472, -0.907146, -0.900838, -0.894547,
      -0.888272, -0.882014, -0.875773, -0.869547, -0.863337,
      -0.857142, -0.850963, -0.844798, -0.838648, -0.832512,
      -0.826390, -0.820282, -0.814187, -0.808106, -0.802038,
      -0.795982, -0.789940, -0.783909, -0.777891, -0.771884,
      -0.765889, -0.759906, -0.753934, -0.747973, -0.742023,
      -0.736084, -0.730155, -0.724237, -0.718328, -0.712429,
      -0.706541, -0.700661, -0.694791, -0.688931, -0.683079,
      -0.677236, -0.671402, -0.665576, -0.659759, -0.653950,
      -0.648149, -0.642356, -0.636570, -0.630793, -0.625022,
      -0.619259, -0.613503, -0.607754, -0.602012, -0.596276,
      -0.590548, -0.584825, -0.579109, -0.573399, -0.567695,
      -0.561997, -0.556305, -0.550618, -0.544937, -0.539262,
      -0.533592, -0.527926, -0.522266, -0.516611, -0.510961,
      -0.505315, -0.499674, -0.494037, -0.488405, -0.482777,
      -0.477153, -0.471533, -0.465917, -0.460305, -0.454697,
      -0.449092, -0.443491, -0.437893, -0.432299, -0.426707,
      -0.421119, -0.415534, -0.409951, -0.404372, -0.398795,
      -0.393221, -0.387649, -0.382080, -0.376513, -0.370949,
      -0.365387, -0.359826, -0.354268, -0.348712, -0.343157,
      -0.337604, -0.332053, -0.326503, -0.320955, -0.315408,
      -0.309863, -0.304318, -0.298775, -0.293233, -0.287692,
      -0.282152, -0.276613, -0.271074, -0.265536, -0.259999,
      -0.254462, -0.248926, -0.243389, -0.237854, -0.232318,
      -0.226783, -0.221247, -0.215712, -0.210176, -0.204641,
      -0.199105, -0.193568, -0.188032, -0.182495, -0.176957,
      -0.171419, -0.165880, -0.160341, -0.154800, -0.149259,
      -0.143717, -0.138173, -0.132629, -0.127083, -0.121537,
      -0.115989, -0.110439, -0.104889, -0.099336, -0.093782,
      -0.088227, -0.082670, -0.077111, -0.071550, -0.065987,
      -0.060423, -0.054856, -0.049288, -0.043717, -0.038144,
      -0.032569, -0.026991, -0.021411, -0.015828, -0.010243,
      -0.004656, 00.000934, 00.006527, 00.012123, 00.017722,
      00.023323, 00.028928, 00.034535, 00.040146, 00.045759,
      00.051376, 00.056997, 00.062620, 00.068247, 00.073877,
      00.079511, 00.085149, 00.090790, 00.096435, 00.102083,
      00.107736, 00.113392, 00.119052, 00.124716, 00.130385,
      00.136057, 00.141734, 00.147414, 00.153100, 00.158789,
      00.164483, 00.170181, 00.175884, 00.181592, 00.187304,
      00.193021, 00.198743, 00.204469, 00.210201, 00.215937,
      00.221678, 00.227425, 00.233177, 00.238933, 00.244696,
      00.250463, 00.256236, 00.262014, 00.267798, 00.273587,
      00.279382, 00.285183, 00.290989, 00.296801, 00.302619,
      00.308443, 00.314273, 00.320109, 00.325951, 00.331799,
      00.337654, 00.343515, 00.349382, 00.355255, 00.361135,
      00.367022, 00.372915, 00.378815, 00.384721, 00.390634,
      00.396554, 00.402481, 00.408415, 00.414356, 00.420304,
      00.426260, 00.432222, 00.438192, 00.444169, 00.450153,
      00.456145, 00.462144, 00.468151, 00.474166, 00.480188,
      00.486218, 00.492256, 00.498302, 00.504356, 00.510418,
      00.516488, 00.522566, 00.528653, 00.534747, 00.540850,
      00.546962, 00.553082, 00.559210, 00.565347, 00.571493,
      00.577648, 00.583811, 00.589983, 00.596164, 00.602355,
      00.608554, 00.614762, 00.620980, 00.627207, 00.633444,
      00.639689, 00.645945, 00.652210, 00.658484, 00.664768,
      00.671062, 00.677366, 00.683680, 00.690004, 00.696338,
      00.702682, 00.709036, 00.715400, 00.721775, 00.728160,
      00.734556, 00.740963, 00.747379, 00.753807, 00.760246,
      00.766695, 00.773155, 00.779627, 00.786109, 00.792603,
      00.799107, 00.805624, 00.812151, 00.818690, 00.825241,
      00.831803, 00.838377, 00.844962, 00.851560, 00.858170,
      00.864791, 00.871425, 00.878071, 00.884729, 00.891399,
      00.898082, 00.904778, 00.911486, 00.918206, 00.924940,
      00.931686, 00.938446, 00.945218, 00.952003, 00.958802,
      00.965614, 00.972439, 00.979278, 00.986130, 00.992996,
      00.999875, 01.006769, 01.013676, 01.020597, 01.027533,
      01.034482, 01.041446, 01.048424, 01.055417, 01.062424,
      01.069446, 01.076482, 01.083534, 01.090600, 01.097681,
      01.104778, 01.111889, 01.119016, 01.126159, 01.133316,
      01.140490, 01.147679, 01.154884, 01.162105, 01.169342,
      01.176595, 01.183864, 01.191149, 01.198451, 01.205770,
      01.213105, 01.220457, 01.227826, 01.235211, 01.242614,
      01.250034, 01.257471, 01.264926, 01.272398, 01.279888,
      01.287395, 01.294921, 01.302464, 01.310026, 01.317605,
      01.325203, 01.332819, 01.340454, 01.348108, 01.355780,
      01.363472, 01.371182, 01.378912, 01.386660, 01.394429,
      01.402216, 01.410024, 01.417851, 01.425698, 01.433565,
      01.441453, 01.449360, 01.457288, 01.465237, 01.473206,
      01.481196, 01.489208, 01.497240, 01.505293, 01.513368,
      01.521465, 01.529583, 01.537723, 01.545885, 01.554068,
      01.562275, 01.570503, 01.578754, 01.587028, 01.595325,
      01.603644, 01.611987, 01.620353, 01.628743, 01.637156,
      01.645593, 01.654053, 01.662538, 01.671047, 01.679581,
      01.688139, 01.696721, 01.705329, 01.713961, 01.722619,
      01.731303, 01.740011, 01.748746, 01.757506, 01.766293,
      01.775106, 01.783945, 01.792810, 01.801703, 01.810623,
      01.819569, 01.828543, 01.837545, 01.846574, 01.855631,
      01.864717, 01.873830, 01.882972, 01.892143, 01.901343,
      01.910572, 01.919830, 01.929117, 01.938434, 01.947781,
      01.957158, 01.966566, 01.976004, 01.985473, 01.994972,
      02.004503, 02.014065, 02.023659, 02.033285, 02.042943,
      02.052633, 02.062355, 02.072110, 02.081899, 02.091720,
      02.101575, 02.111464, 02.121386, 02.131343, 02.141334,
      02.151360, 02.161421, 02.171517, 02.181648, 02.191815,
      02.202018, 02.212257, 02.222533, 02.232845, 02.243195,
      02.253582, 02.264006, 02.274468, 02.284968, 02.295507,
      02.306084, 02.316701, 02.327356, 02.338051, 02.348786,
      02.359562, 02.370377, 02.381234, 02.392131, 02.403070,
      02.414051, 02.425073, 02.436138, 02.447246, 02.458397,
      02.469591, 02.480828, 02.492110, 02.503436, 02.514807,
      02.526222, 02.537684, 02.549190, 02.560743, 02.572343,
      02.583989, 02.595682, 02.607423, 02.619212, 02.631050,
      02.642936, 02.654871, 02.666855, 02.678890, 02.690975,
      02.703110, 02.715297, 02.727535, 02.739825, 02.752168,
      02.764563, 02.777012, 02.789514, 02.802070, 02.814681,
      02.827347, 02.840069, 02.852846, 02.865680, 02.878570,
      02.891518, 02.904524, 02.917588, 02.930712, 02.943894,
      02.957136, 02.970439, 02.983802, 02.997227, 03.010714,
      03.024263, 03.037875, 03.051551, 03.065290, 03.079095,
      03.092965, 03.106900, 03.120902, 03.134971, 03.149107,
      03.163312, 03.177585, 03.191928, 03.206340, 03.220824,
      03.235378, 03.250005, 03.264704, 03.279477, 03.294323,
      03.309244, 03.324240, 03.339312, 03.354461, 03.369687,
      03.384992, 03.400375, 03.415838, 03.431381, 03.447005,
      03.462711, 03.478500, 03.494372, 03.510328, 03.526370,
      03.542497, 03.558711, 03.575012, 03.591402, 03.607881,
      03.624450, 03.641111, 03.657863, 03.674708, 03.691646,
      03.708680, 03.725809, 03.743034, 03.760357, 03.777779,
      03.795300, 03.812921, 03.830645, 03.848470, 03.866400,
      03.884434, 03.902574, 03.920821, 03.939176, 03.957640,
      03.976215, 03.994901, 04.013699, 04.032612, 04.051639,
      04.070783, 04.090045, 04.109425, 04.128925, 04.148547,
      04.168292, 04.188160, 04.208154, 04.228275, 04.248524,
      04.268903, 04.289413, 04.310056, 04.330832, 04.351745,
      04.372794, 04.393982, 04.415310, 04.436781, 04.458395,
      04.480154, 04.502060, 04.524114, 04.546319, 04.568676,
      04.591187, 04.613854, 04.636678, 04.659662, 04.682807,
      04.706116, 04.729590, 04.753231, 04.777041, 04.801024,
      04.825179, 04.849511, 04.874020, 04.898710, 04.923582,
      04.948639, 04.973883, 04.999316, 05.024942, 05.050761,
      05.076778, 05.102993, 05.129411, 05.156034, 05.182864,
      05.209903, 05.237156, 05.264625, 05.292312, 05.320220,
      05.348354, 05.376714, 05.405306, 05.434131, 05.463193,
      05.492496, 05.522042, 05.551836, 05.581880, 05.612178,
      05.642734, 05.673552, 05.704634, 05.735986, 05.767610,
      05.799512, 05.831694, 05.864161, 05.896918, 05.929968,
      05.963316, 05.996967, 06.030925, 06.065194, 06.099780,
      06.134687, 06.169921, 06.205486, 06.241387, 06.277630,
      06.314220, 06.351163, 06.388465, 06.426130, 06.464166,
      06.502578, 06.541371, 06.580553, 06.620130, 06.660109,
      06.700495, 06.741297, 06.782520, 06.824173, 06.866262,
      06.908795, 06.951780, 06.995225, 07.039137, 07.083525,
      07.128398, 07.173764, 07.219632, 07.266011, 07.312910,
      07.360339, 07.408308, 07.456827, 07.505905, 07.555554,
      07.605785, 07.656608, 07.708035, 07.760077, 07.812747,
      07.866057, 07.920019, 07.974647, 08.029953, 08.085952,
      08.142657, 08.200083, 08.258245, 08.317158, 08.376837,
      08.437300, 08.498562, 08.560641, 08.623554, 08.687319,
      08.751955, 08.817481, 08.883916, 08.951282, 09.019600,
      09.088889, 09.159174, 09.230477, 09.302822, 09.376233,
      09.450735, 09.526355, 09.603118, 09.681054, 09.760191,
      09.840558, 09.922186, 10.005107, 10.089353, 10.174959,
      10.261958, 10.350389, 10.440287, 10.531693, 10.624646,
      10.719188, 10.815362, 10.913214, 11.012789, 11.114137,
      11.217307, 11.322352, 11.429325, 11.538283, 11.649285,
      11.762390, 11.877664, 11.995170, 12.114979, 12.237161,
      12.361791, 12.488946, 12.618708, 12.751161, 12.886394,
      13.024498, 13.165570, 13.309711, 13.457026, 13.607625,
      13.761625, 13.919145, 14.080314, 14.245263, 14.414134,
      14.587072, 14.764233, 14.945778, 15.131877, 15.322712,
      15.518470, 15.719353, 15.925570, 16.137345, 16.354912,
      16.578520, 16.808433, 17.044929, 17.288305, 17.538873,
      17.796967, 18.062943, 18.337176, 18.620068, 18.912049,
      19.213574, 19.525133, 19.847249, 20.180480, 20.525429,
      20.882738, 21.253102, 21.637266, 22.036036, 22.450278,
      22.880933, 23.329017, 23.795634, 24.281981, 24.789364,
      25.319207, 25.873062, 26.452634, 27.059789, 27.696581,
      28.365274, 29.068370, 29.808638, 30.589157, 31.413354,
      32.285060, 33.208568, 34.188705, 35.230920, 36.341388,
      37.527131, 38.796172, 40.157721, 41.622399, 43.202525,
      44.912465, 46.769077, 48.792279, 51.005773, 53.437996,
      56.123356, 59.103894
    };
  double X, U, V, RANLAN;
  int I;

  X = gsl_rng_uniform_pos(r);
  U = 1000.0 * X;
  I = U;
  U = U - I;

  if (I >= 70 && I <= 800)
    {
      RANLAN = F[I] + U * (F[I + 1] - F[I]);
    }
  else if (I >= 7 && I <= 980)
    {
      RANLAN = F[I] 
        + U * (F[I + 1] - F[I] 
               - 0.25 * (1 - U) * (F[I + 2] - F[I + 1] - F[I] + F[I - 1]));
    }
  else if (I < 7)
    {
      V = log(X);
      U = 1 / V;
      RANLAN = ((0.99858950 + (3.45213058E1 + 1.70854528E1 * U) * U) /
                (1 + (3.41760202E1 + 4.01244582 * U) * U)) *
               ( -log( -0.91893853 - V) - 1);
    }
  else
    {
      U = 1 - X;
      V = U * U;
      if (X <= 0.999)
        {
          RANLAN = (1.00060006 + 2.63991156E2 * U + 4.37320068E3 * V) /
                   ((1 + 2.57368075E2 * U + 3.41448018E3 * V) * U);
        }
      else
        {
          RANLAN = (1.00001538 + 6.07514119E3 * U + 7.34266409E5 * V) /
                   ((1 + 6.06511919E3 * U + 6.94021044E5 * V) * U);
        }
    }

  return RANLAN;
}

/* The stable Levy probability distributions have the form

   p(x) dx = (1/(2 pi)) \int dt exp(- it x - |c t|^alpha)

   with 0 < alpha <= 2. 

   For alpha = 1, we get the Cauchy distribution
   For alpha = 2, we get the Gaussian distribution with sigma = sqrt(2) c.

   Fromn Chapter 5 of Bratley, Fox and Schrage "A Guide to
   Simulation". The original reference given there is,

   J.M. Chambers, C.L. Mallows and B. W. Stuck. "A method for
   simulating stable random variates". Journal of the American
   Statistical Association, JASA 71 340-344 (1976).

   */

double
gsl_ran_levy (const gsl_rng * r, const double c, const double alpha)
{
  double u, v, t, s;

  u = M_PI * (gsl_rng_uniform_pos (r) - 0.5);

  if (alpha == 1)               /* cauchy case */
    {
      t = tan (u);
      return c * t;
    }

  do
    {
      v = gsl_ran_exponential (r, 1.0);
    }
  while (v == 0);

  if (alpha == 2)             /* gaussian case */
    {
      t = 2 * sin (u) * sqrt(v);
      return c * t;
    }

  /* general case */

  t = sin (alpha * u) / pow (cos (u), 1 / alpha);
  s = pow (cos ((1 - alpha) * u) / v, (1 - alpha) / alpha);

  return c * t * s;
}


/* The following routine for the skew-symmetric case was provided by
   Keith Briggs.

   The stable Levy probability distributions have the form

   2*pi* p(x) dx

     = \int dt exp(mu*i*t-|sigma*t|^alpha*(1-i*beta*sign(t)*tan(pi*alpha/2))) for alpha!=1
     = \int dt exp(mu*i*t-|sigma*t|^alpha*(1+i*beta*sign(t)*2/pi*log(|t|)))   for alpha==1

   with 0<alpha<=2, -1<=beta<=1, sigma>0.

   For beta=0, sigma=c, mu=0, we get gsl_ran_levy above.

   For alpha = 1, beta=0, we get the Lorentz distribution
   For alpha = 2, beta=0, we get the Gaussian distribution

   See A. Weron and R. Weron: Computer simulation of Lévy alpha-stable 
   variables and processes, preprint Technical University of Wroclaw.
   http://www.im.pwr.wroc.pl/~hugo/Publications.html

*/

double
gsl_ran_levy_skew (const gsl_rng * r, const double c, 
                   const double alpha, const double beta)
{
  double V, W, X;

  if (beta == 0)  /* symmetric case */
    {
      return gsl_ran_levy (r, c, alpha);
    }

  V = M_PI * (gsl_rng_uniform_pos (r) - 0.5);

  do
    {
      W = gsl_ran_exponential (r, 1.0);
    }
  while (W == 0);

  if (alpha == 1)
    {
      X = ((M_PI_2 + beta * V) * tan (V) -
           beta * log (M_PI_2 * W * cos (V) / (M_PI_2 + beta * V))) / M_PI_2;
      return c * (X + beta * log (c) / M_PI_2);
    }
  else
    {
      double t = beta * tan (M_PI_2 * alpha);
      double B = atan (t) / alpha;
      double S = pow (1 + t * t, 1/(2 * alpha));

      X = S * sin (alpha * (V + B)) / pow (cos (V), 1 / alpha)
        * pow (cos (V - alpha * (V + B)) / W, (1 - alpha) / alpha);
      return c * X;
    }
}

/* Logarithmic distribution 

   prob(n) =   p^n / (n log(1/(1-p)) for n = 1, 2, 3, ...

   We use Kemp's second accelerated generator, from Luc Devroye's book
   on "Non-Uniform Random Variate Generation", Springer */

unsigned int
gsl_ran_logarithmic (const gsl_rng * r, const double p)
{
  double c = log (1-p) ;

  double v = gsl_rng_uniform_pos (r);
  
  if (v >= p)
    {
      return 1 ;
    }
  else
    {
      double u = gsl_rng_uniform_pos (r);      
      double q = 1 - exp (c * u);

      if (v <= q*q)
        {
          double x = 1 + log(v)/log(q) ;
          return x ;
        }
      else if (v <= q)
        {
          return 2;
        }
      else
        {
          return 1 ;
        }
    }
}

double
gsl_ran_logarithmic_pdf (const unsigned int k, const double p)
{
  if (k == 0)
    {
      return 0 ;
    }
  else 
    {
      double P = pow(p, (double)k) / (double) k / log(1/(1-p)) ;
      return P;
    }
}

/* The logistic distribution has the form,

   p(x) dx = (1/a) exp(-x/a) / (1 + exp(-x/a))^2 dx

   for -infty < x < infty */

double
gsl_ran_logistic (const gsl_rng * r, const double a)
{
  double x, z;

  do
    {
      x = gsl_rng_uniform_pos (r);
    }
  while (x == 1);

  z = log (x / (1 - x));

  return a * z;
}

double
gsl_ran_logistic_pdf (const double x, const double a)
{
  double u = exp (-fabs(x)/a);
  double p = u / (fabs(a) * (1 + u) * (1 + u));
  return p;
}

/* The lognormal distribution has the form 

   p(x) dx = 1/(x * sqrt(2 pi sigma^2)) exp(-(ln(x) - zeta)^2/2 sigma^2) dx

   for x > 0. Lognormal random numbers are the exponentials of
   gaussian random numbers */

double
gsl_ran_lognormal (const gsl_rng * r, const double zeta, const double sigma)
{
  double u, v, r2, normal, z;

  do
    {
      /* choose x,y in uniform square (-1,-1) to (+1,+1) */

      u = -1 + 2 * gsl_rng_uniform (r);
      v = -1 + 2 * gsl_rng_uniform (r);

      /* see if it is in the unit circle */
      r2 = u * u + v * v;
    }
  while (r2 > 1.0 || r2 == 0);

  normal = u * sqrt (-2.0 * log (r2) / r2);

  z =  exp (sigma * normal + zeta);

  return z;
}

double
gsl_ran_lognormal_pdf (const double x, const double zeta, const double sigma)
{
  if (x <= 0)
    {
      return 0 ;
    }
  else
    {
      double u = (log (x) - zeta)/sigma;
      double p = 1 / (x * fabs(sigma) * sqrt (2 * M_PI)) * exp (-(u * u) /2);
      return p;
    }
}

/* The multinomial distribution has the form

                                      N!           n_1  n_2      n_K
   prob(n_1, n_2, ... n_K) = -------------------- p_1  p_2  ... p_K
                             (n_1! n_2! ... n_K!) 

   where n_1, n_2, ... n_K are nonnegative integers, sum_{k=1,K} n_k = N,
   and p = (p_1, p_2, ..., p_K) is a probability distribution. 

   Random variates are generated using the conditional binomial method.
   This scales well with N and does not require a setup step.

   Ref: 
   C.S. David, The computer generation of multinomial random variates,
   Comp. Stat. Data Anal. 16 (1993) 205-217
*/

void
gsl_ran_multinomial (const gsl_rng * r, const size_t K,
                     const unsigned int N, const double p[], unsigned int n[])
{
  size_t k;
  double norm = 0.0;
  double sum_p = 0.0;

  unsigned int sum_n = 0;

  /* p[k] may contain non-negative weights that do not sum to 1.0.
   * Even a probability distribution will not exactly sum to 1.0
   * due to rounding errors. 
   */

  for (k = 0; k < K; k++)
    {
      norm += p[k];
    }

  for (k = 0; k < K; k++)
    {
      if (p[k] > 0.0)
        {
          n[k] = gsl_ran_binomial (r, p[k] / (norm - sum_p), N - sum_n);
        }
      else
        {
          n[k] = 0;
        }

      sum_p += p[k];
      sum_n += n[k];
    }

}


double
gsl_ran_multinomial_pdf (const size_t K,
                         const double p[], const unsigned int n[])
{
  return exp (gsl_ran_multinomial_lnpdf (K, p, n));
}


double
gsl_ran_multinomial_lnpdf (const size_t K,
                           const double p[], const unsigned int n[])
{
  size_t k;
  unsigned int N = 0;
  double log_pdf = 0.0;
  double norm = 0.0;

  for (k = 0; k < K; k++)
    {
      N += n[k];
    }

  for (k = 0; k < K; k++)
    {
      norm += p[k];
    }

  log_pdf = gsl_sf_lnfact (N);

  for (k = 0; k < K; k++)
    {
      /* Handle case where n[k]==0 and p[k]==0 */

      if (n[k] > 0) 
        {
          log_pdf += log (p[k] / norm) * n[k] - gsl_sf_lnfact (n[k]);
        }
    }

  return log_pdf;
}

/* The binomial distribution has the form,

   f(x) =  n!/(x!(n-x)!) * p^x (1-p)^(n-x) for integer 0 <= x <= n
        =  0                               otherwise

   This implementation follows the public domain ranlib function
   "ignbin", the bulk of which is the BTPE (Binomial Triangle
   Parallelogram Exponential) algorithm introduced in
   Kachitvichyanukul and Schmeiser[1].  It has been translated to use
   modern C coding standards.

   If n is small and/or p is near 0 or near 1 (specifically, if
   n*min(p,1-p) < SMALL_MEAN), then a different algorithm, called
   BINV, is used which has an average runtime that scales linearly
   with n*min(p,1-p).

   But for larger problems, the BTPE algorithm takes the form of two
   functions b(x) and t(x) -- "bottom" and "top" -- for which b(x) <
   f(x)/f(M) < t(x), with M = floor(n*p+p).  b(x) defines a triangular
   region, and t(x) includes a parallelogram and two tails.  Details
   (including a nice drawing) are in the paper.

   [1] Kachitvichyanukul, V. and Schmeiser, B. W.  Binomial Random
   Variate Generation.  Communications of the ACM, 31, 2 (February,
   1988) 216.

   Note, Bruce Schmeiser (personal communication) points out that if
   you want very fast binomial deviates, and you are happy with
   approximate results, and/or n and n*p are both large, then you can
   just use gaussian estimates: mean=n*p, variance=n*p*(1-p).

   This implementation by James Theiler, April 2003, after obtaining
   permission -- and some good advice -- from Drs. Kachitvichyanukul
   and Schmeiser to use their code as a starting point, and then doing
   a little bit of tweaking.

   Additional polishing for GSL coding standards by Brian Gough.  */

#define SMALL_MEAN 14           /* If n*p < SMALL_MEAN then use BINV
                                   algorithm. The ranlib
                                   implementation used cutoff=30; but
                                   on my computer 14 works better */

#define BINV_CUTOFF 110         /* In BINV, do not permit ix too large */

#define FAR_FROM_MEAN 20        /* If ix-n*p is larger than this, then
                                   use the "squeeze" algorithm.
                                   Ranlib used 20, and this seems to
                                   be the best choice on my machine as
                                   well */

#define LNFACT(x) gsl_sf_lnfact(x)

inline static double
Stirling (double y1)
{
  double y2 = y1 * y1;
  double s =
    (13860.0 -
     (462.0 - (132.0 - (99.0 - 140.0 / y2) / y2) / y2) / y2) / y1 / 166320.0;
  return s;
}

unsigned int
gsl_ran_binomial_tpe (const gsl_rng * rng, double p, unsigned int n)
{
  return gsl_ran_binomial (rng, p, n);
}

unsigned int
gsl_ran_binomial (const gsl_rng * rng, double p, unsigned int n)
{
  int ix;                       /* return value */
  int flipped = 0;
  double q, s, np;

  if (n == 0)
    return 0;

  if (p > 0.5)
    {
      p = 1.0 - p;              /* work with small p */
      flipped = 1;
    }

  q = 1 - p;
  s = p / q;
  np = n * p;

  /* Inverse cdf logic for small mean (BINV in K+S) */

  if (np < SMALL_MEAN)
    {
      double f0 = gsl_pow_uint (q, n);   /* f(x), starting with x=0 */

      while (1)
        {
          /* This while(1) loop will almost certainly only loop once; but
           * if u=1 to within a few epsilons of machine precision, then it
           * is possible for roundoff to prevent the main loop over ix to
           * achieve its proper value.  following the ranlib implementation,
           * we introduce a check for that situation, and when it occurs,
           * we just try again.
           */

          double f = f0;
          double u = gsl_rng_uniform (rng);

          for (ix = 0; ix <= BINV_CUTOFF; ++ix)
            {
              if (u < f)
                goto Finish;
              u -= f;
              /* Use recursion f(x+1) = f(x)*[(n-x)/(x+1)]*[p/(1-p)] */
              f *= s * (n - ix) / (ix + 1);
            }

          /* It should be the case that the 'goto Finish' was encountered
           * before this point was ever reached.  But if we have reached
           * this point, then roundoff has prevented u from decreasing
           * all the way to zero.  This can happen only if the initial u
           * was very nearly equal to 1, which is a rare situation.  In
           * that rare situation, we just try again.
           *
           * Note, following the ranlib implementation, we loop ix only to
           * a hardcoded value of SMALL_MEAN_LARGE_N=110; we could have
           * looped to n, and 99.99...% of the time it won't matter.  This
           * choice, I think is a little more robust against the rare
           * roundoff error.  If n>LARGE_N, then it is technically
           * possible for ix>LARGE_N, but it is astronomically rare, and
           * if ix is that large, it is more likely due to roundoff than
           * probability, so better to nip it at LARGE_N than to take a
           * chance that roundoff will somehow conspire to produce an even
           * larger (and more improbable) ix.  If n<LARGE_N, then once
           * ix=n, f=0, and the loop will continue until ix=LARGE_N.
           */
        }
    }
  else
    {
      /* For n >= SMALL_MEAN, we invoke the BTPE algorithm */

      int k;

      double ffm = np + p;      /* ffm = n*p+p             */
      int m = (int) ffm;        /* m = int floor[n*p+p]    */
      double fm = m;            /* fm = double m;          */
      double xm = fm + 0.5;     /* xm = half integer mean (tip of triangle)  */
      double npq = np * q;      /* npq = n*p*q            */

      /* Compute cumulative area of tri, para, exp tails */

      /* p1: radius of triangle region; since height=1, also: area of region */
      /* p2: p1 + area of parallelogram region */
      /* p3: p2 + area of left tail */
      /* p4: p3 + area of right tail */
      /* pi/p4: probability of i'th area (i=1,2,3,4) */

      /* Note: magic numbers 2.195, 4.6, 0.134, 20.5, 15.3 */
      /* These magic numbers are not adjustable...at least not easily! */

      double p1 = floor (2.195 * sqrt (npq) - 4.6 * q) + 0.5;

      /* xl, xr: left and right edges of triangle */
      double xl = xm - p1;
      double xr = xm + p1;

      /* Parameter of exponential tails */
      /* Left tail:  t(x) = c*exp(-lambda_l*[xl - (x+0.5)]) */
      /* Right tail: t(x) = c*exp(-lambda_r*[(x+0.5) - xr]) */

      double c = 0.134 + 20.5 / (15.3 + fm);
      double p2 = p1 * (1.0 + c + c);

      double al = (ffm - xl) / (ffm - xl * p);
      double lambda_l = al * (1.0 + 0.5 * al);
      double ar = (xr - ffm) / (xr * q);
      double lambda_r = ar * (1.0 + 0.5 * ar);
      double p3 = p2 + c / lambda_l;
      double p4 = p3 + c / lambda_r;

      double var, accept;
      double u, v;              /* random variates */

    TryAgain:

      /* generate random variates, u specifies which region: Tri, Par, Tail */
      u = gsl_rng_uniform (rng) * p4;
      v = gsl_rng_uniform (rng);

      if (u <= p1)
        {
          /* Triangular region */
          ix = (int) (xm - p1 * v + u);
          goto Finish;
        }
      else if (u <= p2)
        {
          /* Parallelogram region */
          double x = xl + (u - p1) / c;
          v = v * c + 1.0 - fabs (x - xm) / p1;
          if (v > 1.0 || v <= 0.0)
            goto TryAgain;
          ix = (int) x;
        }
      else if (u <= p3)
        {
          /* Left tail */
          ix = (int) (xl + log (v) / lambda_l);
          if (ix < 0)
            goto TryAgain;
          v *= ((u - p2) * lambda_l);
        }
      else
        {
          /* Right tail */
          ix = (int) (xr - log (v) / lambda_r);
          if (ix > (double) n)
            goto TryAgain;
          v *= ((u - p3) * lambda_r);
        }

      /* At this point, the goal is to test whether v <= f(x)/f(m) 
       *
       *  v <= f(x)/f(m) = (m!(n-m)! / (x!(n-x)!)) * (p/q)^{x-m}
       *
       */

      /* Here is a direct test using logarithms.  It is a little
       * slower than the various "squeezing" computations below, but
       * if things are working, it should give exactly the same answer
       * (given the same random number seed).  */

#ifdef DIRECT
      var = log (v);

      accept =
        LNFACT (m) + LNFACT (n - m) - LNFACT (ix) - LNFACT (n - ix)
        + (ix - m) * log (p / q);

#else /* SQUEEZE METHOD */

      /* More efficient determination of whether v < f(x)/f(M) */

      k = abs (ix - m);

      if (k <= FAR_FROM_MEAN)
        {
          /* 
           * If ix near m (ie, |ix-m|<FAR_FROM_MEAN), then do
           * explicit evaluation using recursion relation for f(x)
           */
          double g = (n + 1) * s;
          double f = 1.0;

          var = v;

          if (m < ix)
            {
              int i;
              for (i = m + 1; i <= ix; i++)
                {
                  f *= (g / i - s);
                }
            }
          else if (m > ix)
            {
              int i;
              for (i = ix + 1; i <= m; i++)
                {
                  f /= (g / i - s);
                }
            }

          accept = f;
        }
      else
        {
          /* If ix is far from the mean m: k=ABS(ix-m) large */

          var = log (v);

          if (k < npq / 2 - 1)
            {
              /* "Squeeze" using upper and lower bounds on
               * log(f(x)) The squeeze condition was derived
               * under the condition k < npq/2-1 */
              double amaxp =
                k / npq * ((k * (k / 3.0 + 0.625) + (1.0 / 6.0)) / npq + 0.5);
              double ynorm = -(k * k / (2.0 * npq));
              if (var < ynorm - amaxp)
                goto Finish;
              if (var > ynorm + amaxp)
                goto TryAgain;
            }

          /* Now, again: do the test log(v) vs. log f(x)/f(M) */

#if USE_EXACT
          /* This is equivalent to the above, but is a little (~20%) slower */
          /* There are five log's vs three above, maybe that's it? */

          accept = LNFACT (m) + LNFACT (n - m)
            - LNFACT (ix) - LNFACT (n - ix) + (ix - m) * log (p / q);

#else /* USE STIRLING */
          /* The "#define Stirling" above corresponds to the first five
           * terms in asymptoic formula for
           * log Gamma (y) - (y-0.5)log(y) + y - 0.5 log(2*pi);
           * See Abramowitz and Stegun, eq 6.1.40
           */

          /* Note below: two Stirling's are added, and two are
           * subtracted.  In both K+S, and in the ranlib
           * implementation, all four are added.  I (jt) believe that
           * is a mistake -- this has been confirmed by personal
           * correspondence w/ Dr. Kachitvichyanukul.  Note, however,
           * the corrections are so small, that I couldn't find an
           * example where it made a difference that could be
           * observed, let alone tested.  In fact, define'ing Stirling
           * to be zero gave identical results!!  In practice, alv is
           * O(1), ranging 0 to -10 or so, while the Stirling
           * correction is typically O(10^{-5}) ...setting the
           * correction to zero gives about a 2% performance boost;
           * might as well keep it just to be pendantic.  */

          {
            double x1 = ix + 1.0;
            double w1 = n - ix + 1.0;
            double f1 = fm + 1.0;
            double z1 = n + 1.0 - fm;

            accept = xm * log (f1 / x1) + (n - m + 0.5) * log (z1 / w1)
              + (ix - m) * log (w1 * p / (x1 * q))
              + Stirling (f1) + Stirling (z1) - Stirling (x1) - Stirling (w1);
          }
#endif
#endif
        }


      if (var <= accept)
        {
          goto Finish;
        }
      else
        {
          goto TryAgain;
        }
    }

Finish:

  return (flipped) ? (n - ix) : (unsigned int)ix;
}

double gsl_pow_int(double x, int n)
{
  unsigned int un;

  if(n < 0) {
    x = 1.0/x;
    un = -n;
  } else {
    un = n;
  }

  return gsl_pow_uint(x, un);
}

double gsl_pow_uint(double x, unsigned int n)
{
  double value = 1.0;

  /* repeated squaring method 
   * returns 0.0^0 = 1.0, so continuous in x
   */
  do {
     if(n & 1) value *= x;  /* for n odd */
     n >>= 1;
     x *= x;
  } while (n);

  return value;
}

/* The negative binomial distribution has the form,

   prob(k) =  Gamma(n + k)/(Gamma(n) Gamma(k + 1))  p^n (1-p)^k 

   for k = 0, 1, ... . Note that n does not have to be an integer.

   This is the Leger's algorithm (given in the answers in Knuth) */

unsigned int
gsl_ran_negative_binomial (const gsl_rng * r, double p, double n)
{
  double X = gsl_ran_gamma (r, n, 1.0) ;
  unsigned int k = gsl_ran_poisson (r, X*(1-p)/p) ;
  return k ;
}

double
gsl_ran_negative_binomial_pdf (const unsigned int k, const double p, double n)
{
  double P;

  double f = gsl_sf_lngamma (k + n) ;
  double a = gsl_sf_lngamma (n) ;
  double b = gsl_sf_lngamma (k + 1.0) ;

  P = exp(f - a - b + n * log(p) + k * log1p(-p));

  return P;
}

/* The poisson distribution has the form

   p(n) = (mu^n / n!) exp(-mu) 

   for n = 0, 1, 2, ... . The method used here is the one from Knuth. */

unsigned int
gsl_ran_poisson (const gsl_rng * r, double mu)
{
  double emu;
  double prod = 1.0;
  unsigned int k = 0;

  while (mu > 10)
    {
      unsigned int m = mu * (7.0 / 8.0);

      double X = gsl_ran_gamma_int (r, m);

      if (X >= mu)
        {
          return k + gsl_ran_binomial (r, mu / X, m - 1);
        }
      else
        {
          k += m;
          mu -= X; 
        }
    }

  /* This following method works well when mu is small */

  emu = exp (-mu);

  do
    {
      prod *= gsl_rng_uniform (r);
      k++;
    }
  while (prod > emu);

  return k - 1;

}

void
gsl_ran_poisson_array (const gsl_rng * r, size_t n, unsigned int array[],
                       double mu)
{
  size_t i;

  for (i = 0; i < n; i++)
    {
      array[i] = gsl_ran_poisson (r, mu);
    }

  return;
}

double
gsl_ran_poisson_pdf (const unsigned int k, const double mu)
{
  double p;
  double lf = gsl_sf_lnfact (k); 

  p = exp (log (mu) * k - lf - mu);
  return p;
}

/* The Pareto distribution has the form,

   p(x) dx = (a/b) / (x/b)^(a+1) dx     for x >= b

 */

double
gsl_ran_pareto (const gsl_rng * r, double a, const double b)
{
  double x = gsl_rng_uniform_pos (r);

  double z = pow (x, -1 / a);

  return b * z;
}

double
gsl_ran_pareto_pdf (const double x, const double a, const double b)
{
  if (x >= b)
    {
      double p = (a/b) / pow (x/b, a + 1);
      return p;
    }
  else
    {
      return 0;
    }
}

/* The Pascal distribution is a negative binomial with valued integer n

   prob(k) =  (n - 1 + k)!/(n!(k - 1)!) *  p^n (1-p)^k for k = 0, 1, ..., n

   */

unsigned int
gsl_ran_pascal (const gsl_rng * r, double p, unsigned int n)
{
  /* This is a separate interface for the pascal distribution so that
     it can be optimized differently from the negative binomial in
     future.
     
     e.g. if n < 10 it might be faster to generate the Pascal
     distributions as the sum of geometric variates directly.  */
  
  unsigned int k = gsl_ran_negative_binomial (r, p, (double) n);
  return k;
}

double
gsl_ran_pascal_pdf (const unsigned int k, const double p, unsigned int n)
{
  double P = gsl_ran_negative_binomial_pdf (k, p, (double) n);
  return P;
}

/* The Rayleigh distribution has the form

   p(x) dx = (x / sigma^2) exp(-x^2/(2 sigma^2)) dx

   for x = 0 ... +infty */

double
gsl_ran_rayleigh (const gsl_rng * r, const double sigma)
{
  double u = gsl_rng_uniform_pos (r);

  return sigma * sqrt(-2.0 * log (u));
}

double
gsl_ran_rayleigh_pdf (const double x, const double sigma)
{
  if (x < 0)
    {
      return 0 ;
    }
  else
    {
      double u = x / sigma ;
      double p = (u / sigma) * exp(-u * u / 2.0) ;
      
      return p;
    }
}

/* The Rayleigh tail distribution has the form

   p(x) dx = (x / sigma^2) exp((a^2 - x^2)/(2 sigma^2)) dx

   for x = a ... +infty */

double
gsl_ran_rayleigh_tail (const gsl_rng * r, const double a, const double sigma)
{
  double u = gsl_rng_uniform_pos (r);

  return sqrt(a * a - 2.0 * sigma * sigma * log (u));
}

double
gsl_ran_rayleigh_tail_pdf (const double x, const double a, const double sigma)
{
  if (x < a)
    {
      return 0 ;
    }
  else
    {
      double u = x / sigma ;
      double v = a / sigma ;

      double p = (u / sigma) * exp((v + u) * (v - u) / 2.0) ;
      
      return p;
    }
}

void
gsl_ran_dir_2d (const gsl_rng * r, double *x, double *y)
{
  /* This method avoids trig, but it does take an average of 8/pi =
   * 2.55 calls to the RNG, instead of one for the direct
   * trigonometric method.  */

  double u, v, s;
  do
    {
      u = -1 + 2 * gsl_rng_uniform (r);
      v = -1 + 2 * gsl_rng_uniform (r);
      s = u * u + v * v;
    }
  while (s > 1.0 || s == 0.0);

  /* This is the Von Neumann trick. See Knuth, v2, 3rd ed, p140
   * (exercise 23).  Note, no sin, cos, or sqrt !  */

  *x = (u * u - v * v) / s;
  *y = 2 * u * v / s;

  /* Here is the more straightforward approach, 
   *     s = sqrt (s);  *x = u / s;  *y = v / s;
   * It has fewer total operations, but one of them is a sqrt */
}

void
gsl_ran_dir_2d_trig_method (const gsl_rng * r, double *x, double *y)
{
  /* This is the obvious solution... */
  /* It ain't clever, but since sin/cos are often hardware accelerated,
   * it can be faster -- it is on my home Pentium -- than von Neumann's
   * solution, or slower -- as it is on my Sun Sparc 20 at work
   */
  double t = 6.2831853071795864 * gsl_rng_uniform (r);          /* 2*PI */
  *x = cos (t);
  *y = sin (t);
}

void
gsl_ran_dir_3d (const gsl_rng * r, double *x, double *y, double *z)
{
  double s, a;

  /* This is a variant of the algorithm for computing a random point
   * on the unit sphere; the algorithm is suggested in Knuth, v2,
   * 3rd ed, p136; and attributed to Robert E Knop, CACM, 13 (1970),
   * 326.
   */

  /* Begin with the polar method for getting x,y inside a unit circle
   */
  do
    {
      *x = -1 + 2 * gsl_rng_uniform (r);
      *y = -1 + 2 * gsl_rng_uniform (r);
      s = (*x) * (*x) + (*y) * (*y);
    }
  while (s > 1.0);

  *z = -1 + 2 * s;              /* z uniformly distributed from -1 to 1 */
  a = 2 * sqrt (1 - s);         /* factor to adjust x,y so that x^2+y^2
                                 * is equal to 1-z^2 */
  *x *= a;
  *y *= a;
}

void
gsl_ran_dir_nd (const gsl_rng * r, size_t n, double *x)
{
  double d;
  size_t i;
  /* See Knuth, v2, 3rd ed, p135-136.  The method is attributed to
   * G. W. Brown, in Modern Mathematics for the Engineer (1956).
   * The idea is that gaussians G(x) have the property that
   * G(x)G(y)G(z)G(...) is radially symmetric, a function only
   * r = sqrt(x^2+y^2+...)
   */
  d = 0;
  do
    {
      for (i = 0; i < n; ++i)
        {
          x[i] = gsl_ran_gaussian (r, 1.0);
          d += x[i] * x[i];
        }
    }
  while (d == 0);
  d = sqrt (d);
  for (i = 0; i < n; ++i)
    {
      x[i] /= d;
    }
}

/* The t-distribution has the form

   p(x) dx = (Gamma((nu + 1)/2)/(sqrt(pi nu) Gamma(nu/2))
   * (1 + (x^2)/nu)^-((nu + 1)/2) dx

   The method used here is the one described in Knuth */

double
gsl_ran_tdist (const gsl_rng * r, const double nu)
{
  if (nu <= 2)
    {
      double Y1 = gsl_ran_ugaussian (r);
      double Y2 = gsl_ran_chisq (r, nu);

      double t = Y1 / sqrt (Y2 / nu);

      return t;
    }
  else
    {
      double Y1, Y2, Z, t;
      do
        {
          Y1 = gsl_ran_ugaussian (r);
          Y2 = gsl_ran_exponential (r, 1 / (nu/2 - 1));

          Z = Y1 * Y1 / (nu - 2);
        }
      while (1 - Z < 0 || exp (-Y2 - Z) > (1 - Z));

      /* Note that there is a typo in Knuth's formula, the line below
         is taken from the original paper of Marsaglia, Mathematics of
         Computation, 34 (1980), p 234-256 */

      t = Y1 / sqrt ((1 - 2 / nu) * (1 - Z));
      return t;
    }
}

double
gsl_ran_tdist_pdf (const double x, const double nu)
{
  double p;

  double lg1 = gsl_sf_lngamma (nu / 2);
  double lg2 = gsl_sf_lngamma ((nu + 1) / 2);

  p = ((exp (lg2 - lg1) / sqrt (M_PI * nu)) 
       * pow ((1 + x * x / nu), -(nu + 1) / 2));
  return p;
}

/* The Weibull distribution has the form,

   p(x) dx = (b/a) (x/a)^(b-1) exp(-(x/a)^b) dx

 */

double
gsl_ran_weibull (const gsl_rng * r, const double a, const double b)
{
  double x = gsl_rng_uniform_pos (r);

  double z = pow (-log (x), 1 / b);

  return a * z;
}

double
gsl_ran_weibull_pdf (const double x, const double a, const double b)
{
  if (x < 0)
    {
      return 0 ;
    }
  else if (x == 0)
    {
      if (b == 1)
        return 1/a ;
      else
        return 0 ;
    }
  else if (b == 1)
    {
      return exp(-x/a)/a ;
    }
  else
    {
      double p = (b/a) * exp (-pow (x/a, b) + (b - 1) * log (x/a));
      return p;
    }
}

/* Inline swap and copy functions for moving objects around */

static inline 
void swap (void * base, size_t size, size_t i, size_t j)
{
  register char * a = size * i + (char *) base ;
  register char * b = size * j + (char *) base ;
  register size_t s = size ;

  if (i == j)
    return ;
  
  do                                            
    {                                           
      char tmp = *a;                            
      *a++ = *b;                                
      *b++ = tmp;                               
    } 
  while (--s > 0);                              
}

static inline void 
copy (void * dest, size_t i, void * src, size_t j, size_t size)
{
  register char * a = size * i + (char *) dest ;
  register char * b = size * j + (char *) src ;
  register size_t s = size ;
  
  do                                            
    {                                           
      *a++ = *b++;                              
    } 
  while (--s > 0);                              
}

/* Randomly permute (shuffle) N indices

   Supply an array x[N] with nmemb members, each of size size and on
   return it will be shuffled into a random order.  The algorithm is
   from Knuth, SemiNumerical Algorithms, v2, p139, who cites Moses and
   Oakford, and Durstenfeld */

void
gsl_ran_shuffle (const gsl_rng * r, void * base, size_t n, size_t size)
{
  size_t i ;

  for (i = n - 1; i > 0; i--)
    {
      size_t j = gsl_rng_uniform_int(r, i+1); /* originally (i + 1) * gsl_rng_uniform (r) */

      swap (base, size, i, j) ;
    }
}

int
gsl_ran_choose (const gsl_rng * r, void * dest, size_t k, void * src, 
                 size_t n, size_t size)
{
  size_t i, j = 0;

  /* Choose k out of n items, return an array x[] of the k items.
     These items will prevserve the relative order of the original
     input -- you can use shuffle() to randomize the output if you
     wish */

  if (k > n)
    {
      GSL_ERROR ("k is greater than n, cannot sample more than n items",
                 GSL_EINVAL) ;
    }

  for (i = 0; i < n && j < k; i++)
    {
      if ((n - i) * gsl_rng_uniform (r) < k - j)
        {
          copy (dest, j, src, i, size) ;
          j++ ;
        }
    }

  return GSL_SUCCESS;
}

void
gsl_ran_sample (const gsl_rng * r, void * dest, size_t k, void * src, 
                size_t n, size_t size)
{
  size_t i, j = 0;

  /* Choose k out of n items, with replacement */

  for (i = 0; i < k; i++)
    {
      j = gsl_rng_uniform_int (r, n);  /* originally n * gsl_rng_uniform (r) */

      copy (dest, i, src, j, size) ;
    }
}

