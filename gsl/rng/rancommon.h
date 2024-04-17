#include <math.h>
#include <limits.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define INLINE_DECL inline  /* use C99 inline */
#define INLINE_FUN inline
#define GSL_VAR extern
#define GSL_ERROR_VAL(reason, gsl_errno, value) \
       do { \
       gsl_error (reason, __FILE__, __LINE__, gsl_errno) ; \
       return value ; \
       } while (0)
/* GSL_ERROR: call the error handler, and return the error code */

#define GSL_ERROR(reason, gsl_errno) \
       do { \
       gsl_error (reason, __FILE__, __LINE__, gsl_errno) ; \
       return gsl_errno ; \
       } while (0)
#define RETURN_IF_NULL(x) if (!x) { return ; }

#define GSL_DBL_EPSILON        2.2204460492503131e-16

/* The maximum x such that gamma(x) is not
 * considered an overflow.
 */
#define GSL_SF_GAMMA_XMAX  171.0

/* The maximum n such that gsl_sf_fact(n) does not give an overflow. */
#define GSL_SF_FACT_NMAX 170

/* The maximum n such that gsl_sf_doublefact(n) does not give an overflow. */
#define GSL_SF_DOUBLEFACT_NMAX 297

#define LogRootTwoPi_  0.9189385332046727418

#define GSL_SIGN(x)    ((x) >= 0.0 ? 1 : -1)

#define GSL_IS_ODD(n)  ((n) & 1)
#define GSL_IS_EVEN(n) (!(GSL_IS_ODD(n)))
#define GSL_MAX(a,b) ((a) > (b) ? (a) : (b))
#define GSL_MIN(a,b) ((a) < (b) ? (a) : (b))
#define GSL_MAX_INT(a,b)   GSL_MAX(a,b)
#define GSL_MIN_INT(a,b)   GSL_MIN(a,b)
#define GSL_MAX_DBL(a,b)   GSL_MAX(a,b)
#define GSL_MIN_DBL(a,b)   GSL_MIN(a,b)
#define GSL_MAX_LDBL(a,b)  GSL_MAX(a,b)
#define GSL_MIN_LDBL(a,b)  GSL_MIN(a,b)
#define M_SQRTPI   1.77245385090551602729816748334      /* sqrt(pi) */
#define OVERFLOW_ERROR(result) do { (result)->val = GSL_POSINF; (result)->err = GSL_POSINF; GSL_ERROR ("overflow", GSL_EOVRFLW); } while(0)

#define UNDERFLOW_ERROR(result) do { (result)->val = 0.0; (result)->err = GSL_DBL_MIN; GSL_ERROR ("underflow", GSL_EUNDRFLW); } while(0)

#define INTERNAL_OVERFLOW_ERROR(result) do { (result)->val = GSL_POSINF; (result)->err = GSL_POSINF; return GSL_EOVRFLW; } while(0)

#define INTERNAL_UNDERFLOW_ERROR(result) do { (result)->val = 0.0; (result)->err = GSL_DBL_MIN; return GSL_EUNDRFLW; } while(0)

#define DOMAIN_ERROR(result) do { (result)->val = GSL_NAN; (result)->err = GSL_NAN; GSL_ERROR ("domain error", GSL_EDOM); } while(0)

#define DOMAIN_ERROR_MSG(msg, result) do { (result)->val = GSL_NAN; (result)->err = GSL_NAN; GSL_ERROR ((msg), GSL_EDOM); } while(0)

#define DOMAIN_ERROR_E10(result) do { (result)->val = GSL_NAN; (result)->err = GSL_NAN; (result)->e10 = 0 ; GSL_ERROR ("domain error", GSL_EDOM); } while(0)

#define OVERFLOW_ERROR_E10(result) do { (result)->val = GSL_POSINF; (result)->err = GSL_POSINF; (result)->e10 = 0; GSL_ERROR ("overflow", GSL_EOVRFLW); } while(0)
#define UNDERFLOW_ERROR_E10(result) do { (result)->val = 0.0; (result)->err = GSL_DBL_MIN; (result)->e10 = 0; GSL_ERROR ("underflow", GSL_EUNDRFLW); } while(0)


#define OVERFLOW_ERROR_2(r1,r2) do { (r1)->val = GSL_POSINF; (r1)->err = GSL_POSINF; (r2)->val = GSL_POSINF ; (r2)->err=GSL_POSINF; GSL_ERROR ("overflow", GSL_EOVRFLW); } while(0)

#define UNDERFLOW_ERROR_2(r1,r2) do { (r1)->val = 0.0; (r1)->err = GSL_DBL_MIN; (r2)->val = 0.0 ; (r2)->err = GSL_DBL_MIN; GSL_ERROR ("underflow", GSL_EUNDRFLW); } while(0)

#define DOMAIN_ERROR_2(r1,r2) do { (r1)->val = GSL_NAN; (r1)->err = GSL_NAN;  (r2)->val = GSL_NAN; (r2)->err = GSL_NAN;  GSL_ERROR ("domain error", GSL_EDOM); } while(0)

/* evaluate a function discarding the status value in a modifiable way */

#define EVAL_RESULT(fn) \
   gsl_sf_result result; \
   int status = fn; \
   if (status != GSL_SUCCESS) { \
     GSL_ERROR_VAL(#fn, status, result.val); \
   } ; \
   return result.val;

#define EVAL_DOUBLE(fn) \
   int status = fn; \
   if (status != GSL_SUCCESS) { \
     GSL_ERROR_VAL(#fn, status, result); \
   } ; \
   return result;

#define MAXITER_ERROR(result) do { (result)->val = GSL_NAN; (result)->err = GSL_NAN; GSL_ERROR ("too many iterations error", GSL_EMAXITER); } while(0)
# define GSL_POSINF INFINITY
# define GSL_NEGINF (-INFINITY)
# define GSL_NAN NAN
#define M_LNPI     1.14472988584940017414342735135      /* ln(pi) */
#define GSL_DBL_MAX        1.7976931348623157e+308
#define GSL_DBL_MIN        2.2250738585072014e-308
#define GSL_ERROR_SELECT_2(a,b)       ((a) != GSL_SUCCESS ? (a) : ((b) != GSL_SUCCESS ? (b) : GSL_SUCCESS))
#define GSL_ERROR_SELECT_3(a,b,c)     ((a) != GSL_SUCCESS ? (a) : GSL_ERROR_SELECT_2(b,c))
#define GSL_ROOT4_DBL_EPSILON  1.2207031250000000e-04
#define CHECK_UNDERFLOW(r) if (fabs((r)->val) < GSL_DBL_MIN) GSL_ERROR("underflow", GSL_EUNDRFLW);
#define GSL_LOG_DBL_MIN   (-7.0839641853226408e+02)
#define GSL_LOG_DBL_MAX    7.0978271289338397e+02
#define M_EULER    0.57721566490153286060651209008      /* Euler constant */
#define GSL_SQRT_DBL_MAX   1.3407807929942596e+154
#define GSL_SQRT_DBL_MIN   1.4916681462400413e-154

typedef void gsl_stream_handler_t (const char * label, const char * file,
                                   int line, const char * reason);
typedef void gsl_error_handler_t (const char * reason, const char * file,
                                  int line, int gsl_errno);

typedef struct
  {
    const char *name;
    unsigned long int max;
    unsigned long int min;
    size_t size;
    void (*set) (void *state, unsigned long int seed);
    unsigned long int (*get) (void *state);
    double (*get_double) (void *state);
  }
gsl_rng_type;

typedef struct
  {
    const gsl_rng_type * type;
    void *state;
  }
gsl_rng;

enum { 
  GSL_SUCCESS  = 0, 
  GSL_FAILURE  = -1,
  GSL_CONTINUE = -2,  /* iteration has not converged */
  GSL_EDOM     = 1,   /* input domain error, e.g sqrt(-1) */
  GSL_ERANGE   = 2,   /* output range error, e.g. exp(1e100) */
  GSL_EFAULT   = 3,   /* invalid pointer */
  GSL_EINVAL   = 4,   /* invalid argument supplied by user */
  GSL_EFAILED  = 5,   /* generic failure */
  GSL_EFACTOR  = 6,   /* factorization failed */
  GSL_ESANITY  = 7,   /* sanity check failed - shouldn't happen */
  GSL_ENOMEM   = 8,   /* malloc failed */
  GSL_EBADFUNC = 9,   /* problem with user-supplied function */
  GSL_ERUNAWAY = 10,  /* iterative process is out of control */
  GSL_EMAXITER = 11,  /* exceeded max number of iterations */
  GSL_EZERODIV = 12,  /* tried to divide by zero */
  GSL_EBADTOL  = 13,  /* user specified an invalid tolerance */
  GSL_ETOL     = 14,  /* failed to reach the specified tolerance */
  GSL_EUNDRFLW = 15,  /* underflow */
  GSL_EOVRFLW  = 16,  /* overflow  */
  GSL_ELOSS    = 17,  /* loss of accuracy */
  GSL_EROUND   = 18,  /* failed because of roundoff error */
  GSL_EBADLEN  = 19,  /* matrix, vector lengths are not conformant */
  GSL_ENOTSQR  = 20,  /* matrix not square */
  GSL_ESING    = 21,  /* apparent singularity detected */
  GSL_EDIVERGE = 22,  /* integral or series is divergent */
  GSL_EUNSUP   = 23,  /* requested feature is not supported by the hardware */
  GSL_EUNIMPL  = 24,  /* requested feature not (yet) implemented */
  GSL_ECACHE   = 25,  /* cache limit exceeded */
  GSL_ETABLE   = 26,  /* table limit exceeded */
  GSL_ENOPROG  = 27,  /* iteration is not making progress towards solution */
  GSL_ENOPROGJ = 28,  /* jacobian evaluations are not improving the solution */
  GSL_ETOLF    = 29,  /* cannot reach the specified tolerance in F */
  GSL_ETOLX    = 30,  /* cannot reach the specified tolerance in X */
  GSL_ETOLG    = 31,  /* cannot reach the specified tolerance in gradient */
  GSL_EOF      = 32   /* end of file */
} ;

double gsl_ran_gamma(const gsl_rng * r, const double a, const double b);
double gsl_ran_beta(const gsl_rng * r, const double a, const double b);




