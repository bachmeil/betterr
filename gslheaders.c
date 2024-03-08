#include <stddef.h>
#define GSL_VAR extern

double gsl_stats_mean(const double data[], size_t stride, size_t n);

typedef struct {
    const char *name;
    unsigned long int max;
    unsigned long int min;
    size_t size;
    void (*set) (void *state, unsigned long int seed);
    unsigned long int (*get) (void *state);
    double (*get_double) (void *state);
} gsl_rng_type;

typedef struct {
    const gsl_rng_type * type;
    void *state;
} gsl_rng;

typedef struct {
  size_t size;
  double * data;
} gsl_block;

typedef struct {
  size_t size;
  size_t stride;
  double * data;
  gsl_block * block;
  int owner;
} gsl_vector;

// Warning: Not compatible with R's matrix type
// Stores data in row-major order
typedef struct {
  size_t size1; // rows
  size_t size2; // columns
  size_t tda; // Number of columns of storage allocated
  double * data;
  gsl_block * block;
  int owner; // 0 if this is a slice, not freed
} gsl_matrix;


unsigned long int gsl_rng_get(const gsl_rng *r);
double gsl_rng_uniform(const gsl_rng *r);
double gsl_rng_uniform_pos(const gsl_rng *r);
unsigned long int gsl_rng_uniform_int(const gsl_rng *r, unsigned long int n);

unsigned int gsl_ran_bernoulli (const gsl_rng * r, double p);
double gsl_ran_bernoulli_pdf (const unsigned int k, double p);

double gsl_ran_beta (const gsl_rng * r, const double a, const double b);
double gsl_ran_beta_pdf (const double x, const double a, const double b);

unsigned int gsl_ran_binomial (const gsl_rng * r, double p, unsigned int n);
unsigned int gsl_ran_binomial_knuth (const gsl_rng * r, double p, unsigned int n);
unsigned int gsl_ran_binomial_tpe (const gsl_rng * r, double p, unsigned int n);
double gsl_ran_binomial_pdf (const unsigned int k, const double p, const unsigned int n);

double gsl_ran_exponential (const gsl_rng * r, const double mu);
double gsl_ran_exponential_pdf (const double x, const double mu);

double gsl_ran_exppow (const gsl_rng * r, const double a, const double b);
double gsl_ran_exppow_pdf (const double x, const double a, const double b);

double gsl_ran_cauchy (const gsl_rng * r, const double a);
double gsl_ran_cauchy_pdf (const double x, const double a);

double gsl_ran_chisq (const gsl_rng * r, const double nu);
double gsl_ran_chisq_pdf (const double x, const double nu);

void gsl_ran_dirichlet (const gsl_rng * r, const size_t K, const double alpha[], double theta[]);
double gsl_ran_dirichlet_pdf (const size_t K, const double alpha[], const double theta[]);
double gsl_ran_dirichlet_lnpdf (const size_t K, const double alpha[], const double theta[]);

double gsl_ran_erlang (const gsl_rng * r, const double a, const double n);
double gsl_ran_erlang_pdf (const double x, const double a, const double n);

double gsl_ran_fdist (const gsl_rng * r, const double nu1, const double nu2);
double gsl_ran_fdist_pdf (const double x, const double nu1, const double nu2);

double gsl_ran_flat (const gsl_rng * r, const double a, const double b);
double gsl_ran_flat_pdf (double x, const double a, const double b);

double gsl_ran_gamma (const gsl_rng * r, const double a, const double b);
double gsl_ran_gamma_int (const gsl_rng * r, const unsigned int a);
double gsl_ran_gamma_pdf (const double x, const double a, const double b);
double gsl_ran_gamma_mt (const gsl_rng * r, const double a, const double b);
double gsl_ran_gamma_knuth (const gsl_rng * r, const double a, const double b);

double gsl_ran_gaussian (const gsl_rng * r, const double sigma);
double gsl_ran_gaussian_ratio_method (const gsl_rng * r, const double sigma);
double gsl_ran_gaussian_ziggurat (const gsl_rng * r, const double sigma);
double gsl_ran_gaussian_pdf (const double x, const double sigma);

double gsl_ran_ugaussian (const gsl_rng * r);
double gsl_ran_ugaussian_ratio_method (const gsl_rng * r);
double gsl_ran_ugaussian_pdf (const double x);

double gsl_ran_gaussian_tail (const gsl_rng * r, const double a, const double sigma);
double gsl_ran_gaussian_tail_pdf (const double x, const double a, const double sigma);

double gsl_ran_ugaussian_tail (const gsl_rng * r, const double a);
double gsl_ran_ugaussian_tail_pdf (const double x, const double a);

void gsl_ran_bivariate_gaussian (const gsl_rng * r, double sigma_x, double sigma_y, double rho, double *x, double *y);
double gsl_ran_bivariate_gaussian_pdf (const double x, const double y, const double sigma_x, const double sigma_y, const double rho);

int gsl_ran_multivariate_gaussian (const gsl_rng * r, const gsl_vector * mu, const gsl_matrix * L, gsl_vector * result);
int gsl_ran_multivariate_gaussian_log_pdf (const gsl_vector * x,
                                           const gsl_vector * mu,
                                           const gsl_matrix * L,
                                           double * result,
                                           gsl_vector * work);
int gsl_ran_multivariate_gaussian_pdf (const gsl_vector * x,
                                       const gsl_vector * mu,
                                       const gsl_matrix * L,
                                       double * result,
                                       gsl_vector * work);
int gsl_ran_multivariate_gaussian_mean (const gsl_matrix * X, gsl_vector * mu_hat);
int gsl_ran_multivariate_gaussian_vcov (const gsl_matrix * X, gsl_matrix * sigma_hat);

int gsl_ran_wishart (const gsl_rng * r,
                     const double df,
                     const gsl_matrix * L,
                     gsl_matrix * result,
                     gsl_matrix * work);
int gsl_ran_wishart_log_pdf (const gsl_matrix * X,
                             const gsl_matrix * L_X,
                             const double df,
                             const gsl_matrix * L,
                             double * result,
                             gsl_matrix * work);
int gsl_ran_wishart_pdf (const gsl_matrix * X,
                         const gsl_matrix * L_X,
                         const double df,
                         const gsl_matrix * L,
                         double * result,
                         gsl_matrix * work);

double gsl_ran_landau (const gsl_rng * r);
double gsl_ran_landau_pdf (const double x);

unsigned int gsl_ran_geometric (const gsl_rng * r, const double p);
double gsl_ran_geometric_pdf (const unsigned int k, const double p);

unsigned int gsl_ran_hypergeometric (const gsl_rng * r, unsigned int n1, unsigned int n2, unsigned int t);
double gsl_ran_hypergeometric_pdf (const unsigned int k, const unsigned int n1, const unsigned int n2, unsigned int t);

double gsl_ran_gumbel1 (const gsl_rng * r, const double a, const double b);
double gsl_ran_gumbel1_pdf (const double x, const double a, const double b);

double gsl_ran_gumbel2 (const gsl_rng * r, const double a, const double b);
double gsl_ran_gumbel2_pdf (const double x, const double a, const double b);

double gsl_ran_logistic (const gsl_rng * r, const double a);
double gsl_ran_logistic_pdf (const double x, const double a);

double gsl_ran_lognormal (const gsl_rng * r, const double zeta, const double sigma);
double gsl_ran_lognormal_pdf (const double x, const double zeta, const double sigma);

unsigned int gsl_ran_logarithmic (const gsl_rng * r, const double p);
double gsl_ran_logarithmic_pdf (const unsigned int k, const double p);

void gsl_ran_multinomial (const gsl_rng * r, const size_t K,
                          const unsigned int N, const double p[],
                          unsigned int n[] );
double gsl_ran_multinomial_pdf (const size_t K,
                                const double p[], const unsigned int n[] );
double gsl_ran_multinomial_lnpdf (const size_t K,
                           const double p[], const unsigned int n[] );


unsigned int gsl_ran_negative_binomial (const gsl_rng * r, double p, double n);
double gsl_ran_negative_binomial_pdf (const unsigned int k, const double p, double n);

unsigned int gsl_ran_pascal (const gsl_rng * r, double p, unsigned int n);
double gsl_ran_pascal_pdf (const unsigned int k, const double p, unsigned int n);

double gsl_ran_pareto (const gsl_rng * r, double a, const double b);
double gsl_ran_pareto_pdf (const double x, const double a, const double b);

unsigned int gsl_ran_poisson (const gsl_rng * r, double mu);
void gsl_ran_poisson_array (const gsl_rng * r, size_t n, unsigned int array[],
                            double mu);
double gsl_ran_poisson_pdf (const unsigned int k, const double mu);

double gsl_ran_rayleigh (const gsl_rng * r, const double sigma);
double gsl_ran_rayleigh_pdf (const double x, const double sigma);

double gsl_ran_rayleigh_tail (const gsl_rng * r, const double a, const double sigma);
double gsl_ran_rayleigh_tail_pdf (const double x, const double a, const double sigma);

double gsl_ran_tdist (const gsl_rng * r, const double nu);
double gsl_ran_tdist_pdf (const double x, const double nu);

double gsl_ran_laplace (const gsl_rng * r, const double a);
double gsl_ran_laplace_pdf (const double x, const double a);

double gsl_ran_levy (const gsl_rng * r, const double c, const double alpha);
double gsl_ran_levy_skew (const gsl_rng * r, const double c, const double alpha, const double beta);

double gsl_ran_weibull (const gsl_rng * r, const double a, const double b);
double gsl_ran_weibull_pdf (const double x, const double a, const double b);

void gsl_ran_dir_2d (const gsl_rng * r, double * x, double * y);
void gsl_ran_dir_2d_trig_method (const gsl_rng * r, double * x, double * y);
void gsl_ran_dir_3d (const gsl_rng * r, double * x, double * y, double * z);
void gsl_ran_dir_nd (const gsl_rng * r, size_t n, double * x);

void gsl_ran_shuffle (const gsl_rng * r, void * base, size_t nmembm, size_t size);
int gsl_ran_choose (const gsl_rng * r, void * dest, size_t k, void * src, size_t n, size_t size) ;
void gsl_ran_sample (const gsl_rng * r, void * dest, size_t k, void * src, size_t n, size_t size) ;


typedef struct {                /* struct for Walker algorithm */
    size_t K;
    size_t *A;
    double *F;
} gsl_ran_discrete_t;

gsl_ran_discrete_t * gsl_ran_discrete_preproc (size_t K, const double *P);
void gsl_ran_discrete_free(gsl_ran_discrete_t *g);
size_t gsl_ran_discrete (const gsl_rng *r, const gsl_ran_discrete_t *g);
double gsl_ran_discrete_pdf (size_t k, const gsl_ran_discrete_t *g);

gsl_rng *gsl_rng_alloc (const gsl_rng_type * T);
void gsl_rng_free (gsl_rng * r);
void gsl_rng_set (const gsl_rng * r, unsigned long int seed);

GSL_VAR const gsl_rng_type *gsl_rng_borosh13;
GSL_VAR const gsl_rng_type *gsl_rng_coveyou;
GSL_VAR const gsl_rng_type *gsl_rng_cmrg;
GSL_VAR const gsl_rng_type *gsl_rng_fishman18;
GSL_VAR const gsl_rng_type *gsl_rng_fishman20;
GSL_VAR const gsl_rng_type *gsl_rng_fishman2x;
GSL_VAR const gsl_rng_type *gsl_rng_gfsr4;
GSL_VAR const gsl_rng_type *gsl_rng_knuthran;
GSL_VAR const gsl_rng_type *gsl_rng_knuthran2;
GSL_VAR const gsl_rng_type *gsl_rng_knuthran2002;
GSL_VAR const gsl_rng_type *gsl_rng_lecuyer21;
GSL_VAR const gsl_rng_type *gsl_rng_minstd;
GSL_VAR const gsl_rng_type *gsl_rng_mrg;
GSL_VAR const gsl_rng_type *gsl_rng_mt19937;
GSL_VAR const gsl_rng_type *gsl_rng_mt19937_1999;
GSL_VAR const gsl_rng_type *gsl_rng_mt19937_1998;
GSL_VAR const gsl_rng_type *gsl_rng_r250;
GSL_VAR const gsl_rng_type *gsl_rng_ran0;
GSL_VAR const gsl_rng_type *gsl_rng_ran1;
GSL_VAR const gsl_rng_type *gsl_rng_ran2;
GSL_VAR const gsl_rng_type *gsl_rng_ran3;
GSL_VAR const gsl_rng_type *gsl_rng_rand;
GSL_VAR const gsl_rng_type *gsl_rng_rand48;
GSL_VAR const gsl_rng_type *gsl_rng_random128_bsd;
GSL_VAR const gsl_rng_type *gsl_rng_random128_glibc2;
GSL_VAR const gsl_rng_type *gsl_rng_random128_libc5;
GSL_VAR const gsl_rng_type *gsl_rng_random256_bsd;
GSL_VAR const gsl_rng_type *gsl_rng_random256_glibc2;
GSL_VAR const gsl_rng_type *gsl_rng_random256_libc5;
GSL_VAR const gsl_rng_type *gsl_rng_random32_bsd;
GSL_VAR const gsl_rng_type *gsl_rng_random32_glibc2;
GSL_VAR const gsl_rng_type *gsl_rng_random32_libc5;
GSL_VAR const gsl_rng_type *gsl_rng_random64_bsd;
GSL_VAR const gsl_rng_type *gsl_rng_random64_glibc2;
GSL_VAR const gsl_rng_type *gsl_rng_random64_libc5;
GSL_VAR const gsl_rng_type *gsl_rng_random8_bsd;
GSL_VAR const gsl_rng_type *gsl_rng_random8_glibc2;
GSL_VAR const gsl_rng_type *gsl_rng_random8_libc5;
GSL_VAR const gsl_rng_type *gsl_rng_random_bsd;
GSL_VAR const gsl_rng_type *gsl_rng_random_glibc2;
GSL_VAR const gsl_rng_type *gsl_rng_random_libc5;
GSL_VAR const gsl_rng_type *gsl_rng_randu;
GSL_VAR const gsl_rng_type *gsl_rng_ranf;
GSL_VAR const gsl_rng_type *gsl_rng_ranlux;
GSL_VAR const gsl_rng_type *gsl_rng_ranlux389;
GSL_VAR const gsl_rng_type *gsl_rng_ranlxd1;
GSL_VAR const gsl_rng_type *gsl_rng_ranlxd2;
GSL_VAR const gsl_rng_type *gsl_rng_ranlxs0;
GSL_VAR const gsl_rng_type *gsl_rng_ranlxs1;
GSL_VAR const gsl_rng_type *gsl_rng_ranlxs2;
GSL_VAR const gsl_rng_type *gsl_rng_ranmar;
GSL_VAR const gsl_rng_type *gsl_rng_slatec;
GSL_VAR const gsl_rng_type *gsl_rng_taus;
GSL_VAR const gsl_rng_type *gsl_rng_taus2;
GSL_VAR const gsl_rng_type *gsl_rng_taus113;
GSL_VAR const gsl_rng_type *gsl_rng_transputer;
GSL_VAR const gsl_rng_type *gsl_rng_tt800;
GSL_VAR const gsl_rng_type *gsl_rng_uni;
GSL_VAR const gsl_rng_type *gsl_rng_uni32;
GSL_VAR const gsl_rng_type *gsl_rng_vax;
GSL_VAR const gsl_rng_type *gsl_rng_waterman14;
GSL_VAR const gsl_rng_type *gsl_rng_zuf;

int gsl_matrix_swap_rows(gsl_matrix *m, size_t i, size_t j);
void gsl_matrix_free(gsl_matrix *m);
