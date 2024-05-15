#include <common.h>
__import gretltypes;

static int gretl_digits = 6;
/* record of state, and communication of state with outside world */
static int compiling;    /* boolean: are we compiling a function currently? */
static int fn_executing; /* depth of function call stack */
static char gretl_errmsg[ERRLEN];
static char gretl_warnmsg[ERRLEN];
static int gretl_errno;
static int gretl_warnnum;
/* copy of a given model's missmask, or NULL */
static char *refmask;
static int blas_mnk_min = -1;
static int omp_n_threads;
/* Below: setting of the maximal value of K = the shared inner
   dimension in matrix multiplication for use of SIMD. Also
   setting of the minimum value of M x N for doing matrix
   addition and subtraction via SIMD. If these variables are
   set to -1 that disables SIMD by default (though the user
   can change that via the "set" command).
*/
static int simd_k_max = 8;   /* 2014-03-07: was -1 */
static int simd_mn_min = 16; /* 2014-03-07: was -1 */
static int error_printed;
static int alarm_set;
static struct global_vars_ {
    gint8 gretl_debug;
    gint8 gretl_assert;
    gint8 datacols;
    gint8 plot_collect;
    gint8 R_functions;
    gint8 R_lib;
    gint8 loglevel;
    gint8 logstamp;
    gint8 csv_digits;
    gint8 hac_missvals;
    int gmp_bits;
// This used to be globals but it conflicted
} globalvars = {0, 0, 5, 0, 0, 1, 2, 0, UNSET_INT, HAC_ES, 256};
/* the current set of state variables */
/* Find n independent "small" Mersenne Twisters with period 2^521-1;
   set the one corresponding to @self as the one to use
*/

static guint32 dcmt_seed;
static int use_dcmt = 0;
static int n_states;
static GPtrArray *state_stack;
static int state_idx = -1;
static int gretl_model_count = 0;
static const char *name_ok =
    "abcdefghijklmnopqrstuvwxyz"
    "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    "0123456789_";
static int omp_mnk_min = -1;
static char blas_core[64];
static char blas_parallel[32];
static char blas_version[32];
static int blas_variant;

static void (*OB_set_num_threads) (int);
static int (*OB_get_num_threads) (void);
static void (*BLIS_set_num_threads) (int);
static int (*BLIS_get_num_threads) (void);
static void (*BLIS_init) (void);
static void (*BLIS_finalize) (void);
static void (*MKL_finalize) (void);
static void (*MKL_domain_set_num_threads) (int, int);
static int (*MKL_domain_get_max_threads) (int);
static void (*FLAME_init) (void);
static int (*FLAME_initialized) (void);
static void (*FLAME_finalize) (void);
typedef struct flag_match_ {
	gretlopt o;
	char c;
} flag_match;
flag_match flag_matches[] = {
    { OPT_A, 'a' },
    { OPT_B, 'b' },
    { OPT_C, 'c' },
    { OPT_D, 'd' },
    { OPT_E, 'e' },
    { OPT_F, 'f' },
    { OPT_G, 'g' },
    { OPT_H, 'h' },
    { OPT_I, 'i' },
    { OPT_J, 'j' },
    { OPT_K, 'k' },
    { OPT_L, 'l' },
    { OPT_M, 'm' },
    { OPT_N, 'n' },
    { OPT_O, 'o' },
    { OPT_P, 'p' },
    { OPT_Q, 'q' },
    { OPT_R, 'r' },
    { OPT_S, 's' },
    { OPT_T, 't' },
    { OPT_U, 'u' },
    { OPT_V, 'v' },
    { OPT_W, 'w' },
    { OPT_X, 'x' },
    { OPT_Y, 'y' },
    { OPT_Z, 'z' },
    { 0L,   '\0' }
};
double default_nls_toler;
/* for testing purposes */
int stata_sa;
int IGLS;
char glsmat[MAXLEN];
int n_stored_opts;
/* Central accounting for error in matrix allocation */
int gretl_matrix_err;
/* the current set of state variables */
set_state *state;
/* An efficient means of allocating temporary storage for lapack
   operations: this should be used _only_ for temporary allocations
   that would ordinarily be freed before returning from the function
   in question.  In this mode we keep the chunk around for future use,
   expanding it as needed.
*/
void *lapack_mem_chunk;
size_t lapack_mem_sz;
double eq_tol = DEFAULT_EQTOL;
stored_opt * optinfo;

