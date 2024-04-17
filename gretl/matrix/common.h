#include <complex.h>
#include <ctype.h>
#include <errno.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifdef unix
	#include <sys/param.h>
	#include <unistd.h>
#endif
#ifdef _WIN32
	#include <stdint.h>
#endif
#ifdef _WIN64
	#include <stdint.h>
#endif
#ifdef _WIN64
	#define __locale_data
#endif

#define GRETL_VERSION "2023c"
#define LIBGRETL_REVISION  0
#define _(String)  ((char *) String)
#define N_(String) String
#define MAXLINE 131072 /* maximum length of command line */
#define MAXLABEL   128 /* maximum length of descriptive labels for variables */
#define MAXLEN     512 /* max length of regular "long" strings */
#define MAXDISP     32 /* max length of "display names" for variables */
#define VNAMELEN    32 /* space allocated for var names (including termination) */
#define OBSLEN      16 /* space allocated for obs strings (including termination) */
#ifndef M_PI
# define M_PI 3.1415926535897932384626432
#endif
#ifndef M_2PI
# define M_2PI 6.2831853071795864769252864
#endif
#define SQRT_2_PI     2.506628274631000502415765284811
#define LN_2_PI       1.837877066409345483560659472811
#define LN_SQRT_2_PI  0.918938533204672741780329736406
#define PMAX_NOT_AVAILABLE 666
#define NADBL NAN
#define na(x) (isnan(x) || isinf(x))
typedef isinf __builtin_isinf_sign;
#define LISTSEP (-100)

/* numbers smaller than the given limit will print as zero */
#define screen_zero(x)  ((fabs(x) > 1.0e-13)? x : 0.0)
#define Z_COLS_BORROWED 2
#define dset_zcols_borrowed(d) (d->auxiliary == Z_COLS_BORROWED)
#define RESAMPLED ((char *) 0xdeadbeef)
#define free_datainfo(p) do { if (p != NULL) { clear_datainfo(p, 0); free(p); } \
                            } while (0);
#define libset_boolvar(k) (k < STATE_FLAG_MAX || k==R_FUNCTIONS || \
			   k==R_LIB || k==LOGSTAMP)
#define libset_double(k) (k > STATE_INT_MAX && k < STATE_FLOAT_MAX)
#define libset_int(k) ((k > STATE_FLAG_MAX && k < STATE_INT_MAX) || \
		       (k > STATE_VARS_MAX && k < NS_INT_MAX))

#define libset_small_int(k) (k < STATE_SMALL_INT_MAX || \
			     (k > STATE_VARS_MAX && k < NS_SMALL_INT_MAX))

#define coded_intvar(k) (k == GARCH_VCV || \
			 k == GARCH_ALT_VCV || \
			 k == ARMA_VCV || \
			 k == HAC_LAG || \
			 k == HAC_KERNEL || \
                         k == HC_VERSION || \
			 k == PANEL_ROBUST || \
			 k == USE_QR || \
			 k == VECM_NORM || \
			 k == GRETL_OPTIM || \
			 k == MAX_VERBOSE || \
			 k == WILDBOOT_DIST || \
			 k == QUANTILE_TYPE || \
			 k == GRETL_ASSERT || \
			 k == PLOT_COLLECT || \
			 k == LOGLEVEL || \
			 k == HAC_MISSVALS)
#define INTS_OFFSET (1 + log2(STATE_FLAG_MAX))
/* Comment on 'TINY': It's the minimum value for 'test' (see below)
   that libgretl's Cholesky decomposition routine will accept before
   rejecting a data matrix as too highly collinear.  If you set it too
   high, data sets for which Cholesky could produce reasonable
   estimates will be rejected.  If you set it too low (and 100 *
   DBL_EPSILON is definitely too low), gretl will produce more or less
   worthless coefficient estimates when given highly collinear data.
   Before making a permanent change to the value of TINY, check how
   gretl does on the NIST reference data sets for linear regression
   and ensure you're not getting any garbage results.  The current
   value enables us to get decent results on the NIST nonlinear
   regression test suite; it might be a bit too low for some contexts.
*/

#define TINY      8.0e-09 /* was 2.1e-09 (last changed 2007/01/20) */
#define SMALL     2.0e-08 /* threshold for printing a warning for collinearity */
#define YBARZERO  0.5e-14 /* threshold for treating mean of dependent
			     variable as effectively zero */
#define ESSZERO   1.0e-22 /* threshold for considering a tiny error-sum-of-
			     squares value to be effectively zero */
#define G_N_ELEMENTS(arr)		(sizeof (arr) / sizeof ((arr)[0]))
#define SVD_SMIN 1.0e-9
#define QR_RCOND_MIN  1.0e-14
#define QR_RCOND_WARN 1.0e-07
#define gretl_cmatrix_set(m,i,j,x) ((m)->z[(j)*(m)->rows+(i)]=x)

typedef struct gretl_matrix_ gretl_vector;

/**
 * gretl_matrix_get:
 * @m: matrix.
 * @i: row.
 * @j: column.
 *
 * Returns: the @i, @j element of @m.
 */

#define gretl_matrix_get(m,i,j) (m->val[(j)*m->rows+(i)])
#define gretl_cmatrix_get(m,i,j) (m->z[(j)*m->rows+(i)])
/**
 * gretl_vector_get:
 * @v: vector.
 * @i: index.
 *
 * Returns: element @i of @v.
 */

#define gretl_vector_get(v,i) (v->val[i])

/**
 * gretl_matrix_set:
 * @m: matrix.
 * @i: row.
 * @j: column.
 * @x: value to set.
 *
 * Sets the @i, @j element of @m to @x.
 */

#define gretl_matrix_set(m,i,j,x) ((m)->val[(j)*(m)->rows+(i)]=x)
#define gretl_cmatrix_set(m,i,j,x) ((m)->z[(j)*(m)->rows+(i)]=x)
/**
 * gretl_vector_set:
 * @v: vector.
 * @i: index.
 * @x: value to set.
 *
 * Sets element @i of @v to @x.
 */

#define gretl_vector_set(v,i,x) ((v)->val[i]=x)

/**
 * gretl_matrix_cols:
 * @m: matrix to query.
 *
 * Returns: the number of columns in @m.
 */

#define gretl_matrix_cols(m) ((m == NULL)? 0 : m->cols)

/**
 * gretl_matrix_rows:
 * @m: matrix to query.
 *
 * Returns: the number of rows in @m.
 */

#define gretl_matrix_rows(m) ((m == NULL)? 0 : m->rows)

/**
 * gretl_vector_get_length:
 * @v: vector to examine.
 *
 * Returns: the length of vector @v (without regard to whether
 * it is a row or column vector).
 */

#define gretl_vector_get_length(v) ((v == NULL)? 0 : \
				    ((v)->cols == 1)? (v)->rows :	\
				    ((v)->rows == 1)? (v)->cols : 0)

/**
 * gretl_vector_alloc:
 * @i: number of columns.
 *
 * Returns: a new #gretl_vector with @i columns.
 */

#define gretl_vector_alloc(i) gretl_matrix_alloc(1,(i))

/**
 * gretl_column_vector_alloc:
 * @i: number of rows.
 *
 * Returns: a new column gretl_vector with @i rows.
 */

#define gretl_column_vector_alloc(i) gretl_matrix_alloc((i),1)

/**
 * gretl_vector_free:
 * @v: %gretl_vector to free.
 *
 * Frees the vector @v and its associated storage.
 */

#define gretl_vector_free(v) gretl_matrix_free(v)

/**
 * gretl_matrix_is_scalar:
 * @m: matrix to test.
 *
 * Returns: 1 if @m is 1 x 1, else 0.
 */

#define gretl_matrix_is_scalar(m) ((m) != NULL && \
				   (m)->is_complex == 0 && \
                                   (m)->rows == 1 && \
                                   (m)->cols == 1)

#define gretl_is_null_matrix(m) (m == NULL || m->rows == 0 || m->cols == 0)

#define gretl_is_complex(m) (m != NULL && m->is_complex == 1)

/**
 * sample_size:
 * @p: pointer to dataset.
 *
 * Retrieves the length of the current sample range.
 */
#define sample_size(p) ((p == NULL)? 0 : (p->t2 - p->t1 + 1))
# define UTF_WIDTH(s, w) w
#define NAMETRUNC 18
#define floateq(x, y)  (fabs((x) - (y)) < DBL_EPSILON)
#define floatneq(x, y) (fabs((x) - (y)) > DBL_EPSILON)
#define floatgt(x, y)  ((x) - (y) > DBL_EPSILON)
#define floatlt(x, y)  ((y) - (x) > DBL_EPSILON)

#define ok_int(x) (x <= (double) INT_MAX && x >= (double) INT_MIN)
#define dataset_is_panel(p) (p != NULL && p->structure == STACKED_TIME_SERIES)
#define ERRLEN 2048
/**
 * dataset_is_cross_section:
 * @p: pointer to dataset.
 *
 * Attempt to determine whether a dataset contains cross-sectional
 * data (1) or not (0).
 */
#define dataset_is_cross_section(p) (p != NULL && p->structure == CROSS_SECTION)
/**
 * dataset_is_time_series:
 * @p: pointer to dataset.
 *
 * Attempt to determine whether a dataset contains time series
 * data (1) or not (0).
 */
#define dataset_is_time_series(p) (p != NULL && (p->structure == TIME_SERIES || \
						 p->structure == SPECIAL_TIME_SERIES))

/**
 * dataset_is_seasonal:
 * @p: pointer to dataset.
 *
 * Attempt to determine whether a dataset contains seasonal time series
 * data (1) or not (0).
 */
#define dataset_is_seasonal(p) (p != NULL && (p->structure == TIME_SERIES || \
                                p->structure == SPECIAL_TIME_SERIES) && \
                                p->pd > 1)

/**
 * custom_time_series:
 * @p: pointer to dataset.
 *
 * Attempt to determine whether a dataset contains time series
 * data with custom (non-standard) frequency (1) or not (0).
 */
#define custom_time_series(p) (p != NULL && p->structure == SPECIAL_TIME_SERIES)

/**
 * dataset_is_daily:
 * @p: pointer to dataset.
 *
 * Attempt to determine whether a dataset contains daily time series
 * data (1) or not (0).
 */
#define dataset_is_daily(p) (p != NULL && p->structure == TIME_SERIES \
                             && (p->pd == 5 || p->pd == 6 || p->pd == 7))

/**
 * dataset_is_incomplete_daily:
 * @p: pointer to dataset.
 *
 * Attempt to determine whether a dataset contains daily on an
 * incomplete calendar (1) or not (0).
 */
#define dataset_is_incomplete_daily(p) (p != NULL && p->structure == TIME_SERIES \
					&& (p->pd == 5 || p->pd == 6 || p->pd == 7) \
					&& p->markers == DAILY_DATE_STRINGS)

/**
 * dataset_is_weekly:
 * @p: pointer to dataset.
 *
 * Attempt to determine whether a dataset contains weekly time series
 * data (1) or not (0).
 */
#define dataset_is_weekly(p) (p != NULL && p->structure == TIME_SERIES \
                              && p->pd == 52)

/**
 * dataset_is_hourly:
 * @p: pointer to dataset.
 *
 * Attempt to determine whether a dataset contains hourly time series
 * data (1) or not (0).
 */
#define dataset_is_hourly(p) (p != NULL && p->structure == TIME_SERIES \
                              && p->pd == 24)

/**
 * dataset_is_decennial:
 * @p: pointer to dataset.
 *
 * Attempt to determine whether a dataset contains decennial time series
 * data (1) or not (0).
 */
#define dataset_is_decennial(p) (p != NULL && p->structure == TIME_SERIES \
                                 && p->pd == 10)

/**
 * dated_daily_data:
 * @p: pointer to dataset.
 *
 * Attempt to determine whether a dataset contains dated daily time series
 * data (1) or not (0).
 */
#define dated_daily_data(p) (p != NULL && p->structure == TIME_SERIES \
                             && (p->pd == 5 || p->pd == 6 || p->pd == 7) \
                             && p->sd0 > 100000)

/**
 * undated_daily_data:
 * @p: pointer to dataset.
 *
 * Attempt to determine whether a dataset is daily but does not contain
 * any date information data (1) or not (0).
 */
#define undated_daily_data(p) (p != NULL && p->structure == TIME_SERIES \
                               && (p->pd == 5 || p->pd == 6 || p->pd == 7) \
                               && p->sd0 == 1)

/**
 * dated_weekly_data:
 * @p: pointer to dataset.
 *
 * Attempt to determine whether a dataset contains dated weekly 
 * time series data (1) or not (0).
 */
#define dated_weekly_data(p) (p != NULL && p->structure == TIME_SERIES \
                              && p->pd == 52 && \
                              p->sd0 > 100000)

/**
 * calendar_data:
 * @p: pointer to dataset.
 *
 * Attempt to determine whether a dataset uses calendar
 * dates for observation strings (1) or not (0).
 */
#define calendar_data(p) (p != NULL && p->structure == TIME_SERIES && \
                          (p->pd == 5 || p->pd == 6 || p->pd == 7 || p->pd == 52) && \
			  (p->sd0 > 100000 || strchr(p->stobs, '-')))

/**
 * quarterly_or_monthly:
 * @p: pointer to dataset.
 *
 * Attempt to determine whether a dataset is a quarterly
 * or monthly time series (1), or something else (0).
 */
#define quarterly_or_monthly(p) (p != NULL && p->structure == TIME_SERIES && \
                                 (p->pd == 4 || p->pd == 12))

/**
 * annual_data:
 * @p: pointer to dataset.
 *
 * Attempt to determine whether a dataset is an annual
 * time series (1), or something else (0).
 */
#define annual_data(p) (p != NULL && p->structure == TIME_SERIES && \
			p->pd == 1)

/**
 * decennial_data:
 * @p: pointer to dataset.
 *
 * Attempt to determine whether a dataset is a decemmial
 * time series (1), or something else (0).
 */
#define decennial_data(p) (p != NULL && p->structure == TIME_SERIES && \
			   p->pd == 10 && p->sd0 > 1000)

/**
 * dataset_is_panel:
 * @p: pointer to dataset.
 *
 * Attempt to determine whether a dataset contains panel
 * data (1) or not (0).
 */
#define dataset_is_panel(p) (p != NULL && p->structure == STACKED_TIME_SERIES)

/**
 * dataset_is_seasonal_panel:
 * @p: pointer to dataset.
 *
 * Attempt to determine whether a dataset contains panel
 * data with a seasonal time-series dimension (1) or not (0).
 */
#define dataset_is_seasonal_panel(p) (p != NULL && \
				      p->structure == STACKED_TIME_SERIES && \
				      p->panel_pd > 1)

/**
 * dataset_has_markers:
 * @p: pointer to dataset.
 *
 * Determine whether a dataset has observation marker strings (1)
 * or not (0).
 */
#define dataset_has_markers(p) (p != NULL && p->markers && p->S != NULL)

/**
 * dataset_has_panel_time:
 * @p: pointer to dataset.
 *
 * Determine whether a panel dataset has information on its time
 * dimension recorded (1) or not (0).
 */
#define dataset_has_panel_time(p) (p != NULL && \
				   p->structure == STACKED_TIME_SERIES && \
				   p->panel_pd > 0 && p->panel_sd0 > 0.0)

/**
 * sample_size:
 * @p: pointer to dataset.
 *
 * Retrieves the length of the current sample range.
 */
#define sample_size(p) ((p == NULL)? 0 : (p->t2 - p->t1 + 1))

/**
 * dset_get_data:
 * @d: pointer to dataset.
 * @i: index number of variable.
 * @t: observation number.
 *
 * Gets the value of series @i at observation @t.
 */
#define dset_get_data(d,i,t) (d->Z[i][t])

/**
 * dset_set_data:
 * @d: pointer to dataset.
 * @i: index number of variable.
 * @t: observation number.
 * @x: value to set.
 *
 * Sets the value of series @i at observation @t.
 */
#define dset_set_data(d,i,t,x) (d->Z[i][t]=x)
/* for values that really want a non-negative integer */
#define UNSET_INT -9
#define is_unset(i) (i == UNSET_INT)
#define NEEDS_TS    "needs-time-series-data"
#define NEEDS_QM    "needs-qm-data"
#define NEEDS_PANEL "needs-panel-data"
#define NO_DATA_OK  "no-data-ok"
#define FN_NAMELEN 32
#define ok_function_return_type(r) (r == GRETL_TYPE_DOUBLE || \
				    r == GRETL_TYPE_SERIES || \
				    r == GRETL_TYPE_MATRIX || \
				    r == GRETL_TYPE_LIST ||   \
				    r == GRETL_TYPE_STRING || \
				    r == GRETL_TYPE_BUNDLE || \
				    r == GRETL_TYPE_STRINGS ||  \
				    r == GRETL_TYPE_MATRICES || \
				    r == GRETL_TYPE_BUNDLES ||	\
				    r == GRETL_TYPE_LISTS ||	\
				    r == GRETL_TYPE_ARRAYS ||	\
				    r == GRETL_TYPE_VOID || \
				    r == GRETL_TYPE_NUMERIC)
#define RCOND_MIN  1.0e-14
#define RCOND_WARN 1.0e-07
#define gls_rho(p) gretl_model_get_double_default(p, "rho_gls", 0.0)
#define RS_RCOND_MIN 1.0e-15
#define gretl_is_vector(v) (v->rows == 1 || v->cols == 1)
#define matrix_is_scalar(m) (m->rows == 1 && m->cols == 1)
#define mdx(a,i,j) ((j)*a->rows+(i))
#define matrix_transp_get(m,i,j) (m->val[(i)*m->rows+(j)])
#define matrix_transp_set(m,i,j,x) (m->val[(i)*m->rows+(j)]=x)
#define cmatrix_transp_get(m,i,j) (m->z[(i)*m->rows+(j)])
#define cmatrix_transp_set(m,i,j,x) (m->z[(i)*m->rows+(j)]=x)
#define INFO_INVALID 0xdeadbeef
#define is_block_matrix(m) (m->info == (matrix_info *) INFO_INVALID)
#define is_one_by_one(m) (m->rows == 1 && m->cols == 1)
#define no_metadata(m) (m->info == NULL || is_block_matrix(m))
#define R_DIAG_MIN 1.0e-8
#define CI_TIME (-2)
#define CI_UNIT (-3)
#define by_time_and_unit(ci) (ci->dcvar0 < -1 && ci->dcvar1 < -1)
#define by_time_or_unit(ci) (ci->dcvar0 < -1 || ci->dcvar1 < -1)
#define by_time(ci) (ci->dcvar0 == CI_TIME)
#define by_unit(ci) (ci->dcvar0 == CI_UNIT)
#define generic(ci) (ci->dcvar0 > 0)
#define two_way(ci) (ci->dcvar0 != -1 && ci->dcvar1 != -1)
#define is_present(v) (v != -1)
/* commands for which --vcv (OPT_O) is applicable */
#define vcv_opt_ok(c) (MODEL_COMMAND(c) || c == ADD || c == OMIT)

/* commands for which --window (OPT_W) is applicable */
#define window_opt_ok(c) (MODEL_COMMAND(c) || c == VAR || c == VECM)

/* commands which support the $result accessor */
#define yields_result(c) (c == CORR || c == FREQ || c == SUMMARY)

/* commands for which --quiet (OPT_Q) is applicable */
#define quiet_opt_ok(c) (MODEL_COMMAND(c) ||    \
                         yields_result(c) ||    \
                         c == ADD ||            \
                         c == ADF ||            \
                         c == ANOVA ||          \
                         c == APPEND ||         \
			 c == BDS ||		\
                         c == BKW ||            \
                         c == COEFFSUM ||       \
                         c == CHOW ||           \
                         c == COINT2 ||         \
                         c == CORRGM ||         \
                         c == CUSUM ||          \
                         c == DATA ||           \
                         c == DIFFTEST ||       \
                         c == ESTIMATE ||       \
                         c == FCAST ||          \
                         c == FOREIGN ||        \
                         c == FRACTINT ||       \
                         c == FREQ ||           \
                         c == KPSS ||           \
                         c == MAKEPKG ||        \
                         c == MODTEST ||        \
                         c == LEVERAGE ||       \
                         c == LEVINLIN ||       \
                         c == LOOP ||           \
                         c == MAHAL ||          \
                         c == NORMTEST ||       \
                         c == OLS ||            \
                         c == OMIT ||           \
                         c == OPEN ||           \
			 c == PANSPEC ||	\
                         c == PKG ||            \
                         c == QLRTEST ||        \
                         c == RENAME ||         \
                         c == RESET ||          \
                         c == RESTRICT ||       \
                         c == RMPLOT ||         \
                         c == SMPL ||           \
                         c == SYSTEM ||         \
                         c == VAR ||            \
                         c == VECM ||           \
                         c == VIF ||            \
                         c == XCORRGM ||        \
                         c == XTAB)

/* --output (OPT_U) as attached to GNUPLOT */
#define plot_output_opt_ok(c) (c == GNUPLOT ||	\
			       c == PLOT ||	\
			       c == BXPLOT ||	\
			       c == HFPLOT ||	\
			       c == PANPLOT ||	\
			       c == QQPLOT ||	\
			       c == KDPLOT ||   \
			       c == RMPLOT ||	\
			       c == SCATTERS || \
			       c == GRIDPLOT)

/* --plot (OPT_U) as attached to CORR */
#define cmd_plot_opt_ok(c) (c == CORR ||	\
			    c == CORRGM ||	\
			    c == CUSUM ||	\
			    c == FCAST ||	\
			    c == FREQ ||	\
			    c == HURST ||	\
			    c == LEVERAGE ||	\
			    c == PERGM ||	\
			    c == QLRTEST ||	\
			    c == XCORRGM)

/* --outbuf (OPT_b) as attached to GNUPLOT */
#define plot_outbuf_opt_ok(c) (plot_output_opt_ok(c) || \
			       cmd_plot_opt_ok(c))

#define MODEL_COMMAND(c) (c == AR || \
                          c == AR1 || \
                          c == ARCH || \
                          c == ARMA || \
			  c == DPANEL ||   \
			  c == DURATION || \
                          c == GARCH || \
                          c == GMM || \
		          c == HECKIT || \
                          c == HSK || \
                          c == INTREG || \
                          c == IVREG || \
                          c == LAD || \
                          c == LOGISTIC || \
                          c == LOGIT || \
			  c == MIDASREG || \
                          c == MLE || \
                          c == MPOLS || \
			  c == NEGBIN || \
                          c == NLS || \
                          c == OLS || \
                          c == PANEL || \
                          c == POISSON || \
                          c == PROBIT || \
                          c == BIPROBIT || \
                          c == QUANTREG || \
                          c == TOBIT || \
                          c == WLS)

#define AR_MODEL(c) (c == AR || \
		     c == AR1 || \
                     c == ARMA || \
                     c == GARCH)

#define SIMPLE_AR_MODEL(c) (c == AR || c == AR1)

#define ML_ESTIMATOR(c) (c == ARMA || \
			 c == DURATION || \
                         c == GARCH || \
                         c == HECKIT || \
                         c == LOGIT || \
                         c == MLE || \
			 c == NEGBIN ||	\
                         c == POISSON || \
                         c == PROBIT || \
			 c == BIPROBIT || \
                         c == TOBIT)

#define LIMDEP(c) (c == LOGIT || \
                   c == PROBIT || \
                   c == TOBIT || \
                   c == INTREG)

#define COUNT_MODEL(c) (c == POISSON || c == NEGBIN)

#define LSQ_MODEL(c) (c == AR1 || \
                      c == HSK || \
                      c == OLS || \
                      c == WLS)

#define ASYMPTOTIC_MODEL(c) (c == ARMA || \
			     c == DPANEL ||   \
			     c == DURATION || \
                             c == GARCH || \
                             c == GMM || \
                             c == HECKIT || \
                             c == INTREG || \
                             c == IVREG || \
                             c == LOGIT || \
                             c == MLE || \
			     c == NEGBIN || \
                             c == POISSON || \
                             c == PROBIT || \
                             c == TOBIT || \
                             c == BIPROBIT)

#define EQN_SYSTEM_COMMAND(c) (c == VAR || c == VECM || c == SYSTEM)

/* model where the specification is not based on a list
   of variables */
#define NONLIST_MODEL(c) (c == NLS || c == MLE || c == GMM || c == MIDASREG)

#define is_model_ref_cmd(c) (c == ADD || \
	                     c == ARCH || \
	                     c == CHOW || \
	                     c == CUSUM || \
	                     c == FCAST || \
                             c == LEVERAGE || \
	                     c == MODTEST || \
                             c == OMIT || \
	                     c == RESTRICT || \
                             c == VIF)

#define RQ_SPECIAL_MODEL(m) ((m->ci == LAD || m->ci == QUANTREG) &&	\
                             NULL != gretl_model_get_data(m, "rq_tauvec"))

#define POOLED_MODEL(m) ((m->ci == OLS || m->ci == PANEL) && \
                         gretl_model_get_int(m, "pooled"))
#define masked(m,t) (m != NULL && m[t] == '1')
#define gridlimit 2048
#define panel_index(i,t) (i * panidx.T + t + panidx.offset)
#define panel_missing(p, t) (na(p->pooled->uhat[t]))
#define INIT_SIZE 2048
#define MINREM 1024
#define FEW_VALS 32
#define FEWER_VALS 8
#define CHOL_TINY  8.0e-09
#define CHOL_SMALL 1.0e-08
#define CHOL_RCOND_MIN 1.0e-6

#ifdef unix
#define min MIN
#define max MAX
#endif

#define NREPEAT 100
#define DEFAULT_EQTOL 1.0e-9 /* 2014-08-05: was 1.5e-12 */
#define TOEPLITZ_SMALL 1.0e-20
#define QFORM_SMALL 1.0e-20
#define RS_RCOND_MIN 1.0e-15
#define gretl_matrix_cum(m,i,j,x) (m->val[(j)*m->rows+(i)]+=x)
#define gretl_st_result(c,i,j,x,m)                      \
    do {                                                \
        if (m==GRETL_MOD_CUMULATE) {                    \
            c->val[(j)*c->rows+(i)]+=x;                 \
            if (i!=j) c->val[(i)*c->rows+(j)]+=x;       \
        } else if (m==GRETL_MOD_DECREMENT) {            \
            c->val[(j)*c->rows+(i)]-=x;                 \
            if (i!=j) c->val[(i)*c->rows+(j)]-=x;       \
        } else {                                        \
            gretl_matrix_set(c,i,j,x);                  \
            gretl_matrix_set(c,j,i,x);                  \
        }                                               \
    } while (0);
#define true_null_matrix(a) (a->rows == 0 && a->cols == 0)
#define complete_obs(x,y,i) (!na(x[i]) && !na(y[i]))
#define has_colnames(m) (m != NULL && !is_block_matrix(m) && \
                         m->info != NULL && m->info->colnames != NULL)

#define has_rownames(m) (m != NULL && !is_block_matrix(m) && \
                         m->info != NULL && m->info->rownames != NULL)


typedef struct global_vars_ global_vars;
typedef struct parser_ GENERATOR;
typedef struct node NODE;

typedef char   gchar;
typedef short  gshort;
typedef long   glong;
typedef int    gint;
typedef gint   gboolean;

typedef unsigned char   guchar;
typedef unsigned short  gushort;
typedef unsigned long   gulong;
typedef unsigned int    guint;

typedef float   gfloat;
typedef double  gdouble;

typedef int32_t integer;
typedef int logical;
typedef float real;
typedef double doublereal;
typedef int flag;
typedef int ftnlen;
typedef logical (*L_fp)();
typedef struct _cmplx {
    double r;
    double i;
} cmplx;

typedef signed char gint8;
typedef signed short gint16;
typedef signed int gint32;
typedef signed long gint64;
typedef unsigned char guint8;
typedef unsigned short guint16;
typedef unsigned int guint32;
typedef unsigned long guint64;
typedef void* gpointer;
typedef struct GPtrArray_ {
  gpointer* pdata;
  guint len;
} GPtrArray;
typedef struct DATASET_ DATASET;
typedef struct matrix_info_ matrix_info;
/* The following apparatus is used for

   (a) setting and retrieving parameters associated with
       command options, as in --opt=val, and

   (b) storing options for a specified command via the
       "setopt" command (with or without parameters).
*/
typedef struct stored_opt_ stored_opt;
