#include <common.h>
__import globals;
__import gretltypes;
__import common;
__import dataset;

/**
 * get_optval_int:
 * @ci: gretl command index.
 * @opt: gretl option value.
 * @err: location to receive error code.
 *
 * Returns: the integer ancillary value currently
 * associated with option @opt for command @ci, if any,
 * otherwise 0. A non-zero value is written to @err if
 * such a value is required for the option in question
 * but is not present.
 */

//~ int get_optval_int (int ci, gretlopt opt, int *err)
//~ {
    //~ stored_opt *so = matching_stored_opt(ci, opt);
    //~ int status = option_parm_status(ci, opt);
    //~ int ret = 0;

    //~ if (so != NULL && so->val != NULL) {
        //~ ret = gretl_int_from_string(so->val, err);
        //~ if (*err) {
            //~ ret = generate_int(so->val, NULL, err);
        //~ }
        //~ if (*err) {
            //~ gretl_errmsg_sprintf(_("%s: invalid option argument"), so->val);
            //~ *err = E_INVARG;
        //~ }
    //~ } else if (status == 2 && err != NULL) {
        //~ const char *longopt = get_longopt(ci, opt);

        //~ gretl_errmsg_sprintf(_("The option '--%s' requires a parameter"),
                             //~ longopt);
        //~ *err = E_DATA;
    //~ }

    //~ return ret;
//~ }

/* Below: This is used as a one-way mapping from the long form to the
   index (e.g. OPT_Q), so a given index can have more than one
   long-form counterpart depending on the context.  The last field
   indicates whether the given option does not accept (0), accepts
   (1), or requires (2) an associated parameter value.
*/

struct gretl_option gretl_opts[] = {
    { ADD,      OPT_B, "both", 0 },
    { ADD,      OPT_I, "silent", 0 },
    { ADD,      OPT_L, "lm", 0 },
    { ADF,      OPT_N, "nc", 0 },
    { ADF,      OPT_C, "c", 0 },
    { ADF,      OPT_D, "seasonals", 0 },
    { ADF,      OPT_R, "ctt", 0 },
    { ADF,      OPT_T, "ct", 0 },
    { ADF,      OPT_V, "verbose", 1 },
    { ADF,      OPT_F, "difference", 0 },
    { ADF,      OPT_E, "test-down", 1 },
    { ADF,      OPT_G, "gls", 0 },
    { ADF,      OPT_U, "perron-qu", 0 },
    { AR1,      OPT_B, "no-corc", 0 },
    { AR1,      OPT_H, "hilu", 0 },
    { AR1,      OPT_P, "pwe", 0 },
    { AR1,      OPT_L, "loose", 0 },
    { APPEND,   OPT_A, "all-cols", 0 },
    { APPEND,   OPT_T, "time-series", 0 },
    { APPEND,   OPT_R, "rowoffset", 2 },
    { APPEND,   OPT_C, "coloffset", 2 },
    { APPEND,   OPT_F, "fixed-cols", 2 },
    { APPEND,   OPT_M, "rowmask", 2 },
    { APPEND,   OPT_L, "cols", 2 },
    { APPEND,   OPT_S, "sheet", 2 },
    { APPEND,   OPT_V, "verbose", 0 },
    { APPEND,   OPT_U, "update-overlap", 0 },
    { APPEND,   OPT_X, "fixed-sample", 0 },
    { APPEND,   OPT_K, "frompkg", 2 },
    { ARMA,     OPT_A, "as154", 0 },
    { ARMA,     OPT_C, "conditional", 0 },
    { ARMA,     OPT_E, "save-ehat", 0 },
    { ARMA,     OPT_G, "opg", 0 },
    { ARMA,     OPT_H, "hessian", 0 },
    { ARMA,     OPT_K, "kalman", 0 },
    { ARMA,     OPT_L, "lbfgs", 0 },
    { ARMA,     OPT_N, "nc", 0 },
    { ARMA,     OPT_V, "verbose", 0 },
    { ARMA,     OPT_X, "x-13arima", 0 },
    { ARMA,     OPT_X, "x-12-arima", 0 }, /* compatibility alias */
    { ARMA,     OPT_Y, "y-diff-only", 0 },
    { ARMA,     OPT_R, "robust", 0 },
    { ARMA,     OPT_S, "stdx", 0 },
    { ARMA,     OPT_Z, "lagselect", 0 },
    { BDS,      OPT_B, "boot", 2 },
    { BDS,      OPT_C, "corr1", 2 },
    { BDS,      OPT_S, "sdcrit", 2 },
    { BDS,      OPT_X, "matrix", 2 },
    { BIPROBIT, OPT_G, "opg", 0 },
    { BIPROBIT, OPT_R, "robust", 0 },
    { BIPROBIT, OPT_V, "verbose", 0 },
    { BIPROBIT, OPT_C, "cluster", 2 },
    { BIPROBIT, OPT_X, "save-xbeta", 0 },
    { BXPLOT,   OPT_O, "notches", 0 },
    { BXPLOT,   OPT_K, "tweaks", 2 },
    { BXPLOT,   OPT_L, "outliers", 2 },
    { BXPLOT,   OPT_P, "panel", 0 },
    { BXPLOT,   OPT_X, "matrix", 2 },
    { BXPLOT,   OPT_Z, "factorized", 0 },
    { BXPLOT,   OPT_B, "whiskerbars", 0 },
    { CHOW,     OPT_D, "dummy", 0 },
    { CHOW,     OPT_L, "limit-to", 2 },
    { CLEAR,    OPT_D, "dataset", 0 },
    { CLEAR,    OPT_F, "functions", 0 },
    { COINT,    OPT_D, "seasonals", 0 },
    { COINT,    OPT_E, "test-down", 1 },
    { COINT,    OPT_N, "nc", 0 },
    { COINT,    OPT_R, "ctt", 0 },
    { COINT,    OPT_S, "skip-df", 0 },
    { COINT,    OPT_T, "ct", 0 },
    { COINT,    OPT_V, "verbose", 0 },
    { COINT,    OPT_I, "silent", 0 },
    { COINT2,   OPT_A, "crt", 0 },
    { COINT2,   OPT_D, "seasonals", 0 },
    { COINT2,   OPT_N, "nc", 0 },
    { COINT2,   OPT_R, "rc", 0 },
    { COINT2,   OPT_C, "uc", 0 },
    { COINT2,   OPT_T, "ct", 0 },
    { COINT2,   OPT_S, "silent", 0 },
    { COINT2,   OPT_V, "verbose", 0 },
    { COINT2,   OPT_Y, "asy", 0 },
    { CORR,     OPT_K, "kendall", 0 },
    { CORR,     OPT_S, "spearman", 0 },
    { CORR,     OPT_N, "uniform", 0 },
    { CORR,     OPT_V, "verbose", 0 },
    { CORR,     OPT_U, "plot", 2 },
    { CORR,     OPT_X, "matrix", 2 },
    { CORR,     OPT_T, "triangle", 0 },
    { CORRGM,   OPT_B, "bartlett", 0 },
    { CUSUM,    OPT_R, "squares", 0 },
    { DATA,     OPT_C, "compact", 2 },
    { DATA,     OPT_O, "odbc", 0 },
    { DATA,     OPT_N, "name", 2 },
    { DATA,     OPT_V, "verbose", 0 },
    { DATA,     OPT_F, "no-align", 0 },
    { DATAMOD,  OPT_P, "preserve", 0 },
    { DATAMOD,  OPT_T, "panel-time", 0 },
    { DATAMOD,  OPT_W, "weekstart", 2 },
    { DATAMOD,  OPT_R, "repday", 2 },
    { DELEET,   OPT_D, "db", 0 },
    { DELEET,   OPT_F, "force", 0 },
    { DELEET,   OPT_L, "list", 0 },
    { DELEET,   OPT_T, "type", 2 },
    { DIFFTEST, OPT_G, "sign", 0 },
    { DIFFTEST, OPT_R, "rank-sum", 0 },
    { DIFFTEST, OPT_I, "signed-rank", 0 },
    { DIFFTEST, OPT_V, "verbose", 0 },
    { DISCRETE, OPT_R, "reverse", 0 },
    { DPANEL,   OPT_A, "asymptotic", 0 },
    { DPANEL,   OPT_D, "time-dummies", 1 },
    { DPANEL,   OPT_K, "keep-extra", 0 },
    { DPANEL,   OPT_L, "system", 0 },
    { DPANEL,   OPT_T, "two-step", 0 },
    { DPANEL,   OPT_V, "verbose", 0 },
    { DPANEL,   OPT_X, "dpdstyle", 0 },
    { DPANEL,   OPT_C, "collapse", 0 },
    { DUMMIFY,  OPT_F, "drop-first", 0 },
    { DUMMIFY,  OPT_L, "drop-last", 0 },
    { DURATION, OPT_B, "weibull", 0 },
    { DURATION, OPT_E, "exponential", 0 },
    { DURATION, OPT_L, "loglogistic", 0 },
    { DURATION, OPT_Z, "lognormal", 0 },
    { DURATION, OPT_M, "medians", 0 },
    { DURATION, OPT_G, "opg", 0 },
    { DURATION, OPT_R, "robust", 0 },
    { DURATION, OPT_C, "cluster", 2 },
    { DURATION, OPT_V, "verbose", 0 },
    { EQNPRINT, OPT_O, "complete", 0 },
    { EQNPRINT, OPT_T, "t-ratios", 0 },
    { EQNPRINT, OPT_U, "output", 2 },
    { TABPRINT, OPT_O, "complete", 0 },
    { TABPRINT, OPT_C, "csv", 0 },
    { TABPRINT, OPT_R, "rtf", 0 },
    { TABPRINT, OPT_T, "format", 2 },
    { TABPRINT, OPT_U, "output", 2 },
    { EQUATION, OPT_M, "multi", 0 },
    { ESTIMATE, OPT_I, "iterate", 0 },
    { ESTIMATE, OPT_M, "geomean", 0 },
    { ESTIMATE, OPT_N, "no-df-corr", 0 },
    { ESTIMATE, OPT_U, "unrestrict-init", 0 },
    { ESTIMATE, OPT_V, "verbose", 0 },
    { ESTIMATE, OPT_W, "window", 0 },
    { FCAST,    OPT_L, "all-probs", 0 },
    { FCAST,    OPT_D, "dynamic", 0 },
    { FCAST,    OPT_M, "mean-y", 0 },
    { FCAST,    OPT_N, "no-stats", 0 },
    { FCAST,    OPT_T, "stats-only", 0 },
    { FCAST,    OPT_S, "static", 0 },
    { FCAST,    OPT_R, "recursive", 0 },
    { FCAST,    OPT_R, "rolling", 0 }, /* legacy alias */
    { FCAST,    OPT_O, "out-of-sample", 0 },
    { FCAST,    OPT_I, "integrate", 0 },
    { FOREIGN,  OPT_D, "send-data", 1 },
    { FOREIGN,  OPT_V, "verbose", 0 },
    { FOREIGN,  OPT_F, "frame", 0 },
    { FOREIGN,  OPT_N, "no-compile", 0 },
    { FOREIGN,  OPT_I, "io-funcs", 2 },
    { FRACTINT, OPT_G, "gph", 0 },
    { FRACTINT, OPT_A, "all", 0 },
    { FREQ,     OPT_G, "show-plot", 0 }, /* legacy */
    { FREQ,     OPT_O, "gamma", 0 },
    { FREQ,     OPT_S, "silent", 0 },
    { FREQ,     OPT_Z, "normal", 0 },
    { FREQ,     OPT_N, "nbins", 2 },
    { FREQ,     OPT_M, "min", 2 },
    { FREQ,     OPT_W, "binwidth", 2 },
    { FREQ,     OPT_X, "matrix", 2 },
    { FREQ,     OPT_K, "tweaks", 2 },
    { GARCH,    OPT_A, "arma-init", 0 },
    { GARCH,    OPT_F, "fcp", 0 },
    { GARCH,    OPT_N, "nc", 0 },
    { GARCH,    OPT_R, "robust", 0 },
    { GARCH,    OPT_V, "verbose", 0 },
    { GARCH,    OPT_Z, "stdresid", 0 },
    { GMM,      OPT_I, "iterate", 0 },
    { GMM,      OPT_L, "lbfgs", 0 },
    { GMM,      OPT_T, "two-step", 0 },
    { GMM,      OPT_V, "verbose", 0 },
    { GNUPLOT,  OPT_I, "input", 2 },
    { GNUPLOT,  OPT_O, "with-lines", 1 },
    { GNUPLOT,  OPT_F, "fit", 2 },
    { GNUPLOT,  OPT_K, "tweaks", 2 },
    { GNUPLOT,  OPT_M, "with-impulses", 1 },
    { GNUPLOT,  OPT_P, "with-lp", 1 },
    { GNUPLOT,  OPT_B, "with-boxes", 1 },
    { GNUPLOT,  OPT_Q, "with-steps", 1 },
    { GNUPLOT,  OPT_T, "time-series", 0 },
    { GNUPLOT,  OPT_Z, "dummy", 0 },
    { GNUPLOT,  OPT_C, "control", 0 },
    { GNUPLOT,  OPT_U, "output", 2 },
    { GNUPLOT,  OPT_b, "outbuf", 2 },
    { GNUPLOT,  OPT_b, "buffer", 2 }, /* compatibility alias */
    { GNUPLOT,  OPT_i, "inbuf", 2 },
    { GNUPLOT,  OPT_Y, "single-yaxis", 0 },
    { GNUPLOT,  OPT_X, "matrix", 2 },
    { GNUPLOT,  OPT_N, "band", 2 },
    { GNUPLOT,  OPT_J, "band-style", 2 },
    { GNUPLOT,  OPT_a, "bands", 2 },
    { GNUPLOT,  OPT_W, "font", 2 },
    { GNUPLOT,  OPT_L, "ylogscale", 1 },
    { GRAPHPG,  OPT_M, "monochrome", 0 },
    { GRAPHPG,  OPT_O, "output", 2 },
    { GRIDPLOT, OPT_F, "fontsize", 2 },
    { GRIDPLOT, OPT_W, "width", 2 },
    { GRIDPLOT, OPT_H, "height", 2 },
    { GRIDPLOT, OPT_R, "rows", 2 },
    { GRIDPLOT, OPT_C, "cols", 2 },
    { GRIDPLOT, OPT_L, "layout", 2 },
    { GRIDPLOT, OPT_U, "output", 2 },
    { HECKIT,   OPT_M, "ml", 0 },
    { HECKIT,   OPT_G, "opg", 0 },
    { HECKIT,   OPT_R, "robust", 0 },
    { HECKIT,   OPT_C, "cluster", 2 },
    { HECKIT,   OPT_T, "two-step", 0 },
    { HECKIT,   OPT_V, "verbose", 0 },
    { HELP,     OPT_F, "func", 0 },
    { HFPLOT,   OPT_O, "with-lines", 0 },
    { HFPLOT,   OPT_T, "time-series", 0 },
    { HSK,      OPT_N, "no-squares", 0 },
    { INCLUDE,  OPT_F, "force", 0 },
    { INTREG,   OPT_G, "opg", 0 },
    { INTREG,   OPT_R, "robust", 0 },
    { INTREG,   OPT_C, "cluster", 2 },
    { INTREG,   OPT_V, "verbose", 0 },
    { IVREG,    OPT_G, "gmm", 0 },
    { IVREG,    OPT_I, "iterate", 0 },
    { IVREG,    OPT_L, "liml", 0 },
    { IVREG,    OPT_N, "no-df-corr", 0 },
    { IVREG,    OPT_R, "robust", 0 },
    { IVREG,    OPT_S, "save", 0 },
    { IVREG,    OPT_T, "two-step", 0 },
    { IVREG,    OPT_H, "weights", 2 },
    { IVREG,    OPT_X, "no-tests", 0 },
    { IVREG,    OPT_C, "cluster", 2 },
    { JOIN,     OPT_I, "ikey", 2 },
    { JOIN,     OPT_O, "okey", 2 },
    { JOIN,     OPT_F, "filter", 2 },
    { JOIN,     OPT_A, "aggr", 2 },
    { JOIN,     OPT_D, "data", 2 },
    { JOIN,     OPT_T, "tconv-fmt", 2 },
    { JOIN,     OPT_K, "tkey", 2 },
    { JOIN,     OPT_X, "tconvert", 2 },
    { JOIN,     OPT_H, "no-header", 0 },
    { JOIN,     OPT_V, "verbose", 0 },
    { JOIN,     OPT_P, "pd", 2 }, /* undocumented: is it wanted? */
    { JOIN,     OPT_R, "frompkg", 2 },
    { KDPLOT,   OPT_O, "alt", 0 },
    { KDPLOT,   OPT_S, "scale", 2},
    { KPSS,     OPT_T, "trend", 0 },
    { KPSS,     OPT_D, "seasonals", 0 },
    { KPSS,     OPT_V, "verbose", 0 },
    { KPSS,     OPT_F, "difference", 0 },
    { LAGS,     OPT_L, "bylag", 0 },
    { LEVERAGE, OPT_S, "save", 0 },
    { LEVERAGE, OPT_O, "overwrite", 0 },
    { LEVINLIN, OPT_N, "nc", 0 },
    { LEVINLIN, OPT_T, "ct", 0 },
    { LEVINLIN, OPT_V, "verbose", 0 },
    { MAKEPKG,  OPT_D, "dtd", 2 },
    { MAKEPKG,  OPT_I, "index", 0 },
    { MAKEPKG,  OPT_T, "translations", 0 },
    { MARKERS,  OPT_D, "delete", 0 },
    { MARKERS,  OPT_F, "from-file", 2 },
    { MARKERS,  OPT_T, "to-file", 2 },
    { MARKERS,  OPT_A, "to-array", 2 },
    { MARKERS,  OPT_R, "from-array", 2 },
    { MARKERS,  OPT_S, "from-series", 2 },
    { MODTEST,  OPT_A, "autocorr", 0 },
    { MODTEST,  OPT_B, "breusch-pagan", 0 },
    { MODTEST,  OPT_C, "comfac", 0 },
    { MODTEST,  OPT_D, "xdepend", 0 },
    { MODTEST,  OPT_H, "arch", 0 },
    { MODTEST,  OPT_L, "logs", 0 },
    { MODTEST,  OPT_N, "normality", 0 },
    { MODTEST,  OPT_S, "squares", 0 },
    { MODTEST,  OPT_P, "panel", 0 },
    { MODTEST,  OPT_R, "robust", 0 },
    { MODTEST,  OPT_W, "white", 0 },
    { MODTEST,  OPT_X, "white-nocross", 0 },
    { MODTEST,  OPT_I, "silent", 0 },
    { MODTEST,  OPT_U, "univariate", 0 },
    { MPI,      OPT_D, "send-data", 1 },
    { MPI,      OPT_F, "send-functions", 0 },
    { MPI,      OPT_L, "local", 0 },
    { MPI,      OPT_N, "np", 2 },
    { MPI,      OPT_T, "omp-threads", 2 },
    { MPI,      OPT_Q, "quiet", 0},
    { MPI,      OPT_V, "verbose", 1 },
    { MPI,      OPT_S, "single-rng", 0},
    { LABELS,   OPT_D, "delete", 0 },
    { LABELS,   OPT_F, "from-file", 2 },
    { LABELS,   OPT_T, "to-file", 2 },
    { LABELS,   OPT_A, "from-array", 2 },
    { LABELS,   OPT_R, "to-array", 2 },
    { LAD,      OPT_N, "no-vcv", 0 },
    { LOGISTIC, OPT_M, "ymax", 2 },
    { LOGISTIC, OPT_R, "robust", 0 },
    { LOGISTIC, OPT_C, "cluster", 2 },
    { LOGISTIC, OPT_F, "fixed-effects", 0 },
    { LOGIT,    OPT_M, "multinomial", 0 },
    { LOGIT,    OPT_P, "p-values", 0 },
    { LOGIT,    OPT_R, "robust", 0 },
    { LOGIT,    OPT_C, "cluster", 2 },
    { LOGIT,    OPT_V, "verbose", 0 },
    { LOGIT,    OPT_S, "estrella", 0 },
    { LOOP,     OPT_D, "decr", 0 },
    { LOOP,     OPT_P, "progressive", 0 },
    { LOOP,     OPT_V, "verbose", 0 },
    { MAHAL,    OPT_S, "save", 0 },
    { MAHAL,    OPT_V, "vcv", 0 },
    { MEANTEST, OPT_O, "unequal-vars", 0 },
    { MIDASREG, OPT_R, "robust", 0 },
    { MIDASREG, OPT_V, "verbose", 0 },
    { MIDASREG, OPT_L, "levenberg", 0 },
    { MIDASREG, OPT_P, "print-spec", 0 },
    { MIDASREG, OPT_B, "breaktest", 0 },
    { MIDASREG, OPT_C, "clamp-beta", 0 }, /* legacy */
    { MLE,      OPT_A, "auxiliary", 0 },
    { MLE,      OPT_G, "opg", 0 },
    { MLE,      OPT_H, "hessian", 0 },
    { MLE,      OPT_S, "no-gradient-check", 0 },
    { MLE,      OPT_L, "lbfgs", 0 },
    { MLE,      OPT_N, "numerical", 0 },
    { MLE,      OPT_R, "robust", 1 },
    { MLE,      OPT_C, "cluster", 2 },
    { MLE,      OPT_V, "verbose", 0 },
    { MODPRINT, OPT_A, "addstats", 2 },
    { MODPRINT, OPT_O, "output", 2 },
    { MODPRINT, OPT_C, "complete", 0 },
    { MODELTAB, OPT_O, "output", 2 },
    { MODELTAB, OPT_C, "complete", 0 },
    { MPOLS,    OPT_S, "simple-print", 0 },
    { NEGBIN,   OPT_G, "opg", 0 },
    { NEGBIN,   OPT_M, "model1", 0 },
    { NEGBIN,   OPT_R, "robust", 0 },
    { NEGBIN,   OPT_C, "cluster", 2 },
    { NEGBIN,   OPT_V, "verbose", 0 },
    { NLS,      OPT_N, "numerical", 0 },
    { NLS,      OPT_R, "robust", 0 },
    { NLS,      OPT_V, "verbose", 0 },
    { NLS,      OPT_S, "no-gradient-check", 0 },
    { NORMTEST, OPT_A, "all", 0 },
    { NORMTEST, OPT_D, "dhansen", 0 },
    { NORMTEST, OPT_W, "swilk", 0 },
    { NORMTEST, OPT_J, "jbera", 0 },
    { NORMTEST, OPT_L, "lillie", 0 },
    { NULLDATA, OPT_N, "no-index", 0 },
    { NULLDATA, OPT_P, "preserve", 0 },
    { OLS,      OPT_F, "print-final", 0 },
    { OLS,      OPT_J, "jackknife", 0 },
    { OLS,      OPT_N, "no-df-corr", 0 },
    { OLS,      OPT_O, "vcv", 0 },
    { OLS,      OPT_R, "robust", 0 },
    { OLS,      OPT_Q, "quiet", 0 }, /* note: for the sake of documentation */
    { OLS,      OPT_S, "simple-print", 0 },
    { OLS,      OPT_V, "anova", 0 },
    { OLS,      OPT_C, "cluster", 2 },
    { OLS,      OPT_W, "window", 0 },
    { OMIT,     OPT_A, "auto", 1 },
    { OMIT,     OPT_B, "both", 0 },
    { OMIT,     OPT_X, "chi-square", 0 },
    { OMIT,     OPT_I, "silent", 0 },
    { OMIT,     OPT_W, "test-only", 0 },
    { OMIT,     OPT_T, "trend", 0 },      /* omit auto-trend: VAR only */
    { OMIT,     OPT_E, "seasonals", 0 },  /* omit auto-seasonals: VAR only */
    { OPEN,     OPT_A, "all-cols", 0 },
    { OPEN,     OPT_B, "progress-bar", 0 },
    { OPEN,     OPT_D, "drop-empty", 0 },
    { OPEN,     OPT_E, "select", 2 },
    { OPEN,     OPT_F, "fixed-cols", 2 },
    { OPEN,     OPT_O, "odbc", 0 },
    { OPEN,     OPT_P, "preserve", 0 },
    { OPEN,     OPT_R, "rowoffset", 2 },
    { OPEN,     OPT_C, "coloffset", 2 },
    { OPEN,     OPT_S, "sheet", 2 },
    { OPEN,     OPT_W, "www", 0 },
    { OPEN,     OPT_L, "cols", 2 },
    { OPEN,     OPT_M, "rowmask", 2 },
    { OPEN,     OPT_V, "verbose", 0 },
    { OPEN,     OPT_K, "frompkg", 2 },
    { OPEN,     OPT_H, "no-header", 0 },
    { OPEN,     OPT_I, "ignore-quotes", 0 },
    { OPEN,     OPT_U, "bundle", 2 },
    { OUTFILE,  OPT_A, "append", 0 },
    { OUTFILE,  OPT_Q, "quiet", 0 },
    { OUTFILE,  OPT_B, "buffer", 1 }, /* note: 1 is for backward compat */
    { OUTFILE,  OPT_T, "tempfile", 2 },
    { PANEL,    OPT_B, "between", 0 },
    { PANEL,    OPT_C, "cluster", 2 },
    { PANEL,    OPT_D, "time-dummies", 1 },
    { PANEL,    OPT_E, "nerlove", 0 },
    { PANEL,    OPT_F, "fixed-effects", 0 },
    { PANEL,    OPT_I, "iterate", 0 },
    { PANEL,    OPT_M, "matrix-diff", 0 },
    { PANEL,    OPT_N, "no-df-corr", 0 },
    { PANEL,    OPT_P, "pooled", 0 },
    { PANEL,    OPT_R, "robust", 0 },
    { PANEL,    OPT_U, "random-effects", 0 },
    { PANEL,    OPT_V, "verbose", 0 },
    { PANEL,    OPT_H, "unit-weights", 0 },
    { PANEL,    OPT_X, "unbalanced", 1 },
    { PANPLOT,  OPT_M, "means", 0 },
    { PANPLOT,  OPT_V, "overlay", 0 },
    { PANPLOT,  OPT_S, "sequence", 0 },
    { PANPLOT,  OPT_D, "grid", 0 },
    { PANPLOT,  OPT_A, "stack", 0 },
    { PANPLOT,  OPT_B, "boxplots", 0 },
    { PANPLOT,  OPT_C, "boxplot", 0 },
    { PANPLOT,  OPT_Y, "single-yaxis", 0 },
    { PANSPEC,  OPT_M, "matrix-diff", 0 },
    { PANSPEC,  OPT_N, "nerlove", 0 },
    { POISSON,  OPT_R, "robust", 0 },
    { POISSON,  OPT_C, "cluster", 2 },
    { POISSON,  OPT_V, "verbose", 0 },
    { PCA,      OPT_C, "covariance", 0 },
    { PCA,      OPT_A, "save-all", 0 },
    { PCA,      OPT_O, "save", 1 },
    { PCA,      OPT_Q, "quiet", 0 },
    { PERGM,    OPT_O, "bartlett", 0 },
    { PERGM,    OPT_L, "log", 0 },
    { PERGM,    OPT_R, "radians", 0 },
    { PERGM,    OPT_D, "degrees", 0 },
    { PKG,      OPT_L, "local", 0 },
    { PKG,      OPT_V, "verbose", 0 },
    { PLOT,     OPT_C, "control", 0 },
    { PLOT,     OPT_O, "with-lines", 1 },
    { PLOT,     OPT_F, "fit", 2 },
    { PLOT,     OPT_B, "with-boxes", 1 },
    { PLOT,     OPT_Q, "with-steps", 1 },
    { PLOT,     OPT_M, "with-impulses", 1 },
    { PLOT,     OPT_P, "with-lp", 1 },
    { PLOT,     OPT_T, "time-series", 0 },
    { PLOT,     OPT_Y, "single-yaxis", 0 },
    { PLOT,     OPT_Z, "dummy", 0 },
    { PLOT,     OPT_N, "band", 2 },
    { PLOT,     OPT_J, "band-style", 2 },
    { PLOT,     OPT_W, "font", 2 },
    { PLOT,     OPT_L, "ylogscale", 1 },
    { PRINT,    OPT_O, "byobs", 0 },
    { PRINT,    OPT_L, "list", 0 },
    { PRINT,    OPT_D, "no-dates", 0 },
    { PRINT,    OPT_U, "numeric", 0 },
    { PRINT,    OPT_M, "midas", 0 },
    { PRINT,    OPT_C, "complex", 0 },
    { PRINT,    OPT_T, "tree", 0 },
    { PRINT,    OPT_R, "range", 2 },
    { PRINT,    OPT_X, "data-only", 0 },
    { PROBIT,   OPT_P, "p-values", 0 },
    { PROBIT,   OPT_R, "robust", 0 },
    { PROBIT,   OPT_C, "cluster", 2 },
    { PROBIT,   OPT_V, "verbose", 0 },
    { PROBIT,   OPT_E, "random-effects", 0 },
    { PROBIT,   OPT_G, "quadpoints", 2 },
    { PROBIT,   OPT_B, "bootstrap", 1 },
    { PROBIT,   OPT_S, "estrella", 0 },
    { QLRTEST,  OPT_L, "limit-to", 2 },
    { QQPLOT,   OPT_R, "raw", 0 },
    { QQPLOT,   OPT_Z, "z-scores", 0 },
    { QUANTREG, OPT_I, "intervals", 1 },
    { QUANTREG, OPT_N, "no-df-corr", 0 },
    { QUANTREG, OPT_R, "robust", 0 },
    { QUIT,     OPT_X, "exit", 0 },
    { RESET,    OPT_C, "cubes-only", 0 },
    { RESET,    OPT_R, "squares-only", 0 },
    { RESET,    OPT_I, "silent", 0 },
    { RESTRICT, OPT_B, "bootstrap", 1 },
    { RESTRICT, OPT_F, "full", 0 },
    { RESTRICT, OPT_J, "jitter", 0 },
    { RESTRICT, OPT_V, "verbose", 0 },
    { RESTRICT, OPT_L, "lbfgs", 0 },
    { RESTRICT, OPT_N, "no-scaling", 0 },
    { RESTRICT, OPT_S, "silent", 0 },
    { RESTRICT, OPT_W, "wald", 0 },
    { RMPLOT,   OPT_T, "trim", 0 },
    { RUNS,     OPT_D, "difference", 0 },
    { RUNS,     OPT_E, "equal", 0 },
    { SCATTERS, OPT_O, "with-lines", 0 },
    { SCATTERS, OPT_T, "time-series", 0 },
    { SCATTERS, OPT_U, "output", 2 },
    { SCATTERS, OPT_X, "matrix", 2 },
    { SCATTERS, OPT_K, "tweaks", 2 },
    { SET,      OPT_F, "from-file", 2 },
    { SET,      OPT_T, "to-file", 2 },
    { SETINFO,  OPT_C, "continuous", 0 },
    { SETINFO,  OPT_D, "discrete", 0 },
    { SETINFO,  OPT_I, "description", 2 },
    { SETINFO,  OPT_G, "graph-name", 2 },
    { SETINFO,  OPT_M, "midas", 0 },
    { SETINFO,  OPT_F, "coded", 0 },
    { SETINFO,  OPT_N, "numeric", 0 },
    { SETOBS,   OPT_C, "stacked-cross-section", 0 },
    { SETOBS,   OPT_P, "panel-vars", 0 },
    { SETOBS,   OPT_R, "restructure", 0 },
    { SETOBS,   OPT_S, "stacked-time-series", 0 },
    { SETOBS,   OPT_T, "time-series", 0 },
    { SETOBS,   OPT_X, "cross-section", 0 },
    { SETOBS,   OPT_N, "special-time-series", 0 },
    { SETOBS,   OPT_G, "panel-groups", 0 },
    { SETOBS,   OPT_I, "panel-time", 0 },
    { SMPL,     OPT_A, "no-all-missing", 0 },
    { SMPL,     OPT_B, "preserve-panel", 0 },
    { SMPL,     OPT_B, "balanced", 0 }, /* alias */
    { SMPL,     OPT_C, "contiguous", 0 },
    { SMPL,     OPT_D, "dates", 0 },
    { SMPL,     OPT_F, "full", 0 },
    { SMPL,     OPT_O, "dummy", 0 },
    { SMPL,     OPT_M, "no-missing", 0 },
    { SMPL,     OPT_N, "random", 0 },
    { SMPL,     OPT_P, "replace", 0 },
    { SMPL,     OPT_R, "restrict", 0 },
    { SMPL,     OPT_T, "permanent", 0 },
    { SMPL,     OPT_U, "unit", 0 },
    { SMPL,     OPT_X, "time", 0 },
    { SPEARMAN, OPT_V, "verbose", 0 },
    { SQUARE,   OPT_O, "cross", 0 },
    { STDIZE,   OPT_C, "center-only", 0 },
    { STDIZE,   OPT_N, "no-df-corr", 0 },
    { STORE,    OPT_A, "matrix", 2 },
    { STORE,    OPT_D, "database", 0 },
    { STORE,    OPT_E, "comment", 2 },
    { STORE,    OPT_F, "overwrite", 0 },
    { STORE,    OPT_G, "dat", 0 },
    { STORE,    OPT_I, "decimal-comma", 0 },
    { STORE,    OPT_J, "jmulti", 0 },
    { STORE,    OPT_L, "lcnames", 0 },
    { STORE,    OPT_M, "gnu-octave", 0 },
    { STORE,    OPT_N, "no-header", 0 },
    { STORE,    OPT_P, "preserve-strvals", 0 },
    { STORE,    OPT_R, "gnu-R", 0 },
    { STORE,    OPT_X, "omit-obs", 0 },
    { STORE,    OPT_Z, "gzipped", 1 },
    { SUMMARY,  OPT_B, "by", 2 },
    { SUMMARY,  OPT_S, "simple", 0 },
    { SUMMARY,  OPT_W, "weights", 2 },
    { SUMMARY,  OPT_X, "matrix", 2 },
    { SYSTEM,   OPT_I, "iterate", 0 },
    { SYSTEM,   OPT_V, "verbose", 0 },
    { SYSTEM,   OPT_R, "robust", 0 },
    { SYSTEM,   OPT_N, "no-df-corr", 0 },
    { TEXTPLOT, OPT_O, "one-scale", 0 },
    { TEXTPLOT, OPT_S, "time-series", 0 },
    { TEXTPLOT, OPT_T, "tall", 0 },
    { TOBIT,    OPT_L, "llimit", 2 },
    { TOBIT,    OPT_M, "rlimit", 2 },
    { TOBIT,    OPT_G, "opg", 0 },
    { TOBIT,    OPT_R, "robust", 0 },
    { TOBIT,    OPT_C, "cluster", 2 },
    { TOBIT,    OPT_V, "verbose", 0 },
    { VAR,      OPT_D, "seasonals", 0 },
    { VAR,      OPT_F, "variance-decomp", 0 },
    { VAR,      OPT_H, "robust-hac", 0 },
    { VAR,      OPT_I, "impulse-responses", 0 },
    { VAR,      OPT_L, "lagselect", 0 },
    { VAR,      OPT_N, "nc", 0 },
    { VAR,      OPT_R, "robust", 0 },
    { VAR,      OPT_T, "trend", 0 },
    { VAR,      OPT_S, "silent", 0 },
    { VAR,      OPT_M, "minlag", 2 },
    { VARLIST,  OPT_A, "accessors", 0 },
    { VARLIST,  OPT_S, "scalars", 0 },
    { VARLIST,  OPT_T, "type", 2 },
    { VARLIST,  OPT_D, "debug" },
    { VECM,     OPT_A, "crt", 0 },
    { VECM,     OPT_D, "seasonals", 0 },
    { VECM,     OPT_F, "variance-decomp", 0 },
    { VECM,     OPT_I, "impulse-responses", 1 },
    { VECM,     OPT_N, "nc", 0 },
    { VECM,     OPT_R, "rc", 0 },
    { VECM,     OPT_C, "uc", 0 },
    { VECM,     OPT_T, "ct", 0 },
    { VECM,     OPT_V, "verbose", 0 },
    { VECM,     OPT_S, "silent", 0 },
    { WLS,      OPT_R, "robust", 0 },
    { WLS,      OPT_C, "cluster", 2 },
    { WLS,      OPT_Z, "allow-zeros", 0 },
    { XTAB,     OPT_C, "column", 0 },
    { XTAB,     OPT_X, "matrix", 2 },
    { XTAB,     OPT_R, "row", 0 },
    { XTAB,     OPT_Z, "zeros", 0 },
    { XTAB,     OPT_T, "tex", 1 },
    { XTAB,     OPT_N, "no-totals", 0 },
    { XTAB,     OPT_E, "equal", 0 },
    { XTAB,     OPT_F, "no-fisher", 0 },
    { 0,        0L,    NULL, 0 }
};

/* retrieve an integer result directly */

//~ int generate_int (const char *s, DATASET *dset, int *err)
//~ {
    //~ double x = generate_scalar_full(s, dset, NULL, err);
    //~ int ret = -1;

    //~ if (!*err) {
	//~ ret = gretl_int_from_double(x, err);
    //~ }

    //~ return ret;
//~ }



//~ const char *get_longopt (int ci, gretlopt opt)
//~ {
    //~ int i, got_ci = 0;

    //~ for (i=0; gretl_opts[i].ci; i++) {
        //~ if (gretl_opts[i].ci == ci) {
            //~ if (gretl_opts[i].o == opt) {
                //~ return gretl_opts[i].longopt;
            //~ }
            //~ got_ci = 1;
        //~ } else if (got_ci) {
            //~ break;
        //~ }
    //~ }

    //~ return "??";
//~ }

/* retrieve a scalar result directly */
//~ double generate_scalar_full (const char *s, DATASET *dset,
				    //~ PRN *prn, int *err)
//~ {
    //~ parser p;
    //~ double x = NADBL;

    //~ *err = realgen(s, &p, dset, prn, P_PRIV | P_ANON, NUM);

    //~ if (!*err) {
	//~ if (p.ret->t == MAT) {
	    //~ gretl_matrix *m = p.ret->v.m;

	    //~ if (gretl_matrix_is_scalar(m)) {
		//~ x = p.ret->v.m->val[0];
	    //~ } else if (!gretl_is_null_matrix(m)) {
		//~ fprintf(stderr, "generate_scalar: got %d x %d matrix\n",
			//~ m->rows, m->cols);
		//~ *err = E_TYPES;
	    //~ }
	//~ } else if (p.ret->t == NUM) {
	    //~ x = p.ret->v.xval;
	//~ } else {
	    //~ *err = E_TYPES;
	//~ }
    //~ } else if (*err == 1) {
	//~ *err = E_PARSE;
    //~ }

    //~ gen_cleanup(&p);

    //~ return x;
//~ }

//~ int realgen (const char *s, parser *p, DATASET *dset, PRN *prn,
             //~ int flags, int targtype)
//~ {
//~ #if LHDEBUG || EDEBUG || AUX_NODES_DEBUG
    //~ fprintf(stderr, "\n*** realgen: task = %s, depth %d, priv %d\n",
            //~ (flags & P_COMPILE)? "compile" : (flags & P_EXEC)? "exec" : "normal",
            //~ gretl_function_depth(), (flags & P_PRIV)? 1 : 0);
    //~ if (s != NULL) {
        //~ fprintf(stderr, "targ=%d (%s), input='%s'\n", targtype,
                //~ (targtype < PUNCT_MAX)? gretl_type_get_name(targtype) :
                //~ getsymb(targtype), s);
    //~ }
//~ #endif

    //~ if (flags & P_EXEC) {
//~ #if EDEBUG
        //~ fprintf(stderr, "*** printing p->tree (before reinit)\n");
        //~ print_tree(p->tree, p, 0, 0);
//~ #endif
        //~ parser_reinit(p, dset, prn);
        //~ if (p->err) {
            //~ fprintf(stderr, "error in parser_reinit\n");
            //~ goto gen_finish;
        //~ } else if (p->op == INC || p->op == DEC) {
            //~ /* more or less a no-op: the work is done by
               //~ save_generated_var()
            //~ */
            //~ goto gen_finish;
        //~ } else {
            //~ goto starteval;
        //~ }
    //~ } else {
        //~ int done = 0;

        //~ parser_init(p, s, dset, prn, flags, targtype, &done);
        //~ if (p->err) {
            //~ if (gretl_function_depth() == 0) {
                //~ errmsg(p->err, prn);
            //~ }
            //~ goto gen_finish;
        //~ } else if (done) {
            //~ goto gen_finish;
        //~ }
    //~ }

//~ #if EDEBUG
    //~ fprintf(stderr, "after parser %s, p->err = %d (decl? %s)\n",
            //~ (flags & P_EXEC)? "reinit" : "init", p->err,
            //~ (p->flags & P_DECL)? "yes" : "no");
//~ #endif

    //~ if (p->flags & P_DECL) {
        //~ /* check validity of declaration(s) */
        //~ decl_check(p, flags);
        //~ goto gen_finish;
    //~ }

    //~ if (p->op == INC || p->op == DEC) {
        //~ /* implemented via save_generated_var() */
        //~ goto gen_finish;
    //~ }

    //~ /* fire up the lexer */
    //~ lex(p);
    //~ if (p->err) {
//~ #if EDEBUG
        //~ fprintf(stderr, "realgen %p ('%s'): got on lex() error %d\n",
                //~ (void *) p, s, p->err);
//~ #endif
        //~ goto gen_finish;
    //~ }

    //~ /* build the syntax tree */
    //~ p->tree = expr(p);
    //~ if (p->err) {
        //~ goto gen_finish;
    //~ }

//~ #if EDEBUG
    //~ if (p->tree != NULL) {
        //~ fprintf(stderr, "realgen: p->tree at %p, type %d (%s)\n", (void *) p->tree,
                //~ p->tree->t, getsymb(p->tree->t));
    //~ }
    //~ if (p->ch == '\0') {
        //~ fprintf(stderr, " p->ch = NUL, p->sym = %d\n", p->sym);
    //~ } else {
        //~ fprintf(stderr, " p->ch = '%c', p->sym = %d\n", p->ch, p->sym);
    //~ }
//~ #endif

    //~ if (p->sym != EOT || p->ch != 0) {
        //~ int c = p->ch;

        //~ if (c == ' ') {
            //~ c = 0;
        //~ } else if (c != 0) {
            //~ parser_ungetc(p);
            //~ c = p->ch;
        //~ }
        //~ context_error(c, p, "realgen");
        //~ goto gen_finish;
    //~ }

    //~ if (flags & P_NOEXEC) {
        //~ /* we're done at this point */
        //~ goto gen_finish;
    //~ }

    //~ if (!p->err) {
        //~ /* set P_UFRET here if relevant */
        //~ maybe_set_return_flags(p);
    //~ }

 //~ starteval:

//~ #if EDEBUG
    //~ if (flags & P_EXEC) {
        //~ fprintf(stderr, "*** printing p->tree (about to start eval)\n");
        //~ print_tree(p->tree, p, 0, 0);
    //~ }
//~ #endif

    //~ if (autoreg(p)) {
        //~ /* e.g. y = b*y(-1) : evaluate dynamically */
        //~ double *y = p->dset->Z[p->lh.vnum];
        //~ const double *x;
        //~ int t, initted = 0;

        //~ for (t=p->dset->t1; t<=p->dset->t2 && !p->err; t++) {
            //~ /* initialize for this observation */
            //~ p->obs = t;
            //~ if (dataset_is_panel(p->dset) && t % p->dset->pd == 0) {
                //~ initted = 0;
            //~ }
            //~ p->ret = eval(p->tree, p);
            //~ if (p->ret != NULL && p->ret->t == SERIES) {
                //~ x = p->ret->v.xvec;
//~ #if EDEBUG
                //~ autoreg_genr_report(x, y, initted, p);
//~ #endif
                //~ if (!initted && na(x[t])) {
                    //~ ; /* don't overwrite initializer */
                //~ } else {
                    //~ if (p->op == B_ASN) {
                        //~ y[t] = x[t];
                    //~ } else {
                        //~ y[t] = xy_calc(y[t], x[t], p->op, SERIES, p);
                    //~ }
                    //~ initted = 1;
                //~ }
            //~ } else {
                //~ autoreg_error(p, t);
            //~ }
            //~ if (t == p->dset->t1) {
                //~ p->flags &= ~P_START;
            //~ }
        //~ }
    //~ } else {
        //~ /* standard non-dynamic evaluation */
        //~ p->ret = eval(p->tree, p);
    //~ }

    //~ if (p->flags & P_EXEC) {
        //~ p->callcount += 1;
    //~ }

 //~ gen_finish:

    //~ if (p->errprn != NULL) {
        //~ /* Pick and forward any error message that may not be
           //~ seen if realgen was invoked with a NULL value for
           //~ the printer @prn.
        //~ */
        //~ const char *buf = gretl_print_get_buffer(p->errprn);

        //~ if (buf != NULL && *buf != '\0') {
            //~ gretl_errmsg_set(buf);
        //~ }
        //~ gretl_print_destroy(p->errprn);
        //~ p->errprn = NULL;
    //~ }

//~ #if EDEBUG
    //~ fprintf(stderr, "realgen: at finish, err = %d\n", p->err);
//~ # if EDEBUG > 1
    //~ printnode(p->ret, p, 0);
    //~ pputc(prn, '\n');
//~ # endif
//~ #endif

    //~ return p->err;
//~ }

