#include <common.h>
__import glib;

/**
 * gretl_matrix:
 * @rows: number of rows in matrix
 * @cols: number of columns
 * @val: flat array of double-precision values
 *
 * The basic libgretl matrix type; #gretl_vector is an alias
 * that can be used for matrices with @rows or @cols = 1.
 */
typedef struct gretl_matrix_ {
    int rows;
    int cols;
    double *val;
    double _Complex *z; /* was "complex" */
    int is_complex;
    /*< private >*/
    matrix_info *info;
} gretl_matrix;


/**
 * SECTION:gretl_matrix
 * @short_description: construct and manipulate matrices
 * @title: Matrices
 * @include: libgretl.h
 *
 * Libgretl implements most of the matrix functionality that is
 * likely to be required in econometric calculation.  For basics
 * such as decomposition and inversion we use LAPACK as the
 * underlying engine.
 *
 * To get yourself a gretl matrix, use gretl_matrix_alloc() or
 * one of the more specialized constructors; to free such a
 * matrix use gretl_matrix_free().
 */

typedef struct gretl_matrix_block_ {
    int n;
    double *val;
    gretl_matrix **matrix;
} gretl_matrix_block;

struct int_limits {
    SetKey key;
    int min;
    int max;
};

/* matrix metadata struct, not allocated by default */

struct matrix_info_ {
    int t1;
    int t2;
    char **colnames;
    char **rownames;
};

typedef struct matrix_subspec_ matrix_subspec;

union msel {
    int range[2];
    gretl_matrix *m;
    char *str;
};

struct matrix_subspec_ {
    int checked;
    SelType ltype, rtype;
    union msel lsel, rsel;
    int *rslice;
    int *cslice;
};

struct named_val {
    double x;
    const char *s;
};

typedef struct set_state_ set_state;

struct set_state_ {
    int flags;
    /* small integer values */
    gint8 optim;                /* code for preferred optimizer */
    gint8 vecm_norm;            /* VECM beta normalization */
    gint8 garch_vcv;            /* GARCH vcv variant */
    gint8 garch_alt_vcv;        /* GARCH vcv variant, robust estimation */
    gint8 arma_vcv;             /* ARMA vcv variant */
    gint8 wildboot_d;           /* distribution for wild bootstrap */
    gint8 fdjac_qual;           /* quality of "fdjac" function */
    gint8 use_qr;               /* off, on or pivot */
    gint8 max_verbose;          /* optimizer verbosity level */
    gint8 hc_version;           /* HCCME version */
    gint8 panel_robust;         /* panel robust vcv estimator */
    gint8 hac_kernel;           /* HAC kernel type */
    gint8 auto_hac_lag;         /* HAC automatic lag-length formula */
    gint8 user_hac_lag;         /* fixed user-set HAC lag length */
    gint8 lbfgs_mem;            /* memory span for L-BFGS-B */
    gint8 quantile_type;        /* Formula for computing quantiles */
    /* potentially larger integers */
    int horizon;                /* for VAR impulse responses */
    int bootrep;                /* bootstrap replications */
    int loop_maxiter;           /* max no. of iterations in non-for loops */
    int bfgs_maxiter;           /* max iterations, BFGS */
    int bfgs_verbskip;          /* BFGS: show one in n iterations  */
    int boot_iters;             /* max iterations, IRF bootstrap */
    int bhhh_maxiter;           /* max iterations, BHHH */
    int rq_maxiter;             /* max iterations for quantreg, simplex */
    int gmm_maxiter;            /* max iterations for iterated GMM */
    /* floating-point values */
    double conv_huge;           /* conventional value for $huge */
    double nls_toler;           /* NLS convergence criterion */
    double bfgs_toler;          /* convergence tolerance, BFGS */
    double bfgs_maxgrad;        /* max acceptable gradient norm, BFGS */
    double bhhh_toler;          /* convergence tolerance, BHHH */
    double qs_bandwidth;        /* bandwidth for QS HAC kernel */
    double nadarwat_trim;       /* multiple of h to use in nadarwat() for trimming */
    /* strings */
    char csv_write_na[8];       /* representation of NA in CSV output */
    char csv_read_na[8];        /* representation of NA in CSV input */
    /* matrices */
    gretl_matrix *initvals;     /* parameter initializer */
    gretl_matrix *initcurv;     /* initial curvature matrix for BFGS */
    gretl_matrix *matmask;      /* mask for series -> matrix conversion */
};

/* struct to hold model results */
typedef struct MODEL_ {
    int ID;                      /* ID number for model */
    int refcount;                /* for saving/deleting */
    int ci;                      /* "command index" -- estimation method */
    gretlopt opt;                /* record of options */
    int t1, t2, nobs;            /* starting observation, ending
                                    observation, and number of obs */
    char *submask;               /* keep track of sub-sample in force
                                    when model was estimated */
    char *missmask;              /* missing observations mask */
    SAMPLE smpl;                 /* numeric start and end of current sample
                                    when model was estimated */
    int full_n;                  /* full length of dataset on estimation */
    int ncoeff, dfn, dfd;        /* number of coefficents; degrees of
                                    freedom in numerator and denominator */
    int *list;                   /* list of variables by ID number */
    int ifc;                     /* = 1 if the equation includes a constant,
				    else = 0 */
    int nwt;                     /* ID number of the weight variable (WLS) */
    int aux;                     /* code representing the sort of
				    auxiliary regression this is (or not) */
    double *coeff;               /* array of coefficient estimates */
    double *sderr;               /* array of estimated std. errors */
    double *uhat;                /* regression residuals */
    double *yhat;                /* fitted values from regression */
    double *xpx;                 /* X'X matrix, in packed form */
    double *vcv;                 /* VCV matrix for coefficient estimates */
    double ess, tss;             /* Error and Total Sums of Squares */
    double sigma;                /* Standard error of regression */
    double rsq, adjrsq;          /* Unadjusted and adjusted R^2 */     
    double fstt;                 /* overall F-statistic */
    double chisq;                /* overall chi-square statistic */
    double lnL;                  /* log-likelihood */
    double ybar, sdy;            /* mean and std. dev. of dependent var. */
    double criterion[C_MAX];     /* array of model selection statistics */
    double dw, rho;              /* Durbin-Watson stat. and estimated 1st
				    order autocorrelation coefficient */
    ARINFO *arinfo;              /* pointer to struct to hold special info for 
				    autoregressive model */ 
    int errcode;                 /* Error code in case of failure */
    char *name;                  /* for use in GUI */
    char *depvar;                /* name of dependent var in special cases */
    int nparams;                 /* number of named model parameters */
    char **params;               /* for named model parameters */
    gint64 esttime;              /* time of estimation */
    int ntests;                  /* number of attached test results */
    ModelTest *tests;            /* attached hypothesis test results */
    DATASET *dataset;            /* for handling models estimated on a
				    sub-sampled portion of the dataset */
    int n_data_items;            /* number of extra data items */
    model_data_item **data_items; /* pointer to additional data */
} MODEL;

typedef struct FITRESID_ {
    int model_ID;   /* ID of model on which forecast is based */
    int asymp;      /* 0/1 flag for asymptotic estimator */
    int std;        /* 0/1 flag for standardized residuals */
    int model_t1;   /* start of model estimation range */
    int method;     /* one of the ForecastMethod options */
    double *actual; /* array of values of dependent variable */
    double *fitted; /* array of fitted values */
    double *resid;  /* array of residuals */
    double *sderr;  /* array of forecast standard errors (or NULL) */
    double sigma;   /* standard error of regression */
    double alpha;   /* for confidence intervals */
    int pmax;       /* if positive, suggested number of decimal places
                       for use in printing */
    int df;         /* degrees of freedom for model */
    int t0;         /* start of pre-forecast data range */
    int t1;         /* start of forecast range */
    int t2;         /* end of forecast range */
    int k;          /* number of steps ahead (method = FC_KSTEP only) */
    int nobs;       /* length of the arrays actual, fitted, resid */
    char depvar[VNAMELEN]; /* name of dependent variable */
} FITRESID;

typedef struct ARINFO_ {
    int *arlist;          /* list of autoregressive lags */
    double *rho;          /* array of autoreg. coeffs. */
    double *sderr;        /* and their standard errors */
} ARINFO;

typedef struct VCVInfo_ VCVInfo;

struct VCVInfo_ {
    int vmaj;        /* general type of VCV (see VCVMajorType) */
    int vmin;        /* variant of general type */
    int order;       /* for use with HAC */
    VCVFlags flags;  /* includes prewhitening */
    double bw;       /* for use with QS HAC kernel */
    char *cv1;       /* name of (first) cluster variable */
    char *cv2;       /* name of second cluster var */
};

typedef struct {} GHashTable;

typedef struct DATASET_ { 
    int v;              /* number of data series, including const */
    int n;              /* number of observations */
    int pd;             /* periodicity or frequency of data */
    int structure;      /* time series, cross section or whatever */
    double sd0;         /* floating-point representation of stobs */
    int t1, t2;         /* start and end of current sample */
    char stobs[OBSLEN];  /* string representation of starting obs */
    char endobs[OBSLEN]; /* string representation of ending obs */
    double **Z;         /* two-dimensional data array */
    char **varname;     /* array of names of series */
    VARINFO **varinfo;  /* array of metadata per series */
    char markers;       /* NO_MARKERS (0), REGULAR MARKERS or DAILY_DATE_STRINGS */
    char modflag;       /* binary flag for dataset modified or not */
    char **S;           /* to hold observation markers */
    char *descrip;      /* to hold info on data sources, etc. */
    char *submask;      /* subsampling mask */
    char *restriction;  /* record of sub-sampling restriction */
    char *padmask;      /* record of padding to re-balance panel data */
    char *mapfile;      /* name of associated map (polygons) file, if any */
    unsigned int rseed; /* resampling seed */
    int auxiliary;      /* 0 for regular dataset, 1 for auxiliary dataset */
    char *pangrps;      /* panel-only: name of series holding group names */
    int panel_pd;       /* panel-only: panel time-series frequency */
    double panel_sd0;   /* panel-only: time-series start */
} DATASET;

typedef struct series_table_ {
    int n_strs;       /* number of strings in table */
    char **strs;      /* saved strings */
    GHashTable *ht;   /* hash table for quick lookup */
    int flags;        /* status flags (above) */
} series_table;

typedef struct Summary_ {
    gretlopt opt;
    int n;
    int weight_var;
    int *misscount;
    int *list;
    double *stats;
    double *mean;
    double *median;
    double *sd;
    double *skew; 
    double *xkurt;
    double *low;
    double *high;
    double *cv;
    double *perc05;
    double *perc95;
    double *iqr;
    double sw;
    double sb;
} Summary;

typedef struct FreqDist_ {
    char varname[VNAMELEN];  /* for ID purposes */
    int discrete;            /* 1 if series contains integers */
    int strvals;             /* 1 if series is string-valued */
    int dist;                /* code for theoretical distribution */
    int numbins;             /* number of bins or intervals */
    double xbar, sdx;        /* mean and std dev of variable */
    double *midpt;           /* array of midpoints of intervals */
    double *endpt;           /* array of endpoints of intervals */
    char **S;                /* array of string values */
    int *f;                  /* frequencies */
    double test;             /* Chi-squared statistic for testing for a
                                Gaussian distribution, or z statistic
			        for testing for Gamma dist. */
    int n;                   /* included observation */
    int t1, t2;              /* sample limits */
} FreqDist;

typedef struct Xtab_ {
    char rvarname[VNAMELEN]; /* name of rows series */
    char cvarname[VNAMELEN]; /* name of cols series */
    char **Sr;               /* array of row string values */
    char **Sc;               /* array of column string values */
    int rows, cols;          /* dimensions of table */
    double *rval, *cval;     /* row and column numeric values */
    int *rtotal, *ctotal;    /* marginal totals */
    int **f;                 /* array of frequencies */
    int n, missing;          /* observation counts */
    int t1, t2;              /* sample limits */
    int rstrs;               /* row string-valued flag */
    int cstrs;               /* col string-valued flag */
} Xtab;

/**
 * SECTION:gretl_prn
 * @short_description: gretl printing struct
 * @title: PRN
 * @include: libgretl.h
 *
 * Most libgretl functions that print output call for a
 * pointer-to-PRN as one of their arguments. The PRN type is
 * an opaque structure, constructed and manipulated by the
 * functions listed here. It is used as a means of generalizing
 * the printing operation, which may be to a regular file, a
 * temporary file or a buffer.
 *
 * To get hold of a PRN use gretl_print_new() or one of the more
 * specific constructors, and to free it use gretl_print_destroy().
 * If you want to use a PRN dirctly for printing in your own code, use
 * the functions pprintf(), pputs() and pputc(). These are
 * counterparts to the standard C functions fprintf, fputs and
 * fputc, but note that with pputs and pputc the PRN argument
 * must be given first (unlike fputs and fputc in which the
 * FILE argument goes last).
 *
 * Note that whenever a PRN appears as a function parameter
 * in libgretl it is OK to give a NULL argument: in that case
 * pprintf(), pputs() and pputc() are no-ops.
 */

typedef struct PRN_ {
    FILE *fp;          /* file to which to print, or NULL */
    gzFile fz;         /* gzipped file target, or NULL */
    char *buf;         /* buffer to which to print, or NULL */
    size_t bufsize;    /* allocated size of buffer */
    size_t blen;       /* string length of buffer */
    int savepos;       /* saved position in stream or buffer */
    GArray *fplist;    /* stack for use with output redirection */
    PrnFormat format;  /* plain, TeX, RTF */
    gint8 fixed;       /* non-zero for fixed-size buffer */
    gint8 gbuf;        /* non-zero for buffer obtained via GLib */
    guint8 nlcount;    /* count of trailing newlines */
    char delim;        /* CSV field delimiter */
    char *fname;       /* temp file name, or NULL */
} PRN;

struct fpinfo_ {
    FILE *fp;      /* stream to which we're printing */
    int level;     /* level of redirection */
    gchar *fname;  /* name of file or NULL */
    gchar *strvar; /* associated string variable or NULL */
};

typedef struct fpinfo_ fpinfo;

typedef struct {} gzFile;

typedef struct gzFile_ {} gzFile;

typedef struct VARINFO_ {
    char *label;
    char display_name[MAXDISP];
    char parent[VNAMELEN];
    VarFlags flags;
    char compact_method;
    gint64 mtime;
    short transform;    /* note: command index of transform */
    short lag;
    short stack_level;
    short midas_period;
    char midas_freq;
    short orig_pd;
    series_table *st;
} VARINFO;

typedef struct fncall_ {
    ufunc *fun;      /* the function called */
    int argc;        /* argument count */
    int orig_v;      /* number of series defined on entry */
    fn_arg *args;    /* argument array */
    int *ptrvars;    /* list of pointer arguments */
    int *listvars;   /* list of series included in a list argument */
    GList *lists;    /* list of names of list arguments */
    char *retname;   /* name of return value (or dummy string) */
    GretlType rtype; /* return type (when not fixed in advance) */
    obsinfo obs;     /* sample info */
    FCFlags flags;   /* indicators for recursive call, etc */
    fn_line *line;   /* currently executing line */
    linegen *lgen;   /* array of per-line genrs */
    int n_lgen;      /* number of linegens */
} fncall;

/* structure representing a user-defined function */

typedef struct ufunc_ {
    char name[FN_NAMELEN]; /* identifier */
    fnpkg *pkg;            /* pointer to parent package, or NULL */
    fncall *call;          /* pointer to current call, or NULL */
    GList *calls;          /* for use in recursion */
    int pkg_role;          /* printer, plotter, etc. */
    UfunAttrs flags;       /* private, plugin, etc. */
    int line_idx;          /* current line index (compiling) */
    int n_lines;           /* number of lines of code */
    fn_line *lines;        /* array of lines of code */
    int n_params;          /* number of parameters */
    fn_param *params;      /* parameter info array */
    int rettype;           /* return type (if any) */
} ufunc;

typedef struct fn_arg_ {
    char type;            /* argument type */
    char shifted;         /* level was shifted for execution */
    char *upname;         /* name of supplied arg at caller level */
    user_var *uvar;       /* reference to "parent", if any */
    union {
        int idnum;        /* named series arg (series ID) */
        double x;         /* scalar arg */
        double *px;       /* anonymous series arg */
        gretl_matrix *m;  /* matrix arg */
        char *str;        /* string arg */
        int *list;        /* list arg */
        gretl_bundle *b;  /* anonymous bundle pointer */
        gretl_array *a;   /* array argument */
    } val;
} fn_arg;

/* structure representing sample information at start of
   a function call */

typedef struct obsinfo_ {
    int structure;      /* time-series, etc. */
    int pd;             /* data frequency */
    int t1, t2;         /* starting and ending observations */
    int added;          /* number of observations added within function */
    char stobs[OBSLEN]; /* string representation of starting obs */
    int panel_pd;       /* panel time frequency, if applicable */
    double panel_sd0;   /* panel time starting point */
} obsinfo;

/* local symbol for structure representing a line of a
   user-defined function */
typedef struct stmt_ fn_line;

/* Note: the FC_PRESERVE flag is set when we want to allow for repeated
   calls to a given function. In the non-recursive case we just reuse
   the stored fncall struct. In the recursive case we need to ensure
   there's a distinct call struct for each depth of function execution,
   and the function then acquires a GList of fncalls.
*/

struct linegen_ {
    int idx;
    GENERATOR *genr;
};

typedef struct linegen_ linegen;

/* structure representing a function package */

typedef struct fnpkg_ {
    char name[FN_NAMELEN]; /* package name */
    char *fname;      /* filename */
    char *author;     /* author's name */
    char *email;      /* author's email address */
    char *version;    /* package version string */
    char *date;       /* last revision date */
    char *descrip;    /* package description */
    char *help;       /* package help text */
    char *gui_help;   /* GUI-specific help (optional) */
    char *Rdeps;      /* R dependencies (if any) */
    char *sample;     /* sample caller script */
    char *help_fname;     /* filename: package help text */
    char *gui_help_fname; /* filename: GUI-specific help text */
    char *sample_fname;   /* filename: sample caller script */
    char *tags;       /* tag string(s) */
    char *label;      /* for use in GUI menus */
    char *mpath;      /* menu path in GUI */
    int minver;       /* minimum required gretl version */
    char uses_subdir; /* lives in subdirectory (0/1) */
    char prechecked;  /* already checked for data requirement */
    char data_access; /* wants access to full data range */
    DataReq dreq;     /* data requirement */
    int modelreq;     /* required model type, if applicable */
    ufunc **pub;      /* pointers to public interfaces */
    ufunc **priv;     /* pointers to private functions */
    int n_pub;        /* number of public functions */
    int n_priv;       /* number of private functions */
    char overrides;   /* number of overrides of built-in functions */
    char **datafiles; /* names of packaged data files */
    char **depends;   /* names of dependencies */
    char *provider;   /* name of "provider" package, if applicable */
    int n_files;      /* number of data files */
    int n_depends;    /* number of dependencies */
    void *editor;     /* for GUI use */
} fnpkg;

/* structure representing a parameter of a user-defined function */

typedef struct fn_param_ {
    char *name;     /* the name of the parameter */
    char type;      /* its type */
    char *descrip;  /* its description */
    char **labels;  /* value labels, if applicable */
    int nlabels;    /* number of value labels */
    char flags;     /* additional information (e.g. "const" flag) */
    double deflt;   /* default value */
    double min;     /* minimum value (scalar parameters only) */
    double max;     /* maximum value (scalar parameters only) */
    double step;    /* step increment (scalars only) */
} fn_param;

/**
 * gretl_array:
 *
 * An opaque type; use the relevant accessor functions.
 */

typedef struct gretl_array_ {
    GretlType type;  /* type of data */
    int n;           /* number of elements */
    void **data;     /* actual data array */
    double *mdata;   /* for matrix block */
} gretl_array;

typedef struct parser_ {
    const char *input; /* complete input string */
    const char *point; /* remaining unprocessed input */
    const char *rhs;   /* for use in labelling */
    DATASET *dset;     /* convenience pointer to dataset */
    PRN *prn;          /* for printing messages */
    PRN *errprn;       /* for storing error message in case @prn is NULL */
    genflags flags;    /* various attributes (see @genflags above) */
    int targ;          /* target type */
    int op;            /* assignment operator (possibly inflected) */
    lhinfo lh;  /* left-hand side info */
    NODE *lhtree;      /* LHS syntax tree, if needed */
    NODE *lhres;       /* result of eval() on @lhtree */
    NODE *tree;        /* RHS syntax tree */
    NODE *ret;         /* result of eval() on @tree */
    /* below: parser state variables */
    NODE *aux;         /* convenience pointer to current auxiliary node */
    int callcount;
    int dset_n;
    int obs;
    int sym;
    int upsym;
    int ch;
    double xval;
    int idnum;
    char *idstr;
    void *data;
    int err;
} parser;

typedef struct lhinfo_ {
    int t;                 /* type of pre-existing LHS variable, if any */
    char name[VNAMELEN];   /* name of LHS variable */
    char *label;           /* descriptive string for series */
    series_table *stab;    /* holds string values for series */
    int vnum;              /* ID number of pre-existing LHS series */
    user_var *uv;          /* address of pre-existing LHS variable */
    char *expr;            /* expression on left */
    GretlType gtype;       /* gretl type of LHS array, if any, or
			      of LHS bundle member */
    gretl_matrix *mret;    /* matrix output (possibly under bundle or array) */
} lhinfo;

typedef struct node_ {
    gint16 t;        /* type identifier */
    guint8 flags;    /* AUX_NODE etc., see above */
    int vnum;        /* associated series ID number */
    char *vname;     /* associated variable name */
    user_var *uv;    /* associated named variable */
    union val v;     /* value (of whatever type) */
    NODE *L, *M, *R; /* up to three child nodes */
    NODE *aux;       /* auxiliary (result) node */
    NODE *parent;    /* parent node (or NULL) */
    int refcount;    /* reference counter, used by aux nodes */
} node;

struct branchn {
    int n_nodes;
    NODE **n;
};

union val {
    struct branchn bn;
    int idnum;
    char *str;
    double xval;
    double *xvec;
    int *ivec;
    gretl_matrix *m;
    matrix_subspec *mspec;
    gretl_bundle *b;
    gretl_array *a;
    void *ptr;
};

typedef struct ModelTest_ {
    int type;
    int order;
    char *param;
    unsigned char teststat;
    int dfn;
    double dfd;
    double value;
    double pvalue;
    double crit;
    double alpha;
    gretlopt opt;
} ModelTest;

typedef struct model_data_item_ {
    char *key;
    void *ptr;
    int type;
    size_t size;
    void (*destructor) (void *);
} model_data_item;

typedef struct VMatrix_ {
    int ci;
    int dim;
    int t1, t2, n;
    char **names;
    double *vec;
    double *xbar;
    double *ssx;
    int *list;
    int missing;
} VMatrix;

typedef struct SAMPLE_ {
    int t1;
    int t2;
    unsigned int rseed;
} SAMPLE;

typedef enum {
    C_AIC,
    C_BIC,
    C_HQC,
    C_MAX
} ModelSelCriteria;

typedef enum {
    CLEAR_FULL,           /* fully clear the dataset */
    CLEAR_SUBSAMPLE       /* dataset is sub-sampled: clear partially */
} DataClearCode;

typedef enum {
    DS_NONE,
    DS_ADDOBS,
    DS_COMPACT,
    DS_EXPAND,
    DS_TRANSPOSE,
    DS_SORTBY,
    DS_DSORTBY,
    DS_RESAMPLE,
    DS_CLEAR,
    DS_RENUMBER,
    DS_INSOBS,
    DS_PAD_DAILY,
    DS_UNPAD_DAILY
} DatasetOp;

typedef enum {
    DS_COPY_VALUES,
    DS_GRAB_VALUES
} DataCopyFlag;

enum ts_codes {
    CROSS_SECTION,
    TIME_SERIES,
    STACKED_TIME_SERIES,
    STACKED_CROSS_SECTION,
    PANEL_UNKNOWN,
    PANEL_SIDE_BY_SIDE,
    SPECIAL_TIME_SERIES,
    STRUCTURE_UNKNOWN
};

typedef enum {
    OPT_NONE = 0,
    OPT_A = 1 <<  0,
    OPT_B = 1 <<  1,
    OPT_C = 1 <<  2,
    OPT_D = 1 <<  3,
    OPT_E = 1 <<  4,
    OPT_F = 1 <<  5,
    OPT_G = 1 <<  6,
    OPT_H = 1 <<  7,
    OPT_I = 1 <<  8,
    OPT_J = 1 <<  9,
    OPT_K = 1 << 10,
    OPT_L = 1 << 11,
    OPT_M = 1 << 12,
    OPT_N = 1 << 13,
    OPT_O = 1 << 14,
    OPT_P = 1 << 15,
    OPT_Q = 1 << 16,
    OPT_R = 1 << 17,
    OPT_S = 1 << 18,
    OPT_T = 1 << 19,
    OPT_U = 1 << 20,
    OPT_V = 1 << 21,
    OPT_W = 1 << 22,
    OPT_X = 1 << 23,
    OPT_Y = 1 << 24,
    OPT_Z = 1 << 25,
    OPT_a = 1 << 26,
    OPT_b = 1 << 27,
    OPT_i = 1 << 28,
    OPT_UNSET = 1 << 30
} gretlopt;

typedef enum {
    ADD = 1,
    ADF,
    ANOVA,
    APPEND,
    AR,  
    AR1,
    ARCH,
    ARMA,
    BDS,
    BIPROBIT,
    BKW,
    BREAK,
    BXPLOT,
    CHOW,
    CLEAR,
    COEFFSUM,
    COINT,
    COINT2,
    CONTINUE,
    CORR,     
    CORRGM,   
    CUSUM,
    DATA,
    DATAMOD,
    DELEET,
    DIFF,
    DIFFTEST,
    DISCRETE,
    DPANEL,
    DUMMIFY,
    DURATION,
    ELIF,
    ELSE,
    END,
    ENDIF,
    ENDLOOP,
    EQNPRINT, 
    EQUATION,
    ESTIMATE,
    EVAL,
    FCAST,
    FLUSH,
    FOREIGN,
    FRACTINT,
    FREQ, 
    FUNC,
    FUNCERR,
    GARCH,
    GENR,  
    GMM,
    GNUPLOT,
    GPBUILD,
    GRAPHPG,
    GRIDPLOT,
    HECKIT,
    HELP,
    HFPLOT,
    HSK,
    HURST,
    IF,
    INCLUDE,
    INFO,
    INTREG,
    JOIN,
    KDPLOT,
    KPSS,
    LABELS, 
    LAD,
    LAGS,    
    LDIFF,
    LEVERAGE,
    LEVINLIN,
    LOGISTIC,
    LOGIT,
    LOGS,
    LOOP,
    MAHAL,
    MAKEPKG,
    MARKERS,
    MEANTEST,
    MIDASREG,
    MLE,
    MODELTAB,
    MODPRINT,
    MODTEST,
    MPI,
    MPOLS,
    NEGBIN,
    NLS,
    NORMTEST,
    NULLDATA,
    OLS,     
    OMIT,
    OPEN,
    ORTHDEV,
    OUTFILE,
    PANEL,
    PANPLOT,
    PANSPEC,
    PCA,
    PERGM,
    PLOT,    
    POISSON,
    PRINT, 
    PRINTF,
    PROBIT,
    PVAL, 
    QUANTREG,
    QLRTEST,
    QQPLOT,
    QUIT,
    RENAME,
    RESET,
    RESTRICT,
    RMPLOT,
    RUN,
    RUNS,
    SCATTERS,
    SDIFF,
    SET,
    SETINFO,
    SETOBS,
    SETOPT,
    SETMISS,
    SHELL,  
    SMPL,
    SPEARMAN,
    SQUARE,
    STDIZE,
    STORE, 
    SUMMARY,
    SYSTEM,
    TABPRINT,
    TEXTPLOT,
    TOBIT,
    IVREG,
    VAR,
    VARLIST,
    VARTEST,
    VECM,
    VIF,
    WLS,
    XCORRGM,
    XTAB,
    FUNCRET,
    CATCH,
    PKG,
    NC
} GretlCmdIndex;

typedef enum {
    NORM_PHILLIPS,
    NORM_DIAG,
    NORM_FIRST,
    NORM_NONE,
    NORM_MAX
} VECMnorm;

typedef enum {
    OPTIM_AUTO,
    OPTIM_BFGS,
    OPTIM_NEWTON,
    OPTIM_MAX
} OptimCode;

typedef enum {
    HAC_REFUSE,
    HAC_ES,
    HAC_AM
} HACopts;

typedef enum {
    ARELLANO,
    BECK_KATZ,
    DRISCOLL_KRAAY
} PanelRobust;

enum {
    AUTO_LAG_STOCK_WATSON,
    AUTO_LAG_WOOLDRIDGE,
    AUTO_LAG_NEWEYWEST
};

enum {
    CAT_BEHAVE = 1,
    CAT_NUMERIC,
    CAT_RNG,
    CAT_ROBUST,
    CAT_TS,
    CAT_SPECIAL
};

typedef enum {
    SV_ALL,
    SV_INT,
    SV_DOUBLE
} SVType;

struct CoeffIntervals_ {
    int asy;
    int ncoeff;
    double alpha;
    double t;
    char **names;
    double *coeff;
    double *maxerr;
    int df;
    gretlopt opt;
};

typedef enum {
    VCV_CLASSICAL,
    VCV_HC,
    VCV_HAC,
    VCV_ML,
    VCV_PANEL,
    VCV_RQ,
    VCV_CLUSTER
} VCVMajorType;

typedef enum {
    ML_UNSET,
    ML_HESSIAN,
    ML_IM,
    ML_OP,
    ML_QML,
    ML_BW,
    ML_HAC,
    ML_VCVMAX
} MLVCVType;

typedef enum {
    KERNEL_BARTLETT,
    KERNEL_PARZEN,
    KERNEL_QS,
    KERNEL_MAX
} HACKernel;

/**
 * ModelAuxCode:
 * @AUX_NONE: not an auxiliary regression
 * @AUX_SQ: nonlinearity test (squared terms)
 * @AUX_LOG: nonlinearity test (log terms)
 * @AUX_CHOW: Chow test
 * @AUX_ADD: LM test regression for added variables
 * @AUX_AR: autocorrelation test
 * @AUX_ARCH: ARCH test
 * @AUX_WHITE: heteroskedasticity (White's test)
 * @AUX_COINT: cointegration test
 * @AUX_DF: Dickey-Fuller test
 * @AUX_ADF: augmented Dickey-Fuller test
 * @AUX_KPSS: KPSS unit-root test
 * @AUX_OMIT: unused
 * @AUX_RESET: Ramsey's RESET
 * @AUX_SYS: single equation from multivariate system
 * @AUX_VAR: single equation from VAR system
 * @AUX_VECM: single equation from VECM system
 * @AUX_JOHANSEN: Johansen cointegration test
 * @AUX_GROUPWISE: test for groupwise heteroskedasticity
 * @AUX_HET_1: Pesaran-Taylor HET_1 test
 * @AUX_BP: Breusch-Pagan heteroskedastcity test
 * @AUX_AUX: auxiliary regression not otherwise specified
 * @AUX_COMFAC: common factor test
 * @AUX_BIPROB: biprobit initializer
 *
 * Symbolic names to keep track of auxiliary regression models,
 * which are estimated either for the purpose of carrying out
 * some sort of diagnostic test or which form part of a
 * multi-equation system.
 */

typedef enum {
    AUX_NONE,
    AUX_SQ,
    AUX_LOG,
    AUX_CHOW,
    AUX_ADD,
    AUX_AR,
    AUX_ARCH,
    AUX_WHITE,
    AUX_COINT,
    AUX_DF,
    AUX_ADF,
    AUX_KPSS,
    AUX_OMIT,
    AUX_RESET,
    AUX_SYS,
    AUX_VAR,
    AUX_VECM,
    AUX_JOHANSEN,
    AUX_GROUPWISE,
    AUX_HET_1,
    AUX_BP,
    AUX_AUX,
    AUX_COMFAC,
    AUX_BIPROB,
} ModelAuxCode;

typedef enum {
    UVAR_ADD = 1,
    UVAR_DELETE
} UvarAction;

typedef enum {
    UV_PRIVATE = 1 << 0,
    UV_SHELL   = 1 << 1,
    UV_MAIN    = 1 << 2,
    UV_NODECL  = 1 << 3,
    UV_NOREPL  = 1 << 4
} UVFlags;

/* structure representing a call to a user-defined function */

typedef enum {
    FC_RECURSING = 1 << 0, /* call to f() has a call to f() in its ancestry */
    FC_PREV_MSGS = 1 << 1, /* record of "messages" setting before call */
    FC_PREV_ECHO = 1 << 2, /* record of "echo" setting before call */
    FC_PRESERVE  = 1 << 3  /* function call should be preserved */
} FCFlags;

typedef enum {
    UFUN_PRIVATE   = 1 << 0, /* is private to a package */
    UFUN_NOPRINT   = 1 << 1, /* offers no printed output */
    UFUN_MENU_ONLY = 1 << 2, /* is GUI-only */
    UFUN_USES_SET  = 1 << 3, /* includes the "set" command */
    UFUN_HAS_FLOW  = 1 << 4  /* includes flow-control (ifs, loops) */
} UfunAttrs;

typedef enum {
    GRETL_PRINT_STDOUT,
    GRETL_PRINT_STDERR,
    GRETL_PRINT_FILE,
    GRETL_PRINT_BUFFER,
    GRETL_PRINT_TEMPFILE,
    GRETL_PRINT_STREAM,
    GRETL_PRINT_GZFILE
} PrnType;

typedef enum {
    P_DISCARD = 1 <<  0, /* compute and print, don't save */
    P_START   = 1 <<  1, /* first round of evaluation */
    P_AUTOREG = 1 <<  2, /* expression is autoregressive */
    P_DECL    = 1 <<  3, /* statement is actually a declaration */
    P_PRIV    = 1 <<  4, /* generating a "private" or internal var */
    P_COMPILE = 1 <<  5, /* compiling the parse tree */
    P_EXEC    = 1 <<  6, /* evaluating a compiled tree */
    P_NATEST  = 1 <<  7, /* testing for NAs in expression */
    P_UFRET   = 1 <<  8, /* returning value generated by user function */
    P_QUIET   = 1 <<  9, /* don't print any messages or labels */
    P_GETSTR  = 1 << 10, /* state: flag acceptance of plain strings */
    P_MMASK   = 1 << 11, /* genr result is masked matrix */
    P_SLICING = 1 << 12, /* state: calculating object slice (temporary) */
    P_LAGPRSE = 1 << 13, /* state: parsing lag spec (temporary) */
    P_DELTAN  = 1 << 14, /* flag for change in series length */
    P_CATCH   = 1 << 15, /* "catch" is in force */
    P_NODECL  = 1 << 16, /* type of result was not specified */
    P_LISTDEF = 1 << 17, /* expression defines a list */
    P_ANON    = 1 << 18, /* generating an anonymous object */
    P_VOID    = 1 << 19, /* function call, no assignment */
    P_NOEXEC  = 1 << 20, /* just compile, don't evaluate */
    P_MSAVE   = 1 << 21, /* trying for reuse of an aux matrix */
    P_OBSVAL  = 1 << 22, /* generating value of observation in series */
    P_ALIASED = 1 << 23, /* state: handling aliased object (temporary) */
    P_AND     = 1 << 24, /* state: working on right-hand term of B_AND */
    P_OR      = 1 << 25, /* state: working on right-hand term of B_OR */
    P_STACK   = 1 << 26, /* executing stack() */
    P_ALTINP  = 1 << 27, /* the input string has been substituted */
    P_OBJQRY  = 1 << 28, /* querying the existence of an object */
    P_PRNLIST = 1 << 29  /* defining a list for "print" */
} genflags;

typedef enum {
    FN_NEEDS_DATA,  /* needs some (any) sort of dataset */
    FN_NEEDS_TS,    /* function requires time-series data */
    FN_NEEDS_QM,    /* function requires quarterly or monthly data */
    FN_NEEDS_PANEL, /* function requires panel data */
    FN_NODATA_OK    /* function does not require a dataset */
} DataReq;

typedef enum {
    BUNDLE_PLAIN,
    BUNDLE_KALMAN
} BundleType;

typedef enum {
    SEL_NULL,    /* nothing supplied */
    SEL_RANGE,   /* integer range p:q provided */
    SEL_MATRIX,  /* selection matrix provided */
    SEL_ALL,     /* comma-separated blank */
    SEL_DIAG,    /* the "diag" dummy constant */
    SEL_UPPER,   /* the "upper" dummy constant */
    SEL_LOWER,   /* the "lower" dummy constant */
    SEL_REAL,    /* the "real" dummy constant */
    SEL_IMAG,    /* the "imag" dummy constant */
    SEL_ELEMENT, /* derived: selection is a single element */
    SEL_CONTIG,  /* derived: selection is contiguous */
    SEL_EXCL,    /* single exclusion (negative index) */
    SEL_SINGLE,  /* derived: degenerate range + null */
    SEL_STR      /* for use with bundles only */
} SelType;

typedef enum {
    GRETL_TEST_ADD,
    GRETL_TEST_ARCH,
    GRETL_TEST_AUTOCORR,
    GRETL_TEST_CHOW,
    GRETL_TEST_CUSUM,
    GRETL_TEST_QLR,
    GRETL_TEST_GROUPWISE,
    GRETL_TEST_LOGS,
    GRETL_TEST_NORMAL,
    GRETL_TEST_OMIT,
    GRETL_TEST_RESET,
    GRETL_TEST_SQUARES,
    GRETL_TEST_WHITES,
    GRETL_TEST_SARGAN,
    GRETL_TEST_IV_HAUSMAN,
    GRETL_TEST_PANEL_HAUSMAN,
    GRETL_TEST_PANEL_F,
    GRETL_TEST_PANEL_BP,
    GRETL_TEST_PANEL_TIMEDUM,
    GRETL_TEST_PANEL_AR,
    GRETL_TEST_HET_1,
    GRETL_TEST_BP,
    GRETL_TEST_CHOWDUM,
    GRETL_TEST_COMFAC,
    GRETL_TEST_INDEP,
    GRETL_TEST_RE,
    GRETL_TEST_WITHIN_F,
    GRETL_TEST_PANEL_WELCH,
    GRETL_TEST_RE_WALD,
    GRETL_TEST_XDEPEND,
    GRETL_TEST_MAX
} ModelTestType;

enum {
    VCV_SIMPLE,
    VCV_ROBUST,
    VCV_XPX
};

typedef enum {
    ARMA_X12A  = 1 << 0, /* using X-12-ARIMA (or X-13) to generate estimates */
    ARMA_EXACT = 1 << 1, /* using exact ML */
    ARMA_LS    = 1 << 2, /* using conditional ML, and O/NLS == CML */
    ARMA_OLS   = 1 << 3  /* OLS == MLE */
} ArmaFlags;

typedef enum {
    PANEL_HAC,  /* clustered by individual (Arellano) */
    PANEL_BK,   /* Beck-Katz PCSE */
    PANEL_TIME, /* clustered by period */
    PANEL_DK,   /* Driscoll-Kraay SCC */
    PANEL_BOTH  /* clustered by both unit and period */
} PanelVCVType;

typedef enum {
    RQ_ASY,   /* asymptotic */
    RQ_NID    /* sandwich */
} RQVCVType;

typedef enum {
    HAC_PREWHITEN = 1
} VCVFlags;

enum {
    SESSION_CLEAR_ALL,
    SESSION_CLEAR_DATASET
};

enum {
    BLAS_UNKNOWN,
    BLAS_NETLIB,
    BLAS_ATLAS,
    BLAS_OPENBLAS,
    BLAS_MKL,
    BLAS_VECLIB,
    BLAS_BLIS
};

enum gretl_warning_codes {
    W_GRADIENT = 1,
    W_GENMISS,     /* 2 */
    W_GENNAN,      /* 3 */
    W_MAX          /* 4 */
};

typedef enum {
    COLNAMES = 1 << 0,
    ROWNAMES = 1 << 1,
    REVERSED = 1 << 2
} NameFlags;

enum test_stats {
    GRETL_STAT_NONE,
    GRETL_STAT_NORMAL_CHISQ,
    GRETL_STAT_LM,
    GRETL_STAT_F,
    GRETL_STAT_LMF,
    GRETL_STAT_HARVEY_COLLIER,
    GRETL_STAT_RESET,
    GRETL_STAT_LR,
    GRETL_STAT_WALD_CHISQ,
    GRETL_STAT_SUP_WALD,
    GRETL_STAT_Z,
    GRETL_STAT_STUDENT,
    GRETL_STAT_LB_CHISQ,
    GRETL_STAT_WF
};

typedef enum {
    GRETL_MOD_NONE = 0,
    GRETL_MOD_TRANSPOSE,
    GRETL_MOD_SQUARE,
    GRETL_MOD_CUMULATE,
    GRETL_MOD_DECREMENT,
    GRETL_MOD_CTRANSP
} GretlMatrixMod;

typedef enum {
    CONF_NONE = 0,
    CONF_ELEMENTS,
    CONF_A_COLVEC,
    CONF_B_COLVEC,
    CONF_A_ROWVEC,
    CONF_B_ROWVEC,
    CONF_A_SCALAR,
    CONF_B_SCALAR,
    CONF_AC_BR,
    CONF_AR_BC
} ConfType;

typedef enum {
    V_SUM,
    V_PROD,
    V_MEAN
} GretlVecStat;

typedef enum {
    GRETL_MATRIX_SQUARE = 1,
    GRETL_MATRIX_LOWER_TRIANGULAR,
    GRETL_MATRIX_UPPER_TRIANGULAR,
    GRETL_MATRIX_SYMMETRIC,
    GRETL_MATRIX_DIAGONAL,
    GRETL_MATRIX_IDENTITY,
    GRETL_MATRIX_SCALAR,
} GretlMatrixStructure;

typedef enum {
    OP_EQ  = '=',
    OP_GT  = '>',
    OP_LT  = '<',
    OP_NEQ = 21,
    OP_GTE = 22,
    OP_LTE = 23
} GretlOp;

typedef enum {
    NO_MARKERS = 0,
    REGULAR_MARKERS,
    DAILY_DATE_STRINGS
} DatasetMarkerType;

typedef enum {
    VAR_DISCRETE   = 1 << 0,
    VAR_HIDDEN     = 1 << 1,
    VAR_GENERATED  = 1 << 2,
    VAR_LISTARG    = 1 << 3,
    VAR_TIMECOL    = 1 << 4,
    VAR_HFANCHOR   = 1 << 5,
    VAR_CODED      = 1 << 6
} VarFlags;

typedef enum {
    GRETL_FORMAT_TXT       = 1 << 0,
    GRETL_FORMAT_TEX       = 1 << 1,
    GRETL_FORMAT_DOC       = 1 << 2,
    GRETL_FORMAT_RTF       = 1 << 3,
    GRETL_FORMAT_RTF_TXT   = 1 << 4,
    GRETL_FORMAT_EQN       = 1 << 5,
    GRETL_FORMAT_SELECTION = 1 << 6,
    GRETL_FORMAT_CSV       = 1 << 7,
    GRETL_FORMAT_TAB       = 1 << 8,
    GRETL_FORMAT_MODELTAB  = 1 << 9,
    GRETL_FORMAT_LANDSCAPE = 1 << 10,
    GRETL_FORMAT_HAS_MINUS = 1 << 11,
    GRETL_FORMAT_XML       = 1 << 12
} PrnFormat;

typedef enum {
    GRETL_TYPE_NONE,
    GRETL_TYPE_BOOL,
    GRETL_TYPE_INT,
    GRETL_TYPE_UNSIGNED,
    GRETL_TYPE_OBS,
    GRETL_TYPE_LIST,
    GRETL_TYPE_DOUBLE,
    GRETL_TYPE_INT_ARRAY,
    GRETL_TYPE_DOUBLE_ARRAY,
    GRETL_TYPE_STRING,
    GRETL_TYPE_CMPLX_ARRAY,
    GRETL_TYPE_SERIES,
    GRETL_TYPE_MATRIX,
    GRETL_TYPE_STRUCT,
    GRETL_TYPE_SCALAR_REF,
    GRETL_TYPE_SERIES_REF,
    GRETL_TYPE_MATRIX_REF,
    GRETL_TYPE_STRING_REF,
    GRETL_TYPE_LIST_REF,
    GRETL_TYPE_USERIES,
    GRETL_TYPE_DATE,
    GRETL_TYPE_BUNDLE,
    GRETL_TYPE_BUNDLE_REF,
    GRETL_TYPE_ARRAY,
    GRETL_TYPE_ARRAY_REF,
    GRETL_TYPE_STRINGS,
    GRETL_TYPE_MATRICES,
    GRETL_TYPE_BUNDLES,
    GRETL_TYPE_LISTS,
    GRETL_TYPE_ARRAYS,
    GRETL_TYPE_STRINGS_REF,
    GRETL_TYPE_MATRICES_REF,
    GRETL_TYPE_BUNDLES_REF,
    GRETL_TYPE_LISTS_REF,
    GRETL_TYPE_ARRAYS_REF,
    GRETL_TYPE_VOID,
    GRETL_TYPE_NUMERIC,
    GRETL_TYPE_ANY
} GretlType;

typedef enum {
    USE_CWD         = 1 << 0,  /* store: use current dir as default */
    ECHO_ON         = 1 << 1,  /* echoing commands or not */
    MSGS_ON         = 1 << 2,  /* emitting non-error messages or not */
    FORCE_DECPOINT  = 1 << 3,  /* override locale decimal character */
    USE_SVD         = 1 << 4,  /* SVD decomposition is matrix OLS default */
    PREWHITEN       = 1 << 5,  /* HAC pre-whitening? */
    FORCE_HC        = 1 << 6,  /* don't use HAC for time series */
    USE_LBFGS       = 1 << 7,  /* prefer LBFGS to BFGS? */
    SHELL_OK        = 1 << 8,  /* "shell" facility is approved? */
    WARNINGS        = 1 << 9,  /* print numerical warning messages */
    SKIP_MISSING    = 1 << 10, /* skip NAs when building matrix from series */
    BFGS_RSTEP      = 1 << 11, /* use Richardson in BFGS numerical gradient */
    ROBUST_Z        = 1 << 12, /* use z- not t-score with HCCM/HAC */
    MWRITE_G        = 1 << 13, /* use %g format with mwrite() */
    MPI_USE_SMT     = 1 << 14, /* MPI: use hyperthreads by default */
    STATE_FLAG_MAX  = 1 << 15, /* separator */
    /* state small int (but non-boolean) vars */
    GRETL_OPTIM,
    VECM_NORM,
    GARCH_VCV,
    GARCH_ALT_VCV,
    ARMA_VCV,
    WILDBOOT_DIST,
    FDJAC_QUAL,
    USE_QR,
    MAX_VERBOSE,
    HC_VERSION,
    PANEL_ROBUST,
    HAC_KERNEL,
    HAC_LAG,
    USER_HAC_LAG,
    LBFGS_MEM,
    QUANTILE_TYPE,
    STATE_SMALL_INT_MAX, /* separator: start state int vars */
    HORIZON,
    BOOTREP,
    LOOP_MAXITER,
    BFGS_MAXITER,
    BFGS_VERBSKIP,
    BOOT_ITERS,
    BHHH_MAXITER,
    RQ_MAXITER,
    GMM_MAXITER,
    STATE_INT_MAX, /* separator: end state int vars */
    CONV_HUGE,
    NLS_TOLER,
    BFGS_TOLER,
    BFGS_MAXGRAD,
    BHHH_TOLER,
    QS_BANDWIDTH,
    NADARWAT_TRIM,
    STATE_FLOAT_MAX, /* separator: end state floats */
    CSV_WRITE_NA,
    CSV_READ_NA,
    INITVALS,
    INITCURV,
    MATMASK,
    STATE_VARS_MAX, /* separator */
    /* non-state vars follow */
    GRETL_DEBUG,
    GRETL_ASSERT,
    DATACOLS,
    PLOT_COLLECT,
    R_FUNCTIONS,
    R_LIB,
    LOGLEVEL,
    LOGSTAMP,
    CSV_DIGITS,
    HAC_MISSVALS,
    NS_SMALL_INT_MAX, /* separator */
    GMP_BITS,
    NS_MAX, /* separator */
    BLAS_MNK_MIN,
    OMP_MNK_MIN,
    OMP_N_THREADS,
    SIMD_K_MAX,
    SIMD_MN_MIN,
    USE_DCMT,
    NS_INT_MAX, /* separator */
    SEED,
    CSV_DELIM,
    STOPWATCH,
    VERBOSE,
    SV_WORKDIR,
    SV_LOGFILE,
    GRAPH_THEME,
    DISP_DIGITS,
    SETVAR_MAX /* sentinel */
} SetKey;

typedef enum {
    E_DATA = 2,
    E_SINGULAR,     /* 3 */
    E_DF,           /* 4 */
    E_ZERO,         /* 5 */
    E_TSS,          /* 6 */
    E_ESS,          /* 7 */
    E_NOTIMP,       /* 8 */
    E_UNSPEC,       /* 9 */
    E_PDWRONG,     /* 10 */
    E_FOPEN,       /* 11 */
    E_ALLOC,       /* 12 */
    E_EQN,         /* 13 */
    E_UNKVAR,      /* 14 */
    E_ARGS,        /* 15 */
    E_OLSONLY,     /* 16 */
    E_INVARG,      /* 17 */
    E_PARSE,       /* 18 */
    E_NOVARS,      /* 19 */
    E_NOOMIT,      /* 20 */
    E_NOADD,       /* 21 */
    E_ADDDUP,      /* 22 */
    E_LOGS,        /* 23 */
    E_SQUARES,     /* 24 */
    E_LAGS,        /* 25 */
    E_SQRT,        /* 26 */
    E_HIGH,        /* 27 */
    E_OBS,         /* 28 */
    E_NOCONST,     /* 29 */
    E_BADSTAT,     /* 30 */
    E_NOMERGE,     /* 31 */
    E_NOCONV,      /* 32 */
    E_CANCEL,      /* 33 */
    E_MISSDATA,    /* 34 */
    E_NAN,         /* 35 */
    E_NONCONF,     /* 36 */
    E_TYPES,       /* 37 */
    E_BADOPT,      /* 38 */
    E_NOIDENT,     /* 39 */
    E_EXTERNAL,    /* 40 */
    E_TOOLONG,     /* 41 */
    E_NODATA,      /* 42 */
    E_NOTPD,       /* 43 */
    E_JACOBIAN,    /* 44 */
    E_TOOFEW,      /* 45 */
    E_FNEST,       /* 46 */
    E_FUNCERR,     /* 47 : error set by function writer */
    E_STOP,        /* 48 : user aborted execution */
    E_BADCATCH,    /* 49 : "catch" used where it's not valid */
    E_CMPLX,       /* 50 : complex arguments/operands not supported */
    E_MIXED,       /* 51 : mixed complex/real operands not supported */
    E_DEPENDS,     /* 52 : gfn dependencies not met */
    E_DB_DUP,      /* 53 : duplicate vars found when saving to database */
    E_OK,          /* 54 : not really an error */
    E_MAX          /* 55 */
} GretlError; 

/**
 * bundled_item:
 *
 * An item of data within a gretl_bundle. This is an
 * opaque type; use the relevant accessor functions.
 */

typedef struct bundled_item_ {
    GretlType type;
    int size;
    gpointer data;
    char *note;
    char *name; /* pointer to associated key */
} bundled_item;

/**
 * gretl_bundle:
 *
 * An opaque type; use the relevant accessor functions.
 */

typedef struct gretl_bundle_ {
    BundleType type; /* see enum in gretl_bundle.h */
    GHashTable *ht;  /* holds key/value pairs */
    char *creator;   /* name of function that built the bundle */
    void *data;      /* holds pointer to struct for some uses */
} gretl_bundle;

typedef struct distmap_ {
	DistCode code;
	const char * s;
} distmap;

struct stored_opt_ {
    int ci;       /* index of the associated command */
    gretlopt opt; /* option flag */
    char *val;    /* option parameter value, or NULL */
    int flags;    /* may include OPT_SETOPT, OPT_PERSIST */
    int fd;       /* "function depth" at which registered */
};

typedef struct user_var_ {
    GretlType type;
    int level;
    UVFlags flags;
    char name[VNAMELEN];
    void *ptr;
} user_var;

typedef enum {
    D_NONE = 0,
    D_UNIFORM,
    D_UDISCRT,
    D_NORMAL,
    D_STUDENT,
    D_CHISQ,
    D_SNEDECOR,
    D_BINOMIAL,
    D_POISSON,
    D_EXPON,
    D_WEIBULL,
    D_GAMMA,
    D_GED,
    D_LAPLACE,
    D_BETA,
    D_DW,
    D_BINORM,
    D_JOHANSEN,
    D_BETABIN,
    D_NC_CHISQ,
    D_NC_F,
    D_NC_T,
    D_LOGISTIC,
    D_DIRICHLET
} DistCode;

enum {
    DROP_NORMAL,
    DROP_SPECIAL
};

struct gretl_option {
    int ci;              /* command index (gives context) */
    gretlopt o;          /* index of integer type */
    const char *longopt; /* "--"-style string representation of option */
    char parminfo;       /* 0 = option can never take a parameter,
                            1 = option may take a parameter,
                            2 = option requires a parameter
                         */
};
