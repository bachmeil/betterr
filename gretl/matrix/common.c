#include <common.h>
__import gretltypes;
__import globals;
__import glib;

int libset_get_int (SetKey key)
{
    void *valp;

    if (check_for_state()) {
	return 0;
    }

    valp = setkey_get_target(key, SV_INT);

    if (valp != NULL) {
#if SVDEBUG
	fprintf(stderr, "libset_get_int: valp %p\n", valp);
#endif
	if (libset_small_int(key)) {
	    return *(gint8 *) valp;
	} else {
	    return *(int *) valp;
	}
    //~ } else if (key == BLAS_MNK_MIN) {
	//~ return get_blas_mnk_min();
    //~ } else if (key == OMP_N_THREADS) {
	//~ return get_omp_n_threads();
    //~ } else if (key == OMP_MNK_MIN) {
	//~ return get_omp_mnk_min();
    //~ } else if (key == SIMD_K_MAX) {
	//~ return get_simd_k_max();
    //~ } else if (key == SIMD_MN_MIN) {
	//~ return get_simd_mn_min();
    //~ } else {
	//~ fprintf(stderr, "libset_get_int: unrecognized "
		//~ "key %d\n", key);
	return 0;
    }
}

/* check_for_state() returns non-zero if the program options
   state is not readable */

int check_for_state (void)
{
    //~ if (state == NULL) {
	//~ return libset_init();
    //~ } else {
//~ #if PDEBUG > 1
	//~ fprintf(stderr, "check_for_state: state = %p\n", (void *) state);
//~ #endif
	return 0;
    //~ }
}

void *setkey_get_target (SetKey key, SVType t)
{
    //~ int i = INTS_OFFSET + key - GRETL_OPTIM;
    //~ setvar *sv = &setvars[i];

    //~ if (sv->key != key) {
	//~ fprintf(stderr, "*** internal error, looking for %s, found %s ***\n",
		//~ setkey_get_name(key), sv->name);
	//~ return NULL;
    //~ } else if ((t == SV_INT && !libset_int(key)) ||
	       //~ (t == SV_DOUBLE && !libset_double(key))) {
	//~ fprintf(stderr, "*** type mismatch in setkey_get_target for %s ***\n",
		//~ sv->name);
	//~ return NULL;
    //~ } else {
	//~ return setvar_get_target(sv);
    //~ }
}



/**
 * series_adjust_sample:
 * @x: series to be checked for missing values.
 * @t1: on entry, initial start of sample range; on exit,
 *      start of sample range adjusted for missing values.
 * @t2: on entry, initial end of sample range; on exit, end
 *      of sample range adjusted for missing values.
 *
 * Adjusts @t1 and @t2 so as to drop any leading or trailing
 * missing observations.
 *
 * Returns: E_MISSDATA if interior missing values were found
 * within the (possibly adjusted) sample range, otherwise 0.
 */

int series_adjust_sample (const double *x, int *t1, int *t2)
{
    int t, t1min = *t1, t2max = *t2;
    int err = 0;

    for (t=t1min; t<t2max; t++) {
	if (na(x[t])) t1min++;
	else break;
    }

    for (t=t2max; t>t1min; t--) {
	if (na(x[t])) t2max--;
	else break;
    }

    for (t=t1min; t<=t2max; t++) {
	if (na(x[t])) {
	    err = E_MISSDATA;
	    break;
	}
    }

    *t1 = t1min;
    *t2 = t2max;

    return err;
}

/**
 * gretl_isconst:
 * @x: data series to examine.
 * @t1: starting observation.
 * @t2: ending observation.
 *
 * Check whether series @x is constant over the
 * given sample range (aside from any missing values).
 *
 * Returns: 1 if the variable is constant, otherwise 0.
 */

int gretl_isconst (int t1, int t2, const double *x)
{
    int t, ret = 1;

    while (na(x[t1]) && t1 <= t2) {
        t1++;
    }

    if (t1 >= t2) {
        return 0;
    }

    for (t=t1+1; t<=t2; t++) {
        if (na(x[t])) {
            continue;
        }
        if (floatneq(x[t], x[t1])) {
            ret = 0;
            break;
        }
    }

    return ret;
}

/**
 * gretl_errmsg_set:
 * @str: an error message.
 *
 * If %gretl_errmsg is currently blank, copy the given string into
 * the message space; or if the error message is not blank but
 * sufficient space remains, append @str to the message.
 */

void gretl_errmsg_set (const char *str)
{
#if EDEBUG
    fprintf(stderr, "gretl_errmsg_set: '%s'\n", str);
#endif

    if (alarm_set && *gretl_errmsg != '\0') {
	/* leave the current error message in place */
	return;
    }

    if (*gretl_errmsg == '\0') {
	strncat(gretl_errmsg, str, ERRLEN - 1);
    } else if (strcmp(gretl_errmsg, str)) {
	/* should we do the following? */
	int n = strlen(gretl_errmsg);
	int m = strlen(str);

	if (n + m + 2 < ERRLEN) {
	    strcat(gretl_errmsg, "\n");
	    strcat(gretl_errmsg, str);
	}
    }

#if EDEBUG
    fprintf(stderr, "gretl_errmsg now: '%s'\n", gretl_errmsg);
#endif
}

/**
 * gretl_compare_doubles:
 * @a: pointer to first element to compare.
 * @b: pointer to second element to compare.
 *
 * Comparison function for use with qsort.  Sorts doubles in
 * ascending order.
 *
 * Returns: appropriate value for qsort.
 */

int gretl_compare_doubles (const void *a, const void *b)
{
    const double *da = (const double *) a;
    const double *db = (const double *) b;

    if (isnan(*da) || isnan(*db)) {
        if (!isnan(*da)) {
            return -1;
        } else if (!isnan(*db)) {
            return 1;
        } else {
            return 0;
        }
    } else {
        return (*da > *db) - (*da < *db);
    }
}

/**
 * strings_array_free:
 * @strs: array of allocated strings.
 * @nstrs: number of strings in array.
 *
 * Frees each allocated string in @strs, then frees @strs itself.
 * Checks that @strs is not NULL before proceeding.
 */

void strings_array_free (char **strs, int nstrs)
{
    int i;

    if (strs != NULL) {
	for (i=0; i<nstrs; i++) {
	    free(strs[i]);
	}
	free(strs);
    }
}

void clear_gretl_matrix_err (void)
{
    gretl_matrix_err = 0;
}

/**
 * pputc:
 * @prn: gretl printing struct.
 * @c: character to print
 *
 * Returns: the number of bytes printed, or -1 on memory allocation
 * failure.
 */

int pputc (PRN *prn, int c)
{
    if (prn == NULL || prn->fixed) {
	return 0;
    }

    if (c == '\n') prn->nlcount += 1;
    else if (c != '\0') prn->nlcount = 0;

    if (prn->fp != NULL) {
	fputc(c, prn->fp);
	return 1;
    } 
    //~ else if (prn->fz != NULL) {
	//~ char s[2] = {c, '\0'};

	//~ return gzputs(prn->fz, s);
    //~ }

    if (prn->buf == NULL) {
	return 0;
    }

    if (prn->bufsize - prn->blen < MINREM) {
	//~ if (realloc_prn_buffer(prn, 0)) {
	    //~ return -1;
	//~ }
    }

#if PRN_DEBUG > 1
    fprintf(stderr, "pputc: adding char at %p\n", (void *) (prn->buf + prn->blen));
#endif

    prn->buf[prn->blen] = c;
    prn->buf[prn->blen + 1] = '\0';
    prn->blen += 1;

    return 1;
}

/**
 * pprintf:
 * @prn: gretl printing struct.
 * @format: as in the C library's printf().
 * @Varargs: arguments to be printed.
 *
 * Multi-purpose printing function: can output to stream, to buffer
 * or to nowhere (silently discarding the output), depending on
 * how @prn was initialized.  Note that it's preferable to use
 * pputs() for large chunks of fixed text.
 *
 * Returns: the number of bytes printed, or -1 on memory allocation
 * failure.
 */

int pprintf (PRN *prn, const char *format, ...)
{
    va_list args;
    int rem, plen = 0;

    if (prn == NULL || prn->fixed) {
	return 0;
    }

    if (format != NULL && *format != '\0') {
	//~ prn->nlcount = get_nl_count(format);
    }

    if (prn->fp != NULL) {
	/* printing to stream: straightforward */
	va_start(args, format);
	plen = vfprintf(prn->fp, format, args);
	va_end(args);
	return plen;
    } 
  //~ else if (prn->fz != NULL) {
	//~ gchar *tmp = NULL;

	//~ va_start(args, format);
	//~ plen = g_vasprintf(&tmp, format, args);
	//~ va_end(args);
	//~ gzputs(prn->fz, tmp);
	//~ g_free(tmp);
	//~ return plen;
    //~ }

    if (strncmp(format, "@init", 5) == 0) {
	//~ return pprintf_init(prn);
    }

    if (prn->buf == NULL) {
	return 0;
    }

    rem = prn->bufsize - prn->blen;
    if (rem < MINREM) {
	//~ if (realloc_prn_buffer(prn, 0)) {
	    //~ return -1;
	//~ }
    }

#if PRN_DEBUG > 1
    fprintf(stderr, "printing at %p\n", (void *) (prn->buf + prn->blen));
#endif

    /* printing to buffer: be careful not to overrun */
    rem = prn->bufsize - prn->blen - 1;
    va_start(args, format);
    plen = vsnprintf(prn->buf + prn->blen, rem, format, args);
    va_end(args);

    if (plen >= rem) {
	/* buffer not big enough: try again */
	size_t newsize = prn->bufsize + plen + 1024;

	//~ if (realloc_prn_buffer(prn, newsize)) {
	    //~ return -1;
	//~ }
	rem = prn->bufsize - prn->blen - 1;
	va_start(args, format);
	plen = vsnprintf(prn->buf + prn->blen, rem, format, args);
	va_end(args);
    }

    if (plen > 0) {
	prn->blen += plen;
    }

    return plen;
}

/**
 * gretl_matrix_print_to_prn:
 * @m: matrix to print.
 * @msg: accompanying message text (or NULL if no message is wanted).
 * @prn: pointer to gretl printing struct.
 *
 * Prints the matrix @m to @prn.
 */

void
gretl_matrix_print_to_prn (const gretl_matrix *m,
			   const char *msg, PRN *prn)
{
    //~ real_matrix_print_to_prn(m, msg, 0, NULL, NULL,
			     //~ -1, -1, prn);
}

/**
 * pputs:
 * @prn: gretl printing struct.
 * @s: constant string to print.
 *
 * Returns: the number of bytes printed, or -1 on memory allocation
 * failure.
 */

int pputs (PRN *prn, const char *s)
{
    int slen, rem;

    if (prn == NULL || prn->fixed) {
	return 0;
    }

    slen = strlen(s);
    if (slen > 0) {
	//~ prn->nlcount = get_nl_count(s);
    }

    if (prn->fp != NULL) {
	fputs(s, prn->fp);
	return slen;
    } 
  //~ else if (prn->fz != NULL) {
	//~ return gzputs(prn->fz, s);
    //~ }

    if (prn->buf == NULL) {
	return 0;
    }

    rem = prn->bufsize - prn->blen;

    while (rem < MINREM || rem <= slen) {
	//~ if (realloc_prn_buffer(prn, 0)) {
	    //~ return -1;
	//~ }
	rem = prn->bufsize - prn->blen;
    }

#if PRN_DEBUG > 1
    fprintf(stderr, "pputs: bufsize=%d, blen=%d, rem=%d, copying %d bytes\n",
	    (int) prn->bufsize, (int) prn->blen, rem, slen);
#endif

    strcpy(prn->buf + prn->blen, s);
    prn->blen += slen;

    return slen;
}

/**
 * gretl_errmsg_sprintf:
 * @fmt: format string.
 * @...: arguments, as to sprintf.
 *
 * Append a formatted message to the current gretl
 * error message.
 */

void gretl_errmsg_sprintf (const char *fmt, ...)
{
#if EDEBUG
    fprintf(stderr, "gretl_errmsg_sprintf: fmt='%s'\n", fmt);
#endif

    if (*gretl_errmsg == '\0') {
	va_list ap;

	va_start(ap, fmt);
	vsnprintf(gretl_errmsg, ERRLEN, fmt, ap);
	va_end(ap);
    } else if (strstr(gretl_errmsg, "*** error in fun") &&
	       strstr(fmt, "*** error in fun")) {
	/* don't print more than one "error in function" 
	   message, as this gets confusing 
	*/
	;
    } else {
	/* find the number of characters left */
	int len0 = strlen(gretl_errmsg);
	int n = ERRLEN - len0 - 2;

	if (n > 31) {
	    char tmp[ERRLEN];
	    va_list ap;

	    *tmp = '\0';
	    va_start(ap, fmt);
	    vsnprintf(tmp, n, fmt, ap);
	    va_end(ap);

	    if (gretl_errmsg[len0 - 1] != '\n') {
		strcat(gretl_errmsg, "\n");
	    }
	    strcat(gretl_errmsg, tmp);
	} 
    }
}

/**
 * vge_gamma_test:
 * @x: data series.
 * @t1: start of sample range.
 * @t2: end of sample range.
 * @err: location to receive error code.
 *
 * Performs the test described by J. A. Villaseñor and E. González-Estrada
 * (Statistics and Probability Letters, 96 (2015) pp. 281–286) for the null
 * hypothesis that @x is gamma-distributed over the range  @t1 to @t2.
 *
 * For the sake of compatibility with the gamma_test() function in the R
 * package named "goft" we divide by n-1 in computing the covariance
 * term @sxz, although this is not recommended by the authors of the above-
 * noted publication.
 *
 * Returns: the z value for the test, or #NADBL on error.
 */

//~ double vge_gamma_test (const double *x, int t1, int t2, int *err)
//~ {
    //~ int t, n = t2 - t1 + 1;
    //~ double xbar = 0, zbar = 0;
    //~ double xc, zc, sxz, s2;
    //~ double a, V, Vstar;

    //~ for (t=t1; t<=t2; t++) {
        //~ if (x[t] <= 0) {
            //~ gretl_errmsg_set(_("Non-positive values encountered"));
            //~ *err = E_DATA;
            //~ return NADBL;
        //~ } else if (na(x[t])) {
	    //~ n--;
	//~ } else {
            //~ xbar += x[t];
            //~ zbar += log(x[t]);
        //~ }
    //~ }

    //~ if (n < 30) {
        //~ /* minimum obs? */
        //~ *err = E_TOOFEW;
        //~ return NADBL;
    //~ }

    //~ xbar /= n;
    //~ zbar /= n;

    //~ sxz = s2 = 0;
    //~ for (t=t1; t<=t2; t++) {
        //~ if (!na(x[t])) {
            //~ xc = x[t] - xbar;
            //~ zc = log(x[t]) - zbar;
            //~ sxz += xc * zc;
            //~ s2 += xc * xc;
        //~ }
    //~ }

    //~ sxz /= (n-1);
    //~ s2 /= (n-1);
    //~ V = s2 / (xbar * sxz);
    //~ a = xbar / sxz;
    //~ Vstar = sqrt(n*a) * (V-1);

    //~ return fabs(Vstar) / sqrt(2.0);
//~ }



int dist_code_from_string (const char *s)
{
  distmap dmap[] = {
	{ D_UNIFORM,  "u" },
	{ D_UDISCRT,  "i" },
	{ D_NORMAL,   "z" },
	{ D_STUDENT,  "t" },
	{ D_CHISQ,    "x" },
	{ D_SNEDECOR, "f" },
	{ D_BINOMIAL, "b" },
	{ D_POISSON,  "p" },
	{ D_EXPON,    "exp" },
	{ D_WEIBULL,  "w" },
	{ D_GAMMA,    "g" },
	{ D_GED,      "e" },
	{ D_LAPLACE,  "l" },
	{ D_BETA,     "beta" },
	{ D_DW,       "d" },
	{ D_BINORM,   "D" },
	{ D_JOHANSEN, "J" },
	{ D_BETABIN,  "bb" },
	{ D_NC_CHISQ, "ncx" },
	{ D_NC_F,     "ncf" },
	{ D_NC_T,     "nct" },
	{ D_LOGISTIC, "s" },
	{ D_DIRICHLET, "dir" },
	{ D_NONE,     NULL }
    };
    char test[8];
    int i;

    if (!strcmp(s, "D")) {
	/* special: case counts for bivariate normal */
	return D_BINORM;
    }

    /* otherwise we'll ignore case */
    for (i=0; i<8 && s[i]; i++) {
	test[i] = tolower(s[i]);
    }
    test[i] = '\0';

    for (i=0; dmap[i].code; i++) {
	if (!strcmp(test, dmap[i].s)) {
	    return dmap[i].code;
	}
    }

    /* backward compatibility */
    if (!strcmp(test, "n")) {
	return D_NORMAL;
    } else if (!strcmp(test, "c")) {
	return D_CHISQ;
    } else if (!strcmp(test, "lgt")) {
	return D_LOGISTIC;
    }

    return D_NONE;
}

/**
 * series_get_string_table:
 * @dset: pointer to dataset.
 * @i: index number of series.
 *
 * Returns: the string table attched to series @i or NULL if
 * there is no such table.
 */

series_table *series_get_string_table (const DATASET *dset, int i)
{
    if (dset != NULL && i > 0 && i < dset->v) {
	return dset->varinfo[i]->st;
    } else {
	return NULL;
    }
}

/**
 * series_table_get_strings:
 * @st: a gretl series table.
 * @n_strs: location to receive the number of strings, or NULL.
 *
 * Returns: the array of strings associated with @st. These
 * should not be modified in any way.
 */

char **series_table_get_strings (series_table *st, int *n_strs)
{
    if (st != NULL) {
	if (n_strs != NULL) {
	    *n_strs = st->n_strs;
	}
	return st->strs;
    } else {
	return NULL;
    }
}

/**
 * count_distinct_values:
 * @x: sorted array of doubles.
 * @n: number of elements in array.
 *
 * Returns: the number of distinct values in array @x,
 * provided that @x is already sorted.
 */

int count_distinct_values (const double *x, int n)
{
    int i, c = 1;

    for (i=1; i<n; i++) {
        if (x[i] != x[i-1]) {
            c++;
        }
    }

    return c;
}

/**
 * is_string_valued:
 * @dset: pointer to dataset.
 * @i: index number of series.
 *
 * Returns: 1 if series @i has a table of string values
 * (that is, a mapping from numerical values to associated
 * string values), otherwise 0.
 */

int is_string_valued (const DATASET *dset, int i)
{
    if (dset != NULL && i > 0 && i < dset->v) {
	return dset->varinfo[i]->st != NULL;
    } else {
	return 0;
    }
}

/**
 * series_is_discrete:
 * @dset: pointer to dataset.
 * @i: index number of series.
 *
 * Returns: non-zero iff series @i should be treated as discrete.
 */

//~ int series_is_discrete (const DATASET *dset, int i)
//~ {
    //~ return dset->varinfo[i]->flags & VAR_DISCRETE;
//~ }

int few_vals (int t1, int t2, const double *x, double *ratio)
{
    double test[FEW_VALS];
    int match;
    int i, t, n = 0, nv = 0;

    for (t=t1; t<=t2; t++) {
        if (!na(x[t])) {
            match = 0;
            for (i=0; i<nv; i++) {
                if (x[t] == test[i]) {
                    /* we've already seen this value */
                    match = 1;
                    break;
                }
            }
            if (!match) {
                /* x[t] is a "new" value */
                if (nv == FEW_VALS) {
                    /* hit the limit */
                    nv++;
                    break;
                }
                test[nv++] = x[t];
            }
            n++;
        }
    }

    /* ratio of distinct values to observations */
    *ratio = nv / (double) n;

    return nv;
}

/**
 * gretl_isdiscrete:
 * @x: data series to examine.
 * @t1: starting observation.
 * @t2: ending observation.
 *
 * Checks the variable @x over the range @t1 to @t2 for discreteness.
 * This is a heuristic whose components are (a) whether the values
 * are "fairly round" (multiples of 0.25) or not, and, if test (a) is
 * passed, (b) whether the variable takes on only "few" distinct
 * values.
 *
 * Returns: 0 if test (a) is not passed or the number of distinct values
 * is > 32; else 1 if the number of distinct values is <= 32; else 2 if
 * the number of distinct values is <= 8.  A return of 1 is supposed
 * to indicate that it's "reasonable" to treat @x as discrete, while
 * a return of 2 indicates that it's probably unreasonable _not_ to
 * treat @x as discrete for the purpose of drawing up a frequency
 * distribution.
 */

int gretl_isdiscrete (int t1, int t2, const double *x)
{
    int t, n = 0, disc = 1;
    int allints = 1;
    double r = 0;

    for (t=t1; t<=t2; t++) {
        if (na(x[t])) {
            continue;
        }
        n++;
        if (!ok_int(x[t])) {
            allints = disc = 0;
            break;
        }
        r = x[t] - floor(x[t]);
        if (allints && r != 0) {
            allints = 0;
        }
        if (r != 0.0 && r != 0.25 && r != 0.5 && r != 0.75) {
            disc = 0;
            break;
        }
    }

    if (n == 0) {
        disc = 0;
    }

    if (disc) {
        n = few_vals(t1, t2, x, &r);
        if (allints) {
            if (n <= FEW_VALS && r < 0.2) {
                disc = 2;
            } else {
                disc = (n <= FEWER_VALS)? 2 : 1;
            }
        } else if (n > FEW_VALS) {
            disc = 0;
        } else if (r > 0.9 && n > 30) {
            /* somewhat arbitrary: but if r (= ratio of distinct
               values to number of cases) is "too high", and the
               number of observations is not tiny, perhaps we should
               not automatically take the var as discrete
            */
            disc = 0;
        } else if (n <= FEWER_VALS) {
            disc = 2;
        }
    }

    return disc;
}

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

// Shouldn't be needed
int get_optval_int (int ci, gretlopt opt, int *err)
{
    //~ stored_opt *so = matching_stored_opt(ci, opt);
    //~ int status = option_parm_status(ci, opt);
    int ret = 0;

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

    return ret;
}

// Shouldn't be needed
stored_opt *matching_stored_opt (int ci, gretlopt opt)
{
    //~ int i, fd = gretl_function_depth();

//~ #if OPTDEBUG
    //~ fprintf(stderr, "matching_stored_opt? ci=%d, fd=%d, opt=%d, n_stored=%d\n",
            //~ ci, fd, opt, n_stored_opts);
//~ #endif

    //~ for (i=0; i<n_stored_opts; i++) {
	//~ stored_opt *so = &optinfo[i];

        //~ if (so->ci == ci && so->opt == opt && so->fd == fd) {
            //~ return so;
        //~ }
    //~ }

    return NULL;
}

/* for a given (@ci, @opt) pair, determine its status with regard
   to a parameter value: 0 = not allowed, 1 = optional, 2 = required
*/

int option_parm_status (int ci, gretlopt opt)
{
    //~ int i, got_ci = 0;

    //~ if (opt == OPT_U) {
	//~ if (plot_output_opt_ok(ci)) {
	    //~ ci = GNUPLOT;
	//~ } else if (cmd_plot_opt_ok(ci)) {
	    //~ ci = CORR;
	//~ }
    //~ } else if (opt == OPT_b && plot_outbuf_opt_ok(ci)) {
        //~ ci = GNUPLOT;
    //~ }

    //~ for (i=0; gretl_opts[i].ci != 0; i++) {
        //~ if (gretl_opts[i].ci == ci) {
            //~ if (gretl_opts[i].o == opt) {
                //~ return gretl_opts[i].parminfo;
            //~ }
            //~ got_ci = 1;
        //~ } else if (got_ci) {
            //~ break;
        //~ }
    //~ }

    return 0;
}

/**
 * gretl_int_from_string:
 * @s: string to examine.
 * @err: location to receive error code.
 *
 * If @s is a valid string representation of an integer,
 * return that integer, otherwise if @s is the name of a
 * scalar variable, return the value of that variable,
 * provided it can be converted to an integer, otherwise
 * set the content of @err to a non-zero value.
 *
 * Returns: integer value.
 */

int gretl_int_from_string (const char *s, int *err)
{
    char *test;
    int n = 0;

    if (s == NULL || *s == '\0') {
        *err = E_DATA;
        return 0;
    }

    if (isalpha(*s)) {
        double x = get_scalar_value_by_name(s, err);

        if (!*err) {
            n = gretl_int_from_double(x, err);
        }
        return n;
    }

    errno = 0;
    n = strtol(s, &test, 10);

    if (*test != '\0' || errno == ERANGE) {
        *err = E_DATA;
        errno = 0;
        return 0;
    }

    return n;
}



/* more "permissive" than gretl_scalar_get_value(): allows
   for @name being the identifier for a 1 x 1 matrix, or
   bundle.member
*/

double get_scalar_value_by_name (const char *name, int *err)
{
    double ret = NADBL;
    //~ user_var *u;

    //~ if (strchr(name, '.')) {
        //~ ret = maybe_get_bundled_scalar(name, err);
        //~ goto bailout;
    //~ }

    //~ u = get_user_var_by_name(name);

    //~ if (u != NULL) {
        //~ if (u->type == GRETL_TYPE_DOUBLE) {
            //~ ret = *(double *) u->ptr;
        //~ } else if (u->type == GRETL_TYPE_MATRIX) {
            //~ gretl_matrix *m = u->ptr;

            //~ if (gretl_matrix_is_scalar(m)) {
                //~ ret = m->val[0];
            //~ } else {
                //~ *err = E_TYPES;
            //~ }
        //~ } else {
            //~ *err = E_TYPES;
        //~ }
    //~ } else {
        //~ ret = get_const_by_name(name, err);
    //~ }

 //~ bailout:

    //~ if (*err) {
        //~ gretl_errmsg_sprintf(_("'%s': not a scalar"), name);
    //~ }

    return ret;
}

void gretl_push_c_numeric_locale (void)
{
    return;
}

void gretl_pop_c_numeric_locale (void)
{
    return;
}

// This shouldn't be needed, so I didn't import all the associated crap
double generate_scalar (const char *s, DATASET *dset, int *err) {return NADBL;}

/**
 * gretl_double_from_string:
 * @s: string to examine.
 * @err: location to receive error code.
 *
 * If @s is a valid string representation of a double,
 * return its value, otherwise if @s is the name of a
 * scalar variable, return the value of that variable,
 * otherwise set the content of @err to a non-zero value.
 *
 * Returns: double value.
 */

double gretl_double_from_string (const char *s, int *err)
{
    char *test;
    double x;

    if (s == NULL || *s == '\0') {
        *err = E_DATA;
        return NADBL;
    }

    if (isalpha(*s)) {
        return get_scalar_value_by_name(s, err);
    }

    gretl_push_c_numeric_locale();
    errno = 0;
    x = strtod(s, &test);
    gretl_pop_c_numeric_locale();

    if (*test != '\0' || errno == ERANGE) {
        *err = E_DATA;
        errno = 0;
        return NADBL;
    }

    return x;
}

/**
 * get_optval_double:
 * @ci: gretl command index.
 * @opt: gretl option value.
 * @err: location to receive error code.
 *
 * Returns: the double-precision ancillary value currently
 * associated with option @opt for command @ci, if any,
 * otherwise #NADBL. If @opt is an active option for
 * @ci but the parameter for this option cannot be
 * interpreted as a numerical value, E_INVARG is written
 * into @err.
 */

double get_optval_double (int ci, gretlopt opt, int *err)
{
    stored_opt *so = matching_stored_opt(ci, opt);
    double ret = NADBL;

    if (so != NULL && so->val != NULL) {
        ret = gretl_double_from_string(so->val, err);
        if (err) {
            ret = generate_scalar(so->val, NULL, err);
        }
        if (*err) {
            gretl_errmsg_sprintf(_("%s: invalid option argument"), so->val);
            *err = E_INVARG;
        }
    }

    return ret;
}

//~ /**
 //~ * gretl_matrix_set_colnames:
 //~ * @m: target matrix.
 //~ * @S: array of strings.
 //~ *
 //~ * Sets an array of strings on @m which can be retrieved
 //~ * using gretl_matrix_get_colnames(). Note that @S must
 //~ * contain as many strings as @m has columns. The matrix
 //~ * takes ownership of @S, which should be allocated and
 //~ * not subsequently touched by the caller.
 //~ *
 //~ * Returns: 0 on success, non-zero code on error.
 //~ */

//~ int gretl_matrix_set_colnames (gretl_matrix *m, char **S)
//~ {
    //~ if (m == NULL) {
        //~ return E_DATA;
    //~ } else if (is_block_matrix(m)) {
        //~ return matrix_block_error("gretl_matrix_set_colnames");
    //~ } else if (S != NULL && m->info == NULL &&
               //~ gretl_matrix_add_info(m)) {
        //~ return E_ALLOC;
    //~ }

    //~ if (m->info != NULL) {
        //~ if (m->info->colnames != NULL) {
            //~ strings_array_free(m->info->colnames, m->cols);
        //~ }
        //~ m->info->colnames = S;
    //~ }

    //~ return 0;
//~ }

/**
 * gretl_matrix_set_rownames:
 * @m: target matrix.
 * @S: array of strings.
 *
 * Sets an array of strings on @m which can be retrieved
 * using gretl_matrix_get_rownames(). Note that @S must
 * contain as many strings as @m has rows. The matrix
 * takes ownership of @S, which should be allocated and
 * not subsequently touched by the caller.
 *
 * Returns: 0 on success, non-zero code on error.
 */

//~ int gretl_matrix_set_rownames (gretl_matrix *m, char **S)
//~ {
    //~ if (m == NULL) {
        //~ return E_DATA;
    //~ } else if (is_block_matrix(m)) {
        //~ return matrix_block_error("gretl_matrix_set_rownames");
    //~ } else if (S != NULL && m->info == NULL &&
               //~ gretl_matrix_add_info(m)) {
        //~ return E_ALLOC;
    //~ }

    //~ if (m->info != NULL) {
        //~ if (m->info->rownames != NULL) {
            //~ strings_array_free(m->info->rownames, m->rows);
        //~ }
        //~ m->info->rownames = S;
    //~ }

    //~ return 0;
//~ }

//~ DATASET *create_auxiliary_dataset (int nvar, int nobs, gretlopt opt)
//~ {
    //~ DATASET *dset = real_create_new_dataset(nvar, nobs, opt);

    //~ if (dset != NULL) {
	//~ if (opt & OPT_B) {
	    //~ dset->auxiliary = Z_COLS_BORROWED;
	//~ } else {
	    //~ dset->auxiliary = 1;
	//~ }
    //~ }

    //~ return dset;
//~ }

/**
 * destroy_dataset:
 * @dset: pointer to dataset.
 *
 * Frees all resources associated with @dset.
 */

//~ void destroy_dataset (DATASET *dset)
//~ {
    //~ if (dset != NULL) {
	//~ free_Z(dset);
	//~ clear_datainfo(dset, CLEAR_FULL);
	//~ free(dset);
    //~ }
//~ }

/**
 * gretl_matrix_get_colnames:
 * @m: matrix
 *
 * Returns: The array of strings set on @m using
 * gretl_matrix_set_colnames(), or NULL if no such
 * strings have been set. The returned array will
 * contain as many strings as @m has columns.
 */

//~ const char **gretl_matrix_get_colnames (const gretl_matrix *m)
//~ {
    //~ if (has_colnames(m)) {
        //~ return (const char **) m->info->colnames;
    //~ } else {
        //~ return NULL;
    //~ }
//~ }

/**
 * gretl_matrix_get_rownames:
 * @m: matrix
 *
 * Returns:The array of strings set on @m using
 * gretl_matrix_set_rownames(), or NULL if no such
 * strings have been set. The returned array will
 * contain as many strings as @m has rows.
 */

//~ const char **gretl_matrix_get_rownames (const gretl_matrix *m)
//~ {
    //~ if (has_rownames(m)) {
        //~ return (const char **) m->info->rownames;
    //~ } else {
        //~ return NULL;
    //~ }
//~ }

/**
 * gretl_matrix_xtab:
 * @x: data array.
 * @y: data array.
 * @n: length of the two arrays.
 * @err: location to receive error code.
 *
 * Computes the cross tabulation of the values contained in the
 * arrays @x (by row) and @y (by column). These should generally
 * be discrete values otherwise the cross-tabulation may be
 * very large and uninformative.
 *
 * Returns: the generated matrix, or NULL on failure.
 */

//~ gretl_matrix *gretl_matrix_xtab (const double *x,
                                 //~ const double *y,
				 //~ int n, int *err)
//~ {
    //~ gretl_matrix *tab = NULL;
    //~ gretl_matrix *vx = NULL;
    //~ gretl_matrix *vy = NULL;
    //~ int complete;
    //~ int i, j, n_ok;

    //~ *err = 0;

    //~ n_ok = ok_xy_count(x, y, n, err);
    //~ if (*err) {
        //~ return NULL;
    //~ }

    //~ complete = (n_ok == n);

    //~ if (complete) {
	//~ /* no need to dodge missing values */
	//~ vx = gretl_matrix_values(x, n, OPT_S, err);
	//~ if (!*err) {
	    //~ vy = gretl_matrix_values(y, n, OPT_S, err);
	//~ }
    //~ } else {
	//~ double *tmp = malloc(2 * n_ok * sizeof *tmp);
	//~ double *okx, *oky;

	//~ if (tmp == NULL) {
	    //~ *err = E_ALLOC;
	//~ } else {
	    //~ okx = tmp;
	    //~ oky = tmp + n_ok;
	    //~ for (i=0, j=0; i<n; i++) {
		//~ if (complete_obs(x, y, i)) {
		    //~ okx[j] = x[i];
		    //~ oky[j] = y[i];
		    //~ j++;
		//~ }
	    //~ }
	    //~ vx = gretl_matrix_values(okx, n_ok, OPT_S, err);
	    //~ if (!*err) {
		//~ vy = gretl_matrix_values(oky, n_ok, OPT_S, err);
	    //~ }
	    //~ free(tmp);
	//~ }
    //~ }

    //~ if (!*err) {
	//~ tab = gretl_zero_matrix_new(vx->rows, vy->rows);
	//~ if (tab == NULL) {
	    //~ *err = E_ALLOC;
	//~ }
    //~ }

    //~ if (!*err) {
	//~ double **X = doubles_array_new(n_ok, 2);

	//~ if (X == NULL) {
	    //~ *err = E_ALLOC;
	//~ } else {
	    //~ for (i=0, j=0; i<n; i++) {
		//~ if (complete || complete_obs(x, y, i)) {
		    //~ X[j][0] = x[i];
		    //~ X[j][1] = y[i];
		    //~ j++;
		//~ }
	    //~ }
        //~ }
	//~ make_matrix_xtab(X, n_ok, vx, vy, tab);
	//~ doubles_array_free(X, n_ok);
    //~ }

    //~ gretl_matrix_free(vx);
    //~ gretl_matrix_free(vy);

    //~ return tab;
//~ }

/**
 * series_table_get_string:
 * @st: a gretl series table.
 * @val: the numerical value to look up.
 *
 * Returns: the string associated with @val in the
 * given series table, or NULL in case there is no match.
 */

//~ const char *series_table_get_string (series_table *st, double val)
//~ {
    //~ const char *ret = NULL;

    //~ if (!na(val)) {
	//~ int k = (int) lrint(val);

	//~ if (k > 0 && k <= st->n_strs) {
	    //~ ret = st->strs[k-1];
	//~ }
    //~ }

    //~ return ret;
//~ }

//~ static gretl_matrix *gretl_filled_matrix_new (int r, int c,
                                              //~ double val)
//~ {
    //~ gretl_matrix *m = NULL;

    //~ if (r < 0 || c < 0) {
        //~ return NULL;
    //~ } else if (r == 0 || c == 0) {
        //~ m = gretl_null_matrix_new();
        //~ if (m != NULL) {
            //~ m->rows = r;
            //~ m->cols = c;
        //~ }
    //~ } else {
        //~ int i, n = r * c;

        //~ m = gretl_matrix_alloc(r, c);
        //~ if (m != NULL) {
            //~ if (val == 0.0) {
                //~ memset(m->val, 0, n * sizeof *m->val);
            //~ } else {
                //~ for (i=0; i<n; i++) {
                    //~ m->val[i] = val;
                //~ }
            //~ }
        //~ }
    //~ }

    //~ return m;
//~ }

/**
 * gretl_zero_matrix_new:
 * @r: desired number of rows in the matrix.
 * @c: desired number of columns in the matrix.
 *
 * Returns: pointer to a newly allocated zero matrix, or NULL
 * on failure.
 */

//~ gretl_matrix *gretl_zero_matrix_new (int r, int c)
//~ {
    //~ return gretl_filled_matrix_new(r, c, 0.0);
//~ }

//~ static gretl_matrix *
//~ real_gretl_matrix_values (const double *x, int n,
			  //~ gretlopt opt, int *n_vals,
			  //~ int *missvals, int *err)
//~ {
    //~ gretl_matrix *v = NULL;
    //~ double *sorted = NULL;
    //~ double last;
    //~ int m = 0;
    //~ int i, k;

    //~ sorted = malloc(n * sizeof *sorted);
    //~ if (sorted == NULL) {
        //~ *err = E_ALLOC;
        //~ return NULL;
    //~ }

    //~ k = 0;
    //~ for (i=0; i<n; i++) {
        //~ if (!na(x[i])) {
            //~ sorted[k++] = x[i];
        //~ }
    //~ }

    //~ if (k == 0) {
	//~ if (n_vals == NULL) {
	    //~ v = gretl_null_matrix_new();
	//~ }
	//~ goto bailout;
    //~ } else if (missvals != NULL) {
	//~ *missvals = n - k;
    //~ }

    //~ qsort(sorted, k, sizeof *sorted, gretl_compare_doubles);
    //~ m = count_distinct_values(sorted, k);

    //~ if (n_vals != NULL) {
	//~ /* the caller just wants the count */
	//~ *n_vals = m;
	//~ goto bailout;
    //~ }

    //~ v = gretl_column_vector_alloc(m);
    //~ if (v == NULL) {
        //~ *err = E_ALLOC;
        //~ goto bailout;
    //~ }

    //~ if (opt & OPT_S) {
        //~ /* sorted */
        //~ v->val[0] = last = sorted[0];
        //~ for (i=1, m=1; i<k; i++) {
            //~ if (sorted[i] != last) {
                //~ last = sorted[i];
                //~ v->val[m++] = sorted[i];
            //~ }
        //~ }
    //~ } else {
        //~ /* unsorted */
        //~ int j, add;

        //~ for (i=0, m=0; i<n; i++) {
            //~ if (!na(x[i])) {
                //~ add = 1;
                //~ for (j=0; j<m; j++) {
                    //~ if (v->val[j] == x[i]) {
                        //~ add = 0;
                        //~ break;
                    //~ }
                //~ }
                //~ if (add) {
                    //~ v->val[m++] = x[i];
                //~ }
            //~ }
        //~ }
    //~ }

 //~ bailout:

    //~ free(sorted);

    //~ return v;
//~ }

/**
 * gretl_matrix_values:
 * @x: array to process.
 * @n: length of array.
 * @opt: if OPT_S the array of values will be sorted, otherwise
 * given in order of occurrence.
 * @err: location to receive error code.
 *
 * Returns: an allocated matrix containing the distinct
 * values in array @x, skipping any missing values, or
 * NULL on failure.
 */

//~ gretl_matrix *gretl_matrix_values (const double *x, int n,
                                   //~ gretlopt opt, int *err)
//~ {
    //~ return real_gretl_matrix_values(x, n, opt, NULL, NULL, err);
//~ }

/* retrieve a matrix result directly */

//~ gretl_matrix *generate_matrix (const char *s, DATASET *dset,
			       //~ int *err)
//~ {
    //~ gretl_matrix *m = NULL;
    //~ parser p;

    //~ *err = realgen(s, &p, dset, NULL, P_PRIV | P_ANON, MAT);

    //~ if (!*err) {
	//~ NODE *n = p.ret;

	//~ if (n->t == MAT) {
	    //~ if (n->flags & TMP_NODE) {
		//~ /* steal the generated matrix */
		//~ m = n->v.m;
		//~ n->v.m = NULL;
	    //~ } else {
		//~ m = gretl_matrix_copy(n->v.m);
		//~ if (m == NULL) {
		    //~ *err = E_ALLOC;
		//~ }
	    //~ }
	//~ } else if (n->t == NUM) {
	    //~ if (na(n->v.xval)) {
		//~ *err = E_NAN;
	    //~ } else {
		//~ m = gretl_matrix_alloc(1, 1);
		//~ if (m == NULL) {
		    //~ *err = E_ALLOC;
		//~ } else {
		    //~ m->val[0] = n->v.xval;
		//~ }
	    //~ }
	//~ } else {
	    //~ *err = E_TYPES;
	//~ }
    //~ } else if (*err == 1) {
	//~ *err = E_PARSE;
    //~ }

    //~ gen_cleanup(&p);

    //~ return m;
//~ }

/**
 * gretl_list_new:
 * @nterms: the maximum number of elements to be stored in the list.
 *
 * Creates a newly allocated list with space for @nterms elements,
 * besides the leading element, which in a gretl list always
 * holds a count of the number of elements that follow.  This
 * leading element is initialized appropriately.  For example, if
 * @nterms = 4, space for 5 integers is allocated and the first
 * element of the array is set to 4.  The other elements of
 * the list are initialized to 0.
 *
 * Returns: the newly allocated list, or NULL on failure.
 */

int *gretl_list_new (int nterms)
{
    int *list = NULL;
    int i;

    if (nterms < 0) {
	return NULL;
    }

    list = malloc((nterms + 1) * sizeof *list);

    if (list != NULL) {
	list[0] = nterms;
	for (i=1; i<=nterms; i++) {
	    list[i] = 0;
	}
    }

    return list;
}

/**
 * gretl_list_separator_position:
 * @list: an array of integer variable ID numbers, the first element
 * of which holds a count of the number of elements following.
 *
 * Returns: if @list contains the separator for compound
 * lists, #LISTSEP, the position in @list at which this is found,
 * else 0.  The search begins at position 1.
 */

int gretl_list_separator_position (const int *list)
{
    int i;

    if (list != NULL) {
	for (i=1; i<=list[0]; i++) {
	    if (list[i] == LISTSEP) {
		return i;
	    }
	}
    }

    return 0;
}

/**
 * get_matrix_by_name:
 * @name: name of the matrix.
 *
 * Looks up a user-defined matrix by name.
 *
 * Returns: pointer to matrix, or %NULL if not found.
 */

//~ gretl_matrix *get_matrix_by_name (const char *name)
//~ {
    //~ gretl_matrix *ret = NULL;

    //~ if (name != NULL && *name != '\0') {
	//~ user_var *u;

	//~ u = get_user_var_of_type_by_name(name, GRETL_TYPE_MATRIX);
	//~ if (u != NULL) {
	    //~ ret = user_var_get_value(u);
	//~ }
    //~ }

    //~ return ret;
//~ }


/**
 * gretl_bundle_copy:
 * @bundle: gretl bundle to be copied.
 * @err: location to receive error code.
 *
 * Returns: a "deep copy" of @bundle (all the items in @bundle
 * are themselves copied), or NULL on failure.
 */

gretl_bundle *gretl_bundle_copy (const gretl_bundle *bundle, int *err)
{
    gretl_bundle *bcpy = NULL;

    //~ if (bundle == NULL) {
        //~ *err = E_DATA;
    //~ } else {
        //~ if (bundle->type == BUNDLE_KALMAN) {
            //~ bcpy = kalman_bundle_copy(bundle, err);
        //~ } else {
            //~ bcpy = gretl_bundle_new();
            //~ if (bcpy == NULL) {
                //~ *err = E_ALLOC;
            //~ }
        //~ }
        //~ if (!*err) {
            //~ g_hash_table_foreach(bundle->ht, copy_bundled_item, bcpy);
        //~ }
    //~ }

    return bcpy;
}

/**
 * gretl_bundle_destroy:
 * @bundle: bundle to destroy.
 *
 * Frees all contents of @bundle as well as the pointer itself.
 */

void gretl_bundle_destroy (gretl_bundle *bundle)
{
    //~ if (bundle != NULL) {
        //~ if (bundle->ht != NULL) {
            //~ g_hash_table_destroy(bundle->ht);
        //~ }
        //~ free(bundle->creator);
        //~ if (bundle->type == BUNDLE_KALMAN) {
            //~ kalman_free(bundle->data);
        //~ }
        //~ free(bundle);
    //~ }
}


/**
 * gretl_compare_strings:
 * @a: pointer to first element to compare.
 * @b: pointer to second element to compare.
 *
 * Comparison function for use with qsort.  Sorts strings in
 * alphabetical order.
 *
 * Returns: appropriate value for qsort.
 */

int gretl_compare_strings (const void *a, const void *b)
{
    const char **sa = (const char **) a;
    const char **sb = (const char **) b;

    return g_strcmp0(*sa, *sb);
}    

/**
 * strings_array_sort:
 * @pS: location of array of strings.
 * @n: location of the number of strings in the array.
 * @opt: may contain %OPT_U to trim the sorted array
 * so that it contains only unique entries.
 *
 * Sorts an array of strings in ascending lexicographical
 * order. If %OPT_U is given, @n holds the number of unique
 * strings on exit. It is assumed that storage for the
 * strings array was obtained via strings_array_new() or
 * a similar libgretl function.
 *
 * Returns: 0 on success, non-zero on error.
 */

int strings_array_sort (char ***pS, int *n, gretlopt opt)
{
    char **S;
    int ns;

    if (pS == NULL || n == NULL) {
	return E_DATA;
    }

    S = *pS;
    ns = *n;

    qsort(S, ns, sizeof *S, gretl_compare_strings);

    if (opt & OPT_U) {
	int i, j, m = ns;

	for (i=0; i<m-1; i++) {
	    if (!strcmp(S[i], S[i+1])) {
		free(S[i+1]);
		for (j=i+1; j<m-1; j++) {
		    S[j] = S[j+1];
		}
		S[m-1] = NULL;
		i--;
		m--;
	    }
	}
	if (m < ns) {
	    char **tmp = realloc(S, m * sizeof *S);

	    if (tmp != NULL) {
		*pS = tmp;
	    }
	    *n = m;
	}
    }

    return 0;
}

/**
 * gretl_strdup_printf:
 * @format: as in printf().
 * @Varargs: arguments to be printed.
 *
 * Print the arguments according to @format.
 *
 * Returns: allocated result of the printing, or NULL on failure.
 */

char *gretl_strdup_printf (const char *format, ...)
{
    va_list args;
    char *buf = NULL;
    int len;

#ifdef HAVE_VASPRINTF
    va_start(args, format);
    len = vasprintf(&buf, format, args);
    va_end(args);
    if (len < 0) {
	buf = NULL;
    }
#else
    int bsize = 2048;

    buf = malloc(bsize);
    if (buf == NULL) {
	return NULL;
    }

    memset(buf, 0, 1);

    va_start(args, format);
    len = vsnprintf(buf, bsize, format, args);
    va_end(args);

    if (len >= bsize) {
	fputs("gretl_strdup_printf warning: string was truncated\n",
	      stderr);
    }
#endif

    return buf;
}

/**
 * get_optval_string:
 * @ci: gretl command index.
 * @opt: gretl option value.
 *
 * Returns: the ancillary string value currently
 * associated with option @opt for command @ci, if any,
 * otherwise %NULL.
 */

const char *get_optval_string (int ci, gretlopt opt)
{
    stored_opt *so = matching_stored_opt(ci, opt);

    return (so != NULL)? so->val : NULL;
}

/**
 * print_xtab:
 * @tab: gretl cross-tabulation struct.
 * @dset: pointer to dataset, or NULL.
 * @opt: may contain %OPT_R to print row percentages, %OPT_C
 * to print column percentages, %OPT_Z to display zero entries,
 * %OPT_T to print as TeX (LaTeX), %OPT_N to omit marginal
 * totals, %OPT_B to bold-face counts where the row and
 * column values are equal (TeX only), %OPT_S to record
 * Pearson test result.
 * @prn: gretl printing struct.
 *
 * Print crosstab to @prn.
 */

//~ void print_xtab (const Xtab *tab, const DATASET *dset,
		 //~ gretlopt opt, PRN *prn)
//~ {
    //~ int done = 0;

    //~ if (opt & OPT_T) {
	//~ /* --tex : are we printing to file? */
	//~ const char *fname = get_optval_string(XTAB, OPT_T);

	//~ if (fname != NULL && *fname != '\0') {
	    //~ int err = 0;
	    //~ PRN *xprn;

	    //~ gretl_maybe_switch_dir(fname);
	    //~ xprn = gretl_print_new_with_filename(fname, &err);
	    //~ if (err) {
		//~ pprintf(prn, _("Couldn't write to %s"), fname);
		//~ pputc(prn, '\n');
	    //~ } else {
		//~ real_print_xtab(tab, dset, opt, xprn);
		//~ gretl_print_destroy(xprn);
	    //~ }
	    //~ done = 1;
	//~ }
    //~ }

    //~ if (!done) {
	//~ real_print_xtab(tab, dset, opt, prn);
    //~ }
//~ }

//~ static void real_print_xtab (const Xtab *tab, const DATASET *dset,
			     //~ gretlopt opt, PRN *prn)
//~ {
    //~ double x, y, cj, ri;
    //~ int n5 = 0;
    //~ double ymin = 1.0e-7;
    //~ double pearson = 0.0;
    //~ double pval = NADBL;
    //~ char lbl[64];
    //~ int dashlen = 0;
    //~ int rlen = 0;
    //~ int clen = 0;
    //~ int cw, tw;
    //~ int stdw = 6;
    //~ int totals = 1;
    //~ int tex = 0;
    //~ int bold = 0;
    //~ int i, j;

    //~ if (opt & OPT_T) {
	//~ /* LaTeX output */
	//~ tex = 1;
	//~ if ((opt & OPT_E) && !(tab->rstrs || tab->cstrs)) {
	    //~ /* bold-face equal values */
	    //~ bold = 1;
	//~ }
	//~ /* don't bother with chi-square */
	//~ pearson = NADBL;
    //~ }

    //~ if (opt & OPT_N) {
	//~ totals = 0;
    //~ }

    //~ if (*tab->rvarname != '\0' && *tab->cvarname != '\0') {
	//~ pputc(prn, '\n');
	//~ if (tex) {
	    //~ tex_xtab_heading(tab, prn);
	//~ } else {
	    //~ pprintf(prn, _("Cross-tabulation of %s (rows) against %s (columns)"),
		    //~ tab->rvarname, tab->cvarname);
	    //~ pputs(prn, "\n\n");
	//~ }
    //~ } else {
	//~ pputc(prn, '\n');
    //~ }

    //~ if (tex) {
	//~ int nc = tab->cols + 1 + totals;

	//~ pputs(prn, "\\begin{tabular}{");
	//~ pputs(prn, "r|");
	//~ for (j=1; j<nc; j++) {
	    //~ pputc(prn, 'r');
	//~ }
	//~ pputs(prn, "}\n");
    //~ }

    //~ if (tex) {
	//~ pputs(prn, "     & ");
    //~ }

    //~ if (tab->Sr != NULL) {
	//~ rlen = 2 + row_strlen(tab);
	//~ if (rlen > 16) {
	    //~ rlen = 16;
	//~ } else if (rlen < 7) {
	    //~ rlen = 7;
	//~ }
    //~ } else {
	//~ rlen = 7;
    //~ }

    //~ /* allow for BIG integers */
    //~ get_xtab_col_widths(tab, stdw, &cw, &tw);

    //~ if (tab->Sc != NULL) {
	//~ clen = 2 + col_strlen(tab);
	//~ if (clen > 12) {
	    //~ clen = 12;
	//~ } else if (clen < cw) {
	    //~ clen = 6;
	//~ }
    //~ } else {
	//~ clen = cw;
    //~ }

    //~ if (totals) {
        //~ dashlen = clen + 1 + tab->cols*cw + tw;
    //~ }

    //~ bufspace(rlen, prn);

    //~ /* header row: column labels */

    //~ for (j=0; j<tab->cols; j++) {
	//~ if (tab->Sc != NULL) {
	    //~ *lbl = '\0';
	    //~ strncat(lbl, tab->Sc[j], 63);
	    //~ gretl_utf8_truncate(lbl, clen-2);
	    //~ if (tex) {
		//~ pputs(prn, lbl);
	    //~ } else {
		//~ bufspace(clen - g_utf8_strlen(lbl, -1), prn);
		//~ pputs(prn, lbl);
	    //~ }
	//~ } else {
	    //~ cj = tab->cval[j];
	    //~ if (tex) {
		//~ pprintf(prn, "%4g", cj);
	    //~ } else {
		//~ pprintf(prn, "[%*g]", cw-2, cj);
	    //~ }
	//~ }
	//~ if (tex) {
	    //~ if (j < tab->cols - 1 || totals) {
		//~ pputs(prn, " & ");
	    //~ } else {
		//~ pputs(prn, "\\\\ \\hline\n");
	    //~ }
	//~ }
    //~ }

    //~ if (totals) {
	//~ if (tex) {
	    //~ pprintf(prn,"$\\Sigma$\\\\ \\hline\n");
	//~ } else {
	    //~ bufspace(2 + (tw - stdw), prn);
	    //~ pprintf(prn,"%s\n", _("TOT."));
            //~ dashline(dashlen, prn);
	//~ }
    //~ } else if (!tex) {
	//~ pputc(prn, '\n');
    //~ }

    //~ /* body of table */

    //~ for (i=0; i<tab->rows; i++) {
	//~ if (tab->rtotal[i] == 0) {
	    //~ continue;
	//~ }
	//~ if (tab->Sr != NULL) {
	    //~ *lbl = '\0';
	    //~ strncat(lbl, tab->Sr[i], 63);
	    //~ gretl_utf8_truncate(lbl, rlen-2);
	    //~ pputs(prn, lbl);
	    //~ if (!tex) {
		//~ bufspace(rlen - g_utf8_strlen(lbl, -1), prn);
	    //~ }
	//~ } else {
	    //~ ri = tab->rval[i];
	    //~ if (tex) {
		//~ pprintf(prn, "%4g", ri);
	    //~ } else {
		//~ pprintf(prn, "[%4g] ", ri);
	    //~ }
	//~ }
	//~ if (tex) {
	    //~ pputs(prn, " & ");
	//~ }
	//~ /* row counts */
	//~ for (j=0; j<tab->cols; j++) {
	    //~ if (tab->ctotal[j] > 0) {
		//~ if (tab->f[i][j] || (opt & OPT_Z)) {
		    //~ if (opt & (OPT_C | OPT_R)) {
			//~ if (opt & OPT_C) {
			    //~ x = 100.0 * tab->f[i][j] / tab->ctotal[j];
			//~ } else {
			    //~ x = 100.0 * tab->f[i][j] / tab->rtotal[i];
			//~ }
			//~ if (tex) {
			    //~ /* FIXME strvals! */
			    //~ if (bold && tab->cval[j] == tab->rval[i]) {
				//~ pprintf(prn, "\\textbf{%.1f%%%%}", x);
			    //~ } else {
				//~ pprintf(prn, "%*.1f%%%%", clen, x);
			    //~ }
			//~ } else {
			    //~ pprintf(prn, "%*.1f%%", clen-1, x);
			//~ }
		    //~ } else {
			//~ if (bold && tab->cval[j] == tab->rval[i]) {
			    //~ pprintf(prn, "\\textbf{%d} ", tab->f[i][j]);
			//~ } else {
			    //~ pprintf(prn, "%*d", clen, tab->f[i][j]);
			//~ }
		    //~ }
		//~ } else if (!tex) {
		    //~ bufspace(clen, prn);
		//~ }
		//~ if (tex && (totals || j < tab->cols-1)) {
		    //~ pputs(prn, "& ");
		//~ }
	    //~ }
	    //~ if (!na(pearson)) {
		//~ /* cumulate chi-square */
		//~ y = ((double) tab->rtotal[i] * tab->ctotal[j]) / tab->n;
		//~ if (y < ymin) {
		    //~ pearson = NADBL;
		//~ } else {
		    //~ x = (double) tab->f[i][j] - y;
		    //~ pearson += x * x / y;
		    //~ if (y >= 5) {
			//~ n5++;
		    //~ }
		//~ }
	    //~ }
	//~ }
	//~ if (totals) {
	    //~ /* row totals */
	    //~ if (opt & OPT_C) {
		//~ x = 100.0 * tab->rtotal[i] / tab->n;
		//~ if (tex) {
		    //~ pprintf(prn, "%5.1f%%%%", x);
		//~ } else {
		    //~ pprintf(prn, "%5.1f%%", x);
		//~ }
	    //~ } else {
		//~ pprintf(prn, "%*d", tw, tab->rtotal[i]);
	    //~ }
	//~ }
	//~ /* terminate row */
	//~ if (tex) {
	    //~ if (totals && i == tab->rows-1) {
		//~ pputs(prn, "\\\\ [2pt]\n");
	    //~ } else {
		//~ pputs(prn, "\\\\\n");
	    //~ }
	//~ } else {
	    //~ pputc(prn, '\n');
	//~ }
    //~ }

    //~ /* footer row */

    //~ if (totals) {
	//~ /* column totals */
	//~ if (tex) {
	    //~ pputs(prn, "$\\Sigma$ & ");
	//~ } else {
            //~ if (dashlen > 0) {
                //~ dashline(dashlen, prn);
            //~ } else {
                //~ pputc(prn, '\n');
            //~ }
	    //~ pprintf(prn, "%-*s", rlen, _("TOTAL"));
	//~ }
	//~ for (j=0; j<tab->cols; j++) {
	    //~ if (opt & OPT_R) {
		//~ x = 100.0 * tab->ctotal[j] / tab->n;
		//~ if (tex) {
		    //~ pprintf(prn, "%*.1f%%%%", clen, x);
		//~ } else {
		    //~ pprintf(prn, "%*.1f%%", clen-1, x);
		//~ }
	    //~ } else {
		//~ pprintf(prn, "%*d", clen, tab->ctotal[j]);
	    //~ }
	    //~ if (tex) {
		//~ pputs(prn, "& ");
	    //~ }
	//~ }
	//~ pprintf(prn, "%*d\n", tw, tab->n);
    //~ }

    //~ if (tex) {
	//~ pputs(prn, "\\end{tabular}\n");
	//~ return;
    //~ }

    //~ /* additional information, if applicable */

    //~ if (tab->missing) {
	//~ pputc(prn, '\n');
	//~ pprintf(prn, _("%d missing values"), tab->missing);
	//~ pputc(prn, '\n');
    //~ }

    //~ if (na(pearson)) {
	//~ pputc(prn, '\n');
	//~ pprintf(prn, _("Pearson chi-square test not computed: some "
		       //~ "expected frequencies were less\n"
		       //~ "than %g\n"), ymin);
    //~ } else {
	//~ double n5p = (double) n5 / (tab->rows * tab->cols);
	//~ int df = (tab->rows - 1) * (tab->cols - 1);

	//~ pval = chisq_cdf_comp(df, pearson);
	//~ if (na(pval)) {
	    //~ pearson = NADBL;
	//~ } else {
	    //~ pputc(prn, '\n');
	    //~ pprintf(prn, _("Pearson chi-square test = %g (%d df, p-value = %g)"),
		    //~ pearson, df, pval);
	    //~ pputc(prn, '\n');
	    //~ if (!tex && n5p < 0.80) {
		//~ /* xgettext:no-c-format */
		//~ pputs(prn, _("Warning: Less than of 80% of cells had expected "
			     //~ "values of 5 or greater.\n"));
	    //~ }
	//~ }
    //~ }

    //~ if (opt & OPT_S) {
	//~ /* saving Pearson test result */
	//~ record_test_result(pearson, pval);
    //~ }

    //~ if (!(opt & OPT_F) && !tex && tab->rows == 2 && tab->cols == 2) {
	//~ fishers_exact_test(tab, prn);
    //~ }
//~ }

/**
 * graphyx:
 * @y: y-axis data.
 * @x: x-axis data.
 * @n: number of observations.
 * @yname: y-axis label.
 * @xname: x-axis label.
 * @prn: gretl printing struct.
 *
 * Generates a simple ascii scatter-plot of @y against @x and 
 * prints the plot to @prn.
 *
 * Returns: 0 on successful completion, error code on error.
 */

//~ int graphyx (const double *y, const double *x, int n,
	     //~ const char *yname, const char *xname, 
	     //~ PRN *prn)
//~ {
    //~ return graphyzx(y, NULL, x, n, yname, NULL, xname, 
		    //~ OPT_NONE, prn);
//~ }

/**
 * gretl_error_clear:
 *
 * Blank out any previously recorded error message.
 */

int gretl_error_clear (void)
{
#if EDEBUG
    fprintf(stderr, "gretl_error_clear\n");
#endif
    if (!alarm_set) {
	*gretl_errmsg = '\0';
    }
    error_printed = 0;
    errno = 0;
    
    return 0;
}

const char *series_get_graph_name (const DATASET *dset, int i)
{
    const char *ret = dset->varname[i];

    if (dset->varinfo != NULL && dset->varinfo[i] != NULL) {
	if (dset->varinfo[i]->display_name[0] != '\0') {
	    ret = dset->varinfo[i]->display_name;
	}
    }

    return ret;
}

void bufspace (int n, PRN *prn)
{
    while (n-- > 0) {
	pputc(prn, ' ');
    }
}

/**
 * list_adjust_sample:
 * @list: list of variables to be tested for missing values,
 * or %NULL to test all series.
 * @t1: on entry, initial start of sample range; on exit,
 *      start of sample range adjusted for missing values.
 * @t2: on entry, initial end of sample range; on exit, end
 *      of sample range adjusted for missing values.
 * @dset: dataset struct.
 * @nmiss: location to receive number of missing values within
 * (possibly adjusted) sample range.
 *
 * Drops leading or trailing observations from the sample range
 * initially given by the values in @t1 and @t2 if missing values
 * are found for any of the variables given in @list.
 *
 * If @nmiss is non-NULL it receives the number of missing values
 * inside the (possibly reduced) sample range, otherwise it is
 * considered an error if there are any such missing values.
 *
 * Returns: 0 on success or %E_MISSDATA or error.
 */

int list_adjust_sample (const int *list, int *t1, int *t2,
			const DATASET *dset, int *nmiss)
{
    int i, t, t1min = *t1, t2max = *t2;
    int k, vi, missing, err = 0;

    if (list != NULL) {
	k = list[0];
    } else {
	/* check all series */
	k = dset->v - 1;
    }

    /* advance start of sample range to skip missing obs? */
    for (t=t1min; t<t2max; t++) {
	missing = 0;
	for (i=1; i<=k; i++) {
	    vi = list == NULL ? i : list[i];
	    if (vi > 0 && vi != LISTSEP) {
		if (na(dset->Z[vi][t])) {
		    missing = 1;
		    break;
		}
	    }
	}
	if (missing) {
	    t1min++;
	} else {
	    break;
	}
    }

    /* retard end of sample range to skip missing obs? */
    for (t=t2max; t>t1min; t--) {
	missing = 0;
	for (i=1; i<=k; i++) {
	    vi = list == NULL ? i : list[i];
	    if (vi > 0 && vi != LISTSEP) {
		if (na(dset->Z[vi][t])) {
		    missing = 1;
		    break;
		}
	    }
	}
	if (missing) {
	    t2max--;
	} else {
	    break;
	}
    }

    if (nmiss != NULL) {
	*nmiss = 0;
    }

    /* check for missing values within remaining range */
    for (t=t1min; t<=t2max && !err; t++) {
	missing = 0;
	for (i=1; i<=k; i++) {
	    vi = list == NULL ? i : list[i];
	    if (vi > 0 && vi != LISTSEP) {
		if (na(dset->Z[vi][t])) {
		    if (nmiss == NULL) {
			err = E_MISSDATA;
		    } else {
			*nmiss += 1;
		    }
		    break;
		}
	    }
	}
    }

    *t1 = t1min;
    *t2 = t2max;

    return err;
}

/**
 * gretl_unit_matrix_new:
 * @r: desired number of rows in the matrix.
 * @c: desired number of columns in the matrix.
 *
 * Returns: pointer to a newly allocated matrix, all
 * of whose elements equal 1, or NULL on failure.
 */

//~ gretl_matrix *gretl_unit_matrix_new (int r, int c)
//~ {
    //~ return gretl_filled_matrix_new(r, c, 1.0);
//~ }

/**
 * gretl_matrix_raise:
 * @m: matrix to operate on.
 * @x: exponent.
 *
 * Raises each element of @m to the power @x.
 */

//~ void gretl_matrix_raise (gretl_matrix *m, double x)
//~ {
    //~ if (!gretl_is_null_matrix(m)) {
        //~ int i, n = m->rows * m->cols;

        //~ for (i=0; i<n; i++) {
            //~ m->val[i] = pow(m->val[i], x);
        //~ }
    //~ }
//~ }

/**
 * gretl_vector_mean:
 * @v: input vector.
 *
 * Returns: the arithmetic mean of the elements of @v, or
 * #NADBL on failure.
 */

double gretl_vector_mean (const gretl_vector *v)
{
    double num = 0.0;
    int i, n, den = 0;

    if (gretl_is_null_matrix(v)) {
        return NADBL;
    }

    n = gretl_vector_get_length(v);
    if (n == 0) {
        return NADBL;
    }

    for (i=0; i<n; i++) {
        if (!na(v->val[i])) {
            num += v->val[i];
            den++;
        }
    }

    return (den > 0)? (num / den) : NADBL;
}

/**
 * gretl_vector_from_series:
 * @x: series from data array.
 * @t1: starting observation.
 * @t2: ending observation.
 *
 * Returns: a newly allocated gretl_vector containing the values
 * of the given data series for the given range, or NULL on failure.
 * Any missing values in the input array are preserved as NaNs in
 * the output.
 */

//~ gretl_vector *gretl_vector_from_series (const double *x,
					//~ int t1, int t2)
//~ {
    //~ gretl_matrix *v = NULL;
    //~ int n = t2 - t1 + 1;

    //~ if (n > 0) {
	//~ size_t sz = n * sizeof *x;

	//~ v = gretl_column_vector_alloc(n);
	//~ if (v != NULL) {
	    //~ memcpy(v->val, x + t1, sz);
	    //~ gretl_matrix_set_t1(v, t1);
	    //~ gretl_matrix_set_t2(v, t2);
	//~ }
    //~ }

    //~ return v;
//~ }

/**
 * gretl_fix_exponent:
 * @s: string representation of floating-point number.
 *
 * Some C libraries (e.g. MS) print an "extra" zero in the exponent
 * when using scientific notation, e.g. "1.45E-002".  This function
 * checks for this and cuts it out if need be.
 *
 * Returns: the corrected numeric string.
 */

char *gretl_fix_exponent (char *s)
{
    char *p;
    int n;

    if ((p = strstr(s, "+00")) || (p = strstr(s, "-00"))) {
	if (*(p+3)) {
	    memmove(p+1, p+2, strlen(p+1));
	}
    }

    n = strlen(s);
    if (s[n-1] == '.' || s[n-1] == ',') {
	/* delete trailing junk */
	s[n-1] = '\0';
    }

    return s;
}

/**
 * gretl_matrix_init:
 * @m: matrix to be initialized.
 *
 * Initializes @m to be zero by zero with NULL data.
 */

gretl_matrix *gretl_matrix_init (gretl_matrix *m)
{
    m->rows = m->cols = 0;
    m->val = NULL;
    m->info = NULL;
    m->is_complex = 0;
    m->z = NULL;
    return m;
}

/**
 * ntolabel:
 * @datestr: char array to which date is to be printed.
 * @t: zero-based observation number.
 * @dset: data information struct.
 *
 * Prints to @datestr (which must be at least #OBSLEN bytes)
 * the calendar representation of observation number @t.
 *
 * Returns: the observation string.
 */

//~ char *ntolabel (char *datestr, int t, const DATASET *dset)
//~ {
    //~ double x;

//~ #if 0
    //~ fprintf(stderr, "ntolabel: t=%d, pd=%d, sd0=%g, incoming stobs='%s'\n",
	    //~ t, dset->pd, dset->sd0, dset->stobs);
    //~ fprintf(stderr, " calendar_data(dset) %d\n", calendar_data(dset));
//~ #endif

    //~ if (calendar_data(dset)) {
	//~ /* handles both daily and dated weekly data */
	//~ if (dataset_has_markers(dset)) {
	    //~ strcpy(datestr, dset->S[t]);
	    //~ if (strchr(datestr, '/')) {
		//~ gretl_charsub(datestr, '/', '-');
	    //~ }
	//~ } else {
	    //~ calendar_date_string(datestr, t, dset);
	//~ }
	//~ return datestr;
    //~ } else if (dataset_is_daily(dset) ||
	       //~ dataset_is_weekly(dset)) {
	//~ /* undated time series */
	//~ x = date_as_double(t, 1, dset->sd0);
	//~ sprintf(datestr, "%d", (int) x);
	//~ return datestr;
    //~ } else if (dataset_is_decennial(dset)) {
	//~ x = dset->sd0 + 10 * t;
	//~ sprintf(datestr, "%d", (int) x);
	//~ return datestr;
    //~ } else if (dataset_is_panel(dset)) {
	//~ panel_obs(datestr, t, dset);
	//~ return datestr;
    //~ }

    //~ x = date_as_double(t, dset->pd, dset->sd0);

    //~ if (dset->pd == 1) {
        //~ sprintf(datestr, "%d", (int) x);
    //~ } else {
	//~ int pdp = dset->pd;
	//~ short len = 1;
	//~ char fmt[10];

	//~ while ((pdp = pdp / 10)) len++;
	//~ sprintf(fmt, "%%.%df", len);
	//~ sprintf(datestr, fmt, x);
	//~ colonize_obs(datestr);
    //~ }

    //~ return datestr;
//~ }

char *maybe_trim_varname (char *targ, const char *src)
{
    int srclen = strlen(src);

    if (srclen < NAMETRUNC) {
	strcpy(targ, src);
    } else {
	const char *p = strrchr(src, '_');

	*targ = '\0';

	if (p != NULL && isdigit(*(p+1)) && strlen(p) < 4) {
	    /* preserve lag identifier? */
	    int snip = srclen - NAMETRUNC + 2;
	    int fore = p - src;

	    strncat(targ, src, fore - snip);
	    strcat(targ, "~");
	    strcat(targ, p);
	} else {
	    strncat(targ, src, NAMETRUNC - 2);
	    strcat(targ, "~");
	}
    }

    return targ;
}

int get_gretl_digits (void)
{
    return gretl_digits;
}

/**
 * gretl_list_copy:
 * @src: an array of integers, the first element of which holds
 * a count of the number of elements following.
 *
 * Returns: an allocated copy @src (or NULL if @src is NULL).
 */

int *gretl_list_copy (const int *src)
{
    int *targ = NULL;

    if (src != NULL) {
	int n = src[0] + 1;

	targ = malloc(n * sizeof *targ);
	if (targ != NULL) {
	    memcpy(targ, src, n * sizeof *targ);
	}
    }

    return targ;
}

/**
 * doubles_array_new:
 * @m: number of sub-arrays.
 * @n: length of each sub-array.
 *
 * Allocates a 2-dimensional array of doubles, that is,
 * @m arrays each containing @n elements.  If @n is
 * zero the sub-arrays are just set to %NULL.
 *
 * Returns: the allocated array, or %NULL on failure.
 */

double **doubles_array_new (int m, int n)
{
    double **X;
    int i;

    if (m == 0) {
        return NULL;
    }

    X = malloc(m * sizeof *X);
    if (X == NULL) {
        return X;
    }

    for (i=0; i<m; i++) {
        X[i] = NULL;
    }

    if (n > 0) {
        for (i=0; i<m; i++) {
            X[i] = malloc(n * sizeof **X);
            if (X[i] == NULL) {
                doubles_array_free(X, m);
                X = NULL;
                break;
            }
        }
    }

    return X;
}

/**
 * doubles_array_free:
 * @X: 2-dimensional array of doubles.
 * @m: number of sub-arrays.
 *
 * Frees a 2-dimensional array of doubles, first freeing
 * each sub-array.
 */

void doubles_array_free (double **X, int m)
{
    if (X != NULL) {
        int i;

        for (i=0; i<m; i++) {
            free(X[i]);
        }
        free(X);
    }
}

/**
 * gretl_list_delete_at_pos:
 * @list: an array of integers, the first element of which holds
 * a count of the number of elements following.
 * @pos: position at which to delete list element.
 *
 * Deletes the element at position @pos from @list and moves any
 * remaining elements forward.  Decrements the value of the first,
 * counter, element of @list.
 *
 * Returns: 0 on success, 1 on error.
 */

int gretl_list_delete_at_pos (int *list, int pos)
{
    int i, err = 0;

    if (pos < 1 || pos > list[0]) {
	err = 1;
    } else {
	for (i=pos; i<list[0]; i++) {
	    list[i] = list[i + 1];
	}

	list[list[0]] = 0;
	list[0] -= 1;
    }

    return err;
}

int count_selection (const char *s, int n)
{
    int i, c = 0;

    for (i=0; i<n; i++) {
	if (s[i] != 0) c++;
    }

    return c;
}

/**
 * gretl_int_from_double:
 * @x: double-precision floating point value
 * @err: location to receive error code.
 *
 * Returns: the value of @x converted to an integer, if
 * possible. Otherwise returns -1 with @err set to a
 * non-zero value. Note that it is considered an
 * error if @x is "too far" from the nearest integer;
 * it must be "almost integral", with tolerance 0.001.
 */

int gretl_int_from_double (double x, int *err)
{
    int k = -1;

    if (na(x) || fabs(x) > INT_MAX || fabs(x - nearbyint(x)) > 0.001) {
        *err = E_INVARG;
    } else {
        k = (int) lrint(x);
    }

    return k;
}

/**
 * ijton:
 * @i: row number (0-based)
 * @j: column number (0-based)
 * @nrows: number of rows (and columns) in symmetric matrix.
 *
 * Given a (row, column) reference into a symmetric 2-dimensional
 * matrix A, finds the index into a 1-dimensional array x
 * composed of the non-redundant (lower) elements of A.
 *
 * E.g. for the 3 x 3 case with 6 non-redundant elements, 0 to 5,
 *
 *    A(0,0) = x[0]  A(0,1) = x[1]  A(0,2) = x[2]
 *    A(1,0) = x[1]  A(1,1) = x[3]  A(1,2) = x[4]
 *    A(2,0) = x[2]  A(2,1) = x[4]  A(2,1) = x[5]
 *
 * Returns: 0-based index into flat array.
 */

int ijton (int i, int j, int nrows)
{
    if (i > j) {
        int tmp = i;

        i = j;
        j = tmp;
    }

    return nrows * i + j - i - ((i - 1) * i / 2);
}

/**
 * gretl_inverse_compare_doubles:
 * @a: pointer to first element to compare.
 * @b: pointer to second element to compare.
 *
 * Comparison function for use with qsort.  Sorts doubles in
 * descending order.
 *
 * Returns: appropriate value for qsort.
 */

int gretl_inverse_compare_doubles (const void *a, const void *b)
{
    const double *da = (const double *) a;
    const double *db = (const double *) b;

    if (isnan(*da) || isnan(*db)) {
        if (!isnan(*da)) {
            return 1;
        } else if (!isnan(*db)) {
            return -1;
        } else {
            return 0;
        }
    } else {
        return (*da < *db) - (*da > *db);
    }
}

/**
 * strings_array_new:
 * @nstrs: number of strings in array.
 *
 * Allocates storage for @nstrs strings and initializes all
 * to NULL.
 *
 * Returns: the allocated array, or NULL on failure.
 */

char **strings_array_new (int nstrs)
{
    char **s;
    int i;

    if (nstrs <= 0) {
	return NULL;
    }

    s = malloc(nstrs * sizeof *s);
    if (s != NULL) {
	for (i=0; i<nstrs; i++) {
	    s[i] = NULL;
	}
    }

    return s;
}

/**
 * gretl_strdup:
 * @src: the string to duplicate.
 *
 * Returns: an allocated copy of @src, or NULL on error.
 */

char *gretl_strdup (const char *src)
{
    char *targ = NULL;

    if (src != NULL) {
	size_t n = strlen(src) + 1;

	targ = calloc(n, 1);
	if (targ != NULL) {
	    memcpy(targ, src, n);
	}
    }

    return targ;
}

/**
 * strings_array_dup:
 * @strs: array of strings to be copied.
 * @n: number of strings in array.
 *
 * Returns: an allocated copy of @strs, or NULL on failure.
 */

char **strings_array_dup (char **strs, int n)
{
    char **S = NULL;
    int i, err = 0;

    if (n <= 0 || strs == NULL) {
	return NULL;
    }

    S = strings_array_new(n);
    if (S == NULL) return NULL;

    for (i=0; i<n; i++) {
	if (strs[i] == NULL) {
	    S[i] = NULL;
	} else {
	    S[i] = gretl_strdup(strs[i]);
	    if (S[i] == NULL) {
		err = E_ALLOC;
		break;
	    }
	}
    }

    if (err) {
	strings_array_free(S, n);
	S = NULL;
    }

    return S;
}

/**
 * series_table_get_string:
 * @st: a gretl series table.
 * @val: the numerical value to look up.
 *
 * Returns: the string associated with @val in the
 * given series table, or NULL in case there is no match.
 */

const char *series_table_get_string (series_table *st, double val)
{
    const char *ret = NULL;

    if (!na(val)) {
	int k = (int) lrint(val);

	if (k > 0 && k <= st->n_strs) {
	    ret = st->strs[k-1];
	}
    }

    return ret;
}


/**
 * print_time:
 * @s: string into which to print: must be at least 48 bytes.
 *
 * Returns: @s, which will contain a string representation of the
 * current date and time, in the format YYYY-mm-dd H:M.
 */

char *print_time (char *s)
{
    time_t now = time(NULL);
    struct tm *local;

    local = localtime(&now);
    strftime(s, 47, "%Y-%m-%d %H:%M", local);

    return s;
}





