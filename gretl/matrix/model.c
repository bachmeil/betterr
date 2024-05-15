#include <common.h>
__import common;
__import globals;
__import dataset;
__import glib;
__import matrix;
__import gretltypes;

/**
 * gretl_model_init:
 * @pmod: pointer to model.
 * @dset: pointer to dataset.
 *
 * Initializes a gretl #MODEL, including setting its pointer members
 * to %NULL. This initialization should be done if the caller has
 * declared a #MODEL struct directly, rather than obtaining a pointer to
 * #MODEL using gretl_model_new() (in which case the initialization is
 * done automatically).
 */

void gretl_model_init (MODEL *pmod, const DATASET *dset)
{
    if (pmod == NULL) return;

#if MDEBUG
    fprintf(stderr, "gretl_model_init: pmod at %p\n", (void *) pmod);
#endif

    pmod->ID = 0;
    pmod->refcount = 0;
    pmod->ci = 0;
    pmod->opt = OPT_NONE;
    pmod->full_n = 0;
    pmod->t1 = pmod->t2 = 0;
    pmod->nobs = 0;

    if (dset != NULL) {
	pmod->smpl.t1 = dset->t1;
	pmod->smpl.t2 = dset->t2;
	pmod->smpl.rseed = dset->rseed;
    } else {
	pmod->smpl.t1 = 0;
	pmod->smpl.t2 = 0;
	pmod->smpl.rseed = 0;
    }

    pmod->ncoeff = 0;
    pmod->ntests = 0;
    pmod->nparams = 0;
    pmod->errcode = 0;
    pmod->ifc = 0;
    pmod->nwt = 0;
    pmod->aux = AUX_NONE;
    pmod->esttime = 0;

    model_stats_init(pmod);

    gretl_model_init_pointers(pmod);
    pmod->n_data_items = 0;
}

/**
 * gretl_model_set_int:
 * @pmod: pointer to #MODEL.
 * @key: key string, used in retrieval.
 * @val: integer value to set.
 *
 * Records an integer value on a model: the value can be retrieved
 * later using gretl_model_get_int(), using the appropriate @key.
 *
 * Returns: 0 on success, 1 on failure.
 */

int gretl_model_set_int (MODEL *pmod, const char *key, int val)
{
    int *valp;
    int err;

    /* if value is already set, reset it */
    valp = gretl_model_get_data(pmod, key);
    if (valp != NULL) {
	*valp = val;
	return 0;
    }

    valp = malloc(sizeof *valp);
    if (valp == NULL) return 1;

    *valp = val;

    err = gretl_model_set_data(pmod, key, valp, GRETL_TYPE_INT,
			       sizeof(int));
    if (err) {
	free(valp);
    }

    return err;
}

/**
 * gretl_model_set_double:
 * @pmod: pointer to model.
 * @key: key string, used in retrieval.
 * @val: double-precision value to set.
 *
 * Records a floating-point value on @pmod: the value can be
 * retrieved later using gretl_model_get_double() with the
 * appropriate @key.
 *
 * Returns: 0 on success, 1 on failure.
 */

int gretl_model_set_double (MODEL *pmod, const char *key, double val)
{
    double *valp;
    int err;

    /* if value is already set, reset it */
    valp = gretl_model_get_data(pmod, key);
    if (valp != NULL) {
	*valp = val;
	return 0;
    }

    valp = malloc(sizeof *valp);
    if (valp == NULL) return 1;

    *valp = val;

    err = gretl_model_set_data(pmod, key, valp, GRETL_TYPE_DOUBLE,
			       sizeof(double));
    if (err) free(valp);

    return err;
}

/**
 * gretl_model_get_int:
 * @pmod: pointer to model.
 * @key: key string.
 *
 * Returns: the integer value identified by @key, or 0 on failure.
 */

int gretl_model_get_int (const MODEL *pmod, const char *key)
{
    int *valp = NULL;
    int i;

    if (pmod == NULL) {
	return 0;
    }

    for (i=0; i<pmod->n_data_items; i++) {
	if (pmod->data_items[i]->type != GRETL_TYPE_INT) {
	    continue;
	}
	if (!strcmp(key, pmod->data_items[i]->key)) {
	    valp = (int *) pmod->data_items[i]->ptr;
	    return *valp;
	}
    }

    return 0;
}

/**
 * gretl_model_get_double:
 * @pmod: pointer to model.
 * @key: key string.
 *
 * Returns: the double-precision value identified by @key, or
 * #NADBL on failure.
 */

double gretl_model_get_double (const MODEL *pmod, const char *key)
{
    double *valp = NULL;
    int i;

    if (pmod == NULL) {
	return NADBL;
    }

    for (i=0; i<pmod->n_data_items; i++) {
	if (pmod->data_items[i]->type != GRETL_TYPE_DOUBLE) {
	    continue;
	}
	if (!strcmp(key, pmod->data_items[i]->key)) {
	    valp = (double *) pmod->data_items[i]->ptr;
	    return *valp;
	}
    }

    return NADBL;
}

int mle_criteria (MODEL *pmod, int addk)
{
    int err = 0;

    if (na(pmod->lnL)) {
	pmod->criterion[C_AIC] = NADBL;
	pmod->criterion[C_BIC] = NADBL;
	pmod->criterion[C_HQC] = NADBL;
	err = 1;
    } else {
	int k = pmod->ncoeff + addk;
	int n = pmod->nobs;

	pmod->criterion[C_AIC] = -2.0 * pmod->lnL + 2.0 * k;
	pmod->criterion[C_BIC] = -2.0 * pmod->lnL + k * log(n);
	pmod->criterion[C_HQC] = -2.0 * pmod->lnL + 2 * k * log(log(n));
    }

    return err;
}

void set_model_id (MODEL *pmod, gretlopt opt)
{
    /* always record model estimation time */
     pmod->esttime = gretl_monotonic_time();

    /* Experimental: limit setting of sequential model
       number to "visible" models estimated in main
       script or session (not in functions).
    */
    if (opt & OPT_A) {
	/* An auxiliary model? Likely, but OPT_A has
	   special meaning for a few estimators
	*/
	if (pmod->ci != DPANEL && pmod->ci != GARCH &&
	    pmod->ci != ARMA) {
	    return;
	}
    }

    if ((opt & OPT_Q) && !(opt & OPT_W)) {
	/* --quiet and not --window, so "invisible" */
	return;
    } else if (gretl_function_depth() > 0) {
	/* model was estimated inside a function */
	return;
    }

    if (pmod->errcode == 0) {
	pmod->ID = ++gretl_model_count;
    }
}

int model_missing (const MODEL *pmod, int t)
{
    if (pmod->missmask != NULL) {
	return pmod->missmask[t] == '1';
    } else {
	return 0;
    }
}

/**
 * gretl_model_set_list_as_data:
 * @pmod: pointer to #MODEL.
 * @key: key string, used in retrieval.
 * @list: list to attach.
 *
 * Attaches @list to @pmod as data, recoverable via the key @key
 * using gretl_model_get_list().  Note that the model takes
 * ownership of the supplied list.
 *
 * Returns: 0 on success, 1 on failure.
 */

int gretl_model_set_list_as_data (MODEL *pmod, const char *key, int *list)
{
    size_t size = (list[0] + 1) * sizeof *list;

    return gretl_model_set_data_with_destructor(pmod, key, (void *) list,
						GRETL_TYPE_LIST, size,
						NULL);
}

/**
 * gretl_model_destroy_data_item:
 * @pmod: pointer to model.
 * @key: key string.
 *
 * Looks up the data pointer, attached to @pmod, that is
 * identified by @key, and if a pointer is found, frees
 * it (or applies the destructor function that was set for
 * the item, if any) and removes it from the model's list of
 * data items.  If you want to remove the item from the
 * model's list without freeing the underlying data pointer,
 * use gretl_model_detach_data_item().
 *
 * Returns: 0 on success, 1 on failure (pointer not found).
 */

int gretl_model_destroy_data_item (MODEL *pmod, const char *key)
{
    return discard_model_data_item(pmod, key, 1);
}

static int discard_model_data_item (MODEL *pmod, const char *key,
				    int free_data)
{
    model_data_item *junk = NULL;
    int i, targ = 0;
    int err = 0;

    for (i=0; i<pmod->n_data_items; i++) {
	if (!strcmp(key, pmod->data_items[i]->key)) {
	    junk = pmod->data_items[i];
	    targ = i;
	    break;
	}
    }

    if (junk == NULL) {
	err = 1;
    } else {
	int n_items = pmod->n_data_items - 1;

	if (n_items == 0) {
	    free(pmod->data_items);
	    pmod->data_items = NULL;
	} else {
	    model_data_item **items;

	    for (i=targ; i<n_items; i++) {
		pmod->data_items[i] = pmod->data_items[i+1];
	    }
	    items = realloc(pmod->data_items, n_items * sizeof *items);
	    if (items != NULL) {
		pmod->data_items = items;
	    }
	}

	pmod->n_data_items -= 1;

	if (free_data) {
	    /* deep free the data item */
	    free_model_data_item(junk);
	} else {
	    /* just free the item, not the actual data */
	    free(junk->key);
	    free(junk);
	}
    }

    return err;
}

/**
 * clear_model:
 * @pmod: pointer to model.
 *
 * Clears a gretl #MODEL, freeing all allocated storage and setting
 * pointer members to %NULL.  Also frees any data pointers attached
 * via gretl_model_set_data().  The model pointer itself is not
 * freed, so this function may be called on a #MODEL which has been
 * declared directly by the caller; in that case the caller should
 * pass the address of the #MODEL in question.
 */

void clear_model (MODEL *pmod)
{
    if (pmod != NULL) {
#if MDEBUG
	fprintf(stderr, "clear model: model at %p\n", (void *) pmod);
#endif

#if MDEBUG > 1
	debug_print_model_info(pmod, "Doing clear_model");
#endif
	if (pmod->list != NULL) free(pmod->list);
	if (pmod->missmask != NULL) free(pmod->missmask);
	if (pmod->coeff != NULL) free(pmod->coeff);
	if (pmod->sderr != NULL) free(pmod->sderr);
	if (pmod->yhat != NULL) free(pmod->yhat);
	if (pmod->uhat != NULL) free(pmod->uhat);
	if (pmod->xpx != NULL) free(pmod->xpx);
	if (pmod->vcv != NULL) free(pmod->vcv);
	if (pmod->name != NULL) free(pmod->name);
	if (pmod->depvar != NULL) free(pmod->depvar);

	if (pmod->submask != NULL) {
	    free_subsample_mask(pmod->submask);
	}

	if (pmod->arinfo != NULL) {
	    clear_ar_info(pmod);
	}
	if (pmod->params != NULL) {
	    strings_array_free(pmod->params, pmod->nparams);
	}

	destroy_dataset(pmod->dataset);
	gretl_model_destroy_tests(pmod);
	destroy_all_data_items(pmod);
    }

    /* this may be redundant */
#if MDEBUG
    fprintf(stderr, "clear_model, calling gretl_model_init\n");
#endif
#if 1
    gretl_model_init(pmod, NULL);
#endif
}

/**
 * gretl_model_set_string_as_data:
 * @pmod: pointer to #MODEL.
 * @key: key string, used in retrieval.
 * @str: string to attach.
 *
 * Attaches @str to @pmod as data, recoverable via the key @key
 * using gretl_model_get_data().
 *
 * Returns: 0 on success, 1 on failure.
 */

int gretl_model_set_string_as_data (MODEL *pmod, const char *key, char *str)
{
    size_t size = strlen(str) + 1;

    return gretl_model_set_data_with_destructor(pmod, key, (void *) str,
						GRETL_TYPE_STRING, size,
						NULL);
}

static void model_stats_init (MODEL *pmod)
{
    int i;

    pmod->ess = pmod->tss = NADBL;
    pmod->sigma = NADBL;
    pmod->rsq = pmod->adjrsq = NADBL;
    pmod->fstt = pmod->chisq = NADBL;
    pmod->lnL = NADBL;
    pmod->ybar = pmod->sdy = NADBL;
    pmod->dw = pmod->rho = NADBL;

    for (i=0; i<C_MAX; i++) {
	pmod->criterion[i] = NADBL;
    }
}

static void gretl_model_init_pointers (MODEL *pmod)
{
    pmod->list = NULL;
    pmod->submask = NULL;
    pmod->missmask = NULL;
    pmod->coeff = NULL;
    pmod->sderr = NULL;
    pmod->yhat = NULL;
    pmod->uhat = NULL;
    pmod->xpx = NULL;
    pmod->vcv = NULL;
    pmod->arinfo = NULL;
    pmod->name = NULL;
    pmod->depvar = NULL;
    pmod->params = NULL;
    pmod->tests = NULL;
    pmod->dataset = NULL;
    pmod->data_items = NULL;
}

/**
 * gretl_model_get_data:
 * @pmod: pointer to model.
 * @key: key string.
 *
 * Returns: the data pointer identified by @key, or %NULL on failure.
 */

void *gretl_model_get_data (const MODEL *pmod, const char *key)
{
    return gretl_model_get_data_full(pmod, key, NULL, NULL, NULL);
}

/**
 * gretl_model_set_data:
 * @pmod: pointer to #MODEL.
 * @key: key string for data, used in retrieval.
 * @ptr: data-pointer to be attached to model.
 * @type: type of the data to set.
 * @size: size of data in bytes.
 *
 * Attaches data to @pmod: the data can be retrieved later using
 * gretl_model_get_data().  Note that the data are not "physically"
 * copied to the model; simply, @ptr is recorded on the model.
 * This means that the data referenced by the pointer now in
 * effect belong to @pmod.  The data pointer will be freed when
 * @pmod is cleared with clear_model().  If the data has deep
 * structure that requires special treatment on freeing, use
 * gretl_model_set_data_with_destructor() instead.
 *
 * The @size is needed in case the model is copied with
 * copy_model(), in which case the target of the copying
 * operation receives a newly allocated copy of the data in
 * question.
 *
 * Returns: 0 on success, 1 on failure.
 */

int gretl_model_set_data (MODEL *pmod, const char *key, void *ptr,
			  GretlType type, size_t size)
{
    return gretl_model_set_data_with_destructor(pmod, key, ptr, type,
						size, NULL);
}

gint64 gretl_monotonic_time (void)
{
    return g_get_monotonic_time();
}

/**
 * gretl_model_set_data_with_destructor:
 * @pmod: pointer to #MODEL.
 * @key: key string for data, used in retrieval.
 * @ptr: data-pointer to be attached to model.
 * @type: type of data to set.
 * @size: size of data in bytes.
 * @destructor: pointer to function that should be used to free
 * the data-pointer in question.
 *
 * Attaches data to @pmod: the data can be retrieved later using
 * gretl_model_get_data().  Note that the data are not "physically"
 * copied to the model; simply, @ptr is recorded on the model.
 * This means that the data referenced by the pointer now in
 * effect belong to @pmod.  When @pmod is cleared with clear_model(),
 * @destructor will be invoked with @ptr as its single argument.
 * If a simple "free" is OK for freeing the data, you can use
 * gretl_model_set_data() instead.
 *
 * The @size is needed in case the model is copied with
 * copy_model(), in which case the target of the copying
 * operation receives a newly allocated copy of the data in
 * question.
 *
 * Returns: 0 on success, 1 on failure.
 */

int gretl_model_set_data_with_destructor (MODEL *pmod, const char *key, void *ptr,
					  GretlType type, size_t size,
					  void (*destructor) (void *))
{
    model_data_item **items;
    model_data_item *item;
    int i, n;

    for (i=0; i<pmod->n_data_items; i++) {
	item = pmod->data_items[i];
	if (!strcmp(key, item->key)) {
	    /* there's a pre-existing item of this name */
	    if (item->destructor != NULL) {
		(*item->destructor)(item->ptr);
	    } else {
		free(item->ptr);
	    }
	    item->type = type;
	    item->ptr = ptr;
	    item->size = size;
	    item->destructor = destructor;
	    /* handled */
	    return 0;
	}
    }

    n = pmod->n_data_items + 1;

    items = realloc(pmod->data_items, n * sizeof *items);
    if (items == NULL) {
	return 1;
    }

    pmod->data_items = items;

    item = create_data_item(key, ptr, type, size, destructor);
    if (item == NULL) {
	return 1;
    }

    pmod->data_items[n - 1] = item;
    pmod->n_data_items += 1;

    return 0;
}

static void free_model_data_item (model_data_item *item)
{
    if (item->destructor != NULL) {
	(*item->destructor)(item->ptr);
    } else {
	free(item->ptr);
    }
    free(item->key);
    free(item);
}

static void clear_ar_info (MODEL *pmod)
{
    if (pmod->arinfo->arlist) {
	free(pmod->arinfo->arlist);
    }
    if (pmod->arinfo->rho) {
	free(pmod->arinfo->rho);
    }
    if (pmod->arinfo->sderr) {
	free(pmod->arinfo->sderr);
    }

    free(pmod->arinfo);
    pmod->arinfo = NULL;
}

/**
 * gretl_model_destroy_tests:
 * @pmod: pointer to model.
 *
 * Clears any hypothesis test structs that have been attached
 * to @pmod.
 */

void gretl_model_destroy_tests (MODEL *pmod)
{
    if (pmod != NULL && pmod->ntests > 0) {
	int i;

	for (i=0; i<pmod->ntests; i++) {
	    if (pmod->tests[i].param != NULL) {
		free(pmod->tests[i].param);
	    }
	}
	free(pmod->tests);
	pmod->tests = NULL;
	pmod->ntests = 0;
    }
}

static void destroy_all_data_items (MODEL *pmod)
{
    int i;

    if (pmod->n_data_items == 0) {
	return;
    }

#ifdef HAVE_X12A
    maybe_delete_x12_file(pmod);
#endif

    for (i=0; i<pmod->n_data_items; i++) {
	free_model_data_item(pmod->data_items[i]);
    }

    free(pmod->data_items);
    pmod->data_items = NULL;
}

/**
 * gretl_model_get_data_full:
 * @pmod: pointer to model.
 * @key: key string.
 * @copied: location to receive flag indicating whether the
 * return value is an allocated copy of the original data.
 * @type: location to receive data type.
 * @sz: location to receive the size of the data.
 *
 * Returns: the data pointer identified by @key, or %NULL on failure.
 * If a non-zero value is written to @copied this indicates that the
 * return value is a copy of the original (and therefore it is the
 * caller's responsibility to free the data when it is no longer
 * required).
 */

void *gretl_model_get_data_full (const MODEL *pmod, const char *key,
				 GretlType *type, int *copied,
				 size_t *sz)
{
    void *ret = NULL;
    GretlType itype = 0;
    size_t isize = 0;
    int alloced = 0;
    int i, found = 0;

    if (pmod == NULL) {
	return NULL;
    }

    for (i=0; i<pmod->n_data_items; i++) {
	if (!strcmp(key, pmod->data_items[i]->key)) {
	    ret = pmod->data_items[i]->ptr;
	    itype = pmod->data_items[i]->type;
	    isize = pmod->data_items[i]->size;
	    found = 1;
	    break;
	}
    }

    if (!found && pmod->tests != NULL) {
	const ModelTest *test;
	const char *tkey;
	gretl_bundle *b;

	for (i=0; i<pmod->ntests; i++) {
	    test = &pmod->tests[i];
	    tkey = test_type_key(test->type);
	    if (tkey != NULL && !strcmp(key, tkey)) {
		b = bundlize_test(test);
		if (b != NULL) {
		    ret = b;
		    itype = GRETL_TYPE_BUNDLE;
		    alloced = 1;
		}
		break;
	    }
	}
    }

    if (ret != NULL) {
	if (type != NULL) {
	    *type = itype;
	}
	if (sz != NULL) {
	    *sz = isize;
	}
	if (copied != NULL) {
	    *copied = alloced;
	}
    }

    return ret;
}

static model_data_item *create_data_item (const char *key, void *ptr,
					  GretlType type, size_t size,
					  void (*destructor) (void *))
{
    model_data_item *item = malloc(sizeof *item);

    if (item != NULL) {
	item->key = gretl_strdup(key);
	if (item->key == NULL) {
	    free(item);
	    item = NULL;
	} else {
	    item->ptr = ptr;
	    item->type = type;
	    item->size = size;
	    item->destructor = destructor;
	}
    }

    return item;
}

static const char *test_type_key (ModelTestType t)
{
    if (t == GRETL_TEST_ADD) {
	return "add_test";
    } else if (t == GRETL_TEST_ARCH) {
	return "arch_test";
    } else if (t == GRETL_TEST_AUTOCORR) {
	return "autocorr_test";
    } else if (t == GRETL_TEST_CHOW ||
	       t == GRETL_TEST_CHOWDUM) {
	return "chow_test";
    } else if (t == GRETL_TEST_CUSUM) {
	return "cusum_test";
    } else if (t == GRETL_TEST_QLR) {
	return "qlr_test";
    } else if (t == GRETL_TEST_GROUPWISE) {
	return "grpwise_test";
    } else if (t == GRETL_TEST_LOGS) {
	return "logs_test";
    } else if (t == GRETL_TEST_NORMAL) {
	return "normality_test";
    } else if (t == GRETL_TEST_OMIT) {
	return "omit_test";
    } else if (t == GRETL_TEST_RESET) {
	return "reset_test";
    } else if (t == GRETL_TEST_SQUARES) {
	return "squares_test";
    } else if (t == GRETL_TEST_WHITES) {
	return "whites_test";
    } else if (t == GRETL_TEST_SARGAN) {
	return "sargan_test";
    } else if (t == GRETL_TEST_IV_HAUSMAN ||
	       t == GRETL_TEST_PANEL_HAUSMAN) {
	return "hausman_test";
    } else if (t == GRETL_TEST_PANEL_F ||
	       t == GRETL_TEST_PANEL_WELCH) {
	return "fixed_effects_F";
    } else if (t == GRETL_TEST_PANEL_BP ||
	       t == GRETL_TEST_BP) {
	return "bp_test";
    } else if (t == GRETL_TEST_PANEL_TIMEDUM) {
	return "timedum_test";
    } else if (t == GRETL_TEST_HET_1) {
	return "het1_test";
    } else if (t == GRETL_TEST_COMFAC) {
	return "comfac_test";
    } else if (t == GRETL_TEST_INDEP) {
	return "independence_test";
    } else if (t == GRETL_TEST_RE) {
	return "rho_test";
    } else if (t == GRETL_TEST_WITHIN_F) {
	return "within_F";
    } else if (t == GRETL_TEST_RE_WALD) {
	return "re_wald_test";
    } else if (t == GRETL_TEST_XDEPEND) {
	return "cross_sectional_dependence_test";
    } else if (t == GRETL_TEST_PANEL_AR) {
	return "panel_ar_test";
    } else {
	fprintf(stderr, "test_type_key(): type %d has no key!\n", t);
	return NULL;
    }
}

static gretl_bundle *bundlize_test (const ModelTest *src)
{
    gretl_bundle *b = gretl_bundle_new();

    if (b == NULL) {
	return NULL;
    }

    if (src->param != NULL && *src->param != '\0') {
	gretl_bundle_set_string(b, "param", src->param);
    }

    if (src->dfn > 0) {
	gretl_bundle_set_scalar(b, "dfn", (double) src->dfn);
    }
    if (src->dfd > 0) {
	gretl_bundle_set_scalar(b, "dfd", src->dfd);
    }
    if (src->order > 0) {
	gretl_bundle_set_scalar(b, "order", src->order);
    }
    if (!na(src->value)) {
	gretl_bundle_set_scalar(b, "test", src->value);
    }
    if (!na(src->pvalue)) {
	gretl_bundle_set_scalar(b, "pvalue", src->pvalue);
    }
    if (!na(src->crit)) {
	gretl_bundle_set_scalar(b, "crit", src->crit);
	gretl_bundle_set_scalar(b, "alpha", src->alpha);
    }

    return b;
}

gretl_bundle *gretl_bundle_new (void) { return NULL; }

/**
 * gretl_bundle_new:
 *
 * Returns: a newly allocated, empty gretl bundle.
 */

//~ gretl_bundle *gretl_bundle_new (void)
//~ {
    //~ gretl_bundle *b = malloc(sizeof *b);

    //~ if (b != NULL) {
        //~ b->type = BUNDLE_PLAIN;
        //~ b->ht = g_hash_table_new_full(g_str_hash, g_str_equal,
                                      //~ NULL, bundle_item_destroy);
        //~ b->creator = NULL;
        //~ b->data = NULL;
    //~ }

    //~ return b;
//~ }

/**
 * gretl_bundle_set_string:
 * @bundle: target bundle.
 * @key: name of key to create or replace.
 * @str: the string to set.
 *
 * Sets @str as a member of @bundle under the name @key.
 * If @key is already present in the bundle the original
 * value is replaced and destroyed.
 *
 * Returns: 0 on success, error code on error.
 */

int gretl_bundle_set_string (gretl_bundle *bundle, const char *key,
                             const char *str)
{
    return gretl_bundle_set_data(bundle, key, (void *) str,
                                 GRETL_TYPE_STRING, 0);
}

/**
 * gretl_bundle_set_scalar:
 * @bundle: target bundle.
 * @key: name of key to create or replace.
 * @val: the value to set.
 *
 * Sets @val as a member of @bundle under the name @key.
 * If @key is already present in the bundle the original
 * value is replaced and destroyed.
 *
 * Returns: 0 on success, error code on error.
 */

int gretl_bundle_set_scalar (gretl_bundle *bundle, const char *key,
                             double val)
{
    return gretl_bundle_set_data(bundle, key, &val,
                                 GRETL_TYPE_DOUBLE, 0);
}

static int real_bundle_set_data (gretl_bundle *b, const char *key,
                                 void *ptr, GretlType type,
                                 int size, int copy,
                                 const char *note) {
	return 0;
}

//~ static int real_bundle_set_data (gretl_bundle *b, const char *key,
                                 //~ void *ptr, GretlType type,
                                 //~ int size, int copy,
                                 //~ const char *note)
//~ {
    //~ int err = 0, done = 0;

    //~ if (key == NULL || key[0] == '\0') {
        //~ gretl_errmsg_sprintf("real_bundle_set_data: missing key string");
        //~ return E_DATA;
    //~ }

    //~ if (b->type == BUNDLE_KALMAN) {
        //~ done = maybe_set_kalman_element(b->data, key,
                                        //~ ptr, type, copy,
                                        //~ &err);
    //~ }

    //~ if (!done && !err) {
        //~ bundled_item *item = g_hash_table_lookup(b->ht, key);
        //~ int replace = 0;

        //~ if (item != NULL) {
            //~ replace = 1;
            //~ if (item->type == type) {
                //~ /* we can take a shortcut */
                //~ return bundled_item_replace_data(item, ptr, type, size, copy);
            //~ }
        //~ }

        //~ item = bundled_item_new(type, ptr, size, copy, note, &err);

        //~ if (!err) {
            //~ item->name = g_strdup(key);
            //~ if (replace) {
                //~ g_hash_table_replace(b->ht, item->name, item);
            //~ } else {
                //~ g_hash_table_insert(b->ht, item->name, item);
            //~ }
        //~ }
    //~ }

    //~ return err;
//~ }

/**
 * gretl_bundle_set_data:
 * @bundle: target bundle.
 * @key: name of key to create or replace.
 * @ptr: data pointer.
 * @type: type of data.
 * @size: if @type == GRETL_TYPE_SERIES, the length of
 * the series, otherwise 0.
 *
 * Sets the data type and pointer to be associated with @key in
 * the bundle given by @name. If @key is already present in
 * the bundle's hash table the original value is replaced
 * and destroyed. The content of @ptr is copied into the
 * bundle; compare gretl_bundle_donate_data().
 *
 * Returns: 0 on success, error code on error.
 */

int gretl_bundle_set_data (gretl_bundle *bundle, const char *key,
                           void *ptr, GretlType type, int size)
{
    int err;

    if (bundle == NULL) {
        err = E_UNKVAR;
    } else {
        err = real_bundle_set_data(bundle, key, ptr, type,
                                   size, 1, NULL);
    }

    return err;
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
 * gretl_model_get_double_default:
 * @pmod: pointer to model.
 * @key: key string.
 * @deflt: default value
 *
 * Returns: the double-precision value identified by @key, or
 * @deflt if there is no such value.
 */

double gretl_model_get_double_default (const MODEL *pmod,
				       const char *key,
				       double deflt)
{
    double *valp = NULL;
    double ret = deflt;
    int i;

    if (pmod == NULL) {
	return NADBL;
    }

    for (i=0; i<pmod->n_data_items; i++) {
	if (pmod->data_items[i]->type != GRETL_TYPE_DOUBLE) {
	    continue;
	}
	if (!strcmp(key, pmod->data_items[i]->key)) {
	    valp = (double *) pmod->data_items[i]->ptr;
	    ret = *valp;
	    break;
	}
    }

    return ret;
}

/**
 * gretl_model_set_hac_vcv_info:
 * @pmod: pointer to model.
 * @kern: kernel type.
 * @order: lag order.
 * @flags: bitflags.
 * @bw: QS bandwidth, or #NADBL if not applicable.
 *
 * Returns: 0 on success, 1 on failure.
 */

int gretl_model_set_hac_vcv_info (MODEL *pmod, int kern,
				  int order, int flags,
				  double bw)
{
    return gretl_model_set_full_vcv_info(pmod, VCV_HAC, kern,
					 order, flags, bw,
					 NULL, NULL);
}

/**
 * gretl_model_set_full_vcv_info:
 *
 * Returns: 0 on success, 1 on failure.
 */

static int
gretl_model_set_full_vcv_info (MODEL *pmod, int vmaj, int vmin,
			       int order, int flags, double bw,
			       const char *cv1, const char *cv2)
{
    VCVInfo *vi;
    int prev = 0;
    int err = 0;

    vi = gretl_model_get_data(pmod, "vcv_info");

    if (vi == NULL) {
	vi = vcv_info_new();
	if (vi == NULL) {
	    return E_ALLOC;
	}
    } else {
	prev = 1;
	free(vi->cv1);
	free(vi->cv2);
	vi->cv1 = vi->cv2 = NULL;
    }

    vi->vmaj = vmaj;
    vi->vmin = vmin;
    vi->order = order;
    vi->flags = flags;
    vi->bw = bw;

    if (cv1 != NULL) {
	vi->cv1 = gretl_strdup(cv1);
    }
    if (cv2 != NULL) {
	vi->cv2 = gretl_strdup(cv2);
    }

    if (!prev) {
	err = gretl_model_set_data_with_destructor(pmod, "vcv_info", vi,
						   GRETL_TYPE_STRUCT,
						   sizeof *vi,
						   vcv_info_free);
    }

    return err;
}

/**
 * gretl_model_set_vcv_info:
 * @pmod: pointer to model.
 * @vmaj: top-level VCV type.
 * @vmin: variant under @vmaj, if applicable.
 *
 * Returns: 0 on success, 1 on failure.
 */

int gretl_model_set_vcv_info (MODEL *pmod, int vmaj, int vmin)
{
    return gretl_model_set_full_vcv_info(pmod, vmaj, vmin,
					 0, 0, 0, NULL, NULL);
}

/**
 * gretl_model_write_vcv:
 * @pmod: pointer to model.
 * @V: covariance matrix.
 *
 * Write the covariance matrix @V into the model @pmod, using the
 * special packed format that is required by the MODEL struct,
 * and set the standard errors to the square root of the diagonal
 * elements of this matrix.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int gretl_model_write_vcv (MODEL *pmod, const gretl_matrix *V)
{
    int i, j, k, n;
    double x, *tmp;
    int err = 0;

    if (gretl_is_null_matrix(V)) {
	return 0; /* no-op */
    }

    if (V->cols != V->rows) {
	return E_NONCONF;
    }

    k = V->rows;
    n = (k * k + k) / 2;

    /* reallocate model vcv in case it's wrongly sized */
    tmp = realloc(pmod->vcv, n * sizeof *tmp);
    if (tmp == NULL) {
	err = E_ALLOC;
    } else {
	pmod->vcv = tmp;
    }

    /* same for standard errors array */
    if (!err) {
	tmp = realloc(pmod->sderr, k * sizeof *tmp);
	if (tmp == NULL) {
	    err = E_ALLOC;
	} else {
	    pmod->sderr = tmp;
	}
    }

    if (!err) {
	int restricted, idx = 0;

	restricted = gretl_model_get_int(pmod, "restricted");

	for (i=0; i<k; i++) {
	    for (j=i; j<k; j++) {
		x = gretl_matrix_get(V, i, j);
		pmod->vcv[idx++] = x;
		if (i == j) {
		    pmod->sderr[i] = vcv_get_se(x, restricted);
		}
	    }
	}
    }

    return err;
}

/**
 * gretl_model_set_cluster_vcv_info:
 * @pmod: pointer to model.
 * @cv1: name of (first) cluster variable.
 * @cv2: name of second cluster variable or NULL.
 *
 * Returns: 0 on success, 1 on failure.
 */

int gretl_model_set_cluster_vcv_info (MODEL *pmod,
				      const char *cv1,
				      const char *cv2)
{
    return gretl_model_set_full_vcv_info(pmod, VCV_CLUSTER,
					 0, 0, 0, 0,
					 cv1, cv2);
}

/**
 * gretl_model_set_hac_order:
 * @pmod: pointer to model.
 * @order: lag order.
 *
 * Returns: 0 on success, 1 on failure.
 */

int gretl_model_set_hac_order (MODEL *pmod, int order)
{
    VCVInfo *vi = gretl_model_get_data(pmod, "vcv_info");

    if (vi != NULL) {
	vi->order = order;
	return 0;
    } else {
	return E_DATA;
    }
}

static VCVInfo *vcv_info_new (void)
{
    VCVInfo *vi;

    vi = malloc(sizeof *vi);

    if (vi != NULL) {
	vi->vmaj = vi->vmin = 0;
	vi->order = vi->flags = 0;
	vi->bw = NADBL;
	vi->cv1 = NULL;
	vi->cv2 = NULL;
    }

    return vi;
}

static void vcv_info_free (void *data)
{
    VCVInfo *vi = data;

    if (vi != NULL) {
	free(vi->cv1);
	free(vi->cv2);
	free(vi);
    }
}

static double vcv_get_se (double vii, int restricted)
{
    if (restricted) {
	vii = fabs(vii) < 1.0e-17 ? 0.0 : vii;
    }

    return (na(vii) || vii < 0)? NADBL : sqrt(vii);
}

/**
 * gretl_model_set_matrix_as_data:
 * @pmod: pointer to #MODEL.
 * @key: key string, used in retrieval.
 * @m: matrix to attach.
 *
 * Attaches @m to @pmod as data, recoverable via the key @key
 * using gretl_model_get_data().
 *
 * Returns: 0 on success, 1 on failure.
 */

int gretl_model_set_matrix_as_data (MODEL *pmod, const char *key,
				    gretl_matrix *m)
{
    return gretl_model_set_data_with_destructor(pmod, key, (void*) m,
						GRETL_TYPE_MATRIX, 0,
						matrix_free_callback);
}

/**
 * gretl_model_add_y_median:
 * @pmod: pointer to target model.
 * @y: array containing the dependent variable.
 *
 * Calculates the median of @y using the valid observations
 * with the model's sample range and attaches the median
 * to the model as data under the key %ymedian.
 *
 * Returns: 0 on success or error code on error.
 */

int gretl_model_add_y_median (MODEL *pmod, const double *y)
{
    int T = pmod->t2 - pmod->t1 + 1;
    double *sy, m;
    int t, n, ok, n2p;

    sy = malloc(T * sizeof *sy);

    if (sy == NULL) {
	return E_ALLOC;
    }

    n = 0;
    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (pmod->uhat != NULL) {
	    ok = !na(pmod->uhat[t]);
	} else {
	    ok = !model_missing(pmod, t);
	}
	if (ok) {
	    sy[n++] = y[t];
	}
    }

    if (n == 0) {
	free(sy);
	return E_DATA;
    }

    qsort(sy, n, sizeof *sy, gretl_compare_doubles);

    n2p = (T = n / 2) + 1;
    m = (n % 2)? sy[n2p - 1] : 0.5 * (sy[T - 1] + sy[n2p - 1]);

    gretl_model_set_double(pmod, "ymedian", m);

    free(sy);

    return 0;
}

//~ int gretl_model_add_normality_test (MODEL *pmod, double X2)
//~ {
    //~ ModelTest *test = model_test_new(GRETL_TEST_NORMAL);
    //~ int err = 0;

    //~ if (test != NULL) {
        //~ model_test_set_teststat(test, GRETL_STAT_NORMAL_CHISQ);
        //~ model_test_set_dfn(test, 2);
        //~ model_test_set_value(test, X2);
        //~ model_test_set_pvalue(test, chisq_cdf_comp(2, X2));
        //~ maybe_add_test_to_model(pmod, test);
    //~ } else {
        //~ err = E_ALLOC;
    //~ }

    //~ return err;
//~ }

/**
 * gretl_model_new_vcv:
 * @pmod: pointer to model.
 * @nelem: pointer to receive number of elements in
 * the packed array, or %NULL;
 *
 * Allocates space for a packed coefficient covariance matrix
 * in @pmod (if such space is not already allocated).  Sets
 * all entries in the array to zero.
 *
 * Returns: 0 on success, %E_ALLOC on error.
 */

int gretl_model_new_vcv (MODEL *pmod, int *nelem)
{
    int nv = pmod->ncoeff;
    int nxpx = (nv * nv + nv) / 2;
    int i, err = 0;

    if (pmod->vcv == NULL) {
	/* not already allocated */
	pmod->vcv = malloc(nxpx * sizeof *pmod->vcv);
	if (pmod->vcv == NULL) {
	    err = E_ALLOC;
	}
    }

    if (pmod->vcv != NULL) {
	for (i=0; i<nxpx; i++) {
	    pmod->vcv[i] = 0.0;
	}
	if (nelem != NULL) {
	    *nelem = nxpx;
	}
    }

    return err;
}

/**
 * model_has_missing_obs:
 * @pmod: pointer to model.
 *
 * Returns: 1 if there are missing observations in the
 * model's sample range, otherwise 0.
 */

int model_has_missing_obs (const MODEL *pmod)
{
    int sample = pmod->t2 - pmod->t1 + 1;

    return (pmod->nobs < sample);
}

static void matrix_free_callback (void *p)
{
    gretl_matrix_free((struct gretl_matrix_ *) p);
}

/**
 * model_test_new:
 * @ttype: type of test to add.
 *
 * Returns: new #ModelTest pointer, or %NULL on failure.
 */

ModelTest *model_test_new (ModelTestType ttype)
{
    ModelTest *test = malloc(sizeof *test);

    if (test != NULL) {
	gretl_test_init(test, ttype);
    }

    return test;
}

/**
 * maybe_add_test_to_model:
 * @pmod: pointer to model.
 * @test: model test to be added.
 *
 * Adds a #ModelTest to @pmod, if the test in question has
 * not already been performed and recorded.  Note that this
 * function takes care of freeing @test.
 *
 * Returns: 1 if the test was added, otherwise 0.
 */

int maybe_add_test_to_model (MODEL *pmod, ModelTest *test)
{
    int i, done = 0, add = 0;

    if (test == NULL || test->teststat == GRETL_STAT_NONE) {
	return 0;
    }

    for (i=0; i<pmod->ntests; i++) {
	if (!model_tests_differ(test, &pmod->tests[i])) {
	    done = 1;
	}
    }

    if (!done) {
	int n = pmod->ntests + 1;
	ModelTest *tests;

	tests = realloc(pmod->tests, n * sizeof *tests);

	if (tests != NULL) {
	    pmod->tests = tests;
	    copy_test(&pmod->tests[n-1], test);
	    pmod->ntests += 1;
	    add = 1;
	}
    }

    free(test->param);
    free(test);

    return add;
}

void model_test_set_teststat (ModelTest *test, unsigned char ts)
{
    test->teststat = ts;
}

void model_test_set_order (ModelTest *test, int order)
{
    test->order = order;
}

void model_test_set_dfn (ModelTest *test, int df)
{
    test->dfn = df;
}

void model_test_set_dfd (ModelTest *test, double df)
{
    test->dfd = df;
}

void model_test_set_value (ModelTest *test, double val)
{
    test->value = val;
}

void model_test_set_pvalue (ModelTest *test, double pval)
{
    test->pvalue = pval;
}

void model_test_set_opt (ModelTest *test, gretlopt opt)
{
    test->opt = opt;
}

void model_test_set_crit_and_alpha (ModelTest *test,
				    double crit,
				    double alpha)
{
    test->crit = crit;
    test->alpha = alpha;
}

void model_test_set_param (ModelTest *test, const char *s)
{
    test->param = gretl_strdup(s);
}

void model_test_set_allocated_param (ModelTest *test, char *s)
{
    test->param = s;
}

/**
 * chisq_cdf_comp:
 * @df: degrees of freedom.
 * @x: the cutoff point in the distribution.
 *
 * Returns: the integral from @x to infinity of the chi-square
 * distribution with @df degrees of freedom, or #NADBL
 * on failure.
 */

//~ double chisq_cdf_comp (double df, double x)
//~ {
    //~ double p = chdtrc(df, x);

    //~ if (get_cephes_errno() == CEPHES_DOMAIN) {
	//~ p = NADBL;
    //~ }

    //~ return p;
//~ }

static void gretl_test_init (ModelTest *test, ModelTestType ttype)
{
    test->type = ttype;
    test->order = 0;
    test->param = NULL;
    test->teststat = GRETL_STAT_NONE;
    test->dfn = 0;
    test->dfd = 0;
    test->value = test->pvalue = NADBL;
    test->crit = test->alpha = NADBL;
    test->opt = OPT_NONE;
}

static int model_tests_differ (ModelTest *mt1, ModelTest *mt2)
{
    int ret = 0;

    if (mt1->type != mt2->type) {
	ret = 1;
    } else if (mt1->order != mt2->order) {
	ret = 1;
    } else if (mt1->teststat != mt2->teststat) {
	ret = 1;
    } else if (test_params_differ(mt1->param, mt2->param)) {
	ret = 1;
    } else if (testvals_differ(mt1->value, mt2->value)) {
	ret = 1;
    }

    return ret;
}

static void copy_test (ModelTest *targ, const ModelTest *src)
{
    targ->type = src->type;

    if (src->param != NULL && *src->param != '\0') {
	targ->param = gretl_strdup(src->param);
    } else {
	targ->param = NULL;
    }

    targ->teststat = src->teststat;
    targ->dfn = src->dfn;
    targ->dfd = src->dfd;
    targ->order = src->order;
    targ->value = src->value;
    targ->pvalue = src->pvalue;
    targ->crit = src->crit;
    targ->alpha = src->alpha;
    targ->opt = src->opt;
}

static int test_params_differ (const char *p, const char *s)
{
    int ret = 0;

    if (s == NULL && p != NULL) {
	ret = 1;
    } else if (p == NULL && s != NULL) {
	ret = 1;
    } else if (p != NULL && s != NULL) {
	ret = strcmp(p, s);
    }

    return ret;
}

static int testvals_differ (double x, double y)
{
    double eq_tol = 1.0e-10;
    double reldiff;

    if (x == 0.0) {
	reldiff = fabs(y);
    } else if (y == 0.0) {
	reldiff = fabs(x);
    } else if (x > y) {
	reldiff = fabs((x - y) / y);
    } else {
	reldiff = fabs((y - x) / x);
    }

    return reldiff > eq_tol;
}
