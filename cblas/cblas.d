/* Code borrowed from
 * https://github.com/DlangScience/cblas
 * This file is a dump of most of the code from there
 * The license remains the original Boost 1.0 found in that repo
 */
module cblas.cblas;

///
alias sdsdot = cblas_sdsdot;
///
alias dsdot = cblas_dsdot;
///
alias dot = cblas_sdot;
///
alias dot = cblas_ddot;
///
alias dotu = cblas_cdotu;
///
alias dotc = cblas_cdotc;
///
alias dotu = cblas_zdotu;
///
alias dotc = cblas_zdotc;
///
alias dotu_sub = cblas_cdotu_sub;
///
alias dotc_sub = cblas_cdotc_sub;
///
alias dotu_sub = cblas_zdotu_sub;
///
alias dotc_sub = cblas_zdotc_sub;
///
alias asum = cblas_sasum;
///
alias asum = cblas_dasum;
///
alias asum = cblas_scasum;
///
alias asum = cblas_dzasum;
///
alias nrm2 = cblas_snrm2;
///
alias nrm2 = cblas_dnrm2;
///
alias nrm2 = cblas_scnrm2;
deprecated("use 'nrm2' alias instead") alias cnrm2 = cblas_scnrm2;
///
alias nrm2 = cblas_dznrm2;
deprecated("use 'nrm2' alias instead") alias znrm2 = cblas_dznrm2;
///
alias iamax = cblas_isamax;
///
alias iamax = cblas_idamax;
///
alias iamax = cblas_icamax;
///
alias iamax = cblas_izamax;
///
alias axpy = cblas_saxpy;
///
alias axpy = cblas_daxpy;
///
alias axpy = cblas_caxpy;
///
alias axpy = cblas_zaxpy;
///
alias copy = cblas_scopy;
///
alias copy = cblas_dcopy;
///
alias copy = cblas_ccopy;
///
alias copy = cblas_zcopy;
///
alias swap = cblas_sswap;
///
alias swap = cblas_dswap;
///
alias swap = cblas_cswap;
///
alias swap = cblas_zswap;
///
alias rot = cblas_srot;
///
alias rot = cblas_drot;
///
alias rotg = cblas_srotg;
///
alias rotg = cblas_drotg;
///
alias rotm = cblas_srotm;
///
alias rotm = cblas_drotm;
///
alias rotmg = cblas_srotmg;
///
alias rotmg = cblas_drotmg;
///
alias scal = cblas_sscal;
///
alias scal = cblas_dscal;
///
alias scal = cblas_cscal;
///
alias scal = cblas_zscal;
///
alias sscal = cblas_csscal;
///
alias dscal = cblas_zdscal;
///
alias gemv = cblas_sgemv;
///
alias gemv = cblas_dgemv;
///
alias gemv = cblas_cgemv;
///
alias gemv = cblas_zgemv;
///
alias ger = cblas_sger;
///
alias ger = cblas_dger;
///
alias geru = cblas_cgeru;
///
alias gerc = cblas_cgerc;
///
alias geru = cblas_zgeru;
///
alias gerc = cblas_zgerc;
///
alias trsv = cblas_strsv;
///
alias trsv = cblas_dtrsv;
///
alias trsv = cblas_ctrsv;
///
alias trsv = cblas_ztrsv;
///
alias trmv = cblas_strmv;
///
alias trmv = cblas_dtrmv;
///
alias trmv = cblas_ctrmv;
///
alias trmv = cblas_ztrmv;
///
alias syr = cblas_ssyr;
///
alias syr = cblas_dsyr;
///
alias her = cblas_cher;
///
alias her = cblas_zher;
///
alias syr2 = cblas_ssyr2;
///
alias syr2 = cblas_dsyr2;
///
alias her2 = cblas_cher2;
///
alias her2 = cblas_zher2;
///
alias gbmv = cblas_sgbmv;
///
alias gbmv = cblas_dgbmv;
///
alias gbmv = cblas_cgbmv;
///
alias gbmv = cblas_zgbmv;
///
alias sbmv = cblas_ssbmv;
///
alias sbmv = cblas_dsbmv;
///
alias tbmv = cblas_stbmv;
///
alias tbmv = cblas_dtbmv;
///
alias tbmv = cblas_ctbmv;
///
alias tbmv = cblas_ztbmv;
///
alias tbsv = cblas_stbsv;
///
alias tbsv = cblas_dtbsv;
///
alias tbsv = cblas_ctbsv;
///
alias tbsv = cblas_ztbsv;
///
alias tpmv = cblas_stpmv;
///
alias tpmv = cblas_dtpmv;
///
alias tpmv = cblas_ctpmv;
///
alias tpmv = cblas_ztpmv;
///
alias tpsv = cblas_stpsv;
///
alias tpsv = cblas_dtpsv;
///
alias tpsv = cblas_ctpsv;
///
alias tpsv = cblas_ztpsv;
///
alias symv = cblas_ssymv;
///
alias symv = cblas_dsymv;
///
alias hemv = cblas_chemv;
///
alias hemv = cblas_zhemv;
///
alias spmv = cblas_sspmv;
///
alias spmv = cblas_dspmv;
///
alias spr = cblas_sspr;
///
alias spr = cblas_dspr;
///
alias hpr = cblas_chpr;
///
alias hpr = cblas_zhpr;
///
alias spr2 = cblas_sspr2;
///
alias spr2 = cblas_dspr2;
///
alias hpr2 = cblas_chpr2;
///
alias hpr2 = cblas_zhpr2;
///
alias hbmv = cblas_chbmv;
///
alias hbmv = cblas_zhbmv;
///
alias hpmv = cblas_chpmv;
///
alias hpmv = cblas_zhpmv;
///
alias gemm = cblas_sgemm;
///
alias gemm = cblas_dgemm;
///
alias gemm = cblas_cgemm;
///
alias gemm = cblas_zgemm;
///
alias symm = cblas_ssymm;
///
alias symm = cblas_dsymm;
///
alias symm = cblas_csymm;
///
alias symm = cblas_zsymm;
///
alias syrk = cblas_ssyrk;
///
alias syrk = cblas_dsyrk;
///
alias syrk = cblas_csyrk;
///
alias syrk = cblas_zsyrk;
///
alias syr2k = cblas_ssyr2k;
///
alias syr2k = cblas_dsyr2k;
///
alias syr2k = cblas_csyr2k;
///
alias syr2k = cblas_zsyr2k;
///
alias trmm = cblas_strmm;
///
alias trmm = cblas_dtrmm;
///
alias trmm = cblas_ctrmm;
///
alias trmm = cblas_ztrmm;
///
alias trsm = cblas_strsm;
///
alias trsm = cblas_dtrsm;
///
alias trsm = cblas_ctrsm;
///
alias trsm = cblas_ztrsm;
///
alias hemm = cblas_chemm;
///
alias hemm = cblas_zhemm;
///
alias herk = cblas_cherk;
///
alias herk = cblas_zherk;
///
alias her2k = cblas_cher2k;
///
alias her2k = cblas_zher2k;
///
alias xerbla = cblas_xerbla;

version(CBLAS_STD_COMPLEX)
{
	import std.complex : Complex;
}
else
{
	import mir.complex : Complex;
}

alias _cfloat = Complex!float;
alias _cdouble = Complex!double;

extern(C) @system @nogc nothrow {

  version(OPENBLAS)
  {
    enum OPENBLAS_SEQUENTIAL = 0;
    /* OpenBLAS is compiled using normal threading model */
    enum OPENBLAS_THREAD = 1;
    /* OpenBLAS is compiled using OpenMP threading model */
    enum OPENBLAS_OPENMP = 2;
  }

  ///
  alias CBLAS_INDEX = size_t;

  version(BLASNATIVEINT)
    ///
    alias blasint = ptrdiff_t;
  else
    ///
    alias blasint = int;

  ///
  enum CBLAS_LAYOUT {
    ///
    RowMajor=101,
    ///
    ColMajor=102
  };

  ///
  alias Layout = CBLAS_LAYOUT;

  alias CBLAS_ORDER = CBLAS_LAYOUT;

  alias Order = CBLAS_ORDER;

  ///
  enum CBLAS_TRANSPOSE {
    ///
    NoTrans=111,
    ///
    Trans=112,
    ///
    ConjTrans=113,
    ///
    ConjNoTrans=114
  };

  ///
  alias Transpose = CBLAS_TRANSPOSE;

  ///
  enum CBLAS_UPLO {
    ///
    Upper=121,
    ///
    Lower=122
  };

  ///
  alias Uplo = CBLAS_UPLO;

  ///
  enum CBLAS_DIAG {
    ///
    NonUnit=131,
    ///
    Unit=132
  };

  ///
  alias Diag = CBLAS_DIAG;

  ///
  enum CBLAS_SIDE {
    ///
    Left=141,
    ///
    Right=142
  };

  ///
  alias Side = CBLAS_SIDE;

  ///
  enum CBLAS_OFFSET {
    ///
    RowOffset=171,
    ///
    ColOffset=172,
    ///
    FixOffset=173
  };

  ///
  alias Offset = CBLAS_OFFSET;

  version(OPENBLAS)
  {
    /// Set the number of threads on runtime.
    void set_num_threads(int num_threads);
    ///ditto
    void goto_set_num_threads(int num_threads);

    /// Get the build configure on runtime.
    char* get_config();

    /// Get the parallelization type which is used by OpenBLAS
    int get_parallel();
  }

  ///
  void cblas_xerbla(blasint p, char *rout, char *form, ...);

  pure:

  ///
  float  cblas_sdsdot(blasint n, float alpha, scope const(float)* x, blasint incx, scope const(float)* y, blasint incy);
  ///
  double cblas_dsdot (blasint n, scope const(float)* x, blasint incx, scope const(float)* y, blasint incy);
  ///
  float  cblas_sdot(blasint n, scope const(float)*  x, blasint incx, scope const(float)*  y, blasint incy);
  ///
  double cblas_ddot(blasint n, scope const(double)* x, blasint incx, scope const(double)* y, blasint incy);

  ///
  _cfloat  cblas_cdotu(blasint n, scope const(_cfloat)*  x, blasint incx, scope const(_cfloat)*  y, blasint incy);
  ///
  _cfloat  cblas_cdotc(blasint n, scope const(_cfloat)*  x, blasint incx, scope const(_cfloat)*  y, blasint incy);
  ///
  _cdouble cblas_zdotu(blasint n, scope const(_cdouble)* x, blasint incx, scope const(_cdouble)* y, blasint incy);
  ///
  _cdouble cblas_zdotc(blasint n, scope const(_cdouble)* x, blasint incx, scope const(_cdouble)* y, blasint incy);

  ///
  void  cblas_cdotu_sub(blasint n, scope const(_cfloat)*  x, blasint incx, scope const(_cfloat)*  y, blasint incy, _cfloat  *ret);
  ///
  void  cblas_cdotc_sub(blasint n, scope const(_cfloat)*  x, blasint incx, scope const(_cfloat)*  y, blasint incy, _cfloat  *ret);
  ///
  void  cblas_zdotu_sub(blasint n, scope const(_cdouble)* x, blasint incx, scope const(_cdouble)* y, blasint incy, _cdouble *ret);
  ///
  void  cblas_zdotc_sub(blasint n, scope const(_cdouble)* x, blasint incx, scope const(_cdouble)* y, blasint incy, _cdouble *ret);

  ///
  float  cblas_sasum (blasint n, scope const(float)*  x, blasint incx);
  ///
  double cblas_dasum (blasint n, scope const(double)* x, blasint incx);
  ///
  float  cblas_scasum(blasint n, scope const(_cfloat)*  x, blasint incx);
  ///
  double cblas_dzasum(blasint n, scope const(_cdouble)* x, blasint incx);

  ///
  float  cblas_snrm2 (blasint n, scope const(float)*  X, blasint incX);
  ///
  double cblas_dnrm2 (blasint n, scope const(double)* X, blasint incX);
  ///
  float  cblas_scnrm2(blasint n, scope const(_cfloat)*  X, blasint incX);
  ///
  double cblas_dznrm2(blasint n, scope const(_cdouble)* X, blasint incX);

  ///
  CBLAS_INDEX cblas_isamax(blasint n, scope const(float)*  x, blasint incx);
  ///
  CBLAS_INDEX cblas_idamax(blasint n, scope const(double)* x, blasint incx);
  ///
  CBLAS_INDEX cblas_icamax(blasint n, scope const(_cfloat)*  x, blasint incx);
  ///
  CBLAS_INDEX cblas_izamax(blasint n, scope const(_cdouble)* x, blasint incx);

  ///
  void cblas_saxpy(blasint n, float alpha, scope const(float)* x, blasint incx, float *y, blasint incy);
  ///
  void cblas_daxpy(blasint n, double alpha, scope const(double)* x, blasint incx, double *y, blasint incy);
  ///
  void cblas_caxpy(blasint n, scope const ref _cfloat alpha, scope const(_cfloat)* x, blasint incx, _cfloat *y, blasint incy);
  ///
  void cblas_zaxpy(blasint n, scope const ref _cdouble alpha, scope const(_cdouble)* x, blasint incx, _cdouble *y, blasint incy);

  ///
  void cblas_scopy(blasint n, scope const(float)* x, blasint incx, float *y, blasint incy);
  ///
  void cblas_dcopy(blasint n, scope const(double)* x, blasint incx, double *y, blasint incy);
  ///
  void cblas_ccopy(blasint n, scope const(_cfloat)* x, blasint incx, _cfloat *y, blasint incy);
  ///
  void cblas_zcopy(blasint n, scope const(_cdouble)* x, blasint incx, _cdouble *y, blasint incy);

  ///
  void cblas_sswap(blasint n, float *x, blasint incx, float *y, blasint incy);
  ///
  void cblas_dswap(blasint n, double *x, blasint incx, double *y, blasint incy);
  ///
  void cblas_cswap(blasint n, _cfloat *x, blasint incx, _cfloat *y, blasint incy);
  ///
  void cblas_zswap(blasint n, _cdouble *x, blasint incx, _cdouble *y, blasint incy);

  ///
  void cblas_srot(blasint n, float *X, blasint incX, float *Y, blasint incY, float c, float s);
  ///
  void cblas_drot(blasint n, double *X, blasint incX, double *Y, blasint incY, double c, double  s);

  void cblas_srotg(float *a, float *b, float *c, float *s);
  void cblas_drotg(double *a, double *b, double *c, double *s);

  ///
  void cblas_srotm(blasint n, float *X, blasint incX, float *Y, blasint incY, scope const(float)* P);
  ///
  void cblas_drotm(blasint n, double *X, blasint incX, double *Y, blasint incY, scope const(double)* P);

  void cblas_srotmg(float *d1, float *d2, float *b1, float b2, float *P);
  void cblas_drotmg(double *d1, double *d2, double *b1, double b2, double *P);

  ///
  void cblas_sscal(blasint n, float alpha, float *X, blasint incX);
  ///
  void cblas_dscal(blasint n, double alpha, double *X, blasint incX);
  ///
  void cblas_cscal(blasint n, scope const ref _cfloat alpha, _cfloat *X, blasint incX);
  ///
  void cblas_zscal(blasint n, scope const ref _cdouble alpha, _cdouble *X, blasint incX);
  ///
  void cblas_csscal(blasint n, _cfloat alpha, _cfloat *X, blasint incX);
  ///
  void cblas_zdscal(blasint n, _cdouble alpha, _cdouble *X, blasint incX);

  ///
  void cblas_sgemv(CBLAS_ORDER order,  CBLAS_TRANSPOSE trans,  blasint m, blasint n, float alpha, scope const(float)*  a, blasint lda,  scope const(float)*  x, blasint incx,  float beta,  float  *y, blasint incy);
  ///
  void cblas_dgemv(CBLAS_ORDER order,  CBLAS_TRANSPOSE trans,  blasint m, blasint n, double alpha, scope const(double)*  a, blasint lda,  scope const(double)*  x, blasint incx,  double beta,  double  *y, blasint incy);
  ///
  void cblas_cgemv(CBLAS_ORDER order,  CBLAS_TRANSPOSE trans,  blasint m, blasint n, scope const ref _cfloat alpha, scope const(_cfloat)*  a, blasint lda,  scope const(_cfloat)*  x, blasint incx,  scope const ref _cfloat beta,  _cfloat  *y, blasint incy);
  ///
  void cblas_zgemv(CBLAS_ORDER order,  CBLAS_TRANSPOSE trans,  blasint m, blasint n, scope const ref _cdouble alpha, scope const(_cdouble)*  a, blasint lda,  scope const(_cdouble)*  x, blasint incx,  scope const ref _cdouble beta,  _cdouble  *y, blasint incy);

  ///
  void cblas_sger (CBLAS_ORDER order, blasint M, blasint n, float   alpha, scope const(float)*  X, blasint incX, scope const(float)*  Y, blasint incY, float  *A, blasint lda);
  ///
  void cblas_dger (CBLAS_ORDER order, blasint M, blasint n, double  alpha, scope const(double)* X, blasint incX, scope const(double)* Y, blasint incY, double *A, blasint lda);
  ///
  void cblas_cgeru(CBLAS_ORDER order, blasint M, blasint n, scope const ref _cfloat  alpha, scope const(_cfloat)*  X, blasint incX, scope const(_cfloat)*  Y, blasint incY, _cfloat  *A, blasint lda);
  ///
  void cblas_cgerc(CBLAS_ORDER order, blasint M, blasint n, scope const ref _cfloat  alpha, scope const(_cfloat)*  X, blasint incX, scope const(_cfloat)*  Y, blasint incY, _cfloat  *A, blasint lda);
  ///
  void cblas_zgeru(CBLAS_ORDER order, blasint M, blasint n, scope const ref _cdouble alpha, scope const(_cdouble)* X, blasint incX, scope const(_cdouble)* Y, blasint incY, _cdouble *A, blasint lda);
  ///
  void cblas_zgerc(CBLAS_ORDER order, blasint M, blasint n, scope const ref _cdouble alpha, scope const(_cdouble)* X, blasint incX, scope const(_cdouble)* Y, blasint incY, _cdouble *A, blasint lda);

  ///
  void cblas_strsv(CBLAS_ORDER order, CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag, blasint n, scope const(float)* A, blasint lda, float *X, blasint incX);
  ///
  void cblas_dtrsv(CBLAS_ORDER order, CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag, blasint n, scope const(double)* A, blasint lda, double *X, blasint incX);
  ///
  void cblas_ctrsv(CBLAS_ORDER order, CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag, blasint n, scope const(_cfloat)* A, blasint lda, _cfloat *X, blasint incX);
  ///
  void cblas_ztrsv(CBLAS_ORDER order, CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag, blasint n, scope const(_cdouble)* A, blasint lda, _cdouble *X, blasint incX);

  ///
  void cblas_strmv(CBLAS_ORDER order, CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag, blasint n, scope const(float)* A, blasint lda, float *X, blasint incX);
  ///
  void cblas_dtrmv(CBLAS_ORDER order, CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag, blasint n, scope const(double)* A, blasint lda, double *X, blasint incX);
  ///
  void cblas_ctrmv(CBLAS_ORDER order, CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag, blasint n, scope const(_cfloat)* A, blasint lda, _cfloat *X, blasint incX);
  ///
  void cblas_ztrmv(CBLAS_ORDER order, CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag, blasint n, scope const(_cdouble)* A, blasint lda, _cdouble *X, blasint incX);

  ///
  void cblas_ssyr(CBLAS_ORDER order, CBLAS_UPLO Uplo, blasint n, float alpha, scope const(float)* X, blasint incX, float *A, blasint lda);
  ///
  void cblas_dsyr(CBLAS_ORDER order, CBLAS_UPLO Uplo, blasint n, double alpha, scope const(double)* X, blasint incX, double *A, blasint lda);
  ///
  void cblas_cher(CBLAS_ORDER order, CBLAS_UPLO Uplo, blasint n, _cfloat alpha, scope const(_cfloat)* X, blasint incX, _cfloat *A, blasint lda);
  ///
  void cblas_zher(CBLAS_ORDER order, CBLAS_UPLO Uplo, blasint n, _cdouble alpha, scope const(_cdouble)* X, blasint incX, _cdouble *A, blasint lda);

  ///
  void cblas_ssyr2(CBLAS_ORDER order, CBLAS_UPLO Uplo, blasint n, float alpha, scope const(float)* X, blasint incX, scope const(float)* Y, blasint incY, float *A, blasint lda);
  ///
  void cblas_dsyr2(CBLAS_ORDER order, CBLAS_UPLO Uplo, blasint n, double alpha, scope const(double)* X, blasint incX, scope const(double)* Y, blasint incY, double *A, blasint lda);
  ///
  void cblas_cher2(CBLAS_ORDER order, CBLAS_UPLO Uplo, blasint n, scope const ref _cfloat alpha, scope const(_cfloat)* X, blasint incX, scope const(_cfloat)* Y, blasint incY, _cfloat *A, blasint lda);
  ///
  void cblas_zher2(CBLAS_ORDER order, CBLAS_UPLO Uplo, blasint n, scope const ref _cdouble alpha, scope const(_cdouble)* X, blasint incX, scope const(_cdouble)* Y, blasint incY, _cdouble *A, blasint lda);

  ///
  void cblas_sgbmv(CBLAS_ORDER order, CBLAS_TRANSPOSE TransA, blasint M, blasint n, blasint KL, blasint KU, float alpha, scope const(float)* A, blasint lda, scope const(float)* X, blasint incX, float beta, float *Y, blasint incY);
  ///
  void cblas_dgbmv(CBLAS_ORDER order, CBLAS_TRANSPOSE TransA, blasint M, blasint n, blasint KL, blasint KU, double alpha, scope const(double)* A, blasint lda, scope const(double)* X, blasint incX, double beta, double *Y, blasint incY);
  ///
  void cblas_cgbmv(CBLAS_ORDER order, CBLAS_TRANSPOSE TransA, blasint M, blasint n, blasint KL, blasint KU, scope const ref _cfloat alpha, scope const(_cfloat)* A, blasint lda, scope const(_cfloat)* X, blasint incX, scope const ref _cfloat beta, _cfloat *Y, blasint incY);
  ///
  void cblas_zgbmv(CBLAS_ORDER order, CBLAS_TRANSPOSE TransA, blasint M, blasint n, blasint KL, blasint KU, scope const ref _cdouble alpha, scope const(_cdouble)* A, blasint lda, scope const(_cdouble)* X, blasint incX, scope const ref _cdouble beta, _cdouble *Y, blasint incY);

  ///
  void cblas_ssbmv(CBLAS_ORDER order, CBLAS_UPLO Uplo, blasint n, blasint K, float alpha, scope const(float)* A, blasint lda, scope const(float)* X, blasint incX, float beta, float *Y, blasint incY);
  ///
  void cblas_dsbmv(CBLAS_ORDER order, CBLAS_UPLO Uplo, blasint n, blasint K, double alpha, scope const(double)* A, blasint lda, scope const(double)* X, blasint incX, double beta, double *Y, blasint incY);


  ///
  void cblas_stbmv(CBLAS_ORDER order, CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag, blasint n, blasint K, scope const(float)* A, blasint lda, float *X, blasint incX);
  ///
  void cblas_dtbmv(CBLAS_ORDER order, CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag, blasint n, blasint K, scope const(double)* A, blasint lda, double *X, blasint incX);
  ///
  void cblas_ctbmv(CBLAS_ORDER order, CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag, blasint n, blasint K, scope const(_cfloat)* A, blasint lda, _cfloat *X, blasint incX);
  ///
  void cblas_ztbmv(CBLAS_ORDER order, CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag, blasint n, blasint K, scope const(_cdouble)* A, blasint lda, _cdouble *X, blasint incX);

  ///
  void cblas_stbsv(CBLAS_ORDER order, CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag, blasint n, blasint K, scope const(float)* A, blasint lda, float *X, blasint incX);
  ///
  void cblas_dtbsv(CBLAS_ORDER order, CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag, blasint n, blasint K, scope const(double)* A, blasint lda, double *X, blasint incX);
  ///
  void cblas_ctbsv(CBLAS_ORDER order, CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag, blasint n, blasint K, scope const(_cfloat)* A, blasint lda, _cfloat *X, blasint incX);
  ///
  void cblas_ztbsv(CBLAS_ORDER order, CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag, blasint n, blasint K, scope const(_cdouble)* A, blasint lda, _cdouble *X, blasint incX);

  ///
  void cblas_stpmv(CBLAS_ORDER order, CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag, blasint n, scope const(float)* Ap, float *X, blasint incX);
  ///
  void cblas_dtpmv(CBLAS_ORDER order, CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag, blasint n, scope const(double)* Ap, double *X, blasint incX);
  ///
  void cblas_ctpmv(CBLAS_ORDER order, CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag, blasint n, scope const(_cfloat)* Ap, _cfloat *X, blasint incX);
  ///
  void cblas_ztpmv(CBLAS_ORDER order, CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag, blasint n, scope const(_cdouble)* Ap, _cdouble *X, blasint incX);

  ///
  void cblas_stpsv(CBLAS_ORDER order, CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag, blasint n, scope const(float)* Ap, float *X, blasint incX);
  ///
  void cblas_dtpsv(CBLAS_ORDER order, CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag, blasint n, scope const(double)* Ap, double *X, blasint incX);
  ///
  void cblas_ctpsv(CBLAS_ORDER order, CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag, blasint n, scope const(_cfloat)* Ap, _cfloat *X, blasint incX);
  ///
  void cblas_ztpsv(CBLAS_ORDER order, CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag, blasint n, scope const(_cdouble)* Ap, _cdouble *X, blasint incX);

  ///
  void cblas_ssymv(CBLAS_ORDER order, CBLAS_UPLO Uplo, blasint n, float alpha, scope const(float)* A, blasint lda, scope const(float)* X, blasint incX, float beta, float *Y, blasint incY);
  ///
  void cblas_dsymv(CBLAS_ORDER order, CBLAS_UPLO Uplo, blasint n, double alpha, scope const(double)* A, blasint lda, scope const(double)* X, blasint incX, double beta, double *Y, blasint incY);
  ///
  void cblas_chemv(CBLAS_ORDER order, CBLAS_UPLO Uplo, blasint n, scope const ref _cfloat alpha, scope const(_cfloat)* A, blasint lda, scope const(_cfloat)* X, blasint incX, scope const ref _cfloat beta, _cfloat *Y, blasint incY);
  ///
  void cblas_zhemv(CBLAS_ORDER order, CBLAS_UPLO Uplo, blasint n, scope const ref _cdouble alpha, scope const(_cdouble)* A, blasint lda, scope const(_cdouble)* X, blasint incX, scope const ref _cdouble beta, _cdouble *Y, blasint incY);


  ///
  void cblas_sspmv(CBLAS_ORDER order, CBLAS_UPLO Uplo, blasint n, float alpha, scope const(float)* Ap, scope const(float)* X, blasint incX, float beta, float *Y, blasint incY);
  ///
  void cblas_dspmv(CBLAS_ORDER order, CBLAS_UPLO Uplo, blasint n, double alpha, scope const(double)* Ap, scope const(double)* X, blasint incX, double beta, double *Y, blasint incY);

  ///
  void cblas_sspr(CBLAS_ORDER order, CBLAS_UPLO Uplo, blasint n, float alpha, scope const(float)* X, blasint incX, float *Ap);
  ///
  void cblas_dspr(CBLAS_ORDER order, CBLAS_UPLO Uplo, blasint n, double alpha, scope const(double)* X, blasint incX, double *Ap);

  ///
  void cblas_chpr(CBLAS_ORDER order, CBLAS_UPLO Uplo, blasint n, float alpha, scope const(float)* X, blasint incX, float *A);
  ///
  void cblas_zhpr(CBLAS_ORDER order, CBLAS_UPLO Uplo, blasint n, double alpha, scope const(double)* X, blasint incX, double *A);

  ///
  void cblas_sspr2(CBLAS_ORDER order, CBLAS_UPLO Uplo, blasint n, float alpha, scope const(float)* X, blasint incX, scope const(float)* Y, blasint incY, float *A);
  ///
  void cblas_dspr2(CBLAS_ORDER order, CBLAS_UPLO Uplo, blasint n, double alpha, scope const(double)* X, blasint incX, scope const(double)* Y, blasint incY, double *A);
  ///
  void cblas_chpr2(CBLAS_ORDER order, CBLAS_UPLO Uplo, blasint n, scope const ref _cfloat alpha, scope const(_cfloat)* X, blasint incX, scope const(_cfloat)* Y, blasint incY, _cfloat *Ap);
  ///
  void cblas_zhpr2(CBLAS_ORDER order, CBLAS_UPLO Uplo, blasint n, scope const ref _cdouble alpha, scope const(_cdouble)* X, blasint incX, scope const(_cdouble)* Y, blasint incY, _cdouble *Ap);

  ///
  void cblas_chbmv(CBLAS_ORDER order, CBLAS_UPLO Uplo, blasint n, blasint K, scope const ref _cfloat alpha, scope const(float)* A, blasint lda, scope const(float)* X, blasint incX, _cfloat beta, float *Y, blasint incY);
  ///
  void cblas_zhbmv(CBLAS_ORDER order, CBLAS_UPLO Uplo, blasint n, blasint K, scope const ref _cdouble alpha, scope const(double)* A, blasint lda, scope const(double)* X, blasint incX, _cdouble beta, double *Y, blasint incY);

  ///
  void cblas_chpmv(CBLAS_ORDER order, CBLAS_UPLO Uplo, blasint n, scope const ref _cfloat alpha, scope const(float)* Ap, scope const(float)* X, blasint incX, _cfloat beta, float *Y, blasint incY);
  ///
  void cblas_zhpmv(CBLAS_ORDER order, CBLAS_UPLO Uplo, blasint n, scope const ref _cdouble alpha, scope const(double)* Ap, scope const(double)* X, blasint incX, _cdouble beta, double *Y, blasint incY);

  ///
  void cblas_sgemm(CBLAS_ORDER order, CBLAS_TRANSPOSE TransA, CBLAS_TRANSPOSE TransB, blasint M, blasint n, blasint K, float alpha, scope const(float)* A, blasint lda, scope const(float)* B, blasint ldb, float beta, float *C, blasint ldc);
  ///
  void cblas_dgemm(CBLAS_ORDER order, CBLAS_TRANSPOSE TransA, CBLAS_TRANSPOSE TransB, blasint M, blasint n, blasint K, double alpha, scope const(double)* A, blasint lda, scope const(double)* B, blasint ldb, double beta, double *C, blasint ldc);
  ///
  void cblas_cgemm(CBLAS_ORDER order, CBLAS_TRANSPOSE TransA, CBLAS_TRANSPOSE TransB, blasint M, blasint n, blasint K, scope const ref _cfloat alpha, scope const(_cfloat)* A, blasint lda, scope const(_cfloat)* B, blasint ldb, scope const ref _cfloat beta, _cfloat *C, blasint ldc);
  ///
  void cblas_zgemm(CBLAS_ORDER order, CBLAS_TRANSPOSE TransA, CBLAS_TRANSPOSE TransB, blasint M, blasint n, blasint K, scope const ref _cdouble alpha, scope const(_cdouble)* A, blasint lda, scope const(_cdouble)* B, blasint ldb, scope const ref _cdouble beta, _cdouble *C, blasint ldc);

  ///
  void cblas_ssymm(CBLAS_ORDER order, CBLAS_SIDE Side, CBLAS_UPLO Uplo, blasint M, blasint n, float alpha, scope const(float)* A, blasint lda, scope const(float)* B, blasint ldb, float beta, float *C, blasint ldc);
  ///
  void cblas_dsymm(CBLAS_ORDER order, CBLAS_SIDE Side, CBLAS_UPLO Uplo, blasint M, blasint n, double alpha, scope const(double)* A, blasint lda, scope const(double)* B, blasint ldb, double beta, double *C, blasint ldc);
  ///
  void cblas_csymm(CBLAS_ORDER order, CBLAS_SIDE Side, CBLAS_UPLO Uplo, blasint M, blasint n, scope const ref _cfloat alpha, scope const(_cfloat)* A, blasint lda, scope const(_cfloat)* B, blasint ldb, scope const ref _cfloat beta, _cfloat *C, blasint ldc);
  ///
  void cblas_zsymm(CBLAS_ORDER order, CBLAS_SIDE Side, CBLAS_UPLO Uplo, blasint M, blasint n, scope const ref _cdouble alpha, scope const(_cdouble)* A, blasint lda, scope const(_cdouble)* B, blasint ldb, scope const ref _cdouble beta, _cdouble *C, blasint ldc);

  ///
  void cblas_ssyrk(CBLAS_ORDER order, CBLAS_UPLO Uplo, CBLAS_TRANSPOSE Trans, blasint n, blasint K, float alpha, scope const(float)* A, blasint lda, float beta, float *C, blasint ldc);
  ///
  void cblas_dsyrk(CBLAS_ORDER order, CBLAS_UPLO Uplo, CBLAS_TRANSPOSE Trans, blasint n, blasint K, double alpha, scope const(double)* A, blasint lda, double beta, double *C, blasint ldc);
  ///
  void cblas_csyrk(CBLAS_ORDER order, CBLAS_UPLO Uplo, CBLAS_TRANSPOSE Trans, blasint n, blasint K, scope const ref _cfloat alpha, scope const(_cfloat)* A, blasint lda, scope const ref _cfloat beta, _cfloat *C, blasint ldc);
  ///
  void cblas_zsyrk(CBLAS_ORDER order, CBLAS_UPLO Uplo, CBLAS_TRANSPOSE Trans, blasint n, blasint K, scope const ref _cdouble alpha, scope const(_cdouble)* A, blasint lda, scope const ref _cdouble beta, _cdouble *C, blasint ldc);

  ///
  void cblas_ssyr2k(CBLAS_ORDER order, CBLAS_UPLO Uplo, CBLAS_TRANSPOSE Trans, blasint n, blasint K, float alpha, scope const(float)* A, blasint lda, scope const(float)* B, blasint ldb, float beta, float *C, blasint ldc);
  ///
  void cblas_dsyr2k(CBLAS_ORDER order, CBLAS_UPLO Uplo, CBLAS_TRANSPOSE Trans, blasint n, blasint K, double alpha, scope const(double)* A, blasint lda, scope const(double)* B, blasint ldb, double beta, double *C, blasint ldc);
  ///
  void cblas_csyr2k(CBLAS_ORDER order, CBLAS_UPLO Uplo, CBLAS_TRANSPOSE Trans, blasint n, blasint K, scope const ref _cfloat alpha, scope const(_cfloat)* A, blasint lda, scope const(_cfloat)* B, blasint ldb, scope const ref _cfloat beta, _cfloat *C, blasint ldc);
  ///
  void cblas_zsyr2k(CBLAS_ORDER order, CBLAS_UPLO Uplo, CBLAS_TRANSPOSE Trans, blasint n, blasint K, scope const ref _cdouble alpha, scope const(_cdouble)* A, blasint lda, scope const(_cdouble)* B, blasint ldb, scope const ref _cdouble beta, _cdouble *C, blasint ldc);

  ///
  void cblas_strmm(CBLAS_ORDER order, CBLAS_SIDE Side, CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag, blasint M, blasint n, float alpha, scope const(float)* A, blasint lda, float *B, blasint ldb);
  ///
  void cblas_dtrmm(CBLAS_ORDER order, CBLAS_SIDE Side, CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag, blasint M, blasint n, double alpha, scope const(double)* A, blasint lda, double *B, blasint ldb);
  ///
  void cblas_ctrmm(CBLAS_ORDER order, CBLAS_SIDE Side, CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag, blasint M, blasint n, scope const ref _cfloat alpha, scope const(_cfloat)* A, blasint lda, _cfloat *B, blasint ldb);
  ///
  void cblas_ztrmm(CBLAS_ORDER order, CBLAS_SIDE Side, CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag, blasint M, blasint n, scope const ref _cdouble alpha, scope const(_cdouble)* A, blasint lda, _cdouble *B, blasint ldb);

  ///
  void cblas_strsm(CBLAS_ORDER order, CBLAS_SIDE Side, CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag, blasint M, blasint n, float alpha, scope const(float)* A, blasint lda, float *B, blasint ldb);
  ///
  void cblas_dtrsm(CBLAS_ORDER order, CBLAS_SIDE Side, CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag, blasint M, blasint n, double alpha, scope const(double)* A, blasint lda, double *B, blasint ldb);
  ///
  void cblas_ctrsm(CBLAS_ORDER order, CBLAS_SIDE Side, CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag, blasint M, blasint n, scope const ref _cfloat alpha, scope const(_cfloat)* A, blasint lda, _cfloat *B, blasint ldb);
  ///
  void cblas_ztrsm(CBLAS_ORDER order, CBLAS_SIDE Side, CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag, blasint M, blasint n, scope const ref _cdouble alpha, scope const(_cdouble)* A, blasint lda, _cdouble *B, blasint ldb);

  ///
  void cblas_chemm(CBLAS_ORDER order, CBLAS_SIDE Side, CBLAS_UPLO Uplo, blasint M, blasint n, scope const ref _cfloat alpha, scope const(float)* A, blasint lda, scope const(float)* B, blasint ldb, _cfloat beta, float *C, blasint ldc);
  ///
  void cblas_zhemm(CBLAS_ORDER order, CBLAS_SIDE Side, CBLAS_UPLO Uplo, blasint M, blasint n, scope const ref _cdouble alpha, scope const(double)* A, blasint lda, scope const(double)* B, blasint ldb, _cdouble beta, double *C, blasint ldc);

  ///
  void cblas_cherk(CBLAS_ORDER order, CBLAS_UPLO Uplo, CBLAS_TRANSPOSE Trans, blasint n, blasint K, float alpha, scope const(float)* A, blasint lda, float beta, float *C, blasint ldc);
  ///
  void cblas_zherk(CBLAS_ORDER order, CBLAS_UPLO Uplo, CBLAS_TRANSPOSE Trans, blasint n, blasint K, double alpha, scope const(double)* A, blasint lda, double beta, double *C, blasint ldc);

  ///
  void cblas_cher2k(CBLAS_ORDER order, CBLAS_UPLO Uplo, CBLAS_TRANSPOSE Trans, blasint n, blasint K, scope const ref _cfloat alpha, scope const(float)* A, blasint lda, scope const(float)* B, blasint ldb, float beta, float *C, blasint ldc);
  ///
  void cblas_zher2k(CBLAS_ORDER order, CBLAS_UPLO Uplo, CBLAS_TRANSPOSE Trans, blasint n, blasint K, scope const ref _cdouble alpha, scope const(double)* A, blasint lda, scope const(double)* B, blasint ldb, double beta, double *C, blasint ldc);
}
