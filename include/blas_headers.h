#ifndef __BLAS_HEADERS_H__
#define __BLAS_HEADERS_H__

extern void dgeqrf_(int *m, int *n, double *a, int *lda, double *tau,
		double *work, int *lwork, int *info);
extern void dorgqr_(int *m, int *n, int *k, double *a, int *lda, double *tau,
		double *work, int *lwork, int *info);
extern void dsyevd_(char *jobz, char *uplo, int *n, double *a, int *lda,
		double *w, double *work, int *lwork,
		int *liwork, int *info);
extern void dgesvd_(char *jobu, char *jobvt, int *m, int *n, double *a,
		int *lda, double *s, double *u, int *ldu, double *vt,
		int *ldvt, double *work, int *lwork, int *info);
extern void dgemm_(char *transa, char *transb, int *m, int *n, int *k,
		double *alpha, double *a, int *lda, double *b, int *ldb,
		double *beta, double *c, int *ldc);

#endif
