#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "blas_headers.h"
#include "full_matrices.h"
#include "rand_gen.h"

void print_matrix(FILE *fp, full_matrix *mat)
{
	int i;
	int j;
	int m;
	int n;
	m = mat->m;
	n = mat->n;

	for (i = 0; i < m; i++) {
		for (j = 0; j < n; j++) {
			fprintf(fp, "%lf,", mat->val[i + m*j]);
		}
		fprintf(fp, "\n");
	}

}

full_matrix *init_zero_matrix(int m, int n)
{
	full_matrix *mat;
	mat = (full_matrix *)malloc(sizeof(full_matrix));
	mat->m = m;
	mat->n = n;
	mat->val = (double *)calloc(m * n, sizeof(double));
	return mat;
}

void destroy_matrix(full_matrix *mat)
{
	free(mat->val);
	free(mat);
}

void symmetrize_matrix(full_matrix *mat)
{
	int i, j;

	assert(mat->m == mat->n);

	for (i = 0; i < mat->m; i++) {
		for (j = 0; j < i; j++) {
			mat->val[i + (mat->m) * j] = mat->val[j + (mat->m) * i];
		}
	}

}

void fill_rand_matrix(full_matrix *mat)
{
	int m, n;

	m = mat->m;
	n = mat->n;

	fill_rand_vec(m * n, mat->val);

}

int calculate_qr_scratch(full_matrix *mat)
{

	int lwork = -1;
	int info;
	double tau;
	double work;
	int m;
	int n;

	m = mat->m;
	n = mat->n;
	dgeqrf_(&m, &n, mat->val, &m, &tau, &work, &lwork, &info);
	return ((int) work);
}

void compute_qr_q(full_matrix *mat, qr_scratch *qr)
{
	int m;
	int n;
	int info;

	m = mat->m;
	n = mat->n;


	dgeqrf_(&m, &n, mat->val, &m, qr->tau, qr->work, &(qr->lwork), &info);
	dorgqr_(&m, &n, &n, mat->val, &m, qr->tau, qr->work, &(qr->lwork),
			&info);

}

qr_scratch *init_qr_scratch(full_matrix *mat)
{
	int ltau;
	qr_scratch *qr;

	qr = (qr_scratch *)malloc(sizeof(qr_scratch));

	qr->lwork = calculate_qr_scratch(mat) + 1;
	qr->work = (double *)malloc((qr->lwork)*sizeof(double));

	if (mat->m > mat->n) {
		ltau = (mat->m) + 1;
	} else {
		ltau = (mat->n) + 1;
	}

	qr->tau = (double *)malloc(ltau * sizeof(double));

	return qr;
}

void destroy_qr_scratch(qr_scratch *qr)
{
	free(qr->work);
	free(qr->tau);
	free(qr);
}

int compute_svd_scratch(full_matrix *mat)
{
	char jobu;
	char jobvt;
	double s;
	double work;
	int lwork;
	int info;
	double vt;
	double u;
	int ldu = 1;
	int ldvt = 1;


	jobu = 'N';
	jobvt = 'N';
	lwork = -1;

	dgesvd_(&jobu, &jobvt, &(mat->m), &(mat->n), mat->val, &(mat->m),
			&s, &u, &ldu, &vt, &ldvt, &work, &lwork, &info);
	return (int) work;
}

svd_scratch *init_svd_scratch(full_matrix *mat)
{
	svd_scratch *svd;

	svd = (svd_scratch *)malloc(sizeof(svd_scratch));
	svd->lwork = compute_svd_scratch(mat);
	svd->work = (double *)malloc((svd->lwork)*sizeof(double));

	return svd;
}

void destroy_svd_scratch(svd_scratch *svd)
{
	free(svd->work);
	free(svd);
}

void compute_svd_vals(full_matrix *mat, svd_scratch *svd, double *svals)
{
	char jobu = 'N';
	char jobvt = 'N';
	int ldu = 1;
	int ldvt = 1;
	int info;
	double u;
	double vt;

	dgesvd_(&jobu, &jobvt, &(mat->m), &(mat->n), mat->val, &(mat->m),
			svals, &u, &ldu, &vt, &ldvt, svd->work, &(svd->lwork),
			&info);

}

void multiply_matrices(full_matrix *a, char transa, full_matrix *b, char
		transb, full_matrix *out)
{
	int m;
	int n;
	int k;
	int lda;
	int ldb;
	int ldout;
	double alpha = 1.0;
	double beta = 0.0;

	if (transa == 'N') {
		m = a->m;
		k = a->n;
		lda = m;
	} else if (transa == 'T') {
		m = a->n;
		k = a->m;
		lda = k;
	}

	if (transb == 'N') {
		if (k != b->m) {
			fprintf(stdout, "ERROR: multiply_matrices: dimension"
					" mismatch with rows of b operand\n");
			exit(-1);
		}
		n = b->n;
		ldb = b->m;
	} else if (transb == 'T') {
		if (k != b->n) {
			fprintf(stdout, "ERROR: multiply_matrices: dimension"
					" mismatch with rows of b operand\n");
			exit(-1);
		}
		n = b->m;
		ldb = b->m;
	}

	ldout = out->m;

	if (((out->m) != m) || ((out->n) != n)) {
		fprintf(stdout, "ERROR: multiply_matrices: dimension"
				" mismatch with rows/columns of output\n");
		exit(-1);
	}

	dgemm_(&transa, &transb, &m, &n, &k, &alpha, a->val, &lda, b->val, &ldb,
			&beta, out->val, &ldout);


}























