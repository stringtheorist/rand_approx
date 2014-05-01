#ifndef __FULL_MATRICES_H__
#define __FULL_MATRICES_H__
typedef struct
{
	int m;
	int n;
	double *val; 
} full_matrix;

typedef struct
{
	double *work;
	double *tau;
	int lwork;
} qr_scratch;

typedef struct
{
	double *work;
	int lwork;
} svd_scratch;

void symmetrize_matrix(full_matrix *mat);
full_matrix *init_zero_matrix(int m, int n);
void destroy_matrix(full_matrix *mat);
void fill_rand_matrix(full_matrix *mat);
void print_matrix(FILE *fp, full_matrix *mat);
int calculate_qr_scratch(full_matrix *mat);
void destroy_qr_scratch(qr_scratch *qr);
qr_scratch *init_qr_scratch(full_matrix *mat);
void compute_qr_q(full_matrix *mat, qr_scratch *qr);
int compute_svd_scratch(full_matrix *mat);
svd_scratch *init_svd_scratch(full_matrix *mat);
void destroy_svd_scratch(svd_scratch *svd);
void compute_svd_vals(full_matrix *mat, svd_scratch *svd, double *svals);
void multiply_matrices(full_matrix *a, char transa, full_matrix *b,
		char transb, full_matrix *out);

#endif
