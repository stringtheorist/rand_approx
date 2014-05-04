#ifndef __RAND_SVD_H__
#define __RAND_SVD_H__

typedef struct {
	svd_scratch *svd;
	full_matrix *Bt;
	int p;
} svd_finder;

void compute_svals_rand(svd_finder *svd_f, csr_mat *At, range_finder *r,
		double *svals);
svd_finder *init_svd_finder(int n, int p);
void destroy_svd_finder(svd_finder *svd_f);

#endif
