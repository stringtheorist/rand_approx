#ifndef __RAND_RANGE_FINDER_H__
#define __RAND_RANGE_FINDER_H__

typedef struct {
	full_matrix *q_mat;
	full_matrix *omega_mat;
	int p;
	int q;
	qr_scratch *qr;
} range_finder;


void range_find(range_finder *r, CSRMat *A, CSRMat *At);
range_finder *init_range_finder(int n, int p, int q);
void destroy_range_finder(range_finder *r);
#endif
