#include <stdlib.h>
#include <stdio.h>
#include "wallclock.h"
#include "global_vars.h"
#include "arch.h"
#include "csr_mat.h"
#include "full_matrices.h"
#include "rand_range_finder.h"

extern timer_struct rand_timer;

void range_find(range_finder *r, csr_mat *A, csr_mat *At)
{
	int p;
	int q;
	int i;
	double time1, time2;

	p = r->p;
	q = r->q;

	//Q = A*omega

	WALLCLOCK(time1);

	csr_mat_mult_vec_block(A, r->omega_mat->val, r->q_mat->val, p);

	for (i = 0; i<q; i++) {
		//At multiply
		csr_mat_mult_vec_block(At, r->q_mat->val, r->omega_mat->val, p);

		//A multiply
		csr_mat_mult_vec_block(A, r->omega_mat->val, r->q_mat->val, p);

	}

	WALLCLOCK(time2);

	rand_timer.time_matvec += time2 - time1;

	WALLCLOCK(time1);
	compute_qr_q(r->q_mat, r->qr);
	WALLCLOCK(time2);

	rand_timer.time_qr += time2 - time1;
}

range_finder *init_range_finder(int n, int p, int q)
{
	range_finder *r;
	double time1, time2;

	r = (range_finder *) malloc(sizeof(range_finder));

	r->p = p;
	r->q = q;
	r->omega_mat = init_zero_matrix(n, p);
	r->q_mat = init_zero_matrix(n, p);

	WALLCLOCK(time1);
	fill_rand_matrix(r->omega_mat);
	WALLCLOCK(time2);

	rand_timer.time_rand += time2 - time1;
	r->qr = init_qr_scratch(r->q_mat);
	return r;
}

void destroy_range_finder(range_finder *r)
{
	destroy_matrix(r->omega_mat);
	destroy_matrix(r->q_mat);
	destroy_qr_scratch(r->qr);
	free(r);
}





