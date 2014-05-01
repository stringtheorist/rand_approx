#include <stdio.h>
#include <stdlib.h>
#include "global_vars.h"
#include "arch.h"
#include "CSRMat.h"
#include "full_matrices.h"
#include "rand_range_finder.h"
#include "rand_svd.h"
#include "wallclock.h"

extern timer_struct rand_timer;

void compute_svals_rand(svd_finder *svd_f, CSRMat *At, range_finder *r,
		double *svals)
{
	double time1, time2;
	/*Compute B' = A'*Q*/
	WALLCLOCK(time1);
	CSRMatMultVecBlock(At, r->q_mat->val, svd_f->Bt->val, svd_f->p);
	WALLCLOCK(time2);
	rand_timer.time_matvec += (time2 - time1);

	/*Compute SVs of B*/
	WALLCLOCK(time1);
	compute_svd_vals(svd_f->Bt, svd_f->svd, svals);
	WALLCLOCK(time2);
	rand_timer.time_svd += (time2 - time1);
}

svd_finder *init_svd_finder(int n, int p)
{
	svd_finder *svd_f;

	svd_f = (svd_finder *)malloc(sizeof(svd_finder));
	svd_f->p = p;
	svd_f->Bt = init_zero_matrix(n, p);
	svd_f->svd = init_svd_scratch(svd_f->Bt);
	return svd_f;
}

void destroy_svd_finder(svd_finder *svd_f)
{
	destroy_matrix(svd_f->Bt);
	destroy_svd_scratch(svd_f->svd);
	free(svd_f);
}

