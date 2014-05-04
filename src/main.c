#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <assert.h>
#include "global_vars.h"
#include "arch.h"
#include "csr_mat.h"
#include <omp.h>
#include "wallclock.h"
#include "rand_gen.h"
#include "blas_headers.h"
#include "full_matrices.h"
#include "rand_range_finder.h"
#include "rand_svd.h"

timer_struct rand_timer;

int numthreads;


void print_usage(char **argv)
{
	printf("Usage: %s matrix-file matrix-trans-file p q"
			"sing-vals-file\n", argv[0]);
}

int main(int argc, char **argv)
{
	int p;
	int q;
	csr_mat *A;
	csr_mat *At;
	range_finder *r;
	svd_finder *svd_f;
	double *svals;
	int i;
	FILE *fp;
	double time1, time2;

	if (argc < 7) {
		print_usage(argv);
		return 0;
	}
	srand(time(NULL));
	numthreads = atoi(argv[6]);
	omp_set_num_threads(numthreads);

	rand_timer.time_matvec = 0.0;
	rand_timer.time_dgemm = 0.0;
	rand_timer.time_svd = 0.0;
	rand_timer.time_qr = 0.0;
	rand_timer.time_rand = 0.0;
	rand_timer.time_total = 0.0;

	printf("Reading sparse matrix ... ");

	A = csr_mat_create_mm(argv[1]);

	assert(A!=NULL);

	printf("done\n");

	printf("Reading sparse matrix transpose ... ");

	At = csr_mat_create_mm(argv[2]);

	assert(At!=NULL);

	printf("done\n");

	p = atoi(argv[3]);
	q = atoi(argv[4]);

	printf("samples = %d\n", p);
	printf("pow. it = %d\n", q);

	r = init_range_finder(A->n, p, q);
	svd_f = init_svd_finder(A->n, p);
	svals = (double *)malloc((A->n)*sizeof(double));

	WALLCLOCK(time1);
	range_find(r, A, At);

	compute_svals_rand(svd_f, At, r, svals);
	WALLCLOCK(time2);
	rand_timer.time_total += time2 - time1;

	printf("matvec = %lf\n", rand_timer.time_matvec);
	printf("qr = %lf\n", rand_timer.time_qr);
	printf("svd = %lf\n", rand_timer.time_svd);
	printf("rand = %lf\n", rand_timer.time_rand);
	printf("total = %lf\n", rand_timer.time_total);

	fp = fopen(argv[5], "w");
	assert(fp!=NULL);

	for (i = 0; i < p; i++) {
		fprintf(fp, "%1.15e\n", svals[i]);
	}
	fclose(fp);

	free(svals);
	csr_mat_destroy(A);
	csr_mat_destroy(At);
	destroy_range_finder(r);
	destroy_svd_finder(svd_f);
	return 0;
}

