#include <stdio.h>
#include <stdlib.h>
#include <stdlib.h>
#include <stdio.h>
#include "mmio.h"
/*#include "mkl.h"*/ /* mkl_cspblas_dcsrgemv*/
#include "arch.h"
#include "csr_mat.h"
#include "wallclock.h"
#include <omp.h>

#define MAX(a,b) ((a)>(b) ? (a) : (b))

#define LINESIZE 100
#define MAXN     2000000

typedef enum symm_t {
	symmetric,
	unsymmetric,
} symm_t;

csr_mat *csr_mat_create_empty(int n, int nnz)
{
	csr_mat *mat = (csr_mat *) malloc(sizeof(csr_mat));

	mat->n = n;
	mat->ia = (mw_index *) malloc((n+1)*sizeof(mw_index));
	mat->ja = (mw_index *) malloc(nnz*sizeof(mw_index));
	mat->a =  (double *) malloc(nnz*sizeof(double));

	return mat;
}

csr_mat *csr_mat_read(const char *filename)
{
	csr_mat *mat = (csr_mat *) malloc(sizeof(csr_mat));

	FILE *fp = fopen(filename, "rb");
	if (fp == NULL)
	{
		printf("could not open file: %s\n", filename);
		return NULL;
	}

	int n, nnz;
	fread(&n,   sizeof(int), 1, fp);
	fread(&nnz, sizeof(int), 1, fp);

	mat->n = n;
	mat->ia = (mw_index *) malloc((mat->n+1)*sizeof(mw_index));
	mat->ja = (mw_index *) malloc(nnz*sizeof(mw_index));
	mat->a = (double *) malloc(nnz*sizeof(double));

	fread(mat->ia, sizeof(int), n+1, fp);
	fread(mat->ja, sizeof(int), nnz, fp);
	fread(mat->a,  sizeof(double), nnz, fp);
	fclose(fp);

	return mat;
}

/* input is 1-based
 * will allocate space for matrix, therefore this is a create function */
csr_mat *csr_mat_create(const char *filename)
{
	csr_mat *mat = (csr_mat *) malloc(sizeof(csr_mat));

	char line[LINESIZE];
	int indi, indj;
	double val;
	int nnz = 0;
	int maxi = 0, maxj = 0;

	FILE *fp = fopen(filename, "r");
	if (fp == NULL)
	{
		printf("could not open file: %s\n", filename);
		return NULL;
	}

	int *counts = (int *) calloc(MAXN, sizeof(int));

	// first pass to count number of nonzeros per row
	while (fgets(line, LINESIZE, fp) != NULL)
	{
		sscanf(line, "%d %d %lf", &indi, &indj, &val);
		maxi = MAX(maxi, indi);
		maxj = MAX(maxj, indj);
		counts[indi-1]++;
		nnz++;
	}

	if (maxi > MAXN || maxj > MAXN)
	{
		printf("matrix is too large\n");
		free(counts);
		return NULL;
	}

	mat->n = MAX(maxi, maxj);  // assume matrix is square
	mat->ia = (mw_index *) malloc((mat->n+1)*sizeof(mw_index));
	mat->ja = (mw_index *) malloc(nnz*sizeof(mw_index));
	mat->a = (double *) malloc(nnz*sizeof(double));

	// set ia pointers
	int i;
	mat->ia[0] = 0;
	for (i=0; i<mat->n; i++)
	{
		mat->ia[i+1] = mat->ia[i] + counts[i];
		counts[i] = mat->ia[i];
		// use counts array to point to empty space in ja array
	}

	// second pass
	rewind(fp);
	while (fgets(line, LINESIZE, fp) != NULL)
	{
		sscanf(line, "%d %d %lf", &indi, &indj, &val);
		int k = counts[indi-1]++;
		mat->ja[k] = indj-1;
		mat->a[k] = val;
	}

	fclose(fp);
	free(counts);

	return mat;
}

csr_mat *csr_mat_create_mm(const char *filename)
{
	symm_t symmetry_code;
	csr_mat *mat;
	int ret_code;
	MM_typecode matcode;
	FILE *f;
	int M, N, nnz, nnz_;
	int off_diag;
	int i, *I, *J;
	int k;
	double *val;
	int *counts;

	if ((f = fopen(filename, "r")) == NULL) {
		fprintf(stderr, "CSRMatCreateMM: Error: File read error.\n");
		return NULL;
	}

	if(mm_read_banner(f, &matcode) != 0) {
		fprintf(stderr, "CSRMatCreateMM: Error: Could not process MM"
				"banner.\n");
		return NULL;
	}

	if (mm_is_complex(matcode)) {
		fprintf(stderr, "CSRMatCreateMM: Error: Complex matrices not"
				"supported.\n");
		return NULL;
	}
	if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nnz_)) != 0) {
		fprintf(stderr, "CSRMatCreateMM: Error: Can't read matrix"
				"size.\n");
		return NULL;
	}

	if (mm_is_symmetric(matcode)) {
		symmetry_code = symmetric;
#if 0
		fprintf(stdout, "Reading symmetric matrix with:\n");
		fprintf(stdout, "nnz = %d\n", nnz_);
		fprintf(stdout, "n = %d\n", N);
#endif
	} else {
		symmetry_code = unsymmetric;
#if 0
		fprintf(stdout, "Reading unsymmetric matrix with:\n");
		fprintf(stdout, "nnz = %d\n", nnz_);
		fprintf(stdout, "n = %d\n", N);
		fprintf(stdout, "m = %d\n", M); 
#endif
	}

	/*Because we don't do anything special for symmetric matrices*/
	/*Irritatingly MM stores only the upper/lower triangle*/
	if (symmetry_code == symmetric) {
		nnz = (2*nnz_) - N;
	} else {
		nnz = nnz_;
	}


	mat = (csr_mat *)malloc(sizeof(csr_mat));

	I = (int *) malloc(nnz * sizeof(int));
	J = (int *) malloc(nnz * sizeof(int));
	val = (double *)malloc(nnz * sizeof(double));

	for (i = 0; i < nnz_; i++) {
		fscanf(f, "%d %d %lg\n", &I[i], &J[i], &val[i]);
		I[i]--;
		J[i]--;
	}

	fclose(f);

	/*More annoying nonsense to be done because only unique entries
	 * are stored by MM for symmetric matrices */

	if (symmetry_code == symmetric) {
		off_diag = 0;
		for (i = 0; i < nnz_; i++) {
			if (I[i] != J[i]) {
				I[off_diag + nnz_] = J[i];
				J[off_diag + nnz_] = I[i];
				val[off_diag + nnz_] = val[i];
				off_diag++;
			}
		}
	}

	/*Now fill in the corresponding CSRMat structure*/
	mat->n = N;
	mat->nnz = nnz;
	mat->ia = (mw_index *) malloc(((mat->n) + 1) * sizeof(mw_index));
	mat->ja = (mw_index *) malloc(nnz * sizeof(mw_index));
	mat->a = (double *) malloc(nnz * sizeof(double));
	counts = (int *) calloc(mat->n,  sizeof(int));

	/*Set counts*/
	for (i = 0; i < nnz; i++) {
		counts[I[i]]++;
	}

	/* set ia pointers */
	mat->ia[0] = 0;
	for (i=0; i<mat->n; i++)
	{
		mat->ia[i+1] = mat->ia[i] + counts[i];
		counts[i] = mat->ia[i];
		/* use counts array to point to empty space in ja array*/
	}
	for (i = 0; i < nnz; i++) {
		k = counts[I[i]]++;
		mat->ja[k] = J[i];
		mat->a[k] = val[i];
	}

	free(I);
	free(J);
	free(counts);
	free(val);
	return mat;

}

void csr_mat_destroy(csr_mat *mat)
{
	free(mat->ia);
	free(mat->ja);
	free(mat->a);
	free(mat);
}

void csr_mat_dump(const csr_mat *mat)
{
	const mw_index *ia = mat->ia;
	const mw_index *ja = mat->ja;
	const double *a = mat->a;
	int i, j;

	for (i=0; i<mat->n; i++)
	{
		for (j=ia[i]; j<ia[i+1]; j++)
			printf("%d %d %f\n", i+1, (int)ja[j]+1, a[j]);
	}
}

void csr_mat_mult_vec(const csr_mat *mat, const double *x, double *y)
{
	const mw_index *ia = mat->ia;
	const mw_index *ja = mat->ja;
	const double *a = mat->a;
	int i, j;
	double t;

#pragma omp parallel for private(i,j,t)
	for (i=0; i<mat->n; i++)
	{
		t = 0.;
		for (j=ia[i]; j<ia[i+1]; j++)
			t += a[j]*x[ja[j]];
		y[i] = t;
	}
}

void csr_mat_mult_vec_block(const csr_mat *mat, const double *x, double *y, int
		blk_sz)
{
	int i;
	for (i = 0; i < blk_sz; i++) {
		csr_mat_mult_vec(mat, &(x[i*(mat->n)]), &(y[i*(mat->n)]));
	}
}


#if 0 
void CSRMatMultVec_instru(const CSRMat *mat, const double *x, double *y)
{
	int numthreads;

	const mw_index *ia = mat->ia;
	const mw_index *ja = mat->ja;
	const double *a = mat->a;
	int i, j;
	double t;
	int id;
	double time0[256], time1[256], time2, time3;
	numthreads = NUMTHREADS;

	WALLCLOCK(time2);
#pragma omp parallel num_threads(numthreads) private(i,j,t,id)
	{
		id = omp_get_thread_num();
		WALLCLOCK(time0[id]);

		//#pragma omp for schedule(dynamic,4096)
#pragma omp for
		for (i=0; i<mat->n; i++)
		{
			t = 0.;
			for (j=ia[i]; j<ia[i+1]; j++)
				t += a[j]*x[ja[j]];
			y[i] = t;
		}
		WALLCLOCK(time1[id]);
	}
	WALLCLOCK(time3);

	double min = 1.e100;
	double max = 0.;
	double ave = 0.;
	for (id=0; id<numthreads; id++)
	{
		double t = time1[id] - time0[id];
		if (t > max) max = t;
		if (t < min) min = t;
		ave += t/numthreads;
		// printf("thread %d  time %f  numiter %d  numthreads %d\n", id, time1[id]-time0[id], numiter, numthreads);
	}
	printf("min %f max %f ave %f\n", min, max, ave);
	printf("total time %f\n", time3-time2);
	printf("rate %f (peak is 163 GB/s)\n", (8.+4.) / 1.e9 * mat->ia[mat->n] / (time3-time2));
}
#endif 

void csr_mat_trans(csr_mat *ain, csr_mat *aout)
{
    int i, j, k, next, n;
    mw_index *ia;
    mw_index *ja;
    mw_index *iao;
    mw_index *jao;
    double *a;
    double *ao;

    n = ain->n;
    ia = ain->ia;
    ja = ain->ja;
    iao = aout->ia;
    jao = aout->ja;
    a = ain->a;
    ao = aout->a;

    /* initialize */
    for (i=0; i<=n; i++)
        iao[i] = 0;

    /* compute lengths of rows of transpose of A*/
    for (i=0; i<n; i++)
        for (k=ia[i]; k<ia[i+1]; k++)
            iao[ja[k]+1]++;

    /* compute pointers from lengths*/
    iao[0] = 0;
    for (i=0; i<n; i++)
        iao[i+1] = iao[i] + iao[i+1];

    /* now do the actual copying*/
    for (i=0; i<n; i++)
    {
        for (k=ia[i]; k<ia[i+1]; k++)
        {
            j = ja[k];
            next = iao[j];
            ao[next] = a[k];
            jao[next] = i;
            iao[j] = next+1;
        }
    }

    /* reshift iao into ia*/
    for (i=n-1; i>=0; i--)
        iao[i+1] = iao[i];
    iao[0] = 0;
}




#if 0
// y = alpha A x + beta y
void xx(double alpha, const CSRMat *mat, 
		const double *x, double beta, double *y)
{
	char transa = 'N';
	char matdesra[6] = "G  C";
	int n;

	mkl_dcsrmv(&transa, &n, &n, &alpha,
			matdescra, mat->a, mat->ia, mat->ja, // undone: need two arrays...
			x, &beta, y);

	/*
	   void mkl_dcsrmv(char *transa, MKL_INT *m, MKL_INT *k, double *alpha, 
	   char *matdescra, double *val, MKL_INT *indx, MKL_INT *pntrb, MKL_INT *pntre, 
	   double *x, double *beta, double *y);
	   */
}
#endif

#if 0
void CSRMatMultVec_MKL(const CSRMat *mat, const double *x, double *y)
{
	mkl_cspblas_dcsrgemv("N", (int *)&mat->n, 
			(double *)mat->a, (int *)mat->ia, (int *)mat->ja, (double *)x, y);
}
#endif 
