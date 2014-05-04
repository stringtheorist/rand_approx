#ifndef _CSRMAT_H_
#define _CSRMAT_H_

typedef struct
{
    int n;
    int nnz; 
    mwIndex *ia;
    mwIndex *ja;
    double *a;
}
csr_mat;

void csr_mat_trans(csr_mat *ai, csr_mat *ao);

csr_mat *csr_mat_create_empty(int n, int nnz);

csr_mat *csr_mat_create(const char *filename);

csr_mat *csr_mat_create_mm(const char *filename);

csr_mat *csr_mat_read(const char *filename);

void csr_mat_destroy(csr_mat *mat);

void csr_mat_dump(const csr_mat *mat);

void csr_mat_mult_vec(const csr_mat *mat, const double *x, double *y);

//void CSRMatMultVec_instru(const CSRMat *mat, const double *x, double *y);

//void CSRMatMultVec_MKL(const CSRMat *mat, const double *x, double *y);

void csr_mat_mult_vec_block(const csr_mat *mat, const double *x, double *y, int
		blk_sz);

#endif
