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
CSRMat;

void CSRMatTrans(CSRMat *ai, CSRMat *ao);

CSRMat *CSRMatCreateEmpty(int n, int nnz);

CSRMat *CSRMatCreate(const char *filename);

CSRMat *CSRMatCreateMM(const char *filename);

CSRMat *CSRMatRead(const char *filename);

void CSRMatDestroy(CSRMat *mat);

void CSRMatDump(const CSRMat *mat);

void CSRMatMultVec(const CSRMat *mat, const double *x, double *y);

//void CSRMatMultVec_instru(const CSRMat *mat, const double *x, double *y);

//void CSRMatMultVec_MKL(const CSRMat *mat, const double *x, double *y);

void CSRMatMultVecBlock(const CSRMat *mat, const double *x, double *y, int
		blk_sz);

#endif
