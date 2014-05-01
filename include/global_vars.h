#ifndef __GLOBAL_VARS_H__
#define __GLOBAL_VARS_H__

typedef struct {
	double time_dgemm;
	double time_matvec;
	double time_svd;
	double time_total;
	double time_rand;
	double time_qr;
} timer_struct;

#endif
