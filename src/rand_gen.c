#include <math.h>
#include <stdlib.h>
#include "global_vars.h"
#include "rand_gen.h"

extern timer_struct rand_timer;

void fill_rand_vec(int n, double *vec)
{
	int i;
	for (i = 0; i < n; i++) 
	{
		vec[i] = gauss_rand(0.0, 1.0); 
	}

}

double unirand(double a, double b)
{
	return rand() / (RAND_MAX + 1.0) * (b - a) + a;
}

/* Following implementation adapted from one available at
 * http://www.taygeta.com/random/boxmuller.html

   (c) Copyright 1994, Everett F. Carter Jr.
   Permission is granted by the author to use
   this software for any
   application provided this
   copyright
   notice
   is
   preserved.
*/
double gauss_rand(double m, double s)	/* normal random variate generator */
{				        /* mean m, standard deviation s */
	double x1, x2, w, y1;
	static double y2;
	static int use_last = 0;

	if (use_last)		        /* use value from previous call */
	{
		y1 = y2;
		use_last = 0;
	}
	else
	{
		do {
			x1 = 2.0 * unirand(0.0, 1.0) - 1.0;
			x2 = 2.0 * unirand(0.0, 1.0) - 1.0;
			w = x1 * x1 + x2 * x2;
		} while ( w >= 1.0 );

		w = sqrt( (-2.0 * log( w ) ) / w );
		y1 = x1 * w;
		y2 = x2 * w;
		use_last = 1;
	}

	return( m + y1 * s );
}
