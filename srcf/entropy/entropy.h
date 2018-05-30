/*
 * entropy.h
 *
 *  Created on: 30 May 2018
 *      Author: Omar
 */

#ifndef ENTROPY_ENTROPY_H_
#define ENTROPY_ENTROPY_H_

double entropy(const int *d, int nsamples, int nvars, int c, bool *v);
double entropyFast(const int *d, int nsamples, int nvars, int c,int *vnew);



#endif /* ENTROPY_ENTROPY_H_ */
