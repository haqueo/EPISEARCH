/*
 * fullSearch.h
 *
 *  Created on: 3 May 2018
 *      Author: Omar
 */

#ifndef FULLSEARCH_H_
#define FULLSEARCH_H_


double entropy(const int *d, int nsamples, int nvars, int c, bool *v);
double entropyFast(const int *d, int nsamples, int nvars, int c, bool*v);


#endif /* FULLSEARCH_H_ */
