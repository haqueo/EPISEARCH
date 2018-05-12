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
double * calculateMeasures(int p1, double Hp1, int p2, double Hp2, int p3,
		int cl, const int *d,int nsamples, int nvars, int c, double Hcl, double Hp1cl);

#endif /* FULLSEARCH_H_ */
