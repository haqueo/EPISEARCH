/*
 * fullSearch.h
 *
 *  Created on: 3 May 2018
 *      Author: Omar
 */
#include <string>
#include <iostream>
#ifndef FULLSEARCH_H_
#define FULLSEARCH_H_


double entropy(const int *d, int nsamples, int nvars, int c, bool *v);
double entropyFast(const int *d, int nsamples, int nvars, int c,int *vnew);
double * calculateMeasures(int p1, double Hp1, int p2, double Hp2, int p3,
		int cl, const int *d,int nsamples, int nvars, int c, double Hcl, double Hp1cl);
std::vector<int> readData(std::string filename, int nrows, int nvars);
void runFullSearch(std::string filename, std::string outputFilename, int nsamples, int nvars,
		int c = 0);

#endif /* FULLSEARCH_H_ */
