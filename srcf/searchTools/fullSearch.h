/*
 * fullSearch.h
 *
 *  Created on: 3 May 2018
 *      Author: Omar
 */
#include <string>
#include <iostream>
#include <vector>
#ifndef FULLSEARCH_H_
#define FULLSEARCH_H_



double * calculateMeasures(int p1, double Hp1, int p2, double Hp2, int p3,
		int cl, const int *d,int nsamples, int nvars, int c, double Hcl, double Hp1cl);
std::vector<int> readData(std::string filename, int nrows, int nvars);
void runFullSearch(std::string filename, std::string outputFilename, int nsamples, int nvars,
		int c, bool printall, double assocLevel);
void runFullSearchIndexes(std::string filename, std::string outputFilename, int nsamples, int nvars,
		int c, int istart, int iend, int jstart, int jend, int kstart, int kend,
		bool printall, double assocLevel);
std::vector<int> readIndices(std::string filename, int n);
void runSearchVIFilter(std::string filename, std::string outputFilename, int nsamples, int nvars,
		int c, int istart, int iend, int jstart, int jend, int kstart,int kend, double epsilon,
		bool printall, double assocLevel, std::vector<double> viDists);
void runSearchClusterFilter(std::string filename, std::string outputFilename, int nsamples, int nvars,
		int c, int istart, int iend, int jstart, int jend, int kstart,int kend,
		bool printall, double assocLevel, std::vector<int> clusterIDs, std::vector<int> clusterSizes);

void runPAM(std::string datafile, std::string outputClusters, int nsamples, int nvars, size_t k);
#endif /* FULLSEARCH_H_ */
