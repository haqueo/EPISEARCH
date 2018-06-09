/*
 * utilities.h
 *
 *  Created on: 6 Jun 2018
 *      Author: Omar
 */
#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#ifndef UTILITIES_UTILITIES_H_
#define UTILITIES_UTILITIES_H_


std::vector<int> readData(std::string filename, int nrows, int nvars);
std::vector<int> readIndices(std::string filename, int n);
std::vector<double> readviDistances(std::string filename,int nvars);
std::vector<int> readClusterIDs(std::string filename, int nvars);
std::vector<int> readClusterSizes(std::string filename, int numClusters);
int getIndexNextCluster(int thisClusterID, const std::vector<int> &clusterIDs, const std::vector<int> &clusterSizes,
		int numClusters, int nvars);
void writeStatistics(double measure0sum, double measure0squaredsum,
		double measure1sum, double measure1squaredsum, double measure2sum,
		double measure2squaredsum, double measure3sum,
		double measure3squaredsum, double measure4sum,
		double measure4squaredsum, int iterations);
void printMeasures(bool printall, const std::ofstream& myfile, int i, int j, int k,
		double assocLevel, double Hcl, double* measureArray);

#endif /* UTILITIES_UTILITIES_H_ */
