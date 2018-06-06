/*
 * utilities.h
 *
 *  Created on: 6 Jun 2018
 *      Author: Omar
 */
#include <string>
#include <iostream>
#include <vector>
#ifndef UTILITIES_UTILITIES_H_
#define UTILITIES_UTILITIES_H_


std::vector<int> readData(std::string filename, int nrows, int nvars);
std::vector<int> readIndices(std::string filename, int n);
std::vector<double> readviDistances(std::string filename);
std::vector<int> readClusterIDs(std::string filename);
std::vector<int> readClusterSizes(std::string filename);
int getIndexNextCluster(int thisClusterID, std::vector<int> clusterIDs, std::vector<int> clusterSizes,
		int numClusters, int nvars);

#endif /* UTILITIES_UTILITIES_H_ */
