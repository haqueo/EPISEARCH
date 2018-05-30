//============================================================================
// Name        : fullSearch.cpp
// Author      : Omar Haque, but almost all of the entropy code is by Patrick E.
// 				 Meyer found in the infotheo package.
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
//using namespace std;
#include <algorithm>
#include <ctime>
#include <cmath>
#include <vector>
#include <map>
#include <fstream>
using namespace std;
#include <string>
#include <sstream>
#include <iostream>
#include <math.h>




double digamma(double z) {
      if(z<=0) return 0;
      double zp, zpr, zprs, digam = 0;
         zp = z;
      while(zp < 30) {
             zpr = 1/zp;
             digam -= zpr;
             zp++;
      }
         zpr = 1/zp;
         zprs = zpr * zpr;
         digam += log(zp)+zpr*(-0.5+zpr*(-1.0/12.0+zprs*(1.0/120.0-zprs/252.0)));
      return digam;
}

double entropy_empirical(std::map< std::vector<int> ,int > frequencies, int nb_samples) {
      double e = 0;
      for (std::map< std::vector<int> ,int>::const_iterator iter = frequencies.begin(); iter != frequencies.end(); ++iter)
            e -= iter->second * log((double)iter->second);
      return log((double)nb_samples) + e/nb_samples;
}

double entropy_miller_madow(std::map< std::vector<int> ,int > frequencies, int nb_samples) {
      return entropy_empirical(frequencies,nb_samples) + (int(frequencies.size())-1)/(2.0*nb_samples);
}

double entropy_dirichlet(std::map< std::vector<int> ,int > frequencies, int nb_samples, double beta) {
      double e = 0;
      for (std::map< std::vector<int> ,int>::const_iterator iter = frequencies.begin(); iter != frequencies.end(); ++iter)
            e+=(iter->second+beta)*(digamma(nb_samples+(frequencies.size()*beta)+1)-digamma(iter->second+beta+1));
      return e/(nb_samples+(frequencies.size()*beta));
}

double entropy_shrink(std::map< std::vector<int> ,int > frequencies, int nb_samples)
{
      double w = 0;
      int p = frequencies.size(), n2 = nb_samples*nb_samples;
      double lambda, beta;
      for (std::map< std::vector<int> ,int>::const_iterator iter = frequencies.begin(); iter != frequencies.end(); ++iter)
            w += iter->second*iter->second;
         lambda = p*(n2 - w)/((nb_samples-1)*(w*p - n2));
      if(lambda >= 1)
        return -log(1.0/p);
      else {
            beta = (lambda/(1-lambda))*nb_samples/frequencies.size();
        return entropy_dirichlet(frequencies, nb_samples, beta);
      }
}

double entropyFast(const int *d, int nsamples, int nvars, int c,int *vnew){

// H(d) using estimator c
	std::map< std::vector<int> ,int > freq;
	std::vector<int> sel;
	int j;
	int nsamples_ok = 0;
	double H = 0;
	for(int s = 0; s < nsamples; ++s) {

		sel.clear();

		for(int i = 0; i < 4; i++){
			j = vnew[i];
			if (j != -1){
				sel.push_back(d[s+j*nsamples]);
			} else{
				break;
			}

		}

		freq[sel]++;
		nsamples_ok++;

	}
	if( c == 0 ) //empirical
		H = entropy_empirical(freq,nsamples_ok);
	else if( c == 1 ) //miller-madow
		H = entropy_miller_madow(freq,nsamples_ok);
	else if( c == 2 ) //dirichlet Schurmann-Grassberger
		H = entropy_dirichlet(freq,nsamples_ok, 1/freq.size());
	else if( c == 3 ) // shrink
		H = entropy_shrink(freq,nsamples_ok);
	return H;
}

double entropy(const int *d, int nsamples, int nvars, int c, bool *v) {
// H(d) using estimator c
	std::map< std::vector<int> ,int > freq;
	std::vector<int> sel;
	bool ok = true;
	int nsamples_ok = 0;
	double H = 0;
	for(int s = 0; s < nsamples; ++s) {
		ok = true;
		sel.clear();
		for(int i = 0; i < nvars; ++i) {
			if(v[i]){
				if(d[s+i*nsamples]!=500) // change this later TODO
					sel.push_back(d[s+i*nsamples]);
				else
					ok = false;
			}
		}
		if(ok) {
			freq[sel]++;
			nsamples_ok++;
		}
	}
	if( c == 0 ) //empirical
		H = entropy_empirical(freq,nsamples_ok);
	else if( c == 1 ) //miller-madow
		H = entropy_miller_madow(freq,nsamples_ok);
	else if( c == 2 ) //dirichlet Schurmann-Grassberger
		H = entropy_dirichlet(freq,nsamples_ok, 1/freq.size());
	else if( c == 3 ) // shrink
		H = entropy_shrink(freq,nsamples_ok);
	return H;
}

double multiinformation(const int *d, int nsamples, int nvars, int c) {
		 bool *sel = new bool[nvars];
		 double sum = 0;
		 for( int i=0; i<nvars; ++i )
			sel[i] = false;
		 for(int i=0;i<nvars; ++i) {
			sel[i] = true;
			sum += entropy(d, nsamples, nvars, c, sel);
			sel[i] = false;
		 }
		 for( int i=0; i<nvars; ++i )
			sel[i] = true;
		 sum -= entropy(d, nsamples, nvars, c, sel);
		 return sum;
}



double interaction(const int *d, int nsamples, int nvars, int c) {
	double sum = 0; //
	int sign = 1;//
	bool *sel= new bool[nvars];//
	for(int i = 0; i < nvars; ++i)//
		sel[i] = false; //

	std::vector<int> indx;
	int n = nvars;
	int j=1;
	int k=n;
	for(int twk=j;twk<=k;twk++){
		int r=twk;
		bool done=true;
		for(int iwk=0;iwk<r;iwk++)indx.push_back(iwk);
		while(done){
			done=false;
			for(int owk=0;owk<r;owk++) //
				sel[indx[owk]]=true;	//
			sum += sign*entropy(d, nsamples, nvars, c, sel); //
			for(int owk=0;owk<r;owk++) //
				sel[indx[owk]]=false;	//
			for(int iwk=r-1;iwk>=0;iwk--){
				if(indx[iwk]<=(n-1)-(r-iwk)){
					indx[iwk]++;
					for(int swk=iwk+1;swk<r;swk++)
						indx[swk]=indx[swk-1]+1;
					iwk=-1;
					done=true;
				}
			}
		}
		sign = -sign; //
		indx.clear();
	}
	return sum;
}

double *calculateMeasures(int p1, double Hp1, int p2, double Hp2, int p3, int cl, const int *d,
		int nsamples, int nvars, int c,double Hcl, double Hp1cl){
	static double final[5];
	int selTempNew[4] = {-1,-1,-1,-1};

	////// calculation of base entropies and joint entropies

	// H(P1)
	// this is entropyp1

	// H(P1,C)
	// this is entropyp1cl

	// H(P2)
	// this is entropyp2

	// H(P2,C)
	selTempNew[0] = p2;
	selTempNew[1] = cl;
	double Hp2cl = entropyFast(d,nsamples,nvars,0,selTempNew);

	// H(P3,C)

	selTempNew[0] = p3;
	double Hp3cl = entropyFast(d,nsamples,nvars,0,selTempNew);

	// H(P3)

	selTempNew[1] = -1;
	double Hp3 = entropyFast(d,nsamples, nvars,0,selTempNew);

	// H(P1,P2)

	selTempNew[0] = p1;
	selTempNew[1] = p2;
	double Hp1p2 = entropyFast(d,nsamples, nvars,0,selTempNew);

	//H(P1,P2,C)

	selTempNew[2] = cl;
	double Hp1p2cl = entropyFast(d,nsamples,nvars,0,selTempNew);

	// H(P1,P3)

	selTempNew[2] = -1;
	selTempNew[1] = p3;

	double Hp1p3 = entropyFast(d,nsamples,nvars,0,selTempNew);

	// H(P1,P3,C)

	selTempNew[2] = cl;
	double Hp1p3cl = entropyFast(d,nsamples,nvars,0,selTempNew);

	// H(P2,P3)
	selTempNew[0] = p2;
	selTempNew[2] = -1;

	double Hp2p3 = entropyFast(d,nsamples,nvars,0,selTempNew);

	// H(P2,P3,cl)
	selTempNew[2] = cl;
	double Hp2p3cl = entropyFast(d,nsamples,nvars,0,selTempNew);

	// H(P1,P2,P3)
	selTempNew[2] = p1;

	double Hp1p2p3 = entropyFast(d,nsamples,nvars,0,selTempNew);

	// H(P1,P2,P3,C)
	selTempNew[3] = cl;
	double Hp1p2p3cl = entropyFast(d,nsamples,nvars,0,selTempNew);

	///////////// calculation of compound measures
	double Ip1p2p3 = Hcl - Hp1p2p3cl + Hp1p2p3;
	double Ip1 = Hcl - Hp1cl + Hp1;
	double Ip2 = Hcl - Hp2cl + Hp2;
	double Ip3 = Hcl - Hp3cl + Hp3;
	double Ip1p2 = Hcl - Hp1p2cl + Hp1p2;
	double Ip1p3 = Hcl - Hp1p3cl + Hp1p3;
	double Ip2p3 = Hcl - Hp2p3cl + Hp2p3;
	double IGp1p2 = Ip1p2 - Ip1 - Ip2;
	double IGp1p3 = Ip1p3 - Ip1 - Ip3;
	double IGp2p3 = Ip2p3 - Ip2 - Ip3;

	///////////// calculation of VI
	double VIp1p2 = 2*Hp1p2 - Hp1 - Hp2;
	double VIp1p3 = 2*Hp1p3 - Hp1 - Hp3;
	double VIp2p3 = 2*Hp2p3 - Hp2 - Hp3;
	double VIp1cl = 2*Hp1cl - Hp1 - Hcl;
	double VIp2cl = 2*Hp2cl - Hp2 - Hcl;
	double VIp3cl = 2*Hp3cl - Hp3 - Hcl;


	double IG_strict = Ip1p2p3 - max(IGp1p2,0.0) - max(IGp1p3,0.0) - max(IGp2p3,0.0) - Ip1 - Ip2 - Ip3;
	double IG_alt = Ip1p2p3 - IGp1p2 - IGp1p3 - IGp2p3 - Ip1 - Ip2 - Ip3;
	double IG_new = Ip1p2p3;
	double VIp = VIp1p2 + VIp1p3 + VIp2p3;
	double VIc = VIp1cl + VIp2cl + VIp3cl;

	final[0] = IG_strict;
	final[1] = IG_alt;
	final[2] = IG_new;
	final[3] = VIp;
	final[4] = VIc;
	return final;
}



std::vector<int> readData(std::string filename, int nrows, int nvars){
	/*
	 * Format data should be: nrows x nvars csv of integers 0,1 (minor or major alleles)
	 * the Class type is the last of these column included in nvars. So its index is
	 * nvars-1, i.e. the last column.
	 *
	 *
	 * Output: the data in a 1d std::vector, arranged by columns first. So c1^T, c2^T,...,cN^T
	 */

	const int datasize = nvars*nrows;
	static std::vector<int> data(datasize);
	static std::vector<int> dataReal(datasize);

    std::fstream myfile(filename.c_str(), std::ios_base::in);

    int a;

    for(int i =0; i< datasize;i++){
    	myfile >> a;
    	data[i] = a;
    }

    for(int s = 0; s < nrows; s++){
    	for (int i = 0; i<nvars;i++){
    		dataReal[s+i*nrows] = data[i+s*nvars];
    	}

    }

    // this reads in the data as is. I actually need to transpose it.
	return dataReal;

}

std::vector<int> readIndices(std::string filename, int n){

	std::ifstream infile(filename.c_str());
	std::string s;
	std::vector<int> indexes;

   //skip N lines
   for(int i = 0; i < n; ++i){
	   std::getline(infile, s);
   }
   // s is now the n'th line

   std::replace( s.begin(), s.end(), ',', ' ' );

   stringstream ss(s);

   copy(istream_iterator<int>(ss), istream_iterator<int>(), back_inserter(indexes));

   return indexes;
}

void runFullSearch(std::string filename, std::string outputFilename, int nsamples, int nvars,
		int c, bool printall, double assocLevel){
	std::vector<int> d = readData(filename, nsamples, nvars);

	const int* p = d.data();

	int v[4] = {-1,-1,-1,-1};
	v[0] = nvars-1;

	double Hcl = entropyFast(p,nsamples,nvars,c,v);
	v[0] = -1;

	double measure0sum = 0.0;
	double measure0squaredsum = 0.0;
	double measure1sum = 0.0;
	double measure1squaredsum = 0.0;
	double measure2sum = 0.0;
	double measure2squaredsum = 0.0;
	double measure3sum = 0.0;
	double measure3squaredsum = 0.0;
	double measure4sum = 0.0;
	double measure4squaredsum = 0.0;
	int iterations = 0;



	ofstream myfile (outputFilename.c_str());

	if (myfile.is_open()){

	for (int i = 0; i < nvars-3;i++){
		v[0] = i;
		double Hp1 = entropyFast(p,nsamples,nvars,c,v);
		v[1] = nvars-1;
		double Hp1cl = entropyFast(p,nsamples,nvars,c,v);

		v[0] = -1;
		v[1] = -1;

		for (int j = i+1;j < nvars-2;j++){

			v[0] = j;
			double Hp2 = entropyFast(p,nsamples,nvars,c,v);
			v[0] = -1;

			for(int k = j+1; k < nvars-1; k++){
				double *measureArray = calculateMeasures(i, Hp1, j, Hp2, k, nvars-1, p,
						nsamples, nvars, c, Hcl, Hp1cl);

				if (printall){
				myfile << "(" << i << "," << j << ","<< k << "): " << *(measureArray+0) <<
						"	"<< *(measureArray+1) <<"	"<< *(measureArray+3) <<
						"	"<< *(measureArray+4) << "\n";
				} else if (*(measureArray+0) > assocLevel*Hcl){

				myfile << "(" << i << "," << j << ","<< k << "): " << *(measureArray+0) <<
										"	"<< *(measureArray+1) <<"	"<< *(measureArray+3) <<
										"	"<< *(measureArray+4) << "\n";
				}

				measure0sum+=*(measureArray+0);
				measure0squaredsum+=pow(*(measureArray+0),2);
				measure1sum+=*(measureArray+1);
				measure1squaredsum+=pow(*(measureArray+1),2);
				measure2sum+=*(measureArray+2);
				measure2squaredsum+=pow(*(measureArray+2),2);
				measure3sum+=*(measureArray+3);
				measure3squaredsum+=pow(*(measureArray+3),2);
				measure4sum+=*(measureArray+4);
				measure4squaredsum+=pow(*(measureArray+4),2);
				iterations++;
			}
		}
	}

	myfile.close();
	} else {
		cout << "Unable to open file";
	}



	ofstream myfile2 ("means.txt");

	if (myfile2.is_open()){
		myfile2 << "averageIGstrict: " << measure0sum/iterations << std::endl;
		myfile2 << "stddevIGstrict: " << sqrt(measure0squaredsum/iterations - (pow(measure0sum/iterations,2))) << std::endl;

		myfile2 << "averageIGalt: " << measure1sum/iterations << std::endl;
		myfile2 << "stddevIGalt: " << sqrt(measure1squaredsum/iterations - (pow(measure1sum/iterations,2))) << std::endl;

		myfile2 << "averageVIp: " << measure3sum/iterations << std::endl;
		myfile2 << "stddevVIp: " << sqrt(measure3squaredsum/iterations - (pow(measure3sum/iterations,2))) << std::endl;


		myfile2 << "averageVIc: " << measure4sum/iterations << std::endl;
		myfile2 << "stddevVIc: " << sqrt(measure4squaredsum/iterations - (pow(measure4sum/iterations,2))) << std::endl;

		myfile2 << "iterations:" << iterations << std::endl;

		myfile2.close();
	}



}

void runFullSearchIndexes(std::string filename, std::string outputFilename, int nsamples, int nvars,
		int c, int istart, int iend, int jstart, int jend, int kstart, int kend,
		bool printall, double assocLevel){

	/* (istart,jstart,jend) to (iend,jend,kend) INCLUSIVE BOTH SIDES */


	double measure0sum = 0.0;
	double measure0squaredsum = 0.0;
	double measure1sum = 0.0;
	double measure1squaredsum = 0.0;
	double measure2sum = 0.0;
	double measure2squaredsum = 0.0;
	double measure3sum = 0.0;
	double measure3squaredsum = 0.0;
	double measure4sum = 0.0;
	double measure4squaredsum = 0.0;
	int iterations = 0;


	std::vector<int> d = readData(filename, nsamples, nvars);
	// read in the array of indexes through filenameIndexes

	const int* p = d.data();

	int v[4] = {-1,-1,-1,-1};
	v[0] = nvars-1;

	double Hcl = entropyFast(p,nsamples,nvars,c,v);
	v[0] = -1;

	ofstream myfile (outputFilename.c_str());

	bool Done = false;

	if (myfile.is_open()){

	int jstartreal;
	int kstartreal;

	for (int i = istart; !Done; i++){
		v[0] = i;
		double Hp1 = entropyFast(p,nsamples,nvars,c,v);
		v[1] = nvars-1;
		double Hp1cl = entropyFast(p,nsamples,nvars,c,v);

		v[0] = -1;
		v[1] = -1;

		if (i == istart){
			jstartreal = jstart;
		} else {
			jstartreal = i+1;
		}


		for (int j = jstartreal; (j < nvars-2) && !Done; j++){

			v[0] = j;
			double Hp2 = entropyFast(p,nsamples,nvars,c,v);
			v[0] = -1;

			if (j == jstart){
				kstartreal = kstart;
			} else {
				kstartreal = j+1;
			}

			for(int k = kstartreal; (k < nvars-1) && !Done; k++){
				double *measureArray = calculateMeasures(i, Hp1, j, Hp2, k, nvars-1, p,
						nsamples, nvars, c, Hcl, Hp1cl);

				if (printall){
				myfile << "(" << i << "," << j << ","<< k << "): " << *(measureArray+0) <<
						"	"<< *(measureArray+1) <<"	"<< *(measureArray+3) <<
						"	"<< *(measureArray+4) << "\n";
				} else if (*(measureArray+0) > assocLevel*Hcl){
				myfile << "(" << i << "," << j << ","<< k << "): " << *(measureArray+0) <<
										"	"<< *(measureArray+1) <<"	"<< *(measureArray+3) <<
										"	"<< *(measureArray+4) << "\n";
				}

				measure0sum+=*(measureArray+0);
				measure0squaredsum+=pow(*(measureArray+0),2);
				measure1sum+=*(measureArray+1);
				measure1squaredsum+=pow(*(measureArray+1),2);
				measure2sum+=*(measureArray+2);
				measure2squaredsum+=pow(*(measureArray+2),2);
				measure3sum+=*(measureArray+3);
				measure3squaredsum+=pow(*(measureArray+3),2);
				measure4sum+=*(measureArray+4);
				measure4squaredsum+=pow(*(measureArray+4),2);
				iterations++;

				Done = (i >= iend) && (j >= jend) && (k >= kend);
			}
		}
	}

	myfile.close();
	} else {
		cout << "Unable to open file";
	}


	ofstream myfile2 ("means.txt");

	if (myfile2.is_open()){
		myfile2 << "averageIGstrict: " << measure0sum/iterations << std::endl;
				myfile2 << "stddevIGstrict: " << sqrt(measure0squaredsum/iterations - (pow(measure0sum/iterations,2))) << std::endl;

				myfile2 << "averageIGalt: " << measure1sum/iterations << std::endl;
				myfile2 << "stddevIGalt: " << sqrt(measure1squaredsum/iterations - (pow(measure1sum/iterations,2))) << std::endl;

				myfile2 << "averageVIp: " << measure3sum/iterations << std::endl;
				myfile2 << "stddevVIp: " << sqrt(measure3squaredsum/iterations - (pow(measure3sum/iterations,2))) << std::endl;


				myfile2 << "averageVIc: " << measure4sum/iterations << std::endl;
				myfile2 << "stddevVIc: " << sqrt(measure4squaredsum/iterations - (pow(measure4sum/iterations,2))) << std::endl;



			myfile2.close();
		}





}

