//============================================================================
// Name        : fullSearch.cpp
// Author      : Omar Haque
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
using namespace std;




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

double entropyFast(const int *d, int nsamples, int nvars, int c, bool*v){

// H(d) using estimator c
	std::map< std::vector<int> ,int > freq;
	std::vector<int> sel;
	int nsamples_ok = 0;
	double H = 0;
	for(int s = 0; s < nsamples; ++s) {



		sel.clear();
		for(int i = 0; i < nvars; ++i) {
			if(v[i]){
					sel.push_back(d[s+i*nsamples]);
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
				if(d[s+i*nsamples]!=500) // change this later
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


