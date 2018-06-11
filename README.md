# EPISEARCH

EPISEARCH is a fast, high-performance C++ command line application to search for epistatic triples of sites in genetic datasets. EPISEARCH runs 98% faster than the equivalent R code even without parallelisation.

## Binaries / Installation

EPISEARCH can only be used as a command line tool. Binaries for Mac OS and Linux are available from: https://sourceforge.net/projects/episearch/files/

In order to build the program locally just run: `$make all`

Requirements are C++11 and OpenMP

## Usage

To see the format for the required input files, see the examples folder.

Basic usage is: `$ ./EPISEARCH [OPTIONS]`

The options are given here:

## Main options: 

* -inputfile TEXT input file for analysis [required]
* -outputfile TEXT location and name for output file [required]
* -nsamples INT number of rows to the dataset [required]
* -nvars INT number of columns including the phenotypic class as the last [required]
* -numThread INT number of threads to paralellise over using OpenMP [default: 1]
* -indexfile TEXT index file for naive parallelisation (see examples)
* -nindex INT which line of the index file to use
* -printall y,n whether to print data for all triples iterated over 
* -assoclevel FLOAT if -printall is false, only print for sites with an association greater than -assoclevel (e.g. 0.05 for 5%)

## Extra options for filtering:

* -VIfilter y,n whether to use the Variation of Information filter
* -VIfilterfile TEXT text file for VI distances (see examples)
* -epsilon FLOAT Only search triples for VI distance greater than epsilon
* -clusterfilter y,n whether to use a clustering as a filter
* -clusterfilterIDsfile TEXT the file of cluster IDs for sites (see example)
* -numClusters INT number of clusters used

## Basic example:
The following codes searches through all the triples in genetic.txt and returns all the triples with scores such that the strict information gain is greater than 5%.
```
./EPISEARCH -inputfile ./genetic.txt -outputfile ./output.txt -nsamples 400 -nvars 4000 -assoclevel 0.05 -printall n
```
See the examples folder for more in-depth examples using the filters to reduce the search space.
