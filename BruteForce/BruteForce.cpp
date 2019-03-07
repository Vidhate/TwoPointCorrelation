#include <iostream>
#include <string>
#include <cstring>
#include <vector>
#include <algorithm>
#include <fstream>
#include <cmath>
#include <random>
#include <chrono>
using namespace std;

//string path="./../../../DLA Mock Catalogue/";
string path="./../Data/";
vector<string> fileNames{"DLAhost_snap49_r1_b1 (copy).txt","DLAhost_snap49_r1_b1alpha (copy).txt","DLAhost_snap49_r1_b1T10 (copy).txt"};
double boxsize=150.0;	// Size of the simulation box in MPc/h
double numPart=1024.0;	// The number of particles used for the simulation will be cube of this number.
double forceRes=(boxsize/numPart)/30;	// Force resolution used in the simulation is 1/30th of mean inter-particle spacing.
double minLength=10.0*forceRes;	// The smallest bin should start from this length. minLength~ 5 to 10 times forceRes.
double maxLength=boxsize/4.0; /*The maximum length upto which bins should be split is 1/5th to 1/10th the size of the box for the following reasons:
- The boxes are periodic and any separation more than L/2 will require accessing other periodic boxes beyond the current one.
- Of the available perturbation modes in the box, the K_min(biggest wavelength)= 2*pi/boxlength. There are modes below this but they are constrained by the size of the simulation box. Hence, another maxlength cutoff on length greater than K_min ~ L/(5 to 10)
Being very flexible with this, I take the factor 4 for this run.*/

/*
Correlation structure holds the bincenter values and the correlation value at those bins. These can be plotted directly. These will be output in a file  for wasier plotting using python tools.
*/
struct Correlation{
	vector<double> centeredBinValues;
	vector<double> correlationValue_BF;		// Found using Brute Force method
};

/*
Splits a string "s" on the given delimiter "ch".
Returns an array of all individual words after splitting.
Reference: https://www.oreilly.com/library/view/c-cookbook/0596007612/ch04s07.html
*/
vector<string> split(string s, char ch){
	string::size_type i = 0;
	string::size_type j = s.find(ch);
	vector<string> afterSplit;

	while (j != string::npos) {
		afterSplit.push_back(s.substr(i, j-i));
		i = ++j;
		j = s.find(ch, j);

		if (j == string::npos)
			afterSplit.push_back(s.substr(i, s.length()));
	}

	return afterSplit;
}

/*
readFile takes in an index number as parameter. This index number references one of the file names from "fileNames".
The data from the given input file is then read into a 2D vector of doubles, sorted and returned.

References:
https://stackoverflow.com/questions/13035674/how-to-read-line-by-line-or-a-whole-text-file-at-once
*/
vector<vector<double> > readFile(int index=-1,string inpFileName=""){
	string dataFile;
	if(index==-1 && inpFileName.length()!=0){
		dataFile=inpFileName;
	}

	if(index!=-1 && inpFileName.length()==0){
		if(index>=fileNames.size()){
			cout<<"File names dictionary doesn't have this file\n";
			exit(0);
		}
		dataFile=path+fileNames.at(index);
	}

	if(index==-1 && inpFileName.length()==0){
		cout<<"Please give an input file\n";
		exit(0);
	}

	if(index!=-1 && inpFileName.length()!=0){
		cout<<"Ambiguity on input file. Please give either index from \"fileNames\" or directly the name of the file. Not both. \n";
		exit(0);
	}

	vector<vector<double> > pos;
	vector<double> buffer;
	vector<string> tokens;	// To receive input of split word from thesplitting function

	cout<<"Reading file : "<<dataFile<<endl;
	ifstream file(dataFile);
	string line;
	while(getline(file,line)){
		if(line[0]!='#'){
			buffer.clear();
			tokens.clear();
			tokens=split(line,' ');
			// Positions are in the index 1,2,3 of tokens. Known from the data.
			buffer.push_back(stod(tokens.at(1)));
			buffer.push_back(stod(tokens.at(2)));
			buffer.push_back(stod(tokens.at(3)));
			pos.push_back(buffer);
		}
	}
	file.close();
//	sort(pos.begin(),pos.end());
	cout<<"File read successfully"<<endl;
	return pos;
}

// Returns the smaller of 2 doubles
double min(double a, double b){
	if(a<b)return a;
	else return b;
}

/*
This function finds the absolute euclidean distance or separation between two points p1 and p2 in N-Dimensional space with considerations of a periodic box.
N will be found from the length of points p1 and p2.
Returns the separation between points.
Formula for 1D: dl=min(abs(x2-x1),L-abs(x2-x1)) 	| where L is the length of the periodic box.
Available in 2 modes - periodic | normal
*/
double separation(vector<double> p1, vector<double> p2){
	if(p1.size()!=p2.size() || p1.size()==0 || p2.size()==0){
		cout<<"Inconsistent dimensions of 2 points while calculating separation\n";
		exit(0);
	}

	double p_distance=0.0;
	for(int i=0;i<p1.size();i++){
		double dl=abs(p2.at(i)-p1.at(i));
		p_distance+=pow(min(dl,boxsize-dl),2);
	}

/*
	for(int i=0;i<p1.size();i++){
		p_distance+=pow(p2.at(i)-p1.at(i),2);
	}	
*/
	return sqrt(p_distance);
}

/*
If no input in given on bin width then a default number of bins are considered: nBins=maxLength/2 | dR=2.0
Bins can be split linearly or logarithmically according to the parameter passed: "lin"=linear | "log"=logarithmic. Default is logarithmic.
*/
vector<double> createBins(double dR=2.0, string mode="log", int numOfBins=30){
	if(mode!="log" && mode!="lin"){
		cout<<"The mode selected for splitting the bins in unavailable or wrong. Please recheck the spelling/availability\n";
		exit(0);
	}

	if(dR>maxLength){
		cout<<"Bin width larger than available box. Please rectify and run again\n";
		exit(0);
	}

	vector<double> bins;
	
	if(mode=="lin"){
		int nBins=(maxLength-minLength)/dR;
		double mark=minLength;
		for(int i=0;i<nBins;i++){
			bins.push_back(mark);
			mark+=dR;
		}
	}

	if(mode=="log"){
		double limit=log10(maxLength/minLength);
		double mark=0.0;
		double step=limit/double(numOfBins-1);		// To create n bins between [0,x] the step should be x/(n-1)
		for(int i=0;i<numOfBins;i++){
			bins.push_back(minLength*pow(10,mark));
			mark+=step;
		}
	}

	return bins;
}

/*
Takes as parameter the separation between 2 points and the available bins.
Returns the bin index in which this particular pair will be binned using binary search.
*/
int findBin(double sep, vector<double> bins){
	int lo=0;
	int hi=bins.size()-1;
	int mid;
	if(sep>bins.at(hi)){return -1;}
	while(lo<=hi){
		mid=(lo+hi)/2;
		if(sep>=bins[mid] && sep<bins[mid+1]){
			return mid;
		}else{
			if(sep<bins[mid])
				hi=mid-1;
			else lo=mid+1;
		}
	}
	return -1;
}

/*
Function to populate the given positions vector with random positions (3D). Used to find 2-point correlations.
Reference:
https://en.cppreference.com/w/cpp/numeric/random/uniform_real_distribution
*/
void randomDistribution(vector<vector<double> > &posRR, long long n){
	random_device rd;
	mt19937 gen(rd());
	vector<double> bufferRR;
	uniform_real_distribution<> dis(0.0,boxsize);
	for(long long i=0;i<n;i++){
		bufferRR.clear();
		for(int j=0;j<3;j++){
			bufferRR.push_back(dis(gen));
		}
		posRR.push_back(bufferRR);
	}
//	sort(posRR.begin(),posRR.end());
}


/* This function adds image of the dataset around the original dataset to simulate periodicity at a certain scale
*/
/*
vector<vector<double> > periodize(vector<vector<double> > pos, double scale=(boxsize/2.0)){
	
	vector<vector<double> >  new_pos;
	for(long i=0;i<pos.size();i++){
		new_pos.push_back(pos.at(i));
	}
	// Expanding in x
	long len=new_pos.size();
	vector<double> x;
	for(long i=0;i<len;i++){
		x=new_pos.at(i);
		if(x.at(0)>=scale){
			vector<double> dummy{x.at(0)-boxsize,x.at(1),x.at(2)};
			new_pos.push_back(dummy);
		}
		if(x.at(0)<=scale){
			vector<double> dummy{x.at(0)+boxsize,x.at(1),x.at(2)};
			new_pos.push_back(dummy);
		}
	}
	
	// Expanding in y
	len=new_pos.size();
	for(long i=0;i<len;i++){
		x=new_pos.at(i);
		if(x.at(1)>=scale){
			vector<double> dummy{x.at(0),x.at(1)-boxsize,x.at(2)};
			new_pos.push_back(dummy);
		}
		if(x.at(1)<=scale){
			vector<double> dummy{x.at(0),x.at(1)+boxsize,x.at(2)};
			new_pos.push_back(dummy);
		}
	}

	// Expanding in z
	len=new_pos.size();
	for(long i=0;i<len;i++){
		x=new_pos.at(i);
		if(x.at(2)>=scale){
			vector<double> dummy{x.at(0),x.at(1),x.at(2)-boxsize};
			new_pos.push_back(dummy);
		}
		if(x.at(2)<=scale){
			vector<double> dummy{x.at(0),x.at(1),x.at(2)+boxsize};
			new_pos.push_back(dummy);
		}
	}

	return new_pos;
}

*/

/*
Doing counts using brute force pair counting.
input required: DDpositions, RRpositions, bins, (empty vector for) DD Counts, (empty vector for) RR Counts
output give: Populates DD counts and RR counts and returns DD/RR - 1 (the correlation function Xi(r))
Reference for timing the code:
https://www.geeksforgeeks.org/measure-execution-time-function-cpp/
*/
vector<double> countsBruteForce(vector<vector<double> > positions, vector<vector<double> > positionsRR, 
	vector<double> bins,vector<long long> &countsDD, vector<long long> &countsRR){

	cout<<"Starting pair counts"<<endl;
	auto startBF=chrono::high_resolution_clock::now();

	long long n=positions.size();
	double sepDD,sepRR;
	for(long long i=0;i<n;i++){
		for(long long j=0;j<n;j++){
			sepDD=separation(positions.at(i),positions.at(j));
			sepRR=separation(positionsRR.at(i),positionsRR.at(j));
			int bDD=findBin(sepDD,bins);
			if(bDD!=-1)
				countsDD.at(bDD)++;
			int bRR=findBin(sepRR,bins);
			if(bRR!=-1)
				countsRR.at(bRR)++;
		}
	}

	vector<double> DDbyRR;
	double ratio;
	for(int i=0;i<countsDD.size();i++){
		ratio=double(countsDD.at(i)+1)/double(countsRR.at(i)+1);
		DDbyRR.push_back(ratio-1.0);
	}

	auto stopBF=chrono::high_resolution_clock::now();
	auto duration=chrono::duration_cast<chrono::seconds>(stopBF - startBF);

	cout<<"Brute Force pair counts done successfully\n";
	cout<<"Time taken to accomplish : "<<duration.count()<<" seconds\n";

	return DDbyRR;
}

/*
This function calculates 2 point correlation for the given distribution of positions of N points in a cubic box(presumption).
*/
//Correlation twoPtCorrelation(vector<vector<double> > positions, double dR=2.0){		// Original
Correlation twoPtCorrelation(vector<vector<double> > positions, vector<vector<double> > positionsRR, double dR=2.0){			// Inserted
	cout<<"Starting 2-pt\n";
	Correlation result;
	long long n=positions.size();
	vector<double> bins;
	bins=createBins(dR);
	int m=bins.size();
	cout<<"Bins created successfully\n";

	vector<long long> countsDD(m-1,0);
	vector<long long> countsRR(m-1,0);

	// Genarting random positions
/*	vector<vector<double> > positionsRR;			// Original
	positionsRR.clear();
	randomDistribution(positionsRR,n);
	cout<<"Random distribution generated successfully\n";*/

	// Pair counting Brute Force and Computing 2 point correlation as ordered ratio of (countsDD/countsRR)-1.0
	result.correlationValue_BF=countsBruteForce(positions,positionsRR,bins,countsDD,countsRR);

	// Centering bins
	vector<double> centeredBins;
	for(int i=0;i<m-1;i++)
		centeredBins.push_back((bins.at(i)+bins.at(i+1))/2.0);
	result.centeredBinValues=centeredBins;

	return result;
}

/*
To write the Correlation result into a file called "OUT "+input_file_name.
References:
https://stackoverflow.com/questions/478075/creating-files-in-c
Note: Always put \n for starting newline. "endl" will cause unkown characters while reading output from (let's say) Python
*/
void writeResult(Correlation result, int index){
	string outFileName="BF "+fileNames[index];
	ofstream outfile(outFileName);
	outfile<<"# BF 2 point correlation for input file : "<<path+fileNames[index]<<"\n";
	outfile<<"# Separation[R]	Correlation Value[Xi(R)]\n";
	for(int i=0;i<result.centeredBinValues.size();i++){
		outfile<<result.centeredBinValues.at(i)<<" "<<result.correlationValue_BF.at(i)<<"\n";
	}
	outfile.close();
	cout<<"Successfully written to Output file\n\n";
}


int main(){
	Correlation result[fileNames.size()];
	for(int i=0;i<fileNames.size();i++){
		vector<vector<double> > positions;
		vector<vector<double> > positionsRR;
		positions=readFile(i);
		positionsRR=readFile(-1,(path+"RANDOM "+fileNames[i]));	// inserted
		result[i]=twoPtCorrelation(positions,positionsRR);	// inserted
//		result[i]=twoPtCorrelation(positions);				// original
		writeResult(result[i],i);
	}
	return 0;
}
