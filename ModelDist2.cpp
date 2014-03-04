#include <unistd.h>
#include <ios>
#include <omp.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <bitset>
#include <sys/resource.h>

using namespace std;
 double pi = 3.14159L;
////Call is BitHashCompare Parent Mutant firstpassfile hashsize

/////////////////////////

const vector<string> Split(const string& line, const char delim) {
    vector<string> tokens;
    stringstream lineStream(line);
    string token;
    while ( getline(lineStream, token, delim) )
        tokens.push_back(token);
    return tokens;
}


void process_mem_usage(double& vm_usage, double& resident_set, double& MAXvm, double& MAXrss)
{
   using std::ios_base;
   using std::ifstream;
   using std::string;

   vm_usage     = 0.0;
   resident_set = 0.0;

   // 'file' stat seems to give the most reliable results
   //
   ifstream stat_stream("/proc/self/stat",ios_base::in);

   // dummy vars for leading entries in stat that we don't care about
   //
   string pid, comm, state, ppid, pgrp, session, tty_nr;
   string tpgid, flags, minflt, cminflt, majflt, cmajflt;
   string utime, stime, cutime, cstime, priority, nice;
   string O, itrealvalue, starttime;

   // the two fields we want
   //
   unsigned long vsize;
   long rss;

   stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
               >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
               >> utime >> stime >> cutime >> cstime >> priority >> nice
               >> O >> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest

   stat_stream.close();

   long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
   vm_usage     = vsize / 1024.0;
   resident_set = rss * page_size_kb;
        if (vm_usage > MAXvm){MAXvm = vm_usage;}
        if (resident_set > MAXrss){MAXrss = resident_set;}
}

double norm(double x, double mu, double sigma)
{
	return ( 1 / (sqrt(2* pi * pow(sigma, 2)))  ) * exp( -( pow(x-mu, 2) / (2*pow(sigma, 2))));
}
double calcError(vector<long>& histo, vector< double>& ModelSums, long Reads, int ReadLength, int HS)
{
	int eH = 10;
	
	double Pe = histo[1]/((double) Reads * (double) ReadLength * (double) eH);
	cout << "Pe = " << histo[1] << "/((double)" <<  Reads << " * (double) " << ReadLength <<" * (double) " << eH << ");" << " = " << Pe <<  endl;
	//cout << "HS = " << HS <<endl;
	
//	cout << "Error Rate = " << Pe << endl;

	long TotalMode = 0;
	for(int i = 1; i < ModelSums.size(); i++)
	{
		TotalMode+=ModelSums[i];
	}
	for(int i = 1; i < 10; i++)
	{
		double SC = 99;
		double sum = 0;
		sum +=  (pow(Pe, i) * ( ((SC-i)*(1-Pe))/(SC-(i-1)))) * Reads * eH * ReadLength;
//		cout << "old method = " << sum <<endl;
		for (SC = 1; SC < ModelSums.size(); SC++)
		{
//			cout << "     SC = " << SC << " = " << pow(Pe, i) * ( ((SC-i)*(1-Pe))/(SC-(i-1))) << endl << " * " << ModelSums[SC] << endl;	
			sum +=  (pow(Pe, i) * ( ((SC-i)*(1-Pe))/(SC-(i-1)))) * ModelSums[SC] ;/// eH * ReadLength;
		}
//		cout <<"    " << i << "\t" << sum << endl;
	} 
}	

double testModel(long SC, double stdev, double factor, vector<long>& histo)
{
	vector< vector<double> > dist;
        vector< vector<double> > prob;
        vector<double> placeHolder;


        dist.push_back(placeHolder);
        prob.push_back(placeHolder);
        for (long i = 1; i < histo.size(); i++)
        {
                vector<double> d;
                vector<double> p;

                double total = 0;
                d.push_back(0);
                p.push_back(0);
                for (long j = 1; j < histo.size()/SC; j++)
                {
                        d.push_back(norm(i, SC*j, stdev + ((j-1) * factor)));
                        total+=d[j];
                }
                dist.push_back(d);
                for (long j = 1; j<d.size(); j++)
                {
                        p.push_back(d[j]/total);
                }
                prob.push_back(p);
        }

	vector<double> RC;
        RC.push_back(0);
        RC.push_back( histo[SC] / norm(SC, SC, stdev));
        for(long a = 2; a < histo.size()/SC; a++)
        {
                double Reads = 0;
		Reads = (histo[SC*a] / norm(SC*a, SC*a, stdev+((a-1)*factor))) * prob[SC*a][a];
                //for (long b = 1; b < prob.size(); b++)
                //{
                //        Reads += (histo[b] * prob[b][a]);
                //}
                RC.push_back(Reads);
        }

	vector< vector<double> > model;	
	vector<double> ModelSums;
	ModelSums.push_back(0);
 	for (long i = 1; i < histo.size(); i++)
        {
                vector<double> m;
                double sum = 0;
                m.push_back(0);
                for (long j = 1; j < histo.size()/SC; j++)
                {
                        m.push_back( norm(i, SC*j, stdev + ((j-1) * factor)) * RC[j]);
                        sum+=m[j];
                }
                model.push_back(m);
                ModelSums.push_back(sum);
        }

        double SumSquares = 0;
        for (long i = SC-stdev; i < SC*5; i++)
        {
                SumSquares += pow(histo[i]-ModelSums[i], 2);
        }
        //cout << "Sum of squares is " << SumSquares << " S = " << stdev << " f = " << factor <<endl;
        return SumSquares;

}	
int main (int argc, char *argv[])
{
	 cout << "Call is histoFile HS Threads" << endl;
	//for(long i = 1; i < 20; i++)
	//cout << i << " = " << norm(i, 10 , 5) << endl;	
	
	ifstream HistoFile;		
	HistoFile.open (argv[1]);
	if ( HistoFile.is_open())
	{ cout << "Parent File open - " << argv[1] << endl;}	//cout << "##File Opend\n";
	else
	{
		cout << "Error, HistoFile could not be opened";
		return 0;
	}
	

	ofstream ProbFile;
	string path = argv[1];
	path += ".prob";
	ProbFile.open(path.c_str());
	if (ProbFile.is_open())
	{}else{cout << "Error, Prob file could not be opened"; return 0;}	
	
	ofstream ModelFile;
        path = argv[1];
        path += ".model";
        ModelFile.open(path.c_str());
        if (ModelFile.is_open())
        {}else{cout << "Error, Model file could not be opened"; return 0;}

	ofstream DistFile;
        path = argv[1];
        path += ".dist";
        DistFile.open(path.c_str());
        if (DistFile.is_open())
        {}else{cout << "Error, Model file could not be opened"; return 0;}

	double SC = -1;
	double SCvalue = -1;
	double stdi = -1;
	double stdvalue = -1;
	long value = -1;
	long HistoSum = 0;
	int HS = -1;
	string par = argv[2];
	HS=atoi(par.c_str());

	int Inflection = -1;
	string line;
	line = argv[3];
	int Threads = atoi(line.c_str());
	long factor = 1;
	vector<string> temp;
	vector<long> histo;
	histo.push_back(0);
	getline(HistoFile, line);
	getline(HistoFile, line);
	temp = Split(line, '\t');
	value = atol(temp[1].c_str());
	histo.push_back(value);
	long last = value;
	
	bool PastInflection = false;
	long i = 1;
	while (getline(HistoFile, line) )
	{
	//	cout << line << endl;
		i++;
		temp = Split(line, '\t');
		value = atol(temp[1].c_str());
		histo.push_back(value);
		HistoSum+=value;
		if (value-last>0 && PastInflection == false)
		{
			Inflection = i-1;
			PastInflection = true;
		}
		if (PastInflection)
		{
			if (SCvalue < value)
			{
				SCvalue = value;
				SC = i;
			}
		}
	//	cout << histo[i] << endl;
		last = histo[i]; 
	}
	for( i = 0; i < 10 ;i++)
	{
		cout << "I = " << i << " \t " << histo[i] << endl;
	}
	cout << "SC = " << SC << " vlaue = " << SCvalue << endl;
	stdvalue = SCvalue * exp(-.5);
	for (i = SC; i<histo.size(); i++)
	{
		if (histo[i]-stdvalue<0)
		{stdi = i; break;}
	}
	double stdev = i-SC;
	cout << "stdi = " << i << " stdev = " << stdev << endl;	
	

	double bestSQ = testModel( SC,  stdev, factor,  histo);
	double bestS = stdev;
	double bestF = factor;
	double bestSC = SC;
	#pragma omp parallel for num_threads(Threads)
	for (int f = 1; f <=20; f+=1)
	{
		for (double s = stdev-2; s<=stdev+2; s+=.25)	
		{
			for( double m = SC-2; m<=SC+2; m+=.5) 
			{
				double SQ =  testModel( m,  s, f,  histo);
				#pragma omp critical
				{
					if (SQ <= bestSQ)	
					{
						bestSQ = SQ;
						bestS = s;
						bestF = f;
							bestSC = SC;
					}
				}
			}
		}
	}

	cout << "Best Model is SC = " << bestSC << " StdDev = " << bestS << " F = " << bestF << " with SSQ = " << bestSQ << endl;
	vector< vector<double> > model;
	vector< vector<double> > dist;
        vector< vector<double> > prob;
	vector<double> placeHolder;
	model.push_back(placeHolder);
	dist.push_back(placeHolder);
        prob.push_back(placeHolder);

	stdev = bestS;
	factor = bestF;
	SC = bestSC;

	for (i = 1; i < histo.size(); i++)
	{
		vector<double> d;
		vector<double> p;

		double total = 0;
		d.push_back(0);
		p.push_back(0);
		for (long j = 1; j < histo.size()/SC; j++)
		{
			d.push_back(norm(i, SC*j, stdev + ((j-1) * factor)));
			total+=d[j];
		}
		dist.push_back(d);
		
		for (long j = 1; j<d.size(); j++)
		{
			p.push_back(d[j]/total);
		}
		
		prob.push_back(p);
	}
	
	vector<double> RC;
        RC.push_back(0);
 	RC.push_back( histo[SC] / norm(SC, SC, stdev));
        for(long a = 2; a < histo.size()/SC; a++)
        {
                double Reads = 0;
                 Reads = (histo[SC*a] / norm(SC*a, SC*a, stdev+((a-1)*factor))) * prob[SC*a][a];
		//for (long b = 1; b < prob.size(); b++)
                //{
                //        Reads += (histo[b] * prob[b][a]);
                //}
                RC.push_back(Reads);
        }
	cout << RC.size() << endl;
	//for (i = 1; i<50; i++)
	//{cout << "i = " << i << ", RC = " <<  RC[i] << "\n" << endl;}

        vector<double> ModelSums;
        ModelSums.push_back(0);
        for (long i = 1; i < histo.size(); i++)
        {
                vector<double> m;
                double sum = 0;
                m.push_back(0);
                for (long j = 1; j < histo.size()/SC; j++)
                {
                        m.push_back( norm(i, SC*j, stdev + ((j-1) * factor)) * RC[j]);
                        sum+=m[j];
                }
                model.push_back(m);
                ModelSums.push_back(sum);
        }
	for (long i = 1; i < SC; i++)	
	{
		prob[i][1] = model[i][1]/histo[i];	
	}
	
    	if (SC - (5*stdev) > Inflection)
	{
		ModelFile <<  3 << endl << SC - (5*stdev) << endl;
	}
	else
	{
		ModelFile <<  3 << endl << Inflection << endl;
	}
	ModelFile << HistoSum << endl;
	ProbFile << HistoSum << endl;
	DistFile << HistoSum << endl;
	for (long CopyNumber = 1; CopyNumber<model[1].size(); CopyNumber++)
	{
		long LocalSum = 0;
		for (long KmerCount = 1; KmerCount < histo.size(); KmerCount++)
		{
			 LocalSum += model[KmerCount][CopyNumber];	
		}
		ModelFile << ((double)LocalSum) / ((double)HistoSum) <<"\t";
		ProbFile <<  ((float)LocalSum) / ((float)HistoSum) <<"\t";
	        DistFile <<  ((float)LocalSum) / ((float)HistoSum) <<"\t";
	}
	ModelFile << endl;
	ProbFile << endl;
	DistFile << endl;
	ModelFile << "K count\tRawCount\tModelSum\t1x\t2x\t3x" << endl;
	for (long KmerCount = 1; KmerCount < SC*5; KmerCount++)
	{
		ModelFile << KmerCount << "\t" << histo[KmerCount] ;
		 ModelFile << '\t' << ModelSums[KmerCount];
		for (long CopyNumber = 1; CopyNumber< 10; CopyNumber++)
		{
			ModelFile << '\t' << model[KmerCount][CopyNumber];
		}
		ModelFile<< endl;
	}	
	ModelFile.close();
	ProbFile << SC << endl;	
	for (long c = 1; c < prob.size(); c++)
        {
                ProbFile << c ;
                for (long k = 1; k<prob[1].size(); k++)
                {
                        ProbFile << '\t' << prob[c][k];
                }
                ProbFile<< endl;
        }
	 DistFile << SC << endl;
        for (long c = 1; c < dist.size(); c++)
        {
                DistFile << c ;
                for (long k = 1; k<dist[1].size(); k++)
                {
                        DistFile << '\t' << dist[c][k];
                }
                DistFile<< endl;
        }

	long GenomeSize = 0;
	//for(int i =1; i<RC.size(); i++)
	//{
	//	cout << i << "x = " << RC[i] << " hashes  = " << (RC[i] * 18) / (i) << endl;	
	//	GenomeSize += ((RC[i]*18)/ (i));
	//}	
	long TotalGoodKmers;	
	for (int i = 1; i < ModelSums.size(); i++)
	{
		if (ModelSums[i] > 0 )
		{TotalGoodKmers += ModelSums[i]*i;}
		//cout << "i = " << i << " Model = " <<  ModelSums[i] << " TK = " << TotalGoodKmers << endl;
	}
	
	GenomeSize =  (TotalGoodKmers/(100-HS+1) * 100)/SC;
	
	cout << "GenomeSize = " << GenomeSize << endl;
	cout << "Inflection point = " << Inflection << endl;
	cout << "Recomended RUFUS cutoff = " << SC - (5*stdev) << endl;
	cout << "-1std = " << SC - (1*stdev) << "-2std = " << SC - (2*stdev) <<"-3std = " << SC - (3*stdev) <<"-4std = " << SC - (4*stdev) << endl;
	double bam = calcError(histo, ModelSums, GenomeSize * SC / 100 , 100 , HS);
	ProbFile.close();
        DistFile.close();
        ModelFile.close();
}







































