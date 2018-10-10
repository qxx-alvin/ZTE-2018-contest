#include <iostream>
#include <vector>
#include <stdio.h>
#include <algorithm>
#include <numeric>
#include <time.h>
#include <fstream>
using namespace std;


extern int nLink;
extern int nReq;
extern int nPath;

extern vector<int> linkCapa;			// capacity of all links, including reverse links
extern vector<int> reqBw;				// bandwidth demand of all requests
extern vector<vector<vector<int>>> reqPaths;			// alternative paths for all requests
			

int nChromos = 10;					// size of population
double CopyPercentageLow = 0.2;			// direct copy to next generation
double CopyPercentageHigh = 0.2;
int nCross;// = nChromos * (1 - CopyPercentage);
int nCopy;// = nChromos - nCross;
//double MutationProb = 0.1;
//int nMutation = nChromos * MutationProb;
int migIt = 50000;			// migration every 50 iterations
double migRateLow = 0.2;
double migRateHigh = 0.2;
//int nMig = nChromos * migRate;


vector<vector<int>> Chromos;               // population
vector<vector<int>> NewChromos;

int ItNum = 40000;

vector<double> Fitness;					// fitness for each chromo
vector<double> NewFitness;				// fitness for each new chromo
vector<double> MaxFit;					// max fitness of each iteration
vector<double> AvFit;					// average fitness of each iteration
vector<double> ChooseProb;				// cumulative choose probability of each chromo

void gaIt(int itNum);
void initPopulation();
double calcOneFitness(vector<int> &chromo);
void calcFitness();
void calcChooseProb();
void cross(int num); 
void mutation(double av);
void copy();
inline int randInt(int a, int b);
int rws();
void migration(int it);
void adaCross(double avFit, double maxFit);
void adaMutation(double avFit, double maxFit);
void adaCopy();
void calcNewFitness();

void ga()
{
//	srand(time(NULL));
	srand(6);
	// GA iteration

	gaIt(ItNum);

	ofstream out("AvFit.txt");
	for (int i = 0; i < ItNum; ++i)
		out << AvFit[i] << " " << MaxFit[i] << endl;
	out.close();
}

void gaIt(int itNum)
{
	initPopulation();
	Fitness.reserve(nChromos);
	NewFitness.reserve(nChromos);
	MaxFit.reserve(itNum);
	ChooseProb.reserve(nChromos);

	for (int i = 0; i < nChromos; ++i)
	{
		NewFitness.push_back(0);
		Fitness.push_back(0);
		ChooseProb.push_back(0);
	}
	for (int i = 0; i < itNum; ++i)
	{
		MaxFit.push_back(0);
		AvFit.push_back(0);
	}

	for (int i = 0; i < itNum; ++i)
	{
		

//		int nCro = nCross;
		calcFitness();
		AvFit[i] = accumulate(Fitness.begin(), Fitness.end(), 0.0) / nChromos;
		auto maxFitPosition = max_element(Fitness.begin(), Fitness.end());
		MaxFit[i] = *maxFitPosition;
		if (i % 500 == 0)
			cout << i << ": " << 1 - MaxFit[i] << endl;
		if (i > 0 && MaxFit[i] < MaxFit[i - 1])
			cout << " bad" << endl;
//		cout << i << ": " << 1 - AvFit[i] << endl;

		double copyPer = CopyPercentageHigh - (CopyPercentageHigh - CopyPercentageLow) * (i * 1.0 / ItNum);//(AvFit[i] / MaxFit[i]);
		nCopy = nChromos * copyPer;
		nCross = nChromos - nCopy;

		calcChooseProb();
		NewChromos.clear();
		NewChromos.reserve(nChromos);
		if (i % migIt == 0)
		{
			migration(i);
//			nCro -= nMig;
		}
		cross(nCross);
		copy();

		//calcNewFitness();
		//double newAv = accumulate(NewFitness.begin(), NewFitness.end(), 0.0) / nChromos;
		//double newMax = *max_element(NewFitness.begin(), NewFitness.end());
		mutation(AvFit[i]);
		//adaMutation();

//		mutation();
		//adaCross(AvFit[i], MaxFit[i]);
		//adaCopy();

		//calcNewFitness();
		//double newAv = accumulate(NewFitness.begin(), NewFitness.end(), 0.0) / nChromos;
		//double newMax = *max_element(NewFitness.begin(), NewFitness.end());
		//adaMutation(newAv, newMax);

		Chromos = NewChromos;
	}
}


void initPopulation()
{
	Chromos.reserve(nChromos);
	for (int i = 0; i < nChromos; ++i)
	{
		vector<int> chromo(nReq);
		for (int j = 0; j < nReq; ++j)
			chromo[j] = randInt(0, nPath);
		Chromos.push_back(chromo);
	}
}


void calcFitness()
{
	int nBadChromo = 0;
	for (int i = 0; i < nChromos; ++i)
		Fitness[i] = calcOneFitness(Chromos[i]);
//	cout << "bad: " << nBadChromo << endl;
}

void calcNewFitness()
{
	int nBadChromo = 0;
	for (int i = 0; i < nChromos; ++i)
		NewFitness[i] = calcOneFitness(NewChromos[i]);
}


double calcOneFitness(vector<int> &chromo)
{
	vector<double> linkUsed(nLink, 0);

	double fit = 0;
	bool badChromo = false;

	for (int j = 0; j < nReq; ++j)
	{
		for (int k = 0; k < reqPaths[j][chromo[j]].size(); ++k)
		{
			int idx = reqPaths[j][chromo[j]][k];
			linkUsed[idx] += reqBw[j];
			if (linkUsed[idx] > linkCapa[idx])
			{
				badChromo = true;
				break;
			}
		}
		if (badChromo)
			break;
	}

	if (!badChromo)
	{
		for (int j = 0; j < nLink; ++j)
			linkUsed[j] /= linkCapa[j];
		fit = 1 - *max_element(linkUsed.begin(), linkUsed.end());
	}

	return fit;
}

void calcChooseProb()
{
	double fitSum = accumulate(Fitness.begin(), Fitness.end(), 0.0);
	ChooseProb[0] = Fitness[0] / fitSum;
	for (int i = 1; i < nChromos; ++i)
		ChooseProb[i] = ChooseProb[i - 1] + Fitness[i] / fitSum;
}


void adaCross(double avFit, double maxFit)
{
	int nNew = 0;
	double K1 = 1, K2 = 0.5, K3 = 1, K4 = 0.5;
	while (nNew < nCross)
	{
		// choose a pair of parents
		int faInd = randInt(0, nChromos);
		int moInd = randInt(0, nChromos);

		// determine cross probability
		double pc;
		int maxInd = Fitness[faInd] > Fitness[moInd] ? faInd : moInd;
		if (Fitness[maxInd] < avFit)
			pc = K3;
		else
			pc = K1 * (maxFit - Fitness[maxInd]) / (maxFit - avFit);

		// copy
		//if (pc == 0)
		//{
		//	NewChromos.push_back(Chromos[maxInd]);
		//	nNew++;
		//}
//		else

		if (rand() / (RAND_MAX + 1.0) >= pc)
			continue;

		vector<int> father = Chromos[faInd];
		vector<int> mother = Chromos[moInd];

		int crossIdx = randInt(1, nReq);

		vector<int> child;
		child.reserve(nReq);
		child.insert(child.end(), father.begin(), father.begin() + crossIdx);
		child.insert(child.end(), mother.begin() + crossIdx, mother.end());

		NewChromos.push_back(child);
		nNew++;

	}
}

void cross(int nCro)
{
	for (int i = 0; i < nCro; ++i)
	{
		// choose a pair of parents
		vector<int> father = Chromos[rws()];
		vector<int> mother = Chromos[rws()];

		int crossIdx = randInt(1, nReq);

		vector<int> child;
		child.reserve(nReq);
		child.insert(child.end(), father.begin(), father.begin() + crossIdx);
		child.insert(child.end(), mother.begin() + crossIdx, mother.end());

		//double fit = calcOneFitness(child);
		//if (fit < calcOneFitness(father) && fit < calcOneFitness(mother))		// not a good child
		//{
		//	i--;
		//	continue;
		//}
		NewChromos.push_back(child);
	}
}

void adaMutation(double avFit, double maxFit)
{
	double K1 = 1, K2 = 0.5, K3 = 1, K4 = 0.5;
	for (int i = 0; i < nChromos; ++i)
	{
		// determine mutation prob
		double pm;
		double fit = NewFitness[i];
		if (fit < avFit)
			pm = K4;
		else
			pm = K2 * (maxFit - fit) / (maxFit - avFit);

		if (pm == 0)
			pm = 0.01;
		if (rand() / (RAND_MAX + 1.0) >= pm)
			continue;
		
		int req = randInt(0, nReq);			// a random request
		int path = randInt(0, nPath);		// a random path
		NewChromos[i][req] = path;
	}
}

void mutation(double av)
{
	for (int i = 0; i < nCross; ++i)
	{
		//if (i == nCross)
		//	continue;
		vector<int> chromo;
//		int idx = randInt(0, nCross);		// a random chromosome
		int idx = i;
		//if (idx == nCross)
		//	continue;
		int req = randInt(0, nReq);			// a random request
		int path = randInt(0, nPath);		// a random path
		NewChromos[idx][req] = path;
		//chromo = NewChromos[idx];
		//chromo[req] = path;
		//if (calcOneFitness(chromo) > calcOneFitness(NewChromos[idx]))
		//	NewChromos[idx] = chromo;
	}
}


void copy()
{
	vector<int> visited(nChromos, 0);

	for (int i = 0; i < nCopy; ++i)		// choose the best nCopy chromosome
	{
		int maxInd = 0;
		while (visited[maxInd])
			maxInd++;

		for (int j = maxInd + 1; j < nChromos; ++j)
		{
			if (!visited[j] && Fitness[maxInd] < Fitness[j])
				maxInd = j;
		}
		NewChromos.push_back(Chromos[maxInd]);
	}
}


void adaCopy()
{
	int nCur = NewChromos.size();
//	cout << "nCur: " << nCur << endl;
	for (int i = 0; i < nChromos - nCur; ++i)		// choose the best nCopy chromosome
	{
		int ind = i;
		for (int j = i + 1; j < nChromos; ++j)
		{
			if (Fitness[i] < Fitness[j])
			{
				double tmp = Fitness[j];
				Fitness[j] = Fitness[i];
				Fitness[i] = tmp;
				ind = j;			// record the index
			}
		}
		NewChromos.push_back(Chromos[ind]);
	}
}

void migration(int it)
{
	int nMig = nChromos * (migRateLow + (migRateHigh - migRateLow) * (it * 1.0 / ItNum));// (AvFit[it] / MaxFit[it]));
	//vector<double> FitCopy;
	//FitCopy.assign(Fitness.begin(), Fitness.end());
	vector<int> visited(nChromos, 0);

	for (int i = 0; i < nMig; ++i)		// choose the worst nMig chromosome
	{
		int minInd = 0;
		while (visited[minInd])
			minInd++;

		for (int j = minInd + 1; j < nChromos; ++j)
		{
			if (!visited[j] && Fitness[minInd] > Fitness[j])
				minInd = j;
		}
		// new chromosome
		vector<int> chromo(nReq);
		for (int j = 0; j < nReq; ++j)
			chromo[j] = randInt(0, nPath);

		Chromos[minInd] = chromo;
		Fitness[minInd] = calcOneFitness(chromo);
		visited[minInd] = 1;
	}
}


int rws()
{
	double r = rand() / (RAND_MAX + 1.0);
	for (int i = 0; i < nChromos; ++i)
		if (ChooseProb[i] >= r && Fitness[i] != 0)
			return i;
}


inline int randInt(int a, int b)
{
	return rand() / (RAND_MAX + 1.0) * (b - a) + a;
}