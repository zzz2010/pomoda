#pragma once
#include "flmatrix.h"


class BGModel :
	public FlMatrix
{
public:
	BGModel(void);
	BGModel(string bgfile, int markovOrder=3);
	BGModel(map<string, double>& instSet);
	~BGModel(void);
public:
	int Order;
	FlMatrix Markov0;
	FlMatrix Markov1;
	FlMatrix Markov2;	
	double getStringProb(string pattern);
	void MixMotifModel(void* m,double factor=0.5);
	void MixBGModel(BGModel& bg,double factor=0.5);
	double getPWMThreshold(FlMatrix* pwm,double FDR=0.01);
	void GenBGData(int taglength);
};
