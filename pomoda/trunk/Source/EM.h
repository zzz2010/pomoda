#pragma once
#include "stdafx.h"
#include "MotifModel.h"

#define FOR(i,n) for(i=0; i<n; i++)
class EM
{
public:
	EM(PARAM* setting);
	MotifModel* BackModel;
	void LoadSeqFile(string filename);
	vector<MotifModel*> LoadSeedModels(vector<MotifModel*> candidates,int outMotifNum);
	vector<long> EMIterate(MotifModel* seed);
	vector<string> DATASET;
	vector<double> Prior;
	double** ErasingFactor;
	PARAM* Setting;

};

