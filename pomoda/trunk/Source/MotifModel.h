#pragma once
#include "HashEngine.h"
#include "flmatrix.h"

#include "intmatrix.h"
#include "stdafx.h"

#define FOR(i,n) for(i=0; i<n; i++)

class MotifModel :
	public FlMatrix
{
public:
	MotifModel(void);
	MotifModel(HashEngine* , int maxlen,PARAM* setting);
	~MotifModel(void);
	double BinPvalue;
	int Efflen;
	PARAM* Setting;
	// a batch of instances are considered belonged to this motif
	vector<string> InstanceSet;
	// Add Instance to  instanceSet and update the PWM
	void AddInstance(string instanceString);
	void Merge(MotifModel* P);
	void MergeList(vector<MotifModel*> simList);
	// remove the lowest likelood score instance
	void RemoveWeakInstance(void);
	// compute serveral scores of this motif model
	void ComputeScore(double sampleratio, double & CDScore, double & ORScore, double & BindingRegion, double & CNSVScore, double & DiffScore);
	// MCMC method generate instance from this motif model
	map<string,double> GenerateInstanceFromPWM(double sampleratio,bool noThreshold=false);
	map<string,double> GenerateInstanceFromPWMPQ(double sampleratio,bool noThreshold=true);
	vector<int> columnOrder(int col);
	void resetPWM();
	HashEngine* SearchEngine;
	// Compute center distributioin score
static	double ComputeCDScore(double* histgram);
	// Compute Intensity Score
static	double ComputeINTScore(double* Histogram);
	// Compute OverRepresented Score
	double ComputeORScore(map<string,double>& InstSet, vector<int> count,bool mflag=false);
	double getWindowSize(map<string,double>& InstSet,double threshold);
	// Compute Conservation Score
	double ComputeCNSVScore(void);
	//gen random instance
	pair<string,double> randomInstance(void);



	// get summary score of this motif model
	double GetMixedScore(void);
	double CDScore;
	double ORScore;
	double Threshold;
	double BindingRegion;
	double CNSVScore;
	double DiffScore;
	//double Ratio1;
	//double Ratio2;
	double SeqPvalue;




	
		//statistics var
	double cdavg,cdvar,intavg,intvar,oravg,orvar;
	bool switchFlag;
	int head,tail;
	double deltaScore;

	int scoreIndex;
	map<int,double> badmove;
	vector<int> instSeq;
	vector<string> BGInstanceSet;
	void PWMRefinement();
	string Consensus;
	string seed;
	// update the PWM by random search , and runid specify the column will be updated
	void divergeSeedPart(vector<int>* hotsites_p=NULL);
	double updateModel(int runid);
	double GetMixedScore(double CDScore, double ORScore,double INTScore, double CNSVScore, double DiffScore);
	double getLogProb(char* instance);
	void InitializePWMofInstanceSet();
	std::string get_consensus(double pseudo_weight = 0, double genome_fra = 0.5, double genome_frc = 0.5 );
	static vector<string> NCloseSet(string);
	static double SimilarityScore(vector<VAL>& plist1,vector<VAL>& plist2,int len1,int len2,int& commonCount,double alterP=0);
	static double SimilarityScore2(int len,int NumOfTry,int countCommon,int aln, int alnScore,int mismatch=1);
	static int Alignment(string& refmotif, string& newmotif,int& bestScore,int mismatch=1000);
	static int AlignmentRC(string& refmotif, string& newmotif, int& bestScore,int mismatch=1000);
	vector<MotifModel*> getSeedMotifs( int SeedNum,int len,int mismatch=1);
	map<string, double> MisMatchClosedSet(string pattern, int mis);
	void printPWM(string name);
	 double  AlignmentScore(string& refmotif, string& newmotif);
	  double AlignmentScoreRC(string& refmotif, string& newmotif);
	vector<VAL> getMatchPos();

	static int AlignmentPWM(MotifModel* refmotif,MotifModel* newmotif,double& bestScore,int & bestoverlap);
	static int AlignmentPWMRC(MotifModel* refmotif,MotifModel* newmotif,double& bestScore);
	MotifModel RC();
	double FittingCurve(double &snr,double& prbE,double& prbN);
	void printMatchPos(string name,vector<VAL>& list);
	void printRealPos(string name);
	static string Hash2ACGT(int hash,int len)
	{
		string ACGT="ACGT";
		int k;
		string pa="";
			FOR(k,len)
			{
				pa.push_back( ACGT[hash%4]);
				hash>>=2;
			}	
			return pa;
	};

	int Length()
	{
		return X()-head-tail;
	}

	vector<long int> POSLIST;
	int lastupdateCol;
	void MarkPos();
	static string printAlignment(string pa, int aln, int reflen)
	{
		string temp="";
		int i;
		FOR(i,reflen-1+aln)
			temp.push_back('N');
		FOR(i,pa.size())
			temp.push_back(pa[i]);
		int left=2*reflen-1-aln-pa.size();
		FOR(i,left)
			temp.push_back('N');
		return temp;
			
	};

	double entropy2(int col)
	{
		int k;
		double entropyT=0;
		FOR(k,4)
			entropyT-=g(col,k)*log(g(col,k))/log(2.0);
		return entropyT;
	}
		static void convert_int_to_fl_matrix (IntMatrix* i1, FlMatrix* f1) 
		{
			double sum;
			int i,j;
			FOR(i,i1->X()){
				sum = 0;
				FOR(j,4) {
					f1->s(i1->g(i,j),i,j);
					sum+=i1->g(i,j);
				}
				FOR(j,4) {
					f1->d(sum,i,j);
				}
			}
		};
private:
	double entropy(int col);
};
