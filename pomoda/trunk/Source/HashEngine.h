#pragma once
#include <iostream>
#include <stdlib.h>
#include <map>
#include "MiscUtilities.h"
#include "MemManager.h"
#include "TextConverter.h"
#include <math.h>
#include <time.h>
#include "ISearchEngine.h"

using namespace std;

class HashEngine :
	public ISearchEngine
{
public:
	HashEngine(void);
	HashEngine(const char* filename, int hashlen=8,string wfile="");
	~HashEngine(void);
	virtual int searchPattern( char* pattern, int mismatches, vector<VAL>& position);
virtual	 VAL getTotalLength();
	 // The size of the text
	 VAL TotalLength;
	 char * gPattern;
	 string weightfile;;
	 int Hashlen;
	 vector<VAL>** HashIndex;
	 double BGProb[4];
	 double CDProb[4];
	 vector<double> SeqWeighting;
	 VAL buildIndex(const char* textfile);
	 // use hash1 hash2 hash3 form the base hash address and then use the range to specify the hash adress range
	 void JoinMerge(vector<VAL>& position,int Range, int hash1, int hash2=-1, int hash3=-1);
	 // use the most conservative part as a anchor, go back the hashtext to check the satisfication
	 void SmartShift(vector<VAL>& position, char* pattern);
	 // Coding Text into a Char array
	 char* CharText;
	 int SeqLen;
	 int searchPatternNRC( char* pattern,int pos, int mismatches, vector<VAL>& position);
	 string getSite(VAL pos, int len);
	 HashEngine* Clone();
	 void LoadWeightFile(string filename);
	 double getSeqWeight(int seqid);
	 bool EnableWeight;

};
