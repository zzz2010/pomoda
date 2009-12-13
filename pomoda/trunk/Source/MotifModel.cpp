#include "stdafx.h"
#include "MotifModel.h"
#include "MiscUtilities.h"
#include "HashEngine.h"



#define PWMThreshold 0.001


	int SEQLEN=4000;
	bool LargeDataFlag=true;

int BINNUMBER=100;
int MAXSEQNUM=17440;
int maxActRLen=5000;
double CDThreshold=0;
	double bgfold=5;
	double wininc=1.1;
	
VAL LIBSIZE;
double ENTROPY_Threshold=1.9;
double tolerance=0.1;


struct ScoreObjStruct
{
	double cdscore;
	double orscore;
	double BindingRegion;
}; 
typedef ScoreObjStruct ScoreObj;




MotifModel::MotifModel(void)
: CDScore(0)
, ORScore(0)
, BindingRegion(0)
, CNSVScore(0)
, DiffScore(0),FlMatrix(16,4)
{
	cdavg=cdvar=intavg=intvar=oravg=orvar=0;
	scoreIndex=0;
	LIBSIZE=MAXSEQNUM*SEQLEN;
}

MotifModel::MotifModel(HashEngine* engine, int motiflen,PARAM* setting)
: CDScore(0)
, ORScore(0)
, BindingRegion(setting->startwSize)
, CNSVScore(0)
, DiffScore(0),FlMatrix(motiflen,4)
{
	cdavg=cdvar=intavg=intvar=oravg=orvar=0;
	ORScore=1;
	scoreIndex=0;
	Setting=setting;
	Efflen=Setting->seedlength;
	SearchEngine=engine;
	head=tail=(X()-Setting->seedlength)/2;
	
	initialise(0.25);
	LIBSIZE=engine->getTotalLength();
	
	SEQLEN=engine->SeqLen;
	MAXSEQNUM=LIBSIZE/SEQLEN;
	BinPvalue=0;
	if(SEQLEN<2000)
	{
       LargeDataFlag=false;
	}

	
	switchFlag=false;
	deltaScore=-MINSCORE;
	lastupdateCol=-1;
	badmove.clear();
	CDScore=ORScore=CNSVScore=DiffScore=0;
	BindingRegion=Setting->startwSize;
	BINNUMBER= ceil((double)SEQLEN/setting->resolution);
	
	
}

MotifModel::~MotifModel(void)
{
	
}

void MotifModel::InitializePWMofInstanceSet()
{
	if(InstanceSet.size()==0)
		return;
	int width=InstanceSet[0].size();
	initialise(0.25);
	int i,j,k;
	double basecount=PWMThreshold*InstanceSet.size();
	double* sumAll=new double [width];
	FOR(i,width)
	{
		sumAll[i]=PWMThreshold*4;	
		double sum[4];
		FOR(j,4)
			sum[j]=PWMThreshold;
		FOR(j,InstanceSet.size())
		{
			int col=acgt(InstanceSet[j][i]);
			if(col>3||col<0)
			{
				FOR(col,4)
					sum[col]+=PWMThreshold;
				sumAll[i]+=PWMThreshold*4;
				continue;
			}
			else
			{
				double weight=1;
				if(j<instSeq.size())
					weight=SearchEngine->getSeqWeight(instSeq[j]);
			sum[col]+=weight;
			sumAll[i]+=weight;
			}
		}

		FOR(j,4)
		{
			double temp=sum[j]/sumAll[i];
			if(sumAll[i]==basecount*4)
				temp=0.25;
			s(temp,(X()-width)/2+i,j);
		}
	}
	//this->print();
	delete [] sumAll;
}


 double MotifModel::SimilarityScore(vector<VAL>& plist1,vector<VAL>& plist2,int len1,int len2,int& commonCount,double alterP)
	{	int i,j,k;
		i=j=0;
		double countIntersect=0;
		double countUnion=0;
		int overlap=min(len1,len2);
		while(i<plist1.size()&&j<plist2.size())
		{
				int seqnum1=plist1[i]/SEQLEN;
				int pos1=plist1[i]%SEQLEN;
				int seqnum2=plist2[j]/SEQLEN;
				int pos2=plist2[j]%SEQLEN;
				if(seqnum1==seqnum2&&(pos2-pos1<=(len1-overlap)&&pos1-pos2<=(len2-overlap)))//abs((int)(pos2-pos1))<len
				{
					countIntersect+=1;
					countUnion-=1;
			
				}
				if(plist1[i]<plist2[j])
					i++;
				else if(plist1[i]>plist2[j])
					j++;
				else
				{
					i++;
					j++;
				}
		}

		commonCount=countIntersect;
		countUnion+=plist1.size()+plist2.size();
		double P=plist2.size()*len1/(double)LIBSIZE;
		P=max(P,alterP);
		P=min(0.5,P);
		double lamda=plist1.size()*P;
		//int minsize=min(plist1.size(),plist2.size());
		double score2=0;
	
		if(commonCount>plist1.size())
			commonCount=plist1.size();
		{
			score2+=binominalCDF(P,commonCount,plist1.size());
		}
		score2=1-score2;

		return (double)score2; 
	}

 //this function only work in 1-mismatch, and equal length
 double MotifModel::SimilarityScore2(int len,int NumOfTry,int countCommon,int aln, int alnScore,int mismatch)
	{	
		int i,j,k;
		i=j=0;
		double P=0;//
		int overlap=len+aln;
		if(aln>0)
			overlap=len-aln;
		int mis=overlap-alnScore;
		int left=len-overlap;
		if(mis>2*mismatch)
			return -1;
		if(mismatch==0)
		{
			P=1/(pow(4.0,left));
		}
		else if(mis==2*mismatch)
		{
			P=1/(3*pow(4.0,left));
		}
		else if(mis==1*mismatch)
		{
			double po=pow(4.0,left);
			P=5/(6*po);
		}
		else 
		{ //0-mismatch
			double po=pow(4.0,left);
			P=(4*left+overlap)/(6*po);
		}
		double lamda=NumOfTry*P;

		double score2=0;
		//score2=(countCommon-lamda)/(sqrt((double)NumOfTry*P*(1-P)));
		FOR(i,countCommon)
		{
			double temp=binominal(P,i,NumOfTry);
			score2+=temp;
			
		}
		score2=1-score2;

		return (double)score2; 
	}

  int MotifModel::AlignmentRC(string& refmotif, string& newmotif,int& bestScore,int mismatch)
  {
	   bestScore=-1000;
	   int aln1=Alignment(refmotif,newmotif,bestScore,mismatch);
	  string RC=reverseString(newmotif);
	  int bestScore2=-1000;
	  int aln2=Alignment(refmotif,RC,bestScore2,mismatch);
	  if(bestScore2>bestScore)
	  {
		  bestScore=bestScore2;
		newmotif=RC;
		return aln2;
	  }
	  
	  return aln1;
  }

 int  MotifModel::Alignment(string& refmotif, string& newmotif,int &alnScore,int mismatch)
 {
		int bestscore=-1000;
		int bestaln=0;
		int size=newmotif.size();
		int size2=refmotif.size();
		int minsize=min(refmotif.size(),newmotif.size());
		for(int i=0-size+1;i<size2;i++)
		{
			int curScore=0;
			 int pp = 0;
            if (i > 0)
                pp = i;
            int j;
			for(j=pp;j-pp<minsize;j++)
			{
				 if(j>=size2||(j-i)>=size)
					      break;
				 if(j>=0)
				{
					if(newmotif[j-i]==refmotif[j])
						curScore+=1;
					else if (refmotif[j]=='N'||newmotif[j-i]=='N')
						curScore=curScore;
					else
						curScore-=1;
				}
			}
			int overlap=j;
                        ////i>0
            if (i > 0)
                  overlap = j - i;
			if(bestscore<curScore&&(overlap-curScore)<=2*mismatch)
			{
				bestscore=curScore;
				bestaln=i;
			}
		}
		alnScore=bestscore;
		return bestaln;
 }

  double MotifModel::AlignmentScoreRC(string& refmotif, string& newmotif)
  {
	  double bestScore=-1000;
	   bestScore=AlignmentScore(refmotif,newmotif);
	  string RC=reverseString(newmotif);
	  double bestScore2=-1000;
	  bestScore2=AlignmentScore(refmotif,RC);
	  if(bestScore2>bestScore)
	  {
		  bestScore=bestScore2;
		newmotif=RC;
	
	  }
	  
	  return bestScore;
  }

 double  MotifModel::AlignmentScore(string& refmotif, string& newmotif)
 {
		int bestscore=-1000;
		int bestaln=0;
		int size=newmotif.size();
		int size2=refmotif.size();
		int minsize=min(refmotif.size(),newmotif.size());
		for(int i=0-size+1;i<size2;i++)
		{
			int curScore=0;
			for(int j=i;j-i<minsize;j++)
			{
				if(j>=size2)
					break;
				if(j>=0)
				{
					if(newmotif[j-i]==refmotif[j])
						curScore+=1;
					else if(refmotif[j]=='N'||newmotif[j-i]=='N')
						curScore=curScore;

					else
						curScore-=1;
				}
			}
			int overlap=minsize-i;
			if(i<0)
				overlap=minsize+i;

			if(bestscore<curScore)
			{
				bestscore=curScore;
				bestaln=i;
			}
		}

		return (double)bestscore/minsize;
 }


void MotifModel::printPWM(string name)
{
	int i,j,k;
	name+=".pwm";
	name="./result/"+name;
	 FILE* outstream=fopen(name.c_str(), "w+");
	
	 if(outstream==NULL)
	 {
		 fclose(outstream);
		 return;
	 }


	 fprintf(outstream,"#%s\n",get_consensus().c_str());
		fprintf(outstream,"#%f\t%f\t%f\t%f\n",GetMixedScore(),CDScore,BindingRegion,ORScore);
	 fprintf(outstream, "%s\n", "PO\tA\tC\tG\tT");
	 FOR(i,X()) {
		 fprintf(outstream, "%d", i+1);
        FOR(j,4) 
		{
			 fprintf(outstream, "\t%f", g(i,j));
		}

		 fprintf(outstream, "\n");

	 }
	  fclose(outstream);
}

double MotifModel::getLogProb(char* instance)
{
	int i;
		double score=0;
		FOR(i,X())
		{
			int a=acgt(*(instance+i));
			if(a>3||a<0)
			{
				score=MINSCORE;
				break;
			}
			double sss=g(i,a);
			if(sss==0)
			{
				score=MINSCORE;
				break;
			}
			if(sss==0.25)
				continue;
			score+=log(sss);
		}
	return score;
}


void MotifModel::divergeSeedPart(vector<int>* hotsites_p)
{

	int i,j,k,q;
	switchFlag=true;
	POSLIST.clear();
	get_consensus(0);
	int seedpos=Consensus.find(seed);
	if(seedpos==-1)
		return;

	double* divCnt=new double[Setting->seedlength];
	FOR(i,Setting->seedlength)
		divCnt[i]=0;
	InstanceSet.clear();
	map<string,double> inst=GenerateInstanceFromPWMPQ(0.8);
			map<string,double>::iterator Iter1=inst.begin();
string ACGT="ACGT";
int lastsize=0;
int nonConPos=0;
	while(Iter1!=inst.end())
	{
		 int centCNT,bgCNT;
				 centCNT=bgCNT=0;
		//search the non-revised one first
		char instemp2[64];
				strcpy(instemp2,Iter1->first.c_str());
		int rcindex=SearchEngine->searchPattern(instemp2,0,POSLIST);	//search for each non-mutate instance
				for(q=lastsize;q<POSLIST.size();q++)
				 {
					int pos=POSLIST[q]%SEQLEN;
					double bias=abs(pos-SEQLEN/2);
					if(bias<BindingRegion/2)
					{
						
						centCNT++;
					}
					else
						bgCNT++;
				 }
			   	double cdd=(double)(centCNT/BindingRegion);
				 double bgg=(double)(bgCNT/(SEQLEN-BindingRegion));
				 double pvalue=binominalTail(BindingRegion/SEQLEN,centCNT,bgCNT+centCNT);
				 if(pvalue>0.05)// reject if not significant
				 {
					 // erase the position list
					 POSLIST.erase(POSLIST.begin()+lastsize,POSLIST.end());
					 Iter1++;
					 continue;
				 }
		//mark the reverse complement to be negative
		for(q=rcindex;q<POSLIST.size();q++)
			POSLIST[q]=0-POSLIST[q];
		lastsize=POSLIST.size();
		InstanceSet.push_back(instemp2);
		FOR(i,Setting->seedlength)
		{
			int col=(X()-Setting->seedlength)/2+i;
			bool ncFlag=false;
			col=seedpos+i;
			FOR(j,4)
			{
				char instemp[64];
				strcpy(instemp,Iter1->first.c_str());
				if(ACGT[j]!=instemp[col])
				instemp[col]=ACGT[j];
				else
					continue;//skip the non-revised one

				 rcindex=SearchEngine->searchPattern(instemp,0,POSLIST); //search one-mismatch instance
			
				 centCNT=bgCNT=0;
				 for(q=lastsize;q<POSLIST.size();q++)
				 {
					int pos=POSLIST[q]%SEQLEN;
					double bias=abs(pos-SEQLEN/2);
					if(bias<BindingRegion/2)
						centCNT++;
					else
						bgCNT++;
				 }
				 double cdd=(double)(centCNT/BindingRegion);
				 double bgg=(double)(bgCNT/(SEQLEN-BindingRegion));
				  double pvalue=binominalTail(BindingRegion/SEQLEN,centCNT,bgCNT+centCNT);
				  double st=cdd/bgg;
				  if(pvalue>0.05||(st<ORScore))// reject if not significant or have smaller uncorrected score
				 {
					 POSLIST.erase(POSLIST.begin()+lastsize,POSLIST.end());
					 continue;
				 }
				for(q=rcindex;q<POSLIST.size();q++)
					POSLIST[q]=0-POSLIST[q];
				lastsize=POSLIST.size();
				ncFlag=true;
				InstanceSet.push_back(instemp);
				// remove background bias
				divCnt[i]+=SearchEngine->BGProb[j];
			}
			if(ncFlag)
			{
			nonConPos++;
			
			}
			
		}
		Iter1++;

	}
	
	nonConPos/=inst.size()*2;
	if(nonConPos<1)
		nonConPos=1;
	//a heuristics to determine the degree of divergence
	FOR(i,Setting->seedlength)
		divCnt[i]/=inst.size();
	

		MotifModel tempmodel(SearchEngine,X(),Setting);

	int width=X();
	Consensus=get_consensus(0);



	double newScore=0;
	double oldScore=GetMixedScore();
	
	
int len=Length();

int ct=len/2;
int lastpos=-1;
          // collect the instances selected above in the center region to form PWM
			FOR(i,POSLIST.size())
			{
						long int wpos=POSLIST[i];

						bool rc=(wpos<0);
						VAL upos=wpos;
						if(rc)
							upos=0-wpos;

				if(SearchEngine->CharText[upos+ct]=='X'||abs(upos-lastpos)<ct*2)
					continue;
						int seqnum=upos/SEQLEN;
						int pos=upos%SEQLEN;
						lastpos=upos;
				int windowsize=BindingRegion;
				double bias=abs(pos-SEQLEN/2);

				if(bias<windowsize/2)
				{
					
						string pa;
						if(rc)
						{
							if(hotsites_p!=NULL)
								hotsites_p->push_back(upos-tail);
							pa=SearchEngine->getSite(upos,len);//(upos-tail,X());
							
							pa=reverseString(pa);
						}
						else
						{
							if(hotsites_p!=NULL)
								hotsites_p->push_back(upos-head);
							pa=SearchEngine->getSite(upos,len);//(upos-head,X());
					
						}
					
					tempmodel.AddInstance(pa);
					if(SearchEngine->EnableWeight)
					tempmodel.instSeq.push_back(seqnum);
				}

			}

				

				if(tempmodel.InstanceSet.size()==0)
				{
					return ;
				}
				
				tempmodel.InitializePWMofInstanceSet();
		
				tempmodel.get_consensus(0);

					FOR(i,Setting->seedlength)
					{
						int col=seedpos+i;
						int id=acgt(Consensus[col]);
						col=tempmodel.head+seedpos+i;
						if(id>3||id<0)
							continue;
						tempmodel.d(1+divCnt[i],col,id); // diverge the orignal seed element by the divergence degree
					
						double sum=0;
						double sum2=0;
						FOR(j,4)
						{
							tempmodel.d(SearchEngine->CDProb[j],col,j);
						}


						FOR(j,4)
						{
							sum+=tempmodel.g(col,j);
					
						}
						FOR(j,4)
						{
							tempmodel.d(sum,col,j);
												
						}


					}
			

				POSLIST.clear();
				int effLen=0;
				int nCnt=0;
				Consensus=tempmodel.get_consensus(0);

				// check the result PWM, if so few nongap position, then reject this PWM
				FOR(i,Setting->seedlength)
				{  
					if(Consensus[seedpos+i]=='A'||Consensus[seedpos+i]=='C'||Consensus[seedpos+i]=='G'||Consensus[seedpos+i]=='T')
						effLen++;
					if(Consensus[seedpos+i]=='N')
						nCnt++;
				}
				if(effLen-nCnt<Setting->seedlength-3||nCnt>1)//
				{
					ORScore=1;
					return;
				}

			//write back to this
		FOR(i,X())
		{
			FOR(j,4)
				s(tempmodel.g(i,j),i,j);
		}
			tempmodel.InstanceSet.clear();
				double ors;
			//compute the score again for the final PWM
			tempmodel.ComputeScore(0.8,CDScore,ORScore,BindingRegion,CNSVScore,DiffScore);
	
			delete []divCnt;
	
}

//mark the position by 'X'
void MotifModel::MarkPos()
{
		int i,j,k;
	vector<VAL> poslist;
	double threshold=log(PWMThreshold);
	 HashEngine* SearchEngine2=(HashEngine*)SearchEngine;
	 map<string,double> insts=GenerateInstanceFromPWMPQ(0.8);
	
	 map<string,double>::iterator Iter;
	double totalCount=0;
	for(Iter=insts.begin();Iter!=insts.end();Iter++)
	{	
		char instemp[64];
		string randomIns=Iter->first;
		strcpy(instemp,randomIns.c_str());
		int rcindex=SearchEngine->searchPattern(instemp,0,poslist);	

	}
	int len=insts.begin()->first.size();
	
	FOR(i,poslist.size())
	{
		k=poslist[i];
		k=abs(k);
		int pos=k%SEQLEN;
		if(abs(pos-SEQLEN/2)<BindingRegion)//(SEQLEN/bgfold)/2
		{
			FOR(j,len)
			 {
				*(SearchEngine2->CharText+k+j)='X';

			 }
		}
	}

	 POSLIST= poslist;

}

vector<VAL> MotifModel::getMatchPos()
{
		int i,j,k;

	vector<VAL> poslist;
	int maxRange=1<<(2*Length());
	if(maxRange<2000)
		maxRange=2000;
	if(maxRange>SEQLEN)
		maxRange=SEQLEN;

	double bgocc=0;

	double threshold=log(PWMThreshold);
	 HashEngine* SearchEngine2=(HashEngine*)SearchEngine;

	 {
		 map<string,double> insts=GenerateInstanceFromPWMPQ(0.8);
		 map<string,double>::iterator Iter;
		double totalCount=0;
		
		for(Iter=insts.begin();Iter!=insts.end();Iter++)
		{	
			char instemp[64];
			string randomIns=Iter->first;

			if(Iter->second<0.05)
			{
				continue;
			}

			strcpy(instemp,randomIns.c_str());
			SearchEngine->searchPattern(instemp,0,poslist);	
		
			if(!LargeDataFlag)
			{
				int len=Iter->first.size();
				double temp=1;
				FOR(j,len)
				{
					if(Iter->first[j]=='N')
						continue;
					int a=acgt(Iter->first[j]);
					if(a>-1)
					temp*=SearchEngine->BGProb[a];

				}
				bgocc+=temp;
			}
		}
	 }


	 POSLIST= poslist;

	 sort(POSLIST.begin(),POSLIST.end());
	 int BINNUM2=SEQLEN/(BindingRegion/2);
	 int centerCnt=0;
	 int bgCnt=0;//new int[BINNUM2];

	 int lastseq=-1;
	 int lastseqL=-1;
	 if(POSLIST.size()==0)
	 {
		 SeqPvalue=-MINSCORE;
		 return POSLIST;
	 }
	 FOR(i,POSLIST.size())
	 {
		 int pos=POSLIST[i]%SEQLEN;
		 int seqnumL=POSLIST[i]/SEQLEN;
		 if(seqnumL!=lastseqL)
		 {
			 lastseq=-1;
			 lastseqL=seqnumL;
		 }
		 int bias=abs(pos-SEQLEN/2);
		 int seqnum=pos/(BindingRegion/2);
		 if(seqnum==lastseq)
			 continue;
		 lastseq=seqnum;
		 double step=1;
		if(SearchEngine->EnableWeight)
			step*=SearchEngine->getSeqWeight(seqnum);
		 if(bias<BindingRegion/2)
		 {
			 centerCnt+=step;
		 }
		 else if(bias<maxRange/2&&bias>=BindingRegion)
			 bgCnt+=step;
		 
	 }
	 int realMAXSEQNUM=0;
	 if(SearchEngine->EnableWeight)
		FOR(i,MAXSEQNUM)
		{
			realMAXSEQNUM+=SearchEngine->getSeqWeight(i);
		}

	 else
		 realMAXSEQNUM=MAXSEQNUM;

	 BINNUM2=BINNUM2*maxRange/SEQLEN;
	
	 double ratio=(double)( bgCnt)/((BINNUM2-4)*realMAXSEQNUM);  //+centerCnt
	 if(ratio>1)
		 ratio=1;
	 if(!LargeDataFlag)
	 {
		 bgCnt=bgocc*((BINNUM2-4)*realMAXSEQNUM)*BindingRegion/2;
		 ratio=bgocc*BindingRegion/2;
	 }
	 SeqPvalue=(double)(centerCnt/2)/((double)bgCnt/(BINNUM2-4));//
	  SeqPvalue=binominalTail(ratio, centerCnt,realMAXSEQNUM*2);	
	 // CDScore=centerCnt;
	  return poslist;

}

void MotifModel::printRealPos(string name)
{
		int i,j,k,q;
	POSLIST.clear();
	name+=".pos";
	name= Setting->outputDIR+"/"+name;
	 FILE* outstream=fopen(name.c_str(), "w+");
	int startpos=Consensus.find_first_not_of("N");
	 if(outstream==NULL)
	 {
		 fclose(outstream);
		 return;
	 }

	map<string,double> inst=GenerateInstanceFromPWMPQ(0.8);
			map<string,double>::iterator Iter1=inst.begin();
		string ACGT="ACGT";
		int lastsize=0;
		int nonConPos=0;
	while(Iter1!=inst.end())
	{
		char instemp2[64];
				strcpy(instemp2,Iter1->first.c_str());
		int rcindex=SearchEngine->searchPattern(instemp2,0,POSLIST);	

		for(q=rcindex;q<POSLIST.size();q++)
			POSLIST[q]=0-POSLIST[q];
		lastsize=POSLIST.size();
		Iter1++;
	}
	FOR(i,POSLIST.size())
	{
		k=abs(POSLIST[i]);
		{
			int seq=k/SEQLEN;
			int pos=k%SEQLEN-SEQLEN/2+startpos;
			
			//if(abs(pos)<(BindingRegion/2))	
			 fprintf(outstream,"%d\t%d\n",seq,pos);
		
		}
	}

	 fclose(outstream);
}
void MotifModel::printMatchPos(string name,vector<VAL>& list)
{
	if(list.size()==0)
		list=getMatchPos();
	int i,j,k;
	name+=".pos";
	name="./MatchPos/"+name;
	 FILE* outstream=fopen(name.c_str(), "w+");
	int startpos=Consensus.find_first_not_of("N");
	 if(outstream==NULL)
	 {
		 fclose(outstream);
		 return;
	 }
  double seqweight=1;
	 FOR(i,list.size())
	{
		k=list[i];

		{
			int seq=k/SEQLEN;
			int pos=k%SEQLEN-SEQLEN/2+startpos;
			 fprintf(outstream,"%d,%d\n",seq,pos);
		
		}
	}
	  fclose(outstream);
}


// enumerate all 'N' elements
	map<string, double> MotifModel::MisMatchClosedSet(string pattern, int mis)
	{
		map<string, double> ret;
		VAL combineNUM=C(pattern.size(),mis);
		double prob=1/(pow((double)4,(double)mis)*combineNUM);
		int i,j,k;
		FOR(i,combineNUM)
		{
			string temp="";
			FOR(j,pattern.size())
				temp.push_back(pattern[j]);
			int p1,left;
			left=i;
			FOR(k,mis)
			{
				p1=left%pattern.size();
				p1+=mis-k;
				temp[p1]='N';
				left/=pattern.size();
			}

			vector<string> nlist=NCloseSet(temp);
			FOR(k, nlist.size())
			{
				ret[nlist[k]]=prob;
			}

		}

		return ret;
	}


	//get a list of seed PWM with good score
 vector<MotifModel*> MotifModel::getSeedMotifs( int SeedNum,int len, int mismatch)
{
	head=(X()-len)/2;
	tail=head; 
	int i,j,k,q;
		int ii;
	cout<<"SeedNum: "<<SeedNum<<"\t len: "<<len<<"\t mismatch: "<<mismatch<<endl;
	LIBSIZE=SearchEngine->getTotalLength();
	string ACGT="ACGT";
	vector<MotifModel*>* Results=new vector<MotifModel*>();
	int loopnum=1<<(2*len);
	
	vector<double> cdscoreList;
	vector<ScoreObj> scoreList;
	
	int BINSIZE2=SEQLEN/2/(BINNUMBER)+1;
	vector<int> patternlist;
	int maxScoreIndex=0;
	double maxScore=MINSCORE;
	double testscore=0;
	ScoreObj testObj;
	int ttt=0;
	int maxRange=1<<(2*Setting->seedlength);
	if(maxRange<2000)
		maxRange=2000;
	
	cdavg=cdvar=intavg=intvar=oravg=orvar=0;orvar=MINSCORE;
	//enumurate all possible
	FOR(i,loopnum)
	{
		string pattern="";
		int hash=i;
		//ignore the reverseComplement
		if(getReverseComplementHashing(hash,len)<i)
		{
			pattern=Hash2ACGT(hash,len);
			pattern.insert(3,string("NN"));
		}
		else
		pattern=Hash2ACGT(hash,len);

		map<string,double> instset;
		instset[pattern]=1;
		double windowsize=1000;//


		char instemp[64];
		strcpy(instemp,pattern.c_str());//
		vector<VAL> positions;
		SearchEngine->searchPattern(instemp,mismatch,positions);
		
		
		double wininc=1.2;
		if(!LargeDataFlag)
			bgfold=1;
		int numwin=log((double)(SEQLEN/Setting->startwSize)/(bgfold))/log(wininc);
		double * CDCounters=new double[numwin];
		double * BGCounters=new double[numwin];

			double* CDHist=new double[BINNUMBER];
	
			FOR(ii,BINNUMBER)
				CDHist[ii]=0;
			int BINSIZE=SEQLEN/BINNUMBER+1;
	
		FOR(ii,numwin)
		{
			CDCounters[ii]=0;
			BGCounters[ii]=0;
		}
		double curScore=-1;
		int coverSeq=0;
		int lastseq=-1;
		
			sort(positions.begin(),positions.end());
		int jj;
		int lastpos=-1;
			FOR(j,positions.size())
			{
				int pos=positions[j]%SEQLEN;
				int seq=positions[j]/SEQLEN;
				if(abs(positions[j]-lastpos)<Setting->seedlength) //zzz
				{
					lastpos=positions[j];
					continue;
				}
				lastpos=positions[j];
				double bias=abs(pos-SEQLEN/2);
					double step=1;
				if(SearchEngine->EnableWeight)
					step*=SearchEngine->getSeqWeight(seq);
		
				if(LargeDataFlag)
				{
					int tt=abs(pos-SEQLEN/2);
					//if(tt<maxActRLen)
					CDHist[tt/BINSIZE2]+=step;
				}
				if(bias>maxRange/2)
					continue;
				if(!LargeDataFlag)
				FOR(jj,numwin)
				{
					int windowsize=(int)Setting->startwSize*pow( wininc,jj);
					
					if(bias<=windowsize/2)
					{
						CDCounters[jj]+=step;					
					}

				}
				if(seq!=lastseq)
				{
					lastseq=seq;
					coverSeq++;
				}
			}
			double maxbgcnt=(double)positions.size()/SEQLEN;
	
			if(coverSeq>(MAXSEQNUM*0.05))// zzz
			{
				int bestwindowId=0;
				
				double bestscore=MINSCORE;
				double bestscore2=MINSCORE;
				int windowsize=0;
				if(LargeDataFlag)
				{
					double sum=0;
					FOR(ii,BINNUMBER)
						sum+=CDHist[ii];
					FOR(ii,BINNUMBER)
						CDHist[ii]=CDHist[ii]/sum;

								double bglen=min(SEQLEN,maxRange);
						int realBINNUMBER=(int)((double)(bglen/SEQLEN)*BINNUMBER);
						FOR(ii,realBINNUMBER/bgfold)
						{
							double cdcnter=0;
							double bgcnter=0;
							
							double score=0;
							FOR(j,ii+1)
								cdcnter+=CDHist[j];
							for(j=ii+1;j< realBINNUMBER;j++)
								bgcnter+=CDHist[j];
							double avgvalue=(double)(bgcnter/( realBINNUMBER-ii-1));
							double ssr=0;
							//for(j=ii+1;j<realBINNUMBER;j++)
							//		ssr+=(CDHist[j]-avgvalue)*(CDHist[j]-avgvalue);
							//ssr=sqrt( (double)(ssr/(realBINNUMBER-ii-1)));
							bgcnter=ssr+avgvalue; // 
							score=(double)((cdcnter/(ii+1))/bgcnter);
							if(score>=bestscore&&cdcnter*sum>Setting->min_supp_ratio*MAXSEQNUM&&bgcnter*sum>bgfold)
							{
								bestwindowId=ii;		
								bestscore=score;
								bestscore2=cdcnter*sum;
							}
						
						}
						windowsize=(bestwindowId+1)*(SEQLEN/BINNUMBER);

				}
				else
				{
					double bgocc=0;
				
					{
						int len=pattern.size();
						double temp=1;
						FOR(j,len)
						{
							if(pattern[j]=='N')
								continue;
							int a=acgt(pattern[j]);
							if(a>-1)
							temp*=SearchEngine->BGProb[a];

						}
						bgocc+=temp;
					}
				FOR(ii,numwin)
				{			
					 windowsize=(Setting->startwSize*pow(wininc,ii));
					 double bglen=min(SEQLEN,maxRange);
					double bgfold=(double)bglen/windowsize;
					 double ratio=(double)(bgfold-1)/bgfold;

					 double score=((CDCounters[ii]+1)/windowsize)/(bgocc*MAXSEQNUM); 

					if(score>bestscore)
					{
						bestwindowId=ii;		
						bestscore=score;
					}
				}

				 windowsize=Setting->startwSize*pow(wininc,bestwindowId);

				}
				//the number of seed depend on the output number
			    if(bestscore<1)
					continue;
				ScoreObj tempobj;
				if(!LargeDataFlag)
				tempobj.cdscore=CDCounters[bestwindowId];
				tempobj.BindingRegion=windowsize;//IntScore;
				tempobj.orscore=bestscore;//still odd ratio

				if(LargeDataFlag)
					tempobj.cdscore=bestscore2;
		
				curScore=bestscore;

				scoreList.push_back(tempobj);
				patternlist.push_back(i);
				cdscoreList.push_back(curScore);
				positions.clear();
		}
		else
		{
				delete[] CDCounters;
				delete[] BGCounters;
				delete[] CDHist;
					continue;
		}
		
		if(maxScore<curScore)
		{
			maxScore=curScore;
			maxScoreIndex=patternlist.size()-1;
		}
		delete[] CDCounters;
		delete[] BGCounters;
		delete[] CDHist;

		
	}
	//caculate the distance matrix
	double r1,r2,br1,br2,step;
	int bestRank=1000000;

	//statistics var
	cdavg=cdavg/scoreList.size();

	intavg=intavg/scoreList.size();

	oravg=oravg/scoreList.size();

	maxScore=MINSCORE;


	r1=r2=0;


	cout<<"Rank: "<<bestRank<<endl;
	int motifid=1;
	FOR(i,SeedNum)
	{
		if(cdscoreList[maxScoreIndex]==MINSCORE)
		{
			cout<<"early stop!"<<endl;
			break;
		}
		MotifModel* im=new MotifModel(SearchEngine,X(),Setting);

		string fn="/AlignmentResult/"+Hash2ACGT(patternlist[ maxScoreIndex],len)+".fa";
		
		string topPattern=Hash2ACGT(patternlist[ maxScoreIndex],len);
		int hash=patternlist[ maxScoreIndex];
		if(getReverseComplementHashing(hash,len)<hash)
		{
			topPattern.insert(3,string("NN"));
		}
		cout<<"top:"<<motifid-1<<"\t"<<topPattern<<":"<<cdscoreList[maxScoreIndex]<<"\t"<<scoreList[maxScoreIndex].cdscore<<"\t"<<scoreList[maxScoreIndex].BindingRegion<<"\t"<<scoreList[maxScoreIndex].orscore<<endl;

		im->AddInstance(printAlignment(topPattern,0,len));
		vector<double> simScorelist;
		vector<int> rank;
		int count=0;
		cdscoreList[maxScoreIndex]=MINSCORE;

		im->CDScore=scoreList[maxScoreIndex].cdscore;
		im->ORScore=scoreList[maxScoreIndex].orscore;
		im->BindingRegion=scoreList[maxScoreIndex].BindingRegion;
		im->cdavg=cdavg;
		im->cdvar=cdvar;
		im->intavg=intavg;
		im->intvar=intvar;
		im->oravg=oravg;
		im->orvar=orvar;
		im->seed=topPattern;
		im->Setting=Setting;
		Results->push_back(im);
		


		maxScore=MINSCORE;
		motifid++;
	FOR(j,cdscoreList.size())
		{
			//merge the similar k-mer
			if(cdscoreList[j]==MINSCORE)
				continue;
			//string tempPattern=Hash2ACGT(patternlist[j],len);
			//int hash=patternlist[j];
			//if(getReverseComplementHashing(hash,len)<hash)
			//{
			//	tempPattern.insert(3,string("NN"));
			//}
			//int alnScore=0;
			//int aln=AlignmentRC(topPattern,tempPattern,alnScore,1);
			//if(aln==0)
			//{
			//	MotifModel* imMerge=new MotifModel(SearchEngine,X(),Setting);
			//	imMerge->AddInstance(printAlignment(tempPattern,0,len));
			//}
			//else 

			if(maxScore<cdscoreList[j])
			{
				maxScore=cdscoreList[j];
				maxScoreIndex=j;
			}
		}

		im->InitializePWMofInstanceSet();
		
		cout<<im->get_consensus(0)<<endl;
		
	}


	return *Results;

}

// Add Instance to  instanceSet and update the PWM
void MotifModel::AddInstance(string instanceString)
{
	InstanceSet.push_back(instanceString);
}

void MotifModel::MergeList(vector<MotifModel*> simList)
{
	int i,j;
	MotifModel* temp=new MotifModel(this->SearchEngine,this->X(),this->Setting);
	temp->head=this->head;
	temp->tail=this->tail;
	FOR(i,this->X())
	{
		//copy first
		FOR(j,4)
			temp->s(this->g(i,j),i,j);
		//boost up by weight
		FOR(j,4)
			this->m(this->ORScore-1,i,j);
	}
	temp->get_consensus(0);
	//use temp to align all
	FOR(i,simList.size())
	{
			MotifModel* p1=temp;
		
		MotifModel* p2=simList[i];
		
		 double divg = 10;
            double temp;
			int overlap=0;
			int bestol=overlap;
          int bestaln=0;
            int aln=AlignmentPWM(p1, p2, temp,overlap);
            divg = temp;
			bestaln=aln;
			bestol=overlap;
			MotifModel* p2RC=&p2->RC();
            aln = AlignmentPWM(p1, p2RC, temp,overlap);
            if (temp < divg)
			{
                divg = temp;
				bestaln=aln;
				bestol=overlap;
				int pp=0;
				int i=bestaln;
                int j;
				if (i > 0)
				{	pp = i;
				    if(p2->Length()>bestol)
					{
						//bestol=p2->Length();
						bestol=p1->X()-p1->head-pp;//to right end of wide pwm p1
					}
				}
				else
				{
					//pp=i;
					//bestol=p2->Length();
					pp= -p1->head;
					bestol=p1->head+i+p2->Length(); //to the left end of wide pwm p1
					if((i-pp)>p2->head)
					{
						pp=i-p2->head;
						bestol=p2->head+p2->Length();
					}
				}
                
	            for(j=pp;j-pp<bestol;j++)
	            {
			         double sump = 0;
                    for (int p = 0; p < 4; p++)
                    {
						this->a(p2RC->g(j-i+p2RC->head,p)*(p2->ORScore-1),j+this->head,p);
                    }
             
	            }
			}
			else
			{
				int pp=0;
				int i=bestaln;
                
				 int j;
				if (i > 0)
				{	pp = i;
				    if(p2->Length()>bestol)
					{
						//bestol=p2->Length();
						bestol=p1->X()-p1->head-pp;//to right end of wide pwm p1
					}
				}
				else
				{
					//pp=i;
					//bestol=p2->Length();
					
					pp= -p1->head;
					bestol=p1->head+i+p2->Length(); //to the left end of wide pwm p1
					if((i-pp)>p2->head)
					{
						pp=i-p2->head;
						bestol=p2->head+p2->Length();
					}
				}
				for(j=pp;j-pp<bestol;j++)
	            {
			         double sump = 0;
                    for (int p = 0; p < 4; p++)
                    {
						this->a(p2->g(j-i+p2->head,p)*(p2->ORScore-1),j+this->head,p);
                    }
             
	            }

			}    
	}

	//normalized
	FOR(i,this->X())
	{double sum=0;
		FOR(j,4)
			sum+=this->g(i,j);
		FOR(j,4)
			this->d(sum,i,j);
	}

	this->get_consensus(0);


}

// merge two PWM using PWM algniment
void MotifModel::Merge(MotifModel* P)
{
			MotifModel* p1=this;

		string c1=this->get_consensus(0);
		string c2=P->get_consensus(0);
		
		MotifModel* p2=P;
		 double divg = 10;
            double temp;
			int overlap=0;
			int bestol=overlap;
          int bestaln=0;
            int aln=AlignmentPWM(p1, p2, temp,overlap);
            divg = temp;
			bestaln=aln;
			bestol=overlap;
			MotifModel* p2RC=&p2->RC();
            aln = AlignmentPWM(p1, p2RC, temp,overlap);
            if (temp < divg)
			{
                divg = temp;
				bestaln=aln;
				bestol=overlap;
				int pp=0;
				int i=bestaln;
                int j;
				if (i > 0)
				{	pp = i;
				    if(p2->Length()>bestol)
						bestol=p2->Length();
				}
				else
				{
					pp=i;
					bestol=p2->Length()+i;
				}
                
	            for(j=pp;j-pp<bestol;j++)
	            {
			         double sump = 0;
                    for (int p = 0; p < 4; p++)
                    {
						p1->a(p2RC->g(j-i+p2RC->head,p)*(p2->ORScore-1),j+p1->head,p);
                    }
             
	            }
			}
			else
			{
				int pp=0;
				int i=bestaln;
                
				 int j;
				if (i > 0)
				{	pp = i;
				    if(p2->Length()>bestol)
						bestol=p2->Length();
				}
				else
				{
					pp=i;
					bestol=p2->Length()+i;
				}
				for(j=pp;j-pp<bestol;j++)
	            {
			         double sump = 0;
                    for (int p = 0; p < 4; p++)
                    {
						p1->a(p2->g(j-i+p2->head,p)*(p2->ORScore-1),j+p1->head,p);
                    }
             
	            }

			}    
           
}

// compute serveral scores of this motif model
void MotifModel::ComputeScore(double sampleratio, double & CDScore, double & ORScore, double & BindingRegion, double & CNSVScore, double & DiffScore)
{
	int i,j,k;
	//CDScore=ORScore=1; 
	if(!LargeDataFlag)
		bgfold=1;
	int numwin=log((double)(SEQLEN/Setting->startwSize)/(bgfold))/log(wininc);
	double * CDCounters=new double[numwin];
	double * BGCounters=new double[numwin];

	FOR(i,numwin)
	{
		CDCounters[i]=0;
		BGCounters[i]=0;
	}
	map<string,double> insts=GenerateInstanceFromPWMPQ(sampleratio);
	
	if(insts.size()==0)
	{
		CDScore=ORScore=CNSVScore=MINSCORE;
		this->CDScore=this->ORScore=MINSCORE;
			delete [] CDCounters;
			delete [] BGCounters;

		return;
	}

	if(DEBUG)
	cout<<insts.size()<<endl;
	vector<VAL> positions;
	map<string,double>::iterator Iter;
	double totalCount=0;
	double BGcount=0;
	double CDcount=0;
	vector<int> countlist;
	double sumNorm=0;
	POSLIST.clear();
	double coverSeq=0;
	int lastposSize=0;

	double* CDHist=new double[BINNUMBER];
	FOR(i,BINNUMBER)
		CDHist[i]=0;
	int BINSIZE=SEQLEN/BINNUMBER+1;
	

	int maxRange=1<<(2*Length());
	if(maxRange<2000)
		maxRange=2000;
	if(maxRange>SEQLEN)
	maxRange=SEQLEN;
int BINSIZE2=SEQLEN/2/(BINNUMBER)+1;
	int maxCover=0;
	int countNonX=0;
	for(Iter=insts.begin();Iter!=insts.end();Iter++)
	{	
		//coverSeq=0; 
		char instemp[64];
		string randomIns=Iter->first;
		strcpy(instemp,randomIns.c_str());
		int rcindex=SearchEngine->searchPattern(instemp,0,positions);	
		
		int lastseq=-1;
		int lastpos=-1;
		
		double probw=(double)(positions.size()-lastposSize+1)/(LIBSIZE);

		
			
		
		int totalsize=positions.size();
		double step=Iter->second/probw;
		int ct=Length()/2;
		
			for(i=lastposSize;i<totalsize;i++)
			{
				
				long int wpos=positions[i];
				

				
				int seqnum=wpos/SEQLEN;
				int pos=wpos%SEQLEN;
				if(seqnum!=lastseq)
				{
					lastseq=seqnum;
					coverSeq+=1;	
					lastpos=-1;
				}
				if(SearchEngine->CharText[wpos+ct]=='X'||wpos-lastpos<ct*3)//zzz
				{
					lastpos=wpos;
					continue;
				}
				lastpos=wpos;
				double seqweight=1;
				if(SearchEngine->EnableWeight)
					seqweight=SearchEngine->getSeqWeight(seqnum);

				if(LargeDataFlag)
				{
					int tt=abs(pos-SEQLEN/2);
						
						  CDHist[tt/BINSIZE2]+=step*seqweight;
				}
				/////////
				//posflag[pos]=1;
				sumNorm+=step*seqweight;		


				countNonX++;

				double bias=abs(pos-SEQLEN/2);	
				if(bias>maxRange/2)
					continue;
				if(i<rcindex)
				{
				   POSLIST.push_back(wpos);
				}
				else
				{
				   POSLIST.push_back(0-wpos);
				}
				
				if(!LargeDataFlag)
				FOR(j,numwin)
				{
					int windowsize=Setting->startwSize*pow( wininc,j);
					
					if(bias<=windowsize/2)
					{
						CDCounters[j]+=step*seqweight;
					
					}

				}

				
				
				
				
			}
			if(maxCover<coverSeq)
				maxCover=coverSeq;
			
			lastposSize=totalsize;

	}

	double nonblkRatio=1;//(double)(nonBLK)/SEQLEN;

	int bestwindowId=-1;
	double bestscore=MINSCORE;
	int windowsize=0;
	if(positions.size()!=0)
		sumNorm=sumNorm/(double)positions.size();
	else
		sumNorm=1;
	
	double maxbgcnt=(double)positions.size()/SEQLEN;

	double expbgcnt=pow(0.25,Efflen)*2*LIBSIZE/BINNUMBER/Length()*insts.size();

	if(maxCover>Setting->min_supp_ratio *MAXSEQNUM)//zzz
	{
		
		if(LargeDataFlag)
		{
			//double sum=0;
			//FOR(i,BINNUMBER)
			//	sum+=CDHist[i];
			FOR(i,BINNUMBER)
				CDHist[i]=CDHist[i]/sumNorm;
			double bglen=min(SEQLEN,maxRange);
			double bestcdcnt,bestbgcnt;
			int realBINNUMBER=(int)((double)(bglen/SEQLEN)*BINNUMBER);
			FOR(i,realBINNUMBER/bgfold)
			{
				double cdcnter=0;
				double bgcnter=0;
				
				double score=0;
				FOR(j,i+1)
					cdcnter+=CDHist[j];
				for(j=i+1;j<realBINNUMBER;j++)
					bgcnter+=CDHist[j];
				//only average on bg
				double avgvalue=(double)(bgcnter/(realBINNUMBER-i-1));
				double ssr=0;
				for(j=i+1;j<realBINNUMBER-1;j++)
						ssr+=(CDHist[j+1]-CDHist[j])*(CDHist[j+1]-CDHist[j]);
				ssr=sqrt( (double)(ssr/(realBINNUMBER-i-2)));
				bgcnter=ssr+avgvalue; //1.95996*
				if(cdcnter<Setting->min_supp_ratio*MAXSEQNUM||bgcnter<expbgcnt)//zzz
					continue;
				score=(double)((cdcnter/(i+1))/bgcnter);
				bgcnter-=ssr;  //1.95996*
				bgcnter*=(realBINNUMBER-i-1);
				windowsize=(i+1)*(SEQLEN/BINNUMBER);
				double bgfold=(double)bglen/windowsize;
				if(ORScore<1)
					ORScore=1;
				 double ratio=(double)(ORScore)/(ORScore+bgfold-1);//(double)(bgfold-1)/bgfold;//
		//zzz	//double pvalue=binominalTail(ratio,cdcnter,(bgcnter+cdcnter))*(CNSVScore-DiffScore);//*pow((double)4,Length()-Setting->seedlength);
				// if(pvalue>0.01)//BGCounters[i]<0.05*(CDCounters[i]+BGCounters[i])||
				//	 continue;
				if(score>=bestscore)
				{
					bestwindowId=i;		
					bestscore=score;
					CDScore=cdcnter;
					bestbgcnt=bgcnter;
					bestcdcnt=cdcnter;
				}
			
			}
			if(bestscore!=MINSCORE)
			{
			windowsize=(bestwindowId+1)*(SEQLEN/BINNUMBER);
			double ratio=((double)(realBINNUMBER-bestwindowId-1)/realBINNUMBER);

			//BinPvalue=binominalCDF(ratio,bestbgcnt,(bestbgcnt+bestcdcnt))*maxRange;
			}
		}
		else
		{
			double bgocc=0;
			for(Iter=insts.begin();Iter!=insts.end();Iter++)
			{
				int len=Iter->first.size();
				double temp=1;
				FOR(j,len)
				{
					if(Iter->first[j]=='N')
						continue;
					int a=acgt(Iter->first[j]);
					if(a>-1)
					temp*=SearchEngine->BGProb[a];

				}
				bgocc+=temp;
			}
			FOR(i,numwin)
			{
				CDCounters[i]/=sumNorm;		
				if(CDCounters[i]<Setting->min_supp_ratio*MAXSEQNUM) //zzz
					continue;			
				 windowsize=Setting->startwSize*pow(wininc,i);
				  double bglen=min(SEQLEN,maxRange);
				double bgfold=(double)bglen/windowsize;
				if(ORScore<1)
					ORScore=1;
				 windowsize=Setting->startwSize*pow(wininc,i);
				 double ratio=(double)(bgfold-1)/(ORScore+bgfold-1);
				double score=((CDCounters[i])/ windowsize)/(bgocc*MAXSEQNUM); 
				if(score>=bestscore)
				{
					bestwindowId=i;		
					bestscore=score;
				}
			}
			 windowsize=Setting->startwSize*pow(wininc,bestwindowId);
		}
	}
	
	if(bestscore==MINSCORE)   
	{
		CDScore=ORScore=CNSVScore=MINSCORE;
		this->CDScore=this->ORScore=MINSCORE;
	}
	else
	{
		if(!LargeDataFlag)
		CDScore=CDCounters[bestwindowId];
		BindingRegion=windowsize; 
		ORScore=bestscore;

		this->ORScore=ORScore;
		this->CDScore=CDScore;
		this->BindingRegion=BindingRegion;
		//double cCap=CDCounters[bestwindowId]/(CDCounters[bestwindowId]+BGCounters[bestwindowId]);
		//	double c0=(double)windowsize/SEQLEN;
		//if(LargeDataFlag)
		//	CDScore=cCap*log(cCap/c0)+(1-cCap)*log((1-cCap)/(1-c0));
	}


	//delete []CDHistogram;
	//delete []INTHistogram;
	delete [] CDCounters;
	delete [] BGCounters;
	delete [] CDHist;

	
}

//return a instance and its probabilty
pair<string,double> MotifModel::randomInstance()
{
	int i,j,k;
	string ACGT="ACGT";
	int totalN=0;
	double ratio=1;
		string instance="";
		FOR(j,X())
		{
			double entropy=0;
			
			entropy=this->entropy2(j);

			if(entropy>ENTROPY_Threshold)
			{
				instance.push_back('N');
				totalN++;
				
			}
			else
			{
				double pt=(double)(rand()%1000)/(double)1000;
				double sum=0;
				FOR(k,4)
				{
					sum+=g(j,k);
					if(sum>pt)
					{
						instance.push_back(ACGT[k]);
						ratio*=g(j,k);
						break;
					}
				}
			}

			
		}

		//trimme N in the head and tail

		head=tail=0;
		FOR(j,X())
		{
			if(instance[j]=='N')
				head++;
			else
				break;
		}
		FOR(j,X())
		{
			if(instance[X()-j-1]=='N')
				tail++;
			else
				break;
		}
		pair<string,double> ret;
		ret.first=(instance.substr(head,instance.size()-tail-head));
		//vector<string> nset=NCloseSet(ret.first);
		//ret.first=nset[rand()%nset.size()];
		ret.second=ratio;//*pow(0.25,totalN-tail-head);
	
		return ret;
	
}

//set a PWM to trivial PWM
void MotifModel::resetPWM()
{
	int i,j;
	int start=(X()-Setting->seedlength)/2;
	FOR(i,X())
	{
		if(i>=start&&i<start+Setting->seedlength)
			continue;
		FOR(j,4)
			s(0.25,i,j);
	}
	lastupdateCol=-1;
	badmove.clear();
	CDScore=ORScore=CNSVScore=DiffScore=0;

}

// MCMC method generate instance from this motif model
map<string,double> MotifModel::GenerateInstanceFromPWM(double sampleratio,bool noThreshold)
{
	int i,j,k;
	string ACGT="ACGT";
	map<string,double> retSet;
	double sumRatio=0;
	while(sumRatio<sampleratio)
	{
		pair<string,double> temp=randomInstance();
		//ignore the low probability
		
			
		//if(temp.second<1/(double)(X()*4))
		//	continue;
		map<string,double>::iterator iter =retSet.find(temp.first);
		if(iter==retSet.end())
		{
			sumRatio+=temp.second;
			if(temp.second<PWMThreshold&&!noThreshold)
			continue;
		
			retSet[temp.first]=temp.second;
		}
	}
	map<string,double>::iterator Iter;
	for(Iter=retSet.begin();Iter!=retSet.end();Iter++)
	{
		Iter->second=Iter->second/sumRatio;
	}
	

	return retSet;	
}


// Exact Priority Queue method generate instance from this motif model
map<string,double> MotifModel::GenerateInstanceFromPWMPQ(double sampleratio,bool noThreshold)
{
	int i,j,k;
	string ACGT="ACGT";
	 int totalN=0;
	double ratio=1;
	vector<int> effIndex;
	map<string,double> retSet;
	map<double,unsigned long long int> PQ;
	vector< vector<int> > orderedCol;
	set<unsigned long long int> hashtable;
	double sumRatio=0;
	string instance="";

	FOR(j,X())
	{
		double entropy=0;
		
		entropy=this->entropy2(j);

		if(entropy<ENTROPY_Threshold)
		{
			effIndex.push_back(j);	
			orderedCol.push_back(columnOrder(j));
		}
	}
	int effLen=effIndex.size();
	head=effIndex[0];

	tail=X()-effIndex[effLen-1]-1;
	unsigned long long int tophashcode=0;
	double probv=1;
	Efflen=effLen;
	FOR(i, effLen)
	{
		tophashcode<<=2;
		tophashcode+=orderedCol[i][0];
		probv*=g(effIndex[i],orderedCol[i][0]);	
	}
	const double smallscale=0.0000000000001;
	PQ[sumRatio*smallscale-probv]=tophashcode;//make the key little different
	hashtable.insert(tophashcode);
	int maxsize=1<<(int)(effLen-Setting->seedlength+2);
	maxsize=min(maxsize,MAXSEQNUM);
	while(sumRatio<sampleratio)
	{
		if(PQ.size()==0)
			break;
		sumRatio+=0-PQ.begin()->first;
		tophashcode=PQ.begin()->second;
		vector<unsigned long long int> tempptn;
		//decode
			string topptn="";
			FOR(k, effLen)
			{
				if(k>0)
				{
					FOR(i,effIndex[effLen-k]-effIndex[effLen-k-1]-1)
						topptn="N"+topptn;

				}
				int acgtid=tophashcode%4;
				topptn=ACGT[acgtid]+topptn;
				tempptn.push_back(acgtid);
				tophashcode>>=2;
			}
			probv=0-PQ.begin()->first;
			retSet[topptn]=probv;
			//hashtable.insert(topptn);// add to hashtable
			PQ.erase(PQ.begin()); //dequeue
			if(sumRatio>sampleratio||(probv<PWMThreshold&&!noThreshold)||(retSet.size()>=maxsize))//&&!noThreshold
				break;
		//mutate 
			FOR(i,effLen)
			{
				//left to right 
				int dd=effLen-i-1;
				int key=tempptn[dd];
				FOR(j,4)
					if(key==orderedCol[i][j])
						break;
				if(j<3)
				{
					unsigned long long int temphashcode=0;
					double tempprobv=1;
					bool nonsenseFlag=false;
					//left to right 
					FOR(k,effLen)
					{
						int cc=effLen-k-1;
						temphashcode<<=2;
						if(k!=i)
						{
							temphashcode+=tempptn[cc];
							double ti=g(effIndex[k],tempptn[cc]);
							if(!noThreshold&&ti<0.01)
							{
									nonsenseFlag=true;
									break;
							}

								tempprobv*=ti;
						}
						else
						{
							int mutcode=orderedCol[k][j+1];
							temphashcode+=mutcode;
							double ti=g(effIndex[k],mutcode);
							/////omit the element with prob less than 0.01
							if(!noThreshold&&ti<0.01)
							{
									nonsenseFlag=true;
									break;
							}

								tempprobv*=ti;
						}
					}
					if(hashtable.find(temphashcode)!=hashtable.end()||nonsenseFlag)
						continue;
					hashtable.insert(temphashcode);
					PQ[hashtable.size()*smallscale-tempprobv]=temphashcode;
				}
			}
				
	}
	map<string, double>::iterator iter1=retSet.begin();
	while(iter1!=retSet.end())
	{
		iter1->second/=sumRatio;
		iter1++;
	}

	return retSet;	
}

//sort the column element from high prob to low prob
vector<int> MotifModel::columnOrder(int col)
{
		int i,j,k;
		vector<int> ret;
		map<double,int> sortlist;
		double smallvalue=0.000000000001;
		FOR(i,4)
		{
			sortlist[smallvalue*i-g(col,i)]=i;
		}
		map<double,int>::iterator iter=sortlist.begin();
		while(iter!=sortlist.end())
		{
			ret.push_back(iter->second);
			iter++;
		}
		return ret;
}


// Compute center distributioin score
//double MotifModel::ComputeCDScore(double* histgramBIN)
//{
//	int binsize=SEQLEN/BINNUMBER+1;
//	
//	double sum=0;
//		
//		for (int i = 0; i <  BINNUMBER; ++i)
//		{
//				sum+=histgramBIN[i];
//		}
//		for (int i = 1; i < BINNUMBER; ++i)
//		{
//			histgramBIN[i]+=histgramBIN[i-1];
//		}
//		bool flag=true;
//		
//		double sum2=0;
//		double tems1=0;
//		double tems2=0;
//		for (int i = 0; i < BINNUMBER; ++i)
//		{			
//			double h=(double)(i+1)/BINNUMBER;
//			double temp=(h-(double)histgramBIN[i]/sum);
//			if(i<=BINNUMBER/2)
//				sum2+=temp;
//			else
//				sum2-=temp;
//
//		}
//		double score=sum2;
//		return score;
//}



// Compute Intensity Score
//double MotifModel::ComputeBindingRegion(double* histgramBIN)
//{
//		int binsize=MAXSEQNUM/BINNUMBER+1;
//	    double sum=0;
//		
//		for (int i = 0; i <  BINNUMBER; ++i)
//		{
//				sum+=histgramBIN[i];
//		}
//		for (int i = 1; i < BINNUMBER; ++i)
//		{
//			histgramBIN[i]+=histgramBIN[i-1];
//		}
//		bool flag=true;
//		
//		double sum2=0;
//		double tems1=0;
//		double tems2=0;
//		for (int i = 0; i < BINNUMBER; ++i)
//		{			
//			double h=(double)(i+1)/BINNUMBER;
//			double temp=(h-(double)histgramBIN[i]/sum);
//				sum2-=temp;
//		}
//		double score=sum2;
//		return abs(score);
//
//}

vector<string> MotifModel::NCloseSet(std::string pattern)
{
	string ACGT="ACGT";
	int i,j,k;
	vector<int> pos;
	vector<string> retSet;
	FOR(i,pattern.size())
	{
		if(pattern[i]=='N')
		{
			pos.push_back(i);
		}
	}

	if(pos.size()==0)
	{
		retSet.push_back(pattern);
		return retSet;
	}
	int loop=1<<(2*pos.size());

	FOR(i,loop)
	{
		string temp=pattern;
		FOR(j,pos.size())
		{
			
			int code=(i>>(j*2))%4;
			temp[pos[j]]=ACGT[code];

		}
		retSet.push_back(temp);
	}

	

	return retSet;
}

// Compute OverRepresented Score
//double MotifModel::ComputeORScore(map<string,double>& InstSet, vector<int> countlist,bool mflag)
//{
//	double ProfM=0;
//	double markovScore=0;
//
//	map<string,double>::iterator Iter=InstSet.begin();
//	int Setting->seedlength=Iter->first.size();
//	double score=0;
//	double score2=0;
//	int i=0;
//	double delta=(double)1/Setting->seedlength;
//	while(Iter!=InstSet.end())
//	{
//		if(countlist[i]==0)
//		{
//			Iter++;
//			i++;
//			continue;
//		}
//		//ProfM+=bgmodel.getStringProb(Iter->first)*(Iter->second);
//		ProfM=bgmodel.getStringProb(Iter->first);
//		double expect=MAXSEQNUM*ProfM;//(LIBSIZE-Setting->seedlength)*ProfM;
//		if(mflag)
//		markovScore+=log(ProfM)*(Iter->second);
//		//score+=log((countlist[i]+delta)/(expect+delta))*(Iter->second);
//		score+=(countlist[i])*(Iter->second);
//		score2+=(expect)*(Iter->second);
//		//double P=(double)countlist[i]/LIBSIZE;
//		//score+=countlist[i]*(log(P/(1-P))-log(ProfM)+ComputeCNSVScore())*(Iter->second);
//		Iter++;
//		i++;
//	}
//	
//	//double score=count/ProfM-(LIBSIZE-Setting->seedlength);
//	//score*=abs(score);
//	//score+=1;
//	//score/=((LIBSIZE-Setting->seedlength+1)/ProfM+((LIBSIZE-Setting->seedlength*2+2)*(LIBSIZE-Setting->seedlength*2+1)-(LIBSIZE-Setting->seedlength+1)*(LIBSIZE-Setting->seedlength+1)));
//
//	 //score=(count+1)/((LIBSIZE-Setting->seedlength)*ProfM+1);
//	return (score)/(score2); //exp(score-markovScore);
//}

// Compute Conservation Score
double MotifModel::ComputeCNSVScore(void)
{
	double score=0;
	double sum=0;
	for (int i = head; i <X()-tail; ++i)
	{	
		double sum=0;
		for(int j=0;j<4;j++)
		{
			sum+=g(i,j);
		}

		for(int j=0;j<4;j++)
		{
			double temp=g(i,j)/sum+0.0000000000000001;
			score+=temp*log(temp);
		}

	}

	score=score;
	return score;
}

// get summary score of this motif model
double MotifModel::GetMixedScore(void)
{
	return GetMixedScore( CDScore, ORScore,BindingRegion, CNSVScore,  DiffScore);
}


string MotifModel::get_consensus(double pseudo_weight , double genome_fra , double genome_frc)
{
	int i,j,k,l;
    string consensus;
    string consensus_pattern;

    double total_freq = 0;
    double matrix_value = 0;
    double matrix_weight = 0;
    double nuc_freq = 0;
    double prior_prob = 0;
    
    for(i = 0; i < X(); ++i) {
	total_freq = 0;
	consensus_pattern = "";

	FOR(j,Y()) {
	    FOR(k,Z()) {
		FOR(l,W()) {
		    total_freq += double(g(i,j,k,l));
		}
	    }
	}

	FOR(j,Y()) {
	    nuc_freq = 0;
    	    FOR(k,Z()) {
		FOR(l,W()) {
		    nuc_freq += double(g(i,j,k,l));
		}
	    }

	    if ((i == 0) || (i == 3)) { prior_prob = genome_fra / 2; }
	    if ((i == 1) || (i == 2)) { prior_prob = genome_frc / 2; }
	    matrix_value = ((nuc_freq / total_freq) + (pseudo_weight * prior_prob)) / (1 + pseudo_weight);
	    matrix_weight = log(matrix_value/ prior_prob);
	
	    if((matrix_weight > 0) && (j == 0)) { consensus_pattern = consensus_pattern + "A"; }
	    if((matrix_weight > 0) && (j == 1)) { consensus_pattern = consensus_pattern + "C"; }
	    if((matrix_weight > 0) && (j == 2)) { consensus_pattern = consensus_pattern + "G"; }
	    if((matrix_weight > 0) && (j == 3)) { consensus_pattern = consensus_pattern + "T"; }		    
	}

	if(entropy2(i)>ENTROPY_Threshold)
		{ consensus = consensus + "N"; }	
	else if (consensus_pattern == "A" ) { consensus = consensus + "A"; }
	else if (consensus_pattern == "C" )	{ consensus = consensus + "C"; }
	else if (consensus_pattern == "G" )	{ consensus = consensus + "G"; }
	else if (consensus_pattern == "T" )	{ consensus = consensus + "T"; }
	else if (consensus_pattern == "AC" )	{ consensus = consensus + "M"; }
	else if (consensus_pattern == "AG" )	{ consensus = consensus + "R"; }
	else if (consensus_pattern == "AT" )	{ consensus = consensus + "W"; }
	else if (consensus_pattern == "CG" )	{ consensus = consensus + "S"; }
	else if (consensus_pattern == "CT" )	{ consensus = consensus + "Y"; }
	else if (consensus_pattern == "GT" )	{ consensus = consensus + "K"; }
	else if (consensus_pattern == "ACG" )	{ consensus = consensus + "V"; }
	else if (consensus_pattern == "ACT" )	{ consensus = consensus + "H"; }
	else if (consensus_pattern == "AGT" )	{ consensus = consensus + "D"; }
	else if (consensus_pattern == "CGT" )	{ consensus = consensus + "B"; }
	else					{ consensus = consensus + "N"; }
    }

   Consensus=consensus;
		head=tail=0;
	FOR(j,X())
	{
		if(Consensus[j]=='N')
			head++;
		else
			break;
	}
	FOR(j,X())
	{
		int k=X()-j-1;
		if(Consensus[k]=='N')
			tail++;
		else
			break;
	}
	Consensus=(consensus.substr(head,X()-head-tail));
		return Consensus;
};

void MotifModel::PWMRefinement()
{
	string temp=get_consensus(0);
	int t=3;
	int i,j;
	for(i=1;i<temp.size();i++)
	{
		if(temp[i]!='N')
			break;
	}
	
	for(j=1;j<temp.size();j++)
	{
		if(temp[temp.size()-j-1]!='N')
			break;
	}

	if(j>t)
		tail+=j;
	
	if(i>t)
		head+=i;
	if(j>t||i>t)
	{
		FOR(i,X())
		{
			if(i<head||i>(X()-tail) )
			{
				FOR(j,4)
				{
					s(0.25,i,j);
				}
			}
		}
		DiffScore=Setting->seedlength;
		CNSVScore=Setting->seedlength;
		
	}
	ComputeScore(0.8,CDScore,ORScore,BindingRegion,CNSVScore,DiffScore);
	int len=Length();
	int ct=len/2;
	InstanceSet.clear();
	FOR(i,POSLIST.size())
	{
		long int wpos=POSLIST[i];

		bool rc=(wpos<0);
		VAL upos=wpos;
		if(rc)
		{
			upos=0-wpos;
			POSLIST[i]=upos;
		}
		int seqnum=upos/SEQLEN;
		int pos=upos%SEQLEN;
		if(SearchEngine->CharText[upos+ct]=='X')
					continue;
			int windowsize=BindingRegion;
			double bias=abs(pos-SEQLEN/2);
			if(bias<windowsize/2)
			{
					string pa;
					if(rc)
					{
						pa=SearchEngine->getSite(upos,len);
			
						pa=reverseString(pa);
					}
					else
					{
						pa=SearchEngine->getSite(upos,len);
				
					}
				
					InstanceSet.push_back(pa);
			}

	}
	
	InitializePWMofInstanceSet();

}

// update the PWM by random search , and runid specify the column will be updated
double MotifModel::updateModel(int runid)
{

	//_CrtDumpMemoryLeaks();
	bool reviseFlag=false;

	if(runid==0)
	{
		DiffScore=Setting->seedlength;
		CNSVScore=Setting->seedlength;
		if(this->Length()>8)
			this->print();
		ComputeScore(0.8,CDScore,ORScore,BindingRegion,CNSVScore,DiffScore);
		bgfold=5; //zzz
		return 0;
	}
	MotifModel tempmodel(SearchEngine,X(),Setting);
	MotifModel tempBGmodel(SearchEngine,X(),Setting);
	int width=X();


	int col=(width+(runid%2)-1)/2+((runid%2)-0.5)*runid;
	

	int i,j,k,q;
	double newScore=0;
	double oldScore=GetMixedScore();
	
Consensus=get_consensus(0);

int rcbias=Length()/2+head;
double binomailP=(double)POSLIST.size()/LIBSIZE;
double expX=(double)POSLIST.size()/Length();
int Xcount=0;
if(badmove.size()>0)
{
	tempmodel.InstanceSet=InstanceSet;
	tempBGmodel.InstanceSet=BGInstanceSet;

}
else
{

	FOR(i,POSLIST.size())
	{
		long int wpos=POSLIST[i];

		bool rc=(wpos<0);
		VAL upos=wpos;
		if(rc)
		{
			upos=0-wpos;
		
		}
		int seqnum=upos/SEQLEN;
		int pos=upos%SEQLEN;
		
			int windowsize=BindingRegion;
			double bias=abs(pos-SEQLEN/2);
			
			if(bias>SEQLEN/bgfold&&LargeDataFlag)//if(bias>windowsize/2)
			{
				//if(bias<(bgfold/2+1)*windowsize||LargeDataFlag)//(bgfold/2)*windowsize)
				if(tempBGmodel.InstanceSet.size()>bgfold*tempmodel.InstanceSet.size())
					continue;
				{
					string pa;
					if(rc)
					{
						pa=SearchEngine->getSite(upos-tail,X());					
						pa=reverseString(pa);
					}
					else
					{
						pa=SearchEngine->getSite(upos-head,X());
			
					}

					tempBGmodel.AddInstance(pa);
					if(SearchEngine->EnableWeight)
					tempBGmodel.instSeq.push_back(seqnum);
				}
			}
			else if(bias<windowsize/2)
			{
					string pa;
					if(rc)
					{
						pa=SearchEngine->getSite(upos-tail,X());
			
						pa=reverseString(pa);
					}
					else
					{
						pa=SearchEngine->getSite(upos-head,X());
				
					}
				
				tempmodel.AddInstance(pa);
				if(SearchEngine->EnableWeight)
				tempmodel.instSeq.push_back(seqnum);
			}

	}
	InstanceSet=tempmodel.InstanceSet;
		BGInstanceSet=tempBGmodel.InstanceSet;
}
	if(LargeDataFlag)
		if(tempBGmodel.InstanceSet.size()==0||tempmodel.InstanceSet.size()==0)
		{
			if(tempBGmodel.InstanceSet.size()==0)
				bgfold*=2;
			return -1;
		}

	tempmodel.head=0;
	tempmodel.InitializePWMofInstanceSet();
	tempBGmodel.head=0;
	if(LargeDataFlag)
	tempBGmodel.InitializePWMofInstanceSet();
	//tempBGmodel.initialise(0.25);//zzz
	//tempmodel.print();
	//tempBGmodel.print();

	//if(tempBGmodel.InstanceSet.size()<bgfold*4*minCDSupport*MAXSEQNUM)//  get_consensus().size()>Setting->seedlength)
	if(LargeDataFlag)
	FOR(i,X())
	{
		double sumbg=0;
		//FOR(j,4)
		//{
		//	tempmodel.d(SearchEngine->CDProb[j],i,j);
		//	sumbg+=tempmodel.g(i,j);
		//}
		//FOR(j,4)
		//	tempmodel.d(sumbg,i,j);

		FOR(j,4)
		{
			double maxprob=max(SearchEngine->BGProb[j],SearchEngine->CDProb[j]);
			if(tempBGmodel.g(i,j)<0.5&&tempBGmodel.g(i,j)>maxprob )//*1.2
			{
			/*	if(tempmodel.g(i,j)<tempBGmodel.g(i,j))
					tempBGmodel.s(tempmodel.g(i,j),i,j);*///
				tempBGmodel.s(maxprob,i,j);
			}
			else 
				if(tempmodel.g(i,j)<PWMThreshold)
				tempmodel.s(tempBGmodel.g(i,j),i,j);
			else if(tempBGmodel.g(i,j)<PWMThreshold)
			{
				tempBGmodel.s(tempmodel.g(i,j),i,j);
			}
			
			
			
		}

	}
	else
	{
		FOR(i,X())
		{
			FOR(j,4)
				tempBGmodel.s(SearchEngine->BGProb[j],i,j);
		}
	}


	//if(Consensus.find_first_of("GACTCA")!=string::npos)
	//{
	//	tempmodel.print();
	//tempBGmodel.print();
	//}
	vector<int> effectIndex;
	double maxGain=MINSCORE;
	int maxIndex=0;
	int maxIndexJ=0;

	FOR(i,X())
	{
		if(entropy2(i)<ENTROPY_Threshold)
		{
			//continue;
			effectIndex.push_back(i);
			double sum1=0;
			double sum2=0;
			double sum3=0;
			double sum4=0;
			double sumXPc,sumXPbg;
			if(lastupdateCol==-1)
				continue;
			if(i==lastupdateCol)
			continue;
			bool flagjump=false;

			FOR(j,4)
			{
				if(tempmodel.g(i,j)>0.9)
					flagjump=true;
				sum1+=tempmodel.g(i,j);
				sum2+=tempmodel.g(i,j)/g(i,j);
			}
			if(flagjump)
				continue;
			sumXPc=sum1/sum2;
			FOR(j,4)
			{
				sum3+=tempBGmodel.g(i,j);
				sum4+=tempBGmodel.g(i,j)/g(i,j);
			}
			sumXPbg=sum3/sum4;
			FOR(j,4)
			{
				if(badmove.find(i*100+j)!=badmove.end())
					continue;
				double temp=tempmodel.g(i,j)/(g(i,j)*sum2);
				temp=temp/(tempBGmodel.g(i,j)/(g(i,j)*sum4));
				temp=temp*sumXPbg/sumXPc;
				temp*=SearchEngine->BGProb[j]/SearchEngine->CDProb[j];
				//if(tempmodel.InstanceSet.size()*(1/(g(i,j)*sum2))>minCDSupport*MAXSEQNUM)
					if(maxGain<temp&&g(i,j)<(1-PWMThreshold))
					{
						maxGain=temp;
						maxIndex=i;
						maxIndexJ=j;
					}
			}

		}
		else //if(tempmodel.entropy2(i)<ENTROPY_Threshold+0.08)
		{
			
			if(i==lastupdateCol)
				continue;
			double infogain=0;//tempBGmodel.entropy2(i)-tempmodel.entropy2(i);
			FOR(j,4)
			{

			double Pvalue=binominalTail(SearchEngine->CDProb[j],tempmodel.InstanceSet.size()*tempmodel.g(i,j),tempmodel.InstanceSet.size());
		
				if(Pvalue>1.e-3)
					continue;

				if(badmove.find(i*100+j)!=badmove.end())
					continue;
				double temp=tempmodel.g(i,j)/tempBGmodel.g(i,j);
				temp*=SearchEngine->BGProb[j]/SearchEngine->CDProb[j];
				//if(i==19)
				//	cout<<tempmodel.g(i,j)<<","<<tempBGmodel.g(i,j)<<endl;;
				infogain=temp;
				//if(tempmodel.InstanceSet.size()*tempmodel.g(i,j)>minCDSupport*MAXSEQNUM)
					if(maxGain<infogain)//&&tempmodel.g(i,j)>0.5
					{
						maxGain=infogain;
						maxIndex=i;
						maxIndexJ=j;
					}
			}
			
		}
	}

	double cds,ints,ors,cns,divs;
	ors=ORScore;
	divs=effectIndex.size();
	//tempBGmodel.print();
	//cout<<endl;
	//tempmodel.print();
	double gainlimit=4;

	
	if(abs(maxIndex)>X()||abs(maxIndexJ)>4)
		return MINSCORE;
					
				
	if(maxGain>1)
	{
		double sum=0;


		if(entropy2(maxIndex)>=ENTROPY_Threshold) //new column
		{
				effectIndex.push_back(maxIndex);

			sort(effectIndex.begin(),effectIndex.end());
			k=0;
			
			FOR(j,4)
			{
				double temp=tempmodel.g(maxIndex,j)/tempBGmodel.g(maxIndex,j);
				//temp*=SearchEngine->BGProb[j]/SearchEngine->CDProb[j];
				map<int,double>::iterator iiter=badmove.find(maxIndex*100+j);
				//if(temp<1)
				//{
				//	
				//}else
				if(temp<1)//temp<1|| ||iiter!=badmove.end()
				{
					
					if(iiter==badmove.end())
					{
						double smallvalue=PWMThreshold;
					tempmodel.s(smallvalue,maxIndex,j);
						sum+=smallvalue;
					}
					else
						//sum+=tempmodel.g(maxIndex,j);
						sum+=1;
				}
				else
				{
					if(temp>gainlimit)
						temp=gainlimit;
					//sum+=tempmodel.g(maxIndex,j)*temp;
					//tempmodel.m(temp,maxIndex,j);
					sum+=temp-1;
					tempmodel.s(temp-1,maxIndex,j);
					
				}
			}
		}
		else //old column
		{
				double sum1=0;
				double sum2=0;
				double sum3=0;
				double sum4=0;
				double sumXPc,sumXPbg;
				reviseFlag=true;
				FOR(j,4)
				{
					sum1+=tempmodel.g(maxIndex,j);
					sum2+=tempmodel.g(maxIndex,j)/g(maxIndex,j);
				}
				sumXPc=sum1/sum2;
				FOR(j,4)
				{
					sum3+=tempBGmodel.g(maxIndex,j);
					sum4+=tempBGmodel.g(maxIndex,j)/g(maxIndex,j);
				}
				sumXPbg=sum3/sum4;
				FOR(j,4)
				{
					double temp=tempmodel.g(maxIndex,j)/(g(maxIndex,j)*sum2);
					temp=temp/(tempBGmodel.g(maxIndex,j)/(g(maxIndex,j)*sum4));
					temp=temp*sumXPbg/sumXPc;
					temp*=SearchEngine->BGProb[j]/SearchEngine->CDProb[j];
					map<int,double>::iterator iiter=badmove.find(maxIndex*100+j);
					
					if(temp<1)//  ||iiter!=badmove.end()
					{
						double smallvalue=PWMThreshold;
						if(iiter==badmove.end())
						{
							tempmodel.s(smallvalue,maxIndex,j);
							sum+=smallvalue;
						}
						else
						{
								//sum+=g(maxIndex,j);
								//tempmodel.s(g(maxIndex,j),maxIndex,j);

							sum+=tempmodel.g(maxIndex,j);
							
						}
					}
					else
					{
						if(temp>gainlimit)
						temp=gainlimit;
						//sum+=g(maxIndex,j)*temp;
						//tempmodel.s(g(maxIndex,j)*temp,maxIndex,j);

						//sum+=tempmodel.g(maxIndex,j)*temp;
						//tempmodel.m(temp,maxIndex,j);
								sum+=temp-1;
							tempmodel.s(temp-1,maxIndex,j);  //the same rule as extend
							
					}
				}
			
		}
		//normalize
		FOR(j,4)
		{
			tempmodel.m(1/sum,maxIndex,j);
		}
	}
	else
	{
		return MINSCORE;
	}
	
	k=0;
	FOR(i,X())
	{
		if(k<effectIndex.size()&&i==effectIndex[k])
		{
			k++;
			if(i!=maxIndex)
			{
				FOR(j,4)
				tempmodel.s(g(i,j),i,j);
			}
			continue;
		}
		FOR(j,4)
		{
			tempmodel.s(0.25,i,j);
		}
	}
	cns=effectIndex.size();
	//tempmodel.print();
	
	tempmodel.ComputeScore(0.8,cds,ors,ints,cns,divs);

	if(tempmodel.GetMixedScore()>GetMixedScore())
	{

		//cout<<maxIndex<<","<<Pvalue<<endl;
		if(reviseFlag)
			cout<<maxGain<<","<<maxIndex<<","<<maxIndexJ<<endl;
		FOR(i,effectIndex.size())
		{
			k=effectIndex[i];
			FOR(j,4)
			{
				s(tempmodel.g(k,j),k,j);
			}

		}
		
		POSLIST=tempmodel.POSLIST;
		badmove.clear();
		//if(!switchFlag&&ors>2*ORScore)
		//{
		//	switchFlag=true;
		//	//if(LargeDataFlag)
		//	//	this->divergeSeedPart(0.01);

		//	this->divergeSeedPart(0.01);
		//}
		newScore=ORScore=ors;
		CDScore=cds;
		BindingRegion=ints;
		BinPvalue=tempmodel.BinPvalue;
	
		
	}
	else
	{
		badmove[maxIndex*100+maxIndexJ]=ors;
		lastupdateCol=maxIndex;
		return -1;
	}


		if(DEBUG)
			cout<<tempmodel.get_consensus(0)<<endl;

	return newScore-oldScore;
}

double MotifModel::GetMixedScore(double CDScore,double ORScore, double BindingRegion, double CNSVScore, double DiffScore)
{
		return ORScore;

}

double MotifModel::entropy(int col)
{int k;
		double entropy=0;
		vector<double> sorted;
		FOR(k,4)
			sorted.push_back(g(col,k));
		//entropy-=g(col,k)*log(g(col,k))/log(2.0);

		sort(sorted.begin(),sorted.end());
		double ndis=0;
		FOR(k,4)
			ndis+=(sorted[k]-0.25)*(sorted[k]-0.25);

		double parray[]={0,0.333,0.333,0.333};
		double dis=0;
		FOR(k,4)
			dis+=(sorted[k]-parray[k])*(sorted[k]-parray[k]);
		if(dis>ndis)
		{
			dis=0;
			double parray2[]={0,0,0.5,0.5};
			FOR(k,4)
			dis+=(sorted[k]-parray2[k])*(sorted[k]-parray2[k]);

			if(dis>ndis)
			{
				dis=0;
				double parray3[]={0,0,0,1};
				FOR(k,4)
				dis+=(sorted[k]-parray3[k])*(sorted[k]-parray3[k]);
				if(dis>ndis)
					return 2;
			}
		}
	return entropy;
}



	int MotifModel::AlignmentPWM(MotifModel* refmotif,MotifModel* newmotif,double& bestscore,int &bestoverlap)
	{
				           bestscore=10;
		            int bestaln=0;
                bestoverlap = 0;
				int size=newmotif->Length();
		            int size2=refmotif->Length();
               int      alnScore=10;
			   
		            int minsize=min(refmotif->Length(),newmotif->Length());
					int minoverlap=5;
					//move newmotif from left to rigth referred to refmotif
		            for(int i=0-size+minoverlap;i<size2-minoverlap+1;i++)
		            {
			            float  curScore=0;
                        int pp = 0;
                        if (i > 0)
                            pp = i;
                        int j;
			            for(j=pp;j-pp<minsize;j++)
			            {
				            if(j>=size2||(j-i)>=size)
					            break;
				            if(j>=0)
				            {
                                double sump = 0;
                                for (int p = 0; p < 4; p++)
                                {
									double temp=refmotif->g(j+refmotif->head,p)-newmotif->g(j-i+newmotif->head,p);
                                    temp*=(temp);
                                    sump += (float)temp;
                                }
                                curScore +=(float) sqrt( sump);
				            }
			            }
                        //size2=minsize
			            int overlap=j;
                        ////i>0
                        if (i > 0)
                            overlap = j - i;
						
                        curScore /= overlap * (float)sqrt(2.0);
						if(bestscore>curScore&&(overlap>=6||(overlap>=minoverlap )))////7
			            {
				            bestscore=curScore;
				            bestaln=i;
                            bestoverlap = overlap;
			            }
		            }
		            alnScore=bestscore;
		            return bestaln;
	
	}
	 int MotifModel::AlignmentPWMRC(MotifModel* refmotif,MotifModel* newmotif,double& bestScore)
	{
		MotifModel* p1=refmotif;
		refmotif->get_consensus(0);
		newmotif->get_consensus(0);
		int overlap=0;
		MotifModel* p2=newmotif;
		 double divg = 10;
            double temp;
          int bestaln=0;
            int aln=AlignmentPWM(p1, p2, temp,overlap);
            divg = temp;
			bestaln=aln;
			MotifModel p1RC=p1->RC();
            aln = AlignmentPWM(&p1RC, p2, temp,overlap);
            if (temp < divg)
			{
                divg = temp;
				bestaln=aln;
			}
          
           bestScore= divg;
			return bestaln;
	
	}
	MotifModel MotifModel::RC()
	{
		int width=X();
		MotifModel rc(SearchEngine,width,Setting);
		rc.head=tail;
		rc.ORScore=this->ORScore;
		rc.tail=head;
            for (int i = 0; i < width; i++)
            {
                rc.s( g(i, 3),width - i - 1, 0) ;
                rc.s( g(i, 2),width - i - 1, 1)  ;
                rc.s(g(i, 1),width - i - 1, 2)   ;
                rc.s(g(i, 0),width - i - 1, 3)  ;
            }

			rc.get_consensus(0);
            return rc;
	
	}


	double MotifModel::FittingCurve(double &snr,double& prbE,double& prbN)
	{
		int i,j,k;
	double* CDHist=new double[BINNUMBER];
	FOR(i,BINNUMBER)
		CDHist[i]=0;	 

	

	int BINSIZE2=maxActRLen/(BINNUMBER)+1;
  double seqweight=1;
  int startpos=Consensus.find_first_not_of("N");
	 FOR(i,POSLIST.size())
	{
		k=abs(POSLIST[i]);

		{
			int seq=k/SEQLEN;
			int pos=k%SEQLEN-SEQLEN/2+startpos;
			
			//if(abs(pos)<(bgfold*BindingRegion/2))

			if(SearchEngine->EnableWeight)
					seqweight=SearchEngine->getSeqWeight(seq);
			int tt=abs(pos);
			if(tt<maxActRLen)
				CDHist[tt/BINSIZE2]+=seqweight;

			
		
		}
	}
	
	 double sum=0;
			FOR(i,BINNUMBER)
				sum+=CDHist[i];
			FOR(i,BINNUMBER)
				CDHist[i]=CDHist[i]/sum;
			
			double avgError=fitParameter(CDHist,BINNUMBER,snr,prbE,prbN);

			/*	string name=Consensus+".his";
				name="./MatchPos/"+name;
				 FILE* outstream=fopen(name.c_str(), "w+");
			
				 if(outstream==NULL)
				 {
					 fclose(outstream);
					
				 }
				 else
				 {
					 FOR(i,BINNUMBER)
					 fprintf(outstream,"%f\n",CDHist[i]);

					  fclose(outstream);
				 
				 }*/
			delete[] CDHist;
			return avgError;
			

	}