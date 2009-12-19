#include "stdafx.h"
#include "EM.h"
#include <iostream>
#include <map>


void EM::LoadSeqFile(string filename)
{
		ifstream openfile(filename.c_str());
		BackModel=new MotifModel();
	if(openfile)
	{	int i,j,k;
		double bgprob[4];
		bgprob[0]=bgprob[1]=bgprob[2]=bgprob[3]=0;
		const int seqlen=10000; //assume max seqlen 10000
		int seqlen2=100;
		int totallen=0;
		while(!openfile.eof())
		{
			
			char content[seqlen+1];
		
			openfile.getline(content,seqlen+1);
			int len=strlen(content);
			if(len==0){
				continue;
			}
			if(content[0]=='>')
			{
				continue;
			}
			string ss=string(content,len);
			DATASET.push_back(ss);
			Prior.push_back(0.5);
			totallen+=ss.size();
			FOR(k,ss.size())
			{
				int id=acgt(ss[k]);
				if(id<0||id>3)
				{
					totallen--;
					continue;
				}
				bgprob[id]+=1;
			}
		}
		FOR(i,4)
			bgprob[i]=bgprob[i]/totallen;
		FOR(i,BackModel->X())
		{
			FOR(j,4)
				BackModel->s(bgprob[j],i,j);
		}
	
	}

}

	vector<long> EM::EMIterate(MotifModel* seed)
	{
		const int maxSeqlen=10000;
		double lamda=0.5; //(seed->ORScore-1)/seed->ORScore;
		BackModel->head=0;
		
		vector<long> sitepos;
		int i,j,k,iter;
		int motiflen=seed->Length();

		BackModel->tail= BackModel->X()-motiflen;
		seed->POSLIST.clear();
		seed->InstanceSet.clear();
		seed->BGInstanceSet.clear();
		seed->instWeight.clear();
		vector<map<int,double> > bgPROBmap;

		double sumWeight;
		//build a bg prob cache
		FOR(i,DATASET.size())
		{
			bgPROBmap.push_back(map<int,double>());
		}
		//iteration EM
		FOR(iter,100)
		{
			seed->POSLIST.clear();
			seed->InstanceSet.clear();
			seed->BGInstanceSet.clear();
			seed->instWeight.clear();

		//assume each sequence contain at most one occurrence
			 sumWeight=0;
		// lamda=0.9;//fix lamda
			FOR(i,DATASET.size())
			{
				double maxscore=MINSCORE;
				string matchsite;
				int matchpos=0;
				bool rc=false;
				string sequence=DATASET[i];
				double normalFactor=0;
				int nSplit=1;
				double adhocSum=0;
				FOR(j,sequence.size()-motiflen+1)
				{
					
					string site=sequence.substr(j,motiflen);
					double pv1=seed->getLogProb(site.c_str());
					string rcsite=reverseString(site);
					double pv2=seed->getLogProb(rcsite.c_str());
					if(pv1<pv2)
					{
						pv1=pv2;
					}
					if(pv1==MINSCORE)
					{
						nSplit++;
						j=j+motiflen-1;
						if(adhocSum>normalFactor)
						normalFactor=adhocSum;
						adhocSum=0;
						continue;
					}
					double bgpv=0;
					if(iter==0)
					{

						bgpv=exp(BackModel->getLogProb(site.c_str()));
						bgPROBmap[i][j]=bgpv;
					}
					else
						bgpv=bgPROBmap[i][j];
					double Epv1=exp(pv1);
				double	Zj=lamda*Epv1/(lamda*(Epv1+(1-Epv1)*bgpv)+(1-lamda)*bgpv) ;//(lamda*Epv1+(1-lamda)*bgpv);
				
				if(ErasingFactor.find(j+maxSeqlen*i)!=ErasingFactor.end())
					Zj=Zj*(1-ErasingFactor[j+maxSeqlen*i]);
				adhocSum+=Zj;
				double Sx=Epv1/bgpv;
				if(Sx<(1-lamda)/lamda)
					continue;
					if(Zj>maxscore)
					{
						maxscore=Zj;
						matchpos=j;
						
						if(pv1==pv2)
						{
				
							rc=true;
							matchsite=rcsite;
							
						}
						else
						{
					
							rc=false;
							matchsite=site;
						}
					}
				}
				if(maxscore<0)
					continue;
				seed->POSLIST.push_back(matchpos+maxSeqlen*i); 
				
				//if(normalFactor>nSplit)
				//	maxscore=maxscore/(normalFactor/nSplit);
				/**********no normalization***********/
				//if (normalFactor>1)
				//	maxscore=maxscore/(normalFactor/1);
				/**********no normalization***********/
				//add prior
				maxscore=maxscore*Prior[i];
	
				
				seed->InstanceSet.push_back(matchsite);
				seed->instWeight.push_back(maxscore);
				sumWeight+=maxscore;	
			}
			//update lamda
			lamda=sumWeight/DATASET.size();

			//copy backup currenct matrix
			MotifModel tempmodel(seed->SearchEngine,seed->X(),seed->Setting);
			tempmodel.head=seed->head;
			tempmodel.tail=seed->tail;
			FOR(i,seed->X())
				FOR(j,4)
				tempmodel.s(seed->g(i,j),i,j);
			//form new PWM
			seed->InitializePWMofInstanceSet(seed->instWeight);

			double algnScore=1000;
			int overlap=0;
			seed->AlignmentPWM(seed,&tempmodel,algnScore,overlap);
			cout<<seed->get_consensus(-1)<<endl;
			//seed->print();
			if(algnScore<0.001||lamda<=0||lamda>=1)
				break;
		}
		
		seed->ORScore=lamda/(1-lamda);
		seed->CDScore=sumWeight;
		double sumlog=0;
		FOR(i,motiflen)
			FOR(j,4)
			{
				double v=seed->g(i+seed->head,j);
				sumlog+=v*log(v);
			}
		double Evalue=lamda*(log(lamda)+sumlog);
		double bgsumlog=0;
		FOR(j,4)
		{
			double v=BackModel->g(0,j);
			bgsumlog+=v*log(v);
		}
		Evalue+=(1-lamda)*(log(1-lamda)+motiflen*bgsumlog);
		seed->SeqPvalue=exp(Evalue);
		cout<<seed->get_consensus(0)<<"\t"<<seed->ORScore<<"\t"<<seed->CDScore<<"\t"<<seed->SeqPvalue<<endl;

		//update erasing factor and prior
		
		FOR(i,seed->POSLIST.size())
		{
			int pos=seed->POSLIST[i];
			int seq=pos/maxSeqlen;
			double Px=seed->instWeight[i]/Prior[seq];//get back the conditional probablity
			FOR(j,motiflen)
			{
				if(ErasingFactor.find(pos+j)!=ErasingFactor.end())
				{
					ErasingFactor[pos+j]=ErasingFactor[pos+j]+(1-ErasingFactor[pos+j])*Px;
				}
				else
					ErasingFactor[pos+j]=Px;
			}
			//update Prior =P(X|theta)+(1-P(X|theta))*Prior
			
			Prior[seq]=Px+(1-Px)*Prior[seq];
		}
		//noramlize Prior
		double sumPrior=0;
		FOR(i,Prior.size())
		{
			sumPrior+=Prior[i];
		}
		sumPrior=sumPrior/(0.5*DATASET.size());
		FOR(i,Prior.size())
			Prior[i]=Prior[i]/sumPrior;
		


		return sitepos;
	}


vector<MotifModel*> EM::LoadSeedModels(vector<MotifModel*> candidates, int outMotifNum)
{
	vector<MotifModel*> retList;
	int i,j,k;
	const int maxSeqlen=10000;
	vector< vector < long> > poslist;
	i=0;
	double * occCounter=new double[candidates.size()];
	FOR(i,candidates.size())
	{
		//extend the length of short motif
		if(candidates[i]->Length()<10)
		{
			candidates[i]->head=max(candidates[i]->head-2,0);
			candidates[i]->tail=max(candidates[i]->tail-2,0);
		}
		//EM refine the model
		vector < long > temp=EMIterate(candidates[i]);
		poslist.push_back(temp);
		occCounter[i]=candidates[i]->CDScore;
	}

	vector<long> markedPos;
	//select max mutual cover set
	FOR(i,outMotifNum)
	{
		int maxCount=0;
		int index=0;
		FOR(j,candidates.size())
		{
			if(candidates[j]->Length()<7)
				occCounter[j]=0;

			if(occCounter[j]>maxCount)
			{
				index=j;
				maxCount=occCounter[j];
			}
		}
		//mark this motif in later loop
		occCounter[index]=0;
		// add this motif position to marked Pos
		FOR(j,candidates[index]->POSLIST.size())
		{
			if(candidates[index]->POSLIST[j]!=-1)
			{
				if(candidates[index]->instWeight[j]>0.05)
				markedPos.push_back(candidates[index]->POSLIST[j]);
				else
				candidates[index]->POSLIST[j]=-1;
			}
		}
		//clean up the rest position
		FOR(j,candidates.size())
		{
			if(occCounter[j]!=0)
			{
				int x, y;
				x=y=0;
				vector <long> mark=candidates[index]->POSLIST;
				int marklen=candidates[index]->Length();
				vector <long> rest=candidates[j]->POSLIST;
				int restlen=candidates[j]->Length();
				while(x<mark.size()&&y<rest.size())
				{
					int pos1=mark[x];
					int pos2=rest[y];
					if(pos1<pos2)
					{
						if(pos1+marklen>pos2)
						{
							//do clean up

							occCounter[j]=occCounter[j]-candidates[j]->instWeight[y];
							candidates[j]->POSLIST[y]=-1;
							x++;
							y++;
						}
						else
							x++;
					}
					else
					{
						if(pos2+restlen>pos1)
						{
							//do clean up
						
							candidates[j]->POSLIST[y]=-1;
			            	occCounter[j]=occCounter[j]-candidates[j]->instWeight[y];
							x++;
							y++;
						}
						else
							y++;
					}
				}
			}
		}
		retList.push_back(candidates[index]);
	}
	delete[] occCounter;
	return retList;
}