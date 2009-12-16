#include "EM.h"
#include "stdafx.h"
#include <iostream>


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
			totallen+=ss.size();
			FOR(k,ss.size())
			{
				int id=acgt(ss[k]);
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
		BackModel->print();
	}

}

	vector<long> EM::EMIterate(MotifModel* seed)
	{
		const int maxSeqlen=10000;
		double lamda=0.5;
		BackModel->head=0;
		
		vector<long> sitepos;
		int i,j,k,iter;
		int motiflen=seed->Length();

		BackModel->tail= BackModel->X()-motiflen;
		seed->POSLIST.clear();
		seed->InstanceSet.clear();
		seed->BGInstanceSet.clear();
		vector<double> weightlist;
		//iteration EM
		FOR(iter,1000)
		{
			seed->POSLIST.clear();
			seed->InstanceSet.clear();
			seed->BGInstanceSet.clear();
			weightlist.clear();
		//assume each sequence contain one occurrence
			double sumWeight=0;
			FOR(i,DATASET.size())
			{
				double maxscore=MINSCORE;
				int matchpos=0;
				bool rc=false;
				string sequence=DATASET[i];
				FOR(j,sequence.size()-motiflen+1)
				{
					string site=sequence.substr(j,motiflen);
					double pv1=seed->getLogProb(site.c_str());
					double pv2=seed->getLogProb(reverseString(site).c_str());
					if(pv1<pv2)
					{
						pv1=pv2;
					}
					if(pv1==MINSCORE)
						continue;
					double bgpv=exp(BackModel->getLogProb(site.c_str()));
					double Epv1=exp(pv1);
					Epv1=lamda*Epv1/(lamda*pv1+(1-lamda)*bgpv);
					if(Epv1>maxscore)
					{
						maxscore=Epv1;
						matchpos=j;
						if(pv1==pv2)
							rc=true;
					}
				}
				seed->POSLIST.push_back(matchpos+maxSeqlen*i); 
				string matchsite=sequence.substr(matchpos,motiflen);
				if(rc)
					matchsite=reverseString(matchsite);
				seed->InstanceSet.push_back(matchsite);
				weightlist.push_back(maxscore);
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
			seed->InitializePWMofInstanceSet(weightlist);
			double algnScore=1000;
			int overlap=0;
			seed->AlignmentPWM(seed,&tempmodel,algnScore,overlap);
			seed->print();
			if(algnScore<0.001)
				break;
		}


		return sitepos;
	}


vector<MotifModel*> EM::LoadSeedModels(vector<MotifModel*> candidates, int outMotifNum)
{
	vector<MotifModel*> retList;
	int i,j,k;
	const int maxSeqlen=10000;
	vector< vector < long> > poslist;
	FOR(i,candidates.size())
	{
		//EM refine the model
		vector < long > temp=EMIterate(candidates[i]);
		poslist.push_back(temp);
	}

	vector<long> markedPos;
	//select max mutual cover set
	FOR(i,outMotifNum)
	{
		int maxCount=0;
		int index=0;
		FOR(j,candidates.size())
		{
			if(candidates[j]->CDScore>maxCount)
			{
				index=j;
			}
		}
		//mark this motif in later loop
		candidates[index]->CDScore=0;
		// add this motif position to marked Pos
		FOR(j,candidates[index]->POSLIST.size())
		{
			if(candidates[index]->POSLIST[j]!=-1)
			markedPos.push_back(candidates[index]->POSLIST[j]);
		}
		//clean up the rest position
		FOR(j,candidates.size())
		{
			if(candidates[j]->CDScore!=0)
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
							int seq=pos2/maxSeqlen;
							int pos=pos2%maxSeqlen;
							string site=DATASET[seq].substr(pos,restlen);
							double pv1=candidates[j]->getLogProb(site.c_str());
							double pv2=candidates[j]->getLogProb(reverseString(site).c_str());
							if(pv1<pv2)
								pv1=pv2;
							candidates[j]->CDScore=candidates[j]->CDScore-exp(pv1);
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
								int seq=pos2/maxSeqlen;
							int pos=pos2%maxSeqlen;
							string site=DATASET[seq].substr(pos,restlen);
							candidates[j]->POSLIST[y]=-1;
							double pv1=candidates[j]->getLogProb(site.c_str());
							double pv2=candidates[j]->getLogProb(reverseString(site).c_str());
							if(pv1<pv2)
								pv1=pv2;
							candidates[j]->CDScore=candidates[j]->CDScore-exp(pv1);
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
	return retList;
}