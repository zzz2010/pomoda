#include "stdafx.h"
#include "EM.h"
#include <iostream>
#include <map>


EM::EM(PARAM* setting)
{
	Setting=setting;
}

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
		ErasingFactor=new double*[DATASET.size()];
		FOR(i,DATASET.size())
		{
			ErasingFactor[i]=(double *)calloc(DATASET[i].size(),sizeof(double));
		}
		
	}


}

	vector<long> EM::EMIterate(MotifModel* seed)
	{
		const int maxSeqlen=10000;
		double lamda=0.5;//(seed->ORScore-1)/seed->ORScore;
		BackModel->head=0;
		
		vector<long> sitepos;
		int i,j,k,iter;
		int motiflen=seed->Length();
	
		BackModel->tail= BackModel->X()-motiflen;
		seed->POSLIST.clear();
		seed->InstanceSet.clear();
		seed->BGInstanceSet.clear();
		seed->instWeight.clear();
		double ** bgPROBmap;
		//build a bg prob cache
		bgPROBmap=new double*[DATASET.size()];
		FOR(i,DATASET.size())
		{
			bgPROBmap[i]=(double *)calloc(DATASET[i].size(),sizeof(double));
		}
		vector<double> ErasingList;

		double sumWeight;
			double Evalue=0;
		double lastEvalue=MINSCORE;

		//iteration EM
		FOR(iter,100)
		{
			Evalue=0;
			seed->POSLIST.clear();
			seed->InstanceSet.clear();
			seed->BGInstanceSet.clear();
			seed->instWeight.clear();
			ErasingList.clear();
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
					double erasing=1;
					FOR(k,motiflen)
						erasing*=(1-ErasingFactor[i][j+k]);
					double posT=1.0/(sequence.size()-motiflen+1);
					//Epv1*=erasing;
				double	Zj=lamda*posT*Epv1/(lamda*Epv1*posT+(1-lamda*posT)*bgpv) ;//(lamda*Epv1+(1-lamda)*bgpv);
				
				//if(ErasingFactor.find(j+maxSeqlen*i)!=ErasingFactor.end())
				
				adhocSum+=Zj;
				
				//if(Zj<0.5)
				//	continue;
				if(Zj>=0.5)
				ErasingList.push_back(j+maxSeqlen*i+Zj);
				
				Zj=Zj*erasing; //time the erasing factor after putting into the list

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
				//maxscore=maxscore*Prior[i];
	
				
				seed->InstanceSet.push_back(matchsite);
				seed->instWeight.push_back(maxscore);
				sumWeight+=maxscore;
			}
			//update lamda
			lamda=(double)sumWeight/DATASET.size(); //
			if(lamda<Setting->min_supp_ratio)
				lamda=2*Setting->min_supp_ratio; //double it to search

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
	
			//seed->print();
				double sumlog=0;
				FOR(i,motiflen)
					FOR(j,4)
					{
						double v=seed->g(i+seed->head,j);
						sumlog+=v*log(v);
					}
			
				if(lamda!=0)
					Evalue=lamda*(log(lamda)+sumlog);
				double bgsumlog=0;
				FOR(j,4)
				{
					double v=BackModel->g(0,j);
					bgsumlog+=v*log(v);
				}
				if(lamda!=1)
					Evalue+=(1-lamda)*(log(1-lamda)+motiflen*bgsumlog);
				cout<<seed->get_consensus(-1)<<"\t"<<sumWeight<<"\t"<<Evalue<<endl;
				if(algnScore<0.01||lamda<=0||lamda>=1)
					break;
				//restore the good one
				//if(lastEvalue>Evalue)
				//{
				//		FOR(i,seed->X())
				//			FOR(j,4)
				//				seed->s(tempmodel.g(i,j),i,j);

				//		break;
				//}
				lastEvalue=Evalue;
		}
		
		seed->ORScore=lamda/(1-lamda);
		seed->CDScore=sumWeight;

		seed->SeqPvalue=exp(Evalue);

		
		//only consider Zvalue>0.5
		FOR(i,seed->POSLIST.size())
		{
			if(seed->instWeight[i]<0.5)// suppose the prior is 0.5
			{
				seed->CDScore-=seed->instWeight[i];
				seed->instWeight[i]=0;
				seed->POSLIST[i]=-1;

			}
		}

		cout<<seed->get_consensus(0)<<"\t"<<seed->ORScore<<"\t"<<seed->CDScore<<"\t"<<seed->SeqPvalue<<endl;
		//update erasing factor and prior
		
		FOR(i,ErasingList.size())
		{
		
			int pos=floor(ErasingList[i]);
			int seq=pos/maxSeqlen;
			double Px=ErasingList[i]-pos;//get back the conditional probablity
			pos=pos%maxSeqlen;
			//Px=1-pow(1-Px,1.0/motiflen);
			FOR(j,motiflen)
			{
					ErasingFactor[seq][pos+j]=ErasingFactor[seq][pos+j]+(1-ErasingFactor[seq][pos+j])*Px;
			}
			//update Prior =P(X|theta)+(1-P(X|theta))*Prior			
			//Prior[seq]=Px+(1-Px)*Prior[seq];
		}
		//noramlize Prior: prior is not good just using main motifs
		//double sumPrior=0;
		//FOR(i,Prior.size())
		//{
		//	sumPrior+=Prior[i];
		//}
		//sumPrior=sumPrior/(0.5*DATASET.size());
		//FOR(i,Prior.size())
		//	Prior[i]=Prior[i]/sumPrior;
		
		//free memory
		FOR(i,DATASET.size())
		{
			free(bgPROBmap[i]);
		}
		delete[] bgPROBmap;

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
		//if(candidates[i]->Length()<10)
		{
			candidates[i]->head=max(candidates[i]->head-2,0);
			candidates[i]->tail=max(candidates[i]->tail-2,0);
		}
		//EM refine the model
		vector < long > temp=EMIterate(candidates[i]);
		poslist.push_back(temp);
		occCounter[i]=candidates[i]->CDScore;
		candidates[i]->deltaScore=0; //record the overlap ratio
	}

	vector<long> markedPos;
	//select max mutual cover set
	FOR(i,min(outMotifNum,candidates.size()))
	{
		int maxCount=0;
		int index=0;
		FOR(j,candidates.size())
		{
			//if(candidates[j]->Length()<7)
			//	occCounter[j]=0;

			if(occCounter[j]>maxCount)//  &&candidates[j]->POSLIST.size()>minSupport  &&candidates[j]->deltaScore<(Setting->olThresh*candidates[j]->CDScore)
			{
				index=j;
				maxCount=occCounter[j];
			}
		}

		// no more good motif
		if(maxCount==0)
			break;
		//mark this motif in later loop
		occCounter[index]=0;
		// add this motif position to marked Pos
		FOR(j,candidates[index]->POSLIST.size())
		{
			if(candidates[index]->POSLIST[j]!=-1)
			{
				if(candidates[index]->instWeight[j]>0.05) //dumy check
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
							// increase the overlap counter
							candidates[j]->deltaScore+=candidates[j]->instWeight[y];
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
							// increase the overlap counter
							candidates[j]->deltaScore+=candidates[j]->instWeight[y];
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