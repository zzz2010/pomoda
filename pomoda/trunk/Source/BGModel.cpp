#include "stdafx.h"
#include "BGModel.h"
#include "MotifModel.h"



BGModel::BGModel(string bgfile, int MarkovOrder):FlMatrix(4,4,4,4,"BackGround Model"),Markov0(4),Markov1(4,4),Markov2(4,4,4)
{
	Markov0.initialise(0);
	Markov1.initialise(0);
	Markov2.initialise(0);
	cout<<"read bg file:" <<bgfile<<endl;
	ifstream openfile(bgfile.c_str());
	this->initialise(0.0000000001);
	if(openfile)
	{
		int seqid=0;
		while(openfile.good())
		{
			string line="";
			openfile>>line;
			if(line=="")
				continue;
			if(line[0]=='>')
			{
	
				continue;
			}
			else
			{			
					string temp=line;
					for(int j=0;j<temp.size()-3;j++)
					{
						int dd[4];
						int k;
						bool flag=false;
						FOR(k,4)
						{
						dd[k]=acgt(temp[j+k]);
						if(dd[k]<0||dd[k]>3)
							flag=true;
						}
						if(flag)
							continue;
						
						a(1,dd[0],dd[1],dd[2],dd[3]);
						Markov0.a(1,dd[0]);
						Markov1.a(1,dd[0],dd[1]);
						Markov2.a(1,dd[0],dd[1],dd[2]);

						
					}
				
			}
		}


		int i,j,k,q;
		double sum0=0;
		
		FOR(i,X())
		{
			sum0+=Markov0.g(i);
			double sum1=0;
				FOR(j,Y())
				{
					sum1+=Markov1.g(i,j);
					double sum2=0;
					FOR(k,Z())
					{
						sum2+=Markov2.g(i,j,k);
						double sum3=0;
						FOR(q,W())
						{
							sum3+=g(i,j,k,q);
						}

						FOR(q,W())
						{
							
							d(sum3,i,j,k,q);
							
						}

					}
					FOR(k,Z())
						Markov2.d(sum2,i,j,k);
				}
				FOR(j,Y())
				Markov1.d(sum1,i,j);
				
		}

		FOR(i,X())
		{
			
			Markov0.d(sum0,i);
	
			
		}

		openfile.close();
	}
	else
		cout<<"BG file error!"<<endl;

	
}


BGModel::BGModel(void):FlMatrix(4,4,4,4,"BackGround Model"),Markov0(4),Markov1(4,4),Markov2(4,4,4)
{
	Markov0.initialise(0.25);
	Markov1.initialise((double)1/16);
	Markov2.initialise((double)1/64);
	this->initialise((double)1/256);
}

BGModel::BGModel(map<string, double>& instSet):FlMatrix(4,4,4,4,"BackGround Model"),Markov0(4),Markov1(4,4),Markov2(4,4,4)
{
		Markov0.initialise(0.0000000001);
	Markov1.initialise(0.0000000001);
	Markov2.initialise(0.0000000001);

	this->initialise(0.0000000001);
	map<string, double>::iterator ITER=instSet.begin();
		int seqid=0;
		while(ITER!=instSet.end())
		{		
				string temp=ITER->first;
					for(int j=0;j<temp.size()-3;j++)
					{
						int dd[4];
						int k;
						bool flag=false;
						FOR(k,4)
						{
						dd[k]=acgt(temp[j+k]);
						if(dd[k]<0||dd[k]>3)
							flag=true;
						}
						if(flag)
							continue;
						
						a(ITER->second,dd[0],dd[1],dd[2],dd[3]);
						Markov0.a(ITER->second,dd[0]);
						Markov1.a(ITER->second,dd[0],dd[1]);
						Markov2.a(ITER->second,dd[0],dd[1],dd[2]);					
					}
			ITER++;
		}


		int i,j,k,q;
		double sum0=0;
		
		FOR(i,X())
		{
			sum0+=Markov0.g(i);
			double sum1=0;
				FOR(j,Y())
				{
					sum1+=Markov1.g(i,j);
					double sum2=0;
					FOR(k,Z())
					{
						sum2+=Markov2.g(i,j,k);
						double sum3=0;
						FOR(q,W())
						{
							sum3+=g(i,j,k,q);
						}

						FOR(q,W())
						{
							
							d(sum3,i,j,k,q);
							
						}

					}
					FOR(k,Z())
						Markov2.d(sum2,i,j,k);
				}
				FOR(j,Y())
				Markov1.d(sum1,i,j);
				
		}

		FOR(i,X())
		{
			
			Markov0.d(sum0,i);
	
			
		}
		

}

BGModel::~BGModel(void)
{
	
	//Markov0.~FlMatrix();
	//Markov1.~FlMatrix();
	//Markov2.~FlMatrix();
}

//double BGModel::getStringProb(string pattern)
//{
//	int d1=acgt(pattern[0]);
//	int d2=acgt(pattern[1]);
//	int d3=acgt(pattern[2]);
//	int i,j,k,q;
//	double prob1=0;
//	double prob2=0;
//	double prob3=0;
//	
//	prob1=Markov0.g(d1);
//	prob2=Markov1.g(d1,d2);
//	prob3=Markov2.g(d1,d2,d3);
//	
//	double PROB=log(prob1)+log(prob2)+log(prob3); 
//	for(i=3;i<pattern.size();i++)
//	{
//		PROB+=log(g(acgt(pattern[i-3]),acgt(pattern[i-2]),acgt(pattern[i-1]),acgt(pattern[i])));
//	}
//	
//	double ret=exp(PROB)+0.000001;
//	return ret;
//}

double BGModel::getStringProb(string pattern)
{
	if(pattern.size()==1)
	{
		int d1=acgt(pattern[0]);
		if(NOTACGT(d1))
			return 1;
		return Markov0.g(d1);
	}
	if(pattern.size()==2)
	{
				int d1=acgt(pattern[0]);
				int d2=acgt(pattern[1]);	
				if(NOTACGT(d2))
					return getStringProb(pattern.substr(0,1));
				if(NOTACGT(d1))
					return getStringProb(pattern.substr(1,1)); //getStringProb(pattern.substr(0,pattern.size()-1))*
		return getStringProb(pattern.substr(0,pattern.size()-1))*Markov1.g(d1,d2);
	}
	if(pattern.size()==3)
	{
				int d1=acgt(pattern[0]);
				int d2=acgt(pattern[1]);
				int d3=acgt(pattern[2]);
				if(NOTACGT(d3))
					return getStringProb(pattern.substr(0,2));
				if(NOTACGT(d1))
				{
					return getStringProb(pattern.substr(1,2));
				}
				if(NOTACGT(d2))
				{
					int i;
					double ret=0;
					//FOR(i,4)
					//	ret+=Markov2.g(d1,i,d3);
					ret=getStringProb(pattern.substr(2,1));
					
					return getStringProb(pattern.substr(0,1))*ret;
				}
		return getStringProb(pattern.substr(0,pattern.size()-1))*Markov2.g(d1,d2,d3);
	}
	if(pattern.size()==4)
	{
				int d1=acgt(pattern[0]);
				int d2=acgt(pattern[1]);
				int d3=acgt(pattern[2]);
				int d4=acgt(pattern[3]);
				if(NOTACGT(d4))
					return getStringProb(pattern.substr(0,pattern.size()-1));
				if(NOTACGT(d1))
					return getStringProb(pattern.substr(1,3));
				if(NOTACGT(d2))
				{
					double ret=0;
					int i,j;
					//FOR(i,4)
					//{
					//	if(NOTACGT(d3))
					//	{
					//		FOR(j,4)
					//			ret+=g(d1,i,j,d4);

					//	}
					//	else
					//	{
					//		ret+=g(d1,i,d3,d4);
					//	}
					//}
					
					return getStringProb(pattern.substr(0,1))*getStringProb(pattern.substr(2,2));
				}
				if(NOTACGT(d3))
				{
					double ret=0;
					int j;
					//FOR(j,4)
					//	ret+=g(d1,d2,j,d4);
					return getStringProb(pattern.substr(0,2))*getStringProb(pattern.substr(3,1));
				}
		return	getStringProb(pattern.substr(0,pattern.size()-1))*g(d1,d2,d3,d4);
	}
	else
	{
				int size=pattern.size();
				int d1=acgt(pattern[size-4]);
				int d2=acgt(pattern[size-3]);
				int d3=acgt(pattern[size-2]);
				int d4=acgt(pattern[size-1]);
				if(NOTACGT(d4))
					return getStringProb(pattern.substr(0,pattern.size()-1))*1;
				if(NOTACGT(d1))
				{
					double ret=0;
					int i,j,k;
					return getStringProb(pattern.substr(0,pattern.size()-4))*getStringProb(pattern.substr(pattern.size()-3,3));
				}
				if(NOTACGT(d2))
				{
					double ret=0;
					int i,j;
					//FOR(i,4)
					//{
					//	if(NOTACGT(d3))
					//	{
					//		FOR(j,4)
					//			ret+=g(d1,i,j,d4);

					//	}
					//	else
					//	{
					//		ret+=g(d1,i,d3,d4);
					//	}
					//}
					return getStringProb(pattern.substr(0,pattern.size()-3))*getStringProb(pattern.substr(pattern.size()-2,2));
				}
				if(NOTACGT(d3))
				{
					double ret=0;
					int j;
					//FOR(j,4)
					//	ret+=g(d1,d2,j,d4);
					return getStringProb(pattern.substr(0,pattern.size()-2))*getStringProb(pattern.substr(pattern.size()-1,1));
				}
		return	getStringProb(pattern.substr(0,pattern.size()-1))*g(d1,d2,d3,d4);
		
	}
}



void BGModel::GenBGData(int taglength)
{
	int i,j,k,q;
	string content="TACG";
	string ACGT="ACGT";
	FOR(i,(taglength*1000)-4)
	{
		
				int d1=acgt(content[i]);
				int d2=acgt(content[i+1]);
				int d3=acgt(content[i+2]);
				double sum=0;
				double p=nextDouble;
				FOR(j,4)
				{
					sum+=g(d1,d2,d3,j);
					if(p<sum)
						break;
				}
				content.push_back(ACGT[j]);
	
	}
	string filename="bg"+Int2String(taglength)+".Seq";
	ofstream writefile(filename.c_str());
	writefile<<content<<endl;
	writefile.close();
	cout<<"Finish Generate the Data"<<endl;
}

	void BGModel::MixMotifModel(void* M,double factor)
	{
		MotifModel * m=(MotifModel*) M;
		map<string, double> insts=m->GenerateInstanceFromPWM(0.5);
		BGModel tempbg(insts);
		this->MixBGModel(tempbg,factor);
		
	}

	void BGModel::MixBGModel(BGModel& bg,double factor)
	{
				int i,j,k,q;

		
		FOR(i,X())
		{
			
				FOR(j,Y())
				{
					if(bg.Markov1.g(i,j)!=0.25)
					Markov1.s( Markov1.g(i,j)*(1-factor)+factor*bg.Markov1.g(i,j),i,j);


					FOR(k,Z())
					{
						if(bg.Markov2.g(i,j,k)!=0.25)
						Markov2.s(Markov2.g(i,j,k)*(1-factor)+factor*bg.Markov2.g(i,j,k),i,j,k);
	
						FOR(q,W())
						{
							if(bg.g(i,j,k,q)!=0.25)
							s(g(i,j,k,q)*(1-factor)+factor*bg.g(i,j,k,q),i,j,k,q);
						
						}

					

					}
					
				}
				
				
		}

		//FOR(i,X())
		//{
		//	
		//	if(bg.Markov0.g(i)!=0.25)
		//	Markov0.s(Markov0.g(i)*(1-factor)+factor*bg.Markov0.g(i),i);
		//}
	}

double BGModel::getPWMThreshold(FlMatrix* pwm,double FDR)
{
	int i,j,k,q;
	string content="TACG";
	string ACGT="ACGT";
	vector<double> score;
	FOR(i,(1000000))
	{
		
			int d1=acgt(content[i]);
			int d2=acgt(content[i+1]);
			int d3=acgt(content[i+2]);
			double sum=0;
			double p=nextDouble;
			FOR(j,4)
			{
				sum+=g(d1,d2,d3,j);
				if(p<sum)
					break;
			}
			content.push_back(ACGT[j]);
			if(i> pwm->X())
			{
			string temp=content.substr(i- pwm->X(), pwm->X());
			double sc=0;
			FOR(j,temp.size())
			{
				double sss=pwm->g(j,acgt(temp[j]));
				if(sss==0)
				{
					sc=MINSCORE;
					break;
				}
				sc+=log(sss);
			}
			score.push_back(sc);
			}
	}
	sort(score.begin(),score.end());
	int totlen=1000000-pwm->X();
	return score[totlen-totlen*FDR-1];
}