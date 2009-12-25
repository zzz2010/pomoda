#include "stdafx.h"
#include "HashEngine.h"
#include "MotifModel.h"
#include "MemManager.h"

#define FOR(i,n) for(i=0; i<n; i++)
#define bgfold 5

string ACGT="ACGT";


HashEngine::HashEngine(void):ISearchEngine()
, TotalLength(0)
{

}

HashEngine::HashEngine(const char* filename,int hashlen,string wfile):ISearchEngine()
, TotalLength(0),Hashlen(hashlen)
{
	int addrange=1<<(Hashlen*2);

	HashIndex=new vector<VAL>* [addrange];
	
	int i;
	FOR(i,addrange)
	{
		HashIndex[i]=new vector<VAL>();
		HashIndex[i]->reserve(1000);

	}
	 EnableWeight=true;
	 weightfile=wfile;
	TotalLength=buildIndex(filename);
	

}


HashEngine::~HashEngine(void)
{
	int i;
	int addrange=1<<(Hashlen*2);
	FOR(i,addrange)
		delete HashIndex[i];


	delete [] HashIndex;
	delete []CharText;
	//delete []BGProb;
	//delete []CDProb;
}

 string HashEngine::getSite(VAL pos, int len)
 {
	 int i;
	 string ret;
	 ret.reserve(len);
	FOR(i,len)
		ret.push_back(*(CharText+pos+i));

	return ret;
	
 }


// use hash1 hash2 hash3 form the base hash address and then use the range to specify the hash adress range
void HashEngine::JoinMerge(vector<VAL>& position,int Range, int hash1, int hash2, int hash3)
{
	int i,j,k,q;
	vector<VAL> tempMerge;
	if(hash3!=-1)
	{
		FOR(i,Range)
		{
			j=k=q=0;
			int size1=(*HashIndex[hash1]).size();
			int size2=(*HashIndex[hash2]).size();
			int size3=(*HashIndex[hash3+i]).size();
			VAL pos1,pos2,pos3;
			while(j<size1&&k<size2&&q<size3)
			{
				pos1=(*HashIndex[hash1])[j];
				pos2=(*HashIndex[hash2])[k];
				pos3=(*HashIndex[hash3+i])[q];
				if(pos1==(pos2-Hashlen)&&pos2==(pos3-Hashlen))
				{
					position.push_back(pos1);
					j++;
					k++;
					q++;
				}
				if(pos1<(pos2-Hashlen))
					j++;
				if(pos2<(pos3-Hashlen))
					k++;
				if(pos3<(pos2+Hashlen))
					q++;
			}
		}
	}
	else if(hash2!=-1)
	{
		FOR(i,Range)
		{
			j=k=q=0;
			int size1=(*HashIndex[hash1]).size();
			int size2=(*HashIndex[hash2+i]).size();
			VAL pos1,pos2;
			while(j<size1&&k<size2)
			{
				pos1=(*HashIndex[hash1])[j];
				pos2=(*HashIndex[hash2+i])[k];
			
				if(pos1==(pos2-Hashlen))
				{
					position.push_back(pos1);
					j++;
					k++;
				
				}
				if(pos1<(pos2-Hashlen))
					j++;
				if(pos2<(pos1+Hashlen))
					k++;
			}
		}
	}
	else
	{
		//if(hash1>65536||hash1<0)
		//	return;
		FOR(i,Range)
		{
			j=k=q=0;
			int size1=(*HashIndex[hash1+i]).size();
			
			VAL pos1,pos2;
			while(j<size1)
			{
				pos1=(*HashIndex[hash1+i])[j];	
				{
					position.push_back(pos1);
					j++;			
				}
				
			}
		}
	}
}


int HashEngine::searchPattern( char* pattern, int mismatches, vector<VAL>& position)
{
	gPattern=pattern;
	searchPatternNRC(pattern,0,mismatches,position);
	int ret=position.size();
		//reverse searching
	string temp=reverseString(pattern);
	if(temp==pattern)
		return ret;
	char buffer[32];
	strcpy(buffer,temp.c_str());
	gPattern=buffer;
	searchPatternNRC(buffer,0,mismatches,position);
	return ret;
}

int HashEngine::searchPatternNRC( char* pattern, int pos, int mismatches, vector<VAL>& position)
{
	int len=strlen(pattern);
	
	//handle mismatch
	if(mismatches>0)
	{
		string ACGT="ACGT";
		int i,j;		
		for(i=pos;i<len;i++)
		{		
			if(pattern[i]=='N')
				continue;
			FOR(j,3)
			{
				char* temp=new char[len+1];
				strcpy(temp,pattern);
				//temp[i]='N';
				temp[i]=ACGT[(acgt(temp[i])+j+1)%4];
				searchPatternNRC(temp,pos+1,mismatches-1,position);
				delete[] temp;
			}
		}
		searchPatternNRC(pattern,pos+1,0,position);
	}
	else //0 mismatch
	{
		//handle 'N'
		if(strlen(pattern)>Hashlen)
		{
			SmartShift(position,pattern);

			return position.size();
		}
		vector<string> searchSet=MotifModel::NCloseSet(pattern);
		if(searchSet.size()==1)
		{
			int left,range;
			int hash[3];
			int i;
			FOR(i,3)
				hash[i]=-1;	
			string pa=searchSet[0];
			if(pa.size()>3*Hashlen)
				pa=pa.substr(0,3*Hashlen);
			
			FOR(i,3)
			{
				
				if(pa.size()<=(i+1)*Hashlen)
					break;
				hash[i]=getHashing(pa.c_str(),i*Hashlen,Hashlen);
				
			}	
			int restlen=pa.size()%Hashlen;
			if(restlen==0)
				restlen=Hashlen;
			hash[i]=getHashing(pa.c_str(),i*Hashlen,restlen);
			left=(Hashlen-pa.size()%Hashlen)%Hashlen;

			range=1<<(2*left);
			hash[i]<<=(2*left);
			JoinMerge(position,range,hash[0],hash[1],hash[2]);
			//reverse searching
			//pa=reverseString(pa);
			//FOR(i,3)
			//{
			//	hash[i]=-1;	
			//	if(pa.size()<=(i+1)*Hashlen)
			//		break;
			//	hash[i]=getHashing(pa.c_str(),i*Hashlen,Hashlen);
			//}			
			//hash[i]=getHashing(pa.c_str(),i*Hashlen,pa.size()%Hashlen);
			//hash[i]<<=(2*left);
			//JoinMerge(position,range,hash[0],hash[1],hash[2]);
			return position.size();
		}

		// launch queries
		int i;
		FOR(i,searchSet.size())
		{
			char* temp=new char[len+1];
			strcpy(temp,searchSet[i].c_str());
			searchPatternNRC(temp,pos+1,0,position);
			delete[] temp;
		}
	}
	return position.size();

}


VAL HashEngine::getTotalLength()
 {
	return TotalLength;
 }


//VAL HashEngine::buildIndex(const char* textfile)
//{
//	VAL start=setStartTime();
//	int maxHash=1<<(2*Hashlen);
//	ifstream openfile(textfile);
//	VAL pos=-1;
//	if(openfile)
//	{	int i,j,k;
//		
//		const int seqlen=1000;
//
//		while(openfile.good())
//		{
//			char content[seqlen+1];
//			openfile.get(content,seqlen+1);
//			int len=strlen(content);
//			FOR(i,len)
//			{
//				pos++;
//				//if(i<Hashlen-1)
//				//	continue;
//				if(i>seqlen-Hashlen)
//					continue;
//
//				int hash=getHashing(content,i,Hashlen);
//				
//				if(hash>=0&&hash<(maxHash))
//				{
//					HashIndex[hash]->push_back(pos); //already sort! since reading sequentially
//					//cout<<hash<<"\t"<<HashIndex[hash].size()<<endl;
//				}
//			}
//		}
//		cout<<"Finish Buidling Index on "<<textfile<<"\t"<<getElapsedTime(start)<<" seconds"<<endl;
//		//FOR(i,65536)
//		//	cout<<HashIndex[i]->size()<<endl;
//		openfile.close();
//	}
//	 TotalLength=pos;
//	CharText=new char[TotalLength+1];
//	 FILE *fp;
//	fp=fopen(textfile,"rb");
//	fread(CharText,sizeof(char),TotalLength,fp);
//	 fclose(fp);
//	return pos;
//}

void HashEngine::LoadWeightFile(string filename)
{
	ifstream openfile(filename.c_str());
	if(openfile)
	{	int i,j,k;
		string line="";
		while(openfile.good())
		{
			openfile>>line;
			int pos=line.find_last_of(':');
			if(pos!=string::npos)
			{
				string weigstr=line.substr(pos+1,line.size()-pos-1);
				SeqWeighting.push_back(atof(weigstr.c_str()));
			}
			//else
			//{
			//	SeqWeighting.push_back(1.0);
			//}
		}
		if(SeqWeighting.size()>0)
		cout<<"Load Weight File Successfully!"<<endl;

	}
	else
		EnableWeight=false;
		
}

double HashEngine::getSeqWeight(int seqid)
{
	if(!EnableWeight||SeqWeighting.size()==0)
		return 1;
	if(seqid>=SeqWeighting.size()||seqid<0)
		return SeqWeighting[SeqWeighting.size()-1];
	else
		return SeqWeighting[seqid];
}

VAL HashEngine::buildIndex(const char* textfile)
{
	VAL start=setStartTime();
	int maxHash=1<<(2*Hashlen);
	ifstream openfile(textfile);
	VAL pos=-1;
	if(openfile)
	{	int i,j,k;
		FOR(i,4)
		{
			BGProb[i]=0;
			CDProb[i]=0;
		}
		const int seqlen=10000; //assume max seqlen 10000
		int seqlen2=100;
		SeqLen=0;
		while(!openfile.eof())
		{
			
			char content[seqlen+1];
			string ss;
			//getline(openfile, ss);
			//char * cont = ss.c_str();
			//strcpy(content,ss.c_str());
			openfile.getline(content,seqlen+1);
			int len=strlen(content);
			if(len==0){
				continue;
			}
			if(content[0]=='>')
			{
				if(SeqLen!=0&&SeqLen>seqlen2)
					seqlen2=SeqLen;
				SeqLen=0;
				continue;
			}
			
			SeqLen+=len;

			FOR(i,len)
			{
				pos++;
				//if(i<Hashlen-1)
				//	continue;
				if(i>len-Hashlen)
					continue;

				int hash=getHashing(content,i,Hashlen);
			
				if(hash>=0&&hash<(maxHash))
				{
					HashIndex[hash]->push_back(pos); //already sort! since reading sequentially
					//cout<<hash<<"\t"<<HashIndex[hash].size()<<endl;
					int acgt=hash%4;
					int bias=abs(SeqLen-len+i-seqlen2/2);
					if(bias<150)
						CDProb[acgt]+=1;
					else if(bias>seqlen2/bgfold)
						BGProb[acgt]+=1;
					else if(seqlen2<2000)
						BGProb[acgt]+=1;
				}
			}

			
		}

		double sumpro=0;
		double sumpro2=0;
		FOR(i,4)
		{
			sumpro+=BGProb[i];
			sumpro2+=CDProb[i];
		}

		FOR(i,4)
		{
			BGProb[i]/=sumpro;
			CDProb[i]/=sumpro2;
		}
		cout<<"Finish Buidling Index on "<<textfile<<"\t"<<getElapsedTime(start)<<" seconds"<<endl;
		
		//FOR(i,65536)
		//	cout<<HashIndex[i]->size()<<endl;
		TotalLength=pos+1;
		SeqLen=seqlen2;
		cout<<"SEQLEN:"<<SeqLen<<"\t"<<"MAXSEQNUM:"<<TotalLength/SeqLen<<endl;
		CharText=new char[TotalLength+1];
		int pos2=-1;
		openfile.close();
		ifstream openfile2(textfile);

		while(openfile2.good())
		{
			char content[seqlen+1];
			openfile2.getline(content,seqlen+1);
			if(content[0]=='>')
				continue;
			int len=strlen(content);
			strcpy(CharText+pos2+1,content);
			pos2+=len;
		}

		openfile2.close();
			if(pos!=pos2)
		cout<<"Hash index building error: POS!=POS2"<<endl;
	}
	
	//int ss=weightfile.find_last_of("/");
	//int ss2=weightfile.find_last_of("\\");
	//if(ss2>ss)
	//	ss=ss2;
	//weightfile
	
	if( EnableWeight)
	LoadWeightFile(weightfile);
	return pos;
}


// use the most conservative part as a anchor, go back the hashtext to check the satisfication
void HashEngine::SmartShift(vector<VAL>& position, char* pattern)
{
	int len=strlen(pattern);
	int i,j,k,q;
	int maxaln=0;
	int bestScore=Hashlen;
	int curScore=0;
	vector<VAL> anchors;
	//find out the conservativest part
	FOR(i,Hashlen)
	{
		if((pattern[i])=='N'||(pattern[i])=='n')
		{
			bestScore--;
		}
	}
	FOR(i,len-Hashlen)
	{
		if(i==0)
			curScore=bestScore;
		else
		{
			if((pattern[i])=='N'||(pattern[i])=='n')
			{
				curScore++; //one go out
			}
			if((pattern[i+Hashlen-1])=='N'||(pattern[i+Hashlen-1])=='n')
			{
				curScore--; //one go out
			}
			if(curScore>bestScore)
			{
				maxaln=i+1;
				bestScore=curScore;
			}

		}
	}
	
	char* temp=new char[Hashlen+1];
	strncpy(temp,&pattern[maxaln],min(len-maxaln,Hashlen));
	temp[Hashlen]='\0';
	searchPatternNRC(temp,0,0,anchors);
	delete [] temp;
	vector<int> mask;
	
	FOR(i,len)
	{

		if(i>=maxaln&&i<(maxaln+Hashlen))
			continue;
		if((pattern[i])=='N'||(pattern[i])=='n')
			continue;
		
		mask.push_back(i);
	}

	int size=anchors.size();
	FOR(i,size)
	{
		int pos=anchors[i];
		bool flag=true;
		FOR(j,mask.size())
		{
			int pos2=pos-(maxaln-mask[j]);
			if(acgt(CharText[pos2])!=acgt(pattern[mask[j]]))
			{
				flag=false;
				break;
			}
		}

		if(flag)
			position.push_back(pos-maxaln);
	}
	

}


 HashEngine* HashEngine::Clone()
 {
	HashEngine* newEng=new HashEngine();
	newEng->HashIndex=this->HashIndex;
	newEng->CharText=new char[this->TotalLength];
	strcpy(newEng->CharText,this->CharText);
	newEng->Hashlen=this->Hashlen;
	newEng->gPattern=this->gPattern;
	newEng->SeqLen=this->SeqLen;
	newEng->TotalLength=this->TotalLength;
	newEng->SeqWeighting=SeqWeighting;
	newEng->EnableWeight=this->EnableWeight;
	int i;
	FOR(i,4)
	{
	newEng->CDProb[i]=this->CDProb[i];
	newEng->BGProb[i]=this->BGProb[i];
	}
	return newEng;
 
 }