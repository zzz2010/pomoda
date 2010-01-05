/*************************************
* Pomoda: Peak Oriented Motif Discovery Algorithm 
* http://compbio.ddns.nus.edu.sg/~chipseq/Pomoda
*
* Function: generates position weight matrixes (PWM) for motifs concentrate around the ChIP-Seq peaks
* 
* Syntax:
* POMODA -i inputFasta -o outputDIR [-w weightFile, -ot overlapThreshold -pdt PWM_Divergence_Threshold
* -ratio minSupportRatio, -FDR p-valueCutoff, -maxlen maxMotifLength, -seedlen SeedLength,-rs resolution_bp, -mbr min_binding_range
* -n numberOfMotifs]
* 
* Comments:
* inputFasta: input data, using fasta format, recommended length per sequence is 10k
* outputDIR: the directory for output of one task, different task should have different outputDIR, otherwise, the output position file will be replaced.
* weightFile: the file store the peak intensity, each line corresponding to one sequence in inputFasta, e.g., ":3251" in the 1st line, says the first sequence intensity is 3251, [optional]
* overlapThreshold: the binding site overlap percentage cutoff for filtering the redundant candidate motif. The higher value, the less motifs will be filtered. [default: 0.02]
* PWM_Divergence_Threshold: the dissimilarity cuttof for filtering the redundant candidate motif. The higher value, the more motifs will be filtered. [default:0.18]
* minSupportRatio: the minimal proportion of sequences contains the reported moitf [default:0.05]. 
* FDR: the p-value cuttoff for testing uniformilty of motif distribution. [default:1.e-7]
* maxMotifLength: the maximum length of reported motif. [default:26]
* SeedLength: the length of k-mer using in seed selection phase. [default:6]
* rs: the resolution level of binding region estimation (unit:bp).  [default:100]
* min_binding_range: the minimal size of binding region. [default:100]
* numberOfMotifs: the number of output motifs. [default:20]
*
* Author: Zhizhuo Zhang   10/10/2009
* Computational Biology Lab, School of Computing, National University of Singapore (NUS)
* Email: zzz2010@gmail.com
* 
* 
* If you use this program in your research, please cite:
* Zhizhuo et.al 2009. Pomoda: Peak Oriented Motif Discovery Algorithm 
* BMC Bioinformatics (in review).
**************************************/
#include "stdafx.h"
#include "MotifModel.h"
#include <sstream>
#include "HashEngine.h"
#include "EM.h"

//string seqFile="select_m_10k.fa";//"snber_1000_sorted.fasta";//"select_H_sapiens_UCSC_hg18_10000.fas";//"E2F_Q4V1200.fasta";//"iterSim/0U.fasta";//"testset/select_MAKEW14S1A2V200L10000N2000U.fasta";//"select_MAKEW7S1A1V1200L10000N3000U.fasta";//"P$CBT_01V200.fasta";//"testset/simdata1.fas";//"select_p300.fasta";//"select_ctcf.fasta";//"snber_1000_sorted.fasta";//"select_MAKEW8S1A1V50.fasta";//"select_m_10k.fa";//"select_i10k.fasta";//"seq_ctcf_sort_intensity.fa";//"select_m_Q1_2000.faW16S1A1V100.fa";//"select_cmyc.fasta";//"select_m_Q1_2000.fa";///"select_m_Q1_2000.faW8S1A1V400.fa";//"select_2fold_im_m.fa.SeqRev";//"select_i_Q1.fa.SeqRev";//"select_2fold_im_i.fa.SeqRev";//"select_m_Q1.fa.SeqRev";//"bg1000.SeqW12S3A4V200.Seq";//"bg1000.SeqW16S1A1V100.Seq";//"select_2fold_im_m.fa.SeqRev";//"MCF7E2.txt.SeqRev";//"seq500_mcf7e2_promoter.fa.SeqRev";//"bg1000.SeqW8S1A1V100.Seq";//

/***********************************************/
PARAM *read_parameters (int nargs, char **argv)
/***********************************************/
{
PARAM *p;
int iarg=1, score=0, i;
static char* syntax = "Pomoda -i inputFasta  [-o outputDIR -w weightFile -ot overlapThreshold -pdt PWM_Divergence_Threshold -ratio minSupportRatio, -FDR p-valueCutoff, -maxlen maxMotifLength, -seedlen SeedLength,-rs Resolution_bp, -mbr min_binding_range, -n numberOfMotifs]";


if(nargs < 3){ printf("%s\n",syntax); exit(0); }
p = new param_st();//(PARAM*)calloc(1,sizeof(PARAM));
p->seedlength= 6;
p->min_supp_ratio=0.05;
p->FDRthresh=1.e-2;
p->olThresh=0.02;
p->pdThresh=0.18;
p->outputDIR=string(".");
p->max_motif_length=26;
p->N_motif=20;
p->resolution=100;
p->startwSize=100;

while(iarg < nargs){
	if(!strcmp(argv[iarg],"-i") && iarg < nargs-1){ p->inputFile=argv[++iarg]; }
	else if(!strcmp(argv[iarg],"-o") && iarg < nargs-1) p->outputDIR=argv[++iarg];
	else if(!strcmp(argv[iarg],"-w") && iarg < nargs-1) p->weightFile=argv[++iarg];
	else if(!strcmp(argv[iarg],"-ot") && iarg < nargs-1) sscanf(argv[++iarg],"%f",&p->olThresh);
	else if(!strcmp(argv[iarg],"-pdt") && iarg < nargs-1) sscanf(argv[++iarg],"%f",&p->pdThresh);
	else if(!strcmp(argv[iarg],"-ratio") && iarg < nargs-1) sscanf(argv[++iarg],"%f",&p->min_supp_ratio);
	else if(!strcmp(argv[iarg],"-FDR") && iarg < nargs-1) sscanf(argv[++iarg],"%f",&p->FDRthresh);
	else if(!strcmp(argv[iarg],"-maxlen") && iarg < nargs-1) sscanf(argv[++iarg],"%d",&p->max_motif_length);
	else if(!strcmp(argv[iarg],"-n") && iarg < nargs-1) sscanf(argv[++iarg],"%d",&p->N_motif);
	else if(!strcmp(argv[iarg],"-seedlen") && iarg < nargs-1) sscanf(argv[++iarg],"%d",&p->seedlength);
	else if(!strcmp(argv[iarg],"-mbr") && iarg < nargs-1) sscanf(argv[++iarg],"%d",&p->startwSize);
	else if(!strcmp(argv[iarg],"-rs") && iarg < nargs-1) sscanf(argv[++iarg],"%d",&p->resolution);
	else{
		printf("ERROR: Wrong option %s\n", argv[iarg]);
		exit(0); 
	}
	++iarg;
}


if(!p->inputFile[0]){ printf("ERROR: %s\n",syntax); exit(0); }


if(p->N_motif<1){ printf("WARNING: Wrong number of output motifs. Using default N=20\n"); p->N_motif=20; }
if(p->min_supp_ratio >1){ printf("WARNING: Minimum ratio should be <1. Using default=0.05\n"); p->min_supp_ratio=0.05; }
return(p);
}

void test(PARAM * setting)
{
		//setting->inputFile="p53_4k.fa";
		
	

		HashEngine engine(setting->inputFile.c_str(),8,setting->weightFile);
		char* q="tctggctcctcatgtgacctctaaaa";
		vector<long> posl;
		engine.searchPatternNRC(q,0,0,posl);
		string g=engine.getSite(61,13);
		MotifModel MM2(&engine,setting->max_motif_length,setting);
		//MM2.AddInstance("GGTCACNSTGAC");//
		MM2.AddInstance("CTTGGCNNNNNNNCA"); //
		MM2.InitializePWMofInstanceSet();
		cout<<MM2.get_consensus()<<endl; 
		MM2.ORScore=1.4;
		MM2.get_consensus(0);
		MotifModel MM3(&engine,setting->max_motif_length,setting);
		//MM3.AddInstance("ACAAAC");  //CAGGACANSNTGNC   GGCAAACAGNNA
		MM3.AddInstance("TGGCACCC");
		//MM3.AddInstance("ACCACGTGC");
		MM3.InitializePWMofInstanceSet();
		MM3.ORScore=2;
		cout<<MM3.get_consensus(0)<<endl;
		vector<MotifModel*> simlist;
		simlist.push_back(&MM2);
		simlist.push_back(&MM3);

	//EM emOpt(setting);
	//emOpt.LoadSeqFile("./result/tber.fa_hotsites.fa");
	//simlist=emOpt.LoadSeedModels(simlist,setting->N_motif);

		MotifModel MM(&engine,setting->max_motif_length,setting);
		MM.seed="TGCCAAG"; //ATGCCC
		MM.AddInstance(MM.seed);
		MM.InitializePWMofInstanceSet();
		cout<<MM.get_consensus(0)<<endl;
		MM.ORScore=2.49284;
		 double temp;
          int bestaln=0;
		  int overlap;
            int aln=MM.AlignmentPWM(&MM, &MM2, temp,overlap);
			MM.MergeList(simlist);

		cout<<	MM.get_consensus(0);
		MM.print();




		double iterCnt=0;
		double improveCount=0;
		iterCnt=improveCount=0;
		double cns=0;
	MM.ComputeScore(0.8,cns,cns,cns,cns,cns);
		for(int j=0;j<MM.X();j++)
		{			
			iterCnt++;
			try
			{
				double improve=MM.updateModel(j);
				if(improve==-1)
				{
					if((iterCnt-improveCount>6))//(MMinst->GetMixedScore())
					break;
				}
				else if(improve!=MINSCORE)
					improveCount++;
				else
					break;

				cout<<MM.get_consensus(0)<<"\t"<<MM.GetMixedScore()<<endl;
			}
			catch( char * str ) {
		     cout << "Exception raised: " << str << '\n';
			 break;
			}

		}
		MM.divergeSeedPart();
			cout<<MM.get_consensus(0)<<"\t"<<MM.GetMixedScore()<<endl;
		getchar();
}

void runAll(PARAM * setting)
{

	int split=setting->inputFile.find_last_of("/\\");
	vector<int> hotsites;
	
	if(split<0)
		split=0;

	string outfilename=setting->outputDIR+"/"+setting->inputFile.substr(split+1)+".pomoda";
	string outfilename2=setting->outputDIR+"/"+setting->inputFile.substr(split+1)+"_hotsites.fa";
	string outfilename3=setting->outputDIR+"/"+setting->inputFile.substr(split+1)+".match";
	cout<<outfilename<<endl;
	ofstream resultout(outfilename.c_str());
	ofstream resultout2(outfilename2.c_str());
	ofstream resultout3(outfilename3.c_str());
	resultout<<"Top "<<setting->N_motif<<" motifs"<<"\t"<<"Window\tScore\tSegment\tCon\tDeg"<<endl;
	int i,j,k;
	VAL start=setStartTime();
	HashEngine engine(setting->inputFile.c_str(),8,setting->weightFile);
	
	HashEngine* engine2=engine.Clone();
	int seqnum=engine.getTotalLength()/engine.SeqLen+1;
	if((setting->min_supp_ratio*seqnum)<50)
		setting->min_supp_ratio=(double)50.0/seqnum;
	cout<<setting->min_supp_ratio<<endl;
	
	MotifModel MM(&engine,setting->max_motif_length,setting);

////////////////////////////////////////////////
	MM.switchFlag=true;
	
	vector<MotifModel*> SeedList=MM.getSeedMotifs(6*setting->N_motif*setting->N_motif,setting->seedlength,setting->min_supp_ratio);
	double markThreshold=SeedList[0]->GetMixedScore();
	cout<<"markThreshold: "<<markThreshold<<endl;
	map<double,MotifModel*> sortlist;
		vector<MotifModel*> markMotifs;
			vector< vector<VAL> > postLists;

	map<MotifModel*,int> filterList;
	map<MotifModel*,vector<MotifModel*> > filterMaps;
		double iterCnt=0;
		double improveCount=0;
		
			
	FOR(i,SeedList.size())
	{		
		cout<<i<<endl;
		MotifModel* MMinst=SeedList[i];
		
		iterCnt=improveCount=0;
		for(int j=0;j<MM.X();j++)
		{			
			iterCnt++;
			try
			{
				double improve=MMinst->updateModel(j);
				if(improve==-1)
				{
					if((iterCnt-improveCount>3))//(MMinst->GetMixedScore())
					break;
				}
				else if(improve!=MINSCORE)
					improveCount++;
				else
					break;
			}
			catch( char * str ) {
		     cout << "Exception raised: " << str << '\n';
			 break;
			}

		}
			
		//	bool ksflag=!KSTest(MMinst->POSLIST,engine.SeqLen);

			if(MMinst->ORScore<1) //||ksflag
			{
				MMinst->POSLIST.clear();
				continue;
			}
			MMinst->PWMRefinement();


		double smallrandom=0.0000000001*sortlist.size();
		if(!isnan(MMinst->GetMixedScore()))
		{
		
				cout<<MMinst->get_consensus()<<endl;
			if(MMinst->GetMixedScore()>markThreshold)
			{				
				MMinst->SearchEngine=engine2;
				MMinst->divergeSeedPart(&hotsites);
				MMinst->SearchEngine=&engine;// need to mark
				if(MMinst->GetMixedScore()>markThreshold&&MMinst->get_consensus().size()>setting->seedlength+1)
				{
					cout<<"Mark!"<<endl;
					MMinst->MarkPos();
				}
					MMinst->SearchEngine=engine2;
				
			}
			else if(MMinst->get_consensus().size()>setting->seedlength+1&&!MMinst->switchFlag)
			{
				MMinst->SearchEngine=&engine;
				MMinst->divergeSeedPart(&hotsites);
				MMinst->SearchEngine=engine2;
			
			}
			else
			{
				int i;
				#ifdef hotsite
				if(MMinst->GetMixedScore()>1)
				FOR(i,MMinst->POSLIST.size())
				{
					long int wpos=MMinst->POSLIST[i];
							bool rc=(wpos<0);
							VAL upos=wpos;
							if(rc)
								upos=0-wpos;
							int seqnum=upos/engine.SeqLen;
							int pos=upos%engine.SeqLen;
					int windowsize=MMinst->BindingRegion;
					double bias=abs(pos-engine.SeqLen/2);
					if(bias<windowsize/2)
					{					
							string pa;
							if(rc)
							{
							
									hotsites.push_back(upos-MMinst->tail);
							}
							else
							{
								
									hotsites.push_back(upos-MMinst->head);
							}
					}

				}
				#endif
				MMinst->InstanceSet.clear();
				
			}


			if(MMinst->GetMixedScore()>1)
			{
				
				sortlist[0-MMinst->GetMixedScore()+ smallrandom]=MMinst;
		

						cout<<MMinst->get_consensus()<<endl;
				cout<<MMinst->GetMixedScore()<<" windowsize:"<<MMinst->BindingRegion<<endl;
						cout<<MMinst->CDScore<<"\t";
				cout<<MMinst->BindingRegion<<"\t";
				//cout<<MMinst->SeqPvalue<<"\t";
				cout<<MMinst->ORScore<<endl;
			}
		}
	}
	cout<<"Elapse Time: "<<getElapsedTime(start)<<endl;
	iterCnt=iterCnt/SeedList.size();
	improveCount=improveCount/SeedList.size();

	
	
#ifdef hotsite
	//handle hotsites
	sort(hotsites.begin(),hotsites.end());
	stringstream headerSS(stringstream::in | stringstream::out);
	stringstream SequenSS(stringstream::in | stringstream::out);
	int lastseq=-1;
	FOR(i,hotsites.size())
	{
		int pos=hotsites[i];
		int pos2=pos;
		FOR(j,hotsites.size()-i-1)
		{
			if(pos2+MM.X()<hotsites[i+j+1])
				break;
			pos2=hotsites[i+j+1];
		}
		if(j==0)
			continue;
		int len=pos2-pos+MM.X();
		string site=engine2->getSite(pos,len);
		int seqnum=pos/engine.SeqLen;
		if(seqnum!=lastseq)
		{
			lastseq=seqnum;
			if(lastseq!=0)
			{	
			resultout2<<">"<<headerSS.str()<<endl;
			resultout2<<SequenSS.str()<<endl;	
		
			headerSS.str("");
			SequenSS.str("");
			}
		}
		headerSS<<"SEQ"<<seqnum<<":"<<pos%engine.SeqLen<<"-"<<pos%engine.SeqLen+len<<"\t"<<j+1<<"\t";
		SequenSS<<site<<"X"; //split different hot sites
		i=i+j;
	}
	resultout2<<">"<<headerSS.str()<<endl;
	resultout2<<SequenSS.str()<<endl;
	resultout2.close();
#endif

	//return ;
	map<double,MotifModel*>::iterator ITER;

	do{
		 ITER=sortlist.begin();
		 markMotifs.clear();
	i=0;
	while(ITER!=sortlist.end())//&&i<setting->N_motif
	{
		MotifModel* MMinst=ITER->second;

		int SEQLEN=engine2->SeqLen;
		MMinst->SearchEngine=engine2;
		
		

		
		vector<VAL> templist=MMinst->getMatchPos();
		
		int j;
		double minscore=0;
		
		int filterId=-1;

		double maxalnscore=-MINSCORE;

		minscore=0;//////
		double threshold=(double)engine.SeqLen/engine.TotalLength;
		if( (MMinst->SeqPvalue>setting->FDRthresh ||MMinst->get_consensus(0).size()<setting->seedlength))// 0.000000001  0.000001  0.00000001 
		{
			ITER++;
			filterList[MMinst]=6000;
			continue;
		}

		int mergeid=0;
		FOR(j,postLists.size())
		{

			int comcount=0;
				try
				{
					double alnscore;
					double bestas;
					int alnn=MMinst->AlignmentPWMRC(MMinst,markMotifs[j],bestas);
					alnscore=(double)bestas;//min(MMinst->Consensus.size(),markMotifs[j]->Consensus.size());
					double alterP=1/(pow(4.0,MMinst->Length()*(1-alnscore)));
					double score=MMinst->SimilarityScore(postLists[j],templist,markMotifs[j]->Consensus.size(),MMinst->Consensus.size(),comcount,alterP);
				score=(double)comcount/min(templist.size(),postLists[j].size());////////
					//if(markMotifs[j]->ORScore>markThreshold)
					//{
					//	if(maxalnscore<0.18)
					//	{
					//		filterId=1000+maxalnscore*100;
					//		maxalnscore=0.01;					
					//		MMinst->seed+="-"+markMotifs[j]->seed;
					//		//merge
					//		mergeid=j;
					//		//filterMaps[markMotifs[j]].push_back(MMinst);
					//	}
					//}
					if(minscore<score)/////
					{
						minscore=score;
						if(minscore>setting->olThresh)//<0.0001
						{
							filterId=2000+minscore*100;
							MMinst->seed+="-"+markMotifs[j]->seed;

							//merge
							mergeid=j;
							//filterMaps[markMotifs[j]].push_back(MMinst);
						}
					}
					if(maxalnscore>alnscore)
					{
						maxalnscore=alnscore;
						if(maxalnscore<setting->pdThresh)
						{
							filterId=1000+maxalnscore*100;
							MMinst->seed+="-"+markMotifs[j]->seed;

							//merge
							mergeid=j;
							//filterMaps[markMotifs[j]].push_back(MMinst);
						}
					}
				}
				catch( char * str ) 
				{
					 cout << "Exception raised: " << str << '\n';
					 break;
				}

		}


		if(minscore>setting->olThresh||maxalnscore<setting->pdThresh) //minscore<0.00000001
		{
			ITER++;
			filterList[MMinst]=filterId;
			int mainId=filterId%1000;	
			filterMaps[markMotifs[mergeid]].push_back(MMinst);
			if(DEBUG)
			cout<<"filter:"<<MMinst->get_consensus(0)<<"\t"<<(0-ITER->first)<<endl;
			
			continue;
		}



/**************************second pass filtering ****************************/
			//map<MotifModel*,int>::iterator ITER22=filterList.begin();
			//maxalnscore=-MINSCORE;minscore=0;
			//while(ITER22!=filterList.end())
			//{
			//	if(ITER22->first->ORScore<1.5||ITER22->second>=6000)//
			//	{
			//		ITER22++;
			//		continue;				
			//	}

			//	double alnscore;
			//		double bestas;
			//		
			//		int alnn=MMinst->AlignmentPWMRC(MMinst,ITER22->first,bestas);
			//		alnscore=(double)bestas;
			//	int comcount=0;
			//	double score=MMinst->SimilarityScore(ITER22->first->POSLIST,templist,ITER22->first->Consensus.size(),MMinst->Consensus.size(),comcount);
			//	score=(double)comcount/min(templist.size(),ITER22->first->POSLIST.size());////////

			//	if(minscore<score)/////
			//		{
			//			minscore=score;
			//			if(minscore>0.4)//<0.0001
			//			{
			//				filterId=4000+minscore*100;
			//				MMinst->seed+="-"+ITER22->first->seed;
			//				//merge
			//				filterMaps[ITER22->first].push_back(MMinst);
			//				break;
			//			}
			//		}
			//		if(maxalnscore>alnscore)
			//		{
			//			maxalnscore=alnscore;
			//			if(maxalnscore<0.14)
			//			{
			//				filterId=5000+maxalnscore*100;
			//				MMinst->seed+="-"+ITER22->first->seed;
			//				//merge
			//				filterMaps[ITER22->first].push_back(MMinst);
			//				break;
			//			}
			//		}
			//	ITER22++;
			//}

			//if(minscore>0.4||maxalnscore<0.17) //minscore<0.00000001
			//{
			//	ITER++;
			//	filterList[MMinst]=filterId;
			//	if(DEBUG)
			//	cout<<"filter:"<<MMinst->get_consensus(0)<<"\t"<<(0-ITER->first)<<endl;
			//	//delete MMinst;
			//	continue;
			//}

/**************************second pass filtering ****************************/

		double snr,prbE,prbN;

		if(MMinst->ORScore>1)
		{
		postLists.push_back(MMinst->POSLIST);
		
		markMotifs.push_back(MMinst);
			cout<<"Motif "<<i<<":";
			
			cout<<MMinst->seed<<":"<<MMinst->get_consensus(0)<<endl;
		cout<<(0-ITER->first) <<"\t";
		cout<<MMinst->CDScore<<"\t";
		cout<<MMinst->BindingRegion<<"\t";
		cout<<MMinst->SeqPvalue<<"\t";
		cout<<MMinst->ORScore<<endl;
		string name="Motif"+Int2String(i);
		i++;
		}
		 //MMinst->printPWM(name);
		// MMinst->printRealPos(name);
		ITER++;
		

		filterMaps[MMinst]=vector<MotifModel*>();
		}
		SeedList.clear();


		

		break;
		//if(markMotifs.size()==sortlist.size())

			
		/***************Merge the filtered Motifs*************************/
		FOR(k,markMotifs.size())
		{
				MotifModel* MMinst=markMotifs[k];
				MMinst->MergeList(filterMaps[MMinst]);
				SeedList.push_back(MMinst);
				MMinst->PWMRefinement();
		}
		/***************Merge the filtered Motifs*************************/
		filterMaps.clear();
		postLists.clear();
		sortlist.clear();
		
	
///***************Generalization Phase*****************/
	//EM emOpt(setting);
	//emOpt.LoadSeqFile(outfilename2);//("./result/stber300.fa");
	//markMotifs=emOpt.LoadSeedModels(markMotifs,setting->N_motif);
	sortlist.clear();
	FOR(k,markMotifs.size())
	{
		MotifModel* MMinst=markMotifs[k];
		MMinst->SearchEngine=engine2;
		MMinst->ComputeScore(0.8,MMinst->CDScore,MMinst->ORScore,MMinst->BindingRegion,MMinst->CNSVScore,MMinst->DiffScore);
		if(MMinst->ORScore==MINSCORE)
			MMinst->ComputeScore(0.99,MMinst->CDScore,MMinst->ORScore,MMinst->BindingRegion,MMinst->CNSVScore,MMinst->DiffScore);
		if(MMinst->ORScore<1)
		{
			cout<<"oop! "<<MMinst->Consensus<<MMinst->ORScore<<endl;
			continue;
		}
		double smallrandom=0.0000000001*sortlist.size();
		sortlist[0-MMinst->ORScore+smallrandom]=MMinst;
	}
		
		
		markMotifs.clear();
		
	}while(true);

 //    ITER=sortlist.begin();
	//	k=0;
	//while(ITER!=sortlist.end())
	//{
	//	if(0-ITER->first<1)
	//		break;
	//	MotifModel* MMinst=ITER->second;
	//	vector<VAL> templist=MMinst->POSLIST;
	//	
	//	int j;
	//	double minscore=0;
 //       double maxalnscore=-MINSCORE;
	//	int filterId=-1;
	//	FOR(j,markMotifs.size())
	//	{

	//		int comcount=0;
	//			try
	//			{
	//				double alnscore;
	//				double bestas;
	//				int alnn=MMinst->AlignmentPWMRC(MMinst,markMotifs[j],bestas);
	//				alnscore=(double)bestas;//min(MMinst->Consensus.size(),markMotifs[j]->Consensus.size());
	//				double alterP=1/(pow(4.0,MMinst->Length()*(1-alnscore)));
	//				double score=MMinst->SimilarityScore(postLists[j],templist,markMotifs[j]->Consensus.size(),MMinst->Consensus.size(),comcount,alterP);
	//			score=(double)comcount/min(templist.size(),postLists[j].size());////////

	//				if(minscore<score)/////
	//				{
	//					minscore=score;
	//					if(minscore>setting->olThresh)//<0.0001
	//					{
	//						filterId=9000+minscore*100;
	//						MMinst->seed+="-"+markMotifs[j]->seed;

	//						//merge
	//						filterMaps[markMotifs[j]].push_back(MMinst);
	//					}
	//				}
	//				if(maxalnscore>alnscore)
	//				{
	//					maxalnscore=alnscore;
	//					if(maxalnscore<setting->pdThresh)
	//					{
	//						filterId=8000+maxalnscore*100;
	//						MMinst->seed+="-"+markMotifs[j]->seed;
	//					}
	//				}
	//			}
	//			catch( char * str ) 
	//			{
	//				 cout << "Exception raised: " << str << '\n';
	//				 break;
	//			}

	//	}


	//	if(minscore>setting->olThresh||maxalnscore<setting->pdThresh) //minscore<0.00000001
	//	{
	//		ITER++;
	//		filterList[MMinst]=filterId;
	//		int mainId=filterId%1000;	
	//		if(DEBUG)
	//		cout<<"filter:"<<MMinst->get_consensus(0)<<"\t"<<(0-ITER->first)<<endl;
	//		
	//		continue;
	//	}
	//	markMotifs.push_back(MMinst);
	//	ITER++;
	//}
///***************Generalization Phase*****************/

/***************output distance array*************************/
	string outfilename4=setting->outputDIR+"/"+setting->inputFile.substr(split+1)+".dist";
	string outfilename5=setting->outputDIR+"/"+setting->inputFile.substr(split+1)+".name";
	ofstream resultout4(outfilename4.c_str());
	ofstream resultout5(outfilename5.c_str());
	FOR(k,markMotifs.size()-1)
	{
		MotifModel* MMinst1=markMotifs[k];
		for(j=k+1;j<markMotifs.size();j++)
		{
			MotifModel* MMinst2=markMotifs[j];
				int comcount=0;
					double alnscore;
					double bestas;
					int alnn=MMinst1->AlignmentPWMRC(MMinst1,MMinst2,bestas);
					alnscore=(double)bestas;//min(MMinst->Consensus.size(),markMotifs[j]->Consensus.size());
					double alterP=1/(pow(4.0,MMinst1->Length()*(1-alnscore)));
					double score=MMinst1->SimilarityScore(MMinst1->POSLIST,MMinst2->POSLIST,MMinst1->Consensus.size(),MMinst2->Consensus.size(),comcount,alterP);
				    score=(double)comcount/min(MMinst1->POSLIST.size(),MMinst2->POSLIST.size());////////
					resultout4<<1-score<<"\t"<<alnscore<<endl;
		}
		resultout5<<MMinst1->Consensus<<"\t"<<MMinst1->ORScore<<endl;
	
	}
	resultout5<<markMotifs[k]->Consensus<<"\t"<<markMotifs[k]->ORScore<<endl;
	resultout4.close();
	resultout5.close();
/***************output distance array*************************/


		FOR(k,markMotifs.size())
		{
				MotifModel* MMinst=markMotifs[k];
				//cout<<"top "<<k<<":\t"<<MMinst->get_consensus(0)<<endl;
				////MMinst->PWMRefinement();
				//cout<<"top "<<k<<":\t"<<MMinst->get_consensus(0)<<endl;
			//	MMinst->MergeList(filterMaps[MMinst]);
			stringstream bindingsites(stringstream::in | stringstream::out);
#ifdef hotsite
			FOR(j,MMinst->POSLIST.size())
			{
				int SEQLEN=10000;
				long int wpos=MMinst->POSLIST[j];
				if(wpos==-1)
					continue;
				bool rc=(wpos<0);
				VAL upos=wpos;
				
				if(rc)
					upos=0-wpos;
				int seqnum=upos/SEQLEN;
				int pos=upos%SEQLEN;
				int windowsize=MMinst->BindingRegion;
				double bias=abs(pos-SEQLEN/2);
				int motiflen=MMinst->Length();
				//if(bias<=windowsize/2)
				{
					bindingsites<<"Motif_"<<k<<"\t"<<seqnum<<"\t"; //start from 1
					//if(rc)
					//	bindingsites<<"b";
					//else
					//	bindingsites<<"f";
					bindingsites<<pos<<"\t";
					bindingsites<<pos+motiflen<<endl;		
					//bindingsites<<MMinst->SearchEngine->getSite(upos,MMinst->Length())<<endl;								
				}
			}
#endif
			
			resultout<<"DE	Motif_"<<k<<"\t";
			resultout<<MMinst->BindingRegion<<"\t";
			resultout<<MMinst->ORScore<<"\t";
			resultout<<MMinst->CDScore<<"\t";
			resultout<<MMinst->Consensus<<"\t";
			resultout<<reverseString( MMinst->Consensus)<<endl;
			resultout<<"PO	A	C	G	T"<<endl;
			for(int kkk=MMinst->head;kkk<MMinst->X()-MMinst->tail;kkk++)
			{
				resultout<<kkk-MMinst->head+1<<"\t";
				resultout<<MMinst->g(kkk,0)<<"\t";
				resultout<<MMinst->g(kkk,1)<<"\t";
				resultout<<MMinst->g(kkk,2)<<"\t";
				resultout<<MMinst->g(kkk,3)<<endl;

			}
			resultout<<"XX"<<endl;
			resultout<<endl;
			resultout3<<bindingsites.str();
			resultout<<endl;

	}

	map<MotifModel*,int>::iterator ITER2=filterList.begin();
	cout<<"Filter List:"<<endl;
	while(ITER2!=filterList.end())
	{
		cout<<ITER2->first->seed<<":"<<ITER2->first->get_consensus(0)<<"\t"<<(ITER2->first->GetMixedScore())<<"\tfilter"<<ITER2->second<<endl;
		ITER2++;
	}
	cout<<"Elapse Time: "<<getElapsedTime(start)<<endl;
	resultout<<"Elapse Time: "<<getElapsedTime(start)<<endl;
	resultout.close();
	resultout3.close();

	

	
}



int main(int argc, char* argv[])
{
	//std::ofstream log("result/log.txt");
 //   std::streambuf *oldbuf = std::cout.rdbuf(log.rdbuf());

	PARAM * setting=read_parameters (argc, argv);
	unsigned long int aa=1;
	cout<<sizeof(aa)*8<<endl;


	//test(setting);
	runAll(setting);

	string tag;
	//cin>>tag;
	return 1;
	
	
}