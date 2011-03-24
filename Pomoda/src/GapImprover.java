import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map.Entry;
import java.util.TreeMap;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.biojava.bio.dist.DistributionTools;
import org.biojava.bio.dist.UniformDistribution;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;

import cern.colt.Arrays;
import auc.AUCCalculator;
import auc.Confusion;


public class GapImprover {

	/**
	 * @param args
	 */
	public String outputPrefix="./";
	public String inputFasta;
	public String ctrlFasta="";
	public boolean OOPS=false; //only one dependence per sequence
	public boolean OOPG=false; //only one occurrence per sequence
	public boolean removeBG=false; //false:uniform BG assume
	public String bgmodelFile="";
	LinearEngine SearchEngine;
	public double sampling_ratio=1;
	public double FDR=0.01;
	public double entropyThresh=1;
	public int FlankLen=0;
	public int max_gaplen=12;
	
	public BGModel background;
	public GapImprover(Pomoda motiffinder)
	{
		SearchEngine=motiffinder.SearchEngine2;
		sampling_ratio=motiffinder.sampling_ratio;
		FDR=motiffinder.FDR;
		background=motiffinder.background;
		
	}
	
	public GapImprover()
	{
		
	}
	
	public double KL_Divergence_empirical(List<String> sites,PWM motif,int start)
	{
		
		double KLsum=0;
		HashMap<String, Integer> sitecount=new HashMap<String, Integer>();
		Iterator<String> iter=sites.iterator();
		int totalCount=0;
		while(iter.hasNext())
		{
			String temp=iter.next();
			temp=temp.substring(start,start+motif.core_motiflen);
		
			if(temp.contains("N"))
				continue;
			if(sitecount.containsKey(temp))
			{
				sitecount.put(temp, sitecount.get(temp)+1);
			}
			else
			{
				sitecount.put(temp, 1);
			}
			totalCount++;
		}
		
		for(String key:sitecount.keySet())
		{
			int count=sitecount.get(key);
			double p=(double)count/totalCount;
			double motif_logP=motif.scoreWeightMatrix(key);
			if(Double.isInfinite(motif_logP))
				continue;
		
			KLsum+=p*(Math.log(p)-motif_logP);
		}
			
			
		return KLsum;
		
	}
	
	public void initialize()
	{

		SearchEngine=new LinearEngine(4);
		SearchEngine.build_index(this.inputFasta);
	
		background=new BGModel();
		File file=null;
		int bg_markov_order=3;
		if(ctrlFasta.isEmpty())
		{
			bg_markov_order=1;
			file= new File(inputFasta+".bgobj");
		}
		else
			file= new File(ctrlFasta+".bgobj");
		
		if(file.exists()||!bgmodelFile.isEmpty())
		{
			if(bgmodelFile.isEmpty())
			background.LoadModel(file.getAbsolutePath());
			else
				background.LoadModel(bgmodelFile);
		}
		else
		{
			if(ctrlFasta.isEmpty())
			{
		     background.BuildModel(inputFasta, bg_markov_order+1); //3-order bg
		     background.SaveModel(inputFasta+".bgobj");
			}
			else
			{
			     background.BuildModel(ctrlFasta, bg_markov_order+1); //3-order bg
			     background.SaveModel(ctrlFasta+".bgobj");
				}
				
		}

		
	}
	
	
	public HashMap<HashSet<Integer>,HashMap<String,Double>> FindBest(List<GapBGModelingThread> list)
	{

		HashMap<HashSet<Integer>,HashMap<String,Double>> Dmap=new HashMap<HashSet<Integer>,HashMap<String,Double>>();
		Iterator<GapBGModelingThread> iter3=list.iterator();
		int start=-1;
		double minKL=Double.MAX_VALUE;
		double baseScore=0;
		double bestScore=0;
		GapBGModelingThread bestThread=null;
		HashSet<Integer> bestDgroups=null;
		try {

		while(iter3.hasNext())
		{
			GapBGModelingThread t1=iter3.next();	
			t1.join();
				if(t1.depend_Pos.size()==0)
				{
					System.out.println(t1.toString());
					baseScore=t1.KL_Divergence;
					
					
				}
		}
		iter3=list.iterator();	
		//filter the negative threads
		ArrayList<GapBGModelingThread> positiveThread=new ArrayList<GapBGModelingThread>(list.size()/2);
		
		 LinkedList< HashSet<Integer> > queue=new LinkedList< HashSet<Integer> >();
		for(GapBGModelingThread t2:list)
		{
			if(t2.KL_Divergence<baseScore)
			{
				t2.KL_Divergence=baseScore-t2.KL_Divergence;
				positiveThread.add(t2);
				if(t2.KL_Divergence>bestScore)
				{
					bestDgroups=new HashSet<Integer>();
					bestDgroups.add(positiveThread.size()-1);
					bestScore=t2.KL_Divergence;
				}
				
			}
		}
		
		// build graph
		
		HashSet[] adjgraph=new HashSet[positiveThread.size()];
		for (int i = 0; i < positiveThread.size()-1; i++) {
			HashSet<Integer> s1=positiveThread.get(i).depend_Pos;
			for (int j = i+1; j < positiveThread.size(); j++) 
			{
				HashSet<Integer> s2=positiveThread.get(j).depend_Pos;
				boolean overlap=false;
				HashSet<Integer> intersect=new HashSet<Integer>(s2);
				intersect.retainAll(s1);
				if(intersect.size()>0)
					overlap=true;
				
				if(overlap==false)
				{
					if(adjgraph[i]==null)
						adjgraph[i]=new HashSet<Integer>();
					if(adjgraph[j]==null)
						adjgraph[j]=new HashSet<Integer>();
					adjgraph[i].add(j);
					adjgraph[j].add(i);
					HashSet<Integer> set=new HashSet<Integer>();
					set.add(i);
					set.add(j);
					queue.add(set);
				}
			}
		}
		
		//enumerate all possible number of group
		 HashSet<Integer> topElm=null;
		 
		while((topElm=queue.poll())!=null)
		{
			double currscore=0;
			Iterator<Integer> iter=topElm.iterator();
			HashSet<Integer> commonThirdPoint=null;
			
			while(iter.hasNext())
			{
				int id=iter.next();
				currscore+=positiveThread.get(id).KL_Divergence;
				if(commonThirdPoint==null)
				{
					commonThirdPoint=new HashSet<Integer>(adjgraph[id]);
				}
				else
				{
					commonThirdPoint.retainAll(adjgraph[id]);
				}
				
			}
			//push queue
			iter=commonThirdPoint.iterator();
			while(iter.hasNext())
			{
				HashSet<Integer> newCombine=new HashSet<Integer>(topElm);
				newCombine.add(iter.next());
				queue.add(newCombine);
			}
			if(currscore>bestScore)
			{
				bestDgroups=topElm;
				bestScore=currscore;
			}
			
		}
		
		if(bestDgroups!=null)
		for(Integer id : bestDgroups)
		{
			System.out.println("best:"+positiveThread.get(id).toString());
			Dmap.put(positiveThread.get(id).depend_Pos,positiveThread.get(id).DprobMap);
		}
				
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return Dmap;
	}
	public GapPWM fillDependency(PWM motif)
	{
		GapPWM gapPWM=null;
		String Consensus=motif.Consensus(true);
		ArrayList<Integer> gapstart=new ArrayList<Integer>(motif.core_motiflen/2);
		ArrayList<Integer> gapend=new ArrayList<Integer>(motif.core_motiflen/2);
		
		//detect gap range
		int start=-1;
		for (int i = motif.head-FlankLen; i < motif.head+motif.core_motiflen+FlankLen; i++) {
			double entropy=0;
			if(i<0||i>=motif.columns())
				entropy=2;
			else
				entropy=DistributionTools.totalEntropy(motif.getColumn(i)) ;//
			if(entropy>entropyThresh&&(max_gaplen>(i-motif.head+FlankLen-start)||start==-1))
			{
				if(start==-1)
				{
					start=i-motif.head+FlankLen;
				}
			}
			else
			{
				if(start!=-1)
				{
				gapstart.add(start);
				gapend.add(i-motif.head+FlankLen);
				}
				start=-1;
//				if(i-motif.head+FlankLen>18)
//				break;
			}
		}
		if(start!=-1)
		{
			gapstart.add(start);
			gapend.add(motif.core_motiflen+2*FlankLen);
		}
		
		//debug
	//	ArrayList<Double> snull = new ArrayList<Double>(),sbest =new ArrayList<Double>();
		
		//get a set of instance strings
		LinkedList<String> sites=new LinkedList<String>();
		if(!OOPS)
		{
			double pwmThresh=motif.getThresh(sampling_ratio, FDR, background);
			LinkedList<FastaLocation> falocs=SearchEngine.searchPattern(motif, pwmThresh);
			Iterator<FastaLocation> iter=falocs.iterator();
			while(iter.hasNext())
			{
				FastaLocation currloc=iter.next();
				String site=SearchEngine.getSite(currloc.getSeqId(), currloc.getSeqPos()-FlankLen, motif.core_motiflen+2*FlankLen);
				if(currloc.ReverseStrand)
					site=common.getReverseCompletementString(site);
				if(site!=null)
					sites.add(site);		
			}
		}
		else
		{
			 LinkedList<FastaLocation> falocs =null;
			 if(removeBG)
				 falocs=SearchEngine.searchPattern(motif, 0); //enable background in searching
			 else
				 falocs=SearchEngine.searchPattern(motif, Double.NEGATIVE_INFINITY);
        	 Iterator<FastaLocation> iter=falocs.iterator();
        	 int lastseq=-1;
        	 double seqcount=0;
        	 double maxseq_score=	Double.NEGATIVE_INFINITY;
        	 FastaLocation max_currloc=null;
        	 while(iter.hasNext())
        	 {
        		 FastaLocation currloc=iter.next();
        		 if(lastseq!=currloc.getSeqId())
        		 {
        			 seqcount+=1;
        			
        			 if(lastseq!=-1)
        			 {
     				String site=SearchEngine.getSite(max_currloc.getSeqId(), max_currloc.getSeqPos()-FlankLen, motif.core_motiflen+2*FlankLen);
    				if(max_currloc.ReverseStrand)
    					site=common.getReverseCompletementString(site);
    				if(site!=null)
    					sites.add(site);
    				
    				//debug
    				//snull.add(max_currloc.Score+Math.log(0.25));
        			 }
        			 lastseq=currloc.getSeqId(); 
        			 maxseq_score=currloc.Score;
        			 max_currloc=currloc;
        		 }
        		 if(maxseq_score<currloc.Score)
        		 {
        			 maxseq_score=currloc.Score;
        			 max_currloc=currloc;
        		 }
        	 }
        	 //last seq
 			String site=SearchEngine.getSite(max_currloc.getSeqId(), max_currloc.getSeqPos()-FlankLen, motif.core_motiflen+2*FlankLen);
			if(max_currloc.ReverseStrand)
				site=common.getReverseCompletementString(site);
			if(site!=null)
			sites.add(site);
			//debug
			//snull.add(max_currloc.Score+Math.log(0.25));
		}
		try {
		//find the best dependency modeling in each gap region
		LinkedList<GapBGModelingThread> threadPool=new LinkedList<GapBGModelingThread>();
		HashMap<HashSet<Integer>,HashMap<String,Double>> Dmap=new HashMap<HashSet<Integer>,HashMap<String,Double>>();
		if(sites.size()>0)
		{
		
		for (int i = 0; i <gapstart.size(); i++) {
			int gstart=gapstart.get(i);
			int gend=gapend.get(i);
			if(gend-gstart<2)
				continue;
			int combinNum=1<<(gend-gstart);
			for (int j = 0; j < combinNum; j++) {
				HashSet<Integer> dpos=new HashSet<Integer>();
				int bcode=j;
				for (int j2 = 0; j2 < gend-gstart; j2++) {
					if(bcode%2==0)
						dpos.add(j2+gstart);
					bcode>>=1;
				}
				if(dpos.size()==1)
					continue;
				GapBGModelingThread t1=null;
				if(removeBG)
					t1=new GapBGModelingThread(gstart, gend, sites, dpos,background);//null mean not considering BG
				else
					t1=new GapBGModelingThread(gstart, gend, sites, dpos,null);//null mean not considering BG
				t1.run();
				threadPool.add(t1);
			}
			if(!OOPG)
			{
			Dmap.putAll( FindBest(threadPool.subList(0, threadPool.size())));
			threadPool.clear();
			}
		}
		
		}
	
		Iterator<GapBGModelingThread> iter3=threadPool.iterator();
		 start=-1;
		double minKL=Double.MAX_VALUE;
		GapBGModelingThread bestThread=null;
		
		GapBGModelingThread[] Pos_BestThread=new GapBGModelingThread[motif.core_motiflen+2*FlankLen];
		
		if(OOPG)
		while(iter3.hasNext())
		{
			GapBGModelingThread t1=iter3.next();	
				t1.join();
				if(t1.depend_Pos.size()==0)
				{
					System.out.println(t1.toString());
					//snull=t1.debuglist;
				}
				if(t1.gapStart!=start)
				{
					start=t1.gapStart;
					minKL=t1.KL_Divergence;
					if(bestThread!=null)
					{
						System.out.println("best:"+bestThread.toString());
						if(bestThread.depend_Pos.size()>1)
						{
							Dmap.put(bestThread.depend_Pos, bestThread.DprobMap);
							
						}
					}
					bestThread=t1;
				}
				else 
				{
					if(t1.KL_Divergence<minKL)
					{
						minKL=t1.KL_Divergence;
						bestThread=t1;
					}
					
				}
				
				Iterator<Integer> iter4=t1.depend_Pos.iterator();
				while(iter4.hasNext())
				{
					int posId=iter4.next();
					if(Pos_BestThread[posId]!=null)
					{
						if(Pos_BestThread[posId].KL_Divergence>t1.KL_Divergence)
							Pos_BestThread[posId]=t1;
					}
					else
					{
						Pos_BestThread[posId]=t1;
					}
				}
				if(t1.depend_Pos.size()==0)
				{
					for (int i = t1.gapStart; i < t1.gapEnd; i++) {
						if(Pos_BestThread[i].KL_Divergence>t1.KL_Divergence)
							Pos_BestThread[i]=t1;
						
					}
				}
		}
		if(bestThread!=null&&bestThread.depend_Pos.size()>1)
		{
			Dmap.put(bestThread.depend_Pos, bestThread.DprobMap);
		System.out.println("best:"+bestThread.toString());
		//sbest=bestThread.debuglist;
		}

		
//		if(!OOPG)
//		{
//		//fill in multi-dependency in the same gap region	
//		for (int i = 0; i < Pos_BestThread.length; i++) {
//			if(Pos_BestThread[i]!=null&&Pos_BestThread[i].depend_Pos.size()>1)
//			{
//				boolean maxCover=true;
//				for(Integer ii:Pos_BestThread[i].depend_Pos)
//				{
//					if(Pos_BestThread[ii].hashCode()!=Pos_BestThread[i].hashCode())
//					{
//						maxCover=false;
//						break;
//					}
//				}
//				if(maxCover)
//				{
//				Dmap.put(Pos_BestThread[i].depend_Pos, Pos_BestThread[i].DprobMap);
//				System.out.println(Pos_BestThread[i]);
//				}
//			}
//		}
//		}
	
		gapPWM=GapPWM.createGapPWM(motif.subPWM( motif.head,motif.head+motif.core_motiflen), Dmap,FlankLen);
		if(gapPWM.core_motiflen!=motif.core_motiflen&&sites.size()>0)
		{
			motif=new PWM(sites.toArray(new String[1]));
			gapPWM=GapPWM.createGapPWM(motif.subPWM(FlankLen,motif.columns()-FlankLen), Dmap,FlankLen);
			motif=motif.subPWM(gapPWM.head, gapPWM.head+gapPWM.core_motiflen);
		}
		
		//debug
//		for(String site : sites)
//		{
//			double gs=gapPWM.scoreWeightMatrix(site.substring(2,19));
//			sbest.add(gs);
//		}
//		int scnt=0;
//		for (int i = 0; i < sbest.size(); i++) {
//			if(sbest.get(i)>snull.get(i))
//				scnt+=1;
//		}
		
		
//		String testStr="agagaagaagaaagaaagaaagaagaaaggaagaaagaaagaaagaaagaaagaaagaaagaaagaaagaaagaaagaaagaaagaaaagaaagaaagaa";
//		testStr=common.getReverseCompletementString(testStr);
//		for (int i = 0; i < testStr.length()-gapPWM.core_motiflen; i++) {
//			String temp=testStr.substring(i, i+gapPWM.core_motiflen);
//			double score1=gapPWM.scoreWeightMatrix(temp);
//			double score2=motif.scoreWeightMatrix(temp);
//			System.out.println(score1+"\t"+score2);
//		}
		System.out.println(KL_Divergence_empirical(sites, motif,gapPWM.head) +"\t"+KL_Divergence_empirical(sites, gapPWM,gapPWM.head));
		
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IllegalAlphabetException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IllegalSymbolException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		
		return gapPWM;
	}
	
	
	public double AUCtest(PWM motif)
	{
      
		
         LinearEngine BGSearch=new LinearEngine(6);
         
         Iterator<String> iter2=SearchEngine.ForwardStrand.iterator();
         background.r.setSeed(0);
         while(iter2.hasNext())
         {
        	 int len=iter2.next().length();
        	 String bgstr="";
        	 if(removeBG)
        	 {
	        	 KeyValuePair<Double, String> bgstr_p=background.generateRandomSequence(len);
	        	 bgstr=bgstr_p.value;
        	 }
        	 else
        	 {
	        	 UniformDistribution ud=new UniformDistribution(DNATools.getDNA());
	        	 bgstr=DistributionTools.generateSymbolList(ud, len).seqString();
        	 }
        	 BGSearch.ForwardStrand.add(bgstr);	 
        	 BGSearch.TotalLen+=bgstr.length();
         }

     	TreeMap<Double,Integer> Sorted_labels=new TreeMap<Double,Integer>();
         
        	 LinkedList<FastaLocation> falocs =SearchEngine.searchPattern(motif, Double.NEGATIVE_INFINITY);
        	 Iterator<FastaLocation> iter=falocs.iterator();
        	 int lastseq=-1;
        	 double seqcount=0;
        	 double maxseq_score=	Double.NEGATIVE_INFINITY;
        	 while(iter.hasNext())
        	 {
        		 FastaLocation currloc=iter.next();
        		 if(lastseq!=currloc.getSeqId())
        		 {
        			 seqcount+=1;
        			 if(lastseq!=-1)
        			 {
        				 Sorted_labels.put(maxseq_score+seqcount*common.DoubleMinNormal, 1);
        			 }
        				 lastseq=currloc.getSeqId();
        			 maxseq_score=currloc.Score;
        		 }
        		 if(maxseq_score<currloc.Score)
        		 {
        			 maxseq_score=currloc.Score;

        		 }
        	 }
        	 Sorted_labels.put(maxseq_score+seqcount*common.DoubleMinNormal, 1);
        	 
        	 //bg sequences
        	 falocs =BGSearch.searchPattern(motif, Double.NEGATIVE_INFINITY);
        	 iter=falocs.iterator();
        	lastseq=-1;
        	 seqcount=0;
        	maxseq_score=Double.NEGATIVE_INFINITY;
	       	 while(iter.hasNext())
	    	 {
	    		 FastaLocation currloc=iter.next();

	    		 if(lastseq!=currloc.getSeqId())
	    		 {
	    			 seqcount+=1;
	    			
	    			 if(lastseq!=-1)
	    			 {
	    			 Sorted_labels.put(maxseq_score-seqcount*common.DoubleMinNormal, 0);
	    			
	    			 }
	    			 maxseq_score=currloc.Score;
	    			 lastseq=currloc.getSeqId();
	    		 }
	    		 if(maxseq_score<currloc.Score)
	    		 {
	    			 maxseq_score=currloc.Score;
	    		 }
	    		 
	    	 }
	       	Sorted_labels.put(maxseq_score-seqcount*common.DoubleMinNormal, 0);
	       	 int[]  labels=new int[Sorted_labels.size()];
	       	double[]  scores=new double[Sorted_labels.size()];
	       	 int ii=0;
	       	 int one=0;
	       	 for(Double key:Sorted_labels.descendingKeySet())
	       	 {
	       		 labels[ii]=Sorted_labels.get(key);
	       		 if(labels[ii]==1)
	       			 one++;
	       		 scores[ii]=key;
	       		        ii++;
	       	 }
        	 
        	 //AUCcalc.addROCPoint(fp,(double)seqcount/SearchEngine.getSeqNum());
	       	Confusion AUCcalc=AUCCalculator.readArrays(labels, scores);	
	       	System.out.println(scores[scores.length-1]);
         double AUCscore=AUCcalc.calculateAUCROC();
         
         return AUCscore;
		
	}
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		Options options = new Options();
		options.addOption("i", true, "input fasta file");
		options.addOption("pwm", true, "input PWM file");
		options.addOption("c", true, "control fasta file");
		options.addOption("bgmodel", true, "background model file");
		options.addOption("prefix", true, "output directory");
		options.addOption("flank", true, "the numboer of flanking positions around PWM to include(default 0)");
		options.addOption("ratio",true, "sampling ratio (default 1)");
		options.addOption("thresh",true, "minimum entropy threshold for considering a position as a gap(default 1)");
		options.addOption("oops",false,"whether assuming only one occurrence per sequence (default false)");
		options.addOption("oopg",false,"whether assuming only one dependence per gap region (default false)");
		options.addOption("rmbg",false,"whether considering the background probility in learning and evaluating the dependences (default false)");
		options.addOption("maxlen",true,"maxmimum length of gap (default 12)");
		options.addOption("FDR",true,"fasle positive rate");
		String inputPWM;
		CommandLineParser parser = new GnuParser();
		GapImprover GImprover=new GapImprover();
		
		try {
			CommandLine cmd = parser.parse( options, args);
			if(cmd.hasOption("i"))
			{
				GImprover.inputFasta=cmd.getOptionValue("i");
			}
			else
			{
				throw new ParseException("no input fasta file");
			}
			if(cmd.hasOption("pwm"))
			{
				inputPWM=cmd.getOptionValue("pwm");
			}
			else
			{
				throw new ParseException("no input pwm file");
			}
			if(cmd.hasOption("c"))
			{
				GImprover.ctrlFasta=cmd.getOptionValue("c");
			}
			if(cmd.hasOption("bgmodel"))
			{
				GImprover.bgmodelFile=cmd.getOptionValue("bgmodel");
			}
			if(cmd.hasOption("prefix"))
			{
				GImprover.outputPrefix=cmd.getOptionValue("prefix");
			}
			if(cmd.hasOption("flank"))
			{
				GImprover.FlankLen=Integer.parseInt(cmd.getOptionValue("flank"));
			}

			if(cmd.hasOption("ratio"))
			{
				GImprover.sampling_ratio=Double.parseDouble( cmd.getOptionValue("ratio"));
			}
			if(cmd.hasOption("thresh"))
			{
				GImprover.entropyThresh=Double.parseDouble(cmd.getOptionValue("thresh"));
			}
			if(cmd.hasOption("maxlen"))
			{
				GImprover.max_gaplen=Integer.parseInt(cmd.getOptionValue("maxlen"));
			}
			if(cmd.hasOption("FDR"))
			{
				GImprover.FDR=Double.parseDouble(cmd.getOptionValue("FDR"));
			}
			if(cmd.hasOption("oops"))
			{
				GImprover.OOPS=true;
			}
			if(cmd.hasOption("oopg"))
			{
				GImprover.OOPG=true;
			}
			if(cmd.hasOption("rmbg"))
			{
				GImprover.removeBG=true;
			}

		} catch (ParseException e) {
			// TODO Auto-generated catch block
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp( "GapImprover", options );
			return;
		}
		
		GImprover.initialize();
		LinkedList<PWM> pwmlist=common.LoadPWMFromFile(inputPWM);
		Iterator<PWM> iter=pwmlist.iterator();
		//LinkedList<GapPWM> improvedPWMs=new LinkedList<GapPWM>();
		 File directory=new File(inputPWM);
		 String outdir=directory.getParent();
		File file = new File(outdir+"/GPimprover_sorted.dpwm"); 
		TreeMap<Double, PWM> sortedPWMs=new TreeMap<Double, PWM>();
		try {
			BufferedWriter writer= new BufferedWriter(new FileWriter(file));
			while(iter.hasNext())
			{
				PWM rawpwm=iter.next();
				System.out.println(rawpwm.Consensus(true));
				if(GImprover.removeBG)
					GImprover.SearchEngine.EnableBackground(GImprover.background);
				GapPWM gpwm=GImprover.fillDependency(rawpwm);
				gpwm.Name="GPimpover_"+rawpwm.Name;
				GImprover.SearchEngine.DisableBackground();
				GImprover.AUCtest(rawpwm);
				double score=GImprover.AUCtest(gpwm);
				gpwm.Score=score;
				sortedPWMs.put(score, gpwm);
				
			}
			for(Double key:sortedPWMs.descendingKeySet())
			{
				writer.write(sortedPWMs.get(key).toString());
			}
			writer.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		
	}

}