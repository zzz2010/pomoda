import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
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
import auc.AUCCalculator;
import auc.Confusion;


public class GapImprover {

	/**
	 * @param args
	 */
	public String outputPrefix="./";
	public String inputFasta;
	public String ctrlFasta="";
	public boolean OOPS=true; //only one dependence per sequence
	public boolean OOPG=false; //only one occurrence per sequence
	public boolean removeBG=true; //false:uniform BG assume
	public String bgmodelFile="";
	LinearEngine SearchEngine;
	public double sampling_ratio=1;
	public double FDR=0.01;
	public double entropyThresh=1;
	public int max_gaplen=8;
	
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
	
	public double KL_Divergence_empirical(List<String> sites,PWM motif)
	{
		
		double KLsum=0;
		HashMap<String, Integer> sitecount=new HashMap<String, Integer>();
		Iterator<String> iter=sites.iterator();
		int totalCount=0;
		while(iter.hasNext())
		{
			String temp=iter.next();
			temp=temp.substring((temp.length()-motif.core_motiflen)/2,(temp.length()-motif.core_motiflen)/2+motif.core_motiflen);
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
		SearchEngine.build_index(this.inputFasta,1000);
	
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
	public GapPWM fillDependency(PWM motif)
	{
		GapPWM gapPWM=null;
		String Consensus=motif.Consensus(true);
		ArrayList<Integer> gapstart=new ArrayList<Integer>(motif.core_motiflen/2);
		ArrayList<Integer> gapend=new ArrayList<Integer>(motif.core_motiflen/2);
		int FlankLen=2;//Math.min(Math.min(motif.head, motif.tail), 2);
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
			}
		}
		if(start!=-1)
		{
			gapstart.add(start);
			gapend.add(motif.core_motiflen+2*FlankLen);
		}
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
			 LinkedList<FastaLocation> falocs =SearchEngine.searchPattern(motif, Double.NEGATIVE_INFINITY);
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
		}
		try {
		//find the best dependency modeling in each gap region
		LinkedList<GapBGModelingThread> threadPool=new LinkedList<GapBGModelingThread>();
		if(sites.size()>0)
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
		}
		
		
	
		Iterator<GapBGModelingThread> iter3=threadPool.iterator();
		 start=-1;
		double minKL=Double.MAX_VALUE;
		GapBGModelingThread bestThread=null;
		HashMap<HashSet<Integer>,HashMap<String,Double>> Dmap=new HashMap<HashSet<Integer>,HashMap<String,Double>>();
		GapBGModelingThread[] Pos_BestThread=new GapBGModelingThread[motif.core_motiflen+2*FlankLen];
		
		while(iter3.hasNext())
		{
			GapBGModelingThread t1=iter3.next();	
				t1.join();
				if(t1.depend_Pos.size()==0)
					System.out.println(t1.toString());
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
		}
		if(!OOPG)
		{
		//fill in multi-dependency in the same gap region	
		for (int i = 0; i < Pos_BestThread.length; i++) {
			if(Pos_BestThread[i]!=null&&Pos_BestThread[i].depend_Pos.size()>1)
			{
				boolean maxCover=true;
				for(Integer ii:Pos_BestThread[i].depend_Pos)
				{
					if(Pos_BestThread[ii].hashCode()!=Pos_BestThread[i].hashCode())
					{
						maxCover=false;
						break;
					}
				}
				if(maxCover)
				{
				Dmap.put(Pos_BestThread[i].depend_Pos, Pos_BestThread[i].DprobMap);
				System.out.println(Pos_BestThread[i]);
				}
			}
		}
		}
	
		gapPWM=GapPWM.createGapPWM(motif.subPWM( motif.head,motif.head+motif.core_motiflen), Dmap,FlankLen);
		if(gapPWM.core_motiflen!=motif.core_motiflen&&sites.size()>0)
		{
			motif=new PWM(sites.toArray(new String[1]));
			gapPWM=GapPWM.createGapPWM(motif.subPWM(FlankLen,motif.columns()-FlankLen), Dmap,FlankLen);
			motif=motif.subPWM(gapPWM.head, gapPWM.head+gapPWM.core_motiflen);
		}
//		String testStr="agagaagaagaaagaaagaaagaagaaaggaagaaagaaagaaagaaagaaagaaagaaagaaagaaagaaagaaagaaagaaagaaaagaaagaaagaa";
//		testStr=common.getReverseCompletementString(testStr);
//		for (int i = 0; i < testStr.length()-gapPWM.core_motiflen; i++) {
//			String temp=testStr.substring(i, i+gapPWM.core_motiflen);
//			double score1=gapPWM.scoreWeightMatrix(temp);
//			double score2=motif.scoreWeightMatrix(temp);
//			System.out.println(score1+"\t"+score2);
//		}
		System.out.println(KL_Divergence_empirical(sites, motif) +"\t"+KL_Divergence_empirical(sites, gapPWM));
		
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
		options.addOption("ratio",true, "sampling ratio (default 1)");
		options.addOption("thresh",true, "minimum entropy threshold for considering a position as a gap(default 0.5)");
		options.addOption("maxlen",true,"maxmimum length of gap (default 8)");
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

		} catch (ParseException e) {
			// TODO Auto-generated catch block
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp( "GapImprover", options );
			return;
		}
		
		GImprover.initialize();
		LinkedList<PWM> pwmlist=common.LoadPWMFromFile(inputPWM);
		Iterator<PWM> iter=pwmlist.iterator();
		LinkedList<GapPWM> improvedPWMs=new LinkedList<GapPWM>();
		while(iter.hasNext())
		{
			PWM rawpwm=iter.next();
			System.out.println(rawpwm.Consensus(true));
			GapPWM gpwm=GImprover.fillDependency(rawpwm);
			
			GImprover.AUCtest(rawpwm);
			GImprover.AUCtest(gpwm);
			
		}
		
	}

}
