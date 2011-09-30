import java.awt.Color;
import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map.Entry;
import java.util.Random;
import java.util.TreeMap;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.lang.ArrayUtils;
import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dist.DistributionTools;
import org.biojava.bio.dist.SimpleDistribution;
import org.biojava.bio.dist.UniformDistribution;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Symbol;
import org.biojava.utils.ChangeVetoException;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.ApplicationFrame;
import org.jfree.ui.RectangleInsets;
import org.jfree.ui.RefineryUtilities;
import org.jfree.util.ArrayUtilities;

import umontreal.iro.lecuyer.probdist.HypergeometricDist;

import auc.AUCCalculator;
import auc.Confusion;


public class PWMevaluator {

	/**
	 * @param args
	 */
	public String outputPrefix="./";
	public String inputFasta;
	public String ctrlFasta="";
//	public boolean OOPS=true; //only one dependence per sequence
	public boolean removeBG=false; //false:uniform BG assume
	public String bgmodelFile="";
	public int BGFold=1;
	LinearEngine SearchEngine;
	LinearEngine BGSearchEngine;
	public double sampling_ratio=1;
	public double FDR=0.001;
	public double DivergenceThresh=0.24;
	
	public int resolution=10;
	
	public BGModel background=null;
	
	HashMap<String,XYSeries> ROCdata=new HashMap<String, XYSeries>();
	private boolean PBMflag=false;
	
	public PWMevaluator(Pomoda motiffinder)
	{
		//super("");

		SearchEngine=motiffinder.SearchEngine2;
		//sampling_ratio=motiffinder.sampling_ratio;
		FDR=motiffinder.FDR;
		background=motiffinder.background;
		removeBG=true;
		resolution=motiffinder.resolution;
		
	}
	
	public PWMevaluator()
	{//super("");
		
	}
	
	public static PWM shufflePWM(PWM motif)
	{
		int start=motif.head;
		int end=motif.columns()-motif.tail;
		Distribution[] dists=new Distribution[end-start];
		for (int i = start; i < end; i++) {
			FiniteAlphabet alpb=(FiniteAlphabet) motif.getColumn(i).getAlphabet();
			dists[i-start]=new SimpleDistribution(alpb)  ;
			Iterator<Symbol> iter= alpb.iterator();
			while(iter.hasNext())
			{
				Symbol sym=iter.next();
				try {
					dists[i-start].setWeight(sym, motif.getColumn(i).getWeight(sym));
				} catch (IllegalSymbolException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				} catch (ChangeVetoException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			
		}
		ArrayList<Distribution> randist=new ArrayList<Distribution>();
		for (int i = 0; i < dists.length; i++) {
			randist.add(dists[i]);
		}
		Collections.shuffle(randist);
		PWM sub=null;
		try {
			sub = new PWM(randist.toArray(dists));
			sub.core_motiflen=motif.core_motiflen-Math.max(0, start-motif.head)-Math.max(0, motif.columns()-end-motif.tail);
			sub.head=Math.max(0, motif.head-start);
			sub.tail=Math.max(0, motif.tail-(motif.columns()-end));
			sub.Name=motif.Name+"_rand";
			sub.Score=motif.Score;
			sub.strand_en=motif.strand_en;
			sub.pos_en=motif.pos_en;
			sub.peakrank_en=motif.peakrank_en;
			sub.pos_prior=(ArrayList<Double>)motif.pos_prior.clone();
			sub.peakrank_prior=(ArrayList<Double>)motif.peakrank_prior.clone();
			sub.strand_plus_prior=motif.strand_plus_prior;
		} catch (IllegalAlphabetException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		return sub;
	}
	
	public static TreeMap<String, Double> getKmerCount(List<String> seqs,List<Double> weightList,int len)
	{
		TreeMap<String, Double> KmerCount=new TreeMap<String, Double>();
		TreeMap<String, Integer> KmerCount2=new TreeMap<String, Integer>();
		int i=0;
		Iterator<String> seqsIter=seqs.iterator();
		while(seqsIter.hasNext())
		{
			double weight=1.0;
			if(weightList!=null)
				weight=weightList.get(i);
			String seqstr=seqsIter.next();
			for (int j = 0; j < seqstr.length()-len+1; j++) {
				String kmer=seqstr.substring(j,len);
				String rckmer=common.getReverseCompletementString(kmer);
				if(KmerCount.containsKey(kmer))
				{
					KmerCount.put(kmer, KmerCount.get(kmer)+weight);
					KmerCount2.put(kmer, KmerCount2.get(kmer)+1);
				}
				else
				{
					KmerCount.put(kmer, weight);
					KmerCount2.put(kmer, 1);
				}
				if(!rckmer.equalsIgnoreCase(kmer))
				{
					if(KmerCount.containsKey(rckmer))
					{
						KmerCount.put(rckmer, KmerCount.get(rckmer)+weight);
						KmerCount2.put(rckmer, KmerCount2.get(rckmer)+1);
					}
					else
					{
						KmerCount.put(rckmer, weight);
						KmerCount2.put(rckmer, 1);
					}
					
				}
			}
				
				i++;
		}
		
		//normalizing
		for(String key:KmerCount.keySet())
		{
			KmerCount.put(key, KmerCount.get(key)/KmerCount2.get(key));
		}
		
		return KmerCount;
		
	}
	
	public static int comparePositionList(List<FastaLocation> Falocs,String ansfile,int windowsize)
	{
		int maxseqlen=10000;
		ArrayList<Integer> anslist=new ArrayList<Integer>(1000);
		ArrayList<Integer> testlist=new ArrayList<Integer>(Falocs.size());
		try {
			BufferedReader readbuffer = new BufferedReader(new FileReader(ansfile));
			String strRead;
			while ((strRead=readbuffer.readLine())!=null){
				String splitarray[] = strRead.split("\t");
					int seqid=Integer.parseInt(splitarray[splitarray.length-2]);
					int pos=Integer.parseInt(splitarray[splitarray.length-1]);
					anslist.add(seqid*maxseqlen+pos);
				}
			Iterator<FastaLocation> iter=Falocs.iterator();
			while(iter.hasNext())
			{
				FastaLocation curr=iter.next();
				testlist.add(curr.getSeqId()*maxseqlen+curr.getSeqPos());
			}
			SortingThread st1=new SortingThread(anslist);
			SortingThread st2=new SortingThread(testlist);
			st1.run();
			
			st2.run();
			OverlappingThread ot=new OverlappingThread(st1.getResult(), st2.getResult(), windowsize);
			ot.run();
			readbuffer.close();
			return ot.getResult().size();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		
		return 0;
	}
	
	public void initialize()
	{
		common.initialize();
		SearchEngine=new LinearEngine(4);
		if(this.PBMflag)
		{
			SearchEngine.buildPBM_index(inputFasta, 100000,true);
		}
		else
			SearchEngine.build_index(this.inputFasta);
		if(background==null)
		{
			background=new BGModel();
			int bg_markov_order=0;
			if(!ctrlFasta.isEmpty())
			{
				BGSearchEngine=new LinearEngine(4);
				 background.BuildModel(this.ctrlFasta, 5); //5-order bg
				BGSearchEngine.build_index(this.ctrlFasta);
				removeBG=true;
			}
			
			else if(!bgmodelFile.isEmpty())
			{
				removeBG=true;
					background.LoadModel(bgmodelFile);
			}
			
			else
			{
				background=BGModel.CreateUniform();
			}

		}
		
	}
	
	public ArrayList<PWM> similarPWMs(PWM query, List<PWM> library,int K , double[] divergences)
	{
		ArrayList<PWM> similarList=new ArrayList<PWM>(K);
		Iterator<PWM> iter=library.iterator();
		TreeMap<Double, PWM> sortedPWMs=new TreeMap<Double, PWM>();
		while(iter.hasNext())
		{
			PWM p2=iter.next();
			double divergence=common.PWM_Divergence(query,p2 );
			sortedPWMs.put(divergence+(sortedPWMs.size()%K)*common.DoubleMinNormal, p2);
		}
		
		for(Double key:sortedPWMs.keySet())
		{
			divergences[similarList.size()]=key;
			if(key>DivergenceThresh)
				continue;
			similarList.add(sortedPWMs.get(key));
			
			if(similarList.size()==K)
				break;
		}
		
		return similarList;
	}
	
	
	public HashMap<String,Double> calc_mutliscore(PWM motif,LinearEngine bg_search)
	{
		///////////////compute SN,PPV,PC,ASP,CC////////////
		HashMap<String,Double> retscores=new HashMap<String,Double>();
		
		if(motif==null)
			return retscores;
		BGModel uniform_bg=new  BGModel();
		uniform_bg.BuildModel(new String[]{"ACGT"}, 1);
		double thresh=motif.getThresh(0.9999, FDR, uniform_bg,true);
		SearchThread.bestonly=true;
		 LinearEngine SearEngine=null;
		
			 SearEngine=this.SearchEngine;

			 LinearEngine BGSearch=null;
		if(bg_search==null)
		{
		  BGSearch=new LinearEngine(6);
         Iterator<String> iter2=SearchEngine.ForwardStrand.iterator();
         background.r.setSeed(0);
         while(iter2.hasNext())
         {
        	 int len=iter2.next().length();
        	 for (int i = 0; i < BGFold; i++) {		
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
	        	 BGSearch.accSeqLen.add( BGSearch.TotalLen);
        	 }
         }
		}
		else
		{
			BGSearch=bg_search;
		}

     	TreeMap<Double,Integer> Sorted_labels=new TreeMap<Double,Integer>();
     	double lamda=(double)SearchEngine.getSeqNum()/SearEngine.TotalLen/2;
        // SearchThread.recordSiteThreshold=Math.log((1-lamda)/lamda)+motif.core_motiflen*Math.log(0.25);
         motif.matchsite.clear();
        	 LinkedList<FastaLocation> falocs =SearchEngine.searchPattern(motif, thresh);
      double TP=falocs.size();
      double FN=SearchEngine.ForwardStrand.size()-TP;
        	 //bg sequences
        	 falocs =BGSearch.searchPattern(motif, thresh);
     double FP=falocs.size();
     double TN=BGSearch.ForwardStrand.size()-FP;
    double HGscore=HypergeometricDist.cdf((int)FP+1, BGSearch.ForwardStrand.size()+2, this.SearchEngine.ForwardStrand.size()+2, (int)TP+1) ;
    retscores.put("HG", HGscore);
    double SN=TP/(TP+FN);
    double SPC=TN/BGSearch.ForwardStrand.size();
    double PPV=TP/(TP+FP);
    double ASP=(SN+PPV)/2;
    double CC=(TP*TN-FP*FN)/Math.sqrt((TP+FN)*(TN+FP)*(TP+FP)*(TN*FN));
    retscores.put("SN", SN);
    retscores.put("SPC", SPC);
    retscores.put("PPV", PPV);
    retscores.put("ASP", ASP);
    retscores.put("CC", CC);
    return retscores;

	}
	
	public double calcAUC(PWM motif,LinearEngine bg_search)
	{
		if(motif==null)
			return Double.NEGATIVE_INFINITY;
		SearchThread.bestonly=true;
		 LinearEngine SearEngine=null;
		
			 SearEngine=this.SearchEngine;

			 LinearEngine BGSearch=null;
		if(bg_search==null)
		{
		  BGSearch=new LinearEngine(6);
         Iterator<String> iter2=SearchEngine.ForwardStrand.iterator();
         background.r.setSeed(0);
         while(iter2.hasNext())
         {
        	 int len=iter2.next().length();
        	 for (int i = 0; i < BGFold; i++) {

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
	        	 BGSearch.accSeqLen.add( BGSearch.TotalLen);
        	 }
         }
		}
		else
		{
			BGSearch=bg_search;
		}

     	TreeMap<Double,Integer> Sorted_labels=new TreeMap<Double,Integer>();
     	double lamda=(double)SearchEngine.getSeqNum()/SearEngine.TotalLen/2;
        SearchThread.recordSiteThreshold=Math.log((1-lamda)/lamda)+motif.core_motiflen*Math.log(0.25);
        motif.matchsite.clear();
        	 LinkedList<FastaLocation> falocs =SearchEngine.searchPattern(motif, Double.NEGATIVE_INFINITY);
        SearchThread.recordSiteThreshold=Double.POSITIVE_INFINITY;
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
        				 Sorted_labels.put(maxseq_score-seqcount*common.DoubleMinNormal, 1);
        				// System.err.println(maxseq_score);
        			 }
        				 lastseq=currloc.getSeqId();
        			 maxseq_score=currloc.Score;
        		 }
        		 if(maxseq_score<currloc.Score)
        		 {
        			 maxseq_score=currloc.Score;

        		 }
        	 }
        	 Sorted_labels.put(maxseq_score-seqcount*common.DoubleMinNormal, 1);
        	 
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
	    			 Sorted_labels.put(maxseq_score+seqcount*common.DoubleMinNormal, 0);
	    			
	    			 }
	    			 maxseq_score=currloc.Score;
	    			 lastseq=currloc.getSeqId();
	    		 }
	    		 if(maxseq_score<currloc.Score)
	    		 {
	    			 maxseq_score=currloc.Score;
	    		 }
	    		 
	    	 }
	       	Sorted_labels.put(maxseq_score+seqcount*common.DoubleMinNormal, 0);
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
	       	XYSeries series1 = new XYSeries(motif.Name);
	       	int poscount=0;
	       	int skip=labels.length/50;
	       	for (int i = 0; i < labels.length; i++) {
				if(labels[i]==1)
					poscount++;
				if(skip!=0&&i%skip==0)
				{
				series1.add((double)(i+1-poscount)/(labels.length-one), (double)(poscount)/one);
				//System.err.println(labels[i]+"\t"+scores[i]);
				}
			}
	       	ROCdata.put(motif.Name, series1);
	       	
         double AUCscore=AUCcalc.calculateAUCROC();
         SearchThread.bestonly=false;
         return AUCscore;

	}
	
	public double HyperGeometricScore(PWM motif, LinearEngine bg_search)
	{
		double HGscore=0;
		 LinearEngine SearEngine=null;
			
		 SearEngine=this.SearchEngine;
		double thresh=motif.getThresh(1, 0.0001, background,false);
		double lamda=(double)SearchEngine.getSeqNum()/SearEngine.TotalLen/2;
		 SearchThread.recordSiteThreshold=thresh;//Math.log((1-lamda)/lamda)+motif.core_motiflen*Math.log(0.25);
	        motif.matchsite.clear();
		LinkedList<FastaLocation> falocs= this.SearchEngine.searchPattern(motif, thresh);
		motif.Score=thresh;
	      SearchThread.recordSiteThreshold=Double.POSITIVE_INFINITY;
	      
			 LinearEngine BGSearch=null;
				if(bg_search==null)
				{
				  BGSearch=new LinearEngine(6);
		         Iterator<String> iter2=SearchEngine.ForwardStrand.iterator();
		         background.r.setSeed(0);
		         while(iter2.hasNext())
		         {
		        	 int len=iter2.next().length();
		        	 for (int i = 0; i < BGFold; i++) {

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
			        	 BGSearch.accSeqLen.add( BGSearch.TotalLen);
		        	 }
		         }
				}
				else
				{
					BGSearch=bg_search;
				}
		LinkedList<FastaLocation> falocs_bg= BGSearch.searchPattern(motif, thresh);
		
//		HGscore=HypergeometricDist.cdf(falocs_bg.size()+1, BGSearch.ForwardStrand.size()+2, this.SearchEngine.ForwardStrand.size()+2, falocs.size()+1) ;
		double p=(double)(falocs_bg.size()+1)/BGSearch.TotalLen;
		HGscore=(falocs.size()-this.SearchEngine.TotalLen*p)/Math.sqrt(this.SearchEngine.TotalLen*p*(1-p));
		return HGscore;
	}
	
	
	
	//compute rank correlation
	public static double CalcSpearmanRankCorr(PWM motif, TreeMap<String,Double> kmercount)
	{
		ArrayList<String> rankedKmer=new ArrayList<String>(kmercount.size());
		ArrayList<Double> PWMScore=new ArrayList<Double>(kmercount.size());
		HashSet<String> filterRC=new HashSet<String>(kmercount.size());
		for (Entry<String, Double> entry  : common.entriesSortedByValues(kmercount)) {
			String rckmer=common.getReverseCompletementString(entry.getKey());
			if(filterRC.contains(rckmer))
				continue;
			filterRC.add(entry.getKey());
		    rankedKmer.add(entry.getKey());
		    double score=motif.scoreWeightMatrix(entry.getKey());
		    double score2=motif.scoreWeightMatrix(rckmer);
		    if(score2>score)
		    	score=score2;
		    PWMScore.add(score);
		}
		
		
		 ArrayList<Double> sortedValues = new ArrayList<Double>(PWMScore);
		    Collections.sort(sortedValues);
		 
		    double[] ranks = new double[sortedValues.size()];
		    
		    double[] ranks_null = new double[sortedValues.size()];
		 
		    for (int i=0; i<PWMScore.size(); i++)
		    {
		        ranks[i]=((double)sortedValues.indexOf(PWMScore.get(i)));
		        ranks_null[i]=((double)i);
		    }
		 
		double spearman=common.getPearsonCorrelation(ranks, ranks_null);
		
		return spearman;
	}
	
	public double CalcCorrelation(PWM motif)
	{
		TreeMap<Double,Integer> Sorted_labels=new TreeMap<Double,Integer>();
        
   	 LinkedList<FastaLocation> falocs =SearchEngine.searchPattern(motif, Double.NEGATIVE_INFINITY);
   	 Iterator<FastaLocation> iter=falocs.iterator();
   	 int lastseq=-1;
   	 double seqcount=0;
   	 double maxseq_score=	Double.NEGATIVE_INFINITY;
   	 double[] scores=new double[SearchEngine.seqWeighting.size()];
   	 int seqid=0;
   	 while(iter.hasNext())
   	 {
   		 FastaLocation currloc=iter.next();
   		 if(lastseq!=currloc.getSeqId())
   		 {
   			 seqcount+=1;
   			 if(lastseq!=-1)
   			 {
   				 Sorted_labels.put(maxseq_score+seqcount*common.DoubleMinNormal, 1);
   				scores[seqid]=Math.exp(maxseq_score);//Math.exp
   				seqid++;
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
		scores[seqid]=Math.exp(maxseq_score); //Math.exp
			
		double corr=common.getPearsonCorrelation(scores,ArrayUtils.toPrimitive(SearchEngine.seqWeighting.toArray(new Double[1])));
   	
    return corr;
	}
	
	
	
	public double calcAUC(PWM motif,ArrayList<Double[]> DnaseLib, LinkedList<String> sequences)
	{
		 LinearEngine SearEngine=null;
		 boolean noPWMscore=true;
		 
		 if(sequences==null)
			 SearEngine=this.SearchEngine;
		 else
		 {
			 SearEngine=new LinearEngine(6);
			 SearEngine.ForwardStrand=sequences;
		 }
		 LinearEngine BGSearch=new LinearEngine(6);
         Iterator<String> iter2=SearchEngine.ForwardStrand.iterator();
         background.r.setSeed(0);
         while(iter2.hasNext())
         {
        	 int len=iter2.next().length();
        	 for (int i = 0; i < BGFold; i++) {
	
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
		        	 BGSearch.accSeqLen.add( BGSearch.TotalLen);
        	 }
         }

     	TreeMap<Double,Integer> Sorted_labels=new TreeMap<Double,Integer>();
         
        	 LinkedList<FastaLocation> falocs =SearchEngine.searchPattern(motif, Double.NEGATIVE_INFINITY);
        	 Iterator<FastaLocation> iter=falocs.iterator();
        	 int lastseq=-1;
        	 double seqcount=0;
        	 double maxseq_score=	Double.NEGATIVE_INFINITY;
        	 FastaLocation max_loc=null;
        	 int motiflen=motif.core_motiflen;
        	 while(iter.hasNext())
        	 {
        		 FastaLocation currloc=iter.next();
        		 if(lastseq!=currloc.getSeqId())
        		 {
        			 seqcount+=1;
        			 if(lastseq!=-1)
        			 {
        				 Double[] temp_dnase2=null;
        				 double logDnaseProb=0;
     					if(motif.Dnase_prob!=null)
     					{
     						int posbin=(int)(motif.pos_prior.size()*((max_loc.getSeqPos()+motiflen/2)%max_loc.getSeqLen()/(double)max_loc.getSeqLen()));
     						temp_dnase2=new Double[motif.Dnase_prob.size()];
     						Arrays.fill(temp_dnase2, new Double(0));
     						Double[] dnaseseq=DnaseLib.get(max_loc.getSeqId());
     						for (int i = 0; i <motif.Dnase_prob.size(); i++) {
     							if(max_loc.ReverseStrand)
     							{
     								temp_dnase2[i]=dnaseseq[posbin+motif.Dnase_prob.size()-i-1];
     							}
     							else
     							{
     								temp_dnase2[i]=dnaseseq[posbin+i];
     							}
     						}
     						
     						logDnaseProb=motif.calcLogDnaseProb(temp_dnase2,0);
     					}
     					if(noPWMscore)
     						maxseq_score=0;
        				 Sorted_labels.put(maxseq_score+logDnaseProb-seqcount*common.DoubleMinNormal, 1);
        			 }
        				 lastseq=currloc.getSeqId();
        			 maxseq_score=currloc.Score;
        			 max_loc=currloc;
        		 }
        		 if(maxseq_score<currloc.Score)
        		 {
        			 maxseq_score=currloc.Score;
        			 max_loc=currloc;

        		 }
        	 }
				if(noPWMscore)
						maxseq_score=0;
        	 Sorted_labels.put(maxseq_score-seqcount*common.DoubleMinNormal, 1);
        	 
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
	    				 Double[] temp_dnase2=null;
        				 double logDnaseProb=0;
     					if(motif.Dnase_prob!=null)
     					{
     						int posbin=(int)(motif.pos_prior.size()*((max_loc.getSeqPos()+motiflen/2)%max_loc.getSeqLen()/(double)max_loc.getSeqLen()));
     						temp_dnase2=new Double[motif.Dnase_prob.size()];
     						Arrays.fill(temp_dnase2, new Double(0));
     						Double[] dnaseseq=DnaseLib.get(max_loc.getSeqId());
     						Random r=new Random(common.randomseed);
     						for (int i = 0; i <motif.Dnase_prob.size(); i++) {
     							{
     								temp_dnase2[i]= dnaseseq[r.nextInt(dnaseseq.length)];//posbin+i
     							}
     						}
     						
     						logDnaseProb=motif.calcLogDnaseProb(temp_dnase2,0);
     					}
     					if(noPWMscore)
     						maxseq_score=0;
	    			 Sorted_labels.put(maxseq_score+logDnaseProb+seqcount*common.DoubleMinNormal, 0);
	    			
	    			 }
	    			 maxseq_score=currloc.Score;
	    			 max_loc=currloc;
	    			 lastseq=currloc.getSeqId();
	    		 }
	    		 if(maxseq_score<currloc.Score)
	    		 {
	    			 maxseq_score=currloc.Score;
	    			 max_loc=currloc;
	    		 }
	    		 
	    	 }
				if(noPWMscore)
						maxseq_score=0;
	       	Sorted_labels.put(maxseq_score+seqcount*common.DoubleMinNormal, 0);
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
	       	XYSeries series1 = new XYSeries(motif.Name);
	       	int poscount=0;
	       	int skip=labels.length/50;
	       	for (int i = 0; i < labels.length; i++) {
				if(labels[i]==1)
					poscount++;
				if(i%skip==0)
				series1.add((double)(i+1-poscount)/(labels.length-one), (double)(poscount)/one);
			}
	       	ROCdata.put(motif.Name, series1);
	       	
         double AUCscore=AUCcalc.calculateAUCROC();
         
         return AUCscore;

	}
	
	public  void DrawROC(String pngfile)
	{
		  XYSeriesCollection dataset = new XYSeriesCollection();
		 for(String name: ROCdata.keySet())
		 {
			 dataset.addSeries(ROCdata.get(name));
		 }
		 JFreeChart chart = ChartFactory.createXYLineChart(
	                "ROC curve", // chart title
	                "False Positive Rate", // x axis label
	                "True Positive Rate", // y axis label
	                dataset, // data
	                PlotOrientation.VERTICAL,
	                true, // include legend
	                true, // tooltips
	                false // urls
	                );
	// NOW DO SOME OPTIONAL CUSTOMISATION OF THE CHART...
	        chart.setBackgroundPaint(Color.white);
	// get a reference to the plot for further customisation...
	        XYPlot plot = (XYPlot) chart.getPlot();
	        plot.setBackgroundPaint(Color.lightGray);
	        plot.setAxisOffset(new RectangleInsets(5.0, 5.0, 5.0, 5.0));
	        plot.setDomainGridlinePaint(Color.white);
	        plot.setRangeGridlinePaint(Color.white);
	        XYLineAndShapeRenderer renderer = (XYLineAndShapeRenderer) plot.getRenderer();
	        renderer.setShapesVisible(true);
	        renderer.setShapesFilled(true);
	// change the auto tick unit selection to integer units only...
	        NumberAxis rangeAxis = (NumberAxis) plot.getRangeAxis();
	        rangeAxis.setStandardTickUnits(NumberAxis.createStandardTickUnits());
	        
	        ChartPanel chartPanel = new ChartPanel(chart);
	        chartPanel.setPreferredSize(new java.awt.Dimension(800, 600));
	      //  setContentPane(chartPanel);
	        try {
				ChartUtilities.saveChartAsPNG(new File(pngfile), chart, 800, 600);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}

	}
	
	
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		Options options = new Options();
		options.addOption("i", true, "input fasta file (or PBM format file with -pbm)");
		options.addOption("pwm", true, "input PWM file");
		options.addOption("N", true, "the number of PWM want to evaluate");
		options.addOption("c", true, "control fasta file");
		options.addOption("dpwm", false,"compare dpwm and its pwm version");
		options.addOption("convert", false, "convert input PWM file to the transfac format");
		options.addOption("roc", false, "compute AUC and draw ROC curve for the given pwm file");
		options.addOption("pbm", false, "input file is PBM format, and will compute signal correlation");
		options.addOption("multiscore", false, "compute SN,PPV,PC,ASP,CC for the given pwm file");
		options.addOption("match", true, "find similar motifs in known PWM library (path to the library, e.g., jaspar.pwm)");
		options.addOption("bgmodel", true, "background model file");
		options.addOption("bgfold", true, "if background model is set, the fold change between negative sequence and positive sequences (default 1, Integer)");
		options.addOption("markov", true, "use markov model of the control sequences rather than directly control sequences");
		options.addOption("prefix", true, "output directory");
		options.addOption("ratio",true, "sampling ratio (default 1)");
		options.addOption("genrand",false, "shuffle the colums of input PWMs and evaluate");
		options.addOption("thresh",true, "minimum PWM divergence threshold for considering a match known motif(default 0.24)");
		options.addOption("FDR",true,"fasle positive rate");
		String inputPWM;
		CommandLineParser parser = new GnuParser();
		PWMevaluator evaluator=new PWMevaluator();
		boolean rocflag=false;
		boolean dpwmflag=false;
		boolean genrand=false;
		boolean multiscore_flag=false;
		boolean convertflag=false;
		
		LinkedList<PWM> PWMLibrary=null;
		int topN=1000000;
		try {
			CommandLine cmd = parser.parse( options, args);
			if(cmd.hasOption("i"))
			{
				evaluator.inputFasta=cmd.getOptionValue("i");
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
				evaluator.ctrlFasta=cmd.getOptionValue("c");
				if(cmd.hasOption("markov"))
				{
					common.initialize();
					int bgorder=Integer.parseInt(cmd.getOptionValue("markov"))+1;
					evaluator.removeBG=true;
					evaluator.background=new BGModel();
					evaluator.background.BuildModel(evaluator.ctrlFasta, bgorder);
				}
			}
			if(cmd.hasOption("N"))
			{
				topN=Integer.parseInt(cmd.getOptionValue("N"));
			}
			if(cmd.hasOption("bgmodel"))
			{
				evaluator.bgmodelFile=cmd.getOptionValue("bgmodel");
			}
			if(cmd.hasOption("bgfold"))
			{
				evaluator.BGFold=Integer.parseInt(cmd.getOptionValue("bgfold"));
			}
			if(cmd.hasOption("match"))
			{
				PWMLibrary=common.LoadPWMFromFile(cmd.getOptionValue("match"));
			}
			if(cmd.hasOption("roc"))
			{
				rocflag=true;
			}
			if(cmd.hasOption("dpwm"))
			{
				dpwmflag=true;
			}
			if(cmd.hasOption("pbm"))
			{
				evaluator.PBMflag=true;
			}
			if(cmd.hasOption("multiscore"))
			{
				multiscore_flag=true;
			}
			if(cmd.hasOption("convert"))
			{
				convertflag =true;
			}
			if(cmd.hasOption("prefix"))
			{
				evaluator.outputPrefix=cmd.getOptionValue("prefix");
			}

			if(cmd.hasOption("ratio"))
			{
				evaluator.sampling_ratio=Double.parseDouble( cmd.getOptionValue("ratio"));
			}
			if(cmd.hasOption("thresh"))
			{
				evaluator.DivergenceThresh=Double.parseDouble(cmd.getOptionValue("thresh"));
			}
			if(cmd.hasOption("FDR"))
			{
				evaluator.FDR=Double.parseDouble(cmd.getOptionValue("FDR"));
			}
			if(cmd.hasOption("genrand"))
			{
				genrand=true;
			}

		} catch (ParseException e) {
			// TODO Auto-generated catch block
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp( "PWMevaluator", options );
			return;
		} 
		
		
		 
		try {
			
			if(genrand)
			{
				File randfile = new File(inputPWM+"_rand.pwm");
				BufferedWriter randwriter = new BufferedWriter(new FileWriter(randfile));
				List<PWM> pwmlist=common.LoadPWMFromFile(inputPWM);
				Iterator<PWM> iter=pwmlist.iterator();
				while(iter.hasNext())
				{
					PWM p1=iter.next();
					randwriter.write(shufflePWM(p1).toString()+"\n");
					
				}
				randwriter.close();
				inputPWM=inputPWM+"_rand.pwm";
			}
			File file = new File(inputPWM+"_eval.txt");
			BufferedWriter writer = new BufferedWriter(new FileWriter(file));
		
		if(rocflag)
		{
		evaluator.initialize();
		
		List<PWM> pwmlist=common.LoadPWMFromFile(inputPWM);
		if(pwmlist.size()>topN)
			pwmlist= pwmlist.subList(0, topN);
		Iterator<PWM> iter=pwmlist.iterator();
		writer.write("AUC Result:\n");
		TreeMap<Double,PWM> sortedPWMs=new TreeMap<Double,PWM>();
		while(iter.hasNext())
		{
			PWM p1=iter.next();
			if(dpwmflag&&p1 instanceof GapPWM )
			{
				pwmlist.add(p1.Clone());
			}
			double auc=0;
			if(evaluator.SearchEngine.TotalLen/evaluator.SearchEngine.getSeqNum()<500)
				auc=evaluator.calcAUC(p1, evaluator.BGSearchEngine);
			else
				auc=evaluator.HyperGeometricScore(p1, evaluator.BGSearchEngine);
			p1.Score=auc;
			writer.write(p1.Name+"\t"+auc+"\n");
			sortedPWMs.put(auc, p1);
		}
		
		if(convertflag)
		{
			File file2 = new File(inputPWM+"_sorted.pwm"); 
			BufferedWriter writer2= new BufferedWriter(new FileWriter(file2));
			for(Double key:sortedPWMs.descendingKeySet())
			{
				PWM motif=sortedPWMs.get(key);
				writer2.write(motif.toString());
			}
			writer2.close();
			convertflag=false;
		}
		evaluator.DrawROC(inputPWM+"_roc.png");
		//evaluator.pack();
	     //   RefineryUtilities.centerFrameOnScreen(evaluator);
	      //  evaluator.setVisible(true);
		}

		
		
		if(multiscore_flag)
		{
			evaluator.initialize();
			
			List<PWM> pwmlist=common.LoadPWMFromFile(inputPWM);
			if(pwmlist.size()>topN)
				pwmlist= pwmlist.subList(0, topN);
			Iterator<PWM> iter=pwmlist.iterator();
			writer.write("MultiScore Result:\n");
			while(iter.hasNext())
			{
				PWM p1=iter.next();
				if(dpwmflag&&p1 instanceof GapPWM )
				{
					pwmlist.add(p1.Clone());
				}
				HashMap<String,Double> multiscore=evaluator.calc_mutliscore(p1, evaluator.BGSearchEngine);
				writer.write(p1.Name);
				String scoreresult="";
				for(String scorename:multiscore.keySet())
				{
					scoreresult+="\t"+scorename+"|";
					scoreresult+=multiscore.get(scorename).toString();
				}
				writer.write(scoreresult+"\n");
			}
			
		}
		
		if(evaluator.PBMflag)
		{
			evaluator.initialize();
			
			List<PWM> pwmlist=common.LoadPWMFromFile(inputPWM);
			if(pwmlist.size()>topN)
				pwmlist= pwmlist.subList(0, topN);
			Iterator<PWM> iter=pwmlist.iterator();
			writer.write("Correlation Result:\n");
			while(iter.hasNext())
			{
				PWM p1=iter.next();
				if(dpwmflag&&p1 instanceof GapPWM )
				{
					pwmlist.add(p1.Clone());
				}
				double corr=0;
					corr=evaluator.CalcCorrelation(p1);
				p1.Score=corr;
				writer.write(p1.Name+"\t"+corr+"\n");
				System.out.println("Pearson Correlation PBM:"+corr);
			}
			//compute the rank correlation
			iter=pwmlist.iterator();
			while(iter.hasNext())
			{
				PWM p1=iter.next();
				TreeMap<String,Double> Kmercount=getKmerCount(evaluator.SearchEngine.ForwardStrand, evaluator.SearchEngine.seqWeighting, p1.columns());
				double corr=0;
				corr=CalcSpearmanRankCorr(p1, Kmercount);
				p1.Score=corr;
				writer.write(p1.Name+"\t"+corr+"\n");
				System.out.println("Spearman Correlation PBM:"+corr);
				
			}
		}
		
		
		if(convertflag)
		{
			File file2 = new File(inputPWM+"_sorted.pwm"); 
			BufferedWriter writer2= new BufferedWriter(new FileWriter(file2));
			LinkedList<PWM> pwmlist=common.LoadPWMFromFile(inputPWM);
			for(PWM p1:pwmlist)
			{
				writer2.write(p1.toString());
			}
			writer2.close();
			convertflag=false;
		}
	    if(PWMLibrary!=null)
	    {
	    	writer.write("Similar Known Motifs:\n");
	    	LinkedList<PWM> pwmlist=common.LoadPWMFromFile(inputPWM);
			Iterator<PWM> iter=pwmlist.iterator();
			
			while(iter.hasNext())
			{
				PWM p1=iter.next();
				int K=200;
				double [] pwmdivergences=new double[K];
			    ArrayList<PWM> similars=evaluator.similarPWMs(p1, PWMLibrary, K,pwmdivergences);
			    writer.write(p1.Name);
			    for (int i = 0; i < similars.size(); i++) {
			    	writer.write("\t"+similars.get(i).Name+"|"+pwmdivergences[i]);
				}
				writer.write("\n");
			}
	    }
	    writer.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

}
