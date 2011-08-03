/**
 * @author zhizhuo zhang
 * zzz2010@gmail.com
 */
import java.awt.Color;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Random;
import java.util.Set;
import java.util.TreeMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.lang.ArrayUtils;
import org.biojava.bio.BioException;
import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dist.DistributionTools;
import org.biojava.bio.dp.SimpleWeightMatrix;
import org.biojava.bio.dp.WeightMatrix;
import org.biojava.bio.dp.WeightMatrixAnnotator;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.Feature;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.symbol.Alignment;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.SimpleAlignment;
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
import org.jfree.ui.RectangleInsets;
import org.pr.clustering.hierarchical.LinkageCriterion;

import umontreal.iro.lecuyer.probdist.BinomialDist;
import umontreal.iro.lecuyer.probdist.ChiSquareDistQuick;
import umontreal.iro.lecuyer.probdist.NegativeBinomialDist;
import umontreal.iro.lecuyer.probdist.NormalDist;
import umontreal.iro.lecuyer.probdistmulti.DirichletDist;
import umontreal.iro.lecuyer.util.Num;

import cern.jet.random.Binomial;
import cern.jet.random.engine.RandomEngine;
import edu.stanford.rsl.jpop.AdditiveFunctionAssembler;
import edu.stanford.rsl.jpop.FunctionOptimizer;
import edu.stanford.rsl.jpop.SimpleFunctionController;
import edu.stanford.rsl.jpop.FunctionOptimizer.OptimizationMode;





public class Pomoda {

	public String outputPrefix="./";
	public String inputFasta;
	public int max_iterNum=50;
	public int MIN_SAMPLENUM=1000;
	public int MAX_P=2;
	public String ctrlFasta="";
	public int BGorder=2;
	public String bgmodelFile="";
	public int seedlen=5;
	public boolean debug=false;
	public boolean OOPS=false;
	public int resolution=20;
	public int min_motiflen=7;
	public int ending_windowsize=600;
	public double FDR=0.01;
	public int max_motiflen=55;
	public int max_threadNum=6;
	public int num_motif=5;
	public double sampling_ratio=0.01;
	public double min_support_ratio=0.1;
	public boolean maskflag=true;
	
	public LinearEngine SearchEngine2;
	LinearEngine BGSearch=null;
	public BGModel background;
	
	public String linkage="WPGMA";
	public ArrayList<Double> pos_prior;
	public NegativeBinomialDist dnaseBG=null;
	public ArrayList<Double[]> DnaseLib=null;
	public int DnaseWindow=1;
	
	public void initialize()
	{
		//build hash index
		common.initialize();
		
//		SearchEngine=new HashEngine(seedlen);
//		SearchEngine.build_index(this.inputFasta);
		
		
		

		SearchEngine2=new LinearEngine(max_threadNum);
		SearchEngine2.build_index(this.inputFasta);
		if(SearchEngine2.num_thread>SearchEngine2.TotalLen/400000)
			SearchEngine2.num_thread=SearchEngine2.TotalLen/400000+1;
		SearchEngine2.build_KmerHitList(seedlen);
		

		
//		if(SearchEngine_Test())
//			System.out.println("SearchEngine_Test : pass");
			
//		double[][] sim=SearchEngine2.seq_similarity();
//		double maxsim=0;
//		int maxi=0,maxj=0;
//		for (int i = 0; i < sim.length; i++) {
//			for (int j = i; j < sim.length; j++) {
//				if(maxsim<sim[i][j])
//				{
//					maxi=i;
//					maxj=j;
//					maxsim=sim[i][j];
//				}
//			}
//		}
		
		if(ctrlFasta.equalsIgnoreCase(inputFasta))
			ctrlFasta="";
		background=new BGModel();
		File file=null;
		int bg_markov_order=BGorder-1;
		if(ctrlFasta.isEmpty())
		{
			bg_markov_order=1;
			file= new File(inputFasta+".bg");
		}
		else
		{
			file= new File(ctrlFasta+".bg");
			
		}
		
		if(!bgmodelFile.isEmpty())
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
		     background.SaveModel(inputFasta+".bg");
			}
			else
			{
				BGSearch=new LinearEngine(max_threadNum);
				BGSearch.build_index(ctrlFasta, SearchEngine2.ForwardStrand);
				bg_markov_order=seedlen-1;
				background.EnablePWMBG=true;
				background.kmerlen=seedlen;
				background.flanklen=(this.max_motiflen-seedlen)/2;
				background.BuildModel(BGSearch.ForwardStrand.toArray(new String[1]), bg_markov_order+1);
			    // background.BuildModel(ctrlFasta, bg_markov_order+1); //3-order bg
			     background.SaveModel(ctrlFasta+".bg");
			}
				
		}
//		if(BGModel_Test())
//			System.out.println("BGModel_Test : pass");
		pos_prior=new ArrayList<Double>(SearchEngine2.getTotalLength()/SearchEngine2.getSeqNum()/this.resolution);
	//	SearchEngine2.EnableBackground(background);
		
		if(DnaseLib!=null)
		{
			DnaseWindow=Integer.MAX_VALUE;
			Iterator<String> iter=SearchEngine2.ForwardStrand.iterator();
			int seqid=0;
			while(iter.hasNext())
			{
				int seqlen=iter.next().length();
				Double[] dnaseArr=DnaseLib.get(seqid);
				DnaseWindow=Math.min(DnaseWindow, (dnaseArr.length-seqlen)/2);
				
				Double[] m_dnaseArr=new Double[dnaseArr.length/resolution];
				
				
				Arrays.fill(m_dnaseArr, new Double(0));
				for (int i = 0; i < dnaseArr.length; i++) {
					if(i/resolution<m_dnaseArr.length)

						m_dnaseArr[i/resolution]+=dnaseArr[i];

				}
				DnaseLib.set(seqid, m_dnaseArr);
				seqid++;
				if(DnaseWindow<=0)
				{
					System.out.println("Dnase Data Error: Line"+seqid);
					DnaseLib=null;
					break;
				}
			}
			
			DnaseWindow/=resolution;
//			double mean=0;
//			double variance=0;
//			double count=0;
//			ArrayList<Integer> stat=new ArrayList<Integer>(SearchEngine2.TotalLen);
//			for (int i = 0; i < DnaseLib.size(); i++) {
//				  Double[] arr=DnaseLib.get(i);
//				  for (int j = 0; j < arr.length-2*DnaseWindow; j++) {
//					double sum=0;
//					for (int k = 0; k < 2*DnaseWindow; k++) {
//						sum+=arr[k+j];
//					}
//					
//					mean+=sum;
//					stat.add((int)sum);
//					count++;
//					variance+=sum*sum;
//				}
//			}
//			mean/=count;
//			variance/=count;
//			variance-=mean*mean;
//			
//			double p,r;
//			p=mean/variance;
//			r=p*mean/(1-p);
			//double[] paras=NegativeBinomialDist.getMLE(ArrayUtils.toPrimitive(stat.toArray(new Integer[1])),stat.size());
			dnaseBG=new NegativeBinomialDist(1,0.5);//(r, p);
		}
		
		PWM.bg_prob=new double[]{Math.exp(background.Get_LOGPROB("A")),Math.exp(background.Get_LOGPROB("C")),Math.exp(background.Get_LOGPROB("G")),Math.exp(background.Get_LOGPROB("T"))};
		
////		
//		if(GAP_Test())
//			System.out.println("PWM_Test : pass");
//	    System.exit(1);
		
	}
	

	private boolean BGModel_Test()
	{
		boolean pass=true;
		double score1=0;
		String testpattern="ACGTAC";
		score1=background.Get_LOGPROB(testpattern);
		System.out.println(score1);
		background.SaveModel("bg.obj");
		BGModel bgtest=new BGModel();
		bgtest.LoadModel("bg.obj");
		double score2=bgtest.Get_LOGPROB(testpattern);
		if(score2!=score1)
		{
			pass=false;
			System.out.println(score2);
		}
		
		
		return pass;
	}
	
	
	
	
	

	private boolean GAP_Test()
	{
		boolean pass=true;
		try {
			
			LinkedList<PWM> AR=new LinkedList<PWM>();
			
			AR.add(new PWM(new String[]{"NNNNNNACANNTGTNNNNN"}));
			AR.add(new PWM(new String[]{"NNNGNACANNTGTNCNNN"}));
			AR.add(new PWM(new String[]{"NNNNNNACANNNTGTNNNNN"}));
			AR.add(new PWM(new String[]{"NNNGNACANNNTGTNCNNN"}));
			AR.add(new PWM(new String[]{"NNNNNNACANNNNTGTNNNNN"}));
			AR.add(new PWM(new String[]{"NNNGNACANNNNTGTNCNNN"}));
			AR.add(new PWM(new String[]{"NNNNNNACANNNNNTGTNNNNN"}));
			AR.add(new PWM(new String[]{"NNNGNACANNNNNTGTNCNNN"}));
			AR.add(new PWM(new String[]{"NNNNNNACANNNNNNTGTNNNNN"}));
			AR.add(new PWM(new String[]{"NNNGNACANNNNNNTGTNCNNN"}));

			
			GapImprover gimprover=new GapImprover(this);
			GapPWM gPWM=gimprover.fillDependency(AR.get(3));
			double gthresh=gPWM.getThresh(sampling_ratio, FDR, background);
			LinkedList<FastaLocation> falocs=SearchEngine2.searchPattern(gPWM, gthresh);
			int dd=falocs.size();
			
		} catch (IllegalAlphabetException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IllegalSymbolException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		return pass;
	}
	
	
	public HashSet<String> MismatchSet(String pattern, int mismatch)
	{
		HashSet<String> scannedSet=new HashSet<String>();
		char[] ACGT=new char[]{'A','C','G','T'};
		scannedSet.add(pattern);
		scannedSet.add(common.getReverseCompletementString(pattern));
		for (int i = 0; i < mismatch; i++) {
			HashSet<String> tempset=new HashSet<String>(); 
			for(String patt:scannedSet)
			for (int j = 0; j < patt.length(); j++) {
				StringBuilder sb=new StringBuilder(patt);
				for (int k = 0; k < 4; k++) {
					sb.setCharAt(j, ACGT[k]);
					tempset.add(sb.toString());
				}
				sb.setCharAt(j, patt.charAt(j));
			}
			
			scannedSet.addAll(tempset);
		}
		return scannedSet;
	}
	private boolean MatrixAnnotatorTest()
	{
		

	    Map map = new HashMap();
	    try {
			map.put("seq0", DNATools.createDNA("aggag"));
		    map.put("seq1", DNATools.createDNA("aggaa"));
		    map.put("seq2", DNATools.createDNA("aggag"));
		    map.put("seq3", DNATools.createDNA("aagag"));
		    Alignment align = new SimpleAlignment(map);
		 
		    //make a Distribution[] of the motif
		    Distribution[] dists =
		        DistributionTools.distOverAlignment(align, false, 0.01);
		 
		    //make a Weight Matrix
		    WeightMatrix matrix = new SimpleWeightMatrix(dists);
		 
		    //the sequence to score against
		    Sequence seq = DNATools.createDNASequence("aaagcctaggaagaggagctgat","seq");
		 
		    //annotate the sequence with the weight matrix using a low threshold (0.1)
		    WeightMatrixAnnotator wma = new WeightMatrixAnnotator(matrix, 0.1);
		    seq = wma.annotate(seq);
		 
		    //output match information
		    for (Iterator it = seq.features(); it.hasNext(); ) {
		      Feature f = (Feature)it.next();
		      Location loc = f.getLocation();
		      System.out.println("Match at " + loc.getMin()+"-"+loc.getMax());
		      System.out.println("\tscore : "+f.getAnnotation().getProperty("score"));
		    }
		} catch (IllegalSymbolException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IllegalAlphabetException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (ChangeVetoException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (BioException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	    return true;
	}
	

	
	public void Masking(PWM motif)
	{
		int masklen=0;
		double log_thresh=motif.getThresh(1.0, 0.0001, background);
		String consensus_core=motif.Consensus(true);
		
		
		//double logN025=log025*(motif.head+motif.tail);
		int motiflen=consensus_core.length();
		LinkedList<String> MatchSite=new LinkedList<String>();
		
		//update the loglik matrix
		SearchEngine2.DisableBackground();
		LinkedList<FastaLocation> Falocs=SearchEngine2.searchPattern(motif, log_thresh);
		Iterator<FastaLocation> iter2=Falocs.iterator();
		int count=0;
		int lastseq=-1;
	
		String X_Str="";
		for (int i = 0; i < motiflen; i++) {
			X_Str+="N";
		}
		while(iter2.hasNext())
		{
			//masklen+=seedlen;
			FastaLocation currloc=iter2.next();
			String Site1=SearchEngine2.getSite(currloc.getSeqId(), currloc.getSeqPos(), motiflen);
			
			String rep=SearchEngine2.ForwardStrand.get(currloc.getSeqId()).replace(Site1, X_Str);
			SearchEngine2.ForwardStrand.set(currloc.getSeqId(), rep);

			
		}

		SearchEngine2.TotalLen-=masklen;
		
	}
	
	public ArrayList<PWM>	 getSeedMotifs3(int max_num_Seeds) {
		ArrayList<PWM>	SeedMotifs=new	ArrayList<PWM>(max_num_Seeds);
		int maxkmerlen=(int) (Math.log(2*SearchEngine2.TotalLen)/Math.log(4));
		if(maxkmerlen<seedlen)
			return SeedMotifs;
		HashMap<String,Double> PatternLib=new HashMap<String, Double>((int) Math.pow(maxkmerlen, 4));
		Iterator<String> iter=SearchEngine2.ForwardStrand.iterator();
		while(iter.hasNext())
		{
			String seq=iter.next();
			for (int i = 0; i < seq.length()-maxkmerlen; i++) {
				String patt=seq.substring(i, i+maxkmerlen).toUpperCase();
				String revpatt=common.getReverseCompletementString(patt);
				if(PatternLib.containsKey(patt))
				{
					PatternLib.put(patt, PatternLib.get(patt)+1);
				}
				else if(PatternLib.containsKey(revpatt))
				{
					patt=revpatt;
					PatternLib.put(patt, PatternLib.get(patt)+1);
				}
				else
				{
					PatternLib.put(patt, 1.0);
				}
			}
		}
		
		TreeMap<Double,String> sortPattLib=new TreeMap<Double, String>();
		double minscore=0;
		for(Entry<String,Double> pair:PatternLib.entrySet())
		{
			String pattern=pair.getKey();
			double logprob_bg=Math.log(Math.exp(background.Get_LOGPROB(pattern))+Math.exp(background.Get_LOGPROB(common.getReverseCompletementString(pattern)) )) ;
			double oddratio=PatternLib.get(pattern)/Math.exp(logprob_bg);
//			if(Double.isNaN(oddratio)||oddratio==0)
//				break;
			sortPattLib.put(oddratio,pattern);
			if(max_num_Seeds<sortPattLib.size())
			{
				minscore=sortPattLib.firstKey();
			}
		}

			//just take top ones
		    for (Map.Entry<Double,String> pair : sortPattLib.descendingMap().entrySet()) {
		    	if(SeedMotifs.size()==max_num_Seeds)
		    		break;
			    String toppattern=pair.getValue();
			    System.out.println(toppattern+"\t"+pair.getKey());
			    StringBuffer sb=new StringBuffer(toppattern);
			    for (int j = 0; j < (max_motiflen-this.min_motiflen)/2; j++) {
			    	sb.insert(0, '-');
			    	sb.append('-');					
				}
			    
			    String[] seqs={sb.toString()};
			    try {
					PWM a=new PWM(seqs);
					a.Score=pair.getKey();
					SeedMotifs.add(a);
					//System.out.println(a.Score);
				} catch (IllegalAlphabetException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				} catch (IllegalSymbolException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
		    }
		return SeedMotifs;
	}
	public ArrayList<PWM>	 getSeedMotifs2() {
		int max_num_Seeds=num_motif*num_motif;
		ArrayList<PWM>	SeedMotifs=new	ArrayList<PWM>(max_num_Seeds);
		TreeMap<Double,String> overrepKmer=new TreeMap<Double, String>();
		TreeMap<String,Double> PatternLib=new TreeMap<String, Double>();
		for(String patt:SearchEngine2.KmerHitList.keySet())
		{
			double orratio=SearchEngine2.KmerHitList.get(patt).size()/Math.exp(background.Get_LOGPROB(patt))/SearchEngine2.TotalLen;
			if(orratio>1)
			overrepKmer.put(orratio, patt);
		}
		for(Double orratio:overrepKmer.descendingKeySet())
		{
			String patt=overrepKmer.get(orratio);
			HashMap<Integer,Double> pattScore=new HashMap<Integer, Double>();
			for(Entry<Double, String> pair: overrepKmer.entrySet())
			{
				String patt2=pair.getValue();
				
				int mismatch=0;
				int pattCode=0;
				for (int i = 0; i < patt2.length(); i++) {
					if(patt2.charAt(i)!=patt.charAt(i))
					{
						pattCode+=i*Math.pow(patt.length(),mismatch);
						mismatch++;
						if(mismatch>patt.length()-seedlen)
							break;
					}
				}
				if(mismatch<=patt.length()-seedlen)
				{
					if(pattScore.containsKey(pattCode))
					{
						pattScore.put(pattCode, pattScore.get(pattCode)+orratio);
					}
					else
					{
						pattScore.put(pattCode, orratio);
					}
				}
			}
			//find the best gap pattern
			int bestPattCode=0;
			double bestscore=0;
			for(Entry<Integer,Double> pair:pattScore.entrySet())
			{
				if(pair.getValue()>bestscore)
				{
					bestPattCode=pair.getKey();
					bestscore=pair.getValue();
				}
			}
			StringBuilder sb=new StringBuilder(patt);
			for (int i = 0; i < patt.length()-seedlen; i++) {
				int pos=bestPattCode%patt.length();
				sb.setCharAt(pos,'N');
				bestPattCode-=pos;
				bestPattCode/=patt.length();
			}
			PatternLib.put(sb.toString(), bestscore+orratio);
		}
		 ValueComparator bvc =  new ValueComparator();
			List<Map.Entry<String,Double>> mappingList=new ArrayList<Map.Entry<String,Double>>(PatternLib.entrySet());
			Collections.sort(mappingList, bvc);
			//just take top ones
		    for (Map.Entry<String,Double> pair : mappingList) {
		    	if(SeedMotifs.size()==max_num_Seeds)
		    		break;
			    String toppattern=pair.getKey();
			    System.out.println(toppattern+"\t"+pair.getValue());
			    StringBuffer sb=new StringBuffer(toppattern);
			    for (int j = 0; j < (max_motiflen-this.min_motiflen)/2; j++) {
			    	sb.insert(0, '-');
			    	sb.append('-');					
				}
			    
			    String[] seqs={sb.toString()};
			    try {
					PWM a=new PWM(seqs);
					a.Score=pair.getValue();
					SeedMotifs.add(a);
					//System.out.println(a.Score);
				} catch (IllegalAlphabetException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				} catch (IllegalSymbolException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
		    }
		return SeedMotifs;
	}
	
	public ArrayList<PWM>	 getSeedMotifs() {
		int max_num_Seeds=num_motif*num_motif;
		ArrayList<PWM>	SeedMotifs=new	ArrayList<PWM>(max_num_Seeds);
		int LIBSIZE=SearchEngine2.getTotalLength();
		String ACGT="ACGT";
		int loopnum=1<<(2*seedlen);
		int maxRange=1<<(2*seedlen);
		 ValueComparator bvc =  new ValueComparator();
		TreeMap<String,Double> seedScores=new TreeMap<String,Double>();
		for (String pattern:SearchEngine2.KmerHitList.keySet()) {
			if(seedScores.containsKey(common.getReverseCompletementString(pattern)))
				continue;
//			String pattern="";
//			int hash=i;
//			//ignore the reverseComplement
//			if(common.getReverseComplementHashing(hash,seedlen)<i)
//			{
//				continue;
//			}
//			else
//			pattern=common.Hash2ACGT(hash,seedlen);
			
			//LinkedList<Integer> positionlist=SearchEngine.searchPattern(pattern,0);
			LinkedList<FastaLocation> LocList=new LinkedList<FastaLocation>(); //SearchEngine.Int2Location(positionlist);

			int mismatch=0;
//			if(pattern.length()>=6)
//				mismatch=1;
//			if(pattern.length()>=8)
//				mismatch=2;
//			if(pattern.length()>=10)
//				mismatch=3;

			HashSet<String> candPatt=MismatchSet(pattern,mismatch);
			candPatt.addAll(MismatchSet(common.getReverseCompletementString(pattern),mismatch));
			double prob_bg=0;
			for(String candp:candPatt)
			{
				LinkedList<FastaLocation> templist=SearchEngine2.KmerHitList.get(candp);
				if(templist!=null)
				{
					LocList.addAll(templist);
					prob_bg+=Math.exp(background.Get_LOGPROB(candp));
				}
			}
			
			

		
			if((SearchEngine2.TotalLen*prob_bg)>LocList.size())
				continue;
//			double pvalue=BinomialDist.cdf(SearchEngine2.TotalLen, prob_bg,LocList.size());
//			if(pvalue<1-FDR)
//				continue;
			
			double score=0;
			int lastpos=-1;
			int facount_nonoverlap=0;
			for(FastaLocation fal:LocList)
			{
				if(fal.getMin()-lastpos>seedlen)
					facount_nonoverlap++;
				lastpos=fal.getMin();
			}
//			int overlapcount=LocList.size()-facount_nonoverlap;
//			double avglen=(SearchEngine2.TotalLen-LocList.size()*seedlen)/LocList.size();
//			double p=(double)max_motiflen/avglen;
//			double pv_dimer=BinomialDist.cdf(facount_nonoverlap,p,overlapcount);

			//if(pos_prior.size()==0&&!OOPS)
				score=facount_nonoverlap/(SearchEngine2.TotalLen*prob_bg); //positionlist.size()*(-common.DoubleMinNormal*seedlen-logprob_bg);//sum loglik ,-0.037267253272904234 is from pseudo count
//			if(pattern.equalsIgnoreCase("TTCCC"))
//				pattern="GGGAA";
			if(score<=1)
				continue;
				//else
//				score=sumLLR(LocList,-common.DoubleMinNormal*seedlen,logprob_bg);
//			score=BinomialDist.cdf(SearchEngine2.TotalLen,prob_bg,facount_nonoverlap);
			seedScores.put(pattern, score);
		}
		//sort by score
		List<Map.Entry<String,Double>> mappingList=new ArrayList<Map.Entry<String,Double>>(seedScores.entrySet());
		Collections.sort(mappingList, bvc);
		
		//just take top ones
		    for (Map.Entry<String,Double> pair : mappingList) {
		    	if(SeedMotifs.size()==max_num_Seeds)
		    		break;
			    String toppattern=pair.getKey();
			    System.out.println(toppattern+"\t"+pair.getValue());
			    StringBuffer sb=new StringBuffer(toppattern);
			    for (int j = 0; j < (max_motiflen-seedlen)/2; j++) {
			    	sb.insert(0, '-');
			    	sb.append('-');					
				}
			    
			    String[] seqs={sb.toString()};
			    try {
					PWM a=new PWM(seqs);
					a.Score=pair.getValue();
					SeedMotifs.add(a);
					//System.out.println(a.Score);
				} catch (IllegalAlphabetException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				} catch (IllegalSymbolException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
		    }
		 
		
		return	SeedMotifs;
	}
	

	public void DrawDistribution(ArrayList<Double> dist,String pngfile)
	{
		if(dist==null)
			return;
		  XYSeriesCollection dataset = new XYSeriesCollection();
		  XYSeries series1 = new XYSeries("");
		  for (int i = 0; i <dist.size(); i++) {
			double x=this.resolution*i-dist.size()*this.resolution/2;
			series1.add(x, dist.get(i));
		}
			 dataset.addSeries(series1);
		
		 JFreeChart chart = ChartFactory.createXYLineChart(
	                "Distribution curve", // chart title
	                "Position", // x axis label
	                "Probability", // y axis label
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
	        plot.setBackgroundPaint(Color.white);
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
	        	if(pngfile!=null)
				ChartUtilities.saveChartAsPNG(new File(pngfile), chart, 800, 600);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}

	}
	
	public double FindPrior(PWM motif,LinkedList<FastaLocation> Falocs)
	{
		double Prior_EZ=0;
		double [][]m_matrix=new double [motif.core_motiflen][4];
		int num_priorbin=SearchEngine2.getTotalLength()/SearchEngine2.getSeqNum()/this.resolution;
		double lognullprior=Math.log(1.0/num_priorbin);
		
		Collections.sort(Falocs,new FaScoreGreaterThan());
		
		Iterator<FastaLocation> iter=Falocs.iterator();
		ArrayList<String> sites=new ArrayList<String>(Falocs.size());
		ArrayList<Double> logbg_score=new ArrayList<Double>(Falocs.size());
		double bestscore=Double.NEGATIVE_INFINITY;
		PWM bestPWM=motif.Clone();
		double lamda=0;
		ArrayList<Double> probRatio = new ArrayList<Double>(Falocs.size());
		BGModel motifBG=new BGModel();
		TreeMap<String, Double> bgstrSet=new TreeMap<String, Double>();
		
	
		while(iter.hasNext())
		{
			
			FastaLocation currloc=iter.next();
			String site="";
			//forward site
			if(currloc.ReverseStrand)
				site=SearchEngine2.getSite(currloc.getSeqId(), currloc.getSeqPos(),motif.core_motiflen);
			else
				site=SearchEngine2.getSite(currloc.getSeqId(), currloc.getSeqPos(),motif.core_motiflen);
			//check masking
			if(site.indexOf('X')>-1)
				continue;
			
			//assume only one N for line break
			StringBuffer sb=new StringBuffer(site);
			for (int i = 0; i < site.length(); i++) {
				if(site.charAt(i)=='N'&& i<=site.length()/2)
				{
					for (int j = 0; j < i; j++) {
						sb.setCharAt(j, 'N');
					}
					continue;
				}
				else if(site.charAt(i)=='N')
				{
					for (int j = i+1; j < site.length(); j++) {
						sb.setCharAt(j, 'N');
					}
					break;
				}
				
			}
			site=sb.toString();

			if(currloc.ReverseStrand)
			{
				//reverse site
				site=common.getReverseCompletementString(site);
			}
			//double logprob_BG=background.Get_LOGPROB(site.substring((site.length()-motiflen)/2, motiflen));
//
			bgstrSet.put(site, 1.0);
			double logprob_BG=background.Get_LOGPROB(site);

			logbg_score.add(logprob_BG);//  //
			sites.add(site);
		}
//		motifBG.BuildModel(bgstrSet, 3);
//		try {
//			PWM pwmBG=new PWM(sites.toArray(new String[1]));
//
//		iter=Falocs.iterator();
//		while(iter.hasNext())
//		{
//
//			FastaLocation currloc=iter.next();
//			String site="";
//			//forward site
//			if(currloc.ReverseStrand)
//				site=SearchEngine2.getSite(currloc.getSeqId(), currloc.getSeqPos(),motif.core_motiflen);
//			else
//				site=SearchEngine2.getSite(currloc.getSeqId(), currloc.getSeqPos(),motif.core_motiflen);
//			//check masking
//			if(site.indexOf('X')>-1)
//				continue;
//			
//			//assume only one N for line break
//			StringBuffer sb=new StringBuffer(site);
//			for (int i = 0; i < site.length(); i++) {
//				if(site.charAt(i)=='N'&& i<=site.length()/2)
//				{
//					for (int j = 0; j < i; j++) {
//						sb.setCharAt(j, 'N');
//					}
//					continue;
//				}
//				else if(site.charAt(i)=='N')
//				{
//					for (int j = i+1; j < site.length(); j++) {
//						sb.setCharAt(j, 'N');
//					}
//					break;
//				}
//				
//			}
//			site=sb.toString();
//
//			if(currloc.ReverseStrand)
//			{
//				//reverse site
//				site=common.getReverseCompletementString(site);
//			}
//			//double logprob_BG=background.Get_LOGPROB(site.substring((site.length()-motiflen)/2, motiflen));
//			double logprob_BG=pwmBG.scoreWeightMatrix(site);//  motifBG.Get_LOGPROB(site);
//			logbg_score.add(logprob_BG);//  //
//			probRatio.add(Math.exp(currloc.Score-logprob_BG));
//			
//		}
//		
//		} catch (IllegalAlphabetException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		} catch (IllegalSymbolException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		}
//		FunctionOptimizer opt=new FunctionOptimizer(1);
//
//		FindPrior FPfun=new FindPrior(probRatio, logbg_score);
//		opt.setInitialX(new double[]{0.5});
//		opt.setOptimizationMode(OptimizationMode.Hessian);
//		opt.setFunctionController(new SimpleFunctionController());
////		opt.setFunctionController(new ParallelFunctionController());
//		opt.setFunctionAssembler(new AdditiveFunctionAssembler());
////		opt.optimizeFunction(FPfun);
////		lamda=opt.getOptimum()[0];
//		double bestlamda=0;
//		double minnegllk=Double.MAX_VALUE;
//		for (lamda=0.01;lamda < 1; lamda+=0.01) {
//			double negllk=FPfun.evaluate(new double[]{lamda}, 1);
//			if(negllk<minnegllk)
//			{
//				minnegllk=negllk;
//				bestlamda=lamda;
//			}
//		}
		
		
		int seqcount=0;
		Iterator<FastaLocation> iter2=Falocs.iterator();
		while(iter2.hasNext())
		{
			seqcount++;
			Prior_EZ=(double)seqcount/SearchEngine2.TotalLen;
			FastaLocation currloc=iter2.next();
			String site=sites.get(seqcount-1);

			for (int i = 0; i < motif.core_motiflen; i++) {
				int symid=common.acgt[site.charAt(i)];
				if(symid<4)
					m_matrix[i][symid]+=1;
				else
				{
					for (int j = 0; j < 4; j++) {
						m_matrix[i][j]+=0.25;
					}
				}	
			}
			
			for (int i = 0; i < m_matrix.length; i++) {
				motif.setWeights(i+motif.head,common.Normalize(m_matrix[i].clone()));
			}
			
			Iterator<FastaLocation> iter3=Falocs.iterator();
			int seqcount2=0;
			double loglikscore=0;
			while(iter3.hasNext())
			{
				seqcount2++;
				FastaLocation currloc2=iter3.next();
				String site2=sites.get(seqcount2-1);
				int prior_bin2=(int)(num_priorbin*((currloc2.getSeqPos()+motif.core_motiflen/2)%currloc2.getSeqLen()/(double)currloc2.getSeqLen()));
				int rankbin2=num_priorbin*currloc2.getSeqId()/SearchEngine2.getSeqNum();
				double logprior2=0;
				double logprob_theta2=motif.scoreWeightMatrix(site2)-logbg_score.get(seqcount2-1);//  //
				if(motif.pos_prior.size()!=0&&motif.pos_en)
					logprior2=Math.log(motif.pos_prior.get(prior_bin2)+common.DoubleMinNormal)-lognullprior;
				if(motif.strand_en)
				{
					if(currloc.ReverseStrand)
						logprior2+=Math.log((1-motif.strand_plus_prior)*2);
					else
					logprior2+=Math.log(motif.strand_plus_prior*2);
				}
				if(motif.peakrank_prior.size()!=0&&motif.peakrank_en)
					logprior2+=Math.log(motif.peakrank_prior.get(rankbin2)+common.DoubleMinNormal)-lognullprior;
				
				double logDnaseProb2=0;
				//logprior=0;
				double loglik2=logprob_theta2+logprior2+logDnaseProb2+Math.log(Prior_EZ/(1-Prior_EZ));
				if(loglik2>10)
					loglik2=10;
				double prob_theta2=Math.exp(loglik2)/(Math.exp(loglik2)+1);//Math.exp(currloc.Score);
				loglikscore+=prob_theta2*loglik2;
			}
			
			if(loglikscore>bestscore)
			{
				bestPWM=motif.Clone();
				bestscore=loglikscore;
				lamda=Prior_EZ;
			}
			
		}
		
		for (int i = 0; i < m_matrix.length; i++) {
			motif.setWeights(i+motif.head,common.Normalize(bestPWM.m_matrix[i]));
		}
		
		return lamda;
	}
	
	public PWM Relax_Seed_3(PWM motif)
	{
		double log025=Math.log(0.25);
		double switchvalue=1;//*SearchEngine2.TotalLen*Math.pow(0.25, motif.core_motiflen);
		//relax the conserved column 
		String consensus=motif.Consensus(false);
		double main_prop=0.504166667;//Math.pow(sampling_ratio*9/(seedlen*seedlen-2*seedlen+4), 1.0/(seedlen-2)); //2 mismatch
		if(motif.core_motiflen<10)
			main_prop= 0.704166667;
			//Math.pow(3*sampling_ratio/(seedlen-2), 1.0/(seedlen-1));// 1 mismatch
			
			//Math.pow(sampling_ratio*9/(seedlen*seedlen-2*seedlen+4), 1.0/(seedlen-2)); //2 mismatch
		for (int i = 0; i < consensus.length(); i++) {
			if(consensus.charAt(i)=='N')
				continue;
			boolean conserved=false;
			for (int j = 0; j < 4; j++) {
				 if(motif.m_matrix[i][j]>0.9999)
				 {
					 conserved=true;
					 break;
				 }
			}
			if(conserved)
			{
				for (int j = 0; j < 4; j++) {
					 if(motif.m_matrix[i][j]>main_prop)
					 {
					    motif.setWeight(i, j, main_prop);
						
					 }
					 else
						motif.setWeight(i, j, (1-main_prop)/3);
				}
				
			}
		}

		//EM full site iteration
		
		double bestscore=motif.Score;
		double lastscore=Double.NEGATIVE_INFINITY;
		int num_priorbin=SearchEngine2.getTotalLength()/SearchEngine2.getSeqNum()/this.resolution;
		if(motif.core_motiflen/2>this.resolution)
			num_priorbin=SearchEngine2.getTotalLength()/SearchEngine2.getSeqNum()/(motif.core_motiflen/2);
		
		if(motif.pos_prior.size()==0)
		{
			for (int i = 0; i <num_priorbin ; i++) {
				motif.pos_prior.add(1.0/num_priorbin);
			}
		}
		
		
		HashSet<Integer> stateCodes=new HashSet<Integer>();
	
		int flankingLen=0;
		int orighead=motif.head;
		int origtail=motif.tail;

		 double    Prior_EZ=(double)SearchEngine2.getSeqNum()*0.5/SearchEngine2.TotalLen;//motif.inst_coverage/Falocs.size();//1-SearchEngine.TotalLen*prior_fp/Falocs.size();
		   Prior_EZ=Math.min(Prior_EZ, motif.inst_coverage/SearchEngine2.TotalLen);	
		   
		int iter_count=0;
		PWM bestPWM=motif.Clone();
		bestPWM.Score=0;//ignore previous score
		double sitesperSeq=0;
		LinkedList<FastaLocation> Falocs=null;
//		if(OOPS)
//		{
//			SamplingThread_PS.background=this.background;
//			Falocs=SearchEngine2.samplingPattern_PS(motif,Math.max(MIN_SAMPLENUM,(int)(SearchEngine2.TotalLen*sampling_ratio)));
//		}
//		else
//		{
			SamplingThread.background=this.background;
			Falocs=SearchEngine2.samplingPattern(motif,Math.max(MIN_SAMPLENUM,(int)(SearchEngine2.TotalLen*sampling_ratio)));
//		}
	
			
		
//			Falocs=SearchEngine2.searchPattern(motif, motif.core_motiflen*Math.log(0.25)+Math.log(2.0));

			 
		int newhead=motif.head;
		int newtail=motif.tail;
		if(motif.head+flankingLen+motif.core_motiflen>=motif.columns()||motif.head-flankingLen<0)
			flankingLen=0;
		
		//Prior_EZ=Math.min(0.5,(double)SearchEngine.getSeqNum()/Falocs.size());
		///////////////build newBG for iterations////////////////////
		BGModel motifBG=new BGModel();
		double [][] bgprob=new double[motif.core_motiflen][4];
		Iterator<FastaLocation> iter=Falocs.iterator();
		ArrayList<FastaLocation> filtered_Falocs=new ArrayList<FastaLocation>(Falocs.size());
		TreeMap<String, Double> bgstrSet=new TreeMap<String, Double>();
		int last=-1;
		int lastpos=-1;
		double trueweight=0;
		int seqcount=0;
		double total_sampleweight=0;
		double [] single_bgprob=new double[4];
		
		double [] peakrank_renorm=new double[num_priorbin];
		double [] strand_renorm=new double[2];
		double [] pos_renorm=new double[num_priorbin];
		int truecount=0;
		double duplicateFactor=1;
		int scount=0;
		double prob_fsum=0;
		double prob_rsum=0;
		ArrayList<String> Siteslist=new ArrayList<String>(Falocs.size());
		while(iter.hasNext())
		{
			FastaLocation currloc=iter.next();
			if(currloc.getSeqId()!=last)
			{
				seqcount++;
				last=currloc.getSeqId();
			}
			String site="";
			//forward site
			if(currloc.ReverseStrand)
				site=SearchEngine2.getSite(currloc.getSeqId(), currloc.getSeqPos()-(origtail-newtail),motif.core_motiflen);
			else
				site=SearchEngine2.getSite(currloc.getSeqId(), currloc.getSeqPos()-(orighead-newhead),motif.core_motiflen);
//			check masking
//			if(site.indexOf('X')>-1)
//				continue;
			
			
//			if(site.equalsIgnoreCase(""))
//				continue;

			double probtheta=Math.exp(currloc.Score);
//			int rankbin=num_priorbin*currloc.getSeqId()/SearchEngine2.getSeqNum();
//			int prior_bin=(int)(num_priorbin*((currloc.getSeqPos()+motif.core_motiflen/2)%currloc.getSeqLen()/(double)currloc.getSeqLen()));
//			peakrank_renorm[rankbin]+=1;//probtheta;
//			pos_renorm[prior_bin]+=1;//probtheta;
			if(currloc.ReverseStrand)
				strand_renorm[1]+=1;
			else
				strand_renorm[0]+=1;
			//assume only one N for line break
			StringBuffer sb=new StringBuffer(site);
			for (int i = 0; i < site.length(); i++) {
				if(site.charAt(i)=='N'&& i<=site.length()/2)
				{
					for (int j = 0; j < i; j++) {
						sb.setCharAt(j, 'N');
					}
					continue;
				}
				else if(site.charAt(i)=='N')
				{
					for (int j = i+1; j < site.length(); j++) {
						sb.setCharAt(j, 'N');
					}
					break;
				}
				
			}
			site=sb.toString();
			if(currloc.ReverseStrand)
			{
				//reverse site
				site=common.getReverseCompletementString(site);
			}

			double entropy=common.SeqComplexity(2,site);
			if(entropy<1)
				continue;
			if(site.length()!=motif.core_motiflen)
				continue;
//			if(lastpos==currloc.getMin())
//			{
//				filtered_Falocs.get(filtered_Falocs.size()-1).Score+=currloc.Score;
//				continue;
//			}
			
			Siteslist.add(site);
			filtered_Falocs.add(currloc);
			
			double sampleweight=1.0/probtheta;

			boolean truesite=false;
			boolean overlapflag=Math.abs(currloc.getMin()-lastpos)<motif.core_motiflen;
			if((Character.isLowerCase( site.charAt(0))&&Character.isLowerCase( site.charAt(site.length()-1))))//
			{
				
				truecount++;
				truesite=true;
				//if(!site.matches("A|C|G|T"))
				trueweight+=sampleweight;
				if(currloc.ReverseStrand)
				{
					scount++;
					prob_rsum+=1.0/sampleweight;
				}
				else
					prob_fsum+=1.0/sampleweight;
			
			}

			lastpos=currloc.getMin();
			for (int i = 0; i < site.length(); i++) {
				int symid=common.acgt[site.charAt(i)];
				if(symid<4)
				{
//				if(!truesite&&!overlapflag)
				single_bgprob[symid]+=(1-probtheta)*sampleweight;
				bgprob[i][symid]+=(1-probtheta)*sampleweight;
				}
				else
				{
					for (int k = 0; k < 4; k++) {
						bgprob[i][k]+=0.25*sampleweight;
					}
				}
			}
			bgstrSet.put(site, sampleweight);
			if(sampleweight>switchvalue)
				sampleweight=switchvalue;
			total_sampleweight+=sampleweight;
			
		}
		peakrank_renorm=common.Normalize(peakrank_renorm);
		pos_renorm=common.Normalize(pos_renorm);
		strand_renorm=common.Normalize(strand_renorm);
		
		double prior_true=trueweight/total_sampleweight;
		duplicateFactor=(double)filtered_Falocs.size()/bgstrSet.size();
//		if(OOPS)
//			duplicateFactor=(double)filtered_Falocs.size()/seqcount;
		int bgorder=1;
		/****************manual initialize********************/
		motifBG.order=1;
		single_bgprob[0]+=single_bgprob[3];
		single_bgprob[1]+=single_bgprob[2];
		single_bgprob[2]=single_bgprob[1];
		single_bgprob[3]=single_bgprob[0];
		single_bgprob=common.Normalize(single_bgprob);
		motifBG.conditionProb=new HashMap<String, Double>(4);
		motifBG.conditionProb.put("A", single_bgprob[0]);
		motifBG.conditionProb.put("C", single_bgprob[1]);
		motifBG.conditionProb.put("G", single_bgprob[2]);
		motifBG.conditionProb.put("T", single_bgprob[3]);
		
//		motifBG.conditionProb.put("A", 0.25);
//		motifBG.conditionProb.put("C", 0.25);
//		motifBG.conditionProb.put("G", 0.25);
//		motifBG.conditionProb.put("T", 0.25);
		
		Arrays.fill(single_bgprob, 0);
		/****************manual initialize********************/
		
//		motifBG.BuildModel(bgstrSet, bgorder);
		bgstrSet.clear();
		PWM PWMBG=null;
		try {
			PWMBG=new PWM(bgprob);
			common.fill2DArray(bgprob, 0);
		} catch (IllegalAlphabetException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IllegalSymbolException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		
		//int truepos=PWMevaluator.comparePositionList(Falocs, "D:\\eclipse\\data\\batchsim\\ESR1.ans", motif.core_motiflen);
		//Prior_EZ=(double)truepos/Falocs.size();
		//double prior_fp=motif.getFDR(log_thresh,background);
   
    //Prior_EZ=Prior_EZ*SearchEngine2.TotalLen*seqcount/SearchEngine2.getSeqNum()/filtered_Falocs.size();
    
   RandomEngine rand=RandomEngine.makeDefault();
    double max_Prior_EZ=(double)SearchEngine2.getSeqNum()*MAX_P/SearchEngine2.TotalLen;
        if(Prior_EZ<0)
			Prior_EZ=FDR;
		//Prior_EZ*=0.9;
		double prior_gamma=(Prior_EZ*Falocs.size())/seqcount;//(double)truepos/seqcount;//
		if(prior_gamma>1)
			prior_gamma=0.9999;
		///////////////build newBG for iterations////////////////////
		
		////////////////OOPS correction/////////////////
			
		////////////////OOPS correction/////////////////
		
		int motiflen=motif.core_motiflen+flankingLen*2; 
//		Prior_EZ=(double)truepos/Falocs.size();

		do
		{
			
			iter_count++;
			
//			common.print2DArray(PWMBG.m_matrix);
			
			String consensus_core=motif.Consensus(false).substring(motif.head,motif.head+motif.core_motiflen);
			System.out.println(consensus_core+"\t"+String.valueOf(bestscore));
			
			bestscore=0;
			double[] temp_prior=new double[num_priorbin];
			double[] temp_peakrank=new double[num_priorbin];	
			double[] temp_strand=new double[2];	
			TreeMap<Double, Integer> Zvalue=new TreeMap<Double, Integer>();
			double overlap_loglik=0;
			double overlap_expLLR=0;
			double overlap_prob=0;
			
			double matchcount=0;
			int overlap_pos=-motiflen;
			double lognullprior=Math.log(1.0/num_priorbin);
			//make sure head, tail the same in the iteration
			motif.head=newhead;
			motif.tail=newtail;
			double maxsamplweight=0;
			//double log_thresh=motif.getThresh(sampling_ratio, 2*FDR, background)- motiflen*log025;
			//double log_thresh=Math.log(1-Prior_EZ)-Math.log(Prior_EZ);
			
			double [][]m_matrix=new double [motiflen][4];
			double [][] sumexpLLR=new double[motiflen][4];
			
			//update the loglik matrix
			//LinkedList<FastaLocation> Falocs=SearchEngine2.searchPattern(motif, log_thresh);
			System.out.println("number of occurrences: "+String.valueOf(Falocs.size()));
			//System.out.println(motifBG.conditionProb);
					Iterator<FastaLocation> iter2=filtered_Falocs.iterator();
					Iterator<String> iter_site=Siteslist.iterator();
					int count=0;
					double match_seqCount=0;
					int lastseq=-1;
					lastpos=-1;
					double lastpwmweight=0;
					double lastprob_theta=0;
					String lastsite="";
					double max_seqloglik=0;
                    double matchsitecount_seq=0;
                    double matchsiteweight_seq=0;
					String max_seqsite="";
					while(iter2.hasNext())
					{
						FastaLocation currloc=iter2.next();
						String site=iter_site.next();
						String revsite=common.getReverseCompletementString(site);
						
						//double logprob_theta=currloc.Score;//include the bg log_prob in the score
						double logpwm=motif.scoreWeightMatrix(site);//-Math.log(2.0);
						double logbg=PWMBG.scoreWeightMatrix(site);//motifBG.Get_LOGPROB(site);//motifBG.Get_LOGPROB(site)+Math.log(2.0);;//
						logbg=Math.max(logbg, background.Get_LOGPROB(site));
						double llc=Math.log(Math.exp(logpwm)*Prior_EZ+(1-Prior_EZ)*Math.exp(logbg))/10000;
						double logprob_theta=logpwm-logbg;//  //
					//	double logprob_theta=motif.scoreWeightMatrix(site)-SearchEngine2.BGscoreMap.get(motif.core_motiflen).get(currloc.getSeqId()).get(currloc.getSeqPos());

					
						//double logprob_BG=background.Get_LOGPROB(site.substring((site.length()-motiflen)/2, motiflen));
						double logprior=0;
                   ///////////////////////////// extra feature integration ////////////////////////////
						int prior_bin=(int)(num_priorbin*((currloc.getSeqPos()+motiflen/2)%currloc.getSeqLen()/(double)currloc.getSeqLen()));
						int rankbin=num_priorbin*currloc.getSeqId()/SearchEngine2.getSeqNum();
						if(motif.pos_prior.size()!=0&&motif.pos_en)
							logprior=Math.log(motif.pos_prior.get(prior_bin)+Double.MIN_NORMAL)-lognullprior;
						if(motif.strand_en)
						{
							if(currloc.ReverseStrand)
								logprior+=Math.log((1-motif.strand_plus_prior)*2);
							else
								logprior+=Math.log(motif.strand_plus_prior*2);
						}
						if(motif.peakrank_prior.size()!=0&&motif.peakrank_en)
							logprior+=Math.log(motif.peakrank_prior.get(rankbin)+Double.MIN_NORMAL)-lognullprior;
						
						double logDnaseProb=0;
                    ///////////////////////////// extra feature integration ////////////////////////////
						//logprior=0;
						double loglik=logprob_theta+logprior+logDnaseProb+Math.log(Prior_EZ/(1-Prior_EZ));
						if(loglik>10)
							loglik=10;
						double prob_theta=Math.exp(loglik)/(Math.exp(loglik)+1);//Math.exp(currloc.Score);
						
						//double prob_theta_only=Math.exp(logprob_theta+logprior)/(Math.exp(logprob_theta+logprior)+1);
						if(Double.isNaN(prob_theta))
							prob_theta=1;//upper flow
						double sampleWeight=1.0/Math.exp(currloc.Score); //re-weighting
					//	sampleWeight=1;
//						double we=prob_theta*sampleWeight;
//						if(we>10)
//							System.out.println(site+" "+we+" "+sampleWeight);
						
//						bgstrSet.put(site, (1-prob_theta)*sampleWeight);
						

						
						temp_peakrank[rankbin]+=prob_theta/(peakrank_renorm[rankbin]*num_priorbin);
						temp_prior[prior_bin]+=prob_theta/(pos_renorm[prior_bin]*num_priorbin);//make smaller
						if(currloc.ReverseStrand)
						{						
							temp_strand[1]+=prob_theta*0.5/strand_renorm[1];
						}
						else
						{						
							temp_strand[0]+=prob_theta*0.5/strand_renorm[0];
						}
						if(OOPS)
							loglik-=common.DoubleMinNormal*Math.abs(currloc.getSeqLen()/2-currloc.getSeqPos()-motiflen/2); //add small bias to center
//						if(!OOPS)
							bestscore+=llc*sampleWeight;//(loglik-Math.log(prob_theta))*sampleWeight/total_sampleweight;
					
						if(Double.isInfinite(bestscore))
							break;

						Zvalue.put(1-prob_theta, count);
						count++;
						
						if(sampleWeight>switchvalue)
						sampleWeight=switchvalue;
						
						if(!OOPS)
						{
							matchcount+=prob_theta*sampleWeight;
							double pwmweight=prob_theta*sampleWeight;
							if(maxsamplweight<pwmweight)
								maxsamplweight=pwmweight;


								if((currloc.getMin()-lastpos)<site.length()&&prob_theta>(lastprob_theta+common.DoubleMinNormal))
								{
									for (int i = 0; i < lastsite.length(); i++) {
										int symid=common.acgt[lastsite.charAt(i)];
										if(symid>3)
										{
											continue;
										}
										
										m_matrix[i][symid]-=lastpwmweight*0.5;
									}
								}
								

							if((currloc.getMin()-lastpos)>=site.length()||prob_theta>(lastprob_theta+common.DoubleMinNormal))
							{
								
									for (int i = 0; i < site.length(); i++) {
										int symid=common.acgt[site.charAt(i)];
										if(symid>3)
										{
											continue;
										}	
										
										m_matrix[i][symid]+=prob_theta*sampleWeight;//prob_theta;//
										single_bgprob[symid]+=(1-prob_theta)*sampleWeight;//(1-prob_theta)*
										bgprob[i][symid]+=(1-prob_theta)*sampleWeight;	//(1-prob_theta)*
										
									}
									lastpos=currloc.getMin();
									lastsite=site;
									lastprob_theta=prob_theta;
									lastpwmweight=pwmweight;
								}

						}
						if(OOPS&&currloc.getSeqId()!=lastseq)
						{
							//NegBinFunction.plogis(max_seqloglik);
							lastseq=currloc.getSeqId();
//							if(matchsitecount_seq>currloc.getSeqLen())
//								matchsitecount_seq=currloc.getSeqLen();
							double seqlen=(currloc.getSeqLen()-motif.core_motiflen);
							prior_gamma=Prior_EZ*matchsitecount_seq;//
							if(prior_gamma>1)
								prior_gamma=0.9999;
							sitesperSeq=0;
							if(matchsitecount_seq>0)
							{
								 double renomalizefactor=matchsitecount_seq/matchsiteweight_seq;
								for (int i = 0; i < motiflen; i++) {
									double sumexpLLRallsymid=0;
									 for (int symid = 0; symid < 4; symid++) 
									 {
										 sumexpLLRallsymid+=sumexpLLR[i][symid];
									 }
   								    for (int symid = 0; symid < 4; symid++) 
									{
									//m_matrix[i][symid]+=max_count_matrix[i][symid]/sitesperSeq;
   									
   									 double temp=sumexpLLR[i][symid]*Prior_EZ/((1-prior_gamma)/renomalizefactor+sumexpLLRallsymid*Prior_EZ);
											m_matrix[i][symid]+=temp;
											sitesperSeq+=temp;
									}
								
								}
								
								
							}
							
//							if(matchsitecount_seq>0)
//								bestscore+= max_seqloglik/matchsitecount_seq;

							max_seqloglik=0;	
							sitesperSeq=0;
							overlap_prob=0;
							overlap_pos=-motiflen;
							common.fill2DArray(sumexpLLR,0);
							matchsitecount_seq=0;
							matchsiteweight_seq=0;
						}
						
						if(OOPS)
						{
							double overlapfactor=1;
							if(currloc.getSeqPos()-overlap_pos<motif.core_motiflen)
							{
								
								if(overlap_prob<prob_theta)
								{
									overlapfactor=1;
									for (int i = 0; i < lastsite.length(); i++) {
										int symid=common.acgt[lastsite.charAt(i)];									
										if(symid>3)
										{
											for (int j = 0; j < 4; j++) {
												sumexpLLR[i][j]-=0.25*overlap_expLLR*overlapfactor; //re-weighting
											}
											continue;
										}
									sumexpLLR[i][symid]-=overlap_expLLR*overlapfactor; //re-weighting
									}
								}
								else
									overlapfactor=0;
							}

							
							matchsitecount_seq++;//=sampleWeight;//
							matchsiteweight_seq+=sampleWeight;
							double expLLR=Math.exp(loglik-Math.log(Prior_EZ/(1-Prior_EZ)));
							
								max_seqloglik+=prob_theta*loglik*sampleWeight; //re-weighting
								if(Double.isNaN(max_seqloglik))
									break;
								for (int i = 0; i < site.length(); i++) {
									int symid=common.acgt[site.charAt(i)];									
									if(symid>3)
									{
										for (int j = 0; j < 4; j++) {
											sumexpLLR[i][j]+=0.25*expLLR*sampleWeight*overlapfactor; //re-weighting
											bgprob[i][j]+=0.25*(1-prob_theta)*sampleWeight;
										}
										continue;
									}
									
								single_bgprob[symid]+=(1-prob_theta)*sampleWeight;
								bgprob[i][symid]+=(1-prob_theta)*sampleWeight;
								sumexpLLR[i][symid]+=expLLR*sampleWeight*overlapfactor; //re-weighting
								
								}	
								sitesperSeq+=prob_theta;

								if(currloc.getSeqPos()-overlap_pos>motif.core_motiflen||overlap_prob<prob_theta)
								{
									overlap_prob=prob_theta;
									overlap_loglik=loglik;
									overlap_expLLR=expLLR*sampleWeight;
									lastsite=site;
									overlap_pos=currloc.getSeqPos();
								}

						}
						
						

					}
					
					
					//last instance for OOPS
					if(OOPS)
					{
						double seqlen=(SearchEngine2.TotalLen/SearchEngine2.getSeqNum()-motif.core_motiflen);
						prior_gamma=Prior_EZ*matchsitecount_seq;//
						 double renomalizefactor=matchsitecount_seq/matchsiteweight_seq;
						if(prior_gamma>1)
							prior_gamma=0.9999;
						for (int i = 0; i < motiflen; i++) {
							double sumexpLLRallsymid=0;
							
							for (int symid = 0; symid < 4; symid++) 
							{
								sumexpLLRallsymid+=sumexpLLR[i][symid];
							}
							
							for (int symid = 0; symid < 4; symid++) 
							{
								
							  //m_matrix[i][symid]+=max_count_matrix[i][symid]/sitesperSeq+common.DoubleMinNormal;// pesudo count
								double temp=sumexpLLR[i][symid]*Prior_EZ/((1-prior_gamma)/renomalizefactor+sumexpLLRallsymid*Prior_EZ);
								m_matrix[i][symid]+=temp;
								sitesperSeq+=temp;
								
							}
						}

//							bestscore+= max_seqloglik/matchsitecount_seq;
					}
					
					match_seqCount=0;
					for (int i = 0; i < 4; i++) {
						match_seqCount+=m_matrix[m_matrix.length/2][i];
					}
					//total_sampleweight may introduce more fluctuation
					if(OOPS)
						Prior_EZ=match_seqCount/total_sampleweight;//filtered_Falocs.size();//(SearchEngine2.TotalLen*(double)seqcount/SearchEngine2.getSeqNum());//SearchEngine2.TotalLen;
					else
						Prior_EZ=(double)matchcount/total_sampleweight;//match_seqCount/total_sampleweight;//SearchEngine2.TotalLen;  //total_sampleweight;
					
					
//					if(Prior_EZ>max_Prior_EZ)
//						Prior_EZ=max_Prior_EZ;
					System.out.println(Prior_EZ);
					//System.out.println(match_seqCount/total_sampleweight);
					if(prior_gamma>1)
						prior_gamma=0.9999;
					

//					motifBG.BuildModel(bgstrSet, bgorder);		
					/****************manual initialize********************/
					single_bgprob[0]+=single_bgprob[3];
					single_bgprob[1]+=single_bgprob[2];
					single_bgprob[2]=single_bgprob[1];
					single_bgprob[3]=single_bgprob[0];
					single_bgprob=common.Normalize(single_bgprob);
					motifBG.conditionProb.put("A", single_bgprob[0]);
					motifBG.conditionProb.put("C", single_bgprob[1]);
					motifBG.conditionProb.put("G", single_bgprob[2]);
					motifBG.conditionProb.put("T", single_bgprob[3]);
					
//					motifBG.conditionProb.put("A", 0.25);
//					motifBG.conditionProb.put("C", 0.25);
//					motifBG.conditionProb.put("G", 0.25);
//					motifBG.conditionProb.put("T", 0.25);
					Arrays.fill(single_bgprob, 0);
					/****************manual initialize********************/
					try {
						PWMBG=new PWM(bgprob);
						common.fill2DArray(bgprob, 0);
					} catch (IllegalAlphabetException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					} catch (IllegalSymbolException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
					
					
					
					
					if(Double.isNaN(bestscore) )//||Math.abs(lastscore-bestscore)<FDR||match_seqCount<min_support_ratio*SearchEngine2.getSeqNum()
						break;
					else
					{

						motif.Score=bestscore;//-Math.log(SearchEngine2.getSeqNum())*2*motif.core_motiflen;
						lastscore=bestscore;

					

						  bestPWM=motif.Clone();
						for (int i = 0; i < m_matrix.length; i++) {
							motif.setWeights(i+motif.head,common.Normalize(m_matrix[i]));
						}
						double divergence=common.PWM_Divergence(bestPWM,motif);
						bestPWM.Prior_EZ=Prior_EZ;//*filtered_Falocs.size()/SearchEngine2.TotalLen;
						if(OOPS)
							bestPWM.Prior_EZ=Prior_EZ*filtered_Falocs.size()/(seqcount*SearchEngine2.TotalLen/SearchEngine2.getSeqNum());
						if(divergence<0.001)
							break;
						
//						if(OOPS)
//						{
//							double sumprob=temp_strand[0]+temp_strand[1];
//							if(sumprob>seqcount)
//							{
//								duplicateFactor=sumprob/seqcount;
//							}
//
//						}
						//determine whether pos_prior is significant needed
						double chistat=0;
						double sumPrior=0;
						for (int i = 0; i < temp_prior.length; i++) {
//							temp_peakrank[i]/=duplicateFactor;
							sumPrior+=temp_prior[i];
						}
						//ignore the first and last bin
						double avgE=(sumPrior-temp_prior[0]-temp_prior[temp_prior.length-1])/(temp_prior.length-2);
							for (int j = 1; j < temp_prior.length-1; j++) {
								double temp=temp_prior[j]-avgE;
								chistat+=temp*temp/avgE;
							}
							double invFDR=ChiSquareDistQuick.inverseF(temp_prior.length-3, 1-FDR);
							
							if(chistat>invFDR)
							{
								if(!motif.pos_en)
								System.out.println("position bias on...");
								motif.pos_en=true;
		
							}
							else
							{
								if(motif.pos_en)
								System.out.println("position bias off...");
								motif.pos_en=false;
							}
							
							
							
							motif.pos_prior.clear();
							for (int i = 0; i < temp_prior.length; i++) {
								temp_prior[i]/=sumPrior;
								motif.pos_prior.add(temp_prior[i]);
							}
							//determine whether strand_prior is significant needed
							
//							temp_strand[0]/=duplicateFactor;
//							temp_strand[1]/=duplicateFactor;
							double X=Math.max(temp_strand[0], temp_strand[1])+1;
							Binomial binomial=new Binomial((int)(temp_strand[0]+temp_strand[1])+2,0.5,rand);
							double pvalue_strand=1.0-binomial.cdf((int)X);
							if(pvalue_strand<this.FDR)
							{
								if(!motif.strand_en)
								System.out.println("strand bias on...");
								motif.strand_plus_prior=temp_strand[0]/(temp_strand[0]+temp_strand[1]);
								motif.strand_en=true;
							}
							else
							{
								if(motif.strand_en)
								System.out.println("strand bias off...");
								motif.strand_en=false;
							}
						
							//determine whether peakrank_prior is significant needed, and ignore the final one
							chistat=0;
							sumPrior=0;
							for (int i = 0; i < temp_peakrank.length; i++) {
//								temp_peakrank[i]/=duplicateFactor;
								sumPrior+=temp_peakrank[i];
							}
							avgE=(sumPrior-temp_peakrank[temp_peakrank.length-1])/(temp_peakrank.length-1);
							
								for (int j = 0; j < temp_peakrank.length-1; j++) {
									double temp=temp_peakrank[j]-avgE;
									chistat+=temp*temp/avgE;
								}
								invFDR=ChiSquareDistQuick.inverseF(temp_peakrank.length-2, 1-FDR);
								if(chistat>invFDR)
								{
									if(!motif.peakrank_en)
									System.out.println("rank bias on...");
									motif.peakrank_en=true;
								}
								else
								{
									if(motif.peakrank_en)
									System.out.println("rank bias off...");
									motif.peakrank_en=false;
								}
								motif.peakrank_prior.clear();
								for (int i = 0; i < temp_peakrank.length; i++) {
									temp_peakrank[i]/=sumPrior;
									motif.peakrank_prior.add(temp_peakrank[i]);
								}					
					}
					//DrawDistribution(motif.peakrank_prior,"temppeakrank.png");
					//not allow to grow in the iterations
					flankingLen=0;
//					if(match_seqCount<1)
//						break;

					System.out.println("number of occurred sequences: "+String.valueOf(match_seqCount));

			}while(iter_count<=max_iterNum);
			
		if(DnaseLib!=null)
			DrawDistribution(motif.Dnase_prob,"Dnase_plot.png");
		
		return bestPWM;
	}

	
	
	
	public PWM Column_Replacement_2(PWM motif)
	{
		int maxiter=50;
		int minmotiflen=min_motiflen;
		int maxmotiflen=(this.max_motiflen+seedlen)/2;
		double log025=Math.log(0.25);
		if(this.pos_prior.size()>0)
			motif.pos_prior=(ArrayList<Double>) this.pos_prior.clone();
		double bestscore=motif.Score;
		HashSet<Integer> extendedCols=new HashSet<Integer>();
		HashSet<Integer> stateCodes=new HashSet<Integer>();
		int num_priorbin=SearchEngine2.getTotalLength()/SearchEngine2.getSeqNum()/this.resolution;
		double Prior_EZ=0.5;
		int inst_hash=-1;
		String seedstring=motif.Consensus(true);

		//double thresh=motif.getThresh(1.0, this.FDR, this.background);
		LinkedList<FastaLocation> origFalocs=new LinkedList<FastaLocation>();//SearchEngine.searchPattern(motif, thresh);
		LinkedList<FastaLocation> tempfalocs=(SearchEngine2.KmerHitList.get(seedstring));
		if(seedstring.length()!=SearchEngine2.KmerHitList.keySet().iterator().next().length())
			tempfalocs=SearchEngine2.searchPattern(seedstring, 0);
		////////////////////////////////////Insert the reverse complement///////////////////////////////////////////
		if(!seedstring.equalsIgnoreCase(common.getReverseCompletementString(seedstring)))
		{
			LinkedList<FastaLocation> tempfalocs2=SearchEngine2.KmerHitList.get(common.getReverseCompletementString(seedstring));
			if(tempfalocs==null)
			{
				Iterator<FastaLocation> iter2=tempfalocs2.iterator();
				while(iter2.hasNext())
				{
					FastaLocation pos2=iter2.next();
					pos2.ReverseStrand=true;
					origFalocs.add(pos2);
				}
			}
			else if(tempfalocs2==null)
				origFalocs.addAll(tempfalocs);
			else
			{
				Iterator<FastaLocation> iter1=tempfalocs.iterator();
				Iterator<FastaLocation> iter2=tempfalocs2.iterator();
				FastaLocation pos1 = iter1.next();
				FastaLocation pos2 = iter2.next();
				while(iter2.hasNext()||iter1.hasNext())
				{
					if(pos1.getMin()<pos2.getMin())
					{
						origFalocs.add(pos1);
						if(iter1.hasNext())
						pos1=iter1.next();
					}
					else 
					{
						pos2.ReverseStrand=true;
						origFalocs.add(pos2);
						if(iter2.hasNext())
						pos2=iter2.next();
					}
					if(!iter2.hasNext())
					{
						while(iter1.hasNext())
						{
							pos1=iter1.next();
							origFalocs.add(pos1);
						}
							
					}
					if(!iter1.hasNext())
					{
						while(iter2.hasNext())
						{
							pos2=iter2.next();
							pos2.ReverseStrand=true;
							origFalocs.add(pos2);
						}
							
					}
					
				}
			}
		
		}
		else
		{
			origFalocs.addAll(tempfalocs);
		}
		////////////////////////////////////Insert the reverse complement///////////////////////////////////////////
		
		if(origFalocs.size()<10)
			return null;
//		int truepos=PWMevaluator.comparePositionList(origFalocs, "E:\\eclipse\\data\\batchsim\\vali_Klf4.ans", motif.core_motiflen*2);
		Prior_EZ=1- (Math.exp(background.Get_LOGPROB(motif.Consensus(true)))+Math.exp(background.Get_LOGPROB(common.getReverseCompletementString(motif.Consensus(true)))))*SearchEngine2.TotalLen/origFalocs.size();
		if(Prior_EZ<0)
			Prior_EZ=FDR;

		///heuristic, want to smaller starting point
		//Prior_EZ*=0.9;
		///////////////build newBG for iterations////////////////////
		PWM single_logprob_bg_matrix=null;
		BGModel motifBG=new BGModel();
		
		double seedscore=motif.scoreWeightMatrix(seedstring);
		Iterator<FastaLocation> iter=origFalocs.iterator();
		TreeMap<String, Double> bgstrSet=new TreeMap<String, Double>();
		String[] allsites=new String[origFalocs.size()];
		int last=-1;
		int sitecount=0;
		double[][] single_logprob_bg_=new double[motif.columns()][4];
//		LinearEngine BGsearch=new LinearEngine(2);
//		BGsearch.build_index("E:/eclipse/data/promoter_bg.fa");
		PWM bgmodel=null;
		if(background.kmerPWMBG!=null)
		bgmodel=background.kmerPWMBG.get(motif.Consensus(true));
		int shifthead=0,shifttail=0;
		if(bgmodel!=null)
		{
			shifthead=(bgmodel.columns()-motif.columns())/2;
			shifttail=(bgmodel.columns()-motif.columns())/2;
			
			if(shifthead<0)
			{
				bgmodel=null;
				background.EnablePWMBG=false;
			}
		}
//		while(iter.hasNext())
//		{
//			FastaLocation currloc=iter.next();
//
//			String site=SearchEngine.getSite(currloc.getMin()-(motif.columns()-seedlen)/2, motif.columns());
//			if(currloc.ReverseStrand)
//				site=common.getReverseCompletementString(site);
//			allsites[sitecount]=site;
//			sitecount++;
//		}
//		try {
//			PWM bgmodel=new PWM(allsites);
		if(bgmodel!=null)
		{
			for (int i = 0; i < motif.columns(); i++) {
				for (int j = 0; j < 4; j++) {
					single_logprob_bg_[i][j]=bgmodel.log_matrix[shifthead+i][j];
				}
				
			}
		}
		else
		{
///////////////build newBG for iterations////////////////////

			double[][] tempM=common.forwardFiltering(this.background, (motif.columns()-seedstring.length())/2, seedstring);
			for (int i = (motif.columns()-seedstring.length())/2+seedstring.length(); i < motif.columns(); i++) {
				for (int j = 0; j < 4; j++) {
					single_logprob_bg_[i][j]=Math.log(tempM[i- (motif.columns()-seedstring.length())/2-seedstring.length()][j]);
				}
			}
			tempM=common.backwardSmoothing(this.background, (motif.columns()-seedstring.length())/2, seedstring);
			for (int i = (motif.columns()-seedstring.length())/2-1; i >=0; i--) {
				for (int j = 0; j < 4; j++) {
					single_logprob_bg_[i][j]=Math.log(tempM[(motif.columns()-seedstring.length())/2-1-i][j]);
				}
			}
		}
//			common.print2DArray(bgmodel.m_matrix);
//		} catch (IllegalAlphabetException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		} catch (IllegalSymbolException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		}
		int bgorder=1;
		motifBG.order=background.order;
		motifBG.conditionProb=(HashMap<String, Double>)background.conditionProb.clone(); //.BuildModel(bgstrSet, bgorder);
		bgstrSet.clear();
		//Prior_EZ=(double)truepos/origFalocs.size();  //Math.min(0.5, Math.max(0, (double)seqcount/origFalocs.size()) );
		double prior_gamma=0;//(double)truepos/seqcount;//
		///////////////build newBG for iterations////////////////////
		PWM bestPWM=null;
		int itercount=0;

		double bgseed_ignore2=0;
		if(bgmodel!=null)			
		{
					bgmodel.head=motif.head+shifthead;
					bgmodel.tail=motif.tail+shifttail;
					bgseed_ignore2=bgmodel.scoreWeightMatrix(seedstring)-seedscore;
					
		}
		double bgseed_ignore3=this.background.Get_LOGPROB(seedstring)-seedscore;
		
		do
		{

			itercount++;
		int truecount=0;
		int overlapcount=0;
		
		
		motif.pos_en=true;
		motif.peakrank_en=true;
		motif.strand_en=true;
		
//			double [][] A=new double[motif.columns()][4];
//			double [][] B=new double[motif.columns()][4];
		String consensus_core=motif.Consensus(false).substring(motif.head,motif.head+motif.core_motiflen);
		double bgseed_ignore=motifBG.Get_LOGPROB(seedstring)-seedscore;
		
		if(bgmodel!=null)
		{
			bgmodel.head=motif.head+shifthead;
			bgmodel.tail=motif.tail+shifttail;

		}
		
		System.out.println(consensus_core+"\t"+String.valueOf(bestscore));
		bestscore=Double.NEGATIVE_INFINITY;
		double[] temp_prior=new double[num_priorbin];	
		double[] temp_peakrank=new double[num_priorbin];	
		double[] temp_strand=new double[2];	
		double lognullprior=Math.log(1.0/num_priorbin);
		
		LinkedList<FastaLocation> Falocs=origFalocs;
		motif.inst_coverage=Falocs.size()*Prior_EZ;
		if(motif.inst_coverage<1)
			break;
		System.out.println(Prior_EZ);
		
		if(bestPWM!=null)
		{
			PWM tempMotif=motif.Clone();
			double div=common.PWM_Divergence(bestPWM, tempMotif);
			if(div<0.001&&motif.core_motiflen>min_motiflen&&bestPWM.core_motiflen==tempMotif.core_motiflen)
				return motif;
			bestPWM=tempMotif;
		}
		else
		bestPWM=motif.Clone();
		
//		if(stateCodes.contains(consensus_core.hashCode())||motif.core_motiflen>maxmotiflen)
//			return motif;
//		else
//			stateCodes.add(consensus_core.hashCode());
			
	
		double [][] sumLogBG_matrix=new double[motif.columns()][4];

		double[][] count_matrix=new double[motif.columns()][4];
		double[][] bgprob=new double[motif.columns()][4];
		common.fill2DArray(count_matrix, common.DoubleMinNormal);
		String consensus=motif.Consensus(false);
		//double logN025=log025*(motif.head+motif.tail);
		int motiflen=consensus_core.length();
		int avergeSeqlen=(SearchEngine2.TotalLen/SearchEngine2.getSeqNum());
	//	double[] single_logprob_bg=new double [4];
		double llr=1;

		
		double overlap_logBG=0;
		double overlap_prob=0;
		double overlap_expLLR=0;
		String lastsite="";
		int matchsitecount_seq=0;
		int overlap_pos=-motiflen;
		String ACGT="ACGT";
		//get single sym bgprob
//		for (int i = 0; i < single_logprob_bg.length; i++) {
//			single_logprob_bg[i]=motifBG.Get_LOGPROB(ACGT.substring(i, i+1));
//		}
		//update the loglik matrix
		
			
			
				Iterator<FastaLocation> iter2=Falocs.iterator();
				int count=0;
				int lastseq=-1;
				double [][] logBG_matrix=new double[motif.columns()][4];
				double [][] sumexpLLR=new double[motif.columns()][4];
				if(OOPS)
				{
					common.fill2DArray(logBG_matrix,Double.MIN_VALUE);
					common.fill2DArray(sumexpLLR,0);
				}
				double sitesperSeq=0;
				while(iter2.hasNext())
				{
					FastaLocation currloc=iter2.next();
				
					String site="";
					//forward site
					if(currloc.ReverseStrand)
					{
						site=SearchEngine2.getSite(currloc.getSeqId(), currloc.getSeqPos()-(motif.columns()-seedstring.length())/2, motif.columns());
						site=common.getReverseCompletementString(site);
					}
					else
						site=SearchEngine2.getSite(currloc.getSeqId(),currloc.getSeqPos()-(motif.columns()-seedstring.length())/2, motif.columns());
					if(site.charAt(site.length()/2)=='X')
						continue;
					
					if(Character.isLowerCase(site.charAt(site.length()/2)))
						truecount++;
					if(site.length()!=motif.columns())
						continue;
//					if(site.equalsIgnoreCase(""))
//						continue;
					//assume only one N for line break
					StringBuffer sb=new StringBuffer(site);
					for (int i = 0; i < site.length(); i++) {
						if(site.charAt(i)=='N'&& i<=site.length()/2)
						{
							for (int j = 0; j < i; j++) {
								sb.setCharAt(j, 'N');
							}
							continue;
						}
						else if(site.charAt(i)=='N')
						{
							for (int j = i+1; j < site.length(); j++) {
								sb.setCharAt(j, 'N');
							}
							break;
						}
						
					}
					site=sb.toString();
					String motifsite=site.substring(motif.head,motif.head+motiflen);
					double logprob_theta=motif.scoreWeightMatrix(motifsite);
					double logprob_BG=0;
					if(background.EnablePWMBG&&background.kmerlen==seedstring.length()&&bgmodel!=null)
					{
//						if(itercount<6)
						logprob_BG=bgmodel.scoreWeightMatrix(motifsite)-bgseed_ignore2;
//						else
//							logprob_BG=this.background.Get_LOGPROB(motifsite)-bgseed_ignore3;
					}
					else
						logprob_BG=motifBG.Get_LOGPROB(motifsite)-bgseed_ignore;

					int posbin=(int)(num_priorbin*((currloc.getSeqPos()+motiflen/2)%currloc.getSeqLen()/(double)currloc.getSeqLen()));
					int rankbin=num_priorbin*currloc.getSeqId()/SearchEngine2.getSeqNum();
					double logprior=0;
					if(motif.pos_prior.size()!=0&&motif.pos_en)
						logprior=Math.log(motif.pos_prior.get(posbin)+common.DoubleMinNormal)-lognullprior;
					if(motif.strand_en)
					{
						if(currloc.ReverseStrand)
							logprior+=Math.log((1-motif.strand_plus_prior)*2);
						else
						logprior+=Math.log(motif.strand_plus_prior*2);
					}
					if(motif.peakrank_prior.size()!=0&&motif.peakrank_en)
						logprior+=Math.log(motif.peakrank_prior.get(rankbin)+common.DoubleMinNormal)-lognullprior;
					
					double logDnaseProb=0;
					
					
					double loglik=logprob_theta+logprior-logprob_BG+logDnaseProb+Math.log(Prior_EZ/(1-Prior_EZ));
					if(loglik>10)
						loglik=10;
					double prob_theta=Math.exp(loglik)/(Math.exp(loglik)+1);//Math.exp(currloc.Score);				
					if(Double.isNaN(prob_theta))
						prob_theta=1;
					llr=Math.exp(loglik);
				//	bgstrSet.put(site.substring(0,(motif.columns()-seedlen)/2)+site.substring((motif.columns()-seedlen)/2+seedlen), 1-prob_theta);
					temp_peakrank[rankbin]+=prob_theta;
					temp_prior[posbin]+=prob_theta;
					if(currloc.ReverseStrand)
					{
						temp_strand[1]+=prob_theta;   
					}
					else
					{
						temp_strand[0]+=prob_theta;
					}

					if(OOPS)
						loglik-=common.DoubleMinNormal*Math.abs(currloc.getSeqLen()/2-currloc.getSeqPos()-motiflen/2); //add small bias to center
					if(!OOPS)
						bestscore+=loglik;
					

					if(OOPS&&currloc.getSeqId()!=lastseq)
					{
						if(lastseq!=-1)
							for (int i = 0; i < motif.columns(); i++) {	
								double sumexpLLRallsymid=0;
								 for (int symid = 0; symid < 4; symid++)
								 {
									 sumexpLLRallsymid+=sumexpLLR[i][symid];
								 }
	                         for (int symid = 0; symid < 4; symid++) {
	                        	 if(sitesperSeq<1)
	                        		 sitesperSeq=1;
	                        	// if(logBG_matrix[i][symid]!=Double.MIN_VALUE)
	                        	 {
								//	sumLogBG_matrix[i][symid]+=logBG_matrix[i][symid];

								//	count_matrix[i][symid]+=max_count_matrix[i][symid]/sitesperSeq;
									 prior_gamma=1-Math.pow(1-Prior_EZ, matchsitecount_seq);
									 double temp=sumexpLLR[i][symid]*Prior_EZ/((1-prior_gamma)+sumexpLLRallsymid*Prior_EZ);
						
									count_matrix[i][symid]+=temp;
									bgprob[i][symid]+=matchsitecount_seq-temp;
	                        	 }
							}

							}
						
						lastseq=currloc.getSeqId();
						sitesperSeq=0;
						lastsite="";
						overlap_pos=-motiflen;
						overlap_prob=0;
						matchsitecount_seq=0;
						common.fill2DArray(logBG_matrix,Double.MIN_VALUE);
						common.fill2DArray(sumexpLLR,0);
					}
					if(OOPS)
					{
						matchsitecount_seq++;
						 double expLLR=Math.exp(loglik-Math.log(Prior_EZ/(1-Prior_EZ)));
						if(Math.abs(currloc.getSeqPos()-overlap_pos)<maxmotiflen)
						{// overlap last site
							overlapcount++;
							if(expLLR>overlap_expLLR+common.DoubleMinNormal)
							{//the current site better than last site, remove the effect of last site
								for (int i = 0; i < site.length(); i++) {
									int symid=common.acgt[site.charAt(i)];
									if(symid<4)
									{
										sumexpLLR[i][symid]+=expLLR;

									}
									else
									{// new line seperator or N occurs, equal probability
										for (int j = 0; j < 4; j++) {
											sumexpLLR[i][j]+=expLLR*Math.exp(single_logprob_bg_[i][j]);
										
										}
									}
									if(lastsite.length()>0)
									{// remove the effect of last site
										int lastsymid=common.acgt[lastsite.charAt(i)];
										if(lastsymid<4)
										{
											sumexpLLR[i][lastsymid]-=overlap_expLLR;
											if(sumexpLLR[i][lastsymid]<0)
												sumexpLLR[i][lastsymid]=0; 
										}
										else
										{
											for (int j = 0; j < 4; j++) {
												sumexpLLR[i][j]-=overlap_expLLR*Math.exp(single_logprob_bg_[i][j]);
											}
										}
									}
								
								}	
								sitesperSeq+=prob_theta-overlap_prob;
								overlap_prob=prob_theta;
								overlap_expLLR=expLLR;
								overlap_logBG=logprob_BG;
								lastsite=site;
								overlap_pos=currloc.getSeqPos(); // update last site position
							}
							//else just ignore the low score site
						}								
						else
						{
							for (int i = 0; i < site.length(); i++) {
								int symid=common.acgt[site.charAt(i)];
								if(symid<4)
								{
									sumexpLLR[i][symid]+=expLLR;

								}
								else
								{
									for (int j = 0; j < 4; j++) {
										sumexpLLR[i][j]+=expLLR*Math.exp(single_logprob_bg_[i][j]);
							
									}
								}
							}
							sitesperSeq+=prob_theta;
							overlap_prob=prob_theta;
							overlap_expLLR=expLLR;
							overlap_logBG=logprob_BG;
							lastsite=site;
							overlap_pos=currloc.getSeqPos();
						}
						

						continue;
					}

					//////////////////////////////TCM/////////////////////////////
					for (int i = 0; i < site.length(); i++) {
						int symid=common.acgt[site.charAt(i)];
						if(symid>3)
						{
							for (int j = 0; j < 4; j++) {
								count_matrix[i][j]+=prob_theta*0.25;
							}
							
							continue; //meet new line separator
						}
						count_matrix[i][symid]+=prob_theta;

						if(consensus.charAt(i)!='N')
							continue;
						sumLogBG_matrix[i][symid]+=single_logprob_bg_[i][symid];//single_logprob_bg[symid];
						bgprob[i][symid]+=1-prob_theta;
					}
					//////////////////////////////TCM/////////////////////////////
			}
	
				

		
		//select the best column replacement
		double maxloglik=Double.NEGATIVE_INFINITY;
		
		int numbestSym=4;
		int bestCol=-1;
		ArrayList<Integer> bestSym=new ArrayList<Integer>(4);
		//use binomial p-value to decide stop extention
		RandomEngine rand=RandomEngine.makeDefault();
         // Binomial binomial=new new Binomial(R)
		
		
		int num_col_cand=0;
		//compute the optimal column assignments
		
		double [][] optimalcols=null;
		if(motif.core_motiflen==seedstring.length())
		{
//			OptExtParameter optimzer=new OptExtParameter(A, B);
//			optimzer.count_matrix=count_matrix;
//			optimzer.log_bg_matrix=single_logprob_bg_;
//			optimzer.llr=llr;
//			optimzer.find_parameter();
//			optimalcols=optimzer.optimalcol;
			double[] KL_Div=new double[motif.columns()];
			optimalcols=new double[motif.columns()][4];
			for (int k = 0; k < motif.columns(); k++)
			{
				if(consensus.charAt(k)!='N')
					continue;
				num_col_cand++;
				double[] p1=common.Normalize(count_matrix[k]);
				for (int u = 0; u < 4; u++) {
					if(p1[u]>0)
					{//KL_Div[k]+=Math.max(KL_Div[k], (p1[u]/Math.exp(single_logprob_bg_[k][u])));//
						KL_Div[k]+=p1[u]*(Math.log(p1[u])-single_logprob_bg_[k][u]);//Math.pow((p1[u]-Math.exp(single_logprob_bg_[k][u])),2);//
					}
					optimalcols[k][u]=p1[u]/Math.exp(single_logprob_bg_[k][u]);
				}
				optimalcols[k]=common.Normalize(optimalcols[k]);
			}
			double maxlh=Double.NEGATIVE_INFINITY;
			for (int k = 0; k < motif.columns(); k++) {
				int extralen=0;
				if(k<motif.head)
					extralen=motif.head-k;
				if(k>motif.columns()-motif.tail-1)
					extralen=k-(motif.columns()-motif.tail-1);
				if((extralen+motif.core_motiflen)>maxmotiflen)
					continue;
				if(maxlh<KL_Div[k]&&KL_Div[k]!=0)
				{
					maxlh=KL_Div[k];
					bestCol=k;
					
				}
				
			}
			
//			common.print2DArray(optimalcols); 
//			common.print2DArray(bgmodel.m_matrix);
//			common.print2DArray(count_matrix);
		}
		else
		{
			optimalcols=new double[count_matrix.length][4];
			for (int i = 0; i < optimalcols.length; i++) {
				double sumCount=0;
				for (int j = 0; j < 4; j++) {
					optimalcols[i][j]=count_matrix[i][j];///Math.exp(single_logprob_bg_[i][j])
					sumCount+=optimalcols[i][j];
				}
				for (int j = 0; j < 4; j++) {
					optimalcols[i][j]/=sumCount;
				}
			}
		}
		
		
//		common.print2DArray(count_matrix);
			if(bestCol==-1)
			{
			for (int i = 0; i < motif.columns(); i++) {
				if(consensus.charAt(i)!='N')
					continue;
				num_col_cand++;
				TreeMap<Double,Integer> orderSym=new TreeMap<Double,Integer>();
				double temploglik=0;
				double sumtemp=0;
				//sort the symid by count
				for (int j = 0; j < 4; j++) {
					double temp= count_matrix[i][j]/bgprob[i][j];
					orderSym.put(temp-common.DoubleMinNormal*j, j);
					sumtemp+=temp;
				}
				//System.out.println(sumtemp);
				//only consider best two sym
				int c=0;
				double sumcount=0;
				for(Double key: orderSym.descendingKeySet()) {
				
					if(c>=numbestSym)
						break;
					int symid=orderSym.get(key);
					//extra base, the only difference is  
					double logdiff=(Math.log(optimalcols[i][symid])-single_logprob_bg_[i][symid]);
					temploglik+=count_matrix[i][symid]*(logdiff); //key*(1+boundaryLoss);
//					temploglik=Math.max(temploglik, count_matrix[i][symid]/bgprob[i][symid]);
					
					sumcount+=count_matrix[i][symid];
					c++;			
				}

				//temploglik=temploglik-Falocs.size()*Math.log(0.25);//-common.lnEntropy(optimalcols[i])
				if(temploglik>maxloglik)
				{
					////////////////////make sure bestcol not exceed maxlen///////////////////////
					int extralen=0;
					if(i<motif.head)
						extralen=motif.head-i;
					if(i>motif.columns()-motif.tail-1)
						extralen=i-(motif.columns()-motif.tail-1);
					if((extralen+motif.core_motiflen)>maxmotiflen)
						continue;
					////////////////////make sure bestcol not exceed maxlen///////////////////////
					maxloglik=temploglik;
					bestCol=i;
					bestSym.clear();
					for(Double key: orderSym.descendingKeySet()) {
						int col=orderSym.get(key);
						if(bestSym.size()==numbestSym)
							break;
						bestSym.add(col);
						
					}										
				}
				
			}
			}
			
			
			//also maximize loglik for update other column

			for (int ecol = motif.head; ecol < motif.columns()-motif.tail; ecol++) 
			{
				//int ecol=iter1.next();
				motif.setWeights(ecol,common.Normalize(count_matrix[ecol]));
				
			}
			double max_sumCount=0;
			for (int j = 0; j < 4; j++) {
				max_sumCount+=count_matrix[count_matrix.length/2][j];
			}
			
			Prior_EZ=max_sumCount/Falocs.size();
			if(bestCol==-1)
				continue;




			double [] repColumnValue=new double[4];
			Arrays.fill(repColumnValue, common.DoubleMinNormal);
			

			


			//for (int i = 1; i <= bestSym.size(); i++)
			repColumnValue=optimalcols[bestCol];		
			
//			for (int i = 0; i < bgmodel.columns(); i++) {
//				bgmodel.setWeights(i, common.Normalize(bgprob[i]));
//			}

			
		//	motifBG.BuildModel(bgstrSet, bgorder);
			



			double boundaryLoss=0.25*Math.min(Math.abs(bestCol-motif.head), Math.abs(motiflen+motif.head-bestCol-1))/(double)avergeSeqlen;
			if((motif.head<bestCol&&bestCol<motiflen+motif.head)||OOPS)
				boundaryLoss=0;

			
//			if(max_sumNorm<=bestscore||max_sumCount<SearchEngine.getSeqNum()*min_support_ratio)
//				break;
//			else
				bestscore=maxloglik;
				motif.Score=bestscore;





			//determine whether pos_prior is significant needed
			double chistat=0;
			double sumPrior=0;
			for (int i = 0; i < temp_prior.length; i++) {
				sumPrior+=temp_prior[i];
			}
			//ignore the first and last bin
			double avgE=(sumPrior-temp_prior[0]-temp_prior[temp_prior.length-1])/(temp_prior.length-2);
				for (int j = 1; j < temp_prior.length-1; j++) {
					double temp=temp_prior[j]-avgE;
					chistat+=temp*temp/avgE;
				}
				double invFDR=ChiSquareDistQuick.inverseF(temp_prior.length-3, 1-FDR);
				
				if(chistat>invFDR)
				{
					if(!motif.pos_en)
					System.out.println("position bias on...");
					motif.pos_en=true;
				}
				motif.pos_prior.clear();
				for (int i = 0; i < temp_prior.length; i++) {
					temp_prior[i]/=sumPrior;
					motif.pos_prior.add(temp_prior[i]);
				}
				
			//determine whether strand_prior is significant needed
				double X=Math.max(temp_strand[0], temp_strand[1])+1;
				Binomial binomial=new Binomial((int)(temp_strand[0]+temp_strand[1])+2,0.5,rand);
				double pvalue_strand=1.0-binomial.cdf((int)X);
				if(pvalue_strand<this.FDR)
				{
					if(!motif.strand_en)
					System.out.println("strand bias on...");
					motif.strand_plus_prior=temp_strand[0]/(temp_strand[0]+temp_strand[1]);
					motif.strand_en=true;
				}
				
				//determine whether peakrank_prior is significant needed
				chistat=0;
				sumPrior=0;
				for (int i = 0; i < temp_peakrank.length; i++) {
					sumPrior+=temp_peakrank[i];
				}
				double avgcount=(sumPrior-temp_peakrank[temp_peakrank.length-1])/(temp_peakrank.length-1);
					for (int j = 0; j < temp_peakrank.length-1; j++) {
						double temp=temp_peakrank[j]-avgcount;
						chistat+=temp*temp/avgcount;
					}
					invFDR=ChiSquareDistQuick.inverseF(temp_peakrank.length-2, 1-FDR);
					if(chistat>invFDR)
					{
						if(!motif.peakrank_en)
						System.out.println("rank bias on...");
						motif.peakrank_en=true;
					}
					motif.peakrank_prior.clear();
					for (int i = 0; i < temp_peakrank.length; i++) {
						temp_peakrank[i]/=sumPrior;
						motif.peakrank_prior.add(temp_peakrank[i]);
					}
			

			
			motif.Prior_EZ=Prior_EZ;
			
			motif.inst_coverage=Prior_EZ*Falocs.size();
			//update motif column value, if significant
			chistat=0;
			double total=0;
			for (int i = 0; i < 4; i++) {
				total+=count_matrix[bestCol][i];
			}
			double rescale=1.5;
			for (int j = 0; j < 4; j++) {
					double Ei=total*Math.exp(single_logprob_bg_[bestCol][j])*rescale;
					double temp=count_matrix[bestCol][j]*rescale-Ei;
					chistat+=temp*temp/Ei;
				}
		
			invFDR=ChiSquareDistQuick.inverseF(3, 1-FDR);//(motif.columns()-seedstring.length())
			int extralen=0;
			if(bestCol<motif.head)
				extralen=motif.head-bestCol;
			if(bestCol>motif.columns()-motif.tail-1)
				extralen=bestCol-(motif.columns()-motif.tail-1);
			
			if(extralen+motif.core_motiflen>maxmotiflen)
				continue;
			if(chistat>invFDR)//num_col_cand*  pvalue<this.FDR
			{
			   motif.setWeights(bestCol, repColumnValue);
			   extendedCols.add(bestCol);
			}

			//when sample size is small, then ostrich policy let it extend
			if(total<10||motif.core_motiflen<=minmotiflen||extralen==0)
			{

				if(extralen+motif.core_motiflen<=(motif.columns()+seedstring.length())/2)
				{
					   motif.setWeights(bestCol, repColumnValue);
					   extendedCols.add(bestCol);
				}
			}


			
			
		}while(maxiter>itercount);
		return bestPWM;
	
		//return motif;
	}
	
	
	
	
	double sumLLR(PWM motif)
	{
		String consensus_core=motif.Consensus(true);
		int num_priorbin=SearchEngine2.getTotalLength()/SearchEngine2.getSeqNum()/this.resolution;
		int motiflen=consensus_core.length(); 
		
		//cant pass the FDR test
		if(motiflen==seedlen)
			return -10000;
		
		double[] temp_prior=new double[num_priorbin];
		double[] temp_peakrank=new double[num_priorbin];	
		double[] temp_strand=new double[2];	
		double overlap_loglik=0;
		double overlap_prob=0;
		double sitesperSeq=0;
		int overlap_pos=-motiflen;
		double lognullprior=Math.log(1.0/num_priorbin);
		double score=0;
		double Prior_EZ= (1-motif.inst_FDR)*SearchEngine2.getSeqNum()/SearchEngine2.TotalLen;
		
		double log_thresh=Double.NEGATIVE_INFINITY; //Math.log(1-Prior_EZ)-Math.log(Prior_EZ);
		//SearchEngine2.EnableBackground(background);
		log_thresh=motif.getThresh(0.99, SearchEngine2.getSeqNum()/SearchEngine2.TotalLen, background);
		SearchThread.bestonly=true;
		LinkedList<FastaLocation> Falocs=SearchEngine2.searchPattern(motif, log_thresh);
		SearchThread.bestonly=false;
		Iterator<FastaLocation> iter2=Falocs.iterator();
		int count=0;
		double [][]m_matrix=new double [motiflen][4];
		double [][] max_count_matrix=new double[motif.columns()][4];
		double match_seqCount=0;
		int lastseq=-1;
		double max_seqloglik=0;
		
	
		while(iter2.hasNext())
		{
			FastaLocation currloc=iter2.next();
//			double X_2=2*currloc.Score;
//			double v=(motif.core_motiflen-background.order)*3;
//			double z_test=Math.pow(X_2/v, 1.0/3)-(1-2/(9*v));
//			if(z_test>0)
				score+=1+Math.exp(currloc.Score);
			
//			String site=SearchEngine2.getSite(currloc.getSeqId(), currloc.getSeqPos(), currloc.getSeqLen());
//			double logprob_theta=currloc.Score-background.Get_LOGPROB(site);//  //
//			//	double logprob_theta=motif.scoreWeightMatrix(site)-SearchEngine2.BGscoreMap.get(motif.core_motiflen).get(currloc.getSeqId()).get(currloc.getSeqPos());
//
//				
//				//double logprob_BG=background.Get_LOGPROB(site.substring((site.length()-motiflen)/2, motiflen));
//				double logprior=0;
//				int prior_bin=(int)(num_priorbin*((currloc.getSeqPos()+motiflen/2)%currloc.getSeqLen()/(double)currloc.getSeqLen()));
//				int rankbin=num_priorbin*currloc.getSeqId()/SearchEngine.getSeqNum();
//				if(motif.pos_prior.size()!=0&&motif.pos_en)
//					logprior=Math.log(motif.pos_prior.get(prior_bin)+common.DoubleMinNormal)-lognullprior;
//				if(motif.strand_en)
//				{
//					if(currloc.ReverseStrand)
//						logprior+=Math.log((1-motif.strand_plus_prior)*2);
//					else
//					logprior+=Math.log(motif.strand_plus_prior*2);
//				}
//				if(motif.peakrank_prior.size()!=0&&motif.peakrank_en)
//					logprior+=Math.log(motif.peakrank_prior.get(rankbin)+common.DoubleMinNormal)-lognullprior;
//				double logDnaseProb=0;
//				double loglik=logprob_theta+logprior+logDnaseProb+Math.log(Prior_EZ/(1-Prior_EZ));
//				if(loglik>10)
//					loglik=10;
//				double prob_theta=Math.exp(loglik)/(Math.exp(loglik)+1);//Math.exp(currloc.Score);
//
//				if(OOPS)
//					loglik-=common.DoubleMinNormal*Math.abs(currloc.getSeqLen()/2-currloc.getSeqPos()-motiflen/2); //add small bias to center
//				if(!OOPS)
//					score+=loglik;
//				if(OOPS&&currloc.getSeqId()!=lastseq)
//				{
//					if(sitesperSeq<1)
//               		 sitesperSeq=1;
//					if(sitesperSeq>1)
//						match_seqCount+=1;
//					else
//						match_seqCount+=sitesperSeq;
//					//NegBinFunction.plogis(max_seqloglik);
//					lastseq=currloc.getSeqId();
//			
//					score+= max_seqloglik;
//
//
//					max_seqloglik=0;	
//					sitesperSeq=0;
//					overlap_prob=0;
//					overlap_pos=-motiflen;
//					common.fill2DArray(max_count_matrix,0);
//				}
//				
//				if(OOPS)
//				{
//					if(Math.abs(currloc.getSeqPos()-overlap_pos)<motiflen)
//					{
//						if(prob_theta>overlap_prob)
//						{
//							max_seqloglik+=prob_theta*loglik-overlap_prob*overlap_loglik;	
//							sitesperSeq+=prob_theta-overlap_prob;
//							overlap_prob=prob_theta;
//							overlap_loglik=loglik;						
//							overlap_pos=currloc.getSeqPos();
//						}
//						//else just ignore the low score site
//					}								
//					else
//					{
//						max_seqloglik+=prob_theta*loglik;
//
//						sitesperSeq+=1;
//						overlap_prob=prob_theta;
//						overlap_loglik=loglik;
//						
//						overlap_pos=currloc.getSeqPos();
//					}
//					
//			
//					
//
//				}

			
		}
		
		//last instance for OOPS
//		if(OOPS)
//		{
//			if(sitesperSeq<1)
//         		 sitesperSeq=1;
//			if(sitesperSeq>1)
//				match_seqCount+=1;
//			else
//				match_seqCount+=sitesperSeq;
//
//				score+= max_seqloglik/sitesperSeq;
//		}
		
		
		if(Double.isNaN(score))
			score=0;
		return score;
	}
	
	//for fix PWM score and BG score, but there is position prior
	double sumLLR(LinkedList<FastaLocation> LocList,double pwmloglik,double bgloglik)
	{
		double score=0;
		
       Iterator<FastaLocation> iter = LocList.iterator();
		int num_priorbin=SearchEngine2.getTotalLength()/SearchEngine2.getSeqNum()/this.resolution;
		double lognullprior=Math.log(1.0/num_priorbin);
		double max_seqloglik=0;
		int lastseq=-1;
		int lastpos=-seedlen;
		double[] loglik_seq=null;
		if(OOPS)
		{
			loglik_seq=new double[SearchEngine2.getSeqNum()+1];
		}
		while(iter.hasNext())
		{
			FastaLocation currloc=iter.next();
			double logprior=0;
			if(pos_prior.size()!=0)
				logprior=pos_prior.get( pos_prior.size()*currloc.getSeqPos()/currloc.getSeqLen())-lognullprior;
			
			double loglik=pwmloglik+logprior-bgloglik;
			
			if(OOPS)
			loglik-=common.DoubleMinNormal*Math.abs(currloc.getSeqLen()/2-currloc.getSeqPos()-seedlen/2); //add small bias to center
			

			
			if(OOPS)
			{
				int seqid=currloc.getSeqId();
				loglik_seq[seqid]=Math.max(loglik, loglik_seq[seqid]);
				continue;
			}
			
			score+=loglik;
		}
		if(OOPS)
		{
			for (int i = 0; i < loglik_seq.length; i++) {
				score+=loglik_seq[i];
			}
		}
		

		return score;
	}
	
	//for fix PWM score but there is 'N' inside the pattern
	double sumLLR(LinkedList<FastaLocation> LocList,double pwmloglik,int motiflen)
	{
		double score=0;
		int count=0;
		int num_priorbin=SearchEngine2.getTotalLength()/SearchEngine2.getSeqNum()/this.resolution;
		double lognullprior=Math.log(1.0/num_priorbin);
		Iterator<FastaLocation> iter2=LocList.iterator();
		while(iter2.hasNext())
		{
			FastaLocation currloc=iter2.next();
			String site="";
			//forward site
			site=SearchEngine2.getSite(currloc.getSeqId(),currloc.getSeqPos(),motiflen);
			if(site.equalsIgnoreCase(""))
				continue;
			if(site.indexOf("N")!=-1)
				continue;
			double logprob_BG=background.Get_LOGPROB(site);
			if(count>=SearchEngine2.forwardCount)
			{
				//reverse site
				site=common.getReverseCompletementString(site);
				
			}	
			double logprior=0;
			if(pos_prior.size()!=0)
				logprior=pos_prior.get( pos_prior.size()*currloc.getSeqPos()/currloc.getSeqLen())-lognullprior;
			double loglik=pwmloglik+logprior-logprob_BG;
			score+=loglik;

	}
		
		return score;
		
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// read parameters from command line
		Options options = new Options();
		options.addOption("i", true, "input fasta file");
		options.addOption("c", true, "control fasta file");
		options.addOption("bgmodel", true, "background model file");
		options.addOption("dnase", true, "dnase data file");
		options.addOption("prefix", true, "output directory");
		options.addOption("seedlen", true, "kmer seed motif length (default 5)");
		options.addOption("ratio",true, "sampling ratio (default 0.8)");
		options.addOption("n",true, "number of motifs in final report (default 5)");
		options.addOption("rs",true, "counting resolution (default 40 bp)");
		options.addOption("supp",true, "minimum support ratio, the percentage of peaks contains motif (default 0.05)");
		options.addOption("maxlen",true,"maxmimum length of the motif (default 30)");
		options.addOption("minlen",true,"minimum length of the motif (default 7)");
		options.addOption("maxw",true, "maximum size of motif binding region (default 600bp)");
		options.addOption("mask",false,"whether marking the top motif location in order to find co-motif");
		options.addOption("oops",false,"whether assuming only one occurrence per sequence");
		options.addOption("FDR",true,"fasle positive rate (default 0.01)");
		options.addOption("clust",true,"linkage type of hierachical clustering:"+Arrays.toString(LinkageCriterion.values()) );
		
		CommandLineParser parser = new GnuParser();
		Pomoda motifFinder=new Pomoda();
		try {
			CommandLine cmd = parser.parse( options, args);
			if(cmd.hasOption("i"))
			{
				motifFinder.inputFasta=cmd.getOptionValue("i");
			}
			else
			{
				throw new ParseException("no input fasta file");
			}
			if(cmd.hasOption("c"))
			{
				motifFinder.ctrlFasta=cmd.getOptionValue("c");
			}
			if(cmd.hasOption("bgmodel"))
			{
				motifFinder.bgmodelFile=cmd.getOptionValue("bgmodel");
			}
			if(cmd.hasOption("dnase"))
			{
				String dnasefile=cmd.getOptionValue("dnase");
				motifFinder.DnaseLib=common.ReadDelimitedFile("\t", dnasefile);
				//the window size is determined by the sequence length and dnase array length
				
			}
			if(cmd.hasOption("prefix"))
			{
				motifFinder.outputPrefix=cmd.getOptionValue("prefix");
			}
			if(cmd.hasOption("seedlen"))
			{
				motifFinder.seedlen=Integer.parseInt(cmd.getOptionValue("seedlen"));
			}
			if(cmd.hasOption("ratio"))
			{
				motifFinder.sampling_ratio=Double.parseDouble( cmd.getOptionValue("ratio"));
			}
			if(cmd.hasOption("n"))
			{
				motifFinder.num_motif=Integer.parseInt(cmd.getOptionValue("n"));
			}
			if(cmd.hasOption("rs"))
			{
				motifFinder.resolution=Integer.parseInt(cmd.getOptionValue("rs"));
			}
			if(cmd.hasOption("supp"))
			{
				motifFinder.min_support_ratio=Double.parseDouble(cmd.getOptionValue("supp"));
			}
			if(cmd.hasOption("maxlen"))
			{
				motifFinder.max_motiflen=Integer.parseInt(cmd.getOptionValue("maxlen"))*2-motifFinder.seedlen;
			}
			if(cmd.hasOption("minlen"))
			{
				motifFinder.min_motiflen=Integer.parseInt( cmd.getOptionValue("minlen"));
			}
			if(cmd.hasOption("maxw"))
			{
				motifFinder.ending_windowsize=Integer.parseInt( cmd.getOptionValue("maxw"));
			}
			if(cmd.hasOption("mask"))
			{
				motifFinder.maskflag=true;
			}
			else
			{
				motifFinder.maskflag=false;
			}
			if(cmd.hasOption("oops"))
			{
				motifFinder.OOPS=true;
			}

			if(cmd.hasOption("FDR"))
			{
				motifFinder.FDR=Double.parseDouble(cmd.getOptionValue("FDR"));
			}
			if(cmd.hasOption("clust"))
			{
				motifFinder.linkage=(cmd.getOptionValue("clust"));
			}
	
			System.out.println(Arrays.toString(args) );
		} catch (ParseException e) {
			// TODO Auto-generated catch block
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp( "JPomoda", options );
			return;
		}
		
		
		long fullstart = System.currentTimeMillis();
		long start = System.currentTimeMillis();
		//initialize Pomoda 
		motifFinder.initialize();
		long end = System.currentTimeMillis();
		System.out.println("Initialization time was "+(end-start)/1000+" seconds.");

		 System.currentTimeMillis();
		//get seed motifs
		ArrayList<PWM>  seedPWMs=motifFinder.getSeedMotifs();
//		ArrayList<PWM>  seedPWMs2=motifFinder.getSeedMotifs3(10);
		//LinkedList<PWM>  seedPWMs=new LinkedList<PWM>();
		
		double topseed_Score=0;
		
		File file = new File(motifFinder.outputPrefix+"jpomoda_raw.pwm"); 
		try {
			
//			seedPWMs.clear();
//		//	seedPWMs.addAll(common.LoadPWMFromFile("D:\\eclipse\\data\\test.pwm").subList(0, 1));
//		//	double llrscore2=motifFinder.sumLLR(seedPWMs.get(0));
//			seedPWMs.add(new PWM(new String[]{"NNNNNNNNNNNNNNNNNNNNNNNNNCGAAANNNNNNNNNNNNNNNNNNNNNNNNN"}));
//			seedPWMs.add(new PWM(new String[]{"NNNNNNNNNNNNNNNNNNNNNNNNNCGCGCNNNNNNNNNNNNNNNNNNNNNNNNN"}));
//			seedPWMs.add(new PWM(new String[]{"NNNNNNNNNNNNNNNNNNNNNNNNNTCGCGNNNNNNNNNNNNNNNNNNNNNNNNN"}));
//			seedPWMs.add(new PWM(new String[]{"NNNNNNNNNNNNNNNNTCACANNNNNNNNNNNNNNNN"}));
//			seedPWMs.add(new PWM(new String[]{"NNNNNNNNNNNNNNNNNNNNNNNNNACTACNNNNNNNNNNNNNNNNNNNNNNNNN"}));
//			seedPWMs.add(new PWM(new String[]{"NNNNNNNNNNNNNNNNNNNNNNNNNGTTAANNNNNNNNNNNNNNNNNNNNNNNNN"}));
			BufferedWriter writer = new BufferedWriter(new FileWriter(file));
			TreeMap<Double, PWM> sortedPWMs=new TreeMap<Double, PWM>();
		//extend and refine motifs
		for (int i = 0; i < seedPWMs.size(); i++) {
			PWM motif=seedPWMs.get(i);
			motif.Score=1;
			try
			{
					int num_priorbin=motifFinder.SearchEngine2.getTotalLength()/motifFinder.SearchEngine2.getSeqNum()/motifFinder.resolution;
					
					if(motif.pos_prior.size()==0)
					{
						for (int i1 = 0; i1 <num_priorbin ; i1++) {
							motif.pos_prior.add(1.0/num_priorbin);
							motif.peakrank_prior.add(1.0/num_priorbin);
						}
					}
					
		
//			if(motifFinder.DnaseLib!=null)
//			{
//				motif.DnaseBG=motifFinder.dnaseBG;
//				motif.DnaseFG=new NegativeBinomialDist(motifFinder.dnaseBG.getGamma(), motifFinder.dnaseBG.getP());
//				motif.Dnase_prob=new ArrayList<Double>(motifFinder.DnaseWindow*2);
//				for (int j = 0; j < motifFinder.DnaseWindow*2; j++) {
//					motif.Dnase_prob.add( 1.0/(motifFinder.DnaseWindow*2));
//				}
//			}
			
			System.out.println("Extending...");
			seedPWMs.set(i,motifFinder.Column_Replacement_2(motif));
			if(seedPWMs.get(i)==null||seedPWMs.get(i).core_motiflen<motifFinder.min_motiflen)
			{
				seedPWMs.remove(i);
				continue;
			}
			System.out.println("Relaxing...");
			seedPWMs.set(i, motifFinder.Relax_Seed_3(seedPWMs.get(i)));
			if(seedPWMs.get(i)==null)
			{
				seedPWMs.remove(i);
				continue;
			}
			seedPWMs.get(i).Name="Motif"+String.valueOf(i+1);
			// to make different length comparable ,need to consider the instance coverage
//			seedPWMs.get(i).Score=seedPWMs.get(i).Score/seedPWMs.get(i).inst_coverage;//corrected score
			if(motifFinder.maskflag&&seedPWMs.get(i).peakrank_en&&seedPWMs.get(i).pos_en)// (seedPWMs.get(i).Prior_EZ*motifFinder.SearchEngine2.TotalLen)>(motifFinder.SearchEngine2.getSeqNum()))
			{
				//do something to mark the locations in SearchEngine
				
				PWM mainmotif=seedPWMs.get(i).trim();
				if(mainmotif==null)
					continue;
				
				if(mainmotif.core_motiflen>=motifFinder.min_motiflen)
				{
					System.out.println("Masking... "+mainmotif.Consensus(true));
						motifFinder.Masking(mainmotif);
				}
				topseed_Score=seedPWMs.get(i).Score;
				//motifFinder.SearchEngine2.EnableBackground(motifFinder.background);
			}		
			writer.write(seedPWMs.get(i).toString());
			}
			catch (Exception ex)
			{
				ex.printStackTrace();
				continue;
			}

		}
		end = System.currentTimeMillis();
		System.out.println("Find "+seedPWMs.size()+" motifs time was "+(end-start)/1000+" seconds.");
		
		writer.close();
		start = System.currentTimeMillis();
		//restore the unmask fasta
		if(motifFinder.maskflag)
			motifFinder.SearchEngine2.build_index(motifFinder.inputFasta);
		PWMevaluator evaluator=new PWMevaluator(motifFinder);
		//using full LLR score
//		motifFinder.SearchEngine2.EnableBackground(motifFinder.background);
		
		 ArrayList<AUCComputeThread> threadpool=new ArrayList<AUCComputeThread>(seedPWMs.size());
//			ExecutorService executor = Executors.newFixedThreadPool(Math.min(motifFinder.max_threadNum,seedPWMs.size()));
			
			if(motifFinder.ctrlFasta!="")
			{
				motifFinder.BGSearch=new LinearEngine(4);
				motifFinder.BGSearch.build_index(motifFinder.ctrlFasta);
				System.gc();
			}
///////////////////////////////////////evaluate different motif in parallel///////////////////////////////////////////		 
		for (int i = 0; i < seedPWMs.size(); i++) {
			AUCComputeThread t1=new AUCComputeThread(evaluator, seedPWMs.get(i).trim(), motifFinder.BGSearch);
//			if( motifFinder.BGSearch!=null)
			if(motifFinder.maskflag)
				t1.ZscoreFlag=true;
				t1.run();
//			else
//				executor.execute(t1);
			threadpool.add(t1);
		}
		
		// Wait until all threads are finish
//		while (!executor.isTerminated()) {
//			Thread.sleep(3000);
//			executor.shutdown();
//		}
		
		  for (int i = 0; i < threadpool.size(); i++) {
				double llrscore=threadpool.get(i).getResult();
				System.out.println(seedPWMs.get(i).Consensus(true)+" AUC:"+ llrscore);
				seedPWMs.get(i).Score=llrscore+seedPWMs.get(i).Prior_EZ/10;
				if(llrscore>0.5)
				sortedPWMs.put(llrscore+seedPWMs.get(i).Prior_EZ/10, seedPWMs.get(i));
		  }
///////////////////////////////////////evaluate different motif in parallel///////////////////////////////////////////			
		seedPWMs.clear();
		for (Double key : sortedPWMs.descendingKeySet()) {
			seedPWMs.add(sortedPWMs.get(key));
		}
		
//		motifFinder.SearchEngine2.DisableBackground();
		end = System.currentTimeMillis();
		
		System.out.println("Evaluation time was "+(end-start)/1000+" seconds.");
		evaluator.SearchEngine.DisableBackground();
		file = new File(motifFinder.outputPrefix+"jpomoda_clust.pwm"); 
		writer = new BufferedWriter(new FileWriter(file));
		//clustering motif, re-initialize
		start = System.currentTimeMillis();
		System.out.println("Clustering...");
		sortedPWMs.clear();
		PWMcluster clustering=new PWMcluster(motifFinder);
	
			ArrayList<PWM>  clusterPWMs=clustering.Clustering_Fast(seedPWMs,motifFinder.num_motif);
	
			for(PWM pwm:clusterPWMs)
			{
				sortedPWMs.put(pwm.Score, pwm); //desc order
			}
			System.out.println("Clustered Motifs:");
			int c=0;
//			GapImprover gimprover=new GapImprover(motifFinder);
			
			for(Double key:sortedPWMs.descendingKeySet())
			{
				sortedPWMs.get(key).Name="Motif_clust"+String.valueOf(c+1);
				c++;
				System.out.println(sortedPWMs.get(key).Consensus(true)+'\t'+sortedPWMs.get(key).Score);
				writer.write(sortedPWMs.get(key).toString());
				//if(sortedPWMs.get(key).pos_en)
				motifFinder.DrawDistribution(sortedPWMs.get(key).pos_prior,motifFinder.outputPrefix+"/"+sortedPWMs.get(key).Name+"_posdist.png");
				//if(sortedPWMs.get(key).peakrank_en)
				motifFinder.DrawDistribution(sortedPWMs.get(key).peakrank_prior,motifFinder.outputPrefix+"/"+sortedPWMs.get(key).Name+"_rankdist.png");
	
				if(motifFinder.DnaseLib==null)
					System.out.println("PWM AUC:"+ evaluator.calcAUC(sortedPWMs.get(key),null));
				else
					System.out.println("PWM AUC:"+ evaluator.calcAUC(sortedPWMs.get(key),motifFinder.DnaseLib,null));

			}
			
			writer.close();
			end = System.currentTimeMillis();
			System.out.println("Clustering time was "+(end-start)/1000+" seconds.");
			
			long fullend = System.currentTimeMillis();
			
			System.out.println("Total time was "+(fullend-fullstart)/1000+" seconds.");
		} 
		catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			
			
//		} catch (IllegalAlphabetException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		} catch (IllegalSymbolException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
			
			
		} 
		}
		

		

	

}
