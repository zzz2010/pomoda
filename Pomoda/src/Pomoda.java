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
import java.util.TreeMap;
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

import umontreal.iro.lecuyer.probdist.NegativeBinomialDist;
import umontreal.iro.lecuyer.probdist.NormalDist;
import umontreal.iro.lecuyer.probdistmulti.DirichletDist;
import umontreal.iro.lecuyer.util.Num;

import cern.jet.random.Binomial;
import cern.jet.random.engine.RandomEngine;





public class Pomoda {

	public String outputPrefix="./";
	public String inputFasta;
	public int max_iterNum=50;
	public String ctrlFasta="";
	public String bgmodelFile="";
	public int seedlen=5;
	public boolean debug=false;
	public boolean OOPS=false;
	public int resolution=10;
	public int starting_windowsize=200;
	public int ending_windowsize=600;
	public double FDR=0.001;
	public int max_motiflen=31;
	public int num_motif=5;
	public double sampling_ratio=1;
	public double min_support_ratio=0.1;
	public boolean maskflag=true;
	public HashEngine SearchEngine;
	public LinearEngine SearchEngine2;
	public BGModel background;
	
	public String linkage="WPGMA";
	public ArrayList<Double> pos_prior;
	public NegativeBinomialDist dnaseBG=null;
	public ArrayList<Double[]> DnaseLib=null;
	public int DnaseWindow=1;
	
	public void initialize()
	{
		//build hash index
		SearchEngine=new HashEngine(5);
		SearchEngine.build_index(this.inputFasta);
		SearchEngine2=new LinearEngine(6);
		SearchEngine2.build_index(this.inputFasta);
//		if(SearchEngine_Test())
//			System.out.println("SearchEngine_Test : pass");
			
		
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
//		if(BGModel_Test())
//			System.out.println("BGModel_Test : pass");
		pos_prior=new ArrayList<Double>(SearchEngine.getTotalLength()/SearchEngine.getSeqNum()/this.resolution);
		SearchEngine2.EnableBackground(background);
		
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
	
	
	private boolean PWM_Test()
	{
		boolean pass=true;
		try {
			//"ATTAAA","AATTAA","AAATTA","AAAATT","AAAAAT"
			PWM testpwm=new PWM(new String[]{"ATTAAA","AATTAA","AAATTA","AAAATT","TAAAAT"});
			BGModel bgtest=new BGModel();
			bgtest.BuildModel(new String[]{"ATTATT","ATTATT","ATTATT","ATTATT","ATTATT","ATTATT","ATTATT","ATTATT","ATTATT","ATTATT","ATTATT"}, 4);
			String testpatt="ATTATT";
			String cns=testpwm.Consensus(true);
			testpwm.print();
			double score1=Math.exp(testpwm.scoreWeightMatrix(testpatt));
			double score2=Math.exp(bgtest.Get_LOGPROB(testpatt));
			if(score2<score1)
			{
				pass=false;
			}
			LinkedList<PWM> AR=new LinkedList<PWM>();
			
			AR.add(new PWM(new String[]{"NNNNNNNNNNGTACANNNNNNNNNN"}));
			AR.add(new PWM(new String[]{"NNNNNNNNNGNACANNNNNNNNN"}));
			AR.add(new PWM(new String[]{"NNNNNNNGNACANNNNGNNNNNNN"}));
			AR.add(new PWM(new String[]{"NNNNNNGNACANNNNGTNNNNNN"}));
			AR.add(new PWM(new String[]{"NNNNNNGNACANNNNGTNCNNNNNN"}));
			AR.add(new PWM(new String[]{"NNNNNGNACANNNTGTNCNNNNN"}));
			
			for (int i = 0; i < AR.size(); i++) {
			//	AR.get(i).print();
			     this.Column_Replacement_(AR.get(i));
			     this.Relax_Seed_(AR.get(i));
			}
				
			
		} catch (IllegalAlphabetException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IllegalSymbolException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
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
	
	private boolean SearchEngine_Test()
	{
		boolean pass=true;
		//Exact Test
		String testpattern="ACGTAC";
		String Ctestpattern="GTACGT";
		LinkedList<Integer> poslist=SearchEngine.searchPattern(testpattern, 0);
		Iterator<Integer> iter=poslist.iterator();
		while(iter.hasNext())
		{
			String sitestring=SearchEngine.getSite(iter.next(), 6);
			if(sitestring.equalsIgnoreCase(testpattern)||sitestring.equalsIgnoreCase(Ctestpattern))
				continue;
			else
			{
				pass=false;
				System.out.println("Exact Match Test Fail...");
				System.out.println(sitestring);
				break;
			}
			
		}
		//Mismatch Test
		testpattern="ACGNTAC";
		Ctestpattern="GTANCGT";
		poslist=SearchEngine.searchPattern(testpattern, 0);
		 iter=poslist.iterator();
		while(iter.hasNext())
		{
			String sitestring=SearchEngine.getSite(iter.next(), 7);
			if(sitestring.matches(testpattern.replace("N", "\\w"))||sitestring.matches(Ctestpattern.replace("N", "\\w")))
				continue;
			else
			{
				pass=false;
				System.out.println("Mismatch Match Test Fail...");
				System.out.println(sitestring);
				break;
			}
			
		}
		
		return pass;
	}
	
	public void Masking(PWM motif)
	{
		
		double log_thresh=motif.getThresh(sampling_ratio, FDR, background);
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
		StringBuffer sb=new StringBuffer(SearchEngine.CharText);
		String X_Str="";
		for (int i = 0; i < motiflen; i++) {
			X_Str+="X";
		}
		while(iter2.hasNext())
		{
			FastaLocation currloc=iter2.next();
			String Site1=SearchEngine2.getSite(currloc.getSeqId(), currloc.getSeqPos(), motiflen);
			String Site2=SearchEngine.getSite(currloc.getMin()+currloc.getSeqId(),motiflen);
			String rep=SearchEngine2.ForwardStrand.get(currloc.getSeqId()).replace(Site1, X_Str);
			SearchEngine2.ForwardStrand.set(currloc.getSeqId(), rep);
		//	rep=SearchEngine2.ReverseStrand.get(currloc.getSeqId()).replace(common.getReverseCompletementString(Site1) , X_Str);
		//	SearchEngine2.ReverseStrand.set(currloc.getSeqId(), rep);
			//hash engine...
			sb.replace(currloc.getMin()+currloc.getSeqId(), currloc.getMin()+currloc.getSeqId()+motiflen, X_Str);
		}
		SearchEngine.CharText=sb.toString();
		
	}
	
	
	public ArrayList<PWM>	 getSeedMotifs() {
		int max_num_Seeds=num_motif*num_motif;
		ArrayList<PWM>	SeedMotifs=new	ArrayList<PWM>(max_num_Seeds);
		int LIBSIZE=SearchEngine.getTotalLength();
		String ACGT="ACGT";
		int loopnum=1<<(2*seedlen);
		int maxRange=1<<(2*seedlen);
		 ValueComparator bvc =  new ValueComparator();
		TreeMap<Integer,Double> seedScores=new TreeMap<Integer,Double>();
		for (int i = 0; i < loopnum; i++) {
			String pattern="";
			int hash=i;
			//ignore the reverseComplement
			if(common.getReverseComplementHashing(hash,seedlen)<i)
			{
				continue;
			}
			else
			pattern=common.Hash2ACGT(hash,seedlen);
			
			LinkedList<Integer> positionlist=SearchEngine.searchPattern(pattern,0);
			LinkedList<FastaLocation> LocList=SearchEngine.Int2Location(positionlist);
			double logprob_bg=Math.max(background.Get_LOGPROB(pattern), background.Get_LOGPROB(common.getReverseCompletementString(pattern) )) ;

			double score=0;
			if(pos_prior.size()==0&&!OOPS)
				score=positionlist.size()*(-common.DoubleMinNormal*seedlen-logprob_bg);//sum loglik ,-0.037267253272904234 is from pseudo count
			else
				score=sumLLR(LocList,-common.DoubleMinNormal*seedlen,logprob_bg);
			seedScores.put(hash, score);
		}
		//sort by score
		List<Map.Entry<Integer,Double>> mappingList=new ArrayList<Map.Entry<Integer,Double>>(seedScores.entrySet());
		Collections.sort(mappingList, bvc);
		
		//just take top ones
		    for (Map.Entry<Integer,Double> pair : mappingList) {
		    	if(SeedMotifs.size()==max_num_Seeds)
		    		break;
			    String toppattern=common.Hash2ACGT(pair.getKey(),seedlen);
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

	public PWM Relax_Seed_2(PWM motif)
	{
		double log025=Math.log(0.25);
		//relax the conserved column 
		String consensus=motif.Consensus(false);
		double main_prop=Math.pow(sampling_ratio*9/(seedlen*seedlen-2*seedlen+4), 1.0/(seedlen-2)); //2 mismatch
			
			//Math.pow(3*sampling_ratio/(seedlen-2), 1.0/(seedlen-1));// 1 mismatch
			
			//Math.pow(sampling_ratio*9/(seedlen*seedlen-2*seedlen+4), 1.0/(seedlen-2)); //2 mismatch
		for (int i = 0; i < consensus.length(); i++) {
			if(consensus.charAt(i)=='N')
				continue;
			boolean conserved=false;
			for (int j = 0; j < 4; j++) {
				 if(motif.m_matrix[i][j]>(1-common.DoubleMinNormal))
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
		double Prior_EZ=(double)SearchEngine.getSeqNum()/(SearchEngine.getTotalLength()/motif.core_motiflen);
		double bestscore=motif.Score;
		double lastscore=0;
		int num_priorbin=SearchEngine.getTotalLength()/SearchEngine.getSeqNum()/this.resolution;
		
		if(motif.pos_prior.size()==0)
		{
			for (int i = 0; i <num_priorbin ; i++) {
				motif.pos_prior.add(1.0/num_priorbin);
			}
		}
		
		
		HashSet<Integer> stateCodes=new HashSet<Integer>();
	
		int flankingLen=0;
		if(motif.head+flankingLen+motif.core_motiflen>=motif.columns()||motif.head-flankingLen<0)
			flankingLen=0;
		int iter_count=0;
		PWM bestPWM=motif.Clone();
		double sitesperSeq=0;
		double log_thresh=Math.log(1-Prior_EZ)-Math.log(Prior_EZ);
		LinkedList<FastaLocation> Falocs=SearchEngine2.searchPattern(motif, log_thresh);
		
		do
		{

			iter_count++;
			String consensus_core=motif.Consensus(true);
			
			int motiflen=consensus_core.length(); 
			System.out.println(consensus_core+"\t"+String.valueOf(bestscore));
			bestscore=0;
			double[] temp_prior=new double[num_priorbin];

			double overlap_loglik=0;
			double overlap_prob=0;
			
			int overlap_pos=-motiflen;
			double lognullprior=Math.log(1.0/num_priorbin);
	
			
			//double log_thresh=motif.getThresh(sampling_ratio, 2*FDR, background)- motiflen*log025;
			//double log_thresh=Math.log(1-Prior_EZ)-Math.log(Prior_EZ);
			
			double [][]m_matrix=new double [motiflen+flankingLen*2][4];
			double [][] max_count_matrix=new double[motif.columns()][4];
			//update the loglik matrix
			//LinkedList<FastaLocation> Falocs=SearchEngine2.searchPattern(motif, log_thresh);
			System.out.println("number of occurrences: "+String.valueOf(Falocs.size()));

					Iterator<FastaLocation> iter2=Falocs.iterator();
					int count=0;
					double match_seqCount=0;
					int lastseq=-1;
					String lastsite="";
					double max_seqloglik=0;

					String max_seqsite="";
					while(iter2.hasNext())
					{
						FastaLocation currloc=iter2.next();
						String site="";
						//forward site
						site=SearchEngine2.getSite(currloc.getSeqId(), currloc.getSeqPos()-flankingLen,motiflen+flankingLen*2);
						//check masking
						if(site.indexOf('X')>-1)
							continue;
						
//						if(site.equalsIgnoreCase(""))
//							continue;
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

						if(site.length()!=motiflen+flankingLen*2)
							continue;
						
						//double logprob_theta=currloc.Score;//include the bg log_prob in the score
						double logprob_theta=motif.scoreWeightMatrix(site)-SearchEngine2.BGscoreMap.get(motif.core_motiflen).get(currloc.getSeqId()).get(currloc.getSeqPos());

						
						//double logprob_BG=background.Get_LOGPROB(site.substring((site.length()-motiflen)/2, motiflen));
						double logprior=0;
						int prior_bin=(int)(num_priorbin*((currloc.getSeqPos()+motiflen/2)%currloc.getSeqLen()/(double)currloc.getSeqLen()));
						if(motif.pos_prior.size()!=0)
							logprior=Math.log(motif.pos_prior.get(prior_bin))-lognullprior;
						double logDnaseProb=0;

						double loglik=logprob_theta+logprior+logDnaseProb+Math.log(Prior_EZ/(1-Prior_EZ));
						if(loglik>200)
							loglik=200;
						double prob_theta=Math.exp(loglik)/(Math.exp(loglik)+1);//Math.exp(currloc.Score);
						//double prob_theta_only=Math.exp(logprob_theta+logprior)/(Math.exp(logprob_theta+logprior)+1);
						if(Double.isNaN(prob_theta))
							prob_theta=1;//upper flow

						
						temp_prior[prior_bin]+=prob_theta;//make smaller
						if(OOPS)
							loglik-=common.DoubleMinNormal*Math.abs(currloc.getSeqLen()/2-currloc.getSeqPos()-motiflen/2); //add small bias to center
						if(!OOPS)
							bestscore+=loglik;
						if(OOPS&&currloc.getSeqId()!=lastseq)
						{
							if(sitesperSeq<1)
                       		 sitesperSeq=1;
							if(sitesperSeq>1)
								match_seqCount+=1;
							else
								match_seqCount+=sitesperSeq;
							//NegBinFunction.plogis(max_seqloglik);
							lastseq=currloc.getSeqId();
							
							//System.out.println(max_seqsite.toUpperCase());
							if(max_seqsite.length()==motiflen+flankingLen*2)
							{
								for (int i = 0; i < motiflen+flankingLen*2; i++) {
   								 for (int symid = 0; symid < 4; symid++) 
									{
									m_matrix[i][symid]+=max_count_matrix[i][symid]/sitesperSeq;
									}
								
								}
								
								
							}

							bestscore+= max_seqloglik/sitesperSeq;

							max_seqloglik=0;	
							sitesperSeq=0;
							overlap_pos=-motiflen;
							common.fill2DArray(max_count_matrix,0);
						}
						
						if(OOPS)
						{
							if(Math.abs(currloc.getSeqPos()-overlap_pos)<motiflen)
							{
								if(prob_theta>overlap_prob)
								{
									max_seqloglik+=prob_theta*loglik-overlap_prob*overlap_loglik;
									for (int i = 0; i < site.length(); i++) {
										int symid=common.acgt(site.charAt(i));
										if(symid<4)
											max_count_matrix[i][symid]+=prob_theta;
									if(lastsite.length()>0)
									{
										int lastsymid=common.acgt(lastsite.charAt(i));
										if(lastsymid<4)
											max_count_matrix[i][lastsymid]-=overlap_prob;
										else
											for (int j = 0; j < 4; j++) {
												max_count_matrix[i][j]-=0.25*overlap_prob;
											}
									}
									
									}	
									sitesperSeq+=prob_theta-overlap_prob;
									overlap_prob=prob_theta;
									overlap_loglik=loglik;
									lastsite=site;
								}
								//else just ignore the low score site
							}								
							else
							{
								max_seqloglik+=prob_theta*loglik;
								
								for (int i = 0; i < site.length(); i++) {
									int symid=common.acgt(site.charAt(i));

									if(symid>3)
									{
										for (int j = 0; j < 4; j++) {
											max_count_matrix[i][j]+=0.25*overlap_prob;
										}
										continue;
									}
								max_count_matrix[i][symid]+=prob_theta;
								
								}	
								sitesperSeq+=prob_theta;
								overlap_prob=prob_theta;
								overlap_loglik=loglik;
								lastsite=site;
							}
							overlap_pos=currloc.getSeqPos();

						}
						
						
						if(!OOPS)
						{
						for (int i = 0; i < motiflen+flankingLen*2; i++) {
							int symid=common.acgt(site.charAt(i));
							if(symid>3)
								continue;
							m_matrix[i][symid]+=prob_theta;//prob_theta;
						}
						

						}
					}
					
					
					//last instance for OOPS
					if(OOPS)
					{
						if(sitesperSeq<1)
                     		 sitesperSeq=1;
						if(sitesperSeq>1)
							match_seqCount+=1;
						else
							match_seqCount+=sitesperSeq;
						for (int i = 0; i < motiflen+flankingLen*2; i++) {
							for (int symid = 0; symid < 4; symid++) 
							{
							m_matrix[i][symid]+=max_count_matrix[i][symid]/sitesperSeq+common.DoubleMinNormal;// pesudo count
							}
						}

							bestscore+= max_seqloglik/sitesperSeq;
					}
					
					Prior_EZ=match_seqCount/Falocs.size();
//					if(lastscore>=bestscore||match_seqCount<min_support_ratio*SearchEngine2.getSeqNum())
//						break;
//					else
					{
						if(flankingLen==0)
						{
//							if(lastscore!=1)
//								motif.Score*=bestscore/lastscore;
//							else
						motif.Score=bestscore;//-Math.log(SearchEngine2.getSeqNum())*2*motif.core_motiflen;
						lastscore=bestscore;
						}
						
						if(Falocs.size()==0|| stateCodes.contains(consensus_core.hashCode()))
							return bestPWM;
						else
							stateCodes.add(consensus_core.hashCode());
					
						if(bestPWM.Score<motif.Score)
						  bestPWM=motif.Clone();
						for (int i = 0; i < m_matrix.length; i++) {
							motif.setWeights(i+motif.head-flankingLen,common.Normalize(m_matrix[i]));
						}
					
					
						motif.pos_prior.clear();
						temp_prior=common.Normalize(temp_prior);
						for (int i = 0; i < temp_prior.length; i++) {
							motif.pos_prior.add(temp_prior[i]);
						}
						

						
					}
				
					//not allow to grow in the iterations
					flankingLen=0;
					
					System.out.println("number of occurred sequences: "+String.valueOf(match_seqCount));
			}while(motif.core_motiflen<max_motiflen&&iter_count<=max_iterNum);
			
		if(DnaseLib!=null)
			DrawDistribution(motif.Dnase_prob,"Dnase_plot.png");
		return bestPWM;
	}
	
	public PWM Relax_Seed_(PWM motif)
	{
		double log025=Math.log(0.25);
		//relax the conserved column 
		String consensus=motif.Consensus(false);
		double main_prop=Math.pow(sampling_ratio*9/(seedlen*seedlen-2*seedlen+4), 1.0/(seedlen-2)); //2 mismatch
			
			//Math.pow(3*sampling_ratio/(seedlen-2), 1.0/(seedlen-1));// 1 mismatch
			
			//Math.pow(sampling_ratio*9/(seedlen*seedlen-2*seedlen+4), 1.0/(seedlen-2)); //2 mismatch
		for (int i = 0; i < consensus.length(); i++) {
			if(consensus.charAt(i)=='N')
				continue;
			boolean conserved=false;
			for (int j = 0; j < 4; j++) {
				 if(motif.m_matrix[i][j]>(1-common.DoubleMinNormal))
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
		double Prior_EZ=0.5;
		double bestscore=motif.Score;
		double lastscore=1;
		int num_priorbin=SearchEngine.getTotalLength()/SearchEngine.getSeqNum()/this.resolution;
		
		if(motif.pos_prior.size()==0)
		{
			for (int i = 0; i <num_priorbin ; i++) {
				motif.pos_prior.add(1.0/num_priorbin);
			}
		}
		
		
		HashSet<Integer> stateCodes=new HashSet<Integer>();
	
		int flankingLen=1;
		if(motif.head+flankingLen+motif.core_motiflen>=motif.columns()||motif.head-flankingLen<0)
			flankingLen=0;
		int iter_count=0;
		PWM bestPWM=motif.Clone();
		do
		{
			 try {
					System.setErr(new PrintStream(new FileOutputStream("system_err.txt")));
				} catch (FileNotFoundException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			iter_count++;
			String consensus_core=motif.Consensus(true);
			
			int motiflen=consensus_core.length(); 
			System.out.println(consensus_core+"\t"+String.valueOf(bestscore));
			bestscore=0;
			double[] temp_prior=new double[num_priorbin];
//			double[] temp_dnase=new double[2*DnaseWindow];
//			Double[] max_temp_dnase2=null;
//			double MultiNomConfidence=0;
//			ArrayList<Integer> Rl_stat=new ArrayList<Integer>(SearchEngine2.getSeqNum());
//			ArrayList<Double> Ez_stat=new ArrayList<Double>(SearchEngine2.getSeqNum());
//			ArrayList<Double> PWMscore=new ArrayList<Double>(SearchEngine2.getSeqNum());
//			ArrayList<Double> NBscore=new ArrayList<Double>(SearchEngine2.getSeqNum());
//			ArrayList<Double> MNscore=new ArrayList<Double>(SearchEngine2.getSeqNum());
//			double sumpl=0;
//			double sumRlpl=0;
			
			double lognullprior=Math.log(1.0/num_priorbin);
	
			
			//double log_thresh=motif.getThresh(sampling_ratio, 2*FDR, background)- motiflen*log025;
			double log_thresh=Math.log(1-Prior_EZ)-Math.log(Prior_EZ);
			
			double [][]m_matrix=new double [motiflen+flankingLen*2][4];
			
			//update the loglik matrix
			LinkedList<FastaLocation> Falocs=SearchEngine2.searchPattern(motif, log_thresh);
			System.out.println("number of occurrences: "+String.valueOf(Falocs.size()));
			Prior_EZ=Math.min(motif.core_motiflen*(double)SearchEngine2.getSeqNum()/SearchEngine2.getTotalLength(),Prior_EZ);
			
//			if(Prior_EZ<0.5)
//				Prior_EZ=0.5;
					Iterator<FastaLocation> iter2=Falocs.iterator();
					int count=0;
					int match_seqCount=0;
					int lastseq=-1;
					double max_seqloglik=0;
					double max_logDnaseMNProb=0;
					double max_seqprob_theta=0;
					double max_seqprob_theta_only=0;
					String max_seqsite="";
					while(iter2.hasNext())
					{
						FastaLocation currloc=iter2.next();
						String site="";
						//forward site
						site=SearchEngine2.getSite(currloc.getSeqId(), currloc.getSeqPos()-flankingLen,motiflen+flankingLen*2);
						//check masking
						if(site.indexOf('X')>-1)
							continue;
						
//						if(site.equalsIgnoreCase(""))
//							continue;
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

						if(site.length()!=motiflen+flankingLen*2)
							continue;
						
						double logprob_theta=currloc.Score;//include the bg log_prob in the score
						
						
						
						//double logprob_BG=background.Get_LOGPROB(site.substring((site.length()-motiflen)/2, motiflen));
						double logprior=0;
						int prior_bin=(int)(num_priorbin*((currloc.getSeqPos()+motiflen/2)%currloc.getSeqLen()/(double)currloc.getSeqLen()));
						if(motif.pos_prior.size()!=0)
							logprior=motif.pos_prior.get(prior_bin)-lognullprior;
						double logDnaseProb=0;
//					double	logDnaseNBProb=0,logDnaseMNProb=0;
//						Double[] temp_dnase2=null;
//						if(motif.Dnase_prob!=null)
//						{
//							temp_dnase2=new Double[2*DnaseWindow];
//							Arrays.fill(temp_dnase2, new Double(0));
//							Double[] dnaseseq=DnaseLib.get(currloc.getSeqId());
//							for (int i = 0; i <2*DnaseWindow ; i++) {
//								if(currloc.ReverseStrand)
//								{
//									temp_dnase2[i]=dnaseseq[prior_bin+motif.Dnase_prob.size()-i-1];
//								}
//								else
//								{
//									temp_dnase2[i]=dnaseseq[prior_bin+i];
//								}
//							}
//							logDnaseNBProb=motif.calcLogDnaseNegBinProb(temp_dnase2, 0);
//							logDnaseMNProb=motif.calcLogDnaseMultiNomProb(temp_dnase2, 0);
//							logDnaseProb=logDnaseNBProb+logDnaseMNProb;
//						}
						//double loglik=logprob_theta+logprior-logprob_BG;
						double loglik=logprob_theta+logprior+logDnaseProb+Math.log(Prior_EZ/(1-Prior_EZ));
						if(loglik>200)
							loglik=200;
						double prob_theta=loglik;//Math.exp(loglik)/(Math.exp(loglik)+1);//Math.exp(currloc.Score);
						//double prob_theta_only=Math.exp(logprob_theta+logprior)/(Math.exp(logprob_theta+logprior)+1);
						if(Double.isNaN(prob_theta))
							prob_theta=1;//upper flow
//						if(motif.Dnase_prob!=null)
//						{
//							sumpl+=prob_theta;
//							Double min = (Double) Collections.min(Arrays.asList(temp_dnase2));
//							System.err.println( common.Array2String(temp_dnase2,'\t'));
//							double Rl=0;
//							for (int i = 0; i <2*DnaseWindow ; i++) {
//								double temp=temp_dnase2[i]-min;
//								Rl+=temp;
//								temp=prob_theta*temp;
//								
//								sumRlpl+=temp;
//							}
//							for (int i = 0; i <2*DnaseWindow ; i++) {
//								temp_dnase[i]+=prob_theta*(temp_dnase2[i]-min);
//							}
//							Rl_stat.add((int)Rl);
//							Ez_stat.add(prob_theta);
//							PWMscore.add(logprob_theta+logprior);
//							NBscore.add(logDnaseNBProb/motif.NegBinConfidence);
//							MNscore.add(logDnaseMNProb/motif.MultiNomConfidence);
//							if(logDnaseMNProb>0)
//								MultiNomConfidence+=prob_theta_only;
//							else
//								MultiNomConfidence+=1-prob_theta_only;
//							
//						}
						
						temp_prior[prior_bin]+=prob_theta;//make smaller
						if(OOPS)
							loglik-=common.DoubleMinNormal*Math.abs(currloc.getSeqLen()/2-currloc.getSeqPos()-motiflen/2); //add small bias to center
						if(!OOPS)
							bestscore+=loglik;
						if(OOPS&&currloc.getSeqId()!=lastseq)
						{
							match_seqCount++;
							bestscore+= max_seqloglik;//NegBinFunction.plogis(max_seqloglik);
							lastseq=currloc.getSeqId();
							max_seqloglik=0;
							//System.out.println(max_seqsite.toUpperCase());
							if(max_seqsite.length()==motiflen+flankingLen*2)
							{
								for (int i = 0; i < motiflen+flankingLen*2; i++) {
									int symid=common.acgt(max_seqsite.charAt(i));
									if(symid>3)
										{
											for (int j = 0; j < 4; j++) {
												m_matrix[i][j]+=0.25;
											}
										
										}
									else
									m_matrix[i][symid]+=max_seqprob_theta;
								}
								
								
//								if(motif.Dnase_prob!=null)
//								{
//									sumpl+=max_seqprob_theta;
//									
//									double Rl=0;
//									for (int i = 0; i <2*DnaseWindow ; i++) {
//										double temp=max_temp_dnase2[i];
//										Rl+=temp;
//										temp=prob_theta*temp;
//										temp_dnase[i]+=temp;
//										sumRlpl+=temp;
//									}
//									Rl_stat.add((int)Rl);
//									Ez_stat.add(max_seqprob_theta_only);
//									if(max_logDnaseMNProb>0)
//										MultiNomConfidence+=max_seqprob_theta_only;
//									else
//										MultiNomConfidence+=1-max_seqprob_theta_only;
//									
//								}
							}
//							max_temp_dnase2=temp_dnase2;
						}
						
						if(max_seqloglik<loglik)
						{
							max_seqloglik=loglik;
//							max_logDnaseMNProb=logDnaseMNProb;
							max_seqprob_theta=prob_theta;
							//max_seqprob_theta_only=prob_theta_only;
							max_seqsite=site;
//							max_temp_dnase2=temp_dnase2;
//							if(prob_theta>1000000)
//								continue;
						}
						
						if(!OOPS)
						{
						for (int i = 0; i < motiflen+flankingLen*2; i++) {
							int symid=common.acgt(site.charAt(i));
							if(symid>3)
								continue;
							m_matrix[i][symid]+=prob_theta;//prob_theta;
						}
						

						}
					}
					
					
					//last instance for OOPS
					if(OOPS)
					{
						match_seqCount++;
						bestscore+=max_seqloglik; //NegBinFunction.plogis(max_seqloglik);
						
						max_seqloglik=0;
						//System.out.println(max_seqsite.toUpperCase());
						if(max_seqsite.length()==motiflen+flankingLen*2)
						for (int i = 0; i < motiflen+flankingLen*2; i++) {
							int symid=common.acgt(max_seqsite.charAt(i));
							if(symid>3)
								{
									for (int j = 0; j < 4; j++) {
										m_matrix[i][j]+=max_seqprob_theta;
									}
								
								}
							else
							m_matrix[i][symid]+=max_seqprob_theta;
						}
					}
					
					if(lastscore>=bestscore||match_seqCount<min_support_ratio*SearchEngine2.getSeqNum())
						break;
					else
					{
						if(flankingLen==0)
						{
//							if(lastscore!=1)
//								motif.Score*=bestscore/lastscore;
//							else
						motif.Score=bestscore;//-Math.log(SearchEngine2.getSeqNum())*2*motif.core_motiflen;
						lastscore=bestscore;
						}
						
							
							
						
							
						
						if(Falocs.size()==0|| stateCodes.contains(Falocs.hashCode()))
							return bestPWM;
						else
							stateCodes.add(Falocs.hashCode());
						
						//motif=new PWM((String[])(MatchSite.toArray(new  String[1])));
						if(bestPWM.Score<motif.Score)
						bestPWM=motif.Clone();
						for (int i = 0; i < m_matrix.length; i++) {
							motif.setWeights(i+motif.head-flankingLen,common.Normalize(m_matrix[i]));
						}
					
					
						motif.pos_prior.clear();
						temp_prior=common.Normalize(temp_prior);
						for (int i = 0; i < temp_prior.length; i++) {
							motif.pos_prior.add(temp_prior[i]);
						}
						
//						if(DnaseLib!=null)
//						{
//							motif.Dnase_prob.clear();
//							double sumRl=0;
//							for (int i = 0; i < temp_dnase.length; i++) {
//								sumRl+=temp_dnase[i];
//							}
//							temp_dnase=common.Normalize(temp_dnase);
//							for (int i = 0; i < temp_dnase.length; i++) {
//								motif.Dnase_prob.add(temp_dnase[i]);
//							}
//							//double newP=motif.DnaseFG.getGamma()*sumpl/(motif.DnaseFG.getGamma()*sumpl+sumRlpl);
//							Prior_EZ=sumpl/Ez_stat.size();
//							if(Prior_EZ>0.9999)
//								Prior_EZ=0.9999;
//							
//							 NegBinFunction solver=new NegBinFunction(Rl_stat, Ez_stat, 0);
//							 double[] paras=null;
//							 if(motif.DnaseFG.getGamma()==motif.DnaseBG.getGamma())
//								 paras=solver.run(null);
//							 else
//							 {
//								 double[] preParas=new double[]{Math.log(motif.DnaseBG.getGamma()),NegBinFunction.qlogis(motif.DnaseBG.getP()),Math.log(motif.DnaseFG.getGamma()),NegBinFunction.qlogis(motif.DnaseFG.getP())};
//								 System.out.println(Arrays.toString(preParas));
//								 paras=solver.run(preParas);								 
//							 }
//							 motif.DnaseBG=new NegativeBinomialDist(Math.exp( paras[0]),NegBinFunction.plogis( paras[1]));
//							 motif.DnaseFG=new NegativeBinomialDist(Math.exp( paras[2]),NegBinFunction.plogis( paras[3]));
//							// LogitFunction logitSolver=new LogitFunction(Ez_stat, PWMscore, NBscore, MNscore);
//							// double[] beta=logitSolver.run(new double[]{motif.NegBinConfidence,motif.MultiNomConfidence});
//							 
//							 motif.NegBinConfidence=0;//beta[0];//solver.FeatureConfidence;
//							 motif.MultiNomConfidence=1;//beta[1];//2*MultiNomConfidence/Ez_stat.size()-1;
//							 System.out.println(motif.DnaseFG+"\t"+motif.DnaseFG.getMean());
//								
//							 System.out.println(motif.DnaseBG+"\t"+motif.DnaseBG.getMean());
//								
//							 System.out.println(Prior_EZ);
//							
//							 System.out.println( motif.NegBinConfidence);
//							 System.out.println( motif.MultiNomConfidence);
//							 DrawDistribution(motif.Dnase_prob,"Dnase_plot.png");
//							 
//						}
						
					}
				
					//not allow to grow in the iterations
					flankingLen=0;
					
					System.out.println("number of occurred sequences: "+String.valueOf(match_seqCount));
			}while(motif.core_motiflen<max_motiflen&&iter_count<=max_iterNum);
			
		if(DnaseLib!=null)
			DrawDistribution(motif.Dnase_prob,"Dnase_plot.png");
		return bestPWM;
	}
	
	
	public PWM Column_Replacement_2(PWM motif)
	{
		double log025=Math.log(0.25);
		if(this.pos_prior.size()>0)
			motif.pos_prior=(ArrayList<Double>) this.pos_prior.clone();
		double bestscore=motif.Score;
		HashSet<Integer> extendedCols=new HashSet<Integer>();
		HashSet<Integer> stateCodes=new HashSet<Integer>();
		int num_priorbin=SearchEngine.getTotalLength()/SearchEngine.getSeqNum()/this.resolution;
		double Prior_EZ=0.5;
		int inst_hash=-1;
		double thresh=motif.getThresh(this.sampling_ratio, this.FDR, this.background);
		LinkedList<FastaLocation> origFalocs=SearchEngine.searchPattern(motif, thresh);
		Prior_EZ=1- Math.exp(background.Get_LOGPROB(motif.Consensus(true)))*SearchEngine.TotalLen/origFalocs.size();
		
		do
		{

			String consensus_core=motif.Consensus(true);
		
		System.out.println(consensus_core+"\t"+String.valueOf(bestscore));
		bestscore=0;
		double[] temp_prior=new double[num_priorbin];		
		double lognullprior=Math.log(1.0/num_priorbin);
		
		LinkedList<FastaLocation> Falocs=origFalocs;
		System.out.println(Falocs.size());
		if(Falocs.size()<min_support_ratio*SearchEngine.getSeqNum() ||stateCodes.contains(consensus_core.hashCode()))
			return motif;
		else
			stateCodes.add(consensus_core.hashCode());
			
	
		double [][] loglik_matrix=new double[motif.columns()][4];

		double[][] count_matrix=new double[motif.columns()][4];
		common.fill2DArray(count_matrix, common.DoubleMinNormal);
		String consensus=motif.Consensus(false);
		//double logN025=log025*(motif.head+motif.tail);
		int motiflen=consensus_core.length();
		int avergeSeqlen=(SearchEngine.TotalLen/SearchEngine.SeqNum);
		double[] single_logprob_bg=new double [4];
		double overlap_loglik=0;
		double overlap_prob=0;
		
		int overlap_pos=-motiflen;
		String ACGT="ACGT";
		//get single sym bgprob
		for (int i = 0; i < single_logprob_bg.length; i++) {
			single_logprob_bg[i]=background.Get_LOGPROB(ACGT.substring(i, i+1));
		}
		//update the loglik matrix
		
			
			
				Iterator<FastaLocation> iter2=Falocs.iterator();
				int count=0;
				int lastseq=-1;
				double [][] max_loglik_matrix=new double[motif.columns()][4];
				double [][] max_count_matrix=new double[motif.columns()][4];
				if(OOPS)
				{
					common.fill2DArray(max_loglik_matrix,Double.MIN_VALUE);
					common.fill2DArray(max_count_matrix,0);
				}
				double sitesperSeq=0;
				while(iter2.hasNext())
				{
					FastaLocation currloc=iter2.next();
				
					String site="";
					//forward site
					if(currloc.ReverseStrand)
						site=SearchEngine.getSite(currloc.getMin()-(motif.columns()-seedlen)/2, motif.columns());
					else
						site=SearchEngine.getSite(currloc.getMin()-(motif.columns()-seedlen)/2, motif.columns());
					if(site.indexOf('X')>-1)
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
					double logprob_theta=motif.scoreWeightMatrix(site.substring(motif.head,motif.head+motiflen));
					double logprob_BG=background.Get_LOGPROB(site.substring(motif.head,motif.head+motiflen));

					if(currloc.ReverseStrand)
					{
						//reverse site
						site=common.getReverseCompletementString(site);			
					}	
					int posbin=(int)(num_priorbin*((currloc.getSeqPos()+motiflen/2)%currloc.getSeqLen()/(double)currloc.getSeqLen()));
					double logprior=0;
					if(motif.pos_prior.size()!=0)
						logprior=Math.log(motif.pos_prior.get(posbin))-lognullprior;
					double logDnaseProb=0;

					double loglik=logprob_theta+logprior-logprob_BG+logDnaseProb+Math.log(Prior_EZ/(1-Prior_EZ));
					if(loglik>200)
						loglik=200;
					double prob_theta=Math.exp(loglik)/(Math.exp(loglik)+1);//Math.exp(currloc.Score);				
					if(Double.isNaN(prob_theta))
						prob_theta=1;
					
					temp_prior[posbin]+=prob_theta;

					if(OOPS)
						loglik-=common.DoubleMinNormal*Math.abs(currloc.getSeqLen()/2-currloc.getSeqPos()-motiflen/2); //add small bias to center
					if(!OOPS)
						bestscore+=loglik;
					if(OOPS&&currloc.getSeqId()!=lastseq)
					{
						if(lastseq!=-1)
							for (int i = 0; i < motif.columns(); i++) {							
	                         for (int symid = 0; symid < 4; symid++) {
	                        	 if(sitesperSeq<1)
	                        		 sitesperSeq=1;
	                        	 if(max_loglik_matrix[i][symid]!=Double.MIN_VALUE)
	                        	 {
									loglik_matrix[i][symid]+=max_loglik_matrix[i][symid]/sitesperSeq;

									count_matrix[i][symid]+=max_count_matrix[i][symid]/sitesperSeq;
	                        	 }
							}

							}
						
						lastseq=currloc.getSeqId();
						sitesperSeq=0;
						overlap_pos=-motiflen;
						common.fill2DArray(max_loglik_matrix,Double.MIN_VALUE);
						common.fill2DArray(max_count_matrix,0);
					}
					if(OOPS)
					{
						if(Math.abs(currloc.getSeqPos()-overlap_pos)<motiflen)
						{
							if(prob_theta>overlap_prob)
							{
								for (int i = 0; i < site.length(); i++) {
									int symid=common.acgt(site.charAt(i));
									if(symid>3)
										continue; //meet new line separator
								max_count_matrix[i][symid]+=prob_theta-overlap_prob;
								if(max_loglik_matrix[i][symid]==Double.MIN_VALUE)
									max_loglik_matrix[i][symid]=prob_theta*loglik-overlap_prob*overlap_loglik;
								else
									max_loglik_matrix[i][symid]+=prob_theta*loglik-overlap_prob*overlap_loglik;
								
								}	
								sitesperSeq+=prob_theta-overlap_prob;
								overlap_prob=prob_theta;
								overlap_loglik=loglik;
							}
							//else just ignore the low score site
						}								
						else
						{
							for (int i = 0; i < site.length(); i++) {
								int symid=common.acgt(site.charAt(i));
								if(symid>3)
								{
									for (int j = 0; j < 4; j++) {
										max_count_matrix[i][j]+=0.25;
									}
									continue;
								}
								
							max_count_matrix[i][symid]+=prob_theta;
							if(max_loglik_matrix[i][symid]==Double.MIN_VALUE)
								max_loglik_matrix[i][symid]=prob_theta*(loglik);
							else
								max_loglik_matrix[i][symid]+=prob_theta*(loglik);
							
							}	
							sitesperSeq+=prob_theta;
							overlap_prob=prob_theta;
							overlap_loglik=loglik;
						}
						overlap_pos=currloc.getSeqPos();

						continue;
					}

					
					for (int i = 0; i < site.length(); i++) {
						int symid=common.acgt(site.charAt(i));
						if(symid>3)
							continue; //meet new line separator
						count_matrix[i][symid]+=prob_theta;
						if(consensus.charAt(i)!='N')
							continue;
						loglik_matrix[i][symid]+=loglik;
					}
			}
	
		
		
		//select the best column replacement
		double maxloglik=0;
		
		int numbestSym=4;
		int bestCol=-1;
		ArrayList<Integer> bestSym=new ArrayList<Integer>(4);
		//use binomial p-value to decide stop extention
		RandomEngine rand=RandomEngine.makeDefault();
         // Binomial binomial=new new Binomial(R)
		
		
		int num_col_cand=0;
		//compute the optimal column assignments
		double [][] optimalcols=new double[count_matrix.length][4];
		for (int i = 0; i < optimalcols.length; i++) {
			double sumCount=0;
			for (int j = 0; j < 4; j++) {
				optimalcols[i][j]=count_matrix[i][j];
				sumCount+=optimalcols[i][j];
			}
			for (int j = 0; j < 4; j++) {
				optimalcols[i][j]/=sumCount;
			}
		}
		
		
			for (int i = 0; i < motif.columns(); i++) {
				if(consensus.charAt(i)!='N')
					continue;
				num_col_cand++;
				TreeMap<Double,Integer> orderSym=new TreeMap<Double,Integer>();
				double temploglik=0;
				double sumtemp=0;
				for (int j = 0; j < 4; j++) {
					double temp=loglik_matrix[i][j];
					orderSym.put(temp-common.DoubleMinNormal*j, j);
					sumtemp+=temp;
				}
				//-single_logprob_bg[symid]
				double boundaryLoss=0.25*Math.min(Math.abs(i-motif.head), Math.abs(motiflen+motif.head-i-1))/(double)avergeSeqlen;
				if((motif.head<i&&i<motiflen+motif.head)||OOPS)
					boundaryLoss=0;
				
				//System.out.println(sumtemp);
				//only consider best two sym
				int c=0;
				double sumcount=0;
				for(Double key: orderSym.descendingKeySet()) {
				
					if(key<0||c>=numbestSym)
						break;
					int symid=orderSym.get(key);
					temploglik+=count_matrix[i][symid]*(Math.log(optimalcols[i][symid])-single_logprob_bg[symid]); //key*(1+boundaryLoss);
					sumcount+=count_matrix[i][symid];
					c++;			
				}
				temploglik=temploglik-common.lnEntropy(optimalcols[i]);
				if(temploglik>maxloglik)
				{
					maxloglik=temploglik;
					bestCol=i;
					bestSym.clear();
					for(Double key: orderSym.descendingKeySet()) {
						int col=orderSym.get(key);
						if(key<0)//||bestSym.size()==numbestSym
							break;
						bestSym.add(col);
						
					}
					
					
				}
				
			}
			
			if(bestCol==-1)
				break;

			double X=count_matrix[bestCol][bestSym.get(0)];
			double total=0;
			for (int j = 0; j < 4; j++) {
				total+=count_matrix[bestCol][j];
			}
			Binomial binomial=new Binomial((int)total,Math.exp(single_logprob_bg[bestSym.get(0)]),rand);
			double pvalue=1.0-binomial.cdf((int)X);
			if(pvalue>this.FDR/(num_col_cand*4))
				break;
			double [] repColumnValue=new double[4];
			Arrays.fill(repColumnValue, common.DoubleMinNormal);
			
			double max_sumNorm=0;
			double max_sumCount=0;
			double max_symCount=0;

			for (int i = 1; i <= bestSym.size(); i++) {
				double sumNorm=0;
				
				double sumcount=0;
				for (int j = 0; j < i; j++)
					sumcount+=count_matrix[bestCol][bestSym.get(j)];
				for (int j = 0; j < i; j++) {
					double temp=loglik_matrix[bestCol][bestSym.get(j)];
						sumNorm+=(temp)+count_matrix[bestCol][bestSym.get(j)]*Math.log(count_matrix[bestCol][bestSym.get(j)]/sumcount);
					
				}

				if(sumNorm>max_sumNorm)
				{
					for (int j = 0; j < i; j++) 
					{
						repColumnValue[bestSym.get(j)]=count_matrix[bestCol][bestSym.get(j)]/sumcount;
					}
					repColumnValue=common.Normalize(repColumnValue);
					max_sumNorm=sumNorm;
					max_symCount=i;
				}
				
				
			}
			
			for (int j = 0; j < 4; j++) {
				max_sumCount+=count_matrix[bestCol][j];
			}
			Prior_EZ=max_sumCount/Falocs.size();

			if(debug)
			for (int j = 0; j < count_matrix.length; j++) {
				System.out.println(Arrays.toString(loglik_matrix[j]));
			}
			
			System.out.println(max_sumNorm);
			double boundaryLoss=0.25*Math.min(Math.abs(bestCol-motif.head), Math.abs(motiflen+motif.head-bestCol-1))/(double)avergeSeqlen;
			if((motif.head<bestCol&&bestCol<motiflen+motif.head)||OOPS)
				boundaryLoss=0;
			max_sumNorm*=(1+boundaryLoss);
			
			if(max_sumNorm<=bestscore||max_sumCount<SearchEngine.getSeqNum()*min_support_ratio)
				break;
			else
			{
				motif.Score=max_sumNorm;
				bestscore=max_sumNorm;
			}


			//update motif column value
			motif.setWeights(bestCol, repColumnValue);
			double sumPrior=0;
			for (int i = 0; i < temp_prior.length; i++) {
				sumPrior+=temp_prior[i];
			}
			//also maximize loglik for update other column
			Iterator<Integer> iter1=extendedCols.iterator();
			while(iter1.hasNext())
			{
				int ecol=iter1.next();
				motif.setWeights(ecol, optimalcols[ecol]);
				
			}
			
			motif.inst_FDR=1-Prior_EZ;
			extendedCols.add(bestCol);
			motif.pos_prior.clear();
			for (int i = 0; i < temp_prior.length; i++) {
				temp_prior[i]/=sumPrior;
				motif.pos_prior.add(temp_prior[i]);
			}

			if(debug)
				motif.print();
			
			
		}while(true);
	
		return motif;
	}
	
	
	public PWM Column_Replacement_(PWM motif)
	{
		double log025=Math.log(0.25);
		if(this.pos_prior.size()>0)
			motif.pos_prior=(ArrayList<Double>) this.pos_prior.clone();
		double bestscore=motif.Score;
		HashSet<Integer> extendedCols=new HashSet<Integer>();
		HashSet<Integer> stateCodes=new HashSet<Integer>();
		int num_priorbin=SearchEngine.getTotalLength()/SearchEngine.getSeqNum()/this.resolution;
		double Prior_EZ=0.5;
		int inst_hash=-1;
//		if(motif.Dnase_prob!=null)
//			motif.Log_Fab_ratio=Num.lnGamma(motif.Dnase_prob.size());
		do
		{

			String consensus_core=motif.Consensus(true);
		
		System.out.println(consensus_core+"\t"+String.valueOf(bestscore));
		bestscore=0;
		double[] temp_prior=new double[num_priorbin];
//		double MultiNomConfidence=0;
//		double[] temp_dnase=new double[2*DnaseWindow];
//		ArrayList<Integer> Rl_stat=new ArrayList<Integer>(SearchEngine2.getSeqNum());
//		ArrayList<Double> Ez_stat=new ArrayList<Double>(SearchEngine2.getSeqNum());
//		double sumpl=0;
//		double sumRlpl=0;
		
		double lognullprior=Math.log(1.0/num_priorbin);
		double thresh=motif.getThresh(this.sampling_ratio, this.FDR, this.background);
		
		LinkedList<FastaLocation> Falocs=SearchEngine.searchPattern(motif, thresh);
		System.out.println(Falocs.size());
		if(Falocs.size()<min_support_ratio*SearchEngine.getSeqNum() ||stateCodes.contains(Falocs.hashCode()))
			return motif;
		else
			stateCodes.add(Falocs.hashCode());
			
	
		double [][] loglik_matrix=new double[motif.columns()][4];
		
		double[][] count_matrix=new double[motif.columns()][4];
		common.fill2DArray(count_matrix, common.DoubleMinNormal);
		String consensus=motif.Consensus(false);
		//double logN025=log025*(motif.head+motif.tail);
		int motiflen=consensus_core.length();
		int avergeSeqlen=(SearchEngine.TotalLen/SearchEngine.SeqNum);
		double[] single_logprob_bg=new double [4];
		String ACGT="ACGT";
		//get single sym bgprob
		for (int i = 0; i < single_logprob_bg.length; i++) {
			single_logprob_bg[i]=background.Get_LOGPROB(ACGT.substring(i, i+1));
		}
		//update the loglik matrix
		
			
			
				Iterator<FastaLocation> iter2=Falocs.iterator();
				int count=0;
				int lastseq=-1;
				double [][] max_loglik_matrix=new double[motif.columns()][4];
				if(OOPS)
					common.fill2DArray(max_loglik_matrix,Double.MIN_VALUE);
				while(iter2.hasNext())
				{
					FastaLocation currloc=iter2.next();
					double logprob_theta=currloc.Score;
					String site="";
					//forward site
					if(currloc.ReverseStrand)
						site=SearchEngine.getSite(currloc.getMin()-motif.tail, motif.columns());
					else
						site=SearchEngine.getSite(currloc.getMin()-motif.head, motif.columns());
					if(site.indexOf('X')>-1)
						continue;
					
					double logprob_BG=background.Get_LOGPROB(site.substring(motif.head,motif.head+motiflen));
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

					if(currloc.ReverseStrand)
					{
						//reverse site
						site=common.getReverseCompletementString(site);
						
					}	
					int posbin=(int)(num_priorbin*((currloc.getSeqPos()+motiflen/2)%currloc.getSeqLen()/(double)currloc.getSeqLen()));
					double logprior=0;
					if(motif.pos_prior.size()!=0)
						logprior=motif.pos_prior.get(posbin)-lognullprior;
					double logDnaseProb=0;
//					double logDnaseNBProb=0;
//					double logDnaseMNProb=0;
//					Double[] temp_dnase2=null;
//					if(motif.Dnase_prob!=null)
//					{
//						temp_dnase2=new Double[2*DnaseWindow];
//						Arrays.fill(temp_dnase2, new Double(0));
//						Double[] dnaseseq=DnaseLib.get(currloc.getSeqId());
//						for (int i = 0; i <2*DnaseWindow ; i++) {
//							if(currloc.ReverseStrand)
//							{
//								temp_dnase2[i]=dnaseseq[posbin+motif.Dnase_prob.size()-i-1];
//							}
//							else
//							{
//								temp_dnase2[i]=dnaseseq[posbin+i];
//							}
//						}
//						logDnaseNBProb=motif.calcLogDnaseNegBinProb(temp_dnase2, 0);
//						logDnaseMNProb=motif.calcLogDnaseMultiNomProb(temp_dnase2, 0);
//						logDnaseProb=logDnaseNBProb+logDnaseMNProb;
//					}
					double loglik=logprob_theta+logprior-logprob_BG+logDnaseProb+Math.log(Prior_EZ/(1-Prior_EZ));
					if(loglik>200)
						loglik=200;
					double prob_theta=Math.exp(loglik)/(Math.exp(loglik)+1);//Math.exp(currloc.Score);
					//double prob_theta_only=Math.exp(logprob_theta+logprior)/(Math.exp(logprob_theta+logprior)+1);
//					if(logDnaseMNProb>0)
//						MultiNomConfidence+=prob_theta_only;
//					else
//						MultiNomConfidence+=1-prob_theta_only;
					
					if(Double.isNaN(prob_theta))
						prob_theta=1;
					
					temp_prior[posbin]+=prob_theta;
//					if(motif.Dnase_prob!=null)
//					{
//						sumpl+=prob_theta;
//						Double min = (Double) Collections.min(Arrays.asList(temp_dnase2));
//						double Rl=0;
//					
//						for (int i = 0; i <2*DnaseWindow ; i++) {
//							double temp=temp_dnase2[i]-min;
//
//							Rl+=temp;
//							temp=prob_theta*temp;
//							temp_dnase[i]+=loglik*(temp_dnase2[i]-min);
//							sumRlpl+=temp;
//						}
//						Rl_stat.add((int)Rl);
//						Ez_stat.add(prob_theta);
//						
//					}
					if(OOPS)
						loglik-=common.DoubleMinNormal*Math.abs(currloc.getSeqLen()/2-currloc.getSeqPos()-motiflen/2); //add small bias to center
					if(!OOPS)
						bestscore+=loglik;
					if(OOPS&&currloc.getSeqId()!=lastseq)
					{
						if(lastseq!=-1)
							for (int i = 0; i < site.length(); i++) {							
	                         for (int symid = 0; symid < 4; symid++) {
	                        	 if(max_loglik_matrix[i][symid]!=Double.MIN_VALUE)
	                        	 {
									loglik_matrix[i][symid]+=max_loglik_matrix[i][symid]-single_logprob_bg[symid];
									count_matrix[i][symid]+=1;
	                        	 }
							}

							}
						bestscore+=loglik;
						lastseq=currloc.getSeqId();
						common.fill2DArray(max_loglik_matrix,Double.MIN_VALUE);
					}
					if(OOPS)
					{
						for (int i = 0; i < site.length(); i++) {
//							if(consensus.charAt(i)!='N')
//								continue;
							int symid=common.acgt(site.charAt(i));
							if(symid>3)
								continue; //meet new line separator
							if(loglik>max_loglik_matrix[i][symid])
								max_loglik_matrix[i][symid]=loglik;
							
						}
					
						continue;
					}
//					if(loglik<0)
//						System.out.println(site);
					
					for (int i = 0; i < site.length(); i++) {
						int symid=common.acgt(site.charAt(i));
						if(symid>3)
							continue; //meet new line separator
						count_matrix[i][symid]+=1;
						if(consensus.charAt(i)!='N')
							continue;
						loglik_matrix[i][symid]+=loglik-single_logprob_bg[symid];
					}
			}
	
		
		
		//select the best column replacement
		double maxloglik=Double.MIN_VALUE;
		
		int numbestSym=1;
		int bestCol=-1;
		ArrayList<Integer> bestSym=new ArrayList<Integer>(4);
		//use binomial p-value to decide stop extention
		RandomEngine rand=RandomEngine.makeDefault();
         // Binomial binomial=new new Binomial(R)
		
		
		int num_col_cand=0;
			for (int i = 0; i < motif.columns(); i++) {
				if(consensus.charAt(i)!='N')
					continue;
				num_col_cand++;
				TreeMap<Double,Integer> orderSym=new TreeMap<Double,Integer>();
				double temploglik=0;
				double sumtemp=0;
				for (int j = 0; j < 4; j++) {
					double temp=loglik_matrix[i][j];
					orderSym.put(temp-common.DoubleMinNormal*j, j);
					sumtemp+=temp;
				}
				
				double boundaryLoss=0.25*Math.min(Math.abs(i-motif.head), Math.abs(motiflen+motif.head-i-1))/(double)avergeSeqlen;
				if((motif.head<i&&i<motiflen+motif.head)||OOPS)
					boundaryLoss=0;
				
				//System.out.println(sumtemp);
				//only consider best two sym
				int c=0;
				for(Double key: orderSym.descendingKeySet()) {
				
					if(key<0||c>=numbestSym)
						break;
					temploglik+=key*(1+boundaryLoss);
					c++;
					
				}
				
				if(temploglik>maxloglik)
				{
					maxloglik=temploglik;
					bestCol=i;
					bestSym.clear();
					for(Double key: orderSym.descendingKeySet()) {
						int col=orderSym.get(key);
						if(key<0)//||bestSym.size()==numbestSym
							break;
						bestSym.add(col);
						
					}
					
					
				}
				
			}
			
			if(bestCol==-1)
				break;
			double X=count_matrix[bestCol][bestSym.get(0)];
			Binomial binomial=new Binomial(Falocs.size(),Math.exp(single_logprob_bg[bestSym.get(0)]),rand);
			double pvalue=1.0-binomial.cdf((int)X);
			if(pvalue>this.FDR/(num_col_cand*4))
				break;
			
			double [] repColumnValue=new double[4];
			Arrays.fill(repColumnValue, common.DoubleMinNormal);
			
			double max_sumNorm=0;
			double max_sumCount=0;
			double max_symCount=0;

			for (int i = 1; i <= bestSym.size(); i++) {
				double sumNorm=0;
				double sumcount=0;
				for (int j = 0; j < i; j++)
					sumcount+=count_matrix[bestCol][bestSym.get(j)];
				for (int j = 0; j < i; j++) {
					double temp=loglik_matrix[bestCol][bestSym.get(j)];
						sumNorm+=(temp)+count_matrix[bestCol][bestSym.get(j)]*Math.log(count_matrix[bestCol][bestSym.get(j)]/sumcount);
					
				}

				if(sumNorm>max_sumNorm)
				{
					for (int j = 0; j < i; j++) 
					{
						repColumnValue[bestSym.get(j)]=count_matrix[bestCol][bestSym.get(j)]/sumcount;
					}
					max_sumNorm=sumNorm;
					max_sumCount=sumcount;
					max_symCount=i;
				}
				
				
			}

			if(debug)
			for (int j = 0; j < count_matrix.length; j++) {
				System.out.println(Arrays.toString(loglik_matrix[j]));
			}
			
			System.out.println(max_sumNorm);
			double boundaryLoss=0.25*Math.min(Math.abs(bestCol-motif.head), Math.abs(motiflen+motif.head-bestCol-1))/(double)avergeSeqlen;
			if((motif.head<bestCol&&bestCol<motiflen+motif.head)||OOPS)
				boundaryLoss=0;
			max_sumNorm*=(1+boundaryLoss);
			
			if(max_sumNorm<=bestscore||max_sumCount<SearchEngine.getSeqNum()*min_support_ratio)
				break;
			else
			{
				motif.Score*=max_sumNorm/bestscore;
				bestscore=max_sumNorm;
			}


			//update motif column value
			motif.setWeights(bestCol, repColumnValue);
			double sumPrior=0;
			for (int i = 0; i < temp_prior.length; i++) {
				sumPrior+=temp_prior[i];
			}
			//also maximize loglik for update other column
			Iterator<Integer> iter1=extendedCols.iterator();
			while(iter1.hasNext())
			{
				int ecol=iter1.next();
				motif.setWeights(ecol, common.Normalize(count_matrix[ecol]));
				
			}
			
			
			extendedCols.add(bestCol);
			motif.pos_prior.clear();
			for (int i = 0; i < temp_prior.length; i++) {
				temp_prior[i]/=sumPrior;
				motif.pos_prior.add(temp_prior[i]);
			}
//			if(DnaseLib!=null)
//			{
//				motif.Dnase_prob.clear();
//				temp_dnase=common.Normalize(temp_dnase);
//				for (int i = 0; i < temp_dnase.length; i++) {
//					motif.Dnase_prob.add(temp_dnase[i]);
//				}
//				Prior_EZ=sumpl/Ez_stat.size();
//				if(Prior_EZ>0.9999)
//					Prior_EZ=0.9999;
//				
//				 NegBinFunction solver=new NegBinFunction(Rl_stat, Ez_stat, 0);
//				 double[] paras=null;
//				 if(motif.DnaseFG.getGamma()==motif.DnaseBG.getGamma())
//					 paras=solver.run(null);
//				 else
//				 {
//					 double[] preParas=new double[]{Math.log(motif.DnaseBG.getGamma()),NegBinFunction.qlogis(motif.DnaseBG.getP()),Math.log(motif.DnaseFG.getGamma()),NegBinFunction.qlogis(motif.DnaseFG.getP())};
//					 System.out.println(Arrays.toString(preParas));
//					 paras=solver.run(preParas);
//				 }
//				 motif.DnaseBG=new NegativeBinomialDist(Math.exp( paras[0]),NegBinFunction.plogis( paras[1]));
//				 motif.DnaseFG=new NegativeBinomialDist(Math.exp( paras[2]),NegBinFunction.plogis( paras[3]));
//				 motif.NegBinConfidence=0;//solver.FeatureConfidence;
//				 motif.MultiNomConfidence=1;//2*MultiNomConfidence/Ez_stat.size()-1;
////				 motif.samplesize=Ez_stat.size();
////			
////				 double[] ones=new double[motif.Dnase_prob.size()];
////				 Arrays.fill(ones, 1.0);
////				 double sum=0;
////				 motif.Log_Fab_ratio=0;
////					for (int i = 0; i <motif.Dnase_prob.size(); i++) {
////						double alpha=motif.Dnase_prob.get(i) *motif.samplesize;
////						 motif.Log_Fab_ratio-=Num.lnGamma(alpha);
////						sum+=alpha;
////					}
////					 motif.Log_Fab_ratio+=Num.lnGamma(sum);
//					
//				DrawDistribution(motif.Dnase_prob,"Dnase_plot.png");
//				System.out.println(motif.DnaseFG+"\t"+motif.DnaseFG.getMean());
//				System.out.println(motif.DnaseBG+"\t"+motif.DnaseBG.getMean());
//				System.out.println(Prior_EZ);
//				System.out.println(solver.FeatureConfidence);
//				System.out.println( motif.MultiNomConfidence);
//			}
			if(debug)
				motif.print();
			
			
		}while(true);
	
		return motif;
	}
	
	
	
	
	//for fix PWM score and BG score, but there is position prior
	double sumLLR(LinkedList<FastaLocation> LocList,double pwmloglik,double bgloglik)
	{
		double score=0;
		
       Iterator<FastaLocation> iter=LocList.iterator();
		int num_priorbin=SearchEngine.getTotalLength()/SearchEngine.getSeqNum()/this.resolution;
		double lognullprior=Math.log(1.0/num_priorbin);
		double max_seqloglik=0;
		int lastseq=-1;
		int lastpos=-seedlen;
		while(iter.hasNext())
		{
			FastaLocation currloc=iter.next();
			double logprior=0;
			if(pos_prior.size()!=0)
				logprior=pos_prior.get( pos_prior.size()*currloc.getSeqPos()/currloc.getSeqLen())-lognullprior;
			
			double loglik=pwmloglik+logprior-bgloglik;
			if(OOPS)
			loglik-=common.DoubleMinNormal*Math.abs(currloc.getSeqLen()/2-currloc.getSeqPos()-seedlen/2); //add small bias to center
			if(OOPS&&currloc.getSeqId()!=lastseq)
			{
				if(lastseq!=-1)
					score+=max_seqloglik;
				lastseq=currloc.getSeqId();
				lastpos=-seedlen;
				max_seqloglik=0;
			}
			if(OOPS)
			{
				//if(loglik>max_seqloglik)
				if(Math.abs(lastpos-currloc.getSeqPos())>seedlen)
				max_seqloglik+=loglik;
				
				lastpos=currloc.getSeqPos();
				continue;
			}
			
			score+=loglik;
		}
		if(OOPS)
			score+=max_seqloglik;
		

		return score;
	}
	
	//for fix PWM score but there is 'N' inside the pattern
	double sumLLR(LinkedList<FastaLocation> LocList,double pwmloglik,int motiflen)
	{
		double score=0;
		int count=0;
		int num_priorbin=SearchEngine.getTotalLength()/SearchEngine.getSeqNum()/this.resolution;
		double lognullprior=Math.log(1.0/num_priorbin);
		Iterator<FastaLocation> iter2=LocList.iterator();
		while(iter2.hasNext())
		{
			FastaLocation currloc=iter2.next();
			String site="";
			//forward site
			site=SearchEngine.getSite(currloc.getMin(),motiflen);
			if(site.equalsIgnoreCase(""))
				continue;
			if(site.indexOf("N")!=-1)
				continue;
			double logprob_BG=background.Get_LOGPROB(site);
			if(count>=SearchEngine.forwardCount)
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
		options.addOption("minw",true,"minimum size of motif binding region (default 200bp)");
		options.addOption("maxw",true, "maximum size of motif binding region (default 600bp)");
		options.addOption("mask",false,"whether marking the top motif location in order to find co-motif");
		options.addOption("oops",false,"whether assuming only one occurrence per sequence");
		options.addOption("FDR",true,"fasle positive rate (default 0.05)");
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
				motifFinder.max_motiflen=Integer.parseInt(cmd.getOptionValue("maxlen"));
			}
			if(cmd.hasOption("minw"))
			{
				motifFinder.starting_windowsize=Integer.parseInt( cmd.getOptionValue("minw"));
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
	
			
		} catch (ParseException e) {
			// TODO Auto-generated catch block
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp( "JPomoda", options );
			return;
		}
		
		//initialize Pomoda 
		motifFinder.initialize();

		
		//get seed motifs
		//ArrayList<PWM>  seedPWMs=motifFinder.getSeedMotifs();
		ArrayList<PWM>  seedPWMs=new ArrayList<PWM>();
		
		double topseed_Score=0;
		
		File file = new File(motifFinder.outputPrefix+"jpomoda_raw.pwm"); 
		try {
			seedPWMs.add(new PWM(new String[]{"NNNNNNNNNNNTGACCNNNNNNNNNNN"}));
			BufferedWriter writer = new BufferedWriter(new FileWriter(file));
			TreeMap<Double, PWM> sortedPWMs=new TreeMap<Double, PWM>();
		//extend and refine motifs
		for (int i = 0; i < seedPWMs.size(); i++) {
			PWM motif=seedPWMs.get(i);
			motif.Score=1;
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
			System.out.println("Relaxing...");
			seedPWMs.set(i, motifFinder.Relax_Seed_2(seedPWMs.get(i)));

			seedPWMs.get(i).Name="Motif"+String.valueOf(i+1);
			// to make different length comparable ,need to consider the instance coverage
//			seedPWMs.get(i).Score=seedPWMs.get(i).Score/seedPWMs.get(i).inst_coverage;//corrected score
			if(motifFinder.maskflag&&seedPWMs.get(i).Score>(topseed_Score))
			{
				//do something to mark the locations in SearchEngine
				System.out.println("Masking...");
				motifFinder.Masking(seedPWMs.get(i));
				topseed_Score=seedPWMs.get(i).Score;
				motifFinder.SearchEngine2.EnableBackground(motifFinder.background);
			}
			
			writer.write(seedPWMs.get(i).toString());
		
		
		}

		writer.close();
		
		file = new File(motifFinder.outputPrefix+"jpomoda_clust.pwm"); 
		writer = new BufferedWriter(new FileWriter(file));
		//clustering motif, re-initialize
		System.out.println("Clustering...");
		motifFinder.SearchEngine2.build_index(motifFinder.inputFasta);
		sortedPWMs.clear();
		PWMcluster clustering=new PWMcluster(motifFinder);
	
			ArrayList<PWM>  clusterPWMs=clustering.Clustering_(seedPWMs,motifFinder.num_motif);
	
			for(PWM pwm:clusterPWMs)
			{
				sortedPWMs.put(pwm.Score, pwm); //desc order
			}
			System.out.println("Clustered Motifs:");
			int c=0;
//			GapImprover gimprover=new GapImprover(motifFinder);
			PWMevaluator evaluator=new PWMevaluator(motifFinder);
			for(Double key:sortedPWMs.descendingKeySet())
			{
				sortedPWMs.get(key).Name="Motif_clust"+String.valueOf(c+1);
				c++;
				System.out.println(sortedPWMs.get(key).Consensus(true)+'\t'+sortedPWMs.get(key).Score);
				writer.write(sortedPWMs.get(key).toString());
				motifFinder.DrawDistribution(sortedPWMs.get(key).pos_prior,motifFinder.outputPrefix+"/"+sortedPWMs.get(key).Name+"_dist.png");
//				GapPWM gpwm=gimprover.fillDependency(sortedPWMs.get(key));
//				
				if(motifFinder.DnaseLib==null)
					System.out.println("PWM AUC:"+ evaluator.calcAUC(sortedPWMs.get(key),null));
				else
					System.out.println("PWM AUC:"+ evaluator.calcAUC(sortedPWMs.get(key),motifFinder.DnaseLib,null));
//				System.out.println("Gap improved PWM AUC:"+gimprover.AUCtest(gpwm));
			}
			
			writer.close();
		
		} 
		catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IllegalAlphabetException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IllegalSymbolException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		}
		

		

	

}
