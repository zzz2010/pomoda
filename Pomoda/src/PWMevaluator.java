import java.awt.Color;
import java.awt.image.BufferedImage;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Random;
import java.util.TreeMap;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.lang.ArrayUtils;
import org.biojava.bio.dist.DistributionTools;
import org.biojava.bio.dist.UniformDistribution;
import org.biojava.bio.seq.DNATools;
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

import auc.AUCCalculator;
import auc.Confusion;


public class PWMevaluator {

	/**
	 * @param args
	 */
	public String outputPrefix="./";
	public String inputFasta;
	public String ctrlFasta="";
	public boolean OOPS=true; //only one dependence per sequence
	public boolean removeBG=true; //false:uniform BG assume
	public String bgmodelFile="";
	LinearEngine SearchEngine;
	public double sampling_ratio=1;
	public double FDR=0.01;
	public double entropyThresh=1;
	
	public int resolution=10;
	
	public BGModel background;
	
	HashMap<String,XYSeries> ROCdata=new HashMap<String, XYSeries>();
	
	public PWMevaluator(Pomoda motiffinder)
	{
		//super("");

		SearchEngine=motiffinder.SearchEngine2;
		sampling_ratio=motiffinder.sampling_ratio;
		FDR=motiffinder.FDR;
		background=motiffinder.background;
		resolution=motiffinder.resolution;
		
	}
	
	public PWMevaluator()
	{//super("");
		
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
		     background.BuildModel(inputFasta, bg_markov_order+1); //1-order bg
		     background.SaveModel(inputFasta+".bgobj");
			}
			else
			{
			     background.BuildModel(ctrlFasta, bg_markov_order+1); //3-order bg
			     background.SaveModel(ctrlFasta+".bgobj");
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
			similarList.add(sortedPWMs.get(key));
			
			if(similarList.size()==K)
				break;
		}
		
		return similarList;
	}
	
	public double calcAUC(PWM motif,LinkedList<String> sequences)
	{
		 LinearEngine SearEngine=null;
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
     						Random r=new Random();
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
		options.addOption("i", true, "input fasta file");
		options.addOption("pwm", true, "input PWM file");
		options.addOption("c", true, "control fasta file");
		options.addOption("convert", false, "convert input PWM file to the transfac format");
		options.addOption("roc", false, "compute AUC and draw ROC curve for the given pwm file");
		options.addOption("match", true, "find similar motifs in known PWM library (path to the library, e.g., jaspar.pwm)");
		options.addOption("bgmodel", true, "background model file");
		options.addOption("prefix", true, "output directory");
		options.addOption("ratio",true, "sampling ratio (default 1)");
		options.addOption("thresh",true, "minimum entropy threshold for considering a position as a gap(default 0.5)");
		options.addOption("FDR",true,"fasle positive rate");
		String inputPWM;
		CommandLineParser parser = new GnuParser();
		PWMevaluator evaluator=new PWMevaluator();
		boolean rocflag=false;
		boolean convertflag=false;
		LinkedList<PWM> PWMLibrary=null;
		
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
			}
			if(cmd.hasOption("bgmodel"))
			{
				evaluator.bgmodelFile=cmd.getOptionValue("bgmodel");
			}
			if(cmd.hasOption("match"))
			{
				PWMLibrary=common.LoadPWMFromFile(cmd.getOptionValue("match"));
			}
			if(cmd.hasOption("roc"))
			{
				rocflag=true;
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
				evaluator.entropyThresh=Double.parseDouble(cmd.getOptionValue("thresh"));
			}
			if(cmd.hasOption("FDR"))
			{
				evaluator.FDR=Double.parseDouble(cmd.getOptionValue("FDR"));
			}

		} catch (ParseException e) {
			// TODO Auto-generated catch block
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp( "PWMevaluator", options );
			return;
		}
		File file = new File(inputPWM+"_eval.txt"); 
		try {
			BufferedWriter writer = new BufferedWriter(new FileWriter(file));
		
		if(rocflag)
		{
		evaluator.initialize();
		
		LinkedList<PWM> pwmlist=common.LoadPWMFromFile(inputPWM);
		Iterator<PWM> iter=pwmlist.iterator();
		writer.write("AUC Result:\n");
		TreeMap<Double,PWM> sortedPWMs=new TreeMap<Double,PWM>();
		while(iter.hasNext())
		{
			PWM p1=iter.next();
			double auc=evaluator.calcAUC(p1, null);
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
				writer2.write(sortedPWMs.get(key).toString());
			}
			writer2.close();
			convertflag=false;
		}
		evaluator.DrawROC(inputPWM+"_roc.png");
		//evaluator.pack();
	     //   RefineryUtilities.centerFrameOnScreen(evaluator);
	      //  evaluator.setVisible(true);
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
				int K=5;
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
