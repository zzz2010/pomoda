import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.util.ArrayList;
import java.util.HashMap;
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
import org.jfree.data.xy.XYSeries;
import org.w3c.dom.stylesheets.LinkStyle;

import auc.AUCCalculator;
import auc.Confusion;

import cern.colt.function.DoubleDoubleFunction;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix1DProcedure;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.jet.math.PlusMult;


public class MultPWMScorer {
	public String outputPrefix="./";
	public String inputFasta;
	public String ctrlFasta="";
	//public boolean removeBG=false; //false:uniform BG assume
	public String bgmodelFile="";
	public int BGFold=1;
	LinearEngine SearchEngine;
	LinearEngine BGSearchEngine;
	public boolean StrandSpec=false;
	public boolean linearOrder=true;
	public static double maxFPdraw=1;
	public int resolution=10;
	
	public BGModel background=null;
	static List<PWM> pwmlist=null;
	static HashMap<String,XYSeries> ROCdata=new HashMap<String, XYSeries>();
	private boolean PBMflag=false;
	public static int max_motif_span=30;
	

	
	public MultPWMScorer()
	{//super("");
		
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
		
		if(StrandSpec)
		{
			SearchEngine.singleStrand=true;
		}
		if(background==null)
		{
			background=new BGModel();
			int bg_markov_order=0;
			if(!ctrlFasta.isEmpty())
			{
				BGSearchEngine=new LinearEngine(4);
				 background.BuildModel(this.ctrlFasta, 4); //4-order bg
				BGSearchEngine.build_index(this.ctrlFasta);
				//removeBG=true;
			}
			
			else if(!bgmodelFile.isEmpty())
			{
				//removeBG=true;
					background.LoadModel(bgmodelFile);
			}
			
			else
			{
				background=BGModel.CreateUniform();
			}

		}
		
		//add uniform bg for both search engine, ensure the score can be positive
		SearchEngine.EnableBackground(BGModel.CreateUniform());
		BGSearchEngine.EnableBackground(BGModel.CreateUniform());
		
	}
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		Options options = new Options();
		options.addOption("i", true, "input fasta file (or PBM format file with -pbm)");
		options.addOption("pwm", true, "input PWM file");
		options.addOption("prefix", true, "output directory");
		options.addOption("c", true, "control fasta file");
		options.addOption("N", true, "the number of PWM want to evaluate");
		options.addOption("bgmodel", true, "background model file");
		options.addOption("match", true, "find similar motifs in known PWM library (path to the library, e.g., jaspar.pwm)");
		options.addOption("convert", false, "convert input PWM file to the transfac format");
		options.addOption("roc", false, "compute AUC and draw ROC curve for the given pwm file");
		options.addOption("corr", false, "compute rank correlation between sequence and PWM score for the given pwm file");
		options.addOption("pbm", false, "input file is PBM format, and will compute signal correlation");
		options.addOption("markov", true, "use markov model of the control sequences rather than directly control sequences");
		options.addOption("maxROCfp",true, "the x axis scale in ROC curve drawing (default 1)");
		
		CommandLineParser parser = new GnuParser();
		MultPWMScorer evaluator=new MultPWMScorer();
		boolean rocflag=false;
		boolean dpwmflag=false;
		boolean corrflag=false;
		boolean genrand=false;
		boolean multiscore_flag=false;
		boolean convertflag=false;
		String inputPWM;
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
					//evaluator.removeBG=true;
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
			if(cmd.hasOption("corr"))
			{
				corrflag=true;
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
			
			if(cmd.hasOption("maxROCfp"))
			{
				evaluator.maxFPdraw=Double.parseDouble( cmd.getOptionValue("maxROCfp"));
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
		PWM.infothresh=0;  //without cutting the flanking position
		try {
			ArrayList<PWM> restorepwmlist=null;
			evaluator.initialize();
			pwmlist=common.LoadPWMFromFile(inputPWM);
			if(pwmlist.size()>topN)
				pwmlist= pwmlist.subList(0, topN);
			
			max_motif_span=0;
			Iterator<PWM> iter=pwmlist.iterator();
			ArrayList<DenseDoubleMatrix2D> ProfileCollection=new ArrayList<DenseDoubleMatrix2D>();
			ArrayList<DenseDoubleMatrix2D> BG_rofileCollection=new ArrayList<DenseDoubleMatrix2D>();
			while(iter.hasNext())
			{
				PWM motif=iter.next();
				max_motif_span+=motif.columns();
				DenseDoubleMatrix2D scoreProfiles=evaluator.getSinglePWMScoreProfile(motif,evaluator.SearchEngine);
				ProfileCollection.add(scoreProfiles);
				DenseDoubleMatrix2D scoreProfiles2=evaluator.getSinglePWMScoreProfile(motif,evaluator.BGSearchEngine);
				BG_rofileCollection.add(scoreProfiles2);
			}
			DenseDoubleMatrix2D finalProfile=evaluator.combineDifferentProfiles(ProfileCollection,BG_rofileCollection);
			
			
			//draw heatmap
			DrawUtil.drawHeatMap(finalProfile, evaluator.outputPrefix+"heatmap.png");
			//draw SignalAroundCenter
			DrawUtil.drawSignalAroundPeakCurve(finalProfile, evaluator.outputPrefix+"curve.png");
			//draw ROC
			DrawUtil.DrawROC(ROCdata, evaluator.outputPrefix+"ROC.png");
			
			//save profile to file
			 FileOutputStream fileOut =
   		         new FileOutputStream(evaluator.outputPrefix+"outputProfile.obj");
   		         ObjectOutputStream out =
   		                            new ObjectOutputStream(fileOut);
			out.writeObject(finalProfile);
	         out.close();
	          fileOut.close();
			
		}
		catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
	
	DenseDoubleMatrix2D combineDifferentProfiles(ArrayList<DenseDoubleMatrix2D> ProfileCollection,ArrayList<DenseDoubleMatrix2D> BG_rofileCollection)
	{
		if(!linearOrder)
		{
			//take max in the window size
			for (int i = 0; i < ProfileCollection.size(); i++) {
				ProfileCollection.set(i, maxWinSize(ProfileCollection.get(i),max_motif_span));
				BG_rofileCollection.set(i, maxWinSize(BG_rofileCollection.get(i),max_motif_span));
			}
		}
		
		int bestBarCode=(int) (Math.pow(2, ProfileCollection.size())-1);
		DenseDoubleMatrix2D bestProfile=null;
		double bestScore=Double.NEGATIVE_INFINITY;
		for (int i = 1; i < Math.pow(2, ProfileCollection.size()); i++) {
			int barcode=i;
			ArrayList<DenseDoubleMatrix2D> tempPColl=new ArrayList<DenseDoubleMatrix2D>();
			ArrayList<DenseDoubleMatrix2D> tempBGPColl=new ArrayList<DenseDoubleMatrix2D>();
			ArrayList<Integer> motiflen=new ArrayList<Integer>();
			for (int j = 0; j < ProfileCollection.size(); j++) {
				if(barcode%2==1)
				{
					tempPColl.add(ProfileCollection.get(j));
					tempBGPColl.add(BG_rofileCollection.get(j));
					motiflen.add(pwmlist.get(j).columns());
				}
				barcode>>=1;
			}
			DenseDoubleMatrix2D mixProfile=null;
			DenseDoubleMatrix2D mixBGProfile=null;
			if(linearOrder)
			{
				mixProfile=maxWinSize(linearOrderMix(tempPColl,motiflen),max_motif_span);
				mixBGProfile=maxWinSize(linearOrderMix(tempBGPColl,motiflen),max_motif_span);
			}
			else //normal add up
			{
				mixProfile=tempPColl.get(0);
				mixBGProfile=tempBGPColl.get(0);
				PlusMult Plus=PlusMult.plusMult(1);
				for (int j = 1; j < tempPColl.size(); j++) {
					mixProfile=(DenseDoubleMatrix2D) mixProfile.assign(tempPColl.get(j), Plus);
					mixBGProfile=(DenseDoubleMatrix2D) mixBGProfile.assign(tempBGPColl.get(j), Plus);
				}
			}
			double score=computeProfileScore(mixProfile,mixBGProfile,String.valueOf(i));
			if(score>bestScore)
			{
				bestProfile=mixProfile;
				bestBarCode=i;
				bestScore=score;
			}
		}
		String barcode_Str=Integer.toBinaryString(bestBarCode);
		System.out.println("Best combination: "+barcode_Str);
		return bestProfile;
	}
	
	
	static double computeProfileScore(DenseDoubleMatrix2D mixProfile, DenseDoubleMatrix2D mixBGProfile,String name)
	{
		double score=0;
		TreeMap<Double,Integer> Sorted_labels=new TreeMap<Double,Integer>();
		for (int i = 0; i < mixProfile.rows(); i++) {
			Sorted_labels.put(max(mixProfile.viewRow(i))-i*common.DoubleMinNormal, 1);
		}
		for (int i = 0; i < mixBGProfile.rows(); i++) {
			Sorted_labels.put(max(mixBGProfile.viewRow(i))+i*common.DoubleMinNormal, 0);
		}
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
		       	Confusion AUCcalc=AUCCalculator.readArrays(labels, scores);	
		       	XYSeries series1 = new XYSeries(name);
		       	int poscount=0;
		       	int skip=labels.length/50;
		       	for (int i = 0; i < labels.length; i++) {
					if(labels[i]==1)
						poscount++;
					double fp=(double)(i+1-poscount)/(labels.length-one);
					if(fp>maxFPdraw)
						break;
					if(i%skip==0)
					series1.add(fp, (double)(poscount)/one);
					
		       	}
		       	ROCdata.put(name, series1);
		       	
	          score=AUCcalc.calculateAUCROC();
		return score;
	}
	
	static DenseDoubleMatrix2D linearOrderMix(ArrayList<DenseDoubleMatrix2D> ProfileCollection,ArrayList<Integer> motiflen)
	{
		DenseDoubleMatrix2D retMat=(DenseDoubleMatrix2D) ProfileCollection.get(0).copy();
		DenseDoubleMatrix2D RretMat=(DenseDoubleMatrix2D) ProfileCollection.get(0).copy();
		DenseDoubleMatrix2D FretMat=(DenseDoubleMatrix2D) ProfileCollection.get(0).copy();
		int forwardGap=motiflen.get(0);
		int reverseGap=1;
		for (int i = 1; i < ProfileCollection.size(); i++) {
			DenseDoubleMatrix2D temp=ProfileCollection.get(i);
			for (int j = 0; j < temp.rows(); j++) {
				DoubleMatrix1D row = temp.viewRow(j);
				for (int j2 = 0; j2 < row.size(); j2++) {
					
					//forward strand
					double FVal=Double.NEGATIVE_INFINITY;
					if(j2+forwardGap<row.size())
					FVal=row.get(j2+forwardGap);
					FretMat.set(j, j2,FVal+FretMat.get(j, j2));
					//reverse strand
					double RVal=Double.NEGATIVE_INFINITY;
					if(j2-reverseGap>-1)
						RVal=row.get(j2-reverseGap);
					RretMat.set(j, j2,RVal+RretMat.get(j, j2));
				}
				
			}
			forwardGap+=motiflen.get(i);
			reverseGap+=motiflen.get(i);
		}
		
		if(ProfileCollection.size()>1)
			retMat=(DenseDoubleMatrix2D) FretMat.assign(RretMat,new maxFunction());

		return retMat;
	}
	
	static double max(DoubleMatrix1D row)
	{
		double maxV=Double.NEGATIVE_INFINITY;
		for (int i = 0; i < row.size(); i++) {
			double tmp=row.getQuick(i);
			if(maxV<tmp)
				maxV=tmp;
		}
		return maxV;
	}
	static DenseDoubleMatrix2D  maxWinSize(DenseDoubleMatrix2D mat,int winsize)
	{
		DenseDoubleMatrix2D ret=(DenseDoubleMatrix2D) mat.copy();
		for (int i = 0; i < mat.rows(); i++) {
			DoubleMatrix1D row = mat.viewRow(i);		
			for (int j = 0; j < mat.columns(); j++) {
				double maxScore=max(row.viewPart(Math.max(0,j-winsize), Math.min(mat.columns()-1,j+winsize)-Math.max(0,j-winsize)));
				ret.set(i, j, maxScore);
			}
		}
		return ret;
	}
	
	DenseDoubleMatrix2D getSinglePWMScoreProfile(PWM motif,LinearEngine engine)
	{
		int maxlen=0;
		for (int i = 0; i < engine.accSeqLen.size()-1; i++) {
			int len=engine.accSeqLen.get(i+1)-engine.accSeqLen.get(i);
			if(maxlen<len)
				maxlen=len;
		}
		DenseDoubleMatrix2D scoreProfile=new DenseDoubleMatrix2D(engine.ForwardStrand.size(), maxlen);
		 LinkedList<FastaLocation> falocs=engine.searchPattern(motif, Double.NEGATIVE_INFINITY);
		 Iterator<FastaLocation> iter=falocs.iterator();
		 while(iter.hasNext())
    	 {
    		 FastaLocation currloc=iter.next();
    		 scoreProfile.set(currloc.getSeqId(),currloc.getSeqPos(),currloc.Score);
    	 }
		return scoreProfile;
		
	}

}
