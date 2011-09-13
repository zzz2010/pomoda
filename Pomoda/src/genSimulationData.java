import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Iterator;
import java.util.LinkedList;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.biojava.bio.dist.DistributionTools;


public class genSimulationData {
	
	static boolean PBMFlag=false;
	public static void generateSimulatedData(LinkedList<PWM> pwmlist, int N, int len, String outputFasta, BGModel background, double sampling_ratio)
	{
		File file2 = new File(outputFasta); 
		File file1 = new File(outputFasta.replace(".fa", ".ans").replace("_PBM.txt", "_PBM.ans")); 
		try {
			BufferedWriter writer2= new BufferedWriter(new FileWriter(file2));
			BufferedWriter writer1= new BufferedWriter(new FileWriter(file1));
	         background.r.setSeed(common.randomseed);

			for (int i = 0; i < N; i++) {
	        	 String bgstr="";
		        	 KeyValuePair<Double, String> bgstr_p=background.generateRandomSequence(len);
		        	 bgstr=bgstr_p.value;

	        	 Iterator<PWM> iter=pwmlist.iterator();
	        	 double bestscore=-100000;
	        	 while(iter.hasNext())
	        	 {
	        		 boolean depflag=false;
	        		 PWM motif=iter.next();
	        		 if( motif instanceof GapPWM)
	        			 depflag=true;
	        		 
	        		 if(background.r.nextDouble()<sampling_ratio)
	        		 {
	        			 String site="";
	        			 site=motif.get_randomSite().toLowerCase();
	        			 double score=motif.scoreWeightMatrix(site);
	        			 if(score>bestscore)
	        				 bestscore=score;
	        			 //no reverse complement
//	        			 //reverse strand
//	        			 if(background.r.nextDouble()<0.5)
//	        				 site=common.getReverseCompletementString(site);
	        			 
	        			 int pos=background.r.nextInt(len-motif.columns());
	        			 bgstr=bgstr.substring(0, pos)+site+bgstr.substring(pos+site.length());
	        			 writer1.write(motif.Name+"\t"+site+"\t"+i+"\t"+pos+"\n");
	        		 }
	        	 }
	        	 if(!PBMFlag)
	        	 {
	        		 writer2.write(">"+i+"\n");
	        		 writer2.write(bgstr+"\n");
	        	 }
	        	 else
	        	 {
	        		 writer2.write(Math.exp(bestscore)*10000+"\t"+bgstr+"\n");
	        	 }
	        	 
			}
			writer2.close();
			writer1.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		int topN=2;
		Options options = new Options();
		options.addOption("i", true, "input fasta file");
		options.addOption("pwm", true, "input PWM file");
		options.addOption("pbm", false, "output PBM format file");
		options.addOption("N", true, "the number of PWMs used to mask");
		options.addOption("c", true, "control fasta file");
		options.addOption("bgmodel", true, "background model file");
		options.addOption("prefix", true, "output directory");
		options.addOption("FDR",true,"fasle positive rate");
		options.addOption("simlen",true,"turn on the simulation function: -N number of sequences, -FDR the percentage contain motif");
		String inputPWM;
		int simlen=-1;
		CommandLineParser parser = new GnuParser();
		PWMevaluator dumy=new PWMevaluator();
		try {
			CommandLine cmd = parser.parse( options, args);

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
				dumy.ctrlFasta=cmd.getOptionValue("c");
			}
			if(cmd.hasOption("pbm"))
			{
				PBMFlag=true;
			}
			if(cmd.hasOption("N"))
			{
				topN=Integer.parseInt(cmd.getOptionValue("N"));
			}
			if(cmd.hasOption("simlen"))
			{
				simlen=Integer.parseInt(cmd.getOptionValue("simlen"));
			}
			if(cmd.hasOption("bgmodel"))
			{
				dumy.bgmodelFile=cmd.getOptionValue("bgmodel");
			}
			if(cmd.hasOption("prefix"))
			{
				dumy.outputPrefix=cmd.getOptionValue("prefix");
			}
			if(cmd.hasOption("FDR"))
			{
				dumy.FDR=Double.parseDouble(cmd.getOptionValue("FDR"));
			}
		} catch (ParseException e) {
			// TODO Auto-generated catch block
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp( "FastaMask", options );
			return;
		}
		
		if(dumy.background==null)
		{
			dumy.background=new BGModel();
			File file=null;
			int bg_markov_order=0;
			if(!dumy.ctrlFasta.isEmpty())
			{
				dumy.BGSearchEngine=new LinearEngine(4);
				dumy.BGSearchEngine.build_index(dumy.ctrlFasta);
				dumy.removeBG=true;
			}
			
			if(!dumy.bgmodelFile.isEmpty())
			{
				dumy.removeBG=true;
					dumy.background.LoadModel(dumy.bgmodelFile);
			}

		}
		LinkedList<PWM> pwmlist=common.LoadPWMFromFile(inputPWM);
		String basename= (new File(inputPWM)).getName();
		if(simlen>0)
		{// turn on simulated data generation, turn off the fasta masking
			String outputFasta=dumy.outputPrefix+"/"+basename.split("\\.")[0]+topN+"_"+simlen+"_"+dumy.FDR+".fa";
			if(PBMFlag)
			{
				outputFasta=outputFasta.replace(".fa", "_PBM.txt");
			}
			generateSimulatedData(pwmlist, topN, simlen, outputFasta, dumy.background, dumy.FDR);
			return ;
		}
		
		for (int i = 0; i < Math.min(topN, pwmlist.size()); i++) {
			PWM motif=pwmlist.get(i);
			double log_thresh=motif.getThresh(1,dumy.FDR, dumy.background,true);
			String consensus_core=motif.Consensus(true);
			
			
			//double logN025=log025*(motif.head+motif.tail);
			int motiflen=consensus_core.length();
			LinkedList<String> MatchSite=new LinkedList<String>();
			
			//update the loglik matrix
			dumy.SearchEngine.DisableBackground();
			LinkedList<FastaLocation> Falocs=dumy.SearchEngine.searchPattern(motif, log_thresh);
			Iterator<FastaLocation> iter2=Falocs.iterator();
			int count=0;
			int lastseq=-1;
			
			String X_Str="";
			for (int j = 0; j < motiflen; j++) {
				X_Str+="N";
			}
			while(iter2.hasNext())
			{
				FastaLocation currloc=iter2.next();
				String Site1=dumy.SearchEngine.getSite(currloc.getSeqId(), currloc.getSeqPos(), motiflen);		
				String rep=dumy.SearchEngine.ForwardStrand.get(currloc.getSeqId()).replace(Site1, X_Str);
				dumy.SearchEngine.ForwardStrand.set(currloc.getSeqId(), rep);
			}
		
		}
		
		File file = new File(dumy.inputFasta+"_"+basename.substring(0,3)+"mask.fa"); 
		try {
			BufferedWriter writer = new BufferedWriter(new FileWriter(file));
			Iterator<String> iter=dumy.SearchEngine.ForwardStrand.iterator();
			int ii=0;
			while(iter.hasNext())
			{
				writer.write(">mask"+ii+"\n");
				ii++;
				writer.write(iter.next().toUpperCase()+"\n");
			}
			writer.close();
		}
		catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		

	}




}
