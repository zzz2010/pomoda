import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.TreeMap;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;


public class FastaMask {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		int topN=2;
		Options options = new Options();
		options.addOption("i", true, "input fasta file");
		options.addOption("pwm", true, "input PWM file");
		options.addOption("N", true, "the number of PWMs used to mask");
		options.addOption("c", true, "control fasta file");
		options.addOption("bgmodel", true, "background model file");
		options.addOption("prefix", true, "output directory");
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
			if(cmd.hasOption("N"))
			{
				topN=Integer.parseInt(cmd.getOptionValue("N"));
			}
			if(cmd.hasOption("bgmodel"))
			{
				GImprover.bgmodelFile=cmd.getOptionValue("bgmodel");
			}
			if(cmd.hasOption("prefix"))
			{
				GImprover.outputPrefix=cmd.getOptionValue("prefix");
			}
			if(cmd.hasOption("FDR"))
			{
				GImprover.FDR=Double.parseDouble(cmd.getOptionValue("FDR"));
			}
		} catch (ParseException e) {
			// TODO Auto-generated catch block
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp( "FastaMask", options );
			return;
		}
		
		GImprover.initialize();
		LinkedList<PWM> pwmlist=common.LoadPWMFromFile(inputPWM);
		for (int i = 0; i < Math.min(topN, pwmlist.size()); i++) {
			PWM motif=pwmlist.get(i);
			double log_thresh=motif.getThresh(1,GImprover.FDR, GImprover.background);
			String consensus_core=motif.Consensus(true);
			
			
			//double logN025=log025*(motif.head+motif.tail);
			int motiflen=consensus_core.length();
			LinkedList<String> MatchSite=new LinkedList<String>();
			
			//update the loglik matrix
			GImprover.SearchEngine.DisableBackground();
			LinkedList<FastaLocation> Falocs=GImprover.SearchEngine.searchPattern(motif, log_thresh);
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
				String Site1=GImprover.SearchEngine.getSite(currloc.getSeqId(), currloc.getSeqPos(), motiflen);		
				String rep=GImprover.SearchEngine.ForwardStrand.get(currloc.getSeqId()).replace(Site1, X_Str);
				GImprover.SearchEngine.ForwardStrand.set(currloc.getSeqId(), rep);
			}
		
		}
		String basename= (new File(inputPWM)).getName();
		File file = new File(GImprover.inputFasta+"_"+basename.substring(0,3)+"mask.fa"); 
		try {
			BufferedWriter writer = new BufferedWriter(new FileWriter(file));
			Iterator<String> iter=GImprover.SearchEngine.ForwardStrand.iterator();
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
