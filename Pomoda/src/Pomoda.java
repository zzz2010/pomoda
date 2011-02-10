import java.util.Iterator;
import java.util.LinkedList;

import org.apache.commons.cli.*;



public class Pomoda {





	public String outputPrefix="./";
	public String inputFasta;
	public int seedlen=5;
	public int resolution=40;
	public int starting_windowsize=200;
	public int ending_windowsize=600;
	public double FDR=0.01;
	public int max_motiflen=30;
	public int num_motif=5;
	public double sampling_ratio=0.8;
	public double min_support_ratio=0.05;
	public boolean markflag=true;
	public ISearchEngine SearchEngine;
	
	public void initialize()
	{
		//build hash index
		SearchEngine=new HashEngine(5);
		SearchEngine.build_index(this.inputFasta);
		SearchEngine_Test();
	}
	
	private boolean SearchEngine_Test()
	{
		boolean pass=true;
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
				System.out.print(sitestring);
				break;
			}
			
		}
		
		return pass;
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// read parameters from command line
		Options options = new Options();
		options.addOption("i", true, "input fasta file");
		options.addOption("prefix", true, "output directory");
		options.addOption("seedlen", true, "kmer seed motif length (default 5)");
		options.addOption("ratio",true, "sampling ratio (default 0.8)");
		options.addOption("n",true, "number of motifs in final report (default 5)");
		options.addOption("rs",true, "counting resolution (default 40 bp)");
		options.addOption("supp",true, "minimum support ratio, the percentage of peaks contains motif (default 0.05)");
		options.addOption("maxlen",true,"maxmimum length of the motif (default 30)");
		options.addOption("minw",true,"minimum size of motif binding region (default 200bp)");
		options.addOption("maxw",true, "maximum size of motif binding region (default 600bp)");
		options.addOption("mark",true,"whether marking the top motif location in order to find co-motif");
		options.addOption("FDR",true,"fasle positive rate");
		
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
			if(cmd.hasOption("mark"))
			{
				motifFinder.markflag=true;
			}
			else
			{
				motifFinder.markflag=false;
			}
			if(cmd.hasOption("FDR"))
			{
				motifFinder.FDR=Double.parseDouble(cmd.getOptionValue("FDR"));
			}
	
			
		} catch (ParseException e) {
			// TODO Auto-generated catch block
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp( "Pomoda", options );
			return;
		}
		
		//initialize Pomoda 
		motifFinder.initialize();

		

	}

}
