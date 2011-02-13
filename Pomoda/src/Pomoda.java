/**
 * @author zhizhuo zhang
 * zzz2010@gmail.com
 */
import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.biojava.bio.dp.ScoreType;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;



public class Pomoda {





	public String outputPrefix="./";
	public String inputFasta;
	public int seedlen=5;
	public int resolution=20;
	public int starting_windowsize=200;
	public int ending_windowsize=600;
	public double FDR=0.01;
	public int max_motiflen=15;
	public int num_motif=5;
	public double sampling_ratio=0.8;
	public double min_support_ratio=0.05;
	public boolean markflag=true;
	public HashEngine SearchEngine;
	public BGModel background;
	public ArrayList<Double> pos_prior;
	
	public void initialize()
	{
		//build hash index
		SearchEngine=new HashEngine(5);
		SearchEngine.build_index(this.inputFasta);
//		if(SearchEngine_Test())
//			System.out.println("SearchEngine_Test : pass");
			
		
		background=new BGModel();
		File file = new File(inputFasta+".bgobj");
		if(file.exists())
		{
			background.LoadModel(inputFasta+".bgobj");
			background.order=4;
		}
		else
		{
		     background.BuildModel(inputFasta, 4); //3-order bg
		     background.SaveModel(inputFasta+".bgobj");
		}
//		if(BGModel_Test())
//			System.out.println("BGModel_Test : pass");
		
//		if(PWM_Test())
//			System.out.println("PWM_Test : pass");
			
		pos_prior=new ArrayList<Double>(SearchEngine.getTotalLength()/SearchEngine.getSeqNum()/this.resolution);
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
			PWM testpwm=new PWM(new String[]{"ATTAAA","AATTAA","AAATTA","AAAATT","AAAAAT"});
			BGModel bgtest=new BGModel();
			bgtest.BuildModel(new String[]{"ATTATT","ATTATT","ATTATT","ATTATT","ATTATT","ATTATT","ATTATT","ATTATT","ATTATT","ATTATT","ATTATT"}, 4);
			String testpatt="ATTATT";
			double score1=Math.exp(testpwm.scoreWeightMatrix(testpatt, ScoreType.PROBABILITY));
			double score2=Math.exp(bgtest.Get_LOGPROB(testpatt));
			if(score2<score1)
			{
				pass=false;
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
	
	
	public ArrayList<PWM>	 getSeedMotifs() {
		
		ArrayList<PWM>	SeedMotifs=new	ArrayList<PWM>(num_motif*num_motif);
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
			if(pos_prior.size()==0)
				score=positionlist.size()*(-0.037267253272904234-logprob_bg);//sum loglik ,-0.037267253272904234 is from pseudo count
			else
				score=CenterDistributionScore(LocList,0,logprob_bg);
			seedScores.put(hash, score);
		}
		//sort by score
		List<Map.Entry<Integer,Double>> mappingList=new ArrayList<Map.Entry<Integer,Double>>(seedScores.entrySet());
		Collections.sort(mappingList, bvc);
		int max_num_Seeds=num_motif*num_motif;
		//just take top ones
		    for (Map.Entry<Integer,Double> pair : mappingList) {
		    	if(SeedMotifs.size()==max_num_Seeds)
		    		break;
			    String toppattern=common.Hash2ACGT(pair.getKey(),seedlen);
			    //System.out.println(pair.getValue());
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
	
	public PWM Column_Replacement(PWM motif)
	{
		double log025=Math.log(0.25);
		
		double bestscore=motif.Score;
		
		int num_priorbin=SearchEngine.getTotalLength()/SearchEngine.getSeqNum()/this.resolution;
		do
		{
			String consensus_core=motif.Censensus(true);
		System.out.println(consensus_core+"\t"+String.valueOf(bestscore));
		
		double[] temp_prior=new double[num_priorbin];
		double lognullprior=Math.log(1.0/num_priorbin);
		LinkedList<Map.Entry<Double,String>> inst=motif.GenerateInstanceFromPWMPQ(this.sampling_ratio, this.FDR, this.background);
		if(inst.size()==0)
			return motif;
		
		Iterator<Map.Entry<Double,String>> iter=inst.iterator();
		double [][] loglik_matrix=new double[motif.columns()][4];
		String consensus=motif.Censensus(false);
		//double logN025=log025*(motif.head+motif.tail);
		int motiflen=consensus_core.length();
		//update the loglik matrix
		while(iter.hasNext())
		{
			Map.Entry<Double,String> patternEntry=iter.next();
			double logprob_theta=patternEntry.getKey();
			LinkedList<Integer> locs=SearchEngine.searchPattern(patternEntry.getValue(), 0);

				LinkedList<FastaLocation> Falocs=SearchEngine.Int2Location(locs);
				Iterator<FastaLocation> iter2=Falocs.iterator();
				int count=0;
				while(iter2.hasNext())
				{
					FastaLocation currloc=iter2.next();
					String site="";
					//forward site
					site=SearchEngine.getSite(currloc.getMin()-motif.head, motif.columns());
					if(site.equalsIgnoreCase(""))
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
					double logprob_BG=background.Get_LOGPROB(site.replace("N", ""));
					if(count>=SearchEngine.forwardCount)
					{
						//reverse site
						site=common.getReverseCompletementString(site);
						
					}	
					double logprior=0;
					if(pos_prior.size()!=0)
						logprior=pos_prior.get( pos_prior.size()*(currloc.getSeqPos()+motiflen/2)%currloc.getSeqLen()/currloc.getSeqLen())-lognullprior;
					double loglik=logprob_theta+logprior-logprob_BG;
//					if(loglik<0)
//						System.out.println(site);
					temp_prior[num_priorbin*(currloc.getSeqPos()+motiflen/2)%currloc.getSeqLen()/currloc.getSeqLen()]+=loglik/10000;//make smaller
					for (int i = 0; i < site.length(); i++) {
						if(consensus.charAt(i)!='N')
							continue;
						int symid=common.acgt(site.charAt(i));
						if(symid>3)
							continue; //meet new line separator
						loglik_matrix[i][symid]+=loglik;
					}
			}
	
		}
		
		//select the best column replacement
		double maxloglik=Double.MIN_VALUE;
		
		int numbestSym=2;
		int bestCol=-1;
		ArrayList<Integer> bestSym=new ArrayList<Integer>(numbestSym);
			for (int i = 0; i < motif.columns(); i++) {
				if(consensus.charAt(i)!='N')
					continue;
				TreeMap<Double,Integer> orderSym=new TreeMap<Double,Integer>();
				double temploglik=0;
				double sumtemp=0;
				for (int j = 0; j < 4; j++) {
					double temp=loglik_matrix[i][j];
					orderSym.put(temp-Double.MIN_NORMAL*j, j);
					sumtemp+=temp;
				}
				//System.out.println(sumtemp);
				//only consider best two sym
				
				for(Double key: orderSym.descendingKeySet()) {
					//int col=orderSym.get(key);
					if(key<0)
						break;
					temploglik+=key;
					
				}
				
				if(temploglik>maxloglik)
				{
					maxloglik=temploglik;
					bestCol=i;
					bestSym.clear();
					for(Double key: orderSym.descendingKeySet()) {
						int col=orderSym.get(key);
						if(key<0||bestSym.size()==numbestSym)
							break;
						bestSym.add(col);
						
					}
					
					
				}
				
			}
			if(bestCol==-1)
				break;
			double [] repColumnValue=new double[4];
			double sumNorm=0;
			for (int j = 0; j < bestSym.size(); j++) {
				double temp=loglik_matrix[bestCol][bestSym.get(j)];
				if(temp>0)
				{
					repColumnValue[j]=(temp);
					sumNorm+=repColumnValue[j];
				}
				
			}
			//normalizing the column
			for (int j = 0; j < 4; j++) {				
					repColumnValue[j]/=sumNorm;				
			}
			
			
			System.out.println(sumNorm);
			if(sumNorm<=bestscore)
				break;
			else
			{
				bestscore=sumNorm;
				
			}


			//update motif column value
			motif.setWeights(bestCol, repColumnValue);
			double sumPrior=0;
			for (int i = 0; i < temp_prior.length; i++) {
				sumPrior+=temp_prior[i];
			}
			pos_prior.clear();
			for (int i = 0; i < temp_prior.length; i++) {
				temp_prior[i]/=sumPrior;
				pos_prior.add(temp_prior[i]);
			}
			
			
		}while(true);
		
		return motif;
	}
	
	//for fix PWM score and BG score, but there is position prior
	double CenterDistributionScore(LinkedList<FastaLocation> LocList,double pwmloglik,double bgloglik)
	{
		double score=0;
		
       Iterator<FastaLocation> iter=LocList.iterator();
		int num_priorbin=SearchEngine.getTotalLength()/SearchEngine.getSeqNum()/this.resolution;
		double lognullprior=Math.log(1.0/num_priorbin);
		while(iter.hasNext())
		{
			FastaLocation currloc=iter.next();
			double logprior=0;
			if(pos_prior.size()!=0)
				logprior=pos_prior.get( pos_prior.size()*currloc.getSeqPos()/currloc.getSeqLen())-lognullprior;
			double loglik=pwmloglik+logprior-bgloglik;
			score+=loglik;
		}
		
		return score;
	}
	
	//for fix PWM score but there is 'N' inside the pattern
	double CenterDistributionScore(LinkedList<FastaLocation> LocList,double pwmloglik,int motiflen)
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
		
		//get seed motifs
		ArrayList<PWM>  seedPWMs=motifFinder.getSeedMotifs();
		
		//extend and refine motifs
		for (int i = 0; i < seedPWMs.size(); i++) {
			seedPWMs.set(i,motifFinder.Column_Replacement(seedPWMs.get(i)));
			if(motifFinder.markflag)
			{
				//do something to mark the locations in SearchEngine
				
			}
		}
		

		

	}

}
