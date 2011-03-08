/**
 * @author zhizhuo zhang
 * zzz2010@gmail.com
 */
import java.io.BufferedWriter;
import java.io.File;
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
import java.util.Map;
import java.util.TreeMap;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
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
import org.pr.clustering.hierarchical.LinkageCriterion;

import cern.jet.random.Binomial;
import cern.jet.random.engine.RandomEngine;





public class Pomoda {





	public String outputPrefix="./";
	public String inputFasta;
	public String ctrlFasta="";
	public String bgmodelFile="";
	public int seedlen=5;
	public boolean debug=false;
	public boolean OOPS=false;
	public int resolution=10;
	public int starting_windowsize=200;
	public int ending_windowsize=600;
	public double FDR=0.001;
	public int max_motiflen=25;
	public int num_motif=5;
	public double sampling_ratio=1;
	public double min_support_ratio=0.1;
	public boolean maskflag=true;
	public HashEngine SearchEngine;
	public LinearEngine SearchEngine2;
	public BGModel background;
	public String linkage="WPGMA";
	public ArrayList<Double> pos_prior;
	
	public void initialize()
	{
		//build hash index
		SearchEngine=new HashEngine(5);
		SearchEngine.build_index(this.inputFasta);
		SearchEngine2=new LinearEngine(4);
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
//			for (int i = 0; i < AR.size(); i++) {
//			//	AR.get(i).print();
//
//			     this.Relax_Seed_(AR.get(i));
//			}
			
			
//			LinkedList<String> sites=new LinkedList<String>();
//			String gappattern="ACANNNTGT";
//			int Nlen=3;
//			int Nstart=3;
//			int gapmerSize=1<<(Nlen*2);
//			double[] gapmerCount=new double[gapmerSize];
//			LinkedList<FastaLocation> falocs=SearchEngine2.searchPattern(gappattern, 0);
//			Iterator<FastaLocation> iter=falocs.iterator();
//			while(iter.hasNext())
//			{
//				FastaLocation currloc=iter.next();
//				String site=SearchEngine2.getSite(currloc.getSeqId(), currloc.getSeqPos(), gappattern.length());
//				if(currloc.ReverseStrand)
//					site=common.getReverseCompletementString(site);
//				int hash=common.getHashing(site, Nstart, Nlen);
//				gapmerCount[hash]+=1;
//				sites.add(site.substring(Nstart, Nstart+Nlen));		
//			}
//			PWM model=new PWM(sites.toArray(new String[1]));
//			for (int i = 0; i < gapmerSize; i++) {
//				String gapstr=common.Hash2ACGT(i, Nlen);
//				double log_p=model.scoreWeightMatrix(gapstr,ScoreType.PROBABILITY);
//				double fold=gapmerCount[i]/(falocs.size()*Math.exp(log_p));
//				System.out.println(gapstr+'\t'+fold+'\t'+gapmerCount[i]/falocs.size());
//			}
			
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
				score=CenterDistributionScore(LocList,-common.DoubleMinNormal*seedlen,logprob_bg);
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
	

	

	
	public PWM Relax_Seed_(PWM motif)
	{
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
		double bestscore=motif.Score;
		double lastscore=1;
		int num_priorbin=SearchEngine.getTotalLength()/SearchEngine.getSeqNum()/this.resolution;
	
		HashSet<Integer> stateCodes=new HashSet<Integer>();
	
		int flankingLen=2;
		if(motif.head+flankingLen+motif.core_motiflen>=motif.columns()||motif.head-flankingLen<0)
			flankingLen=0;
		do
		{
			String consensus_core=motif.Consensus(true);
			int motiflen=consensus_core.length(); 
			System.out.println(consensus_core+"\t"+String.valueOf(bestscore));
			bestscore=0;
			double[] temp_prior=new double[num_priorbin];
			double lognullprior=Math.log(1.0/num_priorbin);
	
			
			//double log_thresh=motif.getThresh(sampling_ratio, FDR, background);
			double log_thresh=0;
			SearchEngine2.EnableBackground(background);
			double [][]m_matrix=new double [motiflen+flankingLen*2][4];
			
			//update the loglik matrix
			LinkedList<FastaLocation> Falocs=SearchEngine2.searchPattern(motif, log_thresh);
			System.out.println("number of occurrences: "+String.valueOf(Falocs.size()));
					Iterator<FastaLocation> iter2=Falocs.iterator();
					int count=0;
					int match_seqCount=0;
					int lastseq=-1;
					double max_seqloglik=0;
					double max_seqprob_theta=0;
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
						//double loglik=logprob_theta+logprior-logprob_BG;
						double loglik=logprob_theta+logprior;
						double prob_theta=Math.exp(logprior)/10000;//Math.exp(currloc.Score);
						temp_prior[prior_bin]+=prob_theta;//make smaller
						if(OOPS)
							loglik-=common.DoubleMinNormal*Math.abs(currloc.getSeqLen()/2-currloc.getSeqPos()-motiflen/2); //add small bias to center
						if(!OOPS)
							bestscore+=loglik;
						if(OOPS&&currloc.getSeqId()!=lastseq)
						{
							match_seqCount++;
							bestscore+=max_seqloglik;
							lastseq=currloc.getSeqId();
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
						
						if(max_seqloglik<loglik)
						{
							max_seqloglik=loglik;
							max_seqprob_theta=prob_theta;
							max_seqsite=site;
							if(prob_theta>1000000)
								continue;
						}
						
						if(!OOPS)
						for (int i = 0; i < motiflen+flankingLen*2; i++) {
							int symid=common.acgt(site.charAt(i));
							if(symid>3)
								continue;
							m_matrix[i][symid]+=prob_theta;//prob_theta;
						}
					}
					
					
					//last instance for OOPS
					if(OOPS)
					{
						match_seqCount++;
						bestscore+=max_seqloglik;
						
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
						motif.Score=bestscore;
						lastscore=bestscore;
						}
						
							
							
						
							
						
						if(Falocs.size()==0|| stateCodes.contains(Falocs.hashCode()))
							return motif;
						else
							stateCodes.add(Falocs.hashCode());
						
						//motif=new PWM((String[])(MatchSite.toArray(new  String[1])));

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
			}while(motif.core_motiflen<max_motiflen);
			

		
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
		
		int inst_hash=-1;
		do
		{
		
			String consensus_core=motif.Consensus(true);
		
		System.out.println(consensus_core+"\t"+String.valueOf(bestscore));
		bestscore=0;
		double[] temp_prior=new double[num_priorbin];
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
					
					double logprior=0;
					if(motif.pos_prior.size()!=0)
						logprior=motif.pos_prior.get((int)( motif.pos_prior.size()*((currloc.getSeqPos()+motiflen/2)%currloc.getSeqLen()/(double)currloc.getSeqLen())))-lognullprior;
					double loglik=logprob_theta+logprior-logprob_BG;
					temp_prior[(int)(num_priorbin*((currloc.getSeqPos()+motiflen/2)%currloc.getSeqLen()/(double)currloc.getSeqLen()))]+=Math.exp(loglik)/10000;//make smaller
					if(OOPS)
						loglik-=common.DoubleMinNormal*Math.abs(currloc.getSeqLen()/2-currloc.getSeqPos()-motiflen/2); //add small bias to center
					if(!OOPS)
						bestscore+=loglik;
					if(OOPS&&currloc.getSeqId()!=lastseq)
					{
						if(lastseq!=-1)
							for (int i = 0; i < site.length(); i++) {
//								if(consensus.charAt(i)!='N')
//									continue;
								
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
				if(motif.head<i&&i<motiflen+motif.head)
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
			if(motif.head<bestCol&&bestCol<motiflen+motif.head)
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
			if(debug)
				motif.print();
			
			
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
		double max_seqloglik=Double.MIN_VALUE;
		int lastseq=-1;
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
				max_seqloglik=Double.MIN_VALUE;
			}
			if(OOPS&&loglik>max_seqloglik)
			{
				max_seqloglik=loglik;
				continue;
			}
			
			score+=loglik;
		}
		if(OOPS)
			score+=max_seqloglik;
		
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
		options.addOption("c", true, "control fasta file");
		options.addOption("bgmodel", true, "background model file");
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
		options.addOption("FDR",true,"fasle positive rate");
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
			formatter.printHelp( "Pomoda", options );
			return;
		}
		
		//initialize Pomoda 
		motifFinder.initialize();

		
		//get seed motifs
		ArrayList<PWM>  seedPWMs=motifFinder.getSeedMotifs();
		double topseed_Score=seedPWMs.get(0).Score;
		topseed_Score=0;//zzz
		File file = new File(motifFinder.outputPrefix+"pomoda_raw.pwm"); 
		try {
			BufferedWriter writer = new BufferedWriter(new FileWriter(file));
			TreeMap<Double, PWM> sortedPWMs=new TreeMap<Double, PWM>();
		//extend and refine motifs
		for (int i = 0; i < seedPWMs.size(); i++) {
			System.out.println("Extending...");
			seedPWMs.set(i,motifFinder.Column_Replacement_(seedPWMs.get(i)));
			System.out.println("Relaxing...");
			seedPWMs.set(i, motifFinder.Relax_Seed_(seedPWMs.get(i)));

			seedPWMs.get(i).Name="Motif"+String.valueOf(i+1);
			// to make different length comparable ,need to consider the instance coverage
//			seedPWMs.get(i).Score=seedPWMs.get(i).Score/seedPWMs.get(i).inst_coverage;//corrected score
			if(motifFinder.maskflag&&seedPWMs.get(i).Score>(topseed_Score))
			{
				//do something to mark the locations in SearchEngine
				System.out.println("Masking...");
				motifFinder.Masking(seedPWMs.get(i));
				topseed_Score=seedPWMs.get(i).Score;
				
			}
			//Runtime.getRuntime().gc();
			writer.write(seedPWMs.get(i).toString());
		
		
		}

		writer.close();
		
		file = new File(motifFinder.outputPrefix+"pomoda_clust.pwm"); 
		writer = new BufferedWriter(new FileWriter(file));
		//clustering motif, re-initialize
		System.out.println("Clustering...");
		motifFinder.SearchEngine2.build_index(motifFinder.inputFasta);
		sortedPWMs.clear();
		PWMcluster clustering=new PWMcluster(motifFinder);
		
			ArrayList<PWM>  clusterPWMs=clustering.Clustering(seedPWMs,motifFinder.num_motif);
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
				
//				GapPWM gpwm=gimprover.fillDependency(sortedPWMs.get(key));
//				
//				System.out.println("PWM AUC:"+gimprover.AUCtest(sortedPWMs.get(key)));
//				System.out.println("Gap improved PWM AUC:"+gimprover.AUCtest(gpwm));
			}
			
			writer.close();
		
		} 
		catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		}
		

		

	

}
