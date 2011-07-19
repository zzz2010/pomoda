import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.TreeMap;

import org.biojava.bio.BioException;
import org.biojava.bio.dp.IllegalTransitionException;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.seq.io.SymbolTokenization;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojavax.bio.seq.RichSequence;


public class LinearEngine {

	LinkedList<String> ForwardStrand;
	//LinkedList<String> ReverseStrand;
	int num_thread;
	BGModel background=null;
	public HashMap<String,LinkedList<FastaLocation>> KmerHitList=null;
	public HashMap<Integer,HashMap<Integer,ArrayList<Double>>> BGscoreMap=null;
    public int forwardCount; //the first n entry in the searchPattern is forward.
	public int TotalLen=0;
	public ArrayList<Integer> accSeqLen;
	public LinearEngine(int num_thread) {
		if(num_thread<1)
			num_thread=1;
		this.num_thread=num_thread;
		ForwardStrand=new LinkedList<String>();
		//ReverseStrand=new LinkedList<String>();
		accSeqLen=new ArrayList<Integer>(10000);
	}
	
	
	public void EnableBackground(BGModel bg)
	{
		background=bg;
		BGscoreMap=new HashMap<Integer,HashMap<Integer,ArrayList<Double>>>();
	}
	public void DisableBackground()
	{
		background=null;
	}
	
	public void build_KmerHitList(int kmerlen)
	{
		KmerHitList=new HashMap<String, LinkedList<FastaLocation>>((int)Math.pow(4, kmerlen));
		Iterator<String> iter=ForwardStrand.iterator();
		int seqid=0;
		int pos=0;
		while(iter.hasNext())
		{
			String seq=iter.next();
			for (int i = 0; i < seq.length()-kmerlen; i++) {
				String kmer=seq.substring(i, i+kmerlen).toUpperCase();
				
				
				if(!KmerHitList.containsKey(kmer))
				{
					KmerHitList.put(kmer, new LinkedList<FastaLocation>());
				}
				FastaLocation fa=new FastaLocation(pos+i, seqid, i, seq.length());
				KmerHitList.get(kmer).add(fa);
			
				
			
			}
			seqid++;
			pos+=seq.length();
		}
	}
	public void build_index(String inputfile,int maxSeq) {
		accSeqLen.clear();
		// TODO Auto-generated method stub
		 try {
			 ForwardStrand.clear();
			// ReverseStrand.clear();
			 //Database to hold the training set
			      BufferedReader br = new BufferedReader(new FileReader(inputfile));
   		          SymbolTokenization toke = AlphabetManager.alphabetForName("DNA").getTokenization("token");
   		          SequenceIterator seqi = RichSequence.IOTools.readFasta(br, toke,null);
   		          
   		          TotalLen=0;
   		       accSeqLen.add(0);
   			      while (seqi.hasNext()) {
   			    	  
   			    	  Sequence seq=seqi.nextSequence();
   			    	  String seqstr=seq.seqString().replace("N", "");
   				     ForwardStrand.add(seqstr);
   				  //ReverseStrand.add(common.getReverseCompletementString( seqstr));
   				if(ForwardStrand.size()>maxSeq)
   				    	 return;
   				TotalLen+=seqstr.length();
   				accSeqLen.add(TotalLen);
   			      }
   	
			     
			
		} catch (IllegalSymbolException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IllegalTransitionException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IllegalAlphabetException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IllegalArgumentException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (BioException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}
	
	public void build_index(String inputfile, LinkedList<String> FilteredSeqs)
	{
		int maxctrlseq=Math.max(1000,FilteredSeqs.size()*2);
		
		int idlen=100;
		HashSet<Integer> filtercode=new HashSet<Integer>(FilteredSeqs.size());
		for(String seq :FilteredSeqs)
		{
			for (int i = 0; i < seq.length()-idlen; i++) {
				filtercode.add(seq.substring(i, i+idlen).hashCode());
			}
		}
		 try {
			 ForwardStrand.clear();
			 accSeqLen.clear();
			// ReverseStrand.clear();
			 //Database to hold the training set
			      BufferedReader br = new BufferedReader(new FileReader(inputfile));
			      int filternum=0;
			      String Line="";
			      String seqstr="";
			      TotalLen=0;
	   		       accSeqLen.add(0);
			     while( (Line=br.readLine())!=null)
			     {
			    	 if(Line.length()>0)
			    	 {
			    		 if(Line.charAt(0)=='>')
			    		 {
			    			 if(seqstr!="")
			    			 {
			    			 seqstr=seqstr.replace("N", "");
			    			 boolean filterflag=false;
			    			 for (int i = 0; i < seqstr.length()-idlen; i++) {
			    					if(filtercode.contains(seqstr.substring(i, i+idlen).hashCode()))
			    					{
			    						filterflag=true;
			    						break;
			    					}
			    				}
			    			 if(!filterflag&&ForwardStrand.size()<=maxctrlseq)
			    			 {
				    			 ForwardStrand.add(seqstr);
				    			 TotalLen+=seqstr.length();
				    			 accSeqLen.add(TotalLen);
			    			 }
			    			 else
			    				 filternum++;
			    			 seqstr="";
			    			 }
			    			 
			    		 }
			    		 else
			    		 {
			    			 seqstr+=Line.trim();
			    		 }
			    	 }
			     }

   		          System.out.println("Filter Sequences: "+filternum);
   		    
   	
			     
			
		} catch (IllegalArgumentException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public void build_index(String inputfile) {
		// TODO Auto-generated method stub
		 try {
			 ForwardStrand.clear();
			 accSeqLen.clear();
			// ReverseStrand.clear();
			 //Database to hold the training set
			      BufferedReader br = new BufferedReader(new FileReader(inputfile));
			      String Line="";
			      String seqstr="";
			      TotalLen=0;
	   		       accSeqLen.add(0);
			     while( (Line=br.readLine())!=null)
			     {
			    	 if(Line.length()>0)
			    	 {
			    		 if(Line.charAt(0)=='>')
			    		 {
			    			 if(seqstr!="")
			    			 {
			    			 seqstr=seqstr.replace("N", "");
			    			 ForwardStrand.add(seqstr);
			    			 TotalLen+=seqstr.length();
			    			 accSeqLen.add(TotalLen);
			    			 seqstr="";
			    			 }
			    			 
			    		 }
			    		 else
			    		 {
			    			 seqstr+=Line.trim();
			    		 }
			    	 }
			     }

   		          
   		    
   	
			     
			
		} catch (IllegalArgumentException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

	public LinkedList<FastaLocation> searchPattern(String pattern, int mismatch) {
		// TODO Auto-generated method stub
		if(num_thread>ForwardStrand.size())
			num_thread=1;
	    int workSize=ForwardStrand.size()/num_thread+1;
	    Iterator<String> iter=ForwardStrand.iterator();
	    LinkedList<FastaLocation> search_result=new LinkedList<FastaLocation>();
	    int count=0;
	    ArrayList<SearchThread> threadpool=new ArrayList<SearchThread>(num_thread);
	    //Forward search
	    for (int i = 0; i < num_thread; i++) {
	    	SearchThread t1 = new SearchThread(pattern, mismatch, ForwardStrand.subList(i*workSize,Math.min(ForwardStrand.size(),(i+1)*workSize )),i*workSize,accSeqLen);
			if(background!=null)
			{
				t1.bgmodel=background;
				t1.BGscoreMap=this.BGscoreMap;
			}
	    	t1.start();
			threadpool.add(t1);
		}

	    
	    //join and wait
		try {
	    for (int i = 0; i < threadpool.size(); i++) {
	
				threadpool.get(i).join();
				if(i==num_thread)
					forwardCount=search_result.size();
			search_result.addAll(threadpool.get(i).getResult());
			

		}
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	    
		return search_result;
	}
	public LinkedList<FastaLocation> samplingPattern_PS(PWM pattern, int Num_sample) {
		if(num_thread>ForwardStrand.size())
			num_thread=1;
		// TODO Auto-generated method stub
	    int workSize=ForwardStrand.size()/num_thread+1;
	    Iterator<String> iter=ForwardStrand.iterator();
	    LinkedList<FastaLocation> search_result=new LinkedList<FastaLocation>();
	    int count=0;
	    ArrayList<SamplingThread_PS> threadpool=new ArrayList<SamplingThread_PS>(num_thread);

	    //Forward search
	    for (int i = 0; i < num_thread; i++) {
	    	SamplingThread_PS t1 = new SamplingThread_PS(pattern, Num_sample/num_thread, ForwardStrand.subList(i*workSize,Math.min(ForwardStrand.size(),(i+1)*workSize )),i*workSize,accSeqLen);

	    	t1.start();
			threadpool.add(t1);
		}

		try {
	    for (int i = 0; i < threadpool.size(); i++) {
				threadpool.get(i).join();
				if(i==num_thread)
					forwardCount=search_result.size();
			search_result.addAll(threadpool.get(i).getResult());
		}
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	    
		return search_result;
	}
	public LinkedList<FastaLocation> samplingPattern(PWM pattern, int Num_sample) {
		if(num_thread>ForwardStrand.size())
			num_thread=1;
		// TODO Auto-generated method stub
	    int workSize=ForwardStrand.size()/num_thread+1;
	    Iterator<String> iter=ForwardStrand.iterator();
	    LinkedList<FastaLocation> search_result=new LinkedList<FastaLocation>();
	    int count=0;
	    ArrayList<SamplingThread> threadpool=new ArrayList<SamplingThread>(num_thread);

	    //Forward search
	    for (int i = 0; i < num_thread; i++) {
	    	SamplingThread t1 = new SamplingThread(pattern, Num_sample/num_thread, ForwardStrand.subList(i*workSize,Math.min(ForwardStrand.size(),(i+1)*workSize )),i*workSize,accSeqLen);

	    	t1.start();
			threadpool.add(t1);
		}

		try {
	    for (int i = 0; i < threadpool.size(); i++) {
				threadpool.get(i).join();
				if(i==num_thread)
					forwardCount=search_result.size();
			search_result.addAll(threadpool.get(i).getResult());
		}
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	    
		return search_result;
	}
	
	public LinkedList<FastaLocation> searchPattern(PWM pattern, double thresh) {
		// TODO Auto-generated method stub
		if(num_thread>ForwardStrand.size())
			num_thread=1;
	    int workSize=ForwardStrand.size()/num_thread+1;
	    Iterator<String> iter=ForwardStrand.iterator();
	    LinkedList<FastaLocation> search_result=new LinkedList<FastaLocation>();
	    int count=0;
	    ArrayList<SearchThread> threadpool=new ArrayList<SearchThread>(num_thread);
	    if(background!=null)
		{
			if(BGscoreMap!=null)
			{
				if(BGscoreMap.containsKey(pattern.core_motiflen))
				{
					for (int i = 0; i < ForwardStrand.size(); i++) {
						if(!BGscoreMap.get(pattern.core_motiflen).containsKey(i))
						{
							BGscoreMap.get(pattern.core_motiflen).put(i, new  ArrayList<Double>() );
						}
					}
				}
				else
				{
					BGscoreMap.put(pattern.core_motiflen,new HashMap<Integer, ArrayList<Double>>());
					for (int i = 0; i < ForwardStrand.size(); i++) {
						BGscoreMap.get(pattern.core_motiflen).put(i, new  ArrayList<Double>() );
					}
				}
			}
			
		}
	    //Forward search
	    for (int i = 0; i < num_thread; i++) {
	    	SearchThread t1 = new SearchThread(pattern, thresh, ForwardStrand.subList(i*workSize,Math.min(ForwardStrand.size(),(i+1)*workSize )),i*workSize,accSeqLen);
	    	if(background!=null)
	    	{
				t1.bgmodel=background;
				t1.BGscoreMap=this.BGscoreMap;
				
	    	}
	    	t1.start();
			threadpool.add(t1);
		}
//	    //Reverse search
//	    for (int i = 0; i < num_thread; i++) {
//	    	SearchThread t1 = new SearchThread(pattern, thresh, ReverseStrand.subList(i*workSize,Math.min(ReverseStrand.size(),(i+1)*workSize )));
//			t1.start();
//			threadpool.add(t1);
//		}
	    
	    //join and wait
		try {
	    for (int i = 0; i < threadpool.size(); i++) {
				threadpool.get(i).join();
				if(i==num_thread)
					forwardCount=search_result.size();
			search_result.addAll(threadpool.get(i).getResult());
			if(SearchThread.recordSiteThreshold<Double.POSITIVE_INFINITY)
				pattern.matchsite.addAll(threadpool.get(i).matchsite);
		    
		}
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	    
		return search_result;
	}



	public String getSite(int seqid, int location, int len) {
		// TODO Auto-generated method stub
		String temp=ForwardStrand.get(seqid);
		 String onsite=temp.substring(Math.max(0, location),  Math.min(temp.length(), location+len));
		 int head=-location;
		 int tail=location+len-temp.length();
		 if(head>0)
			 for (int i = 0; i < head; i++) {
				onsite="N"+onsite;
			}
		 if(tail>0)
			 for (int i = 0; i < tail; i++) {
					onsite=onsite+"N";
				}
		 return onsite;
	}

	public int getTotalLength() {
		// TODO Auto-generated method stub
		return TotalLen;
	}

	public int getSeqNum() {
		// TODO Auto-generated method stub
		return ForwardStrand.size();
	}
	
	
	public double[][] seq_similarity()
	{
		int test=common.longestSubstr("AAACAATTT", "AAAGTTATT");
		double[][] sim_arr=new double[getSeqNum()][getSeqNum()];
		Iterator<String> iter=ForwardStrand.iterator();
		int seqid=0;
		 TreeMap<String, Integer> StringMap=new TreeMap<String, Integer>();
		 int wordlen=18;
		 int seqnum=getSeqNum();
		 while(iter.hasNext())
		 {
			 String seq=iter.next();
			 for (int i = 0; i < seq.length()-wordlen; i++) {
			       String word=seq.substring(i,i+wordlen);
			       if(StringMap.containsKey(word))
			       {
			    	   if(StringMap.get(word)!=seqid)
			    		   StringMap.put(word,StringMap.get(word)*seqnum+seqid);
			    	   i+=wordlen;
			       }
			       else
			       {
			    	   StringMap.put(word, seqid);
			       }
			}
			 seqid++;
		 }
		 //decode process
		 for(Integer seqgroup : StringMap.values())
		 {
			 if(seqgroup<getSeqNum())
				 continue;
			 int round=(int)Math.ceil(Math.log(seqgroup)/Math.log(seqnum));
			 ArrayList<Integer> seqlist=new ArrayList<Integer>(round);
			 for (int i = 0; i < round; i++) {
				 int rest=seqgroup%seqnum;
				 seqlist.add(rest);
				 seqgroup=(seqgroup-rest)/seqnum;
			}
			 System.out.print(round+":");
			 for (int i = 0; i <seqlist.size()-1; i++) {
				 int ii=seqlist.get(i);
				 String seq1=ForwardStrand.get(ii);
				 for (int j = i+1; j < seqlist.size(); j++)
				 {
					 int jj=seqlist.get(j);
					 if(sim_arr[ii][jj]==0)
					 {
						 String seq2=ForwardStrand.get(jj);
						 int lcs=common.longestSubstr(seq1, seq2);
						 int lcs2=common.longestSubstr(common.getReverseCompletementString(seq1), seq2);
						 lcs=Math.max(lcs, lcs2);
						 System.out.println(lcs);
						 sim_arr[ii][jj]=lcs;
						 sim_arr[jj][ii]=lcs;
					 }
				 }
		 }
		 }
		
		return sim_arr;
	}
	
	public PWM getPWM_fromMatch(PWM pattern)
	{
		PWM ret=null;
		LinkedList<FastaLocation> origFalocs=this.searchPattern(pattern, pattern.core_motiflen*Math.log(0.25));
		Iterator<FastaLocation> iter=origFalocs.iterator();
		String[] allsites=new String[origFalocs.size()];
		int sitecount=0;
		while(iter.hasNext())
		{
			FastaLocation currloc=iter.next();

			String site=this.getSite(currloc.getSeqId(),currloc.getSeqPos()-pattern.head, pattern.columns());
			if(currloc.ReverseStrand)
				site=common.getReverseCompletementString(site);
			allsites[sitecount]=site;
			sitecount++;
		}
		
			try {
				ret=new PWM(allsites);
			} catch (IllegalAlphabetException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (IllegalSymbolException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}

		return ret;
	}

}
