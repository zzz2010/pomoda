import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;

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
	public HashMap<Integer,HashMap<Integer,ArrayList<Double>>> BGscoreMap=null;
    public int forwardCount; //the first n entry in the searchPattern is forward.
	public int TotalLen=0;
	public ArrayList<Integer> accSeqLen;
	public LinearEngine(int num_thread) {
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
	
	
	public void build_index(String inputfile) {
		// TODO Auto-generated method stub
		 try {
			 ForwardStrand.clear();
			 accSeqLen.clear();
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

	public LinkedList<FastaLocation> searchPattern(String pattern, int mismatch) {
		// TODO Auto-generated method stub
	    int workSize=ForwardStrand.size()/num_thread+1;
	    Iterator<String> iter=ForwardStrand.iterator();
	    LinkedList<FastaLocation> search_result=new LinkedList<FastaLocation>();
	    int count=0;
	    ArrayList<SearchThread> threadpool=new ArrayList<SearchThread>(num_thread);
	    //Forward search
	    for (int i = 0; i < num_thread; i++) {
	    	SearchThread t1 = new SearchThread(pattern, mismatch, ForwardStrand.subList(i*workSize,Math.min(ForwardStrand.size(),(i+1)*workSize ) ),i*workSize,accSeqLen);
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
	public LinkedList<FastaLocation> samplingPattern(PWM pattern, int Num_sample) {
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

}
