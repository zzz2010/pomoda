import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

import org.biojava.utils.ChangeVetoException;


public class SearchThread extends Thread  {

	static boolean bestonly=false;  //only get the best occurrence for each sequence
	static double recordSiteThreshold=Double.POSITIVE_INFINITY;
	private LinkedList<FastaLocation> result;
	public LinkedList<Integer> matchsite=null;
	PWM motif;
	double thresh;
	String pattern;
	int startSeqId;
	int mismatch;
	int dbsize;
	public HashMap<Integer,HashMap<Integer,ArrayList<Double>>> BGscoreMap;
	BGModel bgmodel=null;
	public List<String> db; //in case want to run multi db, before get result
	boolean PWMflag=false;
	public ArrayList<Integer> accSeqLen;
	public SearchThread(PWM motif,double thresh, List<String> db,int startseqNum,ArrayList<Integer> SeqLenArr)
	{
		startSeqId= startseqNum;
		this.motif=motif;
		this.thresh=thresh;
		this.db=db;
		this.result=new LinkedList<FastaLocation>();
		this.PWMflag=true;
		accSeqLen=SeqLenArr;
	}
	
	public SearchThread(String pattern,int mismatch, List<String> db,int startseqNum,ArrayList<Integer> SeqLenArr)
	{
		this.pattern=pattern;
		this.mismatch=mismatch;
		this.db=db;
		result=new LinkedList<FastaLocation>();
		this.PWMflag=false;
		accSeqLen=SeqLenArr;
		
	}
	@Override
	public void run() {
		// TODO Auto-generated method stub
		if(PWMflag)
		{
			int pos=accSeqLen.get(startSeqId);
			boolean bg_buff_ready=true;
			if(Double.POSITIVE_INFINITY>recordSiteThreshold)
			      matchsite=new LinkedList<Integer>();
			Iterator<String> iter=db.iterator();
			int seqid=startSeqId;
			ArrayList<FastaLocation> tempbestpos=new ArrayList<FastaLocation>(1000);
			while(iter.hasNext())
			{

				String seq=iter.next();
				double bestscore=Double.NEGATIVE_INFINITY;
				if(bestonly)
					tempbestpos.clear();
				FastaLocation bestpos=null;
				try {

					for (int i = 0; i < seq.length()-motif.core_motiflen; i++) {
						String temp=seq.substring(i,i+motif.core_motiflen);
						double score=motif.scoreWeightMatrix(temp);
						double bgscore=0;
						if(bgmodel!=null)
						{
							if(BGscoreMap.get(motif.core_motiflen).get(seqid).size()>i)
							{
								bgscore=BGscoreMap.get(motif.core_motiflen).get(seqid).get(i);
							}
							else
							{
								bgscore=bgmodel.Get_LOGPROB(temp);
								BGscoreMap.get(motif.core_motiflen).get(seqid).add(bgscore);
							}
						}
						
						boolean reverse=false;
						double score2=motif.scoreWeightMatrix(common.getReverseCompletementString(temp));
						
						if(score2>score)
						{
							score=score2;
							reverse=true;
						}
						
						score-=bgscore;

						if(score>recordSiteThreshold)
						{
							matchsite.add(pos+i);
						}
						
						if(score>thresh)
						{
							int addpos=pos+i;
							
						FastaLocation fapos=new FastaLocation(addpos,seqid , i, seq.length());
						//fapos.seq=temp;
	
							fapos.ReverseStrand=reverse;
					      fapos.Score=score;
					      if(bestonly)
					      {
					    	  if(bestscore<score)
					    	  {
					    		  bestpos=fapos;
					    		  bestscore=score;
					    		  tempbestpos.clear();
					    	  }
					    	  else if(bestscore==score)
					    	  {
					    		  tempbestpos.add(fapos);
					    	  }
					      }
					      else
					    	  result.add(fapos);
						}
						
					}
					if(bestonly&&bestpos!=null)
					{
						result.add(bestpos);
						result.addAll(tempbestpos);
					}					 
				    pos+=seq.length();
				    seqid++;
				} catch (ChangeVetoException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		}
		else
		{
			int pos=accSeqLen.get(startSeqId);
			
			Iterator<String> iter=db.iterator();
			int seqid=startSeqId;
while(iter.hasNext())
{

			String seq=iter.next().toUpperCase();
			
			for (int i = 0; i < seq.length()-pattern.length()+1; i++) {
				int num_mismatch=0;
				boolean reverse=false;
				for (int j = 0; j < pattern.length(); j++) {
					if(pattern.charAt(j)=='N'|| pattern.charAt(j)==seq.charAt(i+j))
					{
						continue;
					}
					else
					{
						num_mismatch++;
						if(num_mismatch>mismatch)
							break;
					}
				}
				
				
				//reverse strand
				String rcseq=common.getReverseCompletementString(seq.substring(i,i+pattern.length()));
				int num_mismatch2=0;
				for (int j = 0; j < pattern.length(); j++) {
					if(pattern.charAt(j)=='N'|| pattern.charAt(j)==rcseq.charAt(j))
					{
						continue;
					}
					else
					{
						num_mismatch2++;
						if(num_mismatch2>mismatch)
							break;
					}
				}
				if(num_mismatch2<num_mismatch)
				{
					reverse=true;
					num_mismatch=num_mismatch2;
				}
				
				
				if(num_mismatch<=mismatch)
				{
					FastaLocation fapos=new FastaLocation(pos+i,seqid , i, seq.length());
				    fapos.Score=num_mismatch;  
				    fapos.ReverseStrand=reverse;
					result.add(fapos);
				}
					

			}
		    pos+=seq.length();
		    seqid++;
			
}
			
		}
			

	}
	
	public LinkedList<FastaLocation> getResult()
	{
		return result;
	}
	

}
