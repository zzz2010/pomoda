import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Random;

import org.biojava.utils.ChangeVetoException;


public class SamplingThread_PS extends Thread {
	//static boolean bestonly=false;  //only get the best occurrence for each sequence
	private LinkedList<FastaLocation> result;
	PWM motif;
	int samplenum;
	String pattern;
	Random rand=new Random(1);
	int startSeqId;
	int mismatch;
	int dbsize;
	public HashMap<Integer,HashMap<Integer,ArrayList<Double>>> BGscoreMap;
	BGModel bgmodel=null;
	public List<String> db; //in case want to run multi db, before get result
	boolean PWMflag=false;
	public ArrayList<Integer> accSeqLen;
	public SamplingThread_PS(PWM motif,int samplenum, List<String> db,int startseqNum,ArrayList<Integer> SeqLenArr)
	{
		startSeqId= startseqNum;
		this.motif=motif;
		this.samplenum=samplenum;
		this.db=db;
		this.result=new LinkedList<FastaLocation>();
		this.PWMflag=true;
		accSeqLen=SeqLenArr;
	}
	
	public SamplingThread_PS(String pattern,int mismatch, List<String> db,int startseqNum,ArrayList<Integer> SeqLenArr)
	{
		this.pattern=pattern;
		this.mismatch=mismatch;
		this.db=db;
		result=new LinkedList<FastaLocation>();
		this.PWMflag=false;
		accSeqLen=SeqLenArr;
		
	}
	public void run2() {
		// TODO Auto-generated method stub
	
			int pos=accSeqLen.get(startSeqId);
			boolean bg_buff_ready=true;
			int totallen_db=accSeqLen.get(startSeqId+db.size())-accSeqLen.get(startSeqId);
			int numberSample_PS=(int) Math.ceil((double)samplenum/db.size());
			int[] selids=new int[numberSample_PS];
			Iterator<String> iter=db.iterator();
			int seqid=startSeqId;
			while(iter.hasNext())
			{

				String seq=iter.next();
				double bestscore=Double.NEGATIVE_INFINITY;
				double[] accProb=new double[seq.length()-motif.core_motiflen];
				boolean[] posStrand=new boolean[seq.length()-motif.core_motiflen];
				
				try {

					for (int i = 0; i < seq.length()-motif.core_motiflen; i++) {
						String temp=seq.substring(i,i+motif.core_motiflen);
						double score=motif.scoreWeightMatrix(temp);
						
						
						boolean reverse=false;
						double score2=motif.scoreWeightMatrix(common.getReverseCompletementString(temp));
						
						
						if(score2>score||(score2==score&&rand.nextBoolean()))
						{
							score=score2;
							reverse=true;
						}
						
						accProb[i]=Math.exp(score);
						posStrand[i]=reverse;						
					}
					accProb=common.Normalize(accProb);
					accProb=common.prefixsum(accProb);
					
					for (int i = 0; i < numberSample_PS; i++) {
						double pnt=rand.nextDouble();
						int sel=Math.abs(Arrays.binarySearch(accProb, pnt))-1;
						selids[i]=sel;
					}
					Arrays.sort(selids);
					for (int i = 0; i < numberSample_PS; i++) {
						int sel=selids[i];
						int addpos=pos+sel;
						FastaLocation fapos=new FastaLocation(addpos,seqid , sel, seq.length());
						  fapos.ReverseStrand=posStrand[sel];
						  if(i>0)
							  fapos.Score=accProb[i]-accProb[i-1];
						  else
							  fapos.Score=accProb[i];
						  
						  fapos.Score=Math.log(fapos.Score*numberSample_PS);

					    	  result.add(fapos);
					}
				    pos+=seq.length();
				    seqid++;
				} catch (ChangeVetoException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		
			

	}
	@Override
	public void run() {
		// TODO Auto-generated method stub
	
			int pos=accSeqLen.get(startSeqId);
			boolean bg_buff_ready=true;
			int totallen_db=accSeqLen.get(startSeqId+db.size())-accSeqLen.get(startSeqId);
			double threshold=1.0/samplenum;
			int numberSample_PS=(int) Math.ceil((double)samplenum/db.size());
			int[] selids=new int[samplenum];
			Iterator<String> iter=db.iterator();
			int seqid=startSeqId;
			double nulllog=motif.core_motiflen*Math.log(0.25);
			double[] accProb=new double[totallen_db];
			boolean[] posStrand=new boolean[totallen_db];
			while(iter.hasNext())
			{

				String seq=iter.next();
				double bestscore=Double.NEGATIVE_INFINITY;

				
				try {

					for (int i = 0; i < seq.length()-motif.core_motiflen; i++) {
						String temp=seq.substring(i,i+motif.core_motiflen);
						double score=motif.scoreWeightMatrix(temp);
						
						
						boolean reverse=false;
						double score2=motif.scoreWeightMatrix(common.getReverseCompletementString(temp));
						
						
						if(score2>score||(score2==score&&rand.nextBoolean()))
						{
							score=score2;
							reverse=true;
						}
						
						double Prior_EZ=0.0012468827930174563;
						double loglik=score-nulllog-Math.log(2.0)+Math.log(Prior_EZ/(1-Prior_EZ));
						
						double prob_theta=Math.exp(loglik)/(Math.exp(loglik)+1);
						score=Math.log(prob_theta);
						
						
						if(score==Double.NEGATIVE_INFINITY)
							accProb[pos+i-accSeqLen.get(startSeqId)]=0;
						else
							accProb[pos+i-accSeqLen.get(startSeqId)]=Math.exp(score);
//						if(accProb[pos+i-accSeqLen.get(startSeqId)]>common.DoubleMinNormal)
//							accProb[pos+i-accSeqLen.get(startSeqId)]=accProb[pos+i-accSeqLen.get(startSeqId)];
						posStrand[pos+i-accSeqLen.get(startSeqId)]=reverse;						
					}

				    pos+=seq.length();
				    seqid++;
				} catch (ChangeVetoException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		
			accProb=common.Normalize(accProb);
			double maxprob=0;
			for (int rr = 0; rr < accProb.length; rr++) {
				if(accProb[rr]>maxprob)
					maxprob=accProb[rr];
			}
			accProb=common.prefixsum(accProb);
			
			for (int i = 0; i < samplenum; i++) {
				double pnt=rand.nextDouble();
				int sel=Math.abs(Arrays.binarySearch(accProb, pnt))-1;
				selids[i]=sel;
			}
			Arrays.sort(selids);
			for (int i = 0; i < samplenum; i++) {
				int sel=selids[i];
				int addpos=sel+accSeqLen.get(startSeqId);
				seqid=Arrays.binarySearch(accSeqLen.toArray(new Integer[1]),addpos);
				if(seqid<0)
				 seqid=Math.abs(seqid)-2;

				FastaLocation fapos=new FastaLocation(addpos,seqid , addpos-accSeqLen.get(seqid), db.get(seqid).length());

				  fapos.ReverseStrand=posStrand[addpos];
				  if(addpos>0)
					  fapos.Score=accProb[addpos]-accProb[addpos-1];
				  else
					  fapos.Score=accProb[addpos];
				  if(fapos.Score<(maxprob*threshold))
					  continue;
				  fapos.Score=Math.log(fapos.Score*samplenum);
			    	  result.add(fapos);
			}

	}
	
	public LinkedList<FastaLocation> getResult()
	{
		return result;
	}
	



}
