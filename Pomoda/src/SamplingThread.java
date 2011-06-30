	import java.util.ArrayList;
	import java.util.HashMap;
	import java.util.Iterator;
	import java.util.LinkedList;
	import java.util.List;
import java.util.Random;

import org.biojava.utils.ChangeVetoException;

	public class SamplingThread  extends Thread {



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
		public SamplingThread(PWM motif,int samplenum, List<String> db,int startseqNum,ArrayList<Integer> SeqLenArr)
		{
			startSeqId= startseqNum;
			this.motif=motif;
			this.samplenum=samplenum;
			this.db=db;
			this.result=new LinkedList<FastaLocation>();
			this.PWMflag=true;
			accSeqLen=SeqLenArr;
		}
		
		public SamplingThread(String pattern,int mismatch, List<String> db,int startseqNum,ArrayList<Integer> SeqLenArr)
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
				int totallen_db=accSeqLen.get(startSeqId+db.size())-accSeqLen.get(startSeqId);
				double enhancefactor=motif.core_motiflen*Math.log(4)+Math.log((double)samplenum/totallen_db);
				//enhancefactor=Double.MAX_VALUE; //sample all site
				double nulllog=motif.core_motiflen*Math.log(0.25);
				Iterator<String> iter=db.iterator();
				int seqid=startSeqId;
				while(iter.hasNext())
				{

					String seq=iter.next();
					double bestscore=Double.NEGATIVE_INFINITY;

					FastaLocation bestpos=null;
					try {
						double minprob=Math.log(1.0/6);
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
//							double Prior_EZ=0.0012468827930174563;
//							double loglik=score-nulllog-Math.log(2.0)+Math.log(Prior_EZ/(1-Prior_EZ));
//							
//							double prob_theta=Math.exp(loglik)/(Math.exp(loglik)+1);
//							score=Math.log(prob_theta*20);
							
							
							score=enhancefactor+score;
						
//							if(minprob>score)
//								continue;
							if(score>0)
								score=0;

							if(score>Math.log(rand.nextDouble()))
							{
								int addpos=pos+i;
								
							FastaLocation fapos=new FastaLocation(addpos,seqid , i, seq.length());
							//fapos.seq=temp;
		
							  fapos.ReverseStrand=reverse;
						      fapos.Score=score;//Math.log(Math.exp(score)+minprob);

						    	  result.add(fapos);
							}
							
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
