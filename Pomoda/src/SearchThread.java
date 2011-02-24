import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

import org.biojava.bio.BioException;
import org.biojava.bio.dp.WeightMatrixAnnotator;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.Feature;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.utils.ChangeVetoException;


public class SearchThread extends Thread  {

	private LinkedList<FastaLocation> result;
	PWM motif;
	double thresh;
	String pattern;
	int startSeqId;
	int mismatch;
	int dbsize;
	public List<String> db; //in case want to run multi db, before get result
	boolean PWMflag=false;
	public SearchThread(PWM motif,double thresh, List<String> db,int startseqNum)
	{
		startSeqId= startseqNum;
		this.motif=motif;
		this.thresh=thresh;
		this.db=db;
		this.result=new LinkedList<FastaLocation>();
		this.PWMflag=true;
		
	}
	
	public SearchThread(String pattern,int mismatch, List<String> db)
	{
		this.pattern=pattern;
		this.mismatch=mismatch;
		this.db=db;
		result=new LinkedList<FastaLocation>();
		this.PWMflag=false;
		
	}
	@Override
	public void run() {
		// TODO Auto-generated method stub
		if(PWMflag)
		{
			int pos=0;
		
			Iterator<String> iter=db.iterator();
			int seqid=startSeqId;
			
			while(iter.hasNext())
			{

				String seq=iter.next();
				try {

					for (int i = 0; i < seq.length()-motif.core_motiflen; i++) {
						String temp=seq.substring(i,i+motif.core_motiflen);
						double score=motif.scoreWeightMatrix(temp);
						boolean reverse=false;
						double score2=motif.scoreWeightMatrix(common.getReverseCompletementString(temp));
						if(score2>score)
						{
							score=score2;
							reverse=true;
						}
						if(score>thresh)
						{
							int addpos=pos+i;
							
						FastaLocation fapos=new FastaLocation(addpos,seqid , i, seq.length());
						fapos.seq=temp;
						if(reverse)
							fapos.ReverseStrand=true;
					      fapos.Score=score;
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
			int pos=0;
			
			Iterator<String> iter=db.iterator();
			int seqid=0;
while(iter.hasNext())
{

			String seq=iter.next();
			
			for (int i = 0; i < seq.length()-pattern.length()+1; i++) {
				int num_mismatch=0;
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
				if(num_mismatch<=mismatch)
				{
					FastaLocation fapos=new FastaLocation(pos+i,seqid , i, seq.length());
				    fapos.Score=num_mismatch;  
					result.add(fapos);
				}
					
			    pos+=seq.length();
			    seqid++;
			}
			
}
			
		}
			

	}
	
	public LinkedList<FastaLocation> getResult()
	{
		return result;
	}
	

}
