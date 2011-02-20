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
	int mismatch;
	int dbsize;
	public List<Sequence> db; //in case want to run multi db, before get result
	boolean PWMflag=false;
	public SearchThread(PWM motif,double thresh, List<Sequence> db)
	{
		this.motif=motif;
		this.thresh=thresh;
		this.db=db;
		this.result=new LinkedList<FastaLocation>();
		this.PWMflag=true;
		
	}
	
	public SearchThread(String pattern,int mismatch, List<Sequence> db)
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
		
			Iterator<Sequence> iter=db.iterator();
			while(iter.hasNext())
			{

				Sequence seq=iter.next();
				try {
					String seq2 = seq.seqString();
					for (int i = 0; i < seq2.length()-motif.core_motiflen; i++) {
						String temp=seq2.substring(i,i+motif.core_motiflen);
						double score=motif.scoreWeightMatrix(temp);
						double score2=motif.scoreWeightMatrix(common.getReverseCompletementString(temp));
						if(score2>score)
							score=score2;
						if(score>thresh)
						{
						FastaLocation fapos=new FastaLocation(pos+i,Integer.parseInt(seq.getName()) , i, seq.length());
					      fapos.Score=score;
					      result.add(fapos);
						}
						
					}
					 
				    pos+=seq.length();
				} catch (ChangeVetoException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		}
		else
		{
			int pos=0;
			SymbolList sym_pattern;
			try {
				sym_pattern = DNATools.createDNA(pattern);
				Iterator<Sequence> iter=db.iterator();
	
			while(iter.hasNext())
			{

				Sequence seq=iter.next();
				
				for (int i = 0; i < seq.length()-pattern.length()+1; i++) {
					int num_mismatch=0;
					for (int j = 0; j < pattern.length(); j++) {
						if(sym_pattern.symbolAt(j).getMatches().contains(seq.symbolAt(i+j)))
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
						FastaLocation fapos=new FastaLocation(pos+i,Integer.parseInt(seq.getName()) , i, seq.length());
					    fapos.Score=num_mismatch;  
						result.add(fapos);
					}
						
					
				}
				
			}
			} catch (IllegalSymbolException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			
		}
			

	}
	
	public LinkedList<FastaLocation> getResult()
	{
		return result;
	}
	

}
