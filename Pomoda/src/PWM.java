/**
 * @author zhizhuo zhang
 * zzz2010@gmail.com
 */

import java.util.*;
import java.util.Map.Entry;

import java.util.Iterator;

import org.biojava.bio.alignment.*;
import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dist.DistributionTools;
import org.biojava.bio.dp.DP;
import org.biojava.bio.dp.ScoreType;
import org.biojava.bio.dp.SimpleWeightMatrix;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.*;

public class PWM extends SimpleWeightMatrix {

	double[][]  m_matrix;
	public int head;
	public int tail;
	public double Score;
	public ArrayList<Double>pos_prior;
	public PWM(Distribution[] arg0) throws IllegalAlphabetException {
		super(arg0);
		// TODO Auto-generated constructor stub
		
		head=-1;
		tail=-1;
		pos_prior=new ArrayList<Double>();
		Score=Double.MIN_VALUE;
	}
	
	
	
	public static final Distribution[] alignment2Distribution(String[] alignments) throws IllegalSymbolException, IllegalAlphabetException
	{
		Map<String, SymbolList> map = new HashMap<String, SymbolList>();
		for (int i = 0; i < alignments.length; i++) {
			map.put(String.valueOf(i),DNATools.createDNA( alignments[i])) ;
		}
		Alignment align = new SimpleAlignment(map);
	    //make a Distribution[] of the motif
	    Distribution[] dists =
	        DistributionTools.distOverAlignment(align, false, common.DoubleMinNormal);
	    
	    
	    
	    return	dists;
			
	}
	
	public void setWeights(int col, double[] weights)
	{
		for (int i = 0; i < weights.length; i++) {
			m_matrix[col][i]=weights[i];
		}
		Distribution di=this.getColumn(col);
		SymbolList sla;
		double weight=0;
		
		try {
			sla =DNATools.createDNA("ACGT");
			for (int i = 1; i <= weights.length; i++) {
				this.getColumn(col).setWeight(sla.symbolAt(i), weights[i-1]);
			}
		} catch (IllegalSymbolException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IndexOutOfBoundsException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
	
	public void setWeight(int col, int symid,double weight)
	{
		m_matrix[col][symid]=weight;
		Distribution di=this.getColumn(col);
		SymbolList sla;

			try {
				sla =DNATools.createDNA("ACGT");
				this.getColumn(col).setWeight(sla.symbolAt(symid), weight);
			} catch (IllegalSymbolException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			
		
		
	}
	
	
	public double getWeight(int col, int symid)
	{

		return m_matrix[col][symid];
		
	}

	public PWM(String[] alignments) throws IllegalAlphabetException, IllegalSymbolException  {

	    this(alignment2Distribution(alignments));
		
	    
		SymbolList sla;
		try {
			sla = DNATools.createDNA("ACGT");
			m_matrix=new double[this.columns()][4];
			for (int i = 0; i < this.columns(); i++) {
				Distribution di=this.getColumn(i);
				for (int j = 1; j<= 4; j++)
				{  
					double weight=di.getWeight(sla.symbolAt(j));
			        m_matrix[i][j-1]=weight;
						
				}
				
		
				
			}
		} catch (IllegalSymbolException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		

	    
	    
	}
	
	
	public double scoreWeightMatrix( String seq, ScoreType scoreType)
	{
		try {
			return DP.scoreWeightMatrix(this, DNATools.createDNA(seq), scoreType, 1);
		} catch (IllegalSymbolException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return Double.MIN_VALUE;
		}
		
	}
	
	
	private double getBackgroundLogProb(String pattern,BGModel model)
	{
		double logprob=0;
		char[] ACGT=new char[]{'A','C','G','T'};
		StringBuffer sb=new StringBuffer(pattern);
		//fill N with Random
		Random rand=new Random();
		int Ncount=0;
		for (int i = 0; i < pattern.length(); i++) {
			if(pattern.charAt(i)=='N')
			{
				sb.setCharAt(i,  ACGT[rand.nextInt(4)]);
				Ncount++;
			}
			
		}
		logprob=model.Get_LOGPROB(sb.toString())+Ncount*Math.log(4);
		
		
		return logprob;
	}
	
	
	//sort symid based on the weight: desc order.
	int[] columnOrder(int col)
	{
		int [] order=new int[4];
		TreeMap<Double,Integer> column=new TreeMap<Double,Integer>();
		

			double smallvalue=common.DoubleMinNormal;
			for (int j = 1; j<= 4; j++)
			{  
				double weight=m_matrix[col][j-1];
				column.put(smallvalue*j-weight, j-1);
			}
			Iterator<Integer> iter= column.values().iterator();
			int c=0;
			 while(iter.hasNext())
				 {order[c]=iter.next();
				 c++;
				        
				 }
				 
			
		

		
		return order;
	}
	
	
	// Exact Priority Queue method generate instance from this motif model
	public LinkedList<Map.Entry<Double,String>> GenerateInstanceFromPWMPQ(double sampleratio,double FDRthresh,BGModel bgmodel)
	{
		LinkedList<Map.Entry<Double,String>> InstanceSet=new LinkedList<Map.Entry<Double,String>>();
		TreeMap<Double,String> PQ=new TreeMap<Double,String>();
		double sumFDR=0;
		double sumProb=0;
		double smallscale=0.0000000000001;
		char[] ACGT=new char[]{'A','C','G','T'};
		String consensus=this.Censensus(false);
        ArrayList<Integer>  effIndex=new ArrayList<Integer>(consensus.length());
        ArrayList< int[]> orderedCol=new ArrayList< int[]> (consensus.length());
      
        for (int i = 0; i < consensus.length(); i++) {
			if(consensus.charAt(i)!='N')
			{
				effIndex.add(i);
				orderedCol.add(columnOrder(i));
			}
        }
        int effLen=effIndex.size();
    	if(effLen==0)
    		return InstanceSet;
    	int head=effIndex.get(0);
    	int tail=consensus.length()-effIndex.get(effLen-1)-1;
    	//get first one into PQ
    	StringBuffer sb=new StringBuffer("");
    	double logprob=0;
    	double log025=Math.log(0.25);
    	for (int i = 0; i < effLen; i++) {
			int symid=orderedCol.get(i)[0];
			sb.append(ACGT[symid]);
			logprob+=Math.log(getWeight(effIndex.get(i),symid));
			if(i!=effLen-1)
			{
				int gapsize=effIndex.get(i+1)-effIndex.get(i)-1;
				for (int j = 0; j < gapsize; j++) {
					sb.append('N');	
					//logprob+=log025;
				}
			}
		}
    	PQ.put(logprob, sb.toString());
    	Entry<Double,String> topEntry=PQ.pollLastEntry();
    	sumFDR+=Math.exp(getBackgroundLogProb(topEntry.getValue(),bgmodel));
    	sumProb+=Math.exp(logprob);
    	
		while(sumFDR<FDRthresh )
		{
			InstanceSet.push(topEntry);
			if(sumProb>sampleratio)
				break;
			//mutate one mismatch
			for (int z = 0; z < effLen; z++) {
				logprob=topEntry.getKey();
				int changingIndex=effIndex.get(z);
				sb=new StringBuffer(topEntry.getValue());
				
					int orig_symid=common.acgt( sb.charAt(changingIndex-head));
					int symid=-1;
					int[] orders=orderedCol.get(z);
					int i=0;
					for ( i= 0; i < orders.length; i++) {
					     if(orders[i]==orig_symid)
					    	 break;
					}
					if(i==(orders.length-1))
						continue;
					symid=orders[i+1];
					if(getWeight(changingIndex, symid)<0.05)
						continue;
						
											
					sb.setCharAt(changingIndex-head,(ACGT[symid]));
					logprob+=Math.log(getWeight(changingIndex,symid));
					logprob-=Math.log(getWeight(changingIndex,orig_symid));
					
				PQ.put(logprob, sb.toString());	 
				
			}
			// get pop out the best one from PQ
			topEntry=PQ.pollLastEntry();
			if(topEntry==null)
				break;
	    	sumFDR+=Math.exp(getBackgroundLogProb(topEntry.getValue(),bgmodel));
	    	sumProb+=Math.exp(logprob);						
		}
		
		
		
	   return InstanceSet;	
	}
	
	
	//print matrix
	public void print()
	{
		
		for (int i = 0; i < this.columns(); i++) {
			StringBuffer sb=new StringBuffer("");
			for (int j = 0; j< 4; j++)
			{  
				double weight=m_matrix[i][j];
					sb.append(String.valueOf(weight));
					sb.append('\t');
					
			}
			System.out.println(sb.toString());

			
		}
		
	}
	
	public String Censensus(boolean trim)
	{
		String consensus="";
		String  ACGT="ACGT";
		for (int i = 0; i < this.columns(); i++) {
			StringBuffer sb=new StringBuffer("");
			for (int j = 1; j<= 4; j++)
			{  
				double weight=m_matrix[i][j-1];
				if(weight>0.3)
					sb.append(ACGT.charAt(j-1));
					
			}
			String consensus_pattern=sb.toString();
			if(consensus_pattern.equalsIgnoreCase( "" )) {consensus = consensus + "N"; }
			else if (consensus_pattern.equalsIgnoreCase( "A") ) { consensus = consensus + "A"; }
			else if (consensus_pattern.equalsIgnoreCase("C")  )	{ consensus = consensus + "C"; }
			else if (consensus_pattern.equalsIgnoreCase( "G" ))	{ consensus = consensus + "G"; }
			else if (consensus_pattern.equalsIgnoreCase( "T" ))	{ consensus = consensus + "T"; }
			else if (consensus_pattern.equalsIgnoreCase( "AC" ))	{ consensus = consensus + "M"; }
			else if (consensus_pattern.equalsIgnoreCase( "AG" ))	{ consensus = consensus + "R"; }
			else if (consensus_pattern.equalsIgnoreCase( "AT" ))	{ consensus = consensus + "W"; }
			else if (consensus_pattern.equalsIgnoreCase( "CG" ))	{ consensus = consensus + "S"; }
			else if (consensus_pattern.equalsIgnoreCase( "CT" ))	{ consensus = consensus + "Y"; }
			else if (consensus_pattern.equalsIgnoreCase( "GT" ))	{ consensus = consensus + "K"; }
			else if (consensus_pattern.equalsIgnoreCase( "ACG" ))	{ consensus = consensus + "V"; }
			else if (consensus_pattern.equalsIgnoreCase( "ACT" ))	{ consensus = consensus + "H"; }
			else if (consensus_pattern.equalsIgnoreCase( "AGT" ))	{ consensus = consensus + "D"; }
			else if (consensus_pattern.equalsIgnoreCase( "CGT" ))	{ consensus = consensus + "B"; }
			
		}

		int start,end;
		start=end=-1;
		for (int i = 0; i < consensus.length(); i++) {
			if(consensus.charAt(i)!='N'& start==-1)
				start=i;
			if(consensus.charAt(consensus.length()-i-1)!='N'& end==-1)
				end=consensus.length()-i;
		}
		head=start;
		tail=consensus.length()-end;
		if(trim)
		{
		consensus=consensus.substring(start,end);
		}
		return consensus;
	}
	
	
	

}
