import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.TreeMap;

import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;


public class GapModelingThread extends Thread {

	List<String> Sites;
	public int gapStart;
	public int gapEnd;
	public HashSet<Integer> depend_Pos;
	public PWM gapPWM;
	public double KL_improve;
	public HashMap<String,Double> DprobMap;
	
	public GapModelingThread(int gapstart, int gapend, List<String> sites,HashSet<Integer> depend_pos)
	{
		this.gapEnd=gapend;
		this.gapStart=gapstart;
		this.Sites=sites;
		this.depend_Pos=depend_pos;
		KL_improve=common.DoubleMinNormal;
		DprobMap=new HashMap<String,Double>();
	}
	
	@Override
	public String toString()
	{
		String ret="";
		ret=String.valueOf(gapStart)+"-"+String.valueOf(gapEnd)+":"+Arrays.toString(depend_Pos.toArray(new Integer[1]))+'\t'+KL_improve; 
		
		return ret;
	}
	
	@Override
	public void run() {
		String[] gapstr=new String[ Sites.size()];
		Iterator<String> iter=Sites.iterator();
		int i=0;
		while(iter.hasNext())
		{
			String temp=iter.next();
			gapstr[i]=temp.substring(gapStart, gapEnd);
		
			 i++;
		}
		
		try {
			gapPWM=new PWM(gapstr);
			
			//build Kmer model
			int dmerSize=1<<(depend_Pos.size() *2);
			int gapmerSize=1<<((gapEnd-gapStart)*2);
			double[] gapmerCount=new double[gapmerSize];
			double[] dmerCount=new double[dmerSize];
			
			for (int j = 0; j < gapstr.length; j++) {
				String dmer="";
				Iterator<Integer> iter2=depend_Pos.iterator();
				if(gapstr[j].contains("N"))
					continue;
				int whash=common.getHashing(gapstr[j], 0, gapstr[j].length());
				gapmerCount[whash]+=1;
				while(iter2.hasNext())
				{
					int dpos=iter2.next()-gapStart;
					dmer+=gapstr[j].charAt(dpos);
				}
				if(depend_Pos.size()>1)
				{
				int hash=common.getHashing(dmer, 0, depend_Pos.size());
				dmerCount[hash]+=1;
				}
			}
			gapmerCount=	common.Normalize(gapmerCount);
			int num_top=4* depend_Pos.size()-1;
			if(depend_Pos.size()>1)
			{
				dmerCount=common.Normalize(dmerCount);
				TreeMap<Double,Integer> sorted_column=new TreeMap<Double,Integer>();
				for (int j = 0; j< dmerSize; j++)
				{  
					double weight=dmerCount[j]-(j%num_top)*common.DoubleMinNormal;
					sorted_column.put(weight, j);
				}
			
			double sumprob=0;
			for(Double key:sorted_column.descendingKeySet())
			{
				DprobMap.put(common.Hash2ACGT(sorted_column.get(key), depend_Pos.size()), key);
				sumprob+=key;
				if(DprobMap.size()==(num_top))
					break;
			}
		
			DprobMap.put("N", (1-sumprob)/(dmerSize-4*depend_Pos.size()+1));
			}
			//compute KL-divergence
			
			double sum_plogp_p=0;
			for (int j = 0; j < gapmerSize; j++) {
				double p=gapmerCount[j];
				if(p==0.0)
					continue;
				
				String gapstrj=common.Hash2ACGT(j, gapEnd-gapStart);
				double logq=0;
				if(depend_Pos.size()>1)
				{
					String dmer="";
					Iterator<Integer> iter2=depend_Pos.iterator();
					while(iter2.hasNext())
					{
						int dpos=iter2.next()-gapStart;
						dmer+=gapstrj.charAt(dpos);
					}
					
					if(DprobMap.containsKey(dmer))
						logq+=Math.log(DprobMap.get(dmer));
					else
						logq+=Math.log(DprobMap.get("N"));
				}
				for (int k = 0; k < gapstrj.length(); k++) {
					int symid=common.acgt(gapstrj.charAt(k));
					if(depend_Pos.size()>1&&depend_Pos.contains(gapStart+k))
						continue;
					logq+=gapPWM.getLogWeight(k, symid);
					
				}
				sum_plogp_p+=p*(Math.log(p) -logq);
				
			}
			
			KL_improve=sum_plogp_p;
			
		} catch (IllegalAlphabetException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IllegalSymbolException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
}
