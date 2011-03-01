import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.TreeMap;

import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;


public class GapModelingThread extends Thread {

	List<String> Sites;
	int gapStart;
	int gapEnd;
	HashSet<Integer> depend_Pos;
	PWM gapPWM;
	double KL_improve;
	HashMap<String,Double> DprobMap;
	
	public GapModelingThread(int gapstart, int gapend, List<String> sites,HashSet<Integer> depend_pos)
	{
		this.gapEnd=gapend;
		this.gapStart=gapstart;
		this.Sites=sites;
		this.depend_Pos=depend_pos;
		KL_improve=common.DoubleMinNormal;
		DprobMap=new HashMap<String,Double>();
	}
	
	public void run() {
		String[] gapstr=new String[ Sites.size()];
		Iterator<String> iter=Sites.iterator();
		int i=0;
		while(iter.hasNext())
		{
			gapstr[i]=iter.next().substring(gapStart, gapEnd);
			i++;
		}
		try {
			gapPWM=new PWM(gapstr);
			//build Kmer model
			int dmerSize=1<<(depend_Pos.size() *2);
			double[] dmerCount=new double[dmerSize];
			for (int j = 0; j < gapstr.length; j++) {
				String dmer="";
				Iterator<Integer> iter2=depend_Pos.iterator();
				while(iter2.hasNext())
				{
					int dpos=iter2.next()-gapStart;
					dmer+=gapstr[j].charAt(dpos);
				}
				int hash=common.getHashing(dmer, 0, depend_Pos.size());
				dmerCount[hash]+=1;
			}
			common.Normalize(dmerCount);
			TreeMap<Double,Integer> sorted_column=new TreeMap<Double,Integer>();
			for (int j = 0; j< dmerSize; j++)
			{  
				double weight=dmerCount[j];
				sorted_column.put(weight, j-1);
			}
			double sumprob=0;
			for(Double key:sorted_column.descendingKeySet())
			{
				DprobMap.put(common.Hash2ACGT(sorted_column.get(key), depend_Pos.size()), key);
				sumprob+=key;
				if(DprobMap.size()==(4* depend_Pos.size()-1))
					break;
			}
			DprobMap.put("N", (1-sumprob)/(dmerSize-4*depend_Pos.size()+1));
			
			//compute KL-divergence
			
			
			
		} catch (IllegalAlphabetException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IllegalSymbolException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
}
