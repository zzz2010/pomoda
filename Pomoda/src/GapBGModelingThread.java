import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.TreeMap;

import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;


public class GapBGModelingThread extends Thread {

	List<String> Sites;
	public int gapStart;
	public int gapEnd;
	public double lamda=1;
	public HashSet<Integer> depend_Pos;
	public PWM gapPWM;
	public double KL_Divergence;
	public HashMap<String,Double> DprobMap;
	public BGModel background;
	
	public GapBGModelingThread(int gapstart, int gapend, List<String> sites,HashSet<Integer> depend_pos,BGModel bg)
	{
		this.gapEnd=gapend;
		this.gapStart=gapstart;
		this.Sites=sites;
		this.depend_Pos=depend_pos;
		KL_Divergence=common.DoubleMinNormal;
		DprobMap=new HashMap<String,Double>();
		background=bg;
	}
	
	double findPriorLamda(double[] Pob,double[] Pbg,double[] Pt)
	{
		double lamda=0.5;
		
		double sqerror=Double.MAX_VALUE;
		double lamdachange=Double.MAX_VALUE;
		
		//sort the Pob , and get index of sorted array
		TreeMap<Double, List<Integer>> map = new TreeMap<Double, List<Integer>>();
		for(int i = 0; i < Pob.length; i++) {
		    List<Integer> ind = map.get(Pob[i]);
		    if(ind == null){
		        ind = new ArrayList<Integer>();
		        map.put(Pob[i], ind);
		    }
		    ind.add(i);
		}

		// Now flatten the list
		List<Integer> indices = new ArrayList<Integer>();
		for(List<Integer> arr : map.values()) {
		    indices.addAll(arr);
		}
		double A,B;
		A=B=0;
		//minimize sum_(Pbg*lamda-Pob)^2
		for (int i = 0; i < indices.size()/2; i++) {
			A+=Pob[indices.get(i)]*Pbg[indices.get(i)];
			B+=Pbg[indices.get(i)];
		}
		lamda=1-A/B;
		
		for (int i = 0; i < Pt.length; i++) {
			Pt[i]=(Pob[i]-(1-lamda)*Pbg[i])/lamda;
			if(Pt[i]<0)
				Pt[i]=0;
		}
		common.Normalize(Pt);
		return lamda;
	}
	
	double findPriorLamda2(double[] Pob,double[] Pbg,double[] Pt)
	{
		double lamda;
		
		double sqerror=Double.MAX_VALUE;
		double lamdachange=Double.MAX_VALUE;
		double lastlamda=lamda=0.5;
		while(sqerror>common.DoubleMinNormal&&lamdachange>common.DoubleMinNormal)
		{
			for (int i = 0; i < Pt.length; i++) {
				Pt[i]=(Pob[i]-(1-lamda)*Pbg[i])/lamda;
				if(Pt[i]<0)
					Pt[i]=0;
			}
			common.Normalize(Pt);
			sqerror=0;
			double A,B;
			B=A=0;
			for (int i = 0; i < Pt.length; i++) {
				double err=Pob[i]-(lamda*Pt[i]+(1-lamda)*Pbg[i]);
				sqerror+=err*err;
				A+=(Pob[i]-Pbg[i])*(Pt[i]-Pbg[i]);
				B+=(Pt[i]-Pbg[i])*(Pt[i]-Pbg[i]);
			}
			//minimize (B*lamda-A)^2
			lamda=A/B;
			if(lamda>1)
				lamda=1;
			if(lamda<0)
				lamda=0;
			lamdachange=Math.abs(lastlamda-lamda);
			lastlamda=lamda;
			
		}
		
		
		return lamda;
	}
	
	@Override
	public String toString()
	{
		String ret="";
		ret=String.valueOf(gapStart)+"-"+String.valueOf(gapEnd)+":"+Arrays.toString(depend_Pos.toArray(new Integer[1]))+'\t'+KL_Divergence+'\t'+lamda; 
		
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
			double[] Pt = new double[dmerSize];
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
			common.Normalize(gapmerCount);
		
			int num_top=4* depend_Pos.size()-1;
			if(depend_Pos.size()>1)
			{
				common.Normalize(dmerCount);
				
				if(background!=null)
				{
				Iterator<Integer> iiter=depend_Pos.iterator();
				int[] gapN=new int[depend_Pos.size()-1];
				int last=iiter.next();
				int ii=0;
				while(iiter.hasNext())
				{
					int curr=iiter.next();
					gapN[ii]=curr-last-1;
					ii++;
					last=curr;
				}
				double[] Pbg=new double[dmerSize];
				for (int j = 0; j< dmerSize; j++)
				{
					String dmer=(common.Hash2ACGT(j, depend_Pos.size()));
					String gapdmer=String.valueOf(dmer.charAt(0));
					for (int k = 0; k < gapN.length; k++) {
						if(gapN[k]>0)
							gapdmer+="N";
						gapdmer+=String.valueOf(dmer.charAt(k+1));
					}
					Pbg[j]=Math.exp(background.Get_LOGPROB(gapdmer));
				}
				 lamda=findPriorLamda(dmerCount,Pbg,Pt);
				 
				}
				else
				{
					for (int j = 0; j< dmerSize; j++)
						Pt[j]=dmerCount[j];
				}
				TreeMap<Double,Integer> sorted_column=new TreeMap<Double,Integer>();
				for (int j = 0; j< dmerSize; j++)
				{  
					double weight=Pt[j]-(j%num_top)*common.DoubleMinNormal;
					if(weight<0)
						weight=0;
					sorted_column.put(weight, j);
				}
			
			double sumprob=0;
			//sort the dmer by desc prob£¬ and take top 4d-1 as the dependency model 
			for(Double key:sorted_column.descendingKeySet())
			{
				DprobMap.put(common.Hash2ACGT(sorted_column.get(key), depend_Pos.size()), key);
				sumprob+=key;

				if(DprobMap.size()==(num_top))
					break;
			}
		
			DprobMap.put("N", (1-sumprob)/(dmerSize-4*depend_Pos.size()+1));
//				if((1-sumprob)<=(num_top*num_top)*common.DoubleMinNormal)
//				{
//					KL_Divergence=Double.MAX_VALUE; //overfit!
//					return;
//				}
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
					
					int di=common.getHashing(dmer, 0, dmer.length()); 

					p=p*Pt[di]/dmerCount[di];
				}
				if(p==0.0)
					continue;
				for (int k = 0; k < gapstrj.length(); k++) {
					int symid=common.acgt(gapstrj.charAt(k));
					if(depend_Pos.size()>1&&depend_Pos.contains(gapStart+k))
						continue;
					logq+=gapPWM.getLogWeight(k, symid);
					
				}
				sum_plogp_p+=p*(Math.log(p) -logq);
				
			}
			
			KL_Divergence=sum_plogp_p;
			
		} catch (IllegalAlphabetException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IllegalSymbolException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
}
