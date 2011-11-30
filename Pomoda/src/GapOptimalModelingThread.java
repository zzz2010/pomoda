import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.TreeMap;

import org.apache.commons.lang.ArrayUtils;
import org.jfree.util.ArrayUtilities;

import umontreal.iro.lecuyer.probdist.ChiSquareDist;




public class GapOptimalModelingThread extends GapBGModelingThread {

	
	static boolean ParaRecyclFlag=false;
	TreeMap<Double,ConstrainBlock> Pending=null;
	
	public GapOptimalModelingThread(HashMap<String, Double> gapmerCount,
			PWM sitePWM, HashSet<Integer> depend_pos, BGModel bg) {
		super(gapmerCount, sitePWM, depend_pos, bg);
		// TODO Auto-generated constructor stub
	}
	
	
	static TreeMap<Double,ConstrainBlock> getConstrainBlockQueue(double[] dmerProb, int colid)
	{
		
		
		
		return null;
		
	}


	
	
	
	static double[] findOptimalKconstrainEntries(double[] sorteddmerProb, int k)
	{
		double[] bound=new double[2];
		//KL_j=sum_{i=j}^{j+k}{P_i*log(P_i)}-k*m_j*log(m_j)
		//m_{j+1}=(k*m_j-P_j+P_{j+k+1})/k
		//KL_{j+1}=KL_j+P_{j+k+1}*log(P_{j+k+1})-P_j*log(P_j)+k*m_j*log(m_j)-k*m_{j+1}*log(m_{j+1})
		double last_KL=0;
		double last_m=0;
		double minKL=1000;
		int minIndex=0;
		for (int i = 0; i < sorteddmerProb.length-k+1; i++) {
			double KL=0;
			double m=0;
			if(i==0)
			{
				for (int j = 0; j < k; j++) {
					KL+=sorteddmerProb[j]*Math.log(sorteddmerProb[j]);
					m+=sorteddmerProb[j];
				}
				m/=k;
				KL-=k*m*Math.log(m);
				last_KL=KL;
				last_m=m;
				minKL=KL;
			}
			else
			{
				m=(last_m*k-sorteddmerProb[i-1]+sorteddmerProb[i+k-1])/k;
				KL=last_KL+sorteddmerProb[i+k-1]*Math.log(sorteddmerProb[i+k-1])-sorteddmerProb[i-1]*Math.log(sorteddmerProb[i-1])+k*last_m*Math.log(last_m)-k*m*Math.log(m);
				last_KL=KL;
				last_m=m;
				if(KL<minKL)
				{
					minKL=KL;
					minIndex=i;
				}
			}
		}
		
		bound[0]=sorteddmerProb[minIndex];
		bound[1]=sorteddmerProb[minIndex+k-1];
		return bound;
	}
	
	
	
	
	public void run() {
		
		
		
		//build Kmer model
		int dmerSize=1<<(depend_Pos.size() *2);
		
		
		dmerCount=new double[dmerSize];
		double[] Pt = new double[dmerSize];
		Integer[] sorteddpos=depend_Pos.toArray(new Integer[1]);
		Arrays.sort(sorteddpos);
		for (String gapstrj:gapmerCount.keySet()) {
			double weight=gapmerCount.get(gapstrj);
			String dmer="";
			
			if(sorteddpos[0]!=null)
			for (int k = 0; k < sorteddpos.length; k++) {
				int dpos=sorteddpos[k]-gapStart;
				dmer+=gapstrj.charAt(dpos);
			}
			if(dmer.indexOf('N')>-1)
				continue;
			if(depend_Pos.size()>1)
			{
			int hash=common.getHashing(dmer, 0, depend_Pos.size());
			dmerCount[hash]+=weight;
			}
			
			
		}
	
		//the number of free parameter is 3d not 4d-1
		int num_top=ParaNum*depend_Pos.size();//
		if(depend_Pos.size()>1)
		{
			dmerCount=common.Normalize(dmerCount);

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
				double weight=Pt[j]*(1-j*common.DoubleMinNormal);
				if(weight<0)
					weight=0;
				sorted_column.put(weight, j);
			}
		
		double sumprob=0;
		
		double[] sorteddmerProb=ArrayUtils.toPrimitive(sorted_column.keySet().toArray(new Double[1]));
		double[] bounds=findOptimalKconstrainEntries(sorteddmerProb, sorteddmerProb.length-num_top);
		if(depend_Pos.contains(5)&&depend_Pos.contains(0)&&depend_Pos.size()==2)
			sumprob=0;
		
		//sort the dmer by desc prob£¬ and take top 3d as the dependency model 
		for(Double key:sorted_column.descendingKeySet())
		{
			if(key>bounds[1]||key<bounds[0])
			{	
				DprobMap.put(common.Hash2ACGT(sorted_column.get(key), depend_Pos.size()), key);
				sumprob+=key;
			}

			if(DprobMap.size()==(num_top))
				break;
		}
	
		DprobMap.put("N", (1-sumprob)/(dmerSize-num_top));

		}
		
		//compute KL-divergence
		
		double sum_plogp_p=0;
		//the reason to use gapmer instead of dmer here, to measure, the thing is comparable among threads
		for (String gapstrj:gapmerCount.keySet()) {
			double p=gapmerCount.get(gapstrj);
			if(p==0.0)
				continue;
			
			
			double logq=0;
			if(depend_Pos.size()>1)//have dependency positions
			{
				String dmer="";
				//assume only one dependency group
				if(sorteddpos[0]!=null)
					for (int k = 0; k < sorteddpos.length; k++) {
						int dpos=sorteddpos[k]-gapStart;
						dmer+=gapstrj.charAt(dpos);
					}

				
				if(DprobMap.containsKey(dmer))
					logq+=Math.log(DprobMap.get(dmer));
				else
					logq+=Math.log(DprobMap.get("N"));
				
				int di=common.getHashing(dmer, 0, dmer.length()); 
				
				//compute the real prob, when considering bgprob inside the observed prob
				p=p*Pt[di]/dmerCount[di];

			}
			if(p==0.0)
				continue;
			//consider the non-dependency positions
			for (int k = 0; k < gapstrj.length(); k++) {
				int symid=common.acgt(gapstrj.charAt(k));
				if(depend_Pos.size()>1&&depend_Pos.contains(gapStart+k))
					continue;
				logq+=gapPWM.getLogWeight(k, symid);
				
			}
			sum_plogp_p+=p*(Math.log(p) -logq);
			//debuglist.add(Math.floor( Sites.size()*p)+Math.exp(logq));
		}
		
		KL_Divergence=sum_plogp_p;
		
		if(depend_Pos.size()==0)
		{
			KL_scorethresh=KL_Divergence;
		}
		if(KL_scorethresh<KL_Divergence)
		{
			//no need to consider
			//depend_Pos.clear();
			DprobMap.clear();

		}
		
		chisqPvalue=1-chisqTest(dmerCount,gapmerCount.size(),gapPWM,sorteddpos);
		//parameter recycling
		if(depend_Pos.size()>1&&chisqPvalue<0.05&&ParaRecyclFlag)
		{
			PendingParas=new TreeMap<Double, KeyValuePair<Integer,String>>();
			TreeMap<Double,KeyValuePair<Integer, String>> paraslist=GapImprover.getFreeParaQueue( dmerCount,-this.hashCode());
			int k=0;  //take |dmerCount|-topNum
			for(Double key:paraslist.keySet())
			{
				k++;
				if(paraslist.size()-PendingParas.size()==num_top)
					break;
				PendingParas.put(key, paraslist.get(key));
			}
			
		}
  
		
}
	
	

	 public static void main(String[] args) 
	{
		int dim=4;
		double[] testarray=new double[]{0.01,0.15,0.15,0.19,0.2,0.3};
		double[] bounds=findOptimalKconstrainEntries(testarray,3);
		System.out.println(Arrays.toString(bounds));
		
		
	}
	 
	 
	  

	
}
