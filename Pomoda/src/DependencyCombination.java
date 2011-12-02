import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.jgrapht.alg.BronKerboschCliqueFinder;
import org.jgrapht.graph.DefaultEdge;
import org.jgrapht.graph.SimpleGraph;


public class DependencyCombination {
	//allow multi dep-group in the same region, need to try different combinations
	public static HashMap<HashSet<Integer>,HashMap<String,Double>> FindBest_1(List<GapBGModelingThread> list)
	{

		HashMap<HashSet<Integer>,HashMap<String,Double>> Dmap=new HashMap<HashSet<Integer>,HashMap<String,Double>>();
		Iterator<GapBGModelingThread> iter3=list.iterator();
		int start=-1;
		double minKL=Double.MAX_VALUE;
		double baseScore=0;
		double bestScore=0;
		GapBGModelingThread bestThread=null;
		HashSet<Integer> bestDgroups=new HashSet<Integer>();
		try {

		while(iter3.hasNext())
		{
			GapBGModelingThread t1=iter3.next();	
			t1.join();
				if(t1.depend_Pos.size()==0)
				{
					System.out.println(t1.toString());
					baseScore=t1.KL_Divergence;
					//PWM case
					
				}
				else //debug
					System.out.println(t1.toString());
		}
		
		
		iter3=list.iterator();	
		//filter the negative threads, only consider the dep-group better than independence case
		ArrayList<GapBGModelingThread> positiveThread=new ArrayList<GapBGModelingThread>(list.size()/2);
		
		 LinkedList< HashSet<Integer> > queue=new LinkedList< HashSet<Integer> >();
		for(GapBGModelingThread t2:list)
		{
			
			if(baseScore>t2.KL_Divergence&&t2.chisqPvalue<0.05)
			{
				//I reuse the field KL_Divergence as a score, not the KL_Divergence meaning any more
				t2.KL_Divergence=baseScore-t2.KL_Divergence;
				positiveThread.add(t2);
				//higher the better now, KL_Divergence is score here
				if(t2.KL_Divergence>bestScore)
				{
					bestDgroups.clear();
					bestDgroups.add(positiveThread.size()-1);
					bestScore=t2.KL_Divergence;
				}
				
			}
		}
		//here: bestDgroups is the best only-1 dep-group
		
		// build graph
		SimpleGraph<GapBGModelingThread, DefaultEdge> graph=new SimpleGraph<GapBGModelingThread, DefaultEdge>(DefaultEdge.class);
		
		for (int i = 0; i < positiveThread.size(); i++) {
			GapBGModelingThread t1=positiveThread.get(i);
			t1.lamda=i;
			graph.addVertex(t1);
		}
		for (int i = 0; i < positiveThread.size()-1; i++) {
			GapBGModelingThread t1=positiveThread.get(i);
			HashSet<Integer> s1=t1.depend_Pos;
			for (int j = i+1; j < positiveThread.size(); j++) 
			{
				GapBGModelingThread t2=positiveThread.get(j);
				HashSet<Integer> s2=t2.depend_Pos;
				boolean overlap=false;
				HashSet<Integer> intersect=new HashSet<Integer>(s2);
				intersect.retainAll(s1);
				if(intersect.size()>0)
					overlap=true;
				
				//flexible combination
				if(overlap==false)
				{
					
					graph.addEdge(t1, t2);
					
				}
			}
		}
		//to here: queue contain all "2 dep-groups" set
		
		//enumerate all possible number of dep-group combinations, maximum clique
		BronKerboschCliqueFinder<GapBGModelingThread, DefaultEdge> cliquefinder=new BronKerboschCliqueFinder<GapBGModelingThread, DefaultEdge>(graph);
		Collection<Set<GapBGModelingThread>> cliques=cliquefinder.getAllMaximalCliques();
		System.out.println("Clique Number: "+cliques.size());
		for( Set<GapBGModelingThread> clique:cliques)
		{
			double score=0;
			HashSet<Integer> Dgroups=new HashSet<Integer>();
			for(GapBGModelingThread t : clique)
			{
				score+=t.KL_Divergence;
				Dgroups.add((int)t.lamda);
			}
			if(score>bestScore)
			{
				bestDgroups.clear();
				bestDgroups=Dgroups;
				bestScore=score;
			}
		}
		
		double descSum=0;
		if(bestDgroups.size()>0)
		for(Integer id : bestDgroups)
		{
			System.out.println("max KL desc:"+positiveThread.get(id).toString());
			descSum+=positiveThread.get(id).KL_Divergence;
			Dmap.put(positiveThread.get(id).depend_Pos,positiveThread.get(id).DprobMap);
		}
				

		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		
		return Dmap;
	}
	
	//DeltaKLMatrix the first column are zeros, meaning donate 0 parameter, no KL decrease
	public static double[] DP_computeMinKLdesc(double[][]DeltaKLMatrix)
	{
		int maxDonateNum=DeltaKLMatrix.length*(DeltaKLMatrix[0].length-1);
		double[] minKLdesc_dnoate_i=new double[maxDonateNum+1];
		double[][]DPtable=new double[DeltaKLMatrix.length][maxDonateNum+1];
		
		common.fill2DArray(DPtable, Double.MAX_VALUE);
		//initialize first column DP
		for (int i = 0; i < DeltaKLMatrix[0].length; i++) {
			DPtable[0][i]=DeltaKLMatrix[0][i];
		}
		
		
		for (int i = 1; i < DeltaKLMatrix.length; i++) { //for each position
			for (int donateNum = 0; donateNum < Math.min((i+1)*(DeltaKLMatrix[0].length-1)+1, maxDonateNum+1); donateNum++) {//for each possible donate number
				double minKLdesc=DPtable[i-1][donateNum]; // no need to donate for position i
			
				for (int j = Math.max(1,donateNum-i*(DeltaKLMatrix[0].length-1)); j <= Math.min(DeltaKLMatrix[0].length-1,donateNum); j++) {
					double KLdesc=DeltaKLMatrix[i][j]+DPtable[i-1][donateNum-j];
					if(minKLdesc>KLdesc)
					{
						minKLdesc=KLdesc;
					}
				}
				DPtable[i][donateNum]=minKLdesc;
			}
		}

		//start from 0 donation
		for (int i = 0; i < DPtable[0].length; i++) {
			minKLdesc_dnoate_i[i]=DPtable[DPtable.length-1][i];
		}
		
		//the first element in return is the min KL decrease when donate zero parameter
		return minKLdesc_dnoate_i;
		
	}
	
	public static HashMap<HashSet<Integer>,HashMap<String,Double>> FindBestCombination(List<GapOptimalModelingThread> list,
			HashMap<Integer,ArrayList<ConstrainBlock>> ConservedCBList,HashMap<Integer,ArrayList<ConstrainBlock>> DiverseCBList,double[][] m_matrix, int[] translate)
	{
		
		HashMap<HashSet<Integer>,HashMap<String,Double>> Dmap=new HashMap<HashSet<Integer>,HashMap<String,Double>>();
		Iterator<GapOptimalModelingThread> iter3=list.iterator();
		int start=-1;
		double minKL=Double.MAX_VALUE;
		double baseScore=0;
		double bestScore=0;
		GapBGModelingThread bestThread=null;
		HashSet<Integer> bestDgroups=new HashSet<Integer>();
		try {

		while(iter3.hasNext())
		{
			GapBGModelingThread t1=iter3.next();	
			t1.join();
				if(t1.depend_Pos.size()==0)
				{
					System.out.println(t1.toString());
					baseScore=t1.KL_Divergence;
					//PWM case
					
				}
				else //debug
					System.out.println(t1.toString());
		}
		
		
		iter3=list.iterator();	
		//filter the negative threads, only consider the dep-group better than independence case
		ArrayList<GapOptimalModelingThread> positiveThread=new ArrayList<GapOptimalModelingThread>(list.size()/2);
		
		 LinkedList< HashSet<Integer> > queue=new LinkedList< HashSet<Integer> >();
		for(GapOptimalModelingThread t2:list)
		{
			
			if(baseScore>t2.KL_Divergence&&t2.chisqPvalue<0.05)
			{
				//I reuse the field KL_Divergence as a score, not the KL_Divergence meaning any more
				t2.KL_Divergence=baseScore-t2.KL_Divergence;
				positiveThread.add(t2);
				//higher the better now, KL_Divergence is score here
				if(t2.KL_Divergence>bestScore)
				{
					bestDgroups.clear();
					bestDgroups.add(positiveThread.size()-1);
					bestScore=t2.KL_Divergence;
				}
				
			}
		}
		//here: bestDgroups is the best only-1 dep-group
		
		// build graph
		SimpleGraph<GapOptimalModelingThread, DefaultEdge> graph=new SimpleGraph<GapOptimalModelingThread, DefaultEdge>(DefaultEdge.class);
		
		for (int i = 0; i < positiveThread.size(); i++) {
			GapOptimalModelingThread t1=positiveThread.get(i);
			t1.lamda=i;
			graph.addVertex(t1);
		}
		for (int i = 0; i < positiveThread.size()-1; i++) {
			GapOptimalModelingThread t1=positiveThread.get(i);
			HashSet<Integer> s1=t1.depend_Pos;
			for (int j = i+1; j < positiveThread.size(); j++) 
			{
				GapOptimalModelingThread t2=positiveThread.get(j);
				HashSet<Integer> s2=t2.depend_Pos;
				boolean overlap=false;
				HashSet<Integer> intersect=new HashSet<Integer>(s2);
				intersect.retainAll(s1);
				if(intersect.size()>0)
					overlap=true;
				
				//flexible combination
				if(overlap==false)
				{		
					graph.addEdge(t1, t2);					
				}
			}
		}
		
		double [][] ConservedDeltaKL=new double[ConservedCBList.size()][4];

		int colid=0;
		for(ArrayList<ConstrainBlock> cblist:ConservedCBList.values())
		{
			ConservedDeltaKL[colid][0]=0;
			for (int donateNum = 1; donateNum < ConservedDeltaKL.length; donateNum++)
			{
			    ConservedDeltaKL[colid][donateNum]=cblist.get(donateNum-1).KL;
			}
			
			colid++;					
		}
		double[] min_conversedKL_descr=DP_computeMinKLdesc(ConservedDeltaKL);
		
		//enumerate all possible number of dep-group combinations, maximum clique
		BronKerboschCliqueFinder<GapOptimalModelingThread, DefaultEdge> cliquefinder=new BronKerboschCliqueFinder<GapOptimalModelingThread, DefaultEdge>(graph);
		Collection<Set<GapOptimalModelingThread>> cliques=cliquefinder.getAllMaximalCliques();
		System.out.println("Clique Number: "+cliques.size());
		for( Set<GapOptimalModelingThread> clique:cliques)
		{
				// compute the total KL for each clique, using parameter recycling
			HashSet<Integer> Dpos=new HashSet<Integer>();
			for(GapBGModelingThread t : clique)
			{		
				Dpos.addAll(t.depend_Pos);
			}
			
			//////////////////////////////////compute donate recourses///////////////////////
			double [][] IndColumnsDeltaKL=new double[DiverseCBList.size()-Dpos.size()][4];
			colid=0;
			for(Map.Entry<Integer, ArrayList<ConstrainBlock>>  elm:DiverseCBList.entrySet())
			{
				if(Dpos.contains(elm.getKey()))
					continue;
				IndColumnsDeltaKL[colid][0]=0;
				for (int donateNum = 1; donateNum <4; donateNum++)
				{
					IndColumnsDeltaKL[colid][donateNum]=elm.getValue().get(donateNum-1).KL;
				}			
				colid++;					
			}
			double[] min_DiversedKL_descr=DP_computeMinKLdesc(IndColumnsDeltaKL);	
			int totalDonateNumber=min_DiversedKL_descr.length+min_conversedKL_descr.length-1;
			double[] donateArray=new double[totalDonateNumber];
			//combine conserved and diverse
			double[] V1,V2;
			if(min_conversedKL_descr.length<min_DiversedKL_descr.length)
			{
				V1=min_conversedKL_descr; //short
				V2=min_DiversedKL_descr;  //long
			}
			else
			{
				V2=min_conversedKL_descr;
				V1=min_DiversedKL_descr;
			}
			for (int i = 0; i < donateArray.length; i++) {
				if(i==0)
				{
					donateArray[i]=0;
				}
				else
				{
					double mindeltaKL=Double.MAX_VALUE;
					for (int j = Math.max(0,i-V2.length+1); j < Math.min(i+1, V1.length); j++) {
						double deltaKL=V1[j]+V2[i-j];
						if(deltaKL<mindeltaKL)
						{
							mindeltaKL=deltaKL;
						}
					}
					donateArray[i]=mindeltaKL;
				}
			}
			//////////////////////////////////compute donate recourses///////////////////////
			System.out.println(Arrays.toString(donateArray));
			
		}
		
		double descSum=0;
		if(bestDgroups.size()>0)
			for(Integer id : bestDgroups)
			{
				System.out.println("max KL desc:"+positiveThread.get(id).toString());
				descSum+=positiveThread.get(id).KL_Divergence;
				Dmap.put(positiveThread.get(id).depend_Pos,positiveThread.get(id).DprobMap);
			}
				

		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		
		return Dmap;
		
		
		
	}
	
	
	 public static void main(String[] args) 
	{
		int dim=4;
		double[][] testarray=new double[3][4];
//		testarray[0]=new double[]{0,0,0};
//		testarray[1]=new double[]{0.1,1.1,2.1};
//		testarray[2]=new double[]{0.2,1.21,2.21};
//		testarray[3]=new double[]{0.3,1.32,2.32};
		
		testarray[0]=new double[]{0,0.1,0.2,0.3};
		testarray[1]=new double[]{0,1.1,1.21,1.32};
		testarray[2]=new double[]{0,2.1,2.21,2.32};
		double[] bounds=DP_computeMinKLdesc(testarray);
		System.out.println(Arrays.toString(bounds));
		
		
	}
	 
	
	
}
