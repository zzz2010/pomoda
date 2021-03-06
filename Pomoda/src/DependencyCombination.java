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

import EDU.oswego.cs.dl.util.concurrent.LinkedQueue;
import EDU.oswego.cs.dl.util.concurrent.PooledExecutor;


public class DependencyCombination extends Thread{
	
	static int threadNum=4;
	
	Set<GapOptimalModelingThread> clique;
	Set<GapOptimalModelingThread> bestclique;
	ArrayList<GapOptimalModelingThread> negThreads;
	double score=0;
	static double bestscore_allthreads=0;
	double bestscore=0;
	HashMap<Integer,ArrayList<ConstrainBlock>> DiverseCBList;
	 int[][] DivAction;
	 int[] DivIndAction;
	 int bestRecycNum=0;
	 double maxDeltaKL=0;
	 int[][] DepAction;
	 double[] min_conversedKL_descr;

	public DependencyCombination(Set<GapOptimalModelingThread> clique_,	HashMap<Integer,ArrayList<ConstrainBlock>> DiverseCBList_,double[]  min_conversedKL_descr_,ArrayList<GapOptimalModelingThread> negThreads_)
	{
		min_conversedKL_descr=min_conversedKL_descr_;
		clique=clique_;
		DiverseCBList=DiverseCBList_;
		negThreads=negThreads_;
		
	}
	
	
	static int cliqueScore2(Set<GapOptimalModelingThread> clique,	HashMap<Integer,ArrayList<ConstrainBlock>> DiverseCBList,double[]  min_conversedKL_descr, int[][] DivAction,int[] DivIndAction,int[][] DepAction)
	{
		
		HashSet<Integer> Dpos=new HashSet<Integer>();
		
		for(GapOptimalModelingThread t : clique)
		{		
			Dpos.addAll(t.depend_Pos);
		}
		
		//////////////////////////////////compute donate recourses///////////////////////
		double [][] IndColumnsDeltaKL=new double[DiverseCBList.size()-Dpos.size()][4];
		double[] donateArray=new double[1];

		if(DiverseCBList.size()-Dpos.size()>0)
		{
		int colid=0;
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
		
		 int maxDivIndDonateNum=IndColumnsDeltaKL.length*(IndColumnsDeltaKL[0].length-1);
		
		double[] min_DiversedKL_descr=DP_computeMinKLdesc(IndColumnsDeltaKL,DivAction);	
		int totalDonateNumber=min_DiversedKL_descr.length+min_conversedKL_descr.length-1;
		donateArray=new double[totalDonateNumber];
		//combine conserved and diverse
		double[] V1,V2;
		
		boolean swapFlag=false;
		if(min_conversedKL_descr.length<min_DiversedKL_descr.length)
		{
			V1=min_conversedKL_descr; //short
			V2=min_DiversedKL_descr;  //long
			swapFlag=true;
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
						DivIndAction[i]=j;
						if(swapFlag)
							DivIndAction[i]=i-j;
					}
				}
				donateArray[i]=mindeltaKL;
			}
		}
		}
		else
		{
			donateArray=min_conversedKL_descr;
		}
		//////////////////////////////////compute request recourses///////////////////////
		double[] maxKLincr=null;
		int maxDonateNum=donateArray.length+3*Dpos.size();
		//extend the donateArray for the 3*Dpos.size() free parameters in the front
		double[] newdonateArray=new double[maxDonateNum];
		for (int i = 3*Dpos.size()+1; i < newdonateArray.length; i++) {
			newdonateArray[i]=donateArray[i-3*Dpos.size()-1];
		}
		donateArray=newdonateArray;
		
		
		double score=0;
		HashSet<Integer> Dgroups=new HashSet<Integer>();
		int tid=0;
		for(GapOptimalModelingThread t : clique)
		{
			Dgroups.add((int)t.lamda);
			score+=t.KL_Divergence;
			double[] KLvec=new double[Math.min(maxDonateNum, t.PendingConstrainBlocks.size()) ];
			for (int i = 0; i < KLvec.length; i++) {
				KLvec[i]=t.PendingConstrainBlocks.get(i).KL;
			}
			if(maxKLincr==null)
			{
				maxKLincr=KLvec;
				for(int i = 0; i < KLvec.length; i++)
					DepAction[0][i]=i;
			}
			else
			{
				maxKLincr=DP_mergeMaxKLincr(KLvec, maxKLincr,DepAction[tid]);
			}
			tid++;
		}
		//find the maximum improvement
		double maxDeltaKL=0;
		int bestRecycNum=0;
		for (int i = 0; i < Math.min(maxDonateNum, maxKLincr.length); i++) {
			double DKL=maxKLincr[i]-donateArray[i]-i*common.DoubleMinNormal; //give penalty for small different.
			
			if(DKL>maxDeltaKL)
			{
				maxDeltaKL=DKL;
				bestRecycNum=i;
				
			}
			
		}
		
		return bestRecycNum;

	}
	
	static double cliqueScore(Set<GapOptimalModelingThread> clique,	HashMap<Integer,ArrayList<ConstrainBlock>> DiverseCBList,double[]  min_conversedKL_descr, double[] out_BenifitArr)
	{
		
		HashSet<Integer> Dpos=new HashSet<Integer>();
		
		for(GapOptimalModelingThread t : clique)
		{		
			Dpos.addAll(t.depend_Pos);
		}
		
		//////////////////////////////////compute donate recourses///////////////////////
		double [][] IndColumnsDeltaKL=new double[DiverseCBList.size()-Dpos.size()][4];
		double[] donateArray=new double[1];
		 int[][] DivAction=new int[1][1];
		 int[] DivIndAction=new int[1];
		if(DiverseCBList.size()-Dpos.size()>0)
		{
		int colid=0;
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
		
		 int maxDivIndDonateNum=IndColumnsDeltaKL.length*(IndColumnsDeltaKL[0].length-1);
		 DivAction=new int[IndColumnsDeltaKL.length][maxDivIndDonateNum+1];
		double[] min_DiversedKL_descr=DP_computeMinKLdesc(IndColumnsDeltaKL,DivAction);	
		int totalDonateNumber=min_DiversedKL_descr.length+min_conversedKL_descr.length-1;
		donateArray=new double[totalDonateNumber];
		//combine conserved and diverse
		double[] V1,V2;
		 DivIndAction=new int[totalDonateNumber];
		boolean swapFlag=false;
		if(min_conversedKL_descr.length<min_DiversedKL_descr.length)
		{
			V1=min_conversedKL_descr; //short
			V2=min_DiversedKL_descr;  //long
			swapFlag=true;
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
						DivIndAction[i]=j;
						if(swapFlag)
							DivIndAction[i]=i-j;
					}
				}
				donateArray[i]=mindeltaKL;
			}
		}
		}
		else
		{
			donateArray=min_conversedKL_descr;
		}
		//////////////////////////////////compute request recourses///////////////////////
		double[] maxKLincr=null;
		int maxDonateNum=donateArray.length+3*Dpos.size();
		//extend the donateArray for the 3*Dpos.size() free parameters in the front
		double[] newdonateArray=new double[maxDonateNum];
		for (int i = 3*Dpos.size()+1; i < newdonateArray.length; i++) {
			newdonateArray[i]=donateArray[i-3*Dpos.size()-1];
		}
		donateArray=newdonateArray;
		
		int[][] DepAction=new int[clique.size()][maxDonateNum];
		double score=0;
		HashSet<Integer> Dgroups=new HashSet<Integer>();
		int tid=0;
		for(GapOptimalModelingThread t : clique)
		{
			Dgroups.add((int)t.lamda);
			score+=t.KL_Divergence;
			double[] KLvec=new double[Math.min(maxDonateNum, t.PendingConstrainBlocks.size()) ];
			for (int i = 0; i < KLvec.length; i++) {
				KLvec[i]=t.PendingConstrainBlocks.get(i).KL;
			}
			if(maxKLincr==null)
			{
				maxKLincr=KLvec;
				for(int i = 0; i < KLvec.length; i++)
					DepAction[0][i]=i;
			}
			else
			{
				maxKLincr=DP_mergeMaxKLincr(KLvec, maxKLincr,DepAction[tid]);
			}
			tid++;
		}
		//find the maximum improvement
		double maxDeltaKL=0;
		int bestRecycNum=0;
		for (int i = 0; i < Math.min(maxDonateNum, maxKLincr.length); i++) {
			double DKL=maxKLincr[i]-donateArray[i]-i*common.DoubleMinNormal; //give penalty for small different.
			out_BenifitArr[i]=DKL;
			if(DKL>maxDeltaKL)
			{
				maxDeltaKL=DKL;
				bestRecycNum=i;
				
			}
			
		}
		
		return score+maxDeltaKL;

	}
	
	
	public void run() {
		
		HashSet<Integer> Dpos=new HashSet<Integer>();
		int pruningCount=0;
		double noParaBasicScore=0;
		for(GapOptimalModelingThread t : clique)
		{		
			Dpos.addAll(t.depend_Pos);
			noParaBasicScore=t.KL_Divergence;
		}
		double[] basicBenifitArr=new double[3*DiverseCBList.size()+min_conversedKL_descr.length];
		double basicCliquescore=cliqueScore(clique, DiverseCBList, min_conversedKL_descr,basicBenifitArr);
		bestscore=basicCliquescore;
		bestclique=clique;
		maxDeltaKL=bestscore-noParaBasicScore;
		HashSet<Set<GapOptimalModelingThread>> triedClique=new HashSet<Set<GapOptimalModelingThread>>(negThreads.size());
		LinkedList< Set<GapOptimalModelingThread> > queue=new LinkedList< Set<GapOptimalModelingThread> >();
		queue.add(new HashSet<GapOptimalModelingThread>(clique));
		Set<GapOptimalModelingThread> topElm=null;
		while((topElm=queue.poll())!=null)
		{
			Dpos.clear();
			for(GapOptimalModelingThread t : topElm)
			{		
				Dpos.addAll(t.depend_Pos);
			}

			for(GapOptimalModelingThread t1 : negThreads)
			{
				HashSet<Integer> s2=t1.depend_Pos;
				boolean overlap=false;
				HashSet<Integer> intersect=new HashSet<Integer>(s2);
				intersect.retainAll(Dpos);
				if(intersect.size()>0)
					overlap=true;			
				//flexible combination
				if(overlap==false)
				{
					double upperbound=basicCliquescore+t1.PendingConstrainBlocks.get(t1.PendingConstrainBlocks.size()-1).KL+t1.KL_Divergence;
					HashSet<GapOptimalModelingThread> newclique=new HashSet<GapOptimalModelingThread>(topElm);
					newclique.add(t1);
					if(triedClique.contains(newclique))
						continue;
					else
						triedClique.add(newclique);
					
					queue.add(newclique);
					
					double noParaScore=noParaBasicScore+t1.KL_Divergence;
					topElm.removeAll(clique);
					for(GapOptimalModelingThread t2:topElm)//estimate upperbound for all negative threads
					{
						upperbound+=t2.PendingConstrainBlocks.get(t2.PendingConstrainBlocks.size()-1).KL+t2.KL_Divergence;
						noParaScore+=t2.KL_Divergence;
					}
					if(upperbound<bestscore_allthreads)
					{
						pruningCount++;
						continue;
					}
					double[] BenifitArr=new double[3*DiverseCBList.size()+min_conversedKL_descr.length];
					double newcliqueScore=cliqueScore(newclique, DiverseCBList, min_conversedKL_descr,BenifitArr);
					if(newcliqueScore>bestscore)
					{
						bestscore=newcliqueScore;
						bestclique=newclique;
						maxDeltaKL=bestscore-noParaScore;
						if(newcliqueScore>bestscore_allthreads)
						{
							bestscore_allthreads=newcliqueScore;
						}
					}	
					

				}			
			}
		}
		
	//	 System.out.println("BreathSearch Count:"+pruningCount);
	}

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
			
			if(baseScore>t2.KL_Divergence&&t2.chisqPvalue<GapImprover.FDR)
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
	
	public static double[] DP_mergeMaxKLincr(double[]V1,double[]V2,int[] V1Action)
	{
		double[] mergedV=new double[V1Action.length];
		for (int i = 0; i < mergedV.length; i++) {

				double maxdeltaKL=Double.NEGATIVE_INFINITY;
				for (int j = Math.max(0,i-V2.length+1); j < Math.min(i+1, V1.length); j++) {
					double deltaKL=V1[j]+V2[i-j];
					if(deltaKL>maxdeltaKL)
					{
						maxdeltaKL=deltaKL;
						V1Action[i]=j;
					}
				}
				mergedV[i]=maxdeltaKL;
			}
		
		
		
		
		return mergedV;
	}
	
	//DeltaKLMatrix the first column are zeros, meaning donate 0 parameter, no KL decrease
	public static double[] DP_computeMinKLdesc(double[][]DeltaKLMatrix,int[][] actionMatrix)
	{
		if(DeltaKLMatrix.length==0)
			return new double[1];
		int maxDonateNum=DeltaKLMatrix.length*(DeltaKLMatrix[0].length-1);
		double[] minKLdesc_dnoate_i=new double[maxDonateNum+1];
		double[][]DPtable=new double[DeltaKLMatrix.length][maxDonateNum+1];
		
		common.fill2DArray(DPtable, Double.MAX_VALUE);
		//initialize first column DP
		for (int i = 0; i < DeltaKLMatrix[0].length; i++) {
			DPtable[0][i]=DeltaKLMatrix[0][i];
			actionMatrix[0][i]=i;
		}
		
		
		
		for (int i = 1; i < DeltaKLMatrix.length; i++) { //for each position
			for (int donateNum = 0; donateNum < Math.min((i+1)*(DeltaKLMatrix[0].length-1)+1, maxDonateNum+1); donateNum++) {//for each possible donate number
				double minKLdesc=DPtable[i-1][donateNum]; // no need to donate for position i
			
				for (int j = Math.max(1,donateNum-i*(DeltaKLMatrix[0].length-1)); j <= Math.min(DeltaKLMatrix[0].length-1,donateNum); j++) {
					double KLdesc=DeltaKLMatrix[i][j]+DPtable[i-1][donateNum-j];
					if(minKLdesc>KLdesc)
					{
						actionMatrix[i][donateNum]=j;
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
	
	
	
	public static int[] backtracking_DP(int[][]Actions, int finalRequestNum)
	{
		if(Actions.length==0)
			return new int[1];
		int[] actArray=new int[Actions.length];
		int restRequestNum=finalRequestNum-Actions[Actions.length-1][finalRequestNum];
		actArray[actArray.length-1]=Actions[Actions.length-1][finalRequestNum];
		for (int i = 1; i < Actions.length; i++) {
			actArray[actArray.length-i-1]=Actions[Actions.length-i-1][restRequestNum];
			restRequestNum-=actArray[actArray.length-i-1];
		}
		
		
		return actArray;
	}
	
	//no donator for dep-group
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
		Set<GapOptimalModelingThread> bestDgroups=new HashSet<GapOptimalModelingThread>();
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
			
			if(baseScore>t2.KL_Divergence&&t2.chisqPvalue<GapImprover.FDR)
			{
				//I reuse the field KL_Divergence as a score, not the KL_Divergence meaning any more
				t2.KL_Divergence=baseScore-t2.KL_Divergence;
				positiveThread.add(t2);
				//higher the better now, KL_Divergence is score here
				if(t2.KL_Divergence>bestScore)
				{
					bestDgroups.clear();
					bestDgroups.add(t2);
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
		double recyclingEnhance=0;
		int colid=0;
		for(ArrayList<ConstrainBlock> cblist:ConservedCBList.values())
		{
			ConservedDeltaKL[colid][0]=0;
			for (int donateNum = 1; donateNum < ConservedDeltaKL[0].length; donateNum++)
			{
			    ConservedDeltaKL[colid][donateNum]=cblist.get(donateNum-1).KL;
			}
			
			colid++;					
		}
		
		int[][] bestConAction=new int[1][1],bestDivIndAction=new int[1][1],BestDepAction = new int[1][1];
		int[] bestDivIndFinalAction=new int[1];
		int bestReqNum=0;
		int maxConvDonateNum=0;
		if(ConservedCBList.size()>0)
		maxConvDonateNum=ConservedDeltaKL.length*(ConservedDeltaKL[0].length-1);
		bestConAction=new int[ConservedDeltaKL.length][maxConvDonateNum+1];
		double[] min_conversedKL_descr=DP_computeMinKLdesc(ConservedDeltaKL,bestConAction);
		
		//enumerate all possible number of dep-group combinations, maximum clique
		BronKerboschCliqueFinder<GapOptimalModelingThread, DefaultEdge> cliquefinder=new BronKerboschCliqueFinder<GapOptimalModelingThread, DefaultEdge>(graph);
		Collection<Set<GapOptimalModelingThread>> cliques=cliquefinder.getAllMaximalCliques();
		System.out.println("Clique Number: "+cliques.size());
		
		// no dependency group found
		if(positiveThread.size()==0)
			return Dmap;
		//when only have one positive thread, will report no clique
		if(cliques.size()==0)
		{
			HashSet<GapOptimalModelingThread> hs=new HashSet<GapOptimalModelingThread>();
			hs.add(positiveThread.get(0));
			cliques.add(hs);
		}
		
		for( Set<GapOptimalModelingThread> clique:cliques)
		{
				// compute the total KL for each clique, using parameter recycling


		
			
			HashSet<Integer> Dpos=new HashSet<Integer>();
			
			for(GapOptimalModelingThread t : clique)
			{		
				Dpos.addAll(t.depend_Pos);
			}
			
			//////////////////////////////////compute donate recourses///////////////////////
			double [][] IndColumnsDeltaKL=new double[DiverseCBList.size()-Dpos.size()][4];
			double[] donateArray=new double[1];
			 int[][] DivAction=new int[1][1];
			 int[] DivIndAction=new int[1];
			if(DiverseCBList.size()-Dpos.size()>0)
			{
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
			
			 int maxDivIndDonateNum=IndColumnsDeltaKL.length*(IndColumnsDeltaKL[0].length-1);
			 DivAction=new int[IndColumnsDeltaKL.length][maxDivIndDonateNum+1];
			double[] min_DiversedKL_descr=DP_computeMinKLdesc(IndColumnsDeltaKL,DivAction);	
			int totalDonateNumber=min_DiversedKL_descr.length+min_conversedKL_descr.length-1;
			donateArray=new double[totalDonateNumber];
			//combine conserved and diverse
			double[] V1,V2;
			 DivIndAction=new int[totalDonateNumber];
			boolean swapFlag=false;
			if(min_conversedKL_descr.length<min_DiversedKL_descr.length)
			{
				V1=min_conversedKL_descr; //short
				V2=min_DiversedKL_descr;  //long
				swapFlag=true;
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
							DivIndAction[i]=j;
							if(swapFlag)
								DivIndAction[i]=i-j;
						}
					}
					donateArray[i]=mindeltaKL;
				}
			}
			}
			else
			{
				donateArray=min_conversedKL_descr;
			}
			//////////////////////////////////compute request recourses///////////////////////
			double[] maxKLincr=null;
			int[][] DepAction=new int[clique.size()][donateArray.length];
			double score=0;
			HashSet<Integer> Dgroups=new HashSet<Integer>();
			int tid=0;
			for(GapOptimalModelingThread t : clique)
			{
				Dgroups.add((int)t.lamda);
				score+=t.KL_Divergence;
				double[] KLvec=new double[Math.min(donateArray.length, t.PendingConstrainBlocks.size()+1) ];
				for (int i = 1; i < KLvec.length; i++) {
					KLvec[i]=t.PendingConstrainBlocks.get(i-1).KL;
				}
				if(maxKLincr==null)
				{
					maxKLincr=KLvec;
					for(int i = 0; i < KLvec.length; i++)
						DepAction[0][i]=i;
				}
				else
				{
					maxKLincr=DP_mergeMaxKLincr(KLvec, maxKLincr,DepAction[tid]);
				}
				tid++;
			}
			//find the maximum improvement
			double maxDeltaKL=0;
			int bestRecycNum=0;
			for (int i = 0; i < Math.min(donateArray.length, maxKLincr.length); i++) {
				double DKL=maxKLincr[i]-donateArray[i]-i*common.DoubleMinNormal; //give penalty for small different.
				if(DKL>maxDeltaKL)
				{
					maxDeltaKL=DKL;
					bestRecycNum=i;
					
				}
				
			}
			score+=maxDeltaKL;
			if(score>bestScore)
			{
				bestDgroups.clear();
				bestDgroups=clique;
				recyclingEnhance=maxDeltaKL;
				bestScore=score;
				bestDivIndAction=DivAction;
				bestDivIndFinalAction=DivIndAction;
				BestDepAction=DepAction;
				bestReqNum=bestRecycNum;
			}
			}
			
		
		
		///////////////////////////////decode action//////////////////////////
		if(bestDgroups.size()>0&&bestReqNum>0) 
		{
			int[] bt_depAction=backtracking_DP(BestDepAction, bestReqNum);
			HashSet<Integer> tPosSet=new HashSet<Integer>();
			int depId=0;
			for(GapOptimalModelingThread t : bestDgroups)
			{
				tPosSet.addAll(t.depend_Pos);
				if(bt_depAction[depId]>0)
				{
					ConstrainBlock CB=t.PendingConstrainBlocks.get(bt_depAction[depId]-1);
					t.updateParaNum(t.dmerCount,new double[]{CB.lowerbound,CB.upperbound});
				}
				depId++;
			}
			//change Independent column
			int totaldivcolReq=0;
			if(bestDivIndFinalAction.length>bestReqNum)
				totaldivcolReq=bestDivIndFinalAction[bestReqNum];
			//conserved bases
			int[] bt_conAction=backtracking_DP(bestConAction,bestReqNum-totaldivcolReq);
			int cid=-1;
			for(Map.Entry<Integer, ArrayList<ConstrainBlock>>  elm:ConservedCBList.entrySet())
			{
				cid++;
				if(bt_conAction[cid]==0)
					continue;
				int orignalColumn=elm.getKey();
				HashMap<String,Double> dprobMap=new HashMap<String,Double>();
				ConstrainBlock  cb=elm.getValue().get(bt_conAction[cid]-1);
				double sumprob=0;
				int taken=0;
				HashSet<Integer> dPos=new HashSet<Integer>();
				dPos.add(orignalColumn);
				for (int i = 0; i < 4; i++) {
					double prob=m_matrix[orignalColumn][i];
					if(prob>(cb.upperbound+common.DoubleMinNormal)||prob<cb.lowerbound)
					{
						sumprob+=m_matrix[orignalColumn][i];
						dprobMap.put(common.Hash2ACGT(i, 1), m_matrix[orignalColumn][i]);
						taken++;
					}
				}
				dprobMap.put("N",(1-sumprob)/(4-taken));
				System.out.println("recycling :"+dPos);
				Dmap.put(dPos, dprobMap);
				
			}
			//diverse bases
			cid=-1;
			int[] bt_divAction=backtracking_DP(bestDivIndAction, totaldivcolReq);
			for(Map.Entry<Integer, ArrayList<ConstrainBlock>>  elm:DiverseCBList.entrySet())
			{
				int transColumn=elm.getKey();
				if(tPosSet.contains(transColumn))
				{
					continue;
				}
				cid++;
	
				if(bt_divAction[cid]==0)
					continue;
				
				
				int orignalColumn=translate[transColumn];
				
				HashMap<String,Double> dprobMap=new HashMap<String,Double>();
				ConstrainBlock  cb=elm.getValue().get(bt_divAction[cid]-1);
				double sumprob=0;
				int taken=0;
				HashSet<Integer> dPos=new HashSet<Integer>();
				dPos.add(orignalColumn);
				for (int i = 0; i < 4; i++) {
					double prob=m_matrix[orignalColumn][i];
					if(prob>cb.upperbound+common.DoubleMinNormal||prob<cb.lowerbound)
					{
						sumprob+=m_matrix[orignalColumn][i];
						dprobMap.put(common.Hash2ACGT(i, 1), m_matrix[orignalColumn][i]);
						taken++;
					}
				}
				dprobMap.put("N",(1-sumprob)/(4-taken));
				System.out.println("recycling :"+dPos);
				Dmap.put(dPos, dprobMap);
				
			}
			System.out.println("recycling enhance KL:"+recyclingEnhance);
		}
		
		
		
		double descSum=0;
		if(bestDgroups.size()>0)
			for(GapOptimalModelingThread t : bestDgroups)
			{
				System.out.println("max KL desc:"+t.toString());
				descSum+=t.KL_Divergence;
				Dmap.put(t.depend_Pos,t.DprobMap);
			}
				

		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		
		return Dmap;
		
		
		
	}
	
	//dep-group can be donator
	public static HashMap<HashSet<Integer>,HashMap<String,Double>> FindBestCombination2(List<GapOptimalModelingThread> list,
			HashMap<Integer,ArrayList<ConstrainBlock>> ConservedCBList,HashMap<Integer,ArrayList<ConstrainBlock>> DiverseCBList,double[][] m_matrix, int[] translate)
	{
		
		HashMap<HashSet<Integer>,HashMap<String,Double>> Dmap=new HashMap<HashSet<Integer>,HashMap<String,Double>>();
		Iterator<GapOptimalModelingThread> iter3=list.iterator();
		int start=-1;
		double minKL=Double.MAX_VALUE;
		double baseScore=0;
		double bestScore=0;
		GapBGModelingThread bestThread=null;
		Set<GapOptimalModelingThread> bestDgroups=new HashSet<GapOptimalModelingThread>();
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
		
		
		for(GapOptimalModelingThread t2:list)
		{
			
			if(baseScore>t2.KL_Divergence&&t2.chisqPvalue<GapImprover.FDR)
			{
				//I reuse the field KL_Divergence as a score, not the KL_Divergence meaning any more
				t2.KL_Divergence=baseScore-t2.KL_Divergence;
				positiveThread.add(t2);
				//higher the better now, KL_Divergence is score here
				if(t2.KL_Divergence>bestScore)
				{
					bestDgroups.clear();
					bestDgroups.add(t2);
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
		double recyclingEnhance=0;
		int colid=0;
		for(ArrayList<ConstrainBlock> cblist:ConservedCBList.values())
		{
			ConservedDeltaKL[colid][0]=0;
			for (int donateNum = 1; donateNum < ConservedDeltaKL[0].length; donateNum++)
			{
			    ConservedDeltaKL[colid][donateNum]=cblist.get(donateNum-1).KL;
			}
			
			colid++;					
		}
		
		int[][] bestConAction=new int[1][1],bestDivIndAction=new int[1][1],BestDepAction = new int[1][1];
		int[] bestDivIndFinalAction=new int[1];
		int bestReqNum=0;
		int maxConvDonateNum=0;
		if(ConservedCBList.size()>0)
		maxConvDonateNum=ConservedDeltaKL.length*(ConservedDeltaKL[0].length-1);
		bestConAction=new int[ConservedDeltaKL.length][maxConvDonateNum+1];
		double[] min_conversedKL_descr=DP_computeMinKLdesc(ConservedDeltaKL,bestConAction);
		
		//enumerate all possible number of dep-group combinations, maximum clique
		BronKerboschCliqueFinder<GapOptimalModelingThread, DefaultEdge> cliquefinder=new BronKerboschCliqueFinder<GapOptimalModelingThread, DefaultEdge>(graph);
		Collection<Set<GapOptimalModelingThread>> cliques=cliquefinder.getAllMaximalCliques();
		System.out.println("Clique Number: "+cliques.size());
		// no dependency group found
		if(positiveThread.size()==0)
			return Dmap;
		//when only have one positive thread, will report no clique
		if(cliques.size()==0)
		{
			HashSet<GapOptimalModelingThread> hs=new HashSet<GapOptimalModelingThread>();
			hs.add(positiveThread.get(0));
			cliques.add(hs);
		}
		
		
		//Threads Pool
		LinkedList<DependencyCombination> threadPool=new LinkedList<DependencyCombination>();
		PooledExecutor executor = new PooledExecutor(new LinkedQueue());
		executor.setMinimumPoolSize(1);
		executor.setKeepAliveTime(-1);
		if(cliques.size()>100&&threadNum>1)
			executor.createThreads(threadNum);
		for( Set<GapOptimalModelingThread> clique:cliques)
		{
				// compute the total KL for each clique, using parameter recycling

			
			HashSet<Integer> Dpos=new HashSet<Integer>();
			
			for(GapOptimalModelingThread t : clique)
			{		
				Dpos.addAll(t.depend_Pos);
			}
			
			//////////////////////////////////compute donate recourses///////////////////////
			double [][] IndColumnsDeltaKL=new double[DiverseCBList.size()-Dpos.size()][4];
			double[] donateArray=new double[1];
			 int[][] DivAction=new int[1][1];
			 int[] DivIndAction=new int[1];
			if(DiverseCBList.size()-Dpos.size()>0)
			{
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
			
			 int maxDivIndDonateNum=IndColumnsDeltaKL.length*(IndColumnsDeltaKL[0].length-1);
			 DivAction=new int[IndColumnsDeltaKL.length][maxDivIndDonateNum+1];
			double[] min_DiversedKL_descr=DP_computeMinKLdesc(IndColumnsDeltaKL,DivAction);	
			int totalDonateNumber=min_DiversedKL_descr.length+min_conversedKL_descr.length-1;
			donateArray=new double[totalDonateNumber];
			//combine conserved and diverse
			double[] V1,V2;
			 DivIndAction=new int[totalDonateNumber];
			boolean swapFlag=false;
			if(min_conversedKL_descr.length<min_DiversedKL_descr.length)
			{
				V1=min_conversedKL_descr; //short
				V2=min_DiversedKL_descr;  //long
				swapFlag=true;
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
							DivIndAction[i]=j;
							if(swapFlag)
								DivIndAction[i]=i-j;
						}
					}
					donateArray[i]=mindeltaKL;
				}
			}
			}
			else
			{
				donateArray=min_conversedKL_descr;
			}
			//////////////////////////////////compute request recourses///////////////////////
			double[] maxKLincr=null;
			int maxDonateNum=donateArray.length+3*Dpos.size();
			//extend the donateArray for the 3*Dpos.size() free parameters in the front
			double[] newdonateArray=new double[maxDonateNum];
			for (int i = 3*Dpos.size()+1; i < newdonateArray.length; i++) {
				newdonateArray[i]=donateArray[i-3*Dpos.size()-1];
			}
			donateArray=newdonateArray;
			
			int[][] DepAction=new int[clique.size()][maxDonateNum];
			double score=0;
			HashSet<Integer> Dgroups=new HashSet<Integer>();
			int tid=0;
			for(GapOptimalModelingThread t : clique)
			{
				Dgroups.add((int)t.lamda);
				score+=t.KL_Divergence;
				double[] KLvec=new double[Math.min(maxDonateNum, t.PendingConstrainBlocks.size()) ];
				for (int i = 0; i < KLvec.length; i++) {
					KLvec[i]=t.PendingConstrainBlocks.get(i).KL;
				}
				if(maxKLincr==null)
				{
					maxKLincr=KLvec;
					for(int i = 0; i < KLvec.length; i++)
						DepAction[0][i]=i;
				}
				else
				{
					maxKLincr=DP_mergeMaxKLincr(KLvec, maxKLincr,DepAction[tid]);
				}
				tid++;
			}
			//find the maximum improvement
			double maxDeltaKL=0;
			int bestRecycNum=0;
			for (int i = 0; i < Math.min(maxDonateNum, maxKLincr.length); i++) {
				double DKL=maxKLincr[i]-donateArray[i]-i*common.DoubleMinNormal; //give penalty for small different.
				if(DKL>maxDeltaKL)
				{
					maxDeltaKL=DKL;
					bestRecycNum=i;
					
				}
				
			}
			score+=maxDeltaKL;
			if(score>bestScore)
			{
				bestDgroups.clear();
				bestDgroups=clique;
				recyclingEnhance=maxDeltaKL;
				bestScore=score;
				bestDivIndAction=DivAction;
				bestDivIndFinalAction=DivIndAction;
				BestDepAction=DepAction;
				bestReqNum=bestRecycNum;
			}
		}
			

		
		///////////////////////////////decode action//////////////////////////
		if(bestDgroups.size()>0&&bestReqNum>0) 
		{
			int[] bt_depAction=backtracking_DP(BestDepAction, bestReqNum);
			HashSet<Integer> tPosSet=new HashSet<Integer>();
			int depId=0;
			for(GapOptimalModelingThread t : bestDgroups)
			{
				tPosSet.addAll(t.depend_Pos);
				if(bt_depAction[depId]>0)
				{
					ConstrainBlock CB=t.PendingConstrainBlocks.get(bt_depAction[depId]);
					t.updateParaNum(t.dmerCount,new double[]{CB.lowerbound,CB.upperbound});
				}
				depId++;
			}
			//change Independent column
			int totaldivcolReq=0;
			if(bestDivIndFinalAction.length>bestReqNum-3*tPosSet.size())
				totaldivcolReq=bestDivIndFinalAction[bestReqNum-3*tPosSet.size()];
			//conserved bases
			int[] bt_conAction=backtracking_DP(bestConAction,bestReqNum-totaldivcolReq-3*tPosSet.size());
			int cid=-1;
			for(Map.Entry<Integer, ArrayList<ConstrainBlock>>  elm:ConservedCBList.entrySet())
			{
				cid++;
				if(bt_conAction[cid]==0)
					continue;
				int orignalColumn=elm.getKey();
				HashMap<String,Double> dprobMap=new HashMap<String,Double>();
				ConstrainBlock  cb=elm.getValue().get(bt_conAction[cid]-1);
				double sumprob=0;
				int taken=0;
				HashSet<Integer> dPos=new HashSet<Integer>();
				dPos.add(orignalColumn);
				for (int i = 0; i < 4; i++) {
					double prob=m_matrix[orignalColumn][i];
					if(prob>(cb.upperbound+common.DoubleMinNormal)||prob<cb.lowerbound)
					{
						if(m_matrix[orignalColumn][i]==1)
							m_matrix[orignalColumn][i]=1-common.DoubleMinNormal;
						sumprob+=m_matrix[orignalColumn][i];
						dprobMap.put(common.Hash2ACGT(i, 1), m_matrix[orignalColumn][i]);
						taken++;
					}
				}
				dprobMap.put("N",(1-sumprob)/(4-taken));
				System.out.println("recycling :"+dPos);
				Dmap.put(dPos, dprobMap);
				
			}
			//diverse bases
			cid=-1;
			int[] bt_divAction=backtracking_DP(bestDivIndAction, totaldivcolReq);
			for(Map.Entry<Integer, ArrayList<ConstrainBlock>>  elm:DiverseCBList.entrySet())
			{
				int transColumn=elm.getKey();
				if(tPosSet.contains(transColumn))
				{
					continue;
				}
				cid++;
	
				if(bt_divAction[cid]==0)
					continue;
				
				
				int orignalColumn=translate[transColumn];
				
				HashMap<String,Double> dprobMap=new HashMap<String,Double>();
				ConstrainBlock  cb=elm.getValue().get(bt_divAction[cid]-1);
				double sumprob=0;
				int taken=0;
				HashSet<Integer> dPos=new HashSet<Integer>();
				dPos.add(orignalColumn);
				for (int i = 0; i < 4; i++) {
					double prob=m_matrix[orignalColumn][i];
					if(prob>cb.upperbound+common.DoubleMinNormal||prob<cb.lowerbound)
					{
						if(m_matrix[orignalColumn][i]==1)
							m_matrix[orignalColumn][i]=1-+common.DoubleMinNormal;
						sumprob+=m_matrix[orignalColumn][i];
						dprobMap.put(common.Hash2ACGT(i, 1), m_matrix[orignalColumn][i]);
						taken++;
					}
				}
				dprobMap.put("N",(1-sumprob)/(4-taken));
				System.out.println("recycling :"+dPos);
				Dmap.put(dPos, dprobMap);				
			}
			System.out.println("recycling enhance KL:"+recyclingEnhance);
		}
		
		
		
		double descSum=0;
		if(bestDgroups.size()>0)
			for(GapOptimalModelingThread t : bestDgroups)
			{
				System.out.println("max KL desc:"+t.toString());
				descSum+=t.KL_Divergence;
				Dmap.put(t.depend_Pos,t.DprobMap);
			}
				

		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		
		return Dmap;
		
		
		
	}
	
	//allow recyle the negative thread, but must based on maximum clique group
	public static HashMap<HashSet<Integer>,HashMap<String,Double>> FindBestCombination3(List<GapOptimalModelingThread> list,
			HashMap<Integer,ArrayList<ConstrainBlock>> ConservedCBList,HashMap<Integer,ArrayList<ConstrainBlock>> DiverseCBList,double[][] m_matrix, int[] translate)
	{
		
		HashMap<HashSet<Integer>,HashMap<String,Double>> Dmap=new HashMap<HashSet<Integer>,HashMap<String,Double>>();
		Iterator<GapOptimalModelingThread> iter3=list.iterator();
		int start=-1;
		double minKL=Double.MAX_VALUE;
		double baseScore=0;
		double bestScore=0;
		GapBGModelingThread bestThread=null;
		Set<GapOptimalModelingThread> bestDgroups=new HashSet<GapOptimalModelingThread>();
		try {

		while(iter3.hasNext())
		{
			GapBGModelingThread t1=iter3.next();	
//			if(t1.depend_Pos.size()==2&&t1.depend_Pos.contains(0)&&t1.depend_Pos.contains(1))
//				bestScore=0;
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
		ArrayList<GapOptimalModelingThread> negativeThread=new ArrayList<GapOptimalModelingThread>(list.size()/2);
		
		 LinkedList< HashSet<Integer> > queue=new LinkedList< HashSet<Integer> >();
		for(GapOptimalModelingThread t2:list)
		{
			if(t2.chisqPvalue<GapImprover.FDR)
			{
				t2.KL_Divergence=baseScore-t2.KL_Divergence;
				if(t2.KL_Divergence>0)
				{
					//I reuse the field KL_Divergence as a score, not the KL_Divergence meaning any more					
					positiveThread.add(t2);
					//higher the better now, KL_Divergence is score here
					if(t2.KL_Divergence>bestScore)
					{
						bestDgroups.clear();
						bestDgroups.add(t2);
						bestScore=t2.KL_Divergence;
					}
					
				}
				else if(t2.depend_Pos.size()>0)
				{						
					negativeThread.add(t2);
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
		double recyclingEnhance=0;
		int colid=0;
		for(ArrayList<ConstrainBlock> cblist:ConservedCBList.values())
		{
			ConservedDeltaKL[colid][0]=0;
			for (int donateNum = 1; donateNum < ConservedDeltaKL[0].length; donateNum++)
			{
			    ConservedDeltaKL[colid][donateNum]=cblist.get(donateNum-1).KL;
			}
			
			colid++;					
		}
		
		int[][] bestConAction=new int[1][1],bestDivIndAction=new int[1][1],BestDepAction = new int[1][1];
		int[] bestDivIndFinalAction=new int[1];
		int bestReqNum=0;
		int maxConvDonateNum=0;
		if(ConservedCBList.size()>0)
			maxConvDonateNum=ConservedDeltaKL.length*(ConservedDeltaKL[0].length-1);
		bestConAction=new int[ConservedDeltaKL.length][maxConvDonateNum+1];
		double[] min_conversedKL_descr=DP_computeMinKLdesc(ConservedDeltaKL,bestConAction);
		
		//enumerate all possible number of dep-group combinations, maximum clique
		BronKerboschCliqueFinder<GapOptimalModelingThread, DefaultEdge> cliquefinder=new BronKerboschCliqueFinder<GapOptimalModelingThread, DefaultEdge>(graph);
		Collection<Set<GapOptimalModelingThread>> cliques=cliquefinder.getAllMaximalCliques();
		System.out.println("Clique Number: "+cliques.size());
		
		System.out.println("number of positive:"+positiveThread.size());
		// no dependency group found
		if(positiveThread.size()==0)
			return Dmap;
		
		//when only have one positive thread, will report no clique
		if(cliques.size()==0)
		{
			HashSet<GapOptimalModelingThread> hs=new HashSet<GapOptimalModelingThread>();
			hs.add(positiveThread.get(0));
			cliques.add(hs);
		}
		
		
		//Threads Pool
		LinkedList<DependencyCombination> threadPool=new LinkedList<DependencyCombination>();
		PooledExecutor executor = new PooledExecutor(new LinkedQueue());
		executor.setMinimumPoolSize(threadNum);
		executor.setKeepAliveTime(1000 * 60*5);
		System.out.println("Number of Negative Threads: "+negativeThread.size());
		for( Set<GapOptimalModelingThread> clique:cliques)
		{
				// compute the total KL for each clique, using parameter recycling
				DependencyCombination t11=new DependencyCombination(clique, DiverseCBList, min_conversedKL_descr,negativeThread);
				if(negativeThread.size()>1000)
					executor.execute(t11);
				else
					t11.run();
				threadPool.add(t11);
		}
			
		//fetch threading result
		
			System.out.println("summarizing multiple threads result...");
			executor.shutdownAfterProcessingCurrentlyQueuedTasks();
			executor.awaitTerminationAfterShutdown();
			for(DependencyCombination t11:threadPool)
			{
				if(t11.bestscore>bestScore)
				{
					
					bestDgroups=t11.bestclique;
					recyclingEnhance=t11.maxDeltaKL;
					bestScore=t11.bestscore;
//					bestDivIndAction=t11.DivAction;
//					bestDivIndFinalAction=t11.DivIndAction;
//					BestDepAction=t11.DepAction;
//					bestReqNum=t11.bestRecycNum;
				}
			}
		
			int divIndColNum=DiverseCBList.size();
			HashSet<Integer> tPosSet=new HashSet<Integer>();
			for(GapOptimalModelingThread t : bestDgroups)
			{
				tPosSet.addAll(t.depend_Pos);
			}
			divIndColNum-=tPosSet.size();
			int totalDonateNumber=3*divIndColNum+min_conversedKL_descr.length;
			bestDivIndAction=new int[divIndColNum][divIndColNum*3+1];
			bestDivIndFinalAction=new int[totalDonateNumber];
			int maxDonateNum=totalDonateNumber +3*tPosSet.size();
			BestDepAction=new int[bestDgroups.size()][maxDonateNum];
			bestReqNum=cliqueScore2(bestDgroups, DiverseCBList, min_conversedKL_descr, bestDivIndAction, bestDivIndFinalAction, BestDepAction);
		///////////////////////////////decode action//////////////////////////
		if(bestDgroups.size()>0&&bestReqNum>0) 
		{
			int[] bt_depAction=backtracking_DP(BestDepAction, bestReqNum);
			
			int depId=0;
			for(GapOptimalModelingThread t : bestDgroups)
			{
				
				if(bt_depAction[depId]>0)
				{
					ConstrainBlock CB=t.PendingConstrainBlocks.get(bt_depAction[depId]);
					double oldscore=t.KL_Divergence;
					t.updateParaNum(t.dmerCount,new double[]{CB.lowerbound,CB.upperbound});
					if(oldscore<0)
					{
						double newscore=oldscore+CB.KL;
						System.out.println("negative Thread Improved:"+oldscore+" to "+newscore);
					}
				}
				depId++;
			}
			//change Independent column
			int totaldivcolReq=0;
			if(bestDivIndFinalAction.length>bestReqNum-3*tPosSet.size())
				totaldivcolReq=bestDivIndFinalAction[bestReqNum-3*tPosSet.size()];
			//conserved bases
			int[] bt_conAction=backtracking_DP(bestConAction,bestReqNum-totaldivcolReq-3*tPosSet.size());
			int cid=-1;
			for(Map.Entry<Integer, ArrayList<ConstrainBlock>>  elm:ConservedCBList.entrySet())
			{
				cid++;
				if(bt_conAction[cid]==0)
					continue;
				int orignalColumn=elm.getKey();
				HashMap<String,Double> dprobMap=new HashMap<String,Double>();
				ConstrainBlock  cb=elm.getValue().get(bt_conAction[cid]-1);
				double sumprob=0;
				int taken=0;
				HashSet<Integer> dPos=new HashSet<Integer>();
				dPos.add(orignalColumn);
				for (int i = 0; i < 4; i++) {
					double prob=m_matrix[orignalColumn][i];
					if(prob>(cb.upperbound+common.DoubleMinNormal)||prob<cb.lowerbound)
					{
						if(m_matrix[orignalColumn][i]==1)
							m_matrix[orignalColumn][i]=1-+common.DoubleMinNormal;
						sumprob+=m_matrix[orignalColumn][i];
						dprobMap.put(common.Hash2ACGT(i, 1), m_matrix[orignalColumn][i]);
						taken++;
					}
				}
				dprobMap.put("N",(1-sumprob)/(4-taken));
				System.out.println("recycling :"+dPos);
				Dmap.put(dPos, dprobMap);
				
			}
			//diverse bases
			cid=-1;
			int[] bt_divAction=backtracking_DP(bestDivIndAction, totaldivcolReq);
			for(Map.Entry<Integer, ArrayList<ConstrainBlock>>  elm:DiverseCBList.entrySet())
			{
				int transColumn=elm.getKey();
				if(tPosSet.contains(transColumn))
				{
					continue;
				}
				cid++;
	
				if(bt_divAction[cid]==0)
					continue;
				
				
				int orignalColumn=translate[transColumn];
				
				HashMap<String,Double> dprobMap=new HashMap<String,Double>();
				ConstrainBlock  cb=elm.getValue().get(bt_divAction[cid]-1);
				double sumprob=0;
				int taken=0;
				HashSet<Integer> dPos=new HashSet<Integer>();
				dPos.add(orignalColumn);
				for (int i = 0; i < 4; i++) {
					double prob=m_matrix[orignalColumn][i];
					if(prob>cb.upperbound+common.DoubleMinNormal||prob<cb.lowerbound)
					{
						if(m_matrix[orignalColumn][i]==1)
							m_matrix[orignalColumn][i]=1-common.DoubleMinNormal;
						sumprob+=m_matrix[orignalColumn][i];
						dprobMap.put(common.Hash2ACGT(i, 1), m_matrix[orignalColumn][i]);
						taken++;
					}
				}
				dprobMap.put("N",(1-sumprob)/(4-taken));
				System.out.println("recycling :"+dPos);
				Dmap.put(dPos, dprobMap);				
			}
			System.out.println("recycling enhance KL:"+recyclingEnhance);
		}
		
		
		
		double descSum=0;
		if(bestDgroups.size()>0)
			for(GapOptimalModelingThread t : bestDgroups)
			{
				System.out.println("max KL desc:"+t.toString());
				descSum+=t.KL_Divergence;
				Dmap.put(t.depend_Pos,t.DprobMap);
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
		int[][] action=new int[testarray.length][(testarray.length*(testarray[0].length))+1];
		double[] bounds=DP_computeMinKLdesc(testarray,action);
		System.out.println(Arrays.toString(bounds));
		
		
		
		double[] V1=new double[]{0,0.1,0.2,0.3};
		double[] V2=new double[]{0,0.01,0.22,0.3};
		int[] action2=new int[V1.length+V2.length-1];
		bounds=DP_mergeMaxKLincr(V1, V2,action2);
		System.out.println(Arrays.toString(bounds));
		
		
	}
	 
	
	
}
