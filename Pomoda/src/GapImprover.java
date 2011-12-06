import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map.Entry;
import java.util.Random;
import java.util.TreeMap;


import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.lang.ArrayUtils;
import org.biojava.bio.dist.DistributionTools;
import org.biojava.bio.dist.PairDistribution;
import org.biojava.bio.dist.UniformDistribution;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.jgrapht.graph.SimpleGraph;



import cern.colt.Arrays;
import EDU.oswego.cs.dl.util.concurrent.BoundedBuffer;
import EDU.oswego.cs.dl.util.concurrent.Executor;
import EDU.oswego.cs.dl.util.concurrent.LinkedQueue;
import EDU.oswego.cs.dl.util.concurrent.PooledExecutor;
import auc.AUCCalculator;
import auc.Confusion;


public class GapImprover {

	/**
	 * @param args
	 */
	public String outputPrefix="./";
	public String inputFasta;
	public String ctrlFasta="";
	public boolean OOPS=false; //only one dependence per sequence
	public boolean OOPG=false; //only one occurrence per sequence
	public boolean removeBG=false; //false:uniform BG assume
	public boolean PBMflag=false;
	public boolean is_bsitedata=false;
	public String bgmodelFile="";
	LinearEngine SearchEngine;
	public double sampling_ratio=1;
	public double FDR=0.01;
	public double entropyThresh=1;
	public int FlankLen=0;
	public int threadNum=4;
	public int max_gaplen=8;
	public double KLthresh=0.05;
	
	public BGModel background;
	public GapImprover(Pomoda motiffinder)
	{
		SearchEngine=motiffinder.SearchEngine2;
		sampling_ratio=motiffinder.sampling_ratio;
		FDR=motiffinder.FDR;
		background=motiffinder.background;
		
	}
	
	public GapImprover()
	{
		
	}
	
	public double KL_Divergence_empirical(List<String> sites,PWM motif,int start)
	{
		if(sites==null||sites.size()==0)
			return 10000;
		double KLsum=0;
		HashMap<String, Integer> sitecount=new HashMap<String, Integer>();
		Iterator<String> iter=sites.iterator();
		int totalCount=0;
		while(iter.hasNext())
		{
			String temp=iter.next();
			temp=temp.substring(start,start+motif.core_motiflen);
		
			if(temp.contains("N"))
				continue;
			if(sitecount.containsKey(temp))
			{
				sitecount.put(temp, sitecount.get(temp)+1);
			}
			else
			{
				sitecount.put(temp, 1);
			}
			totalCount++;
		}
		
		for(String key:sitecount.keySet())
		{
			int count=sitecount.get(key);
			double p=(double)count/totalCount;
			double motif_logP=0;
			if(key.length()==motif.core_motiflen&&key.length()>1)
				motif_logP=motif.scoreWeightMatrix(key);
			if(Double.isInfinite(motif_logP))
				continue;
		
			KLsum+=p*(Math.log(p)-motif_logP);
		}
			
			
		return KLsum;
		
	}

	public double KL_Divergence_empirical(List<String> sites,List<Double> SiteWeight,PWM motif,int start)
	{
		if(sites==null||sites.size()==0)
			return 10000;
		double KLsum=0;
		HashMap<String, Double> sitecount=new HashMap<String, Double>();
		Iterator<String> iter=sites.iterator();
		Iterator<Double> iter2=SiteWeight.iterator();
		double totalCount=0;
		while(iter.hasNext())
		{
			String temp=iter.next();
			double weight=iter2.next();
			temp=temp.substring(start,start+motif.core_motiflen);
		
			if(temp.contains("N"))
				continue;
			if(sitecount.containsKey(temp))
			{
				sitecount.put(temp, sitecount.get(temp)+weight);
			}
			else
			{
				sitecount.put(temp, weight);
			}
			totalCount+=weight;
		}
		
		for(String key:sitecount.keySet())
		{
			Double count=sitecount.get(key);
			double p=(double)count/totalCount;
			double motif_logP=0;
			if(key.length()==motif.core_motiflen&&key.length()>1)
				motif_logP=motif.scoreWeightMatrix(key);
			if(Double.isInfinite(motif_logP))
				continue;
		
			KLsum+=p*(Math.log(p)-motif_logP);
		}
			
			
		return KLsum;
		
	}

	public void initialize()
	{
		common.initialize();
		SearchEngine=new LinearEngine(4);
		if(this.PBMflag)
		{
			if(is_bsitedata)
				SearchEngine.buildPBM_index(inputFasta, 100000,false);
			else
				SearchEngine.buildPBM_index(inputFasta, 1000,true);
		}
		else
			SearchEngine.build_index(this.inputFasta);
	
		SearchEngine.num_thread=this.threadNum;
		background=new BGModel();
		File file=null;
		int bg_markov_order=3;

		
		if(!bgmodelFile.isEmpty())
		{
				background.LoadModel(bgmodelFile);
		}
		else if(!ctrlFasta.isEmpty())
		{

			     background.BuildModel(ctrlFasta, bg_markov_order+1); //3-order bg
			     background.SaveModel(ctrlFasta+".bg");
				
		}
		else
		{
			background=BGModel.CreateUniform();
		}

		
	}
	
	
	
	
	//allow multi dep-group in the same region, need to try different combinations
	public HashMap<HashSet<Integer>,HashMap<String,Double>> FindBest(List<GapBGModelingThread> list)
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
			
			if((baseScore-t2.KL_Divergence)>KLthresh&&t2.chisqPvalue<0.05)
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
		
		HashSet[] adjgraph=new HashSet[positiveThread.size()];
		for (int i = 0; i < positiveThread.size()-1; i++) {
			HashSet<Integer> s1=positiveThread.get(i).depend_Pos;
			for (int j = i+1; j < positiveThread.size(); j++) 
			{
				HashSet<Integer> s2=positiveThread.get(j).depend_Pos;
				boolean overlap=false;
				HashSet<Integer> intersect=new HashSet<Integer>(s2);
				intersect.retainAll(s1);
				if(intersect.size()>0)
					overlap=true;
				
				//flexible combination
				if(overlap==false)
				{
					if(adjgraph[i]==null)
						adjgraph[i]=new HashSet<Integer>();
					if(adjgraph[j]==null)
						adjgraph[j]=new HashSet<Integer>();
					adjgraph[i].add(j);
					adjgraph[j].add(i);
					HashSet<Integer> set=new HashSet<Integer>();
					set.add(i);
					set.add(j);
					queue.add(set);
				}
			}
		}
		//to here: queue contain all "2 dep-groups" set
		
		//enumerate all possible number of dep-group combinations
		 HashSet<Integer> topElm=null;
		 //breath first search, each child node is the flexible dep-group for all the parents.
		while((topElm=queue.poll())!=null)
		{
			double currscore=0;
			Iterator<Integer> iter=topElm.iterator();
			HashSet<Integer> commonThirdPoint=null;
			//all id (other dep-group) in commonThirdPoint can be use to form a combine with topElm(set of dep-group)
			while(iter.hasNext())
			{
				int id=iter.next();
				currscore+=positiveThread.get(id).KL_Divergence;
				if(commonThirdPoint==null)
				{
					
					commonThirdPoint=new HashSet<Integer>(adjgraph[id]);
				}
				else
				{
					commonThirdPoint.retainAll(adjgraph[id]);
				}
				
			}
			//push queue
			iter=commonThirdPoint.iterator();
			while(iter.hasNext())
			{
				HashSet<Integer> newCombine=new HashSet<Integer>(topElm);
				newCombine.add(iter.next());
				queue.add(newCombine);
			}
			if(currscore>bestScore)
			{
				bestDgroups=topElm;
				bestScore=currscore;
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
				
		if(descSum<baseScore*KLthresh)
		{
			System.out.println("No dependency found!");
			return new HashMap<HashSet<Integer>,HashMap<String,Double>>();
		}
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		
		return Dmap;
	}
	
	
	
	
	// enable parameters recycling
	public HashMap<HashSet<Integer>,HashMap<String,Double>> FindBest2(List<GapBGModelingThread> list,TreeMap<Double,KeyValuePair<Integer, String>> FreeParaQueue,double[][] m_matrix, int[] translate)
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
			
			if((baseScore-t2.KL_Divergence)>KLthresh&&t2.chisqPvalue<0.05)
			{
				//I reuse the field KL_Divergence as a score, not the KL_Divergence meaning any more
				t2.KL_Divergence=baseScore-t2.KL_Divergence;
				
				//parameter recycling
				double recyclingEnhance=0;
				Iterator<Double> penditer=t2.PendingParas.descendingKeySet().iterator();
				Iterator<Double> freeParaiter=FreeParaQueue.keySet().iterator();
				double pend=penditer.next();
				double free=10000;
				if(FreeParaQueue.keySet().size()>0)
					free=freeParaiter.next();
				while(freeParaiter.hasNext()&&t2.depend_Pos.contains(FreeParaQueue.get(free).key))
					free=freeParaiter.next();
				while(pend>free+10*common.DoubleMinNormal)
				{

					recyclingEnhance+=pend-free;
					  
					if(!penditer.hasNext()||!freeParaiter.hasNext()) 
						break;
					  free=freeParaiter.next();
					while(freeParaiter.hasNext()&&t2.depend_Pos.contains(FreeParaQueue.get(free).key))
							free=freeParaiter.next();
					
					  pend=penditer.next();
				}
				
				positiveThread.add(t2);
				//higher the better now, KL_Divergence is score here
				if(t2.KL_Divergence+recyclingEnhance>bestScore)
				{
					bestDgroups.clear();
					bestDgroups.add(positiveThread.size()-1);
					bestScore=t2.KL_Divergence+recyclingEnhance;
				}
				
			}
		}
		//here: bestDgroups is the best only-1 dep-group
		
		// build graph
		
		HashSet[] adjgraph=new HashSet[positiveThread.size()];
		for (int i = 0; i < positiveThread.size()-1; i++) {
			HashSet<Integer> s1=positiveThread.get(i).depend_Pos;
			for (int j = i+1; j < positiveThread.size(); j++) 
			{
				HashSet<Integer> s2=positiveThread.get(j).depend_Pos;
				boolean overlap=false;
				HashSet<Integer> intersect=new HashSet<Integer>(s2);
				intersect.retainAll(s1);
				if(intersect.size()>0)
					overlap=true;
				
				//flexible combination
				if(overlap==false)
				{
					if(adjgraph[i]==null)
						adjgraph[i]=new HashSet<Integer>();
					if(adjgraph[j]==null)
						adjgraph[j]=new HashSet<Integer>();
					adjgraph[i].add(j);
					adjgraph[j].add(i);
					HashSet<Integer> set=new HashSet<Integer>();
					set.add(i);
					set.add(j);
					queue.add(set);
				}
			}
		}
		//to here: queue contain all "2 dep-groups" set
		
		//enumerate all possible number of dep-group combinations
		 HashSet<Integer> topElm=null;
		 //breath first search, each child node is the flexible dep-group for all the parents.
		while((topElm=queue.poll())!=null)
		{
			double currscore=0;
			Iterator<Integer> iter=topElm.iterator();
			HashSet<Integer> commonThirdPoint=null;
			//all id (other dep-group) in commonThirdPoint can be use to form a combine with topElm(set of dep-group)
			while(iter.hasNext())
			{
				int id=iter.next();
				currscore+=positiveThread.get(id).KL_Divergence;
				if(commonThirdPoint==null)
				{
					
					commonThirdPoint=new HashSet<Integer>(adjgraph[id]);
				}
				else
				{
					commonThirdPoint.retainAll(adjgraph[id]);
				}
				
			}
			//push queue
			iter=commonThirdPoint.iterator();
			while(iter.hasNext())
			{
				HashSet<Integer> newCombine=new HashSet<Integer>(topElm);
				newCombine.add(iter.next());
				queue.add(newCombine);
			}
//////////////////////////parameters recycling////////////////////////
			double recyclingEnhance=0;
			TreeMap<Double,KeyValuePair<Integer, String>> PendingParas=new TreeMap<Double,KeyValuePair<Integer, String>>();
			HashSet<Integer> All_depend_Pos=new HashSet<Integer>();
			for(Integer tid:topElm)
			{
				All_depend_Pos.addAll(positiveThread.get(tid).depend_Pos);
				PendingParas.putAll(positiveThread.get(tid).PendingParas);
			}
			Iterator<Double> penditer=PendingParas.descendingKeySet().iterator();
			Iterator<Double> freeParaiter=FreeParaQueue.keySet().iterator();
			double pend=penditer.next();
			double free=10000;
			if(FreeParaQueue.keySet().size()>0)
			{
				free=freeParaiter.next();
				while(freeParaiter.hasNext()&&All_depend_Pos.contains(FreeParaQueue.get(free).key))
					free=freeParaiter.next();
			}
			while(pend>free+10*common.DoubleMinNormal)
			{
				recyclingEnhance+=pend-free;
				if(!penditer.hasNext()||!freeParaiter.hasNext())
					break;
				  pend=penditer.next();
				  free=freeParaiter.next();
				while(freeParaiter.hasNext()&&All_depend_Pos.contains(FreeParaQueue.get(free).key))
						free=freeParaiter.next();
			}
			
//////////////////////////parameters recycling////////////////////////			
			
			if(currscore+recyclingEnhance>bestScore)
			{
				bestDgroups=topElm;
				bestScore=currscore+recyclingEnhance;
			}
			
		}
		
		double descSum=0;
		if(bestDgroups.size()>0)
		{
//////////////////////////parameters recycling////////////////////////
			HashMap<Integer,HashSet<Integer>> delectedPara=new HashMap<Integer,HashSet<Integer>>();
			double recyclingEnhance=0;
			TreeMap<Double,KeyValuePair<Integer, String>> PendingParas=new TreeMap<Double,KeyValuePair<Integer, String>>();
			HashMap<Integer,Integer> hashcode2ThreadId=new HashMap<Integer, Integer>();
			HashSet<Integer> All_depend_Pos=new HashSet<Integer>();
			for(Integer tid:bestDgroups)
			{
				PendingParas.putAll(positiveThread.get(tid).PendingParas);
				All_depend_Pos.addAll(positiveThread.get(tid).depend_Pos);
				hashcode2ThreadId.put(positiveThread.get(tid).hashCode(), tid);
			}
			Iterator<Double> penditer=PendingParas.descendingKeySet().iterator();
			Iterator<Double> freeParaiter=FreeParaQueue.keySet().iterator();
			double pend=penditer.next();
			double free=10000;
			if(FreeParaQueue.keySet().size()>0)
			{
				free=freeParaiter.next();
				
				while(freeParaiter.hasNext()&&All_depend_Pos.contains(FreeParaQueue.get(free).key))
					free=freeParaiter.next();
			}
			int paraRecyclNum=0;
			while(pend>free+10*common.DoubleMinNormal)
			{
				String dmer=PendingParas.get(pend).value;
				int tid=hashcode2ThreadId.get(-PendingParas.get(pend).key);
				positiveThread.get(tid).DprobMap.put(dmer, positiveThread.get(tid).dmerCount[common.getHashing(dmer, 0, dmer.length())]);
				positiveThread.get(tid).KL_Divergence+=recyclingEnhance;
				
				KeyValuePair<Integer, String> fpair=FreeParaQueue.get(free);
				int colid=fpair.key;
				if(colid>=0)
					colid=translate[colid];
				else
					colid=0-colid-1;
				if(!delectedPara.containsKey(colid))
				{
					delectedPara.put(colid, new HashSet<Integer>());
				}
				delectedPara.get(colid).add(common.getHashing(fpair.getValue(),0,fpair.getValue().length()));
				paraRecyclNum++;
				recyclingEnhance+=pend-free;
				
				if(!penditer.hasNext()||!freeParaiter.hasNext())
					break;
				  pend=penditer.next();
				  free=freeParaiter.next();
				while(freeParaiter.hasNext()&&All_depend_Pos.contains(FreeParaQueue.get(free).key))
						free=freeParaiter.next();
			}
			//recompute "N" entry
			for(Integer id : bestDgroups)
			{
				double sum=0;
				int dmerlen=1;
				for(String Key:positiveThread.get( id ).DprobMap.keySet())
				{
					if(Key.equalsIgnoreCase("N"))
						continue;
					sum+=positiveThread.get( id ).DprobMap.get(Key);
					dmerlen=Key.length();
				}
				if(sum>1)
					sum=1;
				positiveThread.get(id ).DprobMap.put("N", (1-sum)/(Math.pow(4, dmerlen)-positiveThread.get(id ).DprobMap.size()+1));
			}
			
//			FreeParaQueue.clear();
//			FreeParaQueue.putAll(delectedPara);
			for(Integer col : delectedPara.keySet())
			{
				HashSet<Integer> dPos=new HashSet<Integer>();
				dPos.add(col);
				HashMap<String,Double> dprobMap=new HashMap<String,Double>();
				HashSet<Integer> deleteds=delectedPara.get(col);
				double minprob=1;
				int minsymid=0;
				double sumProb=0;
				for (int i = 0; i < 4; i++) {
					if(!deleteds.contains(i))
					{
						sumProb+=m_matrix[col][i];
						if(m_matrix[col][i]<minprob)
						{
							minprob=m_matrix[col][i];
							minsymid=i;
						}
					}
				}
				
				double Nprob=0.25;
				if(minprob<1)
				{
					sumProb-=minprob;
					Nprob=(1-sumProb)/(deleteds.size()+1);
					deleteds.add(minsymid);
				}
				dprobMap.put("N", Nprob);
				for (int i = 0; i < 4; i++) {
					if(!deleteds.contains(i))
					{
						dprobMap.put(common.Hash2ACGT(i, 1), m_matrix[col][i]);
					}
				}
				System.out.println("recycling :"+dPos);
				
				Dmap.put(dPos, dprobMap);
			}
			System.out.println("recycling enhance KL:"+recyclingEnhance);
			
//////////////////////////parameters recycling////////////////////////	
			for(Integer id : bestDgroups)
			{
				System.out.println("max KL desc:"+positiveThread.get(id).toString());

				descSum+=positiveThread.get(id).KL_Divergence;
				Dmap.put(positiveThread.get(id).depend_Pos,positiveThread.get(id).DprobMap);
			}
		}
				
//		if(bestScore<baseScore*KLthresh)
//		{
//			System.out.println("No dependency found!");
//			return new HashMap<HashSet<Integer>,HashMap<String,Double>>();
//		}
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		
		return Dmap;
	}
	
	public GapPWM refineGapPWM(GapPWM motif)
	{
		if(is_bsitedata)
			return motif;
		//get a set of instance strings, assume they are all real binding site
		LinkedList<String> sites=new LinkedList<String>();
		ArrayList<Double> siteWeight=null;
		//need to ignore the flanklen in the refine query sites
		int oldflank=FlankLen;
		int oldhead=motif.head;
		int oldtail=motif.tail;
		int oldlen=motif.core_motiflen;
		motif.core_motiflen=motif.columns();
		motif.head=0;
		motif.tail=0;
		FlankLen=0;
		KeyValuePair<LinkedList<String>, ArrayList<Double>> retpair=querySites(motif);
		FlankLen=oldflank;
		sites=retpair.key;
		siteWeight=retpair.value;
		//recover head tail, corelen
		motif.core_motiflen=oldlen;
		motif.head=oldhead;
		motif.tail=oldtail;
		
			if(sites==null||sites.size()<3)
				return motif;
			System.out.println("Number of refined Sites:"+sites.size());
			PWM dataPWM=null;
			
				try {
					
					if(SearchEngine.seqWeighting!=null)
					dataPWM=PWM.createPWM(sites.toArray(new String[1]), siteWeight.toArray(new Double[1]));
					else
						dataPWM=new PWM(sites.toArray(new String[1]));
				} catch (IllegalAlphabetException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				} catch (IllegalSymbolException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}

			
			//consider the whole motif length	
			HashMap<String, Double> sitecountMap=getSitemerFrequency(sites, siteWeight);
			GapBGModelingThread.KL_scorethresh=Double.MAX_VALUE;
			//compute the new probalitiy with sites data
			for(Integer dgroupId:motif.Dgroup_DmerProb.keySet())
			{
				HashSet<Integer> dpos=new HashSet<Integer>();
				for (int i = 0; i < motif.GroupId.length; i++) {
					if(motif.GroupId[i]==dgroupId)
						dpos.add(i);
				}
				if(dpos.size()==1)
					continue;
				GapBGModelingThread t1=null;
				if(removeBG)
					t1=new GapBGModelingThread(sitecountMap,dataPWM,dpos,background);//(gstart, gend, sites, dpos,background,siteWeight));//null mean not considering BG
				else
					 t1=new GapBGModelingThread(sitecountMap,dataPWM,dpos,null);
				
				t1.run();
				/////////////////parameters recycling/////////////////
				int num_top=motif.Dgroup_DmerProb.get(dgroupId).size()-1;
				TreeMap<Double,Integer> sorted_column=new TreeMap<Double,Integer>();
				for (int j = 0; j< t1.dmerCount.length; j++)
				{  
					double weight=t1.dmerCount[j]*(1-(j%num_top)*common.DoubleMinNormal);
					if(weight<0)
						weight=0;
					sorted_column.put(weight, j);
				}
				//sort the dmer by desc prob£¬ and take top 3d as the dependency model 
				int k=0;
				double sumprob=0;
				t1.DprobMap.clear();
				for(Double key:sorted_column.descendingKeySet())
				{
					k++;
						t1.DprobMap.put(common.Hash2ACGT(sorted_column.get(key), dpos.size()), key);
					sumprob+=key;

					if(t1.DprobMap.size()==(num_top))//
						break;
				}
				if(sumprob>1)
					sumprob=1;
				t1.DprobMap.put("N", (1-sumprob)/(t1.dmerCount.length-num_top));
				
                /////////////////parameters recycling/////////////////
				motif.Dgroup_DmerProb.put(dgroupId,t1.DprobMap);
			}
			
			//dont revise PWM¡¡content
//			for (int i = 0; i < dataPWM.columns(); i++) {
//				motif.setWeights(i, dataPWM.m_matrix[i]);
//			}
			
			
		
	
		
		
		
		return motif;
	}
	
	
	public void cleanupThread(LinkedList<GapBGModelingThread> threadsPool)
	{
		
	
	         Iterator<GapBGModelingThread> iter=threadsPool.iterator();
	         while(iter.hasNext())
	         {
	        	 GapBGModelingThread t1=iter.next();
	        		double thresh=GapBGModelingThread.KL_scorethresh-KLthresh*t1.depend_Pos.size();
	        	 	if(t1.KL_Divergence>thresh)
	        	 		iter.remove();
	         }
	}
	
	
	//this one dont make assumption conserved base will break the interaction regions
	public GapPWM fillDependency1(PWM motif)
	{
		GapPWM gapPWM=null;
		if(!is_bsitedata)
		motif.Consensus(true);

		
		//debug
	//	ArrayList<Double> snull = new ArrayList<Double>(),sbest =new ArrayList<Double>();
		
		//get a set of instance strings, assume they are all real binding site
		LinkedList<String> sites=new LinkedList<String>();
		ArrayList<Double> siteWeight=null;
		if(SearchEngine.seqWeighting!=null)
			siteWeight=new ArrayList<Double>(SearchEngine.ForwardStrand.size());
		if(!is_bsitedata)
		{
			KeyValuePair<LinkedList<String>, ArrayList<Double>> retpair=querySites(motif);
			sites=retpair.key;
			siteWeight=retpair.value;
		}
		else
		{
			//directly use the sequences as binding sites
			sites=SearchEngine.ForwardStrand;
			if(SearchEngine.seqWeighting!=null)
				siteWeight=SearchEngine.seqWeighting;
			
		}
		
		System.out.println("Number of Sites:"+sites.size());
		if(sites.size()<3)
			return GapPWM.createGapPWM(motif, new HashMap<HashSet<Integer>, HashMap<String,Double>>(),0);;
		PWM dataPWM=null;
		try {
			if(SearchEngine.seqWeighting!=null)
				dataPWM=PWM.createPWM(sites.toArray(new String[1]), siteWeight.toArray(new Double[1]));
			else
				dataPWM=new PWM(sites.toArray(new String[1]));
		} catch (IllegalAlphabetException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		} catch (IllegalSymbolException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		if(dataPWM==null)
			dataPWM=motif;
		if(motif==null||FlankLen>0)
			motif=dataPWM;
	
		
		//detect  conserved bases
		HashSet<Integer> conBases=new HashSet<Integer>();
		int start=-1;
		for (int i = motif.head; i < motif.head+motif.core_motiflen; i++) {
			double entropy=0;
			if(i<0||i>=motif.columns())
				entropy=2;
			else
				entropy=DistributionTools.totalEntropy(motif.getColumn(i-motif.head)) ;//use motif to determine conserved base
			if(entropy<entropyThresh)
			{
				
					start=i-motif.head;
					conBases.add(start);
			}
		}
		motif.Prior_EZ=conBases.size();
		System.out.println("Conserved Bases:"+conBases);
		if(conBases.size()>(motif.columns()-2))
			return GapPWM.createGapPWM(motif, new HashMap<HashSet<Integer>, HashMap<String,Double>>(),0);;
			
/******************************************************************************************/		
		//transate : remove the Conserved Columns and do again.
		
		LinkedList<String> sites2=new LinkedList<String>();
		for(String site:sites)
		{
			StringBuilder sb=new StringBuilder();
			for (int i = 0; i < site.length(); i++) {
				if(!conBases.contains(i))
				{
					sb.append(site.charAt(i));
				}
			}
			sites2.add(sb.toString());
		}
		PWM dataPWM2=null;
		HashSet<Integer> conBases2=new HashSet<Integer>(conBases);
		conBases.clear();
		try {
			if(SearchEngine.seqWeighting!=null)
				dataPWM2=PWM.createPWM(sites2.toArray(new String[1]), siteWeight.toArray(new Double[1]));
			else
				dataPWM2=new PWM(sites2.toArray(new String[1]));
			
			PWM temp=dataPWM;
			dataPWM=dataPWM2;
			dataPWM2=temp;
			LinkedList<String> tempsites=null;
			tempsites=sites;
			sites=sites2;
			sites2=tempsites;
			
		} catch (IllegalAlphabetException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		} catch (IllegalSymbolException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		int[] translateTB=new int[dataPWM.columns()];
		int Dj=0;
		for (int i = 0; i < translateTB.length; i++) {
			while(conBases2.contains(Dj))
			{
				Dj++;
			}
			translateTB[i]=Dj;
			Dj++;
		}
/******************************************************************************************/	
		
		
		try {
		//find the best dependency modeling in each gap region
		LinkedList<GapBGModelingThread> threadPool=new LinkedList<GapBGModelingThread>();
		HashMap<HashSet<Integer>,HashMap<String,Double>> Dmap=new HashMap<HashSet<Integer>,HashMap<String,Double>>();
		if(sites.size()>0)
		{
			
		//consider the whole motif length	
			HashMap<String, Double> sitecountMap=getSitemerFrequency(sites, siteWeight);
			int combinNum=(1<<Math.min(max_gaplen-1,dataPWM.columns()))*Math.max(1, dataPWM.columns()-max_gaplen+2); //+2 due to need a final to represent start point not in dep-group
		//Threads Pool
			PooledExecutor executor = new PooledExecutor(new LinkedQueue());
			executor.setMinimumPoolSize(threadNum);
			executor.setKeepAliveTime(-1);
			
			//add non-dep group thread
			HashSet<Integer> dpos=new HashSet<Integer>();
			GapBGModelingThread t1=null;
			if(removeBG)
				t1=new GapBGModelingThread(sitecountMap,dataPWM,dpos,background);//(gstart, gend, sites, dpos,background,siteWeight));//null mean not considering BG
			else
				 t1=new GapBGModelingThread(sitecountMap,dataPWM,dpos,null);//(gstart, gend, sites, dpos,null,siteWeight));//null mean not considering BG
			if(combinNum<=100||threadNum==1)
				t1.run();
			else
			{
			executor.execute(t1);
			}
			threadPool.add(t1);
			
			
		 {
			int gstart=0;
			int gend=dataPWM.columns();
			//each block assume the start base must in the dep-group,so only max_gaplen-1 binary length to consider
		
			
			int shift=0;//move the max_gaplen along the positions
			int blocksize=1<<(max_gaplen-1);
			for (int j = 0; j < combinNum; j++) {
				 dpos=new HashSet<Integer>();
				int bcode=j%blocksize;
				boolean convflag=false;
				shift=j/blocksize;
				if(shift!=(dataPWM.columns()-max_gaplen+1)&&dataPWM.columns()>=max_gaplen)
				{
					if(!conBases.contains(shift))
					{
						dpos.add(shift);
					}
					else 
						continue;
				}
				else
					shift--;//the final block
				
				int blocklen=Math.min(max_gaplen-1,dataPWM.columns());
				for (int j2 = 0; j2 < blocklen; j2++) {
					if(bcode%2==0)
					{
						dpos.add(j2+shift+1);
						if(conBases.contains(j2+shift+1))
						{
							convflag=true;
							break;
						}
					}
					bcode>>=1;
				}
				if(dpos.size()<2||convflag)
					continue;
				 t1=null;
				if(removeBG)
					t1=new GapBGModelingThread(sitecountMap,dataPWM,dpos,background);//(gstart, gend, sites, dpos,background,siteWeight);//null mean not considering BG
				else
					t1=new GapBGModelingThread(sitecountMap,dataPWM,dpos,null);//(gstart, gend, sites, dpos,null,siteWeight);//null mean not considering BG
				if(combinNum<=100||threadNum==1)
					t1.run();
				else
				{
				executor.execute(t1);
				}
	
				threadPool.add(t1);
//				if(threadPool.size()>10000&&threadPool.size()%10000==0)
//				{
//					Thread.sleep(100);
//					System.out.println("before:"+threadPool.size());
//					cleanupThread(threadPool);
//					System.out.println("after:"+threadPool.size());
//				}
			}

			
//			 Wait until all threads are finish
			if(combinNum>100&&threadNum!=1)
			{
			executor.shutdownAfterProcessingCurrentlyQueuedTasks();
			executor.awaitTerminationAfterShutdown();
			}
			if(!OOPG)
			{
				
			Dmap.putAll( FindBest(threadPool.subList(0, threadPool.size())));
			System.out.println("Final:"+threadPool.size());
			threadPool.clear();
			}
		}
		
		}
	
		Iterator<GapBGModelingThread> iter3=threadPool.iterator();
		 start=-1;
		double minKL=Double.MAX_VALUE;
		GapBGModelingThread bestThread=null;
		
		GapBGModelingThread[] Pos_BestThread=new GapBGModelingThread[motif.core_motiflen+2*FlankLen];
		
		if(OOPG)
		while(iter3.hasNext())
		{
			//find different flexible combinations in the same gap region(not allow two dependency group share the same position)
			GapBGModelingThread t1=iter3.next();	
				t1.join();
				if(t1.depend_Pos.size()==0)
				{
					//print PWM case, no dependancy 
					System.out.println(t1.toString());
					//snull=t1.debuglist;
				}
				if(t1.gapStart!=start)
				{
					start=t1.gapStart;
					minKL=t1.KL_Divergence;
					if(bestThread!=null)
					{
						System.out.println("best:"+bestThread.toString());
						if(bestThread.depend_Pos.size()>1)
						{
							Dmap.put(bestThread.depend_Pos, bestThread.DprobMap);
							
						}
					}
					bestThread=t1;
				}
				else 
				{
					if(t1.KL_Divergence<minKL)
					{
						minKL=t1.KL_Divergence;
						bestThread=t1;
					}
					
				}
				
				Iterator<Integer> iter4=t1.depend_Pos.iterator();
				while(iter4.hasNext())
				{
					int posId=iter4.next();
					if(Pos_BestThread[posId]!=null)
					{
						//as only one dep-group in the region, only check whether the given one is the best 
						if(Pos_BestThread[posId].KL_Divergence>t1.KL_Divergence)
							Pos_BestThread[posId]=t1;
					}
					else
					{
						Pos_BestThread[posId]=t1;
					}
				}
				if(t1.depend_Pos.size()==0)
				{
					for (int i = t1.gapStart; i < t1.gapEnd; i++) {
						if(Pos_BestThread[i].KL_Divergence>t1.KL_Divergence)
							Pos_BestThread[i]=t1;
						
					}
				}
		}
		if(bestThread!=null&&bestThread.depend_Pos.size()>1)
		{
			Dmap.put(bestThread.depend_Pos, bestThread.DprobMap);
		System.out.println("best:"+bestThread.toString());
		//sbest=bestThread.debuglist;
		}

		
//		if(!OOPG)
//		{
//		//fill in multi-dependency in the same gap region	
//		for (int i = 0; i < Pos_BestThread.length; i++) {
//			if(Pos_BestThread[i]!=null&&Pos_BestThread[i].depend_Pos.size()>1)
//			{
//				boolean maxCover=true;
//				for(Integer ii:Pos_BestThread[i].depend_Pos)
//				{
//					if(Pos_BestThread[ii].hashCode()!=Pos_BestThread[i].hashCode())
//					{
//						maxCover=false;
//						break;
//					}
//				}
//				if(maxCover)
//				{
//				Dmap.put(Pos_BestThread[i].depend_Pos, Pos_BestThread[i].DprobMap);
//				System.out.println(Pos_BestThread[i]);
//				}
//			}
//		}
//		}

/******************************************************************************************/	
		//translate back to orignial columns ids
		HashMap<HashSet<Integer>, HashMap<String, Double>> Dmap2=new HashMap<HashSet<Integer>, HashMap<String,Double>>();

		for(HashSet<Integer> keyset : Dmap.keySet())
		{
			HashSet<Integer> keyset2= new HashSet<Integer>();
			for(Integer dj : keyset)
			{
				keyset2.add(translateTB[dj]);
			}
			System.out.println("Translate:"+keyset.toString()+" To "+keyset2.toString());
			Dmap2.put(keyset2, Dmap.get(keyset));
			
		}
		Dmap=Dmap2;
		dataPWM=dataPWM2;
		sites=sites2;
/******************************************************************************************/			
		
		if(FlankLen==0)
			gapPWM=GapPWM.createGapPWM(motif.subPWM( motif.head,motif.head+motif.core_motiflen), Dmap,FlankLen);
		else
			gapPWM=GapPWM.createGapPWM(dataPWM.subPWM(FlankLen,dataPWM.columns()-FlankLen), Dmap,FlankLen);
//		if(gapPWM.core_motiflen!=motif.core_motiflen&&sites.size()>0)
//		{
//			motif=new PWM(sites.toArray(new String[1]));
//			gapPWM=GapPWM.createGapPWM(motif.subPWM(FlankLen,motif.columns()-FlankLen), Dmap,FlankLen);
//			motif=motif.subPWM(gapPWM.head, gapPWM.head+gapPWM.core_motiflen);
//		}
		dataPWM.head=gapPWM.head;
		dataPWM.core_motiflen=gapPWM.core_motiflen;
		dataPWM.tail=dataPWM.columns()-dataPWM.core_motiflen-dataPWM.head;
		
		if(siteWeight!=null)
		{
			System.out.println("orginal KL:"+KL_Divergence_empirical(sites,siteWeight, dataPWM,gapPWM.head) +"\timproved KL:"+KL_Divergence_empirical(sites,siteWeight, gapPWM,gapPWM.head));
		}
		else
		System.out.println("orginal KL:"+KL_Divergence_empirical(sites, dataPWM,gapPWM.head) +"\timproved KL:"+KL_Divergence_empirical(sites, gapPWM,gapPWM.head));
		
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} 
		
		return gapPWM;
	}
	
	
	
	
	//this one dont make assumption conserved base will break the interaction regions, optimal free parameter selection
	public GapPWM fillDependency1_1(PWM motif)
	{
		GapPWM gapPWM=null;
		if(!is_bsitedata)
		motif.Consensus(true);

		
		//debug
	//	ArrayList<Double> snull = new ArrayList<Double>(),sbest =new ArrayList<Double>();
		
		//get a set of instance strings, assume they are all real binding site
		LinkedList<String> sites=new LinkedList<String>();
		ArrayList<Double> siteWeight=null;
		if(SearchEngine.seqWeighting!=null)
			siteWeight=new ArrayList<Double>(SearchEngine.ForwardStrand.size());
		
		if(!is_bsitedata)
		{
			KeyValuePair<LinkedList<String>, ArrayList<Double>> retpair=querySites(motif);
			sites=retpair.key;
			siteWeight=retpair.value;
		}
		else
		{
			//directly use the sequences as binding sites
			sites=SearchEngine.ForwardStrand;
			if(SearchEngine.seqWeighting!=null)
				siteWeight=SearchEngine.seqWeighting;
			
		}
		
		System.out.println("Number of Sites:"+sites.size());
		if(sites.size()<3)
			return GapPWM.createGapPWM(motif, new HashMap<HashSet<Integer>, HashMap<String,Double>>(),0);;
		PWM dataPWM=null;
		try {
			if(SearchEngine.seqWeighting!=null)
				dataPWM=PWM.createPWM(sites.toArray(new String[1]), siteWeight.toArray(new Double[1]));
			else
				dataPWM=new PWM(sites.toArray(new String[1]));
		} catch (IllegalAlphabetException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		} catch (IllegalSymbolException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		if(dataPWM==null)
			dataPWM=motif;
		if(motif==null||FlankLen>0)
			motif=dataPWM;
	
		
		//detect  conserved bases
		HashSet<Integer> conBases=new HashSet<Integer>();
		int start=-1;
		for (int i = motif.head; i < motif.head+motif.core_motiflen; i++) {
			double entropy=0;
			if(i<0||i>=motif.columns())
				entropy=2;
			else
				entropy=DistributionTools.totalEntropy(motif.getColumn(i-motif.head)) ;//use motif to determine conserved base
			if(entropy<entropyThresh)
			{
				
					start=i-motif.head;
					conBases.add(start);
			}
		}
		motif.Prior_EZ=conBases.size();
		System.out.println("Conserved Bases:"+conBases);
		if(conBases.size()>(motif.columns()-2))
			return GapPWM.createGapPWM(motif, new HashMap<HashSet<Integer>, HashMap<String,Double>>(),0);;
			
/******************************************************************************************/		
		//transate : remove the Conserved Columns and do again.
		
		LinkedList<String> sites2=new LinkedList<String>();
		for(String site:sites)
		{
			StringBuilder sb=new StringBuilder();
			for (int i = 0; i < site.length(); i++) {
				if(!conBases.contains(i))
				{
					sb.append(site.charAt(i));
				}
			}
			sites2.add(sb.toString());
		}
		PWM dataPWM2=null;
		HashSet<Integer> conBases2=new HashSet<Integer>(conBases);
		conBases.clear();
		try {
			if(SearchEngine.seqWeighting!=null)
				dataPWM2=PWM.createPWM(sites2.toArray(new String[1]), siteWeight.toArray(new Double[1]));
			else
				dataPWM2=new PWM(sites2.toArray(new String[1]));
			
			PWM temp=dataPWM;
			dataPWM=dataPWM2;
			dataPWM2=temp;
			LinkedList<String> tempsites=null;
			tempsites=sites;
			sites=sites2;
			sites2=tempsites;
			
		} catch (IllegalAlphabetException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		} catch (IllegalSymbolException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		int[] translateTB=new int[dataPWM.columns()];
		int Dj=0;
		for (int i = 0; i < translateTB.length; i++) {
			while(conBases2.contains(Dj))
			{
				Dj++;
			}
			translateTB[i]=Dj;
			Dj++;
		}
/******************************************************************************************/	
		
		
		try {
		//find the best dependency modeling in each gap region
		LinkedList<GapBGModelingThread> threadPool=new LinkedList<GapBGModelingThread>();
		HashMap<HashSet<Integer>,HashMap<String,Double>> Dmap=new HashMap<HashSet<Integer>,HashMap<String,Double>>();
		if(sites.size()>0)
		{
			
		//consider the whole motif length	
			HashMap<String, Double> sitecountMap=getSitemerFrequency(sites, siteWeight);
			int combinNum=(1<<Math.min(max_gaplen-1,dataPWM.columns()))*Math.max(1, dataPWM.columns()-max_gaplen+2); //+2 due to need a final to represent start point not in dep-group
		//Threads Pool
			PooledExecutor executor = new PooledExecutor(new LinkedQueue());
			executor.setMinimumPoolSize(threadNum);
			executor.setKeepAliveTime(-1);
			
			//add non-dep group thread
			HashSet<Integer> dpos=new HashSet<Integer>();
			GapBGModelingThread t1=null;
			if(removeBG)
				t1=new GapOptimalModelingThread(sitecountMap,dataPWM,dpos,background);//(gstart, gend, sites, dpos,background,siteWeight));//null mean not considering BG
			else
				 t1=new GapOptimalModelingThread(sitecountMap,dataPWM,dpos,null);//(gstart, gend, sites, dpos,null,siteWeight));//null mean not considering BG
			if(combinNum<=100||threadNum==1)
				t1.run();
			else
			{
			executor.execute(t1);
			}
			threadPool.add(t1);
			
			
		 {
			int gstart=0;
			int gend=dataPWM.columns();
			//each block assume the start base must in the dep-group,so only max_gaplen-1 binary length to consider
		
			
			int shift=0;//move the max_gaplen along the positions
			int blocksize=1<<(max_gaplen-1);
			for (int j = 0; j < combinNum; j++) {
				 dpos=new HashSet<Integer>();
				int bcode=j%blocksize;
				boolean convflag=false;
				shift=j/blocksize;
				if(shift!=(dataPWM.columns()-max_gaplen+1)&&dataPWM.columns()>=max_gaplen)
				{
					if(!conBases.contains(shift))
					{
						dpos.add(shift);
					}
					else 
						continue;
				}
				else
					shift--;//the final block
				
				int blocklen=Math.min(max_gaplen-1,dataPWM.columns());
				for (int j2 = 0; j2 < blocklen; j2++) {
					if(bcode%2==0)
					{
						dpos.add(j2+shift+1);
						if(conBases.contains(j2+shift+1))
						{
							convflag=true;
							break;
						}
					}
					bcode>>=1;
				}
				if(dpos.size()<2||convflag)
					continue;
				 t1=null;
				if(removeBG)
					t1=new GapOptimalModelingThread(sitecountMap,dataPWM,dpos,background);//(gstart, gend, sites, dpos,background,siteWeight);//null mean not considering BG
				else
					t1=new GapOptimalModelingThread(sitecountMap,dataPWM,dpos,null);//(gstart, gend, sites, dpos,null,siteWeight);//null mean not considering BG
				if(combinNum<=100||threadNum==1)
					t1.run();
				else
				{
				executor.execute(t1);
				}
	
				threadPool.add(t1);
//				if(threadPool.size()>10000&&threadPool.size()%10000==0)
//				{
//					Thread.sleep(100);
//					System.out.println("before:"+threadPool.size());
//					cleanupThread(threadPool);
//					System.out.println("after:"+threadPool.size());
//				}
			}

			
//			 Wait until all threads are finish
			if(combinNum>100&&threadNum!=1)
			{
			executor.shutdownAfterProcessingCurrentlyQueuedTasks();
			executor.awaitTerminationAfterShutdown();
			}
			if(!OOPG)
			{
				
		//	Dmap.putAll( FindBest(threadPool.subList(0, threadPool.size())));
			Dmap.putAll( DependencyCombination.FindBest_1(threadPool.subList(0, threadPool.size())));
			System.out.println("Final:"+threadPool.size());
			threadPool.clear();
			}
		}
		
		}
	
		Iterator<GapBGModelingThread> iter3=threadPool.iterator();
		 start=-1;
		double minKL=Double.MAX_VALUE;
		GapBGModelingThread bestThread=null;
		
		GapBGModelingThread[] Pos_BestThread=new GapBGModelingThread[motif.core_motiflen+2*FlankLen];
		
		if(OOPG)
		while(iter3.hasNext())
		{
			//find different flexible combinations in the same gap region(not allow two dependency group share the same position)
			GapBGModelingThread t1=iter3.next();	
				t1.join();
				if(t1.depend_Pos.size()==0)
				{
					//print PWM case, no dependancy 
					System.out.println(t1.toString());
					//snull=t1.debuglist;
				}
				if(t1.gapStart!=start)
				{
					start=t1.gapStart;
					minKL=t1.KL_Divergence;
					if(bestThread!=null)
					{
						System.out.println("best:"+bestThread.toString());
						if(bestThread.depend_Pos.size()>1)
						{
							Dmap.put(bestThread.depend_Pos, bestThread.DprobMap);
							
						}
					}
					bestThread=t1;
				}
				else 
				{
					if(t1.KL_Divergence<minKL)
					{
						minKL=t1.KL_Divergence;
						bestThread=t1;
					}
					
				}
				
				Iterator<Integer> iter4=t1.depend_Pos.iterator();
				while(iter4.hasNext())
				{
					int posId=iter4.next();
					if(Pos_BestThread[posId]!=null)
					{
						//as only one dep-group in the region, only check whether the given one is the best 
						if(Pos_BestThread[posId].KL_Divergence>t1.KL_Divergence)
							Pos_BestThread[posId]=t1;
					}
					else
					{
						Pos_BestThread[posId]=t1;
					}
				}
				if(t1.depend_Pos.size()==0)
				{
					for (int i = t1.gapStart; i < t1.gapEnd; i++) {
						if(Pos_BestThread[i].KL_Divergence>t1.KL_Divergence)
							Pos_BestThread[i]=t1;
						
					}
				}
		}
		if(bestThread!=null&&bestThread.depend_Pos.size()>1)
		{
			Dmap.put(bestThread.depend_Pos, bestThread.DprobMap);
		System.out.println("best:"+bestThread.toString());
		//sbest=bestThread.debuglist;
		}

		
//		if(!OOPG)
//		{
//		//fill in multi-dependency in the same gap region	
//		for (int i = 0; i < Pos_BestThread.length; i++) {
//			if(Pos_BestThread[i]!=null&&Pos_BestThread[i].depend_Pos.size()>1)
//			{
//				boolean maxCover=true;
//				for(Integer ii:Pos_BestThread[i].depend_Pos)
//				{
//					if(Pos_BestThread[ii].hashCode()!=Pos_BestThread[i].hashCode())
//					{
//						maxCover=false;
//						break;
//					}
//				}
//				if(maxCover)
//				{
//				Dmap.put(Pos_BestThread[i].depend_Pos, Pos_BestThread[i].DprobMap);
//				System.out.println(Pos_BestThread[i]);
//				}
//			}
//		}
//		}

/******************************************************************************************/	
		//translate back to orignial columns ids
		HashMap<HashSet<Integer>, HashMap<String, Double>> Dmap2=new HashMap<HashSet<Integer>, HashMap<String,Double>>();

		for(HashSet<Integer> keyset : Dmap.keySet())
		{
			HashSet<Integer> keyset2= new HashSet<Integer>();
			for(Integer dj : keyset)
			{
				keyset2.add(translateTB[dj]);
			}
			System.out.println("Translate:"+keyset.toString()+" To "+keyset2.toString());
			Dmap2.put(keyset2, Dmap.get(keyset));
			
		}
		Dmap=Dmap2;
		dataPWM=dataPWM2;
		sites=sites2;
/******************************************************************************************/			
		
		if(FlankLen==0)
			gapPWM=GapPWM.createGapPWM(motif.subPWM( motif.head,motif.head+motif.core_motiflen), Dmap,FlankLen);
		else
			gapPWM=GapPWM.createGapPWM(dataPWM.subPWM(FlankLen,dataPWM.columns()-FlankLen), Dmap,FlankLen);
//		if(gapPWM.core_motiflen!=motif.core_motiflen&&sites.size()>0)
//		{
//			motif=new PWM(sites.toArray(new String[1]));
//			gapPWM=GapPWM.createGapPWM(motif.subPWM(FlankLen,motif.columns()-FlankLen), Dmap,FlankLen);
//			motif=motif.subPWM(gapPWM.head, gapPWM.head+gapPWM.core_motiflen);
//		}
		dataPWM.head=gapPWM.head;
		dataPWM.core_motiflen=gapPWM.core_motiflen;
		dataPWM.tail=dataPWM.columns()-dataPWM.core_motiflen-dataPWM.head;
		
		if(siteWeight!=null)
		{
			System.out.println("orginal KL:"+KL_Divergence_empirical(sites,siteWeight, dataPWM,gapPWM.head) +"\timproved KL:"+KL_Divergence_empirical(sites,siteWeight, gapPWM,gapPWM.head));
		}
		else
		System.out.println("orginal KL:"+KL_Divergence_empirical(sites, dataPWM,gapPWM.head) +"\timproved KL:"+KL_Divergence_empirical(sites, gapPWM,gapPWM.head));
		
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} 
		
		return gapPWM;
	}
	
	
	//enable the parameter recycling , optimal mode
	public GapPWM fillDependency2_1(PWM motif)
	{
		GapPWM gapPWM=null;
		if(!is_bsitedata)
		motif.Consensus(true);
		//the maximum number parameters of recycling is min(3n-3k,4^k-3k), n is the number of columns, k is the maximum number of dependency column in one group
		// when 3n-3k=4^k-3k => k=log_4(3n)
		GapOptimalModelingThread.ParaRecyclNum= (int)Math.ceil(3*motif.columns()- 3*Math.log( 3*motif.columns())/Math.log(4));
		
		//get a set of instance strings, assume they are all real binding site
		LinkedList<String> sites=null;
		ArrayList<Double> siteWeight=null;

		if(!is_bsitedata)
		{
			KeyValuePair<LinkedList<String>, ArrayList<Double>> retpair=querySites(motif);
			sites=retpair.key;
			siteWeight=retpair.value;
		}
		else
		{
			//directly use the sequences as binding sites
			sites=SearchEngine.ForwardStrand;
			if(SearchEngine.seqWeighting!=null)
				siteWeight=SearchEngine.seqWeighting;
			
		}
		
		System.out.println("Number of Sites:"+sites.size());
		if(sites.size()<3)
			return GapPWM.createGapPWM(motif, new HashMap<HashSet<Integer>, HashMap<String,Double>>(),0);;
		PWM dataPWM=null;
		try {
			if(SearchEngine.seqWeighting!=null)
				dataPWM=PWM.createPWM(sites.toArray(new String[1]), siteWeight.toArray(new Double[1]));
			else
				dataPWM=new PWM(sites.toArray(new String[1]));
		} catch (IllegalAlphabetException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		} catch (IllegalSymbolException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		if(dataPWM==null)
			dataPWM=motif;
		if(motif==null||FlankLen>0)
			motif=dataPWM;
	
		//sorted Free parameter queue
		TreeMap<Double,KeyValuePair<Integer, String>> FreeParaQueue=new TreeMap<Double,KeyValuePair<Integer, String>>();
		HashMap<Integer,ArrayList<ConstrainBlock>> ConservedCBList=new HashMap<Integer,ArrayList<ConstrainBlock>>();
		HashMap<Integer,ArrayList<ConstrainBlock>> DiverseCBList=new HashMap<Integer,ArrayList<ConstrainBlock>>();
		//detect  conserved bases
		HashSet<Integer> conBases=new HashSet<Integer>();
		int start=-1;
		for (int i = motif.head; i < motif.head+motif.core_motiflen; i++) {
			double entropy=0;
			if(i<0||i>=motif.columns())
				entropy=2;
			else
				entropy=DistributionTools.totalEntropy(motif.getColumn(i-motif.head)) ;//use motif to determine conserved base
			if(entropy<entropyThresh)
			{				
					start=i-motif.head;
					conBases.add(start);	
					//here untranslated column use negative colid, also recycle non-conserved bases later
					ConservedCBList.put(start, GapOptimalModelingThread.getConstrainReverseBlockQueueDeltaKL(motif.m_matrix[start], 3)) ;
					
			}
		}
		motif.Prior_EZ=conBases.size();
		System.out.println("Conserved Bases:"+conBases);
		if(conBases.size()>(motif.columns()-2))
			return GapPWM.createGapPWM(motif, new HashMap<HashSet<Integer>, HashMap<String,Double>>(),0);
			
/******************************************************************************************/		
		//transate : remove the Conserved Columns and do again.
		
		LinkedList<String> sites2=new LinkedList<String>();
		for(String site:sites)
		{
			StringBuilder sb=new StringBuilder();
			for (int i = 0; i < site.length(); i++) {
				if(!conBases.contains(i))
				{
					sb.append(site.charAt(i));
				}
			}
			sites2.add(sb.toString());
		}
		PWM dataPWM2=null;
		HashSet<Integer> conBases2=new HashSet<Integer>(conBases);
		conBases.clear();
		try {
			if(SearchEngine.seqWeighting!=null)
				dataPWM2=PWM.createPWM(sites2.toArray(new String[1]), siteWeight.toArray(new Double[1]));
			else
				dataPWM2=new PWM(sites2.toArray(new String[1]));
			
			PWM temp=dataPWM;
			dataPWM=dataPWM2;
			dataPWM2=temp;
			LinkedList<String> tempsites=null;
			tempsites=sites;
			sites=sites2;
			sites2=tempsites;
			
		} catch (IllegalAlphabetException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		} catch (IllegalSymbolException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		int[] translateTB=new int[dataPWM.columns()];
		int Dj=0;
		for (int i = 0; i < translateTB.length; i++) {
			while(conBases2.contains(Dj))
			{
				Dj++;
			}
			translateTB[i]=Dj;
			Dj++;
		}
		
		// recycle non-conserved bases using positive colid
		for (int i = 0; i <  translateTB.length; i++) {
			DiverseCBList.put(i, GapOptimalModelingThread.getConstrainReverseBlockQueueDeltaKL(motif.m_matrix[translateTB[i]], 3)) ;
		}

/******************************************************************************************/	

		
		try {
		//find the best dependency modeling in each gap region
		LinkedList<GapOptimalModelingThread> threadPool=new LinkedList<GapOptimalModelingThread>();
		HashMap<HashSet<Integer>,HashMap<String,Double>> Dmap=new HashMap<HashSet<Integer>,HashMap<String,Double>>();
		if(sites.size()>0)
		{
			
		//consider the whole motif length	
			HashMap<String, Double> sitecountMap=getSitemerFrequency(sites, siteWeight);
			int combinNum=(1<<Math.min(max_gaplen-1,dataPWM.columns()))*Math.max(1, dataPWM.columns()-max_gaplen+2); //+2 due to need a final to represent start point not in dep-group
		//Threads Pool
			PooledExecutor executor = new PooledExecutor(new LinkedQueue());
			executor.setMinimumPoolSize(threadNum);
			executor.setKeepAliveTime(-1);
			
			//add non-dep group thread
			HashSet<Integer> dpos=new HashSet<Integer>();
			GapOptimalModelingThread t1=null;
			if(removeBG)
				t1=new GapOptimalModelingThread(sitecountMap,dataPWM,dpos,background);//(gstart, gend, sites, dpos,background,siteWeight));//null mean not considering BG
			else
				 t1=new GapOptimalModelingThread(sitecountMap,dataPWM,dpos,null);//(gstart, gend, sites, dpos,null,siteWeight));//null mean not considering BG
			if(combinNum<=100||threadNum==1)
				t1.run();
			else
			{
				executor.execute(t1);
			}
			threadPool.add(t1);
			
			
		 {
			int gstart=0;
			int gend=dataPWM.columns();
			//each block assume the start base must in the dep-group,so only max_gaplen-1 binary length to consider
		
			
			int shift=0;//move the max_gaplen along the positions
			int blocksize=1<<(max_gaplen-1);
			for (int j = 0; j < combinNum; j++) {
				 dpos=new HashSet<Integer>();
				int bcode=j%blocksize;
				boolean convflag=false;
				shift=j/blocksize;
				if(shift!=(dataPWM.columns()-max_gaplen+1)&&dataPWM.columns()>=max_gaplen)
				{
					if(!conBases.contains(shift))
					{
						dpos.add(shift);
					}
					else 
						continue;
				}
				else
					shift--;//the final block
				
				int blocklen=Math.min(max_gaplen-1,dataPWM.columns());
				for (int j2 = 0; j2 < blocklen; j2++) {
					if(bcode%2==0)
					{
						dpos.add(j2+shift+1);
						if(conBases.contains(j2+shift+1))
						{
							convflag=true;
							break;
						}
					}
					bcode>>=1;
				}
				if(dpos.size()<2||convflag)
					continue;
				 t1=null;
				if(removeBG)
					t1=new GapOptimalModelingThread(sitecountMap,dataPWM,dpos,background);//(gstart, gend, sites, dpos,background,siteWeight);//null mean not considering BG
				else
					t1=new GapOptimalModelingThread(sitecountMap,dataPWM,dpos,null);//(gstart, gend, sites, dpos,null,siteWeight);//null mean not considering BG
				if(combinNum<=100||threadNum==1)
					t1.run();
				else
				{
				executor.execute(t1);
				}
	
				threadPool.add(t1);
			}

			
//			 Wait until all threads are finish
			if(combinNum>100&&threadNum!=1)
			{
			executor.shutdownAfterProcessingCurrentlyQueuedTasks();
			executor.awaitTerminationAfterShutdown();
			}
			if(!OOPG)
			{
				
			Dmap.putAll( DependencyCombination.FindBestCombination(threadPool.subList(0, threadPool.size()),ConservedCBList,DiverseCBList,motif.m_matrix,translateTB));
			System.out.println("Final:"+threadPool.size());
			threadPool.clear();
			}
		}
		
		}
	
		Iterator<GapOptimalModelingThread> iter3=threadPool.iterator();
		 start=-1;
		double minKL=Double.MAX_VALUE;
		GapOptimalModelingThread bestThread=null;
		
		GapOptimalModelingThread[] Pos_BestThread=new GapOptimalModelingThread[motif.core_motiflen+2*FlankLen];
		
		if(OOPG)
		while(iter3.hasNext())
		{
			//find different flexible combinations in the same gap region(not allow two dependency group share the same position)
			GapOptimalModelingThread t1=iter3.next();	
				t1.join();
				if(t1.depend_Pos.size()==0)
				{
					//print PWM case, no dependancy 
					System.out.println(t1.toString());
					//snull=t1.debuglist;
				}
				if(t1.gapStart!=start)
				{
					start=t1.gapStart;
					minKL=t1.KL_Divergence;
					if(bestThread!=null)
					{
						System.out.println("best:"+bestThread.toString());
						if(bestThread.depend_Pos.size()>1)
						{
							Dmap.put(bestThread.depend_Pos, bestThread.DprobMap);
							
						}
					}
					bestThread=t1;
				}
				else 
				{
					if(t1.KL_Divergence<minKL)
					{
						minKL=t1.KL_Divergence;
						bestThread=t1;
					}
					
				}
				
				Iterator<Integer> iter4=t1.depend_Pos.iterator();
				while(iter4.hasNext())
				{
					int posId=iter4.next();
					if(Pos_BestThread[posId]!=null)
					{
						//as only one dep-group in the region, only check whether the given one is the best 
						if(Pos_BestThread[posId].KL_Divergence>t1.KL_Divergence)
							Pos_BestThread[posId]=t1;
					}
					else
					{
						Pos_BestThread[posId]=t1;
					}
				}
				if(t1.depend_Pos.size()==0)
				{
					for (int i = t1.gapStart; i < t1.gapEnd; i++) {
						if(Pos_BestThread[i].KL_Divergence>t1.KL_Divergence)
							Pos_BestThread[i]=t1;
						
					}
				}
		}
		if(bestThread!=null&&bestThread.depend_Pos.size()>1)
		{
			Dmap.put(bestThread.depend_Pos, bestThread.DprobMap);
		System.out.println("best:"+bestThread.toString());
		//sbest=bestThread.debuglist;
		}

		

/******************************************************************************************/	
		//translate back to orignial columns ids
		HashMap<HashSet<Integer>, HashMap<String, Double>> Dmap2=new HashMap<HashSet<Integer>, HashMap<String,Double>>();

		for(HashSet<Integer> keyset : Dmap.keySet())
		{
			if(keyset.size()==1) //no need to translate for the recycled single based columns
			{
				Dmap2.put(keyset, Dmap.get(keyset));
				continue;
			}
			HashSet<Integer> keyset2= new HashSet<Integer>();
			
			for(Integer dj : keyset)
			{
				keyset2.add(translateTB[dj]);
			}
			System.out.println("Translate:"+keyset.toString()+" To "+keyset2.toString());
			Dmap2.put(keyset2, Dmap.get(keyset));
			
		}
		Dmap=Dmap2;
		dataPWM=dataPWM2;
		sites=sites2;
/******************************************************************************************/			
		
		//if(FlankLen==0)
			gapPWM=GapPWM.createGapPWM(motif.subPWM( motif.head,motif.head+motif.core_motiflen), Dmap,0);
		//else
			//gapPWM=GapPWM.createGapPWM(dataPWM.subPWM(FlankLen,dataPWM.columns()-FlankLen), Dmap,FlankLen);
//		if(gapPWM.core_motiflen!=motif.core_motiflen&&sites.size()>0)
//		{
//			motif=new PWM(sites.toArray(new String[1]));
//			gapPWM=GapPWM.createGapPWM(motif.subPWM(FlankLen,motif.columns()-FlankLen), Dmap,FlankLen);
//			motif=motif.subPWM(gapPWM.head, gapPWM.head+gapPWM.core_motiflen);
//		}
		dataPWM.head=gapPWM.head;
		dataPWM.core_motiflen=gapPWM.core_motiflen;
		dataPWM.tail=dataPWM.columns()-dataPWM.core_motiflen-dataPWM.head;
		
		if(siteWeight!=null)
		{
			System.out.println("orginal KL:"+KL_Divergence_empirical(sites,siteWeight, dataPWM,gapPWM.head) +"\timproved KL:"+KL_Divergence_empirical(sites,siteWeight, gapPWM,gapPWM.head));
		}
		else
		System.out.println("orginal KL:"+KL_Divergence_empirical(sites, dataPWM,gapPWM.head) +"\timproved KL:"+KL_Divergence_empirical(sites, gapPWM,gapPWM.head));
		
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} 
		
		return gapPWM;
	}
	
	public KeyValuePair<LinkedList<String>, ArrayList<Double>> querySites(PWM motif)
	{
		//get a set of instance strings, assume they are all real binding site
		LinkedList<String> sites=new LinkedList<String>();
		
		ArrayList<Double> siteWeight=null;
		if(SearchEngine.seqWeighting!=null)
			siteWeight=new ArrayList<Double>(SearchEngine.ForwardStrand.size());
		if(!OOPS)
		{
			//double pwmThresh=motif.getThresh(sampling_ratio, FDR, background,false);
			SearchEngine.EnableBackground(background);
			LinkedList<FastaLocation> falocs=SearchEngine.searchPattern(motif,2.99573227); //
			Iterator<FastaLocation> iter=falocs.iterator();
			 int lastseq=-1;
			 double highthresh=9.21034037;//;
        	 double seqcount=0;
        	 double maxseq_score=	Double.NEGATIVE_INFINITY;
        	 FastaLocation max_currloc=null;
        	 //take the best occurrences above 1.5 and all sites above 3.0
        	 if(falocs.size()<10)
        		 return null;
        	 while(iter.hasNext())
        	 {
        		 FastaLocation currloc=iter.next();
        		 if(lastseq!=currloc.getSeqId())
        		 {
        			 seqcount+=1;
        			
        			 if(lastseq!=-1&&maxseq_score<=highthresh)
        			 {
	     				String site=SearchEngine.getSite(max_currloc.getSeqId(), max_currloc.getSeqPos()-FlankLen, motif.core_motiflen+2*FlankLen);
	    				if(max_currloc.ReverseStrand)
	    					site=common.getReverseCompletementString(site);
	    				if(site!=null)
	    				{
	    					sites.add(site.toUpperCase());
	    					if(SearchEngine.seqWeighting!=null)
	    						siteWeight.add(SearchEngine.seqWeighting.get(max_currloc.getSeqId()));
	    				}
    				
    				//debug
    				//snull.add(max_currloc.Score+Math.log(0.25));
        			 }
        			 lastseq=currloc.getSeqId(); 
        			 maxseq_score=currloc.Score;
        			 max_currloc=currloc;
        		 }
        		 if(maxseq_score<currloc.Score)
        		 {
        			 maxseq_score=currloc.Score;
        			 max_currloc=currloc;
        		 }
        		 if(currloc.Score>highthresh)
        		 {
	     				String site=SearchEngine.getSite(currloc.getSeqId(), currloc.getSeqPos()-FlankLen, motif.core_motiflen+2*FlankLen);
	    				if(currloc.ReverseStrand)
	    					site=common.getReverseCompletementString(site);
	    				if(site!=null)
	    				{
	    					sites.add(site.toUpperCase());
	    					if(SearchEngine.seqWeighting!=null)
	    						siteWeight.add(SearchEngine.seqWeighting.get(currloc.getSeqId()));
	    				} 
        		 }
        	 }
        	 //last seq
 			String site=SearchEngine.getSite(max_currloc.getSeqId(), max_currloc.getSeqPos()-FlankLen, motif.core_motiflen+2*FlankLen);
			if(max_currloc.ReverseStrand)
				site=common.getReverseCompletementString(site);
			if(site!=null)
			{
			sites.add(site.toUpperCase());
			if(SearchEngine.seqWeighting!=null)
				siteWeight.add(SearchEngine.seqWeighting.get(max_currloc.getSeqId()));
			}
			//debug
			//snull.add(max_currloc.Score+Math.log(0.25));
			
			System.out.println("number of covered sequences: "+seqcount);
			
		}
		else
		{
			 LinkedList<FastaLocation> falocs =null;
			 if(removeBG)
			 {
				 SearchEngine.EnableBackground(background);
				 falocs=SearchEngine.searchPattern(motif, 0); //enable background in searching
			 }
			 else
				 falocs=SearchEngine.searchPattern(motif, Double.NEGATIVE_INFINITY);
        	 Iterator<FastaLocation> iter=falocs.iterator();
        	 int lastseq=-1;
        	 double seqcount=0;
        	 double maxseq_score=	Double.NEGATIVE_INFINITY;
        	 FastaLocation max_currloc=null;
        	 //take the best occurrences
        	 while(iter.hasNext())
        	 {
        		 FastaLocation currloc=iter.next();
        		 if(lastseq!=currloc.getSeqId())
        		 {
        			 seqcount+=1;
        			
        			 if(lastseq!=-1)
        			 {
     				String site=SearchEngine.getSite(max_currloc.getSeqId(), max_currloc.getSeqPos()-FlankLen, motif.core_motiflen+2*FlankLen);
    				if(max_currloc.ReverseStrand)
    					site=common.getReverseCompletementString(site);
    				if(site!=null)
    				{
    					sites.add(site.toUpperCase());
    					if(SearchEngine.seqWeighting!=null)
    						siteWeight.add(SearchEngine.seqWeighting.get(max_currloc.getSeqId()));
    				}
    				
    				//debug
    				//snull.add(max_currloc.Score+Math.log(0.25));
        			 }
        			 lastseq=currloc.getSeqId(); 
        			 maxseq_score=currloc.Score;
        			 max_currloc=currloc;
        		 }
        		 if(maxseq_score<currloc.Score)
        		 {
        			 maxseq_score=currloc.Score;
        			 max_currloc=currloc;
        		 }
        	 }
        	 //last seq
 			String site=SearchEngine.getSite(max_currloc.getSeqId(), max_currloc.getSeqPos()-FlankLen, motif.core_motiflen+2*FlankLen);
			if(max_currloc.ReverseStrand)
				site=common.getReverseCompletementString(site);
			if(site!=null)
			{
			sites.add(site.toUpperCase());
			if(SearchEngine.seqWeighting!=null)
				siteWeight.add(SearchEngine.seqWeighting.get(max_currloc.getSeqId()));
			}
			//debug
			//snull.add(max_currloc.Score+Math.log(0.25));
		}
		
		KeyValuePair<LinkedList<String>, ArrayList<Double>> pair=new KeyValuePair<LinkedList<String>, ArrayList<Double>>(sites, siteWeight);
		return pair;
	}
	
	static TreeMap<Double,KeyValuePair<Integer, String>> getFreeParaQueue(double[] probArr, int colid)
	{
		TreeMap<Double,KeyValuePair<Integer, String>> queue=new TreeMap<Double,KeyValuePair<Integer, String>>();
		int kmerlen=(int)Math.round(Math.log(probArr.length)/Math.log(4));
		/* E_k=l_klog(l_k/(|P|-k))+\sum_{i=1}^{k}{P_ilogP_i}  # entropy
		 * l_k=\sum_{i=k+1}^{|P|}{P_i}
		 * delta_k=E_{k+1}-E_k=(l_k-P_k)log((l_k-P_k)/(|P|-k-1))-l_klog(l_k/(|P|-k))+P_klog(P_k)
		 * 
		 * 
		 */
		int num_top=100;
		TreeMap<Double,Integer> sorted_column=new TreeMap<Double,Integer>();
		for (int j = 0; j< probArr.length; j++)
		{  
			double weight=probArr[j]*10000*(1-(j%probArr.length)*common.DoubleMinNormal);
	
			sorted_column.put(weight, j); //weight just for ordering 
		}
	
	double sumprob=0;
	
	double l_k=0;
	int k=probArr.length;
	for(Double key:sorted_column.keySet())
	{
		k--;
		int symid=sorted_column.get(key);
		String kmer=common.Hash2ACGT(symid, kmerlen);
		l_k+=probArr[symid];
		double P_k=probArr[symid];
		double delta_k=-l_k*Math.log(l_k/(probArr.length-k))+P_k*Math.log(P_k);
		if(k<probArr.length-1) // the first one always 0, as least vaule is the same as the rest average of removing one 
		{
			delta_k+=(l_k-P_k)*Math.log((l_k-P_k)/(probArr.length-k-1));
		queue.put(delta_k+((k+colid)%num_top)*common.DoubleMinNormal, new KeyValuePair<Integer, String>(colid,kmer));
		}
	}
		
		
		return queue;
	}
	
	
	//enable the parameter recycling
	public GapPWM fillDependency2(PWM motif)
	{
		GapPWM gapPWM=null;
		if(!is_bsitedata)
		motif.Consensus(true);

		//get a set of instance strings, assume they are all real binding site
		LinkedList<String> sites=null;
		ArrayList<Double> siteWeight=null;

		if(!is_bsitedata)
		{
			KeyValuePair<LinkedList<String>, ArrayList<Double>> retpair=querySites(motif);
			sites=retpair.key;
			siteWeight=retpair.value;
		}
		else
		{
			//directly use the sequences as binding sites
			sites=SearchEngine.ForwardStrand;
			if(SearchEngine.seqWeighting!=null)
				siteWeight=SearchEngine.seqWeighting;
			
		}
		
		System.out.println("Number of Sites:"+sites.size());
		if(sites.size()<3)
			return GapPWM.createGapPWM(motif, new HashMap<HashSet<Integer>, HashMap<String,Double>>(),0);;
		PWM dataPWM=null;
		try {
			if(SearchEngine.seqWeighting!=null)
				dataPWM=PWM.createPWM(sites.toArray(new String[1]), siteWeight.toArray(new Double[1]));
			else
				dataPWM=new PWM(sites.toArray(new String[1]));
		} catch (IllegalAlphabetException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		} catch (IllegalSymbolException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		if(dataPWM==null)
			dataPWM=motif;
		if(motif==null||FlankLen>0)
			motif=dataPWM;
	
		//sorted Free parameter queue
		TreeMap<Double,KeyValuePair<Integer, String>> FreeParaQueue=new TreeMap<Double,KeyValuePair<Integer, String>>();
		//detect  conserved bases
		HashSet<Integer> conBases=new HashSet<Integer>();
		int start=-1;
		for (int i = motif.head; i < motif.head+motif.core_motiflen; i++) {
			double entropy=0;
			if(i<0||i>=motif.columns())
				entropy=2;
			else
				entropy=DistributionTools.totalEntropy(motif.getColumn(i-motif.head)) ;//use motif to determine conserved base
			if(entropy<entropyThresh)
			{				
					start=i-motif.head;
					conBases.add(start);	
					
					//here untranslated column use negative colid, also recycle non-conserved bases later
					TreeMap<Double,KeyValuePair<Integer, String>> paraslist=getFreeParaQueue( motif.m_matrix[start],-start-1);
					for(Double key : paraslist.keySet())
					{
						FreeParaQueue.put(key, paraslist.get(key));
					}
					
			}

		}
		motif.Prior_EZ=conBases.size();
		System.out.println("Conserved Bases:"+conBases);
		if(conBases.size()>(motif.columns()-2))
			return GapPWM.createGapPWM(motif, new HashMap<HashSet<Integer>, HashMap<String,Double>>(),0);
			
/******************************************************************************************/		
		//transate : remove the Conserved Columns and do again.
		
		LinkedList<String> sites2=new LinkedList<String>();
		for(String site:sites)
		{
			StringBuilder sb=new StringBuilder();
			for (int i = 0; i < site.length(); i++) {
				if(!conBases.contains(i))
				{
					sb.append(site.charAt(i));
				}
			}
			sites2.add(sb.toString());
		}
		PWM dataPWM2=null;
		HashSet<Integer> conBases2=new HashSet<Integer>(conBases);
		conBases.clear();
		try {
			if(SearchEngine.seqWeighting!=null)
				dataPWM2=PWM.createPWM(sites2.toArray(new String[1]), siteWeight.toArray(new Double[1]));
			else
				dataPWM2=new PWM(sites2.toArray(new String[1]));
			
			PWM temp=dataPWM;
			dataPWM=dataPWM2;
			dataPWM2=temp;
			LinkedList<String> tempsites=null;
			tempsites=sites;
			sites=sites2;
			sites2=tempsites;
			
		} catch (IllegalAlphabetException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		} catch (IllegalSymbolException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		int[] translateTB=new int[dataPWM.columns()];
		int Dj=0;
		for (int i = 0; i < translateTB.length; i++) {
			while(conBases2.contains(Dj))
			{
				Dj++;
			}
			translateTB[i]=Dj;
			Dj++;
		}
		
		// recycle non-conserved bases using positive colid
		for (int i = 0; i <  translateTB.length; i++) {
			TreeMap<Double,KeyValuePair<Integer, String>> paraslist=getFreeParaQueue( motif.m_matrix[translateTB[i]],i);
			for(Double key : paraslist.keySet())
			{
				FreeParaQueue.put(key, paraslist.get(key));
			}
		}

/******************************************************************************************/	
		
		
		try {
		//find the best dependency modeling in each gap region
		LinkedList<GapBGModelingThread> threadPool=new LinkedList<GapBGModelingThread>();
		HashMap<HashSet<Integer>,HashMap<String,Double>> Dmap=new HashMap<HashSet<Integer>,HashMap<String,Double>>();
		if(sites.size()>0)
		{
			
		//consider the whole motif length	
			HashMap<String, Double> sitecountMap=getSitemerFrequency(sites, siteWeight);
			int combinNum=(1<<Math.min(max_gaplen-1,dataPWM.columns()))*Math.max(1, dataPWM.columns()-max_gaplen+2); //+2 due to need a final to represent start point not in dep-group
		//Threads Pool
			PooledExecutor executor = new PooledExecutor(new LinkedQueue());
			executor.setMinimumPoolSize(threadNum);
			executor.setKeepAliveTime(-1);
			
			//add non-dep group thread
			HashSet<Integer> dpos=new HashSet<Integer>();
			GapBGModelingThread t1=null;
			if(removeBG)
				t1=new GapBGModelingThread(sitecountMap,dataPWM,dpos,background);//(gstart, gend, sites, dpos,background,siteWeight));//null mean not considering BG
			else
				 t1=new GapBGModelingThread(sitecountMap,dataPWM,dpos,null);//(gstart, gend, sites, dpos,null,siteWeight));//null mean not considering BG
			if(combinNum<=100||threadNum==1)
				t1.run();
			else
			{
			executor.execute(t1);
			}
			threadPool.add(t1);
			
			
		 {
			int gstart=0;
			int gend=dataPWM.columns();
			//each block assume the start base must in the dep-group,so only max_gaplen-1 binary length to consider
		
			
			int shift=0;//move the max_gaplen along the positions
			int blocksize=1<<(max_gaplen-1);
			for (int j = 0; j < combinNum; j++) {
				 dpos=new HashSet<Integer>();
				int bcode=j%blocksize;
				boolean convflag=false;
				shift=j/blocksize;
				if(shift!=(dataPWM.columns()-max_gaplen+1)&&dataPWM.columns()>=max_gaplen)
				{
					if(!conBases.contains(shift))
					{
						dpos.add(shift);
					}
					else 
						continue;
				}
				else
					shift--;//the final block
				
				int blocklen=Math.min(max_gaplen-1,dataPWM.columns());
				for (int j2 = 0; j2 < blocklen; j2++) {
					if(bcode%2==0)
					{
						dpos.add(j2+shift+1);
						if(conBases.contains(j2+shift+1))
						{
							convflag=true;
							break;
						}
					}
					bcode>>=1;
				}
				if(dpos.size()<2||convflag)
					continue;
				 t1=null;
				if(removeBG)
					t1=new GapBGModelingThread(sitecountMap,dataPWM,dpos,background);//(gstart, gend, sites, dpos,background,siteWeight);//null mean not considering BG
				else
					t1=new GapBGModelingThread(sitecountMap,dataPWM,dpos,null);//(gstart, gend, sites, dpos,null,siteWeight);//null mean not considering BG
				if(combinNum<=100||threadNum==1)
					t1.run();
				else
				{
				executor.execute(t1);
				}
	
				threadPool.add(t1);
			}

			
//			 Wait until all threads are finish
			if(combinNum>100&&threadNum!=1)
			{
			executor.shutdownAfterProcessingCurrentlyQueuedTasks();
			executor.awaitTerminationAfterShutdown();
			}
			if(!OOPG)
			{
				
			Dmap.putAll( FindBest2(threadPool.subList(0, threadPool.size()),FreeParaQueue,motif.m_matrix,translateTB));
			System.out.println("Final:"+threadPool.size());
			threadPool.clear();
			}
		}
		
		}
	
		Iterator<GapBGModelingThread> iter3=threadPool.iterator();
		 start=-1;
		double minKL=Double.MAX_VALUE;
		GapBGModelingThread bestThread=null;
		
		GapBGModelingThread[] Pos_BestThread=new GapBGModelingThread[motif.core_motiflen+2*FlankLen];
		
		if(OOPG)
		while(iter3.hasNext())
		{
			//find different flexible combinations in the same gap region(not allow two dependency group share the same position)
			GapBGModelingThread t1=iter3.next();	
				t1.join();
				if(t1.depend_Pos.size()==0)
				{
					//print PWM case, no dependancy 
					System.out.println(t1.toString());
					//snull=t1.debuglist;
				}
				if(t1.gapStart!=start)
				{
					start=t1.gapStart;
					minKL=t1.KL_Divergence;
					if(bestThread!=null)
					{
						System.out.println("best:"+bestThread.toString());
						if(bestThread.depend_Pos.size()>1)
						{
							Dmap.put(bestThread.depend_Pos, bestThread.DprobMap);
							
						}
					}
					bestThread=t1;
				}
				else 
				{
					if(t1.KL_Divergence<minKL)
					{
						minKL=t1.KL_Divergence;
						bestThread=t1;
					}
					
				}
				
				Iterator<Integer> iter4=t1.depend_Pos.iterator();
				while(iter4.hasNext())
				{
					int posId=iter4.next();
					if(Pos_BestThread[posId]!=null)
					{
						//as only one dep-group in the region, only check whether the given one is the best 
						if(Pos_BestThread[posId].KL_Divergence>t1.KL_Divergence)
							Pos_BestThread[posId]=t1;
					}
					else
					{
						Pos_BestThread[posId]=t1;
					}
				}
				if(t1.depend_Pos.size()==0)
				{
					for (int i = t1.gapStart; i < t1.gapEnd; i++) {
						if(Pos_BestThread[i].KL_Divergence>t1.KL_Divergence)
							Pos_BestThread[i]=t1;
						
					}
				}
		}
		if(bestThread!=null&&bestThread.depend_Pos.size()>1)
		{
			Dmap.put(bestThread.depend_Pos, bestThread.DprobMap);
		System.out.println("best:"+bestThread.toString());
		//sbest=bestThread.debuglist;
		}

		

/******************************************************************************************/	
		//translate back to orignial columns ids
		HashMap<HashSet<Integer>, HashMap<String, Double>> Dmap2=new HashMap<HashSet<Integer>, HashMap<String,Double>>();

		for(HashSet<Integer> keyset : Dmap.keySet())
		{
			if(keyset.size()==1) //no need to translate for the recycled single based columns
			{
				Dmap2.put(keyset, Dmap.get(keyset));
				continue;
			}
			HashSet<Integer> keyset2= new HashSet<Integer>();
			
			for(Integer dj : keyset)
			{
				keyset2.add(translateTB[dj]);
			}
			System.out.println("Translate:"+keyset.toString()+" To "+keyset2.toString());
			Dmap2.put(keyset2, Dmap.get(keyset));
			
		}
		Dmap=Dmap2;
		dataPWM=dataPWM2;
		sites=sites2;
/******************************************************************************************/			
		
		//if(FlankLen==0)
			gapPWM=GapPWM.createGapPWM(motif.subPWM( motif.head,motif.head+motif.core_motiflen), Dmap,0);
		//else
			//gapPWM=GapPWM.createGapPWM(dataPWM.subPWM(FlankLen,dataPWM.columns()-FlankLen), Dmap,FlankLen);
//		if(gapPWM.core_motiflen!=motif.core_motiflen&&sites.size()>0)
//		{
//			motif=new PWM(sites.toArray(new String[1]));
//			gapPWM=GapPWM.createGapPWM(motif.subPWM(FlankLen,motif.columns()-FlankLen), Dmap,FlankLen);
//			motif=motif.subPWM(gapPWM.head, gapPWM.head+gapPWM.core_motiflen);
//		}
		dataPWM.head=gapPWM.head;
		dataPWM.core_motiflen=gapPWM.core_motiflen;
		dataPWM.tail=dataPWM.columns()-dataPWM.core_motiflen-dataPWM.head;
		
		if(siteWeight!=null)
		{
			System.out.println("orginal KL:"+KL_Divergence_empirical(sites,siteWeight, dataPWM,gapPWM.head) +"\timproved KL:"+KL_Divergence_empirical(sites,siteWeight, gapPWM,gapPWM.head));
		}
		else
		System.out.println("orginal KL:"+KL_Divergence_empirical(sites, dataPWM,gapPWM.head) +"\timproved KL:"+KL_Divergence_empirical(sites, gapPWM,gapPWM.head));
		
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} 
		
		return gapPWM;
	}
	
	public double CorrelationTest(PWM motif)
	{
		TreeMap<Double,Integer> Sorted_labels=new TreeMap<Double,Integer>();
        
   	 LinkedList<FastaLocation> falocs =SearchEngine.searchPattern(motif, Double.NEGATIVE_INFINITY);
   	 Iterator<FastaLocation> iter=falocs.iterator();
   	 int lastseq=-1;
   	 double seqcount=0;
   	 double maxseq_score=	Double.NEGATIVE_INFINITY;
   	 double[] scores=new double[SearchEngine.seqWeighting.size()];
   	 int seqid=0;
   	 while(iter.hasNext())
   	 {
   		 FastaLocation currloc=iter.next();
   		 if(lastseq!=currloc.getSeqId())
   		 {
   			 seqcount+=1;
   			 if(lastseq!=-1)
   			 {
   				 Sorted_labels.put(maxseq_score+seqcount*common.DoubleMinNormal, 1);
   				scores[seqid]=Math.exp(maxseq_score);//Math.exp
   				seqid++;
   			 }
   				 lastseq=currloc.getSeqId();
   			 maxseq_score=currloc.Score;
   		 }
   		 if(maxseq_score<currloc.Score)
   		 {
   			 maxseq_score=currloc.Score;
   		 }
   	 }
   	Sorted_labels.put(maxseq_score+seqcount*common.DoubleMinNormal, 1);
		scores[seqid]=Math.exp(maxseq_score); //Math.exp
			
		double corr=common.getPearsonCorrelation(scores,ArrayUtils.toPrimitive(SearchEngine.seqWeighting.toArray(new Double[1])));
   	
    return corr;
	}
	
	
	static public HashMap<String,Double> getSitemerFrequency(List<String> Sites,ArrayList<Double> seqWeighting)
	{
		String[] gapstr=new String[ Sites.size()];
		int gapStart=0;
		int gapEnd=Sites.get(0).length();
		Iterator<String> iter=Sites.iterator();
		int i=0;
		while(iter.hasNext())
		{
			String temp=iter.next();
			gapstr[i]=temp.substring(gapStart, gapEnd);
		
			 i++;
		}
		
	
			//build Kmer model
			HashMap<String,Double> gapmerCount=new HashMap<String,Double>(gapstr.length);

			double sumWeight=0;
			for (int j = 0; j < gapstr.length; j++) {
				String dmer="";
				double weight=1;
				if(seqWeighting!=null)
					weight=seqWeighting.get(j);
				
				//Iterator<Integer> iter2=depend_Pos.iterator();
				if(gapstr[j].contains("N"))
					continue;
				if(gapmerCount.containsKey(gapstr[j]))
					gapmerCount.put(gapstr[j], gapmerCount.get(gapstr[j])+weight) ;
				else
					gapmerCount.put(gapstr[j],weight);
				sumWeight+=weight;

			}
			//Normalize gapmerCount
			for (String key:gapmerCount.keySet()) {
				gapmerCount.put(key, gapmerCount.get(key)/sumWeight);
			}
			
			return gapmerCount;
	}
	
	public double AUCtest(PWM motif)
	{
      
		
         LinearEngine BGSearch=new LinearEngine(6);
         BGSearch.accSeqLen.add(0);
         Iterator<String> iter2=SearchEngine.ForwardStrand.iterator();
         background.r.setSeed(0);
         while(iter2.hasNext())
         {
        	 int len=iter2.next().length();
        	 String bgstr="";
        	 if(removeBG)
        	 {
	        	 KeyValuePair<Double, String> bgstr_p=background.generateRandomSequence(len);
	        	 bgstr=bgstr_p.value;
        	 }
        	 else
        	 {
	        	 UniformDistribution ud=new UniformDistribution(DNATools.getDNA());
	        	 bgstr=DistributionTools.generateSymbolList(ud, len).seqString();
        	 }
        	 BGSearch.ForwardStrand.add(bgstr);	 
        	 BGSearch.TotalLen+=bgstr.length();
        	 BGSearch.accSeqLen.add(BGSearch.TotalLen);
         }

     	TreeMap<Double,Integer> Sorted_labels=new TreeMap<Double,Integer>();
         
        	 LinkedList<FastaLocation> falocs =SearchEngine.searchPattern(motif, Double.NEGATIVE_INFINITY);
        	 Iterator<FastaLocation> iter=falocs.iterator();
        	 int lastseq=-1;
        	 double seqcount=0;
        	 double maxseq_score=	Double.NEGATIVE_INFINITY;
        	 while(iter.hasNext())
        	 {
        		 FastaLocation currloc=iter.next();
        		 if(lastseq!=currloc.getSeqId())
        		 {
        			 seqcount+=1;
        			 if(lastseq!=-1)
        			 {
        				 Sorted_labels.put(maxseq_score+seqcount*common.DoubleMinNormal, 1);
        			 }
        				 lastseq=currloc.getSeqId();
        			 maxseq_score=currloc.Score;
        		 }
        		 if(maxseq_score<currloc.Score)
        		 {
        			 maxseq_score=currloc.Score;

        		 }
        	 }
        	 Sorted_labels.put(maxseq_score+seqcount*common.DoubleMinNormal, 1);
        	 
        	 //bg sequences
        	 falocs =BGSearch.searchPattern(motif, Double.NEGATIVE_INFINITY);
        	 iter=falocs.iterator();
        	lastseq=-1;
        	 seqcount=0;
        	maxseq_score=Double.NEGATIVE_INFINITY;
	       	 while(iter.hasNext())
	    	 {
	    		 FastaLocation currloc=iter.next();

	    		 if(lastseq!=currloc.getSeqId())
	    		 {
	    			 seqcount+=1;
	    			
	    			 if(lastseq!=-1)
	    			 {
	    			 Sorted_labels.put(maxseq_score-seqcount*common.DoubleMinNormal, 0);
	    			
	    			 }
	    			 maxseq_score=currloc.Score;
	    			 lastseq=currloc.getSeqId();
	    		 }
	    		 if(maxseq_score<currloc.Score)
	    		 {
	    			 maxseq_score=currloc.Score;
	    		 }
	    		 
	    	 }
	       	Sorted_labels.put(maxseq_score-seqcount*common.DoubleMinNormal, 0);
	       	 int[]  labels=new int[Sorted_labels.size()];
	       	double[]  scores=new double[Sorted_labels.size()];
	       	 int ii=0;
	       	 int one=0;
	       	 for(Double key:Sorted_labels.descendingKeySet())
	       	 {
	       		 labels[ii]=Sorted_labels.get(key);
	       		 if(labels[ii]==1)
	       			 one++;
	       		 scores[ii]=key;
	       		        ii++;
	       	 }
        	 
        	 //AUCcalc.addROCPoint(fp,(double)seqcount/SearchEngine.getSeqNum());
	       	Confusion AUCcalc=AUCCalculator.readArrays(labels, scores);	
//	       	System.out.println(scores[scores.length-1]);
         double AUCscore=AUCcalc.calculateAUCROC();
         
         return AUCscore;
		
	}
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		Options options = new Options();
		options.addOption("i", true, "input file");
		options.addOption("pbm", false, "input file is PBM format (default is fasta format)");
		options.addOption("pwm", true, "input PWM file");
		options.addOption("c", true, "control fasta file");
		options.addOption("version", true, "the version of transformer (1,2,3...)");
		options.addOption("bgmodel", true, "background model file");
		options.addOption("site", false, "input sequences are binding sites, no need to provide PWM");
		options.addOption("prefix", true, "output directory");
		options.addOption("flank", true, "the number of flanking positions around PWM to include(default 0)");
		options.addOption("ratio",true, "sampling ratio (default 1)");
		options.addOption("threadnum",true, "number of threads to use(default :4)");
		options.addOption("paranum",true, "free parameter number for each position(default :3)");
		options.addOption("thresh",true, "minimum entropy threshold for considering a position as a gap(default 1)");
		options.addOption("KLthresh",true, "at least percentage of descrease to claim as a dependency (default 0.01)");
		options.addOption("oops",false,"whether assuming only one occurrence per sequence (default false)");
		options.addOption("oopg",false,"whether assuming only one dependence per gap region (default false)");
		options.addOption("rmbg",false,"whether considering the background probility in learning and evaluating the dependences (default false)");
		options.addOption("maxlen",true,"maxmimum length of gap (default 8)");
		options.addOption("FDR",true,"fasle positive rate");
		String inputPWM="";
		CommandLineParser parser = new GnuParser();
		GapImprover GImprover=new GapImprover();
		double version=2;
		long fullstart = System.currentTimeMillis();
		try {
			CommandLine cmd = parser.parse( options, args);
			if(cmd.hasOption("i"))
			{
				GImprover.inputFasta=cmd.getOptionValue("i");
			}
			else
			{
				throw new ParseException("no input sequence file");
			}
			if(cmd.hasOption("pwm"))
			{
				inputPWM=cmd.getOptionValue("pwm");
			}
			else if(!cmd.hasOption("site"))
			{
				throw new ParseException("no input pwm file");
			}
			if(cmd.hasOption("pbm"))
			{
				GImprover.PBMflag=true;
			}
			if(cmd.hasOption("site"))
			{
				GImprover.is_bsitedata=true;
			}
			
			if(cmd.hasOption("threadnum"))
			{
				GImprover.threadNum=Integer.parseInt(cmd.getOptionValue("threadnum"));
			}
			if(cmd.hasOption("version"))
			{
				version=Double.parseDouble(cmd.getOptionValue("version"));
			}
			if(cmd.hasOption("paranum"))
			{
				GapBGModelingThread.ParaNum=Integer.parseInt(cmd.getOptionValue("paranum"));
			}

			if(cmd.hasOption("c"))
			{
				GImprover.ctrlFasta=cmd.getOptionValue("c");
			}
			if(cmd.hasOption("bgmodel"))
			{
				GImprover.bgmodelFile=cmd.getOptionValue("bgmodel");
			}
			if(cmd.hasOption("prefix"))
			{
				GImprover.outputPrefix=cmd.getOptionValue("prefix");
			}
			if(cmd.hasOption("flank"))
			{
				GImprover.FlankLen=Integer.parseInt(cmd.getOptionValue("flank"));
			}

			if(cmd.hasOption("ratio"))
			{
				GImprover.sampling_ratio=Double.parseDouble( cmd.getOptionValue("ratio"));
			}
			if(cmd.hasOption("thresh"))
			{
				GImprover.entropyThresh=Double.parseDouble(cmd.getOptionValue("thresh"));
			}
			if(cmd.hasOption("KLthresh"))
			{
				GImprover.KLthresh=Double.parseDouble(cmd.getOptionValue("KLthresh"));
			}
			if(cmd.hasOption("maxlen"))
			{
				GImprover.max_gaplen=Integer.parseInt(cmd.getOptionValue("maxlen"));
			}
			if(cmd.hasOption("FDR"))
			{
				GImprover.FDR=Double.parseDouble(cmd.getOptionValue("FDR"));
			}
			if(cmd.hasOption("oops"))
			{
				GImprover.OOPS=true;
			}
			if(cmd.hasOption("oopg"))
			{
				GImprover.OOPG=true;
			}
			if(cmd.hasOption("rmbg"))
			{
				GImprover.removeBG=true;
			}
			System.out.println(Arrays.toString(args) );
		} catch (ParseException e) {
			// TODO Auto-generated catch block
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp( "GapImprover", options );
			return;
		}
		int oldflank=GImprover.FlankLen;
		GImprover.initialize();
		PWM.infothresh=0;
		LinkedList<PWM> pwmlist=null;
		if(!GImprover.is_bsitedata)
		{
			pwmlist=common.LoadPWMFromFile(inputPWM);
			Iterator<PWM> iter=pwmlist.iterator();
			LinkedList<GapPWM> improvedPWMs=new LinkedList<GapPWM>();
			 File directory=new File(inputPWM);
			 String outdir=directory.getParent();
			 if(GImprover.outputPrefix!="./")
				 outdir=GImprover.outputPrefix;
			File file = new File(outdir+"/GPimprover_sorted.dpwm"); 
			TreeMap<Double, PWM> sortedPWMs=new TreeMap<Double, PWM>();
			try {
				BufferedWriter writer= new BufferedWriter(new FileWriter(file));
				while(iter.hasNext())
				{
					PWM rawpwm=iter.next();
					System.out.println(rawpwm.Name+" : "+rawpwm.Consensus(true));
					GImprover.FlankLen=oldflank;
			
					if(GImprover.removeBG)
						GImprover.SearchEngine.EnableBackground(GImprover.background);
					
					String oldname=rawpwm.Name;
					int lastcount=0;
					int currentcount=10;
					int loopcount=0;
					
					GImprover.FlankLen=0;
					//prior PWM refine loop
					while((currentcount-lastcount)>1&&loopcount<30)
					{
						lastcount=currentcount;
						loopcount++;
						KeyValuePair<LinkedList<String>, ArrayList<Double>> retpair=GImprover.querySites(rawpwm);
						if(retpair==null)
						{
							System.out.println("Bad Motif: "+oldname);
							rawpwm=null;
							break;
						}
						currentcount=retpair.key.size();
						try {
							rawpwm=new PWM(retpair.key.toArray(new String[1]));
						} catch (IllegalAlphabetException e) {
							// TODO Auto-generated catch block
							e.printStackTrace();
						} catch (IllegalSymbolException e) {
							// TODO Auto-generated catch block
							e.printStackTrace();
						}
					}
					GImprover.FlankLen=oldflank;
					if(rawpwm!=null)
					{
						GapPWM gpwm=null;
						if(version==2)
						{
							gpwm=GImprover.fillDependency2(rawpwm);
							GImprover.FlankLen=0;
							gpwm=GImprover.fillDependency2(gpwm);
						}
						if(version==1)
						{
							gpwm=GImprover.fillDependency1(rawpwm);
							GImprover.FlankLen=0;
							gpwm=GImprover.fillDependency1(gpwm);
						}
						if(version==1.1)
						{
							gpwm=GImprover.fillDependency1_1(rawpwm);
							GImprover.FlankLen=0;
							gpwm=GImprover.fillDependency1_1(gpwm);
						}
						if(version==2.1)
						{
							gpwm=GImprover.fillDependency2_1(rawpwm);
							GImprover.FlankLen=0;
							gpwm=GImprover.fillDependency2_1(gpwm);
						}
	
						if(rawpwm.Prior_EZ<5)//the conserved bases number less than 5
							gpwm=GImprover.refineGapPWM(gpwm);
						gpwm.Name="GPimpover_"+oldname;
						GImprover.SearchEngine.DisableBackground();
						System.out.print("Original Motif:");
						GImprover.AUCtest(rawpwm);
						System.out.print("Improved Motif:");
						double score=GImprover.AUCtest(gpwm);
						gpwm.Score=score;
						sortedPWMs.put(score+gpwm.Name.hashCode()*common.DoubleMinNormal, gpwm);
						if(GImprover.PBMflag)
						{
							double corr1=GImprover.CorrelationTest(rawpwm);
							System.out.println("Original Motif Correlation with Signal:"+corr1);
							
							double corr2=GImprover.CorrelationTest(gpwm);
							System.out.println("Improved Motif Correlation with Signal:"+corr2);
						}
					}
					
				}
				System.out.println("Writing the DPWM file: "+sortedPWMs.size());
				for(Double key:sortedPWMs.descendingKeySet())
				{
					writer.write(sortedPWMs.get(key).toString());
					
					//output the site PWM 
					PWM dataPWM=sortedPWMs.get(key).subPWM(sortedPWMs.get(key).head, sortedPWMs.get(key).head+sortedPWMs.get(key).core_motiflen);
					dataPWM.Name="sitePWM";
					writer.write(dataPWM.toString());
				}
				writer.close();
				System.out.println("Finish Writing the DPWM file");
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		else
		{
			 File directory=new File(GImprover.inputFasta);
			 String outdir=directory.getParent();
			 if(GImprover.outputPrefix!="./")
				 outdir=GImprover.outputPrefix;
			File file = new File(outdir+"/GPimprover_sorted.dpwm"); 
			try {
				BufferedWriter writer= new BufferedWriter(new FileWriter(file));
				System.out.println("binding sites mode");
				if(GImprover.removeBG)
					GImprover.SearchEngine.EnableBackground(GImprover.background);
				GapPWM gpwm=GImprover.fillDependency2(null);
				gpwm.Name="GPimpover_sitePWM";
				GImprover.SearchEngine.DisableBackground();
				writer.write(gpwm.toString());
				
				
				PWM dataPWM=new PWM(GImprover.SearchEngine.ForwardStrand.toArray(new String[1]));
				dataPWM.Name="sitePWM";
				writer.write(dataPWM.toString());
				System.out.print("Original Motif:");
				GImprover.AUCtest(dataPWM);
				System.out.print("Improved Motif:");
				double score=GImprover.AUCtest(gpwm);
				gpwm.Score=score;
				
				if(GImprover.PBMflag)
				{
					double corr1=GImprover.CorrelationTest(dataPWM);
					System.out.println("Original Motif Correlation with Signal:"+corr1);
					
					double corr2=GImprover.CorrelationTest(gpwm);
					System.out.println("Improved Motif Correlation with Signal:"+corr2);
				}
				writer.close();
				
				
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (IllegalAlphabetException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (IllegalSymbolException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		
		long fullend = System.currentTimeMillis();
		
		System.out.println("Total time was "+(fullend-fullstart)/1000+" seconds.");
		
		System.exit(0);
		
		
	}

}
