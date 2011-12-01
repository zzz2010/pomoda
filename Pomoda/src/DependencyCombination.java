import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
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
	
	
	
}
