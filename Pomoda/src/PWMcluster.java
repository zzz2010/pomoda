import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

import org.apache.commons.cli.Options;
import org.pr.clustering.hierarchical.Hierarchical;
import org.pr.clustering.hierarchical.LinkageCriterion;


public class PWMcluster {

	LinearEngine SearchEngine;
	double threshold;
	public double sampling_ratio=0.8;
	public double FDR=0.01;
	public LinkageCriterion linkage;
	public BGModel background;
	public PWMcluster(String faFile,double thresh)
	{
		SearchEngine=new LinearEngine(6);
		SearchEngine.build_index(faFile);
		threshold=thresh;
	}
	
	public PWMcluster(Pomoda motiffinder)
	{
		SearchEngine=motiffinder.SearchEngine2;
		sampling_ratio=motiffinder.sampling_ratio;
		FDR=motiffinder.FDR;
		background=motiffinder.background;
		linkage=LinkageCriterion.valueOf(motiffinder.linkage);
		System.out.println("LinkageCriterion: "+motiffinder.linkage);
		
	}
	
	
	public ArrayList<PWM> Clustering(List<PWM> rawPwms,int num_cluster)
	{
		ArrayList<PWM> clusterMoitfs=new ArrayList<PWM>(num_cluster);
		ArrayList<Thread> threadpool=new ArrayList<Thread>(rawPwms.size()*rawPwms.size());
		SearchEngine.DisableBackground();
		for (int i = 0; i <rawPwms.size(); i++) {
			PWM rawpwm=rawPwms.get(i);
			System.out.println(rawpwm.Consensus(true)+'\t'+rawpwm.Score);
			double thresh=rawpwm.getThresh(sampling_ratio, FDR, background);
			LinkedList<FastaLocation> falocs=SearchEngine.searchPattern(rawpwm, thresh);
			ArrayList<Integer> pos=new ArrayList<Integer>(falocs.size());
			Iterator<FastaLocation> iter=falocs.iterator();
			while(iter.hasNext())
			{
				pos.add((iter.next().getMin()+rawpwm.columns()/2));
			}
			SortingThread t1=new SortingThread(pos);
			t1.start();
			threadpool.add(t1);
			
		}
		try {
				ArrayList<LinkedList<Integer>> PosSet=new ArrayList<LinkedList<Integer>>(rawPwms.size());
				for (int i = 0; i < threadpool.size(); i++) {
				  SortingThread t1=(SortingThread)threadpool.get(i);
				  t1.join();
				PosSet.add((LinkedList<Integer>)t1.getResult());
				}
				threadpool.clear();
				for (int i = 0; i < rawPwms.size()-1; i++) {
					for (int j = i+1; j < rawPwms.size(); j++) {
						OverlappingThread t2=new OverlappingThread(PosSet.get(i), PosSet.get(j), 10);
						t2.start();
						t2.setName(String.valueOf(i*rawPwms.size()+j));
						threadpool.add(t2);
					}
				}
				double[][] dist=new double[rawPwms.size()][rawPwms.size()];
				for (int i = 0; i < threadpool.size(); i++) {
					OverlappingThread t2=(OverlappingThread)threadpool.get(i);
					int pairid=Integer.parseInt(t2.getName());
					int row=pairid/rawPwms.size();
					int col=pairid%rawPwms.size();						
					t2.join();
					double temp=t2.getResult().size()/(double)Math.min(PosSet.get(row).size(), PosSet.get(col).size());
					dist[row][col]=1-temp; //distance
					dist[col][row]=1-temp;
					
				}
				
//				for (int i = 0; i < dist.length; i++) {
//					StringBuffer sb=new StringBuffer("");
//					for (int j = 0; j< dist.length; j++)
//					{  
//						double weight=dist[i][j];
//							sb.append(String.valueOf(weight));
//							sb.append('\t');
//							
//					}
//					System.out.println(sb.toString());
//
//					
//				}
				
				Hierarchical clustering= new Hierarchical(dist,linkage);
				List<Integer> clusterlabed=clustering.partition(num_cluster);
				for (int i = 0; i <num_cluster; i++)
				{
					clusterMoitfs.add(null);
				}
				for (int i = 0; i <clusterlabed.size(); i++)
				{
					int cid=clusterlabed.get(i);
					
					if(cid>=clusterMoitfs.size())
					{
						int origsize=clusterMoitfs.size();
						for (int j = 0; j < cid+1-origsize; j++) {
							clusterMoitfs.add(null);
						}
					}
					if(clusterMoitfs.get(cid)==null)
					{
						clusterMoitfs.set(cid, rawPwms.get(i));
					}
					else if(clusterMoitfs.get(cid).Score<rawPwms.get(i).Score)
					{
						clusterMoitfs.set(cid, rawPwms.get(i));
					}
				}

				
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		
		
		return clusterMoitfs;
	}
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		Options options = new Options();
		options.addOption("i", true, "input pwm file");
		options.addOption("fa", true, "input fasta file[optional]");
		options.addOption("N", true, "number of cluster motifs[default is 5]");
		

	}

}
