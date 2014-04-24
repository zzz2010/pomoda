import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.TreeMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.pr.clustering.hierarchical.Hierarchical;
import org.pr.clustering.hierarchical.LinkageCriterion;

import umontreal.iro.lecuyer.probdist.HypergeometricDist;


public class PWMcluster {

	LinearEngine SearchEngine;
	public double sampling_ratio=1;
	public double FDR=0.001;
	public LinkageCriterion linkage=LinkageCriterion.WPGMA;
	public BGModel background;
	private String bgmodelFile="";
	public boolean featureFirst=false;
	public double overlapThresh=0.1;
	public String outputPrefix="./";
	public String inputFasta=null;
	public String ctrlFasta="";
	public int min_motiflen=7;

	public PWMcluster()
	{
		
	}
	
	public PWMcluster(SEME motiffinder)
	{
		SearchEngine=motiffinder.SearchEngine2;
		//sampling_ratio=motiffinder.sampling_ratio;
		//FDR=motiffinder.FDR;
		background=motiffinder.background;
		overlapThresh=motiffinder.overlapthresh;
		min_motiflen=motiffinder.min_motiflen;
		linkage=LinkageCriterion.valueOf(motiffinder.linkage);
		System.out.println("LinkageCriterion: "+motiffinder.linkage);
		
	}
	
	public void initialize()
	{
		common.initialize();
		SearchEngine=new LinearEngine(6);
		if(this.inputFasta!=null)
		SearchEngine.build_index(this.inputFasta);
		else
		{
			String  seqstr=BGModel.CreateUniform().generateRandomSequence(100000).getValue();
			SearchEngine.ForwardStrand.add(seqstr);
			SearchEngine.SeqNames.add("random");
			SearchEngine.TotalLen+=seqstr.length();
			SearchEngine. accSeqLen.add(SearchEngine.TotalLen);
		}
	
		background=new BGModel();
		File file=null;
		int bg_markov_order=5;
		if(ctrlFasta.isEmpty())
		{
			bg_markov_order=1;
			file= new File(inputFasta+".bg");
		}
		else
			file= new File(ctrlFasta+".bg");
		
		if(file.exists()||!bgmodelFile.isEmpty())
		{
			if(bgmodelFile.isEmpty())
			background.LoadModel(file.getAbsolutePath());
			else
				background.LoadModel(bgmodelFile);
		}
		else
		{
			if(ctrlFasta.isEmpty())
			{
				if(inputFasta!=null)
				{
		     background.BuildModel(inputFasta, bg_markov_order+1); //1-order bg
		     background.SaveModel(inputFasta+".bg");
				}
				else
				{
					background=background.CreateUniform();
				}
			}
			else
			{
			     background.BuildModel(ctrlFasta, bg_markov_order+1); //3-order bg
			     background.SaveModel(ctrlFasta+".bg");
			}				
		}
		
	}
	
	
	public double[][] PairScore(List<PWM> rawPwms)
	{
		double[][] overlapRatio=new double[rawPwms.size()][rawPwms.size()];
		

		ArrayList<PWM> sortedlist=new ArrayList<PWM>(rawPwms.size());
		ArrayList<LinkedList<Integer>> PosSet=new ArrayList<LinkedList<Integer>>(rawPwms.size());
		//sort the positions
		ExecutorService executor = Executors.newFixedThreadPool(6);
		ArrayList<Thread> threadpool=new ArrayList<Thread>(rawPwms.size());
		SearchEngine.DisableBackground();
		int threadid=0;
		HashMap<Integer,Integer> pwm_tid=new HashMap<Integer,Integer>();
		int pwmid=-1;
		for(PWM rawpwm:rawPwms)
		{
			pwmid++;
			sortedlist.add(rawpwm);
			if(rawpwm.core_motiflen<6)
				continue;
			System.out.println(rawpwm.Consensus(true)+'\t'+rawpwm.Score);
			double thresh=rawpwm.getThresh(sampling_ratio, FDR, background,false);

			LinkedList<FastaLocation> falocs=SearchEngine.searchPattern(rawpwm, thresh);
			ArrayList<Integer> pos=new ArrayList<Integer>(falocs.size());
			rawpwm.matchsite=new LinkedList<Integer>();
			for(FastaLocation floc:falocs)
			{
				rawpwm.matchsite.add(floc.getMin());
			}
			Iterator<FastaLocation> iter=falocs.iterator();
			while(iter.hasNext())
			{
				FastaLocation temp=iter.next();
				pos.add((temp.getMin()+rawpwm.columns()/2));
			}
			SortingThread t1=new SortingThread(pos);
			executor.execute(t1);
			//t1.start();
			threadpool.add(t1);
			pwm_tid.put(pwmid, threadid);
			threadid+=1;

		}
		
		try
		{
			executor.shutdown();
			// Wait until all threads are finish
			while (!executor.isTerminated()) {
				Thread.sleep(3000);
			}
			for (int i = 0; i < threadpool.size(); i++) {
				  SortingThread t1=(SortingThread)threadpool.get(i);
				
				PosSet.add((LinkedList<Integer>)t1.getResult());
				}
			
			int id=0;
			for (int i = 0; i < rawPwms.size()-1; i++) {
				if(!pwm_tid.containsKey(i))
					continue;
				for (int j = i+1; j < rawPwms.size(); j++) {
					if(!pwm_tid.containsKey(j))
						continue;
					int overlaplen=Math.max(rawPwms.get(i).core_motiflen, rawPwms.get(j).core_motiflen)/2;
					OverlappingThread t2=new OverlappingThread(PosSet.get(pwm_tid.get(i)), PosSet.get(pwm_tid.get(j)), overlaplen);
					t2.run();
					double temp=t2.getResult().size()/(double)Math.min(PosSet.get(i).size()+1, PosSet.get(j).size()+1);
					overlapRatio[i][j]=temp;
					overlapRatio[j][i]=overlapRatio[i][j];
				}
			}

		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return overlapRatio;
	}
	
	public ArrayList<PWM> Clustering_(List<PWM> rawPwms,int num_cluster)
	{
		
		ArrayList<PWM> clusterMoitfs=new ArrayList<PWM>(num_cluster);
		TreeMap<Double, PWM> sortedPWMs=new TreeMap<Double, PWM>();
		for (int i = 0; i <rawPwms.size(); i++) {
			double key=rawPwms.get(i).Score;
			if(rawPwms.get(i).pos_en||rawPwms.get(i).peakrank_en)
				key+=10000;
			sortedPWMs.put(key, rawPwms.get(i));
		}
		ArrayList<PWM> sortedlist=new ArrayList<PWM>(rawPwms.size());
		ArrayList<LinkedList<Integer>> PosSet=new ArrayList<LinkedList<Integer>>(rawPwms.size());
		//sort the positions
		ExecutorService executor = Executors.newFixedThreadPool(6);
		ArrayList<Thread> threadpool=new ArrayList<Thread>(rawPwms.size());
		SearchEngine.DisableBackground();
		for(Double key:sortedPWMs.descendingKeySet())
		{
			
			PWM rawpwm=sortedPWMs.get(key);
			sortedlist.add(rawpwm);
			if(rawpwm.core_motiflen<min_motiflen)
				continue;
			System.out.println(rawpwm.Consensus(true)+'\t'+rawpwm.Score);
			double thresh=rawpwm.getThresh(sampling_ratio, FDR, background,false);

			LinkedList<FastaLocation> falocs=SearchEngine.searchPattern(rawpwm, thresh);
			ArrayList<Integer> pos=new ArrayList<Integer>(falocs.size());
			rawpwm.matchsite=new LinkedList<Integer>();
			for(FastaLocation floc:falocs)
			{
				rawpwm.matchsite.add(floc.getMin());
			}
			Iterator<FastaLocation> iter=falocs.iterator();
			while(iter.hasNext())
			{
				FastaLocation temp=iter.next();
				pos.add((temp.getMin()+rawpwm.columns()/2));
			}
			SortingThread t1=new SortingThread(pos);
			executor.execute(t1);
			//t1.start();
			threadpool.add(t1);

		}
		
		try
		{
			executor.shutdown();
			// Wait until all threads are finish
			while (!executor.isTerminated()) {
				Thread.sleep(3000);
			}
			for (int i = 0; i < threadpool.size(); i++) {
				  SortingThread t1=(SortingThread)threadpool.get(i);
				
				PosSet.add((LinkedList<Integer>)t1.getResult());
				}
			
			int id=0;
			ArrayList<Integer> clusterMoitfsId=new ArrayList<Integer>(num_cluster);
			for(Double key:sortedPWMs.descendingKeySet())
			{
			PWM	motif=sortedPWMs.get(key);
			
				if(clusterMoitfsId.size()==0)
				{
					clusterMoitfsId.add(id);	
					clusterMoitfs.add(motif);
				}
				else
				{
					boolean newclass=true;
					for (int i = 0; i < clusterMoitfsId.size(); i++) {
						int row=clusterMoitfsId.get(i);
						int col=id;	
						int overlaplen=Math.max(rawPwms.get(row).core_motiflen, rawPwms.get(col).core_motiflen)/2;
						OverlappingThread t2=new OverlappingThread(PosSet.get(clusterMoitfsId.get(i)), PosSet.get(id), overlaplen);
						t2.run();

						double temp=t2.getResult().size()/(double)Math.min(PosSet.get(row).size()+1, PosSet.get(col).size()+1);
						int l=Math.max(SearchEngine.TotalLen/(rawPwms.get(row).core_motiflen),SearchEngine.TotalLen/rawPwms.get(col).core_motiflen);
						int m=PosSet.get(row).size();
						int k=PosSet.get(col).size();
						int x=t2.getResult().size();
						HypergeometricDist dist=new HypergeometricDist(m,l,k);
						double temp2=dist.cdf(x);
						if(temp>overlapThresh)//temp2>0.5||temp>overlapThresh
						{
							System.out.println("-"+motif.Consensus(true)+"\t"+clusterMoitfs.get(i).Consensus(true)+"\t"+temp);
							newclass=false;
							break;
						}
					}
					if(newclass)
					{
						clusterMoitfsId.add(id);	
						clusterMoitfs.add(motif);
					}
				}
				
				id++;
				if(clusterMoitfs.size()==num_cluster)
					break;
			}
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
			if(clusterMoitfs.size()<num_cluster&&overlapThresh<1)
			{
				overlapThresh*=1.2;
				System.out.println("adjust overlap ratio:"+overlapThresh);
				return Clustering_Fast(sortedlist,num_cluster);
			}
			return clusterMoitfs;
	}
	
	public ArrayList<PWM> Clustering_Fast(List<PWM> sortedPWMs,int num_cluster)
	{
		ArrayList<PWM> clusterMoitfs=new ArrayList<PWM>(num_cluster);
	
		
		int id=0;
		ArrayList<Integer> clusterMoitfsId=new ArrayList<Integer>(num_cluster);
		for(PWM motif:sortedPWMs)
		{
			if(motif.core_motiflen<6)
				continue;
			if(clusterMoitfsId.size()==0)
			{
				clusterMoitfsId.add(id);	
				clusterMoitfs.add(motif);
			}
			else
			{
				boolean newclass=true;
				for (int i = 0; i < clusterMoitfsId.size(); i++) {

					int row=clusterMoitfsId.get(i);
					int col=id;	
					int overlaplen=Math.max(sortedPWMs.get(row).core_motiflen, sortedPWMs.get(col).core_motiflen)/2;
					OverlappingThread t2=new OverlappingThread(sortedPWMs.get(row).matchsite ,sortedPWMs.get(col).matchsite,overlaplen);
					t2.run();
					double temp=t2.getResult().size()/(double)Math.min(sortedPWMs.get(row).matchsite.size()+1, sortedPWMs.get(col).matchsite.size()+1);
					int l=Math.max(SearchEngine.TotalLen/(sortedPWMs.get(row).core_motiflen),SearchEngine.TotalLen/sortedPWMs.get(col).core_motiflen);
					int m=sortedPWMs.get(row).matchsite.size();
					int k=sortedPWMs.get(col).matchsite.size();
					int x=t2.getResult().size();
					//HypergeometricDist dist=new HypergeometricDist(m,l,k);
					//double temp2=dist.cdf(x);
					if(temp>overlapThresh)//temp2>0.5||temp>overlapThresh
					{
						System.out.println("-"+motif.Consensus(true)+"\t"+clusterMoitfs.get(i).Consensus(true)+"\t"+temp);
						newclass=false;
						break;
					}
				}
				if(newclass)
				{
					clusterMoitfsId.add(id);	
					clusterMoitfs.add(motif);
				}
			}
			
			id++;
			if(clusterMoitfs.size()==num_cluster)
				break;
		}
		
		if(clusterMoitfs.size()<num_cluster&&overlapThresh<1)
		{
			overlapThresh*=1.2;
			System.out.println("adjust overlap ratio:"+overlapThresh);
			return Clustering_Fast(sortedPWMs,num_cluster);
		}
		
		return clusterMoitfs;
	}
	
	public ArrayList<PWM> Clustering(List<PWM> rawPwms,int num_cluster)
	{
		ArrayList<PWM> clusterMoitfs=new ArrayList<PWM>(num_cluster);
		ArrayList<Thread> threadpool=new ArrayList<Thread>(rawPwms.size()*rawPwms.size());
		SearchEngine.DisableBackground();
		 ExecutorService executor = Executors.newFixedThreadPool(6);
		for (int i = 0; i <rawPwms.size(); i++) {
			PWM rawpwm=rawPwms.get(i);
			System.out.println(rawpwm.Consensus(true)+'\t'+rawpwm.Score);
			double thresh=rawpwm.getThresh(sampling_ratio, FDR, background,false);
			LinkedList<FastaLocation> falocs=SearchEngine.searchPattern(rawpwm, thresh);
			ArrayList<Integer> pos=new ArrayList<Integer>(falocs.size());
			Iterator<FastaLocation> iter=falocs.iterator();
			while(iter.hasNext())
			{
				FastaLocation temp=iter.next();
				pos.add((temp.getMin()+rawpwm.columns()/2));
//				String ss=SearchEngine.getSite(pos.get(pos.size()-1)/400, pos.get(pos.size()-1)%400-rawpwm.columns()/2, rawpwm.columns());
//				//String ss=SearchEngine.getSite(temp.getSeqId(), temp.getSeqPos(), rawpwm.columns());
//				double score=rawpwm.scoreWeightMatrix(ss);
//				score=Math.max(score, rawpwm.scoreWeightMatrix(common.getReverseCompletementString(ss)));
//				if(score<thresh)
//					score=0;
			}
			SortingThread t1=new SortingThread(pos);
			executor.execute(t1);
			//t1.start();
			threadpool.add(t1);
			
		}
		try {
			executor.shutdown();
			// Wait until all threads are finish
			while (!executor.isTerminated()) {
				Thread.sleep(3000);
			}
				ArrayList<LinkedList<Integer>> PosSet=new ArrayList<LinkedList<Integer>>(rawPwms.size());
				for (int i = 0; i < threadpool.size(); i++) {
				  SortingThread t1=(SortingThread)threadpool.get(i);
				
				PosSet.add((LinkedList<Integer>)t1.getResult());
				}
				executor = Executors.newFixedThreadPool(6);
				threadpool.clear();
				for (int i = 0; i < rawPwms.size()-1; i++) {
					for (int j = i+1; j < rawPwms.size(); j++) {
						OverlappingThread t2=new OverlappingThread(PosSet.get(i), PosSet.get(j), 5);
						executor.execute(t2);
						t2.setName(String.valueOf(i*rawPwms.size()+j));
						threadpool.add(t2);
					}
					
				}
				executor.shutdown();
				// Wait until all threads are finish
				while (!executor.isTerminated()) {
					Thread.sleep(3000);
				}
				double[][] dist=new double[rawPwms.size()][rawPwms.size()];
				for (int i = 0; i < threadpool.size(); i++) {
					OverlappingThread t2=(OverlappingThread)threadpool.get(i);
					int pairid=Integer.parseInt(t2.getName());
					int row=pairid/rawPwms.size();
					int col=pairid%rawPwms.size();						
					//t2.join();
					double temp=t2.getResult().size()/(double)Math.min(PosSet.get(row).size()+1, PosSet.get(col).size()+1);
//					for(int abspos : t2.getResult())
//					{
//						String site=SearchEngine.getSite(abspos/400, abspos%400-15, 30).toLowerCase();
//						double maxscore=-10000;
//						int maxpos=0;
//						for (int j = 0; j < site.length()-rawPwms.get(0).core_motiflen; j++) {
//							double score=rawPwms.get(0).scoreWeightMatrix(site.substring(j, j+rawPwms.get(0).core_motiflen));
//							if(score>maxscore)
//							{
//								maxpos=j;
//								maxscore=score;
//							}
//						}
//						site=site.replaceFirst(site.substring(maxpos, maxpos+rawPwms.get(0).core_motiflen), site.substring(maxpos, maxpos+rawPwms.get(0).core_motiflen).toUpperCase());
//						
//						 maxscore=-10000;
//						 maxpos=0;
//						for (int j = 0; j < site.length()-rawPwms.get(1).core_motiflen; j++) {
//							double score=rawPwms.get(1).scoreWeightMatrix(site.substring(j, j+rawPwms.get(1).core_motiflen));
//							if(score>maxscore)
//							{
//								maxpos=j;
//								maxscore=score;
//							}
//						}
//						site=site.replaceFirst(site.substring(maxpos, maxpos+rawPwms.get(1).core_motiflen), site.substring(maxpos, maxpos+rawPwms.get(1).core_motiflen).toUpperCase());
//						
//						System.err.println(site+"\t"+site.substring(maxpos, maxpos+rawPwms.get(1).core_motiflen));
//						
//					}
					dist[row][col]=1-temp; //distance
					dist[col][row]=1-temp;
					
				}
				
				
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
					System.out.println(rawPwms.get(i).Consensus(true)+"\t"+cid);
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
		options.addOption("i", true, "input fasta file");
		options.addOption("pwm", true, "input PWM file");
		options.addOption("pairscore_only", false, "only output pairscore");
		options.addOption("c", true, "control fasta file");
		options.addOption("convert", false, "convert input PWM file to the transfac format");
		options.addOption("clust",true,"linkage type of hierachical clustering:"+Arrays.toString(LinkageCriterion.values()) );
		options.addOption("match", true, "find similar motifs in known PWM library (path to the library, e.g., jaspar.pwm)");
		options.addOption("bgmodel", true, "background model file");
		options.addOption("prefix", true, "output directory");
		options.addOption("thresh", true, "overlapping threshold to detemine whether cluster (default 0.2)");
		options.addOption("ratio",true, "sampling ratio (default 1)");
		options.addOption("FDR",true,"fasle positive rate");
		options.addOption("N", true, "number of cluster motifs[default is 5]");
		String inputPWM;
		CommandLineParser parser = new GnuParser();

		PWMcluster clustering=new PWMcluster();
		boolean convertflag=false;
		LinkedList<PWM> PWMLibrary=null;
		boolean pairscore_only=false;
		int num_cluster=5;
		try {
			CommandLine cmd = parser.parse( options, args);
			if(cmd.hasOption("i"))
			{
				clustering.inputFasta=cmd.getOptionValue("i");
			}
			if(cmd.hasOption("pwm"))
			{
				inputPWM=cmd.getOptionValue("pwm");
			}
			else
			{
				throw new ParseException("no input pwm file");
			}
			if(cmd.hasOption("c"))
			{
				clustering.ctrlFasta=cmd.getOptionValue("c");
			}
			if(cmd.hasOption("bgmodel"))
			{
				clustering.bgmodelFile=cmd.getOptionValue("bgmodel");
			}
			if(cmd.hasOption("match"))
			{
				PWMLibrary=common.LoadPWMFromFile(cmd.getOptionValue("match"));
			}

			if(cmd.hasOption("convert"))
			{
				convertflag =true;
			}
			if(cmd.hasOption("pairscore_only"))
			{
				pairscore_only =true;
			}
			if(cmd.hasOption("prefix"))
			{
				clustering.outputPrefix=cmd.getOptionValue("prefix");
			}
			if(cmd.hasOption("ratio"))
			{
				clustering.sampling_ratio=Double.parseDouble( cmd.getOptionValue("ratio"));
			}
			if(cmd.hasOption("thresh"))
			{
				clustering.overlapThresh=Double.parseDouble( cmd.getOptionValue("thresh"));
			}
			if(cmd.hasOption("N"))
			{
				num_cluster=Integer.parseInt( cmd.getOptionValue("N"));
			}
			if(cmd.hasOption("FDR"))
			{
				clustering.FDR=Double.parseDouble(cmd.getOptionValue("FDR"));
			}
			if(cmd.hasOption("clust"))
			{
				clustering.linkage=LinkageCriterion.valueOf((cmd.getOptionValue("clust")));
			}

		} catch (ParseException e) {
			// TODO Auto-generated catch block
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp( "PWMcluster", options );
			return;
		}
		
	   
	   try {
		   clustering.initialize();

		List<PWM> pwmlist=common.LoadPWMFromFile(inputPWM).subList(0,10);
		if(pairscore_only)
		{
			File file2 = new File(inputPWM+".pairscore"); 
			double[][] overlapRatio=clustering.PairScore(pwmlist);
			PrintWriter writer2=new PrintWriter(file2);
			String header="";
			for (int i = 0; i < pwmlist.size(); i++) {
				header+="\t"+pwmlist.get(i).Name;
			}
			writer2.println(header);
			for (int i = 0; i < pwmlist.size(); i++) {
				String line=pwmlist.get(i).Name;
				for (int j = 0; j < pwmlist.size(); j++) {
					line+="\t"+overlapRatio[i][j];
				}
				writer2.println(line);
			}
			writer2.close();
			return;
		}
		
	    File file = new File(inputPWM+"_clust.pwm"); 
		BufferedWriter writer = new BufferedWriter(new FileWriter(file));
		ArrayList<PWM>  clusterPWMs=clustering.Clustering_(pwmlist, num_cluster);
		TreeMap<Double, PWM> sortedPWMs=new TreeMap<Double, PWM>();
		for(PWM pwm:clusterPWMs)
		{
			sortedPWMs.put(pwm.Score, pwm); //desc order
		}
		int c=0;
		for(Double key:sortedPWMs.descendingKeySet())
		{
			sortedPWMs.get(key).Name="Motif_clust"+String.valueOf(c+1);
			c++;
			System.out.println(sortedPWMs.get(key).Name+":"+sortedPWMs.get(key).Consensus(true)+'\t'+sortedPWMs.get(key).Score);
			writer.write(sortedPWMs.get(key).toString());
		}
		writer.close();
	} catch (IOException e) {
		// TODO Auto-generated catch block
		e.printStackTrace();
	}

	}

}
