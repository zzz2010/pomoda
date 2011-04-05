import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.TreeMap;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.pr.clustering.hierarchical.Hierarchical;
import org.pr.clustering.hierarchical.LinkageCriterion;


public class PWMcluster {

	LinearEngine SearchEngine;
	public double sampling_ratio=0.8;
	public double FDR=0.01;
	public LinkageCriterion linkage=LinkageCriterion.WPGMA;
	public BGModel background;
	private String bgmodelFile="";
	public String outputPrefix="./";
	public String inputFasta;
	public String ctrlFasta="";

	public PWMcluster()
	{
		
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
	
	public void initialize()
	{

		SearchEngine=new LinearEngine(6);
		SearchEngine.build_index(this.inputFasta);
	
		background=new BGModel();
		File file=null;
		int bg_markov_order=3;
		if(ctrlFasta.isEmpty())
		{
			bg_markov_order=1;
			file= new File(inputFasta+".bgobj");
		}
		else
			file= new File(ctrlFasta+".bgobj");
		
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
		     background.BuildModel(inputFasta, bg_markov_order+1); //1-order bg
		     background.SaveModel(inputFasta+".bgobj");
			}
			else
			{
			     background.BuildModel(ctrlFasta, bg_markov_order+1); //3-order bg
			     background.SaveModel(ctrlFasta+".bgobj");
			}				
		}
		
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
						t2.run();
						t2.setName(String.valueOf(i*rawPwms.size()+j));
						threadpool.add(t2);
					}
					Runtime.getRuntime().gc();
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
		options.addOption("i", true, "input fasta file");
		options.addOption("pwm", true, "input PWM file");
		options.addOption("c", true, "control fasta file");
		options.addOption("convert", false, "convert input PWM file to the transfac format");
		options.addOption("clust",true,"linkage type of hierachical clustering:"+Arrays.toString(LinkageCriterion.values()) );
		options.addOption("match", true, "find similar motifs in known PWM library (path to the library, e.g., jaspar.pwm)");
		options.addOption("bgmodel", true, "background model file");
		options.addOption("prefix", true, "output directory");
		options.addOption("ratio",true, "sampling ratio (default 1)");
		options.addOption("FDR",true,"fasle positive rate");
		options.addOption("N", true, "number of cluster motifs[default is 5]");
		String inputPWM;
		CommandLineParser parser = new GnuParser();

		PWMcluster clustering=new PWMcluster();
		boolean convertflag=false;
		LinkedList<PWM> PWMLibrary=null;
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
			if(cmd.hasOption("prefix"))
			{
				clustering.outputPrefix=cmd.getOptionValue("prefix");
			}

			if(cmd.hasOption("ratio"))
			{
				clustering.sampling_ratio=Double.parseDouble( cmd.getOptionValue("ratio"));
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
		
	   File file = new File(inputPWM+"_clust.pwm"); 
	   try {
		   clustering.initialize();
		BufferedWriter writer = new BufferedWriter(new FileWriter(file));
		LinkedList<PWM> pwmlist=common.LoadPWMFromFile(inputPWM);
		ArrayList<PWM>  clusterPWMs=clustering.Clustering(pwmlist, num_cluster);
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
			System.out.println(sortedPWMs.get(key).Consensus(true)+'\t'+sortedPWMs.get(key).Score);
			writer.write(sortedPWMs.get(key).toString());
		}
		writer.close();
	} catch (IOException e) {
		// TODO Auto-generated catch block
		e.printStackTrace();
	}

	}

}
