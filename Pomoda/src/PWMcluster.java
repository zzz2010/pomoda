import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

import org.apache.commons.cli.Options;


public class PWMcluster {

	LinearEngine SearchEngine;
	double threshold;
	public double sampling_ratio=0.8;
	public double FDR=0.01;
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
		
		
	}
	
	
	public ArrayList<PWM> Clustering(List<PWM> rawPwms,int num_cluster)
	{
		ArrayList<PWM> clusterMoitfs=new ArrayList<PWM>(num_cluster);
		for (int i = 0; i < clusterMoitfs.size(); i++) {
			PWM rawpwm=rawPwms.get(i);
			double thresh=rawpwm.getThresh(sampling_ratio, FDR, background);
			LinkedList<FastaLocation> falocs=SearchEngine.searchPattern(rawpwm, thresh);
			
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
