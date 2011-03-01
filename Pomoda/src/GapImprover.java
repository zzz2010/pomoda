import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedList;

import org.biojava.bio.dist.DistributionTools;
import org.pr.clustering.hierarchical.LinkageCriterion;


public class GapImprover {

	/**
	 * @param args
	 */
	
	LinearEngine SearchEngine;
	public double sampling_ratio=0.8;
	public double FDR=0.01;
	public double entropyThresh=0.5;
	public BGModel background;
	public GapImprover(Pomoda motiffinder)
	{
		SearchEngine=motiffinder.SearchEngine2;
		sampling_ratio=motiffinder.sampling_ratio;
		FDR=motiffinder.FDR;
		background=motiffinder.background;
		
	}
	
	public GapPWM fillDependency(PWM motif)
	{
		String Consensus=motif.Consensus(true);
		ArrayList<Integer> gapstart=new ArrayList<Integer>(motif.core_motiflen/2);
		ArrayList<Integer> gapend=new ArrayList<Integer>(motif.core_motiflen/2);
		int FlankLen=Math.min(Math.min(motif.head, motif.tail), 2);
		//detect gap range
		int start=-1;
		for (int i = motif.head-FlankLen; i < motif.head+motif.core_motiflen+FlankLen; i++) {
			double entropy=DistributionTools.totalEntropy(motif.getColumn(i)) ;
			if(entropy>entropyThresh)
			{
				if(start==-1)
				{
					start=i-motif.head-FlankLen;
				}
			}
			else
			{
				if(start!=-1)
				{
					gapstart.add(start);
					gapend.add(i-motif.head-FlankLen);
					start=-1;
				}
			}
		}
		//get a set of instance strings
		LinkedList<String> sites=new LinkedList<String>();
		double pwmThresh=motif.getThresh(sampling_ratio, FDR, background);
		LinkedList<FastaLocation> falocs=SearchEngine.searchPattern(motif, pwmThresh);
		Iterator<FastaLocation> iter=falocs.iterator();
		while(iter.hasNext())
		{
			FastaLocation currloc=iter.next();
			String site=SearchEngine.getSite(currloc.getSeqId(), currloc.getSeqPos()-FlankLen, motif.core_motiflen+2*FlankLen);
			if(currloc.ReverseStrand)
				site=common.getReverseCompletementString(site);
			sites.add(site);		
		}
		
		
		return null;
	}
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}
