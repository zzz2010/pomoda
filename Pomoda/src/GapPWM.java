import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Map;

import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dist.UniformDistribution;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;


public class GapPWM extends PWM {
	HashMap<Integer,HashMap<String,Double>> Dgroup_DmerProb;
	int[] GroupId;
	int Num_Group;
	
	public GapPWM(String[] alignments) throws IllegalAlphabetException,
			IllegalSymbolException {
		super(alignments);
		// TODO Auto-generated constructor stub
	}
	public GapPWM(Distribution[] dists) throws IllegalAlphabetException
	{
		super(dists);
	}
	

	public static GapPWM createGapPWM(PWM pwm,HashMap<HashSet<Integer>,HashMap<String,Double>> Dmap, int FlankLen )
	{
		GapPWM ret=null;
		Distribution[] dists=new Distribution[pwm.columns()+2*FlankLen];
		for (int i = 0; i < dists.length; i++) {
			if(i<FlankLen)
				dists[i]=new UniformDistribution(DNATools.getDNA());
			else if((pwm.columns()+FlankLen)<=i)
				dists[i]=new UniformDistribution(DNATools.getDNA());
			else
				dists[i]=pwm.getColumn(i-FlankLen);
		}
		try {
			ret=new GapPWM(dists);
			ret.core_motiflen=pwm.core_motiflen;
			ret.head=pwm.head;
			ret.tail=pwm.tail;
			ret.Name=pwm.Name;
			ret.pos_prior=(ArrayList<Double>) pwm.pos_prior.clone();
			int gid=1;
			ret.GroupId=new int[pwm.columns()+2*FlankLen];
			ret.Dgroup_DmerProb=new HashMap<Integer,HashMap<String,Double>>();
			for(HashSet<Integer> key:Dmap.keySet())
			{
				Iterator<Integer>iter=key.iterator();
				while(iter.hasNext())
				{
					int dpos=iter.next();
					ret.GroupId[dpos]=gid;					
					ret.head=Math.min(dpos, pwm.head+FlankLen);
					ret.tail=Math.min(ret.columns()-dpos-1,pwm.tail+FlankLen);
				}
				if(Dmap.get(key).size()>0)
				{
				ret.Dgroup_DmerProb.put(gid, Dmap.get(key));
				gid++;
				}
			}
			ret.Num_Group=gid-1;
			ret.core_motiflen=ret.columns()-ret.head-ret.tail;
		} catch (IllegalAlphabetException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return ret;
		
	}
	
	
	public double getThresh(double sampleratio,double FDRthresh,BGModel bgmodel)
	{
		
		double log025=Math.log(0.25);
	
		//number sampling
		int num_sampl=100000;

			ArrayList<Double> scorelist=new ArrayList<Double>(num_sampl);
			
			int count=0;
			//double sumfdr=0;
			//double sumProb=0;
			while(count<num_sampl)
			{
			 KeyValuePair<Double, String>	sample=bgmodel.generateRandomSequence(core_motiflen);
				//sumfdr+=sample.getKey();
				double score=scoreWeightMatrix(sample.getValue());
				double score2=scoreWeightMatrix(common.getReverseCompletementString(sample.getValue()));
				if(score<score2)
					score=score2;
				scorelist.add(score);
				count++;
				//sumProb+=Math.exp(score);
				
			}
			Collections.sort(scorelist);
			double thresh=scorelist.get((int)Math.floor(scorelist.size()*(1-FDRthresh)));
			if(thresh==scorelist.get(scorelist.size()-1))
				thresh-=common.DoubleMinNormal;
			return thresh;
			
	
		
		
	}
	
	
	//only consider the core-part, ignore flanking , log score
	public double scoreWeightMatrix( String seq)
	{
		double score=0;
		String[] dmers=new String[Num_Group];
		for (int i = 0; i < Num_Group; i++) {
			dmers[i]="";
		}
		
		
		int len=Math.min(core_motiflen, seq.length());
         for (int i = 0; i < len; i++) {
        	 if(common.acgt(seq.charAt(i))>3)
        		 return Double.NEGATIVE_INFINITY;
        	 if(GroupId[head+i]==0)
        		 score+=log_matrix[head+i][common.acgt(seq.charAt(i))];
        	 else
        		 dmers[GroupId[head+i]-1]+=seq.charAt(i);
		}
         for (int i = 0; i < Num_Group; i++) {
        	 dmers[i]=dmers[i].toUpperCase();
			if(Dgroup_DmerProb.get(i+1).containsKey(dmers[i]))
			{
				score+=Math.log(Dgroup_DmerProb.get(i+1).get(dmers[i]));
			}
			else
			{
				score+=Math.log(Dgroup_DmerProb.get(i+1).get("N"));
			}
		}
         
         return score;
		
	}

}
