import java.io.BufferedReader;
import java.io.IOException;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dist.DistributionFactory;
import org.biojava.bio.dist.UniformDistribution;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.utils.ChangeVetoException;


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
	
	
	public String toString()
	{
		StringBuffer TransStr=new StringBuffer("");
		String consensus=Consensus(false);
		TransStr.append("DE\t"+Name+"\t"+consensus+"\t"+String.valueOf(this.Score)+"\n");
		TransStr.append("PO\tA\tC\tG\tT\n");
		for (int i = 0; i < this.columns(); i++) {
			TransStr.append(i+1);
			TransStr.append('\t');
			for (int j = 0; j< 4; j++)
			{  
				double weight=m_matrix[i][j];
				TransStr.append(String.valueOf(weight));
				TransStr.append('\t');
					
			}
			TransStr.append(consensus.charAt(i)+"\n");
	
		}
		TransStr.append("XXD\n");
		for (int i = 0; i < Num_Group; i++) {
			StringBuffer dpos_str=new StringBuffer("");
			for (int j = 0; j < GroupId.length; j++) {
				if(GroupId[j]==i+1)
				{
					dpos_str.append(String.valueOf(j));
					dpos_str.append("-");
				}
			}
			dpos_str.setCharAt(dpos_str.length()-1, ':');
			HashMap<String, Double> dprob=Dgroup_DmerProb.get(i+1);
			for(String key : dprob.keySet())
			{
				dpos_str.append(key+"|"+dprob.get(key)+"\t");
			}
			dpos_str.setCharAt(dpos_str.length()-1, '\n');
			TransStr.append(dpos_str.toString());
		}
		TransStr.append("XXX\n");
		
		return TransStr.toString();
	}

	public static GapPWM parseTransfac(String transfcontent)
	{
		String str;
		BufferedReader reader = new BufferedReader(	new StringReader(transfcontent));
		ArrayList<Distribution> dists=new ArrayList<Distribution>();
		HashMap<HashSet<Integer>,HashMap<String,Double>> Dmap=new HashMap<HashSet<Integer>, HashMap<String,Double>>();
		GapPWM pwm=null;
		String pwmName="";
		boolean pwmrow=true;
		try {
			while ((str = reader.readLine()) != null) {
				if(str.startsWith("DE"))
				{
					String[] elms=str.split("\t| ");
					if(elms.length>1)
					pwmName=elms[1];
				}
				else if(str.startsWith("PO"))
				{
					continue;
				}
				else if(str.startsWith("XX"))
				{
					pwmrow=false;
				}
				else if(str.startsWith("XXX"))
				{
					break;
				}
				else if(str.length()>2 &&pwmrow==true)
				{
					String[] elms=str.split("\t| ");
					Distribution di= DistributionFactory.DEFAULT.createDistribution(DNATools.getDNA());
					if(elms.length<5)
						break;
					double [] count=new double[4];
					for (int i = 1; i <= 4; i++) {
						count[i-1]=Double.parseDouble(elms[i])+common.DoubleMinNormal;
					}
					common.Normalize(count);
				
						di.setWeight(DNATools.a(), count[0]);
						di.setWeight(DNATools.c(), count[1]);
						di.setWeight(DNATools.g(), count[2]);
						di.setWeight(DNATools.t(), count[3]);
					dists.add(di);
					//
				}
				else if(str.length()>2 &&pwmrow==false)
				{
					String[] elms=str.split(":");
					String[] dpos=elms[0].split("-");
					HashSet<Integer> dpos_set=new HashSet<Integer>();
					HashMap<String, Double> dprob_set=new HashMap<String, Double>();
					for (int i = 0; i < dpos.length; i++) {
						dpos_set.add(Integer.parseInt(dpos[i]));
					}
					
					String[] dprob=elms[1].split("\t");
					for (int i = 0; i < dprob.length; i++) {
						String[] comp=dprob[i].split("\\|");
						dprob_set.put(comp[0], Double.parseDouble(comp[1]));
					}
					Dmap.put(dpos_set, dprob_set);
				}
			}
			
			
			pwm=createGapPWM(new PWM(dists.toArray(new Distribution[1])),Dmap,0);
			pwm.Name=pwmName;
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IllegalAlphabetException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IllegalSymbolException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (ChangeVetoException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		return pwm;
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
			ret.head=pwm.head+FlankLen;
			ret.tail=pwm.tail+FlankLen;
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
					ret.head=Math.min(dpos, ret.head);
					ret.tail=Math.min(ret.columns()-dpos-1,ret.tail);
					
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
	
	
	@Override
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
	@Override
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
