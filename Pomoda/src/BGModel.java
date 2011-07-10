import java.io.*;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Random;

import org.biojava.bio.BioException;
import org.biojava.bio.dp.IllegalTransitionException;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.seq.io.*;
import org.biojava.bio.symbol.*;
import org.biojavax.bio.seq.RichSequence;

import com.sun.corba.se.impl.encoding.OSFCodeSetRegistry.Entry;



public class BGModel implements Serializable{
	
	
	HashMap<String, Double> conditionProb;
	HashMap<String, PWM> kmerPWMBG=null;
	boolean EnablePWMBG=false;
	int kmerlen=5;
	int flanklen=25;
	 Random r=new Random();
	public int order;
	
	//return probability
	public KeyValuePair<Double, String> generateRandomSequence(int len)
	{
		double score=1;
		String seq="";
		String ACGT="ACGT";
		StringBuffer sb=new StringBuffer("");
		for (int i = 0; i < len; i++) {
			sb.append('N');
			double cut=r.nextDouble();
			double sumprob=0;
			double tempprob=0;
			for (int j = 0; j < 4; j++) {
				sb.setCharAt(i, ACGT.charAt(j));
				String temp=sb.toString();
				if(order<temp.length())
					temp=temp.substring(temp.length()-order);
				tempprob=conditionProb.get(temp);
				sumprob+=tempprob;
				if(sumprob>cut)
				{
					break;
				}
			}
			score*=tempprob;
		}
		seq=sb.toString();
		
		KeyValuePair<Double, String> result =new KeyValuePair<Double, String>(score, seq);
		return result;
	}
	
	//compute logprob of a sequence, support gap-sequence
	public double Get_LOGPROB(String seq)
	{
		seq=seq.toUpperCase();
		if(seq.contains("X"))
			return Double.POSITIVE_INFINITY;
		double logprob=0;
		int start=0;
		int i;
		while(start<seq.length())
		{
			while(start<seq.length()&&seq.charAt(start)=='N')
				start++;
			for (i = 0; i < order&&i<(seq.length()-start); i++) {
				if(seq.charAt(start+i)=='N')
				{
					logprob+=Math.log(0.25);
					break;
				}
				logprob+=Math.log( conditionProb.get(seq.substring(start,start+i+1)));

			}
			if(start+i==seq.length())
				break;
			if(start+i<seq.length()&& seq.charAt(start+i)=='N')
			{
				start=start+i+1;
					logprob+=Math.log(0.25);
				continue;
			}
			
			for (i = start+1; i < (seq.length()-order+1); i++) {
				start=i+order-1;
				if(seq.charAt(start)=='N')
				{
					logprob+=Math.log(0.25);
					break;
				}
				logprob+=Math.log( conditionProb.get(seq.substring(i,i+order)));
			}
			start++;
		}
		//prob=path.getScore();
	    return logprob;
		
	}
	
	

	public void BuildModel(String[] Seqs, int order)
	{
		this.order=order;
		initializeConditionProb();
		 for (int i = 0; i < Seqs.length; i++) {
			addCount(Seqs[i]);
			addCount(common.getReverseCompletementString(Seqs[i]));
		}
		 
		 normalizeConditionProb();
	}
	
	public void BuildModel(Map<String,Double> Seqs, int order)
	{
		this.order=order;
		initializeConditionProb();
		for(String seq : Seqs.keySet())
		{
			addCount(seq,Seqs.get(seq));
			addCount(common.getReverseCompletementString(seq),Seqs.get(seq));
		}
		 
		 normalizeConditionProb();
		
	}
	
	private void addCount(String seq)
	{
		seq=seq.toUpperCase().replace("N", "");
		for (int i = 0; i < seq.length()-order+1; i++) {
			for (int j = 0; j < order; j++) {
				String segment=seq.substring(i, i+j+1);
//			   if(conditionProb.containsKey(segment))
//			   {
				   conditionProb.put(segment, conditionProb.get(segment)+1);
//			   }
//			   else
//			   {
//				   conditionProb.put(segment, 1.0);
//			   }
			}
			
			if(EnablePWMBG&&seq.length()>(2*flanklen+kmerlen))
			{
				String kmer=seq.substring(i, i+kmerlen);
				String seq1="";
				if(i<flanklen)
				{
					String N="";
					for (int j = 0; j < (flanklen-i); j++) {
						N+="N";
					}
					seq1=N+seq.substring(0,i+kmerlen+flanklen);
				}
				else if((seq.length()-i-kmerlen)<flanklen)
				{
					String N="";
					for (int j = 0; j < (flanklen+i-seq.length()+kmerlen); j++) {
						N+="N";
					}
					seq1=seq.substring(i-flanklen,seq.length())+N;
				}
				else
				{
					seq1=seq.substring(i-flanklen,i+kmerlen+flanklen);
				}
				if(kmerPWMBG.containsKey(kmer))
				{
					PWM temp=kmerPWMBG.get(kmer);
					for (int j = 0; j < temp.columns(); j++) {
						int symid=common.acgt(seq1.charAt(j));
						if(symid>3)
						{
							for (int k = 0; k < 4; k++) {
								temp.m_matrix[j][k]+=0.25;
								
							}
						}
						else
						{
							temp.m_matrix[j][symid]+=1;
						}
					}
					//kmerPWMBG.put(kmer, temp);
				}
				else
				{
					try {
						PWM temp=new PWM(new String[]{seq1});
						kmerPWMBG.put(kmer, temp);
					} catch (IllegalAlphabetException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					} catch (IllegalSymbolException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
				
				
			}
			
		}
		
	}
	private void addCount(String seq,double prob)
	{
		seq=seq.toUpperCase().replace("N", "");
		seq=seq.toUpperCase().replace("X", "");
		for (int i = 0; i < seq.length()-order+1; i++) {
			for (int j = 0; j < order; j++) {
				String segment=seq.substring(i, i+j+1);
				   conditionProb.put(segment, conditionProb.get(segment)+prob);
			}
			
		}
		
	}
	
	private void normalizeConditionProb()
	{
		double totalsample=conditionProb.get("A")+conditionProb.get("C")+conditionProb.get("G")+conditionProb.get("T");
		
		for (int i = 0; i < order; i++) {
			int total=1<<(2*(order-i));;
			for (int j = 0; j < total; j++) {
				String pattern=common.Hash2ACGT(j,order-i);
				if(i==order-1)
					conditionProb.put(pattern, conditionProb.get(pattern)/totalsample); 
				else
				conditionProb.put(pattern, conditionProb.get(pattern)/conditionProb.get(pattern.substring(0, pattern.length()-1)));
			}
		}
		if(EnablePWMBG)
		{
			Iterator<java.util.Map.Entry<String, PWM>> iter=kmerPWMBG.entrySet().iterator();
			while(iter.hasNext())
			{
				java.util.Map.Entry<String, PWM> curr=iter.next();
				String kmer=curr.getKey();
				PWM temp=curr.getValue();
				for (int i = 0; i < temp.columns(); i++) {
					temp.setWeights(i, common.Normalize(temp.m_matrix[i]));
				}
			}
		}
		
	}
	
	private void initializeConditionProb()
	{	
		if(EnablePWMBG)
		{
			kmerPWMBG=new HashMap<String, PWM>();
		}
		conditionProb=new HashMap<String, Double>();
		for (int i = 0; i < order; i++) {
			int total=1<<(2*(i+1));;
			for (int j = 0; j < total; j++) {
				String pattern=common.Hash2ACGT(j,i+1);
				conditionProb.put(pattern, 1.0/(total)); //pseudo count 1/total
			}
			
		}
		
	}
	

	
	
	public void BuildModel(String inputFasta, int order)
	{
		 try {
			 this.order=order;
			 initializeConditionProb();
			 //Database to hold the training set
			      BufferedReader br = new BufferedReader(new FileReader(inputFasta));
   		          SymbolTokenization toke = AlphabetManager.alphabetForName("DNA").getTokenization("token");
   		          SequenceIterator seqi = RichSequence.IOTools.readFasta(br, toke,null);
   			      while (seqi.hasNext()) {
   			    	  Sequence seq=seqi.nextSequence();
   				     addCount(seq.seqString());
   				  addCount(common.getReverseCompletementString(seq.seqString()));
   			      }
   			  normalizeConditionProb();
			     
			
		} catch (IllegalSymbolException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IllegalTransitionException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IllegalAlphabetException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IllegalArgumentException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (BioException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
	
	public void SaveModel(String filename)
	{
		

		try {
			FileOutputStream fo=new FileOutputStream(filename);
	        ObjectOutputStream oos;
			oos = new ObjectOutputStream(fo);
			oos.writeObject(this);
			oos.close();			
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}  
        
        
	}
	public void LoadModel(String filename)
	{
		try {
			FileInputStream fi=new FileInputStream(filename);
			ObjectInputStream si=new ObjectInputStream(fi);
			BGModel loaded=(BGModel)si.readObject();
			conditionProb=loaded.conditionProb;
			order=loaded.order;
			kmerPWMBG=loaded.kmerPWMBG;
			kmerlen=loaded.kmerlen;
			flanklen=loaded.flanklen;
			EnablePWMBG=loaded.EnablePWMBG;
			si.close();
			
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (ClassNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IllegalArgumentException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} 
		
	}
	

}
