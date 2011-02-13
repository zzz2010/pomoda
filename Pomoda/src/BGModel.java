import java.io.*;
import java.util.HashMap;

import org.biojava.bio.BioException;
import org.biojava.bio.dp.*;
import org.biojava.bio.dist.DistributionFactory;
import org.biojava.bio.dp.DP;
import org.biojava.bio.dp.DPFactory;
import org.biojava.bio.dp.IllegalTransitionException;
import org.biojava.bio.dp.ProfileHMM;
import org.biojava.bio.dp.SimpleModelTrainer;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.seq.db.*;
import org.biojava.bio.seq.io.*;
import org.biojava.bio.symbol.*;
import org.biojavax.bio.seq.RichSequence;


public class BGModel implements Serializable{
	
	
	HashMap<String, Double> conditionProb;
	public int order;
	
	//compute logprob of a sequence, support gap-sequence
	public double Get_LOGPROB(String seq)
	{
		double logprob=0;
		int start=0;
		int i;
		while(start<seq.length())
		{
			while(seq.charAt(start)=='N')
				start++;
			for (i = 0; i < order&&i<(seq.length()-start); i++) {
				if(seq.charAt(start+i)=='N')
					break;
				logprob+=Math.log( conditionProb.get(seq.substring(start,start+i+1)));

			}
			if(start+i<seq.length()&& seq.charAt(start+i)=='N')
			{
				start=start+i+1;
				break;
			}
			
			for (i = start+1; i < (seq.length()-order+1); i++) {
				start=i+order-1;
				if(seq.charAt(start)=='N')
					break;
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
	
	private void addCount(String seq)
	{
		seq=seq.toUpperCase();
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
		
	}
	
	private void initializeConditionProb()
	{
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
   				     addCount(seqi.nextSequence().seqString());
   				  addCount(common.getReverseCompletementString(seqi.nextSequence().seqString()));
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
