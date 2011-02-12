import java.io.*;

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


public class BGModel {
	
	DP dp;
	public double Get_LOGPROB(String seq)
	{
	    //decode the most likely state path and produce an 'odds' score
	    StatePath path = null;
		try {
			SymbolList test = DNATools.createDNA(seq);
			
		    /*
		     * put the test sequence in an array, an array is used because for pairwise
		     * alignments using an HMM there would need to be two SymbolLists in the 
		     * array
		     */
		 
		    SymbolList[] sla = {test};
			path = dp.viterbi(sla, ScoreType.PROBABILITY);
		} catch (IllegalSymbolException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IllegalArgumentException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IllegalAlphabetException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IllegalTransitionException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		double prob=path.getScore();
	    return prob;
		
	}
	
	
	
	public void BuildModel(String inputFasta, int order)
	{
		 try {
			ProfileHMM hmm = new ProfileHMM(DNATools.getDNA(),
			         order,
			         DistributionFactory.DEFAULT,
			         DistributionFactory.DEFAULT,
			         "DNA Markov Model");
			dp = DPFactory.DEFAULT.createDP(hmm);
			ModelTrainer mt = new SimpleModelTrainer();
			mt.registerModel(hmm);
			mt.setNullModelWeight(1.0);
			mt.train();
			BaumWelchTrainer bwt = new BaumWelchTrainer(dp);
			    //anonymous implementation of the stopping criteria interface to stop after 20 iterations
			StoppingCriteria stopper = new StoppingCriteria(){
			      public boolean isTrainingComplete(TrainingAlgorithm ta){
			    	  System.out.println(ta.getCurrentScore());
			        return (ta.getCycle() > 1);
			      }
			    };
			    
			    /*
			     * optimize the dp matrix to reflect the training set in db using a null model
			     * weight of 1.0 and the Stopping criteria defined above.
			     */
			    
			 //Database to hold the training set
			      BufferedReader br = new BufferedReader(new FileReader(inputFasta));
			      //get a SequenceDB of all sequences in the file
			      SequenceDB db = new HashSequenceDB();
   		          SymbolTokenization toke = AlphabetManager.alphabetForName("DNA").getTokenization("token");
   		          SequenceIterator seqi = RichSequence.IOTools.readFasta(br, toke,null);
   			      while (seqi.hasNext()) {
   				      db.addSequence(seqi.nextSequence());
   			      }
			    
			    
			 bwt.train(db,1.0,stopper);
			 

			 
			
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
			oos.writeObject(dp.getModel());
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
			dp=DPFactory.DEFAULT.createDP((MarkovModel)si.readObject());
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
		} catch (BioException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} 
		
	}
	

}
