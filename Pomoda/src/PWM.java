/**
 * @author zhizhuo zhang
 * zzz2010@gmail.com
 */

import java.io.BufferedReader;
import java.io.IOException;
import java.io.StringReader;
import java.util.*;
import java.util.Map.Entry;

import java.util.Iterator;

import org.apache.commons.lang.ArrayUtils;
import org.biojava.bio.dist.Count;
import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dist.DistributionFactory;
import org.biojava.bio.dist.DistributionTools;
import org.biojava.bio.dist.IndexedCount;
import org.biojava.bio.dist.SimpleDistribution;
import org.biojava.bio.dp.DP;
import org.biojava.bio.dp.ScoreType;
import org.biojava.bio.dp.SimpleWeightMatrix;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.*;
import org.biojava.utils.ChangeVetoException;

import umontreal.iro.lecuyer.probdist.NegativeBinomialDist;
import umontreal.iro.lecuyer.probdist.NormalDist;
import umontreal.iro.lecuyer.util.Num;

public class PWM extends SimpleWeightMatrix {

	static double[] bg_prob=new double[]{0.25,0.25,0.25,0.25};
	double[][]  m_matrix;
	double[][]  log_matrix;
	public int head;
	public String Name="";
	LinkedList<Integer> matchsite=new LinkedList<Integer>();
	boolean strand_en=false;
	boolean pos_en=false;
	boolean peakrank_en=false;
	public double strand_plus_prior=0.5;
	public ArrayList<Double>peakrank_prior=new ArrayList<Double>();
	public int core_motiflen;
	public int tail;
	public double inst_coverage=1;
	public double inst_FDR=1;
	double infothresh=0.4;
	public double Prior_EZ=0;
	//public ArrayList<Double> debuglist=new ArrayList<Double>(); 
	public double Score;
	public ArrayList<Double>pos_prior;
	public ArrayList<Double>Dnase_prob=null;
	public NegativeBinomialDist DnaseBG=null;
	public NegativeBinomialDist DnaseFG=null;
	
	public PWM(Distribution[] arg0) throws IllegalAlphabetException {
		super(arg0);
		// TODO Auto-generated constructor stub
		
		head=0;
		tail=0;
		core_motiflen=arg0.length;
		pos_prior=new ArrayList<Double>();
		Score=Double.MIN_VALUE;
		SymbolList sla;
		try {
			sla = DNATools.createDNA("ACGT");
			m_matrix=new double[this.columns()][4];
			log_matrix=new double[this.columns()][4];
			for (int i = 0; i < this.columns(); i++) {
				Distribution di=this.getColumn(i);
				for (int j = 1; j<= 4; j++)
				{  
					double weight=di.getWeight(sla.symbolAt(j));
			        m_matrix[i][j-1]=weight;
			        log_matrix[i][j-1]=Math.log(weight);
						
				}
				
		
				
			}
		} catch (IllegalSymbolException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	
	
//	public double calcLogDnaseProb(Double[] data, int start)
//	{
//		double sum=0;
//		for (int i = start; i < start+Dnase_prob.size(); i++) {
//			sum+=data[i];
//		}
//		double p1=sum*Math.log(1-DnaseFG.getP())+DnaseFG.getGamma()*Math.log(DnaseFG.getP())+Num.lnGamma(sum+DnaseFG.getGamma())-Num.lnGamma(DnaseFG.getGamma());
//		p1-=sum*Math.log(1-DnaseBG.getP())+DnaseBG.getGamma()*Math.log(DnaseBG.getP())+Num.lnGamma(sum+DnaseBG.getGamma())-Num.lnGamma(DnaseBG.getGamma());
//		double p2=0;
//		double[] data_1=common.Normalize(ArrayUtils.toPrimitive( data));
//		for (int i = start; i <start+Dnase_prob.size(); i++) {
//			p2+=Math.log((Dnase_prob.get(i-start)*Dnase_prob.size()))*data_1[i];
//		}
//		if(Double.isInfinite(p2+p1)||Double.isNaN(p2+p1))
//			return 0;
//		
//		return p1;
//	}
	
	double NegBinConfidence=0.0001;
	double MultiNomConfidence=0.0001;
	public double calcLogDnaseProb(Double[] data, int start)
	{

		double p1= calcLogDnaseNegBinProb(data, start);
		double p2=calcLogDnaseMultiNomProb(data, start);

	
		
		if(Double.isInfinite(p2+p1)||Double.isNaN(p2+p1))
			return 0;
		
		return p1+p2;
	}
	
	public double calcLogDnaseNegBinProb(Double[] data, int start)
	{
		double sum=0;
		Double min = (Double) Collections.min(Arrays.asList(data));
		for (int i = start; i < start+Dnase_prob.size(); i++) {
			sum+=data[i]-min;
		}
		double p1=sum*Math.log(1-DnaseFG.getP())+DnaseFG.getGamma()*Math.log(DnaseFG.getP())+Num.lnGamma(sum+DnaseFG.getGamma())-Num.lnGamma(DnaseFG.getGamma());
		p1-=sum*Math.log(1-DnaseBG.getP())+DnaseBG.getGamma()*Math.log(DnaseBG.getP())+Num.lnGamma(sum+DnaseBG.getGamma())-Num.lnGamma(DnaseBG.getGamma());
		
		p1=p1*NegBinConfidence;
		return p1;
	}
	
	
//	double samplesize=1;
//	double Log_Fab_ratio=0;
	public double calcLogDnaseMultiNomProb(Double[] data, int start)
	{
		double p2=0;
//		Double min = (Double) Collections.min(Arrays.asList(data));
		double[] data_1=common.Normalize(ArrayUtils.toPrimitive( data));
		
		Double min = (Double) Collections.min(Arrays.asList(ArrayUtils.toObject(data_1)));
		
		for (int i = start; i <start+Dnase_prob.size(); i++) {
			p2+=Math.log((Dnase_prob.get(i-start)*Dnase_prob.size()))*(data_1[i]-min);//(data[i]-min);//
		}
//		p2=Log_Fab_ratio;
//		int sum=0;
//		for (int i = start; i <start+Dnase_prob.size(); i++) { 
//			p2+=Num.lnGamma(Dnase_prob.get(i-start)*samplesize+data[i]-min);
//			sum+=Dnase_prob.get(i-start)*samplesize+data[i]-min;
//		}
//		p2-=Num.lnGamma(sum);
//		p2+=sum*Math.log(Dnase_prob.size());
	
		p2=p2*MultiNomConfidence;
		return p2-Math.log(Dnase_prob.size());
	}
	
	public static final Distribution[] countMatrix2Distribution(double[][] count_matrix) throws IllegalSymbolException, IllegalAlphabetException
	{

	    //make a Distribution[] of the motif
	    Distribution[] dists =new Distribution[count_matrix.length];
	    FiniteAlphabet alpha = DNATools.getDNA();
	    for (int i = 0; i < dists.length; i++) {
	    	Count c = new IndexedCount(alpha);
	        c.increaseCount(DNATools.a(),count_matrix[i][0]);
	        c.increaseCount(DNATools.c(),count_matrix[i][1]);
	        c.increaseCount(DNATools.g(),count_matrix[i][2]);
	        c.increaseCount(DNATools.t(),count_matrix[i][3]);
	    	dists[i]=DistributionTools.countToDistribution(c);
		}
	    
	    return	dists;
			
	}
	public static final Distribution[] alignment2Distribution(String[] alignments) throws IllegalSymbolException, IllegalAlphabetException
	{
		Map<String, SymbolList> map = new HashMap<String, SymbolList>();
		for (int i = 0; i < alignments.length; i++) {
			map.put(String.valueOf(i),DNATools.createDNA( alignments[i])) ;
		}
		Alignment align = new SimpleAlignment(map);
	    //make a Distribution[] of the motif
	    Distribution[] dists =
	        DistributionTools.distOverAlignment(align, false, common.DoubleMinNormal);
	    
	    
	    
	    return	dists;
			
	}
	
	public void setWeights(int col, double[] weights)
	{
		for (int i = 0; i < weights.length; i++) {
			m_matrix[col][i]=weights[i];
			log_matrix[col][i]=Math.log(weights[i]);
			
		}
		Distribution di=this.getColumn(col);
		SymbolList sla;
		double weight=0;
		if(head>col)
		{
			core_motiflen+=head-col;
			head=col;
		}
		else if(tail>columns()-col-1)
		{
			core_motiflen+=tail-(columns()-col-1);
			tail=columns()-col-1;
		}
		
		try {
			sla =DNATools.createDNA("ACGT");
			for (int i = 1; i <= weights.length; i++) {
				this.getColumn(col).setWeight(sla.symbolAt(i), weights[i-1]);
			}
		} catch (IllegalSymbolException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IndexOutOfBoundsException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
	
	public void setWeight(int col, int symid,double weight)
	{
		m_matrix[col][symid]=weight;
		log_matrix[col][symid]=Math.log(weight);
		Distribution di=this.getColumn(col);
		SymbolList sla;

			try {
				sla =DNATools.createDNA("ACGT");
				//start from 1
				this.getColumn(col).setWeight(sla.symbolAt(symid+1), weight);
			} catch (IllegalSymbolException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			
		
		
	}
	
	
	public double getWeight(int col, int symid)
	{

		return m_matrix[col][symid];
		
	}
	
	public double getLogWeight(int col, int symid)
	{

		return log_matrix[col][symid];
		
	}
	


	public PWM(String[] alignments) throws IllegalAlphabetException, IllegalSymbolException  {

	    this(alignment2Distribution(alignments));
		
	    
	}

	public PWM(double[][] count_matrix)throws IllegalAlphabetException, IllegalSymbolException  {

	    this(countMatrix2Distribution(count_matrix));
		
	}

	public PWM ReverseComplement()
     {
		  PWM rc =null;
		 try {
		 ArrayList<Distribution> dists=new ArrayList<Distribution>();
		 for (int i = 0; i < columns(); i++) {
			 Distribution di;
			
				di = DistributionFactory.DEFAULT.createDistribution(DNATools.getDNA());
		
             Distribution dii=this.getColumn(columns()-i-1);
			 di.setWeight(DNATools.a(),dii.getWeight(DNATools.t()));
				di.setWeight(DNATools.c(), dii.getWeight(DNATools.g()));
				di.setWeight(DNATools.g(), dii.getWeight(DNATools.c()));
				di.setWeight(DNATools.t(), dii.getWeight(DNATools.a()));
			dists.add(di);
		}
			rc= new PWM(dists.toArray(new Distribution[1]));
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
       
	
         rc.Name = this.Name+"_RC";
         return rc;
     }
	
	public String MaximumConsensus()
	{
		StringBuilder sb=new StringBuilder();
		char[] ACGT=new char[]{'A','C','G','T'};
		for (int i = 0; i < core_motiflen; i++) {
			double maxv=0;
			int maxI=0;
			for (int j = 0; j < 4; j++) {
				if(m_matrix[i+head][j]>maxv)
				{
					maxv=m_matrix[i+head][j];
					maxI=j;
				}
			}
			sb.append(ACGT[maxI]);
			
		}
		
		return sb.toString();
	}
	 
	  
	
	public double scoreWeightMatrix( String seq, ScoreType scoreType)
	{
		try {
			return DP.scoreWeightMatrix(this, DNATools.createDNA(seq), scoreType, 1);
		} catch (IllegalSymbolException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return Double.MIN_VALUE;
		}
		
	}
	
	//only consider the core-part, ignore flanking , log score
	public double scoreWeightMatrix( String seq)
	{
		double score=0;
		int len=Math.min(columns()-head-tail, seq.length());
         for (int i = 0; i < len; i++) {
        	 if(common.acgt(seq.charAt(i))>3)
        		 score+=Math.log(0.25);
        	 else
        		 score+=log_matrix[head+i][common.acgt(seq.charAt(i))];
		}
         
         return score;
		
	}
	

	
	public PWM subPWM(int start, int end)
	{
		
		Distribution[] dists=new Distribution[end-start];
		for (int i = start; i < end; i++) {
			dists[i-start]=getColumn(i);
		}
		PWM sub=null;
		try {
			sub = new PWM(dists);
			sub.core_motiflen=core_motiflen-Math.max(0, start-head)-Math.max(0, this.columns()-end-tail);
			sub.head=Math.max(0, head-start);
			sub.tail=Math.max(0, tail-(this.columns()-end));
			sub.Name=sub.Name;
			sub.Score=this.Score;
			sub.pos_prior=(ArrayList<Double>)pos_prior.clone();
		} catch (IllegalAlphabetException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		return sub;
	}
	
	public PWM Clone()
	{
		int start=this.head;
		int end=this.columns()-this.tail;
		Distribution[] dists=new Distribution[end-start];
		for (int i = start; i < end; i++) {
			FiniteAlphabet alpb=(FiniteAlphabet) getColumn(i).getAlphabet();
			dists[i-start]=new SimpleDistribution(alpb)  ;
			Iterator<Symbol> iter= alpb.iterator();
			while(iter.hasNext())
			{
				Symbol sym=iter.next();
				try {
					dists[i-start].setWeight(sym, getColumn(i).getWeight(sym));
				} catch (IllegalSymbolException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				} catch (ChangeVetoException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			
		}
		PWM sub=null;
		try {
			sub = new PWM(dists);
			sub.core_motiflen=core_motiflen-Math.max(0, start-head)-Math.max(0, this.columns()-end-tail);
			sub.head=Math.max(0, head-start);
			sub.tail=Math.max(0, tail-(this.columns()-end));
			sub.Name=sub.Name;
			sub.Score=this.Score;
			sub.strand_en=strand_en;
			sub.pos_en=pos_en;
			sub.peakrank_en=peakrank_en;
			sub.pos_prior=(ArrayList<Double>)pos_prior.clone();
			sub.peakrank_prior=(ArrayList<Double>)peakrank_prior.clone();
			sub.strand_plus_prior=strand_plus_prior;
		} catch (IllegalAlphabetException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		return sub;
	}
	
	
	protected double getBackgroundLogProb(String pattern,BGModel model)
	{
		double logprob=0;
		char[] ACGT=new char[]{'A','C','G','T'};
		StringBuffer sb=new StringBuffer(pattern);
		//fill N with Random
//		Random rand=new Random();
		int Ncount=0;
		for (int i = 0; i < pattern.length(); i++) {
			if(pattern.charAt(i)=='N')
			{
				Ncount++;
			}
			
		}
		logprob=model.Get_LOGPROB(sb.toString())+Ncount*Math.log(4);
		
		
		return logprob;
	}
	
	public double getFDR(double logthresh,BGModel bgmodel)
	{
		int num_sampl=100000;
		ArrayList<Double> scorelist=new ArrayList<Double>(num_sampl);
		String Consensus=this.Consensus(true);
		int count=0;
		//double sumfdr=0;
		//double sumProb=0;
		int fpcount=0;
		while(count<num_sampl)
		{
		 KeyValuePair<Double, String>	sample=bgmodel.generateRandomSequence(Consensus.length());
			//sumfdr+=sample.getKey();
			double score=scoreWeightMatrix(sample.getValue());
			double score2=scoreWeightMatrix(common.getReverseCompletementString(sample.getValue()));
			if(score<score2)
				score=score2;
			score=(score-bgmodel.Get_LOGPROB(sample.getValue()));
			if(score>logthresh)
				fpcount++;
			count++;
			//sumProb+=Math.exp(score);
			
		}
	
		return (double)fpcount/num_sampl;
	}
	public double getThresh(double sampleratio,double FDRthresh,BGModel bgmodel, boolean strict)
	{
		String Consensus=this.Consensus(true);
		//compute the possible path
		double num_path=1;
		double log025=Math.log(0.25);
		int N_num=0;
		for (int i = 0; i < Consensus.length(); i++) {
			if(Consensus.charAt(i)=='N')
			{
				num_path*=4;
				N_num++;
			}
			else 
			{
				int symid=common.acgt(Consensus.charAt(i));
				if(symid<4)
				{
					if(m_matrix[head+i][symid]>0.99)
					num_path*=1;
					else if(m_matrix[head+i][symid]>0.7)
						num_path*=2;
					else
						num_path*=3;
						
				}
				else if(symid>10)
				{
					num_path*=3;
				}
				else
					num_path*=2;
			}
		}
		
		//number sampling
		int num_sampl=100000;
		boolean samplflag=false;
		
		if((num_sampl)>num_path||bgmodel.conditionProb==null)//zzz
		{
			if(bgmodel.conditionProb==null)
			{	
				bgmodel=BGModel.CreateUniform();
				
			}
			LinkedList<Map.Entry<Double,String>> inst=GenerateInstanceFromPWMPQ(sampleratio, FDRthresh, bgmodel);
			
			if(inst.size()>0)
			{
				if(strict)
				//strictly in FDR
					return inst.getFirst().getKey()+log025*N_num+common.DoubleMinNormal;
				else
					return inst.getFirst().getKey()+log025*N_num;
			}
			else
				samplflag=true;
		}
		else
			samplflag=true;
		if(samplflag)
		{
			ArrayList<Double> scorelist=new ArrayList<Double>(num_sampl);
			
			int count=0;
			//double sumfdr=0;
			//double sumProb=0;
			while(count<num_sampl)
			{
			 KeyValuePair<Double, String>	sample=bgmodel.generateRandomSequence(Consensus.length());
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
			//strictly in FDR
			double thresh=scorelist.get((int)Math.floor(num_sampl*(1-FDRthresh)));

			if(strict)
				thresh+=common.DoubleMinNormal;
			return thresh;
			
			
		}
		return Double.MIN_VALUE;
		
		
	}
	
	
	public static PWM parseTransfac(String transfcontent)
	{
		String str;
		BufferedReader reader = new BufferedReader(	new StringReader(transfcontent));
		ArrayList<Distribution> dists=new ArrayList<Distribution>();
		PWM pwm=null;
		String pwmName="";
		double score=0;
		double prior_EZ=0;
		boolean rankbias=false;
		boolean posbias=false;
		try {
			while ((str = reader.readLine()) != null) {
				if(str.startsWith("DE"))
				{
					String[] elms=str.split("\t| ");
					if(elms.length>1)
					pwmName=elms[1];
					if(elms.length>3)
					{
						try
						{
							if(str.indexOf("pos|")>0)
								posbias=true;
							if(str.indexOf("rank|")>0)
								rankbias=true;
							prior_EZ=Double.valueOf(elms[elms.length-1]);
							if(prior_EZ<0||prior_EZ>1)
								prior_EZ=0;
							
						}
						catch(Exception e)
						{
							prior_EZ=0;
						}
						try
						{
							score =Double.valueOf(elms[3]);
							
						}
						catch(Exception e)
						{
							score=0;
						}
					}
				}
				else if(str.startsWith("P"))
				{
					continue;
				}
				else if(str.length()>2)
				{
					
					String[] elms=str.split("\t| ");
					elms=common.trim_StringArray(elms);
					Distribution di= DistributionFactory.DEFAULT.createDistribution(DNATools.getDNA());
					if(elms.length<5)
						break;
					double [] count=new double[4];
					
					for (int i = 1; i <= 4; i++) {

						count[i-1]=Double.parseDouble(elms[i])+common.DoubleMinNormal;
					}
					count=common.Normalize(count);
				
						di.setWeight(DNATools.a(), count[0]);
						di.setWeight(DNATools.c(), count[1]);
						di.setWeight(DNATools.g(), count[2]);
						di.setWeight(DNATools.t(), count[3]);
					dists.add(di);
					//
				}
			}
			pwm=new PWM(dists.toArray(new Distribution[1]));
			pwm.Name=pwmName;
			pwm.Score=score;
			pwm.Prior_EZ=prior_EZ;
			pwm.pos_en=posbias;
			pwm.peakrank_en=rankbias;
			reader.close();
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
	
	
	//sort symid based on the weight: desc order.
	int[] columnOrder(int col)
	{
		int [] order=new int[4];
		TreeMap<Double,Integer> column=new TreeMap<Double,Integer>();
		

			double smallvalue=common.DoubleMinNormal;
			for (int j = 1; j<= 4; j++)
			{  
				double weight=m_matrix[col][j-1];
				column.put(smallvalue*j-weight, j-1);
			}
			Iterator<Integer> iter= column.values().iterator();
			int c=0;
			 while(iter.hasNext())
				 {order[c]=iter.next();
				 c++;
				        
				 }
				 
			
		

		
		return order;
	}
	
	
	// Exact Priority Queue method generate instance from this motif model
	public LinkedList<Map.Entry<Double,String>> GenerateInstanceFromPWMPQ(double sampleratio,double FDRthresh,BGModel bgmodel)
	{
		LinkedList<Map.Entry<Double,String>> InstanceSet=new LinkedList<Map.Entry<Double,String>>();
		TreeMap<Double,String> PQ=new TreeMap<Double,String>();
		double sumFDR=0;
		double sumProb=0;
		char[] ACGT=new char[]{'A','C','G','T'};
		String consensus=this.Consensus(false);
        ArrayList<Integer>  effIndex=new ArrayList<Integer>(consensus.length());
        ArrayList< int[]> orderedCol=new ArrayList< int[]> (consensus.length());
      
        for (int i = 0; i < consensus.length(); i++) {
			if(consensus.charAt(i)!='N')
			{
				effIndex.add(i);
				orderedCol.add(columnOrder(i));
			}
        }
        int effLen=effIndex.size();
    	if(effLen==0)
    		return InstanceSet;
    	int head=effIndex.get(0);
    	int tail=consensus.length()-effIndex.get(effLen-1)-1;
    	//get first one into PQ
    	StringBuffer sb=new StringBuffer("");
    	double logprob=0;
    	double log025=Math.log(0.25);
    	for (int i = 0; i < effLen; i++) {
			int symid=orderedCol.get(i)[0];
			sb.append(ACGT[symid]);
			logprob+=getLogWeight(effIndex.get(i),symid);
			if(i!=effLen-1)
			{
				int gapsize=effIndex.get(i+1)-effIndex.get(i)-1;
				for (int j = 0; j < gapsize; j++) {
					sb.append('N');	
					//logprob+=log025;
				}
			}
		}
    	PQ.put(logprob, sb.toString());
    	Entry<Double,String> topEntry=PQ.pollLastEntry();
    	
    	sumFDR+=Math.exp(getBackgroundLogProb(topEntry.getValue(),bgmodel));
    	//consider reverse complement
    	sumFDR+=Math.exp(getBackgroundLogProb(common.getReverseCompletementString(topEntry.getValue()),bgmodel));
    	sumProb+=Math.exp(logprob);
    	
		while(sumFDR<FDRthresh||InstanceSet.size()==0 )
		{
		
			InstanceSet.push(topEntry);
			if(sumProb>sampleratio)
				break;
			//mutate one mismatch
			for (int z = 0; z < effLen; z++) {
				logprob=topEntry.getKey();
				int changingIndex=effIndex.get(z);
				sb=new StringBuffer(topEntry.getValue());
				
					int orig_symid=common.acgt( sb.charAt(changingIndex-head));
					int symid=-1;
					int[] orders=orderedCol.get(z);
					int i=0;
					for ( i= 0; i < orders.length; i++) {
					     if(orders[i]==orig_symid)
					    	 break;
					}
					if(i==(orders.length-1))
						continue;
					symid=orders[i+1];
					if(getWeight(changingIndex, symid)<0.05)
						continue;
						
											
					sb.setCharAt(changingIndex-head,(ACGT[symid]));
					logprob+=getLogWeight(changingIndex,symid);
					logprob-=getLogWeight(changingIndex,orig_symid);
					
				PQ.put(logprob+z*common.DoubleMinNormal, sb.toString());	 
				
			}
			// get pop out the best one from PQ
			topEntry=PQ.pollLastEntry();
			if(topEntry==null)
				break;
	    	sumFDR+=Math.exp(getBackgroundLogProb(topEntry.getValue(),bgmodel));
	    	sumProb+=Math.exp(topEntry.getKey());						
		}
		
		
		inst_coverage=sumProb; //denote how many percentage of sample above thresh
		inst_FDR=sumFDR;
	   return InstanceSet;	
	}
	
	
	//print matrix
	public void print()
	{
		
		for (int i = 0; i < this.columns(); i++) {
			StringBuffer sb=new StringBuffer("");
			for (int j = 0; j< 4; j++)
			{  
				double weight=m_matrix[i][j];
					sb.append(String.valueOf(weight));
					sb.append('\t');
					
			}
			System.out.println(sb.toString());

			
		}
		
	}
	
	public PWM trim()
	{
		int newhead=head;
		int newtail=tail;
		double [] nullarr=new double[]{0.25,0.25,0.25,0.25};
		for (int i = head; i < this.columns()-tail; i++) {
		    //calculate the information content
		    double info = DistributionTools.bitsOfInformation(this.getColumn(i));
		    if(info<=infothresh)
		    {
		    	newhead=i+1;
		    	setWeights(i,nullarr);
		    }
		    else
		    	break;
		}
		if(newhead>=this.columns())
			return null;
		for (int i = this.columns()-tail-1; i > head; i--) {
		    //calculate the information content
		    double info = DistributionTools.bitsOfInformation(this.getColumn(i));
		    if(info<infothresh)
		    {
		    	newtail=this.columns()-i;
		    	setWeights(i,nullarr);
		    }
		    else
		    	break;
		}
		head=newhead;
		tail=newtail;
		core_motiflen=columns()-head-tail;
		return this;
	}
	
	@Override
	public String toString()
	{
		StringBuffer TransStr=new StringBuffer("");
		String consensus=Consensus(true);
		int newhead=head;
		int newtail=tail;
		if(consensus.length()>8)
		{
			for (int i = head; i < this.columns()-tail; i++) {
			    //calculate the information content
			    double info = DistributionTools.bitsOfInformation(this.getColumn(i));
			    if(info<infothresh)
			    	newhead=i+1;
			    else
			    	break;
			}
			if(newhead>=this.columns())
				return "";
			for (int i = this.columns()-tail-1; i > head; i--) {
			    //calculate the information content
			    double info = DistributionTools.bitsOfInformation(this.getColumn(i));
			    if(info<infothresh)
			    	newtail=this.columns()-i;
			    else
			    	break;
			}
		}
		String priorlist="";
		if(pos_en)
			priorlist+="pos|";
		if(strand_en)
			priorlist+="strand"+strand_plus_prior+"|";
		if(peakrank_en)
			priorlist+="peakrank|";
		TransStr.append("DE\t"+Name+"\t"+consensus.substring(newhead-head,consensus.length()-(newtail-tail))+"\t"+String.valueOf(this.Score)+"\t"+priorlist+"\t"+(Prior_EZ)+"\n");
		TransStr.append("PO\tA\tC\tG\tT\n");
		for (int i = newhead; i < this.columns()-newtail; i++) {

			TransStr.append(i-newhead+1);
			TransStr.append('\t');
			for (int j = 0; j< 4; j++)
			{  
				double weight=m_matrix[i][j];
				TransStr.append(String.valueOf(weight));
				TransStr.append('\t');
					
			}
			TransStr.append(consensus.charAt(i-head)+"\n");
	
		}
		TransStr.append("XX\n");
		
		
		return TransStr.toString();
	}
	
	
	public String Consensus(boolean trim)
	{
		String consensus="";
		String  ACGT="ACGT";
		for (int i = 0; i < this.columns(); i++) {
			StringBuffer sb=new StringBuffer("");
			for (int j = 1; j<= 4; j++)
			{  
				double weight=m_matrix[i][j-1];
				if(weight>(bg_prob[j-1]+0.01)&&weight!=0.25)//side effect control extending length
					sb.append(ACGT.charAt(j-1));
					
			}
			String consensus_pattern=sb.toString();
			if(consensus_pattern.equalsIgnoreCase( "" )) {consensus = consensus + "N"; }
			else if (consensus_pattern.equalsIgnoreCase( "A") ) { consensus = consensus + "A"; }
			else if (consensus_pattern.equalsIgnoreCase("C")  )	{ consensus = consensus + "C"; }
			else if (consensus_pattern.equalsIgnoreCase( "G" ))	{ consensus = consensus + "G"; }
			else if (consensus_pattern.equalsIgnoreCase( "T" ))	{ consensus = consensus + "T"; }
			else if (consensus_pattern.equalsIgnoreCase( "AC" ))	{ consensus = consensus + "M"; }
			else if (consensus_pattern.equalsIgnoreCase( "AG" ))	{ consensus = consensus + "R"; }
			else if (consensus_pattern.equalsIgnoreCase( "AT" ))	{ consensus = consensus + "W"; }
			else if (consensus_pattern.equalsIgnoreCase( "CG" ))	{ consensus = consensus + "S"; }
			else if (consensus_pattern.equalsIgnoreCase( "CT" ))	{ consensus = consensus + "Y"; }
			else if (consensus_pattern.equalsIgnoreCase( "GT" ))	{ consensus = consensus + "K"; }
			else if (consensus_pattern.equalsIgnoreCase( "ACG" ))	{ consensus = consensus + "V"; }
			else if (consensus_pattern.equalsIgnoreCase( "ACT" ))	{ consensus = consensus + "H"; }
			else if (consensus_pattern.equalsIgnoreCase( "AGT" ))	{ consensus = consensus + "D"; }
			else if (consensus_pattern.equalsIgnoreCase( "CGT" ))	{ consensus = consensus + "B"; }
			
		}

		int start,end;
		start=-1;
			end=0;
		for (int i = 0; i < consensus.length(); i++) {
			if(consensus.charAt(i)!='N'& start==-1)
				start=i;
			if(consensus.charAt(consensus.length()-i-1)!='N'& end==0)
				end=consensus.length()-i;
		}

		if(trim)
		{
			if(start>-1)
		   consensus=consensus.substring(start,end);
			else
			{
				 consensus=consensus.substring(0,core_motiflen);
			}
			if(start>-1)
			{
			head=start;
			tail=columns()-end;
			core_motiflen=end-start;
			}
		}

	
		return consensus;
	}
	
	
	

}
