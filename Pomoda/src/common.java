import java.awt.Color;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.LinkedList;
import java.util.Map;

import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dist.DistributionFactory;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.utils.ChangeVetoException;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.RectangleInsets;

/**
 * @author zhizhuo zhang
 * zzz2010@gmail.com
 */
public class common {
	
	
	static double DoubleMinNormal=0.00000000000001;
	
	public static ArrayList<Double[]> ReadDelimitedFile(String sep, String file)
	{
		ArrayList<Double[]> ret=new ArrayList<Double[]>();
		
		try {
			BufferedReader readbuffer = new BufferedReader(new FileReader(file));
			String strRead;
			while ((strRead=readbuffer.readLine())!=null){
				String splitarray[] = strRead.split(sep);
				Double[] arr=new Double[splitarray.length];
				for (int i = 0; i < arr.length; i++) {
					arr[i]=Double.parseDouble(splitarray[i]);
				}
				ret.add(arr);
				}

				
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		
		return ret;
		
	}
	
    public  static LinkedList<PWM> LoadPWMFromFile(String file)
    {
    	LinkedList<PWM> retlist = new LinkedList<PWM>();
        if (file.endsWith(".pomoda") )
        	retlist = TransfacHandler(file);
        else if (file.endsWith(".traw"))
            retlist = TrawlerHandler(file);
        else if (file.endsWith(".dpwm"))
            retlist = GapPWMHandler(file);
        else if (file.endsWith(".pwm"))
            retlist = TransfacHandler(file);
        else if (file.endsWith(".meme"))
            retlist = MEMEHandler(file);
        else if (file.indexOf("meme.txt") != -1)
            retlist = MEMEHandler(file);
       else if (file.endsWith(".wee"))
            retlist = WeederHandler(file);
       else if (file.endsWith(".cisf") )
           retlist = CisFinderHandler(file);
       else if (file.indexOf("cisfinder") != -1)
           retlist = CisFinderHandler(file);
       else if (file.indexOf("thetafinal") != -1)
           retlist = HMSHandler(file);
       else if (file.indexOf("ChIPMunk") != -1)
           retlist = ChIPMunkHandler(file);

            return retlist;
       
    }
    
    
    
    public static LinkedList<PWM> ChIPMunkHandler(String file) {
    	 LinkedList<PWM> retlist = new LinkedList<PWM>();
         String line = "";
       
         try {
         	  BufferedReader sr = new BufferedReader(new FileReader(new File(file)));
         	  
         		while ((line = sr.readLine()) != null)
     			{  
         	if(line.indexOf("KDIC|")!=-1)
         	{
         	  double[][] m_matrix=null;
         	  int row=0;
         	 int w = 5;
         	 String nname = "ChIPMunk_" + retlist.size();
 			while ((line = sr.readLine()) != null)
 			{
 			   
 			  
 			  String[] comps = line.trim().split("[ \\|\t]+");
 			  if(comps.length<w)
 				  break;
 			  w=comps.length-1;
 			  if(m_matrix==null)
 				  m_matrix=new double[w][4];
 			   for (int i = 0; i < w; i++) {
				m_matrix[i][row]=Double.parseDouble(comps[i+1])+DoubleMinNormal;
			}
 			    row++;
 			    if(row==4)
 			    	break;
 			}
 			  ArrayList<Distribution> dists=new ArrayList<Distribution>();
 			for (int i = 0; i < m_matrix.length; i++) {
				double [] col=m_matrix[i];
 				  common.Normalize(col);
		            Distribution di= DistributionFactory.DEFAULT.createDistribution(DNATools.getDNA());
		             di.setWeight(DNATools.a(), col[0]);
						di.setWeight(DNATools.c(), col[1]);
						di.setWeight(DNATools.g(), col[2]);
						di.setWeight(DNATools.t(), col[3]);
					dists.add(di);
			}
 			  PWM candidate = new PWM(dists.toArray(new Distribution[1]));
 			         candidate.Name = nname;
 			         retlist.add(candidate);

         	}
     			}
 			
 			sr.close();
 		} catch (NumberFormatException e) {
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
 		} catch (IOException e) {
 			// TODO Auto-generated catch block
 			e.printStackTrace();
 		}
         

         return retlist;
	}



	public static LinkedList<PWM> HMSHandler(String file) {
    	 LinkedList<PWM> retlist = new LinkedList<PWM>();
         String line = "";
       
         try {
         	  BufferedReader sr = new BufferedReader(new FileReader(new File(file)));
         	  double[][] m_matrix=null;
         	  int row=0;
         	 int w = 5;
         	
 			while ((line = sr.readLine()) != null)
 			{
 			   
 				 String nname = "HMS_" + retlist.size();
 			  String[] comps = line.trim().split("[ |\t]+");
 			 if(row%4==0)
 			 {
 				 w=5;
 				m_matrix=null;
 				 
 			 }
 			  if(comps.length<w)
 				  break;
 			  w=comps.length;
 			  if(m_matrix==null)
 				  m_matrix=new double[w][4];
 			   for (int i = 0; i < comps.length; i++) {
				m_matrix[i][row%4]=Double.parseDouble(comps[i])+DoubleMinNormal;
			}
 			    row++;
 			    if(row%4==0)
 			    {
 		 			  ArrayList<Distribution> dists=new ArrayList<Distribution>();
 		  			for (int i = 0; i < m_matrix.length; i++) {
 		 				double [] col=m_matrix[i];
 		  				  common.Normalize(col);
 		 		            Distribution di= DistributionFactory.DEFAULT.createDistribution(DNATools.getDNA());
 		 		             di.setWeight(DNATools.a(), col[0]);
 		 						di.setWeight(DNATools.c(), col[1]);
 		 						di.setWeight(DNATools.g(), col[2]);
 		 						di.setWeight(DNATools.t(), col[3]);
 		 					dists.add(di);
 		 			}
 		  			  PWM candidate = new PWM(dists.toArray(new Distribution[1]));
 		  			         candidate.Name = nname;
 		  			         retlist.add(candidate);
 			    }
 			}


 			    
 			
 			sr.close();
 		} catch (NumberFormatException e) {
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
 		} catch (IOException e) {
 			// TODO Auto-generated catch block
 			e.printStackTrace();
 		}
         

         return retlist;
	}



	public static LinkedList<PWM> CisFinderHandler(String file) {
    	 LinkedList<PWM> retlist = new LinkedList<PWM>();
         String line = "";
       
         try {
         	  BufferedReader sr = new BufferedReader(new FileReader(new File(file)));
 			while ((line = sr.readLine()) != null)
 			{
 			    String nname = "Cisfinder_" + retlist.size();
 			    if (line.indexOf(">") != -1)
 			    {
 			       
 			        ArrayList<Distribution> dists=new ArrayList<Distribution>();
 			        while ((line = sr.readLine()) != null)
 			        {
 			            String[] comps = line.trim().split("[ |\t]+");
 			            if (comps.length < 4)
 			                break;
 			           double[] col = new double[4];
 			            float sum = 0;
 			            for (int i = 0; i < 4; i++)
 			            {
 			                double tt = Double.parseDouble(comps[i+1])+DoubleMinNormal;
 			                sum += tt;
 			                col[i]=tt;
 			            }
 			            common.Normalize(col);
 			            Distribution di= DistributionFactory.DEFAULT.createDistribution(DNATools.getDNA());
 			             di.setWeight(DNATools.a(), col[0]);
 							di.setWeight(DNATools.c(), col[1]);
 							di.setWeight(DNATools.g(), col[2]);
 							di.setWeight(DNATools.t(), col[3]);
 						dists.add(di);
 			        }

 			        int w = dists.size();
 			         if (w < 5)
 			             continue;
 			         PWM candidate = new PWM(dists.toArray(new Distribution[1]));
 			         candidate.Name = nname;
 			         retlist.add(candidate);

 			    }
 			}
 			sr.close();
 		} catch (NumberFormatException e) {
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
 		} catch (IOException e) {
 			// TODO Auto-generated catch block
 			e.printStackTrace();
 		}
         

         return retlist;
	}



	public static LinkedList<PWM> MEMEHandler(String file) {
    	 LinkedList<PWM> retlist = new LinkedList<PWM>();
         String line = "";
       
         try {
         	  BufferedReader sr = new BufferedReader(new FileReader(new File(file)));
 			while ((line = sr.readLine()) != null)
 			{
 			    String nname = "MEME_" + retlist.size();
 			    if (line.indexOf("position-specific probability matrix") != -1)
 			    {
 			       
 			        sr.readLine();
 			       sr.readLine();
 			        ArrayList<Distribution> dists=new ArrayList<Distribution>();
 			        while ((line = sr.readLine()) != null)
 			        {
 			            String[] comps = line.trim().split("[ |\t]+");
 			            if (comps.length < 4)
 			                break;
 			           double[] col = new double[4];
 			            float sum = 0;
 			            for (int i = 0; i < 4; i++)
 			            {
 			                double tt = Double.parseDouble(comps[i])+DoubleMinNormal;
 			                sum += tt;
 			                col[i]=tt;
 			            }
 			            common.Normalize(col);
 			            Distribution di= DistributionFactory.DEFAULT.createDistribution(DNATools.getDNA());
 			             di.setWeight(DNATools.a(), col[0]);
 							di.setWeight(DNATools.c(), col[1]);
 							di.setWeight(DNATools.g(), col[2]);
 							di.setWeight(DNATools.t(), col[3]);
 						dists.add(di);
 			        }

 			        int w = dists.size();
 			         if (w < 5)
 			             continue;
 			         PWM candidate = new PWM(dists.toArray(new Distribution[1]));
 			         candidate.Name = nname;
 			         retlist.add(candidate);

 			    }
 			}
 			sr.close();
 		} catch (NumberFormatException e) {
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
 		} catch (IOException e) {
 			// TODO Auto-generated catch block
 			e.printStackTrace();
 		}
         

         return retlist;
	}



	public static double PWM_Divergence(PWM p1, PWM p2)
    {
        double divg = 10;
               
        AlignmentResult result1=Alignment(p1, p2);
        divg =result1.alnScore;
        AlignmentResult result2= Alignment(p1.ReverseComplement(), p2);
       
        if (result2.alnScore < divg)
            divg = result2.alnScore;
        
        return divg;
    }
    
    static AlignmentResult Alignment(PWM refmotif, PWM newmotif)
    {
    	AlignmentResult ret=new AlignmentResult();
    	 float alnScore;
    	int bestoverlap;
                   float bestscore=10;
                   int bestaln=0;
           bestoverlap = 0;
                   int size=newmotif.core_motiflen;
                   int size2=refmotif.core_motiflen;
            alnScore=10;
                   int minsize=Math.min(refmotif.core_motiflen,newmotif.core_motiflen);
                   for(int i=0-size+1;i<size2;i++)
                   {
                           float  curScore=0;
               int pp = 0;
               if (i > 0)
                   pp = i;
               int j;
                           for(j=pp;j-pp<minsize;j++)
                           {
                                   if(j>=size2||(j-i)>=size)
                                           break;
                                   if(j>=0)
                                   {
                       double sump = 0;
                       for (int p = 0; p < 4; p++)
                       {
                           double temp=refmotif.m_matrix[j][p]-newmotif.m_matrix[j-i][p];
                           temp*=(temp);
                           sump += (float)temp;
                       }
                       curScore +=(float) Math.sqrt( sump);
                                   }
                           }
               //size2=minsize
                           int overlap=j;
               ////i>0
               if (i > 0)
                   overlap = j - i;
             
               curScore /= overlap * (float)Math.sqrt(2);
                           if(bestscore>curScore&&(overlap>=7||overlap>=minsize-2))
                           {
                                   bestscore=curScore;
                                   bestaln=i;
                   bestoverlap = overlap;
                           }
                   }
                   alnScore=bestscore;
                   
                   ret.bestaln=(bestaln);
                   ret.alnScore=(alnScore);
                   ret.bestoverlap=(bestoverlap);
                   return ret;
    }
    
    
	static LinkedList<PWM> WeederHandler(String file)
    {
		 LinkedList<PWM> retlist = new LinkedList<PWM>();
        String line = "";
      
        try {
        	  BufferedReader sr = new BufferedReader(new FileReader(new File(file)));
			while ((line = sr.readLine()) != null)
			{
			    String nname = "Weeder_" + retlist.size();
			    if (line.indexOf("All Occurrences") != -1)
			    {
			       
			        sr.readLine();
			        ArrayList<Distribution> dists=new ArrayList<Distribution>();
			        while ((line = sr.readLine()) != null)
			        {
			            String[] comps = line.split("[ |\t]+");
			            if (comps.length < 4)
			                break;
			           double[] col = new double[4];
			            float sum = 0;
			            for (int i = 0; i < 4; i++)
			            {
			                double tt = Double.parseDouble(comps[i + 1])+DoubleMinNormal;
			                sum += tt;
			                col[i]=tt;
			            }
			            common.Normalize(col);
			            Distribution di= DistributionFactory.DEFAULT.createDistribution(DNATools.getDNA());
			             di.setWeight(DNATools.a(), col[0]);
							di.setWeight(DNATools.c(), col[1]);
							di.setWeight(DNATools.g(), col[2]);
							di.setWeight(DNATools.t(), col[3]);
						dists.add(di);
			        }

			        int w = dists.size();
			         if (w < 5)
			             continue;
			         PWM candidate = new PWM(dists.toArray(new Distribution[1]));
			         candidate.Name = nname;
			         retlist.add(candidate);

			    }
			}
			sr.close();
		} catch (NumberFormatException e) {
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
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
        

        return retlist;
    }
	
	
	
	
	 static LinkedList<PWM> TrawlerHandler(String file)
     {
		 LinkedList<PWM> retlist = new LinkedList<PWM>();
         String line = "";
      
      
         try {
        	   BufferedReader sr = new BufferedReader(new FileReader(new File(file)));
			while ((line = sr.readLine()) != null)
			 {
			    
			     while (line.indexOf(">family") != -1)
			     {
			         String nname = "trawler"+line.substring(7);
			       
			         ArrayList<Distribution> dists=new ArrayList<Distribution>();
			         while ((line = sr.readLine()) != null)
			         {
			             String[] comps = line.split("\t");
			             if (comps.length < 4)
			                 break;
			             ArrayList<Double> col = new ArrayList<Double>();
			             float sum = 0;
			             for (int i = 0; i < 4; i++)
			             {
			                 Double tt = Double.parseDouble(comps[i])+DoubleMinNormal;
			                 sum += tt;
			                 col.add(tt);
			             }
			             for (int i = 0; i < 4; i++)
			                 col.set(i, col.get(i)/sum);

			     
			             Distribution di= DistributionFactory.DEFAULT.createDistribution(DNATools.getDNA());
			             di.setWeight(DNATools.a(), col.get(0));
							di.setWeight(DNATools.c(), col.get(1));
							di.setWeight(DNATools.g(), col.get(2));
							di.setWeight(DNATools.t(), col.get(3));
						dists.add(di);
			         }
			         
			         int w = dists.size();
			         if (w < 5)
			             continue;
			         PWM candidate = new PWM(dists.toArray(new Distribution[1]));
			         candidate.Name = nname;
			         retlist.add(candidate);
			         if(line==null)
			        	 break;
			     }
			 }
			 sr.close();
		} catch (NumberFormatException e) {
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
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
        

         return retlist;
     }
	
	 
	 public static LinkedList<PWM> GapPWMHandler(String pwmfile)
		{
			LinkedList<PWM> pwmlist=new LinkedList<PWM>();
			File file = new File(pwmfile);
			try {
				BufferedReader reader = new BufferedReader(new FileReader(file));
				String text = null;
				String content="";
				while ((text = reader.readLine()) != null) {
				    if(text.startsWith("DE"))
				    {
				    	content="";
				    }				    
				    else if(text.startsWith("XXX"))
				    {
				    	pwmlist.add(GapPWM.parseTransfac(content));
				    }
				    else if(text.startsWith("XX")&&!text.startsWith("XXD"))
				    {
				    	pwmlist.add(PWM.parseTransfac(content));
				    }
				    content+=text+'\n';
				}
			} catch (FileNotFoundException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			
			return pwmlist;
		}
	 
	public static LinkedList<PWM> TransfacHandler(String pwmfile)
	{
		LinkedList<PWM> pwmlist=new LinkedList<PWM>();
		File file = new File(pwmfile);
		try {
			BufferedReader reader = new BufferedReader(new FileReader(file));
			String text = null;
			String content="";
			while ((text = reader.readLine()) != null) {
			    if(text.startsWith("DE"))
			    {
			    	content="";
			    }
			    else if(text.startsWith("XX"))
			    {
			    	pwmlist.add(PWM.parseTransfac(content));
			    }
			    content+=text+'\n';
			}
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		return pwmlist;
	}
	
	public static int CountPositive(double[] arr)
	{
		int count=0;
		for (int i = 0; i < arr.length; i++) {
			if(arr[i]>0)
				count++;
		}
		return count;
	}
	
	public static double Zscore(double sumcount,double p, double X)
	{
		return (X-sumcount*p)/Math.sqrt(sumcount*p*(1-p));
		
	}
	
	
	   // return integer nearest to x
	   static long nint(double x) {
	      if (x < 0.0) return (long) Math.ceil(x - 0.5);
	      return (long) Math.floor(x + 0.5);
	   }

	   // return log n!
	   static double logFactorial(int n) {
	      double ans = 0.0;
	      for (int i = 1; i <= n; i++)
	         ans += Math.log(i);
	      return ans;
	   }

	   // return the binomial coefficient n choose k.
	   static long binomial(int n, int k) {
	      return nint(Math.exp(logFactorial(n) - logFactorial(k) - logFactorial(n-k)));
	   }

	static String Hash2ACGT(int hash,int len)
	{
		String ACGT="ACGT";
		int k;
		String kmer="";
			for(k=0;k<len;k++)
			{
				kmer=ACGT.charAt(hash%4)+kmer;
			
				hash>>=2;
			}	
			return kmer;
	};

	

	
	static int getReverseComplementHashing(int hash, int len)
	{
		int ret=0;
		for(int i = 0; i < len; ++i)
		{
			ret<<=2;
			int temp=hash%4;
			hash>>=2;
			ret+=3-temp;
			
			
		}

		return ret;
	}
	
	static char reverseC(char w)
	{
		switch(w)
		{
		case 'a': return 't';
		case 'c': return 'g';
		case 'g': return 'c';
		case 't': return 'a';
		case 'A': return 'T';
		case 'C': return 'G';
		case 'G': return 'C';
		case 'T': return 'A';
		case 'N': return 'N';
		case 'n': return 'N';
		case 'R': return 'Y';
		case 'Y': return 'R';
		case 'W': return 'S';
		case 'S': return 'W';
		case 'V': return 'B';
		case 'B': return 'V';
		case 'K': return 'M';
		case 'M': return 'K';
		case 'H': return 'D';
		case 'D': return 'H';
		
		}
		return 'N';
	}
	static double[] Normalize(double[] Arr)
	{
		double sum=0;
		for (int i = 0; i < Arr.length; i++) {
			sum+=Arr[i];
		}
		if(sum==0)
		{
			for (int i = 0; i < Arr.length; i++) {
				Arr[i]=1.0/ Arr.length;
			}
			return Arr;
		}
		for (int i = 0; i < Arr.length; i++) {
			Arr[i]/=sum;
		}
		return Arr;
		
	}
	
	static String getReverseCompletementString(String Tag)
	{
			String temp="";
			for(int i=0;i<Tag.length();i++)
				temp=temp.concat(String.valueOf((reverseC(Tag.charAt(Tag.length()-1-i)))));
			//temp.push_back(reverseC(Tag[i]));
			return temp;
	}
	static void fill2DArray(double[][] max_loglik_matrix,double val)
	{
		for (int i = 0; i < max_loglik_matrix.length; i++) {
			Arrays.fill(max_loglik_matrix[i], val);
		}
	}
	
	
	static void print2DArray(double[][] max_loglik_matrix)
	{
		for (int i = 0; i < max_loglik_matrix.length; i++) {
			System.out.println( Arrays.toString(max_loglik_matrix[i]));
		}
	}
	
	static String Array2String(Object[] arr , char sep)
	{
		StringBuffer sb=new StringBuffer("");
		for (int i = 0; i < arr.length; i++) {
			if(i>0)
				sb.append(sep);
			sb.append(arr[i]);
		}
		
		return sb.toString();
	}

	public static String replaceCharAt(String s, int pos, char c) {
	   return s.substring(0,pos) + c + s.substring(pos+1);
	}
	
	static  int acgt(char w)
	{
		switch(w)
		{
		case 'a': return 0;
		case 'c': return 1;
		case 'g': return 2;
		case 't': return 3;
		case 'A': return 0;
		case 'C': return 1;
		case 'G': return 2;
		case 'T': return 3;
		case 'R': return 4;
		case 'Y': return 5;
		case 'K': return 6;
		case 'M': return 7;
		case 'S': return 8;
		case 'W': return 9;
		case 'B': return 10;
		case 'D': return 11;
		case 'H': return 12;
		case 'V': return 13;
		case 'N': return 14;
		case 'n': return 14;
		case 'X': return 15;
			//RYKMSWBDHVN
		default: return -1;
		}

	}
	
	static Integer getHashing(String source, int start, int len)
	{
			int i=0;
			int ret=0;
			
			ret = common.acgt( source.charAt(start));
			if(ret<0)
				return -1;
			for (i = start + 1; i < start + len; i++)
			{ 
				ret<<=2;
				int temp=common.acgt(source.charAt(i));
				if(temp<0||temp>3)
					return -1;
				ret+=temp;			
			}
				return ret;

	}

}

class  AlignmentResult
{
	  int bestaln;
      double alnScore;
      int bestoverlap;
}

class ValueComparator implements Comparator<Map.Entry<Integer,Double>> {

	
	  @Override
	public int compare(Map.Entry<Integer,Double> e1, Map.Entry<Integer,Double> e2) {
	        if (e1.getValue() < e2.getValue()){
	            return 1;
	        } else if (e1.getValue() == e2.getValue()) {
	            return 0;
	        } else {
	            return -1;
	        }
	    }


	}
