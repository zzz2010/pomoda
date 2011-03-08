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

/**
 * @author zhizhuo zhang
 * zzz2010@gmail.com
 */
public class common {
	
	
	static double DoubleMinNormal=0.00000000000001;
	
	
    public  static List<PWM> LoadPWMFromFile(String file, int topk)
    {
        List<PWM> retlist = new List<PWM>();
        if (file.indexOf(".pomoda") != -1)
        	retlist = TransfacHandler(file);
        else if (file.indexOf(".traw") != -1)
            retlist = TrawlerHandler(file);
        else if (file.indexOf(".trans") != -1)
            retlist = TransfacHandler(file);
        else if (file.indexOf(".admout") != -1)
            retlist = AmadeusHandler(file);
        else if (file.indexOf(".wee") != -1)
            retlist = WeederHandler(file);

       if(retlist.Count>topk)
       retlist.RemoveRange(topk,retlist.Count-topk);
       return retlist;
       
    }
    
	static List<PWM> WeederHandler(string file)
    {
        List<PWM> retlist = new List<PWM>();
        string line = "";
        StreamReader sr = new StreamReader(file);
        while ((line = sr.ReadLine()) != null)
        {
            string nname = "Motif_" + retlist.Count.ToString();
            if (line.IndexOf("All Occurrences") != -1)
            {
               
                sr.ReadLine();
                List<List<float>> temlist = new List<List<float>>();
                while ((line = sr.ReadLine()) != null)
                {
                    string[] comps = line.Split(new string[] { " " },StringSplitOptions.RemoveEmptyEntries);
                    if (comps.Length < 4)
                        break;
                    List<float> col = new List<float>();
                    float sum = 0;
                    for (int i = 0; i < 4; i++)
                    {
                        float tt = float.Parse(comps[i + 2]);
                        sum += tt;
                        col.Add(tt);
                    }
                    for (int i = 0; i < 4; i++)
                        col[i] /= sum;
                    temlist.Add(col);
                }

                int w = temlist.Count;
                if (w < 5)
                    continue;
                PWM candidate = new PWM(w);
                candidate.Name = nname;
                for (int i = 0; i < w; i++)
                {
                    for (int j = 0; j < 4; j++)
                        candidate.matrix[i, j] = temlist[i][j];
                }
                retlist.Add(candidate);

            }
        }
        sr.Close();

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
			         String nname = line.substring(1);
			       
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
			                 Double tt = Double.parseDouble(comps[i]);
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
			         if (line == null)
			             break;

			         int w = dists.size();
			         if (w < 5)
			             continue;
			         PWM candidate = new PWM(dists.toArray(new Distribution[1]));
			         candidate.Name = nname;
			  
			      
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
