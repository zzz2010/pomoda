import java.util.Comparator;
import java.util.Map;

/**
 * @author zhizhuo zhang
 * zzz2010@gmail.com
 */
public class common {
	
	static String Hash2ACGT(int hash,int len)
	{
		String ACGT="ACGT";
		int k;
		StringBuffer pa=new StringBuffer("");
			for(k=0;k<len;k++)
			{
				pa.append( ACGT.charAt(hash%4));
				hash>>=2;
			}	
			return pa.toString();
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
	static String getReverseCompletementString(String Tag)
	{
			String temp="";
			for(int i=0;i<Tag.length();i++)
				temp=temp.concat(String.valueOf((reverseC(Tag.charAt(Tag.length()-1-i)))));
			//temp.push_back(reverseC(Tag[i]));
			return temp;
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

}

class ValueComparator implements Comparator<Map.Entry<Integer,Double>> {

	
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
