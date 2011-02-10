
public class common {
	
	
	
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
	static String reverseString(String Tag)
	{
			String temp="";
			for(int i=0;i<Tag.length();i++)
			temp.concat(String.valueOf((reverseC(Tag.charAt(Tag.length()-1-i)))));
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
