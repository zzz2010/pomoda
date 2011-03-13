/**
 * @author zhizhuo zhang
 * zzz2010@gmail.com
 */
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.LinkedList;


public class HashEngine implements ISearchEngine {

	public HashEngine(int hashlen) {
		super();
		this.Hashlen = hashlen;
		
		int addrange=1<<(hashlen*2);
		HashIndex= new ArrayList< LinkedList<Integer>>(addrange);
		for(int i=0;i<addrange;i++)
			HashIndex.add(new LinkedList<Integer>());
		
		BGProb=new double[4];
		CDProb=new double[4];
	}

	private int Hashlen;
	private ArrayList< LinkedList<Integer>> HashIndex;
	public String CharText; 
	public int SeqNum;
	public int TotalLen;
	public ArrayList<Integer> accSeqLen;
	public double[] BGProb;
	public double[] CDProb;
	public int forwardCount; //the first n entry in the searchPattern is forward.
	
	
	

	@Override
	public void build_index(String inputfile) {
		int maxHash=1<<(2*Hashlen);
		CharText="";
		int i,j,k;
		SeqNum=0;
		TotalLen=0;
		accSeqLen=new ArrayList<Integer>(10000);
		 try {
			BufferedReader input =  new BufferedReader(new FileReader(inputfile));
			String line = null;
	        while (( line = input.readLine()) != null){
	            if(line.startsWith(">"))
	            {
	            	SeqNum++;
	            	accSeqLen.add(TotalLen);
	            }
	            else
	            {
	            	
	            	
	            	for(i=0;i<line.length()-Hashlen;i++)
	            	{
	            		int hash=common.getHashing(line,i,Hashlen);
	    				if(hash>=0&&hash<(maxHash))
	    				{
	    					HashIndex.get(hash).add(TotalLen+i);
	    					
	    				}
	            	
	            	}
	            	line=line.toUpperCase().replace("N", "");
	            	TotalLen+=line.length();
	            	CharText=CharText.concat(line).concat("N");
	            	TotalLen+=1;
	            	
	            }
	            
	          
	          }
	        accSeqLen.add(TotalLen);
			
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

	private String gPattern;
	LinkedList<Integer> SmartShift(String pattern)
	{
		LinkedList<Integer> ret=new LinkedList<Integer>();
		int len= pattern.length();
		int i,j,k,q;
		int maxaln=0;
		int bestScore=Hashlen;
		int curScore=0;
		LinkedList<Integer> anchors = new LinkedList<Integer>() ;
		//find out the conservativest part
		for(i=0;i<Hashlen;i++)
		{
			if((pattern.charAt(i))=='N'||(pattern.charAt(i))=='n')
			{
				bestScore--;
			}
		}
		for(i=0;i<len-Hashlen;i++)
		{
			if(i==0)
				curScore=bestScore;
			else
			{
				if((pattern.charAt(i))=='N'||(pattern.charAt(i))=='n')
				{
					curScore++; //one go out
				}
				if((pattern.charAt(i+Hashlen-1))=='N'||(pattern.charAt(i+Hashlen-1))=='n')
				{
					curScore--; //one go out
				}
				if(curScore>bestScore)
				{
					maxaln=i+1;
					bestScore=curScore;
				}

			}
		}
		
		String temp=pattern.substring(maxaln, maxaln+Math.min(len-maxaln,Hashlen));
		
		
		anchors.addAll(searchPatternNRC(temp,0,0));
		
		LinkedList<Integer> mask = new LinkedList<Integer>();
		
		for(i=0;i<len;i++)
		{

			if(i>=maxaln&&i<(maxaln+Hashlen))
				continue;
			if((pattern.charAt(i))=='N'||(pattern.charAt(i))=='n')
				continue;
			
			mask.add(i);
		}

		int size=anchors.size();
		Iterator<Integer> iter=anchors.iterator();
		for(i=0;i<size;i++)
		{
			int pos=iter.next();
			boolean flag=true;
			Iterator<Integer> iter2=mask.iterator();
			for(j=0;j<mask.size();j++)
			{
				int markiter=iter2.next();
				int pos2=pos-(maxaln-markiter);
				if(pos2<0||common.acgt(CharText.charAt(pos2))!=common.acgt(pattern.charAt(markiter)))
				{
					flag=false;
					break;
				}
			}

			if(flag)
				ret.add(pos-maxaln);
		}
		
		return ret;
	
	}
	
	LinkedList<String> NCloseSet(String pattern)
	{
		String ACGT="ACGT";
		int i,j,k;
		LinkedList<Integer> pos = new LinkedList<Integer>();
		LinkedList<String> retSet = new LinkedList<String>();
		for(i=0;i<pattern.length();i++)
		{
			if(pattern.charAt(i)=='N')
			{
				pos.add(i);
			}
		}

		if(pos.size()==0)
		{
			retSet.add(pattern);
			return retSet;
		}
		int loop=1<<(2*pos.size());

		for(i=0;i<loop;i++)
		{
			String temp=pattern;
			Iterator<Integer> iterator = pos.iterator(); 
			for(j=0;j<pos.size();j++)
			{
				
				int code=(i>>(j*2))%4;
				temp=common.replaceCharAt(temp, iterator.next(), ACGT.charAt(code));
				

			}
			retSet.add(temp);
		}

		

		return retSet;
	
	}
	
	LinkedList<Integer> searchPatternNRC( String pattern, int pos, int mismatches)
	{
		LinkedList<Integer> ret=new LinkedList<Integer>();
		int len=pattern.length();
		if(mismatches>0)
		{
			String ACGT="ACGT";
			int i,j;		
			for(i=pos;i<len;i++)
			{		
				if(pattern.charAt(i)=='N')
					continue;
				for(j=0;j<3;i++)
				{
					String temp=new String(pattern);
					temp=common.replaceCharAt(temp, i, ACGT.charAt((common.acgt(temp.charAt(i))+j+1)%4));
					
					ret.addAll(searchPatternNRC(temp,pos+1,mismatches-1));
					
				}
			}
			ret.addAll(searchPatternNRC(pattern,pos+1,0));
		}
		else //0 mismatch
		{
			//handle 'N'
			if(pattern.length()>Hashlen)
			{
				ret.addAll(SmartShift(pattern));

				return ret;
			}
			LinkedList<String> searchSet=NCloseSet(pattern);
			if(searchSet.size()==1)
			{
				int left,range;
				int[] hash=new int[3];
				int i;
				for(i=0;i<3;i++)
					hash[i]=-1;	
				String pa=searchSet.getFirst();
				if(pa.length()>3*Hashlen)
					pa=pa.substring(0,3*Hashlen);
				
				for(i=0;i<3;i++)
				{
					
					if(pa.length()<=(i+1)*Hashlen)
						break;
					hash[i]=common.getHashing(pa,i*Hashlen,Hashlen);
					
				}	
				int restlen=pa.length()%Hashlen;
				if(restlen==0)
					restlen=Hashlen;
				hash[i]=common.getHashing(pa,i*Hashlen,restlen);
				left=(Hashlen-pa.length()%Hashlen)%Hashlen;

				range=1<<(2*left);
				hash[i]<<=(2*left);
				ret.addAll( JoinMerge(range,hash[0],hash[1],hash[2]));

				return ret;
			}
			int i;
			Iterator<String> iterator=searchSet.iterator();
			for(i=0;i<searchSet.size();i++)
			{
				String temp=iterator.next();
			
				ret.addAll(searchPatternNRC(temp,pos+1,0));
				
			}
		}
		return ret;
		
	}
		
		
		LinkedList<Integer> JoinMerge(int Range, int hash1, int hash2, int hash3)
		{
			LinkedList<Integer> ret=new LinkedList<Integer>();
			int i;
			if(hash3!=-1)
			{
				for(i=0;i<Range;i++)
				{
					
					int pos1,pos2,pos3;
					Iterator<Integer> iter1=HashIndex.get(hash1).iterator();
					Iterator<Integer> iter2=HashIndex.get(hash2).iterator();
					Iterator<Integer> iter3=HashIndex.get(hash3+i).iterator();
					pos1=iter1.next();
					pos2=iter2.next();
					pos3=iter3.next();
					while(iter1.hasNext()&&iter2.hasNext()&&iter3.hasNext())
					{


						if(pos1==(pos2-Hashlen)&&pos2==(pos3-Hashlen))
						{
							ret.add(pos1);
							pos1=iter1.next();
							pos2=iter2.next();
							pos3=iter3.next();
						}
						if(pos1<(pos2-Hashlen))
							pos1=iter1.next();
						if(pos2<(pos3-Hashlen))
							pos2=iter2.next();
						if(pos3<(pos2+Hashlen))
							pos3=iter3.next();
					}
				}
			}
			else if(hash2!=-1)
			{
				for(i=0;i<Range;i++)
				{
					
					int pos1,pos2;
					Iterator<Integer> iter1=HashIndex.get(hash1).iterator();
					Iterator<Integer> iter2=HashIndex.get(hash2+i).iterator();
					pos1=iter1.next();
					pos2=iter2.next();
					while(iter1.hasNext()&&iter2.hasNext())
					{
					
						if(pos1==(pos2-Hashlen))
						{
							ret.add(pos1);
							pos1=iter1.next();
							pos2=iter2.next();
						
						}
						if(pos1<(pos2-Hashlen))
							pos1=iter1.next();
						if(pos2<(pos1+Hashlen))
							pos2=iter2.next();
					}
				}
			}
			else
			{
				//if(hash1>65536||hash1<0)
				//	return;
				for(i=0;i<Range;i++)
				{
					
				
					
					int pos1;
					Iterator<Integer> iter1=HashIndex.get(hash1+i).iterator();
					if(!iter1.hasNext())
						continue;
					pos1=iter1.next();
					while(iter1.hasNext())
					{
						
						ret.add(pos1);
						pos1=iter1.next();
						
					}
				}
			}
		
			return ret;
		}
	@Override
	public LinkedList<Integer> searchPattern(String pattern, int mismatch) {
		LinkedList<Integer> ret=new LinkedList<Integer>();
		gPattern=pattern;
		ret.addAll( searchPatternNRC(pattern,0,mismatch));
		
			//reverse searching
		String temp=common.getReverseCompletementString(pattern);
		forwardCount=ret.size();
		if(temp.equals(pattern))
			return ret;
		gPattern=temp;
		
		ret.addAll( searchPatternNRC(temp,0,mismatch));
		return ret;
	}
	
	public LinkedList<FastaLocation> searchPattern(PWM motif, double thresh) {
		LinkedList<FastaLocation> ret=new LinkedList<FastaLocation>();
		//find most converse Hashlen part
		String Consensus=motif.Consensus(true);
		double[] columnCons=new double[Consensus.length()];
		int beststart=-1;
		double minCon_score=Double.MAX_VALUE;
		double Con_score=0;
		for (int i = 0; i < Consensus.length(); i++) {
			int symid=common.acgt(Consensus.charAt(i));
			if(symid<4)
			{
				columnCons[i]=1/motif.m_matrix[motif.head+i][symid];
				
			}
			else
				columnCons[i]=4;
			
			Con_score+=columnCons[i];
			if(i>Hashlen-2)
			{
				if(Con_score<minCon_score)
				{
					beststart=i-Hashlen+1;
					minCon_score=Con_score;
				}
				if(i-Hashlen+1>=0)
					Con_score-=columnCons[i-Hashlen+1];
			}
			
		}
		
		Integer[] seqlenlist=accSeqLen.toArray(new Integer[1]);
		//forward strand
		String pattern=Consensus.substring(beststart,beststart+Hashlen);
		LinkedList<Integer> poslist=searchPatternNRC(pattern,0,0);
		Iterator<Integer> iter=poslist.iterator();
		while(iter.hasNext())
		{
			int pos=iter.next();
			String site=this.getSite(pos-beststart, Consensus.length());
			double score=motif.scoreWeightMatrix(site);
			if(score>thresh)
			{
				pos=pos-beststart;
				int seqNum=0;
				//use binary search
				seqNum=Arrays.binarySearch(seqlenlist, pos);
				if(seqNum<0)
					seqNum=-seqNum-2;
				int seqLen=accSeqLen.get(seqNum+1)-accSeqLen.get(seqNum);
				FastaLocation fapos=new FastaLocation(pos, seqNum, pos-accSeqLen.get(seqNum), seqLen);
				fapos.Score=score;
				ret.add(fapos);
			}
		}
		forwardCount=ret.size();
		
		if(Consensus.equalsIgnoreCase(common.getReverseCompletementString(Consensus)))
			return ret;
		//reverse strand
		pattern=common.getReverseCompletementString(pattern);
		poslist=searchPatternNRC(pattern,0,0);
		iter=poslist.iterator();
		while(iter.hasNext())
		{
			int pos=iter.next();
			String site=this.getSite(pos-(Consensus.length()-beststart-Hashlen), Consensus.length());
			site=common.getReverseCompletementString(site);
			double score=motif.scoreWeightMatrix(site);
			if(score>thresh)
			{
				pos=pos-(Consensus.length()-beststart-Hashlen);
				int seqNum=0;
				
				//use binary search
				seqNum=Arrays.binarySearch(seqlenlist, pos);
				if(seqNum<0)
					seqNum=-seqNum-2;
				int seqLen=accSeqLen.get(seqNum+1)-accSeqLen.get(seqNum);
				FastaLocation fapos=new FastaLocation(pos, seqNum, pos-accSeqLen.get(seqNum), seqLen);
				fapos.Score=score;
				fapos.ReverseStrand=true;
			
				ret.add(fapos);
			}
		}
		
		
		return ret;
	}
	
	@Override
	public LinkedList<FastaLocation> Int2Location(LinkedList<Integer> input)
	{
		LinkedList<FastaLocation> ret=new LinkedList<FastaLocation>();
		Iterator<Integer> iter=input.iterator();
		
		while(iter.hasNext())
		{
			int seqNum=0;
			int pos=iter.next();
			while(pos>accSeqLen.get(seqNum+1))
				seqNum++;
			int seqLen=accSeqLen.get(seqNum+1)-accSeqLen.get(seqNum);
			
			FastaLocation a=new FastaLocation(pos, seqNum, pos-accSeqLen.get(seqNum), seqLen);

			ret.add(a);
			
		}
		
		
		return ret;
	}

	@Override
	public String getSite(int location, int len) {
		  if(location<0)
		  {
			  String ret="";
			  for (location = 0; location < 0; location++) {
				ret+="N";
			}
			  ret+=CharText.substring(location, location+len);
		  }
		if((location+len)>=CharText.length())
			return "X";
		return CharText.substring(location, location+len);
	}
	
	


	@Override
	public int getTotalLength() {
		// TODO Auto-generated method stub
		return TotalLen;
	}

	@Override
	public int getSeqNum() {
		// TODO Auto-generated method stub
		return SeqNum;
	}

}
