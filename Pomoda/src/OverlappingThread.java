import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;


public class OverlappingThread extends Thread {
	
	List<Integer> list1;
	List<Integer> list2;
	int windowsize=20;
	LinkedList<Integer> result;
	//assume the list is sorted
	public OverlappingThread(List<Integer> list1,List<Integer> list2 , int windowsize)
	{	
		super();
		this.list1=list1;
		this.list2=list2;
		this.windowsize=windowsize;
		result=new LinkedList<Integer>();
	}
	
	@Override
	public void run() {
		Iterator<Integer> iter1=list1.iterator();
		Iterator<Integer> iter2=list2.iterator();
		if(list1.size()==0||list2.size()==0)
			return;
		Integer pos1 = iter1.next();
		Integer pos2 = iter2.next();
		while(iter1.hasNext()&&iter2.hasNext())
		{
		
			if(Math.abs(pos1-pos2)<=windowsize)
			{
				result.add((pos1+pos2)/2);
				pos1=iter1.next();
				pos2=iter2.next();
			
			}
			else if(pos1<(pos2-windowsize))
				pos1=iter1.next();
			else if(pos2<(pos1-windowsize))
				pos2=iter2.next();
		}
		
	}
	
	//return the overlapping point
	public LinkedList<Integer> getResult()
	{
		return result;
	}
	

}
