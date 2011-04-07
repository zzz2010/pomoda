import java.util.Collections;
import java.util.LinkedList;
import java.util.List;


public class SortingThread extends Thread {
	
	public SortingThread(List<Integer> sortedlist) {
		super();
		this.sortedlist = sortedlist;
	}

	List<Integer> sortedlist;
	
	@Override
	public void run()
	{
		Collections.sort(sortedlist); 
		LinkedList<Integer> nodulist=new LinkedList<Integer>();
		
		//remove duplicate position
		int lastpos=Integer.MIN_VALUE;
		for (int i = 0; i < sortedlist.size(); i++) {
			int currpos=sortedlist.get(i);
			if(lastpos<currpos-5)
			{
				nodulist.add(currpos);
			}
			lastpos=currpos;
		}
		this.sortedlist=nodulist;
	}
	
	public List<Integer> getResult()
	{
		return sortedlist;
	}
	

}
