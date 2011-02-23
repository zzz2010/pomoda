import java.util.Collections;
import java.util.List;


public class SortingThread extends Thread {
	
	public SortingThread(List<Integer> sortedlist) {
		super();
		this.sortedlist = sortedlist;
	}

	List<Integer> sortedlist;
	
	public void run()
	{
		Collections.sort(sortedlist); 
	}
	
	public List<Integer> getResult()
	{
		return sortedlist;
	}
	

}
