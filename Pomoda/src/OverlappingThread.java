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
		this.list1=list1;
		this.list2=list2;
		this.windowsize=windowsize;
	}
	
	public void run() {
		Iterator<Integer> iter1=list1.iterator();
		Iterator<Integer> iter2=list2.iterator();
		
		
	}
	
	//return the overlapping point
	public LinkedList<Integer> getResult()
	{
		return result;
	}
	

}
