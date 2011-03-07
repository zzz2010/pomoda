/**
 * @author zhizhuo zhang
 * zzz2010@gmail.com
 */
import java.util.*;
public interface ISearchEngine {
	
 void	build_index(String inputfile);
 LinkedList<Integer> searchPattern(String pattern, int mismatch);
 LinkedList<FastaLocation> Int2Location(LinkedList<Integer> input);
 String getSite(int location,int len);
 int getTotalLength();
 int getSeqNum();
 

}
