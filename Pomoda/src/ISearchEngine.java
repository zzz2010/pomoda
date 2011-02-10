import java.util.*;
public interface ISearchEngine {
	
 void	build_index(String inputfile);
 LinkedList<Integer> searchPattern(String pattern, int mismatch);
 String getSite(int location,int len);

}
