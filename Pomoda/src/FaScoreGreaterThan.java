import java.util.Comparator;


public class FaScoreGreaterThan implements Comparator<FastaLocation> {

	@Override
	public int compare(FastaLocation arg0, FastaLocation arg1) {
		// TODO Auto-generated method stub
		if(arg0.Score<arg1.Score)
			return 1;
		else if(arg0.Score==arg1.Score)
			return 0;
		else
			return -1;
	}

}
