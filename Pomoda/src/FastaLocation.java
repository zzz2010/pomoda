import org.biojava.bio.symbol.PointLocation;


public class FastaLocation extends PointLocation {

	public double Score;
	public String seq="";
	public boolean ReverseStrand=false;
	public FastaLocation(int point, int seqId, int seqPos, int seqLen) {
		super(point);
		this.setSeqId(seqId);
		this.setSeqPos(seqPos);
		this.setSeqLen(seqLen);
	}
	public FastaLocation(int point) {
		super(point);
		// TODO Auto-generated constructor stub
	}
	public void setSeqId(int seqId) {
		this.seqId = seqId;
	}
	public int getSeqId() {
		return seqId;
	}
	public void setSeqPos(int seqPos) {
		this.seqPos = seqPos;
	}
	public int getSeqPos() {
		return seqPos;
	}
	public void setSeqLen(int seqLen) {
		this.seqLen = seqLen;
	}
	public int getSeqLen() {
		return seqLen;
	}
	private int seqId;
	private int seqPos;
	private int seqLen;

}


	