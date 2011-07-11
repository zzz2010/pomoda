import java.util.LinkedList;


public class AUCComputeThread extends Thread {

	public AUCComputeThread(PWMevaluator evaluator, PWM motif,
			LinearEngine bGSearch) {
		super();
		Evaluator = evaluator;
		this.motif = motif;
		this.sequences = bGSearch;
	}
	PWMevaluator Evaluator;
	PWM motif;
	LinearEngine sequences;
	double AUCresult=0;
	@Override
	public void run() {
		// TODO Auto-generated method stub
		super.run();
		AUCresult=Evaluator.calcAUC(motif,sequences);
	}
	public double getResult()
	{
		return AUCresult;
	}
	
	

}
