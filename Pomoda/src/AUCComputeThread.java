import java.util.LinkedList;


public class AUCComputeThread extends Thread {

	public AUCComputeThread(PWMevaluator evaluator, PWM motif,
			LinkedList<String> sequences) {
		super();
		Evaluator = evaluator;
		this.motif = motif;
		this.sequences = sequences;
	}
	PWMevaluator Evaluator;
	PWM motif;
	LinkedList<String> sequences;
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
