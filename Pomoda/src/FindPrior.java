import java.util.ArrayList;


import edu.stanford.rsl.jpop.HessianOptimizableFunction;



public class FindPrior implements HessianOptimizableFunction {

	ArrayList<Double> ProbRatio;
	ArrayList<Double> LogBG;
	
	public FindPrior(ArrayList<Double> probRatio, ArrayList<Double> logBG) {
		super();
		ProbRatio = probRatio;
		LogBG = logBG;
	}

	@Override
	public double[] gradient(double[] x, int block) {
		double df=0;
		double lamda=x[0];
		if(lamda>0.9999)
		{
			lamda=0.9999;
			x[0]=0.9999;
		}
		for (int i = 0; i < ProbRatio.size(); i++) {
			df-=(ProbRatio.get(i)-1)/(1+lamda*(ProbRatio.get(i)-1));
		}
		return new double[]{df};
	}

	@Override
	public void setNumberOfProcessingBlocks(int number) {
		// TODO Auto-generated method stub

	}

	@Override
	public int getNumberOfProcessingBlocks() {
		// TODO Auto-generated method stub
		return 1;
	}

	@Override
	public double evaluate(double[] x, int block) {
		// TODO Auto-generated method stub
		double f=0;
		double lamda=x[0];
		if(lamda>0.9999)
		{
			lamda=0.9999;
			x[0]=0.9999;
		}
		for (int i = 0; i < ProbRatio.size(); i++) {
			f-=Math.log(lamda*ProbRatio.get(i)+1-lamda)+LogBG.get(i);
		}
		return f;
	}

	@Override
	public double[][] hessian(double[] x, int block) {
		double df=0;
		double lamda=x[0];
		if(lamda>0.9999)
		{
			lamda=0.9999;
			x[0]=0.9999;
		}
		for (int i = 0; i < ProbRatio.size(); i++) {
			double temp=(ProbRatio.get(i)-1)/(1+lamda*(ProbRatio.get(i)-1));
			df-=-temp*temp;
		}
		
		double[][] d2x=new double[1][1];
		d2x[0][0]=df;
		return d2x;
	}



}
