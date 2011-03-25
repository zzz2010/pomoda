import java.util.ArrayList;

import umontreal.iro.lecuyer.util.Num;

import edu.stanford.rsl.jpop.*;
import edu.stanford.rsl.jpop.FunctionOptimizer.OptimizationMode;


public class NegBinFunction implements GradientOptimizableFunction {

	
	ArrayList<Integer> Rl;
	ArrayList<Double> Ez;
	double tao;
	public NegBinFunction(ArrayList<Integer> Rl_stat,ArrayList<Double> Ez_stat, double P)
	{
		Rl=Rl_stat;
		Ez=Ez_stat;
		tao=P;		
	}

	public double run()
	{
		double alpha=0;
		double myM=0;
		double myV=0;
		double sumREz=0;
		double sumEz=0;
		double sumR2Ez=0;
		for (int i = 0; i < Rl.size(); i++) {
			sumREz+=Rl.get(i)*Ez.get(i);
			sumEz+=Ez.get(i);
		//	sumR2Ez+=Rl.get(i)*Rl.get(i)*Ez.get(i);
		}
		myM=(sumREz/sumEz);

		alpha=myM*tao/(1-tao);
		
		int[] iflag=new int[] {0} ;
		double[] x=new double[]{ alpha};
		int[] term=new int[2];
		FunctionOptimizer opt=new FunctionOptimizer(1);
		opt.setInitialX(x);
		opt.setOptimizationMode(OptimizationMode.Gradient);
		opt.optimizeFunction(this);
		x=opt.getOptimum();
		

		return x[0];
		
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
		double alpha=x[0];
//		return alpha*alpha-alpha+1;
		double sum=0;
		for (int i = 0; i < Rl.size(); i++) {
			sum+=Ez.get(i)*(-Num.lnBeta(alpha, Rl.get(i))+alpha*Math.log(tao)+Rl.get(i)*Math.log(1-tao));
		}
	    
		return -sum;
	}
	@Override
	public double[] gradient(double[] x, int block) {
		double alpha=x[0];
//		 return new double[]{2*alpha-1};
		 
			double sum=0;
			for (int i = 0; i < Rl.size(); i++) {
				sum+=Ez.get(i)*(Num.digamma(Rl.get(i)+alpha)-Num.digamma(alpha)+Math.log(tao));
			}
				
			return new double[]{-sum};
	}

}
