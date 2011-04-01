import java.util.ArrayList;

import edu.stanford.rsl.jpop.AdditiveFunctionAssembler;
import edu.stanford.rsl.jpop.FunctionOptimizer;
import edu.stanford.rsl.jpop.GradientOptimizableFunction;
import edu.stanford.rsl.jpop.SimpleFunctionController;
import edu.stanford.rsl.jpop.FunctionOptimizer.OptimizationMode;


public class LogitFunction implements GradientOptimizableFunction {
	ArrayList<Double> Ez;
	ArrayList<Double> PWMscore;

	ArrayList<Double> NBscore;
	ArrayList<Double> MNscore;
	
	
	public LogitFunction(ArrayList<Double> ez, ArrayList<Double> pWMscore,
			ArrayList<Double> nBscore, ArrayList<Double> mNscore) {
		super();
		Ez = ez;
		PWMscore = pWMscore;
		NBscore = nBscore;
		MNscore = mNscore;
	}
	public double[] run(double[] preLogitpara)
	{
		

		//compute FG model first
		FunctionOptimizer opt=new FunctionOptimizer(2);
		opt.setInitialX(preLogitpara);
		opt.setOptimizationMode(OptimizationMode.Gradient);
		opt.setFunctionController(new SimpleFunctionController());
//		opt.setFunctionController(new ParallelFunctionController());
		opt.setFunctionAssembler(new AdditiveFunctionAssembler());
		opt.optimizeFunction(this);
		double[] x=opt.getOptimum();
		if(x[0]<0.00001||Double.isNaN(x[0]))
			x[0]=0.00001;
		if(x[1]<0.00001||Double.isNaN(x[1]))
			x[1]=0.00001;
		return x;
		
		
	}
	@Override
	public void setNumberOfProcessingBlocks(int number) {
		// TODO Auto-generated method stub

	}

	@Override
	public int getNumberOfProcessingBlocks() {

		return 1;
	}

	@Override
	public double evaluate(double[] x, int block) {
		double sum=0;
//		if(Double.isNaN(x[0])||Double.isNaN(x[1]))
//			return 1e+300;
		if(x[0]<0.00001||Double.isNaN(x[0]))
			x[0]=0.00001;
		if(x[1]<0.00001||Double.isNaN(x[1]))
			x[1]=0.00001;
		if(x[0]>100)
			x[0]=100;
		if(x[1]>100)
			x[1]=100;
		for (int i = 0; i < Ez.size(); i++) {
			double temp=PWMscore.get(i)+NBscore.get(i)*x[0]+MNscore.get(i)*x[1];
			sum+=-(Ez.get(i)*( temp)-Math.log(1+Math.exp(temp)));
			if(Double.isNaN(sum)||Double.isInfinite(sum))
				return 1e+300;
		}
		
		return sum;
	}

	@Override
	public double[] gradient(double[] x, int block) {
		// TODO Auto-generated method stub
		if(x[0]<0.00001||Double.isNaN(x[0]))
			x[0]=0.00001;
		if(x[1]<0.00001||Double.isNaN(x[1]))
			x[1]=0.00001;
		if(x[0]>100)
			x[0]=100;
		if(x[1]>100)
			x[1]=100;
		double[] sum=new double[2];
		for (int i = 0; i < Ez.size(); i++) {
			double temp= NegBinFunction.plogis(PWMscore.get(i)+NBscore.get(i)*x[0]+MNscore.get(i)*x[1]);
			sum[0]+=-(Ez.get(i)-( temp))*NBscore.get(i);
			sum[1]+=-(Ez.get(i)-( temp))*MNscore.get(i);
		}

		return sum;
	}

}
