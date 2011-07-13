import edu.stanford.rsl.jpop.GradientOptimizableFunction;


public class EEM_auxfunction implements GradientOptimizableFunction {

	public EEM_auxfunction(double[] a, double[] b) {
		super();
		A = a;
		B = b;
	}

	double[] A;
	double[] B;
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
		double y=x[0];
		double f=0;
		for (int i = 0; i < 4; i++) {
			f+=Math.exp((y-A[i]-B[i])/A[i]);
		}
		f=(f-1)*(f-1);
//		if(Double.isNaN(f))
//		{	x[0]=0;
//		f=0;
//		}
		
		return f;
	}

	@Override
	public double[] gradient(double[] x, int block) {
		// TODO Auto-generated method stub
		double f=0;
		double y=x[0];
		for (int i = 0; i < 4; i++) {
			f+=Math.exp((y-A[i]-B[i])/A[i]);
		}
//		if(Double.isNaN(f))
//		{
//			x[0]=0;
//			f=0;
//		}
		double df=2*(f-1);
		double sum=0;
		for (int i = 0; i < 4; i++) {
			sum+=1.0/(A[i]*Math.exp((A[i] + B[i] - y)/A[i]));
		}
		df=df*sum;
		return new double[]{df};
	}

}
