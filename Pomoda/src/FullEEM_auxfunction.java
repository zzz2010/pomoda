import edu.stanford.rsl.jpop.GradientOptimizableFunction;


public class FullEEM_auxfunction implements GradientOptimizableFunction {

	public FullEEM_auxfunction(double[] count_matrix, double llr, double[] b) {
		super();
		this.count_matrix = count_matrix;
		this.llr = llr;
		this.b = b;
	}

	double[] count_matrix;
	double llr;
	double[] b;
	
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
		double sumX=0;
		
		for (int i = 0; i < 4; i++) {
			if(x[i]>1)
				return 100000000;
			if(x[i]<0)
				return 100000000;
			double ratio=llr*x[i]/b[i];
			double Z=ratio/(ratio+1);
			f+=-count_matrix[i]*Z*Math.log(ratio);
			sumX+=x[i];
		}
		f+=x[4]*(sumX-1);
		return f;
	}

	@Override
	public double[] gradient(double[] x, int block) {
		double[] df=new double[5];
		double sumX=0;
		for (int i = 0; i < 4; i++) {
			if(x[i]>1)
				x[i]=0.999999;
			if(x[i]<0)
				x[i]=0.00000001;
			double ratio=llr/b[i];
			df[i]=-count_matrix[i]*(ratio+ratio*Math.log(ratio*x[i])+ratio*ratio*x[i])/Math.pow((1+ratio*x[i]),2)+x[4];
			sumX+=x[i];
		}
		df[4]=(sumX-1);
		return df;
	}

}
