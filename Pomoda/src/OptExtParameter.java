import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.io.OutputStream;
import edu.stanford.rsl.jpop.AdditiveFunctionAssembler;
import edu.stanford.rsl.jpop.FunctionOptimizer;
import edu.stanford.rsl.jpop.SimpleFunctionController;
import edu.stanford.rsl.jpop.FunctionOptimizer.OptimizationMode;


public class OptExtParameter  {

public OptExtParameter(double[][] a, double[][] b) {
		super();
		A = a;
		B = b;
		optimalcol=new double[A.length][4];
		colgoodness=new double[A.length];
	}
double[][] A; //sum{llr/b}
double[][] B; //sum{llr/b*log(llr/b)}
double[][] optimalcol;
double[]  colgoodness;
double[][] count_matrix;
double[][] log_bg_matrix;
double llr;

public void find_parameter2(double llr)
{
	for (int i = 0; i < A.length; i++) {
		if(A[i][0]==0)
			continue;
		double[] b=new double[4];
		for (int j = 0; j < b.length; j++) {
			b[j]=Math.exp(log_bg_matrix[i][j]);
		}
		FullEEM_auxfunction func=new FullEEM_auxfunction(count_matrix[i], llr, b);
		FunctionOptimizer opt=new FunctionOptimizer(5);
		opt.setInitialX(new double[]{0.25,0.25,0.25,0.25,1});
		opt.setOptimizationMode(OptimizationMode.Gradient);
		opt.setFunctionController(new SimpleFunctionController());
	//		opt.setFunctionController(new ParallelFunctionController());
		opt.setFunctionAssembler(new AdditiveFunctionAssembler());
		opt.optimizeFunction(func);
		double[] x=opt.getOptimum();
		for (int j = 0; j < 4; j++) {
			optimalcol[i][j]=x[j];
		}
	}
}

public void find_parameter()
{
    FileOutputStream fos;
    PrintStream original = System.out;
	System.setOut(new PrintStream(new OutputStream() {
	    public void write(int b) {
	        //DO NOTHING
	    }
	}));
		for (int i = 0; i < A.length; i++) {
			//compute FG model first
			
			if(A[i][0]==0)
				continue;
			EEM_auxfunction func=new EEM_auxfunction(A[i], B[i]);
			FunctionOptimizer opt=new FunctionOptimizer(1);
			
			opt.setInitialX(new double[]{guess_intiX(A,B,i)});
			opt.setOptimizationMode(OptimizationMode.Gradient);
			opt.setFunctionController(new SimpleFunctionController());
		//		opt.setFunctionController(new ParallelFunctionController());
			opt.setFunctionAssembler(new AdditiveFunctionAssembler());
			opt.optimizeFunction(func);
			double[] x=opt.getOptimum();
			optimalcol[i]=solveOptColumn(x[0],A,B,i);
			
			double goodness=0;
			for (int j = 0; j < 4; j++) {
				double ratio=llr*optimalcol[i][j]/Math.exp(log_bg_matrix[i][j]);
				double Z=ratio/(1+ratio);
				goodness+=count_matrix[i][j]*Z*Math.log(ratio);
					//B[i][j]*optimalcol[i][j]+A[i][j]*optimalcol[i][j]*Math.log(optimalcol[i][j]);
			}
			colgoodness[i]=goodness;
		}
	
	 System.setOut(original);
}

double guess_intiX(double[][] A,double[][] B,int col)
{
	double x=0;
	double coef=0;
	double ApB=0;
	for (int i = 0; i < 4; i++) {
		coef+=1.0/A[col][i];
		ApB+=(A[col][i]+B[col][i])/A[col][i];
	}
	x=(ApB+4*Math.log(0.25))/coef;
	return x;
}

double[] solveOptColumn(double lamda, double[][] A,double[][] B, int col)
{
	double[] optcol=new double[4];
	for (int i = 0; i < 4; i++) {
		optcol[i]=Math.exp((lamda-A[col][i]-B[col][i])/A[col][i]);
	}
	
	
	return optcol;
}

}
