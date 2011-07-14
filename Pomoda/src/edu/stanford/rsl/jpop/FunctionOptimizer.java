package edu.stanford.rsl.jpop;

import java.util.Arrays;

import edu.stanford.rsl.jpop.fortran.*;

public class FunctionOptimizer extends Object {

	public static enum OptimizationMode {Function, Gradient, Hessian};
	/**
	 *Termination codes
	 *                      ITRMCD =  0:  Optimal solution found
	 *                      ITRMCD =  1:  Terminated with gradient small,
	 *                                  X is probably optimal
	 *                      ITRMCD =  2:  Terminated with stepsize small,
	 *                                  X is probably optimal
	 *                      ITRMCD =  3:  Lower point cannot be found,
	 *                                  X is probably optimal
	 *                      ITRMCD =  4:  Iteration limit (150) exceeded
	 *                      ITRMCD =  5:  Too many large steps,
	 *                                  function may be unbounded
	 */
	public static enum TerminationCode {OptimalSolution, GradientSmall, StepsizeSmall, LowestPointFound, IterationLimit, TooManyLargeStepsProbableUnboundFunction};

	private int dimension;
	private double [] initialX;
	private double [] workSpaceX;
	private double [] functionValueAtX;
	private double [] gradientWorkspace;
	private int [] terminationCode;
	private double [] [] hessianWorkSpace;
	private double [] diagonalWorkspace;
	private UncminForJava uncmin;
	private FunctionController controller;
	private OptimizationMode optimizationMode = OptimizationMode.Function;

	public FunctionOptimizer (int dimension){
		this();
		setDimension(dimension);
	}

	public FunctionOptimizer(){
		controller = new ParallelFunctionController();
		controller.setAssembler(new AdditiveFunctionAssembler());
		uncmin = new UncminForJava(controller);
		setDimension(5);
	}

	/**
	 * Sets a new FunctionAssembler for the optimization process. Default setting is AdditiveFunctionAssembler.
	 * @param assembler
	 */
	public void setFunctionAssembler(FunctionAssembler assembler){
		controller.setAssembler(assembler);
	}

	/**
	 * Sets a new FunctionController. Defaults to ParallelFunctionController. In rare cases the overhead caused by the parallel execution may reduce the actual execution time. In this case use a SimpleFunctionController instead.
	 * @param controller
	 */
	public void setFunctionController(FunctionController controller){
		this.controller = controller;
		uncmin = new UncminForJava(controller);
	}

	/**
	 * Method will optimize the given function using the current settings of the FunctionOptimizer. Note that if you set the OptimizationMode to Gradient or Hessian, the OptimizableFunction must implement GradientOptimizableFunction or HessianOptmizableFunction respectively.
	 * @param function
	 * @return
	 */
	@SuppressWarnings("deprecation")
	public double [] optimizeFunction(OptimizableFunction function){
		// Declarations only needed for Fortran Code
		//
		// Those declared of length 2 should always be declared of length 2.
		// Those declared of length 5 should be declared of length narg+1
		// where narg is the number of arguments over which you are
		// optimizing.

		double typsiz[] = new double[dimension+1];
		double fscale[] = new double[2];
		int method[] = new int[2];
		int iexp[] = new int[2];
		int msg[] = new int[2];
		int ndigit[] = new int[2];
		int itnlim[] = new int[2];
		int iagflg[] = new int[2];
		int iahflg[] = new int[2];
		double dlt[] = new double[2];
		double gradtl[] = new double[2];
		double stepmx[] = new double[2];
		double steptl[] = new double[2];	
		if (optimizationMode == OptimizationMode.Function) {
			uncmin.optimizeFunction0(dimension,initialX,function,workSpaceX,functionValueAtX,gradientWorkspace,terminationCode,hessianWorkSpace,diagonalWorkspace);
		} else if (optimizationMode == OptimizationMode.Gradient) {

			UncminForJava.initialize(dimension,initialX,typsiz,fscale,method,iexp,
					msg,ndigit,itnlim,iagflg,iahflg,
					dlt,gradtl,stepmx,steptl);
			iagflg[1] = 1;
			iahflg[1] = 0;
			iexp[1] = 0;
			uncmin.optimizeFunction7(dimension,initialX,function,typsiz,fscale,method,iexp,
					msg,ndigit,itnlim,iagflg,iahflg,
					dlt,gradtl,stepmx,steptl,
					workSpaceX,functionValueAtX,gradientWorkspace,terminationCode,hessianWorkSpace,diagonalWorkspace);
		} else if (optimizationMode == OptimizationMode.Hessian) {
			UncminForJava.initialize(dimension,initialX,typsiz,fscale,method,iexp,
					msg,ndigit,itnlim,iagflg,iahflg,
					dlt,gradtl,stepmx,steptl);
			iagflg[1] = 1;
			iahflg[1] = 1;
			iexp[1] = 0;
			uncmin.optimizeFunction7(dimension,initialX,function,typsiz,fscale,method,iexp,
					msg,ndigit,itnlim,iagflg,iahflg,
					dlt,gradtl,stepmx,steptl,
					workSpaceX,functionValueAtX,gradientWorkspace,terminationCode,hessianWorkSpace,diagonalWorkspace);
		}
		return getOptimum();
	}

	public void setDimension(int dimension){
		this.dimension = dimension;
		initialX = new double [dimension];
		workSpaceX = new double [dimension+1];
		functionValueAtX  = new double [2];
		gradientWorkspace =  new double [dimension +1];
		terminationCode = new int[2];
		hessianWorkSpace = new double [dimension+1][dimension+1];
		diagonalWorkspace = new double [dimension+1];
	}

	public void setInitialX(double [] x){
		if (x.length > dimension) {
			initialX = Arrays.copyOfRange(x, 0, dimension);
		} else {
			initialX = x;
		}
	}

	public double [] getOptimum(){
		return Arrays.copyOfRange(workSpaceX, 1, workSpaceX.length);
	}

	public double [] getGradientAtOptimum(){
		return Arrays.copyOfRange(gradientWorkspace, 1, gradientWorkspace.length);
	}

	public double getFunctionAtOptimum(){
		return functionValueAtX[1];
	}

	public double [][] getHessianAtOptimum(){
		double [] [] hessian = new double [hessianWorkSpace.length-1][hessianWorkSpace.length-1];
		for (int i = 0; i < hessian.length; i++){
			System.arraycopy(hessianWorkSpace[i+1], 1, hessian[i], 0, hessian[i].length);
		}
		return hessian;
	}

	/**
	 * @return the optimizationMode
	 */
	public OptimizationMode getOptimizationMode() {
		return optimizationMode;
	}

	/**
	 * @param optimizationMode the optimizationMode to set
	 */
	public void setOptimizationMode(OptimizationMode optimizationMode) {
		this.optimizationMode = optimizationMode;
	}


	/**
	 * returns the TerminationCode of the Optimization process. It should be checked to determine, whether the optimization was successful.
	 * @return the TermininationCode
	 */
	public TerminationCode getTerminationCode(){
		/**
		 *Termination codes
		 *                      ITRMCD =  0:  Optimal solution found
		 *                      ITRMCD =  1:  Terminated with gradient small,
		 *                                  X is probably optimal
		 *                      ITRMCD =  2:  Terminated with stepsize small,
		 *                                  X is probably optimal
		 *                      ITRMCD =  3:  Lower point cannot be found,
		 *                                  X is probably optimal
		 *                      ITRMCD =  4:  Iteration limit (150) exceeded
		 *                      ITRMCD =  5:  Too many large steps,
		 *                                  function may be unbounded
		 */
		switch (terminationCode[1]){
		case 0: return TerminationCode.OptimalSolution;
		case 1: return TerminationCode.GradientSmall;
		case 2: return TerminationCode.StepsizeSmall;
		case 3: return TerminationCode.LowestPointFound;
		case 4: return TerminationCode.IterationLimit;
		case 5: return TerminationCode.TooManyLargeStepsProbableUnboundFunction;
		default:
			return null;
		}
	}

}
