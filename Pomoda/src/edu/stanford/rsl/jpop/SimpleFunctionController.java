package edu.stanford.rsl.jpop;

/**
 * Implements a single threaded FunctionController. Class is mainly implemented for testing and debugging.
 * @author akmaier
 *
 */
public class SimpleFunctionController extends FunctionController {

	@Override
	public double evaluate(OptimizableFunction function, double[] x) {
		int threads = function.getNumberOfProcessingBlocks();
		double [] results = new double [threads];
		for (int i = 0; i < threads; i++){
			results[i] = function.evaluate(x, i);
		}
		return assembler.assembleEvaluationBlocks(results);
	}

	@Override
	public double[] gradient(GradientOptimizableFunction function, double[] x) {
		int threads = function.getNumberOfProcessingBlocks();
		double [][] results = new double [threads][];
		for (int i = 0; i < threads; i++){
			results[i] = function.gradient(x, i);
		}
		return assembler.assembleGradientBlocks(results);
	}

	@Override
	public double[][] hessian(HessianOptimizableFunction function, double[] x) {
		int threads = function.getNumberOfProcessingBlocks();
		double [][][] results = new double [threads][][];
		for (int i = 0; i < threads; i++){
			results[i] = function.hessian(x, i);
		}
		return assembler.assembleHessianBlocks(results);
	}


}
