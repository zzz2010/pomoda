package edu.stanford.rsl.jpop;

import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import edu.stanford.rsl.jpop.FunctionOptimizer.OptimizationMode;


/**
 * Class to process the gradient evaluations in parallel to speed up computationally expensive function / gradient / Hessian evaluations.
 * @author akmaier
 *
 */
public class ParallelFunctionController extends FunctionController {

	protected boolean debug = true;
	
	@Override
	public double evaluate(OptimizableFunction function, double[] x) {
		BlockThread [] blocks = initializeBlocks(x, function, OptimizationMode.Function);
		try {
			runParallel(blocks);
		} catch (InterruptedException e) {
			e.printStackTrace();
		}
		double [] results = new double [blocks.length];
		for (int i = 0; i < blocks.length; i++){
			results[i] = (Double) blocks[i].getResult();
		}
		return assembler.assembleEvaluationBlocks(results);
	}

	@Override
	public double[] gradient(GradientOptimizableFunction function, double[] x) {
		BlockThread [] blocks = initializeBlocks(x, function, OptimizationMode.Gradient);
		try {
			runParallel(blocks);
		} catch (InterruptedException e) {
			e.printStackTrace();
		}
		double [][] results = new double [blocks.length][];
		for (int i = 0; i < blocks.length; i++){
			results[i] = (double []) blocks[i].getResult();
		}
		return assembler.assembleGradientBlocks(results);
	}

	@Override
	public double[][] hessian(HessianOptimizableFunction function, double[] x) {
		BlockThread [] blocks = initializeBlocks(x, function, OptimizationMode.Hessian);
		try {
			runParallel(blocks);
		} catch (InterruptedException e) {
			e.printStackTrace();
		}
		double [][][] results = new double [blocks.length][][];
		for (int i = 0; i < blocks.length; i++){
			results[i] = (double [][]) blocks[i].getResult();
		}
		return assembler.assembleHessianBlocks(results);
	}
	
	/**
	 * Creates and Prepares the blocks for execution
	 * @param x the point to evaluate
	 * @param function the function
	 * @param mode the evaluation Mode
	 * @return the prepared blocks
	 */
	protected BlockThread [] initializeBlocks(double [] x, OptimizableFunction function, OptimizationMode mode){
		int threads = function.getNumberOfProcessingBlocks();
		BlockThread [] blocks = new BlockThread[function.getNumberOfProcessingBlocks()];
		for (int i = 0; i < threads; i++){
			blocks[i] = new BlockThread();
			blocks[i].setX(x);
			blocks[i].setBlock(i);
			blocks[i].setFunction(function);
			blocks[i].setMode(mode);
		}
		return blocks;
	}
	
	/**
	 * executes the parallel processing.
	 * @param threads
	 * @throws InterruptedException 
	 */
	protected void runParallel(BlockThread [] threads) throws InterruptedException{
		ExecutorService e = Executors.newFixedThreadPool(threads.length);
		CountDownLatch latch = new CountDownLatch(threads.length);
		for (int i = 0; i < threads.length; i++){
			if (debug) System.out.println("Running " + i + " of " + threads.length);
			threads[i].setLatch(latch);
			e.submit(threads[i]);
		}
		latch.await();
		e.shutdownNow();
	}

}
