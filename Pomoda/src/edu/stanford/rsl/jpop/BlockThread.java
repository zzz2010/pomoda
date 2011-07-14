package edu.stanford.rsl.jpop;

import java.util.Arrays;
import java.util.concurrent.CountDownLatch;

import edu.stanford.rsl.jpop.FunctionOptimizer.OptimizationMode;

/**
 * Class to run the block evaluation in an independent Thread.
 * @author akmaier
 *
 */
public class BlockThread extends Thread {
	
	private OptimizationMode mode;
	private double [] x;
	private int block;
	private OptimizableFunction function;
	private Object result;
	private CountDownLatch latch;

	/**
	 * @return the mode
	 */
	public OptimizationMode getMode() {
		return mode;
	}
	
	/**
	 * @param mode the mode to set
	 */
	public void setMode(OptimizationMode mode) {
		this.mode = mode;
	}
	
	/**
	 * @return the x
	 */
	public double[] getX() {
		return x;
	}

	/**
	 * @param x the x to set
	 */
	public void setX(double[] x) {
		this.x = Arrays.copyOf(x, x.length);
	}

	/**
	 * @return the block
	 */
	public int getBlock() {
		return block;
	}

	/**
	 * @param block the block to set
	 */
	public void setBlock(int block) {
		this.block = block;
	}

	/**
	 * @return the function
	 */
	public OptimizableFunction getFunction() {
		return function;
	}
	
	/**
	 * @param function the function to set
	 */
	public void setFunction(OptimizableFunction function) {
		this.function = function;
	}

	/**
	 * @return the result
	 */
	public Object getResult() {
		return result;
	}

	/**
	 * @param latch the latch to set
	 */
	public void setLatch(CountDownLatch latch) {
		this.latch = latch;
	}

	/**
	 * @return the latch
	 */
	public CountDownLatch getLatch() {
		return latch;
	}

	@Override
	public void run(){
		if (mode == OptimizationMode.Function){
			result = new Double(function.evaluate(x, block));
		} else if (mode == OptimizationMode.Gradient) {
			result = ((GradientOptimizableFunction)function).gradient(x, block);
		} else if (mode == OptimizationMode.Hessian) {
			result = ((HessianOptimizableFunction) function).hessian(x, block);
		}
		latch.countDown();
	}

	
}
