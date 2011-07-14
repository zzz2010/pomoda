package edu.stanford.rsl.jpop;

public interface OptimizableFunction {
	
	/**
	 * Sets the number of parallel processing blocks. This number should optimally equal to the number of available processors in the current machine. 
	 * Default value should be 1. 
	 * @param number
	 */
	public void setNumberOfProcessingBlocks(int number);

	/**
	 * returns the number of parallel processing blocks.
	 * @return
	 */
	public int getNumberOfProcessingBlocks();
	
	/**
	 * Evaluates the function at position x.<BR> 
	 * (Note that x is a Fortran Array which starts at 1.)
	 * @param x the position
	 * @param block the block identifier. First block is 0. block is < getNumberOfProcessingBlocks().
	 * @return the function value at x
	 */
	public double evaluate(double x[], int block);
}
