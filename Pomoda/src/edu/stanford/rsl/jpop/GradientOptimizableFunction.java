package edu.stanford.rsl.jpop;

public interface GradientOptimizableFunction extends OptimizableFunction {


	/**
	 * Computes the Gradient at position x.
	 * (Note that x is a Fortran Array which starts at 1.)
	 * @param x the position
	 * @param block the block identifier. First block is 0. block is < getNumberOfProcessingBlocks().
	 * @return the gradient at x. (In Fortran Style)
	 */
	double [] gradient(double x[], int block);



}
