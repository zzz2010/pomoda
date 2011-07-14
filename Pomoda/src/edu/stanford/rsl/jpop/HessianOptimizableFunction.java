package edu.stanford.rsl.jpop;

public interface HessianOptimizableFunction extends GradientOptimizableFunction {

	/**
	 * Computes the Hessian matrix at position x.
	 * (Note that x is a Fortran Array which starts at 1.)
	 * @param x
	 * @param block the block identifier. First block is 0. block is < getNumberOfProcessingBlocks().
	 * @return the Hessian at x. (In Fortran Style)
	 */
	double [][] hessian(double x[], int block);
	
}
