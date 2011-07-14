package edu.stanford.rsl.jpop;

/**
 * This interface provides methods to assemble the information after the block-wise processing of the evaluate, gradient, and hessian.
 * @author akmaier
 *
 */
public interface FunctionAssembler {
	
	/**
	 * This method is provided with an array of block results and assembles them into a single function evaluation result.
	 * @param blockResults the results of the processing blocks.
	 * @return the evaluation result
	 */
	public double assembleEvaluationBlocks (double [] blockResults);
	
	/**
	 * This method is provided with an array of gradients which were obtained by the different processing blocks. The methods assembles the data into a single gradient vector. The first index of the gradient array is supposed to refer to the block number. 
	 * @param blockResults the block results
	 * @return the assembled gradient vector
	 */
	public double [] assembleGradientBlocks (double [] [] blockResults);
	
	
	/**
	 * This method is provided with an array of hessians as computed fron the parallel procesing blocks. The first index of the array should refer to the block index.
	 * @param blockResults the block results.
	 * @return the assembled Hessian.
	 */
	public double [] [] assembleHessianBlocks (double [] [] [] blockResults);
	
}
