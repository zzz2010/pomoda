package edu.stanford.rsl.jpop;

/**
 * Implements an additive FunctionAssembler meaning that it computes the sum element by element of each block result.
 * This FunctionAssembler can be used as an example for alternative FunctionAssembler implementations.
 * 
 * @author akmaier
 *
 */
public class AdditiveFunctionAssembler implements FunctionAssembler {

	@Override
	public double assembleEvaluationBlocks(double[] blockResults) {
		// all block results must have the same size
		double revan = 0;
		for (int i = 0; i < blockResults.length; i++){
			revan += blockResults[i];
		}
		return revan;
	}

	@Override
	public double[] assembleGradientBlocks(double[][] blockResults) {
		// all block results must have the same size
		double revan []  = new double [blockResults[0].length];
		for (int i = 0; i < blockResults.length; i++){ // for all blocks
			for (int j = 0; j < blockResults[0].length; j++) { // for all elements in the gradient vector.
				revan[j] += blockResults[i][j];
			}
		}
		return revan;
	}

	@Override
	public double[][] assembleHessianBlocks(double[][][] blockResults) {
		// all block results must have the same size
		double revan [][]  = new double [blockResults[0].length][blockResults[0][0].length];
		for (int i = 0; i < blockResults.length; i++){ // for all blocks
			for (int j = 0; j < blockResults[0].length; j++) { // for all elements in the hessian.
				for (int k = 0; k < blockResults[0][0].length; k++) {
					revan[j][k] += blockResults[i][j][k];
				}
			}
		}
		return revan;
	}

}
