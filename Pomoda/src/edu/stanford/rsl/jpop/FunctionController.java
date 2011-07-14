package edu.stanford.rsl.jpop;

public abstract class FunctionController {
	
	protected FunctionAssembler assembler;
	
	/**
	 * @return the assembler
	 */
	public FunctionAssembler getAssembler() {
		return assembler;
	}

	/**
	 * @param assembler the assembler to set
	 */
	public void setAssembler(FunctionAssembler assembler) {
		this.assembler = assembler;
	}

	/**
	 * Evaluates the given Function and returns the result.
	 * @param function
	 * @return
	 */
	public abstract double evaluate(OptimizableFunction function, double [] x);
	
	/**
	 * Evaluates the given Function's Gradient and returns the result.
	 * @param function
	 * @return
	 */
	public abstract double [] gradient (GradientOptimizableFunction function, double [] x);
	
	/**
	 * Evaluates the given Function's Hessian and returns the result.
	 * @param function
	 * @return
	 */
	public abstract double [] [] hessian (HessianOptimizableFunction function, double [] x);
}
