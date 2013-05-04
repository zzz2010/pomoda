import cern.colt.function.DoubleDoubleFunction;


public class maxFunction implements DoubleDoubleFunction {

	@Override
	public double apply(double arg0, double arg1) {
		// TODO Auto-generated method stub
		return Math.max(arg0, arg1);
	}

}
