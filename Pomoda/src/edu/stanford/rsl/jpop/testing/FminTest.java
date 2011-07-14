package edu.stanford.rsl.jpop.testing;

import edu.stanford.rsl.jpop.SimpleFunctionOptimizer;
import edu.stanford.rsl.jpop.SimpleOptimizableFunction;
import edu.stanford.rsl.jpop.utils.UserUtil;

/**
 *
 *This class tests the Fmin class.
 *
 *@author Steve Verrill
 *@version .5 --- March 25, 1998
 *
 *Modified to fit new API
 *
 *@author Andreas Maier
 *@version 0.7 --- April 27, 2010
 *
 */


public class FminTest extends Object implements SimpleOptimizableFunction{

	int id_f_to_min;
	double c,d,e;

	FminTest(int idtemp, double ctemp, double dtemp, double etemp) {

		id_f_to_min = idtemp;
		c = ctemp;
		d = dtemp;
		e = etemp;

	}

	public static void main (String args[]) {

		int another;
		int idtemp;
		double ctemp,dtemp,etemp;
		double a,b,tol,xmin;

		ctemp = dtemp = etemp = 0.0;

		another = 1;

		while (another == 1) { 
			try{

				idtemp = UserUtil.queryInt("\nWhat function do you " +
						"want to minimize?\n\n" +
						"1 -- (x - c)(x - d)\n" +
						"2 -- (x - c)(x - d)(x - e)\n" +
						"3 -- sin(x)\n\n",1);

				if (idtemp == 1) {

					ctemp = UserUtil.queryDouble("\nWhat is the c value?  ", 1);
					dtemp = UserUtil.queryDouble("\nWhat is the d value?  ", 1);

				} else if (idtemp == 2) {

					ctemp = UserUtil.queryDouble("\nWhat is the c value?  ", 1);
					dtemp = UserUtil.queryDouble("\nWhat is the d value?  ", 1);
					etemp = UserUtil.queryDouble("\nWhat is the e value?  ", 1);

				}

				FminTest fmintest = new FminTest(idtemp,ctemp,dtemp,etemp);

				a = UserUtil.queryDouble("\nWhat is the a value?  ", 1);
				b = UserUtil.queryDouble("\nWhat is the b value?  ", 1);
				tol = UserUtil.queryDouble("\nWhat is the tol value?  ", 1);

				SimpleFunctionOptimizer opti = new SimpleFunctionOptimizer();
				opti.setAbsoluteTolerance(tol);
				opti.setLeftEndPoint(a);
				opti.setRightEndPoint(b);
				xmin = opti.minimize(fmintest);

				System.out.print("\nThe xmin value is " + xmin + "\n");      

				another = UserUtil.queryInt("\nAnother test?" +
						"   0 - no   1 - yes\n\n", 0);
			}catch (Exception e){
				e.printStackTrace();
			}
		}

		System.out.print("\n");

	}


	public double evaluate(double x) {

		double f;

		if (id_f_to_min == 1) {

			f = (x - c)*(x - d);

		} else if (id_f_to_min == 2) {

			f = (x - c)*(x - d)*(x - e);

		} else {

			f = Math.sin(x);

		}

		return f;         

	}


}
