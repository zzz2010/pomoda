package edu.stanford.rsl.jpop.testing;

import edu.stanford.rsl.jpop.SimpleFunctionOptimizer;
import edu.stanford.rsl.jpop.SimpleOptimizableFunction;
import edu.stanford.rsl.jpop.utils.UserUtil;

/**
*
*This class tests the Fzero class.
*
*@author Steve Verrill
*@version .5 --- April 18, 2001
*
*Modified to fit to new API.
*April 25, 2010 
*@author akmaier
*
*/


public class FzeroTest extends Object implements SimpleOptimizableFunction {

   int id_f_to_zero;
   double d,e,f;

   FzeroTest(int idtemp, double dtemp, double etemp, double ftemp) {

      id_f_to_zero = idtemp;
      d = dtemp;
      e = etemp;
      f = ftemp;

   }

   public static void main (String args[]) {

      int another;
      int idtemp;
      double dtemp,etemp,ftemp;

      double b[] = new double[2];
      double c[] = new double[2];
      double r,re,ae;

      int iflag[] = new int[2];

      dtemp = etemp = ftemp = 0.0;

      another = 1;

      while (another == 1) { 
    	 try{
         idtemp = UserUtil.queryInt("\nFor which function do you " +
         "want to find zeros?\n\n" +
         "1 -- (x - d)(x - e)\n" +
         "2 -- (x - d)(x - e)(x - f)\n" +
         "3 -- sin(x)\n\n", 1);

         if (idtemp == 1) {

            dtemp = UserUtil.queryDouble("\nWhat is the d value?  ", 1);
            etemp = UserUtil.queryDouble("\nWhat is the e value?  ", 1);

         } else if (idtemp == 2) {

            dtemp = UserUtil.queryDouble("\nWhat is the d value?  ", 1);
            etemp = UserUtil.queryDouble("\nWhat is the e value?  ", 1);
            ftemp = UserUtil.queryDouble("\nWhat is the f value?  ", 1);

         }

         FzeroTest fzerotest = new FzeroTest(idtemp,dtemp,etemp,ftemp);

         b[1] = UserUtil.queryDouble("\nWhat is the b value?  ", 1);
         c[1] = UserUtil.queryDouble("\nWhat is the c value?  ", 1);
         r = (b[1] + c[1])/2.0;
         re = UserUtil.queryDouble("\nWhat is the re value?  ", 1);
         ae = UserUtil.queryDouble("\nWhat is the ae value?  ", 1);

         
         SimpleFunctionOptimizer opti = new SimpleFunctionOptimizer();
         opti.setAbsoluteTolerance(ae);
         opti.setRelativeTolerance(re);
         opti.setLeftEndPoint(b[1]);
         opti.setRightEndPoint(c[1]);
         opti.setUseInitialGuess(true);
         opti.setInitialGuess(r);
         opti.findRoot(fzerotest);

         System.out.print("\nThe b value is " + b[1] + "\n");
         System.out.print("\nThe iflag value is " + iflag[1] + "\n");

         another = UserUtil.queryInt("\nAnother test?" +
         "   0 - no   1 - yes\n\n", 0);
    	 } catch (Exception e){
    		 e.printStackTrace();
    	 }
      }

      System.out.print("\n");

   }


   public double evaluate(double x) {

      double ff;

      if (id_f_to_zero == 1) {

         ff = (x - d)*(x - e);

      } else if (id_f_to_zero == 2) {

         ff = (x - d)*(x - e)*(x - f);

      } else {

         ff = Math.sin(x);

      }

      return ff;         

   }


}




