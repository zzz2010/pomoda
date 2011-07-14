/*
    Fmin.java copyright claim:

    This software is based on the public domain fmin routine.
    The FORTRAN version can be found at

    www.netlib.org

    This software was translated from the FORTRAN version
    to Java by a US government employee on official time.  
    Thus this software is also in the public domain.


    The translator's mail address is:

    Steve Verrill 
    USDA Forest Products Laboratory
    1 Gifford Pinchot Drive
    Madison, Wisconsin
    53705


    The translator's e-mail address is:

    steve@ws13.fpl.fs.fed.us


 ***********************************************************************

DISCLAIMER OF WARRANTIES:

THIS SOFTWARE IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND. 
THE TRANSLATOR DOES NOT WARRANT, GUARANTEE OR MAKE ANY REPRESENTATIONS 
REGARDING THE SOFTWARE OR DOCUMENTATION IN TERMS OF THEIR CORRECTNESS, 
RELIABILITY, CURRENTNESS, OR OTHERWISE. THE ENTIRE RISK AS TO 
THE RESULTS AND PERFORMANCE OF THE SOFTWARE IS ASSUMED BY YOU. 
IN NO CASE WILL ANY PARTY INVOLVED WITH THE CREATION OR DISTRIBUTION 
OF THE SOFTWARE BE LIABLE FOR ANY DAMAGE THAT MAY RESULT FROM THE USE 
OF THIS SOFTWARE.

Sorry about that.

 ***********************************************************************


History:

Date        Translator        Changes

3/24/98     Steve Verrill     Translated

 */


package edu.stanford.rsl.jpop;

import edu.stanford.rsl.jpop.fortran.Blas_f77;


/**
 * This class was is based on the Fortran-Java translation of Steve Verrill of fzero and fmin.
 * The code was adopted and altered to fit into a hierarchical object-oriented program structure.
 * 
 *@author Andreas Maier
 *@version 0.7 --- April 27, 2010
 * 
 */

public class SimpleFunctionOptimizer extends Object {

	protected double leftEndPoint = Double.MIN_VALUE;
	protected double rightEndPoint = Double.MAX_VALUE;
	protected double absoluteTolerance = 3.0e-8;
	protected double relativeTolerance = 3.0e-4;
	protected boolean printWarnings = true;
	protected boolean useInitialGuess = false;
	protected double initialGuess = 0;

	/**
	 * @return the leftEndPoint of the search interval.
	 */
	public double getLeftEndPoint() {
		return leftEndPoint;
	}

	/**
	 * Sets the left end point of the search interval.<br>
	 * Default is Double.MIN_VALUE.<br>
	 * @param leftEndPoint the leftEndPoint to set
	 */
	public void setLeftEndPoint(double leftEndPoint) {
		this.leftEndPoint = leftEndPoint;
	}

	/**
	 * @return the rightEndPoint of the search interval.
	 */
	public double getRightEndPoint() {
		return rightEndPoint;
	}

	/**
	 * Sets the right end point of the search interval.<br>
	 * Defaults to Double.MAX_VALUE.<BR>
	 * @param rightEndPoint the rightEndPoint to set
	 */
	public void setRightEndPoint(double rightEndPoint) {
		this.rightEndPoint = rightEndPoint;
	}

	/**
	 * @return the absoluteTolerance
	 */
	public double getAbsoluteTolerance() {
		return absoluteTolerance;
	}

	/**
	 * Sets the absolute tolerance of the search.<BR>
	 * Defaults to 3.0E-08. 1.1102e-16 is machine precision.<br>
	 * @param absoluteTolerance the absoluteTolerance to set
	 */
	public void setAbsoluteTolerance(double absoluteTolerance) {
		this.absoluteTolerance = absoluteTolerance;
	}

	/**
	 * @return the relativeTolerance
	 */
	public double getRelativeTolerance() {
		return relativeTolerance;
	}

	/**
	 * Sets the relative tolerance.<BR>
	 * Defaults to 3.0E-04. 1.1102e-16 is machine precision.<br>
	 * @param relativeTolerance the relativeTolerance to set
	 */
	public void setRelativeTolerance(double relativeTolerance) {
		this.relativeTolerance = relativeTolerance;
	}

	/**
	 * @return the printWarnings
	 */
	public boolean isPrintWarnings() {
		return printWarnings;
	}

	/**
	 * Turns warnings on or off.<br>
	 * Default is on.<BR>
	 * @param printWarnings the printWarnings to set
	 */
	public void setPrintWarnings(boolean printWarnings) {
		this.printWarnings = printWarnings;
	}

	/**
	 * @return the useInitialGuess
	 */
	public boolean isUseInitialGuess() {
		return useInitialGuess;
	}

	/**
	 * Turns the use of an initial guess on or off.<BR>Default is off.<br>
	 * @param useInitialGuess the useInitialGuess to set
	 */
	public void setUseInitialGuess(boolean useInitialGuess) {
		this.useInitialGuess = useInitialGuess;
	}

	/**
	 * @return the initialGuess
	 */
	public double getInitialGuess() {
		return initialGuess;
	}

	/**
	 * An initial guess of a zero can be set. <br>
	 * The value will only be used if useInitialGuess is true.<br><br>
	 * A guess of a zero of the function could help in
	 * speeding up convergence. If evaluate(leftEndPoint) and evaluate(initialGuess) have
	 * opposite signs, a root will be found in the interval
	 * (leftEndPoint,inititalGuess); if not, but evaluate(initialGuess) and evaluate(rightEndPoint) have opposite
	 * signs, a root will be found in the interval (initialGuess,rightEndPoint);
	 * otherwise, the interval (leftEndPoint,rightEndPoint) will be searched for a
	 * possible root.  
	 * @param initialGuess the initialGuess to set
	 */
	public void setInitialGuess(double initialGuess) {
		this.initialGuess = initialGuess;
	}

	/**
	 * Minimizes a SimpleOptimizableFunction. Searches between leftEndPoint and rightEndPoint with a tolerance of absoluteTolerance. Parameters must be set using the respective getters and setters.
	 * @param function the function to minimize.
	 * @return the position of the minimum.
	 */
	public double minimize(SimpleOptimizableFunction function){
		return fmin(leftEndPoint, rightEndPoint, function, absoluteTolerance);
	}

	/**
	 *
	 *<p>
	 *This method performs a 1-dimensional minimization.
	 *It implements Brent's method which combines a golden-section
	 *search and parabolic interpolation.  The introductory comments from
	 *the FORTRAN version are provided below.
	 *
	 *This method is a translation from FORTRAN to Java of the Netlib 
	 *function fmin.  In the Netlib listing no author is given.
	 *
	 *Translated by Steve Verrill, March 24, 1998.
	 *
	 *@param  leftEndPoint         Left endpoint of initial interval
	 *@param  rightEndPoint         Right endpoint of initial interval
	 *@param  function  A class that defines a method, evaluate(),
	 *                  to minimize.  The class must implement
	 *                  the SimpleOptimizableFunction interface.
	 *                  evaluate must have one
	 *                  double valued argument.
	 *@param  toleranceFactor       Desired length of the interval in which
	 *                  the minimum will be determined to lie
	 *                  (This should be greater than, roughly, 3.0e-8.)
	 *@return the minimal value;
	 */
	@Deprecated
	public static double fmin (double leftEndPoint, double rightEndPoint, SimpleOptimizableFunction function,
			double toleranceFactor) {

		/*

Here is a copy of the Netlib documentation:

c
c      An approximation x to the point where f attains a minimum on
c  the interval (ax,bx) is determined.
c
c  input..
c
c  ax    left endpoint of initial interval
c  bx    right endpoint of initial interval
c  f     function subprogram which evaluates f(x) for any x
c        in the interval (ax,bx)
c  tol   desired length of the interval of uncertainty of the final
c        result (.ge.0.)
c
c  output..
c
c  fmin  abcissa approximating the point where  f  attains a
c        minimum
c
c     The method used is a combination of golden section search and
c  successive parabolic interpolation.  Convergence is never much slower
c  than that for a Fibonacci search.  If f has a continuous second
c  derivative which is positive at the minimum (which is not at ax or
c  bx), then convergence is superlinear, and usually of the order of
c  about 1.324.
c     The function f is never evaluated at two points closer together
c  than eps*abs(fmin)+(tol/3), where eps is approximately the square
c  root of the relative machine precision.  If f is a unimodal
c  function and the computed values of f are always unimodal when
c  separated by at least eps*abs(x)+(tol/3), then fmin approximates
c  the abcissa of the global minimum of f on the interval (ax,bx) with
c  an error less than 3*eps*abs(fmin)+tol.  If f is not unimodal,
c  then fmin may approximate a local, but perhaps non-global, minimum to
c  the same accuracy.
c      This function subprogram is a slightly modified version of the
c  Algol 60 procedure localmin given in Richard Brent, Algorithms For
c  Minimization Without Derivatives, Prentice-Hall, Inc. (1973).
c

		 */

		double c,d,e,eps,xm,p,q,r,tol1,t2,
		u,v,w,fu,fv,fw,fx,x,tol3;

		c = .5*(3.0 - Math.sqrt(5.0));
		d = 0.0;

		// 1.1102e-16 is machine precision

		eps = 1.2e-16;
		tol1 = eps + 1.0;
		eps = Math.sqrt(eps);

		v = leftEndPoint + c*(rightEndPoint-leftEndPoint);
		w = v;
		x = v;
		e = 0.0;
		fx = function.evaluate(x);
		fv = fx;
		fw = fx;
		tol3 = toleranceFactor/3.0;

		xm = .5*(leftEndPoint + rightEndPoint);
		tol1 = eps*Math.abs(x) + tol3;
		t2 = 2.0*tol1;

		// main loop

		while (Math.abs(x-xm) > (t2 - .5*(rightEndPoint-leftEndPoint))) {

			p = q = r = 0.0;

			if (Math.abs(e) > tol1) {

				// fit the parabola
				r = (x-w)*(fx-fv);
				q = (x-v)*(fx-fw);
				p = (x-v)*q - (x-w)*r;
				q = 2.0*(q-r);
				if (q > 0.0) {
					p = -p;
				} else {
					q = -q;
				}
				r = e;
				e = d;         
			}
			if ((Math.abs(p) < Math.abs(.5*q*r)) &&
					(p > q*(leftEndPoint-x)) &&
					(p < q*(rightEndPoint-x))) {
				// a parabolic interpolation step
				d = p/q;
				u = x+d;
				// f must not be evaluated too close to a or b
				if (((u-leftEndPoint) < t2) || ((rightEndPoint-u) < t2)) {
					d = tol1;
					if (x >= xm) d = -d;
				}
			} else {
				// a golden-section step
				if (x < xm) {
					e = rightEndPoint-x;
				} else {
					e = leftEndPoint-x;
				}
				d = c*e;
			}
			// f must not be evaluated too close to x
			if (Math.abs(d) >= tol1) {
				u = x+d;
			} else {
				if (d > 0.0) {
					u = x + tol1;
				} else {
					u = x - tol1;
				}
			}
			fu = function.evaluate(u);
			// Update a, b, v, w, and x
			if (fx <= fu) {
				if (u < x) {
					leftEndPoint = u;
				} else {
					rightEndPoint = u;
				}
			}
			if (fu <= fx) {
				if (u < x) {
					rightEndPoint = x;
				} else {
					leftEndPoint = x;
				}
				v = w;
				fv = fw;
				w = x;
				fw = fx;
				x = u;
				fx = fu;
				xm = .5*(leftEndPoint + rightEndPoint);
				tol1 = eps*Math.abs(x) + tol3;
				t2 = 2.0*tol1;
			} else {
				if ((fu <= fw) || (w == x)) {
					v = w;
					fv = fw;
					w = u;
					fw = fu;
					xm = .5*(leftEndPoint + rightEndPoint);
					tol1 = eps*Math.abs(x) + tol3;
					t2 = 2.0*tol1;
				} else if ((fu > fv) && (v != x) && (v != w)) {
					xm = .5*(leftEndPoint + rightEndPoint);
					tol1 = eps*Math.abs(x) + tol3;
					t2 = 2.0*tol1;
				} else {
					v = u;
					fv = fu;
					xm = .5*(leftEndPoint + rightEndPoint);
					tol1 = eps*Math.abs(x) + tol3;
					t2 = 2.0*tol1;
				}
			}
		}
		return x;
	}   

	/**
	 * Method which prints all warnings. <BR>
	 * May be overloaded in a sub class of SimpleFunctionOptimizer to redirect the warnings.
	 * @param warning the warning.
	 */
	protected static void warn(String warning){
		System.out.println("Warning: " + warning);
	}

	/**
	 * Method to find a root of a function between leftEndPoint and rightEndPoint.<BR> May use an initial guess to speed to convergence.<BR><BR>Accuracy can be adjusted using relativeTolerance and absoluteTolerance. Use setters of this class to adjust the parameters or to disable the printing of warnings.<BR>
	 * @param function
	 * @return
	 */
	public double findRoot(SimpleOptimizableFunction function){
		double [] b = {0, leftEndPoint};
		double [] c = {0, rightEndPoint};
		int [] iflag = {0 ,0};
		if (!useInitialGuess) initialGuess = leftEndPoint;
		fzero(function, b, c, initialGuess,relativeTolerance, absoluteTolerance, iflag);
		if (printWarnings) {
			switch (iflag[1]){
			case 3:
				warn("Absolute value of zero position is greater than either of the absolute values of the end points! i.e. " +
						" Math.abs(evaluate(findRoot(function))) > " +
						" Math.max(Math.abs(evaluate(leftEndPoint))," +
				"        Math.abs(evaluate(rightEndPoint)))");
				break;
			case 4:
				warn("No change in sign of evaluate(x) was found although the " +
						" interval (leftEndPoint, rightEndPoint) collapsed to the requested tolerance. " +
						" You should examine this case and decide whether " +
						" findRoot(function) is near a local minimum of evaluate(x), or it is near a" +
				" zero of even multiplicity, or neither of these.");
				break;
			case 5:
				warn("Too many (> 500) function evaluations used.");
				break;
			}	
		}
		return b[1];
	}

	/**
	 *
	 *<p>
	 *This method searches for a zero of a function f(x) between
	 *the given values b and c until the width of the interval
	 *(b,c) has collapsed to within a tolerance specified by
	 *the stopping criterion, Math.abs(b-c) <= 2.0*(rw*Math.abs(b)+ae).
	 *The method used is an efficient combination of bisection
	 *and the secant rule.
	 *The introductory comments from
	 *the FORTRAN version are provided below.
	 *
	 *This method is a translation from FORTRAN to Java of the Netlib
	 *(actually it is in the SLATEC library which is available at Netlib)
	 *function dfzero.  In the FORTRAN code, L.F. Shampine and H.A. Watts
	 *are listed as the authors.
	 *
	 *Translated by Steve Verrill, April 17, 2001.
	 *</p>
	 *@param  function   A class that defines a method, evaluate,
	 *                   that returns a value that is to be zeroed.
	 *                   The class must implement
	 *                   the SimpleOptimizableFunction interface. 
	 *                   evaluate must have one
	 *                   double valued argument.
	 *@param  b[1]       One end point of the initial interval.  The value returned
	 *                   for b[1] is usually the better approximation to a
	 *                   zero of evaluate.
	 *@param  c[1]       The other end point of the initial interval.
	 *@param  r          Initial guess of a zero.  A (better) guess of a zero 
	 *                   of evaluate() which could help in
	 *                   speeding up convergence.  If evaluate(b) and evaluate(r) have
	 *                   opposite signs, a root will be found in the interval
	 *                   (b,r); if not, but evaluate(r) and evaluate(c) have opposite
	 *                   signs, a root will be found in the interval (r,c);
	 *                   otherwise, the interval (b,c) will be searched for a
	 *                   possible root.  When no better guess is known, it is
	 *                   recommended that r be set to b or c; because if r is
	 *                   not interior to the interval (b,c), it will be ignored.
	 *@param re          Relative error used for rw in the stopping criterion.
	 *                   If the requested re is less than machine precision,
	 *                   then rw is set to approximately machine precision.
	 *@param ae          Absolute error used in the stopping criterion.  If the
	 *                   given interval (b,c) contains the origin, then a
	 *                   nonzero value should be chosen for AE.
	 *
	 *@param iflag[1]   A status code.  User must check iflag[1] after each call.
	 *                  Control returns to the user from dfzero in all cases.
	 *
	 *                  1:  b is within the requested tolerance of a zero.
	 *                      The interval (b,c) collapsed to the requested
	 *                      tolerance, the function changes sign in (b,c), and
	 *                      evaluate(x) decreased in magnitude as (b,c) collapsed.
	 *
	 *                  2:  evaluate(b) = 0.  However, the interval (b,c) may not have
	 *                      collapsed to the requested tolerance.
	 *
	 *                  3:  b may be near a singular point of evaluate(x).
	 *                      The interval (b,c) collapsed to the requested tol-
	 *                      erance and the function changes sign in (b,c), but
	 *                      evaluate(x) increased in magnitude as (b,c) collapsed,i.e.
	 *                      Math.abs(evaluate(b out)) > 
	 *                      Math.max(Math.abs(evaluate(b in)),
	 *                               Math.abs(evaluate(c in)))
	 *
	 *                  4:  No change in sign of evaluate(x) was found although the
	 *                      interval (b,c) collapsed to the requested tolerance.
	 *                      The user must examine this case and decide whether
	 *                      b is near a local minimum of evaluate(x), or b is near a
	 *                      zero of even multiplicity, or neither of these.
	 *
	 *                  5:  Too many (> 500) function evaluations used.
	 *
	 */
	@Deprecated
	public static void fzero (SimpleOptimizableFunction function, double b[], double c[], 
			double r, double re, double ae, int iflag[]) {
		/*
	Here is a copy of the Netlib documentation:

	      SUBROUTINE DFZERO(F,B,C,R,RE,AE,IFLAG)
	C***BEGIN PROLOGUE  DFZERO
	C***DATE WRITTEN   700901   (YYMMDD)
	C***REVISION DATE  861211   (YYMMDD)
	C***CATEGORY NO.  F1B
	C***KEYWORDS  LIBRARY=SLATEC,TYPE=DOUBLE PRECISION(FZERO-S DFZERO-D),
	C             BISECTION,NONLINEAR,NONLINEAR EQUATIONS,ROOTS,ZEROES,
	C             ZEROS
	C***AUTHOR  SHAMPINE,L.F.,SNLA
	C           WATTS,H.A.,SNLA
	C***PURPOSE  Search for a zero of a function F(X) in a given
	C            interval (B,C).  It is designed primarily for problems
	C            where F(B) and F(C) have opposite signs.
	C***DESCRIPTION
	C
	C       **** Double Precision version of FZERO ****
	C
	C     Based on a method by T J Dekker
	C     written by L F Shampine and H A Watts
	C
	C            DFZERO searches for a zero of a function F(X) between
	C            the given values B and C until the width of the interval
	C            (B,C) has collapsed to within a tolerance specified by
	C            the stopping criterion, DABS(B-C) .LE. 2.*(RW*DABS(B)+AE).
	C            The method used is an efficient combination of bisection
	C            and the secant rule.
	C
	C     Description Of Arguments
	C
	C     F,B,C,R,RE and AE are DOUBLE PRECISION input parameters.
	C     B and C are DOUBLE PRECISION output parameters and IFLAG (flagged
	C        by an * below).
	C
	C        F     - Name of the DOUBLE PRECISION valued external function.
	C                This name must be in an EXTERNAL statement in the
	C                calling program.  F must be a function of one double
	C                precision argument.
	C
	C       *B     - One end of the interval (B,C).  The value returned for
	C                B usually is the better approximation to a zero of F.
	C
	C       *C     - The other end of the interval (B,C)
	C
	C        R     - A (better) guess of a zero of F which could help in
	C                speeding up convergence.  If F(B) and F(R) have
	C                opposite signs, a root will be found in the interval
	C                (B,R); if not, but F(R) and F(C) have opposite
	C                signs, a root will be found in the interval (R,C);
	C                otherwise, the interval (B,C) will be searched for a
	C                possible root.  When no better guess is known, it is
	C                recommended that r be set to B or C; because if R is
	C                not interior to the interval (B,C), it will be ignored.
	C
	C        RE    - Relative error used for RW in the stopping criterion.
	C                If the requested RE is less than machine precision,
	C                then RW is set to approximately machine precision.
	C
	C        AE    - Absolute error used in the stopping criterion.  If the
	C                given interval (B,C) contains the origin, then a
	C                nonzero value should be chosen for AE.
	C
	C       *IFLAG - A status code.  User must check IFLAG after each call.
	C                Control returns to the user from FZERO in all cases.
	C                XERROR does not process diagnostics in these cases.
	C
	C                1  B is within the requested tolerance of a zero.
	C                   The interval (B,C) collapsed to the requested
	C                   tolerance, the function changes sign in (B,C), and
	C                   F(X) decreased in magnitude as (B,C) collapsed.
	C
	C                2  F(B) = 0.  However, the interval (B,C) may not have
	C                   collapsed to the requested tolerance.
	C
	C                3  B may be near a singular point of F(X).
	C                   The interval (B,C) collapsed to the requested tol-
	C                   erance and the function changes sign in (B,C), but
	C                   F(X) increased in magnitude as (B,C) collapsed,i.e.
	C                     abs(F(B out)) .GT. max(abs(F(B in)),abs(F(C in)))
	C
	C                4  No change in sign of F(X) was found although the
	C                   interval (B,C) collapsed to the requested tolerance.
	C                   The user must examine this case and decide whether
	C                   B is near a local minimum of F(X), or B is near a
	C                   zero of even multiplicity, or neither of these.
	C
	C                5  Too many (.GT. 500) function evaluations used.
	C***REFERENCES  L. F. SHAMPINE AND H. A. WATTS, *FZERO, A ROOT-SOLVING
	C                 CODE*, SC-TM-70-631, SEPTEMBER 1970.
	C               T. J. DEKKER, *FINDING A ZERO BY MEANS OF SUCCESSIVE
	C                 LINEAR INTERPOLATION*, 'CONSTRUCTIVE ASPECTS OF THE
	C                 FUNDAMENTAL THEOREM OF ALGEBRA', EDITED BY B. DEJON
	C                 P. HENRICI, 1969.
	C***ROUTINES CALLED  D1MACH
	C***END PROLOGUE  DFZERO
	c
		 */
		int ic,kount;
		double a,acbs,acmb,aw,cmb,er,fa,fb,fc,fx,fz,p,q,rw,
		t,tol,z;
		// 1.1102e-16 is machine precision
		er = 2.0*1.2e-16;
		// Initialize
		z = r;
		if (r <= Math.min(b[1],c[1]) || r >= Math.max(b[1],c[1])) {
			z = c[1];
		}
		rw = Math.max(re,er);
		aw = Math.max(ae,0.0);
		ic = 0;
		t = z;
		fz = function.evaluate(t);
		fc = fz;
		t = b[1];
		fb = function.evaluate(t);
		kount = 2;
		if (Blas_f77.sign_f77(1.0,fz) != Blas_f77.sign_f77(1.0,fb)) {
			c[1] = z;
		} else {
			if (z != c[1]) {
				t = c[1];
				fc = function.evaluate(t);
				kount = 3;
				if (Blas_f77.sign_f77(1.0,fz) != Blas_f77.sign_f77(1.0,fc)) {
					b[1] = z;
					fb = fz;
				}
			}
		}
		a = c[1];
		fa = fc;
		acbs = Math.abs(b[1] - c[1]);
		fx = Math.max(Math.abs(fb),Math.abs(fc));
		// Main loop
		while (true) {
			if (Math.abs(fc) < Math.abs(fb)) {
				// Perform interchange
				a = b[1];
				fa = fb;
				b[1] = c[1];
				fb = fc;
				c[1] = a;
				fc = fa;
			}
			cmb = 0.5*(c[1] - b[1]);
			acmb = Math.abs(cmb);
			tol = rw*Math.abs(b[1]) + aw;
			// Test stopping criterion and function count.
			if (acmb <= tol) {
				// Finished.  Process the results for the proper setting
				// of the flag.
				if (Blas_f77.sign_f77(1.0,fb) == Blas_f77.sign_f77(1.0,fc)) {
					iflag[1] = 4;
					return;
				}
				if (Math.abs(fb) > fx) {
					iflag[1] = 3;
					return;
				}
				iflag[1] = 1;
				return;
			}

			if (fb == 0.0) {
				iflag[1] = 2;
				return;
			}

			if (kount >= 500) {
				iflag[1] = 5;
				return;
			}


			/*
			 *              CALCULATE NEW ITERATE IMPLICITLY AS B+P/Q
			 *              WHERE WE ARRANGE P .GE. 0.
			 *              THE IMPLICIT FORM IS USED TO PREVENT OVERFLOW.
			 */

			p = (b[1] - a)*fb;
			q = fa - fb;
			if (p < 0.0) {
				p = -p;
				q = -q;
			}

			/*
			 *              UPDATE A AND CHECK FOR SATISFACTORY REDUCTION
			 *              IN THE SIZE OF THE BRACKETING INTERVAL.
			 *              IF NOT, PERFORM BISECTION.
			 */

			a = b[1];
			fa = fb;
			ic++;
			if (ic >= 4 && 8.0*acmb >= acbs) {
				// Use bisection
				b[1] = .5*(c[1] + b[1]);
			} else {
				if (ic >= 4) {
					ic = 0;
					acbs = acmb;
				}
				// Test for too small a change.
				if (p <= Math.abs(q)*tol) {
					// Increment by the tolerance.
					b[1] += Blas_f77.sign_f77(tol,cmb);
				} else {
					// Root ought to be between b[1] and (c[1] + b[1])/2.
					if (p >= cmb*q) {
						// Use bisection.
						b[1] = 0.5*(c[1] + b[1]);
					} else {
						// Use the secant rule.
						b[1] += p/q;
					}
				}
			}
			// Have completed the computation for a new iterate b[1].
			t = b[1];
			fb = function.evaluate(t);
			kount++;
			// Decide whether the next step is an interpolation
			// or an extrapolation.
			if (Blas_f77.sign_f77(1.0,fb) == Blas_f77.sign_f77(1.0,fc)) {
				c[1] = a;
				fc = fa;
			}
		}
	}

}
