/*************************************************************************
Copyright (c) 2007-2008, Sergey Bochkanov (ALGLIB project).

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

- Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.

- Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer listed
  in this license in the documentation and/or other materials
  provided with the distribution.

- Neither the name of the copyright holders nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*************************************************************************/

#ifndef _lbfgs_h
#define _lbfgs_h

#include "ap.h"

struct lbfgsstate
{
    int n;
    int m;
    double epsg;
    double epsf;
    double epsx;
    int maxits;
    int flags;
    int nfev;
    int mcstage;
    int k;
    int q;
    int p;
    ap::real_1d_array rho;
    ap::real_2d_array y;
    ap::real_2d_array s;
    ap::real_1d_array theta;
    ap::real_1d_array d;
    double stp;
    ap::real_1d_array work;
    double fold;
    double gammak;
    ap::real_1d_array x;
    double f;
    ap::real_1d_array g;
    bool xupdated;
    ap::rcommstate rstate;
    int repiterationscount;
    int repnfev;
    int repterminationtype;
    bool brackt;
    bool stage1;
    int infoc;
    double dg;
    double dgm;
    double dginit;
    double dgtest;
    double dgx;
    double dgxm;
    double dgy;
    double dgym;
    double finit;
    double ftest1;
    double fm;
    double fx;
    double fxm;
    double fy;
    double fym;
    double stx;
    double sty;
    double stmin;
    double stmax;
    double width;
    double width1;
    double xtrapf;
};
struct lbfgsreport
{
    int iterationscount;
    int nfev;
    int terminationtype;
};


/*************************************************************************
        LIMITED MEMORY BFGS METHOD FOR LARGE SCALE OPTIMIZATION

The subroutine minimizes function F(x) of N arguments by  using  a  quasi-
Newton method (LBFGS scheme) which is optimized to use  a  minimum  amount
of memory.

The subroutine generates the approximation of an inverse Hessian matrix by
using information about the last M steps of the algorithm  (instead of N).
It lessens a required amount of memory from a value  of  order  N^2  to  a
value of order 2*N*M.

Input parameters:
    N   -   problem dimension. N>0
    M   -   number of corrections in the BFGS scheme of Hessian
            approximation update. Recommended value:  3<=M<=7. The smaller
            value causes worse convergence, the bigger will  not  cause  a
            considerably better convergence, but will cause a fall in  the
            performance. M<=N.
    X   -   initial solution approximation, array[0..N-1].
    EpsG -  positive number which  defines  a  precision  of  search.  The
            subroutine finishes its work if the condition ||G|| < EpsG  is
            satisfied, where ||.|| means Euclidian norm, G - gradient, X -
            current approximation.
    EpsF -  positive number which  defines  a  precision  of  search.  The
            subroutine finishes its work if on iteration  number  k+1  the
            condition |F(k+1)-F(k)| <= EpsF*max{|F(k)|, |F(k+1)|, 1}    is
            satisfied.
    EpsX -  positive number which  defines  a  precision  of  search.  The
            subroutine finishes its work if on iteration number k+1    the
            condition |X(k+1)-X(k)| <= EpsX is fulfilled.
    MaxIts- maximum number of iterations. If MaxIts=0, the number of
            iterations is unlimited.
    Flags - additional settings:
            * Flags = 0     means no additional settings
            * Flags = 1     "do not allocate memory". used when solving
                            a many subsequent tasks with  same N/M  values.
                            First  call MUST  be without this flag bit set,
                            subsequent calls of MinLBFGS with same LBFGSState
                            structure can set Flags to 1.

Output parameters:
    State - structure used for reverse communication.
    
See also MinLBFGSIteration, MinLBFGSResults

  -- ALGLIB --
     Copyright 14.11.2007 by Bochkanov Sergey
*************************************************************************/
void minlbfgs(const int& n,
     const int& m,
     ap::real_1d_array& x,
     const double& epsg,
     const double& epsf,
     const double& epsx,
     const int& maxits,
     int flags,
     lbfgsstate& state);


/*************************************************************************
One L-BFGS iteration

Called after initialization with MinLBFGS.
See HTML documentation for examples.

Input parameters:
    State   -   structure which stores algorithm state between calls and
                which is used for reverse communication. Must be initialized
                with MinLBFGS.

If suborutine returned False, iterative proces has converged.

If subroutine returned True, caller should calculate function value
State.F an gradient State.G[0..N-1] at State.X[0..N-1] and call
MinLBFGSIteration again.

  -- ALGLIB --
     Copyright 20.04.2009 by Bochkanov Sergey
*************************************************************************/
bool minlbfgsiteration(lbfgsstate& state);


/*************************************************************************
L-BFGS algorithm results

Called after MinLBFGSIteration returned False.

Input parameters:
    State   -   algorithm state (used by MinLBFGSIteration).

Output parameters:
    X       -   array[0..N-1], solution
    Rep     -   optimization report:
                * Rep.TerminationType completetion code:
                    * -2    rounding errors prevent further improvement.
                            X contains best point found.
                    * -1    incorrect parameters were specified
                    *  1    relative function improvement is no more than
                            EpsF.
                    *  2    relative step is no more than EpsX.
                    *  4    gradient norm is no more than EpsG
                    *  5    MaxIts steps was taken
                * Rep.IterationsCount contains iterations count
                * NFEV countains number of function calculations

  -- ALGLIB --
     Copyright 14.11.2007 by Bochkanov Sergey
*************************************************************************/
void minlbfgsresults(const lbfgsstate& state,
     ap::real_1d_array& x,
     lbfgsreport& rep);


#endif
