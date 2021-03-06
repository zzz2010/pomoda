/*************************************************************************
Copyright (c) 2009, Sergey Bochkanov (ALGLIB project).

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

#ifndef _minlm_h
#define _minlm_h

#include "ap.h"

#include "blas.h"
#include "trinverse.h"
#include "cholesky.h"
#include "spdsolve.h"
#include "lbfgs.h"


struct lmstate
{
    bool wrongparams;
    int n;
    int m;
    double epsf;
    double epsx;
    int maxits;
    int flags;
    int usermode;
    ap::real_1d_array x;
    double f;
    ap::real_1d_array fi;
    ap::real_2d_array j;
    ap::real_2d_array h;
    ap::real_1d_array g;
    bool needf;
    bool needfg;
    bool needfgh;
    bool needfij;
    bool xupdated;
    lbfgsstate internalstate;
    lbfgsreport internalrep;
    ap::real_1d_array xprec;
    ap::real_1d_array xbase;
    ap::real_1d_array xdir;
    ap::real_1d_array gbase;
    ap::real_2d_array rawmodel;
    ap::real_2d_array model;
    ap::real_1d_array work;
    ap::rcommstate rstate;
    int repiterationscount;
    int repterminationtype;
    int repnfunc;
    int repnjac;
    int repngrad;
    int repnhess;
    int repncholesky;
};
struct lmreport
{
    int iterationscount;
    int terminationtype;
    int nfunc;
    int njac;
    int ngrad;
    int nhess;
    int ncholesky;
};


/*************************************************************************
    LEVENBERG-MARQUARDT-LIKE METHOD FOR NON-LINEAR OPTIMIZATION

Optimization using function gradient and Hessian.  Algorithm -  Levenberg-
Marquardt   modification   with   L-BFGS   pre-optimization  and  internal
pre-conditioned L-BFGS optimization after each Levenberg-Marquardt step.

Function F has general form (not "sum-of-squares"):

    F = F(x[0], ..., x[n-1])

EXAMPLE

See HTML-documentation.

INPUT PARAMETERS:
    N       -   dimension, N>1
    X       -   initial solution, array[0..N-1]
    EpsF    -   stopping criterion. Algorithm stops if
                |F(k+1)-F(k)| <= EpsF*max{|F(k)|, |F(k+1)|, 1}
    EpsX    -   stopping criterion. Algorithm stops if
                |X(k+1)-X(k)| <= EpsX*(1+|X(k)|)
    MaxIts  -   stopping criterion. Algorithm stops after MaxIts iterations.
                MaxIts=0 means no stopping criterion.

�������� ���������:
    State   -   structure which stores algorithm state between subsequent
                calls of MinLMIteration. Used for reverse communication.
                This structure should be passed to MinLMIteration subroutine.

See also MinLMIteration, MinLMResults.

NOTE

Passing EpsF=0, EpsX=0 and MaxIts=0 (simultaneously) will lead to automatic
stopping criterion selection (small EpsX).

  -- ALGLIB --
     Copyright 30.03.2009 by Bochkanov Sergey
*************************************************************************/
void minlmfgh(const int& n,
     const ap::real_1d_array& x,
     const double& epsf,
     const double& epsx,
     const int& maxits,
     lmstate& state);


/*************************************************************************
    LEVENBERG-MARQUARDT-LIKE METHOD FOR NON-LINEAR OPTIMIZATION

Optimization using function gradient and Jacobian.  Algorithm -  Levenberg-
Marquardt   modification   with   L-BFGS   pre-optimization  and  internal
pre-conditioned L-BFGS optimization after each Levenberg-Marquardt step.

Function F is represented as sum of squares:

    F = f[0]^2(x[0],...,x[n-1]) + ... + f[m-1]^2(x[0],...,x[n-1])

EXAMPLE

See HTML-documentation.

INPUT PARAMETERS:
    N       -   dimension, N>1
    M       -   number of functions f[i]
    X       -   initial solution, array[0..N-1]
    EpsF    -   stopping criterion. Algorithm stops if
                |F(k+1)-F(k)| <= EpsF*max{|F(k)|, |F(k+1)|, 1}
    EpsX    -   stopping criterion. Algorithm stops if
                |X(k+1)-X(k)| <= EpsX*(1+|X(k)|)
    MaxIts  -   stopping criterion. Algorithm stops after MaxIts iterations.
                MaxIts=0 means no stopping criterion.

�������� ���������:
    State   -   structure which stores algorithm state between subsequent
                calls of MinLMIteration. Used for reverse communication.
                This structure should be passed to MinLMIteration subroutine.

See also MinLMIteration, MinLMResults.

NOTE

Passing EpsF=0, EpsX=0 and MaxIts=0 (simultaneously) will lead to automatic
stopping criterion selection (small EpsX).

  -- ALGLIB --
     Copyright 30.03.2009 by Bochkanov Sergey
*************************************************************************/
void minlmfgj(const int& n,
     const int& m,
     const ap::real_1d_array& x,
     const double& epsf,
     const double& epsx,
     const int& maxits,
     lmstate& state);


/*************************************************************************
    CLASSIC LEVENBERG-MARQUARDT METHOD FOR NON-LINEAR OPTIMIZATION

Optimization using Jacobi matrix. Algorithm  -  classic Levenberg-Marquardt
method.

Function F is represented as sum of squares:

    F = f[0]^2(x[0],...,x[n-1]) + ... + f[m-1]^2(x[0],...,x[n-1])

EXAMPLE

See HTML-documentation.

INPUT PARAMETERS:
    N       -   dimension, N>1
    M       -   number of functions f[i]
    X       -   initial solution, array[0..N-1]
    EpsF    -   stopping criterion. Algorithm stops if
                |F(k+1)-F(k)| <= EpsF*max{|F(k)|, |F(k+1)|, 1}
    EpsX    -   stopping criterion. Algorithm stops if
                |X(k+1)-X(k)| <= EpsX*(1+|X(k)|)
    MaxIts  -   stopping criterion. Algorithm stops after MaxIts iterations.
                MaxIts=0 means no stopping criterion.

�������� ���������:
    State   -   structure which stores algorithm state between subsequent
                calls of MinLMIteration. Used for reverse communication.
                This structure should be passed to MinLMIteration subroutine.

See also MinLMIteration, MinLMResults.

NOTE

Passing EpsF=0, EpsX=0 and MaxIts=0 (simultaneously) will lead to automatic
stopping criterion selection (small EpsX).

  -- ALGLIB --
     Copyright 30.03.2009 by Bochkanov Sergey
*************************************************************************/
void minlmfj(const int& n,
     const int& m,
     const ap::real_1d_array& x,
     const double& epsf,
     const double& epsx,
     const int& maxits,
     lmstate& state);


/*************************************************************************
One Levenberg-Marquardt iteration.

Called after inialization of State structure with MinLMXXX subroutine.
See HTML docs for examples.

Input parameters:
    State   -   structure which stores algorithm state between subsequent
                calls and which is used for reverse communication. Must be
                initialized with MinLMXXX call first.

If subroutine returned False, iterative algorithm has converged.

If subroutine returned True, then:
* if State.NeedF=True,      -   function value F at State.X[0..N-1]
                                is required
* if State.NeedFG=True      -   function value F and gradient G
                                are required
* if State.NeedFiJ=True     -   function vector f[i] and Jacobi matrix J
                                are required
* if State.NeedFGH=True     -   function value F, gradient G and Hesian H
                                are required

One and only one of this fields can be set at time.

Results are stored:
* function value            -   in LMState.F
* gradient                  -   in LMState.G[0..N-1]
* Jacobi matrix             -   in LMState.J[0..M-1,0..N-1]
* Hessian                   -   in LMState.H[0..N-1,0..N-1]

  -- ALGLIB --
     Copyright 10.03.2009 by Bochkanov Sergey
*************************************************************************/
bool minlmiteration(lmstate& state);


/*************************************************************************
Levenberg-Marquardt algorithm results

Called after MinLMIteration returned False.

Input parameters:
    State   -   algorithm state (used by MinLMIteration).

Output parameters:
    X       -   array[0..N-1], solution
    Rep     -   optimization report:
                * Rep.TerminationType completetion code:
                    * -1    incorrect parameters were specified
                    *  1    relative function improvement is no more than
                            EpsF.
                    *  2    relative step is no more than EpsX.
                    *  4    gradient norm is no more than EpsG
                    *  5    MaxIts steps was taken
                * Rep.IterationsCount contains iterations count
                * Rep.NFunc     - number of function calculations
                * Rep.NJac      - number of Jacobi matrix calculations
                * Rep.NGrad     - number of gradient calculations
                * Rep.NHess     - number of Hessian calculations
                * Rep.NCholesky - number of Cholesky decomposition calculations

  -- ALGLIB --
     Copyright 10.03.2009 by Bochkanov Sergey
*************************************************************************/
void minlmresults(const lmstate& state, ap::real_1d_array& x, lmreport& rep);


#endif
