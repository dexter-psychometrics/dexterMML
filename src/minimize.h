#ifndef DXM_MINIMIZE_
#define DXM_MINIMIZE_

#include <RcppArmadillo.h>
#include "shared.h"

/*****************************************************************************************************************
* This implementation of nlm is a partial translation from R source code which is a translation of fortran uncmin 
* which was written by Dennis and Schnabel
*
* dtrsl (used by nlm code) is a cpp translation from fortran code in linpack 
* originally developed by Jack Dongarra, Jim Bunch, Cleve Moler and Pete Stewart.
* and published under the MIT license
*
* license information for the original R code follows below
******************************************************************************************************************/


/*
possible error codes:
0: no problem
other: problem
*/

// simple 1 dimensional NR minimization for convex functions
// since nlm does not work well for 1D problems and dichotomous rasch should be a convex function
template <class T>
void D1min(arma::vec& pars, const double gtol, int &iter, double &fret, T &funcd, int& err)
{
	if(pars.n_elem != 1)
	{
		err = 1; // problem is not 1d
		return;
	}
	const int max_iter = 200;
	const double min_step = 1e-10, max_step=5;
	double step;
	
	err = 0;
	arma::vec g(1);
	arma::mat h(1,1);

	for(iter=1; iter<=max_iter; iter++)
	{
		funcd.df(pars,g);
		funcd.hess(pars,h);
	
		step = g[0]/h[0];
		if(std::abs(step)>max_step)
			step = std::copysign(max_step, step);
				
		pars[0] -= step;
		
		if(std::abs(g[0]) < gtol || std::abs(step) < min_step )
			break;
		
	}
	if(iter==max_iter)
		err=2;
	
}

/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1999-2017   The R Core Team
 *  Copyright (C) 2003-2017   The R Foundation
 *  Copyright (C) 1997-1999   Saikat DebRoy
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 */

/* ../appl/uncmin.f
   -- translated by f2c (version of 1 June 1993 23:00:00).
   -- and hand edited by Saikat DebRoy
   */

/*--- The Dennis + Schnabel Minimizer -- used by R's  nlm() ---*/


namespace RNLM{


template <class T>
static void R_lnsrch(const arma::vec& x, const double f, const arma::vec& g, arma::vec& p, arma::vec& xpls,
       double& fpls, T& fcn, bool& mxtake, int& iretcd,
       const double stepmx, const double steptl, const arma::vec& sx)
{
/* Find a next newton iterate by line search. 

 * PARAMETERS :

 *	x(n)	    --> old iterate:	x[k-1]
 *	f	     	--> function value at old iterate, f(x)
 *	g(n)	    --> gradient at old iterate, g(x), or approximate
 *	p(n)	    --> non-zero newton step
 *	xpls(n)	    <--	new iterate x[k]
 *	fpls	    <--	function value at new iterate, f(xpls)
 *	fcn	        --> name of subroutine to evaluate function
 *	iretcd	    <--	return code
 *	mxtake	    <-- boolean flag indicating step of maximum length used
 *	stepmx	    --> maximum allowable step size
 *	steptl	    --> relative step size at which successive iterates considered close enough to terminate algorithm
 *	sx(n)	    --> diagonal scaling matrix for x

 *	internal variables

 *	sln		 newton length
 *	rln		 relative length of newton step
*/
	const int n = x.n_elem;
    int i;
    bool firstback = true;
    double disc;
    double a3, b;
    double t1, t2, t3, lambda, tlmbda, rmnlmb;
    double scl, rln, sln, slp;
    double temp1;
    double pfpls = 0, plmbda = 0; /* -Wall */

    temp1 = 0.;
    for (i = 0; i < n; ++i)
		temp1 += sx[i] * sx[i] * p[i] * p[i];
    sln = std::sqrt(temp1);
    if (sln > stepmx) 
	{
		/*	newton step longer than maximum allowed */
		scl = stepmx / sln;
		//F77_CALL(dscal)(&n, &scl, p, &one); //scales a vector by a constant.
		p *= scl;
		sln = stepmx;
    }
    //slp = F77_CALL(ddot)(&n, g, &one, p, &one); //dot product of two vectors
    slp = arma::dot(g,p);
	rln = 0.;
    for (i = 0; i < n; ++i) 
	{
		temp1 = std::abs(p[i])/ std::max(std::abs(x[i]), 1.0/sx[i]);
		if(rln < temp1) rln = temp1;
    }	
	
    rmnlmb = steptl / rln;
    lambda = 1.0;

    /*	check if new iterate satisfactory.  generate new lambda if necessary. */

    mxtake = false;
    iretcd = 2;
	int ls_iter=0;
	const int max_ls_iter=200;
    do 
	{
		xpls = x + lambda * p;
		fpls = fcn(xpls);
		//(*fcn)(n, xpls, fpls, state);
		if (fpls <= f + slp * 1e-4 * lambda) 
		{ 
			/* solution found */
			iretcd = 0;
			if (lambda == 1. && sln > stepmx * .99) mxtake = true;
			return;
		}
		/* else : solution not (yet) found */

		/* First find a point with a finite value */

		if (lambda < rmnlmb) 
		{
			/* no satisfactory xpls found sufficiently distinct from x */
			iretcd = 1;
			return;
		}
		else 
		{ 
			/*	calculate new lambda */
			if (!std::isfinite(fpls))
			{
				lambda *= 0.1;
				firstback = true;
			}
			else 
			{
				if (firstback) 
				{ 
					/*	first backtrack: quadratic fit */
					tlmbda = -lambda * slp / ((fpls - f - slp) * 2.);
					firstback = false;
				}
				else 
				{ 
					/*	all subsequent backtracks: cubic fit */
					t1 = fpls - f - lambda * slp;
					t2 = pfpls - f - plmbda * slp;
					t3 = 1. / (lambda - plmbda);
					a3 = 3. * t3 * (t1 / (lambda * lambda)
						   - t2 / (plmbda * plmbda));
					b = t3 * (t2 * lambda / (plmbda * plmbda)
						  - t1 * plmbda / (lambda * lambda));
					disc = b * b - a3 * slp;
					if (disc > b * b) /* only one positive critical point, must be minimum */
						tlmbda = (-b + ((a3 < 0) ? -std::sqrt(disc) : std::sqrt(disc))) /a3;
					else /* both critical points positive, first is minimum */
						tlmbda = (-b + ((a3 < 0) ? std::sqrt(disc) : -std::sqrt(disc))) /a3;

					if (tlmbda > lambda * .5)
						tlmbda = lambda * .5;
				}
				plmbda = lambda;
				pfpls = fpls;
				if (tlmbda < lambda * .1)
					lambda *= 0.1;
				else
					lambda = tlmbda;
			}
		}
		ls_iter++;
		if(ls_iter>max_ls_iter)
			break; // with iretcd=2
    } while(iretcd > 1);
} 

static void choldc(arma::mat& a, const double diagmx, const double tol, double& addmax)
{
	int nr = a.n_cols;
	int n=nr;

    double tmp1, tmp2;
    int i, j, k;
    double aminl, offmax, amnlsq;
    double sum;

    addmax = 0.0;
    aminl = sqrt(diagmx * tol);
    amnlsq = aminl * aminl;

    /*	form row i of L */

    for (i = 0; i < n; ++i) 
	{
		// A[i,j] := * || find i,j element of lower triangular matrix L
		for (j = 0; j < i; ++j) 
		{
			sum = 0.;
			for (k = 0; k < j; ++k)
				sum += a[i + k * nr] * a[j + k * nr];
			a[i + j * nr] = (a[i + j * nr] - sum) / a[j + j * nr];
		}

		// A[i,i] := * || find diagonal elements of L
		sum = 0.;
		for (k = 0; k < i; ++k)
			sum += a[i + k * nr] * a[i + k * nr];

		tmp1 = a[i + i * nr] - sum;
		if (tmp1 >= amnlsq) { // normal Cholesky
			a[i + i * nr] = std::sqrt(tmp1);
		}
		else { // augment diagonal of L
			/*	find maximum off-diagonal element in row */
			offmax = 0.;
			for (j = 0; j < i; ++j) {
			if(offmax < (tmp2 = fabs(a[i + j * nr])))
				offmax = tmp2;
			}
			if (offmax <= amnlsq) offmax = amnlsq;

			/* add to diagonal element to
			 * allow cholesky decomposition to continue */
			a[i + i * nr] = std::sqrt(offmax);
			if(addmax < (tmp2 = offmax - tmp1)) addmax = tmp2;
		}
    }
} /* choldc */


static void chlhsn(arma::mat& a, const arma::vec& sx, arma::vec& udiag, const double epsm)
{
	int nr = a.n_cols;
	int n = nr;

    int i, j;
    double evmin, evmax;
    double addmax, diagmn, diagmx, offmax, offrow, posmax;
    double sdd, amu, tol, tmp;


    /*	scale hessian */
    /*	pre- and post- multiply "a" by inv(sx) */

    for (j = 0; j < n; ++j)
		for (i = j; i < n; ++i)
			a[i + j * nr] /= sx[i] * sx[j];

    /*	step1
     *	-----
     *	note:  if a different tolerance is desired throughout this
     *	algorithm, change tolerance here: */

    tol = std::sqrt(epsm);

    diagmx = a[0];
    diagmn = a[0];
    if (n > 1) {
		for (i = 1; i < n; ++i) 
		{
			tmp = a[i + i * nr];
			if(diagmn > tmp)	diagmn = tmp;
			if(diagmx < tmp)	diagmx = tmp;
		}
    }
    posmax = std::max(diagmx, 0.0);

    if (diagmn <= posmax * tol) 
	{
		amu = tol * (posmax - diagmn) - diagmn;
		if (amu == 0.) {
			/*	find largest off-diagonal element of a */
			offmax = 0.;
			for (i = 1; i < n; ++i) {
			for (j = 0; j < i; ++j)
				if (offmax < (tmp = std::fabs(a[i + j * nr]))) offmax = tmp;
			}

			if (offmax == 0.)
			amu = 1.;
			else
			amu = offmax * (tol + 1.);
		}
		/*	a=a + mu*i */
		for (i = 0; i < n; ++i)
			a[i + i * nr] += amu;

		diagmx += amu;
    }

    /*	copy lower triangular part of "a" to upper triangular part */
    /*	and diagonal of "a" to udiag */
    for (i = 0; i < n; ++i) 
	{
		udiag[i] = a[i + i * nr];
		for (j = 0; j < i; ++j)
			a[j + i * nr] = a[i + j * nr];
    }
    choldc(a, diagmx, tol, addmax);


    /*	step3
     *	if addmax=0, "a" was positive definite going into step 2,
     *	the LL+ decomposition has been done, and we return.
     *	otherwise, addmax>0.  perturb "a" so that it is safely
     *	diagonally dominant and find LL+ decomposition */

    if (addmax > 0.0) 
	{
		/*	restore original "a" (lower triangular part and diagonal) */

		for (i = 0; i < n; ++i) {
			a[i + i * nr] = udiag[i];
			for (j = 0; j < i; ++j)
			a[i + j * nr] = a[j + i * nr];
		}

		/*	find sdd such that a+sdd*i is safely positive definite */
		/*	note:  evmin<0 since a is not positive definite; */

		evmin = 0.;
		evmax = a[0];
		for (i = 0; i < n; ++i) {
			offrow = 0.;
			for (j = 0; j < i; ++j)
			offrow += std::fabs(a[i + j * nr]);
			for (j = i+1; j < n; ++j)
			offrow += std::fabs(a[j + i * nr]);
			tmp = a[i + i * nr] - offrow;
			if(evmin > tmp) evmin = tmp;
			tmp = a[i + i * nr] + offrow;
			if(evmax < tmp) evmax = tmp;
		}
		sdd = tol * (evmax - evmin) - evmin;

		/*	perturb "a" and decompose again */

		amu = std::min(sdd, addmax);
		for (i = 0; i < n; ++i) 
		{
			a[i + i * nr] += amu;
			udiag[i] = a[i + i * nr];
		}

		/*	 "a" now guaranteed safely positive definite */

		choldc(a, 0.0, tol, addmax);
    }
    /*	unscale hessian and cholesky decomposition matrix */

    for (j = 0; j < n; ++j) 
	{
		for (i = j; i < n; ++i)
			a[i + j * nr] *= sx[i];
		for (i = 0; i < j; ++i)
			a[i + j * nr] *= sx[i] * sx[j];
		udiag[j] *= sx[j] * sx[j];
    }
} /* chlhsn */




static void dtrsl(const arma::mat& t, arma::vec& b, const int job, int& info)
{
	const int n = t.n_cols;

    double ddot,temp;
    int j,jj;

	if(arma::any(t.diag()==0.0))
		return;
    info = 0;

	/*
	in R:
		b[n]=b[n]/a[n,n]
		for(jj in 2:n)
		{
		  j = n-jj+1
		  b[j] = b[j] - sum(a[(j+1):n,j] * b[(j+1):n])
		  b[j] = b[j]/a[j,j]
		}
	*/
	if(job==10)
	{
		b[n-1] = b[n-1]/t.at(n-1,n-1);
		if(n<2) return;
		for(jj=1;jj<n;jj++)
		{
			j = n-jj-1;
			ddot=0;
			for(int i=j+1; i<n;i++)
				ddot += t.at(i,j) * b[i];
			b[j] -= ddot;
			b[j] = b[j]/t.at(j,j);		
		}
	}

	/*
	in R:
		b[1]=b[1]/a[1,1]
		n=2
		for(j in 2:n)
		{
		  temp = -b[j-1]
		  l = n-j
		  b[j:(l+j)] = temp * a[j:(l+j),j-1] + b[j:(l+j)] 
		  b[j] = b[j]/a[j,j]
		}	
	*/
	
	if(job == 0)
	{
		b[0] = b[0]/t.at(0,0);
		if(n<2) return;
		for(j=1; j<n; j++)
		{
			temp = -b[j-1];
			for(jj=j; jj<n; jj++)
				b[jj] = temp * t.at(jj,j-1) + b[jj];
			b[j] = b[j] / t.at(j,j);		
		}	
	}
}


static void lltslv(arma::mat& a, arma::vec& x,  const arma:: vec& b)
{
    int job = 0, info;
	x = b;
	dtrsl(a, x, job, info);
    job = 10;
    dtrsl(a, x, job, info);
} 

static int opt_stop(const arma::vec& xpls, double fpls, const arma::vec& gpls, const arma::vec& x,
	     const int itncnt, int& icscmx, const double gradtl, const double steptl,
	     const arma::vec& sx, const double fscale, const int itnlim,
	     const int iretcd, const bool mxtake)
{
/* Unconstrained minimization stopping criteria :

 * Find whether the algorithm should terminate, due to any
 * of the following:
 *	1) problem solved within user tolerance
 *	2) convergence within user tolerance
 *	3) iteration limit reached
 *	4) divergence or too restrictive maximum step (stepmx) suspected

 * ARGUMENTS :

 *	xpls(n)	     --> new iterate x[k]
 *	fpls	     --> function value at new iterate f(xpls)
 *	gpls(n)	     --> gradient at new iterate, g(xpls), or approximate
 *	x(n)	     --> old iterate x[k-1]
 *	itncnt	     --> current iteration k
 *	icscmx	    <--> number consecutive steps >= stepmx
 *			 		[retain value between successive calls]
 *	gradtl	     --> tolerance at which relative gradient considered close
 *			 			enough to zero to terminate algorithm
 *	steptl	     --> relative step size at which successive iterates
 *			 			considered close enough to terminate algorithm
 *	sx(n)	     --> diagonal scaling matrix for x
 *	fscale	     --> estimate of scale of objective function
 *	itnlim	     --> maximum number of allowable iterations
 *	iretcd	     --> return code
 *	mxtake	     --> boolean flag indicating step of maximum length used
 *
 * VALUE :
 *	`itrmcd' : termination code
 */

	const int n = x.n_elem;
    int i, jtrmcd;
    double d, relgrd, relstp, rgx, rsx;

    /*	last global step failed to locate a point lower than x */
    if (iretcd == 1)
		return 3;

    /* else : */

    /* find direction in which relative gradient maximum. */

    /* check whether within tolerance */
    d = std::max(std::abs(fpls), fscale);
    rgx = 0.;
    for (i = 0; i < n; ++i) 
	{
		relgrd = std::abs(gpls[i]) * std::max(std::abs(xpls[i]), 1./sx[i]) / d;
		if(rgx < relgrd) 
			rgx = relgrd;
    }
    jtrmcd = 1;
    if (rgx > gradtl) 
	{
		if (itncnt == 0)
			return 0;

		/* find direction in which relative stepsize maximum */
		/* check whether within tolerance. */
		rsx = 0.;
		for (i = 0; i < n; ++i) 
		{
			relstp = std::abs(xpls[i] - x[i]) / std::max(std::abs(xpls[i]), 1./sx[i]);
			if(rsx < relstp) 
				rsx = relstp;
		}
		jtrmcd = 2;
		if (rsx > steptl) 
		{ 
			/*	check iteration limit */
			jtrmcd = 4;
			if (itncnt < itnlim) 
			{
				/*	check number of consecutive steps \ stepmx */
				if (!mxtake) {
					icscmx = 0; 
					return 0;
				} 
				else 
				{
					icscmx++;
					if (icscmx < 5) 
						return 0;
					jtrmcd = 5;
				}
			}
		}
    }
    return jtrmcd;
} /* opt_stop */



static void optchk(const arma::vec& x, arma::vec& typsiz, arma::vec& sx, double& fscale,
				   const double gradtl,  double& dlt,
				   double& stepmx,
				   int& msg)
{
/* Check input for reasonableness.
 * Return *msg in {-1,-2,..,-7}	 if something is wrong

 * PARAMETERS :

 *	x(n)	     --> on entry, estimate to root of fcn
 *	typsiz(n)   <--> typical size of each component of x
 *	sx(n)	    <--	 diagonal scaling matrix for x
 *	fscale	    <--> estimate of scale of objective function fcn
 *	gradtl	     --> tolerance at which gradient considered close
 *			 			enough to zero to terminate algorithm
 *	dlt	        <--> trust region radius
 *	stepmx	    <--> maximum step size
 *	msg	    	<--> message and error code
 */
    const int n = x.n_elem;
	int i;
    double stpsiz;

    /*	check dimension of problem */
    if (n == 1 && msg % 2 == 0) {
		msg = -2;  
		return;
    }

    /*	compute scale matrix */
    for (i = 0; i < n; ++i) 
	{
		if (typsiz[i] == 0.)
			typsiz[i] = 1.;
		else if (typsiz[i] < 0.)
			typsiz[i] = -typsiz[i];
		sx[i] = 1. / typsiz[i];
    }

    /*	compute default maximum step size if not provided */
    if (stepmx <= 0.) 
	{
		stpsiz = 0.;
		for (i = 0; i < n; ++i)
			stpsiz += SQR(x[i]) * SQR(sx[i]);
		stepmx = 1000. * std::max(std::sqrt(stpsiz), 1.0);
    }

    /*	check function scale */
    if (fscale == 0.)
		fscale = 1.0;
    else if (fscale < 0.)
		fscale = -(fscale);

    /*	check gradient tolerance */
    if (gradtl < 0.) 
	{
		msg = -3;  
		return;
    }

    /*	check trust region radius */
    if (dlt <= 0.) 
	{
		dlt = -1.0;
    } 
	else if (dlt > stepmx)
	{
		dlt = stepmx;
    }
    return;
} /* optchk */


static void optdrv_end(arma::vec& xpls, const arma::vec& x, arma::vec& gpls, const arma::vec& g, double& fpls, const double f, int& msg, const int itrmcd)
{
    int i;
	const int n = xpls.n_elem;
    /*	termination :
	reset xpls,fpls,gpls,  if previous iterate solution */
    if (itrmcd == 3) 
	{
		fpls = f;
		for (i = 0; i < n; ++i) {
			xpls[i] = x[i];
			gpls[i] = g[i];
		}
    }
    msg = 0;
} /* optdrv_end */

// length of param must be >=2
template <class T>
void R_nlm(arma::vec& x, T& fcn, 
       arma::vec& typsiz, double fscale, 
       int& msg,  const int itnlim, 
       double dlt, double gradtl, double stepmx, double steptl,
       arma::vec& xpls, double& fpls, arma::vec& gpls, int& itrmcd,
       int& itncnt)
{
/* Driver for non-linear optimization problem  

 * PARAMETERS :

 *	x(n)	 	--> on entry: estimate to a root of fcn
 *	fcn	     	--> name of subroutine to evaluate optimization function
 *	typsiz(n)   --> typical size for each component of x
 *	fscale	    --> estimate of scale of objective function
 *	msg	   	    --> on output: ( < 0) error code; =0 no error
 *	itnlim	    --> maximum number of allowable iterations
 *	dlt	        --> trust region radius
 *	gradtl	    --> tolerance at which gradient considered close
 *					 enough to zero to terminate algorithm
 *	stepmx	     --> maximum allowable step size
 *	steptl	     --> relative step size at which successive iterates
 *			 			considered close enough to terminate algorithm
 *	xpls(n)	    <--> on exit:  xpls is local minimum
 *	fpls	    <--> on exit:  function value at solution, xpls
 *	gpls(n)	    <--> on exit:  gradient at solution xpls
 *	itrmcd	    <--	 termination code
 *	itncnt	     current iteration, k  {{was `internal'}}


 *	internal variables

 *	analtl		 tolerance for comparison of estimated and
 *			 analytical gradients and hessians
 *	epsm		 machine epsilon
 *	f		 function value: fcn(x)
 *	rnf		 relative noise in optimization function fcn.
 *			      noise=10.**(-ndigit)
 */
    const int n = x.n_elem, ndigit=15;
	bool mxtake = false;
    int iretcd = 0, icscmx = 0;
    double epsm, f, analtl;
    double rnf;
	
	arma::mat a(n,n);
	arma::vec g(n), p(n,arma::fill::zeros), sx(n), udiag(n);

    itncnt = 0;
    epsm = std::numeric_limits<double>::epsilon();
    optchk(x, typsiz, sx, fscale, gradtl, dlt, stepmx, msg);
    if (msg < 0) return;

    rnf = std::pow(10.0, -(double)ndigit);	
    rnf = std::max(rnf, epsm);
    analtl = std::sqrt(rnf);
    analtl = std::max(0.1, analtl);

    /*	evaluate fcn(x) */
	f = fcn(x);
	fcn.df(x, g);
	
    iretcd = -1;
	xpls.zeros();
    itrmcd = opt_stop(x, f, g, xpls, itncnt, icscmx, gradtl, steptl, sx, 
						fscale, itnlim, iretcd, false);
    if (itrmcd != 0) 
	{
		optdrv_end(xpls, x, gpls, g, fpls, f, msg, itrmcd);
		xpls=x;
		fpls=f;
		return;
    }

	/* analytic hessian */
	fcn.hess(x,a);
    
    /* THE Iterations : */

    while(true)
	{
		itncnt++;
		chlhsn(a, sx, udiag,epsm);
		lltslv(a, p, -g);

		R_lnsrch(x, f, g, p, xpls, fpls, fcn, mxtake,
			   iretcd, stepmx, steptl, sx);
	

		/*	calculate step for output */
		p = xpls-x;

		/*	calculate gradient at xpls */
		fcn.df(xpls,gpls);
		

		/*	check whether stopping criteria satisfied */
		itrmcd = opt_stop(xpls, fpls, gpls, x, itncnt, icscmx,
				   gradtl, steptl, sx, fscale, itnlim, iretcd, mxtake);

		if(itrmcd != 0) break;

		/*	evaluate hessian at xpls */
		fcn.hess(xpls,a);
		//a(0,1)=0;

		/*	x <-- xpls  and	 g <-- gpls  and  f <-- fpls */
		f = fpls;
		x = xpls;
		g = gpls;
    } /* END while(1) */

	optdrv_end(xpls, x, gpls, g, fpls, f, msg, itrmcd);
} /* optdrv */
}

// simple wrapper
// set some parameters to good defaults (from R and trial/error)
template <class T>
void nlm(arma::vec& p, const double gtol, int &iter, double &fret, T &funcd, int& err)
{
	const int n = p.n_elem;
	const int itnlim = 200;
	int itrmcd = 0;
	arma::vec xpls(n), gpls(n);
	arma::vec typsiz = funcd.typical_size;
	err=0;
	double fscale=1; // estimate of size of f at minimum, leave at R default for now
	double dlt = 1; //trust region radius
	/*R: max(1000 * sqrt(sum((p/typsize)^2)), 1000)*/
	//p.print("p:");
	//typsiz.print("typsiz:");
	//(p/typsiz).print("div:");
	
	double stepmx = std::max(10 * std::sqrt(arma::accu(arma::square(p/typsiz))),10.0);
	double steptl = 1e-6; //R has 1e-6
	
	RNLM::R_nlm(p, funcd,
				   typsiz, fscale, 
				   err,  itnlim, 
				   dlt, gtol, stepmx, steptl,
				   xpls, fret, gpls, itrmcd,
				   iter);
	p = xpls; 
}



#endif