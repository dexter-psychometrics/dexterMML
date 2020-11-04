#ifndef DXM_MINIMIZE_
#define DXM_MINIMIZE_

#include <RcppArmadillo.h>

// dfpmin is borrowed from numerical recipes, does not have a publishable license


/*
possible error codes:
0: no problem
1: Roundoff problem in line_search.
2: max iter in line_search reached
3: too many iterations in dfpmin
4: too many iterations in NR
5: solve(hessian, gradient) not succesfull, hessian may not be invertible
*/

template <class T>
void line_search(const arma::vec& xold, const double fold, const arma::vec& g, arma::vec& p,
			arma::vec& x, double &f, const double stpmax, bool &check, T &func, int& err) 
{
	const double ALF=1.0e-4, TOLX=std::numeric_limits<double>::epsilon();
	double a,alam,alam2=0.0,alamin,b,disc,f2=0.0;
	double rhs1,rhs2,tmplam;

	check=false;
	const double sum = std::sqrt(arma::accu(arma::square(p)));
	if (sum > stpmax) 
		p *= stpmax/sum;

	const double slope = arma::accu(g % p);

	if (slope >= 0)
	{
		err = 1;
		return;
	}
	const double test = arma::max(arma::abs(p/arma::clamp(xold,-1.0,1.0)));

	alamin=TOLX/test;
	alam=1.0;
	for (int iter=0; iter<200; iter++) {
		x = xold + alam * p;

		f=func(x);
		if (alam < alamin) {
			x = xold;
			check=true;
			return;
		} else if (f <= fold+ALF*alam*slope) return;
		else {
			if (alam == 1.0)
				tmplam = -slope/(2.0*(f-fold-slope));
			else {
				rhs1=f-fold-alam*slope;
				rhs2=f2-fold-alam2*slope;
				a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
				b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
				if (a == 0.0) tmplam = -slope/(2.0*b);
				else {
					disc=b*b-3.0*a*slope;
					if (disc < 0.0) tmplam=0.5*alam;
					else if (b <= 0.0) tmplam=(-b+std::sqrt(disc))/(3.0*a);
					else tmplam=-slope/(b+std::sqrt(disc));
				}
				if (tmplam>0.5*alam)
					tmplam=0.5*alam;
			}
		}
		alam2=alam;
		f2 = f;
		alam=std::max(tmplam,0.1*alam);
	}
	err = 2;
}


template <class T>
void dfpmin(arma::vec& p, const double gtol, int &iter, double &fret, T &funcd, int& err)
{
	const int ITMAX=200;
	const double EPS=std::numeric_limits<double>::epsilon();
	const double TOLX=4*EPS,STPMX=0.1; //STPMX needs to be considerably smaller than in NR
	bool check;
	double den,fac,fad,fae,fp,sumdg,sumxi,test;
	int n=p.n_elem;
	arma::vec dg(n,arma::fill::zeros), g(n), hdg(n,arma::fill::zeros), pnew(n,arma::fill::zeros), xi(n);
	arma::mat hessinv(n,n,arma::fill::zeros);
	err = 0;
	fp = funcd(p);
	funcd.df(p,g);
	hessinv.diag().ones();
	xi = -g;

	const double stpmax = STPMX * std::max(std::sqrt(arma::accu(arma::square(p))), (double)n);
	
	for (int its=0;its<ITMAX;its++) {
		iter=its;
		line_search(p,fp,g,xi,pnew,fret,stpmax,check,funcd, err);
		if(err>0) return;
		fp = fret;
		xi = pnew-p;
		p = pnew;
		test = arma::max(arma::abs( xi / arma::clamp(p,-1.0,1.0)));
		if (test < TOLX)
			return;
		dg = g;
		funcd.df(p,g);
				
		den = std::max(std::abs(fret),1.0);
		test = arma::max(arma::abs(g % arma::clamp(p,-1.0,1.0)))/den;
		if (test < gtol)
			return;
			
		dg = g-dg;
		hdg.zeros();
		for (int i=0;i<n;i++)
			for (int j=0;j<n;j++) hdg[i] += hessinv.at(i,j)*dg[j];

		fac = arma::accu(dg % xi);
		fae = arma::accu(dg % hdg);
		sumdg = arma::accu(arma::square(dg));
		sumxi = arma::accu(arma::square(xi));

		if (fac > std::sqrt(EPS*sumdg*sumxi)) {
			fac = 1/fac;
			fad = 1/fae;
			dg = fac * xi - fad * hdg;			
			for (int i=0;i<n;i++) {
				for (int j=i;j<n;j++) {
					hessinv.at(i,j) += fac*xi[i]*xi[j]-fad*hdg[i]*hdg[j]+fae*dg[i]*dg[j];
					hessinv.at(j,i) = hessinv.at(i,j);
				}
			}
		}
		xi.zeros();
		for (int i=0;i<n;i++)
			for (int j=0;j<n;j++) 
				xi[i] -= hessinv.at(i,j)*g[j];
	}
	err=3;
}



template <class T>
void NRmin(arma::vec& p, const double gtol, int &iter, double &fret, T &funcd, int& err)
{
	const int max_iter=50, npar=p.n_elem;

	err=0;	
	arma::vec g(npar), delta(npar,arma::fill::zeros);
	arma::mat h(npar,npar);
	
	
	for(iter=0; iter<max_iter; iter++)
	{
		funcd.df(p, g);
		funcd.hess(p, h);
		
		bool slvd = arma::solve(delta, h, g);
		if(!slvd || arma::any(arma::abs(delta)>10))
		{
			err=5;
			return;
		}		
		
		p += delta;
		if(arma::max(arma::abs(g)) < gtol)
		{
			fret = funcd(p);
			return;
		}
	}
	err=4;
}


#endif
