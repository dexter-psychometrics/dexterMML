#ifndef DXM_MINIMIZE_
#define DXM_MINIMIZE_

#include <RcppArmadillo.h>

// borrowed from numerical recipes
// if published, we'll have to bake our own our borrow it from somewhere with a better license
// we can also go for full NR since the hessian is not that hard to derive in most cases but num recp argues against that

// stopgap measure with std::min, should find constrained optim line search
template <class T>
void lnsrch(const arma::vec& xold, const double fold, const arma::vec& g, arma::vec& p,
			arma::vec& x, double &f, const double stpmax, bool &check, T &func) 
{
	const double ALF=1.0e-4, TOLX=std::numeric_limits<double>::epsilon();
	double a,alam,alam2=0.0,alamin,b,disc,f2=0.0;
	double rhs1,rhs2,slope=0.0,sum=0.0,temp,test,tmplam;
	int i,n=xold.n_elem;
	check=false;
	for (i=0;i<n;i++) sum += p[i]*p[i];
	sum=std::sqrt(sum);
	if (sum > stpmax)
		for (i=0;i<n;i++)
			p[i] *= stpmax/sum;
	for (i=0;i<n;i++)
		slope += g[i]*p[i];
	if (slope >= 0.0) Rcpp::stop("Roundoff problem in lnsrch.");
	test=0.0;
	for (i=0;i<n;i++) {
		temp=std::abs(p[i])/std::max(std::abs(xold[i]),1.0);
		if (temp > test) test=temp;
	}
	alamin=TOLX/test;
	alam=1.0;
	for (int iter=0; iter<200; iter++) {
		
		for (i=0;i<n;i++) x[i]=xold[i]+alam*p[i];

		f=func(x);
		if (alam < alamin) {
			for (i=0;i<n;i++) x[i]=xold[i];
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
	Rcpp::stop("max iter in lnsrch reached");
}


template <class T>
void dfpmin(arma::vec& p, const double gtol, int &iter, double &fret, T &funcd)
{
	const int ITMAX=200;
	const double EPS=std::numeric_limits<double>::epsilon();
	const double TOLX=4*EPS,STPMX=0.1; //STPMX needs to be considerably smaller than in NR
	bool check;
	double den,fac,fad,fae,fp,stpmax,sum=0.0,sumdg,sumxi,temp,test;
	int n=p.n_elem;
	arma::vec dg(n,arma::fill::zeros),g(n,arma::fill::zeros),hdg(n,arma::fill::zeros),pnew(n,arma::fill::zeros),xi(n,arma::fill::zeros);
	arma::mat hessin(n,n,arma::fill::zeros);
	fp=funcd(p);
	funcd.df(p,g);
	for (int i=0;i<n;i++) {
		for (int j=0;j<n;j++) hessin.at(i,j)=0.0;
		hessin.at(i,i)=1.0;
		xi[i] = -g[i];
		sum += p[i]*p[i];
	}
	stpmax=STPMX*std::max(std::sqrt(sum),(double)n);
	for (int its=0;its<ITMAX;its++) {
		iter=its;
		lnsrch(p,fp,g,xi,pnew,fret,stpmax,check,funcd);
		fp=fret;
		for (int i=0;i<n;i++) {
			xi[i]=pnew[i]-p[i];
			p[i]=pnew[i];
		}
		test=0.0;
		for (int i=0;i<n;i++) {
			temp=std::abs(xi[i])/std::max(std::abs(p[i]),1.0);
			if (temp > test) test=temp;
		}
		if (test < TOLX)
			return;
		for (int i=0;i<n;i++) dg[i]=g[i];
		funcd.df(p,g);
		test=0.0;
		den=std::max(std::abs(fret),1.0);
		for (int i=0;i<n;i++) {
			temp=std::abs(g[i])*std::max(std::abs(p[i]),1.0)/den;
			if (temp > test) test=temp;
		}
		if (test < gtol)
			return;
		for (int i=0;i<n;i++)
			dg[i]=g[i]-dg[i];
		for (int i=0;i<n;i++) {
			hdg[i]=0.0;
			for (int j=0;j<n;j++) hdg[i] += hessin.at(i,j)*dg[j];
		}
		fac=fae=sumdg=sumxi=0.0;
		for (int i=0;i<n;i++) {
			fac += dg[i]*xi[i];
			fae += dg[i]*hdg[i];
			sumdg += dg[i] * dg[i];
			sumxi += xi[i] * xi[i];
		}
		if (fac > std::sqrt(EPS*sumdg*sumxi)) {
			fac=1.0/fac;
			fad=1.0/fae;
			for (int i=0;i<n;i++) dg[i]=fac*xi[i]-fad*hdg[i];
			for (int i=0;i<n;i++) {
				for (int j=i;j<n;j++) {
					hessin.at(i,j) += fac*xi[i]*xi[j]
						-fad*hdg[i]*hdg[j]+fae*dg[i]*dg[j];
					hessin.at(j,i)=hessin.at(i,j);
				}
			}
		}
		for (int i=0;i<n;i++) {
			xi[i]=0.0;
			for (int j=0;j<n;j++) xi[i] -= hessin.at(i,j)*g[j];
		}
	}
	Rcpp::stop("too many iterations in dfpmin");
}





#endif
