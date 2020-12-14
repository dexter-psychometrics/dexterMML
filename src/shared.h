#ifndef DXM_SHARED_
#define DXM_SHARED_

#include <queue>
#include <RcppArmadillo.h>

inline double SQR(const double v){ return v*v; };
inline long double SQR(const long double v){ return v*v; };
inline int SQR(const int v){ return v*v; };
inline double CUB(const double v){ return v*v*v; };
inline long double CUB(const long double v){ return v*v*v; };
inline int CUB(const int v){ return v*v*v; };

inline int kron(const int a, const int b){return a==b ? 1 : 0; };

arma::vec gaussian_pts(const double mu, const double s, const arma::vec& theta);

arma::field<arma::mat>& field_plus(arma::field<arma::mat>& a, const arma::field<arma::mat>& b); 

arma::field<arma::mat> field_init(const arma::field<arma::mat>& orig);

arma::mat mat_init(const arma::mat& orig);

arma::vec vec_init(const arma::vec& orig);

#pragma omp declare reduction( + : arma::field<arma::mat> : field_plus(omp_out, omp_in)) \
initializer( omp_priv = field_init(omp_orig) )

#pragma omp declare reduction( + : arma::mat : omp_out += omp_in ) \
initializer( omp_priv = mat_init(omp_orig) )

#pragma omp declare reduction( + : arma::vec : omp_out += omp_in ) \
initializer( omp_priv = vec_init(omp_orig) )

struct progress
{
	int max_iter, perc, ckInterrupt,calls, bar_width;
	bool draw_bar, draw_any;
	std::string fmt;

	progress(const int max_iter_, const int char_width, const int ckInterrupt_=5)
	{
		perc=0;
		max_iter = max_iter_;
		ckInterrupt = ckInterrupt_;
		bar_width = char_width - 8;
		draw_bar = bar_width >= 5;
		draw_any = char_width >= 4;
		if(draw_bar)
		{
			fmt = std::string("\r|%-") + std::to_string(bar_width) + std::string("s|% 3i%%");			
		} 
		draw();
	}
	void update(const int iter)
	{
		int p = std::min(int(std::round((100.0*iter)/max_iter)),100);
		if(p!=perc)
		{
			perc=p;
			draw();
		}	
		calls++;
		if(calls % ckInterrupt == 0)
			Rcpp::checkUserInterrupt();
	}
	void draw()
	{
		if(draw_bar)
		{
			int bw = std::round(bar_width * perc/100.0);
			Rprintf(fmt.c_str(), std::string(bw,'=').c_str(), perc);
		}
		else if(draw_any) Rprintf("\r% 4i%%", perc);	
	}	
	void close()
	{
		if(draw_any)
		{
			if(perc<100)
			{
				perc=100;
				draw();
			}
			Rprintf("\n");
		}
	}
};

struct progress_est : progress
{
	double target;
	int max_max_iter;	
	bool geo;
	std::queue<double> hist; 
	
	progress_est(const int max_iter_, const int char_width, const double target_=.0001, 
					const int ckInterrupt_=5, bool geometric=true) //geo false is linear progress
		: progress(max_iter_, char_width, ckInterrupt_)
	{
		max_max_iter = max_iter_;
		geo=geometric;
		if(geo)
			target = std::log(target_);
		else
			target = target;
	}
	void update(const double dif, const int iter)
	{
		if(iter>=8)
		{
			if(hist.size() == 5)
				hist.pop();
			if(geo)
				hist.push(std::log(dif));
			else
				hist.push(dif);
		}
		if(hist.size()>=2)
		{
			double a = (hist.back() - hist.front())/hist.size();
			if((target - hist.back())/a > 0)
				max_iter = std::min((int)(std::ceil((target - hist.back())/a) + iter), max_max_iter);
			
		} 
		else max_iter = max_max_iter/2;
		progress::update(iter);
	}
};



#endif