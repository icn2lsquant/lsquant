#ifndef  CHEBYSHEV_COEFFICIENTS
#define CHEBYSHEV_COEFFICIENTS

#include <complex>		/* for std::vector mostly class*/

using namespace std;

inline
double delta_chebF(double x, const double m)
{
	const double f0 =  sqrt(1.0 - x*x ) ;
	const double fm =  cos(m*acos(x));
	return fm/f0/ M_PI;
};

inline
complex<double> greenR_chebF(double x, const double m)
{
	const complex<double> I(0.0,1.0);
	const double f0 =  sqrt(1.0 - x*x ) ;
	const complex<double> fm =  pow( x - I*f0, m );
	return -I*fm/f0;
};

inline
complex<double> DgreenR_chebF(double x, const double m)
{
	const complex<double> I(0.0,1.0);
	const double f0 =  sqrt(1.0 - x*x ) ;
	const complex<double> fm =  pow( x - I*f0 , m )*(x + I*m*f0);
	return -I*fm/f0/f0/f0;
};

inline 
complex<double> GammaAC_chebF(const int m0,const int m1, double w, double x)
{
	const complex<double> I(0.0,1.0);
	const double gamma1m0 = delta_chebF(x,m0)*M_PI;
	const double gamma1m1 = delta_chebF(x,m1)*M_PI;
	const complex<double> gamma2p = I*greenR_chebF(x+w,m1); 
	const complex<double> gamma2m = I*greenR_chebF(x-w,m0) ;
	return gamma1m0*std::conj(gamma2p)-gamma1m1*(gamma2m);
}





#endif
