#include "special_functions.hpp"

double chebyshev::besselJ(const int n, const double x) {
  
  switch (n) {
  case 0:
    return gsl_sf_bessel_Jn(0,x);
//    return j0(x);
//    return std::sph_bessel ( 0,x);
    break;
  case 1:
    return gsl_sf_bessel_Jn(1,x);
 //   return j1(x);
 //   return std::sph_bessel ( 1,x);
    break;
  default:
    return gsl_sf_bessel_Jn(n,x);
//    return jn(n, x);
//    return std::sph_bessel ( n,x);
  }

  return -123456;
}






/*#include "special_functions.hpp"

double chebyshev::besselJ(const int n, const double x) {
  
  switch (n) {
  case 0:
    return j0(x);
//    return std::sph_bessel ( 0,x);
    break;
  case 1:
    return j1(x);
 //   return std::sph_bessel ( 1,x);
    break;
  default:
    return jn(n, x);
//    return std::sph_bessel ( n,x);
  }

  return -123456;
}*/
