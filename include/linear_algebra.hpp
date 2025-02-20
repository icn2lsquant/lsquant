

#ifndef ALGEBRA_FUNCTION
#define ALGEBRA_FUNCTION
#include <complex>
#include <vector>
#include <cassert>

using namespace std;

#define MKL_Complex16 complex<double>
#include "mkl.h"

namespace linalg
{

// VECTOR-VECTOR FUNCTIONS
void scal(const int dim, complex<double> a, complex<double> *x);

void scal(const complex<double>& a, vector< complex<double> >& x);

void axpy(const int dim, complex<double> a, const complex<double> *x, complex<double> *y);

void axpy(const complex<double>& a, const vector< complex<double> >& x, vector< complex<double> >& y);

void copy(const int dim, const complex<double> *x, complex<double> *y);

void copy(const vector< complex<double> >&x,vector< complex<double> >& y);

complex<double> vdot(const vector< complex<double> >& x,const vector< complex<double> >& y);

complex<double> vdot(const int dim, const complex<double> *x, complex<double> *y);

double nrm2(const vector< complex<double> >& x);

double nrm2(const int dim, const complex<double> *x);


void batch_vdot(const int dim,const int batchSize,const complex<double>* leftb,const complex<double>* rightb,complex<double>* output);

} // namespace linalg

#endif
