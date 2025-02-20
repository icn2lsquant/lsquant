#include "chebyshev_moments.hpp"


void chebyshev::Moments::SetInitVectors( const Moments::vector_t& T0 )
{
	assert( T0.size() == this->SystemSize() );
	const auto dim = this->SystemSize();

	if( this->Chebyshev0().size()!= dim )
		this->Chebyshev0() = Moments::vector_t(dim,Moments::value_t(0)); 

	if( this->Chebyshev1().size()!= dim )
		this->Chebyshev1() = Moments::vector_t(dim,Moments::value_t(0)); 
	//From now on this-> will be discarded in Chebyshev0() and Chebyshev1()

	linalg::copy ( T0, this->Chebyshev0() );
	this->Hamiltonian().Multiply( this->Chebyshev0(), this->Chebyshev1() );
};


void chebyshev::Moments::SetInitVectors( SparseMatrixType &OP ,const Moments::vector_t& T0 )
{
	const auto dim = this->SystemSize();
	assert( OP.rank() == this->SystemSize() && T0.size() == this->SystemSize() );

	if( this->Chebyshev0().size()!= dim )
		this->Chebyshev0() = Moments::vector_t(dim,Moments::value_t(0)); 

	if( this->Chebyshev1().size()!= dim )
		this->Chebyshev1() = Moments::vector_t(dim,Moments::value_t(0)); 
	//From now on this-> will be discarded in Chebyshev0() and Chebyshev1()

	linalg::copy ( T0, this->Chebyshev1() );
	OP.Multiply( this->Chebyshev1(), this->Chebyshev0() );
	this->Hamiltonian().Multiply( this->Chebyshev0(), this->Chebyshev1() );

	return ;
};



int chebyshev::Moments::Iterate( )
{
	this->Hamiltonian().Multiply(2.0,this->Chebyshev1(),-1.0,this->Chebyshev0());
	this->Chebyshev0().swap(this->Chebyshev1());
	return 0;
};

	//light functions
int chebyshev::Moments::JacksonKernelMomCutOff( const double broad )
{
	assert( broad >0 );
	const double eta   =  2.0*broad/1000/this->BandWidth();
	return ceil(M_PI/eta);
};
	
//light functions
double chebyshev::Moments::JacksonKernel(const double m,  const double Mom )
{
	const double
	phi_J = M_PI/(double)(Mom+1.0);
	return ( (Mom-m+1)*cos( phi_J*m )+ sin(phi_J*m)/tan(phi_J) )*phi_J/M_PI;
};

