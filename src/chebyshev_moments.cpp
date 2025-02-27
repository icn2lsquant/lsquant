#include "chebyshev_moments.hpp"



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
	OP.Multiply( 1.0, this->Chebyshev1(), 0.0, this->Chebyshev0() );
	this->Hamiltonian().Multiply( 1.0, this->Chebyshev0(), 0.0, this->Chebyshev1() );

	return ;
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
	this->Hamiltonian().Multiply( 1.0, this->Chebyshev0(), 0.0, this->Chebyshev1() );
};



void chebyshev::Moments1D_nonOrth::SetInitVectors_nonOrthogonal( Moments::vector_t& T0 )
{
  

	assert( T0.size() == this->SystemSize() );
	const auto dim = this->SystemSize();

	if( this->Chebyshev0().size()!= dim )
		this->Chebyshev0() = Moments::vector_t(dim,Moments::value_t(0)); 

	if( this->Chebyshev1().size()!= dim )
		this->Chebyshev1() = Moments::vector_t(dim,Moments::value_t(0)); 
	//From now on this-> will be discarded in Chebyshev0() and Chebyshev1()

	
        vector_t tmp_(dim,0.0);


	
	linalg::copy ( T0, this->Chebyshev0() );	
	this->Hamiltonian().Multiply( 1.0, this->Chebyshev0(), 0.0, tmp_ );
        linalg::orthogonalize(*S_,  tmp_, this->Chebyshev1());

	
	double b =  this->ShiftFactor();
	
 #pragma omp parallel for
	for(std::size_t i=0; i< dim;i++)	  
	  Chebyshev1()[i] += b * Chebyshev0()[i];

};


void chebyshev::Moments1D_nonOrth::SetInitVectors_nonOrthogonal( SparseMatrixType &OP , Moments::vector_t& T0 )
{
	const auto dim = this->SystemSize();
	assert( OP.rank() == this->SystemSize() && T0.size() == this->SystemSize() );

	if( this->Chebyshev0().size()!= dim )
		this->Chebyshev0() = Moments::vector_t(dim,Moments::value_t(0)); 

	if( this->Chebyshev1().size()!= dim )
		this->Chebyshev1() = Moments::vector_t(dim,Moments::value_t(0)); 
	//From now on this-> will be discarded in Chebyshev0() and Chebyshev1()


	vector_t tmp_(dim,0.0);



	linalg::copy ( T0, this->Chebyshev1() );
        OP.Multiply( 1.0, this->Chebyshev1(), 0.0, this->Chebyshev0() );

	
	this->Hamiltonian().Multiply( 1.0, this->Chebyshev0(), 0.0, tmp_ );
        linalg::orthogonalize(*S_, tmp_ , this->Chebyshev1());


	double b =  this->ShiftFactor();
	
 #pragma omp parallel for
	for(std::size_t i=0; i< dim;i++)	  
	  Chebyshev1()[i] += b * Chebyshev0()[i];
	




	return ;
};


void chebyshev::Moments1D_nonOrth::set_Preconditioner(){
  /* -------------------This isnt working with long int eigen matrices??----------------------// 


  Eigen::IncompleteCholesky<std::complex<double>, Eigen::Lower, Eigen::NaturalOrdering<int>> ichol(S_->eigen_matrix());
  std::cout <<"Incomplete Cholesky solver status:  "<< ichol.info() << std::endl;
  Eigen::SparseMatrix<std::complex<double>> L = ichol.matrixL();
  Eigen::VectorXcd Sichol = ichol.scalingS();
  Eigen::MatrixXcd D = Sichol.asDiagonal().inverse();
  Eigen::SparseMatrix<std::complex<double>> Dinv = D.sparseView();
  Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> P = ichol.permutationP();

  //Eigen::SparseMatrix<std::complex<double>> sparse_L_chol = Dinv * L ;
  //Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor> sparse_S_approx = sparse_L_chol * sparse_L_chol.adjoint() ;


  Sin_sparse_approx_.set_eigen_matrix(sparse_S_approx);
  */

  /*
  //TEST
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<std::complex<double>>> inverter;
  
  Eigen::SparseMatrix<std::complex<double>> I(S_->eigen_matrix().rows(),S_->eigen_matrix().rows());
  I.setIdentity();
  inverter.compute(S_->eigen_matrix());
  Eigen::SparseMatrix<std::complex<double>> inverse_matrix = inverter.solve(I);

  std::cout<<sparse_S_approx.block(0,0,10,10)<<std::endl<<std::endl<<S_->eigen_matrix().block(0,0,10,10)<<std::endl<<std::endl;
  

  
  // -------------------This isnt working with long int eigen matrices??----------------------// 
  */
};



double chebyshev::Moments1D_nonOrth::Iterate_nonOrthogonal_test( SparseMatrixType &orth_Ham )
{

        double error=0;  
        std::size_t dim = Chebyshev0().size();
        vector_t tmp_( dim, 0.0 ),
	  tmp_2_( dim, 0.0),
	  compare( dim, 0.0 );

	double b =  this->ShiftFactor();

	
	this->Hamiltonian().Multiply( this->Chebyshev1(),tmp_);
       	linalg::orthogonalize(*S_, tmp_, tmp_2_);

	
	orth_Ham.Multiply(this->Chebyshev1(), compare);
	
	
	#pragma omp parallel for
	for(std::size_t i=0; i< dim;i++){
	  std::complex<double> tmp_2 = tmp_2_[i] + b * Chebyshev1()[i];
	  
	  Chebyshev0()[i] = 2.0 * tmp_2 - Chebyshev0()[i];

	  error += std::abs( ( compare[i] - tmp_2 ) / compare[i] );
	}
	
	this->Chebyshev0().swap(this->Chebyshev1());

	

	return error/dim;
};


int chebyshev::Moments1D_nonOrth::Iterate_nonOrthogonal( )
{

        std::size_t dim = Chebyshev0().size();
        vector_t tmp_( dim, 0.0 ),
	         tmp_2_( dim, 0.0);


	
	this->Hamiltonian().Multiply( this->Chebyshev1(),tmp_);
	
       	linalg::orthogonalize(*S_, tmp_, tmp_2_);
	
	#pragma omp parallel for
	for(std::size_t i=0; i< dim;i++){
	  Chebyshev0()[i] = 2.0 * tmp_2_[i] - Chebyshev0()[i];
	}
	
	this->Chebyshev0().swap(this->Chebyshev1());

	

	return 0;
};



int chebyshev::Moments::Iterate( )
{
	this->Hamiltonian().Multiply(2.0,this->Chebyshev1(),-1.0,this->Chebyshev0());
	this->Chebyshev0().swap(this->Chebyshev1());
	return 0;
};
