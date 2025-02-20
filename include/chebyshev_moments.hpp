// With contributions made by Angel D. Prieto S.
#ifndef CHEBYSHEV_MOMENTS
#define CHEBYSHEV_MOMENTS


#include <complex>
#include <vector>
#include <string>
#include <array>

#include "sparse_matrix.hpp" //contain SparseMatrixType
#include <cassert>			 //needed for assert
#include <fstream>   		 //For ifstream and ofstream
#include <limits>    		 //Needed for dbl::digits10
#include "linear_algebra.hpp"
#include "vector_list.hpp"
#include "special_functions.hpp"
#include "chebyshev_coefficients.hpp"

namespace chebyshev 
{
	const double CUTOFF = 0.99;

class Moments
{
	public:
	typedef std::complex<double>  value_t;
	typedef std::vector< value_t > vector_t;

	//default constructor
	Moments():
	_pNHAM(0),system_label(""),system_size(0),
	band_width(0),band_center(0){};

	//GETTERS
	void getMomentsParams( Moments& mom)
	{
		this->SetHamiltonian( mom.Hamiltonian() ) ; 
		this->SystemLabel( mom.SystemLabel());
		this->BandWidth( mom.BandWidth() );
		this->BandCenter( mom.BandCenter() );
	};	
	
	inline
	size_t SystemSize() const { return system_size; };

	inline
	string SystemLabel() const { return system_label; };

	inline
	double BandWidth() const { return band_width; };

	inline
	double HalfWidth() const { return BandWidth()/2.0; };

	inline
	double BandCenter() const { return band_center; };

	inline
	double ScaleFactor() const { return chebyshev::CUTOFF/HalfWidth(); };

	inline
	double ShiftFactor() const { return -BandCenter()/HalfWidth()*chebyshev::CUTOFF; };

	inline 
	vector_t& MomentVector() { return mu ;}

	inline
	value_t& MomentVector(const size_t i){return  mu[i]; };

	inline
	Moments::vector_t& Chebyshev0(){ return ChebV0; } 

	inline
	Moments::vector_t& Chebyshev1(){ return ChebV1; } 

	
	inline
	SparseMatrixType& Hamiltonian()
	{ 
		return *_pNHAM; 
	};



	//SETTERS
	inline
	void SetHamiltonian( SparseMatrixType& NHAM )
	{ 
		if ( this->SystemSize() == 0 ) //Use the rank of the hamiltonian as system size
			this->SystemSize( NHAM.rank() );	
		assert( NHAM.rank() == this->SystemSize()  );
		_pNHAM = &NHAM; 
	};


	//Heavy functions
	int  Rescale2ChebyshevDomain()
	{
		this->Hamiltonian().Rescale(this->ScaleFactor(),this->ShiftFactor());
		return 0;
	};


	inline
	void SetAndRescaleHamiltonian(SparseMatrixType& NHAM)
	{ 
		this->SetHamiltonian(NHAM );
		this->Rescale2ChebyshevDomain();
	};


	inline
	void SystemSize(const int dim)  { system_size = dim; };

	inline
	void SystemLabel(string label)  { system_label = label; };

	inline
	void BandWidth( const double x)  { band_width = x; };

	inline
	void BandCenter(const double x) { band_center = x; };

	inline 
	void MomentVector(const vector_t _mu ) { mu= _mu;}

	void SetInitVectors( const vector_t& T0 );

	void SetInitVectors( SparseMatrixType &OP ,const vector_t& T0 );


	inline
	double Rescale2ChebyshevDomain(const double energ)
	{ 
		return (energ - this->BandCenter() )/this->HalfWidth(); 
	};

	int Iterate( );

	//light functions
    int JacksonKernelMomCutOff( const double broad );
	
	//light functions
    double JacksonKernel(const double m,  const double Mom );


	private:
	SparseMatrixType* _pNHAM;

	Moments::vector_t ChebV0,ChebV1,OPV;
	std::string system_label;
	size_t system_size;
	double band_width,band_center;
	vector_t mu;	
};


class Moments1D: public Moments
{
	public: 

	Moments1D():_numMoms(0){};

	Moments1D(const size_t m0):_numMoms(m0){ this->MomentVector( Moments::vector_t(_numMoms, 0.0) ); };

	Moments1D( std::string momfilename );

	//GETTERS
	inline
	size_t MomentNumber() const { return _numMoms; };

	inline
	size_t HighestMomentNumber() const { return  _numMoms; };

	//SETTERS

	//OPERATORS
	inline
	Moments::value_t& operator()(const size_t m0) { return this->MomentVector(m0); };

	void MomentNumber(const size_t _numMoms );

	void saveIn(std::string filename);

	//Transformation
	void ApplyJacksonKernel( const double broad );

	// Input/Output
	void Print();

	private:
	size_t _numMoms;
};


class Moments2D: public Moments
{
	public:
	Moments2D():numMoms({0,0}){};


	Moments2D(const size_t m0,const size_t m1):numMoms({m0,m1}){ this->MomentVector( Moments::vector_t(numMoms[1]*numMoms[0], 0.0) );    };

	Moments2D( std::string momfilename );


	Moments2D( Moments2D mom, const size_t m0,const size_t m1 )
	{
		this->getMomentsParams(mom);
		this->numMoms={m0,m1};
		this->MomentVector( Moments::vector_t(numMoms[1]*numMoms[0], 0.0) );    
	};

	//GETTERS

	array<int, 2> MomentNumber() const { return numMoms; };

	int HighestMomentNumber(const int i) const { return numMoms[i]; };

	inline
	int HighestMomentNumber() const { return  (numMoms[1] > numMoms[0]) ? numMoms[1] : numMoms[0]; };

	//SETTERS
	void MomentNumber(const int mom0, const int mom1 );

	//OPERATORS
	inline
	Moments::value_t& operator()(const int m0,const int m1) { return this->MomentVector(m0*numMoms[1] + m1 ); };

	//Transformation
	void ApplyJacksonKernel( const double b0, const double b1 );

	// Input/Output 
	//COSTFUL FUNCTIONS
	void saveIn(std::string filename);

	
	void AddSubMatrix( Moments2D& sub , const int mL, const int mR)
	{
		for(int m0=0; m0<sub.HighestMomentNumber(0); m0++)
		for(int m1=0; m1<sub.HighestMomentNumber(1); m1++)
			this->operator()(mL+m0,mR+m1) += sub(m0,m1);
	} 

	void InsertSubMatrix( Moments2D& sub , const int mL, const int mR)
	{
		for(int m0=0; m0<sub.HighestMomentNumber(0); m0++)
		for(int m1=0; m1<sub.HighestMomentNumber(1); m1++)
			this->operator()(mL+m0,mR+m1) = sub(m0,m1);
	} 
	
	

	void Print();

	protected:
	array<int, 2> numMoms;
};

  class MomentsTD : public Moments
  {
	  private:
		SparseMatrixType* _pNCON;
	  public:
	  const double HBAR = 0.6582119624 ;//planck constant in eV.fs
	 


	  MomentsTD():
	  _numMoms(1), _maxTimeStep(1),
	  _timeStep(0),
	  _dt(0)
	  {};
	  
	  MomentsTD( const size_t m, const size_t n ): 
	  _numMoms(m), _maxTimeStep(n),
	  _timeStep(0),
	  _dt(0)
	  { this->MomentVector( Moments::vector_t(m*n, 0.0) );    };
	  
	  MomentsTD( std::string momfilename );
	  
	  //GETTERS
	  inline
	  size_t MomentNumber() const { return _numMoms;};
	  
	  inline
	  size_t HighestMomentNumber() const { return _numMoms;};
	  
	  inline
	  size_t CurrentTimeStep() const { return _timeStep;};

	  inline 
	  size_t MaxTimeStep() const { return _maxTimeStep; };
	  
	  inline
	  double TimeDiff() const   { return _dt; };
	  	  
	  inline
	  double ChebyshevFreq() const   { return HalfWidth()/chebyshev::CUTOFF/HBAR; };

	  inline
	  double ChebyshevFreq_0() const   { return BandCenter()/HBAR; };
	  
	  inline
	  Moments::vector_t& Chebyshevx0(){ return ChebVx0; } 

	  inline
	  Moments::vector_t& Chebyshevx1(){ return ChebVx1; } 

	  inline
	  Moments::vector_t& DeltaPhi(){ return DelPhi; } 

	  inline
	  Moments::vector_t& MSDWF(){ return WFGaussian; } 

	  void SetConmuter( SparseMatrixType& NCON )
	  { 
		if ( this->SystemSize() == 0 ) //Use the rank of the hamiltonian as system size
			this->SystemSize( NCON.rank() );	
		assert( NCON.rank() == this->SystemSize()  );
		_pNCON = &NCON; 
	  };

	  int  RescaleCon2ChebyshevDomain()
	  {
		const auto I = Moments::value_t(0, 1);
		this->Conmuter().Rescale(this->ScaleFactor(),this->ShiftFactor());
		std::cout<<"Warning only works for symmetric boundaries in energy"<<std::endl;
		return 0;
          };

	  int ResetGaussian()
	  {
		const auto dim = this->SystemSize();
		if( this->MSDWF().size()!= dim )
			this->MSDWF() = Moments::vector_t(dim,Moments::value_t(0)); 
		linalg::scal(0.0,this->MSDWF());
		return 0;
	  };
 
	  void SetAndRescaleConmuter(SparseMatrixType& NCON)
	  { 
		this->SetConmuter(NCON );
		this->RescaleCon2ChebyshevDomain();
	  };



	  SparseMatrixType& Conmuter()
	  { 
		return *_pNCON; 
  	  };



	  //SETTERS
	
	  void resetCurrentTm(){
		  const auto dim = this->SystemSize();
		  if( this->DeltaPhi().size()!= dim )
			  this->DeltaPhi() = Moments::vector_t(dim,Moments::value_t(0)); 
		  currentTm=0;
	  }

	  void MomentNumber(const size_t mom);
	  
	  void MaxTimeStep(const  size_t maxTimeStep )  {  _maxTimeStep = maxTimeStep; };

	  inline
	  void IncreaseTimeStep(){ _timeStep++; };

	  inline
	  void ResetTime(){ _timeStep=0; };
	  
	  inline
	  void TimeDiff(const double dt ) { _dt = dt; };
	  
	  
	  int Evolve(  vector_t& Phi);

	  int Conv_Evolve(  vector_t& Phi);

	  int Deltaiter(double E);

	  int MSD_Evolve(  vector_t& Phi);

	  //OPERATORS
	  inline
	  Moments::value_t& operator()(const size_t m, const size_t n)
	  {
		  return this->MomentVector( m*MaxTimeStep() + n );
	  };
	  
	  //Transformation
	  void ApplyJacksonKernel( const double broad );
	  
	  //COSTFUL FUNCTIONS
	  void saveIn(std::string filename);
	  
	  void Print();

  private:
    int currentTm;
    Moments::vector_t ChebVx0,ChebVx1,WFGaussian,DelPhi;
    size_t _numMoms, _maxTimeStep, _timeStep;
    double _dt;
  };









  class MomentsLocal : public Moments
  {
	  private:
		SparseMatrixType* _pNCON;
	  public:
	  const double HBAR = 0.6582119624 ;//planck constant in eV.fs
	 


	  MomentsLocal():
	  _numMoms(1), _Norb(1)
	  {};
	  
	  MomentsLocal( const size_t m, const size_t n ): 
	  _numMoms(m), _Norb(n)
	  { this->MomentVector( Moments::vector_t(m*n, 0.0) );    };
	  
//	  MomentsLocal( std::string momfilename );

	  MomentsLocal( std::string momfilename,bool binflag);


	  //GETTERS
	  inline
	  size_t MomentNumber() const { return _numMoms;};
	  
	  inline
	  size_t HighestMomentNumber() const { return _numMoms;};
	  
	  inline
	  size_t CurrentOrbital() const { return _currentOrbital;};

	  inline 
	  size_t NumberOfOrbitals() const { return _Norb; };
	  
	  inline
	  double ChebyshevFreq() const   { return HalfWidth()/chebyshev::CUTOFF/HBAR; };

	  inline
	  double ChebyshevFreq_0() const   { return BandCenter()/HBAR; };
	  
	  inline
	  Moments::vector_t& DeltaPhi(){ return DelPhi; } 





	  //SETTERS
	
	  void resetCurrentTm(){
		  const auto dim = this->SystemSize();
		  if( this->DeltaPhi().size()!= dim )
			  this->DeltaPhi() = Moments::vector_t(dim,Moments::value_t(0)); 
		  currentTm=0;
	  }

	  void MomentNumber(const size_t mom);
	  
	  void NumberOfOrbitals(const  size_t Norb )  {  _Norb = Norb; };

	  inline
	  void IncreaseCurrentOrbital(){ _currentOrbital++; };

	  inline
	  void ResetOrbital(){ _currentOrbital=0; };
	  
	  
	  int Deltaiter(double E);

	  //OPERATORS
	  inline
	  Moments::value_t& operator()(const size_t m, const size_t n)
	  {	  //size_t otherm=1073;
//	  	  for (size_t i=0;i<NumberOfOrbitals()*HighestMomentNumber();i++){
//			  size_t idx=otherm*NumberOfOrbitals()+i;
//			  std::cout<<idx<<std::endl;
//			  std::cout<<" "<<i<<" "<<idx<<" "<<MomentVector(idx)<<std::endl;
//		  }	  
//		  if (m>1072){
//		  std::cout<<m*NumberOfOrbitals()+n<<std::endl;
//		  }
		  return this->MomentVector( m*NumberOfOrbitals() + n );
	  };
	  
	  //Transformation
	  void ApplyJacksonKernel( const double broad );
	  
	  //COSTFUL FUNCTIONS
	  void saveIn(std::string filename);

	  void saveInBin(std::string filename);


	  void Print();

  private:
    int currentTm;
    Moments::vector_t WFGaussian,DelPhi;
    size_t _numMoms, _Norb, _currentOrbital;
  };








class MomentsAC: public Moments2D
{
 public:
	MomentsAC( std::string momfilename );
	//SETTERS
	bool Freq(std::string frequenciesfilename);
	void Freq(vector<double> freq_){freq=freq_;}
	void HighestFreq(double freq){highestFreq=freq;}
	void NumFreq(int numFreq_){numFreq=numFreq_;}
	//GETTERS
	double CurrentFreq(){return freq[currentFreq];}
	int NumFreq(){return numFreq;}
	double HighestFreq(){return highestFreq;}
	
	
	double IterateFreq(){return currentFreq++;}
	void RescaleFreq(){for (int i=0;i<this->NumFreq();i++){freq[i]*=this->ScaleFactor();}}

 private:
	std::vector< double > freq;
	int currentFreq,numFreq;
	double highestFreq;
};



























class Vectors : public Moments
{
	public: 
	typedef VectorList< Moments::value_t > vectorList_t;


	int NumberOfVectors() const
	{ return numVecs;}


	int SetNumberOfVectors( const int x)
	{ numVecs = x; return 0;}


	Vectors():Chebmu(0,0){};
	
	Vectors(const size_t nMoms,const size_t dim ):Chebmu(nMoms,dim) {  };

	Vectors( Moments1D& mom ): Chebmu(mom.HighestMomentNumber(), mom.SystemSize() )
	{ 
		this->getMomentsParams(mom);
	};
	
	Vectors( Moments2D& mom ): Chebmu(mom.HighestMomentNumber(), mom.SystemSize() )
	{ 
		this->getMomentsParams(mom);
	};

	
	Vectors( Moments2D& mom, size_t i  ): Chebmu(mom.HighestMomentNumber(), mom.SystemSize() )
	{ 
		this->getMomentsParams(mom);
	};
	
	
   
    int CreateVectorSet()
    {
		try
		{
			const int vec_size  = this->SystemSize();
			const int list_size = this->NumberOfVectors();
			Chebmu(list_size, vec_size );
		}
		catch (...)
		{ std::cerr<<"Failed to initilize the vector list."<<std::endl;}

		
		return 0;
	}
    
    
    
    Vectors( MomentsTD& mom ): Chebmu(mom.HighestMomentNumber(), mom.SystemSize() )
	{ 
		this->getMomentsParams(mom);
	};

	void getMomentsParams( Moments& mom)
	{
		this->SetHamiltonian( mom.Hamiltonian() ) ; 
		this->SystemLabel( mom.SystemLabel());
		this->BandWidth( mom.BandWidth() );
		this->BandCenter( mom.BandCenter() );
		if ( mom.Chebyshev0().size() == this->SystemSize() )
			this->SetInitVectors( mom.Chebyshev0() );
	}

	inline
	size_t Size() const
	{
		return  (long unsigned int)this->SystemSize()*
				(long unsigned int)this->NumberOfVectors();
	}

	inline 
	double SizeInGB() const
	{
		const double vec_size  = this->SystemSize();
		const double list_size = this->NumberOfVectors();
		return  sizeof(value_t)*vec_size*list_size*pow(2.0,-30.0);
	}

	inline
	size_t HighestMomentNumber() const { return  this->Chebmu.ListSize(); };


	inline
	vectorList_t& List() { return this->Chebmu; };

	inline
	Moments::vector_t& Vector(const size_t m0) { return this->Chebmu.ListElem(m0); };


	inline
	Moments::value_t& operator()(const size_t m0) { return this->Chebmu(m0,0); };

	int IterateAll( );

	int EvolveAll( const double DeltaT, const double Omega0);

	int Multiply( SparseMatrixType &OP );


	double MemoryConsumptionInGB();


	private:
	Moments::vector_t OPV;
	vectorList_t Chebmu;	
	int numVecs;
};

};

#endif 


