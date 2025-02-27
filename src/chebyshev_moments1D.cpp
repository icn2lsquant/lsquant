#include "chebyshev_moments.hpp"

void chebyshev::Moments1D::Print()
{
	std::cout<<"\n\nCHEBYSHEV 1D MOMENTS INFO"<<std::endl;
	std::cout<<"\tSYSTEM:\t\t\t"<<this->SystemLabel()<<std::endl;
	if( this-> SystemSize() > 0 )
		std::cout<<"\tSIZE:\t\t\t"<<this-> SystemSize()<<std::endl;

	std::cout<<"\tMOMENTS SIZE:\t\t"<<this->HighestMomentNumber()<<std::endl;
	std::cout<<"\tSCALE FACTOR:\t\t"<<this->ScaleFactor()<<std::endl;
	std::cout<<"\tSHIFT FACTOR:\t\t"<<this->ShiftFactor()<<std::endl;
	std::cout<<"\tENERGY SPECTRUM:\t("
			 <<-this->HalfWidth()+this->BandCenter()<<" , "
			 << this->HalfWidth()+this->BandCenter()<<")"<<std::endl<<std::endl;

};

chebyshev::Moments1D::Moments1D( std::string momfilename )
{
	//Check if the input_momfile have the right extension 
	std::size_t ext_pos = string( momfilename ).find(".chebmom1D"); 
	if( ext_pos == string::npos )
	{ std::cerr<<"The first argument does not seem to be a valid .chebmom1D file"<<std::endl; assert(false);}

	//if it does, use it to get the extension
	this->SystemLabel( momfilename.substr(0,ext_pos) ); 

	//and then try to open the file
	std::ifstream momfile( momfilename.c_str() );
	assert( momfile.is_open() );

	//if succesful, read the header
	int ibuff; double dbuff;
	momfile>>ibuff; this->SystemSize(ibuff);
	momfile>>dbuff; this->BandWidth(dbuff);
	momfile>>dbuff; this->BandCenter(dbuff);

	//create the moment array and read the data
	
	momfile>>this->_numMoms;

	this->MomentVector( Moments::vector_t(_numMoms, 0.0) );
	double rmu,imu;
	for( int m0 = 0 ; m0 < _numMoms ; m0++)
	{ 
		momfile>>rmu>>imu;
		this->operator()(m0) = Moments::value_t(rmu,imu);
	}
	momfile.close();
};


void chebyshev::Moments1D::MomentNumber(const size_t numMoms )
{ 
	assert( this->_numMoms <= this->HighestMomentNumber() );
	
	//Copy all previous moments smaller or equal than the new size into new vector;
	chebyshev::Moments1D new_mom( numMoms );
	for( size_t m = 0 ; m < numMoms  ; m++)
      new_mom(m) = this->operator()(m); 

	this->_numMoms = new_mom._numMoms;
	this->MomentVector( new_mom.MomentVector() );
};


void chebyshev::Moments1D::ApplyJacksonKernel( const double broad )
{
  assert( broad > 0);
  const double eta = 2.0*broad/1000/this->BandWidth();
  int maxMom =  ceil(M_PI/eta);
  
  if(  maxMom > this->HighestMomentNumber() )
	maxMom = this->HighestMomentNumber();
  std::cout << "Kernel reduced the number of moments to " << maxMom <<" for a broadening of "<<M_PI/maxMom << std::endl;
  this->MomentNumber(maxMom);

  const double phi_J = M_PI/(double)(maxMom+1.0);
  double g_D_m;

  for( size_t m = 0 ; m < maxMom ; m++)
  {
	  g_D_m = ( (maxMom - m + 1) * cos(phi_J * m) + sin(phi_J * m) /tan(phi_J) ) * phi_J/M_PI;
	  this->operator()(m) *= g_D_m;
	}
}

void chebyshev::Moments1D::ApplyLorentzKernel( const double broad, const double lambda )
{
  assert( broad > 0);
  const double eta = 2.0*broad/1000/this->BandWidth();
  int maxMom =  ceil(lambda/eta);// ceil(M_PI/eta);
  
  if(  maxMom > this->HighestMomentNumber() )
	maxMom = this->HighestMomentNumber();
  std::cout << "Kernel reduced the number of moments to " << maxMom <<" for a broadening of "<<lambda/maxMom << std::endl;
  this->MomentNumber(maxMom);


  double g_D_m;

  for( size_t m = 0 ; m < maxMom ; m++)
  {
    g_D_m = std::sinh( lambda * (1.0 - double(m) / double(maxMom)) ) / std::sinh( lambda ); 
    this->operator()(m) *= g_D_m;
  }
}


	
void chebyshev::Moments1D::saveIn(std::string filename)
{
  typedef std::numeric_limits<double> dbl;
  ofstream outputfile(filename.c_str());
  assert(outputfile.is_open());
  outputfile.precision(dbl::digits10);
  outputfile << this->SystemSize() << " " << this->BandWidth() << " " << this->BandCenter()<<std::endl;
  //Print the number of moments for all directions in a line
  outputfile << this->HighestMomentNumber() << std::endl;

  for ( auto mom : this->MomentVector() )
    outputfile << mom.real() << " " << mom.imag() << std::endl;
  outputfile.close();
};


