#include "chebyshev_moments.hpp"

int chebyshev::MomentsLocal::Deltaiter(double E)
{

	const auto NumMoms = this->HighestMomentNumber();


	/*if( this->Chebyshev0().size()!= dim )
		std::cout<<"warning"<<std::endl;
		this->Chebyshev0() = Moments::vector_t(dim,Moments::value_t(0)); 

	if( this->Chebyshev1().size()!= dim )
		std::cout<<"warning"<<std::endl;
		this->Chebyshev1() = Moments::vector_t(dim,Moments::value_t(0)); 
*/


	//From now on this-> will be discarded in Chebyshev0() and Chebyshev1()


  	double g_D_m=chebyshev::Moments::JacksonKernel(currentTm,  NumMoms );
	double TE=cos(currentTm*acos(E));
	if (currentTm==0){
//  		for ( auto wfcoef : Chebyshev0())
//    			std::cout << wfcoef.real() << " " << wfcoef.imag() << std::endl;
		std::cout<<std::endl<<"First broadening term in step "<<currentTm<<" is "<<TE*delta_chebF(E,currentTm)*g_D_m<<std::endl; 
		linalg::axpy( TE*delta_chebF(E,currentTm)*g_D_m , Chebyshev0(), DeltaPhi());
		currentTm++;
		return 0;
	}

	else if (currentTm==1){
		linalg::axpy( TE*delta_chebF(E,currentTm)*g_D_m , Chebyshev1(), DeltaPhi());
		currentTm++;
		return 0;
	}
	else{
		this->Hamiltonian().Multiply(2.0,  Chebyshev1(), -1.0,  Chebyshev0() );
		Chebyshev0().swap(Chebyshev1());
		linalg::axpy( TE*delta_chebF(E,currentTm)*g_D_m , Chebyshev1(), DeltaPhi());
		currentTm++;
		return 0;
	}
};




chebyshev::MomentsLocal::MomentsLocal( std::string momfilename,bool binflag)
{
if (!binflag){
  //Check if the input_momfile have the right extension 
  std::size_t ext_pos = string( momfilename ).find(".chebmomLocal"); 
  if( ext_pos == string::npos )
    { std::cerr<<"The first argument does not seem to be a valid .chebmomLocal file"<<std::endl; assert(false);}

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
  momfile>>ibuff; this->NumberOfOrbitals(ibuff);

  //create the moment array and read the data
	
  momfile>>this->_numMoms>>this->_Norb;

  this->MomentVector( Moments::vector_t( this->HighestMomentNumber() * this->NumberOfOrbitals(), 0.0) );
  double rmu, imu;
  for( int m = 0; m < this->HighestMomentNumber() ; m++)
  for( int n = 0; n < this->NumberOfOrbitals() ; n++)
  { 
	  momfile>>rmu>>imu;
	  this->operator()(m,n) = Moments::value_t(rmu, imu);
  }
  momfile.close();
}
else{
  //Made to read from binary files
  //Check if the input_momfile have the right extension 
  std:: cout <<"Asking the program to read from binary file"<<std::endl;
  std::size_t ext_pos = string( momfilename ).find(".chebmomLocal.bin"); 
  if( ext_pos == string::npos )
    { std::cerr<<"The first argument does not seem to be a valid .chebmomLocal.bin binary file"<<std::endl; assert(false);}

  //if it does, use it to get the extension
  this->SystemLabel( momfilename.substr(0,ext_pos) ); 

  //and then try to open the file


  


  FILE * momfile = fopen(momfilename.c_str(), "rb");
  if (momfile!=nullptr){
  	size_t tbuff;
	double dbuff;
	fread(&tbuff,sizeof(size_t),1,momfile); this->SystemSize(tbuff);
	fread(&dbuff,sizeof(double),1,momfile); this->BandWidth(dbuff);
	fread(&dbuff,sizeof(double),1,momfile); this->BandCenter(dbuff);
	fread(&tbuff,sizeof(size_t),1,momfile); this->NumberOfOrbitals(tbuff);
	fread(&tbuff,sizeof(size_t),1,momfile); _numMoms=tbuff;
	fread(&tbuff,sizeof(size_t),1,momfile); _Norb=tbuff;

        size_t nvalues=this->NumberOfOrbitals()*this->HighestMomentNumber();
  	this->MomentVector( Moments::vector_t( nvalues, 0.0) );
	if(MomentVector().data()!=nullptr){
		fread(MomentVector().data(),sizeof(value_t),nvalues,momfile);
	}else{
		std::cout<<"Unable to allocate MomentVector"<<std::endl;
	}

  fclose(momfile);

  }else{std::cout<<"Unable to open file"<<std::endl;}
}

};



















void chebyshev::MomentsLocal::saveIn(std::string filename)
{
  typedef std::numeric_limits<double> dbl;
  ofstream outputfile(filename.c_str());
  outputfile.precision(dbl::digits10);
  outputfile << this->SystemSize() << " " << this->BandWidth() << " " << this->BandCenter() << " " 
	     << this->NumberOfOrbitals() << " " << std::endl;
  //Print the number of moments for all directions in a line
  outputfile << this->HighestMomentNumber() << " " << this->NumberOfOrbitals() << " ";

  for ( auto mom : this->MomentVector() )
    outputfile << mom.real() << " " << mom.imag() << std::endl;
  outputfile.close();
};


void chebyshev::MomentsLocal::saveInBin(std::string filename)
{
  std::cout<<"Storing as binary file"<<std::endl;
  FILE * outputfile = fopen(filename.c_str(), "wb");
  size_t nvalues=this->HighestMomentNumber()*this->NumberOfOrbitals();
  size_t tbuff;
  double dbuff;
  if (outputfile!=nullptr){
	tbuff=this->SystemSize();fwrite(&tbuff,sizeof(size_t),1,outputfile);
	dbuff=this->BandWidth();fwrite(&dbuff,sizeof(double),1,outputfile);
	dbuff=this->BandCenter();fwrite(&dbuff,sizeof(double),1,outputfile);
	tbuff=this->NumberOfOrbitals();fwrite(&tbuff,sizeof(size_t),1,outputfile);
	tbuff=this->HighestMomentNumber();fwrite(&tbuff,sizeof(size_t),1,outputfile);
	tbuff=this->NumberOfOrbitals();fwrite(&tbuff,sizeof(size_t),1,outputfile);
	fwrite(this->MomentVector().data(),sizeof(value_t),nvalues,outputfile);
	fclose(outputfile);
  }
};









void chebyshev::MomentsLocal::Print()
{
  std::cout<<"\n\nCHEBYSHEV Local MOMENTS INFO"<<std::endl;
  std::cout<<"\tSYSTEM:\t\t\t"<<this->SystemLabel()<<std::endl;
  if( this-> SystemSize() > 0 )
    std::cout<<"\tSIZE:\t\t\t"<<this-> SystemSize()<<std::endl;

  std::cout<<"\tMOMENTS SIZE:\t\t"<<"("
	   <<this->HighestMomentNumber()<< " x " <<this->NumberOfOrbitals()<<")"<<std::endl;
  std::cout<<"\tSCALE FACTOR:\t\t"<<this->ScaleFactor()<<std::endl;
  std::cout<<"\tSHIFT FACTOR:\t\t"<<this->ShiftFactor()<<std::endl;
  std::cout<<"\tENERGY SPECTRUM:\t("
	   <<-this->HalfWidth()+this->BandCenter()<<" , "
	   << this->HalfWidth()+this->BandCenter()<<")"<<std::endl<<std::endl;
  std::cout<<"\tNUMBER OF LOCAL ORBITALS:\t\t"<<this->NumberOfOrbitals()<<std::endl;
  
};

void chebyshev::MomentsLocal::MomentNumber(const size_t numMoms )
{ 
	assert( numMoms <= this->HighestMomentNumber() );
	const auto numberoforbitals = this->NumberOfOrbitals();
	
	//Copy all previous moments smaller or equal than the new size into new vector;
	chebyshev::MomentsLocal new_mom( numMoms, numberoforbitals);
	for( size_t m = 0 ; m < numMoms  ; m++)
	for( size_t n = 0 ; n < numberoforbitals ; n++)
      new_mom(m, n) = this->operator()(m, n); 

	this->_numMoms = new_mom._numMoms;
	this->MomentVector( new_mom.MomentVector() );
};

void chebyshev::MomentsLocal::ApplyJacksonKernel( const double broad )
{
  assert( broad > 0);
  const double eta = 2.0*broad/1000/this->BandWidth();
  int maxMom =  ceil(M_PI/eta);
  
  if(  maxMom > this->HighestMomentNumber() )
	maxMom = this->HighestMomentNumber();
  std::cout << "Kernel reduced the number of moments to " << maxMom <<" for a broadening of "<<M_PI/maxMom*this->BandWidth() << "eV" << std::endl;
  this->MomentNumber(maxMom);

  const double phi_J = M_PI/(double)(maxMom+1.0);
  double g_D_m;

  for( size_t m = 0 ; m < maxMom ; m++)
  {
	  g_D_m = ( (maxMom - m + 1) * cos(phi_J * m) + sin(phi_J * m) /tan(phi_J) ) * phi_J/M_PI;
	  for( size_t n = 0 ; n < this->NumberOfOrbitals() ; n++) this->operator()(m, n) *= g_D_m;
  }
}
