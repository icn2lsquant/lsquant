#include "chebyshev_moments.hpp"


void chebyshev::Moments2D::MomentNumber(const int mom0, const int mom1 )
{ 
	assert ( mom0<= numMoms[0] && mom1 <= numMoms[1] );
	chebyshev::Moments2D new_mom( mom0, mom1 );
	for( int m0 = 0 ; m0 < mom0 ; m0++)
	for( int m1 = 0 ; m1 < mom1 ; m1++)
		new_mom(m0,m1) = this->operator()(m0,m1); 
	this->numMoms= new_mom.MomentNumber();
	this->MomentVector( new_mom.MomentVector() );
};


chebyshev::Moments2D::Moments2D( std::string momfilename )
{
	//Check if the input_momfile have the right extension 
	std::size_t ext_pos = string( momfilename ).find(".chebmom2D"); 
	if( ext_pos == string::npos )
	{ std::cerr<<"The first argument does not seem to be a valid .chebmom2D file"<<std::endl; assert(false);}

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
	
	momfile>>this->numMoms[0]>>this->numMoms[1];

	this->MomentVector( Moments::vector_t(numMoms[1]*numMoms[0], 0.0) );
	double rmu,imu;
	for( int m0 = 0 ; m0 < numMoms[0] ; m0++)
	for( int m1 = 0 ; m1 < numMoms[1] ; m1++)
	{ 
		momfile>>rmu>>imu;
		this->operator()(m0,m1) = Moments::value_t(rmu,imu);
	}
	momfile.close();
};



void chebyshev::Moments2D::saveIn(std::string filename)
{

    typedef std::numeric_limits<double> dbl;
	ofstream outputfile(filename.c_str());
    outputfile.precision(dbl::digits10);
    outputfile << this->SystemSize() << " " << this->BandWidth() << "  " << this->BandCenter() << std::endl;
    //Print the number of moments for all directions in a line
    for ( auto x : numMoms )
        outputfile << x << " ";
    outputfile << std::endl;

    for ( auto mom : this->MomentVector() )
        outputfile << mom.real() << " " << mom.imag() << std::endl;
    outputfile.close();
};


void chebyshev::Moments2D::Print()
{
	std::cout<<"\n\nCHEBYSHEV 2D MOMENTS INFO"<<std::endl;
	std::cout<<"\tSYSTEM:\t\t\t"<<this->SystemLabel()<<std::endl;
	if( this-> SystemSize() > 0 )
		std::cout<<"\tSIZE:\t\t\t"<<this-> SystemSize()<<std::endl;

	std::cout<<"\tMOMENTS SIZE:\t\t"<<"("<<this->HighestMomentNumber(0)<<" x " <<this->HighestMomentNumber(1)<<")"<<std::endl;
	std::cout<<"\tSCALE FACTOR:\t\t"<<this->ScaleFactor()<<std::endl;
	std::cout<<"\tSHIFT FACTOR:\t\t"<<this->ShiftFactor()<<std::endl;
	std::cout<<"\tENERGY SPECTRUM:\t("
			 <<-this->HalfWidth()+this->BandCenter()<<" , "
			 << this->HalfWidth()+this->BandCenter()<<")"<<std::endl<<std::endl;

};

void chebyshev::Moments2D::ApplyJacksonKernel( const double b0, const double b1 )
{
	assert( b0 >0 && b1>0);
	const double eta0   =  2.0*b0/1000/this->BandWidth();
	const double eta1   =  2.0*b1/1000/this->BandWidth();
		
	int maxMom0=  ceil(M_PI/eta0);
	int maxMom1=  ceil(M_PI/eta1);
	//Removing the negatives to control the overflow
	if(  maxMom0 > numMoms[0] || maxMom0<=0 ) maxMom0 = numMoms[0];
	if(  maxMom1 > numMoms[1] || maxMom1<=0) maxMom1 = numMoms[1];
	std::cout<<"Kernel reduced the number of moments to "<<maxMom0<<" "<<maxMom1<<std::endl;
	this->MomentNumber( maxMom0,maxMom1 ) ;



	const double
	phi_J0 = M_PI/(double)(numMoms[0]+1.0),
	phi_J1 = M_PI/(double)(numMoms[1]+1.0);
		
	double g_D_m0,g_D_m1;
	for( int m0 = 0 ; m0 < numMoms[0] ; m0++)
	{
		g_D_m0=( (numMoms[0]-m0+1)*cos( phi_J0*m0 )+ sin(phi_J0*m0)/tan(phi_J0) )*phi_J0/M_PI;
		for( int m1 = 0 ; m1 < numMoms[1] ; m1++)
		{
			g_D_m1=( (numMoms[1]-m1+1)*cos( phi_J1*m1 )+ sin(phi_J1*m1)/tan(phi_J1) )*phi_J1/M_PI;
			this->operator()(m0,m1) *= g_D_m0*g_D_m1;
		}
	}
}


