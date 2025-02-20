#include "chebyshev_moments.hpp"


chebyshev::MomentsAC::MomentsAC( std::string momfilename )
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


bool chebyshev::MomentsAC::Freq(std::string frequenciesfilename)
{

	std::ifstream freqfile( frequenciesfilename.c_str() );
	assert(freqfile.is_open());

	//read the header
	std::string sbuff;
	freqfile>>sbuff;
	int ibuff;
	double dbuff;
	if (sbuff=="grid")
	{
		double maxFreq=0;
		freqfile>> ibuff;
		this->NumFreq(ibuff);
		std::vector<double> freq_(this->NumFreq());
		for (int i=0;i<this->NumFreq();i++)
		{
			freqfile>>dbuff;
			if (dbuff < 0){std::cout<<"Not allowed negative frequencies"<<std::endl; return false;}
			else if(dbuff>maxFreq){maxFreq=dbuff;}
			

			freq_[i]=dbuff;
		}
		this->HighestFreq(maxFreq);	
		this->Freq(freq_);
		return true;

	}
	else if (sbuff=="range")
	{
		int nfreq;double inifreq,endfreq;
		freqfile>>inifreq;
		freqfile>>endfreq;
		freqfile>>nfreq;
		//check the validity of the data
		assert (endfreq >inifreq);
		assert (nfreq>0);
		if (inifreq < 0||endfreq<0){std::cout<<"Not allowed negative frequencies"<<std::endl; return false;}

		this ->NumFreq(nfreq);
		this ->HighestFreq(endfreq);
		double dfreq=(endfreq-inifreq)/(nfreq-1);
		std::vector<double> freq_(this->NumFreq());
		for (int i=0;i<this->NumFreq();i++)
		{
			freq_[i]=inifreq+dfreq*i;
		}
		this->Freq(freq_);
		return true;
		
	
	}
	else
	{
		std::cout<<"Unable to read the freq file"<< std::endl;
		return false;
	}

};



