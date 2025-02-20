
// C & C++ libraries
#include <iostream>		/* for std::cout mostly */
#include <string>		/* for std::string class */
#include <fstream>		/* for std::ofstream and std::ifstream functions classes*/
#include <vector>		/* for std::vector mostly class*/
#include <complex>		/* for std::vector mostly class*/
#include <omp.h>

#include "chebyshev_solver.hpp"
#include "chebyshev_coefficients.hpp"
#include "chebyshev_moments.hpp"
#include "optical_conductivity_tools.hpp"


void printHelpMessage();
void printWelcomeMessage();







int main(int argc, char *argv[])
{	
	if (argc != 4)
	{
		printHelpMessage();
		return 0;
	}
	else
		printWelcomeMessage();

	const double broadening  = stod(argv[3]);
	const int spinless=true;
	std::array<double,2>  matterBounds = AC_utils::MatterBounds();



	//Read the chebyshev moments from the file
 	chebyshev::MomentsAC mu(argv[1]); 
	if (!mu.Freq(argv[2])){return 0;}
	mu.RescaleFreq();
	//and apply the appropiated kernel
	mu.ApplyJacksonKernel(broadening,broadening);

	int prefactorSpin=1;
	if (spinless ){prefactorSpin*=2;}

	const int num_div = 2*mu.HighestMomentNumber()+1;
	
	const double
	xboundm = matterBounds[0]*mu.ScaleFactor()+mu.ShiftFactor();
	const double
	xboundp = matterBounds[1]*mu.ScaleFactor()+mu.ShiftFactor();
		
	std::vector< double >  energies(num_div,0);
	for( int i=0; i < num_div; i++)
		energies[i] = xboundm + i*(xboundp-xboundm)/(num_div-1) ;
	
	



	std::cout<<"Computing the Optical Conductivity kernel using "<<mu.HighestMomentNumber(0)<<" x "<<mu.HighestMomentNumber(1)<<" moments "<<std::endl;
	std::cout<<"Integrated over a grid normalized grid  ("<<xboundm<<","<<xboundp<<") in steps of "<<energies[1] - energies[0]<< std::endl;
	std::cout<<"The first 10 moments are"<<std::endl;
	for(int m0 =0 ; m0< 1; m0++ )
	for(int m1 =0 ; m1< 10; m1++ )
		std::cout<<mu(m0,m1).real()<<" "<<mu(m0,m1).imag()<<std::endl;
	

	std::string
	outputName  ="OpticalConductivity_"+mu.SystemLabel()+"JACKSON.dat";

	std::cout<<"Saving the data in "<<outputName<<std::endl;
	std::cout<<"PARAMETERS: "<< mu.SystemSize()<<" "<<mu.HalfWidth()<<std::endl;
	std::ofstream outputfile( outputName.c_str() );
	
	outputfile<<mu.NumFreq()<<' '<<num_div-1<<std::endl;
	
	std::vector<complex<double>> kernel(num_div,0.0);
 	
	std::cout<< "He llegado"<<std::endl<<GammaAC_chebF(6,30,.3,.1)<<std::endl;


	for (int nomega =0; nomega<mu.NumFreq();nomega++)
	{
		double current_omega=mu.CurrentFreq();
		std::cout<<"Computation for Freq "<<nomega<<std::endl;
		std::cout<<"Freq = "<< current_omega/mu.ScaleFactor()<<std::endl;
		#pragma omp parallel for
		for( int i=0; i < num_div; i++)
		{	
			kernel[i]=0.0;
			const double energ = energies[i];
			for( int m0 = 0 ; m0 < mu.HighestMomentNumber(0) ; m0++)
			for( int m1 = 0 ; m1 < mu.HighestMomentNumber(1) ; m1++)
				kernel[i] +=  GammaAC_chebF(m0,m1,current_omega,energ)*mu(m0,m1)  ;
			kernel[i] *=  mu.SystemSize()*mu.ScaleFactor()*mu.ScaleFactor();
		}
	
		complex <double> acc = 0;
		for( int i=0; i < num_div-1; i++)
		{
			const double energ  = energies[i];
			const double denerg = energies[i+1]-energies[i];
			acc +=kernel[i]*denerg;
			outputfile<<current_omega/mu.ScaleFactor()<<' '<<energ/mu.ScaleFactor() - mu.ShiftFactor() <<" "<<acc.real()<<' '<<acc.imag() <<std::endl;
		}
		mu.IterateFreq();
	}

	outputfile.close();

	std::cout<<"The program finished succesfully."<<std::endl;
return 0;
}
	

void printHelpMessage()
{
	std::cout << "The program should be called with the following options: moments_filename Freq_file broadening(meV)" << std::endl
			  << std::endl;
	std::cout << "moments_filename will be used to look for .chebmom2D file" << std::endl;
	std::cout << "Freq_file is the file that stores the Frequencies. (See documentation for format)" << std::endl;
	std::cout << "broadening in (meV) will define the broadening of the delta functions" << std::endl;
};

void printWelcomeMessage()
{
	std::cout << "WELCOME: This program will compute the chebyshev sum of the kubo-bastin formula for non equlibrium properties" << std::endl;
};
