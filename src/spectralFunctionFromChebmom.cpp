		
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

void printHelpMessage();
void printWelcomeMessage();

int main(int argc, char *argv[])
{	
	if (argc != 3)
	{
		printHelpMessage();
		return 0;
	}
	else
		printWelcomeMessage();

	const double broadening  = stod(argv[2]);


	//Read the chebyshev moments from the file
 	chebyshev::Moments1D mu(argv[1]); 
	//and apply the appropiated kernel
	mu.ApplyJacksonKernel(broadening);

	const int num_div = 30*mu.HighestMomentNumber();
	
	const double
	xbound = chebyshev::CUTOFF;
		
	std::vector< double >  energies(num_div,0);
	for( int i=0; i < num_div; i++)
		energies[i] = -xbound + i*(2*xbound)/(num_div-1) ;



	std::cout<<"Computing the Spectral function using "<<mu.HighestMomentNumber()<<" moments "<<std::endl;
	std::cout<<"Over a grid normalized grid  ("<<-xbound<<","<<xbound<<") in steps of "<<energies[1] - energies[0]<< std::endl;
	std::cout<<"The first 10 moments are"<<std::endl;
	for(int m0 =0 ; m0< 1; m0++ )
		std::cout<<mu(m0).real()<<" "<<mu(m0).imag()<<std::endl;
	

	std::string
	outputName  ="mean"+mu.SystemLabel()+"JACKSON.dat";

	std::cout<<"Saving the data in "<<outputName<<std::endl;
	std::cout<<"PARAMETERS: "<< mu.SystemSize()<<" "<<mu.HalfWidth()<<std::endl;
	std::ofstream outputfile( outputName.c_str() );

	std::vector<double> output(num_div,0.0);
	#pragma omp parallel for
	for( int i=0; i < num_div; i++)
	{
		output[i] = 0;
		const double energ = energies[i];
		for( int m0 = 0 ; m0 < mu.HighestMomentNumber() ; m0++)
			output[i] += delta_chebF(energ,m0)*mu(m0).real() ;
		output[i] *= mu.SystemSize()*mu.ScaleFactor();
	}

	for( int i=0; i < num_div; i++)
		outputfile<<energies[i]/mu.ScaleFactor() + mu.ShiftFactor() <<" "<< output[i] <<std::endl;

	outputfile.close();

	std::cout<<"The program finished succesfully."<<std::endl;
return 0;
}
	

void printHelpMessage()
{
	std::cout << "The program should be called with the following options: moments_filename broadening(meV)" << std::endl
			  << std::endl;
	std::cout << "moments_filename will be used to look for .chebmom1D file" << std::endl;
	std::cout << "broadening in (meV) will define the broadening of the delta functions" << std::endl;
};

void printWelcomeMessage()
{
	std::cout << "WELCOME: This program will compute the chebyshev sum of the spectral function" << std::endl;
};
