// C & C++ libraries
#include <iostream> /* for std::cout mostly */
#include <string>   /* for std::string class */
#include <fstream>  /* for std::ofstream and std::ifstream functions classes*/
#include <stdlib.h>
#include <chrono>


#include "kpm_noneqop.hpp" //Message functions
#include "chebyshev_moments.hpp"
#include "sparse_matrix.hpp"
#include "quantum_states.hpp"
#include "chebyshev_solver.hpp"

namespace spectral
{
	void printHelpMessage();
	void printWelcomeMessage();
};


int main(int argc, char *argv[])
{
	if ( !(argc == 4 || argc == 5 ) )
	{
		spectral::printHelpMessage();
		return 0;
	}
	else
		spectral::printWelcomeMessage();
	
	const std::string
		LABEL  = argv[1],
		S_OP   = argv[2],
		S_NMOM = argv[3];

	const int numMoms   = atoi(S_NMOM.c_str() );
	chebyshev::Moments1D chebMoms( numMoms ); //load number of moments

	const int NOPS = 2;
	SparseMatrixType OP[NOPS];
	OP[0].SetID("HAM");
	OP[1].SetID( S_OP );
	
	// Build the operators from Files
	SparseMatrixBuilder builder;
	std::array<double,2> spectral_bounds;	
	
	for (int i = 0; i < NOPS; i++)
	if( OP[i].isIdentity() )
			std::cout<<"The operator "<<OP[i].ID()<<" is treated as the identity"<<std::endl;
	else
	{
		std::string input = "operators/" + LABEL + "." + OP[i].ID() + ".CSR";
		builder.setSparseMatrix(&OP[i]);
		builder.BuildOPFromCSRFile(input);
		
		if( i == 0 ) //is hamiltonian, obtain automatically the energy bounds
		 spectral_bounds = chebyshev::utility::SpectralBounds(OP[0]);
	}
	
	//CONFIGURE THE CHEBYSHEV MOMENTS
	chebMoms.SystemLabel(LABEL);
	chebMoms.BandWidth ( (spectral_bounds[1]-spectral_bounds[0])*1.0);
	chebMoms.BandCenter( (spectral_bounds[1]+spectral_bounds[0])*0.5);
	chebMoms.SetAndRescaleHamiltonian(OP[0]);
	chebMoms.Print();


	//Compute the chebyshev expansion table
	qstates::generator gen;
	if( argc == 5)
		gen  = qstates::LoadStateFile(argv[4]);

		
	//gen.NumberOfStates(std::stoi(argv[4]));




	chebyshev::SpectralMoments(OP[1],chebMoms, gen);

	auto prefix="SpectralOp"+OP[1].ID();
	if( OP[1].isIdentity() )
		prefix="SpectralOp"+OP[1].ID();
	std::string outputfilename=prefix+LABEL+"KPM_M"+S_NMOM+"_state"+gen.StateLabel()+".chebmom1D";	

	std::cout<<"Saving the moments in  "<<outputfilename<<std::endl;
	chebMoms.saveIn(outputfilename);
	std::cout<<"End of program"<<std::endl;
	return 0;
};


	void spectral::printHelpMessage()
	{
		std::cout << "The program should be called with the following options: Label Op numMom " << std::endl
				  << std::endl;
		std::cout << "Label will be used to look for Label.Ham, Label.Op" << std::endl;
		std::cout << "Op will be used to located the  matrix file of the operators for the spectral " << std::endl;
		std::cout << "numMom will be used to set the number of moments in the chebyshev table" << std::endl;
	};

	
	void spectral::printWelcomeMessage()
	{
		std::cout << "WELCOME: This program will compute a list needed for expanding the spectral function in Chebyshev polynomials" << std::endl;
	};



