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


namespace LDOS
{
	void printHelpMessage();
	void printWelcomeMessage();
}


int main(int argc, char *argv[])
{
	bool saveInBin=true;


	if ( !(argc == 5 || argc == 6 ) )
	{
		LDOS::printHelpMessage();
		return 0;
	}
	else
		LDOS::printWelcomeMessage();
	
	const std::string
		LABEL   = argv[1],
		S_OP    = argv[2],
		S_NMOM  = argv[3],
		S_NSTAT = argv[4];

	const int numMoms   = atoi(S_NMOM.c_str()  );
	const int numStates = atoi(S_NSTAT.c_str() );







	
	SparseMatrixType OP[2];
	OP[0].SetID("HAM");
	OP[1].SetID(S_OP);
	
	// Build the operators from Files
	SparseMatrixBuilder builder;
	std::array<double,2> spectral_bounds;	
	for (int i = 0; i < 2; i++)
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
	
	std::cout<<"Building the Moments table"<<std::endl;
	std::cout<<"Dimensions are going to be ("<<numMoms<<"x"<<OP[0].rank()<<")"<<std::endl;
	chebyshev::MomentsLocal chebMoms(numMoms, OP[0].rank()); //load number of moments
	std::cout<<"Moments table succesfully built"<<std::endl;

	//CONFIGURE THE CHEBYSHEV MOMENTS
	chebMoms.SystemLabel(LABEL);
	chebMoms.BandWidth ( (spectral_bounds[1]-spectral_bounds[0])*1.0);
	chebMoms.BandCenter( (spectral_bounds[1]+spectral_bounds[0])*0.5);
	chebMoms.SetAndRescaleHamiltonian(OP[0]);
	chebMoms.Print();

	//Compute the chebyshev expansion table
	qstates::generator gen;
	if( argc == 6)
		gen  = qstates::LoadStateFile(argv[5]);
	else gen.NumberOfStates(numStates);

	chebyshev::ComputeDeltaPhi(OP[1],chebMoms, gen );

	if (saveInBin){
	std::string outputfilename="LOCALSpectral"+S_OP+LABEL+"KPM_M"+S_NMOM+"_state_"+S_NSTAT+"_"+gen.StateLabel()+".chebmomLocal.bin";	
	std::cout<<"Saving convergence data in "<<outputfilename<<std::endl;
	chebMoms.saveInBin(outputfilename);
	}else{
	std::string outputfilename="LOCALSpectral"+S_OP+LABEL+"KPM_M"+S_NMOM+"_state"+gen.StateLabel()+".chebmomLocal";	
	std::cout<<"Saving convergence data in "<<outputfilename<<std::endl;
	chebMoms.saveIn(outputfilename);
	}
	std::cout<<"End of program"<<std::endl;
	return 0;
}

void LDOS::printHelpMessage()
	{
		std::cout << "The program should be called with the following options: Label Op numMom numStates (optional) num_states_file" << std::endl
				  << std::endl;
		std::cout << "Label will be used to look for Label.Ham, Label.Op " << std::endl;
		std::cout << "Op will be used to located the sparse matrix file of two operators for the projection evolution" << std::endl;
		std::cout << "numMom will be used to set the number of moments in the chebyshev table" << std::endl;
		std::cout << "numStates will be used to set the number of random phases for the calculation" << std::endl;
	};

	inline
void LDOS::printWelcomeMessage()
	{
		std::cout << "WELCOME: This program will compute a table needed for expanding the Local Density of States in Chebyshev polynomialms" << std::endl;
	};
