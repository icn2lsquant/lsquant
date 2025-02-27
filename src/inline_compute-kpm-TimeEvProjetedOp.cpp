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


namespace time_evolution
{
	void printHelpMessage();
	void printWelcomeMessage();
}


int main(int argc, char *argv[])
{
	if ( !(argc == 7 || argc == 8 ) )
	{
		time_evolution::printHelpMessage();
		return 0;
	}
	else
		time_evolution::printWelcomeMessage();
	
	const std::string
		LABEL  = argv[1],
		S_OP   = argv[2],
		S_OPPRJ= argv[3],
		S_NMOM = argv[4],
		S_NTIME= argv[5],
		S_TMAX = argv[6];

	const int numMoms = atoi(S_NMOM.c_str() );
	const int numTimes= atoi(S_NTIME.c_str() );
	const double tmax = stod(S_TMAX );

	chebyshev::MomentsTD chebMoms(numMoms, numTimes); //load number of moments

	
	SparseMatrixType OP[3];
	OP[0].SetID("HAM");
	OP[1].SetID(S_OP);
	OP[2].SetID(S_OPPRJ);
	
	// Build the operators from Files
	SparseMatrixBuilder builder;
	std::array<double,2> spectral_bounds;	
	for (int i = 0; i < 3; i++)
	{
		std::string input = "operators/" + LABEL + "." + OP[i].ID() + ".CSR";
		builder.setSparseMatrix(&OP[i]);
		builder.BuildOPFromCSRFile(input);
	
		if( i == 0 ) //is hamiltonian
		//Obtain automatically the energy bounds
		 spectral_bounds = chebyshev::utility::SpectralBounds(OP[0]);
	};
	//CONFIGURE THE CHEBYSHEV MOMENTS
	chebMoms.SystemLabel(LABEL);
	chebMoms.BandWidth ( (spectral_bounds[1]-spectral_bounds[0])*1.0);
	chebMoms.BandCenter( (spectral_bounds[1]+spectral_bounds[0])*0.5);
	chebMoms.TimeDiff( tmax/(numTimes-1) );
	chebMoms.SetAndRescaleHamiltonian(OP[0]);
	chebMoms.Print();

	//Compute the chebyshev expansion table
	qstates::generator gen;
	if( argc == 8)
		gen  = qstates::LoadStateFile(argv[7]);

	chebyshev::TimeEvolvedProjectedOperator( OP[1], OP[2], chebMoms, gen);

	
	std::string outputfilename="TimeProj"+S_OP+"-"+S_OPPRJ+LABEL+"KPM_M"+S_NMOM+"_state"+gen.StateLabel()+".chebmomTD";	
	std::cout<<"Saving convergence data in "<<outputfilename<<std::endl;
	chebMoms.saveIn(outputfilename);
	std::cout<<"End of program"<<std::endl;
	return 0;
}

void time_evolution::printHelpMessage()
	{
		std::cout << "The program should be called with the following options: Label Op PrOp numMom numTimeSteps MaxTime (optional) num_states_file" << std::endl
				  << std::endl;
		std::cout << "Label will be used to look for Label.Ham, Label.Op and Label.PrOp" << std::endl;
		std::cout << "Op and PrOp will be used to located the sparse matrix file of two operators for the projection evolution" << std::endl;
		std::cout << "numMom will be used to set the number of moments in the chebyshev table" << std::endl;
		std::cout << "numTimeSteps  will be used to set the number of timesteps in the chebyshev table" << std::endl;
		std::cout << "TimeMax  will be set the maximum time where the correlation will be evaluted " << std::endl;
	};

	inline
void time_evolution::printWelcomeMessage()
	{
		std::cout << "WELCOME: This program will compute a table needed for expanding the correlation function in Chebyshev polynomialms" << std::endl;
	};
