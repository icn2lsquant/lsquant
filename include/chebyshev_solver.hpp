#ifndef CHEBYSHEV_SOLVER
#define CHEBYSHEV_SOLVER

// C & C++ libraries
#include <cassert>   //for assert
#include <array>

#include <vector>    //for std::vector mostly class
#include <numeric>   //for std::accumulate *
#include <algorithm> //for std::max_elem
#include <complex>   ///for std::complex
#include <fstream>   //For ofstream
#include <limits>    //For getting machine precision limit
#include "sparse_matrix.hpp"
#include "chebyshev_moments.hpp"
#include "chebyshev_coefficients.hpp"
#include "linear_algebra.hpp"
#include <omp.h>
#include <chrono>
#include "quantum_states.hpp"
#include "kpm_noneqop.hpp" //Get Batch function
#include "special_functions.hpp"

namespace chebyshev
{
	typedef std::complex<double> value_t;	
	typedef std::vector<value_t> vector_t ;


	namespace utility
	{
		std::array<double,2> SpectralBounds( SparseMatrixType& HAM);
	};



	int ComputeMomTable( chebyshev::Vectors &chevVecL, chebyshev::Vectors& chevVecR , chebyshev::Moments2D& sub);

	int CorrelationExpansionMoments( 	const vector_t& PhiR, const vector_t& PhiL,
										SparseMatrixType &OPL,
										SparseMatrixType &OPR,  
										chebyshev::Vectors &chevVecL,
										chebyshev::Vectors &chevVecR,
										chebyshev::Moments2D &chebMoms
										);

	
	int CorrelationExpansionMoments( SparseMatrixType &OPL, SparseMatrixType &OPR,  chebyshev::Moments2D &chebMoms, qstates::generator& gen );

	int TimeDependentCorrelations( SparseMatrixType &OPL, SparseMatrixType &OPR,  chebyshev::MomentsTD &chebMoms, qstates::generator& gen);

	int SpectralMoments(SparseMatrixType &OP,  chebyshev::Moments1D &chebMoms, qstates::generator& gen);

        int TimeEvolvedProjectedOperator(SparseMatrixType &OP, SparseMatrixType &OPPRJ,  chebyshev::MomentsTD &chebMoms, qstates::generator& gen  );

        int TimeEvolvedOperator(SparseMatrixType &OP,  chebyshev::MomentsTD &chebMoms, qstates::generator& gen  );

	int MeanSquareDisplacement(chebyshev::MomentsTD &chebMoms, qstates::generator& gen  );

	int TimeEvolvedProjectedOperatorWF(SparseMatrixType &OP, SparseMatrixType &OPRJ,  chebyshev::MomentsTD &chebMoms, qstates::generator& gen  );
    
      	int TimeEvolvedOperatorWithWF(chebyshev::MomentsTD &chebMoms, qstates::generator& gen  );

	int ComputeDeltaPhi( SparseMatrixType &OP,  chebyshev::MomentsLocal &chebMoms, qstates::generator& gen);

	int ComputeDeltaPhi2( SparseMatrixType &OP,  chebyshev::MomentsTD &chebMoms, qstates::generator& gen ,double Energy);

}; // namespace chebyshev

#endif
