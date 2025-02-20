#include "sparse_matrix.hpp"

bool Sparse::OPERATOR_FromCSRFile(const std::string input, int &dim, vector<int> &columns, vector<int> &rowIndex, vector<complex<double> > &values)
{
  std::cout << "\nReading the CSR file located at: " << input << std::endl;
  //OPEN MATRIX FILE
  std::ifstream matrix_file(input.c_str());
  if (!matrix_file.is_open())
  {
    std::cerr << "ERROR: Can not open the matrix file";
    return false;
  };

  //READ DIMENSION OF THE MATRIX
  int nnz;
  matrix_file >> dim >> nnz;

  //CREATE ARRAYS TO STORE THE MATRIX
  values   = vector<complex<double> >(nnz);
  columns  = vector<int>(nnz); 
  rowIndex = vector<int>(dim + 1);

  //READ VALUES
  double rev, imv;
  for (int i = 0; i < nnz; i++)
  {
    matrix_file >> rev >> imv;
    values[i] = complex<double>(rev, imv);
  }

  //READ COLUMNS
  int col;
  for (int i = 0; i < nnz; i++)
  {
    matrix_file >> col;
    columns[i] = col;
  }

  //READ ROW_INDEX_ARRAY
  int rowIdx;
  for (int i = 0; i < dim + 1; i++)
  {
    matrix_file >> rowIdx;
    rowIndex[i] = rowIdx;
  }
  matrix_file.close();

  std::cout << "FINISH READING.Status: SUCCED" << std::endl;
  return true;
};
