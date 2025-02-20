#ifndef CHEB_VEC
#define CHEB_VEC

#include <array>
#include <vector>
#include <complex>
#include <string>
#include <cassert>
#include <iostream>		/* for std::cout mostly */
#include <fstream>		/* for std::ofstream and std::ifstream functions classes*/

using namespace std;


template <class T>
class VectorList 
{
	public:
	typedef vector< T > vector_t;
	typedef vector< vector_t > list_t;
	
	//Default constructor. Constructs an empty container with a default-constructed allocator.
	VectorList( ){};
	
	VectorList( size_t list_count, size_t count, const T value = T(0) ):
	_list_size(list_count),
	_vector_size(count),
	_data( list_t(list_count) )
	{
		 for( auto& elem : this->List() )
			elem = vector_t(count); 
	};


	inline
	size_t VectorSize() const { return _vector_size; }
	
	inline
	size_t ListSize() const { return _list_size; }
	
	inline 
	list_t& List(){ return this->_data; };

	inline 
	vector_t& ListElem(const size_t i){ return this->_data[i]; };
	
	inline 
	T& operator()(const size_t i, const size_t j ){return this->_data[i][j];} 
        

	private:
	size_t _list_size, _vector_size ;
	vector< vector< T> > _data;

};


#endif 
