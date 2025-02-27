
namespace AC_utils
{
std::array<double,2> MatterBounds()
{
	double highest=0, lowest=0;
	std::ifstream bounds_file( "AC_MatterBounds" );
	if( bounds_file.is_open() )
	{
		bounds_file>>lowest>>highest;
		bounds_file.close();
	}
	else
	{
		lowest = -100;
		highest = 100;
	}
	return { lowest, highest};
};







};





