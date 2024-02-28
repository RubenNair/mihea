/*-=-=-=-=-=-=-=-=-=-=-=-=-=-= Section Includes -=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
#include "gomea/src/fitness/benchmarks-discrete.hpp"
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
namespace gomea{
namespace fitness{

using namespace gomea;

oneMax_t::oneMax_t( int number_of_variables ) : GBOFitnessFunction_t<char>(number_of_variables)
{
	this->name = "OneMax function";
	this->vtr = number_of_variables;
	this->use_vtr = true;
	this->optimization_mode = opt_mode::MAX;
	this->initialize();
}
		
int oneMax_t::getNumberOfSubfunctions()
{
	return number_of_variables;
}
		
vec_t<int> oneMax_t::inputsToSubfunction( int subfunction_index )
{
	vec_t<int> vec;
	vec.push_back(subfunction_index);
	return vec;
}
		
double oneMax_t::subfunction( int subfunction_index, vec_t<char> &variables )
{
	if(optimization_mode == opt_mode::MAX)
		return( variables[subfunction_index] == 1 ? 1.0 : 0.0 );
	else
		return( variables[subfunction_index] == 1 ? 0.0 : 1.0 );
}

}}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
