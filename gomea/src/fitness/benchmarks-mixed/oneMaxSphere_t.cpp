/*-=-=-=-=-=-=-=-=-=-=-=-=-=-= Section Includes -=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
#include "gomea/src/fitness/benchmarks-mixed.hpp"
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
namespace gomea{
namespace fitness{

using namespace gomea;

oneMaxSphere_t::oneMaxSphere_t( int number_of_variables, int number_of_c_variables ) : GBOFitnessFunction_t<int>(number_of_variables)
{
	this->name = "F1: OneMax-Sphere function";
    this->number_of_c_variables = number_of_c_variables;
	this->vtr = 10e-10;
	this->use_vtr = true;
	this->optimization_mode = opt_mode::MIN;
	this->initialize();
}

double oneMaxSphere_t::getLowerRangeBound( int dimension )
{
	return( -1e308 );
}
		
double oneMaxSphere_t::getUpperRangeBound( int dimension )
{
	return( 1e308 );
}
		
int oneMaxSphere_t::getNumberOfSubfunctions()
{
	return number_of_variables;
}
		
vec_t<int> oneMaxSphere_t::inputsToSubfunction( int subfunction_index )
{
	vec_t<int> vec;
	vec.push_back(subfunction_index);
	return vec;
}
		
double oneMaxSphere_t::subfunction( int subfunction_index, vec_t<int> &variables )
{
	if(optimization_mode == opt_mode::MAX)
		return( variables[subfunction_index] == 1 ? 1.0 : 0.0 );
	else
		return( variables[subfunction_index] == 1 ? 0.0 : 1.0 );

    
}

double oneMaxSphere_t::discrete_subfunction(int subfunction_index, vec_t<int> &variables)
{
    if(optimization_mode == opt_mode::MAX)
		return( variables[subfunction_index] == 1 ? 1.0 : 0.0 );
	else
		return( variables[subfunction_index] == 1 ? 0.0 : 1.0 );
}

double oneMaxSphere_t::continuous_subfunction(int subfunction_index, vec_t<double> &c_variables)
{
    return c_variables[subfunction_index] * c_variables[subfunction_index];
}

void oneMaxSphere_t::evaluationFunction( solution_t<int> *solution )
{
    // solution = static_cast<solution_mixed*>(solution);
	solution_mixed *solution_mix = static_cast<solution_mixed*>(solution);
	// solution_m->getNumberOfCVariables();
    solution_mix->initFitnessBuffers(getNumberOfFitnessBuffers());
	solution_mix->clearFitnessBuffers();
	for( int i = 0; i < solution_mix->getNumberOfVariables(); i++ )
	{
		int buffer_index = this->getIndexOfFitnessBuffer(i);
		double fsub = discrete_subfunction(i, solution_mix->variables);
		solution_mix->addToFitnessBuffer(buffer_index, fsub);
	}

    for( int i = 0; i < solution_mix->getNumberOfCVariables(); i++ )
    {
        int buffer_index = this->getIndexOfFitnessBuffer(i);
        double fsub = continuous_subfunction(i, solution_mix->c_variables);
        solution_mix->addToFitnessBuffer(buffer_index, fsub);
    }

	for( int i = 0; i < this->number_of_objectives; i++ )
	{
		double ffitness = objectiveFunction( i, solution_mix );
		solution_mix->setObjectiveValue(ffitness);
	}
	double fcons = constraintFunction(solution);
	solution_mix->setConstraintValue(fcons);

	this->full_number_of_evaluations++;
	this->number_of_evaluations++;
}

double oneMaxSphere_t::objectiveFunction( int objective_index, solution_t<int> *solution )
{
    return objectiveFunction(objective_index,solution->fitness_buffers);
}

double oneMaxSphere_t::objectiveFunction( int objective_index, vec_t<double> &fitness_buffers )
{
    return fitness_buffers[objective_index];
}

double oneMaxSphere_t::constraintFunction( solution_t<int> *solution )
{
    return 0.0;
}


}}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
