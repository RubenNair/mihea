/*-=-=-=-=-=-=-=-=-=-=-=-=-=-= Section Includes -=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
#include "gomea/src/fitness/benchmarks-mixed.hpp"
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
namespace gomea{
namespace fitness{

using namespace gomea;

DT5Sphere_t::DT5Sphere_t( int number_of_variables, int number_of_c_variables ) : GBOFitnessFunction_t<char>(number_of_variables)
{
	this->name = "F3: DeceptiveTrap5-Sphere function";
    this->number_of_c_variables = number_of_c_variables;
	this->vtr = 10e-10;
	this->use_vtr = true;
	this->optimization_mode = opt_mode::MIN;
	this->initialize();
}

double DT5Sphere_t::getLowerRangeBound( int dimension )
{
	return( -10 );
}
		
double DT5Sphere_t::getUpperRangeBound( int dimension )
{
	return( 10 );
}
		
int DT5Sphere_t::getNumberOfSubfunctions()
{
	return number_of_variables;
}
		
vec_t<int> DT5Sphere_t::inputsToSubfunction( int subfunction_index )
{
	vec_t<int> vec;
	vec.push_back(subfunction_index);
	return vec;
}
		
double DT5Sphere_t::subfunction( int subfunction_index, vec_t<char> &variables )
{
	if(optimization_mode == opt_mode::MAX)
		return( variables[subfunction_index] == '\001' ? 1.0 : 0.0 );
	else
		return( variables[subfunction_index] == '\001' ? 0.0 : 1.0 );

    
}

double DT5Sphere_t::discrete_subfunction(int subfunction_index, vec_t<char> &variables)
{
    int sum = 0;
	for(int i = 0; i < 5; i++) {
		sum += variables[(subfunction_index*5) + i] == '\001' ? 1 : 0;
	}

    if(optimization_mode == opt_mode::MAX)
		return( sum == 5 ? 1.0 : (4 - sum) / 5.0 );
	else
		return( sum == 5 ? 0.0 : 1 - ((4 - sum) / 5.0) );
}

double DT5Sphere_t::continuous_subfunction(int subfunction_index, vec_t<double> &c_variables)
{
    return c_variables[subfunction_index] * c_variables[subfunction_index];
}

void DT5Sphere_t::evaluationFunction( solution_t<char> *solution )
{
    // solution = static_cast<solution_mixed*>(solution);
	solution_mixed *solution_mix = static_cast<solution_mixed*>(solution);
	// solution_m->getNumberOfCVariables();
    solution_mix->initFitnessBuffers(getNumberOfFitnessBuffers());
	solution_mix->clearFitnessBuffers();
	
    solution_mix->d_objective_value = 0.0;
    solution_mix->c_objective_value = 0.0;

	for( int i = 0; i < solution_mix->getNumberOfVariables() / 5; i++ )
	{
		int buffer_index = this->getIndexOfFitnessBuffer(i);
		double fsub = discrete_subfunction(i, solution_mix->variables);
		solution_mix->addToFitnessBuffer(buffer_index, fsub);
		solution_mix->d_objective_value += fsub;
	}

    for( int i = 0; i < solution_mix->getNumberOfCVariables(); i++ )
    {
        int buffer_index = this->getIndexOfFitnessBuffer(i);
        double fsub = continuous_subfunction(i, solution_mix->c_variables);
        solution_mix->addToFitnessBuffer(buffer_index, fsub);
		solution_mix->c_objective_value += fsub;
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

// template<class T>
// void DT5Sphere_t::partialEvaluationFunction( solution_t<T> *parent, partial_solution_t<T> *solution)
// {
//     evaluatePartialSolution(parent, solution);
// }

// template<class T>
// void DT5Sphere_t::evaluatePartialSolution( solution_t<T> *parent, partial_solution_t<T> *solution)
// {
//     parent = static_cast<solution_mixed*>(parent);
//     solution->initFitnessBuffers(getNumberOfFitnessBuffers());
// 	solution->resetFitnessBuffers();

// 	std::set<int> touched_subfunctions;
// 	//if( dependent_subfunctions.size() == 0 )
// 	{
// 		for( int ind : solution->touched_indices )
// 		{
// 			assert( this->subfunction_dependency_map[ind].size() > 0 );
// 			touched_subfunctions.insert(this->subfunction_dependency_map[ind].begin(), this->subfunction_dependency_map[ind].end());
// 		}
// 	}
// 	//else
// 		//touched_subfunctions = dependent_subfunctions;
	
// 	double objective_value_delta = 0.0;
// 	// Calculate sum of touched subfunctions for parent
// 	for( int subfunction_index : touched_subfunctions )
// 	{
// 		double subf_result = subfunction( subfunction_index, parent->variables );
// 		objective_value_delta -= subf_result; 
// 		int buffer_index = this->getIndexOfFitnessBuffer(subfunction_index);
// 		solution->subtractFromFitnessBuffer( buffer_index, subf_result );
// 	}
	
// 	// Create backup of parent variables before modification
// 	vec_t<T> partial_backup = parent->getCopyOfVariables( solution->touched_indices );

// 	// Insert variables of partial solution and then calculate sum of touched subfunctions for offspring
// 	parent->insertVariables( solution->touched_variables, solution->touched_indices );
// 	for( int subfunction_index : touched_subfunctions ) 
// 	{
// 		double subf_result = subfunction( subfunction_index, parent->variables );
// 		//solution->partial_objective_values[subfunction_index] = subf_result;
// 		objective_value_delta += subf_result;
// 		int buffer_index = this->getIndexOfFitnessBuffer(subfunction_index);
// 		solution->addToFitnessBuffer( buffer_index, subf_result );
// 	}

// 	// Return parent variables to original state
// 	parent->insertVariables(partial_backup, solution->touched_indices);

// 	// Update fitness of partial solution
// 	//solution->setObjectiveValue(parent->getObjectiveValue() + objective_value_delta);
// 	//solution->setConstraintValue(parent->getConstraintValue());

// 	// Add parent buffer for final result of buffer	
// 	vec_t<double> parent_buffers = parent->fitness_buffers;
// 	for( size_t i = 0; i < parent_buffers.size(); i++ )
// 		solution->addToFitnessBuffer(i,parent_buffers[i]);

// 	// Apply function to calculate objective value from fitness buffers
// 	for( int i = 0; i < this->number_of_objectives; i++ )
// 	{
// 		double ffitness = objectiveFunction( i, solution );
// 		solution->setObjectiveValue(ffitness);
// 	}
// 	double fcons = constraintFunction(solution);
// 	solution->setConstraintValue(fcons);

// 	this->full_number_of_evaluations++;
// 	this->number_of_evaluations += touched_subfunctions.size() / (double) this->getNumberOfSubfunctions();
// }    

double DT5Sphere_t::objectiveFunction( int objective_index, solution_t<char> *solution )
{
    return objectiveFunction(objective_index,solution->fitness_buffers);
}

double DT5Sphere_t::objectiveFunction( int objective_index, vec_t<double> &fitness_buffers )
{
    return fitness_buffers[objective_index];
}

double DT5Sphere_t::constraintFunction( solution_t<char> *solution )
{
    return 0.0;
}

// template<class T>
// void DT5Sphere_t::evaluatePartialSolutionBlackBox(solution_t<T> *parent, partial_solution_t<T> *solution)
// {
//     evaluatePartialSolution(parent, solution);
// }



}}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
