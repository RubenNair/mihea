// WIP: This is a copy of the DT5Sphere_t class. It is not yet fully adapted to the BN case.

/*-=-=-=-=-=-=-=-=-=-=-=-=-=-= Section Includes -=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
#include "gomea/src/fitness/benchmarks-mixed.hpp"
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
namespace gomea{
namespace fitness{

using namespace gomea;

BNStructureLearning::BNStructureLearning( int number_of_variables, int number_of_c_variables, int problem_index, const shared_ptr<DataStructure<double>> &data, size_t max_number_of_parents, size_t max_number_of_discretizations, bool transformCVariables) : GBOFitnessFunction_t<int>(number_of_variables)
{
	this->name = "Bayesian Network structure learning using density to calculate fitness";
    this->number_of_c_variables = number_of_c_variables;
	this->vtr = 1e308;
	this->use_vtr = false;
	this->optimization_mode = opt_mode::MIN; // Since minimization is hardcoded in some spots, easiest way right now is to negate fitness and still minimize. (In original code, maximization)
	this->transformCVariables = transformCVariables;
	this->initialize();
	// TODO: initialize density here -> will need more / different parameters to this constructor to do so.
	this->density = new Density(problem_index, name, data, max_number_of_parents, max_number_of_discretizations);
}

double BNStructureLearning::getLowerRangeBound( int dimension )
{
	return( this->transformCVariables ? 1.0 : 0.0 );
}
		
double BNStructureLearning::getUpperRangeBound( int dimension )
{
	return( this->transformCVariables ? INFINITY : 0.99 );
}
		
int BNStructureLearning::getNumberOfSubfunctions()
{
	return number_of_variables;
}
		
vec_t<int> BNStructureLearning::inputsToSubfunction( int subfunction_index )
{
	vec_t<int> vec;
	vec.push_back(subfunction_index);
	return vec;
}
		
double BNStructureLearning::subfunction( int subfunction_index, vec_t<int> &variables )
{
	cout << "ERROR: subfunction() should not be called for BN structure learning." << endl;
	return -1;    
}

double BNStructureLearning::discrete_subfunction(int subfunction_index, vec_t<int> &variables)
{
    cout << "ERROR: discrete_subfunction() should not be called for BN structure learning." << endl;
	return -1;
}

double BNStructureLearning::continuous_subfunction(int subfunction_index, vec_t<double> &c_variables)
{
    cout << "ERROR: continuous_subfunction() should not be called for BN structure learning." << endl;
	return -1;
}

void BNStructureLearning::evaluationFunction( solution_t<int> *solution )
{
    // TODO: calculate fitness using density class. Probably need to cast the solution to solution_so type.
	solution_BN *casted_solution = dynamic_cast<solution_BN *>(solution);
	// Make sure boundaries are updated before evaluation (TODO not sure if this should be here, but it works for now)
	casted_solution->updateBoundaries();
	density->computeFitness(*casted_solution);

	// Copy calculated fitness to fitness buffer so the rest of the code can handle it. 
	casted_solution->initFitnessBuffers(getNumberOfFitnessBuffers());
	casted_solution->clearFitnessBuffers();
	// TODO for now assumed only single objective, but maybe weird way to write it like this.
	for( int i = 0; i < this->number_of_objectives; i++ )
	{
		casted_solution->addToFitnessBuffer(i, casted_solution->getFitness());
		double ffitness = objectiveFunction( i, casted_solution );
		casted_solution->setObjectiveValue(ffitness);
	}
	double fcons = constraintFunction(solution);
	casted_solution->setConstraintValue(fcons);
	// density->computeFitness(*solution);
}

const Density *BNStructureLearning::getDensity() const
{
	return this->density;
}

// template<class T>
// void BNStructureLearning::partialEvaluationFunction( solution_t<T> *parent, partial_solution_t<T> *solution)
// {
//     evaluatePartialSolution(parent, solution);
// }

// template<class T>
// void BNStructureLearning::evaluatePartialSolution( solution_t<T> *parent, partial_solution_t<T> *solution)
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

double BNStructureLearning::objectiveFunction( int objective_index, solution_t<int> *solution )
{
    return objectiveFunction(objective_index,solution->fitness_buffers);
}

double BNStructureLearning::objectiveFunction( int objective_index, vec_t<double> &fitness_buffers )
{
    return fitness_buffers[objective_index];
}

double BNStructureLearning::constraintFunction( solution_t<int> *solution )
{
    return 0.0;
}

// template<class T>
// void BNStructureLearning::evaluatePartialSolutionBlackBox(solution_t<T> *parent, partial_solution_t<T> *solution)
// {
//     evaluatePartialSolution(parent, solution);
// }



}}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
