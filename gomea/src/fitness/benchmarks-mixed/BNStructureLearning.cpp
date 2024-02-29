
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
	this->density = new Density(problem_index, name, data, max_number_of_parents, max_number_of_discretizations);
}

BNStructureLearning::~BNStructureLearning()
{
	delete this->density;
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
	solution_BN *casted_solution = dynamic_cast<solution_BN *>(solution);
	// Make sure boundaries are updated before evaluation
	casted_solution->updateBoundaries();
	density->computeFitness(*casted_solution);

	// Copy calculated fitness to fitness buffer so the rest of the code can handle it. 
	casted_solution->initFitnessBuffers(getNumberOfFitnessBuffers());
	casted_solution->clearFitnessBuffers();
	// NOTE: for now assumed only single objective
	for( int i = 0; i < this->number_of_objectives; i++ )
	{
		casted_solution->addToFitnessBuffer(i, casted_solution->getFitness());
		double ffitness = objectiveFunction( i, casted_solution );
		casted_solution->setObjectiveValue(ffitness);
	}
	double fcons = constraintFunction(solution);
	casted_solution->setConstraintValue(fcons);

	// Update timestamp and number of evaluations of evaluated solution (number of evaluations up to this solution, not including this one)
	clock_t currentTime = clock();
	casted_solution->setTimeStamp(currentTime);

	casted_solution->setNumberOfEvaluations(this->full_number_of_evaluations);
}

const Density *BNStructureLearning::getDensity() const
{
	return this->density;
}

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


}}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
