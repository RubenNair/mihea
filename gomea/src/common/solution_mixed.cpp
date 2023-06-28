#include "gomea/src/common/solution_mixed.hpp"

namespace gomea{

solution_mixed::solution_mixed( int number_of_variables_, size_t alphabetSize_, int number_of_c_variables_ ) : solution_t(number_of_variables_, alphabetSize_), c_variables(vec_t<double>(number_of_c_variables_)) {}

solution_mixed::solution_mixed( vec_t<char> &variables, vec_t<double> &c_variables ) : solution_t(variables), c_variables(vec_t<double>(c_variables)) {}

solution_mixed::solution_mixed( size_t numberOfVariables_, size_t alphabetSize_, size_t numberOfCVariables_, fitness_t<char> *problemInstance_ ) : solution_mixed(numberOfVariables_, alphabetSize_, numberOfCVariables_)
{
	// solution_t(numberOfVariables_, alphabetSize_);
    fill(c_variables.begin(), c_variables.end(), 0);
    this->problemInstance = problemInstance_;
}

void solution_mixed::randomInit(std::mt19937 *rng)
{
	for (int i = 0; i < getNumberOfVariables(); ++i)
	{
		variables[i] = (*rng)() % getAlphabetSize();
	}

    for (int i = 0; i < getNumberOfCVariables(); ++i) 
    {
        c_variables[i] = problemInstance->getLowerRangeBound(i) + ((*rng)() / (double)(*rng).max()) * (problemInstance->getUpperRangeBound(i) - problemInstance->getLowerRangeBound(i));
    }
}

int solution_mixed::getNumberOfCVariables() const
{
	return c_variables.size();
}

void solution_mixed::insertCVariables( const vec_t<double> &vars_to_insert )
{
	assert( vars_to_insert.size() == c_variables.size() );
	for( size_t i = 0; i < c_variables.size(); i++ )
	{
		c_variables[i] = vars_to_insert[i];
	}
}

void solution_mixed::insertCVariables( vec_t<double> vars_to_insert, vec_t<int> indices_to_insert )
{
	for( size_t i = 0; i < indices_to_insert.size(); i++ )
	{
		int ind = indices_to_insert[i];
		c_variables[ind] = vars_to_insert[i];
	}
}


void solution_mixed::insertSolution( solution_mixed *solution )
{
	insertVariables(solution->variables);
    insertCVariables(solution->c_variables);
	setObjectiveValues( solution->getObjectiveValues() );
	setConstraintValue( solution->getConstraintValue() );
	setFitnessBuffers( solution->fitness_buffers );

	d_objective_value = solution->d_objective_value;
	c_objective_value = solution->c_objective_value;
}

void solution_mixed::print()
{
	for( size_t i = 0; i < variables.size(); i++ )
		printf("%c ", variables[i]);
    for( size_t i = 0; i < c_variables.size(); i++ )
	{
        printf("%6.5f ",(double)c_variables[i]);
	}
	printf("\n");
}


}
