#pragma once

#include "gomea/src/common/gomea_defs.hpp"
#include "gomea/src/common/partial_solution.hpp"
#include "gomea/src/mixed_integer/Config.hpp"
#include "gomea/src/common/solution.hpp"

namespace gomea{

using gomea::fitness::fitness_t;


class solution_mixed : public solution_t<char>
{

    // laat variables discrete, voeg real_variables toe als doubles
	public:
    vec_t<double> c_variables;
	double d_objective_value;
	double c_objective_value;
    solution_mixed( int number_of_variables_, size_t alphabetSize_, int number_of_c_variables_ );
    solution_mixed (vec_t<char> &variables, vec_t<double> &c_variables);
    solution_mixed( size_t numberOfVariables_, size_t alphabetSize_, size_t numberOfCVariables_, fitness_t<char> *problemInstance_ );
    
    bool operator==(const solution_mixed &solutionB)
    {
        for (int i = 0; i < getNumberOfVariables(); ++i)
        {
            if (this->variables[i] != solutionB.variables[i])
                return false;
        }
        for (int i = 0; i < getNumberOfCVariables(); ++i)
        {
            if (this->c_variables[i] != solutionB.c_variables[i])
                return false;
        }
        return true;
    }

	bool operator!=(const solution_mixed &solutionB)
	{
		return !(*this == solutionB);
	}

	friend std::ostream &operator<<(std::ostream &out, const solution_mixed &solution)
	{
		for (int i = 0; i < solution.getNumberOfVariables(); ++i)
			out << +solution.variables[i];

		for(int i = 0; i < solution.getNumberOfCVariables(); ++i)
			out << +solution.c_variables[i];
		out << " | " << solution.getObjectiveValue();
		return out;
	}

	int getNumberOfCVariables() const;
	void randomInit(std::mt19937 *rng);
	void insertCVariables( const vec_t<double> &vars_to_insert );
	void insertCVariables( vec_t<double> vars_to_insert, vec_t<int> indices_to_insert );
	void insertSolution( solution_mixed *solution );
	void print();

	fitness_t<char> *problemInstance;
// {
// 	for (int i = 0; i < getNumberOfCVariables(); ++i)
// 	{
// 		cVariables[i] = problemInstance->getLowerRangeBound(i) + ((*rng)() / (double)(*rng).max()) * (problemInstance->getUpperRangeBound(i) - problemInstance->getLowerRangeBound(i));
// 	}
// }

	
	// 	solution_t( int number_of_variables );
	// 	solution_t( vec_t<char> &variables );
	// 	solution_t( size_t numberOfVariables_, size_t alphabetSize_ );
	// 	solution_t(size_t numberOfDVariables_, size_t numberOfCVariables_, size_t alphabetSize_, fitness_t<xhar> *problemInstance_);
	// 	// TODO RUBEN: need a solution_t initializer for continuous; will need to pass 

	// 	bool operator==(const solution_t<T> &solutionB)
	// 	{
	// 		for (int i = 0; i < getNumberOfVariables(); ++i)
	// 		{
	// 			if (this->variables[i] != solutionB.variables[i])
	// 				return false;
	// 		}
	// 		return true;
	// 	}

	// 	friend std::ostream &operator<<(std::ostream &out, const solution_t<T> &solution)
	// 	{
	// 		for (int i = 0; i < solution.getNumberOfVariables(); ++i)
	// 			out << +solution.variables[i];
	// 		out << " | " << solution.getObjectiveValue();
	// 		return out;
	// 	}
		
	// 	void randomInit(std::mt19937 *rng);

	// 	void initMemory(int number_of_objectives, int number_of_fitness_buffers);
	// 	void initObjectiveValues(int number_of_objectives);
	// 	void initFitnessBuffers(int number_of_fitness_buffers);

	// 	int getNumberOfVariables() const;
	// 	int getNumberOfDVariables() const;
	// 	int getNumberOfCVariables() const;
	// 	int getNumberOfObjectives() const;
	// 	double getObjectiveValue( int objective_value_index = 0 ) const;
	// 	const vec_t<double> getObjectiveValues() const;
	// 	double getPartialObjectiveValue( int subfunction_index ) const;
	// 	double getConstraintValue() const;
	// 	double getPartialConstraintValue( int subfunction_index ) const;

	// 	void setObjectiveValue( double v );
	// 	void setObjectiveValue( int objective_value_index, double v );
	// 	void setObjectiveValues( const vec_t<double> &v );
	// 	void setConstraintValue( double v );
	// 	void setPartialObjectiveValue( int subfunction_index, double v );
	// 	void setPartialConstraintValue( int subfunction_index, double v );

	// 	double getFitnessBuffer( int buffer_index ) const;
	// 	const vec_t<double> getFitnessBuffers() const;
	// 	void addToFitnessBuffer( int buffer_index, double partial_fitness );
	// 	void subtractFromFitnessBuffer( int buffer_index, double partial_fitness );
	// 	void setFitnessBuffers( const vec_t<double> &buffers );
	// 	void clearFitnessBuffers();

	// 	partial_solution_t<T> getPartialCopy( const vec_t<int> &variable_indices ) const;
	// 	const vec_t<T> getCopyOfVariables( const vec_t<int> &variable_indices = vec_t<int>()) const;
	// 	void insertVariables( const vec_t<T> &vars_to_insert );
	// 	void insertVariables(vec_t<T> vars_to_insert, vec_t<int> indices_to_insert);
	// 	void insertSolution( solution_t<T> *solution );
	// 	void insertPartialSolution( partial_solution_t<T> *solution );

	// 	void print();

	// 	size_t getAlphabetSize();
		
	// 	vec_t<T> variables;
	// 	vec_t<double> fitness_buffers;
	// 	vec_t<T> dVariables;
	// 	bool is_MI_solution = false;
	// 	fitness_t *problemInstance;

	// private:
	// 	vec_t<double> objective_values;
	// 	double constraint_value;
		
	// 	vec_t<double> partial_objective_values;
	// 	vec_t<double> partial_constraint_values;

	// 	size_t alphabetSize = 0;
};

}
