#pragma once

#include "gomea/src/common/gomea_defs.hpp"
#include "gomea/src/common/partial_solution.hpp"
#include "gomea/src/mixed_integer/Config.hpp"
#include "gomea/src/common/solution.hpp"

namespace gomea{

using gomea::fitness::fitness_t;


class solution_mixed : public solution_t<int>
{

	public:
    vec_t<double> c_variables;
    solution_mixed( int number_of_variables_, size_t alphabetSize_, int number_of_c_variables_ );
    solution_mixed (vec_t<int> &variables, vec_t<double> &c_variables);
    solution_mixed( size_t numberOfVariables_, size_t alphabetSize_, size_t numberOfCVariables_, fitness_t<int> *problemInstance_ );
    solution_mixed( const solution_mixed &other );
	virtual ~solution_mixed() = default;

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
	virtual void insertSolution( solution_mixed *solution );
	void print();

	virtual solution_mixed *clone();

	fitness_t<int> *problemInstance;

	
	clock_t getTimeStamp() const;
	void setTimeStamp(clock_t timeStamp);

	size_t getNumberOfFullEvaluations() const;
	void setNumberOfFullEvaluations(size_t numberOfFullEvaluations);


	clock_t timeStamp;                      // The time stamp when a solution was evaluated
	size_t numberOfFullEvaluations;         // The number of full evaluations executed to get the solution
};

}
