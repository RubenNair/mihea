#pragma once

#include "gomea/src/fitness/fitness.hpp"
#include "gomea/src/fitness/gbo_fitness.hpp"
#include "gomea/src/fitness/bbo_fitness.hpp"
// #include "gomea/src/common/solution.hpp"
#include "gomea/src/common/solution_mixed.hpp"
#include "gomea/src/common/partial_solution.hpp"
#include "gomea/src/common/gomea_defs.hpp"
#include "gomea/src/utils/tools.hpp"
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

namespace gomea{
namespace fitness{

class oneMaxSphere_t: public GBOFitnessFunction_t<char>
{
	public:
		oneMaxSphere_t( int number_of_variables, int number_of_c_variables );
		int getNumberOfSubfunctions(); 
		vec_t<int> inputsToSubfunction( int subfunction_index );
		int number_of_c_variables;
		double getLowerRangeBound( int dimension );
		double getUpperRangeBound( int dimension );
		
	private:
		double subfunction( int subfunction_index, vec_t<char> &variables );
		double discrete_subfunction(int subfunction_index, vec_t<char> &variables);
		double continuous_subfunction(int subfunction_index, vec_t<double> &c_variables);
		// void evaluationFunction( solution_mixed *solution );
		// double objectiveFunction( int objective_index, solution_mixed *solution );
		// double objectiveFunction( int objective_index, vec_t<double> &fitness_buffers );
		// double constraintFunction( solution_mixed *solution );
		// from gbo_fitness.hpp
		void evaluationFunction( solution_t<char> *solution );
		// void partialEvaluationFunction( solution_t<char> *parent, partial_solution_t<char> *solution );
		double objectiveFunction( int objective_index, solution_t<char> *solution );
		double objectiveFunction( int objective_index, vec_t<double> &fitness_buffers );
		double constraintFunction( solution_t<char> *solution );

		// from fitness.hpp;
		// void evaluate( solution_t<T> *solution ); // Actually not necessary to override this function, I think
		// void evaluatePartialSolution( solution_t<char> *parent, partial_solution_t<char> *solution );
		// bool betterFitness( solution_t<T> *sol_x, solution_t<T> *sol_y ); // Also not necessary to override this function, I think
		// void evaluatePartialSolutionBlackBox( solution_t<char> *parent, partial_solution_t<char> *solution );

};

class oneMaxRotEllip_t: public GBOFitnessFunction_t<char>
{
	public:
		oneMaxRotEllip_t( int number_of_variables, int number_of_c_variables );
		int getNumberOfSubfunctions(); 
		vec_t<int> inputsToSubfunction( int subfunction_index );
		int number_of_c_variables;
		double getLowerRangeBound( int dimension );
		double getUpperRangeBound( int dimension );
		
	private:
		double subfunction( int subfunction_index, vec_t<char> &variables );
		double discrete_subfunction(int subfunction_index, vec_t<char> &variables);
		double continuous_subfunction(int subfunction_index, vec_t<double> &c_variables);
		// void evaluationFunction( solution_mixed *solution );
		// double objectiveFunction( int objective_index, solution_mixed *solution );
		// double objectiveFunction( int objective_index, vec_t<double> &fitness_buffers );
		// double constraintFunction( solution_mixed *solution );
		// from gbo_fitness.hpp
		void evaluationFunction( solution_t<char> *solution );
		// void partialEvaluationFunction( solution_t<char> *parent, partial_solution_t<char> *solution );
		double objectiveFunction( int objective_index, solution_t<char> *solution );
		double objectiveFunction( int objective_index, vec_t<double> &fitness_buffers );
		double constraintFunction( solution_t<char> *solution );
		vec_t<double> rotateParameters( vec_t<double> &parameters);
		// from fitness.hpp;
		// void evaluate( solution_t<T> *solution ); // Actually not necessary to override this function, I think
		// void evaluatePartialSolution( solution_t<char> *parent, partial_solution_t<char> *solution );
		// bool betterFitness( solution_t<T> *sol_x, solution_t<T> *sol_y ); // Also not necessary to override this function, I think
		// void evaluatePartialSolutionBlackBox( solution_t<char> *parent, partial_solution_t<char> *solution );


		double **initializeObjectiveRotationMatrix( double rotation_angle, int rotation_block_size );
		double *rotateVariables( double *variables, int num_variables, double **rotation_matrix );
		
};

class DT5Sphere_t: public GBOFitnessFunction_t<char>
{
	public:
		DT5Sphere_t( int number_of_variables, int number_of_c_variables );
		int getNumberOfSubfunctions(); 
		vec_t<int> inputsToSubfunction( int subfunction_index );
		int number_of_c_variables;
		double getLowerRangeBound( int dimension );
		double getUpperRangeBound( int dimension );
		
	private:
		double subfunction( int subfunction_index, vec_t<char> &variables );
		double discrete_subfunction(int subfunction_index, vec_t<char> &variables);
		double continuous_subfunction(int subfunction_index, vec_t<double> &c_variables);
		// void evaluationFunction( solution_mixed *solution );
		// double objectiveFunction( int objective_index, solution_mixed *solution );
		// double objectiveFunction( int objective_index, vec_t<double> &fitness_buffers );
		// double constraintFunction( solution_mixed *solution );
		// from gbo_fitness.hpp
		void evaluationFunction( solution_t<char> *solution );
		// void partialEvaluationFunction( solution_t<char> *parent, partial_solution_t<char> *solution );
		double objectiveFunction( int objective_index, solution_t<char> *solution );
		double objectiveFunction( int objective_index, vec_t<double> &fitness_buffers );
		double constraintFunction( solution_t<char> *solution );

		// from fitness.hpp;
		// void evaluate( solution_t<T> *solution ); // Actually not necessary to override this function, I think
		// void evaluatePartialSolution( solution_t<char> *parent, partial_solution_t<char> *solution );
		// bool betterFitness( solution_t<T> *sol_x, solution_t<T> *sol_y ); // Also not necessary to override this function, I think
		// void evaluatePartialSolutionBlackBox( solution_t<char> *parent, partial_solution_t<char> *solution );

};

class DT5RotEllip_t: public GBOFitnessFunction_t<char>
{
	public:
		DT5RotEllip_t( int number_of_variables, int number_of_c_variables );
		int getNumberOfSubfunctions(); 
		vec_t<int> inputsToSubfunction( int subfunction_index );
		int number_of_c_variables;
		double getLowerRangeBound( int dimension );
		double getUpperRangeBound( int dimension );
		
	private:
		double subfunction( int subfunction_index, vec_t<char> &variables );
		double discrete_subfunction(int subfunction_index, vec_t<char> &variables);
		double continuous_subfunction(int subfunction_index, vec_t<double> &c_variables);
		// void evaluationFunction( solution_mixed *solution );
		// double objectiveFunction( int objective_index, solution_mixed *solution );
		// double objectiveFunction( int objective_index, vec_t<double> &fitness_buffers );
		// double constraintFunction( solution_mixed *solution );
		// from gbo_fitness.hpp
		void evaluationFunction( solution_t<char> *solution );
		// void partialEvaluationFunction( solution_t<char> *parent, partial_solution_t<char> *solution );
		double objectiveFunction( int objective_index, solution_t<char> *solution );
		double objectiveFunction( int objective_index, vec_t<double> &fitness_buffers );
		double constraintFunction( solution_t<char> *solution );
		vec_t<double> rotateParameters( vec_t<double> &parameters);
		// from fitness.hpp;
		// void evaluate( solution_t<T> *solution ); // Actually not necessary to override this function, I think
		// void evaluatePartialSolution( solution_t<char> *parent, partial_solution_t<char> *solution );
		// bool betterFitness( solution_t<T> *sol_x, solution_t<T> *sol_y ); // Also not necessary to override this function, I think
		// void evaluatePartialSolutionBlackBox( solution_t<char> *parent, partial_solution_t<char> *solution );


		double **initializeObjectiveRotationMatrix( double rotation_angle, int rotation_block_size );
		double *rotateVariables( double *variables, int num_variables, double **rotation_matrix );
		
};

// class deceptiveTrap_t: public GBOFitnessFunction_t<char>
// {
// 	public:
// 		deceptiveTrap_t( int number_of_variables, int trap_size );
// 		int getNumberOfSubfunctions(); 
// 		vec_t<int> inputsToSubfunction( int subfunction_index );
		
// 	private:
// 		int trap_size;
// 		double subfunction( int subfunction_index, vec_t<char> &variables );
// };

// class deceptiveTrapBBO_t: public BBOFitnessFunction_t<char>
// {
// 	public:
// 		deceptiveTrapBBO_t( int number_of_variables, int trap_size );
// 		double objectiveFunction( int objective_index, vec_t<char> &variables );
		
// 	private:
// 		int trap_size;
// };

// class maxCut_t: public GBOFitnessFunction_t<char>
// {
// 	public:
// 		maxCut_t( std::string input_file, std::string vtr_file );
// 		int getNumberOfSubfunctions(); 
// 		vec_t<int> inputsToSubfunction( int subfunction_index );
		
// 	private:
// 		vec_t<vec_t<int>> edges;
// 		vec_t<double> edge_weights;

// 		void readInputFile( std::string input_file );
// 		void readVTRFile( std::string input_file );
// 		double subfunction( int subfunction_index, vec_t<char> &variables );
// };

// class maxCutBBO_t: public BBOFitnessFunction_t<char>
// {
// 	public:
// 		maxCutBBO_t( std::string input_file, std::string vtr_file );
// 		double objectiveFunction( int objective_index, vec_t<char> &variables );
		
// 	private:
// 		vec_t<vec_t<int>> edges;
// 		vec_t<double> edge_weights;

// 		void readInputFile( std::string input_file );
// 		void readVTRFile( std::string input_file );
// };

}}
