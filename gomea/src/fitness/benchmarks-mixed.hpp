#pragma once

#include "gomea/src/fitness/fitness.hpp"
#include "gomea/src/fitness/gbo_fitness.hpp"
#include "gomea/src/fitness/bbo_fitness.hpp"
#include "gomea/src/common/solution_mixed.hpp"
#include "gomea/src/common/partial_solution.hpp"
#include "gomea/src/common/gomea_defs.hpp"
#include "gomea/src/utils/tools.hpp"
#include "gomea/src/fitness/density.h"
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

namespace gomea{
namespace fitness{

class oneMaxSphere_t: public GBOFitnessFunction_t<int>
{
	public:
		oneMaxSphere_t( int number_of_variables, int number_of_c_variables );
		int getNumberOfSubfunctions(); 
		vec_t<int> inputsToSubfunction( int subfunction_index );
		int number_of_c_variables;
		double getLowerRangeBound( int dimension );
		double getUpperRangeBound( int dimension );
		
	private:
		double subfunction( int subfunction_index, vec_t<int> &variables );
		double discrete_subfunction(int subfunction_index, vec_t<int> &variables);
		double continuous_subfunction(int subfunction_index, vec_t<double> &c_variables);
		void evaluationFunction( solution_t<int> *solution );
		double objectiveFunction( int objective_index, solution_t<int> *solution );
		double objectiveFunction( int objective_index, vec_t<double> &fitness_buffers );
		double constraintFunction( solution_t<int> *solution );

};

class oneMaxRotEllip_t: public GBOFitnessFunction_t<int>
{
	public:
		oneMaxRotEllip_t( int number_of_variables, int number_of_c_variables );
		int getNumberOfSubfunctions(); 
		vec_t<int> inputsToSubfunction( int subfunction_index );
		int number_of_c_variables;
		double getLowerRangeBound( int dimension );
		double getUpperRangeBound( int dimension );
		
	private:
		double subfunction( int subfunction_index, vec_t<int> &variables );
		double discrete_subfunction(int subfunction_index, vec_t<int> &variables);
		double continuous_subfunction(int subfunction_index, vec_t<double> &c_variables);
		void evaluationFunction( solution_t<int> *solution );
		double objectiveFunction( int objective_index, solution_t<int> *solution );
		double objectiveFunction( int objective_index, vec_t<double> &fitness_buffers );
		double constraintFunction( solution_t<int> *solution );
		vec_t<double> rotateParameters( vec_t<double> &parameters);


		double **initializeObjectiveRotationMatrix( double rotation_angle, int rotation_block_size );
		double *rotateVariables( double *variables, int num_variables, double **rotation_matrix );
		
};

class DT5Sphere_t: public GBOFitnessFunction_t<int>
{
	public:
		DT5Sphere_t( int number_of_variables, int number_of_c_variables );
		int getNumberOfSubfunctions(); 
		vec_t<int> inputsToSubfunction( int subfunction_index );
		int number_of_c_variables;
		double getLowerRangeBound( int dimension );
		double getUpperRangeBound( int dimension );
		double getInitUpperRangeBound( int dimension );
		double getInitLowerRangeBound( int dimension );
		
	private:
		double subfunction( int subfunction_index, vec_t<int> &variables );
		double discrete_subfunction(int subfunction_index, vec_t<int> &variables);
		double continuous_subfunction(int subfunction_index, vec_t<double> &c_variables);
		void evaluationFunction( solution_t<int> *solution );
		double objectiveFunction( int objective_index, solution_t<int> *solution );
		double objectiveFunction( int objective_index, vec_t<double> &fitness_buffers );
		double constraintFunction( solution_t<int> *solution );

};

class DT5RotEllip_t: public GBOFitnessFunction_t<int>
{
	public:
		DT5RotEllip_t( int number_of_variables, int number_of_c_variables );
		int getNumberOfSubfunctions(); 
		vec_t<int> inputsToSubfunction( int subfunction_index );
		int number_of_c_variables;
		double getLowerRangeBound( int dimension );
		double getUpperRangeBound( int dimension );
		
	private:
		double subfunction( int subfunction_index, vec_t<int> &variables );
		double discrete_subfunction(int subfunction_index, vec_t<int> &variables);
		double continuous_subfunction(int subfunction_index, vec_t<double> &c_variables);
		void evaluationFunction( solution_t<int> *solution );
		double objectiveFunction( int objective_index, solution_t<int> *solution );
		double objectiveFunction( int objective_index, vec_t<double> &fitness_buffers );
		double constraintFunction( solution_t<int> *solution );
		vec_t<double> rotateParameters( vec_t<double> &parameters);

		double **initializeObjectiveRotationMatrix( double rotation_angle, int rotation_block_size );
		double *rotateVariables( double *variables, int num_variables, double **rotation_matrix );
		
};

class DT5BlockRotEllip_t: public GBOFitnessFunction_t<int>
{
	public:
		DT5BlockRotEllip_t( int number_of_variables, int number_of_c_variables, double a = 1.1 );
		int getNumberOfSubfunctions(); 
		vec_t<int> inputsToSubfunction( int subfunction_index );
		int number_of_c_variables;
		double getLowerRangeBound( int dimension );
		double getUpperRangeBound( int dimension );
		
	private:
		double subfunction( int subfunction_index, vec_t<int> &variables );
		double discrete_subfunction(int subfunction_index, vec_t<int> &variables);
		double continuous_subfunction(int subfunction_index, vec_t<int> &variables, vec_t<double> &c_variables);
		void evaluationFunction( solution_t<int> *solution ) override;
		double objectiveFunction( int objective_index, solution_t<int> *solution );
		double objectiveFunction( int objective_index, vec_t<double> &fitness_buffers );
		double constraintFunction( solution_t<int> *solution );
		vec_t<double> rotateParameters( vec_t<double> &parameters);

		double **initializeObjectiveRotationMatrix( double rotation_angle, int rotation_block_size );
		double *rotateVariables( double *variables, int num_variables, double **rotation_matrix );
		vec_t<vec_t<double> > ellipsoid_centres;
		double a;
		int k = 5;

};

class DT5BlockNOTRotEllip_t: public GBOFitnessFunction_t<int>
{
	public:
		DT5BlockNOTRotEllip_t( int number_of_variables, int number_of_c_variables, double a = 1.1 );
		int getNumberOfSubfunctions(); 
		vec_t<int> inputsToSubfunction( int subfunction_index );
		int number_of_c_variables;
		double getLowerRangeBound( int dimension );
		double getUpperRangeBound( int dimension );
		
	private:
		double subfunction( int subfunction_index, vec_t<int> &variables );
		double discrete_subfunction(int subfunction_index, vec_t<int> &variables);
		double continuous_subfunction(int subfunction_index, vec_t<int> &variables, vec_t<double> &c_variables);
		void evaluationFunction( solution_t<int> *solution ) override;
		double objectiveFunction( int objective_index, solution_t<int> *solution );
		double objectiveFunction( int objective_index, vec_t<double> &fitness_buffers );
		double constraintFunction( solution_t<int> *solution );
		vec_t<double> rotateParameters( vec_t<double> &parameters);

		double **initializeObjectiveRotationMatrix( double rotation_angle, int rotation_block_size );
		double *rotateVariables( double *variables, int num_variables, double **rotation_matrix );
		vec_t<vec_t<double> > ellipsoid_centres;
		double a;
		int k = 5;

};

class DT5BlockFullRotEllipWrongExponent_t: public GBOFitnessFunction_t<int>
{
	public:
		DT5BlockFullRotEllipWrongExponent_t( int number_of_variables, int number_of_c_variables, double a = 1.1 );
		int getNumberOfSubfunctions(); 
		vec_t<int> inputsToSubfunction( int subfunction_index );
		int number_of_c_variables;
		double getLowerRangeBound( int dimension );
		double getUpperRangeBound( int dimension );
		
	private:
		double subfunction( int subfunction_index, vec_t<int> &variables );
		double discrete_subfunction(int subfunction_index, vec_t<int> &variables);
		double continuous_subfunction(int subfunction_index, vec_t<int> &variables, vec_t<double> &c_variables);
		void evaluationFunction( solution_t<int> *solution ) override;
		double objectiveFunction( int objective_index, solution_t<int> *solution );
		double objectiveFunction( int objective_index, vec_t<double> &fitness_buffers );
		double constraintFunction( solution_t<int> *solution );
		vec_t<double> rotateParameters( vec_t<double> &parameters);

		double **initializeObjectiveRotationMatrix( double rotation_angle, int rotation_block_size );
		double *rotateVariables( double *variables, int num_variables, double **rotation_matrix );
		vec_t<vec_t<double> > ellipsoid_centres;
		double a;
		int k = 5;

};

class DT5BlockNOTRotEllipWrongExponent_t: public GBOFitnessFunction_t<int>
{
	public:
		DT5BlockNOTRotEllipWrongExponent_t( int number_of_variables, int number_of_c_variables, double a = 1.1 );
		int getNumberOfSubfunctions(); 
		vec_t<int> inputsToSubfunction( int subfunction_index );
		int number_of_c_variables;
		double getLowerRangeBound( int dimension );
		double getUpperRangeBound( int dimension );
		
	private:
		double subfunction( int subfunction_index, vec_t<int> &variables );
		double discrete_subfunction(int subfunction_index, vec_t<int> &variables);
		double continuous_subfunction(int subfunction_index, vec_t<int> &variables, vec_t<double> &c_variables);
		void evaluationFunction( solution_t<int> *solution ) override;
		double objectiveFunction( int objective_index, solution_t<int> *solution );
		double objectiveFunction( int objective_index, vec_t<double> &fitness_buffers );
		double constraintFunction( solution_t<int> *solution );
		vec_t<double> rotateParameters( vec_t<double> &parameters);

		double **initializeObjectiveRotationMatrix( double rotation_angle, int rotation_block_size );
		double *rotateVariables( double *variables, int num_variables, double **rotation_matrix );
		vec_t<vec_t<double> > ellipsoid_centres;
		double a;
		int k = 5;

};

class DT5BlockRotEllipCentersZero_t: public GBOFitnessFunction_t<int>
{
	public:
		DT5BlockRotEllipCentersZero_t( int number_of_variables, int number_of_c_variables, double a = 1.1 );
		int getNumberOfSubfunctions(); 
		vec_t<int> inputsToSubfunction( int subfunction_index );
		int number_of_c_variables;
		double getLowerRangeBound( int dimension );
		double getUpperRangeBound( int dimension );
		
	private:
		double subfunction( int subfunction_index, vec_t<int> &variables );
		double discrete_subfunction(int subfunction_index, vec_t<int> &variables);
		double continuous_subfunction(int subfunction_index, vec_t<int> &variables, vec_t<double> &c_variables);
		void evaluationFunction( solution_t<int> *solution ) override;
		double objectiveFunction( int objective_index, solution_t<int> *solution );
		double objectiveFunction( int objective_index, vec_t<double> &fitness_buffers );
		double constraintFunction( solution_t<int> *solution );
		vec_t<double> rotateParameters( vec_t<double> &parameters);

		double **initializeObjectiveRotationMatrix( double rotation_angle, int rotation_block_size );
		double *rotateVariables( double *variables, int num_variables, double **rotation_matrix );
		vec_t<vec_t<double> > ellipsoid_centres;
		double a;
		int k = 5;

};


class BNStructureLearning: public GBOFitnessFunction_t<int>
{
	public:
		BNStructureLearning( int number_of_variables, int number_of_c_variables );
		BNStructureLearning( int number_of_variables, int number_of_c_variables, int problem_index, const shared_ptr<DataStructure<double>> &data, size_t max_number_of_parents, size_t max_number_of_discretizations, bool transformCVariables = false);
		~BNStructureLearning();
		int getNumberOfSubfunctions(); 
		vec_t<int> inputsToSubfunction( int subfunction_index );
		int number_of_c_variables;
		double getLowerRangeBound( int dimension );
		double getUpperRangeBound( int dimension );
		double getInitUpperRangeBound( int dimension );
		double getInitLowerRangeBound( int dimension );
		const Density *getDensity() const;
		
	private:
		Density *density;
		bool transformCVariables;
		double subfunction( int subfunction_index, vec_t<int> &variables );
		double discrete_subfunction(int subfunction_index, vec_t<int> &variables);
		double continuous_subfunction(int subfunction_index, vec_t<double> &c_variables);
		void evaluationFunction( solution_t<int> *solution );
		double objectiveFunction( int objective_index, solution_t<int> *solution );
		double objectiveFunction( int objective_index, vec_t<double> &fitness_buffers );
		double constraintFunction( solution_t<int> *solution );

};

}}
