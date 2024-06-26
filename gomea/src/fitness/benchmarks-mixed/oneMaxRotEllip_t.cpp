/*-=-=-=-=-=-=-=-=-=-=-=-=-=-= Section Includes -=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
#include "gomea/src/fitness/benchmarks-mixed.hpp"
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
namespace gomea{
namespace fitness{

using namespace gomea;

oneMaxRotEllip_t::oneMaxRotEllip_t( int number_of_variables, int number_of_c_variables ) : GBOFitnessFunction_t<int>(number_of_variables)
{
	this->name = "F2: OneMax-RotatedEllipse function";
    this->number_of_c_variables = number_of_c_variables;
	this->vtr = 10e-10;
	this->use_vtr = true;
	this->optimization_mode = opt_mode::MIN;
	this->initialize();
	rotation_matrix = initializeObjectiveRotationMatrix(45, number_of_c_variables);
}

double oneMaxRotEllip_t::getLowerRangeBound( int dimension )
{
	return( -1e308 );
}
		
double oneMaxRotEllip_t::getUpperRangeBound( int dimension )
{
	return( 1e308 );
}
		
int oneMaxRotEllip_t::getNumberOfSubfunctions()
{
	return number_of_variables;
}
		
vec_t<int> oneMaxRotEllip_t::inputsToSubfunction( int subfunction_index )
{
	vec_t<int> vec;
	vec.push_back(subfunction_index);
	return vec;
}

		
double oneMaxRotEllip_t::subfunction( int subfunction_index, vec_t<int> &variables )
{
	if(optimization_mode == opt_mode::MAX)
		return( variables[subfunction_index] == 1 ? 1.0 : 0.0 );
	else
		return( variables[subfunction_index] == 1 ? 0.0 : 1.0 );

    
}

// Assuming subfunction_index is referring to the group of 5 variables for DT5
double oneMaxRotEllip_t::discrete_subfunction(int subfunction_index, vec_t<int> &variables)
{
	if(optimization_mode == opt_mode::MAX)
		return( variables[subfunction_index] == 1 ? 1.0 : 0.0 );
	else
		return( variables[subfunction_index] == 1 ? 0.0 : 1.0 );
}

double oneMaxRotEllip_t::continuous_subfunction(int subfunction_index, vec_t<double> &c_variables)
{
	return pow( 10.0, 6.0*(((double) (subfunction_index))/((double) (number_of_c_variables-1)))) * (c_variables[subfunction_index]) * (c_variables[subfunction_index]);
}

void oneMaxRotEllip_t::evaluationFunction( solution_t<int> *solution )
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

	// Use rotated variables
	double *rotated_variables_ptr = rotateVariables(solution_mix->c_variables.data(), number_of_c_variables, rotation_matrix);
	vec_t<double> rotated_variables(rotated_variables_ptr, rotated_variables_ptr + number_of_c_variables);

    for( int i = 0; i < solution_mix->getNumberOfCVariables(); i++ )
    {
        int buffer_index = this->getIndexOfFitnessBuffer(i);
        double fsub = continuous_subfunction(i, rotated_variables);
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


double oneMaxRotEllip_t::objectiveFunction( int objective_index, solution_t<int> *solution )
{
    return objectiveFunction(objective_index,solution->fitness_buffers);
}

double oneMaxRotEllip_t::objectiveFunction( int objective_index, vec_t<double> &fitness_buffers )
{
    return fitness_buffers[objective_index];
}

double oneMaxRotEllip_t::constraintFunction( solution_t<int> *solution )
{
    return 0.0;
}


double **oneMaxRotEllip_t::initializeObjectiveRotationMatrix( double rotation_angle, int rotation_block_size )
{
    if( rotation_angle == 0.0 )
        return NULL;

    double **matrix = new double*[rotation_block_size];
    for( int i = 0; i < rotation_block_size; i++ )
        matrix[i] = new double[rotation_block_size];

    double **rotation_matrix = new double*[rotation_block_size];
    for( int i = 0; i < rotation_block_size; i++ )
        rotation_matrix[i] = new double[rotation_block_size];

    /* Initialize the rotation matrix to the identity matrix */
    for( int i = 0; i < rotation_block_size; i++ )
    {
        for( int j = 0; j < rotation_block_size; j++ )
            rotation_matrix[i][j] = 0.0;
        rotation_matrix[i][i] = 1.0;
    }

    /* Construct all rotation matrices (quadratic number) and multiply */
    double theta     = (rotation_angle/180.0)*M_PI;
    double cos_theta = cos( theta );
    double sin_theta = sin( theta );
    for( int index0 = 0; index0 < rotation_block_size-1; index0++ )
    {
        for( int index1 = index0+1; index1 < rotation_block_size; index1++ )
        {
            for( int i = 0; i < rotation_block_size; i++ )
            {
                for( int j = 0; j < rotation_block_size; j++ )
                    matrix[i][j] = 0.0;
                matrix[i][i] = 1.0;
            }
            matrix[index0][index0] = cos_theta;
            matrix[index0][index1] = -sin_theta;
            matrix[index1][index0] = sin_theta;
            matrix[index1][index1] = cos_theta;
	
            double **product = gomea::utils::matrixMatrixMultiplication( matrix, rotation_matrix, rotation_block_size, rotation_block_size, rotation_block_size );
            for( int i = 0; i < rotation_block_size; i++ )
                for( int j = 0; j < rotation_block_size; j++ )
                    rotation_matrix[i][j] = product[i][j];

            for( int i = 0; i < rotation_block_size; i++ )
                delete[] product[i];
            delete[] product;
        }
    }

    for( int i = 0; i < rotation_block_size; i++ )
	{
        delete[] matrix[i];
	}
	delete[] matrix;

	return( rotation_matrix );
}

double *oneMaxRotEllip_t::rotateVariables( double *variables, int num_variables, double **rotation_matrix )
{
	double *rotated_variables = gomea::utils::matrixVectorMultiplication( rotation_matrix, variables, num_variables, num_variables );
    return( rotated_variables );
}

}}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
