/*-=-=-=-=-=-=-=-=-=-=-=-=-=-= Section Includes -=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
#include "gomea/src/fitness/benchmarks-mixed.hpp"
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
namespace gomea{
namespace fitness{

using namespace gomea;

DT5BlockNOTRotEllipWrongExponent_t::DT5BlockNOTRotEllipWrongExponent_t( int number_of_variables, int number_of_c_variables, double a ) : GBOFitnessFunction_t<int>(number_of_variables)
{
	this->name = "F5: DeceptiveTrap5-Block NOT rotated ellipse function (with cross-domain dependencies)";
    this->number_of_c_variables = number_of_c_variables;
	this->vtr = 10e-10;
	this->use_vtr = true;
	this->optimization_mode = opt_mode::MIN;
	this->initialize();
	rotation_matrix = initializeObjectiveRotationMatrix(45, k);

    for(int i = 0; i < pow(2, k); i++)
    {
        ellipsoid_centres.push_back(vec_t<double>());
        for(int j = 0; j < k; j++)
        {
            double val = ((gomea::utils::rng)() / (double)(gomea::utils::rng).max()) * 10 - 5;
            ellipsoid_centres[i].push_back(val);
        }
    }
    this->a = a;
}

double DT5BlockNOTRotEllipWrongExponent_t::getLowerRangeBound( int dimension )
{
	return( -1e308 );
}
		
double DT5BlockNOTRotEllipWrongExponent_t::getUpperRangeBound( int dimension )
{
	return( 1e308 );
}
		
int DT5BlockNOTRotEllipWrongExponent_t::getNumberOfSubfunctions()
{
	return number_of_variables;
}
		
vec_t<int> DT5BlockNOTRotEllipWrongExponent_t::inputsToSubfunction( int subfunction_index )
{
	vec_t<int> vec;
	vec.push_back(subfunction_index);
	return vec;
}

		
double DT5BlockNOTRotEllipWrongExponent_t::subfunction( int subfunction_index, vec_t<int> &variables )
{
	if(optimization_mode == opt_mode::MAX)
		return( variables[subfunction_index] == 1 ? 1.0 : 0.0 );
	else
		return( variables[subfunction_index] == 1 ? 0.0 : 1.0 );

    
}

// Assuming subfunction_index is referring to the group of 5 variables for DT5
double DT5BlockNOTRotEllipWrongExponent_t::discrete_subfunction(int subfunction_index, vec_t<int> &variables)
{
	int sum = 0;
	for(int i = 0; i < k; i++) {
		sum += variables[(subfunction_index*k) + i] == 1 ? 1 : 0;
	}

    if(optimization_mode == opt_mode::MAX)
		return( sum == k ? 1.0 : (k-1 - sum) / (double)k );
	else
		return( sum == k ? 0.0 : 1 - ((k-1 - sum) / (double)k) );
}

double DT5BlockNOTRotEllipWrongExponent_t::continuous_subfunction(int subfunction_index, vec_t<int> &variables, vec_t<double> &c_variables)
{
	// Convert k discrete (binary) variables to a decimal value to index the ellipsoid centres
    int decimal_value = 0;
    for(int i = 0; i < k; i++) {
        decimal_value += variables[(subfunction_index*k) + i] == 1 ? pow(2, k-1-i) : 0;
    }
    vec_t<double> ellipsoid_centre = ellipsoid_centres[decimal_value];

    double res = 0.0;
    for(int i = 0; i < k; i++) {
        res += pow(10, 6 * ((double)(i+subfunction_index*k) / (double)(number_of_c_variables - 1))) * pow(c_variables[subfunction_index * k + i] - ellipsoid_centre[i], 2);
    }
    return res;
}

void DT5BlockNOTRotEllipWrongExponent_t::evaluationFunction( solution_t<int> *solution )
{
    // solution = static_cast<solution_mixed*>(solution);
	solution_mixed *solution_mix = static_cast<solution_mixed*>(solution);

    // For this function, there has to be an equal number of discrete and continuous variables
    assert(solution_mix->getNumberOfVariables() == solution_mix->getNumberOfCVariables());
	
    solution_mix->initFitnessBuffers(getNumberOfFitnessBuffers());
	solution_mix->clearFitnessBuffers();

	for( int i = 0; i < solution_mix->getNumberOfVariables() / k; i++ )
	{
		int buffer_index = this->getIndexOfFitnessBuffer(i);
		double ftrap_sub = discrete_subfunction(i, solution_mix->variables);
        double fsub_ellipse = continuous_subfunction(i, solution_mix->variables, solution_mix->c_variables);
        double discrete_part = 1 + pow(10, a) * ftrap_sub;
        double continuous_part = 1 + fsub_ellipse;
        double fsub = discrete_part * continuous_part;

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

double DT5BlockNOTRotEllipWrongExponent_t::objectiveFunction( int objective_index, solution_t<int> *solution )
{
    return objectiveFunction(objective_index,solution->fitness_buffers);
}

double DT5BlockNOTRotEllipWrongExponent_t::objectiveFunction( int objective_index, vec_t<double> &fitness_buffers )
{
    return fitness_buffers[objective_index];
}

double DT5BlockNOTRotEllipWrongExponent_t::constraintFunction( solution_t<int> *solution )
{
    return 0.0;
}


double **DT5BlockNOTRotEllipWrongExponent_t::initializeObjectiveRotationMatrix( double rotation_angle, int rotation_block_size )
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

double *DT5BlockNOTRotEllipWrongExponent_t::rotateVariables( double *variables, int num_variables, double **rotation_matrix )
{
	double *rotated_variables = gomea::utils::matrixVectorMultiplication( rotation_matrix, variables, num_variables, num_variables );
    return( rotated_variables );
}

}}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
