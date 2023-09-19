/*-=-=-=-=-=-=-=-=-=-=-=-=-=-= Section Includes -=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
#include "gomea/src/fitness/benchmarks-mixed.hpp"
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
namespace gomea{
namespace fitness{

using namespace gomea;

DT5BlockRotEllipBBO_t::DT5BlockRotEllipBBO_t( int number_of_variables, int number_of_c_variables, double a ) : BBOFitnessFunction_t<char>(number_of_variables)
{
	this->name = "(BBO) F5: DeceptiveTrap5-BlockRotatedEllipse function (with cross-domain dependencies)";
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

double DT5BlockRotEllipBBO_t::getLowerRangeBound( int dimension )
{
	return( -1e308 );
}
		
double DT5BlockRotEllipBBO_t::getUpperRangeBound( int dimension )
{
	return( 1e308 );
}

// Assuming subfunction_index is referring to the group of 5 variables for DT5
double DT5BlockRotEllipBBO_t::discrete_subfunction(int subfunction_index, vec_t<char> &variables)
{
	int sum = 0;
	for(int i = 0; i < k; i++) {
		sum += variables[(subfunction_index*k) + i] == '\001' ? 1 : 0;
	}

    if(optimization_mode == opt_mode::MAX)
		return( sum == k ? 1.0 : (k-1 - sum) / (double)k );
	else
		return( sum == k ? 0.0 : 1 - ((k-1 - sum) / (double)k) );
}

double DT5BlockRotEllipBBO_t::continuous_subfunction(int subfunction_index, vec_t<char> &variables, vec_t<double> &c_variables)
{
	// Convert k discrete (binary) variables to a decimal value to index the ellipsoid centres
    int decimal_value = 0;
    for(int i = 0; i < k; i++) {
        decimal_value += variables[(subfunction_index*k) + i] == '\001' ? pow(2, k-1-i) : 0;
    }
    vec_t<double> ellipsoid_centre = ellipsoid_centres[decimal_value];

    // Use rotated variables
    vec_t<double> sliced_vars(c_variables.begin() + subfunction_index * k, c_variables.begin() + (subfunction_index + 1) * k);
	double *rotated_variables_ptr = rotateVariables(sliced_vars.data(), k, rotation_matrix);
	vec_t<double> rotated_variables(rotated_variables_ptr, rotated_variables_ptr + k);

    double res = 0.0;
    for(int i = 0; i < k; i++) {
        res += pow(10, 6 * ((double)i / (double)(k - 1))) * pow(rotated_variables[i] - ellipsoid_centre[i], 2);
    }
    return res;
}

void DT5BlockRotEllipBBO_t::evaluationFunction( solution_t<char> *solution )
{
    // solution = static_cast<solution_mixed*>(solution);
	solution_mixed *solution_mix = static_cast<solution_mixed*>(solution);

    // For this function, there has to be an equal number of discrete and continuous variables
    assert(solution_mix->getNumberOfVariables() == solution_mix->getNumberOfCVariables());
	

	double fitness_value = 0.0;
	for( int i = 0; i < solution_mix->getNumberOfVariables() / k; i++ )
	{
		double ftrap_sub = discrete_subfunction(i, solution_mix->variables);
        double fsub_ellipse = continuous_subfunction(i, solution_mix->variables, solution_mix->c_variables);
        double discrete_part = 1 + pow(10, a) * ftrap_sub;
        double continuous_part = 1 + fsub_ellipse;
        fitness_value += discrete_part * continuous_part;
	}

	// Assuming single-objective
	solution_mix->setObjectiveValue(fitness_value);

	double fcons = constraintFunction(solution);
	solution_mix->setConstraintValue(fcons);

	this->full_number_of_evaluations++;
	this->number_of_evaluations++;
}

double DT5BlockRotEllipBBO_t::objectiveFunction( int objective_index, solution_t<char> *solution )
{
	solution_mixed *solution_mix = static_cast<solution_mixed*>(solution);
	// For this function, there has to be an equal number of discrete and continuous variables
    assert(solution_mix->getNumberOfVariables() == solution_mix->getNumberOfCVariables());

	double fitness_total = 0.0;

	for( int i = 0; i < solution_mix->getNumberOfVariables() / k; i++ )
	{
		double ftrap_sub = discrete_subfunction(i, solution_mix->variables);
        double fsub_ellipse = continuous_subfunction(i, solution_mix->variables, solution_mix->c_variables);
        double discrete_part = 1 + pow(10, a) * ftrap_sub;
        double continuous_part = 1 + fsub_ellipse;
        fitness_total += discrete_part * continuous_part;
	}

	return fitness_total;
}

double DT5BlockRotEllipBBO_t::objectiveFunction( int objective_index, vec_t<char> &variables )
{
	cerr << "Wrong objectiveFunction used: DT5BlockRotEllipBBO_t::objectiveFunction( int objective_index, vec_t<char> &variables )\n";
	cerr << "Use: DT5BlockRotEllipBBO_t::objectiveFunction( int objective_index, solution_t<char> *solution )" << endl;
	exit( 0 );
}

double DT5BlockRotEllipBBO_t::constraintFunction( solution_t<char> *solution )
{
	return 0.0;
}

double **DT5BlockRotEllipBBO_t::initializeObjectiveRotationMatrix( double rotation_angle, int rotation_block_size )
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

double *DT5BlockRotEllipBBO_t::rotateVariables( double *variables, int num_variables, double **rotation_matrix )
{
	double *rotated_variables = gomea::utils::matrixVectorMultiplication( rotation_matrix, variables, num_variables, num_variables );
    return( rotated_variables );
}

}}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
