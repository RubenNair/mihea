
using namespace std;

#include "gomea/src/mixed_integer/iamalgam.hpp"

#include "gomea/src/common/solution_BN.hpp"
#include "gomea/src/fitness/benchmarks-mixed.hpp"

namespace gomea{
namespace mixedinteger{

iamalgam::iamalgam()
{
    return;

}

iamalgam::iamalgam(Config *config_): config(config_)
{
    population_size       = config->basePopulationSize;
    problemInstance       = config->fitness;
    number_of_parameters  = config->numberOfcVariables;
    number_of_populations = 1;
    haveNextNextGaussian = false;

}

iamalgam::iamalgam(Config *config_, Solutionset *population_) : iamalgam(config_)
{
  population = population_;
  population_size = population_->size();
}

iamalgam::~iamalgam()
{}


void iamalgam::initialize() {
    number_of_generations = 0;

    if( number_of_starts == 1 )
    number_of_evaluations = 0;

  alpha_AMS = 0.5*tau*(((double) population_size)/((double) (population_size-1)));
  delta_AMS = 2.0;

  initializeMemory();

  initializeDistributionMultipliers();

  computeRanks();

  eta_p = 1.0-exp(-1.2*pow(selection_size,0.31)/pow(number_of_parameters,0.50));
  eta_s = 1.0-exp(-1.1*pow(selection_size,1.20)/pow(number_of_parameters,1.60));
    
}

void iamalgam::initializeMemory() {
    int i, j;
  assert(number_of_populations == 1);

  selection_size = (int) (tau*(population_size));
  
  populations_terminated           = (bool *) malloc( number_of_populations*sizeof( bool ) );
  no_improvement_stretch           = (int *) malloc( number_of_populations*sizeof( int ) );
  ranks                            = (double **) malloc( number_of_populations*sizeof( double * ) );
  
  selections.resize(selection_size);
  for(int i = 0; i < selection_size; i++) 
  {
    if(config->useBN) {
      // If useBN flag is true, assume we're dealing with BNs. Static cast problemInstance to BNStructureLearning
      fitness::BNStructureLearning *BNproblemInstance = dynamic_cast<fitness::BNStructureLearning*>(problemInstance);

      vec_t<double> maxValuesData, minValuesData;
      vec_t<vec_t<double>> data = config->data->getDataMatrix().getRawMatrix();
      std::tie(maxValuesData, minValuesData) = findMaxAndMinValuesInData(data);

      selections[i] = new solution_BN(config->numberOfVariables, config->alphabetSize, config->numberOfcVariables, config->data->getColumnType(), BNproblemInstance->getDensity()->getOriginalData()->getColumnNumberOfClasses(), config->discretization_policy_index, config->maxParents, config->maxInstantiations, problemInstance, maxValuesData, minValuesData, config->lower_user_range, config->upper_user_range, config->data);
      selections[i]->initObjectiveValues(problemInstance->number_of_objectives);
    } else {
      selections[i] = new solution_mixed(config->numberOfVariables, config->alphabetSize, config->numberOfcVariables, problemInstance);
      selections[i]->randomInit(&gomea::utils::rng);
      selections[i]->initObjectiveValues(problemInstance->number_of_objectives);
    }
    if(population->solutions[0]->fitness_buffers.size() > 0)
      {
        selections[i]->initFitnessBuffers(problemInstance->number_of_objectives);
      }
  }

  mean_vectors                     = (double **) malloc( number_of_populations*sizeof( double * ) );
  mean_vectors_previous            = (double **) malloc( number_of_populations*sizeof( double * ) );
  covariance_matrices              = (double ***) malloc( number_of_populations*sizeof( double ** ) );
  generational_covariance_matrices = (double ***) malloc( number_of_populations*sizeof( double ** ) );
  aggregated_covariance_matrices   = (double ***) malloc( number_of_populations*sizeof( double ** ) );
  cholesky_factors_lower_triangle  = (double ***) malloc( number_of_populations*sizeof( double ** ) );
  ams_vectors                      = (double **) malloc( number_of_populations*sizeof( double * ) );

  for( i = 0; i < number_of_populations; i++ )
  {

    populations_terminated[i] = 0;

    no_improvement_stretch[i] = 0;

    ranks[i] = (double *) malloc( population_size*sizeof( double ) );

    mean_vectors[i] = (double *) malloc( number_of_parameters*sizeof( double ) );

    mean_vectors_previous[i] = (double *) malloc( number_of_parameters*sizeof( double ) );

    covariance_matrices[i] = (double **) malloc( number_of_parameters*sizeof( double * ) );
    for( j = 0; j < number_of_parameters; j++ )
      covariance_matrices[i][j] = (double *) malloc( number_of_parameters*sizeof( double ) );
    
    generational_covariance_matrices[i] = (double **) malloc( number_of_parameters*sizeof( double * ) );
    for( j = 0; j < number_of_parameters; j++ )
      generational_covariance_matrices[i][j] = (double *) malloc( number_of_parameters*sizeof( double ) );
    
    aggregated_covariance_matrices[i] = (double **) malloc( number_of_parameters*sizeof( double * ) );
    for( j = 0; j < number_of_parameters; j++ )
      aggregated_covariance_matrices[i][j] = (double *) malloc( number_of_parameters*sizeof( double ) );

    cholesky_factors_lower_triangle[i] = NULL;

    ams_vectors[i] = (double *) malloc( number_of_parameters*sizeof( double ) );
  }

}

void iamalgam::initializeDistributionMultipliers() 
{
  int i;

  distribution_multipliers = (double *) malloc( number_of_populations*sizeof( double ) );
  for( i = 0; i < number_of_populations; i++ )
    distribution_multipliers[i] = 1.0;

  samples_drawn_from_normal = (int *) malloc( number_of_populations*sizeof( int ) );
  out_of_bounds_draws       = (int *) malloc( number_of_populations*sizeof( int ) );

  distribution_multiplier_increase = 1.0/distribution_multiplier_decrease;
}


void iamalgam::computeRanks()
{
  int i;

  for( i = 0; i < number_of_populations; i++ )
    if( !populations_terminated[i] )
      computeRanksForOnePopulation( i );
}

void iamalgam::computeRanksForOnePopulation(int population_index)
{
  int i, *sorted, rank;
  double *objective_values, *constraint_values;

  objective_values  = (double *) malloc( population_size*sizeof( double ) );
  constraint_values = (double *) malloc( population_size*sizeof( double ) );

  for( i = 0; i < population_size; i++)
  {
    objective_values[i] = population->solutions[i]->getObjectiveValue();
    constraint_values[i] = population->solutions[i]->getConstraintValue();
  }

  sorted = mergeSortFitness(  objective_values, constraint_values, population_size );
  
  rank                               = 0;
  ranks[population_index][sorted[0]] = rank;
  for( i = 1; i < population_size; i++ )
  {
    if( population->solutions[sorted[i]]->getObjectiveValue() != population->solutions[sorted[i-1]]->getObjectiveValue() )
      rank++;

    ranks[population_index][sorted[i]] = rank;
  }

  free( sorted );
  free( objective_values );
  free( constraint_values );
}

int *iamalgam::mergeSortFitness(double *objectives, double *constraints, int number_of_solutions)
{
  int i, *sorted, *tosort;

  sorted = (int *) malloc( number_of_solutions * sizeof( int ) );
  tosort = (int *) malloc( number_of_solutions * sizeof( int ) );
  for( i = 0; i < number_of_solutions; i++ )
    tosort[i] = i;

  if( number_of_solutions == 1 )
    sorted[0] = 0;
  else
    mergeSortFitnessWithinBounds( objectives, constraints, sorted, tosort, 0, number_of_solutions-1 );

  free( tosort );

  return( sorted );
}

void iamalgam::mergeSortFitnessWithinBounds(double *objectives, double *constraints, int *sorted, int *tosort, int p, int q)
{
  int r;

  if( p < q )
  {
    r = (p + q) / 2;
    mergeSortFitnessWithinBounds( objectives, constraints, sorted, tosort, p, r );
    mergeSortFitnessWithinBounds( objectives, constraints, sorted, tosort, r+1, q );
    mergeSortFitnessMerge( objectives, constraints, sorted, tosort, p, r+1, q );
  }
}

void iamalgam::mergeSortFitnessMerge(double *objectives, double *constraints, int *sorted, int *tosort, int p, int r, int q)
  {
    int i, j, k, first;

  i = p;
  j = r;
  for( k = p; k <= q; k++ )
  {
    first = 0;
    if( j <= q )
    {
      if( i < r )
      {
        if( betterFitness( objectives[tosort[i]], constraints[tosort[i]],
                           objectives[tosort[j]], constraints[tosort[j]] ) )
          first = 1;
      }
    }
    else
      first = 1;

    if( first )
    {
      sorted[k] = tosort[i];
      i++;
    }
    else
    {
      sorted[k] = tosort[j];
      j++;
    }
  }

  for( k = p; k <= q; k++ )
    tosort[k] = sorted[k];
}

bool iamalgam::betterFitness(double objective_value_x, double constraint_value_x, double objective_value_y, double constraint_value_y)
{
  bool result;


  result = false;

  // NOTE: only focussing on non-constrained problems for now, so ignoring constraint values.
  if( objective_value_x < objective_value_y ) // NOTE: assuming minimization
        result = true;

  return( result );
}

bool iamalgam::checkTerminationConditionForRunOnce() 
{
  bool allTrue;
  int   i;

  if( checkNumberOfEvaluationsTerminationCondition() )
    return( true );
  
  if( use_vtr )
  {
    if( checkVTRTerminationCondition() )
      return( true );
  }

  checkFitnessVarianceTermination();

  checkDistributionMultiplierTerminationCondition();

  allTrue = true;
  for( i = 0; i < number_of_populations; i++ )
  {
    if( !populations_terminated[i] )
    {
      allTrue = false;
      break;
    }
  }

  return( allTrue );
}

bool iamalgam::checkNumberOfEvaluationsTerminationCondition()
{
  if( number_of_evaluations >= maximum_number_of_evaluations )
    return( true );

  return( false );
}

bool iamalgam::checkVTRTerminationCondition()
{
  int population_of_best, index_of_best;

  determineBestSolutionInCurrentPopulations( &population_of_best, &index_of_best );

  if( population->solutions[index_of_best]->getConstraintValue() == 0 && population->solutions[index_of_best]->getObjectiveValue() <= problemInstance->vtr )
    return( true );

  return( false );
}

void iamalgam::checkFitnessVarianceTermination()
{
  int i;

  for( i = 0; i < number_of_populations; i++ )
  {
    if( !populations_terminated[i] )
      if( checkFitnessVarianceTerminationSinglePopulation( i ) )
        populations_terminated[i] = 1;
  }
}

bool iamalgam::checkFitnessVarianceTerminationSinglePopulation( int population_index )
{
  int    i;
  double objective_avg, objective_var;
  
  objective_avg = 0.0;
  for( i = 0; i < population_size; i++ )
    objective_avg  += population->solutions[i]->getObjectiveValue();
  objective_avg = objective_avg / ((double) population_size);

  objective_var = 0.0;
  for( i = 0; i < population_size; i++ )
    objective_var  += (population->solutions[i]->getObjectiveValue()-objective_avg)*(population->solutions[i]->getObjectiveValue()-objective_avg);
  objective_var = objective_var / ((double) population_size);

  if( objective_var <= 0.0 )
    objective_var = 0.0;

  if( objective_var <= fitness_variance_tolerance )
    return( true );

  return( false );
}

void iamalgam::checkDistributionMultiplierTerminationCondition()
{
  int i;

  for( i = 0; i < number_of_populations; i++ )
  {
    if( !populations_terminated[i] )
      if( distribution_multipliers[i] < 1e-10 )
        populations_terminated[i] = 1;
  }
}

void iamalgam::makeSelections()
{
  int i;
  
  for( i = 0; i < number_of_populations; i++ )
    if( !populations_terminated[i] )
      makeSelectionsForOnePopulation( i );
}

void iamalgam::makeSelectionsForOnePopulation(int population_index)
{
  int i, *sorted;

  sorted = mergeSort( ranks[population_index], population_size );

  if( ranks[population_index][sorted[selection_size-1]] == 0 )
    makeSelectionsForOnePopulationUsingDiversityOnRank0( population_index );
  else
  {
    for( i = 0; i < selection_size; i++ )
    {
        delete selections[i];
        selections[i] = population->solutions[sorted[i]]->clone();
    }
  }

  free( sorted );
}

vec_t<int> iamalgam::countBoundaries(vec_t<solution_mixed *> selections)
{
  // cast selections to solution_BN
  vec_t<int> boundaries_count = vec_t<int>(selections[0]->getNumberOfCVariables(), 0);
  for(int i = 0; i < selections.size(); i++)
  {
    solution_BN *casted_solution = (solution_BN*) selections[i];
    if(casted_solution == NULL)
    {
      // Selections does not contain solution_BNs, so abort counting
      boundaries_count = vec_t<int>();
      break;
    }
    int numBoundaries = casted_solution->getBoundaries()[0].size();
    boundaries_count[numBoundaries]++;
  }

  return boundaries_count;
}

void iamalgam::makeSelectionsForOnePopulationUsingDiversityOnRank0(int population_index)
{
  int     i, j, number_of_rank0_solutions, *preselection_indices,
         *selection_indices, index_of_farthest, number_selected_so_far;
  double *nn_distances, distance_of_farthest, value;

  number_of_rank0_solutions = 0;
  for( i = 0; i < population_size; i++ )
  {
    if( ranks[population_index][i] == 0 )
      number_of_rank0_solutions++;
  }

  preselection_indices = (int *) malloc( number_of_rank0_solutions*sizeof( int ) );
  j                    = 0;
  for( i = 0; i < population_size; i++ )
  {
    if( ranks[population_index][i] == 0 )
    {
      preselection_indices[j] = i;
      j++;
    }
  }

  index_of_farthest    = 0;
  distance_of_farthest = population->solutions[preselection_indices[0]]->getObjectiveValue(); //objective_values[population_index][preselection_indices[0]];
  for( i = 1; i < number_of_rank0_solutions; i++ )
  {
    if( population->solutions[preselection_indices[i]]->getObjectiveValue() > distance_of_farthest )
    {
      index_of_farthest    = i;
      distance_of_farthest = population->solutions[preselection_indices[i]]->getObjectiveValue(); //objective_values[population_index][preselection_indices[i]];
    }
  }

  number_selected_so_far                    = 0;
  selection_indices                         = (int *) malloc( selection_size*sizeof( int ) );
  selection_indices[number_selected_so_far] = preselection_indices[index_of_farthest];
  preselection_indices[index_of_farthest]   = preselection_indices[number_of_rank0_solutions-1];
  number_of_rank0_solutions--;
  number_selected_so_far++;

  nn_distances = (double *) malloc( number_of_rank0_solutions*sizeof( double ) );
  for( i = 0; i < number_of_rank0_solutions; i++ )
  {
      nn_distances[i] = distanceInParameterSpace( population->solutions[preselection_indices[i]]->c_variables, population->solutions[selection_indices[number_selected_so_far-1]]->c_variables );
  }

  while( number_selected_so_far < selection_size )
  {
    index_of_farthest    = 0;
    distance_of_farthest = nn_distances[0];
    for( i = 1; i < number_of_rank0_solutions; i++ )
    {
      if( nn_distances[i] > distance_of_farthest )
      {
        index_of_farthest    = i;
        distance_of_farthest = nn_distances[i];
      }
    }
    
    selection_indices[number_selected_so_far] = preselection_indices[index_of_farthest];
    preselection_indices[index_of_farthest]   = preselection_indices[number_of_rank0_solutions-1];
    nn_distances[index_of_farthest]           = nn_distances[number_of_rank0_solutions-1];
    number_of_rank0_solutions--;
    number_selected_so_far++;

    for( i = 0; i < number_of_rank0_solutions; i++ )
    {
      value = distanceInParameterSpace( population->solutions[preselection_indices[i]]->c_variables, population->solutions[selection_indices[number_selected_so_far-1]]->c_variables ); //distanceInParameterSpace( populations[population_index][preselection_indices[i]], populations[population_index][selection_indices[number_selected_so_far-1]] );
      if( value < nn_distances[i] )
        nn_distances[i] = value;
    }
  }

  for( i = 0; i < selection_size; i++ )
  {
    delete selections[i];
    selections[i] = population->solutions[selection_indices[i]]->clone();
  }

  free( nn_distances );
  free( selection_indices );
  free( preselection_indices );
}

double iamalgam::distanceInParameterSpace(vec_t<double> solution_a, vec_t<double> solution_b)
{
  int    i;
  double value, result;
  
  result = 0.0;
  for( i = 0; i < number_of_parameters; i++ )
  {
    value   = solution_b[i] - solution_a[i];
    result += value*value;
  }
  result = sqrt( result );
  
  return( result );
}

int *iamalgam::mergeSort(double *array, int array_size)
{
  int i, *sorted, *tosort;

  sorted = (int *) malloc( array_size * sizeof( int ) );
  tosort = (int *) malloc( array_size * sizeof( int ) );
  for( i = 0; i < array_size; i++ )
    tosort[i] = i;

  if( array_size == 1 )
    sorted[0] = 0;
  else
    mergeSortWithinBounds( array, sorted, tosort, 0, array_size-1 );

  free( tosort );

  return( sorted );
}

void iamalgam::mergeSortWithinBounds( double *array, int *sorted, int *tosort, int p, int q )
{
  int r;

  if( p < q )
  {
    r = (p + q) / 2;
    mergeSortWithinBounds( array, sorted, tosort, p, r );
    mergeSortWithinBounds( array, sorted, tosort, r+1, q );
    mergeSortMerge( array, sorted, tosort, p, r+1, q );
  }
}

void iamalgam::mergeSortMerge( double *array, int *sorted, int *tosort, int p, int r, int q )
{
  int i, j, k, first;

  i = p;
  j = r;
  for( k = p; k <= q; k++ )
  {
    first = 0;
    if( j <= q )
    {
      if( i < r )
      {
        if( array[tosort[i]] < array[tosort[j]] )
          first = 1;
      }
    }
    else
      first = 1;

    if( first )
    {
      sorted[k] = tosort[i];
      i++;
    }
    else
    {
      sorted[k] = tosort[j];
      j++;
    }
  }

  for( k = p; k <= q; k++ )
    tosort[k] = sorted[k];
}

void iamalgam::makePopulations()
{
  estimateParametersAllPopulations();
  
  copyBestSolutionsToPopulations();

  applyDistributionMultipliers();

  generateAndEvaluateNewSolutionsToFillPopulations();

  computeRanks();

  adaptDistributionMultipliers();
}

void iamalgam::estimateParametersAllPopulations()
{
  int i;

  for( i = 0; i < number_of_populations; i++ )
    if( !populations_terminated[i] )
      estimateParameters( i );
}

void iamalgam::estimateParameters(int population_index )
{
  int i, j;

  estimateParametersML( population_index );

  if( number_of_generations == 0 )
  {
    for( i = 0; i < number_of_parameters; i++ )
    {
      for( j = 0; j < number_of_parameters; j++ )
      {
        aggregated_covariance_matrices[population_index][i][j] = generational_covariance_matrices[population_index][i][j];
        if( i != j )
          aggregated_covariance_matrices[population_index][i][j] = 0.0;
      }
    }
  }
  else
  {
    if( number_of_generations == 1 )
    {
      for( i = 0; i < number_of_parameters; i++ )
        ams_vectors[population_index][i] = mean_vectors[population_index][i]-mean_vectors_previous[population_index][i];
    }
    else
    {
      for( i = 0; i < number_of_parameters; i++ )
        ams_vectors[population_index][i] =
          (1.0-eta_p)*ams_vectors[population_index][i] +
          eta_p*(mean_vectors[population_index][i]-mean_vectors_previous[population_index][i]);
    }

    for( i = 0; i < number_of_parameters; i++ )
    {
      for( j = 0; j < number_of_parameters; j++ )
      {
        aggregated_covariance_matrices[population_index][i][j] =
          (1.0-eta_s)*aggregated_covariance_matrices[population_index][i][j] +
          eta_s*generational_covariance_matrices[population_index][i][j];
      }
    }
  }

  for( i = 0; i < number_of_parameters; i++ )
    for( j = 0; j < number_of_parameters; j++ )
      covariance_matrices[population_index][i][j] = aggregated_covariance_matrices[population_index][i][j];
}

void iamalgam::estimateParametersML(int population_index)
{
  estimateMeanVectorML( population_index );

  estimateCovarianceMatrixML( population_index ); 
}

void iamalgam::estimateMeanVectorML( int population_index )
{
  int i, j;

  if( number_of_generations > 0 )
  {
    for( i = 0; i < number_of_parameters; i++ )
      mean_vectors_previous[population_index][i] = mean_vectors[population_index][i];
  }

  for( i = 0; i < number_of_parameters; i++ )
  {
    mean_vectors[population_index][i] = 0.0;

    for( j = 0; j < selection_size; j++ )
    {
      mean_vectors[population_index][i] += selections[j]->c_variables[i];
    }

    mean_vectors[population_index][i] /= (double) selection_size;
  }

  /* Change the focus of the search to the best solution */
  if( distribution_multipliers[population_index] < 1.0 )
    for( i = 0; i < number_of_parameters; i++ )
    {
      mean_vectors[population_index][i] = selections[0]->c_variables[i];
    }
}

void iamalgam::estimateCovarianceMatrixML( int population_index )
{
  int i, j, k;

  /* First do the maximum-likelihood estimate from data */
  for( i = 0; i < number_of_parameters; i++ )
  {
    for( j = i; j < number_of_parameters; j++ )
    {
      generational_covariance_matrices[population_index][i][j] = 0.0;

      for( k = 0; k < selection_size; k++ )
      {
        generational_covariance_matrices[population_index][i][j] += (selections[k]->c_variables[i]-mean_vectors[population_index][i])*(selections[k]->c_variables[j]-mean_vectors[population_index][j]);
      }
      generational_covariance_matrices[population_index][i][j] /= (double) selection_size;
    }
  }

  for( i = 0; i < number_of_parameters; i++ )
    for( j = 0; j < i; j++ )
      generational_covariance_matrices[population_index][i][j] = generational_covariance_matrices[population_index][j][i];
}

void iamalgam::copyBestSolutionsToPopulations()
{
  int i;

  for( i = 0; i < number_of_populations; i++ )
  {
    if( !populations_terminated[i] )
    {
      delete population->solutions[0];
      population->solutions[0] = selections[0]->clone();
    }

  }
}

void iamalgam::applyDistributionMultipliers()
{
  int i, j, k;
  
  for( i = 0; i < number_of_populations; i++ )
  {
    if( !populations_terminated[i] )
    {
      for( j = 0; j < number_of_parameters; j++ )
        for( k = 0; k < number_of_parameters; k++ )
          covariance_matrices[i][j][k] *= distribution_multipliers[i];
    }
  }
}

void iamalgam::generateAndEvaluateNewSolutionsToFillPopulations()
{
  short   out_of_range;
  int     i, j, k, q, number_of_AMS_solutions;
  double *solution, *solution_AMS, shrink_factor;

  solution_AMS = (double *) malloc( number_of_parameters*sizeof( double ) );

  for( i = 0; i < number_of_populations; i++ )
  {
    computeParametersForSampling( i );

    if( !populations_terminated[i] )
    {
      number_of_AMS_solutions      = (int) (alpha_AMS*(population_size-1));
      samples_drawn_from_normal[i] = 0;
      out_of_bounds_draws[i]       = 0;
      q                            = 0;
    

      for( j = 1; j < population_size; j++ )
      {
        solution = generateNewSolution( i );
  
        for( k = 0; k < number_of_parameters; k++ )
          population->solutions[j]->c_variables[k] = solution[k];

        // If we are dealing with a BN, then normalize the values and update the boundaries after changing c_variables
        if(config->useBN)
        {
          solution_BN *casted_sol = dynamic_cast<solution_BN*>(population->solutions[j]);
          if(config->useNormalizedCVars)
          {
            casted_sol->normalize();
          }
          casted_sol->updateBoundaries();
        }
  
        if( (number_of_generations > 0) && (q < number_of_AMS_solutions) )
        {
          out_of_range  = 1;
          shrink_factor = 2;
          while( (out_of_range == 1) && (shrink_factor > 1e-10) )
          {
            shrink_factor *= 0.5;
            out_of_range   = 0;
            for( k = 0; k < number_of_parameters; k++ )
            {
              solution_AMS[k] = solution[k] + shrink_factor*delta_AMS*distribution_multipliers[i]*ams_vectors[i][k];
              if( !isParameterInRangeBounds( solution_AMS[k], k, (config->extraCVarForNumberOfBins && k < config->data->getNumberOfContinuousVariables())))
              {
                out_of_range = 1;
                break;
              }
            }
          }
          if( !out_of_range )
          {
            for( k = 0; k < number_of_parameters; k++ )
              population->solutions[j]->c_variables[k] = solution_AMS[k];
          }
        }

        q++;
  
        free( solution );
      }

      vec_t<solution_t<int> *> casted_pop;
      for(solution_mixed *sol : population->solutions)
      {
        casted_pop.push_back(dynamic_cast<solution_t<int>*>(sol));
      }
      problemInstance->evaluatePopulation(casted_pop, true);
    }
  }

  free( solution_AMS );
}

double iamalgam::averageFitnessPopulation()
{
  double average = 0;
  for(int i = 0; i < population_size; i++)
  {
    average += population->solutions[i]->getObjectiveValue();
  }
  return average/population_size;
}

void iamalgam::computeParametersForSampling(int population_index)
{
  int i;

  if( cholesky_factors_lower_triangle[population_index] )
  {
    for( i = 0; i < number_of_parameters; i++ )
      free( cholesky_factors_lower_triangle[population_index][i] );
    free( cholesky_factors_lower_triangle[population_index] );
  }

  cholesky_factors_lower_triangle[population_index] = choleskyDecomposition( covariance_matrices[population_index], number_of_parameters );
}

double **iamalgam::choleskyDecomposition( double **matrix, int n )
{
  int     i, j, k, info, *ipvt;
  double *a, *work, **result;
  
  a    = (double *) malloc( n*n*sizeof( double ) );
  work = (double *) malloc( n*sizeof( double ) );
  ipvt = (int *) malloc( n*sizeof( int ) );

  k = 0;
  for( i = 0; i < n; i++ )
  {
    for( j = 0; j < n; j++ )
    {
      a[k] = matrix[i][j];
      k++;
    }
    ipvt[i] = 0;
  }

  info = linpackDCHDC( a, n, n, work, ipvt );

  result = matrixNew( n, n );
  if( info != n ) /* Matrix is not positive definite */
  {
    k = 0;
    for( i = 0; i < n; i++ )
    {
      for( j = 0; j < n; j++ )
      {
        result[i][j] = i != j ? 0.0 : sqrt( matrix[i][j] );
        k++;
      }
    }
  }
  else
  {
    k = 0;
    for( i = 0; i < n; i++ )
    {
      for( j = 0; j < n; j++ )
      {
        result[i][j] = i < j ? 0.0 : a[k];
        k++;
      }
    }
  }

  free( ipvt );
  free( work );
  free( a );
  
  return( result );
}

int iamalgam::linpackDCHDC( double a[], int lda, int p, double work[], int ipvt[] )
{
  int    info, j, jp, k, l, maxl, pl, pu;
  double maxdia, temp;

  pl   = 1;
  pu   = 0;
  info = p;
  for( k = 1; k <= p; k++ )
  {
    maxdia = a[k-1+(k-1)*lda];
    maxl   = k;
    if( pl <= k && k < pu )
    {
      for( l = k+1; l <= pu; l++ )
      {
        if( maxdia < a[l-1+(l-1)*lda] )
        {
          maxdia = a[l-1+(l-1)*lda];
          maxl   = l;
        }
      }
    }

    if( maxdia <= 0.0 )
    {
      info = k - 1;

      return( info );
    }

    if( k != maxl )
    {
      blasDSWAP( k-1, a+0+(k-1)*lda, 1, a+0+(maxl-1)*lda, 1 );

      a[maxl-1+(maxl-1)*lda] = a[k-1+(k-1)*lda];
      a[k-1+(k-1)*lda]       = maxdia;
      jp                     = ipvt[maxl-1];
      ipvt[maxl-1]           = ipvt[k-1];
      ipvt[k-1]              = jp;
    }
    work[k-1]        = sqrt( a[k-1+(k-1)*lda] );
    a[k-1+(k-1)*lda] = work[k-1];

    for( j = k+1; j <= p; j++ )
    {
      if( k != maxl )
      {
        if( j < maxl )
        {
          temp                = a[k-1+(j-1)*lda];
          a[k-1+(j-1)*lda]    = a[j-1+(maxl-1)*lda];
          a[j-1+(maxl-1)*lda] = temp;
        }
        else if ( maxl < j )
        {
          temp                = a[k-1+(j-1)*lda];
          a[k-1+(j-1)*lda]    = a[maxl-1+(j-1)*lda];
          a[maxl-1+(j-1)*lda] = temp;
        }
      }
      a[k-1+(j-1)*lda] = a[k-1+(j-1)*lda] / work[k-1];
      work[j-1]        = a[k-1+(j-1)*lda];
      temp             = -a[k-1+(j-1)*lda];

      blasDAXPY( j-k, temp, work+k, 1, a+k+(j-1)*lda, 1 );
    }
  }

  return( info );
}

int iamalgam::blasDSWAP( int n, double *dx, int incx, double *dy, int incy )
{
  double dtmp;
  
  if (n > 0)
  {
    incx *= sizeof( double );
    incy *= sizeof( double );

    dtmp  = (*dx);
    *dx   = (*dy);
    *dy   = dtmp;

    while( (--n) > 0 )
    {
      dx = (double *) ((char *) dx + incx);
      dy = (double *) ((char *) dy + incy);
      dtmp = (*dx); *dx = (*dy); *dy = dtmp;
    }
  }

  return( 0 );
}

int iamalgam::blasDAXPY(int n, double da, double *dx, int incx, double *dy, int incy)
{
  double dtmp0, dtmp, *dx0, *dy0;
  
  if( n > 0 && da != 0. )
  {
    incx *= sizeof(double);
    incy *= sizeof(double);
    *dy  += da * (*dx);

    if( (n & 1) == 0 )
    {
      dx   = (double *) ((char *) dx + incx);
      dy   = (double *) ((char *) dy + incy);
      *dy += da * (*dx);
      --n;
    }
    n = n >> 1;
    while( n > 0 )
    {
      dy0   = (double *) ((char *) dy + incy);
      dy    = (double *) ((char *) dy0 + incy);
      dtmp0 = (*dy0);
      dtmp  = (*dy); 
      dx0   = (double *) ((char *) dx + incx); 
      dx    = (double *) ((char *) dx0 + incx);
      *dy0  = dtmp0 + da * (*dx0);
      *dy   = dtmp + da * (*dx);
      --n;
    }
  }

  return( 0 );
}

double **iamalgam::matrixNew( int n, int m )
{
  int      i;
  double **result;

  result = (double **) malloc( n*( sizeof( double * ) ) );
  for( i = 0; i < n; i++ )
    result[i] = (double *) malloc( m*( sizeof( double ) ) );

  return( result );
}

double *iamalgam::generateNewSolution( int population_index )
{
  short   ready;
  int     i, times_not_in_bounds;
  double *result, *z;

  times_not_in_bounds = -1;
  out_of_bounds_draws[population_index]--;

  ready = 0;
  do
  {
    times_not_in_bounds++;
    samples_drawn_from_normal[population_index]++;
    out_of_bounds_draws[population_index]++;
    if( times_not_in_bounds >= 100 )
    {
      result = (double *) malloc( number_of_parameters*sizeof( double ) );
      for( i = 0; i < number_of_parameters; i++ )
      {
        if(config->extraCVarForNumberOfBins && i < config->data->getNumberOfContinuousVariables())
        {
          result[i] =  fmin( std::lround( ((gomea::utils::rng)() / (double)(gomea::utils::rng).max()) * (config->maxDiscretizations - 1) + 1.5), config->maxDiscretizations);
        }
        else if(config->transformCVariables)
        {
          result[i] = 1.0 / (config->lower_user_range + ((gomea::utils::rng)() / (double)(gomea::utils::rng).max()) * (config->upper_user_range - config->lower_user_range));
        }
        else
        {
          result[i] = config->lower_user_range + ((gomea::utils::rng)() / (double)(gomea::utils::rng).max()) * (config->upper_user_range - config->lower_user_range);
        }
      }
    }
    else
    {
      z = (double *) malloc( number_of_parameters*sizeof( double ) );

      for( i = 0; i < number_of_parameters; i++ )
         z[i] = random1DNormalUnit();

      result = matrixVectorMultiplication( cholesky_factors_lower_triangle[population_index], z, number_of_parameters, number_of_parameters );

      for( i = 0; i < number_of_parameters; i++ )
        result[i] += mean_vectors[population_index][i];

      free( z );
    }
    
    ready = 1;
    for( i = 0; i < number_of_parameters; i++ )
    {
      if( !isParameterInRangeBounds( result[i], i, (config->extraCVarForNumberOfBins && i < config->data->getNumberOfContinuousVariables()) ) )
      {
        ready = 0;
        break;
      }
    }
    if( !ready )
      free( result );
  }
  while( !ready );

  return( result );
}

double iamalgam::randomRealUniform01()
{
  int64_t n26, n27;
  double  result;

  random_seed_changing = (random_seed_changing * 0x5DEECE66DLLU + 0xBLLU) & ((1LLU << 48) - 1);
  n26                  = (int64_t)(random_seed_changing >> (48 - 26));
  random_seed_changing = (random_seed_changing * 0x5DEECE66DLLU + 0xBLLU) & ((1LLU << 48) - 1);
  n27                  = (int64_t)(random_seed_changing >> (48 - 27));
  result               = (((int64_t)n26 << 27) + n27) / ((double) (1LLU << 53));

  return( result );
}

double iamalgam::random1DNormalUnit()
{
  double v1, v2, s, multiplier, value;
  s = 0.0;

  if( haveNextNextGaussian )
  {
    haveNextNextGaussian = false;

    return( nextNextGaussian );
  }
  else
  {
    do
    {
      v1 = 2 * (randomRealUniform01()) - 1;
      v2 = 2 * (randomRealUniform01()) - 1;
      s = v1 * v1 + v2 * v2;
    } while (s >= 1);

    value                = -2 * log(s)/s;
    multiplier           = value <= 0.0 ? 0.0 : sqrt( value );
    nextNextGaussian     = v2 * multiplier;
    haveNextNextGaussian = true;

    return( v1 * multiplier );
  }
}

double *iamalgam::matrixVectorMultiplication( double **matrix, double *vector, int n0, int n1 )
{
  int     i;
  double *result;
  
  result = (double *) malloc( n0*sizeof( double ) );
  for( i = 0; i < n0; i++ )
    result[i] = vectorDotProduct( matrix[i], vector, n1 );
  
  return( result );
}

double iamalgam::vectorDotProduct( double *vector0, double *vector1, int n0 )
{
  int    i;
  double result;
  
  result = 0.0;
  for( i = 0; i < n0; i++ )
    result += vector0[i]*vector1[i];
    
  return( result );
}

bool iamalgam::isParameterInRangeBounds( double parameter, int dimension, bool is_extra_cvar)
{
  if(is_extra_cvar)
  {
    long rounded_parameter = std::lround(parameter);
    if( rounded_parameter < 1 || rounded_parameter > config->maxDiscretizations )
    {
      return( false );
    }
    return( true );
  }
  else
  {
    if( parameter < problemInstance->getLowerRangeBound(dimension) ||
        parameter > problemInstance->getUpperRangeBound(dimension) ||
        isnan( parameter ) )
    {
      return( false );
    }
    
    return( true );
  }
}

void iamalgam::adaptDistributionMultipliers()
{
  int    i;

  for( i = 0; i < number_of_populations; i++ )
  {
    adaptDistributionMultipliersForOnePopulation( i );
  }
}

void iamalgam::adaptDistributionMultipliersForOnePopulation(int i)
{
  short  improvement;
  double st_dev_ratio;

  if( !populations_terminated[i] )
    {
      if( (((double) out_of_bounds_draws[i])/((double) samples_drawn_from_normal[i])) > 0.9 )
        distribution_multipliers[i] *= 0.5;
  
      improvement = generationalImprovementForOnePopulation( i, &st_dev_ratio );
    
      if( improvement )
      {
        no_improvement_stretch[i] = 0;

        if( distribution_multipliers[i] < 1.0 )
          distribution_multipliers[i] = 1.0;
  
        if( st_dev_ratio > st_dev_ratio_threshold )
          distribution_multipliers[i] *= distribution_multiplier_increase;
      }
      else
      {
        if( distribution_multipliers[i] <= 1.0 )
          (no_improvement_stretch[i])++;
  
        if( (distribution_multipliers[i] > 1.0) || (no_improvement_stretch[i] >= maximum_no_improvement_stretch) )
          distribution_multipliers[i] *= distribution_multiplier_decrease;
  
        if( (no_improvement_stretch[i] < maximum_no_improvement_stretch) && (distribution_multipliers[i] < 1.0) )
          distribution_multipliers[i] = 1.0;
      }
    }
}

bool iamalgam::generationalImprovementForOnePopulation( int population_index, double *st_dev_ratio )
{
  int     i, j, index_best_selected, index_best_population,
          number_of_improvements;
  double *average_parameters_of_improvements;

  /* Determine best selected solutions */
  index_best_selected = 0;
  for( i = 0; i < selection_size; i++ )
  {
    if( betterFitness(selections[i]->getObjectiveValue(), selections[i]->getConstraintValue(),
                      selections[index_best_selected]->getObjectiveValue(), selections[index_best_selected]->getConstraintValue()) )
      index_best_selected = i;
  }

  /* Determine best in the population and the average improvement parameters */
  average_parameters_of_improvements = (double *) malloc( number_of_parameters*sizeof( double ) );
  for( i = 0; i < number_of_parameters; i++ )
    average_parameters_of_improvements[i] = 0.0;

  index_best_population   = 0;
  number_of_improvements  = 0;
  for( i = 0; i < population_size; i++ )
  {
    if( betterFitness( population->solutions[i]->getObjectiveValue(), population->solutions[i]->getConstraintValue(),
                       population->solutions[index_best_population]->getObjectiveValue(), population->solutions[index_best_population]->getConstraintValue()) )
      index_best_population = i;

    if( betterFitness( population->solutions[i]->getObjectiveValue(), population->solutions[i]->getConstraintValue(),
                       selections[index_best_selected]->getObjectiveValue(), selections[index_best_selected]->getConstraintValue()) )
    {
      number_of_improvements++;
      for( j = 0; j < number_of_parameters; j++ )
      {
        average_parameters_of_improvements[j] += population->solutions[i]->c_variables[j];
      }
    }
  }

  /* Determine st.dev. ratio */
  *st_dev_ratio = 0.0;
  if( number_of_improvements > 0 )
  {
    for( i = 0; i < number_of_parameters; i++ )
      average_parameters_of_improvements[i] /= (double) number_of_improvements;

    *st_dev_ratio = getStDevRatio( population_index, average_parameters_of_improvements );
  }

  free( average_parameters_of_improvements );


  if( fabs( selections[index_best_selected]->getObjectiveValue() - population->solutions[index_best_population]->getObjectiveValue() ) == 0.0 )
    return( false );

  return( true );
}

double iamalgam::getStDevRatio(int population_index, double *parameters)
{
  int      i;
  double **inverse, result, *x_min_mu, *z;

  inverse = matrixLowerTriangularInverse( cholesky_factors_lower_triangle[population_index], number_of_parameters );

  x_min_mu = (double *) malloc( number_of_parameters*sizeof( double ) );
  
  for( i = 0; i < number_of_parameters; i++ )
    x_min_mu[i] = parameters[i]-mean_vectors[population_index][i];

  z = matrixVectorMultiplication( inverse, x_min_mu, number_of_parameters, number_of_parameters );

  result = 0.0;
  for( i = 0; i < number_of_parameters; i++ )
  {
    if( fabs( z[i] ) > result )
      result = fabs( z[i] );
  }

  free( z );
  free( x_min_mu );
  for( i = 0; i < number_of_parameters; i++ )
    free( inverse[i] );
  free( inverse );

  return( result );
}

double **iamalgam::matrixLowerTriangularInverse( double **matrix, int n )
{
  int     i, j, k;
  double *t, **result;
  
  t = (double *) malloc( n*n*sizeof( double ) );

  k = 0;
  for( i = 0; i < n; i++ )
  {
    for( j = 0; j < n; j++ )
    {
      t[k] = matrix[j][i];
      k++;
    }
  }

  linpackDTRDI( t, n, n );

  result = matrixNew( n, n );
  k = 0;
  for( i = 0; i < n; i++ )
  {
    for( j = 0; j < n; j++ )
    {
      result[j][i] = i > j ? 0.0 : t[k];
      k++;
    }
  }

  free( t );
  
  return( result );
}

int iamalgam::linpackDTRDI( double t[], int ldt, int n )
{
  int    j, k, info;
  double temp;

  info = 0;
  for( k = n; 1 <= k; k-- )
  {
    if ( t[k-1+(k-1)*ldt] == 0.0 )
    {
      info = k;
      break;
    }

    t[k-1+(k-1)*ldt] = 1.0 / t[k-1+(k-1)*ldt];
    temp = -t[k-1+(k-1)*ldt];

    if ( k != n )
    {
      blasDSCAL( n-k, temp, t+k+(k-1)*ldt, 1 );
    }

    for( j = 1; j <= k-1; j++ )
    {
      temp = t[k-1+(j-1)*ldt];
      t[k-1+(j-1)*ldt] = 0.0;
      blasDAXPY( n-k+1, temp, t+k-1+(k-1)*ldt, 1, t+k-1+(j-1)*ldt, 1 );
    }
  }

  return( info );
}

void iamalgam::blasDSCAL( int n, double sa, double x[], int incx )
{
  int i, ix, m;

  if( n <= 0 )
  {
  }
  else if( incx == 1 )
  {
    m = n % 5;

    for( i = 0; i < m; i++ )
    {
      x[i] = sa * x[i];
    }

    for( i = m; i < n; i = i + 5 )
    {
      x[i]   = sa * x[i];
      x[i+1] = sa * x[i+1];
      x[i+2] = sa * x[i+2];
      x[i+3] = sa * x[i+3];
      x[i+4] = sa * x[i+4];
    }
  }
  else
  {
    if( 0 <= incx )
    {
      ix = 0;
    }
    else
    {
      ix = ( - n + 1 ) * incx;
    }

    for( i = 0; i < n; i++ )
    {
      x[ix] = sa * x[ix];
      ix = ix + incx;
    }
  }
}

void iamalgam::determineBestSolutionSoFar()
{
  int i, population_of_best, index_of_best;

  determineBestSolutionInCurrentPopulations( &population_of_best, &index_of_best );

  if( number_of_starts == 1 ||
      betterFitness( population->solutions[index_of_best]->getObjectiveValue(),
                     population->solutions[index_of_best]->getConstraintValue(),
                     best_so_far_objective_value,
                     best_so_far_constraint_value ) )                     
  {
    best_so_far_objective_value  = population->solutions[index_of_best]->getObjectiveValue();
    best_so_far_constraint_value = population->solutions[index_of_best]->getConstraintValue();
    for( i = 0; i < number_of_parameters; i++ )
      best_so_far_solution[i] = population->solutions[index_of_best]->c_variables[i];
  }
}

void iamalgam::determineBestSolutionInCurrentPopulations(int *population_of_best, int *index_of_best)
{
  int i, j;

  (*population_of_best) = 0;
  (*index_of_best)      = 0;
  for( i = 0; i < number_of_populations; i++ )
  {
    for( j = 0; j < population_size; j++ )
    {
      if( betterFitness( population->solutions[j]->getObjectiveValue(), population->solutions[j]->getConstraintValue(),
                         population->solutions[(*index_of_best)]->getObjectiveValue(), population->solutions[(*index_of_best)]->getConstraintValue() ) )
      {
        (*population_of_best) = i;
        (*index_of_best)      = j;
      }
    }
  }
}

void iamalgam::ezilaitini()
{
  ezilaitiniMemory();

  ezilaitiniDistributionMultipliers();
}

void iamalgam::ezilaitiniMemory()
{
  int i, j;

  for( i = 0; i < number_of_populations; i++ )
  {
    free( ranks[i] );

    free( mean_vectors[i] );

    free( mean_vectors_previous[i] );

    for( j = 0; j < number_of_parameters; j++ )
      free( covariance_matrices[i][j] );
    free( covariance_matrices[i] );

    for( j = 0; j < number_of_parameters; j++ )
      free( generational_covariance_matrices[i][j] );
    free( generational_covariance_matrices[i] );
  
    for( j = 0; j < number_of_parameters; j++ )
      free( aggregated_covariance_matrices[i][j] );
    free( aggregated_covariance_matrices[i] );

    if( cholesky_factors_lower_triangle[i] )
    {
      for( j = 0; j < number_of_parameters; j++ )
        free( cholesky_factors_lower_triangle[i][j] );
      free( cholesky_factors_lower_triangle[i] );
    }

    free( ams_vectors[i] );
  }

  free( covariance_matrices );
  free( generational_covariance_matrices );
  free( aggregated_covariance_matrices );
  free( cholesky_factors_lower_triangle );
  free( populations_terminated );
  free( no_improvement_stretch );
  free( ranks );
  free( mean_vectors );
  free( mean_vectors_previous );
  free( ams_vectors );

  // delete selections
  for(solution_mixed *sol : selections)
  {
    delete sol;
  }
}

void iamalgam::ezilaitiniDistributionMultipliers()
{
  free( distribution_multipliers );
  free( samples_drawn_from_normal );
  free( out_of_bounds_draws );
}


// Note: original run function of iAMaLGaM as reference. Not used in MIHEA, but gives nice overview of normal iAMaLGaM flow
void iamalgam::runOnce() 
{
  number_of_starts++;

  initialize();

  while( !checkTerminationConditionForRunOnce())
  {
    makeSelections();

    makePopulations();

    number_of_generations++;
  }

  determineBestSolutionSoFar();

  ezilaitini();


}

// IAMaLGaM learning continuous model step from MIHEA
void iamalgam::learnContinuousModel(int population_index)
{
  // First, make sure there is only 1 population by checking the number_of_populations and population_index
  assert(number_of_populations == 1 && population_index == 0);
  writePopulationBoundaryStatsToFile(config->folder, population->solutions, "GEN " + to_string(number_of_generations) + " - POPULATION");

  computeRanksForOnePopulation(population_index);

  makeSelectionsForOnePopulation(population_index);
  writePopulationToFile(config->folder, selections, "SELECTIONS in iamalgam ----------------------------------", config->logDebugInformation);  
  writePopulationBoundaryStatsToFile(config->folder, selections, "GEN " + to_string(number_of_generations) + " - SELECTIONS");

  estimateParameters(population_index);
  writeMatrixToFile(config->folder, covariance_matrices[population_index], number_of_parameters, number_of_parameters,  "COVARIANCE MATRIX in iamalgam ----------------------------------", config->logDebugInformation);
  writeMatrixToFile(config->folder, aggregated_covariance_matrices[population_index], number_of_parameters, number_of_parameters,  "AGGREGATED COVARIANCE MATRIX in iamalgam ----------------------------------", config->logDebugInformation);
  writeMatrixToFile(config->folder, generational_covariance_matrices[population_index], number_of_parameters, number_of_parameters,  "GENERATIONAL COVARIANCE MATRIX in iamalgam ----------------------------------", config->logDebugInformation);
  writeVectorToFile(config->folder, ams_vectors[population_index], number_of_parameters, "AMS VECTOR in iamalgam ----------------------------------", config->logDebugInformation);
  writeMessageToLogFile(config->folder, "distribution multiplier: " + to_string(distribution_multipliers[population_index]) + "\n", config->logDebugInformation);

  copyBestSolutionsToPopulations(); 

  applyDistributionMultipliers();
}

void iamalgam::generateNewPopulation(int population_index)
{
  short   out_of_range;
  int     j, k, q, number_of_AMS_solutions;
  double *solution, *solution_AMS, shrink_factor;

  solution_AMS = (double *) malloc( number_of_parameters*sizeof( double ) );
  
  computeParametersForSampling(population_index);

  if( !populations_terminated[population_index] )
  {
    number_of_AMS_solutions      = (int) (alpha_AMS*(population_size-1));
    samples_drawn_from_normal[population_index] = 0;
    out_of_bounds_draws[population_index]       = 0;
    q                            = 0;
    for( j = 1; j < population_size; j++ )
    {
      solution = generateNewSolution( population_index );

      for( k = 0; k < number_of_parameters; k++ )
      {
        population->solutions[j]->c_variables[k] = solution[k];
      }

      if( (number_of_generations > 0) && (q < number_of_AMS_solutions) )
      {
        out_of_range  = 1;
        shrink_factor = 2;
        while( (out_of_range == 1) && (shrink_factor > 1e-10) )
        {
          shrink_factor *= 0.5;
          out_of_range   = 0;
          for( k = 0; k < number_of_parameters; k++ )
          {
            solution_AMS[k] = solution[k] + shrink_factor*delta_AMS*distribution_multipliers[population_index]*ams_vectors[population_index][k];
            if( !isParameterInRangeBounds( solution_AMS[k], k ) )
            {
              out_of_range = 1;
              break;
            }
          }
        }
        if( !out_of_range )
        {
          for( k = 0; k < number_of_parameters; k++ )
          {
            population->solutions[j]->c_variables[k] = solution_AMS[k];
          }
        }
      }

      q++;

      free( solution );
    }
  }

  free( solution_AMS );
}


}}