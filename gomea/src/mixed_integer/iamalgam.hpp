#pragma once


#include "gomea/src/mixed_integer/Config.hpp"
#include "gomea/src/mixed_integer/utils.hpp"
#include "gomea/src/mixed_integer/Solutionset.hpp"

namespace gomea{
namespace mixedinteger{

class iamalgam
{
    public:
        Config *config;
        int basePopulationSize;
        fitness_t *problemInstance = NULL;
        

        iamalgam();
        iamalgam(Config *config_);
        iamalgam(Config *config_, Solutionset *population_);
        ~iamalgam();
        
        void ezilaitini();
        void ezilaitiniMemory();
        void ezilaitiniDistributionMultipliers();
        void run();
        
        void runOnce();
        void initialize();
        void mergeSortFitnessWithinBounds( double *objectives, double *constraints, int *sorted, int *tosort, int p, int q );
        void mergeSortFitnessMerge( double *objectives, double *constraints, int *sorted, int *tosort, int p, int r, int q );
        int *mergeSort(double *array, int array_size);
        void mergeSortWithinBounds( double *array, int *sorted, int *tosort, int p, int q );
        void mergeSortMerge( double *array, int *sorted, int *tosort, int p, int r, int q );
        void printVerboseOverview();
        bool checkTerminationConditionForRunOnce(); 
        bool checkNumberOfEvaluationsTerminationCondition();
        bool checkVTRTerminationCondition();
        void determineBestSolutionInCurrentPopulations(int *population_of_best, int *index_of_best); // These pointers are not arrays, just single values
        void checkFitnessVarianceTermination();
        bool checkFitnessVarianceTerminationSinglePopulation(int population_index);
        void checkDistributionMultiplierTerminationCondition();
        void writeGenerationalStatistics();
        bool betterFitness(double objective_value_x, double constraint_value_x, double objective_value_y, double constraint_value_y);
        void writeGenerationalSolutions(bool final);
        void makeSelections();
        void makeSelectionsForOnePopulation(int population_index);
        void makeSelectionsForOnePopulationUsingDiversityOnRank0(int population_index);
        double distanceInParameterSpace(vec_t<double> solution_a, vec_t<double> solution_b);
        void makePopulations();
        void estimateParametersAllPopulations();
        void estimateParameters(int population_index);
        void estimateParametersML(int population_index);
        void estimateMeanVectorML(int population_index);
        void estimateCovarianceMatrixML(int population_index);
        void copyBestSolutionsToPopulations();
        void applyDistributionMultipliers();
        void generateAndEvaluateNewSolutionsToFillPopulations();
        double averageFitnessPopulation();
        void computeParametersForSampling(int population_index);
        double **choleskyDecomposition( double **matrix, int n );
        int linpackDCHDC( double a[], int lda, int p, double work[], int ipvt[] );
        int linpackDTRDI( double t[], int ldt, int n );
        void blasDSCAL( int n, double sa, double x[], int incx );
        int blasDSWAP( int n, double *dx, int incx, double *dy, int incy );
        int blasDAXPY(int n, double da, double *dx, int incx, double *dy, int incy);
        double **matrixNew( int n, int m );
        double *generateNewSolution(int population_index);
        double randomRealUniform01();
        double random1DNormalUnit();
        double *matrixVectorMultiplication( double **matrix, double *vector, int n0, int n1 );
        double vectorDotProduct( double *vector0, double *vector1, int n0 );
        bool isParameterInRangeBounds(double parameter, int dimension, bool is_extra_cvar = false);
        void computeRanks();
        void adaptDistributionMultipliers();
        void adaptDistributionMultipliersForOnePopulation(int i = 0);
        bool generationalImprovementForOnePopulation(int population_index, double *st_dev_ratio); // st_dev_ratio is a pointer to a single value
        double getStDevRatio(int population_index, double *parameters);
        double **matrixLowerTriangularInverse( double **matrix, int n );
        void determineBestSolutionSoFar();

        bool checkTerminationCondition();
        void initializeRandomNumberGenerator();
        void initializeMemory();
        void initializeDistributionMultipliers();
        void initializePopulationsAndFitnessValues();
        void computeRanksForOnePopulation(int population_index);
        int *mergeSortFitness(double *objectives, double *constraints, int number_of_solutions);

        void learnContinuousModel(int population_index = 0);
        void generateNewPopulation(int population_index = 0);

        vec_t<int> countBoundaries(vec_t<solution_mixed *> selections);

        bool write_generational_statistics;
        bool write_generational_solutions;
        int number_of_generations;

        int number_of_starts = 0;
        int number_of_evaluations;
        int maximum_number_of_evaluations;
        double fitness_variance_tolerance;
        int number_of_parameters;
        int number_of_populations = 1; // Should always be 1
        int population_size;
        int selection_size;
        bool use_vtr;
        double tau = 0.35;
        double eta_p;
        double eta_s;
        double alpha_AMS;
        double delta_AMS;
        Solutionset *population; // NOTE: This population is the *offspringPopulation* from Population.cpp! (if -k is not active)
        double *distribution_multipliers;
        double distribution_multiplier_increase;
        double distribution_multiplier_decrease = 0.9;
        int *samples_drawn_from_normal;
        int *out_of_bounds_draws;
        bool *populations_terminated;
        int *no_improvement_stretch;
        double **ranks;
        vec_t<solution_mixed*> selections;
        double **mean_vectors;
        double **mean_vectors_previous;
        double ***covariance_matrices;
        double ***generational_covariance_matrices;
        double ***aggregated_covariance_matrices;
        double ***cholesky_factors_lower_triangle;
        double **ams_vectors;
        double *lower_range_bounds;
        double *upper_range_bounds;
        double *lower_init_ranges;
        double *upper_init_ranges;
        double *best_so_far_solution;
        double best_so_far_objective_value;
        double best_so_far_constraint_value;
        double st_dev_ratio_threshold = 1.0;
        int maximum_no_improvement_stretch = 200;
        bool haveNextNextGaussian;
        double nextNextGaussian;
        int64_t random_seed_changing;

};

}}
