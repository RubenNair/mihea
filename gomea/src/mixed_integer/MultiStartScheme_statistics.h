// //
// // Created by Damy Ha on 27-Sep-22.
// //

// #include <chrono>
// #include <iostream>
// #include <limits>
// #include <memory>
// #include <time.h>


// #include "gomea/src/fitness/fitness_BN.h"

// using namespace std;

// class MultiStartScheme {
// public:
//     MultiStartScheme(size_t algorithm_index,
//                      size_t discretizationPolicy,
//                      const shared_ptr<Fitness> fitnessFunction,
//                      size_t maximum_number_of_evaluations,
//                      long maximum_time_in_seconds,
//                      bool use_vtr,
//                      double vtr,
//                      shared_ptr<mt19937> rng,
//                      bool write_generational_statistics,
//                      bool print_progress,
//                      int runIndex,
//                      string writeDirectory);         // Constructor
//     ~MultiStartScheme();        // Destructor

//     void run();                 // Run the process

//     void executeOptimizers();   // Execute the optimizers

//     // Getters
//     clock_t getStartingTime() const;
//     clock_t getStartingRunTime() const;
//     clock_t getEndTime() const;
//     size_t getCurrentNumberOfInstances();
//     const Archive_so_solution &getElitistSolutions() const;
//     int getReferenceSolutionIndex() const;
//     int getDoPostRunDiscretization() const;

//     // Setter
//     void setStartingTime(clock_t startingTime);                 // Set the starting time (since turning on the algorithm)
//     void setPrintFos(bool printFos);                            // Print GOMEA optimizer's FOS
//     void setLinkageModelIndex(int linkageModelIndex);           // The linkage model of GOMEA
//     void setLocalSearchIndex(int local_search_index);           // The local optimizer index of BN-GOMEA
//     void setReferenceSolutionIndex(int referenceSolutionIndex); // Sets the reference solution index
//     void setDoPostRunDiscretization(int doPostRunDiscretization);


// private:
//     // Scheme
//     size_t algorithm_index;             // The algorithm to use to generate new instances
//     int discretization_policy_index;    // The index of the discretization policy to use
//     bool use_vtr;                       // Use the VTR as a termination threshold
//     double vtr;                         // The fitness value to reach
//     bool terminateImmediately;          // Terminates the program after 1 generation
//     int runIndex;                       // The index of the run

//     // Multi-pop scheme
//     size_t maximum_number_of_instances;                     // The maximum number of instances to run
//     size_t maximum_number_of_evaluations;                   // The maximum number of evaluations to execute
//     long maximum_time_in_seconds;                           // The maximum execution time in seconds
//     size_t number_of_subgenerations_per_instance_factor;    // The difference in generations (factor) between population sizes
//     size_t base_population_size;                            // The starting population size

//     vector<shared_ptr<Optimizer_SO>> optimizers;            // The optimizers solving the same problem
//     size_t number_of_generations;                           // The number of generations of the multi-pop scheme
//     size_t indexMinimumInstance;                            // The minimum optimizer index that corresponds to the GOMEA that is still allowed to run (lower ones should be stopped because of average fitness being lower than that of a higher one).

//     // Fitness function
//     shared_ptr<Fitness> fitnessFunction;            // The fitness function passed by main

//     // Statistics
//     bool print_progress;                                    // Print the progress of the multi start scheme
//     Archive_so_solution archive;                            // The best found solutions
//     bool write_generational_statistics;                     // Write the generational statistics file
//     string writeDirectory;                                  // The directory to write statistics into
//     string pathMultiStartSchemeStatistics;                  // Path to write the multi start scheme statistics to
//     string pathMultiStartSchemeSolutions;                   // Path to write the multi start scheme solutions to
//     shared_ptr<ofstream> streamMultiStartSchemeStatistics;  // Stream to write multi start scheme the statistics to
//     shared_ptr<ofstream> streamMultiStartSchemeSolutions;   // Stream to write multi start scheme solutions to


//     // RNG
//     shared_ptr<mt19937> rng;    // The random number generator

//     // Time
//     clock_t startingTime;       // The starting time
//     clock_t startingRunTime;    // The time when starting the optimization
//     clock_t endTime;            // The end time

//     // BN-GOMEA
//     bool print_FOS;                 // Prints the FOS (Only useful for FOS based algorithms
//     int linkage_model_index;        // The linkage model index of GOMEA
//     int localSearchIndex;           // The local optimizer index of BN-GOMEA

//     // Reference Optimizer
//     int referenceSolutionIndex;     // The reference solution to plot

//     // Other settings
//     int doPostRunDiscretization;    // Performs a discretization

//     // Multi start scheme
//     void initializeNewInstance();                       // Initializes a new optimizer instances
//     void generationalStepAllOptimizers();               // Starts the recursive optimizer execution scheme
//     void generationalStepAllGOMEAsRecursiveFold(size_t indexSmallestInstance, size_t indexLargestInstance);      // The recursive execution method
//     size_t findSmallestActiveInstance();                // Finds the optimizer with the smallest population that is still active
//     bool instanceShouldTerminate(size_t instanceIndex); // Determines if an optimizer instance should terminate because of its performance (with respect to other instances)
//     void determineElitistSolutions();                   // Retrieves the elitist solution from all optimizers and keeps the best solution
//     void executePostRunDiscretization();                // Executes a post run discretization

//     // Optimizer initialization
//     size_t determineNextTargetPopulationSize();                         // Determines the target population size for the next optimizer
//     void initializeReferenceSolutionOptimizer();                        // Initializes a Reference solution optimizer
//     void initializeCBN_GOMEAOptimizer(size_t targetPopulationSize);     // Initializes a CBN-GOMEA optimizer
//     void initializePostRunDiscretizationOptimizer(const vector<shared_ptr<Solution_so>>& solutionsToRediscretize, size_t targetPopulationSize, size_t discretizationIndex, bool overWriteComputationTime);   // Initializes a post run discretization optimizer

//     // Termination
//     bool checkTermination();                            // Checks the termination criterion
//     bool checkTerminateImmediately();                   // Checks if the algorithm should terminate after the first generation
//     bool checkEvaluationBudgetHasBeenExceeded();        // Checks if the valuation budget has been exceeded
//     bool checkTimeBudgetHasBeenExceeded();              // Checks if the time budget has been exceeded
//     bool fitnessValueHasBeenReached();                  // Checks if the vtr has been reached

//     // Time
//     double getRunTime();

//     // Statistics
//     void closeStatisticsStreamsOfAllOptmizers();                        // Closes the statistics file for all optimizers
//     void initialize_multi_start_scheme_statistics_file();               // Initializes the multi start scheme statistics file
//     void initialize_multi_start_scheme_solutions_file();                // Initializes the multi start scheme solutions file
//     shared_ptr<Solution_so> determineSolutionToWrite();                 // Determines the solution to write statistics of from the archive
//     void write_multi_start_scheme_statistics(size_t indexOptimizer);    // Writes statistics in the multi-start scheme statistics file
//     void write_multi_start_scheme_statistics(size_t indexOptimizer, const shared_ptr<Solution_so>& solutionToWrite);   // Writes statistics of a specific solution in the multi-start scheme statistics file

//     // Writes the multi-start scheme statistics of a single solution
//     void write_single_solution_to_multi_start_scheme_statistics(const shared_ptr<Solution_so>& solutionToWrite, const shared_ptr<Optimizer_SO>& optimizerOfSolution);
//     // Writes the multi-start scheme solutions
//     void write_single_solution_to_multi_start_scheme_solutions(const shared_ptr<Solution_so> &solutionToWrite, const shared_ptr<Optimizer_SO> &optimizerOfSolution);

//     // Other functions
//     double determineRunTime();                  // Determines the current run time
//     double getAverageElitistFitness();          // Calculates the average fitness of all best found solutions of the optimizers
//     string convertSolutionNetworkToString(const vector<int> &network);              // Converts a network to a string
//     string convertInstantiationCountToString(const vector<size_t> &instantiations); // Converts a list of instantiation counts to a string
//     string convertBoundariesToString(const vector<vector<double>> &boundaries);     // Converts the boundaries to a string

//     // Print
//     void printOptimizerStatus();        // Prints the status of all optimizers



// };

