// //
// // Created by Damy Ha on 12-Mar-23.
// //

// #include "gomea/src/mixed_integer/MultiStartScheme_statistics.h"

// /**
//  * Writes statistics of the current generation, given the best found solution in the archive.
//  * @param indexOptimizer The index of the optimizer
//  */
// void MultiStartScheme::write_multi_start_scheme_statistics(size_t indexOptimizer) {
//     // Retrieve the best fitting solution from the archive
//     shared_ptr<Solution_so> solution = this->determineSolutionToWrite();

//     // Write the solution
//     this->write_multi_start_scheme_statistics(indexOptimizer, solution);
// }

// /**
//  * Writes statistics of the current generation, given a specific solution
//  * @param indexOptimizer The index of the optimizer
//  * @param solutionToWrite The specific solution to write
//  */
// void MultiStartScheme::write_multi_start_scheme_statistics(size_t indexOptimizer,
//                                                            const shared_ptr<Solution_so>& solutionToWrite) {
//     // Write the solution
//     if (solutionToWrite) {
//         // Prepare variables
//         shared_ptr<Optimizer_SO> optimizer = this->optimizers[indexOptimizer];  // Determine the optimizer that wrote the solution

//         // Write the statistics of the solution
//         if (write_generational_statistics) {
//             this->write_single_solution_to_multi_start_scheme_statistics(solutionToWrite, optimizer);
//         }

//         // Write the solution
//         if (write_generational_statistics) {
//             this->write_single_solution_to_multi_start_scheme_solutions(solutionToWrite, optimizer);
//         }
//     }

// }


// /**
//  * Determines the solution to write to the MSS statistics files
//  * @return The solution to write
//  */
// shared_ptr<Solution_so> MultiStartScheme::determineSolutionToWrite() {
//     // Find the most fitting solution to write to the file
//     shared_ptr<Solution_so> solutionToWrite = nullptr;
//     double eliteAccuracy, eliteSensitivity, eliteAverage_hamming_distance, eliteArc_ratio;
//     for(const auto &solution : archive.getArchive()) {
//         // Check if the solutions was obtained within the time budget
//         double runTime = double(solution->getTimeStamp() - this->startingRunTime) / CLOCKS_PER_SEC;
//         if (maximum_time_in_seconds > 0 and runTime > maximum_time_in_seconds) {
//             continue;
//         }

//         // Check if the solution was obtained within the evaluation budget
//         if (maximum_number_of_evaluations > 0 and solution->getNumberOfEvaluations() > maximum_number_of_evaluations) {
//             continue;
//         }

//         // Check if the solution has a better score than the current best solution
//         if (!solutionToWrite) {
//             solutionToWrite = solution;
//         } else {
//             // Determine the score metric that is most important
//             double accuracy, sensitivity, average_hamming_distance, arc_ratio;
//             this->fitnessFunction->evaluateNetwork(*solution, accuracy, sensitivity, average_hamming_distance, arc_ratio);

//             // Update the best solution
//             if (accuracy > eliteAccuracy) {
//                 solutionToWrite = solution;
//                 eliteAccuracy = accuracy;
//             }
//         }
//     }

//     return solutionToWrite;
// }


// /////////////////////////////////////
// /// Multi-start scheme statistics ///
// /////////////////////////////////////
// /**
//  * Initializes the generational statistics file
//  */
// void MultiStartScheme::initialize_multi_start_scheme_statistics_file() {
//     // Create and open statistics file
//     shared_ptr<ofstream> statisticsStream(new ofstream);
//     // Open
//     this->streamMultiStartSchemeStatistics = statisticsStream;
//     this->streamMultiStartSchemeStatistics->open(this->pathMultiStartSchemeStatistics, ofstream::out | ofstream::trunc);

//     // Add header to statistics file
// #ifdef DEBUG_STATISTICS
//     *statisticsStream
//             << "     Algorithm  Gen  OptInd  TotalEvaluations          Time          "
//                "Best-objective       Average-objective       Accuracy    Sensitivity AvgHammingDistance       "
//                "LogLikelihood       LLDifference    Dist.DistanceKL       ArcRatio      AvgNoDisc        AvgDiscDist        "
//                "MaxDiscDist      Ham.Boundaries\n";
// #else
//     *statisticsStream
//             << "     Algorithm  Gen  OptInd  TotalEvaluations          Time          "
//                "Best-objective       Average-objective\n";

// #endif
// }



// /**
//  * Writes a single solution to the multi start scheme statistics file.
//  * @param solutionToWrite The solution to write to the file
//  * @param optimizerOfSolution The optimizer that created the solution
//  */
// void MultiStartScheme::write_single_solution_to_multi_start_scheme_statistics(const shared_ptr<Solution_so>& solutionToWrite,
//                                                                               const shared_ptr<Optimizer_SO>& optimizerOfSolution) {
//     // Perform check
//     if (!this->streamMultiStartSchemeStatistics->is_open()) {
//         throw runtime_error("Multi Start Scheme statistics file was not open.");
//     }

//     // Load variables
//     string optimizerName = optimizerOfSolution->getOptimizerName();
//     double runTime = double(solutionToWrite->getTimeStamp() - this->startingRunTime) / CLOCKS_PER_SEC;

//     // Write to statistics file
//     *streamMultiStartSchemeStatistics
//             << setw(14) << optimizerName
//             << " " << setw(4) << number_of_generations
//             << " " << setw(7) << optimizerOfSolution->getOptimizerIndex()
//             << " " << setw(17) << solutionToWrite->getNumberOfFullEvaluations()
//             << " " << setw(13) << scientific << setprecision(3) << runTime
//             << " " << setw(23) << scientific << setprecision(16) << solutionToWrite->getFitness()
//             << " " << setw(23) << scientific << setprecision(16) << this->getAverageElitistFitness()
//             << endl;
// }


// ////////////////////////////////////
// /// Multi-start scheme solutions ///
// ////////////////////////////////////
// /**
//  * Initializes the generational solution file
//  */
// void MultiStartScheme::initialize_multi_start_scheme_solutions_file() {
//     // Create and open statistics file
//     shared_ptr<ofstream> statisticsStream(new ofstream);
//     // Open
//     this->streamMultiStartSchemeSolutions = statisticsStream;
//     this->streamMultiStartSchemeSolutions->open(this->pathMultiStartSchemeSolutions, ofstream::out | ofstream::trunc);

//     // Add header to statistics file
//     *statisticsStream
//     << "     Algorithm  Gen  OptInd  TotalEvaluations          Time          "
//        "Best-objective "
//        "Network Discretizations Edges\n";
// }

// /**
//  * Writes a single solution to the multi start scheme solution file.
//  * @param solutionToWrite The solution to write to the file
//  * @param optimizerOfSolution The optimizer that created the solution
//  */
// void MultiStartScheme::write_single_solution_to_multi_start_scheme_solutions(const shared_ptr<Solution_so> &solutionToWrite,
//                                                                              const shared_ptr<Optimizer_SO>& optimizerOfSolution) {
//     // Perform check
//     if (!this->streamMultiStartSchemeSolutions->is_open()) {
//         throw runtime_error("Multi Start Scheme solutions file was not open.");
//     }

//     // Load variables
//     string optimizerName = optimizerOfSolution->getOptimizerName();
//     double runTime = double(solutionToWrite->getTimeStamp() - this->startingRunTime) / CLOCKS_PER_SEC;

//     // Write to statistics file
//     *streamMultiStartSchemeSolutions
//             << setw(14) << optimizerName
//             << "/" << setw(4) << number_of_generations
//             << "/" << setw(7) << optimizerOfSolution->getOptimizerIndex()
//             << "/" << setw(17) << solutionToWrite->getNumberOfFullEvaluations()
//             << "/" << setw(13) << scientific << setprecision(3) << runTime
//             << "/" << setw(23) << scientific << setprecision(16) << solutionToWrite->getFitness()
//             << "/" << convertSolutionNetworkToString(solutionToWrite->getNetworkParameters())
//             << "/" << convertInstantiationCountToString(solutionToWrite->getDiscretizationPolicy()->getNumberOfInstantiations())
//             << "/" << convertBoundariesToString(solutionToWrite->getDiscretizationPolicy()->getBoundaries())
//             << endl;
// }


// /////////////////////////////
// /// Statistical Functions ///
// /////////////////////////////
// /**
//  * Calculates the average best solution fitness
//  * @return The average fitness over all elitist solutions
//  */
// double MultiStartScheme::getAverageElitistFitness() {
//     double sum_fitness = 0.0;
//     for (const auto &optimizer : this->optimizers) {
//         sum_fitness += optimizer->getElitistSolution()->getFitness();
//     }

//     return sum_fitness / (double) this->optimizers.size();
// }

// ///////////////////////
// /// Other functions ///
// ///////////////////////
// /**
//  * Determines the current run time
//  * @return The run time
//  */
// double MultiStartScheme::determineRunTime() {
//     // Determine time that has passed
//     clock_t currentTime = clock();
//     double result = double(currentTime - this->startingRunTime) / CLOCKS_PER_SEC;
//     return result;
// }

// /**
//  * Converts a network to a string
//  * @param network The network to convert
//  * @return The network as a string
//  */
// string MultiStartScheme::convertSolutionNetworkToString(const vector<int> &network) {
//     // Initialize string stream
//     stringstream ss;
//     ss << "[";

//     // Add the elements of the network
//     for (size_t i = 0; i < network.size(); ++i) {
//         ss << network[i];
//         // Skip the comma for the last item
//         if (i != network.size() - 1) {
//             ss << ",";
//         }
//     }

//     ss << "]";
//     return ss.str();
// }

// /**
//  * Converts a list of instantiation counts to a string
//  * @param instantiations The instantiations to write as a string
//  * @return The instantiation count as a string
//  */
// string MultiStartScheme::convertInstantiationCountToString(const vector<size_t> &instantiations) {
//     // Initialize string stream
//     stringstream ss;
//     ss << "[";

//     // Add the elements of the network
//     for (size_t i = 0; i < instantiations.size(); ++i) {
//         ss << instantiations[i];
//         // Skip the comma for the last item
//         if (i != instantiations.size() - 1) {
//             ss << ",";
//         }
//     }

//     ss << "]";
//     return ss.str();

// }

// /**
//  * Converts the boundaries to a string
//  * @param boundaries The boundaries
//  * @return The boundaries as a string
//  */
// string MultiStartScheme::convertBoundariesToString(const vector<vector<double>> &boundaries) {
//     // Initialize string stream
//     stringstream ss;
//     ss << "[";

//     // Go over each node
//     for (size_t i = 0; i < boundaries.size(); ++i) {
//         ss << "[";

//         // Add the values of the boundaries
//         for (size_t j = 0; j < boundaries[i].size(); ++j) {
//             ss << scientific << setprecision(16) << boundaries[i][j];

//             // Skip the last item
//             if (j != boundaries[i].size() - 1) { ss << ","; }
//         }

//         ss << "]";
//         // Skip the last item
//         if (i != boundaries.size() - 1) { ss << ";";}
//     }

//     ss << "]";

//     return ss.str();
// }