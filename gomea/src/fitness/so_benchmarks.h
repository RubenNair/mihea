//
// Created by Damy Ha on 10-Oct-22.
//

#ifndef IMPLEMENTATIONS_SO_BENCHMARKS_H
#define IMPLEMENTATIONS_SO_BENCHMARKS_H

#include <iostream>
#include <fstream>
#include <algorithm>

#include "gomea/src/utils/data_structure.h"


// // Retrieves the objective function given the problem index
// shared_ptr<Fitness> getSOObjectivePointer(int index, int scoreMetricIndex, int runIndex, const shared_ptr<mt19937>& rng, bool initializeData);

// Print functions
void printAllInstalledProblems();

/////////////////////////
/// Fitness functions ///
/////////////////////////
// shared_ptr<Fitness> unimplemented_problem(          int problemIndex, int scoreMetricIndex, int runIndex, const shared_ptr<mt19937>& rng, bool initializeData);
// shared_ptr<Fitness> test(                           int problemIndex, int scoreMetricIndex, int runIndex, const shared_ptr<mt19937>& rng, bool initializeData);
// shared_ptr<Fitness> chain_2_classes(                int problemIndex, int scoreMetricIndex, int runIndex, const shared_ptr<mt19937>& rng, bool initializeData);
// shared_ptr<Fitness> chain_3_classes(                int problemIndex, int scoreMetricIndex, int runIndex, const shared_ptr<mt19937>& rng, bool initializeData);
// shared_ptr<Fitness> fork_2_classes(                 int problemIndex, int scoreMetricIndex, int runIndex, const shared_ptr<mt19937>& rng, bool initializeData);
// shared_ptr<Fitness> fork_3_classes(                 int problemIndex, int scoreMetricIndex, int runIndex, const shared_ptr<mt19937>& rng, bool initializeData);
// shared_ptr<Fitness> collider_2_classes(             int problemIndex, int scoreMetricIndex, int runIndex, const shared_ptr<mt19937>& rng, bool initializeData);
// shared_ptr<Fitness> collider_3_classes(             int problemIndex, int scoreMetricIndex, int runIndex, const shared_ptr<mt19937>& rng, bool initializeData);
// shared_ptr<Fitness> asia(                           int problemIndex, int scoreMetricIndex, int runIndex, const shared_ptr<mt19937>& rng, bool initializeData);
// shared_ptr<Fitness> asia_continuous(                int problemIndex, int scoreMetricIndex, int runIndex, const shared_ptr<mt19937>& rng, bool initializeData);
// shared_ptr<Fitness> insurance(                      int problemIndex, int scoreMetricIndex, int runIndex, const shared_ptr<mt19937>& rng, bool initializeData);
// shared_ptr<Fitness> alarm(                          int problemIndex, int scoreMetricIndex, int runIndex, const shared_ptr<mt19937>& rng, bool initializeData);
// shared_ptr<Fitness> heparii(                        int problemIndex, int scoreMetricIndex, int runIndex, const shared_ptr<mt19937>& rng, bool initializeData);
// shared_ptr<Fitness> independent_5nodes(             int problemIndex, int scoreMetricIndex, int runIndex, const shared_ptr<mt19937>& rng, bool initializeData);
// shared_ptr<Fitness> dependent_5nodes(               int problemIndex, int scoreMetricIndex, int runIndex, const shared_ptr<mt19937>& rng, bool initializeData);
// shared_ptr<Fitness> dependent_5nodes_uniform_range( int problemIndex, int scoreMetricIndex, int runIndex, const shared_ptr<mt19937>& rng, bool initializeData);
// shared_ptr<Fitness> uniform_range_27nodes(          int problemIndex, int scoreMetricIndex, int runIndex, const shared_ptr<mt19937>& rng, bool initializeData);
// shared_ptr<Fitness> uniform_range_37nodes(          int problemIndex, int scoreMetricIndex, int runIndex, const shared_ptr<mt19937>& rng, bool initializeData);
// shared_ptr<Fitness> network_table_7(                int problemIndex, int scoreMetricIndex, int runIndex, const shared_ptr<mt19937>& rng, bool initializeData);
// shared_ptr<Fitness> uniform_range_11nodes(          int problemIndex, int scoreMetricIndex, int runIndex, const shared_ptr<mt19937>& rng, bool initializeData);
// shared_ptr<Fitness> random_uniform_network(int generationMethod, int problemIndex, int scoreMetricIndex, int runIndex, const shared_ptr<mt19937>& rng, bool initializeData);

// Data paths
string determinePathInfo(const string &problemName, int runIndex);
string determinePathData(const string &problemName, int runIndex);
string determinePathOptimalSolution(const string &problemName, int runIndex);
void copyDataFilesToTargetDir(const string &targetDir, const string &problemName, int runIndex);

// Helper functions
bool copyFile(const string &inputPath, const string &outputPath);
// shared_ptr<Fitness> createGenericFitnessFunction(int problemIndex, int scoreMetricIndex, const string& functionName, const string& baseName, const shared_ptr<DataStructure<double>> &data, size_t maxNumberOfParents, size_t maxNumberOfDiscretizations);
shared_ptr<DataStructure<double>> initializeDataFromPath(bool initializeData, const string& baseName, int runIndex);
void retrieveOptimalSolution(const string& pathBestSolution, string& stringOriginalNetwork, vector<vector<double>> &originalBoundaries);
// void setOptimalSolutionToFitnessFunction(const shared_ptr<Fitness>& fitnessFunction, const string& baseName, int runIndex );

vector<int> convertNetworkStringToIntVector(string input);

bool path_exists(string path_to_file); 
#endif //IMPLEMENTATIONS_SO_BENCHMARKS_H
