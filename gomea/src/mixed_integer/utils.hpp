#pragma once

#include <sys/stat.h>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <fstream>
#include <filesystem>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <exception>
#include <cassert>

using namespace std;

// #include "gomea/src/common/solution_mixed.hpp"
#include "gomea/src/common/solution_BN.hpp"
#include "gomea/src/common/gomea_defs.hpp"
#include "gomea/src/fitness/density.h"

namespace gomea{
namespace mixedinteger{

typedef std::chrono::time_point<std::chrono::steady_clock> chtime;

struct archiveRecord
{
  bool isFound = false;
  double value = 0.0;
};

struct hashVector
{ 
    size_t operator()(const vector<int> &vec) const
    { 
        hash <int> hashChar; 
        size_t hash_value = 0;
        for (size_t i = 0; i < vec.size(); ++i) 
            hash_value = hash_value*31 + hashChar(vec[i]); 
        return hash_value; 
    } 
}; 

class solutionsArchive
{
    size_t maxArchiveSize;
public:
    solutionsArchive(size_t maxArchiveSize_): maxArchiveSize(maxArchiveSize_){};
    unordered_map<vector<int>, double, hashVector > archive;
    void checkAlreadyEvaluated(vector<int> &genotype, archiveRecord *result);
    void insertSolution(vector<int> &genotype, double fitness);
};

void prepareFolder(string &folder);
void initElitistFile(string &folder);
void initStatisticsFile(string &folder, bool useBN = false);
void initLogFile(string &folder);
void initBoundaryStatsFile(string &folder);
void writeStatisticsToFile(string &folder, long long numberOfEvaluations, long long time, solution_t<int> *solution, size_t populationSize, bool vtrHit = false);
void writeBNStatisticsToFile(string &folder, long long numberOfEvaluations, long long time, solution_BN *solution, size_t populationSize, bool vtrHit = false);
void writeElitistSolutionToFile(string &folder, long long numberOfEvaluations, long long time, solution_mixed *solution);
void writePopulationToFile(string &folder, vec_t<solution_mixed*> population, string message, bool doLog = true);
void writePopulationBoundaryStatsToFile(string &folder, vec_t<solution_mixed*> population, string message);
void writeBuildingBlocksToFile(string &folder, vec_t<solution_mixed*> population, string message, int k, bool doLog = true);
bool allBuildingBlocksStillExist(vec_t<solution_mixed*> population, int k);
string countBuildingBlocks(vec_t<solution_mixed*> population, int k);
void printPopulation(vec_t<solution_mixed *> &population);
void writeMatrixToFile(string &folder, double **matrix, int rows, int cols, string message, bool doLog = true);
void writeVectorToFile(string &folder, double *vector, int length, string message, bool doLog = true);
void writeMessageToLogFile(string &folder, string message, bool doLog = true);
size_t calculateNumberOfLinks(size_t number_of_nodes);
tuple<vec_t<double>, vec_t<double>> findMaxAndMinValuesInData(vec_t<vec_t<double>> &data);

void writeRunCompletedFile(string &folder, const long long numberOfEvaluations, const clock_t startTime, bool doLog = false);
void copyDataFilesToTargetDir(const string& pathToDataDir, const string &targetDir, const string &problemName, int runIndex);
string determinePathData(const string& pathToDataDir, const string &problemName, int runIndex);
string determinePathInfo(const string& pathToDataDir, const string &problemName, int runIndex);
string determinePathOptimalSolution(const string& pathToDataDir, const string &problemName, int runIndex);
bool copyFile(const string &inputPath, const string &outputPath);
void writeParametersFile(string &folder, Config *config, const Density *fitnessFunction, bool doLog = false);

void write_multi_start_scheme_statistics(string &folder, solution_mixed *elitist, size_t optimizerIndex,
                                            string &optimizerName, size_t number_of_generations, 
                                            clock_t startingRunTime, double avg_elitist_fitness, const vec_t<ColumnDataType> &column_types);
void initialize_multi_start_scheme_statistics_file(string &folder);
void write_single_solution_to_multi_start_scheme_statistics(string &folder, const solution_BN* solutionToWrite, 
                                                            string &optimizerName, size_t optimizerIndex, 
                                                            size_t number_of_generations, clock_t startingRunTime,
                                                            double avg_elitist_fitness);
void initialize_multi_start_scheme_solutions_file(string &folder);
void write_single_solution_to_multi_start_scheme_solutions(string &folder, const solution_BN* solutionToWrite, 
                                                            string &optimizerName, size_t optimizerIndex, 
                                                            size_t number_of_generations, clock_t startingRunTime, const vec_t<ColumnDataType> &column_types);

string convertSolutionNetworkToString(const vector<int> &network);
string convertInstantiationCountToString(const vector<size_t> &instantiations);
string convertBoundariesToString(const vector<vector<double>> &boundaries, const vec_t<ColumnDataType> &column_types);
vector<size_t> getNumberOfInstantiations(const solution_BN *solution); // Calculate the number of instantiations per node based on the amount of boundaries

}}