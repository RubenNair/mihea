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

#include "gomea/src/common/solution_mixed.hpp"
#include "gomea/src/common/gomea_defs.hpp"

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
    size_t operator()(const vector<char> &vec) const
    { 
        hash <char> hashChar; 
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
    unordered_map<vector<char>, double, hashVector > archive;
    void checkAlreadyEvaluated(vector<char> &genotype, archiveRecord *result);
    void insertSolution(vector<char> &genotype, double fitness);
};

void prepareFolder(string &folder);
void initElitistFile(string &folder);
void initStatisticsFile(string &folder);
void initLogFile(string &folder);
void writeStatisticsToFile(string &folder, long long numberOfEvaluations, long long time, solution_t<char> *solution, size_t populationSize, bool vtrHit = false);
void writeElitistSolutionToFile(string &folder, long long numberOfEvaluations, long long time, solution_mixed *solution);
void writePopulationToFile(string &folder, vec_t<solution_mixed*> population, string message, bool doLog = true);
void writeBuildingBlocksToFile(string &folder, vec_t<solution_mixed*> population, string message, int k, bool doLog = true);
bool allBuildingBlocksStillExist(vec_t<solution_mixed*> population, int k);
string countBuildingBlocks(vec_t<solution_mixed*> population, int k);
void printPopulation(vec_t<solution_mixed *> &population);
void writeMatrixToFile(string &folder, double **matrix, int rows, int cols, string message, bool doLog = true);
void writeVectorToFile(string &folder, double *vector, int length, string message, bool doLog = true);
void writeMessageToLogFile(string &folder, string message, bool doLog = true);

}}