#pragma once

#include <string>
#include <cstdlib>
#include <iostream>
#include <random>
#include <chrono>
#include <cassert>
#include <unistd.h>
#include <getopt.h>

using namespace std;

#include "gomea/src/common/linkage_model.hpp"
#include "gomea/src/fitness/fitness.hpp"
#include "gomea/src/utils/data_structure.h"


namespace gomea{
namespace mixedinteger{

typedef gomea::fitness::fitness_t<int> fitness_t;

class Config
{
	void splitString(const string &str, vector<string> &splitted, char delim);
    bool isNumber(const string &s);

public:
	Config();
    ~Config();

    fitness_t *getFitnessClassDiscrete(int problem_index, int number_of_variables);
    void setMethodInitParams(int settingIndex);
    bool parseCommandLine(int argc, char **argv);
    void checkOptions();
    void printUsage();
    void printOverview();
    
	fitness_t *fitness;
	int usePartialEvaluations              = 0,                  
		useParallelGOM		               = 1,                  
		useParallelFOSOrder	               = 0,
		popUpdatesDuringGOM				   = 0,
		fixFOSOrderForPopulation		   = 0,
        AnalyzeFOS                         = 0,
        writeElitists					   = 0,
        printNewElitists                   = 0,
        saveEvaluations                    = 0,
        logDebugInformation                = 0,
        useForcedImprovements              = 0,
        runIndex                           = 0,
        dontUseOffspringPopulation         = 0,
        printHelp                          = 0,
		maximumNumberOfEvaluations		   = -1,
		maximumNumberOfGenerations		   = -1;
	double maximumNumberOfSeconds = -1;
    double     lower_user_range = 0.0,                              /* The initial lower range-bound indicated by the user (same for all dimensions). */
			   upper_user_range = 1.0;                              /* The initial upper range-bound indicated by the user (same for all dimensions). */
    double vtr = 1e+308;
    size_t k = 1, s = 1,   
        FOSIndex = 0;
    double a_value = 1.1;
	int GPUIndex = -1;
	int maximumFOSSetSize = -1;

    string folder = "output"; //"discrete_gomea_output"; // RUBEN changed from "discrete_gomea_output" to "output" 
    //string problemName,
    string FOSName;
    string problemInstancePath = "";

    //long long timelimitMilliseconds = -1,
    bool fix_seed = false;
    long long randomSeed = 0;

    bool useBN = false;
    bool useNormalizedCVars = false;
    bool useOptimalSolution = false;
    bool guaranteedInitSpread = false;
    bool transformCVariables = false;
    bool extraCVarForNumberOfBins = false;
    int discretization_policy_index = 0;
    
    size_t alphabetSize = 2;
    size_t maxArchiveSize = 1000000;
    int maximumNumberOfGOMEAs   = 10, // RUBEN was 100
        maximumNumberOfGAMBITs  = 10, 
        IMSsubgenerationFactor = 4,
        basePopulationSize     = 2,
        numberOfVariables = 10,
        numberOfdVariables = 10,
        numberOfcVariables = 10;
    linkage_config_t *linkage_config = NULL;

    size_t maxDiscretizations = -1;
    size_t maxParents = 6;
    size_t maxInstantiations = 15;
    shared_ptr<DataStructure<double>> data = NULL;

    private:
        int problemIndex = 0;
};

}}