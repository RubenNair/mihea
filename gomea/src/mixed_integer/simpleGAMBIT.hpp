#pragma once

#include <vector>
using namespace std;

#include "gomea/src/mixed_integer/config.hpp"
#include "gomea/src/mixed_integer/iamalgam.hpp"
#include "gomea/src/mixed_integer/Config.hpp"
#include "gomea/src/mixed_integer/Population.hpp"
#include "gomea/src/mixed_integer/shared.hpp"
#include "gomea/src/mixed_integer/gomea.hpp"
#include "gomea/src/fitness/benchmarks-discrete.hpp"
#include "gomea/src/utils/time.hpp"
#include "gomea/src/common/output_statistics.hpp"

namespace gomea{
namespace mixedinteger{

class simpleGAMBIT
{
    public:
    Config *config;
    vector<Population*> GAMBITs;
    sharedInformation *sharedInformationInstance = NULL;

    int IMSsubgenerationFactor;
    int numberOfGAMBITs = 0;
    int maximumNumberOfGAMBITs;
    int numberOfGenerationsGAMBIT = 0;
    int basePopulationSize = 10;
    int currentGAMBITIndex = 0;
    int minimumGAMBITIndex = 0;

    bool isInitialized = false;
    bool hasTerminated = false;

    time_t start_time;
    clock_t clock_start_time;

    fitness_t *problemInstance = NULL;

    int gen;
    double prevGenElitistFitness;
    bool prevGenElitistInitialized;
    


    simpleGAMBIT();
    simpleGAMBIT(Config *config_);
    ~simpleGAMBIT();
    void initialize();

    void initializeNewGAMBIT();

    void ezilaitini();
    void run();
    void runGeneration();
    void runGeneration(int GAMBITIndex);
    void generationalStepAllGAMBITs();
    void GAMBITGenerationalStepAllGAMBITsRecursiveFold(int GAMBITIndexSmallest, int GAMBITIndexBiggest);


    bool checkTermination();
    bool checkEvaluationLimitTerminationCriterion();
    bool checkTimeLimitTerminationCriterion();
    bool checkTerminationGAMBIT(int GAMBITIndex);
    void FindCurrentGAMBIT();

    double getAverageElitistFitness();
};

}}