#pragma once

#include <vector>
using namespace std;

#include "gomea/src/mixed_integer/gomeaIMS.hpp"
#include "gomea/src/mixed_integer/config.hpp"
#include "gomea/src/mixed_integer/iamalgam.hpp"

namespace gomea{
namespace mixedinteger{

class simpleGAMBIT
{
    public:
    gomeaIMS *gomeaIMSInstance;
    iamalgam *iamalgamInstance;
    Config *config;
    vector<Population*> GAMBITs;
    sharedInformation *sharedInformationInstance = NULL;

 
    int numberOfGAMBITs = 0;
    int maximumNumberOfGAMBITs;
    int numberOfGenerationsGAMBIT = 0;
    int basePopulationSize = 10;
    int currentGAMBITIndex = 0;
    int minimumGAMBITIndex = 0;

    bool isInitialized = false;
    bool hasTerminated = false;

    time_t start_time;
    output_statistics_t output;

    fitness_t *problemInstance = NULL;

    


    simpleGAMBIT();
    simpleGAMBIT(Config *config_);
    ~simpleGAMBIT();
    void initialize();

    void initializeNewGAMBIT();

    void ezilaitini();
    void run();

    bool checkTermination();
    bool checkEvaluationLimitTerminationCriterion();
    bool checkTimeLimitTerminationCriterion();
    bool checkTerminationGAMBIT(int GAMBITIndex);
    void FindCurrentGAMBIT();
};

}}