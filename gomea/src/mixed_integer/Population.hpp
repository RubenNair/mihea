#pragma once

#include <cmath>
#include <iostream> 
#include <vector>
using namespace std;

#include "gomea/src/mixed_integer/Config.hpp"
#include "gomea/src/mixed_integer/shared.hpp"
#include "gomea/src/fitness/fitness.hpp"
#include "gomea/src/fitness/benchmarks-discrete.hpp"
#include "gomea/src/fitness/benchmarks-mixed.hpp"
#include "gomea/src/common/solution_BN.hpp"
#include "gomea/src/mixed_integer/Solutionset.hpp"
#include "gomea/src/common/partial_solution.hpp"
#include "gomea/src/common/linkage_model.hpp"
#include "gomea/src/mixed_integer/iamalgam.hpp"

namespace gomea{
namespace mixedinteger{

class Population
{
public:
    Config *config;
    fitness_t *problemInstance;
    sharedInformation *sharedInformationPointer;
    size_t GOMEAIndex;
    size_t populationSize;

    iamalgam *iamalgamInstance;
    Solutionset *population;
    Solutionset *offspringPopulation;
    vec_t<int> noImprovementStretches;

    bool terminated;
    double averageFitness;
    size_t numberOfGenerations;
    
    linkage_model_pt FOSInstance = NULL;

    size_t optimizerIndex;
    string optimizerName;
    double optimizerElitistFitness;
    bool optimizerElitistFitnessInitialized;

    Population(Config *config_, fitness_t *problemInstance_, sharedInformation *sharedInformationPointer_, size_t GOMEAIndex_, size_t populationSize_, size_t optimizerIndex, string optimizerName, linkage_model_pt FOSInstance_ = NULL );
    ~Population();

    friend ostream & operator << (ostream &out, const Population &populationInstance);

    void calculateAverageFitness();
    double getFitnessMean();
    double getFitnessVariance();
    double getConstraintValueMean();
    double getConstraintValueVariance();
    solution_mixed *getBestSolution();
    solution_mixed *getWorstSolution();
    bool allSolutionsAreEqual();
    void makeOffspring();
    void learnDiscreteModel();
    void copyOffspringToPopulation();
    void copyPopulationToOffspring();
    void generateOffspring();
    void generateSingleOffspring(int FOS_index);
    void determineFOSOrder();
    void evaluateAllSolutionsInPopulation();
    void evaluateSolution(solution_mixed *solution);
    void evaluateSolution(solution_mixed *solution, solution_mixed *solutionBefore, vec_t<int> &touchedGenes, double fitnessBefore);
    bool GOM(size_t offspringdIndex);
    bool GOMSingleFOS(size_t offspringIndex, size_t FOSIndex);
    bool FI(size_t offspringdIndex);
    bool FISingleFOS(size_t offspringIndex, size_t FOSIndex);
    void updateElitistAndCheckVTR(solution_mixed *solution);
    void checkTimeLimit();
    
    // iAMaLGaM functions
    void learnContinuousModel();
    void generateDiscretePopulation(int FOS_index);
    void checkDiscretePopulationConvergedNotOptimal();
    bool allContinuousEqual(solution_mixed *a, solution_mixed *b);

};

}}