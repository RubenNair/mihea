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
#include "gomea/src/common/solution_mixed.hpp"
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
    vec_t<solution_mixed*> population;
    vec_t<solution_mixed*> offspringPopulation;
    vec_t<int> noImprovementStretches;

    bool terminated;
    double averageFitness;
    size_t numberOfGenerations;
    
    linkage_model_pt FOSInstance = NULL;

    Population(Config *config_, fitness_t *problemInstance_, sharedInformation *sharedInformationPointer_, size_t GOMEAIndex_, size_t populationSize_, linkage_model_pt FOSInstance_ = NULL );
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
    // TODO RUBEN createRandomSolution from GAMBIT sould be method in this class -> I think this is the randomInit in solution_mixed,
    //      might just need a function here that calls that function for all individuals in population.
    //      Actually, currently randomInit is called for each individual when the population is created, so no longer necessary?
    
    // iAMaLGaM functions
    // void updatePopulation(int population_index, Population *currGAMBIT, bool onlyObjectiveAndConstraints = false);
    void learnContinuousModel();
    void generateDiscretePopulation(int FOS_index);
    void checkDiscretePopulationConvergedNotOptimal();
    void checkForDuplicate(string message);
    void checkIndividualSolvedDiscrete();
    bool allContinuousEqual(solution_mixed *a, solution_mixed *b);
    // void generateNewContinuousPopulation();

};

}}