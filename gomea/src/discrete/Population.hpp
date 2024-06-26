#pragma once

#include <cmath>
#include <iostream> 
#include <vector>
using namespace std;

#include "gomea/src/discrete/Config.hpp"
#include "gomea/src/discrete/shared.hpp"
#include "gomea/src/fitness/fitness.hpp"
#include "gomea/src/fitness/benchmarks-discrete.hpp"
#include "gomea/src/common/solution.hpp"
#include "gomea/src/common/partial_solution.hpp"
#include "gomea/src/common/linkage_model.hpp"

namespace gomea{
namespace discrete{

class Population
{
public:
    Config *config;
    fitness_t *problemInstance;
    sharedInformation *sharedInformationPointer;
    size_t GOMEAIndex;
    size_t populationSize;

    vec_t<solution_t<char>*> population;
    vec_t<solution_t<char>*> offspringPopulation;
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
    solution_t<char> *getBestSolution();
    solution_t<char> *getWorstSolution();
    bool allSolutionsAreEqual();
    void makeOffspring();
    void copyOffspringToPopulation();
    void copyPopulationToOffspring();
    void generateOffspring();
    void generateOffspringSingleFOSElement(int FOS_index);
    void evaluateSolution(solution_t<char> *solution);
    void evaluateSolution(solution_t<char> *solution, solution_t<char> *solutionBefore, vec_t<int> &touchedGenes, double fitnessBefore);
    void determineFOSOrder();
    bool GOM(size_t offspringIndex);
    bool GOMSingleFOSElement(size_t offspringIndex, int FOS_index);
    bool FI(size_t offspringIndex);
    void updateElitistAndCheckVTR(solution_t<char> *solution);
    void checkTimeLimit();
};

}}