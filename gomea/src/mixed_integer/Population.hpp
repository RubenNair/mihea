#pragma once

#include <cmath>
#include <iostream> 
#include <vector>
using namespace std;

#include "gomea/src/mixed_integer/Config.hpp"
#include "gomea/src/mixed_integer/shared.hpp"
#include "gomea/src/fitness/fitness.hpp"
#include "gomea/src/fitness/benchmarks-discrete.hpp"
#include "gomea/src/common/solution.hpp"
#include "gomea/src/common/partial_solution.hpp"
#include "gomea/src/common/linkage_model.hpp"

namespace gomea{
namespace mixedinteger{

class Population
{
public:
    Config *config;
    fitness_t *problemInstance;
    sharedInformation *sharedInformationPointer;
    size_t GOMEAIndex;
    size_t dPopulationSize;
    size_t cPopulationSize;

    vec_t<solution_t<char>*> dPopulation;
    vec_t<solution_t<char>*> offspringdPopulation;
    vec_t<solution_t<double>*> cPopulation;
    vec_t<solution_t<double>*> offspringcPopulation;
    vec_t<int> noImprovementStretches;

    bool terminated;
    double averageFitness;
    size_t numberOfGenerations;
    
    linkage_model_pt FOSInstance = NULL;

    Population(Config *config_, fitness_t *problemInstance_, sharedInformation *sharedInformationPointer_, size_t GOMEAIndex_, size_t dPopulationSize_, linkage_model_pt FOSInstance_ = NULL );
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
    void generateOffspring();
    void evaluateSolution(solution_t<char> *solution);
    void evaluateSolution(solution_t<char> *solution, solution_t<char> *solutionBefore, vec_t<int> &touchedGenes, double fitnessBefore);
    bool GOM(size_t offspringdIndex);
    bool FI(size_t offspringdIndex);
    void updateElitistAndCheckVTR(solution_t<char> *solution);
    void checkTimeLimit();
};

}}