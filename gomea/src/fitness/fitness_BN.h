//
// Created by Damy Ha on 03-Oct-22.
//

#ifndef IMPLEMENTATIONS_FITNESS_H
#define IMPLEMENTATIONS_FITNESS_H

#include <cstddef>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <tuple>
#include <utility>
#include <unordered_set>

#include "gomea/src/common/solution_BN.hpp"
#include "gomea/src/common/gomea_defs.hpp"
// #include "../../../Continuous_Bayesian_Network_GOMEA/include/util.h"

using namespace std;

namespace gomea {

class Fitness_BN {
public:
    // Constructor for non data dependent problems
//    Fitness_BN(int problem_index,
//            string fitness_function_name,
//            size_t number_of_nodes,
//            size_t number_of_links,
//            size_t number_of_nodes_to_discretize,
//            size_t maximum_number_of_parents,
//            size_t max_number_of_discretizations);

    // Constructor for data dependent problems
    Fitness_BN(int problem_index,
            string fitness_function_name,
            const shared_ptr<DataStructure<double>> &data,
            size_t max_number_of_parents,
            size_t max_number_of_discretizations);

    ~Fitness_BN();     // Destructor

    // Virtual functions
    void computeFitness(solution_BN &solution);                             // Computes the fitness (Calls a private function)
    void computeFitnessWithoutCountingEvaluations(solution_BN &solution);   // Computes the fitness without counting it as an evaluation
    void computePartialFitness(solution_BN &solution);                      // Computes the fitness with partial evaluations

    // Functions
    void resetToOriginalData();

    // Initializes the matrices needed to calculate the occurrences
    void initializeCountMatrices(solution_BN &solution, vector<vector<int>> &m_prior, vector<vector<int>> &m_all, vector<int> &max_number_of_combinations_per_variable, vector<vector<vector<int>>> &combinations, vector<vector<vector<int>>> &mijk);
    // Compute the prior and conditional probabilities
    void calculatePriorAndConditionalCounts(solution_BN &solution, vector<vector<int>> &m_prior, vector<vector<int>> &m_all, vector<int> &max_number_of_combinations_per_variable, vector<vector<vector<int>>> &combinations, vector<vector<vector<int>>> &mijk);
    // Initializes matrix m_prior and m_ijk
    void computeSampleValue(int node_index, const solution_BN& solution, int maximum_number_of_parent_combinations_of_i, vector<vector<vector<int>>> &combinations, vector<vector<vector<int>>> &mijk,vector<vector<int>> &p_prior);
    // Used to determine the observed number of combinations of variable for a particular node. See the implementation for more details
    void findAllCombinationsValues(int node_index, const solution_BN &solution, vector<vector<vector<int>>> &combinations, vector<int> &number_of_combinations_per_variable);
    // Used to check if parent values matches a combination vector
    bool sampleMatchesAnyCombination(vector<vector<int>> combinations, size_t combinationsFound, const vector<double>& sample_data, vector<int> &combinationNodeIndices, int numberOfNodesInCombination);
    bool sampleMatchesCombination(vector<int> combination, vector<double> sample_data, vector<int> combinationNodeIndices, int numberOfNodesInCombination);
    int findCombinationIndex(vector<vector<int>> validCombinations, const vector<double>& sample_data, const vector<int>& combinationNodeIndices );

    /// Partial evaluations
    // Compute the nodes that have changed w.r.t the previous solution
    tuple<vector<size_t>, vector<size_t>> computeNodesThatHaveChanged(solution_BN &solution);
    // Reinitialize the matrices needed to calculate occurrences
    void reinitializeCountMatrices(solution_BN &solution, const vector<size_t> &nodesToReassess, vector<vector<int>> &mPrior, vector<vector<int>> &mAll, vector<int> &maxNumberOfCombinationsPerVariable, vector<vector<vector<int>>> &combinations, vector<vector<vector<int>>> &mijk);
    // Recompute occurrences
    void recalculateCounts(solution_BN &solution, const vector<size_t> &nodesToRecount, vector<vector<int>> &mPrior, vector<vector<int>> &mAll, vector<int> &maxNumberOfCombinationsPerVariable, vector<vector<vector<int>>> &combinations, vector<vector<vector<int>>> &mijk);

    // Getters
    const string &getFitnessFunctionName() const;
    const string &getFitnessFunctionBaseName() const;
    const string &getFitnessFunctionType() const;
    size_t getNumberOfNodes() const;
    size_t getNumberOfLinks() const;
    size_t getNumberOfNodesToDiscretize() const;
    const vector<ColumnDataType> &getDataTypeNodes() const;
    size_t getMaximumNumberOfParents() const;
    size_t getMaxNumberOfDiscretizations() const;
    const shared_ptr<DataStructure<double>> &getOriginalData() const;
    const shared_ptr<DataStructure<double>> &getData() const;       // Should only be used for debugging
    size_t getSampleSize();
    size_t getNumberOfFullEvaluations() const;
    size_t getNumberOfEvaluations() const;
    bool getUsePartialEvaluations() const;
    bool hasOptimalNetwork() const;
    bool hasOptimalBoundariesAvailable() const;
    const vector<int> &getOptimalNetwork() const;
    const vector<vector<double>> &getOptimalBoundaries() const;
    const vector<size_t> &getContinuousNodeIndices() const;
    const vector<size_t> &getDiscreteNodeIndices() const;
    const vector<size_t> &getInitialNumberOfInstantiations() const;
    shared_ptr<solution_BN> getGoldenSolution();
    bool getPostProcessing() const;
    const shared_ptr<mt19937> &getRng() const;

    // Setters
    void setFitnessFunctionBaseName(const string &fitnessFunctionBaseName);
    void setData(const shared_ptr<DataStructure<double>> &data);
    void setOptimalNetwork(const vector<int> &optimalSolution);
    void setOptimalBoundaries(const vector<vector<double>> &optimalBoundaries);
    void setUsePartialEvaluations(bool usePartialEvaluations);
    void setPostProcessing(bool postProcessing);
    void setPreviousSolution(const shared_ptr<solution_BN> &previousSolution);      // Only used for testing
    void setRng(const shared_ptr<mt19937> &rng);

    // Printing
    void printOptimalNetwork();
    void printOptimalBoundaries();

    /// Deprecated statistical functions
    // Statistics
    void evaluateNetwork(const solution_BN& solution, double &accuracy, double &sensitivity, double &average_hamming_distance, double &arc_ratio);                          // Evaluates the network
    void evaluateDiscretization(const solution_BN& solution, double &avgProposedNumberOfDiscretizations, double &averageBoundaryDistance, double &maxBoundaryDistance, int &hammingDistanceBoundaries);     // Evaluates the discretization

    // Statistics - Network evaluation
    double calculateAccuracy(size_t true_positive, size_t true_negative, size_t false_positive, size_t false_negative);                 // Calculates the accuracy of the network
    double calculateSensitivity(size_t true_positive, size_t true_negative, size_t false_positive, size_t false_negative);              // Calculates the sensitivity of the network
    double calculateAverageHammingDistance(size_t true_positive, size_t true_negative, size_t false_positive, size_t false_negative);   // Calculates the average hamming distance of the network
    double calculateRatioOfCorrectArcs(size_t true_negative, size_t arc_present);                                                       // Calculates the correct number of arcs

    // Statistics - Discretization evaluation
    double logLikelihoodDifference(solution_BN &solution);                                                              // Calculates the log likelihood difference between the golden network
    double log_likelihood_score(solution_BN &solution);                                                                 // Calculates the log likelihood of a solution
    vector<double> calculateClosestDistances(vector<double> nodeBoundaries, vector<double> nodeOptimalBoundaries);      // Calculates the closest proposed distance w.r.t. the optimal discretization boundary
    vector<double> calculatePercentualDifferenceBoundaries(vector<double> closestDistance, vector<double> optimalBoundary);                 // Normalizes the boundaries as a percentage of error.
    double calculateAverageBoundaryDistance(const vector<vector<double>>& minimumDistancesPerNode);                     // Calculates the average distance per optimal discretization boundary
    double calculateMaxBoundaryDistance(const vector<vector<double>>& minimumDistancesPerNode);                         // Calculates the worst distance of all discretization boundaries
    double calculateAvgNumberOfProposedBoundaries(const solution_BN& solution);                                         // Calculates the average number of boundaries per node
    int calculateHammingDistanceBoundaries(const solution_BN& solution, vector<size_t> optimalBoundaryCount);           // Calculates the hamming distance between the optimal boundary count

    // Statistics - CPD Kullback leibler
    double calculateKLDivergenceScore(solution_BN &solution, double &scoreKL, double &scoreLS);                         // Calculates the divergence and distance scores
    void prepareCPDVariables(solution_BN &solution,
                             vector<vector<vector<double>>> &resultConditionalDensity,
                             shared_ptr<DataStructure<double>> &resultDiscretizedData,
                             vector<vector<vector<int>>> &resultParentCombinations);                                    // Prepares the variables needed for the KL score
    vector<vector<vector<double>>> computeConditionalProbabilityDensity(solution_BN &solution,
                                                                        vector<vector<vector<int>>> &mijk,
                                                                        vector<vector<int>> &m_all,
                                                                        vector<vector<int>> &m_prior,
                                                                        vector<int> &max_number_of_combinations_per_variable);      // Calculates the conditional probability density table
    double calculateKLOverDensity(solution_BN &desiredSolution,
                                  solution_BN &givenSolution,
                                  const vector<vector<vector<double>>> &desiredDensity,
                                  const vector<vector<vector<double>>> &solutionDensity,
                                  const shared_ptr<DataStructure<double>> &dataDesired,
                                  const shared_ptr<DataStructure<double>> &dataSolution,
                                  const vector<vector<vector<int>>> &desiredCombinations,
                                  const vector<vector<vector<int>>> &solutionCombinations);     // Calculates the Kullback Leibler divergence
    double calculateL2overDensities(solution_BN &desiredSolution,
                                    solution_BN &givenSolution,
                                    const vector<vector<vector<double>>> &desiredDensity,
                                    const vector<vector<vector<double>>> &solutionDensity,
                                    const shared_ptr<DataStructure<double>> &dataDesired,
                                    const shared_ptr<DataStructure<double>> &dataSolution,
                                    const vector<vector<vector<int>>> &desiredCombinations,
                                    const vector<vector<vector<int>>> &solutionCombinations);   // Calculates the Least square distane

    // Statistics - CPT evaluation
    double calculateL1DistributionDistance(solution_BN &solution);                                                                                        // Calculates the distance between distributions of the golden and given solution
    vector<vector<vector<double>>> computeCPT(solution_BN &solution);                                                                                   // Computes the conditonal probability tables of every node
    vector<vector<vector<double>>> normalizeToCPT(solution_BN &solution, const vector<vector<int>> &m_all, const vector<vector<vector<int>>> &mijk);    // Normalizes the counting of data to a conditional probability table
    double calculateL1DifferenceDistributions(const vector<double> &givenDistribution, const vector<double> &targetDistribution);                       // Calculates the L1 difference between distributions
    vector<vector<double>> recursiveAppendZeroToDistribution(vector<vector<double>> currentDistributions, size_t targetSize);                           // Adds zeros to the distributions (in all possible combinations)
    vector<vector<double>> removeDuplicateDistributions(const vector<vector<double>>& distributions);                                                   // Removes duplicate distributions

protected:
    // Variables
    int problemIndex;                           // The index of the problem (used for debugging)
    string fitnessFunctionType;                 // The type of the fitness function (e.g. BDeu, MDL)
    string fitnessFunctionName;                 // The name of the fitness function
    string fitnessFunctionBaseName;             // The base name of the fitness function

    // Network variables
    size_t numberOfNodes;                       // The number of nodes in the problem
    size_t numberOfLinks;                       // The number of links in the network structure
    size_t numberOfNodesToDiscretize;           // The number of nodes to discretize
    vector<ColumnDataType> dataTypeNodes;       // The data type of a node (Continuous or Discrete)
    vector<size_t> continuousNodeIndices;       // The nodes to discretize
    vector<size_t> discreteNodeIndices;         // The nodes that are discrete

    // RNG
    shared_ptr<mt19937> rng;                    // Pointer to the random number generator

    // Variables partial evaluations
    shared_ptr<solution_BN> previousSolution;   // The previous solution of which the fitness has been calculated

    // Data
    shared_ptr<DataStructure<double>> originalData;     // Contains the original data
    shared_ptr<DataStructure<double>> data;             // Data to use for computing scores

    // Evaluations
    bool usePartialEvaluations;                 // Use partial evaluations
    size_t numberOfFullEvaluations;             // The number of full evaluations done
    double numberOfEvaluations;                 // The number of (partial) evaluations

    // Limits
    size_t maxNumberOfParents;                  // The maximum number of parents per node
    size_t maxNumberOfDiscretizations;          // The maximum number of discretization per node

    // Golden Network
    bool optimalNetworkAvailable;               // Whether an optimal network is available
    bool optimalBoundariesAvailable;            // Whether the boundaries for the nodes to be discretized is available
    shared_ptr<solution_BN> optimalSolution;    // The optimal solution
    vector<int> optimalNetwork;                 // The optimal network
    vector<vector<double>> optimalBoundaries;   // The optimal boundaries
    double loglikelihood;                       // The log likelihood. Todo: Move this to the solution if possible.

    // Post-processing
    bool postProcessing;                      // Write solution files for post-processing

    /// Fitness function
    // Implementation specific fitness computation
    virtual tuple<double, vector<double>, double> computeFitnessValue(solution_BN &solution);
    // // Implementation specific partial fitness computation
    // virtual tuple<double, vector<double>, double> computePartialFitnessValue(solution_BN &solution, vector<size_t> &nodesToUpdate, vector<size_t> &nodesToRediscretize);

    // /// Partial evaluation functions
    // // Computes the node indices which have changed their (parent) topology
    // virtual void computeNodesThatHaveChangedTopology(solution_BN &solution, vector<size_t> &result);
    // // Computes the nodes that have changed the number of discretizations
    // void computeNodesThatHaveChangedDiscretization(solution_BN &solution,
    //                                                vector<size_t> &resultNodesToUpdate,
    //                                                vector<size_t> &resultNodesToRediscretize);
    // // Saves the partial evaluation variables of the previous solution (in case the variables are reset)
    // virtual void savePartialEvaluationVariablesOfPreviousSolution(vector<size_t> &nodesToUpdate, vector<size_t> &nodesToRediscretize);
    // // Resets the partial evaluation variables to the previous solution's variables
    // virtual void resetPartialEvaluationVariablesToPreviousSolution(vector<size_t> &nodesToUpdate, vector<size_t> &nodesToRediscretize);


    // Discretize data
    void discretizeData(solution_BN &solution);                                                 // Discretizes data
    vector<vector<double>> discretizeData(solution_BN &solution, DataStructure<double> &data);  // Discretizes data using continuous part learned by iAMaLGaM
    vector<double> discretizeData(solution_BN &solution, size_t nodeIndex, int continuousIndex, const vector<double>& dataOfNode);
    vector<size_t> getNumberOfInstantiations(solution_BN &solution);                            // Calculate the number of instantiations per node based on the amount of boundaries
    void rediscretizeNodes(solution_BN &solution, const vector<size_t>& nodesToRediscretize);   // Discretizes data for partial evaluations

    // Helper functions
    void separateNodeType();                    // Determines the nodes to discretize
    vector<double> normalizeBoundaries(vector<double> boundaries, vector<double> continuousData);   // Normalize the boundaries
    double log_likelihood(solution_BN &solution, vector<vector<vector<int>>> &mijk, vector<vector<int>> &m_all, vector<int> &max_number_of_combinations_per_variable);  // The (normalized) log likelihood
};

}

#endif //IMPLEMENTATIONS_FITNESS_H
