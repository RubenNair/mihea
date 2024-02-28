//
// Created by Damy Ha on 15-Dec-22.
//

#include "gomea/src/fitness/fitness_BN.h"
#include "gomea/src/utils/data_structure.h"
// #include "../../../Continuous_Bayesian_Network_GOMEA/include/Random.h"

#include <iostream>
#include <memory>
#include <sys/stat.h>
#include <utility>

using namespace std;

#ifndef IMPLEMENTATIONS_DENSITY_H
#define IMPLEMENTATIONS_DENSITY_H

namespace gomea{

class Density : public Fitness_BN {
public:
    Density(int problem_index,
            string fitness_function_name,
            const shared_ptr<DataStructure<double>> &data,
            size_t max_number_of_parents,
            size_t max_number_of_discretizations);      // Constructor
    ~Density();                                         // Destructor

    // Methods
//    virtual void computePartialFitness(solution_BN &solution);     // Computes the fitness with partial evaluations

private:
    // Variables
    double penaltyFactor;                               // The weight of the penalty

    /// Partial Evaluations
    vector<vector<vector<int>>> combinedCount;          // [node_i][node_i_parents_are_value_j][node_i_is_value_k]: Counts the number of occurrences of node_i taking a value of 'k' when its parents are taking a value (combination) of 'j'.
    vector<vector<int>> countNodePrior;                 // [node_i][node_i_takes_value_k]: Counts the number of occurrences of node_i taking a value of 'k'.
    vector<vector<int>> countParentCombinations;        // [node_i][number_of_combinations_of_i]: The number of times the parent of node_i take a value combination of 'j'
    vector<int> maxNumberOfCombinationsPerVariable;     // [node_i]: The (maximum) number of combinations i's parents can make, i.e. the 'columns' in the probability table of node_i
    vector<vector<vector<int>>> parentCombinations;     // [node_i][number_of_combinations_of_i][node_k]: The number of combinations node i's parents make according to the data (example: theoretically we can make 4 combinations with 2 True/False parents. If 'TT' is massing, we get less than 4 combinations).
    vector<double> scorePerNode;                        // The log likelihood over the densities per node

    /// PartialEvaluations Previous solution
    vector<vector<vector<int>>> previousCombinedCount;      // [node_i][node_i_parents_are_value_j][node_i_is_value_k]: Counts the previous number of occurrences of node_i taking a value of 'k' when its parents are taking a value (combination) of 'j'.
    vector<vector<int>> previousCountNodePrior;             // [node_i][node_i_takes_value_k]: Counts the previous number of occurrences of node_i taking a value of 'k'.
    vector<vector<int>> previousCountParentCombinations;    // [node_i][number_of_combinations_of_i]: The number of times the previous parent of node_i take a value combination of 'j'
    vector<int> previousMaxNumberOfCombinationsPerVariable; // [node_i]: The previous (maximum) number of combinations i's parents can make, i.e. the 'columns' in the probability table of node_i
    vector<vector<vector<int>>> previousParentCombinations; // [node_i][number_of_combinations_of_i][node_k]: The previous number of combinations node i's parents make according to the data (example: theoretically we can make 4 combinations with 2 True/False parents. If 'TT' is massing, we get less than 4 combinations).
    vector<double> previousScorePerNode;                    // The previous log likelihood over the densities per node


    /// Fitness
    // Computes the fitness of a solution
    tuple<double, vector<double>, double> computeFitnessValue(solution_BN &solution);
    // // Computes the fitness of a solution via partial evaluations
    // tuple<double, vector<double>, double> computePartialFitnessValue(solution_BN &solution, vector<size_t> &nodesToUpdate, vector<size_t> &nodesToRediscretize);

    // Calculates the density score
    tuple<double, vector<double>, double> computeDensityScore(solution_BN &solution,
                                                              vector<vector<vector<int>>> &mijk,
                                                              vector<vector<int>> &m_all,
                                                              vector<vector<int>> &m_prior,
                                                              vector<vector<vector<int>>> &combinations,
                                                              vector<int> &max_number_of_combinations_per_variable);

    // Calculates the log likelihood over densities
    vector<double> loglikelihood_density(solution_BN &solution,
                                         vector<vector<vector<int>>> &mijk,
                                         vector<vector<int>> &m_all,
                                         vector<vector<int>> &m_prior,
                                         vector<int> &max_number_of_combinations_per_variable);

    // /// Partial score
    // tuple<double, vector<double>, double> computePartialDensityScore(solution_BN &solution,
    //                                                                  vector<size_t> &nodesToReassess,
    //                                                                  vector<double> &scoreNode,
    //                                                                  vector<vector<vector<int>>> &mijk,
    //                                                                  vector<vector<int>> &m_all,
    //                                                                  vector<vector<int>> &m_prior,
    //                                                                  vector<vector<vector<int>>> &combinations,
    //                                                                  vector<int> &max_number_of_combinations_per_variable);
    double loglikelihood_density_node(solution_BN &solution,
                                      size_t nodeIndex,
                                      vector<vector<vector<int>>> &mijk,
                                      vector<vector<int>> &m_all,
                                      vector<vector<int>> &m_prior,
                                      vector<int> &max_number_of_combinations_per_variable);
    // // Saves the partial evaluation variables of the previous solution (in case the variables are reset)
    // void savePartialEvaluationVariablesOfPreviousSolution(vector<size_t> &nodesToUpdate, vector<size_t> &nodesToRediscretize);
    // // Resets the partial evaluation variables to the previous solution's variables
    // void resetPartialEvaluationVariablesToPreviousSolution(vector<size_t> &nodesToUpdate, vector<size_t> &nodesToRediscretize);
};

}

#endif //IMPLEMENTATIONS_DENSITY_H
