//
// Created by Damy Ha on 15-Dec-22.
//

#include "gomea/src/fitness/density.h"

namespace gomea{

/**
 * Constructor.
 * Everything done in the constructor occurs before the start time is set.
 * @param problem_index The index of this problem (used for debugging)
 * @param fitness_function_name The name of this fitness function/ name of the problem
 * @param data The data used by the BDEU function
 * @param max_number_of_parents The maximum number of parents a node can have
 * @param max_number_of_discretizations The maximum number of discretizations per node (with continuous data)
 */
Density::Density(int problem_index,
                 string fitness_function_name,
                 const shared_ptr<DataStructure<double>> &data,
                 size_t max_number_of_parents,
                 size_t max_number_of_discretizations) :
                 Fitness_BN(problem_index, std::move(fitness_function_name), data, max_number_of_parents, max_number_of_discretizations) {
    // All variables are taken care of in the parent constructor
    this->penaltyFactor = 1.0;
    this->fitnessFunctionType = "Density_" + to_string(this->penaltyFactor);

}

Density::~Density() = default;


/**
 * Computes the BDeu fitness of the solution, given the data.
 * @param solution The solution to calculate the BDeu score of.
 * @return The fitness value, fitness per node, constraint value
 */
tuple<double, vector<double>, double> Density::computeFitnessValue(solution_BN &solution) {
    // Perform the discretization
    if (solution.getNumberOfNodesToDiscretize() > 0) {
        this->discretizeData(solution);
    }

    // Initialize variables
    vector<vector<int>> m_prior;                            // array[node_i][node_i_takes_value_k]: Counts the number of occurrences of node_i taking a value of 'k'.
    vector<vector<int>> m_all;                              // array[node_i][number_of_combinations_of_i]: The number of times the parent of node_i take a value combination of 'j'
    vector<int> max_number_of_combinations_per_variable;    // array[node_i]: The (maximum) number of combinations i's parents can make, i.e. the 'columns' in the probability table of node_i
    vector<vector<vector<int>>> combinations;               // combinations[node_i][number_of_combinations_of_i][node_k]: The number of combinations node i's parents make according to the data (example: theoretically we can make 4 combinations with 2 True/False parents. If 'TT' is massing, we get less than 4 combinations).
    vector<vector<vector<int>>> mijk;                       // mijk[node_i][node_i_parents_are_value_j][node_i_is_value_k]: Counts the number of occurrences of node_i taking a value of 'k' when its parents are taking a value (combination) of 'j'.
    vector<double> scores;                                  // array[node_i]: The score of node i

    // Initialize the matrices and vectors
    this->initializeCountMatrices(solution, m_prior, m_all, max_number_of_combinations_per_variable, combinations, mijk);

    // Compute the prior and conditional counts
    this->calculatePriorAndConditionalCounts(solution, m_prior, m_all, max_number_of_combinations_per_variable, combinations, mijk);

    // Compute and assign the BDEU score
    auto result = this->computeDensityScore(solution, mijk, m_all, m_prior, combinations, max_number_of_combinations_per_variable);

    // Save the results (if partial evaluations are used in combination with a valid discretization strategy)
    if (this->usePartialEvaluations) { // and solution.getDiscretizationPolicy()->getSupportsPartialEvaluations()) {
        this->combinedCount = mijk;
        this->countNodePrior = m_prior;
        this->countParentCombinations = m_all;
        this->maxNumberOfCombinationsPerVariable = max_number_of_combinations_per_variable;
        this->parentCombinations = combinations;
        this->scorePerNode = get<1>(result);
    }

    return result;
}


/**
 * Calculates the Density score of the solution.
 * @param solution The solution to calculate the Density score of
 * @param mijk mijk[node_i][node_i_parents_are_value_j][node_i_is_value_k]: Counts the number of occurrences of node_i taking a value of 'k' when its parents are taking a value (combination) of 'j'.
 * @param m_all array[node_i][number_of_combinations_of_i]: The number of times the parent of node_i take a value combination of 'j'
 * @param m_prior array[node_i][node_i_takes_value_k]: Counts the number of occurrences of node_i taking a value of 'k'.
 * @param combinations combinations[node_i][number_of_combinations_of_i][node_k]: The number of combinations node i's parents make according to the data (example: theoretically we can make 4 combinations with 2 True/False parents. If 'TT' is massing, we get less than 4 combinations).
 * @param max_number_of_combinations_per_variable array[node_i]: The (maximum) number of combinations i's parents can make, i.e. the 'columns' in the probability table of node_i
 * @return The fitness value, fitness per node, constraint value
 */
tuple<double, vector<double>, double> Density::computeDensityScore(solution_BN &solution,
                                                                   vector<vector<vector<int>>> &mijk,
                                                                   vector<vector<int>> &m_all,
                                                                   vector<vector<int>> &m_prior,
                                                                   vector<vector<vector<int>>> &combinations,
                                                                   vector<int> &max_number_of_combinations_per_variable) {
    // Determine variables
    vector<size_t> unique_variables_count_pre_discretization = this->originalData->getNumberOfUniqueValues();

    // Calculate the log likelihood over densities
    vector<double> log_likelihood = loglikelihood_density(solution, mijk, m_all, m_prior, max_number_of_combinations_per_variable);

    // Sum the scores (and invert it as we do maximization)
    double density_score = 0;
    vector<double> density_node(this->numberOfNodes);
    for (size_t node_index = 0; node_index < this->numberOfNodes; ++node_index) {
        density_node[node_index] = log_likelihood[node_index]; // - cost_network[node_index] - cost_discretization[node_index];
        density_score += density_node[node_index];
    }

    // RUBEN: revert scores back, since code is written for minimization
    density_score = -density_score;

    // Set the results
    auto result = make_tuple(density_score, density_node, 0);

    return result;
}

/**
 * Calculates the mutual information score
 * @param solution The solution to calculate the mutual information score of
 * @param mijk mijk[node_i][node_i_parents_are_value_j][node_i_is_value_k]: Counts the number of occurrences of node_i taking a value of 'k' when its parents are taking a value (combination) of 'j'.
 * @param m_all array[node_i][number_of_combinations_of_i]: The number of times the parent of node_i take a value combination of 'j'
 * @param m_prior array[node_i][node_i_takes_value_k]: Counts the number of occurrences of node_i taking a value of 'k'.
 * @param combinations combinations[node_i][number_of_combinations_of_i][node_k]: The number of combinations node i's parents make according to the data (example: theoretically we can make 4 combinations with 2 True/False parents. If 'TT' is massing, we get less than 4 combinations).
 * @param max_number_of_combinations_per_variable array[node_i]: The (maximum) number of combinations i's parents can make, i.e. the 'columns' in the probability table of node_i
 * @return
 */
vector<double> Density::loglikelihood_density(solution_BN &solution,
                                              vector<vector<vector<int>>> &mijk,
                                              vector<vector<int>> &m_all,
                                              vector<vector<int>> &m_prior,
                                              vector<int> &max_number_of_combinations_per_variable) {
    // Retrieve variables
    size_t sampleSize = this->data->getNumberOfDataRows();

    // Go over all nodes
    vector<double> result(this->numberOfNodes);
    int continuousCount = 0;
    for(size_t node_index = 0; node_index < this->numberOfNodes; ++node_index) {
        // Determine variables of node
        ColumnDataType nodeType = this->dataTypeNodes[node_index];

        // Complexity score
        int maxParentCombinations = max_number_of_combinations_per_variable[node_index];
        size_t NumberOfinstantiationsOfNode = this->data->getColumnNumberOfClasses()[node_index];
        double complexityScore = (double) maxParentCombinations * (double) (NumberOfinstantiationsOfNode - 1);

        // Go over all parent combination
        double loglikelihood = 0, sumParents = 0;
        for (size_t parentCombinationIndex = 0; parentCombinationIndex < mijk[node_index].size(); ++parentCombinationIndex) {
            sumParents += m_all[node_index][parentCombinationIndex];    // Used for check

            // Calculate the density
            size_t numberOfVariablesOfNode =  mijk[node_index][parentCombinationIndex].size();
            double sumChild = 0;
            for (size_t valueIndex = 0; valueIndex < numberOfVariablesOfNode; ++ valueIndex){
                // Skip the counts that are zero. This does not raise any problems
                size_t count_combination_j = m_all[node_index][parentCombinationIndex];
                size_t count_j_and_k = mijk[node_index][parentCombinationIndex][valueIndex];
                if (count_j_and_k == 0) {   // We don't check 'count_combination_j' as this should be skipped at all times.
                    continue;
                }

                sumChild += mijk[node_index][parentCombinationIndex][valueIndex];   // Used for check

                // Prepare density estimation of p(X | Pa(x))
                double rangeX;
                if (nodeType == Continuous) {
                    // Normalize the boundaries -> RUBEN: shouldn't be necessary, assuming we're using bin width encoding.
                    const shared_ptr<DiscretizationPolicy>& policy = solution.getDiscretizationPolicy();
                    vector<double> boundaries = solution.getBoundaries()[continuousCount];
                    vector<double> dataColumn = this->originalData->getDataMatrix().getColumn(node_index);
                    vector<double> normalizedBoundaries = this->normalizeBoundaries(boundaries, dataColumn);

                    // Determine the range of the boundary
                    if (valueIndex == 0){
                        rangeX = normalizedBoundaries[valueIndex];
                    } else if (valueIndex == numberOfVariablesOfNode - 1) {
                        rangeX = 1 - normalizedBoundaries[valueIndex - 1];
                    } else {
                        rangeX = normalizedBoundaries[valueIndex] - normalizedBoundaries[valueIndex-1];
                    }
                } else {
                    // Discrete
                    rangeX = 1.0 / (double) numberOfVariablesOfNode;
                }

                // Determine the density of the current node and parent node f(X=x and Pa(X)=y)
                // count_j_and_k and count_combination_j can never be zero.
                double p_ijk = (double) count_j_and_k / (double)  count_combination_j;
                double f_x = (1 / rangeX) * p_ijk; // (1 / (double) discretizationsX) *

                loglikelihood += (((double) count_j_and_k / (double) sampleSize) * log(f_x));
            }
            // Perform check
            assert ( sumChild == m_all[node_index][parentCombinationIndex]);
        }

        // Compute the node score
        result[node_index] = ((double) sampleSize*loglikelihood) - this->penaltyFactor * (complexityScore*log((double) sampleSize)/ 2);

        // Perform check
        assert (sumParents == sampleSize);

        // If the node was continuous, increment continuousCount here
        if (nodeType == Continuous) {
            continuousCount++;
        }
    }


    return result;
}

}