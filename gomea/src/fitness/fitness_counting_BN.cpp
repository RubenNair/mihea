//
// Created by Damy Ha on 14-Dec-22.
//

#include "gomea/src/fitness/fitness_BN.h"

namespace gomea{

/**
 * Initializes the matrices and vectors needed to count occurrences. These matrices and vectors are used to calculate the MDL score
 * @param solution The solution that contains general information of the network
 * @param m_prior array[node_i][node_i_takes_value_k]: Counts the number of occurrences of node_i taking a value of 'k'.
 * @param m_all array[node_i][number_of_combinations_of_i]: The number of times the parent of node_i take a value combination of 'j'
 * @param max_number_of_combinations_per_variable array[node_i]: The (maximum) number of combinations i's parents can make, i.e. the 'columns' in the probability table of node_i
 * @param combinations  combinations[node_i][number_of_combinations_of_i][node_k]: The number of combinations node i's parents make according to the data (example: theoretically we can make 4 combinations with 2 True/False parents. If 'TT' is massing, we get less than 4 combinations).
 * @param mijk mijk[node_i][node_i_parents_are_value_j][node_i_is_value_k]: Counts the number of occurrences of node_i taking a value of 'k' when its parents are taking a value (combination) of 'j'. (reused in 'computesamplevalue'), when the parents
 */
void Fitness_BN::initializeCountMatrices(solution_BN &solution,
                                      vector<vector<int>> &m_prior,
                                      vector<vector<int>> &m_all,
                                      vector<int> &max_number_of_combinations_per_variable,
                                      vector<vector<vector<int>>> &combinations,
                                      vector<vector<vector<int>>> &mijk) {
    // Resize the variables
    m_prior.resize(this->numberOfNodes);
    m_all.resize(this->numberOfNodes);
    max_number_of_combinations_per_variable.resize(this->numberOfNodes, 1);
    combinations.resize(this->numberOfNodes);
    mijk.resize(this->numberOfNodes);

    // Initialize the matrices
    for (int node_i = 0; node_i < this->numberOfNodes; ++node_i) {
        // Retrieve variables
        size_t number_of_parents_of_i = solution.getNumberOfParents()[node_i];

        // Set the number of combinations per variable
        // If a node has 1 parent of Low, Med, High, it gets a probability table of 3 columns
        // If a node has 2 parents of Low, Med, High it gets a probability table of 9 columns
        for (int local_parent_index = 0; local_parent_index < number_of_parents_of_i; ++local_parent_index) {
            // Determine the parent's node index
            size_t parent_node = solution.getParentMatrix()[node_i][local_parent_index];

            // Determine the number of combinations per variable
            int classes_of_parent_node = (int) data->getColumnNumberOfClasses()[parent_node];
            max_number_of_combinations_per_variable[node_i] *= classes_of_parent_node;
        }

        // Initialize the 3D: combinations[node_i][number_of_combinations_of_i][number_of_parents_of_i]
        int number_of_combinations_of_i = max_number_of_combinations_per_variable[node_i];
        vector<vector<int>> comb(number_of_combinations_of_i, vector<int>(number_of_parents_of_i, -1));
        combinations[node_i] = comb;

        // Initialize the mijk matrix: mijk[node_i][number_of_combinations_of_i][number_of_classes_of_i]
        size_t number_of_classes_of_i = solution.getDiscretizedData()->getColumnNumberOfClasses()[node_i];
        vector<vector<int>> ijk(number_of_combinations_of_i, vector<int>(number_of_classes_of_i, 0));
        mijk[node_i] = ijk;

        // Initialize the m_all matrix: m_all[node_i][number_of_combinations_of_i]
        vector<int> vec_combinations = vector<int>(number_of_combinations_of_i, 0);
        m_all[node_i] = vec_combinations;

        // Initialize the m_prior matrix: m_prior[node_i][number_of_classes_of_i]
        vector<int> vec_nodes = vector<int>(number_of_classes_of_i, 0);
        m_prior[node_i] = vec_nodes;
    }
}


/**
 * Computes the prior and conditional counts by going over the data
 * @param solution The solution that contains general information of the network
 * @param m_prior array[node_i][node_i_takes_value_k]: Counts the number of occurrences of node_i taking a value of 'k'.
 * @param m_all array[node_i][number_of_combinations_of_i]: The number of times the parent of node_i take a value combination of 'j'
 * @param max_number_of_combinations_per_variable array[node_i]: The (maximum) number of combinations i's parents can make, i.e. the 'columns' in the probability table of node_i
 * @param combinations  combinations[node_i][number_of_combinations_of_i][node_k]: The number of combinations node i's parents make according to the data (example: theoretically we can make 4 combinations with 2 True/False parents. If 'TT' is massing, we get less than 4 combinations).
 * @param mijk mijk[node_i][node_i_parents_are_value_j][node_i_is_value_k]: Counts the number of occurrences of node_i taking a value of 'k' when its parents are taking a value (combination) of 'j'. (reused in 'computesamplevalue'), when the parents
 */
void Fitness_BN::calculatePriorAndConditionalCounts(solution_BN &solution,
                                                 vector<vector<int>> &m_prior,
                                                 vector<vector<int>> &m_all,
                                                 vector<int> &max_number_of_combinations_per_variable,
                                                 vector<vector<vector<int>>> &combinations,
                                                 vector<vector<vector<int>>> &mijk) {
    for (int node_i = 0; node_i < this->numberOfNodes; ++node_i) {
        // Determine variables
        size_t number_of_classes_of_i = solution.getDiscretizedData()->getColumnNumberOfClasses()[node_i];
        int number_of_combinations_of_i = max_number_of_combinations_per_variable[node_i];

        // Calculate the probabilities
        size_t number_of_parents_of_i = solution.getNumberOfParents()[node_i];
        if (number_of_parents_of_i == 0) {
            // Node i does not have parents, no need to calculate the combination of parent variables
            // Count the "prior probabilities" (mijk is not calculated in 'computeSampleValue' because 'i' does not have parents)
            this->computeSampleValue(node_i, solution, max_number_of_combinations_per_variable[node_i], combinations, mijk, m_prior);

            // Assign the prior probabilities
            for (int value_k = 0; value_k < number_of_classes_of_i; ++value_k) {
                // As node_i does not have parents, set the prior probabilities.
                mijk[node_i][0][value_k] = m_prior[node_i][value_k];

                // Determine the number of times, the parents of 'i' take a value of 'j' by counting over all occurrences where 'i' takes a value of 'k'
                m_all[node_i][0] += mijk[node_i][0][value_k];
            }

        } else {
            // Node_i has children thus its combination of parents must be determined
            this->findAllCombinationsValues(node_i, solution, combinations, max_number_of_combinations_per_variable);

            // Count the "Probabilities"
            this->computeSampleValue(node_i, solution, max_number_of_combinations_per_variable[node_i], combinations, mijk, m_prior);

            // Determine the number of times, the parents of 'i' take a value of 'j' by counting over all occurrences where 'i' takes a value of 'k'
            for (int comb = 0; comb < number_of_combinations_of_i; ++comb) {
                for (int value_k = 0; value_k < number_of_classes_of_i; ++value_k) {
                    m_all[node_i][comb] += mijk[node_i][comb][value_k];
                }
            }
        }
    }
}

/**
 * Method to find the actual parent combinations (over theoretical). Assumes that 'combinations' is initialized with '-1'.
 * combinations[node_i][number_of_combinations_of_i][node_k]: The number of combinations node i's parents make according
 * to the data (example: theoretically we can make 4 combinations with 2 True/False parents. If 'TT' is missing in the data,
 * we get less than 4 combinations).
 * 'combinations' is a 3D matrix, accessible by [node_index]. The remaining vector contains the unique combinations of variables.
 * Example: combinations[node_1] shows the unique combinations of variables, e.g.
 * [[0, 0], [0, 1], [1, 0], [-1, -1]], where we could not find [1, 1].
 * @param node_index The index of the node to determine the combinations of
 * @param solution The solution itself, used to get data from
 * @param combinations A pointer to the 'combinations' matrix that needs to be constructed
 * @param number_of_combinations_per_variable The (alleged) number of variable a node can take. This yields an upper bound on the number of combinations
 */
void Fitness_BN::findAllCombinationsValues(int node_index,
                                        const solution_BN &solution,
                                        vector<vector<vector<int>>> &combinations,
                                        vector<int> &number_of_combinations_per_variable) {
    // Initialize variables
    int combination_index = 0;                  // The number of (unique) combinations found from the data
    size_t numberOfParentsOfi = solution.getNumberOfParents()[node_index];  // The number of parents of node i
    vector<int> parents_of_node_i = solution.getParentMatrix()[node_index]; // The parents of node i


    // Continue to find combinations until we have iterated through all possible combinations or until the samples run out
    for (int local_sample_index = 0; local_sample_index < this->getSampleSize(); ++local_sample_index) {
        // Check if the maximum number of combinations have been found, i.e. all combinations have been found.
        if (combination_index >= number_of_combinations_per_variable[node_index]) {
            break;  // Done with the search
        }

        // Initialize variables
        bool combinationAlreadyProcessed = true;    // Used as a flag when a new combination has been found.

        // Retrieve the sample
        vector<double> sample_data = solution.getDiscretizedData()->getDataMatrix().getRow(local_sample_index);

        // Go over all existing combinations that have been seen until now.
        // If a solution contains a variable that has not been seen yet, a new combination might have been detected.
        // If this holds for all existing combinations, i.e. it differs from all existing combinations, we have found a new combination
        for (int i = 0; i < combination_index; ++i) {
            // Assume that we have already processed the combination, until proven otherwise
            vector<int> current_parent_combination = combinations[node_index][i];

            // Compare the existing combination and the possible new combination for every variable (number of parents)
            combinationAlreadyProcessed = sampleMatchesCombination(current_parent_combination,
                                                                   sample_data,
                                                                   parents_of_node_i,
                                                                   (int) numberOfParentsOfi);

            // Check if the combination matches an existing combination
            if (combinationAlreadyProcessed) {
                // As it matches an existing combination, we can stop checking the other existing combinations and
                // go to the next sample.
                break;
            }
        }

        // Check if the combination has been processed yet
        if (!combinationAlreadyProcessed || combination_index == 0) {
            // Add the combination to the combination matrix (by taking the value of each parent)
            size_t number_of_parents_of_node_i = solution.getNumberOfParents()[node_index];
            for (size_t local_parent_index = 0; local_parent_index < number_of_parents_of_node_i; ++local_parent_index) {
                // Retrieve the parent node to copy from
                int parent_node = solution.getParentMatrix()[node_index][local_parent_index];
                int parent_value = (int) sample_data[parent_node];

                combinations[node_index][combination_index][local_parent_index] = parent_value;
            }

            // As a (new) combination has been found, increase the number of combinations found
            combination_index += 1;
        }
    }
}

/**
 * Initializes 'm_prior' and 'mijk'. These vectors count the number of occurrences in the data.
 * @param node_index The index of the node to compute a sample value for
 * @param solution The solution which stores information such as the parent nodes
 * @param maximum_number_of_parent_combinations_of_i The (maximum) number of combinations i's parents can make, i.e. the 'columns' in the probability table of node_i
 * @param combinations A 3D matrix accessible by [node_index]. The remaining vector contains the unique combinations of variables, e.g. for node 2 [[0, 1], [1, 0], [0, 0]]
 * @param mijk Accessible by: [node_i][node_i_parents_are_value_j][node_i_is_value_k]. Counts the number of occurrences of node_i taking a value of 'k' when its parents are taking a value (combination) of 'j'. (reused in 'computesamplevalue'), when the parents
 * @param p_prior Count the number of occurrences in the data where node_i takes value 'k'
 */
void Fitness_BN::computeSampleValue(int node_index,
                                 const solution_BN &solution,
                                 int maximum_number_of_parent_combinations_of_i,
                                 vector<vector<vector<int>>> &combinations,
                                 vector<vector<vector<int>>> &mijk,
                                 vector<vector<int>> &m_prior) {

    // Initialize variables
    size_t number_of_classes_of_i = solution.getDiscretizedData()->getColumnNumberOfClasses()[node_index];

    // Set the prior probabilities by going over the samples
    for (size_t sampleIndex = 0; sampleIndex < this->getSampleSize(); ++sampleIndex) {
        // Retrieve the random sample
        vector<double> sample_data = solution.getDiscretizedData()->getDataMatrix().getRow(sampleIndex);

        // Value of the node
        int sample_value = (int) sample_data[node_index];

        // Increment the class (given the value of the node)
        m_prior[node_index][sample_value]++;
    }

    // Set the parent probabilities for being in configuration j for all nodes i, while node_i is in state 'k'
    // Only do this for nodes with parents (as parentless nodes will use something else)
    vector<size_t> numberOfParents = solution.getNumberOfParents();
    size_t numberOfParentsOfi = numberOfParents[node_index];
    if (numberOfParentsOfi != 0) {
        // Go over all samples
        for (size_t sampleIndex = 0; sampleIndex < this->getSampleSize(); ++sampleIndex) {
            // Retrieve the random sample
            vector<double> sample_data = solution.getDiscretizedData()->getDataMatrix().getRow(sampleIndex);

            // Prepare variables
            vector<int> parents_of_node_i = solution.getParentMatrix()[node_index];

            // Determine the parent combination index, i.e. the parent combination's index in 'combinations'
            int parent_combination_index = -1;
            for (int combination_index_j = 0; combination_index_j < maximum_number_of_parent_combinations_of_i; ++combination_index_j) {
                // Go over all parent combinations and check if e.g. data=[1,2,0] matches every element of the current combination
                // Determine the current parent combination, of which we want to check if it matches the data
                vector<int> current_parent_combination = combinations[node_index][combination_index_j];

                // Check if this parent combination matches all variables of the data
                bool foundParentIndex = sampleMatchesCombination(current_parent_combination,
                                                                 sample_data,
                                                                 parents_of_node_i,
                                                                 (int) numberOfParentsOfi);

                // Check if the combination has been found.
                if (foundParentIndex) {
                    parent_combination_index = combination_index_j;
                    break;
                }

                if (combination_index_j == maximum_number_of_parent_combinations_of_i - 1){
                    assert (parent_combination_index != -1);
                }
            }

            // Determine the class of i
            int value_of_sample_at_node_i = (int) sample_data[node_index];    // Retrieve the value of node i

            // Check for debugging: Combination exists in list of combinations
            assert (parent_combination_index != -1);
            assert (value_of_sample_at_node_i < number_of_classes_of_i);

            // Increment the count
            mijk[node_index][parent_combination_index][value_of_sample_at_node_i]++;
        }
    }
}

////////////////////////
/// Helper functions ///
////////////////////////
bool Fitness_BN::sampleMatchesAnyCombination(vector<vector<int>> combinations,
                                          size_t combinationsFound,
                                          const vector<double>& sample_data,
                                          vector<int> &combinationNodeIndices,
                                          int numberOfNodesInCombination) {
    // Go over all existing combinations that have been seen until now.
    // Check if any existing combination matches the sample
    bool combinationAlreadyProcessed;
    for (int combinationIndex = 0; combinationIndex < combinationsFound; ++combinationIndex) {
        // Assume that we have already processed the combination, until proven otherwise
        vector<int> current_parent_combination = combinations[combinationIndex];

        // Compare the existing combination and the possible new combination for every variable
        combinationAlreadyProcessed = this->sampleMatchesCombination(current_parent_combination,
                                                                     sample_data,
                                                                     combinationNodeIndices,
                                                                     numberOfNodesInCombination);

        // Check if the combination matches an existing combination
        if (combinationAlreadyProcessed) {
            // As it matches an existing combination, we can stop checking the other existing combinations
            return true;
        }
    }

    return false;
}
/**
 * Given a combination (of e.g. parents or spouses) this method tests equality between the combination and the data
 * e.g. if 2 parents produce the combination "23", it will return if 'sample_data' contains "23" at the given columns
 * (combinationNodeIndices) in the data.
 * @param combination The combination which we want to find in the data
 * @param sample_data The data sample
 * @param combinationNodeIndices  The node indices (column indices) of the 'combination' used to retrieve the
 *  combination in the data
 * @param numberOfNodesInCombination The number the node has
 * @return Equality between the data and the parent node
 */
bool Fitness_BN::sampleMatchesCombination(vector<int> combination,
                                       vector<double> sample_data,
                                       vector<int> combinationNodeIndices,
                                       int numberOfNodesInCombination) {
    // Perform test
    assert(combination.size() == numberOfNodesInCombination);

    bool sampleMatchesParentCombination = true;
    for (int localNodeIndex = 0; localNodeIndex < numberOfNodesInCombination; ++localNodeIndex) {
        // Retrieve the value of the current combination to check
        int valueCombination = combination[localNodeIndex];

        // Retrieve the value of the data
        size_t node_index = combinationNodeIndices[localNodeIndex];
        int valueData = (int) sample_data[node_index];

        // Check if the currently selected combination matches with the data
        if (valueData != valueCombination) {
            // Current combination has at least one difference -> Does not match
            sampleMatchesParentCombination = false;
            break;
        }
    }

    return sampleMatchesParentCombination;
}

/**
 * Finds the matching combination index of the nodes to inspect in the sample.
 * @param validCombinations All valid combinations
 * @param sample_data A sample of the data of which we want to find one of the combinations
 * @param combinationNodeIndices The node indices to find the combination of.
 * @return The index of the combination
 */
int Fitness_BN::findCombinationIndex(vector<vector<int>> validCombinations,
                                  const vector<double>& sample_data,
                                  const vector<int>& combinationNodeIndices ) {
    // Go over all combinations and check if the data matches every element of the current combination
    for (int combinationIndex = 0; combinationIndex < validCombinations.size(); ++combinationIndex) {
        // Retrieve the current combination
        vector<int> currentCombination = validCombinations[combinationIndex];

        // Check if this combination matches all variables in the data
        bool combinationMatches = this->sampleMatchesCombination(currentCombination, sample_data, combinationNodeIndices, (int) combinationNodeIndices.size());
        if (combinationMatches) {
            return combinationIndex;
        }
    }

    return -1;
}

}