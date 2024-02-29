//
// Created by Damy Ha on 27-Sep-22.
//

#include "gomea/src/fitness/fitness_BN.h"

namespace gomea{
/**
 * Constructor fitness function
 * @param problem_index The problem index (used for debugging)
 * @param fitness_function_name The fitness function name
 * @param data The data used in the fitness funciton
 * @param max_number_of_parents The maximum number of parents allowed
 * @param max_number_of_discretizations The maximum number of discretizations
 */
Fitness_BN::Fitness_BN(int problem_index,
                 string fitness_function_name,
                 const shared_ptr<DataStructure<double>> &data,
                 size_t max_number_of_parents,
                 size_t max_number_of_discretizations) {
    // Set variables
    this->problemIndex = problem_index;
    this->fitnessFunctionName = std::move(fitness_function_name);
    this->maxNumberOfParents = max_number_of_parents;
    this->maxNumberOfDiscretizations = max_number_of_discretizations;

    // Initialize variables
    this->optimalNetworkAvailable = false;
    this->numberOfFullEvaluations = 0;
    this->numberOfEvaluations = 0;
    this->usePartialEvaluations = false;
    this->postProcessing = true;

    // Initialize variables dependent on the data
    if (data != nullptr) {
        // Determine variables
        this->originalData = data;
        this->data = data->clone();
        this->numberOfNodes = data->getNumberOfDataColumns();
        this->numberOfLinks = calculateNumberOfLinks(this->numberOfNodes);
        this->numberOfNodesToDiscretize = data->getNumberOfContinuousVariables();
        this->dataTypeNodes = data->getColumnType();

        // Determine the indices of the nodes to discretize
        this->separateNodeType();

        // Check that the maximum number of discretizations is valid
        if (this->numberOfNodesToDiscretize != 0 && max_number_of_discretizations < 1 ) {
            throw runtime_error("The maximum number of discretizations (for continuous variables) must be larger than 1.");
        }

        this->loglikelihood = numeric_limits<double>::max();
    }
}

Fitness_BN::~Fitness_BN() = default;

/**
 * Determines the indices of the nodes that need to be discretized
 */
void Fitness_BN::separateNodeType() {
    this->continuousNodeIndices = this->data->determineContinuousNodeIndices();
    this->discreteNodeIndices = this->data->determineDiscreteNodeIndices();
}

/**
 * Computes the fitness value
 * @param solution
 */
tuple<double, vector<double>, double> Fitness_BN::computeFitnessValue(solution_BN &solution) {
    throw runtime_error("Fitness function not implemented.");
}


/**
 * Discretizes the data based on the solution's discretization policy
 * @param solution The solution (contains the discretization policy)
 */
void Fitness_BN::discretizeData(solution_BN &solution) {
    // Actually, discretize based on continuous variable values.
    vector<vector<double>> discretizedData = discretizeData(solution, *this->originalData);
    vector<size_t> instantiationsPerNode = getNumberOfInstantiations(solution);

    shared_ptr<DataStructure<double>> result = make_shared<DataStructure<double>>(discretizedData,
                                                                                  this->originalData->getColumnNames(),
                                                                                  this->originalData->getColumnType(),  // NOTE: Officially this should be discrete data.
                                                                                  instantiationsPerNode,
                                                                                  this->originalData->getNumberOfUniqueValues());
    solution.setDiscretizedData(result);
}

vector<size_t> Fitness_BN::getNumberOfInstantiations(solution_BN &solution)
{
    vector<size_t> result;
    int continuousIndex = 0;
    for (size_t nodeIndex = 0; nodeIndex < solution.getNumberOfNodes(); ++nodeIndex) {
        if (solution.getNodeDataTypes()[nodeIndex] == Discrete) {
            result.push_back(solution.getNumberOfDiscretizationsperNode()[nodeIndex]);
        } else {
            result.push_back(solution.getBoundaries()[continuousIndex].size() + 1);
            continuousIndex++;
        }
    }

    return result;
}

/**
 * Discretizes the data, if the boundaries have been determined
 * @param nodeIndex The index of the node to discretize (i.e. used to select the policy)
 * @param dataOfNode The data of the node to discretize
 * @return The discretized data
 */
vector<double> Fitness_BN::discretizeData(solution_BN &solution, size_t nodeIndex, int continuousIndex, const vector<double>& dataOfNode) {
    // Check that the boundaries have been determined
    assert(!solution.getBoundaries().empty());                      // Check that the boundaries have been initialized
    assert(!(solution.getBoundaries().size() < continuousIndex));         // Check that this node has potential boundaries
    assert(!solution.getBoundaries()[continuousIndex].empty());           // Check that this node has boundaries
    assert(solution.getNodeDataTypes()[nodeIndex] != Discrete);   // Check that this node is not discrete

    // Initialize result
    vector<double> nodeBoundaries = solution.getBoundaries()[continuousIndex];
    vector<double> result(dataOfNode.size());

    // Go over each data item
    for (size_t dataIndex = 0; dataIndex < dataOfNode.size(); ++dataIndex) {
        // Load the data
        double dataValue = dataOfNode[dataIndex];

        // Try out each boundary
        size_t clusterValue = nodeBoundaries.size();       // Set to the last element, as the value will not be checked
        for (int boundaryIndex = 0; boundaryIndex < nodeBoundaries.size(); ++boundaryIndex) {
            // Check if the boundary is larger than the value
            if (dataValue <= nodeBoundaries[boundaryIndex]) {
                clusterValue = boundaryIndex;
                break;
            }
        }

        result[dataIndex] = (double) clusterValue;
    }

    return result;
}

/**
 * Discretize all data
 * @param data The data to discretize
 * @return All data discretized
 */
vector<vector<double>> Fitness_BN::discretizeData(solution_BN &solution, DataStructure<double> &data) {
    // Perform check
    assert(!solution.c_variables.empty());
    assert(!getInitialNumberOfInstantiations().empty());


    // Initialize results
    size_t numberOfRows = data.getNumberOfDataRows();
    size_t numberOfColumns = data.getNumberOfDataColumns();
    DataMatrix<double> result(numberOfRows, numberOfColumns);
    int continuousIndex = 0;

    // Discretize all nodes
    for(size_t nodeIndex = 0; nodeIndex < numberOfColumns; ++nodeIndex) {
        ColumnDataType dataType = data.getColumnType()[nodeIndex];
        if (dataType == Continuous) {
            // Discretize data
            vector<double> continuousData = data.getDataMatrix().getColumn(nodeIndex);
            vector<double> discretizedData = this->discretizeData(solution, nodeIndex, continuousIndex, continuousData);

            result.setColumn(nodeIndex, discretizedData);
            continuousIndex++;
        } else {
            // Add the already discrete data
            vector<double> targetData = data.getDataMatrix().getColumn(nodeIndex);
            result.setColumn(nodeIndex, targetData);
        }
    }

    return result.getRawMatrix();
}



/**
 * Rediscretize particular nodes
 * @param solution The new solution of which we use the discretization policy
 * @param nodesToRediscretize A vector of nodes to rediscretize
 */
void Fitness_BN::rediscretizeNodes(solution_BN &solution, const vector<size_t>& nodesToRediscretize) {
    cout << "Used rediscretizeNodes function, likely needs to be reworked first" << endl;
    exit(0);
}


/**
 * Computes the fitness of a solution.
 * Executes 'computeFitnessValue()' and takes care of the other variables
 * @param solution The solution to calculate the fitness of
 */
void Fitness_BN::computeFitness(solution_BN &solution) {
    // Computes the fitness
    double score, constraintValue;
    vector<double> nodeScore;
    tie(score, nodeScore, constraintValue) = this->computeFitnessValue(solution);

    // Increment number of evaluations and set the previous solution
    this->numberOfEvaluations++;
    this->numberOfFullEvaluations++;

    // Sets the statistics
    clock_t currentTime = clock();

    solution.setFitness(score);
    solution.setConstraintValue(constraintValue);
    solution.setNumberOfEvaluations(this->numberOfEvaluations);
    solution.setNumberOfFullEvaluations(this->numberOfFullEvaluations);
    solution.setTimeStamp(currentTime);

    if (this->usePartialEvaluations) {
        shared_ptr<solution_BN> prevSolPtr(solution.clone());
        this->previousSolution = prevSolPtr;
    }
}

/**
 * Computes the fitness of a solution without counting it as an evaluation
 * Please be careful with this method
 * @param solution The solution to calculate the fitness of
 */
void Fitness_BN::computeFitnessWithoutCountingEvaluations(solution_BN &solution) {
    // Compute the fitness
    double score, constraintValue;
    vector<double> nodeScore;
    tie(score, nodeScore, constraintValue) = this->computeFitnessValue(solution);

    // Sets the statistics
    solution.setFitness(score);
    solution.setConstraintValue(constraintValue);
    solution.setNumberOfEvaluations(this->numberOfEvaluations);
    solution.setNumberOfFullEvaluations(this->numberOfFullEvaluations);
}

/**
 * Reset the data (used to compute the BDeu score) to the original data
 */
void Fitness_BN::resetToOriginalData(){
    this->data = this->originalData->clone();
}

/**
 * Creates the original solution, if available.
 * @return The optimal solution
 */
shared_ptr<solution_BN> Fitness_BN::getGoldenSolution() {
    // Check if the golden solution is available
    if (this->optimalNetworkAvailable) {
        // Check if the golden solution has already been created
        if (!this->optimalSolution) {
            cout << "Golden solution not yet implemented/adapted" << endl;
            exit(0);
        }
        shared_ptr<solution_BN> result(this->optimalSolution->clone());
        return result;
    } else {
        return nullptr;
    }
}

// Getters
bool Fitness_BN::hasOptimalNetwork() const { return optimalNetworkAvailable; }
bool Fitness_BN::hasOptimalBoundariesAvailable() const { return optimalBoundariesAvailable; }
const vector<int> &Fitness_BN::getOptimalNetwork() const { return optimalNetwork; }
const vector<vector<double>> &Fitness_BN::getOptimalBoundaries() const { return optimalBoundaries; }
size_t Fitness_BN::getNumberOfNodes() const { return numberOfNodes; }
size_t Fitness_BN::getNumberOfLinks() const { return numberOfLinks; }
size_t Fitness_BN::getMaximumNumberOfParents() const { return maxNumberOfParents; }
size_t Fitness_BN::getNumberOfFullEvaluations() const { return numberOfFullEvaluations; }
size_t Fitness_BN::getNumberOfEvaluations() const { return numberOfEvaluations; }
const string &Fitness_BN::getFitnessFunctionName() const { return fitnessFunctionName; }
const string &Fitness_BN::getFitnessFunctionType() const { return fitnessFunctionType;}
size_t Fitness_BN::getNumberOfNodesToDiscretize() const { return numberOfNodesToDiscretize; }
size_t Fitness_BN::getMaxNumberOfDiscretizations() const { return maxNumberOfDiscretizations; }
const shared_ptr<DataStructure<double>> &Fitness_BN::getOriginalData() const { return originalData; }
const shared_ptr<DataStructure<double>> &Fitness_BN::getData() const { return data; }
const vector<ColumnDataType> &Fitness_BN::getDataTypeNodes() const { return dataTypeNodes; }
size_t Fitness_BN::getSampleSize() { return this->data->getNumberOfDataRows(); }
bool Fitness_BN::getUsePartialEvaluations() const { return usePartialEvaluations; }
const vector<size_t> &Fitness_BN::getContinuousNodeIndices() const { return continuousNodeIndices; }
const vector<size_t> &Fitness_BN::getDiscreteNodeIndices() const { return discreteNodeIndices; }
const vector<size_t> &Fitness_BN::getInitialNumberOfInstantiations() const { return originalData->getColumnNumberOfClasses(); }
const string &Fitness_BN::getFitnessFunctionBaseName() const { return fitnessFunctionBaseName; }
bool Fitness_BN::getPostProcessing() const { return postProcessing; }
const shared_ptr<mt19937> &Fitness_BN::getRng() const { return rng; }

// Setters
void Fitness_BN::setData(const shared_ptr<DataStructure<double>> &data) { Fitness_BN::data = data; }
void Fitness_BN::setUsePartialEvaluations(bool usePartialEvaluations) { Fitness_BN::usePartialEvaluations = usePartialEvaluations; }
void Fitness_BN::setFitnessFunctionBaseName(const string &fitnessFunctionBaseName) { Fitness_BN::fitnessFunctionBaseName = fitnessFunctionBaseName; }
void Fitness_BN::setPreviousSolution(const shared_ptr<solution_BN> &previousSolution) { Fitness_BN::previousSolution = previousSolution; }
void Fitness_BN::setRng(const shared_ptr<mt19937> &rng) { Fitness_BN::rng = rng; }
void Fitness_BN::setOptimalNetwork(const vector<int> &optimalSolution) {
    Fitness_BN::optimalNetwork = optimalSolution;
    Fitness_BN::optimalNetworkAvailable = true;
}
void Fitness_BN::setOptimalBoundaries(const vector<vector<double>> &optimalBoundaries) {
    Fitness_BN::optimalBoundaries = optimalBoundaries;
    Fitness_BN::optimalBoundariesAvailable = true;
}
void Fitness_BN::setPostProcessing(bool postProcessing) {
    if (postProcessing) {
        if (this->fitnessFunctionBaseName.empty()) {
            cout << "Post-Processing required, however the baseName is not set." << endl;
            Fitness_BN::postProcessing = false;
        } else {
            Fitness_BN::postProcessing = true;
        }
    } else {
        Fitness_BN::postProcessing = false;
    }
}

/// Printing
void Fitness_BN::printOptimalNetwork() {
    if (optimalNetworkAvailable) {
        for (const auto &value : optimalNetwork) {
            cout << (int) value;
        }
        cout << endl;
    }
}

void Fitness_BN::printOptimalBoundaries() {
    if (optimalBoundariesAvailable) {
        for (size_t nodeIndex = 0; nodeIndex < numberOfNodes; ++nodeIndex) {
            if (this->dataTypeNodes[nodeIndex] == Discrete) {
                cout << ";";
                continue;
            }

            for(size_t boundaryIndex = 0; boundaryIndex < optimalBoundaries[nodeIndex].size(); ++boundaryIndex) {
                cout << optimalBoundaries[nodeIndex][boundaryIndex];

                if (boundaryIndex != optimalBoundaries[nodeIndex].size() - 1) {
                    cout << ", ";
                }
            }

            if (nodeIndex != numberOfNodes - 1) {
                cout << ";";
            }
        }
        cout << endl;
    }
}

/**
 * Remaps a set of boundaries to [0, 1] -> also makes sure that bin widths add up to 1.
 * @param boundaries The boundaries that need to be remapped
 * @param continuousData The continuous data (of which this function needs the min and max value)
 * @return The normalized boundaries
 */
vector<double> Fitness_BN::normalizeBoundaries(vector<double> boundaries, vector<double> continuousData) {
    // Retrieve the min and max values
    double minDataValue = *min_element(continuousData.begin(), continuousData.end());
    double maxDataValue = *max_element(continuousData.begin(), continuousData.end());

    // Todo: Quick hack when data does not span entire range
    double minValue = min(minDataValue, boundaries[0]);
    double maxValue = max(maxDataValue, boundaries.back());

    // Remap every boundary
    vector<double> result(boundaries.size());

    double total = 0.0;
    for (size_t boundaryIndex = 0; boundaryIndex < boundaries.size(); ++boundaryIndex) {
        // Retrieve the current value
        double currentBoundaryValue = boundaries[boundaryIndex];
        
        // Remap to [0, 1]
        result[boundaryIndex] = (1/(maxValue - minValue)) * (currentBoundaryValue - minValue);

        // Perform check
        assert(result[boundaryIndex] >= 0);
        assert(result[boundaryIndex] <= 1);
    }

    return result;
}

}
