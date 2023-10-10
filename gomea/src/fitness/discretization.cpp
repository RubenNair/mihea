//
// Created by Damy Ha on 07-Nov-22.
//

#include <utility>

#include "gomea/src/fitness/discretization.h"

/**
 * Constructor
 * @param dataTypePerNode The data type per node
 */
DiscretizationPolicy::DiscretizationPolicy(vector<ColumnDataType> dataTypePerNode) {
    // Set variables
    this->policyName = "ParentPolicy";
    this->dataTypePerNode = std::move(dataTypePerNode);
    this->needsNetworkStructure = false;
    this->needsNumberOfInstantiationsInAdvance = false;
}

/**
 * Constructor
 * @param dataTypePerNode The data type per node
 * @param desiredNumberOfInstantiations The requested number of instantiations per continuous node
 */
DiscretizationPolicy::DiscretizationPolicy(vector<ColumnDataType> dataTypePerNode, vector<size_t> desiredNumberOfInstantiations) {
    this->policyName = "ParentPolicy";
    this->dataTypePerNode = std::move(dataTypePerNode);
    this->desiredNumberOfInstantiations = std::move(desiredNumberOfInstantiations);
    this->needsNetworkStructure = false;
}

DiscretizationPolicy::~DiscretizationPolicy() = default;

/**
 * Initializes the network structure for discretization policies that need it
 * @param newNumberOfParents The number of parent per node
 * @param newNumberOfChildren The number of children per node
 * @param newChildMatrix The indices of the children per node
 * @param newParentMatrix The indices of the parents per node
 * @param newSpouseMatrix The spouses of each child node
 */
void DiscretizationPolicy::initializeNetworkStructure(vector<size_t> newNumberOfParents,
                                                      vector<size_t> newNumberOfChildren,
                                                      vector<vector<int>> newChildMatrix,
                                                      vector<vector<int>> newParentMatrix,
                                                      vector<vector<vector<size_t>>> newSpouseMatrix) {
    // Set the variables
    this->numberOfParents = std::move(newNumberOfParents);
    this->numberOfChildren = std::move(newNumberOfChildren);
    this->childMatrix = std::move(newChildMatrix);
    this->parentMatrix = std::move(newParentMatrix);
    this->spouseMatrix = std::move(newSpouseMatrix);
}

/**
 * Determines the discretization boundaries of the data for all nodes
 * @param data The data to discretize
 * @param rng The random number generator
 */
void DiscretizationPolicy::determineDiscretizationBoundaries(DataStructure<double> &data, const shared_ptr<mt19937>& rng) {
    throw runtime_error("Discretization method not implemented.");
}

/**
 * Determines the discretization boundaries of the data for the given nodesToRediscretize
 * @param data The data to discretize
 * @param nodesToRediscretize The nodesToRediscretize to determine the boundaries of
 * @param boundariesToKeep The boundaries that should not be redetermined
 * @param rng The random number generator
 */
void DiscretizationPolicy::redetermineSelectBoundaries(DataStructure<double> &data,
                                                       const vector<size_t> &nodesToRediscretize,
                                                       vector<vector<double>> boundariesToKeep,
                                                       const shared_ptr<mt19937>& rng) {
    throw runtime_error("Discretization method for given nodesToRediscretize not implemented.");
}

/**
 * Discretizes the data, if the boundaries have been determined
 * @param nodeIndex The index of the node to discretize (i.e. used to select the policy)
 * @param dataOfNode The data of the node to discretize
 * @return The discretized data
 */
vector<double> DiscretizationPolicy::discretizeData(size_t nodeIndex, const vector<double>& dataOfNode) {
    // Check that the boundaries have been determined
    assert(!this->boundaries.empty());                      // Check that the boundaries have been initialized
    assert(!(this->boundaries.size() < nodeIndex));         // Check that this node has potential boundaries
    assert(!this->boundaries[nodeIndex].empty());           // Check that this node has boundaries
    assert(this->dataTypePerNode[nodeIndex] != Discrete);   // Check that this node is not discrete

    // Initialize result
    vector<double> nodeBoundaries = this->boundaries[nodeIndex];
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
vector<vector<double>> DiscretizationPolicy::discretizeData(DataStructure<double> &data) {
    // Perform check
    assert(!this->boundaries.empty());
//    assert(!this->boundaryEdges.empty());
    assert(!this->numberOfInstantiations.empty());

    // Initialize results
    size_t numberOfRows = data.getNumberOfDataRows();
    size_t numberOfColumns = data.getNumberOfDataColumns();
    DataMatrix<double> result(numberOfRows, numberOfColumns);

    // Discretize all nodes
    for(size_t nodeIndex = 0; nodeIndex < numberOfColumns; ++nodeIndex) {
        ColumnDataType dataType = data.getColumnType()[nodeIndex];
        if (dataType == Continuous) {
            // Discretize data
            vector<double> continuousData = data.getDataMatrix().getColumn(nodeIndex);
            vector<double> discretizedData = this->discretizeData(nodeIndex, continuousData);

            result.setColumn(nodeIndex, discretizedData);
        } else {
            // Add the already discrete data
            vector<double> targetData = data.getDataMatrix().getColumn(nodeIndex);
//            vector<size_t> discreteData = convertToType<double, size_t>(targetData);
            result.setColumn(nodeIndex, targetData);
        }
    }

    return result.getRawMatrix();
}

/**
 * Discretizes continuous data
 * @param data The data to discretize
 * @param edges The edges
 * @return The discretized data
 */
vector<size_t> DiscretizationPolicy::discretizeContinuousData(vector<double> data, const vector<double>& edges) {
    // Initialize result
    size_t sampleSize = data.size();
    vector<size_t> result(sampleSize);

    // Go over each sample value
    for(size_t sampleIndex = 0; sampleIndex < sampleSize; ++sampleIndex) {
        // Determine the sample value
        double sampleValue = data[sampleIndex];

        // Find the correct bin
        size_t binNumber = edges.size();    // If a single boundary does not match, it's the last instantiation
        for (int binIndex = 0; binIndex < edges.size(); ++binIndex) {
            if (sampleValue <= edges[binIndex]) {    // Increment, as we have one starting edge
                binNumber = binIndex;
                break;
            }
        }

        result[sampleIndex] = binNumber;
    }

    return result;
}

// Getters
const string &DiscretizationPolicy::getPolicyName() const { return policyName; }
const vector<ColumnDataType> &DiscretizationPolicy::getDataTypePerNode() const { return dataTypePerNode; }
const vector<size_t> &DiscretizationPolicy::getDesiredNumberOfInstantiations() const { return desiredNumberOfInstantiations; }
const vector<size_t> &DiscretizationPolicy::getNumberOfInstantiations() const { return numberOfInstantiations; }
const vector<vector<double>> &DiscretizationPolicy::getBoundaries() const { return boundaries; }
bool DiscretizationPolicy::getNeedsNumberOfInstantiationsInAdvance() const { return needsNumberOfInstantiationsInAdvance; }
bool DiscretizationPolicy::getNeedsNetworkStructure() const { return needsNetworkStructure; }
bool DiscretizationPolicy::getSupportsPartialEvaluations() const { return supportsPartialEvaluations; }

/**
 * Constructor for cloning
 */
DiscretizationPolicy::DiscretizationPolicy(string policyName,
                                           vector<size_t> desiredNumberOfInstantiations, vector<size_t> numberOfInstantiations,
                                           vector<vector<double>> boundaries,
                                           bool needsNetworkStructure,
                                           vector<size_t> numberOfParents, vector<size_t> numberOfChildren,
                                           vector<vector<int>> childMatrix, vector<vector<int>> parentMatrix,
                                           vector<vector<vector<size_t>>> spouseMatrix) {
    this->desiredNumberOfInstantiations = std::move(desiredNumberOfInstantiations);
    this->numberOfInstantiations = std::move(numberOfInstantiations);
    this->boundaries = std::move(boundaries);
    this->needsNetworkStructure = needsNetworkStructure;
    this->numberOfParents = std::move(numberOfParents);
    this->numberOfChildren = std::move(numberOfChildren);
    this->childMatrix = std::move(childMatrix);
    this->parentMatrix = std::move(parentMatrix);
    this->spouseMatrix = std::move(spouseMatrix);
}

/**
 * Constructor for cloning
 */
DiscretizationPolicy::DiscretizationPolicy(string policyName,
                                           vector<ColumnDataType> dataTypePerNode,
                                           vector<size_t> desiredNumberOfInstantiations,
                                           vector<size_t> numberOfInstantiations,
                                           vector<vector<double>> boundaries,
                                           bool needsNetworkStructure) {
    this->policyName = std::move(policyName);
    this->dataTypePerNode = std::move(dataTypePerNode);
    this->desiredNumberOfInstantiations = std::move(desiredNumberOfInstantiations);
    this->numberOfInstantiations = std::move(numberOfInstantiations);
    this->boundaries = std::move(boundaries);
    this->needsNetworkStructure = needsNetworkStructure;
}
/**
 * Cloning
 * @return An error, as each policy should implement its own cloning method.
 */
shared_ptr<DiscretizationPolicy> DiscretizationPolicy::clone() {
    throw runtime_error("Cloning method of discretization policy not implemented.");
}


