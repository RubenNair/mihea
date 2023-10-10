//
// Created by Damy Ha on 07-Nov-22.
//

#ifndef IMPLEMENTATIONS_DISCRETIZATION_H
#define IMPLEMENTATIONS_DISCRETIZATION_H

#include <limits>
#include <numeric>
#include <random>
#include <set>
#include <tuple>

#include "gomea/src/utils/data_structure.h"
#include "gomea/src/utils/data_matrix.h"
// #include "../../../Continuous_Bayesian_Network_GOMEA/include/Random.h"

using namespace std;

/**
 * Used to compactly return information from the discretization process
 */
struct BoundaryInfo {
    vector<double> boundaries;                  // The boundaries determined by the optimization process
    vector<size_t> edges;                       // The edges of the boundaries
    size_t numberOfInstantiations;              // The number of instantiations, of the boundaries
};

struct BoundariesInfo {
    vector<vector<double>> boundaries;          // The boundaries of all nodes
    vector<vector<size_t>> edges;               // The edges of the boundaries
    vector<size_t> numberOfInstantiations;      // The number of instantiations of all nodes
};

struct SpouseCountsIntervalInfo {
    vector<vector<size_t>> combinedCounts;      // [SpouseCombination][ChildInstantiation] -> count
    vector<size_t> countSpouse;                 // [SpouseCombination] -> count
    vector<size_t> countChild;                  // [ChildInstantiation] -> count
    set<size_t> spouseInstantiations;           // All child instantiations in the interval
    set<size_t> childInstantiations;            // All child instantiations in the interval

    // Assignment operator
    SpouseCountsIntervalInfo& operator=(const SpouseCountsIntervalInfo& other) {
        if (this != &other) {
            combinedCounts = other.combinedCounts;
            countSpouse = other.countSpouse;
            countChild = other.countChild;
            spouseInstantiations = other.spouseInstantiations;
            childInstantiations = other.childInstantiations;
        }
        return *this;
    }
};


class DiscretizationPolicy {
public:
    // Constructor
    DiscretizationPolicy(vector<ColumnDataType> dataTypePerNode);
    DiscretizationPolicy(vector<ColumnDataType> dataTypePerNode, vector<size_t> desiredNumberOfInstantiations);
    ~DiscretizationPolicy();

    // Constructor for cloning
    DiscretizationPolicy(string policyName,
                         vector<ColumnDataType> dataTypePerNode,
                         vector<size_t> desiredNumberOfInstantiations,
                         vector<size_t> numberOfInstantiations,
                         vector<vector<double>> boundaries,
                         bool needsNetworkStructure);
    DiscretizationPolicy(string policyName,
                         vector<size_t> desiredNumberOfInstantiations,
                         vector<size_t> numberOfInstantiations,
                         vector<vector<double>> boundaries,
                         bool needsNetworkStructure,
                         vector<size_t> numberOfParents,
                         vector<size_t> numberOfChildren,
                         vector<vector<int>> childMatrix,
                         vector<vector<int>> parentMatrix,
                         vector<vector<vector<size_t>>> spouseMatrix
                         );
    virtual shared_ptr<DiscretizationPolicy> clone();

    // Additional initialization
    void initializeNetworkStructure(vector<size_t> newNumberOfParents,
                                    vector<size_t> newNumberOfChildren,
                                    vector<vector<int>> newChildMatrix,
                                    vector<vector<int>> newParentMatrix,
                                    vector<vector<vector<size_t>>> newSpouseMatrix);

    // Methods
    // Determines the discretization boundaries
    virtual void determineDiscretizationBoundaries(DataStructure<double> &data, const shared_ptr<mt19937>& rng);
    // Determines the discretization boundaries of the given nodesToRediscretize
    virtual void redetermineSelectBoundaries(DataStructure<double> &data, const vector<size_t> &nodesToRediscretize, vector<vector<double>> boundariesToKeep, const shared_ptr<mt19937>& rng);
    // Discretizes the data of a single node
    vector<double> discretizeData(size_t nodeIndex, const vector<double>& dataOfNode);
    // Discretizes the data of all nodes
    vector<vector<double>> discretizeData(DataStructure<double> &data);

    // Getters
    const string &getPolicyName() const;
    const vector<ColumnDataType> &getDataTypePerNode() const;
    const vector<size_t> &getDesiredNumberOfInstantiations() const;
    const vector<size_t> &getNumberOfInstantiations() const;
    const vector<vector<double>> &getBoundaries() const;
    bool getNeedsNumberOfInstantiationsInAdvance() const;
    bool getNeedsNetworkStructure() const;
    bool getSupportsPartialEvaluations() const;

protected:
    // Variables
    string policyName;                              // The name of the discretization policy
    vector<ColumnDataType> dataTypePerNode;         // The data type per node
    vector<size_t> desiredNumberOfInstantiations;   // The desired number of discretizations (to perform) per node
    vector<size_t> numberOfInstantiations;          // The number of instantiations determined by the discretization policy
    vector<vector<double>> boundaries;              // The boundaries of the discretization method (for a given node)

    // Instantiations: Some policies need to know the number of instantiations, others determine this value on their own
    bool needsNumberOfInstantiationsInAdvance;      // The policy needs the number of discretizations before doing discretization from the solution
    bool supportsPartialEvaluations;                // The policy supports partial evaluations

    // Network structure: Some policies need to know more about the network structure
    bool needsNetworkStructure;                     // The policy needs the underlying network structure
    vector<size_t> numberOfParents;                 // The number of parents of each node
    vector<size_t> numberOfChildren;                // The number of children of each node
    vector<vector<int>> childMatrix;                // Matrix indicating outgoing links (per row)
    vector<vector<int>> parentMatrix;               // Matrix indicating ingoing links (per row)
    vector<vector<vector<size_t>>> spouseMatrix;    // The spouse of each child of node i (can be empty)

    /// Helper functions
    // Determines the continuous node order in a reverse topology
    vector<size_t> determineReverseTopologyOfContinuousNodes();
    // Determine the continuous node order randomly
    vector<size_t> determineRandomTopologyOfContinuousNodes(const shared_ptr<mt19937>& rng);
    // Determines the most instantiations of all nodes (returns 0 when there are no discrete nodes)
    size_t determineMostInstantiationsDiscreteVariables(DataStructure<double> &data);
    // Determines the lambda among the markov blanket
    size_t determineLambdaValue(vector<vector<size_t>> data);
    // Determines the permutation of a sort
    vector<size_t> determinePermutationOfSort(vector<double> continuousData);
    // Finds the index of a given combination
    int findIndexOfCombination(vector<vector<size_t>> &combinations, vector<size_t> &combinationToCheck);
    // Computes the unique parent combinations from parent data
    tuple<vector<size_t>, vector<vector<size_t>>> mapRowDataToCombinationIndex(const DataMatrix<size_t>& rowData);
    // Computes the head indices
    vector<size_t> determineHeads(const vector<double>& sortedContinuousData, size_t numberOfUniqueValues);
    // Computes the tail indices
    vector<size_t> determineTails(const vector<double>& sortedContinuousData, size_t numberOfUniqueValues);
    // Converts unique edges to boundaries
    vector<double> convertEdgesToBoundaries(const vector<size_t>& uniqueEdges, const vector<size_t>& heads, const vector<size_t>& tails, const vector<double>& sortedContinuousData, size_t policyNumberOfInstantiations);

    /// Discretization
    // Discretizes continuous data given edges
    vector<size_t> discretizeContinuousData(vector<double> data, const vector<double>& edges);

    /// Equal width discretization
    // Performs an initial equal width discretization on continuous nodes
    tuple<DataMatrix<size_t>, vector<size_t>> applyEqualWidthOnData(DataStructure<double> data, size_t numberOfBins);
    vector<vector<size_t>> discretizeAllDataOld(DataStructure<double> data, size_t numberOfBins);
    // Divides the data of a node into equal width sized bins
    vector<size_t> equalWidthDiscretization(vector<double> nodeData, size_t numberOfBins);
    // Determines the edges of the equal width discretization method given sorted data
    BoundaryInfo determineEqualDiscretizationEdges(vector<double> sortedContinuousDataNode_i, size_t numberOfBins);

    /// Sorting
    // Sort the continuous data of a node and discretized data based on a specific node's data
    tuple<vector<double>, DataMatrix<size_t>> sortContinuousAndDiscreteData(size_t nodeIndex, DataStructure<double> &continuousData, DataMatrix<size_t> &discretizedData);
};

#endif //IMPLEMENTATIONS_DISCRETIZATION_H
