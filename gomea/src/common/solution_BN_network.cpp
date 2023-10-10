#include "gomea/src/common/solution_BN.hpp"

namespace gomea{

/////////////////////////
/// Network structure ///
/////////////////////////
/**
 * Method to determine information from the nodes of the network.
 * 'output_number_of_parents' contains the number of parents of node 'i'
 * 'output_number_of_children' contains the number of children of node 'i'
 * 'output_parent_matrix' contains the parent nodes (indices) of node 'i' (row). The columns are accessible by: [0, ..., output_number_of_parents[i]].
 * 'output_child_matrix' contains the children (indices) of node 'i' (row). The columns are accessible by: [0, ..., output_number_of_children[i]].
 * Example: output_number_of_children[0] = 2, means that node 0 is the parent of 2 other nodes.
 * Example: output_child_matrix = [[3, 0, 0, 0, 0,], [0, 3, 4, 0, 0, ], ...], given Node 1, 2, ... means
 * that there is a connection from Node 1 to Node 4 (Assuming row 0 is node 1)
 *  and there is a connection from Node 2 to Node 1, 4 and 5
 *  'output_number_of_children' indicates the number of outgoing links per row.
 * @param parametersToProcess The parameter to process. Links can be set to 0 if the number of parents is exceeded.
 * @return A struct with information of the nodes.
 */
NodeInformation solution_BN::findParents(vec_t<int> &parametersToProcess) {
    // Prepare results
    vec_t<size_t> result_number_of_parents(number_of_nodes, 0);
    vec_t<size_t> result_number_of_children(number_of_nodes, 0);
    vec_t<vec_t<int>> result_child_matrix(number_of_nodes);
    vec_t<vec_t<int>> result_parent_matrix(number_of_nodes);

    // Go over all links
    // The loop below tries to go over all combinations of links, e.g. 4 nodes will give A:12, A13, A14, A23, A24, A34
    int node_i = 0, node_j = 0; // See paper where links are defined as going from node i to node j
    for (int i = 0; i < this->number_of_links; ++i) {
        // The scheme excludes link going from node_i to itself
        if (node_i == node_j && node_j < this->number_of_nodes - 1) {
            node_j++;
        } else if (node_i == node_j && node_j == this->number_of_nodes) {
            node_j = node_i + 1;
            node_i++;
        }

        // Process the link value
        if (parametersToProcess[i] == 1) {
            // Process a link from node_i to node_j
            // Check if maximum number of parents has been exceeded for node_i
            if (maximum_number_of_parents > 0 and result_number_of_parents[node_j] >= maximum_number_of_parents) {
                parametersToProcess[i] = 0;
                //dbn.constraint_value=1;
            } else {
                // Set a link in the connection matrix
                // As we know the current number of children/outgoing links (before processing this link),
                // we can use 'result_number_of_children' as an index counter
                result_child_matrix[node_i].push_back(node_j);
                result_parent_matrix[node_j].push_back(node_i);

                // Increment the number of outgoing links
                result_number_of_children[node_i]++;

                // Increment the number of parents
                result_number_of_parents[node_j]++;
            }
        } else if (parametersToProcess[i] == 2) {
            // Process a link from node_j to node_i
            // Check if maximum number of parents has been exceeded for node_j
            if (maximum_number_of_parents > 0 and result_number_of_parents[node_i] >= maximum_number_of_parents) {
                parametersToProcess[i] = 0;
                //dbn.constraint_value=1;
            } else {
                // Set a link in the connection matrix
                result_child_matrix[node_j].push_back(node_i);
                result_parent_matrix[node_i].push_back(node_j);

                // Increment the number of outgoing links
                result_number_of_children[node_j]++;

                // Increment the number of parents
                result_number_of_parents[node_i]++;
            }
        } else {
            assert (parametersToProcess[i] == 0);
        }

        // Increment the node_j until it has visited all higher indexed nodes
        if (node_j < this->number_of_nodes - 1) {
            node_j++;
        } else {
            node_j = node_i + 1;
            node_i++;
        }

        // Just in case
        assert (parametersToProcess[i] == 0 || parametersToProcess[i] == 1 || parametersToProcess[i] == 2);
    }

    // Assign results
    NodeInformation result;
    result.output_number_of_parents = result_number_of_parents;
    result.output_number_of_children = result_number_of_children;
    result.output_parent_matrix = result_parent_matrix;
    result.output_child_matrix = result_child_matrix;

    return result;
}

/**
 * Determines the outgoing links of the adjacency matrix, i.e. some links are missing.
 * The matrix represents a connection between child node_j (row) to parent node_i (column) when it is 1.
 * There is no connection when an entry is 0.
 * @param numberOfChildren A vector containing the number of children node 'i' has.
 * @param parentMatrix A matrix in which row 'i' contains the child node indices of node 'i'
 * @return The adjacency matrix
 */
vec_t<vec_t<int>> solution_BN::findAdjacencyMatrix(vec_t<size_t> numberOfChildren, vec_t<vec_t<int>> parentMatrix) {
    // Initialize results
    vec_t<vec_t<int>> result(this->number_of_nodes, vec_t<int>(this->number_of_nodes, 0));

    // Fill in the adjacency matrix, by going over all nodes
    for(int node_i=0; node_i<this->number_of_nodes; ++node_i) {
        // Go over the children of node_i
        size_t number_of_children_of_i = numberOfChildren[node_i];
        for(int local_index_child=0; local_index_child<number_of_children_of_i; ++local_index_child) {
            int node_j = parentMatrix[node_i][local_index_child];
            result[node_j][node_i] = 1;
        }
    }

    // Assign result
    return result;
}

/**
 * Determines the spouses of the children of the nodes
 */
void solution_BN::determineSpouseNodes(vec_t<size_t> numberOfParents,
                                       vec_t<size_t> numberOfChildren,
                                       vec_t<vec_t<int>> parentMatrix,
                                       vec_t<vec_t<int>> childMatrix) {
    // Initialize result
    vec_t<vec_t<vec_t<size_t>>> result(this->number_of_nodes);   // array[node_i][child_index_j][spouses]

    // Initialize the matrices
    for (size_t nodeIndex = 0; nodeIndex < this->number_of_nodes; ++nodeIndex) {
        // Retrieve variables
        size_t numberOfChildrenNode_i = numberOfChildren[nodeIndex];
        vec_t<int> childrenOfNode_i = childMatrix[nodeIndex];

        // Check if node has children at all, otherwise continue
        if (numberOfChildrenNode_i == 0) { continue; }

        // Go over the children of this node
        vec_t<vec_t<size_t>> spouses(numberOfChildrenNode_i);        // The spouses per child
        for (size_t localChildIndex_j = 0; localChildIndex_j < numberOfChildrenNode_i; ++localChildIndex_j) {
            // Get the spouses of node 'i', with respect to child node 'j'
            int childNodeIndex_j = childrenOfNode_i[localChildIndex_j];             // The node index of child 'j'
            vec_t<int> parentIndices_j = parentMatrix[childNodeIndex_j];    // The parents of node 'j'
            size_t numberOfParentsOf_j = numberOfParents[childNodeIndex_j]; // The number of parents of node 'j'
            vec_t<size_t> spouses_ij = this->determineSpousesOfChildNode(parentIndices_j, numberOfParentsOf_j, (int) nodeIndex);    // We exclude the current node from the list of parents

            // Set the spouses of this child
            spouses[localChildIndex_j] = spouses_ij;
        }

        // Set the result
        result[nodeIndex] = spouses;
    }

    this->spouse_matrix = result;
}

/**
 * Determines the spouses of a node (given a particular child node)
 * @param parentIndicesOfChild The parents of the child node
 * @param numberOfParents The number of parents of the child
 * @param nodeIndexToExclude The parent node index of which the spouses are to be determined
 * @return A list of spouses
 */
vec_t<size_t> solution_BN::determineSpousesOfChildNode(vec_t<int> parentIndicesOfChild,
                                                        size_t numberOfParents,
                                                        size_t nodeIndexToExclude) {
    // As 'parentIndicesOfChild' may (read usually) contains excess zeros (because of the old implementation)
    // we need to manually copy the spouses.
    // Initialize the results
    vec_t<size_t> result;
    result.reserve(parentIndicesOfChild.size() - 1);

    for(size_t localIndex = 0; localIndex < numberOfParents; ++localIndex) {
        // Retrieve the parent index
        size_t parentIndex = (size_t) parentIndicesOfChild[localIndex];

        // Check if the parent index is not the index that needs to be excluded
        if (parentIndex != nodeIndexToExclude) {
            result.push_back(parentIndex);
        }
    }

    assert(result.size() == numberOfParents - 1);   // Check that indeed, one parent is removed
    return result;
}


///////////////////////////
/// Topological sorting ///
///////////////////////////
/**
 * Performs a topological search where links that cause the network to become cyclic are removed.
 * @param parametersToProcess The parameters forming the network. This variable might be changed
 * @param numberOfParents A vector containing the number of parents a node 'i' has.
 * @param adjacencyMatrix The (partial) adjacency matrix of the network
 * // Todo: I think the adjancency matrix should not be processed as it is currently
 */
void solution_BN::topologicalSort(vec_t<int> &parametersToProcess,
                                  vec_t<size_t> numberOfParents,
                                  vec_t<vec_t<int>> adjacencyMatrix) {
    // Prepare variables
    vec_t<int> visited(number_of_nodes, 0);
    vec_t<int> stack(number_of_nodes, 0);
    vec_t<vec_t<int>> temp_adjacency_matrix(std::move(adjacencyMatrix));    // Copy the adjacency matrix

    // Visit nodes without children first as these will not cause a cycle
    for (size_t nodeIndex = 0; nodeIndex < number_of_nodes; ++nodeIndex) {
        if (numberOfParents[nodeIndex] == 0) {
            dfs(nodeIndex, temp_adjacency_matrix, visited, stack, parametersToProcess);
        }
    }

    // Do depth first traversal for all remaining nodes
    for (size_t nodeIndex = 0; nodeIndex < number_of_nodes; ++nodeIndex) {
        // Find non visited nodes
        while (visited[nodeIndex] == 1 && nodeIndex < number_of_nodes - 1) {
            ++nodeIndex;
        }

        // Do the depth first search
        dfs(nodeIndex, temp_adjacency_matrix, visited, stack, parametersToProcess);
    }
}

/**
 * Using the adjacency matrix, a depth first search is executed to detected and remove links that cause
 * the network to become a cyclic network.
 * @param currentNodeIndex The node index to perform the depth first search on.
 * @param tempAdjacentMatrix The temporary adjacency matrix, of which the entries can be changed
 * @param visitedNodes A vector where an element is 1 if node 'i' has been visited
 * @param nodesInStack A vector where an element is 1 if node 'i' is on the stack
 * @param parameters The parameters where links might be removed to break up cycles
 */
void solution_BN::dfs(size_t currentNodeIndex,
                      vec_t<vec_t<int>> &tempAdjacentMatrix,
                      vec_t<int> &visitedNodes,
                      vec_t<int> &nodesInStack,
                      vec_t<int> &preprocessedParameters) {

    // Set node as visited
    visitedNodes[currentNodeIndex] = 1;
    nodesInStack[currentNodeIndex] = 1;

    // Perform a recursive search over all nodes
    for (int u = 0; u < this->number_of_nodes; ++u) {
        // Determine if there is a link form Node U to node I
        int linkFromNodeUToNodeIa = tempAdjacentMatrix[currentNodeIndex][u];

        if (linkFromNodeUToNodeIa == 1 && visitedNodes[u] && nodesInStack[u]) {
            // Link between U and I, but Node U has been visited already and is on the stack
            // This means that there is a cycle: u -> ... -> parent -> currentNodeIndex
            removeEdge((int) currentNodeIndex, u, preprocessedParameters);
        } else if (linkFromNodeUToNodeIa && !visitedNodes[u]) {
            // Link between U and I, but Node U has not been visited yet: Continue search
            dfs(u, tempAdjacentMatrix, visitedNodes, nodesInStack, preprocessedParameters);
        }
    }

    // Reset current node's presence in the stack
    nodesInStack[currentNodeIndex] = 0;
}

/**
 * Delete an edge forming a cycle.
 * @param currentNodeIndex The current index being inspected
 * @param u The last link that completes the cycle
 * @param preprocessedParameters The parameters describing the links
 */
void solution_BN::removeEdge(int currentNodeIndex, int u, vec_t<int> &preprocessedParameters) {
    int pos = 0;
    // Determine which link should be removed
    if (u < currentNodeIndex) {
        // Determine the value in 'preprocessedParameters' related to the link
        for (int i = 0; i < u; i++)
            pos = pos + (this->number_of_nodes - i - 1);
        pos = pos + (currentNodeIndex - u) - 1;

        // Remove the link
        if (preprocessedParameters[pos] == 1) {
            preprocessedParameters[pos] = 0;
//            printf("Remove arc from %d -> %d, position %d\n",currentNodeIndex,u,pos);
        } else {
            // For debugging
            printf("Trying to Remove arc from %d -> %d, position %d but this failed?\n",currentNodeIndex,u,pos);
            assert (preprocessedParameters[pos] == 0 || preprocessedParameters[pos] == 1 || preprocessedParameters[pos] == 2);
        }
    } else {
        // Determine the value in 'preprocessedParameters' related to the link
        for (int i = 0; i < currentNodeIndex; i++)
            pos = pos + (this->number_of_nodes - i - 1);
        pos = pos + (u - currentNodeIndex) - 1;

        // Remove the link
        if (preprocessedParameters[pos] == 2) {
            preprocessedParameters[pos] = 0;
//            printf("Remove arc from %d -> %d, position %d\n",currentNodeIndex,u,pos);
        } else {
            // For debugging
            printf("Trying to Remove arc from %d -> %d, position %d but this failed?\n",currentNodeIndex,u,pos);
            assert (preprocessedParameters[pos] == 0 || preprocessedParameters[pos] == 1 || preprocessedParameters[pos] == 2);
        }
    }
}

////////////////////////
/// Max combinations ///
////////////////////////
/**
 * Determines the max number of parent combinations given the node instantiations
 * @param numberOfInstantiations The instantiations per node (read: max value a node can take)
 * @return The max number of parent combinations that can occur
 */
vec_t<size_t> solution_BN::determineMaxNumberOfParentCombinations(vec_t<size_t> numberOfInstantiations) {
    // Determine the max number of spouse combinations for each node
    vec_t<size_t> result(this->number_of_nodes, 1);
    for (size_t nodeIndex = 0; nodeIndex < this->number_of_nodes; ++nodeIndex) {
        // Set the number of combinations per variable
        // If a node has 1 parent of Low, Med, High, it gets a probability table of 3 columns
        // If a node has 2 parents of Low, Med, High it gets a probability table of 9 columns
        size_t numberOfParents = this->number_of_parents[nodeIndex];
        for (size_t localParentIndex = 0; localParentIndex < numberOfParents; ++localParentIndex) {
            // Determine the parent's node index
            size_t parentNodeIndex = this->parent_matrix[nodeIndex][localParentIndex];

            // Determine the number of combinations per variable
            size_t instantiationsOfParentNode = numberOfInstantiations[parentNodeIndex];
            result[nodeIndex] *= instantiationsOfParentNode;
        }
    }

    return result;
}
/**
 * Determines the maximum number of spouse combination of node i, child j
 * @param numberOfInstantiations The number of instantiations per node
 * @return The number of spouse combinations of node i, child j
 */
vec_t<vec_t<size_t>> solution_BN::determineMaxNumberOfSpouseCombinations(vec_t<size_t> numberOfInstantiations) {
    // Perform check
    assert(!this->spouse_matrix.empty());

    // Determine the number of spouse combinations per node
    vec_t<vec_t<size_t>> result(this->number_of_nodes);      // array[node_i][child_index_j]
    for (size_t nodeIndex = 0; nodeIndex < this->number_of_nodes; ++nodeIndex) {
        // Retrieve variables
        vec_t<vec_t<size_t>> spousesOfChildren = this->spouse_matrix[nodeIndex];
        size_t numberOfChildrenNode_i = spousesOfChildren.size();

        // Check if node has children at all
        if (spousesOfChildren.empty()) { continue; }

        // Go over the (local) children of this node
        vec_t<size_t> maxSpouseCombinationVector(numberOfChildrenNode_i, 1);       // The number of spouse combinations per child node (if equal to 1, there are no spouses. We use this to sum over all child value counts later)
        for (size_t localChildIndex_j = 0; localChildIndex_j < numberOfChildrenNode_i; ++localChildIndex_j) {
            // Determine the maximum number of spouse combinations
            // Get the spouses of node 'i', with respect to child node 'j'
            vec_t<size_t> spouses_ij = this->spouse_matrix[nodeIndex][localChildIndex_j];

            // Go over all spouses to count the maximum number of combination they can make
            for (const auto &spouseNodeIndex : spouses_ij) {
                // Determine the number of combinations per variable
                size_t instantiationsOfSpouseNode = numberOfInstantiations[spouseNodeIndex];
                maxSpouseCombinationVector[localChildIndex_j] *= instantiationsOfSpouseNode;
            }
        }

        // Set the vectors
        result[nodeIndex] = maxSpouseCombinationVector;
    }

    return result;
}

}