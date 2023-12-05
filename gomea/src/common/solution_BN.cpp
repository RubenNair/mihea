// WIP: Create solution class for BN problems by modifying code below, based on solution_so.cpp from other codebase.

#include <cmath>

#include "gomea/src/common/solution_BN.hpp"

#include "gomea/src/fitness/so_benchmarks.h"


using namespace std;


namespace gomea{



// solution_BN::solution_BN( int number_of_variables_, size_t alphabetSize_, int number_of_c_variables_ ) : solution_t(number_of_variables_, alphabetSize_), c_variables(vec_t<double>(number_of_c_variables_)) {}

// solution_BN::solution_BN( vec_t<int> &variables, vec_t<double> &c_variables ) : solution_mixed(variables), c_variables(vec_t<double>(c_variables)) {}

solution_BN::solution_BN( size_t numberOfVariables_, size_t alphabetSize_, size_t numberOfCVariables_, fitness_t<int> *problemInstance_ ) : solution_mixed(numberOfVariables_, alphabetSize_, numberOfCVariables_)
{
	// solution_t(numberOfVariables_, alphabetSize_);
    fill(c_variables.begin(), c_variables.end(), 0);
    this->problemInstance = problemInstance_;
}


solution_BN::solution_BN(size_t numberOfVariables_, 
                         size_t alphabetSize_, 
                         size_t numberOfCVariables_,
                        //  vec_t<int> &unprocessedParameters,
                         vec_t<ColumnDataType> node_data_types,
                         const vec_t<size_t> &initialNumberOfInstantiations,
                         int discretizationPolicyIndex,
                         size_t maximum_number_of_parents,
                         size_t maximum_number_of_instantiations,
                         fitness_t<int> *problemInstance_,
                         vec_t<double> maxValuesData, vec_t<double> minValuesData,
                         double lower_user_range, double upper_user_range,
                         shared_ptr<DataStructure<double>> data,
                         double populationIndexRatio,
                         bool useNormalizedCVars,
                         bool transformCVariables,
                         bool useOptimalSolution,
                         bool guaranteedInitSpread,
                         bool extraCVarForNumberOfBins,
                         string problemInstancePath,
                         int runIndex) : solution_BN(numberOfVariables_, alphabetSize_, numberOfCVariables_, problemInstance_)
{
	// Set parameters
    this->node_data_types = node_data_types;
    this->maximum_number_of_parents = (int) maximum_number_of_parents;

    this->maxValuesData = maxValuesData;
    this->minValuesData = minValuesData;

    this->lower_user_range = lower_user_range;
    this->upper_user_range = upper_user_range;

    this->useNormalizedCVars = useNormalizedCVars;
    this->transformCVariables = transformCVariables;
    this->useOptimalSolution = useOptimalSolution;
    this->guaranteedInitSpread = guaranteedInitSpread;
    this->extraCVarForNumberOfBins = extraCVarForNumberOfBins;
    this->problemInstancePath = problemInstancePath;
    this->runIndex = runIndex;
    

    if(data != NULL) 
    {
        this->data = data;
    }

    // Determine parameters
    this->number_of_nodes = (int) node_data_types.size();
    this->number_of_links = calculateNumberOfLinks(this->number_of_nodes);
    this->number_of_nodes_to_discretize = std::count(node_data_types.begin(), node_data_types.end(), Continuous);

    // Initialize parameters
    this->fitness = numeric_limits<double>::min();              // Set to the lowest double
    this->constraintValue = numeric_limits<double>::max();      // Set to the highest value

    // Perform check
//    assert(number_of_links + number_of_nodes_to_discretize == unprocessedParameters.size());    // Check that all parameters are used

    // Initialize variables
    if(useOptimalSolution)
    {
        optimalInit();
    } else if(populationIndexRatio >= 0)
    {
        randomInit(&gomea::utils::rng, populationIndexRatio);
    } else
    {
        randomInit(&gomea::utils::rng);
    }

    if(transformCVariables)
    {
        execTransformationCVariables();
    }
    
    // if using normalized c_vars, do normalization after potential transformation
    //  (so it is done on the transformed values, similar to how it will be called from iamalgam)
    if(useNormalizedCVars)
    {
        normalize();
    }

    // Process the solution
    NetworkStructure solutionInformation = processParametersSolution(variables, initialNumberOfInstantiations, discretizationPolicyIndex, maximum_number_of_instantiations);

    // Assign the results
    this->variables.resize(solutionInformation.processedParameters.size());
    // this->variables = solutionInformation.processedParameters;
    transform(solutionInformation.processedParameters.begin(), solutionInformation.processedParameters.end(), this->variables.begin(), [](int i){return i;});
    this->number_of_parents = solutionInformation.nodeInfo.output_number_of_parents;
    this->number_of_children = solutionInformation.nodeInfo.output_number_of_children;
    this->parent_matrix = solutionInformation.nodeInfo.output_parent_matrix;
    this->child_matrix = solutionInformation.nodeInfo.output_child_matrix;
    this->spouse_matrix.resize(0);  // Not set as not all algorithms need this
    this->adjacency_matrix = solutionInformation.adjacencyMatrix;
    this->numberOfDiscretizationsperNode = solutionInformation.discretizationPerNode;
    // this->discretizationPolicy = solutionInformation.discretizationPolicy;

    // Store the boundaries
    updateBoundaries();
}

void solution_BN::optimalInit()
{
    assert(problemInstancePath != "");
    string pathOptimalSolution = determinePathOptimalSolution(problemInstancePath, runIndex);

    string stringOptimalSolution;
    vec_t<vec_t<double>> optimalBoundaries;
    retrieveOptimalSolution(pathOptimalSolution, stringOptimalSolution, optimalBoundaries);

    // Convert the string to a vector of ints
    for(int i = 0; i < getNumberOfVariables(); i++)
    {
        variables[i] = stringOptimalSolution[i] - '0';
    }

    // Store the optimalBoundaries in the boundaries variable. Skip rows that have no elements.
    vec_t<vec_t<double>> boundaries(this->number_of_nodes_to_discretize);
    int cVarsCount = 0;
    for(int i = 0; i < optimalBoundaries.size(); i++)
    {
        if(optimalBoundaries[i].size() > 0)
        {
            for(int j = 0; j < optimalBoundaries[i].size(); j++)
            {
                boundaries[cVarsCount].push_back(optimalBoundaries[i][j]);
            }
            cVarsCount++;
        }
    }
    this->boundaries = boundaries;
}

// Only used in situations like initializing offspringpopulation, where the exact values don't matter (just that the arrays are initialized)
void solution_BN::randomInit(std::mt19937 *rng)
{
	for (int i = 0; i < getNumberOfVariables(); ++i)
	{
		variables[i] = (*rng)() % getAlphabetSize();
	}

    // c_variables = {0.2000153015, 0.2000154391, 0.2000154391, 0.2000154391, 0.1999383812, 0, 0, 0};

    // variables = {1, 2, 0, 2, 0, 1};
    // c_variables = {0.1987034297, 0.2003387483, 0.2003387483, 0.2003387483, 0.2002803253, 0, 0, 0, 0,
    //                 0.3343046869, 0.3344275367, 0.3312677764, 0, 0, 0, 0, 0, 0,
    //                 0.4999268203, 0.5000731797, 0, 0, 0, 0, 0, 0, 0};

    for (int i = 0; i < getNumberOfCVariables(); ++i) 
    {
        c_variables[i] += lower_user_range + ((*rng)() / (double)(*rng).max()) * (upper_user_range - lower_user_range); //(*rng)() / (double)(*rng).max() * 0.02 - 0.01; //
		// // TODO RUBEN: hardcoded upper and lower bounds for now, should be read from config maybe?
		// c_variables[i] = -10 + ((*rng)() / (double)(*rng).max()) * 20; 
    }

    if(useNormalizedCVars)
    {
        normalize();
    }
}

/**
 * Normalizes just the c_variables, the first numberOfBins variables for each continuous node.
 * If numberOfBins is -1, all c_variables are normalized.
*/
void solution_BN::normalize(int numberOfBins) 
{
    // If we are using the optimal solution, don't normalize (just to make sure they will always match conditions of the optimum)
    if(useOptimalSolution)
    {
        return;
    }

    // Normalize in [0, 1] range, so if variables were transformed, temporarily transform them back
    // For now, if numberOfBins is not -1, it means we are in init stage of transformedCVariables - guaranteedSpread combo (method 2 init 3)
    // In which case the variables are not yet transformed, and therefore they shouldn't be transformed in this function (will be handled after)
    if(transformCVariables && numberOfBins == -1)
    {
        execTransformationCVariables();
    }

    int maxDiscretizations;
    int offset = 0;
    if(extraCVarForNumberOfBins)
    {
        maxDiscretizations = (this->c_variables.size() / this->number_of_nodes_to_discretize) - 1;
        offset = this->number_of_nodes_to_discretize;
    } else
    {
        maxDiscretizations = this->c_variables.size() / this->number_of_nodes_to_discretize;
    }
    
    int cVarsCount = 0;
    for(int j = 0; j < this->number_of_nodes; j++)
    {
        if(node_data_types[j] == Continuous)
        {
            double sum = 0.0;
            for(int i = 0; i < maxDiscretizations; i++)
            {
                if(numberOfBins != -1 && i >= numberOfBins)
                {
                    break;
                }
                int c_var_index = offset + maxDiscretizations * cVarsCount + i;
                sum += c_variables[c_var_index];
            }
            for(int i = 0; i < maxDiscretizations; i++)
            {
                int c_var_index = offset + maxDiscretizations * cVarsCount + i;
                c_variables[c_var_index] /= sum;
            }
            cVarsCount++;
        }
    }
    // See previous comment about transformation in this function
    if(transformCVariables && numberOfBins == -1)
    {
        execTransformationCVariables();
    }
}

void solution_BN::execTransformationCVariables()
{
    if(transformCVariables)
    {
        for(int i = 0; i < getNumberOfCVariables(); i++)
        {
            if(extraCVarForNumberOfBins && i < this->number_of_nodes_to_discretize)
            {
                continue;
            }
            
            c_variables[i] = 1.0 / c_variables[i];
        }
    }
}

/**
 * random initialization function of the solution that scales the init of c_variables with the index of the solution in the population.
*/
void solution_BN::randomInit(std::mt19937 *rng, double populationIndexRatio)
{
    if(problemInstancePath == "network7")
    {
        // This network has the network structure fixed in the Yi-Chun Chen paper, so fix it here as well.
        variables = {0, 1, 1};
    } else
    {
        for (int i = 0; i < getNumberOfVariables(); ++i)
        {
            variables[i] = (*rng)() % getAlphabetSize();
        }
    }

    int maxDiscretizations = getNumberOfCVariables() / this->number_of_nodes_to_discretize;
    double minUpperRangeBound = 1.0 / maxDiscretizations;
    
    // The number of boundaries per node should be between 2 and maxDiscretizations (at least 2 boundaries, since after normalization the latter will be 1 and disregarded)
    long numberOfBoundariesPerNode = fmin(std::lround(populationIndexRatio * (maxDiscretizations - 1) + 1.5), maxDiscretizations);

    for (int i = 0; i < getNumberOfCVariables(); ++i) 
    {
        if(extraCVarForNumberOfBins)
        {
           // For each continuous node, there is an extra c-variable at the start of the array to indicate the amount of bins for that node.
           // First initialize these variables to the correct value (i.e. the same for all solutions)
              if(i < this->number_of_nodes_to_discretize)
              {
                c_variables[i] = numberOfBoundariesPerNode;
                continue;
              }
        }
        if(guaranteedInitSpread)
        {
            c_variables[i] = lower_user_range + ((*rng)() / (double)(*rng).max()) * (upper_user_range - lower_user_range);
        } else
        {
            // Scale upper bound of c_variables based on the normalized index of the solution in the population
            double newUpperBound = minUpperRangeBound + populationIndexRatio * (upper_user_range - minUpperRangeBound);
            assert(newUpperBound >= lower_user_range);
            c_variables[i] = lower_user_range + ((*rng)() / (double)(*rng).max()) * (newUpperBound - lower_user_range);
            // // TODO RUBEN: hardcoded upper and lower bounds for now, should be read from config maybe?
            // c_variables[i] = -10 + ((*rng)() / (double)(*rng).max()) * 20; 
        }
    }

    if(guaranteedInitSpread)
    {
        normalize(numberOfBoundariesPerNode);
    }
}

/**
 * generate the boundary values based on the continuous variables and data.
*/
void solution_BN::updateBoundaries()
{
    // If we are using the optimal solution, don't update the boundaries (just to make sure they will always be the same as optimum)
    if(useOptimalSolution)
    {
        return;
    }

    // If c_variables are transformed, transform them back to get them in [0, 1] range for calculating the boundaries
    if(transformCVariables)
    {
        execTransformationCVariables();
    }

    if(guaranteedInitSpread)
    {
        updateBoundariesBasedOnNumberOfDataSamples();
    } else
    {
        updateBoundariesBasedOnBinWidths();
    }

    // Transform back if necessary
    if(transformCVariables)
    {
        execTransformationCVariables();
    }
}

void solution_BN::updateBoundariesBasedOnNumberOfDataSamples()
{
    const int data_size = data->getNumberOfDataRows();
    double step_size = 1.0 / data_size;
    int cVarsCount = 0;
    int offset = 0;
    int maxDiscretizations;
    int loopCondition = maxDiscretizations - 1;
    if(extraCVarForNumberOfBins)
    {
        maxDiscretizations = (this->c_variables.size() / this->number_of_nodes_to_discretize) - 1;
        offset = this->number_of_nodes_to_discretize;
    }
    else 
    {
        maxDiscretizations = this->c_variables.size() / this->number_of_nodes_to_discretize;
    }
    
    vec_t<vec_t<double>> boundaries(this->number_of_nodes_to_discretize);
    for(int i = 0; i < this->number_of_nodes; i++)
    {
        if(node_data_types[i] == Continuous)
        {
            double binwidthUsed = 0.0;
            int loopCondition = extraCVarForNumberOfBins ? std::lround(c_variables[cVarsCount]) : loopCondition;

            for(int j = 0; j < loopCondition; j++)
            {
                int c_var_index = offset + maxDiscretizations * cVarsCount + j;
                // If c_variables[c_var_index] is 0, skip (because there shouldn't be a bin with 0 width)
                if(c_variables[c_var_index] == 0) {
                    continue;
                }
                double boundary;
                if(binwidthUsed + c_variables[c_var_index] - 1.0 >= -1e-5) { // Check if adding this binwidth would cause the boundary to exceed 1.0 (rounding to 5-digit precision). Could give issues if stepsize < 1e-5 (i.e. more than 100.000 samples)!
                    break;
                }
                else {
                    vec_t<int> sorted_indices = data->getSortedDataIndices()[i];
                    DataMatrix<double> data_matrix = data->getDataMatrix();
                    int data_samples_up_to_curr_boundary = (binwidthUsed + c_variables[c_var_index]) / step_size;
                    // Only add a boundary if there will be at least one data sample in the new bin
                    if(data_samples_up_to_curr_boundary != (int)(binwidthUsed / step_size) && data_samples_up_to_curr_boundary + 1 < data_size)
                    {
                        boundary = (data_matrix.getElement(sorted_indices[data_samples_up_to_curr_boundary], i) + data_matrix.getElement(sorted_indices[data_samples_up_to_curr_boundary + 1], i)) / 2.0;
                    } else
                    {  
                        // Boundary would create a bin without datapoints in it, or be after the last sample in data, so ignore it instead (still adding to binwidthUsed)? (TODO: figure out if that is correct interpretation)
                        binwidthUsed += c_variables[c_var_index];            
                        continue;
                    }
                } 
                boundaries[cVarsCount].push_back(boundary);
                binwidthUsed += c_variables[c_var_index];
            }
            cVarsCount++;
        }
    }

    this->boundaries = boundaries;
}

void solution_BN::updateBoundariesBasedOnBinWidths()
{
    int cVarsCount = 0;
    vec_t<vec_t<double>> boundaries(this->number_of_nodes_to_discretize);
    int offset = 0;
    int maxDiscretizations;
    int loopCondition = maxDiscretizations - 1;
    if(extraCVarForNumberOfBins)
    {
        maxDiscretizations = (this->c_variables.size() / this->number_of_nodes_to_discretize) - 1;
        offset = this->number_of_nodes_to_discretize;
    }
    else 
    {
        maxDiscretizations = this->c_variables.size() / this->number_of_nodes_to_discretize;
    }
    
    for(int i = 0; i < this->number_of_nodes; i++)
    {
        if(node_data_types[i] == Continuous)
        {
            double binwidthUsed = 0.0;
            double range = maxValuesData[i] - minValuesData[i];
            int loopCondition = extraCVarForNumberOfBins ? std::lround(c_variables[cVarsCount]) : loopCondition;

            for(int j = 0; j < loopCondition; j++)
            {
                int c_var_index = offset + maxDiscretizations * cVarsCount + j;
                // If c_variables[c_var_index] is 0, skip (because there shouldn't be a bin with 0 width)
                if(c_variables[c_var_index] == 0) {
                    continue;
                }
                double boundary;
                if(binwidthUsed + c_variables[c_var_index] - 1.0 >= -1e-5) { // Check if adding this binwidth would cause the boundary to exceed 1.0 (rounding to 5-digit precision)
                    break;
                }
                else {
                    boundary = binwidthUsed + c_variables[c_var_index];
                } 
                boundaries[cVarsCount].push_back(minValuesData[i] + boundary * range);
                binwidthUsed += c_variables[c_var_index];
            }
            cVarsCount++;
        }
    }

    this->boundaries = boundaries;
}


int solution_BN::getNumberOfCVariables() const
{
	return c_variables.size();
}

void solution_BN::insertCVariables( const vec_t<double> &vars_to_insert )
{
	assert( vars_to_insert.size() == c_variables.size() );
	for( size_t i = 0; i < c_variables.size(); i++ )
	{
		c_variables[i] = vars_to_insert[i];
	}
}

void solution_BN::insertCVariables( vec_t<double> vars_to_insert, vec_t<int> indices_to_insert )
{
	for( size_t i = 0; i < indices_to_insert.size(); i++ )
	{
		int ind = indices_to_insert[i];
		c_variables[ind] = vars_to_insert[i];
	}
}


void solution_BN::insertSolution( solution_mixed *solution )
{
    solution_BN *casted_solution = (solution_BN*) solution;
    assert(casted_solution != nullptr); // TODO
	insertVariables(casted_solution->variables);
    insertCVariables(casted_solution->c_variables);
	setObjectiveValues( casted_solution->getObjectiveValues() );
	setConstraintValue( casted_solution->getConstraintValue() );
	setFitnessBuffers( casted_solution->fitness_buffers );

    assert(0); // TODO
    // // Store the boundaries
    // updateBoundaries();
}

void solution_BN::print()
{
	for( size_t i = 0; i < variables.size(); i++ )
		printf("%c ", variables[i]);
    for( size_t i = 0; i < c_variables.size(); i++ )
	{
        printf("%6.5f ",(double)c_variables[i]);
	}
	printf("\n");
}


/**
 * Reinitialize the network variables after saving memory
 */
void solution_BN::reinitializeNetworkVariables() {
    NodeInformation nodeInformation = this->findParents(this->variables);

    this->number_of_parents = nodeInformation.output_number_of_parents;
    this->number_of_children = nodeInformation.output_number_of_children;
    this->child_matrix = nodeInformation.output_child_matrix;
    this->parent_matrix = nodeInformation.output_parent_matrix;
    this->adjacency_matrix = this->findAdjacencyMatrix(this->number_of_children, this->child_matrix);
}

/**
 * Delete the network variables to save memory
 */
void solution_BN::clearNetworkVariables() {
    // Reset the bayesian network structure parameters
    this->number_of_parents.resize(0);
    this->number_of_children.resize(0);
    this->child_matrix.resize(0);
    this->parent_matrix.resize(0);
    this->adjacency_matrix.resize(0);
    this->spouse_matrix.resize(0);   // Matrix indicating the spouses of each child node of node i
}

/**
 * Processes the raw solution.
 * First excess parents, i.e. links that cause nodes to have too many parents are removed.
 * Then links that cause the network to have cycles are removed.
 * The network structure's information: Parents and children per node, adjacency matrix are simultaneously determined.
 * After this, the number of discretizations per node are determined, together with the initialization of the discretization policy itself
 * @param parametersToProcess The solution as an int array
 * @param initialNumberOfInstantiations The initial number of instantiations (for discrete nodes)
 * @param discretizationPolicyIndex The discretization policy index
 * @param maximum_number_of_instantiations The maximum number of instantiations per node
 * @return A struct containing the cleaned up solution, resulting vectors and matrices.
 */
NetworkStructure solution_BN::processParametersSolution(const vec_t<int> &parametersToProcess,
                                                        const vec_t<size_t> &initialNumberOfInstantiations,
                                                        int discretizationPolicyIndex,
                                                        size_t maximum_number_of_instantiations) {
    // Process the parameters related to the network structure optimization
    NetworkStructure network = processNetworkStructure(parametersToProcess);

    // Process the discretization parameters
    network = processDiscretization(network, initialNumberOfInstantiations, discretizationPolicyIndex, maximum_number_of_instantiations);

    return network;
}

/**
 * Method to recalculate network structure (mostly to enforce the maximum number of parents and remove loops).
 * Should be called whenever `variables` is changed.
*/
void solution_BN::reProcessParametersSolution(vector<int> newParameters)
{
    // Process the solution
    NetworkStructure solutionInformation = processNetworkStructure(newParameters);

    // Assign the results
    this->variables.resize(solutionInformation.processedParameters.size());
    // this->variables = solutionInformation.processedParameters;
    transform(solutionInformation.processedParameters.begin(), solutionInformation.processedParameters.end(), this->variables.begin(), [](int i){return i;});
    this->number_of_parents = solutionInformation.nodeInfo.output_number_of_parents;
    this->number_of_children = solutionInformation.nodeInfo.output_number_of_children;
    this->parent_matrix = solutionInformation.nodeInfo.output_parent_matrix;
    this->child_matrix = solutionInformation.nodeInfo.output_child_matrix;
    this->spouse_matrix.resize(0);  // Not set as not all algorithms need this
    this->adjacency_matrix = solutionInformation.adjacencyMatrix;
    // this->numberOfDiscretizationsperNode = solutionInformation.discretizationPerNode;
    // this->discretizationPolicy = solutionInformation.discretizationPolicy;
}

/**
 * Processes the parameters related to the bayesian network structure learning.
 * The maximum number of parents per node is checked and enforced in this method.
 * Loops in the network structure are also removed.
 * @param parametersToProcess The raw solution in vector form.
 * @return The network parameters:
 *  - The solution without exceeding number of parents and loops
 *  - The number of parents and children
 *  - The indices of the parents and children per node
 *  - The adjacency matrix
 */
NetworkStructure solution_BN::processNetworkStructure(vec_t<int> parametersToProcess) {
    // Preprocess the parameters by removing links that cause children to have too many parents
    // and to calculate the number of children/parents each node has + the connection matrices (used in the next steps)
    NodeInformation preProcessedNodeInformation = this->findParents(parametersToProcess);

    // Determine (part of) the adjacency matrix of the preprocessed solution
    // We use the children over the parents, as removing the cyclic links is easier to do using the children.
    vec_t<vec_t<int>> preProcessedAdjacencyMatrix = this->findAdjacencyMatrix(
            preProcessedNodeInformation.output_number_of_children,
            preProcessedNodeInformation.output_child_matrix);

    // Remove cyclic links in the parameters
    // We use the children, as going from parent to child is easier than child to parents.
    this->topologicalSort(parametersToProcess,
                          preProcessedNodeInformation.output_number_of_children,
                          preProcessedAdjacencyMatrix);

    // Recalculate the parents/children + connection matrix of the sanitized solution (without cycles/too many parents)
    NodeInformation nodeInformation = this->findParents(parametersToProcess);

    // Recalculate the adjacency matrix
    vec_t<vec_t<int>> adjacencyMatrix = this->findAdjacencyMatrix(nodeInformation.output_number_of_parents,
                                                                    nodeInformation.output_parent_matrix);

    // Partly fill in the results
    NetworkStructure result;
    result.processedParameters = parametersToProcess;
    result.nodeInfo = nodeInformation;
    result.adjacencyMatrix = adjacencyMatrix;

    return result;
}

/**
 * Determines the number of discretizations per node and initializes the discretization policies.
 * @param result The Network structure after removing loops and limiting the number of parents
 * @param initialNumberOfInstantiations The initial number of instantiations (for discrete nodes)
 * @param discretizationPolicyIndex The discretization policy index
 * @param maximum_number_of_instantiations The maximum number of discretizations
 * @return The network structure including the number of discretizations per node.
 */
NetworkStructure solution_BN::processDiscretization(NetworkStructure result,
                                                    vec_t<size_t> initialNumberOfInstantiations,
                                                    int discretizationPolicyIndex,
                                                    size_t maximum_number_of_instantiations) {

    // For now, no discretization policy (since we will learn discretizations with iAMaLGaM), so just set discretizationsPerNode to parameter value and then return.
    result.discretizationPerNode = initialNumberOfInstantiations; // NOTE: only values for discrete nodes are correct and should be used.
    return result;

    // // Generate an empty policy
    // shared_ptr<DiscretizationPolicy> policySettings = generateDiscretizationPolicyForSettings(discretizationPolicyIndex);

    // // Retrieve the number of instantiations from the solution.
    // // Todo: Make sure that the number of instantiations are set
    // size_t rawSolutionPointerIndex = this->number_of_links - 1; // Keeps track which parameters from the raw solution have been processed
    // vec_t<size_t> desiredNumberOfInstantiationsPerNode;
    // if (policySettings->getNeedsNumberOfInstantiationsInAdvance()) {
    //     desiredNumberOfInstantiationsPerNode = determineNumberOfDiscretizationsFromSolution(result, rawSolutionPointerIndex, initialNumberOfInstantiations);
    // }

    // // Initialize the discretization policy objects
    // shared_ptr<DiscretizationPolicy> policy = generateDiscretizationPolicy(discretizationPolicyIndex, this->node_data_types, desiredNumberOfInstantiationsPerNode);

    // // Add the network to the policy if it needs the policy
    // if (policy->getNeedsNetworkStructure()) {
    //     // Determine the spouse nodes
    //     NodeInformation nodeInformation = result.nodeInfo;
    //     vec_t<size_t> numberOfParents = nodeInformation.output_number_of_parents;
    //     vec_t<size_t> numberOfChildren = nodeInformation.output_number_of_children;
    //     vec_t<vec_t<int>> parentMatrix = nodeInformation.output_parent_matrix;
    //     vec_t<vec_t<int>> childMatrix = nodeInformation.output_child_matrix;
    //     this->determineSpouseNodes(numberOfParents, numberOfChildren, parentMatrix, childMatrix);
    //     // Pass the network structure
    //     policy->initializeNetworkStructure(numberOfParents, numberOfChildren, childMatrix, parentMatrix, this->spouse_matrix);
    // }


    // // Append the results
    // result.discretizationPolicy = policy;
    // result.discretizationPerNode = policy->getNumberOfInstantiations();

    // return result;
}

/**
 * Determines the number of discretizations per node from the raw solution
 * @param networkStructure The network structure that contains the (cleaned up) solution
 * @param rawSolutionPointerIndex The index to the parameter in the solution that has been last processed
 * @param initialNumberOfInstantiations The number of instantiations of the discrete nodes
 * @return The number of instantiations of the continuous nodes (set by the solution)
 */
vec_t<size_t> solution_BN::determineNumberOfDiscretizationsFromSolution(const NetworkStructure &networkStructure,
                                                                         size_t &rawSolutionPointerIndex,
                                                                         const vec_t<size_t> &initialNumberOfInstantiations) {
    // Retrieve the number of discretizations for every continuous node
    vec_t<size_t> result(this->number_of_nodes);
    for(size_t nodeIndex = 0; nodeIndex < this->number_of_nodes; ++nodeIndex) {
        // Determine the data type of the node
        ColumnDataType node_type = this->node_data_types[nodeIndex];
        if (node_type == Continuous) {
            // Increment the index in the parameter vector
            rawSolutionPointerIndex++;

            // Retrieve the number of discretizations
            vec_t<int> parameterVector = networkStructure.processedParameters;
            auto numberOfDiscretizationsInSolution = (size_t) parameterVector[rawSolutionPointerIndex];

            // Set the number of discretizations
            result[nodeIndex] = numberOfDiscretizationsInSolution;
        } else {
            // Use the original (discrete) value
            result[nodeIndex] = initialNumberOfInstantiations[nodeIndex];
        }
    }

    return result;
}


///////////////
/// Getters ///
///////////////
vec_t<int> solution_BN::getParameters() const { return variables; }
int solution_BN::getNumberOfNodes() const { return number_of_nodes; }
size_t solution_BN::getNumberOfLinks() const { return number_of_links; }
size_t solution_BN::getNumberOfNodesToDiscretize() const { return number_of_nodes_to_discretize; }
double solution_BN::getFitness() const { return fitness;}
double solution_BN::getConstraintValue() const { return constraintValue; }
int solution_BN::getMaximumNumberOfParents() const { return maximum_number_of_parents; }
const vec_t<size_t> &solution_BN::getNumberOfParents() const { return number_of_parents; }
const vec_t<size_t> &solution_BN::getNumberOfChildren() const { return number_of_children; }
const vec_t<vec_t<int>> &solution_BN::getParentMatrix() const { return parent_matrix; }
const vec_t<vec_t<int>> &solution_BN::getChildMatrix() const { return child_matrix; }
const vec_t<vec_t<vec_t<size_t>>> &solution_BN::getSpouseMatrix() const { return spouse_matrix; }
const vec_t<vec_t<int>> &solution_BN::getAdjacencyMatrix() const { return adjacency_matrix; }
const vec_t<ColumnDataType> &solution_BN::getNodeDataTypes() const { return node_data_types; }
const vec_t<size_t> &solution_BN::getNumberOfDiscretizationsperNode() const { return numberOfDiscretizationsperNode; }
const shared_ptr<DiscretizationPolicy> &solution_BN::getDiscretizationPolicy() const { return discretizationPolicy; }
vec_t<int> solution_BN::getNetworkParameters() const { vec_t<int> result(variables.begin(), variables.begin() + number_of_links); return result; }
clock_t solution_BN::getTimeStamp() const { return timeStamp; }
size_t solution_BN::getNumberOfFullEvaluations() const { return numberOfFullEvaluations; }
double solution_BN::getNumberOfEvaluations() const { return numberOfEvaluations; }
const vec_t<vec_t<double>> &solution_BN::getBoundaries() const { return boundaries; }
const vec_t<double> &solution_BN::getMaxValuesData() const { return maxValuesData; }
const vec_t<double> &solution_BN::getMinValuesData() const { return minValuesData; }
const shared_ptr<DataStructure<double>> &solution_BN::getDiscretizedData() const { return discretizedData; }



///////////////
/// Setters ///
///////////////
void solution_BN::setDiscretizedData(const shared_ptr<DataStructure<double>> &data) { solution_BN::discretizedData = data; }
void solution_BN::setFitness(double newFitnessValue) { solution_BN::fitness = newFitnessValue; }
void solution_BN::setConstraintValue(double newConstraintValue) { solution_BN::constraintValue = newConstraintValue; }
void solution_BN::setTimeStamp(clock_t timeStamp) { solution_BN::timeStamp = timeStamp; }
void solution_BN::setNumberOfFullEvaluations(size_t numberOfFullEvaluations) { solution_BN::numberOfFullEvaluations = numberOfFullEvaluations; }
void solution_BN::setNumberOfEvaluations(double numberOfEvaluations) { solution_BN::numberOfEvaluations = numberOfEvaluations; }
void solution_BN::setDiscretizationPolicy(const shared_ptr<DiscretizationPolicy> &discretizationPolicy) { solution_BN::discretizationPolicy = discretizationPolicy; }

/////////////////
/// Operators ///
////////////////
/**
 * Compares the fitness and constraint value of this solution and 'rhs'
 * @param rhs The other solution to compare
 * @return true if this solution is better than 'rhs' otherwise false
 */
bool solution_BN::thisSolutionIsBetterThan(const solution_BN &rhs) {
    // This solution is infeasible
    if (constraintValue > 0) {
        if (rhs.constraintValue > 0) {
            // Both are infeasible, the best constraint value is better
            return  constraintValue < rhs.constraintValue;
        } else {
            // This solution is infeasible, but rhs is not
            return false;
        }
    } else {
        if (rhs.constraintValue > 0) {
            // This solution is feasible, rhs is not
            return true;
        } else {
            // Both solutions are feasible
            return fitness > rhs.fitness;
        }
    }
}

/**
 * Determines if the solutions are equal in fitness
 * @param rhs
 * @return
 */
bool solution_BN::thisSolutionHasEqualFitness(const solution_BN &rhs) {
    return ((this->constraintValue == rhs.constraintValue) && (this->fitness == rhs.fitness));
}


bool solution_BN::operator==(const solution_BN &rhs) const {
    // https://stackoverflow.com/questions/16422486/can-i-use-to-compare-two-vectors-i-tried-it-and-seems-to-be-working-fine
    return variables == rhs.variables && c_variables == rhs.c_variables;
}

bool solution_BN::operator!=(const solution_BN &rhs) const {
    return !(rhs == *this);
}

///////////////
/// Cloning ///
///////////////

/**
 * Clones this solution
 * @return A clone of this solution
 */
solution_BN *solution_BN::clone() {
    // shared_ptr<DiscretizationPolicy> copyOfDiscretizationPolicy = this->discretizationPolicy->clone();
    /*solution_BN *result = new solution_BN(this->variables, 
                                          this->fitness_buffers, 
                                          getObjectiveValues(), 
                                          getConstraintValue(), 
                                          getAlphabetSize(), 
                                          this->c_variables, 
                                          this->problemInstance,
                                          this->node_data_types,
                                          this->number_of_nodes,
                                          this->number_of_links,
                                          this->number_of_nodes_to_discretize,
                                          this->maximum_number_of_parents,
                                          this->fitness,
                                          this->constraintValue,
                                          this->timeStamp,
                                          this->numberOfFullEvaluations,
                                          this->numberOfEvaluations,
                                          this->number_of_parents,
                                          this->number_of_children,
                                          this->child_matrix,
                                          this->parent_matrix,
                                          this->adjacency_matrix,
                                          this->spouse_matrix,
                                          this->numberOfDiscretizationsperNode,
                                          this->discretizationPolicy, //copyOfDiscretizationPolicy,
                                          this->boundaries,
                                          this->maxValuesData,
                                          this->minValuesData);
    return result;*/
    return new solution_BN(*this);
}


/**
 * Constructor for cloning things
 */
/*solution_BN::solution_BN(vec_t<int> &variables, vec_t<double> fitness_buffers, vec_t<double> objective_values, double constraint_value, 
    size_t alphabetSize, vec_t<double> &c_variables, fitness_t<int> *problemInstance,
    const vec_t<ColumnDataType> &nodeDataTypes, int numberOfNodes,
    size_t numberOfLinks, size_t numberOfNodesToDiscretize, int maximumNumberOfParents,
    double fitness, double constraintValue,
    clock_t timeStamp, size_t numberOfFullEvaluations, double numberOfEvaluations,
    const vec_t<size_t> &numberOfParents, const vec_t<size_t> &numberOfChildren,
    const vec_t<vec_t<int>> &childMatrix, const vec_t<vec_t<int>> &parentMatrix,
    const vec_t<vec_t<int>> &adjacencyMatrix, const vec_t<vec_t<vec_t<size_t>>> &spouse_matrix,
    const vec_t<size_t> &numberOfDiscretizationsperNode,
    shared_ptr<DiscretizationPolicy> &discretizationPolicy,
    vec_t<vec_t<double>> boundaries, 
    vec_t<double> maxValuesData, vec_t<double> minValuesData) : solution_mixed(variables, fitness_buffers, objective_values, constraint_value, alphabetSize, c_variables, problemInstance),
    node_data_types(nodeDataTypes), number_of_nodes(numberOfNodes), number_of_links(numberOfLinks),
    number_of_nodes_to_discretize(numberOfNodesToDiscretize), maximum_number_of_parents(maximumNumberOfParents),
    fitness(fitness), constraintValue(constraintValue),
    timeStamp(timeStamp), numberOfFullEvaluations(numberOfFullEvaluations), numberOfEvaluations(numberOfEvaluations),
    number_of_parents(numberOfParents), number_of_children(numberOfChildren),
    child_matrix(childMatrix), parent_matrix(parentMatrix),
    adjacency_matrix(adjacencyMatrix), spouse_matrix(spouse_matrix),
    numberOfDiscretizationsperNode(numberOfDiscretizationsperNode), discretizationPolicy(discretizationPolicy),
    boundaries(boundaries), maxValuesData(maxValuesData), minValuesData(minValuesData) {}*/

solution_BN::solution_BN( const solution_BN &other ) : solution_mixed(other),
    node_data_types(other.node_data_types), number_of_nodes(other.number_of_nodes), number_of_links(other.number_of_links),
    number_of_nodes_to_discretize(other.number_of_nodes_to_discretize), maximum_number_of_parents(other.maximum_number_of_parents),
    fitness(other.fitness), constraintValue(other.constraintValue),
    timeStamp(other.timeStamp), numberOfFullEvaluations(other.numberOfFullEvaluations), numberOfEvaluations(other.numberOfEvaluations),
    number_of_parents(other.number_of_parents), number_of_children(other.number_of_children),
    child_matrix(other.child_matrix), parent_matrix(other.parent_matrix),
    adjacency_matrix(other.adjacency_matrix), spouse_matrix(other.spouse_matrix),
    numberOfDiscretizationsperNode(other.numberOfDiscretizationsperNode), discretizationPolicy(other.discretizationPolicy),
    boundaries(other.boundaries), maxValuesData(other.maxValuesData), minValuesData(other.minValuesData), 
    lower_user_range(other.lower_user_range), upper_user_range(other.upper_user_range), data(other.data),
    useNormalizedCVars(other.useNormalizedCVars), transformCVariables(other.transformCVariables), 
    useOptimalSolution(other.useOptimalSolution), guaranteedInitSpread(other.guaranteedInitSpread), extraCVarForNumberOfBins(other.extraCVarForNumberOfBins),
    problemInstancePath(other.problemInstancePath), runIndex(other.runIndex)  {}


/*solution_BN::solution_BN( const solution_BN &other ) : solution_mixed(other.variables, other.fitness_buffers, objective_values, constraint_value, alphabetSize, c_variables, problemInstance)
{
    this->node_data_types = other.nodeDataTypes;
    this->number_of_nodes = other.numberOfNodes;
    this->number_of_links = other.numberOfLinks;
    this->number_of_nodes_to_discretize = other.numberOfNodesToDiscretize;
    this->maximum_number_of_parents = other.maximumNumberOfParents;
    this->fitness = other.fitness;
    this->constraintValue = other.constraintValue;
    this->timeStamp = other.timeStamp;
    this->numberOfFullEvaluations = other.numberOfFullEvaluations; 
    this->numberOfEvaluations = other.numberOfEvaluations;
    this->number_of_parents = other.numberOfParents;
    this->number_of_children other.numberOfChildren;
    this->child_matrix = other.childMatrix;
    this->parent_matrix = other.parentMatrix;
    this->adjacency_matrix = other.adjacencyMatrix;
    this->spouse_matrix = other.spouse_matrix;
    this->numberOfDiscretizationsperNode = other.numberOfDiscretizationsperNode;
    this->discretizationPolicy = other.discretizationPolicy;
    this->boundaries = other.boundaries;
    this->maxValuesData = other.maxValuesData;
    this->minValuesData = other.minValuesData;
}
*/
}


