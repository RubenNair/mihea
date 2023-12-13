#pragma once

#include "gomea/src/common/gomea_defs.hpp"
#include "gomea/src/common/partial_solution.hpp"
// #include "gomea/src/mixed_integer/Config.hpp"
#include "gomea/src/common/solution_mixed.hpp"
#include "gomea/src/utils/data_structure.h"
#include "gomea/src/common/gomea_defs.hpp"
#include "gomea/src/fitness/fitness.hpp"
#include "gomea/src/fitness/discretization.h"

#include <algorithm>

namespace gomea {


using gomea::fitness::fitness_t;

/**
 * Used to store information about the nodes
 */
struct NodeInformation {
    vec_t<size_t> output_number_of_parents;        // Number of parents per node i
    vec_t<size_t> output_number_of_children;       // Number of children per node i
    vec_t<vec_t<int>> output_parent_matrix;       // Per row i,  the (node indices) of the parents of node i
    vec_t<vec_t<int>> output_child_matrix;        // Per row i, the (node) indices of the children of node i
};


/**
 * Used to compactly return information from a processing function
 */
struct NetworkStructure {
    vec_t<int> processedParameters;            // The parameters after removing excess parents and cycles
    NodeInformation nodeInfo;                   // The node information after removing excess parents and cycles
    vec_t<vec_t<int>> adjacencyMatrix;        // The Adjacency matrix after removing excess parents and cycles
    vec_t<size_t> discretizationPerNode;       // The number of discretizations per node
    shared_ptr<DiscretizationPolicy> discretizationPolicy;      // The discretization policy of this solution
};



class solution_BN : public solution_mixed
{
	public:
	// `variables` contains the discrete part of the solution: e.g.: [A_12, A13, A14, ..., A23, A24, ..., A34, ..., | D0, D1, ...]
										// RUBEN -> instead of D0, D1, ... we store the bin width (normalized) for each node (in other variable, `c_variables`). Requires knowing max number of discretizations per node.
    // solution_BN( int number_of_variables_, size_t alphabetSize_, int number_of_c_variables_ );
    // solution_BN (vec_t<int> &variables, vec_t<double> &c_variables);
    solution_BN( size_t numberOfVariables_, size_t alphabetSize_, size_t numberOfCVariables_, fitness_t<int> *problemInstance_ );
    solution_BN( size_t numberOfVariables_, 
                 size_t alphabetSize_, 
                 size_t numberOfCVariables_, 
                //  vec_t<int> &unprocessedParameters,
                 vec_t<ColumnDataType> node_data_types,
                 const vec_t<size_t> &initialNumberOfInstantiations,
                 int discretizationPolicyIndex,
                 size_t maximum_number_of_parents,
                 size_t maximum_number_of_Instantiations,
                 fitness_t<int> *problemInstance_,
                 vec_t<double> maxValuesData, vec_t<double> minValuesData,
                 double lower_user_range, double upper_user_range,
                 shared_ptr<DataStructure<double>> data = NULL,
                 double populationIndexRatio = -1.0,
                 bool useNormalizedCVars = false,
                 bool transformCVariables = false,
                 bool useOptimalSolution = false,
                 bool guaranteedInitSpread = false,
                 bool extraCVarForNumberOfBins = false,
                 bool forceNBoundariesUsed = false,
                 string problemInstancePath = "",
                 int runIndex = 0);
    virtual ~solution_BN() = default;

    void reProcessParametersSolution(vector<int> newParameters);
    
    bool operator==(const solution_BN &solutionB)
    {
        return (this->variables == solutionB.variables && this->c_variables == solutionB.c_variables);
    }

	bool operator!=(const solution_BN &solutionB)
	{
		return !(*this == solutionB);
	}

	friend std::ostream &operator<<(std::ostream &out, const solution_BN &solution)
	{
		for (int i = 0; i < solution.getNumberOfVariables(); ++i)
			out << +solution.variables[i];

		for(int i = 0; i < solution.getNumberOfCVariables(); ++i)
			out << +solution.c_variables[i];
		out << " | " << solution.getObjectiveValue();
		return out;
	}

	int getNumberOfCVariables() const;
	void randomInit(std::mt19937 *rng);
    void randomInit(std::mt19937 *rng, double populationIndexRatio);
    void optimalInit();
    void normalize(int numberOfBins = -1);
    void updateBoundaries();
    tuple<vec_t<double>, vec_t<double>> findMaxAndMinValuesInData();
	void insertCVariables( const vec_t<double> &vars_to_insert );
	void insertCVariables( vec_t<double> vars_to_insert, vec_t<int> indices_to_insert );
	void insertSolution( solution_mixed *solution );
	void print();

	fitness_t<int> *problemInstance;

    // Copy
    solution_BN *clone();

    // Initialize and clear up memory: For large problems, memory can become an issue
    void reinitializeNetworkVariables();
    // Clears the network variables
    void clearNetworkVariables();

    /// Network derivations
    // Determines the spouse nodes of the network
    void determineSpouseNodes(vec_t<size_t> numberOfParents, vec_t<size_t> numberOfChildren, vec_t<vec_t<int>> parentMatrix, vec_t<vec_t<int>> childMatrix);
    // Determines the max number of parent combinations per node
    vec_t<size_t> determineMaxNumberOfParentCombinations(vec_t<size_t> numberOfInstantiations);
    // Determines the number of spouse combinations the children can make
    vec_t<vec_t<size_t>> determineMaxNumberOfSpouseCombinations(vec_t<size_t> numberOfInstantiations);


    // Comparison
    bool thisSolutionIsBetterThan(const solution_BN &rhs);      // Compares the fitness and constraintValue with other solutions
    bool thisSolutionHasEqualFitness(const solution_BN &rhs);   // Compares for equality between fitness and constraint value
    bool operator==(const solution_BN &rhs) const;              // The equality operator
    bool operator!=(const solution_BN &rhs) const;              // The difference operator

    // Getters
    vec_t<int> getParameters() const;
    vec_t<int> getNetworkParameters() const;
    int getNumberOfNodes() const;
    size_t getNumberOfLinks() const;
    size_t getNumberOfNodesToDiscretize() const;
    double getFitness() const;
    double getConstraintValue() const;
    const vec_t<double> &getFitnessPerNode() const;
    int getMaximumNumberOfParents() const;
    const vec_t<size_t> &getNumberOfParents() const;
    const vec_t<size_t> &getNumberOfChildren() const;
    const vec_t<vec_t<int>> &getParentMatrix() const;
    const vec_t<vec_t<int>> &getChildMatrix() const;
    const vec_t<vec_t<vec_t<size_t>>> &getSpouseMatrix() const;
    const vec_t<vec_t<int>> &getAdjacencyMatrix() const;
    const vec_t<ColumnDataType> &getNodeDataTypes() const;
    const vec_t<size_t> &getNumberOfDiscretizationsperNode() const;
    const shared_ptr<DiscretizationPolicy> &getDiscretizationPolicy() const;
    clock_t getTimeStamp() const;
    size_t getNumberOfFullEvaluations() const;
    double getNumberOfEvaluations() const;
    const vec_t<vec_t<double>> &getBoundaries() const;
    const vec_t<double> &getMaxValuesData() const;
    const vec_t<double> &getMinValuesData() const;
    const shared_ptr<DataStructure<double>> &getDiscretizedData() const;

    void setDiscretizedData(const shared_ptr<DataStructure<double>> &data);

    // Setters
    void setFitness(double newFitnessValue);
    void setFitnessPerNode(const vec_t<double> &newFitnessPerNode);
    void setConstraintValue(double newConstraintValue);
    void setTimeStamp(clock_t timeStamp);
    void setNumberOfFullEvaluations(size_t numberOfFullEvaluations);
    void setNumberOfEvaluations(double numberOfEvaluations);
    void setDiscretizationPolicy(const shared_ptr<DiscretizationPolicy> &discretizationPolicy);

    // Constructor for cloning (If you add something here make sure you're not just copying pointers).
    /*solution_BN(vec_t<int> &variables, vec_t<double> fitness_buffers, vec_t<double> objective_values, double constraint_value, 
                size_t alphabetSize, vec_t<double> &c_variables, fitness_t<int> *problemInstance,
                const vec_t<ColumnDataType> &nodeDataTypes, int numberOfNodes,
                size_t numberOfLinks, size_t numberOfNodesToDiscretize, int maximumNumberOfParents, double fitness,
                double constraintValue,
                clock_t timeStamp, size_t numberOfFullEvaluations, double numberOfEvaluations,
                const vec_t<size_t> &numberOfParents, const vec_t<size_t> &numberOfChildren,
                const vec_t<vec_t<int>> &childMatrix, const vec_t<vec_t<int>> &parentMatrix,
                const vec_t<vec_t<int>> &adjacencyMatrix, const vec_t<vec_t<vec_t<size_t>>> &spouse_matrix,
                const vec_t<size_t> &numberOfDiscretizationsperNode,
                shared_ptr<DiscretizationPolicy> &discretizationPolicy,
                vec_t<vec_t<double>> boundaries,
                vec_t<double> maxValuesData, vec_t<double> minValuesData);*/
    solution_BN( const solution_BN &other );

protected:
    shared_ptr<DataStructure<double>> data = NULL; // Store the data, to be able to calculate the boundaries based on c_vars and datapoints.
    shared_ptr<DataStructure<double>> discretizedData = NULL; // Store the discretized data, used for fitness calculation (previously stored in fitness, but gives issues with parallel evaluation)
    vec_t<ColumnDataType> node_data_types; // The type of the node [Discrete/Continuous]
    int number_of_nodes;                    // The number of random variables (nodes)
    size_t number_of_links;                 // The number of links in the bayesian network
    size_t number_of_nodes_to_discretize;   // The number of nodes to discretize 
    int maximum_number_of_parents;          // Maximum number of children a node can have

    // Fitness
    double fitness;                         // The fitness of the solution
    double constraintValue;                 // The constraint value of the solution
    clock_t timeStamp;                      // The time stamp when a solution was evaluated
    size_t numberOfFullEvaluations;         // The number of full evaluations executed to get the solution
	double numberOfEvaluations;             // The number of (partial) evaluations executed to get this solution

    // Init lower and upper bound
    double lower_user_range;
    double upper_user_range;

	// Bayesian network structure
    vec_t<size_t> number_of_parents;       // The number of parents of each random variable
    vec_t<size_t> number_of_children;      // The number of children of each random variable
    vec_t<vec_t<int>> child_matrix;       // Matrix indicating outgoing links (per row). See the impl. 'findParents()'.
    vec_t<vec_t<int>> parent_matrix;      // Matrix indicating ingoing links (per row). See the impl. 'findParents()'.
    vec_t<vec_t<int>> adjacency_matrix;   // The adjacency matrix. See the impl. 'findAdjacencyMatrix()'.
    vec_t<vec_t<vec_t<size_t>>> spouse_matrix;   // Matrix indicating the spouses of each child node of node i

    void updateBoundariesBasedOnBinWidths();
    void updateBoundariesBasedOnNumberOfDataSamples();
    void execTransformationCVariables();

    // Solution processing methods
    NetworkStructure processParametersSolution(const vector<int> &parametersToProcess,
                                               const vector<size_t> &initialNumberOfInstantiations,
                                               int discretizationPolicyIndex,
                                               size_t maximum_number_of_instantiations);                        // Processes the parameters of the solution
    NetworkStructure processNetworkStructure(vector<int> parametersToProcess);                                  // Process the parameters related to the network structure
    NetworkStructure processDiscretization(NetworkStructure result,
                                           vector<size_t> initialNumberOfInstantiations,
                                           int discretizationPolicyIndex,
                                           size_t maximum_number_of_instantiations);                                    // Process the parameters related to the discretization of the continuous data
    vector<size_t> determineNumberOfDiscretizationsFromSolution(const NetworkStructure& networkStructure,
                                                                size_t &rawSolutionPointerIndex,
                                                                const vector<size_t> &initialNumberOfInstantiations);   // Process the number of discretizations of the continuous nodes

    NodeInformation findParents(vector<int> &parametersToProcess);                  // Determines the number of children 'number_of_children' and the indices of the children 'child_matrix'
    vector<vector<int>> findAdjacencyMatrix(vector<size_t> numberOfChildren,
                                            vector<vector<int>> parentMatrix);      // Determines the adjacency matrix 'adjacency_matrix'
    vector<size_t> determineSpousesOfChildNode(vector<int> parentIndicesOfChild,
                                               size_t numberOfParents,
                                               size_t nodeIndexToExclude);          // Determines the spouse node given a child node
    void topologicalSort(vector<int> &parametersToProcess,
                         vector<size_t> numberOfParents,
                         vector<vector<int>> adjacencyMatrix);                      // Does a topological sort that removes links that make the graph cyclic

    // Detecting cyclic graphs
    void dfs(size_t currentNodeIndex,
             vector<vector<int>> &tempAdjacentMatrix,
             vector<int> &visitedNodes,
             vector<int> &nodesInStack,
             vector<int> &preprocessedParameters);                                      // Does a depth first search to detect cycles
    void removeEdge(int currentNodeIndex, int u, vector<int> &preprocessedParameters);  // Removes edges that cause a cycle in the graph

	// Discretization // RUBEN -> TODO: figure out how to make this work with my discretization policy / representation
    vec_t<size_t> numberOfDiscretizationsperNode;                       // The number of discretizations per node
    shared_ptr<DiscretizationPolicy> discretizationPolicy;              // The discretization policy
    vec_t<vec_t<double>> boundaries;                                    // Stores the boundaries of the discretization (per node) TODO RUBEN: Make sure this is updated if continuous variables are updated.
    vec_t<double> maxValuesData, minValuesData;                         // The maximum and minimum values in the data for each node
    bool useNormalizedCVars;                                            // Indicates whether the c_vars are normalized or not (and also if boundaries are based on number of samples or on data ranges (bin widths), respectively)
    bool transformCVariables;                                           // Indicates whether c_vars are transformed for iamalgam (initialized in [0,1], also calculate boundaries in this range, but otherwise map using f(x) = 1/x)
    bool useOptimalSolution;
    bool guaranteedInitSpread;                                          // If true, the c_vars will be initialized such that there is an equal number of solutions in the population for each amount of bins (between 2 and max)
    bool extraCVarForNumberOfBins;                                      // If true, an extra c_var is added to the solution per continuous node to indicate the number of bins (between 2 and max)
    bool forceNBoundariesUsed;                                          // If true, the amount of bins used will be equal to the value of the extra c-variable(s).
    string problemInstancePath = "";                                    // The path to the problem instance
    int runIndex;                                                       
    


// {
// 	for (int i = 0; i < getNumberOfCVariables(); ++i)
// 	{
// 		cVariables[i] = problemInstance->getLowerRangeBound(i) + ((*rng)() / (double)(*rng).max()) * (problemInstance->getUpperRangeBound(i) - problemInstance->getLowerRangeBound(i));
// 	}
// }

	
	// 	solution_t( int number_of_variables );
	// 	solution_t( vec_t<char> &variables );
	// 	solution_t( size_t numberOfVariables_, size_t alphabetSize_ );
	// 	solution_t(size_t numberOfDVariables_, size_t numberOfCVariables_, size_t alphabetSize_, fitness_t<xhar> *problemInstance_);
	// 	// TODO RUBEN: need a solution_t initializer for continuous; will need to pass 

	// 	bool operator==(const solution_t<T> &solutionB)
	// 	{
	// 		for (int i = 0; i < getNumberOfVariables(); ++i)
	// 		{
	// 			if (this->variables[i] != solutionB.variables[i])
	// 				return false;
	// 		}
	// 		return true;
	// 	}

	// 	friend std::ostream &operator<<(std::ostream &out, const solution_t<T> &solution)
	// 	{
	// 		for (int i = 0; i < solution.getNumberOfVariables(); ++i)
	// 			out << +solution.variables[i];
	// 		out << " | " << solution.getObjectiveValue();
	// 		return out;
	// 	}
		
	// 	void randomInit(std::mt19937 *rng);

	// 	void initMemory(int number_of_objectives, int number_of_fitness_buffers);
	// 	void initObjectiveValues(int number_of_objectives);
	// 	void initFitnessBuffers(int number_of_fitness_buffers);

	// 	int getNumberOfVariables() const;
	// 	int getNumberOfDVariables() const;
	// 	int getNumberOfCVariables() const;
	// 	int getNumberOfObjectives() const;
	// 	double getObjectiveValue( int objective_value_index = 0 ) const;
	// 	const vec_t<double> getObjectiveValues() const;
	// 	double getPartialObjectiveValue( int subfunction_index ) const;
	// 	double getConstraintValue() const;
	// 	double getPartialConstraintValue( int subfunction_index ) const;

	// 	void setObjectiveValue( double v );
	// 	void setObjectiveValue( int objective_value_index, double v );
	// 	void setObjectiveValues( const vec_t<double> &v );
	// 	void setConstraintValue( double v );
	// 	void setPartialObjectiveValue( int subfunction_index, double v );
	// 	void setPartialConstraintValue( int subfunction_index, double v );

	// 	double getFitnessBuffer( int buffer_index ) const;
	// 	const vec_t<double> getFitnessBuffers() const;
	// 	void addToFitnessBuffer( int buffer_index, double partial_fitness );
	// 	void subtractFromFitnessBuffer( int buffer_index, double partial_fitness );
	// 	void setFitnessBuffers( const vec_t<double> &buffers );
	// 	void clearFitnessBuffers();

	// 	partial_solution_t<T> getPartialCopy( const vec_t<int> &variable_indices ) const;
	// 	const vec_t<T> getCopyOfVariables( const vec_t<int> &variable_indices = vec_t<int>()) const;
	// 	void insertVariables( const vec_t<T> &vars_to_insert );
	// 	void insertVariables(vec_t<T> vars_to_insert, vec_t<int> indices_to_insert);
	// 	void insertSolution( solution_t<T> *solution );
	// 	void insertPartialSolution( partial_solution_t<T> *solution );

	// 	void print();

	// 	size_t getAlphabetSize();
		
	// 	vec_t<T> variables;
	// 	vec_t<double> fitness_buffers;
	// 	vec_t<T> dVariables;
	// 	bool is_MI_solution = false;
	// 	fitness_t *problemInstance;

	// private:
	// 	vec_t<double> objective_values;
	// 	double constraint_value;
		
	// 	vec_t<double> partial_objective_values;
	// 	vec_t<double> partial_constraint_values;

	// 	size_t alphabetSize = 0;
};

}
