#include "gomea/src/mixed_integer/Population.hpp"

namespace gomea{namespace mixedinteger{

// Saving old way of initializing Population for now, without cPopulation
// Population::Population(Config *config_, fitness_t *problemInstance_, sharedInformation *sharedInformationPointer_, size_t GOMEAIndex_, size_t populationSize_, linkage_model_pt FOSInstance_ ): 
//         config(config_), 
//         problemInstance(problemInstance_),
//         sharedInformationPointer(sharedInformationPointer_),
//         GOMEAIndex(GOMEAIndex_), 
//         populationSize(populationSize_)
// {
//         terminated = false;
//         numberOfGenerations = 0;
//         averageFitness = 0.0;
        
//         dPopulation.resize(populationSize);
//         offspringdPopulation.resize(populationSize);
//         noImprovementStretches.resize(populationSize);
        
//         vec_t<int> allGenes(problemInstance->number_of_variables);
//         iota(allGenes.begin(), allGenes.end(), 0);

//         for (size_t i = 0; i < populationSize; ++i)
//         {
//             noImprovementStretches[i] = 0;

//             dPopulation[i] = new solution_t<char>(problemInstance->number_of_variables, config->alphabetSize);
//             dPopulation[i]->randomInit(&gomea::utils::rng);
//             problemInstance->evaluate(dPopulation[i]);
            
//             offspringdPopulation[i] = new solution_t<char>(problemInstance->number_of_variables, config->alphabetSize);
//             *offspringdPopulation[i] = *dPopulation[i];
//         }
			
// 		if( config->linkage_config != NULL )
// 		{
// 			FOSInstance = linkage_model_t::createFOSInstance( *config->linkage_config, problemInstance->number_of_variables );
//             FOSInstance->initializeDependentSubfunctions( problemInstance->subfunction_dependency_map );
// 		}
// 		else if( FOSInstance_ == NULL )
// 		{
// 			FOSInstance = linkage_model_t::createLinkageTreeFOSInstance(config->FOSIndex, problemInstance->number_of_variables, config->linkage_config->lt_similarity_measure, config->linkage_config->lt_maximum_set_size);
// 		}
// 		else FOSInstance = FOSInstance_;
        
//         #ifdef DEBUG
//             cout << "New Population created! Population #" << GOMEAIndex << " populationSize:" << populationSize << endl;
//             cout << this;
//         #endif
// }

Population::Population(Config *config_, fitness_t *problemInstance_, sharedInformation *sharedInformationPointer_, size_t GOMEAIndex_, size_t populationSize_, linkage_model_pt FOSInstance_ ): 
        config(config_), 
        problemInstance(problemInstance_),
        sharedInformationPointer(sharedInformationPointer_),
        GOMEAIndex(GOMEAIndex_), 
        populationSize(populationSize_)
{
        terminated = false;
        numberOfGenerations = 0;
        averageFitness = 0.0;
        
        population.resize(populationSize);
        offspringPopulation.resize(populationSize);
        noImprovementStretches.resize(populationSize);
        
        vec_t<int> allGenes(problemInstance->number_of_variables); // RUBEN TODO might need to change this property of problemInstance, but depends on how fitness_t is handled
        iota(allGenes.begin(), allGenes.end(), 0);

        for (size_t i = 0; i < populationSize; ++i)
        {
            noImprovementStretches[i] = 0;

            population[i] = new solution_mixed(config->numberOfVariables, config->alphabetSize, config->numberOfcVariables, problemInstance); // TODO RUBEN: now just assumes amount of d and c variables is equal, getting real number of c variables depends on how fitness_t is handled
            population[i]->randomInit(&gomea::utils::rng);
            problemInstance->evaluate(population[i]);
            // updateElitistAndCheckVTR(population[i]); // TODO RUBEN I think this can be removed, since this call is made at the start of learnDiscreteModel (after calculating average fitness)
            
            // problemInstance->evaluate(dPopulation[i]); // RUBEN TODO maybe doesn't work like this anymore, depends on how fitness_t is handled
            
            offspringPopulation[i] = new solution_mixed(config->numberOfVariables, config->alphabetSize, config->numberOfcVariables, problemInstance);
            *offspringPopulation[i] = *population[i];
        }
			
		if( config->linkage_config != NULL && config->numberOfVariables > 0 )
		{
			FOSInstance = linkage_model_t::createFOSInstance( *config->linkage_config, problemInstance->number_of_variables );
            FOSInstance->initializeDependentSubfunctions( problemInstance->subfunction_dependency_map );
		}
		else if( FOSInstance_ == NULL )
		{
			FOSInstance = linkage_model_t::createLinkageTreeFOSInstance(config->FOSIndex, problemInstance->number_of_variables, config->linkage_config->lt_similarity_measure, config->linkage_config->lt_maximum_set_size);
		}
		else FOSInstance = FOSInstance_;
        
        if(config->numberOfcVariables > 0)
        {
            iamalgamInstance = new iamalgam(config_, offspringPopulation);
            iamalgamInstance->initialize();
        }



        #ifdef DEBUG
            cout << "New Population created! Population #" << GOMEAIndex << " populationSize:" << populationSize << endl;
            cout << this;
        #endif
}

Population::~Population()
{
    for (size_t i = 0; i < populationSize; ++i)
    {
        delete population[i];
        // delete offspringPopulation[i];
    }
}

ostream & operator << (ostream &out, const Population &populationInstance)
{
    out << "Generation " << populationInstance.numberOfGenerations << ":" << endl;
    for (size_t i = 0; i < populationInstance.populationSize; ++i)
        out << *populationInstance.population[i] << endl;
    out << endl;
    return out;
}

bool Population::allSolutionsAreEqual()
{
    for (size_t i = 1; i < populationSize; i++)
    {
        if(*population[i] != *population[0])
            return false;
    }
    cout << "[DEBUGGING] All solutions are equal! Population #" << GOMEAIndex << endl;
    return true;
}

void Population::calculateAverageFitness()
{
    averageFitness = 0.0;
    for (size_t i = 0; i < populationSize; ++i)
        averageFitness += population[i]->getObjectiveValue();
    averageFitness /= populationSize;
}

double Population::getFitnessMean()
{
	double objective_avg = 0.0;
	for(size_t i = 0; i < populationSize; i++ )
		objective_avg  += population[i]->getObjectiveValue();
	objective_avg = objective_avg / ((double) populationSize);
	return( objective_avg );
}

double Population::getFitnessVariance()
{
	double objective_avg = getFitnessMean();
	double objective_var = 0.0;
	for(size_t i = 0; i < populationSize; i++ )
		objective_var  += (population[i]->getObjectiveValue()-objective_avg)*(population[i]->getObjectiveValue()-objective_avg);
	objective_var = objective_var / ((double) populationSize);

	if( objective_var <= 0.0 )
		objective_var = 0.0;
	return( objective_var );
}

double Population::getConstraintValueMean()
{
	double constraint_avg = 0.0;
	for(size_t i = 0; i < populationSize; i++ )
		constraint_avg  += population[i]->getConstraintValue();
	constraint_avg = constraint_avg / ((double) populationSize);

	return( constraint_avg );
}

double Population::getConstraintValueVariance()
{
	double constraint_avg = getConstraintValueMean();

	double constraint_var = 0.0;
	for(size_t i = 0; i < populationSize; i++ )
		constraint_var  += (population[i]->getConstraintValue()-constraint_avg)*(population[i]->getConstraintValue()-constraint_avg);
	constraint_var = constraint_var / ((double) populationSize);

	if( constraint_var <= 0.0 )
		constraint_var = 0.0;
	return( constraint_var );
}

solution_mixed *Population::getBestSolution()
{
	int index_best = 0;
	for(size_t j = 1; j < populationSize; j++ )
    {
        if( problemInstance->betterFitness( population[j]->getObjectiveValue(), population[j]->getConstraintValue(), population[index_best]->getObjectiveValue(), population[index_best]->getConstraintValue()) )
		{
			index_best = j;
        }
    }
	return( population[index_best] );
}

solution_mixed *Population::getWorstSolution()
{
	int index_worst = 0;
	for(size_t j = 1; j < populationSize; j++ )
    {
        if( problemInstance->betterFitness( population[index_worst]->getObjectiveValue(), population[index_worst]->getConstraintValue(), population[j]->getObjectiveValue(), population[j]->getConstraintValue()) )
		{
			index_worst = j;
        }
    }
	return( population[index_worst] );
}


void Population::copyOffspringToPopulation()
{
    for(size_t i = 0; i < populationSize; i++)
    {
        // *population[i] = *offspringPopulation[i];
        population[i]->insertSolution(offspringPopulation[i]);
    }
    // *population[0] = sharedInformationPointer->elitist;
}

void Population::copyPopulationToOffspring()
{
    for(size_t i = 0; i < populationSize; i++)
    {
        // *population[i] = *offspringPopulation[i];
        offspringPopulation[i]->insertSolution(population[i]);
    }
    // *population[0] = sharedInformationPointer->elitist;
}

void Population::makeOffspring()
{
    if( numberOfGenerations == 0 )
    {
        for (size_t i = 0; i < populationSize; ++i)
            updateElitistAndCheckVTR(population[i]);
    }

    if( FOSInstance->type == linkage::LINKAGE_TREE )
    {
        if (FOSInstance->is_static)
        {
            if (FOSInstance->size() == 0)
            {
                FOSInstance->learnLinkageTreeFOS(problemInstance->getSimilarityMatrix(FOSInstance->getSimilarityMeasure()), false);
                FOSInstance->initializeDependentSubfunctions( problemInstance->subfunction_dependency_map );
            }
        }
        else
        {
            // TODO RUBEN figure out if this is a valid solution to the problem of matching the function argument of learnLinkageTreeFOS
            vec_t<solution_t<char>*> casted_population;
            for(solution_mixed* sol : population) 
            {
                casted_population.push_back(static_cast<solution_t<char>*>(sol));
            }

            FOSInstance->learnLinkageTreeFOS(casted_population, config->alphabetSize );
            FOSInstance->initializeDependentSubfunctions(problemInstance->subfunction_dependency_map);
        }
    }

    FOSInstance->setCountersToZero();
    if (config->AnalyzeFOS)
    {
        FOSInstance->writeToFileFOS(config->folder, GOMEAIndex, numberOfGenerations);
    }
 
    generateOffspring();

    if (config->AnalyzeFOS)
        FOSInstance->writeFOSStatistics(config->folder, GOMEAIndex, numberOfGenerations);

}

void Population::learnDiscreteModel()
{
    // TODO: maybe update the similarity matrix here first? -> If I understand correctly, similarity matrix is made once and not updated?
            // Seems to be based on variable interaction graph. Just TODO need to double-check that it is initialized correctly (but maybe not in this function)

    calculateAverageFitness();

    if( numberOfGenerations == 0 )
    {
        for (size_t i = 0; i < populationSize; ++i)
            updateElitistAndCheckVTR(population[i]);
    }

    // checkForDuplicate("DISCRETE START LEARN MODEL");
    if( FOSInstance->type == linkage::LINKAGE_TREE )
    {
        if (FOSInstance->is_static)
        {
            if (FOSInstance->size() == 0)
            {
                FOSInstance->learnLinkageTreeFOS(problemInstance->getSimilarityMatrix(FOSInstance->getSimilarityMeasure()), false);
            }
        }
        else
        {
            // TODO RUBEN figure out if this is a valid solution to the problem of matching the function argument of learnLinkageTreeFOS
            vec_t<solution_t<char>*> casted_population;
            for(solution_mixed* sol : population) 
            {
                casted_population.push_back(static_cast<solution_t<char>*>(sol));
            }

            FOSInstance->learnLinkageTreeFOS(casted_population, config->alphabetSize );
        }

        // For mixed-integer problems, the FOS set of all discrete elements is still valid (since it is not the whole solution (continuous part), and thus not a complete copy would be made)
        // So add the FOS set of all discrete elements to the FOSInstance.
        if(config->numberOfcVariables > 0)
        {
            vec_t<int> allDiscreteElements(problemInstance->number_of_variables);
            for (int i = 0; i < problemInstance->number_of_variables; ++i)
            {
                allDiscreteElements[i] = i;
            }
            FOSInstance->FOSStructure.push_back(allDiscreteElements);
        }
        FOSInstance->initializeDependentSubfunctions(problemInstance->subfunction_dependency_map);

    }

    FOSInstance->setCountersToZero();
    // FOSInstance->printFOS();
}

void Population::determineFOSOrder()
{
    assert( !config->useParallelFOSOrder || !config->fixFOSOrderForPopulation );
   	if( config->fixFOSOrderForPopulation )
    	FOSInstance->shuffleFOS(); 
	else if( config->useParallelFOSOrder )
    {
        assert( problemInstance->hasVariableInteractionGraph() );
		FOSInstance->determineParallelFOSOrder(problemInstance->variable_interaction_graph );
    }

    if (!config->useParallelFOSOrder && !config->fixFOSOrderForPopulation)
            FOSInstance->shuffleFOS();
}

void Population::generateOffspring()
{
    vec_t<vec_t<int> > neighbors;
   
	assert( !config->useParallelFOSOrder || !config->fixFOSOrderForPopulation );
   	if( config->fixFOSOrderForPopulation )
    	FOSInstance->shuffleFOS(); 
	else if( config->useParallelFOSOrder )
    {
        assert( problemInstance->hasVariableInteractionGraph() );
		FOSInstance->determineParallelFOSOrder(problemInstance->variable_interaction_graph );
    }
    int FIcount = 0;
    for (size_t i = 0; i < populationSize; i++)
    {
        if (!config->useParallelFOSOrder && !config->fixFOSOrderForPopulation)
            FOSInstance->shuffleFOS();

        solution_mixed backup = *population[i];

        bool solutionHasChanged;
        solutionHasChanged = GOM(i);

        /* Phase 2 (Forced Improvement): optimal mixing with elitist solution */
        if (config->useForcedImprovements)
        {
            if ((!solutionHasChanged) || (noImprovementStretches[i] > (1 + (log(populationSize) / log(10)))))
            {
                FI(i);
                FIcount++;
            }
        }
        
        if (!(offspringPopulation[i]->getObjectiveValue() > population[i]->getObjectiveValue())) // RUBEN doesn't this assume minimization?
            noImprovementStretches[i]++;
        else
            noImprovementStretches[i] = 0;
        
        // RUBEN maybe this is a nicer, more general approach for the above if statement
        // if((problemInstance->optimization_mode == fitness::opt_mode::MIN && !(offspringPopulation[i]->getObjectiveValue() > population[i]->getObjectiveValue()))
        // || (problemInstance->optimization_mode == fitness::opt_mode::MAX && !(offspringPopulation[i]->getObjectiveValue() < population[i]->getObjectiveValue())))
        //     noImprovementStretches[i]++;
        // else
        //     noImprovementStretches[i] = 0;
    }
    cout << "[DEBUGGING] FIcount: " << FIcount << endl;
}

// Generate offspring for a single FOS element
void Population::generateSingleOffspring(int FOS_index)
{
    int FIcount = 0;
    for (size_t i = 0; i < populationSize; i++)
    {
        
        solution_mixed backup = *offspringPopulation[i];

        bool solutionHasChanged;
        // checkForDuplicate("DISCRETE BEFORE GOM");
        solutionHasChanged = GOMSingleFOS(i, FOS_index);
        // checkForDuplicate("DISCRETE AFTER GOM");
        /* Phase 2 (Forced Improvement): optimal mixing with elitist solution */
        if (config->useForcedImprovements)
        {
            if ((!solutionHasChanged) || (noImprovementStretches[i] > (1 + (log(populationSize) / log(10)))))
            {
                FISingleFOS(i, FOS_index);
                FIcount++;
            }
        }
        
        if (!(offspringPopulation[i]->getObjectiveValue() < population[i]->getObjectiveValue())) // RUBEN doesn't this assume minimization?
            noImprovementStretches[i]++;
        else
            noImprovementStretches[i] = 0;
        
        // RUBEN maybe this is a nicer, more general approach for the above if statement (TODO: check if < and > symbols are correct here)
        // if((problemInstance->optimization_mode == fitness::opt_mode::MIN && !(offspringPopulation[i]->getObjectiveValue() < population[i]->getObjectiveValue()))
        // || (problemInstance->optimization_mode == fitness::opt_mode::MAX && !(offspringPopulation[i]->getObjectiveValue() > population[i]->getObjectiveValue())))
        //     noImprovementStretches[i]++;
        // else
        //     noImprovementStretches[i] = 0;
        // checkForDuplicate("DISCRETE END OF LOOP");
    }

    // if(!allBuildingBlocksStillExist(population, 5))
    // {
    //     cout << "[DEBUGGING] Not all building blocks exist anymore after GOMSingleFOS for FOS_index " << FOS_index << ", gen " << numberOfGenerations << endl;
    // }
    // writeBuildingBlocksToFile(config->folder, population, "BUILDING BLOCKS After GOMSingleFOS index " + to_string(FOS_index) + " ----------------------------------", 5);
}

void Population::evaluateAllSolutionsInPopulation()
{
    for (size_t i = 0; i < populationSize; i++)
    {
        // TODO RUBEN: turn this back on in case I want to use this function as evaluation again as well.
        // problemInstance->evaluate(offspringPopulation[i]);
        updateElitistAndCheckVTR(offspringPopulation[i]);
    }
}

bool Population::GOM(size_t offspringIndex)
{
    size_t donorIndex;
    bool solutionHasChanged = false;
    bool thisIsTheElitistSolution = *offspringPopulation[offspringIndex] == sharedInformationPointer->elitist;//(sharedInformationPointer->elitistSolutionGOMEAIndex == GOMEAIndex) && (sharedInformationPointer->elitistSolutionOffspringIndex == offspringIndex);
    
    *offspringPopulation[offspringIndex] = *population[offspringIndex];
            
    vec_t<int> donorIndices(populationSize);
    iota(donorIndices.begin(), donorIndices.end(), 0);

    for (size_t i = 0; i < FOSInstance->size(); i++)
    {
        int ind = FOSInstance->FOSorder[i];

        if (FOSInstance->elementSize(ind) == 0 || (int) FOSInstance->elementSize(ind) == problemInstance->number_of_variables)
            continue;

        bool donorEqualToOffspring = true;
        size_t indicesTried = 0;

        while (donorEqualToOffspring && indicesTried < donorIndices.size())
        {
            int j = gomea::utils::rng() % (donorIndices.size() - indicesTried);
            swap(donorIndices[indicesTried], donorIndices[indicesTried + j]);
            donorIndex = donorIndices[indicesTried];
            indicesTried++;

            if (donorIndex == offspringIndex)
                continue;

            vec_t<char> donorGenes;
            for(size_t j = 0; j < FOSInstance->elementSize(ind); j++)
            {
                int variableFromFOS = FOSInstance->FOSStructure[ind][j];
                //offspringPopulation[offspringIndex]->variables[variableFromFOS] = population[donorIndex]->variables[variableFromFOS];
                donorGenes.push_back(population[donorIndex]->variables[variableFromFOS]);
                if (donorGenes[variableFromFOS] != offspringPopulation[offspringIndex]->variables[variableFromFOS]) // RUBEN: shouldn't donorGenes[j] be donorGenes[variableFromFOS]? -> Made this change myself
                    donorEqualToOffspring = false;
            }
            partial_solution_t<char> *partial_offspring = new partial_solution_t<char>(donorGenes, FOSInstance->FOSStructure[ind]);

            if (!donorEqualToOffspring)
            {
                //evaluateSolution(offspringPopulation[offspringIndex], backup, touchedGenes, backup->getObjectiveValue());
                //problemInstance->evaluatePartialSolution(offspringPopulation[offspringIndex], partial_offspring, FOSInstance->getDependentSubfunctions(ind) );
                problemInstance->evaluatePartialSolution(offspringPopulation[offspringIndex], partial_offspring );

                // accept the change if this solution is not the elitist and the fitness is at least equally good (allows random walk in neutral fitness landscape)
                // however, if this is the elitist solution, only accept strict improvements, to avoid convergence problems
                if ((!thisIsTheElitistSolution && (partial_offspring->getObjectiveValue() >= offspringPopulation[offspringIndex]->getObjectiveValue())) || 
                        (thisIsTheElitistSolution && (partial_offspring->getObjectiveValue() > offspringPopulation[offspringIndex]->getObjectiveValue())))     
                {
                    offspringPopulation[offspringIndex]->insertPartialSolution(partial_offspring);
                    // offspringPopulation[offspringIndex]->variables[variableFromFOS] = population[donorIndex]->variables[variableFromFOS];
                    //*backup = *offspringPopulation[offspringIndex];
                    
                    solutionHasChanged = true;
                    updateElitistAndCheckVTR(offspringPopulation[offspringIndex]);

                    FOSInstance->improvementCounters[ind]++;
                }

                FOSInstance->usageCounters[ind]++;

            }
            delete partial_offspring;

            break;
        }
    }
    return solutionHasChanged;
}

bool Population::GOMSingleFOS(size_t offspringIndex, size_t FOSIndex)
{
    size_t donorIndex;
    bool solutionHasChanged = false;
    bool thisIsTheElitistSolution = *offspringPopulation[offspringIndex] == sharedInformationPointer->elitist;//(sharedInformationPointer->elitistSolutionGOMEAIndex == GOMEAIndex) && (sharedInformationPointer->elitistSolutionOffspringIndex == offspringIndex);

    // *offspringPopulation[offspringIndex] = *population[offspringIndex];
    // offspringPopulation[offspringIndex]->insertSolution(population[offspringIndex]);
            
    vec_t<int> donorIndices(populationSize);
    iota(donorIndices.begin(), donorIndices.end(), 0);

    int ind = FOSInstance->FOSorder[FOSIndex];

    if (FOSInstance->elementSize(ind) == 0) // RUBEN removed check for (int) FOSInstance->elementSize(ind) == problemInstance->number_of_variables, since for MI case, FOS element with all discrete variables is allowed.
        return solutionHasChanged;

    bool donorEqualToOffspring = true;
    size_t indicesTried = 0;

    while (donorEqualToOffspring && indicesTried < donorIndices.size())
    {
        int j = gomea::utils::rng() % (donorIndices.size() - indicesTried);
        swap(donorIndices[indicesTried], donorIndices[indicesTried + j]);
        donorIndex = donorIndices[indicesTried];
        indicesTried++;

        if (donorIndex == offspringIndex)
            continue;

        vec_t<char> donorGenes = vec_t<char>(offspringPopulation[offspringIndex]->variables.size());
        // First copy all original genes to donorGenes
        for(size_t i = 0; i < offspringPopulation[offspringIndex]->variables.size(); i++)
            donorGenes[i] = offspringPopulation[offspringIndex]->variables[i];

        // Then copy the genes from the donor
        for(size_t j = 0; j < FOSInstance->elementSize(ind); j++)
        {
            int variableFromFOS = FOSInstance->FOSStructure[ind][j];
            //offspringPopulation[offspringIndex]->variables[variableFromFOS] = population[donorIndex]->variables[variableFromFOS];
            donorGenes[variableFromFOS] = population[donorIndex]->variables[variableFromFOS];
            if (donorGenes[variableFromFOS] != offspringPopulation[offspringIndex]->variables[variableFromFOS]) // RUBEN: shouldn't donorGenes[j] be donorGenes[variableFromFOS]? -> Made this change myself
                donorEqualToOffspring = false;
        }
        // partial_solution_t<char> *partial_offspring = new partial_solution_t<char>(donorGenes, FOSInstance->FOSStructure[ind]);
        

        if (!donorEqualToOffspring)
        {   
            solution_mixed *offspring = new solution_mixed(donorGenes, offspringPopulation[offspringIndex]->c_variables);
            problemInstance->evaluate(offspring);
            //evaluateSolution(offspringPopulation[offspringIndex], backup, touchedGenes, backup->getObjectiveValue());
            //problemInstance->evaluatePartialSolution(offspringPopulation[offspringIndex], partial_offspring, FOSInstance->getDependentSubfunctions(ind) );
            // problemInstance->evaluatePartialSolution(offspringPopulation[offspringIndex], partial_offspring );

            // accept the change if this solution is not the elitist and the fitness is at least equally good (allows random walk in neutral fitness landscape)
            // however, if this is the elitist solution, only accept strict improvements, to avoid convergence problems
            // RUBEN changed > to <, since I am assuming minimization. TODO: more elegant solution, maybe using optimization_mode of problemInstance
            if ((!thisIsTheElitistSolution && (offspring->getObjectiveValue() <= offspringPopulation[offspringIndex]->getObjectiveValue())) || 
                    (thisIsTheElitistSolution && (offspring->getObjectiveValue() < offspringPopulation[offspringIndex]->getObjectiveValue())))     
            {
                // offspringPopulation[offspringIndex]->insertPartialSolution(partial_offspring);
                // offspringPopulation[offspringIndex]->variables[variableFromFOS] = population[donorIndex]->variables[variableFromFOS];
                //*backup = *offspringPopulation[offspringIndex];
                offspringPopulation[offspringIndex]->insertSolution(offspring);

                solutionHasChanged = true;
                updateElitistAndCheckVTR(offspringPopulation[offspringIndex]);

                FOSInstance->improvementCounters[ind]++;
            }

            FOSInstance->usageCounters[ind]++;

        }

        break;
    }
    return solutionHasChanged;
}

bool Population::FI(size_t offspringIndex)
{
    if (!config->useParallelFOSOrder && !config->fixFOSOrderForPopulation)
        FOSInstance->shuffleFOS();

    bool solutionHasChanged = 0;

    for (size_t i = 0; i < FOSInstance->size(); i++)
    {
        int ind = FOSInstance->FOSorder[i];
        vec_t<char> touchedGenes = vec_t<char>(FOSInstance->elementSize(ind));
        bool donorEqualToOffspring = true;
        for (size_t j = 0; j < FOSInstance->elementSize(ind); j++)
        {
            int variableFromFOS = FOSInstance->FOSStructure[ind][j];
            touchedGenes[j] = sharedInformationPointer->elitist.variables[variableFromFOS];
            if (population[offspringIndex]->variables[variableFromFOS] != touchedGenes[j])
                donorEqualToOffspring = false;
        }
        gomea::partial_solution_t<char> *partial_offspring = new gomea::partial_solution_t<char>(touchedGenes, FOSInstance->FOSStructure[ind]);

        if (!donorEqualToOffspring)
        {
            problemInstance->evaluatePartialSolution(offspringPopulation[offspringIndex], partial_offspring );

            if (partial_offspring->getObjectiveValue() > offspringPopulation[offspringIndex]->getObjectiveValue() ) 
            {
                offspringPopulation[offspringIndex]->insertPartialSolution(partial_offspring);
                updateElitistAndCheckVTR(offspringPopulation[offspringIndex]);
                solutionHasChanged = true;
            }
        }
        delete partial_offspring;
        if (solutionHasChanged)
            break;
    }

    if (!solutionHasChanged)
    {
        *offspringPopulation[offspringIndex] = sharedInformationPointer->elitist;
    }

    return solutionHasChanged;
}

bool Population::FISingleFOS(size_t offspringIndex, size_t FOSIndex)
{
    bool solutionHasChanged = 0;

    int ind = FOSInstance->FOSorder[FOSIndex];
    // vec_t<char> touchedGenes = vec_t<char>(FOSInstance->elementSize(ind));
    bool donorEqualToOffspring = true;

    vec_t<char> touchedGenes = vec_t<char>(offspringPopulation[offspringIndex]->variables.size());
    // First copy all original genes to touchedGenes
    for(size_t i = 0; i < offspringPopulation[offspringIndex]->variables.size(); i++)
        touchedGenes[i] = offspringPopulation[offspringIndex]->variables[i];

    for (size_t j = 0; j < FOSInstance->elementSize(ind); j++)
    {
        int variableFromFOS = FOSInstance->FOSStructure[ind][j];
        touchedGenes[variableFromFOS] = sharedInformationPointer->elitist.variables[variableFromFOS]; // RUBEN: shouldn't touchedGenes[j] be touchedGenes[variableFromFOS]? -> Made this change myself
        if (offspringPopulation[offspringIndex]->variables[variableFromFOS] != touchedGenes[variableFromFOS]) 
            donorEqualToOffspring = false;
    }
    // gomea::partial_solution_t<char> *partial_offspring = new gomea::partial_solution_t<char>(touchedGenes, FOSInstance->FOSStructure[ind]);
    solution_mixed *offspring = new solution_mixed(touchedGenes, offspringPopulation[offspringIndex]->c_variables);

    if (!donorEqualToOffspring)
    {
        problemInstance->evaluate(offspring);
        // problemInstance->evaluatePartialSolution(offspringPopulation[offspringIndex], partial_offspring );

        if (offspring->getObjectiveValue() < offspringPopulation[offspringIndex]->getObjectiveValue() ) // RUBEN this was assuming maximization (>), changed to assume minimization. TODO: more elegant solution, maybe using optimization_mode of problemInstance
        {
            offspringPopulation[offspringIndex] = offspring;
            updateElitistAndCheckVTR(offspringPopulation[offspringIndex]);
            solutionHasChanged = true;
        }
        // else { // TODO: figure out if this should be done; not sure if this would cause issues with the reference to c_variables that shouldn't be deleted
        //     delete offspring;
        // }
    }

    if (!solutionHasChanged)
    {
        *offspringPopulation[offspringIndex] = sharedInformationPointer->elitist;
    }

    return solutionHasChanged;
}

/*void Population::evaluateSolution(solution_mixed *parent, gomea::partial_solution_t<char> *solution ) 
{
    checkTimeLimit();

    //cout << "before eval" << solution -> fitness << endl;
    if (config->usePartialEvaluations && solution != NULL)
    {
        problemInstance->calculateFitnessPartialEvaluations(solution, solutionBefore, touchedGenes, fitnessBefore);
        sharedInformationPointer->numberOfEvaluations += (double)touchedGenes.size() / problemInstance->number_of_variables;
    }
    else
    {
        problemInstance->calculateFitness(solution);
        sharedInformationPointer->numberOfEvaluations += 1;
    }

    updateElitistAndCheckVTR(solution);
}*/

/*void Population::evaluateSolution(solution_mixed *solution, solution_mixed *solutionBefore, vec_t<int> &touchedGenes, double fitnessBefore)
{
    checkTimeLimit();

    // Do the actual evaluation
    archiveRecord searchResult;
    
    if (config->saveEvaluations)
        sharedInformationPointer->evaluatedSolutions->checkAlreadyEvaluated(solution->variables, &searchResult);
    
    if (searchResult.isFound)
        solution->setObjectiveValue(searchResult.value);
    else
    { 
        //cout << "before eval" << solution -> fitness << endl;
        if (config->usePartialEvaluations && solutionBefore != NULL)
        {
            assert(0);
            // TODO
            //problemInstance->evaluatePartialSolution(solution, solutionBefore, touchedGenes, fitnessBefore);
            sharedInformationPointer->numberOfEvaluations += (double)touchedGenes.size() / problemInstance->number_of_variables;
        }
        else
        {
            assert(0);
            // TODO
            //problemInstance->evaluate(solution);
            sharedInformationPointer->numberOfEvaluations += 1;
        }

        if (config->saveEvaluations)
            sharedInformationPointer->evaluatedSolutions->insertSolution(solution->variables, solution->getObjectiveValue());
    }

    updateElitistAndCheckVTR(solution);
}*/

void Population::checkTimeLimit()
{
    if ( config->maximumNumberOfSeconds > 0 && utils::getElapsedTimeSeconds(sharedInformationPointer->startTime) > config->maximumNumberOfSeconds)
    {
        terminated = true;
        throw utils::customException("time");
    }
}

void Population::updateElitistAndCheckVTR(solution_mixed *solution)
{
    // RUBEN TODO: use problemInstance elitist here, to be consistent with iAMaLGaM as well.
    /* Update elitist solution */
    //if (sharedInformationPointer->firstEvaluationEver || (solution->getObjectiveValue() > sharedInformationPointer->elitist.getObjectiveValue()))
    if (sharedInformationPointer->firstEvaluationEver || problemInstance->betterFitness(solution,&sharedInformationPointer->elitist) )
    {
        if(config->printNewElitists)
        {
            cout << "New elitist solution found: " << solution->getObjectiveValue() << ", c_variables: ";
            for(auto const& val : solution->c_variables)
                cout << val << " ";
            cout << "\t variables: ";
            for(int dval : solution->variables)
                cout << dval << " ";
            cout << endl;
        }
        sharedInformationPointer->elitistSolutionHittingTimeMilliseconds = utils::getElapsedTimeMilliseconds(sharedInformationPointer->startTime);
        sharedInformationPointer->elitistSolutionHittingTimeEvaluations = problemInstance->number_of_evaluations;

        sharedInformationPointer->elitist = *solution;
		sharedInformationPointer->elitistFitness = solution->getObjectiveValue();
		sharedInformationPointer->elitistConstraintValue = solution->getConstraintValue();
        
        /* Check the VTR */
        if (problemInstance->use_vtr && solution->getObjectiveValue() <= problemInstance->vtr + 1e-10) // RUBEN: was >=, changed to <= (since I'm assuming minimization, TODO should probably have a robuster solution (based on probleminstance optimization_mode?))
        {
            writeStatisticsToFile(config->folder, sharedInformationPointer->elitistSolutionHittingTimeEvaluations, sharedInformationPointer->elitistSolutionHittingTimeMilliseconds, solution, populationSize);
            writeElitistSolutionToFile(config->folder, sharedInformationPointer->elitistSolutionHittingTimeEvaluations, sharedInformationPointer->elitistSolutionHittingTimeMilliseconds, solution);
            cout << "VTR HIT! (popsize: " << populationSize << ", stats -> (generation #evals): " << numberOfGenerations << " " << sharedInformationPointer->elitistSolutionHittingTimeEvaluations << ")\n";
            terminated = true;
            throw utils::customException("vtr");
        }
    
        //writeStatisticsToFile(config->folder, sharedInformationPointer->elitistSolutionHittingTimeEvaluations, sharedInformationPointer->elitistSolutionHittingTimeMilliseconds, solution);
        //if( config->writeElitists )
			//writeElitistSolutionToFile(config->folder, sharedInformationPointer->elitistSolutionHittingTimeEvaluations, sharedInformationPointer->elitistSolutionHittingTimeMilliseconds, solution);
    }

    sharedInformationPointer->firstEvaluationEver = false;
}

void Population::generateDiscretePopulation(int FOS_index) 
{
    // checkForDuplicate("START OF DISCRETE");
    generateSingleOffspring(FOS_index);
    // checkIndividualSolvedDiscrete();
    // checkForDuplicate("DISCRETE");
}


void Population::learnContinuousModel()
{
    // checkForDuplicate("CONTINUOUS 1");
    iamalgamInstance->learnContinuousModel();
    // checkForDuplicate("CONTINUOUS 2");

    iamalgamInstance->generateAndEvaluateNewSolutionsToFillPopulations();
    // iamalgamInstance->generateNewPopulation();
    writePopulationToFile(config->folder, population, "POPULATION After generating CONTINUOUS population ----------------------------------", config->logDebugInformation);
    // checkForDuplicate("CONTINUOUS 3");
    
    evaluateAllSolutionsInPopulation();
    // checkForDuplicate("CONTINUOUS 4");
    iamalgamInstance->number_of_generations++;
    iamalgamInstance->adaptDistributionMultipliersForOnePopulation();
}


// void Population::generateNewContinuousPopulation()
// {
//     iamalgamInstance->generateNewPopulation();
// }

void Population::checkForDuplicate(string message)
{
    // cout << setprecision(15);
    // for(size_t i = 0; i < populationSize; i++)
    // {
    //     for(size_t j = i+1; j < populationSize; j++)
    //     {
    //         if(population[i]->getObjectiveValue() == population[j]->getObjectiveValue())
    //         {
    //             cout << "[" + message + "] Duplicate found: " << population[i]->getObjectiveValue() << ", " << population[j]->getObjectiveValue() << "\t(" << i << "," << j << ")" << endl;
    //             cout << "[" + message + "] i: \tc_variables: ";
    //             for(double val : population[i]->c_variables)
    //                 cout << val << " ";
    //             cout << "\t variables: ";
    //             for(int dval : population[i]->variables)
    //                 cout << dval << " ";
    //             cout << endl;
    //             cout << "[" + message + "] j: \tc_variables: ";
    //             for(double val : population[j]->c_variables)
    //                 cout << val << " ";
    //             cout << "\t variables: ";
    //             for(int dval : population[j]->variables)
    //                 cout << dval << " ";
    //             cout << endl;
    //             // exit(0);
    //         }
    //         else if(allContinuousEqual(population[i], population[j]))
    //         {
    //             cout << "[ " << message << "] Two solutions with same continuous variables found: " << "(" << i << "," << j << ")\t i: " << population[i]->getObjectiveValue() << ", j: " << population[j]->getObjectiveValue() << endl;
    //         }
    //     }
    // }
}

bool Population::allContinuousEqual(solution_mixed *a, solution_mixed *b)
{
    for(size_t i = 0; i < a->c_variables.size(); i++)
    {
        if(a->c_variables[i] != b->c_variables[i])
            return false;
    }
    return true;
}

void Population::checkIndividualSolvedDiscrete()
{
    // for(size_t i = 0; i < populationSize; i++)
    // {
    //     bool all_one = true;
    //     for(size_t j = 0; j < population[i]->variables.size(); j++)
    //     {
    //         if(population[i]->variables[j] == '\000')
    //         {
    //             all_one = false;
    //             break;
    //         }
    //     }
    //     if(all_one)
    //     {
    //         cout << "Found individual that solved discrete part. i: " << i << ", gen: " << numberOfGenerations << endl;
    //         exit(0);
    //     }
    // }
}

}}